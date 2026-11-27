"""Independent numerical/accounting oracle for warning-only monitor evidence."""
from pathlib import Path
from functools import lru_cache
import hashlib,importlib.util,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
CASES=['stable','level_step','spectral_line','sparse_transients','repeated_transients','persistent_signal']
def load(p):return json.loads(p.read_text(encoding='utf-8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def canonical(x):return json.dumps(x,sort_keys=True,allow_nan=False).encode()
def associated(c,s):return s is not None and max(c['t_start'],s[0])<=min(c['t_end'],s[1]) and max(c['f_min'],s[2])<=min(c['f_max'],s[3])
@lru_cache(None)
def exact_ks_tail(dnum,n=12,m=284):
    total=math.comb(n+m,n);previous=[0]*(m+1)
    for i in range(n+1):
        row=[0]*(m+1)
        for j in range(m+1):
            if abs(m*i-n*j)>=dnum:continue
            row[j]=1 if i==0 and j==0 else (previous[j] if i else 0)+(row[j-1] if j else 0)
        previous=row
    return (total-previous[m])/total

def ks(local,reference):
    points=np.unique(np.r_[local,reference]);a=np.searchsorted(np.sort(local),points,side='right');b=np.searchsorted(np.sort(reference),points,side='right')
    dnum=int(np.max(np.abs(len(reference)*a-len(local)*b)))
    return dnum/(len(local)*len(reference)),exact_ks_tail(dnum,len(local),len(reference))
def feature_oracle(x):
    x=np.asarray(x,dtype=float);x=x-float(np.median(x));a=np.abs(x);sigma=max(1.4826*float(np.median(a)),1e-18)
    window=.5-.5*np.cos(2*np.pi*np.arange(128)/128.)
    frames=np.lib.stride_tricks.sliding_window_view(x,128)[::64]
    fft=np.fft.rfft(frames*window[None,:],axis=1)/np.sum(window)
    psd=np.median(np.abs(fft)**2,axis=0)/np.log(2.);freq=np.fft.rfftfreq(128,1/1024.)
    values=[np.log(sigma),float(np.quantile(a,.99))/sigma]
    for lo,hi in [(20,80),(80,160),(160,280),(280,450)]:values.append(np.log(max(float(np.median(psd[(freq>=lo)&(freq<=hi)])),1e-30)))
    return np.array(values)

def episode_oracle(flags):
    starts=[i for i,v in enumerate(flags) if v and (i==0 or not flags[i-1])]
    out=[]
    for a in starts:
        b=a+1
        while b<len(flags) and flags[b]:b+=1
        out.append({'start_index':284+a,'end_index_exclusive':284+b,'windows':b-a,'duration_s':4.*(b-a)})
    return out

def main():
    out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/regime_warning_monitor'
    m=load(out/'manifest.json');assert m['status']=='COMPLETE' and len(m['cases'])==30
    cfg=m['monitor_config'];assert cfg['local_window_count']==12 and cfg['mode']=='causal_trailing' and cfg['stft_nperseg']==128
    report={'status':'RUNNING','reference_models':0,'detector_calibrations':0,'monitor_records':0,'search_records':0,'candidate_records':0,'sidecar_records':0,'ks_feature_checks':0,'max_ks_p_absolute_error':0.,'input_records':0,'regenerated_inputs':0,'manual_numpy_feature_checks':0,'case_summaries':[],'gates':{},'scope':'Five paired seed units; descriptive fixed horizon, no automatic noise classification or FAP certification'}
    for seed in m['seeds']:
        sd=out/str(seed);ref=load(sd/'monitor_reference.json');cal=load(sd/'detector_calibration.json');freeze=load(sd/'freeze.json')
        assert freeze['target_decisions_measured']==0
        assert all(hashlib.sha256((sd/f).read_bytes()).hexdigest()==h for f,h in freeze['hashes'].items())
        reference=np.array(ref['ref_matrix']);assert reference.shape==(284,6)
        med=np.median(reference,axis=0);mad=1.4826*np.median(np.abs(reference-med[None,:]),axis=0);std=np.std(reference,axis=0)
        scale=np.where(mad>1e-9,mad,np.where(std>1e-9,std,1.))
        np.testing.assert_array_equal(med,ref['ref_medians']);np.testing.assert_array_equal(scale,ref['ref_mad_scales'])
        source=load(sd/'stable/input_features.json')[:284]
        assert np.array_equal(reference,np.array([r['monitor_features'] for r in source]))
        assert freeze['reference_input_hashes']==[r['input_sha256'] for r in source]
        assert cal['n_norm_samples']==85 and cal['n_threshold_samples']==199 and cal['target_fap']==.005
        ss={}
        for name,s in cal['scale_stats'].items():
            a=[r['raw_maxima'][name] for r in source[:85]];loc=float(np.median(a));q=float(np.quantile(a,.95));spread=max(q-loc,1.)
            assert (s['median_null_score'],s['q95_null_score'],s['spread_scale'])==(loc,q,spread);ss[name]=(loc,spread)
        def search_score(r):return max((r['raw_maxima'][name]-loc)/spread for name,(loc,spread) in ss.items())
        null=sorted(search_score(r) for r in source[85:]);assert null==cal['sorted_global_null_scores'] and cal['global_threshold']==null[-1]
        report['reference_models']+=1;report['detector_calibrations']+=1;stable_flags=None
        for case in CASES:
            dest=sd/case;inputs=load(dest/'input_features.json');assess=load(dest/'monitor_assessments.json');search=load(dest/'search_measurements.json');sidecar=load(dest/'warning_sidecar.json');summary=load(dest/'summary.json')
            assert len(inputs)==800 and len(assess)==len(search)==len(sidecar)==516
            assert [r['input_sha256'] for r in inputs[:284]]==[r['input_sha256'] for r in source]
            assert all(r['index']==i for i,r in enumerate(inputs))
            target=inputs[284:];features=np.array([r['monitor_features'] for r in target]);z=(features-med[None,:])/scale[None,:];absmax=np.max(np.abs(z),axis=1)
            cusum=0.;flags=[];eligible=[]
            for i,(raw,a,s,w) in enumerate(zip(target,assess,search,sidecar)):
                assert a['index']==i and a['sample_id']==raw['sample_id'] and not a['future_data_used'] and a['mode']=='causal_trailing'
                assert a['min_latency_windows']==1 and a['min_latency_seconds']==4.
                cusum=max(0.,cusum+(float(absmax[i])-1.8-cfg['cusum_k']))
                np.testing.assert_allclose(a['cusum_stat'],cusum,rtol=1e-12,atol=1e-10)
                lo=max(0,i-11);assert a['diagnostics']['local_span']==[lo,i+1]
                enough=i>=11;eligible.append(enough)
                if not enough:
                    assert a['credibility_status']=='WARMUP_INSUFFICIENT_HISTORY' and a['ks_min_p_bonferroni'] is None and not a['regime_shifted']
                    np.testing.assert_allclose(a['max_feature_z_drift'],absmax[i],atol=1e-12)
                else:
                    local=features[lo:i+1];shift=np.abs(np.median(local,axis=0)-med)/scale
                    pvals=[ks(local[:,k],reference[:,k])[1] for k in range(6)];bonf=min(1.,6*min(pvals))
                    err=abs(bonf-a['ks_min_p_bonferroni']);report['max_ks_p_absolute_error']=max(report['max_ks_p_absolute_error'],err)
                    assert err<5e-12,(seed,case,i,bonf,a['ks_min_p_bonferroni'])
                    ks_trigger=bonf<cfg['ks_alpha'];z_trigger=float(np.max(shift))>=cfg['drift_z_threshold'];cp=cusum>=cfg['cusum_h']
                    assert a['diagnostics']['ks_triggered']==ks_trigger and a['diagnostics']['z_triggered']==z_trigger and a['diagnostics']['cusum_triggered']==cp
                    assert a['change_point_flag']==cp and a['regime_shifted']==bool(ks_trigger or z_trigger or cp)
                    assert a['credibility_status']==('REGIME_SHIFT_DETECTED' if a['regime_shifted'] else 'STATIONARY_CREDIBLE')
                    np.testing.assert_allclose(a['max_feature_z_drift'],max(shift),atol=1e-12)
                    dominant=int(np.argmin(pvals)) if ks_trigger else int(np.argmax(shift))
                    assert a['dominant_shift_feature']==ref['feature_names'][dominant]
                    assert a['diagnostics']['per_feature_z_shift']=={name:round(float(value),3) for name,value in zip(ref['feature_names'],shift)}
                    report['ks_feature_checks']+=6
                flags.append(bool(a['regime_shifted']));report['monitor_records']+=1
                value=search_score(raw);p=(1+sum(x>=value for x in null))/200.
                assert s['index']==raw['index'] and s['sample_id']==raw['sample_id'] and s['score']==value and s['p']==p and s['threshold']==cal['global_threshold'] and s['available']
                assert s['detected']==bool(value>cal['global_threshold'] and p<=.005)==bool(s['candidates'])
                hit=False
                for c in s['candidates']:
                    loc,spread=ss[c['scale_name']];cpv=(1+sum(x>=c['standardized_score'] for x in null))/200.
                    assert c['standardized_score']==(c['raw_anomaly_score']-loc)/spread and cpv==c['empirical_p_value'] and cpv<=.005 and c['standardized_score']>s['threshold']
                    assert 0<=c['t_start']<=c['t_end']<=4 and 0<=c['f_min']<=c['f_max']<=512 and c['n_bins']>=2
                    hit=hit or associated(c,raw['support']);report['candidate_records']+=1
                assert s['associated']==bool(hit)
                assert w['index']==raw['index'] and w['sample_id']==raw['sample_id'] and w['warning_only'] and w['warning']==a['regime_shifted'] and w['monitor_eligible']==enough
                assert w['search_measurement_sha256']==hashlib.sha256(canonical(s)).hexdigest() and not w['background_updated'] and not w['candidates_suppressed']
                report['search_records']+=1;report['sidecar_records']+=1
            preservation=load(dest/'preservation.json');h=hashlib.sha256(canonical(search)).hexdigest()
            assert preservation['search_sha256_before_annotation']==preservation['search_sha256_after_annotation']==h and preservation['no_detector_gating']
            ep=episode_oracle(flags);assert summary['warning_episodes']==ep and summary['warning_windows']==sum(flags) and summary['any_warning']==any(flags)
            assert summary['monitor_eligible_windows']==sum(eligible)==505 and summary['warmup_windows']==11 and summary['target_windows']==516
            assert summary['prechange_warning_windows']==sum(flags[:76]) and summary['warning_active_before_change']==flags[75]
            post=next((i+284 for i,v in enumerate(flags) if i>=76 and v),None);new=next((e['start_index'] for e in ep if e['start_index']>=360),None)
            assert summary['first_post_change_warning_index']==post and summary['first_new_episode_after_change']==new
            assert summary['first_post_change_warning_delay_s']==(4.*(post+1-360) if post is not None else None)
            assert summary['first_new_episode_delay_s']==(4.*(new+1-360) if new is not None else None)
            signal_indices=[i for i,r in enumerate(target) if r['signal_kind'] is not None]
            assert summary['signal_windows']==len(signal_indices) and summary['signal_warning_windows']==sum(flags[i] for i in signal_indices)
            assert summary['signal_warmup_windows']==sum(not eligible[i] for i in signal_indices)
            assert summary['associated_search_signal_windows']==sum(search[i]['associated'] for i in signal_indices)
            assert summary['associated_recoveries_that_warning_gate_would_remove']==sum(search[i]['associated'] and flags[i] for i in signal_indices)
            assert summary['search_noise_alarms']==sum(s['detected'] for r,s in zip(target,search) if r['signal_kind'] is None)
            assert summary['search_noise_windows']==sum(r['signal_kind'] is None for r in target)
            assert summary['detector_available_windows']==516 and not summary['background_updated'] and not summary['candidates_suppressed']
            assert summary['trigger_windows']=={key:sum(bool(a['diagnostics'].get(key,False)) for a in assess) for key in ['ks_triggered','z_triggered','cusum_triggered']}
            if case=='stable':stable_flags=flags
            assert summary['extra_warning_windows_vs_stable']==sum(a and not b for a,b in zip(flags,stable_flags))
            assert summary['suppressed_warning_windows_vs_stable']==sum(b and not a for a,b in zip(flags,stable_flags))
            last=max(signal_indices) if signal_indices else None
            assert summary['extra_warning_windows_after_last_signal']==(sum(flags[i] and not stable_flags[i] for i in range(last+1,516)) if last is not None else 0)
            report['case_summaries'].append({'seed':seed,'case':case,'summary':summary});report['input_records']+=800
        a=sd/'spectral_line';b=sd/'persistent_signal'
        af=load(a/'input_features.json');bf=load(b/'input_features.json')
        assert [r['input_sha256'] for r in af]==[r['input_sha256'] for r in bf]
        for x,y in zip(load(a/'monitor_assessments.json'),load(b/'monitor_assessments.json')):assert {k:v for k,v in x.items() if k!='sample_id'}=={k:v for k,v in y.items() if k!='sample_id'}
        for x,y in zip(load(a/'search_measurements.json'),load(b/'search_measurements.json')):assert {k:v for k,v in x.items() if k not in ['sample_id','associated']}=={k:v for k,v in y.items() if k not in ['sample_id','associated']}
    spec=importlib.util.spec_from_file_location('warning_pilot',ROOT/'experiments/run_regime_warning_monitor.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    assert Path(driver.__file__).read_bytes()==(out/'run_regime_warning_monitor.py').read_bytes()
    for seed in driver.SEEDS:
        for case in driver.CASES:
            saved=load(out/str(seed)/case/'input_features.json')
            for i,((obs,meta),r) in enumerate(zip(driver.generate(seed,case),saved)):
                assert all(r[k]==value for k,value in meta.items()) and driver.key(obs)==r['input_sha256']
                if i in [284,285,286] or (case=='stable' and i in [0,1,2]):
                    np.testing.assert_allclose(feature_oracle(obs.signal),r['monitor_features'],rtol=1e-12,atol=1e-12);report['manual_numpy_feature_checks']+=1
                report['regenerated_inputs']+=1
    assert report['input_records']==report['regenerated_inputs']==24000
    for name,h in {**m['source_sha256'],**m['dependency_sha256']}.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    allcases=report['case_summaries']
    report['gates']={
      'zero_stable_warning_episodes':all(not x['summary']['warning_episodes'] for x in allcases if x['case']=='stable'),
      'all_changes_new_episode_within_48s_without_preexisting_warning':all(not x['summary']['warning_active_before_change'] and x['summary']['first_new_episode_delay_s'] is not None and x['summary']['first_new_episode_delay_s']<=48. for x in allcases if x['case'] in ['level_step','spectral_line']),
      'zero_extra_transient_warning_windows':all(x['summary']['extra_warning_windows_vs_stable']==0 for x in allcases if x['case'] in ['sparse_transients','repeated_transients'])}
    report.update(status='PASS',source_hashes_unchanged=True,production_prefix_checks=m['prefix_checks'],exact_ks_cache_entries=exact_ks_tail.cache_info().currsize,
      noise_origin_classification_supported=False)
    save(out/'independent_verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='case_summaries'},indent=2))
if __name__=='__main__':main()
