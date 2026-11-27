"""Independent saved-data verification of timeline eligibility, exposure and calibration arithmetic."""
from pathlib import Path
import hashlib,importlib.util,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
POLICIES=['frozen_naive','age_only','age_and_state']
CASES=['stable','announced_line','unannounced_line','repeated_signals','calibration_intrusions']
def load(p):return json.loads(p.read_text(encoding='utf-8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def hit(c,s):return s is not None and max(c['t_start'],s[0])<=min(c['t_end'],s[1]) and max(c['f_min'],s[2])<=min(c['f_max'],s[3])
def exposure(rows):
    pairs={(r['sample_id'],r['start_s'],r['end_s']) for r in rows};assert len(pairs)==len(rows)
    ranges=sorted((a,b) for _,a,b in pairs);assert all(math.isfinite(a) and math.isfinite(b) and b>a for a,b in ranges)
    assert all(ranges[j-1][1]<=ranges[j][0] for j in range(1,len(ranges)))
    return sum(b-a for a,b in ranges)
def main():
    out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/background_lifetime'
    m=load(out/'manifest.json');assert m['status']=='COMPLETE' and len(m['cases'])==15
    report={'status':'RUNNING','input_records':0,'calibrations':0,'diagnostic_measurements':0,'policy_records':0,'candidate_records':0,'nonempty_diagnostics_during_abstention':0,'summaries':[],
        'scope':'Three paired continuous-sequence seeds per case; no scientific lifetime/FAR certification'}
    for entry in m['cases']:
        seed=entry['seed'];case=entry['case'];dest=out/f'{seed}_{case}'
        f=load(dest/'input_features.json');rows=load(dest/'timeline.json');measurements=load(dest/'measurements.json');membership=load(dest/'calibration_membership.json')
        assert len(f)==800 and [r['index'] for r in f]==list(range(800)) and len({r['sample_id'] for r in f})==800
        assert len(rows)==516*3 and len(measurements)==516+116
        stats={};nulls={};calibrations={}
        for stage,lo,hi in [('initial',0,284),('replacement',400,684)]:
            cm=membership[stage];assert cm['source_indices']==list(range(lo,hi)) and cm['fit_index']==hi-1 and cm['first_use_index']==hi
            assert cm['source_start_s']==4.*lo and cm['source_end_s']==4.*hi
            assert cm['source_input_hashes']==[r['input_sha256'] for r in f[lo:hi]] and cm['sample_ids']==[r['sample_id'] for r in f[lo:hi]]
            assert {r['external_state'] for r in f[lo:hi]}=={cm['source_state']}
            assert all(r['commanded_background'] for r in f[lo:hi])
            cal=load(dest/f'calibration_{stage}.json');assert cal['n_norm_samples']==85 and cal['n_threshold_samples']==199 and cal['target_fap']==.005
            ss={}
            for name,s in cal['scale_stats'].items():
                raw=[r['raw_maxima'][name] for r in f[lo:lo+85]];med=float(np.median(raw));q=float(np.quantile(raw,.95));spread=max(q-med,1.)
                assert (s['median_null_score'],s['q95_null_score'],s['spread_scale'])==(med,q,spread);ss[name]=(med,spread)
            def score(r,ss=ss):return max((r['raw_maxima'][name]-med)/spread for name,(med,spread) in ss.items())
            null=sorted(score(r) for r in f[lo+85:hi]);assert null==cal['sorted_global_null_scores']
            # With199B andalpha.005, the admissible threshold is exactly the maximum.
            assert 1/200.==.005 and cal['global_threshold']==null[-1]
            fr=load(dest/f'freeze_{stage}.json');assert fr['decisions_using_this_calibration_so_far']==0
            assert fr['fit_index']==hi-1 and fr['earliest_use_index']==hi
            assert fr['calibration_hash']==hashlib.sha256((dest/f'calibration_{stage}.json').read_bytes()).hexdigest()
            stats[stage]=ss;nulls[stage]=null;calibrations[stage]=cal;report['calibrations']+=1
        for mid,r in measurements.items():
            i=r['index'];stage=r['calibration_id'];cm=membership[stage];raw=f[i];ss=stats[stage];null=nulls[stage];cal=calibrations[stage]
            assert mid==f'{i}_{stage}' and cm['fit_index']<i and cm['source_end_s']<=raw['start_s'] and i not in cm['source_indices']
            score=max((raw['raw_maxima'][name]-med)/spread for name,(med,spread) in ss.items());p=(1+sum(x>=score for x in null))/200.
            assert score==r['score'] and p==r['diagnostic_p'] and r['threshold']==cal['global_threshold']
            assert r['detected']==bool(score>r['threshold'] and p<=.005)==bool(r['candidates'])
            associated=False
            for c in r['candidates']:
                med,spread=ss[c['scale_name']];cp=(1+sum(x>=c['standardized_score'] for x in null))/200.
                assert c['standardized_score']==(c['raw_anomaly_score']-med)/spread and c['standardized_score']>r['threshold']
                assert cp==c['empirical_p_value'] and cp<=.005
                assert 0<=c['t_start']<=c['t_end']<=4 and 0<=c['f_min']<=c['f_max']<=512
                associated=associated or hit(c,raw['support']);report['candidate_records']+=1
            assert r['associated']==bool(associated);report['diagnostic_measurements']+=1
        per_policy={p:[] for p in POLICIES}
        for r in rows:
            i=r['index'];policy=r['policy'];raw=f[i];stage='initial' if policy=='frozen_naive' or i<684 else 'replacement';cm=membership[stage]
            assert all(r[k]==raw[k] for k in raw if k!='raw_maxima')
            assert r['diagnostic_calibration_id']==stage and r['measurement_id']==f'{i}_{stage}' and not r['profile_validated']
            value=measurements[r['measurement_id']]
            # Independently derive the complete fixed schedule, not by calling eligibility().
            state_ok=cm['source_state']==raw['external_state']
            age_ok=raw['end_s']-cm['source_start_s']<=2400.
            available=policy=='frozen_naive' or (age_ok and (policy=='age_only' or state_ok))
            reason='NAIVE_REUSE_UNQUALIFIED' if policy=='frozen_naive' else ('STATE_MISMATCH' if policy=='age_and_state' and not state_ok else ('AGE_EXPIRED' if not age_ok else 'ELIGIBLE_UNVALIDATED_POLICY'))
            assert r['available']==available and r['reason']==reason
            assert r['policy_calibration_id']==(stage if available else None)
            assert r['threshold']==(value['threshold'] if available else None) and r['p']==(value['diagnostic_p'] if available else None)
            assert r['detected']==bool(available and value['detected']) and r['associated']==bool(available and value['associated'])
            assert r['replacement_fitted_for_next']==(i==683 and policy!='frozen_naive')
            if not available and value['candidates']:report['nonempty_diagnostics_during_abstention']+=1
            per_policy[policy].append(r);report['policy_records']+=1
        for policy,seq in per_policy.items():
            assert len(seq)==516 and [r['index'] for r in seq]==list(range(284,800))
            summary=entry['summary'][policy];assert summary==load(dest/'summary.json')[policy]
            unavailable=[r for r in seq if not r['available']]
            expected=[] if policy=='frozen_naive' else list(range(360 if policy=='age_and_state' and case=='announced_line' else 600,684))
            assert [r['index'] for r in unavailable]==expected and summary['unavailable_windows']==len(expected)
            assert summary['available_windows']==516-len(expected) and summary['target_windows']==516
            neg=[r for r in seq if r['signal_kind'] is None];an=[r for r in neg if r['available']];un=[r for r in neg if not r['available']]
            assert summary['available_noise_windows']==len(an) and summary['unavailable_noise_windows']==len(un)
            assert summary['available_noise_seconds']==exposure(an) and summary['unavailable_noise_seconds']==exposure(un)
            assert summary['noise_alarms']==sum(r['detected'] for r in neg)
            assert summary['conditional_noise_rate']==(sum(r['detected'] for r in neg)/len(an) if an else None)
            for kind in ['wanted','intrusion']:
                signals=[r for r in seq if r['signal_kind']==kind];s=summary['signals'][kind]
                assert s=={'total':len(signals),'detected':sum(r['detected'] for r in signals),'associated':sum(r['associated'] for r in signals),
                    'missed_available':sum(r['available'] and not r['associated'] for r in signals),'abstained':sum(not r['available'] for r in signals)}
            late=[r['index'] for r in seq if r['signal_kind']=='wanted' and r['index']>=684 and r['associated']]
            assert summary['first_late_associated_wanted_index']==(late[0] if late else None)
            assert summary['availability_restored_index']==(684 if expected else None)
            assert summary['diagnostic_candidates_during_abstention']==sum(len(measurements[r['measurement_id']]['candidates']) for r in unavailable)
            assert not summary['validated_recovery_claim']
            for name,lo,hi in [('before_change',284,360),('change_to_replacement',360,684),('after_replacement',684,800)]:
                part=[r for r in seq if lo<=r['index']<hi];negative=[r for r in part if r['signal_kind'] is None];positive=[r for r in part if r['signal_kind'] is not None]
                assert summary['eras'][name]=={'noise_total':len(negative),'available_noise':sum(r['available'] for r in negative),'noise_alarms':sum(r['detected'] for r in negative),
                    'signal_total':len(positive),'signal_associated':sum(r['associated'] for r in positive),'signal_abstained':sum(not r['available'] for r in positive)}
        report['summaries'].append(entry);report['input_records']+=len(f)
    # Paired line waveforms are identical. Metadata changes must not alter core measurements.
    for seed in m['seeds']:
        a=out/f'{seed}_announced_line';b=out/f'{seed}_unannounced_line'
        af=load(a/'input_features.json');bf=load(b/'input_features.json')
        assert [r['input_sha256'] for r in af]==[r['input_sha256'] for r in bf]
        assert load(a/'measurements.json')==load(b/'measurements.json')
        ar=load(a/'timeline.json');br=load(b/'timeline.json')
        for x,y in zip(ar,br):
            if x['policy'] in ['frozen_naive','age_only']:
                assert all(x[k]==y[k] for k in ['available','threshold','p','detected','associated','diagnostic_calibration_id'])
    # Regeneration checks preservation and fixed schedule, not scientific simulator validity.
    spec=importlib.util.spec_from_file_location('lifetime',ROOT/'experiments/run_background_lifetime.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    assert (out/'run_background_lifetime.py').read_bytes()==Path(driver.__file__).read_bytes()
    regenerated=0
    for seed in driver.SEEDS:
        for case in driver.CASES:
            saved=load(out/f'{seed}_{case}/input_features.json')
            for (obs,meta),r in zip(driver.generate(seed,case),saved):
                assert all(r[k]==v for k,v in meta.items()) and driver.key(obs)==r['input_sha256'];regenerated+=1
    assert regenerated==report['input_records']==12000
    for name,h in {**m['source_sha256'],**m['dependency_sha256']}.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    assert report['nonempty_diagnostics_during_abstention']>0
    report.update(status='PASS',regenerated_inputs=regenerated,source_hashes_unchanged=True,uncached_parity_checks=m['parity_checks'],focused_fixture_checks=load(out/'focused_checks.json')['count'])
    save(out/'independent_verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='summaries'},indent=2))
if __name__=='__main__':main()
