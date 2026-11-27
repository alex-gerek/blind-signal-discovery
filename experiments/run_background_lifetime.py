"""Causal timeline experiment with scheduled source acquisition and explicit expiry costs."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,json,platform,sys,time,traceback,uuid
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"));sys.path.insert(0,str(ROOT/"experiments"));sys.dont_write_bytecode=True
from run_block_tail_guard import CachedCore,replay,key,save
from run_separate_background_control import clean,associated
from bsd.generator import Observation,_scale_to_energy_snr
from bsd.noise import generate_pink_colored
from background_lifetime_controller import POLICIES,eligibility,interval_exposure,focused_checks
SEEDS=[19101,19102,19103]
CASES=["stable","announced_line","unannounced_line","repeated_signals","calibration_intrusions"]
PROTOCOL=ROOT/"docs/background_lifetime_protocol_2026-11-27.md"
def sha(x):return hashlib.sha256(x).hexdigest()
def generate(seed,case):
    rng=np.random.default_rng(np.random.SeedSequence([seed,1]))
    base=generate_pink_colored(800*4096,1024.,rng,sigma=1.,alpha=1.,f_knee=16.)[0]
    phase=np.random.default_rng(np.random.SeedSequence([seed,2])).uniform(0,2*np.pi)
    wanted=list(range(288,400,8))+list(range(688,800,8)) if case=="repeated_signals" else [304,324,344,364,384,704,724,744,764,784]
    schedule={i:("chirp" if j%2==0 else "burst","wanted") for j,i in enumerate(wanted)}
    if case=="calibration_intrusions":schedule.update({i:("chirp" if j%2==0 else "burst","intrusion") for j,i in enumerate([500,540,580,620])})
    for i in range(800):
        n=base[i*4096:(i+1)*4096].copy()
        if case in ["announced_line","unannounced_line"] and i>=360:
            t=np.arange(i*4096,(i+1)*4096)/1024.;n+=3.*np.sin(2*np.pi*110*t+phase)
        base_hash=key(Observation('base',n,1024.,4.));family=None;kind=None;support=None;x=n
        if i in schedule:
            family,kind=schedule[i];fk=0 if family=="chirp" else 1
            r=np.random.default_rng(np.random.SeedSequence([seed,3,i,fk]));s,support=clean(family,r)
            x=n+_scale_to_energy_snr(s,float(np.std(n)),18.)
        ident=f"s{seed}_{case}_{i:04d}";o=Observation(ident,x,1024.,4.)
        telemetry="B" if case=="announced_line" and i>=360 else "A"
        yield o,{'sample_id':ident,'index':i,'start_s':4.*i,'end_s':4.*(i+1),'input_sha256':key(o),'base_noise_sha256':base_hash,
            'external_state':telemetry,'commanded_background':i<284 or 400<=i<684,'family':family,'signal_kind':kind,'support':support,'snr':18. if family else None}

def summarize(rows,measurements):
    result={}
    for policy in POLICIES:
        seq=[r for r in rows if r['policy']==policy]
        noise=[r for r in seq if r['signal_kind'] is None];signals=[r for r in seq if r['signal_kind'] is not None]
        available=[r for r in noise if r['available']];unavailable=[r for r in noise if not r['available']]
        late=[r for r in signals if r['signal_kind']=='wanted' and r['index']>=684]
        restored=next((r['index'] for j,r in enumerate(seq) if j and r['available'] and not seq[j-1]['available']),None)
        eras={}
        for name,lo,hi in [('before_change',284,360),('change_to_replacement',360,684),('after_replacement',684,800)]:
            group=[r for r in seq if lo<=r['index']<hi];neg=[r for r in group if r['signal_kind'] is None];pos=[r for r in group if r['signal_kind'] is not None]
            eras[name]={'noise_total':len(neg),'available_noise':sum(r['available'] for r in neg),'noise_alarms':sum(r['detected'] for r in neg),'signal_total':len(pos),'signal_associated':sum(r['associated'] for r in pos),'signal_abstained':sum(not r['available'] for r in pos)}
        by_kind={}
        for kind in ['wanted','intrusion']:
            group=[r for r in signals if r['signal_kind']==kind]
            by_kind[kind]={'total':len(group),'detected':sum(r['detected'] for r in group),'associated':sum(r['associated'] for r in group),'missed_available':sum(r['available'] and not r['associated'] for r in group),'abstained':sum(not r['available'] for r in group)}
        result[policy]={'target_windows':len(seq),'available_windows':sum(r['available'] for r in seq),'unavailable_windows':sum(not r['available'] for r in seq),
            'available_noise_windows':len(available),'unavailable_noise_windows':len(unavailable),'available_noise_seconds':interval_exposure([(r['sample_id'],r['start_s'],r['end_s']) for r in available]),
            'unavailable_noise_seconds':interval_exposure([(r['sample_id'],r['start_s'],r['end_s']) for r in unavailable]),'noise_alarms':sum(r['detected'] for r in noise),
            'conditional_noise_rate':sum(r['detected'] for r in noise)/len(available) if available else None,'signals':by_kind,
            'availability_restored_index':restored,'first_late_associated_wanted_index':next((r['index'] for r in late if r['associated']),None),
            'diagnostic_candidates_during_abstention':sum(len(measurements[r['measurement_id']]['candidates']) for r in seq if not r['available']),
            'eras':eras,'validated_recovery_claim':False}
    return result

def main():
    assert replay.git('rev-parse','HEAD')=='e21a7cc48ce4122804f86e3397ee3894afa90639'
    assert not replay.git('diff','--','src','app','pyproject.toml')
    global OUT
    ident=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
    OUT=ROOT/'results/background_lifetime'/ident;OUT.mkdir(parents=True,exist_ok=False)
    print('RUN_DIRECTORY',OUT,flush=True)
    dependencies=[Path(__file__),ROOT/'experiments/background_lifetime_controller.py',ROOT/'experiments/run_block_tail_guard.py',ROOT/'experiments/block_tail_guard.py',ROOT/'experiments/run_separate_background_control.py',ROOT/'experiments/run_historical_adaptive_replay.py']
    for p in dependencies:(OUT/p.name).write_bytes(p.read_bytes())
    (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes())
    source=replay.hashes(list((ROOT/'src/bsd').glob('*.py')))
    manifest={'run_id':ident,'status':'RUNNING','source_commit':replay.git('rev-parse','HEAD'),'dirty_paths':replay.git('status','--short'),'source_sha256':source,'dependency_sha256':replay.hashes(dependencies),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':replay.scipy.__version__,'age_limit_s':2400.,'age_limit_validated':False,'seeds':SEEDS,'cases':[],'parity_checks':0}
    save(OUT/'manifest.json',manifest);save(OUT/'focused_checks.json',focused_checks());start=time.monotonic()
    for seed in SEEDS:
        det=CachedCore()
        for case in CASES:
            dest=OUT/f'{seed}_{case}';dest.mkdir(exist_ok=False)
            records=list(generate(seed,case));observations=[o for o,m in records];features=[];rows=[];measurements={};cals={};freeze={};cal_meta={}
            before_parity=det.parity
            for i in [0,1,2,284,285,286,684,685,686]:det.parity_check(observations[i])
            def fit(stage,lo,hi,fit_index):
                part=observations[lo:hi];states={records[i][1]['external_state'] for i in range(lo,hi)};assert len(states)==1
                cal=det.calibrate(part,target_fap=.005,norm_fraction=.30,calibration_seed=seed,allow_unattainable_fap=True)
                assert cal.n_norm_samples==85 and cal.n_threshold_samples==199
                cals[stage]=cal;filename=f'calibration_{stage}.json';save(dest/filename,asdict(cal));freeze[filename]=sha((dest/filename).read_bytes())
                cal_meta[stage]={'calibration_id':stage,'source_indices':list(range(lo,hi)),'sample_ids':[o.sample_id for o in part],'source_input_hashes':[key(o) for o in part],
                    'source_start_s':4.*lo,'source_end_s':4.*hi,'source_state':next(iter(states)),'fit_index':fit_index,'first_use_index':hi,'n_norm':85,'n_threshold':199}
                save(dest/f'freeze_{stage}.json',{'calibration_hash':freeze[filename],'fit_index':fit_index,'earliest_use_index':hi,'decisions_using_this_calibration_so_far':0})
            fit('initial',0,284,283)
            for i,(o,meta) in enumerate(records):
                raw=det._raw_analyze_all_scales(o);features.append({**meta,'raw_maxima':{k:float(v[0]) for k,v in raw.items()}})
                if i>=284:
                    stages=['initial']+(['replacement'] if i>=684 else [])
                    for stage in stages:
                        cal=cals[stage];result=det.analyze(o,cal);candidates=[asdict(c) for c in result.candidates];mid=f'{i}_{stage}'
                        measurements[mid]={'measurement_id':mid,'index':i,'calibration_id':stage,'score':result.max_global_score,'threshold':cal.global_threshold,
                            'diagnostic_p':cal.empirical_p_value(result.max_global_score)[0],'detected':result.detected,'associated':bool(meta['support'] and any(associated(c,meta['support']) for c in candidates)),'candidates':candidates}
                    for policy in POLICIES:
                        stage='initial' if policy=='frozen_naive' or i<684 else 'replacement'
                        cm=cal_meta[stage];cal=cals[stage];mid=f'{i}_{stage}';measurement=measurements[mid]
                        ok,reason=eligibility(policy,meta['start_s'],meta['end_s'],cm['source_start_s'],cm['source_end_s'],cm['source_state'],meta['external_state'])
                        ok=bool(ok and cal.target_fap_attainable and cal.global_threshold is not None)
                        rows.append({**meta,'policy':policy,'measurement_id':mid,'diagnostic_calibration_id':stage,'policy_calibration_id':stage if ok else None,
                            'available':ok,'reason':reason,'threshold':cal.global_threshold if ok else None,'p':measurement['diagnostic_p'] if ok else None,
                            'detected':bool(ok and measurement['detected']),'associated':bool(ok and measurement['associated']),
                            'replacement_fitted_for_next':False,'profile_validated':False})
                if i==683:
                    fit('replacement',400,684,683)
                    for row in rows[-3:]:row['replacement_fitted_for_next']=row['policy']!='frozen_naive'
            assert all(sha((dest/name).read_bytes())==h for name,h in freeze.items())
            save(dest/'input_features.json',features);save(dest/'calibration_membership.json',cal_meta);save(dest/'measurements.json',measurements);save(dest/'timeline.json',rows)
            summary=summarize(rows,measurements);save(dest/'summary.json',summary)
            for policy in POLICIES:
                expected=0 if policy=='frozen_naive' else (324 if policy=='age_and_state' and case=='announced_line' else 84)
                assert summary[policy]['unavailable_windows']==expected
            manifest['cases'].append({'seed':seed,'case':case,'summary':summary,'calibration_hashes_unchanged':True})
            manifest['parity_checks']+=det.parity-before_parity;save(OUT/'manifest.json',manifest)
            print('CASE_COMPLETE',seed,case,[(p,summary[p]['noise_alarms'],summary[p]['signals']['wanted']['associated'],summary[p]['signals']['wanted']['abstained']) for p in POLICIES],flush=True)
    assert all(sha((ROOT/f).read_bytes())==h for f,h in {**source,**manifest['dependency_sha256']}.items())
    manifest.update(status='COMPLETE',elapsed_s=time.monotonic()-start,source_hashes_unchanged=True);save(OUT/'manifest.json',manifest)
    print('RUN_STATUS COMPLETE',flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        if 'OUT' in globals():save(OUT/'failure.json',{'status':'FAILED','traceback':traceback.format_exc()})
        raise
