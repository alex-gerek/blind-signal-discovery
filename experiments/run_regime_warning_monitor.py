"""Frozen default regime monitor, warning-only independent sequence study."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,json,platform,sys,time,traceback,uuid
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'));sys.dont_write_bytecode=True
from run_block_tail_guard import CachedCore,replay,key,save
from run_separate_background_control import clean,associated
from bsd.generator import Observation,_scale_to_energy_snr
from bsd.noise import generate_pink_colored
from bsd.regime_monitor import BackgroundRegimeMonitor,RegimeMonitorConfig,extract_window_features
from warning_monitor_metrics import summarize
SEEDS=list(range(20101,20106))
CASES=['stable','level_step','spectral_line','sparse_transients','repeated_transients','persistent_signal']
PROTOCOL=ROOT/'docs/regime_warning_monitor_protocol_2026-11-27.md'
def sha(x):return hashlib.sha256(x).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,default=replay.norm,allow_nan=False).encode()
def generate(seed,case):
    r=np.random.default_rng(np.random.SeedSequence([seed,1]));base=generate_pink_colored(800*4096,1024.,r,sigma=1.,alpha=1.,f_knee=16.)[0]
    phase=np.random.default_rng(np.random.SeedSequence([seed,2])).uniform(0,2*np.pi)
    wanted=[304,324,344,364,384,704,724,744,764,784] if case=='sparse_transients' else (list(range(288,800,8)) if case=='repeated_transients' else [])
    schedule={i:('chirp' if j%2==0 else 'burst') for j,i in enumerate(wanted)}
    for i in range(800):
        noise=base[i*4096:(i+1)*4096].copy();x=noise;family=None;kind=None;support=None
        if case=='level_step' and i>=360:x=2.*noise
        if case in ['spectral_line','persistent_signal'] and i>=360:
            x=noise+3.*np.sin(2*np.pi*110*np.arange(i*4096,(i+1)*4096)/1024.+phase)
            if case=='persistent_signal':family='tone';kind='persistent_signal';support=[0.,4.,110.,110.]
        if i in schedule:
            family=schedule[i];kind='transient';fk=0 if family=='chirp' else 1
            rng=np.random.default_rng(np.random.SeedSequence([seed,3,i,fk]));s,support=clean(family,rng)
            x=noise+_scale_to_energy_snr(s,float(np.std(noise)),18.)
        ident=f's{seed}_{case}_{i:04d}';o=Observation(ident,x,1024.,4.)
        yield o,{'sample_id':ident,'index':i,'start_s':4.*i,'end_s':4.*(i+1),'input_sha256':key(o),'base_noise_sha256':key(Observation('base',noise,1024.,4.)),
            'noise_change':case in ['level_step','spectral_line'] and i>=360,'family':family,'signal_kind':kind,'support':support,'snr':18. if kind=='transient' else None}

def main():
    assert replay.git('rev-parse','HEAD')=='6f4ccfac7444fb01c2b19785e4236536559ce7dd'
    assert not replay.git('diff','--','src','app','pyproject.toml')
    global OUT
    ident=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
    OUT=ROOT/'results/regime_warning_monitor'/ident;OUT.mkdir(parents=True,exist_ok=False);print('RUN_DIRECTORY',OUT,flush=True)
    dependencies=[Path(__file__),ROOT/'experiments/warning_monitor_metrics.py',ROOT/'experiments/run_block_tail_guard.py',ROOT/'experiments/block_tail_guard.py',ROOT/'experiments/run_separate_background_control.py',ROOT/'experiments/run_historical_adaptive_replay.py']
    for p in dependencies:(OUT/p.name).write_bytes(p.read_bytes())
    (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes())
    cfg=RegimeMonitorConfig();assert cfg.mode=='causal_trailing'
    sources=replay.hashes(list((ROOT/'src/bsd').glob('*.py')))
    manifest={'run_id':ident,'status':'RUNNING','source_commit':replay.git('rev-parse','HEAD'),'dirty_paths':replay.git('status','--short'),'source_sha256':sources,'dependency_sha256':replay.hashes(dependencies),
      'monitor_config':asdict(cfg),'python':platform.python_version(),'numpy':np.__version__,'scipy':replay.scipy.__version__,'seeds':SEEDS,'cases':[],'prefix_checks':0,'candidates_suppressed':False,'reference_updates':0}
    save(OUT/'manifest.json',manifest);start=time.monotonic()
    for seed in SEEDS:
        seed_dest=OUT/str(seed);seed_dest.mkdir(exist_ok=False)
        det=CachedCore();monitor=BackgroundRegimeMonitor(cfg)
        initial=list(generate(seed,'stable'))[:284];ref_obs=[o for o,m in initial]
        ref=monitor.fit_reference(ref_obs);cal=det.calibrate(ref_obs,target_fap=.005,norm_fraction=.30,calibration_seed=seed)
        save(seed_dest/'monitor_reference.json',asdict(ref));save(seed_dest/'detector_calibration.json',asdict(cal))
        frozen={name:sha((seed_dest/name).read_bytes()) for name in ['monitor_reference.json','detector_calibration.json']}
        save(seed_dest/'freeze.json',{'target_decisions_measured':0,'hashes':frozen,'reference_sample_ids':[o.sample_id for o in ref_obs],'reference_input_hashes':[key(o) for o in ref_obs]})
        stable_flags=None
        for case in CASES:
            dest=seed_dest/case;dest.mkdir(exist_ok=False);records=list(generate(seed,case));features=[]
            assert [key(o) for o,m in records[:284]]==[key(o) for o in ref_obs]
            for o,meta in records:
                raw=det._raw_analyze_all_scales(o)
                features.append({**meta,'monitor_features':extract_window_features(o,cfg.freq_bands_hz,cfg.stft_nperseg).as_vector().tolist(),'raw_maxima':{k:float(v[0]) for k,v in raw.items()}})
            target=[o for o,m in records[284:]];target_meta=[m for o,m in records[284:]]
            assessments=[asdict(a) for a in monitor.monitor_sequence(target,ref)]
            for length in [12,76,96]:
                prefix=[asdict(a) for a in monitor.monitor_sequence(target[:length],ref)]
                assert prefix==assessments[:length];manifest['prefix_checks']+=1
            search=[]
            for o,meta in records[284:]:
                result=det.analyze(o,cal);candidates=[asdict(c) for c in result.candidates]
                search.append({'sample_id':o.sample_id,'index':meta['index'],'score':result.max_global_score,'threshold':cal.global_threshold,'p':cal.empirical_p_value(result.max_global_score)[0],
                    'available':bool(cal.target_fap_attainable and cal.global_threshold is not None),'detected':result.detected,'associated':bool(meta['support'] and any(associated(c,meta['support']) for c in candidates)),'candidates':candidates})
            before=sha(canonical(search));save(dest/'search_measurements.json',search)
            sidecar=[{'index':m['index'],'sample_id':m['sample_id'],'warning_only':True,'warning':a['regime_shifted'],'monitor_eligible':a['credibility_status']!='WARMUP_INSUFFICIENT_HISTORY',
                'search_measurement_sha256':sha(canonical(s)),'background_updated':False,'candidates_suppressed':False} for m,a,s in zip(target_meta,assessments,search)]
            assert sha(canonical(search))==before
            flags=[a['regime_shifted'] for a in assessments]
            if case=='stable':stable_flags=flags
            summary=summarize(target_meta,assessments,search,stable_flags)
            save(dest/'input_features.json',features);save(dest/'monitor_assessments.json',assessments);save(dest/'warning_sidecar.json',sidecar)
            save(dest/'preservation.json',{'search_sha256_before_annotation':before,'search_sha256_after_annotation':sha(canonical(search)),'unchanged':True,'no_detector_gating':True})
            save(dest/'summary.json',summary)
            manifest['cases'].append({'seed':seed,'case':case,'summary':summary});save(OUT/'manifest.json',manifest)
            print('CASE_COMPLETE',seed,case,'warning_windows',summary['warning_windows'],'episodes',len(summary['warning_episodes']),'first_new',summary['first_new_episode_after_change'],'gate_losses',summary['associated_recoveries_that_warning_gate_would_remove'],flush=True)
        assert all(sha((seed_dest/name).read_bytes())==h for name,h in frozen.items())
    assert all(sha((ROOT/name).read_bytes())==h for name,h in {**sources,**manifest['dependency_sha256']}.items())
    manifest.update(status='COMPLETE',elapsed_s=time.monotonic()-start,source_hashes_unchanged=True)
    save(OUT/'manifest.json',manifest);print('RUN_STATUS COMPLETE',flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        if 'OUT' in globals():save(OUT/'failure.json',{'status':'FAILED','traceback':traceback.format_exc()})
        raise
