"""Fixed-source synthetic background control, unchanged core and no admission filter."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,json,platform,sys,time,traceback,uuid
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"));sys.path.insert(0,str(ROOT/"experiments"));sys.dont_write_bytecode=True
from run_block_tail_guard import CachedCore,replay,key,save
from bsd.generator import Observation,generate_clean_signal,_scale_to_energy_snr
from bsd.noise import generate_pink_colored
PROTOCOL=ROOT/"docs/separate_background_control_protocol_2026-11-27.md"
ARMS=["stream_contaminated","paired_noise_oracle","separate_matched","separate_old"]
CASES=["stationary","level_step","spectral_line","impulsive_noise"]
SEEDS=[18101,18102,18103]
def sha(x):return hashlib.sha256(x).hexdigest()
def clean(family,rng):
    if family=="chirp":
        x,a,b,c,d,params=generate_clean_signal("linear_chirp",4096,1024.,rng)
        return x,[a,b,c,d]
    x=signal.sosfiltfilt(signal.butter(4,[35.,350.],btype="bandpass",fs=1024.,output="sos"),rng.normal(size=4096))
    mask=np.zeros(4096);mask[1408:1664]=signal.windows.tukey(256,.15)
    return x*mask,[1.375,1.625,35.,350.]
def associated(candidate,support):
    a,b,c,d=support
    return candidate["t_end"]>=a and candidate["t_start"]<=b and candidate["f_max"]>=c and candidate["f_min"]<=d

def generate(seed,ci):
    case=CASES[ci]
    rng={k:np.random.default_rng(np.random.SeedSequence([seed,ci,k])) for k in range(1,7)}
    def noise(r,profile=case):
        x=generate_pink_colored(4096,1024.,r,sigma=1.,alpha=1.,f_knee=16.)[0]
        if profile=="level_step":x=x*2.
        elif profile=="spectral_line":x=x+3.*np.sin(2*np.pi*110*np.arange(4096)/1024.+r.uniform(0,2*np.pi))
        elif profile=="impulsive_noise" and r.random()<.03:
            x=x.copy();x[int(r.integers(4096))]+=float(r.choice([-12.,12.]))
        return x
    accumulated=[noise(rng[1]) for _ in range(265)]
    separate=[noise(rng[2]) for _ in range(265)]
    old=[noise(rng[3],"stationary") for _ in range(265)]
    contaminated=[x.copy() for x in accumulated]
    labels={}
    for i,family in [(100,"chirp"),(140,"burst"),(180,"chirp"),(220,"burst")]:
        x,support=clean(family,rng[5]);sigma=float(np.std(accumulated[i]))
        contaminated[i]+=_scale_to_energy_snr(x,sigma,18.)
        labels[i]={"family":family,"snr":18.,"support":support}
    sources=dict(zip(ARMS,[contaminated,accumulated,separate,old]))
    calibration={}
    for arm in ARMS:
        records=[]
        for i,x in enumerate(sources[arm]):
            ident=f"s{seed}_c{ci}_{arm}_{i:03d}";o=Observation(ident,x,1024.,4.)
            label=labels.get(i,{}) if arm=="stream_contaminated" else {}
            records.append((o,{"sample_id":ident,"part":"calibration","arm":arm,"index":i,"input_sha256":key(o),"underlying_noise_sha256":key(Observation('noise',accumulated[i] if arm in ARMS[:2] else x,1024.,4.)),**label}))
        calibration[arm]=records
    target=[]
    for i in range(200):
        o=Observation(f"s{seed}_c{ci}_noise_{i:03d}",noise(rng[4]),1024.,4.)
        target.append((o,{"sample_id":o.sample_id,"part":"noise","index":i,"input_sha256":key(o),"base_noise_sha256":key(o),"family":None,"snr":None,"unit":None,"support":None}))
    for fi,family in enumerate(["chirp","burst"]):
        for unit in range(10):
            n=noise(rng[4]);x,support=clean(family,rng[6]);sigma=float(np.std(n))
            noise_hash=key(Observation('noise',n,1024.,4.))
            for ai,snr in enumerate([8.,18.]):
                index=fi*20+unit*2+ai;o=Observation(f"s{seed}_c{ci}_signal_{index:03d}",n+_scale_to_energy_snr(x,sigma,snr),1024.,4.)
                target.append((o,{"sample_id":o.sample_id,"part":"signal","index":index,"input_sha256":key(o),"base_noise_sha256":noise_hash,"family":family,"snr":snr,"unit":unit,"support":support}))
    return calibration,target

def main():
    assert replay.git("rev-parse","HEAD")=="d049e2b1c97626b562d234a22bc12db5d51702dd"
    assert not replay.git("diff","--","src","app","pyproject.toml")
    global OUT
    ident=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    OUT=ROOT/"results/separate_background_control"/ident;OUT.mkdir(parents=True,exist_ok=False)
    print("RUN_DIRECTORY",OUT,flush=True)
    dependencies=[Path(__file__),ROOT/"experiments/run_block_tail_guard.py",ROOT/"experiments/block_tail_guard.py",ROOT/"experiments/run_historical_adaptive_replay.py"]
    for path in dependencies:(OUT/path.name).write_bytes(path.read_bytes())
    (OUT/"protocol.md").write_bytes(PROTOCOL.read_bytes())
    sources=replay.hashes(list((ROOT/"src/bsd").glob("*.py")))
    m={"run_id":ident,"status":"RUNNING","source_commit":replay.git("rev-parse","HEAD"),"dirty_paths":replay.git("status","--short"),"source_sha256":sources,"dependency_sha256":replay.hashes(dependencies),"python":platform.python_version(),"numpy":np.__version__,"scipy":replay.scipy.__version__,"seeds":SEEDS,"profiles":CASES,"cases":[],"parity_checks":0,"scope":"Synthetic fixed-source pilot, not autonomous adaptive recovery or real background-purity certification"}
    save(OUT/"manifest.json",m);start=time.monotonic()
    for seed in SEEDS:
        for ci,case in enumerate(CASES):
            dest=OUT/f"{seed}_{case}";dest.mkdir(exist_ok=False)
            det=CachedCore();calibration,target=generate(seed,ci)
            hashes={arm:{meta["input_sha256"] for _,meta in records} for arm,records in calibration.items()}
            base_hashes={arm:{meta["underlying_noise_sha256"] for _,meta in records} for arm,records in calibration.items()}
            th={meta["base_noise_sha256"] for _,meta in target}
            for arm in ARMS:
                assert len(hashes[arm])==len(base_hashes[arm])==265 and not(base_hashes[arm]&th)
            for a,b in [("paired_noise_oracle","separate_matched"),("paired_noise_oracle","separate_old"),("separate_matched","separate_old")]:assert not(base_hashes[a]&base_hashes[b])
            changed=[i for i in range(265) if key(calibration[ARMS[0]][i][0])!=key(calibration[ARMS[1]][i][0])]
            assert changed==[100,140,180,220]
            save(dest/"separation.json",{"cross_source_overlaps":0,"calibration_target_noise_overlaps":0,"oracle_changed_indices":changed,"target_distinct_noise_bases":len(th)})
            for arm in ARMS:
                for o,_ in calibration[arm][:3]:det.parity_check(o)
            calibrations={};frozen={};cal_features=[]
            for arm in ARMS:
                observations=[o for o,_ in calibration[arm]]
                for o,meta in calibration[arm]:
                    raw=det._raw_analyze_all_scales(o)
                    cal_features.append({**meta,"raw_maxima":{k:float(v[0]) for k,v in raw.items()}})
                for alpha in [.01,.005]:
                    cal=det.calibrate(observations,target_fap=alpha,norm_fraction=.25,calibration_seed=seed,allow_unattainable_fap=True)
                    name=f"{arm}_{alpha}";calibrations[name]=cal;filename=f"calibration_{name}.json"
                    save(dest/filename,asdict(cal));frozen[filename]=sha((dest/filename).read_bytes())
            save(dest/"calibration_features.json",cal_features)
            save(dest/"threshold_freeze.json",{"target_decisions_measured":0,"calibration_hashes":frozen})
            for o,_ in target[:3]:det.parity_check(o)
            tf=[]
            for o,meta in target:
                raw=det._raw_analyze_all_scales(o)
                tf.append({**meta,"raw_maxima":{k:float(v[0]) for k,v in raw.items()}})
            save(dest/"target_features.json",tf)
            summaries=[]
            for arm in ARMS:
                for alpha in [.01,.005]:
                    name=f"{arm}_{alpha}";cal=calibrations[name];rows=[]
                    available=bool(cal.target_fap_attainable and cal.global_threshold is not None)
                    for o,meta in target:
                        result=det.analyze(o,cal);candidates=[asdict(c) for c in result.candidates]
                        hit=bool(meta["support"] and any(associated(c,meta["support"]) for c in candidates))
                        rows.append({**meta,"score":result.max_global_score,"threshold":cal.global_threshold,"p":cal.empirical_p_value(result.max_global_score)[0] if available else None,"available":available,"detected":result.detected,"associated":hit,"candidates":candidates})
                    negative=[r for r in rows if r["part"]=="noise"]
                    groups=[]
                    for family in ["chirp","burst"]:
                        for snr in [8.,18.]:
                            group=[r for r in rows if r["family"]==family and r["snr"]==snr]
                            groups.append({"family":family,"snr":snr,"total":len(group),"detected":sum(r["detected"] for r in group),"associated":sum(r["associated"] for r in group),"abstained":sum(not r["available"] for r in group)})
                    summary={"arm":arm,"alpha":alpha,"threshold":cal.global_threshold,"noise_total":200,"available_noise_windows":sum(r["available"] for r in negative),"available_noise_seconds":4.*sum(r["available"] for r in negative),"noise_alarms":sum(r["detected"] for r in negative),"groups":groups}
                    save(dest/f"evaluation_{name}.json",rows);summaries.append(summary)
            assert all(sha((dest/f).read_bytes())==h for f,h in frozen.items())
            entry={"seed":seed,"case":case,"summaries":summaries,"threshold_hashes_unchanged":True}
            save(dest/"summary.json",entry);m["cases"].append(entry);m["parity_checks"]+=det.parity
            save(OUT/"manifest.json",m)
            print("CASE_COMPLETE",seed,case,[(s["arm"],s["noise_alarms"],sum(g["associated"] for g in s["groups"])) for s in summaries if s["alpha"]==.01],flush=True)
    assert all(sha((ROOT/f).read_bytes())==h for f,h in sources.items())
    assert all(sha((ROOT/f).read_bytes())==h for f,h in m["dependency_sha256"].items())
    m.update(status="COMPLETE",source_hashes_unchanged=True,elapsed_s=time.monotonic()-start)
    save(OUT/"manifest.json",m);print("RUN_STATUS COMPLETE",flush=True)
if __name__=="__main__":
    try:main()
    except Exception:
        if "OUT" in globals():save(OUT/"failure.json",{"status":"FAILED","traceback":traceback.format_exc()})
        raise
