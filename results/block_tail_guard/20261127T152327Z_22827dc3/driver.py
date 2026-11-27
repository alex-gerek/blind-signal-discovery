"""Prespecified fixed-schedule block-admission falsification pilot. No manager integration."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,importlib.util,json,platform,sys,time,traceback,uuid
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"));sys.path.insert(0,str(ROOT/"experiments"));sys.dont_write_bytecode=True
# Reuse the existing synthetic-only IO/network guard, without invoking its historical main.
spec=importlib.util.spec_from_file_location("synthetic_guard",ROOT/"experiments/run_historical_adaptive_replay.py")
replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation,generate_clean_signal,_scale_to_energy_snr
from bsd.noise import generate_pink_colored
from block_tail_guard import block_guard,focused_checks
PROTOCOL=ROOT/"docs/block_tail_guard_protocol_2026-11-27.md"
CASES=["pink","calibration_chirps","calibration_bursts","impulsive_noise"]
def sha(x):return hashlib.sha256(x).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,default=replay.norm,allow_nan=False),encoding="utf-8")
def key(o):return sha(np.asarray(o.signal,dtype="<f8").tobytes()+np.asarray([o.sample_rate,o.duration],dtype="<f8").tobytes())
def clean(family,rng):
    if family=="chirp":return generate_clean_signal("linear_chirp",4096,1024.,rng)[0]
    x=signal.sosfiltfilt(signal.butter(4,[35.,350.],btype="bandpass",fs=1024.,output="sos"),rng.normal(size=4096))
    a,b=int(1.375*1024),int(1.625*1024);mask=np.zeros(4096);mask[a:b]=signal.windows.tukey(b-a,.15)
    return x*mask

def generate(seed):
    rc=np.random.default_rng(np.random.SeedSequence([seed,1]));rn=np.random.default_rng(np.random.SeedSequence([seed,2]))
    rs=np.random.default_rng(np.random.SeedSequence([seed,3]));ri=np.random.default_rng(np.random.SeedSequence([seed,4]))
    def noise(rng):return generate_pink_colored(4096,1024.,rng,sigma=1.,alpha=1.,f_knee=16.)[0]
    cal=[noise(rc) for _ in range(265)];neg=[noise(rn) for _ in range(80)]
    injected={f:{i:clean(f,rs) for i in [100,140,180,220]} for f in ["chirp","burst"]}
    units=[(f,k,noise(rn),clean(f,rs)) for f in ["chirp","burst"] for k in range(3)]
    def impulsive(x):
        y=x.copy()
        if ri.random()<.03:y[int(ri.integers(4096))]+=float(ri.choice([-12.,12.]))
        return y
    heavy_cal=[impulsive(x) for x in cal];heavy_neg=[impulsive(x) for x in neg]
    heavy_units=[(f,k,impulsive(n),s) for f,k,n,s in units]
    for case in CASES:
        cs=[x.copy() for x in (heavy_cal if case=="impulsive_noise" else cal)]
        family={"calibration_chirps":"chirp","calibration_bursts":"burst"}.get(case)
        if family:
            for i,s in injected[family].items():cs[i]+=_scale_to_energy_snr(s,1.,18.)
        records=[]
        def add(part,index,x,family=None,snr=None,unit=None):
            ident=f"s{seed}_{case}_{part}_{index:03d}"
            o=Observation(ident,x,1024.,4.)
            records.append((o,{"sample_id":ident,"seed":seed,"case":case,"part":part,"index":index,"family":family,"snr":snr,"unit":unit,"input_sha256":key(o)}))
        for i,x in enumerate(cs):add("calibration",i,x,family if i in [100,140,180,220] else None,18. if family and i in [100,140,180,220] else None)
        for i,x in enumerate(heavy_neg if case=="impulsive_noise" else neg):add("noise",i,x)
        for j,(f,k,n,s) in enumerate(heavy_units if case=="impulsive_noise" else units):
            for ai,snr in enumerate([8.,18.]):add("signal",2*j+ai,n+_scale_to_energy_snr(s,1.,snr),f,snr,k)
        yield case,records

class CachedCore(MultiResolutionDetector):
    def __init__(self):super().__init__();self.cache={};self.parity=0
    def _raw_analyze_all_scales(self,o):
        h=key(o)
        if h not in self.cache:self.cache[h]=super()._raw_analyze_all_scales(o)
        return self.cache[h]
    def parity_check(self,o):
        cached=self._raw_analyze_all_scales(o);direct=super()._raw_analyze_all_scales(o)
        assert json.dumps(cached,default=lambda v:asdict(v) if hasattr(v,"__dataclass_fields__") else replay.norm(v),sort_keys=True)==json.dumps(direct,default=lambda v:asdict(v) if hasattr(v,"__dataclass_fields__") else replay.norm(v),sort_keys=True)
        self.parity+=1

def main():
    assert replay.git("rev-parse","HEAD")=="d049e2b1c97626b562d234a22bc12db5d51702dd"
    assert not replay.git("diff","--","src","app","pyproject.toml")
    ident=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    global OUT
    OUT=ROOT/"results/block_tail_guard"/ident;OUT.mkdir(parents=True,exist_ok=False)
    print("RUN_DIRECTORY",OUT,flush=True)
    (OUT/"protocol.md").write_bytes(PROTOCOL.read_bytes())
    (OUT/"driver.py").write_bytes(Path(__file__).read_bytes())
    (OUT/"guard.py").write_bytes((ROOT/"experiments/block_tail_guard.py").read_bytes())
    sources=replay.hashes(list((ROOT/"src/bsd").glob("*.py")))
    manifest={"run_id":ident,"status":"RUNNING","source_commit":replay.git("rev-parse","HEAD"),"dirty_paths":replay.git("status","--short"),"source_sha256":sources,"python":platform.python_version(),"numpy":np.__version__,"scipy":replay.scipy.__version__,"seeds":[17101,17102,17103],"scope":"Fixed-schedule exploratory pilot, not end-to-end adaptive recovery","cases":[]}
    save(OUT/"manifest.json",manifest);save(OUT/"focused_checks.json",focused_checks())
    det=CachedCore();start=time.monotonic()
    for seed in manifest["seeds"]:
        for case,records in generate(seed):
            dest=OUT/f"{seed}_{case}";dest.mkdir(exist_ok=False)
            cal_obs=[o for o,m in records if m["part"]=="calibration"]
            evaluations=[(o,m) for o,m in records if m["part"]!="calibration"]
            for o in cal_obs[:3]+[o for o,m in evaluations[:3]]+([cal_obs[100]] if case.startswith("calibration_") else []):det.parity_check(o)
            features=[]
            for o,m in records:
                raw=det._raw_analyze_all_scales(o)
                features.append({**m,"raw_maxima":{k:float(v[0]) for k,v in raw.items()}})
            save(dest/"input_features.json",features)
            summaries=[]
            for alpha in [.01,.005]:
                cal=det.calibrate(cal_obs,target_fap=alpha,norm_fraction=.25,calibration_seed=seed,allow_unattainable_fap=True)
                guard=block_guard(cal.sorted_global_null_scores)
                available=guard["admitted"] and cal.target_fap_attainable and cal.global_threshold is not None
                save(dest/f"calibration_{alpha}.json",asdict(cal));rows=[]
                for o,m in evaluations:
                    result=det.analyze(o,cal)
                    p=cal.empirical_p_value(result.max_global_score)[0]
                    diagnostic=[asdict(c) for c in result.candidates]
                    rows.append({**m,"core_score":result.max_global_score,"core_threshold":cal.global_threshold,"core_p":p,"core_detected":result.detected,
                        "diagnostic_candidates":diagnostic,"policy_available":bool(available),"policy_threshold":cal.global_threshold if available else None,
                        "policy_p":p if available else None,"policy_reason":guard["reason"] if not available else "AVAILABLE",
                        "policy_detected":bool(available and result.detected)})
                negatives=[r for r in rows if r["part"]=="noise"];positives=[r for r in rows if r["part"]=="signal"]
                summary={"alpha":alpha,"guard":guard,"policy_available":bool(available),"core_noise_alarms":sum(r["core_detected"] for r in negatives),
                    "policy_noise_alarms":sum(r["policy_detected"] for r in negatives),"available_noise_windows":80 if available else 0,"available_noise_seconds":320. if available else 0.,
                    "conditional_noise_rate":sum(r["policy_detected"] for r in negatives)/80 if available else None,
                    "core_signal_detected":sum(r["core_detected"] for r in positives),"policy_signal_detected":sum(r["policy_detected"] for r in positives),
                    "policy_signal_missed":sum(not r["policy_detected"] for r in positives) if available else 0,"policy_signal_abstained":0 if available else 12,
                    "lost_baseline_detections":sum(r["core_detected"] and not r["policy_detected"] for r in positives)}
                save(dest/f"evaluation_{alpha}.json",rows);summaries.append(summary)
            entry={"seed":seed,"case":case,"summaries":summaries};save(dest/"summary.json",entry)
            manifest["cases"].append(entry);save(OUT/"manifest.json",manifest)
            print("CASE_COMPLETE",seed,case,"admitted",summaries[0]["guard"]["admitted"],"lost",summaries[0]["lost_baseline_detections"],flush=True)
    assert all(replay.sha((ROOT/f).read_bytes())==h for f,h in sources.items())
    manifest.update(status="COMPLETE",source_hashes_unchanged=True,elapsed_s=time.monotonic()-start,cache_uncached_parity_checks=det.parity,unique_cached_observations=len(det.cache))
    save(OUT/"manifest.json",manifest);print("RUN_STATUS COMPLETE",flush=True)
if __name__=="__main__":
    try:main()
    except Exception:
        if "OUT" in globals():save(OUT/"failure.json",{"status":"FAILED","traceback":traceback.format_exc()})
        raise
