"""Replay historical synthetic recipes unchanged; record current decisions and exposure."""
from __future__ import annotations
import ast
import csv
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any
import uuid

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.dont_write_bytecode=True
def guard(event,args):
    if event=="open" and args and isinstance(args[0],(str,bytes,os.PathLike)):
        text=os.path.abspath(os.fsdecode(args[0])).replace("\\","/").lower()
        if "o3_holdout" in text or text.startswith(str(ROOT/"data").replace("\\","/").lower()+"/"):
            raise PermissionError("Reserved/data-directory access blocked by synthetic replay")
        mode=args[1] if len(args)>1 else None
        flags=args[2] if len(args)>2 else 0
        writing=(isinstance(mode,str) and any(c in mode for c in "wax+")) or (isinstance(flags,int) and bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
        if writing and any(text.startswith(str(ROOT/"results"/d).replace("\\","/").lower()+"/") for d in ("v03_phase2","v03_phase2b","v03_phase2b_corrected")):
            raise PermissionError("Historical output mutation blocked")
    if event in ("socket.connect","socket.getaddrinfo"):
        raise PermissionError("Network blocked by synthetic replay")
sys.addaudithook(guard)

import numpy as np
import scipy
from scipy import signal as sp_signal
from bsd.adaptive_regime import CRED_CALIBRATED,CRED_RECALIBRATING,CRED_TRANSITION,SIG_UNCERTAIN
from bsd.adaptive_regime_v2b import (GLITCH_ELEVATED,AdaptiveRegimeManagerV2b,AdaptiveRegimeV2bConfig,
    FrozenReferenceBroadbandChannel,calibrate_separate_search_modes,fit_glitch_reference)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation,_scale_to_energy_snr,generate_clean_signal,make_split_rng
from bsd.noise import generate_noise
from bsd.workstation_runner import _build_rows

BASELINE="043f0eb95be8eeae65b9a5f5bb63015d93559b73"
SCRIPT=ROOT/"experiments/run_v03_phase2b_corrected_benchmark.py"
PROTOCOL=ROOT/"docs/historical_adaptive_replay_protocol_2026-11-26.md"
ALLOW={"_gen_obs_list","_make_chirp_chunks","_wilson_ci90","_run_suite_a_corrected","_run_suite_b_corrected"}
A_NAMES=["stationary_null_200win","abrupt_step_500win","transient_chirps_stationary","glitch_storm_120win"]
B_NAMES=[f"{bg}__{m}" for bg in ("stationary","drifting_1.0_to_1.4x")
         for m in ("80s_sinusoid","80s_multitone_comb","40s_slow_chirp","80s_broadband_stochastic")]
PLAN=A_NAMES+B_NAMES

def git(*a):
    return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
def sha(b):
    return hashlib.sha256(b).hexdigest()
def hashes(paths):
    return {str(p.relative_to(ROOT)).replace("\\","/"):sha(p.read_bytes()) for p in paths}
def norm(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)
def save(p,x):
    p.write_text(json.dumps(x,indent=2,default=norm,allow_nan=False),encoding="utf-8")
def csv_save(p,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def array_hash(obs):
    arr=np.stack([o.signal for o in obs]).astype("<f8",copy=False)
    return {"sha256":sha(arr.tobytes(order="C")),"shape":list(arr.shape),"dtype":"<f8","encoding":"little-endian C-order raw bytes"}
def ast_functions(text):
    tree=ast.parse(text)
    return ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ALLOW],type_ignores=[])

def main():
    if git("rev-parse","HEAD")!=BASELINE or git("diff",BASELINE,"--","src","app","tests","pyproject.toml"):
        raise SystemExit("Unexpected source state; write a new protocol before replay")
    historical_source=git("show","edd1ec9:experiments/run_v03_phase2b_corrected_benchmark.py")
    selected=ast_functions(SCRIPT.read_text(encoding="utf-8"))
    assert {n.name for n in selected.body}==ALLOW
    assert ast.dump(selected)==ast.dump(ast_functions(historical_source))
    assert not git("diff","edd1ec9",BASELINE,"--","src/bsd/generator.py","src/bsd/noise.py")
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    out=ROOT/"results/historical_adaptive_replay"/run_id
    out.mkdir(parents=True,exist_ok=False)
    (out/"protocol.md").write_bytes(PROTOCOL.read_bytes())
    (out/"driver.py").write_bytes(Path(__file__).read_bytes())
    ast_text=ast.unparse(selected)
    (out/"replayed_synthetic_functions.py").write_text(ast_text,encoding="utf-8")
    sources=sorted((ROOT/"src/bsd").glob("*.py"))+[SCRIPT,ROOT/"pyproject.toml"]
    preserved=[ROOT/"results/v03_phase2b_corrected"/name for name in ("results.json","report.md","suite_a_corrected.csv","suite_b_corrected.csv")]
    original_hashes=hashes(preserved)
    manifest={"run_id":run_id,"status":"RUNNING","source_commit":BASELINE,"historical_recipe_commit":"edd1ec9",
        "recipe_ast_sha256":sha(ast.dump(selected).encode()),"source_sha256":hashes(sources),
        "historical_artifact_sha256_before":original_hashes,"driver_sha256":sha(Path(__file__).read_bytes()),
        "protocol_sha256":sha(PROTOCOL.read_bytes()),"dirty_paths":git("status","--porcelain").splitlines(),
        "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
        "o3_and_data_directory_access":"BLOCKED","network":"BLOCKED","cases":[]}
    save(out/"manifest.json",manifest)
    print("RUN_DIRECTORY",out,flush=True)
    call_index=0

    class RecordingManager(AdaptiveRegimeManagerV2b):
        def process_sequence(self,observations,*args,**kwargs):
            nonlocal call_index
            name=PLAN[call_index];call_index+=1
            dest=out/name;dest.mkdir(exist_ok=False)
            started=time.monotonic()
            if name==A_NAMES[1]: sig=[108,125,145]+list(range(330,430,10))
            elif name==A_NAMES[2]:sig=list(range(15,115,5))
            elif name in B_NAMES:sig=list(range(50,60 if "40s_slow_chirp" in name else 70))
            else:sig=[]
            labels={"signal_indices":sig,"post_step_injection_indices":list(range(330,430,10)) if name==A_NAMES[1] else [],
                "step_index":100 if name==A_NAMES[1] else None,
                "impulsive_noise_start_index":60 if name==A_NAMES[3] else None,
                "manager_config":asdict(self.config),"input_array":array_hash(observations)}
            save(dest/"input_and_labels.json",labels)
            print("START",name,flush=True)
            try:
                windows=super().process_sequence(observations,*args,**kwargs)
                candidates,timeline=_build_rows(windows,4.,self.config.policy)
                csv_save(dest/"timeline.csv",timeline)
                csv_save(dest/"candidates.csv",candidates)
                save(dest/"windows.json",[asdict(w) for w in windows])
                save(dest/"outcome.json",{"status":"COMPLETE","n_windows":len(windows),"elapsed_s":time.monotonic()-started})
                manifest["cases"].append({"case":name,"status":"COMPLETE","n_windows":len(windows),"elapsed_s":time.monotonic()-started})
                save(out/"manifest.json",manifest)
                print("DONE",name,f"{time.monotonic()-started:.1f}s",flush=True)
                return windows
            except Exception as exc:
                save(dest/"outcome.json",{"status":"ERROR","type":type(exc).__name__,"error":str(exc)})
                raise

    def record_calibration(det,bb,reference,*args,**kwargs):
        scal=calibrate_separate_search_modes(det,bb,reference,*args,**kwargs)
        save(out/"reference_calibration.json",{"input_array":array_hash(reference),"kwargs":kwargs,
            "transient_null_scores":scal.transient_stft_cal.sorted_global_null_scores,
            "transient_threshold":scal.transient_stft_cal.global_threshold,
            "union_null_scores":scal.union_stft_cal.sorted_global_null_scores,
            "union_threshold":scal.union_stft_cal.global_threshold})
        return scal

    namespace={k:v for k,v in globals().items() if not k.startswith("__")}
    namespace["AdaptiveRegimeManagerV2b"]=RecordingManager
    namespace["calibrate_separate_search_modes"]=record_calibration
    exec(compile(selected,str(SCRIPT),"exec"),namespace)
    try:
        sa=namespace["_run_suite_a_corrected"]()
        save(out/"suite_a_legacy_aggregation_current_outputs.json",sa["rows"])
        sb=namespace["_run_suite_b_corrected"](sa["shared_ref"])
        save(out/"suite_b_legacy_aggregation_current_outputs.json",sb["rows"])
        assert call_index==12
        manifest["status"]="COMPLETE"
    except Exception as exc:
        manifest["status"]="ERROR"
        manifest["error"]={"type":type(exc).__name__,"message":str(exc)}
        raise
    finally:
        manifest["source_hashes_unchanged"]=hashes(sources)==manifest["source_sha256"]
        manifest["historical_hashes_unchanged"]=hashes(preserved)==original_hashes
        manifest["finished_utc"]=datetime.now(timezone.utc).isoformat()
        save(out/"manifest.json",manifest)
    assert manifest["source_hashes_unchanged"] and manifest["historical_hashes_unchanged"]
    print("RUN_STATUS",manifest["status"],flush=True)

if __name__=="__main__":
    main()
