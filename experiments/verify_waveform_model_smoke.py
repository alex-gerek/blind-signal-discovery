"""Independent numerical reconstruction of waveform-smoke evidence."""
import os
os.environ.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
import sys,json,hashlib,gzip,importlib.util
from pathlib import Path
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
PARENT=ROOT/"results/matched_gap2/20261128T115035Z_bf957291"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(x):return hashlib.sha256(np.asarray(x,dtype="<f8").tobytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n",encoding="utf-8")
def main(out):
    rows=json.loads((out/"results.json").read_text(encoding="utf-8"))
    with np.load(out/"background_weights.npz") as z:weights=z["weights"].copy(); savedpsd=z["mean_welch_psd"].copy()
    spec=importlib.util.spec_from_file_location("frozen_generator",ROOT/"experiments/run_blind_two_tone_exploratory.py");old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    backgrounds=json.loads((out/"background_hashes.json").read_text(encoding="utf-8"));psds=[]
    refs=[r for r in json.loads(gzip.decompress((PARENT/"calibration_raw.json.gz").read_bytes())) if r["w"]==1]
    for row,ref in zip(backgrounds,refs):
        x,_=old.make(40101,row["index"],"noise");assert key(x)==row["hash"]==ref["hash"]
        f,p=signal.welch(x,fs=1024,window="hann",nperseg=1024,noverlap=512,detrend="constant",scaling="density");psds.append(p)
    assert len(backgrounds)==len(refs)==284
    psd=np.mean(psds,axis=0);np.testing.assert_array_equal(psd,savedpsd)
    med=np.median(psd[psd>0]);expected=np.sqrt(med/np.maximum(np.interp(np.fft.rfftfreq(4096,1/1024),f,psd),1e-6*med));np.testing.assert_array_equal(expected,weights)
    env=np.zeros(4096);env[1280:2304]=signal.windows.tukey(1024,.15);t=np.arange(4096)/1024
    def whiten(a):return np.fft.irfft(np.fft.rfft(a-np.mean(a,axis=0),axis=0)*(weights[:,None] if a.ndim==2 else weights),n=4096,axis=0)
    def design(fs):return np.column_stack([v for f in fs for v in [env*np.cos(2*np.pi*f*t),env*np.sin(2*np.pi*f*t)]])
    count=0;normal=[]
    for row in rows:
        with np.load(out/(row["id"]+".npz")) as z:
            x=z["waveform"].copy();assert key(x)==row["input_hash"]
            y=whiten(x);np.testing.assert_allclose(y@y,row["Q0"],rtol=1e-11)
            for name in ["M1","M2"]:
                b=design(row[name+"_frequencies"]);a=whiten(b);beta=np.asarray(row[name+"_coefficients_cos_sin"]);pred=b@beta+row[name+"_offset"]
                np.testing.assert_allclose(pred,z[name+"_prediction"],rtol=1e-9,atol=1e-10)
                np.testing.assert_allclose(x-pred,z[name+"_raw_residual"],rtol=1e-9,atol=1e-10)
                r=whiten(x-pred);np.testing.assert_allclose(r@r,row["Q"+name[-1]],rtol=1e-9)
                relative=float(np.linalg.norm(a.T@r)/(np.linalg.norm(a)*np.linalg.norm(y)))
                assert relative<1e-9;normal.append(relative)
                np.testing.assert_allclose(np.mean(x-pred),0,atol=1e-10)
            for start in row["starts"]:
                a=whiten(design(start["frequencies"]));beta=np.linalg.lstsq(a,y,rcond=None)[0];r=y-a@beta
                np.testing.assert_allclose(r@r,start["Q"],rtol=1e-9);count+=1
            assert row["nested_loss_check"] and row["Q2"]<=row["Q1"]+1e-8
            if row["kind"]=="two":
                errors=np.abs(np.sort(row["M2_frequencies"])-np.sort(row["truth_frequencies"]))
                np.testing.assert_allclose(errors,row["evaluation_errors_hz"])
                assert bool(np.all(errors<=.5))==row["two_frequency_match_0_5_hz"]
            else:assert row["two_frequency_match_0_5_hz"] is None
            assert row["tone_count_decision"] is None
    parent=json.loads((PARENT/"artifact_sha256.json").read_text(encoding="utf-8"))["files"]
    for p,h in parent.items():assert sha(ROOT/p)==h,p
    manifest=json.loads((out/"manifest.json").read_text(encoding="utf-8"))
    for p,h in manifest["production_hashes"].items():assert sha(ROOT/p)==h
    report={"status":"PASS","records":len(rows),"background_hash_and_psd_checks":len(backgrounds),"model_reconstructions":2*len(rows),"optimizer_returned_objectives":count,"max_relative_normal_equation_residual":max(normal),"parent_artifacts_unchanged":len(parent),"production_source_files_unchanged":len(manifest["production_hashes"]),"scope":"Numerical implementation checks only; no validation of tone-count selection, global optimization, unknown-envelope robustness or scientific correctness."}
    save(out/"verification.json",report);print(json.dumps(report))
if __name__=="__main__":main(Path(sys.argv[1]))
