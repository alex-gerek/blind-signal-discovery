"""Known-envelope waveform smoke study; no production imports are modified."""
import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import sys, json, gzip, hashlib, time, uuid, platform, subprocess, importlib.util
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import scipy
from scipy import signal, optimize
ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
PARENT = ROOT / "results/matched_gap2/20261128T115035Z_bf957291"
PROTOCOL = ROOT / "docs/waveform_model_smoke_protocol_2026-11-28.md"
IDS = ["b00_two_0_0_1", "b00_two_0_0_0.5", "b00_two_0.5_0_1", "b05_two_0_0_0.25", "b04_two_0.5_0_0.5", "b09_one_0_0_1", "b00_noise_0_0_1", "b09_noise_0_0_1"]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def key(a): return hashlib.sha256(np.asarray(a, dtype="<f8").tobytes()).hexdigest()
def save(p, a): p.write_text(json.dumps(a, indent=2, allow_nan=False) + "\n", encoding="utf-8")
def readgz(p): return json.loads(gzip.decompress(p.read_bytes()))
def git(*a): return subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()

class Fitter:
    def __init__(self, weights):
        self.weights = weights
        self.t = np.arange(4096) / 1024.
        self.env = np.zeros(4096)
        self.env[1280:2304] = signal.windows.tukey(1024, .15)
        self.grid = np.arange(20., 450.01, .5)
        self.calls = 0
        self.rawgrid = self.basis(self.grid)
        self.wgrid = self.whiten(self.rawgrid)
    def whiten(self, a):
        return np.fft.irfft(np.fft.rfft(a - np.mean(a, axis=0), axis=0) * (self.weights[:, None] if a.ndim == 2 else self.weights), n=4096, axis=0)
    def basis(self, fs):
        angle = 2 * np.pi * self.t[:, None] * np.asarray(fs)[None, :]
        return np.stack([self.env[:, None] * np.cos(angle), self.env[:, None] * np.sin(angle)], axis=2).reshape(4096, -1)
    def fit(self, x, fs):
        self.calls += 1
        b = self.basis(fs)
        a = self.whiten(b)
        y = self.whiten(x)
        beta = np.linalg.lstsq(a, y, rcond=None)[0]
        offset = float(np.mean(x - b @ beta))
        prediction = b @ beta + offset
        residual = y - a @ beta
        return residual, beta, offset, prediction, float(np.linalg.cond(a))
    def grid_gain(self, y, fixed=None):
        c = self.wgrid[:, 0::2].copy()
        s = self.wgrid[:, 1::2].copy()
        if fixed is not None:
            a = self.whiten(self.basis([fixed]))
            q = np.linalg.qr(a)[0]
            y = y - q @ (q.T @ y)
            c -= q @ (q.T @ c)
            s -= q @ (q.T @ s)
        cc = np.sum(c*c, axis=0); ss = np.sum(s*s, axis=0); cs = np.sum(c*s, axis=0)
        yc = y @ c; ys = y @ s
        det = cc*ss-cs*cs
        gains = np.full(len(cc), -np.inf)
        ok = det > 1e-12 * np.maximum(cc*ss, 1e-30)
        gains[ok] = (ss[ok]*yc[ok]**2 + cc[ok]*ys[ok]**2 - 2*cs[ok]*yc[ok]*ys[ok])/det[ok]
        return gains
    def maxima(self, gains, n):
        ix = set(signal.find_peaks(gains)[0].tolist()) | {0, len(gains)-1, int(np.argmax(gains))}
        return sorted(ix, key=lambda i: (-gains[i], i))[:n]
    def solve(self, x):
        started = time.perf_counter(); count0 = self.calls
        y = self.whiten(x); q0 = float(y@y)
        gains = self.grid_gain(y); tops = self.maxima(gains, 5)
        m1_candidates = [(float(self.grid[i]), "grid", None) for i in {0, len(self.grid)-1, int(np.argmax(gains))}]
        for i in tops:
            res = optimize.minimize_scalar(lambda f: float(self.fit(x, [f])[0] @ self.fit(x, [f])[0]), bounds=(max(20, self.grid[i]-.5), min(450, self.grid[i]+.5)), method="bounded", options={"maxiter":80,"xatol":1e-7})
            m1_candidates.append((float(res.x), "refined", {"success":bool(res.success),"nfev":int(res.nfev)}))
        f1, _, _ = min(m1_candidates, key=lambda row: float(np.sum(self.fit(x, [row[0]])[0]**2)))
        r1, b1, off1, pred1, cond1 = self.fit(x, [f1]); q1 = float(r1@r1)
        m1_seconds = time.perf_counter()-started
        s2 = time.perf_counter()
        partners = self.maxima(self.grid_gain(y, f1), 6)
        seeds = [(f1, float(self.grid[i])) for i in partners]
        seeds += [(f1-d/2, f1+d/2) for d in [.25,.5,1,2,4,8,16,32,64]]
        from itertools import combinations
        seeds += [tuple(float(self.grid[i]) for i in pair) for pair in combinations(tops[:3],2)]
        unique=[]
        for pair in seeds:
            pair=tuple(sorted(pair))
            if 20 <= pair[0] < pair[1] <=450 and not any(np.allclose(pair,p,rtol=0,atol=1e-9) for p in unique): unique.append(pair)
        assert len(unique)<=18
        starts=[]; best=None
        # Explicit nested baseline: second sinusoid coefficient zero is feasible.
        for pair in unique:
            fit = optimize.least_squares(lambda fs: self.fit(x, fs)[0], pair, bounds=([20,20],[450,450]), jac="2-point", max_nfev=80, ftol=1e-9, xtol=1e-9, gtol=1e-7)
            fs=sorted(float(f) for f in fit.x)
            r,b,off,pred,cond=self.fit(x,fs); loss=float(r@r)
            row={"initial":list(pair),"frequencies":fs,"Q":loss,"success":bool(fit.success),"status":int(fit.status),"message":str(fit.message),"nfev":int(fit.nfev),"condition":cond,"at_band_boundary":any(abs(f-20)<1e-5 or abs(f-450)<1e-5 for f in fs)}
            starts.append(row)
            if best is None or loss<best[0]: best=(loss,fs,b,off,pred,cond,row)
        q2,fs,b2,off2,pred2,cond2,bestrow=best
        nested=q2<=q1+1e-8*max(1,q1)
        result={"Q0":q0,"Q1":q1,"Q2":q2,"gain_Q1_minus_Q2":q1-q2,"M1_frequencies":[f1],"M2_frequencies":fs,"M1_coefficients_cos_sin":b1.tolist(),"M2_coefficients_cos_sin":b2.tolist(),"M1_offset":off1,"M2_offset":off2,"M1_condition":cond1,"M2_condition":cond2,"unreliable_condition":cond2>1e8,"nested_loss_check":nested,"selected_success":bestrow["success"],"M1_candidates":[{"f":f,"origin":origin,"optimizer":opt} for f,origin,opt in m1_candidates],"starts":starts,"M1_seconds":m1_seconds,"M2_seconds":time.perf_counter()-s2,"fit_calls":self.calls-count0,"tone_count_decision":None}
        return result,pred1,pred2

def fixtures(fitter):
    checks=[]
    for fs,beta in [([71.37], [.73,-.24]),([67.13,73.81],[.7,-.2,.3,.4])]:
        x=fitter.basis(fs)@np.asarray(beta)+.31
        r,b,off,pred,cond=fitter.fit(x,fs)
        assert np.max(np.abs(pred-x))<1e-10 and abs(off-.31)<1e-10 and np.max(abs(b-beta))<1e-10
        checks.append("fixed_design_reconstruction_"+str(len(fs)))
        if len(fs)==2:
            rp=fitter.fit(x,fs[::-1])[0]
            assert np.max(abs(r-rp))<1e-10
            assert r@r <= np.sum(fitter.fit(x,fs[:1])[0]**2)+1e-9
            checks += ["permutation_invariance","nested_fixed_design_loss"]
    return {"status":"PASS","checks":checks}

def main():
    out=ROOT/"results/waveform_model_smoke"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
    out.mkdir(parents=True,exist_ok=False);print("RUN_DIRECTORY",out,flush=True)
    (out/"driver.py").write_bytes(Path(__file__).read_bytes());(out/"protocol.md").write_bytes(PROTOCOL.read_bytes())
    manifest={"status":"RUNNING","source_commit":git("rev-parse","HEAD"),"dirty":git("status","--short"),"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform(),"processor":platform.processor(),"threads":{k:os.environ[k] for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]},"production_hashes":{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/"src/bsd").glob("*.py")}}
    save(out/"manifest.json",manifest)
    start=time.perf_counter()
    spec=importlib.util.spec_from_file_location("frozen_blind",ROOT/"experiments/run_blind_two_tone_exploratory.py");old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    refs=[r for r in readgz(PARENT/"calibration_raw.json.gz") if r["w"]==1];assert len(refs)==284
    backgrounds=[];psds=[]
    for i,row in enumerate(refs):
        assert row["id"]==f"cal_{i}"
        x,_=old.make(40101,i,"noise");assert key(x)==row["hash"]
        f,p=signal.welch(x,fs=1024,window="hann",nperseg=1024,noverlap=512,detrend="constant",scaling="density")
        psds.append(p);backgrounds.append({"index":i,"hash":key(x)})
    psd=np.mean(psds,axis=0); med=float(np.median(psd[psd>0])); interp=np.interp(np.fft.rfftfreq(4096,1/1024),f,psd);weights=np.sqrt(med/np.maximum(interp,1e-6*med))
    np.savez_compressed(out/"background_weights.npz",welch_frequencies=f,mean_welch_psd=psd,weights=weights)
    save(out/"background_hashes.json",backgrounds)
    background_seconds=time.perf_counter()-start
    t=time.perf_counter();fitter=Fitter(weights);precompute_seconds=time.perf_counter()-t
    save(out/"fixtures.json",fixtures(fitter))
    evaluation={r["id"]:r for r in readgz(PARENT/"evaluation.json.gz") if r["w"]==1}
    results=[]
    for ident in IDS:
        row=evaluation[ident];path=PARENT/row["trace"]
        with np.load(path) as z:x=z["waveform"].copy()
        assert key(x)==row["meta"]["hash"]
        result,p1,p2=fitter.solve(x)
        # Evaluation labels are attached ONLY after blind fitting.
        m=row["meta"]; truth=[m[k] for k in ["f1","f2"] if m.get(k) is not None]
        errors=[abs(a-b) for a,b in zip(result["M2_frequencies"],sorted(truth))] if m["kind"]=="two" else None
        result.update(id=ident,kind=m["kind"],truth_frequencies=truth,evaluation_errors_hz=errors,two_frequency_match_0_5_hz=all(e<=.5 for e in errors) if errors is not None else None,input_hash=key(x),trace_sha256=sha(path),truth_metadata=m)
        np.savez_compressed(out/(ident+".npz"),waveform=x,M1_prediction=p1,M2_prediction=p2,M1_raw_residual=x-p1,M2_raw_residual=x-p2)
        results.append(result);save(out/"results.json",results)
        print(ident,"M2",result["M2_frequencies"],"gain",round(result["gain_Q1_minus_Q2"],4),"seconds",round(result["M1_seconds"]+result["M2_seconds"],3),flush=True)
    assert all(sha(ROOT/p)==h for p,h in manifest["production_hashes"].items())
    manifest.update(status="COMPLETE",background_seconds=background_seconds,precompute_seconds=precompute_seconds,total_seconds=time.perf_counter()-start,records=len(results),background_records=284)
    save(out/"manifest.json",manifest)
    print("COMPLETE",out,flush=True)
if __name__=="__main__":main()
