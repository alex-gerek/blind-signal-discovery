"""Independent arithmetic/accounting verification of saved whole-block guard pilot."""
from pathlib import Path
import hashlib,importlib.util,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def load(p):
    return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding="utf-8")
def main():
    out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/"results/block_tail_guard"
    m=load(out/"manifest.json");assert m["status"]=="COMPLETE" and len(m["cases"])==12
    report={"status":"RUNNING","input_records":0,"calibrations":0,"evaluation_records":0,"candidate_records":0,"groups":{},"scope":"Descriptive fixed-schedule pilot, three paired seed units, no adaptive recovery inference"}
    for case in m["cases"]:
        dest=out/f"{case['seed']}_{case['case']}"
        features=load(dest/"input_features.json");by_id={r["sample_id"]:r for r in features}
        assert len(features)==len(by_id)==357
        calibration=[r for r in features if r["part"]=="calibration"]
        assert len(calibration)==265
        assert [r["index"] for r in calibration]==list(range(265))
        assert sum(r["family"] is not None for r in calibration)==(4 if case["case"].startswith("calibration_") else 0)
        assert all(r["family"] is None for r in calibration[:66])
        for alpha,summary in zip([.01,.005],case["summaries"]):
            cal=load(dest/f"calibration_{alpha}.json");rows=load(dest/f"evaluation_{alpha}.json")
            assert cal["n_norm_samples"]==66 and cal["n_threshold_samples"]==199 and len(rows)==92
            stats={}
            for name in cal["scale_stats"]:
                a=np.array([r["raw_maxima"][name] for r in calibration[:66]])
                med=float(np.median(a));q=float(np.quantile(a,.95));spread=max(q-med,1.)
                s=cal["scale_stats"][name]
                assert s["median_null_score"]==med and s["q95_null_score"]==q and s["spread_scale"]==spread
                stats[name]=(med,spread)
            def maximum(r):return max((r["raw_maxima"][name]-med)/spread for name,(med,spread) in stats.items())
            null=sorted(maximum(r) for r in calibration[66:])
            assert null==cal["sorted_global_null_scores"]
            # Enumerate allowable tail counts, independently from production rank helper.
            tail=max(k for k in range(200) if k/200.<=alpha)
            assert tail>=1
            threshold=null[200-tail-1];assert threshold==cal["global_threshold"]
            med=float(np.median(null));spread=max(float(np.quantile(null,.95))-med,1.);limit=med+10*spread
            admitted=max(null)<=limit
            assert summary["guard"]=={"admitted":admitted,"reason":"ADMITTED" if admitted else "WHOLE_BLOCK_TAIL_QUARANTINE","median":med,"spread":spread,"limit":limit,"maximum":max(null)}
            assert summary["policy_available"]==admitted
            for r in rows:
                source=by_id[r["sample_id"]]
                assert all(r[k]==source[k] for k in ["seed","case","part","index","family","snr","unit","input_sha256"])
                score=maximum(source);p=(1+sum(x>=score for x in null))/200.
                assert r["core_score"]==score and r["core_p"]==p and r["core_threshold"]==threshold
                detected=score>threshold and p<=alpha
                assert r["core_detected"]==detected==bool(r["diagnostic_candidates"])
                assert r["policy_detected"]==bool(admitted and detected) and r["policy_available"]==admitted
                assert r["policy_p"]==(p if admitted else None)
                assert r["policy_threshold"]==(threshold if admitted else None)
                assert r["policy_reason"]==("AVAILABLE" if admitted else "WHOLE_BLOCK_TAIL_QUARANTINE")
                for candidate in r["diagnostic_candidates"]:
                    cp=(1+sum(x>=candidate["standardized_score"] for x in null))/200.
                    assert candidate["standardized_score"]>threshold and candidate["empirical_p_value"]==cp and cp<=alpha
                    report["candidate_records"]+=1
                report["evaluation_records"]+=1
            negative=[r for r in rows if r["part"]=="noise"];positive=[r for r in rows if r["part"]=="signal"]
            assert len(negative)==len({r["sample_id"] for r in negative})==80 and len(positive)==12
            assert len({(r["family"],r["unit"]) for r in positive})==6
            assert summary["core_noise_alarms"]==sum(r["core_detected"] for r in negative)
            assert summary["policy_noise_alarms"]==sum(r["policy_detected"] for r in negative)
            assert summary["available_noise_windows"]==(80 if admitted else 0)
            assert summary["available_noise_seconds"]==(320. if admitted else 0.)
            assert summary["conditional_noise_rate"]==(sum(r["policy_detected"] for r in negative)/80 if admitted else None)
            assert summary["core_signal_detected"]==sum(r["core_detected"] for r in positive)
            assert summary["policy_signal_detected"]==sum(r["policy_detected"] for r in positive)
            assert summary["policy_signal_abstained"]==(0 if admitted else 12)
            assert summary["policy_signal_missed"]==(sum(not r["policy_detected"] for r in positive) if admitted else 0)
            assert summary["lost_baseline_detections"]==sum(r["core_detected"] and not r["policy_detected"] for r in positive)
            k=f"{case['case']}_fap{alpha}"
            group=report["groups"].setdefault(k,{"blocks":0,"admitted":0,"core_noise_alarms":0,"policy_noise_alarms":0,"available_noise_windows":0,"core_signal_detected":0,"policy_signal_detected":0,"policy_signal_abstained":0,"lost_baseline_detections":0})
            group["blocks"]+=1;group["admitted"]+=int(admitted)
            for key in list(group):
                if key not in ["blocks","admitted"]:group[key]+=summary[key]
            report["calibrations"]+=1
        report["input_records"]+=len(features)
    # Regenerate all inputs from the frozen generator; this checks preservation, not generator validity.
    spec=importlib.util.spec_from_file_location("guard_pilot",ROOT/"experiments/run_block_tail_guard.py")
    pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)
    assert (out/"driver.py").read_bytes()==Path(pilot.__file__).read_bytes()
    regenerated=0
    for seed in m["seeds"]:
        for case,records in pilot.generate(seed):
            saved=load(out/f"{seed}_{case}/input_features.json")
            for (o,meta),r in zip(records,saved):
                assert all(r[k]==v for k,v in meta.items()) and pilot.key(o)==r["input_sha256"]
                regenerated+=1
    assert regenerated==report["input_records"]==4284
    for name,h in m["source_sha256"].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    report["regenerated_inputs"]=regenerated;report["source_hashes_unchanged"]=True
    report["focused_fixture_checks"]=load(out/"focused_checks.json")["count"]
    report["uncached_parity_checks"]=m["cache_uncached_parity_checks"]
    report["criteria"]={
        "clean_pink_all_admitted":all(v["admitted"]==v["blocks"] for k,v in report["groups"].items() if k.startswith("pink_")),
        "contaminated_all_quarantined":all(v["admitted"]==0 for k,v in report["groups"].items() if k.startswith("calibration_")),
        "no_baseline_recovery_lost":all(v["lost_baseline_detections"]==0 for v in report["groups"].values())}
    report["adoption_rejected"]=not all(report["criteria"].values());report["status"]="PASS"
    save(out/"independent_verification.json",report);print(json.dumps(report,indent=2))
if __name__=="__main__":main()
