"""Independent saved-record checks; no detector execution or production imports."""
from pathlib import Path
import hashlib,json,math,sys

ROOT=Path(__file__).resolve().parents[1]
def load(p): return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
def save(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding="utf-8")
def exposure(rows):
    seen={}; intervals=[]
    for ident,a,b in rows:
        if not(math.isfinite(a) and math.isfinite(b)) or b<a: raise ValueError("Invalid exposure interval")
        if ident in seen and seen[ident]!=(a,b): raise ValueError("Conflicting duplicate ID")
        seen[ident]=(a,b)
        if b>a: intervals.append((a,b))
    end=None;total=0.
    for a,b in sorted(intervals):
        if end is None or a>=end: total+=b-a;end=b
        elif b>end: total+=b-end;end=b
    return total

def main():
    out=Path(sys.argv[1]).resolve()
    assert out.parent==ROOT/"results/adaptive_admission_trace"
    m=load(out/"manifest.json");assert m["status"]=="COMPLETE"
    assert exposure([])==0
    assert exposure([("a",0,0)])==0
    assert exposure([("a",0,4),("a",0,4),("b",2,6),("c",8,10)])==8
    for rows in [[("a",1,0)],[("a",0,float("inf"))],[("a",float("nan"),1)],[("a",0,1),("a",0,2)]]:
        try: exposure(rows)
        except ValueError: pass
        else: raise AssertionError("invalid interval accepted")
    report={"scope":"Unchanged synthetic historical reproduction, not independent efficacy validation", "interval_fixture_checks":7,"arms":{},"checks":{"windows":0,"calibrations":0,"members":0,"p_values":0}}
    later=list(range(370,430,10))
    old=ROOT/"results/historical_adaptive_followup/20261126T164211Z_dd2ea691"
    for arm in ["baseline","without_earlier_four_injections"]:
        dest=out/arm
        assert load(dest/"outcome.json")["input_array"]==load(old/arm/"outcome.json")["input_array"]
        w=load(dest/"windows.json");u=load(dest/"used_calibration_timeline.json");traces=load(dest/"calibration_trace.json")
        assert w==load(old/arm/"windows.json") and len(w)==len(u)==500
        old_traces=load(old/arm/"calibration_trace.json")
        assert len(traces)==len(old_traces)
        by_id={t["calibration_id"]:t for t in traces}
        for t,prev in zip(traces,old_traces):
            assert all(t[k]==v for k,v in prev.items())
            cal=load(dest/(t["calibration_id"]+".json"))
            n=t["n_total"];na=t["n_norm"];nb=t["n_threshold"]
            assert n==len(t["sample_ids"])==len(t["input_hashes"])==na+nb
            assert len(set(t["sample_ids"]))==n
            assert cal["n_norm_samples"]==na and cal["n_threshold_samples"]==nb
            assert cal["global_threshold"]==t["global_threshold"]
            null=cal["sorted_global_null_scores"]
            assert len(null)==nb and null==sorted(null) and all(map(math.isfinite,null))
            for entry in t["known_injected_members"]:
                k=t["sample_ids"].index(entry["sample_id"])
                assert entry["partition"]==("NORMALIZATION_A" if k<na else "THRESHOLD_B")
            report["checks"]["calibrations"]+=1;report["checks"]["members"]+=n
        switches=[]
        for i,(r,v) in enumerate(zip(u,w)):
            assert r["window_index"]==v["index"]==i and r["sample_id"]==v["sample_id"]
            assert r["score"]==v["stft_max_score"] and r["threshold"]==v["stft_threshold"]
            assert r["available"]==v["diagnostics"]["dedicated_stft_available"]==v["detection_available"]
            assert r["raw_candidate_count"]==len(v["raw_stft_candidates"])
            t=by_id[r["calibration_id"]];assert t["earliest_use_index"]<=i
            assert t["fit_window_index"] is None or t["fit_window_index"]<i
            cal=load(dest/(r["calibration_id"]+".json"))
            assert r["threshold"]==cal["global_threshold"]
            if i==0 or r["calibration_id"]!=u[i-1]["calibration_id"]:
                switches.append({"index":i,"calibration_id":r["calibration_id"],"threshold":r["threshold"]})
                if i: assert u[i-1]["updated_for_next"] and t["fit_window_index"]==i-1
            if r["available"]:
                null=cal["sorted_global_null_scores"]
                empirical=(1+sum(x>=r["score"] for x in null))/(len(null)+1)
                assert empirical==r["empirical_p"] and r["p_unavailable_reason"] is None
                assert r["raw_nominated"]==bool(r["score"]>r["threshold"] and empirical<=cal["target_fap"])
                report["checks"]["p_values"]+=1
            else:
                assert r["empirical_p"] is None and r["p_unavailable_reason"]=="CALIBRATION_UNAVAILABLE"
                assert not r["raw_nominated"]
            report["checks"]["windows"]+=1
        eligible=[r for r in u if r["available"] and not r["known_signal"]]
        unavailable=[r for r in u if not r["available"] and not r["known_signal"]]
        mature=[(r,v) for r,v in zip(u,w) if not r["known_signal"] and v["regime_id"]>=1 and v["credibility_status"]==w[0]["credibility_status"] and v["glitch_storm_status"]==w[0]["glitch_storm_status"] and r["available"] and r["mature"]]
        report["arms"][arm]={"switches":switches,"available_noise_windows":len(eligible),"unavailable_noise_windows":len(unavailable),
            "available_noise_seconds":exposure([(r["sample_id"],r["start_s"],r["end_s"]) for r in eligible]),
            "unavailable_noise_seconds":exposure([(r["sample_id"],r["start_s"],r["end_s"]) for r in unavailable]),
            "available_noise_alarms":sum(r["post_veto_nominated"] for r in eligible),
            "later_available":sum(u[i]["available"] for i in later),"later_detected":sum(u[i]["post_veto_nominated"] for i in later),
            "raw_candidates_during_unavailability":sum(r["raw_candidate_count"] for r in u if not r["available"]),
            "mature_clean_noise_windows":len(mature),"mature_clean_noise_stft":sum(r["post_veto_nominated"] for r,v in mature),
            "mature_clean_noise_broadband":sum(v["broadband_nominated"] for r,v in mature),"mature_clean_noise_union":sum(v["joint_nominated_post_veto"] for r,v in mature)}
    baseline=load(out/"baseline/used_calibration_timeline.json");oracle=load(out/"without_earlier_four_injections/used_calibration_timeline.json")
    assert all(baseline[i]["score"]==oracle[i]["score"] for i in later)
    assert report["arms"]["baseline"]["later_detected"]==0
    assert report["arms"]["without_earlier_four_injections"]["later_detected"]==5
    assert baseline[367]["threshold"]==163.83940094146237
    assert oracle[367]["threshold"]==1.6580395349128854
    for a in report["arms"].values():
        assert a["later_available"]==6 and a["mature_clean_noise_windows"]==127
        assert a["mature_clean_noise_stft"]==0 and a["mature_clean_noise_broadband"]==a["mature_clean_noise_union"]==127
    for f,h in m["source_sha256"].items(): assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    report["source_hashes_unchanged"]=True;report["status"]="PASS"
    save(out/"independent_verification.json",report)
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
