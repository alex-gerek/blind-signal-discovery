"""Saved-data source-control accounting, independent from core calibration/decision methods."""
from pathlib import Path
import hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
ARMS=["stream_contaminated","paired_noise_oracle","separate_matched","separate_old"]
def load(p):return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding="utf-8")
def associated(c,s):
    return s is not None and max(c["t_start"],s[0])<=min(c["t_end"],s[1]) and max(c["f_min"],s[2])<=min(c["f_max"],s[3])
def main():
    out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/"results/separate_background_control"
    m=load(out/"manifest.json");assert m["status"]=="COMPLETE" and len(m["cases"])==12
    report={"status":"RUNNING","calibrations":0,"evaluation_records":0,"candidate_records":0,"input_records":0,"groups":{},"paired_association":{},"scope":"Descriptive synthetic fixed-source control, three calibration seeds per profile"}
    for case in m["cases"]:
        dest=out/f"{case['seed']}_{case['case']}"
        cf=load(dest/"calibration_features.json");tf=load(dest/"target_features.json")
        assert len(cf)==len({r["sample_id"] for r in cf})==1060
        assert len(tf)==len({r["sample_id"] for r in tf})==240
        by_id={r["sample_id"]:r for r in tf};split={arm:[r for r in cf if r["arm"]==arm] for arm in ARMS}
        for arm,a in split.items():
            assert len(a)==265 and [r["index"] for r in a]==list(range(265))
        assert [i for i,(a,b) in enumerate(zip(split[ARMS[0]],split[ARMS[1]])) if a["input_sha256"]!=b["input_sha256"]]==[100,140,180,220]
        assert all(a["raw_maxima"]==b["raw_maxima"] for a,b in zip(split[ARMS[0]][:66],split[ARMS[1]][:66]))
        bh={arm:{r["underlying_noise_sha256"] for r in a} for arm,a in split.items()}
        target_hashes={r["base_noise_sha256"] for r in tf}
        assert len(target_hashes)==220 and all(not(h&target_hashes) for h in bh.values())
        for a,b in [(ARMS[1],ARMS[2]),(ARMS[1],ARMS[3]),(ARMS[2],ARMS[3])]:assert not(bh[a]&bh[b])
        freeze=load(dest/"threshold_freeze.json");assert freeze["target_decisions_measured"]==0
        assert len(freeze["calibration_hashes"])==8
        assert all(hashlib.sha256((dest/f).read_bytes()).hexdigest()==h for f,h in freeze["calibration_hashes"].items())
        evaluated={}
        for summary in case["summaries"]:
            arm=summary["arm"];alpha=summary["alpha"];name=f"{arm}_{alpha}"
            cal=load(dest/f"calibration_{name}.json");rows=load(dest/f"evaluation_{name}.json")
            a=split[arm];assert len(rows)==240 and cal["n_norm_samples"]==66 and cal["n_threshold_samples"]==199
            stats={}
            for scale,s in cal["scale_stats"].items():
                raw=[r["raw_maxima"][scale] for r in a[:66]];med=float(np.median(raw));q=float(np.quantile(raw,.95));spread=max(q-med,1.)
                assert (s["median_null_score"],s["q95_null_score"],s["spread_scale"])==(med,q,spread)
                stats[scale]=(med,spread)
            def score(r):return max((r["raw_maxima"][k]-med)/spread for k,(med,spread) in stats.items())
            null=sorted(score(r) for r in a[66:]);assert null==cal["sorted_global_null_scores"]
            allowed=max(i for i in range(200) if i/200.<=alpha);assert allowed>=1
            threshold=null[200-allowed-1];assert threshold==cal["global_threshold"]==summary["threshold"]
            for row in rows:
                raw=by_id[row["sample_id"]];assert all(row[k]==v for k,v in raw.items() if k!="raw_maxima")
                value=score(raw);p=(1+sum(x>=value for x in null))/200.
                assert value==row["score"] and p==row["p"] and row["threshold"]==threshold
                assert row["available"] and row["detected"]==bool(value>threshold and p<=alpha)==bool(row["candidates"])
                hit=False
                for c in row["candidates"]:
                    med,spread=stats[c["scale_name"]]
                    assert c["standardized_score"]==(c["raw_anomaly_score"]-med)/spread
                    assert 0<=c["t_start"]<=c["t_end"]<=4 and 0<=c["f_min"]<=c["f_max"]<=512
                    assert c["n_bins"]>=2
                    cp=(1+sum(x>=c["standardized_score"] for x in null))/200.
                    assert c["standardized_score"]>threshold and cp==c["empirical_p_value"] and cp<=alpha
                    hit=hit or associated(c,row["support"]);report["candidate_records"]+=1
                assert row["associated"]==bool(hit)
                report["evaluation_records"]+=1
            negative=[r for r in rows if r["part"]=="noise"]
            assert len(negative)==200 and summary["available_noise_windows"]==200 and summary["available_noise_seconds"]==800.
            assert summary["noise_alarms"]==sum(r["detected"] for r in negative)
            k=f"{case['case']}_{arm}_fap{alpha}"
            group=report["groups"].setdefault(k,{"calibration_seeds":0,"noise_windows":0,"noise_alarms":0,"available_noise_seconds":0.,"signals":{}})
            group["calibration_seeds"]+=1;group["noise_windows"]+=200;group["noise_alarms"]+=summary["noise_alarms"];group["available_noise_seconds"]+=800.
            for sg in summary["groups"]:
                sub=[r for r in rows if r["family"]==sg["family"] and r["snr"]==sg["snr"]]
                assert len(sub)==sg["total"]==10 and sg["abstained"]==0
                assert sg["detected"]==sum(r["detected"] for r in sub) and sg["associated"]==sum(r["associated"] for r in sub)
                sk=f"{sg['family']}_snr{sg['snr']}";g=group["signals"].setdefault(sk,{"total":0,"detected":0,"associated":0})
                for field in g:g[field]+=sg[field]
            evaluated[name]=rows;report["calibrations"]+=1
        for alpha in [.01,.005]:
            original=evaluated[f"{ARMS[0]}_{alpha}"]
            for arm in ARMS[1:]:
                candidate=evaluated[f"{arm}_{alpha}"]
                for a,b in zip(original,candidate):assert a["sample_id"]==b["sample_id"] and a["input_sha256"]==b["input_sha256"]
                for family in ["chirp","burst"]:
                    for snr in [8.,18.]:
                        pairs=[(a,b) for a,b in zip(original,candidate) if a["family"]==family and a["snr"]==snr]
                        k=f"{case['case']}_{arm}_fap{alpha}_{family}_snr{snr}"
                        g=report["paired_association"].setdefault(k,{"gained":0,"lost":0,"pairs":0})
                        g["gained"]+=sum(not a["associated"] and b["associated"] for a,b in pairs)
                        g["lost"]+=sum(a["associated"] and not b["associated"] for a,b in pairs)
                        g["pairs"]+=len(pairs)
        for alpha in [.01,.005]:
            c0=load(dest/f"calibration_{ARMS[0]}_{alpha}.json");c1=load(dest/f"calibration_{ARMS[1]}_{alpha}.json")
            assert c0["scale_stats"]==c1["scale_stats"]
            assert all(a["score"]==b["score"] for a,b in zip(evaluated[f"{ARMS[0]}_{alpha}"],evaluated[f"{ARMS[1]}_{alpha}"]))
        report["input_records"]+=len(cf)+len(tf)
    spec=importlib.util.spec_from_file_location("source_control",ROOT/"experiments/run_separate_background_control.py")
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    assert Path(driver.__file__).read_bytes()==(out/"run_separate_background_control.py").read_bytes()
    regenerated=0
    for seed in driver.SEEDS:
        for ci,case in enumerate(driver.CASES):
            calibration,target=driver.generate(seed,ci);dest=out/f"{seed}_{case}"
            cf=load(dest/"calibration_features.json");tf=load(dest/"target_features.json")
            generated=[item for arm in ARMS for item in calibration[arm]]+target
            for (o,meta),stored in zip(generated,cf+tf):
                assert all(stored[k]==v for k,v in meta.items()) and driver.key(o)==stored["input_sha256"]
                regenerated+=1
    assert regenerated==report["input_records"]==15600
    for name,h in {**m["source_sha256"],**m["dependency_sha256"]}.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    report.update(status="PASS",regenerated_inputs=regenerated,source_hashes_unchanged=True,uncached_parity_checks=m["parity_checks"])
    save(out/"independent_verification.json",report)
    print(json.dumps({k:v for k,v in report.items() if k not in ["groups","paired_association"]},indent=2))
    for case in driver.CASES:
        print(case,"FAP0.01",[(arm,report["groups"][f"{case}_{arm}_fap0.01"]["noise_alarms"],{k:v["associated"] for k,v in report["groups"][f"{case}_{arm}_fap0.01"]["signals"].items()}) for arm in ARMS])
if __name__=="__main__":main()
