"""Independent accounting and deterministic recipe regeneration; no detector rerun."""
from __future__ import annotations
import ast,csv,hashlib,importlib.util,json,sys
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime,timezone
import numpy as np

root=Path(sys.argv[1]).resolve()
repo=root.parents[2]
sys.path.insert(0,str(repo/"src"))
spec=importlib.util.spec_from_file_location("replay_driver",repo/"experiments/run_historical_adaptive_replay.py")
driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
def load(p):
    return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):
    assert not p.exists(),f"Do not overwrite: {p}"
    p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding="utf-8")
manifest=load(root/"manifest.json")
assert manifest["status"]=="COMPLETE" and len(manifest["cases"])==12
for name,h in manifest["source_sha256"].items():assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==h
for name,h in manifest["historical_artifact_sha256_before"].items():assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==h
all_windows={}
summaries=[]
for case in manifest["cases"]:
    name=case["case"];p=root/name
    labels=load(p/"input_and_labels.json");windows=load(p/"windows.json")
    with (p/"timeline.csv").open(encoding="utf-8",newline="") as f:timeline=list(csv.DictReader(f))
    assert len(windows)==len(timeline)==case["n_windows"]
    signals=set(labels["signal_indices"])
    isclean=lambda w:w["credibility_status"]==driver.CRED_CALIBRATED and w["glitch_storm_status"]!="GLITCH_STORM" and w["detection_available"]
    # Use the actual is_clean_calibrated contract (GLITCH_QUIET constant).
    from bsd.adaptive_regime_v2b import GLITCH_QUIET
    isclean=lambda w:w["credibility_status"]==driver.CRED_CALIBRATED and w["glitch_storm_status"]==GLITCH_QUIET and w["detection_available"]
    for w,t in zip(windows,timeline):
        d=w["diagnostics"]
        assert float(t["stft_max_score_full_precision"])==w["stft_max_score"]
        assert (t["detection_available"]=="True")==w["detection_available"]
        if d["dedicated_stft_available"]:
            assert w["stft_threshold"] is not None
            assert w["stft_nominated_raw"]==(w["stft_max_score"]>w["stft_threshold"])
        if d["dedicated_broadband_available"]:
            assert w["broadband_threshold"] is not None
            assert w["broadband_nominated"]==(w["broadband_excess_score"]>w["broadband_threshold"])
        assert (t["active_mode_anomaly_flag"]=="True")==bool(w["detection_available"] and
            (w["stft_nominated_post_veto"] if labels["manager_config"]["policy"]=="transient" else w["broadband_nominated"]))
    channel_defs={
        "dedicated_stft":(lambda w:w["diagnostics"]["dedicated_stft_available"],lambda w:w["stft_nominated_post_veto"]),
        "dedicated_broadband":(lambda w:w["diagnostics"]["dedicated_broadband_available"],lambda w:w["broadband_nominated"]),
        "full_union":(lambda w:w["diagnostics"]["union_mode_status"]=="FULL_UNION_AVAILABLE",lambda w:w["joint_nominated_post_veto"])}
    channels={}
    for ch,(avail,nom) in channel_defs.items():
        noise=[w for w in windows if w["index"] not in signals]
        sig=[w for w in windows if w["index"] in signals]
        an=[w for w in noise if avail(w)];asi=[w for w in sig if avail(w)]
        clean=[w for w in an if isclean(w)]
        mature=[w for w in clean if w["regime_id"]>=1 and w["diagnostics"]["is_mature"]]
        count=lambda seq:sum(bool(nom(w)) for w in seq)
        channels[ch]={"noise_total":len(noise),"noise_available":len(an),"noise_abstained":len(noise)-len(an),
            "noise_fp":count(an),"noise_fraction":count(an)/len(an) if an else None,
            "noise_available_exposure_s":4*len(an),"clean_noise_available":len(clean),"clean_noise_fp":count(clean),
            "mature_clean_noise_available":len(mature),"mature_clean_noise_fp":count(mature),
            "mature_clean_noise_fraction":count(mature)/len(mature) if mature else None,
            "signal_total":len(sig),"signal_available":len(asi),"signal_abstained":len(sig)-len(asi),
            "signal_detected":count(asi),"signal_missed_available":len(asi)-count(asi),
            "signal_efficiency_available":count(asi)/len(asi) if asi else None}
    off=labels["step_index"] or 0
    after=[w for w in windows if w["index"]>=off]
    first=lambda f:next((w["index"]-off for w in after if f(w)),None)
    timing={"first_new_regime_calibrated":first(lambda w:w["regime_id"]>=1 and w["credibility_status"]==driver.CRED_CALIBRATED),
        "first_new_regime_stft_available":first(lambda w:w["regime_id"]>=1 and w["diagnostics"]["dedicated_stft_available"]),
        "first_new_regime_full_union_available":first(lambda w:w["regime_id"]>=1 and w["diagnostics"]["union_mode_status"]=="FULL_UNION_AVAILABLE"),
        "first_mature_available":first(lambda w:w["regime_id"]>=1 and w["diagnostics"]["is_mature"] and w["diagnostics"]["dedicated_stft_available"])}
    post=[]
    for i in labels["post_step_injection_indices"]:
        w=windows[i];av=w["diagnostics"]["dedicated_stft_available"]
        post.append({"window_index":i,"available":av,"clean":isclean(w),"mature":w["diagnostics"]["is_mature"],
            "detected":bool(av and w["stft_nominated_post_veto"]),"score":w["stft_max_score"],"threshold":w["stft_threshold"],
            "diagnostic_candidate_count":len(w["raw_stft_candidates"])})
    summary={"case":name,"n_windows":len(windows),"policy":labels["manager_config"]["policy"],"channels":channels,
        "timing_relative_to_step":timing,"post_step_injections":post,
        "union_status_counts":{s:sum(w["diagnostics"]["union_mode_status"]==s for w in windows) for s in sorted({w["diagnostics"]["union_mode_status"] for w in windows})},
        "clean_windows":sum(isclean(w) for w in windows),"detection_available_windows":sum(w["detection_available"] for w in windows),
        "glitch_warning_windows":sum(w["glitch_storm_status"]!=GLITCH_QUIET for w in windows),
        "regimes":max(w["regime_id"] for w in windows)+1,"glitch_ref_refreshes":windows[-1]["diagnostics"]["glitch_ref_refresh_count"],
        "unavailable_diagnostic_candidate_windows":sum(not w["detection_available"] and bool(w["raw_stft_candidates"]) for w in windows)}
    summaries.append(summary);all_windows[name]=windows

# Regenerate only the original input recipes; feed saved outputs back to aggregators.
# No detector calibration/analyze/process_sequence is executed in this verification.
regenerated=[]
clean_calls=[]
call_index=0
initial=all_windows[driver.A_NAMES[0]][0]["calibration_uncertainty"]
def clean_record(*a,**kw):
    value=driver.generate_clean_signal(*a,**kw)
    clean_calls.append({"t_start_s":float(value[1]),"t_end_s":float(value[2]),
        "f_min_hz":float(value[3]),"f_max_hz":float(value[4]),"parameters":value[5]})
    return value
class SavedOutputManager:
    def __init__(self,config=None):self.config=config
    def fit_continuum_reference(self,reference):return None
    def process_sequence(self,observations,*args,**kwargs):
        global call_index
        name=driver.PLAN[call_index];call_index+=1
        expected=load(root/name/"input_and_labels.json")
        assert driver.array_hash(observations)==expected["input_array"],name
        regenerated.append({"case":name,"input_hash_matched":True,"generated_chirps":list(clean_calls)})
        clean_calls.clear()
        result=[]
        for w in all_windows[name]:
            obj=SimpleNamespace(**w)
            obj.calibration_uncertainty=SimpleNamespace(**w["calibration_uncertainty"])
            obj.is_clean_calibrated=w["credibility_status"]==driver.CRED_CALIBRATED and w["glitch_storm_status"]==GLITCH_QUIET and w["detection_available"]
            result.append(obj)
        return result
def saved_calibration(det,bb,reference,*a,**kw):
    assert driver.array_hash(reference)==load(root/"reference_calibration.json")["input_array"]
    return SimpleNamespace(transient_stft_uncertainty=SimpleNamespace(**initial))
namespace={k:v for k,v in vars(driver).items() if not k.startswith("__")}
namespace.update(AdaptiveRegimeManagerV2b=SavedOutputManager,calibrate_separate_search_modes=saved_calibration,
    fit_glitch_reference=lambda _:None,generate_clean_signal=clean_record)
tree=ast.parse((root/"replayed_synthetic_functions.py").read_text(encoding="utf-8"))
exec(compile(tree,"replayed_synthetic_functions.py","exec"),namespace)
sa=namespace["_run_suite_a_corrected"]()
sb=namespace["_run_suite_b_corrected"](sa["shared_ref"])
assert call_index==12
for label,rows in (("a",sa["rows"]),("b",sb["rows"])):
    assert json.loads(json.dumps(rows,default=driver.norm))==load(root/f"suite_{label}_legacy_aggregation_current_outputs.json")
# Compare only the permitted A/B historical evidence; old O3 audit/real-data fields are not used.
historical=load(repo/"results/v03_phase2b_corrected/results.json")
comparisons=[]
for suite,current in (("suite_a",sa["rows"]),("suite_b",sb["rows"])):
    for old,new in zip(historical[suite],current):
        diffs={}
        for key,value in old.items():
            if key in new:
                a,b=value,new[key]
                equal=(abs(float(a)-float(b))<1e-12 if isinstance(a,(int,float)) and isinstance(b,(int,float)) else a==b)
                if not equal:diffs[key]={"historical":a,"current_legacy_aggregation":b}
        comparisons.append({"suite":suite,"scenario":new.get("scenario",str(new.get("background"))+"__"+str(new.get("morphology"))),"differences":diffs})
save(root/"independent_analysis.json",{"status":"VERIFIED","source_and_historical_hashes_match":True,
    "cases":summaries,"input_regeneration":regenerated,"historical_comparisons":comparisons,
    "limitations":["Historical original input bytes/environment are unavailable; current regenerated inputs match exactly",
        "Legacy aggregates use historical denominator semantics and are not the primary corrected statistics",
        "Single seed per suite; no real data or scientific novelty established"]})

lines=["# Historical adaptive Phase 2b — synthetic replay","",
    f"Run: `{manifest['run_id']}`. Current source `043f0eb`; historical recipe `edd1ec9`.",
    "Twelve original synthetic cases completed without detector tuning. Source and historical artifact hashes match. All twelve current input arrays and the reference regenerated exactly. Original historical array hashes are unavailable; bitwise identity with old execution is not established.",
    "", "## Channel-specific current measurements", "",
    "| Case | Clean / total | STFT noise FP / available | STFT abstained noise | STFT detected / available injections | Full Union noise FP / available |",
    "| --- | --- | --- | --- | --- | --- |"]
for s in summaries:
    a=s["channels"]["dedicated_stft"];u=s["channels"]["full_union"]
    lines.append(f"| {s['case']} | {s['clean_windows']}/{s['n_windows']} | {a['noise_fp']}/{a['noise_available']} | {a['noise_abstained']} | {a['signal_detected']}/{a['signal_available']} (abstained {a['signal_abstained']}) | {u['noise_fp']}/{u['noise_available']} |")
ab=next(s for s in summaries if s["case"]=="abrupt_step_500win")
lines+=["","## Abrupt step: timings and fixed post-step injections","",
    "Timings below are zero-based window-start offsets from the shift at index100; one window is4 seconds. Credibility labels and actual available thresholds are different definitions.",
    "",json.dumps(ab["timing_relative_to_step"],indent=2),"",
    "| Injection index | Available | Clean | Mature | Detected | Score | Used threshold |",
    "| --- | --- | --- | --- | --- | --- | --- |"]
for r in ab["post_step_injections"]:
    lines.append(f"| {r['window_index']} | {r['available']} | {r['clean']} | {r['mature']} | {r['detected']} | {r['score']:.8g} | {r['threshold']} |")
lines+=["","## Historical comparison","",
    "Only nonidentical historical fields follow. These compare the old aggregation recipe on current outputs; corrected available/clean/mature denominators are in independent_analysis.json.",
    "",json.dumps([c for c in comparisons if c["differences"]],indent=2),"",
    "## Interpretation limits","",
    "- This replay is a previously inspected development recipe with one seed per suite, not independent confirmation or a blind search.",
    "- No new physical signal or astrophysical event can be established from synthetic inputs.",
    "- Changed recovery timing, availability and thresholds are software/method observations. They do not establish methodological novelty.",
    "- Zero observed mature false alarms is conditional on the recorded available denominator, not proof of zero population FAR.",
    "- Glitch-storm and drifting-background counts are robustness diagnostics; availability alone does not validate exchangeability or stationary FAP.",
    "- Persistent broadband versus background drift remains ambiguous. Association of missed injections with resonant noise is a hypothesis; per-injection frequency metadata are preserved in independent_analysis.json, not a causal intervention.",
    "- Old A/B results remain preserved. Real-data suiteD, coherence suiteC and the earlier phase2/initial phase2b implementations were not rerun.",
    "- No O3 input/manifest, commit, push, release or Zenodo operation was involved.",""]
assert not (root/"report.md").exists()
(root/"report.md").write_text("\n".join(lines),encoding="utf-8")
print("ANALYSIS_VERIFIED")
for s in summaries:
    if s["case"] in driver.A_NAMES:
        print(s["case"],json.dumps({"channels":s["channels"],"timing":s["timing_relative_to_step"]}))
print("HISTORICAL_FIELDS_CHANGED",json.dumps([c for c in comparisons if c["differences"]]))
