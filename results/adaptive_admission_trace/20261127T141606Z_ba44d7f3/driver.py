"""Measurement-only calibration membership and used-state tracing; no policy changes."""
from __future__ import annotations
import ast,csv,hashlib,importlib.util,json
from dataclasses import asdict
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import sys,time,uuid,platform,traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
spec=importlib.util.spec_from_file_location("adaptive_replay",ROOT/"experiments/run_historical_adaptive_replay.py")
replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
from bsd.adaptive_regime_v2b import GLITCH_QUIET
PARENT=ROOT/"results/historical_adaptive_replay/20261126T163629Z_75c7a756"
PROTOCOL=ROOT/"docs/adaptive_admission_trace_protocol_2026-11-27.md"
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def save(p,x):replay.save(p,x)
def main():
    assert replay.git("rev-parse","HEAD")=="d049e2b1c97626b562d234a22bc12db5d51702dd"
    assert not replay.git("diff",replay.BASELINE,"--","src","app","pyproject.toml")
    ident=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    out=ROOT/"results/adaptive_admission_trace"/ident;out.mkdir(parents=True,exist_ok=False)
    global RUN_OUTPUT
    RUN_OUTPUT=out
    (out/"protocol.md").write_bytes(PROTOCOL.read_bytes());(out/"driver.py").write_bytes(Path(__file__).read_bytes())
    record={"run_id":ident,"source_commit":replay.git("rev-parse","HEAD"),"dirty_paths":replay.git("status","--short"),"python_version":platform.python_version(),"numpy_version":replay.np.__version__,"scipy_version":replay.scipy.__version__,"parent_run":PARENT.name,"status":"RUNNING",
        "source_sha256":load(PARENT/"manifest.json")["source_sha256"],"parent_manifest_sha256":replay.sha((PARENT/"manifest.json").read_bytes()),
        "driver_sha256":replay.sha(Path(__file__).read_bytes()),"protocol_sha256":replay.sha(PROTOCOL.read_bytes()),
        "post_hoc":"Specified after observing mature threshold inflation; single exploratory paired intervention, no detector tuning",
        "removed_signal_indices":[330,340,350,360],"later_paired_injection_indices":[370,380,390,400,410,420],"arms":[]}
    save(out/"manifest.json",record)
    print("RUN_DIRECTORY",out,flush=True)
    state={"targets":{},"raw_shifted":None,"reference":None,"call_index":0}
    initial=load(PARENT/"stationary_null_200win/windows.json")[0]["calibration_uncertainty"]
    def gen_record(*args,**kwargs):
        obs=replay.generate_noise(*args,**kwargs)
        return obs
    class RecipeOnlyManager:
        def __init__(self,config=None):pass
        def fit_continuum_reference(self,ref):return None
        def process_sequence(self,observations,*args,**kwargs):
            name=replay.A_NAMES[state["call_index"]];state["call_index"]+=1
            state["targets"][name]=observations
            expected=load(PARENT/name/"input_and_labels.json")["input_array"]
            assert replay.array_hash(observations)==expected
            rows=[]
            for w in load(PARENT/name/"windows.json"):
                o=SimpleNamespace(**w);o.calibration_uncertainty=SimpleNamespace(**w["calibration_uncertainty"])
                o.is_clean_calibrated=w["credibility_status"]==replay.CRED_CALIBRATED and w["glitch_storm_status"]==GLITCH_QUIET and w["detection_available"]
                rows.append(o)
            return rows
    def fake_cal(det,bb,ref,*args,**kwargs):
        state["reference"]=ref
        assert replay.array_hash(ref)==load(PARENT/"reference_calibration.json")["input_array"]
        return SimpleNamespace(transient_stft_uncertainty=SimpleNamespace(**initial))
    ns={k:v for k,v in vars(replay).items() if not k.startswith("__")}
    ns.update(AdaptiveRegimeManagerV2b=RecipeOnlyManager,calibrate_separate_search_modes=fake_cal,fit_glitch_reference=lambda _:None)
    exec(compile(ast.parse((PARENT/"replayed_synthetic_functions.py").read_text(encoding="utf-8")),"synthetic_recipe","exec"),ns)
    base_gen=ns["_gen_obs_list"]
    def record_gen(*args,**kwargs):
        obs=base_gen(*args,**kwargs)
        if args[-1]=="ab1":state["raw_shifted"]=[o.signal.copy() for o in obs]
        return obs
    ns["_gen_obs_list"]=record_gen
    ns["_run_suite_a_corrected"]()
    reference=state["reference"];baseline=state["targets"]["abrupt_step_500win"]
    control=list(baseline)
    for index in record["removed_signal_indices"]:
        old=baseline[index]
        control[index]=replay.Observation(old.sample_id,state["raw_shifted"][index-100].copy(),old.sample_rate,old.duration)
    # All later evaluated injections and all other inputs are bitwise identical.
    for i in range(500):
        if i not in record["removed_signal_indices"]:assert replay.np.array_equal(baseline[i].signal,control[i].signal)
    truth=set([108,125,145]+list(range(330,430,10)))
    true_ids={baseline[i].sample_id:i for i in truth}
    initial_source=record["source_sha256"]
    assert all(replay.sha((ROOT/f).read_bytes())==h for f,h in initial_source.items())
    (out/"replay_helper.py").write_bytes((ROOT/"experiments/run_historical_adaptive_replay.py").read_bytes())
    for arm,obs in (("baseline",baseline),("without_earlier_four_injections",control)):
        dest=out/arm;dest.mkdir(exist_ok=False);traces=[];used=[];cal_ids={};cal_objects=[]
        class TracedDetector(replay.MultiResolutionDetector):
            def calibrate(self,observations,*args,**kwargs):
                cal=super().calibrate(observations,*args,**kwargs)
                caller=sys._getframe(1)
                fit_index=caller.f_locals.get("i") if caller.f_code.co_name=="process_sequence" else None
                call_id=f"cal_{len(traces):03d}"
                cal_ids[id(cal)]=call_id;cal_objects.append(cal)
                save(dest/(call_id+".json"),asdict(cal))
                ids=[o.sample_id for o in observations]
                entries=[{"window_index":true_ids[s],"sample_id":s,
                    "partition":"NORMALIZATION_A" if k<cal.n_norm_samples else "THRESHOLD_B"}
                    for k,s in enumerate(ids) if s in true_ids and not (arm!="baseline" and true_ids[s] in record["removed_signal_indices"])]
                traces.append({"calibration_id":call_id,"fit_window_index":fit_index,"earliest_use_index":0 if fit_index is None else fit_index+1,"input_hashes":[replay.sha(replay.np.asarray(o.signal,dtype="<f8").tobytes()) for o in observations],"target_fap":cal.target_fap,"n_total":len(ids),"n_norm":cal.n_norm_samples,
                    "n_threshold":cal.n_threshold_samples,"global_threshold":cal.global_threshold,
                    "sample_ids":ids,"known_injected_members":entries})
                return cal
            def analyze(self,observation,calibration):
                result=super().analyze(observation,calibration)
                caller=sys._getframe(1)
                if caller.f_code.co_name=="process_sequence" and caller.f_locals.get("used_stft_cal") is calibration:
                    available=bool(calibration.target_fap_attainable and calibration.global_threshold is not None)
                    used.append({"window_index":caller.f_locals["i"],"sample_id":observation.sample_id,
                        "calibration_id":cal_ids[id(calibration)],"score":result.max_global_score,
                        "threshold":calibration.global_threshold,"empirical_p":calibration.empirical_p_value(result.max_global_score)[0] if available else None,
                        "p_unavailable_reason":None if available else "CALIBRATION_UNAVAILABLE"})
                return result
        det=TracedDetector();bb=replay.FrozenReferenceBroadbandChannel()
        scal=replay.calibrate_separate_search_modes(det,bb,reference,target_fap_transient=.01,
            target_fap_persistent=.01,target_fwer_union=.01,norm_fraction=.30,calibration_seed=3026)
        mgr=replay.AdaptiveRegimeManagerV2b(replay.AdaptiveRegimeV2bConfig(mode="causal_trailing",policy="transient"))
        start=time.monotonic()
        windows=mgr.process_sequence(obs,det,bb,scal,mgr.fit_continuum_reference(reference),replay.fit_glitch_reference(reference))
        serialized=json.loads(json.dumps([asdict(w) for w in windows],default=replay.norm,allow_nan=False))
        if arm=="baseline":assert serialized==load(PARENT/"abrupt_step_500win/windows.json")
        old_arm=ROOT/"results/historical_adaptive_followup/20261126T164211Z_dd2ea691"/arm
        assert serialized==load(old_arm/"windows.json")
        assert len(used)==len(obs)==500
        for row,w in zip(used,serialized):
            assert row["window_index"]==w["index"] and row["sample_id"]==w["sample_id"]
            assert row["score"]==w["stft_max_score"] and row["threshold"]==w["stft_threshold"]
            row.update(start_s=4.*w["index"],end_s=4.*(w["index"]+1),
                available=w["detection_available"],mature=w["diagnostics"]["is_mature"],
                known_signal=w["index"] in truth and not(arm!="baseline" and w["index"] in record["removed_signal_indices"]),
                raw_nominated=w["stft_nominated_raw"],post_veto_nominated=w["stft_nominated_post_veto"],
                raw_candidate_count=len(w["raw_stft_candidates"]),
                union_mode_status=w["diagnostics"]["union_mode_status"],
                union_stft_threshold=w["diagnostics"]["union_stft_threshold"],
                updated_for_next=w["diagnostics"]["calibration_updated_for_next_window"])
        save(dest/"used_calibration_timeline.json",used)
        replay.csv_save(dest/"timeline.csv",replay._build_rows(windows,4.,"transient")[1])
        save(dest/"windows.json",serialized);save(dest/"calibration_trace.json",traces)
        results=[]
        for i in record["later_paired_injection_indices"]:
            w=windows[i]
            results.append({"window_index":i,"available":w.detection_available,"detected":bool(w.detection_available and w.stft_nominated_post_veto),
                "score":w.stft_max_score,"threshold":w.stft_threshold,"mature":w.diagnostics["is_mature"]})
        entry={"arm":arm,"input_array":replay.array_hash(obs),"elapsed_s":time.monotonic()-start,
            "later_injections":results,"later_detected":sum(r["detected"] for r in results),
            "baseline_matches_parent_exactly":True if arm=="baseline" else None,"matches_previous_followup_exactly":True}
        save(dest/"outcome.json",entry);record["arms"].append(entry);save(out/"manifest.json",record)
        print("ARM_COMPLETE",arm,"later_detected",entry["later_detected"],flush=True)
    assert all(replay.sha((ROOT/f).read_bytes())==h for f,h in initial_source.items())
    record["status"]="COMPLETE";record["source_hashes_unchanged"]=True
    save(out/"manifest.json",record)
    print("RUN_STATUS COMPLETE",flush=True)
if __name__=="__main__":
    try:
        main()
    except Exception:
        if "RUN_OUTPUT" in globals():
            save(RUN_OUTPUT/"failure.json",{"status":"FAILED","traceback":traceback.format_exc()})
        raise
