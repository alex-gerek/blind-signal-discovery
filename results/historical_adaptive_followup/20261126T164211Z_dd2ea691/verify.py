"""Verify paired outputs and write bounded mechanistic findings."""
import hashlib,json,sys
from pathlib import Path
from datetime import datetime,timezone
p=Path(sys.argv[1]).resolve()
repo=p.parents[2]
def load(x):return json.loads(x.read_text(encoding="utf-8"),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
m=load(p/"manifest.json");assert m["status"]=="COMPLETE"
for f,h in m["source_sha256"].items():assert hashlib.sha256((repo/f).read_bytes()).hexdigest()==h
arms={}
for arm in ("baseline","without_earlier_four_injections"):
    w=load(p/arm/"windows.json");out=load(p/arm/"outcome.json");tr=load(p/arm/"calibration_trace.json")
    for r in w:
        if r["diagnostics"]["dedicated_stft_available"]:
            assert r["stft_nominated_raw"]==(r["stft_max_score"]>r["stft_threshold"])
    later=[r for r in w if r["index"] in m["later_paired_injection_indices"]]
    assert out["later_detected"]==sum(r["detection_available"] and r["stft_nominated_post_veto"] for r in later)
    cal=[r for r in tr if r["target_fap"]==.01][-1]
    first=next(r["index"] for r in w if r["regime_id"]>0 and r["diagnostics"]["is_mature"])
    sig=set([108,125,145]+list(range(330,430,10)))
    noise=[r for r in w if r["index"]>=first and r["index"] not in sig and r["detection_available"]]
    arms[arm]={"mature_threshold":cal["global_threshold"],"threshold_split_size":cal["n_threshold"],
        "known_injected_threshold_members":[v["window_index"] for v in cal["known_injected_members"] if v["partition"]=="THRESHOLD_B"],
        "mature_start_index":first,"later_detected":out["later_detected"],"later_available":sum(r["detection_available"] for r in later),
        "later_scores":[r["stft_max_score"] for r in later],"mature_noise_false_alarms":sum(r["stft_nominated_post_veto"] for r in noise),
        "mature_noise_available":len(noise),"last_ljung_box_p":w[-1]["calibration_uncertainty"]["ljung_box_p"]}
assert arms["baseline"]["known_injected_threshold_members"]==[330,340,350,360]
assert arms["without_earlier_four_injections"]["known_injected_threshold_members"]==[]
assert arms["baseline"]["later_scores"]==arms["without_earlier_four_injections"]["later_scores"]
assert arms["baseline"]["mature_start_index"]==arms["without_earlier_four_injections"]["mature_start_index"]==367
parent=repo/"results/historical_adaptive_replay"/m["parent_run"]
assert load(p/"baseline/windows.json")==load(parent/"abrupt_step_500win/windows.json")
result={"status":"VERIFIED","verified_at_utc":datetime.now(timezone.utc).isoformat(),"arms":arms,
    "confirmed":"In this paired case the four earlier injected signals entered threshold SplitB and changed the threshold, not the six later scores; later detections changed0/6 to5/6.",
    "limitations":["One post-hoc exploratory intervention on a previously used synthetic seed",
        "No universal prevalence or ready operational admission rule established",
        "Oracle knowledge of injected windows is unavailable on unknown real input",
        "No evidence of a new physical signal or methodological novelty"]}
assert not (p/"independent_verification.json").exists()
(p/"independent_verification.json").write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
a,b=arms.values()
report=f"""# Adaptive calibration contamination — exploratory paired follow-up

Source043f0eb; parent replay{m['parent_run']}. Protocol was saved after observing the parent failure and before this intervention. No detector parameter or source change.

## Confirmed observations

| Measurement | Original injected sequence | Four earlier injections replaced with exact original noise |
| --- | --- | --- |
| Injected members in threshold SplitB |330,340,350,360 |none |
| SplitB size |199 |199 |
| Mature threshold |{a['mature_threshold']:.12g} |{b['mature_threshold']:.12g} |
| Mature first used window |367 |367 |
| Same six later injections detected |0/6 |5/6 |
| Later mature noise false alarms |0/127 |0/127 |
| Calibration Ljung-Box p |{a['last_ljung_box_p']:.8g} |{b['last_ljung_box_p']:.8g} |

All samples outside the four intervention windows were asserted bitwise identical. The baseline's complete500-window outputs matched the parent exactly. The six later scores were EXACTLY identical across arms. Therefore the detection change is explained by changed thresholds, rather than changes in those later waveforms/scores or recovery timing.

The original sequence fed its own known signals into the estimated noise tail. In this constructed example, a conservative contaminated threshold hid five otherwise detectable later injections. Mature/available does not imply recovered sensitivity. The old resonance-only explanation does not account for this controlled difference, although spectral masking may still matter for other misses.

## Scope and practical impact

This is a confirmed mechanism in one post-hoc synthetic paired experiment, not proof of its frequency in real streams or a novel physical discovery. The control used oracle injection labels and is NOT a deployable fix. Simply deleting high-score windows can bias the noise calibration and break false-alarm control; do not adopt that rule without a separate design/protocol.

Keep the baseline and failed0/6 result. Before advertising reliable adaptive recovery, design and independently validate admission/segmentation for recalibration, signal retention and uncertainty. Include stationary/no-signal extremes, recurrent true signals, drift and shifted_ar2, with held-out sensitivity and noise exposure. Preserve all unavailable decisions and candidate diagnostics.

The frozen broadband/Union false alarms in the parent are a separate limitation:127/127 mature noise windows were nominated. They do not establish astrophysical signals or recovered1% Union FWER after a changed background.

No original data, O3, commit, push, release or Zenodo operation.
"""
assert not (p/"report.md").exists();(p/"report.md").write_text(report,encoding="utf-8")
print(json.dumps(result,indent=2))
