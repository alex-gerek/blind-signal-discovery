"""Independent exported-row accounting for the bounded control study."""
import csv, hashlib, json, sys
from pathlib import Path
from datetime import datetime, timezone

root = Path(sys.argv[1])
def load(p):
    return json.loads(p.read_text(encoding="utf-8"), parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
manifest = load(root / "manifest.json")
assert manifest["status"] in ("COMPLETE", "COMPLETED_WITH_ERRORS")
outcomes = []
failures = []
paired = {}
adaptive_paired = {}
regenerated = load(root / "regenerated_input_hashes.json")
assert regenerated["status"] == "ALL_SIX_INPUTS_MATCHED"
for name, expected in manifest["source_sha256"].items():
    assert hashlib.sha256((root.parents[2] / name).read_bytes()).hexdigest() == expected
for case in manifest["case_results"]:
    p = root / case["case"]
    x = load(p / "outcome.json")
    if x["status"] != "COMPLETE":
        failures.append({"case":case["case"],"error":x.get("error")})
        continue
    with (p/"timeline.csv").open(encoding="utf-8",newline="") as f:
        rows = list(csv.DictReader(f))
    is_stationary = case["case"].startswith("seed")
    indices = list(range(201,220,2)) if is_stationary else ([] if x.get("capped") else list(range(351,370,2)))
    negatives = [r for r in rows if (int(r["window_index"])<200 if is_stationary else int(r["window_index"]) not in indices)]
    positives = [r for r in rows if int(r["window_index"]) in indices]
    available = lambda r:r["detection_available"]=="True"
    hit = lambda r:available(r) and r["active_mode_anomaly_flag"]=="True"
    n_available = sum(available(r) for r in negatives)
    fp = sum(hit(r) for r in negatives)
    pos_available = sum(available(r) for r in positives)
    tp = sum(hit(r) for r in positives)
    for r in rows:
        score = float(r["stft_max_score_full_precision"])
        threshold = r["stft_threshold_full_precision"]
        if available(r):
            assert threshold != ""
            assert hit(r) == (score > float(threshold)), (case["case"], r["window_index"])
        else:
            assert not hit(r)
    m = x["metrics"]
    assert (fp,n_available,len(negatives)-n_available,tp,pos_available-tp,len(positives)-pos_available) == (
        m["negative_control_false_alarms_available"],m["negative_control_available_windows"],
        m["negative_control_abstained_windows"],
        m["positive_injections_detected"] if is_stationary else sum(r["detected"] for r in x["injections"]),
        m["positive_injections_missed_available"],m["positive_injections_abstained"])
    assert m["negative_control_available_exposure_s"] == 4*n_available
    assert m["negative_control_far_available"] == (fp/n_available if n_available else None)
    assert m["detection_efficiency_available"] == (tp/pos_available if pos_available else None)
    result = {"case":case["case"],"noise_fp":fp,"noise_available":n_available,
        "noise_abstained":len(negatives)-n_available,"noise_available_s":4*n_available,
        "injections_detected":tp,"injections_available":pos_available,
        "injections_missed":pos_available-tp,"injections_abstained":len(positives)-pos_available}
    if is_stationary:
        result.update(seed=x["seed"],snr=x["snr"],calibration_status=x["calibration"]["status"])
        key = x["seed"]
        current = {"background":x["input_hashes"]["background"],"noise_rows":[
            (r["stft_max_score_full_precision"],r["stft_threshold_full_precision"],
             r["detection_available"],r["active_mode_anomaly_flag"]) for r in negatives]}
        if key in paired:
            assert paired[key]==current, "Paired noise/control mismatch"
        paired[key] = current
    else:
        if not x["capped"]:
            if x["seed"] in adaptive_paired:
                assert adaptive_paired[x["seed"]] == x["input_hashes"], "Paired adaptive input mismatch"
            adaptive_paired[x["seed"]] = x["input_hashes"]
        result.update(first_unavailable_index=x["first_unavailable_index"],
            first_new_regime_available_index=x["first_new_regime_available_index"],
            first_mature_available_index=x["first_mature_available_index"],
            unavailable_windows=x["unavailable_windows"], union_status_counts=x["union_status_counts"],
            compatibility_flags=x["compatibility_flags"],
            unavailable_diagnostic_candidate_windows=x["unavailable_diagnostic_candidate_windows"])
        flags=x["compatibility_flags"]
        if x["capped"]:
            if not (flags["resource_block_reported"] and flags["expected_unavailable_at_end"]):
                failures.append({"case":case["case"],"compatibility":"Expected resource-limited abstention missing"})
        elif not (flags["recovered_new_regime"] and flags["actual_strong_injection_detected"]):
            failures.append({"case":case["case"],"compatibility":"Recovery/detection expectation not met"})
    outcomes.append(result)
stationary = [r for r in outcomes if r["case"].startswith("seed")]
unique_noise = [r for r in stationary if r["snr"]==8]
noise_fp = sum(r["noise_fp"] for r in unique_noise)
noise_n = sum(r["noise_available"] for r in unique_noise)
summary = {"verified_at_utc":datetime.now(timezone.utc).isoformat(),
    "accounting":"Independently recomputed from full-precision exported timeline CSV; paired noise counted once per seed",
    "status":"VERIFIED_WITH_FINDINGS" if failures else "ACCOUNTING_VERIFIED",
    "cases_complete":len(outcomes),"cases_planned":15,"findings":failures,
    "input_regeneration": "6/6 exact array hashes matched; paired adaptive input hashes matched across FAP",
    "stationary_noise":{"false_alarms":noise_fp,"available_windows":noise_n,"available_exposure_s":4*noise_n,
        "abstained_windows":sum(r["noise_abstained"] for r in unique_noise),
        "fraction":noise_fp/noise_n if noise_n else None,
        "scope":"Descriptive pooled fraction across three frozen calibrations; no common-rate iid CI claimed"},
    "stationary_injections":{str(snr):{"detected":sum(r["injections_detected"] for r in stationary if r["snr"]==snr),
        "available":sum(r["injections_available"] for r in stationary if r["snr"]==snr),
        "abstained":sum(r["injections_abstained"] for r in stationary if r["snr"]==snr)} for snr in (8,18)},
    "cases":outcomes}
for name in ("independent_verification.json","report.md","artifact_sha256.json"):
    assert not (root/name).exists(), "Do not overwrite an earlier analysis"
(root/"independent_verification.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
lines=["# Post-audit Python Workstation controls — bounded synthetic study",
    "",f"- Run: `{manifest['run_id']}`; source `{manifest['source_commit']}`.",
    "- The detector implementation was unchanged. The new uncommitted driver, protocol and dirty-path list are recorded in manifest.json.",
    f"- Completed cases: {len(outcomes)}/15. Accounting status: {summary['status']}.",
    "- Inputs are generated synthetic arrays. O3-manifest reads and network access were blocked.",
    "", "## Stationary controls",
    "", "| Seed | SNR | Noise FP / available | Noise abstained | Injections detected / available | Injection abstained | Calibration |",
    "| --- | --- | --- | --- | --- | --- | --- |"]
for r in stationary:
    lines.append(f"| {r['seed']} | {r['snr']} | {r['noise_fp']}/{r['noise_available']} | {r['noise_abstained']} | {r['injections_detected']}/{r['injections_available']} | {r['injections_abstained']} | {r['calibration_status']} |")
lines += ["",f"Paired SNR controls share identical background hashes and identical pre-injection exported noise scores/thresholds/decisions. Count them once: {noise_fp}/{noise_n}, available exposure {4*noise_n} s.",
    "The fraction is descriptive across three frozen calibrations; a pooled common-rate Bernoulli confidence interval is not justified by this design.",
    "", "## Adaptive and resource controls", "",
    "| Case | First unavailable | First new-regime available | First mature available | Unavailable windows | Noise FP / available | Injections detected / available | Abstained injections |",
    "| --- | --- | --- | --- | --- | --- | --- | --- |"]
for r in outcomes:
    if r["case"].startswith("seed"):continue
    lines.append(f"| {r['case']} | {r['first_unavailable_index']} | {r['first_new_regime_available_index']} | {r['first_mature_available_index']} | {r['unavailable_windows']} | {r['noise_fp']}/{r['noise_available']} | {r['injections_detected']}/{r['injections_available']} | {r['injections_abstained']} |")
lines += ["", "Indices are zero-based from the amplitude step at target index 0; multiply by 4 seconds for target-relative window start time. Recovery is recorded when the new threshold is actually used, not when it is prepared for the next window.",
    "The adaptive noise denominator includes available decisions before the monitor recognized the shift, as well as available post-recovery noise. It is not a stationary FAP estimate or a mature-only FAR.",
    "Adaptive inputs and injections are paired across FAP settings: six scenario results are not 60 independent injected signals. At union FWER=0.005 the Sidak per-channel alpha is about 0.00250313, requiring at least 399 Split-B scores: 570 initial windows at the 30% normalization split, or 532 recalibration windows at the 25% split. The 320-reference/370-target design cannot establish full Union recovery; unavailable Union is expected, not evidence of a software defect.",
    "Union status counts and resource flags are retained in independent_verification.json; a recovered dedicated STFT channel must not be described as full Union recovery.",
    "", "## Findings and scope", "",
    json.dumps(failures,ensure_ascii=False),
    "", "- This pilot has three independent calibration realizations, one stationary noise model and one amplitude-step model. It does not replace the 100-seed historical confirmatory study.",
    "- Strong synthetic recovery does not establish sensitivity to weak signals, resonant shifted_ar2 noise, persistent broadband signals, real sensors or gravitational-wave events.",
    "- All misses and unavailable cases are retained. No seed was replaced and no detector parameter was tuned.",
    "- The 15 cases are experiments, not 15 additional pytest tests. Existing saved regression evidence remains 135 passed / 1 deselected at the previously recorded source.",
    "- Headless export checks do not establish native Streamlit/browser interaction or decoder behavior. Manual operator verification remains pending.",
    "- Historical adaptive/coincidence/GW170814 and browser stride-decimation outputs remain historical and require separate planned reruns.",
    "- The initial system-Python import failure is recorded in docs/validation/post_audit_control_import_failure_20261126.json.",
    ""]
(root/"report.md").write_text("\n".join(lines),encoding="utf-8")
inventory = {str(p.relative_to(root)).replace("\\","/"):hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(root.rglob("*")) if p.is_file() and p.name!="artifact_sha256.json"}
(root/"artifact_sha256.json").write_text(json.dumps(inventory,indent=2),encoding="utf-8")
print(json.dumps({k:v for k,v in summary.items() if k!="cases"},indent=2))
