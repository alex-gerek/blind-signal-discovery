"""Bounded, append-only synthetic controls; no detector tuning or external data."""
from __future__ import annotations
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import uuid
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.dont_write_bytecode = True

def guard(event, args):
    if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        if "o3_holdout" in str(args[0]).lower():
            raise PermissionError("Reserved-manifest/data access blocked by control protocol")
    if event in ("socket.connect", "socket.getaddrinfo"):
        raise PermissionError("Network access blocked by control protocol")
sys.addaudithook(guard)

import numpy as np
import scipy
from bsd import workstation_runner as wr
from bsd.adaptive_regime_v2b import (
    AdaptiveRegimeManagerV2b, AdaptiveRegimeV2bConfig,
    FrozenReferenceBroadbandChannel, calibrate_separate_search_modes, fit_glitch_reference,
)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, generate_clean_signal, _scale_to_energy_snr

PROTOCOL = ROOT / "docs/post_audit_workstation_control_protocol_2026-11-26.md"
BASELINE = "043f0eb95be8eeae65b9a5f5bb63015d93559b73"

def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

def digest(data):
    return hashlib.sha256(data).hexdigest()

def array_hash(a):
    a = np.asarray(a, dtype="<f8", order="C")
    return {"sha256": digest(a.tobytes(order="C")), "shape": list(a.shape),
            "dtype": "<f8", "byte_order": "little", "encoding": "C-order raw array bytes"}

def json_write(p, data):
    p.write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")

def csv_write(p, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with p.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def source_hashes():
    paths = sorted((ROOT / "src/bsd").glob("*.py")) + [ROOT / "app/bsd_workstation.py", ROOT / "pyproject.toml"]
    return {str(p.relative_to(ROOT)).replace("\\", "/"): digest(p.read_bytes()) for p in paths}

def validate_windows(windows):
    for w in windows:
        available = bool(w.detection_available)
        threshold = w.stft_threshold
        if available and (threshold is None or not np.isfinite(threshold)):
            raise AssertionError("Available STFT has no finite used threshold")
        if threshold is not None and bool(w.stft_nominated_raw) != bool(w.stft_max_score > threshold):
            raise AssertionError("Raw nomination disagrees with used threshold")
        if not available and w.is_clean_calibrated:
            raise AssertionError("Unavailable window included in clean exposure")
    return {"windows_checked": len(windows), "used_threshold_and_availability_checks": "PASS"}

def stationary(case, out, seed, snr):
    captured = {}
    original = wr.run_workstation_analysis
    def forward(target, background, config=None):
        captured.update({"background": array_hash(background.signal), "target": array_hash(target.signal)})
        return original(target, background, config)
    wr.run_workstation_analysis = forward
    try:
        result = wr.run_synthetic_self_test(
            wr.WorkstationConfig(), n_bg_windows=220, n_neg_control_windows=200,
            n_pos_injections=10, noise_model="pink_colored", signal_family="linear_chirp",
            target_snr=snr, seed=seed)
    finally:
        wr.run_workstation_analysis = original
    check = validate_windows(result.windows)
    gt = result.ground_truth_audit
    assert gt["positive_injections_detected"] + gt["positive_injections_missed_available"] + gt["positive_injections_abstained"] == 10
    wr.export_workstation_run(result, out)
    json_write(out / "array_provenance.json", captured)
    # Export contract checked using independently decoded saved bytes.
    saved = json.loads((out / "run_report.json").read_text(encoding="utf-8"))
    assert saved == json.loads(json.dumps(wr.build_workstation_report(result), allow_nan=False))
    rows = list(csv.DictReader(io.StringIO((out / "timeline.csv").read_text(encoding="utf-8"))))
    assert len(rows) == len(result.windows)
    for row, w in zip(rows, result.windows):
        assert float(row["stft_max_score_full_precision"]) == w.stft_max_score
        assert (row["detection_available"] == "True") == bool(w.detection_available)
        assert (row["active_mode_anomaly_flag"] == "True") == bool(w.detection_available and w.stft_nominated_post_veto)
        if w.stft_threshold is not None:
            assert float(row["stft_threshold_full_precision"]) == w.stft_threshold
    return {"case": case, "seed": seed, "snr": snr, "calibration": result.calibration_audit.to_dict(),
            "metrics": gt, "summary": result.summary_metrics, "input_hashes": captured, "checks": check}

def adaptive(case, out, seed, fap, capped=False):
    rng_ref, rng_noise, rng_inj = [np.random.default_rng(seed + k) for k in range(3)]
    def obs(rng, n, sigma, prefix):
        return [Observation(f"{prefix}_{i}", rng.normal(0, sigma, 4096), 1024., 4.) for i in range(n)]
    ref = obs(rng_ref, 320, 1., "ref")
    target = obs(rng_noise, 370, 4., "shift")
    inj_indices = list(range(351, 370, 2))
    for idx in inj_indices:
        clean, *_ = generate_clean_signal("linear_chirp", 4096, 1024., rng_inj)
        target[idx] = Observation(target[idx].sample_id,
            target[idx].signal + _scale_to_energy_snr(clean, 4., 30.), 1024., 4.)
    if capped:
        target = target[:180]
        inj_indices = []
    hashes = {"background": array_hash(np.stack([o.signal for o in ref])),
              "target": array_hash(np.stack([o.signal for o in target]))}
    det, bb = MultiResolutionDetector(), FrozenReferenceBroadbandChannel()
    scal = calibrate_separate_search_modes(det, bb, ref, target_fap_transient=fap,
        target_fap_persistent=fap, target_fwer_union=fap, calibration_seed=0)
    cfg = AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient",
        provisional_recal_windows=60, mature_recal_windows=200,
        max_recal_buffer_windows=150 if capped else 2000, target_fwer=fap)
    mgr = AdaptiveRegimeManagerV2b(cfg)
    windows = mgr.process_sequence(target, det, bb, scal, mgr.fit_continuum_reference(ref), fit_glitch_reference(ref))
    check = validate_windows(windows)
    candidates, timeline = wr._build_rows(windows, 4., "transient")
    csv_write(out / "timeline.csv", timeline)
    json_write(out / "candidate_records.json", candidates)
    json_write(out / "window_diagnostics.json", [w.diagnostics for w in windows])
    json_write(out / "array_provenance.json", hashes)
    records = [{"window_index": idx, "detected": bool(windows[idx].detection_available and windows[idx].stft_nominated_post_veto),
        "detection_available": bool(windows[idx].detection_available),
        "evaluation_status": "ABSTAINED" if not windows[idx].detection_available else ("DETECTED" if windows[idx].stft_nominated_post_veto else "MISSED")}
        for idx in inj_indices]
    noise = [w for i,w in enumerate(windows) if i not in inj_indices]
    metrics = wr._self_test_availability_metrics(noise, records, 4., "transient")
    first = lambda pred: next((i for i,w in enumerate(windows) if pred(w)), None)
    recovery = first(lambda w: w.regime_id >= 1 and w.detection_available)
    mature = first(lambda w: w.regime_id >= 1 and w.detection_available and w.diagnostics.get("is_mature", False))
    flags = {
        "resource_block_reported": any(w.diagnostics.get("resource_limit_blocked_recovery", False) for w in windows),
        "recovered_new_regime": recovery is not None,
        "actual_strong_injection_detected": any(r["detected"] for r in records),
    }
    if capped:
        flags["expected_unavailable_at_end"] = not windows[-1].detection_available and windows[-1].stft_threshold is None
    return {"case": case, "seed": seed, "fap": fap, "capped": capped, "manager_config": cfg.__dict__,
        "input_hashes": hashes, "checks": check, "metrics": metrics, "injections": records,
        "first_unavailable_index": first(lambda w: not w.detection_available),
        "first_new_regime_available_index": recovery, "first_mature_available_index": mature,
        "unavailable_windows": sum(not w.detection_available for w in windows),
        "available_exposure_s": 4 * sum(w.detection_available for w in windows),
        "unavailable_diagnostic_candidate_windows": sum(not w.detection_available and bool(w.raw_stft_candidates) for w in windows),
        "union_status_counts": {s: sum(r["union_mode_status"] == s for r in timeline) for s in sorted({r["union_mode_status"] for r in timeline})},
        "compatibility_flags": flags}

def main():
    if git("rev-parse", "HEAD") != BASELINE:
        raise SystemExit("Unexpected baseline; specify a new protocol before running another revision")
    if git("diff", BASELINE, "--", "src", "app", "pyproject.toml", "tests"):
        raise SystemExit("Implementation/test modifications are outside this protocol")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = f"{stamp}_{uuid.uuid4().hex[:8]}"
    out = ROOT / "results/post_audit_workstation_controls" / run_id
    out.mkdir(parents=True, exist_ok=False)
    (out / "protocol.md").write_bytes(PROTOCOL.read_bytes())
    (out / "driver.py").write_bytes(Path(__file__).read_bytes())
    hashes = source_hashes()
    manifest = {"run_id": run_id, "started_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": git("rev-parse", "HEAD"), "source_tree": git("rev-parse", "HEAD^{tree}"),
        "uncommitted_paths": git("status", "--porcelain").splitlines(),
        "source_sha256": hashes, "driver_sha256": digest(Path(__file__).read_bytes()),
        "protocol_sha256": digest(PROTOCOL.read_bytes()), "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__, "reserved_access": "BLOCKED", "network": "BLOCKED",
        "status": "RUNNING", "case_results": []}
    json_write(out / "manifest.json", manifest)
    print("RUN_DIRECTORY", out, flush=True)
    plan = [(f"seed{s}_snr{snr}", stationary, {"seed":s,"snr":snr}) for s in (2026,2027,2028) for snr in (8,18)]
    plan += [(f"adaptive_seed{s}_fap{fap}", adaptive, {"seed":s,"fap":fap}) for s in (601,611,621) for fap in (.01,.005)]
    plan += [(f"capped_seed{s}_fap0.005", adaptive, {"seed":s,"fap":.005,"capped":True}) for s in (601,611,621)]
    for case, fn, kwargs in plan:
        dest = out / case
        dest.mkdir(exist_ok=False)
        t0 = time.monotonic()
        try:
            result = fn(case, dest, **kwargs)
            result.update(status="COMPLETE", elapsed_s=time.monotonic()-t0)
        except Exception as exc:
            result = {"case":case, "status":"ERROR", "error_type":type(exc).__name__,
                      "error":str(exc), "traceback":traceback.format_exc(), "elapsed_s":time.monotonic()-t0}
        json_write(dest / "outcome.json", result)
        manifest["case_results"].append({"case":case, "status":result["status"], "elapsed_s":result["elapsed_s"]})
        json_write(out / "manifest.json", manifest)
        print(case, result["status"], f'{result["elapsed_s"]:.1f}s', flush=True)
    manifest["source_hashes_unchanged"] = source_hashes() == hashes
    manifest["status"] = "COMPLETE" if all(r["status"]=="COMPLETE" for r in manifest["case_results"]) and manifest["source_hashes_unchanged"] else "COMPLETED_WITH_ERRORS"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    json_write(out / "manifest.json", manifest)
    print("RUN_STATUS", manifest["status"], flush=True)
    return 0 if manifest["status"]=="COMPLETE" else 1

if __name__ == "__main__":
    raise SystemExit(main())
