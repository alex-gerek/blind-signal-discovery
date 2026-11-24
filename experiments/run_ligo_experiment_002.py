"""Experiment 002 — Real LIGO Data Validation (GW150914, GWOSC Official HDF5).

Evaluates the frozen BSD v0.1 detector (`TemplateFreeDetector`) on authentic
LIGO Hanford (H1) and Livingston (L1) strain measurements from GWOSC without
using the published event timestamp or waveform template during calibration
or candidate extraction.

Two explicitly separated preprocessing pipelines are evaluated:
  1. `frozen_v01_minimal`: numerical scaling + decimation to 1024 Hz only.
  2. `adapted_ligo_whitened`: background median-Welch whitening + 30-380 Hz
     bandpass + decimation to 1024 Hz (physical adaptation layer).

Usage:
    python experiments/run_ligo_experiment_002.py
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bsd.calibration import CalibrationResult, calibrate_detector, clopper_pearson_ci, wilson_ci
from bsd.detector import DetectionResult, DetectorConfig, TemplateFreeDetector
from bsd.ligo_adapter import (
    BackgroundModel,
    GWOSCStrainRecord,
    PaddedSegment,
    candidate_to_gps,
    extract_clean_segments,
    fit_background_model,
    load_gwosc_hdf5,
    preprocess_segment,
)

# Published GWOSC reference values — accessed ONLY during Section 7 final evaluation
PUBLISHED_GW150914_GPS = 1126259462.422
PUBLISHED_GW150914_FREQ_BAND_HZ = (35.0, 250.0)
PHYSICAL_LIGHT_TRAVEL_S = 0.010  # ~10 ms max H1-L1 propagation delay

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 180,
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
)

_DATA_NOTE = "Authentic LIGO O1 strain data from GWOSC (GWTC-1 / GW150914 v3) — BSD Experiment 002"


def _add_footer(fig: plt.Figure, note: str = _DATA_NOTE) -> None:
    fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=8, style="italic", color="#444444")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _run_detector_pipeline_for_site(
    detector_name: str,
    rec_32: GWOSCStrainRecord,
    rec_4096: GWOSCStrainRecord,
    bsd_detector: TemplateFreeDetector,
    mode: str,
    target_fap: float = 0.01,
    n_psd_segments: int = 60,
    n_cal_segments: int = 500,
) -> dict[str, Any]:
    """Run strictly separated background fit, calibration, bg validation, and blind 32s scan."""
    # Exclude +-64 s guard around the entire 32-second target file from all background sets
    exclude_interval = (rec_32.gps_start - 64.0, rec_32.gps_end + 64.0)
    bg_all = extract_clean_segments(
        rec_4096,
        segment_duration_s=4.0,
        pad_s=1.0,
        exclude_gps_interval=exclude_interval,
        require_dq_data_bit=True,
    )
    if len(bg_all) < n_psd_segments + n_cal_segments + 50:
        raise RuntimeError(
            f"{detector_name}: insufficient clean background segments ({len(bg_all)} found)."
        )

    psd_segs = bg_all[:n_psd_segments]
    cal_segs = bg_all[n_psd_segments : n_psd_segments + n_cal_segments]
    val_segs = bg_all[n_psd_segments + n_cal_segments :]

    # 1. Fit background model strictly on psd_segs
    bg_model = fit_background_model(
        detector=detector_name,
        bg_segments=psd_segs,
        mode=mode,
        target_sample_rate=1024.0,
        bandpass_hz=(30.0, 380.0) if mode == "adapted_ligo_whitened" else None,
    )

    # 2. Calibrate frozen BSD v0.1 detector on independent cal_segs
    cal_obs = [
        preprocess_segment(seg, f"{detector_name}_{mode}_cal_{idx:04d}", bg_model)
        for idx, seg in enumerate(cal_segs)
    ]
    cal = calibrate_detector(
        detector=bsd_detector,
        calibration_observations=cal_obs,
        target_fap=target_fap,
        calibration_noise_models=(f"ligo_{detector_name}_{mode}",),
        calibration_seed=0,
    )

    # 3. Independent background validation (empirical FAR + background candidate timestamps)
    val_scores: list[float] = []
    val_fp = 0
    val_candidates_gps: list[dict[str, Any]] = []
    val_segment_detected: list[tuple[float, bool, float]] = []

    for idx, seg in enumerate(val_segs):
        seg_id = f"{detector_name}_{mode}_bgval_{idx:04d}"
        obs = preprocess_segment(seg, seg_id, bg_model)
        res = bsd_detector.analyze(obs, calibration=cal, return_maps=False)
        val_scores.append(res.max_anomaly_score)
        val_segment_detected.append((seg.gps_start, res.detected, res.max_anomaly_score))
        if res.detected:
            val_fp += 1
            for cand in res.candidates:
                val_candidates_gps.append(candidate_to_gps(cand, seg.gps_start, detector_name, seg_id))

    n_val = len(val_segs)
    val_far = val_fp / n_val
    val_far_cp = clopper_pearson_ci(val_fp, n_val)
    val_far_w = wilson_ci(val_fp, n_val)

    # 4. Blind scan of the 32-second target file (uniform 4-s grid with 1-s padding inside 32-s file)
    target_segs = extract_clean_segments(
        rec_32,
        segment_duration_s=4.0,
        pad_s=1.0,
        exclude_gps_interval=None,
        require_dq_data_bit=True,
        step_s=4.0,
    )

    target_windows: list[dict[str, Any]] = []
    target_candidates_gps: list[dict[str, Any]] = []
    target_maps_bundle: list[tuple[PaddedSegment, Any, DetectionResult]] = []

    for idx, seg in enumerate(target_segs):
        seg_id = f"{detector_name}_{mode}_target_{idx:02d}"
        obs = preprocess_segment(seg, seg_id, bg_model)
        res = bsd_detector.analyze(obs, calibration=cal, return_maps=True)
        target_maps_bundle.append((seg, obs, res))
        target_windows.append(
            {
                "segment_id": seg_id,
                "detector": detector_name,
                "mode": mode,
                "gps_start": seg.gps_start,
                "gps_end": seg.gps_start + seg.duration_s,
                "detected": res.detected,
                "max_anomaly_score": res.max_anomaly_score,
                "threshold": cal.threshold,
                "n_candidates": len(res.candidates),
                "input_kurtosis": res.diagnostics["input_kurtosis"],
                "input_crest_factor": res.diagnostics["input_crest_factor"],
            }
        )
        for cand in res.candidates:
            target_candidates_gps.append(candidate_to_gps(cand, seg.gps_start, detector_name, seg_id))

    return {
        "detector": detector_name,
        "mode": mode,
        "bg_model": bg_model,
        "calibration": cal,
        "n_psd_segments": len(psd_segs),
        "n_cal_segments": len(cal_segs),
        "n_val_segments": n_val,
        "val_fp": val_fp,
        "val_far": val_far,
        "val_far_cp_95ci": list(val_far_cp),
        "val_far_wilson_95ci": list(val_far_w),
        "val_statistically_credible": bool(val_far_cp[0] <= target_fap),
        "val_scores": val_scores,
        "val_segment_detected": val_segment_detected,
        "val_candidates_gps": val_candidates_gps,
        "target_windows": target_windows,
        "target_candidates_gps": target_candidates_gps,
        "target_maps_bundle": target_maps_bundle,
    }


def _find_coincidences(
    h1_cands: list[dict[str, Any]],
    l1_cands: list[dict[str, Any]],
    max_dt_s: float = 0.10,
) -> list[dict[str, Any]]:
    """Find H1-L1 candidate pairs whose peak times or intervals fall within max_dt_s."""
    pairs: list[dict[str, Any]] = []
    for c1 in h1_cands:
        for c2 in l1_cands:
            dt_peak = c1["gps_peak"] - c2["gps_peak"]
            # Interval gap (0.0 if intervals overlap)
            gap = max(0.0, max(c1["gps_start"], c2["gps_start"]) - min(c1["gps_end"], c2["gps_end"]))
            if abs(dt_peak) <= max_dt_s or gap <= max_dt_s:
                pairs.append(
                    {
                        "h1_gps_peak": c1["gps_peak"],
                        "l1_gps_peak": c2["gps_peak"],
                        "dt_h1_minus_l1_s": float(dt_peak),
                        "interval_gap_s": float(gap),
                        "h1_score": c1["anomaly_score"],
                        "l1_score": c2["anomaly_score"],
                        "h1_freq_range_hz": [c1["f_min_hz"], c1["f_max_hz"]],
                        "l1_freq_range_hz": [c2["f_min_hz"], c2["f_max_hz"]],
                        "within_10ms_light_travel": bool(abs(dt_peak) <= PHYSICAL_LIGHT_TRAVEL_S),
                    }
                )
    return pairs


def _time_shift_background_coincidences(
    h1_val_cands: list[dict[str, Any]],
    l1_val_cands: list[dict[str, Any]],
    max_dt_s: float = 0.10,
    n_shifts: int = 200,
) -> dict[str, Any]:
    """Estimate accidental H1-L1 coincidence rate via unphysical time shifts on background."""
    # Shifts of +-4s, +-8s, ..., +-400s (well outside 10 ms physical propagation delay)
    shifts = [k * 4.0 for k in range(-n_shifts // 2, n_shifts // 2 + 1) if k != 0][:n_shifts]
    counts: list[int] = []
    for dt_shift in shifts:
        shifted_l1 = [
            {
                **c,
                "gps_start": c["gps_start"] + dt_shift,
                "gps_end": c["gps_end"] + dt_shift,
                "gps_peak": c["gps_peak"] + dt_shift,
            }
            for c in l1_val_cands
        ]
        coinc = _find_coincidences(h1_val_cands, shifted_l1, max_dt_s=max_dt_s)
        counts.append(len(coinc))

    total_coinc = int(sum(counts))
    shifts_with_coinc = int(sum(1 for c in counts if c > 0))
    return {
        "n_time_shifts": len(shifts),
        "coincidence_window_s": max_dt_s,
        "n_h1_bg_candidates": len(h1_val_cands),
        "n_l1_bg_candidates": len(l1_val_cands),
        "total_accidental_coincidences": total_coinc,
        "mean_coincidences_per_shift": float(np.mean(counts)) if counts else 0.0,
        "fraction_shifts_with_ge1_coinc": float(shifts_with_coinc / len(shifts)) if shifts else 0.0,
    }


def _evaluate_against_gw150914(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    """Post-freeze evaluation comparing blind candidates to published GW150914 time."""
    if not candidates:
        return {"recovered": False, "closest_candidate": None, "timing_error_s": None}

    best = min(candidates, key=lambda c: abs(c["gps_peak"] - PUBLISHED_GW150914_GPS))
    err_peak = best["gps_peak"] - PUBLISHED_GW150914_GPS
    contains_gps = (best["gps_start"] <= PUBLISHED_GW150914_GPS <= best["gps_end"])
    near_gps = abs(err_peak) <= 0.25 or contains_gps
    f_lo, f_hi = PUBLISHED_GW150914_FREQ_BAND_HZ
    freq_overlaps = (best["f_max_hz"] >= f_lo) and (best["f_min_hz"] <= f_hi)
    return {
        "recovered": bool(near_gps and freq_overlaps),
        "contains_published_gps": bool(contains_gps),
        "timing_error_s": float(err_peak),
        "freq_overlaps_35_250hz": bool(freq_overlaps),
        "closest_candidate": best,
    }


def _generate_figures(
    rec_h1_32: GWOSCStrainRecord,
    rec_l1_32: GWOSCStrainRecord,
    res_min_h1: dict[str, Any],
    res_min_l1: dict[str, Any],
    res_adp_h1: dict[str, Any],
    res_adp_l1: dict[str, Any],
    coinc_adp: list[dict[str, Any]],
    ts_adp: dict[str, Any],
    fig_dir: Path,
) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []

    # Fig 01: Raw strain for H1 and L1 (32-s target file, NO event timestamp marker)
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 5.2), sharex=True)
    t_h1 = np.arange(rec_h1_32.strain.size) / rec_h1_32.sample_rate
    t_l1 = np.arange(rec_l1_32.strain.size) / rec_l1_32.sample_rate
    axes[0].plot(t_h1, rec_h1_32.strain * 1e18, lw=0.45, color="#1f4e79")
    axes[0].set_title(f"LIGO Hanford (H1) Raw Strain — GPS Start {rec_h1_32.gps_start:.1f} (32 s, 4096 Hz)")
    axes[0].set_ylabel("Strain (x 10^-18)")
    axes[1].plot(t_l1, rec_l1_32.strain * 1e18, lw=0.45, color="#2f855a")
    axes[1].set_title(f"LIGO Livingston (L1) Raw Strain — GPS Start {rec_l1_32.gps_start:.1f} (32 s, 4096 Hz)")
    axes[1].set_ylabel("Strain (x 10^-18)")
    axes[1].set_xlabel(f"Time Since GPS {rec_h1_32.gps_start:.1f} (s)")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig01_raw_strain_h1_l1.png")
    plt.close(fig)
    generated.append("fig01_raw_strain_h1_l1.png")

    # Fig 02: Processed (whitened & 30-380 Hz bandpassed) strain across the 32-s target windows
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 5.2), sharex=True)
    for ax, bundle, label, color in [
        (axes[0], res_adp_h1["target_maps_bundle"], "H1 Adapted (Whitened + 30-380 Hz)", "#1f4e79"),
        (axes[1], res_adp_l1["target_maps_bundle"], "L1 Adapted (Whitened + 30-380 Hz)", "#2f855a"),
    ]:
        for seg, obs, _det in bundle:
            t_rel = (seg.gps_start - rec_h1_32.gps_start) + np.arange(obs.signal.size) / obs.sample_rate
            ax.plot(t_rel, obs.signal, lw=0.5, color=color)
        ax.set_title(f"{label} Strain (1024 Hz, Independent Background PSD)")
        ax.set_ylabel("Whitened Strain (sigma)")
    axes[1].set_xlabel(f"Time Since GPS {rec_h1_32.gps_start:.1f} (s)")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig02_processed_strain_h1_l1.png")
    plt.close(fig)
    generated.append("fig02_processed_strain_h1_l1.png")

    # Select highest-scoring blind window for H1 and L1 in adapted mode (NO event-time peeking!)
    h1_top_bundle = max(res_adp_h1["target_maps_bundle"], key=lambda item: item[2].max_anomaly_score)
    l1_top_bundle = max(res_adp_l1["target_maps_bundle"], key=lambda item: item[2].max_anomaly_score)

    # Fig 03: Spectrograms and excess-power maps (top blind window per detector, NO event marker)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 5.8), sharey=True)
    for row, (name, (seg, _obs, det)) in enumerate([("H1", h1_top_bundle), ("L1", l1_top_bundle)]):
        assert det.tf_maps is not None
        maps = det.tf_maps
        p_db = 10.0 * np.log10(np.maximum(maps.power_spectrogram, 1e-18))
        im0 = axes[row, 0].pcolormesh(maps.times, maps.freqs, p_db, shading="auto", cmap="viridis")
        axes[row, 0].set_title(f"{name} Highest-Score Blind Window (GPS {seg.gps_start:.1f}) — STFT Power (dB)")
        axes[row, 0].set_ylabel("Frequency (Hz)")
        fig.colorbar(im0, ax=axes[row, 0], label="Power (dB)", pad=0.015)

        im1 = axes[row, 1].pcolormesh(
            maps.times, maps.freqs, np.clip(maps.normalized_excess_map, 0.0, 20.0),
            shading="auto", cmap="magma", vmin=0.0, vmax=20.0,
        )
        axes[row, 1].set_title(
            f"{name} Excess Power Z[f,t] (max score={det.max_anomaly_score:.1f}, thr={det.threshold_used:.1f})"
        )
        fig.colorbar(im1, ax=axes[row, 1], label="Z[f,t]", pad=0.015)
    axes[1, 0].set_xlabel("Time Within 4-s Window (s)")
    axes[1, 1].set_xlabel("Time Within 4-s Window (s)")
    fig.suptitle("Blind Highest-Score Windows — Adapted LIGO Pipeline + Frozen BSD v0.1 Detector", y=0.995)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig03_spectrograms_and_excess_maps.png")
    plt.close(fig)
    generated.append("fig03_spectrograms_and_excess_maps.png")

    # Fig 04: Anomaly score timelines across all 32-s target windows (both pipelines)
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 5.8), sharex=True)
    for ax, (title, h1_r, l1_r) in zip(
        axes,
        [
            ("Pipeline A: Frozen v0.1 Minimal (Scale + Decimate Only)", res_min_h1, res_min_l1),
            ("Pipeline B: Adapted LIGO (Background Whitening + 30-380 Hz Bandpass) + Frozen v0.1", res_adp_h1, res_adp_l1),
        ],
    ):
        t_h = [0.5 * (w["gps_start"] + w["gps_end"]) - rec_h1_32.gps_start for w in h1_r["target_windows"]]
        s_h = [w["max_anomaly_score"] for w in h1_r["target_windows"]]
        t_l = [0.5 * (w["gps_start"] + w["gps_end"]) - rec_l1_32.gps_start for w in l1_r["target_windows"]]
        s_l = [w["max_anomaly_score"] for w in l1_r["target_windows"]]

        ax.plot(t_h, s_h, "-o", color="#1f4e79", label=f"H1 max score (1% thr={h1_r['calibration'].threshold:.1f})")
        ax.axhline(h1_r["calibration"].threshold, color="#1f4e79", ls="--", lw=1.2, alpha=0.8)
        ax.plot(t_l, s_l, "-s", color="#2f855a", label=f"L1 max score (1% thr={l1_r['calibration'].threshold:.1f})")
        ax.axhline(l1_r["calibration"].threshold, color="#2f855a", ls=":", lw=1.4, alpha=0.8)
        ax.axvline(
            PUBLISHED_GW150914_GPS - rec_h1_32.gps_start,
            color="#d62728",
            ls="-.",
            lw=1.3,
            label="Published GW150914 GPS (eval reference only)",
        )
        ax.set_title(title)
        ax.set_ylabel("Max Cluster Score")
        ax.legend(loc="upper right", fontsize=8)
    axes[1].set_xlabel(f"Window Center Time Since GPS {rec_h1_32.gps_start:.1f} (s)")
    fig.suptitle("Anomaly Score Timelines Across the 32-Second Target Segment (Final Evaluation Overlay)", y=0.995)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig04_anomaly_score_timelines.png")
    plt.close(fig)
    generated.append("fig04_anomaly_score_timelines.png")

    # Fig 05: Candidate overlays + published GW150914 timestamp (final evaluation figure)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True)
    for ax, name, (seg, _obs, det) in zip(axes, ["H1", "L1"], [h1_top_bundle, l1_top_bundle]):
        assert det.tf_maps is not None
        maps = det.tf_maps
        im = ax.pcolormesh(
            maps.times, maps.freqs, np.clip(maps.normalized_excess_map, 0.0, 20.0),
            shading="auto", cmap="inferno", vmin=0.0, vmax=20.0,
        )
        fig.colorbar(im, ax=ax, label="Z[f,t]", pad=0.015)
        for idx, cand in enumerate(det.candidates):
            ax.add_patch(
                patches.Rectangle(
                    (cand.t_start, cand.f_min),
                    max(0.05, cand.t_end - cand.t_start),
                    max(8.0, cand.f_max - cand.f_min),
                    lw=1.8,
                    edgecolor="#00e5ff",
                    facecolor="none",
                    label="Blind candidate" if idx == 0 else None,
                )
            )
        t_gw_local = PUBLISHED_GW150914_GPS - seg.gps_start
        if 0.0 <= t_gw_local <= seg.duration_s:
            ax.axvline(t_gw_local, color="#ffffff", ls="--", lw=1.5, label="Published GW150914 GPS")
        ax.set_title(f"{name} | GPS {seg.gps_start:.1f} | score={det.max_anomaly_score:.1f} (thr={det.threshold_used:.1f})")
        ax.set_xlabel("Time Within 4-s Window (s)")
        ax.set_ylim(0, 420)
        ax.legend(loc="upper right", framealpha=0.8)
    axes[0].set_ylabel("Frequency (Hz)")
    fig.suptitle("Final Evaluation: Blind Candidate Overlays vs. Published GW150914 Timestamp", y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig05_candidate_overlays_eval.png")
    plt.close(fig)
    generated.append("fig05_candidate_overlays_eval.png")

    # Fig 06: Calibration null distributions and independent background validation FAR
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for name, r_adp, color in [("H1", res_adp_h1, "#1f4e79"), ("L1", res_adp_l1, "#2f855a")]:
        cal = r_adp["calibration"]
        scores = np.clip(cal.sorted_null_scores, 0, np.quantile(cal.sorted_null_scores, 0.995) * 1.25)
        axes[0].hist(scores, bins=35, alpha=0.45, density=True, color=color, label=f"{name} cal null (N={cal.n_calibration_samples})")
        axes[0].axvline(cal.threshold, color=color, lw=1.8, ls="--", label=f"{name} 1% thr={cal.threshold:.1f}")
    axes[0].set_title("Adapted LIGO Calibration Null Distributions")
    axes[0].set_xlabel("Max Cluster Anomaly Score")
    axes[0].set_ylabel("Density")
    axes[0].legend(fontsize=8)

    far_items = [
        ("H1 Minimal", res_min_h1),
        ("L1 Minimal", res_min_l1),
        ("H1 Adapted", res_adp_h1),
        ("L1 Adapted", res_adp_l1),
    ]
    x = np.arange(len(far_items))
    rates = [item[1]["val_far"] * 100.0 for item in far_items]
    lo = [(item[1]["val_far"] - item[1]["val_far_cp_95ci"][0]) * 100.0 for item in far_items]
    hi = [(item[1]["val_far_cp_95ci"][1] - item[1]["val_far"]) * 100.0 for item in far_items]
    colors = ["#2b6cb0" if item[1]["val_statistically_credible"] else "#c53030" for item in far_items]
    bars = axes[1].bar(x, rates, color=colors, alpha=0.85, width=0.55)
    axes[1].errorbar(x, rates, yerr=[lo, hi], fmt="none", ecolor="#111111", capsize=4)
    axes[1].axhline(1.0, color="#111111", ls="--", lw=1.3, label="Target FAP = 1.0%")
    for bar, rate, h_err, (_lbl, r) in zip(bars, rates, hi, far_items):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + h_err + 0.35,
            f"{rate:.2f}%\n({r['val_fp']}/{r['n_val_segments']})",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([it[0] for it in far_items], rotation=12)
    axes[1].set_ylabel("Independent Background FAR (%)")
    axes[1].set_ylim(0, max(max(r + h for r, h in zip(rates, hi)) + 2.5, 6.0))
    axes[1].set_title("Independent LIGO Background Validation FAR (95% Exact CIs)")
    axes[1].legend(loc="upper right")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig06_calibration_and_background_far.png")
    plt.close(fig)
    generated.append("fig06_calibration_and_background_far.png")

    # Fig 07: Cross-detector candidate timing comparison
    fig, ax = plt.subplots(figsize=(9.5, 4.0))
    h1_c = res_adp_h1["target_candidates_gps"]
    l1_c = res_adp_l1["target_candidates_gps"]
    if h1_c:
        ax.scatter(
            [c["gps_peak"] - rec_h1_32.gps_start for c in h1_c],
            [1.0] * len(h1_c),
            s=[min(300, max(40, c["anomaly_score"] * 0.6)) for c in h1_c],
            color="#1f4e79",
            alpha=0.85,
            label=f"H1 blind candidates (n={len(h1_c)})",
        )
    if l1_c:
        ax.scatter(
            [c["gps_peak"] - rec_l1_32.gps_start for c in l1_c],
            [0.0] * len(l1_c),
            s=[min(300, max(40, c["anomaly_score"] * 0.6)) for c in l1_c],
            color="#2f855a",
            alpha=0.85,
            label=f"L1 blind candidates (n={len(l1_c)})",
        )
    ax.axvline(
        PUBLISHED_GW150914_GPS - rec_h1_32.gps_start,
        color="#d62728",
        ls="--",
        lw=1.4,
        label="Published GW150914 GPS (eval reference)",
    )
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["L1 (Livingston)", "H1 (Hanford)"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlim(0.0, 32.0)
    ax.set_xlabel(f"Candidate Peak Time Since GPS {rec_h1_32.gps_start:.1f} (s)")
    ax.set_title(
        f"Cross-Detector Blind Candidate Timing (Adapted Pipeline)\n"
        f"Coincidences (|dt|<=0.1s): {len(coinc_adp)} | "
        f"Bg Time-Shift Rate: {ts_adp['mean_coincidences_per_shift']:.3f}/shift"
    )
    ax.legend(loc="upper right")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _add_footer(fig)
    fig.savefig(fig_dir / "fig07_cross_detector_timing.png")
    plt.close(fig)
    generated.append("fig07_cross_detector_timing.png")

    return generated


def _write_report(results: dict[str, Any], out_path: Path) -> Path:
    prov = results["provenance"]
    min_h1 = results["pipelines"]["frozen_v01_minimal"]["H1"]
    min_l1 = results["pipelines"]["frozen_v01_minimal"]["L1"]
    adp_h1 = results["pipelines"]["adapted_ligo_whitened"]["H1"]
    adp_l1 = results["pipelines"]["adapted_ligo_whitened"]["L1"]
    coinc_min = results["coincidence"]["frozen_v01_minimal"]
    coinc_adp = results["coincidence"]["adapted_ligo_whitened"]
    eval_min = results["gw150914_evaluation"]["frozen_v01_minimal"]
    eval_adp = results["gw150914_evaluation"]["adapted_ligo_whitened"]

    def _fmt_cand_table(cands: list[dict[str, Any]]) -> list[str]:
        if not cands:
            return ["*No candidates exceeded the calibrated 1% FAP threshold.*", ""]
        rows = [
            "| Detector | GPS Peak (s) | Interval [Start, End] (s) | Duration (s) | Freq Range (Hz) | Peak Freq (Hz) | Score | Emp. p-value |",
            "| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
        for c in sorted(cands, key=lambda x: x["anomaly_score"], reverse=True):
            p_str = (
                f"<= {c['empirical_p_value']:.5f}"
                if c["p_value_is_upper_bound"]
                else f"{c['empirical_p_value']:.5f}"
            )
            rows.append(
                f"| {c['detector']} | `{c['gps_peak']:.3f}` | "
                f"`[{c['gps_start']:.3f}, {c['gps_end']:.3f}]` | "
                f"`{c['duration_s']:.3f}` | `[{c['f_min_hz']:.1f}, {c['f_max_hz']:.1f}]` | "
                f"`{c['peak_freq_hz']:.1f}` | `{c['anomaly_score']:.1f}` | `{p_str}` |"
            )
        rows.append("")
        return rows

    lines: list[str] = [
        "# Blind Signal Discovery (BSD) — Experiment 002: Real LIGO Data Validation (GW150914)",
        "",
        "> **Scientific Scope Statement:** This is an external methodological validation of the frozen BSD v0.1 detector on authentic public LIGO O1 strain data from GWOSC. **A statistically unusual event is not automatically a gravitational wave.** Instrumental glitches remain an alternative hypothesis, and astrophysical confirmation belongs strictly to the LIGO Scientific Collaboration and Virgo Collaboration.",
        "",
        "## 1. Objective",
        "",
        "Evaluate whether the frozen BSD v0.1 template-free anomaly detector (`src/bsd/detector.py`, unchanged from `main`) can identify statistically unusual events in real LIGO Hanford (H1) and Livingston (L1) strain measurements — including the independently confirmed GW150914 event — without using any signal template or the published event timestamp during calibration or detection.",
        "",
        "## 2. Official Data Provenance & Dataset Details",
        "",
        "All data were downloaded unchanged from the official Gravitational Wave Open Science Center (GWOSC) GWTC-1-confident GW150914 v3 release (`https://gwosc.org/events/GW150914/`). Raw files are stored under `data/gwosc/` on drive D: (excluded from Git via `.gitignore`).",
        "",
        "| Detector | Duration (s) | Sample Rate (Hz) | GPS Interval | Size (MB) | All Finite? | DQ DATA (%) | SHA-256 (prefix) |",
        "| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]
    for p in prov:
        lines.append(
            f"| {p['detector']} | {p['duration_s']} | {p['sample_rate_hz']:.0f} | "
            f"`[{p['gps_start']:.0f}, {p['gps_end']:.0f}]` | {p['size_bytes']/1e6:.2f} | "
            f"{p['all_finite']} | {p['dq_data_bit_fraction']*100:.1f}% | `{p['sha256'][:16]}...` |"
        )
    lines.extend(
        [
            "",
            "**Why 4096-second datasets were required in addition to the 32-second files:**",
            "The 32-second files contain at most 7 non-overlapping 4-second windows (with 1-second padding), which cannot statistically calibrate a 1% experiment-level false-alarm probability (minimum resolvable empirical p-value would be $1/8 = 12.5\\%$). Using the official 4096-second GWOSC release for the same observation epoch provides $>900$ independent, DQ-verified 4-second background segments per detector after strictly excluding a $\\pm 64\\text{ s}$ guard zone around the entire 32-second target file (`GPS [1126259383, 1126259543]`).",
            "",
            "## 3. Preprocessing & Preservation of Frozen BSD v0.1",
            "",
            "The BSD v0.1 detector (`TemplateFreeDetector` with `nperseg=128`, `noverlap=96`, `window='hann'`, `freq_median_kernel=11`, `time_median_kernel=15`, `bin_seed_threshold=4.0`, `min_cluster_bins=2`) was **not modified**. Because raw LIGO strain has amplitude $h \\sim 10^{-19}$ (which would underflow `_EPS = 1e-18` in squared STFT power), a separate data-adaptation layer (`src/bsd/ligo_adapter.py`) was implemented with two explicitly distinguished pipelines:",
            "",
            "1. **Pipeline A (`frozen_v01_minimal`):** Median subtraction, division by background MAD standard deviation (bringing strain to $O(1)$ without altering spectral shape), and zero-phase anti-aliased decimation `4096 Hz -> 1024 Hz` (`scipy.signal.decimate`, IIR).",
            "2. **Pipeline B (`adapted_ligo_whitened`):** `1e18` numerical scaling, frequency-domain whitening using a robust `average='median'` Welch PSD (`4.0 s` Hann windows) estimated strictly on independent background segments, zero-phase 4th-order Butterworth bandpass `30-380 Hz` (standard aLIGO O1 sensitive band per GWOSC tutorial), and zero-phase decimation `4096 Hz -> 1024 Hz`. All filtering uses `1.0 s` left/right context padding trimmed after filtering to prevent edge transients.",
            "",
            "## 4. Noise Calibration & Independent Background Validation",
            "",
            "For each detector and pipeline, clean 4-second segments outside `GPS [1126259383, 1126259543]` were partitioned into three strictly disjoint sets:",
            f"- **PSD / scale estimation:** `{adp_h1['n_psd_segments']}` segments (`240 s`)",
            f"- **Threshold calibration (target FAP = 1.0%):** `{adp_h1['n_cal_segments']}` segments (`2000 s`) via `calibrate_detector`",
            f"- **Independent background FAR validation:** `{adp_h1['n_val_segments']}` segments (H1) and `{adp_l1['n_val_segments']}` segments (L1)",
            "",
            "| Pipeline | Detector | Calibrated 1% Threshold | Independent Bg FP / N | Empirical Bg FAR | 95% Exact CI | Statistically Credible? |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
            f"| `frozen_v01_minimal` | H1 | `{min_h1['calibration']['threshold']:.2f}` | {min_h1['val_fp']} / {min_h1['n_val_segments']} | {min_h1['val_far']*100:.2f}% | [{min_h1['val_far_cp_95ci'][0]*100:.2f}%, {min_h1['val_far_cp_95ci'][1]*100:.2f}%] | {'Yes' if min_h1['val_statistically_credible'] else 'NO (FAR inflated)'} |",
            f"| `frozen_v01_minimal` | L1 | `{min_l1['calibration']['threshold']:.2f}` | {min_l1['val_fp']} / {min_l1['n_val_segments']} | {min_l1['val_far']*100:.2f}% | [{min_l1['val_far_cp_95ci'][0]*100:.2f}%, {min_l1['val_far_cp_95ci'][1]*100:.2f}%] | {'Yes' if min_l1['val_statistically_credible'] else 'NO (FAR inflated)'} |",
            f"| `adapted_ligo_whitened` | H1 | `{adp_h1['calibration']['threshold']:.2f}` | {adp_h1['val_fp']} / {adp_h1['n_val_segments']} | {adp_h1['val_far']*100:.2f}% | [{adp_h1['val_far_cp_95ci'][0]*100:.2f}%, {adp_h1['val_far_cp_95ci'][1]*100:.2f}%] | {'Yes' if adp_h1['val_statistically_credible'] else 'NO (FAR inflated)'} |",
            f"| `adapted_ligo_whitened` | L1 | `{adp_l1['calibration']['threshold']:.2f}` | {adp_l1['val_fp']} / {adp_l1['n_val_segments']} | {adp_l1['val_far']*100:.2f}% | [{adp_l1['val_far_cp_95ci'][0]*100:.2f}%, {adp_l1['val_far_cp_95ci'][1]*100:.2f}%] | {'Yes' if adp_l1['val_statistically_credible'] else 'NO (FAR inflated)'} |",
            "",
            "## 5. Blind Single-Detector Results (32-Second Target Segment)",
            "",
            "All candidates below were extracted **without** knowledge of the GW150914 timestamp across the uniform 4-second window grid (`GPS 1126259448` to `1126259476`).",
            "",
            "### 5.1 Pipeline A (`frozen_v01_minimal` — Scale + Decimate Only)",
            "",
        ]
    )
    lines.extend(_fmt_cand_table(min_h1["target_candidates_gps"] + min_l1["target_candidates_gps"]))
    lines.extend(
        [
            "### 5.2 Pipeline B (`adapted_ligo_whitened` — Background Whitening + 30-380 Hz Bandpass)",
            "",
        ]
    )
    lines.extend(_fmt_cand_table(adp_h1["target_candidates_gps"] + adp_l1["target_candidates_gps"]))

    lines.extend(
        [
            "## 6. Two-Detector Comparison (H1 vs. L1) & Time-Shift Controls",
            "",
            "The maximum physical gravitational-wave light-travel delay between Hanford (H1) and Livingston (L1) is approximately `10 ms` (`0.010 s`). However, the frozen BSD v0.1 STFT has a frame step of `dt = 31.25 ms` (`nperseg=128, noverlap=96` at `1024 Hz`). **Because `dt = 31.25 ms > 10 ms`, BSD v0.1 cannot resolve sub-10 ms inter-site propagation delays directly.** We therefore test temporal coincidence within `|dt| <= 0.10 s` (~3 STFT frames) and calibrate accidental coincidence rates using 200 unphysical time shifts on the independent background validation set.",
            "",
            f"- **Pipeline A (`frozen_v01_minimal`):** `{len(coinc_min['target_coincidences'])}` target H1-L1 coincidences within `0.10 s`; background time-shift rate = `{coinc_min['time_shift_control']['mean_coincidences_per_shift']:.4f}` coincidences/shift (`{coinc_min['time_shift_control']['total_accidental_coincidences']}` across `{coinc_min['time_shift_control']['n_time_shifts']}` shifts).",
            f"- **Pipeline B (`adapted_ligo_whitened`):** `{len(coinc_adp['target_coincidences'])}` target H1-L1 coincidences within `0.10 s`; background time-shift rate = `{coinc_adp['time_shift_control']['mean_coincidences_per_shift']:.4f}` coincidences/shift (`{coinc_adp['time_shift_control']['total_accidental_coincidences']}` across `{coinc_adp['time_shift_control']['n_time_shifts']}` shifts).",
            "",
        ]
    )
    if coinc_adp["target_coincidences"]:
        lines.append("| H1 GPS Peak | L1 GPS Peak | dt (H1 - L1, s) | Interval Gap (s) | H1 Score | L1 Score | Within 10 ms? |")
        lines.append("| :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        for cp in coinc_adp["target_coincidences"]:
            lines.append(
                f"| `{cp['h1_gps_peak']:.4f}` | `{cp['l1_gps_peak']:.4f}` | `{cp['dt_h1_minus_l1_s']:+.4f}` | "
                f"`{cp['interval_gap_s']:.4f}` | `{cp['h1_score']:.1f}` | `{cp['l1_score']:.1f}` | "
                f"{'Yes' if cp['within_10ms_light_travel'] else 'No (STFT-limited)'} |"
            )
        lines.append("")

    lines.extend(
        [
            "## 7. Independent Evaluation Against Published GW150914 Reference",
            "",
            f"Published GWOSC reference (consulted only after all candidates were frozen): `GPS = {PUBLISHED_GW150914_GPS:.3f}`, frequency sweep `~35-250 Hz`.",
            "",
            "| Pipeline | Detector | Blind Recovery of GW150914? | Timing Error (Peak - Ref, s) | Contains Ref GPS? | Freq Overlaps 35-250 Hz? |",
            "| :--- | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for p_name, ev_pair in [("frozen_v01_minimal", eval_min), ("adapted_ligo_whitened", eval_adp)]:
        for det_name in ("H1", "L1"):
            ev = ev_pair[det_name]
            t_err = f"{ev['timing_error_s']:+.4f}" if ev["timing_error_s"] is not None else "N/A"
            lines.append(
                f"| `{p_name}` | {det_name} | **{'YES' if ev['recovered'] else 'NO'}** | `{t_err}` | "
                f"{ev.get('contains_published_gps', False)} | {ev.get('freq_overlaps_35_250hz', False)} |"
            )

    lines.extend(
        [
            "",
            "### Epistemic Distinction",
            "1. **Blind anomalies:** Candidates in Section 5 were produced solely by threshold exceedances against independent LIGO background noise, with zero access to the GW150914 timestamp or waveform.",
            "2. **Post-hoc evaluation:** Overlay of `GPS = 1126259462.422` in Figures 4, 5, and 7 was applied strictly after detection to measure localization accuracy.",
            "3. **Astrophysical confirmation:** Independent astrophysical identification of GW150914 as a binary black hole merger relies on the LIGO/Virgo collaboration's matched-filter, waveform-consistency, and detector-characterization analyses, not on this single-resolution excess-power prototype.",
            "",
            "## 8. Limitations & Failure Analysis",
            "",
            "1. **Failure of unadapted `frozen_v01_minimal`:** Without external LIGO PSD whitening and `30-380 Hz` bandpass filtering, seismic noise below `25 Hz` (6-8 orders of magnitude stronger than the sensitive band) leaks across the `125 ms` Hann STFT window and swamps the astrophysical signal.",
            "2. **Real LIGO non-Gaussian glitches:** Real interferometer background contains instrumental glitches (heavy tails), causing finite-sample threshold variance and potential FAR drift over hour-long timescales compared to stationary synthetic Gaussian noise.",
            "3. **STFT timing resolution (`31.25 ms`):** Coarser than the `10 ms` H1-L1 light-travel bound, preventing sub-10 ms propagation-delay verification without finer time resolution or cross-correlation.",
            "4. **Finite calibration tail (`N_cal = 500`):** Minimum resolvable empirical p-value is `1 / 501 = 0.00200`; scores exceeding the maximum background sample are strictly reported as upper bounds (`<= 0.00200`).",
            "",
            "## 9. Generated Figures",
            "",
        ]
    )
    for fn in results["figures"]:
        lines.append(f"- `figures/{fn}`")

    lines.extend(
        [
            "",
            "## 10. Reproducibility Commands",
            "",
            "```powershell",
            ".\\.venv\\Scripts\\python.exe experiments\\download_gwosc_data.py",
            ".\\.venv\\Scripts\\python.exe -m pytest -v",
            ".\\.venv\\Scripts\\python.exe experiments\\run_ligo_experiment_002.py",
            "```",
            "",
            "## 11. Recommendations for Experiment 003",
            "",
            "1. **Multi-resolution Q-transform / finer time step:** Reduce STFT hop size to `<= 2 ms` (or add band-limited H1-L1 cross-correlation on coincident candidates) to resolve the `10 ms` physical propagation window.",
            "2. **Glitch-veto diagnostics:** Incorporate GWOSC `CBC_CAT2`/`CBC_CAT3` hardware-injection and data-quality flags alongside a calibrated kurtosis/duration veto to separate short instrumental blips from chirping transients.",
            "3. **Multi-event blind validation:** Evaluate the same frozen adapted pipeline on additional GWTC-1 events (e.g., GW151226, GW170104, GW170814) and null off-source stretches to measure empirical sensitivity across real SNR regimes.",
            "",
        ]
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def run_experiment_002() -> dict[str, Any]:
    t0 = time.perf_counter()
    data_dir = PROJECT_ROOT / "data" / "gwosc"
    out_dir = PROJECT_ROOT / "results" / "ligo_experiment_002"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    prov_path = out_dir / "data_provenance.json"
    provenance = json.loads(prov_path.read_text(encoding="utf-8")) if prov_path.is_file() else []

    print("[BSD-002] Loading official GWOSC H1 and L1 HDF5 files...")
    rec_h1_32 = load_gwosc_hdf5(data_dir / "H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5")
    rec_l1_32 = load_gwosc_hdf5(data_dir / "L-L1_GWOSC_4KHZ_R1-1126259447-32.hdf5")
    rec_h1_4096 = load_gwosc_hdf5(data_dir / "H-H1_GWOSC_4KHZ_R1-1126257415-4096.hdf5")
    rec_l1_4096 = load_gwosc_hdf5(data_dir / "L-L1_GWOSC_4KHZ_R1-1126257415-4096.hdf5")

    # Frozen BSD v0.1 detector (exact default configuration, zero modifications)
    bsd_detector = TemplateFreeDetector(DetectorConfig())

    print("[BSD-002] Running Pipeline A (frozen_v01_minimal) on H1 and L1...")
    res_min_h1 = _run_detector_pipeline_for_site("H1", rec_h1_32, rec_h1_4096, bsd_detector, mode="frozen_v01_minimal")
    res_min_l1 = _run_detector_pipeline_for_site("L1", rec_l1_32, rec_l1_4096, bsd_detector, mode="frozen_v01_minimal")

    print("[BSD-002] Running Pipeline B (adapted_ligo_whitened) on H1 and L1...")
    res_adp_h1 = _run_detector_pipeline_for_site("H1", rec_h1_32, rec_h1_4096, bsd_detector, mode="adapted_ligo_whitened")
    res_adp_l1 = _run_detector_pipeline_for_site("L1", rec_l1_32, rec_l1_4096, bsd_detector, mode="adapted_ligo_whitened")

    # Cross-detector coincidence + 200-shift background control
    coinc_min_pairs = _find_coincidences(res_min_h1["target_candidates_gps"], res_min_l1["target_candidates_gps"], max_dt_s=0.10)
    ts_min = _time_shift_background_coincidences(res_min_h1["val_candidates_gps"], res_min_l1["val_candidates_gps"], max_dt_s=0.10)

    coinc_adp_pairs = _find_coincidences(res_adp_h1["target_candidates_gps"], res_adp_l1["target_candidates_gps"], max_dt_s=0.10)
    ts_adp = _time_shift_background_coincidences(res_adp_h1["val_candidates_gps"], res_adp_l1["val_candidates_gps"], max_dt_s=0.10)

    # Post-freeze evaluation against published GW150914 timestamp
    eval_min = {
        "H1": _evaluate_against_gw150914(res_min_h1["target_candidates_gps"]),
        "L1": _evaluate_against_gw150914(res_min_l1["target_candidates_gps"]),
    }
    eval_adp = {
        "H1": _evaluate_against_gw150914(res_adp_h1["target_candidates_gps"]),
        "L1": _evaluate_against_gw150914(res_adp_l1["target_candidates_gps"]),
    }

    figures = _generate_figures(
        rec_h1_32,
        rec_l1_32,
        res_min_h1,
        res_min_l1,
        res_adp_h1,
        res_adp_l1,
        coinc_adp_pairs,
        ts_adp,
        fig_dir,
    )

    def _serialize_site(r: dict[str, Any]) -> dict[str, Any]:
        return {
            "detector": r["detector"],
            "mode": r["mode"],
            "n_psd_segments": r["n_psd_segments"],
            "n_cal_segments": r["n_cal_segments"],
            "n_val_segments": r["n_val_segments"],
            "calibration": r["calibration"].to_summary_dict(),
            "val_fp": r["val_fp"],
            "val_far": r["val_far"],
            "val_far_cp_95ci": r["val_far_cp_95ci"],
            "val_far_wilson_95ci": r["val_far_wilson_95ci"],
            "val_statistically_credible": r["val_statistically_credible"],
            "target_windows": r["target_windows"],
            "target_candidates_gps": r["target_candidates_gps"],
        }

    elapsed = time.perf_counter() - t0
    results: dict[str, Any] = {
        "metadata": {
            "experiment": "ligo_experiment_002",
            "event_reference": "GW150914",
            "published_gps_eval_only": PUBLISHED_GW150914_GPS,
            "detector_config": bsd_detector.config.to_dict(),
            "runtime_seconds": round(elapsed, 2),
        },
        "provenance": provenance,
        "pipelines": {
            "frozen_v01_minimal": {
                "H1": _serialize_site(res_min_h1),
                "L1": _serialize_site(res_min_l1),
            },
            "adapted_ligo_whitened": {
                "H1": _serialize_site(res_adp_h1),
                "L1": _serialize_site(res_adp_l1),
            },
        },
        "coincidence": {
            "frozen_v01_minimal": {
                "target_coincidences": coinc_min_pairs,
                "time_shift_control": ts_min,
            },
            "adapted_ligo_whitened": {
                "target_coincidences": coinc_adp_pairs,
                "time_shift_control": ts_adp,
            },
        },
        "gw150914_evaluation": {
            "frozen_v01_minimal": eval_min,
            "adapted_ligo_whitened": eval_adp,
        },
        "figures": figures,
    }

    (out_dir / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

    all_cands = (
        res_min_h1["target_candidates_gps"]
        + res_min_l1["target_candidates_gps"]
        + res_adp_h1["target_candidates_gps"]
        + res_adp_l1["target_candidates_gps"]
    )
    _write_csv(out_dir / "target_candidates.csv", all_cands)

    all_windows = (
        res_min_h1["target_windows"]
        + res_min_l1["target_windows"]
        + res_adp_h1["target_windows"]
        + res_adp_l1["target_windows"]
    )
    _write_csv(out_dir / "target_windows.csv", all_windows)

    report_path = _write_report(results, out_dir / "report.md")

    print(f"[BSD-002] Finished in {elapsed:.2f} s.")
    print(
        f"[BSD-002] Minimal H1: thr={res_min_h1['calibration'].threshold:.1f}, "
        f"bg_FAR={res_min_h1['val_far']*100:.2f}%, target_cands={len(res_min_h1['target_candidates_gps'])}, "
        f"GW150914_recovered={eval_min['H1']['recovered']}"
    )
    print(
        f"[BSD-002] Minimal L1: thr={res_min_l1['calibration'].threshold:.1f}, "
        f"bg_FAR={res_min_l1['val_far']*100:.2f}%, target_cands={len(res_min_l1['target_candidates_gps'])}, "
        f"GW150914_recovered={eval_min['L1']['recovered']}"
    )
    print(
        f"[BSD-002] Adapted H1: thr={res_adp_h1['calibration'].threshold:.1f}, "
        f"bg_FAR={res_adp_h1['val_far']*100:.2f}%, target_cands={len(res_adp_h1['target_candidates_gps'])}, "
        f"GW150914_recovered={eval_adp['H1']['recovered']}"
    )
    print(
        f"[BSD-002] Adapted L1: thr={res_adp_l1['calibration'].threshold:.1f}, "
        f"bg_FAR={res_adp_l1['val_far']*100:.2f}%, target_cands={len(res_adp_l1['target_candidates_gps'])}, "
        f"GW150914_recovered={eval_adp['L1']['recovered']}"
    )
    print(f"[BSD-002] Adapted H1-L1 coincidences (<=0.10s): {len(coinc_adp_pairs)} | bg time-shift mean: {ts_adp['mean_coincidences_per_shift']:.4f}")
    print(f"[BSD-002] Report: {report_path}")
    return results


if __name__ == "__main__":
    run_experiment_002()
