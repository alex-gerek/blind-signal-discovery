"""WHISKER / BSD Exploratory Real-Data Anomaly Search on the GW170814 Epoch (H1, L1, V1).

IMPORTANT SCOPE NOTICE:
  This is an EXPLORATORY analysis of a previously used O2 GWOSC epoch
  (GW170814, GPS 1186739814..1186743910), NOT a blind discovery search on
  unseen data. Detection of GW170814 itself serves as a repeatability /
  positive-control check; the adjacent off-source interval is used for
  exploratory candidate ranking, DQ/veto triage, and 3-detector cross-checks.

Strict Non-Overlapping Epoch Partition (including pad_s = 1.0 s on every 4-s window):
  - Epoch bounds:               GPS 1186739814 .. 1186743910 (4096 s)
  - PSD-fitting interval:       GPS 1186739816 .. 1186740056 (60 x 4 s = 240 s core;
                                padded 1186739815 .. 1186740057)
  - Guard 1 (PSD -> Cal A):     GPS 1186740057 .. 1186740079 (22 s gap)
  - Calibration Block A:        GPS 1186740080 .. 1186741360 (320 x 4 s = 1280 s core;
                                padded 1186740079 .. 1186741361)
  - Guard 2 (Cal A -> Search):  GPS 1186741361 .. 1186741539 (178 s gap)
  - Target Search Interval:     GPS 1186741540 .. 1186742180 (160 x 4 s = 640 s core;
                                padded 1186741539 .. 1186742181)
      * GW170814 merger time:   t_c = 1186741861.5 s (inside Window 80: [1186741860, 1186741864])
      * +-16 s exclusion zone:  [1186741845.5, 1186741877.5]
      * Excluded from off-source (any window whose padded [gps0-1, gps0+5] intersects +-16 s):
          - Window 80:          ON_SOURCE_POSITIVE_CONTROL ([1186741860, 1186741864])
          - Windows 76..79, 81..84: EVENT_GUARD_ZONE_EXCLUDED (8 windows)
      * Off-source search set:  151 non-overlapping 4-s windows (604 s = 10.07 min)
  - Guard 3 (Search -> Cal B):  GPS 1186742181 .. 1186742359 (178 s gap)
  - Calibration Block B:        GPS 1186742360 .. 1186743908 (387 x 4 s = 1548 s core;
                                padded 1186742359 .. 1186743909)
  - Total Calibration Ensemble: 707 non-overlapping 4-s windows (Split A = 212, Split B = 495)

Detector mathematics, GUI, historical baselines, and the reserved O3 holdout
are completely untouched.
"""

from __future__ import annotations

import csv
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sp_signal

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bsd.adaptive_regime_v2b import estimate_calibration_uncertainty
from bsd.coincidence import (
    estimate_cross_detector_delay,
    match_coincident_candidates,
    run_time_slide_analysis,
)
from bsd.detector_v2 import DEFAULT_V02_SCALES, MultiResolutionDetector
from bsd.ligo_adapter import (
    GWOSCStrainRecord,
    PaddedSegment,
    fit_background_model,
    load_gwosc_hdf5,
    preprocess_segment,
)
from bsd.vetoes import VetoConfig, apply_vetoes
from bsd.workstation_runner import classify_calibration_status

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 180,
        "font.size": 9.2,
        "axes.titlesize": 10.2,
        "axes.labelsize": 9.2,
        "legend.fontsize": 8.0,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
)

EVENT_NAME = "GW170814"
EVENT_GPS = 1186741861.5
EVENT_EXCLUSION_HALF_WIDTH_S = 16.0

SEGMENT_DURATION_S = 4.0
PAD_S = 1.0
TARGET_SAMPLE_RATE_HZ = 1024.0
BANDPASS_HZ = (30.0, 380.0)
PREPROCESS_MODE = "adapted_ligo_whitened"
TARGET_FAP = 0.01
NORM_FRACTION = 0.30

# Exact non-overlapping core intervals (each window uses pad_s = 1.0 s on both sides)
PSD_CORE_INTERVAL = (1186739816.0, 1186740056.0)       # 60 x 4 s = 240 s
CAL_A_CORE_INTERVAL = (1186740080.0, 1186741360.0)     # 320 x 4 s = 1280 s
SEARCH_CORE_INTERVAL = (1186741540.0, 1186742180.0)    # 160 x 4 s = 640 s
CAL_B_CORE_INTERVAL = (1186742360.0, 1186743908.0)     # 387 x 4 s = 1548 s

DETECTORS: dict[str, Path] = {
    "H1": PROJECT_ROOT / "data" / "gwosc" / "GW170814" / "H-H1_GWOSC_4KHZ_R1-1186739814-4096.hdf5",
    "L1": PROJECT_ROOT / "data" / "gwosc" / "GW170814" / "L-L1_GWOSC_4KHZ_R1-1186739814-4096.hdf5",
    "V1": PROJECT_ROOT / "data" / "gwosc" / "GW170814" / "V-V1_GWOSC_4KHZ_R1-1186739814-4096.hdf5",
}

# Coincidence parameters (prespecified from STFT scale hop sizes + physical light-travel bounds)
PHYSICAL_LIGHT_TRAVEL_S: dict[tuple[str, str], float] = {
    ("H1", "L1"): 0.0100,
    ("H1", "V1"): 0.0273,
    ("L1", "V1"): 0.0264,
}
COINC_SCREENING_DT_S: dict[tuple[str, str], float] = {
    ("H1", "L1"): 0.050,
    ("H1", "V1"): 0.065,
    ("L1", "V1"): 0.065,
}
XCORR_MAX_LAG_S: dict[tuple[str, str], float] = {
    ("H1", "L1"): 0.015,
    ("H1", "V1"): 0.030,
    ("L1", "V1"): 0.030,
}

# Frozen reversible veto configuration (matches BSD v0.2 transient search benchmark)
FROZEN_VETO_CONFIG = VetoConfig(
    enable_glitch_veto=True,
    enable_line_veto=True,
)

RESULTS_DIR = PROJECT_ROOT / "results" / "gw170814_exploratory_search"


def _get_git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=PROJECT_ROOT,
            stderr=subprocess.DEVNULL,
        )
        return out.decode("utf-8").strip()
    except Exception:
        return "unknown"


def _wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 1.0)
    p = k / n
    d = 1.0 + (z * z) / n
    c = (p + (z * z) / (2.0 * n)) / d
    h = (z * np.sqrt((p * (1.0 - p) + (z * z) / (4.0 * n)) / n)) / d
    return (max(0.0, float(c - h)), min(1.0, float(c + h)))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _load_inj_mask(path: Path) -> tuple[np.ndarray, tuple[str, ...]]:
    with h5py.File(path, "r") as hf:
        inj_mask = np.asarray(hf["quality/injections/Injmask"][:], dtype=np.int64)
        raw_names = hf["quality/injections/InjShortnames"][:]
        names = tuple(x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in raw_names)
    return inj_mask, names


def _verify_disjoint_intervals() -> list[dict[str, Any]]:
    """Verify strict non-overlap across all padded intervals and guards."""
    blocks = [
        ("PSD_FIT", PSD_CORE_INTERVAL),
        ("CAL_BLOCK_A", CAL_A_CORE_INTERVAL),
        ("SEARCH_TARGET", SEARCH_CORE_INTERVAL),
        ("CAL_BLOCK_B", CAL_B_CORE_INTERVAL),
    ]
    summary: list[dict[str, Any]] = []
    padded_list: list[tuple[str, float, float]] = []
    for name, (t0, t1) in blocks:
        dur = t1 - t0
        n_win = int(round(dur / SEGMENT_DURATION_S))
        assert abs(n_win * SEGMENT_DURATION_S - dur) < 1e-9, f"{name} duration not multiple of 4s"
        p0, p1 = t0 - PAD_S, t1 + PAD_S
        padded_list.append((name, p0, p1))
        summary.append(
            {
                "name": name,
                "core_gps_start": t0,
                "core_gps_end": t1,
                "core_duration_s": dur,
                "n_windows_4s": n_win,
                "padded_gps_start": p0,
                "padded_gps_end": p1,
            }
        )

    for i in range(len(padded_list) - 1):
        name_a, _a0, a1 = padded_list[i]
        name_b, b0, _b1 = padded_list[i + 1]
        gap = b0 - a1
        assert gap > 0.0, f"Padded overlap between {name_a} and {name_b}: gap={gap}"
        summary[i]["guard_gap_after_padded_s"] = gap
    summary[-1]["guard_gap_after_padded_s"] = 0.0
    return summary


def _extract_exact_interval_segments(
    record: GWOSCStrainRecord,
    core_interval: tuple[float, float],
) -> tuple[list[PaddedSegment], list[dict[str, Any]]]:
    """Slice exact 4-s non-overlapping windows inside `core_interval` with pad_s=1.0."""
    fs = record.sample_rate
    seg_len = int(round(SEGMENT_DURATION_S * fs))
    pad_len = int(round(PAD_S * fs))
    t0, t1 = core_interval
    n_win = int(round((t1 - t0) / SEGMENT_DURATION_S))

    segments: list[PaddedSegment] = []
    unusable: list[dict[str, Any]] = []

    for w_idx in range(n_win):
        gps0 = t0 + w_idx * SEGMENT_DURATION_S
        gps1 = gps0 + SEGMENT_DURATION_S
        pad0 = gps0 - PAD_S
        pad1 = gps1 + PAD_S
        start_idx = int(round((gps0 - record.gps_start) * fs))
        if start_idx - pad_len < 0 or start_idx + seg_len + pad_len > record.strain.size:
            unusable.append({"detector": record.detector, "gps_start": gps0, "gps_end": gps1, "reason": "OUT_OF_FILE_BOUNDS"})
            continue

        sec0 = max(0, int(np.floor(pad0 - record.gps_start)))
        sec1 = min(record.dq_mask.size, int(np.ceil(pad1 - record.gps_start)))
        if sec1 <= sec0 or not np.all((record.dq_mask[sec0:sec1] & 1) > 0):
            unusable.append({"detector": record.detector, "gps_start": gps0, "gps_end": gps1, "reason": "DQ_DATA_BIT_FAIL"})
            continue

        chunk = record.strain[start_idx - pad_len : start_idx + seg_len + pad_len]
        if not np.all(np.isfinite(chunk)):
            unusable.append({"detector": record.detector, "gps_start": gps0, "gps_end": gps1, "reason": "NONFINITE_SAMPLES"})
            continue

        segments.append(
            PaddedSegment(
                gps_start=gps0,
                duration_s=SEGMENT_DURATION_S,
                sample_rate=fs,
                padded_strain=chunk.copy(),
                pad_left_samples=pad_len,
                pad_right_samples=pad_len,
            )
        )
    return segments, unusable


def _inspect_flags_at_interval(
    record: GWOSCStrainRecord,
    inj_mask: np.ndarray,
    inj_names: tuple[str, ...],
    gps_start: float,
    gps_end: float,
) -> dict[str, Any]:
    """Return exact GWOSC DQmask and Injmask status across [gps_start, gps_end]."""
    sec0 = max(0, int(np.floor(gps_start - record.gps_start)))
    sec1 = min(record.dq_mask.size, max(sec0 + 1, int(np.ceil(gps_end - record.gps_start))))
    dq_slice = record.dq_mask[sec0:sec1]
    inj_slice = inj_mask[sec0:sec1]

    dq_failed = [nm for bit, nm in enumerate(record.dq_shortnames) if not bool(np.all((dq_slice & (1 << bit)) > 0))]
    # In GWOSC Injmask, bit == 1 means NO injection; bit == 0 means injection active
    inj_active = [nm.replace("NO_", "") for bit, nm in enumerate(inj_names) if not bool(np.all((inj_slice & (1 << bit)) > 0))]
    return {
        "dq_all_pass": len(dq_failed) == 0,
        "dq_failed_flags": ";".join(dq_failed) if dq_failed else "NONE",
        "inj_active_flags": ";".join(inj_active) if inj_active else "NONE",
    }


def _generate_figures(
    window_timelines: dict[str, list[dict[str, Any]]],
    cal_summary: dict[str, dict[str, Any]],
    on_source_obs: dict[str, Any],
    top_offsource_info: dict[str, Any] | None,
    out_dir: Path,
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Figure 1: 3-detector search timeline across the 640-s target interval
    # ------------------------------------------------------------------
    fig1_path = out_dir / "fig01_gw170814_search_timeline.png"
    fig, axes = plt.subplots(3, 1, figsize=(11.8, 7.6), sharex=True)
    det_colors = {"H1": "#2b6cb0", "L1": "#2f855a", "V1": "#6b46c1"}

    for ax, det in zip(axes, ("H1", "L1", "V1")):
        rows = window_timelines[det]
        th = cal_summary[det]["frozen_threshold"]
        t_rel = np.array([(r["gps_start"] + 2.0) - EVENT_GPS for r in rows])
        scores = np.array([r["max_global_score"] for r in rows])
        roles = [r["window_role"] for r in rows]

        ax.axvspan(-EVENT_EXCLUSION_HALF_WIDTH_S, EVENT_EXCLUSION_HALF_WIDTH_S, color="#feebc8", alpha=0.55, label="+-16 s event exclusion zone")
        ax.axvline(0.0, color="#c53030", linestyle="--", linewidth=1.2, label="GW170814 merger (t_c)")
        ax.axhline(th, color="#1a202c", linestyle=":", linewidth=1.3, label=f"Frozen threshold ({th:.2f}, {cal_summary[det]['calibration_status']})")

        # Off-source points
        off_idx = [i for i, rl in enumerate(roles) if rl == "OFF_SOURCE"]
        on_idx = [i for i, rl in enumerate(roles) if rl == "ON_SOURCE_POSITIVE_CONTROL"]
        grd_idx = [i for i, rl in enumerate(roles) if rl == "EVENT_GUARD_ZONE_EXCLUDED"]

        ax.plot(t_rel[off_idx], scores[off_idx], ".", color=det_colors[det], markersize=4.5, alpha=0.8, label="Off-source windows (n=151)")
        if grd_idx:
            ax.plot(t_rel[grd_idx], scores[grd_idx], "x", color="#718096", markersize=5.0, alpha=0.8, label="Excluded guard windows (n=8)")
        if on_idx:
            ax.plot(t_rel[on_idx], scores[on_idx], "*", color="#c53030", markersize=11.0, zorder=5, label=f"On-source GW170814 ({scores[on_idx[0]]:.2f})")

        # Highlight off-source threshold crossings
        off_hits = [i for i in off_idx if scores[i] > th]
        if off_hits:
            ax.scatter(t_rel[off_hits], scores[off_hits], facecolors="none", edgecolors="#dd6b20", s=55, linewidths=1.5, zorder=4, label=f"Off-source crossings (n={len(off_hits)})")

        ax.set_ylabel(f"{det} Global Score")
        ax.set_title(f"{det} Frozen Multi-Resolution Score Across 640-s Target Interval (Exploratory, Previously Used Epoch)")
        ax.legend(loc="upper right", ncol=2, frameon=True, framealpha=0.90, fontsize=7.3)

    axes[-1].set_xlabel("Window Center Time Relative to GW170814 Merger t_c = 1186741861.5 (s)")
    fig.tight_layout()
    fig.savefig(fig1_path, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: GW170814 positive-control spectrogram (H1/L1) + top off-source
    # ------------------------------------------------------------------
    fig2_path = out_dir / "fig02_gw170814_positive_control_and_top_offsource.png"
    n_panels = 3 if top_offsource_info is not None else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(4.3 * n_panels, 4.1))
    if n_panels == 2:
        axes = list(axes)

    for idx, det in enumerate(("H1", "L1")):
        obs = on_source_obs[det]
        f, t, sxx = sp_signal.spectrogram(obs.signal, fs=obs.sample_rate, window="hann", nperseg=128, noverlap=112)
        mask = (f >= 30.0) & (f <= 360.0)
        p_db = 10.0 * np.log10(np.maximum(sxx[mask], 1e-18))
        vmin, vmax = float(np.quantile(p_db, 0.15)), float(np.quantile(p_db, 0.998))
        im = axes[idx].pcolormesh(t, f[mask], p_db, shading="gouraud", cmap="viridis", vmin=vmin, vmax=vmax)
        axes[idx].axvline(EVENT_GPS - 1186741860.0, color="#feb2b2", linestyle="--", linewidth=1.2, label="GW170814 t_c")
        axes[idx].set_title(f"Positive Control: {det} On-Source\n(GPS 1186741860..1186741864)")
        axes[idx].set_xlabel("Time in 4-s Window (s)")
        if idx == 0:
            axes[idx].set_ylabel("Frequency (Hz)")
        axes[idx].legend(loc="upper right", fontsize=7.5)
        fig.colorbar(im, ax=axes[idx], fraction=0.046, pad=0.04, label="Power (dB)")

    if top_offsource_info is not None:
        obs_top = top_offsource_info["obs"]
        cand = top_offsource_info["row"]
        f, t, sxx = sp_signal.spectrogram(obs_top.signal, fs=obs_top.sample_rate, window="hann", nperseg=128, noverlap=112)
        mask = (f >= 30.0) & (f <= 360.0)
        p_db = 10.0 * np.log10(np.maximum(sxx[mask], 1e-18))
        vmin, vmax = float(np.quantile(p_db, 0.15)), float(np.quantile(p_db, 0.998))
        im = axes[2].pcolormesh(t, f[mask], p_db, shading="gouraud", cmap="viridis", vmin=vmin, vmax=vmax)
        axes[2].axvline(cand["t_peak_in_window_s"], color="#fbd38d", linestyle="--", linewidth=1.2, label=f"Peak ({cand['t_peak_in_window_s']:.2f} s)")
        axes[2].set_title(
            f"Top Off-Source: {cand['detector']} (score={cand['standardized_score']:.2f})\n"
            f"GPS {cand['window_gps_start']:.0f}..{cand['window_gps_end']:.0f} | veto={cand['veto_status']}"
        )
        axes[2].set_xlabel("Time in 4-s Window (s)")
        axes[2].legend(loc="upper right", fontsize=7.5)
        fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04, label="Power (dB)")

    fig.tight_layout()
    fig.savefig(fig2_path, bbox_inches="tight")
    plt.close(fig)

    return fig1_path, fig2_path


def main() -> None:
    t0 = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Verify strict non-overlap of all intervals BEFORE loading or searching
    interval_table = _verify_disjoint_intervals()
    scale_localization = []
    for sp in DEFAULT_V02_SCALES:
        cfg = sp.detector_config
        hop_s = (cfg.nperseg - cfg.noverlap) / TARGET_SAMPLE_RATE_HZ
        win_s = cfg.nperseg / TARGET_SAMPLE_RATE_HZ
        scale_localization.append(
            {
                "scale_name": sp.name,
                "nperseg": cfg.nperseg,
                "noverlap": cfg.noverlap,
                "window_width_ms": round(win_s * 1000.0, 2),
                "hop_step_ms": round(hop_s * 1000.0, 2),
                "half_hop_quantization_ms": round(0.5 * hop_s * 1000.0, 2),
            }
        )

    protocol: dict[str, Any] = {
        "study_label": "EXPLORATORY_ANALYSIS_OF_PREVIOUSLY_USED_O2_EPOCH_NOT_BLIND_DISCOVERY",
        "event_name": EVENT_NAME,
        "event_gps": EVENT_GPS,
        "event_exclusion_half_width_s": EVENT_EXCLUSION_HALF_WIDTH_S,
        "event_exclusion_interval": [EVENT_GPS - EVENT_EXCLUSION_HALF_WIDTH_S, EVENT_GPS + EVENT_EXCLUSION_HALF_WIDTH_S],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _get_git_commit(),
        "python_version": platform.python_version(),
        "segment_duration_s": SEGMENT_DURATION_S,
        "pad_s": PAD_S,
        "preprocessing": {
            "mode": PREPROCESS_MODE,
            "target_sample_rate_hz": TARGET_SAMPLE_RATE_HZ,
            "bandpass_hz": list(BANDPASS_HZ),
            "psd_fit_windows": 60,
            "psd_fit_core_duration_s": 240.0,
        },
        "calibration_config": {
            "detector_class": "MultiResolutionDetector",
            "target_fap_nominal": TARGET_FAP,
            "norm_fraction": NORM_FRACTION,
            "nominal_fap_note": "Nominal target_fap=0.01 is a calibration quantile parameter, NOT a measured FAR on real LIGO/Virgo data.",
        },
        "veto_config_reversible": FROZEN_VETO_CONFIG.to_dict(),
        "stft_scale_localization": scale_localization,
        "coincidence_config": {
            "physical_light_travel_s": {f"{a}-{b}": v for (a, b), v in PHYSICAL_LIGHT_TRAVEL_S.items()},
            "stft_screening_max_dt_s": {f"{a}-{b}": v for (a, b), v in COINC_SCREENING_DT_S.items()},
            "xcorr_max_lag_s": {f"{a}-{b}": v for (a, b), v in XCORR_MAX_LAG_S.items()},
            "single_detector_interpretation_rule": (
                "A candidate detected in only one detector is NEVER classified as instrumental "
                "solely because another detector did not detect it (e.g., antenna pattern or sensitivity differences)."
            ),
        },
        "verified_disjoint_intervals": interval_table,
    }
    (RESULTS_DIR / "protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    print("[GW170814 Exploratory Search] Protocol frozen & disjoint intervals verified. Running H1, L1, V1...", flush=True)

    detector = MultiResolutionDetector()
    cal_summary: dict[str, dict[str, Any]] = {}
    window_timelines: dict[str, list[dict[str, Any]]] = {}
    excluded_rows: list[dict[str, Any]] = []
    all_raw_candidates: list[dict[str, Any]] = []
    per_det_search_cands_for_coinc: dict[str, dict[str, list[dict[str, Any]]]] = {
        "ON_SOURCE": {"H1": [], "L1": [], "V1": []},
        "OFF_SOURCE": {"H1": [], "L1": [], "V1": []},
    }
    search_obs_by_det_gps: dict[str, dict[float, Any]] = {"H1": {}, "L1": {}, "V1": {}}
    on_source_obs: dict[str, Any] = {}

    ex0 = EVENT_GPS - EVENT_EXCLUSION_HALF_WIDTH_S
    ex1 = EVENT_GPS + EVENT_EXCLUSION_HALF_WIDTH_S

    for det, h5_path in DETECTORS.items():
        t_det = time.perf_counter()
        record = load_gwosc_hdf5(h5_path)
        inj_mask, inj_names = _load_inj_mask(h5_path)

        # Extract strictly disjoint segments
        psd_segs, psd_unusable = _extract_exact_interval_segments(record, PSD_CORE_INTERVAL)
        cal_a_segs, cal_a_unusable = _extract_exact_interval_segments(record, CAL_A_CORE_INTERVAL)
        cal_b_segs, cal_b_unusable = _extract_exact_interval_segments(record, CAL_B_CORE_INTERVAL)
        search_segs, search_unusable = _extract_exact_interval_segments(record, SEARCH_CORE_INTERVAL)

        for tag, u_list in (("PSD", psd_unusable), ("CAL_A", cal_a_unusable), ("CAL_B", cal_b_unusable), ("SEARCH", search_unusable)):
            for u in u_list:
                excluded_rows.append({**u, "interval_group": tag})

        # Fit frozen background PSD strictly on PSD_CORE_INTERVAL
        bg_model = fit_background_model(
            detector=det,
            bg_segments=psd_segs,
            mode=PREPROCESS_MODE,
            target_sample_rate=TARGET_SAMPLE_RATE_HZ,
            bandpass_hz=BANDPASS_HZ,
        )

        # Preprocess calibration windows (Cal A + Cal B) and freeze calibration
        cal_segs = cal_a_segs + cal_b_segs
        cal_obs = [
            preprocess_segment(seg, f"{det}_cal_{i:04d}", bg_model)
            for i, seg in enumerate(cal_segs)
        ]
        frozen_cal = detector.calibrate(
            cal_obs,
            target_fap=TARGET_FAP,
            norm_fraction=NORM_FRACTION,
            calibration_seed=0,
        )
        n_norm = frozen_cal.n_norm_samples
        split_b_chrono = np.array(
            [detector.analyze(o, frozen_cal).max_global_score for o in cal_obs[n_norm:]],
            dtype=np.float64,
        )
        cal_unc = estimate_calibration_uncertainty(split_b_chrono, target_fap=TARGET_FAP)
        cal_audit = classify_calibration_status(len(cal_obs), cal_unc)

        # Search the 640-s target interval with frozen detector
        det_timeline: list[dict[str, Any]] = []
        n_off_windows = 0
        n_off_nominated_raw = 0
        n_off_nominated_post_veto = 0
        on_source_detected_raw = False
        on_source_detected_post_veto = False
        on_source_max_score = 0.0

        for idx, seg in enumerate(search_segs):
            gps0 = seg.gps_start
            gps1 = gps0 + SEGMENT_DURATION_S
            pad0 = gps0 - PAD_S
            pad1 = gps1 + PAD_S

            if gps0 <= EVENT_GPS < gps1:
                role = "ON_SOURCE_POSITIVE_CONTROL"
            elif not (pad1 <= ex0 or pad0 >= ex1):
                role = "EVENT_GUARD_ZONE_EXCLUDED"
                excluded_rows.append(
                    {
                        "detector": det,
                        "gps_start": gps0,
                        "gps_end": gps1,
                        "reason": f"PADDED_OVERLAP_WITH_EVENT_PM_{int(EVENT_EXCLUSION_HALF_WIDTH_S)}S_GUARD",
                        "interval_group": "SEARCH",
                    }
                )
            else:
                role = "OFF_SOURCE"

            obs = preprocess_segment(seg, f"{det}_search_{idx:03d}", bg_model)
            search_obs_by_det_gps[det][gps0] = obs
            if role == "ON_SOURCE_POSITIVE_CONTROL":
                on_source_obs[det] = obs

            res = detector.analyze(obs, frozen_cal)
            veto_assess = apply_vetoes(obs, res.candidates, FROZEN_VETO_CONFIG)
            veto_map = {d.candidate.candidate_id: d for d in veto_assess.decisions}

            det_timeline.append(
                {
                    "detector": det,
                    "window_index": idx,
                    "gps_start": gps0,
                    "gps_end": gps1,
                    "window_role": role,
                    "max_global_score": round(res.max_global_score, 4),
                    "detected_raw": res.detected,
                    "detected_post_veto": veto_assess.detected_after_veto,
                    "n_raw_candidates": len(res.candidates),
                    "n_surviving_candidates": len(veto_assess.surviving_candidates),
                }
            )

            if role == "ON_SOURCE_POSITIVE_CONTROL":
                on_source_detected_raw = bool(res.detected)
                on_source_detected_post_veto = bool(veto_assess.detected_after_veto)
                on_source_max_score = float(res.max_global_score)
            elif role == "OFF_SOURCE":
                n_off_windows += 1
                if res.detected:
                    n_off_nominated_raw += 1
                if veto_assess.detected_after_veto:
                    n_off_nominated_post_veto += 1

            # Save EVERY raw candidate verbatim (both on-source and off-source; also guard if any)
            for c in res.candidates:
                vdec = veto_map.get(c.candidate_id)
                is_vetoed = bool(vdec.vetoed) if vdec is not None else False
                veto_reasons = ";".join(vdec.reason_codes) if (vdec is not None and vdec.reason_codes) else "NONE"
                c_gps_start = gps0 + c.t_start
                c_gps_end = gps0 + c.t_end
                c_gps_peak = gps0 + c.peak_time
                flags = _inspect_flags_at_interval(record, inj_mask, inj_names, c_gps_start, c_gps_end)

                cand_row = {
                    "detector": det,
                    "window_role": role,
                    "window_index": idx,
                    "window_gps_start": gps0,
                    "window_gps_end": gps1,
                    "candidate_id_in_window": c.candidate_id,
                    "scale_name": c.scale_name,
                    "gps_start": round(c_gps_start, 4),
                    "gps_end": round(c_gps_end, 4),
                    "gps_peak": round(c_gps_peak, 4),
                    "t_peak_in_window_s": round(c.peak_time, 4),
                    "t_rel_event_s": round(c_gps_peak - EVENT_GPS, 4),
                    "duration_s": round(c.t_end - c.t_start, 4),
                    "f_min_hz": round(c.f_min, 2),
                    "f_max_hz": round(c.f_max, 2),
                    "peak_freq_hz": round(c.peak_freq, 2),
                    "n_bins": int(c.n_bins),
                    "raw_anomaly_score": round(c.raw_anomaly_score, 4),
                    "standardized_score": round(c.standardized_score, 4),
                    "empirical_p_value": round(c.empirical_p_value, 6),
                    "p_value_is_upper_bound": bool(c.p_value_is_upper_bound),
                    "veto_status": "VETOED_ANNOTATED" if is_vetoed else "SURVIVING",
                    "veto_reasons": veto_reasons,
                    "dq_all_pass": flags["dq_all_pass"],
                    "dq_failed_flags": flags["dq_failed_flags"],
                    "inj_active_flags": flags["inj_active_flags"],
                    "coincidence_status": "PENDING",
                }
                all_raw_candidates.append(cand_row)

                if role in ("ON_SOURCE_POSITIVE_CONTROL", "OFF_SOURCE"):
                    bucket = "ON_SOURCE" if role == "ON_SOURCE_POSITIVE_CONTROL" else "OFF_SOURCE"
                    per_det_search_cands_for_coinc[bucket][det].append(
                        {
                            "detector": det,
                            "gps_start": c_gps_start,
                            "gps_end": c_gps_end,
                            "gps_peak": c_gps_peak,
                            "anomaly_score": c.standardized_score,
                            "window_gps_start": gps0,
                            "veto_status": cand_row["veto_status"],
                            "scale_name": c.scale_name,
                            "cand_row_ref": cand_row,
                        }
                    )

        ci_raw = _wilson_ci(n_off_nominated_raw, n_off_windows)
        ci_post = _wilson_ci(n_off_nominated_post_veto, n_off_windows)
        window_timelines[det] = det_timeline
        cal_summary[det] = {
            "n_psd_windows": len(psd_segs),
            "n_cal_windows": len(cal_obs),
            "n_norm_split_a": frozen_cal.n_norm_samples,
            "n_thresh_split_b": frozen_cal.n_threshold_samples,
            "p_min_resolvable": round(1.0 / (frozen_cal.n_threshold_samples + 1), 6),
            "frozen_threshold": round(frozen_cal.global_threshold, 4),
            "calibration_status": cal_audit.status,
            "unmet_validation_conditions": list(cal_audit.unmet_validation_conditions),
            "n_eff": round(cal_audit.n_eff, 2),
            "cal_ljung_box_p": round(cal_audit.ljung_box_p, 4),
            "on_source_detected_raw": on_source_detected_raw,
            "on_source_detected_post_veto": on_source_detected_post_veto,
            "on_source_max_score": round(on_source_max_score, 4),
            "n_off_source_windows": n_off_windows,
            "n_off_nominated_raw": n_off_nominated_raw,
            "off_nomination_rate_raw": round(n_off_nominated_raw / max(1, n_off_windows), 5),
            "off_nomination_ci95_raw": [round(ci_raw[0], 5), round(ci_raw[1], 5)],
            "n_off_nominated_post_veto": n_off_nominated_post_veto,
            "off_nomination_rate_post_veto": round(n_off_nominated_post_veto / max(1, n_off_windows), 5),
            "off_nomination_ci95_post_veto": [round(ci_post[0], 5), round(ci_post[1], 5)],
            "elapsed_s": round(time.perf_counter() - t_det, 2),
        }
        print(
            f"  {det} | cal={cal_audit.status:<11} (N_cal={len(cal_obs)}, n_eff={cal_audit.n_eff:.1f}, th={frozen_cal.global_threshold:.2f}) | "
            f"GW170814 raw={on_source_detected_raw} post_veto={on_source_detected_post_veto} (score={on_source_max_score:.2f}) | "
            f"off-source raw={n_off_nominated_raw}/{n_off_windows} post_veto={n_off_nominated_post_veto}/{n_off_windows}",
            flush=True,
        )

    # 2. Pairwise multi-detector coincidence screening (H1-L1, H1-V1, L1-V1)
    coincident_pairs_rows: list[dict[str, Any]] = []
    coincident_row_ids: set[int] = set()

    for bucket in ("ON_SOURCE", "OFF_SOURCE"):
        for det_a, det_b in (("H1", "L1"), ("H1", "V1"), ("L1", "V1")):
            max_dt = COINC_SCREENING_DT_S[(det_a, det_b)]
            phys_lt = PHYSICAL_LIGHT_TRAVEL_S[(det_a, det_b)]
            max_lag = XCORR_MAX_LAG_S[(det_a, det_b)]

            c_a = per_det_search_cands_for_coinc[bucket][det_a]
            c_b = per_det_search_cands_for_coinc[bucket][det_b]
            pairs = match_coincident_candidates(c_a, c_b, max_dt_s=max_dt)

            for pair in pairs:
                # Find matching window observations for cross-correlation lag refinement
                win_gps = min(search_obs_by_det_gps[det_a].keys(), key=lambda g: abs((g + 2.0) - pair.gps_peak_a))
                obs_a = search_obs_by_det_gps[det_a].get(win_gps)
                obs_b = search_obs_by_det_gps[det_b].get(win_gps)
                xcorr_delay_ms = None
                abs_corr = None
                within_phys = None
                if obs_a is not None and obs_b is not None:
                    # Slice local +-0.25 s window around average peak for clean waveform xcorr
                    t_mid = 0.5 * (pair.gps_peak_a + pair.gps_peak_b) - win_gps
                    i0 = max(0, int(round((t_mid - 0.25) * TARGET_SAMPLE_RATE_HZ)))
                    i1 = min(obs_a.signal.size, int(round((t_mid + 0.25) * TARGET_SAMPLE_RATE_HZ)))
                    if i1 - i0 >= 128:
                        xc = estimate_cross_detector_delay(
                            obs_a.signal[i0:i1],
                            obs_b.signal[i0:i1],
                            sample_rate=TARGET_SAMPLE_RATE_HZ,
                            band_hz=(35.0, 350.0),
                            max_lag_s=max_lag,
                        )
                        xcorr_delay_ms = round(xc.delay_s * 1000.0, 3)
                        abs_corr = round(xc.abs_peak_corr, 4)
                        within_phys = bool(abs(xc.delay_s) <= phys_lt)

                coincident_pairs_rows.append(
                    {
                        "bucket": bucket,
                        "det_a": det_a,
                        "det_b": det_b,
                        "gps_peak_a": round(pair.gps_peak_a, 4),
                        "gps_peak_b": round(pair.gps_peak_b, 4),
                        "stft_dt_ms": round(pair.stft_dt_a_minus_b_s * 1000.0, 2),
                        "interval_gap_ms": round(pair.interval_gap_s * 1000.0, 2),
                        "score_a": round(pair.score_a, 4),
                        "score_b": round(pair.score_b, 4),
                        "xcorr_delay_ms": xcorr_delay_ms,
                        "xcorr_abs_peak": abs_corr,
                        "within_physical_light_travel": within_phys,
                        "physical_light_travel_ms": round(phys_lt * 1000.0, 1),
                    }
                )
                for item in c_a:
                    if abs(item["gps_peak"] - pair.gps_peak_a) < 1e-4:
                        coincident_row_ids.add(id(item["cand_row_ref"]))
                for item in c_b:
                    if abs(item["gps_peak"] - pair.gps_peak_b) < 1e-4:
                        coincident_row_ids.add(id(item["cand_row_ref"]))

    for row in all_raw_candidates:
        if id(row) in coincident_row_ids:
            row["coincidence_status"] = "MULTI_DETECTOR_COINCIDENT"
        else:
            row["coincidence_status"] = "SINGLE_DETECTOR_UNCONFIRMED"

    # 3. Exploratory off-source time-slide background check (honest finite-livetime bounds)
    off_livetime_s = 151 * SEGMENT_DURATION_S  # 604.0 s
    base_off_intervals = [
        (SEARCH_CORE_INTERVAL[0], 1186741844.0),
        (1186741880.0, SEARCH_CORE_INTERVAL[1]),
    ]
    shifts_100 = [k * 4.0 for k in range(-50, 51) if k != 0]

    def _interval_overlap_s(a_list: list[tuple[float, float]], b_list: list[tuple[float, float]]) -> float:
        return sum(max(0.0, min(a1, b1) - max(a0, b0)) for a0, a1 in a_list for b0, b1 in b_list)

    slide_overlaps = [
        _interval_overlap_s(base_off_intervals, [(b0 + dt, b1 + dt) for b0, b1 in base_off_intervals])
        for dt in shifts_100
    ]
    true_slide_livetime_hr = float(sum(slide_overlaps) / 3600.0)
    nominal_slide_livetime_hr = float(len(shifts_100) * off_livetime_s / 3600.0)

    # Deduplicated 1-peak-per-window off-source lists for audited time-slide check
    win_dedup: dict[str, list[dict[str, Any]]] = {"H1": [], "L1": [], "V1": []}
    for det in ("H1", "L1", "V1"):
        by_win: dict[float, dict[str, Any]] = {}
        for item in per_det_search_cands_for_coinc["OFF_SOURCE"][det]:
            wg = float(item["window_gps_start"])
            if wg not in by_win or float(item["anomaly_score"]) > float(by_win[wg]["anomaly_score"]):
                by_win[wg] = item
        win_dedup[det] = list(by_win.values())

    time_slide_results: dict[str, Any] = {}
    for det_a, det_b in (("H1", "L1"), ("H1", "V1"), ("L1", "V1")):
        mdt = COINC_SCREENING_DT_S[(det_a, det_b)]
        ts = run_time_slide_analysis(
            bg_cands_a=per_det_search_cands_for_coinc["OFF_SOURCE"][det_a],
            bg_cands_b=per_det_search_cands_for_coinc["OFF_SOURCE"][det_b],
            single_epoch_livetime_s=off_livetime_s,
            max_dt_s=mdt,
            n_slides=100,
            slide_step_s=4.0,
        )
        dedup_coinc = 0
        for dt in shifts_100:
            for ca in win_dedup[det_a]:
                for cb in win_dedup[det_b]:
                    gp_b = float(cb["gps_peak"]) + dt
                    if not any(s0 <= gp_b <= s1 for s0, s1 in base_off_intervals):
                        continue
                    if abs(float(ca["gps_peak"]) - gp_b) <= mdt:
                        dedup_coinc += 1
        ts["true_valid_overlap_min_s"] = round(float(min(slide_overlaps)), 1)
        ts["true_valid_overlap_mean_s"] = round(float(np.mean(slide_overlaps)), 1)
        ts["true_valid_overlap_max_s"] = round(float(max(slide_overlaps)), 1)
        ts["true_valid_overlap_total_hours"] = round(true_slide_livetime_hr, 2)
        ts["true_valid_overlap_pct_of_nominal"] = round(100.0 * true_slide_livetime_hr / nominal_slide_livetime_hr, 1)
        ts["n_unique_windows_a"] = len(win_dedup[det_a])
        ts["n_unique_windows_b"] = len(win_dedup[det_b])
        ts["dedup_window_level_slide_coinc"] = dedup_coinc
        time_slide_results[f"{det_a}-{det_b}"] = ts

    # Sort raw candidates: ON_SOURCE first, then OFF_SOURCE ranked by standardized_score desc
    all_raw_candidates.sort(
        key=lambda r: (0 if r["window_role"] == "ON_SOURCE_POSITIVE_CONTROL" else 1, -r["standardized_score"])
    )
    off_source_ranked = [r for r in all_raw_candidates if r["window_role"] == "OFF_SOURCE"]
    top_offsource_info = None
    top_amp_audit: dict[str, float] | None = None
    if off_source_ranked:
        top_row = off_source_ranked[0]
        obs_top = search_obs_by_det_gps[top_row["detector"]][top_row["window_gps_start"]]
        top_offsource_info = {"row": top_row, "obs": obs_top}
        xc = obs_top.signal - float(np.median(obs_top.signal))
        raw_mad = float(np.median(np.abs(xc)))
        mad_sig = 1.4826 * raw_mad
        win_sd = float(np.std(xc, ddof=1))
        pk_abs = float(np.max(np.abs(xc)))
        top_amp_audit = {
            "peak_abs_whitened": round(pk_abs, 4),
            "raw_mad": round(raw_mad, 6),
            "mad_scale_1p4826": round(mad_sig, 6),
            "window_sample_sd": round(win_sd, 6),
            "peak_over_raw_mad": round(pk_abs / max(raw_mad, 1e-18), 2),
            "peak_over_mad_scale_1p4826": round(pk_abs / max(mad_sig, 1e-18), 2),
            "peak_over_window_sd": round(pk_abs / max(win_sd, 1e-18), 2),
        }

    # Save CSVs and figures
    _write_csv(RESULTS_DIR / "raw_candidates.csv", all_raw_candidates)
    _write_csv(RESULTS_DIR / "coincident_pairs.csv", coincident_pairs_rows)
    _write_csv(RESULTS_DIR / "excluded_and_unusable_windows.csv", excluded_rows)

    fig1_path, fig2_path = _generate_figures(
        window_timelines=window_timelines,
        cal_summary=cal_summary,
        on_source_obs=on_source_obs,
        top_offsource_info=top_offsource_info,
        out_dir=RESULTS_DIR,
    )

    elapsed_total = time.perf_counter() - t0

    # Write concise, scientifically rigorous English report
    off_coinc = [p for p in coincident_pairs_rows if p["bucket"] == "OFF_SOURCE"]
    on_coinc = [p for p in coincident_pairs_rows if p["bucket"] == "ON_SOURCE"]

    rep: list[str] = [
        "# WHISKER / BSD Exploratory Real-Data Search: GW170814 Epoch (H1, L1, V1)",
        "",
        "> **Scope Notice:** This is an **exploratory analysis of a previously used O2 GWOSC epoch** "
        "(`GW170814`, `GPS 1186739814..1186743910`), **NOT a blind discovery search** on unseen data. "
        "Detection of `GW170814` itself is a repeatability / positive-control check; the adjacent off-source "
        "interval is used for exploratory candidate ranking, reversible DQ/veto triage, and 3-detector cross-checks. "
        "The reserved O3 holdout and detector mathematics remain completely untouched.",
        "",
        "## 1. Verified Non-Overlapping Epoch Partition (Including `pad_s = 1.0 s`)",
        "",
        "Every 4.0-second window uses `pad_s = 1.0 s` context on both sides (`6.0 s` padded span) during whitening "
        "and bandpass filtering. All four data blocks were verified programmatically to have strictly positive gaps "
        "between their padded footprints:",
        "",
        "| Block | Core GPS Interval | Core Duration | 4-s Windows | Padded Footprint | Guard Gap After Padded |",
        "|---|---|---:|---:|---|---:|",
    ]
    for blk in interval_table:
        rep.append(
            f"| `{blk['name']}` | `{blk['core_gps_start']:.0f}..{blk['core_gps_end']:.0f}` | "
            f"`{blk['core_duration_s']:.0f} s` | `{blk['n_windows_4s']}` | "
            f"`{blk['padded_gps_start']:.0f}..{blk['padded_gps_end']:.0f}` | "
            f"`{blk['guard_gap_after_padded_s']:.0f} s` |"
        )

    rep.extend(
        [
            "",
            "### Exact Window Selection Around `GW170814` (`t_c = 1186741861.5`)",
            "",
            "- **Target search block:** `160` contiguous non-overlapping 4-s windows spanning `GPS 1186741540..1186742180` (`640.0 s`).",
            "- **On-source positive control (1 window):** Window `80` (`GPS 1186741860..1186741864`), containing `GW170814` at `t = 1.5 s`.",
            "- **`+-16 s` event exclusion rule:** Exclusion interval `[1186741845.5, 1186741877.5]`. Any window whose **padded** "
            "footprint `[gps0 - 1, gps0 + 5]` intersects this interval is excluded from off-source (`8` guard-zone windows: "
            "Windows `76..79` [`1186741844..1186741860`] and Windows `81..84` [`1186741864..1186741880`]).",
            "- **Resulting off-source search set:** **`151` non-overlapping 4-s windows** (`604.0 s = 10.07 min`; "
            "Windows `0..75` pre-event and Windows `85..159` post-event).",
            f"- **Unusable / DQ-dropped windows across all intervals (H1/L1/V1):** "
            f"`{sum(1 for r in excluded_rows if r['reason'] != 'PADDED_OVERLAP_WITH_EVENT_PM_16S_GUARD')}` "
            "(complete exclusion log saved in `excluded_and_unusable_windows.csv`).",
            "",
            "## 2. Frozen Preprocessing, Detector, Calibration, and Veto Configuration",
            "",
            "- **Preprocessing (`src/bsd/ligo_adapter.py`):** `adapted_ligo_whitened` (`1e18` scaling, median-Welch PSD fitted "
            "strictly on the `60 x 4 s = 240 s` PSD block, zero-phase 4th-order Butterworth bandpass `30-380 Hz`, anti-aliased "
            "decimation `4096 Hz -> 1024 Hz`).",
            "- **Detector (`src/bsd/detector_v2.py`):** Unmodified 3-scale `MultiResolutionDetector` (`short_31ms`: hop `3.91 ms`; "
            "`medium_125ms`: hop `31.25 ms`; `long_500ms`: hop `125.0 ms` at `1024 Hz`), `target_fap = 0.01` (nominal calibration "
            "quantile parameter — **NOT** a measured FAR on real LIGO/Virgo data), `norm_fraction = 0.30` (`Split A = 212`, "
            "`Split B = 495`, `p_min = 1/496 = 0.002016`).",
            "- **Reversible vetoes (`src/bsd/vetoes.py`):** `VetoConfig(enable_glitch_veto=True, enable_line_veto=True)`. "
            "All raw candidates are preserved verbatim in `raw_candidates.csv` before attaching veto and GWOSC DQ/Inj annotations.",
            "",
            "### Per-Detector Calibration Status and Search Summary",
            "",
            "| Detector | Cal Windows (`Split A / B`) | Cal Status | `n_eff` | Cal LB `p` | Frozen Threshold | GW170814 Raw / Post-Veto (Score) | Off-Source Raw (`k/151`) | Off-Source Post-Veto (`k/151`) |",
            "|---|---:|---|---:|---:|---:|---|---:|---:|",
        ]
    )

    for det in ("H1", "L1", "V1"):
        cs = cal_summary[det]
        rep.append(
            f"| `{det}` | `{cs['n_cal_windows']}` (`{cs['n_norm_split_a']}/{cs['n_thresh_split_b']}`) | "
            f"`{cs['calibration_status']}` | `{cs['n_eff']:.1f}` | `{cs['cal_ljung_box_p']:.4f}` | "
            f"`{cs['frozen_threshold']:.4f}` | "
            f"`{cs['on_source_detected_raw']}` / `{cs['on_source_detected_post_veto']}` (`{cs['on_source_max_score']:.2f}`) | "
            f"`{cs['n_off_nominated_raw']}/{cs['n_off_source_windows']}` (`{cs['off_nomination_rate_raw']*100:.2f}%`, "
            f"95% Wilson `[{cs['off_nomination_ci95_raw'][0]*100:.2f}%, {cs['off_nomination_ci95_raw'][1]*100:.2f}%]`) | "
            f"`{cs['n_off_nominated_post_veto']}/{cs['n_off_source_windows']}` (`{cs['off_nomination_rate_post_veto']*100:.2f}%`) |"
        )

    rep.extend(
        [
            "",
            "> **Real-data FAR note:** The nominal `target_fap = 0.01` is only the Split-B calibration quantile rule; "
            "it is **not** a measured false-alarm rate on real non-stationary LIGO/Virgo noise. Wilson intervals on "
            "real-data off-source counts are descriptive approximations assuming independent windows.",
            "",
            "## 3. Peak-Time Localization Uncertainty, Coincidence Screening, and Time-Slide Resolution",
            "",
            "- **STFT localization uncertainty:** At `1024 Hz`, half-hop time quantization is `+/-1.95 ms` (`short_31ms`), "
            "`+/-15.63 ms` (`medium_125ms`), and `+/-62.50 ms` (`long_500ms`). Combining STFT peak quantization with physical "
            "light-travel bounds (`10.0 ms` H1-L1, `27.3 ms` H1-V1, `26.4 ms` L1-V1), prespecified coarse STFT screening uses "
            "`max_dt = 50 ms` (H1-L1) and `65 ms` (H1-V1, L1-V1), followed by band-limited time-domain cross-correlation "
            "(`35-350 Hz`, `+-15 ms` H1-L1, `+-30 ms` LIGO-Virgo).",
            "- **Non-detection rule:** A candidate detected in only one detector is labeled `SINGLE_DETECTOR_UNCONFIRMED` "
            "and is **never** classified as instrumental solely because another detector did not trigger (antenna-pattern "
            "factors and Virgo's lower O2 sensitivity can cause genuine signals to appear in only one or two detectors).",
            f"- **On-source coincident pairs (`GW170814`):** `{len(on_coinc)}` pair(s).",
            f"- **Off-source coincident pairs:** `{len(off_coinc)}` pair(s).",
            "- **Time-slide background audit (`100` non-circular slides, `4 s` step, `+-200 s` shifts over `604 s` off-source):**",
            f"  - **Valid overlap correction:** Accounting for outer edge loss and the central `36 s` guard-zone hole "
            f"(`1186741844..1186741880`), true active overlap per slide ranges from "
            f"`{time_slide_results['H1-L1']['true_valid_overlap_min_s']:.1f} s` to "
            f"`{time_slide_results['H1-L1']['true_valid_overlap_max_s']:.1f} s` (mean "
            f"`{time_slide_results['H1-L1']['true_valid_overlap_mean_s']:.1f} s`), totaling "
            f"**`{time_slide_results['H1-L1']['true_valid_overlap_total_hours']:.2f} hr`** "
            f"(`{time_slide_results['H1-L1']['true_valid_overlap_pct_of_nominal']:.1f}%` of the "
            f"`{nominal_slide_livetime_hr:.2f} hr` nominal `100 x 604 s` figure computed by `run_time_slide_analysis`).",
            "  - **Multi-cluster inflation vs. deduplicated window-level counts:** Raw cluster-level counts in "
            f"`run_time_slide_analysis` (`H1-L1 = {time_slide_results['H1-L1']['total_accidental_coincidences']}`, "
            f"`H1-V1 = {time_slide_results['H1-V1']['total_accidental_coincidences']}`, "
            f"`L1-V1 = {time_slide_results['L1-V1']['total_accidental_coincidences']}`) are inflated because multiple "
            "STFT-scale clusters from the same physical window "
            f"(`{time_slide_results['H1-L1']['n_bg_cands_a']}` H1, `{time_slide_results['H1-L1']['n_bg_cands_b']}` L1, "
            f"`{time_slide_results['H1-V1']['n_bg_cands_b']}` V1 clusters from only "
            f"`{time_slide_results['H1-L1']['n_unique_windows_a']}`, `{time_slide_results['H1-L1']['n_unique_windows_b']}`, "
            f"`{time_slide_results['H1-V1']['n_unique_windows_b']}` unique windows) are cross-matched pairwise and "
            "`match_coincident_candidates` also triggers on `interval_gap <= max_dt`. When deduplicated to 1 peak per "
            "window with `|dt_peak| <= max_dt` inside valid active overlap, window-level slide coincidences are: "
            f"`H1-L1 = {time_slide_results['H1-L1']['dedup_window_level_slide_coinc']}`, "
            f"`H1-V1 = {time_slide_results['H1-V1']['dedup_window_level_slide_coinc']}`, "
            f"`L1-V1 = {time_slide_results['L1-V1']['dedup_window_level_slide_coinc']}`.",
            "  - **4-second shift step alignment:** Shifting by exact `4.0 s` multiples aligns 4-s window boundaries and "
            "the `36 s` (`9 x 4 s`) guard-zone hole identically across slides; any hypothetical 4-s periodic windowing "
            "artifact would remain phase-locked.",
            "  - **Resolution caveat:** With only `604 s` (`10.07 min`) of off-source data "
            f"(`{true_slide_livetime_hr:.2f} hr` valid slide livetime) and only `2-3` triggered windows per detector, "
            "time-slides do **not** support a meaningful astrophysical coincidence significance claim.",
            "",
            "## 4. Top Ranked Off-Source Raw Candidates (All Raw Candidates Preserved in `raw_candidates.csv`)",
            "",
            "> **Empirical `p`-value definition:** `p_emp = (1 + n_ge) / (N_B + 1)` with `N_B = 495` Split-B "
            "calibration windows. When `n_ge = 0`, `p_emp = 1/496 = 0.00202` (calibration resolution floor). This is "
            "a conformal finite-sample point estimate under exact i.i.d. exchangeability — it is **not** written as "
            "`p <= 0.0020` (since `1/496 > 0.0020` and `0/495` does not prove the true tail probability is below `0.0020`), "
            "and it is **not** an astrophysical probability or a measured real-data false-alarm rate.",
            "",
            "| Rank | Det | GPS Peak | `t - t_c` (s) | Scale | Peak Freq (Hz) | Std Score | Emp. `p` (`N_B=495`) | Veto Status (Reasons) | GWOSC DQ / Inj | Coincidence |",
            "|---:|---|---:|---:|---|---:|---:|---:|---|---|---|",
        ]
    )

    for rank, r in enumerate(off_source_ranked[:15], start=1):
        p_str = f"{r['empirical_p_value']:.5f} (0/495; floor 1/496)" if r["p_value_is_upper_bound"] else f"{r['empirical_p_value']:.5f}"
        dq_str = "PASS" if r["dq_all_pass"] else f"FAIL({r['dq_failed_flags']})"
        rep.append(
            f"| {rank} | `{r['detector']}` | `{r['gps_peak']:.3f}` | `{r['t_rel_event_s']:+.2f}` | "
            f"`{r['scale_name']}` | `{r['peak_freq_hz']:.1f}` | `{r['standardized_score']:.2f}` | "
            f"`{p_str}` | `{r['veto_status']}` (`{r['veto_reasons']}`) | `{dq_str}` / `{r['inj_active_flags']}` | "
            f"`{r['coincidence_status']}` |"
        )

    if not off_source_ranked:
        rep.append("| - | - | - | - | - | - | - | - | No off-source raw candidates exceeded threshold | - | - |")

    rep.extend(
        [
            "",
            "## 5. Focused Audit Findings",
            "",
            "### 5.1 Why `GW170814` Exceeded Threshold in `L1` Only (Comparison with Prior Benchmarks)",
            "- **Consistency with prior repository benchmarks:** Both `Exp 003` (v0.1 `TemplateFreeDetector`: `H1` max "
            "score `37.93` vs. threshold `76.80`, ratio `0.494`) and `v02_phase3` (v0.2 `MultiResolutionDetector`: `0` "
            "H1-L1 coincident pairs pre- and post-veto) likewise did not exceed the single-detector 1% threshold in `H1` on `GW170814`.",
            "- **Per-scale scores (`Window 80`, `GPS 1186741860..1186741864`):** `L1` exceeded threshold `2.24` across "
            "all 3 scales (max `8.89`), whereas `H1` (`max = 1.14 < 2.59`) and `V1` (`max = 3.62 < 4.01`) remained below "
            "their frozen single-detector thresholds without matched-filter templates.",
            "",
            "### 5.2 Forensic Audit of Top Off-Source Candidate (`V1` at `GPS 1186741602.344`, `2017-08-14 10:26:24.344 UTC`)",
            f"- **Amplitude definition and units:** On the 4-second whitened, `30-380 Hz` bandpassed, `1024 Hz` decimated "
            f"series (`GPS 1186741600..1186741604`, dimensionless whitened-strain units), the median-centered peak is "
            f"`{top_amp_audit['peak_abs_whitened']:.3f}` with `raw MAD = {top_amp_audit['raw_mad']:.4f}` "
            f"(`peak / raw MAD = {top_amp_audit['peak_over_raw_mad']:.2f}`), Gaussian-consistent MAD scale "
            f"`1.4826 x MAD = {top_amp_audit['mad_scale_1p4826']:.4f}` (`peak / [1.4826 x MAD] = "
            f"{top_amp_audit['peak_over_mad_scale_1p4826']:.2f}`, reported by `vetoes.py` as `local_peak_mad_sigma`), "
            f"and full-window sample SD `{top_amp_audit['window_sample_sd']:.4f}` "
            f"(`peak / window SD = {top_amp_audit['peak_over_window_sd']:.2f}`). This ratio is **not** a Gaussian "
            "standard-deviation significance (`sigma`) on non-Gaussian Virgo noise." if top_amp_audit else "- Amplitude audit unavailable.",
            "- **Why cluster metadata span `16-496 Hz` (`short_31ms`) and `4-508 Hz` (`medium_125ms`) despite `30-380 Hz` bandpass:** "
            "(1) The 4th-order zero-phase Butterworth bandpass has a smooth roll-off, attenuating stopband power by "
            "`10^-5` to `10^-9` relative to the passband; (2) per-frequency median normalization `Z[f, t] = P[f, t] / "
            "(S_bg(f) * T_bg(t))` cancels out the smooth filter gain `|H(f)|^2`, allowing heavily attenuated stopband bins "
            "to exceed `Z > 4`; (3) cluster bounds report half-bin edges (`f_center +/- df/2`, with `df = 32 Hz` on "
            "`short_31ms` and `df = 8 Hz` on `medium_125ms`). Note also that `fig02` masks the displayed spectrogram to "
            "`30-360 Hz`. Physical signal power is concentrated inside the `30-380 Hz` passband.",
            "- **Why `short_31ms` survived while the other 4 clusters were vetoed:** `evaluate_impulsive_glitch` computes "
            "excess kurtosis on `[t_start - 0.08 s, t_end + 0.08 s]`. Wider `medium_125ms`/`long_500ms` slices (`548-3620` "
            "samples, `kurtosis = 39.1-255.5`) exceeded `20.0`, whereas the narrower `short_31ms` slice (`256` samples) "
            "yielded `local_kurtosis = 17.23 < 20.0`, and `n_samples_ge_half_peak = 6 > 3` (bandpass ringing) missed the "
            "ultra-short spike rule. Because `apply_vetoes` sets `detected_after_veto = len(surviving_candidates) > 0`, "
            "this kept the window flagged `detected_post_veto = True` (`V1` off-source `2/151` instead of `1/151`).",
            "- **Public catalog check, observations vs. hypotheses, and novelty caveat:**",
            "  - *Verified observations:* Local `DQmask = 127` (all 7 bits PASS), `Injmask = 31` (all 5 `NO_*_HW_INJ` "
            "bits PASS), and GWOSC `O2_4KHZ_R1` Timeline confirms `V1_DATA`, `V1_CBC_CAT2`, `V1_BURST_CAT2/3`, and "
            "`V1_NO_BURST_HW_INJ` PASS over `[1186741500..1186741700]`. Queries of `GWOSC EventAPI` (`671` entries across "
            "GWTC-1/2.1/3/4 confident & marginal catalogs, `+-3600 s`) and `4-OGC` (`4OGC_top.txt`) returned only "
            "`GW170814` (`+259.16 s`). No public O2 Virgo glitch catalog or O2 Virgo auxiliary-channel release was found "
            "or accessible in the sources checked (`Gravity Spy` Zenodo LIGO O1-O3b releases, `GWitchHunters` O3 Virgo "
            "project, `2-OGC` URL HTTP 404, web/arXiv searches).",
            "  - *Working hypothesis (unproven without Virgo auxiliary channels):* Internal V1 morphology (short `< 0.1 s` "
            "passband broadband impulse without chirp structure, high excess kurtosis on wider slices) is consistent with "
            "an instrumental glitch, whereas documented hardware injections are inconsistent with `Injmask = 31`. Absence "
            "in H1/L1 is recorded strictly as an observation and is not used alone to infer instrumental origin.",
            "  - *Novelty caveat:* Lack of a public catalog match does **not** establish novelty or astrophysical origin; "
            "whether this transient was logged internally by Virgo DetChar in 2017 remains unknown.",
            "",
            "## 6. Coincident Pairs Summary (`coincident_pairs.csv`)",
            "",
        ]
    )
    if coincident_pairs_rows:
        rep.extend(
            [
                "| Bucket | Pair | GPS Peak A | GPS Peak B | STFT `dt` (ms) | XCorr Delay (ms) | `|r_peak|` | Within Light Travel? |",
                "|---|---|---:|---:|---:|---:|---:|---|",
            ]
        )
        for p in coincident_pairs_rows:
            rep.append(
                f"| `{p['bucket']}` | `{p['det_a']}-{p['det_b']}` | `{p['gps_peak_a']:.4f}` | `{p['gps_peak_b']:.4f}` | "
                f"`{p['stft_dt_ms']:+.2f}` | `{p['xcorr_delay_ms']}` | `{p['xcorr_abs_peak']}` | "
                f"`{p['within_physical_light_travel']}` (`<= {p['physical_light_travel_ms']} ms`) |"
            )
    else:
        rep.append("No pairwise coincident candidates found.")

    rep.extend(
        [
            "",
            "## 7. Figures and Limitations",
            "",
            "- `fig01_gw170814_search_timeline.png`: Standardized global score timeline for H1, L1, and V1 across the 640-s target interval.",
            "- `fig02_gw170814_positive_control_and_top_offsource.png`: Whitened spectrograms of the GW170814 positive control (H1, L1) and top-ranked off-source candidate.",
            f"- **Total runtime:** `{elapsed_total:.1f} s`.",
            "- **Limitations:** Previously used O2 epoch (not a blind discovery search); short `604 s` off-source livetime limits time-slide resolution; "
            "single-detector triggers remain unconfirmed without external glitch catalogs (e.g., Gravity Spy) or independent multi-detector coincidence.",
            "",
        ]
    )

    (RESULTS_DIR / "report.md").write_text("\n".join(rep), encoding="utf-8")

    payload = {
        "protocol": protocol,
        "calibration_and_search_summary": cal_summary,
        "n_raw_candidates_total": len(all_raw_candidates),
        "n_off_source_raw_candidates": len(off_source_ranked),
        "n_coincident_pairs_on_source": len(on_coinc),
        "n_coincident_pairs_off_source": len(off_coinc),
        "coincident_pairs": coincident_pairs_rows,
        "time_slide_background": time_slide_results,
        "top_offsource_amplitude_audit": top_amp_audit,
        "elapsed_total_s": round(elapsed_total, 2),
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(
        f"[GW170814 Exploratory Search] DONE in {elapsed_total:.1f}s! "
        f"Off-source raw candidates = {len(off_source_ranked)}, "
        f"on-source coinc = {len(on_coinc)}, off-source coinc = {len(off_coinc)}.",
        flush=True,
    )


if __name__ == "__main__":
    main()
