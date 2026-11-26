"""WHISKER / BSD Expanded Synthetic Negative-Control Experiment (N_eval = 2,000).

Evaluates the THREE existing frozen calibrations (120, 200, and 220 nested
background windows from commit `6a43f08`) on ONE shared, brand-new independent
evaluation set of 2,000 non-overlapping 4-second `pink_colored` windows
(8,000 s = 2.222 hours), drawn from a separate RNG stream:
  `make_split_rng(2026, "neg_control_eval_2000_v1")`.

Explicitly distinguishes:
  1. Frozen Raw Detector Outcome (`frozen_raw_stft`):
     Stateless per-window `MultiResolutionDetector.analyze(obs, frozen_cal)`
     using ONLY each window's 4096 samples and the frozen calibration.
     Trials are strictly i.i.d. Bernoulli given the frozen threshold.
  2. Workstation Adaptive Pipeline Outcome (`workstation_adaptive`):
     Full `AdaptiveRegimeManagerV2b` sequential pass (rolling 12-window
     continuum/glitch monitors + potential recalibration), reporting both
     full-timeline and clean-calibrated surviving outcomes separately.

Note: This experiment estimates performance CONDITIONAL on these three
specific frozen calibrations (`seed = 2026`). It does NOT generalize to the
full calibration procedure across random seeds or to real detector data.
Zero detector math, baselines, GUI, or O3 holdout files are touched.
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
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
EXP_DIR = PROJECT_ROOT / "experiments"
for p in (str(SRC_DIR), str(EXP_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from bsd.adaptive_regime_v2b import estimate_calibration_uncertainty
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, make_split_rng
from bsd.io_loader import LoadedTimeSeries, PreprocessingSummary, slice_into_observations
from bsd.noise import generate_noise
from bsd.workstation_runner import WorkstationConfig, run_workstation_analysis
import run_workstation_calibration_study as cal_study

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 180,
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.2,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
)

EVAL_SEED = 2026
EVAL_STREAM_NAME = "neg_control_eval_2000_v1"
N_EVAL_WINDOWS = 2000
EXPECTED_FROZEN_THRESHOLDS = {120: 2.2568, 200: 2.4250, 220: 2.3089}


def _wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Exact Wilson score confidence interval (default z=1.96 -> 95% CI)."""
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


def _generate_independent_eval_set(n_windows: int = N_EVAL_WINDOWS, fs: float = 1024.0, win_dur: float = 4.0) -> tuple[list[Observation], LoadedTimeSeries]:
    """Generate n_windows strictly independent 4-s pink_colored windows from a new RNG stream."""
    win_len = int(round(win_dur * fs))
    rng_eval = make_split_rng(EVAL_SEED, EVAL_STREAM_NAME)
    chunks: list[np.ndarray] = []
    obs_list: list[Observation] = []
    for i in range(n_windows):
        x, _ = generate_noise("pink_colored", win_len, fs, rng_eval, sigma=1.0)
        chunks.append(x)
        obs_list.append(Observation(f"neg2000_{i:04d}", x, fs, win_dur))

    full_arr = np.concatenate(chunks)
    prep = PreprocessingSummary(
        file_format="synthetic",
        channel_selected="synthetic_1d",
        units="normalized_sigma",
        original_sample_rate_hz=fs,
        effective_sample_rate_hz=fs,
        resampled=False,
        antialias_filter="none",
        effective_bandwidth_hz=(0.0, fs / 2.0),
        analysis_bands_hz=(20.0, 450.0),
        ligo_whitening_applied=False,
        steps=(f"Generated {n_windows} independent pink_colored 4-s windows (seed={EVAL_SEED}, stream={EVAL_STREAM_NAME})",),
    )
    eval_ts = LoadedTimeSeries(
        source_name=f"eval_neg_control_{n_windows}win_pink",
        file_format="synthetic",
        sha256=f"eval_{EVAL_STREAM_NAME}_{EVAL_SEED}_{full_arr.size}",
        n_samples=int(full_arr.size),
        sample_rate_hz=fs,
        duration_s=float(full_arr.size / fs),
        signal=full_arr,
        preprocessing=prep,
        time_interval=(5000.0, 5000.0 + n_windows * win_dur),
        time_reference="SYNTHETIC_CLOCK",
    )
    return obs_list, eval_ts


def _generate_figure(cond_rows: list[dict[str, Any]], fig_path: Path) -> None:
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.4))

    lbls = [f"{r['n_bg_windows']}w\n(th={r['frozen_threshold']:.3f})" for r in cond_rows]
    xx = np.arange(len(cond_rows))

    # Left: Frozen Raw Detector FAR + 95% Wilson CI vs. 1% Target FAP
    rates = np.array([r["frozen_raw_far"] * 100.0 for r in cond_rows])
    lo = np.array([r["frozen_raw_ci95_lo"] * 100.0 for r in cond_rows])
    hi = np.array([r["frozen_raw_ci95_hi"] * 100.0 for r in cond_rows])
    yerr = np.vstack([rates - lo, hi - rates])

    cols = ["#2b6cb0", "#2f855a", "#6b46c1"]
    ax1.bar(xx, rates, width=0.5, color=cols, alpha=0.80, edgecolor="white")
    ax1.errorbar(xx, rates, yerr=yerr, fmt="none", ecolor="#1a202c", elinewidth=1.6, capsize=5)
    ax1.axhline(1.0, color="#d62728", ls="--", lw=1.3, label="Target FAP (1.0%)")
    for i, r in enumerate(cond_rows):
        ax1.text(xx[i], hi[i] + 0.12, f"{r['frozen_raw_fp']}/2000\n({rates[i]:.2f}%)", ha="center", fontsize=8.2)
    ax1.set_xticks(xx)
    ax1.set_xticklabels(lbls)
    ax1.set_ylim(0, max(hi.max() + 0.8, 2.5))
    ax1.set_ylabel("Observed False-Alarm Rate (%) with 95% Wilson CI")
    ax1.set_title("Frozen Raw Detector FAR on 2,000 Independent Windows")
    ax1.legend(loc="upper right")

    # Right: Frozen Raw vs. Adaptive Full-Timeline vs. Adaptive Clean-Calibrated FAR
    w = 0.25
    r_raw = [r["frozen_raw_far"] * 100.0 for r in cond_rows]
    r_afull = [r["adaptive_full_far"] * 100.0 for r in cond_rows]
    r_acln = [r["adaptive_clean_far"] * 100.0 for r in cond_rows]
    ax2.bar(xx - w, r_raw, width=w, color="#1f4e79", label="Frozen Raw Detector (i.i.d.)")
    ax2.bar(xx, r_afull, width=w, color="#dd6b20", label="Adaptive Full-Timeline")
    ax2.bar(xx + w, r_acln, width=w, color="#2f855a", label="Adaptive Clean-Calibrated")
    ax2.axhline(1.0, color="#d62728", ls="--", lw=1.2, label="Target FAP (1.0%)")
    ax2.set_xticks(xx)
    ax2.set_xticklabels(lbls)
    ax2.set_ylim(0, max(max(r_raw), max(r_afull), max(r_acln), 2.0) * 1.35)
    ax2.set_ylabel("Observed False-Alarm Rate (%)")
    ax2.set_title("Raw Frozen vs. Adaptive Workstation Pipeline FAR")
    ax2.legend(loc="upper right", fontsize=7.6)

    fig.suptitle("WHISKER Expanded Negative Control (N_eval = 2,000 Independent pink_colored Windows)", y=0.995)
    fig.tight_layout(rect=[0, 0.01, 1, 0.95])
    fig.savefig(fig_path)
    plt.close(fig)


def main() -> None:
    t0 = time.time()
    out_dir = PROJECT_ROOT / "results" / "workstation_negative_control_2000"
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = WorkstationConfig(
        search_mode="transient",
        processing_mode="causal_trailing",
        window_duration_s=4.0,
        target_fap=0.01,
        enable_glitch_veto=False,
        enable_line_veto=False,
        calibration_seed=0,
    )

    print("[Neg2000] Rebuilding exact frozen 120 / 200 / 220 nested backgrounds (seed=2026)...")
    ds = cal_study._build_nested_datasets(seed=2026, target_snr=13.0)
    bg_map = {120: ds["bg_120_ts"], 200: ds["bg_200_ts"], 220: ds["bg_220_ts"]}

    det = MultiResolutionDetector()
    frozen_cals: dict[int, Any] = {}
    for n_bg, bg_ts in bg_map.items():
        bg_obs = slice_into_observations(bg_ts, window_duration_s=4.0, prefix="bg", min_windows=60)
        cal = det.calibrate(bg_obs, target_fap=cfg.target_fap, norm_fraction=0.30, calibration_seed=cfg.calibration_seed)
        assert abs(cal.global_threshold - EXPECTED_FROZEN_THRESHOLDS[n_bg]) < 1e-4, (
            f"Threshold mismatch for {n_bg}w: got {cal.global_threshold:.4f}, expected {EXPECTED_FROZEN_THRESHOLDS[n_bg]:.4f}"
        )
        frozen_cals[n_bg] = cal
        print(f"  Verified frozen {n_bg}w threshold = {cal.global_threshold:.4f} (matches 6a43f08)")

    print(f"[Neg2000] Generating {N_EVAL_WINDOWS} independent pink_colored evaluation windows (stream={EVAL_STREAM_NAME})...")
    eval_obs, eval_ts = _generate_independent_eval_set(N_EVAL_WINDOWS)

    # Precompute raw per-scale max scores ONCE (strictly stateless per window, mathematically identical to det.analyze)
    print("[Neg2000] Computing stateless per-window 3-scale STFT scores across 2,000 windows...")
    raw_per_win: list[dict[str, float]] = []
    for obs in eval_obs:
        per_scale = det._raw_analyze_all_scales(obs)
        raw_per_win.append({s_name: float(val[0]) for s_name, val in per_scale.items()})

    # Evaluate Frozen Raw Detector across all 3 frozen thresholds on identical 2,000 windows
    frozen_flags: dict[int, np.ndarray] = {}
    frozen_scores: dict[int, np.ndarray] = {}
    eval_dep_audit: dict[int, dict[str, Any]] = {}

    for n_bg, cal in frozen_cals.items():
        scores = np.empty(N_EVAL_WINDOWS, dtype=np.float64)
        for i, rmap in enumerate(raw_per_win):
            scores[i] = max(cal.scale_stats[s.name].standardize(rmap[s.name]) for s in det.scales)
        flags = scores > cal.global_threshold
        frozen_scores[n_bg] = scores
        frozen_flags[n_bg] = flags
        unc_eval = estimate_calibration_uncertainty(scores, target_fap=0.01)
        eval_dep_audit[n_bg] = {
            "eval_lag1_autocorr": round(unc_eval.lag1_autocorr, 4),
            "eval_ljung_box_p": round(unc_eval.ljung_box_p, 4),
            "eval_n_eff": round(unc_eval.n_eff, 1),
        }

    # Evaluate full Workstation Adaptive Pipeline (run_workstation_analysis) on the exact same 2,000 windows
    cond_rows: list[dict[str, Any]] = []
    for n_bg in (120, 200, 220):
        print(f"[Neg2000] Running full Workstation Adaptive Pipeline for {n_bg}w on 2,000 windows...")
        res = run_workstation_analysis(eval_ts, bg_map[n_bg], cfg)
        n_win = len(res.windows)
        assert n_win == N_EVAL_WINDOWS

        # 1. Frozen raw detector counts
        k_raw = int(np.sum(frozen_flags[n_bg]))
        ci_raw = _wilson_ci(k_raw, N_EVAL_WINDOWS)

        # 2. Adaptive full-timeline counts
        k_afull = sum(1 for w in res.windows if w.stft_nominated_post_veto)
        ci_afull = _wilson_ci(k_afull, n_win)

        # 3. Adaptive clean-calibrated counts (eligible windows only)
        clean_wins = [w for w in res.windows if w.is_clean_calibrated]
        n_clean = len(clean_wins)
        k_acln = sum(1 for w in clean_wins if w.stft_nominated_post_veto)
        ci_acln = _wilson_ci(k_acln, n_clean)

        cond_rows.append(
            {
                "n_bg_windows": n_bg,
                "n_thresh_split_b": frozen_cals[n_bg].n_threshold_samples,
                "calibration_status": res.calibration_audit.status,
                "frozen_threshold": round(frozen_cals[n_bg].global_threshold, 4),
                "n_eval_total": N_EVAL_WINDOWS,
                "frozen_raw_fp": k_raw,
                "frozen_raw_far": round(k_raw / N_EVAL_WINDOWS, 6),
                "frozen_raw_ci95_lo": round(ci_raw[0], 6),
                "frozen_raw_ci95_hi": round(ci_raw[1], 6),
                "adaptive_full_fp": k_afull,
                "adaptive_full_far": round(k_afull / n_win, 6),
                "adaptive_full_ci95_lo": round(ci_afull[0], 6),
                "adaptive_full_ci95_hi": round(ci_afull[1], 6),
                "adaptive_n_clean_eligible": n_clean,
                "adaptive_clean_duty_cycle": round(n_clean / n_win, 4),
                "adaptive_clean_fp": k_acln,
                "adaptive_clean_far": round(k_acln / max(1, n_clean), 6),
                "adaptive_clean_ci95_lo": round(ci_acln[0], 6),
                "adaptive_clean_ci95_hi": round(ci_acln[1], 6),
                "adaptive_n_regimes": res.summary_metrics["n_regimes_established"],
                "adaptive_cont_warnings": res.summary_metrics["n_continuum_drift_warnings"],
                "adaptive_glitch_warnings": res.summary_metrics["n_glitch_storm_warnings"],
                **eval_dep_audit[n_bg],
            }
        )

    # Paired overlap & discordance analysis on the exact same 2,000 windows (frozen raw detector)
    f120, f200, f220 = frozen_flags[120], frozen_flags[200], frozen_flags[220]
    paired = {
        "n_windows": N_EVAL_WINDOWS,
        "all_three_true": int(np.sum(f120 & f200 & f220)),
        "only_120": int(np.sum(f120 & (~f200) & (~f220))),
        "only_200": int(np.sum((~f120) & f200 & (~f220))),
        "only_220": int(np.sum((~f120) & (~f200) & f220)),
        "120_and_200_only": int(np.sum(f120 & f200 & (~f220))),
        "120_and_220_only": int(np.sum(f120 & (~f200) & f220)),
        "200_and_220_only": int(np.sum((~f120) & f200 & f220)),
        "any_true": int(np.sum(f120 | f200 | f220)),
        "none_true": int(np.sum((~f120) & (~f200) & (~f220))),
        "120_vs_200_discordant": (int(np.sum(f120 & (~f200))), int(np.sum((~f120) & f200))),
        "120_vs_220_discordant": (int(np.sum(f120 & (~f220))), int(np.sum((~f120) & f220))),
        "200_vs_220_discordant": (int(np.sum(f200 & (~f220))), int(np.sum((~f200) & f220))),
    }

    paired_rows: list[dict[str, Any]] = []
    for i in range(N_EVAL_WINDOWS):
        if f120[i] or f200[i] or f220[i]:
            paired_rows.append(
                {
                    "window_index": i,
                    "score_120": round(float(frozen_scores[120][i]), 4),
                    "fa_120": bool(f120[i]),
                    "score_200": round(float(frozen_scores[200][i]), 4),
                    "fa_200": bool(f200[i]),
                    "score_220": round(float(frozen_scores[220][i]), 4),
                    "fa_220": bool(f220[i]),
                }
            )

    _write_csv(out_dir / "summary.csv", cond_rows)
    _write_csv(out_dir / "paired_false_alarms.csv", paired_rows)
    _generate_figure(cond_rows, out_dir / "figures" / "fig01_expanded_neg_control_2000.png")

    # Write English scientific report
    lines: list[str] = [
        "# WHISKER / BSD Expanded Synthetic Negative-Control Experiment (N_eval = 2,000)",
        "",
        "**Scope:** Conditional evaluation of the three frozen nested calibrations (`120`, `200`, `220` windows from commit `6a43f08`, `seed = 2026`) on **2,000 brand-new independent `pink_colored` evaluation windows** (`8,000 s` = `2.222 hours`).",
        "**Reproducibility:** Evaluation RNG stream = `make_split_rng(2026, 'neg_control_eval_2000_v1')`, strictly disjoint from calibration and prior evaluation streams. Zero thresholds tuned on evaluation data. Zero O3 holdout access.",
        "",
        "## 1. Generator & Pipeline Independence Audit",
        "",
        "- **Noise generator (`src/bsd/noise.py:generate_pink_colored`):** Each 4-second window (`4096` samples at `1024 Hz`) is generated by an independent call drawing a fresh spectral slope `alpha ~ Uniform(0.85, 1.15)` and fresh Gaussian white FFT bins, shaped via per-window `rfft/irfft` and RMS-normalized independently to `sigma = 1.0`. **No filter state, normalization, or memory carries across windows.**",
        "- **Frozen Raw Detector (`MultiResolutionDetector.analyze`):** Strictly stateless per window; therefore, binary false-alarm indicators `I(score_i > threshold)` across the `2,000` evaluation windows satisfy the i.i.d. Bernoulli assumption of the 95% Wilson score interval (conditional on the frozen threshold).",
        "- **Adaptive Workstation Pipeline (`AdaptiveRegimeManagerV2b`):** Uses rolling `12`-window continuum and glitch-storm monitors plus finite-state recalibration. Its window-level decisions can exhibit serial dependence across adjacent windows; accordingly, raw frozen outcomes and adaptive outcomes are reported strictly separately.",
        "- **False-Alarm Definition:** For a single non-overlapping 4-second pure-noise window, a false alarm occurs when `stft_max_score > active_threshold` (strict inequality, `enable_glitch_veto=False`, `enable_line_veto=False`).",
        "",
        "## 2. Primary Results: Frozen Raw Detector vs. Adaptive Workstation Pipeline (N = 2,000)",
        "",
        "| Condition | Status | Frozen Threshold | Frozen Raw FA (`k/N`) | Frozen Raw FAR (95% Wilson CI) | Adaptive Full FA (`k/N`) | Adaptive Clean FA (`k/N_clean`) | Clean Duty | Adaptive Regimes | Eval `Ljung-Box p` |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in cond_rows:
        lines.append(
            f"| `{r['n_bg_windows']}w` (`N_th={r['n_thresh_split_b']}`) | `{r['calibration_status']}` | `{r['frozen_threshold']:.4f}` | "
            f"**`{r['frozen_raw_fp']}/{r['n_eval_total']}`** | **`{r['frozen_raw_far']*100:.2f}%`** `[{r['frozen_raw_ci95_lo']*100:.2f}%, {r['frozen_raw_ci95_hi']*100:.2f}%]` | "
            f"`{r['adaptive_full_fp']}/{r['n_eval_total']}` (`{r['adaptive_full_far']*100:.2f}%`) | "
            f"`{r['adaptive_clean_fp']}/{r['adaptive_n_clean_eligible']}` (`{r['adaptive_clean_far']*100:.2f}%`) | "
            f"`{r['adaptive_clean_duty_cycle']*100:.1f}%` | `{r['adaptive_n_regimes']}` | `{r['eval_ljung_box_p']:.3f}` |"
        )

    lines.extend(
        [
            "",
            "## 3. Paired Comparison Across the Same 2,000 Evaluation Windows (Frozen Raw Detector)",
            "",
            f"- **False alarm under ALL THREE calibrations (`120w & 200w & 220w`):** `{paired['all_three_true']}` windows",
            f"- **False alarm under ANY calibration (`120w | 200w | 220w`):** `{paired['any_true']}` / `{N_EVAL_WINDOWS}` windows (`{paired['none_true']}` clean under all three)",
            f"- **Exclusive single-condition false alarms:** `only_120 = {paired['only_120']}`, `only_200 = {paired['only_200']}`, `only_220 = {paired['only_220']}`",
            f"- **Pairwise-only overlaps:** `120&200 only = {paired['120_and_200_only']}`, `120&220 only = {paired['120_and_220_only']}`, `200&220 only = {paired['200_and_220_only']}`",
            f"- **Discordant pairs `(A only, B only)`:** `120 vs 200 = {paired['120_vs_200_discordant']}`, `120 vs 220 = {paired['120_vs_220_discordant']}`, `200 vs 220 = {paired['200_vs_220_discordant']}`",
            "",
            "## 4. Statistical Interpretation, Adaptive Metrics & Explicit Scope Limitations",
            "",
            "- **Observed rates vs. true rates & Wilson CI limits:** Reported rates are **observed false-alarm rates** on this 2,000-window evaluation set, not exact true population rates. The 95% Wilson score interval assumes i.i.d. Bernoulli trials conditional on the fixed threshold, is an asymptotic normal-approximation interval, and captures only evaluation sample uncertainty (`N_eval = 2,000`), not calibration-threshold variance across seeds.",
            "- **Generator-based independence vs. Ljung–Box diagnostics:** Independence for `Frozen Raw Detector` rests strictly on the generator design (`generate_pink_colored`, no carried state) and stateless raw detector. Non-significant Ljung–Box p-values test only linear lag-1..10 autocorrelation and do **not** prove statistical independence.",
            "- **Specificity & non-nesting:** Paired false-alarm sets are not nested despite ordered scalar thresholds because each condition uses a different Split-A per-scale standardization (`median` and `Q95 - median`), mapping the same raw 3-scale STFT vector to different standardized scores. Differences across `120w`, `200w`, and `220w` are specific to these three frozen calibrations (`seed = 2026`) and this synthetic evaluation set; they do not imply that larger calibration sets monotonically reduce FAR on every seed.",
            "- **Adaptive full vs. clean metrics:** `Adaptive Full` divides all post-veto nominated windows by `2,000`; `Adaptive Clean` divides nominated windows among clean-calibrated windows by `N_clean` (excluded windows from the rolling 12-window Stream B monitors are never counted as true negatives, and serial dependence in Stream B eligibility makes adaptive CIs approximate).",
            "- **Figures:** `figures/fig01_expanded_neg_control_2000.png`",
            "",
        ]
    )
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")

    payload = {
        "experiment": "expanded_synthetic_negative_control_2000",
        "eval_seed": EVAL_SEED,
        "eval_stream_name": EVAL_STREAM_NAME,
        "n_eval_windows": N_EVAL_WINDOWS,
        "config": cfg.to_dict(),
        "o3_holdout_accessed": False,
        "conditions": cond_rows,
        "paired_frozen_raw": paired,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    (out_dir / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\n=== Expanded Negative Control (N_eval = 2,000 Independent Windows) ===")
    for r in cond_rows:
        print(
            f"  {r['n_bg_windows']}w (th={r['frozen_threshold']:.4f}, {r['calibration_status']}): "
            f"FrozenRaw={r['frozen_raw_fp']}/2000 ({r['frozen_raw_far']*100:.2f}%, 95% CI [{r['frozen_raw_ci95_lo']*100:.2f}%, {r['frozen_raw_ci95_hi']*100:.2f}%]) | "
            f"AdaptFull={r['adaptive_full_fp']}/2000 ({r['adaptive_full_far']*100:.2f}%) | "
            f"AdaptClean={r['adaptive_clean_fp']}/{r['adaptive_n_clean_eligible']} ({r['adaptive_clean_far']*100:.2f}%, duty={r['adaptive_clean_duty_cycle']*100:.1f}%, reg={r['adaptive_n_regimes']})"
        )
    print(f"  Paired (Frozen Raw): all_three={paired['all_three_true']}, any={paired['any_true']}, only_120={paired['only_120']}, only_200={paired['only_200']}, only_220={paired['only_220']}")
    print(f"Elapsed: {payload['elapsed_seconds']} s")


if __name__ == "__main__":
    main()
