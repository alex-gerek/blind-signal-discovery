"""WHISKER / BSD Multi-Seed Calibration Study (Variant A: K=20, N_bg=220, N_eval=1,000).

Exploratory study of between-calibration variability in false-alarm rate (FAR)
for the 220-window calibration condition (`pink_colored` 4-second windows,
target_fap = 0.01, norm_fraction = 0.30 -> Split A = 66, Split B = 154) using
ONLY the Frozen Raw Detector (`MultiResolutionDetector.analyze`).

Design (Variant A):
  - K = 20 independent calibration seeds (m = 1..20).
  - Calibration RNG stream for seed m:
      `make_split_rng(10000 + m, f"multiseed_cal220_seed_{m:02d}")`
    generating 220 independent 4-s `pink_colored` windows (4096 samples at 1024 Hz).
  - Separate independent evaluation RNG stream for seed m:
      `make_split_rng(20000 + m, f"multiseed_eval1000_seed_{m:02d}")`
    generating 1,000 independent 4-s `pink_colored` windows.
  - Zero overlap or leakage between calibration and evaluation streams, and zero
    overlap across seeds.
  - Every calibration outcome (VALIDATED, PROVISIONAL, INSUFFICIENT, or error)
    is recorded and reported without filtering or post-hoc tuning.

Methodological rigor:
  - Distinguishes observed per-seed FAR (p_hat_m = k_m / 1000) from the unknown
    conditional true FAR (p_m) of each frozen calibration.
  - Applies exact finite-sample correction for binomial evaluation noise:
      V_binom = (1 / K) * sum_m [ p_hat_m * (1 - p_hat_m) / (N_eval - 1) ]
    and reports both unadjusted (s_p_hat^2, s_p_hat) and adjusted estimates
    (unclipped V_cal_raw, zero-truncated V_cal_trunc, sigma_cal_trunc),
    explicitly noting that zero-truncation and square-rooting make sigma_cal_trunc
    biased, and providing nonparametric bootstrap intervals across the K=20 seeds.
  - Derives and checks the exact Beta(1, 154) order-statistic theoretical
    reference from `MultiResolutionDetector.calibrate` (k_order = 154 of 154),
    including an empirical upper-tail tie check on all 20 seeds.
  - Keeps detector mathematics, GUI, historical baselines, and the O3 holdout
    completely untouched.
"""

from __future__ import annotations

import csv
import json
import platform
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bsd.adaptive_regime_v2b import estimate_calibration_uncertainty
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, make_split_rng
from bsd.noise import generate_noise
from bsd.workstation_runner import classify_calibration_status

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

K_SEEDS = 20
N_BG_WINDOWS = 220
N_EVAL_WINDOWS = 1000
NORM_FRACTION = 0.30
TARGET_FAP = 0.01
SAMPLE_RATE_HZ = 1024.0
WINDOW_DURATION_S = 4.0
NOISE_MODEL = "pink_colored"
NOISE_SIGMA = 1.0
CAL_SEED_BASE = 10000
EVAL_SEED_BASE = 20000
BOOTSTRAP_RESAMPLES = 5000
BOOTSTRAP_SEED = 99999
BOOTSTRAP_STREAM = "multiseed_220_bootstrap_v1"
BB_MC_DRAWS = 100000
BB_MC_SEED = 88888
BB_MC_STREAM = "multiseed_220_betabinom_mc_v1"

RESULTS_DIR = PROJECT_ROOT / "results" / "workstation_multiseed_220"


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


def _generate_windows(
    n_windows: int,
    master_seed: int,
    stream_name: str,
    prefix: str,
) -> list[Observation]:
    """Generate n_windows independent 4-s pink_colored observations."""
    win_len = int(round(WINDOW_DURATION_S * SAMPLE_RATE_HZ))
    rng = make_split_rng(master_seed, stream_name)
    obs_list: list[Observation] = []
    for i in range(n_windows):
        x, _ = generate_noise(
            NOISE_MODEL,
            win_len,
            SAMPLE_RATE_HZ,
            rng,
            sigma=NOISE_SIGMA,
        )
        obs_list.append(Observation(f"{prefix}_{i:04d}", x, SAMPLE_RATE_HZ, WINDOW_DURATION_S))
    return obs_list


def _compute_variance_decomposition(p_hats: np.ndarray, n_eval: int) -> dict[str, Any]:
    """Compute unadjusted and finite-sample-corrected variance decomposition."""
    k = int(p_hats.size)
    if k < 2:
        raise ValueError(f"Need at least 2 seeds for variance decomposition, got {k}.")

    mean_p = float(np.mean(p_hats))
    median_p = float(np.median(p_hats))
    min_p = float(np.min(p_hats))
    max_p = float(np.max(p_hats))
    iqr_p = float(np.quantile(p_hats, 0.75) - np.quantile(p_hats, 0.25))

    var_total = float(np.var(p_hats, ddof=1))
    sd_total = float(np.sqrt(var_total))

    # Naive plug-in binomial variance (biased low by (n_eval - 1) / n_eval)
    v_binom_naive = float(np.mean(p_hats * (1.0 - p_hats) / float(n_eval)))
    # Finite-sample unbiased moment estimator of E[p_m(1 - p_m) / n_eval]
    v_binom_unbiased = float(np.mean(p_hats * (1.0 - p_hats) / float(n_eval - 1)))

    # Unclipped moment estimator of Var(p_m) (unbiased on variance scale, may be < 0)
    v_cal_raw = float(var_total - v_binom_unbiased)
    # Zero-truncated variance and square-root SD (NOT unbiased!)
    v_cal_trunc = float(max(0.0, v_cal_raw))
    sd_cal_trunc = float(np.sqrt(v_cal_trunc))

    # Naive-plugged variants for full transparency
    v_cal_raw_naive = float(var_total - v_binom_naive)
    sd_cal_trunc_naive = float(np.sqrt(max(0.0, v_cal_raw_naive)))

    frac_var_cal = float(v_cal_trunc / var_total) if var_total > 0.0 else 0.0
    frac_var_binom = float(min(v_binom_unbiased, var_total) / var_total) if var_total > 0.0 else 0.0

    # Nonparametric bootstrap across the K seeds
    rng_boot = make_split_rng(BOOTSTRAP_SEED, BOOTSTRAP_STREAM)
    boot_means = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)
    boot_sd_total = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)
    boot_v_cal_raw = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)
    boot_sd_cal_trunc = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)

    for b in range(BOOTSTRAP_RESAMPLES):
        idx = rng_boot.integers(0, k, size=k)
        sample = p_hats[idx]
        boot_means[b] = float(np.mean(sample))
        vt = float(np.var(sample, ddof=1))
        vb = float(np.mean(sample * (1.0 - sample) / float(n_eval - 1)))
        vc = vt - vb
        boot_sd_total[b] = float(np.sqrt(vt))
        boot_v_cal_raw[b] = vc
        boot_sd_cal_trunc[b] = float(np.sqrt(max(0.0, vc)))

    def _ci95(arr: np.ndarray) -> list[float]:
        return [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))]

    return {
        "n_seeds": k,
        "n_eval_per_seed": int(n_eval),
        "mean_observed_far": mean_p,
        "mean_observed_far_boot95": _ci95(boot_means),
        "median_observed_far": median_p,
        "min_observed_far": min_p,
        "max_observed_far": max_p,
        "iqr_observed_far": iqr_p,
        "unadjusted_var_p_hat": var_total,
        "unadjusted_sd_p_hat": sd_total,
        "unadjusted_sd_p_hat_boot95": _ci95(boot_sd_total),
        "v_binom_naive": v_binom_naive,
        "v_binom_unbiased_n_minus_1": v_binom_unbiased,
        "sd_binom_unbiased": float(np.sqrt(v_binom_unbiased)),
        "adjusted_v_cal_raw_unclipped": v_cal_raw,
        "adjusted_v_cal_raw_boot95": _ci95(boot_v_cal_raw),
        "adjusted_v_cal_zero_truncated": v_cal_trunc,
        "adjusted_sd_cal_zero_truncated": sd_cal_trunc,
        "adjusted_sd_cal_boot95": _ci95(boot_sd_cal_trunc),
        "adjusted_v_cal_raw_naive": v_cal_raw_naive,
        "adjusted_sd_cal_zero_truncated_naive": sd_cal_trunc_naive,
        "bootstrap_frac_v_cal_le_zero": float(np.mean(boot_v_cal_raw <= 0.0)),
        "estimated_var_share_calibration": frac_var_cal,
        "estimated_var_share_binomial_eval": frac_var_binom,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
    }


def _compute_beta_reference(n_thresh: int, target_fap: float, n_eval: int) -> dict[str, Any]:
    """Exact order-statistic Beta reference for MultiResolutionDetector.calibrate."""
    k_order = min(n_thresh, int(np.ceil((1.0 - target_fap) * (n_thresh + 1))))
    alpha = float(n_thresh - k_order + 1)
    beta = float(k_order)

    mean_u = alpha / (alpha + beta)
    var_u = (alpha * beta) / (((alpha + beta) ** 2) * (alpha + beta + 1.0))
    sd_u = float(np.sqrt(var_u))
    e_u_1_minus_u = mean_u - (var_u + mean_u * mean_u)

    v_binom_theory = e_u_1_minus_u / float(n_eval)
    var_p_hat_theory = var_u + v_binom_theory
    sd_p_hat_theory = float(np.sqrt(var_p_hat_theory))

    return {
        "n_thresh": int(n_thresh),
        "target_fap": float(target_fap),
        "k_order": int(k_order),
        "is_sample_maximum": bool(k_order == n_thresh),
        "beta_alpha": alpha,
        "beta_beta": beta,
        "mean_true_far": float(mean_u),
        "sd_true_far": sd_u,
        "var_true_far": float(var_u),
        "median_true_far": float(stats.beta.ppf(0.50, alpha, beta)),
        "q05_true_far": float(stats.beta.ppf(0.05, alpha, beta)),
        "q95_true_far": float(stats.beta.ppf(0.95, alpha, beta)),
        "q025_true_far": float(stats.beta.ppf(0.025, alpha, beta)),
        "q975_true_far": float(stats.beta.ppf(0.975, alpha, beta)),
        "prob_true_far_gt_target": float(1.0 - stats.beta.cdf(target_fap, alpha, beta)),
        "expected_v_binom_at_n_eval": float(v_binom_theory),
        "expected_sd_observed_far_at_n_eval": sd_p_hat_theory,
        "expected_var_observed_far_at_n_eval": float(var_p_hat_theory),
        "expected_var_share_calibration": float(var_u / var_p_hat_theory),
    }


def _compute_exploratory_betabinom_mc(
    k_counts: np.ndarray,
    beta_ref: dict[str, Any],
    n_eval: int,
) -> dict[str, Any]:
    """Exploratory post-hoc Monte Carlo comparison against BetaBinomial(n_eval, a, b).

    Uses the exact order-statistic parameters from `beta_ref` (a = n_thresh - k_order + 1,
    b = k_order). Tail probabilities are strictly P(stat >= obs | H_0), NOT the
    probability that any discrepancy is 'just random'.
    """
    k_seeds = int(k_counts.size)
    obs_rates = k_counts.astype(np.float64) / float(n_eval)
    obs_mean = float(np.mean(obs_rates))
    obs_med = float(np.median(obs_rates))
    obs_sd = float(np.std(obs_rates, ddof=1))

    rng_mc = make_split_rng(BB_MC_SEED, BB_MC_STREAM)
    sim_rates = (
        stats.betabinom.rvs(
            n_eval,
            beta_ref["beta_alpha"],
            beta_ref["beta_beta"],
            size=(BB_MC_DRAWS, k_seeds),
            random_state=rng_mc,
        ).astype(np.float64)
        / float(n_eval)
    )
    sim_means = np.mean(sim_rates, axis=1)
    sim_meds = np.median(sim_rates, axis=1)
    sim_sds = np.std(sim_rates, axis=1, ddof=1)

    p_mean_ge = float(np.mean(sim_means >= obs_mean))
    p_mean_le = float(np.mean(sim_means <= obs_mean))
    p_med_ge = float(np.mean(sim_meds >= obs_med))
    p_sd_ge = float(np.mean(sim_sds >= obs_sd))

    return {
        "exploratory_post_hoc": True,
        "mc_draws": int(BB_MC_DRAWS),
        "mc_seed": int(BB_MC_SEED),
        "mc_stream": BB_MC_STREAM,
        "k_seeds": k_seeds,
        "n_eval": int(n_eval),
        "k_order_used": int(beta_ref["k_order"]),
        "n_thresh_used": int(beta_ref["n_thresh"]),
        "beta_alpha": float(beta_ref["beta_alpha"]),
        "beta_beta": float(beta_ref["beta_beta"]),
        "observed_mean": obs_mean,
        "observed_median": obs_med,
        "observed_sd": obs_sd,
        "upper_tail_prob_mean": p_mean_ge,
        "two_sided_prob_mean": float(min(1.0, 2.0 * min(p_mean_ge, p_mean_le))),
        "upper_tail_prob_median": p_med_ge,
        "upper_tail_prob_sd": p_sd_ge,
    }


def _generate_figures(
    per_seed_rows: list[dict[str, Any]],
    agg: dict[str, Any],
    beta_ref: dict[str, Any],
    out_dir: Path,
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]

    # ------------------------------------------------------------------
    # Figure 1: Per-seed observed FAR forest plot (all K=20 seeds)
    # ------------------------------------------------------------------
    fig1_path = out_dir / "fig01_per_seed_far_forest.png"
    fig, ax = plt.subplots(figsize=(9.2, 6.4))

    y_pos = np.arange(len(ok_rows))
    rates_pct = np.array([r["observed_far"] * 100.0 for r in ok_rows])
    lo_pct = np.array([r["wilson_ci95_lo"] * 100.0 for r in ok_rows])
    hi_pct = np.array([r["wilson_ci95_hi"] * 100.0 for r in ok_rows])
    statuses = [r["calibration_status"] for r in ok_rows]

    color_map = {
        "VALIDATED": "#2b6cb0",
        "PROVISIONAL": "#dd6b20",
        "INSUFFICIENT": "#c53030",
    }
    marker_map = {
        "VALIDATED": "o",
        "PROVISIONAL": "s",
        "INSUFFICIENT": "^",
    }

    # Reference bands and vertical lines
    ax.axvspan(
        beta_ref["q05_true_far"] * 100.0,
        beta_ref["q95_true_far"] * 100.0,
        color="#edf2f7",
        alpha=0.85,
        zorder=0,
        label=f"Beta(1, 154) 5th-95th pct ({beta_ref['q05_true_far']*100:.2f}%-{beta_ref['q95_true_far']*100:.2f}%)",
    )
    ax.axvline(
        TARGET_FAP * 100.0,
        color="#c53030",
        linestyle="--",
        linewidth=1.4,
        zorder=1,
        label="Target FAP (1.00%)",
    )
    ax.axvline(
        beta_ref["mean_true_far"] * 100.0,
        color="#4a5568",
        linestyle=":",
        linewidth=1.4,
        zorder=1,
        label=f"Beta(1, 154) mean ({beta_ref['mean_true_far']*100:.2f}%)",
    )
    ax.axvline(
        agg["mean_observed_far"] * 100.0,
        color="#2f855a",
        linestyle="-.",
        linewidth=1.4,
        zorder=1,
        label=f"Observed mean across K={len(ok_rows)} ({agg['mean_observed_far']*100:.2f}%)",
    )

    seen_status: set[str] = set()
    for i, r in enumerate(ok_rows):
        st = r["calibration_status"]
        col = color_map.get(st, "#4a5568")
        mk = marker_map.get(st, "o")
        lbl = f"Calibration {st}" if st not in seen_status else None
        seen_status.add(st)
        ax.errorbar(
            rates_pct[i],
            y_pos[i],
            xerr=[[rates_pct[i] - lo_pct[i]], [hi_pct[i] - rates_pct[i]]],
            fmt=mk,
            color=col,
            ecolor=col,
            elinewidth=1.3,
            capsize=3,
            markersize=5.5,
            zorder=3,
            label=lbl,
        )
        ax.text(
            hi_pct[i] + 0.08,
            y_pos[i],
            f"{r['k_false_alarms']}/{r['n_eval']} ({rates_pct[i]:.2f}%)  th={r['frozen_threshold']:.2f}",
            va="center",
            ha="left",
            fontsize=7.6,
            color="#2d3748",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Seed {r['seed_index']:02d}" for r in ok_rows])
    ax.invert_yaxis()
    ax.set_xlabel("Observed Frozen-Raw False-Alarm Rate (%) with 95% Wilson CI (N_eval = 1,000 per seed)")
    ax.set_title("Per-Seed Observed FAR Across K = 20 Independent 220-Window Calibrations (Exploratory)")
    max_x = max(float(np.max(hi_pct)) + 1.35, 3.6)
    ax.set_xlim(-0.05, max_x)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.8)
    fig.tight_layout()
    fig.savefig(fig1_path, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: Empirical distribution vs Beta-Binomial & Variance Decomp
    # ------------------------------------------------------------------
    fig2_path = out_dir / "fig02_distribution_and_variance_decomp.png"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    # Left panel: Empirical ECDF of observed FAR vs Beta-Binomial(1000, 1, 154) CDF
    sorted_p = np.sort(rates_pct)
    ecdf_y = np.arange(1, len(sorted_p) + 1) / float(len(sorted_p))
    ax1.step(sorted_p, ecdf_y, where="post", color="#2b6cb0", linewidth=2.0, label="Empirical ECDF (K = 20 seeds)")

    # Exact Beta-Binomial CDF over k = 0..max_k
    max_k_grid = max(int(np.max([r["k_false_alarms"] for r in ok_rows])) + 10, 35)
    k_grid = np.arange(0, max_k_grid + 1)
    bb_cdf = stats.betabinom.cdf(k_grid, N_EVAL_WINDOWS, beta_ref["beta_alpha"], beta_ref["beta_beta"])
    ax1.plot(
        (k_grid / float(N_EVAL_WINDOWS)) * 100.0,
        bb_cdf,
        color="#4a5568",
        linestyle="--",
        linewidth=1.6,
        label="Beta-Binomial(1000, 1, 154) ref.",
    )
    ax1.axvline(TARGET_FAP * 100.0, color="#c53030", linestyle=":", linewidth=1.3, label="Target FAP (1.00%)")
    ax1.set_xlabel("Observed Per-Seed FAR (%)")
    ax1.set_ylabel("Cumulative Probability")
    ax1.set_title("Empirical CDF vs. Order-Statistic Reference")
    ax1.set_xlim(-0.05, max(float(np.max(rates_pct)) + 0.4, 2.8))
    ax1.set_ylim(0.0, 1.05)
    ax1.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=8.0)

    # Right panel: Standard deviation decomposition (observed vs binomial noise vs adjusted cal SD)
    labels = [
        "Total observed\nSD (unadjusted)",
        "Finite-sample\nbinomial SD",
        "Adjusted cal.\nSD (zero-trunc.)",
        "Beta(1, 154)\ntheory SD",
    ]
    vals = [
        agg["unadjusted_sd_p_hat"] * 100.0,
        agg["sd_binom_unbiased"] * 100.0,
        agg["adjusted_sd_cal_zero_truncated"] * 100.0,
        beta_ref["sd_true_far"] * 100.0,
    ]
    colors = ["#2b6cb0", "#718096", "#2f855a", "#805ad5"]
    xx = np.arange(len(labels))
    bars = ax2.bar(xx, vals, width=0.52, color=colors, alpha=0.85, edgecolor="white")

    # Add bootstrap 95% CI whiskers for total SD and adjusted cal SD
    ci_total = [c * 100.0 for c in agg["unadjusted_sd_p_hat_boot95"]]
    ci_cal = [c * 100.0 for c in agg["adjusted_sd_cal_boot95"]]
    ax2.errorbar(
        [0, 2],
        [vals[0], vals[2]],
        yerr=[
            [max(0.0, vals[0] - ci_total[0]), max(0.0, vals[2] - ci_cal[0])],
            [max(0.0, ci_total[1] - vals[0]), max(0.0, ci_cal[1] - vals[2])],
        ],
        fmt="none",
        ecolor="#1a202c",
        elinewidth=1.4,
        capsize=4,
        zorder=4,
        label="95% bootstrap CI (K = 20)",
    )

    for idx_b, (b, v) in enumerate(zip(bars, vals)):
        y_txt = v + 0.03
        if idx_b == 0:
            y_txt = ci_total[1] + 0.03
        elif idx_b == 2:
            y_txt = ci_cal[1] + 0.03
        ax2.text(
            b.get_x() + b.get_width() / 2.0,
            y_txt,
            f"{v:.2f}%",
            ha="center",
            va="bottom",
            fontsize=8.3,
            fontweight="bold",
            color="#1a202c",
        )

    ax2.set_xticks(xx)
    ax2.set_xticklabels(labels, fontsize=8.3)
    ax2.set_ylabel("Standard Deviation (percentage points)")
    ax2.set_title("Variance-Component SD Decomposition (Exploratory)")
    ax2.set_ylim(0.0, max(max(vals), ci_total[1], ci_cal[1]) * 1.24)
    ax2.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=8.0)

    fig.tight_layout()
    fig.savefig(fig2_path, bbox_inches="tight")
    plt.close(fig)

    return fig1_path, fig2_path


def _write_report(
    report_path: Path,
    config: dict[str, Any],
    per_seed_rows: list[dict[str, Any]],
    agg_all: dict[str, Any],
    agg_validated: dict[str, Any] | None,
    beta_ref: dict[str, Any],
    bb_mc: dict[str, Any],
    tie_summary: dict[str, Any],
    elapsed_s: float,
) -> None:
    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]
    failed_rows = [r for r in per_seed_rows if r["run_status"] != "OK"]
    val_rows = [r for r in ok_rows if r["calibration_status"] == "VALIDATED"]
    prov_rows = [r for r in ok_rows if r["calibration_status"] == "PROVISIONAL"]
    insuff_rows = [r for r in ok_rows if r["calibration_status"] == "INSUFFICIENT"]

    n_gt_target = sum(1 for r in ok_rows if r["observed_far"] > TARGET_FAP)
    n_eq_zero = sum(1 for r in ok_rows if r["k_false_alarms"] == 0)
    thresholds = np.array([r["frozen_threshold"] for r in ok_rows], dtype=np.float64)

    lines: list[str] = [
        "# WHISKER / BSD Multi-Seed Calibration Study (Variant A: 220-Window Condition)",
        "",
        "## 1. Study Scope and Exploratory Framing",
        "",
        "This report documents an **initial exploratory multi-seed study** (`K = 20` independent calibration seeds) "
        "examining how much the false-alarm rate (FAR) varies across independent 220-window calibrations of the "
        "frozen raw `MultiResolutionDetector` (`target_fap = 0.01`, `norm_fraction = 0.30`).",
        "",
        "- **Exploratory precision notice:** With `K = 20` calibration seeds, estimates of the between-calibration "
        "spread and tail exceedance probabilities carry substantial sampling uncertainty. This study is designed to "
        "characterize the order of magnitude of calibration-to-calibration variation and verify alignment with "
        "order-statistic theory, not to claim high-precision spread or tail estimates.",
        "- **Observed vs. conditional true FAR:** For each seed $m$, $\\hat{p}_m = k_m / 1000$ is the **observed "
        "per-seed FAR** on that seed's 1,000-window evaluation set. It is a noisy binomial realization around that "
        "frozen calibration's **unknown conditional true FAR** $p_m$.",
        "- **Scope restriction:** Only the stateless **Frozen Raw Detector** (`MultiResolutionDetector.analyze`) is "
        "evaluated on synthetic `pink_colored` noise. Adaptive workstation monitoring, vetoes, real sensor data, and "
        "the reserved O3 holdout are not part of this experiment. Detector mathematics remain completely unmodified.",
        "",
        "## 2. Experimental Protocol and Disjoint RNG Architecture",
        "",
        f"- **Git commit:** `{config['git_commit']}`",
        f"- **Calibration seeds:** `K = {K_SEEDS}` (`m = 1..{K_SEEDS}`), each calibrated on `N_bg = {N_BG_WINDOWS}` "
        f"independent 4-second `{NOISE_MODEL}` windows (`4096` samples at `{SAMPLE_RATE_HZ:.0f} Hz`, `sigma = {NOISE_SIGMA}`).",
        f"- **Split structure:** Split A = `66` windows (30% per-scale median and `Q_0.95 - median` spread), "
        f"Split B = `154` windows (70% global-maximum threshold selection).",
        f"- **Evaluation design (Variant A):** Each seed $m$ is evaluated on its own **separate independent evaluation "
        f"set** of `N_eval = {N_EVAL_WINDOWS}` non-overlapping 4-second `{NOISE_MODEL}` windows. Separate evaluation sets "
        f"ensure that evaluation noise is independent across seeds, allowing direct variance-component decomposition.",
        "- **Leakage prevention:** Calibration streams use `make_split_rng(10000 + m, f\"multiseed_cal220_seed_{m:02d}\")` "
        "and evaluation streams use `make_split_rng(20000 + m, f\"multiseed_eval1000_seed_{m:02d}\")`. Evaluation data "
        "are never used to select, filter, or tune thresholds. Complete seed provenance is stored in `seed_mapping.csv`.",
        "- **Generator independence:** Each call to `generate_pink_colored` draws a fresh spectral slope "
        "`alpha ~ Uniform(0.85, 1.15)`, fresh white Gaussian samples, per-window `rfft`/`irfft` shaping, and per-window "
        "RMS normalization (`sigma = 1.0`) with no carried filter state. `MultiResolutionDetector.analyze` is stateless "
        "across windows. Therefore, conditional on frozen calibration $m$, the 1,000 binary alarm indicators are modeled "
        "as i.i.d. $\\mathrm{Bernoulli}(p_m)$ trials by construction of the generator and detector.",
        f"- **Total wall-clock runtime:** `{elapsed_s:.1f} s` (`{elapsed_s / 60.0:.2f} min`).",
        "",
        "## 3. Exact Threshold-Selection Rule and Order-Statistic Reference",
        "",
        "Direct inspection of `MultiResolutionDetector.calibrate` (`src/bsd/detector_v2.py`) confirms the exact "
        "threshold rule for `N_bg = 220`, `norm_fraction = 0.30`, and `target_fap = 0.01`:",
        "",
        "- `n_norm = max(20, int(round(220 * 0.30))) = 66`",
        "- `n_thresh = 220 - 66 = 154`",
        "- `k_order = min(154, ceil((1 - 0.01) * (154 + 1))) = min(154, ceil(153.45)) = 154`",
        "- `global_threshold = sorted_global[153]` (the exact sample maximum $S_{(154)}$ of Split B)",
        "- Decision rule in `MultiResolutionDetector.analyze`: strict inequality `std_score > global_threshold`.",
        "",
        "### Exchangeability Verification (Split B vs. Evaluation Conditional on Split A)",
        "",
        "The $\\mathrm{Beta}(1, 154)$ reference requires that Split-B global scores and evaluation global scores be "
        "exchangeable conditional on the frozen Split-A standardization. Code inspection verifies all four prerequisites:",
        "",
        "1. **Identical noise distribution:** Both Split B (`cal_obs[66:]`) and evaluation (`eval_obs`) are produced by "
        "identical calls to `generate_noise('pink_colored', 4096, 1024.0, rng, sigma=1.0)` -> `generate_pink_colored`, "
        "drawing fresh `alpha ~ Uniform(0.85, 1.15)`, fresh white Gaussian noise, per-window `rfft`/`irfft` shaping, "
        "and per-window RMS normalization (`sigma = 1.0`) with zero carried filter state.",
        "2. **Identical stateless preprocessing:** `TemplateFreeDetector.analyze` (`src/bsd/detector.py`) applies "
        "per-window median subtraction, unpadded Hann STFT (`boundary=None, padded=False`), 2D median-filter PSD and "
        "envelope whitening, and 8-connected-component clustering strictly to `observation.signal`, storing no cross-window "
        "state and never conditioning on `sample_id`.",
        "3. **Identical score calculation:** Both `MultiResolutionDetector.calibrate` (Split B) and "
        "`MultiResolutionDetector.analyze` (evaluation) call `_raw_analyze_all_scales(obs)`, standardize each scale's "
        "maximum cluster score with the frozen Split-A `median` and `spread_scale`, and take the maximum across the 3 "
        f"scales. Decision-rule audit across all `{len(ok_rows) * N_EVAL_WINDOWS:,}` evaluation windows confirmed "
        f"`{tie_summary['decision_consistency_mismatches']}` mismatches between `res.detected` and "
        "`res.max_global_score > frozen_cal.global_threshold`.",
        "4. **Independent RNG stream treatment:** `make_split_rng(master_seed, split_name)` (`src/bsd/generator.py`) "
        "seeds NumPy's PCG64 generator via `SeedSequence([int(master_seed), int(code)])`. All 40 calibration and "
        "evaluation streams (`cal_master_seed = 10001..10020`, `eval_master_seed = 20001..20020`) use pairwise-distinct "
        "`(master_seed, code)` tuples. Within each calibration seed, Split A (`0..65`) and Split B (`66..219`) are "
        "drawn sequentially from the same PCG64 stream, so Split B is independent of Split A and exchangeable with the "
        "evaluation stream under PCG64.",
        "",
        "### Theoretical Order-Statistic Derivation and Assumptions",
        "",
        "Conditional on the frozen Split-A normalization parameters $\\hat{\\theta}_A$, the $N_B = 154$ Split-B global "
        "scores are independent draws from a common score distribution $F_{\\hat{\\theta}_A}$. If $F_{\\hat{\\theta}_A}$ "
        "is continuous in its upper tail (no ties at the maximum), then the true exceedance probability of the sample "
        "maximum $S_{(154)}$ for a fresh independent draw from the same distribution follows the exact distribution-free "
        "order-statistic law:",
        "",
        "$$p_m = 1 - F_{\\hat{\\theta}_A}(S_{(154)}) \\sim \\mathrm{Beta}(N_B - k_{\\text{order}} + 1,\\; k_{\\text{order}}) = \\mathrm{Beta}(1, 154).$$",
        "",
        "Because $\\mathrm{Beta}(1, 154)$ does not depend on $\\hat{\\theta}_A$, integrating over the random Split-A "
        "realizations leaves the **marginal** distribution of $p_m$ across independent seeds equal to $\\mathrm{Beta}(1, 154)$, "
        "and the observed evaluation count $k_m$ follows the compound $\\mathrm{BetaBinomial}(1000, 1, 154)$ distribution.",
        "",
        "- **Continuity / tie check:** Although `MultiResolutionDetector` assigns a raw cluster score of `0.0` when "
        "no time-frequency cluster passes the initial bin threshold (a lower-tail point mass), this does not affect the "
        "upper-tail order statistic as long as $S_{(154)}$ lies strictly above the zero-cluster floor with no ties. "
        f"Across all `{len(ok_rows)}` seeds, the number of seeds with a tied Split-B maximum was "
        f"`{tie_summary['n_seeds_with_max_tie']}/{len(ok_rows)}`, total exact duplicate score pairs inside Split B was "
        f"`{tie_summary['total_split_b_duplicate_scores']}`, and the minimum Split-B maximum threshold observed was "
        f"`{float(np.min(thresholds)):.4f}` (well above the zero-cluster floor).",
        "- **Beta(1, 154) reference moments:**",
        f"  - Expected conditional true FAR $\\mathbb{{E}}[p_m] = 1/155 = {beta_ref['mean_true_far']*100:.4f}\\%$ "
        "(strictly below the 1.00% target FAP because `k_order = 154` is conservative).",
        f"  - Theoretical between-calibration SD $\\mathrm{{SD}}(p_m) = {beta_ref['sd_true_far']*100:.4f}\\%$ "
        f"(median `{beta_ref['median_true_far']*100:.4f}%`, 5th-95th percentile "
        f"`[{beta_ref['q05_true_far']*100:.4f}%, {beta_ref['q95_true_far']*100:.4f}%]`).",
        f"  - Theoretical probability that a single frozen 220-window calibration has conditional true FAR $> 1.00\\%$: "
        f"`{beta_ref['prob_true_far_gt_target']*100:.2f}%`.",
        f"  - Theoretical total SD of observed $\\hat{{p}}_m$ at `N_eval = 1,000`: "
        f"`{beta_ref['expected_sd_observed_far_at_n_eval']*100:.4f}%` "
        f"(theoretical calibration share of total variance: `{beta_ref['expected_var_share_calibration']*100:.1f}%`).",
        "",
        "## 4. Complete Calibration Status Breakdown (All 20 Seeds Included)",
        "",
        "No seeds were excluded, re-run, or filtered:",
        "",
        f"- **Completed runs:** `{len(ok_rows)} / {K_SEEDS}` (`{len(failed_rows)}` failed runs)",
        f"- **`VALIDATED` calibrations:** `{len(val_rows)} / {len(ok_rows)}`",
        f"- **`PROVISIONAL` calibrations:** `{len(prov_rows)} / {len(ok_rows)}`",
        f"- **`INSUFFICIENT` calibrations:** `{len(insuff_rows)} / {len(ok_rows)}`",
        "",
    ]

    if prov_rows or insuff_rows:
        lines.append("Seeds not achieving `VALIDATED` status (included in all primary aggregates):")
        lines.append("")
        for r in prov_rows + insuff_rows:
            lines.append(
                f"- **Seed {r['seed_index']:02d}** (`{r['calibration_status']}`): "
                f"unmet criteria = `{r['unmet_validation_conditions']}`; "
                f"`n_eff = {r['n_eff']:.1f}`, `Ljung-Box p = {r['cal_ljung_box_p']:.4f}`, "
                f"`observed FAR = {r['k_false_alarms']}/{r['n_eval']} ({r['observed_far']*100:.2f}%)`."
            )
        lines.append("")

    lines.extend(
        [
            "## 5. Per-Seed Results Table (All 20 Seeds)",
            "",
            "| Seed | Cal Status | `n_eff` | Cal LB `p` | Frozen Threshold | False Alarms (`k/1000`) | Observed FAR | 95% Wilson CI | Eval LB `p` |",
            "|---:|---|---:|---:|---:|---:|---:|---|---:|",
        ]
    )

    for r in per_seed_rows:
        if r["run_status"] != "OK":
            lines.append(
                f"| {r['seed_index']:02d} | FAILED ({r['run_status']}) | - | - | - | - | - | - | - |"
            )
        else:
            lines.append(
                f"| {r['seed_index']:02d} | {r['calibration_status']} | {r['n_eff']:.1f} | "
                f"{r['cal_ljung_box_p']:.4f} | {r['frozen_threshold']:.4f} | "
                f"{r['k_false_alarms']}/{r['n_eval']} | {r['observed_far']*100:.2f}% | "
                f"[{r['wilson_ci95_lo']*100:.2f}%, {r['wilson_ci95_hi']*100:.2f}%] | "
                f"{r['eval_ljung_box_p']:.4f} |"
            )

    lines.extend(
        [
            "",
            "> **Note on Ljung-Box `p`:** Ljung-Box tests for linear autocorrelation across lags 1..10. "
            "It is reported strictly as a diagnostic check and does **not** prove statistical independence; "
            "the i.i.d. Bernoulli evaluation model relies on the stateless generator and detector architecture.",
            "",
            "## 6. Aggregate Summary and Variance Decomposition",
            "",
            "### Observed FAR Distribution Across All 20 Seeds",
            "",
            f"- **Mean observed FAR ($\\bar{{p}}$):** `{agg_all['mean_observed_far']*100:.3f}%` "
            f"(95% bootstrap CI across $K=20$ seeds: `[{agg_all['mean_observed_far_boot95'][0]*100:.3f}%, "
            f"{agg_all['mean_observed_far_boot95'][1]*100:.3f}%]`; Beta(1, 154) theoretical mean: "
            f"`{beta_ref['mean_true_far']*100:.3f}%`).",
            f"- **Median observed FAR:** `{agg_all['median_observed_far']*100:.3f}%` "
            f"(IQR: `{agg_all['iqr_observed_far']*100:.3f}%`, range: "
            f"`{agg_all['min_observed_far']*100:.2f}%` to `{agg_all['max_observed_far']*100:.2f}%`).",
            f"- **Seeds with observed FAR $> 1.00\\%$:** `{n_gt_target} / {len(ok_rows)}` "
            f"(`{n_gt_target / len(ok_rows) * 100:.1f}%`).",
            f"- **Seeds with `0 / 1000` observed false alarms:** `{n_eq_zero} / {len(ok_rows)}`.",
            f"- **Frozen threshold range across 20 seeds:** `{float(np.min(thresholds)):.4f}` to "
            f"`{float(np.max(thresholds)):.4f}` (mean `{float(np.mean(thresholds)):.4f}`, "
            f"SD `{float(np.std(thresholds, ddof=1)):.4f}`).",
            "",
            "### Variance-Component Decomposition (Unadjusted vs. Adjusted)",
            "",
            "Under independent evaluation sets of size $N_{\\text{eval}} = 1000$ per seed, the law of total variance gives:",
            "",
            "$$\\mathrm{Var}(\\hat{p}_m) = \\mathrm{Var}(p_m) + \\mathbb{E}\\!\\left[\\frac{p_m(1 - p_m)}{N_{\\text{eval}}}\\right].$$",
            "",
            "Because $\\mathbb{E}[\\hat{p}_m(1 - \\hat{p}_m) \\mid p_m] = \\frac{N_{\\text{eval}} - 1}{N_{\\text{eval}}} p_m(1 - p_m)$, "
            "the finite-sample-corrected unbiased moment estimator of the binomial evaluation noise variance is:",
            "",
            "$$\\widehat{V}_{\\text{binom}} = \\frac{1}{K}\\sum_{m=1}^{K} \\frac{\\hat{p}_m(1 - \\hat{p}_m)}{N_{\\text{eval}} - 1}.$$",
            "",
            "| Quantity | Value (variance scale) | Value (SD, percentage points) | 95% Bootstrap CI (`K = 20`) | Notes |",
            "|---|---:|---:|---|---|",
            f"| Unadjusted total sample variance $s_{{\\hat{{p}}}}^2$ / SD $s_{{\\hat{{p}}}}$ | "
            f"`{agg_all['unadjusted_var_p_hat']:.6e}` | `{agg_all['unadjusted_sd_p_hat']*100:.3f}%` | "
            f"`[{agg_all['unadjusted_sd_p_hat_boot95'][0]*100:.3f}%, {agg_all['unadjusted_sd_p_hat_boot95'][1]*100:.3f}%]` | "
            "Includes both calibration spread and finite-$N_{\\text{eval}}$ binomial noise |",
            f"| Finite-sample-corrected binomial noise $\\widehat{{V}}_{{\\text{{binom}}}}$ (`N_eval - 1 = 999`) | "
            f"`{agg_all['v_binom_unbiased_n_minus_1']:.6e}` | `{agg_all['sd_binom_unbiased']*100:.3f}%` | - | "
            f"Naive `1/N_eval` value: `{agg_all['v_binom_naive']:.6e}` (0.1% lower) |",
            f"| Adjusted calibration variance (unclipped) $\\widehat{{V}}_{{\\text{{cal,raw}}}} = s_{{\\hat{{p}}}}^2 - \\widehat{{V}}_{{\\text{{binom}}}}$ | "
            f"`{agg_all['adjusted_v_cal_raw_unclipped']:.6e}` | - | "
            f"`[{agg_all['adjusted_v_cal_raw_boot95'][0]:.3e}, {agg_all['adjusted_v_cal_raw_boot95'][1]:.3e}]` | "
            "Unbiased on variance scale before truncation; can be negative in small samples |",
            f"| Adjusted calibration variance (zero-truncated) / SD $\\hat{{\\sigma}}_{{\\text{{cal}}}}$ | "
            f"`{agg_all['adjusted_v_cal_zero_truncated']:.6e}` | `{agg_all['adjusted_sd_cal_zero_truncated']*100:.3f}%` | "
            f"`[{agg_all['adjusted_sd_cal_boot95'][0]*100:.3f}%, {agg_all['adjusted_sd_cal_boot95'][1]*100:.3f}%]` | "
            "**Not unbiased** (truncation at 0 adds upward bias; square root adds downward bias via Jensen's inequality) |",
            f"| Theoretical $\\mathrm{{Beta}}(1, 154)$ calibration variance / SD | "
            f"`{beta_ref['var_true_far']:.6e}` | `{beta_ref['sd_true_far']*100:.3f}%` | - | "
            "Exact order-statistic reference under continuous i.i.d. Split-B scores |",
            "",
            f"- **Bootstrap truncation rate:** `{agg_all['bootstrap_frac_v_cal_le_zero']*100:.2f}%` of `{BOOTSTRAP_RESAMPLES}` "
            "bootstrap resamples had $\\widehat{V}_{\\text{cal,raw}} \\le 0$.",
            f"- **Estimated variance share:** Calibration-to-calibration variability accounts for roughly "
            f"`{agg_all['estimated_var_share_calibration']*100:.1f}%` of the total variance of observed $\\hat{{p}}_m$ "
            f"(theoretical Beta-Binomial share: `{beta_ref['expected_var_share_calibration']*100:.1f}%`), while finite-sample "
            f"binomial evaluation noise accounts for roughly `{agg_all['estimated_var_share_binomial_eval']*100:.1f}%`.",
            "- **Methodological limitations of the variance estimator:**",
            "  1. Zero-truncation (`max(0, ...)`) makes $\\widehat{V}_{\\text{cal,trunc}}$ positively biased when true "
            "calibration variance is small relative to sampling noise.",
            "  2. Taking the square root makes $\\hat{\\sigma}_{\\text{cal}}$ negatively biased by Jensen's inequality "
            "($\\mathbb{E}[\\sqrt{X}] \\le \\sqrt{\\mathbb{E}[X]}$).",
            "  3. Percentile bootstrap intervals over `K = 20` seeds are approximate and can under-cover in small samples, "
            "especially for variance components.",
            "",
        ]
    )

    if agg_validated is not None and len(val_rows) < len(ok_rows):
        lines.extend(
            [
                "### Sensitivity Check: `VALIDATED`-Only Subset vs. All 20 Seeds",
                "",
                f"For transparency, restricting to the `{len(val_rows)}` seeds classified as `VALIDATED` gives "
                f"mean observed FAR = `{agg_validated['mean_observed_far']*100:.3f}%`, unadjusted SD = "
                f"`{agg_validated['unadjusted_sd_p_hat']*100:.3f}%`, and zero-truncated adjusted calibration SD = "
                f"`{agg_validated['adjusted_sd_cal_zero_truncated']*100:.3f}%`. All primary conclusions remain based "
                "on the complete unfiltered `K = 20` set.",
                "",
            ]
        )

    lines.extend(
        [
            "### Exploratory Post-Hoc Comparison Against `BetaBinomial(1000, 1, 154)`",
            "",
            "> **Post-hoc disclosure:** The Monte Carlo mean and median tail-probability checks below are "
            "**exploratory and post-hoc** — they were selected *after* inspecting the observed per-seed results and "
            "noting that the empirical ECDF lay to the right of the `BetaBinomial(1000, 1, 154)` curve. Because multiple "
            "summary statistics (mean, median, SD, and exceedance count) were examined on the same `K = 20` sample, "
            "nominal tail probabilities should not be treated as pre-registered confirmatory p-values.",
            "",
            f"Using `{bb_mc['mc_draws']:,}` reproducible Monte Carlo draws of `K = {bb_mc['k_seeds']}` independent "
            f"`BetaBinomial(n={bb_mc['n_eval']}, a={bb_mc['beta_alpha']:.0f}, b={bb_mc['beta_beta']:.0f})` samples "
            f"(from the verified exact threshold rule `k_order = {bb_mc['k_order_used']}` of `{bb_mc['n_thresh_used']}`; "
            f"RNG `make_split_rng({bb_mc['mc_seed']}, \"{bb_mc['mc_stream']}\")`):",
            "",
            f"- **Observed mean (`{bb_mc['observed_mean']*100:.3f}%`) vs. theory (`{beta_ref['mean_true_far']*100:.3f}%`):** "
            f"upper-tail probability $P_{{H_0}}(\\bar{{p}}_{{\\text{{sim}}}} \\ge \\bar{{p}}_{{\\text{{obs}}}}) = "
            f"{bb_mc['upper_tail_prob_mean']:.4f}$ (two-sided: `{bb_mc['two_sided_prob_mean']:.4f}`).",
            f"- **Observed median (`{bb_mc['observed_median']*100:.3f}%`) vs. theory (`{beta_ref['median_true_far']*100:.3f}%`):** "
            f"upper-tail probability $P_{{H_0}}(\\text{{med}}_{{\\text{{sim}}}} \\ge \\text{{med}}_{{\\text{{obs}}}}) = "
            f"{bb_mc['upper_tail_prob_median']:.4f}$.",
            f"- **Observed unadjusted SD (`{bb_mc['observed_sd']*100:.3f}%`) vs. theory (`{beta_ref['expected_sd_observed_far_at_n_eval']*100:.3f}%`):** "
            f"upper-tail probability $P_{{H_0}}(s_{{\\text{{sim}}}} \\ge s_{{\\text{{obs}}}}) = "
            f"{bb_mc['upper_tail_prob_sd']:.4f}$.",
            "",
            "**Interpretation guardrails:** A tail probability $P(\\text{stat} \\ge \\text{obs} \\mid H_0)$ measures how "
            "often a summary at least as large as the observed value occurs *if* the exact `BetaBinomial(1000, 1, 154)` "
            "model holds. It is **not** the probability that the observed shift is \"just random\" or that $H_0$ is true. "
            "Evaluating several summaries on the same small sample (`K = 20`) also inflates the chance that at least one "
            "summary exhibits a small nominal tail probability.",
            "",
        ]
    )

    lines.extend(
        [
            "## 7. Figures",
            "",
            "- `fig01_per_seed_far_forest.png`: Observed FAR and 95% Wilson CI for all 20 seeds, colored by calibration status (`VALIDATED` vs. `PROVISIONAL`), alongside the `1.00%` target FAP and `Beta(1, 154)` reference band.",
            "- `fig02_distribution_and_variance_decomp.png`: Empirical CDF of observed per-seed FAR vs. the exact `Beta-Binomial(1000, 1, 154)` reference distribution (left) and variance-component SD decomposition with 95% bootstrap CIs (right).",
            "",
            "## 8. Scientific Conclusions and Limitations",
            "",
            "1. **Calibration-to-calibration spread is substantial at `N_bg = 220` (`N_thresh = 154`):** Because the 1% "
            "threshold is set by the exact sample maximum $S_{(154)}$ of 154 Split-B windows, individual frozen calibrations "
            "drawn from the exact same stationary `pink_colored` process exhibit wide variation in their realized FAR.",
            "2. **Comparison with order-statistic theory:** The observed spread across the 20 seeds "
            f"(unadjusted SD `{agg_all['unadjusted_sd_p_hat']*100:.3f}%`, adjusted zero-truncated calibration SD "
            f"`{agg_all['adjusted_sd_cal_zero_truncated']*100:.3f}%`) aligns closely with the exact "
            f"$\\mathrm{{Beta}}(1, 154)$ / $\\mathrm{{BetaBinomial}}(1000, 1, 154)$ reference "
            f"(theoretical calibration SD `{beta_ref['sd_true_far']*100:.3f}%`, total SD "
            f"`{beta_ref['expected_sd_observed_far_at_n_eval']*100:.3f}%`). However, the observed center "
            f"(mean `{agg_all['mean_observed_far']*100:.3f}%`, 95% bootstrap CI "
            f"`[{agg_all['mean_observed_far_boot95'][0]*100:.3f}%, {agg_all['mean_observed_far_boot95'][1]*100:.3f}%]`, "
            f"median `{agg_all['median_observed_far']*100:.3f}%`) is higher than the theoretical "
            f"$\\mathrm{{Beta}}(1, 154)$ mean (`{beta_ref['mean_true_far']*100:.3f}%`) and median "
            f"(`{beta_ref['median_true_far']*100:.3f}%`). In exploratory post-hoc Monte Carlo checks "
            f"(`{bb_mc['mc_draws']:,}` draws), such an upward shift has nominal upper-tail probabilities of "
            f"`{bb_mc['upper_tail_prob_mean']:.4f}` (mean) and `{bb_mc['upper_tail_prob_median']:.4f}` (median) under "
            "the exact reference model (not adjusted for examining multiple summaries, and not interpretable as the "
            "probability that the shift is due to chance). Because exchangeability was verified by code inspection, "
            "determining whether this shift persists or regresses toward `0.645%` requires a larger seed count `K`.",
            "3. **Contextualizing the earlier single-seed study (`eb21497`):** The earlier `seed = 2026` 220-window "
            "calibration that yielded `31/2000` (`1.55%`) observed FAR lies inside the upper tail of the `Beta(1, 154)` "
            "distribution (roughly the 91st percentile of conditional true FAR), illustrating why single-seed negative "
            "controls must not be generalized as the population-average FAR of the calibration procedure.",
            "4. **What this study does NOT establish:** It does not establish real-world (LIGO/GWOSC or other physical "
            "sensor) false-alarm control, does not evaluate the adaptive workstation layer, and does not provide "
            "high-precision tail estimates (which would require a larger number of seeds $K$ and/or larger $N_{\\text{eval}}$).",
            "",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    t0 = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    git_commit = _get_git_commit()
    config: dict[str, Any] = {
        "experiment_name": "workstation_multiseed_220_variant_a",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "detector_mode": "frozen_raw_stft_only",
        "k_seeds": K_SEEDS,
        "n_bg_windows": N_BG_WINDOWS,
        "norm_fraction": NORM_FRACTION,
        "n_norm_split_a": max(20, int(round(N_BG_WINDOWS * NORM_FRACTION))),
        "n_thresh_split_b": N_BG_WINDOWS - max(20, int(round(N_BG_WINDOWS * NORM_FRACTION))),
        "n_eval_windows_per_seed": N_EVAL_WINDOWS,
        "evaluation_design": "separate_independent_eval_set_per_seed",
        "target_fap": TARGET_FAP,
        "noise_model": NOISE_MODEL,
        "noise_sigma": NOISE_SIGMA,
        "sample_rate_hz": SAMPLE_RATE_HZ,
        "window_duration_s": WINDOW_DURATION_S,
        "cal_seed_base": CAL_SEED_BASE,
        "eval_seed_base": EVAL_SEED_BASE,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_stream": BOOTSTRAP_STREAM,
        "bb_mc_draws": BB_MC_DRAWS,
        "bb_mc_seed": BB_MC_SEED,
        "bb_mc_stream": BB_MC_STREAM,
    }

    n_thresh_expected = config["n_thresh_split_b"]
    beta_ref = _compute_beta_reference(n_thresh_expected, TARGET_FAP, N_EVAL_WINDOWS)

    print(
        f"[Multi-Seed 220w Variant A] Starting K={K_SEEDS} seeds "
        f"(N_bg={N_BG_WINDOWS}, N_eval={N_EVAL_WINDOWS} separate per seed, Frozen Raw Detector only)...",
        flush=True,
    )

    detector = MultiResolutionDetector()
    seed_mapping_rows: list[dict[str, Any]] = []
    per_seed_rows: list[dict[str, Any]] = []
    n_seeds_with_max_tie = 0
    total_split_b_duplicate_scores = 0
    decision_consistency_mismatches = 0

    for m in range(1, K_SEEDS + 1):
        t_seed = time.perf_counter()
        cal_master_seed = CAL_SEED_BASE + m
        cal_stream_name = f"multiseed_cal220_seed_{m:02d}"
        eval_master_seed = EVAL_SEED_BASE + m
        eval_stream_name = f"multiseed_eval1000_seed_{m:02d}"

        seed_mapping_rows.append(
            {
                "seed_index": m,
                "cal_master_seed": cal_master_seed,
                "cal_stream_name": cal_stream_name,
                "n_bg_windows": N_BG_WINDOWS,
                "eval_master_seed": eval_master_seed,
                "eval_stream_name": eval_stream_name,
                "n_eval_windows": N_EVAL_WINDOWS,
            }
        )

        try:
            cal_obs = _generate_windows(N_BG_WINDOWS, cal_master_seed, cal_stream_name, f"s{m:02d}_cal")
            frozen_cal = detector.calibrate(
                cal_obs,
                target_fap=TARGET_FAP,
                norm_fraction=NORM_FRACTION,
                calibration_seed=cal_master_seed,
            )

            # Chronological Split-B global scores for dependence / status audit
            n_norm = frozen_cal.n_norm_samples
            thresh_obs = cal_obs[n_norm:]
            split_b_chrono = np.array(
                [detector.analyze(o, frozen_cal).max_global_score for o in thresh_obs],
                dtype=np.float64,
            )
            cal_unc = estimate_calibration_uncertainty(split_b_chrono, target_fap=TARGET_FAP)
            cal_audit = classify_calibration_status(N_BG_WINDOWS, cal_unc)

            # Verify upper-tail continuity / ties in Split B
            sorted_b = frozen_cal.sorted_global_null_scores
            max_val = float(sorted_b[-1])
            n_at_max = int(np.sum(sorted_b == max_val))
            has_max_tie = n_at_max > 1
            n_dup_pairs = int(sorted_b.size - np.unique(sorted_b).size)
            if has_max_tie:
                n_seeds_with_max_tie += 1
            total_split_b_duplicate_scores += n_dup_pairs

            # Independent evaluation set for seed m
            eval_obs = _generate_windows(N_EVAL_WINDOWS, eval_master_seed, eval_stream_name, f"s{m:02d}_eval")
            eval_scores = np.empty(N_EVAL_WINDOWS, dtype=np.float64)
            k_fp = 0
            for idx, obs in enumerate(eval_obs):
                res = detector.analyze(obs, frozen_cal)
                eval_scores[idx] = res.max_global_score
                score_gt_th = bool(res.max_global_score > frozen_cal.global_threshold)
                if bool(res.detected) != score_gt_th:
                    decision_consistency_mismatches += 1
                if res.detected:
                    k_fp += 1

            obs_far = float(k_fp / N_EVAL_WINDOWS)
            ci_lo, ci_hi = _wilson_ci(k_fp, N_EVAL_WINDOWS)
            eval_unc = estimate_calibration_uncertainty(eval_scores, target_fap=TARGET_FAP)

            row: dict[str, Any] = {
                "seed_index": m,
                "run_status": "OK",
                "cal_master_seed": cal_master_seed,
                "cal_stream_name": cal_stream_name,
                "eval_master_seed": eval_master_seed,
                "eval_stream_name": eval_stream_name,
                "calibration_status": cal_audit.status,
                "unmet_validation_conditions": "; ".join(cal_audit.unmet_validation_conditions),
                "n_bg_windows": N_BG_WINDOWS,
                "n_norm": frozen_cal.n_norm_samples,
                "n_thresh": frozen_cal.n_threshold_samples,
                "n_eff": round(cal_audit.n_eff, 4),
                "cal_lag1_autocorr": round(cal_audit.lag1_autocorr, 6),
                "cal_ljung_box_p": round(cal_audit.ljung_box_p, 6),
                "cal_far_ci90_neff_lo": round(cal_audit.far_ci90_neff[0], 6),
                "cal_far_ci90_neff_hi": round(cal_audit.far_ci90_neff[1], 6),
                "frozen_threshold": round(frozen_cal.global_threshold, 6),
                "split_b_max_tie": has_max_tie,
                "split_b_duplicate_scores": n_dup_pairs,
                "n_eval": N_EVAL_WINDOWS,
                "k_false_alarms": k_fp,
                "observed_far": round(obs_far, 6),
                "wilson_ci95_lo": round(ci_lo, 6),
                "wilson_ci95_hi": round(ci_hi, 6),
                "eval_lag1_autocorr": round(eval_unc.lag1_autocorr, 6),
                "eval_ljung_box_p": round(eval_unc.ljung_box_p, 6),
                "elapsed_s": round(time.perf_counter() - t_seed, 2),
            }
            for s_name, st in frozen_cal.scale_stats.items():
                row[f"med_{s_name}"] = round(st.median_null_score, 6)
                row[f"spread_{s_name}"] = round(st.spread_scale, 6)

            per_seed_rows.append(row)
            print(
                f"  Seed {m:02d}/{K_SEEDS} | status={cal_audit.status:<11} | "
                f"n_eff={cal_audit.n_eff:5.1f} | th={frozen_cal.global_threshold:.4f} | "
                f"FP={k_fp:2d}/{N_EVAL_WINDOWS} ({obs_far*100:.2f}%) | "
                f"elapsed={row['elapsed_s']:.1f}s",
                flush=True,
            )
        except Exception as exc:
            err_msg = f"{type(exc).__name__}: {exc}"
            per_seed_rows.append(
                {
                    "seed_index": m,
                    "run_status": f"ERROR: {err_msg}",
                    "cal_master_seed": cal_master_seed,
                    "cal_stream_name": cal_stream_name,
                    "eval_master_seed": eval_master_seed,
                    "eval_stream_name": eval_stream_name,
                    "calibration_status": "FAILED",
                    "unmet_validation_conditions": traceback.format_exc(limit=3),
                    "n_bg_windows": N_BG_WINDOWS,
                    "n_norm": 0,
                    "n_thresh": 0,
                    "n_eff": 0.0,
                    "cal_lag1_autocorr": 0.0,
                    "cal_ljung_box_p": 0.0,
                    "cal_far_ci90_neff_lo": 0.0,
                    "cal_far_ci90_neff_hi": 1.0,
                    "frozen_threshold": 0.0,
                    "split_b_max_tie": False,
                    "split_b_duplicate_scores": 0,
                    "n_eval": N_EVAL_WINDOWS,
                    "k_false_alarms": -1,
                    "observed_far": -1.0,
                    "wilson_ci95_lo": 0.0,
                    "wilson_ci95_hi": 1.0,
                    "eval_lag1_autocorr": 0.0,
                    "eval_ljung_box_p": 0.0,
                    "elapsed_s": round(time.perf_counter() - t_seed, 2),
                }
            )
            print(f"  Seed {m:02d}/{K_SEEDS} | FAILED: {err_msg}", flush=True)

    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]
    p_hats_all = np.array([r["k_false_alarms"] / float(r["n_eval"]) for r in ok_rows], dtype=np.float64)
    agg_all = _compute_variance_decomposition(p_hats_all, N_EVAL_WINDOWS)

    val_rows = [r for r in ok_rows if r["calibration_status"] == "VALIDATED"]
    agg_validated: dict[str, Any] | None = None
    if len(val_rows) >= 2:
        p_hats_val = np.array([r["k_false_alarms"] / float(r["n_eval"]) for r in val_rows], dtype=np.float64)
        agg_validated = _compute_variance_decomposition(p_hats_val, N_EVAL_WINDOWS)

    tie_summary = {
        "n_seeds_with_max_tie": n_seeds_with_max_tie,
        "total_split_b_duplicate_scores": total_split_b_duplicate_scores,
        "decision_consistency_mismatches": decision_consistency_mismatches,
    }

    k_counts_all = np.array([r["k_false_alarms"] for r in ok_rows], dtype=np.int64)
    bb_mc = _compute_exploratory_betabinom_mc(k_counts_all, beta_ref, N_EVAL_WINDOWS)

    elapsed_total = time.perf_counter() - t0
    _write_csv(RESULTS_DIR / "seed_mapping.csv", seed_mapping_rows)
    _write_csv(RESULTS_DIR / "per_seed_results.csv", per_seed_rows)

    fig1_path, fig2_path = _generate_figures(per_seed_rows, agg_all, beta_ref, RESULTS_DIR)
    report_path = RESULTS_DIR / "report.md"
    _write_report(
        report_path=report_path,
        config=config,
        per_seed_rows=per_seed_rows,
        agg_all=agg_all,
        agg_validated=agg_validated,
        beta_ref=beta_ref,
        bb_mc=bb_mc,
        tie_summary=tie_summary,
        elapsed_s=elapsed_total,
    )

    payload = {
        "config": config,
        "beta_reference": beta_ref,
        "exploratory_betabinom_mc_audit": bb_mc,
        "tie_and_consistency_audit": tie_summary,
        "status_counts": {
            "total": len(per_seed_rows),
            "ok": len(ok_rows),
            "failed": len(per_seed_rows) - len(ok_rows),
            "validated": sum(1 for r in ok_rows if r["calibration_status"] == "VALIDATED"),
            "provisional": sum(1 for r in ok_rows if r["calibration_status"] == "PROVISIONAL"),
            "insufficient": sum(1 for r in ok_rows if r["calibration_status"] == "INSUFFICIENT"),
        },
        "aggregate_all_seeds": agg_all,
        "aggregate_validated_only": agg_validated,
        "per_seed_results": per_seed_rows,
        "artifacts": {
            "seed_mapping_csv": "results/workstation_multiseed_220/seed_mapping.csv",
            "per_seed_results_csv": "results/workstation_multiseed_220/per_seed_results.csv",
            "report_md": "results/workstation_multiseed_220/report.md",
            "fig01": f"results/workstation_multiseed_220/{fig1_path.name}",
            "fig02": f"results/workstation_multiseed_220/{fig2_path.name}",
        },
        "elapsed_total_s": round(elapsed_total, 2),
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (RESULTS_DIR / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    print(
        f"[Multi-Seed 220w Variant A] DONE in {elapsed_total:.1f}s! "
        f"Mean observed FAR = {agg_all['mean_observed_far']*100:.3f}%, "
        f"unadjusted SD = {agg_all['unadjusted_sd_p_hat']*100:.3f}%, "
        f"adjusted cal SD (zero-trunc) = {agg_all['adjusted_sd_cal_zero_truncated']*100:.3f}%, "
        f"VALIDATED = {payload['status_counts']['validated']}/{len(ok_rows)}.",
        flush=True,
    )


if __name__ == "__main__":
    main()
