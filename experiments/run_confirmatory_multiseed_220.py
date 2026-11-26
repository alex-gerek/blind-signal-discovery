"""WHISKER / BSD Confirmatory Multi-Seed Calibration Study (K=100, N_bg=220, N_eval=500).

Separate confirmatory synthetic experiment designed to test whether the mean
observed false-alarm rate (FAR) across independent 220-window calibrations of
the frozen raw `MultiResolutionDetector` (`target_fap = 0.01`, `norm_fraction = 0.30`)
differs from the `Beta(1, 154)` / `BetaBinomial(500, 1, 154)` theoretical reference
(expected FAR mu_0 = 1/155 ~= 0.64516%).

Strictly separate from the exploratory 20-seed study (`commit 8c27877`,
seeds 10001..10020 / 20001..20020), which is NOT reused here.

Prespecified Protocol:
  1. Single primary statistic:
       T = p_bar = (1 / K) * sum_{m=1}^{100} (k_m / 500) = S / 50000,
     where S = sum_{m=1}^{100} k_m is the total false-alarm count across all
     K = 100 independent calibrations and N_eval = 500 independent windows/seed.
     Two-sided equal-tailed significance level alpha = 0.05 (0.025 per tail).
  2. Prespecified Monte Carlo null simulation (estimates, not exact analytical values):
       RNG: make_split_rng(77777, "confirm_220_k100_n500_h0_v1"), M = 200,000 draws
       of S_sim = sum_{m=1}^{100} BetaBinomial(500, 1, 154).
     Discrete boundary and tie rule (inclusive, conservative, no randomized ties):
       - s_lo = largest integer with P_hat_H0(S_sim <= s_lo) <= 0.025
       - s_hi = smallest integer with P_hat_H0(S_sim >= s_hi) <= 0.025
       - Reject H_0 iff S_obs <= s_lo or S_obs >= s_hi.
       - Achieved Type-I error is estimated as
           alpha_hat = P_hat_H0(S_sim <= s_lo) + P_hat_H0(S_sim >= s_hi) <= 0.05.
  3. Independent RNG streams (zero overlap with 2026 or 10001..20020):
       - Calibration m = 1..100: make_split_rng(30000 + m, f"confirm_cal220_seed_{m:03d}")
       - Evaluation  m = 1..100: make_split_rng(40000 + m, f"confirm_eval500_seed_{m:03d}")
  4. Technical-failure handling (prespecified):
       - Every calibration outcome (VALIDATED, PROVISIONAL, INSUFFICIENT, FAILED)
         is recorded without selective dropping or seed replacement.
       - Any infrastructure/OS exception marks the run INCOMPLETE_TECHNICAL_FAILURE;
         any detector exception on valid input is logged verbatim as FAILED.
  5. Detector mathematics, GUI, historical baselines, and O3 holdout are untouched.
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

K_SEEDS = 100
N_BG_WINDOWS = 220
N_EVAL_WINDOWS = 500
NORM_FRACTION = 0.30
TARGET_FAP = 0.01
SAMPLE_RATE_HZ = 1024.0
WINDOW_DURATION_S = 4.0
NOISE_MODEL = "pink_colored"
NOISE_SIGMA = 1.0

CAL_SEED_BASE = 30000
EVAL_SEED_BASE = 40000

ALPHA_TWO_SIDED = 0.05
ALPHA_PER_TAIL = 0.025
H0_MC_DRAWS = 200000
H0_MC_SEED = 77777
H0_MC_STREAM = "confirm_220_k100_n500_h0_v1"

BOOTSTRAP_RESAMPLES = 10000
BOOTSTRAP_SEED = 99901
BOOTSTRAP_STREAM = "confirm_220_k100_n500_boot_v1"

RESULTS_DIR = PROJECT_ROOT / "results" / "workstation_confirmatory_k100"


def _get_git_dirty() -> bool | None:
    try:
        out = subprocess.check_output(["git", "status", "--porcelain"], cwd=PROJECT_ROOT, stderr=subprocess.DEVNULL)
        return bool(out.decode("utf-8").strip())
    except Exception:
        return None


def _compute_two_sided_mc_p_and_se(p_le: float, p_ge: float, mc_draws: int) -> tuple[float, float]:
    """Compute two-sided MC p-value `min(1, 2*min(p_le, p_ge))` and proper scaled SE `2*SE(p_one)`."""
    if not (np.isfinite(p_le) and np.isfinite(p_ge)) or not (0.0 <= p_le <= 1.0 and 0.0 <= p_ge <= 1.0):
        raise ValueError(f"p_le and p_ge must be finite probabilities in [0, 1], got ({p_le}, {p_ge}).")
    if mc_draws <= 0:
        raise ValueError(f"mc_draws must be positive, got {mc_draws}.")
    p_one = float(min(p_le, p_ge))
    p_two = float(min(1.0, 2.0 * p_one))
    se_two = float(2.0 * np.sqrt(max(0.0, p_one * (1.0 - p_one)) / float(mc_draws)))
    return p_two, se_two


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
        "prob_true_far_gt_target": float(1.0 - stats.beta.cdf(target_fap, alpha, beta)),
        "expected_v_binom_at_n_eval": float(v_binom_theory),
        "expected_sd_observed_far_at_n_eval": sd_p_hat_theory,
        "expected_var_observed_far_at_n_eval": float(var_p_hat_theory),
        "expected_var_share_calibration": float(var_u / var_p_hat_theory),
    }


def _prespecify_h0_discrete_rule(beta_ref: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray]:
    """Compute prespecified Monte Carlo estimated critical integers under H_0.

    Returns both the protocol specification dict and the simulated integer sum
    array `sim_sums` (length H0_MC_DRAWS) for exact p-value evaluation.
    """
    total_trials = K_SEEDS * N_EVAL_WINDOWS
    rng_h0 = make_split_rng(H0_MC_SEED, H0_MC_STREAM)

    # Chunked simulation keeps peak memory modest (~160 MB total) while exact
    sim_sums = np.empty(H0_MC_DRAWS, dtype=np.int64)
    chunk = 20000
    for start in range(0, H0_MC_DRAWS, chunk):
        end = min(H0_MC_DRAWS, start + chunk)
        draws = stats.betabinom.rvs(
            N_EVAL_WINDOWS,
            beta_ref["beta_alpha"],
            beta_ref["beta_beta"],
            size=(end - start, K_SEEDS),
            random_state=rng_h0,
        )
        sim_sums[start:end] = np.sum(draws, axis=1, dtype=np.int64)

    max_s = int(np.max(sim_sums))
    counts = np.bincount(sim_sums, minlength=max_s + 1)
    pmf = counts / float(H0_MC_DRAWS)
    cdf = np.cumsum(pmf)
    sf_inclusive = np.cumsum(pmf[::-1])[::-1]  # P(S >= s)

    # Largest integer s_lo with P_hat(S <= s_lo) <= 0.025
    lo_candidates = np.where(cdf <= ALPHA_PER_TAIL)[0]
    s_lo = int(lo_candidates[-1]) if lo_candidates.size > 0 else -1
    p_lo = float(cdf[s_lo]) if s_lo >= 0 else 0.0

    # Smallest integer s_hi with P_hat(S >= s_hi) <= 0.025
    hi_candidates = np.where(sf_inclusive <= ALPHA_PER_TAIL)[0]
    s_hi = int(hi_candidates[0]) if hi_candidates.size > 0 else total_trials + 1
    p_hi = float(sf_inclusive[s_hi]) if s_hi <= max_s else 0.0

    alpha_est = p_lo + p_hi
    alpha_mc_se = float(np.sqrt(alpha_est * (1.0 - alpha_est) / float(H0_MC_DRAWS)))

    spec = {
        "null_model": f"K={K_SEEDS} i.i.d. BetaBinomial(n={N_EVAL_WINDOWS}, a={beta_ref['beta_alpha']:.0f}, b={beta_ref['beta_beta']:.0f})",
        "total_eval_windows": int(total_trials),
        "lattice_step_far": float(1.0 / total_trials),
        "nominal_alpha_two_sided": ALPHA_TWO_SIDED,
        "nominal_alpha_per_tail": ALPHA_PER_TAIL,
        "mc_draws": H0_MC_DRAWS,
        "mc_seed": H0_MC_SEED,
        "mc_stream": H0_MC_STREAM,
        "critical_values_are_mc_estimates": True,
        "boundary_and_tie_rule": (
            "Inclusive conservative discrete rule on integer total false alarms S = sum(k_m): "
            "s_lo = max{s : P_hat(S <= s) <= 0.025}, s_hi = min{s : P_hat(S >= s) <= 0.025}. "
            "Reject H_0 iff S_obs <= s_lo or S_obs >= s_hi. No randomized tie-breaking."
        ),
        "mc_estimated_s_lo_inclusive": s_lo,
        "mc_estimated_s_hi_inclusive": s_hi,
        "mc_estimated_c_lo_far": float(s_lo / total_trials),
        "mc_estimated_c_hi_far": float(s_hi / total_trials),
        "mc_estimated_lower_tail_prob": p_lo,
        "mc_estimated_upper_tail_prob": p_hi,
        "mc_estimated_achieved_type1_error": alpha_est,
        "mc_se_achieved_type1_error": alpha_mc_se,
        "mc_null_mean_far": float(np.mean(sim_sums) / total_trials),
        "mc_null_sd_of_mean_far": float(np.std(sim_sums, ddof=1) / total_trials),
    }
    return spec, sim_sums


def _compute_variance_decomposition(p_hats: np.ndarray, n_eval: int) -> dict[str, Any]:
    """Compute unadjusted and finite-sample-corrected variance decomposition."""
    k = int(p_hats.size)
    mean_p = float(np.mean(p_hats))
    median_p = float(np.median(p_hats))
    min_p = float(np.min(p_hats))
    max_p = float(np.max(p_hats))
    iqr_p = float(np.quantile(p_hats, 0.75) - np.quantile(p_hats, 0.25))

    var_total = float(np.var(p_hats, ddof=1))
    sd_total = float(np.sqrt(var_total))
    se_mean = float(sd_total / np.sqrt(k))
    t_crit = float(stats.t.ppf(0.975, df=k - 1))
    t_ci95 = [mean_p - t_crit * se_mean, mean_p + t_crit * se_mean]

    v_binom_naive = float(np.mean(p_hats * (1.0 - p_hats) / float(n_eval)))
    v_binom_unbiased = float(np.mean(p_hats * (1.0 - p_hats) / float(n_eval - 1)))
    v_cal_raw = float(var_total - v_binom_unbiased)
    v_cal_trunc = float(max(0.0, v_cal_raw))
    sd_cal_trunc = float(np.sqrt(v_cal_trunc))

    rng_boot = make_split_rng(BOOTSTRAP_SEED, BOOTSTRAP_STREAM)
    boot_means = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)
    boot_sd_total = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)
    boot_sd_cal = np.empty(BOOTSTRAP_RESAMPLES, dtype=np.float64)

    for b in range(BOOTSTRAP_RESAMPLES):
        idx = rng_boot.integers(0, k, size=k)
        s = p_hats[idx]
        boot_means[b] = float(np.mean(s))
        vt = float(np.var(s, ddof=1))
        vb = float(np.mean(s * (1.0 - s) / float(n_eval - 1)))
        boot_sd_total[b] = float(np.sqrt(vt))
        boot_sd_cal[b] = float(np.sqrt(max(0.0, vt - vb)))

    def _ci95(arr: np.ndarray) -> list[float]:
        return [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))]

    return {
        "n_seeds": k,
        "n_eval_per_seed": int(n_eval),
        "mean_observed_far": mean_p,
        "se_mean_observed_far": se_mean,
        "mean_observed_far_t_ci95": t_ci95,
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
        "adjusted_v_cal_zero_truncated": v_cal_trunc,
        "adjusted_sd_cal_zero_truncated": sd_cal_trunc,
        "adjusted_sd_cal_boot95": _ci95(boot_sd_cal),
        "estimated_var_share_calibration": float(v_cal_trunc / var_total) if var_total > 0.0 else 0.0,
        "estimated_var_share_binomial_eval": float(min(v_binom_unbiased, var_total) / var_total) if var_total > 0.0 else 0.0,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
    }


def _generate_figures(
    per_seed_rows: list[dict[str, Any]],
    agg: dict[str, Any],
    beta_ref: dict[str, Any],
    h0_spec: dict[str, Any],
    sim_sums: np.ndarray,
    out_dir: Path,
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]
    total_trials = K_SEEDS * N_EVAL_WINDOWS

    # ------------------------------------------------------------------
    # Figure 1: Primary Confirmatory Test (Null distribution of p_bar vs Observed)
    # ------------------------------------------------------------------
    fig1_path = out_dir / "fig01_confirmatory_primary_test.png"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))

    sim_means_pct = (sim_sums.astype(np.float64) / float(total_trials)) * 100.0
    obs_mean_pct = agg["mean_observed_far"] * 100.0
    mu0_pct = beta_ref["mean_true_far"] * 100.0
    lo_pct = h0_spec["mc_estimated_c_lo_far"] * 100.0
    hi_pct = h0_spec["mc_estimated_c_hi_far"] * 100.0

    bins = np.linspace(float(np.quantile(sim_means_pct, 0.0005)), float(np.quantile(sim_means_pct, 0.9995)), 65)
    ax1.hist(sim_means_pct, bins=bins, density=True, color="#cbd5e0", edgecolor="#a0aec0", alpha=0.85, label="MC H0 null (200k draws)")
    ax1.axvspan(bins[0], lo_pct, color="#feb2b2", alpha=0.35, label=f"Reject region (<= {lo_pct:.3f}% or >= {hi_pct:.3f}%)")
    ax1.axvspan(hi_pct, max(bins[-1], obs_mean_pct + 0.08), color="#feb2b2", alpha=0.35)
    ax1.axvline(mu0_pct, color="#4a5568", linestyle="--", linewidth=1.5, label=f"Beta(1, 154) mean ({mu0_pct:.3f}%)")
    ax1.axvline(obs_mean_pct, color="#2b6cb0", linestyle="-", linewidth=2.2, label=f"Observed mean ({obs_mean_pct:.3f}%)")
    ax1.set_xlim(min(bins[0], obs_mean_pct - 0.05), max(bins[-1], obs_mean_pct + 0.05))
    ax1.set_xlabel("Mean Observed FAR Across K = 100 Seeds (%)")
    ax1.set_ylabel("Density")
    ax1.set_title("Primary Confirmatory Test Against BetaBinomial(500, 1, 154)")
    ax1.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.6)

    # Right panel: ECDF of the 100 per-seed observed FARs vs BetaBinomial(500, 1, 154) CDF
    rates_pct = np.sort(np.array([r["observed_far"] * 100.0 for r in ok_rows]))
    ecdf_y = np.arange(1, len(rates_pct) + 1) / float(len(rates_pct))
    ax2.step(rates_pct, ecdf_y, where="post", color="#2b6cb0", linewidth=2.0, label=f"Empirical ECDF (K = {len(ok_rows)})")

    max_k = max(int(np.max([r["k_false_alarms"] for r in ok_rows])) + 5, 22)
    k_grid = np.arange(0, max_k + 1)
    bb_cdf = stats.betabinom.cdf(k_grid, N_EVAL_WINDOWS, beta_ref["beta_alpha"], beta_ref["beta_beta"])
    ax2.plot((k_grid / float(N_EVAL_WINDOWS)) * 100.0, bb_cdf, color="#4a5568", linestyle="--", linewidth=1.6, label="Beta-Binomial(500, 1, 154) CDF")
    ax2.axvline(TARGET_FAP * 100.0, color="#c53030", linestyle=":", linewidth=1.3, label="Target FAP (1.00%)")
    ax2.set_xlabel("Per-Seed Observed FAR (%)")
    ax2.set_ylabel("Cumulative Probability")
    ax2.set_title("Per-Seed Observed FAR ECDF vs. Reference")
    ax2.set_xlim(-0.05, max(float(np.max(rates_pct)) + 0.25, 3.0))
    ax2.set_ylim(0.0, 1.05)
    ax2.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.8)

    fig.tight_layout()
    fig.savefig(fig1_path, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: Per-Seed Observed FAR Scatter by Status & Frozen Threshold
    # ------------------------------------------------------------------
    fig2_path = out_dir / "fig02_per_seed_distribution_and_status.png"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))

    color_map = {"VALIDATED": "#2b6cb0", "PROVISIONAL": "#dd6b20", "INSUFFICIENT": "#c53030"}
    marker_map = {"VALIDATED": "o", "PROVISIONAL": "s", "INSUFFICIENT": "^"}

    for st in ("VALIDATED", "PROVISIONAL", "INSUFFICIENT"):
        sub = [r for r in ok_rows if r["calibration_status"] == st]
        if not sub:
            continue
        xs = [r["seed_index"] for r in sub]
        ys = [r["observed_far"] * 100.0 for r in sub]
        ax1.scatter(xs, ys, s=26, color=color_map[st], marker=marker_map[st], alpha=0.85, label=f"{st} (n={len(sub)})")

    ax1.axhline(TARGET_FAP * 100.0, color="#c53030", linestyle="--", linewidth=1.2, label="Target FAP (1.00%)")
    ax1.axhline(mu0_pct, color="#4a5568", linestyle=":", linewidth=1.3, label=f"Beta(1, 154) mean ({mu0_pct:.3f}%)")
    ax1.axhline(obs_mean_pct, color="#2f855a", linestyle="-.", linewidth=1.3, label=f"Observed mean ({obs_mean_pct:.3f}%)")
    ax1.set_xlabel("Confirmatory Seed Index (m = 1..100)")
    ax1.set_ylabel("Observed Per-Seed FAR (%)")
    ax1.set_title("Observed FAR Across All 100 Independent Calibrations")
    ax1.set_xlim(0, K_SEEDS + 1)
    ax1.set_ylim(-0.1, max(float(np.max(rates_pct)) + 0.4, 3.2))
    ax1.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.5)

    # Right panel: Frozen threshold vs Observed FAR
    for st in ("VALIDATED", "PROVISIONAL", "INSUFFICIENT"):
        sub = [r for r in ok_rows if r["calibration_status"] == st]
        if not sub:
            continue
        xs = [r["frozen_threshold"] for r in sub]
        ys = [r["observed_far"] * 100.0 for r in sub]
        ax2.scatter(xs, ys, s=26, color=color_map[st], marker=marker_map[st], alpha=0.85, label=f"{st} (n={len(sub)})")

    ax2.axhline(TARGET_FAP * 100.0, color="#c53030", linestyle="--", linewidth=1.2, label="Target FAP (1.00%)")
    ax2.set_xlabel("Frozen Global Threshold S_(154)")
    ax2.set_ylabel("Observed Per-Seed FAR (%)")
    ax2.set_title("Frozen Threshold vs. Realized Observed FAR")
    ax2.set_ylim(-0.1, max(float(np.max(rates_pct)) + 0.4, 3.2))
    ax2.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=7.5)

    fig.tight_layout()
    fig.savefig(fig2_path, bbox_inches="tight")
    plt.close(fig)

    return fig1_path, fig2_path


def _write_protocol_files(protocol: dict[str, Any], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    h0 = protocol["h0_discrete_decision_rule"]
    lines = [
        "# Confirmatory Protocol (Prespecified Before Execution)",
        "",
        f"- **Recorded UTC:** `{protocol['timestamp_utc']}`",
        f"- **Git commit:** `{protocol['git_commit']}`",
        f"- **Design:** `K = {K_SEEDS}` independent 220-window calibrations (`Split A = 66`, `Split B = 154`), "
        f"`N_eval = {N_EVAL_WINDOWS}` independent 4-s `{NOISE_MODEL}` evaluation windows per calibration, "
        "Frozen Raw Detector (`MultiResolutionDetector.analyze`) only.",
        "- **Separation from exploratory data:** The earlier 20-seed study (`commit 8c27877`, seeds `10001..10020` / "
        "`20001..20020`) is strictly exploratory and is excluded from all confirmatory calculations.",
        "- **RNG streams:**",
        f"  - Calibration: `make_split_rng({CAL_SEED_BASE} + m, f\"confirm_cal220_seed_{{m:03d}}\")` (`m = 1..{K_SEEDS}`)",
        f"  - Evaluation: `make_split_rng({EVAL_SEED_BASE} + m, f\"confirm_eval{N_EVAL_WINDOWS}_seed_{{m:03d}}\")` (`m = 1..{K_SEEDS}`)",
        f"  - Verified pairwise-distinct `(master_seed, code)` tuples: `{protocol['unique_rng_stream_pairs']} / {2 * K_SEEDS}`.",
        "- **Single primary statistic:** Mean observed FAR across all 100 seeds, "
        "$T = \\bar{p} = \\frac{1}{100}\\sum_{m=1}^{100} \\frac{k_m}{500} = S / 50000$, tested two-sided at $\\alpha = 0.05$.",
        "- **Prespecified Monte Carlo null simulation (estimates, not exact values):**",
        f"  - Model: `{h0['null_model']}`",
        f"  - RNG: `make_split_rng({h0['mc_seed']}, \"{h0['mc_stream']}\")`, `M = {h0['mc_draws']:,}` draws",
        f"  - Boundary/tie rule: `{h0['boundary_and_tie_rule']}`",
        f"  - MC-estimated critical integers: `s_lo = {h0['mc_estimated_s_lo_inclusive']}` (`{h0['mc_estimated_c_lo_far']*100:.4f}%`), "
        f"`s_hi = {h0['mc_estimated_s_hi_inclusive']}` (`{h0['mc_estimated_c_hi_far']*100:.4f}%`)",
        f"  - MC-estimated achieved Type-I error: `{h0['mc_estimated_achieved_type1_error']:.5f}` "
        f"(lower tail `{h0['mc_estimated_lower_tail_prob']:.5f}` + upper tail `{h0['mc_estimated_upper_tail_prob']:.5f}`, "
        f"MC SE `+/- {h0['mc_se_achieved_type1_error']:.5f}`).",
        "- **Power scope caveat:** Pre-study power estimates (`51.9%` at `0.80%`, `72.4%` at `0.85%`, `85.9%` at `0.90%`, "
        "`96.7%` at `0.98%`) apply specifically to the proposed `Beta(1, (1 - mu_1)/mu_1)` alternative model and stated "
        "effect sizes; they are not a general guarantee for arbitrary non-Beta deviations.",
        "- **Technical-failure handling:** Every seed outcome (`VALIDATED`, `PROVISIONAL`, `INSUFFICIENT`, `FAILED`) "
        "is recorded without selective dropping or replacement. Infrastructure exceptions mark the run "
        "`INCOMPLETE_TECHNICAL_FAILURE`; detector exceptions on valid input are logged verbatim as `FAILED`.",
        "",
    ]
    (out_dir / "protocol.md").write_text("\n".join(lines), encoding="utf-8")


def _write_report(
    report_path: Path,
    protocol: dict[str, Any],
    primary_test: dict[str, Any],
    per_seed_rows: list[dict[str, Any]],
    agg_all: dict[str, Any],
    agg_validated: dict[str, Any] | None,
    beta_ref: dict[str, Any],
    tie_summary: dict[str, Any],
    elapsed_s: float,
) -> None:
    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]
    failed_rows = [r for r in per_seed_rows if r["run_status"] != "OK"]
    val_rows = [r for r in ok_rows if r["calibration_status"] == "VALIDATED"]
    prov_rows = [r for r in ok_rows if r["calibration_status"] == "PROVISIONAL"]
    insuff_rows = [r for r in ok_rows if r["calibration_status"] == "INSUFFICIENT"]
    h0 = protocol["h0_discrete_decision_rule"]

    n_gt_target = sum(1 for r in ok_rows if r["observed_far"] > TARGET_FAP)
    n_eq_zero = sum(1 for r in ok_rows if r["k_false_alarms"] == 0)
    thresholds = np.array([r["frozen_threshold"] for r in ok_rows], dtype=np.float64)

    decision_str = primary_test["decision"]
    delta_pp = primary_test["observed_diff_from_mu0"] * 100.0
    delta_boot_lo = (agg_all["mean_observed_far_boot95"][0] - beta_ref["mean_true_far"]) * 100.0
    delta_boot_hi = (agg_all["mean_observed_far_boot95"][1] - beta_ref["mean_true_far"]) * 100.0
    delta_t_lo = (agg_all["mean_observed_far_t_ci95"][0] - beta_ref["mean_true_far"]) * 100.0
    delta_t_hi = (agg_all["mean_observed_far_t_ci95"][1] - beta_ref["mean_true_far"]) * 100.0

    lines: list[str] = [
        "# WHISKER / BSD Confirmatory Multi-Seed Study (K = 100, N_bg = 220, N_eval = 500)",
        "",
        "## 1. Prespecified Confirmatory Objective and Separation from Exploratory Data",
        "",
        "This report documents a **separate confirmatory synthetic experiment** (`K = 100` independent 220-window "
        "calibrations, `N_eval = 500` independent evaluation windows per calibration, `50,000` total evaluation windows) "
        "evaluating the Frozen Raw `MultiResolutionDetector` on synthetic `pink_colored` noise.",
        "",
        "- **Strict separation:** The earlier 20-seed study (`commit 8c27877`, seeds `10001..10020` / `20001..20020`) "
        "is strictly exploratory and was **not** reused in any confirmatory calculation.",
        "- **Prespecified protocol:** Recorded in `protocol.json` and `protocol.md` prior to seed execution.",
        "- **Unmodified boundaries:** Detector mathematics, GUI, historical baselines, and the reserved O3 holdout were "
        "completely untouched.",
        f"- **Total wall-clock runtime:** `{elapsed_s:.1f} s` (`{elapsed_s / 60.0:.2f} min`).",
        "",
        "## 2. Prespecified Primary Test and Monte Carlo Null Specification",
        "",
        "- **Verified exact threshold rule:** For `N_bg = 220`, `norm_fraction = 0.30`, `target_fap = 0.01`, "
        "`MultiResolutionDetector.calibrate` sets `n_norm = 66`, `n_thresh = 154`, and "
        "`k_order = min(154, ceil(0.99 * 155)) = 154` (exact Split-B sample maximum $S_{(154)}$), with strict inequality "
        "`std_score > global_threshold` in `analyze`. Under conditional exchangeability and upper-tail continuity, "
        "$p_m \\sim \\mathrm{Beta}(1, 154)$ ($\\mu_0 = 1/155 = 0.64516\\%$) and $k_m \\sim \\mathrm{BetaBinomial}(500, 1, 154)$.",
        "- **Single primary statistic:** Mean observed FAR across all 100 calibrations, "
        "$T = \\bar{p} = S / 50000$, where $S = \\sum_{m=1}^{100} k_m$ is the integer total false-alarm count.",
        "- **Two-sided discrete decision rule (Monte Carlo estimates, not exact analytical values):**",
        f"  - Null simulation: `M = {h0['mc_draws']:,}` draws via `make_split_rng({h0['mc_seed']}, \"{h0['mc_stream']}\")`.",
        f"  - Inclusive conservative boundaries on integer $S$: `s_lo = {h0['mc_estimated_s_lo_inclusive']}` "
        f"(`{h0['mc_estimated_c_lo_far']*100:.4f}%`) and `s_hi = {h0['mc_estimated_s_hi_inclusive']}` "
        f"(`{h0['mc_estimated_c_hi_far']*100:.4f}%`). Reject $H_0$ iff $S_{{\\text{{obs}}}} \\le {h0['mc_estimated_s_lo_inclusive']}$ "
        f"or $S_{{\\text{{obs}}}} \\ge {h0['mc_estimated_s_hi_inclusive']}$. No randomized tie-breaking.",
        f"  - MC-estimated achieved Type-I error: `{h0['mc_estimated_achieved_type1_error']:.5f}` "
        f"(lower tail `{h0['mc_estimated_lower_tail_prob']:.5f}` + upper tail `{h0['mc_estimated_upper_tail_prob']:.5f}`, "
        f"MC SE `+/- {h0['mc_se_achieved_type1_error']:.5f}`).",
        "- **Power scope caveat:** Pre-study power estimates (`51.9%` at `0.80%`, `72.4%` at `0.85%`, `85.9%` at `0.90%`, "
        "`96.7%` at `0.98%`) apply specifically to the proposed `Beta(1, (1 - mu_1)/mu_1)` alternative model and stated "
        "effect sizes; they are not a general guarantee for arbitrary deviations.",
        "",
        "## 3. Primary Confirmatory Result",
        "",
        "| Primary Quantity | Value |",
        "|---|---|",
        f"| Total evaluation windows (`K x N_eval`) | `{primary_test['total_eval_windows']:,}` (`{K_SEEDS} x {N_EVAL_WINDOWS}`) |",
        f"| Observed total false alarms ($S_{{\\text{{obs}}}}$) | `{primary_test['observed_total_false_alarms']}` |",
        f"| Observed mean FAR ($\\bar{{p}}$) | `{primary_test['observed_mean_far']*100:.4f}%` |",
        f"| Reference $\\mathrm{{Beta}}(1, 154)$ mean ($\\mu_0 = 1/155$) | `{beta_ref['mean_true_far']*100:.4f}%` |",
        f"| Observed difference ($\\bar{{p}} - \\mu_0$) | `{delta_pp:+.4f}` percentage points |",
        f"| 95% nonparametric bootstrap CI for $\\bar{{p}}$ (`{BOOTSTRAP_RESAMPLES:,}` resamples) | "
        f"`[{agg_all['mean_observed_far_boot95'][0]*100:.4f}%, {agg_all['mean_observed_far_boot95'][1]*100:.4f}%]` "
        f"(diff: `[{delta_boot_lo:+.4f}, {delta_boot_hi:+.4f}]` pp) |",
        f"| 95% t-interval for $\\bar{{p}}$ (`df = {K_SEEDS - 1}`) | "
        f"`[{agg_all['mean_observed_far_t_ci95'][0]*100:.4f}%, {agg_all['mean_observed_far_t_ci95'][1]*100:.4f}%]` "
        f"(diff: `[{delta_t_lo:+.4f}, {delta_t_hi:+.4f}]` pp) |",
        f"| MC-estimated critical region ($\\alpha = 0.05$, inclusive) | "
        f"`S <= {h0['mc_estimated_s_lo_inclusive']}` (`<= {h0['mc_estimated_c_lo_far']*100:.4f}%`) or "
        f"`S >= {h0['mc_estimated_s_hi_inclusive']}` (`>= {h0['mc_estimated_c_hi_far']*100:.4f}%`) |",
        f"| MC-estimated lower / upper tail probabilities | "
        f"`P_hat(S <= S_obs) = {primary_test['mc_estimated_p_le']:.5f}`, "
        f"`P_hat(S >= S_obs) = {primary_test['mc_estimated_p_ge']:.5f}` |",
        f"| MC-estimated two-sided p-value | `{primary_test['mc_estimated_p_two_sided']:.5f}` "
        f"(MC SE `+/- {primary_test['mc_se_p_two_sided']:.5f}`) |",
        f"| **Prespecified primary decision** | **`{decision_str}`** |",
        "",
        "## 4. Complete Calibration Status Breakdown and Audit Checks (All 100 Seeds)",
        "",
        "No seeds were dropped, replaced, or filtered:",
        "",
        f"- **Completed runs:** `{len(ok_rows)} / {K_SEEDS}` (`{len(failed_rows)}` failed runs)",
        f"- **`VALIDATED` calibrations:** `{len(val_rows)} / {len(ok_rows)}` (`{len(val_rows)/len(ok_rows)*100:.1f}%`)",
        f"- **`PROVISIONAL` calibrations:** `{len(prov_rows)} / {len(ok_rows)}` (`{len(prov_rows)/len(ok_rows)*100:.1f}%`)",
        f"- **`INSUFFICIENT` calibrations:** `{len(insuff_rows)} / {len(ok_rows)}` (`{len(insuff_rows)/len(ok_rows)*100:.1f}%`)",
        f"- **Split-B continuity / consistency audit:** `{tie_summary['n_seeds_with_max_tie']}/{len(ok_rows)}` seeds "
        f"with a tied Split-B maximum, `{tie_summary['total_split_b_duplicate_scores']}` total Split-B duplicate scores, "
        f"`{tie_summary['decision_consistency_mismatches']}` mismatches between `res.detected` and "
        "`res.max_global_score > frozen_cal.global_threshold` across all `50,000` evaluation windows.",
        f"- **RNG isolation check:** `{protocol['unique_rng_stream_pairs']}/{2 * K_SEEDS}` pairwise-distinct "
        "`(master_seed, code)` tuples across all calibration and evaluation streams.",
        "",
        "Complete per-seed results (all 100 rows, including unmet validation criteria for every non-`VALIDATED` seed) "
        "are saved in `per_seed_results.csv` and `seed_mapping.csv`.",
        "",
        "## 5. Secondary Descriptive Summaries (Not Part of Primary Decision Rule)",
        "",
        f"- **Median observed FAR:** `{agg_all['median_observed_far']*100:.3f}%` "
        f"(IQR `{agg_all['iqr_observed_far']*100:.3f}%`, range `{agg_all['min_observed_far']*100:.2f}%` to "
        f"`{agg_all['max_observed_far']*100:.2f}%`; `Beta(1, 154)` theoretical median: `{beta_ref['median_true_far']*100:.3f}%`).",
        f"- **Seeds with observed FAR $> 1.00\\%$:** `{n_gt_target} / {len(ok_rows)}` (`{n_gt_target/len(ok_rows)*100:.1f}%`).",
        f"- **Seeds with `0 / 500` observed false alarms:** `{n_eq_zero} / {len(ok_rows)}`.",
        f"- **Frozen threshold range:** `{float(np.min(thresholds)):.4f}` to `{float(np.max(thresholds)):.4f}` "
        f"(mean `{float(np.mean(thresholds)):.4f}`, SD `{float(np.std(thresholds, ddof=1)):.4f}`).",
        f"- **Variance decomposition (all 100 seeds):**",
        f"  - Unadjusted total SD $s_{{\\hat{{p}}}}$: `{agg_all['unadjusted_sd_p_hat']*100:.3f}%` "
        f"(95% bootstrap CI `[{agg_all['unadjusted_sd_p_hat_boot95'][0]*100:.3f}%, {agg_all['unadjusted_sd_p_hat_boot95'][1]*100:.3f}%]`; "
        f"theoretical Beta-Binomial total SD: `{beta_ref['expected_sd_observed_far_at_n_eval']*100:.3f}%`).",
        f"  - Finite-sample-corrected binomial noise SD (`N_eval - 1 = 499`): `{agg_all['sd_binom_unbiased']*100:.3f}%` "
        f"($\\widehat{{V}}_{{\\text{{binom}}}} = {agg_all['v_binom_unbiased_n_minus_1']:.6e}$).",
        f"  - Adjusted zero-truncated calibration SD $\\hat{{\\sigma}}_{{\\text{{cal}}}}$ (not unbiased due to zero-truncation "
        f"and Jensen's inequality): `{agg_all['adjusted_sd_cal_zero_truncated']*100:.3f}%` "
        f"(95% bootstrap CI `[{agg_all['adjusted_sd_cal_boot95'][0]*100:.3f}%, {agg_all['adjusted_sd_cal_boot95'][1]*100:.3f}%]`; "
        f"theoretical `Beta(1, 154)` SD: `{beta_ref['sd_true_far']*100:.3f}%`).",
    ]

    if agg_validated is not None and len(val_rows) < len(ok_rows):
        lines.append(
            f"- **`VALIDATED`-only descriptive subset (`n = {len(val_rows)}`):** "
            f"mean observed FAR = `{agg_validated['mean_observed_far']*100:.3f}%`, "
            f"unadjusted SD = `{agg_validated['unadjusted_sd_p_hat']*100:.3f}%`, "
            f"adjusted zero-truncated calibration SD = `{agg_validated['adjusted_sd_cal_zero_truncated']*100:.3f}%` "
            "(reported for transparency only; primary inference uses all 100 seeds)."
        )

    lines.extend(
        [
            "",
            "## 6. Figures",
            "",
            "- `fig01_confirmatory_primary_test.png`: Monte Carlo null distribution of $\\bar{p}$ with inclusive $\\alpha = 0.05$ rejection regions vs. observed mean (left), and empirical ECDF of all 100 per-seed FARs vs. `BetaBinomial(500, 1, 154)` CDF (right).",
            "- `fig02_per_seed_distribution_and_status.png`: Per-seed observed FAR across all 100 seeds colored by calibration status (left), and frozen threshold $S_{(154)}$ vs. observed FAR (right).",
            "",
            "## 7. Interpretation and Limitations",
            "",
        ]
    )

    if decision_str == "FAIL_TO_REJECT_H0":
        lines.extend(
            [
                f"1. **Primary interpretation (`FAIL_TO_REJECT_H0`):** The confirmatory mean observed FAR "
                f"(`{primary_test['observed_mean_far']*100:.4f}%`, difference `{delta_pp:+.4f}` pp from `0.6452%`, "
                f"95% bootstrap CI for difference `[{delta_boot_lo:+.4f}, {delta_boot_hi:+.4f}]` pp; MC-estimated "
                f"two-sided `p = {primary_test['mc_estimated_p_two_sided']:.5f}`) does not fall in the prespecified "
                f"$\\alpha = 0.05$ rejection region.",
                "2. **What this does NOT prove:** Failing to reject $H_0$ is **not** proof that $\\mu = 1/155$ exactly "
                "or proof that any discrepancy has disappeared. The confidence interval quantifies the range of residual "
                "mean shifts compatible with these 100 confirmatory calibrations.",
            ]
        )
    else:
        lines.extend(
            [
                f"1. **Primary interpretation (`{decision_str}`):** The confirmatory mean observed FAR "
                f"(`{primary_test['observed_mean_far']*100:.4f}%`, difference `{delta_pp:+.4f}` pp from `0.6452%`, "
                f"95% bootstrap CI for difference `[{delta_boot_lo:+.4f}, {delta_boot_hi:+.4f}]` pp; MC-estimated "
                f"two-sided `p = {primary_test['mc_estimated_p_two_sided']:.5f}`) falls inside the prespecified "
                "$\\alpha = 0.05$ rejection region, providing confirmatory evidence of a systematic deviation from the "
                "`BetaBinomial(500, 1, 154)` reference under this synthetic protocol.",
                "2. **What this does NOT establish:** Rejecting $H_0$ does not by itself identify the physical or "
                "numerical cause of the shift and does not directly characterize performance on real detector data.",
            ]
        )

    lines.extend(
        [
            "3. **Scope limitations:** Results apply strictly to synthetic stationary `pink_colored` noise evaluated with "
            "the stateless Frozen Raw `MultiResolutionDetector` at `N_bg = 220`, `target_fap = 0.01`, and 4-second windows. "
            "They do not establish false-alarm control on real LIGO/GWOSC or other physical sensor data, nor do they "
            "evaluate the adaptive workstation monitoring layer or vetoes.",
            "",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    t0 = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    git_commit = _get_git_commit()
    n_norm_expected = max(20, int(round(N_BG_WINDOWS * NORM_FRACTION)))
    n_thresh_expected = N_BG_WINDOWS - n_norm_expected
    beta_ref = _compute_beta_reference(n_thresh_expected, TARGET_FAP, N_EVAL_WINDOWS)

    # Verify pairwise-distinct (master_seed, code) tuples BEFORE running
    pairs: set[tuple[int, int]] = set()
    seed_mapping_rows: list[dict[str, Any]] = []
    for m in range(1, K_SEEDS + 1):
        cal_ms = CAL_SEED_BASE + m
        cal_sn = f"confirm_cal220_seed_{m:03d}"
        eval_ms = EVAL_SEED_BASE + m
        eval_sn = f"confirm_eval{N_EVAL_WINDOWS}_seed_{m:03d}"
        for ms, sn in ((cal_ms, cal_sn), (eval_ms, eval_sn)):
            code = 9000 + sum((idx + 1) * ord(ch) for idx, ch in enumerate(sn)) % 100000
            pairs.add((ms, code))
        seed_mapping_rows.append(
            {
                "seed_index": m,
                "cal_master_seed": cal_ms,
                "cal_stream_name": cal_sn,
                "n_bg_windows": N_BG_WINDOWS,
                "eval_master_seed": eval_ms,
                "eval_stream_name": eval_sn,
                "n_eval_windows": N_EVAL_WINDOWS,
            }
        )
    if len(pairs) != 2 * K_SEEDS:
        raise RuntimeError(f"RNG stream collision detected: {len(pairs)} unique pairs != {2 * K_SEEDS}")

    print("[Confirmatory K=100, N_eval=500] Computing prespecified MC H0 decision rule (200,000 draws)...", flush=True)
    h0_spec, sim_sums = _prespecify_h0_discrete_rule(beta_ref)

    protocol: dict[str, Any] = {
        "experiment_name": "workstation_confirmatory_k100_n500",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "git_dirty": _get_git_dirty(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "detector_mode": "frozen_raw_stft_only",
        "k_seeds": K_SEEDS,
        "n_bg_windows": N_BG_WINDOWS,
        "norm_fraction": NORM_FRACTION,
        "n_norm_split_a": n_norm_expected,
        "n_thresh_split_b": n_thresh_expected,
        "n_eval_windows_per_seed": N_EVAL_WINDOWS,
        "target_fap": TARGET_FAP,
        "noise_model": NOISE_MODEL,
        "noise_sigma": NOISE_SIGMA,
        "sample_rate_hz": SAMPLE_RATE_HZ,
        "window_duration_s": WINDOW_DURATION_S,
        "cal_seed_base": CAL_SEED_BASE,
        "eval_seed_base": EVAL_SEED_BASE,
        "unique_rng_stream_pairs": len(pairs),
        "exploratory_study_excluded": "commit 8c27877 (seeds 10001..10020 / 20001..20020)",
        "h0_discrete_decision_rule": h0_spec,
        "bootstrap_resamples": BOOTSTRAP_RESAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_stream": BOOTSTRAP_STREAM,
    }
    _write_protocol_files(protocol, RESULTS_DIR)
    _write_csv(RESULTS_DIR / "seed_mapping.csv", seed_mapping_rows)

    print(
        f"[Confirmatory K=100, N_eval=500] Protocol recorded. "
        f"MC critical integers: S <= {h0_spec['mc_estimated_s_lo_inclusive']} "
        f"({h0_spec['mc_estimated_c_lo_far']*100:.4f}%) or "
        f"S >= {h0_spec['mc_estimated_s_hi_inclusive']} "
        f"({h0_spec['mc_estimated_c_hi_far']*100:.4f}%), "
        f"est. Type-I error = {h0_spec['mc_estimated_achieved_type1_error']:.5f}.",
        flush=True,
    )

    detector = MultiResolutionDetector()
    per_seed_rows: list[dict[str, Any]] = []
    n_seeds_with_max_tie = 0
    total_split_b_duplicate_scores = 0
    decision_consistency_mismatches = 0

    for row_map in seed_mapping_rows:
        m = row_map["seed_index"]
        cal_ms = row_map["cal_master_seed"]
        cal_sn = row_map["cal_stream_name"]
        eval_ms = row_map["eval_master_seed"]
        eval_sn = row_map["eval_stream_name"]
        t_seed = time.perf_counter()

        try:
            cal_obs = _generate_windows(N_BG_WINDOWS, cal_ms, cal_sn, f"c{m:03d}_cal")
            frozen_cal = detector.calibrate(
                cal_obs,
                target_fap=TARGET_FAP,
                norm_fraction=NORM_FRACTION,
                calibration_seed=cal_ms,
            )

            n_norm = frozen_cal.n_norm_samples
            thresh_obs = cal_obs[n_norm:]
            split_b_chrono = np.array(
                [detector.analyze(o, frozen_cal).max_global_score for o in thresh_obs],
                dtype=np.float64,
            )
            cal_unc = estimate_calibration_uncertainty(split_b_chrono, target_fap=TARGET_FAP)
            cal_audit = classify_calibration_status(N_BG_WINDOWS, cal_unc)

            sorted_b = frozen_cal.sorted_global_null_scores
            max_val = float(sorted_b[-1])
            has_max_tie = int(np.sum(sorted_b == max_val)) > 1
            n_dup_pairs = int(sorted_b.size - np.unique(sorted_b).size)
            if has_max_tie:
                n_seeds_with_max_tie += 1
            total_split_b_duplicate_scores += n_dup_pairs

            eval_obs = _generate_windows(N_EVAL_WINDOWS, eval_ms, eval_sn, f"c{m:03d}_eval")
            k_fp = 0
            for obs in eval_obs:
                res = detector.analyze(obs, frozen_cal)
                score_gt_th = bool(res.max_global_score > frozen_cal.global_threshold)
                if bool(res.detected) != score_gt_th:
                    decision_consistency_mismatches += 1
                if res.detected:
                    k_fp += 1

            obs_far = float(k_fp / N_EVAL_WINDOWS)
            ci_lo, ci_hi = _wilson_ci(k_fp, N_EVAL_WINDOWS)

            per_seed_rows.append(
                {
                    "seed_index": m,
                    "run_status": "OK",
                    "cal_master_seed": cal_ms,
                    "cal_stream_name": cal_sn,
                    "eval_master_seed": eval_ms,
                    "eval_stream_name": eval_sn,
                    "calibration_status": cal_audit.status,
                    "unmet_validation_conditions": "; ".join(cal_audit.unmet_validation_conditions),
                    "n_bg_windows": N_BG_WINDOWS,
                    "n_norm": frozen_cal.n_norm_samples,
                    "n_thresh": frozen_cal.n_threshold_samples,
                    "n_eff": round(cal_audit.n_eff, 4),
                    "cal_lag1_autocorr": round(cal_audit.lag1_autocorr, 6),
                    "cal_ljung_box_p": round(cal_audit.ljung_box_p, 6),
                    "cal_far_ci90_neff_hi": round(cal_audit.far_ci90_neff[1], 6),
                    "frozen_threshold": round(frozen_cal.global_threshold, 6),
                    "split_b_max_tie": has_max_tie,
                    "split_b_duplicate_scores": n_dup_pairs,
                    "n_eval": N_EVAL_WINDOWS,
                    "k_false_alarms": k_fp,
                    "observed_far": round(obs_far, 6),
                    "wilson_ci95_lo": round(ci_lo, 6),
                    "wilson_ci95_hi": round(ci_hi, 6),
                    "elapsed_s": round(time.perf_counter() - t_seed, 2),
                }
            )
            if m % 10 == 0 or m == 1:
                running_mean = float(np.mean([r["observed_far"] for r in per_seed_rows if r["run_status"] == "OK"]))
                print(
                    f"  Seed {m:03d}/{K_SEEDS} | status={cal_audit.status:<11} | "
                    f"FP={k_fp:2d}/{N_EVAL_WINDOWS} ({obs_far*100:.2f}%) | "
                    f"running_mean={running_mean*100:.3f}% | elapsed={time.perf_counter()-t0:.0f}s",
                    flush=True,
                )
        except Exception as exc:
            err_msg = f"{type(exc).__name__}: {exc}"
            per_seed_rows.append(
                {
                    "seed_index": m,
                    "run_status": f"ERROR: {err_msg}",
                    "cal_master_seed": cal_ms,
                    "cal_stream_name": cal_sn,
                    "eval_master_seed": eval_ms,
                    "eval_stream_name": eval_sn,
                    "calibration_status": "FAILED",
                    "unmet_validation_conditions": traceback.format_exc(limit=3),
                    "n_bg_windows": N_BG_WINDOWS,
                    "n_norm": 0,
                    "n_thresh": 0,
                    "n_eff": 0.0,
                    "cal_lag1_autocorr": 0.0,
                    "cal_ljung_box_p": 0.0,
                    "cal_far_ci90_neff_hi": 1.0,
                    "frozen_threshold": 0.0,
                    "split_b_max_tie": False,
                    "split_b_duplicate_scores": 0,
                    "n_eval": N_EVAL_WINDOWS,
                    "k_false_alarms": -1,
                    "observed_far": -1.0,
                    "wilson_ci95_lo": 0.0,
                    "wilson_ci95_hi": 1.0,
                    "elapsed_s": round(time.perf_counter() - t_seed, 2),
                }
            )
            print(f"  Seed {m:03d}/{K_SEEDS} | FAILED: {err_msg}", flush=True)

    ok_rows = [r for r in per_seed_rows if r["run_status"] == "OK"]
    failed_rows = [r for r in per_seed_rows if r["run_status"] != "OK"]
    if failed_rows:
        raise RuntimeError(
            f"INCOMPLETE_TECHNICAL_FAILURE: {len(failed_rows)} of {K_SEEDS} seeds failed. "
            "Per prespecified protocol, seeds are not selectively dropped."
        )

    k_counts = np.array([r["k_false_alarms"] for r in ok_rows], dtype=np.int64)
    p_hats = k_counts.astype(np.float64) / float(N_EVAL_WINDOWS)
    s_obs = int(np.sum(k_counts))
    total_trials = K_SEEDS * N_EVAL_WINDOWS
    obs_mean = float(s_obs / total_trials)

    p_le = float(np.mean(sim_sums <= s_obs))
    p_ge = float(np.mean(sim_sums >= s_obs))
    p_two_sided, p_two_sided_se = _compute_two_sided_mc_p_and_se(p_le, p_ge, H0_MC_DRAWS)

    if s_obs <= h0_spec["mc_estimated_s_lo_inclusive"]:
        decision = "REJECT_H0_BELOW_REFERENCE"
    elif s_obs >= h0_spec["mc_estimated_s_hi_inclusive"]:
        decision = "REJECT_H0_ABOVE_REFERENCE"
    else:
        decision = "FAIL_TO_REJECT_H0"

    primary_test = {
        "total_eval_windows": total_trials,
        "observed_total_false_alarms": s_obs,
        "observed_mean_far": obs_mean,
        "reference_mu0": beta_ref["mean_true_far"],
        "observed_diff_from_mu0": float(obs_mean - beta_ref["mean_true_far"]),
        "mc_estimated_p_le": p_le,
        "mc_estimated_p_ge": p_ge,
        "mc_estimated_p_two_sided": p_two_sided,
        "mc_se_p_two_sided": p_two_sided_se,
        "decision": decision,
    }

    agg_all = _compute_variance_decomposition(p_hats, N_EVAL_WINDOWS)
    val_rows = [r for r in ok_rows if r["calibration_status"] == "VALIDATED"]
    agg_validated: dict[str, Any] | None = None
    if len(val_rows) >= 2:
        p_hats_val = np.array([r["k_false_alarms"] / float(N_EVAL_WINDOWS) for r in val_rows], dtype=np.float64)
        agg_validated = _compute_variance_decomposition(p_hats_val, N_EVAL_WINDOWS)

    tie_summary = {
        "n_seeds_with_max_tie": n_seeds_with_max_tie,
        "total_split_b_duplicate_scores": total_split_b_duplicate_scores,
        "decision_consistency_mismatches": decision_consistency_mismatches,
    }

    elapsed_total = time.perf_counter() - t0
    _write_csv(RESULTS_DIR / "per_seed_results.csv", per_seed_rows)
    fig1_path, fig2_path = _generate_figures(per_seed_rows, agg_all, beta_ref, h0_spec, sim_sums, RESULTS_DIR)
    _write_report(
        report_path=RESULTS_DIR / "report.md",
        protocol=protocol,
        primary_test=primary_test,
        per_seed_rows=per_seed_rows,
        agg_all=agg_all,
        agg_validated=agg_validated,
        beta_ref=beta_ref,
        tie_summary=tie_summary,
        elapsed_s=elapsed_total,
    )

    payload = {
        "protocol": protocol,
        "beta_reference": beta_ref,
        "primary_test_result": primary_test,
        "tie_and_consistency_audit": tie_summary,
        "status_counts": {
            "total": len(per_seed_rows),
            "ok": len(ok_rows),
            "failed": len(failed_rows),
            "validated": sum(1 for r in ok_rows if r["calibration_status"] == "VALIDATED"),
            "provisional": sum(1 for r in ok_rows if r["calibration_status"] == "PROVISIONAL"),
            "insufficient": sum(1 for r in ok_rows if r["calibration_status"] == "INSUFFICIENT"),
        },
        "aggregate_all_seeds": agg_all,
        "aggregate_validated_only": agg_validated,
        "per_seed_results": per_seed_rows,
        "artifacts": {
            "protocol_json": "results/workstation_confirmatory_k100/protocol.json",
            "protocol_md": "results/workstation_confirmatory_k100/protocol.md",
            "seed_mapping_csv": "results/workstation_confirmatory_k100/seed_mapping.csv",
            "per_seed_results_csv": "results/workstation_confirmatory_k100/per_seed_results.csv",
            "report_md": "results/workstation_confirmatory_k100/report.md",
            "fig01": f"results/workstation_confirmatory_k100/{fig1_path.name}",
            "fig02": f"results/workstation_confirmatory_k100/{fig2_path.name}",
        },
        "elapsed_total_s": round(elapsed_total, 2),
    }
    (RESULTS_DIR / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(
        f"[Confirmatory K=100, N_eval=500] DONE in {elapsed_total:.1f}s! "
        f"S_obs = {s_obs}/{total_trials} ({obs_mean*100:.4f}%), "
        f"mu0 = {beta_ref['mean_true_far']*100:.4f}%, "
        f"MC p_two_sided = {p_two_sided:.5f}, "
        f"DECISION = {decision}, "
        f"VALIDATED = {payload['status_counts']['validated']}/{K_SEEDS}.",
        flush=True,
    )


if __name__ == "__main__":
    main()
