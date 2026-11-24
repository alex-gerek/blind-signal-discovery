"""Statistical threshold calibration and confidence intervals for BSD.

Calibration is performed strictly on an independent noise-only dataset to
control the experiment-level false-alarm probability (FAP) while accounting
for multiple comparisons across all time-frequency bins and STFT correlations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy import stats

if TYPE_CHECKING:
    from bsd.detector import TemplateFreeDetector
    from bsd.generator import Observation


def wilson_ci(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Compute Wilson score confidence interval for a binomial proportion."""
    if n <= 0:
        return (0.0, 1.0)
    z = float(stats.norm.ppf(0.5 + confidence / 2.0))
    p = k / n
    denom = 1.0 + (z**2) / n
    center = (p + (z**2) / (2.0 * n)) / denom
    half_width = (z / denom) * np.sqrt((p * (1.0 - p) / n) + (z**2) / (4.0 * n**2))
    return (max(0.0, float(center - half_width)), min(1.0, float(center + half_width)))


def clopper_pearson_ci(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Compute exact Clopper-Pearson confidence interval for a binomial proportion."""
    if n <= 0:
        return (0.0, 1.0)
    alpha = 1.0 - confidence
    lower = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2.0, k, n - k + 1))
    upper = 1.0 if k == n else float(stats.beta.ppf(1.0 - alpha / 2.0, k + 1, n - k))
    return (lower, upper)


@dataclass(frozen=True)
class CalibrationResult:
    """Frozen threshold calibration artifact produced from noise-only data."""

    target_fap: float
    n_calibration_samples: int
    threshold: float
    sorted_null_scores: np.ndarray
    calibration_noise_models: tuple[str, ...]
    total_tf_bins_per_sample: int
    naive_per_bin_exp_threshold: float
    bonferroni_per_bin_exp_threshold: float
    empirical_calibration_fap: float
    calibration_seed: int
    quantile_summary: dict[str, float]
    limitations: list[str] = field(default_factory=list)

    def empirical_p_value(self, score: float) -> tuple[float, bool]:
        """Compute finite-sample conformal empirical p-value against calibration null.

        Returns ``(p_value, is_upper_bound)`` where ``is_upper_bound`` is True
        when ``score`` exceeds every observed calibration noise score, meaning
        the true tail probability is bounded only by ``1 / (N_cal + 1)``.
        """
        n = len(self.sorted_null_scores)
        if n == 0:
            return (1.0, False)
        idx = int(np.searchsorted(self.sorted_null_scores, score, side="left"))
        n_ge = n - idx
        p_val = (1.0 + n_ge) / (n + 1.0)
        is_upper_bound = n_ge == 0
        return (float(p_val), is_upper_bound)

    def to_summary_dict(self) -> dict[str, Any]:
        """Serialize calibration summary (excluding raw score vector)."""
        return {
            "target_fap": self.target_fap,
            "n_calibration_samples": self.n_calibration_samples,
            "threshold": self.threshold,
            "calibration_noise_models": list(self.calibration_noise_models),
            "total_tf_bins_per_sample": self.total_tf_bins_per_sample,
            "naive_per_bin_exp_threshold": self.naive_per_bin_exp_threshold,
            "bonferroni_per_bin_exp_threshold": self.bonferroni_per_bin_exp_threshold,
            "empirical_calibration_fap": self.empirical_calibration_fap,
            "calibration_seed": self.calibration_seed,
            "quantile_summary": dict(self.quantile_summary),
            "min_resolvable_p_value": 1.0 / (self.n_calibration_samples + 1.0),
            "limitations": list(self.limitations),
        }


def calibrate_detector(
    detector: TemplateFreeDetector,
    calibration_observations: list[Observation],
    target_fap: float = 0.01,
    calibration_noise_models: tuple[str, ...] = (),
    calibration_seed: int = 0,
) -> CalibrationResult:
    """Calibrate experiment-level detection threshold on noise-only observations.

    Uses the conservative order statistic k = ceil((1 - alpha) * (N + 1)) of
    the experiment-level maximum cluster score T_i = max_C S(C) across each
    noise-only sample. Under exchangeability between calibration and test noise,
    P(T_new > T_(k)) <= alpha automatically accounts for scanning all
    time-frequency bins and their STFT overlap correlations.
    """
    if not (0.0 < target_fap < 1.0):
        raise ValueError(f"target_fap must lie in (0, 1), got {target_fap}.")
    n_cal = len(calibration_observations)
    if n_cal < 20:
        raise ValueError(f"At least 20 calibration samples are required, got {n_cal}.")

    null_scores = np.empty(n_cal, dtype=np.float64)
    total_tf_bins = 0

    for idx, obs in enumerate(calibration_observations):
        res = detector.analyze(obs, calibration=None, return_maps=False)
        null_scores[idx] = res.max_anomaly_score
        if idx == 0:
            total_tf_bins = int(res.diagnostics["total_tf_bins"])

    sorted_scores = np.sort(null_scores)
    # 1-indexed order statistic k = ceil((1 - alpha) * (N + 1)), clamped to N
    k_order = min(n_cal, int(np.ceil((1.0 - target_fap) * (n_cal + 1))))
    threshold = float(sorted_scores[k_order - 1])

    # Fraction of calibration samples strictly exceeding the threshold
    emp_cal_fap = float(np.mean(sorted_scores > threshold))

    # Reference single-bin exponential thresholds to contrast with experiment-level FAP
    naive_per_bin = float(-np.log(target_fap))
    bonferroni_per_bin = float(-np.log(target_fap / max(1, total_tf_bins)))

    quantile_summary = {
        "p50": float(np.quantile(sorted_scores, 0.50)),
        "p90": float(np.quantile(sorted_scores, 0.90)),
        "p95": float(np.quantile(sorted_scores, 0.95)),
        "p99": float(np.quantile(sorted_scores, 0.99)),
        "max": float(sorted_scores[-1]),
    }

    min_p = 1.0 / (n_cal + 1.0)
    expected_tail_count = target_fap * n_cal
    limitations = [
        (
            f"Finite-sample tail resolution: with N_cal={n_cal} and target_fap={target_fap:.3f}, "
            f"only ~{expected_tail_count:.1f} calibration samples lie in the upper tail; "
            f"minimum resolvable empirical p-value is {min_p:.5f}."
        ),
        (
            "Exchangeability assumption: experiment-level FAP control holds only when test "
            "noise follows the same distribution as the calibration noise mixture. Noise-model "
            "shifts (e.g., impulsive interference or unseen colored spectra) can inflate FAR."
        ),
        (
            "Single-observation horizon: FAP is calibrated per individual observation window; "
            "scanning M independent observations without further correction yields an expected "
            "false-alarm count of M * target_fap."
        ),
    ]

    return CalibrationResult(
        target_fap=float(target_fap),
        n_calibration_samples=n_cal,
        threshold=threshold,
        sorted_null_scores=sorted_scores,
        calibration_noise_models=tuple(calibration_noise_models),
        total_tf_bins_per_sample=total_tf_bins,
        naive_per_bin_exp_threshold=naive_per_bin,
        bonferroni_per_bin_exp_threshold=bonferroni_per_bin,
        empirical_calibration_fap=emp_cal_fap,
        calibration_seed=int(calibration_seed),
        quantile_summary=quantile_summary,
        limitations=limitations,
    )
