"""BSD v0.2 Phase 1 — Multi-Resolution Template-Free Anomaly Detector.

Implements a 3-scale STFT bank (31.25 ms, 125 ms, and 500 ms at 1024 Hz)
wrapping the untouched BSD v0.1 `TemplateFreeDetector` per scale, with
split-sample per-scale score standardization and unified experiment-level
1% FAP calibration across all resolutions simultaneously.

Phase 1 intentionally excludes instrumental glitch/line vetoes and
cross-detector coincidence logic (scheduled for subsequent v0.2 phases).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from bsd.calibration import conformal_order_rank
from bsd.detector import CandidateEvent, DetectorConfig, TemplateFreeDetector
from bsd.generator import Observation


@dataclass(frozen=True)
class ScaleSpec:
    """Configuration for one STFT resolution in the multi-resolution bank."""

    name: str
    detector_config: DetectorConfig


DEFAULT_V02_SCALES: tuple[ScaleSpec, ...] = (
    ScaleSpec(
        name="short_31ms",
        detector_config=DetectorConfig(
            nperseg=32,
            noverlap=28,
            window="hann",
            freq_median_kernel=5,
            time_median_kernel=31,
            bin_seed_threshold=4.0,
            min_cluster_bins=2,
        ),
    ),
    ScaleSpec(
        name="medium_125ms",
        detector_config=DetectorConfig(
            nperseg=128,
            noverlap=96,
            window="hann",
            freq_median_kernel=11,
            time_median_kernel=15,
            bin_seed_threshold=4.0,
            min_cluster_bins=2,
        ),
    ),
    ScaleSpec(
        name="long_500ms",
        detector_config=DetectorConfig(
            nperseg=512,
            noverlap=384,
            window="hann",
            freq_median_kernel=21,
            time_median_kernel=7,
            bin_seed_threshold=4.0,
            min_cluster_bins=2,
        ),
    ),
)


@dataclass(frozen=True)
class ScaleNormalizationStats:
    """Robust per-scale null location/scale estimated on calibration split A."""

    scale_name: str
    median_null_score: float
    q95_null_score: float
    spread_scale: float
    n_norm_samples: int

    def standardize(self, raw_score: float) -> float:
        return float((raw_score - self.median_null_score) / self.spread_scale)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MultiResolutionCalibration:
    """Frozen joint calibration across all STFT resolutions."""

    target_fap: float
    n_norm_samples: int
    n_threshold_samples: int
    global_threshold: float | None
    scale_stats: dict[str, ScaleNormalizationStats]
    sorted_global_null_scores: np.ndarray
    empirical_calibration_fap: float
    calibration_seed: int
    quantile_summary: dict[str, float]
    scale_win_fractions_in_null: dict[str, float]
    target_fap_attainable: bool = True
    calibration_state: str = "AVAILABLE"

    def empirical_p_value(self, standardized_score: float) -> tuple[float, bool]:
        """Conformal finite-sample empirical p-value against global null."""
        if not np.isfinite(standardized_score):
            return (1.0, False)
        n = len(self.sorted_global_null_scores)
        if n == 0:
            return (1.0, False)
        idx = int(np.searchsorted(self.sorted_global_null_scores, standardized_score, side="left"))
        n_ge = n - idx
        p_val = (1.0 + n_ge) / (n + 1.0)
        return (float(p_val), n_ge == 0)

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "target_fap": self.target_fap,
            "target_fap_attainable": self.target_fap_attainable,
            "calibration_state": self.calibration_state,
            "n_norm_samples": self.n_norm_samples,
            "n_threshold_samples": self.n_threshold_samples,
            "n_total_calibration_samples": self.n_norm_samples + self.n_threshold_samples,
            "global_threshold": self.global_threshold,
            "empirical_calibration_fap": self.empirical_calibration_fap,
            "min_resolvable_p_value": 1.0 / (self.n_threshold_samples + 1.0),
            "calibration_seed": self.calibration_seed,
            "quantile_summary": dict(self.quantile_summary),
            "scale_stats": {k: v.to_dict() for k, v in self.scale_stats.items()},
            "scale_win_fractions_in_null": dict(self.scale_win_fractions_in_null),
        }


@dataclass(frozen=True)
class MultiResolutionCandidate:
    """A single candidate event from one scale, with raw and standardized scores."""

    candidate_id: int
    scale_name: str
    t_start: float
    t_end: float
    f_min: float
    f_max: float
    peak_time: float
    peak_freq: float
    raw_anomaly_score: float
    standardized_score: float
    peak_excess_power: float
    n_bins: int
    empirical_p_value: float | None
    p_value_is_upper_bound: bool
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MultiResolutionResult:
    """Output of MultiResolutionDetector on one Observation."""

    sample_id: str
    detected: bool
    max_global_score: float
    winning_scale: str | None
    threshold_used: float | None | None
    per_scale_raw_max: dict[str, float]
    per_scale_std_max: dict[str, float]
    candidates: list[MultiResolutionCandidate]
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "detected": self.detected,
            "max_global_score": self.max_global_score,
            "winning_scale": self.winning_scale,
            "threshold_used": self.threshold_used,
            "per_scale_raw_max": dict(self.per_scale_raw_max),
            "per_scale_std_max": dict(self.per_scale_std_max),
            "n_candidates": len(self.candidates),
            "candidates": [c.to_dict() for c in self.candidates],
            "diagnostics": dict(self.diagnostics),
        }


class MultiResolutionDetector:
    """3-scale STFT template-free anomaly detector (BSD v0.2 Phase 1).

    Wraps the frozen BSD v0.1 `TemplateFreeDetector` across multiple STFT
    resolutions without modifying the baseline detector.
    """

    def __init__(self, scales: tuple[ScaleSpec, ...] = DEFAULT_V02_SCALES) -> None:
        if not scales:
            raise ValueError("At least one ScaleSpec is required.")
        names = [s.name for s in scales]
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate scale names: {names}")
        self.scales = scales
        self._detectors: dict[str, TemplateFreeDetector] = {
            s.name: TemplateFreeDetector(s.detector_config) for s in scales
        }

    def _raw_analyze_all_scales(
        self, observation: Observation
    ) -> dict[str, tuple[float, list[CandidateEvent], dict[str, Any]]]:
        if not isinstance(observation, Observation):
            raise TypeError(
                f"MultiResolutionDetector strictly requires an Observation, got {type(observation).__name__}."
            )
        out: dict[str, tuple[float, list[CandidateEvent], dict[str, Any]]] = {}
        for spec in self.scales:
            res = self._detectors[spec.name].analyze(observation, calibration=None, return_maps=False)
            out[spec.name] = (res.max_anomaly_score, res.candidates, res.diagnostics)
        return out

    def calibrate(
        self,
        calibration_observations: list[Observation],
        target_fap: float = 0.01,
        norm_fraction: float = 0.30,
        calibration_seed: int = 0,
        allow_unattainable_fap: bool = False,
    ) -> MultiResolutionCalibration:
        """Split-sample joint calibration across all STFT resolutions.

        Split A (`norm_fraction`) estimates robust per-scale null location
        (median) and upper-tail scale (`Q_0.95 - median`). Split B (remaining
        independent samples) calibrates the joint global-maximum threshold
        via the conservative order statistic `k = ceil((1 - alpha) * (N_B + 1))`.

        When `k > N_B` (i.e. `1 / (N_B + 1) > target_fap`), the requested FAP
        is mathematically unattainable:
          - By default (`allow_unattainable_fap=False`), raises ValueError.
          - If `allow_unattainable_fap=True`, sets `global_threshold = None`
            and `calibration_state = "UNAVAILABLE_FAP_UNATTAINABLE"` so
            detection explicitly abstains (JSON-safe `null`).
        """
        if not np.isfinite(target_fap) or not (0.0 < target_fap < 1.0):
            raise ValueError(f"target_fap must be finite and lie in (0, 1), got {target_fap}.")
        if not np.isfinite(norm_fraction) or not (0.0 < norm_fraction < 1.0):
            raise ValueError(f"norm_fraction must be finite and lie in (0, 1), got {norm_fraction}.")
        n_total = len(calibration_observations)
        if n_total < 60:
            raise ValueError(f"At least 60 calibration observations are required, got {n_total}.")

        n_norm = max(20, int(round(n_total * norm_fraction)))
        n_thresh = n_total - n_norm
        if n_thresh < 40:
            raise ValueError(f"Threshold calibration split too small: {n_thresh}.")

        k_order, fap_attainable, min_required_n_thresh = conformal_order_rank(n_thresh, target_fap)
        if not fap_attainable and not allow_unattainable_fap:
            raise ValueError(
                f"Unattainable target_fap={target_fap:.5g} for n_thresh={n_thresh} Split-B samples "
                f"(n_total={n_total}, norm_fraction={norm_fraction}): minimum resolvable "
                f"p_min=1/({n_thresh}+1)={1.0 / (n_thresh + 1):.5f} > target_fap. "
                f"Need at least {min_required_n_thresh} Split-B samples or target_fap >= {1.0 / (n_thresh + 1):.5f}."
            )

        norm_obs = calibration_observations[:n_norm]
        thresh_obs = calibration_observations[n_norm:]

        # Split A: per-scale null normalization statistics
        raw_norm_scores: dict[str, list[float]] = {s.name: [] for s in self.scales}
        for obs in norm_obs:
            per_scale = self._raw_analyze_all_scales(obs)
            for name, (max_s, _cands, _diag) in per_scale.items():
                raw_norm_scores[name].append(max_s)

        scale_stats: dict[str, ScaleNormalizationStats] = {}
        for s in self.scales:
            arr = np.asarray(raw_norm_scores[s.name], dtype=np.float64)
            med = float(np.median(arr))
            q95 = float(np.quantile(arr, 0.95))
            spread = max(q95 - med, 1.0)
            scale_stats[s.name] = ScaleNormalizationStats(
                scale_name=s.name,
                median_null_score=med,
                q95_null_score=q95,
                spread_scale=spread,
                n_norm_samples=n_norm,
            )

        # Split B: joint global-maximum threshold calibration
        global_scores = np.empty(n_thresh, dtype=np.float64)
        win_counts: dict[str, int] = {s.name: 0 for s in self.scales}

        for idx, obs in enumerate(thresh_obs):
            per_scale = self._raw_analyze_all_scales(obs)
            best_std = -float("inf")
            best_name = self.scales[0].name
            for s in self.scales:
                raw_max = per_scale[s.name][0]
                std_val = scale_stats[s.name].standardize(raw_max)
                if std_val > best_std:
                    best_std = std_val
                    best_name = s.name
            global_scores[idx] = best_std
            win_counts[best_name] += 1

        sorted_global = np.sort(global_scores)
        if fap_attainable and k_order is not None:
            global_threshold: float | None = float(sorted_global[k_order - 1])
            emp_cal_fap = float(np.mean(sorted_global > global_threshold))
            cal_state = "AVAILABLE"
        else:
            global_threshold = None
            emp_cal_fap = 0.0
            cal_state = "UNAVAILABLE_FAP_UNATTAINABLE"

        q_summary = {
            "p50": float(np.quantile(sorted_global, 0.50)),
            "p90": float(np.quantile(sorted_global, 0.90)),
            "p95": float(np.quantile(sorted_global, 0.95)),
            "p99": float(np.quantile(sorted_global, 0.99)),
            "max": float(sorted_global[-1]),
        }
        win_fracs = {k: float(v / n_thresh) for k, v in win_counts.items()}

        return MultiResolutionCalibration(
            target_fap=float(target_fap),
            n_norm_samples=n_norm,
            n_threshold_samples=n_thresh,
            global_threshold=global_threshold,
            scale_stats=scale_stats,
            sorted_global_null_scores=sorted_global,
            empirical_calibration_fap=emp_cal_fap,
            calibration_seed=int(calibration_seed),
            quantile_summary=q_summary,
            scale_win_fractions_in_null=win_fracs,
            target_fap_attainable=fap_attainable,
            calibration_state=cal_state,
        )

    def analyze(
        self,
        observation: Observation,
        calibration: MultiResolutionCalibration,
    ) -> MultiResolutionResult:
        """Run blind multi-resolution detection using a frozen joint calibration."""
        if not isinstance(calibration, MultiResolutionCalibration):
            raise TypeError("analyze requires a frozen MultiResolutionCalibration.")

        per_scale = self._raw_analyze_all_scales(observation)
        per_scale_raw_max: dict[str, float] = {}
        per_scale_std_max: dict[str, float] = {}
        all_accepted: list[MultiResolutionCandidate] = []

        max_global = -float("inf")
        winning_scale: str | None = None
        total_tf_bins = 0
        det_available = bool(calibration.target_fap_attainable and calibration.global_threshold is not None)

        for s in self.scales:
            raw_max, raw_cands, diag = per_scale[s.name]
            stats = calibration.scale_stats[s.name]
            std_max = stats.standardize(raw_max)
            per_scale_raw_max[s.name] = float(raw_max)
            per_scale_std_max[s.name] = float(std_max)
            total_tf_bins += int(diag["total_tf_bins"])

            if std_max > max_global:
                max_global = std_max
                winning_scale = s.name

            if det_available and calibration.global_threshold is not None:
                for c in raw_cands:
                    std_score = stats.standardize(c.anomaly_score)
                    if std_score > calibration.global_threshold:
                        p_val, p_ub = calibration.empirical_p_value(std_score)
                        all_accepted.append(
                            MultiResolutionCandidate(
                                candidate_id=len(all_accepted) + 1,
                                scale_name=s.name,
                                t_start=c.t_start,
                                t_end=c.t_end,
                                f_min=c.f_min,
                                f_max=c.f_max,
                                peak_time=c.peak_time,
                                peak_freq=c.peak_freq,
                                raw_anomaly_score=c.anomaly_score,
                                standardized_score=std_score,
                                peak_excess_power=c.peak_excess_power,
                                n_bins=c.n_bins,
                                empirical_p_value=p_val,
                                p_value_is_upper_bound=p_ub,
                                diagnostics=dict(c.diagnostics),
                            )
                        )

        all_accepted.sort(key=lambda c: c.standardized_score, reverse=True)
        detected = det_available and len(all_accepted) > 0

        return MultiResolutionResult(
            sample_id=observation.sample_id,
            detected=detected,
            max_global_score=float(max_global),
            winning_scale=winning_scale,
            threshold_used=calibration.global_threshold if det_available else None,
            per_scale_raw_max=per_scale_raw_max,
            per_scale_std_max=per_scale_std_max,
            candidates=all_accepted,
            diagnostics={
                "n_scales": len(self.scales),
                "total_tf_bins_across_scales": total_tf_bins,
                "detection_available": det_available,
                "calibration_state": calibration.calibration_state,
            },
        )
