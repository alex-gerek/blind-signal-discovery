"""Template-free baseline anomaly detector for Blind Signal Discovery (BSD).

Methodology (classical signal processing and statistics, no templates, no ML):
1. Robust time-domain DC removal (median subtraction) and outlier diagnostics.
2. STFT power spectrogram using a Hann (or configured) window without zero-pad
   boundary artifacts. Interior frequency bins (excluding DC and Nyquist) are
   retained so every bin follows a common Exponential(1) null distribution under
   locally stationary Gaussian noise.
3. Robust 2D background power estimation:
   - Per-frequency baseline PSD estimated by temporal median divided by ln(2),
     then smoothed across frequency with a median filter so stationary narrowband
     signals are not absorbed into the background estimate.
   - Temporal variance envelope estimated by frequency median of the whitened
     spectrogram divided by ln(2), then smoothed across time with a median filter
     so short broadband pulses are not absorbed into the temporal background.
4. Normalized excess-power map Z[f, t] = P[f, t] / (S_bg(f) * T_bg(t)).
5. 8-connected-component clustering of bins exceeding a seed threshold, scored
   by integrated excess power sum_{(f,t) in C} (Z[f, t] - 1).
6. Comparison against a frozen, independently calibrated experiment-level threshold.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy import ndimage, signal as sp_signal, stats

from bsd.calibration import CalibrationResult
from bsd.generator import Observation

_LN2 = float(np.log(2.0))
_EPS = 1e-18


@dataclass(frozen=True)
class DetectorConfig:
    """Configuration for the template-free excess-power detector."""

    nperseg: int = 128
    noverlap: int = 96
    window: str = "hann"
    freq_median_kernel: int = 11
    time_median_kernel: int = 15
    bin_seed_threshold: float = 4.0
    min_cluster_bins: int = 2

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CandidateEvent:
    """A single anomalous time-frequency region extracted by the detector."""

    candidate_id: int
    t_start: float
    t_end: float
    f_min: float
    f_max: float
    peak_time: float
    peak_freq: float
    anomaly_score: float
    peak_excess_power: float
    n_bins: int
    empirical_p_value: float | None
    p_value_is_upper_bound: bool
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class TimeFrequencyMaps:
    """Intermediate STFT and excess-power arrays for visualization/diagnostics."""

    freqs: np.ndarray
    times: np.ndarray
    power_spectrogram: np.ndarray
    normalized_excess_map: np.ndarray
    background_psd: np.ndarray
    background_time_envelope: np.ndarray


@dataclass(frozen=True)
class DetectionResult:
    """Complete output of the detector on a single Observation."""

    sample_id: str
    detected: bool
    max_anomaly_score: float
    threshold_used: float | None
    candidates: list[CandidateEvent]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    tf_maps: TimeFrequencyMaps | None = None

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "detected": self.detected,
            "max_anomaly_score": self.max_anomaly_score,
            "threshold_used": self.threshold_used,
            "n_candidates": len(self.candidates),
            "candidates": [c.to_dict() for c in self.candidates],
            "diagnostics": dict(self.diagnostics),
        }


class TemplateFreeDetector:
    """Interpretable STFT excess-power anomaly detector without signal templates."""

    def __init__(self, config: DetectorConfig | None = None) -> None:
        self.config = config or DetectorConfig()
        if self.config.freq_median_kernel < 3 or self.config.freq_median_kernel % 2 == 0:
            raise ValueError("freq_median_kernel must be an odd integer >= 3.")
        if self.config.time_median_kernel < 3 or self.config.time_median_kernel % 2 == 0:
            raise ValueError("time_median_kernel must be an odd integer >= 3.")
        if self.config.noverlap >= self.config.nperseg:
            raise ValueError("noverlap must be strictly less than nperseg.")

    def analyze(
        self,
        observation: Observation,
        calibration: CalibrationResult | None = None,
        return_maps: bool = False,
    ) -> DetectionResult:
        """Run blind anomaly detection on a single Observation.

        Strict input isolation: only `Observation` instances are accepted.
        """
        if not isinstance(observation, Observation):
            raise TypeError(
                f"Detector strictly requires an Observation instance, got {type(observation).__name__}."
            )

        x_raw = np.asarray(observation.signal, dtype=np.float64)
        if x_raw.ndim != 1 or x_raw.size < self.config.nperseg:
            raise ValueError(
                f"Observation signal must be 1D with at least {self.config.nperseg} samples."
            )
        if not np.all(np.isfinite(x_raw)):
            raise ValueError(f"Observation '{observation.sample_id}' contains non-finite (NaN/inf) samples.")
        if not np.isfinite(observation.sample_rate) or observation.sample_rate <= 0.0:
            raise ValueError(f"Observation '{observation.sample_id}' has invalid sample_rate={observation.sample_rate}.")

        # 1. Robust time-domain preprocessing & diagnostics
        input_median = float(np.median(x_raw))
        x = x_raw - input_median
        mad = float(np.median(np.abs(x)))
        mad_sigma = max(1.4826 * mad, _EPS)
        input_kurtosis = float(stats.kurtosis(x, fisher=True))
        input_crest = float(np.max(np.abs(x)) / mad_sigma)

        # 2. STFT power spectrogram without zero-padding boundary artifacts
        freqs_all, times, stft_mat = sp_signal.stft(
            x,
            fs=observation.sample_rate,
            window=self.config.window,
            nperseg=self.config.nperseg,
            noverlap=self.config.noverlap,
            boundary=None,
            padded=False,
        )
        # Drop DC (0 Hz) and Nyquist bins (real-valued chi2_1 rather than complex Exp(1))
        freqs = freqs_all[1:-1]
        power = np.abs(stft_mat[1:-1, :]) ** 2
        n_freq, n_time = power.shape

        # 3. Robust 2D background power estimation
        # 3a. Spectral baseline: median over time, then median filter across frequency
        raw_psd = np.median(power, axis=1) / _LN2
        freq_k = min(self.config.freq_median_kernel, n_freq if n_freq % 2 == 1 else n_freq - 1)
        bg_psd = ndimage.median_filter(raw_psd, size=max(3, freq_k), mode="nearest")
        bg_psd = np.maximum(bg_psd, _EPS)

        # 3b. Spectral whitening
        whitened = power / bg_psd[:, None]

        # 3c. Slowly varying temporal noise floor: median over frequency, smoothed in time
        raw_time_env = np.median(whitened, axis=0) / _LN2
        time_k = min(self.config.time_median_kernel, n_time if n_time % 2 == 1 else n_time - 1)
        bg_time = ndimage.median_filter(raw_time_env, size=max(3, time_k), mode="nearest")
        bg_time = np.maximum(bg_time, _EPS)

        # 4. Normalized excess-power map Z[f, t]
        excess_map = whitened / bg_time[None, :]

        # 5. Connected-component clustering above bin seed threshold
        seed_mask = excess_map > self.config.bin_seed_threshold
        structure = np.ones((3, 3), dtype=int)  # 8-connectivity
        labeled, n_raw_clusters = ndimage.label(seed_mask, structure=structure)

        dt = float(times[1] - times[0]) if n_time > 1 else observation.duration
        df = float(freqs[1] - freqs[0]) if n_freq > 1 else observation.sample_rate / self.config.nperseg

        raw_candidates: list[CandidateEvent] = []
        max_score = 0.0

        for label_idx in range(1, n_raw_clusters + 1):
            f_idx, t_idx = np.where(labeled == label_idx)
            n_bins = int(f_idx.size)
            if n_bins < self.config.min_cluster_bins:
                continue

            cluster_z = excess_map[f_idx, t_idx]
            score = float(np.sum(cluster_z - 1.0))
            if score > max_score:
                max_score = score

            peak_pos = int(np.argmax(cluster_z))
            peak_z = float(cluster_z[peak_pos])
            peak_f = float(freqs[f_idx[peak_pos]])
            peak_t = float(times[t_idx[peak_pos]])

            t_min_idx = int(np.min(t_idx))
            t_max_idx = int(np.max(t_idx))
            f_min_idx = int(np.min(f_idx))
            f_max_idx = int(np.max(f_idx))

            t_start = max(0.0, float(times[t_min_idx]) - 0.5 * dt)
            t_end = min(observation.duration, float(times[t_max_idx]) + 0.5 * dt)
            f_min = max(0.0, float(freqs[f_min_idx]) - 0.5 * df)
            f_max = min(observation.sample_rate / 2.0, float(freqs[f_max_idx]) + 0.5 * df)

            boundary_adjacent = bool(
                t_min_idx == 0
                or t_max_idx == (n_time - 1)
                or f_min_idx == 0
                or f_max_idx == (n_freq - 1)
            )

            p_val: float | None = None
            p_is_upper = False
            if calibration is not None:
                p_val, p_is_upper = calibration.empirical_p_value(score)

            cand = CandidateEvent(
                candidate_id=len(raw_candidates) + 1,
                t_start=t_start,
                t_end=t_end,
                f_min=f_min,
                f_max=f_max,
                peak_time=peak_t,
                peak_freq=peak_f,
                anomaly_score=score,
                peak_excess_power=peak_z,
                n_bins=n_bins,
                empirical_p_value=p_val,
                p_value_is_upper_bound=p_is_upper,
                diagnostics={
                    "time_span_s": float(t_end - t_start),
                    "bandwidth_hz": float(f_max - f_min),
                    "mean_excess_power": float(np.mean(cluster_z)),
                    "boundary_adjacent": boundary_adjacent,
                    "input_kurtosis": input_kurtosis,
                    "input_crest_factor": input_crest,
                },
            )
            raw_candidates.append(cand)

        # Sort candidates by descending anomaly score
        raw_candidates.sort(key=lambda c: c.anomaly_score, reverse=True)

        # 6. Apply calibrated experiment-level threshold if provided
        det_available = True
        if calibration is not None:
            if calibration.threshold is None or not calibration.target_fap_attainable:
                det_available = False
                threshold_used: float | None = None
                detected = False
                candidates_out: list[CandidateEvent] = []
            else:
                threshold_used = calibration.threshold
                # Strict '>' matches the finite-sample guarantee P(T_new > T_(k)) <= target_fap
                accepted = [c for c in raw_candidates if c.anomaly_score > calibration.threshold]
                detected = len(accepted) > 0
                candidates_out = accepted
        else:
            threshold_used = None
            detected = False
            candidates_out = raw_candidates

        diagnostics = {
            "n_freq_bins": n_freq,
            "n_time_frames": n_time,
            "total_tf_bins": n_freq * n_time,
            "n_seed_bins": int(np.sum(seed_mask)),
            "n_raw_clusters": int(n_raw_clusters),
            "n_qualified_clusters": len(raw_candidates),
            "input_median": input_median,
            "input_mad_sigma": mad_sigma,
            "input_kurtosis": input_kurtosis,
            "input_crest_factor": input_crest,
            "detection_available": det_available,
        }

        tf_maps = (
            TimeFrequencyMaps(
                freqs=freqs,
                times=times,
                power_spectrogram=power,
                normalized_excess_map=excess_map,
                background_psd=bg_psd,
                background_time_envelope=bg_time,
            )
            if return_maps
            else None
        )

        return DetectionResult(
            sample_id=observation.sample_id,
            detected=detected,
            max_anomaly_score=max_score,
            threshold_used=threshold_used,
            candidates=candidates_out,
            diagnostics=diagnostics,
            tf_maps=tf_maps,
        )
