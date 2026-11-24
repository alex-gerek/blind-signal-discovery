"""BSD v0.3 Phase 1 — Independent Background Regime Monitor.

Monitors statistical stationarity of background noise across sequential
time-series windows (`Observation`) without requiring signal templates and
without modifying or suppressing raw anomaly candidates.

Explicitly separates two operating modes:
  - ``causal_trailing``: strictly causal (uses only past and current windows,
    zero future-data dependence, latency = 1 observation window).
  - ``offline_symmetric``: non-causal retrospective mode (uses symmetric
    past+future local window around index `i`, future-data dependence =
    `local_window // 2` windows).

Regime-shift decisions are reported as transparent credibility metadata
(`STATIONARY_CREDIBLE`, `REGIME_SHIFT_DETECTED`, `WARMUP_INSUFFICIENT_HISTORY`)
and NEVER automatically classify shifts as instrumental noise or delete
detector candidates.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy import signal as sp_signal, stats

from bsd.generator import Observation

STATUS_STATIONARY = "STATIONARY_CREDIBLE"
STATUS_REGIME_SHIFT = "REGIME_SHIFT_DETECTED"
STATUS_WARMUP = "WARMUP_INSUFFICIENT_HISTORY"

DEFAULT_FREQ_BANDS_HZ: tuple[tuple[float, float], ...] = (
    (20.0, 80.0),
    (80.0, 160.0),
    (160.0, 280.0),
    (280.0, 450.0),
)


@dataclass(frozen=True)
class RegimeMonitorConfig:
    """Configuration for BackgroundRegimeMonitor."""

    mode: str = "causal_trailing"  # "causal_trailing" or "offline_symmetric"
    local_window_count: int = 12  # Number of sequential windows in local test block
    ks_alpha: float = 0.005  # Bonferroni-corrected KS significance level
    drift_z_threshold: float = 3.5  # Robust median-shift Z threshold
    cusum_k: float = 0.75  # Page-Hinkley / CUSUM slack parameter
    cusum_h: float = 12.0  # CUSUM alarm threshold on standardized feature shift
    freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ
    stft_nperseg: int = 128

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class WindowBackgroundFeatures:
    """Template-free robust background summary vector for one Observation."""

    sample_id: str
    mad_sigma: float
    tail_q99_ratio: float
    log_band_psd: tuple[float, ...]

    def as_vector(self) -> np.ndarray:
        return np.array(
            [np.log(max(self.mad_sigma, 1e-18)), self.tail_q99_ratio, *self.log_band_psd],
            dtype=np.float64,
        )


@dataclass(frozen=True)
class ReferenceBackgroundState:
    """Frozen reference background feature distribution fitted on calibration data."""

    feature_names: tuple[str, ...]
    ref_matrix: np.ndarray  # shape (N_ref, n_features)
    ref_medians: np.ndarray
    ref_mad_scales: np.ndarray
    n_reference_windows: int

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "n_reference_windows": self.n_reference_windows,
            "feature_names": list(self.feature_names),
            "ref_medians": [float(v) for v in self.ref_medians],
            "ref_mad_scales": [float(v) for v in self.ref_mad_scales],
        }


@dataclass(frozen=True)
class RegimeAssessment:
    """Per-window background stationarity and credibility report."""

    index: int
    sample_id: str
    mode: str
    future_data_used: bool
    min_latency_windows: int
    min_latency_seconds: float
    credibility_status: str
    regime_shifted: bool
    ks_min_p_bonferroni: float | None
    max_feature_z_drift: float
    cusum_stat: float
    change_point_flag: bool
    dominant_shift_feature: str | None
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def extract_window_features(
    observation: Observation,
    freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ,
    stft_nperseg: int = 128,
) -> WindowBackgroundFeatures:
    """Extract robust time- and frequency-domain noise features from one Observation."""
    if not isinstance(observation, Observation):
        raise TypeError(
            f"BackgroundRegimeMonitor strictly requires an Observation, got {type(observation).__name__}."
        )
    x = np.asarray(observation.signal, dtype=np.float64)
    x = x - float(np.median(x))
    abs_x = np.abs(x)
    mad_sigma = max(1.4826 * float(np.median(abs_x)), 1e-18)
    q99 = float(np.quantile(abs_x, 0.99))
    tail_ratio = q99 / mad_sigma

    freqs, _times, zxx = sp_signal.stft(
        x,
        fs=observation.sample_rate,
        window="hann",
        nperseg=min(stft_nperseg, x.size),
        noverlap=min(stft_nperseg, x.size) // 2,
        boundary=None,
        padded=False,
    )
    power = np.abs(zxx) ** 2
    med_psd = np.median(power, axis=1) / float(np.log(2.0))

    band_vals: list[float] = []
    for f_lo, f_hi in freq_bands_hz:
        mask = (freqs >= f_lo) & (freqs <= f_hi)
        if np.any(mask):
            val = float(np.median(med_psd[mask]))
        else:
            val = float(np.median(med_psd))
        band_vals.append(float(np.log(max(val, 1e-30))))

    return WindowBackgroundFeatures(
        sample_id=observation.sample_id,
        mad_sigma=mad_sigma,
        tail_q99_ratio=tail_ratio,
        log_band_psd=tuple(band_vals),
    )


class BackgroundRegimeMonitor:
    """Template-free background stationarity and change-point monitor."""

    def __init__(self, config: RegimeMonitorConfig | None = None) -> None:
        self.config = config or RegimeMonitorConfig()
        if self.config.mode not in ("causal_trailing", "offline_symmetric"):
            raise ValueError(
                f"Unsupported mode {self.config.mode!r}; expected 'causal_trailing' or 'offline_symmetric'."
            )
        if self.config.local_window_count < 4:
            raise ValueError("local_window_count must be >= 4.")

    def _feature_names(self) -> tuple[str, ...]:
        band_names = tuple(f"log_psd_{int(a)}_{int(b)}hz" for a, b in self.config.freq_bands_hz)
        return ("log_mad_sigma", "tail_q99_ratio", *band_names)

    def fit_reference(self, reference_observations: list[Observation]) -> ReferenceBackgroundState:
        """Fit frozen reference background feature distribution on calibration windows."""
        if len(reference_observations) < 20:
            raise ValueError(
                f"At least 20 reference background windows are required, got {len(reference_observations)}."
            )
        feats = [
            extract_window_features(obs, self.config.freq_bands_hz, self.config.stft_nperseg).as_vector()
            for obs in reference_observations
        ]
        mat = np.vstack(feats)
        medians = np.median(mat, axis=0)
        mad_scales = 1.4826 * np.median(np.abs(mat - medians[None, :]), axis=0)
        # Guard against zero dispersion on constant features
        std_fallback = np.std(mat, axis=0)
        scales = np.where(mad_scales > 1e-9, mad_scales, np.where(std_fallback > 1e-9, std_fallback, 1.0))

        return ReferenceBackgroundState(
            feature_names=self._feature_names(),
            ref_matrix=mat,
            ref_medians=medians,
            ref_mad_scales=scales,
            n_reference_windows=len(reference_observations),
        )

    def monitor_sequence(
        self,
        observations: list[Observation],
        reference: ReferenceBackgroundState,
    ) -> list[RegimeAssessment]:
        """Evaluate background regime stationarity across a sequential list of Observations.

        Strictly preserves causality in ``causal_trailing`` mode: assessment at
        index ``i`` depends ONLY on ``observations[0 : i + 1]`` and ``reference``.
        """
        if not observations:
            return []
        cfg = self.config
        w = cfg.local_window_count
        is_symmetric = cfg.mode == "offline_symmetric"
        half_w = w // 2
        n_obs = len(observations)
        win_dur = float(observations[0].duration)
        min_lat_wins = (half_w + 1) if is_symmetric else 1
        min_lat_s = min_lat_wins * win_dur

        feats_mat = np.vstack(
            [
                extract_window_features(obs, cfg.freq_bands_hz, cfg.stft_nperseg).as_vector()
                for obs in observations
            ]
        )
        n_features = feats_mat.shape[1]

        # Per-window standardized L-infinity shift relative to frozen reference
        z_mat = (feats_mat - reference.ref_medians[None, :]) / reference.ref_mad_scales[None, :]
        abs_z_max = np.max(np.abs(z_mat), axis=1)

        # Two-sided CUSUM (Page-Hinkley style) on standardized max-feature shift
        # Strictly causal by construction (depends only on 0..i)
        cusum_vals = np.zeros(n_obs, dtype=np.float64)
        s_pos = 0.0
        for i in range(n_obs):
            # Under stationary Gaussian-like reference, E[max_k |z_k|] ~ 1.6 for K=6;
            # subtract baseline offset + slack k so stationary drift decays to 0
            s_pos = max(0.0, s_pos + (abs_z_max[i] - 1.8 - cfg.cusum_k))
            cusum_vals[i] = s_pos

        assessments: list[RegimeAssessment] = []
        for i, obs in enumerate(observations):
            if is_symmetric:
                lo = max(0, i - half_w)
                hi = min(n_obs, i + half_w + 1)
                has_enough = (hi - lo) >= max(4, half_w + 1)
            else:
                lo = max(0, i - w + 1)
                hi = i + 1
                has_enough = (hi - lo) >= w

            if not has_enough:
                assessments.append(
                    RegimeAssessment(
                        index=i,
                        sample_id=obs.sample_id,
                        mode=cfg.mode,
                        future_data_used=is_symmetric,
                        min_latency_windows=min_lat_wins,
                        min_latency_seconds=min_lat_s,
                        credibility_status=STATUS_WARMUP,
                        regime_shifted=False,
                        ks_min_p_bonferroni=None,
                        max_feature_z_drift=float(abs_z_max[i]),
                        cusum_stat=float(cusum_vals[i]),
                        change_point_flag=False,
                        dominant_shift_feature=None,
                        diagnostics={"local_span": [lo, hi]},
                    )
                )
                continue

            local_block = feats_mat[lo:hi, :]
            local_medians = np.median(local_block, axis=0)
            local_z_shift = np.abs(local_medians - reference.ref_medians) / reference.ref_mad_scales
            dom_idx = int(np.argmax(local_z_shift))
            max_z_drift = float(local_z_shift[dom_idx])

            # Two-sample KS test per feature with Bonferroni correction
            p_vals = [
                float(stats.ks_2samp(local_block[:, k], reference.ref_matrix[:, k]).pvalue)
                for k in range(n_features)
            ]
            min_p_raw = min(p_vals)
            min_p_bonf = min(1.0, n_features * min_p_raw)
            ks_dom_idx = int(np.argmin(p_vals))

            ks_triggered = min_p_bonf < cfg.ks_alpha
            z_triggered = max_z_drift >= cfg.drift_z_threshold
            cp_flag = bool(cusum_vals[i] >= cfg.cusum_h)
            shifted = bool(ks_triggered or z_triggered or cp_flag)

            dom_name = reference.feature_names[ks_dom_idx if ks_triggered else dom_idx]
            status = STATUS_REGIME_SHIFT if shifted else STATUS_STATIONARY

            assessments.append(
                RegimeAssessment(
                    index=i,
                    sample_id=obs.sample_id,
                    mode=cfg.mode,
                    future_data_used=is_symmetric,
                    min_latency_windows=min_lat_wins,
                    min_latency_seconds=min_lat_s,
                    credibility_status=status,
                    regime_shifted=shifted,
                    ks_min_p_bonferroni=min_p_bonf,
                    max_feature_z_drift=max_z_drift,
                    cusum_stat=float(cusum_vals[i]),
                    change_point_flag=cp_flag,
                    dominant_shift_feature=dom_name,
                    diagnostics={
                        "local_span": [lo, hi],
                        "ks_triggered": ks_triggered,
                        "z_triggered": z_triggered,
                        "cusum_triggered": cp_flag,
                        "per_feature_z_shift": {
                            reference.feature_names[k]: round(float(local_z_shift[k]), 3)
                            for k in range(n_features)
                        },
                    },
                )
            )

        return assessments


def summarize_gating_impact(
    assessments: list[RegimeAssessment],
    detected_flags: list[bool],
    is_signal_flags: list[bool],
    window_duration_s: float = 4.0,
) -> dict[str, Any]:
    """Compute mandatory regime-gating metrics without false credit for excluded intervals.

    Excluded (gated-out) background windows are NEVER counted as true negatives.
    End-to-end detection efficiency is reported over ALL target signals both
    before and after gating.
    """
    n = len(assessments)
    if not (len(detected_flags) == n and len(is_signal_flags) == n):
        raise ValueError("assessments, detected_flags, and is_signal_flags must have equal length.")

    gated_mask = [a.regime_shifted for a in assessments]
    n_excluded = sum(gated_mask)
    total_livetime_hr = (n * window_duration_s) / 3600.0
    excluded_livetime_hr = (n_excluded * window_duration_s) / 3600.0
    duty_cycle = 1.0 - (n_excluded / n) if n > 0 else 0.0

    # Background windows (is_signal == False)
    bg_total = sum(1 for s in is_signal_flags if not s)
    bg_fp_before = sum(1 for d, s in zip(detected_flags, is_signal_flags) if (not s) and d)
    # Surviving (non-gated) background windows ONLY as denominator for post-gating FAR
    bg_surviving = sum(1 for g, s in zip(gated_mask, is_signal_flags) if (not s) and (not g))
    bg_fp_after = sum(
        1 for d, g, s in zip(detected_flags, gated_mask, is_signal_flags) if (not s) and (not g) and d
    )

    # Signal windows (is_signal == True)
    sig_total = sum(1 for s in is_signal_flags if s)
    sig_tp_before = sum(1 for d, s in zip(detected_flags, is_signal_flags) if s and d)
    sig_lost_to_gating = sum(
        1 for d, g, s in zip(detected_flags, gated_mask, is_signal_flags) if s and d and g
    )
    sig_tp_after = sum(
        1 for d, g, s in zip(detected_flags, gated_mask, is_signal_flags) if s and d and (not g)
    )

    return {
        "n_windows": n,
        "n_excluded_windows": n_excluded,
        "total_livetime_hours": total_livetime_hr,
        "excluded_livetime_hours": excluded_livetime_hr,
        "effective_duty_cycle": duty_cycle,
        "bg_n_total": bg_total,
        "bg_fp_before_gating": bg_fp_before,
        "bg_far_before_gating": (bg_fp_before / bg_total) if bg_total > 0 else None,
        "bg_n_surviving": bg_surviving,
        "bg_fp_after_gating": bg_fp_after,
        "bg_far_after_gating_on_surviving": (bg_fp_after / bg_surviving) if bg_surviving > 0 else None,
        "sig_n_total": sig_total,
        "sig_tp_before_gating": sig_tp_before,
        "sig_tp_after_gating": sig_tp_after,
        "signals_lost_to_gating": sig_lost_to_gating,
        "efficiency_before_gating": (sig_tp_before / sig_total) if sig_total > 0 else None,
        "efficiency_after_gating": (sig_tp_after / sig_total) if sig_total > 0 else None,
    }
