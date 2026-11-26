"""BSD v0.3 Phase 2b — Dependence-Aware Progressive Calibration, Orthogonal
Glitch/Continuum Diagnostics, and Multiple-Testing-Controlled Broadband Channel.

Preserves frozen v0.1 (`detector.py`), v0.2 (`detector_v2.py`, `vetoes.py`,
`coincidence.py`), v0.3 Phase 1 (`regime_monitor.py`), and v0.3 Phase 2
(`adaptive_regime.py`) completely unmodified.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy import signal as sp_signal, stats

from bsd.adaptive_regime import (
    CRED_CALIBRATED,
    CRED_CALIBRATION_UNAVAILABLE,
    CRED_RECALIBRATING,
    CRED_TRANSITION,
    SIG_EXACT,
    SIG_UNCERTAIN,
    extract_continuum_features,
)
from bsd.calibration import conformal_order_rank
from bsd.detector_v2 import (
    MultiResolutionCalibration,
    MultiResolutionCandidate,
    MultiResolutionDetector,
)
from bsd.generator import Observation
from bsd.regime_monitor import DEFAULT_FREQ_BANDS_HZ, ReferenceBackgroundState
from bsd.vetoes import VetoAssessment, VetoConfig, apply_vetoes

SIG_EMPIRICAL_PROVISIONAL = "EMPIRICAL_PROVISIONAL"
SIG_EMPIRICAL_VALIDATED = "EMPIRICAL_VALIDATED"
SIG_EMPIRICAL_DEPENDENT = "EMPIRICAL_DEPENDENT_UNCERTAIN"
SIG_EMPIRICAL_UNAVAILABLE = "EMPIRICAL_UNAVAILABLE_FAP_UNATTAINABLE"

CONT_STABLE = "STABLE"
CONT_SHIFT = "SHIFT_DETECTED"
CONT_RECAL = "RECALIBRATING"

GLITCH_QUIET = "QUIET"
GLITCH_ELEVATED = "ELEVATED_GLITCH_RATE"


@dataclass(frozen=True)
class CalibrationUncertaintySummary:
    """Finite-sample and temporal-dependence uncertainty audit for a calibration ensemble."""

    n_thresh: int
    n_eff: float
    lag1_autocorr: float
    ljung_box_p: float
    p_min: float
    expected_far_iid: float
    far_ci90_iid: tuple[float, float]
    far_ci90_neff: tuple[float, float]
    dependence_detected: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def estimate_calibration_uncertainty(
    chronological_null_scores: np.ndarray,
    target_fap: float = 0.01,
) -> CalibrationUncertaintySummary:
    """Test linear lag-1..10 autocorrelation and compute nominal i.i.d. + N_eff Beta FAR CIs.

    Note: Ljung-Box and Bartlett N_eff screen for linear autocorrelation at lags 1..10
    only; passing them does NOT prove statistical independence, stationarity, or
    real-data exchangeability.
    """
    x = np.asarray(chronological_null_scores, dtype=np.float64)
    if x.size > 0 and not np.all(np.isfinite(x)):
        raise ValueError("chronological_null_scores must contain only finite values.")
    if not np.isfinite(target_fap) or not (0.0 < target_fap < 1.0):
        raise ValueError(f"target_fap must be finite and in (0, 1), got {target_fap}.")
    n = int(x.size)
    if n < 5:
        return CalibrationUncertaintySummary(
            n_thresh=n,
            n_eff=float(n),
            lag1_autocorr=0.0,
            ljung_box_p=1.0,
            p_min=1.0 / (n + 1.0),
            expected_far_iid=1.0 / (n + 1.0),
            far_ci90_iid=(0.0, 1.0),
            far_ci90_neff=(0.0, 1.0),
            dependence_detected=True,
        )

    xc = x - float(np.mean(x))
    denom = float(np.dot(xc, xc))
    max_lag = max(1, min(10, n // 5))
    rhos: list[float] = []
    if denom > 1e-18:
        for lag in range(1, max_lag + 1):
            r = float(np.dot(xc[:-lag], xc[lag:]) / denom)
            rhos.append(r)
    else:
        rhos = [0.0] * max_lag

    lag1 = rhos[0]
    # Ljung-Box Q statistic
    q_lb = float(n * (n + 2) * sum((r * r) / max(1, n - (k + 1)) for k, r in enumerate(rhos)))
    p_lb = float(1.0 - stats.chi2.cdf(q_lb, df=max_lag))

    # Bartlett-weighted effective sample size
    bartlett_sum = sum((1.0 - (k + 1) / (max_lag + 1.0)) * max(0.0, r) for k, r in enumerate(rhos))
    n_eff = float(np.clip(n / (1.0 + 2.0 * bartlett_sum), 1.0, float(n)))

    # Order statistic k for target_fap under i.i.d. and under n_eff
    k_iid = min(n, max(1, int(np.ceil((1.0 - target_fap) * (n + 1)))))
    a_iid, b_iid = float(n - k_iid + 1), float(k_iid)
    e_far_iid = a_iid / (a_iid + b_iid)
    ci_iid = (
        float(stats.beta.ppf(0.05, a_iid, b_iid)),
        float(stats.beta.ppf(0.95, a_iid, b_iid)),
    )

    # Dependence-aware Beta CI using effective sample size n_eff
    n_eff_int = max(2.0, n_eff)
    k_eff = min(n_eff_int, max(1.0, np.ceil((1.0 - target_fap) * (n_eff_int + 1.0))))
    a_eff = max(0.5, n_eff_int - k_eff + 1.0)
    b_eff = max(0.5, k_eff)
    ci_neff = (
        float(stats.beta.ppf(0.05, a_eff, b_eff)),
        float(stats.beta.ppf(0.95, a_eff, b_eff)),
    )

    return CalibrationUncertaintySummary(
        n_thresh=n,
        n_eff=n_eff,
        lag1_autocorr=lag1,
        ljung_box_p=p_lb,
        p_min=1.0 / (n + 1.0),
        expected_far_iid=e_far_iid,
        far_ci90_iid=ci_iid,
        far_ci90_neff=ci_neff,
        dependence_detected=bool(p_lb < 0.05 or n_eff < 0.75 * n),
    )


@dataclass(frozen=True)
class BroadbandChannelCalibration:
    """Frozen-reference broadband excess-power calibration (matched FWER)."""

    band_medians: np.ndarray
    band_mad_scales: np.ndarray
    threshold: float | None
    alpha_channel: float
    sorted_null_scores: np.ndarray
    uncertainty: CalibrationUncertaintySummary
    target_fap_attainable: bool = True
    calibration_state: str = "AVAILABLE"

    def empirical_p_value(self, score: float) -> tuple[float, bool]:
        if not np.isfinite(score):
            return (1.0, False)
        n = int(self.sorted_null_scores.size)
        if n == 0:
            return (1.0, False)
        idx = int(np.searchsorted(self.sorted_null_scores, score, side="left"))
        n_ge = n - idx
        return (float((1.0 + n_ge) / (n + 1.0)), n_ge == 0)


@dataclass(frozen=True)
class JointDetectorCalibration:
    """Multiple-testing-controlled joint calibration for STFT + Frozen-Reference Broadband."""

    target_fwer: float
    alpha_per_family_sidak: float
    stft_calibration: MultiResolutionCalibration
    broadband_calibration: BroadbandChannelCalibration
    stft_uncertainty: CalibrationUncertaintySummary
    empirical_combined_cal_far: float


class FrozenReferenceBroadbandChannel:
    """External-reference broadband excess-power channel (complements frozen STFT)."""

    def __init__(
        self,
        freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ,
        stft_nperseg: int = 512,
    ) -> None:
        self.freq_bands_hz = freq_bands_hz
        self.stft_nperseg = stft_nperseg

    def _band_vec(self, obs: Observation) -> np.ndarray:
        x = np.asarray(obs.signal, dtype=np.float64)
        if x.size == 0 or not np.all(np.isfinite(x)):
            raise ValueError(f"Observation '{obs.sample_id}' must be non-empty and finite.")
        feat = extract_continuum_features(obs, self.freq_bands_hz, self.stft_nperseg)
        return feat[2:]  # 4 log continuum bands

    def score(self, obs: Observation, band_medians: np.ndarray, band_mad_scales: np.ndarray) -> float:
        """Positive broadband continuum excess chi-squared Z-sum relative to frozen reference R_0."""
        v = self._band_vec(obs)
        z = np.maximum(0.0, (v - band_medians) / np.maximum(band_mad_scales, 1e-9))
        return float(np.sum(z * z))

    def calibrate_with_split(
        self,
        norm_observations: list[Observation],
        thresh_observations: list[Observation],
        alpha_channel: float,
        allow_unattainable_fap: bool = False,
    ) -> BroadbandChannelCalibration:
        if not norm_observations or not thresh_observations:
            raise ValueError("norm_observations and thresh_observations must be non-empty.")
        if not np.isfinite(alpha_channel) or not (0.0 < alpha_channel < 1.0):
            raise ValueError(f"alpha_channel must be finite and in (0, 1), got {alpha_channel}.")

        n = len(thresh_observations)
        k_order, fap_attainable, min_req = conformal_order_rank(n, alpha_channel)
        if not fap_attainable and not allow_unattainable_fap:
            raise ValueError(
                f"Unattainable alpha_channel={alpha_channel:.5g} for n_thresh={n}: "
                f"p_min=1/({n}+1)={1.0 / (n + 1):.5f} > alpha_channel (need >= {min_req} Split-B samples)."
            )

        mat_norm = np.vstack([self._band_vec(o) for o in norm_observations])
        med = np.median(mat_norm, axis=0)
        mad = 1.4826 * np.median(np.abs(mat_norm - med[None, :]), axis=0)
        std_fb = np.std(mat_norm, axis=0)
        scales = np.where(mad > 1e-9, mad, np.where(std_fb > 1e-9, std_fb, 1.0))

        null_scores = np.array([self.score(o, med, scales) for o in thresh_observations], dtype=np.float64)
        unc = estimate_calibration_uncertainty(null_scores, target_fap=alpha_channel)
        sorted_s = np.sort(null_scores)
        if fap_attainable and k_order is not None:
            thresh: float | None = float(sorted_s[k_order - 1])
            cal_state = "AVAILABLE"
        else:
            thresh = None
            cal_state = "UNAVAILABLE_FAP_UNATTAINABLE"
        return BroadbandChannelCalibration(
            band_medians=med,
            band_mad_scales=scales,
            threshold=thresh,
            alpha_channel=alpha_channel,
            sorted_null_scores=sorted_s,
            uncertainty=unc,
            target_fap_attainable=fap_attainable,
            calibration_state=cal_state,
        )


def calibrate_joint_stft_and_broadband(
    detector: MultiResolutionDetector,
    broadband_channel: FrozenReferenceBroadbandChannel,
    calibration_observations: list[Observation],
    target_fwer: float = 0.01,
    norm_fraction: float = 0.30,
    calibration_seed: int = 0,
    allow_unattainable_fap: bool = False,
) -> JointDetectorCalibration:
    """Calibrate STFT bank + frozen-reference broadband channel with Sidak FWER budget.

    Note: The Sidak per-channel budget `1 - sqrt(1 - target_fwer)` is exact under
    channel independence (conservative under positive orthant dependence) and is
    NOT a distribution-free guarantee under arbitrary dependence; Bonferroni
    (`target_fwer / 2`) provides a distribution-free bound.
    """
    n_total = len(calibration_observations)
    if n_total < 60:
        raise ValueError(f"Need >= 60 calibration observations for joint calibration, got {n_total}.")
    if not np.isfinite(target_fwer) or not (0.0 < target_fwer < 1.0):
        raise ValueError(f"target_fwer must be finite and in (0, 1), got {target_fwer}.")
    alpha_each = float(1.0 - np.sqrt(1.0 - target_fwer))

    # Chronological split preserves temporal order for dependence testing on threshold split!
    n_norm = max(20, int(round(n_total * norm_fraction)))
    norm_obs = calibration_observations[:n_norm]
    thresh_obs = calibration_observations[n_norm:]

    stft_cal = detector.calibrate(
        calibration_observations,
        target_fap=alpha_each,
        norm_fraction=norm_fraction,
        calibration_seed=calibration_seed,
        allow_unattainable_fap=allow_unattainable_fap,
    )
    # Chronological STFT null scores on threshold split for temporal dependence audit
    stft_chrono = np.array([detector.analyze(o, stft_cal).max_global_score for o in thresh_obs], dtype=np.float64)
    stft_unc = estimate_calibration_uncertainty(stft_chrono, target_fap=alpha_each)

    bb_cal = broadband_channel.calibrate_with_split(
        norm_obs, thresh_obs, alpha_channel=alpha_each, allow_unattainable_fap=allow_unattainable_fap
    )
    bb_chrono = np.array(
        [broadband_channel.score(o, bb_cal.band_medians, bb_cal.band_mad_scales) for o in thresh_obs],
        dtype=np.float64,
    )
    stft_hit = (
        (stft_chrono > stft_cal.global_threshold)
        if stft_cal.global_threshold is not None
        else np.zeros_like(stft_chrono, dtype=bool)
    )
    bb_hit = (
        (bb_chrono > bb_cal.threshold)
        if bb_cal.threshold is not None
        else np.zeros_like(bb_chrono, dtype=bool)
    )
    comb_fp = np.mean(stft_hit | bb_hit)
    return JointDetectorCalibration(
        target_fwer=target_fwer,
        alpha_per_family_sidak=alpha_each,
        stft_calibration=stft_cal,
        broadband_calibration=bb_cal,
        stft_uncertainty=stft_unc,
        empirical_combined_cal_far=float(comb_fp),
    )


@dataclass(frozen=True)
class SeparateSearchModesCalibration:
    """Separately calibrated Transient Search, Persistent Search, and Sidak Union budgets.

    - `transient_stft_cal`: full `target_fap_transient` (default 1%) on STFT bank
      without penalizing transient sensitivity for the existence of a persistent channel.
    - `persistent_bb_cal`: full `target_fap_persistent` (default 1%) on frozen-reference
      broadband excess channel.
    - `union_stft_cal` & `union_bb_cal`: Sidak-corrected (`1 - sqrt(1 - target_fwer_union)`)
      thresholds used strictly when combining both modes into a single unified discovery claim.
    """

    target_fap_transient: float
    target_fap_persistent: float
    target_fwer_union: float
    alpha_union_sidak: float
    transient_stft_cal: MultiResolutionCalibration
    transient_stft_uncertainty: CalibrationUncertaintySummary
    persistent_bb_cal: BroadbandChannelCalibration
    union_stft_cal: MultiResolutionCalibration
    union_bb_cal: BroadbandChannelCalibration
    empirical_union_cal_far: float


def calibrate_separate_search_modes(
    detector: MultiResolutionDetector,
    broadband_channel: FrozenReferenceBroadbandChannel,
    calibration_observations: list[Observation],
    target_fap_transient: float = 0.01,
    target_fap_persistent: float = 0.01,
    target_fwer_union: float = 0.01,
    norm_fraction: float = 0.30,
    calibration_seed: int = 0,
) -> SeparateSearchModesCalibration:
    """Calibrate Transient Search, Persistent Search, and Sidak Combined Union budgets."""
    n_total = len(calibration_observations)
    if n_total < 60:
        raise ValueError(f"Need >= 60 calibration observations, got {n_total}.")
    alpha_union = float(1.0 - np.sqrt(1.0 - target_fwer_union))

    n_norm = max(20, int(round(n_total * norm_fraction)))
    norm_obs = calibration_observations[:n_norm]
    thresh_obs = calibration_observations[n_norm:]

    # 1. Dedicated Transient Search STFT calibration (full target_fap_transient)
    tr_stft_cal = detector.calibrate(
        calibration_observations,
        target_fap=target_fap_transient,
        norm_fraction=norm_fraction,
        calibration_seed=calibration_seed,
    )
    stft_chrono = np.array([detector.analyze(o, tr_stft_cal).max_global_score for o in thresh_obs], dtype=np.float64)
    tr_unc = estimate_calibration_uncertainty(stft_chrono, target_fap=target_fap_transient)

    # 2. Dedicated Persistent Search broadband calibration (full target_fap_persistent)
    pers_bb_cal = broadband_channel.calibrate_with_split(norm_obs, thresh_obs, alpha_channel=target_fap_persistent)

    # 3. Explicit Sidak-corrected Union calibration (assumes channel independence / POD;
    #    abstain with +inf threshold if Split-B sample size cannot attain alpha_union)
    union_stft_cal = detector.calibrate(
        calibration_observations,
        target_fap=alpha_union,
        norm_fraction=norm_fraction,
        calibration_seed=calibration_seed,
        allow_unattainable_fap=True,
    )
    union_bb_cal = broadband_channel.calibrate_with_split(
        norm_obs, thresh_obs, alpha_channel=alpha_union, allow_unattainable_fap=True
    )
    bb_chrono = np.array(
        [broadband_channel.score(o, union_bb_cal.band_medians, union_bb_cal.band_mad_scales) for o in thresh_obs],
        dtype=np.float64,
    )
    stft_u_hit = (
        (stft_chrono > union_stft_cal.global_threshold)
        if union_stft_cal.global_threshold is not None
        else np.zeros_like(stft_chrono, dtype=bool)
    )
    bb_u_hit = (
        (bb_chrono > union_bb_cal.threshold)
        if union_bb_cal.threshold is not None
        else np.zeros_like(bb_chrono, dtype=bool)
    )
    union_fp = float(np.mean(stft_u_hit | bb_u_hit))

    return SeparateSearchModesCalibration(
        target_fap_transient=target_fap_transient,
        target_fap_persistent=target_fap_persistent,
        target_fwer_union=target_fwer_union,
        alpha_union_sidak=alpha_union,
        transient_stft_cal=tr_stft_cal,
        transient_stft_uncertainty=tr_unc,
        persistent_bb_cal=pers_bb_cal,
        union_stft_cal=union_stft_cal,
        union_bb_cal=union_bb_cal,
        empirical_union_cal_far=union_fp,
    )


@dataclass(frozen=True)
class GlitchReferenceState:
    """Reference baseline for short-timescale broadband impulsive glitch-storm indicator."""

    crest_median: float
    crest_mad: float
    single_window_z_thresh: float = 4.5
    elevated_window_z_thresh: float = 2.8


def is_glitch_free_reference_block(
    observations: list[Observation],
    initial_glitch_ref: GlitchReferenceState,
    max_outlier_fraction: float = 0.05,
    max_mad_inflation_ratio: float = 2.5,
) -> bool:
    """Verify a candidate new-regime block is free of impulsive glitch-storm contamination."""
    if len(observations) < 20:
        return False
    vals = np.array([extract_broadband_impulsive_crest(o) for o in observations], dtype=np.float64)
    med = float(np.median(vals))
    mad = float(1.4826 * np.median(np.abs(vals - med)))
    scale = max(mad, initial_glitch_ref.crest_mad, 1e-6)
    outlier_frac = float(np.mean((vals - med) / scale >= 4.0))
    mad_ratio = mad / max(initial_glitch_ref.crest_mad, 1e-6)
    return bool(outlier_frac <= max_outlier_fraction and mad_ratio <= max_mad_inflation_ratio)


def extract_broadband_impulsive_crest(observation: Observation, nperseg: int = 32) -> float:
    """Short-timescale (31 ms) frequency-median frame power crest ratio.

    Taking the median across all frequency bins in each 31-ms frame makes this
    statistic invariant to narrowband sinusoids while strongly responsive to
    broadband impulsive glitches.
    """
    x = np.asarray(observation.signal, dtype=np.float64)
    x = x - float(np.median(x))
    _f, _t, zxx = sp_signal.stft(
        x,
        fs=observation.sample_rate,
        window="hann",
        nperseg=min(nperseg, x.size),
        noverlap=min(nperseg, x.size) * 3 // 4,
        boundary=None,
        padded=False,
    )
    power = np.abs(zxx[1:-1, :]) ** 2
    frame_med = np.median(power, axis=0)
    base = max(float(np.median(frame_med)), 1e-30)
    peak = max(float(np.max(frame_med)), 1e-30)
    return float(np.log(peak / base))


def fit_glitch_reference(reference_observations: list[Observation]) -> GlitchReferenceState:
    vals = np.array([extract_broadband_impulsive_crest(o) for o in reference_observations], dtype=np.float64)
    med = float(np.median(vals))
    mad = float(1.4826 * np.median(np.abs(vals - med)))
    scale = mad if mad > 1e-9 else max(float(np.std(vals)), 1e-6)
    return GlitchReferenceState(crest_median=med, crest_mad=scale)


@dataclass(frozen=True)
class AdaptiveRegimeV2bConfig:
    """Configuration for Phase 2b AdaptiveRegimeManagerV2b."""

    mode: str = "causal_trailing"  # "causal_trailing" or "offline_symmetric"
    policy: str = "transient"  # "transient" or "persistent"
    local_window_count: int = 12
    ks_alpha: float = 0.001
    drift_z_threshold: float = 4.0
    cusum_k: float = 1.0
    cusum_h: float = 22.0
    shift_confirm_windows: int = 2
    stability_confirm_windows: int = 10
    provisional_recal_windows: int = 60
    mature_recal_windows: int = 200
    max_recal_buffer_windows: int = 2000
    internal_stability_z_max: float = 3.0
    glitch_rate_threshold: float = 0.25  # >= 3 of 12 windows elevated -> GLITCH_STORM
    target_fwer: float = 0.01
    min_exact_neff: float = 160.0
    min_validated_neff: float = 120.0
    freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ
    stft_nperseg: int = 512

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DualOutputWindowV2b:
    """Non-destructive Phase 2b dual-output record with orthogonal diagnostics & raw preservation."""

    index: int
    sample_id: str
    # Stream A: Complete raw measurements & separate nomination decisions (never deleted)
    raw_stft_candidates: tuple[MultiResolutionCandidate, ...]
    surviving_stft_candidates: tuple[MultiResolutionCandidate, ...]
    stft_max_score: float
    stft_threshold: float | None
    stft_nominated_raw: bool
    stft_nominated_post_veto: bool
    broadband_excess_score: float
    broadband_threshold: float | None
    broadband_nominated: bool
    joint_nominated_post_veto: bool
    veto_assessment: VetoAssessment | None
    # Stream B: Orthogonal diagnostics & dependence-aware calibration credibility
    regime_id: int
    credibility_status: str
    continuum_status: str
    glitch_storm_status: str
    warning_causes: tuple[str, ...]
    significance_status: str
    calibration_uncertainty: CalibrationUncertaintySummary
    mode: str
    policy: str
    future_data_used: bool
    min_latency_windows: int
    glitch_crest_z: float
    max_feature_z_drift: float
    diagnostics: dict[str, Any] = field(default_factory=dict)
    detection_available: bool = True

    @property
    def is_clean_calibrated(self) -> bool:
        """True only when continuum is calibrated, no glitch storm is active, AND detection is available."""
        return (
            self.credibility_status == CRED_CALIBRATED
            and self.glitch_storm_status == GLITCH_QUIET
            and self.detection_available
        )


class AdaptiveRegimeManagerV2b:
    """Phase 2b dual-output manager with staged causal calibration & orthogonal diagnostics."""

    def __init__(self, config: AdaptiveRegimeV2bConfig | None = None) -> None:
        self.config = config or AdaptiveRegimeV2bConfig()
        if self.config.mode not in ("causal_trailing", "offline_symmetric"):
            raise ValueError(f"Invalid mode: {self.config.mode!r}")
        if self.config.policy not in ("transient", "persistent"):
            raise ValueError(f"Invalid policy: {self.config.policy!r}")
        for attr in ("max_recal_buffer_windows", "provisional_recal_windows", "mature_recal_windows"):
            val = getattr(self.config, attr)
            if isinstance(val, bool) or not isinstance(val, int) or val < 1:
                raise ValueError(f"{attr} must be a positive integer >= 1, got {val!r}.")

    def fit_continuum_reference(self, reference_observations: list[Observation]) -> ReferenceBackgroundState:
        cfg = self.config
        mat = np.vstack(
            [extract_continuum_features(o, cfg.freq_bands_hz, cfg.stft_nperseg) for o in reference_observations]
        )
        med = np.median(mat, axis=0)
        mad = 1.4826 * np.median(np.abs(mat - med[None, :]), axis=0)
        std_fb = np.std(mat, axis=0)
        scales = np.where(mad > 1e-9, mad, np.where(std_fb > 1e-9, std_fb, 1.0))
        bnames = tuple(f"log_cont_{int(a)}_{int(b)}hz" for a, b in cfg.freq_bands_hz)
        return ReferenceBackgroundState(
            feature_names=("log_continuum_total", "log_continuum_q90_burst", *bnames),
            ref_matrix=mat,
            ref_medians=med,
            ref_mad_scales=scales,
            n_reference_windows=len(reference_observations),
        )

    def _is_internally_stable(self, block: np.ndarray, scale_ref: np.ndarray) -> bool:
        n = block.shape[0]
        if n < 6:
            return False
        mid = n // 2
        m1 = np.median(block[:mid, :], axis=0)
        m2 = np.median(block[mid:, :], axis=0)
        bmod = np.median(block, axis=0)
        bmad = 1.4826 * np.median(np.abs(block - bmod[None, :]), axis=0)
        eff = np.maximum(np.maximum(scale_ref, bmad), 1e-9)
        return bool(float(np.max(np.abs(m2 - m1) / eff)) <= self.config.internal_stability_z_max)

    def _block_matches_reference(self, block: np.ndarray, ref: ReferenceBackgroundState) -> bool:
        nf = block.shape[1]
        loc_med = np.median(block, axis=0)
        z_shift = np.abs(loc_med - ref.ref_medians) / ref.ref_mad_scales
        if float(np.max(z_shift)) >= self.config.drift_z_threshold:
            return False
        p_vals = [float(stats.ks_2samp(block[:, k], ref.ref_matrix[:, k]).pvalue) for k in range(nf)]
        return bool(min(1.0, nf * min(p_vals)) >= self.config.ks_alpha)

    def _classify_significance(
        self,
        regime_id: int,
        unc: CalibrationUncertaintySummary,
        is_mature: bool,
    ) -> str:
        cfg = self.config
        if regime_id == 0 and unc.n_eff >= cfg.min_exact_neff and not unc.dependence_detected:
            return SIG_EXACT
        if is_mature:
            if unc.n_eff >= cfg.min_validated_neff and unc.far_ci90_neff[1] <= 0.025:
                return SIG_EMPIRICAL_VALIDATED
            return SIG_EMPIRICAL_DEPENDENT
        return SIG_EMPIRICAL_PROVISIONAL

    def process_sequence(
        self,
        observations: list[Observation],
        detector: MultiResolutionDetector,
        broadband_channel: FrozenReferenceBroadbandChannel,
        joint_calibration: JointDetectorCalibration | SeparateSearchModesCalibration,
        initial_continuum_ref: ReferenceBackgroundState,
        glitch_ref: GlitchReferenceState,
        veto_config: VetoConfig | None = None,
    ) -> list[DualOutputWindowV2b]:
        if not observations:
            return []
        cfg = self.config
        n_obs = len(observations)
        w = cfg.local_window_count
        is_sym = cfg.mode == "offline_symmetric"
        half_w = w // 2
        min_lat_wins = (half_w + 1) if is_sym else 1

        feats_mat = np.vstack(
            [extract_continuum_features(o, cfg.freq_bands_hz, cfg.stft_nperseg) for o in observations]
        )
        n_features = feats_mat.shape[1]
        # Raw per-window crest feature (purely local to window i; normalized causally by active_glitch_ref)
        raw_crest_arr = np.array([extract_broadband_impulsive_crest(o) for o in observations], dtype=np.float64)

        # Frozen-reference broadband channel NEVER updates its reference R_0!
        is_separate = isinstance(joint_calibration, SeparateSearchModesCalibration)
        if is_separate:
            bb_cal = joint_calibration.persistent_bb_cal
            active_stft_cal = joint_calibration.transient_stft_cal
            active_unc = joint_calibration.transient_stft_uncertainty
            stft_recal_fap = joint_calibration.target_fap_transient
            union_stft_thresh = joint_calibration.union_stft_cal.global_threshold
            union_bb_thresh = joint_calibration.union_bb_cal.threshold
            union_sidak_fap = joint_calibration.alpha_union_sidak
        else:
            bb_cal = joint_calibration.broadband_calibration
            active_stft_cal = joint_calibration.stft_calibration
            active_unc = joint_calibration.stft_uncertainty
            stft_recal_fap = joint_calibration.alpha_per_family_sidak
            union_stft_thresh = active_stft_cal.global_threshold
            union_bb_thresh = bb_cal.threshold
            union_sidak_fap = joint_calibration.alpha_per_family_sidak

        active_ref = initial_continuum_ref
        active_glitch_ref = glitch_ref
        active_is_mature = True
        glitch_ref_refresh_count = 0
        glitch_ref_blocked_count = 0

        # Dynamically compute minimum total windows needed for BOTH dedicated STFT FAP
        # and (in separate-modes) the stricter Sidak Union STFT FAP at norm_fraction=0.25.
        def _req_total_for_fap(target_fap: float) -> int:
            _, _, min_nth = conformal_order_rank(1, target_fap)
            req_nth = max(40, min_nth)
            ntot = max(60, int(np.ceil(req_nth / 0.75)))
            while ntot - max(20, int(round(ntot * 0.25))) < req_nth:
                ntot += 1
            while ntot > 60 and (ntot - 1) - max(20, int(round((ntot - 1) * 0.25))) >= req_nth:
                ntot -= 1
            return ntot

        required_total_for_stft = _req_total_for_fap(stft_recal_fap)
        required_total_for_union = _req_total_for_fap(union_sidak_fap) if is_separate else required_total_for_stft
        effective_mature_windows = max(cfg.mature_recal_windows, required_total_for_stft, required_total_for_union)
        resource_limit_blocked_stft = (
            required_total_for_stft > cfg.max_recal_buffer_windows
            or cfg.provisional_recal_windows > cfg.max_recal_buffer_windows
        )
        resource_limit_blocked_union = is_separate and (
            required_total_for_union > cfg.max_recal_buffer_windows
            or cfg.provisional_recal_windows > cfg.max_recal_buffer_windows
        )
        resource_limit_blocked = resource_limit_blocked_stft or resource_limit_blocked_union

        regime_id = 0
        state = CRED_CALIBRATED
        cusum_s = 0.0
        consec_shift = 0
        regime_start_idx = 0
        regime_accum_indices: list[int] = []
        transition_indices: list[int] = []
        recal_indices: list[int] = []
        outputs: list[DualOutputWindowV2b] = []

        for i, obs in enumerate(observations):
            # Snapshot exact calibration state actually used to evaluate window i
            # (any recalibration triggered at window i takes effect strictly at window i + 1)
            used_stft_cal = active_stft_cal
            used_union_stft_thresh = union_stft_thresh
            used_unc = active_unc
            used_regime_id = regime_id
            used_is_mature = active_is_mature
            used_stft_avail = bool(used_stft_cal.target_fap_attainable and used_stft_cal.global_threshold is not None)
            used_bb_avail = bool(bb_cal.target_fap_attainable and bb_cal.threshold is not None)
            used_union_stft_avail = (
                bool(used_union_stft_thresh is not None and np.isfinite(used_union_stft_thresh))
                if is_separate
                else used_stft_avail
            )
            used_union_bb_avail = (
                bool(union_bb_thresh is not None and np.isfinite(union_bb_thresh))
                if is_separate
                else used_bb_avail
            )
            if used_union_stft_avail and used_union_bb_avail:
                used_union_status = "FULL_UNION_AVAILABLE"
            elif used_union_bb_avail:
                used_union_status = "PARTIAL_BROADBAND_ONLY"
            elif used_union_stft_avail:
                used_union_status = "PARTIAL_STFT_ONLY"
            else:
                used_union_status = "UNION_UNAVAILABLE"
            used_det_avail = used_stft_avail if cfg.policy == "transient" else used_bb_avail
            cal_updated_for_next = False

            # --- Stream A: ALWAYS compute both raw scores & nominations using used_stft_cal ---
            det_res = detector.analyze(obs, used_stft_cal)
            raw_cands = tuple(det_res.candidates)
            if veto_config is not None:
                v_assess: VetoAssessment | None = apply_vetoes(obs, raw_cands, veto_config)
                surv_cands = tuple(v_assess.surviving_candidates)
            else:
                v_assess = None
                surv_cands = raw_cands

            bb_score = broadband_channel.score(obs, bb_cal.band_medians, bb_cal.band_mad_scales)
            bb_nom = bool(bb_cal.threshold is not None and bb_score > bb_cal.threshold)
            stft_raw_nom = len(raw_cands) > 0
            stft_post_nom = len(surv_cands) > 0
            # Explicit Sidak-controlled union nomination when combining both search modes
            if is_separate:
                stft_union_nom = bool(
                    stft_post_nom
                    and used_union_stft_thresh is not None
                    and det_res.max_global_score > used_union_stft_thresh
                )
                bb_union_nom = bool(union_bb_thresh is not None and bb_score > union_bb_thresh)
                joint_nom = bool(stft_union_nom or bb_union_nom)
            else:
                joint_nom = bool(stft_post_nom or bb_nom)

            # --- Stream B Channel 2: Orthogonal Impulsive Glitch-Storm Monitor (causal against active_glitch_ref) ---
            if is_sym:
                g_lo = max(regime_start_idx, i - half_w)
                g_hi = min(n_obs, i + half_w + 1)
            else:
                g_lo = max(regime_start_idx, i - w + 1)
                g_hi = i + 1
            gz_block = (raw_crest_arr[g_lo:g_hi] - active_glitch_ref.crest_median) / active_glitch_ref.crest_mad
            gz_i = float((raw_crest_arr[i] - active_glitch_ref.crest_median) / active_glitch_ref.crest_mad)
            g_elev_rate = float(np.mean(gz_block >= active_glitch_ref.elevated_window_z_thresh)) if gz_block.size else 0.0
            glitch_storm_active = bool(
                gz_i >= active_glitch_ref.single_window_z_thresh
                or (g_hi - g_lo >= 6 and g_elev_rate >= cfg.glitch_rate_threshold)
            )
            glitch_status = GLITCH_ELEVATED if glitch_storm_active else GLITCH_QUIET

            # --- Stream B Channel 1: Spectral Continuum Monitor + Staged Causal Recalibration ---
            z_cur = np.abs(feats_mat[i] - active_ref.ref_medians) / active_ref.ref_mad_scales
            abs_z_i = float(np.max(z_cur))
            max_z_drift = abs_z_i

            if state == CRED_CALIBRATED:
                cusum_s = max(0.0, cusum_s + (abs_z_i - 2.0 - cfg.cusum_k))
                if is_sym:
                    lo = max(regime_start_idx, i - half_w)
                    hi = min(n_obs, i + half_w + 1)
                    enough = (hi - lo) >= max(4, half_w + 1)
                else:
                    lo = max(regime_start_idx, i - w + 1)
                    hi = i + 1
                    enough = (hi - lo) >= w

                shifted = False
                if enough:
                    block = feats_mat[lo:hi, :]
                    loc_med = np.median(block, axis=0)
                    z_shift = np.abs(loc_med - active_ref.ref_medians) / active_ref.ref_mad_scales
                    max_z_drift = float(np.max(z_shift))
                    p_vals = [
                        float(stats.ks_2samp(block[:, k], active_ref.ref_matrix[:, k]).pvalue)
                        for k in range(n_features)
                    ]
                    ks_p = min(1.0, n_features * min(p_vals))
                    shifted = bool(
                        (ks_p < cfg.ks_alpha)
                        or (max_z_drift >= cfg.drift_z_threshold)
                        or (cusum_s >= cfg.cusum_h)
                    )

                consec_shift = (consec_shift + 1) if shifted else 0
                if shifted and consec_shift >= cfg.shift_confirm_windows:
                    if cfg.policy == "persistent":
                        cred_out = CRED_TRANSITION
                        cont_status = CONT_SHIFT
                        sig_out = SIG_UNCERTAIN
                    else:
                        state = CRED_TRANSITION
                        transition_indices = [i]
                        recal_indices = []
                        regime_accum_indices = []
                        cusum_s = 0.0
                        consec_shift = 0
                        cred_out = CRED_TRANSITION
                        cont_status = CONT_SHIFT
                        sig_out = SIG_UNCERTAIN
                else:
                    # Classify window i strictly using the calibration that evaluated window i
                    if used_stft_cal.target_fap_attainable and used_stft_cal.global_threshold is not None:
                        cred_out = CRED_CALIBRATED
                        cont_status = CONT_STABLE
                        sig_out = self._classify_significance(used_regime_id, used_unc, used_is_mature)
                    else:
                        cred_out = CRED_CALIBRATION_UNAVAILABLE
                        cont_status = CONT_RECAL
                        sig_out = SIG_EMPIRICAL_UNAVAILABLE

                    # Progressive Stage 1 -> Stage 2 causal upgrade (strictly past windows <= i; takes effect at i + 1)
                    if regime_id >= 1 and not active_is_mature:
                        regime_accum_indices.append(i)
                        if len(regime_accum_indices) > cfg.max_recal_buffer_windows:
                            regime_accum_indices = regime_accum_indices[-cfg.max_recal_buffer_windows :]
                        n_acc = len(regime_accum_indices)
                        n_norm_acc = max(20, int(round(n_acc * 0.25)))
                        n_thresh_acc = n_acc - n_norm_acc
                        _, acc_stft_ok, _ = conformal_order_rank(max(1, n_thresh_acc), stft_recal_fap)
                        _, acc_union_ok, _ = (
                            conformal_order_rank(max(1, n_thresh_acc), union_sidak_fap)
                            if is_separate
                            else (None, acc_stft_ok, 0)
                        )
                        need_stft_rec = (
                            (not active_stft_cal.target_fap_attainable)
                            and acc_stft_ok
                            and n_thresh_acc >= 40
                            and not resource_limit_blocked_stft
                        )
                        need_union_rec = (
                            is_separate
                            and (union_stft_thresh is None)
                            and acc_union_ok
                            and n_thresh_acc >= 40
                            and not resource_limit_blocked_union
                        )
                        need_mature = (
                            n_acc >= effective_mature_windows
                            and acc_stft_ok
                            and (not is_separate or acc_union_ok)
                            and n_thresh_acc >= 40
                            and not resource_limit_blocked
                        )
                        if need_stft_rec or need_union_rec or need_mature:
                            slice_len = min(effective_mature_windows, cfg.max_recal_buffer_windows)
                            m_obs = [observations[idx] for idx in regime_accum_indices[-slice_len:]]
                            active_ref = self.fit_continuum_reference(m_obs)
                            active_stft_cal = detector.calibrate(
                                m_obs,
                                target_fap=stft_recal_fap,
                                norm_fraction=0.25,
                                calibration_seed=regime_id * 10 + 2,
                                allow_unattainable_fap=True,
                            )
                            if is_separate:
                                union_stft_thresh = detector.calibrate(
                                    m_obs,
                                    target_fap=union_sidak_fap,
                                    norm_fraction=0.25,
                                    calibration_seed=regime_id * 10 + 2,
                                    allow_unattainable_fap=True,
                                ).global_threshold
                            else:
                                union_stft_thresh = active_stft_cal.global_threshold
                            m_chrono = np.array(
                                [
                                    detector.analyze(o, active_stft_cal).max_global_score
                                    for o in m_obs[active_stft_cal.n_norm_samples :]
                                ],
                                dtype=np.float64,
                            )
                            active_unc = estimate_calibration_uncertainty(m_chrono, target_fap=stft_recal_fap)
                            if is_glitch_free_reference_block(m_obs, glitch_ref):
                                active_glitch_ref = fit_glitch_reference(m_obs)
                                glitch_ref_refresh_count += 1
                            else:
                                glitch_ref_blocked_count += 1
                            cal_updated_for_next = True
                            if (
                                n_acc >= effective_mature_windows
                                and active_stft_cal.target_fap_attainable
                                and active_stft_cal.global_threshold is not None
                                and (not is_separate or union_stft_thresh is not None)
                                and not resource_limit_blocked
                            ):
                                active_is_mature = True

            elif state == CRED_TRANSITION:
                transition_indices.append(i)
                cred_out = CRED_TRANSITION
                cont_status = CONT_SHIFT
                sig_out = SIG_UNCERTAIN
                if len(transition_indices) >= cfg.stability_confirm_windows:
                    recent_idx = transition_indices[-cfg.stability_confirm_windows :]
                    block_new = feats_mat[recent_idx, :]
                    if self._block_matches_reference(block_new, active_ref):
                        state = CRED_CALIBRATED
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                    elif self._is_internally_stable(block_new, active_ref.ref_mad_scales):
                        state = CRED_RECALIBRATING
                        recal_indices = list(recent_idx)

            else:  # CRED_RECALIBRATING (Stage 0 -> Stage 1 provisional calibration)
                recal_indices.append(i)
                if len(recal_indices) > cfg.max_recal_buffer_windows:
                    recal_indices = recal_indices[-cfg.max_recal_buffer_windows :]
                if cfg.provisional_recal_windows > cfg.max_recal_buffer_windows:
                    cred_out = CRED_CALIBRATION_UNAVAILABLE
                    cont_status = CONT_RECAL
                    sig_out = SIG_EMPIRICAL_UNAVAILABLE
                    used_det_avail = False
                    used_stft_avail = False
                    used_union_stft_avail = False
                    used_union_status = "PARTIAL_BROADBAND_ONLY" if used_union_bb_avail else "UNION_UNAVAILABLE"
                    used_is_mature = False
                else:
                    cred_out = CRED_RECALIBRATING
                    cont_status = CONT_RECAL
                    sig_out = SIG_UNCERTAIN
                if len(recal_indices) >= cfg.provisional_recal_windows:
                    recal_block = feats_mat[recal_indices, :]
                    if self._block_matches_reference(recal_block, active_ref):
                        state = CRED_CALIBRATED
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                        recal_indices = []
                    elif self._is_internally_stable(recal_block, active_ref.ref_mad_scales):
                        p_obs = [observations[idx] for idx in recal_indices]
                        active_ref = self.fit_continuum_reference(p_obs)
                        active_stft_cal = detector.calibrate(
                            p_obs,
                            target_fap=stft_recal_fap,
                            norm_fraction=0.25,
                            calibration_seed=regime_id * 10 + 1,
                            allow_unattainable_fap=True,
                        )
                        if is_separate:
                            union_stft_thresh = detector.calibrate(
                                p_obs,
                                target_fap=union_sidak_fap,
                                norm_fraction=0.25,
                                calibration_seed=regime_id * 10 + 1,
                                allow_unattainable_fap=True,
                            ).global_threshold
                        else:
                            union_stft_thresh = active_stft_cal.global_threshold
                        p_chrono = np.array(
                            [
                                detector.analyze(o, active_stft_cal).max_global_score
                                for o in p_obs[active_stft_cal.n_norm_samples :]
                            ],
                            dtype=np.float64,
                        )
                        active_unc = estimate_calibration_uncertainty(p_chrono, target_fap=stft_recal_fap)
                        # Correction 2: synchronize glitch_ref ONLY if p_obs is uncontaminated
                        if is_glitch_free_reference_block(p_obs, glitch_ref):
                            active_glitch_ref = fit_glitch_reference(p_obs)
                            glitch_ref_refresh_count += 1
                        else:
                            glitch_ref_blocked_count += 1
                        regime_id += 1
                        active_is_mature = False
                        regime_accum_indices = list(recal_indices[-cfg.max_recal_buffer_windows :])
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                        recal_indices = []
                        state = CRED_CALIBRATED
                        cal_updated_for_next = True
                    else:
                        recal_indices = recal_indices[-cfg.stability_confirm_windows :]

            causes: list[str] = []
            if cont_status != CONT_STABLE:
                causes.append("CONTINUUM_DRIFT")
            if glitch_status != GLITCH_QUIET:
                causes.append("GLITCH_STORM")
            if bb_nom:
                causes.append("SINGLE_CHANNEL_BROADBAND_AMBIGUOUS")

            if not used_det_avail:
                causes.append("CALIBRATION_UNAVAILABLE")
            if is_separate and not used_union_stft_avail:
                causes.append("UNION_STFT_UNAVAILABLE")

            outputs.append(
                DualOutputWindowV2b(
                    index=i,
                    sample_id=obs.sample_id,
                    raw_stft_candidates=raw_cands,
                    surviving_stft_candidates=surv_cands,
                    stft_max_score=det_res.max_global_score,
                    stft_threshold=used_stft_cal.global_threshold,
                    stft_nominated_raw=stft_raw_nom,
                    stft_nominated_post_veto=stft_post_nom,
                    broadband_excess_score=bb_score,
                    broadband_threshold=bb_cal.threshold,
                    broadband_nominated=bb_nom,
                    joint_nominated_post_veto=joint_nom,
                    veto_assessment=v_assess,
                    regime_id=used_regime_id,
                    credibility_status=cred_out,
                    continuum_status=cont_status,
                    glitch_storm_status=glitch_status,
                    warning_causes=tuple(causes),
                    significance_status=sig_out,
                    calibration_uncertainty=used_unc,
                    mode=cfg.mode,
                    policy=cfg.policy,
                    future_data_used=is_sym,
                    min_latency_windows=min_lat_wins,
                    glitch_crest_z=gz_i,
                    max_feature_z_drift=max_z_drift,
                    diagnostics={
                        "is_mature": used_is_mature,
                        "glitch_elev_rate": g_elev_rate,
                        "glitch_ref_refresh_count": glitch_ref_refresh_count,
                        "glitch_ref_blocked_count": glitch_ref_blocked_count,
                        "separate_modes": is_separate,
                        "detection_available": used_det_avail,
                        "dedicated_stft_available": used_stft_avail,
                        "dedicated_broadband_available": used_bb_avail,
                        "union_stft_available": used_union_stft_avail,
                        "union_broadband_available": used_union_bb_avail,
                        "union_stft_threshold": used_union_stft_thresh,
                        "union_mode_status": used_union_status,
                        "stft_calibration_state": used_stft_cal.calibration_state,
                        "stft_n_norm_samples": used_stft_cal.n_norm_samples,
                        "stft_n_threshold_samples": used_stft_cal.n_threshold_samples,
                        "effective_mature_recal_windows": effective_mature_windows,
                        "max_recal_buffer_windows": cfg.max_recal_buffer_windows,
                        "resource_limit_blocked_stft": resource_limit_blocked_stft,
                        "resource_limit_blocked_union": resource_limit_blocked_union,
                        "resource_limit_blocked_recovery": resource_limit_blocked,
                        "calibration_updated_for_next_window": cal_updated_for_next,
                    },
                    detection_available=used_det_avail,
                )
            )

        return outputs
