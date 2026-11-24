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
    CRED_RECALIBRATING,
    CRED_TRANSITION,
    SIG_EXACT,
    SIG_UNCERTAIN,
    extract_continuum_features,
)
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
    """Test temporal dependence and compute nominal i.i.d. + dependence-aware Beta FAR CIs."""
    x = np.asarray(chronological_null_scores, dtype=np.float64)
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
    threshold: float
    alpha_channel: float
    sorted_null_scores: np.ndarray
    uncertainty: CalibrationUncertaintySummary

    def empirical_p_value(self, score: float) -> tuple[float, bool]:
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
    ) -> BroadbandChannelCalibration:
        mat_norm = np.vstack([self._band_vec(o) for o in norm_observations])
        med = np.median(mat_norm, axis=0)
        mad = 1.4826 * np.median(np.abs(mat_norm - med[None, :]), axis=0)
        std_fb = np.std(mat_norm, axis=0)
        scales = np.where(mad > 1e-9, mad, np.where(std_fb > 1e-9, std_fb, 1.0))

        null_scores = np.array([self.score(o, med, scales) for o in thresh_observations], dtype=np.float64)
        unc = estimate_calibration_uncertainty(null_scores, target_fap=alpha_channel)
        sorted_s = np.sort(null_scores)
        n = sorted_s.size
        k = min(n, max(1, int(np.ceil((1.0 - alpha_channel) * (n + 1)))))
        thresh = float(sorted_s[k - 1])
        return BroadbandChannelCalibration(
            band_medians=med,
            band_mad_scales=scales,
            threshold=thresh,
            alpha_channel=alpha_channel,
            sorted_null_scores=sorted_s,
            uncertainty=unc,
        )


def calibrate_joint_stft_and_broadband(
    detector: MultiResolutionDetector,
    broadband_channel: FrozenReferenceBroadbandChannel,
    calibration_observations: list[Observation],
    target_fwer: float = 0.01,
    norm_fraction: float = 0.30,
    calibration_seed: int = 0,
) -> JointDetectorCalibration:
    """Calibrate STFT bank + frozen-reference broadband channel with Sidak FWER control."""
    n_total = len(calibration_observations)
    if n_total < 60:
        raise ValueError(f"Need >= 60 calibration observations for joint calibration, got {n_total}.")
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
    )
    # Chronological STFT null scores on threshold split for temporal dependence audit
    stft_chrono = np.array([detector.analyze(o, stft_cal).max_global_score for o in thresh_obs], dtype=np.float64)
    stft_unc = estimate_calibration_uncertainty(stft_chrono, target_fap=alpha_each)

    bb_cal = broadband_channel.calibrate_with_split(norm_obs, thresh_obs, alpha_channel=alpha_each)
    bb_chrono = np.array(
        [broadband_channel.score(o, bb_cal.band_medians, bb_cal.band_mad_scales) for o in thresh_obs],
        dtype=np.float64,
    )
    comb_fp = np.mean((stft_chrono > stft_cal.global_threshold) | (bb_chrono > bb_cal.threshold))
    return JointDetectorCalibration(
        target_fwer=target_fwer,
        alpha_per_family_sidak=alpha_each,
        stft_calibration=stft_cal,
        broadband_calibration=bb_cal,
        stft_uncertainty=stft_unc,
        empirical_combined_cal_far=float(comb_fp),
    )


@dataclass(frozen=True)
class GlitchReferenceState:
    """Reference baseline for short-timescale broadband impulsive glitch-storm indicator."""

    crest_median: float
    crest_mad: float
    single_window_z_thresh: float = 4.5
    elevated_window_z_thresh: float = 2.8


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
    stft_threshold: float
    stft_nominated_raw: bool
    stft_nominated_post_veto: bool
    broadband_excess_score: float
    broadband_threshold: float
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

    @property
    def is_clean_calibrated(self) -> bool:
        """True only when continuum is calibrated AND no glitch storm is active."""
        return (
            self.credibility_status == CRED_CALIBRATED
            and self.glitch_storm_status == GLITCH_QUIET
        )


class AdaptiveRegimeManagerV2b:
    """Phase 2b dual-output manager with staged causal calibration & orthogonal diagnostics."""

    def __init__(self, config: AdaptiveRegimeV2bConfig | None = None) -> None:
        self.config = config or AdaptiveRegimeV2bConfig()
        if self.config.mode not in ("causal_trailing", "offline_symmetric"):
            raise ValueError(f"Invalid mode: {self.config.mode!r}")
        if self.config.policy not in ("transient", "persistent"):
            raise ValueError(f"Invalid policy: {self.config.policy!r}")

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
        joint_calibration: JointDetectorCalibration,
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
        glitch_z_arr = np.array(
            [
                (extract_broadband_impulsive_crest(o) - glitch_ref.crest_median) / glitch_ref.crest_mad
                for o in observations
            ],
            dtype=np.float64,
        )

        # Frozen-reference broadband channel NEVER updates its reference R_0!
        bb_cal = joint_calibration.broadband_calibration
        active_ref = initial_continuum_ref
        active_stft_cal = joint_calibration.stft_calibration
        active_unc = joint_calibration.stft_uncertainty
        active_is_mature = True

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
            # --- Stream A: ALWAYS compute both raw scores & nominations ---
            det_res = detector.analyze(obs, active_stft_cal)
            raw_cands = tuple(det_res.candidates)
            if veto_config is not None:
                v_assess: VetoAssessment | None = apply_vetoes(obs, raw_cands, veto_config)
                surv_cands = tuple(v_assess.surviving_candidates)
            else:
                v_assess = None
                surv_cands = raw_cands

            bb_score = broadband_channel.score(obs, bb_cal.band_medians, bb_cal.band_mad_scales)
            bb_nom = bool(bb_score > bb_cal.threshold)
            stft_raw_nom = len(raw_cands) > 0
            stft_post_nom = len(surv_cands) > 0
            joint_nom = bool(stft_post_nom or bb_nom)

            # --- Stream B Channel 2: Orthogonal Impulsive Glitch-Storm Monitor ---
            if is_sym:
                g_lo = max(0, i - half_w)
                g_hi = min(n_obs, i + half_w + 1)
            else:
                g_lo = max(0, i - w + 1)
                g_hi = i + 1
            gz_i = float(glitch_z_arr[i])
            g_elev_rate = float(np.mean(glitch_z_arr[g_lo:g_hi] >= glitch_ref.elevated_window_z_thresh))
            glitch_storm_active = bool(
                gz_i >= glitch_ref.single_window_z_thresh
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
                    # Progressive Stage 1 -> Stage 2 causal upgrade (strictly past windows <= i)
                    if regime_id >= 1 and not active_is_mature:
                        regime_accum_indices.append(i)
                        if len(regime_accum_indices) >= cfg.mature_recal_windows:
                            m_obs = [observations[idx] for idx in regime_accum_indices[-cfg.mature_recal_windows :]]
                            active_ref = self.fit_continuum_reference(m_obs)
                            active_stft_cal = detector.calibrate(
                                m_obs,
                                target_fap=joint_calibration.alpha_per_family_sidak,
                                norm_fraction=0.25,
                                calibration_seed=regime_id * 10 + 2,
                            )
                            n_norm_m = max(20, int(round(len(m_obs) * 0.25)))
                            m_chrono = np.array(
                                [detector.analyze(o, active_stft_cal).max_global_score for o in m_obs[n_norm_m:]],
                                dtype=np.float64,
                            )
                            active_unc = estimate_calibration_uncertainty(
                                m_chrono, target_fap=joint_calibration.alpha_per_family_sidak
                            )
                            active_is_mature = True

                    cred_out = CRED_CALIBRATED
                    cont_status = CONT_STABLE
                    sig_out = self._classify_significance(regime_id, active_unc, active_is_mature)

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
                            target_fap=joint_calibration.alpha_per_family_sidak,
                            norm_fraction=0.25,
                            calibration_seed=regime_id * 10 + 1,
                        )
                        n_norm_p = max(15, int(round(len(p_obs) * 0.25)))
                        p_chrono = np.array(
                            [detector.analyze(o, active_stft_cal).max_global_score for o in p_obs[n_norm_p:]],
                            dtype=np.float64,
                        )
                        active_unc = estimate_calibration_uncertainty(
                            p_chrono, target_fap=joint_calibration.alpha_per_family_sidak
                        )
                        regime_id += 1
                        active_is_mature = False
                        regime_accum_indices = list(recal_indices)
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                        recal_indices = []
                        state = CRED_CALIBRATED
                    else:
                        recal_indices = recal_indices[-cfg.stability_confirm_windows :]

            causes: list[str] = []
            if cont_status != CONT_STABLE:
                causes.append("CONTINUUM_DRIFT")
            if glitch_status != GLITCH_QUIET:
                causes.append("GLITCH_STORM")
            if bb_nom:
                causes.append("SINGLE_CHANNEL_BROADBAND_AMBIGUOUS")

            outputs.append(
                DualOutputWindowV2b(
                    index=i,
                    sample_id=obs.sample_id,
                    raw_stft_candidates=raw_cands,
                    surviving_stft_candidates=surv_cands,
                    stft_max_score=det_res.max_global_score,
                    stft_threshold=active_stft_cal.global_threshold,
                    stft_nominated_raw=stft_raw_nom,
                    stft_nominated_post_veto=stft_post_nom,
                    broadband_excess_score=bb_score,
                    broadband_threshold=bb_cal.threshold,
                    broadband_nominated=bb_nom,
                    joint_nominated_post_veto=joint_nom,
                    veto_assessment=v_assess,
                    regime_id=regime_id,
                    credibility_status=cred_out,
                    continuum_status=cont_status,
                    glitch_storm_status=glitch_status,
                    warning_causes=tuple(causes),
                    significance_status=sig_out,
                    calibration_uncertainty=active_unc,
                    mode=cfg.mode,
                    policy=cfg.policy,
                    future_data_used=is_sym,
                    min_latency_windows=min_lat_wins,
                    glitch_crest_z=gz_i,
                    max_feature_z_drift=max_z_drift,
                    diagnostics={"is_mature": active_is_mature, "glitch_elev_rate": g_elev_rate},
                )
            )

        return outputs
