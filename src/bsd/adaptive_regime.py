"""BSD v0.3 Phase 2 — AdaptiveRegimeManager (Dual-Output & Finite Regime Recovery).

Implements:
  1. Spectral-continuum background features (frequency-median-filtered STFT
     continuum PSD + robust tail ratio) treated as a testable hypothesis for
     reducing narrowband persistent-signal absorption.
  2. Dual-output architecture separating:
     - Stream A (Anomaly Stream): 100% of raw candidates preserved verbatim
       (including during transitions and recalibration warmup), with reversible
       veto annotations.
     - Stream B (Credibility Stream): independent regime state, recovery status,
       and strict significance classification (`CALIBRATED_EXACT` ONLY for large
       stationary initial reference ensembles; `EMPIRICAL_FINITE_SAMPLE` after
       local post-transition recalibration; `UNCERTAIN_UNVALIDATED` during
       transitions/warmup).
  3. Three-state finite-state machine with hysteresis (`CALIBRATED_STATIONARY`
     -> `TRANSITION_UNCREDIBLE` -> `RECALIBRATING_WARMUP` -> `CALIBRATED_STATIONARY`)
     under configurable `transient` vs. `persistent` policies and `causal_trailing`
     vs. `offline_symmetric` modes.

Does NOT modify frozen BSD v0.1 (`detector.py`) or v0.2 (`detector_v2.py`,
`vetoes.py`, `coincidence.py`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy import signal as sp_signal, stats

from bsd.calibration import conformal_order_rank
from bsd.detector_v2 import (
    MultiResolutionCalibration,
    MultiResolutionCandidate,
    MultiResolutionDetector,
)
from bsd.generator import Observation
from bsd.regime_monitor import DEFAULT_FREQ_BANDS_HZ, ReferenceBackgroundState
from bsd.vetoes import VetoAssessment, VetoConfig, apply_vetoes

CRED_CALIBRATED = "CALIBRATED_STATIONARY"
CRED_TRANSITION = "TRANSITION_UNCREDIBLE"
CRED_RECALIBRATING = "RECALIBRATING_WARMUP"
CRED_CALIBRATION_UNAVAILABLE = "CALIBRATION_UNAVAILABLE"

SIG_EXACT = "CALIBRATED_EXACT"
SIG_EMPIRICAL = "EMPIRICAL_FINITE_SAMPLE"
SIG_UNCERTAIN = "UNCERTAIN_UNVALIDATED"
SIG_UNAVAILABLE = "CALIBRATION_UNAVAILABLE"


@dataclass(frozen=True)
class AdaptiveRegimeConfig:
    """Configuration for AdaptiveRegimeManager."""

    mode: str = "causal_trailing"  # "causal_trailing" or "offline_symmetric"
    policy: str = "transient"  # "transient" or "persistent"
    transient_duration_horizon_s: float = 10.0  # Experimental boundary parameter
    local_window_count: int = 12
    ks_alpha: float = 0.001
    drift_z_threshold: float = 4.0
    cusum_k: float = 1.0
    cusum_h: float = 22.0
    shift_confirm_windows: int = 2
    stability_confirm_windows: int = 10
    recal_window_count: int = 60  # >=60 for split-sample MultiResolutionDetector calibration
    internal_stability_z_max: float = 3.0
    freq_median_kernel: int = 15
    freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ
    stft_nperseg: int = 512
    target_fap: float = 0.01
    min_exact_cal_samples: int = 200

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DualOutputWindow:
    """Non-destructive dual-output record for one sequential observation window."""

    index: int
    sample_id: str
    # Stream A: Anomaly & Candidate Output (never suppressed)
    raw_candidates: tuple[MultiResolutionCandidate, ...]
    surviving_candidates: tuple[MultiResolutionCandidate, ...]
    detected_raw: bool
    detected_post_veto: bool
    max_global_score: float
    active_threshold: float | None
    veto_assessment: VetoAssessment | None
    # Stream B: Background Stability & Credibility Output
    regime_id: int
    credibility_status: str
    significance_status: str
    min_resolvable_p: float | None
    mode: str
    policy: str
    future_data_used: bool
    min_latency_windows: int
    min_latency_seconds: float
    ks_min_p_bonferroni: float | None
    max_feature_z_drift: float
    cusum_stat: float
    diagnostics: dict[str, Any] = field(default_factory=dict)
    detection_available: bool = True

    @property
    def is_calibrated(self) -> bool:
        return (
            self.credibility_status == CRED_CALIBRATED
            and self.detection_available
            and self.active_threshold is not None
        )


def extract_continuum_features(
    observation: Observation,
    freq_bands_hz: tuple[tuple[float, float], ...] = DEFAULT_FREQ_BANDS_HZ,
    stft_nperseg: int = 512,
    freq_median_kernel: int = 15,
) -> np.ndarray:
    """Extract pure spectral-continuum background features via high-resolution frequency medians.

    Uses a 512-point Hann STFT (df = 2 Hz at 1024 Hz) and takes the median
    across frequency bins within each band and across all frequency bins per
    time frame (no sliding-window smearing of narrowband lines). Narrowband
    persistent tones occupy <5% of bins per band and leave the continuum
    unperturbed, whereas broadband noise shifts and broadband bursts shift the
    median continuum directly.

    Returns a 6-element vector:
      [log_continuum_total_med, log_continuum_q90_burst, log_continuum_band_0..3].
    """
    _ = freq_median_kernel  # Retained in signature for config compatibility
    if not isinstance(observation, Observation):
        raise TypeError(
            f"AdaptiveRegimeManager strictly requires an Observation, got {type(observation).__name__}."
        )
    x = np.asarray(observation.signal, dtype=np.float64)
    x = x - float(np.median(x))

    freqs_all, _times, zxx = sp_signal.stft(
        x,
        fs=observation.sample_rate,
        window="hann",
        nperseg=min(stft_nperseg, x.size),
        noverlap=min(stft_nperseg, x.size) // 2,
        boundary=None,
        padded=False,
    )
    freqs = freqs_all[1:-1]
    power = np.abs(zxx[1:-1, :]) ** 2
    raw_psd = np.median(power, axis=1) / float(np.log(2.0))
    raw_psd = np.maximum(raw_psd, 1e-30)

    log_total = float(np.log(max(float(np.median(raw_psd)), 1e-30)))
    # Median across all frequency bins per time frame -> immune to narrowband lines,
    # sensitive to broadband impulsive bursts.
    frame_cont_med = np.median(power, axis=0)
    log_q90_burst = float(np.log(max(float(np.quantile(frame_cont_med, 0.90)), 1e-30)))

    band_vals: list[float] = []
    for f_lo, f_hi in freq_bands_hz:
        mask = (freqs >= f_lo) & (freqs <= f_hi)
        val = float(np.median(raw_psd[mask])) if np.any(mask) else float(np.median(raw_psd))
        band_vals.append(float(np.log(max(val, 1e-30))))

    return np.array([log_total, log_q90_burst, *band_vals], dtype=np.float64)


class AdaptiveRegimeManager:
    """Dual-output adaptive regime manager with finite post-transition recovery."""

    def __init__(self, config: AdaptiveRegimeConfig | None = None) -> None:
        self.config = config or AdaptiveRegimeConfig()
        if self.config.mode not in ("causal_trailing", "offline_symmetric"):
            raise ValueError(f"Invalid mode: {self.config.mode!r}")
        if self.config.policy not in ("transient", "persistent"):
            raise ValueError(f"Invalid policy: {self.config.policy!r}")
        if self.config.recal_window_count < 60:
            raise ValueError("recal_window_count must be >= 60 for split-sample calibration.")

    def _feature_names(self) -> tuple[str, ...]:
        bnames = tuple(f"log_cont_{int(a)}_{int(b)}hz" for a, b in self.config.freq_bands_hz)
        return ("log_continuum_total", "log_continuum_q90_burst", *bnames)

    def fit_continuum_reference(self, reference_observations: list[Observation]) -> ReferenceBackgroundState:
        """Fit spectral-continuum reference state on stationary background windows."""
        if len(reference_observations) < 20:
            raise ValueError(f"Need >= 20 reference observations, got {len(reference_observations)}.")
        cfg = self.config
        mat = np.vstack(
            [
                extract_continuum_features(
                    o, cfg.freq_bands_hz, cfg.stft_nperseg, cfg.freq_median_kernel
                )
                for o in reference_observations
            ]
        )
        medians = np.median(mat, axis=0)
        mad_scales = 1.4826 * np.median(np.abs(mat - medians[None, :]), axis=0)
        std_fb = np.std(mat, axis=0)
        scales = np.where(mad_scales > 1e-9, mad_scales, np.where(std_fb > 1e-9, std_fb, 1.0))
        return ReferenceBackgroundState(
            feature_names=self._feature_names(),
            ref_matrix=mat,
            ref_medians=medians,
            ref_mad_scales=scales,
            n_reference_windows=len(reference_observations),
        )

    def _is_internally_stable(
        self,
        block: np.ndarray,
        scale_ref: np.ndarray,
    ) -> bool:
        """Check whether two halves of a candidate new-regime block are mutually consistent."""
        n = block.shape[0]
        if n < 6:
            return False
        mid = n // 2
        m1 = np.median(block[:mid, :], axis=0)
        m2 = np.median(block[mid:, :], axis=0)
        block_med = np.median(block, axis=0)
        block_mad = 1.4826 * np.median(np.abs(block - block_med[None, :]), axis=0)
        eff_scale = np.maximum(np.maximum(scale_ref, block_mad), 1e-9)
        z_diff = np.abs(m2 - m1) / eff_scale
        return bool(float(np.max(z_diff)) <= self.config.internal_stability_z_max)

    def _block_matches_reference(
        self,
        block: np.ndarray,
        ref: ReferenceBackgroundState,
    ) -> bool:
        """Return True if a candidate block has returned to agreement with `ref`."""
        n_features = block.shape[1]
        loc_med = np.median(block, axis=0)
        z_shift = np.abs(loc_med - ref.ref_medians) / ref.ref_mad_scales
        if float(np.max(z_shift)) >= self.config.drift_z_threshold:
            return False
        p_vals = [
            float(stats.ks_2samp(block[:, k], ref.ref_matrix[:, k]).pvalue)
            for k in range(n_features)
        ]
        ks_p_bonf = min(1.0, n_features * min(p_vals))
        return bool(ks_p_bonf >= self.config.ks_alpha)

    def process_sequence(
        self,
        observations: list[Observation],
        detector: MultiResolutionDetector,
        initial_calibration: MultiResolutionCalibration,
        initial_reference: ReferenceBackgroundState,
        veto_config: VetoConfig | None = None,
    ) -> list[DualOutputWindow]:
        """Run dual-output detection and regime management across sequential windows."""
        if not observations:
            return []
        cfg = self.config
        n_obs = len(observations)
        w = cfg.local_window_count
        is_sym = cfg.mode == "offline_symmetric"
        half_w = w // 2
        win_dur = float(observations[0].duration)
        min_lat_wins = (half_w + 1) if is_sym else 1
        min_lat_s = min_lat_wins * win_dur

        feats_mat = np.vstack(
            [
                extract_continuum_features(
                    o, cfg.freq_bands_hz, cfg.stft_nperseg, cfg.freq_median_kernel
                )
                for o in observations
            ]
        )
        n_features = feats_mat.shape[1]

        active_ref = initial_reference
        active_cal = initial_calibration
        regime_id = 0
        state = CRED_CALIBRATED
        cusum_s = 0.0
        consec_shift = 0
        regime_start_idx = 0

        transition_indices: list[int] = []
        recal_indices: list[int] = []
        outputs: list[DualOutputWindow] = []

        for i, obs in enumerate(observations):
            # Snapshot exact calibration state actually used to evaluate window i
            # (any recalibration triggered at window i takes effect strictly at window i + 1)
            used_cal = active_cal
            used_regime_id = regime_id
            cal_updated_for_next = False

            # 1. Stream A: ALWAYS run detector and preserve ALL raw candidates
            det_res = detector.analyze(obs, used_cal)
            raw_cands = tuple(det_res.candidates)
            if veto_config is not None:
                v_assess: VetoAssessment | None = apply_vetoes(obs, raw_cands, veto_config)
                surv_cands = tuple(v_assess.surviving_candidates)
            else:
                v_assess = None
                surv_cands = raw_cands

            # 2. Stream B: evaluate background stability against active_ref
            z_cur = np.abs(feats_mat[i] - active_ref.ref_medians) / active_ref.ref_mad_scales
            abs_z_i = float(np.max(z_cur))
            ks_p_bonf: float | None = None
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
                    ks_p_bonf = min(1.0, n_features * min(p_vals))
                    shifted = bool(
                        (ks_p_bonf < cfg.ks_alpha)
                        or (max_z_drift >= cfg.drift_z_threshold)
                        or (cusum_s >= cfg.cusum_h)
                    )

                if shifted:
                    consec_shift += 1
                else:
                    consec_shift = 0

                if shifted and consec_shift >= cfg.shift_confirm_windows:
                    if cfg.policy == "persistent":
                        # Persistent policy: flag uncertainty, NEVER adaptively recalibrate
                        cred_out = CRED_TRANSITION
                        sig_out = SIG_UNCERTAIN
                        min_p_out: float | None = None
                    else:
                        # Transient policy: enter TRANSITION_UNCREDIBLE and begin hysteresis check
                        state = CRED_TRANSITION
                        transition_indices = [i]
                        recal_indices = []
                        cusum_s = 0.0
                        consec_shift = 0
                        cred_out = CRED_TRANSITION
                        sig_out = SIG_UNCERTAIN
                        min_p_out = None
                else:
                    cred_out = CRED_CALIBRATED
                    # CALIBRATED_EXACT only for initial large stationary reference R_0
                    if regime_id == 0 and active_cal.n_threshold_samples >= cfg.min_exact_cal_samples:
                        sig_out = SIG_EXACT
                    else:
                        sig_out = SIG_EMPIRICAL
                    min_p_out = 1.0 / (active_cal.n_threshold_samples + 1.0)

            elif state == CRED_TRANSITION:
                transition_indices.append(i)
                cred_out = CRED_TRANSITION
                sig_out = SIG_UNCERTAIN
                min_p_out = None

                if len(transition_indices) >= cfg.stability_confirm_windows:
                    recent_idx = transition_indices[-cfg.stability_confirm_windows :]
                    block_new = feats_mat[recent_idx, :]
                    if self._block_matches_reference(block_new, active_ref):
                        # Temporary excursion returned to active reference: preserve active_cal
                        state = CRED_CALIBRATED
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                    elif self._is_internally_stable(block_new, active_ref.ref_mad_scales):
                        state = CRED_RECALIBRATING
                        recal_indices = list(recent_idx)

            else:  # state == CRED_RECALIBRATING
                recal_indices.append(i)
                cred_out = CRED_RECALIBRATING
                sig_out = SIG_UNCERTAIN
                min_p_out = None

                if len(recal_indices) >= cfg.recal_window_count:
                    recal_block = feats_mat[recal_indices, :]
                    if self._block_matches_reference(recal_block, active_ref):
                        state = CRED_CALIBRATED
                        regime_start_idx = i + 1
                        cusum_s = 0.0
                        consec_shift = 0
                        transition_indices = []
                        recal_indices = []
                    elif self._is_internally_stable(recal_block, active_ref.ref_mad_scales):
                        n_rec = len(recal_indices)
                        n_norm_rec = max(20, int(round(n_rec * 0.25)))
                        n_thresh_rec = n_rec - n_norm_rec
                        _, fap_ok, _ = conformal_order_rank(max(1, n_thresh_rec), cfg.target_fap)
                        if fap_ok and n_thresh_rec >= 40:
                            recal_obs = [observations[idx] for idx in recal_indices]
                            active_ref = self.fit_continuum_reference(recal_obs)
                            active_cal = detector.calibrate(
                                recal_obs,
                                target_fap=cfg.target_fap,
                                norm_fraction=0.25,
                                calibration_seed=regime_id + 1,
                                allow_unattainable_fap=False,
                            )
                            regime_id += 1
                            regime_start_idx = i + 1
                            cusum_s = 0.0
                            consec_shift = 0
                            transition_indices = []
                            recal_indices = []
                            state = CRED_CALIBRATED
                            cal_updated_for_next = True
                        else:
                            # Keep accumulating clean windows until target_fap becomes attainable
                            cred_out = CRED_CALIBRATION_UNAVAILABLE
                            sig_out = SIG_UNAVAILABLE
                            min_p_out = 1.0 / (n_thresh_rec + 1.0)
                    else:
                        # Regime still drifting: slide recalibration buffer forward
                        recal_indices = recal_indices[-cfg.stability_confirm_windows :]

            det_avail = bool(
                used_cal.target_fap_attainable
                and used_cal.global_threshold is not None
                and cred_out != CRED_CALIBRATION_UNAVAILABLE
            )
            if not det_avail and cred_out == CRED_CALIBRATED:
                cred_out = CRED_CALIBRATION_UNAVAILABLE
                sig_out = SIG_UNAVAILABLE

            outputs.append(
                DualOutputWindow(
                    index=i,
                    sample_id=obs.sample_id,
                    raw_candidates=raw_cands,
                    surviving_candidates=surv_cands,
                    detected_raw=(len(raw_cands) > 0) if det_avail else False,
                    detected_post_veto=(len(surv_cands) > 0) if det_avail else False,
                    max_global_score=det_res.max_global_score,
                    active_threshold=used_cal.global_threshold if det_avail else None,
                    veto_assessment=v_assess,
                    regime_id=used_regime_id,
                    credibility_status=cred_out,
                    significance_status=sig_out,
                    min_resolvable_p=min_p_out,
                    mode=cfg.mode,
                    policy=cfg.policy,
                    future_data_used=is_sym,
                    min_latency_windows=min_lat_wins,
                    min_latency_seconds=min_lat_s,
                    ks_min_p_bonferroni=ks_p_bonf,
                    max_feature_z_drift=max_z_drift,
                    cusum_stat=cusum_s,
                    diagnostics={
                        "active_n_thresh": used_cal.n_threshold_samples,
                        "detection_available": det_avail,
                        "recal_buffer_len": len(recal_indices),
                        "calibration_updated_for_next_window": cal_updated_for_next,
                    },
                    detection_available=det_avail,
                )
            )

        return outputs


def summarize_dual_output_timeline(
    outputs: list[DualOutputWindow],
    is_signal_flags: list[bool],
    onset_idx: int | None = None,
    window_duration_s: float = 4.0,
) -> dict[str, Any]:
    """Summarize full-timeline FAR, efficiency, coverage, and regime recovery time."""
    n = len(outputs)
    if len(is_signal_flags) != n:
        raise ValueError("outputs and is_signal_flags must have equal length.")

    cal_mask = [o.is_calibrated for o in outputs]
    n_cal_win = sum(cal_mask)
    n_excl = n - n_cal_win
    total_hr = (n * window_duration_s) / 3600.0
    excl_hr = (n_excl * window_duration_s) / 3600.0
    duty = n_cal_win / n if n > 0 else 0.0

    # Recovery time after onset_idx: first window >= onset_idx with regime_id >= 1 and is_calibrated
    recovery_windows: int | None = None
    if onset_idx is not None:
        for o in outputs[onset_idx:]:
            if o.regime_id >= 1 and o.is_calibrated:
                recovery_windows = o.index - onset_idx
                break

    # Background metrics (excluded windows NEVER counted as true negatives)
    bg_total = sum(1 for s in is_signal_flags if not s)
    bg_fp_full = sum(1 for o, s in zip(outputs, is_signal_flags) if (not s) and o.detected_post_veto)
    bg_cands_full = sum(len(o.surviving_candidates) for o, s in zip(outputs, is_signal_flags) if not s)

    bg_surv_n = sum(1 for o, s in zip(outputs, is_signal_flags) if (not s) and o.is_calibrated)
    bg_surv_fp = sum(
        1 for o, s in zip(outputs, is_signal_flags) if (not s) and o.is_calibrated and o.detected_post_veto
    )
    bg_surv_cands = sum(
        len(o.surviving_candidates) for o, s in zip(outputs, is_signal_flags) if (not s) and o.is_calibrated
    )
    surv_bg_hr = (bg_surv_n * window_duration_s) / 3600.0

    # Signal metrics (full timeline vs calibrated-only)
    sig_total = sum(1 for s in is_signal_flags if s)
    sig_tp_full_stream_a = sum(1 for o, s in zip(outputs, is_signal_flags) if s and o.detected_post_veto)
    sig_tp_calibrated_only = sum(
        1 for o, s in zip(outputs, is_signal_flags) if s and o.is_calibrated and o.detected_post_veto
    )
    sig_in_uncalibrated = sum(
        1 for o, s in zip(outputs, is_signal_flags) if s and (not o.is_calibrated) and o.detected_post_veto
    )

    n_regimes = max((o.regime_id for o in outputs), default=0) + 1
    cands_in_uncal = sum(len(o.raw_candidates) for o in outputs if not o.is_calibrated)

    return {
        "n_windows": n,
        "n_excluded_windows": n_excl,
        "total_livetime_hours": total_hr,
        "excluded_livetime_hours": excl_hr,
        "effective_duty_cycle": duty,
        "n_regimes_established": n_regimes,
        "regime_recovery_windows": recovery_windows,
        "regime_recovery_seconds": (recovery_windows * window_duration_s) if recovery_windows is not None else None,
        "raw_candidates_preserved_in_uncalibrated": cands_in_uncal,
        "bg_n_total": bg_total,
        "bg_far_full_timeline": (bg_fp_full / bg_total) if bg_total else None,
        "bg_cands_per_hr_full_timeline": (bg_cands_full / ((bg_total * window_duration_s) / 3600.0)) if bg_total else None,
        "bg_n_surviving_calibrated": bg_surv_n,
        "bg_far_calibrated_surviving": (bg_surv_fp / bg_surv_n) if bg_surv_n else None,
        "bg_cands_per_hr_calibrated_surviving": (bg_surv_cands / surv_bg_hr) if surv_bg_hr > 0 else None,
        "sig_n_total": sig_total,
        "sig_eff_full_stream_a": (sig_tp_full_stream_a / sig_total) if sig_total else None,
        "sig_eff_calibrated_only": (sig_tp_calibrated_only / sig_total) if sig_total else None,
        "sig_flagged_uncertain_during_transition": sig_in_uncalibrated,
    }
