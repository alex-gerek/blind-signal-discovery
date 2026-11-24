"""Automated pytest suite for BSD v0.3 Phase 2 (AdaptiveRegimeManager)."""

from __future__ import annotations

import numpy as np
import pytest

from bsd.adaptive_regime import (
    CRED_CALIBRATED,
    SIG_EMPIRICAL,
    SIG_EXACT,
    SIG_UNCERTAIN,
    AdaptiveRegimeConfig,
    AdaptiveRegimeManager,
    summarize_dual_output_timeline,
)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, _scale_to_energy_snr


def _white_obs(rng: np.random.Generator, n: int, sigma: float = 1.0, prefix: str = "w") -> list[Observation]:
    return [Observation(f"{prefix}_{i}", rng.normal(0.0, sigma, size=4096), 1024.0, 4.0) for i in range(n)]


def test_finite_regime_recovery_and_significance_labels() -> None:
    """After an abrupt step-change, transient policy recovers in finite time with EMPIRICAL_FINITE_SAMPLE."""
    rng = np.random.default_rng(501)
    det = MultiResolutionDetector()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    cal = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.33, calibration_seed=501)

    mgr = AdaptiveRegimeManager(
        AdaptiveRegimeConfig(
            mode="causal_trailing",
            policy="transient",
            local_window_count=10,
            stability_confirm_windows=8,
            recal_window_count=60,
            min_exact_cal_samples=60,
        )
    )
    ref_state = mgr.fit_continuum_reference(ref_obs)

    # 40 stationary windows + 100 shifted windows (2.5x variance)
    seq = _white_obs(rng, 40, sigma=1.0, prefix="a") + _white_obs(rng, 100, sigma=2.5, prefix="b")
    out = mgr.process_sequence(seq, det, cal, ref_state)

    # Initial regime R_0 has CALIBRATED_EXACT (since n_thresh=60 >= min_exact_cal_samples=60)
    assert out[15].regime_id == 0
    assert out[15].significance_status == SIG_EXACT

    # Must recover into regime_id >= 1 within <= 75 windows of onset (idx=40)
    summ = summarize_dual_output_timeline(out, [False] * len(seq), onset_idx=40)
    assert summ["regime_recovery_windows"] is not None
    assert summ["regime_recovery_windows"] <= 75
    assert summ["n_regimes_established"] >= 2

    # Post-recovery windows MUST be marked EMPIRICAL_FINITE_SAMPLE, never CALIBRATED_EXACT
    rec_win = next(o for o in out[40:] if o.regime_id >= 1 and o.is_calibrated)
    assert rec_win.significance_status == SIG_EMPIRICAL
    assert rec_win.credibility_status == CRED_CALIBRATED


def test_raw_candidates_preserved_during_transition_and_persistent_policy_no_recal() -> None:
    """Candidates during transitions are never deleted; persistent policy never recalibrates."""
    rng = np.random.default_rng(602)
    det = MultiResolutionDetector()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    cal = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.33, calibration_seed=602)

    mgr_pers = AdaptiveRegimeManager(
        AdaptiveRegimeConfig(mode="causal_trailing", policy="persistent", local_window_count=10)
    )
    ref_state = mgr_pers.fit_continuum_reference(ref_obs)

    # Shifted noise with a strong pulse injected during transition at index 35
    seq = _white_obs(rng, 20, sigma=1.0, prefix="p1") + _white_obs(rng, 60, sigma=3.0, prefix="p2")
    t = np.arange(4096, dtype=np.float64) / 1024.0
    pulse = np.exp(-0.5 * ((t - 2.0) / 0.005) ** 2) * np.cos(2.0 * np.pi * 200.0 * (t - 2.0))
    pulse = _scale_to_energy_snr(pulse, 3.0, 25.0)
    seq[35] = Observation(seq[35].sample_id, seq[35].signal + pulse, 1024.0, 4.0)

    out = mgr_pers.process_sequence(seq, det, cal, ref_state)
    # Persistent policy MUST NEVER increment regime_id (no single-channel recalibration)
    assert all(o.regime_id == 0 for o in out)
    # Candidate at index 35 MUST still be preserved in Stream A even if flagged UNCERTAIN
    assert out[35].detected_raw is True
    assert len(out[35].raw_candidates) >= 1
    assert out[35].significance_status == SIG_UNCERTAIN


def test_continuum_features_preserve_persistent_narrowband_sinusoid() -> None:
    """An 80-s narrowband sinusoid should not trigger false broadband continuum shift."""
    rng = np.random.default_rng(703)
    det = MultiResolutionDetector()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    cal = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.33, calibration_seed=703)

    mgr = AdaptiveRegimeManager(
        AdaptiveRegimeConfig(mode="causal_trailing", policy="transient", local_window_count=12)
    )
    ref_state = mgr.fit_continuum_reference(ref_obs)

    seq = _white_obs(rng, 60, sigma=1.0, prefix="sin")
    t = np.arange(4096, dtype=np.float64) / 1024.0
    is_sig = [False] * 60
    for idx in range(20, 40):
        tone = _scale_to_energy_snr(np.sin(2.0 * np.pi * 180.0 * (t + idx * 4.0)), 1.0, 20.0)
        seq[idx] = Observation(seq[idx].sample_id, seq[idx].signal + tone, 1024.0, 4.0)
        is_sig[idx] = True

    out = mgr.process_sequence(seq, det, cal, ref_state)
    summ = summarize_dual_output_timeline(out, is_sig)
    # Narrowband sinusoid should retain >= 85% calibrated efficiency (vs 30% in Phase 1!)
    assert summ["sig_eff_calibrated_only"] is not None
    assert summ["sig_eff_calibrated_only"] >= 0.85


def test_causal_trailing_strict_causality_and_offline_symmetric_flag() -> None:
    """causal_trailing outputs at index k must be invariant to future windows j > k."""
    rng = np.random.default_rng(804)
    det = MultiResolutionDetector()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    cal = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.33, calibration_seed=804)

    mgr_c = AdaptiveRegimeManager(AdaptiveRegimeConfig(mode="causal_trailing", local_window_count=10))
    ref_state = mgr_c.fit_continuum_reference(ref_obs)

    seq_a = _white_obs(rng, 40, sigma=1.0, prefix="c")
    seq_b = list(seq_a)
    # Strongly perturb all future windows after index 24
    for j in range(25, 40):
        seq_b[j] = Observation(seq_b[j].sample_id, seq_b[j].signal * 10.0, 1024.0, 4.0)

    out_a = mgr_c.process_sequence(seq_a, det, cal, ref_state)
    out_b = mgr_c.process_sequence(seq_b, det, cal, ref_state)
    for k in range(25):
        assert out_a[k].future_data_used is False
        assert out_a[k].credibility_status == out_b[k].credibility_status
        assert out_a[k].regime_id == out_b[k].regime_id
        assert out_a[k].cusum_stat == pytest.approx(out_b[k].cusum_stat)
        assert out_a[k].max_feature_z_drift == pytest.approx(out_b[k].max_feature_z_drift)

    mgr_s = AdaptiveRegimeManager(AdaptiveRegimeConfig(mode="offline_symmetric", local_window_count=10))
    out_s = mgr_s.process_sequence(seq_a, det, cal, ref_state)
    assert out_s[15].future_data_used is True
    assert out_s[15].min_latency_windows == 6


def test_exact_significance_denied_when_initial_reference_small() -> None:
    """CALIBRATED_EXACT must not be assigned when n_threshold_samples < min_exact_cal_samples."""
    rng = np.random.default_rng(905)
    det = MultiResolutionDetector()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    cal = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.33, calibration_seed=905)
    assert cal.n_threshold_samples < 200

    mgr = AdaptiveRegimeManager(AdaptiveRegimeConfig(mode="causal_trailing", min_exact_cal_samples=200))
    ref_state = mgr.fit_continuum_reference(ref_obs)
    seq = _white_obs(rng, 25, sigma=1.0, prefix="sm")
    out = mgr.process_sequence(seq, det, cal, ref_state)
    assert out[15].regime_id == 0
    assert out[15].credibility_status == CRED_CALIBRATED
    assert out[15].significance_status == SIG_EMPIRICAL
