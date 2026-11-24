"""Automated pytest suite for BSD v0.3 Phase 2b."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import signal as sp_signal

from bsd.adaptive_regime_v2b import (
    GLITCH_ELEVATED,
    GLITCH_QUIET,
    AdaptiveRegimeManagerV2b,
    AdaptiveRegimeV2bConfig,
    FrozenReferenceBroadbandChannel,
    calibrate_joint_stft_and_broadband,
    calibrate_separate_search_modes,
    estimate_calibration_uncertainty,
    fit_glitch_reference,
    is_glitch_free_reference_block,
)
from bsd.cross_coherence import (
    COH_ARTIFACT,
    COH_NULL,
    COH_PHYSICAL_EVIDENCE,
    COH_ZERO_LAG_ENV,
    CrossCoherenceConfig,
    evaluate_cross_coherence,
)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, _scale_to_energy_snr
from bsd.noise import generate_noise


def _white_obs(rng: np.random.Generator, n: int, sigma: float = 1.0, prefix: str = "w") -> list[Observation]:
    return [Observation(f"{prefix}_{i}", rng.normal(0.0, sigma, size=4096), 1024.0, 4.0) for i in range(n)]


def test_temporal_dependence_reduces_neff_and_widens_ci() -> None:
    """Correlated AR(1) calibration scores must reduce N_eff and flag dependence."""
    rng = np.random.default_rng(101)
    iid = rng.normal(0.0, 1.0, size=200)
    u_iid = estimate_calibration_uncertainty(iid, target_fap=0.01)
    assert u_iid.n_thresh == 200
    assert u_iid.n_eff >= 140.0
    assert u_iid.dependence_detected is False

    # Strong AR(1) phi=0.85
    ar = np.zeros(200, dtype=np.float64)
    for i in range(1, 200):
        ar[i] = 0.85 * ar[i - 1] + rng.normal(0.0, 0.5)
    u_ar = estimate_calibration_uncertainty(ar, target_fap=0.01)
    assert u_ar.dependence_detected is True
    assert u_ar.n_eff < 0.4 * u_ar.n_thresh
    assert (u_ar.far_ci90_neff[1] - u_ar.far_ci90_neff[0]) > (u_ar.far_ci90_iid[1] - u_ar.far_ci90_iid[0])


def test_joint_sidak_multiple_testing_control_and_broadband_nomination() -> None:
    """Sidak calibration controls combined FWER while broadband channel catches broadband persistent signal."""
    rng = np.random.default_rng(202)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(rng, 100, sigma=1.0, prefix="ref")
    jcal = calibrate_joint_stft_and_broadband(det, bb, ref_obs, target_fwer=0.01, norm_fraction=0.30, calibration_seed=202)

    assert jcal.alpha_per_family_sidak < jcal.target_fwer
    assert abs(1.0 - (1.0 - jcal.alpha_per_family_sidak) ** 2 - jcal.target_fwer) < 1e-9

    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="persistent"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    # 10 windows of stationary noise + 10 windows of 80-300 Hz broadband stochastic signal
    seq = _white_obs(rng, 20, sigma=1.0, prefix="bb")
    sos = sp_signal.butter(4, [80.0, 300.0], btype="bandpass", fs=1024.0, output="sos")
    for idx in range(10, 20):
        s = _scale_to_energy_snr(sp_signal.sosfiltfilt(sos, rng.normal(0.0, 1.0, size=4096)), 1.0, 26.0)
        seq[idx] = Observation(seq[idx].sample_id, seq[idx].signal + s, 1024.0, 4.0)

    out = mgr.process_sequence(seq, det, bb, jcal, cref, gref)
    bb_hits = sum(1 for o in out[10:20] if o.broadband_nominated)
    assert bb_hits >= 8
    assert all("SINGLE_CHANNEL_BROADBAND_AMBIGUOUS" in o.warning_causes for o in out[10:20] if o.broadband_nominated)


def test_orthogonal_glitch_storm_without_recalibration_and_sinusoid_preservation() -> None:
    """Glitch storm fires GLITCH_STORM without recalibrating; 80-s sinusoid fires neither warning."""
    rng = np.random.default_rng(303)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    jcal = calibrate_joint_stft_and_broadband(det, bb, ref_obs, target_fwer=0.02, calibration_seed=303)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    # Glitch storm sequence: 15 white + 20 impulsive
    gl_seq = _white_obs(rng, 15, sigma=1.0, prefix="gw") + [
        Observation(f"gi_{i}", generate_noise("impulsive", 4096, 1024.0, rng, sigma=1.0)[0], 1024.0, 4.0)
        for i in range(20)
    ]
    out_gl = mgr.process_sequence(gl_seq, det, bb, jcal, cref, gref)
    glitch_flagged = sum(1 for o in out_gl[15:] if o.glitch_storm_status == GLITCH_ELEVATED)
    assert glitch_flagged >= 15
    assert all("GLITCH_STORM" in o.warning_causes for o in out_gl[15:] if o.glitch_storm_status == GLITCH_ELEVATED)
    # Glitch storm alone must not force regime recalibration
    assert all(o.regime_id == 0 for o in out_gl)

    # 80-s narrowband sinusoid must trigger neither GLITCH_STORM nor CONTINUUM_DRIFT
    sin_seq = _white_obs(rng, 40, sigma=1.0, prefix="sin")
    t = np.arange(4096, dtype=np.float64) / 1024.0
    for idx in range(15, 35):
        tone = _scale_to_energy_snr(np.sin(2.0 * np.pi * 180.0 * (t + idx * 4.0)), 1.0, 20.0)
        sin_seq[idx] = Observation(sin_seq[idx].sample_id, sin_seq[idx].signal + tone, 1024.0, 4.0)
    out_sin = mgr.process_sequence(sin_seq, det, bb, jcal, cref, gref)
    clean_sin = sum(1 for o in out_sin[15:35] if o.is_clean_calibrated)
    assert clean_sin >= 18


def test_causal_trailing_strict_causality_v2b() -> None:
    """Perturbing future windows j > k must not alter Phase 2b causal_trailing output at k."""
    rng = np.random.default_rng(404)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="ref")
    jcal = calibrate_joint_stft_and_broadband(det, bb, ref_obs, target_fwer=0.02, calibration_seed=404)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    seq_a = _white_obs(rng, 35, sigma=1.0, prefix="c")
    seq_b = list(seq_a)
    for j in range(20, 35):
        seq_b[j] = Observation(seq_b[j].sample_id, seq_b[j].signal * 8.0, 1024.0, 4.0)

    out_a = mgr.process_sequence(seq_a, det, bb, jcal, cref, gref)
    out_b = mgr.process_sequence(seq_b, det, bb, jcal, cref, gref)
    for k in range(20):
        assert out_a[k].future_data_used is False
        assert out_a[k].credibility_status == out_b[k].credibility_status
        assert out_a[k].glitch_storm_status == out_b[k].glitch_storm_status
        assert out_a[k].warning_causes == out_b[k].warning_causes


def test_cross_coherence_five_regimes() -> None:
    """Cross-coherence must separate physical-lag signal from null, zero-lag env noise, and common artifacts."""
    rng = np.random.default_rng(505)
    cfg = CrossCoherenceConfig(n_time_slides=100, alpha_fap=0.01)
    n_win = 10
    n_tot = n_win * 4096
    sos = sp_signal.butter(4, [80.0, 300.0], btype="bandpass", fs=1024.0, output="sos")

    def _to_obs(arr: np.ndarray, pfx: str) -> list[Observation]:
        return [Observation(f"{pfx}_{k}", arr[k * 4096 : (k + 1) * 4096].copy(), 1024.0, 4.0) for k in range(n_win)]

    # 1. Independent noise
    r_null = evaluate_cross_coherence(_to_obs(rng.normal(0, 1, n_tot), "a"), _to_obs(rng.normal(0, 1, n_tot), "b"), cfg)
    assert r_null.classification == COH_NULL

    # 2. Correlated zero-lag broadband environmental noise
    shared_env = sp_signal.sosfiltfilt(sos, rng.normal(0, 1, n_tot)) * 0.6
    r_env = evaluate_cross_coherence(
        _to_obs(rng.normal(0, 1, n_tot) + shared_env, "e1"),
        _to_obs(rng.normal(0, 1, n_tot) + shared_env, "e2"),
        cfg,
    )
    assert r_env.classification == COH_ZERO_LAG_ENV

    # 3. Common narrowband spectral line artifact (180 Hz shared)
    t = np.arange(n_tot, dtype=np.float64) / 1024.0
    line = 0.5 * np.sin(2.0 * np.pi * 180.0 * t)
    r_art = evaluate_cross_coherence(
        _to_obs(rng.normal(0, 1, n_tot) + line, "l1"),
        _to_obs(rng.normal(0, 1, n_tot) + line, "l2"),
        cfg,
    )
    assert r_art.classification == COH_ARTIFACT

    # 4. Physical time-delay (+5 ms = 5 samples at 1024 Hz) shared broadband signal
    sig_bb = sp_signal.sosfiltfilt(sos, rng.normal(0, 1, n_tot)) * 0.6
    sig_bb_del = np.roll(sig_bb, 5)
    r_phys = evaluate_cross_coherence(
        _to_obs(rng.normal(0, 1, n_tot) + sig_bb, "p1"),
        _to_obs(rng.normal(0, 1, n_tot) + sig_bb_del, "p2"),
        cfg,
    )
    assert r_phys.classification == COH_PHYSICAL_EVIDENCE
    assert abs(abs(r_phys.peak_lag_ms) - (5.0 / 1024.0) * 1000.0) <= 1.5


def test_separate_search_modes_preserve_transient_budget_and_control_union() -> None:
    """Transient Search keeps full 1% FAP while Union applies explicit Sidak FWER correction."""
    rng = np.random.default_rng(606)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(rng, 120, sigma=1.0, prefix="sep")
    scal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.01, target_fap_persistent=0.01, target_fwer_union=0.01, calibration_seed=606
    )
    assert scal.transient_stft_cal.target_fap == pytest.approx(0.01)
    assert scal.persistent_bb_cal.alpha_channel == pytest.approx(0.01)
    assert scal.alpha_union_sidak < 0.01
    # Union threshold must be at least as strict as dedicated Transient / Persistent thresholds
    assert scal.union_stft_cal.global_threshold >= scal.transient_stft_cal.global_threshold
    assert scal.union_bb_cal.threshold >= scal.persistent_bb_cal.threshold


def test_glitch_ref_synchronizes_on_clean_transition_and_blocks_contamination() -> None:
    """glitch_ref updates on clean pink->shifted_ar2 transition and refuses impulsive contamination."""
    rng = np.random.default_rng(707)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = [Observation(f"r_{i}", generate_noise("pink_colored", 4096, 1024.0, rng, sigma=1.0)[0], 1024.0, 4.0) for i in range(90)]
    scal = calibrate_separate_search_modes(det, bb, ref_obs, calibration_seed=707)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    # 1. Contamination guard: impulsive glitch storm block must NOT pass clean check
    imp_block = [Observation(f"im_{i}", generate_noise("impulsive", 4096, 1024.0, rng, sigma=1.0)[0], 1024.0, 4.0) for i in range(60)]
    assert is_glitch_free_reference_block(imp_block, gref) is False

    # 2. Clean transition: 30 pink + 110 shifted_ar2 (1.5x) -> recovers R1 AND refreshes glitch_ref
    seq = [
        Observation(f"a_{i}", generate_noise("pink_colored", 4096, 1024.0, rng, sigma=1.0)[0], 1024.0, 4.0) for i in range(30)
    ] + [
        Observation(f"b_{i}", generate_noise("shifted_ar2", 4096, 1024.0, rng, sigma=1.5)[0], 1024.0, 4.0) for i in range(110)
    ]
    out = mgr.process_sequence(seq, det, bb, scal, cref, gref)
    assert out[-1].regime_id >= 1
    assert out[-1].diagnostics["glitch_ref_refresh_count"] >= 1
    # Post-recovery shifted_ar2 windows must be GLITCH_QUIET (not falsely flagged as glitch storm!)
    post_clean_frac = float(np.mean([o.is_clean_calibrated for o in out[105:]]))
    assert post_clean_frac >= 0.85


def test_strict_causality_and_candidate_preservation_with_separate_modes() -> None:
    """causal_trailing stays invariant to future windows and preserves transition candidates."""
    rng = np.random.default_rng(808)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(rng, 90, sigma=1.0, prefix="cr")
    scal = calibrate_separate_search_modes(det, bb, ref_obs, calibration_seed=808)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    seq_a = _white_obs(rng, 25, sigma=1.0, prefix="ca") + _white_obs(rng, 35, sigma=2.5, prefix="cb")
    t = np.arange(4096, dtype=np.float64) / 1024.0
    pulse = _scale_to_energy_snr(np.exp(-0.5 * ((t - 2.0) / 0.005) ** 2) * np.cos(2.0 * np.pi * 200.0 * (t - 2.0)), 2.5, 25.0)
    seq_a[32] = Observation(seq_a[32].sample_id, seq_a[32].signal + pulse, 1024.0, 4.0)

    seq_b = list(seq_a)
    for j in range(40, 60):
        seq_b[j] = Observation(seq_b[j].sample_id, seq_b[j].signal * 6.0, 1024.0, 4.0)

    out_a = mgr.process_sequence(seq_a, det, bb, scal, cref, gref)
    out_b = mgr.process_sequence(seq_b, det, bb, scal, cref, gref)
    for k in range(40):
        assert out_a[k].credibility_status == out_b[k].credibility_status
        assert out_a[k].glitch_storm_status == out_b[k].glitch_storm_status
        assert out_a[k].regime_id == out_b[k].regime_id
    # Raw candidate during transition at idx=32 must be preserved in Stream A
    assert out_a[32].stft_nominated_raw is True
    assert len(out_a[32].raw_stft_candidates) >= 1
