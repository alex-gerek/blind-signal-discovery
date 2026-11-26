"""Regression and edge-case tests for independent Codex audit rounds 1 & 2.

Covers:
  1. Unattainable target FAP in `calibration.py`, `detector_v2.py`,
     `adaptive_regime_v2b.py`, and `workstation_runner.py` (refusal by default,
     explicit `None` / strict RFC 8259 JSON `null` unavailable state when opted in,
     never `Infinity`/`NaN`, and exact floating-point boundary handling).
  2. Both adaptive managers (`AdaptiveRegimeManager` and `AdaptiveRegimeManagerV2b`)
     with insufficient vs. sufficient recalibration buffers, excluding unavailable
     windows from `is_clean_calibrated` / `detection_available`, and verifying
     actual signal detection sensitivity after recovery.
  3. Coincidence matching rejecting 2.8-second peak separation with overlapping
     intervals, zero cross-correlation yielding `delay_identifiable=False` and
     `within_10ms_light_travel=False`, and non-finite input rejection.
  4. Time-slide accounting: merging duplicate/overlapping `active_intervals`,
     misaligned window-start deduplication (by `window_gps_start`/`window_id`,
     NOT absolute GPS grid), zero-overlap / zero-slide / missing-interval
     `coincidence_far_per_hour is None` (never `0.0`), and `n_slides_performed`.
  5. Workstation user-facing attainability refusal, `INSUFFICIENT` status,
     strict JSON `allow_nan=False` export, and qualified UI strings.
  6. Proper scaled two-sided MC SE `2 * SE(p_one)` and non-finite guards.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from bsd.adaptive_regime import (
    CRED_CALIBRATED,
    CRED_CALIBRATION_UNAVAILABLE,
    CRED_RECALIBRATING,
    AdaptiveRegimeConfig,
    AdaptiveRegimeManager,
)
from bsd.adaptive_regime_v2b import (
    AdaptiveRegimeManagerV2b,
    AdaptiveRegimeV2bConfig,
    CalibrationUncertaintySummary,
    FrozenReferenceBroadbandChannel,
    calibrate_separate_search_modes,
    estimate_calibration_uncertainty,
    fit_glitch_reference,
)
from bsd.calibration import calibrate_detector, conformal_order_rank
from bsd.coincidence import (
    deduplicate_candidates_by_window,
    estimate_cross_detector_delay,
    match_coincident_candidates,
    run_time_slide_analysis,
)
from bsd.detector import TemplateFreeDetector
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, _scale_to_energy_snr
from bsd.io_loader import InputValidationError, LoadedTimeSeries, PreprocessingSummary
from bsd.workstation_runner import (
    CAL_STATUS_INSUFFICIENT,
    WorkstationConfig,
    WorkstationRunResult,
    _build_rows,
    classify_calibration_status,
    export_workstation_run,
    run_workstation_analysis,
)


def _white_obs(n: int, seed: int = 123, sigma: float = 1.0, prefix: str = "w") -> list[Observation]:
    rng = np.random.default_rng(seed)
    return [Observation(f"{prefix}_{i}", rng.normal(0.0, sigma, size=4096), 1024.0, 4.0) for i in range(n)]


def test_fap_floating_point_boundaries_and_invalid_observations() -> None:
    """Floating-point boundaries around 1/(N+1) are exact, and non-finite inputs fail before scoring."""
    # Exact rational boundary: n=6, target_fap=1/7 -> p_min=1/7 <= 1/7 (ATTAINABLE!)
    k6, ok6, _ = conformal_order_rank(6, 1.0 / 7.0)
    assert ok6 is True and k6 == 6

    # Exact decimal boundary: n=99, target_fap=0.01 -> ATTAINABLE; 0.01 - 1e-6 -> UNATTAINABLE
    k99, ok99, _ = conformal_order_rank(99, 0.01)
    assert ok99 is True and k99 == 99
    k99_lo, ok99_lo, _ = conformal_order_rank(99, 0.01 - 1e-6)
    assert ok99_lo is False and k99_lo is None
    k98, ok98, _ = conformal_order_rank(98, 0.01)
    assert ok98 is False and k98 is None

    # Invalid observations (NaN / inf / bad sample_rate) rejected BEFORE computing scores
    nan_sig = np.zeros(4096, dtype=np.float64)
    nan_sig[100] = np.nan
    bad_obs = Observation("bad_nan", nan_sig, 1024.0, 4.0)
    with pytest.raises(ValueError, match="non-finite"):
        TemplateFreeDetector().analyze(bad_obs)
    with pytest.raises(ValueError, match="non-finite"):
        MultiResolutionDetector()._raw_analyze_all_scales(bad_obs)
    with pytest.raises(ValueError, match="finite"):
        FrozenReferenceBroadbandChannel()._band_vec(bad_obs)


def test_unattainable_fap_explicit_null_state_and_json_safety() -> None:
    """Unavailable thresholds serialize as JSON null (never Infinity/NaN) and explicitly abstain."""
    obs_60 = _white_obs(60, seed=11)

    det_v1 = TemplateFreeDetector()
    with pytest.raises(ValueError, match="Unattainable target_fap"):
        calibrate_detector(det_v1, obs_60, target_fap=0.01)

    cal_v1 = calibrate_detector(det_v1, obs_60, target_fap=0.01, allow_unattainable_fap=True)
    assert cal_v1.target_fap_attainable is False
    assert cal_v1.threshold is None
    assert cal_v1.calibration_state == "UNAVAILABLE_FAP_UNATTAINABLE"
    # Strict RFC 8259 JSON serialization (allow_nan=False) MUST succeed
    s_v1 = json.dumps(cal_v1.to_summary_dict(), allow_nan=False)
    assert '"threshold": null' in s_v1
    res_v1 = det_v1.analyze(obs_60[0], calibration=cal_v1)
    assert res_v1.detected is False and res_v1.threshold_used is None
    assert res_v1.diagnostics["detection_available"] is False

    det_v2 = MultiResolutionDetector()
    with pytest.raises(ValueError, match="Unattainable target_fap"):
        det_v2.calibrate(obs_60, target_fap=0.01, norm_fraction=0.30)

    cal_v2 = det_v2.calibrate(obs_60, target_fap=0.01, norm_fraction=0.30, allow_unattainable_fap=True)
    assert cal_v2.target_fap_attainable is False
    assert cal_v2.global_threshold is None
    assert cal_v2.calibration_state == "UNAVAILABLE_FAP_UNATTAINABLE"
    s_v2 = json.dumps(cal_v2.to_summary_dict(), allow_nan=False)
    assert '"global_threshold": null' in s_v2
    res_v2 = det_v2.analyze(obs_60[0], cal_v2)
    assert res_v2.detected is False and res_v2.threshold_used is None
    assert res_v2.diagnostics["detection_available"] is False


def test_adaptive_managers_accumulate_until_attainable_and_detect_injected_signal() -> None:
    """Both adaptive managers exclude unavailable windows from clean exposure, accumulate until
    calibration becomes attainable, and genuinely detect an injected signal after recovery."""
    rng = np.random.default_rng(909)
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(90, seed=909, sigma=1.0, prefix="ref")
    cal_v1 = det.calibrate(ref_obs, target_fap=0.02, norm_fraction=0.30, calibration_seed=909)

    # Sequence: 20 R0 windows + 95 shifted R1 windows (sigma=2.5) with a strong burst injected at end
    t = np.arange(4096, dtype=np.float64) / 1024.0
    burst = _scale_to_energy_snr(
        np.exp(-0.5 * ((t - 2.0) / 0.006) ** 2) * np.cos(2.0 * np.pi * 190.0 * (t - 2.0)),
        noise_sigma=2.5,
        target_energy_snr=30.0,
    )

    # 1. AdaptiveRegimeManager (v0.3 Phase 2): target_fap=0.02 needs n_thresh>=49 -> n_total>=66 (at 25% norm)
    #    Case A (insufficient post-shift): only 50 shifted windows -> stays unavailable, never falsely calibrated
    seq_short = _white_obs(20, seed=910, sigma=1.0, prefix="a") + _white_obs(50, seed=911, sigma=2.5, prefix="b")
    mgr1 = AdaptiveRegimeManager(
        AdaptiveRegimeConfig(
            mode="causal_trailing",
            policy="transient",
            target_fap=0.02,
            recal_window_count=60,
            local_window_count=10,
            stability_confirm_windows=8,
        )
    )
    ref_state = mgr1.fit_continuum_reference(ref_obs)
    out_short = mgr1.process_sequence(seq_short, det, cal_v1, ref_state)
    assert out_short[-1].is_calibrated is False

    #    Case B (sufficient post-shift): 100 shifted windows -> accumulates >=66, recovers attainable threshold,
    #    and genuinely detects the injected SNR=30 burst in the final window!
    seq_long = _white_obs(20, seed=912, sigma=1.0, prefix="a") + _white_obs(100, seed=913, sigma=2.5, prefix="b")
    seq_long[-1] = Observation("b_injected", seq_long[-1].signal + burst, 1024.0, 4.0)
    out_long = mgr1.process_sequence(seq_long, det, cal_v1, ref_state)
    assert out_long[-1].regime_id >= 1
    assert out_long[-1].is_calibrated is True
    assert out_long[-1].detection_available is True
    assert out_long[-1].active_threshold is not None and np.isfinite(out_long[-1].active_threshold)
    assert out_long[-1].detected_raw is True

    # 2. AdaptiveRegimeManagerV2b: target_fwer=0.02 -> recal_fap=0.02 needs n_thresh>=49 (n_total>=66)
    scal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.02, target_fap_persistent=0.02, target_fwer_union=0.02, calibration_seed=909
    )
    mgr2 = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing",
            policy="transient",
            target_fwer=0.02,
            provisional_recal_windows=60,  # 60 windows -> n_thresh=45 < 49 -> initially UNAVAILABLE!
            mature_recal_windows=80,
            local_window_count=10,
            stability_confirm_windows=8,
        )
    )
    cref = mgr2.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)
    out_v2b = mgr2.process_sequence(seq_long, det, bb, scal, cref, gref)
    # Any window with unavailable threshold MUST be excluded from is_clean_calibrated & detection_available
    for w in out_v2b:
        if w.stft_threshold is None:
            assert w.detection_available is False
            assert w.is_clean_calibrated is False
            assert w.credibility_status in (CRED_CALIBRATION_UNAVAILABLE, CRED_RECALIBRATING)
    # After accumulating >=66 clean post-shift windows, recovery yields finite threshold & detects burst!
    assert out_v2b[-1].regime_id >= 1
    assert out_v2b[-1].detection_available is True
    assert out_v2b[-1].stft_threshold is not None and np.isfinite(out_v2b[-1].stft_threshold)
    assert out_v2b[-1].stft_nominated_raw is True


def test_coincidence_2p8s_separation_and_zero_correlation_unidentifiable() -> None:
    """2.8 s peak separation is rejected; zero correlation yields unidentifiable delay (NOT physical)."""
    ca = [{"detector": "H1", "gps_start": 1000.0, "gps_end": 1003.5, "gps_peak": 1000.2, "anomaly_score": 6.0}]
    cb = [{"detector": "L1", "gps_start": 1000.5, "gps_end": 1003.8, "gps_peak": 1003.0, "anomaly_score": 6.0}]
    assert abs(ca[0]["gps_peak"] - cb[0]["gps_peak"]) == pytest.approx(2.8)
    assert match_coincident_candidates(ca, cb, max_dt_s=0.05) == []

    # Zero arrays -> delay_identifiable=False, delay_s=None, within_10ms_light_travel=False!
    zeros = np.zeros(512, dtype=np.float64)
    res_zero = estimate_cross_detector_delay(zeros, zeros, sample_rate=1024.0, max_lag_s=0.015)
    assert res_zero.delay_identifiable is False
    assert res_zero.delay_s is None
    assert res_zero.within_10ms_light_travel is False
    assert res_zero.abs_peak_corr == 0.0

    # Non-finite input rejected
    nan_arr = np.zeros(512, dtype=np.float64)
    nan_arr[10] = np.nan
    with pytest.raises(ValueError, match="finite"):
        estimate_cross_detector_delay(nan_arr, zeros, sample_rate=1024.0)


def test_time_slide_misaligned_windows_duplicate_intervals_and_zero_overlap() -> None:
    """Time-slides merge duplicate intervals, deduplicate by window_gps_start on misaligned starts,
    and return FAR=None (UNAVAILABLE) on zero overlap, zero slides, or missing active_intervals."""
    # 1. Duplicate and overlapping active intervals must be merged (never double-counted!)
    dup_intervals = [(1000.0, 1200.0), (1000.0, 1200.0), (1100.0, 1300.0)]  # merges to [(1000.0, 1300.0)]
    clean_intervals = [(1000.0, 1300.0)]
    ts_dup = run_time_slide_analysis([], [], single_epoch_livetime_s=300.0, n_slides=10, slide_step_s=4.0, active_intervals=dup_intervals)
    ts_clean = run_time_slide_analysis([], [], single_epoch_livetime_s=300.0, n_slides=10, slide_step_s=4.0, active_intervals=clean_intervals)
    assert ts_dup["equivalent_time_slide_livetime_hours"] == pytest.approx(ts_clean["equivalent_time_slide_livetime_hours"])
    assert ts_dup["n_slides_performed"] == 10

    # 2. Misaligned window start (1001.5..1005.5): peaks at 1001.8 and 1004.2 share the SAME window_gps_start
    #    and MUST collapse to 1 candidate (not split across an absolute 4-s GPS grid!)
    misaligned_a = [
        {"detector": "L1", "window_gps_start": 1001.5, "gps_start": 1001.6, "gps_end": 1002.0, "gps_peak": 1001.8, "anomaly_score": 12.0},
        {"detector": "L1", "window_gps_start": 1001.5, "gps_start": 1004.0, "gps_end": 1004.5, "gps_peak": 1004.2, "anomaly_score": 8.0},
    ]
    ts_mis = run_time_slide_analysis(
        misaligned_a, [], single_epoch_livetime_s=300.0, n_slides=4, slide_step_s=4.0,
        active_intervals=[(1000.0, 1300.0)], deduplicate_window_s=4.0,
    )
    assert ts_mis["n_bg_cands_a"] == 1

    # 3. Zero overlap (shift 100 s on a 20 s interval) -> FAR is None (UNAVAILABLE), NEVER 0.0!
    ts_no_ov = run_time_slide_analysis(
        [], [], single_epoch_livetime_s=20.0, n_slides=2, slide_step_s=100.0, active_intervals=[(1000.0, 1020.0)]
    )
    assert ts_no_ov["equivalent_time_slide_livetime_hours"] == 0.0
    assert ts_no_ov["coincidence_far_per_hour"] is None
    assert ts_no_ov["far_status"] == "UNAVAILABLE_ZERO_EXPOSURE"

    # 4. Zero slides or missing active_intervals -> FAR is None (UNAVAILABLE), NEVER 0.0!
    ts_zero = run_time_slide_analysis([], [], single_epoch_livetime_s=604.0, n_slides=0)
    assert ts_zero["coincidence_far_per_hour"] is None
    assert ts_zero["far_status"] == "UNAVAILABLE_ZERO_SLIDES"
    ts_no_int = run_time_slide_analysis([], [], single_epoch_livetime_s=604.0, n_slides=10)
    assert ts_no_int["coincidence_far_per_hour"] is None
    assert ts_no_int["far_status"] == "UNAVAILABLE_NO_ACTIVE_INTERVALS"


def test_workstation_rejects_unattainable_fap_and_qualifies_ui_labels() -> None:
    """Workstation refuses unattainable target_fap with InputValidationError and qualifies UI text."""
    rng = np.random.default_rng(77)
    prep = PreprocessingSummary("synthetic", "ch0", "sigma", 1024.0, 1024.0, False, "none", (0.0, 512.0), (20.0, 450.0), False, ())
    bg_ts = LoadedTimeSeries("bg", "synthetic", "sha_bg", 80 * 4096, 1024.0, 320.0, rng.normal(0, 1, 80 * 4096), prep, (0.0, 320.0), "SYN")
    tg_ts = LoadedTimeSeries("tg", "synthetic", "sha_tg", 4 * 4096, 1024.0, 16.0, rng.normal(0, 1, 4 * 4096), prep, (500.0, 516.0), "SYN")

    with pytest.raises(InputValidationError, match="Unattainable target_fap"):
        run_workstation_analysis(tg_ts, bg_ts, WorkstationConfig(target_fap=0.01))

    u = CalibrationUncertaintySummary(70, 65.0, 0.01, 0.80, 1.0 / 71.0, 1.0 / 71.0, (0.001, 0.02), (0.001, 0.02), False)
    ca = classify_calibration_status(100, u, target_fap=0.01)
    assert ca.status == CAL_STATUS_INSUFFICIENT

    with pytest.raises(ValueError, match="finite"):
        estimate_calibration_uncertainty(np.array([1.0, np.nan, 2.0, 3.0, 4.0, 5.0]))

    app_src = (Path(__file__).resolve().parent.parent / "app" / "bsd_workstation.py").read_text(encoding="utf-8")
    assert "Ljung-Box Linear-Lag p" in app_src
    assert "does NOT guarantee real-sensor stationarity" in app_src
    assert "Det.-Avail. Duty Cycle" in app_src


def test_confirmatory_two_sided_mc_se_scaling_and_invalid_inputs() -> None:
    """Two-sided MC p-value `2*p_one` has SE `2*sqrt(p_one*(1-p_one)/M)` and rejects NaN/out-of-range."""
    script_path = Path(__file__).resolve().parent.parent / "experiments" / "run_confirmatory_multiseed_220.py"
    spec = importlib.util.spec_from_file_location("confirmatory_mod", script_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    p_two, se_two = mod._compute_two_sided_mc_p_and_se(0.44642, 0.56436, 200000)
    assert p_two == pytest.approx(0.89284)
    assert se_two == pytest.approx(2.0 * np.sqrt(0.44642 * (1.0 - 0.44642) / 200000.0), rel=1e-9)

    with pytest.raises(ValueError, match="finite probabilities"):
        mod._compute_two_sided_mc_p_and_se(float("nan"), 0.5, 1000)
    with pytest.raises(ValueError, match="positive"):
        mod._compute_two_sided_mc_p_and_se(0.1, 0.9, 0)


def _strong_chirp_obs(seed: int = 502, snr: float = 30.0, sigma: float = 1.0, sample_id: str = "strong_chirp") -> Observation:
    rng = np.random.default_rng(seed)
    noise = rng.normal(0.0, sigma, size=4096)
    t = np.arange(4096) / 1024.0
    mask = (t >= 1.5) & (t <= 2.3)
    t_loc = t[mask] - 1.5
    pulse = np.sin(2.0 * np.pi * (80.0 * t_loc + 80.0 * t_loc**2)) * np.hanning(t_loc.size)
    sig = noise.copy()
    sig[mask] += _scale_to_energy_snr(pulse, noise_sigma=sigma, target_energy_snr=snr)
    return Observation(sample_id, sig, 1024.0, 4.0)


def test_conformal_exact_nextafter_boundary_no_p_gt_alpha() -> None:
    """Exact np.nextafter(1/(N+1), 0) must NEVER permit detection with p_emp > target_fap."""
    obs_60 = _white_obs(60, seed=501)
    strong_sig = _strong_chirp_obs(seed=502, snr=30.0)

    # 1. v0.1 TemplateFreeDetector (N=60 -> p_min = 1/61)
    det_v1 = TemplateFreeDetector()
    alpha_below_v1 = float(np.nextafter(1.0 / 61.0, 0.0))
    alpha_exact_v1 = 1.0 / 61.0
    alpha_above_v1 = float(np.nextafter(1.0 / 61.0, 1.0))

    with pytest.raises(ValueError, match="Unattainable target_fap"):
        calibrate_detector(det_v1, obs_60, target_fap=alpha_below_v1, allow_unattainable_fap=False)
    cal_v1_below = calibrate_detector(det_v1, obs_60, target_fap=alpha_below_v1, allow_unattainable_fap=True)
    assert cal_v1_below.target_fap_attainable is False
    assert cal_v1_below.threshold is None
    res_v1_below = det_v1.analyze(strong_sig, calibration=cal_v1_below, return_maps=False)
    assert res_v1_below.detected is False
    assert res_v1_below.threshold_used is None
    assert res_v1_below.diagnostics["detection_available"] is False

    for alpha_ok in (alpha_exact_v1, alpha_above_v1):
        cal_ok = calibrate_detector(det_v1, obs_60, target_fap=alpha_ok, allow_unattainable_fap=False)
        assert cal_ok.target_fap_attainable is True
        assert cal_ok.threshold is not None
        res_ok = det_v1.analyze(strong_sig, calibration=cal_ok, return_maps=False)
        assert res_ok.detected is True
        p_emp, _ = cal_ok.empirical_p_value(res_ok.max_anomaly_score)
        assert p_emp <= alpha_ok

    # 2. v0.2 MultiResolutionDetector (N_total=60, norm_fraction=20/60 -> N_B=40, p_min=1/41)
    det_v2 = MultiResolutionDetector()
    alpha_below_v2 = float(np.nextafter(1.0 / 41.0, 0.0))
    alpha_exact_v2 = 1.0 / 41.0
    alpha_above_v2 = float(np.nextafter(1.0 / 41.0, 1.0))

    with pytest.raises(ValueError, match="Unattainable target_fap"):
        det_v2.calibrate(obs_60, target_fap=alpha_below_v2, norm_fraction=20.0 / 60.0, allow_unattainable_fap=False)
    cal_v2_below = det_v2.calibrate(obs_60, target_fap=alpha_below_v2, norm_fraction=20.0 / 60.0, allow_unattainable_fap=True)
    assert cal_v2_below.n_threshold_samples == 40
    assert cal_v2_below.target_fap_attainable is False
    assert cal_v2_below.global_threshold is None
    res_v2_below = det_v2.analyze(strong_sig, cal_v2_below)
    assert res_v2_below.detected is False
    assert res_v2_below.threshold_used is None
    assert res_v2_below.diagnostics["detection_available"] is False

    for alpha_ok in (alpha_exact_v2, alpha_above_v2):
        cal_ok = det_v2.calibrate(obs_60, target_fap=alpha_ok, norm_fraction=20.0 / 60.0, allow_unattainable_fap=False)
        assert cal_ok.n_threshold_samples == 40
        assert cal_ok.target_fap_attainable is True
        assert cal_ok.global_threshold is not None
        res_ok = det_v2.analyze(strong_sig, cal_ok)
        assert res_ok.detected is True
        p_emp, _ = cal_ok.empirical_p_value(res_ok.max_global_score)
        assert p_emp <= alpha_ok


def test_v2b_low_fap_0p005_recovery_and_resource_cap() -> None:
    """v2b recovers at FAP=0.005 (needs N_thresh>=199 -> N_total>=265 > mature_recal_windows=200)
    and detects an SNR=30 signal after recovery; resource cap blocks recovery cleanly."""
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(320, seed=601, sigma=1.0, prefix="ref")
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.005, target_fap_persistent=0.005,
        target_fwer_union=0.01, calibration_seed=1,
    )
    assert sep_cal.transient_stft_cal.target_fap_attainable is True
    assert sep_cal.transient_stft_cal.global_threshold is not None

    # 350 shifted-regime windows (amplitude 4.0) + final SNR=30 chirp in the shifted regime
    shifted = _white_obs(350, seed=602, sigma=4.0, prefix="sh")
    injected_obs = _strong_chirp_obs(seed=603, snr=30.0, sigma=4.0, sample_id="v2b_low_fap_inj")
    seq = shifted + [injected_obs]

    mgr = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing", policy="transient",
            provisional_recal_windows=60, mature_recal_windows=200,
        )
    )
    c_ref = mgr.fit_continuum_reference(ref_obs)
    g_ref = fit_glitch_reference(ref_obs)
    out = mgr.process_sequence(seq, det, bb, sep_cal, c_ref, g_ref)

    # Must recover past 200 windows (~265+) with finite threshold, is_mature=True, and DETECT the SNR=30 chirp!
    assert out[-1].regime_id >= 1
    assert out[-1].detection_available is True
    assert out[-1].stft_threshold is not None and np.isfinite(out[-1].stft_threshold)
    assert out[-1].diagnostics["is_mature"] is True
    assert out[-1].diagnostics["effective_mature_recal_windows"] >= 265
    assert out[-1].stft_nominated_raw is True

    # Resource-limit cap smaller than required 265 windows -> stays explicitly CALIBRATION_UNAVAILABLE
    mgr_capped = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing", policy="transient",
            provisional_recal_windows=60, mature_recal_windows=200, max_recal_buffer_windows=150,
        )
    )
    out_cap = mgr_capped.process_sequence(seq[:180], det, bb, sep_cal, c_ref, g_ref)
    assert out_cap[-1].detection_available is False
    assert out_cap[-1].stft_threshold is None
    assert out_cap[-1].credibility_status == CRED_CALIBRATION_UNAVAILABLE
    assert out_cap[-1].diagnostics["is_mature"] is False
    assert out_cap[-1].diagnostics["resource_limit_blocked_recovery"] is True


def test_calibration_switch_window_records_used_threshold_and_applies_next() -> None:
    """At every calibration switch window i, window i records the pre-update threshold/availability
    actually used to evaluate window i, and the new calibration takes effect starting at window i + 1."""
    det = MultiResolutionDetector()
    ref_obs = _white_obs(160, seed=701, sigma=1.0, prefix="ref")
    shifted = _white_obs(160, seed=702, sigma=4.0, prefix="sh")

    # 1. AdaptiveRegimeManager (v1)
    cal0 = det.calibrate(ref_obs, target_fap=0.01, norm_fraction=0.30, calibration_seed=1)
    mgr1 = AdaptiveRegimeManager(AdaptiveRegimeConfig(mode="causal_trailing", policy="transient", recal_window_count=60, target_fap=0.01))
    out1 = mgr1.process_sequence(shifted, det, cal0, mgr1.fit_continuum_reference(ref_obs))
    switch_idxs1 = [w.index for w in out1 if w.diagnostics.get("calibration_updated_for_next_window")]
    assert len(switch_idxs1) >= 1
    si = switch_idxs1[0]
    assert si + 1 < len(out1)
    # Window si was evaluated BEFORE the new calibration took effect (regime 0, recalibrating)
    assert out1[si].regime_id == 0
    assert out1[si].credibility_status in (CRED_RECALIBRATING, CRED_CALIBRATION_UNAVAILABLE)
    # Window si + 1 is the first window evaluated with regime 1 and the new finite threshold
    assert out1[si + 1].regime_id == 1
    assert out1[si + 1].detection_available is True
    assert out1[si + 1].active_threshold is not None and np.isfinite(out1[si + 1].active_threshold)

    # 2. AdaptiveRegimeManagerV2b
    bb = FrozenReferenceBroadbandChannel()
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.02, target_fap_persistent=0.02,
        target_fwer_union=0.04, calibration_seed=1,
    )
    mgr2 = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient", provisional_recal_windows=60, mature_recal_windows=100)
    )
    out2 = mgr2.process_sequence(shifted, det, bb, sep_cal, mgr2.fit_continuum_reference(ref_obs), fit_glitch_reference(ref_obs))
    switch_idxs2 = [w.index for w in out2 if w.diagnostics.get("calibration_updated_for_next_window")]
    assert len(switch_idxs2) >= 2  # Stage 0->1 (provisional) and Stage 1->2 (recovery)
    s_prov = switch_idxs2[0]
    assert out2[s_prov].regime_id == 0
    assert out2[s_prov].stft_threshold == sep_cal.transient_stft_cal.global_threshold
    assert out2[s_prov + 1].regime_id == 1
    assert out2[s_prov + 1].stft_threshold is None  # provisional 60 windows unattainable at FAP=0.02
    assert out2[s_prov + 1].detection_available is False

    s_rec = switch_idxs2[1]
    assert out2[s_rec].stft_threshold is None
    assert out2[s_rec].detection_available is False
    assert out2[s_rec + 1].stft_threshold is not None and np.isfinite(out2[s_rec + 1].stft_threshold)
    assert out2[s_rec + 1].detection_available is True


def test_adaptive_v1_preserves_diagnostic_candidates_during_abstention() -> None:
    """During CRED_CALIBRATION_UNAVAILABLE abstention, AdaptiveRegimeManager preserves Stream A
    diagnostic candidates while setting detected_raw=False, detection_available=False, active_threshold=None."""
    det = MultiResolutionDetector()
    ref_obs = _white_obs(160, seed=711, sigma=1.0, prefix="ref")
    cal0 = det.calibrate(ref_obs, target_fap=0.01, norm_fraction=0.30, calibration_seed=1)

    # Shifted regime with a strong chirp injected right during the unavailable recalibration accumulation
    seq = _white_obs(95, seed=712, sigma=4.0, prefix="sh")
    seq[80] = _strong_chirp_obs(seed=713, snr=30.0, sigma=4.0, sample_id="abstain_inj")

    mgr = AdaptiveRegimeManager(
        AdaptiveRegimeConfig(mode="causal_trailing", policy="transient", recal_window_count=60, target_fap=0.01)
    )
    out = mgr.process_sequence(seq, det, cal0, mgr.fit_continuum_reference(ref_obs))
    w_inj = out[80]
    assert w_inj.credibility_status == CRED_CALIBRATION_UNAVAILABLE
    assert w_inj.detection_available is False
    assert w_inj.active_threshold is None
    assert w_inj.detected_raw is False
    assert w_inj.detected_post_veto is False
    assert len(w_inj.raw_candidates) > 0


def test_workstation_plot_timeline_handles_mixed_and_all_none_thresholds() -> None:
    """Workstation _plot_timeline must not raise TypeError when stft_threshold mixes float and None."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from dataclasses import replace

    app_path = Path(__file__).resolve().parent.parent / "app" / "bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_workstation_app", app_path)
    assert spec is not None and spec.loader is not None
    app_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app_mod)

    # Run a tiny valid workstation analysis, then replace window thresholds with mixed and all-None
    rng = np.random.default_rng(801)
    prep = PreprocessingSummary("synthetic", "ch0", "sigma", 1024.0, 1024.0, False, "none", (0.0, 512.0), (20.0, 450.0), False, ())
    bg_ts = LoadedTimeSeries("bg", "synthetic", "sha_bg", 160 * 4096, 1024.0, 640.0, rng.normal(0, 1, 160 * 4096), prep, (0.0, 640.0), "SYN")
    tg_ts = LoadedTimeSeries("tg", "synthetic", "sha_tg", 4 * 4096, 1024.0, 16.0, rng.normal(0, 1, 4 * 4096), prep, (700.0, 716.0), "SYN")
    res = run_workstation_analysis(tg_ts, bg_ts, WorkstationConfig(target_fap=0.01))

    # Mixed float and None thresholds
    mixed_wins = [
        replace(w, stft_threshold=None, detection_available=False, warning_causes=("CALIBRATION_UNAVAILABLE",))
        if idx in (1, 2) else w
        for idx, w in enumerate(res.windows)
    ]
    res_mixed = replace(res, windows=mixed_wins)
    fig_mixed = app_mod._plot_timeline(res_mixed)
    assert fig_mixed is not None
    plt.close(fig_mixed)

    # All-None thresholds
    all_none_wins = [
        replace(w, stft_threshold=None, detection_available=False, warning_causes=("CALIBRATION_UNAVAILABLE",))
        for w in res.windows
    ]
    res_none = replace(res, windows=all_none_wins)
    fig_none = app_mod._plot_timeline(res_none)
    assert fig_none is not None
    plt.close(fig_none)


def test_foreground_background_window_dedup_symmetry() -> None:
    """One multi-scale transient (3 clusters in same 4-s window) produces 1 foreground event
    and 1 background event under deduplicate_candidates_by_window."""
    c_a_raw = [
        {"detector": "H1", "window_gps_start": 1000.0, "gps_start": 1001.4, "gps_end": 1001.6, "gps_peak": 1001.50, "anomaly_score": 10.0},
        {"detector": "H1", "window_gps_start": 1000.0, "gps_start": 1001.3, "gps_end": 1001.7, "gps_peak": 1001.51, "anomaly_score": 8.0},
        {"detector": "H1", "window_gps_start": 1000.0, "gps_start": 1001.2, "gps_end": 1001.8, "gps_peak": 1001.52, "anomaly_score": 6.0},
    ]
    c_b_raw = [
        {"detector": "L1", "window_gps_start": 1000.0, "gps_start": 1001.4, "gps_end": 1001.6, "gps_peak": 1001.505, "anomaly_score": 9.5},
        {"detector": "L1", "window_gps_start": 1000.0, "gps_start": 1001.3, "gps_end": 1001.7, "gps_peak": 1001.515, "anomaly_score": 7.5},
        {"detector": "L1", "window_gps_start": 1000.0, "gps_start": 1001.2, "gps_end": 1001.8, "gps_peak": 1001.525, "anomaly_score": 5.5},
    ]
    c_a_win = deduplicate_candidates_by_window(c_a_raw, window_s=4.0)
    c_b_win = deduplicate_candidates_by_window(c_b_raw, window_s=4.0)
    assert len(c_a_win) == 1 and len(c_b_win) == 1
    fg_pairs = match_coincident_candidates(c_a_win, c_b_win, max_dt_s=0.05, allow_interval_gap=False)
    assert len(fg_pairs) == 1

    ts = run_time_slide_analysis(
        c_a_raw, c_b_raw, single_epoch_livetime_s=100.0, max_dt_s=0.05,
        n_slides=4, slide_step_s=4.0, active_intervals=[(1000.0, 1100.0)], deduplicate_window_s=4.0,
    )
    assert ts["n_bg_cands_a"] == 1 and ts["n_bg_cands_b"] == 1
    assert ts["total_accidental_coincidences"] == 0

    # Place B exactly 4.0 s earlier so shift +4.0 s aligns them -> exactly 1 window-level slide coincidence (NOT 9!)
    c_b_shifted = [
        {**c, "window_gps_start": 1004.0, "gps_start": c["gps_start"] + 4.0, "gps_end": c["gps_end"] + 4.0, "gps_peak": c["gps_peak"] + 4.0}
        for c in c_b_raw
    ]
    ts_hit = run_time_slide_analysis(
        c_a_raw, c_b_shifted, single_epoch_livetime_s=100.0, max_dt_s=0.05,
        n_slides=2, slide_step_s=4.0, active_intervals=[(1000.0, 1100.0)], deduplicate_window_s=4.0,
    )
    assert ts_hit["n_bg_cands_a"] == 1 and ts_hit["n_bg_cands_b"] == 1
    assert ts_hit["total_accidental_coincidences"] == 1
    assert ts_hit["significance_status"] == "AVAILABLE"
    assert ts_hit["min_resolvable_p_value"] == pytest.approx(1.0 / 3.0)


def test_v2b_union_stft_recovery_and_partial_union_reporting() -> None:
    """v2b recovers BOTH dedicated STFT and stricter Sidak Union STFT (FAP~0.00501 -> >=265 windows),
    reports PARTIAL_BROADBAND_ONLY while Union STFT is unavailable, and detects SNR=30 in Union after recovery."""
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(320, seed=920, sigma=1.0, prefix="ref")
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.01, target_fap_persistent=0.01,
        target_fwer_union=0.01, calibration_seed=1,
    )
    assert sep_cal.union_stft_cal.global_threshold is not None

    shifted = _white_obs(340, seed=921, sigma=4.0, prefix="sh")
    injected = _strong_chirp_obs(seed=922, snr=30.0, sigma=4.0, sample_id="union_inj")
    seq = shifted + [injected]

    mgr = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing", policy="transient",
            provisional_recal_windows=60, mature_recal_windows=200,
        )
    )
    c_ref = mgr.fit_continuum_reference(ref_obs)
    g_ref = fit_glitch_reference(ref_obs)
    out = mgr.process_sequence(seq, det, bb, sep_cal, c_ref, g_ref)

    # Verify intermediate partial state (dedicated STFT available, Union STFT still None) is reported explicitly
    partial_wins = [
        w for w in out
        if w.regime_id >= 1
        and w.diagnostics["dedicated_stft_available"] is True
        and w.diagnostics["union_stft_available"] is False
    ]
    assert len(partial_wins) >= 1
    for pw in partial_wins:
        assert pw.diagnostics["union_mode_status"] == "PARTIAL_BROADBAND_ONLY"
        assert "UNION_STFT_UNAVAILABLE" in pw.warning_causes
        assert pw.diagnostics["is_mature"] is False

    # After >=265 windows, Union STFT recovers, is_mature=True, and SNR=30 chirp is detected in Union!
    last = out[-1]
    assert last.diagnostics["dedicated_stft_available"] is True
    assert last.diagnostics["union_stft_available"] is True
    assert last.diagnostics["union_mode_status"] == "FULL_UNION_AVAILABLE"
    assert last.diagnostics["is_mature"] is True
    assert last.stft_nominated_post_veto is True
    assert last.joint_nominated_post_veto is True

    # Cap at 200 (enough for dedicated STFT 132, NOT enough for Union 265): partial Union explicitly flagged!
    mgr_cap = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing", policy="transient",
            provisional_recal_windows=60, mature_recal_windows=200, max_recal_buffer_windows=200,
        )
    )
    out_cap = mgr_cap.process_sequence(seq[:240], det, bb, sep_cal, c_ref, g_ref)
    assert out_cap[-1].diagnostics["dedicated_stft_available"] is True
    assert out_cap[-1].diagnostics["union_stft_available"] is False
    assert out_cap[-1].diagnostics["union_mode_status"] == "PARTIAL_BROADBAND_ONLY"
    assert out_cap[-1].diagnostics["resource_limit_blocked_union"] is True
    assert out_cap[-1].diagnostics["is_mature"] is False


def test_v2b_provisional_audit_split_matches_calibrator_n_norm_samples() -> None:
    """Provisional 60-window calibration uses n_norm=20, n_thresh=40 -> audit must report N=40 and p_min=1/41 (NOT 45 / 1/46)."""
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(160, seed=930, sigma=1.0, prefix="ref")
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.02, target_fap_persistent=0.02,
        target_fwer_union=0.04, calibration_seed=1,
    )
    shifted = _white_obs(100, seed=931, sigma=4.0, prefix="sh")
    mgr = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient", provisional_recal_windows=60, mature_recal_windows=120)
    )
    out = mgr.process_sequence(shifted, det, bb, sep_cal, mgr.fit_continuum_reference(ref_obs), fit_glitch_reference(ref_obs))
    r1_wins = [w for w in out if w.regime_id == 1]
    assert len(r1_wins) >= 1
    first_r1 = r1_wins[0]
    assert first_r1.diagnostics["stft_n_norm_samples"] == 20
    assert first_r1.diagnostics["stft_n_threshold_samples"] == 40
    assert first_r1.calibration_uncertainty.n_thresh == 40
    assert first_r1.calibration_uncertainty.p_min == pytest.approx(1.0 / 41.0)


def test_v2b_max_recal_buffer_windows_validation_and_insufficient_positive() -> None:
    """Zero, negative, bool, or non-int max_recal_buffer_windows is rejected; small positive limit blocks recovery."""
    for bad_val in (0, -10, True, 12.5):
        with pytest.raises(ValueError, match="max_recal_buffer_windows"):
            AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(max_recal_buffer_windows=bad_val))  # type: ignore[arg-type]

    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(160, seed=940, sigma=1.0, prefix="ref")
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.01, target_fap_persistent=0.01,
        target_fwer_union=0.02, calibration_seed=1,
    )
    shifted = _white_obs(100, seed=941, sigma=4.0, prefix="sh")
    mgr_tiny = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(provisional_recal_windows=60, mature_recal_windows=200, max_recal_buffer_windows=30)
    )
    out_tiny = mgr_tiny.process_sequence(shifted, det, bb, sep_cal, mgr_tiny.fit_continuum_reference(ref_obs), fit_glitch_reference(ref_obs))
    assert out_tiny[-1].regime_id == 0
    assert out_tiny[-1].detection_available is False
    assert out_tiny[-1].credibility_status == CRED_CALIBRATION_UNAVAILABLE
    assert out_tiny[-1].diagnostics["resource_limit_blocked_recovery"] is True
    assert out_tiny[-1].diagnostics["is_mature"] is False


def test_time_slide_separates_technical_resolution_from_unavailable_significance() -> None:
    """Empty or missing active_intervals must return min_resolvable_p_value=None while keeping technical_slide_count_resolution."""
    ts_empty = run_time_slide_analysis(
        [], [], single_epoch_livetime_s=100.0, n_slides=2, slide_step_s=4.0, active_intervals=[]
    )
    assert ts_empty["far_status"] == "UNAVAILABLE_ZERO_EXPOSURE"
    assert ts_empty["coincidence_far_per_hour"] is None
    assert ts_empty["significance_status"] == "UNAVAILABLE_ZERO_EXPOSURE"
    assert ts_empty["min_resolvable_p_value"] is None
    assert ts_empty["technical_slide_count_resolution"] == pytest.approx(1.0 / 3.0)

    ts_none = run_time_slide_analysis([], [], single_epoch_livetime_s=100.0, n_slides=4, slide_step_s=4.0)
    assert ts_none["far_status"] == "UNAVAILABLE_NO_ACTIVE_INTERVALS"
    assert ts_none["significance_status"] == "UNAVAILABLE_NO_ACTIVE_INTERVALS"
    assert ts_none["min_resolvable_p_value"] is None
    assert ts_none["technical_slide_count_resolution"] == pytest.approx(1.0 / 5.0)


def test_benchmark_delay_none_guard_and_gw170814_report_accuracy() -> None:
    """Benchmark guards xc.delay_s is not None on GW150914 curve; GW170814 report accurately describes current dedup method."""
    root = Path(__file__).resolve().parent.parent
    bench_src = (root / "experiments" / "run_v02_phase3_coincidence_benchmark.py").read_text(encoding="utf-8")
    assert "and xc.delay_identifiable\n" in bench_src or "xc.delay_identifiable" in bench_src
    assert "xc.delay_s is not None" in bench_src

    gw_src = (root / "experiments" / "run_gw170814_exploratory_search.py").read_text(encoding="utf-8")
    assert "Window-level deduplicated time-slide counts (current audited method)" in gw_src
    assert "Raw cluster-level counts in" not in gw_src


def test_workstation_export_abstention_active_decision_and_union_states(tmp_path: Path) -> None:
    """End-to-end CSV/JSON export checks:
    1. SNR=30 injection during abstention (max_recal_buffer_windows=30): diagnostic candidates preserved,
       active_mode_anomaly_flag=False, sum(active_mode_anomaly_flag) == n_active_mode_anomaly_windows.
    2. Transient & persistent active decisions match mode and detection_available.
    3. Export & read-back of FULL_UNION_AVAILABLE, PARTIAL_BROADBAND_ONLY, UNION_UNAVAILABLE, and
       calibration-switch window state in both timeline.csv and run_report.json (strict JSON null).
    """
    import csv
    from dataclasses import replace

    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    ref_obs = _white_obs(160, seed=960, sigma=1.0, prefix="ref")
    sep_cal = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.01, target_fap_persistent=0.01,
        target_fwer_union=0.02, calibration_seed=1,
    )
    # 1. Abstention with max_recal_buffer_windows=30 + SNR=30 chirp injected at window 60 (during abstention)
    seq_abstain = _white_obs(80, seed=961, sigma=4.0, prefix="sh")
    seq_abstain[60] = _strong_chirp_obs(seed=962, snr=30.0, sigma=4.0, sample_id="inj_during_abstention")
    mgr_cap = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode="causal_trailing", policy="transient",
            provisional_recal_windows=60, mature_recal_windows=200, max_recal_buffer_windows=30,
        )
    )
    wins_abstain = mgr_cap.process_sequence(
        seq_abstain, det, bb, sep_cal, mgr_cap.fit_continuum_reference(ref_obs), fit_glitch_reference(ref_obs)
    )
    cands_rows, tl_rows = _build_rows(wins_abstain, window_duration_s=4.0, search_mode="transient")
    n_anom_summary = sum(1 for w in wins_abstain if w.detection_available and w.stft_nominated_post_veto)

    # Check window 60: diagnostic candidates preserved, diagnostic nomination True, active decision False!
    w60_row = tl_rows[60]
    assert w60_row["detection_available"] is False
    assert w60_row["stft_nominated_raw"] is True
    assert w60_row["active_mode_diagnostic_nomination"] is True
    assert w60_row["active_mode_anomaly_flag"] is False
    assert w60_row["resource_limit_blocked_stft"] is True
    assert w60_row["resource_limit_blocked_union"] is True
    assert any(c["window_index"] == 60 and c["candidate_type"] == "STFT_DETECTION_CANDIDATE" for c in cands_rows)
    assert sum(1 for r in tl_rows if r["active_mode_anomaly_flag"]) == n_anom_summary

    # Build a WorkstationRunResult and export to tmp_path, then read back CSV and JSON!
    prep = PreprocessingSummary("synthetic", "ch0", "sigma", 1024.0, 1024.0, False, "none", (0.0, 512.0), (20.0, 450.0), False, ())
    bg_ts = LoadedTimeSeries("bg", "synthetic", "sha_bg", 160 * 4096, 1024.0, 640.0, np.concatenate([o.signal for o in ref_obs]), prep, (0.0, 640.0), "SYN")
    tg_ts = LoadedTimeSeries("tg", "synthetic", "sha_tg", 4 * 4096, 1024.0, 16.0, np.concatenate([o.signal for o in seq_abstain[:4]]), prep, (700.0, 716.0), "SYN")
    base_res = run_workstation_analysis(tg_ts, bg_ts, WorkstationConfig(search_mode="transient", target_fap=0.01))

    # Also construct a UNION_UNAVAILABLE window to verify all three Union statuses round-trip cleanly
    unavail_union_win = replace(
        wins_abstain[60],
        stft_threshold=None,
        broadband_threshold=None,
        detection_available=False,
        diagnostics={
            **wins_abstain[60].diagnostics,
            "dedicated_stft_available": False,
            "dedicated_broadband_available": False,
            "union_stft_available": False,
            "union_broadband_available": False,
            "union_stft_threshold": None,
            "union_mode_status": "UNION_UNAVAILABLE",
        },
    )
    mixed_wins = list(wins_abstain) + [replace(unavail_union_win, index=len(wins_abstain), sample_id="union_unavail_win")]
    cands_m, tl_m = _build_rows(mixed_wins, window_duration_s=4.0, search_mode="transient")
    n_anom_m = sum(1 for r in tl_m if r["active_mode_anomaly_flag"])
    res_export = replace(
        base_res,
        windows=mixed_wins,
        candidates_rows=cands_m,
        timeline_rows=tl_m,
        summary_metrics={**base_res.summary_metrics, "n_active_mode_anomaly_windows": n_anom_m},
    )
    paths = export_workstation_run(res_export, tmp_path / "export_abstention")

    # Read back written CSV and JSON
    with paths["timeline_csv"].open("r", encoding="utf-8") as fh:
        csv_rows = list(csv.DictReader(fh))
    report_data = json.loads(paths["run_report_json"].read_text(encoding="utf-8"))
    with paths["candidates_csv"].open("r", encoding="utf-8") as fh:
        csv_cands = list(csv.DictReader(fh))

    required_union_cols = {
        "dedicated_stft_available",
        "dedicated_broadband_available",
        "union_stft_available",
        "union_broadband_available",
        "union_stft_threshold",
        "union_mode_status",
        "resource_limit_blocked_stft",
        "resource_limit_blocked_union",
        "active_mode_diagnostic_nomination",
        "active_mode_anomaly_flag",
        "calibration_updated_for_next_window",
    }
    assert required_union_cols.issubset(set(csv_rows[0].keys()))
    assert sum(1 for r in csv_rows if r["active_mode_anomaly_flag"] == "True") == report_data["summary_metrics"]["n_active_mode_anomaly_windows"]
    assert csv_rows[60]["detection_available"] == "False"
    assert csv_rows[60]["active_mode_diagnostic_nomination"] == "True"
    assert csv_rows[60]["active_mode_anomaly_flag"] == "False"
    assert any(r["window_index"] == "60" for r in csv_cands)

    # Check FULL_UNION_AVAILABLE, PARTIAL_BROADBAND_ONLY, and UNION_UNAVAILABLE in written JSON & CSV
    json_tl = report_data["window_availability_timeline"]
    statuses_json = {r["union_mode_status"] for r in json_tl}
    statuses_csv = {r["union_mode_status"] for r in csv_rows}
    assert {"FULL_UNION_AVAILABLE", "PARTIAL_BROADBAND_ONLY", "UNION_UNAVAILABLE"}.issubset(statuses_json)
    assert statuses_json == statuses_csv
    assert json_tl[-1]["union_mode_status"] == "UNION_UNAVAILABLE"
    assert json_tl[-1]["union_stft_threshold"] is None
    assert json_tl[-1]["stft_threshold"] is None
    assert json_tl[-1]["broadband_threshold"] is None

    # 2. Persistent mode active decision check (uses broadband_nominated gated by detection_available)
    _, tl_pers = _build_rows(mixed_wins, window_duration_s=4.0, search_mode="persistent")
    for w, r in zip(mixed_wins, tl_pers):
        assert r["active_mode_diagnostic_nomination"] is bool(w.broadband_nominated)
        assert r["active_mode_anomaly_flag"] is bool(w.detection_available and w.broadband_nominated)
    assert tl_pers[-1]["active_mode_anomaly_flag"] is False  # detection_available=False -> False!

    # 3. Calibration-switch window export accuracy (pre-update state on window i, post-update on window i + 1)
    sep_cal_sw = calibrate_separate_search_modes(
        det, bb, ref_obs, target_fap_transient=0.02, target_fap_persistent=0.02,
        target_fwer_union=0.04, calibration_seed=1,
    )
    shifted_sw = _white_obs(140, seed=963, sigma=4.0, prefix="sw")
    mgr_sw = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient", provisional_recal_windows=60, mature_recal_windows=100)
    )
    wins_sw = mgr_sw.process_sequence(
        shifted_sw, det, bb, sep_cal_sw, mgr_sw.fit_continuum_reference(ref_obs), fit_glitch_reference(ref_obs)
    )
    _, tl_sw = _build_rows(wins_sw, window_duration_s=4.0, search_mode="transient")
    sw_idxs = [r["window_index"] for r in tl_sw if r["calibration_updated_for_next_window"]]
    assert len(sw_idxs) >= 2
    s_prov = sw_idxs[0]
    assert tl_sw[s_prov]["regime_id"] == 0
    assert tl_sw[s_prov]["union_mode_status"] == "FULL_UNION_AVAILABLE"
    assert tl_sw[s_prov]["union_stft_threshold"] is not None
    assert tl_sw[s_prov + 1]["regime_id"] == 1
    assert tl_sw[s_prov + 1]["stft_threshold"] is None
    assert tl_sw[s_prov + 1]["union_stft_threshold"] is None
    assert tl_sw[s_prov + 1]["union_mode_status"] == "PARTIAL_BROADBAND_ONLY"
    assert tl_sw[s_prov + 1]["detection_available"] is False
    assert tl_sw[s_prov + 1]["active_mode_anomaly_flag"] is False

    s_rec = sw_idxs[1]
    assert tl_sw[s_rec]["stft_threshold"] is None
    assert tl_sw[s_rec]["union_mode_status"] == "PARTIAL_BROADBAND_ONLY"
    assert tl_sw[s_rec + 1]["stft_threshold"] is not None
    assert tl_sw[s_rec + 1]["detection_available"] is True
