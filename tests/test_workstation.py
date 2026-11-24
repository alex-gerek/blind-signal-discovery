"""Automated pytest suite for the BSD Local Measurement Workstation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.io import wavfile

from bsd.adaptive_regime_v2b import CalibrationUncertaintySummary
from bsd.io_loader import (
    InputValidationError,
    O3HoldoutBlockedError,
    assert_not_o3_holdout,
    check_background_independence,
    load_timeseries_file,
)
from bsd.workstation_runner import (
    CAL_STATUS_INSUFFICIENT,
    CAL_STATUS_PROVISIONAL,
    CAL_STATUS_VALIDATED,
    WorkstationConfig,
    classify_calibration_status,
    export_workstation_run,
    run_synthetic_self_test,
    run_workstation_analysis,
)


def test_io_loader_formats_validation_and_resampling(tmp_path: Path) -> None:
    """Validate .wav, .csv, .npy loading, explicit rejections, and opt-in resampling."""
    rng = np.random.default_rng(11)

    # Unsupported extension rejected
    bad_ext = tmp_path / "data.txt"
    bad_ext.write_text("1.0\n2.0\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Unsupported file format"):
        load_timeseries_file(bad_ext, sample_rate_hz=1024.0)

    # Sample rate < 1024 Hz rejected
    low_fs = tmp_path / "low.npy"
    np.save(low_fs, rng.normal(0, 1, 4096))
    with pytest.raises(InputValidationError, match="Unsupported sample rate"):
        load_timeseries_file(low_fs, sample_rate_hz=512.0)

    # NaN values rejected
    nan_arr = rng.normal(0, 1, 4096)
    nan_arr[10] = np.nan
    nan_file = tmp_path / "nan.npy"
    np.save(nan_file, nan_arr)
    with pytest.raises(InputValidationError, match="non-finite"):
        load_timeseries_file(nan_file, sample_rate_hz=1024.0)

    # Stereo WAV requires explicit channel_index; 2048 Hz requires allow_resample=True
    stereo = (rng.normal(0, 0.2, (8192, 2)) * 32767).astype(np.int16)
    wav_path = tmp_path / "stereo_2048.wav"
    wavfile.write(str(wav_path), 2048, stereo)
    with pytest.raises(InputValidationError, match="channels"):
        load_timeseries_file(wav_path, allow_resample=True)
    with pytest.raises(InputValidationError, match="allow_resample|resampling"):
        load_timeseries_file(wav_path, channel_index=0, allow_resample=False)

    ts_wav = load_timeseries_file(wav_path, channel_index=1, allow_resample=True)
    assert ts_wav.sample_rate_hz == pytest.approx(1024.0)
    assert ts_wav.preprocessing.resampled is True
    assert ts_wav.preprocessing.original_sample_rate_hz == pytest.approx(2048.0)
    assert ts_wav.preprocessing.channel_selected == "channel_1"
    assert ts_wav.preprocessing.ligo_whitening_applied is False


def test_o3_holdout_hard_blocked(tmp_path: Path) -> None:
    """Reserved O3 filenames and O3 GPS intervals must be hard-blocked."""
    fake_o3 = tmp_path / "H-H1_GWOSC_4KHZ_GW190521_074359-4096.npy"
    np.save(fake_o3, np.ones(4096))
    with pytest.raises(O3HoldoutBlockedError):
        load_timeseries_file(fake_o3, sample_rate_hz=1024.0)
    with pytest.raises(O3HoldoutBlockedError):
        assert_not_o3_holdout("arbitrary.hdf5", gps_interval=(1242459000.0, 1242460000.0))


def test_background_independence_rejects_same_sha_and_overlap(tmp_path: Path) -> None:
    """Identical SHA-256 and overlapping timestamps must be rejected; unverified files warn."""
    rng = np.random.default_rng(22)
    f1 = tmp_path / "sig1.npy"
    f2 = tmp_path / "sig2.npy"
    np.save(f1, rng.normal(0, 1, 4096))
    np.save(f2, rng.normal(0, 1, 4096))

    ts1 = load_timeseries_file(f1, sample_rate_hz=1024.0)
    ts1_dup = load_timeseries_file(f1, sample_rate_hz=1024.0)
    ts2 = load_timeseries_file(f2, sample_rate_hz=1024.0)

    with pytest.raises(InputValidationError, match="identical SHA-256"):
        check_background_independence(ts1, ts1_dup)

    audit = check_background_independence(ts1, ts2)
    assert audit.status == "UNVERIFIED_EXTERNAL_CHECK_REQUIRED"
    assert len(audit.warnings) >= 1

    # Overlapping CSV time columns rejected
    t_a = np.arange(4096) / 1024.0
    t_b = t_a + 2.0  # overlaps [2.0, 4.0]
    c1 = tmp_path / "a.csv"
    c2 = tmp_path / "b.csv"
    c1.write_text("time,val\n" + "\n".join(f"{t:.6f},{v:.5f}" for t, v in zip(t_a, rng.normal(0, 1, 4096))), encoding="utf-8")
    c2.write_text("time,val\n" + "\n".join(f"{t:.6f},{v:.5f}" for t, v in zip(t_b, rng.normal(0, 1, 4096))), encoding="utf-8")
    tsc1 = load_timeseries_file(c1, csv_value_column="val", csv_time_column="time")
    tsc2 = load_timeseries_file(c2, csv_value_column="val", csv_time_column="time")
    with pytest.raises(InputValidationError, match="overlaps"):
        check_background_independence(tsc1, tsc2)


def test_calibration_status_three_tiers() -> None:
    """VALIDATED must never be assigned solely from window count when statistical criteria fail."""
    u_good = CalibrationUncertaintySummary(200, 185.0, 0.02, 0.72, 0.005, 0.01, (0.002, 0.021), (0.002, 0.022), False)
    assert classify_calibration_status(260, u_good).status == CAL_STATUS_VALIDATED

    # Large window count (260) but temporal dependence detected -> MUST be PROVISIONAL, never VALIDATED!
    u_dep = CalibrationUncertaintySummary(200, 130.0, 0.35, 0.004, 0.005, 0.01, (0.002, 0.021), (0.002, 0.024), True)
    audit_dep = classify_calibration_status(260, u_dep)
    assert audit_dep.status == CAL_STATUS_PROVISIONAL
    assert any("Ljung-Box" in r for r in audit_dep.unmet_validation_conditions)

    # Too few effective samples -> INSUFFICIENT
    u_small = CalibrationUncertaintySummary(60, 30.0, 0.45, 0.001, 0.016, 0.016, (0.001, 0.05), (0.001, 0.09), True)
    assert classify_calibration_status(80, u_small).status == CAL_STATUS_INSUFFICIENT


def test_end_to_end_local_wav_workflow_and_export(tmp_path: Path) -> None:
    """End-to-end local .wav ingestion, calibration, detection, and reproducible CSV/JSON export."""
    rng = np.random.default_rng(33)
    fs = 1024
    # 100 background windows (400 s -> n_thresh=70 -> PROVISIONAL) + 12 target windows with one strong burst in win 5
    bg = rng.normal(0.0, 0.1, size=100 * 4 * fs)
    tgt = rng.normal(0.0, 0.1, size=12 * 4 * fs)
    t_win = np.arange(4 * fs, dtype=np.float64) / fs
    burst = 0.8 * np.exp(-0.5 * ((t_win - 2.0) / 0.01) ** 2) * np.cos(2.0 * np.pi * 190.0 * (t_win - 2.0))
    tgt[5 * 4 * fs : 6 * 4 * fs] += burst

    bg_wav = tmp_path / "bg.wav"
    tgt_wav = tmp_path / "tgt.wav"
    wavfile.write(str(bg_wav), fs, bg.astype(np.float32))
    wavfile.write(str(tgt_wav), fs, tgt.astype(np.float32))

    bg_ts = load_timeseries_file(bg_wav)
    tgt_ts = load_timeseries_file(tgt_wav)
    res = run_workstation_analysis(tgt_ts, bg_ts, WorkstationConfig(search_mode="transient", calibration_seed=33))

    assert res.calibration_audit.status == CAL_STATUS_PROVISIONAL  # 70 windows < 150 n_thresh
    assert res.windows[5].stft_nominated_raw is True
    assert "NOT_ESTABLISHED_BY_SINGLE_CHANNEL_WORKSTATION" in res.scientific_semantics["independent_physical_confirmation"]

    exported = export_workstation_run(res, tmp_path / "out")
    assert exported["candidates_csv"].is_file()
    assert exported["timeline_csv"].is_file()
    payload = json.loads(exported["run_report_json"].read_text(encoding="utf-8"))
    assert payload["run_metadata"]["o3_holdout_accessed"] is False
    assert payload["run_metadata"]["target_preprocessing"]["ligo_whitening_applied"] is False


def test_synthetic_self_test_positive_and_negative_controls() -> None:
    """Synthetic self-test reports negative-control FAR and both detected & missed injections."""
    res = run_synthetic_self_test(
        WorkstationConfig(search_mode="transient", processing_mode="causal_trailing"),
        n_bg_windows=80,
        n_neg_control_windows=20,
        n_pos_injections=4,
        noise_model="pink_colored",
        signal_family="linear_chirp",
        target_snr=18.0,
        seed=44,
    )
    assert res.ground_truth_audit is not None
    gt = res.ground_truth_audit
    assert gt["negative_control_windows"] == 20
    assert gt["positive_injections_total"] == 4
    assert gt["positive_injections_detected"] + gt["positive_injections_missed"] == 4
    assert gt["positive_injections_detected"] >= 3


def test_gui_measurement_clarity_labels() -> None:
    """Verify all six GUI measurement clarity labels on generated figures and app source."""
    import importlib.util
    import matplotlib.pyplot as plt

    app_path = Path(__file__).resolve().parent.parent / "app" / "bsd_workstation.py"
    spec = importlib.util.spec_from_file_location("bsd_workstation_app", app_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    res = run_synthetic_self_test(
        WorkstationConfig(search_mode="transient", processing_mode="causal_trailing"),
        n_bg_windows=80,
        n_neg_control_windows=10,
        n_pos_injections=2,
        seed=55,
    )

    # 1, 2, 6: Waveform units, relative dB spectrogram colorbar, and absolute time interval
    fig_w = mod._plot_window_waveform_and_spec(res, win_idx=5)
    ax_wave, ax_spec, ax_cbar = fig_w.axes[0], fig_w.axes[1], fig_w.axes[2]
    assert "normalized_sigma" in ax_wave.get_ylabel()
    assert "Absolute Time: [20.0, 24.0] s" in ax_wave.get_title()
    assert "Relative STFT Power" in ax_cbar.get_ylabel()
    assert "dB re 1 unit^2" in ax_cbar.get_ylabel()
    plt.close(fig_w)

    # 3: Timeline dimensionless Y-axis and numeric threshold in legend
    fig_t = mod._plot_timeline(res)
    ax_t = fig_t.axes[0]
    assert "dimensionless" in ax_t.get_ylabel()
    legend_texts = [t.get_text() for t in ax_t.get_legend().get_texts()]
    assert any("Active Calibrated Threshold =" in s and "dimensionless" in s for s in legend_texts)
    plt.close(fig_t)

    # 4 & 5: Candidate decision-rule caption and linear (NOT dB) SNR slider label
    src_text = app_path.read_text(encoding="utf-8")
    assert "standardized_score > active_threshold" in src_text
    assert "raw_anomaly_score" in src_text
    assert "linear dimensionless" in src_text and "NOT dB" in src_text
