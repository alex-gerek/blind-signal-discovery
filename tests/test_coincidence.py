"""Automated pytest suite for BSD v0.2 Phase 3 (sub-10 ms timing & time-slides)."""

from __future__ import annotations

import numpy as np
import pytest

from bsd.coincidence import (
    estimate_cross_detector_delay,
    match_coincident_candidates,
    run_time_slide_analysis,
)


def _fractional_delay_signal(x: np.ndarray, delay_s: float, fs: float) -> np.ndarray:
    """Shift 1D signal by exact fractional delay (positive => delayed/later) via FFT."""
    n = x.size
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    X = np.fft.rfft(x)
    phase = np.exp(-2j * np.pi * freqs * delay_s)
    return np.fft.irfft(X * phase, n=n)


@pytest.mark.parametrize("true_delay_ms,polarity", [(6.9, -1.0), (-4.2, 1.0), (0.0, -1.0), (8.5, 1.0)])
def test_cross_detector_delay_sub_millisecond_accuracy(true_delay_ms: float, polarity: float) -> None:
    """Band-limited |r(tau)| recovers sub-10 ms delays and opposite polarity accurately."""
    fs = 1024.0
    t = np.arange(512, dtype=np.float64) / fs
    # Short broadband chirp-like transient centered at t=0.25 s
    env = np.exp(-0.5 * ((t - 0.25) / 0.03) ** 2)
    clean_b = env * np.sin(2.0 * np.pi * (80.0 * t + 300.0 * (t - 0.20) ** 2))
    clean_a = polarity * _fractional_delay_signal(clean_b, true_delay_ms / 1000.0, fs)

    rng = np.random.default_rng(100)
    sig_a = clean_a + rng.normal(0.0, 0.05, size=t.size)
    sig_b = clean_b + rng.normal(0.0, 0.05, size=t.size)

    res = estimate_cross_detector_delay(sig_a, sig_b, sample_rate=fs, band_hz=(35.0, 350.0), max_lag_s=0.015)
    assert abs(res.delay_s * 1000.0 - true_delay_ms) < 0.5
    assert res.within_10ms_light_travel is True
    assert res.abs_peak_corr > 0.85
    assert np.sign(res.signed_peak_corr) == polarity


def test_match_coincident_candidates_and_time_slides() -> None:
    """Time-slide analysis computes exact counts, livetime, and min resolvable p-value."""
    ca = [{"detector": "H1", "gps_start": 100.0, "gps_end": 100.2, "gps_peak": 100.10, "anomaly_score": 10.0}]
    cb = [{"detector": "L1", "gps_start": 100.01, "gps_end": 100.21, "gps_peak": 100.12, "anomaly_score": 9.0}]
    pairs = match_coincident_candidates(ca, cb, max_dt_s=0.05)
    assert len(pairs) == 1
    assert pairs[0].stft_dt_a_minus_b_s == pytest.approx(-0.02)

    ts = run_time_slide_analysis(
        ca,
        cb,
        single_epoch_livetime_s=1600.0,
        max_dt_s=0.05,
        n_slides=200,
        active_intervals=[(0.0, 1600.0)],
    )
    assert ts["n_slides"] == 200
    assert ts["total_accidental_coincidences"] == 0
    assert ts["significance_status"] == "AVAILABLE"
    assert ts["min_resolvable_p_value"] == pytest.approx(1.0 / 201.0)
    assert ts["equivalent_time_slide_livetime_hours"] > 0.0
    assert ts["equivalent_time_slide_livetime_hours"] < (200 * 1600.0 / 3600.0)

    # Without active_intervals, statistical significance and FAR are explicitly unavailable
    ts_no_int = run_time_slide_analysis(ca, cb, single_epoch_livetime_s=1600.0, max_dt_s=0.05, n_slides=200)
    assert ts_no_int["coincidence_far_per_hour"] is None
    assert ts_no_int["min_resolvable_p_value"] is None
    assert ts_no_int["technical_slide_count_resolution"] == pytest.approx(1.0 / 201.0)
