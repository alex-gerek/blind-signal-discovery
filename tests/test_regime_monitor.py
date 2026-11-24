"""Automated pytest suite for BSD v0.3 Phase 1 (BackgroundRegimeMonitor)."""

from __future__ import annotations

import numpy as np
import pytest

from bsd.generator import Observation
from bsd.regime_monitor import (
    STATUS_REGIME_SHIFT,
    STATUS_STATIONARY,
    STATUS_WARMUP,
    BackgroundRegimeMonitor,
    RegimeMonitorConfig,
    summarize_gating_impact,
)


def _make_white_obs(rng: np.random.Generator, n_win: int, sigma: float = 1.0, prefix: str = "w") -> list[Observation]:
    return [
        Observation(f"{prefix}_{i}", rng.normal(0.0, sigma, size=4096), 1024.0, 4.0)
        for i in range(n_win)
    ]


def test_causal_trailing_has_zero_future_dependence() -> None:
    """Modifying future observations j > k must not alter causal_trailing assessment at k."""
    rng = np.random.default_rng(101)
    ref_obs = _make_white_obs(rng, 60, sigma=1.0, prefix="ref")
    seq_a = _make_white_obs(rng, 40, sigma=1.0, prefix="seq")
    # Corrupt future observations indices 25..39 with massive 10x variance jump
    seq_b = list(seq_a[:25]) + [
        Observation(o.sample_id, o.signal * 10.0, o.sample_rate, o.duration) for o in seq_a[25:]
    ]

    mon = BackgroundRegimeMonitor(RegimeMonitorConfig(mode="causal_trailing", local_window_count=10))
    ref = mon.fit_reference(ref_obs)
    out_a = mon.monitor_sequence(seq_a, ref)
    out_b = mon.monitor_sequence(seq_b, ref)

    for idx in range(25):
        assert out_a[idx].future_data_used is False
        assert out_a[idx].min_latency_windows == 1
        assert out_a[idx].regime_shifted == out_b[idx].regime_shifted
        assert out_a[idx].max_feature_z_drift == pytest.approx(out_b[idx].max_feature_z_drift)
        assert out_a[idx].cusum_stat == pytest.approx(out_b[idx].cusum_stat)

    # After the abrupt jump, causal monitor must detect the regime shift
    assert any(a.regime_shifted for a in out_b[28:])


def test_offline_symmetric_documents_future_dependence_and_detects_shift() -> None:
    """offline_symmetric documents future-data use and detects an abrupt shift cleanly."""
    rng = np.random.default_rng(202)
    ref_obs = _make_white_obs(rng, 60, sigma=1.0, prefix="ref")
    seq = _make_white_obs(rng, 25, sigma=1.0, prefix="s1") + _make_white_obs(rng, 25, sigma=2.2, prefix="s2")

    mon = BackgroundRegimeMonitor(RegimeMonitorConfig(mode="offline_symmetric", local_window_count=12))
    ref = mon.fit_reference(ref_obs)
    out = mon.monitor_sequence(seq, ref)

    assert all(a.future_data_used is True for a in out)
    assert all(a.min_latency_windows == 7 for a in out)
    # Early stationary stretch should mostly be STATIONARY_CREDIBLE
    assert sum(a.credibility_status == STATUS_STATIONARY for a in out[:18]) >= 14
    # Late 2.2x variance stretch must be flagged as REGIME_SHIFT_DETECTED
    assert sum(a.credibility_status == STATUS_REGIME_SHIFT for a in out[30:]) >= 18


def test_stationary_noise_preserves_high_duty_cycle() -> None:
    """On stationary in-distribution noise, monitor must maintain high duty cycle."""
    rng = np.random.default_rng(303)
    ref_obs = _make_white_obs(rng, 100, sigma=1.0, prefix="ref")
    test_obs = _make_white_obs(rng, 80, sigma=1.0, prefix="tst")

    mon = BackgroundRegimeMonitor(RegimeMonitorConfig(mode="causal_trailing", local_window_count=12))
    ref = mon.fit_reference(ref_obs)
    out = mon.monitor_sequence(test_obs, ref)

    assert out[0].credibility_status == STATUS_WARMUP
    summary = summarize_gating_impact(out, [False] * 80, [False] * 80)
    assert summary["effective_duty_cycle"] >= 0.90


def test_summarize_gating_impact_no_false_credit() -> None:
    """Excluded background windows must NOT be counted in the surviving FAR denominator."""
    rng = np.random.default_rng(404)
    ref_obs = _make_white_obs(rng, 40, sigma=1.0)
    seq = _make_white_obs(rng, 10, sigma=1.0) + _make_white_obs(rng, 10, sigma=3.0)
    mon = BackgroundRegimeMonitor(RegimeMonitorConfig(mode="offline_symmetric", local_window_count=6))
    out = mon.monitor_sequence(seq, mon.fit_reference(ref_obs))

    # 20 windows: 18 background (1 FP in stationary, 5 FPs in shifted), 2 signals (1 in stationary, 1 in shifted)
    is_sig = [False] * 20
    is_sig[3] = True
    is_sig[17] = True
    detected = [False] * 20
    detected[2] = True   # stationary BG FP
    detected[3] = True   # stationary signal TP
    detected[15] = True  # shifted BG FP
    detected[17] = True  # shifted signal TP (lost to gating)

    s = summarize_gating_impact(out, detected, is_sig)
    assert s["n_excluded_windows"] > 0
    assert s["bg_n_surviving"] < s["bg_n_total"]
    # Post-gating FAR denominator MUST be surviving background windows only, never bg_n_total
    assert s["bg_far_after_gating_on_surviving"] == pytest.approx(
        s["bg_fp_after_gating"] / s["bg_n_surviving"]
    )
    assert s["signals_lost_to_gating"] == 1
    assert s["efficiency_before_gating"] == pytest.approx(1.0)
    assert s["efficiency_after_gating"] == pytest.approx(0.5)
