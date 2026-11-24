"""Automated pytest suite for Experiment 003 multi-event benchmark utilities."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXP_DIR = Path(__file__).resolve().parent.parent / "experiments"
if str(EXP_DIR) not in sys.path:
    sys.path.insert(0, str(EXP_DIR))

import run_ligo_experiment_003 as exp003


def test_partition_background_segments_strict_separation() -> None:
    """Background partitioning enforces disjoint sets and rejects insufficient data."""
    assert exp003.partition_background_segments(150) is None
    n_psd, n_cal, n_val = exp003.partition_background_segments(982)  # type: ignore[misc]
    assert (n_psd, n_cal, n_val) == (60, 500, 422)
    assert n_psd + n_cal + n_val == 982

    # Fallback regime with gaps
    n_psd2, n_cal2, n_val2 = exp003.partition_background_segments(400)  # type: ignore[misc]
    assert n_psd2 + n_cal2 + n_val2 == 400
    assert n_psd2 >= 30 and n_cal2 >= 100 and n_val2 >= 50


def test_find_coincidences_and_time_shift_control() -> None:
    """Coincidence finder detects close pairs and rejects unphysically separated ones."""
    c_h1 = [
        {"detector": "H1", "gps_start": 100.0, "gps_end": 100.3, "gps_peak": 100.20, "anomaly_score": 120.0}
    ]
    c_l1_near = [
        {"detector": "L1", "gps_start": 100.05, "gps_end": 100.32, "gps_peak": 100.23, "anomaly_score": 95.0}
    ]
    c_l1_far = [
        {"detector": "L1", "gps_start": 120.0, "gps_end": 120.3, "gps_peak": 120.20, "anomaly_score": 95.0}
    ]
    pairs_near = exp003.find_coincidences(c_h1, c_l1_near, max_dt_s=0.10)
    assert len(pairs_near) == 1
    assert pairs_near[0]["dt_a_minus_b_s"] == pytest.approx(-0.03)
    assert exp003.find_coincidences(c_h1, c_l1_far, max_dt_s=0.10) == []

    ts = exp003.time_shift_coincidence_control(c_h1, c_l1_near, max_dt_s=0.10, n_shifts=40)
    assert ts["n_shifts"] == 40
    assert ts["total_accidental_coincidences"] == 0


def test_evaluate_detector_against_reference_gps_separates_glitches() -> None:
    """Glitches separated by >0.5 s from published GPS must not count as event recovery."""
    published_gps = 1187008882.4
    glitch_only = [
        {
            "detector": "L1",
            "gps_start": 1187008881.2,
            "gps_end": 1187008881.4,
            "gps_peak": 1187008881.3,
            "anomaly_score": 500.0,
            "empirical_p_value": 0.002,
            "p_value_is_upper_bound": True,
        }
    ]
    res_glitch = exp003.evaluate_detector_against_reference_gps(glitch_only, published_gps, tolerance_s=0.5)
    assert res_glitch["recovered"] is False
    assert res_glitch["n_other_candidates_in_32s"] == 1

    matched_plus_glitch = glitch_only + [
        {
            "detector": "L1",
            "gps_start": 1187008882.2,
            "gps_end": 1187008882.5,
            "gps_peak": 1187008882.38,
            "anomaly_score": 110.0,
            "empirical_p_value": 0.004,
            "p_value_is_upper_bound": False,
        }
    ]
    res_ok = exp003.evaluate_detector_against_reference_gps(matched_plus_glitch, published_gps, tolerance_s=0.5)
    assert res_ok["recovered"] is True
    assert res_ok["n_matched_candidates"] == 1
    assert res_ok["n_other_candidates_in_32s"] == 1
    assert res_ok["best_matched_score"] == pytest.approx(110.0)
