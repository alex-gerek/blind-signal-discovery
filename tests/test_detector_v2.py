"""Automated pytest suite for BSD v0.2 Phase 1 (MultiResolutionDetector)."""

from __future__ import annotations

import numpy as np
import pytest

from bsd.detector_v2 import (
    DEFAULT_V02_SCALES,
    MultiResolutionCalibration,
    MultiResolutionDetector,
)
from bsd.generator import (
    ALL_SIGNAL_FAMILIES,
    generate_noise_dataset,
    generate_sample,
    make_split_rng,
)


def test_default_scales_match_architecture_spec() -> None:
    """Default v0.2 scales must be 31.25 ms, 125 ms, and 500 ms at 1024 Hz."""
    fs = 1024.0
    durations_ms = [round(1000.0 * s.detector_config.nperseg / fs, 2) for s in DEFAULT_V02_SCALES]
    assert durations_ms == [31.25, 125.0, 500.0]


def test_input_isolation_and_calibration_validation() -> None:
    """Detector strictly accepts Observation objects and valid calibration."""
    det = MultiResolutionDetector()
    cal_split = generate_noise_dataset(
        "calibration", 80, 1024.0, 4.0, ("gaussian_white", "pink_colored"), master_seed=19
    )
    cal = det.calibrate(cal_split.observations, target_fap=0.05, norm_fraction=0.35, calibration_seed=19)
    assert isinstance(cal, MultiResolutionCalibration)
    assert cal.n_norm_samples + cal.n_threshold_samples == 80
    assert cal.empirical_calibration_fap <= 0.05 + 1e-9
    assert len(cal.scale_stats) == 3
    assert sum(cal.scale_win_fractions_in_null.values()) == pytest.approx(1.0)

    rng = make_split_rng(19, "development")
    obs, gt = generate_sample("iso_v2", 1024.0, 4.0, "gaussian_white", rng, "sinusoid", 15.0)
    with pytest.raises(TypeError):
        det.analyze(gt, cal)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        det.analyze(obs, None)  # type: ignore[arg-type]


def test_score_consistency_and_reproducibility() -> None:
    """Standardized scores, global maximum, and repeated runs must be exact and consistent."""
    det = MultiResolutionDetector()
    cal_split = generate_noise_dataset(
        "calibration", 80, 1024.0, 4.0, ("gaussian_white", "pink_colored", "nonstationary"), master_seed=23
    )
    cal1 = det.calibrate(cal_split.observations, target_fap=0.05, calibration_seed=23)
    cal2 = det.calibrate(cal_split.observations, target_fap=0.05, calibration_seed=23)
    assert cal1.global_threshold == pytest.approx(cal2.global_threshold)

    rng = make_split_rng(23, "test_signal")
    obs, _gt = generate_sample("cons_0", 1024.0, 4.0, "pink_colored", rng, "linear_chirp", 18.0)
    r1 = det.analyze(obs, cal1)
    r2 = det.analyze(obs, cal2)

    assert r1.max_global_score == pytest.approx(r2.max_global_score)
    assert r1.winning_scale == r2.winning_scale
    assert r1.max_global_score == pytest.approx(max(r1.per_scale_std_max.values()))
    assert r1.winning_scale == max(r1.per_scale_std_max, key=r1.per_scale_std_max.get)  # type: ignore[arg-type]
    if r1.detected:
        assert r1.candidates[0].standardized_score > cal1.global_threshold
        assert r1.candidates[0].standardized_score == pytest.approx(r1.max_global_score)


@pytest.mark.parametrize("family", ALL_SIGNAL_FAMILIES)
def test_multiresolution_detects_all_signal_families(family: str) -> None:
    """All known and withheld signal families at strong SNR=25 are detected blindly."""
    det = MultiResolutionDetector()
    cal_split = generate_noise_dataset(
        "calibration", 90, 1024.0, 4.0, ("gaussian_white", "pink_colored", "nonstationary"), master_seed=31
    )
    cal = det.calibrate(cal_split.observations, target_fap=0.02, calibration_seed=31)

    rng = make_split_rng(31, f"v2_strong_{family}")
    obs, _gt = generate_sample(f"v2_{family}", 1024.0, 4.0, "pink_colored", rng, family, 25.0)
    res = det.analyze(obs, cal)
    assert res.detected is True
    assert len(res.candidates) >= 1
    assert res.candidates[0].empirical_p_value is not None
