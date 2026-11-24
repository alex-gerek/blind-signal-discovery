"""Automated pytest suite for the GWOSC LIGO data adapter (Experiment 002)."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from bsd.detector import CandidateEvent, DetectorConfig, TemplateFreeDetector
from bsd.ligo_adapter import (
    GWOSCStrainRecord,
    candidate_to_gps,
    extract_clean_segments,
    fit_background_model,
    load_gwosc_hdf5,
    preprocess_segment,
)


def _create_mock_gwosc_hdf5(
    path: Path,
    detector: str = "H1",
    gps_start: float = 1126259400.0,
    duration_s: int = 32,
    fs: float = 4096.0,
    inject_nan_at: int | None = None,
    dq_zero_second: int | None = None,
) -> Path:
    """Write a minimal GWOSC-compatible HDF5 file for unit testing."""
    rng = np.random.default_rng(1234)
    n_samples = int(duration_s * fs)
    strain = rng.normal(0.0, 1.0e-19, size=n_samples).astype(np.float64)
    if inject_nan_at is not None:
        strain[inject_nan_at] = np.nan

    dq_mask = np.ones(duration_s, dtype=np.int32)
    if dq_zero_second is not None:
        dq_mask[dq_zero_second] = 0

    with h5py.File(path, "w") as h5:
        grp_s = h5.create_group("strain")
        ds = grp_s.create_dataset("Strain", data=strain)
        ds.attrs["Xstart"] = gps_start
        ds.attrs["Xspacing"] = 1.0 / fs

        grp_q = h5.create_group("quality/simple")
        grp_q.create_dataset("DQmask", data=dq_mask)
        grp_q.create_dataset("DQShortnames", data=np.array([b"DATA"]))

        grp_m = h5.create_group("meta")
        grp_m.create_dataset("Detector", data=detector.encode("utf-8"))
        grp_m.create_dataset("GPSstart", data=gps_start)
        grp_m.create_dataset("Duration", data=duration_s)
    return path


def test_hdf5_loading_and_timestamp_preservation(tmp_path: Path) -> None:
    """HDF5 loader preserves exact GPS timestamps, sample rate, and metadata."""
    h5_path = _create_mock_gwosc_hdf5(tmp_path / "H1_mock.hdf5", detector="H1", gps_start=1126259447.0, duration_s=32)
    rec = load_gwosc_hdf5(h5_path)
    assert isinstance(rec, GWOSCStrainRecord)
    assert rec.detector == "H1"
    assert rec.gps_start == pytest.approx(1126259447.0)
    assert rec.gps_end == pytest.approx(1126259479.0)
    assert rec.sample_rate == pytest.approx(4096.0)
    assert rec.strain.size == 32 * 4096
    assert rec.n_nonfinite == 0


def test_noise_segment_isolation_and_dq_filtering(tmp_path: Path) -> None:
    """Excluded GPS intervals, NaN samples, and DQ=0 seconds are strictly omitted."""
    h5_path = _create_mock_gwosc_hdf5(
        tmp_path / "L1_mock.hdf5",
        detector="L1",
        gps_start=1000.0,
        duration_s=40,
        inject_nan_at=int(30.5 * 4096),
        dq_zero_second=22,
    )
    rec = load_gwosc_hdf5(h5_path)
    assert rec.n_nonfinite == 1

    # Exclude GPS [1008, 1016]
    segs = extract_clean_segments(
        rec,
        segment_duration_s=4.0,
        pad_s=1.0,
        exclude_gps_interval=(1008.0, 1016.0),
        require_dq_data_bit=True,
    )
    assert len(segs) > 0
    for seg in segs:
        pad0 = seg.gps_start - 1.0
        pad1 = seg.gps_start + seg.duration_s + 1.0
        # No overlap with excluded target interval
        assert (pad1 <= 1008.0) or (pad0 >= 1016.0)
        # No overlap with DQ=0 second [1022, 1023]
        assert (pad1 <= 1022.0) or (pad0 >= 1023.0)
        # All samples finite
        assert np.all(np.isfinite(seg.padded_strain))


@pytest.mark.parametrize("mode", ["frozen_v01_minimal", "adapted_ligo_whitened"])
def test_strain_scaling_preprocessing_and_detector_reproducibility(tmp_path: Path, mode: str) -> None:
    """Preprocessing scales ~1e-19 strain cleanly to 1024 Hz and yields reproducible BSD outputs."""
    h5_path = _create_mock_gwosc_hdf5(tmp_path / f"H1_{mode}.hdf5", duration_s=48)
    rec = load_gwosc_hdf5(h5_path)
    segs = extract_clean_segments(rec, segment_duration_s=4.0, pad_s=1.0)
    assert len(segs) >= 6

    bg_model = fit_background_model("H1", segs[:4], mode=mode, target_sample_rate=1024.0)
    obs1 = preprocess_segment(segs[4], "seg_test_a", bg_model)
    obs2 = preprocess_segment(segs[4], "seg_test_a", bg_model)

    assert obs1.sample_rate == pytest.approx(1024.0)
    assert obs1.signal.size == 4096
    assert np.all(np.isfinite(obs1.signal))
    # Scaled strain must be O(1), well above 1e-18 underflow
    assert 1e-3 < float(np.std(obs1.signal)) < 1e4
    assert np.array_equal(obs1.signal, obs2.signal)

    det = TemplateFreeDetector(DetectorConfig())
    r1 = det.analyze(obs1, calibration=None)
    r2 = det.analyze(obs2, calibration=None)
    assert r1.max_anomaly_score == pytest.approx(r2.max_anomaly_score)
    assert len(r1.candidates) == len(r2.candidates)


def test_candidate_timestamp_conversion() -> None:
    """Local candidate times convert accurately to absolute GPS timestamps."""
    cand = CandidateEvent(
        candidate_id=1,
        t_start=1.25,
        t_end=1.50,
        f_min=48.0,
        f_max=224.0,
        peak_time=1.422,
        peak_freq=144.0,
        anomaly_score=180.5,
        peak_excess_power=25.0,
        n_bins=24,
        empirical_p_value=0.002,
        p_value_is_upper_bound=True,
        diagnostics={"boundary_adjacent": False},
    )
    gps_dict = candidate_to_gps(cand, segment_gps_start=1126259461.0, detector="H1", segment_id="seg_01")
    assert gps_dict["gps_start"] == pytest.approx(1126259462.25)
    assert gps_dict["gps_end"] == pytest.approx(1126259462.50)
    assert gps_dict["gps_peak"] == pytest.approx(1126259462.422)
    assert gps_dict["duration_s"] == pytest.approx(0.25)
    assert gps_dict["detector"] == "H1"
