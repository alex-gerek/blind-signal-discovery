"""GWOSC HDF5 loader and physical data-adaptation layer for LIGO strain data.

This module intentionally lives outside the frozen BSD v0.1 detector
(`src/bsd/detector.py`) so the baseline detection algorithm remains completely
unchanged. Two clearly separated preprocessing modes are provided:

1. ``frozen_v01_minimal``:
   Minimal numerical scaling (estimated from independent background) plus
   anti-aliased decimation to 1024 Hz so the frozen BSD v0.1 detector can run
   without floating-point underflow on ~1e-19 LIGO strain values.

2. ``adapted_ligo_whitened``:
   Domain-specific physical adaptation for real LIGO interferometers:
   - 1e18 numerical scaling into O(1) float64 units
   - Frequency-domain whitening using a robust median-Welch PSD estimated
     strictly from independent background-only data
   - Zero-phase 4th-order Butterworth bandpass filter (30-380 Hz)
   - Anti-aliased decimation to 1024 Hz
   All filtering/whitening uses padded context windows when available so the
   4-second evaluation windows remain free of edge transients.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import signal as sp_signal
from scipy.signal.windows import tukey

from bsd.detector import CandidateEvent
from bsd.generator import Observation

STRAIN_NUMERICAL_SCALE = 1.0e18


def _decode_scalar(val: Any) -> Any:
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray) and val.ndim == 0:
        return _decode_scalar(val.item())
    if hasattr(val, "item"):
        try:
            return val.item()
        except Exception:
            return str(val)
    return val


@dataclass(frozen=True)
class GWOSCStrainRecord:
    """In-memory representation of an official GWOSC HDF5 strain file."""

    detector: str
    gps_start: float
    gps_end: float
    sample_rate: float
    strain: np.ndarray
    dq_mask: np.ndarray
    dq_shortnames: tuple[str, ...]
    source_file: str
    n_nonfinite: int
    meta: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class BackgroundModel:
    """Noise characterization estimated strictly from independent background data."""

    detector: str
    mode: str
    source_sample_rate: float
    target_sample_rate: float
    scale_factor: float
    psd_freqs: np.ndarray | None
    psd_values: np.ndarray | None
    bandpass_hz: tuple[float, float] | None
    n_background_segments_used: int
    background_gps_span: tuple[float, float]


@dataclass(frozen=True)
class PaddedSegment:
    """A 4-second target segment paired with optional surrounding context samples."""

    gps_start: float
    duration_s: float
    sample_rate: float
    padded_strain: np.ndarray
    pad_left_samples: int
    pad_right_samples: int


def load_gwosc_hdf5(path: Path) -> GWOSCStrainRecord:
    """Load an official GWOSC HDF5 strain dataset and validate metadata/integrity."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"GWOSC HDF5 file not found: {path}")

    with h5py.File(path, "r") as h5:
        if "strain" not in h5 or "Strain" not in h5["strain"]:
            raise ValueError(f"Missing 'strain/Strain' dataset in {path}")
        strain_ds = h5["strain/Strain"]
        strain = np.asarray(strain_ds[:], dtype=np.float64)
        xspacing = float(strain_ds.attrs["Xspacing"])
        xstart = float(strain_ds.attrs["Xstart"])
        if xspacing <= 0.0:
            raise ValueError(f"Invalid Xspacing={xspacing} in {path}")
        fs = 1.0 / xspacing

        dq_mask = np.asarray(h5["quality/simple/DQmask"][:], dtype=np.int64)
        dq_shortnames = tuple(_decode_scalar(s) for s in h5["quality/simple/DQShortnames"][:])
        meta = {k: _decode_scalar(h5["meta"][k][()]) for k in h5["meta"].keys()}
        detector = str(meta.get("Detector", path.name.split("-")[0]))

    n_nonfinite = int(np.sum(~np.isfinite(strain)))
    gps_end = xstart + (strain.size / fs)
    return GWOSCStrainRecord(
        detector=detector,
        gps_start=xstart,
        gps_end=gps_end,
        sample_rate=fs,
        strain=strain,
        dq_mask=dq_mask,
        dq_shortnames=dq_shortnames,
        source_file=path.name,
        n_nonfinite=n_nonfinite,
        meta=meta,
    )


def extract_clean_segments(
    record: GWOSCStrainRecord,
    segment_duration_s: float = 4.0,
    pad_s: float = 1.0,
    exclude_gps_interval: tuple[float, float] | None = None,
    require_dq_data_bit: bool = True,
    step_s: float | None = None,
) -> list[PaddedSegment]:
    """Slice a GWOSC record into finite, DQ-verified segments with optional padding.

    Any segment whose padded interval overlaps ``exclude_gps_interval`` is
    strictly omitted, preventing background/calibration leakage from target data.
    """
    fs = record.sample_rate
    seg_len = int(round(segment_duration_s * fs))
    pad_len = int(round(pad_s * fs))
    step_len = seg_len if step_s is None else int(round(step_s * fs))
    n_total = record.strain.size

    segments: list[PaddedSegment] = []
    start_idx = pad_len
    while start_idx + seg_len + pad_len <= n_total:
        gps0 = record.gps_start + start_idx / fs
        gps1 = gps0 + segment_duration_s
        pad_gps0 = gps0 - pad_s
        pad_gps1 = gps1 + pad_s

        if exclude_gps_interval is not None:
            ex0, ex1 = exclude_gps_interval
            if not (pad_gps1 <= ex0 or pad_gps0 >= ex1):
                start_idx += step_len
                continue

        # Data-quality check (1 Hz DQmask aligned to gps_start)
        if require_dq_data_bit and record.dq_mask.size > 0:
            sec0 = int(np.floor(pad_gps0 - record.gps_start))
            sec1 = int(np.ceil(pad_gps1 - record.gps_start))
            sec0 = max(0, sec0)
            sec1 = min(record.dq_mask.size, sec1)
            if sec1 <= sec0 or not np.all((record.dq_mask[sec0:sec1] & 1) > 0):
                start_idx += step_len
                continue

        chunk = record.strain[start_idx - pad_len : start_idx + seg_len + pad_len]
        if not np.all(np.isfinite(chunk)):
            start_idx += step_len
            continue

        segments.append(
            PaddedSegment(
                gps_start=gps0,
                duration_s=segment_duration_s,
                sample_rate=fs,
                padded_strain=chunk.copy(),
                pad_left_samples=pad_len,
                pad_right_samples=pad_len,
            )
        )
        start_idx += step_len

    return segments


def fit_background_model(
    detector: str,
    bg_segments: list[PaddedSegment],
    mode: str = "adapted_ligo_whitened",
    target_sample_rate: float = 1024.0,
    bandpass_hz: tuple[float, float] | None = (30.0, 380.0),
) -> BackgroundModel:
    """Estimate numerical scale and (for adapted mode) median-Welch PSD from background."""
    if not bg_segments:
        raise ValueError("bg_segments must be non-empty.")
    source_fs = bg_segments[0].sample_rate

    # Use unpadded central 4-second chunks for background estimation
    core_chunks = [
        seg.padded_strain[
            seg.pad_left_samples : seg.padded_strain.size - seg.pad_right_samples
        ]
        if seg.pad_right_samples > 0
        else seg.padded_strain[seg.pad_left_samples :]
        for seg in bg_segments
    ]
    concat = np.concatenate(core_chunks)
    gps_span = (
        float(min(s.gps_start for s in bg_segments)),
        float(max(s.gps_start + s.duration_s for s in bg_segments)),
    )

    if mode == "frozen_v01_minimal":
        centered = concat - np.median(concat)
        mad_sigma = float(1.4826 * np.median(np.abs(centered)))
        if mad_sigma <= 0.0:
            raise ValueError("Background MAD sigma is non-positive.")
        return BackgroundModel(
            detector=detector,
            mode=mode,
            source_sample_rate=source_fs,
            target_sample_rate=target_sample_rate,
            scale_factor=mad_sigma,
            psd_freqs=None,
            psd_values=None,
            bandpass_hz=None,
            n_background_segments_used=len(bg_segments),
            background_gps_span=gps_span,
        )

    if mode == "adapted_ligo_whitened":
        scaled = (concat - np.median(concat)) * STRAIN_NUMERICAL_SCALE
        nperseg = int(round(4.0 * source_fs))
        freqs, psd = sp_signal.welch(
            scaled,
            fs=source_fs,
            window="hann",
            nperseg=min(nperseg, scaled.size),
            noverlap=min(nperseg, scaled.size) // 2,
            average="median",
        )
        psd = np.maximum(psd, 1e-30)
        return BackgroundModel(
            detector=detector,
            mode=mode,
            source_sample_rate=source_fs,
            target_sample_rate=target_sample_rate,
            scale_factor=STRAIN_NUMERICAL_SCALE,
            psd_freqs=freqs.astype(np.float64),
            psd_values=psd.astype(np.float64),
            bandpass_hz=bandpass_hz,
            n_background_segments_used=len(bg_segments),
            background_gps_span=gps_span,
        )

    raise ValueError(f"Unknown preprocessing mode: {mode!r}")


def _decimate_if_needed(x: np.ndarray, source_fs: float, target_fs: float) -> np.ndarray:
    if np.isclose(source_fs, target_fs):
        return x
    ratio = source_fs / target_fs
    q = int(round(ratio))
    if not np.isclose(ratio, q) or q < 2:
        raise ValueError(f"Unsupported decimation ratio {source_fs} -> {target_fs}")
    return sp_signal.decimate(x, q=q, ftype="iir", zero_phase=True).astype(np.float64)


def preprocess_segment(
    segment: PaddedSegment,
    sample_id: str,
    bg_model: BackgroundModel,
) -> Observation:
    """Convert a raw LIGO strain segment into a BSD Observation using a frozen BackgroundModel."""
    x_raw = np.asarray(segment.padded_strain, dtype=np.float64)
    if not np.all(np.isfinite(x_raw)):
        raise ValueError(f"Non-finite values in segment {sample_id}")

    source_fs = segment.sample_rate
    target_fs = bg_model.target_sample_rate
    q = int(round(source_fs / target_fs))

    if bg_model.mode == "frozen_v01_minimal":
        x = (x_raw - np.median(x_raw)) / bg_model.scale_factor
        x_dec = _decimate_if_needed(x, source_fs, target_fs)
    elif bg_model.mode == "adapted_ligo_whitened":
        assert bg_model.psd_freqs is not None and bg_model.psd_values is not None
        x = (x_raw - np.median(x_raw)) * bg_model.scale_factor
        n = x.size
        # Gentle Tukey taper over outer padding to suppress FFT edge leakage
        taper_alpha = min(0.25, (2.0 * max(segment.pad_left_samples, 16)) / max(n, 1))
        x_win = x * tukey(n, alpha=taper_alpha)

        dt = 1.0 / source_fs
        freqs = np.fft.rfftfreq(n, d=dt)
        hf = np.fft.rfft(x_win)
        interp_psd = np.interp(freqs, bg_model.psd_freqs, bg_model.psd_values)
        asd = np.sqrt(np.maximum(interp_psd, 1e-30))
        white_hf = hf / (asd / np.sqrt(dt / 2.0))
        white_hf[0] = 0.0
        x_white = np.fft.irfft(white_hf, n=n)

        if bg_model.bandpass_hz is not None:
            f_lo, f_hi = bg_model.bandpass_hz
            sos = sp_signal.butter(4, [f_lo, f_hi], btype="bandpass", fs=source_fs, output="sos")
            x_white = sp_signal.sosfiltfilt(sos, x_white)

        x_dec = _decimate_if_needed(x_white, source_fs, target_fs)
    else:
        raise ValueError(f"Unknown mode: {bg_model.mode!r}")

    # Strip left/right padding after filtering/decimation
    trim_left = segment.pad_left_samples // q
    trim_right = segment.pad_right_samples // q
    if trim_right > 0:
        x_core = x_dec[trim_left:-trim_right]
    else:
        x_core = x_dec[trim_left:]

    expected_len = int(round(segment.duration_s * target_fs))
    x_core = x_core[:expected_len]

    return Observation(
        sample_id=sample_id,
        signal=np.asarray(x_core, dtype=np.float64),
        sample_rate=float(target_fs),
        duration=float(x_core.size / target_fs),
    )


def candidate_to_gps(
    candidate: CandidateEvent,
    segment_gps_start: float,
    detector: str,
    segment_id: str,
) -> dict[str, Any]:
    """Convert a segment-relative CandidateEvent into absolute GPS timestamps."""
    return {
        "detector": detector,
        "segment_id": segment_id,
        "candidate_id": candidate.candidate_id,
        "gps_start": float(segment_gps_start + candidate.t_start),
        "gps_end": float(segment_gps_start + candidate.t_end),
        "gps_peak": float(segment_gps_start + candidate.peak_time),
        "duration_s": float(candidate.t_end - candidate.t_start),
        "f_min_hz": float(candidate.f_min),
        "f_max_hz": float(candidate.f_max),
        "peak_freq_hz": float(candidate.peak_freq),
        "anomaly_score": float(candidate.anomaly_score),
        "peak_excess_power": float(candidate.peak_excess_power),
        "n_bins": int(candidate.n_bins),
        "empirical_p_value": candidate.empirical_p_value,
        "p_value_is_upper_bound": candidate.p_value_is_upper_bound,
        "boundary_adjacent": bool(candidate.diagnostics.get("boundary_adjacent", False)),
    }
