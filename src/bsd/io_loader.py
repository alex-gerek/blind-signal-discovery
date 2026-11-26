"""Domain-agnostic, strictly validated time-series ingestion for BSD Workstation.

Enforces:
  1. Explicit format validation (.wav, .csv, .npy, optional GWOSC .hdf5).
  2. Hard blocklist preventing loading, previewing, or calibrating reserved O3
     holdout data (`docs/v03_o3_holdout_manifest.json` event IDs and O3 GPS span).
  3. Complete, auditable preprocessing provenance (format, channel, units,
     original sample rate, anti-alias resampling, effective bandwidth).
  4. Never silently applies LIGO-specific whitening/bandpass to arbitrary files.
  5. Explicit background/target independence checks (SHA-256 identity and
     timestamp/GPS overlap rejection).
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy import signal as sp_signal
from scipy.io import wavfile

from bsd.generator import Observation
from bsd.ligo_adapter import load_gwosc_hdf5

CANONICAL_SAMPLE_RATE_HZ = 1024.0
MIN_SUPPORTED_SAMPLE_RATE_HZ = 1024.0
SUPPORTED_EXTENSIONS = (".wav", ".csv", ".npy", ".hdf5", ".h5")

# Official LIGO/Virgo Observing Run O3 GPS span (2019-04-01 to 2020-03-27 + guard)
O3_GPS_GUARD_MIN = 1238000000.0
O3_GPS_GUARD_MAX = 1270000000.0


class InputValidationError(ValueError):
    """Raised when an input file or measurement configuration is unsupported or invalid."""


class O3HoldoutBlockedError(PermissionError):
    """Raised when an input matches the reserved untouched O3 holdout manifest."""


@dataclass(frozen=True)
class PreprocessingSummary:
    """Explicit record of every data transformation applied during ingestion."""

    file_format: str
    channel_selected: str
    units: str
    original_sample_rate_hz: float
    effective_sample_rate_hz: float
    resampled: bool
    antialias_filter: str
    effective_bandwidth_hz: tuple[float, float]
    analysis_bands_hz: tuple[float, float]
    ligo_whitening_applied: bool
    steps: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class LoadedTimeSeries:
    """Validated, strictly finite 1D time series ready for BSD window slicing."""

    source_name: str
    file_format: str
    sha256: str
    n_samples: int
    sample_rate_hz: float
    duration_s: float
    signal: np.ndarray
    preprocessing: PreprocessingSummary
    time_interval: tuple[float, float] | None = None
    time_reference: str | None = None  # e.g. "GPS", "CSV_TIME_COLUMN", None
    provenance: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class IndependenceAudit:
    """Result of background/target independence verification."""

    passed: bool
    status: str  # "VERIFIED_DISJOINT_TIMESTAMPS" | "UNVERIFIED_EXTERNAL_CHECK_REQUIRED"
    target_sha256: str
    background_sha256: str
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _load_o3_event_ids() -> tuple[str, ...]:
    manifest_path = Path(__file__).resolve().parent.parent.parent / "docs" / "v03_o3_holdout_manifest.json"
    if not manifest_path.is_file():
        return ()
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    return tuple(str(ev["event_id"]) for ev in data.get("events", []))


def assert_not_o3_holdout(source_name: str, gps_interval: tuple[float, float] | None = None) -> None:
    """Block any reserved O3 holdout filename or O3 GPS interval."""
    name_upper = Path(source_name).name.upper()
    for ev_id in _load_o3_event_ids():
        if ev_id.upper() in name_upper:
            raise O3HoldoutBlockedError(
                f"Access blocked: '{source_name}' matches reserved O3 holdout event '{ev_id}' "
                "in docs/v03_o3_holdout_manifest.json. Reserved O3 data must remain untouched."
            )
    if gps_interval is not None:
        g0, g1 = float(gps_interval[0]), float(gps_interval[1])
        if not (g1 < O3_GPS_GUARD_MIN or g0 > O3_GPS_GUARD_MAX):
            raise O3HoldoutBlockedError(
                f"Access blocked: GPS interval [{g0:.1f}, {g1:.1f}] overlaps Observing Run O3 "
                f"[{O3_GPS_GUARD_MIN:.0f}, {O3_GPS_GUARD_MAX:.0f}]. Reserved O3 data must remain untouched."
            )


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _validate_and_prepare_signal(
    raw_signal: np.ndarray,
    original_fs: float,
    allow_resample: bool,
    source_name: str,
    file_format: str,
    channel_selected: str,
    units: str,
    sha256_hex: str,
    time_interval: tuple[float, float] | None = None,
    time_reference: str | None = None,
    provenance: dict[str, Any] | None = None,
) -> LoadedTimeSeries:
    raw = np.asarray(raw_signal)
    if raw.dtype.kind not in "biuf":
        raise InputValidationError(f"Input '{source_name}' must contain real numeric samples; got {raw.dtype}.")
    x = raw.astype(np.float64)
    if x.ndim != 1:
        raise InputValidationError(
            f"Expected 1D time-series after channel selection for '{source_name}', got shape {x.shape}."
        )
    if x.size == 0:
        raise InputValidationError(f"Input '{source_name}' is empty (0 samples).")
    n_nonfinite = int(np.sum(~np.isfinite(x)))
    if n_nonfinite > 0:
        raise InputValidationError(
            f"Input '{source_name}' contains {n_nonfinite} non-finite (NaN/Inf) samples. "
            "Clean or remove gaps before loading."
        )
    if float(np.std(x)) <= 1e-30:
        raise InputValidationError(
            f"Input '{source_name}' is constant (zero variance) and cannot be analyzed."
        )
    if original_fs <= 0.0 or not np.isfinite(original_fs):
        raise InputValidationError(f"Invalid sample rate {original_fs} Hz for '{source_name}'.")
    if original_fs < MIN_SUPPORTED_SAMPLE_RATE_HZ - 1e-6:
        raise InputValidationError(
            f"Unsupported sample rate {original_fs:.2f} Hz for '{source_name}'. "
            f"Existing BSD multi-resolution STFT and continuum bands (20-450 Hz) require "
            f"sample_rate >= {MIN_SUPPORTED_SAMPLE_RATE_HZ:.0f} Hz (Nyquist >= 512 Hz)."
        )

    steps: list[str] = [f"Loaded {file_format} ({x.size} samples at {original_fs:.2f} Hz, channel={channel_selected})"]
    warnings_list: list[str] = []
    resampled = False
    aa_filter = "none"
    eff_fs = float(original_fs)

    if abs(original_fs - CANONICAL_SAMPLE_RATE_HZ) > 1e-3:
        if not allow_resample:
            raise InputValidationError(
                f"Input '{source_name}' has sample rate {original_fs:.2f} Hz != canonical "
                f"{CANONICAL_SAMPLE_RATE_HZ:.0f} Hz. Enable explicit anti-aliased resampling "
                f"to {CANONICAL_SAMPLE_RATE_HZ:.0f} Hz (Nyquist 512 Hz) to proceed."
            )
        frac = Fraction(CANONICAL_SAMPLE_RATE_HZ / original_fs).limit_denominator(1024)
        x = sp_signal.resample_poly(x, up=frac.numerator, down=frac.denominator).astype(np.float64)
        eff_fs = CANONICAL_SAMPLE_RATE_HZ
        resampled = True
        aa_filter = f"scipy.signal.resample_poly(up={frac.numerator}, down={frac.denominator}, Kaiser FIR anti-alias)"
        steps.append(f"Anti-aliased resampling from {original_fs:.2f} Hz -> {eff_fs:.1f} Hz ({aa_filter})")
        warnings_list.append(
            f"Signal resampled from {original_fs:.2f} Hz to {eff_fs:.1f} Hz with anti-alias FIR filtering; "
            "frequencies above 512 Hz are removed."
        )

    prep = PreprocessingSummary(
        file_format=file_format,
        channel_selected=channel_selected,
        units=units,
        original_sample_rate_hz=float(original_fs),
        effective_sample_rate_hz=eff_fs,
        resampled=resampled,
        antialias_filter=aa_filter,
        effective_bandwidth_hz=(0.0, eff_fs / 2.0),
        analysis_bands_hz=(20.0, 450.0),
        ligo_whitening_applied=False,
        steps=tuple(steps),
    )
    return LoadedTimeSeries(
        source_name=source_name,
        file_format=file_format,
        sha256=sha256_hex,
        n_samples=int(x.size),
        sample_rate_hz=eff_fs,
        duration_s=float(x.size / eff_fs),
        signal=x,
        preprocessing=prep,
        time_interval=time_interval,
        time_reference=time_reference,
        provenance=provenance or {},
        warnings=tuple(warnings_list),
    )


def load_timeseries_file(
    path: Path | str,
    *,
    sample_rate_hz: float | None = None,
    channel_index: int | None = None,
    csv_value_column: str | int = 0,
    csv_time_column: str | int | None = None,
    allow_resample: bool = False,
    units: str = "arbitrary",
) -> LoadedTimeSeries:
    """Load and validate a local .wav, .csv, .npy, or GWOSC .hdf5 file."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {p}")

    assert_not_o3_holdout(p.name)
    ext = p.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise InputValidationError(
            f"Unsupported file format '{ext}' for '{p.name}'. "
            f"Supported formats: {', '.join(SUPPORTED_EXTENSIONS)}."
        )
    sha = _sha256_file(p)

    if ext == ".wav":
        try:
            fs_wav, data = wavfile.read(str(p))
        except Exception as exc:
            raise InputValidationError(f"Failed to read WAV file '{p.name}': {exc}") from exc
        arr = np.asarray(data)
        if arr.ndim == 2:
            if channel_index is None:
                raise InputValidationError(
                    f"WAV file '{p.name}' has {arr.shape[1]} channels. "
                    "Explicitly select a single channel_index (0-based)."
                )
            if not (0 <= channel_index < arr.shape[1]):
                raise InputValidationError(
                    f"channel_index={channel_index} out of range for {arr.shape[1]}-channel WAV '{p.name}'."
                )
            ch_label = f"channel_{channel_index}"
            arr = arr[:, channel_index]
        else:
            ch_label = "mono"
        # Convert integer PCM to [-1, 1] float64 explicitly
        if np.issubdtype(arr.dtype, np.integer):
            info = np.iinfo(arr.dtype)
            scale = float(max(abs(info.min), info.max))
            sig = arr.astype(np.float64) / scale
            wav_units = f"PCM_{arr.dtype}_normalized"
        else:
            sig = arr.astype(np.float64)
            wav_units = f"WAV_{arr.dtype}"
        return _validate_and_prepare_signal(
            sig, float(fs_wav), allow_resample, p.name, "wav", ch_label, wav_units, sha
        )

    if ext == ".npy":
        if sample_rate_hz is None:
            raise InputValidationError(
                f"NumPy file '{p.name}' contains no sample-rate metadata. "
                "Supply sample_rate_hz explicitly."
            )
        try:
            arr = np.load(str(p), allow_pickle=False)
        except (OSError, ValueError, EOFError) as exc:
            raise InputValidationError(f"Failed to read safe numeric NPY '{p.name}': {exc}") from exc
        if arr.ndim == 2:
            if channel_index is None:
                raise InputValidationError(
                    f"NumPy file '{p.name}' is 2D with shape {arr.shape}. "
                    "Explicitly select a single channel_index (column)."
                )
            if not (0 <= channel_index < arr.shape[1]):
                raise InputValidationError(f"channel_index={channel_index} out of range for shape {arr.shape}.")
            ch_label = f"col_{channel_index}"
            sig = arr[:, channel_index]
        elif arr.ndim == 1:
            ch_label = "1d_array"
            sig = arr
        else:
            raise InputValidationError(f"Unsupported NumPy array ndim={arr.ndim} in '{p.name}' (expected 1D or 2D).")
        return _validate_and_prepare_signal(
            sig, float(sample_rate_hz), allow_resample, p.name, "npy", ch_label, units, sha
        )

    if ext == ".csv":
        try:
            with p.open("r", encoding="utf-8-sig", newline="") as fh:
                reader = csv.reader(fh, strict=True)
                numbered_rows = [(reader.line_num, row) for row in reader]
        except (OSError, UnicodeError, csv.Error) as exc:
            raise InputValidationError(f"Failed to read CSV '{p.name}': {exc}") from exc
        if not numbered_rows or not numbered_rows[0][1] or not any(cell.strip() for cell in numbered_rows[0][1]):
            raise InputValidationError(f"CSV file '{p.name}' is empty or starts with a blank row.")
        first = numbered_rows[0][1]
        try:
            float(first[0])
            has_header = False
        except ValueError:
            has_header = True
        header = [h.strip() for h in first] if has_header else []
        data_rows = numbered_rows[1:] if has_header else numbered_rows
        if len(data_rows) < 2:
            raise InputValidationError(f"CSV file '{p.name}' has too few data rows ({len(data_rows)}).")

        def _col_idx(spec: str | int) -> int:
            if isinstance(spec, int) and not isinstance(spec, bool) and spec >= 0:
                return spec
            if isinstance(spec, str) and has_header and header.count(spec) == 1:
                return header.index(spec)
            raise InputValidationError(f"CSV column {spec!r} is missing, ambiguous or invalid in '{p.name}'.")

        width = len(first)
        for line_number, row in data_rows:
            if len(row) != width:
                raise InputValidationError(
                    f"CSV row {line_number} in '{p.name}' has {len(row)} columns; expected {width}."
                )

        def _values(index: int) -> np.ndarray:
            if index >= width:
                raise InputValidationError(f"CSV column {index} out of range for {width} columns in '{p.name}'.")
            values: list[float] = []
            for line_number, row in data_rows:
                try:
                    value = float(row[index])
                except ValueError as exc:
                    raise InputValidationError(
                        f"CSV row {line_number}, column {index} in '{p.name}' is not numeric."
                    ) from exc
                if not np.isfinite(value):
                    raise InputValidationError(
                        f"CSV row {line_number}, column {index} in '{p.name}' contains a non-finite sample."
                    )
                values.append(value)
            return np.asarray(values, dtype=np.float64)

        v_idx = _col_idx(csv_value_column)
        vals = _values(v_idx)
        t_interval: tuple[float, float] | None = None
        t_ref: str | None = None
        fs_eff = sample_rate_hz
        if csv_time_column is not None:
            t_idx = _col_idx(csv_time_column)
            if t_idx == v_idx:
                raise InputValidationError("CSV time and value columns must be different.")
            times = _values(t_idx)
            dt = np.diff(times)
            if np.any(dt <= 0) or not np.all(np.isfinite(dt)):
                raise InputValidationError(f"CSV time column in '{p.name}' is not strictly increasing and finite.")
            med_dt = float(np.median(dt))
            if float(np.max(np.abs(dt - med_dt)) / med_dt) > 0.02:
                raise InputValidationError(
                    f"CSV time column in '{p.name}' is non-uniformly sampled (>2% jitter). "
                    "BSD requires uniformly sampled time series."
                )
            fs_eff = 1.0 / float((times[-1] - times[0]) / (times.size - 1))
            if sample_rate_hz is not None and (
                not np.isfinite(sample_rate_hz) or sample_rate_hz <= 0
                or not np.isclose(float(sample_rate_hz), fs_eff, rtol=0.02, atol=1e-6)
            ):
                raise InputValidationError("Supplied CSV sample rate conflicts with the time column.")
            t_interval = (float(times[0]), float(times[-1]))
            t_ref = "CSV_TIME_COLUMN"
        elif fs_eff is None:
            raise InputValidationError(
                f"CSV file '{p.name}' has no csv_time_column specified. "
                "Provide either csv_time_column or explicit sample_rate_hz."
            )
        ch_label = header[v_idx] if has_header else f"col_{v_idx}"
        return _validate_and_prepare_signal(
            vals, float(fs_eff), allow_resample, p.name, "csv", ch_label, units, sha, t_interval, t_ref
        )

    # GWOSC HDF5 (.hdf5 / .h5) — strictly checked against O3 GPS guard!
    rec = load_gwosc_hdf5(p)
    assert_not_o3_holdout(p.name, gps_interval=(rec.gps_start, rec.gps_end))
    if rec.n_nonfinite > 0:
        raise InputValidationError(
            f"GWOSC HDF5 '{p.name}' contains {rec.n_nonfinite} non-finite strain samples. "
            "Use clean finite segments without data-quality gaps."
        )
    # Note: Do NOT silently apply adapted_ligo_whitened here; only scale by 1e18 for float64 stability
    # and document it explicitly.
    return _validate_and_prepare_signal(
        rec.strain * 1.0e18,
        float(rec.sample_rate),
        allow_resample,
        p.name,
        "gwosc_hdf5",
        rec.detector,
        "GWOSC_strain_x1e18_unwhitened",
        sha,
        time_interval=(rec.gps_start, rec.gps_end),
        time_reference="GPS",
        provenance={"detector": rec.detector, "gps_start": rec.gps_start, "gps_end": rec.gps_end},
    )


def check_background_independence(
    target: LoadedTimeSeries,
    background: LoadedTimeSeries,
) -> IndependenceAudit:
    """Reject identical files or overlapping time intervals; warn if independence is unverified."""
    rates = (target.sample_rate_hz, background.sample_rate_hz)
    if not all(np.isfinite(rate) and rate > 0 for rate in rates) or not np.isclose(
        rates[0], rates[1], rtol=1e-6, atol=1e-6
    ):
        raise InputValidationError("Target and background must have matching finite effective sample rates.")
    if target.sha256 == background.sha256:
        raise InputValidationError(
            "Background independence violation: Target and Background inputs have identical SHA-256 "
            f"({target.sha256[:12]}...). Calibration data MUST be strictly independent from target data."
        )

    if (
        target.time_interval is not None
        and background.time_interval is not None
        and target.time_reference == background.time_reference
    ):
        t0, t1 = target.time_interval
        b0, b1 = background.time_interval
        if not (t1 <= b0 or t0 >= b1):
            raise InputValidationError(
                f"Background independence violation: Target interval [{t0:.2f}, {t1:.2f}] overlaps "
                f"Background interval [{b0:.2f}, {b1:.2f}] ({target.time_reference}). "
                "Provide disjoint background and target intervals."
            )
        return IndependenceAudit(
            passed=True,
            status="VERIFIED_DISJOINT_TIMESTAMPS",
            target_sha256=target.sha256,
            background_sha256=background.sha256,
            warnings=(),
        )

    return IndependenceAudit(
        passed=True,
        status="UNVERIFIED_EXTERNAL_CHECK_REQUIRED",
        target_sha256=target.sha256,
        background_sha256=background.sha256,
        warnings=(
            "INDEPENDENCE_UNVERIFIED_FROM_METADATA: Files have distinct SHA-256 hashes, but no "
            "comparable absolute timestamps exist. Confirm externally that target and background "
            "recordings do not overlap in time or share leaked signal energy.",
        ),
    )


def slice_into_observations(
    ts: LoadedTimeSeries,
    window_duration_s: float = 4.0,
    prefix: str = "win",
    min_windows: int = 1,
) -> list[Observation]:
    """Slice a validated LoadedTimeSeries into non-overlapping fixed-duration Observations."""
    if window_duration_s <= 0.0:
        raise InputValidationError(f"window_duration_s must be positive, got {window_duration_s}.")
    win_len = int(round(window_duration_s * ts.sample_rate_hz))
    if win_len < 1024:
        raise InputValidationError(
            f"Window length {win_len} samples ({window_duration_s} s at {ts.sample_rate_hz} Hz) "
            "is too short for 512-sample STFT analysis (minimum 1024 samples)."
        )
    n_win = ts.n_samples // win_len
    if n_win < min_windows:
        raise InputValidationError(
            f"Input '{ts.source_name}' yields only {n_win} complete {window_duration_s:.1f}-s windows "
            f"(minimum required: {min_windows}, i.e., >= {min_windows * window_duration_s:.0f} s)."
        )
    obs_list: list[Observation] = []
    for k in range(n_win):
        seg = ts.signal[k * win_len : (k + 1) * win_len].copy()
        obs_list.append(Observation(f"{prefix}_{k:04d}", seg, ts.sample_rate_hz, window_duration_s))
    return obs_list
