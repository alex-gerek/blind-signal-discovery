"""Synthetic time-series data generator for Blind Signal Discovery (BSD).

Strict ground-truth isolation is enforced by construction:
- `Observation` contains ONLY the raw time-series array and acquisition rate.
- `GroundTruth` stores all injection metadata, signal labels, SNR values, and
  true time-frequency bounds separately for blind evaluation only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy.signal.windows import tukey

from bsd.noise import generate_noise

KNOWN_SIGNAL_FAMILIES: tuple[str, ...] = (
    "sinusoid",
    "broadband_pulse",
    "linear_chirp",
    "intermittent_narrowband",
)

WITHHELD_SIGNAL_FAMILIES: tuple[str, ...] = ("sinusoidal_fm",)

ALL_SIGNAL_FAMILIES: tuple[str, ...] = KNOWN_SIGNAL_FAMILIES + WITHHELD_SIGNAL_FAMILIES

# Fixed integer stream codes guarantee statistically independent RNG streams
# across datasets when combined with a master experiment seed via SeedSequence.
SPLIT_STREAM_CODES: dict[str, int] = {
    "development": 1001,
    "calibration": 2002,
    "test_noise": 3003,
    "test_signal": 4004,
    "test_withheld": 5005,
    "test_shifted_noise": 6006,
    "test_shifted_signal": 7007,
    "robustness_leakage": 8001,
    "robustness_windowing": 8002,
    "robustness_boundary": 8003,
    "robustness_impulsive": 8004,
    "robustness_nonstationary": 8005,
    "robustness_overlapping": 8006,
}


@dataclass(frozen=True)
class Observation:
    """Detector-visible observation containing no ground-truth metadata."""

    sample_id: str
    signal: np.ndarray
    sample_rate: float
    duration: float


@dataclass(frozen=True)
class GroundTruth:
    """Evaluation-only ground truth stored separately from observations."""

    sample_id: str
    is_signal: bool
    signal_family: str | None
    noise_model: str
    target_snr: float | None
    empirical_energy_snr: float | None
    rms_snr: float | None
    peak_snr: float | None
    t_start: float | None
    t_end: float | None
    f_min: float | None
    f_max: float | None
    params: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DatasetSplit:
    """Paired observations and ground truths for a named, independent split."""

    split_name: str
    seed: int
    observations: list[Observation]
    ground_truths: list[GroundTruth]

    def __len__(self) -> int:
        return len(self.observations)


def make_split_rng(master_seed: int, split_name: str) -> np.random.Generator:
    """Create a deterministic, independent RNG for a named dataset split."""
    if split_name in SPLIT_STREAM_CODES:
        code = SPLIT_STREAM_CODES[split_name]
    else:
        # Deterministic fallback hash using simple polynomial rolling code
        code = 9000 + sum((idx + 1) * ord(ch) for idx, ch in enumerate(split_name)) % 100000
    ss = np.random.SeedSequence([int(master_seed), int(code)])
    return np.random.default_rng(ss)


def _scale_to_energy_snr(
    clean_signal: np.ndarray,
    noise_sigma: float,
    target_energy_snr: float,
) -> np.ndarray:
    """Scale clean signal so sqrt(sum(s^2)) / noise_sigma == target_energy_snr."""
    energy_norm = float(np.linalg.norm(clean_signal))
    if energy_norm <= 0.0:
        raise ValueError("Generated clean signal has zero energy.")
    return clean_signal * (target_energy_snr * noise_sigma / energy_norm)


def _generate_sinusoid(
    t: np.ndarray,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    params = params or {}
    nyquist = sample_rate / 2.0
    freq = float(params.get("freq_hz", rng.uniform(0.10 * nyquist, 0.85 * nyquist)))
    phase = float(params.get("phase_rad", rng.uniform(0.0, 2.0 * np.pi)))
    taper_alpha = float(params.get("taper_alpha", 0.05))
    window = tukey(len(t), alpha=taper_alpha) if taper_alpha > 0 else np.ones_like(t)
    sig = window * np.sin(2.0 * np.pi * freq * t + phase)
    df_half = max(4.0, sample_rate / 128.0)
    return (
        sig,
        float(t[0]),
        float(t[-1]),
        max(0.0, freq - df_half),
        min(nyquist, freq + df_half),
        {"freq_hz": freq, "phase_rad": phase, "taper_alpha": taper_alpha},
    )


def _generate_broadband_pulse(
    t: np.ndarray,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    params = params or {}
    duration = float(t[-1] - t[0])
    nyquist = sample_rate / 2.0
    t0 = float(params.get("t_center_s", rng.uniform(0.15 * duration, 0.85 * duration)))
    sigma_t = float(params.get("sigma_t_s", rng.uniform(0.0025, 0.0065)))
    f_center = float(params.get("f_center_hz", rng.uniform(0.28 * nyquist, 0.72 * nyquist)))
    phase = float(params.get("phase_rad", rng.uniform(0.0, 2.0 * np.pi)))

    envelope = np.exp(-0.5 * ((t - t0) / sigma_t) ** 2)
    sig = envelope * np.cos(2.0 * np.pi * f_center * (t - t0) + phase)
    # Approximate 2-sigma temporal and spectral support
    half_bw = 2.0 / (2.0 * np.pi * sigma_t)
    t_start = max(float(t[0]), t0 - 3.0 * sigma_t)
    t_end = min(float(t[-1]), t0 + 3.0 * sigma_t)
    f_min = max(0.0, f_center - half_bw)
    f_max = min(nyquist, f_center + half_bw)
    return (
        sig,
        t_start,
        t_end,
        f_min,
        f_max,
        {
            "t_center_s": t0,
            "sigma_t_s": sigma_t,
            "f_center_hz": f_center,
            "phase_rad": phase,
        },
    )


def _generate_linear_chirp(
    t: np.ndarray,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    params = params or {}
    duration = float(t[-1] - t[0])
    nyquist = sample_rate / 2.0
    chirp_dur = float(params.get("chirp_duration_s", rng.uniform(0.22 * duration, 0.55 * duration)))
    margin = float(params.get("margin_s", 0.12 * duration))
    max_start = max(margin, duration - margin - chirp_dur)
    t_start = float(params.get("t_start_s", rng.uniform(margin, max_start)))
    t_end = min(duration, t_start + chirp_dur)

    if "f_start_hz" in params and "f_end_hz" in params:
        f_start = float(params["f_start_hz"])
        f_end = float(params["f_end_hz"])
    else:
        f_start = float(rng.uniform(0.12 * nyquist, 0.78 * nyquist))
        sweep = float(rng.uniform(0.16 * nyquist, 0.35 * nyquist))
        direction = 1.0 if (f_start + sweep <= 0.88 * nyquist and rng.random() < 0.5) else -1.0
        if f_start + direction * sweep < 0.10 * nyquist:
            direction = 1.0
        f_end = float(np.clip(f_start + direction * sweep, 0.10 * nyquist, 0.88 * nyquist))

    phase0 = float(params.get("phase_rad", rng.uniform(0.0, 2.0 * np.pi)))
    mask = (t >= t_start) & (t <= t_end)
    n_active = int(np.sum(mask))
    sig = np.zeros_like(t)
    if n_active > 4:
        tau = t[mask] - t_start
        actual_dur = max(float(tau[-1]), 1.0 / sample_rate)
        slope = (f_end - f_start) / actual_dur
        inst_phase = 2.0 * np.pi * (f_start * tau + 0.5 * slope * tau**2) + phase0
        sig[mask] = tukey(n_active, alpha=0.15) * np.sin(inst_phase)

    return (
        sig,
        t_start,
        t_end,
        min(f_start, f_end),
        max(f_start, f_end),
        {
            "t_start_s": t_start,
            "t_end_s": t_end,
            "f_start_hz": f_start,
            "f_end_hz": f_end,
            "phase_rad": phase0,
        },
    )


def _generate_intermittent_narrowband(
    t: np.ndarray,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    params = params or {}
    duration = float(t[-1] - t[0])
    nyquist = sample_rate / 2.0
    freq = float(params.get("freq_hz", rng.uniform(0.12 * nyquist, 0.84 * nyquist)))
    n_bursts = int(params.get("n_bursts", rng.integers(2, 5)))
    phase = float(params.get("phase_rad", rng.uniform(0.0, 2.0 * np.pi)))

    sig = np.zeros_like(t)
    slot_width = (0.84 * duration) / n_bursts
    burst_intervals: list[list[float]] = []

    for idx in range(n_bursts):
        slot_start = 0.08 * duration + idx * slot_width
        burst_dur = float(rng.uniform(0.35 * slot_width, 0.65 * slot_width))
        burst_start = float(rng.uniform(slot_start, slot_start + slot_width - burst_dur))
        burst_end = burst_start + burst_dur
        mask = (t >= burst_start) & (t <= burst_end)
        n_active = int(np.sum(mask))
        if n_active > 4:
            sig[mask] += tukey(n_active, alpha=0.25) * np.sin(2.0 * np.pi * freq * t[mask] + phase)
            burst_intervals.append([burst_start, burst_end])

    t_start = min(iv[0] for iv in burst_intervals)
    t_end = max(iv[1] for iv in burst_intervals)
    df_half = max(6.0, sample_rate / 128.0)
    return (
        sig,
        t_start,
        t_end,
        max(0.0, freq - df_half),
        min(nyquist, freq + df_half),
        {
            "freq_hz": freq,
            "n_bursts": n_bursts,
            "burst_intervals_s": burst_intervals,
            "phase_rad": phase,
        },
    )


def _generate_sinusoidal_fm(
    t: np.ndarray,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    """Withheld signal family: time-localized sinusoidal frequency-modulated wave."""
    params = params or {}
    duration = float(t[-1] - t[0])
    nyquist = sample_rate / 2.0
    sig_dur = float(params.get("signal_duration_s", rng.uniform(0.28 * duration, 0.55 * duration)))
    margin = 0.12 * duration
    t_start = float(params.get("t_start_s", rng.uniform(margin, duration - margin - sig_dur)))
    t_end = t_start + sig_dur

    f_carrier = float(params.get("f_carrier_hz", rng.uniform(0.26 * nyquist, 0.72 * nyquist)))
    f_mod = float(params.get("f_mod_hz", rng.uniform(1.5, 4.0)))
    freq_dev = float(params.get("freq_dev_hz", rng.uniform(35.0, 85.0)))
    beta = freq_dev / f_mod
    phi0 = float(rng.uniform(0.0, 2.0 * np.pi))
    phi_m = float(rng.uniform(0.0, 2.0 * np.pi))

    mask = (t >= t_start) & (t <= t_end)
    n_active = int(np.sum(mask))
    sig = np.zeros_like(t)
    if n_active > 4:
        tau = t[mask] - t_start
        inst_phase = 2.0 * np.pi * f_carrier * tau + beta * np.sin(2.0 * np.pi * f_mod * tau + phi_m) + phi0
        sig[mask] = tukey(n_active, alpha=0.15) * np.sin(inst_phase)

    f_min = max(0.0, f_carrier - freq_dev)
    f_max = min(nyquist, f_carrier + freq_dev)
    return (
        sig,
        t_start,
        t_end,
        f_min,
        f_max,
        {
            "t_start_s": t_start,
            "t_end_s": t_end,
            "f_carrier_hz": f_carrier,
            "f_mod_hz": f_mod,
            "freq_dev_hz": freq_dev,
            "modulation_index": beta,
        },
    )


def generate_clean_signal(
    signal_family: str,
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, float, float, float, float, dict[str, Any]]:
    """Generate a unit-unscaled clean signal and its ground-truth TF bounds."""
    t = np.arange(n_samples, dtype=np.float64) / sample_rate
    if signal_family == "sinusoid":
        return _generate_sinusoid(t, sample_rate, rng, params)
    if signal_family == "broadband_pulse":
        return _generate_broadband_pulse(t, sample_rate, rng, params)
    if signal_family == "linear_chirp":
        return _generate_linear_chirp(t, sample_rate, rng, params)
    if signal_family == "intermittent_narrowband":
        return _generate_intermittent_narrowband(t, sample_rate, rng, params)
    if signal_family == "sinusoidal_fm":
        return _generate_sinusoidal_fm(t, sample_rate, rng, params)
    raise ValueError(
        f"Unknown signal family: {signal_family!r}. Expected one of {ALL_SIGNAL_FAMILIES}."
    )


def generate_sample(
    sample_id: str,
    sample_rate: float,
    duration: float,
    noise_model: str,
    rng: np.random.Generator,
    signal_family: str | None = None,
    target_snr: float | None = None,
    noise_sigma: float = 1.0,
    noise_params: dict[str, Any] | None = None,
    signal_params: dict[str, Any] | None = None,
    secondary_signal: tuple[str, float, dict[str, Any] | None] | None = None,
) -> tuple[Observation, GroundTruth]:
    """Generate a single simulated observation and its isolated ground truth.

    Parameters
    ----------
    target_snr:
        Total energy SNR defined as ``sqrt(sum(s^2)) / noise_sigma``.
    secondary_signal:
        Optional ``(family, snr, params)`` tuple to inject a second overlapping
        signal for robustness testing.
    """
    n_samples = int(round(sample_rate * duration))
    noise, used_noise_params = generate_noise(
        noise_model=noise_model,
        n_samples=n_samples,
        sample_rate=sample_rate,
        rng=rng,
        sigma=noise_sigma,
        params=noise_params,
    )

    if signal_family is None:
        obs = Observation(
            sample_id=sample_id,
            signal=noise,
            sample_rate=float(sample_rate),
            duration=float(duration),
        )
        gt = GroundTruth(
            sample_id=sample_id,
            is_signal=False,
            signal_family=None,
            noise_model=noise_model,
            target_snr=None,
            empirical_energy_snr=None,
            rms_snr=None,
            peak_snr=None,
            t_start=None,
            t_end=None,
            f_min=None,
            f_max=None,
            params={"noise": used_noise_params},
        )
        return obs, gt

    if target_snr is None or target_snr <= 0.0:
        raise ValueError("target_snr must be positive when signal_family is provided.")

    clean_raw, t_start, t_end, f_min, f_max, used_sig_params = generate_clean_signal(
        signal_family=signal_family,
        n_samples=n_samples,
        sample_rate=sample_rate,
        rng=rng,
        params=signal_params,
    )
    clean_scaled = _scale_to_energy_snr(clean_raw, noise_sigma=noise_sigma, target_energy_snr=target_snr)

    extra_params: dict[str, Any] = {}
    if secondary_signal is not None:
        sec_family, sec_snr, sec_params = secondary_signal
        sec_raw, s_t0, s_t1, s_f0, s_f1, sec_used = generate_clean_signal(
            signal_family=sec_family,
            n_samples=n_samples,
            sample_rate=sample_rate,
            rng=rng,
            params=sec_params,
        )
        sec_scaled = _scale_to_energy_snr(sec_raw, noise_sigma=noise_sigma, target_energy_snr=sec_snr)
        clean_scaled = clean_scaled + sec_scaled
        t_start = min(t_start, s_t0)
        t_end = max(t_end, s_t1)
        f_min = min(f_min, s_f0)
        f_max = max(f_max, s_f1)
        extra_params["secondary_signal"] = {"family": sec_family, "snr": sec_snr, **sec_used}

    observed_signal = noise + clean_scaled
    emp_noise_rms = float(np.std(noise))
    emp_sig_norm = float(np.linalg.norm(clean_scaled))
    emp_energy_snr = emp_sig_norm / emp_noise_rms if emp_noise_rms > 0 else 0.0
    rms_snr = float(np.sqrt(np.mean(clean_scaled**2))) / emp_noise_rms if emp_noise_rms > 0 else 0.0
    peak_snr = float(np.max(np.abs(clean_scaled))) / emp_noise_rms if emp_noise_rms > 0 else 0.0

    obs = Observation(
        sample_id=sample_id,
        signal=observed_signal.astype(np.float64),
        sample_rate=float(sample_rate),
        duration=float(duration),
    )
    gt = GroundTruth(
        sample_id=sample_id,
        is_signal=True,
        signal_family=signal_family,
        noise_model=noise_model,
        target_snr=float(target_snr),
        empirical_energy_snr=float(emp_energy_snr),
        rms_snr=float(rms_snr),
        peak_snr=float(peak_snr),
        t_start=float(t_start),
        t_end=float(t_end),
        f_min=float(f_min),
        f_max=float(f_max),
        params={"noise": used_noise_params, "signal": used_sig_params, **extra_params},
    )
    return obs, gt


def generate_noise_dataset(
    split_name: str,
    n_samples_count: int,
    sample_rate: float,
    duration: float,
    noise_models: list[str] | tuple[str, ...],
    master_seed: int,
) -> DatasetSplit:
    """Generate a balanced noise-only dataset split with deterministic seeding."""
    rng = make_split_rng(master_seed, split_name)
    observations: list[Observation] = []
    ground_truths: list[GroundTruth] = []

    for idx in range(n_samples_count):
        noise_model = noise_models[idx % len(noise_models)]
        sample_id = f"{split_name}_s{master_seed}_{idx:05d}"
        obs, gt = generate_sample(
            sample_id=sample_id,
            sample_rate=sample_rate,
            duration=duration,
            noise_model=noise_model,
            rng=rng,
            signal_family=None,
        )
        observations.append(obs)
        ground_truths.append(gt)

    return DatasetSplit(
        split_name=split_name,
        seed=master_seed,
        observations=observations,
        ground_truths=ground_truths,
    )


def generate_signal_dataset(
    split_name: str,
    n_samples_count: int,
    sample_rate: float,
    duration: float,
    signal_families: list[str] | tuple[str, ...],
    snr_levels: list[float] | tuple[float, ...],
    noise_models: list[str] | tuple[str, ...],
    master_seed: int,
) -> DatasetSplit:
    """Generate a stratified signal-injected dataset split across families, SNRs, and noise."""
    rng = make_split_rng(master_seed, split_name)
    observations: list[Observation] = []
    ground_truths: list[GroundTruth] = []

    combos = [
        (fam, float(snr), nm)
        for fam in signal_families
        for snr in snr_levels
        for nm in noise_models
    ]

    for idx in range(n_samples_count):
        fam, snr, noise_model = combos[idx % len(combos)]
        sample_id = f"{split_name}_s{master_seed}_{idx:05d}"
        obs, gt = generate_sample(
            sample_id=sample_id,
            sample_rate=sample_rate,
            duration=duration,
            noise_model=noise_model,
            rng=rng,
            signal_family=fam,
            target_snr=snr,
        )
        observations.append(obs)
        ground_truths.append(gt)

    return DatasetSplit(
        split_name=split_name,
        seed=master_seed,
        observations=observations,
        ground_truths=ground_truths,
    )
