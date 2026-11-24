"""Noise models for synthetic time-series generation in Blind Signal Discovery (BSD).

All noise generators are strictly reproducible given a NumPy Generator instance
and return both the simulated noise array and a parameter dictionary for
ground-truth logging (never passed to the detector).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import signal as sp_signal

DEVELOPMENT_NOISE_MODELS: tuple[str, ...] = (
    "gaussian_white",
    "pink_colored",
    "nonstationary",
)

ALL_NOISE_MODELS: tuple[str, ...] = (
    "gaussian_white",
    "pink_colored",
    "nonstationary",
    "impulsive",
    "shifted_ar2",
)


def _normalize_rms(x: np.ndarray, target_sigma: float) -> np.ndarray:
    """Zero-mean and scale array to exact target RMS standard deviation."""
    x_centered = x - np.mean(x)
    rms = float(np.sqrt(np.mean(x_centered**2)))
    if rms <= 0.0:
        return np.zeros_like(x_centered)
    return x_centered * (target_sigma / rms)


def generate_gaussian_white(
    n_samples: int,
    rng: np.random.Generator,
    sigma: float = 1.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate stationary zero-mean i.i.d. Gaussian white noise."""
    noise = rng.normal(loc=0.0, scale=sigma, size=n_samples)
    return noise.astype(np.float64), {"model": "gaussian_white", "sigma": float(sigma)}


def generate_pink_colored(
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    sigma: float = 1.0,
    alpha: float | None = None,
    f_knee: float = 16.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate pink-like 1/f^alpha colored noise with a smooth low-frequency knee.

    Spectral shaping is applied in the real FFT domain so the power spectral
    density scales as 1 / (f^2 + f_knee^2)^(alpha / 2), preventing DC divergence.
    """
    if alpha is None:
        alpha = float(rng.uniform(0.85, 1.15))
    white = rng.normal(loc=0.0, scale=1.0, size=n_samples)
    freqs = np.fft.rfftfreq(n_samples, d=1.0 / sample_rate)
    spectrum = np.fft.rfft(white)
    shaping = 1.0 / np.power(freqs**2 + f_knee**2, alpha / 4.0)
    shaping[0] = 0.0
    colored = np.fft.irfft(spectrum * shaping, n=n_samples)
    colored = _normalize_rms(colored, sigma)
    return colored.astype(np.float64), {
        "model": "pink_colored",
        "sigma": float(sigma),
        "alpha": float(alpha),
        "f_knee_hz": float(f_knee),
    }


def generate_nonstationary(
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    sigma: float = 1.0,
    modulation_depth: float | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate Gaussian noise with smoothly time-varying standard deviation."""
    t = np.arange(n_samples, dtype=np.float64) / sample_rate
    if modulation_depth is None:
        modulation_depth = float(rng.uniform(0.35, 0.70))
    f1 = float(rng.uniform(0.15, 0.45))
    f2 = float(rng.uniform(0.50, 0.90))
    phi1 = float(rng.uniform(0.0, 2.0 * np.pi))
    phi2 = float(rng.uniform(0.0, 2.0 * np.pi))

    envelope = (
        1.0
        + modulation_depth * np.sin(2.0 * np.pi * f1 * t + phi1)
        + 0.35 * modulation_depth * np.cos(2.0 * np.pi * f2 * t + phi2)
    )
    envelope = np.clip(envelope, 0.25, None)
    envelope = envelope / np.sqrt(np.mean(envelope**2))

    base = rng.normal(loc=0.0, scale=1.0, size=n_samples)
    noise = _normalize_rms(base * envelope, sigma)
    return noise.astype(np.float64), {
        "model": "nonstationary",
        "sigma": float(sigma),
        "modulation_depth": float(modulation_depth),
        "f1_hz": f1,
        "f2_hz": f2,
        "min_envelope": float(np.min(envelope)),
        "max_envelope": float(np.max(envelope)),
    }


def generate_impulsive(
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    sigma: float = 1.0,
    n_impulses: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate Gaussian background contaminated by sparse high-amplitude impulses.

    Note: The background Gaussian component has standard deviation `sigma` before
    impulses are added, so impulses represent genuine non-Gaussian outliers.
    """
    noise = rng.normal(loc=0.0, scale=sigma, size=n_samples)
    if n_impulses is None:
        n_impulses = int(rng.integers(1, 5))

    impulse_indices: list[int] = []
    impulse_amplitudes: list[float] = []
    margin = max(16, int(0.02 * n_samples))
    for _ in range(n_impulses):
        idx = int(rng.integers(margin, n_samples - margin))
        width = int(rng.integers(1, 4))
        sign = float(rng.choice([-1.0, 1.0]))
        amp = float(rng.uniform(6.0, 12.0)) * sigma * sign
        end_idx = min(n_samples, idx + width)
        noise[idx:end_idx] += amp
        impulse_indices.append(idx)
        impulse_amplitudes.append(amp)

    return noise.astype(np.float64), {
        "model": "impulsive",
        "sigma": float(sigma),
        "n_impulses": n_impulses,
        "impulse_times_s": [idx / sample_rate for idx in impulse_indices],
        "impulse_amplitudes": impulse_amplitudes,
    }


def generate_shifted_ar2(
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    sigma: float = 1.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate out-of-distribution resonant AR(2) colored noise with mild drift.

    This distribution is withheld from development and calibration to evaluate
    detector robustness under noise-model mismatch.
    """
    f_res = float(rng.uniform(110.0, 340.0))
    pole_radius = float(rng.uniform(0.72, 0.86))
    omega = 2.0 * np.pi * f_res / sample_rate
    a_coeffs = np.array([1.0, -2.0 * pole_radius * np.cos(omega), pole_radius**2])
    b_coeffs = np.array([1.0])

    # Burn-in to remove filter startup transient
    burn_in = 512
    white = rng.normal(loc=0.0, scale=1.0, size=n_samples + burn_in)
    filtered = sp_signal.lfilter(b_coeffs, a_coeffs, white)[burn_in:]

    # Mild slow nonstationary drift
    t = np.arange(n_samples, dtype=np.float64) / sample_rate
    drift_phase = float(rng.uniform(0.0, 2.0 * np.pi))
    envelope = 1.0 + 0.25 * np.sin(2.0 * np.pi * 0.25 * t + drift_phase)
    noise = _normalize_rms(filtered * envelope, sigma)
    return noise.astype(np.float64), {
        "model": "shifted_ar2",
        "sigma": float(sigma),
        "resonant_freq_hz": f_res,
        "pole_radius": pole_radius,
    }


def generate_noise(
    noise_model: str,
    n_samples: int,
    sample_rate: float,
    rng: np.random.Generator,
    sigma: float = 1.0,
    params: dict[str, Any] | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Dispatch noise generation by model name."""
    params = params or {}
    if noise_model == "gaussian_white":
        return generate_gaussian_white(n_samples=n_samples, rng=rng, sigma=sigma)
    if noise_model == "pink_colored":
        return generate_pink_colored(
            n_samples=n_samples,
            sample_rate=sample_rate,
            rng=rng,
            sigma=sigma,
            alpha=params.get("alpha"),
            f_knee=float(params.get("f_knee", 16.0)),
        )
    if noise_model == "nonstationary":
        return generate_nonstationary(
            n_samples=n_samples,
            sample_rate=sample_rate,
            rng=rng,
            sigma=sigma,
            modulation_depth=params.get("modulation_depth"),
        )
    if noise_model == "impulsive":
        return generate_impulsive(
            n_samples=n_samples,
            sample_rate=sample_rate,
            rng=rng,
            sigma=sigma,
            n_impulses=params.get("n_impulses"),
        )
    if noise_model == "shifted_ar2":
        return generate_shifted_ar2(
            n_samples=n_samples,
            sample_rate=sample_rate,
            rng=rng,
            sigma=sigma,
        )
    raise ValueError(f"Unknown noise model: {noise_model!r}. Expected one of {ALL_NOISE_MODELS}.")
