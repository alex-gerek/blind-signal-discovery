"""BSD v0.3 Phase 2b — Multi-Window Cross-Detector Coherence Prototype.

Treats cross-detector coherence strictly as additional corroborating evidence,
NEVER automatic proof of an astrophysical signal. Explicitly distinguishes:
  - INCOHERENT_CONSISTENT_WITH_NULL (consistent with time-slide null)
  - COHERENT_COMMON_ARTIFACT_WARNING (narrowband line or synchronous impulsive glitch)
  - COHERENT_ZERO_LAG_ENV_WARNING (peak lag <= 0.5 ms: correlated environmental noise
    or face-on signal — never automatically classified as astrophysical)
  - COHERENT_PHYSICAL_LAG_EVIDENCE (broadband non-impulsive excess with physical
    inter-detector delay 0.5 ms < |dt| <= 10 ms above time-slide null threshold)
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import signal as sp_signal

from bsd.generator import Observation

COH_NULL = "INCOHERENT_CONSISTENT_WITH_NULL"
COH_ARTIFACT = "COHERENT_COMMON_ARTIFACT_WARNING"
COH_ZERO_LAG_ENV = "COHERENT_ZERO_LAG_ENV_WARNING"
COH_PHYSICAL_EVIDENCE = "COHERENT_PHYSICAL_LAG_EVIDENCE"


@dataclass(frozen=True)
class CrossCoherenceConfig:
    """Configuration for multi-window cross-detector coherence evaluation."""

    band_hz: tuple[float, float] = (50.0, 350.0)
    max_physical_lag_ms: float = 10.0
    zero_lag_guard_ms: float = 0.5
    n_time_slides: int = 120
    min_slide_s: float = 0.25
    max_slide_s: float = 1.75
    alpha_fap: float = 0.01
    welch_nperseg: int = 512
    narrowband_concentration_thresh: float = 0.45
    impulsive_block_concentration_thresh: float = 0.35


@dataclass(frozen=True)
class CrossCoherenceResult:
    """Auditable cross-detector coherence summary with time-slide null control."""

    peak_xcorr_physical: float
    peak_lag_ms: float
    mean_band_msc: float
    top3_bin_msc_fraction: float
    time_slide_p_value: float
    time_slide_threshold: float
    significant: bool
    classification: str
    n_windows: int
    n_time_slides: int
    diagnostics: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _bandpass_whiten(x: np.ndarray, fs: float, band_hz: tuple[float, float]) -> np.ndarray:
    x = x - float(np.median(x))
    sos = sp_signal.butter(4, [band_hz[0], band_hz[1]], btype="bandpass", fs=fs, output="sos")
    xf = sp_signal.sosfiltfilt(sos, x)
    std = float(np.std(xf))
    return xf / max(std, 1e-18)


def _peak_lag_xcorr(x1: np.ndarray, x2: np.ndarray, max_lag_samples: int) -> tuple[float, int]:
    """Normalized FFT cross-correlation peak within [-max_lag_samples, +max_lag_samples]."""
    n = x1.size
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f1 = np.fft.rfft(x1, n=nfft)
    f2 = np.fft.rfft(x2, n=nfft)
    cc = np.fft.irfft(f1 * np.conj(f2), n=nfft) / max(1, n)
    # Positive lags: 0..max_lag, negative lags: nfft-max_lag..nfft-1
    pos = cc[: max_lag_samples + 1]
    neg = cc[nfft - max_lag_samples :]
    windowed = np.concatenate([neg, pos])
    lags = np.arange(-max_lag_samples, max_lag_samples + 1, dtype=int)
    best_i = int(np.argmax(np.abs(windowed)))
    return float(np.abs(windowed[best_i])), int(lags[best_i])


def evaluate_cross_coherence(
    obs_list_1: list[Observation],
    obs_list_2: list[Observation],
    config: CrossCoherenceConfig | None = None,
) -> CrossCoherenceResult:
    """Evaluate multi-window cross-detector coherence against circular time-slide controls."""
    cfg = config or CrossCoherenceConfig()
    if len(obs_list_1) != len(obs_list_2) or not obs_list_1:
        raise ValueError("obs_list_1 and obs_list_2 must be non-empty and of equal length.")

    fs = float(obs_list_1[0].sample_rate)
    x1 = _bandpass_whiten(np.concatenate([o.signal for o in obs_list_1]), fs, cfg.band_hz)
    x2 = _bandpass_whiten(np.concatenate([o.signal for o in obs_list_2]), fs, cfg.band_hz)
    n = x1.size

    max_lag_samp = max(1, int(round((cfg.max_physical_lag_ms / 1000.0) * fs)))
    peak_cc, peak_lag_samp = _peak_lag_xcorr(x1, x2, max_lag_samp)
    peak_lag_ms = (peak_lag_samp / fs) * 1000.0

    # Multi-window magnitude-squared coherence (MSC) in target band
    freqs, msc = sp_signal.coherence(
        x1, x2, fs=fs, window="hann", nperseg=min(cfg.welch_nperseg, n // 4), noverlap=min(cfg.welch_nperseg, n // 4) // 2
    )
    bmask = (freqs >= cfg.band_hz[0]) & (freqs <= cfg.band_hz[1])
    band_msc = msc[bmask] if np.any(bmask) else msc
    mean_msc = float(np.mean(band_msc))
    top3_frac = float(np.sum(np.sort(band_msc)[-3:]) / max(float(np.sum(band_msc)), 1e-18))

    # Time-slide null distribution (non-physical lags)
    slide_offsets_s = np.linspace(cfg.min_slide_s, cfg.max_slide_s, cfg.n_time_slides)
    null_peaks = np.empty(cfg.n_time_slides, dtype=np.float64)
    for k, dt_s in enumerate(slide_offsets_s):
        shift = int(round(dt_s * fs)) % n
        x2_shifted = np.roll(x2, shift)
        p_null, _ = _peak_lag_xcorr(x1, x2_shifted, max_lag_samp)
        null_peaks[k] = p_null

    n_ge = int(np.sum(null_peaks >= peak_cc))
    p_slide = float((1.0 + n_ge) / (cfg.n_time_slides + 1.0))
    k_th = min(cfg.n_time_slides, max(1, int(np.ceil((1.0 - cfg.alpha_fap) * (cfg.n_time_slides + 1)))))
    thresh = float(np.sort(null_peaks)[k_th - 1])
    sig = bool(peak_cc > thresh)

    # Check common-artifact signatures:
    # (a) narrowband line concentration in top-3 MSC bins
    # (b) synchronous impulsive glitch: fraction of total |cross-product| energy
    #     concentrated in the top 1% of 64-sample (~62 ms) time blocks.
    #     (Stationary Gaussian broadband product has ~0.04; synchronous glitches > 0.35.)
    x2_aligned = np.roll(x2, peak_lag_samp)
    prod = np.abs(x1 * x2_aligned)
    blk = 64
    n_blk = prod.size // blk
    blk_means = prod[: n_blk * blk].reshape(n_blk, blk).mean(axis=1)
    top_k = max(1, int(np.ceil(0.01 * n_blk)))
    top1pct_blk_frac = float(np.sum(np.sort(blk_means)[-top_k:]) / max(float(np.sum(blk_means)), 1e-18))

    max_msc = float(np.max(band_msc))
    is_narrowband_line = bool(max_msc >= 0.50 and top3_frac >= cfg.narrowband_concentration_thresh)
    is_impulsive_artifact = bool(top1pct_blk_frac >= cfg.impulsive_block_concentration_thresh)
    is_zero_lag = bool(abs(peak_lag_ms) <= cfg.zero_lag_guard_ms)

    if is_narrowband_line:
        # Shared periodic narrowband line elevates MSC directly even when circular
        # time-slides remain periodic; flag explicitly as common artifact warning.
        cls = COH_ARTIFACT
    elif not sig:
        cls = COH_NULL
    elif is_impulsive_artifact:
        cls = COH_ARTIFACT
    elif is_zero_lag:
        cls = COH_ZERO_LAG_ENV
    else:
        cls = COH_PHYSICAL_EVIDENCE

    return CrossCoherenceResult(
        peak_xcorr_physical=peak_cc,
        peak_lag_ms=peak_lag_ms,
        mean_band_msc=mean_msc,
        top3_bin_msc_fraction=top3_frac,
        time_slide_p_value=p_slide,
        time_slide_threshold=thresh,
        significant=sig,
        classification=cls,
        n_windows=len(obs_list_1),
        n_time_slides=cfg.n_time_slides,
        diagnostics={
            "top1pct_block_energy_frac": top1pct_blk_frac,
            "is_narrowband_line": is_narrowband_line,
            "is_impulsive_artifact": is_impulsive_artifact,
            "is_zero_lag": is_zero_lag,
        },
    )
