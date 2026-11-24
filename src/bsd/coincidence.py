"""BSD v0.2 Phase 3 — Sub-10 ms Cross-Detector Timing & Time-Slide Coincidence.

Provides:
  1. Band-limited normalized time-domain cross-correlation `estimate_cross_detector_delay`
     over a configurable physical lag window (default +-15 ms) with parabolic
     sub-sample interpolation. Evaluates |r(tau)| to accommodate opposite
     antenna-pattern polarity between H1 and L1.
  2. Pairwise candidate coincidence matching (`match_coincident_candidates`)
     combining temporal overlap/proximity with cross-correlation lag estimates.
  3. Strengthened unphysical time-slide background analysis (`run_time_slide_analysis`)
     reporting accidental coincidence counts, rates per slide and per hour of
     equivalent time-slide livetime, and explicit resolution-bounded empirical
     p-values (`p >= 1 / (N_slides + 1)`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import signal as sp_signal

PHYSICAL_H1_L1_LIGHT_TRAVEL_S = 0.010


@dataclass(frozen=True)
class CrossCorrelationResult:
    """Result of band-limited cross-detector lag estimation."""

    delay_s: float  # Estimated t_A - t_B (positive => signal arrives later in A than B)
    abs_peak_corr: float
    signed_peak_corr: float
    within_10ms_light_travel: bool
    band_hz: tuple[float, float]
    max_lag_s: float
    lags_s: np.ndarray | None = None
    corr_curve: np.ndarray | None = None

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "delay_s": self.delay_s,
            "delay_ms": self.delay_s * 1000.0,
            "abs_peak_corr": self.abs_peak_corr,
            "signed_peak_corr": self.signed_peak_corr,
            "within_10ms_light_travel": self.within_10ms_light_travel,
            "band_hz": list(self.band_hz),
            "max_lag_s": self.max_lag_s,
        }


@dataclass(frozen=True)
class CoincidentPair:
    """A coincident candidate pair across two detectors."""

    det_a: str
    det_b: str
    gps_peak_a: float
    gps_peak_b: float
    stft_dt_a_minus_b_s: float
    interval_gap_s: float
    score_a: float
    score_b: float
    xcorr_delay_s: float | None = None
    abs_peak_corr: float | None = None
    signed_peak_corr: float | None = None
    within_10ms_light_travel: bool | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def estimate_cross_detector_delay(
    sig_a: np.ndarray,
    sig_b: np.ndarray,
    sample_rate: float,
    band_hz: tuple[float, float] = (35.0, 350.0),
    max_lag_s: float = 0.015,
    return_curve: bool = False,
) -> CrossCorrelationResult:
    """Estimate inter-detector arrival delay (t_A - t_B) via band-limited |r(tau)|.

    Convention: if ``sig_a(t) = +- sig_b(t - delay)``, then the signal arrives
    ``delay`` seconds later in detector A than in detector B, and the returned
    ``delay_s`` is approximately ``+delay``.
    """
    a = np.asarray(sig_a, dtype=np.float64)
    b = np.asarray(sig_b, dtype=np.float64)
    if a.ndim != 1 or b.ndim != 1 or a.size != b.size:
        raise ValueError("sig_a and sig_b must be 1D arrays of equal length.")
    if a.size < 32:
        raise ValueError("Signals must have at least 32 samples for cross-correlation.")

    a = a - np.mean(a)
    b = b - np.mean(b)

    nyq = 0.5 * sample_rate
    f_lo = max(1.0, min(band_hz[0], nyq * 0.90))
    f_hi = min(nyq * 0.95, max(band_hz[1], f_lo + 10.0))
    if a.size >= 64 and f_hi > f_lo:
        sos = sp_signal.butter(4, [f_lo, f_hi], btype="bandpass", fs=sample_rate, output="sos")
        a = sp_signal.sosfiltfilt(sos, a)
        b = sp_signal.sosfiltfilt(sos, b)

    max_lag_samples = max(1, min(int(round(max_lag_s * sample_rate)), a.size // 2 - 1))
    lags_samples = np.arange(-max_lag_samples, max_lag_samples + 1, dtype=np.int64)
    corr = np.zeros(lags_samples.size, dtype=np.float64)

    for idx, lag in enumerate(lags_samples):
        # Positive lag: compare a[lag:] with b[:-lag] (signal later in A than B)
        if lag > 0:
            va = a[lag:]
            vb = b[:-lag]
        elif lag < 0:
            va = a[:lag]
            vb = b[-lag:]
        else:
            va = a
            vb = b
        denom = float(np.linalg.norm(va) * np.linalg.norm(vb))
        corr[idx] = float(np.dot(va, vb) / denom) if denom > 0.0 else 0.0

    abs_corr = np.abs(corr)
    peak_idx = int(np.argmax(abs_corr))
    signed_peak = float(corr[peak_idx])
    abs_peak = float(abs_corr[peak_idx])

    # Parabolic sub-sample interpolation around peak of |r(tau)|
    frac_offset = 0.0
    if 0 < peak_idx < (abs_corr.size - 1):
        y_m1 = abs_corr[peak_idx - 1]
        y_0 = abs_corr[peak_idx]
        y_p1 = abs_corr[peak_idx + 1]
        denom_p = (y_m1 - 2.0 * y_0 + y_p1)
        if abs(denom_p) > 1e-12:
            frac_offset = float(np.clip(0.5 * (y_m1 - y_p1) / denom_p, -0.5, 0.5))

    delay_samples = float(lags_samples[peak_idx]) + frac_offset
    delay_s = delay_samples / sample_rate
    lags_s = lags_samples / sample_rate

    return CrossCorrelationResult(
        delay_s=delay_s,
        abs_peak_corr=abs_peak,
        signed_peak_corr=signed_peak,
        within_10ms_light_travel=bool(abs(delay_s) <= PHYSICAL_H1_L1_LIGHT_TRAVEL_S),
        band_hz=(f_lo, f_hi),
        max_lag_s=max_lag_s,
        lags_s=lags_s if return_curve else None,
        corr_curve=corr if return_curve else None,
    )


def match_coincident_candidates(
    cands_a: list[dict[str, Any]],
    cands_b: list[dict[str, Any]],
    max_dt_s: float = 0.05,
) -> list[CoincidentPair]:
    """Match candidate pairs across two detectors within ``max_dt_s``."""
    pairs: list[CoincidentPair] = []
    for ca in cands_a:
        for cb in cands_b:
            dt_peak = float(ca["gps_peak"] - cb["gps_peak"])
            gap = float(max(0.0, max(ca["gps_start"], cb["gps_start"]) - min(ca["gps_end"], cb["gps_end"])))
            if abs(dt_peak) <= max_dt_s or gap <= max_dt_s:
                pairs.append(
                    CoincidentPair(
                        det_a=str(ca["detector"]),
                        det_b=str(cb["detector"]),
                        gps_peak_a=float(ca["gps_peak"]),
                        gps_peak_b=float(cb["gps_peak"]),
                        stft_dt_a_minus_b_s=dt_peak,
                        interval_gap_s=gap,
                        score_a=float(ca["anomaly_score"]),
                        score_b=float(cb["anomaly_score"]),
                    )
                )
    return pairs


def run_time_slide_analysis(
    bg_cands_a: list[dict[str, Any]],
    bg_cands_b: list[dict[str, Any]],
    single_epoch_livetime_s: float,
    max_dt_s: float = 0.05,
    n_slides: int = 500,
    slide_step_s: float = 4.0,
) -> dict[str, Any]:
    """Run strengthened unphysical time-slide background coincidence analysis.

    Explicitly bounds statistical significance by the resolution supported by
    the available number of independent slides (`min_resolvable_p = 1/(N+1)`).
    """
    half = n_slides // 2
    shifts = [k * slide_step_s for k in range(-half, half + 1) if k != 0][:n_slides]
    counts: list[int] = []

    for dt_shift in shifts:
        shifted_b = [
            {
                **c,
                "gps_start": float(c["gps_start"]) + dt_shift,
                "gps_end": float(c["gps_end"]) + dt_shift,
                "gps_peak": float(c["gps_peak"]) + dt_shift,
            }
            for c in bg_cands_b
        ]
        pairs = match_coincident_candidates(bg_cands_a, shifted_b, max_dt_s=max_dt_s)
        counts.append(len(pairs))

    n_used = len(shifts)
    total_coinc = int(sum(counts))
    slides_with_coinc = int(sum(1 for c in counts if c > 0))
    equiv_livetime_hr = (n_used * single_epoch_livetime_s) / 3600.0

    return {
        "n_slides": n_used,
        "slide_step_s": slide_step_s,
        "coincidence_window_s": max_dt_s,
        "n_bg_cands_a": len(bg_cands_a),
        "n_bg_cands_b": len(bg_cands_b),
        "single_epoch_livetime_hours": single_epoch_livetime_s / 3600.0,
        "equivalent_time_slide_livetime_hours": equiv_livetime_hr,
        "total_accidental_coincidences": total_coinc,
        "mean_coincidences_per_slide": float(np.mean(counts)) if counts else 0.0,
        "fraction_slides_with_ge1_coinc": float(slides_with_coinc / n_used) if n_used else 0.0,
        "coincidence_far_per_hour": (total_coinc / equiv_livetime_hr) if equiv_livetime_hr > 0 else 0.0,
        "min_resolvable_p_value": 1.0 / (n_used + 1.0),
        "significance_limitation": (
            f"Empirical coincidence significance is strictly bounded by available time-slide "
            f"resolution: p >= {1.0 / (n_used + 1.0):.5f} ({n_used} slides, "
            f"{equiv_livetime_hr:.1f} equiv. hours)."
        ),
    }
