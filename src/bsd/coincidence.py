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

    delay_s: float | None  # Estimated t_A - t_B (None when delay is unidentifiable, e.g., zero correlation)
    abs_peak_corr: float
    signed_peak_corr: float
    within_10ms_light_travel: bool
    band_hz: tuple[float, float]
    max_lag_s: float
    lags_s: np.ndarray | None = None
    corr_curve: np.ndarray | None = None
    delay_identifiable: bool = True

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "delay_s": self.delay_s,
            "delay_ms": (self.delay_s * 1000.0) if self.delay_s is not None else None,
            "abs_peak_corr": self.abs_peak_corr,
            "signed_peak_corr": self.signed_peak_corr,
            "delay_identifiable": self.delay_identifiable,
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
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("sig_a and sig_b must contain only finite values.")
    if not np.isfinite(sample_rate) or sample_rate <= 0.0:
        raise ValueError(f"sample_rate must be finite and positive, got {sample_rate}.")
    if not np.isfinite(max_lag_s) or max_lag_s <= 0.0:
        raise ValueError(f"max_lag_s must be finite and positive, got {max_lag_s}.")

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
    delay_identifiable = bool(np.isfinite(abs_peak) and abs_peak > 0.0)

    if not delay_identifiable:
        lags_s = lags_samples / sample_rate
        return CrossCorrelationResult(
            delay_s=None,
            abs_peak_corr=0.0,
            signed_peak_corr=0.0,
            within_10ms_light_travel=False,
            band_hz=(f_lo, f_hi),
            max_lag_s=max_lag_s,
            lags_s=lags_s if return_curve else None,
            corr_curve=corr if return_curve else None,
            delay_identifiable=False,
        )

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
        delay_identifiable=True,
    )


def match_coincident_candidates(
    cands_a: list[dict[str, Any]],
    cands_b: list[dict[str, Any]],
    max_dt_s: float = 0.05,
    allow_interval_gap: bool = False,
) -> list[CoincidentPair]:
    """Match candidate pairs across two detectors within ``max_dt_s`` of peak time.

    By default (`allow_interval_gap=False`), requires `|gps_peak_a - gps_peak_b| <= max_dt_s`
    so broad overlapping clusters with far-apart peaks are never falsely matched.
    """
    if not np.isfinite(max_dt_s) or max_dt_s < 0.0:
        raise ValueError(f"max_dt_s must be finite and >= 0, got {max_dt_s}.")
    pairs: list[CoincidentPair] = []
    for ca in cands_a:
        for cb in cands_b:
            pa, pb = float(ca["gps_peak"]), float(cb["gps_peak"])
            sa, ea = float(ca["gps_start"]), float(ca["gps_end"])
            sb, eb = float(cb["gps_start"]), float(cb["gps_end"])
            if not all(np.isfinite(v) for v in (pa, pb, sa, ea, sb, eb)):
                raise ValueError("Candidate timestamps must be finite.")
            dt_peak = pa - pb
            gap = float(max(0.0, max(sa, sb) - min(ea, eb)))
            matched = (abs(dt_peak) <= max_dt_s) or (allow_interval_gap and gap <= max_dt_s)
            if matched:
                pairs.append(
                    CoincidentPair(
                        det_a=str(ca["detector"]),
                        det_b=str(cb["detector"]),
                        gps_peak_a=pa,
                        gps_peak_b=pb,
                        stft_dt_a_minus_b_s=dt_peak,
                        interval_gap_s=gap,
                        score_a=float(ca["anomaly_score"]),
                        score_b=float(cb["anomaly_score"]),
                    )
                )
    return pairs


def _normalize_active_intervals(
    active_intervals: list[tuple[float, float]],
) -> list[tuple[float, float]]:
    """Validate and merge overlapping or duplicate (start, end) intervals."""
    cleaned: list[tuple[float, float]] = []
    for item in active_intervals:
        a0, a1 = float(item[0]), float(item[1])
        if not (np.isfinite(a0) and np.isfinite(a1)):
            raise ValueError(f"Active interval bounds must be finite, got ({a0}, {a1}).")
        if a1 <= a0:
            raise ValueError(f"Active interval must satisfy start < end, got ({a0}, {a1}).")
        cleaned.append((a0, a1))
    if not cleaned:
        return []
    cleaned.sort(key=lambda pair: (pair[0], pair[1]))
    merged: list[tuple[float, float]] = [cleaned[0]]
    for a0, a1 in cleaned[1:]:
        prev0, prev1 = merged[-1]
        if a0 <= prev1:
            merged[-1] = (prev0, max(prev1, a1))
        else:
            merged.append((a0, a1))
    return merged


def _interval_list_overlap_s(
    intervals_a: list[tuple[float, float]],
    intervals_b: list[tuple[float, float]],
) -> float:
    total = 0.0
    for a0, a1 in intervals_a:
        for b0, b1 in intervals_b:
            total += max(0.0, min(a1, b1) - max(a0, b0))
    return float(total)


def deduplicate_candidates_by_window(
    cands: list[dict[str, Any]],
    window_s: float = 4.0,
) -> list[dict[str, Any]]:
    """Public window-level event deduplication (one highest-scoring cluster per window/event).

    Uses the original window identifier (`window_id` or `window_gps_start`) whenever
    present so misaligned window starts are NEVER split by an arbitrary GPS grid.
    Falls back to greedy peak-proximity clustering within `window_s` when no window key exists.
    """
    return _deduplicate_candidates_by_window(cands, window_s)


def _deduplicate_candidates_by_window(
    cands: list[dict[str, Any]],
    window_s: float,
) -> list[dict[str, Any]]:
    """Deduplicate multi-scale clusters to one highest-scoring candidate per window/event.

    Uses the original window identifier (`window_id` or `window_gps_start`) whenever
    present so misaligned window starts are NEVER split by an arbitrary GPS grid.
    Falls back to greedy peak-proximity clustering within `window_s` when no window key exists.
    """
    if not cands:
        return []
    has_win_key = all(("window_id" in c or "window_gps_start" in c) for c in cands)
    if has_win_key:
        best_by_key: dict[Any, dict[str, Any]] = {}
        for c in cands:
            gp = float(c["gps_peak"])
            if not np.isfinite(gp):
                raise ValueError("Candidate gps_peak must be finite.")
            key = c["window_id"] if "window_id" in c else round(float(c["window_gps_start"]), 6)
            if key not in best_by_key or float(c["anomaly_score"]) > float(best_by_key[key]["anomaly_score"]):
                best_by_key[key] = c
        return list(best_by_key.values())

    if window_s <= 0.0:
        return list(cands)
    # Greedy peak-proximity deduplication (highest score first) without absolute GPS grid bias
    ordered = sorted(cands, key=lambda c: float(c["anomaly_score"]), reverse=True)
    kept: list[dict[str, Any]] = []
    for c in ordered:
        gp = float(c["gps_peak"])
        if not np.isfinite(gp):
            raise ValueError("Candidate gps_peak must be finite.")
        if all(abs(gp - float(k["gps_peak"])) >= window_s for k in kept):
            kept.append(c)
    return kept


def run_time_slide_analysis(
    bg_cands_a: list[dict[str, Any]],
    bg_cands_b: list[dict[str, Any]],
    single_epoch_livetime_s: float,
    max_dt_s: float = 0.05,
    n_slides: int = 500,
    slide_step_s: float = 4.0,
    active_intervals: list[tuple[float, float]] | None = None,
    deduplicate_window_s: float | None = None,
) -> dict[str, Any]:
    """Run unphysical time-slide background coincidence analysis.

    - Validates and merges overlapping/duplicate `active_intervals` to compute
      true overlapping exposure per slide.
    - When `active_intervals` is None, or when actual overlapping exposure or
      the number of performed slides is zero, `coincidence_far_per_hour` is
      returned as `None` (UNAVAILABLE), never `0.0`.
    - Deduplicates by original window identifier (`window_id` / `window_gps_start`)
      or peak proximity when `deduplicate_window_s` is set.
    """
    if not np.isfinite(single_epoch_livetime_s) or single_epoch_livetime_s <= 0.0:
        raise ValueError(f"single_epoch_livetime_s must be finite and positive, got {single_epoch_livetime_s}.")
    if n_slides < 0:
        raise ValueError(f"n_slides must be >= 0, got {n_slides}.")
    if not np.isfinite(slide_step_s) or slide_step_s <= 0.0:
        raise ValueError(f"slide_step_s must be finite and positive, got {slide_step_s}.")

    merged_intervals: list[tuple[float, float]] | None = (
        _normalize_active_intervals(active_intervals) if active_intervals is not None else None
    )

    eff_a = (
        _deduplicate_candidates_by_window(bg_cands_a, deduplicate_window_s)
        if deduplicate_window_s is not None
        else list(bg_cands_a)
    )
    eff_b = (
        _deduplicate_candidates_by_window(bg_cands_b, deduplicate_window_s)
        if deduplicate_window_s is not None
        else list(bg_cands_b)
    )

    half = n_slides // 2
    shifts = [k * slide_step_s for k in range(-half, half + 1) if k != 0][:n_slides]
    counts: list[int] = []
    overlap_s_list: list[float] = []

    for dt_shift in shifts:
        if merged_intervals is not None:
            shifted_intervals = [(a0 + dt_shift, a1 + dt_shift) for a0, a1 in merged_intervals]
            ov_s = _interval_list_overlap_s(merged_intervals, shifted_intervals)
            overlap_s_list.append(ov_s)
            shifted_b = []
            for c in eff_b:
                gp_s = float(c["gps_peak"]) + dt_shift
                if any(s0 <= gp_s <= s1 for s0, s1 in merged_intervals):
                    shifted_b.append(
                        {
                            **c,
                            "gps_start": float(c["gps_start"]) + dt_shift,
                            "gps_end": float(c["gps_end"]) + dt_shift,
                            "gps_peak": gp_s,
                        }
                    )
        else:
            overlap_s_list.append(single_epoch_livetime_s)
            shifted_b = [
                {
                    **c,
                    "gps_start": float(c["gps_start"]) + dt_shift,
                    "gps_end": float(c["gps_end"]) + dt_shift,
                    "gps_peak": float(c["gps_peak"]) + dt_shift,
                }
                for c in eff_b
            ]
        pairs = match_coincident_candidates(eff_a, shifted_b, max_dt_s=max_dt_s, allow_interval_gap=False)
        counts.append(len(pairs))

    n_used = len(shifts)
    total_coinc = int(sum(counts))
    slides_with_coinc = int(sum(1 for c in counts if c > 0))
    nominal_livetime_hr = (n_used * single_epoch_livetime_s) / 3600.0
    equiv_livetime_hr = float(sum(overlap_s_list) / 3600.0) if overlap_s_list else 0.0

    if n_used == 0:
        far_per_hr: float | None = None
        far_status = "UNAVAILABLE_ZERO_SLIDES"
    elif merged_intervals is None:
        far_per_hr = None
        far_status = "UNAVAILABLE_NO_ACTIVE_INTERVALS"
    elif equiv_livetime_hr <= 0.0:
        far_per_hr = None
        far_status = "UNAVAILABLE_ZERO_EXPOSURE"
    else:
        far_per_hr = float(total_coinc / equiv_livetime_hr)
        far_status = "AVAILABLE"

    tech_slide_res = (1.0 / (n_used + 1.0)) if n_used > 0 else None
    sig_available = bool(far_status == "AVAILABLE" and n_used > 0)
    return {
        "n_slides": n_used,
        "n_slides_performed": n_used,
        "n_slides_requested": int(n_slides),
        "slide_step_s": slide_step_s,
        "coincidence_window_s": max_dt_s,
        "n_bg_cands_a": len(eff_a),
        "n_bg_cands_b": len(eff_b),
        "deduplicated_window_s": deduplicate_window_s,
        "overlap_corrected": bool(merged_intervals is not None),
        "single_epoch_livetime_hours": single_epoch_livetime_s / 3600.0,
        "nominal_time_slide_livetime_hours": nominal_livetime_hr,
        "equivalent_time_slide_livetime_hours": equiv_livetime_hr,
        "total_accidental_coincidences": total_coinc,
        "mean_coincidences_per_slide": float(np.mean(counts)) if counts else 0.0,
        "fraction_slides_with_ge1_coinc": float(slides_with_coinc / n_used) if n_used else None,
        "coincidence_far_per_hour": far_per_hr,
        "far_status": far_status,
        "significance_status": "AVAILABLE" if sig_available else far_status,
        "technical_slide_count_resolution": tech_slide_res,
        "min_resolvable_p_value": tech_slide_res if sig_available else None,
        "significance_limitation": (
            f"Performed {n_used} slide(s); overlap_corrected={merged_intervals is not None}; "
            f"equiv_livetime={equiv_livetime_hr:.2f} hr; far_status={far_status}; "
            f"significance_status={'AVAILABLE' if sig_available else far_status}. "
            "Non-circular slides with diminishing overlap are not equal-exposure i.i.d. trials."
        ),
    }
