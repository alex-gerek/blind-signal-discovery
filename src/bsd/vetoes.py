"""BSD v0.2 Phase 2 — Auditable & Reversible Instrumental Artifact Diagnostics.

Provides two independent, physically interpretable artifact diagnostics:
  1. Impulsive glitch diagnostic (`IMPULSIVE_GLITCH_*`): identifies ultra-short
     (1-3 sample) high-amplitude spikes and extreme non-Gaussian local kurtosis
     in the time-domain waveform around a candidate.
  2. Persistent narrowband line diagnostic (`PERSISTENT_NARROWBAND_LINE`):
     identifies candidates spanning most of the observation window within a
     narrow frequency band (configurable: disabled for general searches where
     stationary tones are targets, enabled for transient searches).

Raw-Candidate Preservation Principle:
  Every raw candidate (score, timestamps, frequency range) is preserved
  verbatim. Veto decisions are attached as reversible annotations with explicit
  reason codes and diagnostic measurements so pre-veto and post-veto performance
  can always be compared side-by-side.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from scipy import stats

from bsd.detector import CandidateEvent
from bsd.detector_v2 import MultiResolutionCandidate
from bsd.generator import Observation

CandidateLike = CandidateEvent | MultiResolutionCandidate

REASON_GLITCH_KURTOSIS = "IMPULSIVE_GLITCH_EXTREME_KURTOSIS"
REASON_GLITCH_SHORT_SPIKE = "IMPULSIVE_GLITCH_ULTRA_SHORT_SPIKE"
REASON_PERSISTENT_LINE = "PERSISTENT_NARROWBAND_LINE"


@dataclass(frozen=True)
class VetoConfig:
    """Configuration for reversible instrumental artifact vetoes."""

    enable_glitch_veto: bool = True
    glitch_kurtosis_threshold: float = 20.0
    glitch_spike_sigma_threshold: float = 6.0
    glitch_max_half_peak_samples: int = 3
    enable_line_veto: bool = False
    line_min_time_coverage: float = 0.75
    line_max_bandwidth_hz: float = 48.0
    local_context_s: float = 0.08

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class VetoDecision:
    """Reversible veto decision attached to a preserved raw candidate."""

    candidate: CandidateLike
    vetoed: bool
    reason_codes: tuple[str, ...]
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate": self.candidate.to_dict(),
            "vetoed": self.vetoed,
            "reason_codes": list(self.reason_codes),
            "diagnostics": dict(self.diagnostics),
        }


@dataclass(frozen=True)
class VetoAssessment:
    """Complete pre-veto and post-veto audit record for one Observation."""

    sample_id: str
    raw_candidates: tuple[CandidateLike, ...]
    decisions: tuple[VetoDecision, ...]
    surviving_candidates: tuple[CandidateLike, ...]
    detected_before_veto: bool
    detected_after_veto: bool
    n_vetoed_glitch: int
    n_vetoed_line: int
    config: VetoConfig

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "detected_before_veto": self.detected_before_veto,
            "detected_after_veto": self.detected_after_veto,
            "n_raw_candidates": len(self.raw_candidates),
            "n_surviving_candidates": len(self.surviving_candidates),
            "n_vetoed_glitch": self.n_vetoed_glitch,
            "n_vetoed_line": self.n_vetoed_line,
            "decisions": [d.to_dict() for d in self.decisions],
            "config": self.config.to_dict(),
        }


def _max_contiguous_true(mask: np.ndarray) -> int:
    if mask.size == 0 or not np.any(mask):
        return 0
    padded = np.concatenate(([False], mask, [False]))
    diffs = np.diff(padded.astype(np.int8))
    starts = np.where(diffs == 1)[0]
    ends = np.where(diffs == -1)[0]
    return int(np.max(ends - starts)) if starts.size else 0


def evaluate_impulsive_glitch(
    observation: Observation,
    candidate: CandidateLike,
    config: VetoConfig,
) -> tuple[list[str], dict[str, Any]]:
    """Evaluate time-domain impulsive-glitch diagnostics around a candidate."""
    x_full = np.asarray(observation.signal, dtype=np.float64)
    x_centered = x_full - float(np.median(x_full))
    global_mad_sigma = max(1.4826 * float(np.median(np.abs(x_centered))), 1e-18)

    fs = observation.sample_rate
    i0 = max(0, int(np.floor((candidate.t_start - config.local_context_s) * fs)))
    i1 = min(x_centered.size, int(np.ceil((candidate.t_end + config.local_context_s) * fs)))
    if i1 - i0 < 16:
        i0 = 0
        i1 = x_centered.size

    seg = x_centered[i0:i1]
    abs_seg = np.abs(seg)
    peak_amp = float(np.max(abs_seg)) if seg.size else 0.0
    peak_mad_sigma = peak_amp / global_mad_sigma
    local_kurtosis = float(stats.kurtosis(seg, fisher=True)) if seg.size >= 8 else 0.0

    # Count samples reaching >= 50% of local peak amplitude (1-3 for Dirac-like spikes)
    half_peak_mask = abs_seg >= (0.5 * peak_amp) if peak_amp > 0 else np.zeros_like(abs_seg, dtype=bool)
    n_ge_half_peak = int(np.sum(half_peak_mask))
    max_contig_5sigma = _max_contiguous_true(abs_seg >= (5.0 * global_mad_sigma))

    reasons: list[str] = []
    if config.enable_glitch_veto:
        if local_kurtosis >= config.glitch_kurtosis_threshold:
            reasons.append(REASON_GLITCH_KURTOSIS)
        if (
            peak_mad_sigma >= config.glitch_spike_sigma_threshold
            and n_ge_half_peak <= config.glitch_max_half_peak_samples
        ):
            reasons.append(REASON_GLITCH_SHORT_SPIKE)

    diag = {
        "local_kurtosis": local_kurtosis,
        "local_peak_mad_sigma": peak_mad_sigma,
        "n_samples_ge_half_peak": n_ge_half_peak,
        "max_contiguous_ge_5sigma": max_contig_5sigma,
        "local_slice_samples": int(seg.size),
    }
    return reasons, diag


def evaluate_persistent_line(
    observation: Observation,
    candidate: CandidateLike,
    config: VetoConfig,
) -> tuple[list[str], dict[str, Any]]:
    """Evaluate whether a candidate is a persistent narrowband instrumental line."""
    cand_dur = max(0.0, float(candidate.t_end - candidate.t_start))
    time_cov = cand_dur / max(float(observation.duration), 1e-9)
    bw_hz = max(0.0, float(candidate.f_max - candidate.f_min))

    reasons: list[str] = []
    if config.enable_line_veto:
        if time_cov >= config.line_min_time_coverage and bw_hz <= config.line_max_bandwidth_hz:
            reasons.append(REASON_PERSISTENT_LINE)

    diag = {
        "time_coverage_fraction": time_cov,
        "bandwidth_hz": bw_hz,
    }
    return reasons, diag


def apply_vetoes(
    observation: Observation,
    candidates: list[CandidateLike] | tuple[CandidateLike, ...],
    config: VetoConfig | None = None,
) -> VetoAssessment:
    """Apply auditable, reversible vetoes while preserving all raw candidates."""
    if not isinstance(observation, Observation):
        raise TypeError(
            f"apply_vetoes strictly requires an Observation, got {type(observation).__name__}."
        )
    cfg = config or VetoConfig()
    raw_tuple = tuple(candidates)
    decisions: list[VetoDecision] = []
    surviving: list[CandidateLike] = []
    n_glitch = 0
    n_line = 0

    for cand in raw_tuple:
        g_reasons, g_diag = evaluate_impulsive_glitch(observation, cand, cfg)
        l_reasons, l_diag = evaluate_persistent_line(observation, cand, cfg)
        all_reasons = tuple(g_reasons + l_reasons)
        vetoed = len(all_reasons) > 0
        if g_reasons:
            n_glitch += 1
        if l_reasons:
            n_line += 1
        decisions.append(
            VetoDecision(
                candidate=cand,
                vetoed=vetoed,
                reason_codes=all_reasons,
                diagnostics={**g_diag, **l_diag},
            )
        )
        if not vetoed:
            surviving.append(cand)

    return VetoAssessment(
        sample_id=observation.sample_id,
        raw_candidates=raw_tuple,
        decisions=tuple(decisions),
        surviving_candidates=tuple(surviving),
        detected_before_veto=len(raw_tuple) > 0,
        detected_after_veto=len(surviving) > 0,
        n_vetoed_glitch=n_glitch,
        n_vetoed_line=n_line,
        config=cfg,
    )
