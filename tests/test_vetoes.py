"""Automated pytest suite for BSD v0.2 Phase 2 (auditable & reversible vetoes)."""

from __future__ import annotations

import numpy as np
import pytest

from bsd.detector import CandidateEvent
from bsd.generator import Observation, generate_sample, make_split_rng
from bsd.vetoes import (
    REASON_GLITCH_KURTOSIS,
    REASON_GLITCH_SHORT_SPIKE,
    REASON_PERSISTENT_LINE,
    VetoConfig,
    apply_vetoes,
)


def _dummy_cand(t0: float, t1: float, f0: float, f1: float, score: float = 120.0) -> CandidateEvent:
    return CandidateEvent(
        candidate_id=1,
        t_start=t0,
        t_end=t1,
        f_min=f0,
        f_max=f1,
        peak_time=0.5 * (t0 + t1),
        peak_freq=0.5 * (f0 + f1),
        anomaly_score=score,
        peak_excess_power=15.0,
        n_bins=12,
        empirical_p_value=0.002,
        p_value_is_upper_bound=True,
        diagnostics={},
    )


def test_raw_candidates_preserved_and_veto_reversible() -> None:
    """Raw candidates must never be altered or deleted; disabling vetoes restores them."""
    rng = np.random.default_rng(7)
    x = rng.normal(0.0, 1.0, size=4096)
    x[2048] += 15.0  # 1-sample impulsive glitch at t=2.0 s
    obs = Observation("glitch_sample", x, 1024.0, 4.0)
    cand = _dummy_cand(1.95, 2.05, 20.0, 480.0, score=250.0)

    assess_on = apply_vetoes(obs, [cand], VetoConfig(enable_glitch_veto=True, enable_line_veto=False))
    assert assess_on.detected_before_veto is True
    assert assess_on.detected_after_veto is False
    assert assess_on.raw_candidates == (cand,)
    assert len(assess_on.surviving_candidates) == 0
    assert assess_on.n_vetoed_glitch == 1
    assert any(
        r in (REASON_GLITCH_KURTOSIS, REASON_GLITCH_SHORT_SPIKE)
        for r in assess_on.decisions[0].reason_codes
    )

    # Reversibility: disabling glitch veto recovers exact same candidate
    assess_off = apply_vetoes(obs, [cand], VetoConfig(enable_glitch_veto=False, enable_line_veto=False))
    assert assess_off.detected_after_veto is True
    assert assess_off.surviving_candidates == (cand,)
    assert assess_off.decisions[0].reason_codes == ()


def test_genuine_linear_chirp_not_vetoed() -> None:
    """A genuine structured linear chirp must survive both glitch and line vetoes."""
    rng = make_split_rng(42, "veto_chirp_test")
    obs, gt = generate_sample("chirp_ok", 1024.0, 4.0, "gaussian_white", rng, "linear_chirp", 18.0)
    assert gt.t_start is not None and gt.t_end is not None and gt.f_min is not None and gt.f_max is not None
    cand = _dummy_cand(gt.t_start, gt.t_end, gt.f_min, gt.f_max, score=300.0)

    assess = apply_vetoes(obs, [cand], VetoConfig(enable_glitch_veto=True, enable_line_veto=True))
    assert assess.detected_after_veto is True
    assert assess.n_vetoed_glitch == 0
    assert assess.n_vetoed_line == 0


def test_persistent_narrowband_line_veto_is_configurable() -> None:
    """Persistent line veto flags full-window narrowband clusters only when enabled."""
    rng = np.random.default_rng(11)
    obs = Observation("line_obs", rng.normal(0.0, 1.0, size=4096), 1024.0, 4.0)
    line_cand = _dummy_cand(0.05, 3.95, 18.0, 34.0, score=4500.0)

    assess_transient_mode = apply_vetoes(
        obs, [line_cand], VetoConfig(enable_glitch_veto=True, enable_line_veto=True)
    )
    assert assess_transient_mode.detected_before_veto is True
    assert assess_transient_mode.detected_after_veto is False
    assert assess_transient_mode.n_vetoed_line == 1
    assert REASON_PERSISTENT_LINE in assess_transient_mode.decisions[0].reason_codes

    assess_general_mode = apply_vetoes(
        obs, [line_cand], VetoConfig(enable_glitch_veto=True, enable_line_veto=False)
    )
    assert assess_general_mode.detected_after_veto is True
    assert assess_general_mode.n_vetoed_line == 0


def test_apply_vetoes_rejects_non_observation() -> None:
    cand = _dummy_cand(1.0, 1.2, 50.0, 150.0)
    with pytest.raises(TypeError):
        apply_vetoes(np.zeros(4096), [cand])  # type: ignore[arg-type]
