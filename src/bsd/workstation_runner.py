"""Headless runner for the BSD Local Measurement Workstation.

Wraps existing BSD v0.3 Phase 2b modules (`AdaptiveRegimeManagerV2b`,
`calibrate_separate_search_modes`, `FrozenReferenceBroadbandChannel`,
`MultiResolutionDetector`, `VetoConfig`) WITHOUT rewriting detector mathematics.

Enforces:
  - Explicit 3-tier calibration status: INSUFFICIENT, PROVISIONAL, VALIDATED
    (never assigning VALIDATED from raw window count alone).
  - Strict scientific output semantics separating measured features, statistical
    anomalies, detection candidates, and independent physical event confirmation.
  - Reproducible synthetic self-test with positive ground-truth injections AND
    pure-noise negative controls, reporting both detected and missed injections.
  - Complete reproducible CSV/JSON export with run metadata & input provenance.
"""

from __future__ import annotations

import csv
import json
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from bsd.adaptive_regime_v2b import (
    AdaptiveRegimeManagerV2b,
    AdaptiveRegimeV2bConfig,
    CalibrationUncertaintySummary,
    DualOutputWindowV2b,
    FrozenReferenceBroadbandChannel,
    calibrate_separate_search_modes,
    fit_glitch_reference,
)
from bsd.calibration import conformal_order_rank
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, _scale_to_energy_snr, generate_clean_signal, make_split_rng
from bsd.io_loader import (
    IndependenceAudit,
    InputValidationError,
    LoadedTimeSeries,
    PreprocessingSummary,
    check_background_independence,
    slice_into_observations,
)
from bsd.noise import generate_noise
from bsd.vetoes import VetoConfig

CAL_STATUS_INSUFFICIENT = "INSUFFICIENT"
CAL_STATUS_PROVISIONAL = "PROVISIONAL"
CAL_STATUS_VALIDATED = "VALIDATED"

MIN_BG_WINDOWS_PROVISIONAL = 60
MIN_THRESH_VALIDATED = 150
MIN_NEFF_VALIDATED = 120.0
MAX_UPPER_FAR_CI_VALIDATED = 0.025


@dataclass(frozen=True)
class WorkstationConfig:
    """User-configurable parameters for a workstation measurement run."""

    search_mode: str = "transient"  # "transient" or "persistent"
    processing_mode: str = "causal_trailing"  # "causal_trailing" or "offline_symmetric"
    window_duration_s: float = 4.0
    target_fap: float = 0.01
    enable_glitch_veto: bool = False  # Off by default for arbitrary non-LIGO sensors
    enable_line_veto: bool = False  # Off by default for arbitrary non-LIGO sensors
    calibration_seed: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CalibrationStatusAudit:
    """Explicit 3-tier calibration assessment (INSUFFICIENT / PROVISIONAL / VALIDATED)."""

    status: str
    n_bg_windows: int
    n_thresh: int
    n_eff: float
    lag1_autocorr: float
    ljung_box_p: float
    p_min: float
    expected_far_iid: float
    far_ci90_iid: tuple[float, float]
    far_ci90_neff: tuple[float, float]
    dependence_detected: bool
    unmet_validation_conditions: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class WorkstationRunResult:
    """Complete auditable output of a workstation run."""

    run_metadata: dict[str, Any]
    scientific_semantics: dict[str, str]
    independence_audit: IndependenceAudit
    calibration_audit: CalibrationStatusAudit
    summary_metrics: dict[str, Any]
    windows: list[DualOutputWindowV2b]
    target_observations: list[Observation]
    candidates_rows: list[dict[str, Any]]
    timeline_rows: list[dict[str, Any]]
    ground_truth_audit: dict[str, Any] | None = None


def _get_git_commit() -> str:
    try:
        root = Path(__file__).resolve().parent.parent.parent
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=root, stderr=subprocess.DEVNULL)
        return out.decode("utf-8").strip()
    except Exception:
        return "unknown"


def _wilson_ci90(k: int, n: int) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 1.0)
    z = 1.6448536
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = (z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)) / denom
    return (max(0.0, float(center - half)), min(1.0, float(center + half)))


def classify_calibration_status(
    n_bg_windows: int,
    uncertainty: CalibrationUncertaintySummary,
    target_fap: float | None = None,
) -> CalibrationStatusAudit:
    """Assign INSUFFICIENT, PROVISIONAL, or VALIDATED using explicit statistical criteria.

    Note: VALIDATED indicates that internal finite-sample and linear lag-1..10
    screening criteria were met on the supplied background ensemble. It does NOT
    prove stationarity, exchangeability, or real-world false-alarm control.
    """
    unmet: list[str] = []
    fap_unattainable = bool(target_fap is not None and np.isfinite(target_fap) and uncertainty.p_min > target_fap)
    if fap_unattainable and target_fap is not None:
        unmet.append(
            f"target_fap={target_fap:.4g} unattainable with n_thresh={uncertainty.n_thresh} "
            f"(p_min={uncertainty.p_min:.5f} > target_fap)"
        )
    if n_bg_windows < MIN_BG_WINDOWS_PROVISIONAL:
        unmet.append(f"n_bg_windows={n_bg_windows} < {MIN_BG_WINDOWS_PROVISIONAL}")
    if uncertainty.n_eff < 45.0:
        unmet.append(f"n_eff={uncertainty.n_eff:.1f} < 45.0 (severe dependence or insufficient sample size)")

    if fap_unattainable or n_bg_windows < MIN_BG_WINDOWS_PROVISIONAL or uncertainty.n_eff < 45.0:
        tier = CAL_STATUS_INSUFFICIENT
    else:
        if uncertainty.n_thresh < MIN_THRESH_VALIDATED:
            unmet.append(f"n_thresh={uncertainty.n_thresh} < {MIN_THRESH_VALIDATED}")
        if uncertainty.n_eff < MIN_NEFF_VALIDATED:
            unmet.append(f"n_eff={uncertainty.n_eff:.1f} < {MIN_NEFF_VALIDATED:.0f}")
        if uncertainty.ljung_box_p < 0.05:
            unmet.append(f"temporal dependence detected (Ljung-Box p={uncertainty.ljung_box_p:.4f} < 0.05)")
        if uncertainty.far_ci90_neff[1] > MAX_UPPER_FAR_CI_VALIDATED:
            unmet.append(
                f"dependence-aware upper 90% FAR CI={uncertainty.far_ci90_neff[1]*100:.2f}% > "
                f"{MAX_UPPER_FAR_CI_VALIDATED*100:.1f}%"
            )
        tier = CAL_STATUS_VALIDATED if not unmet else CAL_STATUS_PROVISIONAL

    return CalibrationStatusAudit(
        status=tier,
        n_bg_windows=n_bg_windows,
        n_thresh=uncertainty.n_thresh,
        n_eff=uncertainty.n_eff,
        lag1_autocorr=uncertainty.lag1_autocorr,
        ljung_box_p=uncertainty.ljung_box_p,
        p_min=uncertainty.p_min,
        expected_far_iid=uncertainty.expected_far_iid,
        far_ci90_iid=uncertainty.far_ci90_iid,
        far_ci90_neff=uncertainty.far_ci90_neff,
        dependence_detected=uncertainty.dependence_detected,
        unmet_validation_conditions=tuple(unmet),
    )


def _scientific_semantics_block() -> dict[str, str]:
    return {
        "measured_features": (
            "Raw time-domain waveform, 3-scale STFT excess power, 512-pt spectral-continuum band powers, "
            "and 31-ms broadband impulsive crest ratio."
        ),
        "statistical_anomalies": (
            "Windows whose raw STFT or frozen-reference broadband score exceeds the empirical calibration "
            "threshold relative to the supplied background dataset."
        ),
        "detection_candidates": (
            "Localized time-frequency clusters nominated in Stream A, preserved 100% verbatim alongside "
            "reversible veto annotations and orthogonal Stream B credibility warnings."
        ),
        "independent_physical_confirmation": (
            "NOT_ESTABLISHED_BY_SINGLE_CHANNEL_WORKSTATION: Statistical anomalies and detection candidates "
            "are never proof of a physical or astrophysical event without independent multi-sensor or "
            "domain-specific confirmation."
        ),
    }


def _build_rows(
    windows: list[DualOutputWindowV2b],
    window_duration_s: float,
    search_mode: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    timeline_rows: list[dict[str, Any]] = []
    candidates_rows: list[dict[str, Any]] = []

    for w in windows:
        t0_win = w.index * window_duration_s
        diag = w.diagnostics
        active_diag_nom = bool(w.stft_nominated_post_veto if search_mode == "transient" else w.broadband_nominated)
        active_decision = bool(w.detection_available and active_diag_nom)

        ded_stft_avail = bool(
            diag.get("dedicated_stft_available", w.stft_threshold is not None and np.isfinite(w.stft_threshold))
        )
        ded_bb_avail = bool(
            diag.get(
                "dedicated_broadband_available",
                w.broadband_threshold is not None and np.isfinite(w.broadband_threshold),
            )
        )
        u_stft_avail = bool(diag.get("union_stft_available", ded_stft_avail))
        u_bb_avail = bool(diag.get("union_broadband_available", ded_bb_avail))
        u_stft_th_raw = diag.get("union_stft_threshold", None)
        u_stft_th = (
            round(float(u_stft_th_raw), 4)
            if (u_stft_th_raw is not None and np.isfinite(u_stft_th_raw))
            else None
        )
        if "union_mode_status" in diag:
            u_mode_status = str(diag["union_mode_status"])
        elif u_stft_avail and u_bb_avail:
            u_mode_status = "FULL_UNION_AVAILABLE"
        elif u_bb_avail:
            u_mode_status = "PARTIAL_BROADBAND_ONLY"
        elif u_stft_avail:
            u_mode_status = "PARTIAL_STFT_ONLY"
        else:
            u_mode_status = "UNION_UNAVAILABLE"

        rl_stft = bool(diag.get("resource_limit_blocked_stft", False))
        rl_union = bool(diag.get("resource_limit_blocked_union", False))
        cal_upd_next = bool(diag.get("calibration_updated_for_next_window", False))

        timeline_rows.append(
            {
                "window_index": w.index,
                "sample_id": w.sample_id,
                "t_window_start_s": round(t0_win, 4),
                "t_window_end_s": round(t0_win + window_duration_s, 4),
                "search_mode": search_mode,
                "stft_max_score": round(w.stft_max_score, 4),
                "stft_threshold": round(w.stft_threshold, 4) if w.stft_threshold is not None else None,
                "stft_nominated_raw": w.stft_nominated_raw,
                "stft_nominated_post_veto": w.stft_nominated_post_veto,
                "broadband_excess_score": round(w.broadband_excess_score, 4),
                "broadband_threshold": round(w.broadband_threshold, 4) if w.broadband_threshold is not None else None,
                "broadband_nominated": w.broadband_nominated,
                "active_mode_diagnostic_nomination": active_diag_nom,
                "active_mode_anomaly_flag": active_decision,
                "union_sidak_nominated": w.joint_nominated_post_veto,
                "dedicated_stft_available": ded_stft_avail,
                "dedicated_broadband_available": ded_bb_avail,
                "union_stft_available": u_stft_avail,
                "union_broadband_available": u_bb_avail,
                "union_stft_threshold": u_stft_th,
                "union_mode_status": u_mode_status,
                "resource_limit_blocked_stft": rl_stft,
                "resource_limit_blocked_union": rl_union,
                "calibration_updated_for_next_window": cal_upd_next,
                "regime_id": w.regime_id,
                "credibility_status": w.credibility_status,
                "continuum_status": w.continuum_status,
                "glitch_storm_status": w.glitch_storm_status,
                "warning_causes": ";".join(w.warning_causes),
                "significance_status": w.significance_status,
                "detection_available": w.detection_available,
                "is_clean_calibrated": w.is_clean_calibrated,
            }
        )

        veto_map: dict[int, str] = {}
        if w.veto_assessment is not None:
            for vc in w.veto_assessment.vetoed_candidates:
                veto_map[vc.candidate.candidate_id] = ";".join(vc.reasons)

        for c in w.raw_stft_candidates:
            v_reason = veto_map.get(c.candidate_id, "")
            candidates_rows.append(
                {
                    "candidate_type": "STFT_DETECTION_CANDIDATE",
                    "window_index": w.index,
                    "sample_id": w.sample_id,
                    "t_abs_start_s": round(t0_win + c.t_start, 4),
                    "t_abs_end_s": round(t0_win + c.t_end, 4),
                    "peak_time_abs_s": round(t0_win + c.peak_time, 4),
                    "f_min_hz": round(c.f_min, 2),
                    "f_max_hz": round(c.f_max, 2),
                    "peak_freq_hz": round(c.peak_freq, 2),
                    "scale_name": c.scale_name,
                    "raw_anomaly_score": round(c.raw_anomaly_score, 4),
                    "standardized_score": round(c.standardized_score, 4),
                    "empirical_p_value": None if c.empirical_p_value is None else round(c.empirical_p_value, 6),
                    "vetoed": bool(v_reason),
                    "veto_reasons": v_reason,
                    "credibility_status": w.credibility_status,
                    "warning_causes": ";".join(w.warning_causes),
                    "physical_confirmation": "UNCONFIRMED_STATISTICAL_CANDIDATE",
                }
            )
        if w.broadband_nominated:
            candidates_rows.append(
                {
                    "candidate_type": "BROADBAND_STATISTICAL_ANOMALY",
                    "window_index": w.index,
                    "sample_id": w.sample_id,
                    "t_abs_start_s": round(t0_win, 4),
                    "t_abs_end_s": round(t0_win + window_duration_s, 4),
                    "peak_time_abs_s": round(t0_win + 0.5 * window_duration_s, 4),
                    "f_min_hz": 20.0,
                    "f_max_hz": 450.0,
                    "peak_freq_hz": None,
                    "scale_name": "frozen_ref_continuum_512",
                    "raw_anomaly_score": round(w.broadband_excess_score, 4),
                    "standardized_score": round(w.broadband_excess_score, 4),
                    "empirical_p_value": None,
                    "vetoed": False,
                    "veto_reasons": "",
                    "credibility_status": w.credibility_status,
                    "warning_causes": ";".join(w.warning_causes),
                    "physical_confirmation": "SINGLE_CHANNEL_AMBIGUOUS_NOT_CONFIRMED",
                }
            )

    return candidates_rows, timeline_rows


def run_workstation_analysis(
    target: LoadedTimeSeries,
    background: LoadedTimeSeries,
    config: WorkstationConfig | None = None,
) -> WorkstationRunResult:
    """Execute a complete, reproducible BSD workstation measurement run."""
    cfg = config or WorkstationConfig()
    if cfg.search_mode not in ("transient", "persistent"):
        raise InputValidationError(f"Invalid search_mode: {cfg.search_mode!r}")
    if cfg.processing_mode not in ("causal_trailing", "offline_symmetric"):
        raise InputValidationError(f"Invalid processing_mode: {cfg.processing_mode!r}")
    if not np.isfinite(cfg.target_fap) or not (0.001 <= cfg.target_fap <= 0.25):
        raise InputValidationError(f"target_fap must be finite and in [0.001, 0.25], got {cfg.target_fap}.")

    indep = check_background_independence(target, background)
    bg_obs = slice_into_observations(
        background, window_duration_s=cfg.window_duration_s, prefix="bg", min_windows=MIN_BG_WINDOWS_PROVISIONAL
    )
    n_norm_split = max(20, int(round(len(bg_obs) * 0.30)))
    n_thresh_split = len(bg_obs) - n_norm_split
    _, fap_ok, min_n_thresh_req = conformal_order_rank(max(1, n_thresh_split), cfg.target_fap)
    if not fap_ok:
        min_bg_req = int(np.ceil(min_n_thresh_req / 0.70))
        raise InputValidationError(
            f"Unattainable target_fap={cfg.target_fap:.4g} for {len(bg_obs)} background windows "
            f"(Split B has n_thresh={n_thresh_split}, minimum resolvable p_min=1/({n_thresh_split}+1)="
            f"{1.0 / (n_thresh_split + 1):.5f} > target_fap). Supply at least {min_bg_req} background windows "
            f"or raise target_fap to >= {1.0 / (n_thresh_split + 1):.5f}."
        )
    tgt_obs = slice_into_observations(
        target, window_duration_s=cfg.window_duration_s, prefix="tgt", min_windows=1
    )

    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    scal = calibrate_separate_search_modes(
        det,
        bb,
        bg_obs,
        target_fap_transient=cfg.target_fap,
        target_fap_persistent=cfg.target_fap,
        target_fwer_union=cfg.target_fap,
        norm_fraction=0.30,
        calibration_seed=cfg.calibration_seed,
    )

    active_unc = (
        scal.transient_stft_uncertainty
        if cfg.search_mode == "transient"
        else scal.persistent_bb_cal.uncertainty
    )
    cal_audit = classify_calibration_status(len(bg_obs), active_unc, target_fap=cfg.target_fap)

    mgr = AdaptiveRegimeManagerV2b(
        AdaptiveRegimeV2bConfig(
            mode=cfg.processing_mode,
            policy=cfg.search_mode,
            target_fwer=cfg.target_fap,
        )
    )
    cref = mgr.fit_continuum_reference(bg_obs)
    gref = fit_glitch_reference(bg_obs)
    vcfg = (
        VetoConfig(enable_glitch_veto=cfg.enable_glitch_veto, enable_line_veto=cfg.enable_line_veto)
        if (cfg.enable_glitch_veto or cfg.enable_line_veto)
        else None
    )

    windows = mgr.process_sequence(tgt_obs, det, bb, scal, cref, gref, veto_config=vcfg)
    cands_rows, tl_rows = _build_rows(windows, cfg.window_duration_s, cfg.search_mode)

    n_win = len(windows)
    n_det_avail = sum(1 for w in windows if w.detection_available)
    n_cal_unavail = n_win - n_det_avail
    n_clean = sum(1 for w in windows if w.is_clean_calibrated)
    n_anom_active = sum(1 for r in tl_rows if r["active_mode_anomaly_flag"])
    union_status_counts: dict[str, int] = {}
    for r in tl_rows:
        st_name = str(r["union_mode_status"])
        union_status_counts[st_name] = union_status_counts.get(st_name, 0) + 1
    run_ts = datetime.now(timezone.utc)
    run_id = run_ts.strftime("bsd_run_%Y%m%dT%H%M%SZ")

    summary_metrics = {
        "n_target_windows": n_win,
        "target_duration_s": round(n_win * cfg.window_duration_s, 2),
        "n_detection_available_windows": n_det_avail,
        "n_calibration_unavailable_windows": n_cal_unavail,
        "detection_available_duty_cycle": round(n_det_avail / max(1, n_win), 4),
        "n_clean_calibrated_windows": n_clean,
        "clean_duty_cycle": round(n_clean / max(1, n_win), 4),
        "excluded_livetime_s": round((n_win - n_clean) * cfg.window_duration_s, 2),
        "n_active_mode_anomaly_windows": n_anom_active,
        "n_raw_stft_candidates": sum(len(w.raw_stft_candidates) for w in windows),
        "n_broadband_anomaly_windows": sum(1 for w in windows if w.broadband_nominated),
        "n_continuum_drift_warnings": sum(1 for w in windows if "CONTINUUM_DRIFT" in w.warning_causes),
        "n_glitch_storm_warnings": sum(1 for w in windows if "GLITCH_STORM" in w.warning_causes),
        "n_regimes_established": max((w.regime_id for w in windows), default=0) + 1,
        "union_mode_status_counts": union_status_counts,
    }

    run_metadata = {
        "run_id": run_id,
        "timestamp_utc": run_ts.isoformat(),
        "bsd_prototype_version": "0.3.0-phase2b-workstation",
        "git_commit": _get_git_commit(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "config": cfg.to_dict(),
        "target_source": target.source_name,
        "target_sha256": target.sha256,
        "target_preprocessing": target.preprocessing.to_dict(),
        "target_warnings": list(target.warnings),
        "background_source": background.source_name,
        "background_sha256": background.sha256,
        "background_preprocessing": background.preprocessing.to_dict(),
        "background_warnings": list(background.warnings),
        "o3_holdout_accessed": False,
        "prototype_disclaimer": (
            "BSD is an evaluated research prototype (v0.3 Phase 2b). Statistical anomalies and "
            "detection candidates are not independently confirmed physical events."
        ),
    }

    return WorkstationRunResult(
        run_metadata=run_metadata,
        scientific_semantics=_scientific_semantics_block(),
        independence_audit=indep,
        calibration_audit=cal_audit,
        summary_metrics=summary_metrics,
        windows=windows,
        target_observations=tgt_obs,
        candidates_rows=cands_rows,
        timeline_rows=tl_rows,
        ground_truth_audit=None,
    )


def run_synthetic_self_test(
    config: WorkstationConfig | None = None,
    *,
    n_bg_windows: int = 120,
    n_neg_control_windows: int = 40,
    n_pos_injections: int = 10,
    noise_model: str = "pink_colored",
    signal_family: str = "linear_chirp",
    target_snr: float = 18.0,
    seed: int = 2026,
) -> WorkstationRunResult:
    """Reproducible synthetic injection self-test with positive ground truth AND negative control."""
    cfg = config or WorkstationConfig()
    rng = make_split_rng(seed, "workstation_self_test")
    fs = 1024.0
    win_len = int(round(cfg.window_duration_s * fs))

    def _mk_noise(n_w: int) -> np.ndarray:
        return np.concatenate([generate_noise(noise_model, win_len, fs, rng, sigma=1.0)[0] for _ in range(n_w)])

    bg_sig = _mk_noise(n_bg_windows)
    # Target = n_neg_control_windows pure-noise negative control + (2 * n_pos_injections) test windows
    n_pos_span = max(n_pos_injections * 2, n_pos_injections + 4)
    n_tgt = n_neg_control_windows + n_pos_span
    tgt_chunks = [generate_noise(noise_model, win_len, fs, rng, sigma=1.0)[0] for _ in range(n_tgt)]

    inj_records: list[dict[str, Any]] = []
    inj_indices = [n_neg_control_windows + 2 * k + 1 for k in range(n_pos_injections)]
    for k, idx in enumerate(inj_indices):
        clean, t0, t1, f0, f1, _ = generate_clean_signal(signal_family, win_len, fs, rng)
        scaled = _scale_to_energy_snr(clean, 1.0, target_snr)
        tgt_chunks[idx] = tgt_chunks[idx] + scaled
        inj_records.append(
            {
                "injection_id": k,
                "window_index": idx,
                "t_abs_start_s": round(idx * cfg.window_duration_s + t0, 3),
                "t_abs_end_s": round(idx * cfg.window_duration_s + t1, 3),
                "f_min_hz": round(f0, 1),
                "f_max_hz": round(f1, 1),
                "signal_family": signal_family,
                "target_snr": target_snr,
            }
        )
    tgt_sig = np.concatenate(tgt_chunks)

    def _wrap_synth(arr: np.ndarray, name: str, t_int: tuple[float, float]) -> LoadedTimeSeries:
        prep = PreprocessingSummary(
            file_format="synthetic",
            channel_selected="synthetic_1d",
            units="normalized_sigma",
            original_sample_rate_hz=fs,
            effective_sample_rate_hz=fs,
            resampled=False,
            antialias_filter="none",
            effective_bandwidth_hz=(0.0, fs / 2.0),
            analysis_bands_hz=(20.0, 450.0),
            ligo_whitening_applied=False,
            steps=(f"Generated synthetic {noise_model} stream (seed={seed})",),
        )
        return LoadedTimeSeries(
            source_name=name,
            file_format="synthetic",
            sha256=f"synthetic_{name}_{seed}_{arr.size}",
            n_samples=int(arr.size),
            sample_rate_hz=fs,
            duration_s=float(arr.size / fs),
            signal=arr,
            preprocessing=prep,
            time_interval=t_int,
            time_reference="SYNTHETIC_CLOCK",
        )

    bg_ts = _wrap_synth(bg_sig, f"synth_bg_{noise_model}", (0.0, bg_sig.size / fs))
    tgt_ts = _wrap_synth(
        tgt_sig,
        f"synth_target_{noise_model}_{signal_family}",
        (bg_sig.size / fs + 100.0, bg_sig.size / fs + 100.0 + tgt_sig.size / fs),
    )

    res = run_workstation_analysis(tgt_ts, bg_ts, cfg)

    # Evaluate negative control (pure noise windows 0..n_neg_control_windows-1)
    neg_wins = res.windows[:n_neg_control_windows]
    neg_fp = sum(
        1
        for w in neg_wins
        if w.detection_available
        and (w.stft_nominated_post_veto if cfg.search_mode == "transient" else w.broadband_nominated)
    )
    neg_ci = _wilson_ci90(neg_fp, len(neg_wins))

    # Evaluate positive injections (explicitly reporting both detected AND missed!)
    detected_injs: list[dict[str, Any]] = []
    missed_injs: list[dict[str, Any]] = []
    for rec in inj_records:
        w = res.windows[rec["window_index"]]
        hit = bool(
            w.detection_available
            and (w.stft_nominated_post_veto if cfg.search_mode == "transient" else w.broadband_nominated)
        )
        rec_out = {
            **rec,
            "detected": hit,
            "stft_max_score": round(w.stft_max_score, 3),
            "stft_threshold": round(w.stft_threshold, 3) if w.stft_threshold is not None else None,
            "detection_available": w.detection_available,
            "credibility_status": w.credibility_status,
            "warning_causes": ";".join(w.warning_causes),
        }
        if hit:
            detected_injs.append(rec_out)
        else:
            missed_injs.append(rec_out)

    pos_ci = _wilson_ci90(len(detected_injs), len(inj_records))
    gt_audit = {
        "seed": seed,
        "noise_model": noise_model,
        "signal_family": signal_family,
        "target_snr": target_snr,
        "negative_control_windows": len(neg_wins),
        "negative_control_false_alarms": neg_fp,
        "negative_control_far": round(neg_fp / max(1, len(neg_wins)), 4),
        "negative_control_far_ci90": [round(neg_ci[0], 4), round(neg_ci[1], 4)],
        "positive_injections_total": len(inj_records),
        "positive_injections_detected": len(detected_injs),
        "positive_injections_missed": len(missed_injs),
        "detection_efficiency": round(len(detected_injs) / max(1, len(inj_records)), 4),
        "detection_efficiency_ci90": [round(pos_ci[0], 4), round(pos_ci[1], 4)],
        "detected_records": detected_injs,
        "missed_records": missed_injs,
    }

    return WorkstationRunResult(
        run_metadata=res.run_metadata,
        scientific_semantics=res.scientific_semantics,
        independence_audit=res.independence_audit,
        calibration_audit=res.calibration_audit,
        summary_metrics=res.summary_metrics,
        windows=res.windows,
        target_observations=res.target_observations,
        candidates_rows=res.candidates_rows,
        timeline_rows=res.timeline_rows,
        ground_truth_audit=gt_audit,
    )


def export_workstation_run(result: WorkstationRunResult, output_dir: Path | str) -> dict[str, Path]:
    """Export reproducible candidates.csv, timeline.csv, and run_report.json."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    cands_path = out / "candidates.csv"
    tl_path = out / "timeline.csv"
    json_path = out / "run_report.json"

    def _dump_csv(path: Path, rows: list[dict[str, Any]], fallback_fields: list[str]) -> None:
        fields = list(rows[0].keys()) if rows else fallback_fields
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, restval="")
            w.writeheader()
            if rows:
                w.writerows(rows)

    _dump_csv(
        cands_path,
        result.candidates_rows,
        ["candidate_type", "window_index", "t_abs_start_s", "t_abs_end_s", "raw_anomaly_score", "physical_confirmation"],
    )
    _dump_csv(tl_path, result.timeline_rows, ["window_index", "sample_id", "stft_max_score", "credibility_status"])

    window_availability_timeline = [
        {
            "window_index": r["window_index"],
            "sample_id": r["sample_id"],
            "regime_id": r["regime_id"],
            "credibility_status": r["credibility_status"],
            "detection_available": r["detection_available"],
            "active_mode_diagnostic_nomination": r["active_mode_diagnostic_nomination"],
            "active_mode_anomaly_flag": r["active_mode_anomaly_flag"],
            "stft_threshold": r["stft_threshold"],
            "broadband_threshold": r["broadband_threshold"],
            "dedicated_stft_available": r["dedicated_stft_available"],
            "dedicated_broadband_available": r["dedicated_broadband_available"],
            "union_stft_available": r["union_stft_available"],
            "union_broadband_available": r["union_broadband_available"],
            "union_stft_threshold": r["union_stft_threshold"],
            "union_mode_status": r["union_mode_status"],
            "resource_limit_blocked_stft": r["resource_limit_blocked_stft"],
            "resource_limit_blocked_union": r["resource_limit_blocked_union"],
            "calibration_updated_for_next_window": r["calibration_updated_for_next_window"],
        }
        for r in result.timeline_rows
    ]
    report_payload = {
        "run_metadata": result.run_metadata,
        "scientific_semantics": result.scientific_semantics,
        "independence_audit": result.independence_audit.to_dict(),
        "calibration_audit": result.calibration_audit.to_dict(),
        "summary_metrics": result.summary_metrics,
        "window_availability_timeline": window_availability_timeline,
        "ground_truth_audit": result.ground_truth_audit,
        "n_exported_candidates": len(result.candidates_rows),
        "n_exported_timeline_windows": len(result.timeline_rows),
    }
    json_path.write_text(json.dumps(report_payload, indent=2, allow_nan=False), encoding="utf-8")
    return {"candidates_csv": cands_path, "timeline_csv": tl_path, "run_report_json": json_path}
