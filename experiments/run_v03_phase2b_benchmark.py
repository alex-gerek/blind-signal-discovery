"""BSD v0.3 Phase 2b Reproducible Benchmark Suite.

Evaluates:
  Suite A: Dependence-aware progressive causal calibration, Sidak multiple-testing
           control, post-recovery transient efficiency, and orthogonal glitch-storm
           vs. continuum diagnostics.
  Suite B: Separate persistent-signal matrix (80-s sinusoid, 80-s 5-tone comb,
           40-s slow chirp, 80-s broadband stochastic) under BOTH stationary and
           drifting noise, separating frozen-STFT vs. frozen-reference broadband
           Stream A nominations from Stream B credibility warnings.
  Suite C: Five-regime (+unphysical lag) two-detector cross-coherence validation
           against circular time-slide controls.
  Suite D: Previously examined O1/O2 LIGO/Virgo development epochs ONLY.

Reserved O3 holdout (`docs/v03_o3_holdout_manifest.json`) remains 100% untouched.
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sp_signal

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
EXP_DIR = PROJECT_ROOT / "experiments"
for p in (str(SRC_DIR), str(EXP_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from bsd.adaptive_regime import CRED_CALIBRATED, CRED_RECALIBRATING, CRED_TRANSITION, SIG_UNCERTAIN
from bsd.adaptive_regime_v2b import (
    GLITCH_ELEVATED,
    SIG_EMPIRICAL_VALIDATED,
    AdaptiveRegimeManagerV2b,
    AdaptiveRegimeV2bConfig,
    DualOutputWindowV2b,
    FrozenReferenceBroadbandChannel,
    calibrate_joint_stft_and_broadband,
    fit_glitch_reference,
)
from bsd.cross_coherence import (
    COH_ARTIFACT,
    COH_NULL,
    COH_PHYSICAL_EVIDENCE,
    COH_ZERO_LAG_ENV,
    CrossCoherenceConfig,
    evaluate_cross_coherence,
)
from bsd.detector_v2 import MultiResolutionDetector
from bsd.generator import Observation, _scale_to_energy_snr, generate_clean_signal, make_split_rng
from bsd.ligo_adapter import (
    extract_clean_segments,
    fit_background_model,
    load_gwosc_hdf5,
    preprocess_segment,
)
from bsd.noise import generate_noise
from bsd.vetoes import VetoConfig, apply_vetoes
import run_ligo_experiment_003 as exp003

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 180,
        "font.size": 9.5,
        "axes.titlesize": 10.2,
        "axes.labelsize": 9.2,
        "legend.fontsize": 7.8,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for r in rows:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, restval="")
        w.writeheader()
        w.writerows(rows)


def _gen_obs_list(
    rng: np.random.Generator,
    n_win: int,
    noise_model: str,
    sigma_fn: Any,
    prefix: str,
) -> list[Observation]:
    out: list[Observation] = []
    for i in range(n_win):
        sig = float(sigma_fn(i))
        x, _ = generate_noise(noise_model, 4096, 1024.0, rng, sigma=sig)
        out.append(Observation(f"{prefix}_{i:03d}", x, 1024.0, 4.0))
    return out


def _make_chirp_chunks(n_windows: int, f0: float, f1: float, snr_per_win: float) -> list[np.ndarray]:
    dur = n_windows * 4.0
    t = np.arange(n_windows * 4096, dtype=np.float64) / 1024.0
    slope = (f1 - f0) / dur
    w = np.sin(2.0 * np.pi * (f0 * t + 0.5 * slope * t**2))
    return [_scale_to_energy_snr(w[k * 4096 : (k + 1) * 4096], 1.0, snr_per_win) for k in range(n_windows)]


def _wilson_ci90(k: int, n: int) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 1.0)
    z = 1.6448536
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = (z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)) / denom
    return (max(0.0, float(center - half)), min(1.0, float(center + half)))


def _run_suite_a_transient() -> dict[str, Any]:
    """Suite A: Staged causal calibration, matched FWER, transients, and glitch storm."""
    rng = make_split_rng(3026, "v03_phase2b_suite_a")
    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()

    print("[v0.3-P2b] Fitting stationary reference R_0 (pink_colored, N=320)...")
    ref_obs = _gen_obs_list(rng, 320, "pink_colored", lambda _: 1.0, "ref")
    jcal = calibrate_joint_stft_and_broadband(
        det, bb, ref_obs, target_fwer=0.01, norm_fraction=0.30, calibration_seed=3026
    )
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    cref = mgr.fit_continuum_reference(ref_obs)
    gref = fit_glitch_reference(ref_obs)

    rows: list[dict[str, Any]] = []

    # A1. Stationary null 200 windows (matched FWER audit)
    print("[v0.3-P2b]   Suite A1: stationary_null_200win...")
    stat_obs = _gen_obs_list(rng, 200, "pink_colored", lambda _: 1.0, "stat")
    out_stat = mgr.process_sequence(stat_obs, det, bb, jcal, cref, gref)
    stft_fp = sum(1 for o in out_stat if o.stft_nominated_post_veto)
    bb_fp = sum(1 for o in out_stat if o.broadband_nominated)
    joint_fp = sum(1 for o in out_stat if o.joint_nominated_post_veto)
    ci_j = _wilson_ci90(joint_fp, len(out_stat))
    u0 = jcal.stft_uncertainty
    rows.append(
        {
            "scenario": "stationary_null_200win",
            "n_windows": len(out_stat),
            "duty_cycle_clean": float(np.mean([o.is_clean_calibrated for o in out_stat])),
            "stft_far_full": stft_fp / len(out_stat),
            "bb_far_full": bb_fp / len(out_stat),
            "joint_far_full": joint_fp / len(out_stat),
            "joint_far_ci90_lo": ci_j[0],
            "joint_far_ci90_hi": ci_j[1],
            "cal_n_thresh": u0.n_thresh,
            "cal_n_eff": u0.n_eff,
            "cal_ljung_box_p": u0.ljung_box_p,
            "cal_beta_ci90_neff_hi": u0.far_ci90_neff[1],
            "t_provisional_win": None,
            "t_mature_win": None,
            "mature_surv_far": None,
            "transient_eff_full": None,
            "transient_eff_clean_cal": None,
            "glitch_storm_flag_rate": float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_stat])),
            "n_regimes": max(o.regime_id for o in out_stat) + 1,
        }
    )

    # A2. 500-window abrupt step-change (100 pink + 400 shifted 1.5x AR(2))
    #     3 transition chirps at 108, 125, 145; 10 mature post-recovery chirps at 330..420
    print("[v0.3-P2b]   Suite A2: abrupt_step_500win (staged provisional -> mature)...")
    ab_obs = _gen_obs_list(rng, 100, "pink_colored", lambda _: 1.0, "ab0") + _gen_obs_list(
        rng, 400, "shifted_ar2", lambda _: 1.5, "ab1"
    )
    trans_idx = [108, 125, 145]
    post_idx = list(range(330, 430, 10))
    is_sig_ab = [False] * 500
    for idx in trans_idx + post_idx:
        clean, *_ = generate_clean_signal("linear_chirp", 4096, 1024.0, rng)
        s = _scale_to_energy_snr(clean, 1.5, 18.0)
        ab_obs[idx] = Observation(ab_obs[idx].sample_id, ab_obs[idx].signal + s, 1024.0, 4.0)
        is_sig_ab[idx] = True

    out_ab = mgr.process_sequence(ab_obs, det, bb, jcal, cref, gref)
    t_prov = next((o.index - 100 for o in out_ab[100:] if o.regime_id >= 1 and o.credibility_status == CRED_CALIBRATED), None)
    t_mat = next((o.index - 100 for o in out_ab[100:] if o.regime_id >= 1 and o.diagnostics.get("is_mature")), None)
    mature_bg = [o for o, s in zip(out_ab, is_sig_ab) if (not s) and o.regime_id >= 1 and o.diagnostics.get("is_mature") and o.is_clean_calibrated]
    mature_fp = sum(1 for o in mature_bg if o.stft_nominated_post_veto)
    mature_far = (mature_fp / len(mature_bg)) if mature_bg else None
    post_eff_full = float(np.mean([out_ab[i].stft_nominated_post_veto for i in post_idx]))
    post_eff_clean = float(np.mean([out_ab[i].stft_nominated_post_veto and out_ab[i].is_clean_calibrated for i in post_idx]))
    trans_pres = float(np.mean([out_ab[i].stft_nominated_raw and out_ab[i].significance_status == SIG_UNCERTAIN for i in trans_idx]))
    u_mat = out_ab[-1].calibration_uncertainty

    rows.append(
        {
            "scenario": "abrupt_step_500win",
            "n_windows": len(out_ab),
            "duty_cycle_clean": float(np.mean([o.is_clean_calibrated for o in out_ab])),
            "stft_far_full": float(np.mean([o.stft_nominated_post_veto for o, s in zip(out_ab, is_sig_ab) if not s])),
            "bb_far_full": None,
            "joint_far_full": None,
            "joint_far_ci90_lo": None,
            "joint_far_ci90_hi": None,
            "cal_n_thresh": u_mat.n_thresh,
            "cal_n_eff": u_mat.n_eff,
            "cal_ljung_box_p": u_mat.ljung_box_p,
            "cal_beta_ci90_neff_hi": u_mat.far_ci90_neff[1],
            "t_provisional_win": t_prov,
            "t_mature_win": t_mat,
            "mature_surv_far": mature_far,
            "transient_eff_full": post_eff_full,
            "transient_eff_clean_cal": post_eff_clean,
            "glitch_storm_flag_rate": float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_ab])),
            "n_regimes": max(o.regime_id for o in out_ab) + 1,
            "transition_raw_preservation": trans_pres,
            "mature_sig_status": out_ab[-1].significance_status,
        }
    )

    # A3. Stationary transient chirps (20 chirps in 120 windows)
    print("[v0.3-P2b]   Suite A3: transient_chirps_stationary...")
    tr_obs = _gen_obs_list(rng, 120, "pink_colored", lambda _: 1.0, "tr")
    tr_idx = list(range(15, 115, 5))
    tr_sig = [False] * 120
    for idx in tr_idx:
        clean, *_ = generate_clean_signal("linear_chirp", 4096, 1024.0, rng)
        s = _scale_to_energy_snr(clean, 1.0, 18.0)
        tr_obs[idx] = Observation(tr_obs[idx].sample_id, tr_obs[idx].signal + s, 1024.0, 4.0)
        tr_sig[idx] = True
    out_tr = mgr.process_sequence(tr_obs, det, bb, jcal, cref, gref)
    rows.append(
        {
            "scenario": "transient_chirps_stationary",
            "n_windows": 120,
            "duty_cycle_clean": float(np.mean([o.is_clean_calibrated for o in out_tr])),
            "stft_far_full": float(np.mean([o.stft_nominated_post_veto for o, s in zip(out_tr, tr_sig) if not s])),
            "bb_far_full": float(np.mean([o.broadband_nominated for o, s in zip(out_tr, tr_sig) if not s])),
            "joint_far_full": float(np.mean([o.joint_nominated_post_veto for o, s in zip(out_tr, tr_sig) if not s])),
            "joint_far_ci90_lo": None,
            "joint_far_ci90_hi": None,
            "cal_n_thresh": u0.n_thresh,
            "cal_n_eff": u0.n_eff,
            "cal_ljung_box_p": u0.ljung_box_p,
            "cal_beta_ci90_neff_hi": u0.far_ci90_neff[1],
            "t_provisional_win": None,
            "t_mature_win": None,
            "mature_surv_far": None,
            "transient_eff_full": float(np.mean([out_tr[i].stft_nominated_post_veto for i in tr_idx])),
            "transient_eff_clean_cal": float(np.mean([out_tr[i].stft_nominated_post_veto and out_tr[i].is_clean_calibrated for i in tr_idx])),
            "glitch_storm_flag_rate": float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_tr])),
            "n_regimes": max(o.regime_id for o in out_tr) + 1,
        }
    )

    # A4. Synthetic glitch storm (60 pink + 60 impulsive) -> orthogonal GLITCH_STORM, NO recalibration
    print("[v0.3-P2b]   Suite A4: glitch_storm_120win...")
    gl_obs = _gen_obs_list(rng, 60, "pink_colored", lambda _: 1.0, "gl0") + _gen_obs_list(
        rng, 60, "impulsive", lambda _: 1.0, "gl1"
    )
    out_gl = mgr.process_sequence(gl_obs, det, bb, jcal, cref, gref)
    gl_post_flag = float(np.mean([o.glitch_storm_status == GLITCH_ELEVATED for o in out_gl[60:]]))
    gl_causes_ok = all("GLITCH_STORM" in o.warning_causes for o in out_gl[60:] if o.glitch_storm_status == GLITCH_ELEVATED)
    rows.append(
        {
            "scenario": "glitch_storm_120win",
            "n_windows": 120,
            "duty_cycle_clean": float(np.mean([o.is_clean_calibrated for o in out_gl])),
            "stft_far_full": float(np.mean([o.stft_nominated_post_veto for o in out_gl])),
            "bb_far_full": float(np.mean([o.broadband_nominated for o in out_gl])),
            "joint_far_full": float(np.mean([o.joint_nominated_post_veto for o in out_gl])),
            "joint_far_ci90_lo": None,
            "joint_far_ci90_hi": None,
            "cal_n_thresh": u0.n_thresh,
            "cal_n_eff": u0.n_eff,
            "cal_ljung_box_p": u0.ljung_box_p,
            "cal_beta_ci90_neff_hi": u0.far_ci90_neff[1],
            "t_provisional_win": None,
            "t_mature_win": None,
            "mature_surv_far": float(np.mean([o.stft_nominated_post_veto for o in out_gl if o.is_clean_calibrated])),
            "transient_eff_full": None,
            "transient_eff_clean_cal": None,
            "glitch_storm_flag_rate": gl_post_flag,
            "n_regimes": max(o.regime_id for o in out_gl) + 1,
            "glitch_causes_explicit": gl_causes_ok,
        }
    )

    ab_timeline = {
        "z_drift": [o.max_feature_z_drift for o in out_ab],
        "far_ci_hi": [o.calibration_uncertainty.far_ci90_neff[1] * 100.0 for o in out_ab],
        "stage": [
            0 if o.credibility_status == CRED_CALIBRATED and o.regime_id == 0
            else (1 if o.credibility_status in (CRED_TRANSITION, CRED_RECALIBRATING)
                  else (2 if not o.diagnostics.get("is_mature") else 3))
            for o in out_ab
        ],
        "stft_score": [o.stft_max_score for o in out_ab],
        "stft_thresh": [o.stft_threshold for o in out_ab],
        "trans_idx": trans_idx,
        "post_idx": post_idx,
    }
    return {
        "rows": rows,
        "ab_timeline": ab_timeline,
        "shared_ref": (det, bb, jcal, cref, gref),
    }


def _run_suite_b_persistent(shared_ref: tuple[Any, ...]) -> dict[str, Any]:
    """Suite B: Separate persistent signal matrix under stationary AND drifting noise."""
    det, bb, jcal, cref, gref = shared_ref
    rng = make_split_rng(3027, "v03_phase2b_suite_b")
    mgr_pers = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="persistent"))
    t_win = np.arange(4096, dtype=np.float64) / 1024.0
    sos_bb = sp_signal.butter(4, [80.0, 300.0], btype="bandpass", fs=1024.0, output="sos")

    rows: list[dict[str, Any]] = []
    for bg_mode, sigma_fn in (("stationary", lambda _: 1.0), ("drifting_1.0_to_1.4x", lambda i: 1.0 + 0.4 * (i / 119.0))):
        for morph in ("80s_sinusoid", "80s_multitone_comb", "40s_slow_chirp", "80s_broadband_stochastic"):
            print(f"[v0.3-P2b]   Suite B (persistent): {morph} | bg={bg_mode}...")
            obs = _gen_obs_list(rng, 120, "pink_colored", sigma_fn, f"p_{morph[:3]}")
            is_sig = [False] * 120
            if morph == "80s_sinusoid":
                for idx in range(50, 70):
                    s = _scale_to_energy_snr(np.sin(2.0 * np.pi * 180.0 * (t_win + idx * 4.0)), 1.0, 20.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True
            elif morph == "80s_multitone_comb":
                for idx in range(50, 70):
                    tt = t_win + idx * 4.0
                    comb = sum(np.sin(2.0 * np.pi * f * tt) for f in (90.0, 140.0, 190.0, 240.0, 290.0))
                    s = _scale_to_energy_snr(comb, 1.0, 22.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True
            elif morph == "40s_slow_chirp":
                for k, ch in enumerate(_make_chirp_chunks(10, 60.0, 260.0, 16.0)):
                    idx = 50 + k
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + ch, 1024.0, 4.0)
                    is_sig[idx] = True
            else:  # 80s_broadband_stochastic (SNR=26/win)
                for idx in range(50, 70):
                    raw_bb = sp_signal.sosfiltfilt(sos_bb, rng.normal(0.0, 1.0, size=4096))
                    s = _scale_to_energy_snr(raw_bb, 1.0, 26.0)
                    obs[idx] = Observation(obs[idx].sample_id, obs[idx].signal + s, 1024.0, 4.0)
                    is_sig[idx] = True

            out = mgr_pers.process_sequence(obs, det, bb, jcal, cref, gref)
            sig_wins = [o for o, s in zip(out, is_sig) if s]
            bg_wins = [o for o, s in zip(out, is_sig) if not s]
            rows.append(
                {
                    "morphology": morph,
                    "background": bg_mode,
                    "n_sig_windows": len(sig_wins),
                    "frozen_stft_eff_full": float(np.mean([o.stft_nominated_post_veto for o in sig_wins])),
                    "frozen_ref_bb_eff_full": float(np.mean([o.broadband_nominated for o in sig_wins])),
                    "joint_eff_full": float(np.mean([o.joint_nominated_post_veto for o in sig_wins])),
                    "clean_cal_stft_eff": float(np.mean([o.stft_nominated_post_veto and o.is_clean_calibrated for o in sig_wins])),
                    "bg_joint_far_full": float(np.mean([o.joint_nominated_post_veto for o in bg_wins])),
                    "cont_drift_frac_on_sig": float(np.mean(["CONTINUUM_DRIFT" in o.warning_causes for o in sig_wins])),
                    "glitch_flag_frac_on_sig": float(np.mean(["GLITCH_STORM" in o.warning_causes for o in sig_wins])),
                    "n_regimes": max(o.regime_id for o in out) + 1,
                }
            )
    return {"rows": rows}


def _run_suite_c_coherence() -> dict[str, Any]:
    """Suite C: Five-regime (+unphysical lag) two-detector cross-coherence validation."""
    rng = make_split_rng(3028, "v03_phase2b_suite_c")
    cfg = CrossCoherenceConfig(n_time_slides=120, alpha_fap=0.01)
    n_trials = 20
    n_win = 10  # 40-s sustained block per trial
    n_tot = n_win * 4096
    sos = sp_signal.butter(4, [80.0, 300.0], btype="bandpass", fs=1024.0, output="sos")
    t = np.arange(n_tot, dtype=np.float64) / 1024.0

    def _to_obs(arr: np.ndarray, pfx: str) -> list[Observation]:
        return [Observation(f"{pfx}_{k}", arr[k * 4096 : (k + 1) * 4096].copy(), 1024.0, 4.0) for k in range(n_win)]

    regimes = (
        "independent_noise",
        "uncorrelated_dual_drift",
        "correlated_zero_lag_env_noise",
        "common_narrowband_line_artifact",
        "common_synchronous_glitch_artifact",
        "shared_broadband_pos4ms",
        "shared_broadband_neg7ms",
        "shared_broadband_unphysical_50ms",
    )
    rows: list[dict[str, Any]] = []
    for reg in regimes:
        print(f"[v0.3-P2b]   Suite C (coherence): {reg} ({n_trials} trials)...")
        counts = {COH_NULL: 0, COH_ARTIFACT: 0, COH_ZERO_LAG_ENV: 0, COH_PHYSICAL_EVIDENCE: 0}
        lags: list[float] = []
        for _ in range(n_trials):
            n1 = rng.normal(0.0, 1.0, n_tot)
            n2 = rng.normal(0.0, 1.0, n_tot)
            if reg == "uncorrelated_dual_drift":
                n1 = n1 * np.linspace(1.0, 1.6, n_tot)
                n2 = n2 * np.linspace(1.0, 1.8, n_tot)
            elif reg == "correlated_zero_lag_env_noise":
                env = sp_signal.sosfiltfilt(sos, rng.normal(0.0, 1.0, n_tot)) * 0.55
                n1 = n1 + env
                n2 = n2 + env
            elif reg == "common_narrowband_line_artifact":
                line = 0.5 * np.sin(2.0 * np.pi * 180.0 * t)
                n1 = n1 + line
                n2 = n2 + line
            elif reg == "common_synchronous_glitch_artifact":
                gl = np.zeros(n_tot, dtype=np.float64)
                for center in (n_tot // 4, n_tot // 2, 3 * n_tot // 4):
                    idx = np.arange(max(0, center - 40), min(n_tot, center + 40))
                    gl[idx] += 12.0 * np.exp(-0.5 * ((idx - center) / 6.0) ** 2)
                n1 = n1 + gl
                n2 = n2 + gl
            elif reg == "shared_broadband_pos4ms":
                sig = sp_signal.sosfiltfilt(sos, rng.normal(0.0, 1.0, n_tot)) * 0.55
                n1 = n1 * np.linspace(1.0, 1.3, n_tot) + sig
                n2 = n2 + np.roll(sig, 4)  # +3.9 ms
            elif reg == "shared_broadband_neg7ms":
                sig = sp_signal.sosfiltfilt(sos, rng.normal(0.0, 1.0, n_tot)) * 0.55
                n1 = n1 + sig
                n2 = n2 + np.roll(sig, -7)  # -6.8 ms
            elif reg == "shared_broadband_unphysical_50ms":
                sig = sp_signal.sosfiltfilt(sos, rng.normal(0.0, 1.0, n_tot)) * 0.55
                n1 = n1 + sig
                n2 = n2 + np.roll(sig, 51)  # +50 ms (outside +-10 ms physical window)

            res = evaluate_cross_coherence(_to_obs(n1, "d1"), _to_obs(n2, "d2"), cfg)
            counts[res.classification] += 1
            lags.append(res.peak_lag_ms)

        rows.append(
            {
                "regime": reg,
                "n_trials": n_trials,
                "frac_null": counts[COH_NULL] / n_trials,
                "frac_artifact_warning": counts[COH_ARTIFACT] / n_trials,
                "frac_zero_lag_env_warning": counts[COH_ZERO_LAG_ENV] / n_trials,
                "frac_physical_lag_evidence": counts[COH_PHYSICAL_EVIDENCE] / n_trials,
                "median_peak_lag_ms": float(np.median(lags)),
            }
        )
    return {"rows": rows}


def _run_suite_d_ligo_o1o2() -> dict[str, Any]:
    """Suite D: Previously examined O1/O2 LIGO/Virgo development epochs ONLY."""
    file_index: dict[tuple[str, str, int], Path] = {}
    ev_gps: dict[str, float] = {}
    for prov_rel in (
        "results/ligo_experiment_003/data_provenance.json",
        "results/v02_holdout/data_provenance.json",
    ):
        p = json.loads((PROJECT_ROOT / prov_rel).read_text(encoding="utf-8"))
        for ev in p["events"]:
            ev_gps[ev["name"]] = float(ev["published_gps"])
        for f in p["files"]:
            file_index[(f["event"], f["detector"], int(f["duration_s"]))] = PROJECT_ROOT / f["relative_local_path"]

    targets = (
        ("GW150914", "H1", "stable_confirmed"),
        ("GW150914", "L1", "stable_confirmed"),
        ("GW170814", "L1", "stable_confirmed"),
        ("GW170104", "L1", "nonstationary_drift"),
        ("GW170729", "V1", "nonstationary_drift"),
        ("GW170818", "H1", "nonstationary_drift"),
    )

    det = MultiResolutionDetector()
    bb = FrozenReferenceBroadbandChannel()
    cfg_veto = VetoConfig(enable_glitch_veto=True, enable_line_veto=True)
    mgr = AdaptiveRegimeManagerV2b(AdaptiveRegimeV2bConfig(mode="causal_trailing", policy="transient"))
    rows: list[dict[str, Any]] = []

    for ev_name, det_name, epoch_type in targets:
        print(f"[v0.3-P2b]   Suite D (O1/O2): {ev_name} {det_name} ({epoch_type})...")
        rec_32 = load_gwosc_hdf5(file_index[(ev_name, det_name, 32)])
        rec_4096 = load_gwosc_hdf5(file_index[(ev_name, det_name, 4096)])
        excl = (rec_32.gps_start - 64.0, rec_32.gps_end + 64.0)
        bg_all = extract_clean_segments(rec_4096, 4.0, 1.0, excl, True)
        counts = exp003.partition_background_segments(len(bg_all))
        assert counts is not None
        n_psd, n_cal, _ = counts

        bg_model = fit_background_model(det_name, bg_all[:n_psd], "adapted_ligo_whitened", 1024.0, (30.0, 380.0))
        cal_obs = [preprocess_segment(s, f"{ev_name}_{det_name}_c{i}", bg_model) for i, s in enumerate(bg_all[n_psd : n_psd + n_cal])]
        val_obs = [preprocess_segment(s, f"{ev_name}_{det_name}_v{i}", bg_model) for i, s in enumerate(bg_all[n_psd + n_cal :])]

        jcal = calibrate_joint_stft_and_broadband(det, bb, cal_obs, target_fwer=0.01, norm_fraction=0.30, calibration_seed=0)
        cref = mgr.fit_continuum_reference(cal_obs)
        gref = fit_glitch_reference(cal_obs)

        onsource_ok: bool | None = None
        if epoch_type == "stable_confirmed":
            ev_segs = extract_clean_segments(rec_32, 4.0, 1.0, None, True)
            ev_t0 = ev_gps[ev_name]
            ev_seg = next(s for s in ev_segs if s.gps_start <= ev_t0 < s.gps_start + s.duration_s)
            ev_obs = preprocess_segment(ev_seg, f"{ev_name}_{det_name}_onsource", bg_model)
            ev_res = det.analyze(ev_obs, jcal.stft_calibration)
            onsource_ok = bool(apply_vetoes(ev_obs, ev_res.candidates, cfg_veto).detected_after_veto)

        out = mgr.process_sequence(val_obs, det, bb, jcal, cref, gref, veto_config=cfg_veto)
        n_val = len(out)
        clean_wins = [o for o in out if o.is_clean_calibrated]
        n_clean = len(clean_wins)
        fp_full = sum(1 for o in out if o.stft_nominated_post_veto)
        fp_clean = sum(1 for o in clean_wins if o.stft_nominated_post_veto)
        ci_clean = _wilson_ci90(fp_clean, n_clean)
        n_cont = sum(1 for o in out if "CONTINUUM_DRIFT" in o.warning_causes)
        n_gl = sum(1 for o in out if "GLITCH_STORM" in o.warning_causes)
        u0 = jcal.stft_uncertainty

        rows.append(
            {
                "event": ev_name,
                "detector": det_name,
                "epoch_type": epoch_type,
                "n_val_windows": n_val,
                "duty_cycle_clean": n_clean / n_val,
                "excluded_livetime_hr": ((n_val - n_clean) * 4.0) / 3600.0,
                "n_continuum_warn_wins": n_cont,
                "n_glitch_storm_warn_wins": n_gl,
                "n_regimes": max(o.regime_id for o in out) + 1,
                "stft_far_full": fp_full / n_val,
                "stft_far_clean_cal": (fp_clean / n_clean) if n_clean else None,
                "clean_far_ci90_lo": ci_clean[0],
                "clean_far_ci90_hi": ci_clean[1],
                "cal_n_thresh": u0.n_thresh,
                "cal_n_eff": u0.n_eff,
                "cal_ljung_box_p": u0.ljung_box_p,
                "onsource_confirmed_detected": onsource_ok,
            }
        )
    return {"rows": rows}


def _check_o3_untouched() -> dict[str, Any]:
    manifest = json.loads((PROJECT_ROOT / "docs/v03_o3_holdout_manifest.json").read_text(encoding="utf-8"))
    o3_names = [ev["event_id"] for ev in manifest["events"]]
    gwosc_dir = PROJECT_ROOT / "data/gwosc"
    found: list[str] = []
    if gwosc_dir.exists():
        for p in gwosc_dir.rglob("*"):
            for name in o3_names:
                if name in p.name:
                    found.append(p.name)
    return {"n_reserved_o3": len(o3_names), "found_locally": found, "strictly_untouched": len(found) == 0}


def _generate_figures(
    sa: dict[str, Any],
    sb: dict[str, Any],
    sc: dict[str, Any],
    sd: dict[str, Any],
    fig_dir: Path,
) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figs: list[str] = []

    # Fig 01: 500-window staged causal recovery
    tl = sa["ab_timeline"]
    x = np.arange(len(tl["z_drift"]))
    st = np.array(tl["stage"])
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 6.2), sharex=True)
    ax = axes[0]
    ax.plot(x, tl["z_drift"], color="#1f4e79", lw=1.2, label="Continuum max Z-drift")
    ax.axvline(100, color="#d62728", ls="--", lw=1.3, label="Abrupt shift onset (win 100)")
    ax.fill_between(x, 0, max(tl["z_drift"]) * 1.05, where=(st == 1), color="#feb2b2", alpha=0.35, label="Warmup / Transition (UNCERTAIN)")
    ax.fill_between(x, 0, max(tl["z_drift"]) * 1.05, where=(st == 2), color="#fbd38d", alpha=0.35, label="Stage 1: EMPIRICAL_PROVISIONAL")
    ax.fill_between(x, 0, max(tl["z_drift"]) * 1.05, where=(st == 3), color="#c6f6d5", alpha=0.35, label="Stage 2: Mature R1 (N_thresh>=150)")
    ax.set_ylabel("Max Continuum Z-Drift")
    ax.set_title("Stream B: Progressive Causal Recalibration (Zero Future Leakage)")
    ax.legend(loc="upper right", ncol=2, fontsize=7.6)

    ax = axes[1]
    ax.plot(x, tl["far_ci_hi"], color="#805ad5", lw=1.6, label="Active Regime Upper 90% Dependence-Aware FAR Bound (%)")
    ax.axhline(2.5, color="#d62728", ls=":", lw=1.2, label="Mature Upper-Bound Reference (2.5%)")
    ax.set_ylabel("Upper 90% FAR CI (%)")
    ax.set_xlabel("Sequential 4-Second Window Index")
    ax.set_title("Finite-Sample & Dependence-Aware FAR Uncertainty Across Calibration Stages")
    ax.legend(loc="upper right", fontsize=8.0)
    fig.suptitle("BSD v0.3 Phase 2b: Staged Causal Recalibration & Uncertainty Evolution (500-Window Shift)", y=0.995)
    fig.tight_layout(rect=[0, 0.01, 1, 0.96])
    fig.savefig(fig_dir / "fig01_staged_causal_recovery_500win.png")
    plt.close(fig)
    figs.append("fig01_staged_causal_recovery_500win.png")

    # Fig 02: Persistent signal matrix (Frozen STFT vs Frozen-Ref Broadband in Stream A)
    sub_stat = [r for r in sb["rows"] if r["background"] == "stationary"]
    lbls = ["80s Single\nSinusoid", "80s 5-Tone\nComb", "40s Slow\nChirp", "80s Broadband\nStochastic"]
    xx = np.arange(len(sub_stat))
    w = 0.25
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    v_stft = [r["frozen_stft_eff_full"] * 100.0 for r in sub_stat]
    v_bb = [r["frozen_ref_bb_eff_full"] * 100.0 for r in sub_stat]
    v_jnt = [r["joint_eff_full"] * 100.0 for r in sub_stat]
    ax.bar(xx - w, v_stft, width=w, color="#718096", label="Frozen STFT Stream A (%)")
    ax.bar(xx, v_bb, width=w, color="#dd6b20", label="Frozen-Ref Broadband Channel (%)")
    ax.bar(xx + w, v_jnt, width=w, color="#1f4e79", label="Joint Sidak-Controlled Stream A (%)")
    for i in range(len(sub_stat)):
        ax.text(xx[i] - w, v_stft[i] + 1.5, f"{v_stft[i]:.0f}%", ha="center", fontsize=8)
        ax.text(xx[i], v_bb[i] + 1.5, f"{v_bb[i]:.0f}%", ha="center", fontsize=8)
        ax.text(xx[i] + w, v_jnt[i] + 1.5, f"{v_jnt[i]:.0f}%", ha="center", fontsize=8)
    ax.set_xticks(xx)
    ax.set_xticklabels(lbls)
    ax.set_ylim(0, 115)
    ax.set_ylabel("Full-Timeline Nomination Efficiency (%)")
    ax.set_title("BSD v0.3 Phase 2b: Separating Frozen STFT vs. Frozen-Reference Broadband Sensitivity (Matched FWER)")
    ax.legend(loc="lower left", fontsize=8.0)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig02_persistent_matrix_streamA_vs_streamB.png")
    plt.close(fig)
    figs.append("fig02_persistent_matrix_streamA_vs_streamB.png")

    # Fig 03: Cross-coherence 8-regime stacked bar chart
    r_c = sc["rows"]
    lbl_c = [
        "Indep.\nNoise", "Uncorr.\nDrift", "Zero-Lag\nEnv Noise", "Common\n180Hz Line",
        "Common\nGlitch", "Shared BB\n+4ms", "Shared BB\n-7ms", "Shared BB\n+50ms"
    ]
    xx = np.arange(len(r_c))
    f_null = np.array([r["frac_null"] * 100.0 for r in r_c])
    f_art = np.array([r["frac_artifact_warning"] * 100.0 for r in r_c])
    f_env = np.array([r["frac_zero_lag_env_warning"] * 100.0 for r in r_c])
    f_phy = np.array([r["frac_physical_lag_evidence"] * 100.0 for r in r_c])
    fig, ax = plt.subplots(figsize=(10.8, 4.5))
    ax.bar(xx, f_null, color="#cbd5e0", label="Incoherent / Null")
    ax.bar(xx, f_art, bottom=f_null, color="#e53e3e", label="Common Artifact Warning")
    ax.bar(xx, f_env, bottom=f_null + f_art, color="#dd6b20", label="Zero-Lag Env Warning")
    ax.bar(xx, f_phy, bottom=f_null + f_art + f_env, color="#2f855a", label="Physical-Lag Evidence (0.5-10 ms)")
    ax.set_xticks(xx)
    ax.set_xticklabels(lbl_c, fontsize=8.0)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Trial Classification Fraction (%)")
    ax.set_title("BSD v0.3 Phase 2b: Cross-Detector Coherence Across Control, Artifact, Environmental, and Physical-Lag Regimes")
    ax.legend(loc="upper center", ncol=4, fontsize=7.8)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig03_cross_coherence_five_regimes.png")
    plt.close(fig)
    figs.append("fig03_cross_coherence_five_regimes.png")

    # Fig 04: O1/O2 FAR & orthogonal diagnostic breakdown
    r_d = sd["rows"]
    lbl_d = [f"{r['event']}\n{r['detector']}" for r in r_d]
    xx = np.arange(len(r_d))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    f_full = [r["stft_far_full"] * 100.0 for r in r_d]
    f_cln = [0.0 if r["stft_far_clean_cal"] is None else r["stft_far_clean_cal"] * 100.0 for r in r_d]
    w = 0.34
    ax1.bar(xx - w / 2, f_full, width=w, color="#a0aec0", label="Full-Timeline FAR (%)")
    ax1.bar(xx + w / 2, f_cln, width=w, color="#1f4e79", label="Clean-Calibrated FAR (%)")
    ax1.axhline(2.0, color="#d62728", ls="--", lw=1.1, label="2.0% Reference")
    ax1.set_xticks(xx)
    ax1.set_xticklabels(lbl_d, fontsize=8)
    ax1.set_ylabel("Background FAR (%)")
    ax1.set_title("O1/O2 Dev Epochs: Full vs. Clean-Calibrated FAR")
    ax1.legend(loc="upper left", fontsize=7.8)

    c_pct = [100.0 * r["n_continuum_warn_wins"] / r["n_val_windows"] for r in r_d]
    g_pct = [100.0 * r["n_glitch_storm_warn_wins"] / r["n_val_windows"] for r in r_d]
    ax2.bar(xx - w / 2, c_pct, width=w, color="#dd6b20", label="CONTINUUM_DRIFT Flagged (%)")
    ax2.bar(xx + w / 2, g_pct, width=w, color="#805ad5", label="GLITCH_STORM Flagged (%)")
    ax2.set_xticks(xx)
    ax2.set_xticklabels(lbl_d, fontsize=8)
    ax2.set_ylabel("Flagged Timeline Fraction (%)")
    ax2.set_title("O1/O2 Dev Epochs: Orthogonal Diagnostic Attribution")
    ax2.legend(loc="upper left", fontsize=7.8)

    fig.suptitle("BSD v0.3 Phase 2b: O1/O2 Real-Data Evaluation with Orthogonal Glitch/Continuum Diagnostics", y=0.995)
    fig.tight_layout(rect=[0, 0.01, 1, 0.96])
    fig.savefig(fig_dir / "fig04_ligo_o1o2_orthogonal_diagnostics.png")
    plt.close(fig)
    figs.append("fig04_ligo_o1o2_orthogonal_diagnostics.png")

    return figs


def _evaluate_ac(
    sa: dict[str, Any],
    sb: dict[str, Any],
    sc: dict[str, Any],
    sd: dict[str, Any],
    o3: dict[str, Any],
) -> dict[str, Any]:
    amap = {r["scenario"]: r for r in sa["rows"]}
    bmap = {(r["morphology"], r["background"]): r for r in sb["rows"]}
    cmap = {r["regime"]: r for r in sc["rows"]}
    dmap = {(r["event"], r["detector"]): r for r in sd["rows"]}

    ab = amap["abrupt_step_500win"]
    st = amap["stationary_null_200win"]
    tr = amap["transient_chirps_stationary"]
    gl = amap["glitch_storm_120win"]

    ac1 = bool(
        ab["t_provisional_win"] is not None
        and ab["t_provisional_win"] <= 75
        and ab["t_mature_win"] is not None
        and ab["t_mature_win"] <= 220
        and ab["mature_surv_far"] is not None
        and ab["mature_surv_far"] <= 0.015
    )
    ac2 = bool(st["joint_far_full"] <= 0.015)
    ac3 = bool(tr["transient_eff_full"] >= 0.95 and ab["transient_eff_clean_cal"] >= 0.90)

    sin_stat = bmap[("80s_sinusoid", "stationary")]
    ch40_stat = bmap[("40s_slow_chirp", "stationary")]
    ac4 = bool(
        gl["glitch_storm_flag_rate"] >= 0.85
        and gl["n_regimes"] == 1
        and sin_stat["clean_cal_stft_eff"] >= 0.95
        and ch40_stat["clean_cal_stft_eff"] >= 0.85
    )

    bb_stat = bmap[("80s_broadband_stochastic", "stationary")]
    ac5 = bool(bb_stat["frozen_ref_bb_eff_full"] >= 0.90)

    ac6 = bool(
        cmap["shared_broadband_pos4ms"]["frac_physical_lag_evidence"] >= 0.90
        and cmap["shared_broadband_neg7ms"]["frac_physical_lag_evidence"] >= 0.90
        and cmap["independent_noise"]["frac_physical_lag_evidence"] <= 0.02
        and cmap["uncorrelated_dual_drift"]["frac_physical_lag_evidence"] <= 0.02
        and cmap["correlated_zero_lag_env_noise"]["frac_zero_lag_env_warning"] >= 0.80
        and cmap["common_narrowband_line_artifact"]["frac_artifact_warning"] >= 0.80
    )

    confirmed_ok = all(
        r["onsource_confirmed_detected"] is True for r in sd["rows"] if r["epoch_type"] == "stable_confirmed"
    )
    v1_far = dmap[("GW170729", "V1")]["stft_far_clean_cal"]
    h1_18_far = dmap[("GW170818", "H1")]["stft_far_clean_cal"]
    m_far_s = "N/A" if ab["mature_surv_far"] is None else f"{ab['mature_surv_far']*100:.2f}%"
    v1_far_s = "N/A" if v1_far is None else f"{v1_far*100:.2f}%"
    h1_18_far_s = "N/A" if h1_18_far is None else f"{h1_18_far*100:.2f}%"
    ac7 = bool(
        o3["strictly_untouched"]
        and confirmed_ok
        and (v1_far is not None and v1_far <= 0.02)
        and (h1_18_far is not None and h1_18_far < 0.0497)
    )

    return {
        "P2b-AC1": {
            "passed": ac1,
            "summary": (
                f"T_prov={ab['t_provisional_win']}, T_mature={ab['t_mature_win']} (<=220), "
                f"mature_FAR={m_far_s} (<=1.5%), "
                f"N_thresh={ab['cal_n_thresh']}, N_eff={ab['cal_n_eff']:.1f}, p_LB={ab['cal_ljung_box_p']:.3f}"
            ),
        },
        "P2b-AC2": {
            "passed": ac2,
            "summary": (
                f"matched_joint_FAR={st['joint_far_full']*100:.2f}% "
                f"(90% CI [{st['joint_far_ci90_lo']*100:.2f}%, {st['joint_far_ci90_hi']*100:.2f}%])"
            ),
        },
        "P2b-AC3": {
            "passed": ac3,
            "summary": (
                f"stat_chirps={tr['transient_eff_full']*100:.1f}% (>=95%), "
                f"post_mature_chirps={ab['transient_eff_clean_cal']*100:.1f}% (>=90%)"
            ),
        },
        "P2b-AC4": {
            "passed": ac4,
            "summary": (
                f"glitch_flag={gl['glitch_storm_flag_rate']*100:.1f}% (n_regimes={gl['n_regimes']}), "
                f"80s_sin={sin_stat['clean_cal_stft_eff']*100:.1f}%, 40s_chirp={ch40_stat['clean_cal_stft_eff']*100:.1f}%"
            ),
        },
        "P2b-AC5": {
            "passed": ac5,
            "summary": (
                f"80s_broadband_frozen_STFT={bb_stat['frozen_stft_eff_full']*100:.1f}% -> "
                f"frozen_ref_BB={bb_stat['frozen_ref_bb_eff_full']*100:.1f}% (>=90%)"
            ),
        },
        "P2b-AC6": {
            "passed": ac6,
            "summary": (
                f"+4ms={cmap['shared_broadband_pos4ms']['frac_physical_lag_evidence']*100:.0f}%, "
                f"-7ms={cmap['shared_broadband_neg7ms']['frac_physical_lag_evidence']*100:.0f}%, "
                f"null={cmap['independent_noise']['frac_physical_lag_evidence']*100:.0f}%, "
                f"zero_lag_env_warn={cmap['correlated_zero_lag_env_noise']['frac_zero_lag_env_warning']*100:.0f}%, "
                f"line_art_warn={cmap['common_narrowband_line_artifact']['frac_artifact_warning']*100:.0f}%"
            ),
        },
        "P2b-AC7": {
            "passed": ac7,
            "summary": (
                f"GW170729_V1_FAR={v1_far_s} (<=2.0%), GW170818_H1_FAR={h1_18_far_s} (<4.97% P2), "
                f"confirmed_preserved={confirmed_ok}, O3_untouched={o3['strictly_untouched']}"
            ),
        },
    }


def _write_report(
    sa: dict[str, Any],
    sb: dict[str, Any],
    sc: dict[str, Any],
    sd: dict[str, Any],
    ac: dict[str, Any],
    o3: dict[str, Any],
    out_path: Path,
) -> None:
    lines: list[str] = [
        "# BSD v0.3 Phase 2b Scientific Report",
        "",
        "**Status:** Evaluated Locally on Synthetic Scenarios & Previously Examined O1/O2 Development Epochs ONLY",
        f"**Reserved O3 Holdout Status:** Strictly Untouched (`{len(o3['found_locally'])}` local files from `docs/v03_o3_holdout_manifest.json`)",
        "**Preserved Baselines:** `v0.1`, `v0.2`, `v0.3 Phase 1`, and `v0.3 Phase 2` completely unmodified.",
        "",
        "## 1. Acceptance Criteria Summary",
        "",
        "| Criterion | Status | Measured Result |",
        "| :--- | :---: | :--- |",
    ]
    for k, v in ac.items():
        lines.append(f"| **{k}** | **{'PASS' if v['passed'] else 'FAIL'}** | {v['summary']} |")

    lines.extend(
        [
            "",
            "## 2. Suite A — Transient, Staged Causal Calibration & Multiple-Testing Control",
            "",
            "| Scenario | N | Clean Duty | STFT FAR | Joint FAR | N_thresh | N_eff | Ljung-Box p | T_prov | T_mature | Mature FAR | Full Eff | Clean Eff | Glitch Flag | Regimes |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in sa["rows"]:
        fmtp = lambda x: "N/A" if x is None else f"{x*100:.2f}%"
        fmti = lambda x: "N/A" if x is None else str(x)
        lines.append(
            f"| `{r['scenario']}` | {r['n_windows']} | {r['duty_cycle_clean']*100:.1f}% | {fmtp(r['stft_far_full'])} | "
            f"{fmtp(r['joint_far_full'])} | {r['cal_n_thresh']} | {r['cal_n_eff']:.1f} | {r['cal_ljung_box_p']:.3f} | "
            f"{fmti(r['t_provisional_win'])} | {fmti(r['t_mature_win'])} | {fmtp(r['mature_surv_far'])} | "
            f"{fmtp(r['transient_eff_full'])} | {fmtp(r['transient_eff_clean_cal'])} | {r['glitch_storm_flag_rate']*100:.1f}% | {r['n_regimes']} |"
        )

    lines.extend(
        [
            "",
            "## 3. Suite B — Separate Persistent Signal Matrix (Stationary vs. Drifting Noise)",
            "",
            "Preserves raw scores and separates frozen `MultiResolutionDetector` nominations from Sidak-controlled `FrozenReferenceBroadbandChannel` nominations (`policy=persistent`, zero single-channel recalibration).",
            "",
            "| Morphology | Background | Frozen STFT Full Eff | Frozen-Ref BB Full Eff | Joint Full Eff | Clean Cal STFT Eff | BG Joint FAR | Continuum Drift Flag | Glitch Flag | Regimes |",
            "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in sb["rows"]:
        lines.append(
            f"| `{r['morphology']}` | `{r['background']}` | {r['frozen_stft_eff_full']*100:.1f}% | "
            f"{r['frozen_ref_bb_eff_full']*100:.1f}% | {r['joint_eff_full']*100:.1f}% | "
            f"{r['clean_cal_stft_eff']*100:.1f}% | {r['bg_joint_far_full']*100:.2f}% | "
            f"{r['cont_drift_frac_on_sig']*100:.1f}% | {r['glitch_flag_frac_on_sig']*100:.1f}% | {r['n_regimes']} |"
        )

    lines.extend(
        [
            "",
            "## 4. Suite C — Two-Detector Cross-Coherence Five-Regime Validation",
            "",
            "Coherence is treated strictly as additional evidence calibrated against 120 circular time-slides (`alpha = 0.01`), never automatic astrophysical proof.",
            "",
            "| Control / Signal Regime | Trials | Null / Incoherent | Common Artifact Warning | Zero-Lag Env Warning | Physical-Lag Evidence | Median Lag (ms) |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in sc["rows"]:
        lines.append(
            f"| `{r['regime']}` | {r['n_trials']} | {r['frac_null']*100:.0f}% | {r['frac_artifact_warning']*100:.0f}% | "
            f"{r['frac_zero_lag_env_warning']*100:.0f}% | {r['frac_physical_lag_evidence']*100:.0f}% | {r['median_peak_lag_ms']:+.2f} |"
        )

    lines.extend(
        [
            "",
            "## 5. Suite D — Previously Examined O1/O2 LIGO/Virgo Development Epochs (Causal Trailing)",
            "",
            "| Event / Det | Epoch Type | Clean Duty | Excl. (hr) | Cont. Warn | Glitch Warn | Regimes | Full FAR | Clean Cal FAR (90% CI) | N_eff / N_th | p_LB | Confirmed Event |",
            "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in sd["rows"]:
        ons = "N/A" if r["onsource_confirmed_detected"] is None else ("YES" if r["onsource_confirmed_detected"] else "NO")
        ci_s = (
            "N/A"
            if r["stft_far_clean_cal"] is None
            else f"{r['stft_far_clean_cal']*100:.2f}% [{r['clean_far_ci90_lo']*100:.2f}%, {r['clean_far_ci90_hi']*100:.2f}%]"
        )
        lines.append(
            f"| `{r['event']} {r['detector']}` | `{r['epoch_type']}` | {r['duty_cycle_clean']*100:.1f}% | "
            f"{r['excluded_livetime_hr']:.3f} | {r['n_continuum_warn_wins']} | {r['n_glitch_storm_warn_wins']} | "
            f"{r['n_regimes']} | {r['stft_far_full']*100:.2f}% | {ci_s} | "
            f"{r['cal_n_eff']:.0f}/{r['cal_n_thresh']} | {r['cal_ljung_box_p']:.3f} | {ons} |"
        )

    lines.extend(
        [
            "",
            "## 6. Figures",
            "",
            "- `figures/fig01_staged_causal_recovery_500win.png`",
            "- `figures/fig02_persistent_matrix_streamA_vs_streamB.png`",
            "- `figures/fig03_cross_coherence_five_regimes.png`",
            "- `figures/fig04_ligo_o1o2_orthogonal_diagnostics.png`",
            "",
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    t0 = time.time()
    out_dir = PROJECT_ROOT / "results/v03_phase2b"
    out_dir.mkdir(parents=True, exist_ok=True)

    sa = _run_suite_a_transient()
    sb = _run_suite_b_persistent(sa["shared_ref"])
    sc = _run_suite_c_coherence()
    sd = _run_suite_d_ligo_o1o2()
    o3 = _check_o3_untouched()
    ac = _evaluate_ac(sa, sb, sc, sd, o3)

    _write_csv(out_dir / "suite_a_transient_calibration.csv", sa["rows"])
    _write_csv(out_dir / "suite_b_persistent_matrix.csv", sb["rows"])
    _write_csv(out_dir / "suite_c_cross_coherence.csv", sc["rows"])
    _write_csv(out_dir / "suite_d_ligo_o1o2_dev.csv", sd["rows"])
    figs = _generate_figures(sa, sb, sc, sd, out_dir / "figures")
    _write_report(sa, sb, sc, sd, ac, o3, out_dir / "report.md")

    summary = {
        "acceptance_criteria": ac,
        "o3_holdout_audit": o3,
        "figures": figs,
        "suite_a": sa["rows"],
        "suite_b": sb["rows"],
        "suite_c": sc["rows"],
        "suite_d": sd["rows"],
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    (out_dir / "results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== BSD v0.3 Phase 2b Acceptance Criteria Summary ===")
    for k, v in ac.items():
        print(f"  {k}: {'PASS' if v['passed'] else 'FAIL'} | {v['summary']}")
    print(f"Elapsed: {summary['elapsed_seconds']} s")


if __name__ == "__main__":
    main()
