"""Independent blind evaluation, localization metrics, and robustness suite for BSD.

Ground-truth metadata is accessed ONLY after the detector has produced its
predictions on raw `Observation` objects using a frozen `CalibrationResult`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from bsd.calibration import CalibrationResult, calibrate_detector, clopper_pearson_ci, wilson_ci
from bsd.detector import CandidateEvent, DetectorConfig, TemplateFreeDetector
from bsd.generator import (
    DatasetSplit,
    GroundTruth,
    generate_sample,
    make_split_rng,
)


@dataclass(frozen=True)
class SampleEvaluation:
    """Evaluation record for one blind test sample."""

    sample_id: str
    split_name: str
    seed: int
    is_signal: bool
    detected: bool
    outcome: str  # TP, FP, TN, FN
    signal_family: str | None
    noise_model: str
    target_snr: float | None
    empirical_energy_snr: float | None
    max_anomaly_score: float
    threshold_used: float
    n_candidates: int
    top_p_value: float | None
    p_value_is_upper_bound: bool
    time_iou: float | None
    freq_iou: float | None
    time_center_error_s: float | None
    freq_center_error_hz: float | None
    tf_overlap: bool | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _interval_iou(a0: float, a1: float, b0: float, b1: float) -> float:
    inter = max(0.0, min(a1, b1) - max(a0, b0))
    union = max(a1, b1) - min(a0, b0)
    if union <= 0.0:
        return 0.0
    return float(inter / union)


def _compute_localization(
    candidates: list[CandidateEvent],
    gt: GroundTruth,
) -> tuple[float | None, float | None, float | None, float | None, bool | None]:
    """Compute time/frequency localization metrics between candidates and ground truth."""
    if not gt.is_signal or not candidates:
        return (None, None, None, None, None)
    assert gt.t_start is not None and gt.t_end is not None
    assert gt.f_min is not None and gt.f_max is not None

    # Use union of detected candidates for multi-burst/fragmented signals,
    # and also check direct TF box overlap on any individual candidate.
    c_t0 = min(c.t_start for c in candidates)
    c_t1 = max(c.t_end for c in candidates)
    c_f0 = min(c.f_min for c in candidates)
    c_f1 = max(c.f_max for c in candidates)

    t_iou = _interval_iou(c_t0, c_t1, gt.t_start, gt.t_end)
    f_iou = _interval_iou(c_f0, c_f1, gt.f_min, gt.f_max)

    gt_tc = 0.5 * (gt.t_start + gt.t_end)
    gt_fc = 0.5 * (gt.f_min + gt.f_max)
    top = candidates[0]
    t_err = abs(top.peak_time - gt_tc)
    f_err = abs(top.peak_freq - gt_fc)

    any_overlap = any(
        (_interval_iou(c.t_start, c.t_end, gt.t_start, gt.t_end) > 0.0)
        and (_interval_iou(c.f_min, c.f_max, gt.f_min, gt.f_max) > 0.0)
        for c in candidates
    )
    return (t_iou, f_iou, float(t_err), float(f_err), any_overlap)


def evaluate_split(
    detector: TemplateFreeDetector,
    split: DatasetSplit,
    calibration: CalibrationResult,
) -> list[SampleEvaluation]:
    """Evaluate a frozen detector and calibration on an independent DatasetSplit."""
    records: list[SampleEvaluation] = []
    for obs, gt in zip(split.observations, split.ground_truths):
        if obs.sample_id != gt.sample_id:
            raise ValueError(f"Mismatched sample IDs: {obs.sample_id} != {gt.sample_id}")

        # Strictly blind call: detector sees only `obs` and frozen `calibration`
        det = detector.analyze(obs, calibration=calibration, return_maps=False)

        if gt.is_signal and det.detected:
            outcome = "TP"
        elif gt.is_signal and not det.detected:
            outcome = "FN"
        elif not gt.is_signal and det.detected:
            outcome = "FP"
        else:
            outcome = "TN"

        t_iou, f_iou, t_err, f_err, tf_overlap = _compute_localization(det.candidates, gt)
        top_p = det.candidates[0].empirical_p_value if det.candidates else None
        top_ub = det.candidates[0].p_value_is_upper_bound if det.candidates else False

        records.append(
            SampleEvaluation(
                sample_id=obs.sample_id,
                split_name=split.split_name,
                seed=split.seed,
                is_signal=gt.is_signal,
                detected=det.detected,
                outcome=outcome,
                signal_family=gt.signal_family,
                noise_model=gt.noise_model,
                target_snr=gt.target_snr,
                empirical_energy_snr=gt.empirical_energy_snr,
                max_anomaly_score=det.max_anomaly_score,
                threshold_used=calibration.threshold,
                n_candidates=len(det.candidates),
                top_p_value=top_p,
                p_value_is_upper_bound=top_ub,
                time_iou=t_iou,
                freq_iou=f_iou,
                time_center_error_s=t_err,
                freq_center_error_hz=f_err,
                tf_overlap=tf_overlap,
            )
        )
    return records


def summarize_records(
    records: list[SampleEvaluation],
    target_fap: float = 0.01,
) -> dict[str, Any]:
    """Compute detection probability, FAR, CIs, precision/recall, and credibility flags."""
    tp = sum(1 for r in records if r.outcome == "TP")
    fp = sum(1 for r in records if r.outcome == "FP")
    tn = sum(1 for r in records if r.outcome == "TN")
    fn = sum(1 for r in records if r.outcome == "FN")

    n_signal = tp + fn
    n_noise = fp + tn

    det_prob = tp / n_signal if n_signal > 0 else None
    det_wilson = wilson_ci(tp, n_signal) if n_signal > 0 else None
    det_cp = clopper_pearson_ci(tp, n_signal) if n_signal > 0 else None

    far = fp / n_noise if n_noise > 0 else None
    far_wilson = wilson_ci(fp, n_noise) if n_noise > 0 else None
    far_cp = clopper_pearson_ci(fp, n_noise) if n_noise > 0 else None

    precision = tp / (tp + fp) if (tp + fp) > 0 else None
    recall = det_prob
    f1 = (
        (2.0 * precision * recall / (precision + recall))
        if (precision is not None and recall is not None and (precision + recall) > 0)
        else None
    )

    tp_records = [r for r in records if r.outcome == "TP"]
    t_ious = [r.time_iou for r in tp_records if r.time_iou is not None]
    f_ious = [r.freq_iou for r in tp_records if r.freq_iou is not None]
    t_errs = [r.time_center_error_s for r in tp_records if r.time_center_error_s is not None]
    f_errs = [r.freq_center_error_hz for r in tp_records if r.freq_center_error_hz is not None]
    overlaps = [r.tf_overlap for r in tp_records if r.tf_overlap is not None]

    # Statistical credibility assessment
    credibility_warnings: list[str] = []
    if n_noise > 0 and n_noise < 100:
        credibility_warnings.append(
            f"Small noise sample size (N={n_noise}) limits precision of 1% FAR estimation."
        )
    if far_cp is not None and far_cp[0] > target_fap:
        credibility_warnings.append(
            f"Empirical FAR ({far:.4f}, 95% exact CI [{far_cp[0]:.4f}, {far_cp[1]:.4f}]) "
            f"significantly exceeds calibrated target FAP ({target_fap:.4f}); "
            "statistical significance estimates are NOT credible under this noise condition."
        )
    if n_signal > 0 and n_signal < 30:
        credibility_warnings.append(
            f"Small signal sample size (N={n_signal}) yields wide detection-rate uncertainty."
        )

    return {
        "n_total": len(records),
        "n_signal": n_signal,
        "n_noise": n_noise,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "detection_probability": det_prob,
        "detection_wilson_95ci": list(det_wilson) if det_wilson else None,
        "detection_clopper_pearson_95ci": list(det_cp) if det_cp else None,
        "false_alarm_rate": far,
        "far_wilson_95ci": list(far_wilson) if far_wilson else None,
        "far_clopper_pearson_95ci": list(far_cp) if far_cp else None,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "localization": {
            "n_localized": len(tp_records),
            "median_time_iou": float(np.median(t_ious)) if t_ious else None,
            "median_freq_iou": float(np.median(f_ious)) if f_ious else None,
            "median_time_center_error_s": float(np.median(t_errs)) if t_errs else None,
            "median_freq_center_error_hz": float(np.median(f_errs)) if f_errs else None,
            "tf_overlap_rate": float(np.mean(overlaps)) if overlaps else None,
        },
        "statistically_credible": len(credibility_warnings) == 0,
        "credibility_warnings": credibility_warnings,
    }


def group_and_summarize(
    records: list[SampleEvaluation],
    key_fn: Any,
    target_fap: float = 0.01,
) -> dict[str, dict[str, Any]]:
    """Group evaluation records by a key function and summarize each group."""
    buckets: dict[str, list[SampleEvaluation]] = {}
    for rec in records:
        k = str(key_fn(rec))
        buckets.setdefault(k, []).append(rec)
    return {k: summarize_records(v, target_fap=target_fap) for k, v in sorted(buckets.items())}


def run_robustness_suite(
    detector: TemplateFreeDetector,
    calibration: CalibrationResult,
    calibration_split: DatasetSplit,
    sample_rate: float,
    duration: float,
    master_seed: int,
    n_per_condition: int = 100,
) -> dict[str, Any]:
    """Run targeted robustness stress tests covering 6 key failure/artifact modes."""
    target_fap = calibration.target_fap
    df = sample_rate / detector.config.nperseg
    results: dict[str, Any] = {}

    # 1. Spectral leakage: exact STFT bin center vs half-bin offset sinusoid at SNR=12
    rng_leak = make_split_rng(master_seed, "robustness_leakage")
    on_bin_recs: list[SampleEvaluation] = []
    off_bin_recs: list[SampleEvaluation] = []
    for idx in range(n_per_condition):
        bin_k = int(rng_leak.integers(10, 50))
        f_on = bin_k * df
        f_off = (bin_k + 0.5) * df
        obs_on, gt_on = generate_sample(
            f"leak_on_{idx}",
            sample_rate,
            duration,
            "gaussian_white",
            rng_leak,
            signal_family="sinusoid",
            target_snr=12.0,
            signal_params={"freq_hz": f_on, "taper_alpha": 0.0},
        )
        obs_off, gt_off = generate_sample(
            f"leak_off_{idx}",
            sample_rate,
            duration,
            "gaussian_white",
            rng_leak,
            signal_family="sinusoid",
            target_snr=12.0,
            signal_params={"freq_hz": f_off, "taper_alpha": 0.0},
        )
        on_bin_recs.extend(
            evaluate_split(
                detector,
                DatasetSplit("robustness_leakage_on", master_seed, [obs_on], [gt_on]),
                calibration,
            )
        )
        off_bin_recs.extend(
            evaluate_split(
                detector,
                DatasetSplit("robustness_leakage_off", master_seed, [obs_off], [gt_off]),
                calibration,
            )
        )
    results["spectral_leakage"] = {
        "on_bin_sinusoid_snr12": summarize_records(on_bin_recs, target_fap=target_fap),
        "off_bin_sinusoid_snr12": summarize_records(off_bin_recs, target_fap=target_fap),
    }

    # 2. Windowing effects: compare Hann, Rectangular (boxcar), and Blackman at matched 1% FAP
    rng_win = make_split_rng(master_seed, "robustness_windowing")
    win_obs = []
    win_gts = []
    for idx in range(n_per_condition):
        bin_k = int(rng_win.integers(10, 50))
        obs_w, gt_w = generate_sample(
            f"win_{idx}",
            sample_rate,
            duration,
            "pink_colored",
            rng_win,
            signal_family="sinusoid",
            target_snr=14.0,
            signal_params={"freq_hz": (bin_k + 0.5) * df, "taper_alpha": 0.0},
        )
        win_obs.append(obs_w)
        win_gts.append(gt_w)
    win_split = DatasetSplit("robustness_windowing", master_seed, win_obs, win_gts)

    window_summary: dict[str, Any] = {}
    for win_name in ("hann", "boxcar", "blackman"):
        w_cfg = DetectorConfig(
            nperseg=detector.config.nperseg,
            noverlap=detector.config.noverlap,
            window=win_name,
            freq_median_kernel=detector.config.freq_median_kernel,
            time_median_kernel=detector.config.time_median_kernel,
            bin_seed_threshold=detector.config.bin_seed_threshold,
            min_cluster_bins=detector.config.min_cluster_bins,
        )
        w_det = TemplateFreeDetector(w_cfg)
        w_cal = (
            calibration
            if win_name == detector.config.window
            else calibrate_detector(
                w_det,
                calibration_split.observations,
                target_fap=target_fap,
                calibration_noise_models=calibration.calibration_noise_models,
                calibration_seed=master_seed,
            )
        )
        w_recs = evaluate_split(w_det, win_split, w_cal)
        window_summary[win_name] = {
            "calibrated_threshold": w_cal.threshold,
            **summarize_records(w_recs, target_fap=target_fap),
        }
    results["windowing_effects"] = window_summary

    # 3. Boundary artifacts: pulses near t=0 / t=T vs interior pulses at SNR=12
    rng_bnd = make_split_rng(master_seed, "robustness_boundary")
    bnd_edge_obs, bnd_edge_gts = [], []
    bnd_int_obs, bnd_int_gts = [], []
    for idx in range(n_per_condition):
        edge_t = 0.06 if (idx % 2 == 0) else (duration - 0.06)
        o_e, g_e = generate_sample(
            f"bnd_edge_{idx}",
            sample_rate,
            duration,
            "gaussian_white",
            rng_bnd,
            signal_family="broadband_pulse",
            target_snr=12.0,
            signal_params={"t_center_s": edge_t, "sigma_t_s": 0.004},
        )
        o_i, g_i = generate_sample(
            f"bnd_int_{idx}",
            sample_rate,
            duration,
            "gaussian_white",
            rng_bnd,
            signal_family="broadband_pulse",
            target_snr=12.0,
            signal_params={"t_center_s": float(rng_bnd.uniform(1.0, duration - 1.0)), "sigma_t_s": 0.004},
        )
        bnd_edge_obs.append(o_e)
        bnd_edge_gts.append(g_e)
        bnd_int_obs.append(o_i)
        bnd_int_gts.append(g_i)
    results["boundary_artifacts"] = {
        "edge_pulse_snr12": summarize_records(
            evaluate_split(detector, DatasetSplit("bnd_edge", master_seed, bnd_edge_obs, bnd_edge_gts), calibration),
            target_fap=target_fap,
        ),
        "interior_pulse_snr12": summarize_records(
            evaluate_split(detector, DatasetSplit("bnd_int", master_seed, bnd_int_obs, bnd_int_gts), calibration),
            target_fap=target_fap,
        ),
    }

    # 4. Impulsive interference: noise-only FAR and signal detection
    rng_imp = make_split_rng(master_seed, "robustness_impulsive")
    imp_n_obs, imp_n_gts = [], []
    imp_s_obs, imp_s_gts = [], []
    for idx in range(n_per_condition):
        o_n, g_n = generate_sample(f"imp_n_{idx}", sample_rate, duration, "impulsive", rng_imp)
        o_s, g_s = generate_sample(
            f"imp_s_{idx}",
            sample_rate,
            duration,
            "impulsive",
            rng_imp,
            signal_family="linear_chirp",
            target_snr=12.0,
        )
        imp_n_obs.append(o_n)
        imp_n_gts.append(g_n)
        imp_s_obs.append(o_s)
        imp_s_gts.append(g_s)
    results["impulsive_interference"] = {
        "noise_only": summarize_records(
            evaluate_split(detector, DatasetSplit("imp_noise", master_seed, imp_n_obs, imp_n_gts), calibration),
            target_fap=target_fap,
        ),
        "linear_chirp_snr12": summarize_records(
            evaluate_split(detector, DatasetSplit("imp_sig", master_seed, imp_s_obs, imp_s_gts), calibration),
            target_fap=target_fap,
        ),
    }

    # 5. Strong nonstationary background (modulation_depth=0.85)
    rng_ns = make_split_rng(master_seed, "robustness_nonstationary")
    ns_n_obs, ns_n_gts = [], []
    ns_s_obs, ns_s_gts = [], []
    for idx in range(n_per_condition):
        o_n, g_n = generate_sample(
            f"ns_n_{idx}",
            sample_rate,
            duration,
            "nonstationary",
            rng_ns,
            noise_params={"modulation_depth": 0.85},
        )
        o_s, g_s = generate_sample(
            f"ns_s_{idx}",
            sample_rate,
            duration,
            "nonstationary",
            rng_ns,
            signal_family="intermittent_narrowband",
            target_snr=12.0,
            noise_params={"modulation_depth": 0.85},
        )
        ns_n_obs.append(o_n)
        ns_n_gts.append(g_n)
        ns_s_obs.append(o_s)
        ns_s_gts.append(g_s)
    results["strong_nonstationary"] = {
        "noise_only": summarize_records(
            evaluate_split(detector, DatasetSplit("ns_noise", master_seed, ns_n_obs, ns_n_gts), calibration),
            target_fap=target_fap,
        ),
        "intermittent_snr12": summarize_records(
            evaluate_split(detector, DatasetSplit("ns_sig", master_seed, ns_s_obs, ns_s_gts), calibration),
            target_fap=target_fap,
        ),
    }

    # 6. Weak and partially overlapping signals
    rng_ov = make_split_rng(master_seed, "robustness_overlapping")
    weak_obs, weak_gts = [], []
    ov_obs, ov_gts = [], []
    for idx in range(n_per_condition):
        o_w, g_w = generate_sample(
            f"weak_{idx}",
            sample_rate,
            duration,
            "pink_colored",
            rng_ov,
            signal_family="linear_chirp",
            target_snr=6.0,
        )
        o_o, g_o = generate_sample(
            f"overlap_{idx}",
            sample_rate,
            duration,
            "pink_colored",
            rng_ov,
            signal_family="linear_chirp",
            target_snr=10.0,
            secondary_signal=("intermittent_narrowband", 10.0, None),
        )
        weak_obs.append(o_w)
        weak_gts.append(g_w)
        ov_obs.append(o_o)
        ov_gts.append(g_o)
    results["weak_and_overlapping"] = {
        "weak_chirp_snr6": summarize_records(
            evaluate_split(detector, DatasetSplit("weak_sig", master_seed, weak_obs, weak_gts), calibration),
            target_fap=target_fap,
        ),
        "overlapping_chirp_plus_intermittent_snr10": summarize_records(
            evaluate_split(detector, DatasetSplit("overlap_sig", master_seed, ov_obs, ov_gts), calibration),
            target_fap=target_fap,
        ),
    }

    return results
