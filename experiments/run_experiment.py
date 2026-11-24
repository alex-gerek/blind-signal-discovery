"""Main reproducible experiment runner for Blind Signal Discovery (BSD) v0.1.

Usage:
    python experiments/run_experiment.py --config configs/smoke.toml
    python experiments/run_experiment.py --config configs/research.toml
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

# Ensure src/ is importable when run directly from repository root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bsd.calibration import CalibrationResult, calibrate_detector, clopper_pearson_ci, wilson_ci
from bsd.detector import DetectionResult, DetectorConfig, TemplateFreeDetector
from bsd.evaluation import (
    SampleEvaluation,
    evaluate_split,
    group_and_summarize,
    run_robustness_suite,
    summarize_records,
)
from bsd.generator import (
    GroundTruth,
    Observation,
    generate_noise_dataset,
    generate_sample,
    generate_signal_dataset,
    make_split_rng,
)
from bsd.reporting import generate_markdown_report
from bsd.visualization import (
    plot_calibration_distribution,
    plot_candidate_overlays,
    plot_confusion_matrix,
    plot_detection_vs_snr,
    plot_failure_examples,
    plot_false_alarm_results,
    plot_spectrograms_and_excess_maps,
    plot_time_domain_examples,
)


def _merge_binary_summaries(summaries: list[dict[str, Any]], target_fap: float) -> dict[str, Any]:
    """Pool raw TP/FP/TN/FN counts across seeds and recompute exact CIs."""
    tp = sum(s["tp"] for s in summaries)
    fp = sum(s["fp"] for s in summaries)
    tn = sum(s["tn"] for s in summaries)
    fn = sum(s["fn"] for s in summaries)
    n_sig = tp + fn
    n_noise = fp + tn

    det_prob = tp / n_sig if n_sig > 0 else None
    det_w = list(wilson_ci(tp, n_sig)) if n_sig > 0 else None
    det_cp = list(clopper_pearson_ci(tp, n_sig)) if n_sig > 0 else None

    far = fp / n_noise if n_noise > 0 else None
    far_w = list(wilson_ci(fp, n_noise)) if n_noise > 0 else None
    far_cp = list(clopper_pearson_ci(fp, n_noise)) if n_noise > 0 else None

    warnings: list[str] = []
    if far_cp is not None and far_cp[0] > target_fap:
        warnings.append(
            f"Pooled FAR ({far:.4f}, 95% exact CI [{far_cp[0]:.4f}, {far_cp[1]:.4f}]) "
            f"exceeds calibrated target FAP ({target_fap:.4f})."
        )

    out: dict[str, Any] = {
        "n_total": tp + fp + tn + fn,
        "n_signal": n_sig,
        "n_noise": n_noise,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "detection_probability": det_prob,
        "detection_wilson_95ci": det_w,
        "detection_clopper_pearson_95ci": det_cp,
        "false_alarm_rate": far,
        "far_wilson_95ci": far_w,
        "far_clopper_pearson_95ci": far_cp,
        "statistically_credible": len(warnings) == 0,
        "credibility_warnings": warnings,
    }
    if "calibrated_threshold" in summaries[0]:
        out["calibrated_threshold"] = float(np.mean([s["calibrated_threshold"] for s in summaries]))
    return out


def _pool_robustness(seed_robustness: list[dict[str, Any]], target_fap: float) -> dict[str, Any]:
    pooled: dict[str, Any] = {}
    for category in seed_robustness[0]:
        pooled[category] = {}
        for subkey in seed_robustness[0][category]:
            pooled[category][subkey] = _merge_binary_summaries(
                [sr[category][subkey] for sr in seed_robustness],
                target_fap=target_fap,
            )
    return pooled


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_experiment(config_path: Path) -> dict[str, Any]:
    t0 = time.perf_counter()
    with config_path.open("rb") as f:
        cfg = tomllib.load(f)

    exp_cfg = cfg["experiment"]
    sig_cfg = cfg["signal"]
    noise_cfg = cfg["noise"]
    det_cfg_raw = cfg["detector"]
    cal_cfg = cfg["calibration"]
    eval_cfg = cfg["evaluation"]

    exp_name: str = exp_cfg["name"]
    seeds: list[int] = [int(s) for s in exp_cfg["seeds"]]
    output_dir = (PROJECT_ROOT / exp_cfg["output_dir"]).resolve()
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    sample_rate = float(sig_cfg["sample_rate"])
    duration = float(sig_cfg["duration"])
    snr_levels = [float(s) for s in sig_cfg["snr_levels"]]
    known_families = list(sig_cfg["known_families"])
    withheld_families = list(sig_cfg["withheld_families"])

    cal_noise_models = tuple(noise_cfg["calibration_models"])
    eval_noise_models = tuple(noise_cfg["evaluation_models"])
    shifted_noise_models = tuple(noise_cfg["shifted_models"])

    det_config = DetectorConfig(
        nperseg=int(det_cfg_raw["nperseg"]),
        noverlap=int(det_cfg_raw["noverlap"]),
        window=str(det_cfg_raw["window"]),
        freq_median_kernel=int(det_cfg_raw["freq_median_kernel"]),
        time_median_kernel=int(det_cfg_raw["time_median_kernel"]),
        bin_seed_threshold=float(det_cfg_raw["bin_seed_threshold"]),
        min_cluster_bins=int(det_cfg_raw["min_cluster_bins"]),
    )
    detector = TemplateFreeDetector(det_config)
    target_fap = float(cal_cfg["target_fap"])
    n_cal = int(cal_cfg["n_calibration_samples"])

    print(f"[BSD] Starting experiment '{exp_name}' across {len(seeds)} seeds: {seeds}")

    all_noise_recs: list[SampleEvaluation] = []
    all_signal_recs: list[SampleEvaluation] = []
    all_withheld_recs: list[SampleEvaluation] = []
    all_shifted_noise_recs: list[SampleEvaluation] = []
    all_shifted_signal_recs: list[SampleEvaluation] = []
    calibrations: list[CalibrationResult] = []
    per_seed_summary: list[dict[str, Any]] = []
    seed_robustness_list: list[dict[str, Any]] = []

    # Visual examples collected from first seed only
    fp_visual_examples: list[tuple[str, Observation, GroundTruth, DetectionResult]] = []
    fn_visual_examples: list[tuple[str, Observation, GroundTruth, DetectionResult]] = []
    first_calibration: CalibrationResult | None = None

    for seed_idx, seed in enumerate(seeds):
        print(f"[BSD] --- Seed {seed} ({seed_idx + 1}/{len(seeds)}) ---")

        # 1. Dedicated calibration split -> freeze threshold
        cal_split = generate_noise_dataset(
            split_name="calibration",
            n_samples_count=n_cal,
            sample_rate=sample_rate,
            duration=duration,
            noise_models=cal_noise_models,
            master_seed=seed,
        )
        cal = calibrate_detector(
            detector=detector,
            calibration_observations=cal_split.observations,
            target_fap=target_fap,
            calibration_noise_models=cal_noise_models,
            calibration_seed=seed,
        )
        calibrations.append(cal)
        if first_calibration is None:
            first_calibration = cal
        print(
            f"[BSD]   Calibrated 1% FAP threshold = {cal.threshold:.2f} "
            f"(N_cal={n_cal}, emp_cal_fap={cal.empirical_calibration_fap*100:.2f}%)"
        )

        # 2. Strictly independent held-out test splits (no tuning after calibration)
        test_noise_split = generate_noise_dataset(
            split_name="test_noise",
            n_samples_count=int(eval_cfg["n_noise_test"]),
            sample_rate=sample_rate,
            duration=duration,
            noise_models=eval_noise_models,
            master_seed=seed,
        )
        test_signal_split = generate_signal_dataset(
            split_name="test_signal",
            n_samples_count=int(eval_cfg["n_signal_test"]),
            sample_rate=sample_rate,
            duration=duration,
            signal_families=known_families,
            snr_levels=snr_levels,
            noise_models=eval_noise_models,
            master_seed=seed,
        )
        test_withheld_split = generate_signal_dataset(
            split_name="test_withheld",
            n_samples_count=int(eval_cfg["n_withheld_test"]),
            sample_rate=sample_rate,
            duration=duration,
            signal_families=withheld_families,
            snr_levels=snr_levels,
            noise_models=eval_noise_models,
            master_seed=seed,
        )
        test_shifted_noise_split = generate_noise_dataset(
            split_name="test_shifted_noise",
            n_samples_count=int(eval_cfg["n_shifted_noise"]),
            sample_rate=sample_rate,
            duration=duration,
            noise_models=shifted_noise_models,
            master_seed=seed,
        )
        test_shifted_signal_split = generate_signal_dataset(
            split_name="test_shifted_signal",
            n_samples_count=int(eval_cfg["n_shifted_signal"]),
            sample_rate=sample_rate,
            duration=duration,
            signal_families=known_families,
            snr_levels=snr_levels,
            noise_models=shifted_noise_models,
            master_seed=seed,
        )

        r_noise = evaluate_split(detector, test_noise_split, cal)
        r_sig = evaluate_split(detector, test_signal_split, cal)
        r_withheld = evaluate_split(detector, test_withheld_split, cal)
        r_sh_noise = evaluate_split(detector, test_shifted_noise_split, cal)
        r_sh_sig = evaluate_split(detector, test_shifted_signal_split, cal)

        all_noise_recs.extend(r_noise)
        all_signal_recs.extend(r_sig)
        all_withheld_recs.extend(r_withheld)
        all_shifted_noise_recs.extend(r_sh_noise)
        all_shifted_signal_recs.extend(r_sh_sig)

        # Collect diagnostic FP and FN visual examples from the first seed
        if seed_idx == 0:
            for rec, obs, gt in zip(r_noise, test_noise_split.observations, test_noise_split.ground_truths):
                if rec.outcome == "FP" and len(fp_visual_examples) < 2:
                    det_maps = detector.analyze(obs, calibration=cal, return_maps=True)
                    fp_visual_examples.append((f"FP ({obs.sample_id})", obs, gt, det_maps))
            # Also capture one impulsive FP example if < 3 FP examples
            rng_imp_ex = make_split_rng(seed, "visual_impulsive_fp")
            for k in range(10):
                if len(fp_visual_examples) >= 3:
                    break
                o_i, g_i = generate_sample(f"imp_ex_{k}", sample_rate, duration, "impulsive", rng_imp_ex)
                d_i = detector.analyze(o_i, calibration=cal, return_maps=True)
                if d_i.detected:
                    fp_visual_examples.append((f"FP Impulsive ({o_i.sample_id})", o_i, g_i, d_i))
                    break

            for rec, obs, gt in zip(r_sig, test_signal_split.observations, test_signal_split.ground_truths):
                if rec.outcome == "FN" and len(fn_visual_examples) < 3:
                    # Prefer moderate SNR missed signals so they are informative
                    if gt.target_snr is not None and gt.target_snr >= 6.0:
                        det_maps = detector.analyze(obs, calibration=cal, return_maps=True)
                        fn_visual_examples.append(
                            (f"FN {gt.signal_family} SNR={gt.target_snr}", obs, gt, det_maps)
                        )

        # 3. Robustness stress suite
        rob_seed = run_robustness_suite(
            detector=detector,
            calibration=cal,
            calibration_split=cal_split,
            sample_rate=sample_rate,
            duration=duration,
            master_seed=seed,
            n_per_condition=int(eval_cfg["n_robustness_per_condition"]),
        )
        seed_robustness_list.append(rob_seed)

        s_noise_sum = summarize_records(r_noise, target_fap=target_fap)
        s_sig_sum = summarize_records(r_sig, target_fap=target_fap)
        s_wh_sum = summarize_records(r_withheld, target_fap=target_fap)
        per_seed_summary.append(
            {
                "seed": seed,
                "threshold": cal.threshold,
                "noise_fp": s_noise_sum["fp"],
                "noise_n": s_noise_sum["n_noise"],
                "noise_far": s_noise_sum["false_alarm_rate"],
                "known_tp": s_sig_sum["tp"],
                "known_n": s_sig_sum["n_signal"],
                "known_recall": s_sig_sum["detection_probability"],
                "withheld_tp": s_wh_sum["tp"],
                "withheld_n": s_wh_sum["n_signal"],
                "withheld_recall": s_wh_sum["detection_probability"],
            }
        )
        print(
            f"[BSD]   Held-out Noise FAR = {s_noise_sum['false_alarm_rate']*100:.2f}% "
            f"({s_noise_sum['fp']}/{s_noise_sum['n_noise']}) | "
            f"Known Recall = {s_sig_sum['detection_probability']*100:.2f}% | "
            f"Withheld FM Recall = {s_wh_sum['detection_probability']*100:.2f}%"
        )

    # Aggregate results across all seeds
    pooled_in_dist = summarize_records(all_noise_recs + all_signal_recs, target_fap=target_fap)
    pooled_noise_only = summarize_records(all_noise_recs, target_fap=target_fap)
    pooled_noise_by_model = group_and_summarize(all_noise_recs, lambda r: r.noise_model, target_fap=target_fap)
    pooled_signal_by_family = group_and_summarize(all_signal_recs, lambda r: r.signal_family, target_fap=target_fap)
    pooled_withheld = summarize_records(all_withheld_recs, target_fap=target_fap)
    pooled_shifted_noise = summarize_records(all_shifted_noise_recs, target_fap=target_fap)
    pooled_shifted_signal = summarize_records(all_shifted_signal_recs, target_fap=target_fap)
    pooled_robustness = _pool_robustness(seed_robustness_list, target_fap=target_fap)

    # Build family x SNR detection curves (known + withheld)
    family_snr_curves: dict[str, dict[float, dict[str, Any]]] = {}
    for fam in known_families + withheld_families:
        src = all_withheld_recs if fam in withheld_families else all_signal_recs
        fam_recs = [r for r in src if r.signal_family == fam]
        family_snr_curves[fam] = {}
        for snr in snr_levels:
            snr_recs = [r for r in fam_recs if r.target_snr == snr]
            if snr_recs:
                family_snr_curves[fam][snr] = summarize_records(snr_recs, target_fap=target_fap)

    thresholds = [c.threshold for c in calibrations]
    assert first_calibration is not None
    calibration_summary = {
        "target_fap": target_fap,
        "n_calibration_per_seed": n_cal,
        "calibration_noise_models": list(cal_noise_models),
        "total_tf_bins_per_sample": first_calibration.total_tf_bins_per_sample,
        "naive_per_bin_exp_threshold": first_calibration.naive_per_bin_exp_threshold,
        "bonferroni_per_bin_exp_threshold": first_calibration.bonferroni_per_bin_exp_threshold,
        "threshold_mean": float(np.mean(thresholds)),
        "threshold_std": float(np.std(thresholds)),
        "threshold_min": float(np.min(thresholds)),
        "threshold_max": float(np.max(thresholds)),
        "min_resolvable_p_value": 1.0 / (n_cal + 1.0),
        "per_seed_calibrations": [c.to_summary_dict() for c in calibrations],
        "limitations": list(first_calibration.limitations),
    }

    # Generate clean illustrative examples for figures 1-3 using a dedicated visual seed
    rng_vis = make_split_rng(seeds[0], "visual_gallery")
    td_examples: list[tuple[str, Observation, GroundTruth]] = []
    tf_examples: list[tuple[str, Observation, GroundTruth, DetectionResult]] = []

    obs_n0, gt_n0 = generate_sample("vis_noise", sample_rate, duration, "pink_colored", rng_vis)
    td_examples.append(("Noise-Only", obs_n0, gt_n0))

    vis_specs = [
        ("Sinusoid (known)", "sinusoid", "pink_colored", 16.0),
        ("Broadband Pulse (known)", "broadband_pulse", "gaussian_white", 16.0),
        ("Linear Chirp (known)", "linear_chirp", "nonstationary", 16.0),
        ("Intermittent Narrowband (known)", "intermittent_narrowband", "pink_colored", 16.0),
        ("Sinusoidal FM (WITHHELD)", "sinusoidal_fm", "nonstationary", 16.0),
    ]
    for title, fam, nm, snr in vis_specs:
        o_v, g_v = generate_sample(f"vis_{fam}", sample_rate, duration, nm, rng_vis, signal_family=fam, target_snr=snr)
        d_v = detector.analyze(o_v, calibration=first_calibration, return_maps=True)
        td_examples.append((title, o_v, g_v))
        tf_examples.append((title, o_v, g_v, d_v))

    # Generate all 9 publication-quality figures
    generated_figures: list[str] = []
    fig_calls = [
        ("fig01_time_domain_examples.png", lambda p: plot_time_domain_examples(td_examples, p)),
        ("fig02_spectrograms_and_excess_maps.png", lambda p: plot_spectrograms_and_excess_maps(tf_examples[:4], p)),
        ("fig03_candidate_overlays.png", lambda p: plot_candidate_overlays(tf_examples[1:], p)),
        ("fig04_detection_vs_snr.png", lambda p: plot_detection_vs_snr(family_snr_curves, target_fap, p)),
        (
            "fig05_false_alarm_results.png",
            lambda p: plot_false_alarm_results(
                {
                    "Pooled In-Dist": pooled_noise_only,
                    "Gaussian White": pooled_noise_by_model["gaussian_white"],
                    "Pink Colored": pooled_noise_by_model["pink_colored"],
                    "Nonstationary": pooled_noise_by_model["nonstationary"],
                    "Shifted AR(2)": pooled_shifted_noise,
                    "Strong Nonstat.": pooled_robustness["strong_nonstationary"]["noise_only"],
                    "Impulsive": pooled_robustness["impulsive_interference"]["noise_only"],
                },
                target_fap,
                p,
            ),
        ),
        ("fig06_confusion_matrix.png", lambda p: plot_confusion_matrix(pooled_in_dist, p)),
        (
            "fig07_false_positive_examples.png",
            lambda p: plot_failure_examples(
                fp_visual_examples,
                "Measured False-Positive Examples (Simulated Noise-Only Observations)",
                p,
            ),
        ),
        (
            "fig08_missed_signal_examples.png",
            lambda p: plot_failure_examples(
                fn_visual_examples,
                "Measured Missed-Signal (False-Negative) Examples (Simulated Injected Signals)",
                p,
            ),
        ),
        ("fig09_calibration_null_distribution.png", lambda p: plot_calibration_distribution(first_calibration, p)),
    ]
    for fname, fn in fig_calls:
        fn(figures_dir / fname)
        generated_figures.append(fname)

    elapsed_s = time.perf_counter() - t0
    results: dict[str, Any] = {
        "metadata": {
            "experiment_name": exp_name,
            "config_path": str(config_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
            "seeds": seeds,
            "sample_rate": sample_rate,
            "duration": duration,
            "snr_levels": snr_levels,
            "known_families": known_families,
            "withheld_families": withheld_families,
            "detector_config": det_config.to_dict(),
            "runtime_seconds": round(elapsed_s, 2),
        },
        "calibration_summary": calibration_summary,
        "per_seed_summary": per_seed_summary,
        "held_out_in_distribution": pooled_in_dist,
        "held_out_noise_only": pooled_noise_only,
        "held_out_noise_by_model": pooled_noise_by_model,
        "held_out_signal_by_family": pooled_signal_by_family,
        "held_out_withheld_family": pooled_withheld,
        "held_out_shifted_noise": pooled_shifted_noise,
        "held_out_shifted_signal": pooled_shifted_signal,
        "family_snr_curves": {
            fam: {str(k): v for k, v in snr_dict.items()}
            for fam, snr_dict in family_snr_curves.items()
        },
        "robustness_pooled": pooled_robustness,
        "figures": generated_figures,
    }

    # Save machine-readable JSON and CSV outputs
    json_path = output_dir / "results.json"
    json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    all_eval_rows = [
        r.to_dict()
        for r in (
            all_noise_recs
            + all_signal_recs
            + all_withheld_recs
            + all_shifted_noise_recs
            + all_shifted_signal_recs
        )
    ]
    _write_csv(output_dir / "per_sample_results.csv", all_eval_rows)

    curve_rows: list[dict[str, Any]] = []
    for fam, snr_dict in family_snr_curves.items():
        for snr, stats_d in snr_dict.items():
            curve_rows.append(
                {
                    "signal_family": fam,
                    "is_withheld": fam in withheld_families,
                    "target_snr": snr,
                    "n_signal": stats_d["n_signal"],
                    "tp": stats_d["tp"],
                    "detection_probability": stats_d["detection_probability"],
                    "wilson_ci_low": stats_d["detection_wilson_95ci"][0],
                    "wilson_ci_high": stats_d["detection_wilson_95ci"][1],
                }
            )
    _write_csv(output_dir / "summary_by_family_snr.csv", curve_rows)

    # Generate English research report
    report_path = generate_markdown_report(results, output_dir / "report.md")

    print(f"[BSD] Completed in {elapsed_s:.2f} s.")
    print(f"[BSD] Pooled Held-Out Noise FAR: {pooled_noise_only['false_alarm_rate']*100:.2f}% ({pooled_noise_only['fp']}/{pooled_noise_only['n_noise']})")
    print(f"[BSD] Pooled Held-Out Known-Signal Recall: {pooled_in_dist['detection_probability']*100:.2f}% ({pooled_in_dist['tp']}/{pooled_in_dist['n_signal']})")
    print(f"[BSD] Pooled Withheld FM Recall: {pooled_withheld['detection_probability']*100:.2f}% ({pooled_withheld['tp']}/{pooled_withheld['n_signal']})")
    print(f"[BSD] Shifted AR(2) Noise FAR: {pooled_shifted_noise['false_alarm_rate']*100:.2f}% | Credible: {pooled_shifted_noise['statistically_credible']}")
    imp_far = pooled_robustness["impulsive_interference"]["noise_only"]
    print(f"[BSD] Impulsive Noise FAR: {imp_far['false_alarm_rate']*100:.2f}% | Credible: {imp_far['statistically_credible']}")
    print(f"[BSD] Report written to: {report_path}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Blind Signal Discovery (BSD) v0.1 experiment.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/smoke.toml"),
        help="Path to TOML experiment configuration file (relative to repo root or absolute).",
    )
    args = parser.parse_args()
    cfg_path = args.config if args.config.is_absolute() else (PROJECT_ROOT / args.config)
    run_experiment(cfg_path)


if __name__ == "__main__":
    main()
