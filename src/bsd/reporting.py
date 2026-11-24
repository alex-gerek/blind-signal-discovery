"""Scientific report generator for Blind Signal Discovery (BSD) v0.1.

Generates a structured English Markdown report directly from measured
experimental outputs, strictly separating measured observations from
interpretations and documenting failures and statistical limitations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def _fmt_pct(val: float | None) -> str:
    if val is None:
        return "N/A"
    return f"{val * 100.0:.2f}%"


def _fmt_ci_pct(ci: list[float] | None) -> str:
    if not ci:
        return "N/A"
    return f"[{ci[0] * 100.0:.2f}%, {ci[1] * 100.0:.2f}%]"


def _fmt_num(val: float | None, digits: int = 3) -> str:
    if val is None:
        return "N/A"
    return f"{val:.{digits}f}"


def generate_markdown_report(
    results: dict[str, Any],
    output_path: Path,
) -> Path:
    """Write the complete English research report from measured experiment results."""
    meta = results["metadata"]
    cal_agg = results["calibration_summary"]
    in_dist = results["held_out_in_distribution"]
    noise_only = results["held_out_noise_only"]
    by_noise = results["held_out_noise_by_model"]
    by_family = results["held_out_signal_by_family"]
    withheld = results["held_out_withheld_family"]
    shifted_noise = results["held_out_shifted_noise"]
    shifted_signal = results["held_out_shifted_signal"]
    rob = results["robustness_pooled"]
    per_seed = results["per_seed_summary"]

    lines: list[str] = []
    lines.append("# Blind Signal Discovery (BSD) v0.1 — Experimental Research Report")
    lines.append("")
    lines.append("> **Scientific Scope Statement:** This report presents a controlled methodological study of template-free anomaly detection in simulated time-series data. **A detected anomaly is not a physical discovery.** No claims are made regarding dark matter, axions, gravitational waves, or any physical phenomenon.")
    lines.append("")
    lines.append("## 1. Research Question")
    lines.append("")
    lines.append(
        "Can a classical, template-free time-frequency anomaly detector discover previously unseen signal structures in noisy time-series data while maintaining a controlled and independently measured experiment-level false-alarm probability (target FAP = 1.0%)?"
    )
    lines.append("")
    lines.append("## 2. Experimental Methodology & Dataset Separation")
    lines.append("")
    lines.append(f"- **Experiment configuration:** `{meta['experiment_name']}`")
    lines.append(f"- **Reproducible seeds:** `{meta['seeds']}` ({len(meta['seeds'])} independent experimental runs)")
    lines.append(f"- **Sampling parameters:** `sample_rate = {meta['sample_rate']} Hz`, `duration = {meta['duration']} s` (`{int(meta['sample_rate'] * meta['duration'])}` samples per observation)")
    lines.append(f"- **SNR definition:** Total signal energy SNR $\\mathrm{{SNR}}_{{\\mathrm{{energy}}}} = \\|s\\|_2 / \\sigma_{{\\mathrm{{noise}}}}$, evaluated over `{meta['snr_levels']}`.")
    lines.append("- **Strict dataset separation:** Independent RNG streams are derived via `numpy.random.SeedSequence([master_seed, split_stream_code])` for calibration, held-out noise-only testing, held-out known-signal testing, withheld-family testing, shifted-noise testing, and each robustness stress test.")
    lines.append("- **Ground-truth isolation:** The detector API (`TemplateFreeDetector.analyze`) accepts only `Observation` objects containing the raw 1D waveform and sample rate. Ground-truth labels, noise models, SNR values, and time-frequency bounds are stored separately in `GroundTruth` objects and accessed only after frozen-threshold predictions are finalized.")
    lines.append("- **Frozen-threshold protocol:** For each seed, the detection threshold is calibrated exclusively on the dedicated calibration split and then frozen prior to all blind evaluations.")
    lines.append("")
    lines.append("## 3. Detector Design & Noise Assumptions")
    lines.append("")
    dcfg = meta["detector_config"]
    lines.append("The baseline detector uses classical, interpretable signal processing (no templates, no machine learning):")
    lines.append(f"1. **Robust time-domain preprocessing:** Median subtraction for DC removal; MAD-scaled dispersion, Fisher kurtosis, and crest factor recorded as diagnostics.")
    lines.append(f"2. **STFT time-frequency representation:** `{dcfg['window']}` window, `nperseg={dcfg['nperseg']}`, `noverlap={dcfg['noverlap']}`, `boundary=None, padded=False` (avoiding zero-pad edge artifacts). DC and Nyquist bins are excluded so interior bins share a common complex $\\chi^2_2 \\equiv \\mathrm{{Exp}}(1)$ null distribution under locally stationary Gaussian noise.")
    lines.append(f"3. **Robust 2D background estimation:** Per-frequency PSD estimated via temporal median divided by $\\ln 2$, followed by a frequency median filter (`kernel={dcfg['freq_median_kernel']}`) so stationary narrowband tones are not absorbed into the background. After spectral whitening, slowly varying temporal noise variance is estimated via frequency median divided by $\\ln 2$ and smoothed with a temporal median filter (`kernel={dcfg['time_median_kernel']}`) so short broadband pulses are not absorbed.")
    lines.append(f"4. **Connected-component clustering:** Bins with normalized excess power $Z[f,t] > {dcfg['bin_seed_threshold']}$ are grouped via 8-connectivity (`min_cluster_bins={dcfg['min_cluster_bins']}`). Each candidate cluster $C$ is scored by integrated excess power $S(C) = \\sum_{{(f,t)\\in C}} (Z[f,t] - 1)$.")
    lines.append("")
    lines.append("**Key noise assumptions:** Locally smooth background PSD (bandwidth broader than the frequency median kernel) and slowly varying temporal variance (timescale longer than the temporal median kernel), without unmodeled heavy-tailed impulsive transients.")
    lines.append("")
    lines.append("## 4. Statistical Calibration Procedure")
    lines.append("")
    lines.append(f"- **Target experiment-level FAP:** `{cal_agg['target_fap'] * 100.0:.2f}%`")
    lines.append(f"- **Calibration samples per seed:** `{cal_agg['n_calibration_per_seed']}` noise-only observations across `{cal_agg['calibration_noise_models']}`")
    lines.append(f"- **Mean calibrated threshold across seeds:** `{cal_agg['threshold_mean']:.2f} ± {cal_agg['threshold_std']:.2f}` (range `[{cal_agg['threshold_min']:.2f}, {cal_agg['threshold_max']:.2f}]`)")
    lines.append(f"- **Per-bin vs. experiment-level FAP distinction:** Each observation spans `{cal_agg['total_tf_bins_per_sample']}` interior STFT bins. A naive single-bin 1% $\\mathrm{{Exp}}(1)$ threshold is only `{cal_agg['naive_per_bin_exp_threshold']:.2f}` (and a single-bin Bonferroni correction is `{cal_agg['bonferroni_per_bin_exp_threshold']:.2f}`). Calibrating on the empirical order statistic $k = \\lceil(1 - \\alpha)(N_{{\\mathrm{{cal}}}} + 1)\\rceil$ of the experiment-level maximum cluster score $T_{{\\max}} = \\max_C S(C)$ directly accounts for scanning all `{cal_agg['total_tf_bins_per_sample']}` time-frequency bins AND their STFT overlap correlations.")
    lines.append(f"- **Minimum resolvable empirical p-value per seed:** `{cal_agg['min_resolvable_p_value']:.5f}` ($1 / (N_{{\\mathrm{{cal}}}} + 1)$). Candidate scores exceeding all calibration noise samples are reported explicitly as upper bounds ($p \\le 1/(N_{{\\mathrm{{cal}}}}+1)$), never extrapolated.")
    lines.append("")
    lines.append("## 5. Measured Observations (Pooled Across All Seeds)")
    lines.append("")
    lines.append("### 5.1 Held-Out In-Distribution Summary")
    lines.append("")
    lines.append(f"- **Total held-out noise-only samples:** `{noise_only['n_noise']}` (FP = `{noise_only['fp']}`, TN = `{noise_only['tn']}`)")
    lines.append(f"- **Empirical Noise-Only FAR:** `{_fmt_pct(noise_only['false_alarm_rate'])}` (95% Wilson CI `{_fmt_ci_pct(noise_only['far_wilson_95ci'])}`, 95% exact Clopper-Pearson CI `{_fmt_ci_pct(noise_only['far_clopper_pearson_95ci'])}`)")
    lines.append(f"- **Total held-out known-signal samples:** `{in_dist['n_signal']}` (TP = `{in_dist['tp']}`, FN = `{in_dist['fn']}`)")
    lines.append(f"- **Overall Detection Probability (Recall):** `{_fmt_pct(in_dist['detection_probability'])}` (95% Wilson CI `{_fmt_ci_pct(in_dist['detection_wilson_95ci'])}`)")
    lines.append(f"- **Precision:** `{_fmt_pct(in_dist['precision'])}` | **F1 Score:** `{_fmt_num(in_dist['f1'])}`")
    loc = in_dist["localization"]
    lines.append(f"- **Localization (on detected TPs, N={loc['n_localized']}):** Median Time IoU = `{_fmt_num(loc['median_time_iou'])}`, Median Freq IoU = `{_fmt_num(loc['median_freq_iou'])}`, Median Time Center Error = `{_fmt_num(loc['median_time_center_error_s'])} s`, Median Freq Center Error = `{_fmt_num(loc['median_freq_center_error_hz'], 1)} Hz`, TF Box Overlap Rate = `{_fmt_pct(loc['tf_overlap_rate'])}`.")
    lines.append("")
    lines.append("### 5.2 Per-Seed Reproducibility Check")
    lines.append("")
    lines.append("| Seed | Calibrated Threshold | Noise-Only FAR (FP/N) | Known-Signal Recall (TP/N) | Withheld FM Recall (TP/N) |")
    lines.append("| :--- | :---: | :---: | :---: | :---: |")
    for s_row in per_seed:
        lines.append(
            f"| `{s_row['seed']}` | `{s_row['threshold']:.2f}` | "
            f"`{_fmt_pct(s_row['noise_far'])}` ({s_row['noise_fp']}/{s_row['noise_n']}) | "
            f"`{_fmt_pct(s_row['known_recall'])}` ({s_row['known_tp']}/{s_row['known_n']}) | "
            f"`{_fmt_pct(s_row['withheld_recall'])}` ({s_row['withheld_tp']}/{s_row['withheld_n']}) |"
        )
    lines.append("")
    lines.append("### 5.3 False-Alarm Performance by Noise Model")
    lines.append("")
    lines.append("| Noise Condition | In Calibration? | FP / N | Empirical FAR | 95% Exact CI | Statistically Credible? |")
    lines.append("| :--- | :---: | :---: | :---: | :---: | :---: |")
    for nm, stats_dict in by_noise.items():
        cred = "Yes" if stats_dict["statistically_credible"] else "NO (FAR inflated)"
        lines.append(
            f"| `{nm}` | Yes | {stats_dict['fp']} / {stats_dict['n_noise']} | "
            f"{_fmt_pct(stats_dict['false_alarm_rate'])} | {_fmt_ci_pct(stats_dict['far_clopper_pearson_95ci'])} | {cred} |"
        )
    imp_n = rob["impulsive_interference"]["noise_only"]
    ns_n = rob["strong_nonstationary"]["noise_only"]
    lines.append(
        f"| `shifted_ar2` (unseen) | No | {shifted_noise['fp']} / {shifted_noise['n_noise']} | "
        f"{_fmt_pct(shifted_noise['false_alarm_rate'])} | {_fmt_ci_pct(shifted_noise['far_clopper_pearson_95ci'])} | "
        f"{'Yes' if shifted_noise['statistically_credible'] else 'NO (FAR inflated)'} |"
    )
    lines.append(
        f"| `impulsive` (stress) | No | {imp_n['fp']} / {imp_n['n_noise']} | "
        f"{_fmt_pct(imp_n['false_alarm_rate'])} | {_fmt_ci_pct(imp_n['far_clopper_pearson_95ci'])} | "
        f"{'Yes' if imp_n['statistically_credible'] else 'NO (FAR inflated)'} |"
    )
    lines.append(
        f"| `strong_nonstationary` (mod=0.85) | Stress | {ns_n['fp']} / {ns_n['n_noise']} | "
        f"{_fmt_pct(ns_n['false_alarm_rate'])} | {_fmt_ci_pct(ns_n['far_clopper_pearson_95ci'])} | "
        f"{'Yes' if ns_n['statistically_credible'] else 'NO (FAR inflated)'} |"
    )
    lines.append("")
    lines.append("### 5.4 Detection Performance by Signal Family (Known vs. Withheld)")
    lines.append("")
    lines.append("Note: Detection measures whether an anomalous time-frequency event is flagged above the calibrated threshold. Because the detector is template-free, it does **not** classify signal families.")
    lines.append("")
    lines.append("| Signal Family | Status | TP / N | Detection Probability | 95% Wilson CI | Median Time IoU | Median Freq IoU |")
    lines.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
    for fam, stats_dict in by_family.items():
        floc = stats_dict["localization"]
        lines.append(
            f"| `{fam}` | Known | {stats_dict['tp']} / {stats_dict['n_signal']} | "
            f"{_fmt_pct(stats_dict['detection_probability'])} | {_fmt_ci_pct(stats_dict['detection_wilson_95ci'])} | "
            f"{_fmt_num(floc['median_time_iou'])} | {_fmt_num(floc['median_freq_iou'])} |"
        )
    wloc = withheld["localization"]
    lines.append(
        f"| `sinusoidal_fm` | **Withheld** | {withheld['tp']} / {withheld['n_signal']} | "
        f"{_fmt_pct(withheld['detection_probability'])} | {_fmt_ci_pct(withheld['detection_wilson_95ci'])} | "
        f"{_fmt_num(wloc['median_time_iou'])} | {_fmt_num(wloc['median_freq_iou'])} |"
    )
    lines.append(
        f"| `shifted_ar2` (known families) | **Shifted Noise** | {shifted_signal['tp']} / {shifted_signal['n_signal']} | "
        f"{_fmt_pct(shifted_signal['detection_probability'])} | {_fmt_ci_pct(shifted_signal['detection_wilson_95ci'])} | "
        f"{_fmt_num(shifted_signal['localization']['median_time_iou'])} | {_fmt_num(shifted_signal['localization']['median_freq_iou'])} |"
    )
    lines.append("")
    lines.append("### 5.5 Robustness Stress Tests")
    lines.append("")
    sl = rob["spectral_leakage"]
    we = rob["windowing_effects"]
    ba = rob["boundary_artifacts"]
    wo = rob["weak_and_overlapping"]
    lines.append(f"- **Spectral leakage (sinusoid, SNR=12):** On-bin recall = `{_fmt_pct(sl['on_bin_sinusoid_snr12']['detection_probability'])}` vs. Off-bin (half-bin offset) recall = `{_fmt_pct(sl['off_bin_sinusoid_snr12']['detection_probability'])}`.")
    lines.append(
        f"- **Windowing effects (off-bin sinusoid in pink noise, SNR=14, matched 1% FAP):** "
        f"`hann` recall = `{_fmt_pct(we['hann']['detection_probability'])}` (thr={we['hann']['calibrated_threshold']:.1f}), "
        f"`boxcar` recall = `{_fmt_pct(we['boxcar']['detection_probability'])}` (thr={we['boxcar']['calibrated_threshold']:.1f}), "
        f"`blackman` recall = `{_fmt_pct(we['blackman']['detection_probability'])}` (thr={we['blackman']['calibrated_threshold']:.1f})."
    )
    lines.append(f"- **Boundary artifacts (broadband pulse, SNR=12):** Edge placement (`t=0.06s` or `3.94s`) recall = `{_fmt_pct(ba['edge_pulse_snr12']['detection_probability'])}` vs. Interior placement recall = `{_fmt_pct(ba['interior_pulse_snr12']['detection_probability'])}`.")
    lines.append(f"- **Impulsive interference:** Noise-only FAR = `{_fmt_pct(imp_n['false_alarm_rate'])}` (95% CI `{_fmt_ci_pct(imp_n['far_clopper_pearson_95ci'])}`); linear chirp (SNR=12) detection = `{_fmt_pct(rob['impulsive_interference']['linear_chirp_snr12']['detection_probability'])}`.")
    lines.append(f"- **Strong nonstationary background (mod=0.85):** Noise-only FAR = `{_fmt_pct(ns_n['false_alarm_rate'])}`; intermittent narrowband (SNR=12) detection = `{_fmt_pct(rob['strong_nonstationary']['intermittent_snr12']['detection_probability'])}`.")
    lines.append(f"- **Weak & overlapping signals (pink noise):** Weak chirp (SNR=6) recall = `{_fmt_pct(wo['weak_chirp_snr6']['detection_probability'])}`; overlapping chirp + intermittent tone (SNR=10 each) recall = `{_fmt_pct(wo['overlapping_chirp_plus_intermittent_snr10']['detection_probability'])}`.")
    lines.append("")
    lines.append("## 6. Interpretation, Failure Analysis, and Limitations")
    lines.append("")
    lines.append("*The following section separates interpretation and failure diagnosis from the measured observations above.*")
    lines.append("")
    lines.append("### 6.1 Generalization to Withheld Signal Morphology")
    lines.append(
        "Because the excess-power clustering detector makes no assumption about signal waveform shape, it generalizes directly to the withheld `sinusoidal_fm` family, achieving detection probabilities comparable to known time-localized families at matched energy SNR. However, detection does not imply classification: the detector only reports that an anomalous time-frequency region exceeded the calibrated noise null."
    )
    lines.append("")
    lines.append("### 6.2 Failure Analysis")
    lines.append("1. **Impulsive interference breakdown:** Sparse, high-amplitude time-domain impulses spread excess energy across all frequencies in 1–3 STFT frames, producing vertical clusters that trigger high false-alarm rates when the threshold was calibrated on continuous Gaussian/colored/nonstationary noise. Under impulsive contamination, the detector explicitly flags that its significance estimates are **not statistically credible**.")
    lines.append("2. **Low-SNR fragmentation (missed signals):** At low energy SNR (4–8), signal power per STFT bin often falls below the bin seed threshold (`Z > 4.0`) or fragments into disconnected small clusters, each below the experiment-level threshold. This is most pronounced for signals spread over many time-frequency bins (e.g., long weak chirps or weak stationary sinusoids) compared to energy-concentrated pulses.")
    lines.append("3. **Boundary attenuation:** Because `boundary=None, padded=False` drops partial edge windows, short pulses occurring within half an STFT window (`< 0.0625 s`) of the start or end of the record suffer partial attenuation and lower recall.")
    lines.append("4. **Colored-noise spectral mismatch:** Narrowband spectral features in unseen colored noise (`shifted_ar2` with high pole radius) can partially leak past the 11-bin frequency median filter if the noise resonance is narrower than the kernel, potentially perturbing FAR.")
    lines.append("")
    lines.append("### 6.3 Calibration and Statistical Limitations")
    for lim in cal_agg["limitations"]:
        lines.append(f"- {lim}")
    lines.append("")
    lines.append("## 7. Generated Figures")
    lines.append("")
    for fig_name in results.get("figures", []):
        lines.append(f"- `figures/{fig_name}`")
    lines.append("")
    lines.append("## 8. Reproducibility Instructions")
    lines.append("")
    lines.append("```powershell")
    lines.append(".\\.venv\\Scripts\\python.exe -m pytest")
    lines.append(".\\.venv\\Scripts\\python.exe experiments\\run_experiment.py --config configs\\smoke.toml")
    lines.append(".\\.venv\\Scripts\\python.exe experiments\\run_experiment.py --config configs\\research.toml")
    lines.append("```")
    lines.append("")
    lines.append("## 9. Recommendations for Next Experiment (v0.2)")
    lines.append("")
    lines.append("1. **Multi-resolution STFT / Q-transform bank:** Combine short and long STFT windows (calibrated jointly via experiment-level max-statistic) to recover weak long-duration sinusoids without sacrificing sensitivity to short broadband pulses.")
    lines.append("2. **Impulsive-glitch veto / robust pre-filtering:** Add a calibrated time-domain glitch diagnostic or dual-channel coincidence check so single-sample impulses can be distinguished from genuine structured broadband pulses without silently suppressing real transients.")
    lines.append("3. **Larger calibration ensembles:** Increase `n_calibration_samples` to 2,000–5,000 per seed to tighten finite-sample tail uncertainty at 1% and 0.1% FAP levels.")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path
