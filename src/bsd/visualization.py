"""Publication-quality scientific visualization for Blind Signal Discovery (BSD).

All figures use English labels, physical/simulation units, clear captions, and
explicitly identify all plotted waveforms as simulated synthetic data.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from bsd.calibration import CalibrationResult
from bsd.detector import DetectionResult
from bsd.generator import GroundTruth, Observation

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 180,
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
)

_SIM_NOTE = "Simulated synthetic data — Blind Signal Discovery (BSD) v0.1"


def _add_footer(fig: plt.Figure, text: str = _SIM_NOTE) -> None:
    fig.text(0.5, 0.005, text, ha="center", va="bottom", fontsize=8, style="italic", color="#444444")


def plot_time_domain_examples(
    examples: list[tuple[str, Observation, GroundTruth]],
    output_path: Path,
) -> Path:
    """Plot representative simulated time-domain waveforms across families."""
    n = len(examples)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.5, 2.5 * nrows), sharex=True)
    axes_flat = np.atleast_1d(axes).ravel()

    for ax, (title, obs, gt) in zip(axes_flat, examples):
        t = np.arange(len(obs.signal)) / obs.sample_rate
        ax.plot(t, obs.signal, lw=0.6, color="#1f4e79")
        subtitle = f"{title} | noise={gt.noise_model}"
        if gt.is_signal and gt.target_snr is not None:
            subtitle += f" | energy SNR={gt.target_snr:.1f}"
        ax.set_title(subtitle)
        ax.set_ylabel("Amplitude (arb.)")
    for ax in axes_flat[n:]:
        ax.axis("off")
    for ax in axes_flat[-ncols:]:
        ax.set_xlabel("Time (s)")

    fig.suptitle("Simulated Time-Domain Examples (Noise-Only and Injected Signal Families)", y=0.995)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_spectrograms_and_excess_maps(
    examples: list[tuple[str, Observation, GroundTruth, DetectionResult]],
    output_path: Path,
) -> Path:
    """Plot raw STFT power spectrograms (dB) alongside normalized excess-power maps."""
    n = len(examples)
    fig, axes = plt.subplots(n, 2, figsize=(11.5, 2.7 * n), sharex=True, sharey=True)
    if n == 1:
        axes = np.array([axes])

    for row, (label, _obs, gt, det) in enumerate(examples):
        assert det.tf_maps is not None
        maps = det.tf_maps
        power_db = 10.0 * np.log10(np.maximum(maps.power_spectrogram, 1e-18))

        im0 = axes[row, 0].pcolormesh(
            maps.times, maps.freqs, power_db, shading="auto", cmap="viridis"
        )
        axes[row, 0].set_title(f"{label} — Raw STFT Spectrogram (dB) [{gt.noise_model}]")
        axes[row, 0].set_ylabel("Frequency (Hz)")
        fig.colorbar(im0, ax=axes[row, 0], label="Power (dB)", pad=0.015)

        im1 = axes[row, 1].pcolormesh(
            maps.times,
            maps.freqs,
            np.clip(maps.normalized_excess_map, 0.0, 20.0),
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=20.0,
        )
        axes[row, 1].set_title(
            f"{label} — Normalized Excess Power Z[f,t] (max score={det.max_anomaly_score:.1f})"
        )
        fig.colorbar(im1, ax=axes[row, 1], label="Excess Power Z", pad=0.015)

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    fig.suptitle("Raw Spectrograms vs. Whitened Excess-Power Anomaly Score Maps", y=0.995)
    fig.tight_layout(rect=[0, 0.02, 1, 0.975])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_candidate_overlays(
    examples: list[tuple[str, Observation, GroundTruth, DetectionResult]],
    output_path: Path,
) -> Path:
    """Plot excess-power maps with detected candidate boxes and ground-truth overlays."""
    n = len(examples)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.5, 3.3 * nrows), sharex=True, sharey=True)
    axes_flat = np.atleast_1d(axes).ravel()

    for ax, (label, _obs, gt, det) in zip(axes_flat, examples):
        assert det.tf_maps is not None
        maps = det.tf_maps
        im = ax.pcolormesh(
            maps.times,
            maps.freqs,
            np.clip(maps.normalized_excess_map, 0.0, 20.0),
            shading="auto",
            cmap="inferno",
            vmin=0.0,
            vmax=20.0,
        )
        fig.colorbar(im, ax=ax, label="Z[f,t]", pad=0.015)

        if gt.is_signal and gt.t_start is not None and gt.f_min is not None:
            gt_rect = patches.Rectangle(
                (gt.t_start, gt.f_min),
                max(0.05, (gt.t_end or gt.t_start) - gt.t_start),
                max(8.0, (gt.f_max or gt.f_min) - gt.f_min),
                linewidth=1.5,
                edgecolor="#ffffff",
                linestyle="--",
                facecolor="none",
                label="Ground truth (simulated)",
            )
            ax.add_patch(gt_rect)

        for idx, cand in enumerate(det.candidates[:3]):
            c_rect = patches.Rectangle(
                (cand.t_start, cand.f_min),
                max(0.05, cand.t_end - cand.t_start),
                max(8.0, cand.f_max - cand.f_min),
                linewidth=1.6,
                edgecolor="#00e5ff",
                linestyle="-",
                facecolor="none",
                label="Detected candidate" if idx == 0 else None,
            )
            ax.add_patch(c_rect)

        ax.set_title(
            f"{label} | SNR={gt.target_snr} | score={det.max_anomaly_score:.1f} | cand={len(det.candidates)}"
        )
        ax.set_ylabel("Frequency (Hz)")
        ax.legend(loc="upper right", framealpha=0.75)

    for ax in axes_flat[n:]:
        ax.axis("off")
    for ax in axes_flat[-ncols:]:
        ax.set_xlabel("Time (s)")

    fig.suptitle("Candidate Event Overlays on Excess-Power Maps (Cyan=Detected, Dashed White=Truth)", y=0.995)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_detection_vs_snr(
    family_snr_curves: dict[str, dict[float, dict[str, Any]]],
    target_fap: float,
    output_path: Path,
) -> Path:
    """Plot detection probability vs. energy SNR across known and withheld families."""
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    styles = {
        "sinusoid": ("-o", "#1f77b4", "Sinusoid (known)"),
        "broadband_pulse": ("-s", "#2ca02c", "Broadband pulse (known)"),
        "linear_chirp": ("-^", "#ff7f0e", "Linear chirp (known)"),
        "intermittent_narrowband": ("-d", "#9467bd", "Intermittent narrowband (known)"),
        "sinusoidal_fm": ("--*", "#d62728", "Sinusoidal FM (WITHHELD)"),
    }

    for family, snr_map in family_snr_curves.items():
        snrs = sorted(snr_map.keys())
        probs = [snr_map[s]["detection_probability"] for s in snrs]
        lo = [snr_map[s]["detection_probability"] - snr_map[s]["detection_wilson_95ci"][0] for s in snrs]
        hi = [snr_map[s]["detection_wilson_95ci"][1] - snr_map[s]["detection_probability"] for s in snrs]
        fmt, color, label = styles.get(family, ("-o", "#333333", family))
        ax.errorbar(
            snrs,
            probs,
            yerr=[lo, hi],
            fmt=fmt,
            color=color,
            capsize=3,
            lw=1.7,
            markersize=6,
            label=label,
        )

    ax.axhline(target_fap, color="#555555", linestyle=":", lw=1.2, label=f"Target FAP ({target_fap*100:.1f}%)")
    ax.set_xlabel("Injected Signal Energy SNR")
    ax.set_ylabel("Empirical Detection Probability")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("Blind Detection Probability vs. Signal Energy SNR (95% Wilson CIs)")
    ax.legend(loc="lower right", framealpha=0.9)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_false_alarm_results(
    far_by_condition: dict[str, dict[str, Any]],
    target_fap: float,
    output_path: Path,
) -> Path:
    """Plot empirical false-alarm rates with 95% exact Clopper-Pearson CIs."""
    labels = list(far_by_condition.keys())
    rates = [far_by_condition[k]["false_alarm_rate"] * 100.0 for k in labels]
    lo = [
        (far_by_condition[k]["false_alarm_rate"] - far_by_condition[k]["far_clopper_pearson_95ci"][0]) * 100.0
        for k in labels
    ]
    hi = [
        (far_by_condition[k]["far_clopper_pearson_95ci"][1] - far_by_condition[k]["false_alarm_rate"]) * 100.0
        for k in labels
    ]
    colors = [
        "#c53030" if far_by_condition[k]["far_clopper_pearson_95ci"][0] > target_fap else "#2b6cb0"
        for k in labels
    ]

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    x = np.arange(len(labels))
    bars = ax.bar(x, rates, color=colors, alpha=0.85, width=0.58)
    ax.errorbar(x, rates, yerr=[lo, hi], fmt="none", ecolor="#111111", capsize=4, lw=1.3)
    ax.axhline(target_fap * 100.0, color="#111111", linestyle="--", lw=1.4, label=f"Target FAP = {target_fap*100:.1f}%")

    for bar, rate, hi_err, k in zip(bars, rates, hi, labels):
        n_val = far_by_condition[k]["n_noise"]
        fp_val = far_by_condition[k]["fp"]
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + hi_err + 1.5,
            f"{rate:.2f}%\n({fp_val}/{n_val})",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Empirical False-Alarm Rate (%)")
    ax.set_title("Empirical False-Alarm Rate Across Calibrated and Shifted Noise Conditions (95% Exact CIs)")
    ax.set_ylim(0, max(max(r + h for r, h in zip(rates, hi)) + 10.0, 8.0))
    ax.legend(loc="upper left")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    _add_footer(fig, "Blue = within/compatible with calibrated FAP; Red = statistically significant FAR inflation. " + _SIM_NOTE)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_confusion_matrix(
    summary: dict[str, Any],
    output_path: Path,
) -> Path:
    """Plot experiment-level confusion matrix for held-out blind evaluation."""
    cm = np.array(
        [
            [summary["tn"], summary["fp"]],
            [summary["fn"], summary["tp"]],
        ],
        dtype=np.int64,
    )
    row_sums = np.maximum(cm.sum(axis=1, keepdims=True), 1)
    cm_norm = cm / row_sums

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    im = ax.imshow(cm_norm, cmap="Blues", vmin=0.0, vmax=1.0)
    fig.colorbar(im, ax=ax, label="Row-Normalized Proportion")

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Predicted: Noise", "Predicted: Anomaly"])
    ax.set_yticklabels(["True: Noise-Only", "True: Signal Injected"])
    ax.grid(False)

    cell_labels = [["TN", "FP"], ["FN", "TP"]]
    for i in range(2):
        for j in range(2):
            val = cm[i, j]
            pct = cm_norm[i, j] * 100.0
            color = "white" if cm_norm[i, j] > 0.55 else "#111111"
            ax.text(j, i, f"{cell_labels[i][j]}\n{val}\n({pct:.1f}%)", ha="center", va="center", color=color, fontsize=10, fontweight="bold")

    ax.set_title("Held-Out Blind Evaluation Confusion Matrix\n(In-Distribution Noise)", pad=10)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_failure_examples(
    examples: list[tuple[str, Observation, GroundTruth, DetectionResult]],
    title: str,
    output_path: Path,
) -> Path:
    """Plot time-domain waveform and excess-power map for FP or FN diagnostic cases."""
    if not examples:
        fig, ax = plt.subplots(figsize=(7.5, 3.0))
        ax.text(0.5, 0.5, "No failure examples occurred in this split.", ha="center", va="center", fontsize=11)
        ax.axis("off")
        fig.suptitle(title)
        fig.savefig(output_path)
        plt.close(fig)
        return output_path

    n = len(examples)
    fig, axes = plt.subplots(n, 2, figsize=(11.5, 2.8 * n))
    if n == 1:
        axes = np.array([axes])

    for row, (label, obs, gt, det) in enumerate(examples):
        assert det.tf_maps is not None
        t = np.arange(len(obs.signal)) / obs.sample_rate
        axes[row, 0].plot(t, obs.signal, lw=0.6, color="#2d3748")
        axes[row, 0].set_title(f"{label} — Time Series ({gt.noise_model})")
        axes[row, 0].set_ylabel("Amplitude (arb.)")

        maps = det.tf_maps
        im = axes[row, 1].pcolormesh(
            maps.times,
            maps.freqs,
            np.clip(maps.normalized_excess_map, 0.0, 20.0),
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=20.0,
        )
        fig.colorbar(im, ax=axes[row, 1], label="Z[f,t]", pad=0.015)

        if gt.is_signal and gt.t_start is not None and gt.f_min is not None:
            axes[row, 1].add_patch(
                patches.Rectangle(
                    (gt.t_start, gt.f_min),
                    max(0.05, (gt.t_end or gt.t_start) - gt.t_start),
                    max(8.0, (gt.f_max or gt.f_min) - gt.f_min),
                    linewidth=1.5,
                    edgecolor="#ffffff",
                    linestyle="--",
                    facecolor="none",
                )
            )
        for cand in det.candidates[:2]:
            axes[row, 1].add_patch(
                patches.Rectangle(
                    (cand.t_start, cand.f_min),
                    max(0.05, cand.t_end - cand.t_start),
                    max(8.0, cand.f_max - cand.f_min),
                    linewidth=1.6,
                    edgecolor="#00e5ff",
                    facecolor="none",
                )
            )

        thr_str = f"{det.threshold_used:.1f}" if det.threshold_used is not None else "N/A"
        axes[row, 1].set_title(
            f"Excess Power Map (score={det.max_anomaly_score:.1f}, thr={thr_str})"
        )
        axes[row, 1].set_ylabel("Frequency (Hz)")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    fig.suptitle(title, y=0.995)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def plot_calibration_distribution(
    calibration: CalibrationResult,
    output_path: Path,
) -> Path:
    """Plot empirical calibration null distribution of experiment-level max scores."""
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    scores = calibration.sorted_null_scores
    ax.hist(scores, bins=40, color="#4a5568", alpha=0.8, edgecolor="white", density=True, label="Calibration noise-only max scores")
    ax.axvline(
        calibration.threshold,
        color="#d62728",
        lw=2.0,
        linestyle="-",
        label=f"Calibrated 1% FAP threshold = {calibration.threshold:.2f}",
    )
    ax.axvline(
        calibration.naive_per_bin_exp_threshold,
        color="#2b6cb0",
        lw=1.5,
        linestyle="--",
        label=f"Naive single-bin Exp(1) 1% threshold = {calibration.naive_per_bin_exp_threshold:.2f}",
    )
    ax.set_xlabel("Experiment-Level Maximum Cluster Anomaly Score")
    ax.set_ylabel("Probability Density")
    ax.set_title(
        f"Noise-Only Calibration Null Distribution (N={calibration.n_calibration_samples}, "
        f"{calibration.total_tf_bins_per_sample} TF bins/sample)"
    )
    ax.legend(loc="upper right")
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    _add_footer(fig)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path
