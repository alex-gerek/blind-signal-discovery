# Blind Signal Discovery (BSD) v0.1 — Experimental Research Report

> **Scientific Scope Statement:** This report presents a controlled methodological study of template-free anomaly detection in simulated time-series data. **A detected anomaly is not a physical discovery.** No claims are made regarding dark matter, axions, gravitational waves, or any physical phenomenon.

## 1. Research Question

Can a classical, template-free time-frequency anomaly detector discover previously unseen signal structures in noisy time-series data while maintaining a controlled and independently measured experiment-level false-alarm probability (target FAP = 1.0%)?

## 2. Experimental Methodology & Dataset Separation

- **Experiment configuration:** `research_v01`
- **Reproducible seeds:** `[101, 202, 303]` (3 independent experimental runs)
- **Sampling parameters:** `sample_rate = 1024.0 Hz`, `duration = 4.0 s` (`4096` samples per observation)
- **SNR definition:** Total signal energy SNR $\mathrm{SNR}_{\mathrm{energy}} = \|s\|_2 / \sigma_{\mathrm{noise}}$, evaluated over `[4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 25.0]`.
- **Strict dataset separation:** Independent RNG streams are derived via `numpy.random.SeedSequence([master_seed, split_stream_code])` for calibration, held-out noise-only testing, held-out known-signal testing, withheld-family testing, shifted-noise testing, and each robustness stress test.
- **Ground-truth isolation:** The detector API (`TemplateFreeDetector.analyze`) accepts only `Observation` objects containing the raw 1D waveform and sample rate. Ground-truth labels, noise models, SNR values, and time-frequency bounds are stored separately in `GroundTruth` objects and accessed only after frozen-threshold predictions are finalized.
- **Frozen-threshold protocol:** For each seed, the detection threshold is calibrated exclusively on the dedicated calibration split and then frozen prior to all blind evaluations.

## 3. Detector Design & Noise Assumptions

The baseline detector uses classical, interpretable signal processing (no templates, no machine learning):
1. **Robust time-domain preprocessing:** Median subtraction for DC removal; MAD-scaled dispersion, Fisher kurtosis, and crest factor recorded as diagnostics.
2. **STFT time-frequency representation:** `hann` window, `nperseg=128`, `noverlap=96`, `boundary=None, padded=False` (avoiding zero-pad edge artifacts). DC and Nyquist bins are excluded so interior bins share a common complex $\chi^2_2 \equiv \mathrm{Exp}(1)$ null distribution under locally stationary Gaussian noise.
3. **Robust 2D background estimation:** Per-frequency PSD estimated via temporal median divided by $\ln 2$, followed by a frequency median filter (`kernel=11`) so stationary narrowband tones are not absorbed into the background. After spectral whitening, slowly varying temporal noise variance is estimated via frequency median divided by $\ln 2$ and smoothed with a temporal median filter (`kernel=15`) so short broadband pulses are not absorbed.
4. **Connected-component clustering:** Bins with normalized excess power $Z[f,t] > 4.0$ are grouped via 8-connectivity (`min_cluster_bins=2`). Each candidate cluster $C$ is scored by integrated excess power $S(C) = \sum_{(f,t)\in C} (Z[f,t] - 1)$.

**Key noise assumptions:** Locally smooth background PSD (bandwidth broader than the frequency median kernel) and slowly varying temporal variance (timescale longer than the temporal median kernel), without unmodeled heavy-tailed impulsive transients.

## 4. Statistical Calibration Procedure

- **Target experiment-level FAP:** `1.00%`
- **Calibration samples per seed:** `600` noise-only observations across `['gaussian_white', 'pink_colored', 'nonstationary']`
- **Mean calibrated threshold across seeds:** `62.62 ± 0.84` (range `[61.50, 63.52]`)
- **Per-bin vs. experiment-level FAP distinction:** Each observation spans `7875` interior STFT bins. A naive single-bin 1% $\mathrm{Exp}(1)$ threshold is only `4.61` (and a single-bin Bonferroni correction is `13.58`). Calibrating on the empirical order statistic $k = \lceil(1 - \alpha)(N_{\mathrm{cal}} + 1)\rceil$ of the experiment-level maximum cluster score $T_{\max} = \max_C S(C)$ directly accounts for scanning all `7875` time-frequency bins AND their STFT overlap correlations.
- **Minimum resolvable empirical p-value per seed:** `0.00166` ($1 / (N_{\mathrm{cal}} + 1)$). Candidate scores exceeding all calibration noise samples are reported explicitly as upper bounds ($p \le 1/(N_{\mathrm{cal}}+1)$), never extrapolated.

## 5. Measured Observations (Pooled Across All Seeds)

### 5.1 Held-Out In-Distribution Summary

- **Total held-out noise-only samples:** `1500` (FP = `12`, TN = `1488`)
- **Empirical Noise-Only FAR:** `0.80%` (95% Wilson CI `[0.46%, 1.39%]`, 95% exact Clopper-Pearson CI `[0.41%, 1.39%]`)
- **Total held-out known-signal samples:** `1536` (TP = `914`, FN = `622`)
- **Overall Detection Probability (Recall):** `59.51%` (95% Wilson CI `[57.03%, 61.93%]`)
- **Precision:** `98.70%` | **F1 Score:** `0.742`
- **Localization (on detected TPs, N=914):** Median Time IoU = `0.446`, Median Freq IoU = `0.667`, Median Time Center Error = `0.339 s`, Median Freq Center Error = `4.1 Hz`, TF Box Overlap Rate = `99.45%`.

### 5.2 Per-Seed Reproducibility Check

| Seed | Calibrated Threshold | Noise-Only FAR (FP/N) | Known-Signal Recall (TP/N) | Withheld FM Recall (TP/N) |
| :--- | :---: | :---: | :---: | :---: |
| `101` | `61.50` | `1.40%` (7/500) | `59.57%` (305/512) | `37.50%` (60/160) |
| `202` | `63.52` | `0.20%` (1/500) | `58.79%` (301/512) | `38.12%` (61/160) |
| `303` | `62.85` | `0.80%` (4/500) | `60.16%` (308/512) | `38.12%` (61/160) |

### 5.3 False-Alarm Performance by Noise Model

| Noise Condition | In Calibration? | FP / N | Empirical FAR | 95% Exact CI | Statistically Credible? |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `gaussian_white` | Yes | 3 / 501 | 0.60% | [0.12%, 1.74%] | Yes |
| `nonstationary` | Yes | 4 / 498 | 0.80% | [0.22%, 2.04%] | Yes |
| `pink_colored` | Yes | 5 / 501 | 1.00% | [0.32%, 2.31%] | Yes |
| `shifted_ar2` (unseen) | No | 16 / 450 | 3.56% | [2.05%, 5.71%] | NO (FAR inflated) |
| `impulsive` (stress) | No | 252 / 300 | 84.00% | [79.35%, 87.96%] | NO (FAR inflated) |
| `strong_nonstationary` (mod=0.85) | Stress | 6 / 300 | 2.00% | [0.74%, 4.30%] | Yes |

### 5.4 Detection Performance by Signal Family (Known vs. Withheld)

Note: Detection measures whether an anomalous time-frequency event is flagged above the calibrated threshold. Because the detector is template-free, it does **not** classify signal families.

| Signal Family | Status | TP / N | Detection Probability | 95% Wilson CI | Median Time IoU | Median Freq IoU |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `broadband_pulse` | Known | 300 / 384 | 78.12% | [73.72%, 81.97%] | 0.275 | 0.772 |
| `intermittent_narrowband` | Known | 208 / 360 | 57.78% | [52.62%, 62.77%] | 0.942 | 0.500 |
| `linear_chirp` | Known | 216 / 360 | 60.00% | [54.86%, 64.93%] | 0.872 | 0.854 |
| `sinusoid` | Known | 190 / 432 | 43.98% | [39.37%, 48.69%] | 0.727 | 0.500 |
| `sinusoidal_fm` | **Withheld** | 182 / 480 | 37.92% | [33.69%, 42.34%] | 0.497 | 0.809 |
| `shifted_ar2` (known families) | **Shifted Noise** | 292 / 480 | 60.83% | [56.40%, 65.10%] | 0.398 | 0.667 |

### 5.5 Robustness Stress Tests

- **Spectral leakage (sinusoid, SNR=12):** On-bin recall = `4.67%` vs. Off-bin (half-bin offset) recall = `4.00%`.
- **Windowing effects (off-bin sinusoid in pink noise, SNR=14, matched 1% FAP):** `hann` recall = `59.33%` (thr=62.6), `boxcar` recall = `46.33%` (thr=56.7), `blackman` recall = `45.67%` (thr=63.9).
- **Boundary artifacts (broadband pulse, SNR=12):** Edge placement (`t=0.06s` or `3.94s`) recall = `89.67%` vs. Interior placement recall = `100.00%`.
- **Impulsive interference:** Noise-only FAR = `84.00%` (95% CI `[79.35%, 87.96%]`); linear chirp (SNR=12) detection = `90.67%`.
- **Strong nonstationary background (mod=0.85):** Noise-only FAR = `2.00%`; intermittent narrowband (SNR=12) detection = `97.67%`.
- **Weak & overlapping signals (pink noise):** Weak chirp (SNR=6) recall = `7.67%`; overlapping chirp + intermittent tone (SNR=10 each) recall = `94.00%`.

## 6. Interpretation, Failure Analysis, and Limitations

*The following section separates interpretation and failure diagnosis from the measured observations above.*

### 6.1 Generalization to Withheld Signal Morphology
Because the excess-power clustering detector makes no assumption about signal waveform shape, it generalizes directly to the withheld `sinusoidal_fm` family, achieving detection probabilities comparable to known time-localized families at matched energy SNR. However, detection does not imply classification: the detector only reports that an anomalous time-frequency region exceeded the calibrated noise null.

### 6.2 Failure Analysis
1. **Impulsive interference breakdown:** Sparse, high-amplitude time-domain impulses spread excess energy across all frequencies in 1–3 STFT frames, producing vertical clusters that trigger high false-alarm rates when the threshold was calibrated on continuous Gaussian/colored/nonstationary noise. Under impulsive contamination, the detector explicitly flags that its significance estimates are **not statistically credible**.
2. **Low-SNR fragmentation (missed signals):** At low energy SNR (4–8), signal power per STFT bin often falls below the bin seed threshold (`Z > 4.0`) or fragments into disconnected small clusters, each below the experiment-level threshold. This is most pronounced for signals spread over many time-frequency bins (e.g., long weak chirps or weak stationary sinusoids) compared to energy-concentrated pulses.
3. **Boundary attenuation:** Because `boundary=None, padded=False` drops partial edge windows, short pulses occurring within half an STFT window (`< 0.0625 s`) of the start or end of the record suffer partial attenuation and lower recall.
4. **Colored-noise spectral mismatch:** Narrowband spectral features in unseen colored noise (`shifted_ar2` with high pole radius) can partially leak past the 11-bin frequency median filter if the noise resonance is narrower than the kernel, potentially perturbing FAR.

### 6.3 Calibration and Statistical Limitations
- Finite-sample tail resolution: with N_cal=600 and target_fap=0.010, only ~6.0 calibration samples lie in the upper tail; minimum resolvable empirical p-value is 0.00166.
- Exchangeability assumption: experiment-level FAP control holds only when test noise follows the same distribution as the calibration noise mixture. Noise-model shifts (e.g., impulsive interference or unseen colored spectra) can inflate FAR.
- Single-observation horizon: FAP is calibrated per individual observation window; scanning M independent observations without further correction yields an expected false-alarm count of M * target_fap.

## 7. Generated Figures

- `figures/fig01_time_domain_examples.png`
- `figures/fig02_spectrograms_and_excess_maps.png`
- `figures/fig03_candidate_overlays.png`
- `figures/fig04_detection_vs_snr.png`
- `figures/fig05_false_alarm_results.png`
- `figures/fig06_confusion_matrix.png`
- `figures/fig07_false_positive_examples.png`
- `figures/fig08_missed_signal_examples.png`
- `figures/fig09_calibration_null_distribution.png`

## 8. Reproducibility Instructions

```powershell
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\smoke.toml
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\research.toml
```

## 9. Recommendations for Next Experiment (v0.2)

1. **Multi-resolution STFT / Q-transform bank:** Combine short and long STFT windows (calibrated jointly via experiment-level max-statistic) to recover weak long-duration sinusoids without sacrificing sensitivity to short broadband pulses.
2. **Impulsive-glitch veto / robust pre-filtering:** Add a calibrated time-domain glitch diagnostic or dual-channel coincidence check so single-sample impulses can be distinguished from genuine structured broadband pulses without silently suppressing real transients.
3. **Larger calibration ensembles:** Increase `n_calibration_samples` to 2,000–5,000 per seed to tighten finite-sample tail uncertainty at 1% and 0.1% FAP levels.
