# Blind Signal Discovery (BSD) v0.1 — Scientific Review Package

> **Review Integrity Notice:** No code, detector parameters, thresholds, or test results in this repository were modified during this audit. All metrics below reflect the frozen v0.1 benchmark (`research_v01`).

---

## 1. Full Research Report (`research_v01`)

The full generated report is stored verbatim at [results/research_v01/report.md](../results/research_v01/report.md). Its measured results across 3 independent seeds (`[101, 202, 303]`) are summarized below together with the complete Family x SNR detection table from [results/research_v01/summary_by_family_snr.csv](../results/research_v01/summary_by_family_snr.csv).

### 1.1 Headline Measured Results (Pooled Across 3 Seeds)

| Metric / Split | Value | 95% CI | Notes |
| :--- | :---: | :---: | :--- |
| Calibration samples per seed ($N_{\text{cal}}$) | 600 | — | Mixture: `gaussian_white`, `pink_colored`, `nonstationary` |
| Calibrated 1% FAP threshold ($\hat{\tau}$) | `62.62 ± 0.84` | range `[61.50, 63.52]` | Per-seed order statistic $k=595$ of 600 |
| Held-out in-distribution noise FAR (`12 / 1500`) | **0.80%** | Exact `[0.41%, 1.39%]` | Credible (consistent with 1.0% target) |
| Held-out `gaussian_white` FAR (`3 / 501`) | 0.60% | Exact `[0.12%, 1.74%]` | Credible |
| Held-out `nonstationary` FAR (`4 / 498`) | 0.80% | Exact `[0.22%, 2.04%]` | Credible |
| Held-out `pink_colored` FAR (`5 / 501`) | 1.00% | Exact `[0.32%, 2.31%]` | Credible |
| Held-out known-signal recall (`914 / 1536`) | **59.51%** | Wilson `[57.03%, 61.93%]` | Pooled across SNR 4–25; Precision `98.70%`, F1 `0.742` |
| Held-out **withheld** `sinusoidal_fm` recall (`182 / 480`) | **37.92%** | Wilson `[33.69%, 42.34%]` | Pooled across SNR 4–25; reaches `94.44%` at SNR 25 |
| Held-out unseen `shifted_ar2` noise FAR (`16 / 450`) | **3.56%** | Exact `[2.05%, 5.71%]` | **Flagged NOT credible** (FAR inflation) |
| Stress `strong_nonstationary` (mod=0.85) FAR (`6 / 300`) | 2.00% | Exact `[0.74%, 4.30%]` | Compatible with 1% within CI |
| Stress `impulsive` noise FAR (`252 / 300`) | **84.00%** | Exact `[79.35%, 87.96%]` | **Flagged NOT credible** (severe breakdown) |

### 1.2 Detection Probability by Family and Energy SNR

| Energy SNR | `sinusoid` (known) | `broadband_pulse` (known) | `linear_chirp` (known) | `intermittent_narrowband` (known) | `sinusoidal_fm` (**WITHHELD**) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 4.0 | 0.0% (0/54) | 14.8% (8/54) | 4.4% (2/45) | 2.2% (1/45) | 3.2% (2/63) |
| 6.0 | 7.4% (4/54) | 53.7% (29/54) | 15.6% (7/45) | 8.9% (4/45) | 4.8% (3/63) |
| 8.0 | 20.4% (11/54) | 78.4% (40/51) | 26.7% (12/45) | 33.3% (15/45) | 9.5% (6/63) |
| 10.0 | 33.3% (18/54) | 95.6% (43/45) | 60.0% (27/45) | 44.4% (20/45) | 22.2% (14/63) |
| 12.0 | 38.9% (21/54) | 100.0% (45/45) | 75.6% (34/45) | 82.2% (37/45) | 33.3% (21/63) |
| 16.0 | 59.3% (32/54) | 100.0% (45/45) | 97.8% (44/45) | 93.3% (42/45) | 63.2% (36/57) |
| 20.0 | 92.6% (50/54) | 100.0% (45/45) | 100.0% (45/45) | 100.0% (45/45) | 90.7% (49/54) |
| 25.0 | 100.0% (54/54) | 100.0% (45/45) | 100.0% (45/45) | 97.8% (44/45) | 94.4% (51/54) |

### 1.3 Localization & Robustness Summary

- **Localization on detected TPs (N=914):** TF box overlap rate `99.45%`, median time IoU `0.446`, median frequency IoU `0.667`, median time center error `0.339 s`, median frequency center error `4.1 Hz`.
- **Spectral leakage (sinusoid, SNR=12, white noise):** On-bin recall `4.67%` vs. half-bin off-bin recall `4.00%`.
- **Windowing (off-bin sinusoid in pink noise, SNR=14, matched 1% FAP):** `hann` `59.33%` (thr `62.6`), `boxcar` `46.33%` (thr `56.7`), `blackman` `45.67%` (thr `63.9`).
- **Boundary artifacts (broadband pulse, SNR=12):** Edge (`0.06 s` / `3.94 s`) recall `89.67%` vs. interior `100.00%`.
- **Weak & overlapping (pink noise):** Weak chirp (SNR=6) `7.67%`; overlapping chirp + intermittent tone (SNR=10 each) `94.00%`.

---

## 2. Detector and Calibration Implementation Details

### 2.1 Detector Pipeline ([src/bsd/detector.py](../src/bsd/detector.py))

1. **Robust time-domain preprocessing (`L143-L148`):**
   - DC offset removed via median subtraction: $x = x_{\text{raw}} - \mathrm{median}(x_{\text{raw}})$.
   - Time-domain diagnostics computed (not used to alter $x$): $\hat{\sigma}_{\text{MAD}} = 1.4826 \cdot \mathrm{median}(|x|)$, Fisher kurtosis, and crest factor $\max|x| / \hat{\sigma}_{\text{MAD}}$.
2. **STFT power spectrogram (`L151-L163`):**
   - `scipy.signal.stft` with `window="hann"`, `nperseg=128` (`125 ms`, $\Delta f = 8\text{ Hz}$), `noverlap=96` (75% overlap, $\Delta t = 31.25\text{ ms}$), `boundary=None, padded=False`.
   - For `4096` samples (`4.0 s` at `1024 Hz`), this produces `65` frequency bins x `125` time frames.
   - DC (`f = 0 Hz`) and Nyquist (`f = 512 Hz`) bins are excluded (`stft_mat[1:-1, :]`), leaving `63 x 125 = 7,875` interior complex bins. Under locally stationary Gaussian noise, interior $|X[f,t]|^2$ follows $\mathrm{Exp}(\lambda)$ ($\chi^2_2$), whereas DC/Nyquist are real-valued ($\chi^2_1$).
3. **2D robust background whitening (`L165-L180`):**
   - **Spectral baseline:** `raw_psd = median_t(P) / ln(2)` (unbiased for $\mathrm{Exp}(\lambda)$ since $\mathrm{median}(\mathrm{Exp}(\lambda)) = \lambda \ln 2$), smoothed across frequency with a 1D median filter (`freq_median_kernel = 11` bins = `88 Hz`, `mode="nearest"`). This prevents stationary narrowband tones (1–3 bins wide) from being absorbed into the background PSD.
   - **Spectral whitening:** $W[f,t] = P[f,t] / \max(S_{\text{bg}}(f), 10^{-18})$.
   - **Temporal variance baseline:** `raw_time_env = median_f(W) / ln(2)`, smoothed across time with a 1D median filter (`time_median_kernel = 15` frames $\approx 0.47\text{ s}$, `mode="nearest"`). This tracks slow nonstationary amplitude modulation without absorbing short broadband pulses (1–3 frames wide).
   - **Normalized excess-power map:** $Z[f,t] = W[f,t] / \max(T_{\text{bg}}(t), 10^{-18})$.
4. **Connected-component clustering & scoring (`L182-L263`):**
   - Binary seed mask: $M[f,t] = \mathbb{I}(Z[f,t] > \theta_{\text{bin}})$ with $\theta_{\text{bin}} = 4.0$ (single-bin exceedance probability $e^{-4} \approx 1.83\%$).
   - 8-connected-component labeling (`scipy.ndimage.label` with `3x3` ones kernel); clusters with fewer than `min_cluster_bins = 2` bins are discarded.
   - Cluster anomaly score (integrated excess power):
     $$S(C) = \sum_{(f,t) \in C} \bigl(Z[f,t] - 1.0\bigr)$$
   - Experiment-level test statistic: $T_{\max} = \max_C S(C)$ (or `0.0` if no cluster has $\ge 2$ bins).
   - Acceptance rule when frozen `CalibrationResult` is supplied (`L267`): strict inequality `c.anomaly_score > calibration.threshold`.

### 2.2 Calibration Implementation ([src/bsd/calibration.py](../src/bsd/calibration.py))

- `calibrate_detector` (`L103-L176`) runs `detector.analyze(obs, calibration=None)` on every observation in the dedicated calibration split, collects $T_1, \dots, T_{N_{\text{cal}}}$, sorts them as $T_{(1)} \le \dots \le T_{(N_{\text{cal}})}$, and sets:
  $$k = \min\!\left(N_{\text{cal}},\ \left\lceil (1 - \alpha)(N_{\text{cal}} + 1) \right\rceil\right), \qquad \hat{\tau} = T_{(k)}.$$
- Empirical p-value for a candidate with score $s$ (`L66-L78`):
  $$p_{\text{emp}}(s) = \frac{1 + \#\{i : T_i \ge s\}}{N_{\text{cal}} + 1},$$
  with `p_value_is_upper_bound = True` whenever $s > T_{(N_{\text{cal}})}$ (never extrapolating parametric tails).

---

## 3. Exact Data Separation and Ground-Truth Isolation Checks

### 3.1 Statistical Independence of Splits ([src/bsd/generator.py](../src/bsd/generator.py#L29-L77))

Every dataset split receives an independent RNG via `numpy.random.SeedSequence([master_seed, split_stream_code])` with disjoint integer codes:
- `development: 1001`, `calibration: 2002`, `test_noise: 3003`, `test_signal: 4004`, `test_withheld: 5005`, `test_shifted_noise: 6006`, `test_shifted_signal: 7007`, `robustness_*: 8001–8006`.
- `SeedSequence` spawns cryptographically-mixed independent PCG64 streams, so calibration, held-out tests, and robustness splits never share RNG state.

### 3.2 Structural Ground-Truth Isolation

1. **Data-type boundary ([src/bsd/generator.py](../src/bsd/generator.py#L46-L73)):**
   - `Observation` is a `@dataclass(frozen=True)` with **only 4 fields**: `sample_id`, `signal`, `sample_rate`, `duration`.
   - `GroundTruth` is a completely separate dataclass holding `is_signal`, `signal_family`, `noise_model`, `target_snr`, `empirical_energy_snr`, `rms_snr`, `peak_snr`, `t_start`, `t_end`, `f_min`, `f_max`, and `params`.
2. **Runtime enforcement ([src/bsd/detector.py](../src/bsd/detector.py#L134-L137)):**
   - `TemplateFreeDetector.analyze` explicitly checks `if not isinstance(observation, Observation): raise TypeError(...)`.
3. **Evaluation call site ([src/bsd/evaluation.py](../src/bsd/evaluation.py#L103-L104)):**
   - `detector.analyze(obs, calibration=calibration)` receives only `obs` and the frozen `calibration`; `gt` is consulted only afterward to assign `TP/FP/TN/FN` and compute localization errors.
4. **Automated unit test verification ([tests/test_bsd.py](../tests/test_bsd.py#L112-L127)):**
   - `test_detector_input_isolation` programmatically asserts that `Observation` has zero forbidden ground-truth fields and that passing `GroundTruth`, `(obs, gt)`, or raw `np.ndarray` to `detector.analyze` raises `TypeError`.

---

## 4. Full False-Alarm Calibration & Multiple-Comparison Methodology

### 4.1 Per-Bin Threshold vs. Experiment-Level FAP

Each `4.0 s` observation contains **$M = 7{,}875$ interior time-frequency bins**, strongly correlated along time (75% STFT overlap) and frequency (Hann main-lobe width $\approx 3$ bins):
- **Naive single-bin 1% $\text{Exp}(1)$ threshold:** $-\ln(0.01) = 4.61$. Applying this across 7,875 bins would trigger false alarms on virtually 100% of noise-only observations.
- **Single-bin Bonferroni threshold:** $-\ln(0.01 / 7875) = 13.58$. While stricter, it only bounds single-pixel exceedances and does not apply to 2D connected-component cluster sums $S(C)$, whose null distribution depends on both the look-elsewhere factor $M$ and the 2D STFT correlation structure.

### 4.2 Distribution-Free Order-Statistic Control

Instead of assuming independence across bins or fitting a parametric extreme-value tail, BSD calibrates directly on the **experiment-level maximum** $T_i = \max_C S(C)$ per noise-only observation. By exchangeability between calibration samples $T_1, \dots, T_{N_{\text{cal}}}$ and an independent test sample $T_{\text{new}}$ from the same noise distribution:
$$\Pr\bigl(T_{\text{new}} > T_{(k)}\bigr) \le \frac{N_{\text{cal}} + 1 - k}{N_{\text{cal}} + 1} \le \alpha.$$
For $N_{\text{cal}} = 600$ and $\alpha = 0.01$:
$$k = \lceil 0.99 \times 601 \rceil = 595 \implies \Pr(T_{\text{new}} > T_{(595)}) = \frac{6}{601} = 0.9983\% \le 1.0\%.$$
This automatically accounts for:
1. All `7,875` time-frequency bins searched per observation (the look-elsewhere effect).
2. All 2D STFT window-overlap correlations and cluster-formation statistics.
3. Finite-sample background-estimation variance from the 2D median filters.

**Empirical verification on held-out in-distribution test data (`1,500` independent observations across 3 seeds):** `12 / 1500 = 0.80%` (95% exact Clopper-Pearson CI `[0.41%, 1.39%]`).

---

## 5. Explanation of the 84% False-Alarm Rate Under Impulsive Interference

The `84.00%` (`252 / 300`, 95% exact CI `[79.35%, 87.96%]`) false-alarm rate on `impulsive` noise is a **direct, mathematically predictable consequence** of two interacting design choices:

1. **Violation of calibration exchangeability:** `impulsive` noise was deliberately excluded from `calibration_models = ["gaussian_white", "pink_colored", "nonstationary"]` in [configs/research.toml](../configs/research.toml#L16) to serve as an out-of-distribution stress test.
2. **Time-frequency physics of a short impulse in STFT:**
   - Each `impulsive` record injects `1` to `4` sparse spikes (`1`–`3` samples wide, amplitude $6\sigma$–$12\sigma$, [src/bsd/noise.py](../src/bsd/noise.py#L112-L132)).
   - A 1–3 sample Dirac-like spike has a nearly flat, **broadband Fourier spectrum** spanning all `63` interior frequency bins across `1`–`4` overlapping STFT frames (duration `125 ms`).
   - Because the impulse lasts $\le 4$ frames, the `15`-frame temporal median filter `bg_time` **completely ignores it** (just as designed for genuine broadband pulses).
   - Consequently, a single $8\sigma$–$12\sigma$ spike elevates $Z[f,t] > 4.0$ across a vertical column of ~20–60 adjacent frequency bins in 1–3 frames, easily accumulating an integrated excess power $S(C) \approx 80\text{–}300 \gg \hat{\tau} \approx 62.6$.
   - Because the detector is **template-free and single-channel**, a 2-sample impulsive glitch is structurally indistinguishable in the STFT magnitude plane from a genuine short broadband pulse without an explicit time-domain kurtosis/glitch veto or multi-channel coincidence check.
   - Crucially, `summarize_records` automatically detects that the 95% exact CI lower bound (`79.35%`) exceeds `target_fap = 1.0%` and flags `statistically_credible = False`.

---

## 6. Reproducibility Audit

### 6.1 Environment & Baseline Commit

- **Baseline benchmark commit:** `ef8d3ee` (`main`)
- **Runtime versions:** `Python 3.14.3`, `NumPy 2.5.5`, `SciPy 1.18.1`, `Matplotlib 3.11.2`, `pytest 9.2.0`.
- **Execution:** Local CPU only; `research.toml` runtime `74.47 s`, `smoke.toml` runtime `18.25 s`, `pytest` runtime `9.08 s`.

### 6.2 Exact Reproduction Commands (from repository root)

```powershell
# 1. Automated test suite (19 tests)
.\.venv\Scripts\python.exe -m pytest -v

# 2. Fast smoke-test configuration
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\smoke.toml

# 3. Full 3-seed research configuration
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\research.toml
```

### 6.3 Exact Configuration Files

- Full research config: [configs/research.toml](../configs/research.toml) (`seeds = [101, 202, 303]`, `sample_rate = 1024.0`, `duration = 4.0`, `snr_levels = [4, 6, 8, 10, 12, 16, 20, 25]`, `n_calibration_samples = 600`, `n_noise_test = 500`, `n_signal_test = 512`, `n_withheld_test = 160`, `n_shifted_noise = 150`, `n_shifted_signal = 160`, `n_robustness_per_condition = 100`).
- Fast smoke config: [configs/smoke.toml](../configs/smoke.toml) (`seeds = [101, 202]`, `n_calibration_samples = 120`, `n_noise_test = 100`, `n_signal_test = 120`).

---

## 7. Methodological Weaknesses and Potential Statistical Errors

An honest pre-publication audit identifies the following methodological weaknesses and statistical caveats in BSD v0.1:

1. **Global energy SNR vs. whitened matched-filter SNR in colored noise:**
   - Signal amplitude is scaled to global $\mathrm{SNR}_{\text{energy}} = \|s\|_2 / \sigma_{\text{noise}}$ ([src/bsd/generator.py](../src/bsd/generator.py#L79-L87)). In `pink_colored` or `shifted_ar2` noise, the local PSD is much higher at low frequencies than at high frequencies. Thus, two narrowband signals with identical global energy SNR have substantially different *whitened* SNRs depending on carrier frequency. This explains why `intermittent_narrowband` had `44/45` (`97.8%`) rather than `100%` at SNR=25.
2. **Non-divisible stratification counts in `research.toml`:**
   - `test_signal` cycles through `4 families x 8 SNRs x 3 noise_models = 96` combinations, but `n_signal_test = 512` is not divisible by `96` (`512 = 5x96 + 32`). Consequently, `sinusoid` (`432` pooled samples) and low-SNR `broadband_pulse` (`54` samples at SNR 4–6 vs. `45` at SNR 10–25) are slightly over-represented in the pooled overall recall (`59.51%`). Similarly, `n_withheld_test = 160` is not divisible by `24` (`160 = 6x24 + 16`), slightly overweighting lower SNRs (`63` samples vs. `54`) in the pooled withheld recall (`37.92%`). Per-family x per-SNR curves in `summary_by_family_snr.csv` are unaffected by this imbalance.
3. **Conditional vs. unconditional confidence intervals when pooling across seeds:**
   - Within each seed $s$, all 500 held-out noise test samples share a single finite-sample estimated threshold $\hat{\tau}_s$ (from $N_{\text{cal}} = 600$, where only ~6 samples populate the 1% upper tail). Therefore, test false alarms within seed $s$ are conditionally independent given $\hat{\tau}_s$, not unconditionally i.i.d. Standard Clopper-Pearson/Wilson CIs on the pooled `12 / 1500` count measure conditional FAR precision given the realized thresholds, slightly understating total uncertainty from finite $N_{\text{cal}}$ (visible in per-seed FAR spread: `1.40%`, `0.20%`, `0.80%`).
4. **Mixture-level vs. per-model calibration:**
   - Threshold $\hat{\tau}$ is calibrated on an equal mixture of `gaussian_white`, `pink_colored`, and `nonstationary` noise. While empirical per-model FARs happened to be well balanced (`0.60%`, `0.80%`, `1.00%`), calibrating on a mixture only guarantees average FAP $\le 1\%$ across the mixture, not strictly $\le 1\%$ on every component if one noise model has a heavier tail.
5. **Shifted AR(2) resonance leakage past the 11-bin frequency median filter:**
   - High pole radii (`0.72–0.86`) in `shifted_ar2` produce spectral peaks narrow enough to be slightly underestimated by the `11`-bin (`88 Hz`) frequency median filter, raising `Z[f,t]` near the resonant frequency and inflating FAR to `3.56%`.
6. **Single-resolution STFT bias against weak long-duration signals:**
   - Using a single `125 ms` window concentrates broadband pulse energy into 1–2 frames (`78.1%` pooled recall) while spreading a `4.0 s` stationary sinusoid across `125` frames (`44.0%` pooled recall). Weak sinusoids (`SNR=6–12`) fragment below `bin_seed_threshold = 4.0`.
7. **Unmodeled impulsive interference & absence of a glitch veto:**
   - As analyzed in Section 5, single-channel excess power cannot separate 1–3 sample impulsive glitches from real broadband pulses, yielding `84.0%` FAR under impulsive contamination.
8. **Windowing robustness test uses a single calibration split per seed:**
   - In `run_robustness_suite`, `boxcar` and `blackman` are calibrated on the same `cal_split` as `hann` (valid for held-out comparison on `win_split`, but `win_split` has only `300` pooled samples across 3 seeds and only tests off-bin sinusoids in pink noise, not all families).
9. **Empirical p-value floor & multiple-observation horizon:**
   - With `N_cal = 600`, the smallest resolvable empirical p-value is `1/601 = 0.00166`; strong detections (`score > 63`) cannot be ranked by p-value without a larger calibration ensemble. Furthermore, 1% FAP applies *per 4-second observation*, not across an entire multi-observation campaign.
