# WHISKER / Blind Signal Discovery (BSD) — Mathematical Measurement Guide & Source Code Audit

**Scope:** Exact mathematical documentation of the current BSD Local Measurement Workstation implementation (`v0.3 Phase 2b`, commits `edd1ec9`–`40c7bca`), verified directly against source code in `src/bsd/` and `app/bsd_workstation.py`.
**Policy:** Zero detector math, calibration, or GUI code was modified during this audit. Reserved O3 holdout data (`docs/v03_o3_holdout_manifest.json`) remains 100% untouched.

---

## 1. Concise Reference Table of All GUI Quantities, Units & Calculations

| GUI Quantity / Label | Symbol | Units / Status | Exact Formula in Code | Source Implementation |
| :--- | :---: | :---: | :--- | :--- |
| **Waveform Amplitude** | $x[n]$ | Synthetic: **dimensionless** (`normalized_sigma`, $\sigma_{\text{noise}} = 1$); WAV PCM: normalized `[-1, 1]`; GWOSC: `strain × 1e18`; CSV/NPY: user-declared | Direct 1D float64 array `obs.signal` | `bsd/io_loader.py`, `bsd/workstation_runner.py`, `app/bsd_workstation.py:_plot_window_waveform_and_spec` |
| **Relative Window Time** | $t_{\text{rel}}$ | Seconds (`s`), $[0, T_w)$ | $t_{\text{rel}}[n] = n / f_s$, where $T_w = 4.0\text{ s}$, $f_s = 1024\text{ Hz}$ | `app/bsd_workstation.py` (`t = np.arange(N) / fs`) |
| **Absolute Timeline Time** | $t_{\text{abs}}$ | Seconds (`s`) from target start | $t_{\text{abs}} = k \cdot T_w + t_{\text{rel}}$ for window index $k$ | `bsd/workstation_runner.py:_build_rows` |
| **Frequency** | $f$ | Hertz (`Hz`), $[0, f_s / 2]$ | STFT frequency grid $f_m = m \cdot (f_s / N_{\text{seg}})$ | `scipy.signal.stft` in `bsd/detector.py`, `app/bsd_workstation.py` |
| **Spectrogram Power (`Power (dB)`)** | $P_{\text{dB}}(t,f)$ | **Relative dB** (`10 log10(|Z|^2)` re $1\text{ unit}^2$, NOT `dB/Hz` or `dB SPL`) | $10 \log_{10}\!\big(\max(|Z[f,t]|^2,\, 10^{-24})\big)$ with `nperseg=128`, `noverlap=96`, Hann window | `app/bsd_workstation.py:_plot_window_waveform_and_spec` |
| **Noise RMS ($\sigma$)** | $\sigma_{\text{noise}}$ | Same units as $x[n]$ (synthetic: `1.0` dimensionless) | Synthetic noise explicitly RMS-normalized: $\sqrt{\frac{1}{N}\sum x[n]^2} = \sigma$ | `bsd/noise.py:_normalize_rms` |
| **Injected Time-Domain Energy SNR** | $\mathrm{SNR}_E$ | **Dimensionless linear root-energy ratio** (NOT power ratio, NOT dB) | $\mathrm{SNR}_E = \frac{\sqrt{\sum_{n=0}^{N-1} s[n]^2}}{\sigma_{\text{noise}}} = \frac{\|s\|_2}{\sigma_{\text{noise}}}$ | `bsd/generator.py:_scale_to_energy_snr` |
| **Per-Bin Normalized Excess Power** | $Z[f,t]$ | **Dimensionless** power ratio | $Z[f,t] = \frac{|S[f,t]|^2}{B_{\text{psd}}[f]\, B_{\text{time}}[t]}$ (2D median-normalized) | `bsd/detector.py:TemplateFreeDetector.analyze` |
| **Raw Cluster Anomaly Score** | $A_{\text{raw}}$ | **Dimensionless** excess-power sum | $A_{\text{raw}}(C) = \sum_{(f,t) \in C} \big(Z[f,t] - 1\big)$ over 8-connected bins with $Z > 4.0$, $|C| \ge 2$ | `bsd/detector.py` (lines 187–204) |
| **Standardized Global Score (`stft_max_score`)** | $S_{\text{global}}$ | **Dimensionless** null-tail standardized units | $S_{\text{global}} = \max_{r \in \text{scales}} \frac{A_{\text{raw},r} - m_r}{\max(q_{95,r} - m_r,\, 1.0)}$ | `bsd/detector_v2.py:ScaleNormalizationStats.standardize` & `.analyze` |
| **Active Threshold (`stft_threshold`)** | $\tau_{\text{global}}$ | **Dimensionless** (exact same scale as $S_{\text{global}}$) | $k$-th order statistic $S_{(k)}$ of Split-B null global scores, $k = \min(N_{\text{th}}, \lceil (1-\alpha)(N_{\text{th}}+1) \rceil)$ | `bsd/detector_v2.py:MultiResolutionDetector.calibrate` |
| **Empirical p-Value (`empirical_p_value`)** | $p_{\text{emp}}$ | **Dimensionless** probability in $[p_{\min}, 1]$ | $p_{\text{emp}}(s) = \frac{1 + \#\{i : S_i \ge s\}}{N_{\text{th}} + 1}$, resolution floor $p_{\min} = \frac{1}{N_{\text{th}} + 1}$ | `bsd/detector_v2.py:MultiResolutionCalibration.empirical_p_value` |
| **Target FAP vs. Observed FAR** | $\alpha$ vs. $\widehat{\mathrm{FAR}}$ | **Dimensionless** probability (often displayed as `%`) | Target $\alpha \in (0,1)$; Observed $\widehat{\mathrm{FAR}} = \frac{N_{\text{false}}}{N_{\text{null}}}$ with 90% Wilson CI | `bsd/workstation_runner.py:_wilson_ci90` |
| **Effective Sample Size (`N_eff`)** | $N_{\text{eff}}$ | Effective independent windows (dimensionless) | $N_{\text{eff}} = \frac{N_{\text{th}}}{1 + 2\sum_{k=1}^{L} (1 - \frac{k}{L+1})\max(\hat{\rho}_k, 0)}$ | `bsd/adaptive_regime_v2b.py:estimate_calibration_uncertainty` |
| **Frozen-Ref Broadband Score** | $E_{\text{bb}}$ | **Dimensionless** positive excess $\chi^2$ sum | $\sum_{b=1}^{4} \max\!\big(0, \frac{\log \hat{S}_b - \mu_{b,R_0}}{\sigma_{b,R_0}}\big)^2$ over 4 continuum bands | `bsd/adaptive_regime_v2b.py:FrozenReferenceBroadbandChannel.score` |

---

## 2. Physical Quantities, Units, Assumptions & Limitations

### 2.1 Waveform Amplitude and Noise Statistics
- **Verified implementation (`src/bsd/noise.py`):** Synthetic noise generators (`gaussian_white`, `pink_colored`, `nonstationary`, `impulsive`, `shifted_ar2`) generate zero-mean arrays and pass them through `_normalize_rms(x, sigma)` so that $\sqrt{\frac{1}{N}\sum_{n=0}^{N-1} x[n]^2} = \sigma$ exactly (`sigma = 1.0` in the self-test). Synthetic waveforms have **no physical units** (they are dimensionless normalized units `normalized_sigma`).
- **Real files (`src/bsd/io_loader.py`):** Integer PCM `.wav` files are scaled to `[-1, 1]` (`PCM_int16_normalized`); GWOSC `.hdf5` files are scaled by $10^{18}$ into `O(1)` float64 units (`GWOSC_strain_x1e18_unwhitened`); `.csv`/`.npy` files retain user-supplied units without silent rescaling.
- **Pink colored noise (`generate_pink_colored`):** Shapes white Gaussian FFT bins by $H(f) = (f^2 + f_{\text{knee}}^2)^{-\alpha/4}$ with $f_{\text{knee}} = 16\text{ Hz}$, $\alpha \sim \mathrm{Uniform}(0.85, 1.15)$ per window, DC bin zeroed, then RMS-normalized to $\sigma = 1.0$.

### 2.2 Time and Frequency Coordinates
- **Window slicing (`slice_into_observations`):** Non-overlapping segments of $T_w = 4.0\text{ s}$ at $f_s = 1024.0\text{ Hz}$ ($N = 4096$ samples/window).
- **Absolute vs. relative time:** Window $k$ covers $[k \cdot T_w,\, (k+1)\cdot T_w)$. Inside window $k$, relative time runs $t_{\text{rel}} \in [0, 4.0)\text{ s}$. Candidate timestamps in `candidates.csv` are converted to absolute time $t_{\text{abs}} = k \cdot T_w + t_{\text{rel}}$.
  - **Worked example (Window 57):** $k = 57 \implies$ window spans **$t_{\text{abs}} \in [228.000\text{ s},\, 232.000\text{ s})$**. The injected chirp spans $t_{\text{rel}} \in [1.3162\text{ s},\, 3.2291\text{ s}]$, which corresponds to **$t_{\text{abs}} \in [229.316\text{ s},\, 231.229\text{ s}]$**.

### 2.3 Spectral Power and dB Conversion
- **Verified implementation (`app/bsd_workstation.py`, lines 62–72):** The GUI spectrogram computes `scipy.signal.stft(x - median(x), fs=1024, window="hann", nperseg=128, noverlap=96, boundary=None, padded=False)` (frequency step $\Delta f = 8\text{ Hz}$, hop $\Delta t = 31.25\text{ ms}$) and plots:
  $$P_{\text{dB}}[f, t] = 10 \log_{10}\!\big(\max(|Z[f,t]|^2,\, 10^{-24})\big)$$
- **Important distinction:** `scipy.signal.stft` returns complex windowed Fourier amplitudes normalized by the window sum (`scaling="spectrum"` default). Therefore $|Z[f,t]|^2$ is **uncalibrated squared spectral amplitude** (`unit^2`), NOT one-sided power spectral density (`unit^2 / Hz`) and NOT acoustic `dB SPL`. Negative dB values simply mean $|Z[f,t]|^2 < 1\text{ unit}^2$.

---

## 3. Mathematical Explanation of All Three GUI Plots (Worked Example: Window 57)

All numerical values below were **verified directly by executing** `run_synthetic_self_test` with `noise_model="pink_colored"`, `signal_family="linear_chirp"`, `n_bg_windows=120`, `n_neg_control_windows=40`, `n_pos_injections=10`, `seed=2026`, `window_duration_s=4.0`, `sample_rate=1024 Hz`, `target_fap=0.01`, `calibration_seed=0`.

1. **Plot 1 — Waveform (`_plot_window_waveform_and_spec`, top panel):**
   - **X-axis:** `Time Within Window (s)` — relative time $t_{\text{rel}} \in [0, 4.0)\text{ s}$ (`4096` samples).
   - **Y-axis:** `Amplitude` — raw time-series $x[n] = n_{57}[n] + s_{57}[n]$ in dimensionless normalized noise units ($\sigma_{\text{noise}} = 1.0$).
   - **Window 57 values:** Noise RMS $= 1.0000$. At `SNR = 12`, peak signal amplitude is `0.4029` (signal RMS `0.1875`) — completely buried below the `1.0` RMS noise floor in the time domain! At `SNR = 18`, peak signal amplitude is `0.6043` (signal RMS `0.2812`) — still visually submerged in the time-domain waveform.
2. **Plot 2 — Spectrogram (`_plot_window_waveform_and_spec`, bottom panel):**
   - **X-axis:** `Time Within Window (s)` ($t_{\text{rel}} \in [0, 4.0)\text{ s}$).
   - **Y-axis:** `Frequency (Hz)` ($[0, 512]\text{ Hz}$).
   - **Colorbar:** `Power (dB)` $= 10 \log_{10}(\max(|Z[f,t]|^2, 10^{-24}))$.
   - **Red bounding boxes:** Drawn for every accepted candidate in `win.raw_stft_candidates` using `(c.t_start, c.f_min)` to `(c.t_end, c.f_max)` in relative window coordinates.
   - **Window 57 behavior:** The injected chirp sweeps linearly from `112.97 Hz` to `210.10 Hz` over `1.316 s` to `3.229 s`. Because pink noise concentrates power at low frequencies and drops off above `100 Hz`, the diagonal chirp track stands out in the `110–210 Hz` band even though it is invisible in the time-domain waveform. At `SNR = 12`, zero boxes are drawn (missed); at `SNR = 18`, 4 multi-scale candidate boxes are drawn, led by `long_500ms` spanning `[1.562 s, 3.312 s] × [131.0 Hz, 207.0 Hz]`.
3. **Plot 3 — Full-Timeline Anomaly Score vs. Threshold (`_plot_timeline`):**
   - **X-axis:** `Window Index` $k \in \{0, 1, \dots, 59\}$ (`0..39` are the 40 pure-noise negative-control windows; `41, 43, 45, ..., 59` are the 10 positive-injection windows, with Window `57` being injection `#8`).
   - **Y-axis:** `STFT Max Global Score` ($S_{\text{global}}$, dimensionless) in `transient` mode, or `Frozen-Ref Broadband Excess Score` ($E_{\text{bb}}$, dimensionless) in `persistent` mode.
   - **Red dashed line:** `Active Empirical Threshold` ($\tau_{\text{global}} = 2.2568$ for both `SNR = 12` and `SNR = 18`).
   - **Shading:** Orange = `CONTINUUM_DRIFT` warning; Purple = `GLITCH_STORM` warning.

---

## 4. Exact SNR Formula & Determinism Audit (`SNR = 12` vs. `SNR = 18`)

### 4.1 Exact Implementation (`src/bsd/generator.py`, lines 105–114)
```python
def _scale_to_energy_snr(clean_signal, noise_sigma, target_energy_snr):
    energy_norm = float(np.linalg.norm(clean_signal))
    return clean_signal * (target_energy_snr * noise_sigma / energy_norm)
```
The exact mathematical formula is:
$$s_{\text{scaled}}[n] = s_{\text{clean}}[n] \cdot \frac{\mathrm{SNR}_E \cdot \sigma_{\text{noise}}}{\|s_{\text{clean}}\|_2}, \qquad \text{so that} \qquad \mathrm{SNR}_E = \frac{\sqrt{\sum_{n=0}^{N-1} s_{\text{scaled}}[n]^2}}{\sigma_{\text{noise}}}$$

- **Type of ratio:** $\mathrm{SNR}_E$ is a **linear root-energy-to-RMS-noise ratio** (dimensionless, **not** in dB and **not** a squared power ratio). In white Gaussian noise of variance $\sigma_{\text{noise}}^2$, $\mathrm{SNR}_E$ equals the exact discrete matched-filter optimal SNR $\rho_{\text{opt}}$. In **colored** noise (`pink_colored` or `shifted_ar2`), $\mathrm{SNR}_E$ is an **unweighted** time-domain energy ratio and can differ substantially from the true colored matched-filter SNR $\sqrt{4\int |H(f)|^2 / S_n(f)\,df}$!
- **Worked numerical example (Window 57):**
  - Unscaled chirp norm: $\|s_{\text{clean}}\|_2 = 29.786218$; $\sigma_{\text{noise}} = 1.0$.
  - At `SNR = 12`: scale factor $= 12.0 / 29.786218 = 0.402871 \implies \|s_{12}\|_2 = 12.000$, energy $\sum s_{12}^2 = 144.0$, peak $= 0.4029$, window RMS $= 12.0 / \sqrt{4096} = 0.1875$.
  - At `SNR = 18`: scale factor $= 18.0 / 29.786218 = 0.604306 \implies \|s_{18}\|_2 = 18.000$, energy $\sum s_{18}^2 = 324.0$, peak $= 0.6043$, window RMS $= 18.0 / \sqrt{4096} = 0.2812$.
  - **Ratio (`18` vs `12`):** Amplitude increases by **$1.50\times$**; signal energy and STFT power increase by **$1.50^2 = 2.25\times$** (**$+3.52\text{ dB}$**).

### 4.2 Verified Determinism Audit (`seed = 2026`)
Tracing `run_synthetic_self_test` (`src/bsd/workstation_runner.py`, lines 350–380) and verifying via `np.array_equal`:
1. **Identical across `SNR = 12` and `SNR = 18`:** Because `rng = make_split_rng(2026, "workstation_self_test")` draws all 120 background windows, all 60 target noise windows, and all 10 unscaled `generate_clean_signal` waveforms in the exact same order **before** multiplying by `target_snr`, the following are **100% bit-for-bit identical** between `SNR = 12` and `SNR = 18`:
   - All 120 background calibration windows, `N_thresh = 84`, `N_eff = 46.81`, `p_LB = 0.3321`, and global threshold **`2.2568`**.
   - All 40 negative-control windows (same `2/40 = 5.0%` false alarms) and underlying target noise in Window 57.
   - The unscaled chirp waveform, start/end times (`1.3162–3.2291 s`), and frequency trajectory (`112.97–210.10 Hz`).
2. **What changes:** Strictly the linear amplitude multiplier on the 10 injected chirps (`1.50×` amplitude, `2.25×` power).

---

## 5. Exact Derivation of the STFT Detection Statistic

### Step 1 — Multi-Resolution STFT Bank (`src/bsd/detector_v2.py`, lines 31–68)
Each 4-second window is analyzed independently at **three Hann-windowed STFT scales** (DC and Nyquist bins dropped):
- `short_31ms`: `nperseg=32`, `noverlap=28` ($\Delta f = 32\text{ Hz}$, $\Delta t = 3.9\text{ ms}$, `freq_k=5`, `time_k=31`)
- `medium_125ms`: `nperseg=128`, `noverlap=96` ($\Delta f = 8\text{ Hz}$, $\Delta t = 31.25\text{ ms}$, `freq_k=11`, `time_k=15`)
- `long_500ms`: `nperseg=512`, `noverlap=384` ($\Delta f = 2\text{ Hz}$, $\Delta t = 125\text{ ms}$, `freq_k=21`, `time_k=7`)

### Step 2 — Within-Window 2D Robust Background Normalization (`src/bsd/detector.py`, lines 169–184)
For each scale, let $P[f, t] = |S[f,t]|^2$. Under stationary Gaussian noise, $P[f,t]$ is exponentially distributed with median $= \ln(2) \cdot \mathbb{E}[P[f,t]]$:
1. Spectral baseline: $\hat{S}_{\text{raw}}[f] = \mathrm{median}_t(P[f,t]) / \ln 2$, smoothed across frequency by a median filter of width `freq_k` to obtain $B_{\text{psd}}[f]$.
2. Spectral whitening: $W[f,t] = P[f,t] / B_{\text{psd}}[f]$.
3. Temporal floor: $\hat{E}_{\text{raw}}[t] = \mathrm{median}_f(W[f,t]) / \ln 2$, smoothed in time by a median filter of width `time_k` to obtain $B_{\text{time}}[t]$.
4. Dimensionless excess-power map: $Z[f,t] = W[f,t] / B_{\text{time}}[t]$ (under stationary Gaussian noise, $Z[f,t] \approx \mathrm{Exp}(1)$ with mean $1$).

### Step 3 — Connected-Component Clustering & Raw Score (`src/bsd/detector.py`, lines 187–204)
- Threshold bins at $Z[f,t] > 4.0$ (`bin_seed_threshold = 4.0`), group into 8-connected components, keep clusters with $|C| \ge 2$ bins (`min_cluster_bins = 2`).
- Each cluster's dimensionless **raw anomaly score** is:
  $$A_{\text{raw}}(C) = \sum_{(f,t) \in C} \big(Z[f,t] - 1.0\big)$$
- Window max raw score at scale $r$: $A_{\text{raw},r} = \max_{C} A_{\text{raw}}(C)$ (or `0.0` if no qualified cluster).

### Step 4 — Split-Sample Per-Scale Standardization & Global Decision (`src/bsd/detector_v2.py`, lines 238–359)
The background windows are split chronologically into **Split A (`30%`, `n_norm`)** and **Split B (`70%`, `n_thresh`)**:
- On **Split A**, for each scale $r$, compute null median $m_r = \mathrm{median}(A_{\text{raw},r})$ and upper-tail spread $s_r = \max(Q_{0.95}(A_{\text{raw},r}) - m_r,\, 1.0)$.
- Standardize every cluster and window max into **dimensionless null-tail units**:
  $$S_r = \frac{A_{\text{raw},r} - m_r}{s_r}, \qquad S_{\text{global}} = \max_{r \in \{\text{short},\,\text{medium},\,\text{long}\}} S_r$$
- On **Split B**, collect $N_{\text{th}}$ null global scores $S_{\text{global}, 1}, \dots, S_{\text{global}, N_{\text{th}}}$, sort them, and set the **calibrated threshold**:
  $$\tau_{\text{global}} = S_{(k)}, \qquad k = \min\!\Big(N_{\text{th}},\, \big\lceil (1 - \alpha)(N_{\text{th}} + 1) \big\rceil\Big)$$
- **Decision rule:** A candidate is nominated if and only if $S_r(C) > \tau_{\text{global}}$ (strict inequality). Its conformal empirical p-value is $p_{\text{emp}}(s) = (1 + \#\{i : S_i \ge s\}) / (N_{\text{th}} + 1)$.

### Step 5 — Direct Comparability of `23.7075` vs. `2.2568` (Verified!)
- **Yes — directly comparable!** In Window 57 at `SNR = 18`, the top candidate (`long_500ms`) had `raw_anomaly_score = 501.1308` (`72` connected bins, `peak_excess_power = 18.81`), which standardized to **`standardized_score = 23.7075`** (`stft_max_score = 23.7075`). The calibrated threshold **`stft_threshold = 2.2568`** is in the **exact same dimensionless standardized scale**. Because `23.7075 > 2.2568` (and exceeds all `84` Split-B null scores), it is nominated with $p_{\text{emp}} = 1 / (84 + 1) = 0.011765$ (`p_value_is_upper_bound = True`, meaning $p \le 0.0118$).
- **Why `SNR = 12 -> 18` caused a 13× score jump (`1.8018 -> 23.7075`):** Because of the hard per-bin seed cutoff ($Z > 4.0$) and 8-connected clustering. At `SNR = 12`, only scattered bins along the chirp ridge crossed $Z > 4.0$, giving $S_{\text{global}} = 1.8018 < 2.2568$ (missed). Increasing signal power by $2.25\times$ pushed the continuous ridge above $Z > 4.0$, merging into a `72`-bin cluster (`raw = 501.13`, `std = 23.71`).

---

## 6. Calibration, Uncertainty & Why a Missed Injection ≠ Signal Absence

- **Target FAP ($\alpha$) vs. Observed FAR ($\widehat{\mathrm{FAR}}$):** $\alpha = 0.01$ is the *design target* false-alarm probability. With finite calibration size $N_{\text{th}} = 84$ (from `120` background windows with `30%` Split A = `36` and `70%` Split B = `84`), $k = \min(84, \lceil 0.99 \times 85 \rceil) = 84$ (the sample maximum!).
- **Empirical p-value resolution vs. true FAR uncertainty:** $p_{\min} = 1/(84 + 1) = 1.18\%$ is strictly the discrete resolution floor of $p_{\text{emp}}$. Under i.i.d. sampling, the true population FAR of $\tau_{\text{global}}$ follows $\mathrm{Beta}(1, 84)$ with 90% CI `[0.06%, 3.50%]`. Accounting for observed lag autocorrelation (`N_eff = 46.81`, `Ljung-Box p = 0.3321`), the dependence-aware 90% Beta FAR CI widens to **`[0.11%, 6.20%]`**, which is completely consistent with the observed `2/40 = 5.0%` false alarms on the negative control!
- **`PROVISIONAL` vs. `VALIDATED` (`classify_calibration_status`):**
  - `INSUFFICIENT`: `n_bg < 60` or `N_eff < 45`.
  - `VALIDATED`: requires ALL four conditions: `N_thresh >= 150`, `N_eff >= 120`, `Ljung-Box p >= 0.05`, AND dependence-aware upper 90% FAR CI `<= 2.5%`.
  - `PROVISIONAL`: anything in between. In our `120`-window self-test (`N_th = 84`, `N_eff = 46.8`, upper CI `6.2%`), status is correctly **`PROVISIONAL`**.
- **Why a missed injection does NOT prove a signal is absent:** Window 57 at `SNR = 12` contained a real `linear_chirp` (`||s||_2 = 12.0`), yet scored `1.8018 < 2.2568`. A non-detection only means the signal's local STFT excess-power clusters did not exceed the calibrated `1%` FAP threshold under the current background noise realization — never that no signal was present.

---

## 7. Interface & Scientific Clarity Audit (Findings & Recommendations — No Code Modified Yet)

Inspecting `app/bsd_workstation.py` revealed **six concrete labeling ambiguities** to improve in a future UI polish commit:

1. **Spectrogram colorbar (`Power (dB)`):** Currently implies calibrated physical power or PSD (`dB/Hz`). **Recommendation:** Relabel to `Relative STFT Power (10 log10 |Z|^2, dB re 1 unit^2)` and note `nperseg=128, Hann`.
2. **Waveform Y-axis (`Amplitude`):** Omits units. **Recommendation:** Dynamically render `Amplitude ({units})` from `run_metadata["target_preprocessing"]["units"]` (e.g., `Amplitude (normalized_sigma)` for synthetic runs, `Amplitude (PCM_int16_normalized)` for WAV).
3. **Timeline Y-axis (`STFT Max Global Score`):** Does not state that the quantity is dimensionless standardized null-tail units directly comparable to the threshold. **Recommendation:** Label `Standardized Global Anomaly Score (dimensionless)` and show the exact numeric threshold in the legend (`Active Threshold = 2.257`).
4. **Candidate table (`raw_anomaly_score` vs. `standardized_score`):** Both columns appear side-by-side without explicitly stating that `standardized_score` (not `raw_anomaly_score`) is what is compared against `stft_threshold`. **Recommendation:** Add a short caption above the candidate table clarifying `Decision compares standardized_score > stft_threshold`.
5. **Self-test slider label (`Injected Time-Domain Energy SNR`):** Users may mistake `12` or `18` for decibels (`12 dB`). **Recommendation:** Clarify `Injected Root-Energy SNR (||s||_2 / sigma, linear dimensionless — NOT dB)`.
6. **Window Inspector time axis (`Time Within Window (s)`):** Does not display the corresponding absolute timeline span `[228.0 s, 232.0 s]` in the panel title. **Recommendation:** Include `Absolute Time: [t0_abs, t1_abs] s` in the window title so users can directly cross-reference `candidates.csv` timestamps.

---

## 8. Beginner-Friendly Glossary

- **STFT (Short-Time Fourier Transform):** Slices a time series into short overlapping frames and measures how much energy is present at each frequency in each frame, producing a 2D time-frequency image (spectrogram).
- **Excess Power ($Z[f,t]$):** How many times stronger a spectrogram pixel is than the local background noise floor ($Z = 1$ is average noise; $Z = 18$ is $18\times$ average noise power).
- **Connected Cluster:** Adjacent bright spectrogram pixels ($Z > 4.0$) touching horizontally, vertically, or diagonally, grouped into a single localized candidate event.
- **Root-Energy SNR ($\|s\|_2 / \sigma$):** A linear, dimensionless measure of total signal strength relative to 1-sample noise RMS. Multiplying amplitude by $1.5\times$ multiplies SNR by $1.5\times$ and power by $2.25\times$ ($+3.52\text{ dB}$).
- **Target FAP vs. Empirical FAR:** Target FAP is the false-alarm rate you *ask for* (`1%`); empirical FAR is the fraction of clean background windows that *actually* cross your finite-sample threshold.
- **Effective Sample Size ($N_{\text{eff}}$):** How many *truly independent* background windows your calibration set is worth after discounting window-to-window temporal correlation.
- **Empirical p-Value ($p_{\text{emp}}$):** The fraction of calibration background windows that scored at least as high as your candidate. With $N_{\text{th}} = 84$, the smallest possible value is $1 / 85 = 0.0118$.
