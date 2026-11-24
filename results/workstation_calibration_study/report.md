# WHISKER / BSD Controlled Calibration Study: 120 vs. 200 vs. 220 Nested Background Windows

**Purpose:** Educational, reproducible controlled experiment extending the original `120 vs. 200` calibration study with a third nested condition (`N_bg = 220`) to test the transition from `PROVISIONAL` to `VALIDATED` calibration status.
**Controlled Design Verification:**
- **Original 120 vs. 200 baseline preserved verbatim:** All `120`-window and `200`-window measurements match the original study bit-for-bit.
- **Strictly nested backgrounds (`np.array_equal == True`):** `bg_200[:120] == bg_120`, `bg_220[:120] == bg_120`, and `bg_220[:200] == bg_200` (adding exactly `20` new independent `pink_colored` windows).
- **Bitwise-identical target:** All `60` target windows (`40` pure-noise negative controls + `10` `linear_chirp` injections at `SNR = 13`, including Injection #8 in Window `57` spanning `t_rel = 1.316–3.229 s`, `112.97–210.10 Hz`) are strictly identical across all three conditions.
- **Safety:** Zero detector math, GUI, historical reports, or reserved O3 holdout files were modified or accessed.

## 1. Side-by-Side Comparison Table (`120` vs. `200` [Original Baseline] vs. `220` [Extension])

| Metric | 120 Windows *(Orig.)* | 200 Windows *(Orig.)* | 220 Windows *(Ext.)* |
| :--- | :---: | :---: | :---: |
| **Total Background Windows (4 s each)** | `120` (480 s) | `200` (800 s) | `220` (880 s) |
| **Split A (Per-Scale Standardization, 30%)** | `36` | `60` | `66` |
| **Split B (Threshold Calibration, 70%)** | `84` | `140` | **`154`** (`>= 150`!) |
| **Calibrated Threshold (`stft_threshold`, dimensionless)** | **`2.2568`** | **`2.4250`** | **`2.3089`** |
| **Empirical p-Value Floor (`p_min = 1/(N_th+1)`)** | `1.18%` (`1/85`) | `0.71%` (`1/141`) | **`0.65%`** (`1/155`) |
| **Effective Sample Size (`N_eff` / `N_th`)** | `46.8 / 84` | `125.4 / 140` | **`141.5 / 154`** |
| **Lag-1 Autocorrelation / Ljung-Box p** | `+0.039` / `0.3321` | `-0.018` / `0.7459` | `+0.004` / `0.6319` |
| **90% Beta FAR CI (i.i.d.)** | `[0.06%, 3.50%]` | `[0.04%, 2.12%]` | `[0.03%, 1.93%]` |
| **90% Beta FAR CI (`N_eff` dependence-aware)** | `[0.11%, 6.20%]` | `[0.04%, 2.36%]` | **`[0.04%, 2.10%]`** |
| **Workstation Calibration Status** | **`PROVISIONAL`** | **`PROVISIONAL`** | **`VALIDATED`** |
| **Injection #8 (Win 57, SNR=13) Score** | **`4.1188`** | **`4.5251`** | **`4.6796`** |
| **Injection #8 Decision (Margin)** | **DETECTED** (`+1.8620`) | **DETECTED** (`+2.1001`) | **DETECTED** (`+2.3707`) |
| **All 10 Injections Detected (SNR=13)** | **`10 / 10`** | **`10 / 10`** | **`10 / 10`** |
| **Negative Control FAR (Same 40 Windows)** | `2/40` (`5.0%`, CI `[1.7%, 14.0%]`) | `2/40` (`5.0%`, CI `[1.7%, 14.0%]`) | `2/40` (`5.0%`, CI `[1.7%, 14.0%]`) |

## 2. Does the Workstation Reach `VALIDATED` at 220 Windows? (Condition-by-Condition Audit)

Result at `N_bg = 220`: **`VALIDATED`**.

| Validation Criterion | Required | 120w | 200w | 220w |
| :--- | :---: | :---: | :---: | :---: |
| `n_thresh` (Split B size) | `>= 150` | `84` (FAIL) | `140` (FAIL) | **`154`** (PASS) |
| `n_eff` (Bartlett effective samples) | `>= 120.0` | `46.8` (FAIL) | `125.4` (PASS) | **`141.5`** (PASS) |
| `Ljung-Box p` (temporal independence) | `>= 0.05` | `0.3321` (PASS) | `0.7459` (PASS) | **`0.6319`** (PASS) |
| Upper 90% Beta FAR CI (`N_eff`) | `<= 2.50%` | `6.20%` (FAIL) | `2.36%` (PASS) | **`2.10%`** (PASS) |

- `120`w unmet conditions: `n_thresh=84 < 150; n_eff=46.8 < 120; dependence-aware upper 90% FAR CI=6.20% > 2.5%`
- `200`w unmet conditions: `n_thresh=140 < 150`
- `220`w unmet conditions: `None`

## 3. All 10 Linear-Chirp Injections (SNR = 13, Bitwise-Identical Target)

| Inj # | Win | Freq Band (Hz) | 120w Score (Dec) | 200w Score (Dec) | 220w Score (Dec) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| #0 | 41 | 269.7–439.2 | 59.5494 (**DET**) | 50.5059 (**DET**) | 40.7074 (**DET**) |
| #1 | 43 | 212.2–331.8 | 48.7792 (**DET**) | 41.3526 (**DET**) | 33.3397 (**DET**) |
| #2 | 45 | 208.5–352.1 | 25.6492 (**DET**) | 21.6950 (**DET**) | 17.5169 (**DET**) |
| #3 | 47 | 169.5–269.6 | 24.8276 (**DET**) | 21.6216 (**DET**) | 22.3596 (**DET**) |
| #4 | 49 | 129.3–288.5 | 33.6980 (**DET**) | 28.5355 (**DET**) | 23.0230 (**DET**) |
| #5 | 51 | 103.0–274.2 | 10.4541 (**DET**) | 8.7811 (**DET**) | 7.1223 (**DET**) |
| #6 | 53 | 79.3–255.6 | 18.7410 (**DET**) | 15.8239 (**DET**) | 12.7911 (**DET**) |
| #7 | 55 | 307.5–409.4 | 68.8991 (**DET**) | 58.4520 (**DET**) | 47.1034 (**DET**) |
| #8 **(Inj #8)** | 57 | 113.0–210.1 | 4.1188 (**DET**) | 4.5251 (**DET**) | 4.6796 (**DET**) |
| #9 | 59 | 66.7–238.6 | 10.4823 (**DET**) | 8.8051 (**DET**) | 7.1416 (**DET**) |

## 4. Split-A Per-Scale Standardization Factors (`median` / `spread = Q95 - median`)

| Scale | 120w Median / Spread | 200w Median / Spread | 220w Median / Spread |
| :--- | :---: | :---: | :---: |
| `short_31ms` | `54.87` / `33.67` | `56.65` / `30.60` | `56.65` / `31.39` |
| `medium_125ms` | `28.69` / `12.07` | `30.16` / `14.20` | `29.20` / `17.64` |
| `long_500ms` | `30.29` / `19.86` | `30.51` / `18.03` | `30.51` / `17.43` |

## 5. Figures & Limitations

- `figures/fig01_background_score_distributions.png`
- `figures/fig02_target_timeline_120_200_220.png`
- `figures/fig03_all_10_injections_and_uncertainty.png`

**Explicit Methodological Clarifications & Limitations:**
- **Internal calibration validation vs. real-world FAR control:** Reaching `VALIDATED` at `N_bg = 220` confirms only that the *internal calibration sample* meets the workstation's statistical adequacy criteria (`n_thresh >= 150`, `n_eff >= 120`, `Ljung-Box p >= 0.05`, dependence-aware upper 90% Beta FAR CI `<= 2.5%`). It is **not** independent proof of `<= 1%` false-alarm control on future real-world or nonstationary data.
- **Unchanged negative-control FAR (`2/40 = 5.0%`, 90% Wilson CI `[1.7%, 14.0%]`):** Across all three nested calibrations (`120`, `200`, and `220` windows), the observed false-alarm rate on the identical 40-window pure-noise negative control remained `2/40` (`5.0%`), above the nominal `1%` target with a wide confidence interval (`2.5%` discrete step per window). Independently verifying `1%` FAR with tight precision requires hundreds of held-out null windows, not `40`.
- **Limitations of `N_eff` and Ljung–Box diagnostics:** A non-significant Ljung–Box p-value (`p >= 0.05`) and high Bartlett effective sample size (`N_eff = 141.5 / 154`) indicate only that no strong *linear* lag-1..10 autocorrelation was detected in the calibration scores — they are **never** mathematical proof of full statistical independence (nonlinear dependence, non-Gaussian glitch clustering, or slow drift longer than 10 windows can still exist).
- **Single-realization scope & split boundary shift:** This study compares one nested synthetic realization (`seed = 2026`, `pink_colored`, `linear_chirp`, `SNR = 13`) to isolate calibration mechanics, not ensemble-averaged detection curves. Because `MultiResolutionDetector.calibrate` uses a `30%/70%` chronological split, changing `N_bg` shifts both Split A (`0..35`, `0..59`, `0..65`) and Split B (`36..119`, `60..199`, `66..219`) simultaneously, causing both the Split-A normalization parameters and the Split-B order-statistic threshold (`2.2568 -> 2.4250 -> 2.3089`) to change together.
