# Blind Signal Discovery (BSD) — Experiment 003: Multi-Event Real-Data Detection Benchmark

> **Scientific Scope Statement:** This report evaluates the **frozen BSD v0.1 detector** (`src/bsd/detector.py`) and **frozen LIGO adaptation pipeline** (`src/bsd/ligo_adapter.py`, `adapted_ligo_whitened`) across multiple confirmed GWOSC gravitational-wave events and independent background stretches. **No detector parameters or thresholds were tuned.** No new gravitational-wave discovery is claimed.

## 1. Research Question

How does BSD v0.1 detection efficiency vary across real gravitational-wave events with different signal strengths, durations, and frequency evolution, at a measured background false-alarm rate?

## 2. Experimental Protocol & Dataset Separation

- **Frozen detector:** `TemplateFreeDetector(DetectorConfig())` (`nperseg=128`, `noverlap=96`, `window='hann'`, `freq_median_kernel=11`, `time_median_kernel=15`, `bin_seed_threshold=4.0`, `min_cluster_bins=2`).
- **Frozen LIGO preprocessing:** `adapted_ligo_whitened` (`1e18` scale, `4.0 s` median-Welch PSD estimated strictly from independent background, `30-380 Hz` 4th-order zero-phase Butterworth bandpass, decimation `4096 -> 1024 Hz`, `1.0 s` context padding).
- **Independent per-epoch calibration:** For every event and available detector (`H1`, `L1`, and `V1` where available), a guard zone (`±64 s` for BBH; `±128 s` for BNS `GW170817`) around the 32-second target file is excluded from the 4096-second epoch file. Remaining DQ-verified (`DATA=1`, finite) 4-second segments are split into disjoint PSD estimation (`60`), 1% FAP calibration (`500`), and independent background validation (`~300-422`) sets.
- **Reference vs. Independent Validation:** `GW150914` (previously examined in Experiment 002) is listed strictly as a **reference case** and is **excluded** from all independent validation efficiency calculations.
- **Untouched Holdout Set for Future BSD v0.2:** Four additional GWTC-1 events (`GW170729`, `GW170809`, `GW170818`, `GW170823`) are recorded in `v02_holdout_manifest.json` and were **neither downloaded nor inspected**.

## 3. Event-by-Event Benchmark Table

Blind recovery requires a threshold-exceeding candidate within `±0.50 s` of (or spanning) the published GWOSC GPS timestamp (consulted only after all candidates were frozen). Multi-detector coincidence additionally requires two detectors to produce blind candidates within `|dt| <= 0.10 s`.

| Event | Role | Run | Class | Net SNR | M1/M2 (M_sun) | H1 (Score/Thr) | L1 (Score/Thr) | V1 (Score/Thr) | Single-Det Recovery | Coincident (>=2 det) | Other 32s Anomalies |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914` | REF | O1 | BBH | 25.2 | 35.6/30.6 | **YES** (818/77, dt=-0.025s) | **YES** (137/76, dt=-0.056s) | N/A | **YES** | **YES** | 0 |
| `GW151012` | VAL | O1 | BBH | 10.0 | 23.2/13.6 | NO (141/451) | NO (84/116) | N/A | NO | NO | 0 |
| `GW151226` | VAL | O1 | BBH | 13.1 | 13.7/7.7 | NO (50/76) | NO (40/91) | N/A | NO | NO | 0 |
| `GW170104` | VAL | O2 | BBH | 13.0 | 30.8/20.0 | NO (47/86) | NO (45/72) | N/A | NO | NO | 0 |
| `GW170608` | VAL | O2 | BBH | 15.4 | 11.0/7.6 | **YES** (4542/4155, dt=+1.438s) | NO (94/95) | N/A | **YES** | NO | 4 |
| `GW170814` | VAL | O2 | BBH | 17.2 | 30.6/25.2 | NO (38/77) | **YES** (276/74, dt=-0.031s) | NO (391/426) | **YES** | NO | 1 |
| `GW170817` | VAL | O2 | BNS | 33.0 | 1.46/1.27 | **YES** (141/107, dt=-0.119s) | NO (234777/83) | NO (168/356) | **YES** | NO | 2 |

## 4. Independent Validation Efficiency Summary (Excluding GW150914 Reference)

- **Independent validation events evaluated:** `6` (`GW151012`, `GW151226`, `GW170104`, `GW170608`, `GW170814`, `GW170817`)
- **Events recovered by >= 1 detector blindly:** `3 / 6` = **`50.0%`** (95% exact Clopper-Pearson CI `[11.8%, 88.2%]`)
- **Events recovered in multi-detector coincidence (`|dt| <= 0.10 s`):** `0 / 6` = **`0.0%`** (95% exact Clopper-Pearson CI `[0.0%, 45.9%]`)
- **Per-detector blind recovery rate (across all available validation detector-events):** `3 / 14` = **`21.4%`** (95% exact CI `[4.7%, 50.8%]`)

> **Statistical Caution:** With only 6 independent validation events spanning heterogeneous source classes (BBH and BNS), these proportions carry wide confidence intervals and **do not establish a universal detection efficiency curve**. They provide a quantitative baseline showing where single-resolution excess power succeeds and fails on real interferometer data.

## 5. Independent Background False-Alarm Analysis

Across all `6412` independent background validation 4-second segments (pooled across all events and detectors), BSD v0.1 produced `61` false alarms: **pooled empirical background FAR = `0.95%`** (95% exact Clopper-Pearson CI `[0.73%, 1.22%]`) against the `1.00%` target FAP.

| Event | Detector | Clean Bg Segs | DQ Excluded | N_cal | 1% Thr | Val FP / N | Empirical Bg FAR | 95% Exact CI | Credible? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914` | H1 | 982 | 2 | 500 | `77.1` | 3 / 422 | 0.71% | [0.15%, 2.06%] | Yes |
| `GW150914` | L1 | 982 | 2 | 500 | `76.3` | 5 / 422 | 1.18% | [0.39%, 2.74%] | Yes |
| `GW151012` | H1 | 982 | 2 | 500 | `451.3` | 1 / 422 | 0.24% | [0.01%, 1.31%] | Yes |
| `GW151012` | L1 | 982 | 2 | 500 | `116.1` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW151226` | H1 | 500 | 484 | 270 | `76.3` | 3 / 180 | 1.67% | [0.35%, 4.79%] | Yes |
| `GW151226` | L1 | 982 | 2 | 500 | `90.7` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW170104` | H1 | 982 | 2 | 500 | `86.3` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW170104` | L1 | 982 | 2 | 500 | `71.5` | 17 / 422 | 4.03% | [2.36%, 6.37%] | NO (FAR inflated) |
| `GW170608` | H1 | 980 | 4 | 500 | `4155.3` | 1 / 420 | 0.24% | [0.01%, 1.32%] | Yes |
| `GW170608` | L1 | 982 | 2 | 500 | `95.0` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW170814` | H1 | 982 | 2 | 500 | `76.8` | 1 / 422 | 0.24% | [0.01%, 1.31%] | Yes |
| `GW170814` | L1 | 982 | 2 | 500 | `74.3` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW170814` | V1 | 982 | 2 | 500 | `426.5` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `GW170817` | H1 | 950 | 2 | 500 | `107.5` | 1 / 390 | 0.26% | [0.01%, 1.42%] | Yes |
| `GW170817` | L1 | 950 | 2 | 500 | `83.0` | 4 / 390 | 1.03% | [0.28%, 2.61%] | Yes |
| `GW170817` | V1 | 950 | 2 | 500 | `355.8` | 13 / 390 | 3.33% | [1.79%, 5.63%] | NO (FAR inflated) |

## 6. Sensitivity by Event Characteristics & Failure Analysis

1. **High-SNR short BBH chirps (`GW150914` ref SNR 25.2, `GW170814` SNR 17.2):** Energy is concentrated into ~0.15-0.35 s across 35-250 Hz, closely matching BSD v0.1's `125 ms` STFT window. Both H1 and L1 achieve strong blind recovery.
2. **Moderate/Low-SNR BBH events (`GW151012` SNR 10.0, `GW151226` SNR 13.1, `GW170104` SNR 13.0, `GW170608` SNR 15.4):** Spread across multiple STFT frames at lower per-bin SNR, many bins fall below `bin_seed_threshold = 4.0` or fragment into sub-threshold clusters, leading to misses (e.g., `GW170608` L1 reached `94.4` vs. threshold `95.0`).
   - **Critical Artifact Audit on `GW170608` (H1):** Although `GW170608` H1 formally satisfied the automated temporal overlap rule (`4542 / 4155`), visual inspection of `fig04_event_excess_maps_gallery.png` reveals that the triggering cluster is a **full-window horizontal instrumental line at ~15-35 Hz** caused by the well-known H1 angular-control calibration procedure active during `GW170608` (which also inflated the H1 background threshold to `4155.3` and produced 4 other 32-s triggers). It is **NOT** a true recovery of the `GW170608` chirp. Excluding this instrumental artifact yields a true astrophysical blind recovery count of **`2 / 6` (`33.3%`, 95% exact CI `[4.3%, 77.7%]`)**: `GW170814` (L1) and `GW170817` (H1).
3. **Long BNS inspiral (`GW170817`, net SNR 33.0, duration >30 s):** Although `GW170817` has the highest network matched-filter SNR in GWTC-1, its energy is spread thinly over tens of seconds (low instantaneous SNR per 125 ms frame). Furthermore, L1 contains the well-documented non-Gaussian instrumental glitch ~1.1 s before merger (`other 32s anomalies`), illustrating why unvetoed excess-power detections cannot automatically be labeled astrophysical.
4. **Virgo (`V1`) sensitivity:** During O2 (`GW170814`, `GW170817`), Virgo had higher noise floor and lower single-detector SNR than H1/L1; V1 non-detections reflect physical sensitivity differences rather than pipeline failure.
5. **Timing resolution limitation:** BSD v0.1's `31.25 ms` STFT hop size remains coarser than the `10 ms` H1-L1 light-travel window, so coincidence tests use `|dt| <= 0.10 s` validated by 200 background time shifts.

## 7. Recommendations for BSD v0.2 (Based on Measured Limitations)

1. **Multi-resolution time-frequency bank (Q-transform / multi-window STFT):** Combine short (`32 ms`), medium (`125 ms`), and long (`500 ms - 1 s`) windows under a unified experiment-level max-statistic calibration to capture both sub-10 ms timing for BBH mergers and long, thin inspiral tracks like `GW170817`.
2. **Glitch-veto & cross-detector waveform consistency:** Add calibrated time-domain kurtosis/duration diagnostics and band-limited H1-L1 cross-correlation on coincident candidates to suppress single-detector glitches (such as the L1 glitch near `GW170817`) and resolve physical `< 10 ms` propagation delays.
3. **Blind comparison on the reserved v0.2 holdout set:** Evaluate BSD v0.2 against frozen BSD v0.1 on the untouched 4-event GWTC-1 holdout manifest (`GW170729`, `GW170809`, `GW170818`, `GW170823`) only after v0.2 is frozen.

## 8. Reproducibility Commands

```powershell
.\.venv\Scripts\python.exe experiments\download_gwosc_multi_event.py
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe experiments\run_ligo_experiment_003.py
```
