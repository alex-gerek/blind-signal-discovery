# Blind Signal Discovery (BSD) — Experiment 002: Real LIGO Data Validation (GW150914)

> **Scientific Scope Statement:** This is an external methodological validation of the frozen BSD v0.1 detector on authentic public LIGO O1 strain data from GWOSC. **A statistically unusual event is not automatically a gravitational wave.** Instrumental glitches remain an alternative hypothesis, and astrophysical confirmation belongs strictly to the LIGO Scientific Collaboration and Virgo Collaboration.

## 1. Objective

Evaluate whether the frozen BSD v0.1 template-free anomaly detector (`src/bsd/detector.py`, unchanged from `main`) can identify statistically unusual events in real LIGO Hanford (H1) and Livingston (L1) strain measurements — including the independently confirmed GW150914 event — without using any signal template or the published event timestamp during calibration or detection.

## 2. Official Data Provenance & Dataset Details

All data were downloaded unchanged from the official Gravitational Wave Open Science Center (GWOSC) GWTC-1-confident GW150914 v3 release (`https://gwosc.org/events/GW150914/`). Raw files are stored under `data/gwosc/` on drive D: (excluded from Git via `.gitignore`).

| Detector | Duration (s) | Sample Rate (Hz) | GPS Interval | Size (MB) | All Finite? | DQ DATA (%) | SHA-256 (prefix) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| H1 | 32 | 4096 | `[1126259447, 1126259479]` | 1.04 | True | 100.0% | `c5ea87beced5094b...` |
| L1 | 32 | 4096 | `[1126259447, 1126259479]` | 1.01 | True | 100.0% | `9d65bde82f7ec035...` |
| H1 | 4096 | 4096 | `[1126257415, 1126261511]` | 130.16 | True | 100.0% | `8ef8b762cff53b42...` |
| L1 | 4096 | 4096 | `[1126257415, 1126261511]` | 125.87 | True | 100.0% | `bdd2dca270ad5205...` |

**Why 4096-second datasets were required in addition to the 32-second files:**
The 32-second files contain at most 7 non-overlapping 4-second windows (with 1-second padding), which cannot statistically calibrate a 1% experiment-level false-alarm probability (minimum resolvable empirical p-value would be $1/8 = 12.5\%$). Using the official 4096-second GWOSC release for the same observation epoch provides $>900$ independent, DQ-verified 4-second background segments per detector after strictly excluding a $\pm 64\text{ s}$ guard zone around the entire 32-second target file (`GPS [1126259383, 1126259543]`).

## 3. Preprocessing & Preservation of Frozen BSD v0.1

The BSD v0.1 detector (`TemplateFreeDetector` with `nperseg=128`, `noverlap=96`, `window='hann'`, `freq_median_kernel=11`, `time_median_kernel=15`, `bin_seed_threshold=4.0`, `min_cluster_bins=2`) was **not modified**. Because raw LIGO strain has amplitude $h \sim 10^{-19}$ (which would underflow `_EPS = 1e-18` in squared STFT power), a separate data-adaptation layer (`src/bsd/ligo_adapter.py`) was implemented with two explicitly distinguished pipelines:

1. **Pipeline A (`frozen_v01_minimal`):** Median subtraction, division by background MAD standard deviation (bringing strain to $O(1)$ without altering spectral shape), and zero-phase anti-aliased decimation `4096 Hz -> 1024 Hz` (`scipy.signal.decimate`, IIR).
2. **Pipeline B (`adapted_ligo_whitened`):** `1e18` numerical scaling, frequency-domain whitening using a robust `average='median'` Welch PSD (`4.0 s` Hann windows) estimated strictly on independent background segments, zero-phase 4th-order Butterworth bandpass `30-380 Hz` (standard aLIGO O1 sensitive band per GWOSC tutorial), and zero-phase decimation `4096 Hz -> 1024 Hz`. All filtering uses `1.0 s` left/right context padding trimmed after filtering to prevent edge transients.

## 4. Noise Calibration & Independent Background Validation

For each detector and pipeline, clean 4-second segments outside `GPS [1126259383, 1126259543]` were partitioned into three strictly disjoint sets:
- **PSD / scale estimation:** `60` segments (`240 s`)
- **Threshold calibration (target FAP = 1.0%):** `500` segments (`2000 s`) via `calibrate_detector`
- **Independent background FAR validation:** `422` segments (H1) and `422` segments (L1)

| Pipeline | Detector | Calibrated 1% Threshold | Independent Bg FP / N | Empirical Bg FAR | 95% Exact CI | Statistically Credible? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `frozen_v01_minimal` | H1 | `9214.09` | 2 / 422 | 0.47% | [0.06%, 1.70%] | Yes |
| `frozen_v01_minimal` | L1 | `20567.77` | 10 / 422 | 2.37% | [1.14%, 4.31%] | NO (FAR inflated) |
| `adapted_ligo_whitened` | H1 | `77.10` | 3 / 422 | 0.71% | [0.15%, 2.06%] | Yes |
| `adapted_ligo_whitened` | L1 | `76.33` | 5 / 422 | 1.18% | [0.39%, 2.74%] | Yes |

## 5. Blind Single-Detector Results (32-Second Target Segment)

All candidates below were extracted **without** knowledge of the GW150914 timestamp across the uniform 4-second window grid (`GPS 1126259448` to `1126259476`).

### 5.1 Pipeline A (`frozen_v01_minimal` — Scale + Decimate Only)

*No candidates exceeded the calibrated 1% FAP threshold.*

### 5.2 Pipeline B (`adapted_ligo_whitened` — Background Whitening + 30-380 Hz Bandpass)

| Detector | GPS Peak (s) | Interval [Start, End] (s) | Duration (s) | Freq Range (Hz) | Peak Freq (Hz) | Score | Emp. p-value |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| H1 | `1126259462.375` | `[1126259462.141, 1126259462.453]` | `0.312` | `[20.0, 284.0]` | `56.0` | `818.0` | `<= 0.00200` |
| L1 | `1126259462.344` | `[1126259462.297, 1126259462.453]` | `0.156` | `[28.0, 116.0]` | `40.0` | `136.6` | `0.00399` |
| L1 | `1126259462.406` | `[1126259462.328, 1126259462.453]` | `0.125` | `[124.0, 188.0]` | `176.0` | `99.5` | `0.00399` |

## 6. Two-Detector Comparison (H1 vs. L1) & Time-Shift Controls

The maximum physical gravitational-wave light-travel delay between Hanford (H1) and Livingston (L1) is approximately `10 ms` (`0.010 s`). However, the frozen BSD v0.1 STFT has a frame step of `dt = 31.25 ms` (`nperseg=128, noverlap=96` at `1024 Hz`). **Because `dt = 31.25 ms > 10 ms`, BSD v0.1 cannot resolve sub-10 ms inter-site propagation delays directly.** We therefore test temporal coincidence within `|dt| <= 0.10 s` (~3 STFT frames) and calibrate accidental coincidence rates using 200 unphysical time shifts on the independent background validation set.

- **Pipeline A (`frozen_v01_minimal`):** `0` target H1-L1 coincidences within `0.10 s`; background time-shift rate = `0.1400` coincidences/shift (`28` across `200` shifts).
- **Pipeline B (`adapted_ligo_whitened`):** `2` target H1-L1 coincidences within `0.10 s`; background time-shift rate = `0.0100` coincidences/shift (`2` across `200` shifts).

| H1 GPS Peak | L1 GPS Peak | dt (H1 - L1, s) | Interval Gap (s) | H1 Score | L1 Score | Within 10 ms? |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `1126259462.3750` | `1126259462.3438` | `+0.0312` | `0.0000` | `818.0` | `136.6` | No (STFT-limited) |
| `1126259462.3750` | `1126259462.4062` | `-0.0312` | `0.0000` | `818.0` | `99.5` | No (STFT-limited) |

## 7. Independent Evaluation Against Published GW150914 Reference

Published GWOSC reference (consulted only after all candidates were frozen): `GPS = 1126259462.422`, frequency sweep `~35-250 Hz`.

| Pipeline | Detector | Blind Recovery of GW150914? | Timing Error (Peak - Ref, s) | Contains Ref GPS? | Freq Overlaps 35-250 Hz? |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `frozen_v01_minimal` | H1 | **NO** | `N/A` | False | False |
| `frozen_v01_minimal` | L1 | **NO** | `N/A` | False | False |
| `adapted_ligo_whitened` | H1 | **YES** | `-0.0470` | True | True |
| `adapted_ligo_whitened` | L1 | **YES** | `-0.0157` | True | True |

### Epistemic Distinction
1. **Blind anomalies:** Candidates in Section 5 were produced solely by threshold exceedances against independent LIGO background noise, with zero access to the GW150914 timestamp or waveform.
2. **Post-hoc evaluation:** Overlay of `GPS = 1126259462.422` in Figures 4, 5, and 7 was applied strictly after detection to measure localization accuracy.
3. **Astrophysical confirmation:** Independent astrophysical identification of GW150914 as a binary black hole merger relies on the LIGO/Virgo collaboration's matched-filter, waveform-consistency, and detector-characterization analyses, not on this single-resolution excess-power prototype.

## 8. Limitations & Failure Analysis

1. **Failure of unadapted `frozen_v01_minimal`:** Without external LIGO PSD whitening and `30-380 Hz` bandpass filtering, seismic noise below `25 Hz` (6-8 orders of magnitude stronger than the sensitive band) leaks across the `125 ms` Hann STFT window and swamps the astrophysical signal.
2. **Real LIGO non-Gaussian glitches:** Real interferometer background contains instrumental glitches (heavy tails), causing finite-sample threshold variance and potential FAR drift over hour-long timescales compared to stationary synthetic Gaussian noise.
3. **STFT timing resolution (`31.25 ms`):** Coarser than the `10 ms` H1-L1 light-travel bound, preventing sub-10 ms propagation-delay verification without finer time resolution or cross-correlation.
4. **Finite calibration tail (`N_cal = 500`):** Minimum resolvable empirical p-value is `1 / 501 = 0.00200`; scores exceeding the maximum background sample are strictly reported as upper bounds (`<= 0.00200`).

## 9. Generated Figures

- `figures/fig01_raw_strain_h1_l1.png`
- `figures/fig02_processed_strain_h1_l1.png`
- `figures/fig03_spectrograms_and_excess_maps.png`
- `figures/fig04_anomaly_score_timelines.png`
- `figures/fig05_candidate_overlays_eval.png`
- `figures/fig06_calibration_and_background_far.png`
- `figures/fig07_cross_detector_timing.png`

## 10. Reproducibility Commands

```powershell
.\.venv\Scripts\python.exe experiments\download_gwosc_data.py
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe experiments\run_ligo_experiment_002.py
```

## 11. Recommendations for Experiment 003

1. **Multi-resolution Q-transform / finer time step:** Reduce STFT hop size to `<= 2 ms` (or add band-limited H1-L1 cross-correlation on coincident candidates) to resolve the `10 ms` physical propagation window.
2. **Glitch-veto diagnostics:** Incorporate GWOSC `CBC_CAT2`/`CBC_CAT3` hardware-injection and data-quality flags alongside a calibrated kurtosis/duration veto to separate short instrumental blips from chirping transients.
3. **Multi-event blind validation:** Evaluate the same frozen adapted pipeline on additional GWTC-1 events (e.g., GW151226, GW170104, GW170814) and null off-source stretches to measure empirical sensitivity across real SNR regimes.
