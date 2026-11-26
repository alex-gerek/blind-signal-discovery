# WHISKER / BSD Exploratory Real-Data Search: GW170814 Epoch (H1, L1, V1)

> **Scope Notice:** This is an **exploratory analysis of a previously used O2 GWOSC epoch** (`GW170814`, `GPS 1186739814..1186743910`), **NOT a blind discovery search** on unseen data. Detection of `GW170814` itself is a repeatability / positive-control check; the adjacent off-source interval is used for exploratory candidate ranking, reversible DQ/veto triage, and 3-detector cross-checks. The reserved O3 holdout and detector mathematics remain completely untouched.

## 1. Verified Non-Overlapping Epoch Partition (Including `pad_s = 1.0 s`)

Every 4.0-second window uses `pad_s = 1.0 s` context on both sides (`6.0 s` padded span) during whitening and bandpass filtering. All four data blocks were verified programmatically to have strictly positive gaps between their padded footprints:

| Block | Core GPS Interval | Core Duration | 4-s Windows | Padded Footprint | Guard Gap After Padded |
|---|---|---:|---:|---|---:|
| `PSD_FIT` | `1186739816..1186740056` | `240 s` | `60` | `1186739815..1186740057` | `22 s` |
| `CAL_BLOCK_A` | `1186740080..1186741360` | `1280 s` | `320` | `1186740079..1186741361` | `178 s` |
| `SEARCH_TARGET` | `1186741540..1186742180` | `640 s` | `160` | `1186741539..1186742181` | `178 s` |
| `CAL_BLOCK_B` | `1186742360..1186743908` | `1548 s` | `387` | `1186742359..1186743909` | `0 s` |

### Exact Window Selection Around `GW170814` (`t_c = 1186741861.5`)

- **Target search block:** `160` contiguous non-overlapping 4-s windows spanning `GPS 1186741540..1186742180` (`640.0 s`).
- **On-source positive control (1 window):** Window `80` (`GPS 1186741860..1186741864`), containing `GW170814` at `t = 1.5 s`.
- **`+-16 s` event exclusion rule:** Exclusion interval `[1186741845.5, 1186741877.5]`. Any window whose **padded** footprint `[gps0 - 1, gps0 + 5]` intersects this interval is excluded from off-source (`8` guard-zone windows: Windows `76..79` [`1186741844..1186741860`] and Windows `81..84` [`1186741864..1186741880`]).
- **Resulting off-source search set:** **`151` non-overlapping 4-s windows** (`604.0 s = 10.07 min`; Windows `0..75` pre-event and Windows `85..159` post-event).
- **Unusable / DQ-dropped windows across all intervals (H1/L1/V1):** `0` (complete exclusion log saved in `excluded_and_unusable_windows.csv`).

## 2. Frozen Preprocessing, Detector, Calibration, and Veto Configuration

- **Preprocessing (`src/bsd/ligo_adapter.py`):** `adapted_ligo_whitened` (`1e18` scaling, median-Welch PSD fitted strictly on the `60 x 4 s = 240 s` PSD block, zero-phase 4th-order Butterworth bandpass `30-380 Hz`, anti-aliased decimation `4096 Hz -> 1024 Hz`).
- **Detector (`src/bsd/detector_v2.py`):** Unmodified 3-scale `MultiResolutionDetector` (`short_31ms`: hop `3.91 ms`; `medium_125ms`: hop `31.25 ms`; `long_500ms`: hop `125.0 ms` at `1024 Hz`), `target_fap = 0.01` (nominal calibration quantile parameter — **NOT** a measured FAR on real LIGO/Virgo data), `norm_fraction = 0.30` (`Split A = 212`, `Split B = 495`, `p_min = 1/496 = 0.002016`).
- **Reversible vetoes (`src/bsd/vetoes.py`):** `VetoConfig(enable_glitch_veto=True, enable_line_veto=True)`. All raw candidates are preserved verbatim in `raw_candidates.csv` before attaching veto and GWOSC DQ/Inj annotations.

### Per-Detector Calibration Status and Search Summary

| Detector | Cal Windows (`Split A / B`) | Cal Status | `n_eff` | Cal LB `p` | Frozen Threshold | GW170814 Raw / Post-Veto (Score) | Off-Source Raw (`k/151`) | Off-Source Post-Veto (`k/151`) |
|---|---:|---|---:|---:|---:|---|---:|---:|
| `H1` | `707` (`212/495`) | `VALIDATED` | `458.6` | `0.5098` | `2.5905` | `False` / `False` (`1.14`) | `2/151` (`1.32%`, 95% Wilson `[0.36%, 4.70%]`) | `1/151` (`0.66%`) |
| `L1` | `707` (`212/495`) | `VALIDATED` | `432.5` | `0.3014` | `2.2401` | `True` / `True` (`8.89`) | `3/151` (`1.99%`, 95% Wilson `[0.68%, 5.68%]`) | `3/151` (`1.99%`) |
| `V1` | `707` (`212/495`) | `VALIDATED` | `473.6` | `1.0000` | `4.0069` | `False` / `False` (`3.62`) | `2/151` (`1.32%`, 95% Wilson `[0.36%, 4.70%]`) | `2/151` (`1.32%`) |

> **Real-data FAR note:** The nominal `target_fap = 0.01` is only the Split-B calibration quantile rule; it is **not** a measured false-alarm rate on real non-stationary LIGO/Virgo noise. Wilson intervals on real-data off-source counts are descriptive approximations assuming independent windows.

## 3. Peak-Time Localization Uncertainty, Coincidence Screening, and Time-Slide Resolution

- **STFT localization uncertainty:** At `1024 Hz`, half-hop time quantization is `+/-1.95 ms` (`short_31ms`), `+/-15.63 ms` (`medium_125ms`), and `+/-62.50 ms` (`long_500ms`). Combining STFT peak quantization with physical light-travel bounds (`10.0 ms` H1-L1, `27.3 ms` H1-V1, `26.4 ms` L1-V1), prespecified coarse STFT screening uses `max_dt = 50 ms` (H1-L1) and `65 ms` (H1-V1, L1-V1), followed by band-limited time-domain cross-correlation (`35-350 Hz`, `+-15 ms` H1-L1, `+-30 ms` LIGO-Virgo).
- **Non-detection rule:** A candidate detected in only one detector is labeled `SINGLE_DETECTOR_UNCONFIRMED` and is **never** classified as instrumental solely because another detector did not trigger (antenna-pattern factors and Virgo's lower O2 sensitivity can cause genuine signals to appear in only one or two detectors).
- **On-source coincident pairs (`GW170814`):** `0` pair(s).
- **Off-source coincident pairs:** `0` pair(s).
- **Time-slide background audit (`100` non-circular slides, `4 s` step, `+-200 s` shifts over `604 s` off-source):**
  - **Valid overlap correction:** Accounting for outer edge loss and the central `36 s` guard-zone hole (`1186741844..1186741880`), true active overlap per slide ranges from `368.0 s` to `596.0 s` (mean `468.9 s`), totaling **`13.02 hr`** (`77.6%` of the `16.78 hr` nominal `100 x 604 s` figure computed by `run_time_slide_analysis`).
  - **Multi-cluster inflation vs. deduplicated window-level counts:** Raw cluster-level counts in `run_time_slide_analysis` (`H1-L1 = 1`, `H1-V1 = 0`, `L1-V1 = 9`) are inflated because multiple STFT-scale clusters from the same physical window (`2` H1, `4` L1, `8` V1 clusters from only `2`, `3`, `2` unique windows) are cross-matched pairwise and `match_coincident_candidates` also triggers on `interval_gap <= max_dt`. When deduplicated to 1 peak per window with `|dt_peak| <= max_dt` inside valid active overlap, window-level slide coincidences are: `H1-L1 = 0`, `H1-V1 = 0`, `L1-V1 = 0`.
  - **4-second shift step alignment:** Shifting by exact `4.0 s` multiples aligns 4-s window boundaries and the `36 s` (`9 x 4 s`) guard-zone hole identically across slides; any hypothetical 4-s periodic windowing artifact would remain phase-locked.
  - **Resolution caveat:** With only `604 s` (`10.07 min`) of off-source data (`13.02 hr` valid slide livetime) and only `2-3` triggered windows per detector, time-slides do **not** support a meaningful astrophysical coincidence significance claim.

## 4. Top Ranked Off-Source Raw Candidates (All Raw Candidates Preserved in `raw_candidates.csv`)

> **Empirical `p`-value definition:** `p_emp = (1 + n_ge) / (N_B + 1)` with `N_B = 495` Split-B calibration windows. When `n_ge = 0`, `p_emp = 1/496 = 0.00202` (calibration resolution floor). This is a conformal finite-sample point estimate under exact i.i.d. exchangeability — it is **not** written as `p <= 0.0020` (since `1/496 > 0.0020` and `0/495` does not prove the true tail probability is below `0.0020`), and it is **not** an astrophysical probability or a measured real-data false-alarm rate.

| Rank | Det | GPS Peak | `t - t_c` (s) | Scale | Peak Freq (Hz) | Std Score | Emp. `p` (`N_B=495`) | Veto Status (Reasons) | GWOSC DQ / Inj | Coincidence |
|---:|---|---:|---:|---|---:|---:|---:|---|---|---|
| 1 | `V1` | `1186741602.344` | `-259.16` | `medium_125ms` | `32.0` | `457.21` | `0.00202 (0/495; floor 1/496)` | `VETOED_ANNOTATED` (`IMPULSIVE_GLITCH_EXTREME_KURTOSIS`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 2 | `V1` | `1186741602.375` | `-259.12` | `long_500ms` | `112.0` | `221.96` | `0.00202 (0/495; floor 1/496)` | `VETOED_ANNOTATED` (`IMPULSIVE_GLITCH_EXTREME_KURTOSIS`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 3 | `V1` | `1186741602.375` | `-259.12` | `long_500ms` | `28.0` | `147.87` | `0.00202 (0/495; floor 1/496)` | `VETOED_ANNOTATED` (`IMPULSIVE_GLITCH_EXTREME_KURTOSIS`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 4 | `V1` | `1186741602.352` | `-259.15` | `short_31ms` | `128.0` | `26.62` | `0.00403` | `SURVIVING` (`NONE`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 5 | `V1` | `1186741549.344` | `-312.16` | `medium_125ms` | `24.0` | `20.20` | `0.00403` | `SURVIVING` (`NONE`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 6 | `V1` | `1186741549.348` | `-312.15` | `short_31ms` | `32.0` | `8.93` | `0.00403` | `SURVIVING` (`NONE`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 7 | `V1` | `1186741549.375` | `-312.12` | `long_500ms` | `30.0` | `7.79` | `0.00403` | `SURVIVING` (`NONE`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 8 | `H1` | `1186742016.941` | `+155.44` | `short_31ms` | `352.0` | `5.46` | `0.00202 (0/495; floor 1/496)` | `VETOED_ANNOTATED` (`IMPULSIVE_GLITCH_ULTRA_SHORT_SPIKE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 9 | `V1` | `1186741602.375` | `-259.12` | `long_500ms` | `368.0` | `4.76` | `0.00806` | `VETOED_ANNOTATED` (`IMPULSIVE_GLITCH_EXTREME_KURTOSIS`) | `PASS` / `NONE` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 10 | `L1` | `1186741743.000` | `-118.50` | `long_500ms` | `36.0` | `3.18` | `0.00403` | `SURVIVING` (`NONE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 11 | `H1` | `1186741812.500` | `-49.00` | `long_500ms` | `440.0` | `3.05` | `0.00605` | `SURVIVING` (`NONE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 12 | `L1` | `1186741742.918` | `-118.58` | `short_31ms` | `32.0` | `2.50` | `0.00806` | `SURVIVING` (`NONE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 13 | `L1` | `1186742170.125` | `+308.62` | `long_500ms` | `502.0` | `2.48` | `0.00806` | `SURVIVING` (`NONE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |
| 14 | `L1` | `1186741712.875` | `-148.62` | `long_500ms` | `44.0` | `2.38` | `0.00806` | `SURVIVING` (`NONE`) | `PASS` / `CW_HW_INJ` | `SINGLE_DETECTOR_UNCONFIRMED` |

## 5. Focused Audit Findings

### 5.1 Why `GW170814` Exceeded Threshold in `L1` Only (Comparison with Prior Benchmarks)
- **Consistency with prior repository benchmarks:** Both `Exp 003` (v0.1 `TemplateFreeDetector`: `H1` max score `37.93` vs. threshold `76.80`, ratio `0.494`) and `v02_phase3` (v0.2 `MultiResolutionDetector`: `0` H1-L1 coincident pairs pre- and post-veto) likewise did not exceed the single-detector 1% threshold in `H1` on `GW170814`.
- **Per-scale scores (`Window 80`, `GPS 1186741860..1186741864`):** `L1` exceeded threshold `2.24` across all 3 scales (max `8.89`), whereas `H1` (`max = 1.14 < 2.59`) and `V1` (`max = 3.62 < 4.01`) remained below their frozen single-detector thresholds without matched-filter templates.

### 5.2 Forensic Audit of Top Off-Source Candidate (`V1` at `GPS 1186741602.344`, `2017-08-14 10:26:24.344 UTC`)
- **Amplitude definition and units:** On the 4-second whitened, `30-380 Hz` bandpassed, `1024 Hz` decimated series (`GPS 1186741600..1186741604`, dimensionless whitened-strain units), the median-centered peak is `32.445` with `raw MAD = 0.1790` (`peak / raw MAD = 181.22`), Gaussian-consistent MAD scale `1.4826 x MAD = 0.2654` (`peak / [1.4826 x MAD] = 122.23`, reported by `vetoes.py` as `local_peak_mad_sigma`), and full-window sample SD `1.2300` (`peak / window SD = 26.38`). This ratio is **not** a Gaussian standard-deviation significance (`sigma`) on non-Gaussian Virgo noise.
- **Why cluster metadata span `16-496 Hz` (`short_31ms`) and `4-508 Hz` (`medium_125ms`) despite `30-380 Hz` bandpass:** (1) The 4th-order zero-phase Butterworth bandpass has a smooth roll-off, attenuating stopband power by `10^-5` to `10^-9` relative to the passband; (2) per-frequency median normalization `Z[f, t] = P[f, t] / (S_bg(f) * T_bg(t))` cancels out the smooth filter gain `|H(f)|^2`, allowing heavily attenuated stopband bins to exceed `Z > 4`; (3) cluster bounds report half-bin edges (`f_center +/- df/2`, with `df = 32 Hz` on `short_31ms` and `df = 8 Hz` on `medium_125ms`). Note also that `fig02` masks the displayed spectrogram to `30-360 Hz`. Physical signal power is concentrated inside the `30-380 Hz` passband.
- **Why `short_31ms` survived while the other 4 clusters were vetoed:** `evaluate_impulsive_glitch` computes excess kurtosis on `[t_start - 0.08 s, t_end + 0.08 s]`. Wider `medium_125ms`/`long_500ms` slices (`548-3620` samples, `kurtosis = 39.1-255.5`) exceeded `20.0`, whereas the narrower `short_31ms` slice (`256` samples) yielded `local_kurtosis = 17.23 < 20.0`, and `n_samples_ge_half_peak = 6 > 3` (bandpass ringing) missed the ultra-short spike rule. Because `apply_vetoes` sets `detected_after_veto = len(surviving_candidates) > 0`, this kept the window flagged `detected_post_veto = True` (`V1` off-source `2/151` instead of `1/151`).
- **Public catalog check, observations vs. hypotheses, and novelty caveat:**
  - *Verified observations:* Local `DQmask = 127` (all 7 bits PASS), `Injmask = 31` (all 5 `NO_*_HW_INJ` bits PASS), and GWOSC `O2_4KHZ_R1` Timeline confirms `V1_DATA`, `V1_CBC_CAT2`, `V1_BURST_CAT2/3`, and `V1_NO_BURST_HW_INJ` PASS over `[1186741500..1186741700]`. Queries of `GWOSC EventAPI` (`671` entries across GWTC-1/2.1/3/4 confident & marginal catalogs, `+-3600 s`) and `4-OGC` (`4OGC_top.txt`) returned only `GW170814` (`+259.16 s`). No public O2 Virgo glitch catalog or O2 Virgo auxiliary-channel release was found or accessible in the sources checked (`Gravity Spy` Zenodo LIGO O1-O3b releases, `GWitchHunters` O3 Virgo project, `2-OGC` URL HTTP 404, web/arXiv searches).
  - *Working hypothesis (unproven without Virgo auxiliary channels):* Internal V1 morphology (short `< 0.1 s` passband broadband impulse without chirp structure, high excess kurtosis on wider slices) is consistent with an instrumental glitch, whereas documented hardware injections are inconsistent with `Injmask = 31`. Absence in H1/L1 is recorded strictly as an observation and is not used alone to infer instrumental origin.
  - *Novelty caveat:* Lack of a public catalog match does **not** establish novelty or astrophysical origin; whether this transient was logged internally by Virgo DetChar in 2017 remains unknown.

## 6. Coincident Pairs Summary (`coincident_pairs.csv`)

No pairwise coincident candidates found.

## 7. Figures and Limitations

- `fig01_gw170814_search_timeline.png`: Standardized global score timeline for H1, L1, and V1 across the 640-s target interval.
- `fig02_gw170814_positive_control_and_top_offsource.png`: Whitened spectrograms of the GW170814 positive control (H1, L1) and top-ranked off-source candidate.
- **Total runtime:** `96.4 s`.
- **Limitations:** Previously used O2 epoch (not a blind discovery search); short `604 s` off-source livetime limits time-slide resolution; single-detector triggers remain unconfirmed without external glitch catalogs (e.g., Gravity Spy) or independent multi-detector coincidence.
