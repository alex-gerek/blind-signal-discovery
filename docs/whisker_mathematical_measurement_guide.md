# WHISKER / Blind Signal Discovery (BSD) — Mathematical Measurement Guide

**Scope:** Python Local Measurement Workstation, software series `0.3.0`, with calibration and adaptive fixes in `986dd46`. Documentation corrected on 2026-11-26 against `e1bd520`. This is a source-based explanation, not a new scientific validation run.

The [original guide](whisker_mathematical_measurement_guide_historical_40c7bca.md) is preserved as a historical record. Its 120-background-window, FAP=1% detections used the pre-audit fallback and must not be interpreted as valid 1% conformal detections.

The [browser demo](workstation/index.html) is a separate JavaScript implementation: one STFT scale, 512 Hz synthetic sampling and 2-second synthetic windows, different random streams and noise generators. The Python formulas and recorded numerical examples below do not establish browser-demo equivalence. See [known limitations](CURRENT_LIMITATIONS.md) before interpreting file measurements or exported results.

## 1. Quantities, coordinates, and units

| Quantity | Definition and interpretation | Source |
| --- | --- | --- |
| Waveform | `Observation.signal`, a 1D time series. Synthetic units are dimensionless; integer PCM WAV is normalized; GWOSC strain is scaled by `1e18`; CSV/NPY units are user-declared. | `src/bsd/io_loader.py`, `src/bsd/workstation_runner.py` |
| Within-window time | `t[n] = n / fs`. Canonical Python analysis uses `fs = 1024 Hz` and non-overlapping 4-second windows (4096 samples). | `slice_into_observations` |
| Exported timeline time | Window `k` starts at `k * 4 s` relative to target start. Fields labelled absolute time in the UI are relative to that target origin; do not assume they contain a GPS or UTC epoch. | `_build_rows` |
| Frequency | STFT frequency grid `f[m] = m * fs / nperseg`, in Hz. Detector maps exclude DC and Nyquist. | `src/bsd/detector.py` |
| Display spectrogram | `10 log10(max(abs(Z_stft)**2, 1e-24))`, relative to one squared input unit. The Python display uses Hann `nperseg=128`, `noverlap=96`, no boundary extension or padding. | `_plot_window_waveform_and_spec` |
| Root-energy SNR | `norm(s, 2) / sigma`, a linear dimensionless signal norm relative to the generator's noise scale. It is not dB or a power ratio. | `_scale_to_energy_snr` |
| Raw cluster score | Sum of `(excess_power - 1)` over a qualifying 8-connected component. | `TemplateFreeDetector.analyze` |
| Global STFT score | Maximum across scale-standardized raw window scores. It and the calibrated STFT threshold have the same dimensionless scale. | `MultiResolutionDetector` |
| Empirical p-value | `(1 + number of null scores >= score) / (N_B + 1)`. | `MultiResolutionCalibration.empirical_p_value` |
| Observed false-alarm fraction | False positives divided by the stated set of noise-only windows. Availability and exclusions must be reported with the denominator. | `run_synthetic_self_test` |

The display is squared spectral amplitude using SciPy's spectrum scaling, not a power spectral density in units squared per Hz, and not acoustic dB SPL. Its `Z_stft` is a Fourier amplitude, distinct from the detector's normalized excess-power map.

### Noise scale and injected SNR

`generate_gaussian_white` draws Gaussian samples with the specified population standard deviation; its individual windows are not forced to have exactly zero sample mean or exact RMS. Other generators, including `pink_colored`, use `_normalize_rms`; consult the specific generator before asserting exact sample normalization.

Python `pink_colored` uses the amplitude shaping factor

$$H(f) = (f^2 + f_{\mathrm{knee}}^2)^{-\alpha/4},$$

with `f_knee = 16 Hz` and per-window `alpha` uniform on `[0.85, 1.15]`, followed by RMS normalization. Browser pink noise currently uses a different spectrum; equal model names and seeds do not give equal data.

Injected signals are scaled as

$$s_{\mathrm{scaled}} = s_{\mathrm{clean}}\frac{\mathrm{SNR}_E\,\sigma}{\|s_{\mathrm{clean}}\|_2}.$$

For a nonzero waveform, increasing root-energy SNR from 12 to 18 multiplies signal amplitude by 1.5 and signal energy by 2.25 (+3.52 dB). In ideal white noise with known variance this norm agrees with the optimal discrete matched-filter SNR. In colored noise it is an unweighted time-domain norm and generally differs from noise-weighted matched-filter SNR. Neither the noise peak amplitude nor a spectrogram score is the injected SNR.

## 2. Python STFT statistic

The default bank uses three Hann STFT scales at 1024 Hz:

| Scale | `nperseg` / `noverlap` | Frequency spacing | Hop |
| --- | --- | --- | --- |
| `short_31ms` | 32 / 28 | 32 Hz | 3.90625 ms |
| `medium_125ms` | 128 / 96 | 8 Hz | 31.25 ms |
| `long_500ms` | 512 / 384 | 2 Hz | 125 ms |

Each scale removes the window median, computes unpadded STFT power, estimates a frequency baseline from temporal medians divided by `ln(2)`, smooths that baseline, and divides by a similarly estimated temporal envelope. The median-to-mean correction is motivated by the exponential power model for Gaussian noise. Estimated baselines, finite windows, overlapping frames and colored noise mean the final map is not guaranteed to consist of independent `Exp(1)` samples.

Bins above the seed threshold (default 4.0) form 8-connected components; components with at least two bins receive raw score

$$A(C) = \sum_{(f,t)\in C}(Z[f,t]-1).$$

On chronological Split A, each scale receives median `m_r` and spread `s_r = max(q95_r - m_r, 1)`. Scores are standardized as `(A - m_r) / s_r`. Split B contains the resulting global maxima across scales. With default `norm_fraction=0.30`, `n_A=max(20, round(0.30*n_total))` and `N_B=n_total-n_A`.

Persistent broadband search uses a frozen-reference sum of squared positive standardized log-band deviations. That construction does not establish a chi-square null distribution. Its threshold is calibrated from background scores. Single-channel broadband excess remains ambiguous with background drift.

## 3. Exact conformal calibration and decisions

For `N_B` finite null scores and finite `0 < alpha < 1`, the implementation finds the largest integer `m` satisfying

$$\frac{m}{N_B+1}\le\alpha.$$

It checks this condition with the same floating-point division used for empirical p-values, without an upward tolerance. When `m >= 1`, the one-based rank is `k = N_B + 1 - m` and the threshold is the `k`-th sorted null score. A detection requires a strict threshold crossing. Ties at the threshold are conservatively rejected. A permitted crossing must have `p_emp <= alpha`; a large score cannot bypass an unattainable target FAP.

If `1/(N_B+1) > alpha`, no finite order-statistic threshold can support the requested empirical p-value. Standalone calibration rejects this by default; explicit unavailable-calibration paths use `None`/JSON `null` and abstain. The Python Workstation rejects an unattainable initial calibration with an explanatory input error. Adaptive managers can later become unavailable and recover as enough background accumulates. Per-window availability and the threshold actually used must accompany every decision.

| Total background windows | Split A | Split B | Minimum empirical p | FAP=1% attainable? |
| ---: | ---: | ---: | ---: | --- |
| 120 | 36 | 84 | `1/85 = 0.0117647...` | No |
| 200 | 60 | 140 | `1/141 = 0.0070922...` | Yes |
| 220 | 66 | 154 | `1/155 = 0.0064516...` | Yes |

At FAP=1%, at least 99 Split-B scores are required. Do not clip an unavailable rank to the sample maximum. The former 120-window example did so and is corrected in the [erratum](AUDIT_ERRATUM_47a5475.md).

The legacy `p_value_is_upper_bound=True` flag means the score exceeds every observed null score and the returned p-value is at the finite-sample resolution floor. It does not prove that the true population tail probability is below that floor. Population FAR remains uncertain after calibration.

Conformal false-alarm statements rely on exchangeability of calibration and future null scores (conditional on the separately fitted score normalization). They are marginal over repeated calibration and evaluation samples; they do not guarantee the same conditional FAR for every frozen calibration realization or establish stationarity on arbitrary sensors.

## 4. Availability, uncertainty, and validation labels

`PROVISIONAL` and `VALIDATED` are internal Python calibration screening labels. `VALIDATED` requires all of: `N_B >= 150`, estimated `N_eff >= 120`, Ljung-Box linear-lag p-value at least 0.05, and dependence-adjusted upper 90% FAR reference interval no larger than 2.5%. `INSUFFICIENT` includes unattainable FAP or too little background/effective sample size. Initial unattainable input is rejected before a normal Workstation run.

`N_eff` uses Bartlett-weighted positive lag correlations and is an estimate, not a measured count of truly independent observations. A nonsignificant Ljung-Box result does not prove independence or exchangeability. For continuous i.i.d. null scores, the tail probability above a random `k`-th calibration order statistic follows `Beta(N_B+1-k, k)` across repeated calibration samples. Ties and model violations require qualifications. Replacing `N_B` with estimated `N_eff` is a heuristic dependence adjustment; its nominal coverage is not guaranteed for arbitrary dependent sensors.

`CALIBRATED_STATIONARY` describes the regime monitor's state relative to its active reference. It is not independent proof of physical stationarity. Dedicated-channel availability, Union availability, diagnostic nominations, and accepted active decisions are separate fields. The Sidak Union budget assumes independence or suitable positive orthant dependence of channels; it is not an unconditional arbitrary-dependence guarantee.

A non-detection with an available detector is a miss for that injection and configuration. An unavailable decision is abstention. Current self-test summaries can include unavailable windows in FAR/efficiency denominators; consult availability fields and [known limitations](CURRENT_LIMITATIONS.md) before interpreting those summaries.

## 5. Recorded Python example: 220 background windows

The [manual verification record](manual_workstation_synthetic_verification_986dd46.md) documents observations and an in-memory control check at `986dd46`: `transient`, `causal_trailing`, `target_fap=0.01`, calibration seed 0, self-test seed 2026, 220 background windows, 80 noise-only target windows, 10 pink-noise linear-chirp injections, vetoes off.

The recorded threshold is approximately 2.888922, with 154 Split-B windows, estimated `N_eff=99.73` and `PROVISIONAL` calibration. Union is unavailable. SNR 18 recovered 10/10 injections; SNR 8 recovered 3/10; both recorded 1/80 noise-only false alarms. At SNR 8, window 4 was a noise false alarm, window 81 a missed injection, and window 83 a detected injection.

These are previously recorded single-seed observations, not newly generated results. No UI export bundle was saved from that manual session. They do not establish general sensitivity or real-data FAR, and cannot serve as browser-demo reference values.

## 6. Plots, exports, and provenance

The waveform shows input amplitude; the display spectrogram shows relative STFT power; the timeline compares scores with the threshold used. Candidate boxes are diagnostic structures and may exist while an accepted decision is unavailable. A box is not a physical event or proof of novelty.

Use `active_mode_anomaly_flag` for the recorded accepted decision and retain `detection_available`, warnings, Union status and provenance. Current exports round some numeric fields, so near a boundary the rounded score and threshold can be equal even when the original values strictly crossed. Strict JSON syntax and `null` thresholds do not by themselves make an export sufficient for exact replay. See [known limitations](CURRENT_LIMITATIONS.md).

For each future check, save the implementation revision, configuration, seeds, input hashes, preprocessing, calibration evidence and full availability accounting in a fresh run location. Preserve earlier results. The [historical guide](whisker_mathematical_measurement_guide_historical_40c7bca.md), [erratum](AUDIT_ERRATUM_47a5475.md), and [research history](RESEARCH_HISTORY.md) explain which previous observations remain historical.
