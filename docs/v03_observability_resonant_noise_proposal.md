# BSD v0.3 Observability Study Proposal: Detection Limits Under Resonant Background Noise

**Status:** Research-Only Proposal (Pre-Implementation Review — Zero Code or Parameter Changes)
**Branch:** `development/v0.3` (following pushed commits `8ed1813` and `edd1ec9`)
**Preserved Baselines:** `v0.1`, `v0.2`, `v0.3 Phase 1`, `v0.3 Phase 2`, and `v0.3 Phase 2b` (`8ed1813`, `edd1ec9`) remain unmodified.
**Reserved O3 Holdout (`docs/v03_o3_holdout_manifest.json`):** Strictly untouched and uninspected (`0` local files).

---

## 1. Central Research Question

In the Phase 2b corrective benchmark (`edd1ec9`, `abrupt_step_500win`), staged causal recalibration into the shifted AR(2) background (`R_1`) achieved `88.0%` clean duty cycle, `0.0%` stale post-recovery glitch flags, and `0.00%` mature surviving FAR (`T_prov = 62`, `T_mature = 201`, `N_eff = 115.0 / 150`). Yet transient detection efficiency for the 10 linear chirps injected into mature `R_1` (`idx = 330..420`, injected at time-domain energy SNR `= 18.0` relative to `sigma = 1.5`) was only **`50.0%` (`5/10`)**, compared to **`100.0%` (`20/20`)** for equivalent chirps in stationary `pink_colored` noise.

**Central Question:**
> *Why does transient detection efficiency fall to `50%` for chirps after recovery into the shifted AR(2) background?*

Importantly, we **do not assume upfront** that spectral frequency masking is the primary cause. Inspection of `generate_shifted_ar2` (`src/bsd/noise.py`) shows that each window has (a) a per-window randomized AR(2) resonance `f_res ~ Uniform(110, 340) Hz` with pole radius `r ~ Uniform(0.72, 0.86)` (peak PSD gain `~13–51x`, `-3 dB` bandwidth `~46–91 Hz`) and (b) a `±25%` within-window amplitude modulation (`0.25 Hz` sinusoid). The `50%` efficiency drop could therefore arise from any combination of four distinct mechanisms that must be separated experimentally.

---

## 2. Four Competing Hypotheses to Distinguish Experimentally

| Hypothesis | Mechanism | Falsifiable Experimental Signature |
| :--- | :--- | :--- |
| **H1 — Fundamental Information Loss / Low Statistical Distinguishability** | Scaling injections by unweighted time-domain energy SNR (`_scale_to_energy_snr(clean, sigma, 18.0)`) ignores the colored PSD $S_n(f)$. A chirp overlapping the AR(2) resonance peak has a drastically lower true whitened optimal matched-filter SNR $\rho_{\text{opt}} = \sqrt{4 \int \frac{|H(f)|^2}{S_n(f)}\,df}$ than an off-resonance chirp of identical time-domain amplitude. | Efficiency collapses primarily as a function of $\rho_{\text{opt}}$ regardless of detector architecture; once $\rho_{\text{opt}}$ is held constant across frequency placements, in-resonance and off-resonance efficiency converge (or remain near the matched-filter / pre-whitened ceiling). |
| **H2 — Current STFT Detector Sensitivity Loss** | `TemplateFreeDetector` (`src/bsd/detector.py`) estimates local background via within-window 2D median filtering (`freq_median_kernel = 5/11/21`, `time_median_kernel = 31/15/7`). A sharp AR(2) resonance peak (`46–91 Hz` wide) and/or `±25%` within-window envelope modulation may cause `bg_psd` to under-whiten the resonance peak (inflating null scores and thresholds) or over-subtract slow chirps. | Even at **matched $\rho_{\text{opt}}$**, `MultiResolutionDetector` efficiency in `shifted_ar2` lags significantly behind (a) a known-PSD pre-whitened STFT control on the exact same windows and (b) `MultiResolutionDetector` on flat white noise at the same $\rho_{\text{opt}}$. |
| **H3 — Calibration & Threshold Effects** | Mature causal recalibration (`N_recal = 200`, `N_thresh = 150`, `N_eff = 115.0`) on a heterogeneous mixture of `shifted_ar2` windows (each with different `f_res` and `r`) may produce a conservative empirical `99%` threshold (`mature_surv_far = 0.00%`) relative to a large independent stationary oracle calibration (`N_thresh = 1000`). | Re-scoring the exact same `abrupt_step_500win` chirps against an independent `N = 1000` oracle `shifted_ar2` calibration (at matched `1%` empirical FAR) significantly increases detection efficiency without changing detector scores. |
| **H4 — Regime-Monitor Coverage / Gating Effects** | Stream B regime transitions, warmup exclusions, or `GLITCH_STORM` flags suppress otherwise-detected Stream A candidates. | Raw Stream A efficiency (`stft_nominated_raw`) or post-veto efficiency (`stft_nominated_post_veto`) exceeds clean-calibrated efficiency (`is_clean_calibrated`). *(Note: in `edd1ec9`, full and clean efficiency were both `50.0%`, providing a direct baseline check that H4 must still be audited quantitatively across all trials.)* |

---

## 3. Controlled Experimental Design

All experiments are **read-only observational diagnostics** using frozen `MultiResolutionDetector`, frozen `VetoConfig`, and frozen `AdaptiveRegimeManagerV2b` (`edd1ec9`). **No algorithms will be modified and no parameters will be tuned.** All data will be newly generated controlled synthetic signals and backgrounds; the reserved O3 holdout remains 100% untouched.

### Experiment 1 — Forensic Audit of the 10 Exact Post-Recovery Chirps from `edd1ec9`
Reproduce the exact `abrupt_step_500win` sequence (seed `3026`, identical to `edd1ec9`) and record a complete per-injection forensic table for all 10 post-recovery chirps (`idx = 330, 340, ..., 420`), directly comparing the **5 detected** vs. **5 missed** injections:
- Injected chirp parameters: `(f_start, f_end, f_center, duration_s, sweep_rate_hz_s)`, time-domain energy SNR (`18.0`).
- Exact window background parameters: `(f_res, pole_radius r, drift_phase)`, exact analytical AR(2) PSD $S_n(f)$ evaluated along the chirp trajectory.
- Exact theoretical whitened optimal SNR $\rho_{\text{opt}}$ and local band-averaged spectral SNR $\mathrm{SNR}_{\text{spec}}$.
- Detector internals: per-scale raw and standardized max scores (`short_31ms`, `medium_125ms`, `long_500ms`), winning scale, `stft_max_score`, active mature `R_1` threshold vs. independent `N = 1000` oracle threshold, veto status, and Stream B credibility/glitch status.
- Visual diagnostic: side-by-side time-frequency trajectory + AR(2) PSD overlays for representative detected vs. missed injections.

### Experiment 2 — Factorial Sweep Over Frequency Placement, Amplitude Control, Duration, and Sweep Rate
To avoid confounding per-window random `f_res`, fix a controlled resonant AR(2) background at `f_res = 200 Hz`, `r = 0.84` (plus a secondary randomized-`f_res` condition matching `generate_shifted_ar2`, and a **zero-modulation** AR(2) ablation removing the `±25%` `0.25 Hz` envelope to isolate envelope modulation from spectral resonance!).

Vary four controlled factors across **$N = 40$ independent Monte Carlo trials per cell** (reporting exact 90% Wilson confidence intervals):
1. **Frequency placement relative to $f_{\text{res}} = 200\text{ Hz}$:**
   - `in_resonance`: chirp centered on `[175, 225] Hz` (inside AR(2) peak)
   - `crossing_resonance`: chirp sweeping `[120, 280] Hz` (crossing the resonance)
   - `off_resonance_low`: chirp in `[50, 100] Hz` (well below resonance)
   - `off_resonance_high`: chirp in `[320, 380] Hz` (well above resonance)
2. **Paired SNR Control (Critical for separating H1 vs. H2):**
   - **Control Mode A — Fixed Time-Domain Energy SNR** (`18.0`, matching `edd1ec9`)
   - **Control Mode B — Fixed Theoretical Whitened Spectral SNR** $\rho_{\text{opt}} \in \{10, 14, 18, 22, 26\}$ (scaling each chirp so $4 \int |H(f)|^2 / S_n(f)\,df = \rho_{\text{opt}}^2$ exactly!)
3. **Duration & Sweep Rate Grid** (evaluated at matched $\rho_{\text{opt}} = 18$):
   - Short fast chirp (`duration = 0.25 s`, `bandwidth = 120 Hz`, `sweep_rate = 480 Hz/s`)
   - Medium chirp (`duration = 1.0 s`, `bandwidth = 80 Hz`, `sweep_rate = 80 Hz/s`)
   - Long slow chirp (`duration = 2.5 s`, `bandwidth = 40 Hz`, `sweep_rate = 16 Hz/s`)

### Experiment 3 — Independent Controls & Four-Way Loss Attribution
Evaluate three independent controls on the exact same injected waveforms:
1. **Flat-Noise Control (`gaussian_white`) and `pink_colored` Control** at identical $\rho_{\text{opt}} \in \{10, 14, 18, 22, 26\}$ (measures intrinsic STFT detector efficiency without resonance).
2. **Pre-Whitened Oracle-PSD Control** (apply the exact known inverse AR(2) whitening filter `lfilter(a_coeffs, b_coeffs, x)` before running frozen `MultiResolutionDetector` calibrated on pre-whitened null windows at `1%` FAP): directly measures how much efficiency is recoverable by ideal whitening vs. fundamental information loss (H1 vs. H2)!
3. **Calibration Sample-Size Control:** Compare mature causal `R_1` calibration (`N_thresh = 150`) vs. large independent oracle `shifted_ar2` calibration (`N_thresh = 1000`) at matched `1%` empirical FAR (quantifies H3 directly), alongside `>= 300` independent null windows per background verifying matched `1%` FAR with 90% Wilson CIs.

---

## 4. Quantitative Metrics & Statistical Uncertainty

Every experimental cell will report:
- **Detection efficiency** (raw Stream A, post-veto, and clean-calibrated) with **90% Wilson confidence intervals**.
- **Null background FAR** per calibration mode (`>= 300` null windows) with **90% Wilson CIs**, plus `N_thresh`, `N_eff`, and Ljung-Box $p_{\text{LB}}$.
- **Empirical ROC / Efficiency-vs-$\rho_{\text{opt}}$ curves** (`SNR_50` and `SNR_90`: whitened $\rho_{\text{opt}}$ required for `50%` and `90%` detection efficiency) comparing `gaussian_white`, `pink_colored`, `shifted_ar2 (fixed)`, `shifted_ar2 (randomized)`, and `pre-whitened oracle`.
- **Quantitative Loss Attribution Budget:** Decompose the observed `50%` efficiency gap (`100% -> 50%` in `edd1ec9`) into explicit percentage-point contributions from **H1 (whitened SNR loss)**, **H2 (STFT local-whitening / envelope sensitivity loss)**, **H3 (finite-sample `R_1` threshold effect)**, and **H4 (Stream B gating / coverage)**.

---

## 5. Proposed Acceptance Criteria for the Observability Study

Because this is a **diagnostic observability study** (not an algorithm-tuning phase), acceptance criteria verify experimental rigor, statistical power, and unambiguous causal attribution rather than forcing a target efficiency:

| ID | Criterion | Measurable Target |
| :--- | :--- | :--- |
| **OBS-AC1** | **Complete Forensic Audit of `edd1ec9` Chirps** | All 10 post-recovery chirps (`5` detected, `5` missed) from `edd1ec9` are reproduced verbatim and tabulated with exact `(f_start, f_end, f_res, r, rho_opt, raw/std scores, R1 vs. oracle thresholds, Stream B status)`. |
| **OBS-AC2** | **Matched-$\rho_{\text{opt}}$ vs. Time-Domain SNR Separation (H1 vs. H2)** | Efficiency-vs-placement sweep (`in_resonance`, `crossing`, `off_low`, `off_high`) completed with `N >= 40` trials/cell and 90% Wilson CIs under BOTH fixed time-domain SNR (`18`) and fixed whitened $\rho_{\text{opt}} \in \{10, 14, 18, 22, 26\}$, plus `SNR_50`/`SNR_90` estimates. |
| **OBS-AC3** | **Duration, Sweep-Rate & Envelope-Modulation Controls** | Complete factorial results across all 3 `(duration, sweep_rate)` regimes and the zero-modulation AR(2) ablation, quantifying whether within-window `±25%` modulation contributes to threshold inflation or score loss. |
| **OBS-AC4** | **Independent Controls & Matched Null FARs** | `gaussian_white`, `pink_colored`, `pre-whitened oracle-PSD`, and `N = 1000` oracle-calibration controls all evaluated with matched null FAR (`<= 1.5%` on `>= 300` null windows, with explicit 90% Wilson CIs and `N_eff`). |
| **OBS-AC5** | **Quantitative Four-Way Attribution Budget** | Explicit quantitative decomposition attributing the `edd1ec9` `50%` drop across H1, H2, H3, and H4, with clear empirical verdict on which hypotheses are supported or falsified. |
| **OBS-AC6** | **Strict Code, Baseline & O3 Holdout Isolation** | Zero modifications to any detector/monitor algorithm or parameter; `0` files from `docs/v03_o3_holdout_manifest.json` accessed; reproducible script + CSVs + figures + English scientific report on drive D:. |
