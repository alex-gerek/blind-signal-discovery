# BSD v0.3 Phase 2b Research Design — Calibration Uncertainty, Orthogonal Diagnostics, and Persistent Broadband Sensitivity

**Status:** Approved in Principle with Four Methodological Amendments (Frozen Pre-Implementation Design)
**Branch:** `development/v0.3` (following preserved Phase 2 commit `96fc2c0`)
**Frozen Baselines:** `v0.1` (`detector.py`), `v0.2` (`detector_v2.py`, `vetoes.py`, `coincidence.py`), `v0.3` Phase 1 (`regime_monitor.py`), and `v0.3` Phase 2 (`adaptive_regime.py`) remain unmodified.
**Reserved O3 Holdout (`docs/v03_o3_holdout_manifest.json`):** Strictly untouched and uninspected (`0` local files).

---

## 1. Motivation & Summary of Preserved Phase 2 Findings (`96fc2c0`)

Phase 2 (`AdaptiveRegimeManager`, commit `96fc2c0`) demonstrated finite recovery (`T_recover = 62` windows, `80.0%` duty cycle), elimination of narrowband 80-s sinusoid absorption (`0.0%` loss vs. `70–90%` in Phase 1), and `100%` raw candidate preservation in Stream A during transitions, while honestly failing three acceptance criteria (`P2-AC1`, `P2-AC3`, `P2-AC5`) and exposing detector-level blindness to persistent broadband signals:

1. **Calibration Sample-Size & Variance Failure (`P2-AC1`, `P2-AC3`, `P2-AC5`):** Recalibrating on `N_recal = 60` (`N_thresh = 45` threshold samples) produced an expected i.i.d. false-alarm rate of `2.17%` with a `90%` Beta interval `[0.11%, 6.44%]` (`40.3%` chance of true FAR `> 2%`), causing both elevated FAR (`2.17%` synthetic, `4.97%` on `GW170818 H1`, `8.84%` on `GW150914 H1` symmetric) and threshold-variance efficiency loss (`80.0%` on post-recovery chirps).
2. **Conflated Glitch vs. Continuum Sensitivity (`P2-AC5`):** Removing the time-domain tail ratio to protect narrowband sinusoids reduced sensitivity to impulsive non-Gaussian glitch storms (`GW170729 V1` surviving FAR rose from `1.86%` in Phase 1 to `4.14%` in Phase 2).
3. **Detector-Level Blindness to Persistent Broadband Signals:** Only `10%` of the 80-second broadband stochastic signal windows survived in **Stream A** — before Stream B credibility assessment ran — because `TemplateFreeDetector` performs within-window 2D local median background subtraction across both time and frequency.

Phase 2b addresses these three problems independently while enforcing four strict methodological guardrails: **calibration independence testing**, **multiple-testing control**, **critical multi-regime cross-coherence validation**, and **raw measurement + nomination preservation**.

---

## 2. Problem 1: Calibration Sample Size, Temporal Independence, and Uncertainty

### 2.1 Distinguishing Empirical p-Value Resolution from True FAR Uncertainty

Let $X_1, \dots, X_{N}$ be $N = N_{\text{thresh}}$ calibration null scores, and let $X_{(k)}$ be the $k$-th order statistic selected for target false-alarm probability $\alpha$.

1. **Finite Empirical p-Value Resolution:** $p_{\min} = 1 / (N + 1)$ is strictly the discrete quantization floor of the conformal empirical p-value $p_{\text{emp}}(s) = (1 + \#\{i : X_i \ge s\}) / (N + 1)$. It is **not** a lower bound on the true false-alarm probability.
2. **True Exceedance Probability Under Independence:** Only when $X_1, \dots, X_N$ are independent and identically distributed (exchangeable) from a continuous stationary null $F$ does $U = 1 - F(X_{(k)}) \sim \mathrm{Beta}(N - k + 1,\, k)$ hold exactly:

| $N_{\text{thresh}}$ (i.i.d.) | Order Stat $k$ ($\alpha = 0.01$) | $p_{\min} = \frac{1}{N+1}$ | $\mathbb{E}[\mathrm{FAR}]$ | 90% i.i.d. Beta CI | $P(\mathrm{FAR} > 2\%)$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **45** *(Phase 2)* | 45 | 2.17% | 2.17% | `[0.11%, 6.44%]` | **40.3%** |
| **100** | 100 | 0.99% | 0.99% | `[0.05%, 2.95%]` | 13.3% |
| **150** | 150 | 0.66% | 0.66% | `[0.03%, 1.98%]` | 4.8% |
| **200** | 199 | 0.50% | 1.00% | `[0.18%, 2.35%]` | 8.9% |
| **300** | 298 | 0.33% | 1.00% | `[0.27%, 2.08%]` | 6.0% |
| **400** | 397 | 0.25% | 1.00% | `[0.34%, 1.93%]` | 4.1% |

### 2.2 Explicit Temporal Dependence Testing & Effective Sample Size (Amendment 1)

Real detector background windows are frequently **temporally correlated** (slow seismic/thermal drift, clustered glitches, or overlapping whitening PSD estimates). When windows are temporally dependent, the nominal i.i.d. Beta interval underestimates true FAR uncertainty, and **$N_{\text{thresh}} = 150$ is NOT universally sufficient**.

Therefore, Phase 2b will explicitly test and account for temporal dependence on every calibration ensemble:
1. **Autocorrelation & Independence Diagnostics:** Compute lag-$1 \dots L$ ($L = \min(10, \lfloor N/5 \rfloor)$) sample autocorrelation $\hat{\rho}_k$ on the chronological calibration scores and perform a Ljung-Box portmanteau test ($p_{\text{LB}}$).
2. **Effective Sample Size ($N_{\text{eff}}$):** Estimate Bartlett-weighted effective sample size:
   $$N_{\text{eff}} = \frac{N}{1 + 2 \sum_{k=1}^{L} \left(1 - \frac{k}{L+1}\right) \max(\hat{\rho}_k,\, 0)}$$
   capped to $[1, N]$.
3. **Dependence-Aware Uncertainty Reporting:** Report BOTH the nominal i.i.d. Beta 90% CI (using $N_{\text{thresh}}$) AND the dependence-aware 90% CI (using $N_{\text{eff}}$ and circular block-bootstrap validation where $p_{\text{LB}} < 0.05$).
4. **No Universal Sample-Size Claim:** Maturity promotion (`EMPIRICAL_PROVISIONAL -> EMPIRICAL_VALIDATED`) and `CALIBRATED_EXACT` eligibility must evaluate **$N_{\text{eff}}$**, not merely raw window count $N_{\text{thresh}}$. If strong temporal dependence reduces $N_{\text{eff}}$ below the target threshold, the regime remains explicitly marked `EMPIRICAL_DEPENDENT_UNCERTAIN` with its widened uncertainty interval.

### 2.3 Progressive Causal Calibration Schedule (Zero Future-Data Leakage)

In `causal_trailing` mode, recalibration uses strictly past windows $\le t$:
- **Stage 0 (`N_new < 60`):** `RECALIBRATING_WARMUP` / `UNCERTAIN_UNVALIDATED`.
- **Stage 1 (`60 <= N_new < 200`):** `EMPIRICAL_PROVISIONAL`, reporting `N_thresh`, `N_eff`, `p_min`, and both nominal & dependence-aware 90% FAR CIs.
- **Stage 2 (`N_new >= 200`, `N_thresh >= 150`):** Causal refresh using strictly past windows `[t_start, t]`. Upgraded to `EMPIRICAL_VALIDATED` ONLY if $N_{\text{eff}} \ge 120$ and the dependence-aware upper 90% FAR bound is $\le 2.5\%$; otherwise remains `EMPIRICAL_DEPENDENT_UNCERTAIN`.
- Every run reports `T_provisional` (latency to Stage 1), `T_mature` (latency to Stage 2), `N_thresh`, `N_eff`, $p_{\text{LB}}$, and explicit FAR confidence intervals.

---

## 3. Problem 2: Independent Glitch and Continuum Diagnostics

To resolve Phase 2's tradeoff between narrowband sinusoid preservation and impulsive glitch-storm rejection without collapsing heterogeneous effects into one opaque score, Stream B exposes **two orthogonal, named diagnostic channels**:

1. **Spectral Continuum Monitor (`continuum_status`):**
   - Uses high-resolution (`stft_nperseg=512`, `df = 2 Hz`) per-band median continuum PSD features.
   - Governs regime-shift detection and staged causal recalibration (`STABLE`, `SHIFT_DETECTED`, `RECALIBRATING`).
2. **Impulsive Glitch-Storm Monitor (`glitch_storm_status`):**
   - Operates on short-timescale (`nperseg=32`, `31.25 ms`) **frequency-median frame power** $m(t_j) = \mathrm{median}_f\, |Z(t_j, f)|^2$: taking the median across all frequency bins in each 31-ms frame makes $m(t_j)$ invariant to narrowband tones while remaining highly sensitive to broadband impulsive glitches!
   - Tracks the per-window broadband impulsive crest ratio $\max_j m(t_j) / \mathrm{median}_j m(t_j)$ and rolling elevated-burst rate against $R_0$.
   - **Never triggers regime recalibration** (non-Gaussian glitch storms are tail contamination, not a new stationary Gaussian regime). Emits `QUIET` or `ELEVATED_GLITCH_RATE` independently.
3. **Explicit Provenance:** Every window emits `warning_causes: tuple[str, ...]` (e.g., `("CONTINUUM_DRIFT",)`, `("GLITCH_STORM",)`, both, or `()`) while preserving `100%` of raw measurements and candidates in Stream A.

---

## 4. Problem 3: Persistent Broadband Sensitivity, Multiple-Testing Control, and Cross-Detector Coherence

### 4.1 Separating Underlying Detector Blindness from Background Credibility (Amendment 4)

In frozen `TemplateFreeDetector.analyze` (`src/bsd/detector.py`, lines 169–175), within-window 2D local median subtraction (`bg_psd = median_filter(raw_psd, freq_k)` and `whitened = power / bg_psd`) cancels any stationary broadband excess whose bandwidth exceeds `(freq_k / 2) * df` (`~21–80 Hz`). Thus, the `10%` Stream A hit rate on `persistent_80s_broadband` (`80–300 Hz`) is a **detector-level nomination failure**, distinct from Stream B credibility gating.

Phase 2b explicitly separates:
- **Existing candidate preservation:** preserving `100%` of candidates nominated by frozen `MultiResolutionDetector` across all regime states.
- **New persistent broadband nomination (`FrozenReferenceBroadbandChannel`):** a parallel excess-power channel comparing each window's continuum band PSD against the **frozen** initial reference $R_0$:
  $$E_{\text{ref}}(t) = \max_{b} \frac{\log \hat{S}_b(t) - \mu_{b,R_0}}{\sigma_{b,R_0}} \quad \text{and} \quad E_{\text{sum}}(t) = \sum_{b} \max\!\left(0,\, \frac{\log \hat{S}_b(t) - \mu_{b,R_0}}{\sigma_{b,R_0}}\right)$$
  Every window preserves **both raw scores** (`stft_max_score`, `broadband_excess_score`) **and separate nomination decisions** (`stft_nominated`, `broadband_nominated`, `joint_nominated`). Single-channel broadband nominations are explicitly labeled `SINGLE_CHANNEL_AMBIGUOUS`.

### 4.2 Mandatory Multiple-Testing Control (Amendment 2)

Adding `FrozenReferenceBroadbandChannel` expands the hypothesis space. Comparing it against an unchanged `1%` threshold alongside `MultiResolutionDetector` would inflate the combined family-wise false-alarm rate toward `~2%`.

Phase 2b enforces **strict joint / Sidak multiple-testing control**:
- For target combined family-wise false-alarm probability $\alpha_{\text{target}}$ (e.g., `0.01`), each of the $M = 2$ detection families (`MultiResolutionDetector` STFT bank and `FrozenReferenceBroadbandChannel`) is calibrated at the Sidak-corrected per-family level:
  $$\alpha_{\text{each}} = 1 - (1 - \alpha_{\text{target}})^{1/2} \approx 0.00501$$
  **AND/OR** calibrated jointly via the empirical max-standardized null distribution across both families on the same split-sample calibration set, guaranteeing **matched combined FAR** ($\le \alpha_{\text{target}}$ within finite-sample CI) on stationary null data.
- Transient (`MultiResolutionDetector`) and persistent (`FrozenReferenceBroadbandChannel` + `policy=persistent`) results are reported in **strictly separate tables** as well as in the matched-FAR joint audit.

### 4.3 Critical Evaluation of Cross-Detector Coherence (Amendment 3)

Cross-detector coherence is treated strictly as **additional corroborating evidence**, never automatic proof of an astrophysical signal.

For synchronized two-detector streams $(x_1(t), x_2(t))$, compute multi-window (`M = 10..20` windows = `40–80 s`) magnitude-squared coherence $\hat{C}_{12}(f)$ and lag-resolved cross-correlation over physical light-travel delays $|\Delta t| \le 10\text{ ms}$, calibrated against **non-physical time-slide controls** ($|\Delta t| \in [0.25\text{ s}, 1.75\text{ s}]$).

Phase 2b will systematically test **five mandatory two-detector control regimes**:
1. **Independent stationary noise** (null FAR verification against time-slides).
2. **Uncorrelated dual-detector broadband noise drift** (verifying coherence stays null when both detectors drift independently).
3. **Correlated environmental noise** (shared zero-lag broadband or narrowband instrumental contamination, e.g., magnetic/power-line coupling -> demonstrates that coherence *can* fire on correlated non-astrophysical noise and must flag `ZERO_LAG_OR_COMMON_ENV_WARNING` rather than claiming astrophysical proof!).
4. **Common artifacts** (synchronous shared impulsive glitches or common narrowband spectral lines -> verifies line/glitch artifact discrimination).
5. **Genuine shared synthetic signals with injected physical time delays** ($\Delta t \in \{+4\text{ ms}, -7\text{ ms}\}$) vs. **unphysical delays** ($\Delta t = 50\text{ ms}$), evaluated under both stationary and drifting noise.

---

## 5. Measurable Acceptance Criteria for Phase 2b

All evaluations use newly generated synthetic scenarios and previously examined O1/O2 development epochs ONLY. Transient and persistent results are reported separately; all acceptance criteria require **matched FARs, explicit confidence intervals, and full-timeline detection efficiency**.

| ID | Criterion | Measurable Target |
| :--- | :--- | :--- |
| **P2b-AC1** | **Dependence-Aware Calibration & Matched FAR** | On 500-window synthetic abrupt step-change (onset win 100), `causal_trailing` achieves `T_provisional <= 75` win, `T_mature <= 220` win (zero future data), reports `N_thresh`, `N_eff`, Ljung-Box $p_{\text{LB}}$, and both nominal & dependence-aware 90% FAR CIs; mature-regime surviving FAR `<= 1.5%` (consistent with matched target `1%`). |
| **P2b-AC2** | **Multiple-Testing Control (Matched Combined FAR)** | Jointly calibrated STFT + Frozen-Reference Broadband detector achieves matched combined FAR `<= 1.5%` (with explicit 90% CI) on stationary null background — zero uncompensated FAR inflation. |
| **P2b-AC3** | **Transient Efficiency (Full Timeline & Post-Recovery)** | Stationary transient chirps: full-timeline and calibrated efficiency `>= 95%`; post-recovery transient chirps in mature `R_1` (`1.5x shifted_ar2`): `>= 90%`. |
| **P2b-AC4** | **Orthogonal Glitch vs. Continuum Diagnostics** | On synthetic glitch storm + narrowband 80-s sinusoid: (a) `glitch_storm_monitor` flags `>= 85%` of glitch-storm windows with `GLITCH_STORM` without triggering recalibration; (b) `continuum_monitor` keeps `< 5%` loss on 80-s sinusoid and `<= 15%` loss on 40-s chirp; (c) `100%` of warnings specify `warning_causes` and `100%` of raw scores & nominations are preserved. |
| **P2b-AC5** | **Persistent Broadband, Multitone, Chirp & Nonstationary Matrix** | Under both stationary and drifting noise, `FrozenReferenceBroadbandChannel` raises Stream A broadband 80-s persistent detection from `10%` (frozen STFT) to `>= 90%` at matched FAR, while explicitly separating frozen-STFT vs. frozen-reference nominations and documenting behavior on multitone combs, slow chirps, and broadband noise drift. |
| **P2b-AC6** | **Five-Regime Cross-Detector Coherence Validation** | Calibrated against time-slide controls (`<= 1%` null FAR): achieves `>= 90%` detection on physical-delay shared 80-s broadband signals (`\Delta t = +4, -7 ms`), `<= 2%` false coherence on independent noise and uncorrelated dual-detector drift, and explicitly demonstrates & flags correlated environmental zero-lag noise and common artifacts without misclassifying them as proven astrophysical signals. |
| **P2b-AC7** | **Real O1/O2 Development Epochs & Strict O3 Isolation** | Across the 6 previously examined O1/O2 epochs (`GW150914 H1/L1`, `GW170814 L1`, `GW170104 L1`, `GW170729 V1`, `GW170818 H1`): report full-timeline FAR, calibrated FAR with 90% CI, `N_eff`, duty cycle, excluded hours, and diagnostic causes (`CONTINUUM_DRIFT` vs. `GLITCH_STORM`), preserving `100%` of confirmed `GW150914 H1/L1` & `GW170814 L1` events; `0` O3 holdout files touched; all v0.1/v0.2/Phase 1/Phase 2 baselines preserved. |

---

## 6. Scientific Integrity & Reporting Policy

- Do not tune parameters solely to pass acceptance thresholds; report every regression, negative finding, and correlated-artifact vulnerability transparently.
- Keep `docs/v03_o3_holdout_manifest.json` strain files 100% untouched.
- Commit and push this amended research design to `development/v0.3` before implementation; keep Phase 2b implementation commits local until explicitly approved for push.
