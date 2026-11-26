# Blind Signal Discovery (BSD) — Research Proposal: Distinguishing Signals from Nonstationary Background Shifts

> **Project Independence:** Standalone BSD research program, completely separate from PHARMAKON DSP, `libpharmakon`, and Periodic Table of Sound.  
> **Proposal Status:** Approved in Principle (Design Proposal + Frozen O3 Holdout Manifest ONLY — no v0.3 code changes have been made).

---

## 1. Formal Status of BSD v0.2 and Holdout Consumption Record

- **BSD v0.2 Status:** **Evaluated Research Prototype — NOT a Validated Successor to `v0.1-baseline`.**
  - On `development/v0.2` (commit `af9ad5a`), v0.2 passed `AC-1` (synthetic 1% FAP), `AC-3` (impulsive glitch FAR `94.3% -> 7.7%` with `88.0%` broadband pulse recall), `AC-4` (LIGO `GW170608 H1` line and `GW170817 L1` glitch vetoed while preserving confirmed chirps), and Phase 3 sub-10 ms timing validation (`0.37–1.08 ms` MAE on injections; `GW150914` `+7.32 ms`, `r = -0.43`).
  - However, v0.2 **FAILED `AC-2`** (synthetic cross-scale penalty: `linear_chirp -3.47 pp`, withheld `sinusoidal_fm -3.65 pp`) and **FAILED `AC-5` on holdout epochs** (pooled holdout background FAP `2.57%` vs. `1.0%` target, driven by severe hourly background drift in `GW170818 H1` at `15.40%` and `GW170729 V1` at `7.57%`).
  - Accordingly, `development/v0.2` remains unmerged; `main` preserves `v0.1-baseline` and Experiments 001–003 intact.
- **Holdout Consumption Record:** The four reserved GWTC-1 holdout events (`GW170729`, `GW170809`, `GW170818`, `GW170823`) have now been **CONSUMED** in the one-shot v0.2 evaluation (`results/v02_holdout/report.md`). All 11 GWTC-1 confident events have been examined and **must never be reused as untouched final evaluation data** for v0.3 or later versions.

---

## 2. Central Research Question

> **How can a template-free anomaly detector distinguish genuine transient or persistent signals from changes in the statistical distribution of background noise?**

In BSD v0.2, calibration assumed **temporal exchangeability** across a 4096-second observing epoch: early segments (`0–2240 s`) calibrated the threshold applied to later segments (`2240–4096 s`). When the interferometer's noise distribution shifted mid-epoch (e.g., `GW170818 H1`), the static threshold lost false-alarm control (`15.40%` FAR), with the `500 ms` Long scale exhibiting the highest sensitivity to slow spectral drift.

---

## 3. Comparative Analysis of Nonstationary Background Approaches

| Approach | Core Mechanism | Strengths | Failure Modes & Risks |
| :--- | :--- | :--- | :--- |
| **1. Rolling / Symmetric Local Calibration** | Estimate PSD and empirical null threshold from a local sliding window `[t - W, t + W]` excluding a central guard `[t - G, t + G]`, updated continuously in time. | Adapts smoothly to slow (minute-to-hour) noise drift; preserves exact local exchangeability when drift timescale $\tau_{\text{drift}} \gg W$. | Finite-sample variance: a short local window $W$ contains too few 4-s segments for a direct 1% empirical quantile (e.g., $W = 400\text{ s} \Rightarrow N = 100$ segments, tail resolved by only 1 sample). |
| **2. Robust Local Noise Models** | Estimate local per-frequency PSD and temporal dispersion via trimmed/median/MAD or Huber-weighted estimators on short local blocks, coupled with semi-parametric tail modeling (e.g., Generalized Pareto / Peaks-Over-Threshold on local exceedances) or scale-normalized residuals. | Handles moderate nonstationarity and sporadic glitches without inflating PSD; higher sample efficiency than purely nonparametric rolling quantiles. | Parametric tail misspecification can understate heavy-tailed glitch rates; aggressive local whitening can distort broad spectral features. |
| **3. Explicit Change-Point Detection & Regime Gating** | Run a dedicated background-stationarity monitor (e.g., two-sample Kolmogorov-Smirnov / Anderson-Darling or Page-Hinkley / PELT on rolling PSD & score distributions) to detect when the noise regime shifts. When a shift is detected, either trigger recalibration on the new stationary regime or explicitly mark the interval `BACKGROUND_REGIME_SHIFT_UNCREDIBLE`. | **Scientifically honest:** directly separates *stationary calibrated operation* from *nonstationary regime transitions*; prevents silent FAR inflation without altering signal whitening. | Cannot produce calibrated detections during active transition intervals; requires minimum dwell time after a change point before a new threshold is statistically valid. |

### 3.1 Explicit Distinction: Offline Symmetric vs. Causal Trailing-Window Calibration

Every v0.3 calibration study **must explicitly declare and separate two distinct operating modes**, never conflating their results:

| Mode | Calibration Window Relative to Target `t` | Future-Data Dependence | Minimum Latency | Intended Use Case |
| :--- | :---: | :---: | :---: | :--- |
| **`offline_symmetric`** | `[t - W, t + W] \ [t - G, t + G]` | **Yes** (uses future background up to `t + W`) | `>= W + G` (e.g., `>= 10–20 min`) | Archival / retrospective offline reanalysis |
| **`causal_trailing`** | `[t - W - G, t - G]` | **None** (strictly past data only) | `<= 1` observation window (`4 s`) | Online / low-latency streaming monitoring |

Both modes must be benchmarked separately under step-change and drifting noise regimes, because `offline_symmetric` can bracket a slow linear drift symmetrically while `causal_trailing` lags behind abrupt upward noise-floor transitions.

### 3.2 Mandatory Regime-Gating Reporting Requirements

Whenever change-point detection or regime gating is active, every benchmark report **must** publish all five metrics:
1. **`excluded_livetime_hours`** — total observation time flagged as uncalibrated/nonstationary and excluded from calibrated detection.
2. **`effective_duty_cycle`** — fraction of total observation time remaining in calibrated operation (`1 - excluded_fraction`).
3. **`signals_lost_to_gating`** — count and fraction of true injected (or confirmed real) signals falling inside gated-out intervals.
4. **`efficiency_before_gating` vs. `efficiency_after_gating`** — end-to-end detection efficiency measured over *all* target signals (not merely those inside surviving livetime).
5. **No-False-Credit Rule:** Excluded (gated-out) background intervals **must NEVER be counted as successful noise rejections (true negatives)** when computing empirical false-alarm probability or false-alarm rate per hour.

### 3.3 Recommended Hybrid Direction (Subject to Review)

Combine **(3) Explicit Change-Point Regime Gating** with **(1) Local Background Estimation** (evaluated separately in `offline_symmetric` and `causal_trailing` modes) and **Multi-Detector Asymmetry Checks**:
- A genuine astrophysical signal is **localized in time** (or coherent across geographically separated detectors), whereas an instrumental background-distribution shift is **persistent across surrounding guard windows** and typically **uncorrelated across independent sites** (H1 vs. L1 vs. V1).

---

## 4. The Signal Absorption Risk

Any adaptive background estimator introduces a fundamental trade-off: **if the adaptation timescale $\tau_{\text{bg}}$ approaches or overlaps the signal duration $\tau_{\text{sig}}$, the background model will absorb the signal and erase it.**

- **Failure modes of over-adaptation:**
  1. *Persistent / continuous signals* (e.g., stationary sinusoids, continuous-wave tones): a rolling PSD or narrowband line tracker will classify the signal as "new stationary background" and whiten it away.
  2. *Slowly evolving chirps / long inspirals* (e.g., BNS `GW170817`, duration $>30\text{ s}$): a short temporal variance tracker or narrow guard zone $G < 30\text{ s}$ will partially absorb the inspiral energy into $T_{\text{bg}}(t)$ or local calibration quantiles.
- **Required safeguards before implementation:**
  1. **Explicit timescale separation:** Declare the target signal-duration horizon $\tau_{\max}$ upfront and enforce a guard zone $G \ge 2 \tau_{\max}$ around the test window so no target energy enters local PSD or threshold estimation.
  2. **Mode separation (Transient vs. Persistent):** Never claim a single adaptive mode is simultaneously optimal for 20 ms bursts and 100 s continuous tones. Evaluate transient-adapted and persistent-sensitive modes separately with explicit absorption tests.
  3. **Direct absorption stress test:** Inject synthetic persistent sinusoids and 20–40 s slowly sweeping chirps directly into nonstationary drifting noise and measure **signal suppression loss** (recall drop relative to oracle stationary background).

---

## 5. New Independent Final Evaluation Strategy (Defined Prior to Algorithm Changes)

Because all 11 GWTC-1 events are now consumed, any v0.3 prototype must be validated against a **brand-new, strictly untouched evaluation hierarchy** frozen *before* writing v0.3 code:

1. **Synthetic Nonstationary Benchmark Suite (New Disjoint Seeds):**
   - Explicitly model three nonstationary failure regimes: (a) step-change variance/PSD regime shifts, (b) slowly drifting AR(2) spectral resonances, and (c) clustered bursty glitch storms.
   - Include a dedicated **Signal Absorption Control** (long 10–30 s chirps and stationary tones injected across regime boundaries).
   - Report disaggregated FAR metrics (per-window FAP, triggered windows/hr, candidate rate/hr, coincident FAR/hr) and explicit `regime_gating_duty_cycle` (fraction of livetime flagged as uncalibrated/nonstationary — preventing a trivial detector from "passing" FAR simply by flagging 50% of data as uncalibrated).
2. **Real-Data Development Set:**
   - Reuse the already-examined GWTC-1 epochs (`GW170818 H1`, `GW170729 V1`, `GW170104 L1`, `GW150914`, `GW170814`, `GW170817`) strictly as **development diagnostics** to test whether nonstationary drift is properly detected or stabilized.
3. **Frozen Untouched Real-Data Final Holdout Set (GWTC-2.1 / GWTC-3 — O3a & O3b):**
   - Pre-registered and frozen **prior to any v0.3 implementation** in [docs/v03_o3_holdout_manifest.json](v03_o3_holdout_manifest.json):
     - **8 untouched O3 events** (`GW190425` BNS, `GW190521_074359` BBH, `GW190708_232457` BBH, `GW191230_180458` BBH, `GW200115_042309` NSBH, `GW200129_065458` BBH, `GW200224_222234` BBH, `GW200225_060421` BBH; GPS `1240215503`–`1266645879`, zero overlap with O1/O2 GWTC-1 data).
     - **4 deterministic blind off-source null intervals** (`1024 s` each, fixed `+1024 s` offset from 4 anchor O3 events with `>= 512 s` catalog clearance — selected algorithmically without examining strain data, anomaly scores, or noise stationarity).
     - Explicit exclusion guard boundaries: `±64 s` for BBH and `±128 s` for BNS/NSBH around target 32-s windows.
   - **Strict Lock Policy:** Do **not** download, inspect, or run on any O3 holdout strain file until v0.3 code, hyperparameters, calibration mode (`offline_symmetric` vs. `causal_trailing`), and regime-gating rules are completely frozen.
   - Evaluate frozen `v0.1-baseline`, frozen `v0.2`, and frozen `v0.3` in a single locked one-shot pass on the O3 holdout set, reporting transient and persistent evaluations separately alongside pre/post-gating efficiency and effective duty cycle.
