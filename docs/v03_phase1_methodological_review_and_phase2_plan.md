# BSD v0.3 Phase 1 — Methodological Review & Phase 2 Architecture Plan

> **Review Scope:** Methodological audit of `BackgroundRegimeMonitor` (`development/v0.3`, commit `2f6dd2f`) focused on **signal preservation under nonstationary background conditions**.  
> **Approval Status:** Approved with 6 methodological amendments (incorporated below).  
> **Code Freeze Notice:** The 8 reserved O3 events + 4 O3 off-source holdout intervals (`docs/v03_o3_holdout_manifest.json`) remain completely untouched; frozen `v0.1` and `v0.2` baselines remain unmodified.

---

## 1. Diagnosis of ~50% Observation-Time Exclusion & Infinite Recovery Time

In the Phase 1 synthetic benchmark (`results/v03_phase1/report.md`), an abrupt step-change at window 60 of 120 (`abrupt_step_change`: pink noise $\to$ 1.5x shifted AR(2)) yielded:
- `causal_trailing`: **50.0% effective duty cycle** (`60 / 120` windows excluded; `onset_delay = 0` windows).
- `offline_symmetric`: **47.5% effective duty cycle** (`63 / 120` windows excluded; `onset_delay = 0` windows).

**Root Cause:** `BackgroundRegimeMonitor` compares every local block against a **permanently frozen initial `ReferenceBackgroundState`** fitted at the start of the epoch. Once the detector enters a new, internally stationary noise regime at window 60:
1. Every subsequent window `60..119` permanently differs from the old initial reference (`max_feature_z_drift >= 3.5` and Bonferroni KS `p < 0.005`).
2. Furthermore, the one-sided CUSUM statistic accumulates continuously without a reset mechanism.
3. Consequently, **the measured regime recovery time in Phase 1 is $\infty$ (unbounded)**: the monitor never recovers after a genuine transition to a new stationary regime, forcing `100%` of post-transition livetime to be excluded (`50%` of a half-and-half record; `52.5%` in `offline_symmetric` because the centered `12`-window block previews future shifted windows up to `6` steps before window 60).

The same mechanism explains `GW170818 H1` (`53.6%` duty cycle in `causal_trailing`, `50.9%` in `offline_symmetric`): once the H1 background shifted mid-epoch, the remainder of the recording stayed permanently flagged as `REGIME_SHIFT_DETECTED`.

---

## 2. Diagnosis of Signal Absorption Losses (80-s Sinusoid: 70–90%; 40-s Chirp: 30–80%)

Phase 1 quantitatively measured severe signal-preservation losses when long, strong physical signals were injected into stationary noise:

| Signal Scenario | Mode | Efficiency Before Gating | Efficiency After Gating | Signals Lost to Gating |
| :--- | :---: | :---: | :---: | :---: |
| `transient_chirps_stationary` (20 short chirps) | `causal_trailing` | `100.0%` | `95.0%` | `1 / 20` (`5%`) |
| `transient_chirps_stationary` (20 short chirps) | `offline_symmetric` | `100.0%` | `95.0%` | `1 / 20` (`5%`) |
| `persistent_80s_sinusoid` (20 consecutive windows) | `causal_trailing` | `100.0%` | `30.0%` | **`14 / 20` (`70%`)** |
| `persistent_80s_sinusoid` (20 consecutive windows) | `offline_symmetric` | `100.0%` | `10.0%` | **`18 / 20` (`90%`)** |
| `slow_40s_chirp` (10 consecutive windows) | `causal_trailing` | `100.0%` | `70.0%` | **`3 / 10` (`30%`)** |
| `slow_40s_chirp` (10 consecutive windows) | `offline_symmetric` | `100.0%` | `20.0%` | **`8 / 10` (`80%`)** |

**Exact Mechanism:**
1. **Window-scale overlap ($\tau_{\text{sig}} \gtrsim W_{\text{local}}$):** The monitor evaluates `local_window_count = 12` windows (`48 s`). An 80-second sinusoid occupies `20` consecutive windows (`100%` of the local block), raising the time-domain dispersion (`log_mad_sigma`) and the `160–280 Hz` band power (`log_psd_160_280hz`) across the entire local block until KS, Z-drift, and CUSUM trigger.
2. **Causal vs. Symmetric Asymmetry:**
   - In `causal_trailing`, the first `6` windows of the 80-s sinusoid survive before the block median and CUSUM cross threshold (`6 / 20 = 30%` preserved, `70%` lost). For the 40-s chirp (`10` windows < `48 s`), only the final `3` windows cross threshold (`30%` lost).
   - In `offline_symmetric`, the centered block `[i - 6, i + 6]` pulls future signal windows into the local statistics *immediately* at signal onset, collapsing survival to `10%` (`90%` lost) on the 80-s sinusoid and `20%` (`80%` lost) on the 40-s chirp.
3. **Lack of Narrowband Immunity in Background Features:** `extract_window_features` computes `log_mad_sigma` directly on the raw time series and `median(med_psd)` inside each band without first stripping narrowband spectral peaks, allowing a strong narrowband tone to shift broadband regime features.

---

## 3. Epistemic Distinction: Background Shift vs. Instrumental Noise

A detected statistical change in background features is **NOT** evidence that a concurrent anomaly candidate is instrumental noise. Formally:
- **What a regime shift establishes:** The empirical distribution of local background features deviates from the frozen calibration reference, so the calibrated false-alarm probability (e.g., 1% FAP) and empirical `p`-value are **no longer statistically valid**.
- **What a regime shift does NOT establish:** It does not prove the physical origin of the shift. A statistically significant shift can be caused by (a) interferometer hardware drift, (b) a clustered glitch storm, **or** (c) a genuine persistent or slowly evolving astrophysical/physical signal.
- **Required epistemic rule:** Never label a regime-flagged candidate as "instrumental noise" and never silently drop it. Preserve every candidate with an explicit `significance_credibility = UNCERTAIN_REGIME_SHIFT` annotation.

---

## 4. Proposed Dual-Output Architecture

BSD v0.3 Phase 2 will decouple **anomaly detection** from **background statistical credibility** via two parallel, non-destructive output streams for every observation window:

1. **Stream A — Anomaly & Candidate Output (Never Suppressed):**
   - `raw_candidates`: every candidate exceeding the active detector threshold (exact timestamps, frequency bounds, raw & standardized scores, winning scale).
   - `veto_decisions`: reversible Phase 2 artifact flags (`glitch_veto`, `line_veto`) with explicit reason codes and diagnostics.
   - `signal_preservation_guarantee`: candidates occurring during `TRANSITION_UNCREDIBLE` or `RECALIBRATING` intervals remain present in Stream A.
2. **Stream B — Background Stability & Credibility Output (Independent Metadata):**
   - `regime_id`: integer identifier of the active stationary background regime (`R_0, R_1, ...`).
   - `credibility_status`: one of `CALIBRATED_STATIONARY`, `TRANSITION_UNCREDIBLE`, `RECALIBRATING_WARMUP`, or `PERSISTENT_NARROWBAND_ANOMALY_SUSPECTED`.
   - `significance_status`:
     - `CALIBRATED_EXACT`: used **ONLY** when the calibration method genuinely justifies exact finite-sample exchangeability (large independent reference ensemble `N_cal >= 200` on an unshifted stationary regime `R_0`).
     - `EMPIRICAL_FINITE_SAMPLE`: used after post-change-point local recalibration with smaller `N_recal`, explicitly bounded by `p >= 1 / (N_recal + 1)`.
     - `UNCERTAIN_UNVALIDATED`: used during `TRANSITION_UNCREDIBLE` and `RECALIBRATING_WARMUP` intervals; candidates are preserved in Stream A without deleting them.
   - `mode_metadata`: `causal_trailing` vs. `offline_symmetric`, `future_data_used`, and exact `latency_seconds`.
   - `accounting`: `excluded_livetime_hours`, `effective_duty_cycle`, `regime_recovery_time_s`, and `candidates_in_uncalibrated_intervals`.

---

## 5. Separate Policies for Transient vs. Persistent Signal Searches

The nominal `10 s` boundary (`transient_duration_horizon_s`) is treated strictly as an **experimental configuration parameter**, not a fundamental physical law; Phase 2 benchmarks must explicitly probe signal behavior *around* this boundary (e.g., `4 s`, `12 s`, `40 s`, and `80 s` durations).

| Policy Dimension | **Transient Search Policy** (configurable $\tau_{\text{sig}} \le \tau_{\text{horizon}}$) | **Persistent / Slow-Signal Policy** ($\tau_{\text{sig}} > \tau_{\text{horizon}}$) |
| :--- | :--- | :--- |
| **Background Features** | Spectral-continuum features (frequency median filter `kernel=11`) — treated as a **testable hypothesis**, not guaranteed absorption immunity | Strictly frozen multi-hour reference + **cross-detector coherence check** (H1 vs. L1/V1) to distinguish single-site noise drift from coherent persistent signals |
| **Post-Change-Point Recalibration** | **Allowed** once a new stationary regime is confirmed (`dwell >= D_confirm` and internal stability verified over `N_recal` windows); note that a stable new regime is **not** proof that the original shift was instrumental noise | **Prohibited** on single-channel rolling windows (prevents absorbing continuous tones or long inspirals) |
| **Treatment of Persistent Narrowband Excess** | Flagged separately (`PERSISTENT_NARROWBAND_ANOMALY_SUSPECTED` / reversible `line_veto`) rather than absorbed into broadband regime shift | Treated as a primary signal candidate; evaluated against frozen reference + multi-detector consistency |

---

## 6. Proposed Phase 2 Architecture (`AdaptiveRegimeManager`)

To achieve **finite regime recovery time** after abrupt transitions *while rigorously testing* persistent-signal absorption, Phase 2 introduces three targeted mechanisms (while keeping `v0.1` and `v0.2` frozen):

1. **Spectral-Continuum Background Features (Hypothesis Under Test):**
   - Before computing `log_band_psd` and broadband continuum energy, apply a frequency-domain median filter (`kernel = 11` bins) to each window's STFT power spectrum.
   - **Explicit Validation Requirement:** Treat narrowband immunity as a **hypothesis to be tested**, not a guarantee. Benchmark against four distinct persistent/slow morphologies: (a) 80-s single-tone sinusoid, (b) 80-s multi-tone comb, (c) 40-s slowly sweeping chirp, and (d) 80-s **broadband persistent** signal (where broadband continuum features are expected to exhibit measurable absorption, which must be reported honestly).
2. **Three-State Regime Finite-State Machine with Hysteresis (Independent Validation):**
   - `CALIBRATED_STATIONARY` $\xrightarrow{\text{KS / Z / CUSUM shift}}$ `TRANSITION_UNCREDIBLE` (every raw candidate preserved in Stream A with `significance_status = UNCERTAIN_UNVALIDATED`).
   - `TRANSITION_UNCREDIBLE` $\xrightarrow{\text{internal stability over } N_{\text{stable}} \text{ windows}}$ `RECALIBRATING_WARMUP` (all raw candidates still preserved with `significance_status = UNCERTAIN_UNVALIDATED`; collect `N_recal` internally stationary windows to fit refreshed reference + local threshold).
   - `RECALIBRATING_WARMUP` $\xrightarrow{N_{\text{recal}} \text{ complete}}$ `CALIBRATED_STATIONARY (regime_id += 1, significance_status = EMPIRICAL_FINITE_SAMPLE)`.
   - Change-point detection, regime recovery, and post-recovery calibration credibility are validated **independently**; achieving a stable new regime does not prove the transition was instrumental noise.
3. **Full-Timeline Reporting:**
   - Report both candidate false-alarm rates (per-window FAP and candidates/hr) and detection efficiency over the **entire observation timeline**, explicitly accounting for excluded time, effective duty cycle, and candidates flagged as uncertain during transitions/warmup.

---

## 7. Quantitative Summary of Phase 1 Metrics & Proposed Phase 2 Acceptance Criteria

### 7.1 Phase 1 Measured Baseline (to Beat in Phase 2)

- **Observation coverage (`abrupt_step_change`):** `50.0%` (`causal`), `47.5%` (`symmetric`); recovery time = $\infty$.
- **Observation coverage (`GW170818 H1`):** `53.6%` (`causal`), `50.9%` (`symmetric`); FAR `15.40% -> 3.10%` (`causal`), `1.40%` (`symmetric`).
- **Signal loss (`persistent_80s_sinusoid`):** `70%` (`causal`), `90%` (`symmetric`).
- **Signal loss (`slow_40s_chirp`):** `30%` (`causal`), `80%` (`symmetric`).
- **Signal loss (`transient_chirps_stationary`):** `5%` (`1/20`).

### 7.2 Predefined Phase 2 Acceptance Criteria (Before Any O3 Holdout Access)

- **P2-AC1 (Finite Regime Recovery & Duty Cycle):** On a 300-window synthetic abrupt step-change record (shift at window 100), `causal_trailing` must achieve finite recovery time $T_{\text{recover}} \le 80$ windows (`320 s`), raising effective duty cycle from `~33%` (Phase 1 static) to `>= 70%` while maintaining empirical FAR `<= 1.5%` on calibrated surviving livetime (no false credit for excluded intervals).
- **P2-AC2 (Persistent Signal Preservation):** Under the narrowband-immune continuum monitor, `persistent_80s_sinusoid` false regime-gating loss must drop from `70–90%` (Phase 1) to **`< 10%`** in `causal_trailing` mode, and `slow_40s_chirp` loss must drop from `30%` to **`<= 15%`**.
- **P2-AC3 (Transient Signal Preservation):** On `transient_chirps_stationary`, end-to-end efficiency after gating must remain **`>= 95%`** (and `>= 90%` on transients embedded in a post-recovery shifted regime).
- **P2-AC4 (Non-Destructive Dual Output):** `100%` of raw candidates (including those during `TRANSITION_UNCREDIBLE`) must be preserved in Stream A with explicit `significance_status = UNCERTAIN_UNVALIDATED` (verified by automated pytest).
- **P2-AC5 (O1/O2 Development Real-Data Validation):** On `GW170818 H1` and `GW170729 V1`, calibrated surviving-livetime FAR must stay **`<= 2.0%`** with improved duty cycle over Phase 1, while preserving confirmed `GW150914 H1/L1` and `GW170814 L1` detections.
- **P2-AC6 (Holdout Isolation):** The 8 reserved O3 events and 4 O3 off-source intervals (`docs/v03_o3_holdout_manifest.json`) remain strictly uninspected and untouched until full v0.3 freeze.
