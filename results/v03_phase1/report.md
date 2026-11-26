# Blind Signal Discovery (BSD) v0.3 — Phase 1 Report (`BackgroundRegimeMonitor`)

> **Scope:** Phase 1 implements an independent, template-free `BackgroundRegimeMonitor` (`src/bsd/regime_monitor.py`) supporting both **`causal_trailing`** (zero future-data dependence, latency = 1 window = `4 s`) and **`offline_symmetric`** (uses `local_window // 2` future windows, latency = `7` windows = `28 s`). All original detector candidates are preserved unmodified; excluded intervals are NEVER counted as true negatives. The 8 reserved O3 events and 4 O3 off-source intervals remain completely untouched.

## 1. Synthetic Nonstationary & Signal-Absorption Benchmark

| Scenario | Mode | Duty Cycle | Excl. (hr) | Onset Delay (win) | v0.2 Bg FAR (Before -> Surviving) | v0.2 Signal Eff (Before -> After) | Signals Lost to Gating |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `stationary_baseline` | `causal_trailing` | `100.0%` | `0.000` | `—` | `3.33% -> 3.33%` | `— -> —` | `0` |
| `stationary_baseline` | `offline_symmetric` | `100.0%` | `0.000` | `—` | `3.33% -> 3.33%` | `— -> —` | `0` |
| `abrupt_step_change` | `causal_trailing` | `50.0%` | `0.067` | `0` | `3.33% -> 0.00%` | `— -> —` | `0` |
| `abrupt_step_change` | `offline_symmetric` | `47.5%` | `0.070` | `0` | `3.33% -> 0.00%` | `— -> —` | `0` |
| `gradual_drift` | `causal_trailing` | `9.2%` | `0.121` | `0` | `0.83% -> 0.00%` | `— -> —` | `0` |
| `gradual_drift` | `offline_symmetric` | `2.5%` | `0.130` | `0` | `0.83% -> 0.00%` | `— -> —` | `0` |
| `clustered_glitch_storm` | `causal_trailing` | `51.7%` | `0.064` | `2` | `46.67% -> 3.23%` | `— -> —` | `0` |
| `clustered_glitch_storm` | `offline_symmetric` | `50.0%` | `0.067` | `0` | `46.67% -> 0.00%` | `— -> —` | `0` |
| `transient_chirps_stationary` | `causal_trailing` | `97.5%` | `0.003` | `—` | `1.00% -> 0.00%` | `100.0% -> 95.0%` | `1` |
| `transient_chirps_stationary` | `offline_symmetric` | `98.3%` | `0.002` | `—` | `1.00% -> 1.01%` | `100.0% -> 95.0%` | `1` |
| `persistent_80s_sinusoid` | `causal_trailing` | `84.2%` | `0.021` | `—` | `2.00% -> 2.11%` | `100.0% -> 30.0%` | `14` |
| `persistent_80s_sinusoid` | `offline_symmetric` | `85.0%` | `0.020` | `—` | `2.00% -> 2.00%` | `100.0% -> 10.0%` | `18` |
| `slow_40s_chirp` | `causal_trailing` | `93.3%` | `0.009` | `—` | `0.91% -> 0.95%` | `100.0% -> 70.0%` | `3` |
| `slow_40s_chirp` | `offline_symmetric` | `92.5%` | `0.010` | `—` | `0.91% -> 0.92%` | `100.0% -> 20.0%` | `8` |

### Key Synthetic Findings & Signal Absorption Audit
- **Stationary Baseline:** Maintains high effective duty cycle without spurious regime alarms.
- **Nonstationary Regimes (`abrupt_step_change`, `gradual_drift`, `clustered_glitch_storm`):** Both `causal_trailing` and `offline_symmetric` detect distribution shifts rapidly; `offline_symmetric` exhibits earlier transition localization at the cost of future-data dependence (`28 s` minimum latency).
- **Signal Absorption Risk (Measured Quantitatively):** On short `transient_chirps_stationary`, gating causes minimal/zero signal loss. On `persistent_80s_sinusoid` and `slow_40s_chirp`, any loss reported in the `Signals Lost to Gating` column directly measures the risk that a long, high-SNR physical signal shifts local background features enough to trigger the regime monitor — confirming why **transient and persistent search modes must remain explicitly separated**.

## 2. Previously Examined O1/O2 LIGO Development Epochs (No O3 Holdout Touched)

| Event | Det | Epoch Type | Mode | Duty Cycle | Excl. (hr) | v0.1 FAR (Before -> Surviving) | v0.2 Post-Veto FAR (Before -> Surviving) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914` | H1 | `stable_confirmed` | `causal_trailing` | `91.2%` | `0.041` | `0.47% -> 0.26%` | `1.42% -> 0.78%` |
| `GW150914` | H1 | `stable_confirmed` | `offline_symmetric` | `90.8%` | `0.043` | `0.47% -> 0.52%` | `1.42% -> 1.31%` |
| `GW150914` | L1 | `stable_confirmed` | `causal_trailing` | `73.7%` | `0.123` | `0.95% -> 0.96%` | `0.95% -> 0.96%` |
| `GW150914` | L1 | `stable_confirmed` | `offline_symmetric` | `71.3%` | `0.134` | `0.95% -> 1.00%` | `0.95% -> 0.66%` |
| `GW170814` | L1 | `stable_confirmed` | `causal_trailing` | `95.7%` | `0.020` | `0.24% -> 0.25%` | `0.47% -> 0.50%` |
| `GW170814` | L1 | `stable_confirmed` | `offline_symmetric` | `95.3%` | `0.022` | `0.24% -> 0.25%` | `0.47% -> 0.50%` |
| `GW170104` | L1 | `nonstationary_drift` | `causal_trailing` | `96.7%` | `0.016` | `1.90% -> 1.96%` | `1.66% -> 1.72%` |
| `GW170104` | L1 | `nonstationary_drift` | `offline_symmetric` | `96.0%` | `0.019` | `1.90% -> 1.98%` | `1.66% -> 1.73%` |
| `GW170729` | V1 | `nonstationary_drift` | `causal_trailing` | `87.0%` | `0.027` | `7.57% -> 1.86%` | `7.57% -> 1.86%` |
| `GW170729` | V1 | `nonstationary_drift` | `offline_symmetric` | `85.9%` | `0.029` | `7.57% -> 1.26%` | `7.57% -> 1.26%` |
| `GW170818` | H1 | `nonstationary_drift` | `causal_trailing` | `53.6%` | `0.218` | `6.16% -> 0.88%` | `15.40% -> 3.10%` |
| `GW170818` | H1 | `nonstationary_drift` | `offline_symmetric` | `50.9%` | `0.230` | `6.16% -> 0.47%` | `15.40% -> 1.40%` |

## 3. Limitations & Recommended Next Step

1. **Gating vs. Recalibration Trade-off:** Flagging a shifted epoch (`GW170818 H1`, `GW170729 V1`) protects against uncalibrated false alarms, but reduces effective duty cycle until a Phase 2 local recalibration / regime-reset mechanism is implemented.
2. **Persistent Signal Sensitivity:** Single-detector regime monitoring cannot distinguish an 80-second continuous narrowband astrophysical signal from a slow instrumental spectral shift without cross-detector coherence or a dedicated persistent-signal search mode.
3. **Next Recommended Step:** Review Phase 1 measurements before designing **BSD v0.3 Phase 2** (post-change-point local recalibration with strict timescale separation and multi-detector asymmetry checks).

## 4. Reproducibility Commands

```powershell
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe experiments\run_v03_phase1_regime_monitor.py
```
