# BSD v0.3 Phase 2b Scientific Report

**Status:** Evaluated Locally on Synthetic Scenarios & Previously Examined O1/O2 Development Epochs ONLY
**Reserved O3 Holdout Status:** Strictly Untouched (`0` local files from `docs/v03_o3_holdout_manifest.json`)
**Preserved Baselines:** `v0.1`, `v0.2`, `v0.3 Phase 1`, and `v0.3 Phase 2` completely unmodified.

## 1. Acceptance Criteria Summary

| Criterion | Status | Measured Result |
| :--- | :---: | :--- |
| **P2b-AC1** | **FAIL** | T_prov=62, T_mature=201 (<=220), mature_FAR=N/A (<=1.5%), N_thresh=150, N_eff=115.0, p_LB=0.433 |
| **P2b-AC2** | **PASS** | matched_joint_FAR=1.00% (90% CI [0.33%, 2.98%]) |
| **P2b-AC3** | **FAIL** | stat_chirps=100.0% (>=95%), post_mature_chirps=0.0% (>=90%) |
| **P2b-AC4** | **PASS** | glitch_flag=100.0% (n_regimes=1), 80s_sin=100.0%, 40s_chirp=90.0% |
| **P2b-AC5** | **PASS** | 80s_broadband_frozen_STFT=0.0% -> frozen_ref_BB=100.0% (>=90%) |
| **P2b-AC6** | **PASS** | +4ms=100%, -7ms=100%, null=0%, zero_lag_env_warn=100%, line_art_warn=100% |
| **P2b-AC7** | **FAIL** | GW170729_V1_FAR=2.58% (<=2.0%), GW170818_H1_FAR=4.30% (<4.97% P2), confirmed_preserved=False, O3_untouched=True |

## 2. Suite A — Transient, Staged Causal Calibration & Multiple-Testing Control

| Scenario | N | Clean Duty | STFT FAR | Joint FAR | N_thresh | N_eff | Ljung-Box p | T_prov | T_mature | Mature FAR | Full Eff | Clean Eff | Glitch Flag | Regimes |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `stationary_null_200win` | 200 | 100.0% | 0.50% | 1.00% | 224 | 217.3 | 0.874 | N/A | N/A | N/A | N/A | N/A | 0.0% | 1 |
| `abrupt_step_500win` | 500 | 20.4% | 0.62% | N/A | 150 | 115.0 | 0.433 | 62 | 201 | N/A | 50.00% | 0.00% | 79.2% | 2 |
| `transient_chirps_stationary` | 120 | 100.0% | 0.00% | 0.00% | 224 | 217.3 | 0.874 | N/A | N/A | N/A | 100.00% | 100.00% | 0.0% | 1 |
| `glitch_storm_120win` | 120 | 50.0% | 47.50% | 50.00% | 224 | 217.3 | 0.874 | N/A | N/A | 0.00% | N/A | N/A | 100.0% | 1 |

## 3. Suite B — Separate Persistent Signal Matrix (Stationary vs. Drifting Noise)

Preserves raw scores and separates frozen `MultiResolutionDetector` nominations from Sidak-controlled `FrozenReferenceBroadbandChannel` nominations (`policy=persistent`, zero single-channel recalibration).

| Morphology | Background | Frozen STFT Full Eff | Frozen-Ref BB Full Eff | Joint Full Eff | Clean Cal STFT Eff | BG Joint FAR | Continuum Drift Flag | Glitch Flag | Regimes |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `80s_sinusoid` | `stationary` | 100.0% | 0.0% | 100.0% | 100.0% | 2.00% | 0.0% | 0.0% | 1 |
| `80s_multitone_comb` | `stationary` | 100.0% | 5.0% | 100.0% | 100.0% | 2.00% | 0.0% | 0.0% | 1 |
| `40s_slow_chirp` | `stationary` | 90.0% | 10.0% | 90.0% | 90.0% | 0.91% | 0.0% | 0.0% | 1 |
| `80s_broadband_stochastic` | `stationary` | 0.0% | 100.0% | 100.0% | 0.0% | 0.00% | 65.0% | 0.0% | 1 |
| `80s_sinusoid` | `drifting_1.0_to_1.4x` | 100.0% | 100.0% | 100.0% | 0.0% | 69.00% | 100.0% | 0.0% | 1 |
| `80s_multitone_comb` | `drifting_1.0_to_1.4x` | 95.0% | 100.0% | 100.0% | 0.0% | 75.00% | 100.0% | 0.0% | 1 |
| `40s_slow_chirp` | `drifting_1.0_to_1.4x` | 70.0% | 100.0% | 100.0% | 0.0% | 72.73% | 100.0% | 0.0% | 1 |
| `80s_broadband_stochastic` | `drifting_1.0_to_1.4x` | 0.0% | 100.0% | 100.0% | 0.0% | 73.00% | 100.0% | 0.0% | 1 |

## 4. Suite C — Two-Detector Cross-Coherence Five-Regime Validation

Coherence is treated strictly as additional evidence calibrated against 120 circular time-slides (`alpha = 0.01`), never automatic astrophysical proof.

| Control / Signal Regime | Trials | Null / Incoherent | Common Artifact Warning | Zero-Lag Env Warning | Physical-Lag Evidence | Median Lag (ms) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `independent_noise` | 20 | 100% | 0% | 0% | 0% | +0.98 |
| `uncorrelated_dual_drift` | 20 | 100% | 0% | 0% | 0% | +3.91 |
| `correlated_zero_lag_env_noise` | 20 | 0% | 0% | 100% | 0% | +0.00 |
| `common_narrowband_line_artifact` | 20 | 0% | 100% | 0% | 0% | +0.00 |
| `common_synchronous_glitch_artifact` | 20 | 100% | 0% | 0% | 0% | -1.95 |
| `shared_broadband_pos4ms` | 20 | 0% | 0% | 0% | 100% | -3.91 |
| `shared_broadband_neg7ms` | 20 | 0% | 0% | 0% | 100% | +6.84 |
| `shared_broadband_unphysical_50ms` | 20 | 100% | 0% | 0% | 0% | +0.00 |

## 5. Suite D — Previously Examined O1/O2 LIGO/Virgo Development Epochs (Causal Trailing)

| Event / Det | Epoch Type | Clean Duty | Excl. (hr) | Cont. Warn | Glitch Warn | Regimes | Full FAR | Clean Cal FAR (90% CI) | N_eff / N_th | p_LB | Confirmed Event |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914 H1` | `stable_confirmed` | 92.9% | 0.033 | 30 | 0 | 1 | 0.00% | 0.00% [0.00%, 0.69%] | 305/350 | 0.948 | YES |
| `GW150914 L1` | `stable_confirmed` | 92.7% | 0.034 | 29 | 2 | 1 | 0.00% | 0.00% [0.00%, 0.69%] | 346/350 | 1.000 | NO |
| `GW170814 L1` | `stable_confirmed` | 97.6% | 0.011 | 10 | 0 | 1 | 0.00% | 0.00% [0.00%, 0.65%] | 339/350 | 1.000 | NO |
| `GW170104 L1` | `nonstationary_drift` | 99.5% | 0.002 | 0 | 2 | 1 | 1.66% | 1.43% [0.74%, 2.74%] | 268/350 | 0.182 | N/A |
| `GW170729 V1` | `nonstationary_drift` | 83.8% | 0.033 | 16 | 27 | 1 | 7.57% | 2.58% [1.16%, 5.62%] | 146/194 | 0.548 | N/A |
| `GW170818 H1` | `nonstationary_drift` | 71.6% | 0.133 | 120 | 0 | 3 | 8.53% | 4.30% [2.76%, 6.67%] | 262/350 | 0.317 | N/A |

## 6. Figures

- `figures/fig01_staged_causal_recovery_500win.png`
- `figures/fig02_persistent_matrix_streamA_vs_streamB.png`
- `figures/fig03_cross_coherence_five_regimes.png`
- `figures/fig04_ligo_o1o2_orthogonal_diagnostics.png`

## 7. Honest Root-Cause Analysis of Failures & Regressions (`P2b-AC1`, `P2b-AC3`, `P2b-AC7`)

1. **`P2b-AC1` & `P2b-AC3` (Unrefreshed `glitch_ref` after regime shift):** Staged causal recalibration succeeded on `abrupt_step_500win` (`T_prov = 62`, `T_mature = 201 <= 220`, `N_thresh = 150`, `N_eff = 115.0`, `p_LB = 0.433`, full-timeline `STFT FAR = 0.62%`), shrinking the 90% dependence-aware FAR bound from `7.3%` down to `2.57%`. However, `GlitchReferenceState` remained frozen to `R_0` (`pink_colored`) after transitioning to `R_1` (`shifted_ar2`), causing `glitch_storm_monitor` to flag `79.2%` of windows (all 400 `shifted_ar2` windows) and leaving zero `clean_calibrated` windows in mature `R_1` (`mature_FAR = N/A`, `clean_eff = 0.0%`).
2. **`P2b-AC3` & `P2b-AC7` (Transient sensitivity regression from Sidak multiple-testing penalty):** Halving the STFT per-family FAP from `1.0%` to `0.501%` to accommodate the persistent broadband channel within a single joint FWER budget caused two confirmed O1/O2 events (`GW150914 L1` and `GW170814 L1`) to fall below the stricter STFT threshold (`Confirmed Event = NO`) and reduced post-recovery synthetic chirp full-timeline efficiency to `50.0%`. This empirically demonstrates why transient (`MultiResolutionDetector` at full `1%` FAP) and persistent (`FrozenReferenceBroadbandChannel` + `CrossCoherence`) searches should be executed as separate, dedicated search modes rather than penalizing transient sensitivity for a persistent channel.
3. **`P2b-AC7` (`GW170729 V1` & `GW170818 H1`):** Orthogonal diagnostics cleanly separated Virgo impulsive glitch storms (`GW170729 V1`: `27` `GLITCH_STORM` windows, FAR reduced from `4.14%` in Phase 2 to `2.58%` [90% CI `1.16%–5.62%`]) from Hanford continuum drift (`GW170818 H1`: `120` `CONTINUUM_DRIFT` windows, `0` `GLITCH_STORM` windows, FAR reduced from `4.97%` to `4.30%` [90% CI `2.76%–6.67%`]), though `GW170729 V1` remained slightly above the `2.0%` target.
4. **Broadband Channel Under Drifting Noise (Suite B):** Under `drifting_1.0_to_1.4x` noise, `FrozenReferenceBroadbandChannel` fired on `69–75%` of background windows, confirming empirically that single-channel frozen-reference broadband excess is physically degenerate with broadband background drift and requires multi-detector cross-coherence (Suite C) for disambiguation.

