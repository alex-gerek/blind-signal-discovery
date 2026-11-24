# Blind Signal Discovery (BSD) v0.2 — Phase 1 Report (`MultiResolutionDetector`)

> **Scope:** Phase 1 implements ONLY the 3-scale STFT bank (`31.25 ms`, `125 ms`, `500 ms` at `1024 Hz`) with split-sample per-scale normalization and joint 1% experiment-level FAP calibration. Instrumental vetoes and sub-10 ms cross-detector coincidence are intentionally excluded from Phase 1. The 4 reserved GWTC-1 holdout events (`GW170729`, `GW170809`, `GW170818`, `GW170823`) remain strictly untouched.

## 1. Implementation Summary

- **Module:** `src/bsd/detector_v2.py` (`MultiResolutionDetector`) wrapping the untouched BSD v0.1 `TemplateFreeDetector` (`src/bsd/detector.py`) at three resolutions: `short_31ms` (`nperseg=32, noverlap=28, dt=3.91 ms, df=32 Hz`), `medium_125ms` (`nperseg=128, noverlap=96, dt=31.25 ms, df=8 Hz`), and `long_500ms` (`nperseg=512, noverlap=384, dt=125 ms, df=2 Hz`).
- **Joint Calibration:** Split A estimates per-scale null location (`median`) and upper-tail spread (`Q_0.95 - median`) to standardize scores before taking the cross-scale maximum $T_{\text{global}} = \max_m \tilde{S}_m$. Split B (`N_thresh = 1000` per seed in synthetic benchmark) calibrates the single joint 1% FAP threshold via the conservative order statistic $k = \lceil 0.99 (N_B + 1) \rceil$.

## 2. Synthetic Head-to-Head Benchmark (3 Seeds, Exact Divisible Stratification)

### 2.1 Disaggregated False-Alarm Performance (Target FAP = 1.00%)

| Noise Condition | Detector | Per-Window FAP (FP/N) | 95% Exact CI | Candidate FAR (/hour) |
| :--- | :--- | :---: | :---: | :---: |
| In-Distribution (`N=1500`) | `v0.1` | `1.00%` (15/1500) | [0.56%, 1.64%] | `9.00` |
| In-Distribution (`N=1500`) | `v0.2 Phase 1` | `0.93%` (14/1500) | [0.51%, 1.56%] | `8.40` |
| Shifted AR(2) (`N=450`) | `v0.1` | `4.22%` (19/450) | [2.56%, 6.52%] | `40.00` |
| Shifted AR(2) (`N=450`) | `v0.2 Phase 1` | `5.11%` (23/450) | [3.27%, 7.57%] | `50.00` |

### 2.2 Signal Detection Sensitivity (Known & Withheld Families, SNR 4–25)

| Signal Family | Status | v0.1 Recall | v0.2 Phase 1 Recall | Change (Δ pp) | Winning Scale Breakdown (Short / Med / Long) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `sinusoid` | Known | `44.21%` (191/432) | `61.81%` (267/432) | **`+17.59 pp`** | `0 / 47 / 220` |
| `broadband_pulse` | Known | `79.86%` (345/432) | `81.94%` (354/432) | **`+2.08 pp`** | `156 / 196 / 2` |
| `linear_chirp` | Known | `59.49%` (257/432) | `56.02%` (242/432) | **`-3.47 pp`** | `0 / 223 / 19` |
| `intermittent_narrowband` | Known | `59.26%` (256/432) | `63.89%` (276/432) | **`+4.63 pp`** | `0 / 162 / 114` |
| `sinusoidal_fm` | **Withheld** | `39.93%` (230/576) | `36.28%` (209/576) | **`-3.65 pp`** | `30 / 171 / 8` |
| **All Known Pooled** | Known | `60.71%` (1049/1728) | `65.91%` (1139/1728) | **`+5.21 pp`** | — |

## 3. Previously Examined LIGO Development Data Comparison (No Holdout Touched)

- **Pooled Independent Background FAR (`N=6412` 4-s windows):** `v0.1` = **`0.56%`** (`36/6412`, `5.33` cands/hr) vs. `v0.2 Phase 1` = **`0.69%`** (`44/6412`, `9.26` cands/hr).

| Event | Det | Net SNR | v0.1 Bg FAR | v0.2 Bg FAR | v0.1 Rec | v0.2 Rec | v0.2 Winning Scale | v0.2 Timing Err (s) | v0.2 Other 32s Cands |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914` | H1 | 25.2 | 0.47% | 1.42% | YES | **YES** | `medium_125ms` | `-0.0250` | 0 |
| `GW150914` | L1 | 25.2 | 0.95% | 1.18% | YES | **YES** | `short_31ms` | `+0.0062` | 0 |
| `GW151012` | H1 | 10.0 | 0.24% | 0.24% | NO | NO | `—` | `N/A` | 0 |
| `GW151012` | L1 | 10.0 | 0.71% | 0.47% | NO | NO | `—` | `N/A` | 0 |
| `GW151226` | H1 | 13.1 | 0.00% | 0.00% | NO | NO | `—` | `N/A` | 0 |
| `GW151226` | L1 | 13.1 | 0.47% | 0.71% | NO | NO | `—` | `N/A` | 0 |
| `GW170104` | H1 | 13.0 | 0.47% | 0.24% | NO | NO | `—` | `N/A` | 0 |
| `GW170104` | L1 | 13.0 | 1.90% | 1.66% | NO | NO | `—` | `N/A` | 0 |
| `GW170608` | H1 | 15.4 | 0.24% | 0.24% | YES | **YES** | `long_500ms` | `-2.2500` | 1 |
| `GW170608` | L1 | 15.4 | 0.00% | 0.00% | NO | NO | `—` | `N/A` | 0 |
| `GW170814` | H1 | 17.2 | 0.24% | 0.24% | NO | NO | `—` | `N/A` | 0 |
| `GW170814` | L1 | 17.2 | 0.24% | 0.47% | YES | **YES** | `medium_125ms` | `-0.0312` | 4 |
| `GW170814` | V1 | 17.2 | 0.00% | 0.00% | NO | NO | `—` | `N/A` | 0 |
| `GW170817` | H1 | 33.0 | 0.00% | 0.00% | NO | NO | `—` | `N/A` | 0 |
| `GW170817` | L1 | 33.0 | 0.51% | 0.51% | NO | **YES** | `long_500ms` | `-1.0250` | 5 |
| `GW170817` | V1 | 33.0 | 2.31% | 3.33% | NO | NO | `—` | `N/A` | 0 |

## 4. Honest Assessment: Improvements, Regressions, and Limitations

- **Measured Improvements:**
  - Pooled known-signal recall rose **`+5.21 pp`** (`60.71% -> 65.91%`) at matched 1% FAP (`0.93%` vs. `1.00%` per-window in-distribution FAR; `8.40` vs. `9.00` cands/hr).
  - Large gain on stationary `sinusoid` (**`+17.59 pp`**, `44.21% -> 61.81%`, driven by `long_500ms`), plus gains on `intermittent_narrowband` (**`+4.63 pp`**) and `broadband_pulse` (**`+2.08 pp`**, driven by `short_31ms`).
  - On `GW150914 L1`, `short_31ms` achieved `+0.0062 s` peak timing error relative to the published reference (pending controlled time-shift injection validation in Phase 3).
- **Measured Regressions (Reported Transparently):**
  - `linear_chirp` regressed by **`-3.47 pp`** (`59.49% -> 56.02%`) and withheld `sinusoidal_fm` regressed by **`-3.65 pp`** (`39.93% -> 36.28%`). Both morphologies were already near-optimally matched to the single `125 ms` scale, so the cross-scale multiple-comparison penalty slightly raised the effective `125 ms` threshold. Thus, **AC-2 from the v0.2 validation plan is NOT yet satisfied for `sinusoidal_fm` and `linear_chirp`**.
  - Unseen `shifted_ar2` FAR increased slightly (`4.22% -> 5.11%`).
  - On real LIGO development data, split-calibration reduced `N_thresh` from `500` to `350`, increasing finite-sample threshold variance; `GW170608 H1` (`timing_err = -2.25 s`) and `GW170817 L1` (`timing_err = -1.025 s`, `5` other triggers) remain unvetoed instrumental artifacts rather than clean astrophysical recoveries.
- **Remaining Limitations (Scheduled for Phase 2 & Phase 3):**
  1. No glitch or persistent-line vetoes yet (`GW170817 L1` glitch and `GW170608 H1` calibration line).
  2. Cross-scale multiple-comparison penalty on `125 ms`-optimal chirps/FM signals needs careful investigation (e.g., tighter scale weighting or larger `N_cal` / smoother per-scale tail standardization).
  3. Sub-10 ms timing accuracy has not yet been validated via controlled time-shift injections (scheduled for Phase 3).

## 5. Recommended Next Step

Proceed to **BSD v0.2 Phase 2 (`src/bsd/vetoes.py`)**: implement auditable, reversible impulsive-glitch and persistent-line vetoes preserving all raw candidates, and report both pre-veto and post-veto performance.

## 6. Reproducibility Commands

```powershell
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe experiments\run_v02_phase1_benchmark.py
```
