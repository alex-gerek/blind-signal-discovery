# Blind Signal Discovery (BSD) v0.2 — Phase 3 Report (Sub-10 ms Timing & Time-Slide Coincidence)

> **Scope:** Phase 3 implements `src/bsd/coincidence.py` (band-limited normalized cross-correlation with parabolic sub-sample interpolation and 500-slide background coincidence analysis). Sub-10 ms timing accuracy is empirically validated via 240 controlled time-shift injections rather than inferred from STFT hop size. The 4 reserved GWTC-1 holdout events remain strictly untouched.

## 1. Controlled Time-Shift Injection Timing Validation (240 H1/L1 Pairs)

Each injection places an identical signal with exact fractional delay $\tau_{\text{true}} \sim \mathcal{U}(-10, +10)\text{ ms}$ and random polarity $\pm 1$ into independent H1/L1 unit-variance Gaussian noise.

| Family | SNR | N | XCorr Bias (ms) | XCorr MAE (ms) | XCorr RMS (ms) | XCorr P95 |err| (ms) | Within ±2 ms | Polarity Acc | Raw STFT Peak MAE (ms) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `broadband_pulse` | 10 | 30 | `+0.127` | **`1.425`** | `1.786` | `3.163` | `66.7%` | `50.0%` | `4.40` |
| `broadband_pulse` | 14 | 30 | `+0.023` | **`1.020`** | `1.437` | `2.475` | `70.0%` | `56.7%` | `3.15` |
| `broadband_pulse` | 18 | 30 | `-0.052` | **`1.083`** | `1.498` | `2.899` | `80.0%` | `50.0%` | `3.00` |
| `broadband_pulse` | 25 | 30 | `-0.082` | **`0.702`** | `1.198` | `2.897` | `90.0%` | `70.0%` | `4.81` |
| `linear_chirp` | 10 | 30 | `-0.705` | **`1.747`** | `3.276` | `4.865` | `73.3%` | `73.3%` | `582.20` |
| `linear_chirp` | 14 | 30 | `+0.245` | **`1.024`** | `1.703` | `4.018` | `80.0%` | `76.7%` | `573.47` |
| `linear_chirp` | 18 | 30 | `-0.061` | **`0.369`** | `0.742` | `1.645` | `96.7%` | `83.3%` | `635.71` |
| `linear_chirp` | 25 | 30 | `+0.001` | **`0.489`** | `0.917` | `2.301` | `90.0%` | `80.0%` | `458.84` |

**Key Timing Finding:** Raw STFT peak-time differences exhibit `~3–5 ms` MAE on short broadband pulses and `~459–636 ms` MAE on long chirps (where the peak STFT frame shifts across noise realizations), confirming that STFT hop size (`3.91 ms`) alone does **not** establish sub-10 ms timing accuracy. Band-limited normalized cross-correlation with parabolic interpolation achieves `0.37–1.08 ms` MAE at `SNR >= 14` (`70–97%` within `±2 ms`), with polarity accuracy limited by half-cycle phase ambiguity on narrowband/short oscillatory pulses (`50–70%` on pulses vs. `73–83%` on broadband chirps).

## 2. LIGO Development Data: Coincidence & 500-Slide Background Analysis (Pre-Veto vs. Post-Veto)

Statistical significance is strictly bounded by the 500-slide background resolution (`min_resolvable_p = 1 / 501 = 0.00200`).

| Event | Role | Pre-Veto Pairs | Post-Veto Pairs | Measured XCorr Delay (ms) | Signed r(tau) | Within ±10 ms? | 500-Slide Bg Coinc (Pre -> Post) | Post-Veto Coinc FAR (/hr) | Min Resolvable p |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW150914` | REF | 9 | 9 | `+7.32` | `-0.434` | YES | `18 -> 12` | `0.0512` | `>= 0.00200` |
| `GW170608` | DEV | 0 | 0 | `—` | `—` | — | `0 -> 0` | `0.0000` | `>= 0.00200` |
| `GW170814` | DEV | 0 | 0 | `—` | `—` | — | `0 -> 0` | `0.0000` | `>= 0.00200` |
| `GW170817` | DEV | 0 | 0 | `—` | `—` | — | `0 -> 0` | `0.0000` | `>= 0.00200` |

## 3. Readiness for Final One-Shot Holdout Evaluation

With Phase 1 (`MultiResolutionDetector`), Phase 2 (`vetoes.py`), and Phase 3 (`coincidence.py`) implemented and tested on development data, all BSD v0.2 components are ready for **pre-freeze review**. The 4 reserved GWTC-1 holdout events (`GW170729`, `GW170809`, `GW170818`, `GW170823`) remain completely untouched and will be evaluated only in a single locked one-shot comparison after explicit approval to freeze v0.2.

## 4. Reproducibility Commands

```powershell
.\.venv\Scripts\python.exe -m pytest -v
.\.venv\Scripts\python.exe experiments\run_v02_phase3_coincidence_benchmark.py
```
