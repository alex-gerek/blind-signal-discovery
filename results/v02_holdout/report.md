# Blind Signal Discovery (BSD) v0.2 — Final One-Shot Holdout Evaluation

> **Protocol Compliance:** This report documents the **single locked one-shot blind comparison** of frozen `v0.1-baseline` vs. frozen `v0.2` on the 4 previously untouched GWTC-1 holdout events (`GW170729`, `GW170809`, `GW170818`, `GW170823`). No code, parameters, or thresholds were modified after inspecting these data. Because `N = 4` low-to-moderate SNR BBH events (`netSNR 10.8–12.4`) is a small sample, this serves as a **limited final holdout check**, not proof of universal improvement.

## 1. Disaggregated Independent Background False-Alarm Performance

Pooled across all holdout observing epochs (`N = 4277` 4-s background segments = `4.752 hours` livetime):

| Pipeline | Per-Window FAP (FP/N) | 95% Exact CI | Triggered Windows/hr | Total Cands | Candidates/hr |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `v0.1 Raw` | `1.61%` (69/4277) | [1.26%, 2.04%] | `14.52` | `88` | `18.52` |
| `v0.2 Pre-Veto` | `2.76%` (118/4277) | [2.29%, 3.29%] | `24.83` | `241` | `50.71` |
| `v0.2 Post-Veto` | `2.57%` (110/4277) | [2.12%, 3.09%] | `23.15` | `218` | `45.87` |

## 2. Event-by-Event Blind Recovery on the 4 Untouched Holdout Events

- **Events recovered (`>= 1` detector):** `v0.1` = **`1/4`**, `v0.2 Pre-Veto` = **`2/4`**, `v0.2 Post-Veto` = **`2/4`**.

| Event | Net SNR | Det | v0.1 Bg FAR | v0.2 Post Bg FAR | v0.1 Rec | v0.2 Pre Rec | v0.2 Post Rec | v0.2 Post Other 32s |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GW170729` | 10.8 | H1 | 0.24% | 0.24% | NO | NO | NO | 0 |
| `GW170729` | 10.8 | L1 | 1.66% | 0.95% | NO | NO | NO | 0 |
| `GW170729` | 10.8 | V1 | 7.57% | 7.57% | NO | NO | NO | 2 |
| `GW170809` | 12.4 | H1 | 0.95% | 0.47% | NO | NO | NO | 0 |
| `GW170809` | 12.4 | L1 | 1.36% | 3.74% | YES | YES | **YES** | 1 |
| `GW170809` | 12.4 | V1 | 0.24% | 0.00% | NO | NO | NO | 0 |
| `GW170818` | 11.3 | H1 | 6.16% | 15.40% | NO | NO | NO | 0 |
| `GW170818` | 11.3 | L1 | 0.24% | 0.24% | NO | NO | NO | 0 |
| `GW170818` | 11.3 | V1 | 0.71% | 0.47% | NO | NO | NO | 0 |
| `GW170823` | 11.5 | H1 | 0.71% | 0.47% | NO | NO | NO | 0 |
| `GW170823` | 11.5 | L1 | 1.18% | 1.90% | NO | YES | **YES** | 1 |

## 3. H1–L1 Coincidence & 500-Slide Background Summary (Post-Veto)

| Event | Post-Veto H1–L1 Pairs | XCorr Delay (ms) | 500-Slide Bg Coinc | Bg Coinc FAR (/hr) | Min Resolvable p |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `GW170729` | 0 | `—` | `3` | `0.0128` | `>= 0.00200` |
| `GW170809` | 0 | `—` | `20` | `0.1224` | `>= 0.00200` |
| `GW170818` | 0 | `—` | `86` | `0.3668` | `>= 0.00200` |
| `GW170823` | 0 | `—` | `7` | `0.0299` | `>= 0.00200` |

## 4. Merge Gate Assessment (`development/v0.2` -> `main`)

| Criterion | Status | Summary |
| :--- | :---: | :--- |
| **AC-1** (Synthetic 1% FAP control) | **PASSED** | `0.93%` per-window FAR (`N=1500`) |
| **AC-2** (Synthetic recall gain without regression) | **NOT PASSED** | `sinusoid +17.6 pp`, `pulse +2.1 pp`, `intermittent +4.6 pp`, but `linear_chirp -3.5 pp` & withheld `sinusoidal_fm -3.6 pp` |
| **AC-3** (Impulsive glitch rejection <15% FAR, pulse >=85%) | **PASSED** | Impulsive FAR `94.3% -> 7.7%`, broadband pulse recall `88.0%` |
| **AC-4** (LIGO artifact suppression without losing confirmed signals) | **PASSED** | `GW170608 H1` line & `GW170817 L1` glitch vetoed; `GW150914 H1/L1` & `GW170814 L1` preserved |
| **AC-5** (Real LIGO background FAR <= 1.5%) | **NOT PASSED (Holdout)** | Dev `0.58%` (`N=6412`, passed), Holdout `2.57%` (`N=4277`, inflated by nonstationary drift in `GW170818 H1` `15.40%` & `GW170729 V1` `7.57%`) |
| **Phase 3** (Sub-10 ms timing via controlled injections) | **PASSED** | `0.37–1.08 ms` MAE at SNR>=14; `GW150914` `+7.32 ms` (`r=-0.43`) |

**Recommendation:** **Do NOT merge `development/v0.2` into `main` automatically.** While v0.2 passed AC-1, AC-3, AC-4, and Phase 3 timing validation (and recovered `2/4` vs. `1/4` holdout events), **AC-2** (`-3.5 pp` / `-3.6 pp` cross-scale penalty on chirps/FM) and **AC-5 on nonstationary holdout epochs** (`2.57%` vs. `1.61%` for v0.1, driven by `GW170818 H1` background drift between early calibration and late validation windows) remain open scientific issues requiring explicit review before any merge decision.
