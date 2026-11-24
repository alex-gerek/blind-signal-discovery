# Blind Signal Discovery (BSD) v0.2 — Architecture & Validation Plan

> **Status:** Approved in Principle — Design & Validation Plan ONLY (no implementation started).  
> **Branch:** `development/v0.2` (branched from `main` at `12c7638`, preserving Experiments 001–003).  
> **Frozen Baseline Reference:** Annotated tag `v0.1-baseline` (`ef8d3ee`).  
> **Reserved Untouched Holdout Set:** `GW170729`, `GW170809`, `GW170818`, `GW170823` (see `results/ligo_experiment_003/v02_holdout_manifest.json`) — **must not be downloaded, inspected, or used during v0.2 development**.

---

## 1. Problem Statement & Measured v0.1 Limitations

Experiments 001–003 established a reproducible, template-free baseline (`v0.1-baseline`) that maintains ~1% experiment-level false-alarm probability on calibrated backgrounds (`0.80%` synthetic in-distribution; `0.95%` across `6,412` real LIGO/Virgo background segments) and blindly recovers high-SNR real gravitational-wave signals (`GW150914`, `GW170814 L1`, `GW170817 H1`).

However, our measurements exposed four concrete, quantitative failure modes that directly motivate BSD v0.2:

1. **Single-resolution STFT sensitivity & timing ceiling (`125 ms` window, `dt = 31.25 ms`):**
   - Weak long-duration signals fragment below `bin_seed_threshold = 4.0` (`sinusoid` pooled recall `43.98%`, withheld `sinusoidal_fm` `37.92%` in Exp 001; misses on `GW151012`, `GW151226`, `GW170104`, `GW170608 L1` in Exp 003).
   - `dt = 31.25 ms` is coarser than the `10 ms` physical light-travel delay between LIGO H1 and L1, preventing direct sub-10 ms coincidence verification.
2. **Impulsive glitch breakdown:**
   - Unmodeled 1–3 sample impulses produce broadband vertical STFT clusters, inflating synthetic FAR to `84.00%` (Exp 001) and producing a `score = 234,777` trigger on the known L1 glitch near `GW170817` (Exp 003).
3. **Persistent narrowband instrumental lines:**
   - During `GW170608`, the H1 angular-control calibration line (`~15–35 Hz`) produced a full-window horizontal cluster (`score = 4,542`) that falsely triggered the automated overlap check and inflated the H1 calibration threshold to `4,155`.
4. **Nonstationary background drift & benchmark statistical refinements:**
   - Real background FAR drifted above 1% on two nonstationary epochs (`GW170104 L1`: `4.03%`; `GW170817 V1`: `3.33%`) and on unseen resonant noise (`shifted_ar2`: `3.56%`).
   - Experiment 001 benchmark counts (`512`, `160`) were not exact multiples of stratification grid sizes (`96`, `24`), and global energy SNR differed from whitened SNR in colored noise.

---

## 2. Proposed BSD v0.2 Architecture (Classical & Template-Free)

BSD v0.2 preserves the core scientific constraint of v0.1: **classical signal processing and distribution-free statistics only — no signal templates, no machine learning, no LLM classification, strict `Observation` / `GroundTruth` isolation.**

All v0.1 code remains untouched as `TemplateFreeDetector` (`src/bsd/detector.py`) so side-by-side head-to-head comparison against `v0.1-baseline` is always possible. v0.2 will be implemented in new modules (`src/bsd/detector_v2.py`, `src/bsd/vetoes.py`, `src/bsd/coincidence.py`).

### 2.1 Multi-Resolution STFT Bank (`MultiResolutionDetector`)

Replace the single `125 ms` STFT with a 3-scale Hann STFT bank (all unpadded, DC/Nyquist excluded, at `1024 Hz`):

| Scale | `nperseg` | Window Duration | `noverlap` | Hop `dt` | Freq Bin `df` | Primary Target Morphology |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **Short** | `32` | `31.25 ms` | `28` | `3.91 ms` | `32.0 Hz` | Short BBH mergers, broadband pulses, sub-10 ms timing |
| **Medium** | `128` | `125.0 ms` | `96` | `31.25 ms` | `8.0 Hz` | Baseline v0.1 scale, moderate chirps & bursts |
| **Long** | `512` | `500.0 ms` | `384` | `125.0 ms` | `2.0 Hz` | Weak sinusoids, narrowband tones, long BNS inspirals (`GW170817`) |

**Unified Multiple-Comparison Calibration Across Scales:**
- Raw cluster sums $S_m(C)$ have different null scales across Short, Medium, and Long banks (due to different bin counts and overlap factors).
- Using a split-calibration approach on dedicated noise-only data:
  1. Fit per-scale robust location/scale normalization (median & upper-quantile scale of null max-scores $T_{m,i}$ on calibration split A) to obtain standardized scores $\tilde{S}_m(C)$.
  2. Define the experiment-level global statistic $T_{\text{global}} = \max_{m \in \{\text{short}, \text{med}, \text{long}\}} \max_C \tilde{S}_m(C)$.
  3. Calibrate a **single** experiment-level 1% FAP threshold $\hat{\tau}_{\text{global}}$ via the conservative order statistic $k = \lceil(1 - \alpha)(N_{\text{cal}} + 1)\rceil$ on independent calibration split B ($N_{\text{cal}} \ge 1000$).
- This guarantees exact finite-sample 1% experiment-level FAP control across **all three resolutions simultaneously** without Bonferroni over-conservatism.

### 2.2 Auditable Artifact Diagnostics & Reversible Vetoes (`src/bsd/vetoes.py`)

**Raw-Candidate Preservation Principle:** Every anomaly candidate exceeding the calibrated threshold **must be preserved verbatim** in the raw candidate list prior to any instrumental veto. Veto decisions are attached as auditable, reversible metadata flags (`veto_applied`, `veto_reasons`, `veto_metrics`) rather than destructively deleting candidates. All synthetic and real-data benchmarks must report detection performance, false-alarm rates, and candidate tables **both before and after vetoes**.

1. **Impulsive Glitch Diagnostic & Veto (`glitch_veto`):**
   - Computes local time-domain kurtosis, spike width (consecutive samples exceeding $5\hat{\sigma}_{\text{MAD}}$), and single-frame broadband spectral flatness in the Short STFT scale.
   - Flags/vetoes ultra-short 1–3 sample time-domain spikes with flat broadband spectra (distinguishing hardware blips from structured pulses or chirps).
   - Must be verified on synthetic broadband pulses to ensure genuine pulses are **not** vetoed.
2. **Persistent Narrowband Line Veto (`line_veto`, optional transient-search mode):**
   - Tracks frequency bins that remain persistently elevated across `> 80%` of the 4-second window *and* appear in background calibration spectra (e.g., `GW170608 H1` angular calibration line).
   - Keeps stationary-signal sensitivity configurable: disabled in general synthetic benchmarks (where `sinusoid` is a valid target), enabled in transient real-data searches (where stationary instrumental lines are known background artifacts) and explicitly documented in both cases.
3. **Adaptive Spectral Whitening Robustness:**
   - Tune/compare frequency-median kernel width on synthetic `shifted_ar2` (without touching real holdout data) to reduce resonance leakage while preserving narrowband signal sensitivity.

### 2.3 Sub-10 ms Cross-Detector Coincidence & Controlled Timing Validation (`src/bsd/coincidence.py`)

For multi-detector evaluations (H1, L1, V1):
- Use the Short STFT scale (`dt = 3.91 ms`) to narrow candidate peak timing uncertainty.
- For any H1–L1 candidate pair within `|dt| <= 0.05 s`, compute band-limited normalized time-domain cross-correlation $r(\tau)$ over $\tau \in [-15, +15]\text{ ms}$ (physical light-travel bound $\pm 10\text{ ms}$ plus margin).
- **Empirical Timing Validation Requirement:** STFT hop size (`3.91 ms`) alone is **not** evidence of sub-10 ms timing accuracy. Sub-10 ms cross-detector timing accuracy must be explicitly validated using **controlled time-shift signal injections** (known synthetic delays spanning `[-10 ms, +10 ms]` injected into independent H1/L1 noise realizations) to measure empirical timing bias, RMS error, and 95% error bounds.
- **Strengthened Time-Slide Background Analysis:** Expand unphysical background time-slide trials across all independent background validation epochs and report exact empirical upper bounds on accidental coincidence rates. **Never claim statistical significance (or false-alarm probability) finer than the resolution supported by the available independent background data and time-slide livetime** (`p >= 1 / (N_slides + 1)`).

### 2.4 Benchmark & Statistical Hygiene Improvements (Experiment 004)

- Exact divisible stratification in synthetic benchmarks (`n_signal_test = 576 = 6 × 96`, `n_withheld_test = 192 = 8 × 24`).
- Report both **global energy SNR** and **whitened matched-filter SNR** in colored noise.
- Increase synthetic calibration ensemble to `N_cal = 1200` per seed (reducing finite-sample threshold variance and lowering empirical p-value floor to `0.00083`).
- **Disaggregated False-Alarm Reporting:** Every benchmark must report three distinct false-alarm metrics separately (both before and after vetoes):
  1. **Per-window false-alarm probability** (fraction of 4-second background windows containing $\ge 1$ candidate).
  2. **Candidate false-alarm rate per hour** (total background candidates divided by total independent background livetime in hours).
  3. **Coincident-event false-alarm statistics** (empirical accidental coincidence count, rate per time-slide / per hour of equivalent time-slide livetime, and minimum resolvable p-value).

---

## 3. Three-Stage Validation Plan & Predefined Merge Acceptance Criteria

`development/v0.2` will **NOT** be merged into `main` until all three stages are completed and the predefined acceptance criteria are evaluated honestly against `v0.1-baseline`.

### Stage 1 — Synthetic Benchmark (Head-to-Head vs. `v0.1-baseline`)
Run both frozen `v0.1-baseline` and `v0.2` on identical 3-seed synthetic datasets (`seeds = [101, 202, 303]`), including known families, withheld `sinusoidal_fm`, unseen `shifted_ar2`, and impulsive stress tests.

- **AC-1 (FAP Control):** Pooled held-out in-distribution noise-only FAR must remain controlled at the 1.0% target (95% exact Clopper-Pearson CI includes `1.0%`, or point estimate $\le 1.2\%$).
- **AC-2 (Sensitivity Gain):** Measurable, statistically significant recall improvement on `sinusoid` (v0.1: `44.0%`) and withheld `sinusoidal_fm` (v0.1: `37.9%`), with no regression $> 2$ percentage points on `broadband_pulse` or `linear_chirp`.
- **AC-3 (Impulsive Glitch Rejection):** With `glitch_veto` enabled, impulsive noise-only FAR must drop from `84.0%` (v0.1) to `< 15.0%`, while genuine synthetic `broadband_pulse` recall at SNR $\ge 12$ stays $\ge 85.0\%$.

### Stage 2 — Development Real-Data Diagnostic Check (Exp 002–003 Epochs)
Evaluate frozen v0.2 on the previously examined GWOSC events (`GW150914`, `GW151012`, `GW151226`, `GW170104`, `GW170608`, `GW170814`, `GW170817`) — **without** using event timestamps during detection or tuning thresholds on target windows.

- **AC-4 (Artifact Suppression):** Eliminate the false `GW170608 H1` calibration-line trigger and cleanly flag the `GW170817 L1` glitch as an impulsive artifact, without losing blind recovery of `GW150914`, `GW170814 L1`, or `GW170817 H1`.
- **AC-5 (Real Background FAR):** Pooled empirical background FAR across all independent validation segments remains consistent with `1.0%` (95% exact CI covers `1.0%`, point estimate $\le 1.5\%$).

### Stage 3 — One-Shot Blind Comparison on the Reserved Holdout Set
**Only after all v0.2 code, hyperparameters, and veto thresholds are frozen:**
- Download and verify the 4 untouched GWTC-1 holdout events (`GW170729`, `GW170809`, `GW170818`, `GW170823`).
- Run **both** `v0.1-baseline` and frozen `v0.2` blindly in a single locked pass (reporting pre-veto and post-veto outputs, plus all three disaggregated false-alarm metrics).
- Report exact side-by-side recovery counts, timing errors, background FARs, and time-slide coincidence statistics.
- **Scope Limitation:** Treat these 4 reserved GWTC-1 events strictly as a **limited final holdout check**, not as sufficient statistical evidence of universal improvement across all astrophysical populations or observing runs.
- **No post-hoc tuning is permitted after inspecting holdout results.**

---

## 4. Merge Gate to `main`

`development/v0.2` may be merged into `main` only when:
1. All automated pytest suites (existing 27 tests + new v0.2 tests, including controlled time-shift timing validation and veto reversibility tests) pass with 0 failures.
2. Acceptance criteria **AC-1 through AC-5** are verified and documented with exact pre-veto and post-veto measurements, disaggregated false-alarm metrics, and background-bounded significance statements.
3. The one-shot blind holdout comparison (`v0.1-baseline` vs. `v0.2`) is completed and reported transparently, including any regressions or missed holdout events.
4. All original `v0.1-baseline` results, failures, scientific reports, and reproducibility records on `main` (Experiments 001–003) remain completely intact and unmodified.
