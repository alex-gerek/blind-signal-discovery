# BSD Local Measurement Workstation — Repository Audit, Implementation Plan & Interface Contract

**Status:** Historical approved design, retained below; implementation status clarified on 2026-11-26. Parts 1–3 describe the original planning state, not an up-to-date feature inventory.
**Original implementation branch:** `development/v0.3`; the implementation has since been merged into `main`. This historical note does not prescribe current branch-protection settings.
**Core Principle:** Wrap existing BSD Python modules (`src/bsd/`) through a clean, minimal local Streamlit GUI without rewriting detector mathematics, without hardcoding LIGO assumptions on arbitrary data, and without touching the reserved O3 holdout (`docs/v03_o3_holdout_manifest.json`).

---

## Current implementation summary (2026-11-26)

The Python loader now supports WAV, CSV, NPY and optional GWOSC-schema HDF5. The default analysis is 1024 Hz with 4-second windows; higher sampling rates require explicit antialiased resampling opt-in. Background size must support the requested FAP after the actual Split A/B allocation: 120 total windows do not support a 1% target. Initial unattainable calibration is rejected. Adaptive availability and diagnostic nominations are recorded separately from accepted decisions.

The browser demo is a separate single-scale implementation, outside the original zero-duplication Python-wrapper contract. Its data-loading, independence and export limitations are open. The Python upload-basename collision, abstention-summary denominators and rounded exports also remain open. See the [current mathematical guide](whisker_mathematical_measurement_guide.md) and [known limitations](CURRENT_LIMITATIONS.md). Historical thresholds and performance numbers below have not been rerun.

## Part 1 — Historical Repository Audit (Capabilities at Planning Time)

### 1. Loadable File Formats Today
- **GWOSC HDF5 (`.hdf5`) ONLY** via `bsd.ligo_adapter.load_gwosc_hdf5(path)` (`src/bsd/ligo_adapter.py`): strictly expects official GWOSC schema (`/strain/Strain` with `Xspacing`, `Xstart`, `/quality/simple/DQmask`, `/meta`).
- **In-memory synthetic streams** via `bsd.generator` and `bsd.noise`.
- **Not yet supported natively:** generic `.csv`, `.npy`/`.npz`, `.wav`, or arbitrary HDF5 files. The workstation will add a thin, strictly validated, domain-agnostic loader (`src/bsd/io_loader.py`) using already-installed `numpy` and `scipy.io.wavfile` (zero exotic audio dependencies).

### 2. Calibration and Detection Functions / Classes
| Module | Classes / Functions | Role & Maturity |
| :--- | :--- | :--- |
| `bsd.detector` + `bsd.calibration` | `TemplateFreeDetector`, `calibrate_detector` | **Frozen v0.1 baseline** (single-scale STFT excess power + conformal calibration). |
| `bsd.detector_v2` | `MultiResolutionDetector`, `.calibrate()`, `.analyze()` | **v0.2 evaluated prototype** (3-scale STFT: `31.25 / 125 / 500 ms` at `1024 Hz`, split-sample joint FAP calibration). |
| `bsd.vetoes` | `VetoConfig`, `apply_vetoes()` | Auditable, reversible glitch/line veto annotations (never deletes raw candidates). |
| `bsd.adaptive_regime_v2b` | `FrozenReferenceBroadbandChannel`, `calibrate_separate_search_modes()`, `fit_glitch_reference()`, `AdaptiveRegimeManagerV2b` | **v0.3 Phase 2b (corrected `edd1ec9`)** — dedicated `Transient Search` (`1%` FAP), dedicated `Persistent Search` (`1%` FAP), explicit Sidak `Combined Union`, dependence-aware Beta uncertainty (`N_eff`, Ljung-Box $p_{\text{LB}}$), orthogonal `continuum_monitor` & `glitch_storm_monitor`. |
| `bsd.coincidence` & `bsd.cross_coherence` | `evaluate_coincidence()`, `evaluate_cross_coherence()` | Two-detector time-slide prototypes (outside single-channel workstation v1 scope). |

### 3. Input Data Structure, Sample Rate, and Metadata Requirements
- **Core detector input:** `Observation(sample_id: str, signal: np.ndarray [1D float64], sample_rate: float, duration: float)`.
- **Sample-rate & window constraints (derived directly from existing modules):**
  - `MultiResolutionDetector` (`nperseg = 32, 128, 512`) and `AdaptiveRegimeManagerV2b` (`stft_nperseg = 512`, `DEFAULT_FREQ_BANDS_HZ = 20–80, 80–160, 160–280, 280–450 Hz`) require **Nyquist frequency `>= 450 Hz`** (canonical **`sample_rate = 1024.0 Hz`**) and **`window_duration_s = 4.0 s`** (`4096` samples/window).
  - **No silent LIGO assumptions:** Arbitrary user files are **never** passed through `adapted_ligo_whitened` (LIGO `30–380 Hz` bandpass). Files at `fs == 1024 Hz` run directly; files with `fs > 1024 Hz` require explicit user opt-in to anti-aliased decimation/resampling to `1024 Hz` (with Nyquist `512 Hz` disclosed); files with `fs < 1024 Hz`, non-finite values, or `< 60` background windows are **rejected with an explicit explanation**.
- **Separate background/calibration dataset requirement:** `MultiResolutionDetector.calibrate` requires **`>= 60` independent background windows** (`>= 240 s` at `4 s/window`); `>= 200` windows (`>= 800 s`) are recommended for `CALIBRATED_EXACT` / tight empirical FAR bounds.

### 4. Required Parameters
- **Data inputs:** Target time-series file + separate Background/Calibration time-series file (or Synthetic Injection Test mode with configurable noise model, signal family, SNR, and RNG seed).
- **Search mode:** `Transient Search` (dedicated `1%` STFT FAP) or `Persistent Search` (dedicated `1%` Frozen-Reference Broadband FAP + `policy=persistent`, zero single-channel recalibration). Optional read-only `Sidak Combined Union` display with explicit multiple-testing disclosure.
- **Processing mode:** `causal_trailing` (zero future data) or `offline_symmetric` (retrospective, explicitly labeled `future_data_used = True`).
- **Statistical & veto settings:** `target_fap` (default `0.01`), `window_duration_s` (`4.0 s`), optional reversible `enable_glitch_veto` / `enable_line_veto` (explicitly labeled as heuristic/experimental for non-LIGO sensors), and deterministic `calibration_seed`.

### 5. Currently Available Outputs
- **Per-window dual-output stream (`DualOutputWindowV2b`):**
  - *Stream A (never suppressed):* raw & post-veto `MultiResolutionCandidate` list (`t_start, t_end, f_min, f_max, peak_time, peak_freq, raw_anomaly_score, standardized_score, empirical_p_value`), `stft_max_score`, `stft_threshold`, `broadband_excess_score`, `broadband_threshold`, separate nominations (`stft_nominated_raw`, `stft_nominated_post_veto`, `broadband_nominated`, `joint_nominated_post_veto`).
  - *Stream B (credibility & diagnostics):* `regime_id`, `credibility_status`, `continuum_status`, `glitch_storm_status`, `warning_causes`, `significance_status`, and `CalibrationUncertaintySummary` (`n_thresh`, `n_eff`, `lag1_autocorr`, `ljung_box_p`, `p_min`, `far_ci90_iid`, `far_ci90_neff`, `dependence_detected`).
- **Visualizations:** Time-domain waveform, STFT spectrogram (`matplotlib`), timeline of Stream A scores vs. threshold, and Stream B diagnostic/regime shading.

### 6. Experimental Components Not Ready for General Production Use
1. **Overall prototype status:** BSD v0.2/v0.3 is an evaluated research prototype, **not** a validated general-purpose detector; the reserved O3 final holdout has not been evaluated.
2. **Single-channel Persistent Broadband Search:** Physically degenerate with broadband background drift (`SINGLE_CHANNEL_BROADBAND_AMBIGUOUS`); cannot prove signal origin without a second detector.
3. **Post-transition recalibration & resonant noise limits:** Provisional `N = 60` recalibration carries wide FAR uncertainty (`[0.11%, 6.44%]`); nonstationary real data (`GW170818 H1`) still shows `~4.7%` surviving FAR; resonant colored noise (`shifted_ar2`) showed `50%` chirp efficiency loss (subject of the pending observability study).
4. **Heuristic vetoes (`vetoes.py`):** Developed on synthetic and LIGO-like bands; must remain optional and reversible on arbitrary sensors.
5. **Two-detector modules (`coincidence.py`, `cross_coherence.py`):** Experimental prototypes excluded from v1 single-channel GUI.

---

## Part 2 — Proposed Minimal Local Streamlit Workstation Design

### Architecture (Zero Detector-Math Duplication)
```text
<project_root>/
├── src/bsd/
│   ├── io_loader.py               # NEW: Strict format/SR/NaN/O3-guard loader (.wav, .csv, .npy, optional GWOSC .hdf5)
│   └── workstation_runner.py      # NEW: Pure-Python headless runner calling existing BSD modules + metadata/export builder
├── app/
│   └── bsd_workstation.py         # NEW: Minimal Streamlit UI (thin presentation layer over workstation_runner.py)
├── tests/
│   └── test_workstation.py        # NEW: Automated tests for io_loader, O3 blocklist, runner, synthetic test, CSV/JSON export
└── launch_workstation.bat         # NEW: One-click Windows launcher using local .venv\Scripts\streamlit.exe
```

### UI Layout (4 Focused Tabs + Left Sidebar)
- **Sidebar (Configuration & Validation):**
  1. Mode selector: **Real/Local File Analysis** vs. **Synthetic Injection Self-Test (Known Ground Truth)**.
  2. Separate file uploaders / local path selectors for **Target Data** and **Background / Calibration Data**, with instant validation card (format, sample rate, duration, window count `N_bg`, finite check, explicit rejection reasons).
  3. **Search Mode** (`Transient Search` vs. `Persistent Search`), **Processing Mode** (`causal_trailing` vs. `offline_symmetric`), `target_fap`, optional reversible vetoes, and `Run Analysis` button.
- **Tab 1 — Calibration Status & Warnings (Always Visible First):**
  - Explicit banner showing `significance_status` (`CALIBRATED_EXACT`, `EMPIRICAL_VALIDATED`, `EMPIRICAL_PROVISIONAL`, `EMPIRICAL_DEPENDENT_UNCERTAIN`, `UNCERTAIN_UNVALIDATED`), `N_thresh`, `N_eff`, Ljung-Box $p_{\text{LB}}$, 90% i.i.d. & dependence-aware FAR CIs, clean duty cycle, excluded hours, and explicit counts of `CONTINUUM_DRIFT`, `GLITCH_STORM`, and `SINGLE_CHANNEL_BROADBAND_AMBIGUOUS` warnings.
  - Permanent **Research Prototype Disclaimer** banner listing experimental caveats.
- **Tab 2 — Waveform & Spectrogram Inspector:**
  - Interactive window selector (`0 .. N_target - 1`) showing time-domain waveform, STFT spectrogram with candidate bounding boxes `(t_start..t_end, f_min..f_max)`, and full-timeline score vs. threshold + Stream B state overlay.
- **Tab 3 — Candidate Table & Ground-Truth Audit:**
  - Filterable table of **100% preserved raw candidates** (timestamp, window index, frequency band, scale, raw & standardized scores, empirical p-value, veto status, Stream B `credibility_status` & `warning_causes`).
  - In **Synthetic Injection Self-Test** mode: side-by-side ground-truth comparison (true injected windows vs. recovered candidates, empirical FAR, detection efficiency).
- **Tab 4 — Reproducible Export & Run Metadata:**
  - One-click download / local save (under `results/workstation_runs/<timestamp_run_id>/`) of:
    - `candidates.csv` (all raw + annotated candidates)
    - `timeline.csv` (per-window Stream A + Stream B summary)
    - `run_report.json` (complete run metadata: BSD package version, Git commit hash, Python/NumPy/SciPy versions, SHA-256 hashes of input files, exact config, calibration uncertainty summary, and warnings).

---

## Part 3 — Proposed Interface Contract

### 1. Input Loader Contract (`bsd.io_loader`)
```python
@dataclass(frozen=True)
class LoadedTimeSeries:
    source_name: str
    file_format: str                  # "wav" | "csv" | "npy" | "gwosc_hdf5" | "synthetic"
    sha256: str
    original_sample_rate_hz: float
    effective_sample_rate_hz: float   # Always 1024.0 Hz after optional user-approved resampling
    resampled: bool
    n_samples: int
    duration_s: float
    signal: np.ndarray                # 1D float64, strictly finite
    warnings: tuple[str, ...]         # e.g., ("Resampled from 4096 Hz -> 1024 Hz; Nyquist=512 Hz",)
```
- **Hard Validation Rules (raise `InputValidationError` with actionable message):**
  - Reject unsupported file extensions.
  - Reject `O3` holdout files (filename or GPS range matching `docs/v03_o3_holdout_manifest.json`).
  - Reject non-finite (`NaN`/`Inf`), empty, or constant-zero signals.
  - Reject `sample_rate < 1024.0 Hz` (cannot satisfy `MultiResolutionDetector` / `450 Hz` continuum bands without rewriting detector math).
  - Reject multi-channel files unless the user explicitly selects a single channel index.
  - Reject background datasets yielding `< 60` non-overlapping `4.0 s` windows.

### 2. Headless Runner Contract (`bsd.workstation_runner`)
```python
@dataclass(frozen=True)
class WorkstationConfig:
    search_mode: Literal["transient", "persistent"] = "transient"
    processing_mode: Literal["causal_trailing", "offline_symmetric"] = "causal_trailing"
    window_duration_s: float = 4.0
    target_fap: float = 0.01
    enable_glitch_veto: bool = False   # Off by default for arbitrary non-LIGO data
    enable_line_veto: bool = False     # Off by default for arbitrary non-LIGO data
    calibration_seed: int = 0
```
- Returns `WorkstationRunResult` containing:
  - `run_metadata`: `run_id`, UTC timestamp, `bsd_version`, `git_commit`, `search_mode`, `processing_mode`, `config`, input SHA-256s, experimental disclaimer flags.
  - `calibration_summary`: `n_bg_windows`, `n_thresh`, `n_eff`, `ljung_box_p`, `far_ci90_iid`, `far_ci90_neff`, `dependence_detected`, `significance_tier`.
  - `windows`: `list[DualOutputWindowV2b]` (verbatim from `AdaptiveRegimeManagerV2b`).
  - `candidates_rows`: flat list of dicts for `candidates.csv`.
  - `timeline_rows`: flat list of dicts for `timeline.csv`.
  - `ground_truth_audit`: `dict | None` (populated only in Synthetic Injection Self-Test mode).

### 3. Dependency & Launcher Plan
- Add optional `gui = ["streamlit>=1.35"]` to `pyproject.toml` and install into the existing local `.venv` (keeping all packages on drive D:).
- Provide `launch_workstation.bat` at project root: activates `.venv`, sets `TEMP`/`TMP` to `.tmp` inside the project directory, and runs `streamlit run app\bsd_workstation.py --server.headless false`.
- Automated tests (`tests/test_workstation.py`) will test `io_loader`, O3 blocklist enforcement, `workstation_runner` (Transient & Persistent modes + Synthetic Injection Self-Test), and CSV/JSON export without requiring a browser.

---

## Part 4 — Incorporated Methodological Amendments

1. **Calibration Status (`INSUFFICIENT` / `PROVISIONAL` / `VALIDATED`):** `VALIDATED` is never assigned from raw window count alone; it requires `n_thresh >= 150`, `n_eff >= 120`, `ljung_box_p >= 0.05`, and dependence-aware upper 90% Beta FAR `<= 2.5%`. Otherwise status is explicitly `PROVISIONAL` or `INSUFFICIENT` with exact unmet conditions listed.
2. **Explicit Preprocessing Provenance:** Every transformation (format, channel selection, units, original sample rate, anti-alias resampling to `1024 Hz`, effective bandwidth) is logged explicitly. LIGO-specific whitening is never applied silently to arbitrary files.
3. **Background Independence Enforcement:** Identical files (SHA-256) and overlapping time intervals (when metadata/GPS exists) are hard-rejected; unverified independence emits an explicit provenance warning.
4. **Scientific Output Semantics:** Outputs strictly distinguish *measured features*, *statistical anomalies*, *detection candidates*, and *independently confirmed physical events* (never claimed by the single-channel workstation).
5. **Reproducible Self-Test:** Synthetic injection self-test includes both positive injections with known ground truth (reporting both detected and missed injections) and a pure-noise negative control split.
6. **Safe Research Operation:** Frozen detector math and historical benchmarks remain untouched; O3 holdout files and O3 GPS ranges are hard-blocked from loading, previewing, calibrating, or analyzing.

---

## Part 5 — Windows Desktop Launch, Shutdown, and Troubleshooting

### Desktop Launch
- Double-click the **BSD Workstation** shortcut on the Windows Desktop (created via `scripts/create_desktop_shortcut.ps1`).
- The shortcut invokes `launch_workstation.bat` with the working directory set to the repository root, uses the repository-stored icon (`app/assets/bsd_workstation.ico`), activates the local `.venv`, and opens `http://127.0.0.1:8501` in the default web browser.
- **Duplicate-instance protection:** If BSD Workstation is already running on port `8501`, clicking the desktop shortcut again opens the existing session in your browser without spawning a second Streamlit server.

### Clean Shutdown
- Close the browser tab, then press `Ctrl+C` in the **BSD Local Measurement Workstation** console window (or close that console window directly). This terminates the local Streamlit server cleanly.

### Troubleshooting
- **Missing `.venv` or `streamlit.exe`:** The launcher pauses and prints the exact command to recreate `.venv` and install GUI dependencies (`pip install -e ".[gui,ligo,dev]"`).
- **Port `8501` occupied by another application:** Close the conflicting local service on port `8501` or edit `PORT=8501` in `launch_workstation.bat`.
- **Recreating the Desktop Shortcut:** Run `powershell -ExecutionPolicy Bypass -File scripts\create_desktop_shortcut.ps1` from the project root. Only the `.lnk` shortcut is placed on the Desktop; project code, `.venv`, icon assets, and data remain in their configured local locations.

