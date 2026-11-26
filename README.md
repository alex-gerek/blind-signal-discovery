# Blind Signal Discovery (BSD)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24557479.svg)](https://doi.org/10.5281/zenodo.24557479)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A standalone, classical signal-processing and statistical research prototype for detecting statistically significant unknown signals in noisy time-series data without signal templates or machine learning.

> **Current Release (`v0.3.0`):** Audited multi-resolution detector, adaptive regime monitoring, local measurement Workstation (`app/bsd_workstation.py`), historical erratum ([docs/AUDIT_ERRATUM_47a5475.md](docs/AUDIT_ERRATUM_47a5475.md)), research history ([docs/RESEARCH_HISTORY.md](docs/RESEARCH_HISTORY.md)), and manual synthetic verification ([docs/manual_workstation_synthetic_verification_986dd46.md](docs/manual_workstation_synthetic_verification_986dd46.md)). Licensed under the [MIT License](LICENSE); see [SECURITY.md](SECURITY.md) for security guidance.

> **Scientific Scope:** This project is a methodological signal-detection study on simulated time-series data. It investigates template-free anomaly detection under calibrated false-alarm control. **A detected anomaly is not a physical discovery.** No claims are made regarding dark matter, axions, gravitational waves, or any physical phenomenon.

## Primary Research Question

> Can a template-free anomaly detector discover previously unseen signal structures in noisy data while maintaining a controlled and independently measured experiment-level false-alarm rate?

## Documentation and Results

- **Full Experimental Research Report (`research_v01`):** [results/research_v01/report.md](results/research_v01/report.md)
- **Scientific Review Package & Statistical Audit:** [docs/scientific_review_package.md](docs/scientific_review_package.md)
- **Machine-Readable Outputs:** [results.json](results/research_v01/results.json) | [summary_by_family_snr.csv](results/research_v01/summary_by_family_snr.csv) | [per_sample_results.csv](results/research_v01/per_sample_results.csv)

### Summary of Measured v0.1 Results (3 Independent Seeds: `101, 202, 303`)

- **Held-out in-distribution noise-only FAR:** `0.80%` (`12 / 1500`; 95% exact Clopper-Pearson CI `[0.41%, 1.39%]`) against a calibrated `1.0%` experiment-level target FAP.
- **Held-out known-signal detection (SNR 4–25):** `59.51%` recall (`914 / 1536`), `98.70%` precision, `99.45%` time-frequency box overlap on detected events.
- **Withheld signal family (`sinusoidal_fm`, unseen during development):** `37.92%` pooled recall (`182 / 480`), rising to `90.7%` at SNR 20 and `94.4%` at SNR 25.
- **Documented failure modes (automatically flagged as statistically non-credible):** `3.56%` FAR under unseen resonant AR(2) noise (`shifted_ar2`) and `84.00%` FAR under unmodeled impulsive interference.

### Selected Figures (Simulated Synthetic Data)

- [Time-Domain Signal Examples](results/research_v01/figures/fig01_time_domain_examples.png)
- [Spectrograms and Excess-Power Maps](results/research_v01/figures/fig02_spectrograms_and_excess_maps.png)
- [Candidate Overlays (Known and Withheld Families)](results/research_v01/figures/fig03_candidate_overlays.png)
- [Detection Probability vs. Energy SNR](results/research_v01/figures/fig04_detection_vs_snr.png)
- [Empirical False-Alarm Rates Across Noise Conditions](results/research_v01/figures/fig05_false_alarm_results.png)
- [Held-Out Confusion Matrix](results/research_v01/figures/fig06_confusion_matrix.png)
- [False-Positive Diagnostic Examples](results/research_v01/figures/fig07_false_positive_examples.png)
- [Missed-Signal (False-Negative) Diagnostic Examples](results/research_v01/figures/fig08_missed_signal_examples.png)
- [Calibration Null Distribution](results/research_v01/figures/fig09_calibration_null_distribution.png)

## Architecture

The pipeline strictly separates synthetic data generation, blind anomaly detection, statistical calibration, and held-out evaluation:

1. **Noise Models (`src/bsd/noise.py`)**: Gaussian white noise, pink-like $1/f^\alpha$ colored noise, smoothly nonstationary noise, impulsive interference, and an unseen resonant AR(2) shifted-noise distribution.
2. **Synthetic Data Generator (`src/bsd/generator.py`)**: Generates isolated `Observation` objects (raw 1D time series and sample rate only) and separate `GroundTruth` records across four known signal families (`sinusoid`, `broadband_pulse`, `linear_chirp`, `intermittent_narrowband`) plus a withheld blind-test-only family (`sinusoidal_fm`).
3. **Template-Free Baseline Detector (`src/bsd/detector.py`)**: Robust median DC removal, unpadded STFT power spectrogram (excluding DC and Nyquist bins), 2D robust median-filter spectral/temporal background whitening, normalized excess-power scoring, and 8-connected-component clustering.
4. **Statistical Calibration (`src/bsd/calibration.py`)**: Calibrates the experiment-level maximum cluster anomaly score on an independent noise-only calibration split using the finite-sample conservative order statistic $k = \lceil (1 - \alpha)(N_{\text{cal}} + 1) \rceil$ targeting $\alpha = 1\%$ experiment-level false-alarm probability (FAP), accounting for multiple comparisons across all time-frequency bins and STFT overlap correlations.
5. **Independent Blind Evaluation (`src/bsd/evaluation.py`)**: Evaluates frozen thresholds across strictly separated splits (in-distribution noise-only, known signal families, withheld signal family, shifted noise, and a 6-part robustness stress suite) with exact Clopper-Pearson and Wilson 95% confidence intervals and automatic statistical credibility checks.
6. **Visualization & Reporting (`src/bsd/visualization.py`, `src/bsd/reporting.py`)**: Produces 9 publication-quality figures, machine-readable JSON/CSV summaries, and a structured Markdown research report.

## Repository Layout

```text
blind-signal-discovery/
├── README.md
├── pyproject.toml
├── .gitignore
├── docs/
│   └── scientific_review_package.md
├── configs/
│   ├── smoke.toml          # Fast smoke-test configuration
│   └── research.toml       # Full 3-seed research configuration (>=500 noise + >=500 signal/seed)
├── src/
│   └── bsd/
│       ├── __init__.py
│       ├── noise.py
│       ├── generator.py
│       ├── detector.py
│       ├── calibration.py
│       ├── evaluation.py
│       ├── visualization.py
│       └── reporting.py
├── experiments/
│   └── run_experiment.py
├── tests/
│   └── test_bsd.py
└── results/
    └── research_v01/       # Generated JSON, CSV, figures, and research report
```

## Setup and Reproducibility

Requirements: Python 3.11+ on local CPU, NumPy, SciPy, Matplotlib, pytest.

```powershell
# Create and activate virtual environment inside the project
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"

# Run automated test suite
.\.venv\Scripts\python.exe -m pytest -v

# Run fast smoke-test configuration
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\smoke.toml

# Run full 3-seed research configuration
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\research.toml
```

All random number streams are derived deterministically via `numpy.random.SeedSequence([master_seed, split_stream_code])`, guaranteeing exact reproducibility and statistical independence across calibration, held-out test, and robustness splits.
