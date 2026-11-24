# Blind Signal Discovery (BSD) — Research Prototype v0.1

A standalone, classical signal-processing and statistical research prototype for detecting statistically significant unknown signals in noisy time-series data without signal templates or machine learning.

> **Scientific Scope:** This project is a methodological signal-detection study on simulated time-series data. It investigates template-free anomaly detection under calibrated false-alarm control. **A detected anomaly is not a physical discovery.** No claims are made regarding dark matter, axions, gravitational waves, or any physical phenomenon.

## Primary Research Question

> Can a template-free anomaly detector discover previously unseen signal structures in noisy data while maintaining a controlled and independently measured experiment-level false-alarm rate?

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

## Key Results

See [results/research_v01/report.md](results/research_v01/report.md) for the complete measured experimental report, failure analysis, calibration limitations, and generated figures.
