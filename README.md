<div align="center">

# Blind Signal Discovery (BSD)

**Finding unknown structure in noisy time series, with a false-alarm rate you can measure.**

[![Recorded Python tests](https://img.shields.io/badge/Python%20tests-87%20reported%20at%20986dd46-64748b.svg)](docs/CURRENT_LIMITATIONS.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-0f9e8e.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-1f8fd0.svg)](pyproject.toml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24557479.svg)](https://doi.org/10.5281/zenodo.24557479)

**[Open the workstation](https://alex-gerek.github.io/blind-signal-discovery/workstation/)** · [Website](https://alex-gerek.github.io/blind-signal-discovery/) · [How it is measured](#what-has-been-measured) · [Mathematical guide](docs/whisker_mathematical_measurement_guide.md) · [Audit erratum](docs/AUDIT_ERRATUM_47a5475.md) · [Research history](docs/RESEARCH_HISTORY.md)

</div>

---

Most signal detectors search for shapes they already know. This prototype asks the complementary question: **can a purely classical, template-free detector find previously unseen time-frequency structure in noisy data while keeping its false-alarm probability calibrated and independently measurable?**

> **Scientific Scope:** BSD is a methodological research prototype on synthetic and locally supplied time-series data. It strictly distinguishes *measured features*, *statistical anomalies*, and *detection candidates*. **A detected anomaly is never proof of a physical or astrophysical discovery.**

---

## Why we built this (and what comes next)

This repository is the mathematical and statistical testing ground first. Before applying template-free anomaly detection and finite-sample conformal calibration to broader acoustic or musical analysis, we verify the detector's behavior under strict synthetic controls, non-stationary noise shifts, and public LIGO/Virgo O2 time-series benchmarks.

Once the mathematical foundations, calibration limits, and failure modes are thoroughly tested and understood here, we aim to explore whether these same signal-processing and statistical tools can inform our acoustic and musical research projects — including **[Periodic Table of Sound](https://github.com/alex-gerek/periodic-table-of-sound)** and **[PHARMAKON](https://github.com/alex-gerek/PHARMAKON)**. Every cross-domain transfer will be tested separately on its own controls rather than assumed.

---

## The idea in one table

| Stage | Question asked of the signal | Method (`src/bsd/`) |
| --- | --- | --- |
| **Whitening** | What stands out above the local spectral and temporal floor? | Median DC removal, unpadded STFT, 2D robust median-filter background normalization (`detector.py`, `detector_v2.py`) |
| **Multi-scale search** | Is the excess power brief (`31 ms`), intermediate (`125 ms`), or sustained (`500 ms`)? | 3-scale `MultiResolutionDetector` with per-scale robust standardization and 8-connected clustering (`detector_v2.py`) |
| **Calibration** | How high must a score be so pure noise exceeds it with probability $\le \alpha$? | Finite-sample conformal order statistic (`conformal_order_rank`), strict `p_emp <= target_fap`, explicit `CALIBRATION_UNAVAILABLE` (`null`) when unattainable (`calibration.py`) |
| **Regime monitoring** | Has the background drifted or entered an impulsive glitch storm? | Orthogonal spectral-continuum + impulsive-crest monitors, staged causal recalibration, explicit per-channel and Union availability (`adaptive_regime_v2b.py`) |
| **Coincidence & vetoes** | Is a candidate reproducible across detectors or explained by a local artifact? | Reversible glitch/line vetoes (`vetoes.py`), peak-time coincidence, band-limited cross-correlation lag, overlap-corrected time slides (`coincidence.py`) |

---

## What has been measured

Every claim requires a control with a known answer, explicit uncertainty, and documented failure modes. The table distinguishes historical measurements from current-source checks; it does not report new reruns. See [known limitations and independent review scope](docs/CURRENT_LIMITATIONS.md).

| Study / Check | Evidence & Measured Result | Scope / Reference |
| --- | --- | --- |
| **Python regression checks** | Author-reported **`87 / 87` at `986dd46`**; independent follow-up at `e1bd520`: **70 passed, 2 deselected** in a selected subset (not a full-suite rerun). Exercised paths include exact conformal `nextafter` bounds, low-FAP & Union STFT recovery, `n_norm_samples` Split-B alignment, `max_recal_buffer_windows` validation, coincidence/time-slide exposure, strict JSON `allow_nan=False` export (scope differs between the recorded suites) | [tests/](tests/) |
| **Confirmatory synthetic negative control** | `316 / 50,000 = 0.6320%` empirical FAR across `K = 100` independent seeds (`500` windows/seed, `N_B = 154`, expected $\mu_0 = 0.6452\%$, `FAIL_TO_REJECT_H0`) | [results/workstation_confirmatory_k100/](results/workstation_confirmatory_k100/report.md) · Two-sided MC SE erratum (`0.00069` → `0.00222319`) in [AUDIT_ERRATUM_47a5475.md](docs/AUDIT_ERRATUM_47a5475.md) |
| **Manual Workstation verification (`986dd46`)** | Visually inspected synthetic self-test (`220` bg, `80` noise-only, `10` `linear_chirp` in `pink_colored`, `PROVISIONAL`, `UNION_UNAVAILABLE`, threshold `2.8889`): **`SNR=18`: `10/10` detected, `1/80` FA**; **`SNR=8`: `3/10` detected, `7/10` missed, `1/80` FA** | [manual_workstation_synthetic_verification_986dd46.md](docs/manual_workstation_synthetic_verification_986dd46.md) |
| **v0.1 baseline synthetic benchmark** | `0.80%` held-out noise FAR (`12 / 1,500`, 95% CI `[0.41%, 1.39%]`); `59.51%` known-family recall; `37.92%` withheld `sinusoidal_fm` recall (`94.4%` at SNR 25); documented failure modes under `shifted_ar2` (`3.56%`) and impulsive glitches (`84.0%`) | [results/research_v01/report.md](results/research_v01/report.md) · [scientific_review_package.md](docs/scientific_review_package.md) |
| **Historical adaptive & LIGO/Virgo O2 records** | Preserved untouched as pre-audit version records; O3 holdout (`docs/v03_o3_holdout_manifest.json`) strictly untouched | [AUDIT_ERRATUM_47a5475.md](docs/AUDIT_ERRATUM_47a5475.md) · [RESEARCH_HISTORY.md](docs/RESEARCH_HISTORY.md) |

> **Historical note:** Python regression checks cover their exercised code paths. They do not verify the separate JavaScript demo or automatically recompute or revalidate pre-audit historical reports under `results/`. See [docs/AUDIT_ERRATUM_47a5475.md](docs/AUDIT_ERRATUM_47a5475.md).

---

## Open the workstation

- **In your browser (no installation required):** **[Open the Browser Workstation](https://alex-gerek.github.io/blind-signal-discovery/workstation/)** — runs template-free STFT excess-power whitening, 8-connected clustering, exact finite-sample conformal calibration (`1/(N+1) <= alpha`), interactive Viridis spectrograms with candidate boxes, local `.wav`/`.csv` upload in browser memory, and CSV/JSON export. **File-input and reproducibility limitations are open:** WAV aliasing, externally unverified independence when timestamps are absent, abstention-summary denominators, and incomplete synthetic/export provenance. Use synthetic data for demonstration and consult [known limitations](docs/CURRENT_LIMITATIONS.md) before drawing conclusions from file measurements.
- **Input checks:** Uploads are stored in separate directories even when names match. Invalid CSV rows are rejected rather than deleted. Browser file runs compare SHA-256 hashes, effective rates and comparable timestamps; absent timestamps require external separation verification. Browser CSV/TXT supports one value column or two columns (`time_seconds,value`), optional header, and no quoted fields; Python offers explicit CSV column selection.
- **Project overview page:** [https://alex-gerek.github.io/blind-signal-discovery/](https://alex-gerek.github.io/blind-signal-discovery/)
- **Full 3-scale audited Python Workstation (`app/bsd_workstation.py`):** For the full 3-scale `MultiResolutionDetector`, adaptive regime monitoring (`AdaptiveRegimeManagerV2b`), and `.hdf5`/`.npy` support locally on `127.0.0.1:8501`:

```powershell
git clone https://github.com/alex-gerek/blind-signal-discovery.git
cd blind-signal-discovery
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[gui,ligo,dev]"
.\.venv\Scripts\streamlit.exe run app\bsd_workstation.py --server.address 127.0.0.1 --server.port 8501
```

See [docs/bsd_workstation_plan_and_contract.md](docs/bsd_workstation_plan_and_contract.md), [docs/manual_workstation_synthetic_verification_986dd46.md](docs/manual_workstation_synthetic_verification_986dd46.md), and [SECURITY.md](SECURITY.md).

---

## Quick start (CLI & tests)

Requirements: Python 3.11+ on local CPU (`numpy`, `scipy`, `matplotlib`, `pytest`, optional `streamlit` and `h5py`).

```powershell
# Run the Python test suite (the current result must be recorded with its commit)
.\.venv\Scripts\python.exe -m pytest -v

# Run fast smoke-test experiment
.\.venv\Scripts\python.exe experiments\run_experiment.py --config configs\smoke.toml
```

Python synthetic generators use named `numpy.random.SeedSequence` streams for repeatable draws. Repeatability does not prove statistical independence or data separation; the fallback mapping for arbitrary split names can collide. The browser demo uses a separate RNG. Record the implementation revision, configuration, seeds and input provenance for each study.

---

## Documentation, Audits, and Related Projects

- **Known Limitations & Independent Review Scope:** [docs/CURRENT_LIMITATIONS.md](docs/CURRENT_LIMITATIONS.md)
- **Mathematical Measurement Guide:** [docs/whisker_mathematical_measurement_guide.md](docs/whisker_mathematical_measurement_guide.md)
- **Workstation Contract & Guide:** [docs/bsd_workstation_plan_and_contract.md](docs/bsd_workstation_plan_and_contract.md)
- **Manual Workstation Verification (`986dd46`):** [docs/manual_workstation_synthetic_verification_986dd46.md](docs/manual_workstation_synthetic_verification_986dd46.md)
- **Independent Audit Erratum (`47a5475`):** [docs/AUDIT_ERRATUM_47a5475.md](docs/AUDIT_ERRATUM_47a5475.md)
- **Research History Log:** [docs/RESEARCH_HISTORY.md](docs/RESEARCH_HISTORY.md)
- **v0.1 Scientific Review Package:** [docs/scientific_review_package.md](docs/scientific_review_package.md)
- **Security Policy:** [SECURITY.md](SECURITY.md)
- **Related Projects:** [Periodic Table of Sound](https://github.com/alex-gerek/periodic-table-of-sound) · [PHARMAKON](https://github.com/alex-gerek/PHARMAKON)

---

## License & Citation

- **License:** [MIT License](LICENSE) — Copyright (c) 2026 Lex Gerek. Third-party public datasets (e.g., GWOSC / LIGO Open Science Center) remain under their respective upstream licenses.
- **Archived release citation:** [CITATION.cff](CITATION.cff) describes **v0.3.0 at `c406e00`**, archived under DOI [`10.5281/zenodo.24557479`](https://doi.org/10.5281/zenodo.24557479).

The current source may include changes after that archived snapshot. State the exact Git commit when citing later source. Version metadata use the `0.3.0` series consistently; no new release or archive is implied by these documentation corrections.
