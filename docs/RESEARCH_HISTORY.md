# Blind Signal Discovery (BSD) — Research History Log (`development/v0.3`)

This log tracks major experimental studies, independent audits, errata, and manual verification records without modifying or overwriting historical artifacts under `results/`.

## Chronological Entries

1. **v0.1 Synthetic Baseline (`results/research_v01/`)**
   - 3-seed synthetic benchmark (`1,500` held-out noise windows; known + withheld `sinusoidal_fm` signal families) and `docs/scientific_review_package.md`.
2. **v0.2 Multi-Resolution Detector, Reversible Vetoes, Coincidence & One-Shot GWTC-1 Holdout (`results/v02_phase1/`, `results/v02_phase2/`, `results/v02_phase3/`, `results/v02_holdout/`)**
   - Development benchmarks and single locked one-shot evaluation on 4 GWTC-1 holdout events. Pre-audit code versions; see `docs/AUDIT_ERRATUM_47a5475.md`.
3. **v0.3 Non-Stationary Adaptive Regime Monitoring & Phase 2b Dual-Output Architecture (`results/v03_phase1/`, `results/v03_phase2/`, `results/v03_phase2b/`, `results/v03_phase2b_corrected/`)**
   - Staged causal recalibration and orthogonal continuum/glitch diagnostics. Pre-audit code versions; recovery, exposure, and sensitivity metrics have not been recomputed; see `docs/AUDIT_ERRATUM_47a5475.md`.
4. **Workstation Calibration & Confirmatory Negative-Control Studies (`results/workstation_calibration_study/`, `results/workstation_negative_control_2000/`, `results/workstation_multiseed_220/`, `results/workstation_confirmatory_k100/`)**
   - Includes the `K = 100, N_eval = 500` (`50,000` windows) confirmatory synthetic run (`316 / 50,000 = 0.6320%`, `FAIL_TO_REJECT_H0`). Arithmetic erratum for two-sided Monte Carlo SE (`0.00069` published -> `0.00222319` corrected) documented in `docs/AUDIT_ERRATUM_47a5475.md`.
5. **Exploratory Real-Data Search on GW170814 Epoch (`results/gw170814_exploratory_search/`, Commit `47a5475`)**
   - Exploratory O2 `H1/L1/V1` repeatability check and off-source candidate audit (`0` coincident pairs). Preserved as historical record; O3 holdout (`docs/v03_o3_holdout_manifest.json`) untouched.
6. **Multi-Round Independent Code & Methodology Audit (`986dd46`, PR #1, `docs/AUDIT_ERRATUM_47a5475.md`)**
   - Fixed exact conformal FAP boundary (`conformal_order_rank`), unattainable-FAP abstention (`threshold = None` / JSON `null`), adaptive & Union STFT recovery, `n_norm_samples` Split-B audit alignment, `max_recal_buffer_windows` validation, calibration-switch window logging, coincidence peak-time matching and `delay_identifiable` guards, overlap-corrected time-slide exposure and unavailable-significance separation, and Workstation CSV/JSON export alignment (`87 / 87` pytest passed).
7. **Manual Workstation Verification on Synthetic Data (`986dd46`, `docs/manual_workstation_synthetic_verification_986dd46.md`)**
   - Manual visual verification in the local BSD Workstation UI by the project owner, cross-checked against an in-memory computational control run on commit `986dd46` (`transient`, `causal_trailing`, `FAP = 0.01`, `220` background windows, `80` noise-only windows, `10` `linear_chirp` injections in `pink_colored` noise, `seed = 2026`, `calibration_seed = 0`, vetoes off; `PROVISIONAL` calibration, `UNION_UNAVAILABLE`, threshold `2.8889`). Recorded both `SNR = 18` (`10/10` detected, `0/10` missed, `1/80` false alarms) and `SNR = 8` (`3/10` detected, `7/10` missed, `1/80` false alarms), plus visually inspected windows `4` (false alarm), `81` (missed injection), and `83` (detected injection). No UI export bundle saved; no historical reports modified.
