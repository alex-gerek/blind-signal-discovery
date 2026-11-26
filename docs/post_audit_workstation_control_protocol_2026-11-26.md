# Post-audit Python Workstation control protocol (2026-11-26)

Status: specified before execution. This is a bounded synthetic control study, not a confirmatory replacement for the historical 100-calibration study or real-data validation.

## Revision and preservation

- Detector baseline: development/v0.3 `043f0eb95be8eeae65b9a5f5bb63015d93559b73`, identical source tree to main merge `5e19ae4`.
- Do not tune detector, calibration, veto, generator or UI code. Record the full commit, source hashes, driver/protocol hashes and uncommitted file list.
- Start with `.venv/Scripts/python.exe -B experiments/run_post_audit_workstation_controls.py`. Each invocation allocates a fresh UTC/UUID directory under `results/post_audit_workstation_controls/`, saves this protocol and driver before execution, and flushes case outcomes as they finish. Failures and interrupted partial runs remain evidence; no replacement seeds.
- No external inputs, network calls or reserved-manifest reads. Install a Python audit hook that rejects reserved-manifest paths and all network socket connections. Synthetic array fingerprints use canonical little-endian float64 C-order bytes, including separate background and target hashes; do not mistake the self-test's legacy synthetic source labels for cryptographic hashes.
- Source changes during execution invalidate the run. The uncommitted driver and documentation are explicitly identified, not represented as a clean released revision.

## A. Stationary self-test matrix

- Seeds: 2026, 2027, 2028. All seeds included.
- Existing `run_synthetic_self_test`: transient, causal_trailing, FAP=0.01, calibration_seed=0, vetoes off; 220 background windows, 200 initial pure-noise controls, 10 linear-chirp injections, pink_colored noise, fs=1024 Hz, 4-second windows.
- Root-energy SNR: 8 and 18. Within a seed these use the same noise and clean-waveform draws with changed signal amplitude; verify hashes/pairing. The paired 200-noise controls are counted ONCE per seed in aggregate, not once per SNR.
- Report individual and pooled available false-alarm fractions, noise exposure and abstentions; per-SNR detected/available/missed/abstained injection counts. Do not pool 3 frozen calibrations as independent Bernoulli trials with one common rate or claim a precision population bound from 600 windows.
- Expected qualitative direction: stronger injections should usually be easier; do not impose monotonicity as a software invariant. All failures remain visible. No minimum recall is selected after seeing results.
- Preserve Workstation JSON, timeline CSV, candidates CSV, full configuration, array hashes, and candidate/decision consistency checks.

## B. Adaptive recovery and resource control

- Base seeds 601, 611, 621; independent NumPy streams seed, seed+1, seed+2 for reference, target noise, and injections. Inputs are repeated identically at FAP=0.01 and 0.005.
- Reference: 320 independent gaussian_white windows, sigma=1. Target: 350 independent gaussian_white windows at sigma=4, then 20 further windows with 10 linear chirps at root-energy SNR=30 relative to sigma=4, in alternating windows.
- Existing MultiResolutionDetector, FrozenReferenceBroadbandChannel, calibrate_separate_search_modes, AdaptiveRegimeManagerV2b; causal_trailing/transient, provisional=60, mature=200, max buffer=2000; per-search FAP and union FWER all equal the tested FAP.
- Resource-control companion at FAP=0.005 only: same reference and first 180 target noise windows, max buffer=150. This companion must expose unavailable calibration/resource blocking; it is an expected unavailable result, not a sensitivity success.
- Save per-window used thresholds, availability, calibration-switch indicators, regime/union states, and candidates. Measure first unavailable window, first available new-regime decision, first mature new-regime decision, availability/exposure, noise false alarms, diagnostic-candidate preservation and actual post-shift injection detections.
- Report both dedicated STFT and Union availability; do not claim Union recovery where only one channel is usable. Recovery timing is measured in windows relative to the shift at index 0, not a universal latency.
- Compatibility expectation: recover a finite dedicated STFT threshold and detect at least one strong injection; report any failure without tuning. Resource-control outputs must not accept unavailable windows as detections.

## Checks and interpretation fixed before execution

- Strict JSON (no NaN/Infinity); exported full-precision decisions match availability and the threshold actually used; unavailable exposure is excluded; counts partition detected/missed/abstained outcomes.
- A crossing after vetoes-off must agree with the active finite threshold. Exports and in-memory records must agree; preserve diagnostics even when decisions abstain.
- Save elapsed time and exceptions per case; continue independent cases after a failure. Status is COMPLETE only when all planned cases finish, with errors and compatibility flags reported separately.
- No general real-world FAP, universal recovery sensitivity, discovery, novelty or astrophysical-origin claim follows.
- Historical adaptive/coincidence studies remain historical. Next scientific rerun priority is adaptive exposure/sensitivity, then valid time-slide exposure and matching event definitions, then the previously used O2 GW170814 epoch. O3 is excluded.

## Manual operator follow-up (not satisfied by this headless study)

Use Python Local Measurement Workstation, transient / causal_trailing, FAP=0.01, calibration_seed=0, 220 background / 200 noise controls / 10 linear_chirp injections, pink_colored, seed=2026. Run SNR=8 and 18 separately. Save each JSON/CSV bundle to its own folder; capture calibration status, available exposure, false alarms and injection outcomes. Inspect one false alarm if present, one detected injection, and one missed or abstained injection if present. Record exact source revision, browser/Streamlit version, screenshots and bundle hashes. Compare against case seed2026_snr8/18; do not reuse the old 80-noise-window manual reference, whose target generation differs.

Native browser WAV decoder/render/latency validation needs a separate browser protocol and remains outstanding.
