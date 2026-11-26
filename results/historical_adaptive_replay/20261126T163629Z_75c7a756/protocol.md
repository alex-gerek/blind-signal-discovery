# Historical adaptive Phase 2b synthetic replay protocol (2026-11-26)

Status: specified before execution. Owner requested replay and assessment of possible new observations. No scientific discovery is assumed.

## Revision, scope and safety

- Current detector source: development/v0.3 043f0eb, tree equal to published main merge 5e19ae4. New driver/docs are uncommitted and must be listed and hashed.
- Historical comparator: results/v03_phase2b_corrected, introduced by edd1ec9. Preserve all earlier results; hash comparator reports/tables before/after.
- Replay ONLY synthetic suites A and B from experiments/run_v03_phase2b_corrected_benchmark.py. Do not call its main, suite D, historical figure/report writers or O3-check helper.
- Compile only five explicitly allowlisted function definitions from the original AST: _gen_obs_list, _make_chirp_chunks, _wilson_ci90, _run_suite_a_corrected, _run_suite_b_corrected. Supply current detector APIs. Manager subclass forwards calls unchanged and records their inputs/outputs; no threshold, generation or decision changes.
- Verify generator/noise source and synthetic function ASTs match edd1ec9. That supports recipe-level comparability, not bytewise identity with historical arrays: the old run lacks saved array hashes and complete environment provenance.
- No external data or O1/O2/O3 input. Audit hook blocks reserved-manifest access, data-directory opens, socket connections and writes to historical result directories. A fresh UTC/UUID directory under results/historical_adaptive_replay holds each run and partial outcomes. Do not overwrite prior results.

## Fixed cases

Suite A uses make_split_rng(3026, "v03_phase2b_suite_a"), chronological original RNG draw order and a shared 320-window pink reference. All search FAPs / union FWER=0.01, calibration_seed=3026, norm_fraction=.30, causal_trailing, transient, default manager configuration; no vetoes.

1. Stationary pink null: 200 windows.
2. Abrupt step: 100 pink windows then 400 shifted_ar2 windows, sigma=1.5. SNR18 chirps relative to sigma1.5 at transitional indices108/125/145 and fixed post-step indices330,340,...420. All are retained.
3. Stationary chirps: 120 pink windows, SNR18 chirps at15,20,...110.
4. Glitch storm: 60 pink then60 impulsive windows; no injected scientific signal. Report the impulsive half separately from ordinary noise.

Suite B uses make_split_rng(3027, "v03_phase2b_suite_b"), the shared reference, persistent policy and the original RNG order: stationary and drifting1.0_to1.4x backgrounds crossed with80s_sinusoid,80s_multitone_comb,40s_slow_chirp,80s_broadband_stochastic (120 windows each; original durations, frequency/SNR recipes and insertion indices). These are eight cases, not additional seeds.

## Hypotheses fixed before execution

- The previous reported staged times (62/201 windows after step), clean duty88%, post-recovery5/10 chirps and zero observed mature false alarms may change with exact attainability, buffering, calibration-switch logging and exposure fixes.
- Current explicit unavailable intervals should alter valid denominators, even where a legacy full-window ratio appears similar.
- The spectral-resonance/weak-local-excess explanation for shifted_ar2 misses is a hypothesis. Preserve per-injection detection, used threshold, time/frequency and noise spectral diagnostics; do not declare this cause proven by aggregate recall.
- Persistent broadband versus noise drift remains a single-channel ambiguity. High anomaly rates under drifting background do not establish astrophysical signals.
- Glitch storm should remain unrecalibrated; verify refresh counts and regime changes rather than status labels alone.
- No seeds, signal placements, SNRs, thresholds or acceptance targets will be selected after looking at outcomes.

## Saved metrics and comparisons

- Exact hashes of reference/final target array bytes (little-endianfloat64 C order), original synthetic function AST and driver/protocol/source hashes; environment and dirty paths.
- Each case: full-precision timeline CSV, candidate CSV, complete per-window records/diagnostics, manager configuration, signal/background labels, outcomes and elapsed time.
- Primary channel-specific counts: false alarms and detections over AVAILABLE windows, abstentions separately. Give total/available/clean/mature exposure and denominators. No eligible exposure -> null.
- Separate dedicated STFT, dedicated broadband and FULL Union availability. Record partial-Union nominations diagnostically; do not present them as full-Union calibrated significance.
- For the abrupt case: record first new-regime credibility label, first available dedicated STFT, first full Union, first mature available decision, all relative to the step at index100. Preserve historical field definitions alongside corrected definitions without treating them as equivalent.
- Preserve historical-style aggregate rows from the unchanged experiment functions with the label LEGACY_AGGREGATION_ON_CURRENT_OUTPUTS; they are for comparison, not the corrected primary statistics.
- Independently rederive counts and used-threshold crossings from saved per-window records. Compare old tables and new legacy-style rows; identify definition changes and altered decisions.
- Save all failures/interruptions. COMPLETE means execution completed, not scientific acceptance.

## Interpretation and next decision

No new physical signal can be discovered in this synthetic-only replay. Classify findings as reproduced behavior, confirmed implementation/aggregation defect, new bounded diagnostic observation, or untested explanation. Any claim of methodological novelty needs separate external comparison and validation. Do not claim a historical result is unchanged without checking its dependent metrics.

A useful next experiment must follow an observed limitation, use a fresh protocol and outputs, and preserve this replay. No automatic tuning, commit, push, tag/release or Zenodo update is included.
