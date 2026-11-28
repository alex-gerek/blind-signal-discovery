# Matched gap-2 Hz research publication — 28 November 2026

The owner approved committing and pushing the completed matched experiment and saved-map diagnosis to development/v0.3. This is evidence preservation, not detector adoption, a main merge, a release or a new Zenodo archive. Earlier local/pending-publication statements describe their historical state before this approval.

- [Experiment findings](MATCHED_GAP2_FINDINGS_2026-11-28.md): four fixed arms, independent calibration/negative controls, all matched cases and false/wrong/no-pair outcomes.
- [Decision-stage diagnosis](MATCHED_GAP2_DIAGNOSIS_2026-11-28.md):24strict successes,5wrong pairs,187no-pair gates; repeated139Hzpeak and temporal averaging observations.
- [Knowledge map](what_is_known.md): BSD-23/24, including limits and conditions for further work.

## Public versus local evidence

The results/matched_gap2/20261128T115035Z_bf957291/ directory contains the frozen protocol, archived driver/dependency, calibration/raw-score data, frozen thresholds, evaluation records, paired contrasts, verification and report. Three large records are saved directly as .json.gz; their scripts read gzip natively. Results/matched_gap2_diagnosis/20261128T130736Z_b2471926/ contains the diagnostic script, all case tables, profile/gate accounting, control decomposition, figures and reviewed report.

The504compressed NPZ traces (~165.12MiB) remain persistent locally under the experiment's traces/ folder and are excluded by the existing *.npz rule. Exact paths and SHA256 are listed in its artifact_sha256.json. They are NOT backed up by this GitHub push. The public checkout retains all numerical decisions and profile diagnostics; map-level verification additionally needs those traces or a separately recorded deterministic reproduction. Do not describe this as a remote backup of every binary artifact.

## Reproduction and scope

The measurement source pin is fc99a63915cc5137976ea493df1fb61fe6d803d6. Use a separate disposable checkout at that revision, copy the archived/new matched driver and frozen protocol into their experiments/ and docs/ locations, and retain the unchanged blind-pilot dependency already present at that pin. The driver deliberately rejects a different HEAD. Run with the project's recorded Python/dependencies into a NEW output folder; its own snapshots, seeds and metadata preserve the trial. A reproduction is a new research step, not permission to overwrite the saved run.

Verification/report scripts accept the run folder. The forensic diagnostic script explicitly identifies the original parent run and requires its local NPZ maps; analysis of a newly reproduced run needs an explicitly recorded parent-path configuration change. Preserve the reviewed narrative reports rather than silently replacing them with generic builders.

The542saved artifact hashes were checked before staging. Scoped arithmetic/map checks passed earlier; no full regression suite or waveform experiment was rerun just for this publication. Production source, historical outputs, unrelated test cleanup and the article draft are excluded. No O3 data/manifest access. The private cross-project register retains the same step summaries locally; it is not exported.
