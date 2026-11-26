# Verification status and historical rerun queue

Updated on 2026-11-26. This index separates saved implementation fixes, regression evidence and scientific claims. Earlier journal entries keep their original status at the time of writing; dated follow-ups supersede those status snapshots.

## Saved implementation revisions

| Step | Saved commit on `development/v0.3` | Scope and evidence |
| --- | --- | --- |
| Documentation / erratum / version alignment | `99fce40` | Corrected current guide, preserved historical guide, fixed study paths and attainability interpretation, aligned `0.3.0` metadata. Documentation/AST/rank/link checks, not a scientific rerun. |
| Upload isolation and input checks | `8053741` | Unique paths, strict real-valued CSV/NPY input, rate/hash/time checks and error handling. 34 distinct selected pytest cases across final runs, including a six-scenario Node wrapper. |
| Abstention accounting and export precision | `8fab987` | Conditional rates/exposure, separate abstentions, compatibility semantics and shared precise CSV/JSON builders. 18 selected pytest cases across final runs and ten Node scenarios. |
| Browser WAV preprocessing | `8c6d3ae` | Fixed-parameter FIR to exact 512 Hz; 0–192 Hz useful resampled band, disclosed offline lookahead/padding, original/decoded rate provenance. Sixteen Node scenarios and three SciPy references; see the [protocol](browser_wav_resampling_protocol_2026-11-26.md). |

These are development-branch commits. At this check, remote `main` remained at `7e72a12` and tag `v0.3.0` at `c406e00`. Pushing development commits does not merge them into main or replace an archived release. The current fixes must not be inferred to exist in the old release archive. GitHub Pages was also checked: its source is `main:/docs`, and its latest built revision is `7e72a12`. The live site therefore still serves the pre-correction implementation until a separately authorized merge/deployment.

## Consolidated regression evidence

The consolidated run at committed implementation **`8c6d3ae` passed 135 pytest cases with 1 deselected**, exit code 0, in **320.25 seconds**. The excluded test is `test_o3_holdout_hard_blocked`. Loader fixtures received only synthetic reserved-event metadata, and a hard block prevented actual reserved-manifest reads.

The [saved machine-readable record](validation/consolidated_verification_20261126T131531Z.json) contains the exact command/selection, pass/deselection counts, console evidence, source/test hashes and Python/NumPy/SciPy/Node versions. Its file hashes were checked against the committed implementation. Hashes describe raw working-tree bytes; Git line-ending conversion can change byte hashes in another checkout.

The 135 pytest cases include the wrapper that executes the 16 Node browser scenarios and the three SciPy FIR reference cases. These are included checks, not 19 additional pytest cases. Node tests use in-memory DOM/decoder fixtures. No failure or unexpected skip was reported. Implementation source stayed committed and unchanged during this step; only the documentation/history update was uncommitted.

This selection does not verify the real reserved-event filename mapping. It is regression verification using synthetic fixtures, not independent performance measurement on real O1/O2/O3 data. The earlier author-reported 87/87 result and independent 70-test subset are historical scope records, not the current test count.

## Historical results: correction versus rerun

| Evidence | Current interpretation / next step |
| --- | --- |
| Confirmatory `K=100`, 50,000 windows | Preserve the saved 316/50,000 result as historical. The arithmetic MC SE correction `0.00069 -> 0.0022231922256071334` needs no detector rerun. It does not establish all-seed current-code equivalence; any claim about current-code performance needs a separately recorded rerun. |
| Calibration 120/200/220 windows and 2,000-window negative control | Split-B sizes were 84/140/154. The 120-window FAP=1% variant was unattainable and must not be described as valid 1% conformal calibration. Recompute affected adaptive/exposure summaries into a fresh output location before comparing new performance. |
| Adaptive v0.3 phase 2 / phase 2b / corrected phase 2b | Recovery, available exposure and sensitivity metrics remain pre-audit. Repeat under an explicit protocol, preserving failures and unavailable decisions, before updating conclusions. |
| v0.2 coincidence / holdout | Historical FAR and significance used older exposure/event definitions. The current old benchmark callers still need valid active intervals for available significance. Define foreground/background rules and rerun separately without reusing held-out events as blind data. |
| GW170814 exploratory search | Previously used O2 epoch, historical zero coincident pairs. Current schema/deduplication changes do not by themselves reprocess the underlying strain. No novelty or astrophysical-origin conclusion is established. |
| Manual Python Workstation check at `986dd46` | Retain the recorded seed/configuration and observed SNR=18 and SNR=8 outcomes. No UI export was saved. A new manual run must get its own provenance bundle; do not invent an old run ID/hash. |
| Browser WAV measurements made with stride decimation | Historical preprocessing can create aliases. Rerun the source files under the documented FIR with saved original/decoded rates and padding disclosures; do not relabel old output as new-filter output. |

See the [audit erratum](AUDIT_ERRATUM_47a5475.md) and [current limitations](CURRENT_LIMITATIONS.md). No historical scientific rerun was performed by the implementation regression steps.

Before rerunning, check the output contract: historical experiment entry points such as `run_v03_phase2_adaptive_benchmark.py` and `run_workstation_calibration_study.py` write to fixed `results/` locations. They must first receive a reviewed fresh-output mechanism or equivalent isolated run setup. Do not launch them against the preserved historical directories. Record the protocol, revision/configuration, input provenance, available/excluded exposure and a unique destination before execution.

## Remaining verification tasks

- Native-browser decoding, rendering and latency with synthetic controls; save browser/version, input hashes and export bundle.
- Independent real-data false-alarm/sensitivity validation, including filter memory, edge effects and actual available exposure. Finite-sample rank attainability alone does not establish exchangeability.
- Browser seed-zero handling, complete synthetic-run provenance and diagnostic-cluster preservation differences.
- Historical exposure-aware reruns and unresolved distribution/independence assumptions; preserve old outputs and issue a dated superseding record.
- Attribution inventory, icon provenance and the owner's intended public scope for the reserved manifest. No O3 manifest contents or input data are inspected in this step.

Zenodo metadata/release operations are deferred. Do not move an existing tag or silently replace a published archive to make it match a later source tree.
