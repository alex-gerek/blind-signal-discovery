# Verification status and historical rerun queue

Updated on 2026-11-26. This index separates saved implementation fixes, regression evidence and scientific claims. Earlier journal entries keep their original status at the time of writing; dated follow-ups supersede those status snapshots.

## Saved implementation revisions

| Step | Saved commit on `development/v0.3` | Scope and evidence |
| --- | --- | --- |
| Documentation / erratum / version alignment | `99fce40` | Corrected current guide, preserved historical guide, fixed study paths and attainability interpretation, aligned `0.3.0` metadata. Documentation/AST/rank/link checks, not a scientific rerun. |
| Upload isolation and input checks | `8053741` | Unique paths, strict real-valued CSV/NPY input, rate/hash/time checks and error handling. 34 distinct selected pytest cases across final runs, including a six-scenario Node wrapper. |
| Abstention accounting and export precision | `8fab987` | Conditional rates/exposure, separate abstentions, compatibility semantics and shared precise CSV/JSON builders. 18 selected pytest cases across final runs and ten Node scenarios. |
| Browser WAV preprocessing | `8c6d3ae` | Fixed-parameter FIR to exact 512 Hz; 0–192 Hz useful resampled band, disclosed offline lookahead/padding, original/decoded rate provenance. Sixteen Node scenarios and three SciPy references; see the [protocol](browser_wav_resampling_protocol_2026-11-26.md). |

Publication follow-up (2026-11-26): the corrections and history commit `043f0eb` were merged through PR #8 into public `main` at `5e19ae4e6fb9ea177738fa4bae8310c11a2633e6`. GitHub Pages successfully built that revision. The `main` and development source trees match. The earlier pre-merge snapshot is preserved in research history. Tag `v0.3.0` remains at `c406e00`; the existing Zenodo archive does not contain these later fixes.

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


## Post-audit synthetic control pilot (2026-11-26, local results)

The [prespecified protocol](post_audit_workstation_control_protocol_2026-11-26.md) ran all 15 planned cases on unchanged detector source at `043f0eb` (same tree as main `5e19ae4`). The driver/protocol and local dirty-path list are saved with the [run manifest](../results/post_audit_workstation_controls/20261126T154934Z_6cbf4e69/manifest.json). No case was replaced or tuned. One initial system-Python import failure preceded the actual study and is separately recorded; the completed study used the existing project .venv.

- Stationary controls: three seeds, paired SNR 8/18, 220 background windows and 200 initial noise controls per seed. Count the shared noise once: **10/600 = 1.6667%**, 2,400 s available noise exposure, zero abstentions. This observed pilot fraction is above the nominal 1% parameter; three frozen calibrations do not establish a population violation or verify 1% real-data FAP control. Report it as observed, without a pooled common-rate Bernoulli confidence interval.
- Injections: SNR 8 **13/30**; SNR 18 **30/30**; zero abstentions. The changed 200-control target generation is not a rerun of the previous 80-control manual demonstration.
- Adaptive amplitude-step controls: dedicated STFT recovered and detected 10/10 strong SNR 30 chirps in each of six paired seed/FAP scenarios. The FAP settings share inputs, so these are not 60 independent injection trials. First new-regime available decision: index 144 at FAP=0.01; index 277 at FAP=0.005. First mature decision: index 277 at FAP=0.01, not reached at FAP=0.005 in this duration. Multiply the zero-based index by 4 s for target-relative start time, not universal recovery latency.
- At union FWER=0.005 the stricter channel requires 399 Split-B scores (570 initial windows at 30% normalization, 532 adaptive windows at 25%). The 320-reference/370-target pilot therefore does not verify full Union recovery; Union remained unavailable. This is expected insufficient data, not a demonstrated software defect.
- Capped companions (buffer=150, FAP=0.005): all three exposed resource blocking and ended unavailable; 108/180 noise windows were excluded, 72 remained available under the earlier calibration. Do not treat the entire companion as zero exposure or as verified zero false alarms.
- Independent verification recomputed counts/decisions/exposure from exported full-precision CSV, matched strict JSON, verified paired input hashes, and exactly regenerated all six stationary input arrays. Source hashes stayed unchanged. See [independent verification](../results/post_audit_workstation_controls/20261126T154934Z_6cbf4e69/independent_verification.json), [complete report](../results/post_audit_workstation_controls/20261126T154934Z_6cbf4e69/report.md), and [artifact hashes](../results/post_audit_workstation_controls/20261126T154934Z_6cbf4e69/artifact_sha256.json).

These are **15 experiment cases, not 15 new pytest passes**. The saved 135-passed/1-deselected regression record still describes the unchanged implementation. The pilot is limited to pink stationary noise and Gaussian amplitude steps; it does not supersede shifted_ar2 adaptive benchmarks, real-data FAR, GW170814, native UI/browser checks or the historical 100-calibration study. No O3 manifest/data, external input or detector tuning was used.

The [next-release checklist](NEXT_RELEASE_CHECKLIST.md) records the owner-confirmed author Lex Gerek, manual/operator follow-up, attribution inventory and release-version/DOI steps. No commit, push, tag, release or Zenodo change was made by this local control-study step.


## Historical adaptive synthetic replay (2026-11-26, local evidence)

Synthetic suitesA/B of run_v03_phase2b_corrected_benchmark.py were replayed without source tuning:12 cases,1,900 target windows, exact current input/reference regeneration, source/old-artifact hashes unchanged. A subsequent exploratory2-arm abrupt-case trace/control exactly reproduced the baseline and isolated threshold contamination. These are experiments, not additional pytest counts. See [findings](HISTORICAL_ADAPTIVE_REPLAY_FINDINGS_2026-11-26.md) for3/10 vs historical5/10, slower recovery, correct exposure, mature0/6 sensitivity and frozen broadband/fullUnion127/127 mature noise nominations.

The earlier blanket statement that all historical adaptive results remain unrerun is superseded ONLY for these syntheticA/B cases. Real-data suiteD, coherence suiteC, other phase2/phase2b versions and coincidence/GW170814 strain reruns remain pending. No physical discovery or validated adaptive admission repair is established.


## GW150914 current-source real-data control (2026-11-26)

At3497211, two genuine O1 H1/L1 files were processed with an independently fitted frozen PSD and separated calibration/target/control blocks. Both frozen/current Workstation targets recover the known event; background nominations4/320 H1,3/320 L1, off-source0/71 each. No injected waveform or adaptive regime update; this is a previously used positive control, not independent discovery or adaptive-recovery validation. One new-driver API failure is preserved separately. [Interpretation/evidence](GW150914_REAL_DATA_CONTROL_2026-11-26.md) includes exact output, mask/warning caveats, local-only derivative hashes and live UI ingestion limitation. Original135/1 regression source/test hashes remain the relevant unchanged implementation record; this experiment is not extra pytest passes.


## GW170814 bounded current-source replay (2026-11-26)

At8f78ca5, exact historical offline partition/vetoes onH1/L1/V1 reproducesL1-only recovery. All480target decisions checked;18 summary and22 candidate comparisons match old saved precision. No time-slide significance or adaptive run. [Known-event map](KNOWN_EVENT_CONTROL_MAP_2026-11-26.md) records both events, current subthreshold miss diagnostics and incomparable calibration-sample counts. This supersedes the unrerun label ONLY for this fixed detection/candidate scope, not historical FAR/significance.


## GW170814 saved-data miss interpretation (2026-11-27)

At30849e8,472 raw cluster normalizations and all global maxima were checked without a detector rerun. [Diagnosis](GW170814_MISS_DIAGNOSIS_2026-11-27.md) confirms that V1windowmax3.619 is a late374Hz cluster, whereas the near-reference maximum is.664. H1 fails after clustering; L1 passes allscales. This corrects interpretation, not primary decisions, physical cause or efficiency. [CodeQL record](CODEQL_TEST_REGEX_REVIEW_2026-11-27.md) documents approved dismissal of two test-only findings. Existing code/regression hashes remain unchanged.
