# Background acquisition and expiry contract — 27 November 2026

**Draft research contract, not implemented in Workstation.** Evidence revision: `e21a7cc48ce4122804f86e3397ee3894afa90639`. This document defines a practical acquisition path, information requirements and reasons to suspend a background's applicability. It does not certify any real background or introduce a new detector, automatic age limit or admission rule.

See the [separate-source control](SEPARATE_BACKGROUND_CONTROL_FINDINGS_2026-11-27.md), [rejected tail guard](BLOCK_TAIL_GUARD_FINDINGS_2026-11-27.md), [initial admission design](adaptive_background_admission_design_2026-11-27.md) and [research history](RESEARCH_HISTORY.md). A [manifest template](background_source_manifest.example.json) makes missing information explicit; it is not a valid acquisition record or an input supported by the current UI.

## 1. Define what the background represents

Before recording, identify the intended signal, what counts as noise, the channel and operating mode. A persistent line may be an unwanted background feature for a transient search and the desired signal for a persistent search. Do not silently relabel it after looking at the detector output. Keep STFT, broadband and Union claims separate.

The first practical route is a **separately recorded interval from the same physical channel**, with the same instrument state and compatible preprocessing. It is not a different detector/channel used as a surrogate simply because the units match. A cross-channel background needs its own transfer model and validation.

Choose the acquisition schedule and external quality masks before inspecting target scores. Examples for later work are a same-channel O1/O2 off-source recording with documented quality/injection flags, or an audio sensor recording under the intended acquisition setup during a controlled period without the deliberately generated test source. Neither proves that all unknown signals are absent. No new physical files, O3 or reserved manifest were opened for this design.

The reference must represent the intended noise population, including its legitimate impulses and extremes. The failed tail guard demonstrates why selecting only quiet windows is not a generic repair. A room-tone reference does not automatically represent every background in a music mixture.

## 2. Acquire, preserve and qualify the source

1. Record original acquisition intervals, clock origin and uncertainty, physical instrument/channel identity, gain/response version, units, sample rate, bandwidth and acquisition-state history. File modification dates and a shared label such as CSV_TIME_COLUMN do not establish a common absolute clock.
2. Retain original bytes and hashes plus data-content identity where files may be repackaged. Log access/license provenance and external quality evidence. A new filename or hash does not make copied samples independent.
3. Apply predeclared masks derived from instrument logs or acquisition controls, retaining excluded intervals and reasons. Do not use this detector's high scores, target event gates, synthetic labels or a retrospectively selected quiet subset to claim a pure operational background. Known injection masks must include a justified guard for signal duration and processing support.
4. Preserve chronological membership of normalization A, threshold B and a separate quality/validation segment. Do not share target samples with these fitted components. Independence and model assumptions must be assessed separately from identity/overlap checks.
5. Freeze the preprocessing/model configuration and source fingerprint: resampling, antialias filter, whitening/PSD fit intervals, padding, effective frequency bands, window grid and code revision. Shared frozen nuisance models must record their independent training data and conditional assumptions; they must not have been fitted on the target.
6. Record actual accepted exposure after masks and gaps. Never replace missing/invalid rows silently or fill a calibration shortage with duplicated windows. Count complete, unique eligible observations, not file length alone.

Temporal separation includes raw processing dependencies, not only the displayed four-second windows. Padding or a fitted PSD may reach into an otherwise disjoint target interval. The required embargo follows from the actual filter/model support and a validated boundary-error tolerance; there is no generic one-second guarantee. Unknown dependencies mean separation is unverified. A causal acquisition must precede the target; an after-target reference can only support an explicitly retrospective study.

## 3. Enough data means more than an attainable threshold

Use the actual split and empirical-p-value rule. For N_B threshold samples, the minimum resolvable p is 1/(N_B+1). Compute the required total from the real rounding and normalization fraction, including separate Union requirements.

| Purpose | Normalization fraction | Exact minimum total | A | B |
| --- | ---: | ---: | ---: | ---: |
| Dedicated FAP0.005, research split |0.25|265|66|199|
| Dedicated FAP0.005, Workstation split |0.30|284|85|199|
| Union FWER0.005, research split, current Sidak family allocation |0.25|532|133|399|
| Union FWER0.005, Workstation split, same allocation |0.30|570|171|399|

These are arithmetic attainability minima, NOT recommended study sizes, effective-independent counts or a proof of operational false-alarm control. They were checked by enumerating totals with the actual round-based split. Union allocation assumptions remain a separate scientific question.

The source-control pilot's matched line-profile reference still produced11/600noise nominations at nominalFAP0.01. Three random calibrations cannot establish a fixed deployment guarantee. Predeclare an independent validation set, desired uncertainty and intended operating horizon before evaluating background suitability. Correlated sequences and shared calibrations cannot be treated as thousands of independent trials.

## 4. When the background stops being supported

Keep separate dimensions: data integrity, temporal separation, state compatibility, source-quality evidence, age coverage, statistical suspicion and calibration computability. A finite threshold does not establish any of the other dimensions.

| Observation | Proposed classification and action |
| --- | --- |
| Invalid samples, duplicated observations, unresolvable processing overlap or inconsistent clocks | Invalid/unverified acquisition; do not activate a claimed calibrated source. Keep diagnostic measurements and reasons. |
| Documented channel, gain, hardware mode, filter, whitening-model or acquisition-response change | Previous source is outside its declared profile from the first affected interval, unless a separately validated invariance/transfer rule covers that change. Do not assume normalization makes every change harmless. |
| Age exceeds a limit validated for this acquisition profile | AGE_EXPIRED: background is outside the tested age range, not proven physically changed. Acquire a new qualified block; do not refresh its age by saving or renaming it. |
| Age limit, clock origin or state coverage is unknown | AGE_OR_STATE_UNVERIFIED. Do not label the background fresh or treat unknown age as zero. Exploratory measurements can remain explicitly unvalidated. |
| Predeclared drift monitor raises suspicion | SUSPECT_DISTRIBUTION. Preserve the anomaly and investigate; it may be signal rather than a noise change. Do not automatically teach it into the background. |
| Monitor does not reject its null | No evidence from that particular monitor. This is not proof of stationarity, source purity or calibrated FAP. |
| Not enough qualified observations for the requested channel/FAP | Explicit channel-specific unavailability; do not substitute a finite threshold, zero FAR or a zero miss count. |

**No universal lifetime is specified.** A background recorded a minute ago can be unsuitable after a state change; an older recording may still be comparable, but that requires evidence. Time alone neither guarantees nor disproves relevance.

Track the background's start/end, its collection span, the gap to the target and the oldest contributing calibration sample. Track fitted noise-model age separately. Looking only at time since the last background sample hides an old A partition or PSD model. A declared validity interval must cover the complete target window and its processing dependencies, including clock uncertainty. Boundary equality must be specified before implementation; half-open acquisition intervals touching without shared dependencies need not overlap.

Numeric age limits must be derived for a specific instrument/profile from independent continuous-sequence tests with realistic state changes and contamination. A chosen maximum age is a policy parameter until validated. The manifest therefore leaves all age limits null and policy_validated false; it cannot assert fresh calibrated significance.

## 5. Monitor without selecting a convenient noise tail

Candidate diagnostics include robust amplitude, bandwise spectral continuum, line occupancy, impulse/crest statistics and dependence. V2b already has continuum/impulse diagnostics and stability checks (`adaptive_regime_v2b.py:587-609`), but their existing thresholds are not automatically validated for a new acquisition contract.

Keep the monitoring reference and its validation role explicit. Repeated KS tests with a nominal per-test p cut do not by themselves control false expiry alarms over a long deployment. Predeclare the monitoring horizon, persistence/latency, multiplicity handling and sequence-level performance targets. Test clean stable stretches, real distribution changes, rare legitimate noise extremes and repeated signals. Do not pick a monitor because it marks the known historical injections.

The rejected whole-block guard is not revived here: neither deleting individual high scores nor choosing whole blocks based on the same tail and then reusing ordinary conformal claims is automatically valid. Reacquisition must follow an acquisition policy, not repeated attempts until the background looks quiet.

## 6. Activation, expiry and recovery records

Each background/model gets immutable IDs and versioned source evidence. Freeze exact A/B membership, fitted state, quality masks, acquisition/processing footprints, policy version and any validated limits before activation.

A fit using data through window i first affects i+1, preserving the existing V2b used-state snapshot (`adaptive_regime_v2b.py:709-717`). A known instrument-state change has its own effective time; it must not wait for a score-based monitor. For a change discovered during or after a window, preserve the threshold and nomination actually used and append an explicit affected-interval/credibility annotation. Do not rewrite history to pretend the decision used a different calibration.

On expiry/suspicion, raw scores and diagnostic candidates remain available with their actual provenance. Official calibrated significance may be unsupported; represent that explicitly, not as p=1 or FAR=0. Do not silently repurpose the current detection_available field: a future implementation needs distinct source-eligibility and statistical-credibility fields alongside mathematical FAP availability.

Recovery requires a newly acquired qualified block, sufficient channel-specific samples, a recorded causal switch and later sensitivity/noise checks. Marking a calibration mature, finite or available is not recovery evidence. Keep unavailable exposure, operational coverage and all signal abstentions beside conditional rates. Frozen broadband/Union references need their own source/expiry policy; an STFT update is not a full Union repair.

## 7. What exists today and what remains to be tested

Current ingestion checks (`io_loader.py:429-480`) compare finite effective rates, file hashes and comparable recorded intervals, or emit an unverified-independence warning. Workstation calls this check before slicing (`workstation_runner.py:459-470`). Those checks do not establish channel/state/units/preprocessing equivalence, background purity or a validated age policy. The timestamp-disjoint status is narrower than scientific independence.

This contract and template add no enforcement, UI controls or detector changes. They must not be cited as a completed automatic freshness fix.

Before implementation, freeze one concrete acquisition profile and an executable validation protocol covering:

- same channel/state with old and recent stable references;
- known gain, processing, spectral and impulse-distribution changes, including changes within a window;
- old/unknown clock epochs, clock uncertainty, padded dependency overlap, duplicates and missing metadata;
- repeated target signals without known acquisition-state change, so a monitor does not simply erase them as background;
- exact expiry boundaries, a newly acquired replacement, delayed release and insufficient-data/resource limits;
- FAP0.005 with the Workstation and Union splits, strict JSON nulls and nonempty diagnostic candidates during unavailability;
- independent sequence-level noise alarms, missed/abstained signals, availability and recovery latency, with predefined margins and confidence scope.

Initial contract outcome: practical source route and expiry reasons are specified. Numeric age/monitor settings and real source purity remain unvalidated. No new waveform experiment, real-data access, O3, production change, release or Zenodo action occurred in this step.
