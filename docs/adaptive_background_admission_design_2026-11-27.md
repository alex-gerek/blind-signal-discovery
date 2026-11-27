# Adaptive background admission: requirements and validation gates — 27 November 2026

Status: design stage only. Source snapshot `4a1f590a424e81a125d593611a7f162cab54af42`, development/v0.3. Owner authorized point 2: define requirements and a validation plan before changing the detector. No new measurements, candidate implementation, production changes, or automatic adoption. This is NOT a frozen executable experiment protocol: exact generators, seeds, sample sizes and a candidate rule must be fixed in a successor protocol before any outcome-driven evaluation.

## 1. Evidence and the problem to solve

[Historical replay and paired diagnostic](HISTORICAL_ADAPTIVE_REPLAY_FINDINGS_2026-11-26.md) established L09 in one constructed sequence: four earlier chirps enter threshold Split B of the 265-window mature calibration (66 normalization / 199 threshold). Threshold 163.83940094146237 is first used at window 367; all six later available injections are missed. Replacing only those four earlier windows with exact original noise leaves later scores and switch timing unchanged, gives threshold 1.6580395349128854 and restores 5/6 detections. This is an oracle diagnostic, not an implementable repair or prevalence estimate.

Current code admits accumulated observations to STFT recalibration in `src/bsd/adaptive_regime_v2b.py:877-885` and `:964-973`. The glitch-free check at `:905-909` and `:992-997` gates refresh of the glitch reference after the STFT calibration has already been fitted. It does not certify that the STFT calibration input is signal-free. Internal continuum stability is also not proof of signal-free background.

L10 is distinct: frozen broadband/full Union nominates 127/127 mature noise windows after the same shift. A repaired STFT admission policy must not be advertised as fixing frozen broadband or Union FWER. See [current limitations](CURRENT_LIMITATIONS.md).

## 2. Required distinctions

Maintain separate fields for input validity, admission decision, calibration computability, calibration credibility, detection availability, and measured sensitivity. An available finite threshold is not evidence of restored sensitivity; low noise-alarm counts alone are not success. A block may be internally stable but have unknown background purity. Proposed labels below are design requirements, not existing API fields.

Every observation needs a durable identity, acquisition interval and input/provenance hash. Admission records must retain the observed measurements, accept/quarantine/reject decision, explicit reason, policy version and decision time. Retain all diagnostic candidates, including during abstention; do not silently delete rejected calibration records or candidate evidence.

Only past observations may affect causal decisions. A calibration fitted through window i first affects i+1; keep the existing used-state snapshot semantics at `adaptive_regime_v2b.py:709-717`. Export used calibration ID, exact used threshold or null, empirical p-value or null with reason, channel availability and next-state update separately. Strict JSON must contain neither NaN nor Infinity. Do not report unavailable significance as p=1 or an unavailable alarm rate as zero.

## 3. Admission options and constraints

Preferred first diagnostic route: explicit independently sourced background blocks, with valid provenance and separation checks, compared with the unchanged accumulated-stream baseline. This establishes a background-source control, not a deployable autonomous adaptive solution. No such external stream may be assumed available in real deployment; synthetic truth/paired noise is inaccessible to the operational policy.

For an automatic candidate, investigate block-level quarantine and delayed reassessment, without discarding individual observations merely because their detector score is high. A new regime may contain legitimate large noise excursions. Reusing the old threshold to delete them can truncate the calibration tail, lower thresholds and increase false alarms. Random splitting after score-based selection does not remove this selection bias. Temporal separation alone also does not establish distributional exchangeability.

No trimming fraction, score cutoff, repeated veto, gap rule, temporal-kernel change or reference update is approved by this document. Before a candidate can run, specify its exact observable-only predicates, latency, block boundaries, release/discard policy, resource cap and statistical target. Freeze these before evaluation. The original and slower temporal-median methods are not additional arms here; detector-core mathematics remains fixed.

## 4. Required falsification cases

A successor protocol must include all cases, with fixed seeds and no case removal:

- Exact historical replay and paired oracle diagnostic as developmental regression evidence; never count their reused data as independent confirmation.
- Fresh no-signal stationary backgrounds with their full natural upper tail retained. Log which noise extremes each policy quarantines or rejects. Include ties, finite extremes, and invalid/non-finite inputs as separate software fixtures.
- Independent streams with a sustained level step, smooth drift, spectral resonance and impulsive glitch storms. Distinguish a changed noise regime from a persistent signal; ambiguous cases may remain unvalidated rather than being called clean.
- Sparse and repeated transient chirps and broadband bursts during accumulation and after recovery, at weak and strong amplitudes. Vary placement relative to calibration blocks and switch boundaries by a predeclared schedule; preserve paired noise bases and lost baseline recoveries.
- Signal-free and contaminated calibration at a fixed schedule, plus a separate end-to-end adaptive comparison. This separates calibration-source effects from changed transition timing, abstention or reference normalization.
- Insufficient admitted data, empty/repeated IDs, overlapping intervals, a stalled quarantine, and resource-cap exhaustion. Test supported FAP settings including 0.005 and 0.01 and channel-specific Union sample requirements from actual split sizes, not a fixed nominal buffer count.
- Legacy manager/export/UI contract fixtures where shared admission semantics would affect them; never claim both managers are fixed after testing only V2b.

## 5. Separation and accounting

Development, calibration, policy validation and final evaluation streams must be distinct. Do not use injection truth, event time/frequency gates, future data, or paired pre-injection noise to make an admission decision. Truth is for scoring and an explicitly labelled oracle comparator only. All arms receive identical recorded streams; preserve independent sequence IDs and amplitude pairing.

Record membership and order of Split A and B for every proposed and activated calibration, normalization parameters, null scores, rank, ties, minimum resolvable p, frozen references, seeds, resource use and exact effective switch window. Admission decisions and resulting null scores must be auditable independently.

Exposure is a union of actual valid acquisition intervals, separately for each channel and availability state. Duplicate IDs or overlaps must not double count time; gaps and unavailable decisions contribute no available exposure. Keep all-time operational coverage alongside conditional rates. Report abstentions separately from available misses, including time to first finite calibration and time to demonstrated sensitivity on later signals. Permanent abstention cannot win by producing zero alarms.

For both baseline and candidate report noise alarms/counts and time, available/all noise exposure, eligible/detected/missed/abstained signals, support-associated recovery, whole-window nominations, quarantined counts, threshold evolution and latency. Keep STFT, broadband and Union separate. Candidate association is not physical identification. No all-window denominator may conceal unavailable decisions.

## 6. Statistical and adoption gates

Specify the target population and common alarm budget before execution. A calibration target FAP is not automatically an achieved conditional or streamwise FAR. The previous common-budget IID mixture experiment does not establish adaptive validity. Its primary improvement criterion failed; original core remains default.

For adaptive sequences, serial dependence, data-dependent admission and threshold updates invalidate a naive independent-window binomial or McNemar analysis. Use independent generated sequences as the replication unit for the operational comparison; fix sequence counts, endpoint definitions, dependence-aware uncertainty calculation and multiplicity before running. Separate a fixed-calibration IID control where exact binomial inference can legitimately be checked. Preserve an inconclusive outcome if the chosen sample cannot resolve the target.

Predeclare one primary trade-off criterion: measurable post-shift signal preservation subject to a validated noise-alarm budget, with an explicit availability/latency requirement. Exact margins, confidence level and sequence count remain UNSET; they require a feasibility/sample-size calculation independent of evaluation outcomes. No significance testing or successful-repair claim is permitted until these are frozen. Secondary amplitude/family/regime groups remain descriptive unless separately adjusted.

Reject adoption if the policy wins only through abstention, suppresses noise extremes without valid out-of-selection calibration, relies on truth/future access, loses required weak-signal sensitivity, exceeds the budget, misreports used calibration, or conceals L10. Even a positive result supports only its stated synthetic generator and policy scope, not arbitrary sensors or astrophysical significance.

## 7. Next executable step and preservation

Next: build measurement-only instrumentation in experiments for exact admission membership, used/next calibration states, exposure and paired outcome replay. First verify unchanged-baseline reproduction, then propose ONE observable-only rule and freeze a successor protocol with exact distributions, counts/seeds, margins and analysis. Do not optimize on the final evaluation data. Production changes need separate review after results.

Every future run receives a new UTC/UUID output directory and a copied pre-run protocol, code/source hashes, dirty-path record, all inputs/provenance, calibrations, decisions, candidates, exclusions, errors and verification artifacts. Never overwrite old results. No O3 or reserved manifest access, new real-data reads, main merge, release, Zenodo change, commit or push in this design step.

Outcome of this step: existing evidence and code paths inspected; requirements documented. No experiment run and no claim that L09/L10 is resolved. README/article-link changes and CodeQL/article drafts remain unrelated.
