# Background lifetime and reacquisition: causal timeline pilot — 27 November 2026

Pre-execution protocol, source e21a7cc48ce4122804f86e3397ee3894afa90639. Experimental controller only; unchanged detector and production managers/UI. This is a synthetic bench acquisition profile, not an astrophysical observing policy or a validated real lifetime.

## Fixed acquisition timeline

Three fresh seeds19101/19102/19103;800nonoverlapping4s windows at1024Hz per case. Generate ONE continuous pink trace of800*4096samples per seed, sigma1 alpha1 knee16 via generate_pink_colored, RNG SeedSequence([seed,1]). This is a jointly generated FFT-shaped process (periodic boundary model), not independent window noise. Partition it into windows without per-window RMS renormalization. No shared samples across acquisition windows; statistical correlation can remain. Reuse the base trace across five paired cases; independent replication unit is the seed, not pixels/windows/cases/policies.

Scheduled background acquisition intervals are windows[0,284) and[400,684), regardless of scores or policy. These correspond to a bench command to turn OFF the deliberate test source; it does not guarantee absence of unknown signals. Target decisions are recorded for every window284..799, including scheduled reacquisition windows. Acquisition/abstention costs remain in exposure; do not evaluate a new calibration on its own acquisition samples.

Initial calibration uses0..283, fitted through283 and first used284. Replacement uses400..683, fitted AFTER recording window683 and first used684. Norm_fraction0.30 yields85A/199B; dedicatedFAP0.005 just attainable. Record complete membership, full calibration/nulls and freeze hashes. No waveform from a future window can enter a calibration used now. Scheduled source membership must not be selected by candidate scores/truth labels.

## Five fixed cases and observable telemetry

1. stable: continuous pink, external stateA throughout.
2. announced_line: add3*sin(2*pi*110*t+phase) from window360 onward, phase drawn once from SeedSequence([seed,2]); time t is the global sample clock. External instrument-state telemetry changesA->B at360 before evaluating that window.
3. unannounced_line: EXACT same waveform as announced_line, but telemetry staysA. The controller must not receive the hidden change label or waveform-truth metadata. This tests the limits of metadata/age rules, not a drift monitor.
4. repeated_signals: stable noise, unchanged stateA; repeated transient probes outside commanded background periods. Scores do not trigger reacquisition or expiry.
5. calibration_intrusions: stable noise/stateA with the usual sparse probes plus four uncommanded transients at500/540/580/620 (all replacementB). The commanded-source-off flag is unchanged. The controller does not know these labels. This challenges the acquisition-purity assumption; do not remove the contaminated windows.

Sparse wanted probes in all cases except repeated_signals at[304,324,344,364,384,704,724,744,764,784], alternating chirp/burst. Repeated case indices list(range(288,400,8))+list(range(688,800,8)), alternating families. Intrusions alternate chirp/burst as above. Existing chirp defaults; burst4096Gaussian, Butterworth4zero-phase35..350Hz, support1.375..1.625s,Tukey0.15. Scale each to linear root-energy SNR18 relative to std(the noise window BEFORE injection). Clean waveforms use SeedSequence([seed,3,index,family_key]), family_key0chirp/1burst. Same sparse probes across cases; noise profile changes amplitude scaling as specified. Save actual support rectangles, input/base hashes and labels used ONLY for evaluation. Desired source is deliberately off in background intervals; intrusions are the explicit unknown-source exception.

## Three fixed policies; no tuning

- frozen_naive: use initial calibration forever, ignore state changes and age. This is the intentionally unsupported reuse comparator; mathematical computability does not imply scientific credibility.
- age_only: initial source may be used only if the complete target window end <= oldest source sample start+2400s. Equality passes. At window600 it expires; unavailable600..683. Use replacement from684; its oldest sample starts1600s, so it does not expire before this run ends. 2400s is a FIXED UNVALIDATED test parameter, not a recommendation derived from outcomes.
- age_and_state: same age/replacement schedule, but source state must equal external telemetry for the window. In announced_line, initial source becomes ineligible360..683. In unannounced_line it must behave identically to age_only: no hidden truth access. Replacement state comes from observed telemetry during its source interval; all source windows must have a single matching state.

Both updating policies use the replacement at684 even in stable cases; never keep whichever calibration gave better evaluation results. Replacement admission is based on prescribed acquisition and observed state only, not scores/purity truth. Intrusions remain admitted to demonstrate any failure. No automatic signal absorption claim can be resolved by a source-off command alone.

Keep raw diagnostic core measurements with their actual calibration IDs, scores, thresholds, empirical p and accepted candidates when a policy abstains. Official policy threshold/p are null, and unavailable decisions are not counted as noise negatives or available misses. Do not evaluate the not-yet-fitted replacement on earlier target windows, even diagnostically. Record used versus next state explicitly at the switch. No Union decision is computed; fixture verifies199B cannot support its current Sidak allocation atFWER0.005.

## Prespecified evaluation

Every policy/seed/case: actual interval-union available/unavailable noise exposure, alarms and conditional rates (null if zero exposure); wanted and intrusion detected/associated/missed/abstained counts; all-time operational coverage; exact expiry and fit/activation indices; false expiry/coverage cost on stable noise; stale-usage alarms before detection of an unannounced change; first later associated wanted probe after replacement, or explicit no observed recovery. Mathematical finite calibration at684 alone is not recovered sensitivity. Association is accepted candidate rectangle overlap, not physical waveform identification.

Keep three eras:before360,360..683,684..799, with actual denominators, and separate full-window nominations from associations. No inferential p-values, pooled-binomial CI, universal TTL or superiority/adoption claim from three seeds. Policies share inputs; the age and state interventions change availability. Zero alarms during abstention must not be called a performance gain. State-only cannot detect an unannounced change; report that limitation honestly.

Expected bookkeeping (not desired scientific outcomes): age-only has84unavailable target windows in every case. Age+state has324inannounced_line and84elsewhere. Late associated recovery can fail for line or contaminated source. Stable false expiry and reacquisition delays must remain visible. Preserve any result contradicting these assertions as an implementation failure, not an excluded seed. No parameter sweep/retry to improve sensitivity.

## Verification and preservation

Focused controller fixtures:exactageboundary/nextafter,unknown/invalidtime,statechange,empty/duplicate/overlapping acquisition intervals,resource/insufficientsamples,FAP0.005andUnionavailability,nullJSONandnonemptydiagnosticcandidatesduringabstention. These are controller/ledger fixtures, not production tests. Independently reconstructAstatistics,Bnull/rank,p-values,alltargetscores/candidateassociation,policyeligibility,nulls,actualexposure,source/targettiming and switch/recoverycounts from saved records. Verify announced/unannouncedwaveforms and all metadata-independent core results identical; only statepolicy gating differs. Regenerate all per-window hashes and compare uncached raw output for first3source/targetwindows percase plus latefirst3windows.

UniqueUTC/UUID outputs; exact pre-run protocol, drivers/helper copies, source/dependency/dirty/environment records, input features/provenance, calibrations and freeze hashes, full timelines/diagnostics/errors, verifier and report. No old result overwrite, target selection, production/test edits, O3/data/manifest/network/real waveform access, merge/release/Zenodo, commit/push. Initialstatus:design fixed before measurement. Real acquisition purity, repeated-testing drift monitor and validated age limits remain separate future work.
