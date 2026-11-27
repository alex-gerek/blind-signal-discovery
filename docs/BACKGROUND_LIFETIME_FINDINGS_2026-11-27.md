# Background lifetime and reacquisition findings — 27 November 2026

[Protocol](background_lifetime_protocol_2026-11-27.md) · [Full report and figures](../results/background_lifetime/20261127T170612Z_189c239c/report.md) · [Independent verification](../results/background_lifetime/20261127T170612Z_189c239c/independent_verification.json) · [Paired contamination](../results/background_lifetime/20261127T170612Z_189c239c/paired_contamination_verification.json).

Sourcee21a7cc,unchangedcore. Three continuous synthetic pink traces, five paired cases and three experimental policies. DedicatedFAP0.005,284-window source (85A/199B), fixed2400soldest-sample age limit. This limit is UNVALIDATED, not a real-instrument recommendation. No production manager/Workstation implementation or scientific lifetime/FAR certification.

## Main measured outcomes

- Age alone misses an unannounced noise-law change. On the same line-case waveforms, external metadata changes policy eligibility but not core measurements. With no state notification, age+state behaves exactly like age-only, retaining many noise alarms before its scheduled expiry.
- Announced change: frozen-naive1299alarms/1518available noise windows; age-only717/1266; age+state3/552. These are different eligible exposures, not equalFAR sensitivity comparisons. Six wanted probes are abstained across the three state-qualified sequences.
- Reacquisition has a cost: state+age unavailable324windows/1296s per announced sequence; age-only84/336s. Gathering284four-second source windows alone takes1136s (18min56s). Under stable noise the age rule still incurs84unavailable windows/sequence.
- Replacement becomes mathematically available684, but line-case late wanted associations are0/15for both updating policies. Availability is not recovered sensitivity. Frozen naive has15/15late associations while also nominating everypost-change noise window; its apparent recovery is misleading.
- Four hidden replacementB intrusions reduce updated late wanted associations15/15->3/15 relative to the clean paired source. All later scores, normalizationA, wanted inputs and eligibility are identical. Replacement thresholds3.357758->99.897055,2.699896->35.741650,2.568230->84.123418. Source-OFF metadata does not prove purity.

Repeated signals outside the commanded source-acquisition periods do not trigger the metadata/age controller or contaminate those clean source blocks; this does not establish general contamination resistance. Three paired seeds, jointcolored windows and known-support association support descriptive case conclusions only. No automatic monitor for unknown changes was tested. No within-window change, real clock uncertainty, actual capped acquisition or live throughput validation; readiness/invalid metadata are only focused controller fixtures. Union recovery is not supported by199B.

## Checks and preservation

Independent saved-data reconstruction passes30calibrations,9480diagnostic measurements,23220policy records,33239candidate records and12000regenerated input-window records.135directuncachedchecks and23controller/ledgerfixtures pass.1734policy records preserve nonempty diagnostic candidate references while official p/threshold are null. These counts are scoped checks, not pytest passes or scientific proof. Causal use, exact expiry boundary, exposure, fit/activation timing and pairedwaveform identity checked.

No measurement failure, parameter sweep, excluded case or rerun. One post-measurement plot-preparation implementation error was saved and corrected without changing measurements. Final graphs visually checked. Source/dependency hashes unchanged; old results retained, no production/test changes, O3/data/manifest/network/real-waveform access or publication. L09/L10 remain open. New stage remains local.

Next decision: validate an observable distribution-change warning on independent sequences, treating a warning as uncertainty rather than permission to absorb a persistent signal into noise. Independently justify source purity and reacquisition conditions. Do not deploy the2400s setting or treat the acquisition contract as implemented.
