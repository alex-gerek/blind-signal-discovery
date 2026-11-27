# Separate-background source control — 27 November 2026

[Protocol](separate_background_control_protocol_2026-11-27.md) · [Full report and figures](../results/separate_background_control/20261127T154355Z_4f134945/report.md) · [Independent verification](../results/separate_background_control/20261127T154355Z_4f134945/independent_verification.json).

Source d049e2b,unchangedcore. Three fresh seeds/fournoiseprofiles/fourcalibrationsources/twoFAPs. Separate MATCHED background generated independently under the exact same noise law as the target, with no injected signals; it retains noise extremes. Signal-free matching is known by synthetic construction, not inferred by an automatic admission rule. Paired pre-injection noise is a diagnostic oracle; independent old-regime noise tests mismatch.

## Observed fixed-source results at FAP0.01

Arms receive the same targets; all decisions available. Noise denominators600perprofile across3calibrationdraws, eachsignalgroup30units. These are descriptive counts, not matched achievedFAR or a pooledbinomialguarantee.

| Noise profile | Noise alarms: contaminated -> separate matched | Chirp SNR18 associations | Burst SNR18 associations |
| --- | --- | --- | --- |
| Stationary pink |0/600 -> 4/600|17/30 -> 30/30|0/30 -> 27/30|
| Post-step level2x |0/600 -> 4/600|23/30 -> 30/30|1/30 -> 29/30|
| Pink +110Hz line |3/600 -> 11/600|2/30 -> 4/30|0/30 -> 0/30|
| Pink +rare impulses |1/600 -> 4/600|23/30 -> 24/30|1/30 -> 1/30|

Association is accepted-candidate rectangle overlap with known support, not physical identification. Full reports retain everyweakSNR8group and FAP0.005; weakbursts remain poorly recovered (stationary0/30evenwithmatchedsource). Do not select only the favorable stronggroups or claim broad sensitivity improvement at equal trueFAR.

The old-regime background in the spectral-line profile yields600/600noisealarms and120/120positive associations atFAP0.01. This is a warning against judging success by positive-window counts alone; it does not establish that every individual positive candidate is caused by noise. Separate source purity does not ensure noise matching. Matched-source11/600inthisprofile also does not establish1%operationalcontrol fromthreecalibrations.

The paired-noise oracle changes exactlyfourBwindows and keeps allA/targetscores identical, isolating threshold contamination in this constructed control. Separate matched replacesAandB, so calibration-sampling variation affects both normalization andthreshold. It neednot reproduceoracleperformance. Signal/background sources and target noise bases have no unintended hash overlap; this is preservation evidence, not a proof of PRNG independence.

## Verification and decision

Independent saved-data checks pass96calibrations,23040arm/alpha/target records,22090candidate records,15600regenerated calibration/target inputs and180directuncachedrawchecks. Frozen calibration hashes unchanged; reconstructAstatistics,Bscores/ranks,p-values,strictdecisions,candidate standardization/bounds/supportassociation,paired gains/losses andactualavailableexposure. No abstention in this experiment. No measurement failure, droppedcase, retuning orproductionchanges; noO3/data/manifest/network/realwaveformaccess. These are scoped researchchecks,notpytestpasses orindependent scientific replication.

Conclusion: suitable separately acquired background can remove a major contamination effect in these synthetic cases without abstaining, but does not solve all sensitivity/noise-mismatch problems. No autonomous acquisition or quarantine-release policy was implemented. L09/L10 remainopen; fixedsource tests do not establish adaptive timing, realbackgroundpurity,Unioncontrol or astrophysical significance.

Next: specify a realistic acquisition/provenance contract (samechannel/state/preprocessing, temporal separation, contamination/staleness risks and what information is actuallyavailable) before an independent sequence-level operational test. Source separation alone is not such a contract. New results remain local pendingreview; no commit/push/release/Zenodo.
