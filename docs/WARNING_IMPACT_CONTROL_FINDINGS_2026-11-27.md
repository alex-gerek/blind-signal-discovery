# Corrected cadence and warning-impact control — 27 November 2026

[Pre-run protocol](warning_impact_control_protocol_2026-11-27.md) · [Full report and figures](../results/warning_impact_control/20261127T175353Z_5cad334d/report.md) · [Numerical checks](../results/warning_impact_control/20261127T175353Z_5cad334d/independent_verification.json) · [Paired impact analysis](../results/warning_impact_control/20261127T175353Z_5cad334d/independent_impact_analysis.json).

Source `cc68b93`; five new seed traces, eight paired cases, unchanged default monitor and frozen core. This corrects the shared-signal protocol inconsistency from the [earlier pilot](REGIME_WARNING_MONITOR_FINDINGS_2026-11-27.md), which remains preserved. No automatic gating, model updates or production changes.

## Corrected comparison

Repeated grid has64injections per seed. Sparse uses13slots from that grid, with families set by global slot parity. Both schedules alternate and all65shared probes acrossfive seeds have exactly the same waveform, support and core output. Warning flags differ at1/65shared probes. A hypothetical warning gate would remove2shared recoveries in sparse context and3in repeated context. This is a clean monitor-history comparison, not a changed core sensitivity. Do not compare raw65/65versus315/320as a causal cadence effect because the extra probes are different examples.

## Warning versus measured detector impact

- Across5160paired targets (stable->2xstep and sparse->2xstepwithprobes), scores,p-values,nomination and support-association decisions are exactly unchanged. Monitor nevertheless warns on2200changed-level windows. The gain warning is not evidence of additional impairment for this normalized statistic.
- Noise-only stationary and2xstep cases each give25/2580search nominations. The spectral-line case gives2204/2580. Similar warning flags therefore accompany very different impacts. These are descriptive controlled cases, not a trained reliability classifier.
- Baseline25/2580~0.969%also exceeds nominalFAP0.005numerically. This study does not establish nominal0.5%control, diagnose the discrepancy or permit a pooledIIDtest on correlated windows/shared calibration draws. Gain invariance is not calibration validation.
- In eachpositiveprofile65/65known-support associations are observed. With the line present,40of65already exist in exact pre-injection noise; only25are new associations and10are positive-only window nominations. Stable and2xstep controls have0pre-existing noise associations and65new ones. Support overlap remains weaker than physical identification even when newly introduced by a synthetic signal.
- Retrospective fixed12-window future-noise counts after episode onsets:stable2/75,level0/60,line60/60. Ranges unioned and onset excluded. This small dependent-sample comparison is not predictive validation or an automatic impact gate.

Stablemonitor counts:21flagged windows/9episodes in4of5sequences; sparse18/8,repeated26/13. No monitorbudget inferred. Both truechanges startanewwarning at360inall5seeds. Persistent-signal/noise-line waveform ambiguity remains exact. A warning indicates a difference from the reference, not its physical cause or permission to absorb it into noise.

## Verification and scope

Independent reconstruction passes5reference models,5calibrations,20640monitor/search/sidecar records each,83536candidate records and32000regenerated input-window records.135manualNumPyFFTfeatures,120prefixcomparisons and121200integer-pathKScomparisons (maxp difference~1e-15). Corrected65shared identities,195positive/noise pairs,10gain comparisons and5future-exposure fixtures checked independently. Candidate/reference/source hashes unchanged. These are scoped experiment checks, not pytest passes or scientific validity proof.

All40cases retained, no measurement/preparation failures, interruption, retuning or exclusion. Figures visually checked. No O3/data/manifest/network/realwaveforms,production/test changes, publication or adoption. Five paired seed units,twochosennoisechanges,SNR18families and a fixed horizon limit inference. L09/L10 remain open. New stage is local.

Decision: retain generic warning annotations and actual search outputs separately. Do not convert every warning into calibration rejection or count every positive association as identification of the injection. An operational impact rule requires a newly frozen application-specific budget and independent validation, including persistent-signal ambiguity and source-purity conditions.
