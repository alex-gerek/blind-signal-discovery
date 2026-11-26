# Known real-event control map (2026-11-26)

This is a small map of previously used controls, not a representative efficiency estimate. Source versions and calibration recipes differ as recorded; no threshold tuning.

| Known event | Site | Event score / threshold | Outcome before veto | Current evidence |
| --- | --- | --- | --- | --- |
| GW150914 |H1 |35.7673 /2.80539 |Detected |[O1 control](GW150914_REAL_DATA_CONTROL_2026-11-26.md) |
| GW150914 |L1 |8.52803 /3.80464 |Detected |[O1 control](GW150914_REAL_DATA_CONTROL_2026-11-26.md) |
| GW170814 |H1 |1.14481 /2.59054 |Missed |[O2 report](../results/gw170814_current_control/20261126T174011Z_e21955d7/report.md) |
| GW170814 |L1 |8.89240 /2.24008 |Detected |[O2 report](../results/gw170814_current_control/20261126T174011Z_e21955d7/report.md) |
| GW170814 |V1 |3.61917 /4.00689 |Missed |[O2 report](../results/gw170814_current_control/20261126T174011Z_e21955d7/report.md) |

## GW170814 bounded replay conclusions

The prior L1-only outcome reproduces on current source8f78ca5. H1/V1 form seed-qualified raw time-frequency clusters, but ALL per-scale standardized maxima stay below the global threshold; these are pre-veto misses, not unavailable decisions or veto losses. Presence of a raw cluster near the reference does not establish that it represents the physical wave.

- H1 empirical p54/496=0.10887097; L1 p1/496=0.00201613; V1 p5/496=0.01008065.
- V1 is just above requestedFAP0.01. Finite rank resolution matters; this is correct available non-detection, not the earlier upward-tolerance defect. Detection decision was available; no threshold adjustment follows.
- Calibrated off-source counts (151windows/site): H1 raw2/post-veto1; L1 raw3/post-veto3; V1 raw2/post-veto2. These are candidate nomination counts, not proven astrophysical significance or pure-noise FAR.
- Coarse foreground peak matches are zero on/off-source, before/after veto. No time-slide or waveform-delay significance was run.
- All18 historical summary checks and all22 candidate coordinate/score/p/veto comparisons match saved historical precision. New full-precision data are retained. Historical original-array bytes are unavailable, so recipe/output agreement is not proof of bitwise historical input identity.

## What the misses do and do not explain

Confirmed location in the pipeline: after raw clustering, at the standardized global threshold. Different instrument response, local noise, STFT resolution and preprocessing may explain different scores physically; none was isolated by a causal experiment in this bounded replay. No post-hoc lowering of thresholds or event selection is justified.

GW170814 used707calibration windows from historicalCalA+CalB (including post-targetCalB), while GW150914 used320pre-target windows. Both are fixed current-code controls, but they must not be pooled as a universal3/5 sensitivity claim. The O2 replay is offline, not an adaptive-recovery validation. L09/L10 remain unresolved.

[Protocol](gw170814_current_control_protocol_2026-11-26.md), [independent checks and per-scale diagnostics](../results/gw170814_current_control/20261126T174011Z_e21955d7/independent_verification.json), [all22 historical candidate comparisons](../results/gw170814_current_control/20261126T174011Z_e21955d7/historical_candidate_comparison.json), [three-site figure](../results/gw170814_current_control/20261126T174011Z_e21955d7/event_comparison.png), [above/below-threshold raw-cluster display](../results/gw170814_current_control/20261126T174011Z_e21955d7/event_raw_cluster_diagnostics.png).

The detector is demonstrably executable and useful for these measured-data control analyses, with explicit successes/misses. This does not establish a new physical discovery, broad real-data validation or completed research. A project article may describe a working research prototype if it retains this scope and failures. Article drafting is deferred as requested. Data attribution:GWOSC/LIGO/Virgo, CC BY4.0; original/local derived waveforms outsidegit. No O3, tuning, main merge, release or Zenodo update.
