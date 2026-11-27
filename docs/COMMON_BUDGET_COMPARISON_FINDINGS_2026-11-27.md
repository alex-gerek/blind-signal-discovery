# Independent common-budget comparison — 27 November 2026

[Full report and figures](../results/common_budget_comparison/20261127T114754Z_eefa4de6/report.md) · [Prespecified protocol](common_budget_comparison_protocol_2026-11-27.md) · [Independent verification](../results/common_budget_comparison/20261127T114754Z_eefa4de6/independent_verification.json).

Frozen baseline versus fixed approximately1s temporal median; all new synthetic data, no GW170814/other waveform input.200normalization records,1800threshold records,3600held-out noise records and600distinct signal units at two paired amplitudes. Independent streams15001..15004, no cross-split noise-hash intersections. Thresholds frozen before held-out generation.

Budget .01false-alarm windows/4s observation (9window-equivalent/hour) applies to the specified IID one-third stationary/smooth/within-window-step mixture. Conservative rankm10,k1791,operatingalpha10/1801,per-arm97.5% tolerance allowance. This is a common upper cap, not exact equality of actual FAR.

Held-out noise: original21/3600=.5833%,exact95%CI[.3614%,.8903%]; prototype18/3600=.5000%,CI[.2966%,.7891%]. Both upper bounds satisfy the1%overall mixture cap. Conditional drift/step upper bounds exceed1%; no per-regime or real continuous-stream assurance. Rates5.25/4.5window-equivalent alarms/hour are observations, not universal operating rates or an astrophysical event FAR.

ONE primary endpoint, weak broadband burstSNR8:33/300original versus34/300prototype,6slow-only/5baseline-only,exactMcNemarp1.0. The prespecified improvement criterion is not met; no equivalence claim. Secondary descriptive groups: weak chirps187->160/300(3gained,30lost),strong chirps300/300both,strong bursts264->289/300(25gained,0lost). Do not promote the favorable strong-burst group into a new primary result.

Independent checks pass:2calibrations,9600arm/window decisions,5296accepted candidates,6800regenerated input records,600paired signal units,126manual map-arm checks; rank/interval/test identities checked separately. No measurement error, interruption, dropped case or post-evaluation retuning. Counts are not pytest passes or independent real-data validation.

Decision: keep original default; archive candidate as a specialized trade-off. One calibration realization, chosen mixture and limited families do not establish general superiority or sensor performance. Any specialized follow-up needs a new fixed protocol and independent evaluation, not extending this sample until significance. No src/app/main/release/Zenodo/O3 changes. Article/CodeQL drafts,DNA-window and L09/L10 remain separate.
