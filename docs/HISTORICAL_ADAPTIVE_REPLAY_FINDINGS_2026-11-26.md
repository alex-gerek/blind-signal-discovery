# Historical adaptive replay: findings and limits (2026-11-26)

Current source043f0eb; historical recipeedd1ec9. All changes here are measurement drivers, new evidence and documentation. Detector/generator/calibration source was not edited.

## Evidence

- [Prespecified replay protocol](historical_adaptive_replay_protocol_2026-11-26.md).
- [Twelve-case replay report](../results/historical_adaptive_replay/20261126T163629Z_75c7a756/report.md) and [independent accounting](../results/historical_adaptive_replay/20261126T163629Z_75c7a756/independent_analysis.json).
- [Exploratory paired follow-up protocol](historical_adaptive_followup_protocol_2026-11-26.md), [report](../results/historical_adaptive_followup/20261126T164211Z_dd2ea691/report.md) and [independent verification](../results/historical_adaptive_followup/20261126T164211Z_dd2ea691/independent_verification.json).

## Confirmed priority finding: adaptive signal contamination can destroy sensitivity

In the historical abrupt-step recipe, earlier known injections330/340/350/360 enter the THRESHOLD_B partition of the265-window mature recalibration (66 normalization +199 threshold windows). The threshold becomes163.83940094146237 and is first used at window367. All six later injected windows370/380/390/400/410/420 are available/mature but missed.

A post-hoc paired intervention replaced ONLY those four earlier injected windows with their exact original pre-injection noise, keeping all other samples (including the six later injections) bitwise identical. Mature timing stays367 and all six later STFT scores are exactly identical. Threshold becomes1.6580395349128854 and5/6 later injections are detected. Both arms observe0/127 false alarms in the same later mature noise subset; that small sample does not establish general false-alarm control.

This confirms the mechanism in this constructed case: calibration admitted earlier signals as noise-tail examples and elevated the threshold. The old explanation that all post-recovery misses are due to spectral resonance is insufficient for the controlled difference. It does not rule out resonance-related masking for other misses.

- Input admission: src/bsd/adaptive_regime_v2b.py, progressive calibration uses m_obs from accumulated windows (around lines877–900).
- Exact membership and thresholds: each follow-up arm's calibration_trace.json.
- Exact later scores, available decisions and parent-baseline equality: independent_verification.json and windows.json.
- Scope: confirmed design limitation/failure to retain sensitivity, not evidence of a floating-point error or a universal failure. This single exploratory intervention does not establish prevalence or methodological novelty.
- Do not implement oracle signal deletion as a fix. Excluding noise extremes can bias calibration. Review background admission/segmentation and false-alarm validity under an explicit protocol before changing detection code.

## Confirmed priority limitation: frozen broadband/Union does not adapt to changed noise

After the abrupt background step, dedicated broadband nominates all387 post-step noise windows; in the mature clean subset it and full Union nominate127/127 noise windows. Dedicated STFT nominates0/127 in that same mature subset, while also missing all6 later injected signals because of its contaminated threshold.

The broadband reference is deliberately frozen: a sustained change in the noise spectrum is also an anomaly relative to the old reference. Availability/FULL_UNION_AVAILABLE does not establish1% FWER on the changed noise process. These flags are not discoveries of physical signals. This quantifies a known drift ambiguity rather than proving a new physical phenomenon.

The persistent drift matrix retains the historical aggregate results:69/100,75/100,80/110 and73/100 Union noise nominations, depending on morphology. Preserve this failure evidence.

## Historical/current numerical comparison

| Abrupt500-window metric | Historical corrected report | Current replay |
| --- | --- | --- |
| Fixed post-step injections detected |5/10 |3/10 |
| First available new-regime STFT, window offset after shift |62 |134 |
| First mature available decision, offset after shift |201 |267 |
| Clean duty |88.0% |73.6% |
| STFT available duty |Not separately revalidated |428/500=85.6% |
| SplitB of final recalibration |150 |199 |
| Final calibration Ljung-Box p |0.43295 |about1.01e-8 |
| STFT noise nominations, legacy full denominator |7/487=1.437% |7/487=1.437% |
| Correct current available-noise denominator |Not recomputed |7/415=1.687% |
| Mature clean noise STFT |0 observed, old exposure |0/127;6/6 later signals missed |
| Glitch-reference refreshes |2 |3 |

One window is4 seconds. Stage labels, used threshold availability, clean exposure and sensitivity are different measurements. The old available time62 used a finite fallback with40 threshold samples; its1% target was unattainable. Current explicit abstention excludes72 noise windows. A matching legacy7/487 ratio must not conceal the changed valid denominator or lost sensitivity.

The other three suiteA aggregate rows and all eight suiteB aggregate rows match historical saved values within1e-12 comparison tolerance. Their per-window historical identity is not established: the old run has no input-array hashes. Current input arrays/reference regenerated exactly; source recipes and generator/noise code matchedd1ec9. Real-data suiteD, coherence suiteC and earlier phase2/initial phase2b were not rerun.

## What was and was not discovered

- Confirmed: a reproducible sensitivity-loss mechanism through contaminated adaptive threshold calibration, a failure hidden by mature/available labels and zero observed noise alarms.
- Confirmed: full-Union availability after a noise shift does not mean the frozen broadband background matches the new regime.
- Reproduced: stationary control1/200 STFT alarms and20/20 stationary transient injections; glitch-contamination guard retains one regime and zero glitch-reference refreshes; persistent drift limitations remain.
- Not established: a new astrophysical event, novelty of an anomaly, a novel detection method, population performance, independent confirmation or a safe automatic repair. No original physical data were examined and no external novelty review was performed.

## Next decision

Before claiming reliable adaptive recovery, design a separate protocol for background admission/segmentation and signal preservation. Evaluate no-signal tail extremes, repeated signals, drift, spectral resonance, available/clean exposure and held-out sensitivity jointly. Keep this baseline and its failures intact. Do not tune or select seeds to improve the reported score.

A research-prototype release may document these limits honestly; release notes must not say all adaptive sensitivity problems are resolved. No automatic code repair, commit/push, release or Zenodo change was performed.
