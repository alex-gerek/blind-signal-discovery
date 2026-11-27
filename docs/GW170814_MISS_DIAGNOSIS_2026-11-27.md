# GW170814 miss diagnosis: confirmed mechanisms and limits (2026-11-27)

[Protocol](gw170814_miss_analysis_protocol_2026-11-27.md). Source30849e8. [Machine-readable decomposition](../results/gw170814_miss_analysis/20261127T061905Z_9a9108ec/diagnostic_results.json). No new detector or calibration runs: read only the already saved records and exact three processed event-window arrays. No O3, new event, retuning or masking.

## Main correction: a window maximum is not necessarily the known event

| Site | Largest standardized score anywhere in4s | Largest score satisfying the original near-reference gate | Frozen threshold | Decision |
| --- | --- | --- | --- | --- |
| H1 |1.144807 |1.144807 |2.590544 |Miss, before veto |
| L1 |8.892396 |8.892396 |2.240084 |Detected; survives veto |
| V1 |3.619174 |0.664365 |4.006891 |Miss, before veto |

The V1 global winner is a long500ms-scale narrowband cluster spanning local2.4375..3.8125s,369..377Hz, peaklocal3.125s/374Hz. Reference local time is1.5s; the peak is1.625s late. It does NOT satisfy the original near_reference Boolean. The strongest V1 medium-scale cluster (z3.411794,peak3.15625s/376Hz) is also away from the reference. The largest near-reference V1 cluster has z0.664365 on the125ms scale.

Therefore the previous statement that V1's window p=5/496=0.01008065 is close toFAP0.01 must not be interpreted as almost recovering GW170814. That rank belongs to the off-time window winner. Lowering the threshold after seeing the known event could nominate that other structure without demonstrating recovery of the event. Its physical origin is not established here.

The fixed diagnostic gate is the PARENT gate: interval intersects local[1.0,2.0]s and35–250Hz. It is loose association, not identification of a gravitational wave or a newly calibrated ROI detector. No new regional p-value/significance is claimed.

## H1: the failure is after seed clustering and before candidate acceptance

H1 has131 seed-qualified clusters in the full window,21 satisfying the reference gate. Its strongest short-scale cluster is near the reference, local1.501953..1.537109s,48..208Hz,20bins; raw score123.224760 becomes z1.144807.

The conversion is z=(raw-background_median)/background_spread. H1's short-scale values are median68.621237,spread47.696722. The fixed z threshold2.590544 corresponds to raw score192.181715;123.224760 is insufficient. Medium-scale H1 raw36.047195 is below its raw boundary88.317264. Its strongest long-scale cluster is off-time/498Hz and even its standardized score is negative; the near-reference long-scale cluster is weaker.

This confirms the software mechanism without asserting that any particular raw cluster is the physical event. H1 was not unavailable, DQ-dropped or vetoed into a miss. Fragmentation/spectral response/noise may help explain its weak cluster statistic physically, but no causal intervention isolated these factors.

## L1: all three scales cross their frozen thresholds

L1 standardized maxima are8.580261(short),8.892396(medium),2.778833(long), all above2.240084. The maxima satisfy the reference gate and accepted candidates survive the original vetoes.

On the medium scale, H1/L1 background medians are31.827591/31.436325 and raw acceptance boundaries88.317264/87.882687, but observed raw maxima are36.047195/255.509778. Thus the directly observed difference is a much stronger clustered statistic in L1, not just a slightly different global threshold. This arithmetic comparison is not proof of the physical origin of that difference and does not justify using one site's calibration for another.

## V1: heterogeneous background and off-time narrowband structure

V1 normalization spreads are175.296(short),103.824(medium),219.478(long), larger than H1/L1's corresponding spreads. A large raw955.333 long-scale score becomes only z3.619174 against that background; the required raw boundary is1040.428.

Saved PSD-fit block median ASD over35–250Hz:
- H1:7.51948e-24 strain/sqrt(Hz).
- L1:7.27526e-24 strain/sqrt(Hz).
- V1:2.12770e-23 strain/sqrt(Hz), about2.83xH1 and2.92xL1.

These are descriptive spectral medians of a separate early PSD block, not event SNR or confidence intervals. They do not establish how much the actual wave projects onto each detector, nor isolate why a given event cluster is weak. The late374Hz feature must not be labelled instrumental/astrophysical solely from its shape.

## Figures and checks

[Score decomposition](../results/gw170814_miss_analysis/20261127T061905Z_9a9108ec/scale_scores.png).
[Full-window winner locations](../results/gw170814_miss_analysis/20261127T061905Z_9a9108ec/full_window_winners.png).
[Saved background ASD](../results/gw170814_miss_analysis/20261127T061905Z_9a9108ec/background_asd.png).
[Raw acceptance-boundary table](../results/gw170814_miss_analysis/20261127T061905Z_9a9108ec/scale_decomposition.csv).

All472 saved raw cluster standardizations were independently recomputed to absolute tolerance1e-12; every per-scale/global maximum agrees with saved summaries. Input/source records and three event-array hashes match. Full-window plots display30–400Hz; gray/white reference gate and red global-winner boxes are diagnostic only. Where the same cluster wins both categories, outlines overlap. Shared visible-bin5/99.5percentile contrast does not alter data/detection.

No physical response causal explanation, new significance, detection criterion or software repair is established. Existing L09/L10 adaptive-calibration limitations are unaffected. Data/figures derive from GWOSC/LIGO/Virgo measurements underCC BY4.0; see the parent report for attribution. Original/local waveform arrays remain outsidegit.

## Research conclusion and next work

The saved replay now supports a more precise map: H1's strongest nearby structure does not reach the accepted background-tail score; L1 crosses on all scales; V1's apparent near-threshold maximum is from a later different-frequency feature, while nearby structures are much weaker.

Before making claims about causes of misses or overall efficiency, predefine any further diagnostic comparisons (window alignment, spectral conditions, fixed resolution) and calibrate/report all variants rather than selecting the best-looking output. This stage changes interpretation, not the detector, and cannot establish physical novelty.
