# Matched 2 Hz interventions — 28 November 2026

Source fc99a63915cc5137976ea493df1fb61fe6d803d6. Research only, unchanged production detector and unchanged diagnostic. [Protocol](matched_gap2_protocol_2026-11-28.md), [full report](../results/matched_gap2/20261128T115035Z_bf957291/report.md), [paired contrasts](../results/matched_gap2/20261128T115035Z_bf957291/paired_transitions.csv), [verification](../results/matched_gap2/20261128T115035Z_bf957291/independent_verification.json).

## New question

The previous5loose/4strict close-pair successes differed in noise, carrier and absolute phase. This follow-up holds those fixed within each of12freshbase units and varies only offset,relative phase and amplitude ratio. Gap2Hz,resultant root-energySNR18,duration1s,onset1.25s are fixed.18variants perbase=216matched cases; these are not216independent noise draws. Independent calibration and negative-control streams are separate; no thresholds or parameters changed after evaluation.

## Main results

Every two-tone case nominates an anomaly under both banks. Only the1splus-within-cluster-peaks arm estimates any pairs:29claims,24correct with both errors<=0.5Hz,25correct at<=1Hz;5/4wrong respectively.187cases have no eligible proposal. Other3arms have0pairclaims. This separates energy anomaly detection from frequency estimation; it does not reconstruct source waveforms.

| Fixed condition, W1D1 | Strictly correct /12base units |
| --- | ---: |
| Equal amplitudes, on-grid, relative phase0 | 12/12 |
| Same bases, equal amplitudes, off-grid, phase0 | 3/12 |
| Equal amplitudes, on-grid, phase90degrees | 4/12 |
| Equal amplitudes, on-grid, phase180degrees | 1/12 |
| Ratio0.5, on-grid, phase0 | 0/12 |
| Same bases, ratio0.5, off-grid, phase0 | 4/12 |

All other cells are retained in the report/heatmap. Off-grid estimates total7strict successes across108variants (12sharedbases); on-grid17/108. The fact that off-grid works in some matched cases corrects any reading of the older on-grid-only successes as a universal requirement. Four unequal-amplitude off-grid successes also refute treating equal amplitudes as necessary. These are interactions under the specified estimator, not a universal phase/amplitude law.

## Controls and limitations

Every arm has0false pair claims on160independent noise and160independent single-tone inputs. For each0/160, the two-sided95%upper bound is2.28%; no1%certification or zero population error claim. W1has2/160noise anomaly nominations; W0has0/160. Matched controls12noise and24single are reported separately and not pooled.

Fixed resultant energy changes individual component amplitudes with phase/ratio. Moving both carriers0.5Hz in colored noise changes their local noise level as well as alignment. Window comparison changes the frequency-median smoothing footprint. Strong fixed energy, one onset and12base draws do not establish robustness to weak mixtures, edge timing, other durations or real sensors.

No main-case failure occurs at anomaly acceptance or pair cutoff; absence of eligible structural/peak proposals is the bottleneck at the recorded stage. This does not isolate the causal share of STFT, normalization, clustering or the peak/valley gate. The next useful step is inspection of matched profiles/gates on the SAMEbase, not another uncontrolled sweep or production change.

Independent checks PASS:2calibrations,504trace maps,2288claimchecks,864matchingchecks,572regeneratedevaluation inputs,132pairedcontrasts.8direct production parity checks passed; old matching fixtures and216generator pairing/energy checks passed before fitting. These are scoped checks, not pytest passes or a scientific correctness proof. Timing and runtime with conditional-subset caveats are in the full report; strict successes have median absolute onset/end box errors0.1875s and midpoint error0s. Figures visually inspected. No run failures, excluded cases, retuning or reruns.

Local gzip records and504NPZmaps retain all outcomes; large binary traces remain gitignored with hashes. Production and previous evidence hashes unchanged. No O3/data/manifest access, scientific data download, commit/push, main merge, release or Zenodo. Results await review.
