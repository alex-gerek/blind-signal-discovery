# Blind two-tone exploratory findings — 27 November 2026

Source: cc68b93fee699307e3314fce1bd4c9fe403c1e6e. Production detector unchanged. Thresholds and methods frozen before held-out evaluation. This is a synthetic exploratory frequency-estimation study, not source reconstruction or music stem separation.

## Results

| Analysis / research diagnostic | Anomaly / 432 two-tone inputs | Pair claims | Correct pairs (both errors <= 1 Hz) | Wrong pair |
| --- | ---: | ---: | ---: | ---: |
| 0.5 s / disjoint clusters | 317 | 49 | 49 | 0 |
| 0.5 s / plus within-cluster peaks | 317 | 77 | 76 | 1 |
| 1 s / disjoint clusters | 307 | 68 | 67 | 1 |
| 1 s / plus within-cluster peaks | 307 | 98 | 97 | 1 |

Each method produced zero false pair claims on 160 independent noise-only and 160 independent one-tone observations. These control populations are separate and methods are paired; do not pool them to inflate a false-rate denominator. The exact two-sided 95% upper bound for 0/160 is approximately 2.28%, so a 1% false-pair budget is not empirically certified. The 0.5 s bank still produced two noise anomaly nominations; anomaly and pair decisions differ.

At 2 Hz separation only 1 s plus within-cluster peaks found correct pairs: 5/72 with 1 Hz tolerance, 4/72 with stricter 0.5 Hz tolerance. At 1 Hz separation all arms had zero correct estimates. The best arm was far more effective at stronger resultant energy: 4/216 correct at SNR 8 versus 93/216 at SNR 18. These are linear root-energy SNRs, not decibels or matched-filter SNR.

The full factorial uses gaps 1, 2, 4, 8, 16, 96 Hz; amplitude ratios 1, 0.5, 0.25; relative phases 0, pi/2, pi; on/off-grid carriers; two onset alignments; and resultant SNR 8/18. Each of the 432 cells has one independent noise/carrier draw. Unlike the historical study, individual components are not each SNR 18: the SUM is energy-normalized. Lower far-pair counts are therefore not a like-for-like regression against historical 100/100.

No frequencies, source count or onset truth are supplied to inference. Exact duplicate estimates are invalid; duplicate proposals are consolidated; at most one pair is returned. Evaluation checks one-to-one matching, both permutations and inclusive tolerance. Any pair on noise or one tone is a false two-tone claim. Calibration uses independent seeds and separate per-configuration scores/thresholds. Pair calibration has sparse eligible proposals; all four pair cutoffs equal the missing-proposal sentinel (-infinity), with explicit strict-exceedance semantics. It is not an unavailable model or a measured absence of false alarms.

## Trade-offs

With within-cluster diagnostics in both banks, moving from 0.5 s to 1 s gained 28 correctly estimated examples but lost 7; 69 succeeded in both. Do not infer universal dominance. Frequency smoothing bandwidth changes with the FFT grid even though the median kernel parameter is held fixed, so window length is not the only effective representation change.

Among correctly estimated pairs (different selected subsets), median absolute onset/end box errors are 0.0625 s for the 0.5 s arms and 0.1875 s for the 1 s arms. These conditional figures are not a paired timing superiority test. Median analysis-bank times were about 17.94 versus 20.07 ms per four-second input on this machine. Additional within-cluster proposal extraction medians were approximately 0.125 versus 0.166 ms; shared re-labeling overhead is not included in those separate timers. No true cold-process, online latency or deployment memory guarantee was measured.

The separate record-edge supplement has 24 truncated pairs. Correct-pair counts were 2, 5, 1, 2 respectively. All four nominated anomalies in 24/24. This exposes a difference between noticing energy and estimating two carriers. Edge clipping reduces observed energy and changes support; no out-of-record waveform is recovered.

## Checks and retained failures

The first attempt failed before calibration/evaluation because an ambiguous-assignment fixture was incorrectly specified; its complete failure folder and driver are retained. Only that fixture changed, and the completed run uses a fresh output folder. No held-out tuning, excluded signals, threshold changes or waveform-study reruns.

Independent verification passed 2 calibrations, 796 pair-calibration records, 1,552 evaluation records/maps, 4,696 proposal sets, 3,104 claim checks and 1,824 matching checks. Further checks reproduced 1,458 unique waveform hashes, checked six repeated input-only inference cases with altered evaluator labels, and checked production parity on ten records. These are scoped checks, not a full pytest suite or evidence of general scientific correctness. Production and historical artifact hashes are unchanged. Figures were visually inspected.

## Evidence

- [Frozen pre-run protocol](blind_two_tone_exploratory_frozen_2026-11-27.md)
- [Full report](../results/blind_two_tone_exploratory/20261127T192047Z_f279641e/report.md)
- [Marginal strata](../results/blind_two_tone_exploratory/20261127T192047Z_f279641e/stratified_counts.csv)
- [Complete factorial cases](../results/blind_two_tone_exploratory/20261127T192047Z_f279641e/factorial_cases.csv)
- [Verification](../results/blind_two_tone_exploratory/20261127T192047Z_f279641e/independent_verification.json)
- [Reproducibility](../results/blind_two_tone_exploratory/20261127T192047Z_f279641e/reproducibility_verification.json)

The 1,552 compressed longest-scale trace files occupy approximately 508.5 MiB in the persistent run folder; existing gitignore excludes NPZ files. A SHA256 index records every local artifact. The archived generator/source and metadata support deterministic regeneration; publication would require an explicit storage decision. No O3 access, production change, commit, push, merge, release or Zenodo action.

## Decision

Retain all arms as research evidence and keep the production detector unchanged. The experiment shows limited blind frequency-pair estimation in some mixtures; most close pairs remain unresolved or unreturned. A further step should inspect successful and failed close-pair cases, independently replicate calibration and broaden single-tone/noise controls before proposing adoption. Two frequency estimates are not two reconstructed source waveforms.
