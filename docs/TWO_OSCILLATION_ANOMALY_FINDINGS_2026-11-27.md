# Two oscillations: detecting an anomaly versus distinguishing structure

27 November 2026. Source revision: cc68b93fee699307e3314fce1bd4c9fe403c1e6e. Research evidence only; no detector change adopted.

## Question and method

Can the unchanged detector flag unfamiliar structure in mixtures without receiving the component waveforms or their frequencies? Five fresh calibration seeds and 100 paired mixture units were used. Calibration and 1,000 held-out noise windows were independent. Ground truth was used only for evaluation. Each component had linear root-energy SNR 18; sums were not renormalized. These are strong synthetic controls, not real-data validation or an equal-resultant-SNR resolution experiment.

## Observations

| Condition | Anomaly nominations | Two disjoint candidate structures in the same scale |
| --- | ---: | ---: |
| Frequency separation 96 Hz | 100/100 | 100/100 |
| Frequency separation 2 Hz | 100/100 | 0/100 |
| Exact cancellation A minus A | 0/100 | 0/100 |
| Near cancellation A minus 0.95A | 0/100 | 0/100 |

Single-component controls and coherent additions were nominated in all 100 units per case. Independent held-out noise produced 2/1,000 nominations. This sample does not certify a universal false-alarm bound. Cases share noise within a unit and must not be pooled as independent trials.

Exact cancellation leaves the same bytes as the paired noise: there is no observable addition to detect. Near cancellation leaves a residual of SNR 0.9; its absence of nominations is a bounded sensitivity result. A plus A and a single component 2A also yield identical input and complete core output, so the record cannot uniquely identify those two source descriptions.

The two-structure criterion describes candidate morphology, not recovered sources. Loose carrier neighborhoods can both intersect one candidate, even in a single-component control. Two fixed frequency gaps do not establish a universal frequency-resolution boundary. Finite-duration leakage, colored noise, phase and resultant energy also matter.

## Verification and evidence

The independent checking script passed five calibrations, 2,000 measurement records, 2,565 candidate records, 3,420 regenerated inputs and 75 direct uncached parity checks. It checks calibration arithmetic, decisions, candidate bounds, morphology, input identities and energy cross terms. These are scoped checks, not a full pytest run or proof of scientific validity. Both figures were visually inspected. No measurement failure, excluded seed, retuning or waveform rerun occurred.

- [Pre-execution protocol](two_oscillation_anomaly_protocol_2026-11-27.md)
- [Full run report](../results/two_oscillation_anomaly/20261127T183750Z_cc70b1c7/report.md)
- [Independent verification](../results/two_oscillation_anomaly/20261127T183750Z_cc70b1c7/independent_verification.json)
- [Counts](../results/two_oscillation_anomaly/20261127T183750Z_cc70b1c7/summary_counts.csv)
- [Count figure](../results/two_oscillation_anomaly/20261127T183750Z_cc70b1c7/mixture_counts.png)
- [Representative maps](../results/two_oscillation_anomaly/20261127T183750Z_cc70b1c7/mixture_maps.png)

## Decision

Continue template-free anomaly research while distinguishing nomination, candidate morphology and source identification. This pilot adds no source separator and makes no physical discovery. The next useful study should predeclare a wider range of unseen structures and noise conditions, freeze the detector and thresholds, and retain false alarms and misses. No production, application, O3, release or Zenodo changes; results remain local, uncommitted and unpushed.


## Follow-up: precise definition and settings

See [structure diagnostic supplement](../results/two_oscillation_structure_review/20261127T185637Z_569b0ebf/report.md) for the exact rule, STFT settings, per-seed thresholds, per-scale counts and limits of the causal explanation. No new waveform experiment or detector change.
