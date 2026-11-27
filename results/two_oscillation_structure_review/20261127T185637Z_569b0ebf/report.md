# What counted as two separate structures?

Retrospective supplement to run `20261127T183750Z_cc70b1c7`. The original report, protocol and artifacts are unchanged. This review examines existing accepted candidates; it does not regenerate waveforms, retune settings or isolate causal contributions experimentally.

## Exact rule

A unit counts as a two-structure success when there exist two different **accepted candidates in the same STFT scale** such that:

1. Each candidate's time bounding interval intersects the evaluation interval [1.25, 2.25] seconds. End-point contact counts.
2. One frequency bounding interval intersects the known A carrier neighborhood [fA - 4, fA + 4] Hz; the other intersects the partner neighborhood [fB - 4, fB + 4] Hz. End-point contact counts here too.
3. Their frequency bounding intervals are strictly disjoint: `a.f_max < b.f_min` or the reverse. Touching frequency boxes do not count as separate.

This is an existence test over scales; there is no requirement to succeed at all three scales. Candidates from different scales cannot be combined to make a pair. Time overlap between the two candidates is not required beyond each touching the evaluation interval. The criterion does not inspect the exact cluster shape inside its box, and is not a source-count estimator or proof that either waveform was recovered.

Carrier neighborhoods are evaluator knowledge: the detector never receives them. At a 2 Hz separation the +/-4 Hz neighborhoods overlap substantially. One candidate intersecting both neighborhoods is therefore weak evidence and does not count as two structures. No peak-fitting, valley-depth test, deconvolution or super-resolution estimator was used.

Code: `experiments/run_two_oscillation_anomaly.py`, functions `intersects` and `diagnostics`; independently recomputed in `experiments/verify_two_oscillation_anomaly.py`, functions `touch` and `two`.

## Exact settings

Sampling: 1,024 Hz. Observation: 4 seconds (4,096 samples). Two finite tones share a one-second Tukey envelope (alpha 0.15), nonzero support within [1.25, 2.25) seconds. A carrier is an integer 96..176 Hz; B is A+96 Hz or A+2 Hz. Phases are independently drawn. Each component has linear root-energy SNR 18 relative to unit noise sigma 1; the mixture is not rescaled. Colored noise: alpha 1, knee 16 Hz. Five seeds 22101..22105, 20 paired units per seed. Relative phase, colored noise and resultant energy affect the measurement.

| Scale | Hann window samples / duration | Overlap samples | Time hop | Frequency grid spacing | Frequency / time median kernels |
| --- | --- | ---: | ---: | ---: | --- |
| short | 32 / 31.25 ms | 28 | 3.90625 ms | 32 Hz | 5 / 31 |
| medium | 128 / 125 ms | 96 | 31.25 ms | 8 Hz | 11 / 15 |
| long | 512 / 500 ms | 384 | 125 ms | 2 Hz | 21 / 7 |

These are `DEFAULT_V02_SCALES` in `src/bsd/detector_v2.py`. STFT uses no boundary extension or padding and its FFT size equals the segment length (no `nfft` override); DC/Nyquist are removed. The four-second observation is not a four-second Fourier window. Grid spacing is not by itself resolving power: the Hann spectral response spreads each tone across neighboring bins, and finite duration/taper also broadens it.

`src/bsd/detector.py`: median subtraction, temporal-median spectral baseline and frequency-median smoothing; frequency-median temporal envelope and time smoothing; normalized excess power Z. Cells with **Z > 4** seed clusters, using eight-neighbor connectivity (diagonals included). At least two cells are required. Cluster raw score is sum(Z - 1). Frequency boxes extend half a frequency bin past the outer cluster-bin centers; time boxes extend half a hop past the outer time centers. Thus box separation is stricter than merely different peaks or disconnected cell sets.

Each seed has 284 independent calibration windows: first 85 normalize per-scale scores by `(raw - median) / max(Q95 - median, 1)`; remaining 199 calibrate the maximum standardized score over all scales. Target window-level FAP is 0.005. Accepted candidates must exceed that frozen global threshold and have empirical p <= 0.005. With 199 null maxima the minimum p is 1/200; the threshold is the largest calibration null maximum. Thresholds differ by calibration seed:

- 22101: 5.86927757893462
- 22102: 2.35554586520168
- 22103: 3.28665942743702
- 22104: 3.84699375852354
- 22105: 2.2556471296611

These are decision thresholds, not the Z > 4 cell-seeding threshold. No veto, adaptive update, warning gate or component template is used in this pilot.

## What the saved candidates actually show

All denominators below are 100 paired units. "Any accepted" refers to any candidate in that scale; the other columns use the carrier/time diagnostic above.

| Mixture | Scale | Any accepted | One box touches both neighborhoods | Two disjoint boxes |
| --- | --- | ---: | ---: | ---: |
| A_plus_B_far | short_31ms | 7 | 2 | 0 |
| A_plus_B_far | medium_125ms | 100 | 0 | 100 |
| A_plus_B_far | long_500ms | 100 | 0 | 100 |
| A_plus_B_close | short_31ms | 94 | 94 | 0 |
| A_plus_B_close | medium_125ms | 100 | 100 | 0 |
| A_plus_B_close | long_500ms | 100 | 100 | 0 |

Consequently the original 100/100 far result is supported independently by both medium and long scales. The close result is not failure to nominate: medium and long nominate 100/100, with one box touching both overlapping neighborhoods in every unit, but no qualifying disjoint pair at any scale.

Concrete first unit (seed 22101): carriers 96 and 192 Hz produce accepted long-scale frequency boxes [93, 101] and [189, 197] Hz. Carriers 96 and 98 Hz instead produce one accepted long-scale box [93, 105] Hz, spanning 1.1875..2.3125 seconds. Its accepted medium-scale box is [84, 108] Hz. All candidate records are preserved in counts.json and the original measurements, including incidental candidates.

## What limited the close case — and what remains unproved

The finest available grid is 2 Hz, so the close tones are only one grid interval apart there, and less than one interval apart at the shorter scales. Hann broadening and finite-duration leakage overlap their representations. Eight-connected above-threshold cells can then become one component; even separate components would fail our diagnostic if their expanded frequency boxes touch or overlap. The observed accepted boxes show that the current representation, clustering, acceptance and evaluation procedure did not yield a qualifying pair.

This is a mechanistically consistent explanation, not an ablation proving that any single setting caused all 100 failures. This review does not inspect every rejected raw cluster or attribute losses between seeding, connectivity, score acceptance and box geometry. Normalization, noise, phases and unequal resultant energies also remain relevant. We have not established that the information is absent from the waveform, that all methods fail at 2 Hz, or that increasing the observation/window length would fix the problem.

A future causal study would need a predeclared comparison separating raw peaks, seeded components, rejected/accepted candidates and box criteria, with independent calibration whenever the detector settings change. It should vary one factor at a time and retain false alarms. No such experiment is claimed or run here.
