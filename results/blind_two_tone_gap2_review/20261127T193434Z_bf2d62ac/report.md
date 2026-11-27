# Five successful gap-2 Hz cases versus 67 unsuccessful cases

Retrospective diagnosis of the frozen W1D1 (1 second plus within-cluster peaks) pilot, not a new experiment or independent validation. Success uses both errors <=1 Hz; strict matching uses <=0.5 Hz. Inference never received true frequencies or labels. Relative phase is at signal onset and is not constant over time for unequal frequencies. SNR means linear root-energy of the summed addition, not either component or decibels.

| Case | Resultant SNR | Weaker/stronger amplitude | Relative phase at onset | True frequencies (Hz) | Estimates (Hz) | Pass 0.5Hz? |
| --- | ---: | ---: | --- | --- | --- | --- |
| 30701_0072 | 8 | 1 | 0 rad | 116, 118 | [116.0, 118.0] | True |
| 30701_0073 | 18 | 1 | 0 rad | 174, 176 | [174.0, 176.0] | True |
| 30701_0075 | 18 | 1 | 0 rad | 116, 118 | [116.0, 118.0] | True |
| 30701_0091 | 18 | 0.25 | 0 rad | 148, 150 | [148.0, 151.0] | False |
| 30701_0097 | 18 | 1 | 1.5708 rad | 104, 106 | [104.0, 106.0] | True |

| Factor | Value | All cases | Successful | Unsuccessful | Strict successful |
| --- | --- | ---: | ---: | ---: | ---: |
| snr | 8 | 36 | 1 | 35 | 1 |
| snr | 18 | 36 | 4 | 32 | 3 |
| ratio | 0.25 | 24 | 1 | 23 | 0 |
| ratio | 0.5 | 24 | 0 | 24 | 0 |
| ratio | 1 | 24 | 4 | 20 | 4 |
| phase | 0 | 24 | 4 | 20 | 3 |
| phase | 1.5707963 | 24 | 1 | 23 | 1 |
| phase | 3.1415927 | 24 | 0 | 24 | 0 |
| offset | 0 | 36 | 5 | 31 | 4 |
| offset | 0.5 | 36 | 0 | 36 | 0 |
| onset | 1.25 | 36 | 3 | 33 | 3 |
| onset | 1.3115234 | 36 | 2 | 34 | 1 |

Four of five loose-tolerance successes have resultant SNR18; one has SNR8. Stronger energy helps descriptively but is insufficient:32/36strong examples still fail pair estimation, although all36nominate an anomaly. Weak examples nominate15/36anomalies and identify one pair.

Four of five have equal amplitudes. The only unequal-amplitude success is ratio0.25, truth148/150Hz estimated148/151Hz. It sits exactly at the loose1Hz tolerance boundary and fails0.5Hz matching. All four strict successes therefore have equal amplitudes. Ratio0.5 gives0/24; the isolated ratio0.25 loose success must not be interpreted as weaker partners helping.

Four loose successes start at relative phase0 and one at pi/2; none at pi. Strict successes are three at0 and one atpi/2. Counts are marginal descriptions of a small factorial matrix, not a causal phase law. Carriers, absolute phases and independent noise draws also differ between examples.

All five loose successes, and all four strict successes, use carriers on both frequency grids (offset0); none of36off-grid examples succeeds. Exact-grid dependence is a central observed limitation of bin-center peak estimates, not proof of a sole cause. Base offset0.5Hz is half the1Hz grid and a quarter of the2Hz grid. Both onset alignments have successes (3/36aligned versus2/36half-hop); onset alignment alone is not sufficient.

Of the67unsuccessful examples,21fail the anomaly gate and46nominate an anomaly without returning a pair. No wrong pair was returned within these72gap2Hz cases. This separates detectability from pair estimation.

Conclusion: the small set of strict successes is concentrated on equal amplitudes and exact-grid carriers, mostly stronger mixtures, with favorable sampled relative phases. These observations motivate matched follow-up controls; they do not isolate causes because independent noise/carrier/absolute-phase draws changed together. Do not tune on these cases and call a repeat independent validation. The historical detector, calibration, thresholds and artifacts remain untouched. No waveform runs, O3 access, commit or push.
