# Fixed temporal-envelope pilot — 27 November 2026

[Full report, paired counts and figures](../results/temporal_envelope_pilot/20261127T101745Z_0fdf2942/report.md) · [Protocol](temporal_envelope_pilot_protocol_2026-11-27.md) · [Independent verification](../results/temporal_envelope_pilot/20261127T101745Z_0fdf2942/independent_verification.json).

One research-only change: centered temporal median kernels257/33/9frames, approximately1s, instead of31/15/7. Original spectral estimation, Z>4/eight-neighbor/two-cell rule unchanged; no gap-linking. Each method has its own full calibration. No parameter sweep or production adoption.

On three fresh stationary-calibrated synthetic realizations, stationary noise alarms5/600->8/600, smooth-drift10/600->15/600, step5/600->12/600. Weak chirp associations19/30->21/30,18/30->20/30,21/30->22/30 respectively. Weak burst associations1/30->2/30,5/30->8/30,4/30->8/30. Strong chirps stay30/30; strong bursts30/30->30/30,26/30->30/30,29/30->30/30. Step weak-chirp decisions include one lost baseline recovery despite a net gain. Paired conditions reuse underlying data; drift is out-of-distribution, not nominalFAP validation. These small counts do not prove a general improvement/degradation.

On the same stationary SNR18burst pairs, short-scale median temporal response falls2.165689->1.110063. Tracking proxy improves for stationary noise but is mixed under drift/steps. A smaller average tracking error is not a false-alarm guarantee.

GW170814 H1 remains missed, L1 retained. H1 same-peak near-reference cluster grows20->43cells and raw123.224760->252.836460, but its own recalibrated raw boundary rises192.181715->290.480625; score2.612768<threshold3.160593, p7/496=.01411290. No old-threshold recovery claim or known-event tuning. Off-source post-veto H1 1->0/L1 3->2 over151/site are nominations, not pure-noise FAR.

All10calibrations,4,960arm/window decisions,1,931accepted candidates,3,120synthetic input records,10,800tracking records and2,160paired-factor records are independently checked, with90manual map-arm checks and the predefined tracking/paired subset. No measurement failures/interruptions or retuning; hashes match. Counts are not pytest passes or independent replications.

Decision: retain as an experimental trade-off, do not replace the detector. Next validation needs a task-specific alarm tolerance and independent evaluation at a common alarm budget, not a larger smoothing/threshold sweep until H1 passes. No O3, new event, production/app/test change, main/release/Zenodo. DNA-window and L09/L10 remain separate.
