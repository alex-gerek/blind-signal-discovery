# Adaptive calibration contamination — exploratory paired follow-up

Source043f0eb; parent replay20261126T163629Z_75c7a756. Protocol was saved after observing the parent failure and before this intervention. No detector parameter or source change.

## Confirmed observations

| Measurement | Original injected sequence | Four earlier injections replaced with exact original noise |
| --- | --- | --- |
| Injected members in threshold SplitB |330,340,350,360 |none |
| SplitB size |199 |199 |
| Mature threshold |163.839400941 |1.65803953491 |
| Mature first used window |367 |367 |
| Same six later injections detected |0/6 |5/6 |
| Later mature noise false alarms |0/127 |0/127 |
| Calibration Ljung-Box p |1.0113175e-08 |0.54505621 |

All samples outside the four intervention windows were asserted bitwise identical. The baseline's complete500-window outputs matched the parent exactly. The six later scores were EXACTLY identical across arms. Therefore the detection change is explained by changed thresholds, rather than changes in those later waveforms/scores or recovery timing.

The original sequence fed its own known signals into the estimated noise tail. In this constructed example, a conservative contaminated threshold hid five otherwise detectable later injections. Mature/available does not imply recovered sensitivity. The old resonance-only explanation does not account for this controlled difference, although spectral masking may still matter for other misses.

## Scope and practical impact

This is a confirmed mechanism in one post-hoc synthetic paired experiment, not proof of its frequency in real streams or a novel physical discovery. The control used oracle injection labels and is NOT a deployable fix. Simply deleting high-score windows can bias the noise calibration and break false-alarm control; do not adopt that rule without a separate design/protocol.

Keep the baseline and failed0/6 result. Before advertising reliable adaptive recovery, design and independently validate admission/segmentation for recalibration, signal retention and uncertainty. Include stationary/no-signal extremes, recurrent true signals, drift and shifted_ar2, with held-out sensitivity and noise exposure. Preserve all unavailable decisions and candidate diagnostics.

The frozen broadband/Union false alarms in the parent are a separate limitation:127/127 mature noise windows were nominated. They do not establish astrophysical signals or recovered1% Union FWER after a changed background.

No original data, O3, commit, push, release or Zenodo operation.
