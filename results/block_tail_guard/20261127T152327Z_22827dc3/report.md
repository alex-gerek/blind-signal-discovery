# Whole-block tail guard: rejected feasibility pilot

Run `20261127T152327Z_22827dc3`, source d049e2b. One fixed experimental admission rule; unchanged detector core. No adaptive manager integration, production change or release. Full [protocol](protocol.md) and [independent verification](independent_verification.json) retained.

## Rule and design

On the preliminary calibration's Split-B scores z, quarantine the entire block iff max(z)>median(z)+10*max(Q95(z)-median(z),1). No individual extremes removed and no signal truth enters admission. Constant10 is an exploratory heuristic fixed before evaluation, not a significance threshold. Block rejection makes experimental STFT decisions unavailable; preliminary core candidates remain diagnostic. Whole-block score selection can still bias calibration, so admission is not a proof of correct FAP.

Three fresh seeds17101..17103 and four paired calibration cases: pink-only, pink plus four chirps, pink plus four bursts, and pink with rare impulsive NOISE. Each block265windows (66A/199B), evaluation80noise plus12paired-amplitude signal records (six base units) percase. Core FAP0.01/0.005 uses the same data/rule; the latter is exactly attainable at199B. Union atFWER0.005 is not attainable at199B and was only checked as a fixture, not evaluated as recovered. Cases and FAPs reuse inputs; counts must not be pooled as independent observations. Six signal base units per seed are reused across cases and amplitudes:18distinct target noise/clean units total, not72independent signal bases. At most three independent seed replications per case.

This is a one-proposed-block fixed-schedule pilot. There is no automatic reacquisition, release from quarantine or recovery-time measurement. A rejected block abstains throughout its evaluation. Results cannot validate an adaptive recovery policy.

## Results at FAP0.01

| Proposed calibration | Blocks admitted | Core signal nominations | Guard signal nominations | Guard signal abstentions | Lost core nominations |
| --- | ---: | ---: | ---: | ---: | ---: |
| Pink-only |3/3|22/36|22/36|0/36|0|
| Four calibration chirps |0/3|4/36|0/36|36/36|4|
| Four calibration bursts |1/3|18/36|6/36|24/36|12|
| Rare impulsive noise |0/3|13/36|0/36|36/36|13|

These are whole-window signal nominations, NOT known-support associations or source identifications. Each36contains paired SNR8/18 records. No confidence/superiority test was performed.

All three clean pink blocks pass, but the guard admits one contaminated burst block (seed17102:maximum11.0379206078953,limit14.199044362397185). It also rejects every impulsive-noise block, although these contain no added calibration signals. In those cases baseline noise alarms6/240 become0alarms with0available noise windows; the conditional rate is NULL, not zero. It loses13baseline signal nominations, rather than improving operational sensitivity. Pink-only noise remains3/240 with240available windows.

At FAP0.005 the same blocks are admitted/rejected. Core -> guard signal counts: pink20->20,calibrationchirps1->0,calibrationbursts12->4,impulsivenoise9->0 (each36paired-amplitude records). This supports the same rejection decision; these are dependent comparisons, not an independent confirmation.

## Decision

Predeclared clean-pink availability criterion passes. Complete contaminated-block quarantine and preservation of baseline recoveries both FAIL. Do not adopt, do not retune multiplier10 on these outcomes, do not call abstention an adaptive sensitivity repair. L09/L10 remain unresolved.

The useful outcome is a counterexample to this particular rule: an extreme-score guard does not reliably distinguish contamination by signals from legitimate impulsive noise, and whole-block quarantine alone cannot restore detection. This does not prove all possible admission methods impossible.

## Verification and limitations

Independent saved-data arithmetic checks pass24calibrations,2208evaluation records,270preliminary candidate records and4284regenerated input records. Reconstructed A median/Q95/spread, B maxima/order ranks, all evaluation scores/tail p-values/strict decisions, guard statistics, null fields, labels/exposure and operational losses.12focused fixtures include exact/nextafter guard boundaries, ties/invalid arrays and FAP/Union attainability;78direct uncached raw-feature parity checks pass. These are scoped experimental checks, not full pytest passes or scientific efficacy proof. Waveform regeneration verifies preservation, not correctness of the generator itself.

Candidates with core-calibration p-values remain explicitly DIAGNOSTIC when the guard abstains; official guard threshold/p are null. This does not repair production UI/export semantics or preserve all subthreshold seed clusters. Candidate checks verify the records produced by unchanged core; they do not establish physical signal correspondence. All input/source hashes retained, source unchanged. No measurement error/interruption, excluded case, retuning, real waveform/O3/data/manifest/network reads. Initial admissibility results, including failures, are retained in full.

Next proposal: study a separately acquired/provenance-qualified background control against accumulated-stream calibration, explicitly distinguishing what is available to an actual instrument from synthetic truth. Before an autonomous candidate, specify the acquisition/admission information it can really use and a release/recovery mechanism. Freeze a new protocol and independent sequence-level operational evaluation; do not enlarge this pilot until it passes.
