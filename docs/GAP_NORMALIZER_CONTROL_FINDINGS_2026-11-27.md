# Gap-link and paired-normalizer controls — 27 November 2026

[Full report and figures](../results/gap_normalizer_controls/20261127T095901Z_b912adbb/report.md) · [Protocol](gap_normalizer_control_protocol_2026-11-27.md) · [Independent verification](../results/gap_normalizer_controls/20261127T095901Z_b912adbb/independent_verification.json).

## Separately calibrated gap prototype: no observed improvement

H1 remains missed and L1 detected after each method's own full calibration. Event raw maxima at all three scales are unchanged; the one-frame link does not enlarge the strongest event components. Off-source post-veto nominations change H1 1->2 and L1 3->5 over151windows/site. These are small developmental real-data counts, not certified noise-only FAR.

On three fresh synthetic realizations, original/prototype both recover17/30SNR8 and30/30SNR18 injections, with identical paired injection decisions. Noise nominations are6/600 versus7/600 (one extra gap-only window). The changed synthetic distribution prevents comparison with earlier13/30Workstation sensitivity as if it were a version improvement. No statistical general-inferiority claim; no prototype adoption or gap-parameter sweep.

## Paired normalizer control: a measured software response

Added signals can raise the locally estimated temporal envelope. For pink-noise broadband bursts the short-scale median T1/T0 is2.086649atSNR18 and4.753310atSNR36; the temporal-only frozen-noise diagnostic satisfies Z_frozen=Z_production*(T1/T0) pointwise. For pink chirps the medians are1.240820/1.312543. Responses vary; long-scale decreases and all360zero-amplitude identity controls are preserved.

This establishes a response to constructed added signals and isolates its division arithmetic. It does not establish the physical cause of H1's miss, a detection-efficiency loss, or a deployable replacement normalizer. The exact noise-only partner is an oracle unavailable on real event-containing data. A future robust envelope estimator needs independent calibration and drift/noise/sensitivity evaluation.

The two tests completed with unchanged production code and preserved original input/source hashes. Verification checks10calibrations,1,960arm/window decisions,362accepted candidates,2,100input-regeneration records,1,440normalizer records and48independent manual map reconstructions. These counts are not pytest passes or independent scientific replications. One supplemental export-permission failure/retry is retained separately. No O3/new real event/detector adoption. DNA-window and adaptive L09/L10 remain separate.
