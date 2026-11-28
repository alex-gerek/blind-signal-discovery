# Saved-map temporal peak persistence — 28 November 2026

Descriptive analysis only. No detector, gate, calibration or production behavior changed. [Protocol](temporal_peak_persistence_protocol_2026-11-28.md), [full report](../results/temporal_peak_persistence/20261128T135237Z_633a21cf/report.md), [frequency pairs and counts](../results/temporal_peak_persistence/20261128T135237Z_633a21cf/pair_frequencies.csv), [case metrics](../results/temporal_peak_persistence/20261128T135237Z_633a21cf/case_metrics.csv).

## What was inspected

252saved1-second maps:216two-tone variants,24matched singles,12matched noise controls sharing12base noise/carrier/global-phase units.25frames/map,hop0.125s. All6300frames scanned across the full frequency axis for local maxima>=4 and pairs with an interior valley<=0.8*weaker peak. No accepted-cluster/anomaly/truth filter in blind extraction. Near2..3Hzpairs are a predeclared structural subset; all separations are also retained. The independent160+160negative maps were not saved and are not included.

Central9frames1.25..2.25s are used equally for signals and controls below. True-carrier matches (both errors<=0.5Hz) are evaluation-only annotations, not inputs to the tracker.

## Near the actual two tones

- Previously correct24cases: carrier-associated pairs occur in2..6central frames,median4.13/24have the SAMEfrequency pair in at least3consecutive frames.
- Previously unreturned187cases:140have at least one carrier-associated frame,median1(range0..6).24/187have an exact3-frame run;28/187if each frequency can shift<=1Hz per hop. These are not newly calibrated detections.
- Example88/90Hz: appears in5/9central frames but longest exact run3. With the same noise/carrier/global phase and second-tone amplitude0.5,88/90appears2/9and holds2frames; the original average returned no pair.

## Identical blind search in controls

For ANYnear2..3Hzpair anywhere on the map:19/24single-tone controls have a pair in at least one central frame;5have exact2-frame runs,1has3. Noise:7/12have at least one pair frame;2have exact2-frame runs,none3in this small set. One-tone carrier-anchored near-pair annotations show7/24with any frame and1with a3-frame exact run.

For unrestricted pair separation,ALL252records have at least one3-frame exact pair somewhere,including12/12noise controls. Permissive<=1Hz-step paths can last even longer in noise. Existence of stable arbitrary pairs is not a two-tone signature. These candidate patterns are not false-alarm rates of an adopted detector; controls are small and correlated across variants.

## Wrong peaks can be stable

Four earlier wrong134/139Hzpairs (true134.5/136.5) appear5,5,2,3of9central frames,with longest exact runs3,3,1,2. Their single-tone control also contains134/139in4central frames,longest2. Noise has139Hzalone above4in2central frames but lacks a qualifying134Hzpartner. Thus persistence does not eliminate the noise/leakage ambiguity.

The earlier wrong86/89pair (true86/88) lasts5central consecutive frames,6over the full map. Its second frequency is still wrong by1Hz. A persistence-only3or5frame rule would keep that error.

## Stability versus jumps

Exact trajectories repeat frequency bins; permissive trajectories allow1Hz per-component steps and may stitch branches. Per-frame strongest-pair selection is different: it can jump between unrelated frequency regions. For near pairs in the24correct group,8of79adjacent strongest-pair transitions move>1Hz in at least one component; among187previous no-pair cases66of236do. These aggregate shared variants,not independent transitions. Largest jumps reflect different selected structures,not physical oscillators moving hundreds ofHz.

One-second windows overlap87.5%. Three consecutive frames span0.25s between centers and1.25s of union support; five span0.5s and1.5s. They are not3or5independent confirmations. This descriptive longest-path study is retrospective,not streaming-latency or source-identity validation.

## Preservation and next decision

Separate verificationPASS:252maps,6300frames,205439retained pair-gate records,3024exact/jitter recurrence checks;529parent artifact hashes unchanged. Full frequencies,times,paths,gate lists,control comparisons and both figures retained. Figures visually inspected. No waveforms or maps regenerated,no tuning,exclusions,O3,commit,push,mainmerge,release orZenodo.

Possible next question: a truth-blind persistence rule that avoids both transient noise and stable wrong companions. It needs its own fixed operating criteria,independent calibration and new controls. No such rule is adopted or validated here; two peaks are not two reconstructed source waveforms or music stems.
