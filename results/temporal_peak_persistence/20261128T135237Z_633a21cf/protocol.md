# Saved-map temporal peak persistence — descriptive protocol

28 November 2026. Owner requested temporal frequency-pair stability on two-tone,one-tone and noise controls WITHOUT detector changes. Consult BSD-23/24: mean spectra can suppress transient peaks, and one repeated139Hzfeature gives wrong pairs. New question is persistence/drift on the SAME saved maps, not a new signal-search benchmark. All viewed inputs are development evidence. No new waveform generation, map reconstruction, threshold fitting, detector deployment, O3, commit or push.

## Input and availability boundary

Read only W1(1-second) NPZ traces from matched run20261128T115035Z_bf957291:216two-tone variants,24matched single-tone controls,12matched noise controls;252maps share12independent noise/carrier/global-phase base units. Complete maps have511interior frequency bins1..511Hz,25frame centers0.5..3.5s,hop0.125s. Parent phase/ratio/grid matrix and energy normalization unchanged. The160independent noise and160independent one-tone negative-control maps were NOT saved; they are not part of this temporal analysis. Analysing them later would require separately recorded exact-hash regeneration/unchanged map computation. Do not represent36matched controls as320independent controls.

## Fixed truth-blind per-frame diagnostics

For EVERY saved frame over the full frequency axis, find strict local maxima of normalized excess power Z using the same scipy find_peaks default as the earlier diagnostic. Keep maxima with height>=4. Enumerate all unordered distinct peak pairs with at least one interior frequency bin (index gap>=2) and an interior valley<=0.8*smaller peak. Pair score=min(two heights); do NOT add an anomaly,accepted-cluster or true-frequency filter. These are descriptive spectral structures, not detector claims. Save every retained peak/pair and the rejected-gate counts; one frame can contain many pairs.

Run two predeclared structural views: ALL separations and NEAR pairs with separation2..3Hz inclusive. Near is a morphological category scanned over ALL carrier positions, not a supplied true frequency or a calibrated search for two known tones. It enables like-for-like comparison of nearby peaks in controls; full-range findings remain visible. Do not choose a category after inspecting results.

Track exact-bin persistence: same sorted frequency pair in consecutive frames, no skipped frames. Track permissive drift persistence: directed links only between adjacent frames if both sorted component frequencies move at most1Hz; choose longest path ending at each pair, ties by smallest sum of absolute component steps then lexicographically smaller previous pair. The path may branch/stitch candidates; this is not physical component identity. Save best exact and best permissive paths, frequencies, per-step shifts and per-frequency span. Select longest path overall, then highest weaker-endpoint height, then lexicographic endpoint. No gaps bridged. A lost pair ends a track.

Also select the highest-score pair PER FRAME (ties by sorted frequencies) and report its jumps across adjacent frames, including large jumps. This answers whether blindly choosing strongest structures hops around; absence breaks the transition. Do not conflate it with the best retrospective longest path, which is selected after viewing the full sequence.

## Evaluation-only annotations and reporting

Only AFTER blind pair enumeration/tracking, annotate two-tone candidates with one-to-one matching to both true carriers within0.5Hz. For single-tone controls annotate pairs containing one peak within0.5Hz of their single true carrier; the other peak is potentially extraneous, not a recovered source. Noise has no carrier annotation. Keep whole-map blind persistence counts for all three groups.

Report full25frames and fixed central9frame centers1.25..2.25s for every group; the latter matches prior analysis and is applied identically to controls. It is an evaluation interval, not an input to blind full-map peak enumeration. For each case record frame counts, exact observed pair frequencies, longest exact/jitter run, longest carrier-associated run and strongest-pair jumps; attach the original24strict-success/5wrong/187no-pair group only after extraction.

Summarize run-length thresholds1,2,3,5descriptively for ALL and NEAR pairs, not as chosen classifier gates. Longest runs from highly overlapping windows are NOT independent confirmations. One-second windows overlap87.5%; a K-frame track has center span(K-1)*0.125s and union support1+(K-1)*0.125s. Frequency constancy of noise peaks can occur from shared samples. Small matched controls do not certify a false-claim probability. Twelve bases,shared variants,chosen amplitudes and known stationary generator limit inference.

For the5previous wrong pairs additionally report persistence of the exact frequencies that were previously returned (134/139 four times,86/89 once), and inspect the same pair on their noise and single-tone maps even when those controls have no accepted cluster. No frequency matching is used to create the general blind tracks.

## Preservation and verification

Fresh UTC/UUID output folder, archived protocol/script, parent hashes, all frame lists/track summaries,caseCSV,frequencyCSV,aggregate tables,plots and report. No maps or parent reports overwritten. Independently check pair gates,tracking recurrences,frequency-step constraints,frame/time arithmetic,counts and parent hash integrity. Retain null/negative outcomes. Figures describe selection rules and truth annotation. Any new persistence-based decision rule would need a separate protocol,independent calibration and new false-claim evaluation; none is adopted here.
