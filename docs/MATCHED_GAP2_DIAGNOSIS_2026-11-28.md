# Matched gap-2Hz decisions: successes, wrong peaks and no-pair gates

28 November 2026. Saved-data analysis only, unchanged methods and thresholds. Parent run20261128T115035Z_bf957291;216variants share12base units. [Full report and figures](../results/matched_gap2_diagnosis/20261128T130736Z_b2471926/report.md), [all24correct cases](../results/matched_gap2_diagnosis/20261128T130736Z_b2471926/correct_pairs.csv), [five wrong pairs](../results/matched_gap2_diagnosis/20261128T130736Z_b2471926/wrong_pairs.csv), [187no-pair cases](../results/matched_gap2_diagnosis/20261128T130736Z_b2471926/no_pair_cases.csv).

## Correct pairs

24strict successes:12equal-amplitude/phase0/on-grid,4equal/phase90/on-grid,1equal/phase180/on-grid,3equal/phase0/off-grid,4ratio0.5/phase0/off-grid.17pairs exactly match both carriers;7have0.5Hzerror on both. Peak heights6.6984..38.5241in mean normalized excess-power units (gate4). Valley depth20.94..50.70%below the smaller peak; required20%. Some are only slightly beyond the dip criterion, not unambiguously distinct physical sources. The CSV retains every true/estimated frequency, peak prominence,height,depth,margin and condition.

## Wrong pairs

Four correlated variants on base04have true134.5/136.5Hz but choose134/139Hz (errors0.5/2.5Hz), with ratio0.5or0.25 and phase0or90. The repeated139Hzpeak already exists as a local maximum in the same-column noise control. Raw-power accounting confirms noise plus signal leakage plus their interference, not a separate injected139Hztone. It is nearby extraneous morphology, not pure rounding and not a distant-band jump.

For base04/ratio0.5/phase0 ONLY, a nearer137Hzmaximum exists but134/137fails the valley gate (ratio0.911), whereas134/139passes (0.546). In the other three errors,137Hzis not a qualifying local peak. The single-tone accepted cluster ends136.5Hz, excluding139Hz; adding the second tone expands the box to140.5Hz in this example, admitting the extra maximum. This is why the single/noise controls need not claim a pair even though their spectra contain that feature.

Fifth wrong case:86/88Hz ->86/89Hz, weaker/stronger amplitude0.25,phase0,on-grid. The1Hzerror passes the old loose rule but fails strict0.5Hz. All5wrong weak peaks are barely above the height gate (4.2586..5.1889). Raising the gate after looking at this separation would be tuning on development data, not a validated fix. Wrong peaks can have deeper valleys than correct pairs.

## No pair: exclusive stopping-point accounting

| First failed mean-spectrum gate | Number |
| --- | ---: |
| Fewer than two local maxima within an accepted cluster | 141 |
| At least two maxima, fewer than two at height>=4 | 8 |
| Two sufficiently high/spaced maxima, valley too shallow | 38 |
| Interior-bin spacing / later calibration or pair cutoff | 0 |

All187have an accepted long-scale anomaly. These are proposal-gate failures, not unavailable calibration or missing energy. 'No second maximum' describes the averaged representation, not proof of absent waveform information.

Important follow-up:100of the141single-peak averages have an individual saved frame with two carrier-associated peaks that satisfy the same height/spacing/valley gates; in the mean those peaks no longer separate. Across all187no-pair cases,140have at least one such frame. This identifies temporal averaging as a stage that can remove peak morphology in these records. Frame annotations used true frequencies AFTER inference; scanning many times can find noise peaks. No new blind frame-scanning method, false-rate calibration or recovered140/187performance is established.

## Checks and decision

216profile reconstructions and all29claim/no-claim decisions match the stored blind diagnostic. Gate counts sum187. Ten saved raw-power decomposition identities checked;529parent artifact hashes unchanged. Both figures visually inspected. No new waveforms, altered gates, exclusions, scientific rerun, O3, production changes, commit or push.

Next question: can a truth-blind diagnostic requiring temporally consistent peaks retain useful evidence while controlling one-tone/noise false claims? It requires a new fixed protocol and independent calibration/evaluation. Do not simply maximize over frames or raise a threshold on these examples and claim an improvement. Preserve the distinction between anomaly nomination, frequency estimation and waveform/source reconstruction.
