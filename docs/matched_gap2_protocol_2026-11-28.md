# Matched 2 Hz comparison — frozen exploratory protocol, 28 November 2026

## Motivation and prior evidence

BSD-20/21 and the gap-2 review:5/72 estimates within1Hz,4/72 within0.5Hz; all stricter successes equal-amplitude and on-grid. Noise, carrier and global phase differed between the old factorial cells, so those associations were not causal effects. Registry overlap O7 covers the question. New information here is matched interventions with the same noise/carrier/global phase, fresh independent calibrations and evaluation streams. Previous data are development evidence and not reused as blind targets. Production and all historical artifacts remain unchanged. No O3/data/manifest access, network, commit/push, main, release or Zenodo.

## Frozen methods

Reuse exactly the archived blind-pilot scales, extraction and D0/D1 proposal functions, unchanged and hash-checked against its saved driver. Four outputs W0D0,W0D1,W1D0,W1D1; primary descriptive contrast W1D1 versus W0D1. W0 longest0.5s, W1 longest1s,hop0.125s; same kernels21/7 (different frequency footprint remains a confound). D1 proposes peaks inside accepted clusters; D0 uses disjoint accepted boxes. No frequency, source-count, amplitude, onset or phase truth enters inference. No new tuning, interpolation or peak rule.

Each W gets its own new284noise calibration (85A/199B,FAP.005). Each W/D gets199newnoise and199newone-tone maximum-pair nulls; strict rank198 atalpha.01 for each group, maximum cutoff acrossgroups. Missing proposal/cutoff=-infinity, serialized null, same explicit rules as the previous pilot. Freeze threshold/hash before evaluation. One calibration draw only; no population guarantee from a small control sample. Different cutoffs across methods are permitted; decisions use each method's own cutoff.

## Inputs and fixed matrix

1024Hz,4s,pink noise sigma1,alpha1,knee16. Main gap2Hz, one-second Tukeyalpha.15 envelope, onset1.25s and resultant root-energySNR18. Strength/duration/onset are FIXED here; weak mixtures and onset effects are not addressed.

12independent base units. Each has one noise draw, one even carrier in80..178Hz, and one global phase. Cross exactly offsets{0,.5}Hz x relative phases{0,pi/2,pi} x weaker/stronger amplitudes{1,.5,.25}:18variants/base,216mainpairs. Same base noise/carrier/global phase in everyvariant; relative phase is defined at onset and evolves for unequal frequencies. Each summed addition is normalized to resultantSNR18. Save component energies/cross terms: phase/ratio changes can change the required component scale. This estimates effects under FIXED RESULTANT ENERGY, not fixed individual component amplitude. Grid intervention shifts BOTH carriers by.5Hz while retaining gap2,noise and phase; colored-noise level changes with frequency, so this is a physical carrier shift relative to the grid, not an isolated FFT-coordinate operation.

Each base also supplies the exact noise-only record and one matched single tone at each offset (12noise,24single controls). Variants sharing a base are not independent; report paired transitions and denominators in BASE UNITS. Never pool216variants as216independent noise draws or use a naive binomial test over them. Additional160independent noise and160single-tone controls form separate held-out false-claim groups, not pooled with matched controls. Single-tone control generator matches the previous broad SNR8/18/36 distribution.

Streams SeedSequence([seed,index]): development40001; anomaly calibration40101; pair-calibrationnoise40201 andsingle40301; evaluationnoise41101 andsingle41102; mainnoise41201 (12units),mainparameters41202. No unintended cross-split hashes. Main variant inputs paired deliberately. Twelve carrier/global-phase draws are small and cannot establish a universal phase/grid law.

## Outcomes and stage diagnosis

At most one frequency pair/method/input; duplicated estimates invalid, exact duplicate proposals consolidated. Primary correctness: one-to-one BOTH frequency errors<=0.5Hz (inclusive); retain historical loose<=1Hz and all wrong/no-pair outcomes. Wrong pair is a claim failing matching on two-tone input; every pair claim on single/noise is false. No claim is not a correct rejection of an observable tone.

Within every phase/ratio cell, compare on/off-grid perbase:bothcorrect,on-only,off-only,neither. Within everygrid/ratio compare phases0vs pi/2 and0vs pi; withingrid/phase compare ratio1vs.5 and1vs.25. Keep everytransition and count only12bases per comparison. These are exploratory contrasts, multiple tests; no confirmatory p-value or adoption criterion. Also compare W methods on EXACT same records, retain gains/losses. Count whether a failure is no anomaly, no accepted longest-scale candidate, no eligible pair or wrong pair; avoid equating these stages.

Timing: perpair onset/end box and midpoint errors as in oldprotocol; same support[1.25,2.25]. Conditional success subsets differ; report n. Runtime bank and proposal timers as before, plus archive I/O separately if measured. Fixed boundary position means no new boundary generalization or streaming-latency claim.

## Preservation and checks

Unique UTC/UUID folder; frozen protocol,driver,dependency snapshots,source hashes,all calibration maxima/candidates,null scores, input metadata/hashes, all inference records and evaluation, pairedtransitionCSV, plots and report. Evaluation longest-scale NPZ holds waveform/complexSTFT/rawpower/background/Z/mask/labels; large traces stay local/gitignored with hashes. JSON is gzip from the outset to avoid GitHub size limits. Full shorter-scale maps not retained.

Before evaluation verify old matching fixtures, exact paired noise/onset/carrier/phase preservation, energy identity and no unintended parameter changes. Afterwards independent saved-data checks reconstruct calibrations,nullcutoffs,claims,matching,allpairedtransitions,regenerate waveform hashes and compare selected direct production outputs. Record errors, exclusions and interruptions; never enlarge the sample or change methods after seeing results. Update the knowledge map with scope/negative findings; no automatic production adoption, source reconstruction or stem-separation claim.
