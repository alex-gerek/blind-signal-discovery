# Two known oscillations, template-free anomaly detection — 27 November 2026

Pre-execution protocol, source cc68b93fee699307e3314fce1bd4c9fe403c1e6e. Scope: anomaly detection of mixtures, not source separation, inverse waveform reconstruction, novel physical discovery or an information-theoretic limit for all methods. The unchanged MultiResolutionDetector receives ONLY mixture observations and an independent background. Component identities, carriers, labels and timing never enter fitting or detection; they are retained for evaluation only. STFT/clustering/background assumptions remain explicit even without an exact waveform template.

## Fixed data and calibration

Five new seeds22101..22105. Perseed three independent RNG streams SeedSequence([seed,1])calibration284pink windows,[seed,2]held-out200noise windows,[seed,3,unit]20paired mixture units.1024Hz,4s,4096samples,noise.generate_pink_colored sigma1 alpha1 knee16. Calibration norm_fraction0.30,85A/199B,targetFAP0.005. Freeze calibrations/hashes before target detection. No score/threshold tuning,veto,adaptive,monitor,Union or real data. Windows generated independently here, unlike prior continuous-sequence studies; amplitude variants/cases share noise perunit.

Eachunit draws one new pinknoise vector,integer carrier f uniform96..176 inclusive,three phasesuniform0..2pi. Identical Tukey0.15envelope supported1.25..2.25s. A=sin(2pi*f*t+phaseA)enveloped then scaled root-energySNR18relative to sigma1. Bfarcarrierf+96,Bclosecarrierf+2,independentphases,each separately scaledSNR18. Scale before summing; NEVER renormalize the mixture. Save all component/mixture energies, dotproducts, actual residualSNR and inputhashes. Any energy-dependent detection difference is not automatically a frequency-resolution effect.

Nine cases perunit: Aalone,Bfaralone,Bclosealone,AplusBfar,AplusBclose,AplusA,exact_cancel(A-A),near_cancel(A-0.95A),equivalent_single(2A). Also analyze the exact noise alone. Components cancel by vector negation, not approximate sin(phi+pi). AplusA and equivalent_single are byte-identical observations under different source descriptions; exact_cancel equals noise alone. These identities must hold in inputs and fullcoreoutputs (exceptsampleID), not merely detection status.

The cancellation case is NO_OBSERVABLE_ADDITION, not a detector miss or abstention: nonzero latentcomponents produce zero net addition. Nearcancel has root-energySNR0.9; inphase has36. Record energy/cross-term identity ||A+B||²=||A||²+||B||²+2<A,B>. Do not present different resultant energies as matchedSNR comparisons or pool pairedcases as independent injections. No progressive knowledge-removal algorithm is implemented: the anomaly detector has no component knowledge from the start.

## Evaluation

Count full-window nominations and candidate intersections with evaluation-only time interval1.25..2.25s and carrier diagnostic bandsf±4Hz,far(f+96)±4,close(f+2)±4. These are diagnostic bands, not strict physical spectral support of a finite tapered tone. One candidate may cover multiple bands; that is NOT evidence that two sources were recovered.

For pair morphology record whether TWO accepted candidates within the SAME STFTscale have disjoint frequency boxes,eachintersects the shared time interval and respective carrierband. Distinct IDs or different scales alone do not establish two resolved structures. Save per-scale outcomes and candidate counts. Even disjoint structures are not recovered waveforms or known physical sources. Evaluate the same bands on pairednoise-only outputs for accidental associations.

Report9cases separately (100units/case),100pairednoise controls and1000independent held-outnoise windows. No aggregateefficiency including exactcancellation; expose its zeroresidual status. Report actual rootenergySNR,per-scale maximal scores,thresholds,p-values,nomination,candidates,association and morphology. Limitedfrequency/duration/SNR choices do not establish arbitraryunknown-signal sensitivity. Passing mathematical tests does not prove scientificcorrectness.

## Checks and preservation

Beforemeasurement fixedarrayfixtures verify coherent/equivalent/cancel identities and energy sums. Independently reconstruct A/Bcalibration,rawscorestandardization,rank,p-values,allnomination/candidateband/morphologydecisions and pairednoise counts from savedrecords. Regenerate allwaveforms/hashes; assert no cross-stream noise reuse exceptexplicitpairing; compare directuncachedrawoutputs on first3calibration,first3held-outnoise andfirstunit all9cases perseed. Identicalinput/coreoutput checks mandatory; no post-hoc search for a resolvable frequencygap or retuning to improve a score.

UniqueUTC/UUID outputs,exactprotocol/scripts/dependency/source/environment/dirtyrecords,fullrawmaxima/calibrations/measurements/parameters/energies/errors,independentchecks andfigures. Preserve priorresearch and known failures. No O3/data/manifest/network/realwaveformaccess,production/app/test changes,commit/push,merge/release/Zenodo. Initialstatus:fixed before measurement.
