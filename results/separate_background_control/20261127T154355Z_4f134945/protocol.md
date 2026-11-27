# Separately acquired background: fixed-source control — 27 November 2026

Pre-execution protocol, source d049e2b. Unchanged MultiResolutionDetector; no tail guard, trimming, new adaptive policy or production changes. This is an exploratory synthetic source-control pilot, not operational certification that a real background can be acquired without signals. No O3/data/manifest/network or real waveforms.

## Question and four frozen arms

Does replacing signal-contaminated accumulated background with a separate, signal-free draw from the SAME noise process preserve sensitivity without abstention? Distinguish distribution matching from mere sample separation.

Each arm has265calibration observations,66normalizationA/199thresholdB, norm_fraction0.25; FAP0.01 and0.005, production rank/strict p rule unchanged. Calibrate all arms before evaluating any target. Store calibrations and freeze hashes; never match observed alarm counts afterward.

1. stream_contaminated: accumulated calibration noise with four known added signals at100/140/180/220, all inB. Chirps at100/180, bursts at140/220, root-energy SNR18 relative to EACH underlying noise window's standard deviation.
2. paired_noise_oracle: exactly the same underlying calibration noise, before those four additions. This is an inaccessible paired diagnostic, not an implementable source.
3. separate_matched:265new signal-free noise observations from an independent stream, same noise distribution/profile as the target. No scores/target data or truth labels used for choosing samples; retain every noise extreme.
4. separate_old:265independent stationary pink observations without the changed profile. In the stationary case this is another matched draw; in other cases it is intentionally mismatched. Separation alone does not establish calibration relevance. Never pool these arms as independent replications.

Oracle differs from stream at exactly four samples and shares allA normalization. Separate matched also replacesA, so it need not reproduce oracle thresholds or improve every signal. Save bothA/B effects; do not attribute differences automatically to contamination removal.

## Fixed independent streams and data

Seeds18101/18102/18103; cases indexed0..3:stationary,level_step,spectral_line,impulsive_noise. For each seed/case, RNG SeedSequence([seed,case_index,stream_id]), with stream1accumulatedcalnoise,2matchedcalnoise,3oldcalnoise,4targetnoise,5calibrationcleanwaveforms,6targetcleanwaveforms. Three independent seed replications percase; calibrations/arms/alphas/amplitudes share observations and are not extra independent units. No existing pilot inputs reused.

1024Hz,4s,4096samples. Base noise.generate_pink_colored sigma1 alpha1 knee16. Stationary unchanged. Level_step:allcalibrationandtargetobservations current regime2*base; oldsource remainsbase. This is a known post-step distribution, not a continuous transition/recovery experiment. Spectral_line:base+3*sin(2*pi*110*t+phase),phaseuniform0..2piindependentlyperrecord; oldsource lacksline. Impulsive_noise:withprobability0.03addone signedamplitude12impulse atuniformindex; oldsource lacksimpulses. No profile renormalization and no selection/filtering of noisy observations. Impulses/line are noise for accounting; that semantic choice is explicit.

Calibrating signals:existing linear_chirp generator default distributions; burst independent4096Gaussian filtered Butterworth4zero-phase35..350Hz,mask1.375..1.625s,Tukey0.15. Evaluation200newnoise observations, plus10independent target noise/clean units perfamily, each atSNR8/18 (40signalrecords,20baseunits/case). Targets identical acrossallfourarms; amplitudes share noise/clean withinunit. Root-energy SNR scales by std(profilednoise),linear notdB/matchedfilterSNR. Save generator support rectangles; count both whole-window nominations and accepted-candidate support association. The latter is a loose known-support criterion, not waveform identification.

All background sources and target noise bases use independent streams; verify no waveform-hash overlaps between accumulated/matched/old/target sources. Pairednoiseoracle intentionally shares accumulated bases; paired amplitudes intentionally share target bases. Hash disjointness detects reuse but is not a proof of distributional independence.

## Measurements, interpretation and falsification

Everycalibration:source provenance,orderedmembers/hashes,A/Bpartition,rawmaxima,fullnormalization/nullscores/threshold/rank and FAP. Before target analysis freeze calibration hashes and record zero target decisions. Each target:raw maxima,score,p,threshold,acceptedcandidates,availability,signal/supportlabels only for evaluation. No admission/quarantine:all attainable calibrations remain available; if not, expose actualnull/abstention and retain failures.

Report perseed/case/arm/alpha:noisealarms/200,available/noise seconds,signal nomination and supportassociation byfamily/SNR,availablemisses/abstentions,and paired gains/losses. No artificial equalization of thresholds, observednoiseFAR or conditional denominators. Frozen-common nominalFAP is not equality of true achieved false-alarm rates. No superiority/adoption p-values or pooled binomial CI over differentcalibrations; descriptive results only. With199B andthreecalibrationdraws, rates can vary substantially; any finalbudget assurance needs a new properly sized independent study.

A matched source failing to improve a group is retained. Mismatch may or maynot degrade measured behavior; donotassumeit must. Even positive results demonstrate only what source separation under known matching can achieve, not an autonomous signal-free admission rule. L09/L10 not closed:adaptive timing,background acquisition purity,continuous drift,Union and operational release/latency are untested.

## Independent verification and preservation

Reconstruct allAmedians/Q95spreads,Bmaximum scores/ranks,thresholds,target scores,empiricalp-values bydirecttailcount,strictnomination and everycandidate/supportassociation fromsavedmeasurements. Check exactfour-sample oracle difference,input/source separation,freezehashes and actual exposure. Regenerate allinputhashes withsavedseeds; directuncachedrawparity first3calibrationobservations perarm and first3targetobservations percase. These checks do not prove scientificgenerators correct.

UniqueUTC/UUID output,protocol/driver/helpercopies,source/dependency/environment/dirtyhashrecords,fulloutputs/errors,independentverification and final report. Neveroverwriteoldruns. Preserve failures and inconclusive outcomes. No tests rewritten, productionintegration, commit/push, main/release/Zenodo. Initialstatus:protocolfixed,executionpending.
