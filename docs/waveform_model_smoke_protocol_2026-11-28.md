# Waveform model fitting: eight known examples, implementation smoke study

28 November 2026. Owner approved trying the proposed one/two-tone waveform comparison on a few known examples. Consult BSD-23..26 and registry O7. New cue:fit phase-preserving raw samples rather than require separate spectral maxima. This is exploratory implementation/optimization/runtime evidence on ALREADY INSPECTED DEVELOPMENT INPUTS,not independent detection performance. No production BSD change, O3/data/manifest,commit/push,mainmerge,releaseorZenodo.

## Fixed sample and assumptions

Eight original matched W1trace waveforms: b00_two_0_0_1,true88/90equal; b00_two_0_0_0.5,true88/90weak; b00_two_0.5_0_1,true88.5/90.5off-grid; b05_two_0_0_0.25,true86/88previouslywrong86/89; b04_two_0.5_0_0.5,true134.5/136.5previouslywrong134/139; b09_one_0_0_1,true102only; b00_noise_0_0_1 and b09_noise_0_0_1. These labels choose smoke cases and score outputs ONLY;they are not optimizer inputs. Negative controls are two shared-base noise records and one single-tone record,not a calibrated false-rate set.

Use all4096rawsamples at1024Hz. The known experimental-design envelope is shared across models:one-second Tukeyalpha.15 onindices1280..2303 (nominal1.25..2.25s). This is a KNOWN-ENVELOPE/ONSET conditional task,not blind transient localization or unknown-envelope real-data inference. Do not fit/test exact support from each case’s truth. Frequencies20..450Hz,amplitudes and phases unknown. Fit0,1or2sinusoids multiplied by that common envelope plus an unrestricted shared constant offset.

## Independent background weighting

Regenerate284original anomaly-calibration noise waveforms (stream40101,unchanged archived blind generator),verify each hash against saved calibration_raw.json.gz. Estimate a common background spectrum by averaging perrecord Welch PSD(Hann1024,overlap512,constant detrend,spectrum density) on these independent background inputs. Interpolate onto4sFFTgrid,positivefloor1e-6*medianpositivePSD;inverse-square-root weights normalized by medianPSD. No target or same-base control data fit the noise weights.

Whiten raw target and EVERYdesign column with the SAME realFFT weights. Mean removal/constant offset handled identically. This is a circulant stationary-noise weighted least-squares approximation,not an exact likelihood or a validated real-sensor noise model. No foreground-trained whitening. A larger weighted fit gain is not automatically statistical significance.

## Frozen optimization recipe

At fixed frequencies solve sine/cosine coefficients jointly by numpy least squares (variable projection);do not successively subtract a fitted first tone and treat residual frequencies as a final answer. Optimize frequencies continuously;no requirement for two STFTpeaks. Precompute a0.5Hzsingle-frequency grid over20..450Hz. M1:refine up to5highest coarse local maxima with bounded scalar minimization in+/-0.5Hz,max80iterations,xatol1e-7;global best kept. Endpoints and global-grid best retained as candidates.

M2:starts from (a)bestM1frequency plus6best conditional coarse-residual joint-fit partners,(b)bestM1center plus symmetric separations{.25,.5,1,2,4,8,16,32,64}Hz,(c)allpairs of3bestcoarseM1maxima;deduplicate,startbound20..450,no exactduplicatefrequencies. At most18starts,deterministic ordering. Refine BOTHfrequencies jointly with bounded scipy least_squares,jac2point,max_nfev80,ftol/xtol1e-9,gtol1e-7. Frequency ordering is only for reporting. Coarse-residual scans are initializers,not final sequential estimates. Save every start,status,cost,evaluation count and returned frequencies,including failed/boundary/ill-conditioned solutions. Local multistart does not prove global optimality.

Accept no calibrated tone-count decisions in this smoke run. Report weightedQ0,Q1,Q2,explainedenergy and additionalgain Q1-Q2,condition numbers and optimizer convergence. A condition number>1e8is flagged as unreliable,not silently removed or interpreted as one tone. No BIC/chi-square cut or empirical significance threshold is improvised on these8examples.

## Checks, evaluation and stopping

Before sample fits,check noiseless1/2-tone design reconstruction at fixed artificial fixture frequencies unrelated to the smoke target carriers,offset recovery,permutation invariance and nested loss consistency. Same-sample2-tone fits should not be worse than the best included M1baseline;retain/report any failure. No frequency/onset/phase truth enters fitting. Evaluate2-tone sample frequencies AFTERfits:one-to-one sortederrors<=0.5Hz,also retain fullerrors. Single/noise M2frequencies are nuisance fits,not true2-tone discoveries or source reconstruction. Optimization/time failures count as failures,no case replacement.

Stop after these8cases;do not expand or retune based on the result. Report initialization/precomputation time,permodel/casewalltime,totalcalls and hardware/libraryversion,single-thread numericalenvironment. Save fresh UTC/UUIDoutputs,script/protocol,backgroundhashes/PSD,weightedlosses,allstarts,frequencies,linearcoefficients,raw/fitted/residualarrays andillustrative combined-waveformplots,source/parenthashes. Independent checks reconstruct objective/predictions and coefficient normal equations,frequency matching and background input hashes. This checks implementation,not the scientific validity of choosing a number of tones.

Next stage only after review:independent noise/single-tone threshold calibration and fresh held-out comparisons of tone count AND frequencyaccuracy,including envelope-mismatch controls and unknown onset. No threshold or production adoption here.
