# Frozen warning-only regime monitor: independent sequence pilot — 27 November 2026

Pre-execution protocol. Source6f4ccfac7444fb01c2b19785e4236536559ce7dd. ONE existing BackgroundRegimeMonitor in causal_trailing mode with unchanged RegimeMonitorConfig defaults: local12,KSalpha0.005,driftZ3.5,CUSUMk0.75,h12,internal offset1.8,STFT128,existingfourbands. No tuning, extra confirmation filter, restart after alarm, recalibration, candidate suppression, detector/app/test changes or automatic adoption. The monitor produces generic distribution-change WARNINGS only, never noise-origin labels.

## Frozen reference, search detector and information limits

Five new seeds20101..20105. Each seed generates one continuous800*4096sample pink trace at1024Hz,sigma1,alpha1,knee16,SeedSequence([seed,1]); split into800four-second observations without per-window RMS renormalization. This is a jointly FFT-shaped periodic-boundary stochastic model, notIIDwindows. Reference0..283 is shared across six paired cases; freeze its284x6feature matrix, medians/scales and hash before target monitoring. Production feature equations are logMAD amplitude,absQ99/MADtail ratio and median logbandPSD20..80/80..160/160..280/280..450Hz.

Fit an unchanged MultiResolutionDetector on exactly the same initial284samples, norm_fraction0.30,85A/199B,FAP0.005. Freeze this calibration; use it unchanged for every target284..799 and every case. This separate search FAP has no equivalence to monitorKSalpha despite the same numeric0.005. No veto/adaptive/Union path is evaluated.

Warnings are annotations attached to untouched search measurements. Save candidate lists and actual core thresholds/p-values before annotation and prove their serialized content is unchanged afterward, including when warnings occur. Do not gate detection, remove exposure or train on target data. Monitorwarmup is not search unavailability.

## Six paired cases and fixed signals

1.stable:basepink,zeroaddedsource,no truechange.
2.level_step:noise2*basefromwindow360through799,noaddedtestsignals.
3.spectral_line:base+3*sin(2*pi*110*t+phase)from360;globaltime t,phaseoncefromSeedSequence([seed,2]);labelled NOISE-law change.
4.sparse_transients:stablebasepluschirp/burstalternating at[304,324,344,364,384,704,724,744,764,784],nochange innoise law.
5.repeated_transients:stablebaseplusalternatingchirp/burst atlist(range(288,800,8)),nochange innoise law.
6.persistent_signal:EXACTLYsameobservationsasspectral_line,butthetoneislabelledanintendedpersistentSIGNAL,notnoisechange. Truthlabelsneverentermonitor. This constructed ambiguity demonstrates that identical measurements cannot reveal different causal labels; it is not a real-event novelty or population-identifiability claim.

Transientroot-energySNR18 relative to std(originalnoisewindow),linear notdB/matchedfilterSNR. Existingdefaultlinear_chirp and fixedburstGaussian4096/Butterworth4zero-phase35..350Hz,Tukey0.15support1.375..1.625s. CleanRNG SeedSequence([seed,3,index,family_key0chirp/1burst]); sameinjectioninput at sharedindices across transient cases. Persistenttone amplitude3,support[0,4]s/110Hz,notSNR18. All candidates/support associations are diagnostics, not physical identification.

## Warning accounting and fixed horizon

Monitorprocesses516targets,with11warmupwindows(global284..294); first eligible full12-window assessment is295,available after its end1184s. Each later assessment uses at most the last12targetwindows and frozen reference. Save local/global indices and actual prefix. Existing min_latency_seconds4 describes an observation duration, not the48s first-full-history startup. No reset or reference replacement.

Report perseed/case: eligible windows,warningwindows,contiguous warning episodes/onsets/durations,anywarning indicator over this fixed horizon,trigger channels(KS/medianZ/CUSUM),dominant features. Stable warnings are false warnings relative to the known stable noise law. Do not call every flaggedwindow a new alarm. Per-feature Bonferroni is not a repeated-time false-warning guarantee; no naive pooled-binomial CI or claim of0.5%monitorFAP.

For truechanges report warnings before360,whetherwarningalreadyactive359,firstpost-changeflag and firstNEWepisodeafter360 separately. Delayseconds uses assessmentwindowEND minus changeSTART; same-windownotification has4sdelay,notzero. An alreadyactivefalseepisode cannot be credited as new detection. Save absence of a newepisode asnull,notzero latency.

Fortransients report warning overlap with signal windows,monitor-warmup signals,extra/suppressed warningwindows relative to pairedstableinput,postsignal warning persistence,and how many existing support-associated search recoveries WOULD be lost if warnings were used as a gate. This is an explicitly counterfactual tally; actual detector decisions/exposure remainunchanged. Signal-associated warnings are not themselves proof of a noise change or false scientific signal.

Persistent-signal and spectral-line warnings/core outputs must be identical despite different evaluation labels. Anywarning there must remain an ambiguity/credibility annotation, not permission to absorb a persistent signal as background. Persistent warnings are not scored as proof that the monitor is wrong about changed input statistics.

Strict exploratory desiderata, fixed beforemeasurement:zero warningepisodes in all5stable sequences; every level/line change starts a NEW warningepisode by371 (enddelay<=48s) with no activeprechange warning; zero additional warningwindows induced by sparse/repeatedtransients relative to stable. Report each separately. Failure rejects an unqualified deployable interpretation; even allpassing is only a small fixed-pilot screen. No significance test, superiority inference or revision of these gates after outcomes. Persistent ambiguity prohibits automatic noise classification regardless of other gates.

## Independent checks and preservation

Reconstruct reference medians/MADscales/fallbacks and allsavedfeature-derived z,CUSUM,localmedians,KSstatistics/minp,warningflags,warmup and episodes. Check KSstatistic by direct empiricalCDF and exact finite-sample two-sided p via independent integer-lattice path counting (12vs284) with caching of possibleDvalues; this validates numerical arithmetic under the KSmodel, not its applicability to dependent windows. If libraryfloatp differs nearzero retain tolerance and separately verify trigger classification; do not silently change productionKS.

Manualfeature oracle uses directHann128/hop64 NumPyFFT(noSciPYSTFTcall),medianDC,quantiles andbandmedians for first3reference and first3target windows perseed/case. Save/compare fullcore raw maxima/calibration/p/acceptedcandidate decisions; annotatedoutput must retain exactcandidate arrays. Productionmonitorprefixchecks at lengths12/76/96 againstfullassessmentprefix prove nofuturedependence for those prefixes; pairedstable/line prechange and full line/persistent identity are mandatory. Regenerate all24000input records,verifyreference hashes unchanged; record5seedsasreplicationunits,not30independentcases.

UniqueUTC/UUID outputs,exactprotocol/driver/dependency/source/environment/dirtyhashrecords,reference/corecalibrationfreeze,fullfeature/core/monitor/annotation/timeline outputs,errors,independentverification andplots/report. No originalresultoverwrite,O3/data/manifest/network/realwaveforms,productionedits,commit/push,main/release/Zenodo. Initialstatus:fixedprotocolbeforemeasurement;noadoptionorrealmonitoringbudgetcertification.
