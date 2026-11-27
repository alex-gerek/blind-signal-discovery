# Corrected cadence and warning-impact control — 27 November 2026

Pre-execution protocol, source cc68b93fee699307e3314fce1bd4c9fe403c1e6e. This SUPERSEDES only the inconsistent sparse-versus-repeated comparison in the earlier warning pilot; do not overwrite it or reinterpret its data as a corrected experiment. ONE unchanged default BackgroundRegimeMonitor and frozen MultiResolutionDetector, warning-only annotations. No new threshold, warning classifier, gating, adaptive fit, production/UI edits or automatic adoption.

## Fresh paired inputs and the corrected schedule

Five new seeds21101..21105. Same declared continuous-pink bench model:800*4096samples,1024Hz,sigma1,alpha1,knee16,SeedSequence([seed,1]),no per-window RMS renormalization. Source0..283, targets284..799, frozen284x6monitor reference and core85A/199B,FAP0.005. Existing RegimeMonitorConfig defaults unchanged. First11targetwindows are monitorwarmup; search availability remains unaffected. Cases/alphas/windows are not independent replicates; seed is replication unit.

REPEATED=list(range(288,800,8)) has64slots. SPARSE=REPEATED[2::5] has13slots. Family at ANY injection position is determined by its repeated-grid slot j=(index-288)//8:chirp if j even,burst if odd. Since sparse slots differ by5(odd), both schedules alternate families AND share identical families at all common positions. CleanRNG SeedSequence([seed,3,index,family_key0chirp/1burst]) is shared acrosscases. Same noise+family+index implies byte-identical sparse/repeated observation, amplitude and known support. Assert subset,alternation,equalinputhashes and coreoutputs before interpreting a cadence effect.

Eight cases, each on the same base trace perseed:

1.stable:base,noaddedsource.
2.level_step:2*basefrom360.
3.spectral_line:base+3*sin(2*pi*110*t+phase)from360,globaltime,phase SeedSequence([seed,2]);noise-label case.
4.sparse_transients:stablebaseplusSPARSEsignals.
5.repeated_transients:stablebaseplusREPEATEDsignals.
6.persistent_signal:EXACT spectral_linewaveform,buttoneisdesiredsignal,notnoise-lawchange. Evaluationlabels only; monitor cannot identify physicalorigin.
7.step_with_probes:level_stepnoiseplusSPARSEsignals.
8.line_with_probes:spectral_linenoiseplusSPARSEsignals.

Transients root-energySNR18 relative to std(the actual profiled noise window),linear notdB. Existingdefaultchirp and fixedburstGaussian4096,Butterworth4zero-phase35..350Hz,mask1.375..1.625s,Tukey0.15. The same clean waveform/support is used across sparse/step-probe/line-probe arms at the same index, but different noise regimes change amplitude scaling. Paired noise-only cases are retained BEFORE injection. Truth/support never enters monitor, calibration or admission; it is only for evaluation.

## Fixed questions and endpoints (descriptive, no superiority test)

A. Corrected cadence: compare warning flags/history on the13shared identical probe windows perseed; core scores,candidates and detection/association decisions MUST be identical there. Count warning disagreements and hypothetical lost recoveries separately, not a change in core sensitivity. Also retain all additional repeated probes and all stable false-warning episodes. A null difference is evidence to retain, not a reason to increase dose or retune.

B. Warning impact on the search statistic: pair stable versuslevel_step and sparse versusstep_with_probes at all516targets; record maxscore difference,empiricalp/nomination disagreement,noise alarms and signal associations. A raw gain warning neednot imply detector impairment. Compare spectral_line versusstable noise alarms with actual denominators; no assumption that a matched nominalFAP gives equal achievedFAR. Do not invent an automatic harmless/dangerous classifier from these outcomes.

C. Injection-specific evidence: compare each SPARSEpositive observation with its exact pre-injection noise counterpart (sparse->stable,step_probes->level_step,line_probes->spectral_line). On BOTH candidate lists evaluate the SAME known support rectangle after detection. Count positive support association,already-present noise support association,new association,lost association,positive-only whole-window nomination and score differences. An association already present in noise is not proof the injected signal was identified. This is an oracle diagnostic of a synthetic benchmark, not a deployable decision rule.

D. Warning consequences: preserve all actual search/candidate outputs and available exposure. Report warnings,contiguous episodes,startup,firstpostchangeflag versus NEWepisode,logicaldelay from changeSTART to assessmentEND,signal-confounding and hypothetical warning-gate losses. For noise-onlycases stable/level/line, evaluate fixed future windows onset+1..onset+12 aftereachwarningepisode; merge overlappingfuture intervals and truncate at799. Report actual future noise denominators/alarms. This retrospective descriptive association is not causal prediction validation; do not include the warning window as a future outcome.

Persistentline/desiredtone raw/core/monitor identity remains a mandatory ambiguity control. Warning is not a noise-origin label or permission to relearn a persistent signal. No recalibration source is selected here; L09/L10 remain open.

## Numerical checks and preservation

Reuse the earlier independent NumPy-Hann feature oracle,exact integer-lattice KS path check,CUSUM/localZ/warmup reconstruction,core A/B/rank/p/candidate checks and prefix comparisons at12/76/96. Regenerate all32000percase input records; manualfeaturechecks first3reference perseed and first3target percase (135checks);120productionprefixchecks. Frozen reference/core/source/dependency hashes and unmodified sidecar search hashes required. Recheck all candidate bounds,post-detection supports and real exposure. Core availability must not change when a monitor warns.

Additional schedule assertions before measurement and full savedinput/core identity on sharedpositions are required. Independent impact verification must reconstruct cadence discrepancies,gain score deltas,pre-injection support comparisons and unionedfuture exposure without invoking the production monitor or using its summary helper. Exact-gain invariance is an observation to verify,not a source of parameter selection; retain any disagreement.

No inferential p-value for efficacy,pooledIIDCI,validated monitorbudget,thresholdsearch or defaultreplacement from5pairedseedunits. This is a fixed small control with two selected noise changes andSNR18families. Freshdata may differ from the prior pilot; do not merge counts or claim an improvement caused by the protocol correction alone.

UniqueUTC/UUIDfolder;copiedpre-runprotocol/scripts/source/environment/dirtyhashes,fullfeatures/core/monitor/sidecar outputs,pairing/futureexposure diagnostics,errors and independentverification. Preserve every seed/case,including failed expectations. No old artifactoverwrite,O3/data/manifest/network/realwaveforms,production/testedits,commit/push,main/release/Zenodo. Initialstatus:protocolfixedbeforemeasurement;noautomaticrulechange.
