"""Summarize verified warning-only monitoring with fixed-horizon and signal ambiguity limits."""
from pathlib import Path
import csv,json,sys
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
v=load(out/'independent_verification.json');assert v['status']=='PASS'
cases=['stable','level_step','spectral_line','sparse_transients','repeated_transients','persistent_signal'];rows=[]
for case in cases:
 entries=[x for x in v['case_summaries'] if x['case']==case];s=[x['summary'] for x in entries]
 rows.append({'case':case,'sequences':len(entries),'eligible_windows':sum(x['monitor_eligible_windows'] for x in s),'warning_windows':sum(x['warning_windows'] for x in s),
  'warning_episodes':sum(len(x['warning_episodes']) for x in s),'sequences_with_warning':sum(x['any_warning'] for x in s),
  'signal_windows':sum(x['signal_windows'] for x in s),'signal_warning_windows':sum(x['signal_warning_windows'] for x in s),'signal_warmup_windows':sum(x['signal_warmup_windows'] for x in s),
  'extra_windows_vs_stable':sum(x['extra_warning_windows_vs_stable'] for x in s),'suppressed_windows_vs_stable':sum(x['suppressed_warning_windows_vs_stable'] for x in s),
  'associated_search_recoveries':sum(x['associated_search_signal_windows'] for x in s),'hypothetical_warning_gate_losses':sum(x['associated_recoveries_that_warning_gate_would_remove'] for x in s),
  'search_noise_alarms':sum(x['search_noise_alarms'] for x in s),'search_noise_windows':sum(x['search_noise_windows'] for x in s),
  'detector_available_windows':sum(x['detector_available_windows'] for x in s)})
with (out/'summary_counts.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'aggregate_summary.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
lines=['# Warning-only regime monitor: frozen-default independent sequence pilot','',f"Run `{out.name}`, source6f4ccfa. [Pre-run protocol](protocol.md), [independent numerical checks](independent_verification.json), [counts](summary_counts.csv).",'',
'## What was tested','',
'ONE existing BackgroundRegimeMonitor with its unchanged default causal configuration:12window local history,KSalpha0.005,medianZ3.5,CUSUMk0.75,h12 and internal offset1.8; no tuning, debounce, reference resets or added classifier. Five fresh seed traces,800four-second observations each, six paired cases:stable,levelstep,line,sparse/repeatedtransients,and a persistent desired signal exactly identical to the line case. Five seeds are the independent sequence units; cases and windows are paired/correlated.',
'', 'Initial284observations freeze the reference and a separate core search calibration (85A/199B,FAP0.005). Targets284..799 provide516decisions. Monitorwarmup284..294 contains11unassessedwindows; firstfull12-windowresult295 is available at its end1184s,48saftertargetstart. The API4slatency field is an observation duration, not a claim of4sfull-history startup.',
'', 'Core search outputs and candidate arrays remain unchanged. The warning sidecar stores exact search-record hashes, never suppresses a candidate or gates exposure and never updates background. The hypothetical losses below show what would happen under a gate, not an action performed here. SearchFAP and monitorKSalpha have different meanings despite the same numeric0.005.',
'', '## Descriptive fixed-horizon accounting','',
'| Case | Sequences with any warning | Warning windows / eligible | Warning episodes | Extra / suppressed warning windows vs stable | Hypothetical lost search associations |',
'| --- | ---: | ---: | ---: | ---: | ---: |']
for r in rows:
 lines.append(f"| {r['case']} | {r['sequences_with_warning']}/5 | {r['warning_windows']}/{r['eligible_windows']} | {r['warning_episodes']} | {r['extra_windows_vs_stable']} / {r['suppressed_windows_vs_stable']} | {r['hypothetical_warning_gate_losses']}/{r['associated_search_recoveries']} |")
lines += ['', 'Warning windows in a continuous episode are NOT separate new alarms. The stable count is false warning relative to this unchanged synthetic noise law. Signal-associated warnings indicate altered measurements, not proof of instrumental noise or a false physical signal. Persistent-signal window counts represent one sustained episode per seed, not hundreds of independent signals.', '',
'### True-change latency','',
'| Case | Seed | Warning active before change | First post-change flag delay(s) | First NEW episode delay(s) |',
'| --- | ---: | --- | ---: | ---: |']
for x in v['case_summaries']:
 if x['case'] in ['level_step','spectral_line']:
  s=x['summary'];lines.append(f"| {x['case']} | {x['seed']} | {s['warning_active_before_change']} | {s['first_post_change_warning_delay_s']} | {s['first_new_episode_delay_s']} |")
lines += ['', 'A same-window warning has4snotification delay measured from changeSTART to assessmentEND. Already active false warnings are not credited as new detection. Null means no observed qualifying episode, not zero latency.', '',
'## Fixed exploratory gates','']
for name,value in v['gates'].items():lines.append(f'- `{name}`: {value}.')
lines += ['', 'These strict pilot desiderata were fixed before outcomes. Even allpassing would not validate long-horizon false-warning control or deployment; failures are retained, not used to retune defaults on this sample. No formal superiority or pooled-binomial inference was performed.', '',
'## Signal and noise ambiguity','',
'The spectral-line and persistent-signal cases have exactly identical inputs, reference, core scores/candidates and warnings. Only evaluation labels differ; support-association fields can differ because those use truth outside the algorithm. The monitor cannot determine which causal label is correct from these identical observations. This is an intentionally constructed ambiguity, not discovery of a new impossibility theorem or a real-event explanation.',
'', 'A generic warning can still be useful: a frozen reference may not describe the current measurements. It cannot authorize treating a persistent feature as noise, deleting it or fitting a new background to it. High apparent search recovery on a persistent line must be read beside its noise-case alarms. No claim that line warnings identify instrumental origin, and no GW170814/H1 causal inference.',
'', 'Sparse/repeated-transient cases share the same underlying stable noise. Extra warning windows and hypothetical lost recovered injections expose confounding, including warnings on nearby noise windows through local history/CUSUM. This pilot uses only two families atSNR18; absence of a conflict in a seed does not certify all signals.',
'', '## Numerical verification and preservation','',
f"Checks pass{v['reference_models']}frozen reference models,{v['detector_calibrations']}search calibrations,{v['monitor_records']}monitor records,{v['search_records']}search records,{v['candidate_records']}candidate records and{v['sidecar_records']}sidecar hashes. Regenerated{v['regenerated_inputs']}input-window records;{v['manual_numpy_feature_checks']}manualHann/NumPyFFTfeature checks and{v['production_prefix_checks']}productionprefix comparisons. The KSoracle independently checks{v['ks_feature_checks']}local/referencefeature comparisons via exactintegerpath counting; maxBonferronipdifference{v['max_ks_p_absolute_error']:.3g}. This checks numerical arithmetic under the KSmodel, not its statistical validity for dependent windows.",
'', 'Reference medians/MADfallbacks,localZ,CUSUM,KSstatistic/p,allflags/episodes,signaloverlap,corecalibration/empiricalp/strictacceptance,candidate standardization/bounds/supportassociation and unchanged warning annotation verified. Prefix lengths12/76/96 and pairedcase identities provide scoped causality/label-blindness evidence, not a formal proof for every possible input.',
'', 'No production/app/test changes, background update, candidate suppression, excluded case, measurement interruption, retuning, O3/data/manifest/network or real-waveform access. Source/dependency/reference hashes unchanged. The pre-execution verifier syntax check caught a keyword/dictionary typo; retained in verification_preparation_error.json and corrected before verification without restarting measurements. These counts are scoped research checks, not full-suite pytest passes or scientific correctness proof.',
'', '## Limits and decision','',
'This fixed34.4min target horizon has only five seed replications and two chosen noise changes. The existing per-feature Bonferroni correction does not control repeated-time warnings over an arbitrary deployment. JointFFTcolored noise, finite reference variation and limited signal families remain. No automatic reference replacement, expiry policy, validated monitorFAP, Union repair or real source-purity validation was tested. L09/L10 remain open.',
'', 'Keep this monitor as an explicitly provisional diagnostic sidecar. Do not adopt its warnings as noise-origin labels or automatic gates/recalibration triggers. Any redesigned rule needs a newly frozen protocol, independent sequences, a declared warning/latency budget and retained persistent-signal ambiguity. Do not optimize the thresholds on this failed/passing pilot until a desirable score appears.',
'', '![Stable and transient warning counts](warning_counts.png)', '', 'All five seeds shown; these are flagged-window counts, not independent alarms or FAP estimates.', '',
'![Identical input and warning under different labels](signal_noise_ambiguity.png)', '', 'Representative firstseed20101; full five-seed waveform/core/monitor identity checked. Labels are evaluation-only.']
(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Saved complete descriptive warning-monitor report and counts.')
