"""Report verified corrected-schedule and paired detector-impact results without inference tuning."""
from pathlib import Path
import csv,json,sys
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
v=load(out/'independent_verification.json');impact=load(out/'independent_impact_analysis.json');assert v['status']==impact['status']=='PASS'
cases=['stable','level_step','spectral_line','sparse_transients','repeated_transients','persistent_signal','step_with_probes','line_with_probes'];rows=[]
for case in cases:
 summaries=[x['summary'] for x in v['case_summaries'] if x['case']==case]
 rows.append({'case':case,'warning_windows':sum(s['warning_windows'] for s in summaries),'eligible_monitor_windows':sum(s['monitor_eligible_windows'] for s in summaries),
  'warning_episodes':sum(len(s['warning_episodes']) for s in summaries),'sequences_with_warning':sum(s['any_warning'] for s in summaries),
  'signal_windows':sum(s['signal_windows'] for s in summaries),'signal_associations':sum(s['associated_search_signal_windows'] for s in summaries),
  'hypothetical_gate_losses':sum(s['associated_recoveries_that_warning_gate_would_remove'] for s in summaries),
  'noise_alarms':sum(s['search_noise_alarms'] for s in summaries),'noise_windows':sum(s['search_noise_windows'] for s in summaries),
  'available_search_windows':sum(s['detector_available_windows'] for s in summaries)})
with (out/'summary_counts.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'aggregate_summary.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
s=impact['summary']
lines=['# Corrected cadence and warning-impact control','',f"Run `{out.name}`, sourcecc68b93. [Pre-run protocol](protocol.md), [numerical verification](independent_verification.json), [independent paired-impact analysis](independent_impact_analysis.json), [counts](summary_counts.csv).",'',
'## Fixed design and correction','',
'Five FRESH seeds21101..21105;800continuous four-second observations per seed, eight paired cases, unchanged default causal BackgroundRegimeMonitor and frozen core calibration (85A/199B,FAP0.005). No automatic reference update, detector gating, source selection or threshold tuning. The5seeds are replication units; windows/cases are dependent and paired.',
'', 'Repeated injections use the64-slot grid288..792every8windows. Sparse injections use every fifth slot starting at grid position2:13slots304..784. Family is global-grid parity, so both schedules alternate AND every shared probe has identical family,waveform,noise,SNR18,support and core output. The old inconsistent study is retained; new counts are not improvements relative to its different seeds/schedule.',
'', 'The monitor has11warmup target windows; search decisions remain available. Warnings annotate the same frozen measurement records and never suppress candidates. The2400sage policy and all previous adaptive machinery are absent from this study.',
'', '## All cases: descriptive counts','',
'| Case | Warnings / monitor-eligible windows | Episodes | Search noise alarms / noise windows | Signal support associations | Hypothetical gate losses |',
'| --- | ---: | ---: | ---: | ---: | ---: |']
for r in rows:
 lines.append(f"| {r['case']} | {r['warning_windows']}/{r['eligible_monitor_windows']} | {r['warning_episodes']} | {r['noise_alarms']}/{r['noise_windows']} | {r['signal_associations']}/{r['signal_windows']} | {r['hypothetical_gate_losses']} |")
lines += ['', 'Zero-denominator signal groups are not efficiency estimates. Warning windows within an episode are not independent alarms; persistent-tone observations represent one sustained source per seed, not hundreds of independent events. Actual search availability is unchanged in every case. No formal superiority/FAP-confidence inference from these five paired sequences.', '',
'## A. Corrected shared-probe cadence control','',
f"All{s['shared_probe_pairs']}shared probe pairs have identical input and core outputs. Warning flags differ on{s['shared_warning_disagreements']}shared positions. Of{s['shared_core_associations']}core support associations, hypothetical warning gating removes{s['shared_sparse_gate_losses']}in the sparse context and{s['shared_repeated_gate_losses']}in the repeated context. These differences, if present, arise from monitor history with the same current observation, not a changed search statistic. A null result is retained; it does not prove history never matters.",
'', 'Full additional repeated-probe counts are in the table; do not compare raw totals65vs320as if dose/denominators were equal.', '',
'## B. Changed measurements versus search impact','',
f"Across5160paired target comparisons (stable->level_step and sparse->step_with_probes), max search-score difference is{s['max_gain_score_difference']}; nomination disagreements{s['gain_nomination_disagreements']}. The same factor-two gain step can trigger persistent warnings while leaving this locally normalized search statistic unchanged. This is a verified property of this constructed gain intervention, not a proof for arbitrary instrument response changes.",
'', 'The line profile can instead strongly increase stale-calibration noise nominations. Similar generic warning flags therefore do not quantify the same effect on detector reliability. No automatic harmless/dangerous warning classifier was trained here; effects are measured with paired synthetic controls unavailable in ordinary deployment.', '',
'## C. Did the injected signal add a new association?','',
'For every positive probe, score the exact pre-injection noise candidate list against the SAME truth support after detection. The algorithm never sees this support. Counts are not physical waveform identification.', '',
'| Positive case | Probe pairs | Positive association | Noise already associated | New association | Lost association | Positive-only nomination |',
'| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
for case,g in s['injection_support'].items():
 lines.append('| '+' | '.join(str(x) for x in [case,g['pairs'],g['positive_association'],g['noise_already_associated'],g['new_association'],g['lost_association'],g['positive_only_nomination']])+' |')
lines += ['', 'These columns are not a disjoint partition: noise-already includes any lost association. A positive candidate also present in noise is insufficient evidence that the injection was identified. This is an oracle diagnostic, not a new detector acceptance rule.', '',
'## D. Fixed future noise after warning episodes','',
'For each onset, observe the NEXT12windows, excluding the onset itself. Overlapping ranges are unioned and end-of-run ranges truncated. This is retrospective descriptive association, not validated prediction or causation.', '',
'| Noise-only case | Unioned future alarms / windows | Actual exposure(s) |',
'| --- | ---: | ---: |']
for case in ['stable','level_step','spectral_line']:
 group=[r for r in impact['future_noise'] if r['case']==case];n=sum(r['future_noise_windows'] for r in group);a=sum(r['future_noise_alarms'] for r in group);t=sum(r['future_exposure_seconds'] for r in group)
 lines.append(f'| {case} | {a}/{n} | {t} |')
lines += ['', 'An empty future set has unavailable conditional rate(null), not zero. Do not credit a warning as successful reliability control merely because it overlaps a known generated change. No warning gating was performed.', '',
'## Ambiguity and limits','',
'Spectral-line and desired-persistent-tone inputs, core outputs and monitor warnings remain exactly identical; only evaluation labels differ. Warning cannot identify physical noise origin or justify absorbing a persistent signal into background. No real waveform, astrophysical novelty, GW170814 causal explanation or real-instrument lifetime claim.',
'', 'This corrected pilot has two chosen noise changes,SNR18families,five jointFFTcolored seed traces and one fixed observation horizon. Noise windows are notIID; source-calibration variation remains. Equal nominalFAP is not equal achievedFAR. Any future automatic impact policy needs a new frozen rule, application budget and independent validation, retaining false alerts, missed/abstained signals and persistent-signal ambiguity.', '',
'## Verification and preservation','',
f"Numerical checks pass{v['reference_models']}reference models,{v['detector_calibrations']}calibrations,{v['monitor_records']}monitor/search/sidecar records each,{v['candidate_records']}candidate records,{v['regenerated_inputs']}regenerated input records,{v['manual_numpy_feature_checks']}manualNumPyFFTfeatures and{v['production_prefix_checks']}prefix comparisons. KSinteger-path oracle checks{v['ks_feature_checks']}comparisons with maxBonferroniperror{v['max_ks_p_absolute_error']:.3g}. Shared65probeidentity and195positive/noise pairs independently reconstructed. Five future-exposure fixtures test duplicates,overlap,empty andtruncation. These are scoped checks,notpytestpasses orproof ofstatistical validity.",
'', 'Initial/source/reference/calibration hashes and all search/candidate payloads preserved. Old variants untouched. No excluded cases, threshold/seed retuning, production/app/test edits, O3/data/manifest/network/realwaveform access, publication or adoption. No measurements rerun to improve a result. Retain any implementation/error records separately from scientific outcomes.', '',
'![Warning and search impact](warning_search_impact.png)', '',
'![Positive versus pre-injection noise support associations](probe_noise_association.png)', '',
'## Decision','',
'Use the monitor as a provisional annotation of changed measurements. Do not equate every warning with invalid detector calibration or every support-associated positive with identification of the added signal. No new detector or automatic recalibration policy is adopted by this study. L09/L10 remain unresolved.']
(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Saved corrected-schedule impact report and all-case counts.')
