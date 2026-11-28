"""Report fixed temporal diagnostics; no new detector decisions."""
from pathlib import Path
import csv,gzip,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);load=lambda p:json.load(gzip.open(p,'rt',encoding='utf-8')) if p.suffix=='.gz' else json.loads(p.read_text(encoding='utf-8'));rs=load(out/'records.json.gz');agg=load(out/'aggregate.json');verify=load(out/'independent_verification.json');assert verify['status']=='PASS';annotations=load(out/'near_carrier_annotations.json');wrong=load(out/'wrong_pair_control_persistence.json')
labels={'two_correct':'Previously correct pair','two_wrong':'Previously wrong pair','two_no_pair':'Previously no pair','one':'One-tone control','noise':'Noise-only control'}
lines=['# Temporal persistence of saved spectral peak pairs','',f"Run `{out.name}`. [Protocol](protocol.md), [verification](independent_verification.json), [all case metrics](case_metrics.csv), [all frequency-pair counts](pair_frequencies.csv), [frame counts](frame_counts.csv).",'',
'## Scope and definitions','',
'No detector changes, new waveforms, new maps, threshold fitting or classifier adoption. The252saved1-second maps comprise216two-tone variants,24matched singles and12matched noise controls, sharing12base noise/carrier/global-phase units. The320independent negative-control maps were not retained and are not included. All observations are now development evidence.','',
'For each of25frames over the ENTIRE1..511Hz retained frequency range, find local maxima with normalized power Z>=4. A descriptive pair requires at least one intervening bin and a valley<=0.8*smaller peak. There is no accepted-cluster or anomaly gate. This prevents noise controls from appearing clean merely because the detector accepted no cluster. True frequencies enter annotations AFTER blind enumeration. These candidate structures are not newly asserted two-tone detections.','',
'Two structural views are retained:all frequency separations,and near pairs2..3Hz apart anywhere on the map. Near is a predeclared morphology class, not truth-carrier knowledge. Tracks require consecutive frames with no gaps. Exact tracks repeat the same two frequency bins. Permissive tracks allow each sorted component to move<=1Hz per hop and select a retrospective longest path; branching can stitch alternatives,not establish physical identity.','',
'Report full25frames and the SAMEcentral9frames1.25..2.25s in every group. The central interval matches the earlier event support but is only evaluation annotation; full blind enumeration does not use it. One-second Hann windows hop0.125s,overlap87.5%. A3-frame run spans0.25s between centers and1.25s of union window support. A5-frame run spans0.5s between centers and1.5s of union support. These are strongly dependent views of shared samples,not3or5independent confirmations.','',
'## 1. Frequencies near the true carriers: evaluation-only description','',
'This section uses the known two-tone carriers ONLY to annotate the already extracted pairs:one-to-one errors<=0.5Hz. It cannot be used as a blind selection rule.','',
'| Original mean-spectrum group | Cases | Median qualifying central frames /9 (range) | At least one such frame | Same frequencies >=3frames | Allowing1Hzsteps >=3frames |','| --- | ---: | --- | ---: | ---: | ---: |']
for group in ['two_correct','two_wrong','two_no_pair']:
 a=[r for r in rs if r['group']==group];ss=[annotations[r['id']] for r in a];vals=[s['frames_with_carrier_pairs'] for s in ss];lines.append(f"| {labels[group]} | {len(a)} | {np.median(vals):g} ({min(vals)}..{max(vals)}) | {sum(v>0 for v in vals)} | {sum(s['exact_max']>=3 for s in ss)} | {sum(s['jitter_max']>=3 for s in ss)} |")
lines+=['',
'Example true88/90Hz (base00,equal amplitudes,phase0,on-grid):the exact88/90pair appears in5of9central frames,but the longest uninterrupted run is3frames. With the same base and weaker/stronger ratio0.5 (the old mean-spectrum method returns no pair),88/90appears in2of9frames and holds2consecutive frames. Presence in several nonconsecutive frames is different from an uninterrupted stable pair.','',
'In the24previously correct examples,carrier-associated pairs appear in2..6central frames,median4. Only13/24have a3-frame same-bin carrier-associated run. Thus requiring persistence would discard some original successes as well as retain evidence from previously unreturned cases.','',
'In187previously unreturned cases,140have at least one carrier-associated qualifying frame;24have an exact3-frame run,28with allowed1Hzsteps. This is a post-hoc measurement of morphology,not28new calibrated detections or recovered sources.','',
'## 2. Blind near-pair comparison with controls','',
'The following table uses ANY2..3Hzpair over the full frequency axis,without true-frequency selection. A control pair may lie far from a single actual tone. All groups use identical gates and central9frame positions.','',
'| Group | Cases | Median frames containing any near pair (range) | Exact run >=2 | Exact run >=3 | Exact run >=5 | Drift-allowed run >=3 |','| --- | ---: | --- | ---: | ---: | ---: | ---: |']
for group in labels:
 a=next(a for a in agg if a['group']==group and a['view']=='near2to3' and a['interval']=='central9');lines.append(f"| {labels[group]} | {a['cases']} | {a['frames_with_pairs_median']:g} ({a['frames_with_pairs_range'][0]}..{a['frames_with_pairs_range'][1]}) | {a['exact_atleast2']} | {a['exact_atleast3']} | {a['exact_atleast5']} | {a['jitter_atleast3']} |")
lines+=['',
'Of24single-tone controls,19show at least one blind near-pair frame;one has a3-frame same-bin run. Of12noise controls,7show at least one near-pair frame,two have2-frame exact runs,none3-frame in this small sample. These are candidate patterns,not false-positive rates of a calibrated new detector. The controls share12bases and cannot validate a low-error operating rule.','',
'For single-tone controls,post-hoc anchoring one peak to the actual single carrier leaves7/24with at least one near pair;one has a3-frame exact run. Thus the near-pair pattern need not be just an unrelated remote feature:even a known single tone can acquire a nearby spectral partner.','',
'## 3. Do not hide the full-separation search','',
'If any frequency separation is allowed,EVERYrecord has a3-frame exact pair somewhere in the central interval,including12/12noise-only maps. One-tone controls have5/24exact runs>=5;noise-only has0/12exact>=5 but4/12permissive drift paths>=5. Unrestricted pair persistence is therefore not an automatic signature of two tones. Comparing truth-annotated signal tracks against only remote-noise exclusions would be misleading.','',
'## 4. Jumping strongest peaks versus retrospective stable tracks','',
'Per-frame strongest-pair selection can switch to unrelated carrier regions; the longest-path diagnostic instead explicitly favors continuity. Neither is physical source tracking. The same-bin and1Hz-drift results must not be combined into one reliability number.','',
'| Group | Adjacent central near-pair strongest transitions | Transitions moving >1Hz in either component | Largest jump (Hz) |','| --- | ---: | ---: | ---: |']
for group in labels:
 jumps=[j for r in rs if r['group']==group for j in r['summaries']['near2to3']['central9']['strongest_pair_jumps']];lines.append(f"| {labels[group]} | {len(jumps)} | {sum(max(j['component_jumps_hz'])>1 for j in jumps)} | {max((max(j['component_jumps_hz']) for j in jumps),default=0)} |")
lines+=['',
'Counts aggregate shared variants and are descriptive,not independent transition trials. Large jumps reflect selection of a DIFFERENT spectral pair,not evidence that an oscillator physically jumped hundreds ofHz. pair_frequencies.csv lists every exact pair and frame count; records.json.gz retains individual timestamps and paths,so stable frequencies and switches can be audited per case.','',
'## 5. The previously wrong pairs can persist','',
'| Wrong case | Previously returned frequencies | Two-tone pair frames /9; longest exact run | Same-base single-tone control | Same-base noise control |','| --- | --- | --- | --- | --- |']
for r in wrong:
 def text(kind):s=r['groups'][kind];return f"{s['pair_frames']}/9; {s['pair_exact_run']} consecutive"
 lines.append(f"| {r['id']} | {r['frequencies']} | {text('two')} | {text('one')} | {text('noise')} |")
lines+=['',
'Base04wrong134/139pair can repeat3consecutive frames;the same pair also appears in4central frames in its single-tone control,with a2-frame run. The139Hzpeak alone reachesZ>=4in2central noise-control frames,while134Hzdoes not supply a qualifying companion there. This is consistent with the saved noise/leakage/interference explanation;pair persistence does not remove that ambiguity.','',
'The wrong86/89pair (true86/88) is MOREpersistent:6consecutive full-map frames,5central frames. Yet the second frequency is1Hzwrong. A rule requiring3or5frames alone would keep this error;no persistence threshold is adopted from these results.','',
'## Verification and preservation','',
f"Independent saved-map check PASS: {verify['maps']}maps,{verify['frames']}frames,{verify['pair_gates']}retained pair-gate records,{verify['summary_recurrences']}track recurrences. The recurrence oracle uses pair-coordinate distance matrices rather than the producer’s neighbor lookup. Path membership,frequency drift,run lengths,time arithmetic,strongest selection and parent hashes are checked. No waveform generation or method tuning.", '',
'The source maps,parent evaluation and529parent artifact hashes remain unchanged. Full lists,case/frequency/frame CSVs,central/full summaries,wrong-pair controls,protocol andscripts are saved separately. Figures are illustrative evaluations;true-carrier lines and display crops do not enter extraction. The320independent negative-control maps are absent,not silently regenerated or counted.','',
'Any future temporally persistent detector needs a separate fixed decision rule,independent calibration,new negative-control evaluation and honest exposure accounting. This descriptive analysis does not certify three-frame acceptance,source reconstruction,music stem separation or a false-claim probability. New results local;no commit/push,production changes,O3,mainmerge,releaseorZenodo.','',
'![Example temporal maps](temporal_examples.png)','',
'![Exact run length distributions](run_length_distributions.png)','']
(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
# Fixed examples: complete lists remain available; no example selection affects algorithms.
ids=['b00_two_0_0_1','b00_two_0_0_0.5','b04_two_0.5_0_0.5','b05_two_0_0_0.25','b04_one_0.5_0_1','b04_noise_0_0_1'];fig,axes=plt.subplots(3,2,figsize=(11,8.8),layout='constrained')
for ax,ident in zip(axes.flat,ids):
 r=next(r for r in rs if r['id']==ident);m=r['meta'];center=m['f1'] if m['f1'] is not None else m['carrier']+.5;lo,hi=center-5,center+8
 for fr in r['frames']:
  peaks=[p for p in fr['above4_peaks'] if lo<=p[0]<=hi]
  if peaks:ax.scatter([fr['time_s']]*len(peaks),[p[0] for p in peaks],s=[8+min(p[1],50) for p in peaks],c='gray',alpha=.5)
  for p in fr['pairs']:
   if lo<=p[0]<=p[1]<=hi:ax.plot([fr['time_s']]*2,p[:2],color='#e69f00',lw=.7,alpha=.5)
 if m['kind']!='noise':ax.axhline(m['f1'],ls='--',color='black',lw=1)
 if m['kind']=='two':ax.axhline(m['f2'],ls='--',color='black',lw=1)
 ax.axvspan(1.25,2.25,color='blue',alpha=.035);ax.set_xlim(.45,3.55);ax.set_ylim(lo,hi);ax.set_title(ident+' / '+r['group'],fontsize=9);ax.set_xlabel('Frame center (s)');ax.set_ylabel('Frequency (Hz)')
fig.suptitle('Saved1s maps: gray=peaks Z>=4, orange=qualifying pair endpoints\nBlack dashed=true carriers (evaluation only); adjacent frames overlap87.5%');fig.savefig(out/'temporal_examples.png',dpi=160);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
for ax,view in zip(axes,['near2to3','all']):
 for group in ['two_correct','two_no_pair','one','noise']:
  values=[r['summaries'][view]['central9']['exact']['length'] for r in rs if r['group']==group];counts=np.bincount(values,minlength=10)/len(values);ax.plot(range(10),100*counts,'o-',label=f'{labels[group]} (n={len(values)})')
 ax.set_xticks(range(10));ax.set_xlabel('Longest same-bin consecutive run /9');ax.set_ylabel('Cases (%)');ax.set_title('Near2..3Hz pairs' if view=='near2to3' else 'All frequency separations');ax.grid(alpha=.2)
axes[0].legend(fontsize=7);fig.suptitle('Descriptive persistence; shared12base units, no calibrated new detector');fig.savefig(out/'run_length_distributions.png',dpi=160);plt.close(fig)
print(json.dumps({'report':str(out/'report.md'),'figures':2,'new_detector':False}))
