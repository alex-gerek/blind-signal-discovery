"""Summarize the frozen pilot; no measurement or retuning."""
from pathlib import Path
from collections import defaultdict
import csv,json,sys
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'));rows=load(out/'evaluation.json');v=load(out/'independent_verification.json');assert v['status']=='PASS' and load(out/'reproducibility_verification.json')['status']=='PASS'
arms=['W0D0','W0D1','W1D0','W1D1'];names=['0.5 s / boxes','0.5 s / +peaks','1 s / boxes','1 s / +peaks']
def ci(k,n):return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]
def summarize(rs,arm):
 n=len(rs);claims=sum(r['decisions'][arm]['claim'] for r in rs);correct=sum(r['decisions'][arm]['correct'] for r in rs)
 return dict(n=n,anomalies=sum(r['result']['anomaly'] for r in rs),pair_claims=claims,correct=correct,strict_correct=sum(r['decisions'][arm]['strict_correct'] for r in rs),wrong_pair=sum(r['decisions'][arm]['wrong_pair'] for r in rs),false_pair=sum(r['decisions'][arm]['false_pair'] for r in rs),ambiguous=sum(r['decisions'][arm]['ambiguous_match'] for r in rs),no_pair=n-claims,correct_ci95=ci(correct,n))
summary={};strata=[];cells=[];timings={};wrong=[];discord=[]
for arm in arms:
 w=int(arm[1]);rs=[r for r in rows if r['w']==w];summary[arm]={}
 for part in ['main','edge']:
  for kind in ['noise','one','two']:
   subset=[r for r in rs if r['part']==part and r['meta']['kind']==kind]
   if subset:summary[arm][part+'_'+kind]=summarize(subset,arm)
 main=[r for r in rs if r['part']=='main' and r['meta']['kind']=='two']
 for factor in ['gap','ratio','phase','offset','onset','snr']:
  for value in sorted({r['meta'][factor] for r in main}):
   subset=[r for r in main if r['meta'][factor]==value];strata.append({'arm':arm,'factor':factor,'value':value,**{k:v for k,v in summarize(subset,arm).items() if not isinstance(v,list)}})
 for r in main:
  dec=r['decisions'][arm];cells.append({'arm':arm,'id':r['id'],**{k:r['meta'][k] for k in ['gap','ratio','phase','offset','onset','snr','f1','f2','global_phase']},**{k:dec[k] for k in ['claim','correct','strict_correct','wrong_pair','ambiguous_match']},'estimated_frequencies':str(dec['frequencies'])})
  if dec['wrong_pair']:wrong.append({'arm':arm,'id':r['id'],'truth':[r['meta']['f1'],r['meta']['f2']],'estimate':dec['frequencies'],'gap':r['meta']['gap'],'ratio':r['meta']['ratio'],'snr':r['meta']['snr']})
 good=[r for r in main if r['decisions'][arm]['correct']];timings[arm]={'correct_pair_n':len(good),'all_pair_claim_n':sum(r['decisions'][arm]['claim'] for r in main)}
 for axis in ['onset_error','end_error','midpoint_error']:
  values=[r['decisions'][arm]['timing'][axis] for r in good];timings[arm][axis]={'signed_median':float(np.median(values)) if values else None,'absolute_median':float(np.median(np.abs(values))) if values else None,'absolute_p95':float(np.quantile(np.abs(values),.95)) if values else None}
for a,b in [('W0D0','W0D1'),('W1D0','W1D1'),('W0D0','W1D0'),('W0D1','W1D1')]:
 aa={r['id']:r['decisions'][a]['correct'] for r in rows if r['part']=='main' and r['meta']['kind']=='two' and a in r['decisions']};bb={r['id']:r['decisions'][b]['correct'] for r in rows if r['part']=='main' and r['meta']['kind']=='two' and b in r['decisions']};assert aa.keys()==bb.keys();discord.append({'a':a,'b':b,'a_only':sum(aa[k] and not bb[k] for k in aa),'b_only':sum(bb[k] and not aa[k] for k in aa),'both':sum(aa[k] and bb[k] for k in aa)})
def write_csv(name,data):
 with (out/name).open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
write_csv('stratified_counts.csv',strata);write_csv('factorial_cases.csv',cells);write_csv('paired_discordances.csv',discord)
runtime={}
for w in [0,1]:
 rs=[r for r in rows if r['w']==w];bank=[r['result']['bank_seconds'] for r in rs];runtime[f'W{w}']={'first_bank_ms':bank[0]*1000,'later_bank_median_ms':float(np.median(bank[1:])*1000),'later_bank_p95_ms':float(np.quantile(bank[1:],.95)*1000),'trace_io_total_seconds':sum(r['trace_seconds'] for r in rs)}
 for d in [0,1]:
  values=[r['result']['methods'][str(d)]['seconds'] for r in rs];runtime[f'W{w}'][f'D{d}_median_ms']=float(np.median(values)*1000);runtime[f'W{w}'][f'D{d}_p95_ms']=float(np.quantile(values,.95)*1000)
size=sum(p.stat().st_size for p in (out/'traces').glob('*.npz'))
result={'summary':summary,'timing_correct_pairs_only':timings,'runtime':runtime,'wrong_pairs':wrong,'paired_discordances':discord,'trace_archive_bytes':size}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
lines=['# Blind two-tone exploratory comparison','',f"Run `{out.name}`. Source `cc68b93`. [Frozen protocol](protocol.md), [independent verification](independent_verification.json), [reproducibility checks](reproducibility_verification.json).",'',
'## Scope and frozen rules','',
'Four research arms: W0=historical bank with 0.5 s longest window; W1=replace longest window with 1 s, keep 0.125 s hop and the other scales. D0=disjoint accepted longest-scale clusters with overlapping time boxes; D1=also search multiple peaks within one accepted longest-scale cluster. No production changes. The 1 s configuration changes frequency smoothing bandwidth as well as window duration; this is not an isolated duration ablation.', '',
'Truth-blind estimators receive only the observation and frozen calibration. D1 uses a cluster-selected mean spectrum, local maxima >=4, an interior-bin valley <=0.8 times the smaller peak, and no true-frequency neighborhoods. Proposals are ranked by the smaller peak height, exact duplicate frequency pairs deduplicated, equal estimates forbidden. At most one pair is returned. Correctness requires one-to-one matching with both errors <=1 Hz; the stricter gap/4 result is retained. Every claim on noise or a single tone is a false two-tone claim. No separate waveforms are reconstructed.', '',
'Independent streams: 284 anomaly-calibration noise observations,199 pair-calibration noise and199 pair-calibration single-tone observations. W configurations fitted separately, all four pair cutoffs frozen before held-out evaluation. In all four arms the chosen pair-null order statistic was missing-proposal (-infinity), serialized as null with explicit semantics. Thus any finite eligible proposal passes the pair cutoff; missing proposals never pass. This is an outcome of sparse null proposals, not unavailable calibration or a measured zero false rate. Stored null distributions and empirical p checks are retained.', '',
'No threshold/gate was changed after evaluation. One independent calibration draw per W; all arms share evaluation inputs for pairing. Historical 96/98 results were used for development and are not claimed as new held-out evidence.', '',
'## Main held-out results','',
'432 two-tone observations:6 gaps x3 relative phases x3 amplitude ratios x2 grid offsets x2 onsets x2 resultant SNRs. Each full factorial cell has ONE independent noise/carrier draw; marginal counts mix other conditions and have limited precision. Additional160noise and160single-tone controls are independent. Two-tone resultant energy is matched at linear root-energy SNR8/18; components are not each SNR18 as in the historical study.', '',
'| Arm | Anomalies on two tones /432 | Pair claims | Correct <=1Hz | Correct gap/4 | Wrong pair | No pair |', '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
for arm in arms:
 s=summary[arm]['main_two'];lines.append(f"| {arm} | {s['anomalies']} | {s['pair_claims']} | {s['correct']} | {s['strict_correct']} | {s['wrong_pair']} | {s['no_pair']} |")
lines+=['','| Gap (Hz), n=72 each | W0D0 | W0D1 | W1D0 | W1D1 |','| --- | ---: | ---: | ---: | ---: |']
for gap in [1,2,4,8,16,96]:lines.append('| '+str(gap)+' | '+' | '.join(str(next(r['correct'] for r in strata if r['arm']==a and r['factor']=='gap' and r['value']==gap)) for a in arms)+' |')
lines+=['','The W1D1 gap2 result is5/72 at1Hz tolerance,4/72 at stricter0.5Hz tolerance. At gap1 no arm recovered a correct pair. Some methods return a wrong pair on a gap1 input; this does not constitute successful separation. All ambiguous matches and errors are retained.','',
'| Arm | Noise anomalies /160 | Noise false pairs /160 | Single-tone anomalies /160 | Single-tone false pairs /160 |','| --- | ---: | ---: | ---: | ---: |']
for arm in arms:
 n=summary[arm]['main_noise'];s=summary[arm]['main_one'];lines.append(f"| {arm} | {n['anomalies']} | {n['false_pair']} | {s['anomalies']} | {s['false_pair']} |")
upper=ci(0,160)[1];lines+=['',f'Zero false pair claims in160inputs has an exact two-sided95% upper bound of{100*upper:.2f}% per method/population. It does not validate a1% false-pair bound conditionally on this calibration, or any real-data guarantee. Do not pool dependent arms to inflate the denominator. W0 also nominated2/160noise windows as anomalies; pair abstention does not erase those anomaly false alarms.','',
'## Stratified results and failures','',
'[stratified_counts.csv](stratified_counts.csv) reports counts separately by gap, ratio, relative phase, grid offset, onset and SNR. [factorial_cases.csv](factorial_cases.csv) retains each complete cell. Offset0 is on both grids; offset0.5Hz is off both (half the1Hz grid, quarter the2Hz grid). Actual frequency remainders are saved in evaluation metadata. [paired_discordances.csv](paired_discordances.csv) preserves wins/losses on identical observations, without a confirmatory multiple-comparison claim.','']
for factor,title in [('ratio','Amplitude ratio'),('phase','Relative phase (radians)'),('offset','Base frequency-grid offset (Hz)'),('onset','Onset (seconds)'),('snr','Resultant root-energy SNR')]:
 lines += [f'### {title}','','| Value (denominator) | W0D0 | W0D1 | W1D0 | W1D1 |','| --- | ---: | ---: | ---: | ---: |']
 for value in sorted({r['value'] for r in strata if r['factor']==factor}):
  ss=[next(r for r in strata if r['arm']==a and r['factor']==factor and r['value']==value) for a in arms];lines.append(f"| {value:.6g} ({ss[0]['n']}) | "+' | '.join(str(s['correct']) for s in ss)+' |')
 lines+=['']
lines+=['Wrong frequency pairs (one example can fail in more than one arm):','']
for r in wrong:lines.append(f"- {r['arm']}, input{r['id']}: truth{r['truth']}, estimate{r['estimate']}; gap{r['gap']}, ratio{r['ratio']}, SNR{r['snr']}.")
lines+=['','## Timing and record boundaries','',
'Main onsets compare exact hop alignment1.25s with1.3115234375s, about half a hop later. Signals span many overlapping windows. Boundary/offset rows above retain lost cases. Correct-pair timing uses the union of contributing component boxes and mean component peak time, compared with envelope support and its midpoint; midpoint is not a true physical intensity peak. The following statistics are CONDITIONAL on correct pairs and use different subsets; they do not prove one method has better timing.','',
'| Arm | Correct-pair n | Median absolute onset error(s) | End error(s) | Midpoint error(s) |','| --- | ---: | ---: | ---: | ---: |']
for arm in arms:
 t=timings[arm];lines.append(f"| {arm} | {t['correct_pair_n']} | {t['onset_error']['absolute_median']:.4f} | {t['end_error']['absolute_median']:.4f} | {t['midpoint_error']['absolute_median']:.4f} |")
lines+=['','Separate record-edge supplement:24pairs start at-.25or3.5s. The full tone is normalized first, then clipped to the four-second record without padding. Observed energy is lower; these combine truncation, SNR loss and boundary effects, so are not a pure timing ablation. No out-of-record signal is inferred.','',
'| Arm | Edge anomalies /24 | Edge claims | Correct pairs |','| --- | ---: | ---: | ---: |']
for arm in arms:
 s=summary[arm]['edge_two'];lines.append(f"| {arm} | {s['anomalies']} | {s['pair_claims']} | {s['correct']} |")
lines+=['','## Runtime and storage','',
'Perf-counter measurements in one process, with numerical thread environment variables fixed to1. First evaluation call is reported separately; it is NOT a fresh-process cold benchmark because calibration already initialized the libraries. Later medians/p95 mix different signal cases and are exploratory measurements on this machine. Calibration time, verification, and artifact tracing/I/O are not online inference.','',
'| W | First evaluation bank(ms) | Later bank median/p95(ms) | D0 median/p95(ms) | D1 median/p95(ms) |','| --- | ---: | --- | --- | --- |']
for w,r in runtime.items():lines.append(f"| {w} | {r['first_bank_ms']:.3f} | {r['later_bank_median_ms']:.3f}/{r['later_bank_p95_ms']:.3f} | {r['D0_median_ms']:.3f}/{r['D0_p95_ms']:.3f} | {r['D1_median_ms']:.3f}/{r['D1_p95_ms']:.3f} |")
lines+=['',
'Bank time includes production STFT/normalization/clustering and collection of raw candidates/maps. D0/D1 times cover proposal extraction; shared extra labeling for trace/within-cluster use is outside these individual timers. Therefore summing columns is not an exact end-to-end single-arm latency. No separate true cold-process or deployment throughput measurement was made. Peak working-set query was unavailable (null); do not invent a memory figure. Detailed scale and I/O timers remain in evaluation.json.','',f'Compressed traces:1552NPZ files, {size/1024**2:.1f}MiB, excluded from git by existing *.npz rules. Artifact hashes index their persistent local paths. Full longest-scale waveform/complex STFT/power/background/Z/mask/labels retained; shorter-scale raw scores/candidates retained but full shorter-scale maps require deterministic regeneration. The four-second method is offline;1swindow support does not imply a measured streaming latency.','',
'## Verification and preserved errors','',
'Initial attempt20261127T192027Z_8e7afefd failed on an incorrectly specified ambiguous-matching fixture before any calibration/evaluation. Its driver, manifest and failure trace remain intact. The fixture was corrected, not the method; successful run uses a fresh folder. No subsequent excluded/invalid/near-zero cases, threshold changes or waveform-study reruns.','',
f"Independent verification PASS: {v['calibrations']}calibrations, {v['pair_calibration_records']}pair-calibration records, {v['evaluation_records']}evaluation records, {v['trace_files']}trace maps, {v['proposal_sets']}proposal sets, {v['claim_checks']}claim checks and {v['matching_checks']}matching checks. Reproduction verifies1458unique inputs and6input-only repeated-inference/label-change checks; production parity checks10records. These are scoped arithmetic/reproducibility checks, not pytest counts or scientific proof.", '',
'Production source and all previously recorded two-oscillation artifacts remain hash-identical. Calibration thresholds remained frozen throughout evaluation. No O3, external scientific data, production/app changes, commit, push, merge, release or Zenodo action.','',
'## Interpretation and next decision','',
'The within-cluster diagnostic adds correct frequency estimates in this bounded synthetic sample. One-second plus within-cluster peaks sometimes estimates a2Hz pair, but most close pairs remain unreturned; there are wrong pair estimates, losses across arms and worse onset/end precision on the conditional subsets. Zero control false claims is underpowered evidence, not a validated low-false-rate instrument. Keep the original detector unchanged. Before adoption, independently replicate calibrations, improve null-control coverage and power, and separate representation/normalization/gating contributions. This does not reconstruct two source waveforms, identify physical sources or demonstrate music stem separation.','',
'![Correct frequency-pair estimates by gap](correct_pairs_by_gap.png)','',
'![Main counts and false-claim uncertainty](pair_claim_summary.png)','']
(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
# Scientific figures from frozen counts only.
fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');x=np.arange(6);width=.19
for i,(arm,name) in enumerate(zip(arms,names)):
 vals=[next(r['correct'] for r in strata if r['arm']==arm and r['factor']=='gap' and r['value']==g) for g in [1,2,4,8,16,96]];bars=ax.bar(x+(i-1.5)*width,vals,width,label=name);ax.bar_label(bars,fontsize=8,padding=2)
ax.set_xticks(x,[1,2,4,8,16,96]);ax.set_ylim(0,40);ax.set_xlabel('True frequency separation (Hz), evaluation only');ax.set_ylabel('Correct estimated pairs /72');ax.set_title('Blind frequency estimates: both errors <=1 Hz\nAt 2 Hz, 1 s +peaks gives5/72; stricter0.5Hz matching gives4/72');ax.legend(fontsize=8,ncol=2);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);fig.savefig(out/'correct_pairs_by_gap.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained');x=np.arange(4);a=axes[0]
a.bar(x,[summary[s]['main_two']['correct'] for s in arms],label='Correct pair');a.bar(x,[summary[s]['main_two']['wrong_pair'] for s in arms],bottom=[summary[s]['main_two']['correct'] for s in arms],color='#d55e00',label='Wrong pair');a.set_xticks(x,arms);a.set_ylabel('Pair claims /432');a.legend(fontsize=8);a.set_title('Main two-tone sample');a=axes[1];a.errorbar(x,np.zeros(4),yerr=[np.zeros(4),np.repeat(upper*100,4)],fmt='o',capsize=5);a.set_xticks(x,arms);a.set_ylim(-.2,3);a.set_ylabel('False two-tone claims (%)');a.set_title('Each control group:0/160\nExact95% interval, each method');fig.suptitle('An anomaly, a frequency-pair estimate and source reconstruction are different');fig.savefig(out/'pair_claim_summary.png',dpi=170);plt.close(fig)
print(json.dumps({'summary':summary,'runtime':runtime,'trace_MiB':size/1024**2,'report':str(out/'report.md')}))
