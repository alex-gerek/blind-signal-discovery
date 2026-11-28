"""Report the matched exploratory run without changing methods or outputs."""
from pathlib import Path
from itertools import product
from collections import Counter
import csv,gzip,json,sys
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1])
def load(p):
 if p.suffix=='.gz':
  with gzip.open(p,'rt',encoding='utf-8') as f:return json.load(f)
 return json.loads(p.read_text(encoding='utf-8'))
assert load(out/'independent_verification.json')['status']=='PASS';rows=load(out/'evaluation.json.gz');arms=['W0D0','W0D1','W1D0','W1D1'];summary={};strata=[];timing={}
def metrics(rs,arm):return {'n':len(rs),'anomaly':sum(r['result']['anomaly'] for r in rs),'claims':sum(r['decisions'][arm]['claim'] for r in rs),'strict':sum(r['decisions'][arm]['strict_correct'] for r in rs),'loose':sum(r['decisions'][arm]['loose_correct'] for r in rs),'wrong_strict':sum(r['decisions'][arm]['wrong_pair'] for r in rs),'wrong_loose':sum(r['decisions'][arm]['wrong_pair_loose'] for r in rs),'false_pairs':sum(r['decisions'][arm]['false_pair'] for r in rs),'stages':dict(Counter(r['decisions'][arm]['stage'] for r in rs))}
for arm in arms:
 rs=[r for r in rows if r['w']==int(arm[1])];summary[arm]={}
 for part,kind in [('matched','two'),('matched','one'),('matched','noise'),('independent','one'),('independent','noise')]:summary[arm][part+'_'+kind]=metrics([r for r in rs if r['meta']['part']==part and r['meta']['kind']==kind],arm)
 main=[r for r in rs if r['meta']['kind']=='two']
 for factor in ['offset','relative_phase','ratio']:
  for value in sorted({r['meta'][factor] for r in main}):
   m=metrics([r for r in main if r['meta'][factor]==value],arm);strata.append({'arm':arm,'factor':factor,'value':value,**{k:v for k,v in m.items() if k!='stages'}})
 successful=[r for r in main if r['decisions'][arm]['strict_correct']];timing[arm]={'n':len(successful),'median_abs':{k:float(np.median([abs(r['decisions'][arm]['timing'][k]) for r in successful])) if successful else None for k in ['onset','end','midpoint']}}
with (out/'stratified_counts.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(strata[0]));w.writeheader();w.writerows(strata)
result={'summary':summary,'timing_conditional_on_strict_success':timing,'note':'216variants share12base noise/carrier/global-phase draws; variant counts are not independent sample sizes'};(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
trans=list(csv.DictReader((out/'paired_transitions.csv').open(encoding='utf-8')))
lines=['# Matched 2 Hz grid, phase and amplitude interventions','',f"Run `{out.name}`; source `{load(out/'manifest.json')['source_commit']}`. [Frozen protocol](protocol.md), [verification](independent_verification.json), [paired contrasts](paired_transitions.csv), [marginal counts](stratified_counts.csv).",'',
'## What is new','',
'BSD-20/21 previously found five loose / four strict gap-2Hz successes, but noise/carrier/global phase changed across cases. Here12fresh base units retain identical noise, even base carrier and global phase across18variants:2grid offsets x3relative phases x3amplitude ratios. Gap2Hz, resultant root-energySNR18, one-second Tukey envelope and onset1.25s fixed. There are216variant cases, not216independent base draws. Both W banks and D diagnostics are unchanged. No source separation, music-stem reconstruction or physical identification is performed.','',
'Offsets0and0.5Hz shift both carriers, preserving gap/noise/phases. Relative phase is at onset, not constant in time. Fixing resultant energy changes individual amplitudes when phase or ratio changes; shifting carrier in colored noise can change its local noise level. These are bounded matched interventions, not isolated universal phase or FFT-grid effects.','',
'New anomaly calibrations (284noise,85A/199B) and pair controls (199noise+199single) use disjoint streams and per-arm thresholds, frozen before any evaluation inference. No tuning, added seeds or exclusions. The same12noise bases supply12matched noise and24matched single controls;160new independent noise and160new single controls are also reported separately.','',
'Primary success: at most one distinct frequency pair, BOTH errors<=0.5Hz under one-to-one matching. Retain secondary<=1Hz. Exact duplicate estimates are invalid; a claim on one/noise is false regardless of how close one estimate is to the tone. Truth is evaluator-only.','',
'## Main matched results','',
'| Arm | Anomalies /216 | Pair claims | Correct <=0.5Hz | Correct <=1Hz | Wrong under0.5Hz | Wrong under1Hz |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
for arm in arms:
 s=summary[arm]['matched_two'];lines.append(f"| {arm} | {s['anomaly']} | {s['claims']} | {s['strict']} | {s['loose']} | {s['wrong_strict']} | {s['wrong_loose']} |")
for factor,title in [('offset','Grid offset (Hz)'),('relative_phase','Relative phase at onset (radians)'),('ratio','Weaker/stronger amplitude ratio')]:
 lines += ['',f'### {title}','','Strict successes; each denominator is a number of VARIANTS sharing12bases, not independent trials.','','| Value (variants) | W0D0 | W0D1 | W1D0 | W1D1 |','| --- | ---: | ---: | ---: | ---: |']
 for value in sorted({r['value'] for r in strata if r['factor']==factor}):
  ss=[next(s for s in strata if s['arm']==arm and s['factor']==factor and s['value']==value) for arm in arms];lines.append(f"| {value:.8g} ({ss[0]['n']}) | "+' | '.join(str(s['strict']) for s in ss)+' |')
lines+=['','## Within-base contrasts','',
'Each paired contrast has12independent BASE UNITS while fixing the other factors. Both=success in both settings; a-only/b-only=transition; neither=failure in both. The complete132contrast rows preserve all arms, phases and ratios. The following grid contrasts focus on W1D1. No confirmatory p-values or best-cell selection.','',
'| Fixed phase / ratio | Both | On-grid only | Off-grid only | Neither |','| --- | ---: | ---: | ---: | ---: |']
for r in trans:
 if r['arm']=='W1D1' and r['factor']=='offset':
  fixed=json.loads(r['fixed']);lines.append(f"| {fixed['relative_phase']:.6g} / {fixed['ratio']:g} | {r['both']} | {r['a_only']} | {r['b_only']} | {r['neither']} |")
lines+=['','## False claims and missed stages','',
'| Arm | Independent noise false pairs /160 | Independent single false pairs /160 | Matched noise false pairs /12 | Matched single false pairs /24 |','| --- | ---: | ---: | ---: | ---: |']
for arm in arms:lines.append('| '+arm+' | '+' | '.join(str(summary[arm][name]['false_pairs']) for name in ['independent_noise','independent_one','matched_noise','matched_one'])+' |')
lines+=['','Independent-control exact95% two-sided intervals (conditional on the frozen calibration/population; methods share inputs):','']
for arm in arms:
 for name in ['independent_noise','independent_one']:
  k=summary[arm][name]['false_pairs'];n=160;lo=beta.ppf(.025,k,n-k+1) if k else 0;hi=beta.ppf(.975,k+1,n-k) if k<n else 1;lines.append(f'- {arm} {name}: {k}/{n}, interval[{100*lo:.3f}%,{100*hi:.3f}%].')
lines+=['','Matched controls are not pooled with the independent controls. Zero of160does not validate1% false-pair control or a real-sensor guarantee. Pair-cutoff sentinel semantics remain as before; null means missing-proposal/-infinity, not a failed calibration. The manifest and frozen thresholds retain actual cutoffs.','']
for arm in arms:lines.append(f"- {arm} main decision stages: {summary[arm]['matched_two']['stages']}. Wrong claimed pairs are retained in addition to these availability/proposal stages.")
lines+=['','## Timing and runtime','',
'Timing on correct frequency estimates only, with different selected subsets; no paired accuracy claim from these medians. Observed support[1.25,2.25]s, midpoint1.75s. Boundary/onset robustness is not tested anew.','',
'| Arm | Strict success n | Median absolute onset/end/midpoint error(s) |','| --- | ---: | --- |']
for arm,t in timing.items():lines.append(f"| {arm} | {t['n']} | {t['median_abs']} |")
for w in [0,1]:
 rs=[r for r in rows if r['w']==w];bank=[r['result']['bank_seconds'] for r in rs];d1=[r['result']['methods']['1']['seconds'] for r in rs];lines+=['',f"W{w}: bank median/p95={1000*np.median(bank):.3f}/{1000*np.quantile(bank,.95):.3f}ms; D1 proposal median/p95={1000*np.median(d1):.3f}/{1000*np.quantile(d1,.95):.3f}ms. Shared re-labeling and archive I/O are not contained in individual proposal timers. These are mixed-case warm-process measurements, not end-to-end streaming latency."]
lines+=['','## Evidence and limitations','',
'All calibration/evaluation records are saved directly as JSON gzip archives.504matched-condition NPZ traces retain waveform,complexSTFT,power,background,Z,mask,labels; independent negative-control full maps were not saved. Large traces remain local/gitignored with hashes. Historical artifacts and production source hashes are unchanged. The source pin, scripts, calibration data and generator metadata permit separate fresh-folder reproduction.','',
f"Verification PASS: {load(out/'independent_verification.json')}. These are scoped record/map checks, not pytest counts or proof of scientific validity.", '',
'No measurement failure, excluded seed, invalid case, parameter adjustment, interruption or waveform-study rerun. Old matching fixtures and all216generator pairing/energy fixtures passed before fitting. These fixtures generate waveforms but make no held-out inference or selection; thresholds were frozen before evaluation decisions. Selected outputs match direct production multi-scale decisions; generator hashes reproduce.','',
'One calibration realization and12base units limit generalization. Strong resultant energy and fixed onset/duration mean the study does not answer weak-signal or variable-duration questions. Multiple contrasts are exploratory; even matched causal intervention results apply only to this specified generator/configuration. Inspect cell-level successes and failures before proposing any different diagnostic. No automatic production adoption, source reconstruction or stem-separation conclusion.','',
'![Strict pair successes by matched setting](matched_gap2_successes.png)','',
'New results remain local pending review; no commit/push,mainmerge,release,Zenodo or O3 access.']
(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
fig,axes=plt.subplots(2,2,figsize=(9,6.5),layout='constrained');phases=[0.,np.pi/2,np.pi];ratios=[1.,.5,.25]
for i,w in enumerate([0,1]):
 for j,offset in enumerate([0.,.5]):
  arm=f'W{w}D1';ax=axes[i,j];matrix=np.array([[sum(r['decisions'][arm]['strict_correct'] for r in rows if r['w']==w and r['meta']['kind']=='two' and r['meta']['offset']==offset and r['meta']['relative_phase']==phase and r['meta']['ratio']==ratio) for ratio in ratios] for phase in phases]);im=ax.imshow(matrix,vmin=0,vmax=12,cmap='Blues')
  for yi,xi in product(range(3),range(3)):ax.text(xi,yi,f'{matrix[yi,xi]}/12',ha='center',va='center',color='white' if matrix[yi,xi]>6 else 'black')
  ax.set_xticks(range(3),['1','0.5','0.25']);ax.set_yticks(range(3),['0','90','180']);ax.set_xlabel('Amplitude ratio');ax.set_ylabel('Relative phase at onset (degrees)');ax.set_title(f'{"0.5" if w==0 else "1"} s +peaks; offset{offset:g}Hz')
fig.suptitle('Matched 2 Hz tones: both errors <=0.5Hz\nEach cell uses the SAME12base units; resultant SNR18 fixed');fig.colorbar(im,ax=axes.ravel().tolist(),label='Strict successes /12');fig.savefig(out/'matched_gap2_successes.png',dpi=170);plt.close(fig)
print(json.dumps(result))
