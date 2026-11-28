"""Saved-record verifier and matched contrasts; no threshold fitting or tuning."""
from pathlib import Path
from itertools import combinations,permutations,product
import csv,gzip,hashlib,importlib.util,json,sys
import numpy as np
from scipy import ndimage,signal
ROOT=Path(__file__).resolve().parents[1]
def load(p):
 if p.suffix=='.gz':
  with gzip.open(p,'rt',encoding='utf-8') as f:return json.load(f)
 return json.loads(p.read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-10)
def main():
 out=Path(sys.argv[1]);manifest=load(out/'manifest.json');assert manifest['status']=='COMPLETE';rows=load(out/'evaluation.json.gz');prs=load(out/'pair_calibration.json.gz');raw=load(out/'calibration_raw.json.gz');freeze=load(out/'frozen_thresholds.json');assert hashlib.sha256((out/'frozen_thresholds.json').read_bytes()).hexdigest()==manifest['freeze_sha256'];assert freeze['heldout_decisions']==0
 checks={'calibrations':0,'trace_maps':0,'claims':0,'matching':0,'regenerated_inputs':0,'matched_bases':12}
 for w in [0,1]:
  cal=load(out/f'calibration_W{w}.json');rs=[r for r in raw if r['w']==w];assert len(rs)==284;stats=cal['scale_stats']
  for name,s in stats.items():
   values=[r['raw'][name]['max'] for r in rs[:85]];equal(s['median_null_score'],np.median(values));equal(s['spread_scale'],max(np.quantile(values,.95)-np.median(values),1.))
  null=sorted(max((r['raw'][name]['max']-s['median_null_score'])/s['spread_scale'] for name,s in stats.items()) for r in rs[85:]);equal(null,cal['sorted_global_null_scores']);assert null[-1]==cal['global_threshold'];checks['calibrations']+=1
  for d in [0,1]:
   cutoff=[]
   for kind in ['noise','one']:
    ns=[r for r in prs if r['w']==w and r['meta']['kind']==kind];assert len(ns)==199;scores=sorted(r['result']['methods'][str(d)]['proposals'][0]['score'] if r['result']['methods'][str(d)]['proposals'] else -np.inf for r in ns);cutoff.append(scores[197])
   cut=freeze['cuts'][f'W{w}D{d}']['cutoff'];assert max(cutoff)==(-np.inf if cut is None else cut)
 for r in rows+prs:
  w=r['w'];result=r['result'];cal=load(out/f'calibration_W{w}.json');stats=cal['scale_stats'];score=max((v['max']-stats[n]['median_null_score'])/stats[n]['spread_scale'] for n,v in result['raw'].items());equal(score,result['score']);assert result['anomaly']==(score>cal['global_threshold'])
  expected=[]
  for name,v in result['raw'].items():
   for c in v['candidates']:
    if (c['anomaly_score']-stats[name]['median_null_score'])/stats[name]['spread_scale']>cal['global_threshold']:expected.append((name,c['anomaly_score'],c['peak_freq'],c['peak_time']))
  assert sorted(expected)==sorted((c['scale'],c['anomaly_score'],c['peak_freq'],c['peak_time']) for c in result['accepted'])
  if 'decisions' not in r:continue
  for arm,dec in r['decisions'].items():
   ps=result['methods'][arm[-1]]['proposals'];assert len(ps)==len({tuple(p['frequencies']) for p in ps});assert ps==sorted(ps,key=lambda p:(-p['score'],p['frequencies'],p['start']));cut=freeze['cuts'][arm]['cutoff'];claim=bool(ps and ps[0]['score']>(-np.inf if cut is None else cut));assert claim==dec['claim'] and dec['false_pair']==(claim and r['meta']['kind']!='two');checks['claims']+=1
   if r['meta']['kind']=='two':
    truth=[r['meta']['f1'],r['meta']['f2']];f=dec['frequencies'];errs=[[abs(a-b) for a,b in zip(f,p)] for p in permutations(truth)] if claim else [[np.inf,np.inf]];assert dec['strict_correct']==any(max(e)<=.5 for e in errs) and dec['loose_correct']==any(max(e)<=1 for e in errs);checks['matching']+=1
  if r['trace']:
   with np.load(out/r['trace']) as a:
    assert hashlib.sha256(a['waveform'].astype('<f8').tobytes()).hexdigest()==r['meta']['hash'];equal(abs(a['stft'])**2,a['power']);psd=np.maximum(ndimage.median_filter(np.median(a['power'],axis=1)/np.log(2),size=21,mode='nearest'),1e-18);tb=np.maximum(ndimage.median_filter(np.median(a['power']/psd[:,None],axis=0)/np.log(2),size=7,mode='nearest'),1e-18);equal(psd,a['psd']);equal(tb,a['time_background']);equal(a['z'],a['power']/psd[:,None]/tb[None,:]);assert np.array_equal(a['mask'],a['z']>4);labels,_=ndimage.label(a['z']>4,np.ones((3,3),int));assert np.array_equal(labels,a['labels']);checks['trace_maps']+=1
 # Regeneration and explicit within-base pairing: inference code never receives these labels.
 spec=importlib.util.spec_from_file_location('matched',ROOT/'experiments/run_matched_gap2.py');driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
 seen=set()
 for r in rows:
  ident=r['id']
  if ident in seen:continue
  seen.add(ident);m=r['meta']
  if m['part']=='matched':x,_=driver.generate(m['base'],m['kind'],m['offset'],m['relative_phase'] or 0,m['ratio'] or 1)
  else:x,_=driver.old.make(m['seed'],m['index'],m['kind'])
  assert driver.old.key(x)==m['hash'];checks['regenerated_inputs']+=1
 for base in range(12):
  group=[r for r in rows if r['w']==0 and r['meta'].get('base')==base];two=[r for r in group if r['meta']['kind']=='two'];assert len(two)==18;assert len({r['meta']['noise_hash'] for r in group})==len({r['meta']['carrier'] for r in group})==len({r['meta']['global_phase'] for r in group})==1
  for r in two:equal(r['meta']['observed_snr'],18);equal(r['meta']['e1']+r['meta']['e2']+r['meta']['cross'],324)
 assert len(rows)==1144 and checks['regenerated_inputs']==572 and checks['trace_maps']==504
 transitions=[];arms=['W0D0','W0D1','W1D0','W1D1']
 for arm in arms:
  w=int(arm[1]);sub=[r for r in rows if r['w']==w and r['meta']['kind']=='two'];lookup={(r['meta']['base'],r['meta']['offset'],r['meta']['relative_phase'],r['meta']['ratio']):r for r in sub};assert len(lookup)==216
  def contrast(factor,lo,hi,fixed):
   counts={'both':0,'a_only':0,'b_only':0,'neither':0}
   for base in range(12):
    values={'base':base,**fixed,factor:lo};a=lookup[tuple(values[k] for k in ['base','offset','relative_phase','ratio'])]['decisions'][arm]['strict_correct'];values[factor]=hi;b=lookup[tuple(values[k] for k in ['base','offset','relative_phase','ratio'])]['decisions'][arm]['strict_correct'];counts['both' if a and b else 'a_only' if a else 'b_only' if b else 'neither']+=1
   assert sum(counts.values())==12;transitions.append({'arm':arm,'factor':factor,'a':lo,'b':hi,'fixed':json.dumps(fixed,sort_keys=True),'bases':12,**counts})
  for phase,ratio in product([0.,np.pi/2,np.pi],[1.,.5,.25]):contrast('offset',0.,.5,{'relative_phase':phase,'ratio':ratio})
  for offset,ratio in product([0.,.5],[1.,.5,.25]):
   for phase in [np.pi/2,np.pi]:contrast('relative_phase',0.,phase,{'offset':offset,'ratio':ratio})
  for offset,phase in product([0.,.5],[0.,np.pi/2,np.pi]):
   for ratio in [.5,.25]:contrast('ratio',1.,ratio,{'offset':offset,'relative_phase':phase})
 with (out/'paired_transitions.csv').open('w',encoding='utf-8',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(transitions[0]));wr.writeheader();wr.writerows(transitions)
 for p,h in {**manifest['source_sha256'],**manifest['historical_sha256']}.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 checks.update(status='PASS',paired_contrasts=len(transitions),source_and_historical_hashes_unchanged=True,thresholds_unchanged=True);(out/'independent_verification.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8');print(json.dumps(checks))
if __name__=='__main__':
 try:main()
 except Exception:
  import traceback
  p=Path(sys.argv[1])/'verification_failure.json'
  if not p.exists():p.write_text(json.dumps({'traceback':traceback.format_exc()},indent=2),encoding='utf-8')
  raise
