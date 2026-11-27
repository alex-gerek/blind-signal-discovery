"""Independent accounting and saved-map verifier for the blind pilot."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,sys
import numpy as np
from scipy import ndimage,signal
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def equal(a,b):np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-10)
def main():
 out=Path(sys.argv[1]).resolve();m=load(out/'manifest.json');assert m['status']=='COMPLETE';rows=load(out/'evaluation.json');prs=load(out/'pair_calibration_records.json');raw=load(out/'anomaly_calibration_raw.json');frozen=load(out/'frozen_pair_thresholds.json');cuts=frozen['cuts'];stats={};count={'calibrations':0,'pair_calibration_records':len(prs),'evaluation_records':len(rows),'trace_files':0,'proposal_sets':0,'claim_checks':0,'matching_checks':0}
 assert hashlib.sha256((out/'frozen_pair_thresholds.json').read_bytes()).hexdigest()==m['threshold_freeze_sha256'] and frozen['heldout_evaluated']==0
 for w in [0,1]:
  cal=load(out/f'anomaly_calibration_W{w}.json')['calibration'];stats[w]=cal['scale_stats'];rs=[r for r in raw if r['w']==w];assert len(rs)==284
  for name,s in stats[w].items():
   x=np.array([r['raw'][name]['max'] for r in rs[:85]]);med=float(np.median(x));spread=max(float(np.quantile(x,.95))-med,1);equal([s['median_null_score'],s['spread_scale']],[med,spread])
  maxima=sorted(max((r['raw'][n]['max']-s['median_null_score'])/s['spread_scale'] for n,s in stats[w].items()) for r in rs[85:]);equal(maxima,cal['sorted_global_null_scores']);assert cal['global_threshold']==maxima[-1];count['calibrations']+=1
  for d in [0,1]:
   combined=[]
   for kind in ['noise','one']:
    nr=[r for r in prs if r['w']==w and r['meta']['kind']==kind];assert len(nr)==199
    scores=sorted(r['result']['methods'][str(d)]['proposals'][0]['score'] if r['result']['methods'][str(d)]['proposals'] else -np.inf for r in nr);v=scores[197];stored=cuts[f'W{w}D{d}']['by_group'][kind];assert (-np.inf if stored is None else stored)==v;combined.append(v)
   cutoff=cuts[f'W{w}D{d}']['cutoff'];assert (-np.inf if cutoff is None else cutoff)==max(combined)
 # Acceptance and full-observation pair accounting on every calibration and evaluation record.
 for r in prs+rows:
  w=r['w'];res=r['result'];cal=load(out/f'anomaly_calibration_W{w}.json')['calibration'];maximum=max((v['max']-stats[w][n]['median_null_score'])/stats[w][n]['spread_scale'] for n,v in res['raw'].items());equal(maximum,res['score']);assert res['anomaly']==(maximum>cal['global_threshold'])
  expected=[]
  for name,v in res['raw'].items():
   s=stats[w][name]
   for c in v['candidates']:
    z=(c['anomaly_score']-s['median_null_score'])/s['spread_scale']
    if z>cal['global_threshold']:expected.append((name,c['peak_freq'],c['peak_time'],c['anomaly_score']))
  assert sorted(expected)==sorted((c['scale'],c['peak_freq'],c['peak_time'],c['anomaly_score']) for c in res['accepted'])
  for d in [0,1]:
   ps=res['methods'][str(d)]['proposals'];assert len({tuple(p['frequencies']) for p in ps})==len(ps)
   assert ps==sorted(ps,key=lambda p:(-p['score'],p['frequencies'],p['start']))
   for p in ps:assert len(p['frequencies'])==2 and p['frequencies'][0]<p['frequencies'][1]
   count['proposal_sets']+=1
  if 'decisions' not in r:continue
  # Reconstruct the saved longest-scale representation, labeling and D0/D1 proposals.
  with np.load(out/r['trace']) as a:
   count['trace_files']+=1;z=a['z'];f=a['freqs'];t=a['times'];labels,n=ndimage.label(z>4,np.ones((3,3),int));assert np.array_equal(labels,a['labels']) and np.array_equal(z>4,a['mask']);assert hashlib.sha256(a['waveform'].astype('<f8').tobytes()).hexdigest()==r['meta']['hash']
   equal(abs(a['stft'])**2,a['power']);psd=ndimage.median_filter(np.median(a['power'],axis=1)/np.log(2),size=21,mode='nearest');psd=np.maximum(psd,1e-18);env=ndimage.median_filter(np.median(a['power']/psd[:,None],axis=0)/np.log(2),size=7,mode='nearest');env=np.maximum(env,1e-18);equal(psd,a['background_psd']);equal(env,a['background_time']);equal(z,a['power']/psd[:,None]/env[None,:])
   cs=[c for c in res['accepted'] if c['scale']==('long_500ms' if w==0 else 'long_1000ms')]
   base=[]
   def add(pairs,fa,fb,h1,h2,start,end,peak,origin):
    if fa!=fb:pairs.append({'frequencies':sorted([float(fa),float(fb)]),'score':float(min(h1,h2)),'start':float(start),'end':float(end),'peak':float(peak),'origin':origin})
   for i,j in combinations(range(len(cs)),2):
    aa,bb=cs[i],cs[j]
    if (aa['f_max']<bb['f_min'] or bb['f_max']<aa['f_min']) and max(aa['t_start'],bb['t_start'])<=min(aa['t_end'],bb['t_end']):add(base,aa['peak_freq'],bb['peak_freq'],aa['peak_excess_power'],bb['peak_excess_power'],min(aa['t_start'],bb['t_start']),max(aa['t_end'],bb['t_end']),(aa['peak_time']+bb['peak_time'])/2,'disjoint')
   extra=[]
   for c in cs:
    lab=labels[np.argmin(abs(f-c['peak_freq'])),np.argmin(abs(t-c['peak_time']))];cols=np.where(np.any(labels==lab,axis=0))[0];profile=z[:,cols].mean(axis=1);pks=signal.find_peaks(profile)[0];pks=[p for p in pks if c['f_min']<=f[p]<=c['f_max'] and profile[p]>=4]
    for i,j in combinations(pks,2):
     if j-i>=2 and profile[i+1:j].min()<=.8*min(profile[i],profile[j]):add(extra,f[i],f[j],profile[i],profile[j],c['t_start'],c['t_end'],c['peak_time'],'within')
   for d in [0,1]:
    ordered=sorted(base+(extra if d else []),key=lambda p:(-p['score'],p['frequencies'],p['start']));unique={}
    for p in ordered:unique.setdefault(tuple(p['frequencies']),p)
    assert list(unique.values())==res['methods'][str(d)]['proposals']
  for arm,dec in r['decisions'].items():
   d=arm[-1];ps=res['methods'][d]['proposals'];p=ps[0] if ps else None;cut=cuts[arm]['cutoff'];claim=bool(p and p['score']>(-np.inf if cut is None else cut));assert dec['claim']==claim;assert dec['frequencies']==(p['frequencies'] if claim else None);assert dec['false_pair']==(claim and r['meta']['kind']!='two');count['claim_checks']+=1
   if p is not None:
    group_ps=[]
    for kind in ['noise','one']:
     ns=frozen['null_scores'][str(r['w'])][d][kind];group_ps.append((1+sum((-np.inf if x is None else x)>=p['score'] for x in ns))/200)
    assert claim==all(x<=.01 for x in group_ps)
   if r['meta']['kind']=='two':
    truth=[r['meta']['f1'],r['meta']['f2']];fs=dec['frequencies'];good=[];strict=[]
    for perm in permutations(truth):
     errors=[abs(a-b) for a,b in zip(fs,perm)] if fs else [np.inf,np.inf];good.append(max(errors)<=1);strict.append(max(errors)<=r['meta']['gap']/4)
    assert dec['correct']==any(good) and dec['strict_correct']==any(strict) and dec['ambiguous_match']==all(good);count['matching_checks']+=1
   if claim and r['meta']['kind']!='noise':
    meta=r['meta'];equal(list(dec['timing'].values()),[p['start']-meta['observed_start'],p['end']-meta['observed_end'],p['peak']-(meta['observed_start']+meta['observed_end'])/2])
 assert len(rows)==1552 and len({r['id'] for r in rows})==776 and count['trace_files']==1552
 # Inference results/proposals were produced without meta: verify exact input-only function signatures.
 import ast
 tree=ast.parse((out/'driver.py').read_text());funcs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
 assert [a.arg for a in funcs['extract'].args.args]==['x','w','cal'] and [a.arg for a in funcs['proposals'].args.args]==['candidates','z','freqs','times','labels','within']
 for name in ['extract','proposals']:
  assert not any(isinstance(n,ast.Name) and n.id in ['meta','truth','params'] for n in ast.walk(funcs[name]))
 count['input_only_interface_checks']=2
 for p,h in {**m['source_sha256'],**m['historical_sha256']}.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 count.update(status='PASS',production_and_historical_hashes_unchanged=True,threshold_hash_unchanged=True)
 (out/'independent_verification.json').write_text(json.dumps(count,indent=2)+'\n',encoding='utf-8');print(json.dumps(count))
if __name__=='__main__':
 try:main()
 except Exception:
  import traceback
  p=Path(sys.argv[1])/'verification_failure.json'
  if not p.exists():p.write_text(json.dumps({'traceback':traceback.format_exc()},indent=2),encoding='utf-8')
  raise
