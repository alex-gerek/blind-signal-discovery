"""Truth-blind per-frame pair persistence on saved maps; descriptive only."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations,product
from collections import Counter,defaultdict
import csv,gzip,hashlib,json,sys,uuid
import numpy as np
from scipy.signal import find_peaks
ROOT=Path(__file__).resolve().parents[1];PARENT=ROOT/'results/matched_gap2/20261128T115035Z_bf957291';PROTOCOL=ROOT/'docs/temporal_peak_persistence_protocol_2026-11-28.md'
def load(p):return json.load(gzip.open(p,'rt',encoding='utf-8')) if p.suffix=='.gz' else json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):
 text=json.dumps(x,indent=2,allow_nan=False)+'\n'
 if p.suffix=='.gz':
  with gzip.open(p,'wt',encoding='utf-8') as f:f.write(text)
 else:p.write_text(text,encoding='utf-8')
def track(frames,times,indices,jitter):
 prev={};best=None;per_frame=[]
 for j in indices:
  current={}
  for p in frames[j]['pairs']:
   key=(p[0],p[1]);choices=[]
   for da,db in (product([-1,0,1],repeat=2) if jitter else [(0,0)]):
    pk=(key[0]+da,key[1]+db)
    if pk in prev:choices.append((prev[pk],abs(da)+abs(db),pk))
   chosen=min(choices,key=lambda c:(-c[0]['length'],c[1],c[2]))[0] if choices else None
   state={'length':1+(chosen['length'] if chosen else 0),'path':(chosen['path'] if chosen else [])+[{'frame':int(j),'time':float(times[j]),'frequencies':[int(key[0]),int(key[1])]}],'endpoint_score':p[5]};current[key]=state
   if best is None or (-state['length'],-state['endpoint_score'],key,j)<(-best['length'],-best['endpoint_score'],tuple(best['path'][-1]['frequencies']),best['path'][-1]['frame']):best=state
  per_frame.append(max((s['length'] for s in current.values()),default=0));prev=current
 if best is None:return {'length':0,'path':[],'center_span_s':0.,'union_window_support_s':0.,'frequency_span_hz':[0,0],'max_component_step_hz':0,'per_frame_longest_ending_run':per_frame}
 path=best['path'];f=np.array([p['frequencies'] for p in path]);best.update(center_span_s=path[-1]['time']-path[0]['time'],union_window_support_s=1+path[-1]['time']-path[0]['time'],frequency_span_hz=np.ptp(f,axis=0).astype(int).tolist(),max_component_step_hz=int(np.abs(np.diff(f,axis=0)).max()) if len(f)>1 else 0,per_frame_longest_ending_run=per_frame);return best

def summarize(frames,times,indices):
 present=[j for j in indices if frames[j]['pairs']];freqs=Counter(tuple(p[:2]) for j in indices for p in frames[j]['pairs']);selected={j:min(frames[j]['pairs'],key=lambda p:(-p[5],p[0],p[1])) for j in present};jumps=[]
 for a,b in zip(indices,indices[1:]):
  if a in selected and b in selected:jumps.append({'from_frame':int(a),'to_frame':int(b),'from_frequencies':selected[a][:2],'to_frequencies':selected[b][:2],'component_jumps_hz':[abs(selected[b][k]-selected[a][k]) for k in [0,1]]})
 return {'frames_total':len(indices),'frames_with_pairs':len(present),'unique_frequency_pairs':len(freqs),'frequencies_and_frame_counts':[{'frequencies':list(k),'frames':v} for k,v in sorted(freqs.items())],'exact':track(frames,times,indices,False),'jitter1':track(frames,times,indices,True),'strongest_per_frame':[{'frame':j,'frequencies':selected[j][:2],'score':selected[j][5]} for j in present],'strongest_pair_jumps':jumps}
def main():
 out=ROOT/'results/temporal_peak_persistence'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True);print('OUTPUT',out,flush=True);(out/'driver.py').write_bytes(Path(__file__).read_bytes());(out/'protocol.md').write_bytes(PROTOCOL.read_bytes());parent=load(PARENT/'artifact_sha256.json')['files'];assert all(sha(ROOT/p)==h for p,h in parent.items());raw=load(PARENT/'evaluation.json.gz');inputs=[r for r in raw if r['w']==1 and r['meta']['part']=='matched'];assert len(inputs)==252
 records=[];case_rows=[];freq_rows=[];frame_rows=[]
 for ri,row in enumerate(inputs):
  m=row['meta'];d=row['decisions']['W1D1'];group='two_correct' if d['strict_correct'] else 'two_wrong' if m['kind']=='two' and d['claim'] else 'two_no_pair' if m['kind']=='two' else m['kind'];frames=[]
  with np.load(PARENT/row['trace']) as archive:
   z=archive['z'];f=archive['freqs'];times=archive['times'];assert z.shape==(511,25) and np.array_equal(np.diff(f),np.ones(510));assert np.allclose(np.diff(times),.125)
   # Full-map inference: no metadata, event boundaries or true frequencies used here.
   for j in range(len(times)):
    y=z[:,j];peaks=find_peaks(y)[0];high=[int(i) for i in peaks if y[i]>=4];pairs=[];spacing_rejected=valley_rejected=0
    for a,b in combinations(high,2):
     if b-a<2:spacing_rejected+=1;continue
     valley=float(y[a+1:b].min())
     if valley>.8*min(y[a],y[b]):valley_rejected+=1;continue
     pairs.append([int(f[a]),int(f[b]),float(y[a]),float(y[b]),valley,float(min(y[a],y[b]))])
    frames.append({'frame':j,'time_s':float(times[j]),'all_local_peak_count':len(peaks),'above4_peaks':[[int(f[i]),float(y[i])] for i in high],'pairs':pairs,'spacing_rejected':spacing_rejected,'valley_rejected':valley_rejected})
  scopes={'all':frames,'near2to3':[dict(fr,pairs=[p for p in fr['pairs'] if 2<=p[1]-p[0]<=3]) for fr in frames]};intervals={'full':list(range(25)),'central9':[j for j,t in enumerate(times) if 1.25<=t<=2.25]};assert len(intervals['central9'])==9
  summaries={}
  for view,frs in scopes.items():
   summaries[view]={}
   for interval,indices in intervals.items():
    s=summarize(frs,times,indices);summaries[view][interval]=s;js=s['strongest_pair_jumps'];largest=max((max(x['component_jumps_hz']) for x in js),default=0)
    case_rows.append({'id':row['id'],'base':m['base'],'group':group,'view':view,'interval':interval,'frames':len(indices),'frames_with_pairs':s['frames_with_pairs'],'exact_max_run':s['exact']['length'],'jitter_max_run':s['jitter1']['length'],'jitter_frequency_span_hz':s['jitter1']['frequency_span_hz'],'max_strongest_pair_jump_hz':largest,'strongest_pair_transitions':len(js),'strongest_transitions_over1hz':sum(max(x['component_jumps_hz'])>1 for x in js)})
    for freqcount in s['frequencies_and_frame_counts']:freq_rows.append({'id':row['id'],'group':group,'view':view,'interval':interval,'frequency1':freqcount['frequencies'][0],'frequency2':freqcount['frequencies'][1],'frame_count':freqcount['frames']})
  # Evaluation only: metadata labels are first used AFTER full-map pair lists/tracks are fixed.
  associated=[]
  for fr in frames:
   if m['kind']=='two':pairs=[p for p in fr['pairs'] if abs(p[0]-m['f1'])<=.5 and abs(p[1]-m['f2'])<=.5]
   elif m['kind']=='one':pairs=[p for p in fr['pairs'] if min(abs(p[0]-m['f1']),abs(p[1]-m['f1']))<=.5]
   else:pairs=[]
   associated.append(dict(fr,pairs=pairs))
  annotation={interval:summarize(associated,times,indices) for interval,indices in intervals.items()}
  prior_selected=None
  if group=='two_wrong':
   fs=d['frequencies'];fsframes=[dict(fr,pairs=[p for p in fr['pairs'] if p[:2]==fs]) for fr in frames];prior_selected={'frequencies':fs,'summaries':{interval:summarize(fsframes,times,idx) for interval,idx in intervals.items()}}
  record={'id':row['id'],'group':group,'meta':m,'trace':row['trace'],'times':times.tolist(),'frames':frames,'summaries':summaries,'carrier_annotation':annotation,'prior_wrong_pair':prior_selected};records.append(record)
  for fr in frames:frame_rows.append({'id':row['id'],'group':group,'frame':fr['frame'],'time_s':fr['time_s'],'peaks_above4':len(fr['above4_peaks']),'pairs_all':len(fr['pairs']),'pairs_near2to3':sum(2<=p[1]-p[0]<=3 for p in fr['pairs'])})
  if ri%60==0:print('ANALYZED',ri+1,252,flush=True)
 for name,rs in [('case_metrics.csv',case_rows),('pair_frequencies.csv',freq_rows),('frame_counts.csv',frame_rows)]:
  with (out/name).open('w',encoding='utf-8',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
 save(out/'records.json.gz',records)
 aggregate=[]
 for group in ['two_correct','two_wrong','two_no_pair','one','noise']:
  for view,interval in product(['all','near2to3'],['full','central9']):
   rs=[r for r in records if r['group']==group];ss=[r['summaries'][view][interval] for r in rs]
   a={'group':group,'view':view,'interval':interval,'cases':len(rs),'frames_with_pairs_median':float(np.median([s['frames_with_pairs'] for s in ss])),'frames_with_pairs_range':[min(s['frames_with_pairs'] for s in ss),max(s['frames_with_pairs'] for s in ss)]}
   for key,label in [('exact','exact'),('jitter1','jitter')]:
    runs=[s[key]['length'] for s in ss];a[label+'_run_median']=float(np.median(runs));a[label+'_run_max']=max(runs)
    for k in [1,2,3,5]:a[label+'_atleast'+str(k)]=sum(x>=k for x in runs)
   aggregate.append(a)
 save(out/'aggregate.json',aggregate)
 assert all(sha(ROOT/p)==h for p,h in parent.items());save(out/'manifest.json',{'status':'COMPLETE','source_commit':__import__('subprocess').check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_artifacts_verified':len(parent),'input_maps':252,'frames':6300,'new_waveforms':0,'new_calibration':False,'protocol_sha256':sha(PROTOCOL),'driver_sha256':sha(Path(__file__)),'source_trace_sha256':{r['trace']:sha(PARENT/r['trace']) for r in inputs},'parent_evaluation_sha256':sha(PARENT/'evaluation.json.gz'),'groups':dict(Counter(r['group'] for r in records))})
 print('COMPLETE',out,flush=True)
 for a in aggregate:
  if a['view']=='near2to3' and a['interval']=='central9':print(json.dumps(a),flush=True)
if __name__=='__main__':main()
