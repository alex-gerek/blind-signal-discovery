"""Independent saved-map pair and temporal recurrence verification."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import gzip,json,hashlib,sys
import numpy as np
from scipy.signal import find_peaks
ROOT=Path(__file__).resolve().parents[1];PARENT=ROOT/'results/matched_gap2/20261128T115035Z_bf957291'
def load(p):return json.load(gzip.open(p,'rt',encoding='utf-8')) if p.suffix=='.gz' else json.loads(p.read_text(encoding='utf-8'))
def oracle(frames,indices,jitter):
 prev=np.empty((0,2),int);lengths=np.empty(0,int);best=0;per=[]
 for j in indices:
  current=np.array([p[:2] for p in frames[j]['pairs']],int).reshape(-1,2)
  if len(current) and len(prev):
   allowed=(np.max(np.abs(current[:,None,:]-prev[None,:,:]),axis=2)<=int(jitter));now=1+np.max(np.where(allowed,lengths[None,:],0),axis=1)
  else:now=np.ones(len(current),int)
  maximum=int(now.max()) if len(now) else 0;per.append(maximum);best=max(best,maximum);prev=current;lengths=now
 return best,per
def checksummary(stored,frames,times,indices):
 assert stored['frames_total']==len(indices) and stored['frames_with_pairs']==sum(bool(frames[j]['pairs']) for j in indices)
 counts=Counter(tuple(p[:2]) for j in indices for p in frames[j]['pairs']);assert counts==Counter({tuple(x['frequencies']):x['frames'] for x in stored['frequencies_and_frame_counts']})
 for jitter,name in [(False,'exact'),(True,'jitter1')]:
  best,per=oracle(frames,indices,jitter);s=stored[name];assert s['length']==best and s['per_frame_longest_ending_run']==per;path=s['path'];assert len(path)==best
  for point in path:assert point['frame'] in indices and times[point['frame']]==point['time'] and tuple(point['frequencies']) in {tuple(p[:2]) for p in frames[point['frame']]['pairs']}
  for a,b in zip(path,path[1:]):assert b['frame']==a['frame']+1 and max(abs(a['frequencies'][k]-b['frequencies'][k]) for k in [0,1])<=int(jitter)
  if path:assert np.isclose(s['center_span_s'],(best-1)*.125) and np.isclose(s['union_window_support_s'],1+(best-1)*.125)
 # Independently verify strongest-per-frame selection and jumps.
 selected={j:min(frames[j]['pairs'],key=lambda p:(-p[5],p[0],p[1])) for j in indices if frames[j]['pairs']}
 assert stored['strongest_per_frame']==[{'frame':j,'frequencies':selected[j][:2],'score':selected[j][5]} for j in selected]
 for jump in stored['strongest_pair_jumps']:
  a,b=jump['from_frame'],jump['to_frame'];assert b==a+1 and a in selected and b in selected;assert jump['component_jumps_hz']==[abs(selected[b][k]-selected[a][k]) for k in [0,1]]
def main():
 out=Path(sys.argv[1]);rs=load(out/'records.json.gz');assert len(rs)==252;checks={'maps':0,'frames':0,'pair_gates':0,'summary_recurrences':0};near_annotations={};wrong_controls=[]
 for r in rs:
  with np.load(PARENT/r['trace']) as a:
   z=a['z'];freq=a['freqs'];times=a['times'];assert len(times)==25
   for j,fr in enumerate(r['frames']):
    indices=find_peaks(z[:,j])[0];high=[int(i) for i in indices if z[i,j]>=4];assert [[int(freq[i]),float(z[i,j])] for i in high]==fr['above4_peaks'];expected=[]
    for i,k in combinations(high,2):
     if k-i<2:continue
     valley=float(z[i+1:k,j].min())
     if valley<=.8*min(z[i,j],z[k,j]):expected.append([int(freq[i]),int(freq[k]),float(z[i,j]),float(z[k,j]),valley,float(min(z[i,j],z[k,j]))])
    assert expected==fr['pairs'];checks['frames']+=1;checks['pair_gates']+=len(expected)
  checks['maps']+=1;intervals={'full':list(range(25)),'central9':list(range(6,15))}
  for view in ['all','near2to3']:
   frames=[dict(fr,pairs=[p for p in fr['pairs'] if view=='all' or 2<=p[1]-p[0]<=3]) for fr in r['frames']]
   for interval,idx in intervals.items():checksummary(r['summaries'][view][interval],frames,times,idx);checks['summary_recurrences']+=2
  m=r['meta'];annotated=[]
  for fr in r['frames']:
   if m['kind']=='two':ps=[p for p in fr['pairs'] if abs(p[0]-m['f1'])<=.5 and abs(p[1]-m['f2'])<=.5]
   elif m['kind']=='one':ps=[p for p in fr['pairs'] if min(abs(p[0]-m['f1']),abs(p[1]-m['f1']))<=.5]
   else:ps=[]
   annotated.append(dict(fr,pairs=ps))
  for interval,idx in intervals.items():checksummary(r['carrier_annotation'][interval],annotated,times,idx);checks['summary_recurrences']+=2
  # Supplemental annotation restricted to the declared near-pair structural class.
  ann=[dict(fr,pairs=[p for p in fr['pairs'] if 2<=p[1]-p[0]<=3]) for fr in annotated];idx=intervals['central9'];near_annotations[r['id']]={'frames_with_carrier_pairs':sum(bool(ann[j]['pairs']) for j in idx),'exact_max':oracle(ann,idx,False)[0],'jitter_max':oracle(ann,idx,True)[0],'frequencies':dict(Counter(str(p[:2]) for j in idx for p in ann[j]['pairs']))}
  if r['prior_wrong_pair']:
   target=r['prior_wrong_pair']['frequencies'];controls=[x for x in rs if x['meta']['base']==m['base'] and (x['meta']['kind']=='noise' or x['meta']['kind']=='one' and x['meta']['offset']==m['offset'])]
   item={'id':r['id'],'frequencies':target,'groups':{}}
   for control in [r]+controls:
    frames=[dict(fr,pairs=[p for p in fr['pairs'] if p[:2]==target]) for fr in control['frames']];above=[j for j in idx if any(p[0]==target[1] for p in control['frames'][j]['above4_peaks'])];item['groups'][control['meta']['kind']]={'pair_frames':sum(bool(frames[j]['pairs']) for j in idx),'pair_exact_run':oracle(frames,idx,False)[0],'second_peak_above4_frames':above}
   wrong_controls.append(item)
 for path,h in load(out/'manifest.json')['source_trace_sha256'].items():assert hashlib.sha256((PARENT/path).read_bytes()).hexdigest()==h
 assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in load(PARENT/'artifact_sha256.json')['files'].items())
 (out/'near_carrier_annotations.json').write_text(json.dumps(near_annotations,indent=2)+'\n',encoding='utf-8');(out/'wrong_pair_control_persistence.json').write_text(json.dumps(wrong_controls,indent=2)+'\n',encoding='utf-8');checks.update(status='PASS',parent_artifact_hashes_unchanged=529,new_waveforms=0);(out/'independent_verification.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8');print(json.dumps(checks))
if __name__=='__main__':main()
