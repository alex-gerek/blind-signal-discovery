"""Research-only one-missing-frame clustering; production code stays unchanged."""
from collections import defaultdict
import numpy as np
from scipy import ndimage
from bsd.detector import CandidateEvent,TemplateFreeDetector
from bsd.detector_v2 import DEFAULT_V02_SCALES,MultiResolutionDetector

DETECTORS={s.name:TemplateFreeDetector(s.detector_config) for s in DEFAULT_V02_SCALES}

def gap_groups(mask):
 labels,n=ndimage.label(mask,structure=np.ones((3,3),dtype=int))
 parent=list(range(n+1))
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 def join(a,b):
  a,b=find(int(a)),find(int(b))
  if a!=b:parent[max(a,b)]=min(a,b)
 nf,nt=mask.shape
 if nt>=3:
  for df in (-1,0,1):
   lo=max(0,-df);hi=min(nf,nf-df)
   a=labels[lo:hi,:-2];b=labels[lo+df:hi+df,2:];active=(a>0)&(b>0)
   if np.any(active):
    for first,second in np.unique(np.column_stack((a[active],b[active])),axis=0):join(first,second)
 groups=defaultdict(list)
 for f,t in np.argwhere(mask):groups[find(int(labels[f,t]))].append((int(f),int(t)))
 return [sorted(g) for _,g in sorted(groups.items())]

def oracle_groups(mask):
 left=set(map(tuple,np.argwhere(mask)));groups=[]
 for start in sorted(left):
  if start not in left:continue
  left.remove(start);stack=[start];group=[]
  while stack:
   f,t=stack.pop();group.append((f,t))
   neighbors=[(f+df,t+dt) for df in (-1,0,1) for dt in (-1,0,1) if df or dt]
   neighbors += [(f+df,t+dt) for df in (-1,0,1) for dt in (-2,2)]
   for q in neighbors:
    if q in left:left.remove(q);stack.append(q)
  groups.append(sorted(group))
 return sorted(groups)

def focused_checks():
 cases=[('empty',[],0),('singleton',[(2,1)],1),('adjacent',[(2,1),(2,2)],1),('diagonal',[(2,1),(3,2)],1),('one_missing',[(2,1),(2,3)],1),('one_missing_diagonal',[(2,1),(3,3)],1),('two_missing',[(2,1),(2,4)],2),('wide_frequency',[(1,1),(3,3)],2),('cycle',[(2,1),(2,3),(3,2),(3,4)],1)]
 for name,coords,count in cases:
  mask=np.zeros((5,7),dtype=bool)
  for q in coords:mask[q]=True
  actual=gap_groups(mask);assert len(actual)==count,name
  assert sorted(actual)==oracle_groups(mask),name
  assert sum(map(len,actual))==len(coords),name
 z=np.array([[4.,np.nextafter(4.,np.inf),3.9]])
 assert np.array_equal(z>4,[[False,True,False]])
 for seed in range(20):
  mask=np.random.default_rng(seed).random((7,12))<.15
  assert sorted(gap_groups(mask))==oracle_groups(mask)
 # Qualification occurs after linking, so two formerly singleton seeds can qualify.
 mask=np.zeros((3,5),bool);mask[1,0]=mask[1,2]=True
 assert sum(len(g)>=2 for g in gap_groups(mask))==1
 return {'status':'PASS','named_topology_cases':len(cases),'random_oracle_masks':20,'strict_seed_boundary':True,'post_link_qualification':True,'no_double_counting':True}

def raw_pair(observation):
 baseline={};gap={}
 for spec in DEFAULT_V02_SCALES:
  r=DETECTORS[spec.name].analyze(observation,calibration=None,return_maps=True);m=r.tf_maps
  baseline[spec.name]=(r.max_anomaly_score,r.candidates,r.diagnostics)
  z=m.normalized_excess_map;groups=gap_groups(z>4);qualified=[];nf,nt=z.shape
  dt=float(m.times[1]-m.times[0]);df=float(m.freqs[1]-m.freqs[0])
  for g in groups:
   if len(g)<2:continue
   a=np.asarray(g);fi=a[:,0];ti=a[:,1];values=z[fi,ti];peak=int(np.argmax(values))
   t0=max(0.,float(m.times[ti.min()]-.5*dt));t1=min(observation.duration,float(m.times[ti.max()]+.5*dt))
   f0=max(0.,float(m.freqs[fi.min()]-.5*df));f1=min(observation.sample_rate/2,float(m.freqs[fi.max()]+.5*df))
   qualified.append(CandidateEvent(candidate_id=len(qualified)+1,t_start=t0,t_end=t1,f_min=f0,f_max=f1,peak_time=float(m.times[ti[peak]]),peak_freq=float(m.freqs[fi[peak]]),anomaly_score=float(np.sum(values-1)),peak_excess_power=float(values[peak]),n_bins=len(g),empirical_p_value=None,p_value_is_upper_bound=False,diagnostics={'time_span_s':t1-t0,'bandwidth_hz':f1-f0,'mean_excess_power':float(np.mean(values)),'boundary_adjacent':bool(fi.min()==0 or fi.max()==nf-1 or ti.min()==0 or ti.max()==nt-1),'input_kurtosis':r.diagnostics['input_kurtosis'],'input_crest_factor':r.diagnostics['input_crest_factor'],'research_gap_link':True}))
  qualified.sort(key=lambda c:c.anomaly_score,reverse=True)
  maxscore=max((c.anomaly_score for c in qualified),default=0.)
  assert maxscore+1e-12>=r.max_anomaly_score
  diag={**r.diagnostics,'n_raw_clusters':len(groups),'n_qualified_clusters':len(qualified),'research_gap_link':True}
  gap[spec.name]=(maxscore,qualified,diag)
 return {'baseline':baseline,'gap':gap}

class TableDetector(MultiResolutionDetector):
 """Reuse production calibration/acceptance on explicit measured raw features."""
 def __init__(self,table):super().__init__();self.table=table
 def _raw_analyze_all_scales(self,observation):return self.table[observation.sample_id]

def score_table(records,arm):
 return {r['sample_id']:{k:(v,[],{}) for k,v in r[arm].items()} for r in records}
