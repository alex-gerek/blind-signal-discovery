"""Research-only ~1s temporal median; original spectral estimate and topology retained."""
from dataclasses import replace
import numpy as np
from scipy import ndimage
from bsd.detector import CandidateEvent,TemplateFreeDetector
from bsd.detector_v2 import DEFAULT_V02_SCALES,MultiResolutionDetector
DETECTORS={s.name:TemplateFreeDetector(s.detector_config) for s in DEFAULT_V02_SCALES}

def slow_kernel(frame_count,hop_s):
 if frame_count<1 or not np.isfinite(hop_s) or hop_s<=0:raise ValueError('Positive frame count and finite positive hop required')
 wanted=2*int(np.floor(.5/hop_s))+1
 available=frame_count if frame_count%2 else frame_count-1
 return max(1,min(wanted,available))

def slow_envelope(raw_time,hop_s):
 raw_time=np.asarray(raw_time,dtype=float)
 if raw_time.ndim!=1 or not raw_time.size or not np.all(np.isfinite(raw_time)):raise ValueError('Nonempty finite 1D temporal envelope required')
 return np.maximum(ndimage.median_filter(raw_time,size=slow_kernel(raw_time.size,hop_s),mode='nearest'),1e-18)

def slow_maps(m):
 white=m.power_spectrogram/m.background_psd[:,None]
 raw=np.median(white,axis=0)/np.log(2)
 env=slow_envelope(raw,float(m.times[1]-m.times[0]))
 return replace(m,background_time_envelope=env,normalized_excess_map=white/env[None,:])

def cluster_maps(m,duration,sample_rate,diag):
 z=m.normalized_excess_map;labels,n=ndimage.label(z>4,structure=np.ones((3,3),dtype=int));nf,nt=z.shape
 dt=float(m.times[1]-m.times[0]);df=float(m.freqs[1]-m.freqs[0]);qualified=[]
 for lid in range(1,n+1):
  fi,ti=np.where(labels==lid)
  if len(fi)<2:continue
  vals=z[fi,ti];peak=int(np.argmax(vals));t0=max(0.,float(m.times[ti.min()]-.5*dt));t1=min(duration,float(m.times[ti.max()]+.5*dt));f0=max(0.,float(m.freqs[fi.min()]-.5*df));f1=min(sample_rate/2,float(m.freqs[fi.max()]+.5*df))
  qualified.append(CandidateEvent(candidate_id=len(qualified)+1,t_start=t0,t_end=t1,f_min=f0,f_max=f1,peak_time=float(m.times[ti[peak]]),peak_freq=float(m.freqs[fi[peak]]),anomaly_score=float(np.sum(vals-1)),peak_excess_power=float(vals[peak]),n_bins=len(fi),empirical_p_value=None,p_value_is_upper_bound=False,diagnostics={'time_span_s':t1-t0,'bandwidth_hz':f1-f0,'mean_excess_power':float(np.mean(vals)),'boundary_adjacent':bool(fi.min()==0 or fi.max()==nf-1 or ti.min()==0 or ti.max()==nt-1),'input_kurtosis':diag['input_kurtosis'],'input_crest_factor':diag['input_crest_factor'],'research_slow_temporal_median':True}))
 qualified.sort(key=lambda c:c.anomaly_score,reverse=True)
 return (max((c.anomaly_score for c in qualified),default=0.),qualified,{**diag,'n_raw_clusters':n,'n_qualified_clusters':len(qualified),'n_seed_bins':int((z>4).sum()),'research_slow_temporal_median':True})

def raw_pair(observation):
 baseline={};slow={};maps={'baseline':{},'slow':{}}
 for spec in DEFAULT_V02_SCALES:
  r=DETECTORS[spec.name].analyze(observation,return_maps=True);candidate=slow_maps(r.tf_maps)
  baseline[spec.name]=(r.max_anomaly_score,r.candidates,r.diagnostics)
  slow[spec.name]=cluster_maps(candidate,observation.duration,observation.sample_rate,r.diagnostics)
  maps['baseline'][spec.name]=r.tf_maps;maps['slow'][spec.name]=candidate
 return {'baseline':baseline,'slow':slow,'maps':maps}

class TableDetector(MultiResolutionDetector):
 def __init__(self,table):super().__init__();self.table=table
 def _raw_analyze_all_scales(self,observation):return self.table[observation.sample_id]

def score_table(records,arm):return {r['sample_id']:{k:(v,[],{}) for k,v in r[arm].items()} for r in records}

def focused_checks():
 assert [slow_kernel(n,h) for n,h in [(1017,4/1024),(125,32/1024),(29,128/1024)]]==[257,33,9]
 assert slow_kernel(8,4/1024)==7 and slow_kernel(1,4/1024)==1
 assert np.array_equal(slow_envelope(np.ones(401)*2,4/1024),np.ones(401)*2)
 assert np.all(slow_envelope(np.zeros(401),4/1024)==1e-18)
 step=np.r_[np.ones(200),np.ones(201)*4];assert np.array_equal(slow_envelope(step,4/1024),step)
 pulse=np.ones(401);pulse[185:216]=4
 assert ndimage.median_filter(pulse,size=31,mode='nearest')[200]==4 and slow_envelope(pulse,4/1024)[200]==1
 wide=np.ones(401);wide[100:301]=4;assert slow_envelope(wide,4/1024)[200]==4
 for a,h in [(np.array([]),.1),(np.array([np.nan]),.1),(np.ones(4),0.)]:
  try:slow_envelope(a,h)
  except ValueError:pass
  else:raise AssertionError('Invalid envelope accepted')
 from bsd.detector import TimeFrequencyMaps
 z=np.ones((3,5));z[1,0]=z[1,2]=5;z[0,4]=4
 m=TimeFrequencyMaps(np.array([32.,64.,96.]),np.arange(5)*.1+.05,z,z,np.ones(3),np.ones(5));diag={'input_kurtosis':0.,'input_crest_factor':1.,'total_tf_bins':15}
 assert not cluster_maps(m,.5,256.,diag)[1]
 z[2,1]=5;res=cluster_maps(m,.5,256.,diag);assert len(res[1])==1 and res[1][0].n_bins==3 and res[0]==12
 # Independent fixed-map equation check at every production scale.
 for spec in DEFAULT_V02_SCALES:
  from bsd.generator import Observation
  x=np.random.default_rng(119).normal(size=4096);m=DETECTORS[spec.name].analyze(Observation('fixture',x,1024.,4.),return_maps=True).tf_maps
  raw=np.median(m.power_spectrogram/m.background_psd[:,None],axis=0)/np.log(2)
  expected=np.maximum(ndimage.median_filter(raw,size={'short_31ms':257,'medium_125ms':33,'long_500ms':9}[spec.name],mode='nearest'),1e-18)
  np.testing.assert_array_equal(slow_maps(m).background_time_envelope,expected)
 return {'status':'PASS','kernel_sizes':[257,33,9],'odd_frame_clipping':True,'constant_and_positive_floor':True,'ideal_step_preserved':True,'short_spike_rejected':True,'wide_spike_retained':True,'invalid_inputs_rejected':3,'original_eight_neighbor_two_cell_topology':True,'strict_seed_boundary':True,'fixed_map_equation_checks':3}
