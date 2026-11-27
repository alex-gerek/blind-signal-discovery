"""Frozen blind two-tone pilot; research-only, no production mutation."""
from pathlib import Path
from dataclasses import asdict,replace
from datetime import datetime,timezone
from itertools import product,combinations
import hashlib,json,os,platform,subprocess,sys,time,traceback,uuid,zipfile
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True
import numpy as np
import scipy
from scipy import signal,ndimage
from bsd.generator import Observation
from bsd.noise import generate_pink_colored
from bsd.detector import TemplateFreeDetector
from bsd.detector_v2 import MultiResolutionDetector,DEFAULT_V02_SCALES,ScaleSpec
PROTOCOL=ROOT/'docs/blind_two_tone_exploratory_frozen_2026-11-27.md'

def norm(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,default=norm,allow_nan=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(x):return hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
def scales(w):return DEFAULT_V02_SCALES if w==0 else DEFAULT_V02_SCALES[:2]+(ScaleSpec('long_1000ms',replace(DEFAULT_V02_SCALES[-1].detector_config,nperseg=1024,noverlap=896)),)
def make(seed,index,kind,params=None):
 r=np.random.default_rng(np.random.SeedSequence([seed,index]));noise=generate_pink_colored(4096,1024.,r,sigma=1.,alpha=1.,f_knee=16.)[0]
 p=dict(params or {});base=float(2*r.integers(40,90)+p.get('offset',float(r.choice([0,.5]))));onset=float(p.get('onset',r.choice([1.25,1.25+63/1024])));snr=float(p.get('snr',r.choice([8.,18.,36.])));phase=float(r.uniform(0,2*np.pi));gap=float(p.get('gap',0));ratio=float(p.get('ratio',0));delta=float(p.get('phase',0))
 # A one-second waveform is generated on the full integer timeline before record-edge clipping.
 start=int(round(onset*1024));idx=np.arange(start,start+1024);local=np.arange(1024)/1024.;env=signal.windows.tukey(1024,.15)
 a=env*np.sin(2*np.pi*base*local+phase);b=ratio*env*np.sin(2*np.pi*(base+gap)*local+phase+delta) if kind=='two' else np.zeros(1024)
 raw=a+b;energy=float(np.dot(raw,raw));assert energy>1e-10
 scale=snr/np.sqrt(energy);a*=scale;b*=scale;addition=np.zeros(4096);valid=(idx>=0)&(idx<4096);addition[idx[valid]]=(a+b)[valid]
 if kind=='noise':addition[:]=0
 x=noise+addition
 meta={'seed':seed,'index':index,'kind':kind,'f1':base if kind!='noise' else None,'f2':base+gap if kind=='two' else None,'gap':gap if kind=='two' else None,'ratio':ratio if kind=='two' else None,'phase':delta if kind=='two' else None,'global_phase':phase,'offset':base%2,'onset':onset,'observed_start':max(0,onset),'observed_end':min(4,onset+1),'snr':snr if kind!='noise' else 0,'observed_snr':float(np.linalg.norm(addition)),'e1':float(a@a) if kind!='noise' else 0.,'e2':float(b@b) if kind!='noise' else 0.,'cross':float(2*a@b) if kind!='noise' else 0.,'hash':key(x),'grid_remainder_2hz':[base%2,(base+gap)%2] if kind=='two' else [base%2],'grid_remainder_1hz':[base%1,(base+gap)%1] if kind=='two' else [base%1]}
 return x,meta

def proposals(candidates,z,freqs,times,labels,within):
 proposals=[]
 def add(a,b,h1,h2,start,end,peak,origin):
  if a==b:return
  f=sorted([float(a),float(b)]);proposals.append({'frequencies':f,'score':float(min(h1,h2)),'start':float(start),'end':float(end),'peak':float(peak),'origin':origin})
 for a,b in combinations(candidates,2):
  if (a['f_max']<b['f_min'] or b['f_max']<a['f_min']) and max(a['t_start'],b['t_start'])<=min(a['t_end'],b['t_end']):add(a['peak_freq'],b['peak_freq'],a['peak_excess_power'],b['peak_excess_power'],min(a['t_start'],b['t_start']),max(a['t_end'],b['t_end']),(a['peak_time']+b['peak_time'])/2,'disjoint')
 profiles=[]
 if within:
  for c in candidates:
   fi=int(np.argmin(abs(freqs-c['peak_freq'])));ti=int(np.argmin(abs(times-c['peak_time'])));lab=int(labels[fi,ti]);assert lab>0
   cols=np.unique(np.where(labels==lab)[1]);profile=np.mean(z[:,cols],axis=1);peaks,_=signal.find_peaks(profile)
   peaks=[int(i) for i in peaks if c['f_min']<=freqs[i]<=c['f_max'] and profile[i]>=4]
   profiles.append({'cluster_peak':[c['peak_freq'],c['peak_time']],'columns':cols.tolist(),'profile':profile.tolist(),'peaks':peaks})
   for i,j in combinations(peaks,2):
    if j-i>=2 and np.min(profile[i+1:j])<=.8*min(profile[i],profile[j]):add(freqs[i],freqs[j],profile[i],profile[j],c['t_start'],c['t_end'],c['peak_time'],'within')
 unique={}
 for p in sorted(proposals,key=lambda p:(-p['score'],p['frequencies'],p['start'])):unique.setdefault(tuple(p['frequencies']),p)
 ordered=sorted(unique.values(),key=lambda p:(-p['score'],p['frequencies'],p['start']))
 return ordered,profiles

def extract(x,w,cal):
 t0=time.perf_counter();raw={};accepted=[];long=None;stages={}
 for spec in scales(w):
  t=time.perf_counter();r=TemplateFreeDetector(spec.detector_config).analyze(Observation('blind',x,1024.,4.),return_maps=spec==scales(w)[-1]);stages[spec.name]=time.perf_counter()-t
  stat=cal.scale_stats[spec.name];raw[spec.name]={'max':r.max_anomaly_score,'candidates':[asdict(c) for c in r.candidates]}
  for c in r.candidates:
   score=stat.standardize(c.anomaly_score)
   if score>cal.global_threshold:accepted.append({**asdict(c),'scale':spec.name,'standardized':score})
  if spec==scales(w)[-1]:long=r
 banktime=time.perf_counter()-t0;score=max(cal.scale_stats[n].standardize(v['max']) for n,v in raw.items());z=long.tf_maps.normalized_excess_map;labels,n=ndimage.label(z>4,structure=np.ones((3,3),int));cands=[c for c in accepted if c['scale']==scales(w)[-1].name]
 methods={}
 for d in [0,1]:
  t=time.perf_counter();p,profiles=proposals(cands,z,long.tf_maps.freqs,long.tf_maps.times,labels,bool(d));methods[d]={'proposals':p,'profiles':profiles,'seconds':time.perf_counter()-t}
 return {'raw':raw,'accepted':accepted,'anomaly':score>cal.global_threshold,'score':score,'per_scale_seconds':stages,'bank_seconds':banktime,'methods':methods},long.tf_maps,labels

def cutoff(values):
 a=sorted(-np.inf if v is None else v for v in values);v=a[int(np.ceil((len(a)+1)*.99))-1];return None if np.isneginf(v) else float(v)
def claims(p,cut):return bool(p and p['score']>(-np.inf if cut is None else cut))
def match(fs,truth,tol):
 if fs is None or len(fs)!=2 or fs[0]==fs[1]:return False,False,None
 errors=[[abs(fs[0]-truth[0]),abs(fs[1]-truth[1])],[abs(fs[0]-truth[1]),abs(fs[1]-truth[0])]];good=[max(e)<=tol for e in errors];best=min(errors,key=lambda e:(max(e),sum(e)))
 return any(good),all(good),best

def fixtures():
 assert match([96,96],[96,98],1)[0] is False
 assert match([98,96],[96,98],1)[0] and match([95,99],[96,98],1)[0]
 assert not match([94,98],[96,98],1)[0] and match([97,98],[96,98],1)[1]
 assert not claims(None,None) and claims({'score':4},None) and not claims({'score':4},4)
 z=np.array([[1,1],[6,6],[1,1],[8,8],[1,1]],float);freq=np.arange(5.);ts=np.array([1.,2.]);labs=np.ones_like(z,int);c={'peak_freq':3.,'peak_time':1.,'f_min':0.,'f_max':4.,'t_start':.5,'t_end':2.5,'peak_excess_power':8.}
 p,_=proposals([c],z,freq,ts,labs,True);assert p[0]['frequencies']==[1.,3.] and not proposals([c],z,freq,ts,labs,False)[0]
 return {'status':'PASS','checks':['duplicate estimates','permutation matching','inclusive tolerance','ambiguous matches','strict score','no proposal','within-cluster double peak','D0 no pair']}

def main():
 global OUT
 assert git('rev-parse','HEAD')=='cc68b93fee699307e3314fce1bd4c9fe403c1e6e';assert not git('diff','--','src','app','pyproject.toml')
 OUT=ROOT/'results/blind_two_tone_exploratory'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);OUT.mkdir(parents=True,exist_ok=False);print('RUN_DIRECTORY',OUT,flush=True)
 (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes());(OUT/'driver.py').write_bytes(Path(__file__).read_bytes());sources={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in (ROOT/'src/bsd').glob('*.py')}
 historical={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for folder in ['two_oscillation_anomaly','two_oscillation_structure_review','two_oscillation_precluster_inspection'] for p in (ROOT/'results'/folder).rglob('*') if p.is_file()}
 manifest={'status':'RUNNING','source':git('rev-parse','HEAD'),'source_sha256':sources,'historical_sha256':historical,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'dirty':git('status','--short'),'driver_sha256':sha(Path(__file__)),'protocol_sha256':sha(PROTOCOL),'runtime_threads':{n:os.environ[n] for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}}
 save(OUT/'manifest.json',manifest);save(OUT/'fixtures.json',fixtures())
 with zipfile.ZipFile(OUT/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
  for name in sources:archive.write(ROOT/name,name)
 # Tiny smoke: development input, no held-out data or method tuning.
 x,_=make(30101,0,'two',{'gap':8,'ratio':1,'phase':0,'snr':18,'offset':0,'onset':1.25})
 for w in [0,1]:assert TemplateFreeDetector(scales(w)[-1].detector_config).analyze(Observation('smoke',x,1024.,4.),return_maps=True).tf_maps is not None
 print('SMOKE PASS',flush=True)
 models={};calibrations={};nulls={};cal_records=[]
 for w in [0,1]:
  class Logged(MultiResolutionDetector):
   def _raw_analyze_all_scales(self,o):
    result=super()._raw_analyze_all_scales(o);cal_records.append({'w':w,'id':o.sample_id,'hash':key(o.signal),'raw':{n:{'max':v[0],'candidates':[asdict(c) for c in v[1]]} for n,v in result.items()}});return result
  model=Logged(scales(w));obs=[Observation(f'cal_{i}',make(30201,i,'noise')[0],1024.,4.) for i in range(284)];start=time.perf_counter();cal=model.calibrate(obs,target_fap=.005,norm_fraction=.3,calibration_seed=30201);elapsed=time.perf_counter()-start
  calibrations[w]=cal;models[w]=MultiResolutionDetector(scales(w));save(OUT/f'anomaly_calibration_W{w}.json',{'calibration':asdict(cal),'seconds':elapsed});print('ANOMALY_CALIBRATED',w,cal.global_threshold,elapsed,flush=True)
 save(OUT/'anomaly_calibration_raw.json',cal_records)
 pair_rows=[];cuts={};parity=[]
 for w in [0,1]:
  nulls[w]={d:{'noise':[],'one':[]} for d in [0,1]}
  for kind,seed in [('noise',30301),('one',30401)]:
   for i in range(199):
    x,meta=make(seed,i,kind);result,_,_=extract(x,w,calibrations[w]);pair_rows.append({'w':w,'meta':meta,'result':result})
    for d in [0,1]:p=result['methods'][d]['proposals'];nulls[w][d][kind].append(p[0]['score'] if p else None)
    if i==0:
     direct=models[w].analyze(Observation('parity',x,1024.,4.),calibrations[w]);assert direct.detected==result['anomaly'] and np.isclose(direct.max_global_score,result['score']) and len(direct.candidates)==len(result['accepted']);parity.append([w,kind,seed,i])
  for d in [0,1]:
   by={kind:cutoff(v) for kind,v in nulls[w][d].items()};finite=[x for x in by.values() if x is not None];cuts[f'W{w}D{d}']={'by_group':by,'cutoff':max(finite) if finite else None,'negative_infinity_is_no_proposal':True}
  print('PAIR_CALIBRATED',w,{k:v for k,v in cuts.items() if k.startswith(f'W{w}')},flush=True)
 save(OUT/'pair_calibration_records.json',pair_rows);save(OUT/'frozen_pair_thresholds.json',{'cuts':cuts,'null_scores':nulls,'heldout_evaluated':0});frozen=sha(OUT/'frozen_pair_thresholds.json');manifest['threshold_freeze_sha256']=frozen;save(OUT/'manifest.json',manifest)
 jobs=[('noise',30501,i,{},'main') for i in range(160)]+[('one',30601,i,{},'main') for i in range(160)]
 matrix=list(product([1,2,4,8,16,96],[0.,np.pi/2,np.pi],[1.,.5,.25],[0.,.5],[1.25,1.25+63/1024],[8.,18.]))
 assert len(matrix)==432
 for i,(gap,phase,ratio,offset,onset,snr) in enumerate(matrix):jobs.append(('two',30701,i,dict(gap=gap,phase=phase,ratio=ratio,offset=offset,onset=onset,snr=snr),'main'))
 for i,(onset,gap,ratio,phase) in enumerate(product([-.25,3.5],[2,8],[1.,.25],[0.,np.pi/2,np.pi])):jobs.append(('two',30801,i,dict(onset=onset,gap=gap,ratio=ratio,phase=phase,offset=0.,snr=18.),'edge'))
 allrows=[];hashes={};trace_dir=OUT/'traces';trace_dir.mkdir();tstart=time.perf_counter()
 for ji,(kind,seed,index,params,part) in enumerate(jobs):
  x,meta=make(seed,index,kind,params);ident=f'{seed}_{index:04d}';hashes[ident]=meta['hash']
  for w in [0,1]:
   result,maps,labels=extract(x,w,calibrations[w]);t=time.perf_counter()
   # Additional phase-aware trace calculation is artifact overhead, not inference runtime.
   f,ts,stft=signal.stft(x-np.median(x),fs=1024.,window='hann',nperseg=scales(w)[-1].detector_config.nperseg,noverlap=scales(w)[-1].detector_config.noverlap,boundary=None,padded=False)
   np.testing.assert_allclose(abs(stft[1:-1])**2,maps.power_spectrogram,rtol=1e-12,atol=1e-15)
   path=trace_dir/f'{ident}_W{w}.npz';np.savez_compressed(path,waveform=x,stft=stft[1:-1],power=maps.power_spectrogram,background_psd=maps.background_psd,background_time=maps.background_time_envelope,z=maps.normalized_excess_map,freqs=maps.freqs,times=maps.times,mask=maps.normalized_excess_map>4,labels=labels);trace_s=time.perf_counter()-t
   decisions={}
   for d in [0,1]:
    arm=f'W{w}D{d}';ps=result['methods'][d]['proposals'];p=ps[0] if ps else None;claim=claims(p,cuts[arm]['cutoff']);freqs=p['frequencies'] if claim else None;correct=ambiguous=False;errors=None;strict=False
    if kind=='two':correct,ambiguous,errors=match(freqs,[meta['f1'],meta['f2']],1.);strict=match(freqs,[meta['f1'],meta['f2']],meta['gap']/4)[0]
    timing=None
    if claim and kind!='noise':timing={'onset_error':p['start']-meta['observed_start'],'end_error':p['end']-meta['observed_end'],'midpoint_error':p['peak']-(meta['observed_start']+meta['observed_end'])/2}
    decisions[arm]={'claim':claim,'frequencies':freqs,'correct':correct,'strict_correct':strict,'ambiguous_match':ambiguous,'frequency_errors':errors,'false_pair':claim and kind!='two','wrong_pair':claim and kind=='two' and not correct,'timing':timing,'proposal':p}
   if ji in [0,160,320]:
    direct=models[w].analyze(Observation('parity',x,1024.,4.),calibrations[w]);assert direct.detected==result['anomaly'] and direct.max_global_score==result['score'] and len(direct.candidates)==len(result['accepted']);parity.append([w,kind,seed,index])
   allrows.append({'id':ident,'part':part,'meta':meta,'w':w,'result':result,'decisions':decisions,'trace':str(path.relative_to(OUT)),'trace_seconds':trace_s})
  if ji%80==0:save(OUT/'progress.json',{'completed_observations':ji+1,'total':len(jobs),'seconds':time.perf_counter()-tstart});print('EVAL_PROGRESS',ji+1,len(jobs),flush=True)
 save(OUT/'evaluation.json',allrows);save(OUT/'parity.json',{'status':'PASS','checks':parity})
 assert frozen==sha(OUT/'frozen_pair_thresholds.json')
 split_hashes={'anomaly':{key(o.signal) for o in obs},'pair':{r['meta']['hash'] for r in pair_rows},'evaluation':set(hashes.values())}
 assert all(not(split_hashes[a]&split_hashes[b]) for a,b in combinations(split_hashes,2))
 assert all(sha(ROOT/p)==h for p,h in {**sources,**historical}.items())
 save(OUT/'integrity.json',{'status':'PASS','no_cross_split_hash_reuse':True,'unique_evaluation_inputs':len(hashes),'source_and_historical_hashes_unchanged':True,'thresholds_unchanged':True})
 manifest.update(status='COMPLETE',evaluation_observations=len(jobs),evaluation_records=len(allrows),evaluation_elapsed_seconds=time.perf_counter()-tstart)
 try:
  import ctypes
  class PMC(ctypes.Structure):_fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong),('PeakWorkingSetSize',ctypes.c_size_t),('WorkingSetSize',ctypes.c_size_t),('QuotaPeakPagedPoolUsage',ctypes.c_size_t),('QuotaPagedPoolUsage',ctypes.c_size_t),('QuotaPeakNonPagedPoolUsage',ctypes.c_size_t),('QuotaNonPagedPoolUsage',ctypes.c_size_t),('PagefileUsage',ctypes.c_size_t),('PeakPagefileUsage',ctypes.c_size_t)]
  pm=PMC();pm.cb=ctypes.sizeof(pm);ok=ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(),ctypes.byref(pm),pm.cb);manifest['peak_working_set_bytes']=int(pm.PeakWorkingSetSize) if ok else None
 except Exception as exc:manifest['memory_metric_unavailable']=str(exc)
 save(OUT/'manifest.json',manifest);print('COMPLETE',OUT,flush=True)
if __name__=='__main__':
 try:main()
 except Exception:
  if 'OUT' in globals():save(OUT/'failure.json',{'status':'FAILED','traceback':traceback.format_exc()})
  raise
