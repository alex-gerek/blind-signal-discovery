"""Fixed approximately-one-second temporal-envelope pilot; no production edits."""
from __future__ import annotations
import csv,hashlib,importlib.util,json,os,platform,subprocess,sys,uuid
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'));sys.dont_write_bytecode=True
PARENT=ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'
PROTOCOL=ROOT/'docs/temporal_envelope_pilot_protocol_2026-11-27.md'
HELPER=ROOT/'experiments/temporal_envelope_prototype.py'
SITES=('H1','L1');ARMS=('baseline','slow')
ALLOWED={str((ROOT/'data/gwosc/GW170814'/f'{s[0]}-{s}_GWOSC_4KHZ_R1-1186739814-4096.hdf5').resolve()).lower().replace('\\','/') for s in SITES}
DATA=str((ROOT/'data').resolve()).lower().replace('\\','/')+'/'
RESULTS=str((ROOT/'results').resolve()).lower().replace('\\','/')+'/'
ACTIVE=None

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
  writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if 'o3_holdout' in p or (p.startswith(DATA) and (p not in ALLOWED or writing)):raise PermissionError('Nonallowlisted data/O3 access blocked')
  if writing and p.startswith(RESULTS) and (ACTIVE is None or not p.startswith(str(ACTIVE).lower().replace('\\','/')+'/')):raise PermissionError('Historical output writes blocked')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
import scipy
from scipy import signal
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bsd.generator import Observation,generate_clean_signal,_scale_to_energy_snr
from bsd.noise import generate_noise
from bsd.detector_v2 import DEFAULT_V02_SCALES
from temporal_envelope_prototype import raw_pair,TableDetector,score_table,focused_checks,DETECTORS
spec=importlib.util.spec_from_file_location('o2_helpers',ROOT/'experiments/run_gw170814_exploratory_search.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arrayhash(x):return hashlib.sha256(np.asarray(x,dtype='<f8',order='C').tobytes()).hexdigest()
def git(*a):return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=encode),encoding='utf-8')
def csvsave(p,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def obs(x,name):return Observation(name,np.asarray(x,dtype=float),1024.,4.)
def rng(*key):return np.random.default_rng(np.random.SeedSequence(key))
def noise(model,generator):return generate_noise(model,4096,1024.,generator,sigma=1.,params={'alpha':1.,'f_knee':16.})[0]
def overlap(c,gate):
 return gate is not None and c.t_start<=gate[1] and c.t_end>=gate[0] and c.f_min<=gate[3] and c.f_max>=gate[2]

def compare_case(dest,calobs,evaluation,real=False):
 dest.mkdir();features=[]
 for i,o in enumerate(calobs):
  pair=raw_pair(o)
  features.append({'sample_id':o.sample_id,'input_sha256':arrayhash(o.signal),**{arm:{k:v[0] for k,v in pair[arm].items()} for arm in ARMS}})
  if (i+1)%200==0:print('CAL_FEATURES',dest.name,i+1,flush=True)
 save(dest/'calibration_raw_features.json',features)
 cals={}
 for arm in ARMS:
  detector=TableDetector(score_table(features,arm));cal=detector.calibrate(calobs,target_fap=.01,norm_fraction=.30,calibration_seed=0);cals[arm]=cal
  save(dest/f'{arm}_calibration.json',cal.to_summary_dict());save(dest/f'{arm}_null_scores.json',cal.sorted_global_null_scores)
 rows=[];candidate_records=[];eventraw={};tracking=[];paired=[];noise_cache={}
 for item in evaluation:
  o=item['observation'];pair=raw_pair(o)
  row={k:v for k,v in item.items() if k not in ('observation','gate','profile','noise_partner')};row['input_sha256']=arrayhash(o.signal);row['evaluation_gate']=item['gate']
  if 'profile' in item:
   for arm in ARMS:
    for name,maps in pair['maps'][arm].items():
     truth=np.interp(maps.times,np.arange(4096)/1024.,item['profile']**2);truth=truth/np.median(truth)
     est=maps.background_time_envelope/np.median(maps.background_time_envelope)
     tracking.append({'sample_id':o.sample_id,'regime':item['regime'],'arm':arm,'scale':name,'log_rmse_proxy':float(np.sqrt(np.mean((np.log(est)-np.log(truth))**2)))})
  if 'noise_partner' in item:
   partner=item['noise_partner']
   if partner.sample_id not in noise_cache:noise_cache[partner.sample_id]=raw_pair(partner)
   for arm in ARMS:
    for name,maps in pair['maps'][arm].items():
     cleanmaps=noise_cache[partner.sample_id]['maps'][arm][name];active=(maps.times>=item['gate'][0])&(maps.times<=item['gate'][1]);assert active.any()
     paired.append({'sample_id':o.sample_id,'regime':item['regime'],'family':item['family'],'snr':item['snr'],'arm':arm,'scale':name,'temporal_ratio':float(np.median((maps.background_time_envelope/cleanmaps.background_time_envelope)[active]))})
  for arm in ARMS:
   cal=cals[arm];result=TableDetector({o.sample_id:pair[arm]}).analyze(o,cal)
   veto=old.apply_vetoes(o,result.candidates,old.FROZEN_VETO_CONFIG) if real else None
   accepted_near=[c for c in result.candidates if overlap(c,item['gate'])]
   nearstd=[cal.scale_stats[name].standardize(c.anomaly_score) for name,(_,cc,_) in pair[arm].items() for c in cc if overlap(c,item['gate'])]
   nearmax=max(nearstd,default=None)
   row[arm]={'score':result.max_global_score,'threshold':cal.global_threshold,'empirical_p':cal.empirical_p_value(result.max_global_score)[0],'detected_window':result.detected,'associated_raw':bool(accepted_near),'near_score':nearmax,'near_p':cal.empirical_p_value(nearmax)[0] if nearmax is not None else None,'post_veto_window':veto.detected_after_veto if veto else None,'associated_post_veto':any(not v.vetoed and overlap(v.candidate,item['gate']) for v in veto.decisions) if veto else None,'per_scale_raw_max':result.per_scale_raw_max}
   vm={v.candidate.candidate_id:v for v in veto.decisions} if veto else {}
   for c in result.candidates:
    v=vm.get(c.candidate_id);candidate_records.append({'sample_id':o.sample_id,'arm':arm,'role':item['role'],**c.to_dict(),'vetoed':v.vetoed if v else None,'veto_reasons':v.reason_codes if v else []})
   if item['role']=='ON_SOURCE':eventraw[arm]={name:[c.to_dict() for c in cc] for name,(_,cc,_) in pair[arm].items()}
  rows.append(row)
 save(dest/'evaluation.json',rows);save(dest/'candidates.json',candidate_records)
 save(dest/'tracking_proxy.json',tracking);save(dest/'paired_normalizer_response.json',paired)
 if real:save(dest/'event_raw_clusters.json',eventraw)
 summary={'case':dest.name,'n_calibration':len(calobs),'n_evaluation':len(rows),'arms':{},'discordant':{}}
 for arm in ARMS:
  summary['arms'][arm]={}
  for role in sorted({r['role'] for r in rows}):
   rr=[r for r in rows if r['role']==role]
   summary['arms'][arm][role]={'windows':len(rr),'window_nominations':sum(r[arm]['detected_window'] for r in rr),'associated_nominations':sum(r[arm]['associated_raw'] for r in rr),'post_veto_window_nominations':sum(r[arm]['post_veto_window'] for r in rr) if real else None,'associated_post_veto':sum(r[arm]['associated_post_veto'] for r in rr) if real else None}
 for role in sorted({r['role'] for r in rows}):
  rr=[r for r in rows if r['role']==role];key='associated_raw' if role.startswith('SIGNAL') or role=='ON_SOURCE' else 'detected_window'
  summary['discordant'][role]={'slow_only':sum(r['slow'][key] and not r['baseline'][key] for r in rr),'baseline_only':sum(r['baseline'][key] and not r['slow'][key] for r in rr)}
 summary['by_regime']={}
 for regime in sorted({r.get('regime','real') for r in rows}):
  selected=[r for r in rows if r.get('regime','real')==regime];d={'arms':{},'discordant':{}}
  for arm in ARMS:
   d['arms'][arm]={}
   for role in sorted({r['role'] for r in selected}):
    rr=[r for r in selected if r['role']==role]
    d['arms'][arm][role]={'windows':len(rr),'window_nominations':sum(r[arm]['detected_window'] for r in rr),'associated_nominations':sum(r[arm]['associated_raw'] for r in rr)}
  for role in sorted({r['role'] for r in selected}):
   rr=[r for r in selected if r['role']==role];key='associated_raw' if role.startswith('SIGNAL') or role=='ON_SOURCE' else 'detected_window'
   d['discordant'][role]={'slow_only':sum(r['slow'][key] and not r['baseline'][key] for r in rr),'baseline_only':sum(r['baseline'][key] and not r['slow'][key] for r in rr)}
  summary['by_regime'][regime]=d
 save(dest/'summary.json',summary)
 return summary,rows,cals

def run_real(out):
 results=[]
 for site in SITES:
  print('REAL_START',site,flush=True)
  expected=load(PARENT/site/'input_provenance.json')['original'];path=old.DETECTORS[site];assert digest(path)==expected['sha256']
  record=old.load_gwosc_hdf5(path);inj,names=old._load_inj_mask(path)
  sequences=[]
  for interval,n in [(old.PSD_CORE_INTERVAL,60),(old.CAL_A_CORE_INTERVAL,320),(old.CAL_B_CORE_INTERVAL,387),(old.SEARCH_CORE_INTERVAL,160)]:
   segs,bad=old._extract_exact_interval_segments(record,interval);assert not bad and len(segs)==n;sequences.append(segs)
  model=old.fit_background_model(site,sequences[0],mode=old.PREPROCESS_MODE,target_sample_rate=1024.,bandpass_hz=old.BANDPASS_HZ)
  calobs=[old.preprocess_segment(s,f'{site}_cal_{i}',model) for i,s in enumerate(sequences[1]+sequences[2])]
  evaluation=[]
  for i,s in enumerate(sequences[3]):
   gps=s.gps_start;on=gps<=old.EVENT_GPS<gps+4
   role='ON_SOURCE' if on else ('OFF_SOURCE' if gps+5<=old.EVENT_GPS-16 or gps-1>=old.EVENT_GPS+16 else 'GUARD')
   evaluation.append({'sample_id':f'{site}_target_{i}','window_index':i,'gps_start':gps,'role':role,'gate':[old.EVENT_GPS-gps-.5,old.EVENT_GPS-gps+.5,35.,250.] if on else None,'observation':old.preprocess_segment(s,f'{site}_target_{i}',model),'quality':old._inspect_flags_at_interval(record,inj,names,gps,gps+4)})
  dest=out/site;summary,rows,cals=compare_case(dest,calobs,evaluation,True)
  save(dest/'input_provenance.json',{'original':expected,'model':model.__dict__})
  par=load(PARENT/site/'summary.json');event=next(r for r in rows if r['role']=='ON_SOURCE')
  for k in ('score','threshold','empirical_p'):assert abs(event['baseline'][k]-par['event'][k])<1e-10
  assert summary['arms']['baseline']['OFF_SOURCE']['windows']==151
  assert summary['arms']['baseline']['OFF_SOURCE']['window_nominations']==par['offsource_raw']
  assert summary['arms']['baseline']['OFF_SOURCE']['post_veto_window_nominations']==par['offsource_post_veto']
  results.append(summary);print('REAL_DONE',site,json.dumps(summary['arms']),flush=True)
 return results


REGIMES=('stationary','smooth_drift','persistent_step')
FAMILIES=('linear_chirp','band_limited_noise_burst')

def profiled(base,seed,kind,i,regime_index):
 t=np.arange(4096)/1024.;g=rng(seed,4,kind,i,regime_index);params={'regime':REGIMES[regime_index]}
 if regime_index==0:envelope=np.ones(4096);scale=1.;x=base.copy()
 elif regime_index==1:
  phase=float(g.uniform(0,2*np.pi));envelope=1+.5*np.sin(2*np.pi*.25*t+phase);x=base*envelope;scale=float(1/np.sqrt(np.mean(x*x)));x=x*scale;params.update({'phase_rad':phase,'modulation_depth':.5,'frequency_hz':.25})
 else:
  when=float(g.uniform(1,3));factor=2. if g.random()<.5 else .5;envelope=np.where(t<when,1.,factor);x=base*envelope;scale=float(1/np.sqrt(np.mean(x*x)));x=x*scale;params.update({'step_time_s':when,'step_factor':factor})
 params.update({'global_scale':scale,'envelope_sha256':arrayhash(envelope)})
 return x,envelope,params

def clean_signal(seed,i,family_index):
 g=rng(seed,3,i,family_index)
 if family_index==0:
  params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
  clean,a,b,lo,hi,metadata=generate_clean_signal('linear_chirp',4096,1024.,g,params)
 else:
  a,b,lo,hi=1.375,1.625,35.,350.;t=np.arange(4096)/1024.;active=(t>=a)&(t<=b)
  full=signal.sosfiltfilt(signal.butter(4,[35,350],fs=1024.,btype='bandpass',output='sos'),g.normal(size=4096))
  clean=np.zeros(4096);clean[active]=full[active]*tukey(int(active.sum()),alpha=.15)
  metadata={'t_start_s':a,'t_end_s':b,'filter':'butter4_sosfiltfilt_35_350Hz','taper_alpha':.15}
 return clean,[a,b,lo,hi],metadata

def run_synthetic(out):
 results=[]
 for seed in (9101,9102,9103):
  name=f'synthetic_{seed}';print('SYNTHETIC_START',seed,flush=True)
  calobs=[obs(noise('pink_colored',rng(seed,0,i)),f'{name}_cal_{i}') for i in range(320)]
  evaluation=[]
  for regime_index,regime in enumerate(REGIMES):
   for i in range(200):
    base=noise('pink_colored',rng(seed,1,i));x,envelope,params=profiled(base,seed,0,i,regime_index);sample=f'{name}_{regime}_noise_{i}'
    evaluation.append({'sample_id':sample,'role':'NOISE','regime':regime,'noise_index':i,'base_noise_sha256':arrayhash(base),'profile_params':params,'profile':envelope,'gate':None,'observation':obs(x,sample)})
   for family_index,family in enumerate(FAMILIES):
    for i in range(10):
     base=noise('pink_colored',rng(seed,2,i,family_index));x,envelope,params=profiled(base,seed,1+family_index,i,regime_index);clean,gate,metadata=clean_signal(seed,i,family_index)
     partner=obs(x,f'{name}_{regime}_{family}_noise_partner_{i}')
     for snr in (8,18):
      added=_scale_to_energy_snr(clean,float(np.std(x)),snr);mixture=x+added;sample=f'{name}_{regime}_{family}_inj_{i}_snr{snr}'
      np.testing.assert_allclose(np.linalg.norm(added)/np.std(x),snr,rtol=1e-12)
      evaluation.append({'sample_id':sample,'role':f'SIGNAL_{family}_SNR{snr}','family':family,'regime':regime,'injection_id':i,'snr':snr,'base_noise_sha256':arrayhash(base),'noise_sha256':arrayhash(x),'clean_unit_sha256':arrayhash(clean),'signal_sha256':arrayhash(added),'profile_params':params,'generation_params':metadata,'gate':gate,'noise_partner':partner,'observation':obs(mixture,sample)})
  summary,_,_=compare_case(out/name,calobs,evaluation,False);results.append(summary)
  print('SYNTHETIC_DONE',seed,json.dumps(summary['by_regime']),flush=True)
 return results

def make_figures(out):
 summaries=load(out/'summaries.json');synth=[s for s in summaries if s['case'].startswith('synthetic')]
 fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
 for ax,regime in zip(axes,REGIMES):
  labels=('NOISE','chirp8','burst8','chirp18','burst18');keys=('NOISE','SIGNAL_linear_chirp_SNR8','SIGNAL_band_limited_noise_burst_SNR8','SIGNAL_linear_chirp_SNR18','SIGNAL_band_limited_noise_burst_SNR18');x=np.arange(5)
  for offset,arm in [(-.18,'baseline'),(.18,'slow')]:
   vals=[]
   for key in keys:
    items=[s['by_regime'][regime]['arms'][arm][key] for s in synth];num=sum(r['window_nominations' if key=='NOISE' else 'associated_nominations'] for r in items);den=sum(r['windows'] for r in items);vals.append(100*num/den)
   ax.bar(x+offset,vals,.36,label=arm)
  ax.set_xticks(x,labels,rotation=35);ax.set_ylim(0,105);ax.set_title(regime);ax.set_ylabel('Observed percentage (different denominators)');ax.grid(axis='y',alpha=.2)
 axes[-1].legend();fig.suptitle('Descriptive pilot totals; noise denominator600, each signal group30; no common-rate CI')
 fig.savefig(out/'synthetic_comparison.png',dpi=150);plt.close(fig)
 tracking=[];paired=[]
 for seed in (9101,9102,9103):tracking.extend(load(out/f'synthetic_{seed}'/'tracking_proxy.json'));paired.extend(load(out/f'synthetic_{seed}'/'paired_normalizer_response.json'))
 groups=[]
 for regime in REGIMES:
  for arm in ARMS:
   for spec in DEFAULT_V02_SCALES:
    rr=[r['log_rmse_proxy'] for r in tracking if r['regime']==regime and r['arm']==arm and r['scale']==spec.name]
    groups.append({'regime':regime,'arm':arm,'scale':spec.name,'negative_records':len(rr),'median_log_rmse_proxy':float(np.median(rr)),'min':min(rr),'max':max(rr)})
 save(out/'tracking_summary.json',groups)
 response=[]
 for regime in REGIMES:
  for family in FAMILIES:
   for snr in (8,18):
    for arm in ARMS:
     for spec in DEFAULT_V02_SCALES:
      rr=[r['temporal_ratio'] for r in paired if r['regime']==regime and r['family']==family and r['snr']==snr and r['arm']==arm and r['scale']==spec.name]
      response.append({'regime':regime,'family':family,'snr':snr,'arm':arm,'scale':spec.name,'paired_records':len(rr),'median_temporal_ratio':float(np.median(rr)),'min':min(rr),'max':max(rr),'T_increased':sum(v>1 for v in rr),'T_decreased':sum(v<1 for v in rr),'T_unchanged':sum(v==1 for v in rr)})
 save(out/'paired_response_summary.json',response)
 fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
 for ax,regime in zip(axes,REGIMES):
  for arm in ARMS:
   rr=[next(g['median_log_rmse_proxy'] for g in groups if g['regime']==regime and g['arm']==arm and g['scale']==s.name) for s in DEFAULT_V02_SCALES]
   ax.plot(['31ms','125ms','500ms'],rr,'o-',label=arm)
  ax.set_title(regime);ax.set_ylabel('Median relative-envelope log-RMSE proxy');ax.grid(alpha=.2)
 axes[-1].legend();fig.suptitle('Tracking descriptor against known input envelope; not unbiased STFT variance')
 fig.savefig(out/'tracking_proxy.png',dpi=150);plt.close(fig)


def main():
 global ACTIVE
 pm=load(PARENT/'manifest.json');sources=pm['source_sha256'];assert all(digest(ROOT/f)==h for f,h in sources.items())
 assert not git('diff','HEAD','--','src','app','docs/workstation/index.html')
 rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
 out=ROOT/'results/temporal_envelope_pilot'/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 for name,p in [('driver.py',Path(__file__)),('temporal_envelope_prototype.py',HELPER),('protocol.md',PROTOCOL)]:(out/name).write_bytes(p.read_bytes())
 parents={str((PARENT/s/'summary.json').relative_to(ROOT)).replace('\\','/'):digest(PARENT/s/'summary.json') for s in SITES}
 m={'run_id':rid,'source_commit':git('rev-parse','HEAD'),'source_sha256':sources,'parent_sha256':parents,'driver_sha256':digest(Path(__file__)),'helper_sha256':digest(HELPER),'protocol_sha256':digest(PROTOCOL),'dirty_paths':git('status','--porcelain').splitlines(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'status':'RUNNING','stages':{}}
 save(out/'manifest.json',m);print('RUN_DIRECTORY',out,flush=True)
 fixtures=focused_checks();save(out/'focused_checks.json',fixtures)
 real=run_real(out);m['stages']['real']='COMPLETE';save(out/'manifest.json',m)
 synthetic=run_synthetic(out);m['stages']['synthetic']='COMPLETE';save(out/'summaries.json',real+synthetic);save(out/'manifest.json',m)
 make_figures(out)
 assert all(digest(ROOT/f)==h for f,h in sources.items()) and all(digest(ROOT/f)==h for f,h in parents.items()) and digest(HELPER)==m['helper_sha256']
 m['status']='COMPLETE';m['source_and_parent_unchanged']=True;save(out/'manifest.json',m);print('COMPLETE',rid,flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  if ACTIVE is not None and (ACTIVE/'manifest.json').exists():
   m=load(ACTIVE/'manifest.json');m['status']='ERROR';m['error']={'type':type(exc).__name__,'message':str(exc)};save(ACTIVE/'manifest.json',m)
  raise
