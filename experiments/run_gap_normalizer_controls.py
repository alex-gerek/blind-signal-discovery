"""Bounded paired gap-link calibration and normalizer controls; no production edits."""
from __future__ import annotations
import csv,hashlib,importlib.util,json,os,platform,subprocess,sys,uuid
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'));sys.dont_write_bytecode=True
PARENT=ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'
PROTOCOL=ROOT/'docs/gap_normalizer_control_protocol_2026-11-27.md'
HELPER=ROOT/'experiments/gap_link_prototype.py'
SITES=('H1','L1');ARMS=('baseline','gap')
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
from gap_link_prototype import raw_pair,TableDetector,score_table,focused_checks,DETECTORS
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
 rows=[];candidate_records=[];eventraw={}
 for item in evaluation:
  o=item['observation'];pair=raw_pair(o)
  row={k:v for k,v in item.items() if k not in ('observation','gate')};row['input_sha256']=arrayhash(o.signal)
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
 if real:save(dest/'event_raw_clusters.json',eventraw)
 summary={'case':dest.name,'n_calibration':len(calobs),'n_evaluation':len(rows),'arms':{},'discordant':{}}
 for arm in ARMS:
  summary['arms'][arm]={}
  for role in sorted({r['role'] for r in rows}):
   rr=[r for r in rows if r['role']==role]
   summary['arms'][arm][role]={'windows':len(rr),'window_nominations':sum(r[arm]['detected_window'] for r in rr),'associated_nominations':sum(r[arm]['associated_raw'] for r in rr),'post_veto_window_nominations':sum(r[arm]['post_veto_window'] for r in rr) if real else None,'associated_post_veto':sum(r[arm]['associated_post_veto'] for r in rr) if real else None}
 for role in sorted({r['role'] for r in rows}):
  rr=[r for r in rows if r['role']==role];key='associated_raw' if role.startswith('SIGNAL') or role=='ON_SOURCE' else 'detected_window'
  summary['discordant'][role]={'gap_only':sum(r['gap'][key] and not r['baseline'][key] for r in rr),'baseline_only':sum(r['baseline'][key] and not r['gap'][key] for r in rr)}
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

def run_synthetic(out):
 results=[]
 for seed in (7301,7302,7303):
  name=f'synthetic_{seed}';print('SYNTHETIC_START',seed,flush=True)
  calobs=[obs(noise('pink_colored',rng(seed,0,i)),f'{name}_cal_{i}') for i in range(320)]
  evaluation=[{'sample_id':f'{name}_noise_{i}','role':'NOISE','gate':None,'observation':obs(noise('pink_colored',rng(seed,1,i)),f'{name}_noise_{i}')} for i in range(200)]
  for i in range(10):
   base=noise('pink_colored',rng(seed,2,i));g=rng(seed,3,i)
   params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
   clean,a,b,lo,hi,metadata=generate_clean_signal('linear_chirp',4096,1024.,g,params)
   for snr in (8,18):
    signal_scaled=_scale_to_energy_snr(clean,float(np.std(base)),snr);x=base+signal_scaled
    np.testing.assert_allclose(np.linalg.norm(signal_scaled)/np.std(base),snr,rtol=1e-12)
    sample=f'{name}_inj_{i}_snr{snr}'
    evaluation.append({'sample_id':sample,'role':f'SIGNAL_SNR{snr}','injection_id':i,'snr':snr,'noise_sha256':arrayhash(base),'clean_unit_sha256':arrayhash(clean),'signal_sha256':arrayhash(signal_scaled),'generation_params':metadata,'gate':[a,b,lo,hi],'observation':obs(x,sample)})
  summary,_,_=compare_case(out/name,calobs,evaluation,False);results.append(summary)
  print('SYNTHETIC_DONE',seed,json.dumps(summary['arms']),flush=True)
 return results

def normalization_metrics(out):
 records=[];inputs=[]
 for seed in (8401,8402,8403):
  print('NORMALIZER_SEED',seed,flush=True)
  for model_index,model in enumerate(('gaussian_white','pink_colored')):
   for i in range(10):
    base=noise(model,rng(seed,model_index,0,i));baseobs=obs(base,'base')
    base_maps={s.name:DETECTORS[s.name].analyze(baseobs,return_maps=True).tf_maps for s in DEFAULT_V02_SCALES}
    for family_index,family in enumerate(('linear_chirp','band_limited_noise_burst')):
     g=rng(seed,model_index,1,i,family_index)
     if family=='linear_chirp':
      params={'t_start_s':1.25,'chirp_duration_s':.35,'f_start_hz':40.,'f_end_hz':240.,'phase_rad':float(g.uniform(0,2*np.pi))}
      clean,a,b,_,_,metadata=generate_clean_signal('linear_chirp',4096,1024.,g,params)
     else:
      a,b=1.375,1.625;t=np.arange(4096)/1024.;active=(t>=a)&(t<=b)
      full=signal.sosfiltfilt(signal.butter(4,[35,350],fs=1024.,btype='bandpass',output='sos'),g.normal(size=4096))
      clean=np.zeros(4096);clean[active]=full[active]*tukey(int(active.sum()),alpha=.15)
      metadata={'t_start_s':a,'t_end_s':b,'filter':'butter4_sosfiltfilt_35_350Hz','taper_alpha':.15}
     unit_id=f'{seed}_{model}_{i}_{family}'
     for snr in (0,8,18,36):
      scaled=np.zeros_like(base) if snr==0 else _scale_to_energy_snr(clean,float(np.std(base)),snr)
      mixture=base+scaled;input_id=unit_id+f'_snr{snr}'
      inputs.append({'id':input_id,'seed':seed,'model':model,'noise_index':i,'family':family,'snr':snr,'noise_sha256':arrayhash(base),'clean_sha256':arrayhash(clean),'signal_sha256':arrayhash(scaled),'mixture_sha256':arrayhash(mixture),'support':[a,b],'generation_params':metadata,'measured_energy_snr':float(np.linalg.norm(scaled)/np.std(base))})
      for spec in DEFAULT_V02_SCALES:
       m0=base_maps[spec.name];m1=DETECTORS[spec.name].analyze(obs(mixture,input_id),return_maps=True).tf_maps
       p=m1.power_spectrogram;s1=m1.background_psd;t1=m1.background_time_envelope;s0=m0.background_psd;t0=m0.background_time_envelope;z=m1.normalized_excess_map
       temp=p/s1[:,None]/t0[None,:];fully=p/s0[:,None]/t0[None,:]
       np.testing.assert_allclose(temp,z*(t1/t0)[None,:],rtol=1e-12,atol=1e-12)
       np.testing.assert_allclose(fully,z*(s1/s0)[:,None]*(t1/t0)[None,:],rtol=1e-12,atol=1e-12)
       frames=(m1.times>=a)&(m1.times<=b);freq=(m1.freqs>=35)&(m1.freqs<=350);assert frames.any() and freq.any()
       roi=freq[:,None]&frames[None,:]
       if snr==0:
        for key in ('power_spectrogram','background_psd','background_time_envelope','normalized_excess_map'):assert np.array_equal(getattr(m0,key),getattr(m1,key))
       records.append({'input_id':input_id,'seed':seed,'noise_model':model,'noise_index':i,'family':family,'snr':snr,'scale':spec.name,'roi_cells':int(roi.sum()),'temporal_ratio':float(np.median((t1/t0)[frames])),'spectral_ratio':float(np.median((s1/s0)[freq])),'temporal_counterfactual_Z_ratio':float(np.median(np.broadcast_to((t1/t0)[None,:],z.shape)[roi])),'full_counterfactual_Z_ratio':float(np.median(np.broadcast_to((s1/s0)[:,None]*(t1/t0)[None,:],z.shape)[roi])),'production_seeds':int((z[roi]>4).sum()),'temporal_frozen_seeds':int((temp[roi]>4).sum()),'fully_frozen_seeds':int((fully[roi]>4).sum()),'production_positive_excess':float(np.maximum(z[roi]-1,0).sum()),'temporal_frozen_positive_excess':float(np.maximum(temp[roi]-1,0).sum()),'fully_frozen_positive_excess':float(np.maximum(fully[roi]-1,0).sum())})
 save(out/'normalizer_inputs.json',inputs);save(out/'normalizer_records.json',records);csvsave(out/'normalizer_records.csv',records)
 aggregates=[]
 for model in ('gaussian_white','pink_colored'):
  for family in ('linear_chirp','band_limited_noise_burst'):
   for snr in (0,8,18,36):
    for spec in DEFAULT_V02_SCALES:
     rr=[r for r in records if r['noise_model']==model and r['family']==family and r['snr']==snr and r['scale']==spec.name]
     vals=[r['temporal_ratio'] for r in rr]
     aggregates.append({'noise_model':model,'family':family,'snr':snr,'scale':spec.name,'paired_records':len(rr),'temporal_ratio_median':float(np.median(vals)),'temporal_ratio_min':min(vals),'temporal_ratio_max':max(vals),'T_increased':sum(v>1 for v in vals),'T_decreased':sum(v<1 for v in vals),'T_unchanged':sum(v==1 for v in vals)})
 save(out/'normalizer_aggregates.json',aggregates)
 return {'base_noise_realizations':60,'mixtures':len(inputs),'scale_records':len(records),'identity_scale_records':sum(r['snr']==0 for r in records)}

def figures(out,summaries):
 fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for ax,site in zip(axes,SITES):
  rr=load(out/site/'evaluation.json');e=next(r for r in rr if r['role']=='ON_SOURCE')
  ax.bar(ARMS,[e[a]['near_score']/e[a]['threshold'] for a in ARMS],color=['steelblue','darkorange'])
  ax.axhline(1,color='red',ls='--');ax.set_title(site);ax.set_ylabel('Near-reference score / own threshold');ax.set_ylim(bottom=0)
 fig.suptitle('Developmental GW170814 comparison; each arm separately calibrated')
 fig.savefig(out/'real_gap_comparison.png',dpi=150);plt.close(fig)
 agg=load(out/'normalizer_aggregates.json');fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,layout='constrained')
 for row,model in enumerate(('gaussian_white','pink_colored')):
  for col,family in enumerate(('linear_chirp','band_limited_noise_burst')):
   ax=axes[row,col]
   for spec in DEFAULT_V02_SCALES:
    rr=[r for r in agg if r['noise_model']==model and r['family']==family and r['scale']==spec.name]
    ax.plot([r['snr'] for r in rr],[r['temporal_ratio_median'] for r in rr],'o-',label=spec.name)
   ax.axhline(1,color='gray',ls=':');ax.set_title(model+' / '+family);ax.set_xlabel('Root-energy SNR (linear)');ax.set_ylabel('Median paired temporal factor T1/T0');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=8);fig.suptitle('Synthetic normalizer response; medians of paired records, not event significance')
 fig.savefig(out/'normalizer_response.png',dpi=150);plt.close(fig)

def main():
 global ACTIVE
 pm=load(PARENT/'manifest.json');sources=pm['source_sha256'];assert all(digest(ROOT/f)==h for f,h in sources.items())
 assert not git('diff','HEAD','--','src','app','docs/workstation/index.html')
 rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
 out=ROOT/'results/gap_normalizer_controls'/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 for name,p in [('driver.py',Path(__file__)),('gap_link_prototype.py',HELPER),('protocol.md',PROTOCOL)]:(out/name).write_bytes(p.read_bytes())
 parents={str((PARENT/s/'summary.json').relative_to(ROOT)).replace('\\','/'):digest(PARENT/s/'summary.json') for s in SITES}
 m={'run_id':rid,'source_commit':git('rev-parse','HEAD'),'source_sha256':sources,'parent_sha256':parents,'driver_sha256':digest(Path(__file__)),'helper_sha256':digest(HELPER),'protocol_sha256':digest(PROTOCOL),'dirty_paths':git('status','--porcelain').splitlines(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'status':'RUNNING','stages':{}}
 save(out/'manifest.json',m);print('RUN_DIRECTORY',out,flush=True)
 fixtures=focused_checks();save(out/'focused_checks.json',fixtures)
 real=run_real(out);m['stages']['real_gap']='COMPLETE';save(out/'manifest.json',m)
 synthetic=run_synthetic(out);m['stages']['synthetic_gap']='COMPLETE';save(out/'manifest.json',m)
 save(out/'gap_summaries.json',real+synthetic)
 m['stages']['normalizer']=normalization_metrics(out);save(out/'manifest.json',m)
 figures(out,real+synthetic)
 assert all(digest(ROOT/f)==h for f,h in sources.items()) and all(digest(ROOT/f)==h for f,h in parents.items())
 assert digest(HELPER)==m['helper_sha256']
 m['status']='COMPLETE';m['source_and_parent_unchanged']=True;save(out/'manifest.json',m);print('COMPLETE',rid,flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  if ACTIVE is not None and (ACTIVE/'manifest.json').exists():
   m=load(ACTIVE/'manifest.json');m['status']='ERROR';m['error']={'type':type(exc).__name__,'message':str(exc)};save(ACTIVE/'manifest.json',m)
  raise
