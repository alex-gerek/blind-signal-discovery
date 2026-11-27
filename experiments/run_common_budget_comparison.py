"""Independent common-budget comparison; frozen methods, disjoint synthetic streams."""
from __future__ import annotations
import hashlib,json,os,platform,subprocess,sys,uuid
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'));sys.dont_write_bytecode=True
DATA=str((ROOT/'data').resolve()).lower().replace('\\','/')+'/'
RESULTS=str((ROOT/'results').resolve()).lower().replace('\\','/')+'/'
ACTIVE=None

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  if p.startswith(DATA) or 'o3_holdout' in p:raise PermissionError('All real/derived waveform and O3 reads blocked')
  mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
  writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if writing and p.startswith(RESULTS) and (ACTIVE is None or not p.startswith(str(ACTIVE).lower().replace('\\','/')+'/')):raise PermissionError('Historical output writes blocked')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
import scipy
from scipy import signal,stats
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bsd.generator import Observation,generate_clean_signal,_scale_to_energy_snr
from bsd.noise import generate_noise
from temporal_envelope_prototype import raw_pair,TableDetector,score_table,focused_checks
ARMS=('baseline','slow');REGIMES=('stationary','smooth_drift','within_window_step');FAMILIES=('linear_chirp','band_limited_noise_burst')
PROTOCOL=ROOT/'docs/common_budget_comparison_protocol_2026-11-27.md';HELPER=ROOT/'experiments/temporal_envelope_prototype.py'

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ah(x):return hashlib.sha256(np.asarray(x,dtype='<f8',order='C').tobytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def rng(*key):return np.random.default_rng(np.random.SeedSequence(key))
def observation(x,name):return Observation(name,x,1024.,4.)
def intersects(c,g):return g is not None and c.t_start<=g[1] and c.t_end>=g[0] and c.f_min<=g[3] and c.f_max>=g[2]

def noise_record(seed,index,family_key=0):
 regime=int(rng(seed,10,index,family_key).integers(3));g=rng(seed,4,index,family_key,regime);t=np.arange(4096)/1024.
 base=generate_noise('pink_colored',4096,1024.,rng(seed,0,index,family_key),sigma=1.,params={'alpha':1.,'f_knee':16.})[0]
 params={}
 if regime==0:envelope=np.ones(4096);scale=1.;x=base.copy()
 elif regime==1:
  phase=float(g.uniform(0,2*np.pi));envelope=1+.5*np.sin(2*np.pi*.25*t+phase);x=base*envelope;scale=float(1/np.sqrt(np.mean(x*x)));x=x*scale;params={'phase_rad':phase,'modulation_depth':.5,'frequency_hz':.25}
 else:
  when=float(g.uniform(1,3));factor=2. if g.random()<.5 else .5;envelope=np.where(t<when,1.,factor);x=base*envelope;scale=float(1/np.sqrt(np.mean(x*x)));x=x*scale;params={'step_time_s':when,'step_factor':factor}
 metadata={'seed':seed,'index':index,'family_key':family_key,'regime':REGIMES[regime],'base_noise_sha256':ah(base),'noise_sha256':ah(x),'envelope_sha256':ah(envelope),'global_scale':scale,'profile_params':params}
 return x,metadata

def clean_record(seed,index,fi):
 g=rng(seed,3,index,fi)
 if fi==0:
  params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
  x,a,b,lo,hi,meta=generate_clean_signal('linear_chirp',4096,1024.,g,params);return x,[a,b,lo,hi],meta
 active=(np.arange(4096)/1024.>=1.375)&(np.arange(4096)/1024.<=1.625)
 full=signal.sosfiltfilt(signal.butter(4,[35,350],btype='bandpass',fs=1024.,output='sos'),g.normal(size=4096));x=np.zeros(4096);x[active]=full[active]*tukey(int(active.sum()),alpha=.15)
 return x,[1.375,1.625,35.,350.],{'t_start_s':1.375,'t_end_s':1.625,'filter':'butter4_sosfiltfilt_35_350Hz','taper_alpha':.15}

def cp(k,n):
 if n==0:return {'count':k,'n':n,'rate':None,'lower95':None,'upper95':None}
 return {'count':k,'n':n,'rate':k/n,'lower95':0. if k==0 else float(stats.beta.ppf(.025,k,n-k+1)),'upper95':1. if k==n else float(stats.beta.ppf(.975,k+1,n-k))}

def exact_paired(new_only,old_only):
 d=new_only+old_only
 return 1. if d==0 else float(stats.binomtest(new_only,d,.5,alternative='two-sided').pvalue)

def evaluate_one(o,pair,cals,gate,metadata,candidates):
 row={**metadata,'sample_id':o.sample_id,'input_sha256':ah(o.signal),'evaluation_gate':gate}
 for arm in ARMS:
  cal=cals[arm];r=TableDetector({o.sample_id:pair[arm]}).analyze(o,cal)
  row[arm]={'score':r.max_global_score,'threshold':cal.global_threshold,'empirical_p':cal.empirical_p_value(r.max_global_score)[0],'operating_alpha':cal.target_fap,'detected_window':r.detected,'associated':any(intersects(c,gate) for c in r.candidates),'per_scale_raw_max':r.per_scale_raw_max}
  for c in r.candidates:candidates.append({'sample_id':o.sample_id,'arm':arm,**c.to_dict()})
 return row

def inference(out,noise_rows,signal_rows,budget):
 noise_results={};conditional={}
 for arm in ARMS:
  v=cp(sum(r[arm]['detected_window'] for r in noise_rows),len(noise_rows));v['false_alarm_windows_per_hour_equivalent']=v['rate']*900;noise_results[arm]=v
  conditional[arm]={regime:cp(sum(r[arm]['detected_window'] for r in noise_rows if r['generation']['regime']==regime),sum(r['generation']['regime']==regime for r in noise_rows)) for regime in REGIMES}
 budget_supported=all(r['upper95']<=budget['budget_alpha'] for r in noise_results.values())
 groups=[]
 for family in FAMILIES:
  for snr in (8,18):
   rr=[r for r in signal_rows if r['family']==family and r['snr']==snr]
   b=sum(r['slow']['associated'] and not r['baseline']['associated'] for r in rr);c=sum(r['baseline']['associated'] and not r['slow']['associated'] for r in rr)
   group={'family':family,'snr':snr,'n':len(rr),'baseline':cp(sum(r['baseline']['associated'] for r in rr),len(rr)),'slow':cp(sum(r['slow']['associated'] for r in rr),len(rr)),'slow_only':b,'baseline_only':c,'primary':family=='band_limited_noise_burst' and snr==8}
   if group['primary']:group['exact_mcnemar_p_two_sided']=exact_paired(b,c)
   group['by_regime']={regime:{arm:cp(sum(r[arm]['associated'] for r in rr if r['generation']['regime']==regime),sum(r['generation']['regime']==regime for r in rr)) for arm in ARMS} for regime in REGIMES}
   groups.append(group)
 primary=next(g for g in groups if g['primary']);evidence=budget_supported and primary['slow']['count']>primary['baseline']['count'] and primary['exact_mcnemar_p_two_sided']<.05
 result={'budget_alpha':budget['budget_alpha'],'operating_alpha':budget['operating_alpha'],'noise_window_results':noise_results,'noise_by_regime':conditional,'budget_supported_by_evaluation_upper_bounds':budget_supported,'noise_discordance':{'slow_only':sum(r['slow']['detected_window'] and not r['baseline']['detected_window'] for r in noise_rows),'baseline_only':sum(r['baseline']['detected_window'] and not r['slow']['detected_window'] for r in noise_rows)},'signal_groups':groups,'primary_improvement_evidence_for_this_generator':evidence,'automatic_adoption':False,'scope':'One fixed calibration realization; IID synthetic mixture; one predeclared primary association endpoint; secondary groups descriptive'}
 save(out/'inference.json',result)
 fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for ax,results,title,limit in [(axes[0],noise_results,'Held-out noise: count / 3600',None),(axes[1],{a:primary[a] for a in ARMS},'Primary: weak burst associations / 300',None)]:
  vals=[results[a]['count'] for a in ARMS];bars=ax.bar(ARMS,vals,color=['steelblue','darkorange'])
  for bar,value in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,value+.3,str(value),ha='center')
  ax.set_ylim(0,max(vals+[1])*1.25+1);ax.set_title(title);ax.grid(axis='y',alpha=.2)
 fig.suptitle('Frozen independent comparison: separate calibration, common conservative budget')
 fig.savefig(out/'common_budget_counts.png',dpi=150);plt.close(fig)
 return result

def main():
 global ACTIVE
 sources=load(ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7/manifest.json')['source_sha256']
 assert all(digest(ROOT/p)==h for p,h in sources.items());assert not git('diff','HEAD','--','src','app','docs/workstation/index.html')
 rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8];out=ROOT/'results/common_budget_comparison'/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 for name,p in [('driver.py',Path(__file__)),('temporal_envelope_prototype.py',HELPER),('protocol.md',PROTOCOL)]:(out/name).write_bytes(p.read_bytes())
 m={'run_id':rid,'source_commit':git('rev-parse','HEAD'),'source_sha256':sources,'driver_sha256':digest(Path(__file__)),'helper_sha256':digest(HELPER),'protocol_sha256':digest(PROTOCOL),'dirty_paths':git('status','--porcelain').splitlines(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'status':'RUNNING','stages':{}}
 save(out/'manifest.json',m);print('RUN_DIRECTORY',out,flush=True)
 fixtures=focused_checks();assert exact_paired(0,0)==1 and exact_paired(5,0)==.0625 and exact_paired(6,0)==.03125
 assert abs(cp(0,3600)['upper95']-(1-.025**(1/3600)))<1e-12
 save(out/'focused_checks.json',{'prototype':fixtures,'exact_paired_boundaries':'PASS','clopper_pearson_zero_boundary':'PASS'})
 n=1800;alpha=.01;delta=.025;tail=max(k for k in range(1,n+1) if stats.binom.cdf(k-1,n,alpha)<=delta)
 budget={'budget_alpha':alpha,'budget_false_alarm_windows_per_hour_equivalent':alpha*900,'normalization_windows':200,'threshold_windows':n,'failure_allowance_per_arm':delta,'tail_rank_m':tail,'order_k_1based':n+1-tail,'operating_alpha':tail/(n+1),'tolerance_tail_mass_upper975':float(stats.beta.ppf(.975,tail,n+1-tail)),'binomial_failure_probability_at_budget':float(stats.binom.cdf(tail-1,n,alpha)),'next_rank_failure_probability':float(stats.binom.cdf(tail,n,alpha)),'assumptions':'IID predefined mixture, independent A/B; strict crossing/ties conservative; joint method assurances by union bound'}
 assert tail==10 and budget['tolerance_tail_mass_upper975']<=alpha;save(out/'budget_rule.json',budget)
 features=[];calobs=[]
 for split,seed,count in [('A',15001,200),('B',15002,1800)]:
  for i in range(count):
   x,gen=noise_record(seed,i);o=observation(x,f'{split}_{i:04d}');pair=raw_pair(o);calobs.append(o)
   features.append({'sample_id':o.sample_id,'split':split,'generation':gen,'input_sha256':ah(x),**{arm:{name:t[0] for name,t in pair[arm].items()} for arm in ARMS}})
   if (i+1)%400==0:print('CALIBRATION_FEATURES',split,i+1,flush=True)
 save(out/'calibration_raw_features.json',features);cals={};frozen={}
 for arm in ARMS:
  cal=TableDetector(score_table(features,arm)).calibrate(calobs,target_fap=budget['operating_alpha'],norm_fraction=.1,calibration_seed=0);assert (cal.n_norm_samples,cal.n_threshold_samples)==(200,1800);assert cal.global_threshold==cal.sorted_global_null_scores[budget['order_k_1based']-1];cals[arm]=cal
  save(out/f'{arm}_calibration.json',cal.to_summary_dict());save(out/f'{arm}_null_scores.json',cal.sorted_global_null_scores.tolist())
  for name in [f'{arm}_calibration.json',f'{arm}_null_scores.json']:frozen[name]=digest(out/name)
 frozen['budget_rule.json']=digest(out/'budget_rule.json');save(out/'threshold_freeze.json',{'timestamp_utc':datetime.now(timezone.utc).isoformat(),'files_sha256':frozen,'heldout_records_generated':0,'thresholds_will_not_be_retuned':True})
 m['stages']['calibration']='FROZEN';m['threshold_freeze_sha256']=digest(out/'threshold_freeze.json');save(out/'manifest.json',m);print('THRESHOLDS_FROZEN',flush=True)
 del calobs
 noise_rows=[];signal_rows=[];candidates=[]
 for i in range(3600):
  x,gen=noise_record(15003,i);o=observation(x,f'noise_{i:04d}');pair=raw_pair(o)
  noise_rows.append(evaluate_one(o,pair,cals,None,{'role':'NOISE','generation':gen},candidates))
  if (i+1)%900==0:print('HELDOUT_NOISE',i+1,flush=True)
 save(out/'heldout_noise.json',noise_rows);m['stages']['noise']='COMPLETE';save(out/'manifest.json',m)
 for fi,family in enumerate(FAMILIES):
  for i in range(300):
   base,gen=noise_record(15004,i,fi);clean,gate,params=clean_record(15004,i,fi)
   for snr in (8,18):
    added=_scale_to_energy_snr(clean,float(np.std(base)),snr);x=base+added;assert abs(np.linalg.norm(added)/np.std(base)-snr)<1e-10;o=observation(x,f'{family}_{i:04d}_snr{snr}');pair=raw_pair(o)
    signal_rows.append(evaluate_one(o,pair,cals,gate,{'role':'SIGNAL','family':family,'unit_index':i,'snr':snr,'generation':gen,'clean_sha256':ah(clean),'signal_sha256':ah(added),'generation_params':params},candidates))
   if (i+1)%100==0:print('HELDOUT_SIGNALS',family,i+1,flush=True)
 save(out/'heldout_signals.json',signal_rows);save(out/'candidates.json',candidates)
 sets={'A':{r['input_sha256'] for r in features if r['split']=='A'},'B':{r['input_sha256'] for r in features if r['split']=='B'},'noise':{r['input_sha256'] for r in noise_rows},'signal_noise_partners':{r['generation']['noise_sha256'] for r in signal_rows}}
 assert [len(sets[k]) for k in ('A','B','noise','signal_noise_partners')]==[200,1800,3600,600]
 names=list(sets);assert all(not sets[a]&sets[b] for i,a in enumerate(names) for b in names[i+1:]);save(out/'split_isolation.json',{'unique_hash_counts':{k:len(v) for k,v in sets.items()},'cross_split_noise_hash_intersections':0,'signal_amplitudes_paired_by_unit':True})
 result=inference(out,noise_rows,signal_rows,budget)
 assert all(digest(out/name)==h for name,h in frozen.items()) and digest(out/'threshold_freeze.json')==m['threshold_freeze_sha256']
 assert all(digest(ROOT/p)==h for p,h in sources.items()) and digest(HELPER)==m['helper_sha256']
 m['stages']['signals']='COMPLETE';m['status']='COMPLETE';m['frozen_files_and_sources_unchanged']=True;save(out/'manifest.json',m)
 print('COMPLETE',rid,'budget_supported',result['budget_supported_by_evaluation_upper_bounds'],'primary_improvement_evidence',result['primary_improvement_evidence_for_this_generator'],flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  if ACTIVE is not None and (ACTIVE/'manifest.json').exists():
   m=load(ACTIVE/'manifest.json');m['status']='ERROR';m['error']={'type':type(exc).__name__,'message':str(exc)};save(ACTIVE/'manifest.json',m)
  raise
