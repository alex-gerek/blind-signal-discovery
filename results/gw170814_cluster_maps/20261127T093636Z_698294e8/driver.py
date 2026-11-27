"""Saved-event map diagnosis with independent normalization/topology checks; no tuning."""
from __future__ import annotations
import csv,hashlib,json,math,os,platform,subprocess,sys,uuid
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True
PARENT=ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'
PROTOCOL=ROOT/'docs/gw170814_cluster_map_protocol_2026-11-27.md'
SITES=('H1','L1')
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
INPUT_INFO={s:load(PARENT/s/'event_array_manifest.json') for s in SITES}
ALLOWED={str((ROOT/m['path']).resolve()).lower().replace('\\','/') for m in INPUT_INFO.values()}
DATA=str((ROOT/'data').resolve()).lower().replace('\\','/')+'/'
RESULTS=str((ROOT/'results').resolve()).lower().replace('\\','/')+'/'
ACTIVE=None

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
  writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if 'o3_holdout' in p or (p.startswith(DATA) and (p not in ALLOWED or writing)):raise PermissionError('Data access outside exact saved-array allowlist')
  if writing and p.startswith(RESULTS) and (ACTIVE is None or not p.startswith(str(ACTIVE).lower().replace('\\','/')+'/')):raise PermissionError('Historical result write blocked')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
import scipy
from scipy import ndimage,signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsd.detector import TemplateFreeDetector
from bsd.detector_v2 import DEFAULT_V02_SCALES
from bsd.generator import Observation

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=encode),encoding='utf-8')
def git(*a):return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def array_info(x):return {'shape':list(x.shape),'dtype':str(x.dtype),'sha256_C_order_bytes':hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()}
def close(a,b):assert math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-12),(a,b)
def cell_region(fi,ti,f,t,j):return (f[fi]>=35)&(f[fi]<=250)&(t[ti]>=j)&(t[ti]<j+1)

def independent_groups(mask):
 remaining=set(map(tuple,np.argwhere(mask)));groups=[]
 for root in map(tuple,np.argwhere(mask)):
  if root not in remaining:continue
  remaining.remove(root);stack=[root];members=[]
  while stack:
   a,b=stack.pop();members.append((a,b))
   for da in (-1,0,1):
    for db in (-1,0,1):
     q=(a+da,b+db)
     if (da or db) and q in remaining:remaining.remove(q);stack.append(q)
  groups.append(sorted(members))
 return groups

def main():
 global ACTIVE
 parent_manifest=load(PARENT/'manifest.json');sources=parent_manifest['source_sha256']
 assert all(digest(ROOT/f)==h for f,h in sources.items())
 assert not git('diff','HEAD','--','src','app','docs/workstation/index.html')
 rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
 out=ROOT/'results/gw170814_cluster_maps'/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 (out/'driver.py').write_bytes(Path(__file__).read_bytes());(out/'protocol.md').write_bytes(PROTOCOL.read_bytes())
 parent_paths=[PARENT/'manifest.json']+[PARENT/s/n for s in SITES for n in ('calibration.json','event_diagnostics.json','event_array_manifest.json','summary.json')]
 parents={str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in parent_paths}
 manifest={'run_id':rid,'source_commit':git('rev-parse','HEAD'),'source_sha256':sources,'parent_sha256':parents,'protocol_sha256':digest(PROTOCOL),'driver_sha256':digest(Path(__file__)),'dirty_paths':git('status','--porcelain').splitlines(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'status':'RUNNING','scope':'Saved two-array map reconstruction; unchanged detector; diagnostic metrics not calibrated'}
 save(out/'manifest.json',manifest);print('RUN_DIRECTORY',out,flush=True)
 summaries=[];components=[];pixels=[];regions=[];views={};checks={'scales':0,'qualified_clusters_reproduced':0,'all_components_graph_checked':0,'normalization_maps_checked':0}
 for site in SITES:
  meta=INPUT_INFO[site];npy=ROOT/meta['path'];assert digest(npy)==meta['sha256']
  x=np.load(npy,allow_pickle=False);assert x.shape==(4096,) and np.all(np.isfinite(x))
  obs=Observation(f'{site}_saved_event',x,1024.,4.)
  parent=load(PARENT/site/'event_diagnostics.json');cal=load(PARENT/site/'calibration.json')
  for spec in DEFAULT_V02_SCALES:
   cfg=spec.detector_config;r=TemplateFreeDetector(cfg).analyze(obs,calibration=None,return_maps=True);maps=r.tf_maps;assert maps is not None
   f=maps.freqs;t=maps.times;z=maps.normalized_excess_map;power=maps.power_spectrogram
   # Separate direct implementation of the published normalization equations.
   ff,tt,zz=signal.stft(x-np.median(x),fs=1024.,window=cfg.window,nperseg=cfg.nperseg,noverlap=cfg.noverlap,boundary=None,padded=False)
   independent_power=abs(zz[1:-1])**2;nf,nt=independent_power.shape
   fk=min(cfg.freq_median_kernel,nf if nf%2 else nf-1);tk=min(cfg.time_median_kernel,nt if nt%2 else nt-1)
   psd=np.maximum(ndimage.median_filter(np.median(independent_power,axis=1)/np.log(2),size=max(3,fk),mode='nearest'),1e-18)
   white=independent_power/psd[:,None]
   env=np.maximum(ndimage.median_filter(np.median(white,axis=0)/np.log(2),size=max(3,tk),mode='nearest'),1e-18)
   for a,b in ((f,ff[1:-1]),(t,tt),(power,independent_power),(maps.background_psd,psd),(maps.background_time_envelope,env),(z,white/env[None,:])):np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-12)
   checks['normalization_maps_checked']+=1
   mask=z>4;labels,n=ndimage.label(mask,structure=np.ones((3,3),dtype=int));groups=independent_groups(mask)
   assert n==len(groups)
   for i,g in enumerate(groups,1):assert set(g)==set(map(tuple,np.argwhere(labels==i)))
   checks['all_components_graph_checked']+=n
   known={c.candidate_id:c for c in r.candidates};previous={c['candidate_id']:c for c in parent['all_seed_qualified_raw_clusters'] if c['scale']==spec.name}
   assert set(known)==set(previous)
   qid=0;local_components=[];dt=t[1]-t[0];df=f[1]-f[0]
   for lid,g in enumerate(groups,1):
    coords=np.asarray(g);fi=coords[:,0];ti=coords[:,1];values=z[fi,ti];score=float(np.sum(values-1));qualified=len(g)>=2
    peak=int(np.argmax(values));a=max(0.,float(t[ti.min()]-.5*dt));b=min(4.,float(t[ti.max()]+.5*dt));lo=max(0.,float(f[fi.min()]-.5*df));hi=min(512.,float(f[fi.max()]+.5*df))
    c={'site':site,'scale':spec.name,'component_id':lid,'qualified':qualified,'n_bins':len(g),'raw_score':score,'t_start':a,'t_end':b,'f_min':lo,'f_max':hi,'peak_time':float(t[ti[peak]]),'peak_freq':float(f[fi[peak]]),'peak_Z':float(values[peak]),'parent_bbox_gate':a<=2 and b>=1 and lo<=250 and hi>=35,'regional_seed_contributions':[],'regional_seed_counts':[]}
    for j in range(4):
     inside=cell_region(fi,ti,f,t,j);c['regional_seed_counts'].append(int(sum(inside)));c['regional_seed_contributions'].append(float(np.sum(values[inside]-1)))
    if qualified:
     qid+=1;obj=known[qid];prev=previous[qid]
     for key,value in [('anomaly_score',score),('t_start',a),('t_end',b),('f_min',lo),('f_max',hi),('peak_time',c['peak_time']),('peak_freq',c['peak_freq']),('n_bins',len(g))]:close(getattr(obj,key),value);close(prev[key],value)
     assert prev['near_reference']==c['parent_bbox_gate'];checks['qualified_clusters_reproduced']+=1
     c['candidate_id']=qid
    else:c['candidate_id']=None
    for (ii,jj),value in zip(g,values):
     pixels.append({'site':site,'scale':spec.name,'freq_index':int(ii),'time_index':int(jj),'frequency_hz':float(f[ii]),'time_s':float(t[jj]),'power':float(power[ii,jj]),'background_psd':float(psd[ii]),'background_time':float(env[jj]),'Z':float(value),'component_id':lid,'component_size':len(g),'qualified':qualified})
    components.append(c);local_components.append(c)
   assert qid==len(known)
   close(r.max_anomaly_score,parent['per_scale_raw_max'][spec.name])
   stats=cal['scale_stats'][spec.name];threshold=cal['global_threshold'];boundary=stats['median_null_score']+threshold*stats['spread_scale']
   associated=[c for c in local_components if c['qualified'] and c['parent_bbox_gate']]
   winner=max(associated,key=lambda c:c['raw_score']);std=(winner['raw_score']-stats['median_null_score'])/stats['spread_scale']
   for j in range(4):
    roi=(f[:,None]>=35)&(f[:,None]<=250)&(t[None,:]>=j)&(t[None,:]<j+1)
    touched=[c for c in local_components if c['regional_seed_counts'][j]>0]
    total=float(np.sum((z-1)[roi&mask]));largest=max((c['regional_seed_contributions'][j] for c in touched),default=0.)
    largestq=max((c['regional_seed_contributions'][j] for c in touched if c['qualified']),default=0.)
    close(total,sum(c['regional_seed_contributions'][j] for c in touched))
    regions.append({'site':site,'scale':spec.name,'region_start_s':j,'region_end_s':j+1,'cell_band_hz':[35,250],'total_cells':int(roi.sum()),'seed_cells':int((roi&mask).sum()),'all_components_touching':len(touched),'qualified_components_touching':sum(c['qualified'] for c in touched),'seed_excess_sum':total,'largest_component_regional_contribution':largest,'largest_qualified_regional_contribution':largestq,'largest_share':largest/total if total else None})
   event_time=(t>=1)&(t<2);outside=~event_time
   summaries.append({'site':site,'scale':spec.name,'raw_clusters_qualified':qid,'all_components':n,'seed_bins':int(mask.sum()),'singleton_components':sum(not c['qualified'] for c in local_components),'near_reference_winner':winner,'near_standardized_score':std,'frozen_global_threshold':threshold,'raw_acceptance_boundary':boundary,'accepted_near_reference':std>threshold,'time_envelope_median_reference':float(np.median(env[event_time])),'time_envelope_median_other':float(np.median(env[outside])),'map_fingerprints':{key:array_info(value) for key,value in [('frequencies',f),('times',t),('power',power),('background_psd',psd),('background_time',env),('Z',z)]}})
   views[(site,spec.name)]=(f,t,z);checks['scales']+=1
   print('MAP_CHECKED',site,spec.name,'near raw',winner['raw_score'],'boundary',boundary,'seeds',int(mask.sum()),'components',n,flush=True)
 assert checks['qualified_clusters_reproduced']==286
 save(out/'component_records.json',components);save(out/'region_metrics.json',regions);save(out/'scale_summary.json',summaries)
 for name,data in [('seed_pixels.csv',pixels),('region_metrics.csv',[{k:v for k,v in row.items() if k!='cell_band_hz'} for row in regions])]:
  with (out/name).open('w',encoding='utf-8',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True,sharey=True,layout='constrained')
 for row,site in enumerate(SITES):
  for col,spec in enumerate(DEFAULT_V02_SCALES):
   ax=axes[row,col];f,t,z=views[(site,spec.name)]
   im=ax.pcolormesh(t,f,np.log2(np.maximum(z,1e-12)),shading='auto',cmap='viridis',vmin=-1,vmax=4)
   ax.contour(t,f,(z>4).astype(float),levels=[.5],colors='white',linewidths=.45)
   ax.add_patch(Rectangle((1,35),1,215,fill=False,edgecolor='red',lw=1.2,ls='--'))
   ax.axvline(1.5,color='black',ls=':',lw=1);ax.set_xlim(.75,2.25);ax.set_ylim(20,300);ax.set_title(f'{site}: {spec.name}')
   if row==1:ax.set_xlabel('Seconds within saved four-second window')
   if col==0:ax.set_ylabel('Frequency [Hz]')
 cb=fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.85,ticks=[-1,0,1,2,3,4]);cb.ax.set_yticklabels(['0.5','1','2','4','8','16']);cb.set_label('Normalized bin power Z (log2 color scale)')
 fig.suptitle('GW170814: fixed Z>4 seed contours, three scales; dashed red = diagnostic gate')
 fig.savefig(out/'normalized_maps.png',dpi=160);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(12,3.7),sharey=True,layout='constrained')
 for ax,spec in zip(axes,DEFAULT_V02_SCALES):
  for site in SITES:
   rr=[r for r in regions if r['site']==site and r['scale']==spec.name]
   ax.plot([r['region_start_s']+.5 for r in rr],[r['largest_share'] for r in rr],'o-',label=site)
  ax.axvspan(1,2,color='gray',alpha=.12);ax.set_title(spec.name);ax.set_ylim(0,1.05);ax.set_xticks([.5,1.5,2.5,3.5]);ax.set_xlabel('Center of diagnostic one-second region [s]');ax.grid(alpha=.2)
 axes[0].set_ylabel('Largest component share of regional seed excess');axes[-1].legend()
 fig.suptitle('Descriptive concentration of Z>4 seed excess in 35–250 Hz; not significance')
 fig.savefig(out/'regional_concentration.png',dpi=150);plt.close(fig)
 assert all(digest(ROOT/f)==h for f,h in sources.items()) and all(digest(ROOT/f)==h for f,h in parents.items())
 assert all(digest(ROOT/m['path'])==m['sha256'] for m in INPUT_INFO.values())
 checks.update({'status':'VERIFIED','seed_pixels_saved':len(pixels),'input_arrays_and_source_unchanged':True,'no_new_detection_rule':True})
 save(out/'independent_checks.json',checks)
 manifest['status']='COMPLETE';manifest['checks']=checks;manifest['input_arrays']=INPUT_INFO;save(out/'manifest.json',manifest)
 print('COMPLETE',json.dumps(checks),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  if ACTIVE is not None and (ACTIVE/'manifest.json').exists():
   m=load(ACTIVE/'manifest.json');m['status']='ERROR';m['error']={'type':type(exc).__name__,'message':str(exc)};save(ACTIVE/'manifest.json',m)
  raise
