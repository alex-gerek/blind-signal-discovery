"""Prespecified, bounded GW170814 grid/context perturbation; no production edits."""
from __future__ import annotations
import csv,hashlib,importlib.util,json,os,platform,subprocess,sys,time,uuid
from dataclasses import asdict
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True
PARENT=ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'
PROTOCOL=ROOT/'docs/gw170814_window_grid_protocol_2026-11-27.md'
BASE='32388bedf6a1222cdb27b5531e3302787014dd28'
SITES=('H1','L1');PHASES=(0,-2,-1,1)
INPUTS={str((ROOT/'data/gwosc/GW170814'/f'{s[0]}-{s}_GWOSC_4KHZ_R1-1186739814-4096.hdf5').resolve()).lower().replace('\\','/') for s in SITES}
DATA=str((ROOT/'data').resolve()).lower().replace('\\','/')+'/'
RESULTS=str((ROOT/'results').resolve()).lower().replace('\\','/')+'/'
ACTIVE=None

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
  writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if 'o3_holdout' in p or (p.startswith(DATA) and (p not in INPUTS or writing)):raise PermissionError('Reserved/nonallowlisted data access blocked')
  if writing and p.startswith(RESULTS) and (ACTIVE is None or not p.startswith(str(ACTIVE).lower().replace('\\','/')+'/')):raise PermissionError('Historical result write blocked')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
import scipy
spec=importlib.util.spec_from_file_location('historical_grid_helpers',ROOT/'experiments/run_gw170814_exploratory_search.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=encode),encoding='utf-8')
def near(c,local):return c.t_start<=local+.5 and c.t_end>=local-.5 and c.f_min<=250 and c.f_max>=35

def main():
 global ACTIVE
 assert git('rev-parse','HEAD')==BASE
 assert not git('diff','HEAD','--','src','app','docs/workstation/index.html')
 pm=load(PARENT/'manifest.json')
 assert all(digest(ROOT/f)==h for f,h in pm['source_sha256'].items())
 rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
 out=ROOT/'results/gw170814_window_grid'/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 (out/'protocol.md').write_bytes(PROTOCOL.read_bytes());(out/'driver.py').write_bytes(Path(__file__).read_bytes())
 original_hashes={str(PARENT.relative_to(ROOT)).replace('\\','/')+'/'+f:digest(PARENT/f) for f in ['manifest.json','H1/summary.json','L1/summary.json','H1/input_provenance.json','L1/input_provenance.json']}
 m={'run_id':rid,'source_commit':BASE,'source_sha256':pm['source_sha256'],'parent_sha256':original_hashes,'protocol_sha256':digest(PROTOCOL),'driver_sha256':digest(Path(__file__)),'dirty_paths':git('status','--porcelain').splitlines(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'phases_s':PHASES,'status':'RUNNING','variants':[]}
 save(out/'manifest.json',m);print('RUN_DIRECTORY',out,flush=True)
 for site in SITES:
  inputpath=old.DETECTORS[site];expected=load(PARENT/site/'input_provenance.json')['original']
  assert digest(inputpath)==expected['sha256']
  record=old.load_gwosc_hdf5(inputpath);inj,names=old._load_inj_mask(inputpath)
  psd,bad=old._extract_exact_interval_segments(record,old.PSD_CORE_INTERVAL);assert not bad and len(psd)==60
  model=old.fit_background_model(site,psd,mode=old.PREPROCESS_MODE,target_sample_rate=1024.,bandpass_hz=old.BANDPASS_HZ)
  sd=out/site;sd.mkdir();save(sd/'input_provenance.json',{'original':expected,'model':asdict(model),'injection_shortnames':names})
  refcal=None
  for phase in PHASES:
   started=time.monotonic();dest=sd/f'phase_{phase:+d}';dest.mkdir();print('VARIANT_START',site,phase,flush=True)
   blocks=[(a+phase,b+phase) for a,b in (old.CAL_A_CORE_INTERVAL,old.SEARCH_CORE_INTERVAL,old.CAL_B_CORE_INTERVAL)]
   footprints=[(old.PSD_CORE_INTERVAL[0]-1,old.PSD_CORE_INTERVAL[1]+1)]+[(a-1,b+1) for a,b in blocks]
   assert all(record.gps_start<=a<b<=record.gps_end for a,b in footprints)
   assert all(footprints[i][1]<footprints[i+1][0] for i in range(3))
   extracted=[]
   for interval,count in zip(blocks,[320,160,387]):
    segs,bad=old._extract_exact_interval_segments(record,interval)
    assert not bad and len(segs)==count,(site,phase,bad)
    extracted.append(segs)
   ca,target,cb=extracted;det=old.MultiResolutionDetector()
   calobs=[old.preprocess_segment(s,f'{site}_{phase}_cal_{i}',model) for i,s in enumerate(ca+cb)]
   cal=det.calibrate(calobs,target_fap=.01,norm_fraction=.30,calibration_seed=0)
   assert (cal.n_norm_samples,cal.n_threshold_samples)==(212,495) and cal.global_threshold is not None
   if phase==0:refcal=cal
   assert refcal is not None
   chrono=np.array([det.analyze(o,cal).max_global_score for o in calobs[212:]])
   audit=old.classify_calibration_status(707,old.estimate_calibration_uncertainty(chrono,target_fap=.01),target_fap=.01)
   save(dest/'calibration.json',cal.to_summary_dict());save(dest/'null_scores.json',cal.sorted_global_null_scores);save(dest/'calibration_chronological_scores.json',chrono);save(dest/'audit.json',audit.to_dict())
   rows=[];candidates=[];event=None
   for i,s in enumerate(target):
    obs=old.preprocess_segment(s,f'{site}_{phase}_target_{i}',model);r=det.analyze(obs,cal)
    veto=old.apply_vetoes(obs,r.candidates,old.FROZEN_VETO_CONFIG);vmap={d.candidate.candidate_id:d for d in veto.decisions}
    gps=s.gps_start;local=old.EVENT_GPS-gps;on=0<=local<4
    role='ON_SOURCE' if on else ('OFF_SOURCE' if gps+5<=old.EVENT_GPS-16 or gps-1>=old.EVENT_GPS+16 else 'GUARD')
    flags=old._inspect_flags_at_interval(record,inj,names,gps,gps+4)
    row={'window_index':i,'gps_start':gps,'role':role,'score':r.max_global_score,'threshold':cal.global_threshold,'empirical_p':cal.empirical_p_value(r.max_global_score)[0],'detected_raw':r.detected,'detected_post_veto':veto.detected_after_veto,**flags}
    rows.append(row)
    for c in r.candidates:
     v=vmap[c.candidate_id];candidates.append({**c.to_dict(),'window_index':i,'window_gps_start':gps,'role':role,'vetoed':v.vetoed,'veto_reasons':v.reason_codes})
    if on:
     reference=det.analyze(obs,refcal);sub=[]
     for scale,(_,clusters,_) in det._raw_analyze_all_scales(obs).items():
      for c in clusters:sub.append({**c.to_dict(),'scale':scale,'near_reference':near(c,local),'standardized_matched':cal.scale_stats[scale].standardize(c.anomaly_score),'standardized_zero_calibration':refcal.scale_stats[scale].standardize(c.anomaly_score)})
     assert sub
     nearclusters=[c for c in sub if c['near_reference']];assert nearclusters
     nearbest=max(nearclusters,key=lambda c:c['standardized_matched'])
     rn=max(c['standardized_zero_calibration'] for c in nearclusters)
     event={**row,'event_local_s':local,'near_score':nearbest['standardized_matched'],'near_empirical_p':cal.empirical_p_value(nearbest['standardized_matched'])[0],'associated_raw':any(near(c,local) for c in r.candidates),'associated_post_veto':any(near(v.candidate,local) and not v.vetoed for v in veto.decisions),'per_scale_raw_max':r.per_scale_raw_max,'per_scale_std_max':r.per_scale_std_max,'winning_scale':r.winning_scale,'zero_calibration':{'score':reference.max_global_score,'near_score':rn,'threshold':refcal.global_threshold,'empirical_p':refcal.empirical_p_value(reference.max_global_score)[0],'near_empirical_p':refcal.empirical_p_value(rn)[0],'associated_raw':any(near(c,local) for c in reference.candidates)},'processed_event_array_sha256':hashlib.sha256(obs.signal.tobytes()).hexdigest(),'processed_array_hash_definition':'C-order float64 signal bytes, excluding NPY header','quality':flags}
     save(dest/'event_clusters.json',sub);save(dest/'event_summary.json',event)
   assert event is not None and len(rows)==160
   off=[r for r in rows if r['role']=='OFF_SOURCE']
   summary={'site':site,'phase_s':phase,'n_calibration':707,'n_target':160,'event':event,'offsource_windows':len(off),'offsource_exposure_s':4*len(off),'offsource_raw':sum(r['detected_raw'] for r in off),'offsource_post_veto':sum(r['detected_post_veto'] for r in off),'audit':audit.to_dict(),'elapsed_s':time.monotonic()-started,'padded_block_footprints':footprints}
   if phase==0:
    parent=load(PARENT/site/'summary.json')
    for k in ('score','threshold','empirical_p'):assert abs(event[k]-parent['event'][k])<1e-10,(site,k)
    assert summary['offsource_raw']==parent['offsource_raw'] and summary['offsource_post_veto']==parent['offsource_post_veto']
    assert event['detected_raw']==parent['event']['detected_raw'] and event['detected_post_veto']==parent['event']['detected_post_veto']
   with (dest/'timeline.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
   save(dest/'candidates.json',candidates);save(dest/'summary.json',summary);m['variants'].append(summary);save(out/'manifest.json',m)
   print('VARIANT_DONE',site,phase,'near',event['near_score'],'threshold',event['threshold'],'associated',event['associated_raw'],'off',summary['offsource_raw'],summary['offsource_post_veto'],len(off),flush=True)
 assert all(digest(ROOT/f)==h for f,h in m['source_sha256'].items())
 assert all(digest(ROOT/f)==h for f,h in original_hashes.items())
 m['status']='COMPLETE';m['original_sources_and_parent_unchanged']=True;save(out/'manifest.json',m)
 print('RUN_COMPLETE',rid,flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  if ACTIVE is not None:
   mp=ACTIVE/'manifest.json'
   if mp.exists():
    m=load(mp);m['status']='ERROR';m['error']={'type':type(exc).__name__,'message':str(exc)};save(mp,m)
  raise
