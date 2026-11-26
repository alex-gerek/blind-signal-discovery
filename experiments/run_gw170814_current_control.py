"""Bounded historical-recipe GW170814 replay with explicit miss diagnostics."""
from __future__ import annotations
import csv,hashlib,importlib.util,json,os,platform,subprocess,sys,time,uuid
from dataclasses import asdict
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"));sys.dont_write_bytecode=True
ALLOWED={str(ROOT/"data/gwosc/GW170814"/f"{prefix}-{site}_GWOSC_4KHZ_R1-1186739814-4096.hdf5").replace("\\","/").lower()
 for prefix,site in (("H","H1"),("L","L1"),("V","V1"))}
DATA=str(ROOT/"data").replace("\\","/").lower()+"/"
DERIVED=str(ROOT/"data/derived/gw170814_current_control").replace("\\","/").lower()+"/"
def guard(event,args):
 if event=="open" and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).replace("\\","/").lower()
  if "o3_holdout" in p:raise PermissionError("Reserved manifest/data blocked")
  mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
  writing=(isinstance(mode,str) and any(c in mode for c in "wax+")) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if p.startswith(DATA) and not p.startswith(DERIVED):
   if p not in ALLOWED or writing:raise PermissionError("Nonallowlisted/original data access blocked")
  if writing and p.startswith(str(ROOT/"results/gw170814_exploratory_search").replace("\\","/").lower()+"/"):
   raise PermissionError("Historical output write blocked")
 if event in ("socket.connect","socket.getaddrinfo"):raise PermissionError("Network blocked")
sys.addaudithook(guard)
import numpy as np
import scipy
spec=importlib.util.spec_from_file_location("historical_gw170814",ROOT/"experiments/run_gw170814_exploratory_search.py")
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
from bsd.coincidence import deduplicate_candidates_by_window,match_coincident_candidates
BASE="8f78ca5b141ea656957dc69aa275daac4d71eeb0"
PROTOCOL=ROOT/"docs/gw170814_current_control_protocol_2026-11-26.md"
ACTIVE=None
def git(*a):return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
def digest(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=encode),encoding="utf-8")
def writecsv(p,rows):
 with p.open("w",encoding="utf-8",newline="") as f:
  fields=list(dict.fromkeys(k for r in rows for k in r));w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sources():
 return {str(p.relative_to(ROOT)).replace("\\","/"):digest(p) for p in list((ROOT/"src/bsd").glob("*.py"))+[ROOT/"experiments/run_gw170814_exploratory_search.py"]}
def main():
 global ACTIVE
 assert git("rev-parse","HEAD")==BASE and not git("diff",BASE,"--","src","app","tests","pyproject.toml")
 rid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
 out=ROOT/"results/gw170814_current_control"/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE=out
 local=ROOT/"data/derived/gw170814_current_control"/rid;local.mkdir(parents=True,exist_ok=False)
 (out/"protocol.md").write_bytes(PROTOCOL.read_bytes());(out/"driver.py").write_bytes(Path(__file__).read_bytes())
 historic=ROOT/"results/gw170814_exploratory_search/results.json"
 m={"run_id":rid,"source_commit":BASE,"source_sha256":sources(),"driver_sha256":digest(Path(__file__)),
  "protocol_sha256":digest(PROTOCOL),"historical_results_sha256":digest(historic),"status":"RUNNING",
  "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
  "scope":"Previously examined O2 offline control, CalB after target; no blind/adaptive/significance claim",
  "dirty_paths":git("status","--porcelain").splitlines(),"sites":[]}
 save(out/"manifest.json",m);print("RUN_DIRECTORY",out,flush=True)
 old._verify_disjoint_intervals()
 prov=json.loads((ROOT/"results/ligo_experiment_003/data_provenance.json").read_text())
 for site,path in old.DETECTORS.items():
  started=time.monotonic();dest=out/site;dest.mkdir()
  entry=next(r for r in prov["files"] if r["event"]=="GW170814" and r["detector"]==site and r["duration_s"]==4096)
  assert digest(path)==entry["sha256"];record=old.load_gwosc_hdf5(path);inj,names=old._load_inj_mask(path)
  intervals=[old.PSD_CORE_INTERVAL,old.CAL_A_CORE_INTERVAL,old.CAL_B_CORE_INTERVAL,old.SEARCH_CORE_INTERVAL]
  seq=[]
  for interval,n in zip(intervals,[60,320,387,160]):
   observations,bad=old._extract_exact_interval_segments(record,interval)
   assert not bad and len(observations)==n,(site,interval,bad);seq.append(observations)
  model=old.fit_background_model(site,seq[0],mode=old.PREPROCESS_MODE,target_sample_rate=1024.,bandpass_hz=old.BANDPASS_HZ)
  det=old.MultiResolutionDetector()
  calobs=[old.preprocess_segment(s,f"{site}_cal_{i:04d}",model) for i,s in enumerate(seq[1]+seq[2])]
  cal=det.calibrate(calobs,target_fap=.01,norm_fraction=.30,calibration_seed=0)
  chrono=np.array([det.analyze(o,cal).max_global_score for o in calobs[cal.n_norm_samples:]])
  audit=old.classify_calibration_status(len(calobs),old.estimate_calibration_uncertainty(chrono,target_fap=.01),target_fap=.01)
  save(dest/"calibration.json",cal.to_summary_dict());save(dest/"null_scores.json",cal.sorted_global_null_scores)
  rows=[];candidates=[]
  print("SITE_START",site,flush=True)
  for i,s in enumerate(seq[3]):
   obs=old.preprocess_segment(s,f"{site}_search_{i:03d}",model);r=det.analyze(obs,cal)
   v=old.apply_vetoes(obs,r.candidates,old.FROZEN_VETO_CONFIG);vm={d.candidate.candidate_id:d for d in v.decisions}
   gps=s.gps_start;on=gps<=old.EVENT_GPS<gps+4
   role="ON_SOURCE" if on else("OFF_SOURCE" if gps+5<=old.EVENT_GPS-16 or gps-1>=old.EVENT_GPS+16 else "GUARD")
   p,_=cal.empirical_p_value(r.max_global_score)
   rows.append({"window_index":i,"gps_start":gps,"role":role,"score":r.max_global_score,"threshold":cal.global_threshold,
    "empirical_p":p,"detected_raw":r.detected,"detected_post_veto":v.detected_after_veto})
   for c in r.candidates:
    d=vm[c.candidate_id]
    candidates.append({**c.to_dict(),"window_index":i,"role":role,"window_gps_start":gps,
     "gps_start":gps+c.t_start,"gps_end":gps+c.t_end,"gps_peak":gps+c.peak_time,"anomaly_score":c.standardized_score,
     "vetoed":d.vetoed,"veto_reasons":d.reason_codes,"flags":old._inspect_flags_at_interval(record,inj,names,gps+c.t_start,gps+c.t_end)})
   if on:
    raw=det._raw_analyze_all_scales(obs)
    sub=[]
    for scale,(score,clusters,diag) in raw.items():
     for c in clusters:
      sub.append({**c.to_dict(),"scale":scale,"standardized_score":cal.scale_stats[scale].standardize(c.anomaly_score),
       "near_reference":c.t_start<=old.EVENT_GPS-gps+.5 and c.t_end>=old.EVENT_GPS-gps-.5 and c.f_min<=250 and c.f_max>=35})
    save(dest/"event_diagnostics.json",{"per_scale_raw_max":r.per_scale_raw_max,"per_scale_std_max":r.per_scale_std_max,
      "winning_scale":r.winning_scale,"all_seed_qualified_raw_clusters":sub,
      "quality":old._inspect_flags_at_interval(record,inj,names,gps,gps+4),"veto":asdict(v)})
    fn=local/f"GW170814_{site}_event_window_whitened.npy";np.save(fn,obs.signal,allow_pickle=False)
    save(dest/"event_array_manifest.json",{"path":str(fn.relative_to(ROOT)).replace("\\","/"),"sha256":digest(fn),"gps_start":gps,"fs":1024.})
  writecsv(dest/"timeline.csv",rows);save(dest/"candidates.json",candidates)
  off=[r for r in rows if r["role"]=="OFF_SOURCE"];event=next(r for r in rows if r["role"]=="ON_SOURCE")
  assert event["window_index"]==80 and len(off)==151
  summary={"site":site,"event":event,"audit":audit.to_dict(),"offsource_windows":len(off),
   "offsource_raw":sum(r["detected_raw"] for r in off),"offsource_post_veto":sum(r["detected_post_veto"] for r in off),"elapsed_s":time.monotonic()-started}
  save(dest/"input_provenance.json",{"original":entry,"model":asdict(model),"injection_shortnames":names,"dq_shortnames":record.dq_shortnames})
  save(dest/"summary.json",summary);m["sites"].append(summary);save(out/"manifest.json",m)
  print("SITE_DONE",site,"detected",event["detected_raw"],"score",event["score"],"threshold",event["threshold"],flush=True)
 pairs=[]
 for role in ("ON_SOURCE","OFF_SOURCE"):
  for a,b in (("H1","L1"),("H1","V1"),("L1","V1")):
   for post in (False,True):
    def select(site):
     return deduplicate_candidates_by_window([c for c in json.loads((out/site/"candidates.json").read_text()) if c["role"]==role and(not post or not c["vetoed"])])
    matches=match_coincident_candidates(select(a),select(b),max_dt_s=.05 if b=="L1" else .065,allow_interval_gap=False)
    pairs.append({"role":role,"pair":a+"-"+b,"post_veto":post,"coarse_matches":len(matches)})
 save(out/"coarse_foreground_matches.json",pairs)
 assert sources()==m["source_sha256"] and digest(historic)==m["historical_results_sha256"]
 m["status"]="COMPLETE";m["source_and_historical_hashes_unchanged"]=True;save(out/"manifest.json",m);print("RUN_STATUS COMPLETE",flush=True)
if __name__=="__main__":
 try:main()
 except Exception as exc:
  if ACTIVE and (ACTIVE/"manifest.json").exists():
   m=json.loads((ACTIVE/"manifest.json").read_text());m["status"]="ERROR";m["error"]={"type":type(exc).__name__,"message":str(exc)};save(ACTIVE/"manifest.json",m)
  raise
