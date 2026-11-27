"""Inspect existing GW170814 records; no detector execution or retuning."""
from __future__ import annotations
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,os,subprocess,sys,uuid
ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/"results/gw170814_current_control/20261126T174011Z_e21955d7"
PROTOCOL=ROOT/"docs/gw170814_miss_analysis_protocol_2026-11-27.md"
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ALLOWED={str((ROOT/load(PARENT/s/"event_array_manifest.json")["path"]).resolve()).lower().replace("\\","/") for s in ("H1","L1","V1")}
DATA=str((ROOT/"data").resolve()).lower().replace("\\","/")+"/"
def guard(event,args):
 if event=="open" and args and isinstance(args[0],(str,bytes,os.PathLike)):
  text=os.path.abspath(os.fsdecode(args[0])).lower().replace("\\","/")
  if "o3_holdout" in text or (text.startswith(DATA) and text not in ALLOWED):raise PermissionError("Data/O3 access outside saved-array allowlist")
 if event in ("socket.connect","socket.getaddrinfo"):raise PermissionError("Network blocked")
sys.addaudithook(guard)
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.signal import stft
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding="utf-8")
def git(*a):return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
def main():
 assert git("rev-parse","HEAD")=="30849e8189c6c9b643492c085a3594fb1280cbf0"
 assert not git("diff","HEAD","--","src","app","tests","pyproject.toml")
 rid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
 out=ROOT/"results/gw170814_miss_analysis"/rid;out.mkdir(parents=True,exist_ok=False)
 (out/"protocol.md").write_bytes(PROTOCOL.read_bytes());(out/"driver.py").write_bytes(Path(__file__).read_bytes())
 files=[PARENT/"independent_verification.json"]
 for site in ("H1","L1","V1"):
  files.extend(PARENT/site/name for name in ("event_diagnostics.json","calibration.json","input_provenance.json","summary.json","event_array_manifest.json"))
 hashes={str(p.relative_to(ROOT)).replace("\\","/"):digest(p) for p in files}
 m={"run_id":rid,"source_commit":git("rev-parse","HEAD"),"parent":str(PARENT.relative_to(ROOT)).replace("\\","/"),
  "status":"RUNNING","input_sha256":hashes,"protocol_sha256":digest(PROTOCOL),"driver_sha256":digest(Path(__file__)),
  "scope":"Saved-record diagnostic analysis only; no detector/calibration runs, no new detection criterion"}
 save(out/"manifest.json",m);print("RUN_DIRECTORY",out,flush=True)
 sites=[];table=[];spectra={};views={};checked=0
 for site in ("H1","L1","V1"):
  d=load(PARENT/site/"event_diagnostics.json");cal=load(PARENT/site/"calibration.json")
  summary=load(PARENT/site/"summary.json");provenance=load(PARENT/site/"input_provenance.json")
  info=load(PARENT/site/"event_array_manifest.json");npy=ROOT/info["path"];assert digest(npy)==info["sha256"]
  signal=np.load(npy,allow_pickle=False)
  f,t,z=stft(signal,1024.,window="hann",nperseg=128,noverlap=96,boundary=None,padded=False)
  db=10*np.log10(np.maximum(abs(z)**2,1e-24))
  views[site]=(f,t,db)
  stats=cal["scale_stats"];rows=[]
  for scale in d["per_scale_std_max"]:
   s=stats[scale];clusters=[c for c in d["all_seed_qualified_raw_clusters"] if c["scale"]==scale]
   for c in clusters:
    calculated=(c["anomaly_score"]-s["median_null_score"])/s["spread_scale"]
    assert abs(calculated-c["standardized_score"])<1e-12;checked+=1
   near=[c for c in clusters if c["near_reference"]]
   best=max(clusters,key=lambda c:c["standardized_score"])
   nearbest=max(near,key=lambda c:c["standardized_score"]) if near else None
   assert best["standardized_score"]==d["per_scale_std_max"][scale]
   row={"site":site,"scale":scale,"background_median_raw":s["median_null_score"],"background_q95_raw":s["q95_null_score"],
    "background_spread_raw":s["spread_scale"],"global_threshold":cal["global_threshold"],
    "raw_score_boundary":s["median_null_score"]+cal["global_threshold"]*s["spread_scale"],
    "all_cluster_count":len(clusters),"near_reference_cluster_count":len(near),
    "best_anywhere":best,"best_near_reference":nearbest}
   rows.append(row)
   table.append({k:row[k] for k in ("site","scale","background_median_raw","background_spread_raw","global_threshold","raw_score_boundary")}|
    {"raw_observed_anywhere":best["anomaly_score"],"z_anywhere":best["standardized_score"],
     "z_near_reference":nearbest["standardized_score"] if nearbest else None})
  winner=max((r["best_anywhere"] for r in rows),key=lambda c:c["standardized_score"])
  nearwinner=max((r["best_near_reference"] for r in rows if r["best_near_reference"]),key=lambda c:c["standardized_score"])
  assert winner["standardized_score"]==summary["event"]["score"]
  event={"site":site,"event_summary":summary["event"],"global_winner":winner,"near_reference_winner":nearwinner,
   "global_winner_peak_offset_s":winner["peak_time"]-1.5,"rows":rows,"quality":d["quality"]}
  sites.append(event)
  model=provenance["model"];assert model["scale_factor"]==1e18
  freq=np.asarray(model["psd_freqs"]);asd=np.sqrt(np.asarray(model["psd_values"]))/model["scale_factor"]
  band=(freq>=35)&(freq<=250)
  spectra[site]={"frequency":freq,"asd":asd}
  event["background_asd_summary"]={"band_hz":[35,250],"units":"strain/sqrt(Hz)","median_over_frequency":float(np.median(asd[band])),
   "q10_over_frequency":float(np.quantile(asd[band],.1)),"q90_over_frequency":float(np.quantile(asd[band],.9)),
   "interpretation":"Spectral descriptors of saved PSD block; not event SNR or a confidence interval"}
 with (out/"scale_decomposition.csv").open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
 save(out/"diagnostic_results.json",{"sites":sites,"raw_clusters_independently_checked":checked,
  "near_reference_definition":"Parent Boolean: interval intersects[1,2]s and35-250Hz; not physical waveform identity",
  "note":"No regional detection threshold/p-value recalibration or changed primary decisions"})
 fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True)
 for ax,site in zip(axes,sites):
  rows=site["rows"];x=np.arange(3)
  ax.bar(x-.18,[r["best_anywhere"]["standardized_score"] for r in rows],.36,label="Maximum anywhere")
  ax.bar(x+.18,[r["best_near_reference"]["standardized_score"] for r in rows],.36,label="Near reference")
  ax.axhline(site["event_summary"]["threshold"],color="red",ls="--",label="Frozen threshold")
  ax.set_xticks(x,["31ms","125ms","500ms"]);ax.set_title(site["site"]);ax.grid(axis="y",alpha=.2)
 axes[0].set_ylabel("Standardized cluster score");axes[-1].legend(fontsize=8)
 fig.suptitle("GW170814: window maximum versus reference-associated clusters")
 fig.tight_layout();fig.savefig(out/"scale_scores.png",dpi=150);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True)
 visible=np.concatenate([v[2][(v[0]>=30)&(v[0]<=400)].ravel() for v in views.values()])
 lo,hi=np.percentile(visible,[5,99.5])
 for ax,site in zip(axes,sites):
  f,t,db=views[site["site"]]
  ax.pcolormesh(t,f,db,shading="auto",cmap="viridis",vmin=lo,vmax=hi)
  for c,color,label in [(site["global_winner"],"red","Global maximum"),(site["near_reference_winner"],"white","Near-reference maximum")]:
   ax.add_patch(Rectangle((c["t_start"],c["f_min"]),c["t_end"]-c["t_start"],c["f_max"]-c["f_min"],
    fill=False,edgecolor=color,lw=1.2,label=label))
  ax.axvspan(1.,2.,alpha=.12,color="white");ax.axvline(1.5,color="black",ls=":")
  ax.set_xlim(0,4);ax.set_ylim(30,400);ax.set_title(site["site"]);ax.set_xlabel("Seconds within event window")
 axes[0].set_ylabel("Hz");axes[-1].legend(fontsize=7,loc="upper left")
 fig.suptitle("Full four-second window: global winner can lie away from the event")
 fig.tight_layout();fig.savefig(out/"full_window_winners.png",dpi=150);plt.close(fig)
 fig,ax=plt.subplots(figsize=(9,4))
 for site,s in spectra.items():ax.semilogy(s["frequency"],s["asd"],label=site,lw=.8)
 ax.set_xlim(30,380);ax.set_xlabel("Hz");ax.set_ylabel("Saved PSD-block ASD [strain/sqrt(Hz)]");ax.legend();ax.grid(alpha=.2)
 ax.set_title("Background spectra fitted before the target; not event SNR")
 fig.tight_layout();fig.savefig(out/"background_asd.png",dpi=150);plt.close(fig)
 assert all(digest(ROOT/f)==h for f,h in hashes.items())
 m["status"]="COMPLETE";m["input_hashes_unchanged"]=True;m["raw_clusters_checked"]=checked;save(out/"manifest.json",m)
 print("COMPLETE",checked,"cluster standardizations checked")
 for s in sites:
  print(s["site"],"global",s["global_winner"]["standardized_score"],"offset",s["global_winner_peak_offset_s"],
   "near",s["near_reference_winner"]["standardized_score"],"ASDmedian",s["background_asd_summary"]["median_over_frequency"])
if __name__=="__main__":main()
