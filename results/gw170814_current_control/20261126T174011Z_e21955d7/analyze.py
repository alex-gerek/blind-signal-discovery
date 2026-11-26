"""Independent verification and small figures for the bounded GW170814 replay."""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import stft
root=Path(sys.argv[1]).resolve();repo=root.parents[2]
def load(p):return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=load(root/"manifest.json");assert m["status"]=="COMPLETE" and len(m["sites"])==3
for f,h in m["source_sha256"].items():assert sha(repo/f)==h
historic=repo/"results/gw170814_exploratory_search/results.json";assert sha(historic)==m["historical_results_sha256"]
old=load(historic)["calibration_and_search_summary"]
checks=[];fields={}
for site in ("H1","L1","V1"):
 p=root/site;s=load(p/"summary.json");c=load(p/"calibration.json");null=np.asarray(load(p/"null_scores.json"))
 with (p/"timeline.csv").open(encoding="utf-8",newline="") as f:rows=list(csv.DictReader(f))
 assert len(rows)==160
 for i,r in enumerate(rows):
  score=float(r["score"]);th=float(r["threshold"]);pv=(1+len(null)-np.searchsorted(null,score,side="left"))/(len(null)+1)
  assert float(r["gps_start"])==1186741540+4*i and th==c["global_threshold"]
  assert float(r["empirical_p"])==pv and (r["detected_raw"]=="True")==(score>th)
 off=[r for r in rows if r["role"]=="OFF_SOURCE"];assert len(off)==151
 assert sum(r["detected_raw"]=="True" for r in off)==s["offsource_raw"]
 assert sum(r["detected_post_veto"]=="True" for r in off)==s["offsource_post_veto"]
 d=load(p/"event_diagnostics.json")
 assert max(d["per_scale_std_max"].values())==s["event"]["score"]
 oldrow=old[site]
 matches={"event_raw":s["event"]["detected_raw"]==oldrow["on_source_detected_raw"],
  "event_post_veto":s["event"]["detected_post_veto"]==oldrow["on_source_detected_post_veto"],
  "score_rounded4":round(s["event"]["score"],4)==oldrow["on_source_max_score"],
  "threshold_rounded4":round(s["event"]["threshold"],4)==oldrow["frozen_threshold"],
  "offsource_raw":s["offsource_raw"]==oldrow["n_off_nominated_raw"],
  "offsource_post_veto":s["offsource_post_veto"]==oldrow["n_off_nominated_post_veto"]}
 checks.append({"site":site,"summary":s,"historical_agreement":matches,"per_scale_std_max":d["per_scale_std_max"],
  "near_reference_seed_clusters":sum(r["near_reference"] for r in d["all_seed_qualified_raw_clusters"]),
  "decision_explanation":"Above frozen threshold" if s["event"]["detected_raw"] else "All per-scale standardized maxima below frozen threshold; pre-veto miss"})
 info=load(p/"event_array_manifest.json");npfile=repo/info["path"];assert sha(npfile)==info["sha256"]
 x=np.load(npfile,allow_pickle=False)
 f,t,z=stft(x,1024.,window="hann",nperseg=128,noverlap=96,boundary=None,padded=False)
 db=10*np.log10(np.maximum(abs(z)**2,1e-24));t=t-1.5
 mask=(f[:,None]>=30)&(f[:,None]<=380)&(t[None,:]>=-.5)&(t[None,:]<=.5)
 fields[site]=(x,f,t,db,db[mask])
fig,axs=plt.subplots(2,3,figsize=(13,6),sharex="col")
vmin,vmax=np.percentile(np.concatenate([v[-1] for v in fields.values()]),[5,99.5])
for col,site in enumerate(("H1","L1","V1")):
 x,f,t,db,_=fields[site];s=checks[col]["summary"]
 axs[0,col].plot(np.arange(4096)/1024.-1.5,x,lw=.6)
 axs[0,col].set_title(f"{site}: {'DETECTED' if s['event']['detected_raw'] else 'MISSED'}; score {s['event']['score']:.2f} / {s['event']['threshold']:.2f}")
 im=axs[1,col].pcolormesh(t,f,db,cmap="viridis",shading="auto",vmin=vmin,vmax=vmax)
 axs[1,col].set_ylim(30,380);axs[1,col].set_xlabel("Seconds from known control time")
 for a in axs[:,col]:a.set_xlim(-.5,.5);a.axvline(0,color="black",ls=":",lw=.7)
 fig.colorbar(im,ax=axs[1,col],label="STFT dB",extend="both")
axs[0,0].set_ylabel("Whitened strain");axs[1,0].set_ylabel("Hz")
fig.suptitle("Real GW170814 control: identical fixed recipe, no threshold tuning")
fig.text(.5,.008,"Display contrast: pooled visible-bin5th/99.5th percentiles; no amplitude rescaling of detection inputs.",ha="center",fontsize=8)
fig.tight_layout(rect=(0,.03,1,.95));assert not (root/"event_comparison.png").exists();fig.savefig(root/"event_comparison.png",dpi=150);plt.close(fig)
record={"status":"VERIFIED","sites":checks,"coarse_foreground_matches":load(root/"coarse_foreground_matches.json"),
 "display_limits":{"vmin_db":float(vmin),"vmax_db":float(vmax),"percentiles":[5,99.5]},
 "scope":"Fixed historical offline control. No time-slide significance, causal adaptation test or physical discovery. Historical original-array hashes unavailable."}
assert not (root/"independent_verification.json").exists();(root/"independent_verification.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
lines=["# GW170814 bounded real-data replay","",f"Source{m['source_commit']}; run{m['run_id']}. Previously used O2 epoch; same historical partition,707calibration/160target windows, fixed detector and reversible vetoes. CalB is after target:offline control, not causal deployment.",
 "","| Site | Event score | Threshold | Raw / post-veto event | Off-source raw / post-veto (N=151) | Historical six checks agree |",
 "| --- | --- | --- | --- | --- | --- |"]
for x in checks:
 s=x["summary"];e=s["event"]
 lines.append(f"| {x['site']} | {e['score']:.8g} | {e['threshold']:.8g} | {e['detected_raw']} / {e['detected_post_veto']} | {s['offsource_raw']} / {s['offsource_post_veto']} | {all(x['historical_agreement'].values())} |")
lines+=["","## Interpretation",
 "- H1/V1 below-threshold decisions are already misses before vetoes; vetoes did not cause those misses. Per-scale maxima and subthreshold seed-qualified clusters are preserved in event_diagnostics.json.",
 "- A numerical threshold gap explains the code decision. Lower instrument response, local noise and STFT morphology are possible physical/statistical causes, not proven causal explanations by this replay.",
 "- Event-window nomination is a control association, not independent confirmation or a new event. Coarse peak matches are saved separately and do not provide astrophysical significance. No time-slide/background significance was run.",
 "- Off-source nominations are real-data candidates, not certified noise-only false alarms. Extreme V1 candidates and all veto reasons are preserved; no novelty claim.",
 "- Compare with GW150914 carefully: its separate protocol used320calibration windows rather than707. These selected previously examined events are not a representative efficiency estimate; do not quote3/5 detector-event recovery as a population rate.",
 "","![Real event comparison](event_comparison.png)","",
 "## Preservation and attribution",
 "All source/historical hashes and full-precision decisions were independently checked; raw data hashes were checked before measurement. Original/processed event arrays stay ignored/local, with exact hashes. Historical array bytes cannot be independently reconstructed from old reports alone.",
 "Data: GWOSC/LIGO Scientific Collaboration/Virgo Collaboration, CC BY4.0. [Event release](https://gwosc.org/events/GW170814/), [acknowledgement/licensing](https://gwosc.org/acknowledgement/). Processing: independent PSD whitening,30–380Hz smooth bandpass and4096->1024Hz antialias decimation; display contrast affects figures only.",
 "No injected waveform, detector tuning, O3/manifest, extra event, main merge, release or Zenodo operation. Detailed scientific article deferred.",""]
assert not (root/"report.md").exists();(root/"report.md").write_text("\n".join(lines),encoding="utf-8")
print("VERIFIED")
for c in checks:print(c["site"],json.dumps(c))
