"""Independent GW150914 output accounting and measured-data figures."""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import stft
from matplotlib.patches import Rectangle
root=Path(sys.argv[1]).resolve();repo=root.parents[2]
def load(p):return json.loads(p.read_text(encoding="utf-8"),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=load(root/"manifest.json");assert m["status"]=="COMPLETE"
for f,h in m["source_sha256"].items():assert sha(repo/f)==h
checks=[];arrays={}
for site in ("H1","L1"):
 p=root/site
 provenance=load(p/"input_provenance.json")
 assert sha(repo/provenance["original"]["relative_local_path"])==provenance["original"]["sha256"]
 files=load(p/"derived_file_manifest.json")
 for name,meta in files.items():
  npy=repo/meta["npy"];cp=repo/meta["csv"]
  assert sha(npy)==meta["npy_sha256"] and sha(cp)==meta["csv_sha256"]
  a=np.load(npy,allow_pickle=False);data=np.loadtxt(cp,delimiter=",",skiprows=1)
  assert np.array_equal(data[:,1],a),(site,name,"CSV float roundtrip")
  assert np.array_equal(data[:,0],m["blocks"][name][0]+np.arange(len(a))/1024.)
  if name=="target":arrays[site]=a
 null=np.asarray(load(p/"frozen_null_scores.json"));cal=load(p/"frozen_calibration.json")
 for block in ("target","control"):
  with (p/f"frozen_{block}_timeline.csv").open(encoding="utf-8",newline="") as f:rows=list(csv.DictReader(f))
  assert len(rows)==(80 if block=="target" else 320)
  for i,r in enumerate(rows):
   score=float(r["score"]);threshold=float(r["threshold"]);pv=(1+len(null)-np.searchsorted(null,score,side="left"))/(len(null)+1)
   assert threshold==cal["global_threshold"] and float(r["empirical_p"])==pv
   assert (r["detected"]=="True")==(score>threshold)
   if r["detected"]=="True":assert pv<=.01
   assert float(r["gps_start"])==m["blocks"][block][0]+4*i
  candidates=load(p/f"frozen_{block}_candidates.json")
  for c in candidates:
   assert c["gps_peak"]==c["window_gps_start"]+c["peak_time"]
 with (p/"workstation_target/timeline.csv").open(encoding="utf-8",newline="") as f:adaptive=list(csv.DictReader(f))
 assert len(adaptive)==80
 for r in adaptive:
  av=r["detection_available"]=="True";score=float(r["active_score_full_precision"])
  threshold=float(r["active_threshold_full_precision"]) if r["active_threshold_full_precision"] else None
  assert (r["active_mode_anomaly_flag"]=="True")==bool(av and threshold is not None and score>threshold)
 frozen=load(p/"summary.json")
 cand=load(p/"frozen_target_candidates.json")
 associated=[c for c in cand if c["window_index"]==39 and c["gps_t_start"]<=m["known_event_gps"]+.5 and c["gps_t_end"]>=m["known_event_gps"]-.5 and c["f_min"]<=250 and c["f_max"]>=35]
 with (p/"quality_windows.csv").open(encoding="utf-8",newline="") as f:quality=list(csv.DictReader(f))
 event_quality=next(q for q in quality if q["block"]=="target" and q["window_index"]=="39")
 top=max(associated,key=lambda c:c["standardized_score"]) if associated else None
 checks.append({"site":site,"strict_json_csv_decisions":"PASS","derived_data_float_and_gps_roundtrip":"PASS",
  "raw_source_hash":"MATCH","summary":frozen,"event_quality":event_quality,"strongest_reference_associated_candidate":top})

# Full target timeline: preserve every candidate and non-detection.
fig,axs=plt.subplots(2,1,figsize=(11,6),sharex=True)
for ax,site in zip(axs,("H1","L1")):
 with (root/site/"frozen_target_timeline.csv").open(encoding="utf-8",newline="") as f:r=list(csv.DictReader(f))
 with (root/site/"workstation_target/timeline.csv").open(encoding="utf-8",newline="") as f:a=list(csv.DictReader(f))
 x=np.arange(80)
 ax.plot(x,[float(v["score"]) for v in r],label="Frozen STFT score",lw=1)
 ax.axhline(float(r[0]["threshold"]),color="red",ls="--",label="Frozen threshold")
 ax.plot(x,[float(v["active_threshold_full_precision"]) if v["active_threshold_full_precision"] else np.nan for v in a],
  color="orange",ls=":",label="Adaptive used threshold")
 for i,v in enumerate(a):
  if v["detection_available"]!="True":ax.axvspan(i-.5,i+.5,color="gray",alpha=.2)
 ax.axvline(39,color="black",ls=":",label="Known event window39")
 ax.set_ylabel(site+" score");ax.grid(alpha=.2);ax.legend(fontsize=8)
axs[-1].set_xlabel("Target window index (4 s); target starts GPS1126259304")
fig.suptitle("Real GW150914 O1 data: fixed preprocessing, no injected waveform")
fig.tight_layout();dest=root/"target_timeline.png";assert not dest.exists();fig.savefig(dest,dpi=160);plt.close(fig)

# Prespecified visualization around the known control time; no detector-window change.
fig,axs=plt.subplots(2,2,figsize=(11,6),sharex="col")
tc=m["known_event_gps"]
for col,site in enumerate(("H1","L1")):
 signal=arrays[site][39*4096:40*4096];t=np.arange(4096)/1024.+1126259460.-tc
 axs[0,col].plot(t,signal,lw=.6);axs[0,col].set_title(site+" measured, whitened/bandpassed")
 axs[0,col].set_ylabel("Whitened strain (dimensionless)")
 f,tt,z=stft(signal,1024.,window="hann",nperseg=128,noverlap=96,boundary=None,padded=False)
 power=10*np.log10(np.maximum(np.abs(z)**2,1e-24))
 im=axs[1,col].pcolormesh(tt+1126259460.-tc,f,power,shading="auto",cmap="viridis")
 for c in load(root/site/"frozen_target_candidates.json"):
  if c["window_index"]==39:
   axs[1,col].add_patch(Rectangle((c["gps_t_start"]-tc,c["f_min"]),c["gps_t_end"]-c["gps_t_start"],
    c["f_max"]-c["f_min"],fill=False,edgecolor="red",lw=.7))
 axs[1,col].set_ylim(30,380);axs[1,col].set_ylabel("Hz");axs[1,col].set_xlabel("Seconds relative to known event GPS")
 for ax in axs[:,col]:ax.set_xlim(-.5,.3);ax.axvline(0,color="white" if ax is axs[1,col] else "black",ls=":",lw=.7)
 fig.colorbar(im,ax=axs[1,col],label="Relative STFT power (dB re squared input unit)")
fig.suptitle("GW150914 control view: real data; red boxes are detector candidates")
fig.tight_layout();dest=root/"event_waveform_spectrogram.png";assert not dest.exists();fig.savefig(dest,dpi=160);plt.close(fig)

record={"status":"VERIFIED","source_hashes":"UNCHANGED","sites":checks,
 "scope":"Known previously used real O1 event/control; no synthetic injection, no independent astrophysical significance",
 "caveats":["Adjacent padded preprocessing windows share context","Empirical p is calibrated rank, not astrophysical probability",
 "Background nominations may contain instrumental or physical transients; not certified noise-only false alarms",
 "Raw/derived waveform files are local outside git; hashes and official provenance are saved",
 "Prepared CSV ingestion in the live UI has not been validated; reserved-manifest guard review remains pending"]}
assert not (root/"independent_verification.json").exists()
(root/"independent_verification.json").write_text(json.dumps(record,indent=2,allow_nan=False),encoding="utf-8")
lines=["# GW150914 current-code real-data control","",
 f"Run `{m['run_id']}`; source `{m['source_commit']}`. Previously used O1 epoch; no injected waveforms, tuning or blind-discovery claim.",
 "", "## Actual outcomes","",
 "| Site | Frozen event score / threshold | Event window nominated | Reference-associated candidates | Adaptive event detected / available | Background-control nominations | Target off-source nominations |",
 "| --- | --- | --- | --- | --- | --- | --- |"]
for r in checks:
 s=r["summary"];e=s["frozen_event_window"];a=s["adaptive_event"]
 lines.append(f"| {r['site']} | {e['score']:.8g} / {e['threshold']:.8g} | {e['detected']} | {s['frozen_reference_associated_candidates']} | {a['detected']} / {a['available']} | {s['frozen_control_nominations']}/{s['frozen_control_windows']} | {s['frozen_offsource_nominations']}/{s['frozen_offsource_windows']} |")
lines+=["","The background fractions are nomination fractions on real data, not an established population FAR. The two sites are not pooled as independent trials. No coincidence significance or new event is inferred.",
 "", "## Measured-data figures","",
 "![Full target timeline](target_timeline.png)","",
 "![Known-event waveform and spectrogram](event_waveform_spectrogram.png)","",
 "The event figure displays a fixed[-0.5,+0.3]s region around the known reference for inspection; detection used the full four-second window. Display30–380Hz is not a hard detector-bin mask. Preprocessing is offline and zero-phase with1s padding; residual edge effects and neighboring shared context remain possible.",
 "", "## Reproducibility and quality","",
 "- Frozen PSD fit/calibration/target/control block footprints are disjoint, including padding; complete prescribed block counts were required.",
 "- Raw file hashes match stored official GWOSC provenance. Source hashes are unchanged. Exact local CSV/NPY hashes are saved; CSV values and GPS times round-trip exactly.",
 "- All frozen score/threshold/empirical-p decisions and adaptive available decisions were independently recomputed from saved files. Candidates retain GPS timing and scale metadata.",
 "- DATA flags were required; all other DQ and hardware-injection bit masks/definitions are preserved in per-site provenance and quality_windows.csv. A data-quality pass does not prove astrophysical origin.",
 "- Direct Workstation target JSON/CSV are saved. Prepared local data files are already whitened; do not preprocess twice. Live UI file ingestion is not claimed verified under the reserved-manifest constraint.",
 "- Attempt1 failed in new-driver calibration serialization before target analysis; it remains preserved separately. The rerun changed API calls only.",
 "", "## Interpretation","",
 "Use event-window nomination, reference-associated candidate morphology/timing and H1/L1 comparison to assess this known positive control. A successful control establishes this configuration's recovery on these selected files, not universal sensitivity or confirmation independent of the published event.",
 "If frozen/adaptive disagree or background candidates appear, preserve them for inspection rather than lowering the threshold or labelling every candidate a discovery.",
 "", "## Data attribution","",
 "Measured strain supplied by the Gravitational Wave Open Science Center and the LIGO Scientific Collaboration/Virgo Collaboration; GWOSC currently serves the LIGO-Virgo-KAGRA community. Original data:CC BY4.0. Transformations: independent PSD whitening,30–380Hz bandpass and4096->1024Hz antialiased decimation; figures and local derivative files were generated for this control.",
 "Sources: [GW150914 release](https://gwosc.org/events/GW150914/), [GWOSC data acknowledgement/licensing](https://gwosc.org/acknowledgement/), R.Abbott et al., [O1/O2 open data, SoftwareX13(2021)100658](https://doi.org/10.1016/j.softx.2021.100658). The project's MIT license does not replace the source data's CC BY terms.",
 "", "No O3 input/manifest or non-allowlisted data were accessed; network was blocked during measurement. Original and large derived waveform files remain outside git.",""]
assert not (root/"report.md").exists();(root/"report.md").write_text("\n".join(lines),encoding="utf-8")
print("INDEPENDENT_VERIFICATION_PASS")
for r in checks:print(r["site"],json.dumps({"summary":r["summary"],"event_quality":r["event_quality"],"strongest_candidate":r["strongest_reference_associated_candidate"]}))
