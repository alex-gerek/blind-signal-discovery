from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import stft
from matplotlib.patches import Rectangle
p=Path(sys.argv[1]).resolve();repo=p.parents[2];m=json.loads((p/"manifest.json").read_text());tc=m["known_event_gps"]
fields={}
for site in ("H1","L1"):
 info=json.loads((p/site/"derived_file_manifest.json").read_text())
 x=np.load(repo/info["target"]["npy"],allow_pickle=False)[39*4096:40*4096]
 f,t,z=stft(x,1024.,window="hann",nperseg=128,noverlap=96,boundary=None,padded=False)
 db=10*np.log10(np.maximum(np.abs(z)**2,1e-24));rel=t+1126259460.-tc
 mask=(f[:,None]>=30)&(f[:,None]<=380)&(rel[None,:]>=-.5)&(rel[None,:]<=.3)
 fields[site]=(x,f,rel,db,db[mask])
values=np.concatenate([v[-1] for v in fields.values()])
lo,hi=np.percentile(values,[5,99.5])
fig,axs=plt.subplots(2,2,figsize=(11,6),sharex="col")
for col,site in enumerate(("H1","L1")):
 x,f,t,db,_=fields[site]
 axs[0,col].plot(np.arange(4096)/1024.+1126259460.-tc,x,lw=.65)
 axs[0,col].set_title(site+" real whitened measurement")
 axs[0,col].set_ylabel("Dimensionless whitened strain")
 image=axs[1,col].pcolormesh(t,f,db,shading="auto",cmap="viridis",vmin=lo,vmax=hi)
 for c in json.loads((p/site/"frozen_target_candidates.json").read_text()):
  if c["window_index"]==39:
   axs[1,col].add_patch(Rectangle((c["gps_t_start"]-tc,c["f_min"]),c["gps_t_end"]-c["gps_t_start"],
    c["f_max"]-c["f_min"],fill=False,edgecolor="#ef5350",lw=.8))
 axs[1,col].set_ylim(30,380);axs[1,col].set_ylabel("Hz");axs[1,col].set_xlabel("Seconds relative to known event")
 for a in axs[:,col]:a.set_xlim(-.5,.3);a.axvline(0,color="black",ls=":",lw=.7)
 fig.colorbar(image,ax=axs[1,col],label="STFT dB re squared input unit",extend="both")
fig.suptitle("GW150914: real data, no injected waveform; shared display contrast")
fig.text(.5,.004,"Display-only limits: pooled visible-bin 5th/99.5th percentiles. Detection used unmodified full windows.",ha="center",fontsize=8)
fig.tight_layout(rect=(0,.025,1,.95))
out=p/"event_waveform_spectrogram_contrast.png";assert not out.exists();fig.savefig(out,dpi=160);plt.close(fig)
out=p/"display_contrast.json";assert not out.exists();out.write_text(json.dumps({"method":"Pooled visible-bin fifth/99.5th percentiles, both sites; display only",
 "vmin_db":float(lo),"vmax_db":float(hi),"time_relative_to_known_event_s":[-.5,.3],"frequency_display_hz":[30,380],
 "original_full_range_plot_preserved":"event_waveform_spectrogram.png"},indent=2),encoding="utf-8")
print(lo,hi)
