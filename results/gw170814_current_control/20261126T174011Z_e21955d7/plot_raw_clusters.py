from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.signal import stft
p=Path(sys.argv[1]).resolve();repo=p.parents[2]
v=json.loads((p/"independent_verification.json").read_text())["display_limits"]
fig,axes=plt.subplots(1,3,figsize=(12,3.8))
for ax,site in zip(axes,("H1","L1","V1")):
 d=json.loads((p/site/"event_diagnostics.json").read_text());m=json.loads((p/site/"event_array_manifest.json").read_text());s=json.loads((p/site/"summary.json").read_text())
 x=np.load(repo/m["path"],allow_pickle=False)
 f,t,z=stft(x,1024.,nperseg=128,noverlap=96,boundary=None,padded=False)
 ax.pcolormesh(t-1.5,f,10*np.log10(np.maximum(abs(z)**2,1e-24)),shading="auto",cmap="viridis",vmin=v["vmin_db"],vmax=v["vmax_db"])
 for scale in d["per_scale_std_max"]:
  candidates=[c for c in d["all_seed_qualified_raw_clusters"] if c["scale"]==scale and c["near_reference"]]
  if candidates:
   c=max(candidates,key=lambda c:c["standardized_score"]);accepted=c["standardized_score"]>s["event"]["threshold"]
   ax.add_patch(Rectangle((c["t_start"]-1.5,c["f_min"]),c["t_end"]-c["t_start"],c["f_max"]-c["f_min"],
    fill=False,edgecolor="red" if accepted else "#dddddd",lw=1.1,ls="-" if accepted else "--"))
 ax.set_xlim(-.5,.5);ax.set_ylim(30,380);ax.set_title(site);ax.set_xlabel("Seconds from reference")
axes[0].set_ylabel("Hz")
fig.suptitle("Highest near-reference raw cluster per scale: red above threshold; gray below")
fig.text(.5,.01,"Diagnostic clusters can be noise; boxes alone do not identify the gravitational wave.",ha="center",fontsize=9)
fig.tight_layout(rect=(0,.05,1,.91))
q=p/"event_raw_cluster_diagnostics.png";assert not q.exists();fig.savefig(q,dpi=150);plt.close(fig)
