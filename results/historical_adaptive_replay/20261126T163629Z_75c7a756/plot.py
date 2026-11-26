import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
p=Path(sys.argv[1])
w=json.loads((p/"abrupt_step_500win/windows.json").read_text())
x=np.arange(len(w));score=[r["stft_max_score"] for r in w]
threshold=np.array([r["stft_threshold"] if r["stft_threshold"] is not None else np.nan for r in w])
stft=np.array([r["diagnostics"]["dedicated_stft_available"] for r in w])
union=np.array([r["diagnostics"]["union_mode_status"]=="FULL_UNION_AVAILABLE" for r in w])
bb=np.array([r["broadband_nominated"] for r in w])
fig,axes=plt.subplots(3,1,figsize=(11,8),sharex=True,gridspec_kw={"height_ratios":[2,1,1]})
axes[0].plot(x,score,label="STFT score",color="#2364aa",lw=.8)
axes[0].plot(x,threshold,label="Used threshold",color="#c62828",lw=1.6)
for start in [162]:
 axes[0].axvspan(start,234,color="gray",alpha=.2,label="STFT unavailable")
axes[0].set_ylabel("Dimensionless score")
axes[0].set_title("Historical adaptive recipe replay: pink noise to shifted AR(2), source 043f0eb")
axes[0].legend(loc="upper left",ncol=3,fontsize=8)
axes[1].step(x,stft.astype(float),where="post",label="Dedicated STFT available")
axes[1].step(x,union.astype(float),where="post",label="Full Union available",ls="--")
axes[1].set_ylim(-.1,1.2);axes[1].set_yticks([0,1]);axes[1].set_ylabel("Available")
axes[1].legend(loc="lower right",fontsize=8)
axes[2].step(x,bb.astype(float),where="post",label="Frozen broadband nomination",color="#7b1fa2",lw=.8)
for i in range(330,430,10):
 hit=w[i]["detection_available"] and w[i]["stft_nominated_post_veto"]
 axes[2].scatter(i,-.12 if hit else -.30,marker="^" if hit else "x",color="green" if hit else "red",s=35)
axes[2].set_ylim(-.45,1.15);axes[2].set_yticks([0,1]);axes[2].set_ylabel("BB flag")
axes[2].set_xlabel("Window index (4 seconds/window); green triangles: detected fixed injections; red x: missed")
axes[2].legend(loc="upper left",fontsize=8)
for ax in axes:
 ax.axvline(100,color="black",ls=":",lw=1);ax.axvline(367,color="#c62828",ls=":",lw=1)
 ax.grid(alpha=.2)
fig.tight_layout()
dest=p/"adaptive_recovery.png"
assert not dest.exists();fig.savefig(dest,dpi=160);plt.close(fig)
print(dest)
