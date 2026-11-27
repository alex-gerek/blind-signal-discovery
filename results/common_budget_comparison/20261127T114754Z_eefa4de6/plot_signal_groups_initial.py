from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent
results=json.loads((out/"inference.json").read_text())["signal_groups"]
fig,ax=plt.subplots(figsize=(10,4.5),layout="constrained")
x=np.arange(len(results))
for offset,arm,label in [(-.18,"baseline","Original"),(.18,"slow","Approx. 1s median")]:
 values=[100*g[arm]["rate"] for g in results]
 lower=[100*(g[arm]["rate"]-g[arm]["lower95"]) for g in results]
 upper=[100*(g[arm]["upper95"]-g[arm]["rate"]) for g in results]
 bars=ax.bar(x+offset,values,.36,label=label,yerr=[lower,upper],capsize=3)
 for bar,g,value in zip(bars,results,values):ax.text(bar.get_x()+bar.get_width()/2,value+4,f"{g[arm]['count']}/300",ha="center",fontsize=9)
labels=[("Chirp" if g["family"]=="linear_chirp" else "Broadband burst")+f"\nSNR {g['snr']}"+(" (primary)" if g["primary"] else "") for g in results]
ax.set_xticks(x,labels);ax.set_ylim(0,109);ax.set_ylabel("Known-support association [%]");ax.grid(axis="y",alpha=.2);ax.legend()
ax.set_title("Held-out synthetic mixture: 95% per-method intervals; secondary groups descriptive")
target=out/"signal_groups.png"
if target.exists():raise FileExistsError("Preserve existing plot; use a fresh output location")
fig.savefig(target,dpi=150)
