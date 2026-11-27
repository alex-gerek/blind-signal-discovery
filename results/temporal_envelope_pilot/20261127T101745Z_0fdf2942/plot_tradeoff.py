from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent
summaries=json.loads((out/"summaries.json").read_text())
synth=[s for s in summaries if s["case"].startswith("synthetic")]
regimes=("stationary","smooth_drift","persistent_step")
fig,axes=plt.subplots(1,2,figsize=(11,4),layout="constrained")
for ax,role,den,title in [(axes[0],"NOISE",600,"Noise-only false-alarm windows / 600"),(axes[1],"SIGNAL_band_limited_noise_burst_SNR8",30,"Weak broadband-burst associations / 30")]:
 for offset,arm,label in [(-.18,"baseline","Original"),(.18,"slow","Approx. 1s median")]:
  key="window_nominations" if role=="NOISE" else "associated_nominations"
  values=[sum(s["by_regime"][regime]["arms"][arm][role][key] for s in synth) for regime in regimes]
  bars=ax.bar(np.arange(3)+offset,values,.36,label=label)
  for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value+.2,f"{value}/{den}",ha="center",fontsize=9)
 ax.set_xticks(np.arange(3),["Stationary","Smooth drift","Level step"]);ax.set_ylim(0,max(values+([15] if role=="NOISE" else [8]))*1.25+1);ax.set_title(title);ax.grid(axis="y",alpha=.2)
axes[0].legend(fontsize=9)
fig.suptitle("Fixed pilot: observed gain and cost; paired regimes, no significance claim")
target=out/"tradeoff_counts.png"
if target.exists():raise FileExistsError("Preserve the earlier plot; choose a fresh output directory")
fig.savefig(target,dpi=150)
