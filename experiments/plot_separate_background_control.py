"""Presentation of all prespecified profiles/groups, no new inference or detector changes."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=Path(sys.argv[1])
v=json.loads((out/"independent_verification.json").read_text())
assert v["status"]=="PASS"
cases=["stationary","level_step","spectral_line","impulsive_noise"]
labels=["Stationary pink","Post-step level (2x)","Pink + 110 Hz line","Pink + rare impulses"]
arms=["stream_contaminated","paired_noise_oracle","separate_matched","separate_old"]
short=["Contaminated stream","Paired noise (oracle)","Separate matched","Separate old regime"]
colors=["#d55e00","#7570b3","#009e73","#999999"]
fig,ax=plt.subplots(1,1,figsize=(10,4.8))
x=np.arange(4);width=.19
for j,arm in enumerate(arms):
 values=[v["groups"][f"{case}_{arm}_fap0.01"]["noise_alarms"] for case in cases]
 bars=ax.bar(x+(j-1.5)*width,values,width,label=short[j],color=colors[j])
 ax.bar_label(bars,padding=3,fontsize=8)
ax.set_xticks(x,labels);ax.set_ylabel("Noise-window nominations / 600")
ax.set_title("Separate-background source control: observed noise alarms (FAP target 0.01)")
ax.set_ylim(0,max(g["noise_alarms"] for k,g in v["groups"].items() if k.endswith('fap0.01'))*1.25+1)
ax.legend(fontsize=8,ncol=2);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.text(.5,.01,"Three calibration seeds per profile; paired arms. Counts are descriptive, not a certified FAR bound.",ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,1));fig.savefig(out/"background_source_noise.png",dpi=180);plt.close(fig)
keys=["chirp_snr8.0","chirp_snr18.0","burst_snr8.0","burst_snr18.0"]
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,case,title in zip(axes.flat,cases,labels):
 values=np.array([[v["groups"][f"{case}_{arm}_fap0.01"]["signals"][k]["associated"] for k in keys] for arm in arms])
 im=ax.imshow(values,vmin=0,vmax=30,cmap='viridis',aspect='auto')
 for i in range(4):
  for j in range(4):ax.text(j,i,f"{values[i,j]}/30",ha='center',va='center',color='white' if values[i,j]<16 else 'black',fontsize=10)
 ax.set_xticks(range(4),["Chirp 8","Chirp 18","Burst 8","Burst 18"]);ax.set_yticks(range(4),short,fontsize=8);ax.set_title(title)
fig.colorbar(im,ax=axes.ravel().tolist(),label='Known-support associations / 30')
fig.suptitle("Same target inputs, four background sources (FAP target 0.01)\nPaired SNR amplitudes; loose support association, not physical identification",fontsize=12)
fig.savefig(out/"background_source_signals.png",dpi=180);plt.close(fig)
print('Saved two descriptive figures without changing measurement files.')
