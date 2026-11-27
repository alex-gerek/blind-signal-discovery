"""Descriptive source-control figures; all prespecified profiles and signal groups."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);v=json.loads((out/"independent_verification.json").read_text())
assert v["status"]=="PASS"
cases=["stationary","level_step","spectral_line","impulsive_noise"]
titles=["Stationary pink","Post-step level (2x)","Pink + 110 Hz line","Pink + rare impulses"]
arms=["stream_contaminated","paired_noise_oracle","separate_matched","separate_old"]
short=["Contaminated","Paired oracle","Separate matched","Separate old"]
colors=["#d55e00","#7570b3","#009e73","#999999"]
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,case,title in zip(axes.flat,cases,titles):
 values=[v["groups"][f"{case}_{arm}_fap0.01"]["noise_alarms"] for arm in arms]
 bars=ax.bar(np.arange(4),values,color=colors);ax.bar_label(bars,padding=3)
 ax.set_xticks(range(4),short,rotation=12,fontsize=8);ax.set_ylabel("Nominations / 600 noise windows")
 ax.set_ylim(0,max(values)*1.2+1);ax.set_title(title);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.suptitle("Observed noise alarms with four calibration sources (FAP target 0.01)\nPanel y-ranges differ. Three calibration seeds; descriptive paired comparisons.",fontsize=12)
fig.savefig(out/"background_source_noise.png",dpi=180);plt.close(fig)
keys=["chirp_snr8.0","chirp_snr18.0","burst_snr8.0","burst_snr18.0"]
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,case,title in zip(axes.flat,cases,titles):
 values=np.array([[v["groups"][f"{case}_{arm}_fap0.01"]["signals"][k]["associated"] for k in keys] for arm in arms])
 im=ax.imshow(values,vmin=0,vmax=30,cmap='viridis',aspect='auto')
 for i in range(4):
  for j in range(4):ax.text(j,i,f"{values[i,j]}/30",ha='center',va='center',color='white' if values[i,j]<16 else 'black',fontsize=10)
 ax.set_xticks(range(4),["Chirp\nSNR 8","Chirp\nSNR 18","Burst\nSNR 8","Burst\nSNR 18"],fontsize=9);ax.set_yticks(range(4),short,fontsize=8);ax.set_title(title)
fig.colorbar(im,ax=axes.ravel().tolist(),label='Known-support associations / 30')
fig.suptitle("Identical targets, different calibration sources (FAP target 0.01)\nRead with the noise-alarm figure: association is not physical identification.",fontsize=12)
fig.savefig(out/"background_source_signals.png",dpi=180);plt.close(fig)
print('Saved two descriptive source-control figures.')
