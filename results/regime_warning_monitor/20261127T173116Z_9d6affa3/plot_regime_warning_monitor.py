"""Descriptive fixed-monitor plots; warning-only, no new scientific inference."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
v=load(out/'independent_verification.json');assert v['status']=='PASS'
seeds=list(range(20101,20106));cases=['stable','sparse_transients','repeated_transients']
fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');x=np.arange(5);width=.25
for j,(case,color) in enumerate(zip(cases,['#666666','#0072b2','#e69f00'])):
 values=[next(c['summary']['warning_windows'] for c in v['case_summaries'] if c['seed']==s and c['case']==case) for s in seeds]
 bars=ax.bar(x+(j-1)*width,values,width,label=case.replace('_',' '),color=color);ax.bar_label(bars,padding=3)
ax.set_xticks(x,[str(s) for s in seeds]);ax.set_ylabel('Flagged windows / 505 eligible windows');ax.set_title('Default monitor: stable noise and true transient signals\nFive paired traces; 4 s/window; flags are not independent alarms')
maxval=max(c['summary']['warning_windows'] for c in v['case_summaries'] if c['case'] in cases);ax.set_ylim(0,maxval*1.3+1);ax.legend();ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.savefig(out/'warning_counts.png',dpi=180);plt.close(fig)
a=load(out/'20101/spectral_line/monitor_assessments.json');b=load(out/'20101/persistent_signal/monitor_assessments.json');x=np.arange(284,800)
fig,axes=plt.subplots(3,1,figsize=(10,7),sharex=True,layout='constrained')
axes[0].plot(x,[v['max_feature_z_drift'] for v in a],label='Noise-line label',color='#0072b2');axes[0].plot(x,[v['max_feature_z_drift'] for v in b],label='Desired-signal label',color='#d55e00',linestyle='--');axes[0].set_ylabel('Feature drift Z');axes[0].legend()
for ax,data,title,color in [(axes[1],a,'Same waveform labelled as a noise change','#0072b2'),(axes[2],b,'Same waveform labelled as a desired persistent signal','#d55e00')]:
 ax.step(x,[int(v['regime_shifted']) for v in data],where='post',color=color);ax.set_ylim(-.1,1.2);ax.set_yticks([0,1]);ax.set_ylabel('Warning');ax.set_title(title,fontsize=10)
for ax in axes:ax.axvline(360,color='#333333',linestyle=':',linewidth=1);ax.grid(alpha=.15)
axes[-1].set_xlabel('Global window index (4 s/window); tone begins at 360')
fig.suptitle('Identical input -> identical warning; physical cause is not identified\nRepresentative seed20101; labels never enter the monitor',fontsize=12)
fig.savefig(out/'signal_noise_ambiguity.png',dpi=180);plt.close(fig)
print('Saved warning-count and signal/noise-ambiguity figures.')
