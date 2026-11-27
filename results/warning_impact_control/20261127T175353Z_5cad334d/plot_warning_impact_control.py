"""Descriptive warning/decision-impact and injection/noise comparison plots."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
assert load(out/'independent_verification.json')['status']=='PASS'
rows=load(out/'aggregate_summary.json');impact=load(out/'independent_impact_analysis.json')
cases=['stable','level_step','spectral_line'];labels=['Stable pink','2x level step','110 Hz line']
fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for ax,key,title,ylabel,color in [(axes[0],'warning_windows','Monitor warnings','Flagged windows / 2525 eligible','#0072b2'),(axes[1],'noise_alarms','Frozen search noise nominations','Nominations / 2580 noise windows','#d55e00')]:
 vals=[next(r[key] for r in rows if r['case']==case) for case in cases]
 bars=ax.bar(range(3),vals,color=color);ax.bar_label(bars,padding=3);ax.set_ylim(0,max(vals)*1.2+1);ax.set_xticks(range(3),labels,fontsize=9);ax.set_title(title);ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.suptitle('Same generic warning, different impact on a normalized search statistic\nFive paired traces; descriptive counts, not independent-window inference',fontsize=11)
fig.savefig(out/'warning_search_impact.png',dpi=180);plt.close(fig)
support=impact['summary']['injection_support'];cases=['sparse_transients','step_with_probes','line_with_probes'];x=np.arange(3);width=.25
fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
for j,key,color,label in [(0,'positive_association','#0072b2','Positive support association'),(1,'noise_already_associated','#999999','Already in pre-injection noise'),(2,'new_association','#009e73','New association after injection')]:
 vals=[support[case][key] for case in cases];bars=ax.bar(x+(j-1)*width,vals,width,color=color,label=label);ax.bar_label(bars,padding=3)
ax.set_xticks(x,['Stable + probes','2x level + probes','Line + probes']);ax.set_ylim(0,78);ax.set_ylabel('Support associations / 65 paired probes');ax.set_title('What changed when the known signal was added?\nOracle noise pairs; support overlap is not physical identification');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.savefig(out/'probe_noise_association.png',dpi=180);plt.close(fig)
print('Saved warning-impact and paired-noise figures.')
