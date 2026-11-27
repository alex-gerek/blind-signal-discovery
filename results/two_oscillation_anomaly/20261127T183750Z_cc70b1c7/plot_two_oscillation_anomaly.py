"""Descriptive anomaly counts and representative frozen-detector maps."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
v=load(out/'independent_verification.json');assert v['status']=='PASS'
cases=['paired_noise','A_alone','B_far_alone','B_close_alone','A_plus_B_far','A_plus_B_close','A_plus_A','exact_cancel','near_cancel','equivalent_single']
labels=['Paired\nnoise','A','B far','B close','A+B far','A+B close','A+A','A-A','A-0.95A','Single 2A'];x=np.arange(len(cases));width=.36
fig,ax=plt.subplots(figsize=(13,4.8),layout='constrained')
for delta,color,label,values in [(-width/2,'#0072b2','Anomaly nomination',[v['cases'][c]['detected'] for c in cases]),(width/2,'#e69f00','Same-scale disjoint pair (diagnostic)',[max(v['cases'][c]['far_disjoint_pair'],v['cases'][c]['close_disjoint_pair']) for c in cases])]:
 bars=ax.bar(x+delta,values,width,color=color,label=label);ax.bar_label(bars,padding=3,fontsize=8)
ax.set_xticks(x,labels,fontsize=9);ax.set_ylim(0,125);ax.set_ylabel('Count / 100 paired units');ax.set_title('Mixture anomaly detection is not source reconstruction\nComponent SNR18; mixtures not renormalized; cancellation has no observable addition')
ax.legend(fontsize=9,ncol=2);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.savefig(out/'mixture_counts.png',dpi=180);plt.close(fig)
cases=['A_alone','A_plus_B_far','A_plus_B_close','exact_cancel'];titles=['One component A','A + B, gap96Hz','A + B, gap2Hz','A - A: noise-only record']
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,case,title in zip(axes.flat,cases,titles):
 a=load(out/'22101'/('map_'+case+'.json'));z=np.array(a['normalized_excess_map']);im=ax.pcolormesh(a['times'],a['freqs'],np.log10(np.maximum(z,1e-12)),vmin=0,vmax=2,cmap='viridis',shading='auto')
 for c in a['accepted_long_candidates']:ax.add_patch(Rectangle((c['t_start'],c['f_min']),c['t_end']-c['t_start'],c['f_max']-c['f_min'],fill=False,edgecolor='#ff5555',linewidth=1.1))
 ax.set_ylim(0,320);ax.set_xlim(0,4);ax.set_title(title);ax.set_xlabel('Time within 4s window(s)');ax.set_ylabel('Frequency(Hz)')
fig.colorbar(im,ax=axes.ravel().tolist(),label='log10 normalized excess Z (display clipped0..2)')
fig.suptitle('Representative unit: long500ms map, accepted candidate boxes\nCarrier neighborhood overlap alone does not establish two sources',fontsize=12)
fig.savefig(out/'mixture_maps.png',dpi=180);plt.close(fig)
print('Saved mixture counts and representative normalized maps.')
