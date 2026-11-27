"""Descriptive timing and contamination figures; no new statistical inference."""
from pathlib import Path
import json,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
assert load(out/'independent_verification.json')['status']=='PASS'
cases=['stable','announced_line','unannounced_line','repeated_signals','calibration_intrusions'];policies=['frozen_naive','age_only','age_and_state'];short=['Frozen naive','Age only','Age + state']
fig,axes=plt.subplots(5,1,figsize=(11,10),sharex=True,layout='constrained')
for ax,case in zip(axes,cases):
 rows=load(out/f'19101_{case}/timeline.json')
 for j,p in enumerate(policies):
  seq=[r for r in rows if r['policy']==p]
  ax.broken_barh([(r['index'],1) for r in seq if r['available']],(j-.24,.48),facecolors='#b2dfdb')
  ax.broken_barh([(r['index'],1) for r in seq if not r['available']],(j-.24,.48),facecolors='#bdbdbd')
  for r in seq:
   if r['signal_kind'] is not None:
    marker,color=('^','#e69f00') if not r['available'] else (('o','#0072b2') if r['associated'] else ('x','#d55e00'))
    ax.scatter(r['index']+.5,j,s=25,marker=marker,c=color,zorder=4)
 ax.set_yticks(range(3),short,fontsize=8);ax.set_title(case.replace('_',' '),loc='left',fontsize=10)
 ax.axvline(684,color='#333333',linestyle='--',linewidth=.8)
 if 'line' in case:ax.axvline(360,color='#7570b3',linestyle=':',linewidth=1)
 ax.set_xlim(284,800);ax.set_ylim(-.5,2.5);ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
axes[-1].set_xlabel('Window index (4 s/window); replacement first used at 684')
legend=[Patch(color='#b2dfdb',label='Available under test rule'),Patch(color='#bdbdbd',label='Unavailable'),Line2D([0],[0],marker='o',color='none',markerfacecolor='#0072b2',markeredgecolor='#0072b2',label='Support association'),Line2D([0],[0],marker='x',color='#d55e00',linestyle='none',label='Available miss'),Line2D([0],[0],marker='^',color='none',markerfacecolor='#e69f00',label='Abstained probe')]
fig.legend(handles=legend,loc='outside lower center',ncol=3,fontsize=8)
fig.suptitle('Background lifetime pilot: seed19101, fixed unvalidated2400s age limit\nA finite replacement at684 does not itself demonstrate restored sensitivity.',fontsize=12)
fig.savefig(out/'eligibility_timeline.png',dpi=180);plt.close(fig)
records=load(out/'paired_contamination_verification.json')['records'];records=[r for r in records if r['policy']=='age_only']
fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained');x=list(range(3));width=.34
for side,(a,b),title,ylabel in [(0,('stable_replacement_threshold','intrusion_replacement_threshold'),'Replacement threshold','Standardized threshold'),(1,('stable_associated','intrusion_associated'),'Later wanted associations','Count / 5')]:
 ax=axes[side]
 for delta,key,color,label in [(-width/2,a,'#009e73','Clean replacement'),(width/2,b,'#d55e00','Four intrusions in B')]:
  vals=[r[key] for r in records];bars=ax.bar([i+delta for i in x],vals,width,color=color,label=label)
  ax.bar_label(bars,labels=[f'{v:.2f}' if side==0 else str(v) for v in vals],padding=3,fontsize=8)
 ax.set_xticks(x,[str(r['seed']) for r in records]);ax.set_title(title);ax.set_ylabel(ylabel);ax.set_ylim(0,(110 if side==0 else 6));ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
axes[0].legend(fontsize=8)
fig.suptitle('Paired replacement-source challenge: later scores and availability unchanged',fontsize=11)
fig.savefig(out/'replacement_contamination.png',dpi=180);plt.close(fig)
print('Saved descriptive timing and contamination figures.')
