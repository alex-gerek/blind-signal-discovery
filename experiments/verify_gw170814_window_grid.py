"""Independent arithmetic/decision checks of the saved window-grid experiment."""
from __future__ import annotations
import csv,hashlib,json,os,sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).replace('\\','/').lower()
  if '/data/' in p or 'o3_holdout' in p:raise PermissionError('All waveform data and O3 blocked for saved-record verification')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def near(c,t):return c['t_start']<=t+.5 and c['t_end']>=t-.5 and c['f_min']<=250 and c['f_max']>=35

def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/gw170814_window_grid'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE' and len(m['variants'])==8
 assert sha(out/'protocol.md')==m['protocol_sha256'] and sha(out/'driver.py')==m['driver_sha256']
 assert all(sha(ROOT/f)==h for f,h in m['source_sha256'].items())
 assert all(sha(ROOT/f)==h for f,h in m['parent_sha256'].items())
 rows=[];checked=0;clusters_checked=0
 for site in ('H1','L1'):
  ref=load(out/site/'phase_+0/calibration.json')
  for phase in (0,-2,-1,1):
   p=out/site/f'phase_{phase:+d}';s=load(p/'summary.json');cal=load(p/'calibration.json');null=load(p/'null_scores.json');chrono=load(p/'calibration_chronological_scores.json');cs=load(p/'candidates.json');raw=load(p/'event_clusters.json');event=s['event'];t=event['event_local_s']
   assert t==1.5-phase and cal['n_norm_samples']==212 and len(null)==495
   assert sorted(chrono)==null
   allowed=int(Fraction.from_float(.01)*496)
   assert cal['global_threshold']==null[495-allowed]
   timeline=list(csv.DictReader((p/'timeline.csv').open(encoding='utf-8')))
   expected_off=[]
   for i,r in enumerate(timeline):
    gps=1186741540+phase+4*i
    assert float(r['gps_start'])==gps
    role='ON_SOURCE' if gps<=1186741861.5<gps+4 else ('OFF_SOURCE' if gps+5<=1186741861.5-16 or gps-1>=1186741861.5+16 else 'GUARD')
    assert r['role']==role
    z=float(r['score']);pv=(1+sum(v>=z for v in null))/496
    assert float(r['empirical_p'])==pv
    assert float(r['threshold'])==cal['global_threshold']
    accepted=[c for c in cs if c['window_index']==i]
    assert (z>cal['global_threshold'])==(pv<=.01)==(r['detected_raw']=='True')==bool(accepted)
    assert (r['detected_post_veto']=='True')==any(not c['vetoed'] for c in accepted)
    if role=='OFF_SOURCE':expected_off.append(r)
    checked+=1
   assert len(timeline)==160 and s['offsource_windows']==len(expected_off)
   assert s['offsource_exposure_s']==4*len(expected_off)
   assert s['offsource_raw']==sum(r['detected_raw']=='True' for r in expected_off)
   assert s['offsource_post_veto']==sum(r['detected_post_veto']=='True' for r in expected_off)
   for c in cs:
    assert c['standardized_score']>cal['global_threshold']
    assert c['empirical_p_value']==(1+sum(v>=c['standardized_score'] for v in null))/496
   for c in raw:
    assert c['near_reference']==near(c,t)
    for key,model in (('standardized_matched',cal),('standardized_zero_calibration',ref)):
     st=model['scale_stats'][c['scale']]
     assert abs((c['anomaly_score']-st['median_null_score'])/st['spread_scale']-c[key])<1e-12
    clusters_checked+=1
   for field in ('score','threshold','empirical_p','detected_raw','detected_post_veto','gps_start'):
    tr=next(r for r in timeline if r['role']=='ON_SOURCE')
    val=tr[field]
    assert event[field]==(val=='True' if field.startswith('detected') else float(val))
   for scale,z in event['per_scale_std_max'].items():assert z==max(c['standardized_matched'] for c in raw if c['scale']==scale)
   assert event['score']==max(c['standardized_matched'] for c in raw)
   assert event['near_score']==max(c['standardized_matched'] for c in raw if near(c,t))
   assert event['zero_calibration']['score']==max(c['standardized_zero_calibration'] for c in raw)
   assert event['zero_calibration']['near_score']==max(c['standardized_zero_calibration'] for c in raw if near(c,t))
   ec=[c for c in cs if c['role']=='ON_SOURCE']
   assert event['associated_raw']==any(near(c,t) for c in ec)
   assert event['associated_post_veto']==any(near(c,t) and not c['vetoed'] for c in ec)
   for field,kind in (('near_score','near_empirical_p'),('score','empirical_p')):
    assert event[kind]==(1+sum(v>=event[field] for v in null))/496
   refnull=load(out/site/'phase_+0/null_scores.json')
   for field,kind in (('near_score','near_empirical_p'),('score','empirical_p')):
    rc=event['zero_calibration'];assert rc[kind]==(1+sum(v>=rc[field] for v in refnull))/496
   assert event['zero_calibration']['associated_raw']==(event['zero_calibration']['near_score']>ref['global_threshold'])
   if phase==0:
    par=load(ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'/site/'summary.json')
    assert all(abs(event[k]-par['event'][k])<1e-10 for k in ('score','threshold','empirical_p'))
    assert s['offsource_raw']==par['offsource_raw'] and s['offsource_post_veto']==par['offsource_post_veto']
   rows.append({'site':site,'phase_s':phase,'event_local_s':t,'near_score':event['near_score'],'threshold':event['threshold'],'near_p':event['near_empirical_p'],'associated_raw':event['associated_raw'],'associated_post_veto':event['associated_post_veto'],'zero_cal_near_score':event['zero_calibration']['near_score'],'zero_cal_threshold':ref['global_threshold'],'zero_cal_associated':event['zero_calibration']['associated_raw'],'offsource_raw':s['offsource_raw'],'offsource_post_veto':s['offsource_post_veto'],'offsource_windows':len(expected_off),'offsource_exposure_s':s['offsource_exposure_s']})
 with (out/'phase_comparison.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 fig,axes=plt.subplots(1,2,figsize=(10,4),sharex=True)
 for ax,site in zip(axes,('H1','L1')):
  rr=sorted((r for r in rows if r['site']==site),key=lambda r:r['phase_s']);x=[r['phase_s'] for r in rr]
  ax.plot(x,[r['near_score']/r['threshold'] for r in rr],'o-',label='Matched-grid calibration')
  ax.plot(x,[r['zero_cal_near_score']/r['zero_cal_threshold'] for r in rr],'s--',label='Original calibration (diagnostic)')
  ax.axhline(1,color='red',ls=':',label='Acceptance boundary; strict >')
  ax.set_title(site);ax.set_xticks(x);ax.set_xlabel('Grid shift relative to parent [s]');ax.grid(alpha=.2)
 axes[0].set_ylabel('Near-reference standardized score / threshold')
 axes[-1].legend(fontsize=8)
 fig.suptitle('GW170814 grid/context diagnostic: fixed methods, all four phases')
 fig.tight_layout();fig.savefig(out/'phase_comparison.png',dpi=150);plt.close(fig)
 record={'status':'VERIFIED','run_id':m['run_id'],'timeline_windows_checked':checked,'event_raw_clusters_checked':clusters_checked,'variants':len(rows),'scope':'Independent saved-record arithmetic, decisions, role/exposure and baseline checks; no waveform reload or detector execution','source_and_parent_hashes_match':True,'results':rows}
 (out/'independent_verification.json').write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
 (out/'verify_saved.py').write_bytes(Path(__file__).read_bytes())
 print(json.dumps(record,indent=2))
if __name__=='__main__':main()
