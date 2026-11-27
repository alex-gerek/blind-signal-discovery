"""Independent saved-record calibration checks, input regeneration and 48-map audit."""
import csv,hashlib,json,math,os,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  if p.startswith(str(ROOT/'data').lower().replace('\\','/')+'/') or 'o3_holdout' in p:raise PermissionError('All original/derived/O3 data blocked in verifier')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
from scipy import signal,ndimage
from scipy.signal.windows import tukey
from bsd.noise import generate_noise
from bsd.generator import generate_clean_signal,_scale_to_energy_snr

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ah(x):return hashlib.sha256(np.asarray(x,dtype='<f8',order='C').tobytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
def rng(*k):return np.random.default_rng(np.random.SeedSequence(k))
def noise(model,g):return generate_noise(model,4096,1024.,g,sigma=1.,params={'alpha':1.,'f_knee':16.})[0]
def overlap(c,g):return g is not None and c['t_start']<=g[1] and c['t_end']>=g[0] and c['f_min']<=g[3] and c['f_max']>=g[2]
def maps(x,n,overlap,fk,tk):
 f,t,h=signal.stft(x-np.median(x),fs=1024.,window='hann',nperseg=n,noverlap=overlap,boundary=None,padded=False)
 p=abs(h[1:-1])**2;nf,nt=p.shape
 fk=min(fk,nf if nf%2 else nf-1);tk=min(tk,nt if nt%2 else nt-1)
 s=np.maximum(ndimage.median_filter(np.median(p,axis=1)/np.log(2),size=max(3,fk),mode='nearest'),1e-18)
 w=p/s[:,None];env=np.maximum(ndimage.median_filter(np.median(w,axis=0)/np.log(2),size=max(3,tk),mode='nearest'),1e-18)
 return f[1:-1],t,p,s,env,w/env[None,:]

def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/gap_normalizer_controls'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE'
 for name,key in [('driver.py','driver_sha256'),('gap_link_prototype.py','helper_sha256'),('protocol.md','protocol_sha256')]:assert sha(out/name)==m[key]
 assert all(sha(ROOT/f)==h for f,h in m['source_sha256'].items()) and all(sha(ROOT/f)==h for f,h in m['parent_sha256'].items())
 countcal=0;countdecision=0;countcand=0;regen=0
 for case in ('H1','L1','synthetic_7301','synthetic_7302','synthetic_7303'):
  p=out/case;features=load(p/'calibration_raw_features.json');rows=load(p/'evaluation.json');candidates=load(p/'candidates.json');summary=load(p/'summary.json');cals={};nulls={}
  for r in features:
   for scale in r['baseline']:assert r['gap'][scale]+1e-12>=r['baseline'][scale]
  for arm in ('baseline','gap'):
   c=load(p/f'{arm}_calibration.json');null=load(p/f'{arm}_null_scores.json');cals[arm]=c;nulls[arm]=null
   n=c['n_norm_samples'];assert n==round(.3*len(features)) and len(null)==len(features)-n
   for scale,st in c['scale_stats'].items():
    a=[r[arm][scale] for r in features[:n]];med=float(np.median(a));q=float(np.quantile(a,.95));spread=max(q-med,1.)
    close(st['median_null_score'],med);close(st['q95_null_score'],q);close(st['spread_scale'],spread)
   vals=[max((r[arm][scale]-st['median_null_score'])/st['spread_scale'] for scale,st in c['scale_stats'].items()) for r in features[n:]]
   np.testing.assert_allclose(sorted(vals),null,rtol=1e-12,atol=1e-12)
   tail=max(i for i in range(len(null)+1) if i/(len(null)+1)<=.01);k=len(null)+1-tail
   assert c['global_threshold']==null[k-1]
   countcal+=1
  bykey=defaultdict(list)
  for c in candidates:bykey[(c['sample_id'],c['arm'])].append(c)
  for r in rows:
   gate=None
   if r['role']=='ON_SOURCE':gate=[1186741861.5-r['gps_start']-.5,1186741861.5-r['gps_start']+.5,35,250]
   elif r['role'].startswith('SIGNAL'):
    g=r['generation_params'];gate=[g['t_start_s'],g['t_end_s'],min(g['f_start_hz'],g['f_end_hz']),max(g['f_start_hz'],g['f_end_hz'])]
   for scale in r['baseline']['per_scale_raw_max']:assert r['gap']['per_scale_raw_max'][scale]+1e-12>=r['baseline']['per_scale_raw_max'][scale]
   for arm in ('baseline','gap'):
    d=r[arm];cal=cals[arm];null=nulls[arm];cc=bykey[(r['sample_id'],arm)]
    calculated=max((d['per_scale_raw_max'][scale]-st['median_null_score'])/st['spread_scale'] for scale,st in cal['scale_stats'].items());close(calculated,d['score'])
    pv=(1+sum(x>=d['score'] for x in null))/(len(null)+1)
    assert pv==d['empirical_p'];assert d['threshold']==cal['global_threshold']
    assert (d['score']>d['threshold'])==(pv<=.01)==d['detected_window']==bool(cc)
    assert d['associated_raw']==any(overlap(c,gate) for c in cc)
    if d['near_score'] is not None:
     assert d['near_p']==(1+sum(x>=d['near_score'] for x in null))/(len(null)+1)
     assert d['associated_raw']==(d['near_score']>d['threshold'])
    if case in ('H1','L1'):
     assert d['post_veto_window']==any(not c['vetoed'] for c in cc)
     assert d['associated_post_veto']==any(not c['vetoed'] and overlap(c,gate) for c in cc)
    for c in cc:
     st=cal['scale_stats'][c['scale_name']];close((c['raw_anomaly_score']-st['median_null_score'])/st['spread_scale'],c['standardized_score'])
     assert c['standardized_score']>d['threshold']
     assert c['empirical_p_value']==(1+sum(x>=c['standardized_score'] for x in null))/(len(null)+1)
     countcand+=1
    countdecision+=1
  for arm in ('baseline','gap'):
   for role,v in summary['arms'][arm].items():
    rr=[r for r in rows if r['role']==role];assert v['windows']==len(rr)
    assert v['window_nominations']==sum(r[arm]['detected_window'] for r in rr)
    assert v['associated_nominations']==sum(r[arm]['associated_raw'] for r in rr)
    if case in ('H1','L1'):assert v['post_veto_window_nominations']==sum(r[arm]['post_veto_window'] for r in rr)
  if case.startswith('synthetic'):
   seed=int(case.split('_')[1])
   for i,r in enumerate(features):assert ah(noise('pink_colored',rng(seed,0,i)))==r['input_sha256'];regen+=1
   for r in rows:
    if r['role']=='NOISE':
     i=int(r['sample_id'].rsplit('_',1)[1]);assert ah(noise('pink_colored',rng(seed,1,i)))==r['input_sha256'];regen+=1
    else:
     i=r['injection_id'];base=noise('pink_colored',rng(seed,2,i));g=rng(seed,3,i)
     params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
     clean,*_=generate_clean_signal('linear_chirp',4096,1024.,g,params);scaled=_scale_to_energy_snr(clean,float(np.std(base)),r['snr'])
     assert ah(base)==r['noise_sha256'] and ah(clean)==r['clean_unit_sha256'] and ah(scaled)==r['signal_sha256'] and ah(base+scaled)==r['input_sha256'];regen+=1
  else:
   par=load(ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'/case/'summary.json');event=next(r for r in rows if r['role']=='ON_SOURCE')
   for key in ('score','threshold','empirical_p'):close(event['baseline'][key],par['event'][key])
   assert summary['arms']['baseline']['OFF_SOURCE']['windows']==151
   assert summary['arms']['baseline']['OFF_SOURCE']['window_nominations']==par['offsource_raw']
   assert summary['arms']['baseline']['OFF_SOURCE']['post_veto_window_nominations']==par['offsource_post_veto']
 inputs=load(out/'normalizer_inputs.json');records=load(out/'normalizer_records.json');lookup={(r['input_id'],r['scale']):r for r in records};mapchecks=0
 for r in records:
  assert 0<=r['production_seeds']<=r['roi_cells'] and 0<=r['temporal_frozen_seeds']<=r['roi_cells'] and 0<=r['fully_frozen_seeds']<=r['roi_cells']
  close(r['temporal_ratio'],r['temporal_counterfactual_Z_ratio'])
  if r['snr']==0:
   assert r['temporal_ratio']==r['spectral_ratio']==r['full_counterfactual_Z_ratio']==1.
   assert r['production_seeds']==r['temporal_frozen_seeds']==r['fully_frozen_seeds']
   close(r['production_positive_excess'],r['temporal_frozen_positive_excess']);close(r['production_positive_excess'],r['fully_frozen_positive_excess'])
 specs=[('short_31ms',32,28,5,31),('medium_125ms',128,96,11,15),('long_500ms',512,384,21,7)]
 for r in inputs:
  model_index=('gaussian_white','pink_colored').index(r['model']);family_index=('linear_chirp','band_limited_noise_burst').index(r['family']);i=r['noise_index'];seed=r['seed'];base=noise(r['model'],rng(seed,model_index,0,i));g=rng(seed,model_index,1,i,family_index)
  if family_index==0:
   params={'t_start_s':1.25,'chirp_duration_s':.35,'f_start_hz':40.,'f_end_hz':240.,'phase_rad':float(g.uniform(0,2*np.pi))}
   clean,*_=generate_clean_signal('linear_chirp',4096,1024.,g,params)
  else:
   t=np.arange(4096)/1024.;active=(t>=1.375)&(t<=1.625);full=signal.sosfiltfilt(signal.butter(4,[35,350],fs=1024.,btype='bandpass',output='sos'),g.normal(size=4096));clean=np.zeros(4096);clean[active]=full[active]*tukey(int(active.sum()),alpha=.15)
  scaled=np.zeros_like(base) if r['snr']==0 else _scale_to_energy_snr(clean,float(np.std(base)),r['snr']);mixture=base+scaled
  assert (ah(base),ah(clean),ah(scaled),ah(mixture))==(r['noise_sha256'],r['clean_sha256'],r['signal_sha256'],r['mixture_sha256']);regen+=1
  close(np.linalg.norm(scaled)/np.std(base),r['snr'])
  if seed==8401 and i==0:
   for name,n,ov,fk,tk in specs:
    f,t,p0,s0,t0,z0=maps(base,n,ov,fk,tk);_,_,p1,s1,t1,z1=maps(mixture,n,ov,fk,tk);a,b=r['support'];frames=(t>=a)&(t<=b);freq=(f>=35)&(f<=350);roi=freq[:,None]&frames[None,:];q=lookup[(r['id'],name)]
    close(float(np.median((t1/t0)[frames])),q['temporal_ratio']);close(float(np.median((s1/s0)[freq])),q['spectral_ratio'])
    for zz,seedkey,excesskey in [(z1,'production_seeds','production_positive_excess'),(p1/s1[:,None]/t0[None,:],'temporal_frozen_seeds','temporal_frozen_positive_excess'),(p1/s0[:,None]/t0[None,:],'fully_frozen_seeds','fully_frozen_positive_excess')]:
     assert int((zz[roi]>4).sum())==q[seedkey];close(float(np.maximum(zz[roi]-1,0).sum()),q[excesskey])
    mapchecks+=1
 aggs=load(out/'normalizer_aggregates.json')
 for a in aggs:
  rr=[r for r in records if all(r[k]==a[k] for k in ('family','snr','scale')) and r['noise_model']==a['noise_model']];values=[r['temporal_ratio'] for r in rr]
  assert len(rr)==a['paired_records']==30;close(float(np.median(values)),a['temporal_ratio_median']);assert min(values)==a['temporal_ratio_min'] and max(values)==a['temporal_ratio_max']
  assert a['T_increased']==sum(x>1 for x in values) and a['T_decreased']==sum(x<1 for x in values) and a['T_unchanged']==sum(x==1 for x in values)
 result={'status':'VERIFIED','calibrations_checked':countcal,'arm_window_decisions_checked':countdecision,'accepted_candidates_checked':countcand,'synthetic_input_arrays_regenerated':regen,'normalizer_scale_records_checked':len(records),'normalizer_aggregates_checked':len(aggs),'normalizer_maps_independently_recomputed':mapchecks,'normalizer_identity_scale_records':sum(r['snr']==0 for r in records),'scope':'Independent saved calibration/decision/summary checks and deterministic synthetic input regeneration; independent manual normalization equations on predefined 48-map subset; no real waveform/O3 reads'}
 (out/'verify_records.py').write_bytes(Path(__file__).read_bytes());(out/'independent_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
