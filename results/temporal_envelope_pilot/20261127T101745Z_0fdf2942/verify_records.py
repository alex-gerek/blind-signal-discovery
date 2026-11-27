"""Independent saved-record and synthetic-map verification for temporal pilot."""
import hashlib,json,math,os,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  if p.startswith(str(ROOT/'data').lower().replace('\\','/')+'/') or 'o3_holdout' in p:raise PermissionError('All waveform/O3 reads blocked in verifier')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
from scipy import signal,ndimage
from scipy.signal.windows import tukey
from bsd.noise import generate_noise
from bsd.generator import generate_clean_signal,_scale_to_energy_snr
ARMS=('baseline','slow');REGIMES=('stationary','smooth_drift','persistent_step');FAMILIES=('linear_chirp','band_limited_noise_burst')
SPECS=[('short_31ms',32,28,5,31,257),('medium_125ms',128,96,11,15,33),('long_500ms',512,384,21,7,9)]

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ah(x):return hashlib.sha256(np.asarray(x,dtype='<f8',order='C').tobytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
def rng(*k):return np.random.default_rng(np.random.SeedSequence(k))
def noise(g):return generate_noise('pink_colored',4096,1024.,g,sigma=1.,params={'alpha':1.,'f_knee':16.})[0]
def overlap(c,g):return g is not None and c['t_start']<=g[1] and c['t_end']>=g[0] and c['f_min']<=g[3] and c['f_max']>=g[2]

def profile(base,seed,kind,i,ri):
 t=np.arange(4096)/1024.;g=rng(seed,4,kind,i,ri)
 if ri==0:return base.copy(),np.ones(4096)
 if ri==1:env=1+.5*np.sin(2*np.pi*.25*t+float(g.uniform(0,2*np.pi)))
 else:
  when=float(g.uniform(1,3));factor=2. if g.random()<.5 else .5;env=np.where(t<when,1.,factor)
 x=base*env;x=x*(1/np.sqrt(np.mean(x*x)));return x,env

def clean(seed,i,fi):
 g=rng(seed,3,i,fi)
 if fi==0:
  params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
  value,a,b,lo,hi,_=generate_clean_signal('linear_chirp',4096,1024.,g,params);return value,[a,b,lo,hi]
 t=np.arange(4096)/1024.;active=(t>=1.375)&(t<=1.625);full=signal.sosfiltfilt(signal.butter(4,[35,350],fs=1024.,btype='bandpass',output='sos'),g.normal(size=4096));value=np.zeros(4096);value[active]=full[active]*tukey(int(active.sum()),alpha=.15)
 return value,[1.375,1.625,35.,350.]

def map_pair(x,spec):
 name,n,ov,fk,tk,sk=spec
 f,t,h=signal.stft(x-np.median(x),fs=1024.,window='hann',nperseg=n,noverlap=ov,boundary=None,padded=False)
 p=abs(h[1:-1])**2;s=np.maximum(ndimage.median_filter(np.median(p,axis=1)/np.log(2),size=fk,mode='nearest'),1e-18);w=p/s[:,None];raw=np.median(w,axis=0)/np.log(2)
 result={}
 for arm,kernel in [('baseline',tk),('slow',sk)]:
  env=np.maximum(ndimage.median_filter(raw,size=kernel,mode='nearest'),1e-18);z=w/env[None,:];labels,ncomp=ndimage.label(z>4,structure=np.ones((3,3),int))
  scores=[float(np.sum(z[labels==j]-1)) for j in range(1,ncomp+1) if np.count_nonzero(labels==j)>=2]
  result[arm]={'T':env,'raw_max':max(scores,default=0.)}
 return t,result

def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/temporal_envelope_pilot'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE'
 for n,k in [('driver.py','driver_sha256'),('temporal_envelope_prototype.py','helper_sha256'),('protocol.md','protocol_sha256')]:assert sha(out/n)==m[k]
 assert all(sha(ROOT/p)==h for p,h in m['source_sha256'].items()) and all(sha(ROOT/p)==h for p,h in m['parent_sha256'].items())
 counts={'calibrations':0,'arm_window_decisions':0,'accepted_candidates':0,'input_regeneration_records':0,'negative_tracking_records':0,'paired_response_records':0,'manual_map_arm_checks':0,'manual_tracking_proxy_checks':0,'manual_paired_ratio_checks':0}
 alltracking=[];allpaired=[]
 for case in ('H1','L1','synthetic_9101','synthetic_9102','synthetic_9103'):
  p=out/case;features=load(p/'calibration_raw_features.json');rows=load(p/'evaluation.json');cc=load(p/'candidates.json');summary=load(p/'summary.json');cals={};nulls={}
  for arm in ARMS:
   cal=load(p/f'{arm}_calibration.json');null=load(p/f'{arm}_null_scores.json');cals[arm]=cal;nulls[arm]=null;n=cal['n_norm_samples'];assert n==round(.3*len(features)) and len(null)==len(features)-n
   for name,s in cal['scale_stats'].items():
    a=[r[arm][name] for r in features[:n]];med=float(np.median(a));q=float(np.quantile(a,.95));close(s['median_null_score'],med);close(s['q95_null_score'],q);close(s['spread_scale'],max(q-med,1.))
   scores=[max((r[arm][name]-s['median_null_score'])/s['spread_scale'] for name,s in cal['scale_stats'].items()) for r in features[n:]]
   np.testing.assert_allclose(sorted(scores),null,rtol=1e-12,atol=1e-12)
   tail=max(i for i in range(len(null)+1) if i/(len(null)+1)<=.01);assert cal['global_threshold']==null[len(null)-tail]
   counts['calibrations']+=1
  candidates=defaultdict(list)
  for c in cc:candidates[(c['sample_id'],c['arm'])].append(c)
  for r in rows:
   for arm in ARMS:
    d=r[arm];cal=cals[arm];null=nulls[arm];accepted=candidates[(r['sample_id'],arm)];gate=r['evaluation_gate']
    score=max((d['per_scale_raw_max'][name]-s['median_null_score'])/s['spread_scale'] for name,s in cal['scale_stats'].items());close(score,d['score'])
    pv=(1+sum(v>=d['score'] for v in null))/(len(null)+1)
    assert pv==d['empirical_p'] and d['threshold']==cal['global_threshold']
    assert (d['score']>d['threshold'])==(pv<=.01)==d['detected_window']==bool(accepted)
    assert d['associated_raw']==any(overlap(c,gate) for c in accepted)
    if d['near_score'] is not None:
     assert d['near_p']==(1+sum(v>=d['near_score'] for v in null))/(len(null)+1)
     assert d['associated_raw']==(d['near_score']>d['threshold'])
    if case in ('H1','L1'):
     assert d['post_veto_window']==any(not c['vetoed'] for c in accepted)
     assert d['associated_post_veto']==any(not c['vetoed'] and overlap(c,gate) for c in accepted)
    for c in accepted:
     s=cal['scale_stats'][c['scale_name']];close((c['raw_anomaly_score']-s['median_null_score'])/s['spread_scale'],c['standardized_score'])
     assert c['standardized_score']>d['threshold'] and c['empirical_p_value']==(1+sum(v>=c['standardized_score'] for v in null))/(len(null)+1)
     counts['accepted_candidates']+=1
    counts['arm_window_decisions']+=1
  for regime,d in summary['by_regime'].items():
   selected=[r for r in rows if r.get('regime','real')==regime]
   for arm in ARMS:
    for role,s in d['arms'][arm].items():
     rr=[r for r in selected if r['role']==role];assert len(rr)==s['windows'];assert s['window_nominations']==sum(r[arm]['detected_window'] for r in rr);assert s['associated_nominations']==sum(r[arm]['associated_raw'] for r in rr)
   for role,s in d['discordant'].items():
    rr=[r for r in selected if r['role']==role];key='associated_raw' if role.startswith('SIGNAL') or role=='ON_SOURCE' else 'detected_window'
    assert s=={'slow_only':sum(r['slow'][key] and not r['baseline'][key] for r in rr),'baseline_only':sum(r['baseline'][key] and not r['slow'][key] for r in rr)}
  if case in ('H1','L1'):
   par=load(ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'/case/'summary.json');event=next(r for r in rows if r['role']=='ON_SOURCE')
   for key in ('score','threshold','empirical_p'):close(event['baseline'][key],par['event'][key])
   assert summary['arms']['baseline']['OFF_SOURCE']['window_nominations']==par['offsource_raw'] and summary['arms']['baseline']['OFF_SOURCE']['post_veto_window_nominations']==par['offsource_post_veto']
   for r in rows:
    gps=1186741540+4*r['window_index'];assert r['gps_start']==gps
    expected='ON_SOURCE' if gps<=1186741861.5<gps+4 else ('OFF_SOURCE' if gps+5<=1186741861.5-16 or gps-1>=1186741861.5+16 else 'GUARD');assert r['role']==expected
   continue
  seed=int(case.split('_')[1]);tracking=load(p/'tracking_proxy.json');paired=load(p/'paired_normalizer_response.json');alltracking+=tracking;allpaired+=paired
  counts['negative_tracking_records']+=len(tracking);counts['paired_response_records']+=len(paired)
  tracklookup={(r['sample_id'],r['arm'],r['scale']):r for r in tracking};pairlookup={(r['sample_id'],r['arm'],r['scale']):r for r in paired}
  for i,r in enumerate(features):assert ah(noise(rng(seed,0,i)))==r['input_sha256'];counts['input_regeneration_records']+=1
  for r in rows:
   ri=REGIMES.index(r['regime'])
   if r['role']=='NOISE':
    i=r['noise_index'];base=noise(rng(seed,1,i));x,env=profile(base,seed,0,i,ri);assert ah(base)==r['base_noise_sha256'] and ah(x)==r['input_sha256'] and ah(env)==r['profile_params']['envelope_sha256'];gate=None
   else:
    i=r['injection_id'];fi=FAMILIES.index(r['family']);base=noise(rng(seed,2,i,fi));partner,env=profile(base,seed,1+fi,i,ri);cleanwave,gate=clean(seed,i,fi);added=_scale_to_energy_snr(cleanwave,float(np.std(partner)),r['snr']);x=partner+added
    assert gate==r['evaluation_gate']
    assert (ah(base),ah(partner),ah(cleanwave),ah(added),ah(x),ah(env))==(r['base_noise_sha256'],r['noise_sha256'],r['clean_unit_sha256'],r['signal_sha256'],r['input_sha256'],r['profile_params']['envelope_sha256'])
    close(np.linalg.norm(added)/np.std(partner),r['snr'])
   counts['input_regeneration_records']+=1
   if seed==9101 and i==0:
    for spec in SPECS:
     name=spec[0];t,manual=map_pair(x,spec)
     for arm in ARMS:
      close(manual[arm]['raw_max'],r[arm]['per_scale_raw_max'][name]);counts['manual_map_arm_checks']+=1
      if gate is None:
       truth=np.interp(t,np.arange(4096)/1024.,env**2);truth/=np.median(truth);est=manual[arm]['T']/np.median(manual[arm]['T']);proxy=float(np.sqrt(np.mean((np.log(est)-np.log(truth))**2)))
       close(proxy,tracklookup[(r['sample_id'],arm,name)]['log_rmse_proxy']);counts['manual_tracking_proxy_checks']+=1
     if gate is not None:
      _,noise_maps=map_pair(partner,spec);frames=(t>=gate[0])&(t<=gate[1]);assert frames.any()
      for arm in ARMS:close(float(np.median((manual[arm]['T']/noise_maps[arm]['T'])[frames])),pairlookup[(r['sample_id'],arm,name)]['temporal_ratio']);counts['manual_paired_ratio_checks']+=1
 for a in load(out/'tracking_summary.json'):
  vals=[r['log_rmse_proxy'] for r in alltracking if all(r[k]==a[k] for k in ('regime','arm','scale'))];assert len(vals)==a['negative_records']==600;close(float(np.median(vals)),a['median_log_rmse_proxy']);assert min(vals)==a['min'] and max(vals)==a['max']
 for a in load(out/'paired_response_summary.json'):
  vals=[r['temporal_ratio'] for r in allpaired if all(r[k]==a[k] for k in ('regime','family','snr','arm','scale'))];assert len(vals)==a['paired_records']==30;close(float(np.median(vals)),a['median_temporal_ratio']);assert min(vals)==a['min'] and max(vals)==a['max']
  assert a['T_increased']==sum(v>1 for v in vals) and a['T_decreased']==sum(v<1 for v in vals) and a['T_unchanged']==sum(v==1 for v in vals)
 result={'status':'VERIFIED',**counts,'scope':'Saved calibration/decision/summary checks; synthetic input regeneration; prespecified seed9101 first noise/injection manual maps and tracking/paired ratios; no real waveform/O3 reads'}
 (out/'verify_records.py').write_bytes(Path(__file__).read_bytes());(out/'independent_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
