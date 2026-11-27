"""Independent common-budget verification; no data files or detector imports."""
import hashlib,json,math,os,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  if p.startswith(str(ROOT/'data').lower().replace('\\','/')+'/') or 'o3_holdout' in p:raise PermissionError('All waveform/O3 reads blocked')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
import numpy as np
from scipy import signal,ndimage,stats
from scipy.signal.windows import tukey
from bsd.noise import generate_noise
from bsd.generator import generate_clean_signal,_scale_to_energy_snr
ARMS=('baseline','slow');REGIMES=('stationary','smooth_drift','within_window_step');FAMILIES=('linear_chirp','band_limited_noise_burst')
SPECS=[('short_31ms',32,28,5,31,257),('medium_125ms',128,96,11,15,33),('long_500ms',512,384,21,7,9)]

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ah(x):return hashlib.sha256(np.asarray(x,dtype='<f8',order='C').tobytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
def rng(*k):return np.random.default_rng(np.random.SeedSequence(k))
def intersects(c,g):return g is not None and c['t_start']<=g[1] and c['t_end']>=g[0] and c['f_min']<=g[3] and c['f_max']>=g[2]

def input_noise(seed,i,fi):
 ri=int(rng(seed,10,i,fi).integers(3));g=rng(seed,4,i,fi,ri);t=np.arange(4096)/1024.;base=generate_noise('pink_colored',4096,1024.,rng(seed,0,i,fi),sigma=1.,params={'alpha':1.,'f_knee':16.})[0]
 if ri==0:envelope=np.ones(4096);x=base.copy();scale=1.
 else:
  if ri==1:envelope=1+.5*np.sin(2*np.pi*.25*t+float(g.uniform(0,2*np.pi)))
  else:
   when=float(g.uniform(1,3));factor=2. if g.random()<.5 else .5;envelope=np.where(t<when,1.,factor)
  x=base*envelope;scale=float(1/np.sqrt(np.mean(x*x)));x=x*scale
 return x,base,envelope,REGIMES[ri],scale

def input_clean(seed,i,fi):
 g=rng(seed,3,i,fi)
 if fi==0:
  params={'t_start_s':float(g.uniform(1.1,1.4)),'chirp_duration_s':float(g.uniform(.2,.5)),'f_start_hz':float(g.uniform(60,140)),'f_end_hz':float(g.uniform(200,340)),'phase_rad':float(g.uniform(0,2*np.pi))}
  x,a,b,lo,hi,_=generate_clean_signal('linear_chirp',4096,1024.,g,params);return x,[a,b,lo,hi]
 t=np.arange(4096)/1024.;active=(t>=1.375)&(t<=1.625);full=signal.sosfiltfilt(signal.butter(4,[35,350],fs=1024.,btype='bandpass',output='sos'),g.normal(size=4096));x=np.zeros(4096);x[active]=full[active]*tukey(int(active.sum()),alpha=.15)
 return x,[1.375,1.625,35.,350.]

def manual_raw(x,spec):
 name,n,ov,fk,tk,sk=spec
 _,_,h=signal.stft(x-np.median(x),fs=1024.,window='hann',nperseg=n,noverlap=ov,boundary=None,padded=False)
 power=abs(h[1:-1])**2;s=np.maximum(ndimage.median_filter(np.median(power,axis=1)/np.log(2),size=fk,mode='nearest'),1e-18);w=power/s[:,None];raw=np.median(w,axis=0)/np.log(2);result={}
 for arm,k in [('baseline',tk),('slow',sk)]:
  env=np.maximum(ndimage.median_filter(raw,size=k,mode='nearest'),1e-18);z=w/env[None,:];labels,ncomp=ndimage.label(z>4,structure=np.ones((3,3),int));scores=[float(np.sum(z[labels==j]-1)) for j in range(1,ncomp+1) if np.count_nonzero(labels==j)>=2];result[arm]=max(scores,default=0.)
 return result

def check_interval(record):
 k,n=record['count'],record['n']
 if n==0:assert record['rate'] is None and record['lower95'] is None and record['upper95'] is None;return
 assert record['rate']==k/n and 0<=record['lower95']<=k/n<=record['upper95']<=1
 if k:close(float(stats.binom.sf(k-1,n,record['lower95'])),.025)
 else:assert record['lower95']==0
 if k<n:close(float(stats.binom.cdf(k,n,record['upper95'])),.025)
 else:assert record['upper95']==1

def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/common_budget_comparison'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE'
 for name,key in [('driver.py','driver_sha256'),('temporal_envelope_prototype.py','helper_sha256'),('protocol.md','protocol_sha256')]:assert sha(out/name)==m[key]
 assert all(sha(ROOT/p)==h for p,h in m['source_sha256'].items())
 freeze=load(out/'threshold_freeze.json');assert freeze['heldout_records_generated']==0 and freeze['thresholds_will_not_be_retuned'];assert sha(out/'threshold_freeze.json')==m['threshold_freeze_sha256'];assert all(sha(out/p)==h for p,h in freeze['files_sha256'].items())
 budget=load(out/'budget_rule.json');n=1800
 mt=max(k for k in range(1,n+1) if stats.beta.ppf(.975,k,n+1-k)<=.01)
 assert mt==budget['tail_rank_m']==10 and budget['order_k_1based']==1791 and budget['operating_alpha']==10/1801
 assert stats.binom.cdf(mt-1,n,.01)<=.025 and stats.binom.cdf(mt,n,.01)>.025
 close(float(stats.beta.ppf(.975,mt,n+1-mt)),budget['tolerance_tail_mass_upper975'])
 features=load(out/'calibration_raw_features.json');noise=load(out/'heldout_noise.json');signals=load(out/'heldout_signals.json');candidates=load(out/'candidates.json');inference=load(out/'inference.json')
 assert len(features)==2000 and len(noise)==3600 and len(signals)==1200
 assert all(r['split']=='A' for r in features[:200]) and all(r['split']=='B' for r in features[200:])
 cals={};nulls={};manualchecks=0;inputchecks=0;decisionchecks=0;candidatechecks=0
 for arm in ARMS:
  cal=load(out/f'{arm}_calibration.json');null=load(out/f'{arm}_null_scores.json');cals[arm]=cal;nulls[arm]=null
  assert (cal['n_norm_samples'],cal['n_threshold_samples'])==(200,1800) and cal['target_fap']==10/1801
  for name,s in cal['scale_stats'].items():
   a=[r[arm][name] for r in features[:200]];median=float(np.median(a));q=float(np.quantile(a,.95));close(s['median_null_score'],median);close(s['q95_null_score'],q);close(s['spread_scale'],max(q-median,1.))
  vals=[max((r[arm][name]-s['median_null_score'])/s['spread_scale'] for name,s in cal['scale_stats'].items()) for r in features[200:]];np.testing.assert_allclose(sorted(vals),null,rtol=1e-12,atol=1e-12)
  assert cal['global_threshold']==null[1790]
 indexed=defaultdict(list)
 for c in candidates:indexed[(c['sample_id'],c['arm'])].append(c)
 for r in noise+signals:
  for arm in ARMS:
   d=r[arm];cal=cals[arm];null=nulls[arm];cc=indexed[(r['sample_id'],arm)];score=max((d['per_scale_raw_max'][name]-s['median_null_score'])/s['spread_scale'] for name,s in cal['scale_stats'].items());close(score,d['score'])
   p=(1+sum(x>=d['score'] for x in null))/1801;assert p==d['empirical_p'] and d['threshold']==cal['global_threshold'] and d['operating_alpha']==10/1801
   assert (d['score']>d['threshold'])==(p<=10/1801)==d['detected_window']==bool(cc);assert d['associated']==any(intersects(c,r['evaluation_gate']) for c in cc)
   for c in cc:
    s=cal['scale_stats'][c['scale_name']];close((c['raw_anomaly_score']-s['median_null_score'])/s['spread_scale'],c['standardized_score']);assert c['standardized_score']>d['threshold'];assert c['empirical_p_value']==(1+sum(x>=c['standardized_score'] for x in null))/1801;candidatechecks+=1
   decisionchecks+=1
 for r in features+noise+signals:
  gen=r['generation'];seed=gen['seed'];i=gen['index'];fi=gen['family_key'];x,base,env,regime,scale=input_noise(seed,i,fi)
  assert (ah(base),ah(x),ah(env),regime)==(gen['base_noise_sha256'],gen['noise_sha256'],gen['envelope_sha256'],gen['regime']);close(scale,gen['global_scale'])
  if r.get('role')=='SIGNAL':
   assert seed==15004 and fi==FAMILIES.index(r['family']) and i==r['unit_index']
   clean,gate=input_clean(seed,i,fi);added=_scale_to_energy_snr(clean,float(np.std(x)),r['snr']);x=x+added
   assert gate==r['evaluation_gate'] and ah(clean)==r['clean_sha256'] and ah(added)==r['signal_sha256'];close(np.linalg.norm(added)/np.std(input_noise(seed,i,fi)[0]),r['snr'])
  else:assert seed==({'A':15001,'B':15002}[r['split']] if 'split' in r else 15003)
  assert ah(x)==r['input_sha256'];inputchecks+=1
  if i<3:
   for spec in SPECS:
    expected=manual_raw(x,spec)
    for arm in ARMS:
     actual=r[arm][spec[0]] if 'split' in r else r[arm]['per_scale_raw_max'][spec[0]];close(expected[arm],actual);manualchecks+=1
 assert inputchecks==6800 and manualchecks==126
 sets={'A':{r['input_sha256'] for r in features[:200]},'B':{r['input_sha256'] for r in features[200:]},'noise':{r['input_sha256'] for r in noise},'signal_noise_partners':{r['generation']['noise_sha256'] for r in signals}}
 assert [len(sets[k]) for k in ('A','B','noise','signal_noise_partners')]==[200,1800,3600,600]
 names=list(sets);assert all(not sets[a]&sets[b] for j,a in enumerate(names) for b in names[j+1:])
 units=defaultdict(list)
 for r in signals:units[(r['family'],r['unit_index'])].append(r)
 assert len(units)==600
 for rr in units.values():assert sorted(r['snr'] for r in rr)==[8,18] and rr[0]['generation']==rr[1]['generation'] and rr[0]['clean_sha256']==rr[1]['clean_sha256']
 for arm in ARMS:
  rec=inference['noise_window_results'][arm];check_interval(rec);assert rec['count']==sum(r[arm]['detected_window'] for r in noise) and rec['n']==3600
  close(rec['false_alarm_windows_per_hour_equivalent'],rec['rate']*900)
  for regime,rec in inference['noise_by_regime'][arm].items():
   check_interval(rec);rr=[r for r in noise if r['generation']['regime']==regime];assert rec['n']==len(rr) and rec['count']==sum(r[arm]['detected_window'] for r in rr)
 assert inference['budget_supported_by_evaluation_upper_bounds']==all(r['upper95']<=.01 for r in inference['noise_window_results'].values())
 assert inference['budget_alpha']==.01 and inference['operating_alpha']==10/1801
 assert inference['noise_discordance']=={'slow_only':sum(r['slow']['detected_window'] and not r['baseline']['detected_window'] for r in noise),'baseline_only':sum(r['baseline']['detected_window'] and not r['slow']['detected_window'] for r in noise)}
 assert sorted((g['family'],g['snr']) for g in inference['signal_groups'])==sorted((f,a) for f in FAMILIES for a in (8,18))
 assert sum(g['primary'] for g in inference['signal_groups'])==1
 for group in inference['signal_groups']:
  rr=[r for r in signals if r['family']==group['family'] and r['snr']==group['snr']];assert len(rr)==group['n']==300
  b=sum(r['slow']['associated'] and not r['baseline']['associated'] for r in rr);c=sum(r['baseline']['associated'] and not r['slow']['associated'] for r in rr);assert (b,c)==(group['slow_only'],group['baseline_only'])
  for arm in ARMS:check_interval(group[arm]);assert group[arm]['count']==sum(r[arm]['associated'] for r in rr)
  for regime,d in group['by_regime'].items():
   subset=[r for r in rr if r['generation']['regime']==regime]
   for arm in ARMS:check_interval(d[arm]);assert d[arm]['n']==len(subset) and d[arm]['count']==sum(r[arm]['associated'] for r in subset)
  if group['primary']:
   assert group['family']=='band_limited_noise_burst' and group['snr']==8
   d=b+c;exact=1. if d==0 else min(1.,2*sum(math.comb(d,j) for j in range(min(b,c)+1))/(2**d));close(exact,group['exact_mcnemar_p_two_sided'])
   flag=inference['budget_supported_by_evaluation_upper_bounds'] and group['slow']['count']>group['baseline']['count'] and exact<.05
   assert flag==inference['primary_improvement_evidence_for_this_generator']
 assert inference['automatic_adoption'] is False
 result={'status':'VERIFIED','calibrations_checked':2,'input_records_regenerated':inputchecks,'manual_map_arm_checks':manualchecks,'arm_window_decisions_checked':decisionchecks,'accepted_candidates_checked':candidatechecks,'threshold_freeze_hashes_unchanged':True,'noise_hash_split_intersections':0,'paired_signal_units_checked':len(units),'budget_rule_verified_by_beta_and_binomial':True,'intervals_verified_by_binomial_tail_inversion':True,'primary_test_verified_by_exact_combinatorial_sum':True,'scope':'Independent saved-record/math/input checks; synthetic maps only; no detector/real waveform/O3/network use; no retuning'}
 (out/'verify_records.py').write_bytes(Path(__file__).read_bytes());(out/'independent_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
