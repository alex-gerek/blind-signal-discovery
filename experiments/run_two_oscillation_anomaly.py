"""Template-free anomaly control on paired two-oscillation mixtures; no inverse recovery."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,json,platform,sys,time,traceback,uuid
import numpy as np
from scipy.signal.windows import tukey
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'));sys.dont_write_bytecode=True
from run_block_tail_guard import CachedCore,replay,key,save
from bsd.generator import Observation,_scale_to_energy_snr
from bsd.noise import generate_pink_colored
from bsd.detector import TemplateFreeDetector
from bsd.detector_v2 import DEFAULT_V02_SCALES
SEEDS=list(range(22101,22106))
CASES=['A_alone','B_far_alone','B_close_alone','A_plus_B_far','A_plus_B_close','A_plus_A','exact_cancel','near_cancel','equivalent_single']
PROTOCOL=ROOT/'docs/two_oscillation_anomaly_protocol_2026-11-27.md'
def sha(x):return hashlib.sha256(x).hexdigest()
def noise(rng):return generate_pink_colored(4096,1024.,rng,sigma=1.,alpha=1.,f_knee=16.)[0]
def observation(ident,x):return Observation(ident,x,1024.,4.)
def backgrounds(seed,stream,n,prefix):
 r=np.random.default_rng(np.random.SeedSequence([seed,stream]))
 return [observation(f's{seed}_{prefix}_{i:03d}',noise(r)) for i in range(n)]
def unit(seed,index):
 r=np.random.default_rng(np.random.SeedSequence([seed,3,index]));n=noise(r);f=int(r.integers(96,177));phases=r.uniform(0,2*np.pi,3)
 t=np.arange(4096)/1024.;envelope=np.zeros(4096);envelope[1280:2304]=tukey(1024,.15)
 waves=[_scale_to_energy_snr(envelope*np.sin(2*np.pi*freq*t+phase),1.,18.) for freq,phase in zip([f,f+96,f+2],phases)]
 a,b,c=waves;zero=np.zeros_like(a)
 components={'A_alone':(a,zero),'B_far_alone':(b,zero),'B_close_alone':(c,zero),'A_plus_B_far':(a,b),'A_plus_B_close':(a,c),'A_plus_A':(a,a),'exact_cancel':(a,-a),'near_cancel':(a,-.95*a),'equivalent_single':(2*a,zero)}
 assert np.array_equal(a+a,2*a) and np.array_equal(a-a,zero)
 params={'unit':index,'f_a_hz':f,'f_far_hz':f+96,'f_close_hz':f+2,'phases_rad':phases.tolist(),'noise_sigma':float(np.std(n))}
 result=[]
 base=observation(f's{seed}_u{index:02d}_noise',n)
 result.append((base,{**params,'case':'paired_noise','input_sha256':key(base),'noise_sha256':key(base),'component1_energy':0.,'component2_energy':0.,'cross_term':0.,'signal_energy':0.,'resultant_snr':0.,'observable_addition':False}))
 for j,case in enumerate(CASES):
  x,y=components[case];addition=x+y;e1=float(np.dot(x,x));e2=float(np.dot(y,y));cross=float(2*np.dot(x,y));energy=float(np.dot(addition,addition))
  np.testing.assert_allclose(energy,e1+e2+cross,rtol=1e-12,atol=1e-10)
  o=observation(f's{seed}_u{index:02d}_v{j:02d}',n+addition)
  result.append((o,{**params,'case':case,'input_sha256':key(o),'noise_sha256':key(base),'component1_energy':e1,'component2_energy':e2,'cross_term':cross,'signal_energy':energy,'resultant_snr':float(np.sqrt(energy)/np.std(n)),'observable_addition':energy>0}))
 assert key(result[0][0])==key(result[1+CASES.index('exact_cancel')][0])
 assert key(result[1+CASES.index('A_plus_A')][0])==key(result[1+CASES.index('equivalent_single')][0])
 return result

def intersects(c,lo,hi):return c['t_end']>=1.25 and c['t_start']<=2.25 and c['f_max']>=lo and c['f_min']<=hi

def diagnostics(candidates,meta):
 bands={'A':(meta['f_a_hz']-4,meta['f_a_hz']+4),'B_far':(meta['f_far_hz']-4,meta['f_far_hz']+4),'B_close':(meta['f_close_hz']-4,meta['f_close_hz']+4)}
 masks={name:[intersects(c,*band) for c in candidates] for name,band in bands.items()}
 def pair(other):
  return any(i!=j and a['scale_name']==b['scale_name'] and (a['f_max']<b['f_min'] or b['f_max']<a['f_min']) and masks['A'][i] and masks[other][j] for i,a in enumerate(candidates) for j,b in enumerate(candidates))
 return {'carrier_band_associations':{name:any(flags) for name,flags in masks.items()},'same_scale_disjoint_far_pair':pair('B_far'),'same_scale_disjoint_close_pair':pair('B_close'),
  'one_candidate_covers_A_and_far':any(a and b for a,b in zip(masks['A'],masks['B_far'])),'one_candidate_covers_A_and_close':any(a and b for a,b in zip(masks['A'],masks['B_close']))}

def main():
 assert replay.git('rev-parse','HEAD')=='cc68b93fee699307e3314fce1bd4c9fe403c1e6e'
 assert not replay.git('diff','--','src','app','pyproject.toml')
 global OUT
 ident=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8];OUT=ROOT/'results/two_oscillation_anomaly'/ident;OUT.mkdir(parents=True,exist_ok=False)
 print('RUN_DIRECTORY',OUT,flush=True)
 deps=[Path(__file__),ROOT/'experiments/run_block_tail_guard.py',ROOT/'experiments/block_tail_guard.py',ROOT/'experiments/run_historical_adaptive_replay.py']
 for p in deps:(OUT/p.name).write_bytes(p.read_bytes())
 (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes())
 sources=replay.hashes(list((ROOT/'src/bsd').glob('*.py')))
 m={'run_id':ident,'status':'RUNNING','source_commit':replay.git('rev-parse','HEAD'),'dirty_paths':replay.git('status','--short'),'source_sha256':sources,'dependency_sha256':replay.hashes(deps),'seeds':SEEDS,
   'python':platform.python_version(),'numpy':np.__version__,'scipy':replay.scipy.__version__,'cases':CASES,'seed_results':[],'raw_parity_checks':0,'scope':'Anomaly and candidate-structure detection, no source reconstruction'}
 save(OUT/'manifest.json',m);start=time.monotonic()
 # Pure vector identities; no evaluation sample or detector score used.
 a=np.arange(17,dtype=float)-8
 assert np.array_equal(a-a,np.zeros_like(a)) and np.array_equal(a+a,2*a)
 save(OUT/'identity_fixtures.json',{'status':'PASS','exact_negation_cancellation':True,'coherent_equivalence':True})
 for seed in SEEDS:
  dest=OUT/str(seed);dest.mkdir(exist_ok=False);det=CachedCore();reference=backgrounds(seed,1,284,'cal');neg=backgrounds(seed,2,200,'heldout')
  cal=det.calibrate(reference,target_fap=.005,norm_fraction=.30,calibration_seed=seed)
  save(dest/'calibration.json',asdict(cal));frozen=sha((dest/'calibration.json').read_bytes());save(dest/'freeze.json',{'target_decisions_measured':0,'calibration_sha256':frozen,'reference_hashes':[key(o) for o in reference]})
  inputs=[];rows=[]
  for o in reference:
   raw=det._raw_analyze_all_scales(o);inputs.append({'sample_id':o.sample_id,'part':'calibration','input_sha256':key(o),'raw_maxima':{k:float(v[0]) for k,v in raw.items()}})
  for o in reference[:3]+neg[:3]:det.parity_check(o)
  def measure(o,meta):
   raw=det._raw_analyze_all_scales(o);inputs.append({'sample_id':o.sample_id,**meta,'raw_maxima':{k:float(v[0]) for k,v in raw.items()}})
   result=det.analyze(o,cal);candidates=[asdict(c) for c in result.candidates]
   row={'sample_id':o.sample_id,**meta,'score':result.max_global_score,'threshold':cal.global_threshold,'p':cal.empirical_p_value(result.max_global_score)[0],'available':True,'detected':result.detected,'candidates':candidates,'per_scale_std_max':result.per_scale_std_max}
   if meta['part']=='unit':row.update(diagnostics(candidates,meta))
   rows.append(row);return row
  for o in neg:measure(o,{'part':'heldout_noise','case':'heldout_noise','input_sha256':key(o)})
  for index in range(20):
   results={}
   for o,meta in unit(seed,index):
    if index==0 and meta['case']!='paired_noise':det.parity_check(o)
    row=measure(o,{'part':'unit',**meta});results[meta['case']]=row
    if seed==SEEDS[0] and index==0 and meta['case'] in ['A_alone','A_plus_B_far','A_plus_B_close','exact_cancel']:
     maps=TemplateFreeDetector(DEFAULT_V02_SCALES[-1].detector_config).analyze(o,return_maps=True).tf_maps
     save(dest/('map_'+meta['case']+'.json'),{'case':meta['case'],'times':maps.times,'freqs':maps.freqs,'normalized_excess_map':maps.normalized_excess_map,'accepted_long_candidates':[c for c in row['candidates'] if c['scale_name']=='long_500ms']})
   for left,right in [('paired_noise','exact_cancel'),('A_plus_A','equivalent_single')]:
    assert results[left]['input_sha256']==results[right]['input_sha256']
    for k in ['score','threshold','p','detected','candidates','per_scale_std_max']:assert results[left][k]==results[right][k]
  save(dest/'input_features.json',inputs);save(dest/'measurements.json',rows)
  cal_hashes={r['input_sha256'] for r in inputs if r['part']=='calibration'};held={r['input_sha256'] for r in inputs if r['part']=='heldout_noise'};paired={r['noise_sha256'] for r in inputs if r['part']=='unit'}
  assert not(cal_hashes&held or cal_hashes&paired or held&paired)
  assert sha((dest/'calibration.json').read_bytes())==frozen
  summary={'seed':seed,'heldout_noise_alarms':sum(r['detected'] for r in rows if r['part']=='heldout_noise'),'heldout_noise_windows':200,
   'cases':{case:{'units':20,'detected':sum(r['detected'] for r in rows if r['case']==case),'far_disjoint_pair':sum(r.get('same_scale_disjoint_far_pair',False) for r in rows if r['case']==case),'close_disjoint_pair':sum(r.get('same_scale_disjoint_close_pair',False) for r in rows if r['case']==case)} for case in ['paired_noise']+CASES}}
  save(dest/'summary.json',summary);m['seed_results'].append(summary);m['raw_parity_checks']+=det.parity;save(OUT/'manifest.json',m)
  print('SEED_COMPLETE',seed,summary,flush=True)
 assert all(sha((ROOT/f).read_bytes())==h for f,h in {**sources,**m['dependency_sha256']}.items())
 m.update(status='COMPLETE',source_hashes_unchanged=True,elapsed_s=time.monotonic()-start);save(OUT/'manifest.json',m);print('RUN_STATUS COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception:
  if 'OUT' in globals():save(OUT/'failure.json',{'status':'FAILED','traceback':traceback.format_exc()})
  raise
