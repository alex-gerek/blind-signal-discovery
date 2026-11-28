"""Matched grid/phase/amplitude interventions; unchanged research estimators."""
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
from itertools import product,combinations
import os,sys,json,gzip,hashlib,importlib.util,subprocess,time,traceback,uuid
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.dont_write_bytecode=True
import numpy as np
from scipy import signal
from bsd.generator import Observation
from bsd.detector_v2 import MultiResolutionDetector
spec=importlib.util.spec_from_file_location('blind',ROOT/'experiments/run_blind_two_tone_exploratory.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
PROTOCOL=ROOT/'docs/matched_gap2_protocol_2026-11-28.md'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):
 text=json.dumps(x,default=old.norm,allow_nan=False,indent=2)+'\n'
 if p.suffix=='.gz':
  with gzip.open(p,'wt',encoding='utf-8') as f:f.write(text)
 else:p.write_text(text,encoding='utf-8')
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def generate(base,kind,offset=0.,phase=0.,ratio=1.):
 n,_=old.make(41201,base,'noise');r=np.random.default_rng(np.random.SeedSequence([41202,base]));carrier=float(2*r.integers(40,90));global_phase=float(r.uniform(0,2*np.pi));t=np.arange(1024)/1024.;envelope=signal.windows.tukey(1024,.15)
 a=envelope*np.sin(2*np.pi*(carrier+offset)*t+global_phase);b=ratio*envelope*np.sin(2*np.pi*(carrier+offset+2)*t+global_phase+phase) if kind=='two' else np.zeros_like(a);energy=float((a+b)@(a+b));assert energy>1e-10;scale=18/np.sqrt(energy);a*=scale;b*=scale;addition=np.zeros(4096);addition[1280:2304]=a+b
 if kind=='noise':addition[:]=0
 meta={'part':'matched','kind':kind,'base':base,'offset':offset,'relative_phase':phase if kind=='two' else None,'ratio':ratio if kind=='two' else None,'carrier':carrier,'global_phase':global_phase,'f1':carrier+offset if kind!='noise' else None,'f2':carrier+offset+2 if kind=='two' else None,'gap':2 if kind=='two' else None,'snr':18 if kind!='noise' else 0,'observed_snr':float(np.linalg.norm(addition)),'noise_hash':old.key(n),'hash':old.key(n+addition),'e1':float(a@a) if kind!='noise' else 0.,'e2':float(b@b) if kind!='noise' else 0.,'cross':float(2*a@b) if kind!='noise' else 0.}
 if kind!='noise':np.testing.assert_allclose(meta['e1']+meta['e2']+meta['cross'],324,rtol=1e-12,atol=1e-10)
 return n+addition,meta

def main():
 global OUT
 assert old.git('rev-parse','HEAD')=='fc99a63915cc5137976ea493df1fb61fe6d803d6';assert not old.git('diff','--','src','app','pyproject.toml')
 assert digest(ROOT/'experiments/run_blind_two_tone_exploratory.py')==load(ROOT/'results/blind_two_tone_exploratory/20261127T192047Z_f279641e/manifest.json')['driver_sha256']
 OUT=ROOT/'results/matched_gap2'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);OUT.mkdir(parents=True);print('RUN_DIRECTORY',OUT,flush=True)
 sources={p.relative_to(ROOT).as_posix():digest(p) for p in (ROOT/'src/bsd').glob('*.py')};history={p.relative_to(ROOT).as_posix():digest(p) for folder in ['blind_two_tone_exploratory','blind_two_tone_gap2_review'] for p in (ROOT/'results'/folder).rglob('*') if p.is_file()}
 m={'status':'RUNNING','source_commit':old.git('rev-parse','HEAD'),'dirty':old.git('status','--short'),'source_sha256':sources,'historical_sha256':history,'python':sys.version,'numpy':np.__version__,'scipy':old.scipy.__version__,'driver_sha256':digest(Path(__file__)),'dependency_sha256':digest(Path(old.__file__))};save(OUT/'manifest.json',m)
 (OUT/'driver.py').write_bytes(Path(__file__).read_bytes());(OUT/'blind_dependency.py').write_bytes(Path(old.__file__).read_bytes());(OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes())
 save(OUT/'fixtures.json',old.fixtures())
 # Verify the full predeclared matched generator before fitting or target evaluation.
 for base in range(12):
  metas=[]
  for offset,phase,ratio in product([0.,.5],[0.,np.pi/2,np.pi],[1.,.5,.25]):
   x,meta=generate(base,'two',offset,phase,ratio);metas.append(meta);assert x.shape==(4096,) and np.isfinite(x).all();np.testing.assert_allclose(meta['observed_snr'],18.,rtol=1e-12)
  assert len({a['noise_hash'] for a in metas})==len({a['carrier'] for a in metas})==len({a['global_phase'] for a in metas})==1
 print('MATCHING_FIXTURES PASS',flush=True)
 rawcal=[];cals={};cal_times={}
 for w in [0,1]:
  class Logged(MultiResolutionDetector):
   def _raw_analyze_all_scales(self,o):
    r=super()._raw_analyze_all_scales(o);rawcal.append({'w':w,'id':o.sample_id,'hash':old.key(o.signal),'raw':{k:{'max':v[0],'candidates':[asdict(c) for c in v[1]]} for k,v in r.items()}});return r
  observations=[Observation(f'cal_{i}',old.make(40101,i,'noise')[0],1024.,4.) for i in range(284)];t=time.perf_counter();cals[w]=Logged(old.scales(w)).calibrate(observations,target_fap=.005,norm_fraction=.3,calibration_seed=40101);cal_times[w]=time.perf_counter()-t;save(OUT/f'calibration_W{w}.json',asdict(cals[w]));print('ANOMALY_CALIBRATED',w,cals[w].global_threshold,flush=True)
 save(OUT/'calibration_raw.json.gz',rawcal)
 nullrows=[];nulls={};cuts={};parity=[]
 for w in [0,1]:
  nulls[w]={d:{'noise':[],'one':[]} for d in [0,1]}
  for kind,seed in [('noise',40201),('one',40301)]:
   for i in range(199):
    x,meta=old.make(seed,i,kind);res,_,_=old.extract(x,w,cals[w]);nullrows.append({'w':w,'meta':meta,'result':res})
    for d in [0,1]:ps=res['methods'][d]['proposals'];nulls[w][d][kind].append(ps[0]['score'] if ps else None)
  for d in [0,1]:
   cut={k:old.cutoff(v) for k,v in nulls[w][d].items()};values=[v for v in cut.values() if v is not None];cuts[f'W{w}D{d}']={'by_group':cut,'cutoff':max(values) if values else None}
  print('PAIR_CALIBRATED',w,flush=True)
 save(OUT/'pair_calibration.json.gz',nullrows);save(OUT/'frozen_thresholds.json',{'cuts':cuts,'nulls':nulls,'heldout_decisions':0,'calibration_hashes':{w:digest(OUT/f'calibration_W{w}.json') for w in [0,1]},'null_means_negative_infinity':True});freeze=digest(OUT/'frozen_thresholds.json')
 jobs=[]
 for kind,seed in [('noise',41101),('one',41102)]:
  for i in range(160):jobs.append(('independent',kind,seed,i,None))
 for base in range(12):
  jobs.append(('matched','noise',None,base,(0.,0.,1.)))
  for offset in [0.,.5]:jobs.append(('matched','one',None,base,(offset,0.,1.)))
  for offset,phase,ratio in product([0.,.5],[0.,np.pi/2,np.pi],[1.,.5,.25]):jobs.append(('matched','two',None,base,(offset,phase,ratio)))
 assert len(jobs)==572
 rows=[];traces=OUT/'traces';traces.mkdir();start=time.perf_counter()
 for j,(part,kind,seed,i,params) in enumerate(jobs):
  if part=='independent':x,meta=old.make(seed,i,kind);meta['part']=part;ident=f'{seed}_{i:03d}'
  else:x,meta=generate(i,kind,*params);ident=f'b{i:02d}_{kind}_{params[0]:g}_{params[1]:.8g}_{params[2]:g}'
  for w in [0,1]:
   result,maps,labels=old.extract(x,w,cals[w]);decisions={}
   for d in [0,1]:
    arm=f'W{w}D{d}';ps=result['methods'][d]['proposals'];p=ps[0] if ps else None;claim=old.claims(p,cuts[arm]['cutoff']);frequencies=p['frequencies'] if claim else None;strict=loose=ambiguous=False;errors=None
    if kind=='two':strict,_,errors=old.match(frequencies,[meta['f1'],meta['f2']],.5);loose,ambiguous,_=old.match(frequencies,[meta['f1'],meta['f2']],1.)
    longest=[c for c in result['accepted'] if c['scale']==old.scales(w)[-1].name]
    stage='claim' if claim else ('no_anomaly' if not result['anomaly'] else ('no_long_candidate' if not longest else ('no_pair_proposal' if not ps else 'pair_cutoff')))
    timing={'onset':p['start']-1.25,'end':p['end']-2.25,'midpoint':p['peak']-1.75} if claim and part=='matched' and kind!='noise' else None
    decisions[arm]={'claim':claim,'frequencies':frequencies,'strict_correct':strict,'loose_correct':loose,'ambiguous':ambiguous,'errors':errors,'false_pair':claim and kind!='two','wrong_pair':claim and kind=='two' and not strict,'wrong_pair_loose':claim and kind=='two' and not loose,'stage':stage,'timing':timing}
   trace=None
   if part=='matched':
    t=time.perf_counter();f,ts,stft=signal.stft(x-np.median(x),fs=1024.,window='hann',nperseg=old.scales(w)[-1].detector_config.nperseg,noverlap=old.scales(w)[-1].detector_config.noverlap,boundary=None,padded=False);np.testing.assert_allclose(abs(stft[1:-1])**2,maps.power_spectrogram,rtol=1e-12,atol=1e-15);path=traces/f'{ident}_W{w}.npz';np.savez_compressed(path,waveform=x,stft=stft[1:-1],power=maps.power_spectrogram,psd=maps.background_psd,time_background=maps.background_time_envelope,z=maps.normalized_excess_map,freqs=maps.freqs,times=maps.times,mask=maps.normalized_excess_map>4,labels=labels);trace=path.relative_to(OUT).as_posix()
   if j in [0,160,320,323]:
    direct=MultiResolutionDetector(old.scales(w)).analyze(Observation('parity',x,1024.,4.),cals[w]);assert direct.detected==result['anomaly'] and direct.max_global_score==result['score'] and len(direct.candidates)==len(result['accepted']);parity.append([j,w])
   rows.append({'id':ident,'w':w,'meta':meta,'result':result,'decisions':decisions,'trace':trace})
  if j%100==0:print('EVALUATED',j+1,len(jobs),flush=True)
 save(OUT/'evaluation.json.gz',rows);save(OUT/'parity.json',{'status':'PASS','checks':parity})
 groups={'cal':{r['hash'] for r in rawcal},'pair':{r['meta']['hash'] for r in nullrows},'evaluation':{r['meta']['hash'] for r in rows}}
 assert all(not(groups[a]&groups[b]) for a,b in combinations(groups,2));assert freeze==digest(OUT/'frozen_thresholds.json')
 assert all(digest(ROOT/p)==h for p,h in {**sources,**history}.items())
 m.update(status='COMPLETE',evaluation_rows=len(rows),evaluation_inputs=len(jobs),calibration_seconds=cal_times,evaluation_elapsed_seconds=time.perf_counter()-start,freeze_sha256=freeze,source_and_history_unchanged=True);save(OUT/'manifest.json',m);print('COMPLETE',OUT,flush=True)
if __name__=='__main__':
 try:main()
 except Exception:
  if 'OUT' in globals():save(OUT/'failure.json',{'status':'FAILED','traceback':traceback.format_exc()})
  raise
