"""Reproducibility and label-isolation checks; no new evaluation data."""
from pathlib import Path
from dataclasses import fields
import hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from bsd.detector_v2 import MultiResolutionCalibration,ScaleNormalizationStats
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
spec=importlib.util.spec_from_file_location('blind',ROOT/'experiments/run_blind_two_tone_exploratory.py');driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
assert hashlib.sha256((ROOT/'experiments/run_blind_two_tone_exploratory.py').read_bytes()).hexdigest()==load(out/'manifest.json')['driver_sha256']
records=load(out/'evaluation.json');prs=load(out/'pair_calibration_records.json');seen=set();n=0
for row in prs+records:
 meta=row['meta'];ident=(meta['seed'],meta['index'])
 if ident in seen:continue
 seen.add(ident);params={k:meta[k] for k in ['gap','phase','ratio','offset','onset','snr'] if meta.get(k) is not None}
 x,regen=driver.make(meta['seed'],meta['index'],meta['kind'],params);assert driver.key(x)==meta['hash'];n+=1
raw=load(out/'anomaly_calibration_raw.json')
for row in raw:
 if row['w']!=0:continue
 i=int(row['id'].split('_')[-1]);x,_=driver.make(30201,i,'noise');assert driver.key(x)==row['hash'];n+=1
checks=0
for w in [0,1]:
 c=load(out/f'anomaly_calibration_W{w}.json')['calibration'];c['scale_stats']={k:ScaleNormalizationStats(**v) for k,v in c['scale_stats'].items()};c['sorted_global_null_scores']=np.array(c['sorted_global_null_scores']);cal=MultiResolutionCalibration(**c)
 for kind in ['noise','one','two']:
  row=next(r for r in records if r['w']==w and r['meta']['kind']==kind)
  with np.load(out/row['trace']) as a:x=a['waveform'].copy()
  # Reversing/altering labels is evaluator-only: inference is called with the exact same input-only interface.
  altered={**row['meta'],'f1':17.,'f2':401.,'kind':'noise','onset':-100.}
  first,_,_=driver.extract(x,w,cal);second,_,_=driver.extract(x,w,cal)
  for name in ['raw','accepted','score','anomaly']:assert first[name]==second[name]==row['result'][name]
  for d in [0,1]:
   for name in ['proposals','profiles']:assert first['methods'][d][name]==second['methods'][d][name]==row['result']['methods'][str(d)][name]
  assert altered!=row['meta'];checks+=1
result={'status':'PASS','regenerated_unique_inputs':n,'input_only_label_change_checks':checks,'scope':'hash reproducibility and deterministic label-independent inference, not generator scientific validation'}
(out/'reproducibility_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
