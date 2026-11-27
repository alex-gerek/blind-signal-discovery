"""Independent calibration, mixture-identity and morphology accounting oracle."""
from pathlib import Path
from itertools import combinations
import hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
CASES=['A_alone','B_far_alone','B_close_alone','A_plus_B_far','A_plus_B_close','A_plus_A','exact_cancel','near_cancel','equivalent_single']
def load(p):return json.loads(p.read_text(encoding='utf-8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def touch(c,lo,hi):return max(c['t_start'],1.25)<=min(c['t_end'],2.25) and max(c['f_min'],lo)<=min(c['f_max'],hi)
def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/two_oscillation_anomaly'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE'
 report={'status':'RUNNING','calibrations':0,'input_records':0,'measurements':0,'candidate_records':0,'regenerated_inputs':0,'heldout_noise_alarms':0,'heldout_noise_windows':0,
  'cases':{case:{'units':0,'detected':0,'extra_nominations_vs_paired_noise':0,'A_band':0,'B_far_band':0,'B_close_band':0,'far_disjoint_pair':0,'close_disjoint_pair':0,'one_covers_far':0,'one_covers_close':0,'resultant_snr_min':1e100,'resultant_snr_max':0.} for case in ['paired_noise']+CASES},
  'scope':'Paired synthetic anomaly/morphology test, no waveform/source reconstruction'}
 for seed in m['seeds']:
  d=out/str(seed);features=load(d/'input_features.json');rows=load(d/'measurements.json');cal=load(d/'calibration.json');freeze=load(d/'freeze.json')
  assert len(features)==684 and len(rows)==400 and len({r['sample_id'] for r in features})==684
  assert hashlib.sha256((d/'calibration.json').read_bytes()).hexdigest()==freeze['calibration_sha256'] and freeze['target_decisions_measured']==0
  source=[r for r in features if r['part']=='calibration'];assert len(source)==284
  assert cal['n_norm_samples']==85 and cal['n_threshold_samples']==199
  stats={}
  for name,s in cal['scale_stats'].items():
   a=[r['raw_maxima'][name] for r in source[:85]];med=float(np.median(a));q=float(np.quantile(a,.95));spread=max(q-med,1.)
   assert (s['median_null_score'],s['q95_null_score'],s['spread_scale'])==(med,q,spread);stats[name]=(med,spread)
  def standard(r):return {k:(r['raw_maxima'][k]-med)/spread for k,(med,spread) in stats.items()}
  null=sorted(max(standard(r).values()) for r in source[85:]);assert null==cal['sorted_global_null_scores'] and cal['global_threshold']==null[-1]
  by_id={r['sample_id']:r for r in features};units={}
  for row in rows:
   raw=by_id[row['sample_id']];z=standard(raw);score=max(z.values());p=(1+sum(x>=score for x in null))/200.
   assert row['score']==score and row['threshold']==cal['global_threshold'] and row['p']==p and row['available']
   assert row['per_scale_std_max']==z and row['detected']==bool(score>row['threshold'] and p<=.005)==bool(row['candidates'])
   for c in row['candidates']:
    med,spread=stats[c['scale_name']];cp=(1+sum(x>=c['standardized_score'] for x in null))/200.
    assert c['standardized_score']==(c['raw_anomaly_score']-med)/spread and cp==c['empirical_p_value'] and cp<=.005 and c['standardized_score']>row['threshold']
    assert 0<=c['t_start']<=c['t_end']<=4 and 0<=c['f_min']<=c['f_max']<=512 and c['n_bins']>=2
    report['candidate_records']+=1
   if row['part']=='heldout_noise':report['heldout_noise_alarms']+=row['detected'];report['heldout_noise_windows']+=1
   else:
    unit=row['unit'];units.setdefault(unit,{})[row['case']]=row
    np.testing.assert_allclose(row['signal_energy'],row['component1_energy']+row['component2_energy']+row['cross_term'],rtol=1e-12,atol=1e-10)
    np.testing.assert_allclose(row['resultant_snr'],np.sqrt(row['signal_energy'])/row['noise_sigma'],rtol=1e-12,atol=1e-12)
    assert row['observable_addition']==(row['signal_energy']>0)
    bands={'A':(row['f_a_hz']-4,row['f_a_hz']+4),'B_far':(row['f_far_hz']-4,row['f_far_hz']+4),'B_close':(row['f_close_hz']-4,row['f_close_hz']+4)}
    flags={name:[touch(c,*band) for c in row['candidates']] for name,band in bands.items()}
    assert row['carrier_band_associations']=={name:any(v) for name,v in flags.items()}
    def two(other):
     for i,j in combinations(range(len(row['candidates'])),2):
      a=row['candidates'][i];b=row['candidates'][j]
      if a['scale_name']==b['scale_name'] and (a['f_max']<b['f_min'] or b['f_max']<a['f_min']):
       if (flags['A'][i] and flags[other][j]) or(flags['A'][j] and flags[other][i]):return True
     return False
    assert row['same_scale_disjoint_far_pair']==two('B_far') and row['same_scale_disjoint_close_pair']==two('B_close')
    assert row['one_candidate_covers_A_and_far']==any(a and b for a,b in zip(flags['A'],flags['B_far']))
    assert row['one_candidate_covers_A_and_close']==any(a and b for a,b in zip(flags['A'],flags['B_close']))
   report['measurements']+=1
  assert len(units)==20
  for index,group in units.items():
   assert set(group)==set(['paired_noise']+CASES)
   for left,right in [('paired_noise','exact_cancel'),('A_plus_A','equivalent_single')]:
    assert group[left]['input_sha256']==group[right]['input_sha256']
    for key in ['score','p','detected','candidates','per_scale_std_max','carrier_band_associations','same_scale_disjoint_far_pair','same_scale_disjoint_close_pair']:assert group[left][key]==group[right][key]
   assert group['exact_cancel']['signal_energy']==0 and group['exact_cancel']['resultant_snr']==0
   np.testing.assert_allclose(group['near_cancel']['resultant_snr'],.9,atol=1e-12)
   np.testing.assert_allclose(group['A_plus_A']['resultant_snr'],36.,atol=1e-12)
   for case,r in group.items():
    g=report['cases'][case];g['units']+=1;g['detected']+=r['detected'];g['extra_nominations_vs_paired_noise']+=r['detected'] and not group['paired_noise']['detected']
    for name,key in [('A','A_band'),('B_far','B_far_band'),('B_close','B_close_band')]:g[key]+=r['carrier_band_associations'][name]
    g['far_disjoint_pair']+=r['same_scale_disjoint_far_pair'];g['close_disjoint_pair']+=r['same_scale_disjoint_close_pair'];g['one_covers_far']+=r['one_candidate_covers_A_and_far'];g['one_covers_close']+=r['one_candidate_covers_A_and_close']
    g['resultant_snr_min']=min(g['resultant_snr_min'],r['resultant_snr']);g['resultant_snr_max']=max(g['resultant_snr_max'],r['resultant_snr'])
  calhash={r['input_sha256'] for r in source};held={r['input_sha256'] for r in features if r['part']=='heldout_noise'};paired={r['noise_sha256'] for r in features if r['part']=='unit'}
  assert not(calhash&held or calhash&paired or held&paired)
  report['input_records']+=len(features);report['calibrations']+=1
 spec=importlib.util.spec_from_file_location('mixture_control',ROOT/'experiments/run_two_oscillation_anomaly.py');driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
 assert Path(driver.__file__).read_bytes()==(out/'run_two_oscillation_anomaly.py').read_bytes()
 for seed in m['seeds']:
  f=load(out/str(seed)/'input_features.json');lookup={r['sample_id']:r for r in f}
  generated=driver.backgrounds(seed,1,284,'cal')+driver.backgrounds(seed,2,200,'heldout')
  for o in generated:assert driver.key(o)==lookup[o.sample_id]['input_sha256'];report['regenerated_inputs']+=1
  for index in range(20):
   for o,meta in driver.unit(seed,index):
    r=lookup[o.sample_id];assert all(r[k]==v for k,v in meta.items()) and driver.key(o)==r['input_sha256'];report['regenerated_inputs']+=1
 assert report['input_records']==report['regenerated_inputs']==3420 and report['measurements']==2000
 for name,h in {**m['source_sha256'],**m['dependency_sha256']}.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
 report.update(status='PASS',source_hashes_unchanged=True,direct_uncached_parity_checks=m['raw_parity_checks'])
 save(out/'independent_verification.json',report)
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
