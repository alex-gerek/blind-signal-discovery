"""Supplemental paired preservation check; no new waveform measurements or inference."""
from pathlib import Path
import json,sys
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
assert load(out/'independent_verification.json')['status']=='PASS'
records=[]
for seed in [19101,19102,19103]:
 a=out/f'{seed}_stable';b=out/f'{seed}_calibration_intrusions'
 af=load(a/'input_features.json');bf=load(b/'input_features.json')
 assert [i for i,(x,y) in enumerate(zip(af,bf)) if x['input_sha256']!=y['input_sha256']]==[500,540,580,620]
 ac=load(a/'calibration_replacement.json');bc=load(b/'calibration_replacement.json')
 assert ac['scale_stats']==bc['scale_stats']
 assert load(a/'calibration_initial.json')==load(b/'calibration_initial.json')
 am=load(a/'measurements.json');bm=load(b/'measurements.json')
 assert all(am[k]['score']==bm[k]['score'] for k in am if am[k]['index']>=684)
 at=load(a/'timeline.json');bt=load(b/'timeline.json')
 assert all(x['available']==y['available'] for x,y in zip(at,bt))
 late=lambda t,p:[r for r in t if r['policy']==p and r['signal_kind']=='wanted' and r['index']>=684]
 for p in ['frozen_naive','age_only','age_and_state']:
  x=late(at,p);y=late(bt,p)
  assert len(x)==len(y)==5
  assert [am[r['measurement_id']]['score'] for r in x]==[bm[r['measurement_id']]['score'] for r in y]
  if p=='frozen_naive':assert [r['associated'] for r in x]==[r['associated'] for r in y]
  records.append({'seed':seed,'policy':p,'later_probes':5,'stable_associated':sum(r['associated'] for r in x),'intrusion_associated':sum(r['associated'] for r in y),
    'stable_replacement_threshold':ac['global_threshold'],'intrusion_replacement_threshold':bc['global_threshold']})
result={'status':'PASS','scope':'Prespecified replacement-contamination challenge; source differs at fourBwindows, later scores/availability unchanged; descriptive paired diagnostic only','records':records}
(out/'paired_contamination_verification.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(result,indent=2))
