"""Supplemental descriptive comparison of planned gain-step and stable core outputs."""
from pathlib import Path
import json,sys
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
assert load(out/'independent_verification.json')['status']=='PASS'
records=[]
for seed in range(20101,20106):
 a=load(out/str(seed)/'stable/search_measurements.json');b=load(out/str(seed)/'level_step/search_measurements.json')
 assert len(a)==len(b)==516
 diffs=[abs(x['score']-y['score']) for x,y in zip(a,b)]
 records.append({'seed':seed,'paired_target_windows':516,'max_absolute_score_difference':max(diffs),'different_search_decisions':sum(x['detected']!=y['detected'] for x,y in zip(a,b)),
  'different_p_values':sum(x['p']!=y['p'] for x,y in zip(a,b))})
result={'status':'COMPLETE','scope':'Supplemental descriptive comparison; raw amplitude changes need not imply changed detector statistic. Not a universal gain invariance proof or a new primary endpoint.','records':records}
(out/'gain_invariance_diagnostic.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(result,indent=2))
