"""Small exact enumeration cross-check for the independent KS lattice oracle."""
from pathlib import Path
from itertools import combinations
import json,sys
from verify_regime_warning_monitor import exact_ks_tail
out=Path(sys.argv[1]);records=[]
for n,m in [(1,1),(1,2),(2,2),(2,3),(3,3)]:
 statistics=[]
 for positions in combinations(range(n+m),n):
  a=b=0;maximum=0;positions=set(positions)
  for k in range(n+m):
   if k in positions:a+=1
   else:b+=1
   maximum=max(maximum,abs(m*a-n*b))
  statistics.append(maximum)
 for boundary in range(n*m+2):
  expected=sum(d>=boundary for d in statistics)/len(statistics)
  actual=exact_ks_tail(boundary,n,m)
  assert actual==expected,(n,m,boundary,actual,expected)
  records.append({'n':n,'m':m,'D_numerator':boundary,'p':actual,'enumerated_paths':len(statistics)})
result={'status':'PASS','checks':len(records),'scope':'Small untied IID null-path enumeration validates numerical KS oracle only, not correlated-window statistical assumptions','records':records}
(out/'ks_oracle_enumeration_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
print('Exact KS oracle matches',len(records),'small enumerated-null boundaries.')
