"""Cross-check saved seed CSV, component JSON and region metrics without waveform access."""
import csv,hashlib,json,math,os,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True

def guard(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
  p=os.path.abspath(os.fsdecode(args[0])).lower().replace('\\','/')
  if '/data/' in p or 'o3_holdout' in p:raise PermissionError('All waveform/O3 reads blocked in saved-record verifier')
 if event in ('socket.connect','socket.getaddrinfo'):raise PermissionError('Network blocked')
sys.addaudithook(guard)
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)

def main():
 out=Path(sys.argv[1]).resolve();assert out.parent==ROOT/'results/gw170814_cluster_maps'
 m=load(out/'manifest.json');assert m['status']=='COMPLETE'
 assert sha(out/'protocol.md')==m['protocol_sha256'] and sha(out/'driver.py')==m['driver_sha256']
 assert all(sha(ROOT/f)==h for f,h in m['source_sha256'].items())
 assert all(sha(ROOT/f)==h for f,h in m['parent_sha256'].items())
 records=load(out/'component_records.json');scales=load(out/'scale_summary.json');regions=load(out/'region_metrics.json')
 pixels=list(csv.DictReader((out/'seed_pixels.csv').open(encoding='utf-8')))
 bycomp=defaultdict(list);byscale=defaultdict(list)
 for r in pixels:
  k=(r['site'],r['scale']);key=(*k,int(r['component_id']));bycomp[key].append(r);byscale[k].append(r)
  z=float(r['Z']);assert z>4
  close(z,float(r['power'])/float(r['background_psd'])/float(r['background_time']))
 for c in records:
  rr=bycomp[(c['site'],c['scale'],c['component_id'])]
  assert len(rr)==c['n_bins'] and c['qualified']==(len(rr)>=2)
  assert all((r['qualified']=='True')==c['qualified'] and int(r['component_size'])==len(rr) for r in rr)
  close(sum(float(r['Z'])-1 for r in rr),c['raw_score'])
  coords={(int(r['freq_index']),int(r['time_index'])) for r in rr};reached={next(iter(coords))};todo=list(reached)
  while todo:
   a,b=todo.pop()
   for da in (-1,0,1):
    for db in (-1,0,1):
     q=(a+da,b+db)
     if q in coords and q not in reached:reached.add(q);todo.append(q)
  assert reached==coords
  for j in range(4):
   inside=[r for r in rr if 35<=float(r['frequency_hz'])<=250 and j<=float(r['time_s'])<j+1]
   assert len(inside)==c['regional_seed_counts'][j]
   close(sum(float(r['Z'])-1 for r in inside),c['regional_seed_contributions'][j])
 for s in scales:
  key=(s['site'],s['scale']);rr=byscale[key];cc=[c for c in records if (c['site'],c['scale'])==key]
  assert len(rr)==s['seed_bins'] and len(cc)==s['all_components']
  coordinate_labels={(int(r['freq_index']),int(r['time_index'])):int(r['component_id']) for r in rr}
  for (a,b),label in coordinate_labels.items():
   for da in (-1,0,1):
    for db in (-1,0,1):assert coordinate_labels.get((a+da,b+db),label)==label
  assert sum(c['qualified'] for c in cc)==s['raw_clusters_qualified']
  winner=max((c for c in cc if c['qualified'] and c['parent_bbox_gate']),key=lambda c:c['raw_score'])
  assert winner==s['near_reference_winner']
  cal=load(ROOT/'results/gw170814_current_control/20261126T174011Z_e21955d7'/s['site']/'calibration.json')
  st=cal['scale_stats'][s['scale']]
  close((winner['raw_score']-st['median_null_score'])/st['spread_scale'],s['near_standardized_score'])
  close(st['median_null_score']+cal['global_threshold']*st['spread_scale'],s['raw_acceptance_boundary'])
 for q in regions:
  key=(q['site'],q['scale']);j=q['region_start_s'];rr=[r for r in byscale[key] if 35<=float(r['frequency_hz'])<=250 and j<=float(r['time_s'])<j+1]
  assert len(rr)==q['seed_cells'];close(sum(float(r['Z'])-1 for r in rr),q['seed_excess_sum'])
  contributions=defaultdict(float)
  for r in rr:contributions[int(r['component_id'])]+=float(r['Z'])-1
  assert len(contributions)==q['all_components_touching']
  largest=max(contributions.values(),default=0.);close(largest,q['largest_component_regional_contribution'])
  if q['seed_excess_sum']:close(largest/q['seed_excess_sum'],q['largest_share'])
 result={'status':'VERIFIED','seed_pixel_rows_checked':len(pixels),'connected_components_checked':len(records),'region_metrics_checked':len(regions),'scales_checked':len(scales),'scope':'Separate saved CSV/JSON cross-check; no detector execution or waveform reads; no new calibrated statistic'}
 (out/'verify_records.py').write_bytes(Path(__file__).read_bytes());(out/'record_crosscheck.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
 print(json.dumps(result))
if __name__=='__main__':main()
