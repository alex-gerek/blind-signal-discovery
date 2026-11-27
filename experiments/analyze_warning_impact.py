"""Independent paired-impact and future-exposure analysis of verified saved outputs."""
from pathlib import Path
import csv,json,sys
out=Path(sys.argv[1]);load=lambda p:json.loads(p.read_text(encoding='utf-8'))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def association(c,s):return max(c['t_start'],s[0])<=min(c['t_end'],s[1]) and max(c['f_min'],s[2])<=min(c['f_max'],s[3])
def future_union(onsets):return sorted({i for onset in onsets for i in range(onset+1,min(onset+13,800)) if i>=284})
def main():
 v=load(out/'independent_verification.json');assert v['status']=='PASS'
 fixtures=[([],[]),([799],[]),([300,300],list(range(301,313))),([300,305],list(range(301,318))),([790,795],list(range(791,800)))]
 for a,b in fixtures:assert future_union(a)==b
 results={'status':'RUNNING','cadence_pairs':[],'gain_pairs':[],'injection_pairs':[],'future_noise':[],'future_exposure_fixture_checks':len(fixtures),
  'scope':'Descriptive paired synthetic impacts, no learned warning classifier or efficacy p-values'}
 repeated=list(range(288,800,8));sparse=repeated[2::5]
 for seed in range(21101,21106):
  sd=out/str(seed);features={};search={};warnings={}
  for case in ['stable','level_step','spectral_line','sparse_transients','repeated_transients','step_with_probes','line_with_probes']:
   features[case]=load(sd/case/'input_features.json');search[case]=load(sd/case/'search_measurements.json');warnings[case]=load(sd/case/'monitor_assessments.json')
  for index in sparse:
   a=features['sparse_transients'][index];b=features['repeated_transients'][index];sa=search['sparse_transients'][index-284];sb=search['repeated_transients'][index-284]
   assert a['input_sha256']==b['input_sha256'] and a['family']==b['family'] and a['support']==b['support']
   assert {k:v for k,v in sa.items() if k!='sample_id'}=={k:v for k,v in sb.items() if k!='sample_id'}
   wa=warnings['sparse_transients'][index-284]['regime_shifted'];wb=warnings['repeated_transients'][index-284]['regime_shifted']
   results['cadence_pairs'].append({'seed':seed,'index':index,'family':a['family'],'core_associated':sa['associated'],'sparse_warning':wa,'repeated_warning':wb,'warning_disagrees':wa!=wb,
     'sparse_gate_loss':bool(sa['associated'] and wa),'repeated_gate_loss':bool(sb['associated'] and wb),'input_and_core_identical':True})
  for baseline,changed in [('stable','level_step'),('sparse_transients','step_with_probes')]:
   a=search[baseline];b=search[changed]
   results['gain_pairs'].append({'seed':seed,'baseline':baseline,'changed':changed,'target_windows':516,
     'max_absolute_score_difference':max(abs(x['score']-y['score']) for x,y in zip(a,b)),
     'different_nomination_decisions':sum(x['detected']!=y['detected'] for x,y in zip(a,b)),
     'different_p_values':sum(x['p']!=y['p'] for x,y in zip(a,b)),
     'different_support_associations':sum(x['associated']!=y['associated'] for x,y in zip(a,b)),
     'baseline_warning_windows':sum(a['regime_shifted'] for a in warnings[baseline]),'changed_warning_windows':sum(a['regime_shifted'] for a in warnings[changed])})
  for positive,negative in [('sparse_transients','stable'),('step_with_probes','level_step'),('line_with_probes','spectral_line')]:
   for index in sparse:
    pf=features[positive][index];nf=features[negative][index];p=search[positive][index-284];n=search[negative][index-284];support=pf['support']
    assert pf['profiled_noise_sha256']==nf['input_sha256'] and pf['base_noise_sha256']==nf['base_noise_sha256']
    ph=any(association(c,support) for c in p['candidates']);nh=any(association(c,support) for c in n['candidates'])
    assert ph==p['associated']
    results['injection_pairs'].append({'seed':seed,'positive_case':positive,'noise_case':negative,'index':index,'family':pf['family'],
     'positive_association':ph,'noise_already_associated':nh,'new_association':bool(ph and not nh),'lost_association':bool(nh and not ph),
     'positive_nominated':p['detected'],'noise_nominated':n['detected'],'positive_only_nomination':bool(p['detected'] and not n['detected']),
     'positive_score':p['score'],'noise_score':n['score'],'score_difference':p['score']-n['score'],
     'warning_positive':warnings[positive][index-284]['regime_shifted'],'warning_noise':warnings[negative][index-284]['regime_shifted']})
  for case in ['stable','level_step','spectral_line']:
   w=warnings[case];flags=[a['regime_shifted'] for a in w];onsets=[i+284 for i,flag in enumerate(flags) if flag and (i==0 or not flags[i-1])]
   indices=future_union(onsets);bitmap=[False]*516
   for onset in onsets:
    for i in range(onset+1,min(onset+13,800)):bitmap[i-284]=True
   assert indices==[i+284 for i,x in enumerate(bitmap) if x]
   neg=search[case];alarms=sum(neg[i-284]['detected'] for i in indices)
   per_episode=[{'onset_index':onset,'future_indices':future_union([onset]),'alarms':sum(neg[i-284]['detected'] for i in future_union([onset]))} for onset in onsets]
   results['future_noise'].append({'seed':seed,'case':case,'warning_onsets':onsets,'unioned_future_indices':indices,'future_noise_windows':len(indices),
    'future_exposure_seconds':4.*len(indices),'future_noise_alarms':alarms,'conditional_rate':alarms/len(indices) if indices else None,'per_episode':per_episode})
 assert len(results['cadence_pairs'])==65 and len(results['injection_pairs'])==195 and len(results['gain_pairs'])==10
 results['summary']={'shared_probe_pairs':65,'shared_warning_disagreements':sum(r['warning_disagrees'] for r in results['cadence_pairs']),
   'shared_core_associations':sum(r['core_associated'] for r in results['cadence_pairs']),
   'shared_sparse_gate_losses':sum(r['sparse_gate_loss'] for r in results['cadence_pairs']),'shared_repeated_gate_losses':sum(r['repeated_gate_loss'] for r in results['cadence_pairs']),
   'max_gain_score_difference':max(r['max_absolute_score_difference'] for r in results['gain_pairs']),
   'gain_nomination_disagreements':sum(r['different_nomination_decisions'] for r in results['gain_pairs']),
   'injection_support':{case:{'pairs':len([r for r in results['injection_pairs'] if r['positive_case']==case]),**{key:sum(r[key] for r in results['injection_pairs'] if r['positive_case']==case) for key in ['positive_association','noise_already_associated','new_association','lost_association','positive_only_nomination']}} for case in ['sparse_transients','step_with_probes','line_with_probes']}}
 results['status']='PASS';save(out/'independent_impact_analysis.json',results)
 for name in ['cadence_pairs','gain_pairs','injection_pairs']:
  with (out/(name+'.csv')).open('w',encoding='utf-8',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(results[name][0]));writer.writeheader();writer.writerows(results[name])
 print(json.dumps(results['summary'],indent=2))
if __name__=='__main__':main()
