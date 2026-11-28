"""Forensic accounting of fixed W1D1 outcomes from saved maps only."""
from pathlib import Path
from itertools import combinations
from collections import Counter
from datetime import datetime,timezone
import csv,gzip,hashlib,json,sys,uuid
import numpy as np
from scipy.signal import find_peaks,peak_prominences
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/'results/matched_gap2/20261128T115035Z_bf957291'
def load(p):return json.load(gzip.open(p,'rt',encoding='utf-8')) if p.suffix=='.gz' else json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def truthpairs(indices,freq,truth):
 left=[int(i) for i in indices if abs(freq[i]-truth[0])<=.5];right=[int(i) for i in indices if abs(freq[i]-truth[1])<=.5]
 return [(i,j) for i in left for j in right if i!=j]
def main():
 out=ROOT/'results/matched_gap2_diagnosis'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True,exist_ok=False);(out/'driver.py').write_bytes(Path(__file__).read_bytes());print('OUTPUT',out,flush=True)
 records=load(PARENT/'evaluation.json.gz');rs=[r for r in records if r['w']==1 and r['meta']['kind']=='two'];assert len(rs)==216
 parent_hashes=load(PARENT/'artifact_sha256.json')['files'];assert all(sha(ROOT/p)==h for p,h in parent_hashes.items())
 cases=[];controls=[];categories=Counter();good=[];wrong=[];no=[];frame_evidence=Counter()
 for row in rs:
  meta=row['meta'];dec=row['decisions']['W1D1'];cs=[c for c in row['result']['accepted'] if c['scale']=='long_1000ms'];saved=row['result']['methods']['1']['profiles'];assert len(saved)==len(cs)
  with np.load(PARENT/row['trace']) as zfile:
   freq=zfile['freqs'];times=zfile['times'];z=zfile['z'];power=zfile['power'];labels=zfile['labels'];trace_stft=zfile['stft'];truth=[meta['f1'],meta['f2']];profiles=[]
   for c,s in zip(cs,saved):
    assert s['cluster_peak']==[c['peak_freq'],c['peak_time']];cols=np.array(s['columns'],int);y=z[:,cols].mean(axis=1);np.testing.assert_allclose(y,s['profile'],rtol=1e-12,atol=1e-12);pk=find_peaks(y)[0];inside=[int(i) for i in pk if c['f_min']<=freq[i]<=c['f_max']];strong=[i for i in inside if y[i]>=4];spaced=[(i,j) for i,j in combinations(strong,2) if j-i>=2];valid=[(i,j) for i,j in spaced if y[i+1:j].min()<=.8*min(y[i],y[j])];assert strong==s['peaks']
    frame_counts=Counter();eligible_frames=[]
    for j in cols:
     # Truth used after inference only, to annotate a fixed +/-0.5Hz frequency neighborhood.
     if not 1.25<=times[j]<=2.25:continue
     rp=find_peaks(power[:,j])[0];zp=find_peaks(z[:,j])[0];zi=[int(i) for i in zp if c['f_min']<=freq[i]<=c['f_max']];zt=[i for i in zi if z[i,j]>=4]
     frame_counts['raw_power_two_truth_peaks']+=bool(truthpairs(rp,freq,truth));frame_counts['normalized_two_truth_peaks']+=bool(truthpairs(zp,freq,truth));frame_counts['normalized_in_box_two_truth_peaks']+=bool(truthpairs(zi,freq,truth));frame_counts['normalized_in_box_two_above4']+=bool(truthpairs(zt,freq,truth))
     pairs=[(i,k) for i,k in truthpairs(zt,freq,truth) if abs(k-i)>=2 and z[min(i,k)+1:max(i,k),j].min()<=.8*min(z[i,j],z[k,j])]
     frame_counts['normalized_in_box_pair_passes_valley']+=bool(pairs)
     if pairs:eligible_frames.append(float(times[j]))
    profiles.append({'cluster_peak':s['cluster_peak'],'box':[c['f_min'],c['f_max']],'columns':cols.tolist(),'inside_local_peaks_hz':[float(freq[i]) for i in inside],'inside_peak_heights':[float(y[i]) for i in inside],'strong_peak_frequencies_hz':[float(freq[i]) for i in strong],'strong_peak_heights':[float(y[i]) for i in strong],'spaced_pair_count':len(spaced),'valid_pair_count':len(valid),'candidate_pair_valley_ratios':[float(y[i+1:j].min()/min(y[i],y[j])) for i,j in spaced],'mean_full_spectrum_has_truth_pair':bool(truthpairs(pk,freq,truth)),'mean_inside_box_has_truth_pair':bool(truthpairs(inside,freq,truth)),'frame_counts':dict(frame_counts),'eligible_frame_times':eligible_frames})
   item={'id':row['id'],'base':meta['base'],'truth_hz':truth,'offset':meta['offset'],'phase_rad':meta['relative_phase'],'amplitude_ratio':meta['ratio'],'anomaly':row['result']['anomaly'],'profiles':profiles,'estimated_hz':dec['frequencies'],'strict_correct':dec['strict_correct'],'loose_correct':dec['loose_correct'],'errors_hz':dec['errors']}
   if dec['claim']:
    proposal=row['result']['methods']['1']['proposals'][0];assert proposal['origin']=='within';pairs=[]
    for c,s in zip(cs,saved):
     y=np.array(s['profile']);idx=[int(np.flatnonzero(freq==f)[0]) for f in proposal['frequencies']]
     if not all(i in s['peaks'] for i in idx):continue
     i,j=idx;valley=float(y[i+1:j].min());heights=[float(y[i]),float(y[j])]
     if min(heights)!=proposal['score'] or valley>.8*min(heights):continue
     prom=peak_prominences(y,np.array(idx))[0];pairs.append({'peak_heights':heights,'peak_prominences':prom.tolist(),'valley_height':valley,'valley_ratio':valley/min(heights),'valley_depth_percent':100*(1-valley/min(heights)),'height_margin_over4':min(heights)-4,'valley_rule_margin':.8-valley/min(heights),'columns':s['columns'],'cluster_box':[c['f_min'],c['f_max']]})
    assert pairs;item['selected_pair_metrics']=pairs[0];item['category']='strict_success' if dec['strict_correct'] else 'wrong_strict';(good if dec['strict_correct'] else wrong).append(item)
    if not dec['strict_correct']:
     # Compare exactly the same frequency bins/time columns on the shared noise and single-tone controls.
     control_rows=[r for r in records if r['w']==1 and r['meta']['part']=='matched' and r['meta'].get('base')==meta['base'] and (r['meta']['kind']=='noise' or r['meta']['kind']=='one' and r['meta']['offset']==meta['offset'])];control={};cols=np.array(pairs[0]['columns'],int);idx=[int(np.flatnonzero(freq==f)[0]) for f in dec['frequencies']]
     for cr in control_rows:
      with np.load(PARENT/cr['trace']) as cf:
       cy=cf['z'][:,cols].mean(axis=1);peaks=find_peaks(cy)[0];vals={'normalized_values_at_selected_frequencies':cy[idx].tolist(),'is_local_peak_at_selected_frequencies':[bool(i in peaks) for i in idx],'control_anomaly':cr['result']['anomaly'],'control_pair_claim':cr['decisions']['W1D1']['claim']}
       if cr['meta']['kind']=='noise':
        stft_noise=cf['stft'][:,cols];stft_signal=trace_stft[:,cols]-stft_noise;decomp=[]
        for fi in idx:
         pn=float(np.mean(abs(stft_noise[fi])**2));ps=float(np.mean(abs(stft_signal[fi])**2));cross=float(np.mean(2*np.real(stft_noise[fi]*np.conj(stft_signal[fi]))));pt=float(np.mean(power[fi,cols]));np.testing.assert_allclose(pn+ps+cross,pt,rtol=1e-10,atol=1e-15);decomp.append({'raw_noise_power':pn,'raw_added_signal_power':ps,'raw_cross_term':cross,'raw_total_power':pt})
        vals['raw_power_decomposition_at_selected_frequencies']=decomp
       control[cr['meta']['kind']]=vals
     assert set(control)=={'noise','one'};item['same_base_controls']=control;controls.append({'id':row['id'],'controls':control})
   else:
    category='no_accepted_profile' if not profiles else 'fewer_than_two_local_peaks' if not any(len(p['inside_local_peaks_hz'])>=2 for p in profiles) else 'fewer_than_two_above4' if not any(len(p['strong_peak_frequencies_hz'])>=2 for p in profiles) else 'no_interior_bin' if not any(p['spaced_pair_count'] for p in profiles) else 'valley_too_shallow' if not any(p['valid_pair_count'] for p in profiles) else 'unexpected_later_rejection'
    item['category']=category;categories[category]+=1;no.append(item)
    for feature in ['raw_power_two_truth_peaks','normalized_two_truth_peaks','normalized_in_box_two_truth_peaks','normalized_in_box_two_above4','normalized_in_box_pair_passes_valley']:
     present=any(p['frame_counts'].get(feature,0)>0 for p in profiles);frame_evidence[feature]+=present;frame_evidence[category+'__'+feature]+=present
    frame_evidence['mean_full_spectrum_has_truth_pair']+=any(p['mean_full_spectrum_has_truth_pair'] for p in profiles);frame_evidence['mean_inside_box_has_truth_pair']+=any(p['mean_inside_box_has_truth_pair'] for p in profiles)
   cases.append(item)
 assert len(good)==24 and len(wrong)==5 and len(no)==187 and sum(categories.values())==187
 assert categories=={'fewer_than_two_local_peaks':141,'fewer_than_two_above4':8,'valley_too_shallow':38}
 # Rebuild all blind proposals from measured profile gates; compare count and selected frequency pair.
 for row,item in zip(rs,cases):
  has=any(p['valid_pair_count'] for p in item['profiles']);assert has==row['decisions']['W1D1']['claim']
 summary={'status':'PASS','strict_successes':len(good),'wrong_strict':len(wrong),'no_pair':len(no),'no_pair_categories':dict(categories),'no_pair_frame_evidence':dict(frame_evidence),'success_peak_height_range':[min(min(x['selected_pair_metrics']['peak_heights']) for x in good),max(max(x['selected_pair_metrics']['peak_heights']) for x in good)],'success_valley_ratio_range':[min(x['selected_pair_metrics']['valley_ratio'] for x in good),max(x['selected_pair_metrics']['valley_ratio'] for x in good)],'success_error_patterns':dict(Counter(str(x['errors_hz']) for x in good)),'parent_manifest_sha256':sha(PARENT/'manifest.json')}
 save(out/'diagnosis.json',{'summary':summary,'cases':cases});save(out/'wrong_peak_controls.json',controls)
 for name,group in [('correct_pairs.csv',good),('wrong_pairs.csv',wrong),('no_pair_cases.csv',no)]:
  flat=[]
  for r in group:
   f={k:r[k] for k in ['id','base','offset','phase_rad','amplitude_ratio','truth_hz','estimated_hz','errors_hz','strict_correct','loose_correct','category']}
   if 'selected_pair_metrics' in r:f.update(r['selected_pair_metrics'])
   flat.append(f)
  with (out/name).open('w',encoding='utf-8',newline='') as target:
   wr=csv.DictWriter(target,fieldnames=list(flat[0]));wr.writeheader();wr.writerows(flat)
 # Representative cases selected by fixed criteria, not inference retuning.
 representatives=[next(x for x in good if x['offset']==0 and x['phase_rad']==0),next(x for x in good if x['offset']==.5),wrong[0],wrong[-1],next(x for x in no if x['category']=='fewer_than_two_local_peaks'),next(x for x in no if x['category']=='valley_too_shallow')]
 fig,axes=plt.subplots(3,2,figsize=(11,9),layout='constrained')
 for ax,item in zip(axes.flat,representatives):
  row=next(r for r in rs if r['id']==item['id']);p=row['result']['methods']['1']['profiles'][0];c=item['profiles'][0]
  with np.load(PARENT/row['trace']) as trace:freq=trace['freqs'];y=np.array(p['profile'])
  f1,f2=item['truth_hz'];mask=(freq>=f1-5)&(freq<=f2+6);ax.plot(freq[mask],y[mask],'.-',label='Cluster mean Z');ax.axhline(4,color='gray',ls=':',label='Peak height gate=4')
  for f in [f1,f2]:ax.axvline(f,color='black',ls='--',lw=1)
  if item['estimated_hz']:
   for f in item['estimated_hz']:i=int(np.flatnonzero(freq==f)[0]);ax.plot(f,y[i],'ro',mfc='none',ms=9)
  ax.set_title(item['id']+'\n'+item['category'],fontsize=9);ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Mean normalized excess Z');ax.set_ylim(bottom=0)
 axes[0,0].legend(fontsize=8);fig.suptitle('Saved-profile diagnosis: dashed=true carriers; red circles=chosen frequencies\nTruth annotations are evaluation only; mean power peaks do not identify physical sources');fig.savefig(out/'representative_profiles.png',dpi=160);plt.close(fig)
 assert all(sha(ROOT/p)==h for p,h in parent_hashes.items());save(out/'verification.json',{'status':'PASS','cases':216,'profiles_reconstructed':sum(len(r['profiles']) for r in cases),'exclusive_categories_total':187,'claims_match_reconstructed_profile_gates':True,'wrong_control_power_identities':10,'parent_artifact_hashes_unchanged':len(parent_hashes)})
 print(json.dumps(summary));print('REPORT_DIRECTORY',out,flush=True)
if __name__=='__main__':main()
