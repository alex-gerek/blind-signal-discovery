"""Four illustrative saved spectra and exact STFT timing geometry."""
from pathlib import Path
from datetime import datetime,timezone
import csv,gzip,hashlib,json,sys,uuid
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];PARENT=ROOT/'results/matched_gap2/20261128T115035Z_bf957291';TEMP=ROOT/'results/temporal_peak_persistence/20261128T135237Z_633a21cf';PROTOCOL=ROOT/'docs/temporal_spectrum_comparisons_protocol_2026-11-28.md'
CASES=[('stable_true','b00_two_0_0_1',[88,90],'True88/90Hz; equal amplitudes'),('brief_true','b00_two_0_0_0.5',[88,90],'True88/90Hz; weaker partner0.5'),('stable_wrong','b05_two_0_0_0.25',[86,89],'True86/88Hz; wrong reference86/89'),('single_companion','b09_one_0_0_1',[102,105],'Only102Hzinjected; reference102/105')]
def load(p):return json.load(gzip.open(p,'rt',encoding='utf-8')) if p.suffix=='.gz' else json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main():
 out=ROOT/'results/temporal_spectrum_comparisons'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True);print('OUTPUT',out,flush=True);(out/'driver.py').write_bytes(Path(__file__).read_bytes());(out/'protocol.md').write_bytes(PROTOCOL.read_bytes());records=load(TEMP/'records.json.gz');parent_hashes=load(PARENT/'artifact_sha256.json')['files'];temporal_hashes=load(TEMP/'artifact_sha256.json')['files'];assert all(sha(ROOT/p)==h for p,h in {**parent_hashes,**temporal_hashes}.items())
 geometry=[]
 for j in range(25):
  start=j*128;end=start+1024;overlap=max(0,min(end,2304)-max(start,1280));center=.5+j*.125;geometry.append({'frame':j,'center_s':center,'start_sample':start,'last_sample_inclusive':end-1,'start_s':start/1024,'end_s_exclusive':end/1024,'tone_overlap_samples':overlap,'tone_overlap_s':overlap/1024,'central9':1.25<=center<=2.25})
 central=[x for x in geometry if x['central9']];assert len(central)==9 and [g['frame'] for g in central]==list(range(6,15));assert [g['tone_overlap_s'] for g in central]==[.5,.625,.75,.875,1,.875,.75,.625,.5];assert sum(g['tone_overlap_samples']>0 for g in geometry)==15
 with (out/'window_geometry.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(geometry[0]));w.writeheader();w.writerows(geometry)
 prepared=[];frames_csv=[];spectrum_csv=[];summary=[]
 for tag,ident,reference,title in CASES:
  r=next(r for r in records if r['id']==ident);m=r['meta'];truth=[m['f1']]+([m['f2']] if m['kind']=='two' else []);low=min(truth+reference)-5;high=max(truth+reference)+5
  with np.load(PARENT/r['trace']) as a:
   assert a['waveform'].shape==(4096,) and np.array_equal(a['times'],np.array([g['center_s'] for g in geometry]));freq=a['freqs'].copy();times=a['times'].copy();z=a['z'].copy();power=a['power'].copy();psd=a['psd'].copy();tb=a['time_background'].copy()
  mask=(freq>=low)&(freq<=high);indices=[int(np.flatnonzero(freq==x)[0]) for x in reference];frame_metrics=[]
  for g in central:
   j=g['frame'];y=z[:,j];p=power[:,j];zp=find_peaks(y)[0];pp=find_peaks(p)[0];i,k=indices;valley=float(np.min(y[i+1:k]));ratio=valley/min(y[i],y[k]);passes=any(q[:2]==reference for q in r['frames'][j]['pairs']);manual=all(idx in zp and y[idx]>=4 for idx in indices) and k-i>=2 and ratio<=.8;assert passes==manual
   metric={'case':tag,'id':ident,**g,'reference_pair_hz':str(reference),'reference_pair_passes_z_gates':passes,'normalized_height1':float(y[i]),'normalized_height2':float(y[k]),'normalized_valley':valley,'valley_over_smaller_peak':ratio,'raw_power1':float(p[i]),'raw_power2':float(p[k]),'both_reference_bins_local_raw_peaks':all(idx in pp for idx in indices),'both_reference_bins_local_z_peaks':all(idx in zp for idx in indices),'local_raw_peaks_hz':str([float(freq[idx]) for idx in pp if mask[idx]]),'local_z_peaks_hz':str([float(freq[idx]) for idx in zp if mask[idx]])};frames_csv.append(metric);frame_metrics.append(metric)
   for fi in np.flatnonzero(mask):spectrum_csv.append({'case':tag,'frame':j,'center_s':g['center_s'],'frequency_hz':float(freq[fi]),'raw_stft_power':float(p[fi]),'normalized_z':float(y[fi]),'spectral_background':float(psd[fi]),'time_background':float(tb[j])})
  passed=[x['frame'] for x in frame_metrics if x['reference_pair_passes_z_gates']];longest=run=0;last=None
  for j in passed:run=run+1 if last is not None and j==last+1 else 1;longest=max(longest,run);last=j
  s={'case':tag,'id':ident,'truth_hz':truth,'reference_pair_hz':reference,'central_pair_frames':passed,'central_pair_times_s':[geometry[j]['center_s'] for j in passed],'central_pair_count':len(passed),'central_longest_exact_run':longest,'local_frequency_range_hz':[low,high],'meta':m,'frame_metrics':frame_metrics};summary.append(s);prepared.append(dict(tag=tag,ident=ident,title=title,truth=truth,reference=reference,mask=mask,freq=freq,times=times,z=z,power=power,frame_metrics=frame_metrics))
 for name,rows in [('frame_metrics.csv',frames_csv),('local_spectral_samples.csv',spectrum_csv)]:
  with (out/name).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 raw_values=np.concatenate([10*np.log10(np.maximum(d['power'][d['mask']],1e-18)).ravel() for d in prepared]);raw_lo=float(np.quantile(raw_values,.02));raw_hi=float(np.quantile(raw_values,.998));z_hi=float(max(np.max(d['z'][d['mask']]) for d in prepared));log_z_hi=float(np.ceil(np.log10(z_hi)*2)/2)
 save(out/'display_settings.json',{'raw_10log10_power_limits':[raw_lo,raw_hi],'raw_limits_rule':'pooled2nd/99.8thpercentiles across four local displays; display-only clipping','normalized_log10_limits':[-1.,log_z_hi],'raw_power_reference':'arbitrary numerical reference1,not PSD/dBFS','spectral_profiles':'linear exact1Hzbins; common y-axis within each case across nineframes; per-case axes differ','full_spectrogram':'all25frames; frequency crop only','reference_pairs':{d['tag']:d['reference'] for d in prepared}})
 # Common-scale heat maps show both raw local power and normalized representation over all25centers.
 fig,axes=plt.subplots(4,2,figsize=(12,10),layout='constrained')
 for row,d in enumerate(prepared):
  for col in [0,1]:
   ax=axes[row,col];freq=d['freq'][d['mask']];values=d['power'][d['mask']] if col==0 else d['z'][d['mask']];values=10*np.log10(np.maximum(values,1e-18)) if col==0 else np.log10(np.maximum(values,1e-12));im=ax.pcolormesh(d['times'],freq,values,shading='nearest',cmap='viridis',vmin=raw_lo if col==0 else -1,vmax=raw_hi if col==0 else log_z_hi)
   for f in d['truth']:ax.axhline(f,color='white',ls='--',lw=1)
   ax.axvline(1.25,color='white',ls=':',lw=1);ax.axvline(2.25,color='white',ls=':',lw=1);ax.set_xlim(.4375,3.5625);ax.set_ylim(freq.min()-.5,freq.max()+.5);ax.set_xlabel('Frame center (s)');ax.set_ylabel('Frequency (Hz)');ax.set_title(d['title']+(' / raw power' if col==0 else ' / normalized Z'),fontsize=9)
  if row==3:
   fig.colorbar(axes[0,0].collections[0],ax=axes[:,0].tolist(),label='10log10 raw STFT power (reference1 arb.unit)');fig.colorbar(axes[0,1].collections[0],ax=axes[:,1].tolist(),label='log10 normalized excess Z')
 fig.suptitle('Four illustrative records: complete local spectra evolve across overlapping1s windows\nDashed horizontal=true carriers; vertical dotted=central9center interval; shared colors per column');fig.savefig(out/'four_cases_spectral_evolution.png',dpi=170);plt.close(fig)
 # Nine individual spectra preserve the whole local shape, not only two peak heights.
 for d in prepared:
  for mode in ['normalized','raw']:
   fig,axes=plt.subplots(3,3,figsize=(11,7.5),layout='constrained');values=d['z'] if mode=='normalized' else d['power'];ymax=float(np.max(values[d['mask']][:,[g['frame'] for g in central]]))*1.08
   for ax,g,metric in zip(axes.flat,central,d['frame_metrics']):
    j=g['frame'];ax.plot(d['freq'][d['mask']],values[d['mask'],j],'.-',color='#0072b2',lw=1.2)
    for f in d['truth']:ax.axvline(f,color='black',ls='--',lw=.8)
    if mode=='normalized':ax.axhline(4,color='gray',ls=':',lw=.8)
    if metric['reference_pair_passes_z_gates']:
     for f in d['reference']:i=int(np.flatnonzero(d['freq']==f)[0]);ax.plot(f,values[i,j],'o',mfc='none',mec='#d55e00',ms=8)
    ax.set_ylim(0,ymax);ax.set_title(f"t={g['center_s']:.3f}s; input[{g['start_s']:.3f},{g['end_s_exclusive']:.3f})\nreference pair {'passes' if metric['reference_pair_passes_z_gates'] else 'does not pass'} Z gates",fontsize=8);ax.set_xlabel('Frequency (Hz)',fontsize=8);ax.set_ylabel('Z' if mode=='normalized' else 'Raw STFT power',fontsize=8);ax.tick_params(labelsize=8)
   fig.suptitle(d['title']+'\n'+('Normalized spectra; gray=height4; red=reference pair passes' if mode=='normalized' else 'Raw spectra; red annotations reuse Z-gate results,not raw-power detection')+'; dashed=true carriers',fontsize=10);fig.savefig(out/(d['tag']+'_'+mode+'_nine_spectra.png'),dpi=165);plt.close(fig)
 # Explicit geometry: centers within pulse do not imply windows within pulse.
 fig,ax=plt.subplots(figsize=(11,6),layout='constrained');ax.add_patch(Rectangle((0,10.25),4,.5,color='#cccccc'));ax.add_patch(Rectangle((1.25,10.25),1,.5,color='#0072b2'));ax.text(.03,10.9,'Full recording:4s noise; added tone(s) only[1.25,2.25)s',fontsize=10)
 for n,g in enumerate(central):
  y=9-n;ax.broken_barh([(g['start_s'],1)],(y-.28,.56),facecolors='#99ccee',edgecolors='#335577');ax.plot(g['center_s'],y,'ko',ms=4);ax.text(g['end_s_exclusive']+.04,y,f"center{g['center_s']:g}; tone overlap{g['tone_overlap_s']:g}s",va='center',fontsize=8)
 ax.axvspan(1.25,2.25,color='#0072b2',alpha=.10);ax.axvline(1.25,color='#0072b2',ls='--');ax.axvline(2.25,color='#0072b2',ls='--');ax.plot([g['center_s'] for g in geometry],[-.15]*25,'.',color='gray');ax.set_yticks(list(range(1,10)),['frame'+str(g['frame']) for g in reversed(central)]);ax.set_xlim(-.04,4.06);ax.set_ylim(-.5,11.3);ax.set_xlabel('Time within4s record (s)');ax.set_title('Why nine central spectra? 1s window,0.125s hop,selected CENTERS1.25..2.25s\nCentral windows cover[0.75,2.75)s altogether; overlap87.5%,not independent observations');ax.grid(axis='x',alpha=.2);fig.savefig(out/'window_geometry.png',dpi=165);plt.close(fig)
 save(out/'case_summary.json',summary);assert all(sha(ROOT/p)==h for p,h in {**parent_hashes,**temporal_hashes}.items());save(out/'verification.json',{'status':'PASS','cases':4,'central_frames_per_case':9,'all_frame_geometry_checked':25,'local_spectral_samples':len(spectrum_csv),'reference_pair_gate_checks':len(frames_csv),'tone_touching_full_windows':15,'central_union_s':[.75,2.75],'parent_artifact_hashes_unchanged':len(parent_hashes),'temporal_artifact_hashes_unchanged':len(temporal_hashes),'new_waveforms':0})
 print(json.dumps({'cases':[{k:s[k] for k in ['case','truth_hz','reference_pair_hz','central_pair_times_s','central_pair_count','central_longest_exact_run']} for s in summary],'display_raw_limits':[raw_lo,raw_hi],'max_z':z_hi}));print('COMPLETE',out,flush=True)
if __name__=='__main__':main()
