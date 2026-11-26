"""Current-source GW150914 O1 positive control with frozen and Workstation decisions."""
from __future__ import annotations
import csv
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib,json,os
from pathlib import Path
import platform,subprocess,sys,time,uuid

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"));sys.dont_write_bytecode=True
RAW_NAMES={"H-H1_GWOSC_4KHZ_R1-1126257415-4096.hdf5","L-L1_GWOSC_4KHZ_R1-1126257415-4096.hdf5"}
DERIVED_PREFIX=str(ROOT/"data/derived/gw150914_workstation").replace("\\","/").lower()+"/"
RAW_PREFIX=str(ROOT/"data").replace("\\","/").lower()+"/"
def guard(event,args):
    if event=="open" and args and isinstance(args[0],(str,bytes,os.PathLike)):
        text=os.path.abspath(os.fsdecode(args[0])).replace("\\","/").lower()
        if "o3_holdout" in text:raise PermissionError("Reserved manifest/data access blocked")
        if text.startswith(RAW_PREFIX) and not text.startswith(DERIVED_PREFIX):
            allowed={str(ROOT/"data/gwosc"/n).replace("\\","/").lower() for n in RAW_NAMES}
            if text not in allowed:raise PermissionError("Non-allowlisted data access blocked")
            mode=args[1] if len(args)>1 else None;flags=args[2] if len(args)>2 else 0
            if (isinstance(mode,str) and any(c in mode for c in "wax+")) or (isinstance(flags,int) and flags&(os.O_WRONLY|os.O_RDWR|os.O_TRUNC|os.O_CREAT)):
                raise PermissionError("Original data writes blocked")
    if event in ("socket.connect","socket.getaddrinfo"):raise PermissionError("Network blocked during measurement")
sys.addaudithook(guard)

import h5py
import numpy as np
import scipy
from bsd.detector_v2 import MultiResolutionDetector
from bsd.ligo_adapter import load_gwosc_hdf5,extract_clean_segments,fit_background_model,preprocess_segment
from bsd.io_loader import LoadedTimeSeries,PreprocessingSummary
from bsd.workstation_runner import WorkstationConfig,run_workstation_analysis,export_workstation_run

BASE="3497211e0f19832a5e1bd268971e9e3cbc0b83d6"
TC=1126259462.422
BLOCKS={"psd":(1126257416.,1126257656.),"calibration":(1126257680.,1126258960.),
    "target":(1126259304.,1126259624.),"control":(1126260104.,1126261384.)}
PROTOCOL=ROOT/"docs/gw150914_current_control_protocol_2026-11-26.md"
def git(*a):return subprocess.check_output(["git",*a],cwd=ROOT,text=True).strip()
def sha(b):return hashlib.sha256(b).hexdigest()
def filehash(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def encode(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x.relative_to(ROOT))
    raise TypeError(type(x).__name__)
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=encode),encoding="utf-8")
def csv_save(p,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sources():
    paths=list((ROOT/"src/bsd").glob("*.py"))+[ROOT/"app/bsd_workstation.py",ROOT/"pyproject.toml"]
    return {str(p.relative_to(ROOT)).replace("\\","/"):filehash(p) for p in paths}
ACTIVE_OUT=None
def main():
    global ACTIVE_OUT
    assert git("rev-parse","HEAD")==BASE
    assert not git("diff",BASE,"--","src","app","tests","pyproject.toml")
    rid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    out=ROOT/"results/gw150914_current_control"/rid;out.mkdir(parents=True,exist_ok=False);ACTIVE_OUT=out
    derived=ROOT/"data/derived/gw150914_workstation"/rid;derived.mkdir(parents=True,exist_ok=False)
    (out/"protocol.md").write_bytes(PROTOCOL.read_bytes());(out/"driver.py").write_bytes(Path(__file__).read_bytes())
    provenance=json.loads((ROOT/"results/ligo_experiment_003/data_provenance.json").read_text(encoding="utf-8"))
    manifest={"run_id":rid,"source_commit":BASE,"status":"RUNNING","source_sha256":sources(),
        "driver_sha256":filehash(Path(__file__)),"protocol_sha256":filehash(PROTOCOL),
        "dirty_paths":git("status","--porcelain").splitlines(),"python":platform.python_version(),
        "numpy":np.__version__,"scipy":scipy.__version__,"h5py":h5py.__version__,
        "known_event_gps":TC,"blocks":BLOCKS,"derived_data_directory":str(derived.relative_to(ROOT)).replace("\\","/"),
        "scope":"Previously used O1 real-data control, not blind discovery; no injected signals",
        "reserved_access":"BLOCKED","network":"BLOCKED","sites":[]}
    save(out/"manifest.json",manifest);print("RUN_DIRECTORY",out,flush=True)
    spans=sorted((a-1,b+1,name) for name,(a,b) in BLOCKS.items())
    assert all(b<c for (a,b,n),(c,d,m) in zip(spans,spans[1:]))
    for site in ("H1","L1"):
        started=time.monotonic();dest=out/site;dest.mkdir()
        data_entry=next(r for r in provenance["files"] if r["event"]=="GW150914" and r["detector"]==site and r["duration_s"]==4096)
        path=ROOT/data_entry["relative_local_path"]
        assert path.name in RAW_NAMES and filehash(path)==data_entry["sha256"]
        rec=load_gwosc_hdf5(path)
        assert rec.detector==site and rec.sample_rate==4096. and rec.n_nonfinite==0
        print("SITE_START",site,"raw_hash_verified",flush=True)
        with h5py.File(path,"r") as f:
            inj=np.asarray(f["quality/injections/Injmask"][:],dtype=np.int64)
            inj_names=[v.decode() if isinstance(v,bytes) else str(v) for v in f["quality/injections/InjShortnames"][:]]
        all_segments=extract_clean_segments(rec,4.,1.,require_dq_data_bit=True)
        blocks={}
        for name,(a,b) in BLOCKS.items():
            seq=[s for s in all_segments if a<=s.gps_start and s.gps_start+4<=b]
            expected=int((b-a)/4);assert len(seq)==expected,(site,name,len(seq),expected)
            assert all(s.gps_start==a+4*i for i,s in enumerate(seq))
            blocks[name]=seq
        model=fit_background_model(site,blocks["psd"],mode="adapted_ligo_whitened",target_sample_rate=1024.,bandpass_hz=(30.,380.))
        observations={name:[preprocess_segment(s,f"{site}_{name}_{i:04d}",model) for i,s in enumerate(seq)]
            for name,seq in blocks.items() if name!="psd"}
        assert all(len(o.signal)==4096 and np.isfinite(o.signal).all() for seq in observations.values() for o in seq)
        quality=[]
        for name,seq in blocks.items():
            for i,s in enumerate(seq):
                first=int(s.gps_start-1-rec.gps_start);last=int(s.gps_start+5-rec.gps_start)
                quality.append({"block":name,"window_index":i,"gps_start":s.gps_start,
                    "dq_mask_and":int(np.bitwise_and.reduce(rec.dq_mask[first:last])),
                    "injection_mask_and":int(np.bitwise_and.reduce(inj[first:last]))})
        csv_save(dest/"quality_windows.csv",quality)
        save(dest/"input_provenance.json",{"original":data_entry,"dq_shortnames":rec.dq_shortnames,"injection_shortnames":inj_names,
            "preprocessing":asdict(model),"block_padded_spans":spans,
            "input_array_sha256":{name:sha(np.stack([o.signal for o in seq]).astype("<f8").tobytes()) for name,seq in observations.items()}})
        det=MultiResolutionDetector();cal=det.calibrate(observations["calibration"],target_fap=.01,norm_fraction=.30,calibration_seed=0)
        save(dest/"frozen_calibration.json",cal.to_summary_dict())
        save(dest/"frozen_null_scores.json",cal.sorted_global_null_scores)
        frozen={}
        for name in ("target","control"):
            rows=[];candidates=[]
            for i,obs in enumerate(observations[name]):
                result=det.analyze(obs,cal)
                p,_=cal.empirical_p_value(result.max_global_score)
                assert result.detected==(result.max_global_score>cal.global_threshold)
                if result.detected:assert p<=.01
                gps=blocks[name][i].gps_start
                rows.append({"window_index":i,"gps_start":gps,"gps_end":gps+4,
                    "score":result.max_global_score,"threshold":cal.global_threshold,"empirical_p":p,"detected":result.detected})
                for c in result.candidates:
                    candidates.append({**c.to_dict(),"window_index":i,"window_gps_start":gps,
                        "gps_t_start":gps+c.t_start,"gps_t_end":gps+c.t_end,"gps_peak":gps+c.peak_time})
            csv_save(dest/f"frozen_{name}_timeline.csv",rows);save(dest/f"frozen_{name}_candidates.json",candidates)
            frozen[name]={"rows":rows,"candidates":candidates}
        # Save physically adapted data locally for later operator inspection; no original arrays in git.
        file_records={}
        arrays={}
        for name in ("calibration","target"):
            a=np.concatenate([o.signal for o in observations[name]])
            arrays[name]=a
            npy=derived/f"GW150914_{site}_{name}_whitened_1024Hz.npy";np.save(npy,a,allow_pickle=False)
            csv_path=derived/f"GW150914_{site}_{name}_whitened_1024Hz.csv"
            times=BLOCKS[name][0]+np.arange(len(a),dtype=np.float64)/1024.
            np.savetxt(csv_path,np.column_stack((times,a)),delimiter=",",header="gps,strain_whitened",comments="",fmt="%.17g")
            file_records[name]={"csv":str(csv_path.relative_to(ROOT)).replace("\\","/"),"csv_sha256":filehash(csv_path),
                "npy":str(npy.relative_to(ROOT)).replace("\\","/"),"npy_sha256":filehash(npy),"samples":len(a)}
            assert np.array_equal(np.load(npy,allow_pickle=False),a)
        save(dest/"derived_file_manifest.json",file_records)
        prep=PreprocessingSummary("derived_real_GWOSC","whitened_1d","dimensionless_whitened_strain",
            4096.,1024.,True,"scipy IIR zero-phase decimate q=4",(30.,380.),(30.,380.),True,
            ("Raw SHA verified","Independent60-window PSD","Windowed1s padded whitening and30-380Hz bandpass","4096to1024Hz","No injected signal"))
        def wrap(name):
            a=arrays[name]
            return LoadedTimeSeries(f"GW150914_{site}_{name}_preprocessed","derived_real_GWOSC",
                file_records[name]["csv_sha256"],len(a),1024.,len(a)/1024.,a,prep,BLOCKS[name],"GPS",
                {"original_file_sha256":data_entry["sha256"],"run_id":rid})
        adaptive=run_workstation_analysis(wrap("target"),wrap("calibration"),WorkstationConfig())
        export_workstation_run(adaptive,dest/"workstation_target")
        assert adaptive.windows[0].stft_threshold==cal.global_threshold
        on=next(i for i,r in enumerate(frozen["target"]["rows"]) if r["gps_start"]<=TC<r["gps_end"])
        assert on==39
        matched=[c for c in frozen["target"]["candidates"] if c["window_index"]==on
            and c["gps_t_start"]<=TC+.5 and c["gps_t_end"]>=TC-.5 and c["f_min"]<=250 and c["f_max"]>=35]
        off=[r for r in frozen["target"]["rows"] if r["gps_end"]+1<=TC-16 or r["gps_start"]-1>=TC+16]
        aw=adaptive.windows[on]
        summary={"site":site,"frozen_threshold":cal.global_threshold,"calibration_status":adaptive.calibration_audit.to_dict(),
            "event_window_index":on,"event_window_gps":frozen["target"]["rows"][on]["gps_start"],
            "frozen_event_window":frozen["target"]["rows"][on],"frozen_reference_associated_candidates":len(matched),
            "frozen_control_nominations":sum(r["detected"] for r in frozen["control"]["rows"]),"frozen_control_windows":len(frozen["control"]["rows"]),
            "frozen_offsource_nominations":sum(r["detected"] for r in off),"frozen_offsource_windows":len(off),
            "adaptive_event":{"available":aw.detection_available,"score":aw.stft_max_score,"threshold":aw.stft_threshold,
                "detected":bool(aw.detection_available and aw.stft_nominated_post_veto),"credibility":aw.credibility_status,"warnings":aw.warning_causes},
            "adaptive_summary":adaptive.summary_metrics,"elapsed_s":time.monotonic()-started}
        save(dest/"summary.json",summary);manifest["sites"].append(summary);save(out/"manifest.json",manifest)
        print("SITE_DONE",site,"frozen_detected",summary["frozen_event_window"]["detected"],"adaptive_detected",summary["adaptive_event"]["detected"],flush=True)
    assert sources()==manifest["source_sha256"]
    manifest["source_hashes_unchanged"]=True;manifest["status"]="COMPLETE";save(out/"manifest.json",manifest)
    print("RUN_STATUS COMPLETE",flush=True)
if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        if ACTIVE_OUT is not None and (ACTIVE_OUT/"manifest.json").exists():
            record=json.loads((ACTIVE_OUT/"manifest.json").read_text(encoding="utf-8"))
            record["status"]="ERROR"
            record["error"]={"type":type(exc).__name__,"message":str(exc)}
            save(ACTIVE_OUT/"manifest.json",record)
        raise
