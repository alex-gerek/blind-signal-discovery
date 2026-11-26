"""Recreate all stationary inputs from recorded seeds without rerunning detection."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"src"))
from bsd.generator import make_split_rng,generate_clean_signal,_scale_to_energy_snr
from bsd.noise import generate_noise
def sha(a):
    return hashlib.sha256(np.asarray(a,dtype="<f8",order="C").tobytes(order="C")).hexdigest()
out=Path(__file__).resolve().parent
checks=[]
for seed in (2026,2027,2028):
    rng=make_split_rng(seed,"workstation_self_test")
    bg=np.concatenate([generate_noise("pink_colored",4096,1024.,rng,sigma=1.)[0] for _ in range(220)])
    noise=[generate_noise("pink_colored",4096,1024.,rng,sigma=1.)[0] for _ in range(220)]
    clean=[generate_clean_signal("linear_chirp",4096,1024.,rng)[0] for _ in range(10)]
    for snr in (8,18):
        chunks=[a.copy() for a in noise]
        for k,c in enumerate(clean): chunks[201+2*k]+=_scale_to_energy_snr(c,1.,snr)
        target=np.concatenate(chunks)
        rec=json.loads((out/f"seed{seed}_snr{snr}"/"array_provenance.json").read_text(encoding="utf-8"))
        assert sha(bg)==rec["background"]["sha256"]
        assert sha(target)==rec["target"]["sha256"]
        checks.append({"seed":seed,"snr":snr,"background_sha256":sha(bg),"target_sha256":sha(target),
            "noise_control_sha256":sha(np.concatenate(noise[:200])),"status":"MATCHED"})
p=out/"regenerated_input_hashes.json"
assert not p.exists()
p.write_text(json.dumps({"status":"ALL_SIX_INPUTS_MATCHED","checks":checks},indent=2),encoding="utf-8")
print("Exact regenerated input hashes: 6/6 MATCHED; paired pre-injection noise identical")
