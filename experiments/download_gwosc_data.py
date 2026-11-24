"""Download and verify official GWOSC HDF5 strain datasets for GW150914.

Preserves original downloaded files unchanged under `data/gwosc/` and records
source URLs, SHA-256 checksums, sampling rates, GPS timestamps, finite-value
checks, and GWOSC metadata in `results/ligo_experiment_002/data_provenance.json`.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path
from typing import Any

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent

OFFICIAL_GWOSC_FILES: tuple[tuple[str, int, str], ...] = (
    (
        "H1",
        32,
        "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5",
    ),
    (
        "L1",
        32,
        "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_4KHZ_R1-1126259447-32.hdf5",
    ),
    (
        "H1",
        4096,
        "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126257415-4096.hdf5",
    ),
    (
        "L1",
        4096,
        "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_4KHZ_R1-1126257415-4096.hdf5",
    ),
)


def _decode_scalar(val: Any) -> Any:
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray) and val.ndim == 0:
        return _decode_scalar(val.item())
    if hasattr(val, "item"):
        try:
            return val.item()
        except Exception:
            return str(val)
    return val


def download_and_verify() -> list[dict[str, Any]]:
    data_dir = PROJECT_ROOT / "data" / "gwosc"
    out_dir = PROJECT_ROOT / "results" / "ligo_experiment_002"
    data_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    provenance: list[dict[str, Any]] = []
    for detector, duration_s, url in OFFICIAL_GWOSC_FILES:
        filename = url.rsplit("/", 1)[-1]
        dest = data_dir / filename
        if not dest.exists():
            print(f"[GWOSC] Downloading {filename} ...")
            req = urllib.request.Request(url, headers={"User-Agent": "BSD-Research-Prototype/0.1"})
            with urllib.request.urlopen(req, timeout=240) as resp, dest.open("wb") as fh:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    fh.write(chunk)
        else:
            print(f"[GWOSC] Using existing immutable file {filename}")

        sha256 = hashlib.sha256(dest.read_bytes()).hexdigest()
        with h5py.File(dest, "r") as h5:
            strain = h5["strain/Strain"][:]
            xspacing = float(h5["strain/Strain"].attrs["Xspacing"])
            xstart = float(h5["strain/Strain"].attrs["Xstart"])
            fs = 1.0 / xspacing
            dq_mask = h5["quality/simple/DQmask"][:]
            dq_shortnames = [_decode_scalar(s) for s in h5["quality/simple/DQShortnames"][:]]
            n_finite = int(np.sum(np.isfinite(strain)))
            meta = {k: _decode_scalar(h5["meta"][k][()]) for k in h5["meta"].keys()}

        finite_strain = strain[np.isfinite(strain)]
        record = {
            "detector": detector,
            "duration_s": duration_s,
            "filename": filename,
            "url": url,
            "size_bytes": int(dest.stat().st_size),
            "sha256": sha256,
            "gps_start": xstart,
            "gps_end": xstart + strain.size / fs,
            "sample_rate_hz": fs,
            "n_samples": int(strain.size),
            "n_finite": n_finite,
            "all_finite": bool(n_finite == strain.size),
            "strain_std": float(np.std(finite_strain)) if finite_strain.size else None,
            "strain_min": float(np.min(finite_strain)) if finite_strain.size else None,
            "strain_max": float(np.max(finite_strain)) if finite_strain.size else None,
            "dq_mask_len": int(dq_mask.size),
            "dq_shortnames": dq_shortnames,
            "dq_data_bit_fraction": float(np.mean((dq_mask & 1) > 0)),
            "meta": meta,
        }
        provenance.append(record)
        print(
            f"[GWOSC] Verified {filename}: {record['size_bytes'] / 1e6:.2f} MB | "
            f"fs={fs:.1f} Hz | GPS=[{xstart:.1f}, {record['gps_end']:.1f}] | "
            f"finite={record['all_finite']} | std={record['strain_std']:.3e} | "
            f"DQ_DATA={record['dq_data_bit_fraction'] * 100:.1f}% | sha256={sha256[:16]}..."
        )

    prov_json = json.dumps(provenance, indent=2)
    (data_dir / "provenance.json").write_text(prov_json, encoding="utf-8")
    (out_dir / "data_provenance.json").write_text(prov_json, encoding="utf-8")
    print(f"[GWOSC] Provenance written to {out_dir / 'data_provenance.json'}")
    return provenance


if __name__ == "__main__":
    download_and_verify()
