"""Download and verify official GWOSC HDF5 datasets for BSD Experiment 003.

Downloads 32-second and 4096-second (4096 Hz) HDF5 strain files for:
- Independent validation events: GW151012, GW151226, GW170104, GW170608,
  GW170814 (with V1), GW170817 (with V1).
- Previously examined reference: GW150914 (reused from data/gwosc/ if present).

Also writes the untouched holdout manifest (`v02_holdout_manifest.json`) for
future BSD v0.2 comparison without downloading or inspecting those events.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent

EXPERIMENT_003_EVENTS: tuple[dict[str, Any], ...] = (
    {
        "name": "GW150914",
        "version": "v3",
        "role": "reference_previously_examined",
        "observing_run": "O1",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1"),
    },
    {
        "name": "GW151012",
        "version": "v3",
        "role": "independent_validation",
        "observing_run": "O1",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1"),
    },
    {
        "name": "GW151226",
        "version": "v2",
        "role": "independent_validation",
        "observing_run": "O1",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1"),
    },
    {
        "name": "GW170104",
        "version": "v2",
        "role": "independent_validation",
        "observing_run": "O2",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1"),
    },
    {
        "name": "GW170608",
        "version": "v3",
        "role": "independent_validation",
        "observing_run": "O2",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1"),
    },
    {
        "name": "GW170814",
        "version": "v3",
        "role": "independent_validation",
        "observing_run": "O2",
        "source_class": "BBH",
        "detectors_requested": ("H1", "L1", "V1"),
    },
    {
        "name": "GW170817",
        "version": "v3",
        "role": "independent_validation",
        "observing_run": "O2",
        "source_class": "BNS",
        "detectors_requested": ("H1", "L1", "V1"),
    },
)

# Strictly untouched GWTC-1 holdout set reserved for future BSD v0.2 comparison
V02_HOLDOUT_EVENTS: tuple[dict[str, Any], ...] = (
    {"name": "GW170729", "version": "v1", "observing_run": "O2", "catalog": "GWTC-1-confident"},
    {"name": "GW170809", "version": "v1", "observing_run": "O2", "catalog": "GWTC-1-confident"},
    {"name": "GW170818", "version": "v1", "observing_run": "O2", "catalog": "GWTC-1-confident"},
    {"name": "GW170823", "version": "v1", "observing_run": "O2", "catalog": "GWTC-1-confident"},
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


def _fetch_event_json(name: str, version: str) -> dict[str, Any]:
    url = f"https://gwosc.org/eventapi/json/GWTC-1-confident/{name}/{version}/"
    req = urllib.request.Request(url, headers={"User-Agent": "BSD-Research-Prototype/0.1"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    ev = list(data["events"].values())[0]
    ev["_api_url"] = url
    return ev


def _select_strain_file(
    strain_entries: list[dict[str, Any]],
    detector: str,
    duration_s: int,
) -> dict[str, Any] | None:
    matches = [
        s
        for s in strain_entries
        if s.get("detector") == detector
        and s.get("format") == "hdf5"
        and int(s.get("sampling_rate", 0)) == 4096
        and int(s.get("duration", 0)) == duration_s
    ]
    if not matches:
        return None
    # Prefer standard R1 release if multiple releases exist
    r1 = [m for m in matches if "_R1-" in m.get("url", "")]
    return r1[0] if r1 else matches[0]


def _download_and_verify_one(
    event_name: str,
    detector: str,
    duration_s: int,
    url: str,
    dest_dir: Path,
) -> dict[str, Any]:
    filename = url.rsplit("/", 1)[-1]
    # Reuse existing GW150914 files in data/gwosc/ if already downloaded
    legacy_path = PROJECT_ROOT / "data" / "gwosc" / filename
    dest = legacy_path if legacy_path.is_file() else (dest_dir / filename)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if not dest.is_file():
        print(f"[GWOSC-003] Downloading {event_name} {detector} {duration_s}s ({filename}) ...")
        req = urllib.request.Request(url, headers={"User-Agent": "BSD-Research-Prototype/0.1"})
        tmp_dest = dest.with_suffix(".hdf5.part")
        with urllib.request.urlopen(req, timeout=300) as resp, tmp_dest.open("wb") as fh:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)
        tmp_dest.replace(dest)
    else:
        print(f"[GWOSC-003] Reusing verified file {event_name} {detector} {duration_s}s ({dest.name})")

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

    finite_vals = strain[np.isfinite(strain)]
    return {
        "event": event_name,
        "detector": detector,
        "duration_s": duration_s,
        "filename": filename,
        "relative_local_path": str(dest.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "url": url,
        "size_bytes": int(dest.stat().st_size),
        "sha256": sha256,
        "gps_start": xstart,
        "gps_end": xstart + strain.size / fs,
        "sample_rate_hz": fs,
        "n_samples": int(strain.size),
        "n_finite": n_finite,
        "all_finite": bool(n_finite == strain.size),
        "strain_std": float(np.std(finite_vals)) if finite_vals.size else None,
        "dq_data_bit_fraction": float(np.mean((dq_mask & 1) > 0)),
        "dq_shortnames": dq_shortnames,
        "meta": meta,
    }


def main() -> None:
    out_dir = PROJECT_ROOT / "results" / "ligo_experiment_003"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Write untouched v0.2 holdout manifest
    holdout_manifest = {
        "purpose": "Reserved untouched real-data evaluation set for future BSD v0.2 comparison.",
        "status": "UNINSPECTED_AND_NOT_DOWNLOADED",
        "policy": (
            "These confirmed GWTC-1 events were deliberately excluded from Experiment 003 "
            "and must not be downloaded, inspected, or tuned against during v0.2 development."
        ),
        "events": list(V02_HOLDOUT_EVENTS),
    }
    (out_dir / "v02_holdout_manifest.json").write_text(
        json.dumps(holdout_manifest, indent=2), encoding="utf-8"
    )
    print(f"[GWOSC-003] Wrote untouched v0.2 holdout manifest ({len(V02_HOLDOUT_EVENTS)} events).")

    # 2. Query official GWOSC API and build download jobs
    event_catalog: list[dict[str, Any]] = []
    jobs: list[tuple[str, str, int, str, Path]] = []

    for spec in EXPERIMENT_003_EVENTS:
        ev_json = _fetch_event_json(spec["name"], spec["version"])
        strain_list = ev_json.get("strain", [])
        avail_dets: list[str] = []
        for det in spec["detectors_requested"]:
            s32 = _select_strain_file(strain_list, det, 32)
            s4096 = _select_strain_file(strain_list, det, 4096)
            if s32 is not None and s4096 is not None:
                avail_dets.append(det)
                dest_dir = PROJECT_ROOT / "data" / "gwosc" / spec["name"]
                jobs.append((spec["name"], det, 32, s32["url"], dest_dir))
                jobs.append((spec["name"], det, 4096, s4096["url"], dest_dir))

        event_catalog.append(
            {
                "name": spec["name"],
                "version": spec["version"],
                "role": spec["role"],
                "observing_run": spec["observing_run"],
                "source_class": spec["source_class"],
                "published_gps": float(ev_json["GPS"]),
                "network_matched_filter_snr": ev_json.get("network_matched_filter_snr"),
                "mass_1_source": ev_json.get("mass_1_source"),
                "mass_2_source": ev_json.get("mass_2_source"),
                "chirp_mass_source": ev_json.get("chirp_mass_source"),
                "api_url": ev_json["_api_url"],
                "detectors_available": avail_dets,
            }
        )

    print(f"[GWOSC-003] Total HDF5 files to verify/download: {len(jobs)}")
    provenance_files: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        future_map = {
            pool.submit(_download_and_verify_one, *job): job for job in jobs
        }
        for fut in as_completed(future_map):
            rec = fut.result()
            provenance_files.append(rec)
            print(
                f"[GWOSC-003] OK {rec['event']} {rec['detector']} {rec['duration_s']}s | "
                f"{rec['size_bytes']/1e6:.1f} MB | finite={rec['all_finite']} | "
                f"DQ_DATA={rec['dq_data_bit_fraction']*100:.1f}%"
            )

    provenance_files.sort(key=lambda r: (r["event"], r["detector"], r["duration_s"]))
    payload = {
        "events": event_catalog,
        "files": provenance_files,
    }
    (out_dir / "data_provenance.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[GWOSC-003] Saved provenance to {out_dir / 'data_provenance.json'}")


if __name__ == "__main__":
    main()
