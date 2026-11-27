# Research checkpoint — 27 November 2026

Owner approved committing and pushing the completed research stages to development/v0.3 after review. Production mathematics, app, unrelated test cleanup and article drafts are excluded. No main merge, release, Zenodo change or O3 access.

## Preserved stages

- [Corrected warning-impact control](WARNING_IMPACT_CONTROL_FINDINGS_2026-11-27.md): warning differences, pure gain invariance, line-noise alarms and negative findings.
- [Two-oscillation anomaly pilot](TWO_OSCILLATION_ANOMALY_FINDINGS_2026-11-27.md): separate anomaly detection, cluster morphology and source ambiguity.
- Saved-map and structure-definition supplements under results/two_oscillation_structure_review/ and results/two_oscillation_precluster_inspection/.
- [Blind two-tone exploratory comparison](BLIND_TWO_TONE_EXPLORATORY_FINDINGS_2026-11-27.md): four arms, successes and wrong pairs, independent controls, calibration, timing, runtime and retained initial fixture failure.
- Five gap-2Hz successes versus67unsuccessful cases under results/blind_two_tone_gap2_review/. All four stricter successes use equal amplitudes and on-grid carriers; descriptive, not causal evidence.

All new waveform studies use pinned source cc68b93fee699307e3314fce1bd4c9fe403c1e6e. To reproduce, use a separate disposable checkout at that revision, copy the corresponding new scripts into experiments/ and their fixed protocol into docs/, then run the driver with python -B. Its source guard deliberately rejects a newer revision. Dependencies from earlier studies already exist at the pin; copied dependencies/source hashes are retained in run manifests. Use a NEW output folder, never overwrite archived results. Reproduction reruns are research steps requiring a new journal entry. Final reviewed reports can contain qualifications beyond the generic report builders.

## Restoring large JSON results

GitHub rejects individual files over100MiB. The completed blind pilot preserves three exact JSON files as lossless .json.gz archives, with original and archive SHA256 in publication_archives.json. Uncompressed originals remain local and are specifically gitignored; no evidence was removed. From the repository root, restore them before running verifiers/report builders:

```python
from pathlib import Path
import gzip, hashlib, json, shutil
run = Path('results/blind_two_tone_exploratory/20261127T192047Z_f279641e')
for item in json.loads((run / 'publication_archives.json').read_text())['archives']:
    source = run / item['archive']
    target = run / item['file']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == item['archive_sha256']
    if not target.exists():
        with gzip.open(source, 'rb') as src, target.open('xb') as dst:
            shutil.copyfileobj(src, dst)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == item['original_sha256']
```

The1,552NPZ trace maps (~508.5MiB) remain in the persistent local traces/ folder, excluded by the existing *.npz rule. Their exact hashes and paths are in artifact_sha256.json. The public checkout contains numerical results, candidates, diagnostic profiles, figures, generator/driver and source snapshots, but NOT those binary maps. Independent map verification requires those local traces or deterministic regeneration with the archived source and parameters; restoring JSON alone does not restore the maps. This is a local preservation checkpoint, not a remotely backed-up trace archive or Zenodo deposit.

## Historical hash notes

The initial two-oscillation artifact index includes the earlier version of docs/TWO_OSCILLATION_ANOMALY_FINDINGS_2026-11-27.md. A later recorded supplement appended a follow-up link to that living findings document. Its old digest therefore describes the earlier document, not the current one; measurements and protocol did not change. Do not silently rewrite the old index. The publication verification records this specific documentation supersession and hashes the current file separately.

## Next session

Review the five close-pair successes and failures as development evidence; decide whether a matched follow-up should isolate grid position, phase, amplitude ratio or normalization. Keep production unchanged. Replicate calibration and broaden null controls before any adoption claim. No source waveforms or music stems have been recovered.
