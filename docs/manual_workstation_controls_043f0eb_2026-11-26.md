# User-operated Python Workstation controls (2026-11-26)

The owner ran six synthetic controls in the actual local Streamlit interface at source 043f0eb. The desktop shortcut points to the correct launcher, but an old server initially retained an outdated imported module and failed to import build_workstation_report. A fresh-process import succeeded. Only the verified old port-8501 Workstation process family was restarted; server health and a fresh AppTest startup passed. This was a running-session issue, not a missing function in the current file. No detector code was changed.

## Settings and protocol adjustment

Transient, causal_trailing, FAP=0.01, calibration_seed=0, vetoes off, 220 background windows, 80 initial noise-control windows, 10 linear_chirp injections in pink_colored noise; seeds 2026/2027/2028 and root-energy SNR 8/18. Each run contains 100 windows: 80 initial controls, 10 injected windows and 10 interleaved noise windows.

The earlier headless protocol used 200 initial controls. The current UI maximum is 80; this adjustment was stated before the manual runs. The manual matrix is not an exact replay of that 200-control pilot, and its counts must remain separate. Settings were not changed to obtain a preferred validation status or detection rate.

## Saved measurements

| Seed | SNR8 detected | SNR18 detected | Initial noise FP / available, counted once | Calibration |
| --- | --- | --- | --- | --- |
| 2026 | 3/10 | 10/10 | 1/80 | PROVISIONAL |
| 2027 | 5/10 | 10/10 | 2/80 | PROVISIONAL |
| 2028 | 4/10 | 10/10 | 2/80 | VALIDATED |

No injection abstentions occurred. Across the three seeds: SNR8 12/30, SNR18 30/30. Shared pre-injection noise is counted once per seed: 5/240 = 2.0833%, 960 s available exposure. These are descriptive measurements from three frozen calibrations, not a verified population FAP or universal sensitivity curve. VALIDATED is internal screening, not a guarantee of exact observed 1% or physical significance.

## Direct verification

Six JSON/CSV export bundles exist locally under ignored results/workstation_runs/. All 600 exported available decisions were compared with their full-precision score/threshold crossing; JSON parsed without nonfinite constants. Counts agree with initial-control and injection accounting. The 80 noise scores, thresholds, availability and decisions match within each paired seed across SNR. File hashes and exact run IDs are preserved in [the verification record](validation/manual_workstation_controls_043f0eb_20261126.json). The legacy synthetic target_sha256 fields are descriptive source labels, not cryptographic input-array fingerprints; the saved file hashes must not be confused with input hashes.

User screenshots additionally show the seed2026/SNR18 window81 chirp and noise false alarm window4. Those screenshot observations do not establish all-input morphology or real-sensor validity. Original screenshots remain conversation attachments; this record does not invent a repository image artifact.

## Glitch warning and exposure distinction

Seed2028/SNR18 has GLITCH_STORM warnings on windows 97, 98 and 99. Windows 97 and 99 contain known injections, 98 is interleaved noise. All remain detection_available=True; two are detected and window98 is not. is_clean_calibrated explicitly excludes glitch-warning windows, yielding clean duty 97%, while detection-available duty remains 100%. This is consistent with current code and exports. The flag is a diagnostic, not proof of an instrumental glitch. A response to the late injection sequence is plausible, but no causal attribution is claimed without a separate controlled comparison. Union remains unavailable in all six runs.

No O3, detector tuning, historical-result overwrite, commit, push, release or Zenodo change was involved. This closes a bounded user-operated synthetic/export check; native browser WAV decoding, broader GUI workflows, adaptive historical reruns and real-data validation remain open.


## Evidence preservation for GitHub (2026-11-26)

On owner authorization, all six selected verified bundles were copied byte-for-byte to [the tracked research archive](../results/manual_workstation_controls_043f0eb_20261126/README.md), with hashes checked against the original verification record. Original ignored local bundles remain untouched. This supersedes the earlier local-only distribution limitation; no new measurement or metadata identity is fabricated.
