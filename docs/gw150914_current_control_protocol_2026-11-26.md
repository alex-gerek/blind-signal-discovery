# GW150914 current-code real-data control protocol (2026-11-26)

Specified before execution. This is a repeatability/control analysis of a previously examined O1 epoch, NOT a blind discovery or independent confirmation of GW150914.

## Source and data

- Current source:3497211e0f19832a5e1bd268971e9e3cbc0b83d6. Detector, preprocessing and calibration code remain unchanged; new protocol/driver are uncommitted and hashed.
- Use ONLY the two existing official GWTC-1/GW150914 v3 O1 4096Hz4096s HDF5 files:
  H-H1_GWOSC_4KHZ_R1-1126257415-4096.hdf5 and L-L1_GWOSC_4KHZ_R1-1126257415-4096.hdf5.
- Verify each raw file against results/ligo_experiment_003/data_provenance.json. Save source URLs, SHA256, dataset metadata and DQ/injection flag definitions. No download, O3/holdout manifest, other event or external input.
- Official event page:https://gwosc.org/events/GW150914/; known reference GPS1126259462.422 from the existing project reference. It is known before this experiment, so no blind-analysis claim.
- Data provider: Gravitational Wave Open Science Center, LIGO Scientific Collaboration and Virgo Collaboration (current GWOSC service also includes KAGRA). Data are CC BY4.0:https://gwosc.org/acknowledgement/. O1/O2 data reference:R.Abbott et al., SoftwareX13(2021)100658, doi:10.1016/j.softx.2021.100658. Original/derived input arrays stay outside git; report source attribution and transformations explicitly.

## Fixed non-overlapping blocks (same GPS for H1/L1)

| Role | Core GPS interval | Four-second windows |
| --- | --- | --- |
| PSD fit |1126257416..1126257656 |60 |
| Calibration |1126257680..1126258960 |320 |
| Target |1126259304..1126259624 |80 |
| Independent background control |1126260104..1126261384 |320 |

Every window uses1s padding on both sides. Verify non-overlap of BLOCK padded footprints, not only core timestamps. Adjacent windows within a block share padding and need not be statistically independent. Keep the raw DATA flag requirement on the whole padded span, preserve all other quality/hardware-injection annotations. Require full expected block counts and contiguous cores; fail rather than silently stitch dropped windows.

## Fixed preprocessing and detection

- Existing adapted_ligo_whitened: scale1e18; median-Welch PSD fitted ONLY to the60-window PSD block; zero-phase4th-order30–380Hz Butterworth; antialiased4096->1024Hz decimation;1s contexts trimmed; no injected template, amplitude adjustment, per-target RMS normalization, new veto or threshold tuning.
- The30–380Hz filter has smooth attenuation; this is not a hard search-bin mask. Existing detector frequency normalization can nominate attenuated out-of-band bins. Preserve metadata and distinguish display limits from searched bins.
- Primary: frozen MultiResolutionDetector, FAP0.01, norm_fraction0.30, calibration_seed0;320cal windows ->96normalization/224threshold scores. Freeze before target/control. Save scores, thresholds, empirical p, candidates, GPS/time/frequency and all non-detections.
- Secondary: actual Python Workstation transient/causal_trailing FAP0.01/calibration_seed0/vetoes off on the SAME preprocessed calibration and80target windows. Compare used availability/thresholds and event-window decision; do not optimize or claim known adaptive limitations resolved.
- Initial on-source reference window is index39, GPS1126259460..1126259464. Report its nomination and whether any nominated candidate overlaps reference +/-0.5s and35–250Hz. Temporal/band overlap is a control association, not independent physical confirmation.
- Off-source target evaluation excludes all windows whose PADDED footprints intersect reference +/-16s. Report the known-event window separately from guard-zone windows. Background-control nominations are not proven false positives: real data can contain other signals/glitches.
- No coincidence significance or astrophysical FAR is estimated from this short two-site check. Single-detector empirical p is not an astrophysical probability.

## Preservation, verification and operator follow-up

Fresh UTC/UUID results/gw150914_current_control directory. Save protocol/driver/source hashes, input/block/DQ metadata, frozen calibration, full-precision frozen target/control CSV and candidate records, current Workstation JSON/CSV, per-site outcomes and figures. Save partial failures; no seed/parameter replacement.

Generate GPS,value CSVs for calibration and target under ignored data/derived/gw150914_workstation/<run-id> (large CC BY data derivatives stay outside git). Preserve exact hashes and a manifest in the result bundle. These files contain already whitened/bandpassed data; do not apply the LIGO preprocessing again. NPY companions preserve exact array bytes if CSV parsing differs. User Workstation ingestion requires a separate read-only review of its reserved-manifest guard before being used under the no-manifest-read constraint; do not bypass that guard silently.

Verify CSV/NPY target/calibration arrays, frozen empirical-p/threshold decision consistency, per-window GPS mapping, availability and exports. Independently check source/raw hashes and report actual outcomes. User can inspect plotted real waveforms/spectrograms first; no new physical discovery follows.
