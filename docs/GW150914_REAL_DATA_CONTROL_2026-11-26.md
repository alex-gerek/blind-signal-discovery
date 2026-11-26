# GW150914: current-source real-data control (2026-11-26)

This is a real O1 measurement control, not a synthetic injection or blind discovery. Current source3497211, unchanged detector/preprocessing. [Protocol](gw150914_current_control_protocol_2026-11-26.md), [complete report](../results/gw150914_current_control/20261126T170822Z_c5b3ad43/report.md), [independent checks](../results/gw150914_current_control/20261126T170822Z_c5b3ad43/independent_verification.json).

## Result

| Site | Frozen score / threshold | Frozen event | Adaptive Workstation event | Separate background nominations |
| --- | --- | --- | --- | --- |
| H1 |35.7673 /2.80539 |Detected |Detected, available |4/320=1.25% |
| L1 |8.52803 /3.80464 |Detected |Detected, available |3/320=0.9375% |

Both identify window39, GPS1126259460..1126259464, containing the known reference1126259462.422. Each has three candidate clusters near the reference time/band. No added waveform, selected seed or post-result threshold adjustment.

The target has0/71 off-source nominations at each site after the fixed padded +/-16s exclusion. Neither this nor the background fractions establish zero population FAR or new-event significance. H1/L1 counts are not pooled.

The adaptive target run established no new regime and did not update calibration: it agrees with the frozen detector at the known event. This is NOT a new real-data adaptive-recovery test or a resolution of L09/L10.

## Timing/frequency and diagnostics

- Strongest H1 candidate: medium125ms scale, GPSpeak1126259462.375, window bounds2.140625..2.453125s, frequency edges20..284Hz.
- Strongest L1 candidate: short31ms scale, GPSpeak1126259462.40625, window bounds2.345703125..2.427734375s, frequency edges16..304Hz.
- These different-scale STFT peaks are not a physical propagation-delay estimate. No coincidence significance or waveform-delay measurement was performed.
- Smooth30–380Hz filtering is not a hard bin mask. Candidate edges outside the passband do not prove physical power at those frequencies; per-frequency normalization and bin edges must be considered.
- L1 event window has a GLITCH_STORM diagnostic while remaining available/detected; clean duty79/80, detection availability80/80. A warning alone does not justify discarding a candidate or calling it instrumental. This observation does not isolate the warning's cause.
- Event padded-span DATA/DQ mask is127 at both sites. H1 injection mask31 passes all five NO_* flags; L1 mask23 has NO_CW_HW_INJ not passing. The absence-of-continuous-wave-hardware-injection flag is therefore not asserted for L1. Other injection flags pass. This is retained provider metadata, not a newly added waveform or proof that the detected transient is a hardware injection.

## Preservation and view

Two existing original GWOSC files and their official provenance hashes were checked; all four block footprints are disjoint including padding. Preprocessing PSD uses only its dedicated block; no target fitting. All prescribed windows were finite and DATA-qualified. Shared within-block padding and nonstationary real noise limit independence/FAP interpretations.

Frozen threshold/empirical-p/GPS decisions, actual Workstation CSV/JSON and derived CSV/NPY arrays were independently checked. The rank floor is1/225, not an astrophysical probability. Both internal calibration labels are VALIDATED; this is internal screening only.

[Full target timeline](../results/gw150914_current_control/20261126T170822Z_c5b3ad43/target_timeline.png).
[Readable event figure](../results/gw150914_current_control/20261126T170822Z_c5b3ad43/event_waveform_spectrogram_contrast.png).
Original full-color-range figure is retained alongside a second display-only contrast figure, whose pooled visible-bin percentile limits are disclosed in display_contrast.json; detector outputs were not changed.

Large derived calibration/target GPS CSV and exact NPY files are LOCAL ONLY under data/derived/gw150914_workstation/20261126T170822Z_c5b3ad43/. Hashes and source URLs are saved in each site's derived_file_manifest.json. Data are already whitened/bandpassed; do not preprocess twice. The current live UI file-ingestion path normally reads the reserved-event manifest even for O1 files; that path was not invoked or bypassed. A separate guard review is needed before owner-operated loading under the strict no-O3-manifest-read constraint. Synthetic UI tests do not establish this workflow.

Attempt1,20261126T170720Z_7311dda3, failed in the new driver before target analysis due to the wrong calibration serialization API. Its snapshots/partial record remain. The successful rerun corrected API calls only and added error recording.

## Interpretation and next step

The fixed current configuration recovered a known genuine real-data control in BOTH sites. This establishes recovery on these selected previously used files, not universal sensitivity, independent event confirmation or an astrophysical discovery. Inspect real background candidates next using their saved GPS/quality annotations; do not classify them automatically as noise or new signals. Adaptive calibration admission/Union drift problems remain open.

Source attribution: measured strain from GWOSC/LIGO Scientific Collaboration/Virgo Collaboration; original data CC BY4.0. [Official event data](https://gwosc.org/events/GW150914/), [GWOSC acknowledgement/licensing](https://gwosc.org/acknowledgement/), [O1/O2 open-data reference](https://doi.org/10.1016/j.softx.2021.100658). Project MIT licensing does not replace third-party data terms. No O3 manifest/input or non-allowlisted original data was accessed.
