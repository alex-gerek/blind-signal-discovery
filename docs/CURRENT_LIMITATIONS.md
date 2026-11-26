# Current limitations and independent review record

Reviewed on 2026-11-26: `development/v0.3` at `e1bd520`, same tree as public `main` at `7e72a12`. This record preserves the findings behind the documentation correction. The issues below are **open** unless explicitly stated otherwise. Documentation and version alignment do not fix them.

## Implementations and verification scope

The Python Local Measurement Workstation wraps the three-scale detector and adaptive managers. The browser demo in `docs/workstation/index.html` is a separate single-scale implementation with different sampling, window length, RNG and noise models. A shared version series or model label is not evidence of numerical equivalence.

The author recorded 87/87 Python tests passing at `986dd46`. An independent review at `e1bd520` ran 70 selected Python tests successfully, with two selected tests deselected, in 261.09 seconds. The selected files were `test_adaptive_regime.py`, `test_phase2b.py`, `test_detector_v2.py`, `test_coincidence.py`, `test_regime_monitor.py`, `test_vetoes.py`, `test_bsd.py`, and `test_codex_audit_fixes.py`. The selection excluded `result_serialization_and_end_to_end` and `workstation_export`; it did not run every repository test file. The harness disabled bytecode/cache/log writes, blocked filesystem mutation and networking, and blocked O3 manifest access. This was not a full 87-test rerun or an end-to-end GUI validation. Console evidence was retained in the review conversation; no raw log artifact was saved.

Browser checks executed the actual inline functions in Node.js with an in-memory DOM and stubbed canvas calls. WAV checks supplied a decoded synthetic array through a stub AudioContext; they establish the decimation error, not real-browser codec behavior. Python upload/aggregation checks used in-memory mocks. None of these checks involved O3 input data or new scientific performance studies.

## Confirmed open issues

| ID | Evidence at reviewed revision | Impact and proposed verification |
| --- | --- | --- |
| L01 | `readAudioOrCsv` uses stride decimation without an antialias filter. A decoded 48 kHz, 400 Hz synthetic tone produced a near-unit-amplitude alias at 110.6383 Hz after decimation to 510.6383 Hz. | Use disclosed antialiased resampling; test known in-band and out-of-band tones. Existing WAV results require reconsideration. |
| L02 | Browser CSV input `1`, `NaN`, `2`, `Infinity`, `3` becomes `[1,2,3]`; nonfinite rows are silently removed. | Time coordinates are compressed. Reject malformed/nonfinite data with row details and verify timestamp regularity. |
| L03 | Browser local-file workflow accepts the same CSV as background and target; a 20-window example completed without an independence warning. Loaded sampling rates are not compared. | Calibration/evaluation separation and frequency interpretation are unverified. Test identical inputs, overlapping intervals and mismatched rates; record unverified provenance explicitly. |
| L04 | With FAP=0.005 and 40 browser background windows, all 12 target windows were unavailable, yet summary displayed FAR=0/10 and 2/2 missed injections. Python `run_synthetic_self_test` also divides false alarms by all negative-control windows; an isolated unavailable-result check returned FAR=0 for 2 unavailable windows. | Report available exposure, conditional false-alarm fraction and abstentions separately. No available exposure must not imply verified zero false-alarm probability. |
| L05 | Streamlit `_save_upload_to_local_tmp` uses basename only. Target and background both named `recording.csv` resolve to the same path; the second save replaces the first before loading. | Save each role under a unique path; test different contents with the same basename. |
| L06 | Browser `Number(seed) || 2026` maps allowed seed 0 to 2026; both produced identical calibration and timeline in the check. Export metadata omits seed, noise model and source revision. | Validate and preserve seed 0; export complete configuration, implementation identity and input provenance. |
| L07 | Browser and Python timeline exports round scores and thresholds. The illustrative accepted crossing `2.88893 > 2.888922` becomes equal displayed/exported values `2.8889` and `2.8889` with decision `true`. | Preserve full-precision decision inputs alongside presentation values; test CSV/JSON round trips near exact boundaries. |
| L08 | v0.2 benchmark callers omit `active_intervals`, so FAR/significance are unavailable. The report table handles nulls, but the pre-correction narrative claimed significance resolution `1/501`. | The documentation-only patch corrects the narrative. Actual overlap exposure and historical scientific metrics still need a separate reviewed rerun. |

The browser currently stores only nominated candidates, including an empty candidate list during unavailable calibration; it does not implement Python's full diagnostic-candidate preservation. Its count-based `VALIDATED_SYNTH_COUNT` export label is not the Python `VALIDATED` screening procedure. Browser pink-noise shaping is pure 1/f, unlike the Python knee/random-slope generator. These differences need their own protocol, tests and clear scope; do not transfer Python reference values or test claims to the browser.

## Assumptions and work not verified

- Exchangeability underlies conformal guarantees. Linear autocorrelation screening, an estimated effective sample size, and a `VALIDATED` label do not prove it for arbitrary sensors. Dependence-adjusted Beta intervals are heuristic.
- Sidak channel budgeting requires channel-dependence assumptions. Time slides of STFT matches do not supply significance for waveform-correlation-confirmed events.
- The fallback named RNG stream mapping can collide: names `ab` and `ca` map to the same stream code for the same master seed. Repeatability alone does not establish independence or data separation.
- Historical adaptive recovery, sensitivity, exposure and time-slide FAR have not been recomputed with the audited code. Preserve original outputs and run any follow-up into a fresh location. See the [erratum](AUDIT_ERRATUM_47a5475.md).
- The previous manual Python verification saved no UI export bundle. It is an observation/control record, not a complete replay artifact.
- A fresh installation, complete interactive GUI behavior and general real-data performance were not verified in this review. No candidate is established as novel or astrophysical.

## Public repository security and redistribution snapshot

At the review time, secret scanning, push protection, Dependabot, Python CodeQL, full-SHA Action pinning, private vulnerability reporting, and main-branch force-push/deletion protection were enabled. No open alerts were returned. Ordinary main-branch pushes did not require a pull request or test status. Generic non-provider secret patterns remained disabled. Forty-five reachable commits, logs of 17 Actions runs and non-excluded text in three Pages artifacts had no matches under the review's secret patterns; this is a bounded scan, not a guarantee that every possible secret is absent.

The reserved O3 manifest was excluded from content inspection; only its filename was observed in Pages artifacts. Public-distribution intent for that manifest remains an owner decision. Excluding a file from Pages would not remove it from the public Git history. No O3 data were used.

MIT covers the project's licensed original material. The repository's general third-party-data disclaimer is not a complete attribution inventory: GWOSC data have CC BY 4.0 requirements, dependency licenses differ, and icon provenance is not documented. See [GWOSC acknowledgement requirements](https://gwosc.org/acknowledgement/). Resolve source attribution and the scope of redistributed third-party material before the next release.

## Correction policy

This documentation patch preserves historical result files, archived release `v0.3.0` at `c406e00`, and existing DOI identity. No Zenodo metadata were changed. Numerical/backend fixes, new browser reference studies, and archival metadata updates are separate steps requiring their own evidence. Keep failures and superseded results visible.
