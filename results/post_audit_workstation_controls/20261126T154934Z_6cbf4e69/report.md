# Post-audit Python Workstation controls — bounded synthetic study

- Run: `20261126T154934Z_6cbf4e69`; source `043f0eb95be8eeae65b9a5f5bb63015d93559b73`.
- The detector implementation was unchanged. The new uncommitted driver, protocol and dirty-path list are recorded in manifest.json.
- Completed cases: 15/15. Accounting status: ACCOUNTING_VERIFIED.
- Inputs are generated synthetic arrays. O3-manifest reads and network access were blocked.

## Stationary controls

| Seed | SNR | Noise FP / available | Noise abstained | Injections detected / available | Injection abstained | Calibration |
| --- | --- | --- | --- | --- | --- | --- |
| 2026 | 8 | 3/200 | 0 | 4/10 | 0 | PROVISIONAL |
| 2026 | 18 | 3/200 | 0 | 10/10 | 0 | PROVISIONAL |
| 2027 | 8 | 5/200 | 0 | 5/10 | 0 | PROVISIONAL |
| 2027 | 18 | 5/200 | 0 | 10/10 | 0 | PROVISIONAL |
| 2028 | 8 | 2/200 | 0 | 4/10 | 0 | VALIDATED |
| 2028 | 18 | 2/200 | 0 | 10/10 | 0 | VALIDATED |

Paired SNR controls share identical background hashes and identical pre-injection exported noise scores/thresholds/decisions. Count them once: 10/600, available exposure 2400 s.
The fraction is descriptive across three frozen calibrations; a pooled common-rate Bernoulli confidence interval is not justified by this design.

## Adaptive and resource controls

| Case | First unavailable | First new-regime available | First mature available | Unavailable windows | Noise FP / available | Injections detected / available | Abstained injections |
| --- | --- | --- | --- | --- | --- | --- | --- |
| adaptive_seed601_fap0.01 | 72 | 144 | 277 | 72 | 6/288 | 10/10 | 0 |
| adaptive_seed601_fap0.005 | 72 | 277 | None | 205 | 1/155 | 10/10 | 0 |
| adaptive_seed611_fap0.01 | 72 | 144 | 277 | 72 | 0/288 | 10/10 | 0 |
| adaptive_seed611_fap0.005 | 72 | 277 | None | 205 | 0/155 | 10/10 | 0 |
| adaptive_seed621_fap0.01 | 72 | 144 | 277 | 72 | 0/288 | 10/10 | 0 |
| adaptive_seed621_fap0.005 | 72 | 277 | None | 205 | 0/155 | 10/10 | 0 |
| capped_seed601_fap0.005 | 72 | None | None | 108 | 0/72 | 0/0 | 0 |
| capped_seed611_fap0.005 | 72 | None | None | 108 | 0/72 | 0/0 | 0 |
| capped_seed621_fap0.005 | 72 | None | None | 108 | 0/72 | 0/0 | 0 |

Indices are zero-based from the amplitude step at target index 0; multiply by 4 seconds for target-relative window start time. Recovery is recorded when the new threshold is actually used, not when it is prepared for the next window.
The adaptive noise denominator includes available decisions before the monitor recognized the shift, as well as available post-recovery noise. It is not a stationary FAP estimate or a mature-only FAR.
Adaptive inputs and injections are paired across FAP settings: six scenario results are not 60 independent injected signals. At union FWER=0.005 the Sidak per-channel alpha is about 0.00250313, requiring at least 399 Split-B scores: 570 initial windows at the 30% normalization split, or 532 recalibration windows at the 25% split. The 320-reference/370-target design cannot establish full Union recovery; unavailable Union is expected, not evidence of a software defect.
Union status counts and resource flags are retained in independent_verification.json; a recovered dedicated STFT channel must not be described as full Union recovery.

## Findings and scope

[]

- This pilot has three independent calibration realizations, one stationary noise model and one amplitude-step model. It does not replace the 100-seed historical confirmatory study.
- Strong synthetic recovery does not establish sensitivity to weak signals, resonant shifted_ar2 noise, persistent broadband signals, real sensors or gravitational-wave events.
- All misses and unavailable cases are retained. No seed was replaced and no detector parameter was tuned.
- The 15 cases are experiments, not 15 additional pytest tests. Existing saved regression evidence remains 135 passed / 1 deselected at the previously recorded source.
- Headless export checks do not establish native Streamlit/browser interaction or decoder behavior. Manual operator verification remains pending.
- Historical adaptive/coincidence/GW170814 and browser stride-decimation outputs remain historical and require separate planned reruns.
- The initial system-Python import failure is recorded in docs/validation/post_audit_control_import_failure_20261126.json.
