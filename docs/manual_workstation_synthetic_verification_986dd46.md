# Manual Workstation Verification on Synthetic Data (`development/v0.3`, Commit `986dd46`)

> **Verification Type:** Manual Workstation verification on synthetic data. Visual inspection of waveforms, spectrograms, candidate bounding boxes, and UI tables was performed by the project owner in the local BSD Workstation interface (`app/bsd_workstation.py`), and numerical outcomes were cross-checked against an in-memory computational control run using identical parameters on commit `986dd46`.

---

## 1. Provenance and Artifact Availability

- **Repository branch / revision:** `development/v0.3`, commit `986dd46`.
- **Saved export bundle (`timeline.csv`, `candidates.csv`, `run_report.json`):** **Not saved** from the manual UI session. No `results/workstation_runs/` export bundle exists for commit `986dd46` (the only existing directory under `results/workstation_runs/` is an unrelated earlier run from commit `40c7bca`).
- **Screenshots:** Visual inspection was performed directly in the local Workstation UI by the project owner; UI screenshots were not committed to the repository.
- **Provenance policy:** Because no UI export bundle was saved to disk, no synthetic `run_id`, SHA-256 file hashes, or UI execution timestamps are recorded or inferred. All reported metrics below reflect the manual UI observations cross-checked against a deterministic in-memory control run (`run_synthetic_self_test`) on commit `986dd46`.

---

## 2. Shared Configuration and Calibration Status

| Parameter | Value |
| --- | --- |
| Search mode | `transient` |
| Processing mode | `causal_trailing` |
| Window duration | `4.0 s` (`1024 Hz`) |
| Target FAP (`target_fap`) | `0.01` |
| Background calibration windows (`n_bg_windows`) | `220` (`N_thresh = 154` Split-B windows) |
| Negative-control noise-only windows (`n_neg_control_windows`) | `80` |
| Positive injections (`n_pos_injections`) | `10` (`100` target windows total) |
| Noise model | `pink_colored` |
| Injected signal family | `linear_chirp` |
| Self-test seed (`seed`) | `2026` |
| Calibration seed (`calibration_seed`) | `0` |
| Heuristic vetoes (`enable_glitch_veto`, `enable_line_veto`) | Off (`False`, `False`) |
| Active STFT threshold (`stft_threshold`) | `2.8889` (`2.888922` in control run) |

### Calibration and Union Availability

- **Calibration status:** **`PROVISIONAL`** (`n_thresh = 154`, `n_eff = 99.73 < 120`, Ljung-Box linear-lag `p = 0.5687`, dependence-aware upper 90% FAR CI `= 2.96% > 2.5%`).
- **Union mode status:** **`UNION_UNAVAILABLE`** across all `100 / 100` target windows (`union_stft_threshold = null`, warning `UNION_STFT_UNAVAILABLE`). At `target_fwer_union = 0.01`, the Sidak union threshold requires `alpha_union_sidak ≈ 0.00501` (`N_thresh >= 199`), which exceeds the `154` Split-B samples available from `220` background windows.

---

## 3. Summary of Synthetic Self-Test Outcomes (`SNR = 18` vs. `SNR = 8`)

Both SNR configurations and all outcome categories (true detections, missed injections, and noise-only false alarms) are recorded below:

| Target Root-Energy SNR | Injections Detected | Injections Missed | False Alarms on 80 Noise-Only Windows | Total Active Triggered Windows |
| ---: | ---: | ---: | ---: | ---: |
| `18` | `10 / 10` | `0 / 10` | `1 / 80` | `11 / 100` |
| `8` | `3 / 10` | `7 / 10` | `1 / 80` | `4 / 100` |

---

## 4. Manually Inspected Windows at `SNR = 8` (Project Owner Visual Review)

| Window Index | Generator Ground Truth | Active Decision | Observed Score vs. Threshold and Spectrogram Notes |
| ---: | --- | --- | --- |
| `4` | Noise-only | False alarm (`True`) | `stft_max_score = 3.8101` (`> 2.8889`); 1 STFT candidate at `0.69–1.44 s` within window (`t_abs = 16.6875–17.4375 s`), `1.0–13.0 Hz`. |
| `81` | Injection present | Missed (`False`) | `stft_max_score = 0.3680` (`< 2.8889`); injected `linear_chirp` at `0.87–2.47 s` within window (`t_abs = 324.873–326.473 s`), `99.2–247.8 Hz`, `0` nominated candidates. |
| `83` | Injection present | Detected (`True`) | `stft_max_score = 8.6062` (`> 2.8889`); injected `linear_chirp` at `1.98–2.96 s` within window (`t_abs = 333.983–334.958 s`), `285.2–385.3 Hz`; `3` STFT candidate bounding boxes visible on spectrogram (`8.6062`, `7.3816`, `7.0643`). |

---

## 5. Observed Conclusions and Scope Limitations

- **Observed sensitivity reduction:** Lowering the injected root-energy SNR from `18` to `8` on this fixed synthetic realization (`seed = 2026`, `pink_colored`, `linear_chirp`) reduced injection recovery from `10 / 10` to `3 / 10` (`7 / 10` missed).
- **Noise-only false alarms occur above threshold:** Window `4` demonstrates that a pure-noise realization can exceed the calibrated `PROVISIONAL` threshold (`3.8101 > 2.8889`) and produce a localized STFT candidate box.
- **Non-detection does not prove absence of signal:** Window `81` demonstrates that an injected `SNR = 8` `linear_chirp` can remain below the detection threshold (`0.3680 < 2.8889`) without producing any candidate.
- **Strict scope limitation:** These two single-seed synthetic runs (`100` target windows each) do **not** establish a general false-alarm rate (FAR), universal detection sensitivity across signal/noise families, or physical origin of any candidate. Historical benchmark reports under `results/` remain untouched.
