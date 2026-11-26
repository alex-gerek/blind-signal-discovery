# Historical adaptive Phase 2b — synthetic replay

Run: `20261126T163629Z_75c7a756`. Current source `043f0eb`; historical recipe `edd1ec9`.
Twelve original synthetic cases completed without detector tuning. Source and historical artifact hashes match. All twelve current input arrays and the reference regenerated exactly. Original historical array hashes are unavailable; bitwise identity with old execution is not established.

## Channel-specific current measurements

| Case | Clean / total | STFT noise FP / available | STFT abstained noise | STFT detected / available injections | Full Union noise FP / available |
| --- | --- | --- | --- | --- | --- |
| stationary_null_200win | 200/200 | 1/200 | 0 | 0/0 (abstained 0) | 2/200 |
| abrupt_step_500win | 368/500 | 7/415 | 72 | 6/13 (abstained 0) | 186/286 |
| transient_chirps_stationary | 120/120 | 2/100 | 0 | 20/20 (abstained 0) | 0/100 |
| glitch_storm_120win | 60/120 | 57/120 | 0 | 0/0 (abstained 0) | 60/120 |
| stationary__80s_sinusoid | 120/120 | 1/100 | 0 | 20/20 (abstained 0) | 2/100 |
| stationary__80s_multitone_comb | 120/120 | 1/100 | 0 | 20/20 (abstained 0) | 2/100 |
| stationary__40s_slow_chirp | 120/120 | 1/110 | 0 | 9/10 (abstained 0) | 1/110 |
| stationary__80s_broadband_stochastic | 103/120 | 0/100 | 0 | 0/20 (abstained 0) | 0/100 |
| drifting_1.0_to_1.4x__80s_sinusoid | 25/120 | 2/100 | 0 | 20/20 (abstained 0) | 69/100 |
| drifting_1.0_to_1.4x__80s_multitone_comb | 28/120 | 4/100 | 0 | 19/20 (abstained 0) | 75/100 |
| drifting_1.0_to_1.4x__40s_slow_chirp | 24/120 | 1/110 | 0 | 9/10 (abstained 0) | 80/110 |
| drifting_1.0_to_1.4x__80s_broadband_stochastic | 29/120 | 0/100 | 0 | 0/20 (abstained 0) | 73/100 |

## Abrupt step: timings and fixed post-step injections

Timings below are zero-based window-start offsets from the shift at index100; one window is4 seconds. Credibility labels and actual available thresholds are different definitions.

{
  "first_new_regime_calibrated": 134,
  "first_new_regime_stft_available": 134,
  "first_new_regime_full_union_available": 267,
  "first_mature_available": 267
}

| Injection index | Available | Clean | Mature | Detected | Score | Used threshold |
| --- | --- | --- | --- | --- | --- | --- |
| 330 | True | True | False | True | 3.4772675 | 3.1360766413763086 |
| 340 | True | True | False | True | 51.49249 | 3.1360766413763086 |
| 350 | True | True | False | True | 30.891303 | 3.1360766413763086 |
| 360 | True | True | False | False | 0.74424472 | 3.1360766413763086 |
| 370 | True | True | True | False | 9.3242841 | 163.83940094146237 |
| 380 | True | True | True | False | 4.1635108 | 163.83940094146237 |
| 390 | True | True | True | False | 25.707155 | 163.83940094146237 |
| 400 | True | True | True | False | 92.839684 | 163.83940094146237 |
| 410 | True | True | True | False | 1.0152212 | 163.83940094146237 |
| 420 | True | True | True | False | 3.5137935 | 163.83940094146237 |

## Historical comparison

Only nonidentical historical fields follow. These compare the old aggregation recipe on current outputs; corrected available/clean/mature denominators are in independent_analysis.json.

[
  {
    "suite": "suite_a",
    "scenario": "abrupt_step_500win",
    "differences": {
      "duty_cycle_clean": {
        "historical": 0.88,
        "current_legacy_aggregation": 0.736
      },
      "cal_n_thresh": {
        "historical": 150,
        "current_legacy_aggregation": 199
      },
      "cal_n_eff": {
        "historical": 115.00275247453128,
        "current_legacy_aggregation": 181.77933125734825
      },
      "cal_ljung_box_p": {
        "historical": 0.43295466795377036,
        "current_legacy_aggregation": 1.0113175408221764e-08
      },
      "cal_beta_ci90_neff_hi": {
        "historical": 0.025758317220151185,
        "current_legacy_aggregation": 0.02386579523757849
      },
      "t_provisional_win": {
        "historical": 62,
        "current_legacy_aggregation": 134
      },
      "t_mature_win": {
        "historical": 201,
        "current_legacy_aggregation": 267
      },
      "transient_eff_full": {
        "historical": 0.5,
        "current_legacy_aggregation": 0.3
      },
      "transient_eff_clean_cal": {
        "historical": 0.5,
        "current_legacy_aggregation": 0.3
      },
      "glitch_ref_refreshes": {
        "historical": 2,
        "current_legacy_aggregation": 3
      }
    }
  }
]

## Interpretation limits

- This replay is a previously inspected development recipe with one seed per suite, not independent confirmation or a blind search.
- No new physical signal or astrophysical event can be established from synthetic inputs.
- Changed recovery timing, availability and thresholds are software/method observations. They do not establish methodological novelty.
- Zero observed mature false alarms is conditional on the recorded available denominator, not proof of zero population FAR.
- Glitch-storm and drifting-background counts are robustness diagnostics; availability alone does not validate exchangeability or stationary FAP.
- Persistent broadband versus background drift remains ambiguous. Association of missed injections with resonant noise is a hypothesis; per-injection frequency metadata are preserved in independent_analysis.json, not a causal intervention.
- Old A/B results remain preserved. Real-data suiteD, coherence suiteC and the earlier phase2/initial phase2b implementations were not rerun.
- No O3 input/manifest, commit, push, release or Zenodo operation was involved.
