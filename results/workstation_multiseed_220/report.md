# WHISKER / BSD Multi-Seed Calibration Study (Variant A: 220-Window Condition)

## 1. Study Scope and Exploratory Framing

This report documents an **initial exploratory multi-seed study** (`K = 20` independent calibration seeds) examining how much the false-alarm rate (FAR) varies across independent 220-window calibrations of the frozen raw `MultiResolutionDetector` (`target_fap = 0.01`, `norm_fraction = 0.30`).

- **Exploratory precision notice:** With `K = 20` calibration seeds, estimates of the between-calibration spread and tail exceedance probabilities carry substantial sampling uncertainty. This study is designed to characterize the order of magnitude of calibration-to-calibration variation and verify alignment with order-statistic theory, not to claim high-precision spread or tail estimates.
- **Observed vs. conditional true FAR:** For each seed $m$, $\hat{p}_m = k_m / 1000$ is the **observed per-seed FAR** on that seed's 1,000-window evaluation set. It is a noisy binomial realization around that frozen calibration's **unknown conditional true FAR** $p_m$.
- **Scope restriction:** Only the stateless **Frozen Raw Detector** (`MultiResolutionDetector.analyze`) is evaluated on synthetic `pink_colored` noise. Adaptive workstation monitoring, vetoes, real sensor data, and the reserved O3 holdout are not part of this experiment. Detector mathematics remain completely unmodified.

## 2. Experimental Protocol and Disjoint RNG Architecture

- **Git commit:** `eb21497`
- **Calibration seeds:** `K = 20` (`m = 1..20`), each calibrated on `N_bg = 220` independent 4-second `pink_colored` windows (`4096` samples at `1024 Hz`, `sigma = 1.0`).
- **Split structure:** Split A = `66` windows (30% per-scale median and `Q_0.95 - median` spread), Split B = `154` windows (70% global-maximum threshold selection).
- **Evaluation design (Variant A):** Each seed $m$ is evaluated on its own **separate independent evaluation set** of `N_eval = 1000` non-overlapping 4-second `pink_colored` windows. Separate evaluation sets ensure that evaluation noise is independent across seeds, allowing direct variance-component decomposition.
- **Leakage prevention:** Calibration streams use `make_split_rng(10000 + m, f"multiseed_cal220_seed_{m:02d}")` and evaluation streams use `make_split_rng(20000 + m, f"multiseed_eval1000_seed_{m:02d}")`. Evaluation data are never used to select, filter, or tune thresholds. Complete seed provenance is stored in `seed_mapping.csv`.
- **Generator independence:** Each call to `generate_pink_colored` draws a fresh spectral slope `alpha ~ Uniform(0.85, 1.15)`, fresh white Gaussian samples, per-window `rfft`/`irfft` shaping, and per-window RMS normalization (`sigma = 1.0`) with no carried filter state. `MultiResolutionDetector.analyze` is stateless across windows. Therefore, conditional on frozen calibration $m$, the 1,000 binary alarm indicators are modeled as i.i.d. $\mathrm{Bernoulli}(p_m)$ trials by construction of the generator and detector.
- **Total wall-clock runtime:** `574.7 s` (`9.58 min`).

## 3. Exact Threshold-Selection Rule and Order-Statistic Reference

Direct inspection of `MultiResolutionDetector.calibrate` (`src/bsd/detector_v2.py`) confirms the exact threshold rule for `N_bg = 220`, `norm_fraction = 0.30`, and `target_fap = 0.01`:

- `n_norm = max(20, int(round(220 * 0.30))) = 66`
- `n_thresh = 220 - 66 = 154`
- `k_order = min(154, ceil((1 - 0.01) * (154 + 1))) = min(154, ceil(153.45)) = 154`
- `global_threshold = sorted_global[153]` (the exact sample maximum $S_{(154)}$ of Split B)
- Decision rule in `MultiResolutionDetector.analyze`: strict inequality `std_score > global_threshold`.

### Exchangeability Verification (Split B vs. Evaluation Conditional on Split A)

The $\mathrm{Beta}(1, 154)$ reference requires that Split-B global scores and evaluation global scores be exchangeable conditional on the frozen Split-A standardization. Code inspection verifies all four prerequisites:

1. **Identical noise distribution:** Both Split B (`cal_obs[66:]`) and evaluation (`eval_obs`) are produced by identical calls to `generate_noise('pink_colored', 4096, 1024.0, rng, sigma=1.0)` -> `generate_pink_colored`, drawing fresh `alpha ~ Uniform(0.85, 1.15)`, fresh white Gaussian noise, per-window `rfft`/`irfft` shaping, and per-window RMS normalization (`sigma = 1.0`) with zero carried filter state.
2. **Identical stateless preprocessing:** `TemplateFreeDetector.analyze` (`src/bsd/detector.py`) applies per-window median subtraction, unpadded Hann STFT (`boundary=None, padded=False`), 2D median-filter PSD and envelope whitening, and 8-connected-component clustering strictly to `observation.signal`, storing no cross-window state and never conditioning on `sample_id`.
3. **Identical score calculation:** Both `MultiResolutionDetector.calibrate` (Split B) and `MultiResolutionDetector.analyze` (evaluation) call `_raw_analyze_all_scales(obs)`, standardize each scale's maximum cluster score with the frozen Split-A `median` and `spread_scale`, and take the maximum across the 3 scales. Decision-rule audit across all `20,000` evaluation windows confirmed `0` mismatches between `res.detected` and `res.max_global_score > frozen_cal.global_threshold`.
4. **Independent RNG stream treatment:** `make_split_rng(master_seed, split_name)` (`src/bsd/generator.py`) seeds NumPy's PCG64 generator via `SeedSequence([int(master_seed), int(code)])`. All 40 calibration and evaluation streams (`cal_master_seed = 10001..10020`, `eval_master_seed = 20001..20020`) use pairwise-distinct `(master_seed, code)` tuples. Within each calibration seed, Split A (`0..65`) and Split B (`66..219`) are drawn sequentially from the same PCG64 stream, so Split B is independent of Split A and exchangeable with the evaluation stream under PCG64.

### Theoretical Order-Statistic Derivation and Assumptions

Conditional on the frozen Split-A normalization parameters $\hat{\theta}_A$, the $N_B = 154$ Split-B global scores are independent draws from a common score distribution $F_{\hat{\theta}_A}$. If $F_{\hat{\theta}_A}$ is continuous in its upper tail (no ties at the maximum), then the true exceedance probability of the sample maximum $S_{(154)}$ for a fresh independent draw from the same distribution follows the exact distribution-free order-statistic law:

$$p_m = 1 - F_{\hat{\theta}_A}(S_{(154)}) \sim \mathrm{Beta}(N_B - k_{\text{order}} + 1,\; k_{\text{order}}) = \mathrm{Beta}(1, 154).$$

Because $\mathrm{Beta}(1, 154)$ does not depend on $\hat{\theta}_A$, integrating over the random Split-A realizations leaves the **marginal** distribution of $p_m$ across independent seeds equal to $\mathrm{Beta}(1, 154)$, and the observed evaluation count $k_m$ follows the compound $\mathrm{BetaBinomial}(1000, 1, 154)$ distribution.

- **Continuity / tie check:** Although `MultiResolutionDetector` assigns a raw cluster score of `0.0` when no time-frequency cluster passes the initial bin threshold (a lower-tail point mass), this does not affect the upper-tail order statistic as long as $S_{(154)}$ lies strictly above the zero-cluster floor with no ties. Across all `20` seeds, the number of seeds with a tied Split-B maximum was `0/20`, total exact duplicate score pairs inside Split B was `0`, and the minimum Split-B maximum threshold observed was `1.7060` (well above the zero-cluster floor).
- **Beta(1, 154) reference moments:**
  - Expected conditional true FAR $\mathbb{E}[p_m] = 1/155 = 0.6452\%$ (strictly below the 1.00% target FAP because `k_order = 154` is conservative).
  - Theoretical between-calibration SD $\mathrm{SD}(p_m) = 0.6410\%$ (median `0.4491%`, 5th-95th percentile `[0.0333%, 1.9265%]`).
  - Theoretical probability that a single frozen 220-window calibration has conditional true FAR $> 1.00\%$: `21.27%`.
  - Theoretical total SD of observed $\hat{p}_m$ at `N_eval = 1,000`: `0.6889%` (theoretical calibration share of total variance: `86.6%`).

## 4. Complete Calibration Status Breakdown (All 20 Seeds Included)

No seeds were excluded, re-run, or filtered:

- **Completed runs:** `20 / 20` (`0` failed runs)
- **`VALIDATED` calibrations:** `12 / 20`
- **`PROVISIONAL` calibrations:** `8 / 20`
- **`INSUFFICIENT` calibrations:** `0 / 20`

Seeds not achieving `VALIDATED` status (included in all primary aggregates):

- **Seed 01** (`PROVISIONAL`): unmet criteria = `n_eff=110.3 < 120; dependence-aware upper 90% FAR CI=2.68% > 2.5%`; `n_eff = 110.3`, `Ljung-Box p = 0.5405`, `observed FAR = 9/1000 (0.90%)`.
- **Seed 03** (`PROVISIONAL`): unmet criteria = `n_eff=106.7 < 120; dependence-aware upper 90% FAR CI=2.77% > 2.5%`; `n_eff = 106.7`, `Ljung-Box p = 0.6150`, `observed FAR = 29/1000 (2.90%)`.
- **Seed 08** (`PROVISIONAL`): unmet criteria = `n_eff=100.1 < 120; temporal dependence detected (Ljung-Box p=0.0013 < 0.05); dependence-aware upper 90% FAR CI=2.95% > 2.5%`; `n_eff = 100.1`, `Ljung-Box p = 0.0013`, `observed FAR = 8/1000 (0.80%)`.
- **Seed 10** (`PROVISIONAL`): unmet criteria = `temporal dependence detected (Ljung-Box p=0.0107 < 0.05)`; `n_eff = 122.7`, `Ljung-Box p = 0.0107`, `observed FAR = 14/1000 (1.40%)`.
- **Seed 12** (`PROVISIONAL`): unmet criteria = `n_eff=107.7 < 120; dependence-aware upper 90% FAR CI=2.74% > 2.5%`; `n_eff = 107.7`, `Ljung-Box p = 0.4819`, `observed FAR = 7/1000 (0.70%)`.
- **Seed 15** (`PROVISIONAL`): unmet criteria = `n_eff=106.1 < 120; dependence-aware upper 90% FAR CI=2.90% > 2.5%`; `n_eff = 106.1`, `Ljung-Box p = 0.1146`, `observed FAR = 10/1000 (1.00%)`.
- **Seed 16** (`PROVISIONAL`): unmet criteria = `n_eff=83.8 < 120; temporal dependence detected (Ljung-Box p=0.0309 < 0.05); dependence-aware upper 90% FAR CI=3.51% > 2.5%`; `n_eff = 83.8`, `Ljung-Box p = 0.0309`, `observed FAR = 26/1000 (2.60%)`.
- **Seed 17** (`PROVISIONAL`): unmet criteria = `n_eff=110.5 < 120; dependence-aware upper 90% FAR CI=2.67% > 2.5%`; `n_eff = 110.5`, `Ljung-Box p = 0.3523`, `observed FAR = 1/1000 (0.10%)`.

## 5. Per-Seed Results Table (All 20 Seeds)

| Seed | Cal Status | `n_eff` | Cal LB `p` | Frozen Threshold | False Alarms (`k/1000`) | Observed FAR | 95% Wilson CI | Eval LB `p` |
|---:|---|---:|---:|---:|---:|---:|---|---:|
| 01 | PROVISIONAL | 110.3 | 0.5405 | 2.7006 | 9/1000 | 0.90% | [0.47%, 1.70%] | 0.6496 |
| 02 | VALIDATED | 137.1 | 0.3798 | 4.1396 | 0/1000 | 0.00% | [0.00%, 0.38%] | 0.7119 |
| 03 | PROVISIONAL | 106.7 | 0.6150 | 1.9412 | 29/1000 | 2.90% | [2.03%, 4.13%] | 0.6161 |
| 04 | VALIDATED | 139.8 | 0.2795 | 2.0950 | 15/1000 | 1.50% | [0.91%, 2.46%] | 0.1392 |
| 05 | VALIDATED | 130.4 | 0.3532 | 2.8500 | 9/1000 | 0.90% | [0.47%, 1.70%] | 0.0974 |
| 06 | VALIDATED | 126.0 | 0.9749 | 2.8387 | 8/1000 | 0.80% | [0.41%, 1.57%] | 0.1002 |
| 07 | VALIDATED | 127.3 | 0.4174 | 4.9347 | 2/1000 | 0.20% | [0.05%, 0.73%] | 0.9016 |
| 08 | PROVISIONAL | 100.1 | 0.0013 | 1.9556 | 8/1000 | 0.80% | [0.41%, 1.57%] | 0.9989 |
| 09 | VALIDATED | 144.2 | 0.8749 | 3.1525 | 9/1000 | 0.90% | [0.47%, 1.70%] | 0.3581 |
| 10 | PROVISIONAL | 122.7 | 0.0107 | 2.6955 | 14/1000 | 1.40% | [0.84%, 2.34%] | 0.2850 |
| 11 | VALIDATED | 124.0 | 0.9188 | 3.2603 | 0/1000 | 0.00% | [0.00%, 0.38%] | 0.6708 |
| 12 | PROVISIONAL | 107.7 | 0.4819 | 2.6223 | 7/1000 | 0.70% | [0.34%, 1.44%] | 0.0247 |
| 13 | VALIDATED | 134.7 | 0.2323 | 2.8257 | 12/1000 | 1.20% | [0.69%, 2.09%] | 0.9127 |
| 14 | VALIDATED | 122.5 | 0.7097 | 2.3900 | 15/1000 | 1.50% | [0.91%, 2.46%] | 0.6489 |
| 15 | PROVISIONAL | 106.1 | 0.1146 | 1.8866 | 10/1000 | 1.00% | [0.54%, 1.83%] | 0.2358 |
| 16 | PROVISIONAL | 83.8 | 0.0309 | 1.7060 | 26/1000 | 2.60% | [1.78%, 3.78%] | 0.6833 |
| 17 | PROVISIONAL | 110.5 | 0.3523 | 4.4257 | 1/1000 | 0.10% | [0.02%, 0.56%] | 0.0025 |
| 18 | VALIDATED | 133.2 | 0.3841 | 2.4349 | 8/1000 | 0.80% | [0.41%, 1.57%] | 0.2005 |
| 19 | VALIDATED | 136.9 | 0.6845 | 2.2193 | 10/1000 | 1.00% | [0.54%, 1.83%] | 0.7234 |
| 20 | VALIDATED | 140.2 | 0.8792 | 5.1275 | 4/1000 | 0.40% | [0.16%, 1.02%] | 0.8560 |

> **Note on Ljung-Box `p`:** Ljung-Box tests for linear autocorrelation across lags 1..10. It is reported strictly as a diagnostic check and does **not** prove statistical independence; the i.i.d. Bernoulli evaluation model relies on the stateless generator and detector architecture.

## 6. Aggregate Summary and Variance Decomposition

### Observed FAR Distribution Across All 20 Seeds

- **Mean observed FAR ($\bar{p}$):** `0.980%` (95% bootstrap CI across $K=20$ seeds: `[0.670%, 1.325%]`; Beta(1, 154) theoretical mean: `0.645%`).
- **Median observed FAR:** `0.900%` (IQR: `0.625%`, range: `0.00%` to `2.90%`).
- **Seeds with observed FAR $> 1.00\%$:** `6 / 20` (`30.0%`).
- **Seeds with `0 / 1000` observed false alarms:** `2 / 20`.
- **Frozen threshold range across 20 seeds:** `1.7060` to `5.1275` (mean `2.9101`, SD `1.0047`).

### Variance-Component Decomposition (Unadjusted vs. Adjusted)

Under independent evaluation sets of size $N_{\text{eval}} = 1000$ per seed, the law of total variance gives:

$$\mathrm{Var}(\hat{p}_m) = \mathrm{Var}(p_m) + \mathbb{E}\!\left[\frac{p_m(1 - p_m)}{N_{\text{eval}}}\right].$$

Because $\mathbb{E}[\hat{p}_m(1 - \hat{p}_m) \mid p_m] = \frac{N_{\text{eval}} - 1}{N_{\text{eval}}} p_m(1 - p_m)$, the finite-sample-corrected unbiased moment estimator of the binomial evaluation noise variance is:

$$\widehat{V}_{\text{binom}} = \frac{1}{K}\sum_{m=1}^{K} \frac{\hat{p}_m(1 - \hat{p}_m)}{N_{\text{eval}} - 1}.$$

| Quantity | Value (variance scale) | Value (SD, percentage points) | 95% Bootstrap CI (`K = 20`) | Notes |
|---|---:|---:|---|---|
| Unadjusted total sample variance $s_{\hat{p}}^2$ / SD $s_{\hat{p}}$ | `5.743158e-05` | `0.758%` | `[0.414%, 0.997%]` | Includes both calibration spread and finite-$N_{\text{eval}}$ binomial noise |
| Finite-sample-corrected binomial noise $\widehat{V}_{\text{binom}}$ (`N_eval - 1 = 999`) | `9.659059e-06` | `0.311%` | - | Naive `1/N_eval` value: `9.649400e-06` (0.1% lower) |
| Adjusted calibration variance (unclipped) $\widehat{V}_{\text{cal,raw}} = s_{\hat{p}}^2 - \widehat{V}_{\text{binom}}$ | `4.777252e-05` | - | `[9.231e-06, 8.805e-05]` | Unbiased on variance scale before truncation; can be negative in small samples |
| Adjusted calibration variance (zero-truncated) / SD $\hat{\sigma}_{\text{cal}}$ | `4.777252e-05` | `0.691%` | `[0.304%, 0.938%]` | **Not unbiased** (truncation at 0 adds upward bias; square root adds downward bias via Jensen's inequality) |
| Theoretical $\mathrm{Beta}(1, 154)$ calibration variance / SD | `4.108968e-05` | `0.641%` | - | Exact order-statistic reference under continuous i.i.d. Split-B scores |

- **Bootstrap truncation rate:** `0.06%` of `5000` bootstrap resamples had $\widehat{V}_{\text{cal,raw}} \le 0$.
- **Estimated variance share:** Calibration-to-calibration variability accounts for roughly `83.2%` of the total variance of observed $\hat{p}_m$ (theoretical Beta-Binomial share: `86.6%`), while finite-sample binomial evaluation noise accounts for roughly `16.8%`.
- **Methodological limitations of the variance estimator:**
  1. Zero-truncation (`max(0, ...)`) makes $\widehat{V}_{\text{cal,trunc}}$ positively biased when true calibration variance is small relative to sampling noise.
  2. Taking the square root makes $\hat{\sigma}_{\text{cal}}$ negatively biased by Jensen's inequality ($\mathbb{E}[\sqrt{X}] \le \sqrt{\mathbb{E}[X]}$).
  3. Percentile bootstrap intervals over `K = 20` seeds are approximate and can under-cover in small samples, especially for variance components.

### Sensitivity Check: `VALIDATED`-Only Subset vs. All 20 Seeds

For transparency, restricting to the `12` seeds classified as `VALIDATED` gives mean observed FAR = `0.767%`, unadjusted SD = `0.521%`, and zero-truncated adjusted calibration SD = `0.442%`. All primary conclusions remain based on the complete unfiltered `K = 20` set.

### Exploratory Post-Hoc Comparison Against `BetaBinomial(1000, 1, 154)`

> **Post-hoc disclosure:** The Monte Carlo mean and median tail-probability checks below are **exploratory and post-hoc** — they were selected *after* inspecting the observed per-seed results and noting that the empirical ECDF lay to the right of the `BetaBinomial(1000, 1, 154)` curve. Because multiple summary statistics (mean, median, SD, and exceedance count) were examined on the same `K = 20` sample, nominal tail probabilities should not be treated as pre-registered confirmatory p-values.

Using `100,000` reproducible Monte Carlo draws of `K = 20` independent `BetaBinomial(n=1000, a=1, b=154)` samples (from the verified exact threshold rule `k_order = 154` of `154`; RNG `make_split_rng(88888, "multiseed_220_betabinom_mc_v1")`):

- **Observed mean (`0.980%`) vs. theory (`0.645%`):** upper-tail probability $P_{H_0}(\bar{p}_{\text{sim}} \ge \bar{p}_{\text{obs}}) = 0.0251$ (two-sided: `0.0501`).
- **Observed median (`0.900%`) vs. theory (`0.449%`):** upper-tail probability $P_{H_0}(\text{med}_{\text{sim}} \ge \text{med}_{\text{obs}}) = 0.0107$.
- **Observed unadjusted SD (`0.758%`) vs. theory (`0.689%`):** upper-tail probability $P_{H_0}(s_{\text{sim}} \ge s_{\text{obs}}) = 0.2730$.

**Interpretation guardrails:** A tail probability $P(\text{stat} \ge \text{obs} \mid H_0)$ measures how often a summary at least as large as the observed value occurs *if* the exact `BetaBinomial(1000, 1, 154)` model holds. It is **not** the probability that the observed shift is "just random" or that $H_0$ is true. Evaluating several summaries on the same small sample (`K = 20`) also inflates the chance that at least one summary exhibits a small nominal tail probability.

## 7. Figures

- `fig01_per_seed_far_forest.png`: Observed FAR and 95% Wilson CI for all 20 seeds, colored by calibration status (`VALIDATED` vs. `PROVISIONAL`), alongside the `1.00%` target FAP and `Beta(1, 154)` reference band.
- `fig02_distribution_and_variance_decomp.png`: Empirical CDF of observed per-seed FAR vs. the exact `Beta-Binomial(1000, 1, 154)` reference distribution (left) and variance-component SD decomposition with 95% bootstrap CIs (right).

## 8. Scientific Conclusions and Limitations

1. **Calibration-to-calibration spread is substantial at `N_bg = 220` (`N_thresh = 154`):** Because the 1% threshold is set by the exact sample maximum $S_{(154)}$ of 154 Split-B windows, individual frozen calibrations drawn from the exact same stationary `pink_colored` process exhibit wide variation in their realized FAR.
2. **Comparison with order-statistic theory:** The observed spread across the 20 seeds (unadjusted SD `0.758%`, adjusted zero-truncated calibration SD `0.691%`) aligns closely with the exact $\mathrm{Beta}(1, 154)$ / $\mathrm{BetaBinomial}(1000, 1, 154)$ reference (theoretical calibration SD `0.641%`, total SD `0.689%`). However, the observed center (mean `0.980%`, 95% bootstrap CI `[0.670%, 1.325%]`, median `0.900%`) is higher than the theoretical $\mathrm{Beta}(1, 154)$ mean (`0.645%`) and median (`0.449%`). In exploratory post-hoc Monte Carlo checks (`100,000` draws), such an upward shift has nominal upper-tail probabilities of `0.0251` (mean) and `0.0107` (median) under the exact reference model (not adjusted for examining multiple summaries, and not interpretable as the probability that the shift is due to chance). Because exchangeability was verified by code inspection, determining whether this shift persists or regresses toward `0.645%` requires a larger seed count `K`.
3. **Contextualizing the earlier single-seed study (`eb21497`):** The earlier `seed = 2026` 220-window calibration that yielded `31/2000` (`1.55%`) observed FAR lies inside the upper tail of the `Beta(1, 154)` distribution (roughly the 91st percentile of conditional true FAR), illustrating why single-seed negative controls must not be generalized as the population-average FAR of the calibration procedure.
4. **What this study does NOT establish:** It does not establish real-world (LIGO/GWOSC or other physical sensor) false-alarm control, does not evaluate the adaptive workstation layer, and does not provide high-precision tail estimates (which would require a larger number of seeds $K$ and/or larger $N_{\text{eval}}$).
