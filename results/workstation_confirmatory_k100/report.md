# WHISKER / BSD Confirmatory Multi-Seed Study (K = 100, N_bg = 220, N_eval = 500)

## 1. Prespecified Confirmatory Objective and Separation from Exploratory Data

This report documents a **separate confirmatory synthetic experiment** (`K = 100` independent 220-window calibrations, `N_eval = 500` independent evaluation windows per calibration, `50,000` total evaluation windows) evaluating the Frozen Raw `MultiResolutionDetector` on synthetic `pink_colored` noise.

- **Strict separation:** The earlier 20-seed study (`commit 8c27877`, seeds `10001..10020` / `20001..20020`) is strictly exploratory and was **not** reused in any confirmatory calculation.
- **Prespecified protocol:** Recorded in `protocol.json` and `protocol.md` prior to seed execution.
- **Unmodified boundaries:** Detector mathematics, GUI, historical baselines, and the reserved O3 holdout were completely untouched.
- **Total wall-clock runtime:** `1732.1 s` (`28.87 min`).

## 2. Prespecified Primary Test and Monte Carlo Null Specification

- **Verified exact threshold rule:** For `N_bg = 220`, `norm_fraction = 0.30`, `target_fap = 0.01`, `MultiResolutionDetector.calibrate` sets `n_norm = 66`, `n_thresh = 154`, and `k_order = min(154, ceil(0.99 * 155)) = 154` (exact Split-B sample maximum $S_{(154)}$), with strict inequality `std_score > global_threshold` in `analyze`. Under conditional exchangeability and upper-tail continuity, $p_m \sim \mathrm{Beta}(1, 154)$ ($\mu_0 = 1/155 = 0.64516\%$) and $k_m \sim \mathrm{BetaBinomial}(500, 1, 154)$.
- **Single primary statistic:** Mean observed FAR across all 100 calibrations, $T = \bar{p} = S / 50000$, where $S = \sum_{m=1}^{100} k_m$ is the integer total false-alarm count.
- **Two-sided discrete decision rule (Monte Carlo estimates, not exact analytical values):**
  - Null simulation: `M = 200,000` draws via `make_split_rng(77777, "confirm_220_k100_n500_h0_v1")`.
  - Inclusive conservative boundaries on integer $S$: `s_lo = 253` (`0.5060%`) and `s_hi = 399` (`0.7980%`). Reject $H_0$ iff $S_{\text{obs}} \le 253$ or $S_{\text{obs}} \ge 399$. No randomized tie-breaking.
  - MC-estimated achieved Type-I error: `0.04821` (lower tail `0.02420` + upper tail `0.02401`, MC SE `+/- 0.00048`).
- **Power scope caveat:** Pre-study power estimates (`51.9%` at `0.80%`, `72.4%` at `0.85%`, `85.9%` at `0.90%`, `96.7%` at `0.98%`) apply specifically to the proposed `Beta(1, (1 - mu_1)/mu_1)` alternative model and stated effect sizes; they are not a general guarantee for arbitrary deviations.

## 3. Primary Confirmatory Result

| Primary Quantity | Value |
|---|---|
| Total evaluation windows (`K x N_eval`) | `50,000` (`100 x 500`) |
| Observed total false alarms ($S_{\text{obs}}$) | `316` |
| Observed mean FAR ($\bar{p}$) | `0.6320%` |
| Reference $\mathrm{Beta}(1, 154)$ mean ($\mu_0 = 1/155$) | `0.6452%` |
| Observed difference ($\bar{p} - \mu_0$) | `-0.0132` percentage points |
| 95% nonparametric bootstrap CI for $\bar{p}$ (`10,000` resamples) | `[0.5180%, 0.7520%]` (diff: `[-0.1272, +0.1069]` pp) |
| 95% t-interval for $\bar{p}$ (`df = 99`) | `[0.5132%, 0.7508%]` (diff: `[-0.1320, +0.1057]` pp) |
| MC-estimated critical region ($\alpha = 0.05$, inclusive) | `S <= 253` (`<= 0.5060%`) or `S >= 399` (`>= 0.7980%`) |
| MC-estimated lower / upper tail probabilities | `P_hat(S <= S_obs) = 0.44642`, `P_hat(S >= S_obs) = 0.56436` |
| MC-estimated two-sided p-value | `0.89284` (MC SE `+/- 0.00069`) |
| **Prespecified primary decision** | **`FAIL_TO_REJECT_H0`** |

## 4. Complete Calibration Status Breakdown and Audit Checks (All 100 Seeds)

No seeds were dropped, replaced, or filtered:

- **Completed runs:** `100 / 100` (`0` failed runs)
- **`VALIDATED` calibrations:** `50 / 100` (`50.0%`)
- **`PROVISIONAL` calibrations:** `50 / 100` (`50.0%`)
- **`INSUFFICIENT` calibrations:** `0 / 100` (`0.0%`)
- **Split-B continuity / consistency audit:** `0/100` seeds with a tied Split-B maximum, `0` total Split-B duplicate scores, `0` mismatches between `res.detected` and `res.max_global_score > frozen_cal.global_threshold` across all `50,000` evaluation windows.
- **RNG isolation check:** `200/200` pairwise-distinct `(master_seed, code)` tuples across all calibration and evaluation streams.

Complete per-seed results (all 100 rows, including unmet validation criteria for every non-`VALIDATED` seed) are saved in `per_seed_results.csv` and `seed_mapping.csv`.

## 5. Secondary Descriptive Summaries (Not Part of Primary Decision Rule)

- **Median observed FAR:** `0.400%` (IQR `0.800%`, range `0.00%` to `2.60%`; `Beta(1, 154)` theoretical median: `0.449%`).
- **Seeds with observed FAR $> 1.00\%$:** `19 / 100` (`19.0%`).
- **Seeds with `0 / 500` observed false alarms:** `20 / 100`.
- **Frozen threshold range:** `1.8266` to `7.8462` (mean `3.3102`, SD `1.1128`).
- **Variance decomposition (all 100 seeds):**
  - Unadjusted total SD $s_{\hat{p}}$: `0.599%` (95% bootstrap CI `[0.496%, 0.692%]`; theoretical Beta-Binomial total SD: `0.734%`).
  - Finite-sample-corrected binomial noise SD (`N_eval - 1 = 499`): `0.354%` ($\widehat{V}_{\text{binom}} = 1.251415e-05$).
  - Adjusted zero-truncated calibration SD $\hat{\sigma}_{\text{cal}}$ (not unbiased due to zero-truncation and Jensen's inequality): `0.483%` (95% bootstrap CI `[0.364%, 0.583%]`; theoretical `Beta(1, 154)` SD: `0.641%`).
- **`VALIDATED`-only descriptive subset (`n = 50`):** mean observed FAR = `0.664%`, unadjusted SD = `0.645%`, adjusted zero-truncated calibration SD = `0.533%` (reported for transparency only; primary inference uses all 100 seeds).

## 6. Figures

- `fig01_confirmatory_primary_test.png`: Monte Carlo null distribution of $\bar{p}$ with inclusive $\alpha = 0.05$ rejection regions vs. observed mean (left), and empirical ECDF of all 100 per-seed FARs vs. `BetaBinomial(500, 1, 154)` CDF (right).
- `fig02_per_seed_distribution_and_status.png`: Per-seed observed FAR across all 100 seeds colored by calibration status (left), and frozen threshold $S_{(154)}$ vs. observed FAR (right).

## 7. Interpretation and Limitations

1. **Primary interpretation (`FAIL_TO_REJECT_H0`):** The confirmatory mean observed FAR (`0.6320%`, difference `-0.0132` pp from `0.6452%`, 95% bootstrap CI for difference `[-0.1272, +0.1069]` pp; MC-estimated two-sided `p = 0.89284`) does not fall in the prespecified $\alpha = 0.05$ rejection region.
2. **What this does NOT prove:** Failing to reject $H_0$ is **not** proof that $\mu = 1/155$ exactly or proof that any discrepancy has disappeared. The confidence interval quantifies the range of residual mean shifts compatible with these 100 confirmatory calibrations.
3. **Scope limitations:** Results apply strictly to synthetic stationary `pink_colored` noise evaluated with the stateless Frozen Raw `MultiResolutionDetector` at `N_bg = 220`, `target_fap = 0.01`, and 4-second windows. They do not establish false-alarm control on real LIGO/GWOSC or other physical sensor data, nor do they evaluate the adaptive workstation monitoring layer or vetoes.
