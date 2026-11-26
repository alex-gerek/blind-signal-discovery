# Confirmatory Protocol (Prespecified Before Execution)

- **Recorded UTC:** `2026-11-25T18:03:08.702799+00:00`
- **Git commit:** `8c27877`
- **Design:** `K = 100` independent 220-window calibrations (`Split A = 66`, `Split B = 154`), `N_eval = 500` independent 4-s `pink_colored` evaluation windows per calibration, Frozen Raw Detector (`MultiResolutionDetector.analyze`) only.
- **Separation from exploratory data:** The earlier 20-seed study (`commit 8c27877`, seeds `10001..10020` / `20001..20020`) is strictly exploratory and is excluded from all confirmatory calculations.
- **RNG streams:**
  - Calibration: `make_split_rng(30000 + m, f"confirm_cal220_seed_{m:03d}")` (`m = 1..100`)
  - Evaluation: `make_split_rng(40000 + m, f"confirm_eval500_seed_{m:03d}")` (`m = 1..100`)
  - Verified pairwise-distinct `(master_seed, code)` tuples: `200 / 200`.
- **Single primary statistic:** Mean observed FAR across all 100 seeds, $T = \bar{p} = \frac{1}{100}\sum_{m=1}^{100} \frac{k_m}{500} = S / 50000$, tested two-sided at $\alpha = 0.05$.
- **Prespecified Monte Carlo null simulation (estimates, not exact values):**
  - Model: `K=100 i.i.d. BetaBinomial(n=500, a=1, b=154)`
  - RNG: `make_split_rng(77777, "confirm_220_k100_n500_h0_v1")`, `M = 200,000` draws
  - Boundary/tie rule: `Inclusive conservative discrete rule on integer total false alarms S = sum(k_m): s_lo = max{s : P_hat(S <= s) <= 0.025}, s_hi = min{s : P_hat(S >= s) <= 0.025}. Reject H_0 iff S_obs <= s_lo or S_obs >= s_hi. No randomized tie-breaking.`
  - MC-estimated critical integers: `s_lo = 253` (`0.5060%`), `s_hi = 399` (`0.7980%`)
  - MC-estimated achieved Type-I error: `0.04821` (lower tail `0.02420` + upper tail `0.02401`, MC SE `+/- 0.00048`).
- **Power scope caveat:** Pre-study power estimates (`51.9%` at `0.80%`, `72.4%` at `0.85%`, `85.9%` at `0.90%`, `96.7%` at `0.98%`) apply specifically to the proposed `Beta(1, (1 - mu_1)/mu_1)` alternative model and stated effect sizes; they are not a general guarantee for arbitrary non-Beta deviations.
- **Technical-failure handling:** Every seed outcome (`VALIDATED`, `PROVISIONAL`, `INSUFFICIENT`, `FAILED`) is recorded without selective dropping or replacement. Infrastructure exceptions mark the run `INCOMPLETE_TECHNICAL_FAILURE`; detector exceptions on valid input are logged verbatim as `FAILED`.
