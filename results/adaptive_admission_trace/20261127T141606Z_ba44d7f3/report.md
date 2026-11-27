# Adaptive calibration tracing: unchanged-baseline reproduction

Run `20261127T141606Z_ba44d7f3`, source `d049e2b`. Measurement-only experimental wrapper; production detector/manager unchanged. This is a deterministic historical reproduction and oracle comparison, not independent scientific confirmation or a deployable admission rule.

## Results

Both input-array hashes match their previous follow-up arms and all 1,000 serialized window records match exactly. At index367 (zero-based), the baseline threshold is163.83940094146237 and the oracle threshold1.6580395349128854. Later indices370/380/390/400/410/420 are available/mature in both arms and have identical scores; recovery is0/6 versus5/6. The oracle removes only four previously injected signals with known pre-injection noise; it is not available in real operation.

Used STFT calibration changes at indices0,162,234,367. The initial threshold is2.8957490413140143, index162 starts unavailable calibration (null), index234 restores a finite threshold3.1360766413763086. Later threshold inflation therefore occurs AFTER finite calibration availability has already recovered. Each new calibration affects only the following window, never the fitting window.

Baseline available noise exposure:415windows/1660s, with7post-veto STFT alarms; unavailable noise72windows/288s. Oracle:419/1676s,7alarms,72/288s unavailable; the four replaced injections now count as noise, so denominators differ by design. Never pool these paired arms as independent trials. Mature clean new-regime noise:127windows,0STFT but127broadband and127fullUnion nominations in both arms. L09/L10 remain open.

## Instrumentation and checks

- Each calibration has an exact ordered member list, signal hashes, actual normalization/threshold partition counts, calling window and earliest permissible use, and full saved normalization/null-score state. Known-injection membership is diagnostic truth only; never used by the manager to exclude observations.
- Direct per-window detector calls carry an unambiguous calibration-call ID, score, used threshold and empirical p or explicit unavailable reason. Saved windows retain full raw candidates and actual Union thresholds/statuses. Per-window records describe the used state; full trace retains separate dedicated/Union calls, not a guessed identity inferred from threshold equality.
- Independent verifier uses saved data only, no production detector imports. Checks16calibrations,3108member occurrences,1000windows and856available empirical p-values by direct null-tail counts. It checks exact archived outputs/inputs, partition labels, strict threshold decisions, causal switches, source hashes, and7interval fixtures (empty, zero-length, duplicate/overlap/gap, reversed, infinity, NaN, conflicting ID).
- There are zero exported raw candidates in unavailable windows in THIS recipe. This equality check cannot establish diagnostic candidate preservation under a nonempty unavailable-candidate case; that fixture is still needed before changing the policy or API. Seven interval fixtures validate the research exposure helper, not production coincidence or UI behavior.
- This is not a pytest suite result, general false-alarm validation, real-data experiment, admission-policy comparison, or independent validation of statistical exchangeability. No new confidence or superiority claims.

No measurement errors, interrupted runs or parameter changes. Before verification, a redundant unused verifier expression was removed; no detector output or measurement changed. No source/app/old result/test edits, O3/data/manifest/network reads or additional real event. Documentation/source drafts unrelated to this stage remain separate.

## Next decision

Trace-only instrumentation reproduces the baseline and is ready to measure a future admission policy. Before implementing one, freeze its observable-only rule, calibration-source separation, quarantine latency/resource behavior, exact synthetic scenarios and independent sequence-level inference. Do not use the four known injection labels as an admission rule. Keep the unchanged baseline and all historical failures. Do not mark L09/L10 repaired.
