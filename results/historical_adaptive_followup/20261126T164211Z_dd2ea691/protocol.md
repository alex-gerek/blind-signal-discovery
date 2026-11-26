# Exploratory adaptive calibration follow-up (2026-11-26)

Specified after the historical replay and before this follow-up. Parent run: 20261126T163629Z_75c7a756, current source043f0eb.

## Observation and question

The exact recipe replay has a mature STFT threshold jump from3.1360766414 to163.8394009415 at the decision after window366. All six fixed injections370..420 then fail despite available mature decisions. This is not adequately explained by the old aggregate statement about spectral resonance alone.

Hypothesis: known earlier injected signals entered the threshold-calibration partition and contributed to the elevated threshold. Inspect membership before asserting it. This is a post-hoc exploratory hypothesis, not prespecified confirmation or detector tuning.

## Fixed paired experiment

1. Deterministically regenerate original suiteA inputs; compare all current input hashes with the parent. Preserve raw shifted_ar2 noise before signal addition.
2. Run only the abrupt500-window case twice with identical current detector/config/reference: baseline, and a control replacing the four injected windows330/340/350/360 with their exact pre-injection noise. Keep transitional injections108/125/145, all later evaluated injections370..420 and every other sample unchanged. No other removals, seed changes or parameter trials.
3. Instrument MultiResolutionDetector.calibrate by a forwarding subclass. Save exact sample-ID membership, normalization/threshold partition, known injected members and returned thresholds for every calibration. Return the calibration object unchanged.
4. Baseline full serialized outputs must match the parent exactly. Assert all control inputs outside the four intervention windows are bitwise identical.
5. Measure detections/availability/thresholds on the same six later injected windows; inspect whether the four removed signals were actually members of threshold SplitB. Do not select different control windows after seeing results.

## Safety, preservation and interpretation

Fresh UTC/UUID results/historical_adaptive_followup directory, protocol/driver/source hashes, per-arm traces and full per-window output. Block external data/O3/network using the parent replay guard. No detector source change, tuning, historical overwrite, commit/push/release or Zenodo operation.

A paired effect would establish a mechanism in this single constructed sequence. It would not establish prevalence, universal failure, physical novelty, a ready fix, or the correctness of oracle removal in an unknown real stream. A separate design review is needed before changing calibration admission rules; excluding true noise extremes can invalidate false-alarm calibration.
