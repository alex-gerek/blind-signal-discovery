# Whole-block tail guard findings — 27 November 2026

[Prespecified pilot](block_tail_guard_protocol_2026-11-27.md) · [Full report](../results/block_tail_guard/20261127T152327Z_22827dc3/report.md) · [Independent verification](../results/block_tail_guard/20261127T152327Z_22827dc3/independent_verification.json).

**Decision: reject this rule; no production integration.** Whole-block quarantine when maximum exceeds median+10*max(Q95-median,1) was fixed before new synthetic evaluation. No individual tail trimming or truth used by the guard. Three seeds/four paired cases, two FAPs,265-window calibrations. Fixed-schedule one-block pilot, not adaptive recovery validation.

AtFAP0.01,all3pink-only blocks admitted; signal nominations22/36unchanged. All3chirp-contaminated blocks rejected,4baseline nominations lost to abstention. Burst-contaminated blocks2/3rejected,oneadmitted;18->6signal nominations. All3impulsive-NOISE blocks rejected;13->0signal nominations and0available noise exposure. Zero alarms under zero exposure is not an improvement. AtFAP0.005 same selection,lossesremain. Signals counted as window nominations, not waveform association; paired amplitudes/cases/FAPs are not independent observations.

Verification:24calibrations,2208evaluation records,270preliminary candidate records,4284regenerated inputs,12focused fixtures,78uncached parity checks. Official policy p/threshold null when unavailable; diagnostic core candidate records retained separately. No formal efficacy test or full-suite rerun, no admission-bias guarantee. No measurement failures/retuning or production/O3/real-data access. L09/L10 remain open.

Next: define an independent background-source control and realistic acquisition/provenance assumptions before another autonomous policy. Do not retune this heuristic on its failed pilot. New work remains local pending review; prior design/article commitd049e2b is on development.
