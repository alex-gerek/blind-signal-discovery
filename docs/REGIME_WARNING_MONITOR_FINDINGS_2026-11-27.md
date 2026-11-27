# Warning-only regime monitor findings — 27 November 2026

[Pre-run protocol](regime_warning_monitor_protocol_2026-11-27.md) · [Full report/figures](../results/regime_warning_monitor/20261127T173116Z_9d6affa3/report.md) · [Independent checks](../results/regime_warning_monitor/20261127T173116Z_9d6affa3/independent_verification.json).

Source `6f4ccfa`, unchanged default `BackgroundRegimeMonitor`, warning-only sidecar. Five new continuous-sequence seeds, six paired cases, frozen reference and search calibration. No background updates, candidate suppression or production integration.

## Measured behavior

- Both chosen changes (level step and persistent spectral line) start a new warning on the first affected four-second observation in all five sequences, with no active warning immediately before the change. This is window-based notification timing, not a live-runtime benchmark or general detection guarantee.
- Stable noise produces five flagged windows in four episodes, across two of five sequences (2525 monitor-eligible windows). Flagged windows are not independent alarms; the nominal KS cutoff does not certify a long-horizon false-warning rate.
- Sparse transients: six flagged windows and five episodes. Repeated transients: seven flagged windows and six episodes. A hypothetical warning gate would remove2/47 and2/304 existing support-associated search recoveries, respectively. No such gate was applied.
- An intended continuous tone and the noise-line case have identical waveforms and warnings. The monitor identifies altered measurements, not physical origin. The persistent case represents five sustained episodes, not2200independent signals. Do not absorb these features into background automatically.
- Supplemental paired gain check: all2580stable/2xlevel-step search scores,p-values and nomination decisions are exactly equal, although monitoring warns. Raw-data change need not imply a changed normalized search statistic; this is a property of this factor-two control, not arbitrary gain-response invariance.

Strict pilot screens for zero stable episodes and zero extra transient warning windows fail; the change-delay screen passes. These exploratory criteria do not prove a mathematical bug or establish a deployment tolerance. Keep the monitor as a provisional diagnostic annotation, not a validated noise classifier, gate or recalibration trigger.

## Protocol qualification

The pre-run text inconsistently requires alternating families within each transient schedule AND identical injected waveforms at shared positions. At344/704/784 these conflict; the script follows within-case alternation and different families occur at those shared positions. [Audit record](../results/regime_warning_monitor/20261127T173116Z_9d6affa3/protocol_inconsistency.json) preserves the issue. Comparisons of each case with its stable base remain descriptive; sparse-versus-repeated differences cannot be attributed solely to cadence. No silent protocol edit, rerun or seed replacement.

## Verification and scope

Checks pass5reference models,5search calibrations,15480monitor/search/sidecar records each,53711candidate records and24000regenerated input-window records.105independent Hann/NumPyFFT feature checks,90production prefix comparisons,90900independent KSfeature comparisons (maxBonferronipdifference~1e-15),and32small exact-enumeration oracle checks. Reference/source hashes unchanged; every candidate payload and search result preserved by hashes. These are scoped checks, not full pytest passes or proof of scientific validity.

One verifier syntax-preparation error was caught before execution and retained separately; measurements were not restarted. Graphs visually checked. No measurement interruptions, retuning, excluded cases, O3/data/manifest/network/real-waveform access, production/test edits or publication. Continuous colored windows are dependent; five paired seeds do not establish a real sensor monitoring budget. L09/L10 remain open. New results are local.

Next: define an application-specific warning budget and distinguish generic measurement changes from changes relevant to detector calibration. Any revised rule needs fresh independent sequences and a consistent signal schedule; preserve persistent-signal ambiguity and all candidates.
