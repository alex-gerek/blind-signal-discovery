# Adaptive admission tracing — 27 November 2026

[Protocol](adaptive_admission_trace_protocol_2026-11-27.md) · [Full report](../results/adaptive_admission_trace/20261127T141606Z_ba44d7f3/report.md) · [Independent checks](../results/adaptive_admission_trace/20261127T141606Z_ba44d7f3/independent_verification.json).

Source d049e2b, unchanged production code. Both500-window arms exactly reproduce archived input hashes and all1000output records. Baseline later recovery0/6 versus oracle5/6, scores identical; thresholds163.83940094146237/1.6580395349128854 first used at367. No automatic admission fix was tested.

Added experimental calibration-call identity, ordered member hashes/A-B partitions/full calibration state and per-window used-state tracing. Switches0/162/234/367: availability returns234, before sensitivity collapses at367. Baseline415available noise windows/1660s,72unavailable/288s;7STFT alarms. Oracle419available noise windows because four injections became noise. Mature clean noise127/127broadband/fullUnion nominations remains unchanged.

Independent saved-data checks pass16calibrations,3108member occurrences,1000window records,856available p-values and7exposure fixtures. These are specific replay checks, not full-suite passes or independent scientific validation. No nonempty diagnostic candidates occur during abstention in this recipe; that preservation fixture remains missing. No measurement errors or retuning. L09/L10 unresolved; source and previous artifacts preserved, no O3 or real-data access.

Next: define one observable-only policy and freeze the successor experiment before evaluating it. Do not reject high scores by default or use oracle injection labels. No production integration, commit/push of this new stage, release or Zenodo action yet.
