# Reproducing the admission studies — 27 November 2026

The studies were measured with source HEAD `d049e2b1c97626b562d234a22bc12db5d51702dd`. Research scripts enforce this source pin deliberately. Running them at a later publication HEAD will stop; this is not a detector regression.

Use a separate disposable clone/check-out at that pinned revision and the Python dependencies recorded in each run manifest. Preserve your active checkout and its uncommitted work. The historical synthetic replay/follow-up artifacts needed by the trace driver already exist at the pin. Do not use or open O3 or its reserved manifest.

Copy the experimental scripts and their protocol documents from the evidence revision into that clone WITHOUT changing its HEAD. At minimum:

- Trace: `experiments/run_adaptive_admission_trace.py`, `experiments/verify_adaptive_admission_trace.py`, and `docs/adaptive_admission_trace_protocol_2026-11-27.md`.
- Tail guard: `experiments/block_tail_guard.py`, `experiments/run_block_tail_guard.py`, `experiments/verify_block_tail_guard.py`, and `docs/block_tail_guard_protocol_2026-11-27.md`.
- Separate source: `experiments/run_separate_background_control.py`, `experiments/verify_separate_background_control.py`, `experiments/plot_separate_background_sources.py`, `experiments/summarize_separate_background_control.py`, and `docs/separate_background_control_protocol_2026-11-27.md`. It also imports CachedCore and the synthetic IO guard through the preceding tail-guard driver; copy that driver and helper too. The old tail rule is imported as a module dependency but is NOT applied in this source-control experiment.

Driver/helper copies and exact byte hashes are retained in each result directory; check them against the run manifest/index. Retain original pre-run protocols rather than silently editing them for a different source. Replays at newer detector source require a new protocol/run and an explicit comparison.

Run each measurement driver with `python -B experiments/<driver>.py`. It creates a fresh UTC/UUID directory and prints RUN_DIRECTORY. Pass that new directory to the corresponding verifier. For the separate-source study, run its plotting and summarizing scripts only after verification. The final archived report also contains a presentation clarification beyond the generic report builder; keep the archived report as the reviewed interpretation. Do not rerun a builder over the original saved run.

These are standalone experiment drivers, not pytest cases. The trace requires exact equality to the archived historical baseline; the other pilots regenerate their fixed synthetic streams. Verification checks arithmetic, hashes and scoped output semantics, not scientific exchangeability or validity for real sensors. The three seed comparisons are not large independent confirmatory studies.

Archived runs:

- `results/adaptive_admission_trace/20261127T141606Z_ba44d7f3/`
- `results/block_tail_guard/20261127T152327Z_22827dc3/`
- `results/separate_background_control/20261127T154355Z_4f134945/`

Do not overwrite these directories or normalize result line endings: artifact hashes refer to exact bytes. New runs must retain failures and all outcomes. Original detector, GUI, unrelated test cleanup and article drafts were excluded from the evidence commit.

## Additional lifetime/reacquisition pilot

This later study uses pinned source `e21a7cc48ce4122804f86e3397ee3894afa90639`, not d049e2b. Use a separate disposable checkout at e21a7cc and copy the new `background_lifetime_controller.py`, `run_background_lifetime.py`, `verify_background_lifetime.py`, `verify_lifetime_contamination_pair.py`, `report_background_lifetime.py`, and `plot_background_lifetime.py` into experiments/, plus `docs/background_lifetime_protocol_2026-11-27.md`. The preceding core-cache, source-control and synthetic guard dependencies already exist at this pin.

Run the measurement driver with python -B, retain its fresh RUN_DIRECTORY, run both lifetime verifiers on that directory, and only then generate reports/plots. Do not overwrite the archived `results/background_lifetime/20261127T170612Z_189c239c/`. The archived final report includes the reviewed presentation-error disclosure and figure caption clarification beyond the generic report builder; keep it as the reviewed interpretation.

The acquisition/expiry contract and example manifest are design documents, not enforcement in Workstation. The 2400s oldest-sample limit is an unvalidated experimental parameter. Preserve the saved plot-preparation error and negative outcomes alongside all successful arithmetic checks. This publication records a completed synthetic pilot, not a validated operational repair.
