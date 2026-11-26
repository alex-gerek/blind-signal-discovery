# Next release checklist after the post-audit controls

Updated 2026-11-26. This is a preparation record, not a published release or a Zenodo metadata update.

## Identity and deployment

- Owner-confirmed author form: **Lex Gerek** (2026-11-26). LICENSE and CITATION.cff already use this form.
- Existing Zenodo v0.3.0 creator metadata still says Alexander Gericke. For the next software release/deposit, use Lex Gerek consistently; do not infer identity information beyond the owner's chosen attribution.
- Public main: merge `5e19ae4`, source tree equal to development `043f0eb`. PR #8 and Pages deployment succeeded. The local documentation/control-study work after this merge is uncommitted and not deployed.
- Existing archive DOI `10.5281/zenodo.24557479` identifies tag v0.3.0 at `c406e00`. Do not move that tag or relabel old archive contents as current code.

## Required evidence before choosing release contents

1. Finish and independently account for the bounded [synthetic protocol](post_audit_workstation_control_protocol_2026-11-26.md). Keep failures, misses, unavailable exposure and partial Union states.
2. Perform the protocol's manual Python Workstation check with saved JSON/CSV exports. Headless checks do not satisfy the human/operator task.
3. Decide whether native browser codec/render/latency verification is required for the declared release scope. Otherwise state it remains unverified.
4. Resolve distribution attribution: third-party data/figures, dependency notices and icon source/license. MIT does not automatically license imported third-party material.
5. Review the local diff, result size and source hashes. Explicitly separate saved regression evidence (135 passed / 1 deselected) from the new pilot's experiment counts.
6. Historical rerun priority: adaptive recovery/exposure/sensitivity with valid denominators; coincidence with measured overlapping active intervals and identical foreground/background event definitions; then a new run of the previously used O2 GW170814 epoch. Preserve old output directories. O3 is excluded.

## Release mechanics, after evidence review

- Choose a new version/tag and exact commit; a patch version such as 0.3.1 is a proposal, not yet assigned.
- Update package/runtime version, citation metadata, README evidence links and release notes together.
- Use Lex Gerek in the next Zenodo creator metadata. Do not place the old version DOI on a new version as if it identified the new code. Preserve the old DOI in the history; update the version DOI once the new deposit identity is actually known.
- Review and save the release contents before merging/publishing them. Confirm the expected commit and completed checks.
- Publish one GitHub Release for the new tag only when ready. The enabled GitHub–Zenodo integration is expected to create the archive; verify actual deposit status, creator, version, license, GitHub revision and files afterward.
- No release, tag, commit, push, history rewrite or Zenodo operation was performed by writing this checklist.


## Evidence update after historical adaptive replay (2026-11-26)

The bounded headless pilot and six owner-operated80-control UI runs are saved. Synthetic historicalA/B replay is now complete, with a separate paired follow-up; see HISTORICAL_ADAPTIVE_REPLAY_FINDINGS_2026-11-26.md. L09 signal contamination and L10 frozen broadband/Union mismatch remain unresolved. Release notes must not advertise generally recovered adaptive sensitivity or changed-regime1% Union control. A repair requires a new design/protocol, not deleting known synthetic signals in production. Real-data historical reruns, native browser codec checks and attribution inventory remain separate. No new release/archive is published by this update.
