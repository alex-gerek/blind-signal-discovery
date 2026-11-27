# Proposed protocol: blind two-tone structure estimation

Status: PROPOSED, not an executed benchmark. Owner requested review before any large experiment. Historical detector, reports and output hashes remain unchanged. No production modifications, O3 access, commit, push, merge or release. This is synthetic measurement research, not stem separation.

## Questions and boundaries

Separate three endpoints: (1) anomaly nomination, (2) a claim of two frequencies with estimates, (3) waveform reconstruction. Only the first two are tested. Two peaks need not mean two physical sources, and no separate source waveforms will be reconstructed. Exact coincident tones can equal one tone or cancel; they are identifiability controls, not recoverable two-frequency targets.

The detector and pair diagnostic receive waveform, sample rate and independently frozen calibration only. Carrier frequencies, number of tones, onset, phase and amplitudes are held in an evaluation-only manifest. Algorithms must scan the full retained frequency/time representation. No true-frequency neighborhoods, event-time crop or truth-based cluster selection. Truth is joined after outputs are frozen. A separate input-only entry point and a permutation test of evaluation labels will check this isolation.

## Comparison arms

Preserve the historical three-scale bank as a reference. Compare a 2 x 2 experimental layout:

- W0: unchanged bank, longest window 512 samples / 0.5 s, overlap 384.
- W1: same short/medium scales; replace longest window with 1024 samples / 1 s, overlap 896. This keeps the longest-scale hop at 128 samples / 0.125 s. Hann, no padding, seed Z>4, minimum two cells, eight-connectivity remain fixed. No zero-padding presented as improved physical resolution.
- D0: propose a frequency pair only from two accepted components with strictly disjoint frequency boxes in the longest scale; use their recorded peak frequencies.
- D1: retain D0 proposals and also examine multiple peaks within one accepted longest-scale component. It does not split production clusters or alter anomaly scores.

This factorial layout distinguishes effects of longer windows and within-cluster diagnostics. Historical uncalibrated box counts remain historical; the new blind D0 estimator and its false-pair calibration are additional research methods.

For the primary W1 contrast keep frequency/time median kernel lengths 21/7, identical to W0. The frequency smoothing footprint therefore changes from about 42 to 21 Hz, although the code parameter is fixed. This is an explicit normalization confound; any gain is a configuration-level result, not proof that window duration alone caused it. A physical-footprint-matched kernel-41 follow-up would be a separately calibrated secondary arm, not added after seeing results. W0/W1 longest-scale frequency grids are 2/1 Hz; Hann resolving power is broader than grid spacing.

## Proposed within-cluster diagnostic and pair score

Before coding, freeze this input-only recipe or amend the protocol prospectively. For each accepted longest-scale component, form the mean Z spectrum over its occupied time columns (all columns selected from the component, not the true event). Search strict local frequency maxima inside that component's frequency extent. Require two maxima of height >=4, an interior sampled minimum <=0.8 times the smaller peak, and at least one interior frequency sample. No two adjacent bins are called resolved peaks. Use bin-center frequencies; no hidden interpolation or true-carrier fitting. Retain every proposal and reason for rejection.

D0 proposals use disjoint accepted components with overlapping time boxes, restricting accidental sequential tones being called simultaneous. D1 inherits these and adds within-component proposals. Proposed scalar pair score is the smaller normalized peak height; for D0 these are the component peak Z values, for D1 the two mean-profile peak values. The differing score distributions are handled by separate calibration. Choose the highest-score proposal across the full observation; ties broken by ascending frequencies then earliest time. Output at most one pair, or abstain. A missing/invalid proposal scores minus infinity. Require an available anomaly decision and at least one accepted longest-scale component. Mean-profile averaging can erase changing ridges and phase beating; record that limitation rather than adjusting it on evaluation data.

Diagnostic gates (height 4, valley 0.8, time-overlap condition) are hypotheses, not validated settings. Implement and unit-check them on development fixtures; change only with an amended protocol and fresh held-out evaluation seeds. Never choose gates using the already inspected 96/98 case and then call that case independent validation.

## Data and splits

Use independent RNG streams and persistent unique run folders. Separate development fixtures, anomaly calibration noise, pair-claim calibration controls, held-out controls and held-out two-tone targets. W0/W1/D0/D1 see the same evaluation waveforms for paired comparisons. Within a split, generate independent observation noise for each example; do not reuse one noise record across different labels as if independent. Randomized generator metadata must be saved with waveform hashes.

Maintain 4-second records at 1024 Hz. Preserve the old one-second Tukey-0.15 tone envelope as the primary duration. Randomize onset uniformly within 0.75..1.75 s, keeping its full support in the record. Additional two-second tones are a separately reported duration control; increased coherent energy must not be confused with improved resolution. Keep a small historical 96/98 reproduction set as development-only, outside held-out metrics.

Proposed held-out pilot: 200 noise-only and 200 single-tone examples, plus 216 two-tone examples from gaps {1,2,4,8,16,96} Hz x relative phases {0,pi/2,pi} x amplitude ratios {1,0.5,0.25} x four independent noise/carrier draws. This is too small for strong low-false-rate conclusions or finely stratified power; a larger confirmatory sample needs a prospective power calculation after the pilot. No pilot run authorized by this document alone.

Randomize base carriers continuously within 80..180 Hz, not just FFT bins, and global phases uniformly. Include weak and strong levels (root-energy SNR 8 and 18) balanced across the 216 cases, with assignments fixed before generation. For each pair, scale the weaker component by the amplitude ratio and rescale the summed addition to the assigned resultant root-energy SNR. Save component energies and cross term; this controls resultant energy across gaps but changes component amplitudes when phase/correlation changes. Exclude no difficult example silently. Zero resultant energy is an explicit non-observable control rather than divided by zero. Additional unrenormalized historical-SNR-18 controls, if run, are separate descriptive results.

One-tone controls match duration, onset, carrier range and resultant energy distribution of two-tone targets. Include off-bin carriers and stronger single tones to expose sidelobe splitting. Keep exact coincident/cancellation cases as separate identifiability controls. Single-tone and noise examples must be independent of calibration and each other.

Use pink noise matching the historical generator in the primary pilot. White noise, spectral lines and changed noise regimes require separate predeclared strata and suitable matched calibration. Primary results cannot establish robustness to these untested regimes or real sensors.

## Independent calibration

For each W configuration, use its own calibration computed from 1,000 independent noise records (A=300, B=700, target anomaly FAP 0.005). W0/W1 may share calibration waveform IDs for pairing, but scores, normalization and thresholds must be recomputed independently. Neither reuses the old threshold. Freeze and hash before scoring evaluation inputs.

For every W/D arm, independently calibrate the full-observation maximum pair score on two additional negative-control sets: 1,000 noise-only and 1,000 one-tone observations. One-tone controls span the declared amplitude/onset/off-bin-carrier distribution. Set each group's 1% strict-exceedance cutoff using the finite-sample order rank ceil((n+1)*(1-0.01)); take the larger of the two cutoffs. Include no-proposal minus-infinity observations in the null distribution. Thresholds must not depend on true two-tone frequencies or held-out results. This controls marginal claims for the specified calibration populations under exchangeability, not worst-case control over all tones, phases or sensors. Predefine unavailable behavior if a finite attainable threshold cannot be obtained. No selection of the best W/D arm on held-out data followed by an uncorrected significance claim.

## Metrics

Report anomaly nominations, available/abstained decisions, any two-frequency claim, correct matched pair, wrong pair and no pair for every arm and stratum. Count false two-tone claims separately on noise-only and single-tone examples; do not report only conditional accuracy among nominated anomalies. Report per-method binomial intervals for independent examples and paired discordances between arms; stratified cells are exploratory with very small denominators. The pilot is not a proof of superiority or a sub-percent false-pair certification.

For correct pair identification, require a one-to-one assignment of the two returned frequency estimates to the two true carriers, each within 1 Hz (fixed for both windows). Report absolute estimation errors, separation error and how many matches are ambiguous because tolerance neighborhoods overlap. At gaps <=2 Hz, add a stricter descriptive tolerance gap/4 to avoid counting arbitrary near-midpoint peaks as accurate recovery; do not replace the predeclared primary score afterward. Exactly coincident tones have no distinct-frequency target. Anomaly detection alone never counts as a successful pair estimate.

## Stage attribution and saved artifacts

For every input/scale save raw power, complex STFT, normalization factors, normalized Z, frequency/time axes, seed mask, full labeled components including size-rejected ones, raw/standardized scores, acceptance/rejection reason, expanded boxes, spectra/peaks/valleys/pair proposals and final pair decision. Preserve the original waveform or a verified deterministic regeneration recipe, source snapshot, dependency versions, calibration hashes and runtime logs. Known frequencies are used only by a separate evaluator to annotate where evidence disappeared: representation, seed connectivity, size rejection, score acceptance, box rule or pair gate.

For the existing seed-22101 map, Z and clustering can be inspected without recomputation. For raw-power/phase/background attribution and other units, regenerate at the pinned source, verify every historical waveform hash, and archive all intermediate representations to a NEW directory. Do not overwrite historical maps or reinterpret that reproduction as independent evaluation.

## Timing and runtime

Keep hop fixed in the primary W contrast. Measure single-thread warm and cold wall time separately for STFT/normalization, clustering, calibration and peak diagnostic; report median/p95 per 4-second observation, total calibration cost, peak memory, machine/library versions and repeated measurement methodology. Calibration cost is not online inference cost.

Without padding, W0 has 29 longest-scale frames and W1 has 25. Their first frame centers are 0.25/0.5 s, with half-window support 0.25/0.5 s. These support differences are arithmetic, not measured event-timing error. Record onset/offset box error, peak-time error and their tolerance definitions. A 1-second window may broaden event localization or dilute short events. Existing full-record normalization and four-second scoring are offline; do not advertise hop duration as streaming latency. Actual causal implementation latency remains untested.

## Execution gates and limitations

First review this proposal. Then implement isolated research scripts and focused fixtures without production changes, perform a tiny smoke run with a runtime cap, and report its errors/results before proposing the frozen pilot. The large calibration sets dominate work even though the evaluation pilot is small; do not launch them implicitly. Never trade an apparent resolution gain for unmeasured false two-tone claims. No reconstruction, stem separation, physical-source identification or real-data benefit is established by this protocol.
