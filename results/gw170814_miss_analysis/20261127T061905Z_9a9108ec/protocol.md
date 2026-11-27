# Saved-data GW170814 miss analysis (2026-11-27)

Specified before deriving the new diagnostic tables/figures. Current clean source30849e8. Parent evidence:results/gw170814_current_control/20261126T174011Z_e21955d7.

Scope: analyze EXISTING saved cluster/calibration/provenance records and the three existing preprocessed event-window NPY files. No original strain reload, new event, detector run, calibration/threshold change, frequency-mask change or O3 access.

Question: where does each miss arise, and does the global winning cluster actually overlap the known event region? Initial inspection found the V1 global maximum at localt=3.125s,f=374Hz versus reference localt=1.5s. Confirm it rather than interpreting a near-threshold window maximum as near-recovery of the physical event.

- Preserve the parent's already specified near_reference Boolean: cluster interval intersects local[1.0,2.0]s and35–250Hz. This is a loose diagnostic association, not waveform identification or a new detection rule.
- Recompute all raw-cluster standardized scores from saved per-scale median/spread and confirm global maxima. Report both highest ANYWHERE and highest near_reference cluster for each scale/site.
- Translate each fixed global threshold back to its per-scale raw score boundary:median+threshold*spread; strict crossing required. Do not borrow another site's thresholds or optimize them.
- Check whether misses precede vetoes and whether any data windows were unavailable/dropped. Preserve all raw clusters and negative values.
- Summarize saved independently fitted background PSDs in35–250Hz as spectral descriptors, not event SNR, confidence intervals or proof of the physical cause. Undo the documented1e18 scale for ASD units; do not claim per-band median ASD alone predicts detection.
- Save full-window spectrograms and distinguish off-time winners from near-reference clusters. Original event arrays remain outsidegit; verify hashes. Any display contrast is shared, disclosed and does not alter measurement values.
- Confirmed: numerical decision mechanism, timing/frequency association and normalization. Hypotheses: instrument response, spectral masking, STFT fragmentation or background nonstationarity as physical causes. No causal intervention or new significance estimate.
- Fresh UTC/UUID output, protocol/driver/source/input hashes and independent arithmetic checks. Read-only originals. Block data reads except the exact three saved event-array paths, O3 and network.
- Also preserve the approved CodeQL dismissal record separately; no code fix or broad application-security claim.
- Save/push results per standing owner preference, without main merge/release/Zenodo.
