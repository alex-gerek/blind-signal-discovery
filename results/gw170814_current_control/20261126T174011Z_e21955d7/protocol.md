# GW170814 bounded current-code real-data replay (2026-11-26)

Specified before execution. One additional known O2 event, three instruments. Previously examined epoch; repeatability/control, not blind discovery, population efficiency or independent confirmation.

- Source8f78ca5b141ea656957dc69aa275daac4d71eeb0, no detector/preprocessing/veto tuning.
- ONLY existing H1/L1/V1 GW170814 GWOSC4KHZ_R1 files, GPS1186739814..1186743910. Check SHA256 against results/ligo_experiment_003/data_provenance.json. Original files outside git.
- Exact historical exploratory partition: PSD1186739816..1186740056 (60windows), CalA1186740080..1186741360 (320), target1186741540..1186742180 (160), CalB1186742360..1186743908 (387). Each4s window has1s context. Verify block footprint separation; stop on missing windows rather than stitching gaps.
- Existing adapted_ligo_whitened: PSD fit to PSDblock ONLY,1e18 scaling, zero-phase30–380Hz Butterworth, antialias4096->1024Hz decimation. No hard frequency-bin mask and no added waveform.
- Frozen MultiResolutionDetector,FAP.01,norm_fraction.30,seed0;707 calibration observations ->212normalization/495threshold. CalB is AFTER the target: explicitly offline analysis, not causal adaptive deployment.
- Known reference1186741861.5; on-sourcewindow80 GPS1186741860..1186741864. Exclude padded overlaps with +/-16s from off-source;151off-source,8guard and1event window.
- Preserve raw decisions, empiricp and candidates. Apply the ORIGINAL reversible VetoConfig(glitch=True,line=True), preserve pre/post counts and reasons.
- Miss diagnosis: save full-precision per-scale raw/standardized maxima and ALL seed-qualified raw clusters on window80, including subthreshold ones. A below-threshold score explains the software decision, not its entire physical cause. Do not assume antenna response, sensitivity or spectral noise is proven causal.
- Coarse foreground coincidence ONLY: deduplicate by window IDs and require peak-time matching (50msH1-L1,65msVirgo pairs), separately raw/post-veto. No time-slide background, waveform-delay inference, FAR/significance, adaptive Workstation or historical calibration variants this step.
- Compare raw/post-veto event flags, scores, thresholds and off-source counts with the saved historical report; explicitly note changes or recipe agreement. Original saved array hashes absent; cannot verify historical bitwise execution.
- Save fresh UTC/UUID directory, input/source/protocol/driver hashes, quality/injection annotations, full timeline/candidate/diagnostic records, figures and independent verification. Original and processed event waveform files stay ignored underdata/derived; hashes saved.
- Block O3/manifest, nonallowlisted raw data, historical writes and network. Import historical helper module ONLY; never call its main/report writers. Record failures.
- Display fixed event +/-0.5s,30–380Hz with common disclosed visible-bin contrast; no detector score/threshold change. Show selected raw below-threshold boxes separately from accepted candidates.
- Data providerGWOSC/LIGO/Virgo, CC BY4.0:https://gwosc.org/events/GW170814/ and https://gwosc.org/acknowledgement/. Project MIT does not supersede data license.
- Save/push the evidence on development per owner standing request; no main merge, tag/release or Zenodo operation. Article drafting deferred.
