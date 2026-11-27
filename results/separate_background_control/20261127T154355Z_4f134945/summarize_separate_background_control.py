"""Summarize verified source-control evidence without rerunning measurements."""
from pathlib import Path
import csv,json,sys
out=Path(sys.argv[1]);v=json.loads((out/"independent_verification.json").read_text());assert v["status"]=="PASS"
arms=["stream_contaminated","paired_noise_oracle","separate_matched","separate_old"]
cases=["stationary","level_step","spectral_line","impulsive_noise"]
keys=["chirp_snr8.0","chirp_snr18.0","burst_snr8.0","burst_snr18.0"]
rows=[]
for alpha in [.01,.005]:
 for case in cases:
  for arm in arms:
   g=v["groups"][f"{case}_{arm}_fap{alpha}"]
   rows.append({"alpha":alpha,"case":case,"arm":arm,"noise_alarms":g["noise_alarms"],"noise_windows":g["noise_windows"],"available_noise_seconds":g["available_noise_seconds"],**{k:g["signals"][k]["associated"] for k in keys}})
with (out/"comparison_counts.csv").open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
lines=["# Separately acquired background: source-control pilot", "",f"Run `{out.name}`, source d049e2b; unchanged detector core. [Protocol](protocol.md), [independent verification](independent_verification.json), [counts CSV](comparison_counts.csv).", "",
"## Design and scope", "",
"Three fresh seeds18101..18103, four noise profiles, four265-window background sources, FAP0.01/0.005. Each source has66normalizationA/199thresholdB. All eight calibrations per seed/profile frozen before target decisions; all targets identical across arms. No tail filtering or quarantine. Every calibration is available at both FAPs; there are no abstentions in this pilot.","",
"Stream_contaminated adds two chirps/two bursts at indices100/140/180/220 inB. Paired_noise_oracle is exactly that background before those additions; inaccessible in real measurements. Separate_matched is a fresh independently sampled, signal-free background from the same generated noise law as the target. Separate_old is another independent stationary pink draw, intentionally mismatched in the other profiles. Signal-free matching is guaranteed by the synthetic construction, NOT discovered by an admission rule.","",
"The paired oracle shares allA and differs in exactly fourB inputs. Separate_matched replacesA/B and therefore includes ordinary calibration-sampling variation; differences from oracle are not evidence of a second physical mechanism. Level_step is a fixed post-step snapshot (2xnoise), not a test of transition timing, continuous drift or automatic adaptation.","",
"Each profile has600negative windows across three calibration draws,2400available seconds per arm, and30units per fixed family/SNR (120paired-amplitude signal records,60distinct noise/clean units). Across profiles,240distinct signal units/480amplitude records. Arms and FAPs reuse targets; paired amplitudes share noise. Do not count repetitions as independent samples. Support association means accepted rectangle intersects the known signal time/frequency rectangle; it is weaker than waveform identification.","",
"## Results", "",
"Counts below are descriptive aggregates over three calibration seeds, not a pooled independent-binomial FAR guarantee or a formal superiority test. Read association counts beside noise alarms: a mismatched background may nominate noise and produce misleading apparent recovery.",""]
for alpha in [.01,.005]:
 lines += [f"### FAP target {alpha}","","| Noise profile | Calibration source | Noise alarms / 600 | Chirp SNR8 / 30 | Chirp SNR18 / 30 | Burst SNR8 / 30 | Burst SNR18 / 30 |","| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
 for row in rows:
  if row["alpha"]==alpha:lines.append("| "+" | ".join(str(row[k]) for k in ["case","arm","noise_alarms",*keys])+" |")
 lines += [""]
lineold=v["groups"]["spectral_line_separate_old_fap0.01"]
lines += ["## Interpretation and limits","",
f"The deliberately old background in the line profile nominates {lineold['noise_alarms']}/600noise windows at nominalFAP0.01. Its positive-window associations therefore cannot be interpreted as successful signal identification. A signal-free independent source is not enough if its noise distribution fails to match the target.","",
"A separate matched source tests whether suitable independently acquired background can avoid signal contamination without selecting quiet tails or losing exposure. It is not an autonomous admission solution. Comparison with the oracle is diagnostic only; sample-to-sample normalization/threshold variation remains. Any claimed operational improvement still needs a common validated alarm budget, more independent calibration/sequence units and application-specific acquisition requirements.","",
"This study does not establish that the chosen profiles represent a real sensor, that any real off-source segment is pure/stationary, or that a second channel can substitute for the target channel. Comparable preprocessing, gain, units, state and contamination risk must be justified for real background acquisition. No new astrophysical result, GW170814/H1 explanation, frozen broadband/Union repair or adaptive recovery inference.","",
"## Independent verification and preservation","",
f"Checks pass: {v['calibrations']}calibrations, {v['evaluation_records']}target/arm/FAP records, {v['candidate_records']}candidate records, {v['regenerated_inputs']}regenerated calibration/target input records, and {v['uncached_parity_checks']}direct uncached raw-feature parity checks. Independent arithmetic reconstructsAstatistics,Bmaxima/ranks,nulls,alltargetscores/empiricalp/strictdecisions,candidate standardization/bounds/support association,paired gains/losses and actual exposure. Frozen hashes and input-source separation match. Generation replay proves preservation, not scientific validity or PRNG independence.","",
"No measurement failures, interrupted runs, case exclusions, score trimming, parameter retuning or production edits. All original artifacts retained; source/dependency hashes unchanged. No O3/data/manifest/network or real-waveform access. These are scoped experiment checks, not full-suite pytest passes. The original detector remains default; L09/L10 stay unresolved.","",
"![Noise-window nominations](background_source_noise.png)","",
"Noise figure: panel y-ranges differ; denominators600perprofile/arm, same nominalFAP0.01, not equal achieved FAR.","",
"![Known-support associations](background_source_signals.png)","",
"Signal figure:30units per family/SNR/profile; amplitudes and arms paired. View alongside noise alarms.","",
"## Next decision","",
"Specify a realistic source-acquisition contract before attempting autonomous admission: what an instrument can record separately, what makes it comparable, what contamination information is externally available, when it is stale and when to abstain. A new sequence-level experiment must test mismatch, availability and actual recovery under a frozen protocol. Do not promote this oracle-backed synthetic pilot into operational certification."]
(out/"report.md").write_text("\n".join(lines)+"\n",encoding='utf-8')
print('Saved descriptive tables and report from independently verified data.')
