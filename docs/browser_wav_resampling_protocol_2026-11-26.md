# Browser WAV resampling: design and validation record

- Step: `browser-wav-fir-2026-11-26`; source baseline `8fab987`, local patch recorded in research history.
- This is a preprocessing correction, not detector tuning or a general scientific validation. Earlier stride-decimated WAV measurements are not comparable without a new run. Historical synthetic Python/browser results are unchanged by this WAV-only path.

## Fixed design specified before testing

The old path selected every `round(decoded_fs/512)`-th sample without low-pass filtering. At a decoded 48 kHz, a 400 Hz tone became a near-unit 110.6383 Hz alias and the resulting rate was 510.6383 Hz.

The new deterministic stage maps the decoded samples to exactly 512 Hz. For source rate `f_s > 512`, it uses a centered rational Kaiser-windowed sinc FIR:

- Kaiser beta: 8.6; cutoff: 224 Hz.
- Useful passband: 0–192 Hz; transition: 192–256 Hz. This intentionally attenuates content near the new Nyquist limit; do not interpret the full 0–256 Hz band as flat.
- Half support: `H=ceil(24*f_s/512)` source samples; `2*H+1` coefficients per rational phase.
- At output index `m`, source position is `m*f_s/512`. For fractional distance `d` from a source sample, the unnormalized coefficient is

$$h(d)=\frac{448}{f_s}\,\operatorname{sinc}\!\left(\frac{448d}{f_s}\right)\frac{I_0\!\left(8.6\sqrt{1-(d/H)^2}\right)}{I_0(8.6)},\quad |d|\le H,$$

with normalized `sinc(x)=sin(pi*x)/(pi*x)`. Coefficients are normalized to unit DC gain separately for each rational phase. Samples outside the input are zero. The centered evaluation compensates nominal FIR delay; the first/last `H/f_s` seconds remain affected by padding.

The centered FIR is offline: it uses up to `H/f_s` seconds of future samples within the same file. Metadata explicitly mark `future_data_used=true`; it is not a causal streaming transform. Each input file is filtered separately, so no target samples enter the background filter. Filtering before window slicing introduces local dependence between adjacent windows; raw window count and an attainable empirical rank do not establish real-data FAP control. Edge-affected samples/windows are disclosed, not automatically excluded by this patch.

Output count is `floor(N*512/f_s)` and the unrepresented trailing duration is less than one output sample. Native decoded 512 Hz data are copied unchanged. Upsampling, nonfinite input, unsupported integer ratios exceeding two million phase coefficients, and insufficient length are rejected explicitly. At most two rate plans are cached. No waveform RMS normalization is applied.

## Original versus decoded sampling

The loader inspects the source RIFF/WAVE format and data chunks before decoding, then records both source-header and decoded rates. Web Audio may resample during `decodeAudioData`; its resampling and Float32 decode precision are separate from this FIR and can depend on the browser. Metadata include the container fields, decoder-rate change flag, exact output rate/count/duration, filter identity/parameters, DC normalization, zero padding and edge support.

Browser support in this patch is mono little-endian RIFF/WAVE accepted by the browser decoder, with original rate at least 512 Hz. Other containers, multichannel inputs, ambiguous/truncated chunks and costly ratios fail explicitly. Use the Python loader for other formats/channel selection. This does not certify every codec/browser combination.

## Acceptance checks and outcomes

The acceptance criteria were recorded in research history before implementation. On synthetic Float32 tones at decoded 44.1/48 kHz:

| Check | Prespecified criterion | Observed outcome |
| --- | --- | --- |
| 40, 100, 180, 192 Hz gain | Within 1% of unity, excluding padded edges | Passed at both rates |
| 256, 280, 400, 1000, 4096 Hz residual | Peak amplitude below `10^(-70/20)` away from edges | Passed at both rates |
| Old 400 Hz alias control | Reproduce the old near-unit alias, then suppress it | Old alias gain >0.99; new interior residual below the same -70 dB amplitude criterion |
| Timing | Center impulse at 0.5 s remains at output index 256 | Passed |
| Duration/immutability | Exact rate/count, truncation <1/512 s, no input mutation | Passed |
| Pipeline/provenance | Mixed original WAV rates yield 512 Hz and retain both rate/filter records | Passed with synthetic WAV/decoder fixtures |
| Independent reference | Absolute difference from matching SciPy FIR output <1e-5 at 1024/44100/48000 Hz | 3 pytest cases passed |

Tone checks trim 32 output samples at each end, beyond the filter's approximately 24-sample half support. They establish these test cases, not an attenuation guarantee over every frequency, input amplitude or codec. Edge artifacts are disclosed rather than treated as validated signal candidates.

`node --test tests/test_browser_input_validation.js` passed 15 scenarios, including the previous input/FAP/abstention/export checks. `tests/test_browser_wav_resampling.py` passed 3 independent reference cases in 2.16 s (`.tmp/wav-resample-reference-20261126T125501345901Z`). The reference constructs a FIR with SciPy `firwin` on the upsampled grid and applies `resample_poly`; global versus per-phase DC normalization is permitted within the stated tolerance. No test thresholds or filter parameters were changed after these successful checks.

The existing decoder-lifecycle test now supplies a valid RIFF fixture and verifies explicit native-rate metadata instead of expecting the old alias-warning flag. Its original mono rejection/context-close intent is retained. Native metadata, malformed containers, low rates and resource-limit errors have explicit checks.

## Limits and next verification

The browser tests execute production JavaScript with an in-memory DOM and synthetic decoder fixtures. They do not verify native browser decoding, rendering, every supported WAV encoding, or real-sensor false-alarm control. Python detector/adaptive/calibration code and historical result files are unchanged; O3 data and manifest contents are not used. A real-browser manual run should save the new provenance bundle and record browser/version, input hashes and any decoder transformations.

## Method references

- [Web Audio decodeAudioData specification](https://webaudio.github.io/web-audio-api/#dom-baseaudiocontext-decodeaudiodata): decoding may include resampling to the context rate.
- [Microsoft WAVEFORMATEX](https://learn.microsoft.com/en-us/windows/win32/api/mmreg/ns-mmreg-waveformatex): channel count, sampling rate, block alignment and sample format metadata.
- [SciPy resample_poly](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample_poly.html): independent zero-phase FIR/polyphase reference and padding semantics. This custom browser filter is not SciPy's default filter.


## Saved verification record and RIFF follow-up

The final combined reference/wrapper run passed 4 pytest cases in 8.86 seconds: three independent SciPy references plus the wrapper running the then-current 15 browser scenarios. Its environment, commands and production inline-script SHA-256 are preserved in [the unique validation record](validation/browser_wav_resampling_20261126T130026Z.json).

One additional RIFF boundary regression was then added without changing the production script: truncated payload/format/data chunks and duplicate format chunks are rejected; odd-length padded ancillary chunks are accepted. The complete Node suite passed **16 scenarios** in 6.2325 seconds. Filter parameters and numerical acceptance criteria were unchanged. The saved production script hash still matches the final working-tree implementation.

These timings describe this test environment, not an interactive browser performance guarantee. FIR work grows with recording length; large files may be better processed by the Python Workstation. Real-browser decoding and UI responsiveness remain manual-verification tasks.
