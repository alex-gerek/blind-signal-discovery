"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const crypto = require("node:crypto").webcrypto;
const html = fs.readFileSync(path.join(__dirname, "../docs/workstation/index.html"), "utf8");
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1].replace(/\nhandleRunMeasurement\(\);\s*$/, "");
function setup() {
  const elements = new Map();
  for (const id of html.matchAll(/id="([^"]+)"/g)) elements.set(id[1], {
    value: "0", files: [], style: {}, disabled: false, innerHTML: "", textContent: "", width: 960, height: 240,
    addEventListener() {}, getContext() { return {clearRect() {}}; }
  });
  for (const tag of html.matchAll(/<input\b([^>]*)>/g)) {
    const id=tag[1].match(/id="([^"]+)"/), value=tag[1].match(/value="([^"]*)"/);
    if (id && value) elements.get(id[1]).value=value[1];
  }
  for (const tag of html.matchAll(/<select[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)) {
    const options=[...tag[2].matchAll(/<option([^>]*)value="([^"]+)"([^>]*)>/g)];
    const selected=options.find(o => /selected/.test(o[1]+o[3])) || options[0];
    elements.get(tag[1]).value=selected[2];
  }
  const tbody = new Map();
  const document = { getElementById: id => elements.get(id) || null,
    querySelector: id => { if (!tbody.has(id)) tbody.set(id, {innerHTML: ""}); return tbody.get(id); },
    querySelectorAll: () => [] };
  const ctx = vm.createContext({console, crypto, document, window: {}, setTimeout: f => {f();return 0;}});
  vm.runInContext(script, ctx);
  vm.runInContext("drawTimeline=()=>{};drawWindowDetail=()=>{};", ctx);
  return {ctx, elements};
}
function file(name, text) { const bytes = new TextEncoder().encode(text); return {name, text: async () => text, arrayBuffer: async () => bytes.buffer}; }
function wavFile(name, rate, samples, channels=1) {
  const buffer=new ArrayBuffer(44+samples.length*4), view=new DataView(buffer);
  function ascii(offset,text) {for(let i=0;i<text.length;i++)view.setUint8(offset+i,text.charCodeAt(i));}
  ascii(0,"RIFF");view.setUint32(4,buffer.byteLength-8,true);ascii(8,"WAVE");ascii(12,"fmt ");view.setUint32(16,16,true);
  view.setUint16(20,3,true);view.setUint16(22,channels,true);view.setUint32(24,rate,true);view.setUint32(28,rate*channels*4,true);
  view.setUint16(32,channels*4,true);view.setUint16(34,32,true);ascii(36,"data");view.setUint32(40,samples.length*4,true);
  samples.forEach((value,i)=>view.setFloat32(44+i*4,value,true));return {name,arrayBuffer:async()=>buffer.slice(0)};
}
function tone(rate,frequency,duration=2,phase=.37) {return Float32Array.from({length:Math.floor(rate*duration)},(_,i)=>Math.sin(2*Math.PI*frequency*i/rate+phase));}
function fittedAmplitude(signal,frequency,rate=512,trim=32) {
  let ss=0,cc=0,sc=0,ys=0,yc=0;
  for(let i=trim;i<signal.length-trim;i++) {const a=2*Math.PI*frequency*i/rate,si=Math.sin(a),co=Math.cos(a);ss+=si*si;cc+=co*co;sc+=si*co;ys+=signal[i]*si;yc+=signal[i]*co;}
  const d=ss*cc-sc*sc;return Math.hypot((ys*cc-yc*sc)/d,(yc*ss-ys*sc)/d);
}
function interiorPeak(signal,trim=32) {let peak=0;for(let i=trim;i<signal.length-trim;i++)peak=Math.max(peak,Math.abs(signal[i]));return peak;}
function pairData(start=0) { return Array.from({length:20*1024}, (_,i) => `${start+i/512},${Math.sin(i*.113)}`).join("\n"); }
function setValues(e, values) { for (const [key,value] of Object.entries(values)) e.get(key).value=String(value); }

test("strict CSV rejects corrupt samples, gaps and invalid time axes", () => {
  const {ctx}=setup();
  for (const text of ["1\nNaN\n2", "1\nInfinity\n2", "1\n\n2", "1,\n2,3", "0,1\n.1,2\n.3,3", "0,1\n0,2", "1\nnope\n2", "NaN\n1\n2"]) {
    assert.throws(() => ctx.parseBrowserCsv(text,512));
  }
  for (const rate of [0,-1,Infinity,NaN]) assert.throws(() => ctx.parseBrowserCsv("1\n2",rate));
  const valid=ctx.parseBrowserCsv("time,value\n0,1\n0.001953125,2\n0.00390625,3\n",512);
  assert.deepEqual(Array.from(valid.signal),[1,2,3]);assert.equal(valid.fs,512);assert.equal(valid.preprocessing.discarded_data_rows,0);
  assert.throws(() => ctx.parseBrowserCsv("0,1\n.001953125,2",1024),/conflicts/);
  assert.throws(() => ctx.parseBrowserCsv("value\n1\n1",512),/constant/);
});

test("pair checks reject leakage, overlap and mixed effective rates", () => {
  const {ctx}=setup(); const a=ctx.parseBrowserCsv("1\n2\n3",512);
  assert.throws(() => ctx.validateBrowserPair(a,a,"same","same"),/identical/);
  assert.throws(() => ctx.validateBrowserPair(a,{...a,fs:1024},"a","b"),/sample rates/);
  assert.equal(ctx.validateBrowserPair(a,a,"a","b").status,"UNVERIFIED_EXTERNAL_CHECK_REQUIRED");
  assert.throws(() => ctx.validateBrowserPair({...a,timeInterval:[0,2]},{...a,timeInterval:[1,3]},"a","b"),/overlap/);
  assert.equal(ctx.validateBrowserPair({...a,timeInterval:[0,2]},{...a,timeInterval:[3,5]},"a","b").status,"VERIFIED_DISJOINT_TIMESTAMPS");
});

test("real file workflow hashes inputs, exports provenance and retains selected window", async () => {
  const {ctx,elements:e}=setup();setValues(e,{workflow:"files",targetFap:.05,csvFs:512,winIdx:0});
  e.get("bgFile").files=[file("same-name.csv",pairData())];e.get("tgtFile").files=[file("same-name.csv",pairData(100))];
  await ctx.handleRunMeasurement();
  const exported=JSON.parse(vm.runInContext("JSON.stringify(buildExportPayload(LAST_RUN))",ctx));
  assert.equal(exported.metadata.independence_audit.status,"VERIFIED_DISJOINT_TIMESTAMPS");
  assert.equal(exported.metadata.input_provenance.target.sha256.length,64);
  assert.equal(exported.summary.n_target_windows,20);assert.equal(e.get("dlJson").disabled,false);
  e.get("winIdx").value="7";await ctx.handleRunMeasurement();assert.equal(e.get("winIdx").value,7);
  // Failure after success must not offer the old run as a new result.
  e.get("tgtFile").files=e.get("bgFile").files;await ctx.handleRunMeasurement();
  assert.equal(vm.runInContext("LAST_RUN",ctx),null);assert.equal(e.get("dlJson").disabled,true);
  assert.match(e.get("statusLine").textContent,/identical SHA-256/);assert.equal(e.get("runBtn").disabled,false);
});

test("WAV decoder failures close AudioContext and native-rate metadata is explicit", async () => {
  const {ctx}=setup();let closes=0;
  ctx.window={AudioContext:class {async decodeAudioData(){return {numberOfChannels:2};}async close(){closes++;}}};
  await assert.rejects(() => ctx.readAudioOrCsv(wavFile("decoder-mismatch.wav",512,[1,2,3,4]),512),/mono/);assert.equal(closes,1);
  ctx.window={AudioContext:class {async decodeAudioData(){return {numberOfChannels:1,sampleRate:512,getChannelData:()=>Float32Array.from([1,2,3,4])};}async close(){closes++;}}};
  const loaded=await ctx.readAudioOrCsv(wavFile("mono.wav",512,[1,2,3,4]),512);assert.equal(loaded.preprocessing.resampled,false);assert.match(loaded.preprocessing.antialias_filter,/not required/);assert.equal(loaded.preprocessing.original_sample_rate_hz,512);assert.equal(closes,2);
});


test("calibration validation terminates on invalid settings and preserves exact FAP boundaries", () => {
  const {ctx}=setup();
  for (const alpha of [0,1,-1,NaN,Infinity,Number.MIN_VALUE,1e-20]) assert.throws(() => ctx.conformalOrderRank(60,alpha));
  for (const n of [0,-1,1.5,NaN,Infinity,1e30]) assert.throws(() => ctx.conformalOrderRank(n,.01));
  function adjacent(x,up) {const view=new DataView(new ArrayBuffer(8));view.setFloat64(0,x);view.setBigUint64(0,view.getBigUint64(0)+(up?1n:-1n));return view.getFloat64(0);}
  for (const n of [40,60,84,99,154]) {
    const boundary=1/(n+1);
    assert.equal(ctx.conformalOrderRank(n,adjacent(boundary,false)).attainable,false);
    for (const alpha of [boundary,adjacent(boundary,true)]) {
      const rank=ctx.conformalOrderRank(n,alpha);assert.equal(rank.attainable,true);
      const nulls=Float64Array.from({length:n},(_,i)=>i);
      assert.ok(ctx.empiricalPValue(nulls,n).p<=alpha);
      assert.ok(ctx.empiricalPValue(nulls,n-1).p>alpha); // tie at the maximum cannot pass
    }
  }
});

test("default synthetic workflow still works and preserves selection", async () => {
  const {ctx,elements:e}=setup();
  await ctx.handleRunMeasurement();
  let result=JSON.parse(vm.runInContext("JSON.stringify(buildExportPayload(LAST_RUN))",ctx));
  assert.equal(result.metadata.workflow,"synth");assert.equal(result.metadata.n_bg_windows,120);
  assert.equal(result.summary.n_target_windows,50);assert.equal(result.summary.n_detection_available_windows,50);
  assert.equal(result.metadata.input_provenance,null);
  e.get("winIdx").value="13";await ctx.handleRunMeasurement();assert.equal(Number(e.get("winIdx").value),13);
  assert.equal(e.get("dlJson").disabled,false);
});


test("unattainable self-test reports abstentions and unavailable rates in UI and JSON", async () => {
  const {ctx,elements:e}=setup();setValues(e,{workflow:"synth",targetFap:.005,nBg:40,nNeg:10,nPos:2});
  await ctx.handleRunMeasurement();
  const payload=JSON.parse(vm.runInContext("JSON.stringify(buildExportPayload(LAST_RUN))",ctx));
  const gt=payload.ground_truth_audit;
  assert.equal(gt.negative_control_available_exposure_s,0);
  assert.equal(gt.negative_control_far_available,null);assert.equal(gt.negative_control_far_available_ci90,null);
  assert.equal(gt.positive_injections_abstained,2);assert.equal(gt.positive_injections_missed_available,0);
  assert.equal(gt.detection_efficiency_available,null);
  assert.equal(gt.negative_control_far,0); // legacy field explicitly retained
  assert.equal(gt.positive_injections_missed,2); // legacy non-detection count
  assert.ok(payload.timeline.every(r => r.evaluation_status==="ABSTAINED" && !r.active_mode_anomaly_flag && r.threshold_full_precision===null));
  assert.match(e.get("gtMetrics").innerHTML,/UNAVAILABLE/);
  assert.match(ctx.document.querySelector("#gtTable tbody").innerHTML,/ABSTAINED/);
  assert.doesNotMatch(ctx.document.querySelector("#gtTable tbody").innerHTML,/MISSED/);
  assert.match(e.get("jumpButtons").innerHTML,/Abstained Injection/);
  assert.doesNotMatch(e.get("jumpButtons").innerHTML,/Missed Injection/);
});

test("mixed/empty availability uses explicit denominators and conditional exposure", () => {
  const {ctx}=setup();
  const windows=[
    {index:0,inj:null,detectionAvailable:false,activeFlag:false},
    {index:1,inj:null,detectionAvailable:true,activeFlag:true},
    {index:2,inj:null,detectionAvailable:true,activeFlag:false},
    {index:3,inj:{id:0},detectionAvailable:true,activeFlag:true},
    {index:4,inj:{id:1},detectionAvailable:true,activeFlag:false},
    {index:5,inj:{id:2},detectionAvailable:false,activeFlag:false}
  ];
  const gt=ctx.summarizeBrowserGroundTruth(windows,2);
  assert.equal(gt.negative_control_far_available,.5);assert.equal(gt.detection_efficiency_available,.5);
  assert.equal(gt.negative_control_available_exposure_s,4);assert.equal(gt.positive_injections_missed_available,1);
  assert.equal(gt.positive_injections_abstained,1);
  assert.deepEqual(Array.from(gt.injection_records,r=>r.evaluation_status),["DETECTED","MISSED","ABSTAINED"]);
  const empty=ctx.summarizeBrowserGroundTruth([],2);assert.equal(empty.negative_control_far_available,null);assert.equal(empty.detection_efficiency_available,null);
  assert.throws(()=>ctx.summarizeBrowserGroundTruth([{inj:null,detectionAvailable:false,activeFlag:true}],2),/require/);
});

test("precise decision inputs survive CSV and JSON while legacy display values stay rounded", () => {
  const {ctx}=setup();const threshold=2.888922;
  const view=new DataView(new ArrayBuffer(8));view.setFloat64(0,threshold);view.setBigUint64(0,view.getBigUint64(0)+1n);const score=view.getFloat64(0);
  const payload=ctx.buildExportPayload({metadata:{},calibration:{threshold},summary:{},windows:[{index:0,role:"REVIEW",maxScore:score,detectionAvailable:true,diagnosticNomination:true,activeFlag:true,nominatedCandidates:[{}]}],candidates:[]});
  const row=JSON.parse(JSON.stringify(payload)).timeline[0];
  assert.equal(row.max_score,row.threshold);assert.ok(row.max_score_full_precision>row.threshold_full_precision);
  assert.equal(row.max_score_full_precision,score);assert.equal(row.threshold_full_precision,threshold);
  const csv=ctx.toCsv(payload.timeline).split("\r\n"),fields=csv[0].split(","),values=csv[1].split(",");
  assert.equal(Number(values[fields.indexOf("max_score_full_precision")]),score);
  assert.ok(Number(values[fields.indexOf("max_score_full_precision")])>Number(values[fields.indexOf("threshold_full_precision")]));
  assert.throws(()=>ctx.buildExportPayload({...{metadata:{},calibration:{threshold},summary:{},candidates:[]},windows:[{...{index:0,role:"INVALID",detectionAvailable:true,diagnosticNomination:false,activeFlag:false,nominatedCandidates:[]},maxScore:NaN}]}),/non-finite/);
  assert.match(ctx.toCsv([],vm.runInContext("BROWSER_CANDIDATE_FIELDS",ctx)),/score_full_precision/);
  assert.equal(ctx.toCsv([{label:'a,"b"\nc'}]),'label\r\n"a,""b""\nc"');
});

test("actual detected candidate exports retain full score, p-value and threshold", async () => {
  const {ctx,elements:e}=setup();e.get("snr").value="30";await ctx.handleRunMeasurement();
  const payload=JSON.parse(vm.runInContext("JSON.stringify(buildExportPayload(LAST_RUN))",ctx));
  assert.ok(payload.candidates.length>0);
  for(const c of payload.candidates) {
    assert.ok(c.score_full_precision>c.threshold_full_precision);
    assert.ok(c.empirical_p_value_full_precision<=payload.calibration.target_fap);
    assert.equal(c.threshold_full_precision,payload.calibration.threshold);
  }
});


test("WAV FIR passband and stopband meet the prespecified tone checks", () => {
  const {ctx}=setup();
  for(const rate of [44100,48000]) {
    for(const frequency of [40,100,180,192]) {
      const output=ctx.resampleBrowserWav(tone(rate,frequency),rate);
      const gain=fittedAmplitude(output.signal,frequency);
      assert.ok(Math.abs(gain-1)<.01,`${rate} Hz / ${frequency} Hz passband gain ${gain}`);
      assert.equal(output.fs,512);assert.equal(output.signal.length,1024);
    }
    for(const frequency of [256,280,400,1000,4096]) {
      const output=ctx.resampleBrowserWav(tone(rate,frequency,2,Math.PI/2),rate);
      const residual=interiorPeak(output.signal);
      assert.ok(residual<10**(-70/20),`${rate} Hz / ${frequency} Hz residual ${residual}`);
    }
  }
});

test("400 Hz alias is suppressed while an in-band tone survives", () => {
  const {ctx}=setup(), rate=48000, input=tone(rate,400);
  const oldStep=Math.round(rate/512), legacy=Float64Array.from({length:Math.floor(input.length/oldStep)},(_,i)=>input[i*oldStep]);
  const oldRate=rate/oldStep, alias=oldRate-400;
  assert.ok(fittedAmplitude(legacy,alias,oldRate)>.99);
  const output=ctx.resampleBrowserWav(input,rate);
  assert.ok(interiorPeak(output.signal)<10**(-70/20));
  assert.ok(fittedAmplitude(ctx.resampleBrowserWav(tone(rate,100),rate).signal,100)>.99);
});

test("resampling preserves time, native samples and duration bounds without mutating input", () => {
  const {ctx}=setup();
  const native=Float64Array.from([1,-1,.125,0,3]), before=Array.from(native);
  const same=ctx.resampleBrowserWav(native,512);assert.deepEqual(Array.from(same.signal),before);assert.notEqual(same.signal,native);
  const impulse=new Float64Array(48000*2);impulse[24000]=1;
  const converted=ctx.resampleBrowserWav(impulse,48000);
  let max=0;for(let i=0;i<converted.signal.length;i++)if(converted.signal[i]>converted.signal[max])max=i;
  assert.equal(max,256);assert.equal(impulse[24000],1);assert.equal(impulse[0],0);
  const input=tone(44100,100,1.234), copy=Float32Array.from(input), out=ctx.resampleBrowserWav(input,44100);
  assert.deepEqual(Array.from(input),Array.from(copy));assert.equal(out.signal.length,Math.floor(input.length*512/44100));
  const loss=input.length/44100-out.signal.length/512;assert.ok(loss>=0 && loss<1/512);
  assert.equal(out.preprocessing.padding,"zero");assert.equal(out.preprocessing.group_delay_compensated,true);
  assert.ok(out.preprocessing.edge_half_support_s>0);assert.equal(out.preprocessing.kaiser_beta,8.6);
  assert.equal(out.preprocessing.processing_mode,"offline_centered_fir");assert.equal(out.preprocessing.future_data_used,true);
});

test("WAV rate/resource/container failures are explicit and decoder provenance is retained", async () => {
  const {ctx}=setup();
  for(const rate of [0,NaN,Infinity,511,48000.5]) assert.throws(()=>ctx.resampleBrowserWav([1,2,3],rate));
  assert.throws(()=>ctx.resampleBrowserWav([1,NaN],48000),/finite/);
  assert.throws(()=>ctx.resampleBrowserWav([1,2],48000),/short/);
  assert.throws(()=>ctx.resampleBrowserWav(new Float32Array(50000),48001),/resource/);
  assert.throws(()=>ctx.inspectWavContainer(new ArrayBuffer(10)),/RIFF/);
  await assert.rejects(()=>ctx.readAudioOrCsv(wavFile("stereo.wav",44100,[1,2,3,4],2),512),/mono/);
  await assert.rejects(()=>ctx.readAudioOrCsv(wavFile("low.wav",128,[1,2,3,4]),512),/at least/);
  let closed=0;
  ctx.window={AudioContext:class {async decodeAudioData(){return {numberOfChannels:1,sampleRate:48000,getChannelData:()=>tone(48000,100)};}async close(){closed++;}}};
  const output=await ctx.readAudioOrCsv(wavFile("different-container-rate.wav",44100,tone(44100,100)),512);
  assert.equal(output.fs,512);assert.equal(output.preprocessing.original_sample_rate_hz,44100);
  assert.equal(output.preprocessing.decoded_sample_rate_hz,48000);assert.equal(output.preprocessing.decoder_resampled,true);
  assert.equal(output.preprocessing.resampler_id,"bsd_kaiser_sinc_v1");assert.equal(output.preprocessing.cutoff_hz,224);
  assert.equal(output.preprocessing.wav_container.bits_per_sample,32);assert.equal(closed,1);
});


test("whole WAV workflow accepts different original rates and exports exact-rate filter provenance", async () => {
  const {ctx,elements:e}=setup();
  ctx.window={AudioContext:class {
    async decodeAudioData(buffer) {
      const view=new DataView(buffer),rate=view.getUint32(24,true),count=view.getUint32(40,true)/4;
      const samples=Float32Array.from({length:count},(_,i)=>view.getFloat32(44+i*4,true));
      return {numberOfChannels:1,sampleRate:rate,getChannelData:()=>samples};
    }
    async close() {}
  }};
  setValues(e,{workflow:"files",targetFap:.05,csvFs:512,winIdx:0});
  e.get("bgFile").files=[wavFile("background.wav",44100,tone(44100,100,40))];
  e.get("tgtFile").files=[wavFile("target.wav",48000,tone(48000,100,2))];
  await ctx.handleRunMeasurement();
  const payload=JSON.parse(vm.runInContext("JSON.stringify(buildExportPayload(LAST_RUN))",ctx));
  assert.equal(payload.metadata.sample_rate_hz,512);assert.equal(payload.metadata.n_bg_windows,20);
  assert.equal(payload.summary.n_target_windows,1);
  assert.equal(payload.metadata.independence_audit.status,"UNVERIFIED_EXTERNAL_CHECK_REQUIRED");
  for(const role of ["background","target"]) {
    const preprocessing=payload.metadata.input_provenance[role].preprocessing;
    assert.equal(preprocessing.sample_rate_hz,512);assert.equal(preprocessing.cutoff_hz,224);
    assert.equal(preprocessing.antialias_filter,"centered Kaiser-windowed sinc FIR");
    assert.match(preprocessing.warnings[0],/zero padding/);
  }
  assert.equal(payload.metadata.input_provenance.background.preprocessing.original_sample_rate_hz,44100);
  assert.equal(payload.metadata.input_provenance.target.preprocessing.original_sample_rate_hz,48000);
  assert.equal(e.get("dlJson").disabled,false);
});


test("RIFF metadata parser validates chunk boundaries and accepts padded ancillary chunks", async () => {
  const {ctx}=setup(), base=await wavFile("base.wav",44100,[1,2,3,4]).arrayBuffer();
  const truncated=base.slice(0,base.byteLength-1);assert.throws(()=>ctx.inspectWavContainer(truncated),/Truncated/);
  const shortFmt=base.slice(0);new DataView(shortFmt).setUint32(16,15,true);assert.throws(()=>ctx.inspectWavContainer(shortFmt),/format chunk/);
  const badData=base.slice(0);new DataView(badData).setUint32(40,100000,true);assert.throws(()=>ctx.inspectWavContainer(badData),/Truncated/);
  const duplicate=new ArrayBuffer(base.byteLength+24), dup=new Uint8Array(duplicate);
  dup.set(new Uint8Array(base,0,36),0);dup.set(new Uint8Array(base,12,24),36);dup.set(new Uint8Array(base,36),60);
  new DataView(duplicate).setUint32(4,duplicate.byteLength-8,true);assert.throws(()=>ctx.inspectWavContainer(duplicate),/Ambiguous/);
  const padded=new ArrayBuffer(base.byteLength+10), bytes=new Uint8Array(padded), view=new DataView(padded);
  bytes.set(new Uint8Array(base,0,12));bytes.set([74,85,78,75],12);view.setUint32(16,1,true);bytes[20]=7;bytes[21]=0;
  bytes.set(new Uint8Array(base,12),22);view.setUint32(4,padded.byteLength-8,true);
  const metadata=ctx.inspectWavContainer(padded);assert.equal(metadata.sample_rate_hz,44100);assert.equal(metadata.data_bytes,16);
});
