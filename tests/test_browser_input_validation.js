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

test("WAV failures close AudioContext and mono preprocessing retains explicit alias warning", async () => {
  const {ctx}=setup();let closes=0;
  ctx.window={AudioContext:class {async decodeAudioData(){return {numberOfChannels:2};}async close(){closes++;}}};
  await assert.rejects(() => ctx.readAudioOrCsv(file("stereo.wav","x"),512),/mono/);assert.equal(closes,1);
  ctx.window={AudioContext:class {async decodeAudioData(){return {numberOfChannels:1,sampleRate:512,getChannelData:()=>Float32Array.from([1,2,3,4])};}async close(){closes++;}}};
  const loaded=await ctx.readAudioOrCsv(file("mono.wav","x"),512);assert.equal(loaded.preprocessing.antialias_filter,"none");assert.match(loaded.preprocessing.warnings[0],/WAV_ALIASING_RISK/);assert.equal(closes,2);
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
