"use strict";
// Test-only adapter: run the production inline resampler on stdin JSON; no DOM or files are written.
const fs=require("node:fs"),path=require("node:path"),vm=require("node:vm");
const html=fs.readFileSync(path.join(__dirname,"../docs/workstation/index.html"),"utf8");
const script=html.match(/<script>([\s\S]*?)<\/script>/)[1].split("// --- Event Listeners ---")[0];
const ctx=vm.createContext({});vm.runInContext(script,ctx);
const input=JSON.parse(fs.readFileSync(0,"utf8"));
const result=ctx.resampleBrowserWav(Float64Array.from(input.signal),input.source_rate,512);
process.stdout.write(JSON.stringify({signal:Array.from(result.signal),fs:result.fs,preprocessing:result.preprocessing}));
