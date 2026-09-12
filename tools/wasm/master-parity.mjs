// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// P10 acceptance: the mastering ABI compiled to wasm renders what the native reference renders.
//
//   node master-parity.mjs build/fcmaster.node.js ../../build/tools/fcore_master 48000 2 short.f32
//
// It runs BOTH paths of the ABI — the streaming one (configure → process in 4096-frame blocks → flush →
// drop the latency) and the search (`fc_master_solve`) — because they fail differently: a mistake in the
// marshalling shows up in both, a mistake in the block loop or the latency drop only in the first, and a
// mistake in the request struct only in the second. The wasm side mirrors `abiRender` in fcore_master.cpp
// call for call, including its 4096-frame block size, so any difference is the TIER and not the harness.
//
// ==================================================================================================
// WHAT "AGREE" MEANS HERE, AND WHY IT IS NOT BIT-FOR-BIT
// ==================================================================================================
// `fcore_master` is the one target in tools/CMakeLists.txt that deliberately does NOT take
// -ffp-contract=off: it is the reference for the C++ API the DESKTOP app links, so it takes the
// repository's declared -ffp-contract=on (law 10). Baseline wasm has no scalar FMA, so this module
// cannot contract whatever it is told. The two are therefore the SAME PROGRAMME AT TWO ROUNDINGS, and
// the repository has already measured that pair on this very code: 203269 of 288000 samples differ,
// worst 2.03e-6 (tools/CMakeLists.txt, at the fcore_master target).
//
// So the acceptance is a TOLERANCE, and it is stated rather than discovered: 1e-5 in sample value,
// which is ~1.5 bits of a 16-bit master and about 5x the contraction difference already on record. A
// marshalling bug does not land inside it — a field written at the wrong offset moves a gain, a
// threshold or a ceiling, and those move samples by tenths, not by millionths. Bit-exactness IS
// available and IS tested, in the place where it means something: `fcore_master selftest` compares the
// ABI against a direct C++ call in ONE binary on ONE machine, where both sides round identically.
//
// The REPORTED NUMBERS are held tighter, and to a different tolerance for a different reason: latency,
// internalBlock and the pass count are integers and must match exactly, while the measured loudness
// figures are doubles computed over the whole programme and carry the same rounding difference
// amplified by the search's own feedback — 1e-3 dB, which is far below anything audible and far above
// what a swapped field would produce.

import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { createRequire } from 'node:module';
import { resolve, join } from 'node:path';
import { tmpdir } from 'node:os';
import { Struct, sizeOf, assertLayoutMatches, statusName, FC_SOLVE_STATUS, FC_CONSTRAINT }
    from './fc-master-layout.mjs';

const [, , modPath, nativePath, srArg, chArg, rawPath] = process.argv;
if (!rawPath) {
    console.error('usage: node master-parity.mjs <module.js> <fcore_master> <sampleRate> <channels> <raw.f32le>');
    process.exit(2);
}
const sr = Number(srArg), ch = Number(chArg);
const BLOCK = 4096;                       // fcore_master's own default, mirrored so the comparison is of tiers
const TOL_AUDIO = 1e-5;
const TOL_DB = 1e-3;

// ── the fixture, interleaved f32le, de-interleaved to the ABI's planar shape ───────────────────────
const raw = readFileSync(rawPath);
const inter = (raw.byteOffset % 4 === 0)
    ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
    : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
const frames = Math.floor(inter.length / ch);
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c)
    for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];
console.log(`fixture: ${frames} frames, ${ch} ch @ ${sr} Hz (${(frames / sr).toFixed(3)} s)`);

const require = createRequire(import.meta.url);
const M = await (require(resolve(modPath)))();
assertLayoutMatches(M);
console.log(`module: ABI v${M._fc_master_abi_version()}, params ${M._fc_master_sizeof_params()} B, `
          + `config ${M._fc_master_sizeof_config()} B, maxChannels ${M._fc_master_max_channels()}, `
          + `maxEqBands ${M._fc_master_max_eq_bands()}`);

// ── heap helpers. NO VIEW IS EVER HELD ACROSS A CALL — see the note in fc-master-layout.mjs ────────
const alloc = (bytes) => { const p = M._malloc(bytes); if (!p) throw new Error(`wasm OOM: ${bytes} B`); return p; };
const ok = (st, what) => { if (st !== 0) throw new Error(`${what}: ${statusName(st)}`); };
function writeF32 (ptr, arr) { M.HEAPF32.set(arr, ptr >>> 2); }
function readF32 (ptr, n) { return new Float32Array(M.HEAPF32.subarray(ptr >>> 2, (ptr >>> 2) + n)); }

// A fresh handle configured from the core's own defaults, which is what the CLI does.
function makeHandle () {
    const cfgP = alloc(sizeOf('fc_master_config'));
    M._fc_master_config_default(cfgP);
    const cfg = new Struct(M, 'fc_master_config', cfgP);
    cfg.set('sampleRate', sr).set('channels', ch);       // the two the defaults leave at zero, by design
    const hP = alloc(4);
    ok(M._fc_master_create(cfgP, hP), 'create');
    const h = new DataView(M.HEAPF32.buffer).getUint32(hP, true);
    M._free(hP); M._free(cfgP);
    return h;
}
function defaultParams () {
    const p = alloc(sizeOf('fc_master_params'));
    M._fc_master_params_default(p);
    return p;
}

// ── path 1: the streaming render, mirroring abiRender() ───────────────────────────────────────────
function wasmRender () {
    const h = makeHandle();
    const prmP = defaultParams();
    const resP = alloc(sizeOf('fc_master_resolved'));
    const res = new Struct(M, 'fc_master_resolved', resP).init();
    ok(M._fc_master_configure(h, prmP, resP), 'configure');
    const D = res.get('latencySamples');
    const geometry = {
        latency: D, internalBlock: res.get('internalBlock'),
        ceiling: res.get('limiterCeilingDbTp'), release: res.get('limiterReleaseMs'),
    };

    const stride = frames + D;
    const stream = new Float32Array(stride * ch);
    const inP = alloc(BLOCK * ch * 4), outP = alloc(BLOCK * ch * 4);
    let written = 0;
    for (let off = 0; off < frames; ) {
        const m = Math.min(BLOCK, frames - off);
        const slice = new Float32Array(m * ch);          // planar with stride m — the ABI's own rule
        for (let c = 0; c < ch; ++c) slice.set(planar.subarray(c * frames + off, c * frames + off + m), c * m);
        writeF32(inP, slice);
        ok(M._fc_master_process(h, inP, outP, m), 'process');
        const got = readF32(outP, m * ch);
        for (let c = 0; c < ch; ++c) stream.set(got.subarray(c * m, c * m + m), c * stride + written);
        written += m; off += m;
    }
    if (D > 0) {
        const tailP = alloc(D * ch * 4), wP = alloc(4);
        ok(M._fc_master_flush(h, tailP, D, wP), 'flush');
        const got = new DataView(M.HEAPF32.buffer).getUint32(wP, true);
        const tail = readF32(tailP, D * ch);
        for (let c = 0; c < ch; ++c) stream.set(tail.subarray(c * D, c * D + got), c * stride + written);
        written += got;
        M._free(tailP); M._free(wP);
    }
    // Drop the latency, exactly as the CLI does.
    const out = new Float32Array(frames * ch);
    for (let c = 0; c < ch; ++c)
        for (let n = 0; n < frames; ++n) out[c * frames + n] = stream[c * stride + n + D];

    const statsP = alloc(sizeOf('fc_master_stats'));
    const stats = new Struct(M, 'fc_master_stats', statsP).init();
    ok(M._fc_master_get_stats(h, statsP), 'get_stats');
    const seen = { framesIn: stats.get('framesIn'), flushed: stats.get('framesFlushed'),
                   nonFinite: stats.get('nonFiniteIn') };
    M._free(statsP); M._free(inP); M._free(outP); M._free(resP); M._free(prmP);
    ok(M._fc_master_destroy(h), 'destroy');
    return { out, geometry, seen };
}

// ── path 2: the loudness search ───────────────────────────────────────────────────────────────────
function wasmSolve (targetLufs, tpDb) {
    const h = makeHandle();
    const prmP = defaultParams();
    const reqP = alloc(sizeOf('fc_loudness_request'));
    M._fc_loudness_request_default(reqP);
    new Struct(M, 'fc_loudness_request', reqP).set('targetLufs', targetLufs).set('maxTruePeakDbTp', tpDb);

    // The budget, before the call that cannot report an exhausted heap as a status (law 11d).
    const needP = alloc(sizeOf('fc_need'));
    const need = new Struct(M, 'fc_need', needP).init();
    ok(M._fc_master_need(h, 0 /* FC_NEED_SOLVE */, frames, needP), 'need');
    const budget = {
        call: Number(need.get('callBytes')), prepare: Number(need.get('solverPrepareBytes')),
        facade: Number(need.get('facadeBytes')), prepared: need.get('solverPrepared'),
    };
    M._free(needP);

    // in and out MUST be distinct buffers: a search reads its input again on every pass.
    const inP = alloc(frames * ch * 4), outP = alloc(frames * ch * 4), solP = alloc(4);
    writeF32(inP, planar);
    ok(M._fc_master_solve(h, prmP, reqP, inP, outP, frames, solP), 'solve');
    const sol = new DataView(M.HEAPF32.buffer).getUint32(solP, true);

    const sumP = alloc(sizeOf('fc_solution_summary'));
    const sum = new Struct(M, 'fc_solution_summary', sumP).init();
    ok(M._fc_solution_summary_get(sol, sumP), 'summary');
    const measP = alloc(sizeOf('fc_measurement'));
    const meas = new Struct(M, 'fc_measurement', measP).init();
    ok(M._fc_solution_measurement(sol, measP), 'measurement');

    const verdict = {
        status: sum.get('status'), binding: sum.get('binding'),
        gain: sum.get('preLimiterGainDb'), ceiling: sum.get('ceilingDbTp'), passes: sum.get('passes'),
        I: meas.get('integratedLufs'), TP: meas.get('truePeakDbTp'),
        LRA: meas.get('loudnessRangeLu'), PLR: meas.get('plrDb'),
        loudnessValid: meas.get('loudnessValid'), lraValid: meas.get('lraValid'),
    };
    const out = readF32(outP, frames * ch);

    ok(M._fc_solution_destroy(sol), 'solution_destroy');
    ok(M._fc_master_destroy(h), 'destroy');
    for (const p of [sumP, measP, inP, outP, solP, reqP, prmP]) M._free(p);
    return { out, verdict, budget };
}

// ── the native side ───────────────────────────────────────────────────────────────────────────────
const tmp = mkdtempSync(join(tmpdir(), 'fcmaster-parity-'));
// spawnSync rather than execFileSync: fcore_master prints its verdict on stdout and its RESOLVED
// GEOMETRY on stderr, and both are part of what is being compared.
const nativeRun = (args, outFile) => {
    const r = spawnSync(resolve(nativePath), args, { encoding: 'utf8' });
    if (r.status !== 0) throw new Error(`${args[0]}: exit ${r.status}\n${r.stderr}`);
    const bytes = readFileSync(outFile);
    const f = (bytes.byteOffset % 4 === 0)
        ? new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4)
        : new Float32Array(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength));
    const pl = new Float32Array(frames * ch);           // interleaved on disk, planar here
    for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) pl[c * frames + i] = f[i * ch + c];
    return { stdout: r.stdout, stderr: r.stderr, audio: pl };
};

// ── comparison ────────────────────────────────────────────────────────────────────────────────────
let failures = 0;
const check = (pass, what, detail = '') => {
    console.log(`  [${pass ? 'ok' : 'FAIL'}] ${what}${detail ? ' — ' + detail : ''}`);
    if (!pass) ++failures;
};
function compareAudio (a, b, what, tol = TOL_AUDIO) {
    let worst = 0, at = -1, differ = 0, peak = 0;
    for (let i = 0; i < a.length; ++i) {
        const d = Math.abs(a[i] - b[i]);
        if (d !== 0) ++differ;
        if (d > worst) { worst = d; at = i; }
        if (Math.abs(b[i]) > peak) peak = Math.abs(b[i]);
    }
    check(worst <= tol, what,
          `worst ${worst.toExponential(3)} at sample ${at}, ${differ}/${a.length} differ at all, `
        + `native peaks at ${peak.toFixed(4)} (tolerance ${tol.toExponential(0)})`);
    return { worst, differ, peak };
}
const near = (a, b, tol, what) =>
    check(Math.abs(a - b) <= tol, what, `wasm ${a} vs native ${b}`);

console.log('\n=== 1. streaming render: configure → process(4096) → flush → drop latency');
const R = wasmRender();
const rNative = nativeRun(['render', String(sr), String(ch), resolve(rawPath), join(tmp, 'r.f32')], join(tmp, 'r.f32'));
const nativeGeo = /latency=(\d+) internalBlock=(\d+) ceiling=(\S+) release=(\S+)/.exec(rNative.stderr);
check(!!nativeGeo, 'the native render reported its resolved geometry');
if (nativeGeo) {
    near(R.geometry.latency, Number(nativeGeo[1]), 0, 'latencySamples');
    near(R.geometry.internalBlock, Number(nativeGeo[2]), 0, 'internalBlock');
    near(R.geometry.ceiling, Number(nativeGeo[3]), 0, 'limiter ceiling dBTP');
    near(R.geometry.release, Number(nativeGeo[4]), 0, 'limiter release ms');
}
check(Number(R.seen.framesIn) === frames, 'stats: framesIn', `${R.seen.framesIn} of ${frames}`);
check(Number(R.seen.nonFinite) === 0, 'stats: no non-finite input', `${R.seen.nonFinite}`);
const rd = compareAudio(R.out, rNative.audio, 'rendered audio agrees within tolerance');
check(rd.peak > 1e-3, 'the render is not silence', `native peaks at ${rd.peak.toFixed(4)}`);

console.log('\n=== 2. the loudness search: target -14 LUFS, true peak -1 dBTP');
const S = wasmSolve(-14, -1);
console.log(`  budget: call ${S.budget.call} B, solver prepare ${S.budget.prepare} B, `
          + `facade ${S.budget.facade} B, prepared=${S.budget.prepared}`);
const sNative = nativeRun(['solve', String(sr), String(ch), resolve(rawPath), join(tmp, 's.f32'),
                           'target=-14', 'tp=-1'], join(tmp, 's.f32'));
const l1 = /status=(\S+) binding=(\S+) gain=(\S+) ceiling=(\S+) passes=(\S+)/.exec(sNative.stdout);
const l2 = /I=(\S+) TP=(\S+) LRA=(\S+) PLR=(\S+) loudnessValid=(\S+) lraValid=(\S+)/.exec(sNative.stdout);
check(!!l1 && !!l2, 'the native solve reported a verdict');
console.log(`  wasm:   status=${FC_SOLVE_STATUS[S.verdict.status]} binding=${FC_CONSTRAINT[S.verdict.binding]} `
          + `gain=${S.verdict.gain} ceiling=${S.verdict.ceiling} passes=${S.verdict.passes}`);
console.log(`          I=${S.verdict.I} TP=${S.verdict.TP} LRA=${S.verdict.LRA} PLR=${S.verdict.PLR}`);
if (l1 && l2) {
    near(S.verdict.status, Number(l1[1]), 0, 'solve status');
    near(S.verdict.binding, Number(l1[2]), 0, 'binding constraint');
    near(S.verdict.passes, Number(l1[5]), 0, 'pass count');
    near(S.verdict.gain, Number(l1[3]), TOL_DB, 'pre-limiter gain dB');
    near(S.verdict.ceiling, Number(l1[4]), TOL_DB, 'ceiling dBTP');
    near(S.verdict.I, Number(l2[1]), TOL_DB, 'integrated LUFS');
    near(S.verdict.TP, Number(l2[2]), TOL_DB, 'true peak dBTP');
    near(S.verdict.PLR, Number(l2[4]), TOL_DB, 'PLR dB');
    check(S.verdict.loudnessValid === Number(l2[5]), 'loudnessValid');
}
check(Math.abs(S.verdict.I - (-14)) <= 1.0, 'the search actually landed near the target',
      `I=${S.verdict.I.toFixed(3)} LUFS for target -14`);
compareAudio(S.out, sNative.audio, 'delivered master agrees within tolerance');

console.log('\n=== 3. LRA, the one measurement that crosses on its own');
{
    const h = makeHandle();
    const inP = alloc(frames * ch * 4), dP = alloc(8);
    writeF32(inP, planar);
    ok(M._fc_master_measure_lra(h, inP, frames, dP), 'measure_lra');
    const lra = new DataView(M.HEAPF32.buffer).getFloat64(dP, true);
    M._free(inP); M._free(dP);
    ok(M._fc_master_destroy(h), 'destroy');
    const n = spawnSync(resolve(nativePath), ['lra', String(sr), String(ch), resolve(rawPath)],
                        { encoding: 'utf8' });
    if (n.status !== 0) throw new Error(`lra: exit ${n.status}\n${n.stderr}`);
    const nat = Number(n.stdout.trim());
    near(lra, nat, TOL_DB, 'loudness range LU');
}

rmSync(tmp, { recursive: true, force: true });
console.log(`\n${failures === 0 ? 'PASS' : 'FAIL'} — ${failures} failure(s)`);
process.exit(failures === 0 ? 0 : 1);
