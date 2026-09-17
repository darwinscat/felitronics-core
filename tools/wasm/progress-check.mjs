// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';
import { performance } from 'node:perf_hooks';
import { Struct, sizeOf, layoutOf, assertLayoutMatches, statusName } from './fc-master-layout.mjs';

const [, , modPath, srArg, chArg, rawPath, repeatsArg] = process.argv;
if (!rawPath) {
    console.error('usage: node progress-check.mjs <fcmaster.node.js> <sampleRate> <channels> <raw.f32le> [repeats]');
    process.exit(2);
}
const sr = Number(srArg), ch = Number(chArg), REPEATS = Number(repeatsArg ?? 3);
const FC_OK = 0, FC_ERR_STATE = 10, FC_ERR_CANCELLED = 15;

const raw = readFileSync(rawPath);
const inter = new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength - raw.byteLength % 4));
const frames = Math.floor(inter.length / ch);
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];
console.log(`programme: ${frames} frames, ${ch} ch @ ${sr} Hz (${(frames / sr).toFixed(2)} s)`);

const M = await (createRequire(import.meta.url)(resolve(modPath)))();
assertLayoutMatches(M);

let failures = 0;
const check = (pass, what, detail = '') => {
    console.log(`  [${pass ? 'ok' : 'FAIL'}] ${what}${detail ? ' — ' + detail : ''}`);
    if (!pass) ++failures;
};
const alloc = bytes => { const p = M._malloc(bytes); if (!p) throw new Error(`wasm OOM: ${bytes} B`); return p; };
const ok = (st, what) => { if (st !== FC_OK) throw new Error(`${what}: ${statusName(st)}`); };
const dv = () => new DataView(M.HEAPF32.buffer);
const bytesAt = (p, n) => Buffer.from(new Uint8Array(M.HEAPF32.buffer, p, n));
const same = (a, b) => a.length === b.length && Buffer.compare(a, b) === 0;
const SCALAR_BYTES = { i32: 4, u32: 4, f32: 4, f64: 8, u64: 8 };
function fieldBytes (name, p) {
    const parts = [];
    const walk = (n, base) => {
        for (const [, f] of layoutOf(n).fields)
            for (let i = 0; i < (f.count || 1); ++i) {
                const at = base + f.offset + i * f.stride;
                if (SCALAR_BYTES[f.type] !== undefined) parts.push(bytesAt(at, SCALAR_BYTES[f.type]));
                else walk(f.type, at);
            }
    };
    walk(name, p);
    return Buffer.concat(parts);
}
const differing = (a, b) => ['audio', 'summary', 'measurement', 'log', 'traces'].filter(k => !same(a.bytes[k], b.bytes[k]));

const live = [];
function destroyHandles () { while (live.length) ok(M._fc_master_destroy(live.pop()), 'destroy'); }
function makeHandle (deliveryRate = 0) {
    const cfgP = alloc(sizeOf('fc_master_config'));
    const cfg = new Struct(M, 'fc_master_config', cfgP).init();
    ok(M._fc_master_config_defaults(cfgP), 'config_defaults');
    cfg.set('sampleRate', sr).set('channels', ch).set('deliveryRate', deliveryRate);
    const hP = alloc(4);
    ok(M._fc_master_create(cfgP, hP), 'create');
    const h = dv().getUint32(hP, true);
    M._free(hP); M._free(cfgP);
    live.push(h);
    return h;
}
const prmP = alloc(sizeOf('fc_master_params'));
new Struct(M, 'fc_master_params', prmP).init();
ok(M._fc_master_params_defaults(prmP), 'params_defaults');
const reqP = alloc(sizeOf('fc_loudness_request'));
new Struct(M, 'fc_loudness_request', reqP).init();
ok(M._fc_loudness_request_defaults(reqP), 'request_defaults');
new Struct(M, 'fc_loudness_request', reqP).set('targetLufs', -14).set('maxTruePeakDbTp', -1);

function latencyOf (h) {
    const p = alloc(4);
    ok(M._fc_master_latency(h, p), 'latency');
    const d = dv().getInt32(p, true);
    M._free(p);
    return d;
}

const inP = alloc(frames * ch * 4), outP = alloc(frames * ch * 4);
M.HEAPF32.set(planar, inP >>> 2);

function withHandler (handler, call) {
    const msgs = [];
    let inCall = false;
    if (handler !== null)
        M.onProgress = msg => { msgs.push({ ...msg, t: performance.now(), inCall }); return handler(msg, msgs.length - 1); };
    else delete M.onProgress;
    inCall = true;
    const t0 = performance.now();
    const st = call();
    const t1 = performance.now();
    inCall = false;
    delete M.onProgress;
    return { st, msgs, t0, t1 };
}

function solveOn (h, handler, inputP = inP, inFrames = frames, outputP = outP, outFrames = frames, delivered = false) {
    const solP = alloc(4);
    dv().setUint32(solP, 0xDEADBEEF, true);
    const r = withHandler(handler, () => delivered
        ? M._fc_master_solve_delivered(h, prmP, reqP, inputP, inFrames, outputP, outFrames, solP)
        : M._fc_master_solve(h, prmP, reqP, inputP, outputP, inFrames, solP));
    r.solHandleSlot = dv().getUint32(solP, true);
    M._free(solP);
    if (r.st !== FC_OK) return r;
    const sol = r.solHandleSlot;
    const sumP = alloc(sizeOf('fc_solution_summary')), measP = alloc(sizeOf('fc_measurement'));
    const sum = new Struct(M, 'fc_solution_summary', sumP).init();
    new Struct(M, 'fc_measurement', measP).init();
    ok(M._fc_solution_summary_get(sol, sumP), 'summary');
    ok(M._fc_solution_measurement(sol, measP), 'measurement');
    const stride = sizeOf('fc_solve_pass'), logP = alloc(32 * stride), wP = alloc(4);
    ok(M._fc_solution_log(sol, logP, 32, wP), 'log');
    const written = dv().getUint32(wP, true);
    const traceParts = [];
    for (const stage of [0, 1]) {
        const bStride = sizeOf('fc_gr_trace_bucket'), bP = alloc(1000 * bStride), bwP = alloc(4);
        ok(M._fc_solution_gr_trace(sol, stage, bP, 1000, bwP), 'gr_trace');
        const bw = dv().getUint32(bwP, true);
        for (let i = 0; i < bw; ++i) traceParts.push(fieldBytes('fc_gr_trace_bucket', bP + i * bStride));
        M._free(bP); M._free(bwP);
    }
    const log = [];
    for (let i = 0; i < written; ++i) {
        const s = new Struct(M, 'fc_solve_pass', logP + i * stride);
        log.push(Object.fromEntries(['gainDb', 'ceilingDb', 'integratedLufs', 'truePeakDbTp', 'plrDb',
                                     'limiterMaxGrDb', 'loudnessRangeLu', 'violated'].map(k => [k, s.get(k)])));
    }
    r.status = sum.get('status'); r.passes = sum.get('passes');
    const logParts = [];
    for (let i = 0; i < written; ++i) logParts.push(fieldBytes('fc_solve_pass', logP + i * stride));
    r.bytes = { audio: bytesAt(outputP, outFrames * ch * 4), summary: fieldBytes('fc_solution_summary', sumP),
                measurement: fieldBytes('fc_measurement', measP), log: Buffer.concat(logParts),
                traces: Buffer.concat(traceParts) };
    r.log = log;
    ok(M._fc_solution_destroy(sol), 'solution_destroy');
    for (const p of [sumP, measP, logP, wP]) M._free(p);
    return r;
}

function lraOn (h, handler, inputP = inP, inFrames = frames) {
    const outD = alloc(8);
    dv().setFloat64(outD, -12345.5, true);
    const r = withHandler(handler, () => M._fc_master_measure_lra(h, inputP, inFrames, outD));
    r.lraBytes = bytesAt(outD, 8);
    r.lra = dv().getFloat64(outD, true);
    M._free(outD);
    return r;
}

function checkMessages (label, r, framesOf, unitsOf) {
    const { msgs } = r;
    check(msgs.length > 0 && msgs.every(m => m.inCall && m.t >= r.t0 && m.t <= r.t1),
          `${label}: every message arrived inside the call`,
          `${msgs.length} messages, first at +${msgs.length ? (msgs[0].t - r.t0).toFixed(1) : '?'} ms, `
          + `last ${msgs.length ? (r.t1 - msgs[msgs.length - 1].t).toFixed(1) : '?'} ms before it returned, call ${(r.t1 - r.t0).toFixed(0)} ms`);
    const stages = [];
    for (const m of msgs) {
        if (m.fraction === 0) stages.push({ stage: m.stage, pass: m.pass, maxPasses: m.maxPasses, msgs: [m] });
        else if (stages.length) stages[stages.length - 1].msgs.push(m);
    }
    let shape = 0, gaps = 0, worstGap = 0, records = 0;
    for (const s of stages) {
        const f = s.msgs.map(m => m.fraction);
        const last = s.msgs[s.msgs.length - 1];
        if (f[f.length - 1] !== 1) ++shape;
        for (let i = 1; i < f.length; ++i) if (!(f[i] > f[i - 1])) ++shape;
        if (!s.msgs.every(m => m.stage === s.stage && m.pass === s.pass && m.maxPasses === s.maxPasses)) ++shape;
        const rendered = s.stage === 'pass' || s.stage === 'final';
        if (s.msgs.some(m => m !== last && m.record !== undefined)) ++shape;
        if (rendered !== (last.record !== undefined)) ++shape;
        if (rendered) ++records;
        const bound = framesOf(s) / 100, units = unitsOf(s);
        for (let i = 1; i < f.length; ++i) {
            const g = (f[i] - f[i - 1]) * units;
            worstGap = Math.max(worstGap, g / framesOf(s));
            if (g > bound + 1e-6 * units) ++gaps;
        }
    }
    check(stages.length > 0 && shape === 0,
          `${label}: every stage opens at exactly 0, rises, closes at exactly 1, and only a render's close carries a record`,
          stages.map(s => `${s.stage}${s.pass ? ` ${s.pass}/${s.maxPasses}` : ''} (${s.msgs.length})`).join(', '));
    check(gaps === 0, `${label}: no two messages of a stage more than 1 % of the programme apart`,
          `widest ${(100 * worstGap).toFixed(3)} % of its frames`);
    return { stages, records };
}

const recordFields = ['gainDb', 'ceilingDb', 'integratedLufs', 'truePeakDbTp', 'plrDb', 'limiterMaxGrDb', 'loudnessRangeLu', 'violated'];
function checkRecords (label, r, stages) {
    const recs = stages.filter(s => s.stage === 'pass' || s.stage === 'final').map(s => s.msgs[s.msgs.length - 1].record);
    const passNumbers = stages.filter(s => s.stage === 'pass' || s.stage === 'final').map(s => s.pass);
    check(recs.length === r.log.length && recs.every((rec, i) => recordFields.every(k => Object.is(rec[k], r.log[i][k]))),
          `${label}: each render's record is its log entry, field for field`, `${recs.length} records, log ${r.log.length}`);
    check(passNumbers.every((p, i) => p === i + 1) && passNumbers.length === r.passes,
          `${label}: renders numbered 1..passes`, `passes ${r.passes}: ${passNumbers.join(',')}`);
    const maxPasses = new Struct(M, 'fc_loudness_request', reqP).get('maxPasses');
    check(stages.every(s => s.stage === 'pass' ? s.maxPasses === maxPasses
                          : s.stage === 'final' ? s.maxPasses === maxPasses + 1 : s.pass === 0 && s.maxPasses === 0),
          `${label}: a search render's bound is maxPasses (${maxPasses}), the final render's one more, other stages 0/0`);
}

function middleOf (msgs) {
    for (let d = 0; d < msgs.length; ++d)
        for (const i of [Math.floor(msgs.length / 2) + d, Math.floor(msgs.length / 2) - d])
            if (i >= 0 && i < msgs.length && msgs[i].fraction > 0 && msgs[i].fraction < 1) return i;
    return -1;
}

function stopAt (index) {
    const seen = { stopT: 0, after: 0 };
    const handler = (msg, i) => {
        if (i === index) { seen.stopT = performance.now(); return false; }
        if (i > index) ++seen.after;
        return true;
    };
    return { handler, seen };
}

function stageDuration (msgs, index) {
    let a = index; while (a > 0 && msgs[a].fraction !== 0) --a;
    let b = index; while (b < msgs.length - 1 && msgs[b].fraction !== 1) ++b;
    return msgs[b].t - msgs[a].t;
}

const D = latencyOf(makeHandle());
destroyHandles();
const lraStage = { framesOf: () => frames, unitsOf: () => frames };
const solveStage = { framesOf: () => frames, unitsOf: s => (s.stage === 'pass' || s.stage === 'final') ? 2 * frames + D : frames };

console.log('\nfc_master_measure_lra');
{
    const h = makeHandle();
    const plain = lraOn(h, null);
    ok(plain.st, 'lra (no handler)');
    const heard = lraOn(h, () => undefined);
    check(heard.st === FC_OK && same(plain.lraBytes, heard.lraBytes), 'the range is the same bits with a handler',
          `${plain.lra} LU, ${heard.lra} LU`);
    checkMessages('lra', heard, lraStage.framesOf, lraStage.unitsOf);

    const at = middleOf(heard.msgs), stop = stopAt(at);
    const stopped = lraOn(h, stop.handler);
    const back = stopped.t1 - stop.seen.stopT, stageMs = stageDuration(heard.msgs, at);
    check(stopped.st === FC_ERR_CANCELLED, 'a handler returning false mid-stage answers FC_ERR_CANCELLED',
          `${statusName(stopped.st)} at message ${at} of ${heard.msgs.length} (fraction ${heard.msgs[at]?.fraction.toFixed(3)})`);
    check(stop.seen.after === 0 && stopped.lra === -12345.5, 'nothing sent after the stop, *out untouched');
    check(back <= 0.01 * stageMs, 'back within 1 % of the stage', `${back.toFixed(2)} ms of a ${stageMs.toFixed(0)} ms stage`);
    const again = lraOn(h, null);
    check(again.st === FC_OK && same(again.lraBytes, plain.lraBytes), 'the same handle then measures the same bits');

    const thrown = lraOn(h, (m, i) => { if (i === 2) throw new Error('thrown by the handler'); });
    check(thrown.st === FC_ERR_CANCELLED && thrown.msgs.length === 3 && thrown.lra === -12345.5,
          'a handler that throws stops the call: FC_ERR_CANCELLED, nothing after, *out untouched',
          `${statusName(thrown.st)}, ${thrown.msgs.length} messages`);
    const after = lraOn(h, null);
    check(after.st === FC_OK && same(after.lraBytes, plain.lraBytes), 'and the handle then measures the same bits');
    destroyHandles();
}

console.log('\nfc_master_solve');
{
    const plain = solveOn(makeHandle(), null);
    ok(plain.st, 'solve (no handler)');
    const heard = solveOn(makeHandle(), () => undefined);
    check(heard.st === FC_OK, 'solved with a handler', `status ${plain.status}, passes ${plain.passes}`);
    check(differing(plain, heard).length === 0,
          'audio, summary, measurement, log and both traces are the same bits with a handler as without',
          `${plain.bytes.audio.length} audio bytes; differing: [${differing(plain, heard)}]`);
    const { stages } = checkMessages('solve', heard, solveStage.framesOf, solveStage.unitsOf);
    checkRecords('solve', heard, stages);

    const at = middleOf(heard.msgs), stop = stopAt(at);
    const h = makeHandle();
    const stopped = solveOn(h, stop.handler);
    const back = stopped.t1 - stop.seen.stopT, stageMs = stageDuration(heard.msgs, at);
    const m = heard.msgs[at];
    check(stopped.st === FC_ERR_CANCELLED, 'a handler returning false mid-search answers FC_ERR_CANCELLED',
          `${statusName(stopped.st)} at message ${at} of ${heard.msgs.length} — ${m.stage} ${m.pass}/${m.maxPasses} at ${m.fraction.toFixed(3)}`);
    check(stop.seen.after === 0 && stopped.solHandleSlot === 0xDEADBEEF, 'nothing sent after the stop, *out_solution untouched');
    check(back <= 0.01 * stageMs, 'back within 1 % of the stage', `${back.toFixed(2)} ms of a ${stageMs.toFixed(0)} ms stage`);
    const oneP = alloc(ch * 4);
    check(M._fc_master_process(h, inP, oneP, 1) === FC_ERR_STATE,
          'the stopped handle refuses to stream until a configure, as after any solve');
    M._free(oneP);
    const again = solveOn(h, null);
    check(again.st === FC_OK && differing(again, plain).length === 0,
          'the same handle then solves to the same bits', `differing: [${again.st === FC_OK ? differing(again, plain) : '-'}]`);
    destroyHandles();
}

const dRate = sr === 44100 ? 48000 : 44100;
console.log(`\ndelivering ${sr} -> ${dRate} (the first 60 s)`);
{
    const n = Math.min(frames, 60 * sr);
    const sliceP = alloc(n * ch * 4);
    for (let c = 0; c < ch; ++c) M.HEAPF32.set(planar.subarray(c * frames, c * frames + n), (sliceP >>> 2) + c * n);
    const h = makeHandle(dRate);
    const dP = alloc(4);
    ok(M._fc_master_delivered_frames(h, n, dP), 'delivered_frames');
    const d = dv().getUint32(dP, true);
    M._free(dP);
    const dOutP = alloc(d * ch * 4);
    const DD = latencyOf(h);
    const framesOf = s => s.stage === 'convert' ? n : d;
    const unitsOf = s => s.stage === 'convert' ? n : (s.stage === 'lra' ? d : 2 * d + DD);

    const lp = lraOn(h, null, sliceP, n), lh = lraOn(h, () => undefined, sliceP, n);
    check(lp.st === FC_OK && lh.st === FC_OK && same(lp.lraBytes, lh.lraBytes), 'lra: the same bits with a handler');
    const ls = checkMessages('lra', lh, framesOf, unitsOf).stages;
    check(ls.map(s => s.stage).join() === 'convert,lra', 'lra: a convert stage, then the range', ls.map(s => s.stage).join());

    const sp = solveOn(makeHandle(dRate), null, sliceP, n, dOutP, d, true);
    const sh = solveOn(makeHandle(dRate), () => undefined, sliceP, n, dOutP, d, true);
    check(sp.st === FC_OK && sh.st === FC_OK && differing(sp, sh).length === 0,
          'solve_delivered: the same bits with a handler', `differing: [${differing(sp, sh)}]`);
    const ss = checkMessages('solve_delivered', sh, framesOf, unitsOf).stages;
    check(ss[0]?.stage === 'convert', 'solve_delivered: the conversion is the first stage');
    checkRecords('solve_delivered', sh, ss);

    const inConvert = sh.msgs.findIndex(m => m.stage === 'convert' && m.fraction >= 0.5);
    const stop = stopAt(inConvert);
    const h2 = makeHandle(dRate);
    const st = solveOn(h2, stop.handler, sliceP, n, dOutP, d, true);
    check(st.st === FC_ERR_CANCELLED && stop.seen.after === 0 && st.solHandleSlot === 0xDEADBEEF,
          'a stop while converting is FC_ERR_CANCELLED, with nothing after it', statusName(st.st));
    const configured = alloc(d * ch * 4);
    check(M._fc_master_render_delivered(h2, sliceP, n, configured, d) === FC_OK,
          'a stop before the first render moved nothing: the handle still renders without a configure');
    M._free(configured);
    const again = solveOn(h2, null, sliceP, n, dOutP, d, true);
    check(again.st === FC_OK && differing(again, sp).length === 0, 'and then solves to the same bits',
          `differing: [${again.st === FC_OK ? differing(again, sp) : '-'}]`);
    M._free(dOutP); M._free(sliceP);
    destroyHandles();
}

console.log(`\ncost: ${REPEATS} alternated runs each, minimum compared`);
{
    let sink = 0;
    const handler = msg => { sink += msg.fraction; };
    const times = { lra: { off: [], on: [] }, solve: { off: [], on: [] } };
    for (let r = 0; r < REPEATS; ++r)
        for (const on of (r % 2 === 0 ? [false, true] : [true, false])) {
            const h = makeHandle();
            const lr = lraOn(h, on ? handler : null);
            ok(lr.st, 'lra'); times.lra[on ? 'on' : 'off'].push(lr.t1 - lr.t0);
            const so = solveOn(h, on ? handler : null);
            ok(so.st, 'solve'); times.solve[on ? 'on' : 'off'].push(so.t1 - so.t0);
            destroyHandles();
        }
    for (const k of ['lra', 'solve']) {
        const off = Math.min(...times[k].off), on = Math.min(...times[k].on);
        const rel = (on - off) / off;
        check(rel <= 0.01, `${k}: with a handler within 1 % of without`,
              `min ${off.toFixed(0)} ms without, ${on.toFixed(0)} ms with (${(100 * rel).toFixed(2)} %); `
              + `all without [${times[k].off.map(x => x.toFixed(0))}], with [${times[k].on.map(x => x.toFixed(0))}]`);
    }
}

console.log(failures === 0 ? '\nALL OK' : `\n${failures} FAILED`);
process.exit(failures === 0 ? 0 : 1);
