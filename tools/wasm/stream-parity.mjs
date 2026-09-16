// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the streaming-surface parity check (P120). Emits EXACTLY what `fcore_measure stream`
// emits, so the test is a diff:
//
//   fcore_measure stream 88200 2 x.f32 --chunk 4096                              > native.txt
//   node stream-parity.mjs build/fcprobe.node.js stream 88200 2 x.f32 --chunk 4096 > wasm.txt
//   diff native.txt wasm.txt
//
// Usage: node stream-parity.mjs <module.js> stream <sampleRate> <channels> <raw.f32le> [--chunk N]
//
// Both sides cut the file into pieces of exactly N frames (the last one shorter) and read between EVERY piece:
// the loudness readings as bit patterns, and the runs decided so far. Here that reading goes through the
// fc_stream_* handles — the planar copy, the packing, `from` — and there through fcore::StreamProbe directly,
// so a diff covers the C ABI as well as the tier. The runs are drained through a buffer of a few runs, so
// `from` advances in short steps and a run is read whole or not at all.
//
// The detached-view rule of parity.mjs holds: no HEAP view is kept across a call into the module.

import { bitsOf } from './blocks-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, mode, srArg, chArg, rawPath, ...rest] = process.argv;
if (!rawPath || mode !== 'stream') {
    console.error('usage: node stream-parity.mjs <module.js> stream <sampleRate> <channels> <raw.f32le> [--chunk N]');
    process.exit(2);
}
// The native tool's option and argument rules, as clips-parity.mjs mirrors them — see the notes there.
const refuse = msg => { console.error(msg); process.exit(2); };
const given = {};
for (let i = 0; i < rest.length; i++) {
    if (rest[i] === '--precise') continue;
    if (rest[i] !== '--chunk' || rest[i] in given || i + 1 >= rest.length) refuse(`bad, unknown or repeated option: ${rest[i]}`);
    given[rest[i]] = rest[++i];
}
const count = (s, name) => { if (!/^[0-9]{1,18}$/.test(s)) refuse(`bad ${name}: ${s}`); return Number(s); };
const rate = s => {
    if (!/^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/.test(s)) refuse('bad sampleRate');
    const v = Number(s);
    if (!Number.isFinite(v) || !(v > 0)) refuse('bad sampleRate');
    return v;
};
const sr = rate(srArg), ch = count(chArg, 'channels');
if (ch < 1 || ch > 16) refuse('bad sampleRate/channels');            // core::kMaxChannels
const chunk = '--chunk' in given ? count(given['--chunk'], '--chunk') : 4096;
if (chunk < 1 || chunk > 0x7FFFFFFF) refuse('--chunk out of range');

const raw = readFileSync(rawPath);
if (raw.byteLength % (4 * ch) !== 0) refuse(`the file is not a whole number of ${ch}-channel float32 frames`);
const inter = new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
const frames = inter.length / ch;
if (frames === 0) refuse('the file is empty');
const planes = [];
for (let c = 0; c < ch; ++c) {
    const p = new Float32Array(frames);
    for (let i = 0; i < frames; ++i) p[i] = inter[i * ch + c];
    planes.push(p);
}

const require = createRequire(import.meta.url);
const M = await require(resolve(modPath))();

const h = M._fc_stream_create(sr, ch);
if (h === 0) refuse('stream.prepare refused');
const fields = M._fc_stream_loudness_fields();
const stride = M._fc_stream_clips_stride();
if (fields !== 5 || stride !== 6) { console.error(`stream v1 is 5 loudness fields and 6 doubles a run, the module says ${fields} and ${stride}`); process.exit(1); }

const fail = msg => { console.error(msg); process.exit(1); };
const piece = Math.min(chunk, frames);
const pIn = M._malloc(piece * ch * 4);
const pL = M._malloc(fields * 8);
const RUNS = 5;                                                      // a short buffer on purpose — see the header
const pR = M._malloc(RUNS * stride * 8);
if (!pIn || !pL || !pR) fail('wasm OOM');

const out = [`# fcore stream v1 sr=${bitsOf(sr)} ch=${ch} chunk=${chunk}`];
let read = 0;
const drain = () => {
    for (;;) {
        const got = M._fc_stream_clips(h, read, pR, RUNS * stride);
        if (got === 0) break;
        const v = new Float64Array(M.HEAPF64.subarray(pR >>> 3, (pR >>> 3) + got * stride));
        for (let k = 0; k < got; ++k) {
            const b = k * stride;
            out.push(`run ${v[b]} ${v[b + 1]} ${bitsOf(v[b + 2])} ${v[b + 3]} ${v[b + 4]} ${v[b + 5]}`);
        }
        read += got;
    }
};
for (let at = 0; at < frames; at += chunk) {
    const m = Math.min(chunk, frames - at);
    // Planar with THIS piece's length as the stride: channel c at pIn + c*m.
    for (let c = 0; c < ch; ++c) M.HEAPF32.set(planes[c].subarray(at, at + m), (pIn >>> 2) + c * m);
    if (M._fc_stream_process(h, pIn, m) !== 1) fail(`fc_stream_process refused at frame ${at}`);
    if (M._fc_stream_loudness(h, pL) !== 1) fail('fc_stream_loudness refused');
    const L = new Float64Array(M.HEAPF64.subarray(pL >>> 3, (pL >>> 3) + fields));
    out.push(`at ${L[3]} m ${bitsOf(L[0])} s ${bitsOf(L[1])} i ${bitsOf(L[2])} dropped ${L[4]} runs ${M._fc_stream_clips_count(h)}`);
    drain();
}
if (M._fc_stream_finish(h) !== 1) fail('fc_stream_finish refused');
out.push(`finish runs ${M._fc_stream_clips_count(h)}`);
drain();
if (M._fc_stream_destroy(h) !== 1) fail('fc_stream_destroy did not know its own handle');
process.stdout.write(out.join('\n') + '\n');
