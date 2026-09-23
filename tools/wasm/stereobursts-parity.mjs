// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the K3c parity check:
//   fcore_measure stereobursts 48000 2 x.f32                            > native.txt
//   node stereobursts-parity.mjs build/fcprobe.node.js 48000 2 x.f32    > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: prepare() allocates, so memory.grow can fire inside
// _run and a view taken before it would then address freed memory.

import { formatStereoBursts } from './stereobursts-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node stereobursts-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
// THE SAME INPUT DOMAIN AS THE NATIVE TOOL, not a narrower or a wider one — a road that accepts what the
// other rejects has already stopped being a parity check. The residual the older harnesses document
// stands: strtod also reads C99 hex floats, so `0x1p16` is a legal rate natively and is refused here.
const count = (s, name) => { if (!/^[0-9]{1,18}$/.test(s)) refuse(`bad ${name}: ${s}`); return Number(s); };
const rateOf = s => {
    if (!/^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/.test(s)) refuse(`bad sampleRate: ${s}`);
    const v = Number(s);
    if (!Number.isFinite(v) || v <= 0) refuse(`bad sampleRate: ${s}`);
    return v;
};
const sr = rateOf(srArg), ch = count(chArg, 'channels');
if (ch < 1 || ch > 16) refuse(`bad channels: ${chArg}`);

const require = createRequire(import.meta.url);
const M = await require(resolve(modPath))();

// A FILE THAT IS NOT THERE IS A REFUSAL, NOT A CRASH. An uncaught readFileSync leaves node exiting 1,
// while the native tool refuses with 2 — two roads disagreeing on a refusal, which is the half of parity a
// diff of two successful runs can never show.
let raw;
try { raw = readFileSync(rawPath); } catch { refuse(`cannot open ${rawPath}`); }
const inter = (raw.byteOffset % 4 === 0)
    ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
    : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
if (raw.byteLength % (4 * ch) !== 0) refuse(`the file is not a whole number of ${ch}-channel float32 frames`);
const frames = inter.length / ch;
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];

// AN EMPTY PROGRAMME IS A LEGAL MEASUREMENT HERE, as it is for the native road, so it must not travel
// through malloc(0) — which returns 0 and reads exactly like an allocation failure. The core takes a null
// pointer with a zero frame count; what it must never take is a null with a non-zero one.
let ptr = 0;
if (planar.length !== 0) {
    ptr = M._malloc(planar.length * 4);
    if (!ptr) refuse('wasm OOM on the input');
    M.HEAPF32.set(planar, ptr >>> 2);
}
// The class's own defaults, which is what the native tool uses with no flags — the two roads must take
// the same parameters or the diff compares two different measurements.
const ok = M._fc_probe_stereobursts_run(ptr, frames, ch, sr) === 1;
if (ptr) M._free(ptr);
// A refused run exits 2, as fcore_measure does. Exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing — and the refusals are half of what parity means.
if (!ok) { process.exit(2); }

const sEvt = M._fc_probe_stereobursts_evt_stride();
const nSc  = M._fc_probe_stereobursts_scalars_len();

const pull = (fn, elems) => {
    if (elems === 0) return new Float64Array(0);
    const p = M._malloc(elems * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, elems);
    const v = new Float64Array(M.HEAPF64.buffer, p, got).slice();
    M._free(p);
    return v;
};

const s = pull((p, c) => M._fc_probe_stereobursts_scalars(p, c), nSc);
if (s.length !== nSc) refuse(`the module published ${s.length} of ${nSc} scalars`);

// `stored` per axis lives at 16 + 9*a + 7; the copier takes ELEMENTS, as every copier in this ABI does,
// so the request is rows * stride and never rows.
const events = [0, 1].map(a => {
    const stored = s[16 + 9 * a + 7];
    if (stored === 0) return new Float64Array(0);
    const p = M._malloc(stored * sEvt * 8);
    if (!p) refuse('wasm OOM');
    const got = M._fc_probe_stereobursts_events(a, p, stored * sEvt);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * sEvt).slice();
    M._free(p);
    return v;
});

process.stdout.write(formatStereoBursts({ ok: true, s, events, sEvt }));
