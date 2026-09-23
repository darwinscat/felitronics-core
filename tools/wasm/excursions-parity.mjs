// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the K10 parity check:
//   fcore_measure excursions 48000 2 x.f32                            > native.txt
//   node excursions-parity.mjs build/fcprobe.node.js 48000 2 x.f32    > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: prepare() allocates, so memory.grow can fire inside
// _run and a view taken before it would then address freed memory.

import { formatExcursions } from './excursions-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node excursions-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
// THE SAME INPUT DOMAIN AS THE NATIVE TOOL, not a narrower or a wider one: `Number()` reads `0x2`, `+2`
// and `1e0` as counts the native tool refuses, and a road that accepts what the other rejects has already
// stopped being a parity check. The one residual is the one the older harnesses document: strtod also
// reads C99 hex floats, so `0x1p16` is a legal rate natively and is refused here — left refused, because
// a refusal is an exit 2 with no output, which is loud.
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

const raw = readFileSync(rawPath);
const inter = (raw.byteOffset % 4 === 0)
    ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
    : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
if (raw.byteLength % (4 * ch) !== 0) refuse(`the file is not a whole number of ${ch}-channel float32 frames`);
const frames = inter.length / ch;
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];

const ptr = M._malloc(planar.length * 4);
if (!ptr) refuse('wasm OOM on the input');
M.HEAPF32.set(planar, ptr >>> 2);
// The ceiling is the class's own default, which is what the native tool uses too — the two roads must
// take the same parameters or the diff compares two different measurements.
const ok = M._fc_probe_excursions_run(ptr, frames, ch, sr, -1.0) === 1;
M._free(ptr);
// A refused run exits 2, as fcore_measure does. Exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing — and the refusals are half of what parity means.
if (!ok) { process.exit(2); }

const sRun     = M._fc_probe_excursions_run_stride();
const nClasses = M._fc_probe_excursions_classes();
const nCrest   = M._fc_probe_excursions_crest_bins();

const pull = (fn, rows, stride) => {
    if (rows === 0) return new Float64Array(0);
    const p = M._malloc(rows * stride * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, rows * stride);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * stride).slice();
    M._free(p);
    return v;
};

const s       = pull((p, c) => M._fc_probe_excursions_scalars(p, c), M._fc_probe_excursions_scalars_len(), 1);
const classes = pull((p, c) => M._fc_probe_excursions_classes_out(p, c), nClasses, 2);
const crest   = pull((p, c) => M._fc_probe_excursions_crest(p, c), nCrest, 3);
const runs    = pull((p, c) => M._fc_probe_excursions_runs(p, c), s[13], sRun);

// 1e9 for "every local maximum", which is what the native road prints beside the 12 dB one.
const ceilMax      = M._fc_probe_excursions_ceiling_maxima();
const ceilDensity  = M._fc_probe_excursions_ceiling_density(1.0e9, 0.2);
const ceilAbove12  = M._fc_probe_excursions_ceiling_density(12.0, 0.2);

process.stdout.write(formatExcursions({
    ok: true, s, runs, classes, crest, ceilMax, ceilDensity, ceilAbove12, sRun, nClasses, nCrest,
}));
