// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the HumDetector parity check (P77):
//   fcore_measure hum 48000 2 x.f32                          > native.txt
//   node hum-parity.mjs build/fcprobe.node.js 48000 2 x.f32   > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: every prepare() here allocates, so memory.grow can
// fire inside _run and a view taken before it would then address freed memory.

import { formatHum } from './hum-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node hum-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
const sr = Number(srArg), ch = Number(chArg);

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
const ok = M._fc_probe_hum_run(ptr, frames, ch, sr) === 1;
M._free(ptr);
if (!ok) { process.stdout.write(''); process.exit(0); }

const sChan    = M._fc_probe_hum_chan_stride();
const sCand    = M._fc_probe_hum_cand_stride();
const sHarm    = M._fc_probe_hum_harm_stride();
const sStretch = M._fc_probe_hum_stretch_stride();

// Each read: allocate, call, COPY OUT, free — the view never outlives the next call into the module.
const pull = (fn, rows, stride) => {
    if (rows === 0) return new Float64Array(0);
    const p = M._malloc(rows * stride * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, rows * stride);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * stride).slice();
    M._free(p);
    return v;
};

const s = pull((p, c) => M._fc_probe_hum_scalars(p, c), M._fc_probe_hum_scalars_len(), 1);
const nCand = s[6], maxHarm = s[7];
const chan    = pull((p, c) => M._fc_probe_hum_chan(p, c),    ch, sChan);
const cand    = pull((p, c) => M._fc_probe_hum_cand(p, c),    ch * nCand, sCand);
const harm    = pull((p, c) => M._fc_probe_hum_harm(p, c),    ch * nCand * maxHarm, sHarm);
// Stretches are capped by the detector's own storage; ask generously and keep what came back.
const stretch = pull((p, c) => M._fc_probe_hum_stretch(p, c), ch * 4096, sStretch);

process.stdout.write(formatHum({ ok, s, chan, cand, harm, stretch, sChan, sCand, sHarm, sStretch }));
