// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the LowEnd parity check (P77):
//   fcore_measure lowend 48000 2 x.f32                          > native.txt
//   node lowend-parity.mjs build/fcprobe.node.js 48000 2 x.f32   > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: every prepare() here allocates, so memory.grow can
// fire inside _run and a view taken before it would then address freed memory.

import { formatLowEnd } from './lowend-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node lowend-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
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
const ok = M._fc_probe_lowend_run(ptr, frames, ch, sr) === 1;
M._free(ptr);
if (!ok) { process.stdout.write(''); process.exit(0); }

const sSeries = M._fc_probe_lowend_series_stride();
const sBand   = M._fc_probe_lowend_band_stride();
const hBins   = M._fc_probe_lowend_hist_bins();

const pull = (fn, rows, stride) => {
    if (rows === 0) return new Float64Array(0);
    const p = M._malloc(rows * stride * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, rows * stride);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * stride).slice();
    M._free(p);
    return v;
};

const s      = pull((p, c) => M._fc_probe_lowend_scalars(p, c), M._fc_probe_lowend_scalars_len(), 1);
const hist   = pull((p, c) => M._fc_probe_lowend_hist(p, c), hBins, 1);
const series = pull((p, c) => M._fc_probe_lowend_series(p, c), s[26], sSeries);
const bands  = pull((p, c) => M._fc_probe_lowend_bands(p, c), s[6], sBand);

// The note name comes from the module; its length does too.
let noteName = '';
const nBytes = M._fc_probe_lowend_note_name(0, 0);
if (nBytes > 0) {
    const np = M._malloc(nBytes);
    if (!np) refuse('wasm OOM on the note name');
    const wrote = M._fc_probe_lowend_note_name(np, nBytes);
    // HEAPU8 is not among the exported views (build.sh declares HEAPF32/HEAPF64), so the bytes are
    // read through the one buffer that IS exported. Copied out before the next call, as everywhere here.
    noteName = Buffer.from(new Uint8Array(M.HEAPF64.buffer, np, wrote)).toString('latin1');
    M._free(np);
}

process.stdout.write(formatLowEnd({ ok, s, hist, series, bands, noteName, sSeries, sBand, hBins }));
