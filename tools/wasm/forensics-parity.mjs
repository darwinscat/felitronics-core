// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the SourceForensics parity check (P77):
//   fcore_measure forensics 48000 2 x.f32                          > native.txt
//   node forensics-parity.mjs build/fcprobe.node.js 48000 2 x.f32   > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: every prepare() here allocates, so memory.grow can
// fire inside _run and a view taken before it would then address freed memory.

import { formatForensics } from './forensics-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node forensics-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
const sr = Number(srArg), ch = Number(chArg);
if (!Number.isFinite(sr) || sr <= 0) refuse(`bad sampleRate: ${srArg}`);
if (!Number.isInteger(ch) || ch < 1 || ch > 16) refuse(`bad channels: ${chArg}`);
// Integrality is not pedantry: `1.5` passes a Number() parse, survives the frame-size check, deinterleaves
// through fractional indices and reaches the module as 1 — a successful measurement of the wrong audio,
// while the native tool exits 2. The refusals have to match or the diff is comparing two different runs.

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
const ok = M._fc_probe_forensics_run(ptr, frames, ch, sr) === 1;
M._free(ptr);
// A refused run exits 2, as fcore_measure does. Exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing — and the refusals are half of what parity means: a byte
// diff of two SUCCESSFUL runs says nothing about the inputs both roads are supposed to reject.
if (!ok) { process.exit(2); }

const sWall    = M._fc_probe_forensics_wall_stride();
const sGrid    = M._fc_probe_forensics_grid_stride();
const kBuckets = M._fc_probe_forensics_khist_buckets();

const pull = (fn, rows, stride) => {
    if (rows === 0) return new Float64Array(0);
    const p = M._malloc(rows * stride * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, rows * stride);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * stride).slice();
    M._free(p);
    return v;
};

const s     = pull((p, c) => M._fc_probe_forensics_scalars(p, c), M._fc_probe_forensics_scalars_len(), 1);
const wall  = pull((p, c) => M._fc_probe_forensics_wall(p, c),  ch + 1, sWall);   // +1: the file's own wall
const grid  = pull((p, c) => M._fc_probe_forensics_grid(p, c),  ch, sGrid);
const khist = pull((p, c) => M._fc_probe_forensics_khist(p, c), ch, kBuckets);

process.stdout.write(formatForensics({ ok, s, wall, grid, khist, sWall, sGrid, kBuckets }));
