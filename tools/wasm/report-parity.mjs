// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the ProgrammeReport parity check (P77). Emits EXACTLY what `fcore_measure report`
// emits, so the test is a diff:
//
//   fcore_measure report 48000 2 x.f32                          > native.txt
//   node report-parity.mjs build/fcprobe.node.js 48000 2 x.f32   > wasm.txt
//   diff native.txt wasm.txt
//
// The detached-view rule of parity.mjs holds here too: no HEAP view is held across a call into the
// module. Every prepare() in these analyzers allocates, so memory.grow can fire inside _run and a view
// taken before it would then address freed memory.

import { formatReport } from './report-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = msg => { console.error(msg); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node report-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
const sr = Number(srArg), ch = Number(chArg);
if (!Number.isFinite(sr) || sr <= 0) refuse(`bad sampleRate: ${srArg}`);
if (!Number.isInteger(ch) || ch < 1) refuse(`bad channels: ${chArg}`);

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

const ok = M._fc_probe_report_run(ptr, frames, ch, sr) === 1;
M._free(ptr);

if (!ok) { process.stdout.write(''); process.exit(0); }

const strideC = M._fc_probe_report_stride_counts();
const strideV = M._fc_probe_report_stride_values();
const nC = M._fc_probe_report_count_rows();
const nV = M._fc_probe_report_value_rows();

// Sized from the MODULE, never from a number this file made up.
const nameBytes = M._fc_probe_report_names(0, 0);
const nPtr = M._malloc(nameBytes);
if (!nPtr) refuse('wasm OOM on the name table');
const wrote = M._fc_probe_report_names(nPtr, nameBytes);
const blob = Buffer.from(M.HEAPU8 ? M.HEAPU8.slice(nPtr, nPtr + wrote) : new Uint8Array(M.HEAPF64.buffer, nPtr, wrote));
M._free(nPtr);
const names = blob.toString('latin1').split('\0').slice(0, -1);

const cPtr = M._malloc(nC * strideC * 8);
const rowsC = M._fc_probe_report_counts(cPtr, nC * strideC);
const counts = new Float64Array(M.HEAPF64.buffer, cPtr, rowsC * strideC).slice();
M._free(cPtr);

const vPtr = M._malloc(nV * strideV * 8);
const rowsV = M._fc_probe_report_values(vPtr, nV * strideV);
const values = new Float64Array(M.HEAPF64.buffer, vPtr, rowsV * strideV).slice();
M._free(vPtr);

process.stdout.write(formatReport({
    ok, sr, ch, samples: M._fc_probe_report_samples(), names, counts, values, strideC, strideV,
}));
