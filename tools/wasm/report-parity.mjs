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
// THE SAME INPUT DOMAIN AS THE NATIVE TOOL, NOT A NARROWER OR WIDER ONE — the policy P71 established and
// P77 did not carry forward until the release round found it. `Number()` is not `parseCount`/`parseRate`:
// it reads `0x2`, `+2` and `1e0` as channel counts the native tool refuses, so those succeeded here and
// exited 2 there. The counts are decimal digits and nothing else, as fcore_measure's parseCount is; the
// rate mirrors strtod's grammar minus JavaScript's own literal forms. ONE RESIDUAL, the same one P71
// documents: strtod also reads C99 hex floats, so `0x1p16` is a legal rate natively and is refused here —
// left refused rather than reimplemented, because a refusal is an exit 2 with no output, which is loud,
// and not the silent divergence this grammar exists to stop.
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

const ok = M._fc_probe_report_run(ptr, frames, ch, sr) === 1;
M._free(ptr);

// A refused run exits 2, as fcore_measure does. Exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing — and the refusals are half of what parity means: a byte
// diff of two SUCCESSFUL runs says nothing about the inputs both roads are supposed to reject.
if (!ok) { process.exit(2); }

const strideC = M._fc_probe_report_stride_counts();
const strideV = M._fc_probe_report_stride_values();
const nC = M._fc_probe_report_count_rows();
const nV = M._fc_probe_report_value_rows();

// Sized from the MODULE, never from a number this file made up.
const nameBytes = M._fc_probe_report_names(0, 0);
const nPtr = M._malloc(nameBytes);
if (!nPtr) refuse('wasm OOM on the name table');
const wrote = M._fc_probe_report_names(nPtr, nameBytes);
// HEAPU8 is NOT among the exported views (build.sh declares HEAPF32/HEAPF64 only), and in the CHECKED
// module merely touching an unexported runtime property aborts — so a `M.HEAPU8 ? … : …` fallback never
// reaches its fallback. The bytes are read through the one buffer that is exported.
const blob = Buffer.from(new Uint8Array(M.HEAPF64.buffer, nPtr, wrote));
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
