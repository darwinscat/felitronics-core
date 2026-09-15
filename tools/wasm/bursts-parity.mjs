// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the BandBursts parity check (P77):
//   fcore_measure bursts 48000 2 x.f32                          > native.txt
//   node bursts-parity.mjs build/fcprobe.node.js 48000 2 x.f32   > wasm.txt
//   diff native.txt wasm.txt
// No HEAP view is held across a call into the module: every prepare() here allocates, so memory.grow can
// fire inside _run and a view taken before it would then address freed memory.

import { formatBursts } from './bursts-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath) refuse('usage: node bursts-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le>');
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
const ok = M._fc_probe_bursts_run(ptr, frames, ch, sr) === 1;
M._free(ptr);
// A refused run exits 2, as fcore_measure does. Exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing — and the refusals are half of what parity means: a byte
// diff of two SUCCESSFUL runs says nothing about the inputs both roads are supposed to reject.
if (!ok) { process.exit(2); }

const strideChan = M._fc_probe_bursts_chan_stride();
const strideEvt  = M._fc_probe_bursts_evt_stride();
const strideBin  = M._fc_probe_bursts_bin_stride();

// Each read: allocate, call, COPY OUT, free — the view never outlives the next call.
const pull = (fn, rows, stride) => {
    if (rows === 0) return new Float64Array(0);
    const p = M._malloc(rows * stride * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, rows * stride);
    const v = new Float64Array(M.HEAPF64.buffer, p, got * stride).slice();
    M._free(p);
    return v;
};

const s = pull((p, c) => M._fc_probe_bursts_scalars(p, c) / 1, M._fc_probe_bursts_scalars_len(), 1);
const chan   = pull((p, c) => M._fc_probe_bursts_chan(p, c),   ch, strideChan);
const events = pull((p, c) => M._fc_probe_bursts_events(p, c), s[24], strideEvt);
// The histograms are sparse; ask for the worst case the detector can hold and keep what came back.
const ioi = pull((p, c) => M._fc_probe_bursts_ioi(p, c), 4096, strideBin);
const lag = pull((p, c) => M._fc_probe_bursts_lag(p, c), 4096, strideBin);

process.stdout.write(formatBursts({ ok, s, chan, events, ioi, lag, strideChan, strideEvt, strideBin }));
