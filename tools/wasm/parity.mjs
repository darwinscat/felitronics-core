// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the P0 parity check. Emits EXACTLY the format `fcore_measure blocks` emits, so the whole
// test is:
//
//   fcore_measure blocks 48000 2 x.f32  > native.txt
//   node parity.mjs build/fcprobe.node.js 48000 2 x.f32 > wasm.txt
//   diff native.txt wasm.txt
//
// Usage: node parity.mjs <module.js> <sampleRate> <channels> <raw.f32le> [--debug]
//
// THE DETACHED-VIEW TRAP, which is the thing this file exists to get right. With ALLOW_MEMORY_GROWTH, a
// `memory.grow` replaces the underlying ArrayBuffer and every existing typed-array view onto it DETACHES:
// writes are silently dropped, reads come back undefined. Memory grows in two places here — the big _malloc,
// and again INSIDE fc_probe_run (the meter's 1.27 MB prepare can cross a growth boundary). So no HEAP view is
// ever held across a call: `Module.HEAPF32` is re-read immediately before every use, results come back as
// return values rather than through a view, and the output buffer is read only after the run has finished.

import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
if (!rawPath) {
    console.error('usage: node parity.mjs <module.js> <sampleRate> <channels> <raw.f32le> [--debug]');
    process.exit(2);
}
const sr = Number(srArg);
const ch = Number(chArg);

const bitsOf = (() => {
    const dv = new DataView(new ArrayBuffer(8));
    return (x) => { dv.setFloat64(0, x, true); return dv.getBigUint64(0, true).toString(16).padStart(16, '0'); };
})();

const require = createRequire(import.meta.url);
const createFcProbe = require(resolve(modPath));
const M = await createFcProbe();

// The file is interleaved f32le, exactly the bytes the native tool reads. De-interleave into planar — a pure
// float permutation, so it cannot move a bit; the point is that both sides measure the same samples.
const raw = readFileSync(rawPath);
// A Float32Array view needs a 4-byte-aligned byteOffset, which readFileSync's pooled Buffer does not
// promise; copy when it is not aligned rather than throwing on a file that would otherwise be fine.
const inter = (raw.byteOffset % 4 === 0)
    ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
    : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
const frames = Math.floor(inter.length / ch);
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c)
    for (let i = 0; i < frames; ++i)
        planar[c * frames + i] = inter[i * ch + c];

const bytes = planar.length * 4;
const ptr = M._malloc(bytes);
if (!ptr) { console.error('wasm OOM allocating ' + bytes + ' bytes'); process.exit(1); }

// Re-read HEAPF32 AFTER the malloc: that malloc may itself have grown memory and detached any earlier view.
M.HEAPF32.set(planar, ptr >>> 2);          // >>> not >>: a pointer above 2 GB is negative under a signed shift

const ok = M._fc_probe_run(ptr, frames, ch, sr);
if (!ok) { console.error('fc_probe_run rejected its arguments'); process.exit(1); }

const nBlocks = M._fc_probe_block_count();
const outPtr = M._malloc(nBlocks * 8);
if (!outPtr) { console.error('wasm OOM allocating the block-energy buffer'); process.exit(1); }
const written = M._fc_probe_block_energies(outPtr, nBlocks);

// HEAPF64 re-read here too, for the same reason — the second _malloc could have grown memory again. And
// COPIED out immediately rather than held as a view: every line below this is another call into the module,
// and a rule that depends on nobody ever reordering the code is not a rule.
const energies = new Float64Array(M.HEAPF64.subarray(outPtr >>> 3, (outPtr >>> 3) + written));

const tp = M._fc_probe_tp_linear();
const dropped = M._fc_probe_dropped();

let out = `# fcore blocks v1 sr=${bitsOf(sr)} ch=${ch} os=${M._fc_probe_os_factor()}x${M._fc_probe_os_taps()} chunk=${M._fc_probe_chunk()}\n`;
out += `tp ${bitsOf(tp)}\n`;
out += `sp ${bitsOf(M._fc_probe_sample_peak())}\n`;
out += `blocks ${written}\n`;
for (let i = 0; i < written; ++i) out += bitsOf(energies[i]) + '\n';
process.stdout.write(out);

if (dropped !== 0) console.error(`warning: ${dropped} gating blocks dropped`);
if (process.argv.includes('--debug'))
    console.error(`[debug] sizeof(long double) in this wasm = ${M._fc_probe_sizeof_longdouble()} bytes, ` +
                  `lufs = ${M._fc_probe_lufs(ptr, frames, ch, sr)}, dbtp = ${M._fc_probe_dbtp(ptr, frames, ch, sr)}`);

M._free(outPtr);
M._free(ptr);
