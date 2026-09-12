// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the waveform-peaks / stereo-band parity check (P59a). Emits EXACTLY what
// `fcore_measure waveform|stereo|needle` emits, so the test is a diff:
//
//   fcore_measure waveform 44100 2 x.f32 --buckets 1100 --mix max   > native.txt
//   node shapes-parity.mjs build/fcprobe.node.js waveform 44100 2 x.f32 --buckets 1100 --mix max > wasm.txt
//   diff native.txt wasm.txt
//
// Usage: node shapes-parity.mjs <module.js> <waveform|stereo|needle> <sampleRate> <channels> <raw.f32le>
//        [--buckets N] [--mix avr|L|R|max] [--columns N] [--from A --to B]
//
// The detached-view rule of parity.mjs holds here too: no HEAP view is held across a call into the module, and
// every output is copied out of the heap before the next call.

import { formatNeedle, formatPeaks, formatStereo } from './shapes-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, mode, srArg, chArg, rawPath, ...rest] = process.argv;
if (!rawPath || !['waveform', 'stereo', 'needle'].includes(mode)) {
    console.error('usage: node shapes-parity.mjs <module.js> <waveform|stereo|needle> <sampleRate> <channels> <raw.f32le> [options]');
    process.exit(2);
}
// The native tool's option rules, so a command line cannot mean one thing there and another here: every argument is
// `--precise` or a known `--name value`, each at most once; counts are whole decimal numbers. Anything else exits 2,
// as `fcore_measure` does.
const refuse = msg => { console.error(msg); process.exit(2); };
const given = {};
for (let i = 0; i < rest.length; i++) {
    if (rest[i] === '--precise') continue;
    if (!['--buckets', '--columns', '--mix', '--from', '--to'].includes(rest[i]) || rest[i] in given || i + 1 >= rest.length)
        refuse(`bad, unknown or repeated option: ${rest[i]}`);
    given[rest[i]] = rest[++i];
}
const count = (s, name) => { if (!/^[0-9]{1,18}$/.test(s)) refuse(`bad ${name}: ${s}`); return Number(s); };
const opt = (name, fallback) => (name in given ? given[name] : fallback);
if (!/^[0-9]+(\.[0-9]+)?$/.test(srArg) || !(Number(srArg) > 0)) refuse('bad sampleRate');
if (!/^[0-9]{1,2}$/.test(chArg) || Number(chArg) < 1) refuse('bad channels');
const sr = Number(srArg), ch = Number(chArg);
const buckets = count(opt('--buckets', '1000'), '--buckets'), columns = count(opt('--columns', '1200'), '--columns');
const mixCode = ['avr', 'L', 'R', 'max'].indexOf(opt('--mix', 'avr'));

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
if (!ptr) { console.error('wasm OOM'); process.exit(1); }
M.HEAPF32.set(planar, ptr >>> 2);

if (mode === 'needle') {
    const from = count(opt('--from', '0'), '--from'), to = count(opt('--to', String(frames)), '--to');
    // Range-checked HERE, not left to the wasm call: a Number past 2^32 wraps to a uint32 on the way in, and the call
    // would then measure a different stretch than the one printed in the header.
    if (!(from <= to && to <= frames)) refuse('fc_probe_needle refused');
    const out = M._malloc(24);
    if (!M._fc_probe_needle(ptr, frames, ch, from, to, out)) { console.error('fc_probe_needle refused'); process.exit(2); }
    const v = new Float64Array(M.HEAPF64.subarray(out >>> 3, (out >>> 3) + 3));
    process.stdout.write(formatNeedle({ ch, frames, from, to, corr: v[0], width: v[1], rms: v[2] }));
    process.exit(0);
}

if (!M._fc_probe_shapes_run(ptr, frames, ch, sr, buckets, mixCode, columns)) { console.error('fc_probe_shapes_run refused'); process.exit(2); }

if (mode === 'waveform') {
    const n = M._fc_probe_waveform_count();
    const p64 = M._malloc(n * 8), p32 = M._malloc(n * 4);
    const w64 = M._fc_probe_waveform_peaks(p64, n), w32 = M._fc_probe_waveform_peaks_f32(p32, n);
    if (w64 !== n || w32 !== n) { console.error('short copy'); process.exit(1); }
    const peaks = new Float64Array(M.HEAPF64.subarray(p64 >>> 3, (p64 >>> 3) + n));
    const peaks32 = new Float32Array(M.HEAPF32.subarray(p32 >>> 2, (p32 >>> 2) + n));
    process.stdout.write(formatPeaks({ sr, ch, frames, mixCode, decim: M._fc_probe_waveform_decimation(),
                                       emitted: M._fc_probe_waveform_emitted(), peaks, peaks32 }));
} else {
    const cols = M._fc_probe_stereo_cols();
    const bufs = [M._malloc(cols * 4), M._malloc(cols * 4), M._malloc(cols * 4)];
    const ok = M._fc_probe_stereo_width(bufs[0], cols) === cols && M._fc_probe_stereo_corr(bufs[1], cols) === cols
            && M._fc_probe_stereo_rms(bufs[2], cols) === cols;
    if (!ok) { console.error('short copy'); process.exit(1); }
    const [width, corr, rms] = bufs.map(b => new Float32Array(M.HEAPF32.subarray(b >>> 2, (b >>> 2) + cols)));
    process.stdout.write(formatStereo({ ch, frames, cols, mono: M._fc_probe_stereo_is_mono(),
                                        maxRms: M._fc_probe_stereo_max_rms(), width, corr, rms }));
}
