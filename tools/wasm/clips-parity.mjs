// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the clipped-runs parity check (P71). Emits EXACTLY what `fcore_measure clips` emits, so
// the test is a diff:
//
//   fcore_measure clips 48000 2 x.f32 --max-runs 500                              > native.txt
//   node clips-parity.mjs build/fcprobe.node.js clips 48000 2 x.f32 --max-runs 500 > wasm.txt
//   diff native.txt wasm.txt
//
// Usage: node clips-parity.mjs <module.js> clips <sampleRate> <channels> <raw.f32le> [--max-runs N] [--chunk N]
//
// `--chunk` IS ACCEPTED AND DELIBERATELY NOT HONOURED HERE, which makes the diff stronger rather than
// sloppier. The module measures the whole buffer in one fc_probe_clips_run call; the native tool cuts the
// stream into process() calls of the requested size. So a diff taken with `--chunk 7` compares a
// seven-sample-at-a-time native run against a one-call wasm run, and law 8a — the report cannot depend on
// where the stream was cut — is then part of what the parity test proves, across two toolchains at once.
// Without this, the two sides would both end up cutting at multiples of fcore::ClipProbe::kChunk and the diff
// would say nothing at all about re-slicing.
//
// The detached-view rule of parity.mjs holds here too: no HEAP view is held across a call into the module,
// and every output is copied out of the heap before the next call.

import { formatClips } from './clips-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, mode, srArg, chArg, rawPath, ...rest] = process.argv;
if (!rawPath || mode !== 'clips') {
    console.error('usage: node clips-parity.mjs <module.js> clips <sampleRate> <channels> <raw.f32le> [--max-runs N] [--chunk N]');
    process.exit(2);
}
// The native tool's option rules, so a command line cannot mean one thing there and another here: every
// argument is `--precise` or a known `--name value`, each at most once; counts are whole decimal numbers.
const refuse = msg => { console.error(msg); process.exit(2); };
const given = {};
for (let i = 0; i < rest.length; i++) {
    if (rest[i] === '--precise') continue;
    if (!['--max-runs', '--chunk'].includes(rest[i]) || rest[i] in given || i + 1 >= rest.length)
        refuse(`bad, unknown or repeated option: ${rest[i]}`);
    given[rest[i]] = rest[++i];
}
// THE SAME DOMAIN AS THE NATIVE TOOL, NOT A NARROWER ONE. Two roads that disagree about which command lines
// are legal is a defect of its own: the diff then fails, or one side quietly measures something else, for a
// reason that has nothing to do with the measurement. So each rule here mirrors the C it stands opposite.
//
//   · the rate: `parseRate` is strtod plus "the whole string was consumed" (fcore_measure.cpp), which accepts
//     `1e3`, `1000.` and `+1000`. `Number()` accepts exactly those and the same rejections, EXCEPT for
//     JavaScript's own literal forms (`0b1`, `0o7`, `1_0`, `Infinity`) and the empty string, which strtod does
//     not read — so those are refused explicitly rather than left to agree by luck. ONE RESIDUAL, measured:
//     strtod also reads C99 hex floats, so `0x1p10` is a legal 1024 Hz natively and is refused here. It is
//     left refused rather than reimplemented: a spelling nobody types, and a refusal on this side is an
//     exit 2 with no output, which is loud — not the silent divergence the bounds below exist to stop. The
//     same goes for strtod's leading whitespace: `" 48000"` measures natively and is refused here.
//   · the counts: `parseCount` is up to 18 decimal digits and nothing else.
//   · the bounds: the width against core::kMaxChannels, the capacity against fcore::ClipProbe::kMaxRuns, the
//     chunk against the int it is narrowed to. WITHOUT THE CAPACITY BOUND HERE, `--max-runs 4294967296` wraps
//     to 0 on the way into a uint32_t parameter and this side reports a successful run of `maxruns=0` while
//     the native tool exits 2 — a silent disagreement, found by the crew's review round.
const count = (s, name) => { if (!/^[0-9]{1,18}$/.test(s)) refuse(`bad ${name}: ${s}`); return Number(s); };
const rate = s => {
    if (!/^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/.test(s)) refuse('bad sampleRate');
    const v = Number(s);
    if (!Number.isFinite(v) || !(v > 0)) refuse('bad sampleRate');
    return v;
};
const sr = rate(srArg), ch = count(chArg, 'channels');
if (ch < 1 || ch > 16) refuse('bad sampleRate/channels');            // core::kMaxChannels
// 65536 is ClipDetectorParams::maxRuns' own default (ClipDetector.h). It is spelled again here rather than
// exported, because the header line carries `maxruns=` on both sides: a drift between the two defaults shows
// up as a failing diff on the very first comparison instead of hiding until a file overflows.
const maxRuns = '--max-runs' in given ? count(given['--max-runs'], '--max-runs') : 65536;
if (maxRuns > 1048576) refuse('clips.prepare refused');              // fcore::ClipProbe::kMaxRuns
if ('--chunk' in given && count(given['--chunk'], '--chunk') > 0x7FFFFFFF) refuse('--chunk out of range');

const require = createRequire(import.meta.url);
const M = await require(resolve(modPath))();

const raw = readFileSync(rawPath);
const inter = (raw.byteOffset % 4 === 0)
    ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
    : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
if (raw.byteLength % (4 * ch) !== 0) refuse(`the file is not a whole number of ${ch}-channel float32 frames`);
const frames = inter.length / ch;
if (frames === 0) refuse('the file is empty');                       // as `fcore_measure clips` refuses it
const planar = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];

const ptr = M._malloc(planar.length * 4);
if (!ptr) { console.error('wasm OOM'); process.exit(1); }
M.HEAPF32.set(planar, ptr >>> 2);

const ok = M._fc_probe_clips_run(ptr, frames, ch, sr, maxRuns) === 1;
if (!ok) { console.error('fc_probe_clips_run refused'); process.exit(2); }

const channels = M._fc_probe_clips_channels();
const stored = M._fc_probe_clips_stored();
const stride = M._fc_probe_clips_stride();

const pPeaks = M._malloc(channels * 8);
if (M._fc_probe_clips_peaks(pPeaks, channels) !== channels) { console.error('short peak copy'); process.exit(1); }
const peaks = new Float64Array(M.HEAPF64.subarray(pPeaks >>> 3, (pPeaks >>> 3) + channels));

let runs = new Float64Array(0);
if (stored > 0) {
    const pRuns = M._malloc(stored * stride * 8);
    if (!pRuns) { console.error('wasm OOM'); process.exit(1); }
    // `cap` is in DOUBLES here, as it is in every other copier of this ABI; the return is in RUNS.
    if (M._fc_probe_clips_runs(pRuns, stored * stride) !== stored) { console.error('short run copy'); process.exit(1); }
    runs = new Float64Array(M.HEAPF64.subarray(pRuns >>> 3, (pRuns >>> 3) + stored * stride));
}

process.stdout.write(formatClips({
    ok,
    sr: M._fc_probe_clips_rate(),
    ch: channels,
    frames: M._fc_probe_clips_samples(),
    maxRuns: M._fc_probe_clips_max_runs(),
    delay: M._fc_probe_clips_delay(),
    count: M._fc_probe_clips_count(),
    stored, complete: M._fc_probe_clips_complete(), peaks, runs, stride,
}));
