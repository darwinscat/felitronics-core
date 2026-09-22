// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of the BandCrest parity check (K1):
//   fcore_measure crest 48000 2 a.f32 --against b.f32                        > native.txt
//   node crest-parity.mjs build/fcprobe.node.js 48000 2 a.f32 --against b.f32 > wasm.txt
//   diff native.txt wasm.txt
// Flags, with the same names, defaults and refusals on both roads: --edge0 --edge1 --edge2 --hop-ms
// --block-hops --floor-db --share-db --against.
//
// No HEAP view is held across a call into the module: every prepare() here allocates, so memory.grow can fire
// inside a _run and a view taken before it would then address freed memory.

import { formatCrest, CREST_BANDS } from './crest-format.mjs';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, srArg, chArg, rawPath] = process.argv;
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || !srArg || !chArg || !rawPath)
    refuse('usage: node crest-parity.mjs <module.js> <sampleRate> <channels> <raw.f32le> [flags]');

// THE SAME INPUT DOMAIN AS THE NATIVE TOOL, not a narrower or a wider one — the policy P71 established.
const count = (s, name) => { if (!/^[0-9]{1,18}$/.test(s)) refuse(`bad ${name}: ${s}`); return Number(s); };
const rateOf = s => {
    if (!/^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/.test(s)) refuse(`bad sampleRate: ${s}`);
    const v = Number(s);
    if (!Number.isFinite(v) || v <= 0) refuse(`bad sampleRate: ${s}`);
    return v;
};
const finiteOf = (s, name) => {
    if (!/^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/.test(s)) refuse(`bad ${name}: ${s}`);
    const v = Number(s);
    if (!Number.isFinite(v)) refuse(`bad ${name}: ${s}`);
    return v;
};
const sr = rateOf(srArg), ch = count(chArg, 'channels');
if (ch < 1 || ch > 16) refuse(`bad channels: ${chArg}`);

// BandCrestParams' documented defaults, repeated because the module publishes no reader for them. The
// repetition is guarded: the scalar block prints all of them, so a drift between these and the C++ ones shows
// up as a diff on the first line of the output, which is what this harness is for.
const P = { '--edge0': 120, '--edge1': 2000, '--edge2': 6000, '--hop-ms': 100,
            '--floor-db': -70, '--share-db': -40 };
let blockHops = 4, against = null;
{
    const rest = process.argv.slice(6);
    for (let i = 0; i < rest.length; ++i) {
        if (Object.prototype.hasOwnProperty.call(P, rest[i])) {
            if (i + 1 >= rest.length) refuse(`crest: ${rest[i]} needs a finite number`);
            P[rest[i]] = finiteOf(rest[i + 1], rest[i]); ++i; continue;
        }
        if (rest[i] === '--block-hops') {
            if (i + 1 >= rest.length || !/^[0-9]{1,18}$/.test(rest[i + 1])) refuse('crest: --block-hops needs 1..64');
            blockHops = Number(rest[i + 1]);
            if (blockHops < 1 || blockHops > 64) refuse('crest: --block-hops needs 1..64');
            ++i; continue;
        }
        if (rest[i] === '--against') {
            if (i + 1 >= rest.length) refuse('crest: --against needs a file');
            against = rest[i + 1]; ++i; continue;
        }
        // AN OPTION NOBODY KNOWS IS A REFUSAL, on this road too and with the same exit.
        if (rest[i].startsWith('--') && rest[i] !== '--precise')
            refuse(`crest: unknown option ${rest[i]}`);
    }
}

const require = createRequire(import.meta.url);
const M = await require(resolve(modPath))();

const planarise = path => {
    const raw = readFileSync(path);
    const inter = (raw.byteOffset % 4 === 0)
        ? new Float32Array(raw.buffer, raw.byteOffset, Math.floor(raw.byteLength / 4))
        : new Float32Array(raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength));
    if (raw.byteLength % (4 * ch) !== 0) refuse(`${path} is not a whole number of ${ch}-channel float32 frames`);
    const frames = inter.length / ch;
    const planar = new Float32Array(frames * ch);
    for (let c = 0; c < ch; ++c) for (let i = 0; i < frames; ++i) planar[c * frames + i] = inter[i * ch + c];
    return { planar, frames };
};

const runSlot = (slot, path) => {
    const { planar, frames } = planarise(path);
    const ptr = M._malloc(Math.max(4, planar.length * 4));
    if (!ptr) refuse('wasm OOM on the input');
    M.HEAPF32.set(planar, ptr >>> 2);
    const ok = M._fc_probe_crest_run_with(slot, ptr, frames, ch, sr,
                                          P['--edge0'], P['--edge1'], P['--edge2'],
                                          P['--hop-ms'], blockHops,
                                          P['--floor-db'], P['--share-db']) === 1;
    M._free(ptr);
    return ok;
};

// A refused run exits 2, as fcore_measure does: exiting 0 with empty output would tell a caller the
// measurement succeeded and produced nothing, and the refusals are half of what parity means.
if (!runSlot(0, rawPath)) process.exit(2);
if (against !== null && !runSlot(1, against)) process.exit(2);

const pull = (fn, doubles) => {
    if (doubles === 0) return new Float64Array(0);
    const p = M._malloc(doubles * 8);
    if (!p) refuse('wasm OOM');
    const got = fn(p, doubles);
    const v = new Float64Array(M.HEAPF64.buffer, p, doubles).slice();
    M._free(p);
    return { v, got };
};

const scal = pull((p, c) => M._fc_probe_crest_scalars(0, p, c), M._fc_probe_crest_scalars_len());
const s = scal.v;
const blockStride = M._fc_probe_crest_block_stride();
const rows = s[12];                                   // blockCount — see crest-format.mjs's index table
const blk = pull((p, c) => M._fc_probe_crest_blocks(0, p, c), rows * blockStride);

let loss = null;
if (against !== null) {
    loss = [];
    const len = M._fc_probe_crest_loss_len();
    for (let band = 0; band < CREST_BANDS; ++band) loss.push(pull((p, c) => M._fc_probe_crest_loss(band, p, c), len).v);
}

process.stdout.write(formatCrest({ ok: true, s, blocks: blk.v ?? blk, blockStride, loss }));
