// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure crest` (K1). It must match that output byte for byte — the diff IS
// the parity test — so every floating-point number crosses as a raw IEEE-754 bit pattern and every count as a
// decimal integer, exactly as the C++ half prints them.
//
// THE SCALAR AND ROW ORDERS ARE THE CONTRACT. `fc_probe_crest_scalars` and `fc_probe_crest_blocks` write them
// in the order the CLI prints them, and the indices below name that order once. The two halves are a pair;
// neither may be edited alone, and a mismatch shows up as a diff rather than as a plausible wrong number.

import { bitsOf } from './blocks-format.mjs';

const S = {
    sr: 0, ch: 1, hop: 2, blockHops: 3, e0: 4, e1: 5, e2: 6, floor: 7, share: 8,
    samples: 9, hopCount: 10, baseHops: 11, blockCount: 12, nonFinite: 13, reason: 14, peak: 15,
    active0: 16,   // one accepted-block count per band, in band order
    programmeMs: 21,   // the programme's level in the gate's units, dBFS
};

// One loss row, in the order `fc_probe_crest_loss` writes it.
const L = {
    blocks: 0, inActive: 1, usable: 2, p50: 3, p95: 4, cvar95: 5, mean: 6, max: 7, p5: 8,
    over1: 9, over3: 10, over6: 11, peakShift: 12, levelShift: 13,
    outSilent: 14, lagBlocks: 15, valid: 16,
};

export const CREST_BANDS = 5;

export function formatCrest({ ok, s, blocks, blockStride, loss }) {
    if (!ok) return '';
    if (blockStride !== 15)
        throw new Error(`crest v1 is 15 doubles per block row, the module says ${blockStride}`);
    const out = [];
    out.push(`# fcore crest v1 sr=${bitsOf(s[S.sr])} ch=${s[S.ch]} hop=${s[S.hop]}`
           + ` blockHops=${s[S.blockHops]} e0=${bitsOf(s[S.e0])} e1=${bitsOf(s[S.e1])} e2=${bitsOf(s[S.e2])}`
           + ` floor=${bitsOf(s[S.floor])} share=${bitsOf(s[S.share])} chunk=8192`);
    out.push(`samples ${s[S.samples]}`);
    out.push(`hops ${s[S.hopCount]} ${s[S.baseHops]} ${s[S.blockCount]} ${s[S.nonFinite]} ${s[S.reason]}`);
    out.push(`peak ${bitsOf(s[S.peak])}`);
    {
        let line = 'active';
        for (let k = 0; k < CREST_BANDS; ++k) line += ` ${s[S.active0 + k]}`;
        line += ` ${bitsOf(s[S.programmeMs])}`;
        out.push(line);
    }
    for (let j = 0; j < s[S.blockCount]; ++j) {
        const b = j * blockStride;
        let line = `b ${j}`;
        for (let k = 0; k < CREST_BANDS; ++k)
            line += ` ${bitsOf(blocks[b + k * 3])} ${bitsOf(blocks[b + k * 3 + 1])} ${blocks[b + k * 3 + 2]}`;
        out.push(line);
    }
    if (loss) {
        for (let band = 0; band < CREST_BANDS; ++band) {
            const r = loss[band];
            out.push(`loss ${band} ${r[L.blocks]} ${r[L.inActive]} ${r[L.usable]}`
                   + ` ${bitsOf(r[L.p50])} ${bitsOf(r[L.p95])} ${bitsOf(r[L.cvar95])} ${bitsOf(r[L.mean])}`
                   + ` ${bitsOf(r[L.max])} ${bitsOf(r[L.p5])}`
                   + ` ${r[L.over1]} ${r[L.over3]} ${r[L.over6]}`
                   + ` ${bitsOf(r[L.peakShift])} ${bitsOf(r[L.levelShift])}`
                   + ` ${r[L.outSilent]} ${r[L.lagBlocks]} ${r[L.valid]}`);
        }
    }
    return out.join('\n') + '\n';
}
