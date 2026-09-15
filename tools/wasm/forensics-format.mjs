// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure forensics`. Byte-for-byte with the C++ half — the diff IS the test.
// Note the wall table has one more row than there are channels: the last is the FILE's aggregate wall.

import { bitsOf } from './blocks-format.mjs';

const S = { sr:0, ch:1, order:2, hop:3, bins:4, percell:5, exempt:6, distinctlimit:7, plateaucells:8,
            floorcells:9, firstParam:10, samples:26, tail:27, frames:28 };

export function formatForensics({ ok, s, wall, grid, khist, sWall, sGrid, kBuckets }) {
    if (!ok) return '';
    if (sWall !== 39 || sGrid !== 22)
        throw new Error(`forensics v1 is 39/22 doubles per row, the module says ${sWall}/${sGrid}`);
    const L = [];
    L.push(`# fcore forensics v1 sr=${bitsOf(s[S.sr])} ch=${s[S.ch]} order=${s[S.order]} hop=${s[S.hop]} bins=${s[S.bins]} percell=${s[S.percell]}`);
    let p = `params exempt=${s[S.exempt]} distinctlimit=${s[S.distinctlimit]} plateaucells=${s[S.plateaucells]} floorcells=${s[S.floorcells]}`;
    for (let i = 0; i < 16; ++i) p += ` ${bitsOf(s[S.firstParam + i])}`;
    L.push(p);
    L.push(`samples ${s[S.samples]} tail ${s[S.tail]} frames ${s[S.frames]}`);
    for (let i = 0; i < wall.length / sWall; ++i) {
        const w = wall.subarray(i * sWall, i * sWall + sWall);
        L.push(`wall ${w[0]} valid=${w[1]} reason=${w[2]} sharp=${w[3]} nearnyq=${w[4]} clipped=${w[5]} trunc=${w[6]} exempted=${w[7]}`
             + ` second=${w[8]}/${w[9]}/${w[10]}/${w[11]} secondreason=${w[12]} empty=${w[13]} emptyreason=${w[14]}`
             + ` frames=${w[15]}/${w[16]}`);
        let v = `wall ${w[0]} values`;
        for (let k = 17; k < 39; ++k) v += ` ${bitsOf(w[k])}`;
        L.push(v);
    }
    for (let i = 0; i < grid.length / sGrid; ++i) {
        const g = grid.subarray(i * sGrid, i * sGrid + sGrid);
        L.push(`grid ${g[0]} valid=${g[1]} reason=${g[2]} k=${g[3]} pcm=${g[4]} outofrange=${g[5]} bits=${g[6]} robustk=${g[7]}`
             + ` robustbits=${g[8]} zerolow24=${g[9]} peak=${bitsOf(g[10])} min=${bitsOf(g[11])} max=${bitsOf(g[12])}`);
        L.push(`grid ${g[0]} nonzero=${g[13]} zero=${g[14]} nonfinite=${g[15]} absent=${g[16]} offgrid=${g[17]}`
             + ` firstoffgrid=${g[18]} firstmaxk=${g[19]} distinct=${g[20]} complete=${g[21]}`);
        let h = `grid ${g[0]} khist`;
        for (let k = 0; k < kBuckets; ++k) h += ` ${khist[i * kBuckets + k]}`;
        L.push(h);
    }
    return L.join('\n') + '\n';
}
