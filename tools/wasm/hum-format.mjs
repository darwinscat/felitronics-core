// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure hum`. Byte-for-byte with the C++ half — the diff IS the test.
// The rows carry their own (channel, candidate, harmonic) coordinates rather than implying them from
// position, so this file GROUPS by those coordinates instead of re-deriving the CLI's loop nesting.

import { bitsOf } from './blocks-format.mjs';

const H = { sr: 0, ch: 1, order: 2, n: 3, hop: 4, bin: 5, cands: 6, maxHarm: 7 };

export function formatHum({ ok, s, chan, cand, harm, stretch, sChan, sCand, sHarm, sStretch }) {
    if (!ok) return '';
    if (sChan !== 21 || sCand !== 19 || sHarm !== 8 || sStretch !== 5)
        throw new Error(`hum v1 is 21/19/8/5 doubles per row, the module says ${sChan}/${sCand}/${sHarm}/${sStretch}`);
    const L = [];
    L.push(`# fcore hum v1 sr=${bitsOf(s[H.sr])} ch=${s[H.ch]} order=${s[H.order]} n=${s[H.n]} hop=${s[H.hop]} bin=${bitsOf(s[H.bin])}`);
    const rows = (a, st) => { const o = []; for (let i = 0; i < a.length / st; ++i) o.push(a.subarray(i * st, i * st + st)); return o; };
    const R = rows(chan, sChan), K = rows(cand, sCand), Hh = rows(harm, sHarm), St = rows(stretch, sStretch);
    for (const r of R) {
        const c = r[0];
        L.push(`ch ${c} valid ${r[1]} reason ${r[2]} mains ${r[3]} base ${r[4]} fobs ${r[5]} fderived ${r[6]}`);
        L.push(`ch ${c} f0 ${bitsOf(r[7])} hz ${bitsOf(r[8])} tone ${bitsOf(r[9])} peakbin ${bitsOf(r[10])} floor ${bitsOf(r[11])} prom ${bitsOf(r[12])}`);
        L.push(`ch ${c} frames ${r[13]} finite ${r[14]} holed ${r[15]} quiet ${r[16]} stretches ${r[17]} stored ${r[18]} complete ${r[19]} tail ${r[20]}`);
        for (const q of K.filter(x => x[0] === c)) {
            const k = q[1];
            L.push(`ch ${c} cand ${bitsOf(q[2])} found ${q[3]} base ${q[4]} f0 ${bitsOf(q[5])}`
                 + ` sobs ${q[6]} soff ${q[7]} fobs ${q[8]} sspread ${bitsOf(q[9])} fspread ${bitsOf(q[10])}`
                 + ` intra ${bitsOf(q[11])} stat ${q[12]} pass ${q[13]} harm ${q[14]} low ${q[15]}`);
            L.push(`ch ${c} cand ${k} window ${q[16]} hz ${bitsOf(q[17])} prom ${bitsOf(q[18])}`);
            for (const h of Hh.filter(x => x[0] === c && x[1] === k))
                L.push(`ch ${c} cand ${k} h ${h[2]} inband ${h[3]} acc ${h[4]} hz ${bitsOf(h[5])} tone ${bitsOf(h[6])} prom ${bitsOf(h[7])}`);
        }
        for (const t of St.filter(x => x[0] === c))
            L.push(`ch ${c} stretch ${t[1]} ${t[2]} ${t[3]} frames ${t[4]}`);
    }
    return L.join('\n') + '\n';
}
