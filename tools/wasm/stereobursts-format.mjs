// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure stereobursts` (K3c). Byte-for-byte with the C++ half — the diff IS
// the test. Every double prints as its 16-hex-digit pattern, for the reason the other formatters give: the
// two languages do not agree on decimal formatting, and a decimal header would break a whole-file diff
// while every measured bit matched.

import { bitsOf } from './blocks-format.mjs';

export function formatStereoBursts({ ok, s, events, sEvt }) {
    if (!ok) return '';
    // THE WIDTHS ARE THE VERSION GATE, as everywhere here: a row that grew on one side and not the other
    // is how a parity contract drifts in silence, so a mismatch is an error and not a shorter line.
    if (sEvt !== 13) throw new Error(`stereobursts v1 is 13 doubles per event, the module says ${sEvt}`);
    if (s.length < 36) throw new Error(`stereobursts v1 has 36 scalars, the module gave ${s.length}`);
    const L = [];
    L.push(`# fcore stereobursts v1 sr=${bitsOf(s[0])} ch=${s[1]} hop=${s[9]} base=${s[10]}`
           + ` lo=${bitsOf(s[3])} hi=${bitsOf(s[4])} enter=${bitsOf(s[7])} exit=${bitsOf(s[8])} chunk=${s[34]}`);
    L.push(`samples ${s[2]} ran ${s[35]}`);
    L.push(`dome ${bitsOf(s[11])} ${bitsOf(s[12])}`);
    L.push(`side ${s[13]} ${s[14]} ${s[15]}`);
    for (let a = 0; a < 2; ++a) {
        // Nine per axis, Mid first, at a fixed stride — so one axis is indexed as a block and Mid's count
        // can never be read against Side's baseline.
        const o = 16 + 9 * a;
        L.push(`axis ${a} ${s[o]} ${s[o + 1]} ${s[o + 2]} ${s[o + 3]}`
               + ` ${s[o + 4]} ${s[o + 5]} ${s[o + 6]} ${s[o + 7]} ${s[o + 8]}`);
        const rows = events[a];
        for (let i = 0; i < rows.length / sEvt; ++i) {
            const e = rows.subarray(i * sEvt, i * sEvt + sEvt);
            L.push(`e ${a} ${e[0]} ${e[1]} ${e[2]} ${e[8]} ${bitsOf(e[3])} ${bitsOf(e[4])} ${bitsOf(e[5])}`
                   + ` ${bitsOf(e[6])} ${bitsOf(e[7])} ${bitsOf(e[9])} ${bitsOf(e[10])} ${e[11]} ${e[12]}`);
        }
    }
    return L.join('\n') + '\n';
}
