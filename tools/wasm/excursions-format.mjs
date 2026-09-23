// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure excursions` (K10). Byte-for-byte with the C++ half — the diff IS
// the test. Every double prints as its 16-hex-digit pattern for the reason the other formatters do: the
// two languages do not agree on decimal formatting, and a decimal header would break a whole-file diff
// while every measured bit matched.

import { bitsOf } from './blocks-format.mjs';

export function formatExcursions({ ok, s, runs, classes, crest, ceilMax, ceilDensity, ceilAbove12,
                                   sRun, nClasses, nCrest }) {
    if (!ok) return '';
    // THE WIDTHS ARE THE VERSION GATE. A row that grew on one side and not the other is how a parity
    // contract drifts in silence, so the mismatch is an error here and not a shorter line.
    if (sRun !== 6) throw new Error(`excursions v1 is 6 doubles per run, the module says ${sRun}`);
    if (nClasses !== 5) throw new Error(`excursions v1 has 5 duration classes, the module says ${nClasses}`);
    if (nCrest !== 30) throw new Error(`excursions v1 has 30 crest bins, the module says ${nCrest}`);
    const L = [];
    L.push(`# fcore excursions v1 sr=${bitsOf(s[0])} ch=${s[1]} ceiling=${bitsOf(s[2])} merge=${bitsOf(s[4])}`);
    L.push(`reason ${s[7]} valid ${s[8]} samples ${s[5]} measuredOs ${s[6]}`);
    L.push(`peaks recon ${bitsOf(s[9])} sample ${bitsOf(s[10])} truepeak ${bitsOf(s[11])}`);
    L.push(`runs ${s[12]} stored ${s[13]} complete ${s[14]} aboveOs ${s[15]}`);
    L.push(`occupancy ${bitsOf(s[16])} dose ${bitsOf(s[17])} maxexcess ${bitsOf(s[18])} p90 ${bitsOf(s[19])}`
           + ` sat ${s[20]} perminute ${bitsOf(s[21])}`);
    L.push('class count dose');
    for (let k = 0; k < nClasses; ++k) L.push(`c ${k} ${classes[2 * k]} ${bitsOf(classes[2 * k + 1])}`);
    L.push('crest lowHz count dose');
    for (let b = 0; b < nCrest; ++b)
        L.push(`k ${b} ${bitsOf(crest[3 * b])} ${crest[3 * b + 1]} ${bitsOf(crest[3 * b + 2])}`);
    L.push('run startOs lengthOs aboveOs peak dose crestHz');
    for (let i = 0; i < runs.length / sRun; ++i) {
        const r = runs.subarray(i * sRun, i * sRun + sRun);
        L.push(`r ${r[0]} ${r[1]} ${r[2]} ${bitsOf(r[3])} ${bitsOf(r[4])} ${bitsOf(r[5])}`);
    }
    L.push(`ceiling maxima ${ceilMax} density ${bitsOf(ceilDensity)} above12 ${bitsOf(ceilAbove12)}`);
    return L.join('\n') + '\n';
}
