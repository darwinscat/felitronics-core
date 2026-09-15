// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure lowend`. Byte-for-byte with the C++ half — the diff IS the test.
// The note's NAME comes from the module, not from a pitch-class table rebuilt here: that table already
// exists in C++ and a second copy is the drift this ABI is shaped to avoid.

import { bitsOf } from './blocks-format.mjs';

export function formatLowEnd({ ok, s, hist, series, bands, noteName, sSeries, sBand, hBins }) {
    if (!ok) return '';
    if (sSeries !== 6 || sBand !== 11)
        throw new Error(`lowend v1 is 6/11 doubles per row, the module says ${sSeries}/${sBand}`);
    const L = [];
    L.push(`# fcore lowend v1 sr=${bitsOf(s[0])} ch=${s[1]} xover=${bitsOf(s[2])} order=${s[3]} hop=${s[4]} block=${s[5]} bands=${s[6]} chunk=${s[7]}`);
    L.push(`reason ${s[8]} ${s[9]}`);
    L.push(`samples ${s[10]} finite ${s[11]} holes ${s[12]} nonfinite ${s[13]} absent ${s[14]} overflow ${s[15]}`);
    const en = ['lowmid', 'lowside', 'highmid', 'highside', 'rawmid', 'rawside', 'lowfrac'];
    for (let i = 0; i < 7; ++i) L.push(`${en[i]} ${bitsOf(s[16 + i])}`);
    L.push(`highfrac ${bitsOf(s[23])} rawfrac ${bitsOf(s[24])}`);
    L.push(`blocks ${s[25]} stored ${s[26]} complete ${s[27]} histsamples ${s[28]}`);
    for (let i = 0; i < hBins; ++i) L.push(`h${String(i).padStart(2, '0')} ${hist[i]}`);
    L.push(`worst ${s[29]} ${bitsOf(s[30])} ${bitsOf(s[31])}`);
    L.push(`peakenergy ${s[32]} ${bitsOf(s[33])} ${bitsOf(s[34])}`);
    L.push(`peakside ${s[35]} ${bitsOf(s[36])} amp ${bitsOf(s[37])} at ${s[38]}`);
    L.push(`frames used ${s[39]} holed ${s[40]} tail ${s[41]} window ${s[42]} underresolved ${s[43]}`);
    L.push('series index samples finite holes midEnergy sideEnergy');
    for (let i = 0; i < series.length / sSeries; ++i) {
        const r = series.subarray(i * sSeries, i * sSeries + sSeries);
        L.push(`s ${r[0]} ${r[1]} ${r[2]} ${r[3]} ${bitsOf(r[4])} ${bitsOf(r[5])}`);
    }
    L.push('band midi centreHz widthHz binsPerBand midEnergy sideEnergy energy density centroidHz centsOffset');
    for (let i = 0; i < bands.length / sBand; ++i) {
        const r = bands.subarray(i * sBand, i * sBand + sBand);
        let line = `b ${r[0]} ${r[1]}`;
        for (let k = 2; k < 11; ++k) line += ` ${bitsOf(r[k])}`;
        L.push(line);
    }
    L.push(`peak ${s[44]} ${s[45]} density ${s[46]} second ${s[47]}`);
    if (s[48] === 1) L.push(`note ${noteName} nominal ${bitsOf(s[49])} centroid ${bitsOf(s[50])} cents ${bitsOf(s[51])} sidefrac ${bitsOf(s[52])}`);
    else             L.push(`note INVALID reason ${s[9]}`);
    L.push(`frame ${bitsOf(s[53])} rangeshare ${bitsOf(s[54])}`);
    L.push(`background ${bitsOf(s[55])} peakenergy ${bitsOf(s[56])} peakwidth ${bitsOf(s[57])} share ${bitsOf(s[58])} total ${bitsOf(s[59])}`);
    return L.join('\n') + '\n';
}
