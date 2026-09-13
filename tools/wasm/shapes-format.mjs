// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The one place the `waveform` / `stereo` / `needle` streams are formatted on the JavaScript side. They must match
// `fcore_measure waveform|stereo|needle` byte for byte — a `diff` of the two IS the native-vs-wasm parity test — so
// every number goes out as a raw IEEE-754 bit pattern (16 hex digits for a double, 8 for a float32), the same rule
// and the same reason as blocks-format.mjs.

import { bitsOf } from './blocks-format.mjs';

const dv = new DataView(new ArrayBuffer(4));
export function bits32Of(x) {
    dv.setFloat32(0, x, true);
    return dv.getUint32(0, true).toString(16).padStart(8, '0');
}

const MIX_NAMES = ['avr', 'L', 'R', 'max'];

export function formatPeaks({ sr, ch, frames, mixCode, decim, emitted, peaks, peaks32 }) {
    const lines = [`# fcore waveform v1 sr=${bitsOf(sr)} ch=${ch} frames=${frames} buckets=${peaks.length} mix=${MIX_NAMES[mixCode]} decim=${decim} emitted=${emitted}`];
    for (let i = 0; i < peaks.length; ++i) lines.push(`${bitsOf(peaks[i])} ${bits32Of(peaks32[i])}`);
    return lines.join('\n') + '\n';
}

export function formatStereo({ ch, frames, cols, mono, maxRms, width, corr, rms }) {
    const lines = [`# fcore stereo v1 ch=${ch} frames=${frames} cols=${cols} mono=${mono ? 1 : 0}`, `maxrms ${bitsOf(maxRms)}`];
    for (let i = 0; i < cols; ++i) lines.push(`${bits32Of(width[i])} ${bits32Of(corr[i])} ${bits32Of(rms[i])}`);
    return lines.join('\n') + '\n';
}

export function formatNeedle({ ch, frames, from, to, corr, width, rms }) {
    return [`# fcore needle v1 ch=${ch} frames=${frames} from=${from} to=${to}`,
            `corr ${bitsOf(corr)}`, `width ${bitsOf(width)}`, `rms ${bitsOf(rms)}`].join('\n') + '\n';
}
