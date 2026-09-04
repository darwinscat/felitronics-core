// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The one place the `blocks` stream is formatted on the JavaScript side, shared by the node parity harness
// and the page. Two copies of a format that a `diff` is the acceptance test for WILL drift — the page was
// already a version behind before this file existed.
//
// The stream must match `fcore_measure blocks` byte for byte. Everything numeric goes out as a raw IEEE-754
// bit pattern, including the sample rate: C's %g and JavaScript's Number-to-string disagree on a fractional
// rate (48000.123456 prints as 48000.12346 at ten significant digits), which would break the whole-file diff
// while every measured bit matched.

const dv = new DataView(new ArrayBuffer(8));

export function bitsOf(x) {
    dv.setFloat64(0, x, true);
    return dv.getBigUint64(0, true).toString(16).padStart(16, '0');
}

// `energies` may be any indexable of doubles; only the first `count` are emitted (pass a smaller count to
// show a preview, as the page does).
export function formatBlocks({ sr, ch, osFactor, osTaps, chunk, tp, sp, energies, count = energies.length }) {
    const lines = [
        `# fcore blocks v2 sr=${bitsOf(sr)} ch=${ch} os=${osFactor}x${osTaps} chunk=${chunk}`,
        `tp ${bitsOf(tp)}`,
        `sp ${bitsOf(sp)}`,
        `blocks ${energies.length}`,
    ];
    for (let i = 0; i < count; ++i) lines.push(bitsOf(energies[i]));
    return lines.join('\n') + '\n';
}
