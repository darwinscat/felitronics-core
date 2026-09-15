// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The one place the `report` stream is formatted on the JavaScript side. It must match
// `fcore_measure report` byte for byte — a `diff` of the two IS the parity test — so every floating-point
// number goes out as a raw IEEE-754 bit pattern, the same rule and the same reason as blocks-format.mjs.
//
// THE NAMES ARE NOT WRITTEN HERE. They come out of the module, produced by the same visitor walk that
// produced the rows, because ProgrammeReport's visitor is the single enumeration of its fields and a
// hand-written copy of that list on this side would be a second one — which is exactly the drift the
// visitor exists to prevent. A field added to the report therefore appears in this output with nothing
// edited here.
//
// EVERY NUMBER IS READ BACK OUT OF THE MEASUREMENT, never taken from the command line: a header assembled
// from the caller's own arithmetic would agree with the native tool's header even if the module had been
// fed half the buffer.

import { bitsOf } from './blocks-format.mjs';

// `names` is the NUL-separated blob, counts first then values. `counts` and `values` are flat Float64Array
// slices of stride 2 and 4. `ok` is the run's own verdict, mirrored from the C++ half: a report that is
// not ok prints NOTHING, because a caller can forget a bool but cannot forget an empty stdout.
export function formatReport({ ok, sr, ch, samples, names, counts, values, strideC, strideV }) {
    if (!ok) return '';
    if (strideC !== 2 || strideV !== 4) throw new Error(`report v1 is 2/4 doubles per row, the module says ${strideC}/${strideV}`);
    const nC = counts.length / strideC, nV = values.length / strideV;
    if (names.length !== nC + nV) throw new Error(`the module gave ${names.length} names for ${nC + nV} rows`);
    const lines = [`# fcore report v1 sr=${bitsOf(sr)} ch=${ch} samples=${samples}`];
    for (let i = 0; i < nC; ++i) {
        const b = i * strideC;
        lines.push(`C ${names[i]} ${counts[b]} ${counts[b + 1]}`);
    }
    for (let i = 0; i < nV; ++i) {
        const b = i * strideV;
        lines.push(`V ${names[nC + i]} ${values[b]} ${values[b + 1]} ${values[b + 2]} ${bitsOf(values[b + 3])}`);
    }
    return lines.join('\n') + '\n';
}
