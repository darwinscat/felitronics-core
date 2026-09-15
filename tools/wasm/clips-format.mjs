// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The one place the `clips` stream is formatted on the JavaScript side. It must match `fcore_measure clips`
// byte for byte — a `diff` of the two IS the native-vs-wasm parity test — so every floating-point number
// goes out as a raw IEEE-754 bit pattern, the same rule and the same reason as blocks-format.mjs. The C++
// half of this format is tools/fcore_clips_format.h; the two are a pair and neither may be edited alone.
//
// WHY BIT PATTERNS AND NOT DECIMALS, for the level in particular. A clipped run's level is the MEAN of its
// samples, so it is not a grid value: a 12-sample run holding eleven codes of 40 and one of 39 on a 1/64
// lattice has level 479/768, whose double is 3fe3f55555555555 and whose float32 is 3fe3f55560000000. That
// pair is the whole argument. A `%.6f` column prints 0.623698 for both, and a comparison that cannot tell a
// double from its float32 rounding cannot tell anything else about the last bits either.
//
// EVERY NUMBER HERE IS READ BACK OUT OF THE MEASUREMENT, never taken from the command line — see the ABI
// section in fc_probe.cpp. A header assembled from the caller's own arithmetic would agree with the native
// tool's header even if the module had been fed half the buffer.

import { bitsOf } from './blocks-format.mjs';

// `runs` is a flat Float64Array of stride*count doubles: start, length, level, channel, sign, evidence.
//
// `ok` IS THE C++ HALF'S RULE, MIRRORED. formatClips() in fcore_clips_format.h returns nothing at all for a
// report that is not ok, because a caller can forget a bool but cannot forget an empty stdout — and this half
// needs the same rule for the same reason, with one extra teeth: fc_probe_clips_stride() answers 6 even when
// the module holds no result, so a consumer that ignored fc_probe_clips_run()'s return value and read the
// getters anyway would assemble `sr=0000000000000000 ch=0 frames=0 … runs 0 stored 0 complete 0` — exactly the
// plausible empty report both halves exist to forbid. Pass the run's own verdict as `ok`.
export function formatClips({ ok, sr, ch, frames, maxRuns, delay, count, stored, complete, peaks, runs, stride }) {
    if (!ok) return '';
    if (stride !== 6) throw new Error(`clips v1 is 6 doubles per run, the module says ${stride}`);
    const lines = [
        `# fcore clips v1 sr=${bitsOf(sr)} ch=${ch} frames=${frames} maxruns=${maxRuns} delay=${delay}`,
        `runs ${count} stored ${stored} complete ${complete ? 1 : 0}`,
    ];
    for (let c = 0; c < peaks.length; ++c) lines.push(`peak ${c} ${bitsOf(peaks[c])}`);
    for (let i = 0; i < stored; ++i) {
        const b = i * stride;
        // start, length, channel, sign and evidence are integers held exactly in a double (every position in
        // this ABI is below 2^30); `level` is the only one whose bits are the point.
        lines.push(`run ${runs[b]} ${runs[b + 1]} ${bitsOf(runs[b + 2])} ${runs[b + 3]} ${runs[b + 4]} ${runs[b + 5]}`);
    }
    return lines.join('\n') + '\n';
}
