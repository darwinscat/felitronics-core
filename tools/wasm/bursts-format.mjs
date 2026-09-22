// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript half of `fcore_measure bursts`. It must match that output byte for byte — the diff IS
// the parity test — so every floating-point number crosses as a raw IEEE-754 bit pattern and every count
// as a decimal integer, exactly as the C++ half prints them.
//
// THE SCALAR ORDER IS THE CONTRACT. fc_probe_bursts_scalars() writes them in the order the CLI prints
// them, and the indices below name that order once. The two halves are a pair; neither may be edited
// alone, and a mismatch shows up as a diff rather than as a plausible wrong number.

import { bitsOf } from './blocks-format.mjs';

const S = {
    sr: 0, ch: 1, hop: 2, base: 3, lo: 4, hi: 5, enter: 6, exit: 7, chunk: 8,
    samples: 9, hopCount: 10, eligible: 11, zeroBaseline: 12, burstHops: 13,
    damaged: 14, overflow: 15, firstNonFinite: 16,
    tailSamples: 17, tailEnergy: 18,
    eventsValid: 19, eventsReason: 20, progValid: 21, progReason: 22,
    eventCount: 23, storedEvents: 24, eventsComplete: 25,
    onsets: 26, intervals: 27, intervalOverflow: 28, modalHops: 29, modalMass: 30,
    onsetsPerSec: 31,
};

export function formatBursts({ ok, s, chan, events, ioi, lag, strideChan, strideEvt, strideBin }) {
    if (!ok) return '';
    if (strideChan !== 4 || strideEvt !== 12 || strideBin !== 2)
        throw new Error(`bursts v1 is 4/12/2 doubles per row, the module says ${strideChan}/${strideEvt}/${strideBin}`);
    const L = [];
    L.push(`# fcore bursts v1 sr=${bitsOf(s[S.sr])} ch=${s[S.ch]} hop=${s[S.hop]} base=${s[S.base]}`
         + ` lo=${bitsOf(s[S.lo])} hi=${bitsOf(s[S.hi])} enter=${bitsOf(s[S.enter])} exit=${bitsOf(s[S.exit])}`
         + ` chunk=${s[S.chunk]}`);
    L.push(`samples ${s[S.samples]}`);
    L.push(`hops ${s[S.hopCount]} ${s[S.eligible]} ${s[S.zeroBaseline]} ${s[S.burstHops]}`);
    L.push(`damage ${s[S.damaged]} ${s[S.overflow]} ${s[S.firstNonFinite]}`);
    L.push(`tail ${s[S.tailSamples]} ${bitsOf(s[S.tailEnergy])}`);
    L.push(`valid ${s[S.eventsValid]} ${s[S.eventsReason]} ${s[S.progValid]} ${s[S.progReason]}`);
    for (let i = 0; i < chan.length / strideChan; ++i) {
        const b = i * strideChan;
        L.push(`chan ${chan[b]} ${bitsOf(chan[b + 1])} ${chan[b + 2]} ${chan[b + 3]}`);
    }
    L.push(`events ${s[S.eventCount]} ${s[S.storedEvents]} ${s[S.eventsComplete]}`);
    for (let i = 0; i < events.length / strideEvt; ++i) {
        const b = i * strideEvt;
        L.push(`e ${events[b]} ${events[b + 1]} ${events[b + 2]} ${events[b + 3]}`
             + ` ${bitsOf(events[b + 4])} ${bitsOf(events[b + 5])} ${bitsOf(events[b + 6])}`
             + ` ${bitsOf(events[b + 7])} ${bitsOf(events[b + 8])}`
             + ` ${events[b + 9]}${events[b + 10]}${events[b + 11]}`);
    }
    L.push(`onsets ${s[S.onsets]} ${s[S.intervals]} ${s[S.intervalOverflow]} ${s[S.modalHops]} ${s[S.modalMass]}`
         + ` ${bitsOf(s[S.onsetsPerSec])}`);
    for (let i = 0; i < ioi.length / strideBin; ++i) L.push(`ioi ${ioi[i * strideBin]} ${ioi[i * strideBin + 1]}`);
    for (let i = 0; i < lag.length / strideBin; ++i) L.push(`lag ${lag[i * strideBin]} ${lag[i * strideBin + 1]}`);
    return L.join('\n') + '\n';
}
