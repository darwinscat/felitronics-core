// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// Writes the interleaved f32le fixture the `clips` parity check runs on: a programme that IS CLAMPED, so the
// comparison has runs to compare rather than two agreeing zeroes. make-fixture.mjs is not that programme —
// its flat stretches are square waves and a full-scale hit that simply stops, which analysis::ClipDetector
// correctly reports as clean (a step from rest is no ramp), so a parity run on it only proves the two sides
// agree that nothing is there.
//
// DETERMINISM BY CONSTRUCTION, the same rule as make-fixture.mjs and for the same reason: not one
// transcendental is called. Math.sin is not required by ECMAScript to be correctly rounded, so a fixture
// built from it would be a different file on a different engine and the byte-identity claim would rest on
// the engine. Everything here is integer arithmetic plus float32 storage — exactly specified by IEEE-754.
//
// WHY IT IS BUILT IN 16-BIT CODE UNITS. The detector's tau is 2q, and q is the finer of the smallest non-zero
// step the stream has shown and the coarsest 2^-k grid every sample lies on. Generating integer codes and
// dividing by 32768 puts every sample exactly on the 2^-15 lattice (k/32768 is exact in float32 for |k| <=
// 32768), so q is pinned at 2^-15 and the bounds are evaluated against the quantum a 16-bit master really
// has. The `0, 1, 0, -1` prelude is there to SHOW the stream a one-code step in the first four samples, so q
// is that fine from the start rather than after the first quiet passage.
//
// WHAT MAKES A RUN FINDABLE HERE, bound by bound (ClipDetector.h:36-83):
//   · the crest is a PARABOLA — v = 1 - u^2 over an integer period — so the samples outside the clamp are far
//     enough below it to pass the smooth-crest bound V;
//   · the noise is added BEFORE the clamp, never after. Before, it raises the flank activity sigma and so
//     lowers the chance bound's ratio, while the clamped top stays EXACTLY flat (every sample is the ceiling
//     code) and is therefore exempt from the rough-turn test. After the clamp it would roughen the top, and
//     one noisy sample above the ceiling kills the candidate outright;
//   · the approach is moving, not a step from rest, so no flat pair sits before the entry;
//   · the ceiling is an ordinary code (29196, about -1 dBFS), not full scale: clipped-then-turned-down is
//     still clipped, and a fixture that clamps at 1.0 would let a reader think the rule is about full scale.
//
// Usage: node make-clip-fixture.mjs <out.f32> [sampleRate=48000] [channels=2] [seconds=4]

import { writeFileSync } from 'node:fs';

const [, , outPath, srArg = '48000', chArg = '2', secArg = '4'] = process.argv;
if (!outPath) { console.error('usage: node make-clip-fixture.mjs <out.f32> [sampleRate] [channels] [seconds]'); process.exit(2); }

const sr = Number(srArg), ch = Number(chArg), seconds = Number(secArg);
const frames = Math.round(sr * seconds);
const CEIL = 29196;                                   // the clamp, in 16-bit codes (about -1 dBFS)

// xorshift32 — integer only, so every sample below is a pure function of the seed.
function rng(seed) {
    let x = seed >>> 0;
    return () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x; };
}

// A parabolic crest train of integer period `p`: alternating positive and negative bumps, apex `amp` codes,
// zero at every half-period junction. `1 - u*u` with u = 2*(i mod p)/p - 1, in integer arithmetic until the
// last division — so the shape is DETERMINISTIC on any engine (an exact integer numerator below 2^53, then
// one IEEE division and one Math.round, both exactly specified). It is not always EXACT: p*p is a power of
// two only for channel 0's period, and the quotient is rounded for the others. Determinism is what the
// fixture's pinned bytes need; exactness was never the claim.
function crest(i, p, amp) {
    const k = ((i % (2 * p)) + 2 * p) % (2 * p);
    const sign = k < p ? 1 : -1;
    const t = k % p;                                   // 0 .. p-1
    const u = 2 * t - p;                               // -p .. p-2, i.e. u/p in [-1, 1)
    return sign * amp * (p * p - u * u) / (p * p);
}

// THE PASSAGES ARE FRACTIONS OF THE FILE, NOT SECONDS INTO IT. They were absolute at first, and the
// consequence was that every fixture CI actually generates (1 and 2 seconds) stopped before the holes, the
// second silence and the clipped tail — three of the seven passages, including the two that carry the cases
// hardest to get right, were unreachable in the whole parity corpus while the generator looked as though it
// produced them. As fractions, every duration exercises every passage.
const SECTIONS = [0.0625, 0.125, 0.40, 0.475, 0.725, 0.75, 0.80];

const out = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) {
    const rand = rng(0x1f123bb5 + c * 0x27d4eb2f);     // a different, fixed seed per channel
    const b = k => Math.round(frames * SECTIONS[k]);
    // Per channel: a different crest period, so the two channels' runs neither align nor have equal lengths
    // — a report that mixed the channels up would still look plausible with one period.
    const P = c === 0 ? 128 : 96;
    for (let i = 0; i < frames; ++i) {
        let code;
        if (i < 4) {
            code = [0, 1, 0, -1][i];                   // the prelude that shows q = 1 code
        } else if (i < b(0)) {
            code = (rand() % 4001) - 2000;             // quiet noise: clean, and it is where q stays pinned
        } else if (i < b(1)) {
            code = 0;                                  // digital silence — no bands, no runs, no q
        } else if (i < b(2)) {
            // HARD CLAMP. Driven 25 % over the ceiling, so the flat top is about 0.447*P samples long.
            code = Math.round(crest(i, P, Math.round(CEIL * 1.25))) + ((rand() % 61) - 30);
        } else if (i < b(3)) {
            code = (rand() % 8001) - 4000;             // a clean passage between the two clamped sections
        } else if (i < b(4)) {
            // A SHALLOW CLAMP, 0.4 % over the ceiling: flat tops of only a handful of samples, which is the
            // length the chance bound actually argues about.
            code = Math.round(crest(i, P * 2, Math.round(CEIL * 1.004))) + ((rand() % 21) - 10);
        } else if (i < b(5)) {
            // THE HOLES. Non-finite samples break a run and are never part of a peak, a DC or a decision;
            // they sit in their own passage, well away from every witness, because an unknown sample on a
            // flank makes the run beside it undecidable rather than wrong.
            const r = rand() % 97;
            code = r === 0 ? 'nan' : r === 1 ? 'inf' : r === 2 ? '-inf' : (rand() % 2001) - 1000;
        } else if (i < b(6)) {
            code = 0;
        } else {
            // THE TAIL. Clamped right up to the last sample, so the final runs are decided only by finish():
            // the decision lags the stream by 20 ms, and a report read without finish() is missing them.
            code = Math.round(crest(i, P, Math.round(CEIL * 1.25))) + ((rand() % 61) - 30);
        }
        let v;
        if (code === 'nan') v = NaN;
        else if (code === 'inf') v = Infinity;
        else if (code === '-inf') v = -Infinity;
        else v = Math.max(-CEIL, Math.min(CEIL, code)) / 32768;    // the clamp, then onto the 2^-15 lattice
        out[i * ch + c] = v;
    }
}

// Written little-endian EXPLICITLY, as make-fixture.mjs is: a TypedArray's backing bytes are in the host's
// order, and the native tool reads f32**le**.
//
// AND THE NaN IS WRITTEN AS A BIT PATTERN, not as a Number. `writeFloatLE(NaN, …)` stores an
// implementation-chosen encoding — ECMAScript's NumericToRawBytes lets an engine pick any NaN — so a fixture
// that carried NaN through a float store would not be the same BYTES on another engine, and the determinism
// claim above would be false exactly where it is hardest to notice. 0x7fc00000 is the canonical quiet NaN.
// (Nothing downstream can tell the difference — the detector only asks isfinite() — but the fixture's bytes
// are a published number, and a published number has to be reproducible.)
const bytes = Buffer.allocUnsafe(out.length * 4);
for (let i = 0; i < out.length; ++i) {
    if (Number.isNaN(out[i])) bytes.writeUInt32LE(0x7fc00000, i * 4);
    else bytes.writeFloatLE(out[i], i * 4);
}
writeFileSync(outPath, bytes);
console.log(`${outPath}: ${frames} frames x ${ch} ch @ ${sr} Hz (${out.byteLength} bytes), clamped at ${CEIL}/32768`);
