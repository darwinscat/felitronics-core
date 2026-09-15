// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// A programme that HAS mains hum in QUIET STRETCHES, because the loudness fixture does not and the parity
// gate was silently proving nothing with it: on make-fixture.mjs output HumDetector finds zero quiet
// frames and returns NoQuietStretch before it reaches its peak search at all, so 47 lines of "nothing
// found" satisfied `test -s` and compared clean without the measurement path ever running.
//
// TRANSCENDENTAL-FREE, like make-fixture.mjs and for the same reason: a fixture built with Math.sin would
// not be the same bytes on every machine that regenerates it, and a fixture that is not reproducible
// cannot anchor a byte-exact comparison. The tones come from the Chebyshev recurrence
// s[n] = c*s[n-1] - s[n-2] with c = 2*cos(2*pi*f/fs) written as a literal, computed once at 40 digits.
//
// Usage: node make-hum-fixture.mjs <out.f32> [sampleRate=48000] [channels=2] [seconds=10]
// The rate is fixed at 48000 for the literals below; another rate is refused rather than silently wrong.

import { writeFileSync } from 'node:fs';

const [, , outPath, srArg = '48000', chArg = '2', secArg = '30'] = process.argv;
if (!outPath) { console.error('usage: node make-hum-fixture.mjs <out.f32> [sampleRate] [channels] [seconds]'); process.exit(2); }
const sr = Number(srArg), ch = Number(chArg), secs = Number(secArg);
if (sr !== 48000) { console.error('this generator carries 48 kHz literals only; regenerate them for another rate'); process.exit(2); }
if (!Number.isInteger(ch) || ch < 1 || ch > 16) { console.error(`bad channels: ${chArg}`); process.exit(2); }

// 2*cos(2*pi*f/48000) and sin(2*pi*f/48000), at 40 digits, for the mains line and three harmonics.
const TONES = [
    { c: 1.9999571633282585,  s1: 0.0065449379673518581, a: 0.012 },   //  50 Hz — the line
    { c: 1.999828655148014,   s1: 0.01308959557134444,   a: 0.006 },   // 100
    { c: 1.9996144809641296,  s1: 0.019633692460628301,  a: 0.003 },   // 150
    { c: 1.9993146499511145,  s1: 0.026176948307873153,  a: 0.0015 },  // 200
];

const frames = Math.round(sr * secs);
const out = new Float32Array(frames * ch);

// One LCG for the bed, so the file is the same bytes wherever it is regenerated.
let st = 0x9E3779B9 >>> 0;
const rnd = () => { st = (Math.imul(st, 1664525) + 1013904223) >>> 0; return st / 4294967296 * 2 - 1; };

// The tone states, advanced once per frame and shared by every channel: hum is common-mode.
const y = TONES.map(t => ({ prev: 0, cur: t.s1 }));

for (let i = 0; i < frames; ++i) {
    let hum = 0;
    for (let k = 0; k < TONES.length; ++k) {
        hum += TONES[k].a * y[k].cur;
        const next = TONES[k].c * y[k].cur - y[k].prev;
        y[k].prev = y[k].cur; y[k].cur = next;
    }
    // Loud 0-6 s and 14-20 s; quiet 6-14 s and 20-30 s. TWO quiet stretches — the detector needs two to
    // call a line stationary — and each is longer than a whole analysis window plus a hop, which is what
    // it takes for a frame to sit ENTIRELY inside one. An earlier draft used 3-second stretches against a
    // 2.73-second window and produced zero quiet frames: the stretch existed and no frame fitted in it.
    const t = i / sr;
    const loud = (t < 6) || (t >= 14 && t < 20);
    const bed = loud ? 0.30 : 0.0004;
    for (let c = 0; c < ch; ++c) out[i * ch + c] = Math.fround(bed * rnd() + hum);
}

writeFileSync(outPath, Buffer.from(out.buffer, out.byteOffset, out.byteLength));
console.error(`${outPath}: ${frames} frames x ${ch} ch @ ${sr} Hz, two quiet stretches with a 50 Hz line`);
