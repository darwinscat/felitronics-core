// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// Writes the interleaved f32le fixture the native↔wasm parity check runs on, so CI can prove bit-exactness
// on GENERATED audio and nothing private ever enters this public repo.
//
// DETERMINISM BY CONSTRUCTION. Not one transcendental is called: no Math.sin, no Math.pow, no Math.random.
// The signal is an integer xorshift32 PRNG plus float32 arithmetic — every operation exactly specified by
// IEEE-754 and ECMAScript — so the bytes are identical on every JS engine and node version. The fixture is
// reproducible, not merely shared between the two sides of one run. (Sharing one file is all the diff
// strictly needs; reproducibility is what lets a number recorded in a report still mean something later.)
// The output is written little-endian explicitly (see the bottom), so the claim holds on a big-endian host
// too rather than resting on one.
//
// Usage: node make-fixture.mjs <out.f32> [sampleRate=48000] [channels=2] [seconds=10]

import { writeFileSync } from 'node:fs';

const [, , outPath, srArg = '48000', chArg = '2', secArg = '10'] = process.argv;
if (!outPath) { console.error('usage: node make-fixture.mjs <out.f32> [sampleRate] [channels] [seconds]'); process.exit(2); }

const sr = Number(srArg), ch = Number(chArg), seconds = Number(secArg);
const frames = Math.round(sr * seconds);

// xorshift32 — integer only, so every sample below is a pure function of the seed.
function rng (seed) {
    let x = seed >>> 0;
    return () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296 * 2 - 1; };
}

const out = new Float32Array(frames * ch);
for (let c = 0; c < ch; ++c) {
    const rand = rng(0x9e3779b9 + c * 0x85ebca6b);       // a different, fixed seed per channel: decorrelated stereo
    const b = t => Math.round(sr * t);                    // section boundary in frames
    let burst = 0.60;
    for (let i = 0; i < frames; ++i) {
        let v;
        if      (i < b(3))  v = rand() * 0.1;             // ~-20 dBFS noise — ordinary programme, many gated blocks
        else if (i < b(4))  v = 0;                        // DIGITAL SILENCE — the absolute gate, and the K-weighting
                                                          //   filter's subnormal tail (finding F1) on the way through
        else if (i < b(7))  v = rand() * 0.4;             // ~-8 dBFS — blocks well above the relative gate
        else if (i < b(9))  v = ((i & 3) < 2 ? 0.7 : -0.7) + rand() * 0.02;   // fs/4 sampled 45 deg off the crests:
                                                          //   every sample sits at A, the reconstruction peaks at
                                                          //   A*sqrt(2) — so TRUE peak genuinely exceeds sample peak
                                                          //   and the compared number comes from the polyphase FIR
                                                          //   rather than from the sample-peak floor. (An alternating
                                                          //   +-A pattern would NOT do this: at exactly Nyquist the
                                                          //   reconstruction max IS A, measured tp == sp.)
        else if (frames - i > 16) { v = ((i & 3) < 2 ? burst : -burst); burst *= 0.9995; }   // a decaying burst...
        else                v = ((i & 3) < 2 ? 0.99 : -0.99);   // ...and then a full-scale hit occupying the LAST 16
                                                          //   SAMPLES, which simply stops. That is the case that made
                                                          //   fcore::Probe under-report true peak (0.0152 for a true
                                                          //   1.0625) until finish() drained the FIR: the oversampler's
                                                          //   group delay is 63.5 oversampled samples, so the tail of
                                                          //   the file never left the filter. A DECAYING ending does
                                                          //   not test this at all — 0.98*0.9995^48000 is 3.7e-11,
                                                          //   about -209 dBFS, i.e. silence. The burst starts LOWER
                                                          //   than this hit on purpose, so the file's true peak is
                                                          //   decided by the final 16 samples and a missing drain
                                                          //   changes the compared number instead of hiding behind
                                                          //   a louder passage earlier on.
        out[i * ch + c] = v;
    }
}

// Written little-endian EXPLICITLY. A TypedArray's backing bytes are in the host's order, so
// `Buffer.from(out.buffer)` would emit big-endian floats on a big-endian node — and the native tool reads
// f32**le**. Naming the order costs one loop and makes the byte-identity claim unconditional.
const bytes = Buffer.allocUnsafe(out.length * 4);
for (let i = 0; i < out.length; ++i) bytes.writeFloatLE(out[i], i * 4);
writeFileSync(outPath, bytes);
console.log(`${outPath}: ${frames} frames x ${ch} ch @ ${sr} Hz (${out.byteLength} bytes)`);
