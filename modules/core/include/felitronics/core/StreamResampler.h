// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

namespace felitronics::core
{

//==============================================================================
// felitronics::core::StreamResampler — streaming arbitrary-ratio resampler (cubic / Catmull-Rom).
// Promoted from OrbitCab, where it runs a rate-locked neural (NAM) model at its native 48 kHz on any
// host rate — host→model on the way in, model→host on the way out. feed() appends input; the
// produce*() calls emit as many output samples as the buffered history allows, with `pos` carrying
// the sub-sample phase across blocks — arbitrary in/out block sizes, no long-term drift. Identity
// ratio passes the signal with a clean 2-sample delay (catmull @ t=0 reads one sample behind), and
// EVERY ratio delays by exactly 2 of this stage's own input samples (buf holds 3 leading zeros and
// pos starts at 1, so output k reads input position k·inPerOut − 2).
//
// 🔴 WHAT THIS KERNEL COSTS — MEASURED (P32), not asserted. A phase-dependent kernel is a LINEAR
// PERIODICALLY TIME-VARYING filter, so its error is two things at once, and they are the same thing
// seen twice: the per-phase gain M(t) has Fourier coefficients H(Ω+2πk), i.e. the "amplitude
// modulation" and the "interpolation images" are one mechanism, not two. Round trip 44.1↔48 kHz
// (the shipped NAM path), coherent carrier / worst phase, in dB:
//     10 k −0.64/−1.16 · 15 k −2.59/−5.14 · 17.64 k −4.17/−9.27 · 20 k −5.48/−14.79
// The composite gain is periodic with EXACTLY 147 output samples at this ratio and visits every
// phase pair, so those are ceilings from the topology, not sampled maxima.
//
// 🔴 AND THE DECIMATING DIRECTION HAS NO ANTI-ALIASING AT ALL. Going 48 → 44.1 this kernel passes a
// tone above the output Nyquist at −3 dB rms and 0.0 dB PEAK — at phase t = 0 the weights are
// (0,1,0,0), a bare sample pick, which attenuates nothing at any frequency. Everything the driven
// stage makes between 22.05 and 24 kHz folds into 20.1–22.05 kHz.
//
// The older note here said "the driven nonlinear stage masks the interpolation images". Measured, it
// is CONDITIONAL and it does not cover the whole error: (a) the carrier droop is not an added
// component at all, so nothing masks it; (b) against the model's OWN aliasing floor the added
// artifacts sit BELOW it on a high-gain capture and up to +9.8 dB ABOVE it on a clean one, from
// 15 kHz up; (c) driving the stage harder makes it WORSE, not better — over a 42 dB input sweep the
// error grew +11.4 dB and the folded-alias term +12.5 dB, because "driven" is exactly what fills
// 22.05–24 kHz for the un-filtered decimation to fold back. Whether to change the kernel is a
// product decision (it buys transparency with latency: a 32-tap polyphase sinc measures
// −0.11 dB flat at 17.64 kHz with no modulation and 16–42 dB of stopband, for +0.066 %RT per mono
// channel and 30.7 host samples of round-trip delay against today's 3.84). Full numbers, the
// two-oracle protocol and the candidate table: docs/STREAM-RESAMPLER-COST.md.
//
// FAMILY SPLIT vs convolution::resampleIr: THAT is the OFFLINE Kaiser windowed-sinc (≥60 dB-class,
// message-thread, allocates) for rate-converting an impulse response on load — an IR is a
// fingerprint and must survive intact. THIS is the cheap STREAMING rate-match for a live signal
// path. The split still holds — but it is a COST decision, not a transparency claim, and the cost
// above is the price. Don't swap them.
//
// 🔴 Fixed capacity, allocated once in reset() (message thread): feed()/produce*() never allocate,
// lock, do IO, or throw on the audio thread. `buf` is a linear scratch holding `len` valid samples
// at the front, compacted by memmove; `pos` is the fractional read position. Header-only so it can
// be unit-tested directly.
//==============================================================================
struct StreamResampler
{
    double inPerOut = 1.0;        // input samples advanced per output (= inRate / outRate)
    double pos      = 1.0;        // fractional read position into buf (>=1: cubic needs buf[i-1])
    std::vector<float> buf;       // capacity fixed in reset(); buf[0..len) valid
    int len = 0;

    void reset (double inRate, double outRate, int capacity)
    {
        inPerOut = inRate / outRate;
        // Contract-violation guards (message thread only — no cost on the valid path):
        if (capacity < 0) capacity = 0;                      // negative would wrap (size_t) into a huge/UB alloc (or, small negatives, a tiny buffer)
        if (capacity > INT_MAX - 8) capacity = INT_MAX - 8;  // keep buf.size() <= INT_MAX so the (int) buf.size() below never wraps negative
        buf.assign ((std::size_t) capacity + 8, 0.0f);       // ALLOC here (message thread) — never in process
        len = 3;                                             // 3 leading history zeros
        pos = 1.0;
    }

    void feed (const float* in, int n)
    {
        if (n <= 0) return;                                  // guard: n==0 is a no-op (bit-identical); reject negative (std::copy of a reversed range is UB)
        const int cap = (int) buf.size();                    // reset() keeps buf.size() <= INT_MAX, so this never wraps
        if (n > cap)                                         // BACKSTOP: block alone can't fit — keep only its newest `cap` samples, drop all history
        {
            in += (n - cap);
            n   = cap;
            len = 0;
            pos = 1.0;
        }
        if (len + n > cap)                                   // backstop: sized so this shouldn't trigger on the valid path
        {
            int drop = len + n - cap;
            if (drop > len) drop = len;                      // never memmove a size_t-underflowing (negative) count
            std::memmove (buf.data(), buf.data() + drop, (std::size_t) (len - drop) * sizeof (float));
            len -= drop; pos -= drop;
            if (pos < 1.0) pos = 1.0;                        // keep the read head in-bounds — a dropped-past head must not let produce read buf[i-1] below buf[0]
        }
        std::copy (in, in + n, buf.data() + len);
        len += n;
    }

    static float catmull (float a, float b, float c, float d, float t)
    {
        const float t2 = t * t, t3 = t2 * t;
        return 0.5f * ((2.0f * b) + (-a + c) * t
                     + (2.0f * a - 5.0f * b + 4.0f * c - d) * t2
                     + (-a + 3.0f * b - 3.0f * c + d) * t3);
    }

    int produceAvailable (float* out, int cap)
    {
        int k = 0;
        while (k < cap)
        {
            const int i = (int) std::floor (pos);
            if (i + 2 >= len) break;                         // need buf[i-1..i+2]
            out[k++] = catmull (buf[(std::size_t) (i - 1)], buf[(std::size_t) i],
                                buf[(std::size_t) (i + 1)], buf[(std::size_t) (i + 2)],
                                (float) (pos - i));
            pos += inPerOut;
        }
        const int keep = (int) std::floor (pos) - 1;         // keep one sample before the read head
        if (keep > 0)
        {
            std::memmove (buf.data(), buf.data() + keep, (std::size_t) (len - keep) * sizeof (float));
            len -= keep; pos -= keep;
        }
        return k;
    }

    void produceExact (float* out, int want)                 // pad with silence on startup underflow
    {
        const int got = produceAvailable (out, want);
        for (int k = got; k < want; ++k) out[k] = 0.0f;
    }
};

} // namespace felitronics::core
