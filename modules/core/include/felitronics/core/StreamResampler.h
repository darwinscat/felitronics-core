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
// felitronics::core::StreamResampler — streaming arbitrary-ratio resampler (polyphase windowed sinc).
// Promoted from OrbitCab, where it runs a rate-locked neural (NAM) model at its native 48 kHz on any
// host rate — host→model on the way in, model→host on the way out. feed() appends input; the
// produce*() calls emit as many output samples as the buffered history allows, with `pos` carrying
// the sub-sample phase across blocks — arbitrary in/out block sizes, no long-term drift.
//
// THE KERNEL. `kTaps = 64` taps of a Kaiser-windowed sinc (β = 8.6), sampled into `kPhases = 512`
// phase rows built once in reset(), with LINEAR interpolation between the two rows bracketing the
// wanted phase. The window is a fixed 64 INPUT samples wide in both directions; what adapts to the
// ratio is the CUTOFF:
//
//     fc = 0.99 · min(1, outRate / inRate)          in units of the INPUT Nyquist
//
// i.e. the kernel always band-limits to 0.99 of whichever Nyquist is lower — its own on the way up,
// the OUTPUT's on the way down, which is what gives the decimating direction a stopband at all. On
// the shipped 44.1 ↔ 48 kHz NAM round trip both legs land on the same 21 829.5 Hz.
//
// 🔴 WHY IT IS NOT THE CATMULL-ROM CUBIC ANY MORE (P32 measured the cubic; P34 replaced it). A
// phase-dependent kernel is a LINEAR PERIODICALLY TIME-VARYING filter, so its error is two things at
// once and they are the same thing seen twice: the per-phase gain M(t) has Fourier coefficients
// H(Ω+2πk), i.e. the "amplitude modulation" and the "interpolation images" are one mechanism. The
// cubic cost, one round trip 44.1↔48 kHz, coherent carrier / worst phase, in dB, was
//     10 k −0.64/−1.16 · 15 k −2.59/−5.14 · 17.64 k −4.17/−9.27 · 20 k −5.48/−14.79
// and its decimating direction had NO stopband at all (−3 dB rms, 0.0 dB sample peak above the output
// Nyquist, because at phase t = 0 its weights were (0,1,0,0) — a bare sample pick). Worse than the
// droop: a DRIVEN nonlinear stage downstream does not mask those images, it DEMODULATES them — a
// 20 kHz tone at −18 dBFS through a high-gain capture came back with a 100 Hz line 14.5 dB LOUDER than
// its own carrier. Full write-up and the protocol: docs/STREAM-RESAMPLER-COST.md.
//
// 🔴 LATENCY — DERIVED FROM THE GEOMETRY, then measured back. Do not guess this number; the previous
// one was a guess ("~3 samples of lookahead per stage") and was 2.16 samples wrong. reset() primes
// `buf` with `kTaps` leading zeros and `pos = kHalf`, so buf[q] holds input sample q − kTaps, output k
// reads centre position kHalf + k·inPerOut, and therefore
//
//     output k is centred on input position  k·inPerOut − kHalf     →  D = kHalf = 32
//
// EVERY stage delays by exactly kHalf (32) of ITS OWN input samples — the group delay of a symmetric
// 64-tap FIR, no more and no less. The 44.1 ↔ 48 round trip is 32 host samples (down) + 32 model
// samples (up) converted to host rate = 32 + 32·(44100/48000) = 61.4000 host samples, 1.392 ms.
// MEASURED back from the carrier phase at 100 Hz and 500 Hz (below the first whole-period wrap, so
// unambiguous): 61.4000 samples, the geometry to four decimals. It was 3.8375 with the cubic.
// The delay is mildly frequency-dependent on top of that — the coherent term's GROUP delay, which no
// single integer can express and no consumer can act on. Callers that must align a dry path (OrbitCab,
// orbit-amp) take this same figure, so it is an audio-alignment number there and not only PDC.
//
// IDENTITY RATIO. At an exactly equal in/out rate the class short-circuits to a pure delay: a 0.99
// cutoff is a real (if gentle) low-pass, and a caller asking for no rate change must not silently get
// one. The delay is still kHalf, so latency stays one formula for every ratio.
//
// FAMILY SPLIT vs convolution::resampleIr: THAT is the OFFLINE Kaiser windowed-sinc (message-thread,
// allocates, its own cutoff criterion) for rate-converting an impulse response on load — an IR is a
// fingerprint and must survive intact. THIS is the streaming rate-match for a live signal path. The
// two are now the same FAMILY of kernel with different budgets; they are still not interchangeable,
// because that one is sized for an offline one-shot and this one for a per-block audio callback.
//
// 🔴 Fixed capacity and the phase table, both allocated once in reset() (message thread):
// feed()/produce*() never allocate, lock, do IO, or throw on the audio thread. `buf` is a linear
// scratch holding `len` valid samples at the front, compacted by memmove; `pos` is the fractional read
// position. Header-only so it can be unit-tested directly.
//==============================================================================
struct StreamResampler
{
    static constexpr int kTaps   = 64;             // window length in INPUT samples (even)
    static constexpr int kHalf   = kTaps / 2;      // 32 — taps at/ahead of the read head, and D
    static constexpr int kBehind = kHalf - 1;      // 31 — taps behind it
    static constexpr int kPhases = 512;            // phase rows; measured error floor −105 dB (−94 at 256)
    static constexpr double kCutoff = 0.99;        // of the lower Nyquist — the product decision (P34)
    static constexpr double kBeta   = 8.6;         // Kaiser β → ≈87 dB stopband at 64 taps
    static constexpr double kPi     = 3.14159265358979323846;   // repo convention: MSVC has no M_PI

    double inPerOut = 1.0;              // input samples advanced per output (= inRate / outRate)
    double pos      = (double) kHalf;   // fractional read position into buf (>= kBehind: needs buf[i-31])
    std::vector<float> buf;             // capacity fixed in reset(); buf[0..len) valid
    int len = 0;
    std::vector<float> tab;             // (kPhases+1) rows x kTaps, normalised per row — built in reset()
    bool identity = false;              // exact 1:1 ratio → pure delay, no filtering (see header note)

    // The delay this stage adds, in ITS OWN INPUT samples — the single source of truth for every
    // consumer that has to align something against it. It is the group delay of the symmetric
    // kTaps-long FIR, and it does NOT depend on the ratio: output k is centred on input k·r − kHalf
    // (derivation in the header block above). A round trip through two stages costs
    // delayInputSamples() of the first stage's input rate plus the same count of the second's, which
    // is why NamStage reports kHalf·(1 + hostSR/modelRunSR) and not twice one number.
    static constexpr double delayInputSamples() noexcept { return (double) kHalf; }

    // Modified Bessel I0, series. Hand-rolled on purpose: std::cyl_bessel_i is not dependably present
    // across the toolchains this repo builds on (MSVC, Apple clang, emscripten), and the table has to
    // come out the same on all of them.
    static double besselI0 (double x) noexcept
    {
        double sum = 1.0, term = 1.0;
        const double h = x * 0.5;
        for (int k = 1; k < 64; ++k)
        {
            term *= (h * h) / ((double) k * (double) k);
            sum  += term;
            if (term < 1.0e-18 * sum) break;
        }
        return sum;
    }

    // One tap of the continuous kernel at argument x (in input samples), cutoff fc (of input Nyquist).
    static double kernelAt (double x, double fc) noexcept
    {
        const double s = (std::fabs (x) < 1.0e-12) ? fc
                                                   : fc * std::sin (kPi * fc * x) / (kPi * fc * x);
        const double u = x / (double) kHalf;
        if (std::fabs (u) >= 1.0) return 0.0;                       // the window is exactly zero at the edge
        return s * besselI0 (kBeta * std::sqrt (1.0 - u * u)) / besselI0 (kBeta);
    }

    void reset (double inRate, double outRate, int capacity)
    {
        // Contract-violation guards (message thread only — no cost on the valid path). A non-finite or
        // non-positive rate would put inf/NaN into inPerOut, and floor(inf) converted to int is UB in
        // the produce loop; fall back to the identity ratio, which is inert rather than undefined.
        if (! (inRate > 0.0) || ! (outRate > 0.0) || ! std::isfinite (inRate) || ! std::isfinite (outRate))
            inRate = outRate = 1.0;

        inPerOut = inRate / outRate;
        // `fabs(d) <= 0` is exact equality without tripping -Wfloat-equal, and (unlike `!(d > 0)`)
        // it is false for a NaN — which the guard above has already excluded anyway.
        identity = (std::fabs (inRate - outRate) <= 0.0);

        if (capacity < 0) capacity = 0;                                          // negative would wrap (size_t) into a huge/UB alloc
        if (capacity > INT_MAX - (kTaps + 8)) capacity = INT_MAX - (kTaps + 8);  // keep buf.size() <= INT_MAX so the (int) buf.size() below never wraps negative
        buf.assign ((std::size_t) capacity + (std::size_t) kTaps + 8, 0.0f);     // ALLOC here (message thread) — never in process
        len = kTaps;                                                             // kTaps leading history zeros
        pos = (double) kHalf;                                                    // → output k centred on input k·r − kHalf

        // The phase table. Row p holds the taps for phase t = p/kPhases; row kPhases is the t → 1 limit,
        // present so the linear interpolation never has to wrap — a wrapped row would need a SHIFTED tap
        // alignment, which is the classic off-by-one of this construction.
        const double fc = kCutoff * std::min (1.0, outRate / inRate);
        tab.assign ((std::size_t) (kPhases + 1) * (std::size_t) kTaps, 0.0f);
        for (int p = 0; p <= kPhases; ++p)
        {
            const double t = (double) p / (double) kPhases;
            double row[kTaps], s = 0.0;
            for (int j = 0; j < kTaps; ++j) { row[j] = kernelAt ((double) (j - kBehind) - t, fc); s += row[j]; }
            // Partition of unity per row: without it the settled DC gain would MODULATE with the phase
            // (the raw sums drift by ~6 ppm here), which is a defect of exactly the kind this kernel
            // exists to remove. Normalising both rows also normalises their linear interpolant.
            const double inv = (std::fabs (s) > 0.0 ? 1.0 / s : 1.0);
            for (int j = 0; j < kTaps; ++j)
                tab[(std::size_t) p * (std::size_t) kTaps + (std::size_t) j] = (float) (row[j] * inv);
        }
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
            pos = (double) kBehind;
        }
        if (n > cap - len)                                   // backstop: sized so this shouldn't trigger on the valid path
        {                                                    // (`len + n > cap` would be signed overflow at a huge capacity)
            int drop = n - (cap - len);
            if (drop > len) drop = len;                      // never memmove a size_t-underflowing (negative) count
            std::memmove (buf.data(), buf.data() + drop, (std::size_t) (len - drop) * sizeof (float));
            len -= drop; pos -= drop;
            if (pos < (double) kBehind) pos = (double) kBehind;   // keep the read head in-bounds — a dropped-past head must not let produce read below buf[0]
        }
        std::copy (in, in + n, buf.data() + len);
        len += n;
    }

    int produceAvailable (float* out, int cap)
    {
        int k = 0;
        while (k < cap)
        {
            const int i = (int) std::floor (pos);
            if (i + kHalf >= len) break;                     // need buf[i-kBehind .. i+kHalf]
            if (i < kBehind) break;                          // …and the history behind it (backstop paths only)

            const float* x = buf.data() + (std::size_t) (i - kBehind);
            if (identity)
            {
                out[k++] = x[kBehind];                       // pure delay: THIS sample, unfiltered
            }
            else
            {
                const double fp = (pos - (double) i) * (double) kPhases;
                int p = (int) fp;
                if (p < 0) p = 0;
                if (p > kPhases - 1) p = kPhases - 1;
                const float a = (float) (fp - (double) p);
                const float* r0 = tab.data() + (std::size_t) p * (std::size_t) kTaps;
                const float* r1 = r0 + kTaps;

                // 🔴 THE SHAPE OF THIS LOOP WAS CHOSEN BY MEASUREMENT, not by which algebra reads better.
                // Blending the two phase rows per TAP and summing once is algebraically the same as two
                // dot products blended by a scalar — sum_j (r0_j + a(r1_j − r0_j))·x_j = d0 + a(d1 − d0)
                // — and the second form looks cheaper: 128 MACs against 64 MACs plus 64 lerps. It is not.
                // Benchmarked on arm64/Apple clang at 64 taps: per-tap lerp 12.77 ns per output (10.0
                // GMAC/s), two accumulators 16.70 ns, a hand-written four-way split 21.21 ns. A float
                // reduction is not associative, so the compiler may not restructure it; ONE dependency
                // chain with a fused lerp vectorises, TWO compete, and splitting it by hand only spends
                // registers. The three are not bit-identical to each other — this is the shipped order.
                float s = 0.0f;
                for (int j = 0; j < kTaps; ++j) s += (r0[j] + a * (r1[j] - r0[j])) * x[j];
                out[k++] = s;
            }
            pos += inPerOut;
        }
        int keep = (int) std::floor (pos) - kBehind;         // keep kBehind samples before the read head
        if (keep > len) keep = len;                          // a huge ratio can step past everything buffered
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
