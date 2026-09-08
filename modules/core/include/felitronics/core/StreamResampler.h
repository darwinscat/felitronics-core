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
// EVERY stage delays by exactly kHalf (32) of ITS OWN input samples. 🔴 NOT "the group delay of a
// symmetric 64-tap FIR" — that would be 31.5, and saying so was wrong even though the number is right.
// The taps sit at offsets −31 … +32 RELATIVE TO `pos`, and the kernel is symmetric about argument
// zero, i.e. about `pos` itself; the 32 is where the priming puts `pos`, not a property of the tap
// count. The 44.1 ↔ 48 round trip is 32 host samples (down) + 32 model samples (up) converted to host
// rate = 32 + 32·(44100/48000) = 61.4000 host samples, 1.392 ms.
// MEASURED back from the carrier phase at 100 Hz and 200 Hz — 500 Hz is NOT usable, its period is
// 88.2 samples so 61.4 wraps and reads −26.8, which is the same delay and a different branch. Result:
// 61.4000, the geometry to four decimals. It was 3.8375 with the cubic.
// The phase delay is FLAT with frequency: swept on one stage at 50 Hz … 19 kHz on both legs it reads
// 32.000000 in every cell. (The cubic's was not, and its header carried a +0.62-sample-at-20-kHz
// caveat; that caveat died with it and is deliberately not restated.) Callers that must align a dry
// path (OrbitCab, orbit-amp, and rigplayer's own dry/wet blend) take this same figure, so it is an
// audio-alignment number there and not only PDC.
//
// IDENTITY RATIO. At an exactly equal in/out rate the class short-circuits to a pure delay: a 0.99
// cutoff is a real (if gentle) low-pass, and a caller asking for no rate change must not silently get
// one. The delay is still kHalf, so latency stays one formula for every ratio, and no phase table is
// built at all (513 rows x 64 floats = 128 KiB that a bit-copy never reads).
// ⚠️ THAT IS A DISCONTINUITY IN THE TRANSFER FUNCTION, AND IT IS EXACT. Rates that differ by one ULP
// take the filtered path; rates that are bit-identical take the copy. A caller that rate-matches
// "nearly equal" clocks (a drifting external device, a 48000.001 host) will therefore hear the band
// edge appear the moment the two stop being equal. That is the right trade for THIS class — the
// alternative is a tolerance nobody can pick — but it is a sharp edge and callers should know.
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

    // The delay ONE stage adds, in ITS OWN INPUT samples — the single source of truth for every
    // consumer that has to align something against it. It is NOT the (N−1)/2 = 31.5 of a symmetric
    // 64-tap FIR: output k is centred on input k·r − kHalf (derivation in the header block above).
    //
    // 🔴 IT TAKES THE RATES, AND TODAY IT IGNORES THEM. That is deliberate and it is the whole point
    // of the signature. The kernel is a FIXED kTaps wide right now, so the delay is kHalf whatever the
    // conversion — but the open plan item (docs/STREAM-RESAMPLER-COST.md §6.7) is to scale kTaps by
    // max(1, inRate/outRate), which makes the delay a function of the ratio and makes the two legs of
    // a round trip DIFFERENT: 192 kHz → 48 kHz would cost 128 input samples going down and 32 model
    // samples coming back, not 32 and 32. Every consumer that had baked "the delay is a constant" into
    // its own arithmetic would then be wrong, silently, exactly as one downstream repository already
    // was. Taking the rates now costs nothing and means that change edits one body.
    static double delayInputSamples ([[maybe_unused]] double inRate,
                                     [[maybe_unused]] double outRate) noexcept
    {
        return (double) kHalf;
    }

    // What a DOWN+UP PAIR costs, in HOST samples — the composition, owned here so nobody restates it.
    //
    // 🔴 THIS IS GEOMETRY, NOT LATENCY. It answers "what would a pair of these cost", and at equal
    // rates it answers 2·kHalf, because a pair really would cost that. Whether a pair is INSTALLED at
    // all is a policy question belonging to the consumer — nam::NamStage installs one only past a
    // 0.5 Hz difference and reports 0 below it, and that gate lives with the policy, in
    // NamStage::rateMatch(). Reading this number as "the latency" is the mistake this split exists to
    // make impossible.
    //
    // The expression order is the shipped one and is kept deliberately: `a + a·h/m`, not the tidier
    // `a·(1 + h/m)`. They agree bit-for-bit on every rate pair measured (60 pairs, zero differences)
    // for a reason that will EXPIRE — kHalf is a power of two, so scaling by it is exact and both
    // spellings carry a single rounding. Make the delay ratio-dependent (the open kTaps item) and D
    // becomes 35, 59, 118…, where integer straddles do exist: the nearest to the audio grid is a
    // 6930 Hz host against a 44.1 kHz model. One spelling, in one place, is the whole defence.
    static double pairDelayHostSamples (double hostSR, double modelRunSR) noexcept
    {
        const double down = delayInputSamples (hostSR, modelRunSR);      // host samples, going down
        const double up   = delayInputSamples (modelRunSR, hostSR);      // MODEL samples, coming back
        // 🔴 THE GROUPING IS THE SHIPPED ONE, LITERALLY. `up * hostSR / modelRunSR` parses as
        // `(up * hostSR) / modelRunSR`, which is what NamStage computed before this extraction
        // (`d + d * hostSR / modelRunSR`). Writing `up * (hostSR / modelRunSR)` is the same value on
        // every audio rate — 60 pairs, zero differences — and is NOT the same expression: it rounds
        // the quotient first, and at finite extremes the two can part. A crew round caught the
        // regrouping; "no number moves" has to mean the arithmetic, not just the answers we sampled.
        return down + up * hostSR / modelRunSR;                          // …converted to host samples
    }

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
        // 🔴 GUARD THE QUOTIENT, NOT ONLY THE OPERANDS. Both rates can be finite, positive and sane on
        // their own and still produce a step this class cannot walk: reset(DBL_MAX, 1, ...) passes every
        // check above and then converts DBL_MAX to int inside the produce loop, which is UB; a step of
        // DBL_MIN/DBL_MAX is zero, so `pos` never advances and one input sample is emitted forever. The
        // bound is deliberately generous — 1e6:1 in either direction is far past any audio use — and it
        // keeps floor(pos) inside int for any buffer this class can hold.
        if (! std::isfinite (inPerOut) || inPerOut <= 0.0 || inPerOut > 1.0e6 || inPerOut < 1.0e-6)
            inPerOut = 1.0;
        // `fabs(d) <= 0` is exact equality without tripping -Wfloat-equal, and (unlike `!(d > 0)`)
        // it is false for a NaN — which the guard above has already excluded anyway.
        identity = (std::fabs (inRate - outRate) <= 0.0);

        if (capacity < 0) capacity = 0;                                          // negative would wrap (size_t) into a huge/UB alloc
        if (capacity > INT_MAX - (kTaps + 8)) capacity = INT_MAX - (kTaps + 8);  // keep buf.size() <= INT_MAX so the (int) buf.size() below never wraps negative

        // 🔴 ORDER MATTERS ON THE FAILURE PATH. len goes to 0 FIRST, so that if either allocation below
        // throws (bad_alloc, and this function is not noexcept), the object is left in the one state
        // that produces nothing rather than one that reads an empty table: produceAvailable's first
        // test is `i + kHalf >= len`, which is true at len = 0 and breaks immediately. NamStage's
        // prepare catches and marks itself unprepared, but a standalone caller may catch and carry on.
        len = 0;
        pos = (double) kHalf;
        buf.assign ((std::size_t) capacity + (std::size_t) kTaps + 8, 0.0f);     // ALLOC here (message thread) — never in process

        // The phase table. Row p holds the taps for phase t = p/kPhases; row kPhases is evaluated at
        // t = 1 and exists so the linear interpolation never has to wrap — a wrapped row would need a
        // SHIFTED tap alignment, which is the classic off-by-one of this construction. (It is exactly
        // that shift: row kPhases[j] = row 0[j-1] by construction, with a fresh tap entering at j = 0.)
        // 🔴 It is NOT "the continuous kernel's t → 1 limit", which is what an earlier version of this
        // comment claimed. kernelAt() forces the window to zero at |u| >= 1, whereas a true Kaiser has
        // a non-zero endpoint I0(0)/I0(beta) there; the coefficient this drops is about -1.12e-5, i.e.
        // -99 dB. So this is a ZERO-ENDED Kaiser variant, deliberately, and the seam is consistent with
        // it — but the two are not the same function and the comment must not say they are.
        if (identity)
        {
            tab.clear();                    // the copy path never reads it — do not pay 128 KiB for it
            tab.shrink_to_fit();
            len = kTaps;                    // kTaps leading history zeros; the delay is kHalf either way
            return;
        }

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
        len = kTaps;                                                             // kTaps leading history zeros
    }

    // noexcept on the hot path is the house contract, and it is honest here: every call below is
    // on floats, raw pointers or memmove. reset() is deliberately NOT noexcept — it allocates.
    void feed (const float* in, int n) noexcept
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

    int produceAvailable (float* out, int cap) noexcept
    {
        int k = 0;
        while (k < cap)
        {
            const int i = (int) std::floor (pos);
            if (i >= len - kHalf) break;                     // need buf[i-kBehind .. i+kHalf]; written this
                                                             // way so a huge i cannot overflow i + kHalf
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

    void produceExact (float* out, int want) noexcept        // pad with silence on startup underflow
    {
        const int got = produceAvailable (out, want);
        for (int k = got; k < want; ++k) out[k] = 0.0f;
    }
};

} // namespace felitronics::core
