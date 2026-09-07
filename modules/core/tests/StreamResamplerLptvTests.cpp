// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THE COST OF A PHASE-DEPENDENT KERNEL — pinned by measurement (P32, docs/STREAM-RESAMPLER-COST.md).
//
// The other three StreamResampler suites cover the CONTRACT (golden values, DSP theory at low f, and
// contract-violation hardening). None of them puts a number on what the kernel does at the top of the
// band, and the theory suite says so out loud: its high-f SNR case is labelled "a design pin", not a
// spec. That gap is exactly what let a one-line justification ("the driven nonlinear stage masks the
// interpolation images") stand in the header for a whole release cycle with nothing behind it.
//
// WHAT THIS SUITE ASSERTS, and why each number is derived rather than fitted:
//
//  1. A phase-dependent kernel is a LINEAR PERIODICALLY TIME-VARYING filter: for a complex exponential
//     the output is the ideal output times a per-sample complex gain M(t) that depends only on the
//     interpolation phase. The class is LINEAR, so running cos and sin separately and combining them
//     as cos + i·sin gives the exact response to the complex exponential — a per-sample complex gain
//     with NO bucketing, NO FFT and no grid of ours to collide with the signal's. (Bucketing by the
//     wrong period is precisely how the round-trip numbers were first published wrong: it averages the
//     very modulation being measured.)
//
//  2. The round-trip gain at 44100/48000 = 147/160 is periodic with EXACTLY 147 output samples, and
//     that is a proof, not an observation: stage 2 advances 160/147 per output, so 147 outputs advance
//     it by exactly 160 — a whole number of stage-1 phase periods, since stage 1 advances 147/160 and
//     repeats every 160. Asserted by comparing the statistics over 147 against those over 23520
//     (= lcm) and by falsifying every proper divisor of 147.
//     🔴 What those 147 pairs are NOT is the full 160x147 product: the cascade traces a LINE through
//     the two stages' phase torus, fixed by the shipped priming. The worst-phase column is still
//     effectively a ceiling (0.02 dB off the worst over all 160 alignments at 17.64 kHz) but the
//     coherent carrier is an interference term that belongs to the priming — see the CLOSURE group.
//
//  3. The DECIMATING direction has no stopband. At phase t = 0 the Catmull-Rom weights are (0,1,0,0) —
//     a bare sample pick, which attenuates nothing at any frequency — so a tone above the output
//     Nyquist survives at −3 dB rms and 0 dB sample peak (an ENVELOPE figure: no single line exceeds
//     −4.67 dB) and folds back as TWO components, 44100−g and g−3900, i.e. across 18.15–22.05 kHz.
//     Asserted as a FACT ABOUT THE SHIPPED KERNEL, so that a future kernel with a real stopband fails
//     loudly and has to update the document instead of silently inheriting its claims.
//
// A liveness precondition runs FIRST: the same instrument on the unity ratio must read 0.00 dB where
// the header promises transparency. An instrument that cannot read zero cannot be trusted to read −4.

#include <felitronics/core/StreamResampler.h>

#include "felitronics_test.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

using felitronics::core::StreamResampler;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
    constexpr double kPi = 3.14159265358979323846;   // repo convention: MSVC has no M_PI

    // The shipped NAM round trip, with NamStage::configureRates' own capacities and call pattern.
    struct RoundTrip
    {
        StreamResampler down, up;
        std::vector<float> mid;
        int maxBlock, maxModelFrames;

        RoundTrip (double hostSR, double modelSR, int block)
        {
            maxBlock       = block;
            maxModelFrames = (int) std::ceil (maxBlock * (modelSR / std::max (8000.0, hostSR))) + 16;
            down.reset (hostSR, modelSR, maxBlock * 2 + 16);
            up  .reset (modelSR, hostSR, maxModelFrames * 2 + 16);
            mid.assign ((std::size_t) maxModelFrames, 0.0f);
        }
        void process (float* io, int n)
        {
            down.feed (io, n);
            const int m = down.produceAvailable (mid.data(), maxModelFrames);
            if (m > 0) up.feed (mid.data(), m);
            up.produceExact (io, n);
        }
    };

    struct Gains { std::vector<std::complex<double>> g; };

    // Per-sample complex gain of the round trip at frequency f. cos and sin are pushed through TWO
    // independent instances and combined; the class is linear, so the pair is the complex exponential.
    Gains roundTripGains (double f, double hostSR, double modelSR, int block, int skip, int keep)
    {
        RoundTrip rc (hostSR, modelSR, block), rs (hostSR, modelSR, block);
        const double W = 2.0 * kPi * f / hostSR;
        const int total = ((skip + keep + block) / block) * block;
        Gains out;
        out.g.reserve ((std::size_t) total);
        std::vector<float> bc ((std::size_t) block), bs ((std::size_t) block);
        for (int off = 0; off < total; off += block)
        {
            for (int i = 0; i < block; ++i)
            {
                bc[(std::size_t) i] = (float) std::cos (W * (off + i));
                bs[(std::size_t) i] = (float) std::sin (W * (off + i));
            }
            rc.process (bc.data(), block);
            rs.process (bs.data(), block);
            for (int i = 0; i < block; ++i)
            {
                const int m = off + i;
                if (m >= skip && (int) out.g.size() < keep)
                    out.g.push_back (std::complex<double> (bc[(std::size_t) i], bs[(std::size_t) i])
                                     * std::polar (1.0, -W * m));
            }
        }
        return out;
    }

    double dbOf (double x) { return 20.0 * std::log10 (std::max (x, 1e-30)); }

    double coherentDb (const Gains& g)
    {
        std::complex<double> s (0.0, 0.0);
        for (auto v : g.g) s += v;
        return dbOf (std::abs (s / (double) g.g.size()));
    }
    double worstDb (const Gains& g)
    {
        double lo = 1e300;
        for (auto v : g.g) lo = std::min (lo, std::abs (v));
        return dbOf (lo);
    }
    double bestDb (const Gains& g)
    {
        double hi = 0.0;
        for (auto v : g.g) hi = std::max (hi, std::abs (v));
        return dbOf (hi);
    }
} // namespace

int main()
{
    const double H = 44100.0, M = 48000.0;
    const int kPeriod = 147;                    // proved below, not assumed
    const int kSkip = 20000, kKeep = kPeriod * 200;

    // ------------------------------------------------------------------------------------------------
    group ("LIVENESS — the instrument reads zero where the header promises transparency");
    {
        const auto unity = roundTripGains (17640.0, 48000.0, 48000.0, 64, kSkip, kKeep);
        approx (coherentDb (unity), 0.0, 0.001, "unity ratio at 17.64 kHz: coherent gain is 0.000 dB");
        approx (worstDb (unity),    0.0, 0.001, "unity ratio at 17.64 kHz: worst phase is 0.000 dB too");
        const auto low = roundTripGains (1000.0, H, M, 64, kSkip, kKeep);
        approx (coherentDb (low), 0.0, 0.01, "44.1<->48 at 1 kHz reads 0.00 dB — the droop is not an offset");
        ok (worstDb (low) > -0.01, "…and there is no modulation to speak of at 1 kHz either");

        // 🔴 SIGN. Every statistic above is a magnitude, so a kernel that inverted the signal would
        // pass all of them — a crew mutation that negated every fed sample did exactly that, and the
        // block-size comparison cannot catch it either, because both sides of that comparison run the
        // same (mutated) code. Only an ABSOLUTE reference sees it. At 100 Hz the round-trip delay of
        // 3.84 samples costs just 0.055 rad, so the complex coherent gain must sit at very nearly +1.
        {
            const auto vlow = roundTripGains (100.0, H, M, 64, kSkip, kKeep);
            std::complex<double> sum (0.0, 0.0);
            for (auto v : vlow.g) sum += v;
            const std::complex<double> c = sum / (double) vlow.g.size();
            std::printf ("      complex coherent gain at 100 Hz: %+.4f %+.4fi\n", c.real(), c.imag());
            ok (c.real() > 0.99, "the round trip is POSITIVE unity at 100 Hz (Re = " + std::to_string (c.real())
                                 + ") — a polarity inversion would read -1 and every magnitude test would miss it");
            ok (std::abs (c.imag()) < 0.06, "…and its phase is just the 3.84-sample delay, nothing else");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("CLOSURE — the round-trip gain has period exactly 147 output samples, for the SHIPPED priming");
    {
        // stage 2 advances 160/147 per output -> 147 outputs advance it by exactly 160, which is a whole
        // number of stage-1 phase periods (stage 1 advances 147/160 and repeats every 160).
        //
        // 🔴 WHAT THIS IS AND IS NOT. The cascade does NOT visit all 160x147 phase pairs: the stage-1
        // phase at the point stage 2 reads is a deterministic function of the stage-2 phase, so the pairs
        // it visits are a LINE through that torus, fixed by the two reset() offsets. The 147 values below
        // are therefore complete and exact FOR THE SHIPPED PRIMING (pos = 1, len = 3, both stages reset
        // together in NamStage::configureRates), not a kernel-wide ceiling. Swept over all 160 integer
        // alignments the coherent carrier at 17.64 kHz spans -3.59 .. -6.83 dB and the worst phase
        // -6.96 .. -9.29 dB — so the shipped -9.27 is 0.02 dB off the alignment-wide worst, while the
        // shipped -4.17 is one point of a wide span. A crew round caught this being over-claimed.
        const auto g = roundTripGains (17640.0, H, M, 64, kSkip, kPeriod * 5);
        double worstDev = 0.0;
        for (std::size_t m = 0; m < g.g.size(); ++m)
            worstDev = std::max (worstDev, std::abs (g.g[m] - g.g[m % (std::size_t) kPeriod]));
        ok (worstDev < 1e-5, "g[m] == g[m mod 147] across five periods (float; worst dev "
                             + std::to_string (worstDev) + ")");
        bool everyDivisorFails = true;
        for (int d : { 1, 3, 7, 21, 49 })
        {
            double dev = 0.0;
            for (int m = 0; m < kPeriod; ++m) dev = std::max (dev, std::abs (g.g[(std::size_t) m] - g.g[(std::size_t) (m % d)]));
            everyDivisorFails = everyDivisorFails && dev > 0.1;
        }
        ok (everyDivisorFails, "no proper divisor of 147 is a period — 147 is MINIMAL, not merely a period");

        const auto shortRun = roundTripGains (17640.0, H, M, 64, kSkip, kPeriod);
        const auto longRun  = roundTripGains (17640.0, H, M, 64, kSkip, 23520);   // lcm(147,160)
        approx (worstDb (shortRun), worstDb (longRun), 0.001,
                "the worst phase over ONE period equals the worst over the full lcm — the ceiling is reached");
        approx (coherentDb (shortRun), coherentDb (longRun), 0.001,
                "…and so does the coherent carrier");
    }

    // ------------------------------------------------------------------------------------------------
    group ("THE TABLE — round-trip carrier and worst phase, 44.1 <-> 48 kHz (docs/STREAM-RESAMPLER-COST.md)");
    {
        struct Row { double f, coherent, worst, best; };
        // Derived independently from the closed form M(t) = sum_j w_j(t) e^{jW(o_j - t)} enumerated over
        // the exact rational phases, then confirmed by this instrument and by a brute-force simulation.
        const Row rows[] = {
            {  5000.0, -0.05,  -0.08, -0.01 },
            {  8000.0, -0.28,  -0.50, -0.05 },
            { 10000.0, -0.64,  -1.16, -0.11 },
            { 12000.0, -1.22,  -2.28, -0.21 },
            { 15000.0, -2.59,  -5.14, -0.41 },
            { 17640.0, -4.17,  -9.27, -0.61 },
            { 19000.0, -4.98, -12.20, -0.69 },
            { 20000.0, -5.48, -14.79, -0.74 },
        };
        for (const auto& r : rows)
        {
            const auto g = roundTripGains (r.f, H, M, 64, kSkip, kKeep);
            std::printf ("      %6.0f Hz  coherent %7.2f  worst %7.2f  best %7.2f\n",
                         r.f, coherentDb (g), worstDb (g), bestDb (g));
            approx (coherentDb (g), r.coherent, 0.02, std::to_string ((int) r.f) + " Hz: coherent carrier");
            approx (worstDb (g),    r.worst,    0.02, std::to_string ((int) r.f) + " Hz: worst-phase ceiling");
            approx (bestDb (g),     r.best,     0.02, std::to_string ((int) r.f) + " Hz: best-phase ceiling");
            ok (bestDb (g) <= 0.0001, std::to_string ((int) r.f) + " Hz: NO phase exceeds unity — "
                                      "sum|w_j(t)| = 1 + t(1-t) <= 1.25 bounds one pass at +1.94 dB, but the "
                                      "coherent kernel cannot amplify a tone at all, so a positive reading "
                                      "means a broken instrument (that is how a +5.7 dB oracle bug was found)");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("BLOCK SIZE — the cost is a property of the kernel, not of the caller's buffer");
    {
        // 🔴 Compare the COMPLEX coherent gain, not its magnitude. A crew mutation that flipped polarity
        // survived a |mean| comparison, and the first attempt to close that added the worst phase —
        // which is min|g|, invariant under g -> -g just as |mean| is, so it closed nothing. The diff
        // round caught that too. The complex mean carries the sign and the phase, so it does not.
        auto coherent = [] (const Gains& g)
        {
            std::complex<double> s (0.0, 0.0);
            for (auto v : g.g) s += v;
            return s / (double) g.g.size();
        };
        const auto ref = roundTripGains (17640.0, H, M, 64, kSkip, kKeep);
        const std::complex<double> cref = coherent (ref);
        for (int b : { 1, 17, 63, 128, 512 })
        {
            const auto g = roundTripGains (17640.0, H, M, b, kSkip, kKeep);
            const std::complex<double> c = coherent (g);
            ok (std::abs (c - cref) < 3.0e-3,
                "block " + std::to_string (b) + ": the COMPLEX carrier matches block 64 (|delta| = "
                + std::to_string (std::abs (c - cref)) + ") — magnitude alone would not see a sign flip");
            approx (worstDb (g), worstDb (ref), 0.02, "block " + std::to_string (b) + ": …and so does the worst phase");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("COUNT — the EXACT number of outputs a feed yields, not a bound on it");
    {
        // The theory suite already bounds the produced-count drift at <= 2 samples, which is the right
        // statement about long-term accumulator drift and the wrong one about the priming condition:
        // a mutation that stops the loop one sample early (i + 3 >= len instead of i + 2 >= len) stays
        // inside that bound at every length and survives. It is not harmless — one fewer sample available
        // per call is one more silent sample produceExact() pads at the head, i.e. a whole sample of
        // extra latency. So assert the count EXACTLY, from the derivation:
        //   after reset() len = 3, pos = 1; feed(N) -> len = N + 3; produce emits while
        //   floor(1 + k·r) + 2 < N + 3, i.e. k < N/r  ->  K = ceil(N·outRate/inRate).
        //
        // ONE SUBTLETY, and it is arithmetic rather than a defect: when N·outRate/inRate is an exact
        // integer the last k sits precisely ON the boundary, and `pos` is a double accumulation of a
        // ratio that is not a binary fraction (48000/44100 = 160/147 is not), so it can land a few ULP
        // under the boundary and emit one extra sample. The first version of this test allowed +1 at
        // every integral boundary, which the diff round pointed out is unearned where the ratio IS a
        // binary fraction (96000/48000 = 2 accumulates exactly). So the expectation is not a closed
        // form with slack — it is the SAME double accumulation the class does, which has no slack at
        // all and still fails a loop that stops one sample early.
        struct RB { double in, out; };
        bool allExact = true;
        for (const RB rb : { RB {44100,48000}, RB {48000,44100}, RB {96000,48000}, RB {96000,44100},
                             RB {22050,48000}, RB {48000,96000} })
            for (int N : { 1, 64, 147, 160, 1000, 4096 })
            {
                StreamResampler r;
                r.reset (rb.in, rb.out, N + 16);
                std::vector<float> x ((std::size_t) N, 0.25f), o ((std::size_t) (N * 8 + 64));
                r.feed (x.data(), N);
                const int k = r.produceAvailable (o.data(), (int) o.size());
                // predict with the class's own accumulator: emit while floor(pos) <= N, pos starts at
                // 1.0 and advances by the same double inPerOut. No compaction happens inside a single
                // produce call at this capacity, so a plain accumulation is exact.
                int want = 0;
                for (double pos = 1.0; std::floor (pos) <= (double) N; pos += rb.in / rb.out) ++want;
                if (k != want)
                {
                    allExact = false;
                    std::printf ("      %6.0f -> %6.0f  N=%d produced %d, derived %d\n", rb.in, rb.out, N, k, want);
                }
            }
        ok (allExact, "one feed of N yields EXACTLY the count the class's own phase accumulator implies,"
                      " at every ratio and length — no slack anywhere");

        // SIGN, at the level of ONE stage. The round-trip sign check in LIVENESS cannot see a kernel
        // that inverts, because both stages invert and the two cancel — which is precisely why a crew
        // mutation that negated every fed sample survived every suite. A single stage cannot cancel
        // with itself: a settled constant must come out as ITSELF (the weights are a partition of
        // unity at every phase), sign included.
        {
            StreamResampler r;
            r.reset (48000.0, 44100.0, 4096);
            std::vector<float> x (2048, 0.25f), o (4096);
            r.feed (x.data(), 2048);
            const int k = r.produceAvailable (o.data(), (int) o.size());
            double worstDev = 0.0;
            for (int i = 16; i < k - 16; ++i) worstDev = std::max (worstDev, std::fabs ((double) o[(std::size_t) i] - 0.25));
            ok (k > 1800 && worstDev < 1.0e-6,
                "one stage passes a settled +0.25 as +0.25 (worst deviation " + std::to_string (worstDev)
                + ") — partition of unity at every phase, sign included");
        }

        // produceExact()'s documented behaviour on a startup underflow is "pad with silence", and NOTHING
        // in the four suites checked it: a mutation that pads with the PREVIOUS sample instead survived
        // every one of them. Silence and hold are audibly different at a stream start (a held sample is a
        // DC step into the model), and this is the call NamStage uses on the output leg.
        {
            StreamResampler r;
            r.reset (44100.0, 48000.0, 512);
            std::vector<float> out (64, 0.5f);
            r.produceExact (out.data(), 64);
            // `!(fabs(v) > 0)` would also be true for a NaN — the diff round caught that going in with
            // the -Wfloat-equal cleanup. `fabs(v) <= 0` is true ONLY for an exact zero (fabs is never
            // negative) and false for NaN, which is the predicate this test actually wants.
            bool allZero = true;
            for (float v : out) allZero = allZero && (std::fabs (v) <= 0.0f);
            ok (allZero, "produceExact on a freshly reset resampler writes SILENCE, not held samples");

            StreamResampler r2;
            r2.reset (44100.0, 48000.0, 512);
            std::vector<float> in (8, 0.75f), out2 (200, -1.0f);
            r2.feed (in.data(), 8);
            r2.produceExact (out2.data(), 200);          // 8 in -> at most 9 out, the rest is padding
            int lastNonZero = -1;
            for (int i = 0; i < 200; ++i) if (! (std::fabs (out2[(std::size_t) i]) <= 0.0f)) lastNonZero = i;
            bool tailSilent = true;
            for (int i = lastNonZero + 1; i < 200; ++i)
                tailSilent = tailSilent && (std::fabs (out2[(std::size_t) i]) <= 0.0f);
            ok (lastNonZero < 12 && tailSilent,
                "…and the tail past what the history can produce is silence too (last non-zero at "
                + std::to_string (lastNonZero) + " of 200)");
        }

        // NEGATIVE KNOWLEDGE, so nobody re-derives it: the first moment of the round-trip impulse
        // response is NOT a usable delay oracle here. A single impulse is not DC, and a decimating stage
        // does not preserve it — measured DC sums of 0.999 (44.1 kHz), 0.743 (88.2) and 2.000 (96), with
        // moments of 3.824 / 4.773 / 6.000 against geometries of 3.8375 / 5.6750 / 6.0000. The 88.2 kHz
        // row is off by 0.9 samples. The delay is measured from the CARRIER PHASE instead — see the
        // "latency is a MEASUREMENT" group in modules/nam/tests/NamStageTests.cpp.
    }

    // ------------------------------------------------------------------------------------------------
    group ("DECIMATION — the 48 -> 44.1 stage has NO stopband: 0 dB peak above the output Nyquist");
    {
        // At phase t = 0 the weights are (0,1,0,0): a bare sample pick, which cannot attenuate anything.
        approx ((double) StreamResampler::catmull (0.0f, 1.0f, 0.0f, 0.0f, 0.0f), 1.0, 1e-7,
                "catmull at t=0 IS the centre sample — the phase that gives the decimation 0 dB");
        for (double f : { 22100.0, 22500.0, 23000.0, 23500.0, 23900.0 })
        {
            StreamResampler r;
            r.reset (M, H, 8192);
            const int block = 86;
            double peak = 0.0, sq = 0.0; long n = 0;
            std::vector<float> b ((std::size_t) block), o (8192);
            for (int off = 0; off < 200000; off += block)
            {
                for (int i = 0; i < block; ++i) b[(std::size_t) i] = (float) std::sin (2.0 * kPi * f * (off + i) / M);
                r.feed (b.data(), block);
                const int got = r.produceAvailable (o.data(), (int) o.size());
                if (off > 40000)
                    for (int k = 0; k < got; ++k)
                    { const double v = o[(std::size_t) k]; peak = std::max (peak, std::abs (v)); sq += v * v; ++n; }
            }
            const double peakDb = dbOf (peak);
            const double rmsDb  = dbOf (std::sqrt (sq / (double) n) * std::sqrt (2.0));
            std::printf ("      %6.0f Hz above the 22.05 kHz Nyquist -> rms %6.2f dB  peak %6.2f dB\n",
                         f, rmsDb, peakDb);
            ok (peakDb > -0.1, std::to_string ((int) f) + " Hz survives the decimation at FULL level "
                               "(sample peak " + std::to_string (peakDb) + " dB) and folds in. The peak is a "
                               "TIME-domain figure and no spectral line exceeds -4.67 dB: the phase returns to "
                               "within ~1e-13 of t=0 every 147 outputs, and there the output simply IS an input "
                               "sample -- measured 0.999657 against an input peak of exactly 1.0");
            approx (rmsDb, -3.05, 0.25, std::to_string ((int) f) + " Hz: ~3 dB of rms rejection, and that is all there is");
        }

        // WHERE it lands, which is not one place. The first write-up of this said "20.1-22.05 kHz",
        // counting only the k = 0 fold; a crew round found the second component, and it is only 1.5-2.5 dB
        // down. A tone at g in (22.05, 24) kHz returns as 44100 - g AND as g - 3900 (the image at
        // 48000 - g, folded), so the damage covers 18.15-22.05 kHz.
        {
            const double f = 23000.0;
            StreamResampler r;
            r.reset (M, H, 8192);
            const int block = 96, N = (int) H * 2;              // 2 s -> a coherent window at 100 Hz spacing
            std::vector<float> b ((std::size_t) block), o (8192), y;
            y.reserve ((std::size_t) N + 30000);
            for (int off = 0; (int) y.size() < N + 20000; off += block)
            {
                for (int i = 0; i < block; ++i) b[(std::size_t) i] = (float) std::sin (2.0 * kPi * f * (off + i) / M);
                r.feed (b.data(), block);
                const int k = r.produceAvailable (o.data(), (int) o.size());
                y.insert (y.end(), o.begin(), o.begin() + k);
            }
            auto lineDb = [&y, N, H] (double a)
            {
                std::complex<double> acc (0.0, 0.0);
                for (int n = 10000; n < 10000 + N; ++n)
                    acc += (double) y[(std::size_t) n] * std::polar (1.0, -2.0 * kPi * a * n / H);
                return dbOf (2.0 * std::abs (acc) / (double) N);
            };
            const double lo = lineDb (H - f);                   // 21100 Hz — the k = 0 fold
            const double hi = lineDb (f - 3900.0);              // 19100 Hz — the image at 48000 - f, folded
            std::printf ("      23 kHz in -> 21100 Hz %.2f dB and 19100 Hz %.2f dB\n", lo, hi);
            approx (lo, -5.33, 0.15, "23 kHz folds to 21.1 kHz at -5.33 dB");
            approx (hi, -7.04, 0.15, "…and to 19.1 kHz at -7.04 dB — the SECOND component, 1.7 dB down, "
                                     "which is why the damaged band is 18.15-22.05 kHz and not just the top slice");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("THE CRITERION — the rate-match against a DRIVEN nonlinearity's OWN aliasing floor");
    {
        // The header used to justify this kernel with "the driven nonlinear stage masks the interpolation
        // images". Made measurable: the rate-match must not add MORE than the nonlinearity folds down by
        // itself, because that floor is what a product running a network at 48 kHz already accepts. Two
        // quantities, neither containing the other, both read off the same tone.
        //
        // The nonlinearity here is a declared stand-in (tanh at a stated drive), not a capture: a capture
        // cannot be shipped in a test, and the real ones are in docs/STREAM-RESAMPLER-COST.md. The point
        // of a stand-in is that its drive is a KNOB, so the claim can be tested across the whole range the
        // word "driven" covers instead of at one preset.
        //
        // Both windows are coherent, so Parseval gives the non-carrier energy exactly, in O(N), with no
        // window function and no leakage: total power minus the carrier's.
        auto nonCarrierDbc = [] (const std::vector<double>& y, int off, int N, double fs, double f0)
        {
            std::complex<double> acc (0.0, 0.0);
            double tot = 0.0;
            for (int n = 0; n < N; ++n)
            {
                const double v = y[(std::size_t) (off + n)];
                tot += v * v;
                acc += v * std::polar (1.0, -2.0 * kPi * f0 * (off + n) / fs);
            }
            const double amp = 2.0 * std::abs (acc) / (double) N;
            const double car = 0.5 * amp * amp;
            return 10.0 * std::log10 (std::max (tot / (double) N - car, 1e-300) / std::max (car, 1e-300));
        };

        const double f0 = 17500.0;                    // on the 63-sample grid at 44.1 k and 96 at 48 k
        const int N = 441 * 200, skip = 20000, block = 64;
        RoundTrip rt (H, M, block);
        std::vector<float> b ((std::size_t) block);
        std::vector<double> y;
        y.reserve ((std::size_t) (skip + N + block));
        for (int o = 0; o < skip + N + block; o += block)
        {
            for (int i = 0; i < block; ++i) b[(std::size_t) i] = (float) std::sin (2.0 * kPi * f0 * (o + i) / H);
            rt.process (b.data(), block);
            for (int i = 0; i < block; ++i) y.push_back ((double) b[(std::size_t) i]);
        }
        const double added = nonCarrierDbc (y, skip, N, H, f0);
        std::printf ("      rate-match alone at 17.5 kHz adds %.2f dBc\n", added);
        approx (added, -8.84, 0.05, "the round trip adds -8.84 dBc of artifacts at 17.5 kHz");

        // …and how hard the nonlinearity has to be driven before its OWN folding reaches that level.
        // The crossing is BISECTED rather than read off the grid: the first version of this test
        // asserted "the first tested drive where the floor wins", which the diff round correctly called
        // an artefact of the grid — the crossing is between 4 and 8, not at 8.
        struct D { double drive, floorDbc; };
        const D drives[] = { D {0.25, -45.80}, D {1.0, -23.47}, D {4.0, -10.44}, D {8.0, -8.24} };
        for (int di = 0; di < 4; ++di)
        {
            const D d = drives[di];
            const int NF = 96 * 1000;                 // 96 = the tone's period at 48 kHz; coherent
            std::vector<double> z ((std::size_t) NF);
            for (int n = 0; n < NF; ++n)
                z[(std::size_t) n] = std::tanh (d.drive * std::sin (2.0 * kPi * f0 * n / M)) / std::tanh (d.drive);
            const double fl = nonCarrierDbc (z, 0, NF, M, f0);
            std::printf ("      tanh drive %5.2f folds %.2f dBc by itself -> rate-match is %+.2f dB against it\n",
                         d.drive, fl, added - fl);
            approx (fl, d.floorDbc, 0.05, "tanh at drive " + std::to_string (d.drive) + " folds its own "
                    + std::to_string (d.floorDbc) + " dBc");
        }
        {
            auto floorAt = [&] (double drive)
            {
                const int NF = 96 * 1000;
                std::vector<double> z ((std::size_t) NF);
                for (int n = 0; n < NF; ++n)
                    z[(std::size_t) n] = std::tanh (drive * std::sin (2.0 * kPi * f0 * n / M)) / std::tanh (drive);
                return nonCarrierDbc (z, 0, NF, M, f0);
            };
            double lo = 1.0, hi = 16.0;
            for (int i = 0; i < 24; ++i) { const double m = 0.5 * (lo + hi); if (floorAt (m) < added) lo = m; else hi = m; }
            const double cross = 0.5 * (lo + hi);
            std::printf ("      crossing: the tanh has to reach drive %.2f before it folds as much as the "
                         "rate-match adds\n", cross);
            approx (cross, 6.23, 0.05,
                    "the nonlinearity has to be driven to tanh(6.2x) — the peak of a unit sine already sits "
                    "at 0.99999, i.e. a near square wave — before its OWN aliasing reaches what the "
                    "rate-match adds. Below that the rate-match is the LOUDER artifact, which is the "
                    "opposite of masking");
        }
    }

    return felitronics::test::report();
}
