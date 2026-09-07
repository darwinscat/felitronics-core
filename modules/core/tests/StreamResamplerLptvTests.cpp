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
//     that is a closure proof, not an observation: stage 2 advances 160/147 per output, so 147 outputs
//     advance it by exactly 160 — a whole number of stage-1 phase periods, since stage 1 advances
//     147/160 and repeats every 160. Every phase pair the topology can make therefore occurs inside
//     one period, so min/max over it are CEILINGS. Asserted by comparing the statistics over 147
//     against those over 23520 (= lcm) and by falsifying every proper divisor of 147.
//
//  3. The DECIMATING direction has no stopband. At phase t = 0 the Catmull-Rom weights are (0,1,0,0) —
//     a bare sample pick, which attenuates nothing at any frequency — so a tone above the output
//     Nyquist survives at 0 dB peak and folds. This is asserted as a FACT ABOUT THE SHIPPED KERNEL, so
//     that a future kernel with a real stopband fails the assertion loudly and has to update the
//     document instead of silently inheriting its claims.
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
    }

    // ------------------------------------------------------------------------------------------------
    group ("CLOSURE — the round-trip gain has period exactly 147 output samples, so min/max are CEILINGS");
    {
        // stage 2 advances 160/147 per output -> 147 outputs advance it by exactly 160, which is a whole
        // number of stage-1 phase periods (stage 1 advances 147/160 and repeats every 160).
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
        struct Row { double f, coherent, worst; };
        // Derived independently from the closed form M(t) = sum_j w_j(t) e^{jW(o_j - t)} enumerated over
        // the exact rational phases, then confirmed by this instrument and by a brute-force simulation.
        const Row rows[] = {
            {  5000.0, -0.05,  -0.08 },
            {  8000.0, -0.28,  -0.50 },
            { 10000.0, -0.64,  -1.16 },
            { 12000.0, -1.22,  -2.28 },
            { 15000.0, -2.59,  -5.14 },
            { 17640.0, -4.17,  -9.27 },
            { 19000.0, -4.98, -12.20 },
            { 20000.0, -5.48, -14.79 },
        };
        for (const auto& r : rows)
        {
            const auto g = roundTripGains (r.f, H, M, 64, kSkip, kKeep);
            std::printf ("      %6.0f Hz  coherent %7.2f  worst %7.2f  best %7.2f\n",
                         r.f, coherentDb (g), worstDb (g), bestDb (g));
            approx (coherentDb (g), r.coherent, 0.02, std::to_string ((int) r.f) + " Hz: coherent carrier");
            approx (worstDb (g),    r.worst,    0.02, std::to_string ((int) r.f) + " Hz: worst-phase ceiling");
            ok (bestDb (g) <= 0.0001, std::to_string ((int) r.f) + " Hz: no phase EXCEEDS unity "
                                      "(a gain above 0 dB would mean the oracle, not the kernel, is wrong)");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("BLOCK SIZE — the cost is a property of the kernel, not of the caller's buffer");
    {
        const double ref = coherentDb (roundTripGains (17640.0, H, M, 64, kSkip, kKeep));
        for (int b : { 1, 17, 63, 128, 512 })
        {
            const auto g = roundTripGains (17640.0, H, M, b, kSkip, kKeep);
            approx (coherentDb (g), ref, 0.02, "block " + std::to_string (b) + ": same carrier as block 64");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("DECIMATION — the 48 -> 44.1 stage has NO stopband: 0 dB peak above the output Nyquist");
    {
        // At phase t = 0 the weights are (0,1,0,0): a bare sample pick, which cannot attenuate anything.
        approx ((double) StreamResampler::catmull (0.0f, 1.0f, 0.0f, 0.0f, 0.0f), 1.0, 1e-7,
                "catmull at t=0 IS the centre sample — the phase that gives the decimation 0 dB");
        for (double f : { 22100.0, 23000.0, 23900.0 })
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
                               "(peak " + std::to_string (peakDb) + " dB) and folds into the audio band");
            approx (rmsDb, -3.05, 0.25, std::to_string ((int) f) + " Hz: ~3 dB of rms rejection, and that is all there is");
        }
    }

    return felitronics::test::report();
}
