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
        const double ref = coherentDb (roundTripGains (17640.0, H, M, 64, kSkip, kKeep));
        for (int b : { 1, 17, 63, 128, 512 })
        {
            const auto g = roundTripGains (17640.0, H, M, b, kSkip, kKeep);
            approx (coherentDb (g), ref, 0.02, "block " + std::to_string (b) + ": same carrier as block 64");
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
        // ONE EXCEPTION, and it is arithmetic rather than a defect: when N·outRate/inRate is an exact
        // integer the last k sits precisely ON the boundary, and `pos` is a double accumulation of a
        // ratio that is not a binary fraction (48000/44100 = 160/147 is not), so the accumulated value
        // can land a few ULP under the boundary and emit one extra sample. Measured on exactly the two
        // rows where the product is integral: 48000->44100 at N=160 (148 against 147) and 22050->48000
        // at N=147 (321 against 320). So: exact everywhere, +1 allowed ONLY on that boundary — which
        // still fails the one-sample-early mutation on every other row.
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
                const double exact = (double) N * rb.out / rb.in;
                const int want = (int) std::ceil (exact);
                const bool onBoundary = std::fabs (exact - std::floor (exact + 0.5)) < 1e-9;
                if (! (k == want || (onBoundary && k == want + 1)))
                {
                    allExact = false;
                    std::printf ("      %6.0f -> %6.0f  N=%d produced %d, derived %d\n", rb.in, rb.out, N, k, want);
                }
            }
        ok (allExact, "one feed of N yields EXACTLY ceil(N·outRate/inRate) outputs at every ratio and length"
                      " (+1 only where the product is integral and `pos` lands an ULP under the boundary)");

        // produceExact()'s documented behaviour on a startup underflow is "pad with silence", and NOTHING
        // in the four suites checked it: a mutation that pads with the PREVIOUS sample instead survived
        // every one of them. Silence and hold are audibly different at a stream start (a held sample is a
        // DC step into the model), and this is the call NamStage uses on the output leg.
        {
            StreamResampler r;
            r.reset (44100.0, 48000.0, 512);
            std::vector<float> out (64, 0.5f);
            r.produceExact (out.data(), 64);
            bool allZero = true;
            for (float v : out) allZero = allZero && ! (std::fabs (v) > 0.0f);
            ok (allZero, "produceExact on a freshly reset resampler writes SILENCE, not held samples");

            StreamResampler r2;
            r2.reset (44100.0, 48000.0, 512);
            std::vector<float> in (8, 0.75f), out2 (200, -1.0f);
            r2.feed (in.data(), 8);
            r2.produceExact (out2.data(), 200);          // 8 in -> at most 9 out, the rest is padding
            int lastNonZero = -1;
            for (int i = 0; i < 200; ++i) if (std::fabs (out2[(std::size_t) i]) > 0.0f) lastNonZero = i;
            bool tailSilent = true;
            for (int i = lastNonZero + 1; i < 200; ++i)
                tailSilent = tailSilent && ! (std::fabs (out2[(std::size_t) i]) > 0.0f);
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
                               "(sample peak " + std::to_string (peakDb) + " dB — an ENVELOPE figure: no single "
                               "spectral line exceeds -4.67 dB, the peak is where the phases align) and folds in");
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
        struct D { double drive, floorDbc; };
        const D drives[] = { D {0.25, -45.80}, D {1.0, -23.47}, D {4.0, -10.44}, D {8.0, -8.24} };
        int firstOver = -1;
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
            if (fl > added && firstOver < 0) firstOver = di;
        }
        ok (firstOver == 3, "the nonlinearity has to be driven all the way to tanh(8x) — a near square "
                              "wave — before its OWN aliasing reaches what the rate-match adds; below that "
                              "the rate-match is the LOUDER artifact, which is the opposite of masking");
    }

    return felitronics::test::report();
}
