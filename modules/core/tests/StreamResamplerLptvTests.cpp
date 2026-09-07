// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// WHAT A PHASE-DEPENDENT KERNEL COSTS — pinned by measurement. P32 measured the Catmull-Rom cubic
// that used to live here; P34 replaced it with a 64-tap polyphase windowed sinc, and this suite is the
// before/after (docs/STREAM-RESAMPLER-COST.md).
//
// 🔴 THE TRAP THIS SUITE IS BUILT AROUND, and it is specific to a kernel that is now TRANSPARENT.
// P32's instrument proved itself by reading 0.00 dB on the unity ratio: an instrument that cannot read
// zero cannot be trusted to read −4. That argument DIES the moment the kernel under test is itself
// flat, because then "reads 0.00" is what a completely broken instrument reads too — a wrong cutoff, a
// wrong tap centre, a mirrored table, an oracle bucketing by the wrong period: every one of them
// reports a beautiful 0.00 dB on a flat kernel. So liveness here runs the OTHER way as well: the same
// instrument, unchanged, must first read the OLD kernel's −4.17 / −9.27 at 17.64 kHz off a reference
// Catmull-Rom copy kept in this file. Only an instrument that still sees the defect may certify it gone.
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
//     repeats every 160. 🔴 The MINIMALITY half of that proof (no proper divisor of 147 is also a
//     period) CANNOT be run on the new kernel: on a flat kernel every divisor is a period to within
//     the noise, so the assertion would pass vacuously. It is run on the reference cubic, where it has
//     a subject, and the flat result is asserted the other way round for the record.
//
//  3. The decimating direction HAS a stopband now: −9.08 dB at 22.1 kHz down to −88.77 at 23.9 kHz,
//     against a flat −3 dB (and 0.0 dB sample PEAK) before. The 0 dB peak is gone with the sample-pick
//     phase that caused it.
//
//  4. The kernel is transparent in the passband but NOT infinitely so, and the tolerances say which.
//     A Kaiser β = 8.6 window has a real passband ripple δ = 10^(−86.7/20) = 4.6e−5, i.e. ±0.0004 dB
//     per stage and ±0.0008 dB round trip. Asserting "≤ 0.0001 dB" would FAIL a correct kernel — the
//     measured best phase is +0.00039 dB at 17.64 kHz. A criterion that fails the right answer is worse
//     than no criterion.

#include <felitronics/core/StreamResampler.h>

#include "felitronics_test.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <vector>

using felitronics::core::StreamResampler;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
    constexpr double kPi = 3.14159265358979323846;   // repo convention: MSVC has no M_PI

    //==============================================================================================
    // THE REFERENCE KERNEL — a copy of the Catmull-Rom class this module shipped up to v0.26.0, kept
    // ONLY so the instrument can be shown to still see the defect that was fixed. It is never the
    // subject of an acceptance claim; it is the calibration weight on the scale.
    //==============================================================================================
    struct CatmullRef
    {
        double inPerOut = 1.0, pos = 1.0;
        std::vector<float> buf;
        int len = 0;

        void reset (double inRate, double outRate, int capacity)
        {
            inPerOut = inRate / outRate;
            buf.assign ((std::size_t) capacity + 8, 0.0f);
            len = 3;
            pos = 1.0;
        }
        void feed (const float* in, int n)
        {
            if (n <= 0) return;
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
                if (i + 2 >= len) break;
                out[k++] = catmull (buf[(std::size_t) (i - 1)], buf[(std::size_t) i],
                                    buf[(std::size_t) (i + 1)], buf[(std::size_t) (i + 2)],
                                    (float) (pos - i));
                pos += inPerOut;
            }
            const int keep = (int) std::floor (pos) - 1;
            if (keep > 0)
            {
                std::memmove (buf.data(), buf.data() + keep, (std::size_t) (len - keep) * sizeof (float));
                len -= keep; pos -= keep;
            }
            return k;
        }
        void produceExact (float* out, int want)
        {
            const int got = produceAvailable (out, want);
            for (int k = got; k < want; ++k) out[k] = 0.0f;
        }
    };

    // The shipped NAM round trip, with NamStage::configureRates' own capacities and call pattern.
    // Templated so the SAME instrument drives the shipped kernel and the reference cubic.
    template <typename R>
    struct RoundTrip
    {
        R down, up;
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
    template <typename R>
    Gains roundTripGains (double f, double hostSR, double modelSR, int block, int skip, int keep)
    {
        RoundTrip<R> rc (hostSR, modelSR, block), rs (hostSR, modelSR, block);
        const double W = 2.0 * kPi * f / hostSR;
        Gains out;
        out.g.reserve ((std::size_t) keep);
        std::vector<float> bc ((std::size_t) block), bs ((std::size_t) block);
        for (int off = 0; (int) out.g.size() < keep; off += block)
        {
            for (int i = 0; i < block; ++i)
            {
                bc[(std::size_t) i] = (float) std::cos (W * (off + i));
                bs[(std::size_t) i] = (float) std::sin (W * (off + i));
            }
            rc.process (bc.data(), block);
            rs.process (bs.data(), block);
            for (int i = 0; i < block && (int) out.g.size() < keep; ++i)
            {
                const int m = off + i;
                if (m >= skip)
                    out.g.push_back (std::complex<double> (bc[(std::size_t) i], bs[(std::size_t) i])
                                     * std::polar (1.0, -W * m));
            }
        }
        return out;
    }

    double dbOf (double x) { return 20.0 * std::log10 (std::max (x, 1e-30)); }

    std::complex<double> coherent (const Gains& g)
    {
        std::complex<double> s (0.0, 0.0);
        for (auto v : g.g) s += v;
        return s / (double) g.g.size();
    }
    double coherentDb (const Gains& g) { return dbOf (std::abs (coherent (g))); }
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

    // Non-carrier energy relative to the carrier, over a COHERENT window: Parseval gives it exactly, in
    // O(N), with no window function and no leakage — total power minus the carrier's.
    double nonCarrierDbc (const std::vector<double>& y, int off, int N, double fs, double f0)
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
    }

    template <typename R>
    double addedDbc (double H, double M, double f0, int block, int skip, int N)
    {
        RoundTrip<R> rt (H, M, block);
        std::vector<float> b ((std::size_t) block);
        std::vector<double> y;
        y.reserve ((std::size_t) (skip + N + block));
        for (int o = 0; o < skip + N + block; o += block)
        {
            for (int i = 0; i < block; ++i) b[(std::size_t) i] = (float) std::sin (2.0 * kPi * f0 * (o + i) / H);
            rt.process (b.data(), block);
            for (int i = 0; i < block; ++i) y.push_back ((double) b[(std::size_t) i]);
        }
        return nonCarrierDbc (y, skip, N, H, f0);
    }
} // namespace

int main()
{
    const double H = 44100.0, M = 48000.0;
    const int kPeriod = 147;                    // proved below, not assumed
    const int kSkip = 20000, kKeep = kPeriod * 200;

    // The Kaiser passband ripple, DERIVED, and the tolerance every "flat" assertion below leans on.
    // β = 8.6 → stopband A = β/0.1102 + 8.7 = 86.7 dB → δ = 10^(−A/20) = 4.6e−5 → 20·log10(1+δ) =
    // 0.0004 dB per stage, 0.0008 dB for a round trip through two.
    const double kRipple = 20.0 * std::log10 (1.0 + std::pow (10.0, -86.7 / 20.0)) * 2.0;

    // ------------------------------------------------------------------------------------------------
    group ("INSTRUMENT LIVENESS — it must read ZERO on a transparent path AND still read the OLD defect");
    {
        // Direction 1, inherited: the unity ratio must read 0.00 where the class promises transparency.
        {
            const auto g = roundTripGains<StreamResampler> (17640.0, 48000.0, 48000.0, 64, 2000, kPeriod * 20);
            approx (coherentDb (g), 0.0, 1.0e-6, "unity ratio: the instrument reads 0.00 dB");
            approx (worstDb (g),    0.0, 1.0e-6, "…at every phase, because there is no phase");
        }

        // 🔴 Direction 2, and this is the one that matters now the kernel is flat. THE SAME instrument,
        // not a variant of it, driven with the reference cubic. If it cannot reproduce P32's published
        // −4.17 / −9.27 at 17.64 kHz, its 0.000 dB on the shipped kernel means nothing whatever.
        {
            const auto g = roundTripGains<CatmullRef> (17640.0, H, M, 64, kSkip, kKeep);
            std::printf ("      reference cubic at 17.64 kHz: coherent %.2f  worst %.2f  best %.2f\n",
                         coherentDb (g), worstDb (g), bestDb (g));
            approx (coherentDb (g), -4.17, 0.02, "the instrument still reads the OLD kernel's carrier droop");
            approx (worstDb (g),    -9.27, 0.02, "…and its worst phase — so a 0.00 dB reading below is EARNED");
            ok (bestDb (g) - worstDb (g) > 8.0,
                "…and it still resolves 8.7 dB of phase modulation, which is the quantity that has to "
                "disappear; an instrument averaging over the wrong period would report a narrow spread "
                "here and then a comfortable zero afterwards");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("BOTH AXES — the whole point of P34, before and after in one table");
    {
        // Coherent carrier / worst phase / best phase, one round trip. The BEFORE column is P32's
        // published measurement, REPRODUCED here by the reference cubic rather than quoted; the AFTER
        // column is the shipped kernel. Tolerances on the after column are the derived Kaiser ripple,
        // not a fitted slack.
        struct Row { double f, oldCoh, oldWorst, newCoh, newWorst, newBest; };
        const Row rows[] = {
            {  5000.0, -0.05,  -0.08, -0.000128, -0.000392, +0.000034 },
            {  8000.0, -0.28,  -0.50, -0.000170, -0.000300, -0.000055 },
            { 10000.0, -0.64,  -1.16, -0.000054, -0.000175, +0.000025 },
            { 12000.0, -1.22,  -2.28, -0.000041, -0.000181, +0.000028 },
            { 15000.0, -2.59,  -5.14, +0.000003, -0.000053, +0.000113 },
            { 17640.0, -4.17,  -9.27, +0.000160, +0.000016, +0.000385 },
            { 19000.0, -4.98, -12.20, +0.000261, +0.000151, +0.000445 },
            { 20000.0, -5.48, -14.79, -0.013301, -0.013469, -0.013113 },
        };
        std::printf ("      f, Hz |    was: coh / worst  |    now: coh / worst / best\n");
        for (const auto& r : rows)
        {
            const auto oldG = roundTripGains<CatmullRef>      (r.f, H, M, 64, kSkip, kKeep);
            const auto newG = roundTripGains<StreamResampler> (r.f, H, M, 64, kSkip, kKeep);
            std::printf ("      %6.0f | %7.2f / %7.2f      | %+9.6f / %+9.6f / %+9.6f\n",
                         r.f, coherentDb (oldG), worstDb (oldG),
                         coherentDb (newG), worstDb (newG), bestDb (newG));

            approx (coherentDb (oldG), r.oldCoh,   0.02, std::to_string ((int) r.f) + " Hz: the cubic's carrier, as P32 published it");
            approx (worstDb    (oldG), r.oldWorst, 0.02, std::to_string ((int) r.f) + " Hz: …and its worst phase");

            approx (coherentDb (newG), r.newCoh,   0.002, std::to_string ((int) r.f) + " Hz: the sinc's carrier");
            approx (worstDb    (newG), r.newWorst, 0.002, std::to_string ((int) r.f) + " Hz: …its worst phase");
            approx (bestDb     (newG), r.newBest,  0.002, std::to_string ((int) r.f) + " Hz: …and its best phase");

            // 🔴 THE ACTUAL CLAIM: the two phase extremes have COLLAPSED ONTO each other. The carrier
            // number alone would be satisfied by a kernel that modulates symmetrically about unity.
            ok (bestDb (newG) - worstDb (newG) <= kRipple + 0.0005,
                std::to_string ((int) r.f) + " Hz: modulation depth is now "
                + std::to_string (bestDb (newG) - worstDb (newG)) + " dB, inside the derived Kaiser "
                "ripple of " + std::to_string (kRipple) + " dB — it was up to 8.7 dB");

            // …and nothing EXCEEDS unity by more than that same derived ripple. A positive reading
            // beyond it means a broken instrument (that is how a +5.7 dB oracle bug was once found) —
            // but the bound must be the REAL ripple, not zero: this kernel legitimately reads +0.00039,
            // and the assertion it replaced (bestDb <= 0.0001) would have failed the correct answer.
            ok (bestDb (newG) <= kRipple,
                std::to_string ((int) r.f) + " Hz: no phase exceeds the derived passband ripple");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("HOST RATE — where 'transparent' holds, and the rate at which it stops");
    {
        // 🔴 TRANSPARENCY IS NOT A PROPERTY OF THE KERNEL ALONE, and every other number in this suite is
        // measured at one host rate. The window is a fixed kTaps INPUT samples, so its transition width
        // in Hz scales with the INPUT rate: at a 192 kHz host the down leg is a 4:1 decimation and those
        // same 64 taps buy a ~17 kHz transition, which starts eating the audio band. A crew round found
        // this by sweeping rates; nothing here would have.
        //
        // The rows below are the honest scope of the claim, and the table is asserted BOTH ways: flat
        // where it is flat, and DROOPING where it droops. Pinning only the good rates would be the same
        // omission again, one level up.
        struct Row { double host; double at20k; double tol; const char* note; };
        const Row rows[] = {
            {  44100.0, -0.0133, 0.002, "the shipped NAM rate: flat to the band edge" },
            {  88200.0, -0.0009, 0.002, "still flat" },
            {  96000.0, -0.0075, 0.002, "still flat" },
            { 176400.0, -0.6147, 0.010, "🔴 NOT flat: a 3.675:1 decimation with a fixed 64 taps" },
            { 192000.0, -0.7908, 0.010, "🔴 NOT flat: a 4:1 decimation with a fixed 64 taps" },
        };
        for (const auto& r : rows)
        {
            const auto g = roundTripGains<StreamResampler> (20000.0, r.host, M, 64, 20000, 20000);
            std::printf ("      host %6.0f Hz, 20 kHz round trip: %+8.4f dB   %s\n",
                         r.host, coherentDb (g), r.note);
            approx (coherentDb (g), r.at20k, r.tol,
                    std::to_string ((int) r.host) + " Hz host: the 20 kHz carrier is what the fixed "
                    "64-tap window buys at THIS ratio");
        }
        // …and the claim stated as a claim, so a future kernel that fixed it would have to come and
        // edit this line rather than silently pass: transparency is asserted for hosts up to 96 kHz.
        {
            const auto g96 = roundTripGains<StreamResampler> (19000.0, 96000.0, M, 64, 20000, 20000);
            ok (std::fabs (coherentDb (g96)) < 0.002,
                "at 96 kHz the round trip is still flat at 19 kHz (" + std::to_string (coherentDb (g96))
                + " dB) — the transparency claim covers 44.1 / 88.2 / 96, and stops there");
            const auto g192 = roundTripGains<StreamResampler> (19000.0, 192000.0, M, 64, 20000, 20000);
            ok (coherentDb (g192) < -0.1,
                "…and at 192 kHz it is NOT (" + std::to_string (coherentDb (g192)) + " dB at 19 kHz), "
                "which is a documented limit of a fixed-length window under a 4:1 decimation, not a bug "
                "to be tuned away here — scaling kTaps with the ratio is a design change with a CPU "
                "cost proportional to it");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("SIGN AND DELAY — magnitudes cannot see a polarity flip, and this kernel is flat");
    {
        // Every statistic above is a magnitude, so a kernel that inverted the signal would pass all of
        // them — a crew mutation that negated every fed sample did exactly that. Only an ABSOLUTE
        // reference sees it. 100 Hz is chosen because 61.4 samples is well under one whole period there
        // (441 samples), so the phase is the delay with no wrap to resolve: at 500 Hz the same
        // measurement reads −26.8, which is 61.4 − 88.2 and equally true.
        const auto g = roundTripGains<StreamResampler> (100.0, H, M, 64, kSkip, kKeep);
        const std::complex<double> c = coherent (g);
        const double W = 2.0 * kPi * 100.0 / H;
        ok (c.real() > 0.6, "the round trip does NOT invert (Re = " + std::to_string (c.real()) + ")");
        approx (std::abs (c), 1.0, 1.0e-4, "…and it does not change the level at 100 Hz");
        approx (-std::arg (c) / W, 61.4, 0.01,
                "the carrier phase back-counts the delay to 61.4 samples — the geometry "
                "kHalf·(1 + 44100/48000) = 32·1.91875, to four decimals");

        const auto g2 = roundTripGains<StreamResampler> (200.0, H, M, 64, kSkip, kKeep);
        approx (-std::arg (coherent (g2)) / (2.0 * kPi * 200.0 / H), 61.4, 0.01,
                "…and the same delay at 200 Hz, i.e. the phase delay is FLAT with frequency — the "
                "cubic's was not, and its header carried a +0.62-sample-at-20-kHz caveat because of it");

        // 🔴 …AND FLATNESS ASSERTED ONLY AT 100 AND 200 Hz IS FLATNESS ASSERTED WHERE IT IS TRIVIAL.
        // The cubic's own frequency dependence was +0.018 samples at 10 kHz and +0.620 at 20 kHz — i.e.
        // entirely in the top octave, exactly where these two frequencies say nothing. The header
        // claims the new kernel is flat; that claim belongs up there. Above ~360 Hz the delay exceeds
        // half a tone period and the phase wraps, so each reading is resolved into the branch nearest
        // the geometry — which is legitimate BECAUSE the integer is independently fixed by the impulse
        // onset in NamStageTests, not because we like 61.4.
        for (double f : { 5000.0, 10000.0, 17640.0, 19000.0 })
        {
            const auto gh = roundTripGains<StreamResampler> (f, H, M, 64, kSkip, kKeep);
            const double Wf = 2.0 * kPi * f / H;
            const double period = 2.0 * kPi / Wf;
            double d = -std::arg (coherent (gh)) / Wf;
            while (d < 61.4 - 0.5 * period) d += period;
            while (d > 61.4 + 0.5 * period) d -= period;
            approx (d, 61.4, 0.01, std::to_string ((int) f) + " Hz: the delay is STILL 61.4 samples — "
                                   "flat into the top octave, which is the only place the cubic's was not");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("CLOSURE — the 147-sample period, and why minimality can only be proved on the cubic");
    {
        // stage 2 advances 160/147 per output -> 147 outputs advance it by exactly 160, which is a whole
        // number of stage-1 phase periods (stage 1 advances 147/160 and repeats every 160).
        const auto shortRun = roundTripGains<CatmullRef> (17640.0, H, M, 64, kSkip, kPeriod);
        const auto longRun  = roundTripGains<CatmullRef> (17640.0, H, M, 64, kSkip, 23520);   // lcm(147,160)
        // Tolerance 1e-4 dB and not 1e-6, for a reason worth writing down rather than re-deriving: the
        // period is EXACT in arithmetic, so the two runs differ only by the rounding of their own
        // stimulus. The long run evaluates cos(W·m) out to m ≈ 43 500, where the argument reduction
        // costs ulps the short run never pays; the residual measured here is 1e-5 dB. That is the
        // instrument's stimulus, not a period violation — a real violation is 8.7 dB away.
        approx (worstDb (shortRun), worstDb (longRun), 1.0e-4, "147 outputs already contain the worst phase");
        approx (bestDb  (shortRun), bestDb  (longRun), 1.0e-4, "…and the best");
        approx (coherentDb (shortRun), coherentDb (longRun), 1.0e-4, "…and the coherent carrier");

        // 🔴 MINIMALITY, on the cubic, because on a flat kernel it is vacuous. Every proper divisor of
        // 147 (1, 3, 7, 21, 49) must FAIL to be a period.
        for (int d : { 1, 3, 7, 21, 49 })
        {
            double dev = 0.0;
            for (int i = 0; i + d < (int) shortRun.g.size(); ++i)
                dev = std::max (dev, std::abs (shortRun.g[(std::size_t) i] - shortRun.g[(std::size_t) (i + d)]));
            ok (dev > 0.1, "…and " + std::to_string (d) + " is NOT a period of it (max deviation "
                           + std::to_string (dev) + ") — so 147 is minimal, not merely sufficient");
        }

        // On the shipped kernel the SAME measurement is asserted the other way round: every divisor IS
        // a period to within the noise, which is exactly what "no modulation" means and exactly why the
        // minimality proof above has to live on the cubic.
        const auto flat = roundTripGains<StreamResampler> (17640.0, H, M, 64, kSkip, kPeriod);
        double flatDev = 0.0;
        for (int i = 0; i + 1 < (int) flat.g.size(); ++i)
            flatDev = std::max (flatDev, std::abs (flat.g[(std::size_t) i] - flat.g[(std::size_t) (i + 1)]));
        ok (flatDev < 1.0e-4, "on the sinc even a period of 1 fits to " + std::to_string (flatDev)
                              + " — the modulation the closure argument was about is gone, so the "
                                "minimality assertion above would pass vacuously here");
    }

    // ------------------------------------------------------------------------------------------------
    group ("BLOCK SIZE — the cost is a property of the kernel, not of the caller's buffer");
    {
        // Compare the COMPLEX coherent gain, not its magnitude. A crew mutation that flipped polarity
        // survived a |mean| comparison, and the worst phase is min|g|, invariant under g -> -g just as
        // |mean| is, so it closed nothing. The complex mean carries the sign and the phase.
        const auto ref = roundTripGains<StreamResampler> (17640.0, H, M, 64, kSkip, kKeep);
        const std::complex<double> cref = coherent (ref);
        for (int b : { 1, 17, 63, 128, 512 })
        {
            const auto g = roundTripGains<StreamResampler> (17640.0, H, M, b, kSkip, kKeep);
            // 3e-3 was inherited from the cubic, where the block-to-block spread was real. On this
            // kernel the measured spread is 1.5e-10 … 1.6e-9, so 3e-3 admitted a 0.001-sample
            // block-dependent slip without a murmur. 1e-6 is still ~600x the measurement and rejects it.
            ok (std::abs (coherent (g) - cref) < 1.0e-6,
                "block " + std::to_string (b) + ": the COMPLEX carrier matches block 64 (|delta| = "
                + std::to_string (std::abs (coherent (g) - cref)) + ") — magnitude alone would not see a sign flip");
            approx (worstDb (g), worstDb (ref), 0.002, "block " + std::to_string (b) + ": …and so does the worst phase");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("COUNT — the EXACT number of outputs a feed yields, not a bound on it");
    {
        // The produced-count FORMULA is unchanged by the kernel swap: the priming moved from
        // (len = 3, pos = 1) to (len = kTaps, pos = kHalf) and the availability test from `i + 2 < len`
        // to `i + kHalf < len` — both by the same amount — so in exact arithmetic
        //   floor(kHalf + k·r) + kHalf < N + kTaps   <=>   floor(k·r) < N   <=>   K = ceil(N/r).
        //
        // 🔴 IN FLOAT IT IS NOT ALWAYS THE SAME NUMBER, and an earlier version of this comment claimed
        // it was. A crew round measured the exception: `pos` now accumulates from 32.0 instead of 1.0,
        // and at lengths that land EXACTLY on an integer boundary that changes which side of it the
        // sum falls on. Measured against the cubic: N = 147 at 44.1->48 gives 161 where the cubic gave
        // 160, and N = 160 at 48->44.1 gives 147 where it gave 148. Every other length in the sweep
        // agrees exactly. This is harmless inside NamStage — the down leg is bounded by
        // maxModelFrames and the up leg is asked for exactly n — but "unchanged" was the wrong word
        // and the difference is a real ULP fact, not a rounding of the prose.
        //
        // A mutation that stops the loop one sample early stays inside the theory suite's ±2 drift
        // bound at every length and survives; asserted EXACTLY against the class's own accumulator
        // (which has no slack at all), it does not.
        //
        // The expectation is not a closed form with slack — it is the SAME double accumulation the class
        // does, which has no slack at all: 48000/44100 = 160/147 is not a binary fraction and can land a
        // few ULP either side of an integral boundary, while 96000/48000 = 2 accumulates exactly.
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
                int want = 0;
                for (double p = (double) StreamResampler::kHalf;
                     std::floor (p) + StreamResampler::kHalf < (double) (N + StreamResampler::kTaps);
                     p += rb.in / rb.out) ++want;
                if (k != want)
                {
                    allExact = false;
                    std::printf ("      %6.0f -> %6.0f  N=%d produced %d, derived %d\n", rb.in, rb.out, N, k, want);
                }
            }
        ok (allExact, "one feed of N yields EXACTLY the count the class's own phase accumulator implies,"
                      " at every ratio and length — no slack anywhere");

        // SIGN, at the level of ONE stage. The round-trip sign check cannot see a kernel that inverts,
        // because both stages invert and the two cancel — which is precisely why a crew mutation that
        // negated every fed sample survived every suite. A single stage cannot cancel with itself: a
        // settled constant must come out as ITSELF, sign included, because every phase row is
        // normalised to a partition of unity.
        {
            StreamResampler r;
            r.reset (48000.0, 44100.0, 4096);
            std::vector<float> x (2048, 0.25f), o (4096);
            r.feed (x.data(), 2048);
            const int k = r.produceAvailable (o.data(), (int) o.size());
            double worstDev = 0.0;
            for (int i = 64; i < k - 64; ++i) worstDev = std::max (worstDev, std::fabs ((double) o[(std::size_t) i] - 0.25));
            ok (k > 1800 && worstDev < 1.0e-6,
                "one stage passes a settled +0.25 as +0.25 (worst deviation " + std::to_string (worstDev)
                + ") — partition of unity at every phase, sign included");
        }

        // produceExact()'s documented behaviour on a startup underflow is "pad with silence", and a
        // mutation that pads with the PREVIOUS sample instead once survived every suite. Silence and a
        // hold are audibly different at a stream start — a held sample is a DC step. NOTE: inside
        // NamStage this branch is structurally unreachable (K_up = ceil(K_down·h/m) >= n for every n),
        // and the longer kernel does NOT change that, because the availability test moved with the
        // priming. This pins the CLASS's contract, which other callers can reach.
        {
            StreamResampler r;
            r.reset (44100.0, 48000.0, 512);
            std::vector<float> out (64, 0.5f);
            r.produceExact (out.data(), 64);
            // `!(fabs(v) > 0)` would also be true for a NaN. `fabs(v) <= 0` is true ONLY for an exact
            // zero (fabs is never negative) and false for NaN, which is the predicate this test wants.
            bool allZero = true;
            for (float v : out) allZero = allZero && (std::fabs (v) <= 0.0f);
            ok (allZero, "produceExact on a freshly reset resampler writes SILENCE, not held samples");

            StreamResampler r2;
            r2.reset (44100.0, 48000.0, 512);
            std::vector<float> in (8, 0.75f), out2 (200, -1.0f);
            r2.feed (in.data(), 8);
            r2.produceExact (out2.data(), 200);          // 8 in -> ceil(8/r) = 9 out, the other 191 are padding
            int lastNonZero = -1;
            for (int i = 0; i < 200; ++i) if (! (std::fabs (out2[(std::size_t) i]) <= 0.0f)) lastNonZero = i;
            bool tailSilent = true;
            for (int i = lastNonZero + 1; i < 200; ++i)
                tailSilent = tailSilent && (std::fabs (out2[(std::size_t) i]) <= 0.0f);
            // NOTE, because it is the kind of thing that gets "fixed" into a wrong assertion: 8 fed
            // samples DO produce 9 outputs even though the 64-tap window is mostly leading zeros. The
            // count depends on the aperture's FAR edge (output k needs input up to k·r), not on the
            // window being full — which is the same inequality as the cubic's and why the produced
            // count did not move. Those 9 are the FIR ramping up, and everything past them is silence.
            ok (lastNonZero < 12 && tailSilent,
                "…and the tail past what the history can produce is silence too (last non-zero at "
                + std::to_string (lastNonZero) + " of 200)");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("DECIMATION — the 48 -> 44.1 stage now HAS a stopband, and this is it");
    {
        // Before: at phase t = 0 the cubic's weights were (0,1,0,0) — a bare sample pick — so a tone
        // above the output Nyquist came through at −3 dB rms and 0.0 dB sample PEAK. Both halves of
        // that are gone: the peak now tracks the rms, because no phase of this kernel is a sample pick.
        struct Row { double f, rms, oldRms; };
        const Row rows[] = {
            { 22100.0,  -9.083, -3.00 },
            { 22500.0, -15.414, -3.00 },
            { 23000.0, -27.454, -3.06 },
            { 23500.0, -47.563, -3.06 },
            { 23900.0, -88.772, -3.07 },
        };
        for (const auto& row : rows)
        {
            StreamResampler r;
            r.reset (M, H, 8192);
            const int block = 86;
            double peak = 0.0, sq = 0.0; long n = 0;
            std::vector<float> b ((std::size_t) block), o (8192);
            for (int off = 0; off < 200000; off += block)
            {
                for (int i = 0; i < block; ++i) b[(std::size_t) i] = (float) std::sin (2.0 * kPi * row.f * (off + i) / M);
                r.feed (b.data(), block);
                const int got = r.produceAvailable (o.data(), (int) o.size());
                if (off > 40000)
                    for (int k = 0; k < got; ++k)
                    { const double v = o[(std::size_t) k]; peak = std::max (peak, std::abs (v)); sq += v * v; ++n; }
            }
            const double peakDb = dbOf (peak);
            const double rmsDb  = dbOf (std::sqrt (sq / (double) n) * std::sqrt (2.0));
            std::printf ("      %6.0f Hz above the 22.05 kHz Nyquist -> rms %7.2f dB  peak %7.2f dB   (cubic: %.1f / 0.0)\n",
                         row.f, rmsDb, peakDb, row.oldRms);
            approx (rmsDb, row.rms, 0.05, std::to_string ((int) row.f) + " Hz: the designed stopband, not a sample pick");
            // 🔴 THIS BOUND USED TO BE `rmsDb + 3.05` AND WAS THEATRE. `rmsDb` is already scaled by
            // sqrt(2), i.e. it IS the amplitude, so for a sinusoidal residual the peak equals it and
            // the honest margin is ~0, not 3. Measured through the same instrument, the REFERENCE CUBIC
            // — whose whole defect is a 0 dB sample peak — passed the old line at 22.1 kHz (margin
            // +0.069 dB) and at 22.5 kHz (+0.020). A bound the thing it exists to catch walks through
            // is not a bound. The four rows below 23.9 kHz measure peak - rms = 0.002 dB, so 0.10 is
            // fifty times the measurement and still catches a sample pick by 3 dB.
            if (row.f < 23800.0)
                ok (peakDb < rmsDb + 0.10,
                    std::to_string ((int) row.f) + " Hz: the sample PEAK tracks the rms to 0.10 dB (peak "
                    + std::to_string (peakDb) + ", rms " + std::to_string (rmsDb) + "). The cubic read "
                    "0.0 dB here because |M(0)| = 1 at every frequency for a (0,1,0,0) weight set; no "
                    "phase of this kernel is a bare sample pick");
            else
                // Deep in the stopband the residual is no longer a sinusoid — it is the sum of a few
                // very small images, so its crest factor is genuinely higher (+2.85 dB measured) and a
                // 0.10 dB rule would fail a correct kernel. What matters at this frequency is only that
                // the peak is nowhere near 0 dB, which is the claim the cubic failed by 80 dB.
                ok (peakDb < -80.0,
                    std::to_string ((int) row.f) + " Hz: the sample peak is " + std::to_string (peakDb)
                    + " dB, not the 0.0 dB the cubic reproduced here; deep in the stopband the residual "
                    "is not a single sinusoid, so its crest is bounded rather than tracked");
        }

        // WHERE it lands. A tone at g in (22.05, 24) kHz used to return as 44100 − g AND as g − 3900,
        // only 1.7 dB apart, so the damage covered 18.15–22.05 kHz rather than just the top slice. Both
        // are attenuated now, and the SECOND one — the folded image at 48000 − g — is gone entirely.
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
            std::printf ("      23 kHz in -> 21100 Hz %.2f dB (cubic -5.33) and 19100 Hz %.2f dB (cubic -7.04)\n", lo, hi);
            approx (lo, -27.45, 0.15, "23 kHz still folds to 21.1 kHz, but 22 dB down on the cubic");
            approx (hi, -91.58, 0.80, "…and the SECOND component, which used to sit 1.7 dB below the first, "
                                      "is 84 dB below where it was — that is what removed the 18.15-22.05 kHz damage");
        }
    }

    // ------------------------------------------------------------------------------------------------
    group ("THE CRITERION — the rate-match against a DRIVEN nonlinearity's OWN aliasing floor");
    {
        // The header used to justify the cubic with "the driven nonlinear stage masks the interpolation
        // images". Made measurable: the rate-match must not add MORE than the nonlinearity folds down by
        // itself. The cubic needed the nonlinearity driven to tanh(6.2x) — a near square wave — before
        // that was true. The same measurement on the new kernel is the acceptance.
        const double f0 = 17500.0;                    // on the 63-sample grid at 44.1 k and 96 at 48 k
        const int N = 441 * 200, skip = 20000, block = 64;

        const double oldAdded = addedDbc<CatmullRef>      (H, M, f0, block, skip, N);
        const double added    = addedDbc<StreamResampler> (H, M, f0, block, skip, N);
        std::printf ("      rate-match alone at 17.5 kHz: cubic %+.2f dBc -> sinc %+.2f dBc  (%.1f dB better)\n",
                     oldAdded, added, oldAdded - added);
        approx (oldAdded, -8.84, 0.05, "the cubic round trip added -8.84 dBc of artifacts at 17.5 kHz");
        ok (added < -100.0, "the sinc round trip adds " + std::to_string (added)
                            + " dBc — below anything a nonlinearity would have to mask, which RETIRES the "
                              "masking argument rather than winning it");

        // …and how hard the nonlinearity had to be driven before its OWN folding reached the cubic's
        // level. Kept because it calibrates the criterion and shows what the new number is measured
        // against: even a barely-driven tanh folds far more than the sinc round trip adds.
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
            std::printf ("      tanh drive %5.2f folds %.2f dBc by itself -> sinc sits %+.2f dB against it\n",
                         d.drive, fl, added - fl);
            approx (fl, d.floorDbc, 0.05, "tanh at drive " + std::to_string (d.drive) + " folds its own "
                    + std::to_string (d.floorDbc) + " dBc");
            ok (added < fl - 40.0,
                "…and the rate-match is at least 40 dB quieter than that floor at drive "
                + std::to_string (d.drive) + " — the cubic was LOUDER than it at every drive below 6.2");
        }
    }

    return felitronics::test::report();
}
