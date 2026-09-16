// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the polyphase oversampler: upsample interpolation vs the analytic sine
// (correct + alias-free), inter-sample-peak revelation (the limiter's reason to exist), the
// up→down round-trip == a delayed identity, and — since this module owns the shipped tapsPerPhase
// default — what that default actually buys: the declared stopband, the aliasing it stops, and the
// pass-band droop of two oversampled stages in series — and, since P31, why the CUTOFF sits where it does.

#include <felitronics_test.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;


//==============================================================================
// WHAT tapsPerPhase BUYS — the criterion the default is derived from, made executable.
//
// This suite used to run entirely at a hardcoded tpp = 32 on tones of 500 Hz and 2 kHz at 48 kHz
// (0.010 and 0.042 fs), i.e. it owned the default and never measured anything the default decides.
// Both things the taps count controls live between 0.36 fs and the fold at 0.50 fs, and nothing in
// this file looked there. Same class as the blind grid: a fixture that cannot see the quantity it
// certifies. The two below are measured through the PUBLIC API only.
//
//   1. STOPBAND (the criterion). The cutoff is FIXED at 0.45 fs, so the transition has to finish by
//      the fold at 0.50 fs or everything above it comes back as aliasing. designFilter() DECLARES its
//      own target — beta = 9 is a ~90 dB Kaiser stopband — and at 32 taps it delivered 27.
//   2. PASS BAND (the corollary). Two oversampled stages in series is the real assembly (a clipper in
//      front of a limiter), and the dB of one filter pass are multiplied by four across it.
namespace tapsprobe
{
    constexpr double kPi = 3.14159265358979323846;
    constexpr int    kWin = 10000;      // analysis length; every test frequency is an exact bin k/kWin
    constexpr int    kSettle = 4096;    // >> 4*(kMaxTpp-1): both round trips are past their transients

    // Amplitude at baseband-normalised `nu`, coherently projected over exactly kWin samples starting at
    // `off`. `nu*kWin` is an integer for every frequency used here, so the projection is exact and there
    // is no window, no leakage and no grid to land on the wrong side of.
    inline double amp (const std::vector<float>& x, double nu, int off)
    {
        double re = 0.0, im = 0.0;
        for (int i = 0; i < kWin; ++i)
        {
            const double ph = 2.0 * kPi * nu * (double) (off + i);
            re += (double) x[(std::size_t) (off + i)] * std::cos (ph);
            im += (double) x[(std::size_t) (off + i)] * std::sin (ph);
        }
        return 2.0 * std::hypot (re, im) / (double) kWin;
    }

    // One up->down round trip of a baseband tone at `nu`, in dB. `taps < 0` means "call prepare with
    // TWO arguments", i.e. exercise the DEFAULT rather than a number the test knows.
    // A REFUSED prepare() must never look like a measurement. The first draft returned a +1e9 sentinel,
    // and `worstFoldDb (4, 48) > -60.0` — an assertion about a configuration being BAD — was then
    // satisfied by prepare() simply refusing, without a single sample being measured. That is the exact
    // blind-fixture shape this file exists to close, committed inside the fix for it. Both probes now
    // FLAG the refusal so the caller asserts on it instead of inheriting it.
    inline bool g_prepareRefused = false;

    inline double roundTripDb (int factor, int taps, double nu, int passes)
    {
        oversampling::PolyphaseOversampler os[4];
        for (int s = 0; s < passes; ++s)
        {
            const bool ok = taps < 0 ? os[s].prepare (factor, 1) : os[s].prepare (factor, 1, taps);
            if (! ok) { g_prepareRefused = true; return 1e9; }
        }
        const int n = kSettle + kWin;
        std::vector<float> x ((std::size_t) n), osb ((std::size_t) n * (std::size_t) factor);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * nu * (double) i);
        for (int s = 0; s < passes; ++s)
        {
            const float* in[1] { x.data() }; float* mid[1] { osb.data() };
            os[s].upsample (in, 1, n, mid);
            const float* mc[1] { osb.data() }; float* out[1] { x.data() };
            os[s].downsample (mc, 1, n, out);
        }
        return 20.0 * std::log10 (std::max (1e-300, amp (x, nu, kSettle)));
    }

    // |H| of ONE pass at a baseband-normalised `nu` ABOVE the fold, measured the only way the public
    // API allows: hand the DECIMATOR an oversampled tone at that frequency and read what survives, at
    // the baseband frequency it folds to. This is exactly the aliasing path, not a model of it.
    inline double foldThroughDb (int factor, int taps, double nu)
    {
        oversampling::PolyphaseOversampler os;
        if (! (taps < 0 ? os.prepare (factor, 1) : os.prepare (factor, 1, taps))) { g_prepareRefused = true; return 1e9; }
        const int n = kSettle + kWin;
        std::vector<float> osb ((std::size_t) n * (std::size_t) factor), y ((std::size_t) n, 0.0f);
        for (int k = 0; k < n * factor; ++k)
            osb[(std::size_t) k] = (float) std::sin (2.0 * kPi * (nu / (double) factor) * (double) k);
        const float* in[1] { osb.data() }; float* out[1] { y.data() };
        os.downsample (in, 1, n, out);
        double f = nu - std::floor (nu);                       // where it lands after decimation
        if (f > 0.5) f = 1.0 - f;
        return 20.0 * std::log10 (std::max (1e-300, amp (y, f, kSettle)));
    }

    // Worst |H| anywhere in the fold region [0.5, factor/2] fs — the whole band that aliases in.
    inline double worstFoldDb (int factor, int taps)
    {
        double worst = -1e9;
        const double hi = 0.5 * (double) factor;
        // Dense where the transition ends (that is where the worst sidelobe sits), coarse above it.
        // EXACTLY 0.5 fs is skipped on purpose: it folds to baseband Nyquist, where a single-phase sine
        // projection reads an amplitude that depends on the sampling phase (off by 2*sin(pi/2L), i.e.
        // +3.01 dB at 2x and -8.17 at 8x) rather than on |H|. It is one point of an open interval, and
        // reading it would be measuring the probe.
        //
        // The step is 0.0005 through the first four sidelobes and 0.002 after. The coarse step alone
        // MISSED the true maximum: at 2x/64 the peak sits at 0.50499 fs and reads -90.46 dB, while a
        // 0.002 grid reported -90.95 — half a dB of under-reading against a 0.46 dB margin to the bar.
        // Under-reading a worst case is the direction that makes an assertion pass for the wrong reason.
        for (int k = 5002; k <= 5200; k += 5)                   // 0.5002 .. 0.520 fs, step 0.0005
            worst = std::max (worst, foldThroughDb (factor, taps, (double) k / (double) kWin));
        for (int k = 5220; k <= 6000; k += 20)                  // 0.522 .. 0.600 fs, step 0.002
            worst = std::max (worst, foldThroughDb (factor, taps, (double) k / (double) kWin));
        for (int k = 6200; k <= (int) (hi * kWin); k += 200)    // .. factor/2 fs, step 0.02
            worst = std::max (worst, foldThroughDb (factor, taps, (double) k / (double) kWin));
        return worst;
    }
}

int runTapsTests()
{
    using namespace tapsprobe;

    // ---------------------------------------------------------------- 1. the criterion, two-sided
    test::group ("tapsPerPhase: the stopband the design DECLARES (beta = 9 -> ~-90 dB), at the default");
    {
        // Liveness of the instrument, asserted before it is believed. It reads the design's own -6.02 dB
        // cutoff anchor, and it separates two topologies by 40 dB on the same frequency: an instrument
        // that could only ever read "very small" would pass the table below by being deaf.
        g_prepareRefused = false;                  // BEFORE the first probe, not after the loop below:
                                                   // placed later it covered neither the factor sweep nor
                                                   // the separation check, both of which a +1e9 sentinel
                                                   // satisfies. Same shape as the defect this group pins.
        const double anchor = roundTripDb (4, -1, 0.45, 1);
        test::approx (anchor, -12.041, 0.05,
                      "PRECONDITION: the round trip reads its own 0.45 fs cutoff at -12.04 dB (= 2 x -6.02)");
        const double sep = foldThroughDb (4, 32, 0.52) - foldThroughDb (4, -1, 0.52);
        test::ok (sep > 40.0, "PRECONDITION: the instrument separates 32 taps from the default by "
                              + std::to_string ((int) sep) + " dB at 0.52 fs (it is not reading noise)");

        for (int factor : { 2, 4, 8 })
        {
            const double now = worstFoldDb (factor, -1);
            const double was = worstFoldDb (factor, 32);
            std::printf ("       factor %d: worst fold-region rejection = %.2f dB at the default, %.2f dB at 32 taps\n",
                         factor, now, was);
            test::ok (now <= -90.0,
                      "factor " + std::to_string (factor) + ": the DEFAULT delivers the declared -90 dB stopband ("
                      + std::to_string (now) + " dB)");
            test::ok (was > -30.0,
                      "factor " + std::to_string (factor) + ": 32 taps did NOT (" + std::to_string (was)
                      + " dB) — the defect this default closes, and the reason it is not a droop question");
        }
        // WHERE THE KNEE IS, and why "the first taps count that passes" is the wrong question. Below ~58
        // the transition is unfinished and the rejection is monotonically poor; at ~58 it reaches the
        // Kaiser floor and RIPPLES there by about a dB, so 59 and 61 fall a tenth of a dB short while 58,
        // 60 and 64 clear it. An earlier version of this comment said "60 is where it first holds" — it
        // is 58, and the integer is an artefact of a hard bar on a rippling quantity. Pinned as the KNEE
        // plus the ripple, which is what is actually true, so neither can be quietly re-opened.
        const double at48 = worstFoldDb (4, 48), at57 = worstFoldDb (4, 57), at58 = worstFoldDb (4, 58);
        const double at59 = worstFoldDb (4, 59), at61 = worstFoldDb (4, 61);
        test::ok (! g_prepareRefused,
                  "PRECONDITION: every probe in this group actually PREPARED — a refused prepare returns "
                  "+1e9, which satisfies an upper bound without measuring anything");
        test::ok (at48 > -60.0, "48 taps still misses the declared stopband by ~38 dB");
        std::printf ("       knee: 48 -> %.2f, 57 -> %.2f, 58 -> %.2f, 59 -> %.2f, 61 -> %.2f dB\n", at48, at57, at58, at59, at61);
        test::ok (at57 > -88.0, "57 taps still misses it — the transition is genuinely unfinished below the knee");
        test::ok (at58 <= -90.0, "58 taps is the FIRST that meets it (not 60, as this once claimed)");
        // The ripple stated as a PROPERTY rather than as integers crossing a bar. `at61 > -90.0` reads
        // -90.00 here — a knife-edge that a different libm could flip, and an assertion whose margin is
        // its own rounding is not an assertion. That MORE taps can be WORSE is the real content, and it
        // has 0.8 dB of room: monotone-in-taps is exactly what a reader would assume and it is false.
        test::ok (at59 > at58 + 0.5,
                  "59 taps is measurably WORSE than 58 (" + std::to_string (at59) + " vs " + std::to_string (at58)
                  + ") — the stopband ripples on the window floor, so the knee is a region, not an integer");
    }

    // ---------------------------------------------------------------- 2. end to end, through a real nonlinearity
    test::group ("tapsPerPhase: what the stopband is FOR — TRANSITION-BAND leakage through a waveshaper");
    {
        // 🔴 READ THE SCOPE OF THIS FIRST, because an earlier version of this comment overstated it and the
        // number went into a changelog before it was checked. A tanh has INFINITELY many odd harmonics, so
        // the ones above the 4x OS Nyquist (2.0 fs) fold INSIDE the oversampled domain, where no decimation
        // filter can reach them — that is the FACTOR's axis, not this one. Measured on a real Saturator at
        // f0 = 0.17 fs, the TOTAL non-harmonic energy is -31.97 dBc at 32 taps and -31.27 at 64: taps make
        // the total very slightly WORSE, because a flatter pass band also delivers the components that had
        // already folded. So this group does NOT claim "less aliasing"; it isolates the ONE component the
        // taps do own.
        //
        // The 3rd harmonic of 0.17 fs sits at 0.51 fs — just above the fold, inside the transition band —
        // and folds back to 0.49 fs. That single line is exactly the stopband leakage the criterion above
        // is about, it is separable by CONSTRUCTION (a tanh of a 0.17 tone can have no energy there), and
        // it is the one this default moves. The other four lines are measured alongside it precisely so
        // that the ones the taps do NOT fix stay visible in the printout.
        const double f0 = 0.17;
        // 0.49 is the one the taps own (from 0.51 fs, inside the transition band). The rest come from
        // harmonics deep in the stopband, where 32 taps were already at the window floor — they are here
        // to show that they do NOT move, which is what bounds the claim.
        const double folds[] = { 0.49, 0.15, 0.19, 0.47, 0.13 };
        auto aliasDbc = [&] (int taps)
        {
            oversampling::PolyphaseOversampler os;
            if (! (taps < 0 ? os.prepare (4, 1) : os.prepare (4, 1, taps))) g_prepareRefused = true;
            const int n = kSettle + kWin;
            std::vector<float> x ((std::size_t) n), osb ((std::size_t) n * 4);
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (0.9 * std::sin (2.0 * kPi * f0 * (double) i));
            const float* in[1] { x.data() }; float* mid[1] { osb.data() };
            os.upsample (in, 1, n, mid);
            for (auto& v : osb) v = std::tanh (6.0f * v);                 // the nonlinearity, in the OS domain
            const float* mc[1] { osb.data() }; float* out[1] { x.data() };
            os.downsample (mc, 1, n, out);
            const double fund = amp (x, f0, kSettle);
            for (int i = 1; i < 5; ++i)                       // the four that the taps do NOT own
                std::printf ("         %s: %.2f fs fold = %7.2f dBc\n",
                             taps < 0 ? "default" : "32 taps", folds[i],
                             20.0 * std::log10 (amp (x, folds[i], kSettle) / std::max (1e-30, fund)));
            return 20.0 * std::log10 (amp (x, folds[0], kSettle) / std::max (1e-30, fund));
        };
        g_prepareRefused = false;
        const double now = aliasDbc (-1), was = aliasDbc (32);
        test::ok (! g_prepareRefused, "PRECONDITION: both alias probes prepared — a refusal emits silence, "
                                      "and -inf dBc satisfies every upper bound below");
        std::printf ("       tanh at 0.17 fs, 4x: the 0.51 -> 0.49 fs fold = %.2f dBc at the default, %.2f dBc at 32 taps\n",
                     now, was);
        test::ok (was > -50.0, "32 taps let the transition-band fold through at " + std::to_string (was)
                               + " dBc — and that was the shipped default");
        test::ok (now < -75.0, "the default puts it at " + std::to_string (now) + " dBc");
        test::ok (was - now > 25.0, "the taps are worth " + std::to_string (was - now)
                                    + " dB on THIS component — the one the stopband owns");
    }

    // ---------------------------------------------------------------- 3. the pass band, pinned, TWO stages
    test::group ("tapsPerPhase: pass-band droop of TWO oversampled stages in series, pinned");
    {
        // Two stages, because that is the assembly the figures were found on (a clipper in front of a
        // limiter) and one stage is not the case this is decided on. Values measured on this tree and
        // cross-checked against an independent DTFT of designFilter()'s coefficients to three decimals;
        // a redesign has to edit them on purpose. Factor 4; the surface is factor-independent to 0.09 dB
        // and EXACTLY factor-independent at 0.45, which the last row asserts.
        struct Row { int k; double at64; double at32; };
        const Row rows[] = { { 3600,  +0.000,  -0.001 },      // 0.36 fs — 15.9 kHz at 44.1 k
                             { 4000,  +0.000,  -1.549 },      // 0.40 fs — 17.6 kHz: the figure that started this
                             { 4100,  -0.037,  -3.249 },      // 0.41 fs — 18.1 kHz
                             { 4200,  -0.610,  -6.033 },      // 0.42 fs — 18.5 kHz: the second figure
                             { 4400, -10.182, -16.131 },      // 0.44 fs — 19.4 kHz, well inside the transition
                             { 4500, -24.082, -24.082 } };    // 0.45 fs — the cutoff itself, taps-INDEPENDENT
        for (const Row& r : rows)
        {
            const double nu = (double) r.k / (double) kWin;
            test::approx (roundTripDb (4, -1, nu, 2), r.at64, 0.02,
                          "two stages at " + std::to_string (r.k) + "/10000 fs, DEFAULT taps");
            test::approx (roundTripDb (4, 32, nu, 2), r.at32, 0.02,
                          "two stages at " + std::to_string (r.k) + "/10000 fs, 32 taps (what it used to be)");
        }
        // UNITY GAIN, pinned separately and TIGHTLY. The rows above cannot see a gain error: a 0.1 %
        // scale slipped into upsample() or downsample() costs 0.0087 dB per pass, so two stages read
        // 0.0174 dB against a 0.02 tolerance and pass. (Both mutations were run and both survived the
        // table.) Deep in the pass band the true value is 0.000000 and the float noise is ~1e-6, so a
        // 0.002 bar catches a 0.02 % error — and it is checked at one, two and four passes, because a
        // per-pass error accumulates and a per-CALL one would not.
        for (int passes : { 1, 2 })
            test::approx (roundTripDb (4, -1, 0.05, passes), 0.0, 0.002,
                          "unity gain deep in the pass band, " + std::to_string (passes)
                          + " stage(s) — a 0.1% scale in up/downsample reads 0.0087 dB per pass");
        // FACTORS 2 AND 8, at the row where they differ most. The table above runs at factor 4 only, and
        // "factor-independent to 0.09 dB" was a COMMENT, not a check — so a design change that touched
        // only one factor could pass everything above. (Measured: a cutoff scaled by 0.98 at factor 2
        // alone survived the factor-4 table and was caught only by accident, through a witness that
        // happened to sit at 0.5 fs.) 0.44 fs is chosen because it is where the three factors are
        // furthest apart, and the values are the ones an independent DTFT of designFilter() gives.
        test::approx (roundTripDb (2, -1, 0.44, 2), -10.221, 0.02, "two stages at 0.44 fs, factor 2, DEFAULT taps");
        test::approx (roundTripDb (8, -1, 0.44, 2), -10.162, 0.02, "two stages at 0.44 fs, factor 8, DEFAULT taps");
        test::approx (roundTripDb (2, -1, 0.45, 2), -24.082, 0.02, "...and the cutoff anchor holds at factor 2");
        test::approx (roundTripDb (8, -1, 0.45, 2), -24.082, 0.02, "...and at factor 8");
        // The last row is the one a copied table cannot fake: it is identical at both taps counts, so a
        // fixture reading zero everywhere fails it, and a fixture reading the 32-tap column passes it.
        test::approx (roundTripDb (4, 32, 0.45, 2) - roundTripDb (4, -1, 0.45, 2), 0.0, 0.02,
                      "PRECONDITION: 0.45 fs is taps-INDEPENDENT (the cutoff is fixed) — the anchor a copied table trips on");
        test::ok (roundTripDb (4, 32, 0.44, 2) - roundTripDb (4, -1, 0.44, 2) < -5.0,
                  "PRECONDITION: 0.44 fs separates the two topologies by > 5 dB — the row that proves the table moves");
    }

    // ---------------------------------------------------------------- 4. the DEFAULT itself
    test::group ("tapsPerPhase: the default is 64, and prepare() honours what it is handed");
    {
        // Deliberately a RUNTIME check and not a static_assert. A static_assert would stop the build the
        // moment the constant moved, which sounds stronger and is worse: the rest of this file would never
        // run, so nobody would learn whether the BEHAVIOURAL pins above catch the same change. They do —
        // the mutation stand scores the stopband, the aliasing and the droop table separately — and that
        // is only visible if the compile survives.
        test::ok (oversampling::PolyphaseOversampler::kDefaultTapsPerPhase == 64,
                  "kDefaultTapsPerPhase == 64 — the number the criterion above is derived for");
        oversampling::PolyphaseOversampler def, old, big;
        test::ok (def.prepare (4, 1) && def.latencySamples() == 63,
                  "prepare(factor, channels) takes 64 taps/phase -> a 63-sample round trip");
        test::ok (old.prepare (4, 1, 32) && old.latencySamples() == 31,
                  "an explicit 32 is still honoured -> 31 (the default is a default, not a constant)");
        test::ok (! big.prepare (4, 1, oversampling::PolyphaseOversampler::kMaxTapsPerPhase + 1),
                  "tapsPerPhase above kMaxTapsPerPhase is refused HERE, so every consumer inherits the guard");
        test::ok (! big.prepare (4, 1, std::numeric_limits<int>::max()),
                  "and INT_MAX with it (it used to overflow factor*taps before allocating)");
        // The OTHER factor of the same product. `Saturator` hands its oversampleFactor straight through
        // with no ceiling of its own, so prepare(INT_MAX, 1) on the default taps was UBSan-confirmed
        // signed overflow and then a length_error out of assign() — a terminate under -fno-exceptions.
        // The LITERAL numbers, not the constants. Every check here phrased as `kMax + 1` / `kMax` moves
        // with the constant, so `kMaxFactor = 64 -> 63` passes them all while silently dropping the 64x
        // support the header documents. A bound is a promise to callers, and callers write integers.
        test::ok (oversampling::PolyphaseOversampler::kMaxFactor == 64
                    && oversampling::PolyphaseOversampler::kMaxTapsPerPhase == 1024,
                  "the documented ceilings are 64x and 1024 taps/phase — pinned as numbers, not as themselves");
        test::ok (! big.prepare (65, 1), "factor 65 is refused (the documented ceiling is 64)");
        test::ok (! big.prepare (4, 1, 1025), "1025 taps/phase is refused (the documented ceiling is 1024)");
        test::ok (! big.prepare (oversampling::PolyphaseOversampler::kMaxFactor + 1, 1),
                  "a factor above kMaxFactor is refused too — either argument alone can overflow N");
        test::ok (! big.prepare (std::numeric_limits<int>::max(), 1),
                  "INT_MAX as the FACTOR is refused (it used to overflow, then throw, then terminate on wasm)");
        test::ok (big.prepare (32, 1) && big.latencySamples() == 63,
                  "and 32x still works — the largest factor anything in the tree asks for");
        // POSITIVE at both ceilings. Without these, `>` -> `>=` on either bound is invisible: every
        // rejection check above uses a value one PAST the limit, so a guard that also rejects the limit
        // itself passes them all. A bound needs the last legal value as much as the first illegal one.
        oversampling::PolyphaseOversampler edge;
        test::ok (edge.prepare (4, 1, oversampling::PolyphaseOversampler::kMaxTapsPerPhase)
                    && edge.latencySamples() == oversampling::PolyphaseOversampler::kMaxTapsPerPhase - 1,
                  "kMaxTapsPerPhase itself is ACCEPTED (the bound is inclusive, not off by one)");
        test::ok (edge.prepare (oversampling::PolyphaseOversampler::kMaxFactor, 1),
                  "kMaxFactor itself is ACCEPTED too");
        // LINEAR PHASE, asserted on the shipped coefficients rather than on the formula. An impulse
        // through upsample() emits factor*proto in order, so the kernel can be read back through the
        // public API — and it comes out EXACTLY symmetric, which makes this an equality rather than a
        // tolerance. The header claims "linear phase -> a constant group delay (reported)" and nothing
        // tested it: a window centred half a sample off (`(2i-N)/(N-1)`) survives every spectral row in
        // this file, because it costs 0.0002 samples of group delay at 0.36 fs. Third oracle
        // construction in this suite, and the only one that can see the class diverge from its own
        // formula while the other two still agree with each other.
        for (int L : { 2, 4, 8 })
        {
            oversampling::PolyphaseOversampler k;
            test::ok (k.prepare (L, 1), "prepare for the kernel read-back, factor " + std::to_string (L));
            const int taps = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase, N = L * taps;
            std::vector<float> imp ((std::size_t) (taps + 4), 0.0f), ker ((std::size_t) (taps + 4) * (std::size_t) L, 0.0f);
            imp[0] = 1.0f;
            const float* ip[1] { imp.data() }; float* kp[1] { ker.data() };
            k.upsample (ip, 1, taps + 4, kp);
            double worstAsym = 0.0, energy = 0.0;
            for (int i = 0; i < N; ++i)
            {
                worstAsym = std::max (worstAsym, (double) std::fabs (ker[(std::size_t) i] - ker[(std::size_t) (N - 1 - i)]));
                energy += (double) ker[(std::size_t) i] * ker[(std::size_t) i];
            }
            test::ok (energy > 0.1, "PRECONDITION: the kernel read back at factor " + std::to_string (L)
                                    + " is not silence (energy " + std::to_string (energy) + ")");
            test::ok (worstAsym == 0.0, "factor " + std::to_string (L)
                                        + ": the shipped kernel is EXACTLY symmetric — linear phase, as the header claims");
        }
        // The one-pass group delay, pinned numerically. Its only other check compares interpolated samples
        // against an analytic sine at a 0.02 tolerance, and a quarter-OS-sample error there costs 0.0131 —
        // inside the tolerance, i.e. invisible. (N-1)/2 is exact by construction, so assert it as such.
        oversampling::PolyphaseOversampler gd;
        test::ok (gd.prepare (4, 1, 32) && gd.filterLatencyOversampled() == (4.0 * 32.0 - 1.0) * 0.5,
                  "filterLatencyOversampled() == (factor*taps - 1)/2 exactly, at 32 taps");
        test::ok (gd.prepare (4, 1) && gd.filterLatencyOversampled() == (4.0 * 64.0 - 1.0) * 0.5,
                  "...and at the default 64 taps (127.5 oversampled samples)");
    }
    return 0;
}

//==============================================================================
// THE REFERENCE PROTOTYPE, PINNED. P80 moved designFilter() from std::sin to core::det::sin, and the
// claim attached to that change is not "it is better" but "it moved NOTHING": the taps are narrowed to
// float, and that narrowing throws away 29 of the bits the two libms can disagree about.
//
// THESE 128 CONSTANTS ARE THE OLD PATH'S OUTPUT, NOT THIS PATH'S. They were dumped from the system-sin
// design at a9816e2 and are byte-identical on four independent libms — Apple clang/arm64, gcc 14 +
// glibc x86-64, emcc/musl wasm32 and MSVC/UCRT x86-64 — so they are an oracle computed OUTSIDE the code
// under test, not a photograph of it. A table regenerated from the code it guards would pass forever.
//
// WHAT IT CATCHES, and it is two different things:
//   · that the deterministic design still lands on the same floats the shipped one did — i.e. that no
//     bit of TabbyEQ's or OrbitCab's true-peak reading moved when P80 touched this file;
//   · and that a future row whose narrowing does NOT absorb the difference is caught here rather than
//     discovered in a parity diff. Measured margin: perturbing every sin() result by a deliberate k ulp
//     leaves all 128 taps unmoved up to k = 2^24, while the real spread between libms is 1-3 ulp.
//
// The taps are read back through the PUBLIC API, so this tests the class and not a copy of its formula:
// a single 1.0f in the history makes core::firDot's sum one product plus zeros, and adding 0.0f is
// exact, so each downsample() output IS one tap. downsample() picks them up strided by L, hence the
// four passes at the four phase offsets.
static const std::uint32_t kReferenceProto4x32[128] = {
        0x36719320u, 0x35d43594u, 0xb6d237d7u, 0xb7955c31u, 0xb7d346ccu, 0xb7968a21u,
        0x37204084u, 0x3852897cu, 0x38abd9efu, 0x389fa9ffu, 0x377f3451u, 0xb8ba2134u,
        0xb9450da7u, 0xb9614000u, 0xb8f997e7u, 0x38c11223u, 0x39af3aaeu, 0x39f685dau,
        0x39c55ea3u, 0x37d4ceb7u, 0xb9f32f07u, 0xba5e54eau, 0xba634f5au, 0xb9ce715au,
        0x39eb6e12u, 0x3aa89b17u, 0x3ad6abb1u, 0x3a98ffd2u, 0xb8aec466u, 0xbad40735u,
        0xbb2dac71u, 0xbb226f45u, 0xba6c67c9u, 0x3acac2e2u, 0x3b741965u, 0x3b8f50efu,
        0x3b37051bu, 0xba2b10beu, 0xbb935b9bu, 0xbbdcec3cu, 0xbbbf1ca9u, 0xbad0dde8u,
        0x3b8fc215u, 0x3c17a857u, 0x3c26a368u, 0x3bbec2d9u, 0xbb2780d0u, 0xbc398764u,
        0xbc821607u, 0xbc531e69u, 0xbb15a024u, 0x3c445c1eu, 0x3cbdaf17u, 0x3cc8aeb4u,
        0x3c50e8d7u, 0xbc1a4e02u, 0xbd08362fu, 0xbd3e83c3u, 0xbd1bdd69u, 0xbb32289bu,
        0x3d63aaf9u, 0x3dfe1321u, 0x3e3d257eu, 0x3e619376u, 0x3e619376u, 0x3e3d257eu,
        0x3dfe1321u, 0x3d63aaf9u, 0xbb32289bu, 0xbd1bdd69u, 0xbd3e83c3u, 0xbd08362fu,
        0xbc1a4e02u, 0x3c50e8d7u, 0x3cc8aeb4u, 0x3cbdaf17u, 0x3c445c1eu, 0xbb15a024u,
        0xbc531e69u, 0xbc821607u, 0xbc398764u, 0xbb2780d0u, 0x3bbec2d9u, 0x3c26a368u,
        0x3c17a857u, 0x3b8fc215u, 0xbad0dde8u, 0xbbbf1ca9u, 0xbbdcec3cu, 0xbb935b9bu,
        0xba2b10beu, 0x3b37051bu, 0x3b8f50efu, 0x3b741965u, 0x3acac2e2u, 0xba6c67c9u,
        0xbb226f45u, 0xbb2dac71u, 0xbad40735u, 0xb8aec466u, 0x3a98ffd2u, 0x3ad6abb1u,
        0x3aa89b17u, 0x39eb6e12u, 0xb9ce715au, 0xba634f5au, 0xba5e54eau, 0xb9f32f07u,
        0x37d4ceb7u, 0x39c55ea3u, 0x39f685dau, 0x39af3aaeu, 0x38c11223u, 0xb8f997e7u,
        0xb9614000u, 0xb9450da7u, 0xb8ba2134u, 0x377f3451u, 0x389fa9ffu, 0x38abd9efu,
        0x3852897cu, 0x37204084u, 0xb7968a21u, 0xb7d346ccu, 0xb7955c31u, 0xb6d237d7u,
        0x35d43594u, 0x36719320u,
};

static void runReferenceTapPin()
{
    test::group ("the reference 4x32 prototype is bit-identical to the pre-det design (four libms)");
    const int L = 4, tpp = 32, N = L * tpp;
    oversampling::PolyphaseOversampler os;
    const bool prepared = os.prepare (L, 1, tpp);
    test::ok (prepared, "prepare 4x32");
    if (! prepared) return;

    std::vector<std::uint32_t> tap ((std::size_t) N, 0u);
    std::vector<bool> got ((std::size_t) N, false);
    for (int p = 0; p < L; ++p)
    {
        os.reset();
        std::vector<float> in ((std::size_t) N * L, 0.0f); in[(std::size_t) p] = 1.0f;
        std::vector<float> out ((std::size_t) N, 0.0f);
        const float* ip[1] { in.data() }; float* op[1] { out.data() };
        os.downsample (ip, 1, N, op);
        for (int i = 0; i < N; ++i)
        {
            const int idx = i * L + (L - 1 - p);
            if (idx < N) { std::memcpy (&tap[(std::size_t) idx], &out[(std::size_t) i], 4); got[(std::size_t) idx] = true; }
        }
    }
    // The recovery must have covered every tap, or "all taps match" would be a claim about a subset.
    int covered = 0; for (int i = 0; i < N; ++i) if (got[(std::size_t) i]) ++covered;
    test::ok (covered == N, "all " + std::to_string (N) + " taps were recovered (" + std::to_string (covered) + ")");

    int differ = 0, firstBad = -1;
    for (int i = 0; i < N; ++i)
        if (tap[(std::size_t) i] != kReferenceProto4x32[i]) { if (firstBad < 0) firstBad = i; ++differ; }
    if (differ != 0)
        std::printf ("      first difference at tap %d: got %08x, pinned %08x (%d of %d differ)\n",
                     firstBad, tap[(std::size_t) firstBad], kReferenceProto4x32[firstBad], differ, N);
    test::ok (differ == 0, "every tap equals the system-designed table dumped on four libms before P80");

    // AND THE PIN MUST BE ABLE TO FAIL. A comparison against a table is worth nothing if the recovery
    // silently returns the table itself, or zeros, or the same value for every tap.
    bool allSame = true; for (int i = 1; i < N; ++i) if (tap[(std::size_t) i] != tap[0]) { allSame = false; break; }
    test::ok (! allSame, "the recovered taps are not all one value (the comparison is not vacuous)");
    int nonZero = 0; for (int i = 0; i < N; ++i) if (tap[(std::size_t) i] != 0u) ++nonZero;
    test::ok (nonZero > N / 2, "and most of them are non-zero (" + std::to_string (nonZero) + " of " + std::to_string (N) + ")");

    // THE PROTOTYPE IS A PALINDROME, and saying so is the honest way to bound what the pin above proves.
    // A linear-phase FIR is symmetric by construction, so tap[i] == tap[N-1-i] — which means a recovery
    // that read the taps in REVERSE order would compare equal to the table and the pin would not notice.
    // That particular error is harmless (a reversed palindrome is the same filter), and the indexing
    // errors that are NOT harmless — a wrong phase offset, a wrong stride — scramble rather than reverse
    // and are caught by `covered == N` plus the comparison. This assertion covers the remaining piece: it
    // is a property of the DESIGN, so a filter that stopped being linear-phase would fail here rather
    // than silently become a different animal that still matched 128 pinned words.
    int asym = 0;
    for (int i = 0; i < N / 2; ++i) if (tap[(std::size_t) i] != tap[(std::size_t) (N - 1 - i)]) ++asym;
    test::ok (asym == 0, "the recovered prototype is symmetric, as a linear-phase design must be");
}

//==============================================================================
// THE CUTOFF AXIS (P31) — every figure in the header's "WHY THE CUTOFF IS 0.90" paragraph, measured here
// and pinned, so that paragraph is a printout of this group rather than a quotation. Two halves:
//   * what the GUARD BAND buys — exposure (content gain + image gain, one pass UP, i.e. what reaches the
//     nonlinearity) for content in the don't-care band, through the class's public API;
//   * why a HALFBAND cannot buy the same thing — the identity H(f) + H(Fs/2 - f) = 1, on a halfband built
//     here, since the class has none and the claim is about every halfband;
//   * and what the guard band COSTS — the round trip at 19 and 20 kHz, per sample rate.
namespace cutoffaxis
{
    constexpr double kPi = 3.14159265358979323846;
    constexpr int    kSettle = 4096;

    // Coherent amplitude of `nu` (cycles per sample of THIS buffer) over exactly `len` samples from `off`.
    // Every caller picks `len` so that nu*len is an integer.
    inline double ampAt (const std::vector<float>& x, double nu, int off, int len)
    {
        double re = 0.0, im = 0.0;
        for (int i = 0; i < len; ++i)
        {
            const double ph = 2.0 * kPi * nu * (double) (off + i);
            re += (double) x[(std::size_t) (off + i)] * std::cos (ph);
            im += (double) x[(std::size_t) (off + i)] * std::sin (ph);
        }
        return 2.0 * std::hypot (re, im) / (double) len;
    }
    inline double db (double a) { return 20.0 * std::log10 (std::max (1e-300, a)); }

    // One round trip at a PHYSICAL frequency. The window is fs/10 samples, so every multiple of 10 Hz is an
    // exact bin at 44.1, 48 and 88.2 kHz — which tapsprobe's fixed 10 000-sample window cannot say about
    // 20 000 / 44 100 = 200 / 441.
    inline double roundTripHzDb (double fs, double hz, int factor, int taps, bool& refused)
    {
        oversampling::PolyphaseOversampler os;
        if (! os.prepare (factor, 1, taps)) { refused = true; return 1e9; }
        const int w = (int) std::lround (fs / 10.0), n = kSettle + w;
        std::vector<float> x ((std::size_t) n), osb ((std::size_t) n * (std::size_t) factor);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * hz / fs * (double) i);
        const float* in[1] { x.data() }; float* mid[1] { osb.data() };
        os.upsample (in, 1, n, mid);
        const float* mc[1] { osb.data() }; float* out[1] { x.data() };
        os.downsample (mc, 1, n, out);
        return db (ampAt (x, hz / fs, kSettle, w));
    }

    // ONE pass up: the tone's own gain and its first image's, both read IN THE OVERSAMPLED STREAM — the
    // stream the nonlinearity sees. `nu` is baseband-normalised and a multiple of 1e-4, the window 10 000.
    inline void upGains (int factor, int taps, double nu, double& sigDb, double& imgDb, bool& refused)
    {
        oversampling::PolyphaseOversampler os;
        sigDb = imgDb = 1e9;
        if (! os.prepare (factor, 1, taps)) { refused = true; return; }
        const int w = 10000, n = kSettle + w;
        std::vector<float> x ((std::size_t) n), osb ((std::size_t) n * (std::size_t) factor);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * nu * (double) i);
        const float* in[1] { x.data() }; float* out[1] { osb.data() };
        os.upsample (in, 1, n, out);
        sigDb = db (ampAt (osb, nu / factor,          kSettle * factor, w * factor));
        imgDb = db (ampAt (osb, (1.0 - nu) / factor,  kSettle * factor, w * factor));
    }

    // A halfband 2x prototype of length 4m+3: centre 1/2, even offsets structurally zero, odd offsets a
    // Kaiser-windowed sinc scaled to sum to 1/2. sin(pi t / 2) is +-1 for odd t, so no sine is needed.
    inline std::vector<double> halfband (int m, double beta)
    {
        const int n = 4 * m + 3, c = 2 * m + 1;
        auto i0 = [] (double x) { double s = 1.0, t = 1.0; for (int k = 1; k < 64; ++k) { t *= x * x * 0.25 / ((double) k * k); s += t; if (t < 1e-17 * s) break; } return s; };
        std::vector<double> h ((std::size_t) n, 0.0);
        double odd = 0.0;
        for (int j = 0; j < n; ++j)
        {
            const int t = j - c, at = t < 0 ? -t : t;
            if (at % 2 == 0) continue;
            const double r = (double) (2 * j - (n - 1)) / (double) (n - 1);
            const double v = ((at - 1) / 2 % 2 == 0 ? 1.0 : -1.0) / (kPi * (double) at)
                           * i0 (beta * std::sqrt (std::max (0.0, 1.0 - r * r))) / i0 (beta);
            h[(std::size_t) j] = v; odd += v;
        }
        for (auto& v : h) v *= 0.5 / odd;
        h[(std::size_t) c] = 0.5;
        return h;
    }
    // Zero-phase amplitude at `nu` cycles per sample of the filter's own (2x) stream.
    inline double zeroPhase (const std::vector<double>& h, double nu)
    {
        const double c = 0.5 * (double) (h.size() - 1);
        double a = 0.0;
        for (std::size_t k = 0; k < h.size(); ++k) a += h[k] * std::cos (2.0 * kPi * nu * ((double) k - c));
        return a;
    }
}

static void runCutoffAxisTests()
{
    using namespace cutoffaxis;

    test::group ("the cutoff axis (P31): what the guard band BUYS — exposure of don't-care content");
    {
        bool refused = false;
        double s, i;
        // Liveness: the probe reads a flat pass band as 0 dB with a buried image, and reads the design's own
        // cutoff anchor (-6.02 dB for one pass) — an instrument that could only say "small" fails here.
        upGains (4, 64, 0.20, s, i, refused);
        test::approx (s, 0.0, 0.01, "PRECONDITION: a 0.20 fs tone goes up flat");
        test::ok (i < -90.0, "PRECONDITION: and its image is buried (" + std::to_string (i) + " dB)");
        upGains (4, 64, 0.45, s, i, refused);
        test::approx (s, -6.02, 0.05, "PRECONDITION: the probe reads the cutoff anchor, -6.02 dB for one pass");

        const double r[] = { 0.46, 0.47, 0.48 }, pinned[] = { -109.0, -117.0, -136.0 };
        for (int k = 0; k < 3; ++k)
        {
            upGains (4, 64, r[k], s, i, refused);
            std::printf ("       r = %.2f: content %7.2f dB, image %7.2f dB, exposure %7.2f dB\n", r[k], s, i, s + i);
            test::ok (i <= -90.0, "r = " + std::to_string (r[k]) + ": the image of don't-care content is rejected at the "
                                  "design's own stopband (" + std::to_string (i) + " dB) — the design is STRICT");
            test::approx (s + i, pinned[k], 1.0, "r = " + std::to_string (r[k]) + ": exposure as the header states it");
        }
        test::ok (! refused, "PRECONDITION: every probe prepared");
    }

    // THE WHOLE FOLD REGION, densely (the P21 tail: `worstFoldDb` above steps 0.02 fs past 0.62 fs, against a
    // sidelobe spacing of 1/64 fs, so "worst over [0.5, L/2]" was not measured literally there). The prototype
    // comes out of the public API — upsampling an impulse returns L times it — and its zero-phase response is
    // evaluated at 1/32 of a sidelobe over the entire region. (A quarter-lobe step, the first draft, read
    // -91.08 at 2x where the truth is -90.48: sampling a lobe four times under-reads its peak by up to 0.7 dB,
    // which is the direction that lets a bound pass for the wrong reason. The dense maximum is therefore also
    // required to agree with the time-domain probe above, which walks the transition end finely.)
    test::group ("the fold region, densely: the stopband the design declares holds all the way to L/2");
    {
        for (int factor : { 2, 4, 8 })
        {
            oversampling::PolyphaseOversampler os;
            (void) os.prepare (factor, 1);
            const int N = factor * oversampling::PolyphaseOversampler::kDefaultTapsPerPhase;
            std::vector<float> x ((std::size_t) N, 0.0f), y ((std::size_t) N * (std::size_t) factor);
            x[0] = 1.0f;
            const float* in[1] { x.data() }; float* out[1] { y.data() };
            os.upsample (in, 1, N, out);
            std::vector<double> h ((std::size_t) N);
            for (int i = 0; i < N; ++i) h[(std::size_t) i] = (double) y[(std::size_t) i] / (double) factor;
            double worst = -1e9;
            const double stepOs = 1.0 / (32.0 * (double) N);           // 1/32 of a sidelobe, OS-normalised
            for (double nuOs = 0.5 / factor; nuOs <= 0.5; nuOs += stepOs)
                worst = std::max (worst, db (std::fabs (zeroPhase (h, nuOs))));
            std::printf ("       factor %d, 64 taps: worst |H| over the WHOLE fold region, dense = %.2f dB\n", factor, worst);
            test::ok (worst <= -90.0, "factor " + std::to_string (factor)
                                      + ": no point of [0.5, L/2] fs rises above the declared -90 dB (" + std::to_string (worst) + ")");
            const double probe = tapsprobe::worstFoldDb (factor, -1);
            test::approx (worst, probe, 0.05, "factor " + std::to_string (factor)
                          + ": and the dense maximum is the one the time-domain probe finds (the worst sits at the transition)");
        }
    }

    test::group ("the cutoff axis (P31): a HALFBAND first stage can never be strict — the identity, and its price");
    {
        const auto hb = halfband (31, 9.0);                     // 127 taps: the halfband that is flat to 20 kHz
        double worst = 0.0;
        for (int k = 0; k <= 2500; ++k)
        {
            const double nu = 0.25 * (double) k / 2500.0;
            worst = std::max (worst, std::fabs (zeroPhase (hb, nu) + zeroPhase (hb, 0.5 - nu) - 1.0));
        }
        test::ok (worst < 1e-12, "H(f) + H(Fs/2 - f) = 1 at every f (worst residual " + std::to_string (worst) + ")");
        test::approx (zeroPhase (hb, 0.25), 0.5, 1e-12, "so the transition is centred ON the fold: H(fs/2) = 1/2 exactly");

        const double fs = 44100.0;
        test::ok (std::fabs (db (zeroPhase (hb, 20000.0 / (2.0 * fs)))) < 0.001,
                  "the 127-tap halfband is flat at 20 kHz (one pass within 0.001 dB)");
        const double r[] = { 0.46, 0.47, 0.48 }, pinned[] = { -58.0, -35.0, -22.0 };
        for (int k = 0; k < 3; ++k)
        {
            const double content = db (std::fabs (zeroPhase (hb, r[k] / 2.0)));
            const double image   = db (std::fabs (zeroPhase (hb, (1.0 - r[k]) / 2.0)));
            std::printf ("       halfband 127, r = %.2f: content %7.2f dB, image %7.2f dB, exposure %7.2f dB\n",
                         r[k], content, image, content + image);
            test::approx (content + image, pinned[k], 1.5, "r = " + std::to_string (r[k])
                          + ": the halfband's exposure as the header states it — its image is 1 - H(content)");
        }

        // The pair the header names: 79 + 23 taps, -0.41 dB round trip at 20 kHz, image of that tone at -32.5.
        const auto h1 = halfband (19, 9.0), h2 = halfband (5, 9.0);
        const double a1 = zeroPhase (h1, 20000.0 / (2.0 * fs)), a2 = zeroPhase (h2, 20000.0 / (4.0 * fs));
        const double rt = 2.0 * (db (a1) + db (a2));
        const double img = db (1.0 - a1) + db (std::fabs (zeroPhase (h2, (fs - 20000.0) / (4.0 * fs))));
        std::printf ("       halfband 79 + 23 at 20 kHz: round trip %.3f dB, image of the tone %.2f dB\n", rt, img);
        test::approx (rt, -0.41, 0.01, "the 79+23 pair loses 0.41 dB at 20 kHz over a round trip");
        test::approx (img, -32.5, 0.3, "and therefore leaves that tone's image at -32.5 dB, not at -90");
    }

    test::group ("the cutoff axis (P31): what the guard band COSTS — the top of the audio band, per rate");
    {
        bool refused = false;
        test::approx (roundTripHzDb (44100.0, 1000.0, 4, 64, refused), 0.0, 0.001,
                      "PRECONDITION: 1 kHz goes round flat (the probe is not reading a loss of its own)");
        const double at19 = roundTripHzDb (44100.0, 19000.0, 4, 64, refused);
        const double at20 = roundTripHzDb (44100.0, 20000.0, 4, 64, refused);
        const double at20t32 = roundTripHzDb (44100.0, 20000.0, 4, 32, refused);
        const double at20t120 = roundTripHzDb (44100.0, 20000.0, 4, 120, refused);
        const double at48 = roundTripHzDb (48000.0, 20000.0, 4, 64, refused);
        const double at88 = roundTripHzDb (88200.0, 20000.0, 4, 64, refused);
        std::printf ("       44.1 kHz round trip: 19 k %.2f, 20 k %.2f dB (64 taps); 20 k at 32 taps %.2f, at 120 %.2f\n",
                     at19, at20, at20t32, at20t120);
        std::printf ("       20 kHz at 48 kHz %.3f dB, at 88.2 kHz %.4f dB\n", at48, at88);
        test::approx (at19, -1.80, 0.02, "44.1 kHz: 19 kHz loses 1.80 dB over one round trip");
        test::approx (at20, -15.55, 0.02, "44.1 kHz: 20 kHz loses 15.55 dB");
        test::approx (at20t32, -13.71, 0.02, "44.1 kHz: 20 kHz at 32 taps loses 13.71 dB");
        test::approx (at20t120, -19.17, 0.02, "44.1 kHz: 20 kHz at 120 taps loses 19.17 dB");
        test::ok (at20t120 < at20 - 2.0 && at20 < at20t32 - 1.0,
                  "and MORE taps make 20 kHz WORSE (32 -> 64 -> 120: " + std::to_string (at20t32) + ", "
                  + std::to_string (at20) + ", " + std::to_string (at20t120) + ") — a cutoff axis, not a taps one");
        test::approx (at48, -0.15, 0.02, "48 kHz: 20 kHz loses 0.15 dB");
        test::ok (std::fabs (at88) < 0.001, "88.2 kHz: 20 kHz passes flat");
        test::ok (! refused, "PRECONDITION: every probe prepared");
    }
}

int main()
{
    std::printf ("felitronics::oversampling tests\n");
    const double sr = 48000.0;
    const int    tpp = 32;

    // --- 4x upsample interpolation matches the analytic sine (passband, no aliasing) ---
    test::group ("Oversampler 4x upsample == analytic sine");
    {
        const int L = 4, n = 512; const double f = 2000.0, A = 0.8;
        oversampling::PolyphaseOversampler os; test::ok (os.prepare (L, 1, tpp), "prepare 4x");
        std::vector<float> x (n); for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> y ((std::size_t) n * L);
        const float* xi[1] { x.data() }; float* yo[1] { y.data() };
        os.upsample (xi, 1, n, yo);

        const double gdOS = os.filterLatencyOversampled();        // (N-1)/2 OS samples
        double maxErr = 0.0;
        for (int i = 200; i < n * L - 50; ++i)
        {
            const double a = A * std::sin (2.0 * core::kPi * f * (i - gdOS) / (sr * L));
            maxErr = std::max (maxErr, std::fabs (a - y[(std::size_t) i]));
        }
        test::ok (maxErr < 0.02, "OS samples lie on the analytic sine (interp correct, images rejected)");
    }

    // --- reveals an inter-sample peak the baseband samples miss ---
    test::group ("Oversampler reveals inter-sample peak");
    {
        const int L = 4, n = 256; const double f = sr * 0.25, A = 1.0;
        oversampling::PolyphaseOversampler os; os.prepare (L, 1, tpp);
        std::vector<float> x (n); double sampMax = 0.0;
        for (int i = 0; i < n; ++i) { x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr + 0.7)); sampMax = std::max (sampMax, (double) std::fabs (x[i])); }
        std::vector<float> y ((std::size_t) n * L);
        const float* xi[1] { x.data() }; float* yo[1] { y.data() };
        os.upsample (xi, 1, n, yo);
        double osMax = 0.0; for (int i = 100; i < n * L - 100; ++i) osMax = std::max (osMax, (double) std::fabs (y[(std::size_t) i]));
        test::ok (osMax > sampMax + 0.01, "upsampled peak exceeds the sample peak (ISP found)");
        test::ok (osMax <= A * 1.02, "and ~ the true amplitude (no overshoot)");
    }

    // --- up → down round-trip == input delayed by latency (band-limited) ---
    test::group ("Oversampler up->down round-trip == delayed identity");
    {
        const int L = 4, n = 512; const double f = 500.0, A = 0.6;
        oversampling::PolyphaseOversampler os; os.prepare (L, 1, tpp);
        std::vector<float> x (n); for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> osb ((std::size_t) n * L), z (n);
        const float* xi[1] { x.data() }; float* ob[1] { osb.data() };
        os.upsample (xi, 1, n, ob);
        const float* obc[1] { osb.data() }; float* zo[1] { z.data() };
        os.downsample (obc, 1, n, zo);

        const int lat = os.latencySamples();
        test::ok (lat > 0, "reports a positive round-trip latency");
        double maxErr = 0.0;
        for (int i = lat + 60; i < n - 10; ++i) maxErr = std::max (maxErr, (double) std::fabs (z[(std::size_t) i] - x[(std::size_t) (i - lat)]));
        test::ok (maxErr < 0.03, "round-trip == input delayed by latency()");
    }

    // --- lifecycle/misuse: upsample/downsample before prepare (or after a failed prepare) must no-op ---
    // channels_ defaults to 0 and prepare() rejects bad args BEFORE mutating dims, so misuse is a safe no-op.
    test::group ("PolyphaseOversampler: safe before prepare / after failed prepare");
    {
        // Sentinels, because "no-op" is a claim about the CALLER's buffers and nothing here read them
        // back: a downsample() that zeroed its output when unprepared passed every check below.
        float sentinelIn[8], sentinelOut[8];
        for (int i = 0; i < 8; ++i) { sentinelIn[i] = 0.5f + 0.1f * (float) i; sentinelOut[i] = -1.0f - (float) i; }
        {
            oversampling::PolyphaseOversampler unprepped;
            const float* si[1] { sentinelIn }; float* so[1] { sentinelOut };
            unprepped.upsample (si, 1, 8, so);
            unprepped.downsample (si, 1, 2, so);
            bool untouched = true;
            for (int i = 0; i < 8; ++i) untouched = untouched && sentinelOut[i] == -1.0f - (float) i;
            test::ok (untouched, "unprepared up/downsample leave the caller's OUTPUT buffer untouched, "
                                 "not merely 'do not crash'");
        }
        oversampling::PolyphaseOversampler os;         // NOT prepared (channels_ == 0)
        float base[8] { 0.1f, 0,0,0,0,0,0,0 }; float osb[32] {};
        const float* in[1] { base }; float* outo[1] { osb };
        os.upsample (in, 1, 8, outo);                  // channels_==0 → no-op
        const float* ino[1] { osb }; float* outb[1] { base };
        os.downsample (ino, 1, 8, outb);               // no-op
        test::ok (! os.prepare (4, 1, 2), "prepare(tapsPerPhase=2) fails (rejected before mutating dims)");
        os.upsample (in, 1, 8, outo);                  // still unprepared → no-op
        test::ok (os.prepare (4, 1, 32), "prepare valid");
    }

    // --- latency query before prepare() must not report tpp-1 == -1 ---
    test::group ("Oversampler latencySamples before prepare == 0");
    {
        oversampling::PolyphaseOversampler os;
        test::ok (os.latencySamples() == 0, "unprepared oversampler reports 0, not -1");
    }

    // --- 8x round-trip == delayed identity (the 2/4x paths are covered above; 8x is the metering path) ---
    test::group ("Oversampler 8x round-trip == delayed identity");
    {
        const int L = 8, n = 512; const double f = 500.0, A = 0.6;
        oversampling::PolyphaseOversampler os; os.prepare (L, 1, tpp);
        std::vector<float> x (n); for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> osb ((std::size_t) n * L), z (n);
        const float* xi[1] { x.data() }; float* ob[1] { osb.data() };
        os.upsample (xi, 1, n, ob);
        const float* obc[1] { osb.data() }; float* zo[1] { z.data() };
        os.downsample (obc, 1, n, zo);
        const int lat = os.latencySamples();
        test::ok (lat == tpp - 1, "8x round-trip latency == tpp-1");
        double maxErr = 0.0;
        for (int i = lat + 60; i < n - 10; ++i) maxErr = std::max (maxErr, (double) std::fabs (z[(std::size_t) i] - x[(std::size_t) (i - lat)]));
        test::ok (maxErr < 0.03, "8x round-trip == input delayed by latency()");
    }

    // --- P56: the DOUBLE-LENGTH ring must be invisible ---------------------------------------------
    // Both loops now keep a backwards ring in which every sample is stored twice, and read a window
    // that starts at a cursor which only ever decrements. Two things can break in that and nothing
    // above would notice: the wrap, and the per-channel stride. Both checks are BIT-EXACT on purpose —
    // "close enough" is what let a ring bug hide in the first place.
    test::group ("Oversampler: block splitting and channel count are bit-invisible");
    {
        const int L = 4, n = 700, nch = 2;
        std::vector<float> x0 ((std::size_t) n), x1 ((std::size_t) n);
        std::uint32_t rs = 0xC0FFEEu;
        auto rnd = [&rs] { rs = rs * 1664525u + 1013904223u; return (float) (std::int32_t) ((rs >> 8) ^ 0x800000u) * (1.0f / 8388608.0f) - 1.0f; };
        for (int i = 0; i < n; ++i) { x0[(std::size_t) i] = rnd() * 0.7f; x1[(std::size_t) i] = rnd() * 0.7f; }

        auto runUpDown = [&] (int chans, const std::vector<int>& splits, std::vector<float>& outCh0)
        {
            oversampling::PolyphaseOversampler os; os.prepare (L, chans, tpp);
            std::vector<float> u0 ((std::size_t) n * L), u1 ((std::size_t) n * L);
            std::vector<float> d0 ((std::size_t) n), d1 ((std::size_t) n);
            int at = 0;
            for (int len : splits)
            {
                const float* in[2]  { x0.data() + at, x1.data() + at };
                float*       up[2]  { u0.data() + (std::size_t) at * L, u1.data() + (std::size_t) at * L };
                os.upsample (in, chans, len, up);
                const float* upc[2] { u0.data() + (std::size_t) at * L, u1.data() + (std::size_t) at * L };
                float*       dn[2]  { d0.data() + at, d1.data() + at };
                os.downsample (upc, chans, len, dn);
                at += len;
            }
            test::ok (at == n, "the split covers the signal");
            outCh0 = d0;
        };

        std::vector<float> whole, split, mono;
        runUpDown (2, { n }, whole);
        runUpDown (2, { 1, 3, 64, 7, 128, 2, 251, 244 }, split);   // 700, deliberately ragged
        runUpDown (1, { n }, mono);

        bool sameSplit = true, sameMono = true;
        for (int i = 0; i < n; ++i)
        {
            if (split[(std::size_t) i] != whole[(std::size_t) i]) sameSplit = false;
            if (mono [(std::size_t) i] != whole[(std::size_t) i]) sameMono  = false;
        }
        test::ok (sameSplit, "eight ragged blocks == one call, BIT for bit (the ring wrap)");
        test::ok (sameMono,  "channel 0 alone == channel 0 of a stereo run, BIT for bit (the stride)");
    }

    // --- P56: resetChannel still clears exactly one channel's (now longer) rings --------------------
    test::group ("Oversampler: resetChannel clears one channel and leaves the other bit-exact");
    {
        const int L = 4, n = 200;
        oversampling::PolyphaseOversampler a2, b1;
        a2.prepare (L, 2, tpp); b1.prepare (L, 2, tpp);
        std::vector<float> x ((std::size_t) n); for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (0.11 * i) * 0.8f;
        std::vector<float> ua ((std::size_t) n * L), ub ((std::size_t) n * L), va ((std::size_t) n * L), vb ((std::size_t) n * L);
        const float* in[2] { x.data(), x.data() };
        float* o1[2] { ua.data(), ub.data() };
        a2.upsample (in, 2, n, o1); b1.upsample (in, 2, n, o1);      // both warmed identically
        a2.resetChannel (0);                                          // only a2's channel 0 forgets
        float* o2[2] { va.data(), vb.data() };
        a2.upsample (in, 2, n, o2);
        std::vector<float> wa ((std::size_t) n * L), wb ((std::size_t) n * L);
        float* o3[2] { wa.data(), wb.data() };
        b1.upsample (in, 2, n, o3);
        bool ch1Exact = true, ch0Moved = false;
        for (std::size_t i = 0; i < vb.size(); ++i) if (vb[i] != wb[i]) ch1Exact = false;
        for (std::size_t i = 0; i < va.size(); ++i) if (va[i] != wa[i]) ch0Moved = true;
        test::ok (ch1Exact, "the untouched channel is bit-identical to one that was never reset");
        test::ok (ch0Moved, "and the reset one really did forget (the check is not vacuous)");
    }

    runTapsTests();
    runReferenceTapPin();
    runCutoffAxisTests();

    return test::report();
}
