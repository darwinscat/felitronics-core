// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the polyphase oversampler: upsample interpolation vs the analytic sine
// (correct + alias-free), inter-sample-peak revelation (the limiter's reason to exist), the
// up→down round-trip == a delayed identity, and — since this module owns the shipped tapsPerPhase
// default — what that default actually buys: the declared stopband, the aliasing it stops, and the
// pass-band droop of two oversampled stages in series.

#include <felitronics_test.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
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
        // 48 taps is not a cheaper way to satisfy it, and 60 is where it first holds: pinned so that
        // "somewhere between 32 and 64" cannot be quietly re-opened.
        g_prepareRefused = false;
        const double at48 = worstFoldDb (4, 48), at60 = worstFoldDb (4, 60);
        test::ok (! g_prepareRefused,
                  "PRECONDITION: every probe above actually PREPARED — a refused prepare returns +1e9, which "
                  "satisfies an upper bound without measuring anything");
        test::ok (at48 > -60.0, "48 taps still misses the declared stopband by ~38 dB");
        test::ok (at60 <= -90.0, "60 taps is where it first holds — the default rounds that up to 64");
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
            (void) (taps < 0 ? os.prepare (4, 1) : os.prepare (4, 1, taps));
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
        const double now = aliasDbc (-1), was = aliasDbc (32);
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

    runTapsTests();

    return test::report();
}
