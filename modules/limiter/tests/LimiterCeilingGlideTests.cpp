// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// TruePeakLimiter's CEILING GLIDE (TruePeakLimiter::kCeilingGlideMs). A lower ceiling falls to its target along a
// one-pole per oversampled sample; a higher one lands at once; the first write of a stream snaps. The groups:
//
//   * THE BOUND HOLDS PER SAMPLE AGAINST THE CEILING IN FORCE: on the limiter's own oversampled grid, every emitted
//     sample is at most 10^(c_i/20), with c_i replayed here from the stated recursion — read through the taps (the
//     linked peak entering the window and the gain reduction applied D samples later), with the peak clipper off
//     and on, single and dual release.
//   * IT LANDS EXACTLY, UP IS INSTANT, AND THE FIRST WRITE SNAPS (prepare -> set equals set -> prepare).
//   * CLOCKED BY SAMPLES (law 8a): ceiling writes at fixed positions render the same bits under any cut.
//   * IT NO LONGER CLICKS: 0 -> -18 dBTP on a -12 dBFS tone, the second difference of the output against the step's.

#include <felitronics_test.h>
#include <felitronics/limiter/TruePeakLimiter.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using limiter::TruePeakLimiter;
using limiter::TruePeakLimiterParams;

namespace
{
constexpr double kFs = 48000.0, kPi = 3.14159265358979323846;
std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }

TruePeakLimiterParams P (double ceil, bool clip = false, bool dual = false)
{
    TruePeakLimiterParams p; p.ceilingDbTp = ceil; p.releaseMs = 60.0; p.peakClip = clip; p.dualRelease = dual;
    return p;
}

std::vector<float> tone (int n, double amp) { std::vector<float> x ((std::size_t) n); for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (amp * std::sin (2 * kPi * 227.3 * i / kFs)); return x; }

struct Write { int at; TruePeakLimiterParams p; };

// Mono through a fresh limiter; `p0` before the first sample; writes at their positions; calls cut by `cut(k)`.
template <class Cut>
std::vector<float> render (const std::vector<float>& x, const TruePeakLimiterParams& p0, const std::vector<Write>& ws, Cut cut,
                           std::vector<float>* gr = nullptr, std::vector<float>* pk = nullptr, bool setBeforePrepare = true)
{
    std::vector<float> y = x;
    TruePeakLimiter lim;
    if (setBeforePrepare) lim.setParams (p0);
    if (! lim.prepare (kFs, 512, 1)) return {};
    if (! setBeforePrepare) lim.setParams (p0);
    const int n = (int) x.size(), F = lim.oversampleFactor();
    if (gr) gr->assign ((std::size_t) n * (std::size_t) F, 0.0f);
    if (pk) pk->assign ((std::size_t) n * (std::size_t) F, 0.0f);
    std::size_t w = 0;
    for (int off = 0, k = 0; off < n; ++k)
    {
        while (w < ws.size() && ws[w].at <= off) lim.setParams (ws[w++].p);
        int m = std::min (cut (k), n - off);
        if (w < ws.size()) m = std::min (m, ws[w].at - off);
        float* io[1] { y.data() + off };
        limiter::TruePeakLimiterTap t;
        if (gr) { t.gainReductionDb = gr->data() + (std::size_t) off * (std::size_t) F; t.linkedPeakLin = pk->data() + (std::size_t) off * (std::size_t) F; t.capacity = m * F; }
        felitronics::test::run (lim.process (io, 1, m, t));
        off += m;
    }
    return y;
}
auto fixed (int b) { return [b] (int) { return b; }; }
auto ragged (std::uint32_t seed)
{
    return [seed] (int k)
    {
        std::uint32_t s = seed + 2654435761u * (std::uint32_t) (k + 1);
        s ^= s >> 13; s *= 0x5bd1e995u; s ^= s >> 15;
        return (s % 7u == 0u) ? 0 : 1 + (int) (s % 900u);
    };
}
} // namespace

static void testBoundHoldsAgainstTheCeilingInForce()
{
    group ("the bound holds per oversampled sample against the ceiling in force, all the way down the glide");
    const int n = 24000, at = 12000 + 17;
    const std::vector<float> x = tone (n, 0.5);
    for (const bool clip : { false, true })
        for (const bool dual : { false, true })
        {
            std::vector<float> gr, pk;
            (void) render (x, P (0.0, clip, dual), { { at, P (-18.0, clip, dual) } }, fixed (256), &gr, &pk);
            TruePeakLimiter probe;
            ok (probe.prepare (kFs, 512, 1), "PRECONDITION: prepare");
            const int F = probe.oversampleFactor(), D = probe.lookaheadSamples() * F;
            // The ceiling in force, replayed from the stated recursion, per oversampled sample.
            const double coef = std::exp (-1.0 / (TruePeakLimiter::kCeilingGlideMs * 0.001 * kFs * (double) F));
            std::vector<double> cs ((std::size_t) (n * F));
            double c = 0.0, worst = -1e9;
            int glideSamples = 0;
            for (int i = 0; i < n * F; ++i)
            {
                const double target = i >= at * F ? -18.0 : 0.0;   // the write lands before os sample at*F steps
                if (c > target) { const double d = coef * (c - target); c = d < 1.0e-9 ? target : target + d; ++glideSamples; }
                cs[(std::size_t) i] = c;
                if (i < D) continue;
                // The emitted sample at i ENTERED the window D samples earlier — clipped, with the clipper on, at the
                // level riding the ceiling in force THEN (default offset +1 dB, no knee): its magnitude is at most
                // min(peak, clip level). The tap's peak is taken before the clip.
                double in = (double) pk[(std::size_t) (i - D)];
                if (clip) in = std::min (in, std::pow (10.0, (cs[(std::size_t) (i - D)] + 1.0) / 20.0));
                const double out = in * std::pow (10.0, (double) gr[(std::size_t) i] / 20.0);
                const double over = 20.0 * std::log10 (std::max (out, 1e-30)) - c;
                worst = std::max (worst, over);
            }
            // 1e-5 dB is a few float ulps: the limiter's gain is a float of a float dB, and its static ceiling has always
            // carried that rounding (the code-review round measured 1e-6 dB over, contraction on and off alike).
            ok (glideSamples > F * 48 && worst <= 1.0e-5, std::string ("clip ") + (clip ? "on" : "off") + ", dual " + (dual ? "on" : "off")
                + ": every emitted sample within its own ceiling (worst " + std::to_string (worst) + " dB over, "
                + std::to_string (glideSamples) + " gliding os samples)");
        }
}

static void testLandsUpInstantFirstWriteSnaps()
{
    group ("it lands exactly, a higher ceiling is instant, and the first write of a stream snaps");
    TruePeakLimiter lim;
    lim.setParams (P (0.0));
    ok (lim.prepare (kFs, 512, 1), "PRECONDITION: prepare");
    std::vector<float> x = tone (48000, 0.5);
    float* io[1] { x.data() };
    felitronics::test::run (lim.process (io, 1, 4800));
    lim.setParams (P (-18.0));
    ok (lim.ceilingNowDbTp() == 0.0 && lim.effectiveCeilingDbTp() == -18.0, "a lower ceiling is a target: nothing moves before a sample");
    float* io2[1] { x.data() + 4800 };
    felitronics::test::run (lim.process (io2, 1, 4));
    ok (lim.ceilingNowDbTp() < 0.0 && lim.ceilingNowDbTp() > -18.0, "it glides (" + std::to_string (lim.ceilingNowDbTp()) + " dBTP after 4 samples)");
    float* io3[1] { x.data() + 4804 };
    felitronics::test::run (lim.process (io3, 1, 4800));
    ok (lim.ceilingNowDbTp() == -18.0, "and lands on the target exactly");
    lim.setParams (P (-6.0));
    ok (lim.ceilingNowDbTp() == -6.0, "a higher ceiling lands at once (the release makes that move smooth)");

    // prepare -> set equals set -> prepare, and a reset mid-glide lands it.
    const std::vector<float> y = tone (12000, 0.5);
    const std::vector<float> a = render (y, P (-9.0, true), {}, fixed (512), nullptr, nullptr, true);
    const std::vector<float> b = render (y, P (-9.0, true), {}, fixed (512), nullptr, nullptr, false);
    long long d = 0; for (std::size_t i = 0; i < a.size(); ++i) d += bits (a[i]) != bits (b[i]);
    ok (! a.empty() && d == 0, "set -> prepare and prepare -> set render the same bits");
    lim.setParams (P (-24.0));
    felitronics::test::run (lim.process (io, 1, 2));
    ok (lim.ceilingNowDbTp() > -24.0, "PRECONDITION: gliding");
    lim.reset();
    ok (lim.ceilingNowDbTp() == -24.0, "reset() lands the glide");
}

// A CLOCK-ONLY CALL SPENDS THE GLIDE: a lower ceiling written during a gap has landed when audio returns after a gap
// longer than the glide (found by the code-review round: it stood frozen at the old ceiling through the whole gap).
static void testGapSpendsTheGlide()
{
    group ("a clock-only gap spends the ceiling glide as audio time");
    TruePeakLimiter lim;
    lim.setParams (P (0.0, true));
    ok (lim.prepare (kFs, 512, 2), "PRECONDITION: prepare");
    std::vector<float> a = tone (4800, 0.5), b = a;
    float* io[2] { a.data(), b.data() };
    felitronics::test::run (lim.process (io, 2, 4800));
    felitronics::test::run (lim.process (nullptr, 0, 100));       // the gap starts
    lim.setParams (P (-18.0, true));
    felitronics::test::run (lim.process (nullptr, 0, 4800));      // 100 ms of gap, the glide is 2 ms
    ok (lim.ceilingNowDbTp() == -18.0, "the ceiling has landed during the gap (" + std::to_string (lim.ceilingNowDbTp()) + ")");
    TruePeakLimiter cut;
    cut.setParams (P (0.0, true));
    ok (cut.prepare (kFs, 512, 2), "PRECONDITION: prepare");
    std::vector<float> c = tone (4800, 0.5), d = c;
    float* io2[2] { c.data(), d.data() };
    felitronics::test::run (cut.process (io2, 2, 4800));
    felitronics::test::run (cut.process (nullptr, 0, 100));
    cut.setParams (P (-18.0, true));
    for (int k = 0; k < 10; ++k) felitronics::test::run (cut.process (nullptr, 0, 3));   // 30 samples of gap, cut
    TruePeakLimiter whole;
    whole.setParams (P (0.0, true));
    ok (whole.prepare (kFs, 512, 2), "PRECONDITION: prepare");
    std::vector<float> e = tone (4800, 0.5), f = e;
    float* io3[2] { e.data(), f.data() };
    felitronics::test::run (whole.process (io3, 2, 4800));
    felitronics::test::run (whole.process (nullptr, 0, 100));
    whole.setParams (P (-18.0, true));
    felitronics::test::run (whole.process (nullptr, 0, 30));
    ok (cut.ceilingNowDbTp() == whole.ceilingNowDbTp() && cut.ceilingNowDbTp() > -18.0,
        "a gap cut in pieces lands the glide where one call does (" + std::to_string (cut.ceilingNowDbTp()) + ")");
}

static void testClockedBySamples()
{
    group ("law 8a — ceiling writes at fixed positions render the same bits under any cut");
    const int n = 20000;
    const std::vector<float> x = tone (n, 0.6);
    const std::vector<Write> ws { { 3001, P (-12.0, true) }, { 3100, P (-20.0, true) }, { 9000, P (-3.0, true) }, { 12345, P (-15.0, true, true) } };
    const std::vector<float> ref = render (x, P (0.0, true), ws, fixed (n));
    for (const auto& cut : { std::function<int (int)> (fixed (1)), std::function<int (int)> (fixed (64)),
                             std::function<int (int)> (ragged (1)), std::function<int (int)> (ragged (2)) })
    {
        const std::vector<float> got = render (x, P (0.0, true), ws, cut);
        long long d = 0; for (std::size_t i = 0; i < ref.size(); ++i) d += bits (ref[i]) != bits (got[i]);
        ok (d == 0, "bit-identical (" + std::to_string (d) + " differ)");
    }
}

static void testNoClick()
{
    group ("it no longer clicks — 0 -> -18 dBTP on a -12 dBFS tone");
    const int n = 24000, at = 12000 + 17;
    const std::vector<float> x = tone (n, 0.25);
    const std::vector<float> y = render (x, P (0.0), { { at, P (-18.0) } }, fixed (128));
    auto maxD2 = [&] (const std::vector<float>& v, int a, int b)
    {
        double m = 0.0;
        for (int i = std::max (a, 2); i < b; ++i) m = std::max (m, std::fabs ((double) v[(std::size_t) i] - 2.0 * v[(std::size_t) i - 1] + v[(std::size_t) i - 2]));
        return m;
    };
    const double edge = maxD2 (y, at, at + 2400);
    std::printf ("    edge %.1f dBFS\n", 20.0 * std::log10 (edge));
    ok (edge < 1.0e-3, "the edge stays under -60 dBFS of Δ² (the step read -29.9 through the chain)");
}

int main()
{
    std::printf ("felitronics::limiter — the ceiling glide\n");
    testBoundHoldsAgainstTheCeilingInForce();
    testLandsUpInstantFirstWriteSnaps();
    testGapSpendsTheGlide();
    testClockedBySamples();
    testNoClick();
    return felitronics::test::report();
}
