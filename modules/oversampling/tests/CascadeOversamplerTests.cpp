// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Self-tests for oversampling::CascadeOversampler (P31) and the Topology switch in oversampling::Oversampler.
//
// THE INSTRUMENT. A cascade of (zero-stuff, filter) stages is, by the noble identities, ONE zero-stuff by F
// followed by one composite filter G; so upsampling a unit impulse returns F*G exactly, and |G| at every
// frequency above fs/2 IS the image rejection of some content below fs/2 (k*fs +- a covers the whole
// range). The decimator is the mirror: one composite D, then keep one sample in F, so an impulse at
// oversampled phase p comes out as D[m*F - p]. Both composites are recovered through the public API only and
// read with a zero-padded FFT — dense enough (8x) to land within a fraction of a lobe of every sidelobe
// peak, and it includes the bin AT fs/2, where the design puts its -90 dB point.

#include <felitronics_test.h>
#include <alloc_counter.h>
#include <felitronics/oversampling/CascadeOversampler.h>
#include <felitronics/oversampling/Oversampler.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using oversampling::CascadeOversampler;

namespace
{
    constexpr double kPi = 3.14159265358979323846;

    void fft (std::vector<std::complex<double>>& a)
    {
        const std::size_t n = a.size();
        for (std::size_t i = 1, j = 0; i < n; ++i)
        {
            std::size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap (a[i], a[j]);
        }
        for (std::size_t len = 2; len <= n; len <<= 1)
        {
            const double ang = -2.0 * kPi / (double) len;
            for (std::size_t i = 0; i < n; i += len)
                for (std::size_t k = 0; k < len / 2; ++k)
                {
                    const std::complex<double> w (std::cos (ang * (double) k), std::sin (ang * (double) k));
                    const auto u = a[i + k], v = a[i + k + len / 2] * w;
                    a[i + k] = u + v;
                    a[i + k + len / 2] = u - v;
                }
        }
    }

    // |FT| of a real sequence on an 8x zero-padded grid; bin b is b * rate / size.
    std::vector<double> magnitude (const std::vector<double>& h, std::size_t& size)
    {
        size = 1;
        while (size < h.size() * 8) size <<= 1;
        std::vector<std::complex<double>> a (size);
        for (std::size_t i = 0; i < h.size(); ++i) a[i] = h[i];
        fft (a);
        std::vector<double> m (size / 2 + 1);
        for (std::size_t i = 0; i < m.size(); ++i) m[i] = std::abs (a[i]);
        return m;
    }

    // The composite interpolation filter G (at F fs), normalised so that G(0) = 1.
    std::vector<double> upComposite (CascadeOversampler& co)
    {
        const int F = co.factor(), n = 2048;
        std::vector<float> x ((std::size_t) n, 0.0f), y ((std::size_t) n * (std::size_t) F);
        x[0] = 1.0f;
        const float* in[1] { x.data() }; float* out[1] { y.data() };
        co.reset();
        co.upsample (in, 1, n, out);
        std::vector<double> g (y.size());
        for (std::size_t i = 0; i < y.size(); ++i) g[i] = (double) y[i] / (double) F;
        return g;
    }

    // The composite decimation filter D (at F fs), recovered phase by phase; index m*F - p (+ F - 1).
    std::vector<double> downComposite (CascadeOversampler& co)
    {
        const int F = co.factor(), n = 1024;
        std::vector<double> d ((std::size_t) n * (std::size_t) F, 0.0);
        for (int p = 0; p < F; ++p)
        {
            std::vector<float> u ((std::size_t) n * (std::size_t) F, 0.0f), y ((std::size_t) n);
            u[(std::size_t) p] = 1.0f;
            const float* in[1] { u.data() }; float* out[1] { y.data() };
            co.reset();
            co.downsample (in, 1, n, out);
            for (int m = 0; m < n; ++m)
            {
                const long idx = (long) m * F - p + (F - 1);
                if (idx >= 0 && idx < (long) d.size()) d[(std::size_t) idx] = (double) y[(std::size_t) m];
            }
        }
        return d;
    }

    struct Reading { double strictDb = -1e9; double edgeDevDb = 0.0; double dc = 0.0; };

    // Worst |H| at and above fs/2 relative to DC, and the worst one-pass deviation over [0, band edge].
    Reading read (const std::vector<double>& h, double fs, int F, double edgeHz)
    {
        std::size_t size = 0;
        const auto m = magnitude (h, size);
        Reading r;
        r.dc = m[0];
        const double binHz = (double) F * fs / (double) size;
        for (std::size_t b = 0; b < m.size(); ++b)
        {
            const double f = (double) b * binHz, db = 20.0 * std::log10 (std::max (1e-300, m[b] / r.dc));
            if (f >= 0.5 * fs) r.strictDb = std::max (r.strictDb, db);
            if (f <= edgeHz)   r.edgeDevDb = std::max (r.edgeDevDb, std::fabs (db));
        }
        return r;
    }

    std::uint64_t fnv (const std::vector<float>& v)
    {
        std::uint64_t h = 1469598103934665603ull;
        for (float f : v)
        {
            std::uint32_t w; std::memcpy (&w, &f, 4);
            for (int k = 0; k < 4; ++k) { h ^= (w >> (8 * k)) & 0xffu; h *= 1099511628211ull; }
        }
        return h;
    }
}

static void runRuleTests()
{
    test::group ("CascadeOversampler: the design rule holds at every rate and factor (strict, flat to the band edge)");
    const double rates[] = { 8000.0, 22050.0, 32000.0, 44100.0, 46000.0, 48000.0, 50000.0, 64000.0,
                             88200.0, 96000.0, 176400.0, 192000.0, 384000.0, 768000.0 };
    double worstStrict = -1e9, worstEdge = 0.0, worstSym = 0.0;
    bool allPrepared = true, allSymmetric = true;
    for (double fs : rates)
        for (int F : { 2, 4, 8, 16 })
        {
            CascadeOversampler co;
            if (! co.prepare (fs, F, 1)) { allPrepared = false; continue; }
            const auto& d = co.design();
            const auto g = upComposite (co);
            const auto dn = downComposite (co);
            const Reading ru = read (g, fs, F, d.bandEdgeHz), rd = read (dn, fs, F, d.bandEdgeHz);
            worstStrict = std::max ({ worstStrict, ru.strictDb, rd.strictDb });
            worstEdge   = std::max ({ worstEdge, ru.edgeDevDb, rd.edgeDevDb });
            // Linear phase: the composite interpolator is a palindrome about its centre.
            std::size_t lo = 0, hi = g.size() - 1;
            while (lo < g.size() && std::fabs (g[lo]) < 1e-12) ++lo;
            while (hi > lo && std::fabs (g[hi]) < 1e-12) --hi;
            double asym = 0.0;
            for (std::size_t i = lo, j = hi; i < j; ++i, --j) asym = std::max (asym, std::fabs (g[i] - g[j]));
            worstSym = std::max (worstSym, asym);
            if (asym > 1e-6) allSymmetric = false;
            if (F == 4 && (fs == 44100.0 || fs == 48000.0 || fs == 88200.0 || fs == 96000.0))
            {
                int mac = 2 * core::firPadLen (d.firstTapsPerPhase) + core::firPadLen (d.taps[0]);
                for (int k = 1; k < d.stages; ++k) mac += (1 << k) * 2 * core::firPadLen (d.taps[k] / 2 + 1);   // up + down
                std::printf ("       %6.0f Hz 4x: t %3d, latency %3d, MAC %3d, strict %.2f dB, edge %.4f dB\n",
                             fs, d.firstTapsPerPhase, d.latencySamples, mac, std::max (ru.strictDb, rd.strictDb),
                             std::max (ru.edgeDevDb, rd.edgeDevDb));
            }
        }
    std::printf ("       over %zu rates x 4 factors: worst strict %.2f dB, worst one-pass deviation to the edge %.5f dB, "
                 "worst asymmetry %.1e\n", sizeof rates / sizeof rates[0], worstStrict, worstEdge, worstSym);
    test::ok (allPrepared, "every rate from 8 kHz to 768 kHz and every factor 2..16 prepares");
    test::ok (worstStrict <= -90.0, "STRICT: nothing at or above fs/2 comes through either path above -90 dB ("
                                    + std::to_string (worstStrict) + ")");
    test::ok (worstStrict > -95.0, "and the reading is not vacuous — the design sits near its bar, not at the float floor");
    test::ok (worstEdge <= 0.005, "FLAT: one pass stays within 0.005 dB up to the band edge (" + std::to_string (worstEdge) + ")");
    test::ok (allSymmetric, "LINEAR PHASE: every composite interpolator is a palindrome");

    // The table the header prints, pinned — a change of rule has to edit it on purpose.
    struct Row { double fs; int t, lat; };
    const Row rows[] = { { 44100.0, 125, 131 }, { 48000.0, 70, 76 }, { 88200.0, 22, 28 }, { 96000.0, 21, 27 },
                         { 32000.0, 125, 131 }, { 192000.0, 16, 22 } };
    for (const auto& r : rows)
    {
        CascadeOversampler::Design d;
        const bool ok = CascadeOversampler::designFor (r.fs, 4, d);
        test::ok (ok && d.firstTapsPerPhase == r.t && d.latencySamples == r.lat,
                  std::to_string ((int) r.fs) + " Hz 4x: t " + std::to_string (d.firstTapsPerPhase) + ", latency "
                  + std::to_string (d.latencySamples) + " (pinned " + std::to_string (r.t) + ", " + std::to_string (r.lat) + ")");
    }
    CascadeOversampler::Design d20, d32;
    (void) CascadeOversampler::designFor (44100.0, 4, d20);
    (void) CascadeOversampler::designFor (32000.0, 4, d32);
    test::approx (d20.bandEdgeHz, 20000.0, 1e-9, "the band edge is 20 kHz at 44.1 kHz");
    test::approx (d32.bandEdgeHz, 32000.0 * 20000.0 / 44100.0, 1e-6, "and 20/44.1 of the rate below it (32 kHz: 14.51 kHz)");

    // The guard the whole class exists for, in the class's own terms: what the FIXED 0.45 fs cut does at
    // 44.1 kHz, next to what this does. Same probe, both classes.
    {
        auto roundTrip = [] (auto&& up, auto&& down, int F, int lat, double hz)
        {
            const double fs = 44100.0;
            const int w = 4410, settle = 4096, n = settle + w + lat;
            std::vector<float> x ((std::size_t) n), u ((std::size_t) n * (std::size_t) F), y ((std::size_t) n);
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * hz / fs * (double) i);
            up (x.data(), n, u.data());
            down (u.data(), n, y.data());
            double re = 0.0, im = 0.0;
            for (int i = 0; i < w; ++i)
            {
                const double ph = 2.0 * kPi * hz / fs * (double) (settle + i);
                re += (double) y[(std::size_t) (settle + i)] * std::cos (ph);
                im += (double) y[(std::size_t) (settle + i)] * std::sin (ph);
            }
            return 20.0 * std::log10 (2.0 * std::hypot (re, im) / (double) w);
        };
        CascadeOversampler co; (void) co.prepare (44100.0, 4, 1);
        oversampling::PolyphaseOversampler po; (void) po.prepare (4, 1);
        auto cu = [&] (const float* x, int n, float* u) { const float* i[1] { x }; float* o[1] { u }; co.upsample (i, 1, n, o); };
        auto cd = [&] (const float* u, int n, float* y) { const float* i[1] { u }; float* o[1] { y }; co.downsample (i, 1, n, o); };
        auto pu = [&] (const float* x, int n, float* u) { const float* i[1] { x }; float* o[1] { u }; po.upsample (i, 1, n, o); };
        auto pd = [&] (const float* u, int n, float* y) { const float* i[1] { u }; float* o[1] { y }; po.downsample (i, 1, n, o); };
        const double c19 = roundTrip (cu, cd, 4, co.latencySamples(), 19000.0), c20 = roundTrip (cu, cd, 4, co.latencySamples(), 20000.0);
        const double p20 = roundTrip (pu, pd, 4, po.latencySamples(), 20000.0);
        std::printf ("       44.1 kHz round trip at 20 kHz: cascade %.4f dB, PolyphaseOversampler %.2f dB\n", c20, p20);
        test::ok (std::fabs (c19) < 0.01 && std::fabs (c20) < 0.01, "44.1 kHz: the cascade passes 19 and 20 kHz flat over a round trip ("
                                                                    + std::to_string (c19) + ", " + std::to_string (c20) + ")");
        test::ok (p20 < -15.0, "where the fixed-cutoff design loses 15.5 dB (the reason the class exists)");
    }
}

static void runLatencyTests()
{
    test::group ("CascadeOversampler: the round trip is a delayed identity at the reported, integer latency");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (int F : { 2, 4, 8, 16, 32, 64 })
        {
            CascadeOversampler co;
            if (! co.prepare (fs, F, 1)) { test::ok (false, "prepare " + std::to_string ((int) fs) + " x" + std::to_string (F)); continue; }
            const int L = co.latencySamples(), n = 2048, at = 700;
            std::vector<float> x ((std::size_t) n, 0.0f), u ((std::size_t) n * (std::size_t) F), y ((std::size_t) n);
            x[(std::size_t) at] = 1.0f;
            const float* i1[1] { x.data() }; float* o1[1] { u.data() };
            co.upsample (i1, 1, n, o1);
            const float* i2[1] { u.data() }; float* o2[1] { y.data() };
            co.downsample (i2, 1, n, o2);
            int pk = 0;
            for (int i = 0; i < n; ++i) if (std::fabs (y[(std::size_t) i]) > std::fabs (y[(std::size_t) pk])) pk = i;
            double asym = 0.0;
            for (int k = 1; k < 600; ++k) asym = std::max (asym, (double) std::fabs (y[(std::size_t) (at + L + k)] - y[(std::size_t) (at + L - k)]));
            const std::string tag = std::to_string ((int) fs) + " Hz x" + std::to_string (F);
            test::ok (pk == at + L, tag + ": the impulse comes back at the reported latency " + std::to_string (L)
                                    + " (peak at +" + std::to_string (pk - at) + ")");
            test::ok (asym < 1e-6, tag + ": and symmetric about it (" + std::to_string (asym) + ")");

            // A band-limited signal, not only an impulse: 1 kHz comes back as itself, delayed by L.
            std::vector<float> s ((std::size_t) n), su ((std::size_t) n * (std::size_t) F), sy ((std::size_t) n);
            for (int i = 0; i < n; ++i) s[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 1000.0 / fs * (double) i));
            co.reset();
            const float* i3[1] { s.data() }; float* o3[1] { su.data() };
            co.upsample (i3, 1, n, o3);
            const float* i4[1] { su.data() }; float* o4[1] { sy.data() };
            co.downsample (i4, 1, n, o4);
            double err = 0.0;
            for (int i = 400; i < n - L - 10; ++i) err = std::max (err, (double) std::fabs (sy[(std::size_t) (i + L)] - s[(std::size_t) i]));
            test::ok (err < 1e-3, tag + ": 1 kHz round-trips as a delayed identity (max error " + std::to_string (err) + ")");
        }
    CascadeOversampler fresh;
    test::ok (fresh.latencySamples() == 0 && fresh.factor() == 0, "unprepared: latency 0 and factor 0");
    test::ok (CascadeOversampler::latencyFor (44100.0, 3) == 0, "latencyFor answers 0 where prepare refuses");
}

static void runRefusalTests()
{
    test::group ("CascadeOversampler: refusals (law 11b) — whole, and touching nothing");
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    struct Case { double fs; int F, ch; const char* why; };
    const Case bad[] = {
        { nan, 4, 1, "NaN rate" }, { inf, 4, 1, "infinite rate" }, { -44100.0, 4, 1, "negative rate" }, { 0.0, 4, 1, "zero rate" },
        { 999.0, 4, 1, "rate below 1 kHz" }, { 3.0e6 + 1.0, 4, 1, "rate above 3 MHz" },
        { 44100.0, 0, 1, "factor 0" }, { 44100.0, 1, 1, "factor 1" }, { 44100.0, 3, 1, "factor 3 (not a power of two)" },
        { 44100.0, 6, 1, "factor 6" }, { 44100.0, 128, 1, "factor 128" }, { 44100.0, -4, 1, "factor -4" },
        { 44100.0, 4, 0, "no channels" }, { 44100.0, 4, core::kMaxChannels + 1, "channels past kMaxChannels (NOT clamped)" },
        { 44100.0, 4, -1, "negative channels" } };
    for (const auto& b : bad)
    {
        CascadeOversampler co;
        test::ok (co.prepare (48000.0, 2, 2), "setup");
        const int latBefore = co.latencySamples();
        CascadeOversampler::Storage st; st.coeffs = 12345;
        const bool budget = CascadeOversampler::storageFor (b.fs, b.F, b.ch, st);
        CascadeOversampler::Design d; d.latencySamples = 777;
        const bool design = CascadeOversampler::designFor (b.fs, b.F, d);
        const bool prep = co.prepare (b.fs, b.F, b.ch);
        test::ok (! prep && ! budget && st.coeffs == 12345, std::string ("refused: ") + b.why + " (prepare and storageFor agree, out untouched)");
        test::ok (co.latencySamples() == latBefore && co.factor() == 2, std::string ("and the refused call left the previous preparation standing: ") + b.why);
        if (b.ch >= 1 && b.ch <= core::kMaxChannels) test::ok (! design && d.latencySamples == 777, std::string ("designFor refuses it too: ") + b.why);
    }
    test::ok (CascadeOversampler {}.prepare (1000.0, 64, core::kMaxChannels), "the edges themselves are accepted (1 kHz, 64x, kMaxChannels)");
    test::ok (CascadeOversampler {}.prepare (3.0e6, 2, 1), "and 3 MHz");

    // Unprepared calls are no-ops, not bad indices.
    CascadeOversampler un;
    std::vector<float> a (64, 0.25f), b (256, 7.0f);
    const float* ia[1] { a.data() }; float* ob[1] { b.data() };
    un.upsample (ia, 1, 64, ob);
    const float* ib[1] { b.data() }; float* oa[1] { a.data() };
    un.downsample (ib, 1, 32, oa);
    un.reset(); un.resetChannel (0); un.resetChannel (-1);
    test::ok (b[0] == 7.0f && a[0] == 0.25f, "unprepared up/down/reset/resetChannel touch nothing");
}

static void runStorageAndRtTests()
{
    namespace alloc = test::alloc;
    test::group ("CascadeOversampler: what prepare() asks the heap for, and nothing after it");
    {
        // 44.1 kHz, 4x, stereo, derived by hand: stage 1 t = 125 -> tppPad 128 (two phases) + nPad 252
        // coefficients, rings 2*128 + 2*252 per channel; stage 2 halfband 27 -> 14 non-centre taps, pad 16,
        // coefficients 16 + 16, rings 2*16 + 2*16 + centre 2*7 per channel; scratch 2 * 64 * 2; cursors 2*2*3.
        CascadeOversampler::Storage st;
        test::ok (CascadeOversampler::storageFor (44100.0, 4, 2, st), "storageFor 44.1 kHz 4x stereo");
        test::ok (st.coeffs == 2 * 128 + 252 + 16 + 16, "coefficients: " + std::to_string (st.coeffs) + " (pinned 540)");
        test::ok (st.rings == 2 * (2 * 128 + 2 * 252 + 2 * 16 + 2 * 16 + 2 * 7), "rings: " + std::to_string (st.rings) + " (pinned 1676)");
        test::ok (st.scratch == 2 * 64 * 2 && st.cursors == 2 * 2 * 3, "scratch 256, cursors 12");
        const long long c0 = alloc::count.load(), b0 = alloc::bytes.load();
        CascadeOversampler co;
        const bool ok = co.prepare (44100.0, 4, 2);
        const long long dc = alloc::count.load() - c0, db = alloc::bytes.load() - b0;
        test::ok (ok, "prepare");
        test::ok (db == (long long) st.bytes(), "prepare() asked for exactly bytes() — " + std::to_string (db) + " of "
                                                + std::to_string (st.bytes()) + " — so no temporary and nothing the budget misses");
        test::ok (dc == 4, "in four allocations (" + std::to_string (dc) + ")");

        CascadeOversampler::Storage wide;
        (void) CascadeOversampler::storageFor (44100.0, 8, 2, wide);
        test::ok (st.fitsWithin (wide) && ! wide.fitsWithin (st), "fitsWithin orders a narrower topology below a wider one");
    }
    test::group ("CascadeOversampler: upsample/downsample allocate nothing, at any n");
    {
        CascadeOversampler co;
        (void) co.prepare (48000.0, 8, 2);
        std::vector<float> l (1000, 0.1f), r (1000, -0.1f), ul (8000), ur (8000), dl (1000), dr (1000);
        const float* in[2] { l.data(), r.data() }; float* up[2] { ul.data(), ur.data() };
        const float* upc[2] { ul.data(), ur.data() }; float* dn[2] { dl.data(), dr.data() };
        const long long c0 = alloc::count.load();
        for (int n : { 1, 63, 64, 65, 1000 })
        {
            co.upsample (in, 2, n, up);
            co.downsample (upc, 2, n, dn);
        }
        test::okNoAlloc (alloc::count.load() == c0, "no allocation across n = 1, 63, 64, 65, 1000");
    }
}

static void runStreamingTests()
{
    test::group ("CascadeOversampler: block splitting and channel count are bit-invisible; resetChannel is per channel");
    const int F = 8, n = 777;
    std::vector<float> x ((std::size_t) n);
    for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (0.4 * std::sin (0.37 * i) + 0.2 * std::sin (2.1 * i + 0.3));
    auto run = [&] (const std::vector<int>& cuts, std::vector<float>& up, std::vector<float>& dn)
    {
        CascadeOversampler co;
        (void) co.prepare (44100.0, F, 1);
        up.assign ((std::size_t) n * F, 0.0f);
        dn.assign ((std::size_t) n, 0.0f);
        int off = 0;
        for (std::size_t k = 0; off < n; ++k)
        {
            const int m = std::min (cuts[k % cuts.size()], n - off);
            const float* i1[1] { x.data() + off }; float* o1[1] { up.data() + (std::ptrdiff_t) off * F };
            co.upsample (i1, 1, m, o1);
            const float* i2[1] { up.data() + (std::ptrdiff_t) off * F }; float* o2[1] { dn.data() + off };
            co.downsample (i2, 1, m, o2);
            off += m;
        }
    };
    std::vector<float> u1, d1, u2, d2, u3, d3;
    run ({ n }, u1, d1);
    run ({ 1, 5, 64, 3, 129, 17 }, u2, d2);
    run ({ 63 }, u3, d3);
    test::ok (u1 == u2 && d1 == d2 && u1 == u3 && d1 == d3, "whole, ragged and 63-sample blocks give the same bits, both ways");

    // Two channels: channel 1 is bit-identical to a mono run of the same signal; resetChannel(0) leaves it so.
    CascadeOversampler mono, st;
    (void) mono.prepare (44100.0, 4, 1);
    (void) st.prepare (44100.0, 4, 2);
    std::vector<float> a (200), b (200);
    for (int i = 0; i < 200; ++i) { a[(std::size_t) i] = (float) std::sin (0.1 * i); b[(std::size_t) i] = (float) std::cos (0.23 * i); }
    std::vector<float> ma (800), sa (800), sb (800);
    const float* im[1] { b.data() }; float* om[1] { ma.data() };
    const float* is[2] { a.data(), b.data() }; float* os[2] { sa.data(), sb.data() };
    mono.upsample (im, 1, 100, om);
    st.upsample (is, 2, 100, os);
    st.resetChannel (0);
    const float* im2[1] { b.data() + 100 }; float* om2[1] { ma.data() + 400 };
    const float* is2[2] { a.data() + 100, b.data() + 100 }; float* os2[2] { sa.data() + 400, sb.data() + 400 };
    mono.upsample (im2, 1, 100, om2);
    st.upsample (is2, 2, 100, os2);
    test::ok (ma == sb, "the untouched channel is bit-identical to a mono stream that never saw a reset");
    CascadeOversampler fresh; (void) fresh.prepare (44100.0, 4, 1);
    std::vector<float> fa (400);
    const float* ifr[1] { a.data() + 100 }; float* ofr[1] { fa.data() };
    fresh.upsample (ifr, 1, 100, ofr);
    test::ok (std::equal (fa.begin(), fa.end(), sa.begin() + 400), "and the reset one equals a fresh channel");
    test::ok (! std::equal (fa.begin(), fa.end(), sa.begin()), "(the comparison is not vacuous)");
}

static void runBitPins()
{
    // Law 2: the composite taps are the same float words on every row. The hashes were taken on arm64
    // Apple clang and must hold on gcc/x86-64, MSVC/x64, Apple clang/x86-64, gcc/arm64 and wasm; a row where
    // the designed doubles or the narrowing differ fails here by name rather than drifting.
    test::group ("CascadeOversampler: the designed filters are the same bits on every row");
    struct Pin { double fs; int F; std::uint64_t hash; };
    const Pin pins[] = { { 44100.0, 4, 0x06a22c34de9944d0ull }, { 48000.0, 8, 0x9cf934f564435a4full },
                         { 96000.0, 2, 0x0893f7c6b2fa244dull }, { 8000.0, 16, 0x641a714e44beaeb5ull } };
    for (const auto& p : pins)
    {
        CascadeOversampler co;
        (void) co.prepare (p.fs, p.F, 1);
        const int n = 600;
        std::vector<float> x ((std::size_t) n, 0.0f), y ((std::size_t) n * (std::size_t) p.F);
        x[0] = 1.0f;
        const float* in[1] { x.data() }; float* out[1] { y.data() };
        co.upsample (in, 1, n, out);
        std::vector<float> dn ((std::size_t) n);
        const float* in2[1] { y.data() }; float* out2[1] { dn.data() };
        co.downsample (in2, 1, n, out2);
        y.insert (y.end(), dn.begin(), dn.end());
        const std::uint64_t h = fnv (y);
        char buf[32]; std::snprintf (buf, sizeof buf, "%016llx", (unsigned long long) h);
        test::ok (h == p.hash, std::to_string ((int) p.fs) + " Hz x" + std::to_string (p.F) + ": impulse up+down hash " + buf);
    }
}

static void runTopologyTests()
{
    test::group ("Oversampler: Kaiser IS PolyphaseOversampler (same refusals, same bits); Cascade IS CascadeOversampler");
    std::vector<float> x (500);
    for (int i = 0; i < 500; ++i) x[(std::size_t) i] = (float) (0.6 * std::sin (0.21 * i) + 0.3 * std::sin (1.9 * i));
    for (int tpp : { 12, 32, 64 })
        for (int F : { 2, 3, 4, 8 })
        {
            oversampling::Oversampler ov;
            oversampling::PolyphaseOversampler po;
            const bool a = ov.prepare (oversampling::Topology::Kaiser, 44100.0, F, 2, tpp);
            const bool b = po.prepare (F, 2, tpp);
            std::vector<float> u1 (500u * (unsigned) F), u2 (500u * (unsigned) F), d1 (500), d2 (500);
            const float* in[1] { x.data() };
            float* o1[1] { u1.data() }; float* o2[1] { u2.data() };
            ov.upsample (in, 1, 500, o1); po.upsample (in, 1, 500, o2);
            const float* c1[1] { u1.data() }; const float* c2[1] { u2.data() };
            float* e1[1] { d1.data() }; float* e2[1] { d2.data() };
            ov.downsample (c1, 1, 500, e1); po.downsample (c2, 1, 500, e2);
            test::ok (a && b && u1 == u2 && d1 == d2 && ov.latencySamples() == po.latencySamples()
                      && ov.latencySamples() == oversampling::Oversampler::latencyFor (oversampling::Topology::Kaiser, 44100.0, F, tpp),
                      "Kaiser x" + std::to_string (F) + " tpp " + std::to_string (tpp) + ": bit-identical, same latency");
        }
    // The Kaiser refusal set is PolyphaseOversampler's, rate ignored — including its channel CLAMP (P55).
    oversampling::Oversampler::Storage s1;
    oversampling::PolyphaseOversampler::Storage s2;
    test::ok (oversampling::Oversampler::storageFor (oversampling::Topology::Kaiser, std::numeric_limits<double>::quiet_NaN(), 4, 40, 64, s1)
              == oversampling::PolyphaseOversampler::storageFor (4, 40, 64, s2), "Kaiser ignores the rate and clamps channels, as before");
    test::ok (! oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 4, 1, 2000, s1),
              "Cascade still range-checks tapsPerPhase, so the refusal set does not shrink with the topology");
    test::ok (! oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 3, 1, 64, s1),
              "Cascade refuses a factor that is not a power of two");

    oversampling::Oversampler ov;
    CascadeOversampler co;
    (void) ov.prepare (oversampling::Topology::Cascade, 48000.0, 4, 1, 64);
    (void) co.prepare (48000.0, 4, 1);
    std::vector<float> u1 (2000), u2 (2000);
    const float* in[1] { x.data() }; float* o1[1] { u1.data() }; float* o2[1] { u2.data() };
    ov.upsample (in, 1, 500, o1); co.upsample (in, 1, 500, o2);
    test::ok (u1 == u2 && ov.latencySamples() == co.latencySamples() && ov.topology() == oversampling::Topology::Cascade,
              "Cascade through the switch is bit-identical to the class");

    // Switching topology on a live object: a refused preparation keeps the old one; a successful one runs the new.
    const int latC = ov.latencySamples();
    test::ok (! ov.prepare (oversampling::Topology::Kaiser, 48000.0, 4, 1, 2), "a refused Kaiser preparation...");
    test::ok (ov.topology() == oversampling::Topology::Cascade && ov.latencySamples() == latC, "...leaves the cascade in place");
    test::ok (ov.prepare (oversampling::Topology::Kaiser, 48000.0, 4, 1, 64) && ov.latencySamples() == 63, "and a successful one switches");
}

int main()
{
    std::printf ("felitronics::oversampling cascade tests\n");
    runRuleTests();
    runLatencyTests();
    runRefusalTests();
    runStorageAndRtTests();
    runStreamingTests();
    runBitPins();
    runTopologyTests();
    return test::report();
}
