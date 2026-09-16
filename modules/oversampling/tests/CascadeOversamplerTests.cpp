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
#include <type_traits>
#include <utility>
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
        const int F = co.factor(), n = F >= 32 ? 512 : 2048;      // the composite spans < 200 base samples
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
        const int F = co.factor(), n = F >= 32 ? 512 : 1024;
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
    bool allPrepared = true, allSymmetric = true, upLegOk = true, macOk = true;
    for (double fs : rates)
        for (int F : { 2, 4, 8, 16, 32, 64 })
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
            // ...and its centre is where the design says the UP leg ends (the offset a consumer would crop
            // an oversampled trace by — the one number here no other check reads).
            if ((long) (lo + hi) != (long) d.upLegTwice) upLegOk = false;
            if (F == 4 && (fs == 44100.0 || fs == 48000.0 || fs == 88200.0 || fs == 96000.0))
            {
                int mac = 2 * core::firPadLen (d.firstTapsPerPhase) + core::firPadLen (d.taps[0]);
                for (int k = 1; k < d.stages; ++k) mac += (1 << k) * 2 * core::firPadLen (d.taps[k] / 2 + 1);   // up + down
                std::printf ("       %6.0f Hz 4x: t %3d, latency %3d, MAC %3d, strict %.2f dB, edge %.4f dB\n",
                             fs, d.firstTapsPerPhase, d.latencySamples, mac, std::max (ru.strictDb, rd.strictDb),
                             std::max (ru.edgeDevDb, rd.edgeDevDb));
                const int wantMac = fs == 44100.0 ? 572 : (fs == 48000.0 ? 348 : 156);
                if (mac != wantMac) macOk = false;
            }
        }
    std::printf ("       over %zu rates x 6 factors: worst strict %.2f dB, worst one-pass deviation to the edge %.5f dB, "
                 "worst asymmetry %.1e\n", sizeof rates / sizeof rates[0], worstStrict, worstEdge, worstSym);
    test::ok (allPrepared, "every rate from 8 kHz to 768 kHz and every factor 2..64 prepares");
    test::ok (worstStrict <= -90.0, "STRICT: nothing at or above fs/2 comes through either path above -90 dB ("
                                    + std::to_string (worstStrict) + ")");
    test::ok (worstStrict > -95.0, "and the reading is not vacuous — the design sits near its bar, not at the float floor");
    test::ok (worstEdge <= 0.005, "FLAT: one pass stays within 0.005 dB up to the band edge (" + std::to_string (worstEdge) + ")");
    test::ok (allSymmetric, "LINEAR PHASE: every composite interpolator is a palindrome");
    test::ok (upLegOk, "and Design::upLegTwice is twice its centre, in top-rate samples, at every rate and factor");
    // The multiply count the header prints: firDot lengths per base sample, up + down, padding included. The
    // fixed-cutoff stage at 4x/64 is 4 * 64 up + 256 down.
    constexpr int kaiserMac = 4 * 64 + core::firPadLen (4 * 64);
    test::ok (macOk && kaiserMac == 512, "4x multiply count: 572 / 348 / 156 / 156 at 44.1 / 48 / 88.2 / 96 kHz, against 512 for the Kaiser stage");

    // EVERY taps count the rule can produce, not only the ones fourteen rates happen to reach (they reach
    // eight). Stage 1 alone (F = 2) is what the rule sizes; the first rate that yields each t is the probe.
    {
        std::vector<std::pair<int, double>> firstRate;
        for (double fs = CascadeOversampler::kMinSampleRate; fs <= CascadeOversampler::kMaxSampleRate; fs *= 1.0005)
        {
            CascadeOversampler::Design d;
            if (! CascadeOversampler::designFor (fs, 2, d)) continue;
            bool seen = false;
            for (const auto& fr : firstRate) seen = seen || fr.first == d.firstTapsPerPhase;
            if (! seen) firstRate.emplace_back (d.firstTapsPerPhase, fs);
        }
        double wS = -1e9, wE = 0.0; int tS = 0, tMin = 1 << 30, tMax = 0;
        for (const auto& [t, fs] : firstRate)
        {
            CascadeOversampler co;
            (void) co.prepare (fs, 2, 1);
            const Reading r = read (upComposite (co), fs, 2, co.design().bandEdgeHz);
            if (r.strictDb > wS) { wS = r.strictDb; tS = t; }
            wE = std::max (wE, r.edgeDevDb);
            tMin = std::min (tMin, t); tMax = std::max (tMax, t);
        }
        std::printf ("       every reachable taps count (%zu, %d..%d): worst strict %.2f dB (t %d), worst edge %.5f dB\n",
                     firstRate.size(), tMin, tMax, wS, tS, wE);
        test::ok (firstRate.size() == 110 && tMin == 16 && tMax == 125, "the rule produces 110 taps counts, 16..125, over 8 kHz..3 MHz");
        test::ok (wS <= -90.0 && wE <= 0.005, "and every one of them is strict and flat to its band edge ("
                                              + std::to_string (wS) + " dB, " + std::to_string (wE) + " dB)");
    }

    // The table the header prints, pinned — a change of rule has to edit it on purpose.
    struct Row { double fs; int t, lat; };
    // The last two rows are the FLOOR (kMinFirstTaps): the rule alone would build 14 and 13 taps per phase
    // there, still under -90 but by 0.7 dB instead of 2.3 — so the floor is a margin, and it is pinned
    // rather than left to a strictness bar it happens not to cross (a mutation stand removed it unseen).
    const Row rows[] = { { 44100.0, 125, 131 }, { 48000.0, 70, 76 }, { 88200.0, 22, 28 }, { 96000.0, 21, 27 },
                         { 32000.0, 125, 131 }, { 192000.0, 16, 22 }, { 384000.0, 16, 22 }, { 768000.0, 16, 22 } };
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
        { std::nextafter (core::kMinSampleRate, 0.0), 4, 1, "rate just below the core's 8 kHz floor" }, { 1000.0, 4, 1, "1 kHz" },
        { 3.0e6 + 1.0, 4, 1, "rate above 3 MHz" },
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
        // ...and its bits: the object that saw the refusal and one that never did produce the same stream.
        CascadeOversampler twin;
        (void) twin.prepare (48000.0, 2, 2);
        std::vector<float> xs (300), a1 (600), a2 (600);
        for (int i = 0; i < 300; ++i) xs[(std::size_t) i] = (float) std::sin (0.17 * i);
        const float* ix1[1] { xs.data() }; float* oo1[1] { a1.data() }; float* oo2[1] { a2.data() };
        co.upsample (ix1, 1, 300, oo1);
        twin.upsample (ix1, 1, 300, oo2);
        test::ok (a1 == a2, std::string ("and the refused call left the coefficients and rings untouched: ") + b.why);
        if (b.ch >= 1 && b.ch <= core::kMaxChannels) test::ok (! design && d.latencySamples == 777, std::string ("designFor refuses it too: ") + b.why);
    }
    test::ok (CascadeOversampler::kMinSampleRate == core::kMinSampleRate, "the floor equals the core's one floor (P51)");
    test::ok (CascadeOversampler {}.prepare (core::kMinSampleRate, 64, core::kMaxChannels), "the edges themselves are accepted (8 kHz, 64x, kMaxChannels)");
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
        // Read the counter BEFORE the call: the message argument is a std::string built in the same argument
        // list, it allocates, and gcc evaluates it first (a refusal check below was red on gcc for exactly this).
        const bool none = alloc::count.load() == c0;
        test::okNoAlloc (none, "no allocation across n = 1, 63, 64, 65, 1000");
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

    // THE DOWN PATH per channel, with DIFFERENT content: every channel of a stereo decimator equals a mono one
    // fed the same samples (a decimator that read channel 0 for every channel passed everything above).
    {
        const int Fd = 8, nd = 300;
        CascadeOversampler st2, m0, m1;
        (void) st2.prepare (48000.0, Fd, 2); (void) m0.prepare (48000.0, Fd, 1); (void) m1.prepare (48000.0, Fd, 1);
        std::vector<float> u0 ((std::size_t) (nd * Fd)), u1 ((std::size_t) (nd * Fd));
        for (int i = 0; i < nd * Fd; ++i) { u0[(std::size_t) i] = (float) std::sin (0.013 * i); u1[(std::size_t) i] = (float) std::cos (0.029 * i + 1.0); }
        std::vector<float> s0 ((std::size_t) nd), s1 ((std::size_t) nd), r0 ((std::size_t) nd), r1 ((std::size_t) nd);
        const float* si[2] { u0.data(), u1.data() }; float* so[2] { s0.data(), s1.data() };
        st2.downsample (si, 2, nd, so);
        const float* a0[1] { u0.data() }; float* b0[1] { r0.data() };
        const float* a1[1] { u1.data() }; float* b1[1] { r1.data() };
        m0.downsample (a0, 1, nd, b0); m1.downsample (a1, 1, nd, b1);
        test::ok (s0 == r0 && s1 == r1 && s0 != s1, "a stereo decimator is two mono ones, channel for channel, on different content");
    }

    // Channels past the prepared count are not touched, both ways (the consumers refuse such calls first;
    // this is the class's own promise, and without the clamp it indexes rings that do not exist).
    {
        CascadeOversampler one;
        (void) one.prepare (44100.0, 4, 1);
        std::vector<float> c0 (64, 0.25f), c1 (64, 0.5f), u0 (256, 0.0f), u1 (256, 9.0f), d0 (64, 0.0f), d1 (64, 9.0f);
        const float* in2[2] { c0.data(), c1.data() }; float* up2[2] { u0.data(), u1.data() };
        one.upsample (in2, 2, 64, up2);
        const float* uin[2] { u0.data(), u0.data() }; float* dn2[2] { d0.data(), d1.data() };
        one.downsample (uin, 2, 64, dn2);
        bool untouched = true;
        for (float v : u1) untouched = untouched && v == 9.0f;
        for (float v : d1) untouched = untouched && v == 9.0f;
        test::ok (untouched, "a call wider than the preparation leaves the extra channel's output untouched, up and down");
    }

    // resetChannel outside [0, channels) is a no-op on a PREPARED object — the other channels continue bit for
    // bit. (An out-of-range index would write outside the rings; the sanitizer row sees that, this sees the rest.)
    {
        CascadeOversampler p2, q2;
        (void) p2.prepare (44100.0, 4, 2); (void) q2.prepare (44100.0, 4, 2);
        std::vector<float> in (200); for (int i = 0; i < 200; ++i) in[(std::size_t) i] = (float) std::sin (0.3 * i);
        std::vector<float> pa (800), pb (800), qa (800), qb (800);
        const float* ii[2] { in.data(), in.data() };
        float* po[2] { pa.data(), pb.data() }; float* qo[2] { qa.data(), qb.data() };
        p2.upsample (ii, 2, 100, po); q2.upsample (ii, 2, 100, qo);
        p2.resetChannel (-1); p2.resetChannel (2); p2.resetChannel (1 << 20);
        const float* ii2[2] { in.data() + 100, in.data() + 100 };
        float* po2[2] { pa.data() + 400, pb.data() + 400 }; float* qo2[2] { qa.data() + 400, qb.data() + 400 };
        p2.upsample (ii2, 2, 100, po2); q2.upsample (ii2, 2, 100, qo2);
        test::ok (pa == qa && pb == qb, "resetChannel(-1), (channels) and (huge) change nothing on a prepared object");
    }
}

static void runBitPins()
{
    // Law 2: the composite taps are the same float words on every row. The hashes were taken on arm64
    // Apple clang and must hold on gcc/x86-64, MSVC/x64, Apple clang/x86-64, gcc/arm64 and wasm; a row where
    // the designed doubles or the narrowing differ fails here by name rather than drifting.
    test::group ("CascadeOversampler: the designed filters are the same bits on every row");
    struct Pin { double fs; int F; std::uint64_t hash; };
    // 32x and 64x reach the last two halfband stages, which nothing at 16x and below runs.
    const Pin pins[] = { { 44100.0, 4, 0x06a22c34de9944d0ull }, { 48000.0, 8, 0x9cf934f564435a4full },
                         { 96000.0, 2, 0x0893f7c6b2fa244dull }, { 8000.0, 16, 0x641a714e44beaeb5ull },
                         { 44100.0, 32, 0xf54ec22e379c6f60ull }, { 48000.0, 64, 0x487c11a907b74918ull } };
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
    test::ok (! oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 4, 1, 2000, s1)
              && ! oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 4, 1, 3, s1)
              && oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 4, 1, 4, s1),
              "Cascade still range-checks tapsPerPhase at BOTH ends (3 and 2000 refused, 4 accepted), so the refusal set does not shrink with the topology");
    test::ok (! oversampling::Oversampler::storageFor (oversampling::Topology::Cascade, 44100.0, 3, 1, 64, s1),
              "Cascade refuses a factor that is not a power of two");
    {
        using oversampling::Topology; using OS = oversampling::Oversampler;
        test::ok (OS::latencyFor (Topology::Kaiser, 44100.0, 4, 64) == 63 && OS::latencyFor (Topology::Kaiser, 44100.0, 1, 64) == 0
                  && OS::latencyFor (Topology::Kaiser, 44100.0, 4, 2000) == 0 && OS::latencyFor (Topology::Kaiser, 44100.0, 4, 3) == 0
                  && OS::latencyFor (Topology::Cascade, 44100.0, 4, 3) == 0 && OS::latencyFor (Topology::Cascade, 44100.0, 4, 2000) == 0
                  && OS::latencyFor (Topology::Cascade, 44100.0, 4, 64) == 131,
                  "latencyFor answers 0 wherever the preparation would be refused, under EITHER topology");
        OS sw; (void) sw.prepare (Topology::Cascade, 44100.0, 4, 1, 64);
        test::ok (! sw.prepare (Topology::Cascade, 44100.0, 4, 1, 2000) && ! sw.prepare (Topology::Cascade, 44100.0, 4, 1, 3)
                  && sw.latencySamples() == 131 && sw.topology() == Topology::Cascade,
                  "the switch's own prepare() refuses a bad tapsPerPhase under the cascade, and touches nothing");
    }

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

    // The default must not pay for the option: the cascade is held on the heap, only when chosen, and the
    // switch stays copyable (the stages that embed it are copied and moved by callers).
    std::printf ("       sizeof: Oversampler %zu, PolyphaseOversampler %zu, CascadeOversampler %zu (heap-held)\n",
                 sizeof (oversampling::Oversampler), sizeof (oversampling::PolyphaseOversampler), sizeof (CascadeOversampler));
    test::ok (sizeof (oversampling::Oversampler) <= sizeof (oversampling::PolyphaseOversampler) + 40,
              "the switch costs the Kaiser path at most 40 bytes of object size (" + std::to_string (sizeof (oversampling::Oversampler))
              + " against " + std::to_string (sizeof (oversampling::PolyphaseOversampler)) + ")");
    test::ok (std::is_copy_constructible_v<oversampling::Oversampler> && std::is_copy_assignable_v<oversampling::Oversampler>,
              "and it is copyable, as PolyphaseOversampler is");
    {
        namespace alloc = test::alloc;
        using oversampling::Topology;
        oversampling::Oversampler sw;
        oversampling::Oversampler::Storage bk, bc;
        (void) oversampling::Oversampler::storageFor (Topology::Kaiser, 44100.0, 4, 2, 64, bk);
        (void) oversampling::Oversampler::storageFor (Topology::Cascade, 44100.0, 4, 2, 64, bc);
        auto asked = [&] (Topology t) { const long long b0 = alloc::bytes.load(); (void) sw.prepare (t, 44100.0, 4, 2, 64); return alloc::bytes.load() - b0; };
        const long long k1 = asked (Topology::Kaiser);
        const long long c1 = asked (Topology::Cascade);
        const long long c2 = asked (Topology::Cascade);
        const long long k2 = asked (Topology::Kaiser);
        const long long c3 = asked (Topology::Cascade);
        test::ok (k1 == (long long) bk.bytes() && bk.heapObjects == 0, "Kaiser asks for exactly its budget, no heap object ("
                  + std::to_string (k1) + " B)");
        test::ok (c1 == (long long) bc.bytes() && bc.heapObjects == 1, "Cascade asks for exactly its budget, the heap object included ("
                  + std::to_string (c1) + " B)");
        test::ok (c2 == 0, "re-preparing the same cascade asks for nothing (" + std::to_string (c2) + " B)");
        // A switch RELEASES the topology it leaves: going back asks for the whole budget again, which it
        // would not if the old buffers (or the heap object) had been kept.
        test::ok (k2 == (long long) bk.bytes(), "switching back to Kaiser asks for its whole budget again — the Kaiser buffers were released ("
                  + std::to_string (k2) + " B)");
        test::ok (c3 == (long long) bc.bytes(), "and back to Cascade asks for the whole cascade budget, heap object included — it was released too ("
                  + std::to_string (c3) + " B)");

        // A REFUSED preparation touches nothing: it allocates nothing on a fresh switch, and a live Kaiser
        // switch that is refused a cascade still runs its Kaiser bits.
        oversampling::Oversampler fresh;
        const long long f0 = alloc::bytes.load();
        const bool refused = ! fresh.prepare (Topology::Cascade, 44100.0, 3, 1, 64);
        const long long refusedBytes = alloc::bytes.load() - f0;      // read before the message string exists
        test::ok (refused && refusedBytes == 0, "a refused cascade preparation allocates nothing (not even the heap object)");
        oversampling::Oversampler live, twin;
        (void) live.prepare (Topology::Kaiser, 44100.0, 4, 1, 64);
        (void) twin.prepare (Topology::Kaiser, 44100.0, 4, 1, 64);
        test::ok (! live.prepare (Topology::Cascade, 44100.0, 3, 1, 64) && live.topology() == Topology::Kaiser && live.latencySamples() == 63,
                  "a refused cascade leaves a live Kaiser switch in place");
        std::vector<float> xi (100, 0.0f), l1 (400), l2 (400);
        for (int i = 0; i < 100; ++i) xi[(std::size_t) i] = (float) std::sin (0.2 * i);
        const float* xip[1] { xi.data() }; float* lo1[1] { l1.data() }; float* lo2[1] { l2.data() };
        live.upsample (xip, 1, 100, lo1); twin.upsample (xip, 1, 100, lo2);
        test::ok (l1 == l2 && l1[50] != 0.0f, "...with its coefficients and rings intact (bit-identical to a twin)");

        // Storage ordering sees the cascade half and the heap object, not only the Kaiser half.
        test::ok (! bc.fitsWithin (bk) && ! bk.fitsWithin (bc), "a cascade budget does not fit a Kaiser one, nor the reverse");
        oversampling::Oversampler::Storage noHeap = bc; noHeap.heapObjects = 0;
        test::ok (noHeap.fitsWithin (bc) && ! bc.fitsWithin (noHeap), "and the heap object is part of the comparison");
        oversampling::Oversampler::Storage c4, c8;
        (void) oversampling::Oversampler::storageFor (Topology::Cascade, 44100.0, 4, 2, 64, c4);
        (void) oversampling::Oversampler::storageFor (Topology::Cascade, 44100.0, 8, 2, 64, c8);
        test::ok (c4.fitsWithin (c8) && ! c8.fitsWithin (c4), "and two CASCADE budgets are ordered by their cascade halves (4x fits 8x, not the reverse)");
        oversampling::Oversampler copy = sw;                     // sw holds the cascade here
        std::vector<float> ci (50, 0.0f), co1 (200), co2 (200);
        for (int i = 0; i < 50; ++i) ci[(std::size_t) i] = (float) std::cos (0.4 * i);
        const float* cip[1] { ci.data() }; float* cop1[1] { co1.data() }; float* cop2[1] { co2.data() };
        copy.upsample (cip, 1, 50, cop1); sw.upsample (cip, 1, 50, cop2);
        test::ok (copy.topology() == Topology::Cascade && copy.latencySamples() == sw.latencySamples() && co1 == co2,
                  "a copy of a cascade switch is a deep copy: same topology, same latency, same bits");

        // reset() and resetChannel() under the cascade reach the CASCADE: after either, the stream equals a
        // fresh switch's (a switch that routed them to its idle Kaiser member passed everything above).
        {
            oversampling::Oversampler dirty, clean;
            (void) dirty.prepare (Topology::Cascade, 48000.0, 4, 2, 64);
            (void) clean.prepare (Topology::Cascade, 48000.0, 4, 2, 64);
            std::vector<float> noise (200);
            for (int i = 0; i < 200; ++i) noise[(std::size_t) i] = (float) std::sin (1.7 * i) * 0.8f;
            std::vector<float> junk (800), o1 (400), o2 (400), o3 (400), o4 (400);
            const float* ni[2] { noise.data(), noise.data() }; float* jo[2] { junk.data(), junk.data() + 400 };
            dirty.upsample (ni, 2, 100, jo);
            dirty.reset();
            const float* ti[2] { noise.data() + 100, noise.data() + 100 };
            float* d2[2] { o1.data(), o2.data() }; float* c2o[2] { o3.data(), o4.data() };
            dirty.upsample (ti, 2, 100, d2); clean.upsample (ti, 2, 100, c2o);
            test::ok (o1 == o3 && o2 == o4, "reset() under the cascade clears the cascade");
            dirty.upsample (ni, 2, 100, jo);
            dirty.resetChannel (0);
            oversampling::Oversampler clean2;
            (void) clean2.prepare (Topology::Cascade, 48000.0, 4, 1, 64);
            float* d3[2] { o1.data(), o2.data() }; float* c3o[1] { o3.data() };
            dirty.upsample (ti, 2, 100, d3); clean2.upsample (ti, 1, 100, c3o);
            test::ok (o1 == o3, "and resetChannel(0) under the cascade clears channel 0 of the cascade");
        }

        // A MOVED-FROM cascade switch keeps its topology tag and loses its object. It must read as unprepared
        // and ignore calls, not dereference the empty vector (it did: SIGSEGV on latencySamples()).
        oversampling::Oversampler src;
        (void) src.prepare (Topology::Cascade, 48000.0, 4, 1, 64);
        oversampling::Oversampler dst = std::move (src);
        std::vector<float> xin (64, 0.5f), xup (256, 7.0f), xdn (64, 7.0f);
        const float* ii[1] { xin.data() }; float* uo[1] { xup.data() };
        const float* ui[1] { xup.data() }; float* dno[1] { xdn.data() };
        src.upsample (ii, 1, 64, uo); src.downsample (ui, 1, 64, dno); src.reset(); src.resetChannel (0);   // NOLINT: use after move is the test
        test::ok (src.latencySamples() == 0 && src.factor() == 0 && xup[0] == 7.0f && xdn[0] == 7.0f,
                  "a moved-from cascade switch reads as unprepared and its calls touch nothing");
        test::ok (dst.latencySamples() == 76 && dst.factor() == 4, "and the moved-to switch is the prepared cascade");
        oversampling::Oversampler again;
        again = std::move (dst);
        test::ok (again.latencySamples() == 76 && dst.latencySamples() == 0, "move assignment behaves the same way");
    }
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
