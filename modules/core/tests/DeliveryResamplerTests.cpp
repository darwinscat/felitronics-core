// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P57 — core::DeliveryResampler, the contract half: the pure planner over all 30 directed delivery
// pairs, law-11 call discipline clause by clause, RT-safety, the exact latency measured back out of the
// audio, output counts against an independent arithmetic oracle, and the 48 -> 44.1 -> 48 round-trip
// NULL. The BAR itself lives in DeliveryResamplerBarTests.cpp.
//
// Every oracle below is INDEPENDENT of the automaton it checks: counts come from floor((T*L - 1)/M) + 1
// composed through the plan, latency from the carrier phase, the flush from feeding zeros, the gap from
// feeding silence, the channels from mono twins. The first version of this suite compared the
// converter against itself — three slicings of one run, a formula recomputed from the plan's own
// fields — and a mutation stand and a crew round found eight broken implementations it passed.

#include <felitronics/core/DeliveryResampler.h>
#include <felitronics_test.h>

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <vector>

static std::atomic<long long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using felitronics::core::DeliveryResampler;
using felitronics::core::exactlyEqual;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;
using felitronics::test::okNoAlloc;

namespace
{
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

    DeliveryResampler::Params params (double in, double out)
    {
        DeliveryResampler::Params p;
        p.inRate = in; p.outRate = out;
        return p;
    }

    // Feed `x` block by block through a fresh converter, then drain. Everything emitted comes back,
    // including the leading `latencyOutputSamples()` of ramp-up.
    std::vector<float> convert (double inR, double outR, const std::vector<float>& x, int block,
                                double& latency, bool& accepted)
    {
        DeliveryResampler r;
        accepted = r.prepare (params (inR, outR), 1, block);
        if (! accepted) return {};
        latency = r.latencyOutputSamples();
        std::vector<float> y;
        std::vector<float> buf ((std::size_t) r.maxOutputFor (block));
        float* op[1] = { buf.data() };
        for (std::size_t i = 0; i < x.size(); i += (std::size_t) block)
        {
            const int n = (int) std::min ((std::size_t) block, x.size() - i);
            const float* ip[1] = { x.data() + i };
            int got = 0;
            if (! r.process (ip, 1, n, op, (int) buf.size(), got)) { accepted = false; return {}; }
            y.insert (y.end(), buf.begin(), buf.begin() + got);
        }
        std::vector<float> fb ((std::size_t) r.maxFlushOutput());
        float* fp[1] = { fb.data() };
        int got = 0;
        if (! r.flush (1, fp, (int) fb.size(), got)) { accepted = false; return {}; }
        y.insert (y.end(), fb.begin(), fb.begin() + got);
        return y;
    }

    double dbfs (double v) { return 20.0 * std::log10 (std::max (1e-300, v)); }

    bool sameBits (const std::vector<float>& a, const std::vector<float>& b)
    {
        return a.size() == b.size() && (a.empty() || std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0);
    }

    // The independent count oracle: after T inputs a zero-primed L:M stage has emitted
    // floor((T*L - 1)/M) + 1 outputs (none at T = 0). Composed through the plan's stages.
    long long cumulativeCount (const DeliveryResampler::Plan& pl, long long T)
    {
        for (int s = 0; s < pl.count; ++s)
            T = (T <= 0) ? 0 : ((T * pl.stage[s].L - 1) / pl.stage[s].M + 1);
        return T;
    }

    // Delay in output samples recovered from a steady tone's phase at f (one exact second analysed).
    double measuredDelay (const std::vector<float>& y, double fsOut, double f, int start)
    {
        const int len = (int) std::llround (fsOut);
        double re = 0.0, im = 0.0;
        for (int i = 0; i < len; ++i)
        {
            const double t = (double) (start + i) / fsOut;
            const double v = (double) y[(std::size_t) (start + i)];
            re += v * std::cos (2.0 * kPi * f * t);
            im -= v * std::sin (2.0 * kPi * f * t);
        }
        // y(t) = sin(2*pi*f*t - phi), phi = 2*pi*f*D/fsOut  ->  re = -(N/2) sin(phi), im = -(N/2) cos(phi)
        double phase = std::atan2 (-re, -im);
        if (phase < 0.0) phase += 2.0 * kPi;
        return phase / (2.0 * kPi * f) * fsOut;
    }
}

//==============================================================================
// 1. THE PLANNER
static void testPlan()
{
    group ("plan: every delivery pair routes, and the route obeys the cost law");
    int worstStages = 0;
    double worstMacs = 0.0;
    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;                                   // exact equality: these are literals
            const auto pl = DeliveryResampler::plan (params (a, b));
            ok (pl.ok, "plan exists");
            ok (! pl.identity, "a rate change is not the identity path");
            ok (pl.count >= 1 && pl.count <= DeliveryResampler::kMaxStages, "stage count within the search's bound");
            ok (pl.stage[0].inRate == (long long) a, "route starts at the input rate");
            ok (pl.stage[pl.count - 1].outRate == (long long) b, "route ends at the output rate");
            for (int s = 0; s + 1 < pl.count; ++s)
                ok (pl.stage[s].outRate == pl.stage[s + 1].inRate, "the route is connected");
            for (int s = 0; s < pl.count; ++s)
            {
                const auto& sp = pl.stage[s];
                ok (sp.tapsPerPhase % 4 == 0, "taps per phase are a multiple of four");
                ok (sp.tapsPerPhase >= 2 * sp.halfLen + 1, "the padded row holds the whole phase");
                ok (sp.passbandHz < sp.stopbandHz && sp.stopbandHz < sp.nyquistHz,
                    "passband < design stop edge < Nyquist: the edge margin is inside");
                ok ((long long) sp.L * sp.inRate == (long long) sp.M * sp.outRate, "L:M is the stage's own ratio");
            }
            worstStages = std::max (worstStages, pl.count);
            worstMacs   = std::max (worstMacs, pl.macsPerSecond);

            DeliveryResampler::StagePlan direct;
            const double band = DeliveryResampler::kDefaultPassbandFraction * 0.5 * std::min (a, b);
            if (DeliveryResampler::makeStage ((long long) a, (long long) b, band, DeliveryResampler::kDesignStopbandDb, direct))
                ok (pl.macsPerSecond <= DeliveryResampler::stageMacsPerSecond (direct) + 1e-9,
                    "the chosen route is no worse than the direct one");
        }
    ok (worstStages == 2, "two stages is the most any delivery pair needs");
    ok (worstMacs < 45.0e6, "the worst route stays under 45 MMAC/s per channel");

    group ("plan: the pair the decomposition cannot help, and the one it does");
    {
        const auto up = DeliveryResampler::plan (params (176400.0, 192000.0));
        const auto dn = DeliveryResampler::plan (params (192000.0, 176400.0));
        ok (up.count == 1 && dn.count == 1, "176.4 <-> 192 is a single direct stage");
        ok (up.stage[0].L == 160 && up.stage[0].M == 147, "176.4 -> 192 is 160:147");
        ok (dn.stage[0].L == 147 && dn.stage[0].M == 160, "192 -> 176.4 is 147:160");
        ok (up.stage[0].passbandHz > 79000.0, "and it delivers its 80 kHz of band, not 20");
    }
    {
        const auto pl = DeliveryResampler::plan (params (192000.0, 44100.0));
        DeliveryResampler::StagePlan direct;
        (void) DeliveryResampler::makeStage (192000, 44100, DeliveryResampler::kDefaultPassbandFraction * 0.5 * 44100.0,
                                             DeliveryResampler::kDesignStopbandDb, direct);
        ok (pl.count >= 2, "192 -> 44.1 is decomposed");
        ok ((long long) direct.L * direct.tapsPerPhase > 8 * pl.coefficients,
            "and the decomposition is far smaller than the single stage it replaced");
    }

    group ("plan: refusals, including the finite value that is not a representable one");
    ok (! DeliveryResampler::plan (params (44100.5, 48000.0)).ok, "a non-integer rate is refused");
    ok (! DeliveryResampler::plan (params (48000.0, 44100.5)).ok, "…in either position");
    ok (! DeliveryResampler::plan (params (0.0, 48000.0)).ok, "zero is refused");
    ok (! DeliveryResampler::plan (params (-48000.0, 48000.0)).ok, "negative is refused");
    ok (! DeliveryResampler::plan (params (std::nan (""), 48000.0)).ok, "NaN is refused");
    ok (! DeliveryResampler::plan (params (std::numeric_limits<double>::infinity(), 48000.0)).ok, "infinity is refused");
    ok (! DeliveryResampler::plan (params (1.0e300, 48000.0)).ok, "a finite rate llround cannot hold is refused");
    ok (! DeliveryResampler::plan (params (100.0, 48000.0)).ok, "an absurdly low rate is refused");
    ok (! DeliveryResampler::plan (params (44099.0, 48000.0)).ok, "a ratio past kMaxRatioTerm is refused");
    for (double pf : { std::nextafter (1.0, 0.0), 0.999999, 0.9999999999 })
    {
        // Each of these used to produce an out-of-range double->int conversion or a signed overflow
        // inside plan(), and the first returned ok with zero coefficients.
        auto p = params (48000.0, 44100.0);
        p.passbandFraction = pf;
        const auto pl = DeliveryResampler::plan (p);
        ok (! pl.ok, "a passband fraction whose prototype exceeds the stage bound is refused");
        DeliveryResampler r;
        ok (! r.prepare (p, 1, 512), "…and so is the prepare()");
    }
    {
        auto p = params (48000.0, 44100.0);
        p.passbandFraction = 0.999;
        const auto pl = DeliveryResampler::plan (p);
        ok (pl.ok && pl.coefficients <= DeliveryResampler::kMaxStageCoefficients * pl.count,
            "a steep but admissible passband plans within the stage bound");
    }
    {
        auto p = params (48000.0, 44100.0);
        p.passbandFraction = 1.0;
        ok (! DeliveryResampler::plan (p).ok, "a passband at the full Nyquist is refused");
        p.passbandFraction = 0.0;
        ok (! DeliveryResampler::plan (p).ok, "a zero passband is refused");
        p.passbandFraction = DeliveryResampler::kDefaultPassbandFraction;
        p.stopbandDb = 1000.0;
        ok (! DeliveryResampler::plan (p).ok, "an absurd stopband target is refused");
    }

    group ("plan: identity, and the proportional passband edge");
    {
        const auto pl = DeliveryResampler::plan (params (48000.0, 48000.0));
        ok (pl.ok && pl.identity && pl.count == 0, "equal rates plan as the identity, with no stages");
        approx (DeliveryResampler::latencyOutputSamples (pl), 0.0, 0.0, "…and no latency");
    }
    for (double other : kRates)
    {
        if (other == 44100.0) continue;
        const auto pl = DeliveryResampler::plan (params (other, 44100.0));
        approx (pl.stage[0].passbandHz, 20000.0, 1e-9, "44.1 in the pair puts the edge on 20 kHz");
    }
    approx (DeliveryResampler::plan (params (88200.0, 96000.0)).stage[0].passbandHz, 40000.0, 1e-9,
            "88.2 <-> 96 delivers 40 kHz of band");
}

//==============================================================================
// 2. LAW 11, clause by clause
static void testContract()
{
    group ("prepare: binding, refuses what it cannot honour, and disarms on refusal");
    {
        DeliveryResampler r;
        ok (! r.isPrepared(), "a fresh object is not prepared");
        ok (! r.prepare (params (48000.0, 44100.0), 0, 512), "zero channels is refused");
        ok (! r.prepare (params (48000.0, 44100.0), felitronics::core::kMaxChannels + 1, 512), "more than kMaxChannels is refused");
        ok (! r.prepare (params (48000.0, 44100.0), 2, 0), "a zero block is refused");
        ok (! r.prepare (params (44100.5, 48000.0), 2, 512), "a refused plan refuses the prepare");
        ok (! r.isPrepared(), "none of that left it prepared");
        ok (r.prepare (params (48000.0, 44100.0), 2, 512), "a valid configuration is accepted");
        ok (! r.prepare (params (48000.0, 44100.0), felitronics::core::kMaxChannels + 1, 512), "a later refusal is a refusal");
        ok (! r.isPrepared(), "…and it disarmed the object rather than keeping the old plan");
    }

    group ("process: the check order and the verdict");
    {
        DeliveryResampler r;
        std::vector<float> in (2048, 0.25f), out (8192, 0.0f);
        const float* ip[2] = { in.data(), in.data() };
        float* op[2] = { out.data(), out.data() };
        int n = -1;
        ok (! r.process (ip, 1, 512, op, 8192, n) && n == 0, "unprepared is refused, with no count");
        ok (r.prepare (params (48000.0, 44100.0), 2, 512), "prepare");
        ok (! r.process (ip, 1, -1, op, 8192, n), "a negative n is malformed");
        ok (! r.process (ip, -1, 512, op, 8192, n), "a negative nch is malformed");
        ok (! r.process (ip, 3, 512, op, 8192, n), "more channels than prepared is refused (11b)");
        ok (r.process (ip, 1, 0, op, 8192, n) && n == 0, "n == 0 is accepted and produces nothing");
        ok (r.process (nullptr, 0, 512, nullptr, 0, n) && n == 0, "nch == 0 is accepted with null planes (11a)");
        ok (! r.process (nullptr, 1, 512, op, 8192, n), "a null input plane array is refused");
        ok (! r.process (ip, 1, 512, nullptr, 8192, n), "a null output plane array is refused");
        ok (! r.process (ip, 1, 512, op, 1, n), "an output buffer too small is REFUSED, not truncated");
        ok (r.process (ip, 1, 2048, op, 8192, n) && n > 0, "a call FOUR TIMES maxBlock is processed (11a: length is a capacity)");
    }

    group ("11(a): a long call is bit-identical to the caller chunking it at maxBlock");
    for (const auto& pr : { std::pair<double, double> { 192000.0, 44100.0 }, { 44100.0, 192000.0 }, { 48000.0, 44100.0 } })
    {
        const int maxB = 256, total = 3 * maxB + 7;
        std::vector<float> x ((std::size_t) total);
        for (int i = 0; i < total; ++i) x[(std::size_t) i] = (float) std::sin (0.013 * i) * 0.7f;
        DeliveryResampler a, b;
        ok (a.prepare (params (pr.first, pr.second), 1, maxB) && b.prepare (params (pr.first, pr.second), 1, maxB), "prepare both");
        std::vector<float> ya ((std::size_t) a.maxOutputFor (total)), yb;
        const float* ip[1] = { x.data() };
        float* op[1] = { ya.data() };
        int na = 0;
        ok (a.process (ip, 1, total, op, (int) ya.size(), na), "the long call is accepted");
        ya.resize ((std::size_t) na);
        std::vector<float> buf ((std::size_t) b.maxOutputFor (maxB));
        float* bp[1] = { buf.data() };
        for (int off = 0; off < total; off += maxB)
        {
            const int len = std::min (maxB, total - off);
            const float* cp[1] = { x.data() + off };
            int got = 0;
            ok (b.process (cp, 1, len, bp, (int) buf.size(), got), "caller-chunked call");
            yb.insert (yb.end(), buf.begin(), buf.begin() + got);
        }
        ok (sameBits (ya, yb), "long call == caller-chunked, bit for bit");
        ok (! a.process (ip, 1, total, op, a.maxOutputFor (total) - 1, na), "…and one sample of output short is refused");
    }

    group ("maxOutputFor is a true bound, and never negative");
    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;
            DeliveryResampler r;
            if (! r.prepare (params (a, b), 1, 733)) { ok (false, "prepare"); continue; }
            std::vector<float> in (733, 0.25f), out ((std::size_t) r.maxOutputFor (733));
            const float* ip[1] = { in.data() };
            float* op[1] = { out.data() };
            bool over = false, produced = false;
            for (int rep = 0; rep < 40; ++rep)
            {
                const int want = 1 + (rep * 37) % 733;
                int got = 0;
                if (! r.process (ip, 1, want, op, (int) out.size(), got) || got > r.maxOutputFor (want)) { over = true; break; }
                if (got > 0) produced = true;
            }
            ok (! over && produced, "no call exceeded its own bound, and the converter ran");
        }
    {
        DeliveryResampler r;
        ok (r.prepare (params (44100.0, 192000.0), 1, 1), "prepare");
        ok (r.maxOutputFor (INT_MAX) == INT_MAX, "a bound past INT_MAX saturates instead of wrapping negative");
        ok (r.maxOutputFor (-5) == 0, "a negative n bounds nothing");
    }

    group ("output counts match the independent oracle floor((T*L - 1)/M) + 1, call by call");
    for (const auto& pr : { std::pair<double, double> { 192000.0, 44100.0 }, { 44100.0, 192000.0 }, { 48000.0, 44100.0 }, { 88200.0, 44100.0 } })
    {
        DeliveryResampler r;
        ok (r.prepare (params (pr.first, pr.second), 1, 4096), "prepare");
        const auto& pl = r.currentPlan();
        std::vector<float> in (4096, 0.5f), out ((std::size_t) r.maxOutputFor (4096));
        const float* ip[1] = { in.data() };
        float* op[1] = { out.data() };
        long long T = 0;
        bool match = true;
        for (int len : { 1, 4096, 3, 4095, 17, 2, 1000 })
        {
            int got = 0;
            if (! r.process (ip, 1, len, op, (int) out.size(), got)) { match = false; break; }
            const long long want = cumulativeCount (pl, T + len) - cumulativeCount (pl, T);
            if (got != want) match = false;
            T += len;
        }
        ok (match, "every call's count equals the arithmetic oracle");
    }
    {
        // sol's worked case: 192 -> 88.2 -> 44.1 over blocks [1, 4096, 3, 4095, 17] is [1, 941, 0, 941, 4].
        const auto pl = DeliveryResampler::plan (params (192000.0, 44100.0));
        const long long c[] = { cumulativeCount (pl, 1), cumulativeCount (pl, 4097), cumulativeCount (pl, 4100),
                                cumulativeCount (pl, 8195), cumulativeCount (pl, 8212) };
        ok (c[0] == 1 && c[1] - c[0] == 941 && c[2] - c[1] == 0 && c[3] - c[2] == 941 && c[4] - c[3] == 4,
            "the oracle reproduces the worked counts [1, 941, 0, 941, 4]");
    }

    group ("flush: equals feeding the drain as zeros, ends the programme, and reset() re-arms it");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 }, { 44100.0, 96000.0 } })
    {
        const int T = 3000;
        std::vector<float> x ((std::size_t) T);
        for (int i = 0; i < T; ++i) x[(std::size_t) i] = (float) (0.6 * std::sin (2.0 * kPi * 997.0 * i / pr.first));
        DeliveryResampler a, b;
        ok (a.prepare (params (pr.first, pr.second), 1, 512) && b.prepare (params (pr.first, pr.second), 1, 512), "prepare both");
        std::vector<float> buf ((std::size_t) std::max (a.maxOutputFor (512), a.maxFlushOutput()));
        float* op[1] = { buf.data() };
        for (int off = 0; off < T; off += 512)
        {
            const int len = std::min (512, T - off);
            const float* ip[1] = { x.data() + off };
            int g1 = 0, g2 = 0;
            (void) a.process (ip, 1, len, op, (int) buf.size(), g1);
            (void) b.process (ip, 1, len, op, (int) buf.size(), g2);
        }
        int fa = 0;
        ok (a.flush (1, op, (int) buf.size(), fa) && fa > 0, "flush drains something");
        std::vector<float> tailA (buf.begin(), buf.begin() + fa), tailB;
        const std::vector<float> zeros (512, 0.0f);
        for (int left = b.flushInputSamples(); left > 0; left -= 512)
        {
            const float* zp[1] = { zeros.data() };
            int got = 0;
            (void) b.process (zp, 1, std::min (512, left), op, (int) buf.size(), got);
            tailB.insert (tailB.end(), buf.begin(), buf.begin() + got);
        }
        ok (sameBits (tailA, tailB), "the flushed tail is exactly the converter fed flushInputSamples() zeros");

        int again = -1, n = -1;
        const float* ip[1] = { x.data() };
        ok (a.flush (1, op, (int) buf.size(), again) && again == 0, "a second flush is accepted and writes nothing");
        ok (! a.process (ip, 1, 10, op, (int) buf.size(), n), "process() after flush() is REFUSED — the programme ended");
        a.reset();
        ok (a.process (ip, 1, 10, op, (int) buf.size(), n), "reset() re-arms process()");
        ok (! a.flush (1, op, 1, n), "a flush with too small a buffer is refused");
    }

    group ("identity: exactly equal rates copy the bits, including the hostile ones");
    {
        DeliveryResampler r;
        ok (r.prepare (params (96000.0, 96000.0), 1, 16), "prepare at an equal rate");
        approx (r.latencyOutputSamples(), 0.0, 0.0, "no latency");
        const std::uint32_t pat[] = { 0x00000001u, 0x80000001u, 0x00000000u, 0x80000000u, 0x7f7fffffu,
                                      0xff7fffffu, 0x7fc12345u, 0x7f800000u, 0xff800000u };
        std::vector<float> in (9), out (9, 1.0f);
        std::memcpy (in.data(), pat, sizeof (pat));
        const float* ip[1] = { in.data() };
        float* op[1] = { out.data() };
        int n = 0;
        ok (r.process (ip, 1, 9, op, 9, n) && n == 9, "nine hostile samples in, nine out");
        ok (std::memcmp (in.data(), out.data(), sizeof (pat)) == 0,
            "denormals, signed zeros, FLT_MAX, a NaN payload and both infinities survive bit for bit");
        int f = -1;
        ok (r.flush (1, op, 9, f) && f == 0, "and a copy has no tail");
    }

    group ("channels: each is exactly its mono twin, whatever the other carries");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 }, { 44100.0, 176400.0 } })
    {
        const int T = 4000;
        std::vector<float> xa ((std::size_t) T), xb ((std::size_t) T);
        for (int i = 0; i < T; ++i)
        {
            xa[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 613.0 * i / pr.first));
            xb[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * kPi * 4001.0 * i / pr.first + 1.0));
        }
        double lat = 0.0; bool k1 = false, k2 = false;
        const auto ma = convert (pr.first, pr.second, xa, 300, lat, k1);
        const auto mb = convert (pr.first, pr.second, xb, 300, lat, k2);
        DeliveryResampler st;
        ok (k1 && k2 && st.prepare (params (pr.first, pr.second), 2, 300), "prepare");
        std::vector<float> ya, yb, b0 ((std::size_t) std::max (st.maxOutputFor (300), st.maxFlushOutput())), b1 (b0.size());
        float* op[2] = { b0.data(), b1.data() };
        for (int off = 0; off < T; off += 300)
        {
            const int len = std::min (300, T - off);
            const float* ip[2] = { xa.data() + off, xb.data() + off };
            int got = 0;
            (void) st.process (ip, 2, len, op, (int) b0.size(), got);
            ya.insert (ya.end(), b0.begin(), b0.begin() + got);
            yb.insert (yb.end(), b1.begin(), b1.begin() + got);
        }
        int got = 0;
        (void) st.flush (2, op, (int) b0.size(), got);
        ya.insert (ya.end(), b0.begin(), b0.begin() + got);
        yb.insert (yb.end(), b1.begin(), b1.begin() + got);
        ok (sameBits (ya, ma), "channel 0 == its mono conversion, bit for bit");
        ok (sameBits (yb, mb), "channel 1 == its mono conversion, bit for bit");
    }

    group ("11a: a stopped channel does not replay, and a gap is spent on the clock");
    {
        // sol's fixture: channel 1 carries DC, stops for a block, then both are silent. Its history must
        // be gone: the reference (fed silence instead of stopped) is exactly zero after 512 > 228 taps.
        DeliveryResampler r;
        ok (r.prepare (params (48000.0, 44100.0), 2, 512), "prepare");
        std::vector<float> zeros (512, 0.0f), ones (512, 1.0f), o0 (1024), o1 (1024);
        float* op[2] = { o0.data(), o1.data() };
        int n = 0;
        const float* c1[2] = { zeros.data(), ones.data() };
        ok (r.process (c1, 2, 512, op, 1024, n), "both channels run");
        const float* c2[1] = { zeros.data() };
        ok (r.process (c2, 1, 512, op, 1024, n), "channel 1 stops");
        const float* c3[2] = { zeros.data(), zeros.data() };
        ok (r.process (c3, 2, 512, op, 1024, n), "channel 1 returns");
        double peak = 0.0;
        for (int i = 0; i < n; ++i) peak = std::max (peak, (double) std::fabs (o1[(std::size_t) i]));
        ok (peak == 0.0, "the returning channel is exactly silent: its memory was dropped, not replayed");
        ok (n == 471, "…and the stream clock ran through the narrower call (471 outputs, as the reference)");
    }
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 }, { 44100.0, 192000.0 } })
    {
        // nch == 0 for G samples must be bit-identical, afterwards, to feeding G samples of silence —
        // which checks the closed-form clock against the push-by-push automaton.
        const int G = 4099, B = 512;
        std::vector<float> x ((std::size_t) B);
        for (int i = 0; i < B; ++i) x[(std::size_t) i] = (float) (0.8 * std::sin (0.21 * i));
        DeliveryResampler a, b;
        ok (a.prepare (params (pr.first, pr.second), 2, 1024) && b.prepare (params (pr.first, pr.second), 2, 1024), "prepare both");
        std::vector<float> p0 ((std::size_t) b.maxOutputFor (G)), p1 (p0.size());
        float* op[2] = { p0.data(), p1.data() };
        const float* ip[2] = { x.data(), x.data() };
        int na = 0, nb = 0;
        (void) a.process (ip, 2, B, op, (int) p0.size(), na);
        (void) b.process (ip, 2, B, op, (int) p0.size(), nb);
        ok (a.process (nullptr, 0, G, nullptr, 0, na) && na == 0, "the gap is accepted and writes nothing");
        const std::vector<float> zeros ((std::size_t) G, 0.0f);
        const float* zp[2] = { zeros.data(), zeros.data() };
        ok (b.process (zp, 2, G, op, (int) p0.size(), nb), "the reference is fed the gap as silence");
        std::vector<float> ya0 ((std::size_t) a.maxOutputFor (B)), ya1 (ya0.size()), yb0 (ya0.size()), yb1 (ya0.size());
        float* oa[2] = { ya0.data(), ya1.data() };
        float* ob[2] = { yb0.data(), yb1.data() };
        (void) a.process (ip, 2, B, oa, (int) ya0.size(), na);
        (void) b.process (ip, 2, B, ob, (int) yb0.size(), nb);
        ya0.resize ((std::size_t) na); ya1.resize ((std::size_t) na);
        yb0.resize ((std::size_t) nb); yb1.resize ((std::size_t) nb);
        ok (na == nb, "after the gap both emit the same count");
        ok (sameBits (ya0, yb0) && sameBits (ya1, yb1), "…and the same samples: a gap is a gap, not a pause");
    }

    group ("reset: clears the history AND the clock");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 } })
    {
        DeliveryResampler used, fresh;
        ok (used.prepare (params (pr.first, pr.second), 1, 600) && fresh.prepare (params (pr.first, pr.second), 1, 600), "prepare both");
        std::vector<float> ones (600, 1.0f), imp (600, 0.0f), ou ((std::size_t) used.maxOutputFor (600)), of (ou.size());
        imp[0] = 1.0f;
        float* up[1] = { ou.data() };
        float* fp[1] = { of.data() };
        int nu = 0, nf = 0;
        const float* onep[1] = { ones.data() };
        (void) used.process (onep, 1, 357, up, (int) ou.size(), nu);             // leave a mid-phase clock
        used.reset();
        const float* ip[1] = { imp.data() };
        (void) used.process (ip, 1, 600, up, (int) ou.size(), nu);
        (void) fresh.process (ip, 1, 600, fp, (int) of.size(), nf);
        ou.resize ((std::size_t) nu); of.resize ((std::size_t) nf);
        ok (sameBits (ou, of), "after reset() the stream is bit-identical to a fresh object's");
    }

    group ("RT-safety: process(), a long call, a gap, flush() and reset() allocate nothing");
    {
        DeliveryResampler r;
        ok (r.prepare (params (192000.0, 44100.0), 2, 1024), "prepare the worst cascade");
        std::vector<float> in (4096, 0.1f), out ((std::size_t) std::max (r.maxOutputFor (4096), r.maxFlushOutput()));
        const float* ip[2] = { in.data(), in.data() };
        float* op[2] = { out.data(), out.data() };
        int n = 0;
        (void) r.process (ip, 2, 1024, op, (int) out.size(), n);
        const long long before = g_allocs.load();
        for (int k = 0; k < 4; ++k) (void) r.process (ip, 2, 1024, op, (int) out.size(), n);
        (void) r.process (ip, 2, 4096, op, (int) out.size(), n);         // chunked
        (void) r.process (ip, 1, 1024, op, (int) out.size(), n);         // falling edge
        (void) r.process (nullptr, 0, 100000, nullptr, 0, n);            // gap, closed-form clock
        (void) r.flush (2, op, (int) out.size(), n);
        r.reset();
        okNoAlloc (g_allocs.load() == before, "none of those allocated");
    }
}

//==============================================================================
// 3. LATENCY — measured back out of the carrier phase, on a single stage AND on a cascade. Recomputing
// the formula from the plan's own fields is a tautology; only the audio can see a mis-composed cascade.
static void testLatency()
{
    group ("latency: measured from the carrier phase equals the exact rational formula");
    for (const auto& pr : { std::pair<double, double> { 48000.0, 44100.0 }, { 192000.0, 44100.0 },
                            { 44100.0, 192000.0 }, { 176400.0, 48000.0 } })
    {
        const double f = 100.0;                                      // period longer than any latency here
        const int n0 = (int) pr.first * 3;
        std::vector<float> x ((std::size_t) n0);
        for (int i = 0; i < n0; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * f * i / pr.first);
        double lat = 0.0; bool okc = false;
        const auto y = convert (pr.first, pr.second, x, 1024, lat, okc);
        ok (okc, "conversion accepted");
        const int start = (int) std::ceil (lat) + (int) pr.second / 2;
        const double d = measuredDelay (y, pr.second, f, start);
        approx (d, lat, 0.01, "the delay the audio carries is the formula's, to a hundredth of a sample");
    }
}

//==============================================================================
// 4. THE ROUND-TRIP NULL — 48 -> 44.1 -> 48 against an ANALYTIC oracle, at the exact rational delay.
static void testRoundTrip()
{
    group ("round trip 48 -> 44.1 -> 48 nulls against the analytic signal");
    static const double F[]  = { 97.0, 440.0, 1000.0, 3170.0, 7000.0, 11000.0, 15000.0, 19000.0 };
    static const double PH[] = { 0.0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 4.9 };
    auto tone = [&] (double t48, double amp)
    {
        double s = 0.0;
        for (int j = 0; j < 8; ++j) s += std::sin (2.0 * kPi * F[j] * t48 / 48000.0 + PH[j]);
        return amp * s / 8.0;
    };
    const int n0 = 48000 * 4;
    const double ramp = 48000.0 * 0.05;                             // taper both ends — fixture rule 1
    auto run = [&] (double amp, double& peak, double& rms, double& dTot) -> bool
    {
        std::vector<float> x ((std::size_t) n0);
        for (int i = 0; i < n0; ++i)
        {
            double w = 1.0;
            if ((double) i < ramp)          w = 0.5 - 0.5 * std::cos (kPi * (double) i / ramp);
            if ((double) i > n0 - 1 - ramp) w = 0.5 - 0.5 * std::cos (kPi * (double) (n0 - 1 - i) / ramp);
            x[(std::size_t) i] = (float) (tone ((double) i, amp) * w);
        }
        double l1 = 0.0, l2 = 0.0; bool ok1 = false, ok2 = false;
        const auto mid = convert (48000.0, 44100.0, x, 1024, l1, ok1);
        const auto rt  = convert (44100.0, 48000.0, mid, 1024, l2, ok2);
        dTot = l1 * (48000.0 / 44100.0) + l2;
        const int guard = 48000 / 2;
        double mx = 0.0, acc = 0.0; long long m = 0;
        for (int i = guard; i < (int) rt.size() - guard && i < n0 - guard; ++i)
        {
            const double e = (double) rt[(std::size_t) i] - tone ((double) i - dTot, amp);
            mx = std::max (mx, std::fabs (e)); acc += e * e; ++m;
        }
        peak = mx; rms = std::sqrt (acc / (double) std::max (1LL, m));
        return ok1 && ok2 && m > 100000;
    };

    double peak = 0.0, rms = 0.0, dTot = 0.0;
    ok (run (0.9, peak, rms, dTot), "both legs accepted, and the null measured a real stretch");
    {
        const auto d1 = DeliveryResampler::plan (params (48000.0, 44100.0));
        const auto d2 = DeliveryResampler::plan (params (44100.0, 48000.0));
        approx (dTot, (double) d1.stage[0].halfLen + (double) d2.stage[0].halfLen * (48000.0 / 44100.0), 1e-9,
                "the round trip's delay is the two stages' own integer delays carried to 48 kHz");
    }
    ok (dbfs (peak) <= -120.0, "round-trip peak residual at or under -120 dBFS (the requirement)");
    ok (dbfs (peak) <= -140.0, "round-trip peak residual at or under -140 dBFS (delivered -143.3)");
    ok (dbfs (rms) <= -150.0, "round-trip RMS residual at or under -150 dBFS (delivered -160.3)");
    std::printf ("      round trip: peak %.2f dBFS, rms %.2f dBFS, delay %.6f samples\n", dbfs (peak), dbfs (rms), dTot);

    // Recorded, NOT asserted: how the residual moves with level cannot say which term dominates it —
    // see "WHAT THE ROUND TRIP MEASURES" in the header. Two earlier assertions of a mechanism both broke.
    double peakQ = 0.0, rmsQ = 0.0, dQ = 0.0;
    (void) run (0.09, peakQ, rmsQ, dQ);
    std::printf ("      at -20 dB:  peak %.2f dBFS, rms %.2f dBFS (for the record only)\n", dbfs (peakQ), dbfs (rmsQ));
}

//==============================================================================
// 5. SLICING AND FILE LENGTH
static void testSlicingAndLength()
{
    group ("the output does not depend on how the input was sliced");
    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;
            const int n0 = 6000;
            std::vector<float> x ((std::size_t) n0);
            for (int i = 0; i < n0; ++i)
                x[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 613.0 * i / a) + 0.3 * std::sin (2.0 * kPi * 4001.0 * i / a));
            double l = 0.0; bool k1 = false, k2 = false, k3 = false;
            const auto y1 = convert (a, b, x, n0, l, k1);
            const auto y2 = convert (a, b, x, 1, l, k2);
            const auto y3 = convert (a, b, x, 257, l, k3);
            ok (k1 && k2 && k3, "all three slicings accepted");
            ok (sameBits (y1, y2) && sameBits (y1, y3), "…and bit-identical, count and samples");
        }

    group ("process + flush emit the whole file, not the whole file minus the group delay");
    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;
            const int T = 20000;
            std::vector<float> x ((std::size_t) T);
            for (int i = 0; i < T; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * kPi * 501.0 * i / a);
            double lat = 0.0; bool okc = false;
            const auto y = convert (a, b, x, 512, lat, okc);
            const double want = lat + std::floor ((double) T * b / a);
            ok (okc && (double) y.size() >= want, "enough output to drop the latency and keep the file");
            // Each round-up in the drain is an extra INPUT sample, worth outRate/inRate outputs; a drain
            // that ran to maxBlock would pad every file with silence and nothing else here would see it.
            ok ((double) y.size() < want + 4.0 + 4.0 * (b / a), "…and the drain does not overrun");
        }
}

//==============================================================================
// 6. NON-FINITE INPUT washes out in bounded time
static void testHygiene()
{
    group ("a non-finite sample does not poison the converter for ever");
    DeliveryResampler r;
    ok (r.prepare (params (48000.0, 44100.0), 1, 256), "prepare");
    std::vector<float> in (256, 0.0f), out (1024, 0.0f);
    const float* ip[1] = { in.data() };
    float* op[1] = { out.data() };
    int n = 0;
    in[100] = std::numeric_limits<float>::quiet_NaN();
    ok (r.process (ip, 1, 256, op, 1024, n), "a NaN-bearing block is still accepted");
    std::fill (in.begin(), in.end(), 0.0f);
    bool clean = false;
    for (int k = 0; k < 64 && ! clean; ++k)
    {
        if (! r.process (ip, 1, 256, op, 1024, n)) { ok (false, "process"); break; }
        clean = true;
        for (int i = 0; i < n; ++i) if (! std::isfinite (out[(std::size_t) i])) clean = false;
    }
    ok (clean, "the NaN washes out of the history in bounded time");
}

int main()
{
    std::printf ("felitronics::core::DeliveryResampler — contract, plan, latency, round trip\n");
    testPlan();
    testContract();
    testLatency();
    testRoundTrip();
    testSlicingAndLength();
    testHygiene();
    return felitronics::test::report();
}
