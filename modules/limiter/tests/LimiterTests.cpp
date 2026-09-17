// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the true-peak limiter: the GUARANTEE (output true-peak, measured by an
// INDEPENDENT oversampler, stays at/below the ceiling even when the input's inter-sample peaks exceed
// it), transparency below the ceiling, latency, and no-allocation-in-process().

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;

// Independent true-peak (dBTP) of a baseband buffer, via a FRESH 8x oversampler (not the limiter's).
static double measureTruePeakDb (const std::vector<float>& x)
{
    oversampling::PolyphaseOversampler m; (void) m.prepare (8, 1, 32);
    std::vector<float> osb (x.size() * 8);
    const float* xi[1] { x.data() }; float* oo[1] { osb.data() };
    m.upsample (xi, 1, (int) x.size(), oo);
    double mx = 0.0;
    for (std::size_t i = x.size(); i + x.size() < osb.size(); ++i) mx = std::max (mx, (double) std::fabs (osb[i]));
    return 20.0 * std::log10 (mx > 1e-9 ? mx : 1e-9);
}

static double rmsTail (const std::vector<float>& v, double frac)
{
    const int from = (int) (v.size() * (1.0 - frac));
    double s = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { s += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

//==============================================================================
// THE DUAL RELEASE (M2), against a reference of a different construction: the gain law recomputed from the limiter's own
// reconstructed-peak tap, its two windows maxima taken offline by van Herk / Gil-Werman (prefix and suffix maxima over
// blocks of the window's length) rather than by the class's deque, the switches of `dualRelease` applied where the
// contract puts them, and the poison flush at every chunk end. Compared with the gain-reduction tap bit for bit.
namespace dual
{
// max over x[i-W+1 .. i], the samples before 0 absent.
std::vector<float> windowMax (const std::vector<float>& x, int W)
{
    const std::size_t n = x.size(), w = (std::size_t) W;
    std::vector<float> g (n), h (n), out (n);
    for (std::size_t i = 0; i < n; ++i) g[i] = (i % w == 0) ? x[i] : std::max (g[i - 1], x[i]);
    for (std::size_t i = n; i-- > 0;) h[i] = (i + 1 == n || (i + 1) % w == 0) ? x[i] : std::max (h[i + 1], x[i]);
    for (std::size_t i = 0; i < n; ++i) out[i] = (i + 1 < w) ? g[i] : std::max (h[i + 1 - w], g[i]);
    return out;
}

float coefFor (double ms, double fs, int F)
{
    const double t = std::max (ms * 0.001 * fs * (double) F, 8.0 * (double) F);
    float c = (float) std::exp (-1.0 / t);
    if (! (c < 1.0f)) c = std::nextafterf (1.0f, 0.0f);
    return c;
}

struct Segment { int baseband; bool dual; };   // a stretch of the stream processed with `dualRelease` as given

struct Run
{
    std::vector<float> gr, peak;               // the class's taps, OVERSAMPLED
    std::vector<float> out;                    // channel 0 of the output
    double ceiling = 0.0;
    int F = 0, look = 0;
    bool ok = false;
};

// Streams `x` (mono) through a limiter in blocks of `block` baseband samples, `dualRelease` switched per segment.
Run render (const std::vector<float>& x, double fs, double ceilingDb, double relMs, double slowMs,
            const std::vector<Segment>& segs, int block)
{
    Run r;
    limiter::TruePeakLimiter lim;
    if (! lim.prepare (fs, block, 1, {})) return r;
    r.F = lim.oversampleFactor(); r.look = lim.lookaheadSamples();
    const int n = (int) x.size();
    r.gr.assign ((std::size_t) n * (std::size_t) r.F, 0.0f); r.peak = r.gr; r.out = x;
    int at = 0;
    for (const Segment& sg : segs)
    {
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = ceilingDb; p.releaseMs = relMs; p.slowReleaseMs = slowMs; p.dualRelease = sg.dual;
        lim.setParams (p);
        r.ceiling = lim.effectiveCeilingDbTp();
        for (int done = 0; done < sg.baseband; )
        {
            const int m = std::min (block, sg.baseband - done);
            float* io[1] { r.out.data() + at };
            limiter::TruePeakLimiterTap tap;
            tap.gainReductionDb = r.gr.data() + (std::size_t) at * (std::size_t) r.F;
            tap.linkedPeakLin   = r.peak.data() + (std::size_t) at * (std::size_t) r.F;
            tap.capacity        = m * r.F;
            if (! lim.process (io, 1, m, tap)) return r;
            at += m; done += m;
        }
    }
    r.ok = at == n;
    return r;
}

// The contract's gain law over the tap's peaks. `chunkOs` is where the class flushes: the ends of its process chunks.
std::vector<float> reference (const Run& r, double fs, double relMs, double slowMs, const std::vector<Segment>& segs, int block)
{
    const std::size_t n = r.peak.size();
    const std::vector<float> smax = windowMax (r.peak, r.look * r.F + 1);
    std::vector<float> raw (n);
    for (std::size_t i = 0; i < n; ++i)
    {
        double v = r.ceiling - core::gainToDb ((double) smax[i]);
        if (v > 0.0) v = 0.0;
        raw[i] = (float) v;
    }
    const int W = limiter::TruePeakLimiter::slowWindowSamplesFor (fs, r.F);
    const std::vector<float> held = windowMax (raw, W);
    const float cF = coefFor (relMs, fs, r.F), cS = coefFor (slowMs, fs, r.F);
    std::vector<float> gr (n);
    float grF = 0.0f, grS = 0.0f;
    bool on = false;
    std::size_t onAt = 0, i = 0;
    for (const Segment& sg : segs)
    {
        if (sg.dual && ! on) { grS = 0.0f; onAt = i; }
        else if (! sg.dual && on) grF = std::min (grF, grS);
        on = sg.dual;
        for (int done = 0; done < sg.baseband; )
        {
            const int m = std::min (block, sg.baseband - done);
            for (std::size_t k = 0; k < (std::size_t) m * (std::size_t) r.F; ++k, ++i)
            {
                grF = std::min (raw[i], grF * cF);
                float g = grF;
                if (on)
                {
                    grS = std::min (i + 1 - onAt >= (std::size_t) W ? held[i] : 0.0f, grS * cS);
                    g = std::min (grF, grS);
                }
                gr[i] = g;
            }
            core::flushPoison (grF);
            core::flushPoison (grS);
            done += m;
        }
    }
    return gr;
}

std::size_t differ (const std::vector<float>& a, const std::vector<float>& b)
{
    std::size_t d = a.size() == b.size() ? 0u : 1u + std::max (a.size(), b.size());
    for (std::size_t i = 0; i < std::min (a.size(), b.size()); ++i) if (std::memcmp (&a[i], &b[i], sizeof (float)) != 0) ++d;
    return d;
}

// A 1 kHz tone at 48 kHz holds a crest in every millisecond of lookahead, so its required reduction does not return to
// 0 dB while it sounds: the one fixture on which the slow envelope's window fills.
void tone (std::vector<float>& x, double fs, double from, double seconds, double amp)
{
    const int a = (int) std::lround (from * fs), b = (int) std::lround ((from + seconds) * fs);
    for (int i = a; i < b && i < (int) x.size(); ++i)
        x[(std::size_t) i] = (float) (amp * std::sin (2.0 * core::kPi * 1000.0 * (double) i / fs));
}
}   // namespace dual

int main()
{
    std::printf ("felitronics::limiter tests\n");
    const double sr = 48000.0;

    // --- THE GUARANTEE: a signal whose inter-sample peaks exceed the ceiling → output true-peak ≤ ceiling ---
    test::group ("Limiter true-peak guarantee (ISP)");
    {
        const int n = 8192; const double f = sr * 0.25, A = 1.0;     // sample peaks ~0.765, TRUE peak ~1.0 (0 dBTP)
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr + 0.7));
        const double inTp = measureTruePeakDb (x);
        test::ok (inTp > -1.0 + 0.3, "input true-peak is genuinely above the ceiling (real ISP case)");

        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 1, { 1.0, 4, 32 });   // topology is a prepare-time thing now
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0;
        lim.setParams (p);
        felitronics::test::run (lim.process (ch, 1, n));

        const double outTp = measureTruePeakDb (y);
        test::ok (outTp <= -1.0 + 0.5, "output true-peak ≤ ceiling (+0.5 dB downsample-ripple margin)");
        test::ok (lim.gainReductionDb() < -0.2, "limiter actually engaged");
    }

    // --- transparency: a signal below the ceiling passes ~unchanged (just the latency) ---
    test::group ("Limiter transparent below ceiling");
    {
        const int n = 4096; const double f = 1000.0, A = 0.3;        // ~-10 dBFS, well below -1 dBTP
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 1, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0;
        lim.setParams (p);
        felitronics::test::run (lim.process (ch, 1, n));
        test::approx (rmsTail (y, 0.3) / rmsTail (x, 0.3), 1.0, 0.05, "amplitude preserved (transparent)");
        test::ok (lim.gainReductionDb() > -0.3, "no meaningful gain reduction below ceiling");
    }

    // --- latency = oversampler round-trip + lookahead (baseband) ---
    test::group ("Limiter latency");
    {
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, 512, 2, { 1.0, 4, 32 });
        const int look = (int) std::lround (1.0 * 0.001 * sr);       // 48 baseband
        test::ok (lim.latencySamples() == (32 - 1) + look, "latency = (tpp-1) + lookahead");
        // ...and at the DEFAULT topology, which is the one a product gets and which nothing here pinned
        // while every call in this file spelled 32 out. 79 -> 111 samples at 48 kHz is a host-visible
        // resynchronisation, so it is stated as a number rather than inherited from a formula.
        limiter::TruePeakLimiter def; (void) def.prepare (sr, 512, 2, {});
        test::ok (limiter::TruePeakLimiterConfig {}.tapsPerPhase == 64,
                  "TruePeakLimiterConfig defaults to 64 taps/phase (was 32 — see PolyphaseOversampler.h)");
        test::ok (def.latencySamples() == 63 + look && def.latencySamples() == 111,
                  "default topology latency = 63 + 48 = 111 samples at 48 kHz (was 79)");
    }

    // --- no allocation during process() ---
    test::group ("Limiter no-alloc in process()");
    {
        const int n = 512;
        std::vector<float> a (n, 0.6f), b (n, 0.6f);
        float* ch[2] { a.data(), b.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 2, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0;
        lim.setParams (p);
        const long long before = alloc::count.load();
        felitronics::test::run (lim.process (ch, 2, n));
        felitronics::test::run (lim.process (ch, 2, n));
        const long long after = alloc::count.load();
        test::okNoAlloc (after == before, "process() performed zero heap allocations");
    }

    // --- lifecycle/misuse: process before/after a FAILED prepare + oversized blocks must not OOB osBuf ---
    test::group ("TruePeakLimiter: reject process before / after failed prepare, and oversized blocks");
    {
        limiter::TruePeakLimiter lim;                                // NOT prepared (maxCh == 0)
        float a[64] {}, b[64] {}; float* io[2] { a, b };
        test::ok (! lim.process (io, 2, 16), "process() before prepare() is REFUSED (law 11)");
        test::ok (! lim.prepare (48000.0, 16, 2, { 1.0, 4, 2 }),      // tapsPerPhase=2 < 4 → rejected
                  "prepare() REPORTS an unusable configuration instead of half-building");
        test::ok (! lim.process (io, 2, 16), "...and after a FAILED prepare too, still reported");
        test::ok (lim.prepare (48000.0, 16, 2, { 1.0, 4, 32 }), "and accepts a valid one");
        felitronics::test::run (lim.process (io, 2, 16));                                     // works
        felitronics::test::run (lim.process (io, 2, 64));                                     // 64 > maxBlock 16 → CHUNKED now, not dropped
        test::ok (true, "no OOB across failed-prepare / oversized-block process (ASan/UBSan is the check)");
    }

    // --- FALSIFICATION: the SlidingMax deque must survive a strictly-decreasing run that fills it ---
    // (insert-before-expire with capacity == W overwrites the head — the current max — on the W+1-th push)
    test::group ("SlidingMax survives a full deque (4,3,2,1,0 @ W=4)");
    {
        limiter::detail::SlidingMax sm; (void) sm.prepare (4);
        const float in[5]   = { 4.0f, 3.0f, 2.0f, 1.0f, 0.0f };
        const float want[5] = { 4.0f, 4.0f, 4.0f, 4.0f, 3.0f };   // window of 4 → last covers {3,2,1,0}
        for (int i = 0; i < 5; ++i)
            test::approx ((double) sm.push (in[i]), (double) want[i], 0.0, "sliding max after push #" + std::to_string (i));
    }

    // --- SlidingMax == brute-force window max over a hostile sequence (long decreasing ramps + noise) ---
    test::group ("SlidingMax == brute-force window max (hostile sequence)");
    {
        const int W = 5;
        limiter::detail::SlidingMax sm; (void) sm.prepare (W);
        unsigned long long s = 42;
        auto rnd = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 65536.0f; };
        std::vector<float> xs;
        for (int i = 0; i < 40; ++i)
        {
            if (i % 4 == 0) { float v = 1.0f + rnd(); for (int k = 0; k < 50; ++k) { xs.push_back (v); v *= 0.98f; } }  // strictly-decreasing run ≫ W
            else            { for (int k = 0; k < 30; ++k) xs.push_back (rnd()); }
        }
        bool allOk = true;
        for (int i = 0; i < (int) xs.size() && allOk; ++i)
        {
            const float got = sm.push (xs[i]);
            float want = 0.0f;
            for (int j = std::max (0, i - (W - 1)); j <= i; ++j) want = std::max (want, xs[j]);
            allOk = (got == want);
        }
        test::ok (allOk, "sliding max exact for every push (W=5, decreasing runs fill the deque)");
    }

    // --- FALSIFICATION: ceiling guarantee on a >20 ms decaying ramp (deque-overflow end-to-end) ---
    test::group ("Limiter ceiling guarantee on a decaying ramp");
    {
        const int n = 8192;
        std::vector<float> x (n, 0.0f);
        const int rampLen = (int) (0.025 * sr);                      // 25 ms strictly-decreasing > the 20 ms window
        for (int i = 0; i < rampLen; ++i) x[i] = 1.0f - 0.7f * (float) i / (float) rampLen;   // 1.0 → 0.3
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 1, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = 20.0 * std::log10 (0.25); p.releaseMs = 0.5;
        lim.setParams (p);
        felitronics::test::run (lim.process (ch, 1, n));
        double mx = 0.0; for (float v : y) mx = std::max (mx, (double) std::fabs (v));
        test::ok (mx <= 0.25 * 1.03, "every output sample ≤ ceiling on the smooth ramp (got max " + std::to_string (mx) + ")");
    }

    // --- FALSIFICATION: release must start after the LOOKAHEAD window, not the fixed 20 ms max window ---
    test::group ("Limiter release follows the lookahead window (no 20 ms hold)");
    {
        const int n = 4096;
        std::vector<float> x (n, 0.0f);
        for (int i = 0; i < 96; ++i)  x[i] = (float) std::sin (2.0 * core::kPi * 0.25 * i + 0.7);            // 2 ms 0 dBFS fs/4 burst
        for (int i = 96; i < n; ++i)  x[i] = (float) (0.25 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));  // then a −12 dB tone
        std::vector<float> y = x; float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 1, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -6.0; p.releaseMs = 1.0;
        lim.setParams (p);
        felitronics::test::run (lim.process (ch, 1, n));
        auto rmsWin = [] (const std::vector<float>& v, int a, int b) {
            double s2 = 0.0; for (int i = a; i < b; ++i) s2 += (double) v[i] * v[i]; return std::sqrt (s2 / std::max (1, b - a)); };
        const double got  = rmsWin (y, (int) (0.010 * sr), (int) (0.020 * sr));   // 10–20 ms: burst long gone at 1 ms lookahead
        const double want = rmsWin (x, (int) (0.010 * sr), (int) (0.020 * sr));
        test::approx (got / want, 1.0, 0.06, "tone recovered a few ms after the burst (gain not held for 20 ms)");
    }

    // --- non-finite params must not poison the stream (house rule: clamp + reject non-finite) ---
    test::group ("Limiter non-finite params rejected");
    {
        const int n = 1024;
        std::vector<float> y (n);
        for (int i = 0; i < n; ++i) y[i] = (float) (0.8 * std::sin (2.0 * core::kPi * 997.0 * i / sr));
        float* ch[1] { y.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 1, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p;
        p.ceilingDbTp = std::numeric_limits<double>::quiet_NaN();
        p.releaseMs   = std::numeric_limits<double>::quiet_NaN();
        lim.setParams (p);
        felitronics::test::run (lim.process (ch, 1, n));
        bool finite = true; for (float v : y) finite &= (bool) std::isfinite (v);
        test::ok (finite, "NaN ceiling/release → finite output");
        test::ok (lim.latencySamples() >= 0, "latency stays non-negative");
        // The lookahead moved into the prepare-time config, so a non-finite one is REFUSED rather than
        // silently swallowed — and a non-finite SAMPLE RATE too: `rate <= 0.0` is false for NaN, which
        // used to walk straight into std::lround(NaN).
        limiter::TruePeakLimiter bad;
        test::ok (! bad.prepare (sr, n, 1, { std::numeric_limits<double>::quiet_NaN(), 4, 32 }),
                  "prepare() refuses a non-finite lookahead");
        test::ok (! bad.prepare (std::numeric_limits<double>::quiet_NaN(), n, 1), "prepare() refuses a NaN sample rate");
        test::ok (! bad.prepare (std::numeric_limits<double>::infinity(), n, 1), "prepare() refuses an infinite sample rate");
        test::ok (! bad.isPrepared(), "and stays unprepared after all three");
    }

    // --- latency query before prepare() must not report the unprepared oversampler's tpp-1 == -1 ---
    test::group ("Limiter latencySamples before prepare == 0");
    {
        limiter::TruePeakLimiter lim;
        test::ok (lim.latencySamples() == 0, "unprepared limiter reports 0 latency, not -1");
    }

    // THE maxBlock CAP HAS TO BOUND THE ALLOCATION IT EXISTS FOR. prepare() clamps maxBlock to
    // kMaxBlock (1 << 20) and then chunks by the clamped value — but the oversampled scratch used to be
    // sized from the UNCLAMPED argument, so the cap bounded nothing. It is reachable by ordinary use,
    // not a hostile one: the header tells an offline caller that sizing maxBlock to a whole file is
    // normal, and a 10-minute stereo file at 48 kHz then asked for ~460 MB per channel. Counting BYTES
    // is what discriminates — the number of allocations is identical either way.
    test::group ("prepare() respects its own maxBlock cap in BYTES, not only in chunking");
    {
        auto bytesFor = [] (int maxBlock) {
            const long long before = alloc::rawBytes.load();
            {
                limiter::TruePeakLimiter lim;
                test::ok (lim.prepare (48000.0, maxBlock, 1, {}), "prepare(maxBlock = " + std::to_string (maxBlock) + ")");
            }
            return alloc::rawBytes.load() - before;
        };
        const long long atCap   = bytesFor (1 << 20);
        const long long overCap = bytesFor (1 << 23);          // eight times the cap
        test::ok (overCap <= atCap + (atCap / 8),
                  "asking for 8x the cap allocates no more than asking for the cap ("
                  + std::to_string (overCap / 1024) + " KB vs " + std::to_string (atCap / 1024) + " KB)");
        // ...and the clamp is still only a scratch size: a call far larger than maxBlock is processed
        // whole, so nothing was traded away for the bound.
        limiter::TruePeakLimiter lim;
        test::ok (lim.prepare (48000.0, 1 << 23, 1, {}), "prepare with an over-cap maxBlock still succeeds");
        std::vector<float> x ((std::size_t) 200000, 0.9f);
        float* io[1] { x.data() };
        limiter::TruePeakLimiterParams pr; pr.ceilingDbTp = -6.0; lim.setParams (pr);
        felitronics::test::run (lim.process (io, 1, (int) x.size()));
        double peak = 0.0;
        bool finite = true;
        for (float v : x) { finite &= (bool) std::isfinite (v); peak = std::max (peak, (double) std::fabs (v)); }
        test::ok (finite && peak > 0.0 && peak < 0.9, "...and a 200000-sample call is still limited, whole");
    }

    // P31: the topology is part of the budget. prepare() must ask for exactly what storageFor() publishes
    // under either oversampler, and the cascade's rings and scratch are not the Kaiser one's.
    test::group ("prepare() asks for exactly its published budget under both oversampler topologies");
    {
        std::uint64_t kaiserBytes = 0, cascadeBytes = 0;
        for (double rate : { 44100.0, 48000.0, 96000.0 })     // the cascade designs FROM the rate: one rate proves one design
            for (auto topo : { oversampling::Topology::Kaiser, oversampling::Topology::Cascade })
            {
                limiter::TruePeakLimiterConfig cfg; cfg.oversampleFactor = 8; cfg.topology = topo;
                limiter::TruePeakLimiter::Storage st;
                const bool okSt = limiter::TruePeakLimiter::storageFor (rate, 1024, 2, cfg, st);
                limiter::TruePeakLimiter lim;
                const long long before = alloc::bytes.load();
                const bool okP = lim.prepare (rate, 1024, 2, cfg);
                const long long got = alloc::bytes.load() - before;
                const bool kaiser = topo == oversampling::Topology::Kaiser;
                test::ok (okSt && okP && got == (long long) st.bytes(),
                          std::to_string ((int) rate) + " Hz " + (kaiser ? "Kaiser" : "Cascade") + " 8x: asked " + std::to_string (got)
                          + " B, published " + std::to_string (st.bytes()) + " B");
                if (rate == 44100.0) (kaiser ? kaiserBytes : cascadeBytes) = st.bytes();
            }
        test::ok (kaiserBytes != cascadeBytes, "and the two budgets really are different topologies ("
                                               + std::to_string (kaiserBytes) + " vs " + std::to_string (cascadeBytes) + " B)");
    }

    // M2 — THE DUAL RELEASE. The reference above, bit for bit, over tone stretches of 20, 45, 60, 200 and 500 ms at two
    // depths, silence between, dense noise, switches on and off mid-stream, and blocks of 1, 7, 480 and 4096 samples.
    test::group ("dual release: the gain law nulls against an offline reference, switches and chunk ends included");
    {
        const int n = (int) (sr * 3.0);
        std::vector<float> x ((std::size_t) n, 0.0f);
        dual::tone (x, sr, 0.05, 0.020, 0.9);
        dual::tone (x, sr, 0.20, 0.045, 0.9);
        dual::tone (x, sr, 0.40, 0.060, 0.9);
        dual::tone (x, sr, 0.60, 0.200, 0.9);
        dual::tone (x, sr, 0.80, 0.500, 0.6);
        unsigned long long seed = 7;
        for (int i = (int) (1.6 * sr); i < (int) (2.2 * sr); ++i)
        {
            seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            x[(std::size_t) i] = 0.9f * ((float) ((seed >> 40) & 0xffff) / 32768.0f - 1.0f);
        }
        dual::tone (x, sr, 2.3, 0.4, 0.9);
        const std::vector<dual::Segment> schedules[] {
            { { n, true } },
            { { (int) (0.70 * sr), true }, { (int) (0.90 * sr), false }, { n - (int) (1.60 * sr), true } },
            { { (int) (0.32 * sr), false }, { (int) (0.61 * sr), true }, { n - (int) (0.93 * sr), false } },
            // switched on 30 ms before a held tone ends, twice: the window counts from the switch
            { { (int) (1.27 * sr), false }, { (int) (0.90 * sr), true }, { (int) (0.50 * sr), false }, { n - (int) (2.67 * sr), true } },
        };
        std::size_t worst = 0;
        bool allRan = true, slowWorked = false;
        for (std::size_t si = 0; si < std::size (schedules); ++si)
            for (int block : { 1, 7, 480, 4096 })
            {
                const dual::Run r = dual::render (x, sr, -6.0, 20.0, 180.0, schedules[si], block);
                allRan = allRan && r.ok;
                if (! r.ok) continue;
                worst = std::max (worst, dual::differ (r.gr, dual::reference (r, sr, 20.0, 180.0, schedules[si], block)));
                if (si == 0 && block == 4096)
                {
                    const dual::Run single = dual::render (x, sr, -6.0, 20.0, 180.0, { { n, false } }, block);
                    slowWorked = single.ok && dual::differ (single.out, r.out) > 0;
                }
            }
        test::ok (allRan, "PRECONDITION: every schedule rendered at every block size");
        test::ok (slowWorked, "PRECONDITION: the slow envelope changes the output of this fixture");
        test::ok (worst == 0, "the gain-reduction tap is the reference's, bit for bit, on 4 schedules x 4 block sizes ("
                              + std::to_string (worst) + " samples differ at worst)");
    }

    // KNOWN ANSWERS. After the tone stops, the required reduction is 0 dB and the applied one decays by exactly one
    // coefficient per oversampled sample: the slow one after a stretch longer than the window, the fast one after a
    // shorter stretch — the window counted in the reduction REQUIRED, not in the tone, and from the switch that turned
    // the dual release on: switched on 30 ms before a 200 ms tone ends, the fast one.
    test::group ("dual release: after a held reduction the slow release, after a short one the fast release");
    {
        const double C = -6.0, fast = 20.0, slow = 180.0;
        const int F = 4;
        const float cF = dual::coefFor (fast, sr, F), cS = dual::coefFor (slow, sr, F);
        const int W = limiter::TruePeakLimiter::slowWindowSamplesFor (sr, F);
        test::ok (W == 2400 * F, "the window is 50 ms at 48 kHz: 2400 baseband samples x 4 = 9600 oversampled");
        struct Case { double seconds, onAt; };
        for (const Case cs : { Case { 0.200, 0.0 }, Case { 0.030, 0.0 }, Case { 0.200, 0.170 } })
        {
            const double seconds = cs.seconds;
            const int n = (int) (sr * (seconds + 0.5)), on = (int) (sr * cs.onAt);
            std::vector<float> x ((std::size_t) n, 0.0f);
            dual::tone (x, sr, 0.0, seconds, 0.9);
            const dual::Run r = on > 0 ? dual::render (x, sr, C, fast, slow, { { on, false }, { n - on, true } }, 512)
                                       : dual::render (x, sr, C, fast, slow, { { n, true } }, 512);
            if (! test::run (r.ok)) continue;
            // the last oversampled sample whose required reduction is below 0 dB, recomputed from the peaks
            const std::vector<float> smax = dual::windowMax (r.peak, r.look * F + 1);
            std::size_t last = 0, run = 0, longest = 0;
            for (std::size_t i = 0; i < smax.size(); ++i)
            {
                const bool reducing = r.ceiling - core::gainToDb ((double) smax[i]) < 0.0;
                if (reducing) { last = i; ++run; longest = std::max (longest, run); } else run = 0;
            }
            const std::size_t t1 = last + 1 + 480, t2 = t1 + 1920;           // 2.5 ms and 12.5 ms into the release
            const double ratio = (double) r.gr[t2] / (double) r.gr[t1];
            const bool held = longest >= (std::size_t) W && cs.onAt == 0.0;
            const double want = std::pow ((double) (held ? cS : cF), (double) (t2 - t1));
            const double other = std::pow ((double) (held ? cF : cS), (double) (t2 - t1));
            const std::string at = std::to_string ((int) (seconds * 1000.0)) + " ms of tone"
                                 + (cs.onAt > 0.0 ? ", dual release on 30 ms before its end" : "") + ": the reduction held for "
                                 + std::to_string ((double) longest / (sr * F) * 1000.0) + " ms";
            test::ok ((longest >= (std::size_t) W) == (seconds > 0.1), "PRECONDITION: " + at);
            test::ok (r.gr[t1] < -0.5f, "PRECONDITION: the release starts deep (" + std::to_string (r.gr[t1]) + " dB)");
            std::printf ("      %s: ratio %.9f, slow^n %.9f, fast^n %.9f\n", at.c_str(), ratio, (double) std::pow ((double) cS, (double) (t2 - t1)),
                         (double) std::pow ((double) cF, (double) (t2 - t1)));
            test::ok (std::fabs (ratio / want - 1.0) < 1.0e-4 && std::fabs (ratio / other - 1.0) > 1.0e-2,
                      at + " — over 1920 oversampled samples the reduction falls by " + std::to_string (ratio) + ", the "
                      + (held ? "slow" : "fast") + " coefficient's " + std::to_string (want) + ", not the other's " + std::to_string (other));
        }
    }

    // THE DEFAULT IS THE SINGLE RELEASE, whatever the slow release holds; and the getters.
    test::group ("dual release: off is the single release whatever slowReleaseMs holds; the readbacks");
    {
        const int n = (int) (sr * 1.0);
        std::vector<float> x ((std::size_t) n, 0.0f);
        dual::tone (x, sr, 0.1, 0.4, 0.9);
        const dual::Run a = dual::render (x, sr, -6.0, 20.0, 180.0, { { n, false } }, 1024);
        const dual::Run b = dual::render (x, sr, -6.0, 20.0, 3.0,   { { n, false } }, 1024);
        test::ok (a.ok && b.ok && dual::differ (a.out, b.out) == 0 && dual::differ (a.gr, b.gr) == 0,
                  "with dualRelease off, a slow release of 180 ms and of 3 ms render the same bits");
        test::ok (limiter::TruePeakLimiterParams {}.dualRelease == false, "dualRelease is off by default");

        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, 64, 1, {});
        limiter::TruePeakLimiterParams p;
        lim.setParams (p);
        test::ok (lim.effectiveSlowReleaseMs() == 0.0, "off: effectiveSlowReleaseMs() is 0");
        p.dualRelease = true; p.slowReleaseMs = 180.0; lim.setParams (p);
        test::ok (std::fabs (lim.effectiveSlowReleaseMs() - 180.0) < 0.5, "on: the slow release it runs, ~180 ms ("
                  + std::to_string (lim.effectiveSlowReleaseMs()) + ")");
        p.slowReleaseMs = 0.0; lim.setParams (p);
        test::approx (lim.effectiveSlowReleaseMs(), 8.0 / sr * 1000.0, 1.0e-3, "on: floored at 8 baseband samples, like the fast release");
        p.slowReleaseMs = std::numeric_limits<double>::quiet_NaN(); lim.setParams (p);
        test::ok (std::fabs (lim.effectiveSlowReleaseMs() - 200.0) < 0.5, "a NaN slow release falls back to the default 200 ms ("
                  + std::to_string (lim.effectiveSlowReleaseMs()) + ")");

        limiter::TruePeakLimiter::Storage st;
        const bool okSt = limiter::TruePeakLimiter::storageFor (96000.0, 64, 2, { 1.0, 8, 64 }, st);
        test::ok (okSt && st.slowWindow.entries == (std::size_t) limiter::TruePeakLimiter::slowWindowSamplesFor (96000.0, 8)
                  && st.slowWindow.entries == 4800u * 8u, "the budget holds the window: 50 ms at 96 kHz x 8 = 38 400 entries");
    }

    // RT: switching the dual release on and off and processing with it allocates nothing.
    test::group ("dual release: no allocation in process() or setParams(), and the slow state flushes to exact zero");
    {
        const int n = 512;
        std::vector<float> a ((std::size_t) n), b ((std::size_t) n);
        for (int i = 0; i < n; ++i) a[(std::size_t) i] = b[(std::size_t) i] = (float) (0.9 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        float* ch[2] { a.data(), b.data() };
        limiter::TruePeakLimiter lim; (void) lim.prepare (sr, n, 2, {});
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -6.0; p.dualRelease = true;
        const long long before = alloc::count.load();
        lim.setParams (p);
        for (int k = 0; k < 40; ++k) felitronics::test::run (lim.process (ch, 2, n));
        p.dualRelease = false; lim.setParams (p);
        felitronics::test::run (lim.process (ch, 2, n));
        p.dualRelease = true; lim.setParams (p);
        felitronics::test::run (lim.process (ch, 2, n));
        test::okNoAlloc (alloc::count.load() == before, "zero heap allocations across switches and blocks");

        // A reduction held long enough for the slow envelope, then 8 s of silence: the slow state decays past the flush
        // floor and reads exactly 0 — without the flush it would still be a few 1e-17 dB.
        for (int k = 0; k < 40; ++k)
        {
            for (int i = 0; i < n; ++i) a[(std::size_t) i] = b[(std::size_t) i] = (float) (0.9 * std::sin (2.0 * core::kPi * 1000.0 * (k * n + i) / sr));
            felitronics::test::run (lim.process (ch, 2, n));
        }
        const double deep = lim.gainReductionDb();
        std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
        for (int k = 0; k < (int) (8.0 * sr) / n; ++k) felitronics::test::run (lim.process (ch, 2, n));
        test::ok (deep < -3.0 && lim.gainReductionDb() == 0.0, "after the held reduction (" + std::to_string (deep)
                  + " dB) and 8 s of silence the applied reduction is exactly 0 (" + std::to_string (lim.gainReductionDb()) + ")");
    }

    return test::report();
}
