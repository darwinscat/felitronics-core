// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::mastering — `MasteringChainParams::compressorMix`, parallel compression (P60).
//
// THE THREE GATES, and where each one is proved:
//
//  1. mix = 1 is the chain before the field, BIT FOR BIT. Proved three ways. Here: against a bare
//     `dynamics::Compressor` driven by hand, and through the narrowing rule (a mix that rounds to 1 IS 1).
//     In MasteringOracleTests: the whole chain against every stage driven by hand (its mix-1 path is the
//     one it always had; P60 added the dry lines beside it). And
//     once, outside the suite, against a render built at the base commit — every output sample and every
//     tap of seven scenarios, identical (the P60 report carries the numbers).
//  2. mix = 0 is the stage's input delayed by exactly the compressor's lookahead, BIT FOR BIT, found from
//     the INPUT rather than from the chain's own aligner — so a delay that is wrong in the aligner and in
//     the latency together still fails. Signed zeros are in the fixture on purpose.
//  3. In between, the output is within one float ulp of the exact blend of the two ends (a black-box
//     property), and its bits are the double arithmetic the header states (the formula recomputed, on
//     purpose — see the test). NOT "correctly rounded": the double sum is rounded to float, and a double
//     sum can sit on a float midpoint the exact value is not on.
//
// And the things that make a dry path go wrong without touching any of the three: a ring that stops
// advancing while nobody listens to it, a reset that forgets it, a bypass that stops blending, a mix change that
// depends on where the caller cut the stream, and a clamp that lets a value through.

#include <felitronics/mastering/MasteringChain.h>
#include <felitronics_test.h>

#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <string>
#include <vector>

// The allocation counter the RT check reads. Counting only — the budgets are pinned in MasteringChainTests.
static std::atomic<long long> g_allocs { 0 };
// No `throw` anywhere in here: the wasm-audio tier builds with -fno-exceptions (it caught the first draft).
void* operator new   (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[] (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new   (std::size_t s, const std::nothrow_t&) noexcept { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[] (std::size_t s, const std::nothrow_t&) noexcept { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
// The OVER-ALIGNED forms too — the EQ engine is built through one, and a counter without them cannot see an
// aligned allocation the mix path might grow (the code-review round). Same allocator as MasteringChainTests.
static void* alignedNew (std::size_t s, std::align_val_t a)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
#if defined(_MSC_VER)
    return _aligned_malloc (s ? s : 1, (std::size_t) a);
#else
    const std::size_t al = (std::size_t) a < sizeof (void*) ? sizeof (void*) : (std::size_t) a;
    void* p = nullptr;
    return posix_memalign (&p, al, s ? s : 1) == 0 ? p : nullptr;
#endif
}
static void alignedFree (void* p) noexcept
{
#if defined(_MSC_VER)
    _aligned_free (p);
#else
    std::free (p);
#endif
}
void* operator new   (std::size_t s, std::align_val_t a) { return alignedNew (s, a); }
void* operator new[] (std::size_t s, std::align_val_t a) { return alignedNew (s, a); }
void  operator delete   (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::okNoAlloc;
using felitronics::test::group;

namespace
{
constexpr double kPi = 3.14159265358979323846;
using Buf = std::vector<std::vector<float>>;

std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }

Buf silence (int nch, int n) { return Buf ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f)); }

std::vector<float*> planes (Buf& b, int off = 0)
{
    std::vector<float*> p;
    p.reserve (b.size());
    for (auto& v : b) p.push_back (v.data() + off);
    return p;
}

long long countBitDiffs (const Buf& a, const Buf& b, int from = 0) noexcept
{
    if (a.size() != b.size()) return -1;
    long long d = 0;
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        if (a[c].size() != b[c].size()) return -1;
        for (std::size_t i = (std::size_t) from; i < a[c].size(); ++i) d += bits (a[c][i]) != bits (b[c][i]);
    }
    return d;
}

// Loud, dynamic, different per channel, never silent for long — so the compressor really works, the dry
// and the compressed paths really differ, and a stale or misaligned dry sample cannot hide behind a zero.
// EXACT -0.0f is planted every 331 samples, so a dry path that re-signs a zero on the way (a `+ 0.0f`, a
// gate, a different copy) shows in the bit comparisons.
Buf programme (int nch, int n, unsigned seed = 4242u)
{
    Buf x = silence (nch, n);
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / 48000.0;
        const double env = 0.2 + 0.8 * std::fabs (std::sin (2.0 * kPi * 2.3 * t));
        for (int c = 0; c < nch; ++c)
        {
            const double f = 97.0 + 31.0 * c;
            float v = (float) (env * (0.6 * std::sin (2.0 * kPi * f * t) + 0.25 * std::sin (2.0 * kPi * 1733.0 * t + c)))
                    + 0.1f * u (rng);
            if ((i + 17 * c) % 331 == 0) v = -0.0f;
            x[(std::size_t) c][(std::size_t) i] = v;
        }
    }
    return x;
}

int countNegZeros (const Buf& x) noexcept
{
    int k = 0;
    for (const auto& v : x) for (float s : v) k += bits (s) == 0x80000000u;
    return k;
}

// Compressor only, so the chain's output IS the compressor stage's output, delayed by the quantum.
mastering::MasteringChainConfig compOnly (int K = 256, double lookMs = 5.0)
{
    mastering::MasteringChainConfig c;
    c.internalBlock = K;
    c.eq = c.monoBass = c.clipper = c.limiter = c.dither = false;
    c.compressor = true;
    c.compressorLookaheadMs = lookMs;
    return c;
}

mastering::MasteringChainParams heavy (double mix = 1.0)
{
    mastering::MasteringChainParams p;
    p.compressor.thresholdDb = -30.0;
    p.compressor.ratio       = 6.0;
    p.compressor.attackMs    = 2.0;
    p.compressor.releaseMs   = 80.0;
    p.compressor.makeupDb    = 6.0;
    p.compressorMix          = mix;
    return p;
}

// Every stage present and working, the key filter on: the topology where the dry path has the most
// neighbours to be confused with.
mastering::MasteringChainConfig fullConfig (int K = 256)
{
    mastering::MasteringChainConfig c;
    c.internalBlock = K;
    c.monoBass = true;
    c.clipper = true;
    c.sidechainHpfHz = 120.0;
    c.compressorLookaheadMs = 3.0;
    return c;
}

mastering::MasteringChainParams voicing (double mix = 1.0)
{
    mastering::MasteringChainParams p = heavy (mix);
    p.inputGainDb = 1.5;
    p.preLimiterGainDb = 2.0;
    p.eqBands[0].on = true;
    p.eqBands[0].type = eq::FilterType::Bell;
    p.eqBands[0].lane (eq::Lane::Stereo).on = true;
    p.eqBands[0].lane (eq::Lane::Stereo).freq = 2800.0;
    p.eqBands[0].lane (eq::Lane::Stereo).Q = 0.9;
    p.eqBands[0].lane (eq::Lane::Stereo).gainDb = 3.0;
    p.monoBass = { true, 110.0f, 0.0f };
    p.clipper.driveDb = 4.0f;
    p.limiter.ceilingDbTp = -1.0;
    p.dither.bits = 24;
    p.dither.autoBlank = true;
    return p;
}

// One render. `changes` is a list of (stream position, params) applied by SPLITTING the caller's calls at
// that position, so a change lands at the same stream position under every partition.
struct Change { int at; mastering::MasteringChainParams p; };
struct Render { Buf y; std::vector<float> gr; Buf pre; bool ok = false; int latency = 0; };

Render render (const mastering::MasteringChainConfig& cfg, const mastering::MasteringChainParams& p0, const Buf& x,
               const std::vector<int>& part, const std::vector<Change>& changes = {}, bool taps = false)
{
    Render r;
    const int nch = (int) x.size(), n = (int) x[0].size();
    mastering::MasteringChain chain;
    r.y = x;
    if (! chain.prepare (48000.0, nch, cfg)) return r;
    chain.setParams (p0);
    r.latency = chain.latencySamples();
    const int K = chain.internalBlock();
    std::vector<float> gr;
    Buf pre;
    std::vector<float*> prePtr;
    if (taps)
    {
        const int cap = n + K;
        gr.assign ((std::size_t) cap, 0.0f);
        pre = silence (nch, cap);
        prePtr = planes (pre);
    }
    std::size_t nextChange = 0;
    int off = 0, framesSoFar = 0;
    std::size_t pi = 0;
    while (off < n)
    {
        int take = std::min (part[pi % part.size()], n - off);
        if (nextChange < changes.size() && changes[nextChange].at > off) take = std::min (take, changes[nextChange].at - off);
        if (nextChange < changes.size() && changes[nextChange].at == off) { chain.setParams (changes[nextChange].p); ++nextChange; continue; }
        auto io = planes (r.y, off);
        bool accepted = false;
        if (taps)
        {
            mastering::MasteringChainTaps t;
            std::vector<float*> preOff;
            for (auto* q : prePtr) preOff.push_back (q + framesSoFar);
            t.compressorGrDb = gr.data() + framesSoFar;
            t.preLimiter = preOff.data();
            t.frameCapacity = n + K - framesSoFar;
            accepted = chain.process (io.data(), nch, take, t);
            framesSoFar += t.framesWritten;
        }
        else accepted = chain.process (io.data(), nch, take);
        if (! accepted) return r;
        off += take;
        ++pi;
    }
    if (taps)
    {
        gr.resize ((std::size_t) framesSoFar);
        for (auto& v : pre) v.resize ((std::size_t) framesSoFar);
        r.gr = std::move (gr);
        r.pre = std::move (pre);
    }
    r.ok = true;
    return r;
}

std::vector<int> fixedPart (int b) { return { b }; }

std::vector<int> primesPart()
{
    return { 1, 2, 3, 5, 7, 11, 13, 17, 23, 31, 61, 127, 251, 509, 1021, 2039 };
}
} // namespace

//==============================================================================
static void testDefaultsAndReadback()
{
    group ("the field — default, read-back, absent stage");

    ok (core::exactlyEqual (mastering::MasteringChainParams {}.compressorMix, 1.0), "the default mix is 1");

    mastering::MasteringChain unprepared;
    ok (core::exactlyEqual (unprepared.resolved().compressorMix, 0.0), "an unprepared chain reads 0, like every resolved field");

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, 2, compOnly()), "PRECONDITION: prepare");
    ok (core::exactlyEqual (chain.resolved().compressorMix, 1.0), "a prepared chain applies the default of 1 before any process()");
    chain.setParams (heavy (0.25));
    Buf x = programme (2, 512);
    { auto q = planes (x); felitronics::test::run (chain.process (q.data(), 2, 512)); }
    ok (core::exactlyEqual (chain.resolved().compressorMix, 0.25), "resolved() reads the applied mix once a quantum has run");

    auto noComp = compOnly();
    noComp.compressor = false;
    mastering::MasteringChain bare;
    ok (bare.prepare (48000.0, 2, noComp), "PRECONDITION: a chain with no compressor");
    bare.setParams (heavy (0.25));
    ok (core::exactlyEqual (bare.resolved().compressorMix, 0.0), "no compressor, no mix: 0, as every absent stage reads");

    // The latency is a property of the topology, never of a per-block parameter.
    int lat[3] {};
    const double mixes[3] { 0.0, 0.5, 1.0 };
    for (int k = 0; k < 3; ++k)
    {
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, 2, fullConfig()), "PRECONDITION: prepare for the latency check");
        c.setParams (voicing (mixes[k]));
        Buf q0 = programme (2, 256);                        // the mix is APPLIED at a quantum, not at setParams()
        { auto q = planes (q0); felitronics::test::run (c.process (q.data(), 2, 256)); }
        ok (core::exactlyEqual (c.resolved().compressorMix, mixes[k]), "PRECONDITION: the mix has been applied");
        lat[k] = c.latencySamples();
    }
    ok (lat[0] == lat[1] && lat[1] == lat[2], "latencySamples() does not move with the mix (" + std::to_string (lat[0]) + ")");
}

//==============================================================================
// GATE 1, locally: the compressor stage at mix 1 is a bare `dynamics::Compressor`, bit for bit.
static void testMixOneIsTheCompressor()
{
    group ("gate 1 — mix = 1 is the compressor alone, bit for bit");

    const int nch = 2, n = 48000, K = 256;
    const auto cfg = compOnly (K, 5.0);
    const Buf x = programme (nch, n);
    ok (countNegZeros (x) > 100, "PRECONDITION: the fixture carries exact -0.0f samples (" + std::to_string (countNegZeros (x)) + ")");

    const Render viaChain = render (cfg, heavy (1.0), x, fixedPart (1021));
    ok (viaChain.ok, "PRECONDITION: the chain renders");

    dynamics::Compressor comp;
    ok (comp.prepare (48000.0, K, nch, std::max (cfg.compressorLookaheadMs, 1.0)), "hand: compressor prepare");
    { auto cp = heavy().compressor; cp.lookaheadMs = cfg.compressorLookaheadMs; comp.setParams (cp); }
    Buf byHand = x;
    for (int off = 0; off + K <= n; off += K)
    {
        auto q = planes (byHand, off);
        felitronics::test::run (comp.process (q.data(), nch, K));
    }
    long long bad = 0, compressed = 0;
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i + K < n; ++i)
        {
            bad += bits (viaChain.y[(std::size_t) c][(std::size_t) (i + K)]) != bits (byHand[(std::size_t) c][(std::size_t) i]);
            compressed += std::fabs (byHand[(std::size_t) c][(std::size_t) i]) > 1.0e-6f
                          && ! core::exactlyEqual (byHand[(std::size_t) c][(std::size_t) i],
                                                   i >= comp.latencySamples() ? x[(std::size_t) c][(std::size_t) (i - comp.latencySamples())] : 0.0f);
        }
    ok (compressed > n / 2, "PRECONDITION: the compressor really compresses (" + std::to_string (compressed) + " samples moved)");
    ok (bad == 0, "the chain at mix 1 is the compressor driven by hand, bit for bit (" + std::to_string (bad) + " differ)");

    // THE NARROWING RULE. The mix is applied as a float, so every double that rounds to 1.0f IS 1, and reads
    // back as 1. (Which path it then takes does not show in the bits, and this does not pretend to pin it:
    // the compressed sample shares the dry sample's sign, so the blend at exactly 1 returns it unchanged —
    // see the header, and the stand's note on the two equivalent mutants.)
    for (const double m : { 1.0 - 1.0e-9, std::nextafter (1.0, 0.0), 1.0 - 1.0e-8 })
    {
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, cfg), "PRECONDITION: prepare");
        const Render r = render (cfg, heavy (m), x, fixedPart (1021));
        ok (r.ok && countBitDiffs (r.y, viaChain.y) == 0,
            "a mix of " + std::to_string (m) + " narrows to 1.0f and renders as mix 1 exactly ("
            + std::to_string (countBitDiffs (r.y, viaChain.y)) + " differ)");
        c.setParams (heavy (m));
        Buf s = silence (nch, K);
        { auto q = planes (s); felitronics::test::run (c.process (q.data(), nch, K)); }
        ok (core::exactlyEqual (c.resolved().compressorMix, 1.0), "...and resolved() reads the 1 it became");
    }
}

//==============================================================================
// GATE 2: the input, found in the INPUT.
static void testMixZeroIsTheDelayedInput()
{
    group ("gate 2 — mix = 0 is the stage's input delayed by the lookahead, bit for bit");

    const int n = 40000;
    struct Row { int nch; int K; double lookMs; };
    const Row rows[] = { { 2, 256, 5.0 }, { 1, 256, 0.0 }, { 1, 8, 1.0 }, { 2, 64, 20.0 }, { 16, 256, 1.0 }, { 3, 4096, 250.0 } };
    for (const Row& row : rows)
    {
        const auto cfg = compOnly (row.K, row.lookMs);
        const Buf x = programme (row.nch, n, 99u + (unsigned) row.nch);
        const Render r0 = render (cfg, heavy (0.0), x, primesPart());
        const Render r1 = render (cfg, heavy (1.0), x, primesPart());
        const std::string tag = "nch=" + std::to_string (row.nch) + " K=" + std::to_string (row.K)
                              + " lookahead=" + std::to_string ((int) row.lookMs) + " ms";
        ok (r0.ok && r1.ok, "PRECONDITION: both render — " + tag);
        const int D = r0.latency;
        ok (D == row.K + (int) std::lround (row.lookMs * 48.0), "PRECONDITION: the latency is K + the lookahead — " + tag);

        long long bad = 0, primed = 0, negZerosKept = 0;
        for (int c = 0; c < row.nch; ++c)
            for (int i = 0; i < n; ++i)
            {
                const float y = r0.y[(std::size_t) c][(std::size_t) i];
                if (i < D) primed += bits (y) != 0u;
                else
                {
                    const float want = x[(std::size_t) c][(std::size_t) (i - D)];
                    bad += bits (y) != bits (want);
                    negZerosKept += bits (want) == 0x80000000u && bits (y) == 0x80000000u;
                }
            }
        ok (bad == 0, "out[i] == in[i - latency], sign of zero included (" + std::to_string (bad) + " differ) — " + tag);
        ok (primed == 0, "the first `latency` samples are +0.0f exactly — " + tag);
        ok (negZerosKept > 0, "PRECONDITION: -0.0f samples crossed the dry path and kept their sign — " + tag);
        ok (countBitDiffs (r0.y, r1.y) > n / 4, "PRECONDITION: mix 1 is NOT the delayed input here, so the check above can fail — " + tag);
    }

    // THE WHOLE CHAIN: mix 0 renders exactly what a bypassed compressor renders, from the first sample —
    // every downstream stage sees the same input. Not across a bypass TOGGLE: entering bypass lets the
    // compressor's gain reduction release over its release time, and mix 0 has no such glide.
    {
        const int nch = 2;
        const Buf x = programme (nch, n, 7u);
        auto pBypass = voicing (1.0);
        pBypass.bypassCompressor = true;
        const Render a = render (fullConfig(), voicing (0.0), x, fixedPart (777), {}, true);
        const Render b = render (fullConfig(), pBypass, x, fixedPart (777), {}, true);
        const Render c = render (fullConfig(), voicing (1.0), x, fixedPart (777), {}, true);
        ok (a.ok && b.ok && c.ok, "PRECONDITION: three full-chain renders with taps");
        ok (countBitDiffs (a.y, b.y) == 0, "full chain: mix 0 is the bypassed compressor, output bit for bit ("
                                           + std::to_string (countBitDiffs (a.y, b.y)) + " differ)");
        ok (countBitDiffs (a.pre, b.pre) == 0, "...and so is the pre-limiter tap");
        ok (countBitDiffs (a.y, c.y) > 1000, "PRECONDITION: and neither is mix 1 (" + std::to_string (countBitDiffs (a.y, c.y)) + " differ)");

        // THE GAIN-REDUCTION TAP DOES NOT MOVE WITH THE MIX. It is the compressed path's gain reduction —
        // at mix 0 it describes a signal nobody hears — and the header says so. Bypass traces 0 dB.
        long long grDiff = 0, grActive = 0, bypassNonZero = 0;
        for (std::size_t i = 0; i < a.gr.size(); ++i)
        {
            grDiff += bits (a.gr[i]) != bits (c.gr[i]);
            grActive += a.gr[i] < -0.5f;
            bypassNonZero += ! core::exactlyEqual (b.gr[i], 0.0f);
        }
        ok (a.gr.size() == c.gr.size() && grDiff == 0, "the compressor's GR trace is identical at mix 0 and mix 1");
        ok (grActive > 1000, "PRECONDITION: and it shows real gain reduction (" + std::to_string (grActive) + " frames past 0.5 dB)");
        ok (bypassNonZero == 0, "while a bypassed compressor traces exactly 0 dB");
    }
}

//==============================================================================
// GATE 3: the correctly rounded convex combination of the two ends.
static void testMixBetween()
{
    group ("gate 3 — between the ends, the blend of the two ends, in the stated arithmetic");

    const int nch = 2, n = 30000;
    const auto cfg = compOnly (256, 5.0);
    const Buf x = programme (nch, n, 31337u);
    const Render dry = render (cfg, heavy (0.0), x, fixedPart (4096));
    const Render wet = render (cfg, heavy (1.0), x, fixedPart (4096));
    ok (dry.ok && wet.ok, "PRECONDITION: both ends render");
    const int D = dry.latency;

    // TWO CHECKS, and neither pretends to be the other.
    //  (a) BLACK-BOX: the mid render is within ONE float ulp of the exact `(1-m)*A + m*B`, where A and B are
    //      the chain's own mix-0 and mix-1 renders. Weights, alignment and where the blend sits all have to
    //      be right for this; the reference carries each product's error term, so it does not share the
    //      chain's rounding. One ulp and not half, on purpose: the chain rounds the DOUBLE sum to float,
    //      and a double sum can land exactly on a float midpoint the exact value is not on — the
    //      code-review round built one at m = 0x3c800006, where the nearest float to the exact blend is
    //      0x3f21c001 and the chain's is 0x3f21c000. "Correctly rounded" would be a false claim.
    //  (b) THE ARITHMETIC ITSELF: the bits are the double form the header states, products in separate
    //      statements. This IS the formula recomputed, and that is the point of it — the form is the law-10
    //      contract (every row this tree builds agrees on it), and the float spelling of the same blend
    //      passes (a) and fails this (the stand measured it 0.93 ulp out).
    // Mixes on both sides of 2^-6, below which the first double product can round.
    for (const double m : { 0.5, 0.25, 0.7, 0.9, 0.999, 0.1, 0.03, 0.015625, 0.01, 1.0e-3, 1.0e-5 })
    {
        const Render mid = render (cfg, heavy (m), x, fixedPart (4096));
        ok (mid.ok, "PRECONDITION: the mid mix renders");
        const float mf = (float) m;
        const double a = 1.0 - (double) mf, b = (double) mf;
        long long inside = 0, samples = 0, sameBits = 0, bothEnds = 0, pastHalf = 0;
        double worstUlps = 0.0;
        for (int c = 0; c < nch; ++c)
            for (int i = D; i < n; ++i)
            {
                const float  df = dry.y[(std::size_t) c][(std::size_t) i];
                const float  wf = wet.y[(std::size_t) c][(std::size_t) i];
                const float  y  = mid.y[(std::size_t) c][(std::size_t) i];
                const double d = df, w = wf;
                const double pd = a * d, epd = std::fma (a, d, -pd);
                const double pw = b * w, epw = std::fma (b, w, -pw);
                const double exact = pd + pw + (epd + epw);
                const float  ay = std::fabs (y);
                const double ulp = std::max ((double) (std::nextafter (ay, std::numeric_limits<float>::infinity()) - ay),
                                             (double) (ay - std::nextafter (ay, 0.0f)));
                const double err = std::fabs ((double) y - exact);
                // SPELLED POSITIVELY, so a NaN output is outside rather than silently at zero ulps.
                const bool within = std::isfinite (y) && ulp > 0.0 && err <= ulp * (1.0 + 1.0e-6);
                inside += within;
                ++samples;
                if (within) worstUlps = std::max (worstUlps, err / ulp);
                pastHalf += within && err > ulp * (0.5 + 1.0e-6);
                const double spd = a * d;
                const double spw = b * w;
                sameBits += bits ((float) (spd + spw)) == bits (y);
                bothEnds += std::fabs (d - w) > 1.0e-4;
            }
        ok (inside == samples, "mix " + std::to_string (m) + " (a): every sample finite and within one ulp of (1-m)*A + m*B (worst "
                               + std::to_string (worstUlps) + " ulp, " + std::to_string (samples - inside) + " outside, "
                               + std::to_string (pastHalf) + " past half an ulp)");
        ok (sameBits == samples, "mix " + std::to_string (m) + " (b): every sample is the stated double arithmetic, bit for bit ("
                                 + std::to_string (samples - sameBits) + " differ)");
        ok (bothEnds > n / 2, "PRECONDITION: the two ends differ on most samples, so the blend is not trivially either");
    }
}

//==============================================================================
static void testBypassAndMix()
{
    group ("the mix keeps blending through a bypass — a steady bypass untouched, an engaging one without a drop");

    // A STEADY BYPASS IS THE SAME AT EVERY MIX. Bypassed from the first sample the compressed path IS the
    // dry path, and the double blend of a sample with itself returns it — so this holds whether or not
    // the chain blends while bypassed, and it is here as the regression, not as the discriminator.
    const int nch = 2, n = 30000;
    const Buf x = programme (nch, n, 5u);
    auto p1 = voicing (1.0);
    p1.bypassCompressor = true;
    const Render ref = render (fullConfig(), p1, x, fixedPart (512));
    ok (ref.ok, "PRECONDITION: the bypassed reference renders");
    for (const double m : { 0.0, 0.3, 0.999, std::numeric_limits<double>::quiet_NaN() })
    {
        auto p = p1;
        p.compressorMix = m;
        const Render r = render (fullConfig(), p, x, fixedPart (512));
        ok (r.ok && countBitDiffs (r.y, ref.y) == 0,
            "a compressor bypassed from the start renders the same at mix " + std::to_string (m) + " as at mix 1");
    }

    // THE DISCRIMINATOR: a bypass ENGAGED mid-stream. The warm bypass releases the gain reduction over the
    // release time, so for that long the compressed path is still compressed. The first version skipped
    // the blend there and the output fell to the compressed path in one sample — at mix 0, where the
    // output had been the dry signal, a 20 dB drop (the diverse-testing round). The contract now:
    //  (a) at mix 0 a bypass toggle changes NOTHING, bit for bit;
    //  (b) at a mid mix, before, during and after the bypass, the output is the stated blend of the mix-0
    //      render and the mix-1 render of the SAME toggle schedule — the release blended, not dropped into.
    {
        const auto cfg = compOnly (256, 5.0);
        const int K = 256, P = 60 * K, Q = 120 * K, total = 180 * K;
        const Buf y = programme (2, total, 77u);
        auto on = [] (double m) { auto p = heavy (m); p.bypassCompressor = true; return p; };
        const Render zeroToggled = render (cfg, heavy (0.0), y, fixedPart (1000), { { P, on (0.0) }, { Q, heavy (0.0) } });
        const Render zeroPlain   = render (cfg, heavy (0.0), y, fixedPart (1000));
        const Render oneToggled  = render (cfg, heavy (1.0), y, fixedPart (1000), { { P, on (1.0) }, { Q, heavy (1.0) } });
        const Render midToggled  = render (cfg, heavy (0.3), y, fixedPart (1000), { { P, on (0.3) }, { Q, heavy (0.3) } });
        ok (zeroToggled.ok && zeroPlain.ok && oneToggled.ok && midToggled.ok, "PRECONDITION: the toggle renders");
        const int D = midToggled.latency;

        long long glide = 0;
        for (int c = 0; c < 2; ++c)
            for (int i = P + K; i < Q + K; ++i)
                glide += bits (oneToggled.y[(std::size_t) c][(std::size_t) i]) != bits (y[(std::size_t) c][(std::size_t) (i - D)]);
        ok (glide > 100, "PRECONDITION: engaging the bypass releases, so the compressed path is not yet the input ("
                         + std::to_string (glide) + " samples)");

        const long long zeroDiff = countBitDiffs (zeroToggled.y, zeroPlain.y);
        ok (zeroDiff == 0, "(a) at mix 0, engaging and lifting the bypass changes nothing (" + std::to_string (zeroDiff) + " differ)");

        const float mf = 0.3f;
        const double a = 1.0 - (double) mf, b = (double) mf;
        long long off = 0;
        for (int c = 0; c < 2; ++c)
            for (int i = D; i < total; ++i)
            {
                const double pd = a * (double) zeroPlain.y[(std::size_t) c][(std::size_t) i];
                const double pw = b * (double) oneToggled.y[(std::size_t) c][(std::size_t) i];
                off += bits ((float) (pd + pw)) != bits (midToggled.y[(std::size_t) c][(std::size_t) i]);
            }
        ok (off == 0, "(b) at mix 0.3 the toggled render is the blend of the mix-0 render and the toggled mix-1 render, "
                      "bypass and release included (" + std::to_string (off) + " differ)");
    }
}

//==============================================================================
// THE RING MUST NOT STOP WHILE NOBODY LISTENS. After any change, the chain is a chain that had the new
// value all along — which is only true if the dry ring and the compressor both kept running.
static void testDryPathStaysWarm()
{
    group ("the dry path stays warm across mix and bypass changes");

    const int nch = 2, n = 48000, K = 256;
    const auto cfg = compOnly (K, 5.0);
    const Buf x = programme (nch, n, 2024u);
    const int P = 40 * K;                                   // a quantum boundary, so the change applies at P

    struct Case { const char* name; mastering::MasteringChainParams before, after; };
    auto bypassed = heavy (0.0);
    bypassed.bypassCompressor = true;
    const Case cases[] = {
        { "mix 1 -> 0",           heavy (1.0), heavy (0.0) },
        { "mix 1 -> 0.5",         heavy (1.0), heavy (0.5) },
        { "mix 0 -> 1",           heavy (0.0), heavy (1.0) },
        { "mix 0.2 -> 0.8",       heavy (0.2), heavy (0.8) },
        { "bypassed -> mix 0",    bypassed,    heavy (0.0) },
    };
    for (const Case& cs : cases)
    {
        const Render changed = render (cfg, cs.before, x, fixedPart (1000), { { P, cs.after } });
        const Render always  = render (cfg, cs.after,  x, fixedPart (1000));
        ok (changed.ok && always.ok, std::string ("PRECONDITION: renders — ") + cs.name);
        // The quantum holding input [P, P+K) leaves the FIFO at output [P+K, P+2K).
        const long long after = countBitDiffs (changed.y, always.y, P + K);
        ok (after == 0, std::string ("from the first changed quantum on, identical to a chain that always had it — ")
                        + cs.name + " (" + std::to_string (after) + " differ)");
        ok (countBitDiffs (changed.y, always.y) > 0 || std::string (cs.name) == "bypassed -> mix 0",
            std::string ("PRECONDITION: and different before it — ") + cs.name);
    }
}

//==============================================================================
static void testBlockInvariance()
{
    group ("block invariance — a mid mix, and a mix change, at any partition");

    const int nch = 2, n = 60000;
    const Buf x = programme (nch, n, 8080u);
    const std::vector<Change> moves { { 20011, voicing (0.35) }, { 41000, voicing (0.0) }, { 50111, voicing (0.8) } };
    const Render ref = render (fullConfig(), voicing (0.6), x, fixedPart (n), moves);
    ok (ref.ok, "PRECONDITION: the whole-programme render");
    struct Part { const char* name; std::vector<int> p; };
    const Part parts[] = {
        { "1 sample", fixedPart (1) }, { "primes", primesPart() }, { "255", fixedPart (255) },
        { "257", fixedPart (257) }, { "4099", fixedPart (4099) },
    };
    for (const Part& part : parts)
    {
        const Render r = render (fullConfig(), voicing (0.6), x, part.p, moves);
        ok (r.ok && countBitDiffs (r.y, ref.y) == 0,
            std::string ("bit-identical at ") + part.name + " (" + std::to_string (countBitDiffs (r.y, ref.y)) + " differ)");
    }
}

//==============================================================================
static void testClampAndNonFinite()
{
    group ("the clamp and the non-finite rule");

    const int nch = 1, n = 12000;
    const auto cfg = compOnly (256, 1.0);
    const Buf x = programme (nch, n, 11u);
    const Render one  = render (cfg, heavy (1.0), x, fixedPart (256));
    const Render zero = render (cfg, heavy (0.0), x, fixedPart (256));
    ok (one.ok && zero.ok && countBitDiffs (one.y, zero.y) > 0, "PRECONDITION: the two ends differ");

    const double inf = std::numeric_limits<double>::infinity();
    struct Row { double m; bool isOne; const char* name; };
    const Row rows[] = {
        { std::numeric_limits<double>::quiet_NaN(), true, "NaN" }, { inf, true, "+inf" }, { -inf, true, "-inf" },
        { 7.0, true, "7" }, { 1.0000001, true, "just above 1" },
        { -3.0, false, "-3" }, { -0.0, false, "-0.0" }, { 1.0e-300, false, "1e-300 (narrows to 0.0f)" },
    };
    for (const Row& row : rows)
    {
        const Render r = render (cfg, heavy (row.m), x, fixedPart (256));
        const long long d = countBitDiffs (r.y, row.isOne ? one.y : zero.y);
        ok (r.ok && d == 0, std::string ("mix ") + row.name + " renders as mix " + (row.isOne ? "1" : "0")
                            + " (" + std::to_string (d) + " differ)");
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, cfg), "PRECONDITION: prepare");
        c.setParams (heavy (row.m));
        Buf s = silence (nch, 256);
        { auto q = planes (s); felitronics::test::run (c.process (q.data(), nch, 256)); }
        const double got = c.resolved().compressorMix;
        ok (std::isfinite (got) && core::exactlyEqual (got, row.isOne ? 1.0 : 0.0) && ! std::signbit (got),
            std::string ("...and resolved() reads a finite, positively signed ") + (row.isOne ? "1" : "0") + " for " + row.name);
    }
}

//==============================================================================
static void testResetPrepareAndOrder()
{
    group ("reset, re-prepare and configure-before-prepare reach the dry path");

    const int nch = 2, n = 30000;
    const Buf x = programme (nch, n, 606u);

    // reset(): a second render after a reset is a fresh chain's render. At a MID mix, so a dry ring the
    // reset forgot is heard (at mix 1 it is not, and the existing reset check in MasteringChainTests runs
    // at mix 1 — which is why this one exists).
    {
        const auto cfg = compOnly (256, 5.0);
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, cfg), "PRECONDITION: prepare");
        c.setParams (heavy (0.4));
        Buf first = x;
        { auto q = planes (first); felitronics::test::run (c.process (q.data(), nch, n)); }
        c.reset();
        Buf second = x;
        { auto q = planes (second); felitronics::test::run (c.process (q.data(), nch, n)); }
        const Render fresh = render (cfg, heavy (0.4), x, fixedPart (n));
        ok (countBitDiffs (second, fresh.y) == 0, "after reset() a mid-mix render is a fresh chain's, bit for bit ("
                                                  + std::to_string (countBitDiffs (second, fresh.y)) + " differ)");
    }

    // Re-prepare to a LONGER lookahead and a different quantum: the aligner is re-sized with the stage,
    // so mix 0 is still exactly the delayed input — found in the input, not asked of the chain.
    {
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, compOnly (256, 1.0)), "PRECONDITION: a short lookahead first");
        c.setParams (heavy (0.0));
        Buf junk = x;
        { auto q = planes (junk); felitronics::test::run (c.process (q.data(), nch, 5000)); }
        const auto longer = compOnly (128, 30.0);
        ok (c.prepare (48000.0, nch, longer), "PRECONDITION: re-prepared at 30 ms and K = 128");
        Buf y = x;
        { auto q = planes (y); felitronics::test::run (c.process (q.data(), nch, n)); }
        const int D = 128 + 1440;
        ok (c.latencySamples() == D, "PRECONDITION: the new latency is K + 30 ms");
        long long bad = 0;
        for (int ch = 0; ch < nch; ++ch)
            for (int i = 0; i < n; ++i)
                bad += bits (y[(std::size_t) ch][(std::size_t) i]) != bits (i < D ? 0.0f : x[(std::size_t) ch][(std::size_t) (i - D)]);
        ok (bad == 0, "after re-preparing, mix 0 is the input delayed by the NEW latency (" + std::to_string (bad) + " differ)");
    }

    // Configure, THEN prepare — the order a C-ABI facade takes. The mix written first is kept.
    {
        const auto cfg = compOnly (256, 5.0);
        mastering::MasteringChain c;
        c.setParams (heavy (0.0));
        ok (c.prepare (48000.0, nch, cfg), "PRECONDITION: prepare after setParams");
        ok (core::exactlyEqual (c.resolved().compressorMix, 0.0), "the mix written before prepare() is the applied one");
        Buf y = x;
        { auto q = planes (y); felitronics::test::run (c.process (q.data(), nch, n)); }
        const Render ref = render (cfg, heavy (0.0), x, fixedPart (n));
        ok (countBitDiffs (y, ref.y) == 0, "and it renders as prepare-then-configure does");
    }
}

//==============================================================================
static void testNoAllocation()
{
    group ("RT — the blend, the copy and the ring allocate nothing");

    const int nch = 2, n = 20000;
    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, fullConfig()), "PRECONDITION: prepare (allocates, by design)");
    chain.setParams (voicing (0.5));
    Buf x = programme (nch, n, 1u);
    float* px[core::kMaxChannels] {};
    auto at = [&] (int off) { for (int c = 0; c < nch; ++c) px[c] = x[(std::size_t) c].data() + off; };
    const auto p0 = voicing (0.0), p1 = voicing (1.0), pm = voicing (0.37);
    at (0);
    felitronics::test::run (chain.process (px, nch, 1024));           // warm-up, outside the counted region
    chain.reset();

    const long long before = g_allocs.load();
    at (0);    felitronics::test::run (chain.process (px, nch, 5000));
    chain.setParams (p0);
    at (5000); felitronics::test::run (chain.process (px, nch, 5000));
    chain.setParams (p1);
    at (10000); felitronics::test::run (chain.process (px, nch, 5000));
    chain.setParams (pm);
    at (15000); felitronics::test::run (chain.process (px, nch, 5000));
    (void) chain.resolved();
    const long long after = g_allocs.load();
    okNoAlloc (after == before, "process() at mix 0.5, 0, 1 and 0.37 with the changes between allocates nothing ("
                                + std::to_string (after - before) + ")");
}

int main()
{
    std::printf ("felitronics::mastering — compressorMix (parallel compression)\n");
    testDefaultsAndReadback();
    testMixOneIsTheCompressor();
    testMixZeroIsTheDelayedInput();
    testMixBetween();
    testBypassAndMix();
    testDryPathStaysWarm();
    testBlockInvariance();
    testClampAndNonFinite();
    testResetPrepareAndOrder();
    testNoAllocation();
    return felitronics::test::report();
}
