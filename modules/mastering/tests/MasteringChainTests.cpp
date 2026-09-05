// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::mastering — the chain's own acceptance.
//
// The bar this suite is written against, and the reason each part of it is shaped the way it is:
//
//  * BLOCK INVARIANCE IS ASSERTED BIT-FOR-BIT, on an ADVERSARIAL partition. Not 64/512/4096 — those are
//    round, and an internal period can hide behind them. One sample at a time, primes, a random
//    partition with a fixed seed, blocks larger than anything prepared, and the whole programme in a
//    single call. Comparison is on the BIT PATTERN, so -0.0f against +0.0f is a failure rather than a
//    pass.
//  * THE FIXTURES ARE CHOSEN TO BREAK IT, from what was measured on the stages underneath (the full
//    write-up is .private/p6-findings.md): a burst into exact digital silence (the EQ / mono-bass
//    per-call flush), a band parked at 0 dB (which keeps filter state sitting on the flush threshold, so
//    the divergence lands mid-tone rather than in a tail), material at -160 dBFS (the compressor's RMS
//    flush floor), and 24-bit dither with auto-blank ON, which is what turned a 1e-15 difference into
//    three LSB of the export.
//  * EVERY CHECK IS MEANT TO BE VERIFIED BY REVERTING SOMETHING. Where a check exists only because a
//    specific mutation would otherwise survive, the mutation is named in a comment.
//
// The tests here use the chain's own arithmetic where that is what is under test. The checks that must
// NOT trust it — latency, alignment, the composition order — live in MasteringOracleTests.cpp, which
// re-derives them from outside.

#include <felitronics/mastering/OfflineRenderer.h>
#include <felitronics_test.h>

#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <vector>

static std::atomic<long long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::okNoAlloc;
using felitronics::test::group;

namespace
{
constexpr double kPi = 3.14159265358979323846;
using Buf = std::vector<std::vector<float>>;

std::uint32_t bits (float f) noexcept { return std::bit_cast<std::uint32_t> (f); }

// Bit-pattern equality, so a sign of zero is a difference. That is not pedantry: the clipper's own
// `mix = 0` bypass normalises -0.0f to +0.0f, which is exactly why the chain aligns it instead.
bool bitEqual (const Buf& a, const Buf& b) noexcept
{
    if (a.size() != b.size()) return false;
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        if (a[c].size() != b[c].size()) return false;
        for (std::size_t i = 0; i < a[c].size(); ++i)
            if (bits (a[c][i]) != bits (b[c][i])) return false;
    }
    return true;
}

int firstDiff (const Buf& a, const Buf& b) noexcept
{
    for (std::size_t c = 0; c < a.size(); ++c)
        for (std::size_t i = 0; i < a[c].size(); ++i)
            if (bits (a[c][i]) != bits (b[c][i])) return (int) i;
    return -1;
}

std::vector<float*> planes (Buf& b, int off = 0)
{
    std::vector<float*> p;
    p.reserve (b.size());
    for (auto& v : b) p.push_back (v.data() + off);
    return p;
}

//==============================================================================
// FIXTURES. Each one exists because something measured on the stages underneath makes it dangerous.
Buf silence (int nch, int n) { return Buf ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f)); }

// A loud burst, then EXACT digital silence. The shape that drives a filter tail through the 1e-15
// per-call flush threshold and, downstream, decides whether the dither's auto-blank fires.
Buf burstThenSilence (int nch, int n, int burst, unsigned seed = 12345u)
{
    Buf x = silence (nch, n);
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int i = 0; i < burst && i < n; ++i)
    {
        const float e = 0.7f * (float) std::sin (2.0 * kPi * 220.0 * i / 48000.0) + 0.2f * u (rng);
        for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) i] = (c == 0 ? e : -0.8f * e + 0.05f * u (rng));
    }
    return x;
}

// Programme that actually drives gain reduction in BOTH dynamic stages, so the check is not measuring
// a chain that happens to be sitting at unity.
Buf loudProgramme (int nch, int n, unsigned seed = 777u)
{
    Buf x = silence (nch, n);
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / 48000.0;
        const double env = 0.35 + 0.65 * std::fabs (std::sin (2.0 * kPi * 1.7 * t));
        const float e = (float) (env * (0.55 * std::sin (2.0 * kPi * 110.0 * t)
                                      + 0.30 * std::sin (2.0 * kPi * 2350.0 * t))) + 0.12f * u (rng);
        for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) i] = (c == 0 ? e : 0.9f * e + 0.05f * u (rng));
        if (i % 9600 == 0) for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) i] = (c == 0 ? 0.99f : -0.99f);
    }
    return x;
}

// -160 dBFS. The level at which the compressor's RMS flush floor lives; a partition-dependent flush
// there was measured at -7.5 dB of gain reduction in one call against 0.00 dB in per-sample calls.
Buf veryQuiet (int nch, int n)
{
    Buf x = silence (nch, n);
    for (int i = 0; i < n; ++i)
    {
        const float e = 1.0e-8f * (float) std::sin (2.0 * kPi * 997.0 * i / 48000.0);
        for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) i] = e;
    }
    return x;
}

//==============================================================================
// PARTITIONS. Deliberately not round.
std::vector<int> partOne (int n)   { return std::vector<int> ((std::size_t) n, 1); }
std::vector<int> partWhole (int n) { return { n }; }

std::vector<int> partFixed (int n, int b)
{
    std::vector<int> p;
    for (int off = 0; off < n; off += b) p.push_back (std::min (b, n - off));
    return p;
}

// Cycles through primes, so no internal period can stay in step with the caller.
std::vector<int> partPrimes (int n)
{
    static const int primes[] = { 1, 2, 3, 5, 7, 11, 13, 17, 23, 31, 61, 127, 251, 509, 1021, 2039 };
    std::vector<int> p;
    std::size_t k = 0;
    for (int off = 0; off < n; )
    {
        const int b = std::min (primes[k % (sizeof primes / sizeof primes[0])], n - off);
        p.push_back (b);
        off += b;
        ++k;
    }
    return p;
}

// Random sizes from a fixed seed, including sizes far larger than anything the chain was prepared for.
std::vector<int> partRandom (int n, unsigned seed, int hi)
{
    std::vector<int> p;
    std::mt19937 rng (seed);
    std::uniform_int_distribution<int> d (1, hi);
    for (int off = 0; off < n; )
    {
        const int b = std::min (d (rng), n - off);
        p.push_back (b);
        off += b;
    }
    return p;
}

//==============================================================================
mastering::MasteringChainParams voicing()
{
    mastering::MasteringChainParams p;
    p.inputGainDb      = 1.5;
    p.preLimiterGainDb = 3.0;

    // Band 0 does real work; band 1 is parked at 0 dB ON PURPOSE — that is the configuration whose
    // filter state sits on the flush threshold, so the block-dependence lands mid-programme instead of
    // in a decaying tail. Without it the invariance checks below are much weaker than they look.
    p.eqBands[0].on = true;
    p.eqBands[0].type = eq::FilterType::Bell;
    p.eqBands[0].lane (eq::Lane::Stereo).on = true;
    p.eqBands[0].lane (eq::Lane::Stereo).freq = 3200.0;
    p.eqBands[0].lane (eq::Lane::Stereo).Q = 0.9;
    p.eqBands[0].lane (eq::Lane::Stereo).gainDb = 2.5;

    p.eqBands[1].on = true;
    p.eqBands[1].type = eq::FilterType::Bell;
    p.eqBands[1].lane (eq::Lane::Stereo).on = true;
    p.eqBands[1].lane (eq::Lane::Stereo).freq = 400.0;
    p.eqBands[1].lane (eq::Lane::Stereo).Q = 1.0;
    p.eqBands[1].lane (eq::Lane::Stereo).gainDb = 0.0;      // parked — see above

    p.monoBass.enabled = true;
    p.monoBass.frequencyHz = 110.0f;
    p.monoBass.lowWidth = 0.0f;

    p.compressor.thresholdDb = -18.0;
    p.compressor.ratio       = 2.5;
    p.compressor.attackMs    = 12.0;
    p.compressor.releaseMs   = 180.0;
    p.compressor.makeupDb    = 1.0;

    p.clipper.driveDb = 4.0f;
    p.clipper.mix     = 1.0f;

    p.limiter.ceilingDbTp = -1.0;
    p.limiter.releaseMs   = 60.0;

    p.dither.bits      = 24;
    p.dither.shaping   = dither::NoiseShaping::Weighted;
    p.dither.autoBlank = true;                              // the amplifier — see the file header
    return p;
}

mastering::MasteringChainConfig fullConfig (int K = 256)
{
    mastering::MasteringChainConfig c;
    c.internalBlock = K;
    c.eq = c.compressor = c.limiter = c.dither = true;
    c.monoBass = true;
    c.clipper  = true;
    c.sidechainHpfHz = 120.0;
    return c;
}

// Stream `x` through a fresh chain using the given partition; return the streaming output.
Buf runPartition (const mastering::MasteringChainConfig& cfg, const mastering::MasteringChainParams& prm,
                  const Buf& x, const std::vector<int>& part, bool& prepOk)
{
    const int nch = (int) x.size();
    mastering::MasteringChain chain;
    prepOk = chain.prepare (48000.0, nch, cfg);
    Buf y = x;
    if (! prepOk) return y;
    chain.setParams (prm);
    int off = 0;
    for (int b : part)
    {
        auto p = planes (y, off);
        if (! chain.process (p.data(), nch, b)) { prepOk = false; return y; }
        off += b;
    }
    return y;
}
} // namespace

//==============================================================================
static void testBlockInvariance()
{
    group ("block invariance — bit-for-bit on an adversarial partition");

    const int nch = 2, n = 60000;
    const auto cfg = fullConfig();
    const auto prm = voicing();

    struct Fixture { const char* name; Buf x; };
    std::vector<Fixture> fixtures;
    fixtures.push_back ({ "burst -> digital silence", burstThenSilence (nch, n, 4000) });
    fixtures.push_back ({ "loud programme (real GR)", loudProgramme (nch, n) });
    fixtures.push_back ({ "-160 dBFS tone",           veryQuiet (nch, n) });
    fixtures.push_back ({ "digital silence",          silence (nch, n) });
    {
        Buf dc = silence (nch, n);
        for (int c = 0; c < nch; ++c) std::fill (dc[(std::size_t) c].begin(), dc[(std::size_t) c].end(), 0.35f);
        fixtures.push_back ({ "DC", std::move (dc) });
    }
    {
        Buf imp = silence (nch, n);
        imp[0][100] = 0.95f;
        imp[1][n - 1] = -0.95f;                            // one at the very LAST sample: the tail witness
        fixtures.push_back ({ "impulses incl. the last sample", std::move (imp) });
    }

    struct Part { const char* name; std::vector<int> p; };
    for (auto& f : fixtures)
    {
        bool okPrep = false;
        const Buf ref = runPartition (cfg, prm, f.x, partWhole (n), okPrep);
        ok (okPrep, std::string ("prepare + whole-programme call: ") + f.name);

        std::vector<Part> parts;
        parts.push_back ({ "1 sample",        partOne (n) });
        parts.push_back ({ "primes",          partPrimes (n) });
        parts.push_back ({ "random(seed 1)",  partRandom (n, 1u, 3000) });
        parts.push_back ({ "random(seed 99)", partRandom (n, 99u, 700) });
        parts.push_back ({ "1021",            partFixed (n, 1021) });
        parts.push_back ({ "K-1 = 255",       partFixed (n, 255) });
        parts.push_back ({ "K = 256",         partFixed (n, 256) });
        parts.push_back ({ "K+1 = 257",       partFixed (n, 257) });
        parts.push_back ({ "20011 (> K)",     partFixed (n, 20011) });

        for (auto& part : parts)
        {
            bool p2 = false;
            const Buf got = runPartition (cfg, prm, f.x, part.p, p2);
            const bool same = p2 && bitEqual (ref, got);
            ok (same, std::string ("bit-identical: ") + f.name + " @ " + part.name
                      + (same ? "" : "  (first differing sample " + std::to_string (firstDiff (ref, got)) + ")"));
        }
    }

    // MUTATION WITNESS. Removing the internal quantum — handing the caller's block straight to the
    // stages — must fail the checks above. Measured before this module existed: 4349 differing samples
    // between block 4096 and block 1 over this same chain, worst 3.576e-07 (three LSB of 24-bit),
    // because the EQ's per-call flush decides whether the dither's auto-blank sees exact zero.
    // Shrinking the quantum must NOT change that the property holds, only the latency.
    for (int K : { 8, 64, 512, 4096 })
    {
        const auto c2 = fullConfig (K);
        const Buf x = burstThenSilence (nch, n, 4000);
        bool a = false, b = false;
        const Buf refK = runPartition (c2, prm, x, partWhole (n), a);
        const Buf gotK = runPartition (c2, prm, x, partPrimes (n), b);
        ok (a && b && bitEqual (refK, gotK), "bit-identical at internalBlock = " + std::to_string (K));
    }
}

//==============================================================================
static void testBypassNull()
{
    group ("bypass null — every stage present, every stage bypassed");

    const int nch = 2, n = 40000;
    const auto cfg = fullConfig();
    auto prm = voicing();
    prm.bypassEq = prm.bypassMonoBass = prm.bypassCompressor = true;
    prm.bypassClipper = prm.bypassLimiter = prm.bypassDither = true;
    prm.inputGainDb = prm.preLimiterGainDb = 0.0;          // a gain is not a bypassable stage; 0 dB is its identity

    Buf x = loudProgramme (nch, n);
    x[0][0]   = -0.0f;                                     // the sign of zero must survive the whole path
    x[1][123] = -0.0f;

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare with every stage present");
    chain.setParams (prm);

    const int D = chain.latencySamples();
    ok (D > 0, "a fully bypassed chain still declares its topology's latency (" + std::to_string (D) + ")");

    // The STREAMING form: the output must be the input delayed by exactly D, bit for bit.
    Buf y = x;
    {
        int off = 0;
        for (int b : partPrimes (n)) { auto p = planes (y, off); chain.process (p.data(), nch, b); off += b; }
    }
    long long bad = 0;
    for (int c = 0; c < nch; ++c)
        for (int i = D; i < n; ++i)
            if (bits (y[(std::size_t) c][(std::size_t) i]) != bits (x[(std::size_t) c][(std::size_t) (i - D)])) ++bad;
    ok (bad == 0, "streaming: out[i] is bit-identical to in[i-D] (" + std::to_string (bad) + " differ)");

    long long primed = 0;
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < D; ++i) if (! core::exactlyEqual (y[(std::size_t) c][(std::size_t) i], 0.0f)) ++primed;
    ok (primed == 0, "streaming: the first D samples are the chain's own priming, exactly zero");

    // The OFFLINE form: an exact identity, same length, no shift, tail included.
    mastering::OfflineRenderer r;
    r.prepare (nch, 997);
    Buf out = silence (nch, n);
    auto ip = planes (const_cast<Buf&> (x));
    auto op = planes (out);
    ok (r.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render() accepts the fully bypassed chain");
    ok (bitEqual (x, out), "offline: a fully bypassed render is a bit-exact IDENTITY — length and alignment");

    // MUTATION WITNESS: dropping the DryAligner on the limiter, or feeding it only while bypassed,
    // shifts this by its 111 samples and the identity fails. Holding the aligner at the WRONG delay
    // fails it too. That is the pair of mutations this single check is aimed at.
}

//==============================================================================
static void testFlushIsProcessOfZeros()
{
    group ("flush is process() over latencySamples() zeros — by construction, and checked");

    const int nch = 2, n = 20000;
    const auto cfg = fullConfig();
    const auto prm = voicing();
    const Buf x = loudProgramme (nch, n);

    mastering::MasteringChain a, b;
    ok (a.prepare (48000.0, nch, cfg) && b.prepare (48000.0, nch, cfg), "prepare two identical chains");
    a.setParams (prm);
    b.setParams (prm);
    const int D = a.latencySamples();

    // (a) stream x, then flush().
    Buf ya = x;
    { auto p = planes (ya); a.process (p.data(), nch, n); }
    Buf ta = silence (nch, D);
    { auto p = planes (ta); ok (a.flush (p.data(), nch, D) == D, "flush() writes exactly latencySamples() frames"); }

    // (b) stream x followed by D zeros, no flush() at all.
    Buf xb = x;
    for (int c = 0; c < nch; ++c) xb[(std::size_t) c].resize ((std::size_t) (n + D), 0.0f);
    { auto p = planes (xb); b.process (p.data(), nch, n + D); }

    long long bad = 0;
    for (int c = 0; c < nch; ++c)
    {
        for (int i = 0; i < n; ++i) if (bits (ya[(std::size_t) c][(std::size_t) i]) != bits (xb[(std::size_t) c][(std::size_t) i])) ++bad;
        for (int i = 0; i < D; ++i) if (bits (ta[(std::size_t) c][(std::size_t) i]) != bits (xb[(std::size_t) c][(std::size_t) (n + i)])) ++bad;
    }
    ok (bad == 0, "flush() output equals processing D zeros, bit for bit (" + std::to_string (bad) + " differ)");

    // Capacity is honoured rather than overrun, and a chain that is not prepared writes nothing.
    Buf small = silence (nch, 4);
    { auto p = planes (small); ok (a.flush (p.data(), nch, 4) == 4, "flush() honours a capacity below the latency"); }
    mastering::MasteringChain unprepared;
    { Buf s = silence (nch, 8); auto p = planes (s); ok (unprepared.flush (p.data(), nch, 8) == 0, "flush() on an unprepared chain writes nothing"); }
}

//==============================================================================
static void testRendererContract()
{
    group ("OfflineRenderer — length, alignment, the tail, in == out, re-runs");

    const int nch = 2;
    const auto cfg = fullConfig();
    const auto prm = voicing();

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare");
    chain.setParams (prm);
    const int D = chain.latencySamples();

    // The renderer's own block size must not change a single bit — the chain re-blocks anyway.
    {
        const int n = 30000;
        const Buf x = loudProgramme (nch, n);
        Buf ref = silence (nch, n), got = silence (nch, n);
        mastering::OfflineRenderer r1, r2;
        r1.prepare (nch, 4096);
        r2.prepare (nch, 1);
        auto ip = planes (const_cast<Buf&> (x));
        { auto op = planes (ref); ok (r1.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render at blockSize 4096"); }
        { auto op = planes (got); ok (r2.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render at blockSize 1"); }
        ok (bitEqual (ref, got), "renderer blockSize does not change the result");

        mastering::OfflineRenderer r3;
        r3.prepare (nch, 65536);
        Buf big = silence (nch, n);
        { auto op = planes (big); r3.render (chain, (const float* const*) ip.data(), op.data(), nch, n); }
        ok (bitEqual (ref, big), "renderer blockSize 65536 does not change the result either");

        // Re-running must reproduce exactly — P7 renders, measures, adjusts one gain, renders again.
        Buf again = silence (nch, n);
        { auto op = planes (again); r1.render (chain, (const float* const*) ip.data(), op.data(), nch, n); }
        ok (bitEqual (ref, again), "a second render of the same input is bit-identical (dither reseeded)");
    }

    // THE TAIL. A click near the end of the programme must come out — the measured defect of the ffmpeg
    // chain this replaces was that the last D samples never left it at all: a click 4 ms before the end
    // vanished, one 5 ms before survived. Here the click is on the LAST sample.
    {
        const int n = 8000;
        Buf x = silence (nch, n);
        for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) (n - 1)] = 0.8f;
        Buf out = silence (nch, n);
        mastering::OfflineRenderer r;
        r.prepare (nch, 512);
        auto ip = planes (x);
        auto op = planes (out);
        ok (r.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render the last-sample click");
        double tail = 0.0;
        for (int c = 0; c < nch; ++c) tail = std::max (tail, (double) std::fabs (out[(std::size_t) c][(std::size_t) (n - 1)]));
        ok (tail > 0.05, "a click on the very LAST sample reaches the output (peak " + std::to_string (tail) + ")");
        // MUTATION WITNESS: drop the drain from render() and this is exactly 0.
    }

    // Short and degenerate lengths, including frames < latency.
    for (int n : { 0, 1, 7, D - 1, D, D + 1 })
    {
        if (n < 0) continue;
        Buf x = silence (nch, std::max (n, 1));
        for (int i = 0; i < n; ++i) for (int c = 0; c < nch; ++c) x[(std::size_t) c][(std::size_t) i] = 0.4f;
        Buf out = silence (nch, std::max (n, 1));
        mastering::OfflineRenderer r;
        r.prepare (nch, 333);
        auto ip = planes (x);
        auto op = planes (out);
        ok (r.render (chain, (const float* const*) ip.data(), op.data(), nch, n),
            "render frames = " + std::to_string (n));
    }

    // in == out.
    {
        const int n = 12000;
        const Buf src = loudProgramme (nch, n);
        Buf ref = silence (nch, n), inplace = src;
        mastering::OfflineRenderer r;
        r.prepare (nch, 640);
        auto ip = planes (const_cast<Buf&> (src));
        { auto op = planes (ref); r.render (chain, (const float* const*) ip.data(), op.data(), nch, n); }
        { auto pp = planes (inplace); r.render (chain, (const float* const*) pp.data(), pp.data(), nch, n); }
        ok (bitEqual (ref, inplace), "render() with in == out gives the same result as out-of-place");
    }
}

//==============================================================================
static void testSequences()
{
    group ("sequences — channel count, reset, failed prepare, oversized block");

    const int nch = 2, n = 6000;
    const auto cfg = fullConfig();
    const auto prm = voicing();

    // A CHANNEL COUNT THAT IS NOT THE PREPARED ONE IS REFUSED, and a refused call is indistinguishable
    // from one never made. That is what makes "stereo -> mono -> stereo" provable rather than a
    // best-effort reset: [A, refused, B] must equal [A, B], bit for bit.
    {
        const Buf x = loudProgramme (nch, n);
        mastering::MasteringChain a, b;
        ok (a.prepare (48000.0, nch, cfg) && b.prepare (48000.0, nch, cfg), "prepare two stereo chains");
        a.setParams (prm);
        b.setParams (prm);

        Buf ya = x, yb = x;
        const int half = n / 2;
        {
            auto p = planes (ya);
            ok (a.process (p.data(), nch, half), "first half accepted");
            // ...a mono call in the middle, which must change nothing at all
            std::vector<float> mono ((std::size_t) 512, 0.5f);
            float* m[1] { mono.data() };
            ok (! a.process (m, 1, 512), "a mono call on a stereo chain is REFUSED");
            ok (core::exactlyEqual (mono[0], 0.5f) && core::exactlyEqual (mono[511], 0.5f), "the refused call left its buffer untouched");
            std::vector<float> three ((std::size_t) 256, 0.25f);
            float* t3[3] { three.data(), three.data(), three.data() };
            ok (! a.process (t3, 3, 256), "a 3-channel call on a stereo chain is REFUSED");
            auto p2 = planes (ya, half);
            ok (a.process (p2.data(), nch, n - half), "second half accepted");
        }
        { auto p = planes (yb); b.process (p.data(), nch, n); }
        ok (bitEqual (ya, yb), "a refused call is indistinguishable from one never made");
    }

    // reset() puts the chain back to a fresh state — including the EQ's parameter smoothers, which
    // EqBand::reset() does NOT re-snap on its own.
    {
        const Buf x = loudProgramme (nch, n);
        mastering::MasteringChain a, fresh;
        ok (a.prepare (48000.0, nch, cfg) && fresh.prepare (48000.0, nch, cfg), "prepare for the reset check");
        a.setParams (prm);
        fresh.setParams (prm);
        Buf warm = x;
        { auto p = planes (warm); a.process (p.data(), nch, n); }   // dirty the state
        a.reset();
        Buf ya = x, yf = x;
        { auto p = planes (ya); a.process (p.data(), nch, n); }
        { auto p = planes (yf); fresh.process (p.data(), nch, n); }
        ok (bitEqual (ya, yf), "after reset() the chain is bit-identical to a freshly prepared one");
        // MUTATION WITNESS: drop `paramsDirty_ = true` from reset() and this fails, because the EQ's
        // smoothers stay wherever the first pass left them.
    }

    // A FAILED prepare leaves the chain unusable and reporting nothing, and does not resurrect a
    // previously good one.
    {
        mastering::MasteringChain c;
        ok (c.prepare (48000.0, nch, cfg), "a good prepare");
        ok (c.latencySamples() > 0, "...reports a latency");
        ok (! c.prepare (std::nan (""), nch, cfg), "prepare(NaN sample rate) is REFUSED");
        ok (! c.isPrepared() && c.latencySamples() == 0, "a failed prepare leaves the chain unprepared, latency 0");
        Buf x = silence (nch, 64);
        for (int ch = 0; ch < nch; ++ch) std::fill (x[(std::size_t) ch].begin(), x[(std::size_t) ch].end(), 0.5f);
        { auto p = planes (x); ok (! c.process (p.data(), nch, 64), "process() on the failed chain is refused"); }
        ok (core::exactlyEqual (x[0][0], 0.5f), "...and left the buffer untouched");

        // The other refusals, each spelled positively so a NaN cannot slip through.
        mastering::MasteringChain d;
        ok (! d.prepare (0.0, nch, cfg), "prepare(sampleRate 0) refused");
        ok (! d.prepare (-48000.0, nch, cfg), "prepare(negative sample rate) refused");
        ok (! d.prepare (48000.0, 0, cfg), "prepare(0 channels) refused");
        ok (! d.prepare (48000.0, core::kMaxChannels + 1, cfg), "prepare(too many channels) refused");
        auto bad = cfg; bad.internalBlock = 1;
        ok (! d.prepare (48000.0, nch, bad), "prepare(internalBlock below the floor) refused");
        bad = cfg; bad.internalBlock = 1 << 20;
        ok (! d.prepare (48000.0, nch, bad), "prepare(internalBlock above the cap) refused");
        bad = cfg; bad.tapsPerPhase = 2;
        ok (! d.prepare (48000.0, nch, bad), "prepare(tapsPerPhase below the stages' floor) refused");
        bad = cfg; bad.compressorLookaheadMs = std::nan ("");
        ok (! d.prepare (48000.0, nch, bad), "prepare(NaN lookahead) refused");
        bad = cfg; bad.sidechainHpfHz = 30000.0;
        ok (! d.prepare (48000.0, nch, bad), "prepare(sidechain corner above Nyquist) refused");

        // MONO-BASS ON A NON-STEREO CHAIN IS REFUSED, not silently ignored: stereo::MonoBass leaves any
        // buffer that is not exactly two channels untouched, so accepting this would report an enabled
        // stage that does nothing.
        auto mb = cfg; mb.monoBass = true;
        ok (! d.prepare (48000.0, 1, mb), "prepare(mono-bass on a mono chain) refused");
        auto nomb = cfg; nomb.monoBass = false;
        ok (d.prepare (48000.0, 1, nomb), "...and the same mono chain prepares once mono-bass is off");
    }
}

//==============================================================================
static void testPoisonGate()
{
    group ("the gate — a bad sample is bit-identical to the sanitised one it stands for");

    const int nch = 2, n = 12000;
    const auto prm = voicing();

    // The strong form, and the one the plan asks for: not "we survived it" (which `isfinite` on the
    // output would show) but "there was no damage at all" — the poisoned stream and the explicitly
    // sanitised stream must be the SAME BITS. Checked with every stage on, and again with every stage
    // bypassed, because with everything bypassed an ungated Inf would otherwise reach the output.
    for (int variant = 0; variant < 2; ++variant)
    {
        auto p = prm;
        if (variant == 1)
        {
            p.bypassEq = p.bypassMonoBass = p.bypassCompressor = true;
            p.bypassClipper = p.bypassLimiter = p.bypassDither = true;
        }
        const char* what = variant == 0 ? "all stages active" : "all stages bypassed";

        Buf clean = loudProgramme (nch, n);
        Buf dirty = clean;
        const float poison[] = { std::numeric_limits<float>::quiet_NaN(),
                                 std::numeric_limits<float>::infinity(),
                                 -std::numeric_limits<float>::infinity(),
                                 3.0e38f, -3.0e38f };
        // Placed at a quantum boundary, one before it, one after it, and in the middle of one — the
        // positions where a guard that ran per call rather than per sample would behave differently.
        const int at[] = { 0, 255, 256, 257, 700, 5000, n - 1 };
        int k = 0;
        for (int idx : at)
        {
            const float bad = poison[(std::size_t) (k % 5)];
            dirty[0][(std::size_t) idx] = bad;
            dirty[1][(std::size_t) idx] = poison[(std::size_t) ((k + 2) % 5)];
            clean[0][(std::size_t) idx] = std::clamp (std::isfinite (bad) ? bad : 0.0f, -1.0e6f, 1.0e6f);
            const float bad2 = poison[(std::size_t) ((k + 2) % 5)];
            clean[1][(std::size_t) idx] = std::clamp (std::isfinite (bad2) ? bad2 : 0.0f, -1.0e6f, 1.0e6f);
            ++k;
        }

        const auto cfg = fullConfig();
        bool a = false, b = false;
        const Buf yClean = runPartition (cfg, p, clean, partPrimes (n), a);
        const Buf yDirty = runPartition (cfg, p, dirty, partPrimes (n), b);
        ok (a && b && bitEqual (yClean, yDirty),
            std::string ("poisoned input is bit-identical to the sanitised input it stands for — ") + what);

        long long nonFinite = 0;
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < n; ++i) if (! std::isfinite (yDirty[(std::size_t) c][(std::size_t) i])) ++nonFinite;
        ok (nonFinite == 0, std::string ("...and nothing non-finite leaves the chain — ") + what);
    }

    // Bit-transparency of the gate itself: ordinary audio takes neither branch.
    {
        const auto cfg = fullConfig();
        auto p = prm;
        Buf x = loudProgramme (nch, 4000);
        x[0][10] = 1.0e6f; x[0][11] = -1.0e6f;             // exactly on the clamp bounds
        bool a = false, b = false;
        const Buf y1 = runPartition (cfg, p, x, partWhole (4000), a);
        const Buf y2 = runPartition (cfg, p, x, partOne (4000), b);
        ok (a && b && bitEqual (y1, y2), "samples exactly on the gate's bounds stay block-invariant");
    }
}

//==============================================================================
// The four checks below exist because a mutation of the module survived the suite without them. Each
// names the mutation it was written to kill.
static void testMutationGaps()
{
    group ("gaps found by mutating the module");

    // ---- (1) A STAGE COMING BACK OUT OF BYPASS MUST NOT REPLAY WHAT IT HELD -------------------
    // Skipping a stage FREEZES its state; the limiter's oversampler history would then be convolved
    // with the fresh signal on re-entry and emit audio from before the gap. That is the defect measured
    // at +19.76 dB over the ceiling in the limiter itself (#119) and at 0.93 in the saturator.
    // MUTATION KILLED: dropping `if (bypassChanged_.limiter) lim_.reset();`.
    {
        const int nch = 2;
        mastering::MasteringChainConfig cfg;
        cfg.internalBlock = 256;
        cfg.eq = cfg.monoBass = cfg.compressor = cfg.clipper = cfg.dither = false;
        cfg.limiter = true;                                // the only stage in the path

        mastering::MasteringChain chain;
        ok (chain.prepare (48000.0, nch, cfg), "prepare a limiter-only chain");
        mastering::MasteringChainParams p;
        p.limiter.ceilingDbTp = -1.0;
        chain.setParams (p);

        Buf loud = silence (nch, 4096);
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < 4096; ++i)
                loud[(std::size_t) c][(std::size_t) i] = 0.95f * (float) std::sin (2.0 * kPi * 300.0 * i / 48000.0);

        Buf a = loud; { auto pl = planes (a); chain.process (pl.data(), nch, 4096); }          // active, loud
        p.bypassLimiter = true; chain.setParams (p);
        Buf b = loud; { auto pl = planes (b); chain.process (pl.data(), nch, 4096); }          // bypassed, loud
        for (int k = 0; k < 4; ++k)                                                             // drain everything
        { Buf z = silence (nch, 4096); auto pl = planes (z); chain.process (pl.data(), nch, 4096); }
        p.bypassLimiter = false; chain.setParams (p);
        Buf out = silence (nch, 4096);
        { auto pl = planes (out); chain.process (pl.data(), nch, 4096); }                       // active again, SILENCE in

        double peak = 0.0;
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < 4096; ++i) peak = std::max (peak, (double) std::fabs (out[(std::size_t) c][(std::size_t) i]));
        ok (peak < 1.0e-9, "a stage returning from bypass emits nothing from before the gap (peak "
                           + std::to_string (peak) + ")");
    }

    // ---- (1b) THE ALIGNER MUST BE WARM THE INSTANT A STAGE IS BYPASSED -------------------------
    // `DryAligner` is documented as needing to be fed EVERY block, active or not: a ring first written
    // at the moment it is first read holds zeros for its whole delay, so a toggle into bypass would
    // punch a hole of exactly the stage's latency into the programme. Checking the null only well AFTER
    // the toggle misses that entirely — this checks it from the FIRST bypassed sample.
    // MUTATIONS KILLED: moving either aligner's advance() inside the bypass branch.
    {
        const int nch = 2, K = 256, T = 8 * K;             // toggle on a quantum boundary
        mastering::MasteringChainConfig cfg;
        cfg.internalBlock = K;
        cfg.eq = cfg.monoBass = cfg.compressor = cfg.dither = false;
        cfg.clipper = cfg.limiter = true;

        mastering::MasteringChain chain;
        ok (chain.prepare (48000.0, nch, cfg), "prepare a clipper+limiter chain");
        mastering::MasteringChainParams p;
        p.clipper.driveDb     = 6.0f;                      // both stages doing real work before the toggle
        p.limiter.ceilingDbTp = -6.0;
        chain.setParams (p);

        const int n = T + 6 * K;
        Buf x = loudProgramme (nch, n, 4242u);
        const Buf src = x;
        { auto pl = planes (x); chain.process (pl.data(), nch, T); }          // active
        p.bypassClipper = p.bypassLimiter = true;
        chain.setParams (p);                                                  // takes effect at sample T
        { auto pl = planes (x, T); chain.process (pl.data(), nch, n - T); }   // bypassed

        const int D = chain.latencySamples();

        // TWO different things have to hold, and only asserting the second is what let the cold-ring
        // mutations through the first time this was written.
        //
        // (a) NO HOLE. The toggle lands on the quantum holding input [T, T+K), whose result leaves over
        //     output [T+K, T+2K). Across [T+K, T+D) the aligners are still emitting what actually flowed
        //     through them — the CLIPPER-PROCESSED audio, not the raw input — so this window cannot be
        //     nulled against the input. What it can be is CONTINUOUS: a ring first written at the moment
        //     it is first read holds zeros for its whole delay, so the failure has a shape, a run of
        //     exact zeros as long as the stage's latency, and that is what is checked.
        int longestZeroRun = 0;
        for (int c = 0; c < nch; ++c)
        {
            int run = 0;
            for (int i = T + K; i < T + D; ++i)
            {
                run = core::exactlyEqual (x[(std::size_t) c][(std::size_t) i], 0.0f) ? run + 1 : 0;
                longestZeroRun = std::max (longestZeroRun, run);
            }
        }
        ok (longestZeroRun < 16, "a bypass toggle punches no hole: longest run of exact zeros across the "
                                 "transition is " + std::to_string (longestZeroRun) + " samples");

        // (b) AND ONCE THE PRE-TOGGLE CONTENT HAS DRAINED, the bypassed chain is a pure delay again.
        long long bad = 0;
        int firstBad = -1;
        for (int c = 0; c < nch; ++c)
            for (int i = T + D; i < n; ++i)
                if (bits (x[(std::size_t) c][(std::size_t) i]) != bits (src[(std::size_t) c][(std::size_t) (i - D)]))
                { ++bad; if (firstBad < 0) firstBad = i - (T + D); }
        ok (bad == 0, "after the transition the bypassed chain is exactly the input delayed by D ("
                      + std::to_string (bad) + " differ, first at +" + std::to_string (firstBad) + ")");
    }

    // ---- (2) reset() MUST RESTORE THE PARAMETER SMOOTHERS, NOT ONLY THE FILTER STATE -----------
    // `EqBand::reset()` clears filter state and deliberately leaves freq/Q/gain where they are, so a
    // reset in the middle of a 30 ms ramp used to resume it. Measured before the fix: 0.51 of
    // difference, full scale, against a freshly prepared chain.
    // MUTATION KILLED: dropping `forceSnap_ = true;` from reset().
    {
        const int nch = 1, n = 8000;
        mastering::MasteringChainConfig cfg;
        cfg.internalBlock = 256;
        cfg.eq = true;
        cfg.monoBass = cfg.compressor = cfg.clipper = cfg.limiter = cfg.dither = false;

        auto band = [] (double g) {
            eq::BandParams b;
            b.on = true; b.type = eq::FilterType::Bell;
            b.lane (eq::Lane::Stereo).on = true;
            b.lane (eq::Lane::Stereo).freq = 1000.0;
            b.lane (eq::Lane::Stereo).Q = 1.0;
            b.lane (eq::Lane::Stereo).gainDb = g;
            return b;
        };
        Buf x = silence (nch, n);
        for (int i = 0; i < n; ++i) x[0][(std::size_t) i] = 0.4f * (float) std::sin (2.0 * kPi * 1000.0 * i / 48000.0);

        mastering::MasteringChainParams A, B;
        A.eqBands[0] = band (0.0);
        B.eqBands[0] = band (9.0);

        mastering::MasteringChain c1, c2;
        ok (c1.prepare (48000.0, nch, cfg) && c2.prepare (48000.0, nch, cfg), "prepare for the mid-ramp reset check");
        c1.setParams (A);
        { Buf y = x; auto pl = planes (y); c1.process (pl.data(), nch, 1024); }
        c1.setParams (B);                                  // starts the ramp
        { Buf y = x; auto pl = planes (y); c1.process (pl.data(), nch, 512); }   // stop ~10 ms into it
        c1.reset();
        Buf y1 = x; { auto pl = planes (y1); c1.process (pl.data(), nch, n); }

        c2.setParams (B);
        Buf y2 = x; { auto pl = planes (y2); c2.process (pl.data(), nch, n); }
        ok (bitEqual (y1, y2), "reset() in the middle of a parameter ramp restores a FRESH chain");
    }

    // ---- (3) THE CHAIN VALIDATES THE SAMPLE RATE ITSELF ----------------------------------------
    // Two of the stages accept a NaN rate (their guards read `sampleRate <= 0.0`, which is false for a
    // NaN) and go on to emit NaN. With a compressor or a limiter in the topology their positively
    // spelled guards mask that, which is why the check has to run on a topology that has NEITHER.
    // MUTATION KILLED: rewriting the chain's own guard as `sampleRate <= 0.0`.
    {
        mastering::MasteringChainConfig cfg;
        cfg.eq = true;
        cfg.monoBass = cfg.compressor = cfg.clipper = cfg.limiter = cfg.dither = false;
        mastering::MasteringChain c;
        ok (! c.prepare (std::nan (""), 2, cfg), "prepare(NaN) is refused with no stage left to catch it");
        ok (! c.prepare (std::numeric_limits<double>::infinity(), 2, cfg), "prepare(inf) refused likewise");
        ok (c.prepare (48000.0, 2, cfg), "...and the same EQ-only topology prepares at a real rate");
    }

    // ---- (4) THE RENDERER'S FORMULA, NULLED AGAINST ITS OWN DEFINITION -------------------------
    // out[n] == y[n + D], where y is the chain's output for the input followed by D ZEROS. Computed
    // here independently of render(), so a drain that feeds the wrong thing fails even though every
    // render()-against-render() comparison still agrees.
    // MUTATION KILLED: draining with the input's last sample instead of zeros.
    {
        const int nch = 2, n = 20000;
        const auto cfg = fullConfig();
        const auto prm = voicing();
        const Buf x = loudProgramme (nch, n);

        mastering::MasteringChain a, b;
        ok (a.prepare (48000.0, nch, cfg) && b.prepare (48000.0, nch, cfg), "prepare two chains for the formula null");
        a.setParams (prm);
        b.setParams (prm);
        const int D = a.latencySamples();

        Buf padded = x;                                    // x followed by D zeros — the definition
        for (int c = 0; c < nch; ++c) padded[(std::size_t) c].resize ((std::size_t) (n + D), 0.0f);
        { auto pl = planes (padded); a.process (pl.data(), nch, n + D); }

        Buf out = silence (nch, n);
        mastering::OfflineRenderer r;
        r.prepare (nch, 512);
        auto ip = planes (const_cast<Buf&> (x));
        auto op = planes (out);
        ok (r.render (b, (const float* const*) ip.data(), op.data(), nch, n), "render for the formula null");

        long long bad = 0;
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < n; ++i)
                if (bits (out[(std::size_t) c][(std::size_t) i]) != bits (padded[(std::size_t) c][(std::size_t) (i + D)])) ++bad;
        ok (bad == 0, "render(x)[n] == process(x followed by D zeros)[n + D], bit for bit ("
                      + std::to_string (bad) + " differ)");
    }
}

//==============================================================================
static void testRtSafety()
{
    group ("RT safety — no allocation anywhere on the streaming path");

    const int nch = 2, n = 30000;
    const auto cfg = fullConfig();
    auto prm = voicing();

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare (this one DOES allocate, by design)");
    chain.setParams (prm);
    Buf x = loudProgramme (nch, n);
    Buf tail = silence (nch, chain.latencySamples());
    const std::vector<int> part = partPrimes (n);

    // EVERYTHING the test itself needs is built here, outside the counted region — the plane arrays are
    // raw stack pointers rather than the `planes()` helper, which allocates a vector per call and would
    // otherwise be measured as the chain's doing. A counted region that includes the harness's own
    // allocations cannot fail for the right reason.
    float* px[core::kMaxChannels] {};
    float* pt[core::kMaxChannels] {};
    std::vector<float> mono ((std::size_t) 64, 0.0f);
    float* pm[1] { mono.data() };
    for (int c = 0; c < nch; ++c) pt[c] = tail[(std::size_t) c].data();
    auto setPlanes = [&] (int off) { for (int c = 0; c < nch; ++c) px[c] = x[(std::size_t) c].data() + off; };

    // Warm up outside the counted region so a first-touch page fault or a lazily-built static cannot be
    // mistaken for an allocation.
    setPlanes (0);
    chain.process (px, nch, 1024);
    chain.reset();

    const long long before = g_allocs.load();
    setPlanes (0);
    chain.process (px, nch, n);                            // one whole-programme call
    {
        int off = 0;
        for (int b : part) { setPlanes (off); chain.process (px, nch, b); off += b; }
    }
    prm.bypassLimiter = true;  chain.setParams (prm);      // a bypass toggle, mid-stream
    setPlanes (0); chain.process (px, nch, 4096);
    prm.bypassLimiter = false; chain.setParams (prm);
    setPlanes (0); chain.process (px, nch, 4096);
    chain.process (pm, 1, 64);                             // a refused call
    chain.flush (pt, nch, (int) tail[0].size());
    chain.reset();
    (void) chain.resolved();
    const long long after = g_allocs.load();

    okNoAlloc (after == before, "process() / flush() / setParams() / reset() / resolved() allocate nothing ("
                                + std::to_string (after - before) + ")");
}

//==============================================================================
static void testResolvedReadback()
{
    group ("resolved() — what the chain actually applied");

    const int nch = 2;
    auto cfg = fullConfig();
    cfg.limiterLookaheadMs = 1.0;
    cfg.oversampleFactor   = 4;
    cfg.tapsPerPhase       = 64;

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepare");
    auto prm = voicing();
    prm.limiter.ceilingDbTp = 1.0e308;                     // absurd on purpose — the limiter clamps it
    prm.monoBass.frequencyHz = 5.0f;                       // below the stage's own floor of 20 Hz
    chain.setParams (prm);
    Buf x = silence (nch, cfg.internalBlock);
    { auto p = planes (x); chain.process (p.data(), nch, cfg.internalBlock); }   // one quantum applies them

    const auto r = chain.resolved();
    ok (r.latencySamples == chain.latencySamples(), "resolved latency agrees with latencySamples()");
    ok (r.internalBlock == cfg.internalBlock, "resolved internal block");
    ok (r.oversampleFactor == 4, "resolved oversample factor");
    approx (r.limiterCeilingDbTp, 60.0, 1e-9, "an absurd ceiling comes back CLAMPED, not echoed");
    approx ((double) r.monoBass.frequencyHz, 20.0, 1e-4, "a mono-bass corner below the floor comes back clamped");
    ok (r.limiterLookahead == 48, "resolved limiter lookahead in samples (1 ms at 48 kHz)");
    ok (r.compressorLookahead == 48, "resolved compressor lookahead in samples");

    // The identity the whole latency contract rests on, stated in one line.
    ok (r.latencySamples == r.internalBlock + r.compressorLookahead + r.clipperLatency + r.limiterLatency,
        "latency == internal quantum + the sum of the present stages' own reported latencies");
}

//==============================================================================
static void testMonoBassParams()
{
    group ("stereo::MonoBassParams — the added type is exactly the three setters");

    // The stage grew a parameter type in this change; the point of the type is that it CANNOT drift
    // from the setters, so that is what gets checked: same object, bit-identical audio.
    const int n = 4000;
    std::vector<float> l1 ((std::size_t) n), r1 ((std::size_t) n), l2, r2;
    for (int i = 0; i < n; ++i)
    {
        l1[(std::size_t) i] = 0.6f * (float) std::sin (2.0 * kPi * 80.0 * i / 48000.0)
                            + 0.3f * (float) std::sin (2.0 * kPi * 900.0 * i / 48000.0);
        r1[(std::size_t) i] = -0.55f * l1[(std::size_t) i];
    }
    l2 = l1; r2 = r1;

    stereo::MonoBass a, b;
    a.prepare (48000.0, n, 2);
    b.prepare (48000.0, n, 2);
    a.setEnabled (true); a.setFrequency (137.0f); a.setLowWidth (0.25f);
    b.setParams ({ true, 137.0f, 0.25f });
    float* pa[2] { l1.data(), r1.data() };
    float* pb[2] { l2.data(), r2.data() };
    a.process (pa, 2, n);
    b.process (pb, 2, n);
    long long bad = 0;
    for (int i = 0; i < n; ++i)
        if (bits (l1[(std::size_t) i]) != bits (l2[(std::size_t) i]) || bits (r1[(std::size_t) i]) != bits (r2[(std::size_t) i])) ++bad;
    ok (bad == 0, "setParams() is bit-identical to the three setters");

    const auto rp = b.params();
    ok (rp.enabled && core::exactlyEqual (rp.frequencyHz, 137.0f) && core::exactlyEqual (rp.lowWidth, 0.25f),
        "params() reads back what was applied");
    b.setFrequency (2.0f);                                 // below kMinFreq
    approx ((double) b.params().frequencyHz, 20.0, 1e-4, "params() reports the CLAMPED frequency, not the request");
    b.setLowWidth (std::nan (""));
    approx ((double) b.params().lowWidth, 0.25, 1e-9, "a non-finite width is rejected, the last good value stands");
}

//==============================================================================
int main()
{
    std::printf ("felitronics::mastering — chain acceptance\n");
    testBlockInvariance();
    testBypassNull();
    testFlushIsProcessOfZeros();
    testRendererContract();
    testSequences();
    testPoisonGate();
    testMutationGaps();
    testRtSafety();
    testResolvedReadback();
    testMonoBassParams();
    return felitronics::test::report();
}
