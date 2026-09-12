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
#include <felitronics/oversampling/PolyphaseOversampler.h>
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

static std::atomic<long long> g_allocs { 0 };
// BYTES too, and EVERY FORM OF `new` (P41). A counter that overrides only `operator new(size_t)` and
// `new[]` is blind to the OVER-ALIGNED form — and `eq::EqEngine` has `alignof` 64, so 331 KiB of what a
// chain asks for, the single largest request it makes, went past such a counter unseen. That is the
// blindness P52 names; a suite that holds a published budget against a counter cannot have it.
//
// The bytes are the CONTAINER's, not what reached the allocator: MSVC's STL on x86/x64 asks for
// `sizeof(void*) + 31` more (one word more under _DEBUG) on every block of 4096 or more — its own
// alignment, which a REQUESTED-bytes budget leaves to the caller — and this takes that back off. The rule
// is calibrated on literal sizes below, so a check fails if this is not the STL it describes rather than
// a byte-for-byte comparison failing with no explanation. NEITHER FILE MODELS ITERATOR DEBUGGING (MSVC's
// `_DEBUG` adds a container-proxy allocation per container, which is not padding and cannot be subtracted
// off a block): under it the calibration fails by name, which is the point of having one.
#if defined(_MSVC_STL_VERSION) && (defined(_M_IX86) || defined(_M_X64))
#  if defined(_DEBUG)
static constexpr std::size_t kStlBigPad = 2 * sizeof (void*) + 31;
#  else
static constexpr std::size_t kStlBigPad = sizeof (void*) + 31;
#  endif
#else
static constexpr std::size_t kStlBigPad = 0;
#endif
static constexpr std::size_t kStlBigBlock = 4096;
static std::atomic<long long> g_bytes { 0 };
static long long containerBytes (std::size_t s) noexcept
{
    return (long long) (kStlBigPad != 0 && s >= kStlBigBlock + kStlBigPad ? s - kStlBigPad : s);
}
static void* countedNew (std::size_t s)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    g_bytes.fetch_add (containerBytes (s), std::memory_order_relaxed);
    return std::malloc (s ? s : 1);
}
// The over-aligned forms go through the platform's aligned allocator; `std::free` is not valid for it, so
// the matching deletes below use the aligned release. Everything is counted the same way.
static void* countedAlignedNew (std::size_t s, std::size_t a)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    // RAW, not `containerBytes`. That correction exists to undo what MSVC's STL adds ON TOP of a container's
    // own request; an over-aligned `new` in this tree is an OBJECT (`eq::EqEngine`, alignof 64) whose size is
    // exactly `sizeof`, and taking the correction off it subtracts bytes nobody ever added. Caught by the
    // `win` row before CI: the engine's 339 072 read as 339 033, and with it 88 of the chain's budget rows.
    g_bytes.fetch_add ((long long) s, std::memory_order_relaxed);
#if defined(_MSC_VER)
    return _aligned_malloc (s ? s : 1, a);
#else
    const std::size_t al = a < sizeof (void*) ? sizeof (void*) : a;
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
void* operator new      (std::size_t s) { return countedNew (s); }
void* operator new[]    (std::size_t s) { return countedNew (s); }
void* operator new      (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, (std::size_t) a); }
void* operator new[]    (std::size_t s, std::align_val_t a) { return countedAlignedNew (s, (std::size_t) a); }
void* operator new      (std::size_t s, const std::nothrow_t&) noexcept { return countedNew (s); }
void* operator new[]    (std::size_t s, const std::nothrow_t&) noexcept { return countedNew (s); }
void* operator new      (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return countedAlignedNew (s, (std::size_t) a); }
void* operator new[]    (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return countedAlignedNew (s, (std::size_t) a); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete   (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t) noexcept { alignedFree (p); }

// One vector's allocation, as the counter sees it — written through `volatile`, so the optimizer cannot
// remove an allocation nobody observes (the lesson P41 F6 paid for).
static long long vectorRequest (std::size_t n)
{
    const long long before = g_bytes.load();
    {
        std::vector<char> v;
        v.assign (n, 0);
        volatile char* sink = v.data();
        sink[0] = 1;
    }
    return g_bytes.load() - before;
}

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

// THE TOPOLOGY AXIS of the P41 budget matrix. Every optional stage is ABSENT on some row and PRESENT on
// others, which the first version of this table did not do: it moved only the clipper, so `eq`,
// `compressor`, `limiter` and `dither` were on in all 48 rows and `monoBass` in none. A budget that added
// the EQ engine's 331 KiB unconditionally — for a chain that never builds one — was green on every row of
// that matrix (found by the diverse-testing round, as a mutation that survived the whole suite).
constexpr int kTopologies = 9;
mastering::MasteringChainConfig topologyFor (int topo)
{
    mastering::MasteringChainConfig c;
    switch (topo)
    {
        case 0: break;                                            // the default
        case 1: c.clipper = true; break;
        case 2: c.sidechainHpfHz = 80.0; break;
        case 3: c.eq = false; break;                              // no engine at all — 331 KiB of budget
        case 4: c.compressor = false; break;                      // no delay bank
        case 5: c.limiter = false; break;                         // no oversampler, no aligner, tap factor 1
        case 6: c.dither = false; break;
        case 7: c.monoBass = true; break;                         // stereo only: refused at widths 1 and 16
        default:                                                  // the legal maximum
            c.clipper = true; c.sidechainHpfHz = 80.0; c.internalBlock = 8192;
            c.oversampleFactor = 16; c.tapsPerPhase = 1024;
            c.compressorLookaheadMs = 250.0; c.limiterLookaheadMs = 20.0;
            break;
    }
    return c;
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
        for (int b : partPrimes (n)) { auto p = planes (y, off); felitronics::test::run (chain.process (p.data(), nch, b)); off += b; }
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
    felitronics::test::run (r.prepare (nch, 997));
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
    { auto p = planes (ya); felitronics::test::run (a.process (p.data(), nch, n)); }
    Buf ta = silence (nch, D);
    { auto p = planes (ta); ok (a.flush (p.data(), nch, D) == D, "flush() writes exactly latencySamples() frames"); }

    // (b) stream x followed by D zeros, no flush() at all.
    Buf xb = x;
    for (int c = 0; c < nch; ++c) xb[(std::size_t) c].resize ((std::size_t) (n + D), 0.0f);
    { auto p = planes (xb); felitronics::test::run (b.process (p.data(), nch, n + D)); }

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
        felitronics::test::run (r1.prepare (nch, 4096));
        felitronics::test::run (r2.prepare (nch, 1));
        auto ip = planes (const_cast<Buf&> (x));
        { auto op = planes (ref); ok (r1.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render at blockSize 4096"); }
        { auto op = planes (got); ok (r2.render (chain, (const float* const*) ip.data(), op.data(), nch, n), "render at blockSize 1"); }
        ok (bitEqual (ref, got), "renderer blockSize does not change the result");

        mastering::OfflineRenderer r3;
        felitronics::test::run (r3.prepare (nch, 65536));
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
        felitronics::test::run (r.prepare (nch, 512));
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
        felitronics::test::run (r.prepare (nch, 333));
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
        felitronics::test::run (r.prepare (nch, 640));
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
        { auto p = planes (yb); felitronics::test::run (b.process (p.data(), nch, n)); }
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
        { auto p = planes (warm); felitronics::test::run (a.process (p.data(), nch, n)); }   // dirty the state
        a.reset();
        Buf ya = x, yf = x;
        { auto p = planes (ya); felitronics::test::run (a.process (p.data(), nch, n)); }
        { auto p = planes (yf); felitronics::test::run (fresh.process (p.data(), nch, n)); }
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

        Buf a = loud; { auto pl = planes (a); felitronics::test::run (chain.process (pl.data(), nch, 4096)); }          // active, loud
        p.bypassLimiter = true; chain.setParams (p);
        Buf b = loud; { auto pl = planes (b); felitronics::test::run (chain.process (pl.data(), nch, 4096)); }          // bypassed, loud
        for (int k = 0; k < 4; ++k)                                                             // drain everything
        { Buf z = silence (nch, 4096); auto pl = planes (z); felitronics::test::run (chain.process (pl.data(), nch, 4096)); }
        p.bypassLimiter = false; chain.setParams (p);
        Buf out = silence (nch, 4096);
        { auto pl = planes (out); felitronics::test::run (chain.process (pl.data(), nch, 4096)); }                       // active again, SILENCE in

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
        { auto pl = planes (x); felitronics::test::run (chain.process (pl.data(), nch, T)); }          // active
        p.bypassClipper = p.bypassLimiter = true;
        chain.setParams (p);                                                  // takes effect at sample T
        { auto pl = planes (x, T); felitronics::test::run (chain.process (pl.data(), nch, n - T)); }   // bypassed

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
        { Buf y = x; auto pl = planes (y); felitronics::test::run (c1.process (pl.data(), nch, 1024)); }
        c1.setParams (B);                                  // starts the ramp
        { Buf y = x; auto pl = planes (y); felitronics::test::run (c1.process (pl.data(), nch, 512)); }   // stop ~10 ms into it
        c1.reset();
        Buf y1 = x; { auto pl = planes (y1); felitronics::test::run (c1.process (pl.data(), nch, n)); }

        c2.setParams (B);
        Buf y2 = x; { auto pl = planes (y2); felitronics::test::run (c2.process (pl.data(), nch, n)); }
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
        { auto pl = planes (padded); felitronics::test::run (a.process (pl.data(), nch, n + D)); }

        Buf out = silence (nch, n);
        mastering::OfflineRenderer r;
        felitronics::test::run (r.prepare (nch, 512));
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
    felitronics::test::run (chain.process (px, nch, 1024));
    chain.reset();

    const long long before = g_allocs.load();
    setPlanes (0);
    felitronics::test::run (chain.process (px, nch, n));                            // one whole-programme call
    {
        int off = 0;
        for (int b : part) { setPlanes (off); felitronics::test::run (chain.process (px, nch, b)); off += b; }
    }
    prm.bypassLimiter = true;  chain.setParams (prm);      // a bypass toggle, mid-stream
    setPlanes (0); felitronics::test::run (chain.process (px, nch, 4096));
    prm.bypassLimiter = false; chain.setParams (prm);
    setPlanes (0); felitronics::test::run (chain.process (px, nch, 4096));
    const bool monoRefused = ! chain.process (pm, 1, 64);   // asserted after the snapshot — test::ok allocates
    chain.flush (pt, nch, (int) tail[0].size());
    chain.reset();
    (void) chain.resolved();
    const long long after = g_allocs.load();
    test::ok (monoRefused, "a mono call on a stereo chain is refused");

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
    { auto p = planes (x); felitronics::test::run (chain.process (p.data(), nch, cfg.internalBlock)); }   // one quantum applies them

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
// A parameter set written BEFORE prepare() has to survive it, and the readback straight after prepare()
// has to describe the prepared chain rather than the previous one. Both were false: `prepare()` ended
// with `pendingParams_ = params_`, which replaced the caller's pending write with the last APPLIED set.
static void testParamsWrittenBeforePrepare()
{
    group ("setParams() before prepare() is KEPT, and the readback is not stale");

    const int nch = 1;
    auto cfg = fullConfig();
    cfg.eq = false; cfg.monoBass = false; cfg.compressor = false; cfg.clipper = false;
    cfg.limiter = false; cfg.dither = false;               // gain nodes only: the effect is unmistakable

    mastering::MasteringChain chain;
    mastering::MasteringChainParams prm;
    prm.inputGainDb = 12.0;
    chain.setParams (prm);                                 // BEFORE prepare
    ok (chain.prepare (48000.0, nch, cfg), "prepare after setParams");
    approx (chain.params().inputGainDb, 12.0, 0.0, "params() reports the set written before prepare()");

    const int n = 4 * cfg.internalBlock;
    Buf x ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.1f));
    { auto p = planes (x); felitronics::test::run (chain.process (p.data(), nch, n)); }
    // PRECONDITION: past the chain's latency, so the sample examined is a processed one and not priming.
    const int probe = cfg.internalBlock + 44;
    ok (probe < n && probe > chain.latencySamples() - 1, "precondition: the probed sample is past the priming");
    approx ((double) x[0][(std::size_t) probe], 0.1 * std::pow (10.0, 12.0 / 20.0), 1e-6,
            "the audio carries the gain that was set before prepare()");

    // And the same value applied the other way round must give the same audio — the two orders are one
    // behaviour, which is what "kept" means.
    mastering::MasteringChain c2;
    ok (c2.prepare (48000.0, nch, cfg), "prepare then setParams");
    c2.setParams (prm);
    Buf y ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.1f));
    { auto p = planes (y); felitronics::test::run (c2.process (p.data(), nch, n)); }
    long long bad = 0;
    for (int i = 0; i < n; ++i) if (bits (x[0][(std::size_t) i]) != bits (y[0][(std::size_t) i])) ++bad;
    ok (bad == 0, "configure-then-prepare is bit-identical to prepare-then-configure");

    // A REFUSED prepare() must not adopt the pending set either — law 11(b) says a refused prepare
    // leaves the object unusable, and an unusable object reporting new parameters is a lie about both.
    mastering::MasteringChain c3;
    mastering::MasteringChainParams other; other.inputGainDb = -30.0;
    c3.setParams (other);
    ok (! c3.prepare (-1.0, nch, cfg), "a bad sample rate is refused");
    ok (! c3.isPrepared(), "a refused prepare leaves the chain unprepared");
    ok (c3.latencySamples() == 0, "a refused prepare reports no latency");

    // RE-PREPARE with a write still pending. The keep-the-pending-set rule has two cases and only one
    // of them is "a fresh object": a chain that has already run, is written to, and is then prepared
    // again for a new stream must carry that write into the new stream rather than the last one it
    // applied. It is the same line that used to discard it.
    mastering::MasteringChain c4;
    ok (c4.prepare (48000.0, nch, cfg), "re-prepare: first preparation");
    { Buf z ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.05f));
      auto p = planes (z); felitronics::test::run (c4.process (p.data(), nch, n)); }
    mastering::MasteringChainParams later; later.inputGainDb = -6.0;
    c4.setParams (later);                                  // pending, never applied
    ok (c4.prepare (44100.0, nch, cfg), "re-prepare: at a new rate, with a write still pending");
    approx (c4.params().inputGainDb, -6.0, 0.0, "re-prepare KEEPS the pending write, not the applied one");
    approx (c4.sampleRate(), 44100.0, 0.0, "re-prepare reports the new rate");
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
    felitronics::test::run (a.prepare (48000.0, n, 2));
    felitronics::test::run (b.prepare (48000.0, n, 2));
    a.setEnabled (true); a.setFrequency (137.0f); a.setLowWidth (0.25f);
    b.setParams ({ true, 137.0f, 0.25f });
    float* pa[2] { l1.data(), r1.data() };
    float* pb[2] { l2.data(), r2.data() };
    felitronics::test::run (a.process (pa, 2, n));
    felitronics::test::run (b.process (pb, 2, n));
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

//==============================================================================
// THE GATE COUNTS WHAT IT SUBSTITUTES. This belongs HERE and not only in the C-ABI suite: a consumer
// that builds the core as a subproject never builds `tools/`, so a test that lives only there is a test
// that consumer does not have. The counter is also the thing that makes the gate's silence visible at
// all, and its shape is load-bearing for a different reason — see the comment at the gate itself, which
// carries the measurement that made an earlier form of it x2.61 slower.
static void testGateCountsItsSubstitutions()
{
    group ("the input gate COUNTS the samples it substitutes");

    const int nch = 2, n = 4 * 256;
    auto cfg = fullConfig();
    cfg.eq = false; cfg.monoBass = false; cfg.compressor = false;
    cfg.clipper = false; cfg.limiter = false; cfg.dither = false;

    mastering::MasteringChain chain;
    ok (chain.prepare (48000.0, nch, cfg), "prepared");
    ok (chain.nonFiniteInputSamples() == 0, "a fresh chain has counted nothing");

    Buf x ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.25f));
    x[0][10] = std::numeric_limits<float>::quiet_NaN();
    x[1][20] = std::numeric_limits<float>::infinity();
    x[1][30] = -std::numeric_limits<float>::infinity();
    { auto p = planes (x); felitronics::test::run (chain.process (p.data(), nch, n)); }
    ok (chain.nonFiniteInputSamples() == 3, "three substitutions, counted exactly");

    bool finite = true;
    for (const auto& ch : x) for (float v : ch) if (! std::isfinite (v)) finite = false;
    ok (finite, "and nothing non-finite left the chain");

    // A FINITE sample outside the gate's range is CLAMPED and deliberately not counted here: the two
    // are different events and one number cannot mean both.
    Buf y ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.25f));
    y[0][40] = 1.0e9f;
    { auto p = planes (y); felitronics::test::run (chain.process (p.data(), nch, n)); }
    ok (chain.nonFiniteInputSamples() == 3, "a clamped-but-finite sample does not move the counter");
    // The chain delays by its internal quantum, so input sample 40 leaves at output sample 40 + K.
    // Reading index 40 would have read the FIFO's previous contents and made this precondition a
    // tautology — which is exactly how a blind check looks from the inside.
    ok (y[0][40 + cfg.internalBlock] == 1.0e6f,
        "PRECONDITION: it really was clamped, so the check above is live");

    chain.reset();
    ok (chain.nonFiniteInputSamples() == 0, "reset() clears it — the count describes ONE stream");
}


//==============================================================================
// P41 — THE DEMAND IS THE ALLOCATION (law 11d). The chain publishes what preparing it will ask the heap
// for, and this holds the published number against a counter that sees every form of `new`, byte for byte,
// over rates x widths x topologies. Every delta is read into a LOCAL before its check: a call's arguments
// are evaluated in an unspecified order and a message string allocates (gcc builds it first — measured).
void testDemand()
{
    group ("P41: the chain's demand is exactly what preparing it allocates");

    // THE COUNTER'S OWN RULE FIRST, on literal sizes rather than on its own constants: around the STL's
    // big-block threshold and far above it, a vector counts as exactly the bytes it asked for. If this
    // STL pads differently from the rule above, THIS fails with its name on it instead of a byte-for-byte
    // comparison failing with no explanation.
    for (const std::size_t nb : { std::size_t { 4095 }, std::size_t { 4096 }, std::size_t { 4097 }, std::size_t { 1 } << 20 })
    {
        const long long got = vectorRequest (nb);
        ok (got == (long long) nb, "the byte counter counts a " + std::to_string (nb) + "-byte vector as "
                                   + std::to_string (nb) + " bytes (read " + std::to_string (got) + ")");
    }
    // AND THAT IT SEES AN OVER-ALIGNED `new` AT ALL — the form `eq::EqEngine` is built through. A counter
    // blind to it reads the default geometry's largest single request as zero and every budget below passes,
    // which is the shape of P52 and is why this check is a PRECONDITION rather than a nicety.
    {
        const long long before = g_bytes.load();
        {
            auto e = std::make_unique<eq::EqEngine>();
            volatile const void* sink = e.get();
            (void) sink;
        }
        const long long got = g_bytes.load() - before;
        ok (got == (long long) eq::EqEngine::objectBytes(),
            "the counter sees the over-aligned `new` the EQ engine is built through: "
            + std::to_string (got) + " B");
    }

    // THE MATRIX. Four rates, three widths, four topologies — the default, the clipper, the key filter,
    // and the legal maximum. `create` on the C ABI carries the same matrix; this is the core's own side.
    const double rates[] = { 44100.0, 48000.0, 96000.0, 192000.0 };
    const int    widths[] = { 1, 2, 16 };
    int rows = 0, bad = 0, badLat = 0, badAdmit = 0, badReprep = 0;
    long long worst = 0;
    for (const double fs : rates)
        for (const int nch : widths)
            for (int topo = 0; topo < kTopologies; ++topo)
            {
                const mastering::MasteringChainConfig cfg = topologyFor (topo);
                ++rows;
                const std::uint64_t budget = mastering::MasteringChain::prepareBytes (fs, nch, cfg);
                const bool admitted = mastering::MasteringChain::admits (fs, nch, cfg);

                auto chain = std::make_unique<mastering::MasteringChain>();   // FRESH: the budget is a fresh object's
                const long long before = g_bytes.load();
                const bool prepared = chain->prepare (fs, nch, cfg);
                const long long got = g_bytes.load() - before;

                // ADMITS IS PREPARE'S OWN VERDICT, reached without allocating. Not "agrees usually".
                if (admitted != prepared) ++badAdmit;
                if (! prepared) continue;
                if (got != (long long) budget) { ++bad; if (got > worst) worst = got; }
                // The latencies the budget used to size the dry aligners are the ones the prepared stages
                // report. A budget that sized an aligner from a second derivation of the latency would
                // pass every byte check on the default topology and fail on a moved lookahead.
                {
                    mastering::MasteringChain::Storage st;
                    if (! mastering::MasteringChain::storageFor (fs, nch, cfg, st)) ++badLat;
                    if (st.latencySamples != chain->latencySamples()) ++badLat;
                    if (st.tapOversampleFactor != chain->tapOversampleFactor()) ++badLat;
                }
                // RE-PREPARING AT THE SAME GEOMETRY COSTS NOTHING — published and measured.
                const std::uint64_t again = chain->reprepareBytes (fs, nch, cfg);
                const long long b2 = g_bytes.load();
                const bool ok2 = chain->prepare (fs, nch, cfg);
                const long long got2 = g_bytes.load() - b2;
                if (! ok2 || again != 0u || got2 != 0) ++badReprep;
            }
    ok (rows == 4 * 3 * kTopologies, "PRECONDITION: the matrix is 4 rates x 3 widths x "
        + std::to_string (kTopologies) + " topologies (" + std::to_string (rows) + " rows)");
    ok (badAdmit == 0, "admits() is prepare()'s own verdict on every row, reached without allocating ("
                       + std::to_string (badAdmit) + " disagreements)");
    ok (bad == 0, "prepareBytes() is what a fresh prepare() allocates, byte for byte, on every row ("
                  + std::to_string (bad) + " rows off, worst " + std::to_string (worst) + " B)");
    ok (badLat == 0, "the latencies and the tap factor the budget computed are the ones the prepared chain "
                     "reports (" + std::to_string (badLat) + " off)");
    ok (badReprep == 0, "re-preparing at the same geometry asks for nothing, and says so ("
                        + std::to_string (badReprep) + " rows off)");

    // AND A RE-PREPARATION THAT GROWS IS NOT FREE. Every row above re-prepares at the geometry the chain
    // already holds, so "the budget is 0" and "the budget is always 0" pass the same rows — a mutation that
    // returned 0 unconditionally survived the whole suite until this. A bound, not an equality: a chain
    // that grows keeps the containers that still fit, so it asks for LESS than a fresh one, and the header
    // says the number is an upper bound in exactly that case.
    {
        mastering::MasteringChainConfig small;
        mastering::MasteringChainConfig big = small;
        big.internalBlock = 4096; big.clipper = true; big.oversampleFactor = 8; big.tapsPerPhase = 256;
        auto chain = std::make_unique<mastering::MasteringChain>();
        ok (chain->prepare (48000.0, 2, small), "PRECONDITION: a chain at the small geometry");
        const std::uint64_t bound = chain->reprepareBytes (96000.0, 2, big);
        const long long before = g_bytes.load();
        const bool grew = chain->prepare (96000.0, 2, big);
        const long long got = g_bytes.load() - before;
        ok (grew, "PRECONDITION: it re-prepares at the bigger one");
        ok (bound > 0u, "a re-preparation that GROWS is not published as free (" + std::to_string (bound) + " B)");
        ok (got > 0 && (std::uint64_t) got <= bound,
            "and what it asks for is inside the published bound (asked " + std::to_string (got)
            + " B, bound " + std::to_string (bound) + " B)");
    }

    // THE TWO WAYS THE RE-PREPARATION BOUND WAS FOUND TO LEAK (the code-review round), each with the input
    // that found it. Both are about a budget that answered for storage it could not see.
    {
        // (a) A MOVED-FROM chain keeps `prepared_` and its scalars and has given its buffers away. The
        // budget used to recompute "what I have" from those scalars, say it fits, and publish 0 for a
        // preparation that really allocated.
        mastering::MasteringChainConfig bare;
        bare.internalBlock = 64;
        bare.eq = bare.compressor = bare.limiter = bare.dither = false;
        mastering::MasteringChain a;
        ok (a.prepare (48000.0, 2, bare), "PRECONDITION: a prepared chain");
        mastering::MasteringChain b (std::move (a));
        const std::uint64_t bound = a.reprepareBytes (48000.0, 2, bare);     // NOLINT: the point is the moved-from state
        const long long before = g_bytes.load();
        const bool again = a.prepare (48000.0, 2, bare);
        const long long got = g_bytes.load() - before;
        ok (again && got > 0, "PRECONDITION: preparing it again really does allocate (" + std::to_string (got) + " B)");
        ok ((std::uint64_t) got <= bound, "a chain that has given its storage away does not publish 0 (bound "
                                          + std::to_string (bound) + " B)");
    }
    {
        // (b) GEOMETRIC GROWTH. Widening the limiter's per-channel scratch bank by one used to `resize`
        // past the capacity, which doubles — so the call held two buffer objects more than the number paid
        // for, and an "upper bound" was not one. The witness is the code-review round's own: a 3-channel
        // chain grown to 4, everything else moving with it.
        mastering::MasteringChainConfig c1;
        c1.eq = c1.compressor = c1.dither = false;
        c1.internalBlock = 16; c1.oversampleFactor = 2; c1.tapsPerPhase = 4; c1.limiterLookaheadMs = 0.0;
        mastering::MasteringChainConfig c2 = c1;
        c2.internalBlock = 17; c2.oversampleFactor = 3; c2.tapsPerPhase = 5;
        auto chain = std::make_unique<mastering::MasteringChain>();
        ok (chain->prepare (100.0, 3, c1), "PRECONDITION: a narrow chain at a low rate");
        const std::uint64_t bound = chain->reprepareBytes (200.0, 4, c2);
        const long long before = g_bytes.load();
        const bool grew = chain->prepare (200.0, 4, c2);
        const long long got = g_bytes.load() - before;
        ok (grew, "PRECONDITION: and it grows in every dimension at once");
        // THE BOUND IS ON WHAT THE CHAIN ASKS ITS CONTAINERS FOR, not on what they then ask the allocator.
        // A container that has to grow applies its OWN growth policy, and that is the one thing law 11d
        // hands to the caller in as many words ("the runtime's growth step ... the caller's margin"). The
        // `win` row is the witness and the reason this is not an inequality: MSVC's `vector::assign` past
        // the capacity grows by half, so a 48-float buffer asked to hold 68 asks the allocator for 72 —
        // twice over in this very sequence, 32 B past a bound libc++ meets exactly. What IS checkable, and
        // what the mutation stand needs, is that a growing re-preparation is not published as free.
        // NARROWLY, about THIS chain: "a growing re-preparation is never free" is false in general, and the
        // fix round produced the sequence — a chain prepared at K = 8192, re-prepared at K = 8, and then
        // asked for K = 64 asks for nothing, because it never gave the big storage up. What holds is that a
        // chain which has never held more than it holds now, and is asked for more, pays.
        ok (bound > 0u && got > 0, "a re-preparation of a chain that never held more, growing, is not free — "
            "and is not published as free (asked " + std::to_string (got) + " B, bound "
            + std::to_string (bound) + " B)");
    }
    {
        // AND ONE WHERE THE CHAIN'S OWN BUFFERS DO NOT MOVE. The row above grows the quantum, so the FIFO
        // sentinel alone answers it and `fitsWithin` is never consulted — measured on the stand: a mutant
        // that returned 0 unconditionally survived, because the sentinel short-circuited it every time.
        // Here the rate, the width and the quantum are unchanged (FIFO identical, `keyBuf_` absent) and only
        // the limiter's topology grows, so nothing but `fitsWithin` can see it.
        mastering::MasteringChainConfig c1;
        mastering::MasteringChainConfig c2 = c1;
        c2.oversampleFactor = 16; c2.tapsPerPhase = 1024; c2.limiterLookaheadMs = 20.0;
        auto chain = std::make_unique<mastering::MasteringChain>();
        ok (chain->prepare (48000.0, 2, c1), "PRECONDITION: a chain at the default topology");
        mastering::MasteringChain::Storage a {}, b {};
        (void) mastering::MasteringChain::storageFor (48000.0, 2, c1, a);
        (void) mastering::MasteringChain::storageFor (48000.0, 2, c2, b);
        ok (a.fifo == b.fifo && a.keyBuf == b.keyBuf, "PRECONDITION: the chain's own buffers do not move");
        const std::uint64_t bound = chain->reprepareBytes (48000.0, 2, c2);
        const long long before = g_bytes.load();
        const bool grew = chain->prepare (48000.0, 2, c2);
        const long long got = g_bytes.load() - before;
        ok (grew && got > 0, "PRECONDITION: growing only the limiter still allocates (" + std::to_string (got) + " B)");
        ok (bound > 0u, "and a stage that grows behind an unchanged FIFO is not published as free");
    }

    // THE STAGES' OWN BUDGETS, ASKED DIRECTLY, at arguments no chain can reach. Three mutations survived
    // the matrix above for this reason alone: the chain caps its quantum at 8192 (so the limiter's own
    // 1 Mi-sample block cap is never exercised), refuses a width past `core::kMaxChannels` (so the
    // oversampler's channel CLAMP is never exercised), and never asks a dry aligner for a capacity under 2.
    // A budget is only pinned where its own clamps can be reached.
    {
        int off = 0;
        auto measure = [&off] (const char* /*what*/, std::uint64_t budget, auto&& build)
        {
            const long long before = g_bytes.load();
            const bool okPrep = build();
            const long long got = g_bytes.load() - before;
            if (! okPrep || got != (long long) budget) ++off;
        };
        {   // the limiter's block cap: a whole-file maxBlock is a normal thing for an offline caller to pass
            limiter::TruePeakLimiterConfig lc;
            limiter::TruePeakLimiter::Storage st;
            const bool okSt = limiter::TruePeakLimiter::storageFor (48000.0, (1 << 20) + 7, 2, lc, st);
            if (! okSt) ++off;
            auto l = std::make_unique<limiter::TruePeakLimiter>();
            measure ("limiter, block past the cap", st.bytes(),
                     [&] { return l->prepare (48000.0, (1 << 20) + 7, 2, lc); });
        }
        {   // the oversampler's channel clamp — law 11(b)'s unfinished application, tracked as P55. The
            // budget must model what prepare() DOES, not what it ought to do.
            oversampling::PolyphaseOversampler::Storage st;
            const bool okSt = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels + 1, 64, st);
            if (! okSt) ++off;
            auto o = std::make_unique<oversampling::PolyphaseOversampler>();
            measure ("oversampler, width past the maximum", st.bytes(),
                     [&] { return o->prepare (4, core::kMaxChannels + 1, 64); });
        }
        {   // the aligner's floors: a capacity under 2 and a block under 1 are RAISED, not refused, and the
            // budget has to raise them too. Two channels, so the raised counts still exceed the seed the
            // aligner's constructor took (a one-channel aligner at the floor is already the seed, and then
            // preparing it asks for nothing — which would make this check pass whatever the floors did).
            const core::DryAligner::Storage st = core::DryAligner::storageFor (2, 1, 0);
            auto a = std::make_unique<core::DryAligner>();
            const long long before = g_bytes.load();
            a->prepare (2, 1, 0);
            const long long got = g_bytes.load() - before;
            if (got != (long long) st.bytes()) ++off;
        }
        ok (off == 0, "each stage's own budget is what preparing it directly allocates, at the arguments its "
                      "OWN clamps live at (" + std::to_string (off) + " off)");
    }

    // AND THE CLAMPS ARE PINNED AGAINST A PROPERTY, NOT AGAINST THE ALLOCATION — which is the one thing the
    // byte-for-byte checks above CANNOT do. `prepare()` now sizes itself THROUGH `storageFor`, so a mutation
    // inside that function moves the budget and the allocation together and the comparison stays green: the
    // stand proved it, with three clamps removed one at a time and the whole suite still passing. What a
    // clamp needs is an independent statement about itself, and idempotence past the bound is the natural
    // one: beyond the clamp the answer must stop moving.
    {
        limiter::TruePeakLimiterConfig lc;
        limiter::TruePeakLimiter::Storage cap1 {}, cap2 {};
        const bool a1 = limiter::TruePeakLimiter::storageFor (48000.0, 1 << 20, 2, lc, cap1);
        const bool a2 = limiter::TruePeakLimiter::storageFor (48000.0, 1 << 24, 2, lc, cap2);
        ok (a1 && a2 && cap1.bytes() == cap2.bytes() && cap1.osBufSamples == cap2.osBufSamples,
            "the limiter's scratch stops growing at its block cap: a 16x bigger block is the same budget ("
            + std::to_string (cap1.bytes()) + " B vs " + std::to_string (cap2.bytes()) + " B)");

        oversampling::PolyphaseOversampler::Storage w1 {}, w2 {};
        const bool b1 = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels, 64, w1);
        const bool b2 = oversampling::PolyphaseOversampler::storageFor (4, core::kMaxChannels + 8, 64, w2);
        ok (b1 && b2 && w1.bytes() == w2.bytes(),
            "the oversampler's width stops at kMaxChannels, which is what its prepare() CLAMPS to (P55: that "
            "clamp should be a refusal, and until it is, the budget has to describe the clamp)");

        // The two smallest published clamps, which no chain reaches and which therefore had no gate at all
        // until the fix round asked for them: a negative delay is a ring of one slot, and a window under one
        // is a window of one. Both are what their `prepare()` does, and both are stated in `storageFor`.
        ok (core::DelayLine::storageFor (-1).samples == 1 && core::DelayLine::storageFor (-1000).samples == 1
            && core::DelayLine::storageFor (0).samples == 1,
            "a negative delay budgets the one slot prepare() gives it");
        ok (limiter::detail::SlidingMax::storageFor (0).entries == 1
            && limiter::detail::SlidingMax::storageFor (-5).entries == 1,
            "a window under one budgets the one entry prepare() gives it");

        const core::DryAligner::Storage f0 = core::DryAligner::storageFor (2, 0, 0);
        const core::DryAligner::Storage f1 = core::DryAligner::storageFor (2, 1, 2);
        ok (f0.bytes() == f1.bytes() && f0.ring == 4 && f0.scratch == 2,
            "the aligner's floors are floors: a capacity of 0 budgets the 2 slots it will be given, and a "
            "block of 0 the 1 sample (" + std::to_string (f0.bytes()) + " B)");
    }

    // THE DELAY BANK LEAVES WHAT THE `assign` IT REPLACED LEFT. The old form built every line fresh, so a
    // re-prepared bank had no tap; the helper re-uses the lines, and a tap that outlived a SHRINKING
    // preparation used to index before the start of the ring (the code-review round found it, and the
    // invariant it breaks predates this work: `prepare(8); setDelay(8); prepare(2)` on a bare DelayLine).
    {
        std::vector<core::DelayLine> bank;
        core::prepareDelayBank (bank, 2, 8);
        bank[0].setDelay (8);
        core::prepareDelayBank (bank, 1, 2);
        ok (bank.size() == 1 && bank[0].capacity() == 2 && bank[0].delay() == 0,
            "a shrunk bank has the capacity asked for and NO tap, as a freshly built one does");
        core::DelayLine d;
        d.prepare (8);
        d.setDelay (8);
        d.prepare (2);
        ok (d.delay() <= d.capacity(), "and a bare line's tap is re-clamped into the capacity it was re-prepared at");
        float last = 0.0f;
        for (int i = 0; i < 16; ++i) last = d.process ((float) (i + 1));
        ok (std::isfinite (last), "PRECONDITION: reading it afterwards stays inside the ring");
    }

    // A REFUSED PREPARATION ALLOCATES NOTHING — the whole point of `admits()` being ahead of the first
    // buffer. Each of these used to cost between 51 288 and 394 456 bytes on its way to `false` at this
    // geometry, and up to 1 668 312 at sixteen channels and an 8192-sample quantum.
    {
        struct Case { const char* name; double fs; int nch; mastering::MasteringChainConfig cfg; };
        std::vector<Case> cases;
        auto push = [&cases] (const char* n, double fs, int nch, mastering::MasteringChainConfig c)
        { cases.push_back (Case { n, fs, nch, c }); };
        mastering::MasteringChainConfig d;
        push ("a 20 Hz rate", 20.0, 2, d);
        push ("a rate of zero", 0.0, 2, d);
        push ("a NaN rate", std::numeric_limits<double>::quiet_NaN(), 2, d);
        { auto c = d; c.compressorLookaheadMs = 300.0; push ("a lookahead past the compressor's 250 ms", 48000.0, 2, c); }
        { auto c = d; c.tapsPerPhase = 2000;           push ("more taps than the oversampler designs", 48000.0, 2, c); }
        { auto c = d; c.oversampleFactor = 32;         push ("a factor past the limiter's 16", 48000.0, 2, c); }
        { auto c = d; c.monoBass = true;               push ("mono-bass on a mono chain", 48000.0, 1, c); }
        { auto c = d; c.internalBlock = 4;             push ("a quantum under the floor", 48000.0, 2, c); }
        int leaked = 0, admitted = 0;
        for (const Case& cs : cases)
        {
            auto chain = std::make_unique<mastering::MasteringChain> ();
            if (mastering::MasteringChain::admits (cs.fs, cs.nch, cs.cfg)) ++admitted;
            const long long before = g_bytes.load();
            const bool prepared = chain->prepare (cs.fs, cs.nch, cs.cfg);
            const long long got = g_bytes.load() - before;
            if (prepared || got != 0) ++leaked;
        }
        ok (admitted == 0, "PRECONDITION: every case above is one admits() refuses");
        ok (leaked == 0, "a refused prepare() allocates nothing, on all " + std::to_string (cases.size())
                         + " refusals (" + std::to_string (leaked) + " leaked)");
    }

    // THE STAGES' OWN STATICS, against prepared stages. This is what makes the chain's aligner sizing and
    // its budget the same arithmetic the stage runs rather than a copy of it.
    {
        int off = 0;
        for (const double fs : { 44100.0, 48000.0, 192000.0 })
            for (const double look : { 0.0, 1.0, 10.0, 100.0, 250.0 })
            {
                dynamics::Compressor c;
                const double maxLook = std::max (look, 1.0);
                if (! c.prepare (fs, 256, 2, maxLook)) { ++off; continue; }
                dynamics::CompressorParams cp; cp.lookaheadMs = look; c.setParams (cp);
                if (c.latencySamples() != dynamics::Compressor::latencyFor (fs, 256, 2, maxLook, look)) ++off;

                limiter::TruePeakLimiterConfig lc; lc.lookaheadMs = look > 20.0 ? 20.0 : look;
                limiter::TruePeakLimiter l;
                if (! l.prepare (fs, 256, 2, lc)) { ++off; continue; }
                if (l.latencySamples() != limiter::TruePeakLimiter::latencyFor (fs, 256, 2, lc)) ++off;
                if (l.oversampleFactor() != limiter::TruePeakLimiter::oversampleFactorFor (lc)) ++off;

                // EVERY FACTOR, the ones BELOW two included: the saturator's own domain has a second half
                // where the oversampler is not built at all and the latency is 0 whatever the taps say. The
                // chain never asks for it (its `admits` requires a factor of 2 or more), so nothing else in
                // this suite reaches it — measured on the stand, a `latencyFor` that ignored the factor
                // entirely and answered `tpp - 1` survived the whole suite.
                for (const int f : { -1, 0, 1, 2, 4, 16 })
                {
                    saturation::Saturator sat;
                    if (! sat.prepare (fs, 256, 2, f, 64)) { ++off; continue; }
                    if (sat.latencySamples() != saturation::Saturator::latencyFor (fs, 256, 2, f, 64)) ++off;
                }
            }
        ok (off == 0, "every stage's latencyFor() is the latency the prepared stage reports ("
                      + std::to_string (off) + " off)");
    }
}

//==============================================================================
// P41 — THE EQ ENGINE IS REUSED, AND THAT IS BIT-IDENTICAL. `prepare()` used to build a second 331 KiB
// engine before releasing the first; it now re-prepares the one it has. The argument that this is safe is
// the engine's own (`EqBand::reset()` is defined as "the state prepare() leaves it in", and the chain
// rewrites every band straight after), but an argument is not a proof: this renders the same programme
// through a chain prepared ONCE and through a chain prepared THREE times, with parameters written in
// between, and nulls them sample by sample. Before the change the second chain got a brand-new engine
// every time, so a difference here is exactly the reuse being wrong.
void testEqEngineReuseIsBitIdentical()
{
    group ("P41: re-preparing reuses the EQ engine and the audio does not move");
    constexpr int nch = 2, n = 4096;
    const Buf x = burstThenSilence (nch, n, 900);

    mastering::MasteringChainParams pa;          // a parameter set that really drives the EQ
    pa.eqBands[0].on = true;  pa.eqBands[0].lane (eq::Lane::Stereo).gainDb = 6.0;
    pa.eqBands[0].lane (eq::Lane::Stereo).freq = 120.0;
    pa.eqBands[3].on = true;  pa.eqBands[3].lane (eq::Lane::Stereo).gainDb = -4.5;
    pa.eqBands[3].lane (eq::Lane::Stereo).freq = 3200.0;
    pa.eqBands[5].on = true;  pa.eqBands[5].type = eq::FilterType::HighShelf;
    pa.eqBands[5].lane (eq::Lane::Stereo).gainDb = 2.5;
    pa.eqBands[5].lane (eq::Lane::Stereo).freq = 8000.0;

    // A DIFFERENT set for the chain that is re-prepared — and it MOVES THE FILTER TYPE on every band it
    // touches, because the type is the one parameter a re-preparation carries into the next one:
    // `EqBand::prepare()` seeds `lastType_` from the params the band is STILL HOLDING, and the first
    // `setParams()` after a preparation takes the uninitialised branch, which does not touch it. So the
    // `typeChanged` edge inside `updateCoeffs()` sees a different previous type in a re-used engine than in
    // a fresh one, and if that edge cleared anything `clearAudioState()` does not already clear, the two
    // would diverge. Without a type change in this fixture the null below would be blind to exactly that.
    mastering::MasteringChainParams pb = pa;
    pb.eqBands[0].type = eq::FilterType::LowShelf;
    pb.eqBands[0].lane (eq::Lane::Stereo).gainDb = -9.0;
    pb.eqBands[3].type = eq::FilterType::Notch;
    pb.eqBands[3].swept = true;
    pb.eqBands[5].type = eq::FilterType::HighPass;
    pb.eqBands[5].bypass = true;
    pb.eqBands[7].on = true;  pb.eqBands[7].type = eq::FilterType::Tilt;
    pb.eqBands[7].lane (eq::Lane::Stereo).gainDb = 3.0;
    pb.eqBands[9].on = true;  pb.eqBands[9].type = eq::FilterType::BandPass;
    pb.eqBands[9].lane (eq::Lane::Mid).on = true;          // a lane the first set never enabled
    pb.eqBands[9].lane (eq::Lane::Mid).freq = 600.0;
    pb.eqBands[9].dyn.on = true;                           // and the dynamic layer, whose seams are state

    for (const int topo : { 0, 1, 2 })
    {
        mastering::MasteringChainConfig cfg;
        if (topo == 1) { cfg.clipper = true; cfg.sidechainHpfHz = 80.0; }
        // TOPOLOGY 2 IS THE EQ ALONE, and it is the one that actually gates the claim: with a compressor, a
        // limiter and a 24-bit dither downstream, a small difference in the EQ's tail can be quantised away
        // before the float bits are compared, so a null over the full chain is weaker than it looks. Here
        // nothing stands between the engine and the comparison. (The diverse-testing round.)
        if (topo == 2) { cfg.compressor = false; cfg.limiter = false; cfg.dither = false; }

        // AND THE FIXTURE IS LIVE: these settings have to MOVE the programme, or the null below would hold
        // just as well over an EQ that did nothing at all.
        {
            Buf lit = x;
            mastering::MasteringChain probe;
            ok (probe.prepare (48000.0, nch, cfg), "PRECONDITION: a chain for topology " + std::to_string (topo));
            mastering::MasteringChainParams flat = pa;
            flat.bypassEq = true;
            probe.setParams (flat);
            { auto q = planes (lit); felitronics::test::run (probe.process (q.data(), nch, n)); }
            Buf eqd = x;
            mastering::MasteringChain probe2;
            ok (probe2.prepare (48000.0, nch, cfg), "PRECONDITION: and its twin");
            probe2.setParams (pa);
            { auto q = planes (eqd); felitronics::test::run (probe2.process (q.data(), nch, n)); }
            ok (! bitEqual (lit, eqd), "PRECONDITION: the EQ settings really move the programme, topology "
                                       + std::to_string (topo));
        }

        Buf once = x;
        mastering::MasteringChain a;
        ok (a.prepare (48000.0, nch, cfg), "PRECONDITION: a chain prepared once");
        a.setParams (pa);
        { auto p = planes (once); felitronics::test::run (a.process (p.data(), nch, n)); }

        Buf thrice = x;
        mastering::MasteringChain b;
        ok (b.prepare (48000.0, nch, cfg), "PRECONDITION: a chain prepared, ...");
        b.setParams (pb);                                   // a parameter epoch that must not survive
        { Buf junk = x; auto p = planes (junk); felitronics::test::run (b.process (p.data(), nch, 512)); }
        ok (b.prepare (48000.0, nch, cfg), "... re-prepared after audio, ...");
        b.setParams (pb);
        { Buf junk = x; auto p = planes (junk); felitronics::test::run (b.process (p.data(), nch, 777)); }
        ok (b.prepare (48000.0, nch, cfg), "... and re-prepared again");
        b.setParams (pa);                                   // the same set the first chain ran
        { auto p = planes (thrice); felitronics::test::run (b.process (p.data(), nch, n)); }

        ok (bitEqual (once, thrice),
            "a chain prepared three times renders bit-identically to one prepared once, topology "
            + std::to_string (topo) + " (first difference at sample " + std::to_string (firstDiff (once, thrice)) + ")");
    }

    // AND THE GEOMETRY MUST MOVE, or the null proves less than it looks. Every preparation above is at the
    // same rate, width and quantum, and a re-used engine that is NEVER RE-PREPARED still answers correctly
    // there — the chain's `applyParams()` and `reset()` restore everything that matters at ONE geometry.
    // Measured: a mutant that reached for the engine and skipped `eq_->prepare()` on re-use left the whole
    // suite green until this row. Here the second chain is re-prepared at a DIFFERENT rate and quantum, so a
    // stale engine is one designed for the wrong sample rate and the null says so.
    {
        constexpr double fsA = 48000.0, fsB = 96000.0;
        mastering::MasteringChainConfig ca;                  // the geometry it starts at
        mastering::MasteringChainConfig cb;                  // and the one it is moved to
        cb.internalBlock = 512;
        const Buf xb = burstThenSilence (nch, n, 900, 777u);

        Buf direct = xb;
        mastering::MasteringChain fresh;
        ok (fresh.prepare (fsB, nch, cb), "PRECONDITION: a chain built straight at the second geometry");
        fresh.setParams (pa);
        { auto q = planes (direct); felitronics::test::run (fresh.process (q.data(), nch, n)); }

        Buf moved = xb;
        mastering::MasteringChain mover;
        ok (mover.prepare (fsA, nch, ca), "PRECONDITION: and one that starts at the first, ...");
        mover.setParams (pb);
        { Buf junk = xb; auto q = planes (junk); felitronics::test::run (mover.process (q.data(), nch, 640)); }
        ok (mover.prepare (fsB, nch, cb), "... then moves to the second");
        mover.setParams (pa);
        { auto q = planes (moved); felitronics::test::run (mover.process (q.data(), nch, n)); }

        ok (bitEqual (direct, moved),
            "a chain MOVED to another rate and quantum renders bit-identically to one built there (first "
            "difference at sample " + std::to_string (firstDiff (direct, moved)) + ")");
    }

    // AND THE ENGINE REALLY IS THE SAME OBJECT — otherwise the null above would be proving nothing about
    // the reuse, only about the engine being deterministic. A re-preparation that rebuilt it would have to
    // ask the heap for `objectBytes()`; this asks for nothing at all.
    {
        mastering::MasteringChainConfig cfg;
        auto chain = std::make_unique<mastering::MasteringChain>();
        ok (chain->prepare (48000.0, 2, cfg), "PRECONDITION: prepared with an EQ");
        const long long before = g_bytes.load();
        const bool again = chain->prepare (48000.0, 2, cfg);
        const long long got = g_bytes.load() - before;
        ok (again && got == 0, "re-preparing asks the heap for nothing at all, engine included (read "
                               + std::to_string (got) + " B)");
    }
}

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
    testParamsWrittenBeforePrepare();
    testMonoBassParams();
    testGateCountsItsSubstitutions();
    testDemand();
    testEqEngineReuseIsBitIdentical();
    return felitronics::test::report();
}
