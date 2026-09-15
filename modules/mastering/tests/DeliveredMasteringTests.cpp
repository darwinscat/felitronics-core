// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P57b — mastering::DeliveredMastering: SRC first, for a render, a loudness search and a range. The oracle for each
// operation is the composition written out by hand from the parts it is made of — `DeliveryConverter` into a buffer
// of this file's own, then the class that does the job at one rate, on a chain prepared at the delivery rate. A
// composition that skipped the conversion, ran the chain at the source rate, rendered over the wrong buffer or
// dropped a frame differs from that by construction. The budgets are held against the bytes operator new saw.

#include <felitronics/mastering/DeliveredMastering.h>
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using namespace felitronics::mastering;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };
constexpr int    kNch = 2;
constexpr int    kBlock = 4096;

// Different content per channel, a slow envelope and transients, so a planar or a level fault shows.
std::vector<float> programme (double fs, long long frames)
{
    std::vector<float> v ((std::size_t) (frames * kNch), 0.0f);
    for (int c = 0; c < kNch; ++c)
        for (long long i = 0; i < frames; ++i)
        {
            const double t = (double) i / fs;
            const double env = 0.3 + 0.5 * (0.5 + 0.5 * std::sin (2.0 * kPi * 0.9 * t));
            double x = env * (0.55 * std::sin (2.0 * kPi * (330.0 + 110.0 * c) * t) + 0.25 * std::sin (2.0 * kPi * 5100.0 * t));
            if ((i % (long long) (fs / 4.0)) < 16) x += 0.5;
            v[(std::size_t) (c * frames + i)] = (float) x;
        }
    return v;
}

struct Planes
{
    const float* in[core::kMaxChannels] {};
    float* out[core::kMaxChannels] {};
};
Planes planes (const std::vector<float>& in, long long inFrames, std::vector<float>& out, long long outFrames)
{
    Planes p;
    for (int c = 0; c < kNch; ++c)
    {
        p.in[c]  = in.data()  + (std::size_t) (c * inFrames);
        p.out[c] = out.data() + (std::size_t) (c * outFrames);
    }
    return p;
}

long long bitDiff (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return (long long) (a.size() + b.size());
    long long n = 0;
    for (std::size_t i = 0; i < a.size(); ++i) if (std::memcmp (&a[i], &b[i], sizeof (float)) != 0) ++n;
    return n;
}

MasteringChainConfig lightConfig()
{
    MasteringChainConfig c;
    c.eq = false; c.monoBass = false; c.compressor = true; c.clipper = false; c.limiter = true; c.dither = true;
    return c;
}

// THE HAND-WRITTEN ORACLE: the programme converted into a buffer of our own.
bool convertByHand (double a, double b, const std::vector<float>& in, long long inFrames, std::vector<float>& conv, long long& d)
{
    d = DeliveryConverter::deliveredFrames (a, b, inFrames);
    if (d < 0) return false;
    conv.assign ((std::size_t) (d * kNch), 0.0f);
    DeliveryConverter dc;
    if (! dc.prepare (a, b, kNch, 1024)) return false;          // a different block than the class under test: free
    Planes p = planes (in, inFrames, conv, d);
    return dc.convert (p.in, kNch, inFrames, p.out, d);
}
}   // namespace

static void testRenderIsTheComposition()
{
    group ("render is conversion, then OfflineRenderer at the delivery rate — on every pair");
    int pairs = 0, differ = 0, lengthOff = 0, silent = 0;
    for (double a : kRates)
        for (double b : kRates)
        {
            const long long n = (long long) (a * 0.37);
            const auto in = programme (a, n);
            const MasteringChainConfig cfg = lightConfig();

            long long d = 0;
            std::vector<float> conv;
            if (! convertByHand (a, b, in, n, conv, d)) { ++differ; continue; }
            MasteringChain ch1; OfflineRenderer r1;
            std::vector<float> oracle (conv.size(), 0.0f);
            {
                if (! r1.prepare (kNch, 777) || ! ch1.prepare (b, kNch, cfg)) { ++differ; continue; }
                Planes p = planes (conv, d, oracle, d);
                if (! r1.render (ch1, p.in, p.out, kNch, (int) d)) { ++differ; continue; }
            }

            MasteringChain ch2; OfflineRenderer r2; DeliveredMastering dm;
            std::vector<float> got ((std::size_t) (d * kNch), 0.0f);
            if (! r2.prepare (kNch, kBlock) || ! ch2.prepare (b, kNch, cfg) || ! dm.prepare (a, b, kNch, kBlock)) { ++differ; continue; }
            Planes p = planes (in, n, got, d);
            if (! dm.render (ch2, r2, p.in, kNch, n, p.out, d)) { ++differ; continue; }
            ++pairs;
            if (DeliveredMastering::deliveredFrames (a, b, n) != d) ++lengthOff;
            if (bitDiff (oracle, got) != 0) ++differ;
            double peak = 0.0;
            for (float v : got) peak = std::max (peak, (double) std::fabs (v));
            if (peak < 0.05) ++silent;
        }
    ok (pairs == 36, "PRECONDITION: all 36 pairs of the six rates rendered (" + std::to_string (pairs) + ")");
    ok (differ == 0, "and every delivered render is the hand composition, bit for bit (" + std::to_string (differ) + " off)");
    ok (lengthOff == 0, "at the converter's length");
    ok (silent == 0, "PRECONDITION: none of them is silence");
}

static void testShortProgrammes()
{
    group ("one, two and 147 frames: the composition holds at lengths a long fixture never visits");
    int off = 0, runs = 0;
    for (const auto& pr : { std::pair<double, double> { 44100.0, 48000.0 }, std::pair<double, double> { 96000.0, 44100.0 },
                            std::pair<double, double> { 48000.0, 48000.0 } })
        for (long long n : { 1LL, 2LL, 147LL })
        {
            const auto in = programme (pr.first, n);
            const MasteringChainConfig cfg = lightConfig();
            long long d = 0;
            std::vector<float> conv;
            if (! convertByHand (pr.first, pr.second, in, n, conv, d)) { ++off; continue; }
            MasteringChain c1, c2; OfflineRenderer r1, r2; DeliveredMastering dm;
            std::vector<float> oracle (conv.size(), 0.0f), got (conv.size(), 0.5f);
            Planes po = planes (conv, d, oracle, d), pg = planes (in, n, got, d);
            const bool ran = r1.prepare (kNch, kBlock) && c1.prepare (pr.second, kNch, cfg) && r1.render (c1, po.in, po.out, kNch, (int) d)
                          && r2.prepare (kNch, kBlock) && c2.prepare (pr.second, kNch, cfg) && dm.prepare (pr.first, pr.second, kNch, kBlock)
                          && dm.render (c2, r2, pg.in, kNch, n, pg.out, d);
            ++runs;
            if (! ran || d < 1 || bitDiff (oracle, got) != 0) ++off;
        }
    ok (runs == 9 && off == 0, "nine short renders are the hand composition, bit for bit (" + std::to_string (off) + " off)");
}

static void testEqualRatesAreThePlainRender()
{
    group ("equal rates: the delivered render IS the plain render");
    const double fs = 48000.0;
    const long long n = 30000;
    const auto in = programme (fs, n);
    MasteringChainConfig cfg;                                  // the default topology, EQ engine and all
    MasteringChain c1, c2; OfflineRenderer r1, r2; DeliveredMastering dm;
    std::vector<float> plain (in.size(), 0.0f), delivered (in.size(), 0.0f);
    const bool built = r1.prepare (kNch, kBlock) && c1.prepare (fs, kNch, cfg)
                    && r2.prepare (kNch, kBlock) && c2.prepare (fs, kNch, cfg) && dm.prepare (fs, fs, kNch, kBlock);
    ok (built, "PRECONDITION: both built");
    Planes p1 = planes (in, n, plain, n), p2 = planes (in, n, delivered, n);
    ok (built && r1.render (c1, p1.in, p1.out, kNch, (int) n) && dm.render (c2, r2, p2.in, kNch, n, p2.out, n),
        "both rendered");
    ok (bitDiff (plain, delivered) == 0, "bit for bit");
}

static void testRenderRefusesAndAllocatesNothing()
{
    group ("render: refusals, block freedom, and no allocation");
    const double a = 44100.0, b = 96000.0;
    const long long n = 22050;
    const auto in = programme (a, n);
    const long long d = DeliveredMastering::deliveredFrames (a, b, n);
    MasteringChain chain, wrongRate; OfflineRenderer r; DeliveredMastering dm, unprepared;
    const MasteringChainConfig cfg = lightConfig();
    ok (r.prepare (kNch, kBlock) && chain.prepare (b, kNch, cfg) && wrongRate.prepare (a, kNch, cfg)
        && dm.prepare (a, b, kNch, kBlock), "PRECONDITION: built");
    std::vector<float> out ((std::size_t) (d * kNch), 0.0f);
    Planes p = planes (in, n, out, d);

    ok (! dm.render (wrongRate, r, p.in, kNch, n, p.out, d),
        "a chain at the SOURCE rate is refused — it would play the converted programme at the wrong speed");
    ok (! dm.render (chain, r, p.in, kNch, n, p.out, d - 1) && ! dm.render (chain, r, p.in, kNch, n, p.out, d + 1),
        "a delivered length one off either way is refused");
    ok (! dm.render (chain, r, p.in, 1, n, p.out, d), "a width that is not the converter's is refused");
    ok (! unprepared.render (chain, r, p.in, kNch, n, p.out, d), "an unprepared converter is refused");

    const long long allocs = alloc::count.load();
    const bool rendered = dm.render (chain, r, p.in, kNch, n, p.out, d);
    const long long asked = alloc::count.load() - allocs;
    ok (rendered && asked == 0, "the render asks the heap for nothing (" + std::to_string (asked) + ")");

    // THE CONVERTER'S BLOCK IS FREE — the class's own and the renderer's.
    int off = 0;
    for (int blk : { 1, 977, 65536 })
    {
        MasteringChain c2; OfflineRenderer r2; DeliveredMastering d2;
        std::vector<float> o2 (out.size(), 0.0f);
        Planes q = planes (in, n, o2, d);
        if (! r2.prepare (kNch, blk) || ! c2.prepare (b, kNch, cfg) || ! d2.prepare (a, b, kNch, blk)
            || ! d2.render (c2, r2, q.in, kNch, n, q.out, d) || bitDiff (out, o2) != 0) ++off;
    }
    ok (off == 0, "blocks 1, 977 and 65536 render the same bits");

    // AN EMPTY PROGRAMME IS LEGAL, AND THE PLANE RULE STILL APPLIES TO IT. The rule used to be skipped at zero
    // frames, and the render then read the output TABLE to hand it to the renderer: a null table was a crash (on
    // 37c95b4, SIGSEGV). In real tables it renders, as `OfflineRenderer` does at zero frames; without them, or with
    // null planes in them, it is refused.
    {
        std::vector<float> none (1, 0.0f);
        Planes real = planes (none, 0, none, 0);
        ok (dm.render (chain, r, real.in, kNch, 0, real.out, 0), "an empty programme in real tables renders");
        const float* nullIn[core::kMaxChannels] {};
        float* nullOut[core::kMaxChannels] {};
        ok (! dm.render (chain, r, nullIn, kNch, 0, nullOut, 0), "an empty programme with null planes is refused");
        ok (! dm.render (chain, r, nullptr, kNch, 0, nullptr, 0), "and one with null tables is refused, not a crash");
    }
}

static void testSolveAndRangeAreTheComposition()
{
    group ("solve and the range: conversion, then the solver at the delivery rate — and the budgets are the bytes");
    const LoudnessRequest req = [] { LoudnessRequest q; q.targetLufs = -14.0; q.maxTruePeakDbTp = -1.0; return q; } ();
    struct Pair { double a, b; };
    for (const Pair pr : { Pair { 44100.0, 48000.0 }, Pair { 96000.0, 44100.0 }, Pair { 48000.0, 48000.0 } })
    {
        const std::string tag = std::to_string ((int) pr.a) + " -> " + std::to_string ((int) pr.b) + ": ";
        const long long n = (long long) (pr.a * 3.5);
        const auto in = programme (pr.a, n);
        const MasteringChainConfig cfg = lightConfig();
        MasteringChainParams params;

        // The oracle.
        long long d = 0;
        std::vector<float> conv;
        const bool converted = convertByHand (pr.a, pr.b, in, n, conv, d);
        MasteringChain c1; OfflineRenderer r1; TargetLoudnessSolver s1;
        const bool b1 = converted && r1.prepare (kNch, kBlock) && c1.prepare (pr.b, kNch, cfg)
                     && s1.prepare (pr.b, kNch, kBlock, c1.internalBlock(), c1.tapOversampleFactor());
        std::vector<float> oracle (conv.size(), 0.0f);
        Planes po = planes (conv, d, oracle, d);
        double lraOracle = -1.0;
        const bool lra1 = b1 && s1.measureInputLoudnessRange (po.in, kNch, (int) d, lraOracle);
        const LoudnessSolution so = b1 ? s1.solve (c1, r1, params, po.in, po.out, kNch, (int) d, req) : LoudnessSolution {};

        // The class.
        MasteringChain c2; OfflineRenderer r2; TargetLoudnessSolver s2; DeliveredMastering dm;
        const bool b2 = r2.prepare (kNch, kBlock) && c2.prepare (pr.b, kNch, cfg) && dm.prepare (pr.a, pr.b, kNch, kBlock)
                     && s2.prepare (pr.b, kNch, kBlock, c2.internalBlock(), c2.tapOversampleFactor());
        ok (b1 && b2, tag + "PRECONDITION: both built");
        std::vector<float> got ((std::size_t) (d * kNch), 0.0f);
        Planes pg = planes (in, n, got, d);

        const long long identity = (pr.a == pr.b) ? 1 : 0;
        const std::uint64_t lraBudget = DeliveredMastering::measureRangeBytes (pr.a, pr.b, kNch, n);
        double lra = -2.0;
        long long before = alloc::bytes.load();
        const bool lra2 = dm.measureInputLoudnessRange (s2, pg.in, kNch, n, lra);
        const long long lraBytes = alloc::bytes.load() - before;
        ok (lra1 && lra2 && lra == lraOracle, tag + "the range is the hand composition's, exactly");
        ok (lraBytes == (long long) lraBudget, tag + "and allocates its budget (" + std::to_string (lraBytes) + " against "
            + std::to_string (lraBudget) + ")");
        ok (lraBudget == TargetLoudnessSolver::measureRangeBytes (pr.b, (int) d)
                         + (identity ? 0u : (std::uint64_t) (kNch * d * 4)),
            tag + "which is the meter at the delivered length plus the converted programme — none at equal rates");

        const std::uint64_t solveBudget = DeliveredMastering::solveBytes (pr.a, pr.b, kNch, n);
        before = alloc::bytes.load();
        const LoudnessSolution sg = dm.solve (s2, c2, r2, params, pg.in, kNch, n, pg.out, d, req);
        const long long solveBytes = alloc::bytes.load() - before;
        ok (sg.status == so.status && sg.preLimiterGainDb == so.preLimiterGainDb && sg.ceilingDbTp == so.ceilingDbTp
            && sg.passes == so.passes && sg.passes > 0, tag + "the search's verdict is the hand composition's");
        ok (bitDiff (oracle, got) == 0, tag + "and so is its delivered audio, bit for bit");
        const long long programmeBytes = identity ? 0 : (long long) kNch * d * 4;
        const long long perPass = (long long) solveBudget - programmeBytes;
        ok (perPass > 0 && solveBytes == (long long) sg.passes * perPass + programmeBytes,
            tag + "the search allocates passes x its meters + the converted programme ("
            + std::to_string (solveBytes) + ")");
    }

    // THE RANGE IS JUDGED ON THE DELIVERED LENGTH, with nothing converted on the way to a refusal.
    {
        const double a = 44100.0, b = 96000.0;
        MasteringChain c; OfflineRenderer r; TargetLoudnessSolver s; DeliveredMastering dm;
        ok (r.prepare (kNch, kBlock) && c.prepare (b, kNch, lightConfig()) && dm.prepare (a, b, kNch, kBlock)
            && s.prepare (b, kNch, kBlock, c.internalBlock(), c.tapOversampleFactor()), "PRECONDITION: 44.1 -> 96 built");
        const auto in = programme (a, 132300);
        std::vector<float> dummy (1, 0.0f);
        Planes p = planes (in, 132300, dummy, 0);
        ok (DeliveredMastering::deliveredFrames (a, b, 132300) == 288000 && DeliveredMastering::deliveredFrames (a, b, 132299) < 288000,
            "PRECONDITION: 132 300 input frames deliver exactly 3 s at 96 kHz, one fewer does not");
        double v = 0.0;
        ok (DeliveredMastering::measureRangeBytes (a, b, kNch, 132300) > 0u && dm.measureInputLoudnessRange (s, p.in, kNch, 132300, v),
            "3 s delivered: budgeted and measured");
        const long long before = alloc::bytes.load();
        const bool refused = ! dm.measureInputLoudnessRange (s, p.in, kNch, 132299, v);
        const long long asked = alloc::bytes.load() - before;
        ok (DeliveredMastering::measureRangeBytes (a, b, kNch, 132299) == 0u && refused && asked == 0,
            "one frame short: budgeted 0, refused, nothing allocated");
    }
}

static void testSolveRefusals()
{
    group ("solve: refusals are the solver's own verdicts, and allocate nothing");
    const double a = 48000.0, b = 44100.0;
    const long long n = 48000;
    const auto in = programme (a, n);
    const long long d = DeliveredMastering::deliveredFrames (a, b, n);
    LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
    MasteringChain chain, wrongRate; OfflineRenderer r; TargetLoudnessSolver s; DeliveredMastering dm, unprepared;
    ok (r.prepare (kNch, kBlock) && chain.prepare (b, kNch, lightConfig()) && wrongRate.prepare (a, kNch, lightConfig())
        && dm.prepare (a, b, kNch, kBlock) && s.prepare (b, kNch, kBlock, chain.internalBlock(), chain.tapOversampleFactor()),
        "PRECONDITION: built");
    std::vector<float> out ((std::size_t) (d * kNch), 0.0f);
    Planes p = planes (in, n, out, d);
    const MasteringChainParams params;
    long long before = alloc::bytes.load();
    const auto st1 = unprepared.solve (s, chain, r, params, p.in, kNch, n, p.out, d, req).status;
    const auto st2 = dm.solve (s, wrongRate, r, params, p.in, kNch, n, p.out, d, req).status;
    const auto st3 = dm.solve (s, chain, r, params, p.in, kNch, n, p.out, d + 1, req).status;
    const long long asked = alloc::bytes.load() - before;
    ok (st1 == MasteringSolveStatus::NotPrepared, "an unprepared converter: NotPrepared");
    ok (st2 == MasteringSolveStatus::InvalidRequest, "a chain at the source rate: InvalidRequest");
    ok (st3 == MasteringSolveStatus::InvalidRequest, "a delivered length that is not the converter's: InvalidRequest");
    ok (asked == 0, "none of them converted anything (" + std::to_string (asked) + " B)");
    ok (dm.solve (s, chain, r, params, p.in, kNch, 0, p.out, 0, req).status == MasteringSolveStatus::InvalidRequest,
        "an empty programme is forwarded, and the solver's own verdict comes back");

    // `TargetLoudnessSolver::admits` IS NOW THE ONE DEFINITION both paths ask, so each clause of it is pinned here —
    // the mutation stand found the true-peak aim's clause unpinned anywhere in the tree.
    MasteringSolveStatus why = MasteringSolveStatus::Solved;
    for (double aim : { std::numeric_limits<double>::quiet_NaN(), -0.01, std::numeric_limits<double>::infinity() })
    {
        LoudnessRequest bad = req; bad.truePeakAimDb = aim;
        const long long b0 = alloc::bytes.load();
        const auto st = dm.solve (s, chain, r, params, p.in, kNch, n, p.out, d, bad).status;
        const long long asked = alloc::bytes.load() - b0;
        ok (! s.admits (chain, r, kNch, (int) d, bad, why) && why == MasteringSolveStatus::InvalidRequest
            && st == MasteringSolveStatus::InvalidRequest && asked == 0,
            "a true-peak aim of " + std::to_string (aim) + ": refused by admits and by the delivered solve, nothing allocated");
    }
    ok (s.admits (chain, r, kNch, (int) d, req, why), "and the good request is admitted");
}

static void testRefusalsWriteAndAllocateNothing()
{
    group ("every refusal comes before the converter writes a sample and before a programme is allocated");
    const double a = 48000.0, b = 96000.0;
    const long long n = 144000;                                 // 3 s in, 6 s delivered
    const auto in = programme (a, n);
    const long long d = DeliveredMastering::deliveredFrames (a, b, n);
    const MasteringChainConfig cfg = lightConfig();
    MasteringChain chain; OfflineRenderer r, unprepared; DeliveredMastering dm;
    MasteringChain mono; MasteringChainConfig monoCfg = cfg;
    TargetLoudnessSolver s, sUnprepared, sWrongRate;
    ok (r.prepare (kNch, kBlock) && chain.prepare (b, kNch, cfg) && mono.prepare (b, 1, monoCfg) && dm.prepare (a, b, kNch, kBlock)
        && s.prepare (b, kNch, kBlock, chain.internalBlock(), chain.tapOversampleFactor())
        && sWrongRate.prepare (a, kNch, kBlock, chain.internalBlock(), chain.tapOversampleFactor()), "PRECONDITION: built");
    LoudnessRequest req; req.targetLufs = -14.0; req.maxTruePeakDbTp = -1.0;
    const MasteringChainParams params;

    const std::vector<float> sentinel ((std::size_t) (d * kNch), 0.125f);
    auto untouchedAndFree = [&] (const char* what, auto&& call)
    {
        std::vector<float> out = sentinel;
        Planes p = planes (in, n, out, d);
        const long long before = alloc::bytes.load();
        const bool accepted = call (p);
        const long long asked = alloc::bytes.load() - before;
        ok (! accepted && out == sentinel && asked == 0,
            std::string (what) + ": refused, the output untouched, nothing allocated (" + std::to_string (asked) + " B)");
    };
    untouchedAndFree ("render with an unprepared renderer",
                      [&] (Planes& p) { return dm.render (chain, unprepared, p.in, kNch, n, p.out, d); });
    untouchedAndFree ("render on a chain one channel narrower than the converter",
                      [&] (Planes& p) { return dm.render (mono, r, p.in, kNch, n, p.out, d); });
    untouchedAndFree ("solve with an unprepared solver",
                      [&] (Planes& p) { return dm.solve (sUnprepared, chain, r, params, p.in, kNch, n, p.out, d, req).status
                                               != MasteringSolveStatus::NotPrepared; });
    untouchedAndFree ("solve with a request the solver refuses (no target)",
                      [&] (Planes& p) { LoudnessRequest bad; return dm.solve (s, chain, r, params, p.in, kNch, n, p.out, d, bad).status
                                                                   != MasteringSolveStatus::InvalidRequest; });
    untouchedAndFree ("solve with a solver at the SOURCE rate",
                      [&] (Planes& p) { return dm.solve (sWrongRate, chain, r, params, p.in, kNch, n, p.out, d, req).status
                                               != MasteringSolveStatus::InvalidRequest; });
    {
        double v = -1.0;
        const long long before = alloc::bytes.load();
        Planes p = planes (in, n, const_cast<std::vector<float>&> (sentinel), d);
        const bool measured = dm.measureInputLoudnessRange (sWrongRate, p.in, kNch, n, v);
        const long long asked = alloc::bytes.load() - before;
        ok (! measured && v == -1.0 && asked == 0,
            "a range measured by a solver at the SOURCE rate is refused — it would meter 96 kHz samples on a 48 kHz grid");
    }

    // OVERLAPPING PLANES: the conversion writes `out` while it still reads `in`.
    {
        std::vector<float> shared ((std::size_t) (d * kNch + n * kNch), 0.0f);
        std::copy (in.begin(), in.end(), shared.begin());
        const float* ip[core::kMaxChannels] {}; float* op[core::kMaxChannels] {};
        for (int c = 0; c < kNch; ++c) { ip[c] = shared.data() + (std::size_t) (c * n); op[c] = shared.data() + 7 + (std::size_t) (c * d); }
        const std::vector<float> before = shared;
        ok (! dm.render (chain, r, ip, kNch, n, op, d) && shared == before, "render over its own input: refused, nothing written");
        ok (dm.solve (s, chain, r, params, ip, kNch, n, op, d, req).status == MasteringSolveStatus::InvalidRequest && shared == before,
            "solve over its own input: InvalidRequest, nothing written");
        // ACROSS CHANNELS ONLY: output plane 1 starts inside input plane 0, and no plane touches its own namesake.
        std::vector<float> in0 (in.begin(), in.begin() + n), in1 (in.begin() + n, in.end());
        std::vector<float> out0 ((std::size_t) d, 0.0f), cross ((std::size_t) (n + d), 0.0f);
        std::copy (in0.begin(), in0.end(), cross.begin());
        const float* ipx[core::kMaxChannels] { cross.data(), in1.data() };
        float* opx[core::kMaxChannels] { out0.data(), cross.data() + 5 };
        const std::vector<float> crossBefore = cross;
        ok (! dm.render (chain, r, ipx, kNch, n, opx, d) && cross == crossBefore,
            "an output plane over ANOTHER channel's input: refused, nothing written");

        // UPSAMPLING: these spans overlap only beyond the first `inFrames` of the longer output plane.
        std::vector<float> longStride ((std::size_t) (10 * n), 0.125f);
        const float* ipl[core::kMaxChannels] { longStride.data() + n + 1, longStride.data() + 4 * n + 4 };
        float* opl[core::kMaxChannels] { longStride.data(), longStride.data() + 7 * n + 8 };
        std::copy (in.begin(), in.begin() + n, const_cast<float*> (ipl[0]));
        std::copy (in.begin() + n, in.end(), const_cast<float*> (ipl[1]));
        const std::vector<float> longBefore = longStride;
        ok (! dm.render (chain, r, ipl, kNch, n, opl, d) && longStride == longBefore,
            "render: overlap only in the longer output stride is refused, nothing written");
        const auto longSol = dm.solve (s, chain, r, params, ipl, kNch, n, opl, d, req);
        ok (longSol.status == MasteringSolveStatus::InvalidRequest && longSol.passes == 0 && longStride == longBefore,
            "solve: the same longer-stride overlap is InvalidRequest before a render");
    }

    // EQUAL RATES read the input in place, so a null input must be refused before it is read.
    {
        MasteringChain c48; OfflineRenderer r48; TargetLoudnessSolver s48; DeliveredMastering same;
        ok (r48.prepare (kNch, kBlock) && c48.prepare (a, kNch, cfg) && same.prepare (a, a, kNch, kBlock)
            && s48.prepare (a, kNch, kBlock, c48.internalBlock(), c48.tapOversampleFactor()), "PRECONDITION: 48 -> 48 built");
        std::vector<float> out ((std::size_t) (n * kNch), 0.0f);
        float* op[core::kMaxChannels] {};
        for (int c = 0; c < kNch; ++c) op[c] = out.data() + (std::size_t) (c * n);
        ok (same.solve (s48, c48, r48, params, nullptr, kNch, n, op, n, req).status == MasteringSolveStatus::InvalidRequest,
            "equal rates, a null input table: InvalidRequest, not a dereference");
        const float* nullPlanes[core::kMaxChannels] {};
        ok (same.solve (s48, c48, r48, params, nullPlanes, kNch, n, op, n, req).status == MasteringSolveStatus::InvalidRequest,
            "and a table of null planes");
        double v = 0.0;
        Planes p = planes (in, n, out, n);
        ok (! same.measureInputLoudnessRange (s48, p.in, 1, n, v), "and a range asked at a width that is not the converter's");
    }

    // THE RANGE ON A DOWNSAMPLE, judged on the delivered length: 287 999 frames at 96 kHz are 2.99999 s of input and
    // exactly 132 300 = 3 s at 44.1 kHz — measurable here, where a plain 96 kHz handle refuses them. That is the
    // contract (the meter sees the delivered programme), pinned so a change to it is a decision rather than a drift.
    {
        MasteringChain c441; OfflineRenderer r441; TargetLoudnessSolver s441; DeliveredMastering down;
        ok (r441.prepare (kNch, kBlock) && c441.prepare (44100.0, kNch, cfg) && down.prepare (96000.0, 44100.0, kNch, kBlock)
            && s441.prepare (44100.0, kNch, kBlock, c441.internalBlock(), c441.tapOversampleFactor()), "PRECONDITION: 96 -> 44.1 built");
        ok (DeliveredMastering::deliveredFrames (96000.0, 44100.0, 287999) == 132300
            && TargetLoudnessSolver::measureRangeBytes (96000.0, 287999) == 0u,
            "PRECONDITION: 287 999 frames deliver 3 s, and read at the source rate are under 3 s");
        const auto in96 = programme (96000.0, 287999);
        std::vector<float> dummy (1, 0.0f);
        Planes p = planes (in96, 287999, dummy, 0);
        double v = 0.0;
        ok (DeliveredMastering::measureRangeBytes (96000.0, 44100.0, kNch, 287999) > 0u
            && down.measureInputLoudnessRange (s441, p.in, kNch, 287999, v), "so the delivered range is budgeted and measured");
    }
}

static void testNonFiniteInputIsGatedBeforeTheConversion()
{
    group ("a non-finite input sample is converted as the chain's gate would sanitise it, and counted once");
    for (const auto& pr : { std::pair<double, double> { 44100.0, 48000.0 }, std::pair<double, double> { 96000.0, 44100.0 },
                            std::pair<double, double> { 48000.0, 48000.0 } })
    {
        const std::string tag = std::to_string ((int) pr.first) + " -> " + std::to_string ((int) pr.second) + ": ";
        const long long n = (long long) pr.first;
        auto clean = programme (pr.first, n);
        auto dirty = clean;
        dirty[500] = std::numeric_limits<float>::quiet_NaN();                  clean[500] = 0.0f;
        dirty[(std::size_t) n + 9000] = -std::numeric_limits<float>::infinity(); clean[(std::size_t) n + 9000] = 0.0f;
        dirty[20000] = 3.0e38f;                                                clean[20000] = 1.0e6f;
        const long long d = DeliveredMastering::deliveredFrames (pr.first, pr.second, n);
        const MasteringChainConfig cfg = lightConfig();
        auto render = [&] (const std::vector<float>& src, std::vector<float>& out, std::uint64_t& bad)
        {
            MasteringChain ch; OfflineRenderer r; DeliveredMastering dm;
            out.assign ((std::size_t) (d * kNch), 0.0f);
            Planes p = planes (src, n, out, d);
            const bool okay = r.prepare (kNch, kBlock) && ch.prepare (pr.second, kNch, cfg) && dm.prepare (pr.first, pr.second, kNch, kBlock)
                           && dm.render (ch, r, p.in, kNch, n, p.out, d);
            bad = dm.nonFiniteInputSamples();
            return okay;
        };
        std::vector<float> oc, od; std::uint64_t bc = 9, bd = 9;
        ok (render (clean, oc, bc) && render (dirty, od, bd), tag + "both rendered");
        std::size_t nonFinite = 0;
        for (float v : od) if (! std::isfinite (v)) ++nonFinite;
        ok (bitDiff (oc, od) == 0 && nonFinite == 0, tag + "NaN, -inf and 3e38 render exactly as 0, 0 and 1e6 (" + std::to_string (bitDiff (oc, od))
            + " samples differ)");
        ok (bc == 0u && bd == 2u, tag + "and the two non-finite samples are counted, the clamped finite one is not");

        // The SAME count from the range measurement, which at equal rates reads the input in place rather than
        // converting it — and there the measurement itself refuses the poisoned programme, but the count is taken
        // first. (Converting, the gate has already replaced both samples, so the range is measured.)
        {
            MasteringChain ch; OfflineRenderer r; TargetLoudnessSolver sv; DeliveredMastering dm;
            const bool built = r.prepare (kNch, kBlock) && ch.prepare (pr.second, kNch, cfg) && dm.prepare (pr.first, pr.second, kNch, kBlock)
                            && sv.prepare (pr.second, kNch, kBlock, ch.internalBlock(), ch.tapOversampleFactor());
            const long long n4 = 4 * n;                           // a range needs 3 s delivered; 1 s would be refused first
            auto dirty4 = programme (pr.first, n4);
            dirty4[500] = std::numeric_limits<float>::quiet_NaN();
            dirty4[(std::size_t) n4 + 9000] = -std::numeric_limits<float>::infinity();
            std::vector<float> dummy (1, 0.0f);
            Planes p = planes (dirty4, n4, dummy, 0);
            double v = 0.0;
            const bool measured = dm.measureInputLoudnessRange (sv, p.in, kNch, n4, v);
            const bool identity = pr.first == pr.second;
            ok (built && measured != identity && dm.nonFiniteInputSamples() == 2u,
                tag + (identity ? "the range measurement counts the same two — and refuses the programme AFTER its count, keeping it"
                                : "the range measurement counts the same two, and measures the gated programme"));
            const float* nullPlanes[core::kMaxChannels] {};
            ok (! dm.measureInputLoudnessRange (sv, nullPlanes, kNch, n4, v) && dm.nonFiniteInputSamples() == 2u,
                tag + "and refuses a table of null planes BEFORE its count, leaving the previous one");
        }

        // THE COUNT IS THE LAST CALL'S THAT REACHED IT. One object through a sequence: every call refused before its
        // count leaves the previous number exactly where it was — the rule `fc_solution_log` keeps for `written` —
        // and the next call that reaches its count replaces it, zero included.
        {
            MasteringChain ch; OfflineRenderer r; TargetLoudnessSolver sv; DeliveredMastering dm;
            const bool built = r.prepare (kNch, kBlock) && ch.prepare (pr.second, kNch, cfg) && dm.prepare (pr.first, pr.second, kNch, kBlock)
                            && sv.prepare (pr.second, kNch, kBlock, ch.internalBlock(), ch.tapOversampleFactor());
            std::vector<float> out ((std::size_t) (d * kNch), 0.0f);
            Planes pd = planes (dirty, n, out, d);
            ok (built && dm.render (ch, r, pd.in, kNch, n, pd.out, d) && dm.nonFiniteInputSamples() == 2u,
                tag + "a render that reaches its count: 2");
            ok (! dm.render (ch, r, pd.in, kNch, n, pd.out, d - 1) && dm.nonFiniteInputSamples() == 2u,
                tag + "a render refused on its length leaves 2");
            Planes over = planes (dirty, n, dirty, d);                // the output over the input
            ok (! dm.render (ch, r, over.in, kNch, n, over.out, d) && dm.nonFiniteInputSamples() == 2u,
                tag + "a render refused on its planes leaves 2");
            Planes pc = planes (clean, n, out, d);
            LoudnessRequest noTarget;
            ok (dm.solve (sv, ch, r, MasteringChainParams {}, pc.in, kNch, n, pc.out, d, noTarget).status
                    == MasteringSolveStatus::InvalidRequest && dm.nonFiniteInputSamples() == 2u,
                tag + "a solve refused before its count leaves 2");
            ok (dm.render (ch, r, pc.in, kNch, n, pc.out, d) && dm.nonFiniteInputSamples() == 0u,
                tag + "and the next render that reaches its count replaces it: 0");
        }
    }
}

static void testCreateBudget()
{
    group ("createBytes: the chain at the delivery rate plus the converter, and 0 exactly where either refuses");
    const MasteringChainConfig cfg;
    ok (DeliveredMastering::createBytes (48000.0, 96000.0, kNch, cfg, kBlock)
            == createBytes (96000.0, kNch, cfg, kBlock) + DeliveryConverter::prepareBytes (48000.0, 96000.0, kNch, kBlock),
        "the sum of the two core budgets, at the DELIVERY rate for the chain");
    ok (DeliveredMastering::createBytes (48000.0, 44100.5, kNch, cfg, kBlock) == 0u, "a rate the resampler cannot plan: 0");
    ok (DeliveredMastering::createBytes (48000.0, 96000.0, 0, cfg, kBlock) == 0u, "a width the chain refuses: 0");

    const std::uint64_t budget = DeliveredMastering::createBytes (48000.0, 96000.0, kNch, MasteringChainConfig {}, kBlock);
    const std::uint64_t chainOnly = createBytes (96000.0, kNch, MasteringChainConfig {}, kBlock);
    const long long before = alloc::bytes.load();
    {
        DeliveredMastering dm;
        const bool okPrep = dm.prepare (48000.0, 96000.0, kNch, kBlock);
        const long long got = alloc::bytes.load() - before;
        ok (okPrep && got == (long long) (budget - chainOnly), "preparing the converter allocates its part of the budget, to the byte ("
            + std::to_string (got) + ")");
    }
}

int main()
{
    std::printf ("felitronics::mastering::DeliveredMastering — SRC first, render / solve / range\n");
    testRenderIsTheComposition();
    testShortProgrammes();
    testEqualRatesAreThePlainRender();
    testRenderRefusesAndAllocatesNothing();
    testSolveAndRangeAreTheComposition();
    testSolveRefusals();
    testRefusalsWriteAndAllocateNothing();
    testNonFiniteInputIsGatedBeforeTheConversion();
    testCreateBudget();
    return felitronics::test::report();
}
