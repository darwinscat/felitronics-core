// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The mastering chain's EQ stage drives `dynamiceq::LaneDynamics`, so `eqBands[].dyn` is live. What is
// checked here, and what each check is measured AGAINST:
//
//   * UNARMED IS THE ENGINE. With every point's dynamics off the chain's EQ stage is
//     `eq::EqEngine::process` over the same quantum, sample for sample — the null that keeps a default
//     render what it was.
//   * ARMED IS THE PRODUCER. An armed point is `eq::EqBand` + `dynamiceq::LaneDynamics` built OUTSIDE
//     this chain and driven on the same signal, sample for sample. The oracle is the producer itself,
//     not a restatement of its formula.
//   * The cut appears where the programme crosses the threshold and is EXACTLY ZERO where it does not.
//   * The KEY is the EQ section's input, so a point ahead of the armed one changes what it filters and
//     not what it detects on.
//   * The five lanes detect on their own domain only.
//   * Block invariance, zero latency, no allocation in process(), and the published memory budget.
//   * A bypass stretch stops the detectors, and the loudness search sees the same signal at every drive.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqEngine.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using felitronics::test::okNoAlloc;

namespace
{

constexpr double kFs  = 48000.0;
constexpr int    kNch = 2;
constexpr int    kK   = 256;                 // the chain's internal quantum
constexpr double kF0  = 3000.0;              // the band, and the tone that excites it

// EQ ONLY, so what the null compares is the EQ stage and not five other stages' arithmetic.
mastering::MasteringChainConfig eqOnlyConfig()
{
    mastering::MasteringChainConfig c;
    c.internalBlock = kK;
    c.eq = true;
    c.monoBass = c.compressor = c.clipper = c.limiter = c.dither = false;
    return c;
}

// One armed Bell point on the Stereo lane, absolute threshold. ABSOLUTE rather than relative because a
// threshold the programme is provably under gives a delta of EXACTLY zero, which is what the
// "no cut in a quiet place" checks below are written in.
eq::BandParams bellPoint (double freq, double Q, double rangeDb, double thrDb, bool on, eq::Lane lane)
{
    eq::BandParams p;
    p.on = true;
    p.type = eq::FilterType::Bell;
    for (auto& l : p.lanes) l.on = false;
    eq::LaneParams& lp = p.lane (lane);
    lp.on = true; lp.freq = freq; lp.Q = Q; lp.gainDb = 0.0;
    p.dyn.on = on;
    p.dyn.rangeDb = rangeDb;
    p.dyn.thrAuto = false;
    p.dyn.thrDb = thrDb;
    return p;
}

//==============================================================================
// PROGRAMME — a tone at the band centre that spends part of its time far above the threshold and the
// rest far below it, so one fixture carries both the "there is a cut" and the "there is none" case.
// `left`/`right`/`invertRight` choose which domain the loud part lands in, which is what the lane checks
// steer with.

struct Plan
{
    double quietDb = -34.0, loudDb = -6.0;
    bool   left = true, right = true;        // which channels carry the loud part
    double invertRight = 1.0;                // -1 puts the whole tone in Side
    // A steady component well below the band, so a fixture can carry loudness the point does not eat.
    // Off by default — the checks above want the band to be the whole programme.
    double bassDb = -1000.0;
};

std::vector<float> programme (const Plan& plan, int frames)
{
    std::vector<float> planar ((std::size_t) frames * 2, 0.0f);
    const double qa = core::dbToGain (plan.quietDb), la = core::dbToGain (plan.loudDb);
    const double ba = plan.bassDb > -999.0 ? core::dbToGain (plan.bassDb) : 0.0;
    for (int i = 0; i < frames; ++i)
    {
        const double t = (double) i / kFs;
        // 150 ms of loud per 500 ms — long enough for the ballistics to settle either side of the edge.
        const bool loud = std::fmod (t, 0.5) >= 0.35;
        const double s = std::sin (2.0 * core::kPi * kF0 * t);
        const double b = ba * std::sin (2.0 * core::kPi * 220.0 * t);
        const double aL = (loud && plan.left)  ? la : qa;
        const double aR = (loud && plan.right) ? la : qa;
        planar[(std::size_t) i]                          = (float) (aL * s + b);
        planar[(std::size_t) frames + (std::size_t) i]   = (float) (plan.invertRight * (aR * s) + b);
    }
    return planar;
}

// `frames` of programme through the chain, plus its latency of flush, with the priming dropped — the
// same formula `OfflineRenderer` uses, written here so this file needs no renderer.
std::vector<float> render (mastering::MasteringChain& chain, const std::vector<float>& in,
                           int frames, int block)
{
    const int D = chain.latencySamples();
    const std::size_t total = (std::size_t) frames + (std::size_t) D;
    std::vector<float> buf (total * (std::size_t) kNch, 0.0f);
    for (int c = 0; c < kNch; ++c)
        std::copy_n (in.data() + (std::size_t) c * (std::size_t) frames, frames,
                     buf.data() + (std::size_t) c * total);

    float* p[core::kMaxChannels] {};
    for (std::size_t off = 0; off < total; )
    {
        const int n = (int) std::min ((std::size_t) block, total - off);
        for (int c = 0; c < kNch; ++c) p[c] = buf.data() + (std::size_t) c * total + off;
        if (! chain.process (p, kNch, n)) return {};
        off += (std::size_t) n;
    }

    std::vector<float> out ((std::size_t) frames * (std::size_t) kNch, 0.0f);
    for (int c = 0; c < kNch; ++c)
        std::copy_n (buf.data() + (std::size_t) c * total + (std::size_t) D, frames,
                     out.data() + (std::size_t) c * (std::size_t) frames);
    return out;
}

std::vector<float> renderWith (const eq::BandParams& band, const std::vector<float>& in, int frames,
                               int block = 1024)
{
    auto chain = std::make_unique<mastering::MasteringChain>();
    if (! chain->prepare (kFs, kNch, eqOnlyConfig())) return {};
    mastering::MasteringChainParams mp;
    mp.eqBands[0] = band;
    chain->setParams (mp);
    return render (*chain, in, frames, block);
}

bool sameBits (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size() || a.empty()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (! core::sameBits (a[i], b[i])) return false;
    return true;
}

// The first `whole` frames of every plane of a planar render, concatenated — EMPTY if the render was
// refused and came back short, so a comparison against a refused render fails rather than indexing past
// the end of it.
std::vector<float> head (const std::vector<float>& v, int frames, int whole)
{
    if (frames < 0 || whole < 0 || whole > frames) return {};
    const std::size_t n = (std::size_t) frames, w = (std::size_t) whole;
    if (v.size() < n * (std::size_t) kNch) return {};
    std::vector<float> out;
    for (std::size_t c = 0; c < (std::size_t) kNch; ++c)
        out.insert (out.end(), v.begin() + (std::ptrdiff_t) (c * n),
                    v.begin() + (std::ptrdiff_t) (c * n + w));
    return out;
}

// How far apart two renders are, in words — "sizes differ" for a REFUSED render, which is a different
// thing from a render that came out wrong and must not be printed as a sample count.
std::string differingSamples (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return "sizes differ: " + std::to_string (a.size())
                                                             + " against " + std::to_string (b.size());
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); ++i) if (! core::sameBits (a[i], b[i])) ++n;
    return std::to_string (n) + " of " + std::to_string (a.size());
}

// The gate the chain applies to every input sample, so a reference built outside it starts from the
// same samples. Bit-transparent for the fixtures here; spelled out rather than assumed.
void gate (std::vector<float>& v)
{
    for (float& x : v) x = std::clamp (std::isfinite (x) ? x : 0.0f, -1.0e6f, 1.0e6f);
}

//==============================================================================
// 1. UNARMED IS THE ENGINE — `eq::EqEngine::process` over the same quantum, sample for sample.

void testUnarmedIsTheEngine()
{
    group ("with no point armed the EQ stage is eq::EqEngine::process, sample for sample");

    constexpr int frames = 24000;
    const std::vector<float> in = programme (Plan {}, frames);

    // Four points, three of them irrelevant to the dynamics, so the null covers the whole band loop
    // rather than band 0 alone.
    mastering::MasteringChainParams mp;
    mp.eqBands[0] = bellPoint (kF0, 2.0, -24.0, -30.0, false, eq::Lane::Stereo);
    mp.eqBands[3] = bellPoint (200.0, 0.8, 0.0, -24.0, false, eq::Lane::Stereo);
    mp.eqBands[3].lane (eq::Lane::Stereo).gainDb = -4.0;
    mp.eqBands[7] = bellPoint (9000.0, 1.2, 0.0, -24.0, false, eq::Lane::Left);
    mp.eqBands[7].lane (eq::Lane::Left).gainDb = 3.0;
    mp.eqBands[11] = bellPoint (700.0, 3.0, 0.0, -24.0, false, eq::Lane::Side);
    mp.eqBands[11].lane (eq::Lane::Side).gainDb = 2.0;

    auto chain = std::make_unique<mastering::MasteringChain>();
    test::run (chain->prepare (kFs, kNch, eqOnlyConfig()));
    chain->setParams (mp);
    const std::vector<float> got = render (*chain, in, frames, 1024);

    // THE REFERENCE: the engine on its own, driven in quanta of exactly K, on the gated input.
    auto engine = std::make_unique<eq::EqEngine>();
    test::run (engine->prepare (kFs, kK, kNch));
    for (int i = 0; i < eq::EqEngine::kMaxBands; ++i) engine->setBand (i, mp.eqBands[i]);
    std::vector<float> ref = in;
    gate (ref);
    float* p[core::kMaxChannels] {};
    for (int off = 0; off + kK <= frames; off += kK)
    {
        for (int c = 0; c < kNch; ++c) p[c] = ref.data() + (std::size_t) c * (std::size_t) frames + off;
        test::run (engine->process (p, kNch, kK));
    }
    const int whole = (frames / kK) * kK;   // the reference only covers whole quanta
    const std::vector<float> gotCut = head (got, frames, whole), refCut = head (ref, frames, whole);
    ok (sameBits (gotCut, refCut), "unarmed: the chain's EQ stage and eq::EqEngine agree on every sample ("
        + (differingSamples (gotCut, refCut)) + " differ)");

    // ...AND NO OTHER `dyn` FIELD MOVES A BIT WHILE `dyn.on` IS OFF.
    mastering::MasteringChainParams loud = mp;
    for (int b : { 0, 3, 7, 11 })
    {
        loud.eqBands[b].dyn.rangeDb = -30.0;
        loud.eqBands[b].dyn.thrDb   = -120.0;
        loud.eqBands[b].dyn.thrAuto = true;
        loud.eqBands[b].dyn.atk     = 0.0;
        loud.eqBands[b].dyn.rel     = 1.0;
    }
    auto other = std::make_unique<mastering::MasteringChain>();
    test::run (other->prepare (kFs, kNch, eqOnlyConfig()));
    other->setParams (loud);
    ok (sameBits (got, render (*other, in, frames, 1024)),
        "dyn.on off: range / threshold / mode / attack / release change no sample");
}

//==============================================================================
// 2. ARMED IS THE PRODUCER — `eq::EqBand` + `dynamiceq::LaneDynamics` built outside this chain.

// The oracle: ONE band and ONE producer, driven in quanta of K on the same signal, with the section
// input taken before the band runs — which for the first point IS the stage's input. Returns the audio
// and, in `delta`, the producer's own per-quantum gain reduction for the Stereo lane.
std::vector<float> oracle (const eq::BandParams& p, const std::vector<float>& in, int frames,
                           eq::Lane lane, std::vector<double>* delta = nullptr)
{
    eq::EqBand band;
    if (! band.prepare (kFs, kNch)) return {};
    band.setParams (p);
    dynamiceq::LaneDynamics dyn;
    if (! dyn.prepare (kFs, kNch)) return {};
    dyn.setParams (p);

    std::vector<float> out = in;
    gate (out);
    std::vector<float> key ((std::size_t) kK * (std::size_t) kNch, 0.0f);
    for (int off = 0; off + kK <= frames; off += kK)
    {
        float*       aud[core::kMaxChannels] {};
        const float* sc[core::kMaxChannels] {};
        for (int c = 0; c < kNch; ++c)
        {
            aud[c] = out.data() + (std::size_t) c * (std::size_t) frames + off;
            float* k = key.data() + (std::size_t) c * (std::size_t) kK;
            std::copy_n (aud[c], kK, k);
            sc[c] = k;
        }
        if (! dyn.processBand (aud, sc, kNch, kK, band)) return {};
        if (delta != nullptr) delta->push_back (dyn.deltaDb (lane));
    }
    return out;
}

void testArmedIsTheProducer()
{
    group ("an armed point is exactly what dynamiceq::LaneDynamics gives on the same signal");

    constexpr int frames = 24000;
    const std::vector<float> in = programme (Plan {}, frames);
    const eq::BandParams p = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);

    const std::vector<float> got = renderWith (p, in, frames);
    std::vector<double> delta;
    const std::vector<float> ref = oracle (p, in, frames, eq::Lane::Stereo, &delta);
    ok (! got.empty() && ! ref.empty(), "both the chain and the oracle rendered");

    const int whole = (frames / kK) * kK;
    const std::vector<float> a = head (got, frames, whole), b = head (ref, frames, whole);
    ok (sameBits (a, b), "armed: the chain and the standalone band + producer agree on every sample ("
        + (differingSamples (a, b)) + " differ)");

    // AND IT ACTUALLY MOVED, so the null above is not two transparent paths agreeing. The programme
    // spends 150 ms of every 500 ms 28 dB above the quiet stretch and the threshold sits between them.
    double deepest = 0.0, shallowest = -1.0e9;
    for (const double d : delta) { deepest = std::min (deepest, d); shallowest = std::max (shallowest, d); }
    ok (deepest < -6.0, "the loud stretch is cut (deepest " + std::to_string (deepest) + " dB)");
    ok (core::exactlyEqual (shallowest, 0.0),
        "and the quiet stretch is cut by EXACTLY nothing (shallowest " + std::to_string (shallowest) + " dB)");

    // The same point with its dynamics off renders differently — the cut is in the audio, not only in
    // the producer's readout.
    eq::BandParams off = p; off.dyn.on = false;
    const std::vector<float> flat = renderWith (off, in, frames);
    ok (! sameBits (got, flat), "the armed render differs from the same point with dyn.on off ("
        + (differingSamples (got, flat)) + " samples)");
}

//==============================================================================
// 2b. THE KEY IS THE SECTION INPUT — a point ahead of the armed one changes the signal it FILTERS and
// not the signal it DETECTS on. A chain that keyed each point on its own input would have the static
// cut below hide the burst from the armed point entirely.

// The oracle for that: `pre` runs statically over the quantum, then `p`'s producer runs on a key taken
// BEFORE `pre` touched anything.
std::vector<float> oracleBehind (const eq::BandParams& pre, const eq::BandParams& p,
                                 const std::vector<float>& in, int frames)
{
    eq::EqBand front, back;
    if (! front.prepare (kFs, kNch) || ! back.prepare (kFs, kNch)) return {};
    front.setParams (pre);
    back.setParams (p);
    dynamiceq::LaneDynamics dyn;
    if (! dyn.prepare (kFs, kNch)) return {};
    dyn.setParams (p);

    std::vector<float> out = in;
    gate (out);
    std::vector<float> key ((std::size_t) kK * (std::size_t) kNch, 0.0f);
    for (int off = 0; off + kK <= frames; off += kK)
    {
        float*       aud[core::kMaxChannels] {};
        const float* sc[core::kMaxChannels] {};
        for (int c = 0; c < kNch; ++c)
        {
            aud[c] = out.data() + (std::size_t) c * (std::size_t) frames + off;
            float* k = key.data() + (std::size_t) c * (std::size_t) kK;
            std::copy_n (aud[c], kK, k);          // the SECTION input, before `front` runs
            sc[c] = k;
        }
        if (! front.processBlock (aud, kNch, kK)) return {};
        if (! dyn.processBand (aud, sc, kNch, kK, back)) return {};
    }
    return out;
}

void testTheKeyIsTheSectionInput()
{
    group ("the detector key is the EQ section's input, not the point's own");

    constexpr int frames = 24000;
    const std::vector<float> in = programme (Plan {}, frames);

    // A deep static cut at the same frequency, AHEAD of the armed point. Its output is 24 dB down in
    // the band, which is under the armed point's threshold — so a point keyed on its own input would
    // never engage.
    eq::BandParams pre = bellPoint (kF0, 2.0, 0.0, -24.0, false, eq::Lane::Stereo);
    pre.lane (eq::Lane::Stereo).gainDb = -24.0;
    const eq::BandParams armed = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);

    auto chain = std::make_unique<mastering::MasteringChain>();
    test::run (chain->prepare (kFs, kNch, eqOnlyConfig()));
    mastering::MasteringChainParams mp;
    mp.eqBands[0] = pre;
    mp.eqBands[3] = armed;
    chain->setParams (mp);
    const std::vector<float> got = render (*chain, in, frames, 1024);
    const std::vector<float> ref = oracleBehind (pre, armed, in, frames);

    const int whole = (frames / kK) * kK;
    const std::vector<float> a = head (got, frames, whole), b = head (ref, frames, whole);
    ok (sameBits (a, b), "a point behind a static cut detects on the section input ("
        + (differingSamples (a, b)) + " differ)");

    // PRECONDITION: the two keys really do disagree here — the same point keyed on the cut signal is a
    // different render, so the null above is not two paths that happen to agree.
    eq::BandParams flat = armed; flat.dyn.on = false;
    mastering::MasteringChainParams noDyn = mp; noDyn.eqBands[3] = flat;
    auto other = std::make_unique<mastering::MasteringChain>();
    test::run (other->prepare (kFs, kNch, eqOnlyConfig()));
    other->setParams (noDyn);
    ok (! sameBits (got, render (*other, in, frames, 1024)),
        "PRECONDITION: the point behind the cut is engaged at all");
}

//==============================================================================
// 3. THE LANES ARE INDEPENDENT — each detects on its own domain and no other.

void testLanesAreIndependent()
{
    group ("the five lanes detect on their own domain");

    constexpr int frames = 24000;

    // (a) the loud stretch is in LEFT only. The Left lane cuts; the Right lane sees a programme that
    //     never crosses its threshold and therefore applies EXACTLY nothing.
    {
        Plan plan; plan.left = true; plan.right = false;
        const std::vector<float> in = programme (plan, frames);
        for (const eq::Lane l : { eq::Lane::Left, eq::Lane::Right })
        {
            const eq::BandParams armed = bellPoint (kF0, 2.0, -18.0, -24.0, true, l);
            eq::BandParams flat = armed; flat.dyn.on = false;
            const std::vector<float> a = renderWith (armed, in, frames);
            const std::vector<float> f = renderWith (flat, in, frames);
            const bool moved = ! sameBits (a, f);
            ok (moved == (l == eq::Lane::Left),
                std::string (l == eq::Lane::Left ? "the Left lane cuts on a left-only burst"
                                                 : "the Right lane renders bit-identically on a left-only burst"));
        }
    }

    // (b) the whole tone is in SIDE (R is L inverted). The Side lane cuts; Mid carries nothing at all.
    {
        Plan plan; plan.invertRight = -1.0;
        const std::vector<float> in = programme (plan, frames);
        for (const eq::Lane l : { eq::Lane::Side, eq::Lane::Mid })
        {
            const eq::BandParams armed = bellPoint (kF0, 2.0, -18.0, -24.0, true, l);
            eq::BandParams flat = armed; flat.dyn.on = false;
            const bool moved = ! sameBits (renderWith (armed, in, frames), renderWith (flat, in, frames));
            ok (moved == (l == eq::Lane::Side),
                std::string (l == eq::Lane::Side ? "the Side lane cuts on a side-only programme"
                                                 : "the Mid lane renders bit-identically on a side-only programme"));
        }
    }

    // (c) the five placements are five different renders on a programme every one of them can see.
    {
        Plan plan; plan.invertRight = 1.0; plan.right = true;
        const std::vector<float> in = programme (plan, frames);
        std::vector<std::vector<float>> r;
        for (const eq::Lane l : { eq::Lane::Stereo, eq::Lane::Left, eq::Lane::Right, eq::Lane::Mid, eq::Lane::Side })
            r.push_back (renderWith (bellPoint (kF0, 2.0, -18.0, -24.0, true, l), in, frames));
        int pairs = 0, same = 0;
        for (std::size_t i = 0; i < r.size(); ++i)
            for (std::size_t j = i + 1; j < r.size(); ++j)
            { ++pairs; if (sameBits (r[i], r[j])) ++same; }
        // Left and Right are the same render on a programme that is identical in both channels, and Mid
        // is the Stereo lane's own signal there — so the distinguishable set is what is asserted, pair
        // by pair, rather than "all ten differ".
        ok (pairs == 10, "PRECONDITION: ten lane pairs");
        ok (! sameBits (r[0], r[4]) && ! sameBits (r[1], r[4]) && ! sameBits (r[3], r[4]),
            "the Side lane is a different render from Stereo, Left and Mid");
        ok (! sameBits (r[0], r[1]), "the Stereo lane is a different render from Left");
        ok (same <= 3, "at most the analytically identical pairs agree (" + std::to_string (same) + " of 10)");
    }
}

//==============================================================================
// 4. BLOCK INVARIANCE — the caller's partition decides nothing, with a point armed.

void testBlockInvariance()
{
    group ("block invariance with a point armed");

    constexpr int frames = 20000;
    const std::vector<float> in = programme (Plan {}, frames);
    const eq::BandParams p = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);

    const std::vector<float> a = renderWith (p, in, frames, 1);
    int bad = 0;
    for (const int block : { 7, 63, 256, 1021, 4096, 65536 })
        if (! sameBits (a, renderWith (p, in, frames, block))) ++bad;
    ok (! a.empty() && bad == 0,
        "blocks 1 / 7 / 63 / 256 / 1021 / 4096 / 65536 give bit-identical output (" + std::to_string (bad) + " off)");
}

//==============================================================================
// 5. ZERO LATENCY — an armed point does not shift the programme and does not move a tap offset.

void testLatency()
{
    group ("the dynamics cost no latency");

    const eq::BandParams armed = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);
    eq::BandParams flat = armed; flat.dyn.on = false;

    mastering::MasteringChainConfig cfg;          // the FULL default topology, where latency is not just K
    auto chain = std::make_unique<mastering::MasteringChain>();
    test::run (chain->prepare (kFs, kNch, cfg));
    mastering::MasteringChainParams mp; mp.eqBands[0] = flat;
    chain->setParams (mp);
    const int latOff = chain->latencySamples();
    const mastering::MasteringChainResolved rOff = chain->resolved();

    mp.eqBands[0] = armed;
    chain->setParams (mp);
    const mastering::MasteringChainResolved rOn = chain->resolved();
    ok (chain->latencySamples() == latOff, "latencySamples() does not move when a point is armed");
    ok (rOn.latencySamples == rOff.latencySamples && rOn.compressorTapOffset == rOff.compressorTapOffset
        && rOn.limiterTapOffset == rOff.limiterTapOffset && rOn.internalBlock == rOff.internalBlock,
        "the declared latency and both tap offsets are unchanged");

    // AND THE AUDIO AGREES WITH THE DECLARATION: an impulse comes out where an unarmed chain puts it.
    constexpr int frames = 4096;
    std::vector<float> imp ((std::size_t) frames * 2, 0.0f);
    imp[100] = imp[(std::size_t) frames + 100] = 0.5f;
    auto find = [&] (const eq::BandParams& band) {
        const std::vector<float> out = renderWith (band, imp, frames);
        if (out.size() < (std::size_t) frames) return (long long) -2;   // a REFUSED render is not an answer
        for (std::size_t i = 0; i < (std::size_t) frames; ++i)
            if (! core::exactlyEqual (out[i], 0.0f)) return (long long) i;
        return (long long) -1;
    };
    const long long onAt = find (armed), offAt = find (flat);
    ok (onAt == 100 && offAt == 100,
        "the impulse leaves at input index 100 either way (armed " + std::to_string (onAt)
        + ", unarmed " + std::to_string (offAt) + ")");
}

//==============================================================================
// 6. RT SAFETY — process() allocates nothing with a point armed.

void testNoAllocation()
{
    group ("RT safety — an armed point allocates nothing in process()");

    auto chain = std::make_unique<mastering::MasteringChain>();
    mastering::MasteringChainConfig cfg;          // the full default topology
    ok (chain->prepare (kFs, kNch, cfg), "prepare (this one DOES allocate, by design)");

    mastering::MasteringChainParams mp;
    mp.eqBands[0] = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);
    mp.eqBands[5] = bellPoint (200.0, 1.0, 12.0, -30.0, true, eq::Lane::Side);
    chain->setParams (mp);

    constexpr int n = 512;
    std::vector<float> buf ((std::size_t) n * 2, 0.0f);
    for (int i = 0; i < n; ++i)
        buf[(std::size_t) i] = buf[(std::size_t) n + (std::size_t) i]
            = (float) (0.4 * std::sin (2.0 * core::kPi * kF0 * (double) i / kFs));
    float* p[2] { buf.data(), buf.data() + n };

    for (int w = 0; w < 4; ++w) test::run (chain->process (p, kNch, n));   // warm up outside the count

    // The toggle is counted too: arming and disarming is a per-block parameter and must not reach the heap.
    const long long before = alloc::count.load();
    bool good = true;
    for (int i = 0; i < 40; ++i)
    {
        mp.eqBands[0].dyn.on = (i % 2) == 0;
        chain->setParams (mp);
        good = chain->process (p, kNch, n) && good;
    }
    const long long after = alloc::count.load();
    okNoAlloc (after == before, "process() with armed points allocates nothing ("
               + std::to_string (after - before) + " allocations)");
    ok (good, "and every one of those calls was accepted");
}

//==============================================================================
// 7. THE BUDGET — the producers are inside the number `prepareBytes` publishes.

void testBudget()
{
    group ("the dynamics producers are inside the published budget");

    mastering::MasteringChainConfig cfg;          // the full default topology, EQ present
    mastering::MasteringChain::Storage st;
    ok (mastering::MasteringChain::storageFor (kFs, kNch, cfg, st), "the geometry is admitted");
    ok (st.dynBands == (std::size_t) eq::EqEngine::kMaxBands,
        "one producer per band is budgeted (" + std::to_string (st.dynBands) + ")");
    std::printf ("      the producers cost %zu B each, %zu B for the bank\n",
                 sizeof (dynamiceq::LaneDynamics), st.dynBands * sizeof (dynamiceq::LaneDynamics));

    mastering::MasteringChainConfig noEq = cfg; noEq.eq = false;
    mastering::MasteringChain::Storage stNoEq;
    ok (mastering::MasteringChain::storageFor (kFs, kNch, noEq, stNoEq), "and so is the same chain without the EQ");
    ok (stNoEq.dynBands == 0, "a chain with no EQ budgets no producer");
    ok (st.bytes() - stNoEq.bytes()
            == eq::EqEngine::objectBytes() + st.eqScratch.bytes()
               + (std::uint64_t) eq::EqEngine::kMaxBands * (std::uint64_t) sizeof (dynamiceq::LaneDynamics),
        "the EQ stage's whole cost is the engine, its scratch and the producers");

    // WHAT PREPARING IT ACTUALLY ASKS THE HEAP FOR, counted.
    const std::uint64_t budget = mastering::MasteringChain::prepareBytes (kFs, kNch, cfg);
    auto chain = std::make_unique<mastering::MasteringChain>();
    const long long before = alloc::bytes.load();
    const bool prepared = chain->prepare (kFs, kNch, cfg);
    const long long got = alloc::bytes.load() - before;
    ok (prepared, "the chain prepares");
    okNoAlloc (got == (long long) budget, "prepareBytes() is what a fresh prepare() allocates ("
               + std::to_string (got) + " B against " + std::to_string (budget) + ")");

    const std::uint64_t again = chain->reprepareBytes (kFs, kNch, cfg);
    const long long b2 = alloc::bytes.load();
    const bool ok2 = chain->prepare (kFs, kNch, cfg);
    const long long got2 = alloc::bytes.load() - b2;
    ok (ok2, "and re-prepares");
    okNoAlloc (again == 0u && got2 == 0, "re-preparing at the same geometry asks for nothing, and says so ("
               + std::to_string (got2) + " B)");
}

//==============================================================================
// 8. THE EDGES — a bypass that stops the stage stops the detectors with it, and a reset is a restart.

void testEdges()
{
    group ("the stage's edges carry the detectors with them");

    constexpr int frames = 24000;
    const std::vector<float> in = programme (Plan {}, frames);
    const eq::BandParams p = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);

    // TWO RENDERS OF THE SAME PROGRAMME THROUGH ONE CHAIN ARE IDENTICAL — `reset()` restarts the
    // producers as well as the filters, which is what makes a loudness search's second pass mean what
    // its first did.
    auto chain = std::make_unique<mastering::MasteringChain>();
    test::run (chain->prepare (kFs, kNch, eqOnlyConfig()));
    mastering::MasteringChainParams mp; mp.eqBands[0] = p;
    chain->setParams (mp);
    const std::vector<float> first = render (*chain, in, frames, 1024);
    chain->reset();
    chain->setParams (mp);
    const std::vector<float> second = render (*chain, in, frames, 1024);
    ok (sameBits (first, second), "two renders through one chain agree after reset() ("
        + (differingSamples (first, second)) + " differ)");

    // A BYPASSED STAGE IS TRANSPARENT, dynamics or not.
    mastering::MasteringChainParams byp = mp; byp.bypassEq = true;
    auto chB = std::make_unique<mastering::MasteringChain>();
    test::run (chB->prepare (kFs, kNch, eqOnlyConfig()));
    chB->setParams (byp);
    std::vector<float> gated = in; gate (gated);
    ok (sameBits (render (*chB, in, frames, 1024), gated), "a bypassed EQ stage passes the input through");

    // A BYPASS STRETCH STOPS THE DETECTORS. Two streams differing only in how loud their FIRST stretch
    // was, both bypassed through the middle one, must agree sample for sample on the third — a detector
    // that merely froze would replay the first stretch's reduction onto the third.
    auto tail = [&] (double firstLoudDb) {
        Plan hot; hot.loudDb = firstLoudDb; hot.quietDb = firstLoudDb - 2.0;
        const std::vector<float> a = programme (hot, 12000);
        Plan cold; cold.loudDb = -40.0; cold.quietDb = -44.0;
        const std::vector<float> c = programme (cold, 12000);

        auto ch = std::make_unique<mastering::MasteringChain>();
        if (! ch->prepare (kFs, kNch, eqOnlyConfig())) return std::vector<float> {};
        mastering::MasteringChainParams on; on.eqBands[0] = p;
        mastering::MasteringChainParams off = on; off.bypassEq = true;
        ch->setParams (on);   (void) render (*ch, a, 12000, 1024);   // 1. loud, active
        ch->setParams (off);  (void) render (*ch, c, 12000, 1024);   // 2. bypassed
        ch->setParams (on);
        return render (*ch, c, 12000, 1024);                          // 3. quiet, active again
    };
    const std::vector<float> afterHot = tail (-4.0), afterCool = tail (-30.0);
    ok (sameBits (afterHot, afterCool),
        "what happened before a bypass stretch does not reach the audio after it ("
        + (differingSamples (afterHot, afterCool)) + " differ)");
}

//==============================================================================
// 9. THE LOUDNESS SEARCH — an armed point is upstream of the node the search moves, so it sees the
// same signal at every drive. That is what keeps `y(g, c) = 10^(c/20) * y(d, 0)` (LoudnessSolver.h)
// and the monotonicity of the limiter's own statistic true with the dynamics live.

// One render through a chain whose limiter is present, with the PRE-LIMITER tap (taken before the gain
// node) and the limiter's gain-reduction trace handed back.
struct Probe
{
    std::vector<float> out, pre, gr;
    bool ok = false;
};

Probe probeChain (const eq::BandParams& band, const std::vector<float>& in, int frames,
                  double gainDb, double ceilingDbTp)
{
    Probe r;
    mastering::MasteringChainConfig cfg;
    cfg.internalBlock = kK;
    cfg.eq = true; cfg.limiter = true;
    cfg.monoBass = cfg.compressor = cfg.clipper = cfg.dither = false;   // dither is after the limiter and
                                                                        // does not scale — outside the law
    auto chain = std::make_unique<mastering::MasteringChain>();
    if (! chain->prepare (kFs, kNch, cfg)) return r;
    mastering::MasteringChainParams mp;
    mp.eqBands[0] = band;
    mp.preLimiterGainDb = gainDb;
    mp.limiter.ceilingDbTp = ceilingDbTp;
    mp.limiter.releaseMs = 100.0;
    chain->setParams (mp);

    const int D = chain->latencySamples();
    const std::size_t total = (std::size_t) frames + (std::size_t) D;
    const std::size_t cap = total + (std::size_t) kK;
    const int os = chain->tapOversampleFactor();

    std::vector<float> buf (total * (std::size_t) kNch, 0.0f);
    for (int c = 0; c < kNch; ++c)
        std::copy_n (in.data() + (std::size_t) c * (std::size_t) frames, frames,
                     buf.data() + (std::size_t) c * total);
    r.pre.assign (cap * (std::size_t) kNch, 0.0f);
    r.gr.assign (cap * (std::size_t) os, 0.0f);

    float* pre[core::kMaxChannels] {};
    float* p[core::kMaxChannels] {};
    for (int c = 0; c < kNch; ++c) { pre[c] = r.pre.data() + (std::size_t) c * cap; p[c] = buf.data() + (std::size_t) c * total; }
    mastering::MasteringChainTaps taps;
    taps.preLimiter    = pre;
    taps.frameCapacity = (int) cap;
    taps.limiterGrDb   = r.gr.data();
    taps.osCapacity    = (int) (cap * (std::size_t) os);
    if (! chain->process (p, kNch, (int) total, taps)) return r;

    r.out.assign ((std::size_t) frames * (std::size_t) kNch, 0.0f);
    for (int c = 0; c < kNch; ++c)
        std::copy_n (buf.data() + (std::size_t) c * total + (std::size_t) D, frames,
                     r.out.data() + (std::size_t) c * (std::size_t) frames);
    r.gr.resize ((std::size_t) taps.osWritten);
    r.ok = true;
    return r;
}

void testTheSearchIsUndisturbed()
{
    group ("the loudness search: an armed point is upstream of the node the search moves");

    constexpr int frames = 16000;
    // A programme the LIMITER can work on: the band the point rides is only part of it, so an 18 dB cut
    // there does not take the loudness with it.
    Plan plan; plan.quietDb = -22.0; plan.loudDb = -4.0; plan.bassDb = -8.0;
    const std::vector<float> in = programme (plan, frames);
    const eq::BandParams p = bellPoint (kF0, 2.0, -18.0, -24.0, true, eq::Lane::Stereo);

    // (a) THE PRE-LIMITER NODE DOES NOT MOVE WITH THE DRIVE, bit for bit. Everything upstream of the
    //     gain node — the dynamic band's detector included — is a constant of the search.
    const Probe base = probeChain (p, in, frames, 0.0, -1.0);
    ok (base.ok, "the probe renders");
    int movedPre = 0;
    double worstGr = 0.0, prevGr = 0.0;
    bool monotone = true;
    for (const double g : { 3.0, 6.0, 12.0, 18.0 })
    {
        const Probe q = probeChain (p, in, frames, g, -1.0);
        if (! q.ok || ! sameBits (base.pre, q.pre)) ++movedPre;
        // (b) AND THE LIMITER'S OWN STATISTIC IS MONOTONE IN THAT DRIVE, which is what the solver's
        //     constraints are read against.
        double worst = 0.0;
        for (const float v : q.gr) worst = std::min (worst, (double) v);
        if (worst > prevGr + 1.0e-9) monotone = false;
        prevGr = worst;
        worstGr = worst;
    }
    ok (movedPre == 0, "the pre-limiter tap is bit-identical at every drive (" + std::to_string (movedPre)
        + " of 4 moved)");
    ok (monotone, "the limiter's worst gain reduction is monotone non-increasing in the drive (deepest "
        + std::to_string (worstGr) + " dB)");
    ok (worstGr < -1.0, "PRECONDITION: the limiter is working at the top of that sweep");

    // (c) THE SCALE LAW, on the 3x3 grid LoudnessSolver.h states it over, with the point armed.
    double law = 0.0, peak = 0.0;
    for (const double c : { -1.0, -3.0, -6.0 })
        for (const double g : { 3.0, 6.0, 12.0 })
        {
            const Probe a = probeChain (p, in, frames, g, c);
            const Probe b = probeChain (p, in, frames, g - c, 0.0);
            if (! a.ok || ! b.ok) { ok (false, "the scale-law grid rendered"); return; }
            const double k = std::pow (10.0, c / 20.0);
            for (std::size_t i = 0; i < a.out.size(); ++i)
            {
                law  = std::max (law, std::fabs ((double) a.out[i] - k * (double) b.out[i]));
                peak = std::max (peak, std::fabs ((double) a.out[i]));
            }
        }
    ok (peak > 0.1, "PRECONDITION: the compared renders are not silence (" + std::to_string (peak) + ")");
    ok (law < 1.0e-5, "y(g,c) = 10^(c/20) y(d,0) holds with the dynamics live (worst "
        + std::to_string (law) + " at peak " + std::to_string (peak) + ")");
    std::printf ("      scale law with a dynamic point: worst %.3e at peak %.3f\n", law, peak);
}

}   // namespace

int main()
{
    testUnarmedIsTheEngine();
    testArmedIsTheProducer();
    testTheKeyIsTheSectionInput();
    testLanesAreIndependent();
    testBlockInvariance();
    testLatency();
    testNoAllocation();
    testBudget();
    testEdges();
    testTheSearchIsUndisturbed();
    return test::report();
}
