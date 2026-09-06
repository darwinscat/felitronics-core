// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THE TONE BANDS' COEFFICIENT RAMP HAS A LENGTH OF ITS OWN.
//
// It used to be the host's block: the same turn of the same knob travelled in 0.33 ms at a block of 16
// and 85 ms at 4096, and in 180 SECONDS on the whole-file call an offline driver makes. `bandRampLength`
// reads the length off the band instead — five ring times of its own poles, floored at 64 samples and
// capped at 2048 — so the move is the same however the caller cuts the stream.
//
// WHAT EACH TEST HERE IS FOR. The block-invariance test is the claim; the length tests are what stops it
// being satisfied by a constant. Every fixture asserts a precondition first, because each of these would
// pass on a dead one: a band that does not filter is trivially invariant, and two lengths that are equal
// because both hit the floor prove nothing about the a2 they were supposed to be read from.
#include <felitronics_test.h>
#include <felitronics/rigplayer/RigPlayer.h>
#include <felitronics/eq/EqBand.h>

#include <string>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;
using namespace felitronics::rigplayer;
using felitronics::test::ok;
using felitronics::test::group;

// Portable: MSVC's <cmath> lacks kPi by default, and the repo's own convention
// (io/tests/WavTests.cpp:25) is a local constant rather than _USE_MATH_DEFINES.
static constexpr double kPi = 3.14159265358979323846;
static constexpr double kFs = 48000.0;
// The rails are TIME, so the suite derives its expectations from the rate exactly as the header does.
static const int kRampMin = (int) (RigPlayer::kBandRampMinSeconds * kFs);
static const int kRampMax = (int) (RigPlayer::kBandRampMaxSeconds * kFs);

// The ring time this file claims the length is read from, computed independently of the header: the
// poles sit at radius sqrt(|a2|) and the impulse response decays as r^n, so -60 dB is at
// ln(1e-3)/ln(r) samples. If the two disagree, one of them is wrong and the test says which.
// Computed here from a1 AND a2, independently of the header. `sqrt(|a2|)` alone is the GEOMETRIC MEAN
// of the two pole magnitudes and is right only for a complex pair; a low-Q design has a REAL pair, and
// the mean then misses the slow root badly — `matched::highShelfQDb(1000, Q 0.05, +12 dB)` has poles
// 0.9908 and -0.2975, mean 5.3 samples against the slow root's 1053. The first version of this helper
// duplicated the header's own mistake and therefore endorsed it.
static double ringSamples (double a1, double a2)
{
    const double disc = a1 * a1 - 4.0 * a2;
    const double rho  = disc >= 0.0 ? (std::fabs (a1) + std::sqrt (disc)) * 0.5 : std::sqrt (std::fabs (a2));
    if (! (rho > 0.0) || ! (rho < 1.0)) return 0.0;
    return std::log (1.0e-3) / std::log (rho);
}

static void testLengthIsReadOffTheBand()
{
    group ("the ramp length is derived from the band's own poles, not chosen");

    // A high shelf at 10 kHz Q 0.8 — the switch position in the reference pack. Barely rings.
    const SectionBiquad fast = designSection (SectionKind::HighShelf, 10000.0, -9.0, 0.8, kFs);
    // A low shelf at 167 Hz Q 0.58 — the reference pack's bass band. Rings ~8 ms.
    const SectionBiquad slow = designSection (SectionKind::LowShelf, 167.0, 6.0, 0.58, kFs);
    // A resonant bell the pack SCHEMA permits (nothing in namz_rig_load.h bounds hz or q beyond
    // "positive") even though no pack ships one: 100 Hz at Q 8 rings for ~176 ms.
    const SectionBiquad ringing = designSection (SectionKind::Peak, 100.0, 9.0, 8.0, kFs);

    // PRECONDITION — the three really do have different memories, or the rest measures nothing.
    const double rFast = ringSamples (fast.a1, fast.a2), rSlow = ringSamples (slow.a1, slow.a2), rRing = ringSamples (ringing.a1, ringing.a2);
    ok (rFast < rSlow && rSlow < rRing,
        "PRECONDITION the three bands ring for different times: " + std::to_string (1000.0 * rFast / kFs) + " / "
        + std::to_string (1000.0 * rSlow / kFs) + " / " + std::to_string (1000.0 * rRing / kFs) + " ms");

    const int lFast = RigPlayer::bandRampLength (fast, fast, kFs);
    const int lSlow = RigPlayer::bandRampLength (slow, slow, kFs);
    const int lRing = RigPlayer::bandRampLength (ringing, ringing, kFs);

    ok (lFast == kRampMin,
        "a band that barely rings takes the FLOOR (" + std::to_string (lFast) + ")");
    ok (lSlow > kRampMin && lSlow < kRampMax,
        "a band that rings ~8 ms takes a length between the rails (" + std::to_string (lSlow) + " samples, "
        + std::to_string (1000.0 * lSlow / kFs) + " ms)");
    ok (lRing == kRampMax,
        "a band that rings 176 ms is capped by the LAG BUDGET, not served (" + std::to_string (lRing) + ")");
    // ...and the middle one really is ONE ring, computed here independently of the header.
    ok (std::fabs ((double) lSlow - rSlow) <= 1.0,
        "the length between the rails IS one ring time (" + std::to_string (lSlow) + " vs "
        + std::to_string (rSlow) + ")");
    // The number a 512-sample host was already producing, which is why this one is the least
    // surprising choice: the reference pack's bass band lands at 9.1 ms against that host's 10.7.
    ok (lSlow > 300 && lSlow < 600, "the reference pack's bass band lands near what a live host gave ("
                                    + std::to_string (1000.0 * lSlow / kFs) + " ms)");

    // The SLOWER of the two endpoints wins — a move from a fast band to a slow one has to respect the
    // slow one, and a fixture that only ever passed the same band twice would not see it.
    ok (RigPlayer::bandRampLength (fast, slow, kFs) == lSlow && RigPlayer::bandRampLength (slow, fast, kFs) == lSlow,
        "the length is set by the SLOWER endpoint, in either direction");

    // A REAL POLE PAIR — the case `sqrt(|a2|)` gets wrong and the suite had none of. A low-Q shelf is
    // the ordinary way to reach it, and packs are full of low-Q shelves.
    {
        const SectionBiquad lowQ = designSection (SectionKind::HighShelf, 1000.0, 12.0, 0.05, kFs);
        ok (lowQ.a1 * lowQ.a1 - 4.0 * lowQ.a2 > 0.0, "PRECONDITION a Q 0.05 shelf really has REAL poles");
        const double mean = std::log (1.0e-3) / std::log (std::sqrt (std::fabs (lowQ.a2)));
        const double dom  = ringSamples (lowQ.a1, lowQ.a2);
        ok (dom > 20.0 * mean, "PRECONDITION the two estimates disagree by more than 20x ("
                               + std::to_string (mean) + " vs " + std::to_string (dom) + ")");
        ok (RigPlayer::bandRampLength (lowQ, lowQ, kFs) > kRampMin,
            "a real-pole band is NOT collapsed to the floor by the geometric mean ("
            + std::to_string (RigPlayer::bandRampLength (lowQ, lowQ, kFs)) + " samples)");
        ok (std::fabs ((double) RigPlayer::bandRampLength (lowQ, lowQ, kFs) - dom) <= 1.0
            || RigPlayer::bandRampLength (lowQ, lowQ, kFs) == kRampMax,
            "...it is the DOMINANT root's ring, or the ceiling");
    }

    // THE RAILS ARE TIME. A ring is a duration — a 167 Hz shelf rings 9.08 ms at every rate — so the
    // rails have to be too, or the same pack gets a quarter of the lag budget on a 192 kHz session.
    {
        const SectionBiquad r96  = designSection (SectionKind::Peak, 100.0, 9.0, 8.0, 96000.0);
        const SectionBiquad r48  = designSection (SectionKind::Peak, 100.0, 9.0, 8.0, 48000.0);
        const int l96 = RigPlayer::bandRampLength (r96, r96, 96000.0);
        const int l48 = RigPlayer::bandRampLength (r48, r48, 48000.0);
        ok (l96 == (int) (RigPlayer::kBandRampMaxSeconds * 96000.0) && l48 == (int) (RigPlayer::kBandRampMaxSeconds * 48000.0),
            "PRECONDITION both rates put this band at the CEILING (" + std::to_string (l48) + " / " + std::to_string (l96) + ")");
        ok (std::fabs (1000.0 * l96 / 96000.0 - 1000.0 * l48 / 48000.0) < 0.05,
            "the ceiling is the same TIME at 48 and 96 kHz (" + std::to_string (1000.0 * l48 / 48000.0) + " vs "
            + std::to_string (1000.0 * l96 / 96000.0) + " ms)");
        const SectionBiquad f96 = designSection (SectionKind::HighShelf, 10000.0, -9.0, 0.8, 96000.0);
        const SectionBiquad f48 = designSection (SectionKind::HighShelf, 10000.0, -9.0, 0.8, 48000.0);
        ok (std::fabs (1000.0 * RigPlayer::bandRampLength (f96, f96, 96000.0) / 96000.0
                     - 1000.0 * RigPlayer::bandRampLength (f48, f48, 48000.0) / 48000.0) < 0.05,
            "...and so is the floor");
    }

    // Degenerate inputs cannot produce a length outside the rails.
    SectionBiquad flat {}; flat.b0 = 1.0;                     // no pole pair at all
    SectionBiquad unstable {}; unstable.b0 = 1.0; unstable.a2 = 4.0;
    ok (RigPlayer::bandRampLength (flat, flat, kFs) == kRampMin, "a2 = 0 takes the floor");
    ok (RigPlayer::bandRampLength (unstable, unstable, kFs) == kRampMax, "|a2| >= 1 takes the ceiling");
}

//==============================================================================
// The claim, driven through the REAL player. The first version of this test rendered the ramp with a
// hand-written replica of the mechanism — which passed while two mutations of the actual `runBands`
// (its length back to the call, and a restarted ramp beginning at the old start instead of where it
// had got to) sailed through. A replica cannot fail for the reasons the original does; this drives
// `RigPlayer::process` and nothing else.
//
// The knob edit lands on a FIXED ABSOLUTE SAMPLE for every slicing — the render is split there, and
// only the cutting inside each half varies. Otherwise the comparison would be between two different
// parameter timelines, which is not what is being claimed.
static namz::rig::Rig oneBandRig()
{
    namz::rig::Rig rig;
    rig.rigId = "ramp-rig"; rig.name = "Ramp"; rig.modeledBy = "the test";
    namz::rig::Stage st;
    st.kind = namz::rig::StageKind::Nam; st.rawKind = "nam"; st.slot = "pedal";
    st.make = "Darwin's Cat"; st.model = "Ramp"; st.gearType = "pedal";
    auto& d = st.device;
    d.family = "Ramp"; d.rigId = "ramp-rig"; d.slot = "pedal";
    namz::rig::Control g; g.name = "gain"; g.role = namz::rig::Role::Gain; g.values = { "150" }; g.sweep = 300;
    d.controls = { g };
    namz::rig::FileEntry f; f.id = "g150"; f.settings = { { "gain", "150" } };
    d.files = { f };
    // ONE knob, ONE band, and it is the reference pack's bass shelf: 167 Hz Q 0.58 travelling
    // [-12, +6] dB. That band is the one whose ramp length lands between the rails, so the test is
    // about the derived length and not about a floor that would hide it.
    namz::rig::Tone bass;
    bass.name = "bass"; bass.sweep = 300; bass.placement = "post"; bass.reference = "150"; bass.defaultValue = "150";
    namz::rig::Section ls; ls.kind = namz::rig::SectionKind::LowShelf; ls.hz = 167.0; ls.q = 0.58;
    ls.dbAtMin = -12.0; ls.dbAtMax = 6.0;
    bass.sections = { ls };
    st.tone = { bass };
    rig.chain = { st };
    return rig;
}

// A NAM that is a unity gain, so what comes out of the player is the band and nothing else.
static std::vector<std::byte> unityModel()
{
    const char* j = R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[1.0],"sample_rate":48000})";
    std::vector<std::byte> v; for (const char* p = j; *p; ++p) v.push_back((std::byte) *p);
    return v;
}

static void renderThroughPlayer (int block, double fromDial, double toDial, int editAt, int n,
                                 std::vector<float>& out)
{
    RigPlayer p;
    // THE SAME maxBlock FOR EVERY ROW. Only the CALL size varies — `maxBlock` is a configuration (it
    // sizes the internal chunking, the convolver's partitioning and the stage's buffers), so varying
    // it too would compare two different players rather than two slicings of one. The first version of
    // this fixture passed `block` here and every row differed from sample zero, before the edit.
    p.prepare (kFs, n, 1);
    p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
    for (int i = 0; i < 64; ++i) p.serviceHere();              // let the model land
    p.setDial ("bass", fromDial);
    out.assign ((std::size_t) n, 0.0f);
    for (int i = 0; i < n; ++i)
        out[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * kPi * 90.0 * (double) i / kFs)
                                      + 0.2 * std::sin (2.0 * kPi * 1400.0 * (double) i / kFs));
    auto span = [&] (int from, int to)
    {
        for (int off = from; off < to; off += block)
        {
            const int count = std::min (block, to - off);
            float* io[1] { out.data() + off };
            p.process (io, 1, count);
            p.serviceHere();
        }
    };
    // The pre-edit span, ALWAYS cut the same way and always in enough calls for everything else to
    // arrive: `nam::BlendLaw` moves its weight by at most 0.25 PER CALL, so a single long call leaves
    // the blend gain still climbing after the edit — which is what the first version of this fixture
    // measured instead of the band (18000 samples differing, worst 0.32, with the knob held still).
    auto settleThen = [&] (int from, int to)
    {
        const int calls = 32;                                  // >> 1/0.25, with room for the gain glides
        const int step  = std::max (1, (to - from) / calls);
        for (int off = from; off < to; off += step)
        {
            const int c = std::min (step, to - off);
            float* io[1] { out.data() + off };
            p.process (io, 1, c);
            p.serviceHere();
        }
    };
    // EVERYTHING BEFORE THE EDIT IS CUT THE SAME WAY, in one call, for every row. That is not tidiness:
    // `RigPlayer` is NOT slicing-invariant as a whole — `rampInto`'s gain glides and `nam::BlendLaw`'s
    // per-CALL weight step are two more mechanisms on the caller's clock, and they are the rest of this
    // task, not this change. Measured before this line existed: with the knob held perfectly still, two
    // slicings of the same stream differed on 23936 of 24000 samples, worst 0.33. Feeding the pre-edit
    // span identically leaves both runs in the SAME state at the edit — the gains settled, the blend
    // law arrived — so what the rows below compare is the band ramp and nothing else.
    settleThen (0, editAt);
    p.setDial ("bass", toDial);                                // the edit, at a fixed absolute sample
    span (editAt, n);
}

static void testTheMoveIsTheSameHoweverItIsCut()
{
    group ("a declared length makes the move independent of the host block — through RigPlayer::process");
    const int n = 24000, editAt = 6000;

    std::vector<float> ref; renderThroughPlayer (n, 0.0, 300.0, editAt, n, ref);   // whole file, one call

    // PRECONDITION 1 — the player is ALIVE: the render is not silence and not the input.
    {
        double peak = 0.0;
        for (float v : ref) peak = std::fmax (peak, std::fabs ((double) v));
        ok (peak > 0.05, "PRECONDITION the player produces signal (peak " + std::to_string (peak) + ")");
    }
    // PRECONDITION 2 — the knob move is REAL: the same render with the knob left alone differs.
    {
        std::vector<float> still; renderThroughPlayer (n, 0.0, 0.0, editAt, n, still);
        double worst = 0.0;
        for (int i = 0; i < n; ++i) worst = std::fmax (worst, std::fabs ((double) ref[(std::size_t) i] - (double) still[(std::size_t) i]));
        ok (worst > 0.02, "PRECONDITION the knob move changes the signal by " + std::to_string (worst));
    }
    // PRECONDITION 3 — THE PLAYER IS BLOCK-INVARIANT WITH NOTHING MOVING. Everything below compares two
    // slicings and blames the difference on the band ramp; if the player already differs with the knob
    // held still, that blame is misplaced and the rows below measure someone else's mechanism.
    {
        std::vector<float> a, b;
        renderThroughPlayer (n, 0.0, 0.0, editAt, n, a);
        renderThroughPlayer (64, 0.0, 0.0, editAt, n, b);
        long long diff = 0; double worst = 0.0; long long first = -1;
        for (int i = 0; i < n; ++i)
            if (a[(std::size_t) i] != b[(std::size_t) i])
            { ++diff; if (first < 0) first = i; worst = std::fmax (worst, std::fabs ((double) a[(std::size_t) i] - (double) b[(std::size_t) i])); }
        ok (diff == 0, "PRECONDITION with the knob STILL, and the pre-edit span cut identically, the "
                       "player is slicing-invariant — so the rows below can blame the band ramp ("
                       + std::to_string (diff) + " differ, first@" + std::to_string (first)
                       + ", worst " + std::to_string (worst) + ")");
    }

    // PRECONDITION 4 — the ramp spans several of the small blocks below, or "carries across calls" is untested.
    {
        const SectionBiquad a = designSection (SectionKind::LowShelf, 167.0, -12.0, 0.58, kFs);
        const SectionBiquad b = designSection (SectionKind::LowShelf, 167.0,   6.0, 0.58, kFs);
        const int L = RigPlayer::bandRampLength (a, b, kFs);
        ok (L > 128, "PRECONDITION this band's ramp is " + std::to_string (L) + " samples, longer than the small blocks");
    }

    for (int block : { 1, 7, 16, 64, 100, 256, 1024, 4096 })
    {
        std::vector<float> y; renderThroughPlayer (block, 0.0, 300.0, editAt, n, y);
        long long diff = 0; double worst = 0.0;
        for (int i = 0; i < n; ++i)
            if (y[(std::size_t) i] != ref[(std::size_t) i])
            { ++diff; worst = std::fmax (worst, std::fabs ((double) y[(std::size_t) i] - (double) ref[(std::size_t) i])); }
        ok (diff == 0, "block " + std::to_string (block) + " is bit-identical to the whole-file call ("
                       + std::to_string (diff) + " differ, worst " + std::to_string (worst) + ")");
    }

    // A SECOND EDIT ARRIVING MID-RAMP — the case where a restarted ramp must begin where the first one
    // had got to, not where it started. Same claim: the slicing must not decide it.
    {
        auto twoEdits = [&] (int block, std::vector<float>& o)
        {
            RigPlayer p;
            p.prepare (kFs, n, 1);                             // configuration fixed; only the cutting varies
            p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
            for (int i = 0; i < 64; ++i) p.serviceHere();
            p.setDial ("bass", 0.0);
            o.assign ((std::size_t) n, 0.0f);
            for (int i = 0; i < n; ++i)
                o[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * kPi * 90.0 * (double) i / kFs));
            auto span = [&] (int from, int to) {
                for (int off = from; off < to; off += block)
                { const int c = std::min (block, to - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
            };
            for (int off = 0; off < 6000; off += 6000 / 32)     // pre-edit: 32 calls, always the same
            { const int c = std::min (6000 / 32, 6000 - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
            p.setDial ("bass", 300.0);
            span (6000, 6000 + 120);          // 120 samples in — well inside the ramp
            p.setDial ("bass", 90.0);         // ...and the hand moves again
            span (6000 + 120, n);
        };
        std::vector<float> r2; twoEdits (n, r2);
        for (int block : { 1, 16, 64, 256, 1024 })
        {
            std::vector<float> y; twoEdits (block, y);
            long long diff = 0;
            for (int i = 0; i < n; ++i) if (y[(std::size_t) i] != r2[(std::size_t) i]) ++diff;
            ok (diff == 0, "a second edit 120 samples into the ramp: block " + std::to_string (block)
                           + " bit-identical (" + std::to_string (diff) + " differ)");
        }
    }
}

//==============================================================================
// A RAMP RESTARTED MID-FLIGHT MUST BEGIN WHERE IT HAD GOT TO. Comparing slicings cannot see this —
// starting the new ramp from the OLD start is equally wrong at every block size, so every row agrees
// with every other. Nor can a step metric: the filter STATE stays continuous across a coefficient
// jump, so the output has no discontinuity to find (measured — the step at the edit was inside the
// settled signal's own, mutation and all).
//
// WHAT DOES SEE IT is a second edit that sends the knob BACK. Take the dial from 0 to 300, interrupt
// 27 % of the way along, and send it back to 0: the correct ramp starts at ~81 and travels DOWN, while
// one that restarts from the old start is already at 0 and does not move at all. The two are then a
// whole third of the knob apart, and two STATIC renders — the dial held at 0, and held at 81 — say
// which one happened without reproducing the mechanism.
static void testARestartedRampDoesNotJump()
{
    group ("a second edit mid-ramp continues from where the first had got to");
    const int n = 12000, firstEdit = 4000;
    const int L = RigPlayer::bandRampLength (designSection (SectionKind::LowShelf, 167.0, -12.0, 0.58, kFs),
                                             designSection (SectionKind::LowShelf, 167.0,   6.0, 0.58, kFs), kFs);
    const int interrupt = firstEdit + L * 27 / 100;            // 27 % along, so "where it had got to" is far from both ends

    auto play = [&] (double startDial, bool edits, double firstTo, double secondTo, std::vector<float>& o)
    {
        RigPlayer p;
        p.prepare (kFs, n, 1);
        p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
        for (int i = 0; i < 64; ++i) p.serviceHere();
        p.setDial ("bass", startDial);
        o.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i)
            o[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 90.0 * (double) i / kFs));
        auto call = [&] (int from, int to)
        {
            for (int off = from; off < to; off += 64)
            { const int c = std::min (64, to - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
        };
        if (! edits) { call (0, n); return; }
        call (0, firstEdit);
        p.setDial ("bass", firstTo);
        call (firstEdit, interrupt);
        p.setDial ("bass", secondTo);
        call (interrupt, n);
    };

    // The two STATIC references: the dial held where the interrupted ramp had got to, and held at 0.
    const double partWay = 300.0 * (double) (interrupt - firstEdit) / (double) L;
    std::vector<float> atPartWay, atZero, real;
    play (partWay, false, 0, 0, atPartWay);
    play (0.0,     false, 0, 0, atZero);
    play (0.0,     true, 300.0, 0.0, real);

    auto rms = [] (const std::vector<float>& y, int from, int to)
    {
        double s = 0.0; for (int i = from; i < to; ++i) s += (double) y[(std::size_t) i] * y[(std::size_t) i];
        return std::sqrt (s / (double) (to - from));
    };
    const int w0 = interrupt + 2, w1 = interrupt + 24;         // right after the interrupt, inside the new ramp
    const double rPart = rms (atPartWay, w0, w1), rZero = rms (atZero, w0, w1), rReal = rms (real, w0, w1);

    // PRECONDITION — the two references are far apart, or "closer to one of them" says nothing.
    ok (std::fabs (rPart - rZero) > 0.02 * std::max (rPart, rZero),
        "PRECONDITION the two static references differ (" + std::to_string (rPart) + " vs " + std::to_string (rZero) + ")");
    ok (std::fabs (rReal - rPart) < std::fabs (rReal - rZero),
        "just after the interrupt the band is where the FIRST ramp had got to, not back at its start ("
        + std::to_string (rReal) + " against " + std::to_string (rPart) + " / " + std::to_string (rZero) + ")");
}

//==============================================================================
// THE RAMP NULLS AN INDEPENDENT ORACLE, THROUGH THE REAL PLAYER, ON BOTH CHANNELS.
//
// This is the test the rest of the suite turned out not to be. The diff-pass consilium ran three
// mutations against it and all three PASSED: (A) `bandLen_ = 0`, i.e. no ramp at all, just a jump with
// continuous filter state; (B) `takeBands` re-initialising every band on every publish, i.e. reset plus
// jump; (C) the ramp cursor shared across channels, so channel 1 skips the ramp — invisible because
// every other fixture here is MONO. Block-invariance is satisfied by any edit landing on a fixed
// absolute sample, a jump included, so it cannot see (A) or (B); and the "restarted ramp" check reads
// an RMS window where a jump's continuous state sits on the right side of the comparison anyway.
//
// What sees all three is a reference NULL: an oracle built from the primitives — the same designSection
// endpoints, bandRampLength, and the interpolation the header specifies — fed the SAME input through
// one eq::Biquad, with the PLAYER as the thing under test. The replica is legitimate here precisely
// because it is the reference and not the DUT.
static void testTheRampNullsAnOracle()
{
    group ("the ramp nulls an independent oracle through RigPlayer::process, on BOTH channels");
    const int n = 20000, edit1 = 6000, edit2 = 6000 + 120;
    const SectionBiquad cA = designSection (SectionKind::LowShelf, 167.0, -12.0, 0.58, kFs);   // dial 0
    const SectionBiquad cB = designSection (SectionKind::LowShelf, 167.0,   6.0, 0.58, kFs);   // dial 300
    const SectionBiquad cC = designSection (SectionKind::LowShelf, 167.0,  -4.8, 0.58, kFs);   // dial 90

    std::vector<float> in ((std::size_t) n);
    for (int i = 0; i < n; ++i)
        in[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * kPi * 90.0 * (double) i / kFs)
                                     + 0.2 * std::sin (2.0 * kPi * 1400.0 * (double) i / kFs));

    RigPlayer p;
    p.prepare (kFs, n, 2);                                     // STEREO — the only fixture here that is
    p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
    for (int i = 0; i < 64; ++i) p.serviceHere();
    p.setDial ("bass", 0.0);
    std::vector<float> L = in, R = in;
    auto run = [&] (int from, int to, int block) {
        for (int off = from; off < to; off += block)
        { const int c = std::min (block, to - off); float* io[2] { L.data() + off, R.data() + off }; p.process (io, 2, c); p.serviceHere(); }
    };
    run (0, edit1, 6000 / 32);                                 // settle the gains and the blend law
    p.setDial ("bass", 300.0); run (edit1, edit2, 64);
    p.setDial ("bass",  90.0); run (edit2, n, 64);

    auto lerp = [] (const SectionBiquad& f, const SectionBiquad& t, double u) {
        SectionBiquad c;
        c.b0 = f.b0 + (t.b0 - f.b0) * u; c.b1 = f.b1 + (t.b1 - f.b1) * u; c.b2 = f.b2 + (t.b2 - f.b2) * u;
        c.a1 = f.a1 + (t.a1 - f.a1) * u; c.a2 = f.a2 + (t.a2 - f.a2) * u; return c; };
    eq::Biquad bq; bq.setCoeffs (cA); bq.reset();
    std::vector<float> o ((std::size_t) n);
    SectionBiquad cur = cA, from = cA, to = cA; int pos = 0, len = 0;
    for (int i = 0; i < n; ++i)
    {
        if (i == edit1) { from = cur; to = cB; len = RigPlayer::bandRampLength (from, to, kFs); pos = 0; }
        if (i == edit2) { from = cur; to = cC; len = RigPlayer::bandRampLength (from, to, kFs); pos = 0; }
        if (pos < len) { ++pos; cur = lerp (from, to, (double) pos / (double) len); } else cur = to;
        bq.setCoeffs (cur);
        o[(std::size_t) i] = bq.processSample (in[(std::size_t) i]);
    }

    auto worstIn = [&] (const std::vector<float>& y, int a, int b) {
        double w = 0.0; for (int i = a; i < b; ++i) w = std::fmax (w, std::fabs ((double) y[(std::size_t) i] - (double) o[(std::size_t) i])); return w; };
    ok (RigPlayer::bandRampLength (cA, cB, kFs) > 128,
        "PRECONDITION the ramp is longer than the blocks (" + std::to_string (RigPlayer::bandRampLength (cA, cB, kFs)) + ")");
    const double pre = worstIn (L, edit1 - 1000, edit1);
    ok (pre < 1e-5, "PRECONDITION before the edit the player IS the band: worst " + std::to_string (pre));
    ok (worstIn (L, edit1, n) < 1e-5, "the LEFT channel nulls the oracle across both ramps (worst " + std::to_string (worstIn (L, edit1, n)) + ")");
    ok (worstIn (R, edit1, n) < 1e-5, "the RIGHT channel nulls it too (worst " + std::to_string (worstIn (R, edit1, n)) + ")");
    {   // PRECONDITION — the null has teeth: a jump would miss the oracle by a lot.
        eq::Biquad j; j.setCoeffs (cA); j.reset(); double wj = 0.0;
        for (int i = 0; i < n; ++i)
        {
            if (i == edit1) j.setCoeffs (cB);
            if (i == edit2) j.setCoeffs (cC);
            const double y = j.processSample (in[(std::size_t) i]);
            if (i >= edit1) wj = std::fmax (wj, std::fabs (y - (double) o[(std::size_t) i]));
        }
        ok (wj > 1e-2, "PRECONDITION a jump would miss the oracle by " + std::to_string (wj) + " — the null has teeth");
    }
}

//==============================================================================
// WHERE THE FLUSH SITS, not merely whether it exists. A version of this change ran the segment loop
// AFTER all the audio — the right number of flushes, every one of them on the same final state, which
// is the per-call defect wearing the grid's clothes. Deleting the flush is caught by the tail reaching
// exact zero; MISPLACING it is not, because at the end of a call the state is the same either way.
// What separates them is a call that spans many periods: with the flush inside, a whole-file call
// zeroes the tail on the grid exactly as a sliced one does; with it outside, the whole-file call
// carries the tail to the end untouched.
static void testTheFlushIsInsideTheCall()
{
    group ("the flush happens INSIDE the call, on the grid — not after all the audio");
    const int n = 30000;
    auto render = [&] (int block, std::vector<float>& o)
    {
        RigPlayer p;
        p.prepare (kFs, n, 1);
        p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
        for (int i = 0; i < 64; ++i) p.serviceHere();
        p.setDial ("bass", 300.0);
        o.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < 1200; ++i)                          // tone, then digital silence
            o[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 110.0 * (double) i / kFs));
        // The first 2000 samples ALWAYS in 32 equal calls: `nam::BlendLaw` moves its weight by at most
        // 0.25 per CALL, so an arm that starts with one long call leaves the blend gain climbing for the
        // whole render and this test would measure that instead of the flush (measured: 3228 samples).
        for (int off = 0; off < 2000; off += 2000 / 32)
        { const int c = std::min (2000 / 32, 2000 - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
        for (int off = 2000; off < n; off += block)
        { const int c = std::min (block, n - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
    };
    std::vector<float> whole, sliced;
    render (n, whole);          // the tail in ONE call spanning 437 grid periods
    render (100, sliced);       // ...against many, at a size that is not a multiple of the period

    // PRECONDITION 1 — the tail reaches EXACT zero in the sliced arm, i.e. the flush is doing something.
    long long zeros = 0;
    for (int i = 5000; i < n; ++i) if (sliced[(std::size_t) i] == 0.0f) ++zeros;
    ok (zeros > (n - 5000) / 2, "PRECONDITION the sliced tail reaches exact zero (" + std::to_string (zeros) + ")");
    // PRECONDITION 2 — the tone is real, so the tail is a decay and not silence from the start.
    double peak = 0.0;
    for (int i = 0; i < 1200; ++i) peak = std::fmax (peak, std::fabs ((double) sliced[(std::size_t) i]));
    ok (peak > 0.05, "PRECONDITION the tone is there (" + std::to_string (peak) + ")");

    long long diff = 0;
    for (int i = 0; i < n; ++i) if (whole[(std::size_t) i] != sliced[(std::size_t) i]) ++diff;
    ok (diff == 0, "a call spanning 437 grid periods flushes inside itself, exactly as 100-sample calls do ("
                   + std::to_string (diff) + " differ)");
}

//==============================================================================
// prepare() IS A RESTART, WHICH MEANS THE BANDS RESTART TOO. Cancelling a ramp in flight is not enough:
// `bandTo_` and `bandCur_` survive, `rebuildBands` republishes the same COUNT, and `takeBands`
// re-initialises only bands BEYOND that count — so without retiring them the first block after a
// restart GLIDES from the previous stream's coefficients. Measured before the fix: 5951 samples
// differing from a fresh player, worst 0.54.
static void testPrepareRestartsTheBands()
{
    group ("prepare() retires the bands — the next stream does not glide from the last one's");
    const int n = 16000;
    auto render = [&] (bool viaRestart, std::vector<float>& o)
    {
        RigPlayer p;
        p.prepare (kFs, n, 1);
        p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
        for (int i = 0; i < 64; ++i) p.serviceHere();
        if (viaRestart)
        {
            p.setDial ("bass", 0.0);                            // the OLD stream's position...
            std::vector<float> w (256, 0.2f);
            float* io[1] { w.data() };
            p.process (io, 1, 256);                             // ...actually applied
            p.prepare (kFs, n, 1);                              // and now: a new stream
            p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
            for (int i = 0; i < 64; ++i) p.serviceHere();
        }
        p.setDial ("bass", 300.0);
        o.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i)
            o[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * kPi * 90.0 * (double) i / kFs));
        for (int off = 0; off < n; off += 256)
        { const int c = std::min (256, n - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
    };
    std::vector<float> fresh, restarted;
    render (false, fresh);
    render (true,  restarted);
    // PRECONDITION — the two dial positions are far apart, so a stale glide between them would show.
    {
        // Measure the RESPONSE, not b0: the two shelf designs differ by 18 dB down low while their b0
        // values sit close together, so a b0 comparison is a precondition about the wrong quantity —
        // it failed here while the designs were 18 dB apart.
        const SectionBiquad lo = designSection (SectionKind::LowShelf, 167.0, -12.0, 0.58, kFs);
        const SectionBiquad hi = designSection (SectionKind::LowShelf, 167.0,   6.0, 0.58, kFs);
        const double dLo = sectionMagnitudeDb (lo, 40.0, kFs), dHi = sectionMagnitudeDb (hi, 40.0, kFs);
        ok (std::fabs (dHi - dLo) > 12.0, "PRECONDITION dial 0 and dial 300 are " + std::to_string (dHi - dLo)
                                          + " dB apart at 40 Hz — a stale glide between them would show");
    }
    long long diff = 0; double worst = 0.0;
    for (int i = 0; i < n; ++i)
        if (fresh[(std::size_t) i] != restarted[(std::size_t) i])
        { ++diff; worst = std::fmax (worst, std::fabs ((double) fresh[(std::size_t) i] - (double) restarted[(std::size_t) i])); }
    ok (diff == 0, "a player restarted after running at another dial position renders like a fresh one ("
                   + std::to_string (diff) + " differ, worst " + std::to_string (worst) + ")");
}

//==============================================================================
// prepare() RE-ANCHORS THE AUDIO-TIME GRID. Without it a re-prepared player inherits the phase the old
// stream left, so the law-8 flush lands somewhere else and two renders of the same programme differ.
// Invisible to any fixture whose tail never reaches the flush threshold, so this one has a tail.
static void testPrepareReAnchorsTheGrid()
{
    group ("prepare() re-anchors the grid — a re-prepared player renders like a fresh one");
    const int n = 20000, warm = 37;                 // deliberately NOT a multiple of the grid period

    auto render = [&] (bool warmFirst, std::vector<float>& o)
    {
        RigPlayer p;
        p.prepare (kFs, n, 1);
        p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
        for (int i = 0; i < 64; ++i) p.serviceHere();
        p.setDial ("bass", 300.0);
        if (warmFirst)                              // leave the grid at phase 37, then restart
        {
            std::vector<float> w ((std::size_t) warm, 0.1f);
            float* io[1] { w.data() };
            p.process (io, 1, warm);
            p.prepare (kFs, n, 1);
            p.load (oneBandRig(), [] (const std::string&) { return unityModel(); });
            for (int i = 0; i < 64; ++i) p.serviceHere();
            p.setDial ("bass", 300.0);
        }
        o.assign ((std::size_t) n, 0.0f);
        for (int i = 0; i < 1000; ++i)              // tone, then digital silence: the tail is the point
            o[(std::size_t) i] = (float) (0.5 * std::sin (2.0 * kPi * 120.0 * (double) i / kFs));
        // 100, NOT 512. The grid period is 64, so a call size that is a MULTIPLE of it puts the same
        // number of boundaries in every call whatever the phase — floor((0+512)/64) and
        // floor((37+512)/64) are both 8 — and the phase shift this test exists to detect becomes
        // invisible. The first version of this fixture used 512 and passed with the re-anchor deleted.
        for (int off = 0; off < n; off += 100)
        { const int c = std::min (100, n - off); float* io[1] { o.data() + off }; p.process (io, 1, c); p.serviceHere(); }
    };

    std::vector<float> fresh, reused;
    render (false, fresh);
    render (true,  reused);
    // PRECONDITION — the tail reaches EXACT zero, which only happens because the flush ran; without it
    // the phase would be unobservable and this test would pass on anything.
    long long zeros = 0;
    for (int i = 4000; i < n; ++i) if (fresh[(std::size_t) i] == 0.0f) ++zeros;
    ok (zeros > (n - 4000) / 2, "PRECONDITION the tail reaches exact zero (" + std::to_string (zeros) + ")");
    long long diff = 0;
    for (int i = 0; i < n; ++i) if (fresh[(std::size_t) i] != reused[(std::size_t) i]) ++diff;
    ok (diff == 0, "a player used for " + std::to_string (warm)
                   + " samples and re-prepared renders identically (" + std::to_string (diff) + " differ)");
}

//==============================================================================
// And the thing the length is FOR: a fast move overshoots what the same move done slowly produces.
// Measured on this mechanism at +14.88 dB for an instant jump; the floor has to hold that down.
static void testTheFloorHoldsDownTheOvershoot()
{
    group ("the floor holds down the transient a faster move adds — on a band that RAMPS");
    // ⚠️ NOT the reference pack's switch. That switch's two positions carry ZERO and ONE section, so
    // moving between them changes the BAND COUNT, and `takeBands` starts an appearing band AT its
    // coefficients with no ramp at all — an appearance is a hard step, and giving it one is a separate
    // policy (P26 follow-up), not something this floor covers. The first version of this test claimed
    // it did, on a fixture that kept one biquad alive and interpolated identity -> shelf: a transition
    // the player does not perform. What follows is the same shape as a band that DOES ramp, and it is
    // labelled as the synthetic it is.
    const SectionBiquad from = designSection (SectionKind::HighShelf, 10000.0,  0.0, 0.8, kFs);
    const SectionBiquad to   = designSection (SectionKind::HighShelf, 10000.0, -9.0, 0.8, kFs);
    const double probe = 14000.0;

    auto peakOver = [&] (int L)
    {
        eq::Biquad bq; bq.setCoeffs (from); bq.reset();
        const double w = 2.0 * kPi * probe / kFs;
        for (int i = 0; i < 48000; ++i) bq.processSample ((float) (0.5 * std::sin (w * (double) i)));
        double pk = 0.0, ipk = 0.0;
        for (int i = 0; i < 24000; ++i)
        {
            SectionBiquad c = to;
            if (i < L)
            {
                const double t = (double) (i + 1) / (double) L;
                c.b0 = from.b0 + (to.b0 - from.b0) * t; c.b1 = from.b1 + (to.b1 - from.b1) * t;
                c.b2 = from.b2 + (to.b2 - from.b2) * t; c.a1 = from.a1 + (to.a1 - from.a1) * t;
                c.a2 = from.a2 + (to.a2 - from.a2) * t;
            }
            bq.setCoeffs (c);
            const double ph = w * (double) (48000 + i);
            const double y  = (double) bq.processSample ((float) (0.5 * std::sin (ph)));
            const std::complex<double> H = eq::evalCoeffs (c, w);
            pk  = std::fmax (pk,  std::fabs (y));
            ipk = std::fmax (ipk, std::fabs (0.5 * std::abs (H) * std::sin (ph + std::arg (H))));
        }
        return 20.0 * std::log10 (pk / std::max (1e-12, ipk));
    };

    const double instant = peakOver (1);
    const double atFloor = peakOver (kRampMin);
    // PRECONDITION — the instant move DOES overshoot, or the floor is holding down nothing.
    ok (instant > 1.0, "PRECONDITION an instant move overshoots the ideal one by "
                       + std::to_string (instant) + " dB");
    ok (atFloor < 0.2, "at the floor the overshoot is gone (" + std::to_string (atFloor) + " dB)");
}

int main()
{
    std::printf ("felitronics::rigplayer — the tone bands' declared coefficient ramp\n");
    testLengthIsReadOffTheBand();
    testTheMoveIsTheSameHoweverItIsCut();
    testARestartedRampDoesNotJump();
    testTheRampNullsAnOracle();
    testTheFlushIsInsideTheCall();
    testPrepareRestartsTheBands();
    testPrepareReAnchorsTheGrid();
    testTheFloorHoldsDownTheOvershoot();
    return felitronics::test::report();
}
