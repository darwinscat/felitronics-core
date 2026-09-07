// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// LAW 11c — "A PAUSE IS SILENCE". One suite for the whole law, across every module that owns an address
// of it, because the law is one claim and splitting it per module would let a stage answer it its own way.
//
// THE CLAIM, and the only one worth testing: a call with `nch == 0, n > 0` leaves the stage's SHARED,
// one-per-instance state exactly where `n` samples of digital silence at a live width would have left it.
// Bit for bit, at the same call boundaries — the flush cadence of every stage here is once per call, so a
// gap cut into three pieces is compared against a silence cut into the same three pieces and never
// against one long one.
//
// WHY THE REFERENCE RUN IS BUILT THE WAY IT IS. The oracle is a SECOND INSTANCE fed real zeros at a live
// width, driven through the stage's public entry point — not a re-derivation of the recurrence, which
// would share the implementation's own mistakes. Where the stage has per-channel memory UPSTREAM of the
// shared detector (the sidechain SVFs of DeEsser/DynamicEqBand, the sidechain HPF of NoiseGate, the ST
// probe of LaneDynamics) the two runs would disagree for a reason that has nothing to do with this law:
// the gap DROPS those columns (law 11a) and the silence RINGS them down. So every fixture here charges
// the shared state through a channel whose per-channel path is then brought to rest, and asserts that
// precondition out loud before comparing. A fixture that skipped it would be measuring 11a, not 11c.
//
// AND THE PRECONDITIONS ARE ASSERTED, not assumed, in both directions: the shared observable must have a
// non-zero RANGE across the pause (or the comparison is between two constants and passes for a frozen
// stage too), and the per-channel path must be at rest (or the comparison fails for the wrong reason).
// Two fixtures were thrown away getting this right: one put the charge on the channel that never stops,
// and one probed the return at a level loud enough to open the gate on both runs. Both read 0.00 dB.

#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/TransientShaper.h>
#include <felitronics/dynamics/NoiseGate.h>
#include <felitronics/dynamics/GainReductionPath.h>
#include <felitronics/deesser/DeEsser.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/poweramp/PowerAmpStage.h>
#include <felitronics/multiband/MultibandCompressor.h>
#include <felitronics/core/Math.h>
#include <felitronics/eq/EqBand.h>

#include <felitronics_test.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using felitronics::test::run;

namespace
{
constexpr double kFs = 48000.0;

bool bitsEqual (float a, float b) noexcept { return core::sameBits (a, b); }
bool bitsEqual (double a, double b) noexcept { return std::memcmp (&a, &b, sizeof (double)) == 0; }

// The gap lengths. Deliberately NOT multiples of 16 (DynamicEqBand's control period), of 64 (the block
// sizes the older gap tests use) or of the lookahead: a length aligned to the mechanism's own period is
// the third recorded form of a blind fixture, and the counter shortcut in DynamicEqBand is exactly what
// such a length would hide.
const int kGaps[] = { 1, 2, 7, 15, 16, 17, 31, 63, 64, 65, 127, 480, 1000, 4801, 48000 };

// The collapse's own list adds the interval BETWEEN the two events that bound it: with a 5 ms Rms window
// the detector level crosses `kGainToDbFloor` at about 13 200 samples and the recurrence parks at about
// 23 609, and phase 2 is the stretch in between. `kGaps` jumps 4 801 -> 48 000 straight over it, so a
// silent loop that stopped one sample early was invisible: in Peak mode the level reaches zero in a
// single step, and by 48 000 both runs are long since parked on the same word.
const int kGapsCollapse[] = { 1, 2, 7, 15, 16, 17, 31, 63, 64, 65, 127, 480, 1000, 4801,
                              14000, 16000, 20000, 23000, 48000 };

// BRING A STAGE'S PER-CHANNEL PATH TO EXACT REST by feeding it real digital silence until its output is
// bit-zero on every sample of a whole block. It has to be MEASURED, not counted: the fixture's first
// version drained a fixed four blocks, which is enough for a band that CUTS and not for the same band
// BOOSTING — a +18 dB bell rings down about x80 per block, so it needs nine, and at four the residue was
// 1.1e-5 and produced a 1.4e-7 difference on the return that looked exactly like a defect in the code.
// Returns false if rest was not reached inside `maxBlocks`, and every caller asserts that.
template <class Stage>
bool drainToRest (Stage& a, Stage& b, int blockSize, int maxBlocks = 64)
{
    for (int k = 0; k < maxBlocks; ++k)
    {
        std::vector<float> z1 ((std::size_t) blockSize, 0.0f), z2 ((std::size_t) blockSize, 0.0f);
        float* io[2] = { z1.data(), z2.data() };
        std::vector<float> w1 ((std::size_t) blockSize, 0.0f), w2 ((std::size_t) blockSize, 0.0f);
        float* jo[2] = { w1.data(), w2.data() };
        if (! (a.process (io, 2, blockSize) && b.process (jo, 2, blockSize))) return false;
        bool rest = true;
        for (int i = 0; i < blockSize; ++i)
            rest = rest && bitsEqual (z1[(std::size_t) i], 0.0f) && bitsEqual (z2[(std::size_t) i], 0.0f)
                        && bitsEqual (w1[(std::size_t) i], 0.0f) && bitsEqual (w2[(std::size_t) i], 0.0f);
        if (rest && k > 0) return true;          // k > 0: one more block AFTER the first silent one
    }
    return false;
}

// A tone that charges a detector: full scale, so the shared state has somewhere to travel FROM.
void fillTone (std::vector<float>& v, int startSample, double hz, float amp)
{
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = amp * (float) std::sin (2.0 * core::kPi * hz * (double) (startSample + (int) i) / kFs);
}
} // namespace

//==============================================================================
// 1. The premise the collapse stands on. `GainReductionPath::advanceSilence` stops taking the log the
//    moment the detector level reaches `core::gainToDb`'s floor, because from there the conversion —
//    and the curve behind it — returns the same bits for every smaller level. Lower that floor and the
//    collapse is reading a level the clamp no longer hides. A static_assert cannot say this: `std::log10`
//    is not portably constexpr, so the constant could be pinned and the function still walk away from it.
static void floorPremise()
{
    group ("law 11c premise — gainToDb's floor is where the collapse thinks it is");
    const double at = core::gainToDb (core::kGainToDbFloor);
    ok (bitsEqual (core::gainToDb (0.0), at), "gainToDb(0) is the floor value, bit for bit");
    ok (bitsEqual (core::gainToDb (-1.0), at), "a negative level clamps to the floor too");
    ok (bitsEqual (core::gainToDb (1.0e-30), at), "every level under the floor gives the SAME bits");
    ok (bitsEqual (core::gainToDb (core::kGainToDbFloor * 0.5), at), "...including one just under it");
    ok (! bitsEqual (core::gainToDb (core::kGainToDbFloor * 2.0), at), "and a level ABOVE it does not");
    ok (std::fabs (at - (-240.0)) < 1.0e-9, "the floor is -240 dB, which is what the headers say");
}

//==============================================================================
// 2. The collapse against an honest loop. `advanceSilence` is the implementation; the oracle here is a
//    SECOND CONSTRUCTION — `n` plain calls of the same public entry point, with no phase split and no
//    fixed-point exit. If the two-phase shortcut ever stops being the same arithmetic, this is what says so.
static void collapseAgainstHonestLoop()
{
    group ("law 11c — the collapsed silence equals the honest per-sample loop, bit for bit");
    const dynamics::Mode modes[] = { dynamics::Mode::DownCompress, dynamics::Mode::UpCompress, dynamics::Mode::DownExpand };
    const dynamics::Detector dets[] = { dynamics::Detector::Peak, dynamics::Detector::Rms };
    int cases = 0;
    for (dynamics::Mode m : modes)
        for (dynamics::Detector d : dets)
            for (double thr : { -30.0, -300.0 })       // -300 makes DIGITAL SILENCE an ACTIVE sample
                for (int gap : kGapsCollapse)
                {
                    dynamics::GainReductionParams p;
                    p.mode = m; p.detector = d; p.thresholdDb = thr; p.ratio = 4.0; p.kneeDb = 6.0;
                    p.rangeDb = 24.0; p.attackMs = 1.0; p.releaseMs = 250.0; p.rmsWindowMs = 5.0;
                    dynamics::GainReductionPath fast, honest;
                    fast.prepare (kFs); honest.prepare (kFs);
                    fast.setParams (p); honest.setParams (p);
                    for (int i = 0; i < 2000; ++i) { const float x = 0.9f * (float) std::sin (2.0 * core::kPi * 300.0 * i / kFs); (void) fast.processSample (x); (void) honest.processSample (x); }
                    const float chargedLevel = fast.detectorLevel();
                    fast.advanceSilence (gap);
                    for (int i = 0; i < gap; ++i) (void) honest.processSample (0.0f);
                    ++cases;
                    if (! (bitsEqual (fast.valueDb(), honest.valueDb()) && bitsEqual (fast.detectorLevel(), honest.detectorLevel())))
                        ok (false, "collapse == honest loop, gap " + std::to_string (gap));
                    // ...and the state actually MOVED across the pause, or the comparison above is between
                    // two numbers that never travelled (blind-fixture form 5). The DETECTOR LEVEL is the
                    // honest witness here and the gain reduction is not: at `thresholdDb = -300` digital
                    // silence is an ACTIVE sample, so a railed follower legitimately sits still while the
                    // level underneath it falls.
                    if (gap >= 4801) ok (! bitsEqual (fast.detectorLevel(), chargedLevel),
                                         "the detector level moved across the pause (precondition)");
                }
    ok (true, "collapse == honest loop over " + std::to_string (cases) + " configurations");
}

//==============================================================================
// 3. Past the settling horizon a longer pause changes nothing — the fixed-point exit, stated as an
//    observable rather than as a timing. This is what makes a ten-hour gap cost what a settled one costs.
static void fixedPointIsIdempotent()
{
    group ("law 11c — past the horizon, a longer pause is the same state");
    dynamics::GainReductionParams p;
    p.detector = dynamics::Detector::Peak; p.thresholdDb = -30.0; p.ratio = 4.0;
    p.attackMs = 1.0; p.releaseMs = 50.0;
    dynamics::GainReductionPath a, b;
    a.prepare (kFs); b.prepare (kFs); a.setParams (p); b.setParams (p);
    for (int i = 0; i < 2000; ++i) { const float x = 0.9f * (float) std::sin (2.0 * core::kPi * 300.0 * i / kFs); (void) a.processSample (x); (void) b.processSample (x); }
    const float charged = a.valueDb();
    a.advanceSilence (2000000);
    b.advanceSilence (2000000000);                     // a thousand times longer: same answer, same cost
    ok (bitsEqual (a.valueDb(), b.valueDb()), "2e6 and 2e9 samples of pause leave the same gain reduction");
    ok (bitsEqual (a.detectorLevel(), b.detectorLevel()), "...and the same detector level");
    // ...and the pause DID something. Without this the group passes against an `advanceSilence` with an
    // empty body: two untouched instances hold the same charged state and agree perfectly.
    ok (std::fabs (a.valueDb()) < 0.001f && std::fabs (charged) > 1.0f,
        "precondition: the charged " + std::to_string (charged) + " dB actually decayed across the pause");
}

//==============================================================================
// 4. THE INVARIANT ITSELF, per stage: a zero-width gap against real digital silence at a live width,
//    at the SAME call boundaries, compared on the stage's own observable AND on the audio of the return.
//    Every fixture states its two preconditions.

// --- Compressor ---------------------------------------------------------------------------------
static void compressorInvariant()
{
    group ("law 11c — dynamics::Compressor: a gap equals silence of the same length");
    for (dynamics::Mode m : { dynamics::Mode::DownCompress, dynamics::Mode::UpCompress, dynamics::Mode::DownExpand })
      for (dynamics::LinkMode lk : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
        for (int gap : kGaps)
        {
            const int B = 128;
            dynamics::CompressorParams p;
            p.mode = m; p.link = lk; p.thresholdDb = -40.0; p.ratio = 4.0; p.kneeDb = 0.0;
            p.rangeDb = 24.0; p.attackMs = 1.0; p.releaseMs = 250.0; p.lookaheadMs = 0.0;
            dynamics::Compressor A, Bc;
            if (! (A.prepare (kFs, B, 2) && Bc.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
            A.setParams (p); Bc.setParams (p);
            std::vector<float> l (B), r (B);
            for (int k = 0; k < 24; ++k)
            {
                fillTone (l, k * B, 300.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                fillTone (l, k * B, 300.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bc.process (jo, 2, B));
            }
            const double charged = A.gainReductionDb();
            ok (std::fabs (charged) > 1.0, "precondition: the compressor is holding real gain reduction");
            // lookahead 0 => there is NO per-channel memory to be at rest; stated rather than assumed.
            ok (A.latencySamples() == 0, "precondition: no lookahead line, so no per-channel memory in play");

            for (int off = 0; off < gap; )
            {
                const int n = std::min (B, gap - off);
                float* io[2] = { nullptr, nullptr };  run (A.process (io, 0, n));
                std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
                float* jo[2] = { z1.data(), z2.data() }; run (Bc.process (jo, 2, n));   // the SAME boundaries
                off += n;
            }
            if (! bitsEqual ((float) A.gainReductionDb(), (float) Bc.gainReductionDb()))
                ok (false, "gap == silence, gap " + std::to_string (gap));
            if (gap >= 4801) ok (! bitsEqual ((float) A.gainReductionDb(), (float) charged),
                                 "precondition: the state travelled across the pause");
            // and the RETURN is bit-identical, which is the thing a listener hears
            std::vector<float> al (B), ar (B), bl (B), br (B);
            fillTone (al, 0, 300.0, 0.01f); ar = al; bl = al; br = al;
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B)); run (Bc.process (bi, 2, B));
            bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
            if (! same) ok (false, "the return after the gap is bit-identical to the return after silence, gap " + std::to_string (gap));
        }
    ok (true, "Compressor: gap == silence over every mode x link x gap length");
}

// --- TransientShaper ----------------------------------------------------------------------------
static void shaperInvariant()
{
    group ("law 11c — dynamics::TransientShaper: a gap equals silence of the same length");
    for (int gap : kGaps)
    {
        const int B = 128;
        dynamics::TransientShaperParams p;
        // NOT the defaults: they are inert (0 dB attack AND sustain), which is the fifth blind form.
        p.attackDb = 12.0; p.sustainDb = -9.0; p.threshold = 0.05; p.gainSmoothMs = 1.0;
        p.fastReleaseMs = 20.0; p.slowReleaseMs = 150.0;
        dynamics::TransientShaper A, Bs;
        if (! (A.prepare (kFs, B, 2) && Bs.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
        A.setParams (p); Bs.setParams (p);
        std::vector<float> l (B), r (B);
        for (int k = 0; k < 24; ++k)
        {
            for (int i = 0; i < B; ++i) { const int t = k * B + i; const float e = (t % 1200 < 40) ? 1.0f : 0.02f;
                                          l[(std::size_t) i] = e * (float) std::sin (2.0 * core::kPi * 300.0 * t / kFs); }
            r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
            for (int i = 0; i < B; ++i) { const int t = k * B + i; const float e = (t % 1200 < 40) ? 1.0f : 0.02f;
                                          l[(std::size_t) i] = e * (float) std::sin (2.0 * core::kPi * 300.0 * t / kFs); }
            r = l; float* jo[2] = { l.data(), r.data() }; run (Bs.process (jo, 2, B));
        }
        for (int off = 0; off < gap; )
        {
            const int n = std::min (B, gap - off);
            float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
            std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
            float* jo[2] = { z1.data(), z2.data() }; run (Bs.process (jo, 2, n));
            off += n;
        }
        std::vector<float> al (B), ar (B), bl (B), br (B), dry (B);
        fillTone (al, 0, 300.0, 0.2f); ar = al; bl = al; br = al; dry = al;   // the ACTUAL input, kept
        float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
        run (A.process (ai, 2, B)); run (Bs.process (bi, 2, B));
        bool same = true, moved = false;
        // `moved` compares against the buffer that was actually fed, not against a re-derivation of it:
        // `fillTone` narrows a float sine, and re-deriving it in double and narrowing afterwards differs
        // in 23 of 128 samples — enough to satisfy "it is shaping" at mix 0, where it is shaping nothing.
        for (int i = 0; i < B; ++i) { same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
                                      moved = moved || ! bitsEqual (bl[(std::size_t) i], dry[(std::size_t) i]); }
        if (! same) ok (false, "the return after the gap is bit-identical to the return after silence, gap " + std::to_string (gap));
        if (gap == kGaps[0]) ok (moved, "precondition: the shaper is actually shaping the return");
    }
    ok (true, "TransientShaper: gap == silence over every gap length");
}

// --- NoiseGate ----------------------------------------------------------------------------------
static void gateInvariant()
{
    group ("law 11c — dynamics::NoiseGate: a gap equals silence of the same length");
    // THE GATE'S OWN GAP LIST. `kGaps` is blind here: with a 20 ms detector release and a -40 dB
    // threshold the gate closes at about 4 858 samples and reaches the floor by about 19 300, so every
    // kGaps point sits either side of the only window where `coreGain` carries information — 28 of the
    // 30 configurations passed against a gate that was never given the gap at all. These four are inside
    // that window, and they are what makes a one-sample-early exit from the silent loop visible.
    const int gaps[] = { 1, 7, 17, 63, 127, 480, 4801, 5200, 6000, 8000, 12000, 19000, 48000 };
    for (int gap : gaps)
        for (bool on : { true, false })
        {
            const int B = 128;
            dynamics::NoiseGate::Config cfg;
            cfg.floorDb = -90.0f; cfg.envAttackMs = 1.0f; cfg.envReleaseMs = 20.0f;
            cfg.holdMs = 5.0f; cfg.closeMs = 300.0f;
            // A HIGH sidechain high-pass corner, and it is the fixture's hardest constraint rather than a
            // detail. The lanes have to reach EXACT rest before the pause (they sit upstream of the
            // shared detector, and the gap drops them while silence rings them down) — but at the shipped
            // 75 Hz that one-pole needs ~3 200 samples to underflow, which is 67 ms, and the gate's own
            // 20 ms release closes long before then. At 4 kHz the lane empties inside a single block and
            // the gate is still open, which is the only state in which both halves of this test are true
            // at once. The shipped corner is exercised by the 11a group below, where it belongs.
            cfg.sidechainHpHz = 4000.0f;
            dynamics::NoiseGate A, Bg;
            if (! (A.prepare (kFs, B, 2) && Bg.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
            A.setConfig (cfg); Bg.setConfig (cfg); A.seedEnabled (true); Bg.seedEnabled (true);
            std::vector<float> l (B), r (B);
            for (int k = 0; k < 24; ++k)
            {
                fillTone (l, k * B, 500.0, 0.8f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B, true, -40.0f));
                fillTone (l, k * B, 500.0, 0.8f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
            }
            ok (A.currentCoreGain() > 0.5f, "precondition: the gate is OPEN before the pause");
            // ...AND THE LANES BROUGHT TO REST. The sidechain high-pass is per-channel and sits UPSTREAM
            // of the shared detector: the width-2 run rings it down from h0 ~ 0.117 while the width-0 run
            // DROPS it, so without this the two disagree about law 11a and not about law 11c — measured,
            // -0.119 dB on `currentCoreGain` at a 5 200-sample gap. The blind grid above is what hid it.
            {
                for (int k = 0; k < 2; ++k)
                {
                    std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f);
                    float* io[2] = { z1.data(), z2.data() }; run (A.process (io, 2, B, true, -40.0f));
                    std::vector<float> w1 ((std::size_t) B, 0.0f), w2 ((std::size_t) B, 0.0f);
                    float* jo[2] = { w1.data(), w2.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
                }
                ok (A.currentCoreGain() > 0.5f, "precondition: the gate is STILL open after the drain");
            }
            for (int off = 0; off < gap; )
            {
                const int n = std::min (B, gap - off);
                float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n, on, -40.0f));
                std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
                float* jo[2] = { z1.data(), z2.data() }; run (Bg.process (jo, 2, n, on, -40.0f));
                off += n;
            }
            if (! (bitsEqual (A.currentCoreGain(), Bg.currentCoreGain()) && bitsEqual (A.currentGain(), Bg.currentGain())))
                ok (false, "gap == silence, gap " + std::to_string (gap) + (on ? " on" : " off"));
            // The RETURN probe is QUIET — below the open threshold — so only a gate that opened on a
            // ghost shows. A loud probe opens both and reads 0.00 dB whatever the state is.
            std::vector<float> al (B, 0.002f), ar (B, 0.002f), bl (B, 0.002f), br (B, 0.002f);
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B, on, -40.0f)); run (Bg.process (bi, 2, B, on, -40.0f));
            bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
            if (! same) ok (false, "the quiet return is bit-identical, gap " + std::to_string (gap));
            if (gap >= 19000 && on) ok (A.currentCoreGain() < 0.5f, "precondition: a long pause CLOSED the gate — the state moved");
            if (gap == 6000 && on) ok (A.currentCoreGain() > 0.001f && A.currentCoreGain() < 0.999f,
                                       "precondition: a 6 000-sample gap lands MID-CLOSE, where the observable is alive");
        }
    ok (true, "NoiseGate: gap == silence over every gap length and both enable states");
}

// --- DynamicEqBand ------------------------------------------------------------------------------
static void dynamicEqBandInvariant()
{
    group ("law 11c — dynamiceq::DynamicEqBand: a gap equals silence of the same length");
    for (dynamiceq::DynamicEqMode m : { dynamiceq::DynamicEqMode::CutWhenLoud, dynamiceq::DynamicEqMode::BoostWhenLoud, dynamiceq::DynamicEqMode::BoostWhenQuiet })
      for (dynamics::LinkMode lk : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
        for (int K : { 1, 4, 16 })
          for (int gap : kGaps)
          {
              const int B = 128;
              dynamiceq::DynamicEqBandParams p;
              p.mode = m; p.link = lk; p.coeffUpdatePeriod = K;
              p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
              p.attackMs = 1.0; p.releaseMs = 250.0;
              dynamiceq::DynamicEqBand A, Bd;
              if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
              A.setParams (p); Bd.setParams (p);
              // The charge level is PER MODE. A loud tone moves the two "when loud" modes and leaves
              // BoostWhenQuiet at zero — an upward band has nothing to lift when the signal is already
              // over the threshold — so a single amplitude would have made a third of this sweep a
              // comparison between two zeros. The precondition below is what caught it.
              const float amp = (m == dynamiceq::DynamicEqMode::BoostWhenQuiet) ? 0.004f : 0.9f;
              std::vector<float> l (B), r (B);
              for (int k = 0; k < 24; ++k)
              {
                  fillTone (l, k * B, 3000.0, amp); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                  fillTone (l, k * B, 3000.0, amp); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
              }
              // BRING THE PER-CHANNEL SIDECHAIN AND AUDIO FILTERS TO REST. They sit around the shared
              // detector, and the gap DROPS them (law 11a) while the silence RINGS THEM DOWN — without
              // this the two runs would disagree about law 11a, not about law 11c.
              if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest before the pause");
              const double charged = A.dynamicDeltaDb();
              if (gap == kGaps[0]) ok (std::fabs (charged) > 0.5, "precondition: the band is holding a real dynamic delta");
              for (int off = 0; off < gap; )
              {
                  const int n = std::min (B, gap - off);
                  float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
                  std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
                  float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, n));
                  off += n;
              }
              if (! bitsEqual (A.dynamicDeltaDb(), Bd.dynamicDeltaDb()))
                  ok (false, "gap == silence, K " + std::to_string (K) + " gap " + std::to_string (gap));
              std::vector<float> al (B), ar (B), bl (B), br (B);
              fillTone (al, 0, 3000.0, 0.05f); ar = al; bl = al; br = al;
              float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
              run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
              bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
              if (! same) ok (false, "the return is bit-identical, K " + std::to_string (K) + " gap " + std::to_string (gap));
          }
    ok (true, "DynamicEqBand: gap == silence over every mode x link x control period x gap length");
}

// THE COUNTER SHORTCUT, on its own and off phase. `coeffUpdatePeriod` is lowered mid-stream so `ksamp_`
// is left OUTSIDE [0, K), which is the one state in which the closed form and the loop disagree.
static void controlCounterPhase()
{
    group ("law 11c — DynamicEqBand: the control counter survives a pause, off phase and after a K change");
    for (int gap : { 1, 3, 7, 15, 16, 17, 33, 1000 })
        for (int K2 : { 1, 3, 4, 16 })
        {
            const int B = 64;
            dynamiceq::DynamicEqBandParams p;
            p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
            // A FAST release on purpose: with a slow one the follower is still moving for the whole gap
            // and the honest loop consumes all of it, so the closed-form remainder — the thing this group
            // is named after — never runs and deleting it would leave the group passing.
            p.attackMs = 1.0; p.releaseMs = 0.4; p.coeffUpdatePeriod = 16;
            dynamiceq::DynamicEqBand A, Bd;
            if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
            A.setParams (p); Bd.setParams (p);
            std::vector<float> l (B), r (B);
            for (int k = 0; k < 24; ++k)
            {
                fillTone (l, k * B, 3000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                fillTone (l, k * B, 3000.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
            }
            if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest before the pause");
            // leave ksamp_ off phase zero with a 5-sample call, then LOWER K under it
            std::vector<float> s1 (5, 0.0f), s2 (5, 0.0f);
            { float* io[2] = { s1.data(), s2.data() }; run (A.process (io, 2, 5)); }
            { std::vector<float> t1 (5, 0.0f), t2 (5, 0.0f); float* jo[2] = { t1.data(), t2.data() }; run (Bd.process (jo, 2, 5)); }
            p.coeffUpdatePeriod = K2; A.setParams (p); Bd.setParams (p);
            { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, gap)); }
            { std::vector<float> z1 ((std::size_t) gap, 0.0f), z2 ((std::size_t) gap, 0.0f); float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, gap)); }
            std::vector<float> al (B), ar (B), bl (B), br (B);
            fillTone (al, 0, 3000.0, 0.9f); ar = al; bl = al; br = al;
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
            bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
            if (! same) ok (false, "counter phase survives, gap " + std::to_string (gap) + " K " + std::to_string (K2));
        }
    ok (true, "DynamicEqBand: the control counter's phase after a pause matches silence, off phase and after a K change");
}

// --- DeEsser ------------------------------------------------------------------------------------
static void deEsserInvariant()
{
    group ("law 11c — deesser::DeEsser: a gap equals silence, in BOTH topologies");
    for (deesser::DeEsserMode m : { deesser::DeEsserMode::DynamicEq, deesser::DeEsserMode::SplitBand })
      for (bool listen : { false, true })
        for (dynamics::LinkMode lk : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
          for (int gap : kGaps)
          {
              const int B = 128;
              deesser::DeEsserParams p;
              p.mode = m; p.listen = listen; p.link = lk;
              p.fc = 7000.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 12.0;
              p.attackMs = 1.0; p.releaseMs = 250.0;
              deesser::DeEsser A, Bd;
              if (! (A.prepare (kFs, B, 2) && Bd.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
              A.setParams (p); Bd.setParams (p);
              std::vector<float> l (B), r (B);
              for (int k = 0; k < 24; ++k)
              {
                  fillTone (l, k * B, 7000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                  fillTone (l, k * B, 7000.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
              }
              if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest before the pause");
              const double charged = A.gainReductionDb();
              if (gap == kGaps[0]) ok (std::fabs (charged) > 0.3, "precondition: the de-esser is holding real gain reduction");
              for (int off = 0; off < gap; )
              {
                  const int n = std::min (B, gap - off);
                  float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
                  std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
                  float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, n));
                  off += n;
              }
              if (! bitsEqual (A.gainReductionDb(), Bd.gainReductionDb()))
                  ok (false, "gap == silence, mode " + std::to_string ((int) m) + " listen " + std::to_string ((int) listen) + " gap " + std::to_string (gap));
              std::vector<float> al (B), ar (B), bl (B), br (B);
              fillTone (al, 0, 7000.0, 0.05f); ar = al; bl = al; br = al;
              float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
              run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
              bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
              if (! same) ok (false, "the return is bit-identical, mode " + std::to_string ((int) m) + " gap " + std::to_string (gap));
          }
    ok (true, "DeEsser: gap == silence in both topologies, with listen on and off, both links");
}

// --- PowerAmpStage (the SEVENTH address) --------------------------------------------------------
static void powerAmpInvariant()
{
    group ("law 11c — poweramp::PowerAmpStage: the shared sag supply spends the pause");
    for (int gap : kGaps)
    {
        const int B = 128;
        // The observable is the RETURN AUDIO, because the sag supply has no public getter — and the
        // return is the thing a listener hears anyway. `C` is a COLD instance that never saw the loud
        // tone: it is the precondition, and without it this fixture could not tell a live sag rail from
        // an inert one (the `Voicing` defaults are all zero, so a stage with an unfilled voicing does
        // not move a single bit — the fifth blind form, and this repository has already paid for it).
        poweramp::PowerAmpStage A, Bp, C;
        poweramp::Voicing v;
        // A 20-SECOND sag recovery, which is not musical and is the point: the drain to exact rest takes
        // 512 blocks (~1.4 s), and a musical 150 ms recovery would have released the supply completely
        // before the pause even started — the fixture would then compare two rested rails and pass
        // against a stage that freezes the supply. A state defect is tested at settings where the state
        // is still VISIBLE; the musical setting belongs to a test about sound, not about this.
        v.sagMaxDroop = 0.35f; v.sagFastMs = 3.0f; v.sagRecoveryMs = 20000.0f; v.driveScale = 1.0f;
        poweramp::Params pp; pp.driveDb = 18.0f; pp.sag = 1.0f; pp.outputDb = 0.0f;
        A.prepare (kFs, B); Bp.prepare (kFs, B); C.prepare (kFs, B);
        A.setParams (pp, v); Bp.setParams (pp, v); C.setParams (pp, v);
        std::vector<float> l (B), r (B);
        for (int k = 0; k < 24; ++k)
        {
            fillTone (l, k * B, 120.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
            fillTone (l, k * B, 120.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bp.process (jo, 2, B));
        }
        // TWO OBSERVABLES, and the audio one is the load-bearing half. This stage's per-channel path DOES
        // reach exact rest on silence, but slowly: measured, its residue is 5.1e-08 after 64 silent
        // blocks, 8.0e-22 after 256 and exactly zero after 512 — so a 256-block drain would have compared
        // law 11a's drop against a ring-down that had not finished, and read a defect that is not there.
        // Past rest the return is bit-comparable, and it is the only thing that can see the thirteen
        // block-rate GLIDES: `sagDroop()` cannot, because it reads the supply and not the drive.
        if (! drainToRest (A, Bp, B, 1024)) ok (false, "precondition: the per-channel path reaches exact rest before the pause");
        // ...AND A GLIDE IN FLIGHT. The thirteen block-rate smoothers are snapped on the first block and
        // never moved again unless a parameter changes, so a fixture that sets the params once and then
        // pauses would pass against an implementation that freezes the glides — which is half of what
        // this stage's law-11c defect was. Move Drive and Output right before the gap so both runs enter
        // it mid-transition.
        poweramp::Params moved = pp; moved.driveDb = 3.0f; moved.outputDb = -8.0f; moved.sag = 0.3f;
        A.setParams (moved, v); Bp.setParams (moved, v); C.setParams (moved, v);
        const float chargedDroop = A.sagDroop();
        for (int off = 0; off < gap; )
        {
            const int n = std::min (B, gap - off);
            float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
            std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
            float* jo[2] = { z1.data(), z2.data() }; run (Bp.process (jo, 2, n));
            off += n;
        }
        if (! bitsEqual (A.sagDroop(), Bp.sagDroop()))
            ok (false, "the sag supply after a gap equals the sag supply after silence, gap " + std::to_string (gap));
        {
            std::vector<float> al ((std::size_t) B), ar ((std::size_t) B), bl ((std::size_t) B), br ((std::size_t) B);
            fillTone (al, 0, 120.0, 0.9f); ar = al; bl = al; br = al;
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B)); run (Bp.process (bi, 2, B));
            bool same = true, live = false;
            for (int i = 0; i < B; ++i) { same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
                                          live = live || std::fabs ((double) bl[(std::size_t) i]) > 1.0e-3; }
            if (! same) ok (false, "the return after a gap is bit-identical to the return after silence, gap " + std::to_string (gap));
            if (gap == kGaps[0]) ok (live, "precondition: the return actually carries audio");
        }
        if (gap <= 480) ok (! bitsEqual (Bp.sagDroop(), C.sagDroop()),
                            "precondition: the sag supply is genuinely charged — a cold stage answers differently");
        if (gap >= 48000) ok (Bp.sagDroop() < chargedDroop,
                            "...and a second of pause has visibly recovered the rail, exactly as silence does");
    }
    ok (true, "PowerAmpStage: the sag supply and its glides spend a pause exactly as silence does");
}

// --- LaneDynamics -------------------------------------------------------------------------------
// Its oracle is NOT a width-1 call: from a hot stereo state, width 1 leaves the Stereo probe's column 0
// alive while a width-0 falling edge drops every column. The honest reference is a width-2 run over real
// zeros whose probe columns have already been rung down.
static void laneDynamicsInvariant()
{
    group ("law 11c — dynamiceq::LaneDynamics: a gap runs the lanes on silence instead of disengaging");
    for (bool thrAuto : { false, true })
        for (int gap : { 7, 15, 16, 17, 63, 1000, 4801 })
        {
            const int B = 128;
            eq::EqBand bandA, bandB;
            eq::BandParams bp; bp.on = true; bp.type = eq::FilterType::Bell;
            bp.lane (eq::Lane::Stereo).on = true;
            bp.lane (eq::Lane::Stereo).freq = 3000.0;
            bp.lane (eq::Lane::Stereo).Q = 2.0;
            bp.lane (eq::Lane::Stereo).gainDb = 0.0;
            bp.dyn.on = true; bp.dyn.rangeDb = -12.0; bp.dyn.atk = 0.2; bp.dyn.rel = 0.6;
            bp.dyn.thrAuto = thrAuto; bp.dyn.thrDb = -40.0;
            if (! (bandA.prepare (kFs, 2) && bandB.prepare (kFs, 2))) { ok (false, "band prepare"); return; }
            bandA.setParams (bp); bandB.setParams (bp);
            dynamiceq::LaneDynamics A, Bl;
            if (! (A.prepare (kFs, 2) && Bl.prepare (kFs, 2))) { ok (false, "prepare"); return; }
            A.setParams (bp); Bl.setParams (bp);
            // CHARGED WITH BURSTS, not a steady tone. In `thrAuto` the estimator learns whatever is
            // stationary and the lane then correctly IDLES on it — a steady tone leaves the delta at
            // zero, and this sweep would have compared two zeros for half its cases. (The module's own
            // LaneDynamicsTests says the same thing in its own words; the precondition below is what
            // made the fixture admit it.) 250 ms of +12 dB burst per second, the same shape.
            std::vector<float> a1 (B), a2 (B), s1 (B), s2 (B);
            long phase = 0;
            for (int k = 0; k < 400; ++k)
            {
                for (int i = 0; i < B; ++i, ++phase)
                {
                    const double t = (double) phase / kFs;
                    const double amp = core::dbToGain (-24.0 + (std::fmod (t, 1.0) > 0.75 ? 12.0 : 0.0));
                    s1[(std::size_t) i] = (float) (amp * std::sin (2.0 * core::kPi * 3000.0 * t));
                }
                s2 = s1; a1 = s1; a2 = s1;
                float* aud[2] = { a1.data(), a2.data() }; const float* sc[2] = { s1.data(), s2.data() };
                run (A.processBand (aud, sc, 2, B, bandA));
                a1 = s1; a2 = s1;
                float* au2[2] = { a1.data(), a2.data() }; const float* sc2[2] = { s1.data(), s2.data() };
                run (Bl.processBand (au2, sc2, 2, B, bandB));
            }
            // ...and END ON A BURST. The training above stops wherever the phase happens to land, and it
            // landed between bursts: the lane had correctly idled and the delta was -0.029 dB, so the
            // comparison would have been between two numbers that had already arrived. A final 100 ms
            // of the raised level leaves the follower mid-duck without re-teaching the estimator, whose
            // averaging constant is two seconds.
            for (int k = 0; k < 40; ++k)
            {
                for (int i = 0; i < B; ++i, ++phase)
                {
                    const double t = (double) phase / kFs;
                    s1[(std::size_t) i] = (float) (core::dbToGain (-12.0) * std::sin (2.0 * core::kPi * 3000.0 * t));
                }
                s2 = s1; a1 = s1; a2 = s1;
                float* aud[2] = { a1.data(), a2.data() }; const float* sc[2] = { s1.data(), s2.data() };
                run (A.processBand (aud, sc, 2, B, bandA));
                a1 = s1; a2 = s1;
                float* au2[2] = { a1.data(), a2.data() }; const float* sc2[2] = { s1.data(), s2.data() };
                run (Bl.processBand (au2, sc2, 2, B, bandB));
            }
            // Ring the ST probe columns and the band's own filters down with REAL silence, both runs,
            // until the output is bit-zero — measured, not counted, for the reason drainToRest states.
            bool atRest = false;
            for (int k = 0; k < 64 && ! atRest; ++k)
            {
                std::fill (s1.begin(), s1.end(), 0.0f); s2 = s1; a1 = s1; a2 = s1;
                float* aud[2] = { a1.data(), a2.data() }; const float* sc[2] = { s1.data(), s2.data() };
                run (A.processBand (aud, sc, 2, B, bandA));
                std::vector<float> b1 (B, 0.0f), b2 (B, 0.0f), t1 (B, 0.0f), t2 (B, 0.0f);
                float* au2[2] = { b1.data(), b2.data() }; const float* sc2[2] = { t1.data(), t2.data() };
                run (Bl.processBand (au2, sc2, 2, B, bandB));
                bool rest = true;
                for (int i = 0; i < B; ++i) rest = rest && bitsEqual (a1[(std::size_t) i], 0.0f) && bitsEqual (b1[(std::size_t) i], 0.0f);
                if (rest && k > 0) atRest = true;
            }
            if (! atRest) ok (false, "precondition: the probe columns reach rest before the pause");
            const double charged = A.deltaDb (eq::Lane::Stereo);
            if (gap == 7) ok (std::fabs (charged) > 0.2, "precondition: the Stereo lane is holding a real delta (thrAuto "
                              + std::to_string ((int) thrAuto) + ", delta " + std::to_string (charged) + ")");
            {
                std::vector<float> ga1 ((std::size_t) gap, 0.0f), ga2 ((std::size_t) gap, 0.0f);
                float* aud[1] = { nullptr };
                run (A.processBand (aud, nullptr, 0, gap, bandA));
                std::vector<float> gb1 ((std::size_t) gap, 0.0f), gb2 ((std::size_t) gap, 0.0f);
                float* au2[2] = { gb1.data(), gb2.data() }; const float* sc2[2] = { gb1.data(), gb2.data() };
                run (Bl.processBand (au2, sc2, 2, gap, bandB));
            }
            if (! bitsEqual (A.deltaDb (eq::Lane::Stereo), Bl.deltaDb (eq::Lane::Stereo)))
                ok (false, "gap == silence on the Stereo lane, thrAuto " + std::to_string ((int) thrAuto) + " gap " + std::to_string (gap));
            if (gap >= 4801) ok (! bitsEqual (A.deltaDb (eq::Lane::Stereo), charged),
                                 "precondition: the lane's delta travelled across the pause");
        }
    ok (true, "LaneDynamics: a gap advances the lanes on silence, both threshold modes");
}

//==============================================================================
// 5. THE COMPOSITE. `MultibandProcessor` forwards a gap to every band — and used to forward it TWICE to a
//    BYPASSED band, once on the falling edge and once on the zero-width branch. While a band merely froze
//    that was invisible; under law 11c it is a DOUBLE CLOCK. A bypassed band's ballistics must spend the
//    gap exactly once.
static void multibandDoesNotDoubleClock()
{
    group ("law 11c — MultibandProcessor clocks a bypassed band ONCE per gap, not twice");
    const int B = 128, gap = 4801;
    multiband::MultibandCompressor<3> A, Bm;
    if (! (A.prepare (kFs, B, 2) && Bm.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    dynamics::CompressorParams cp;
    cp.thresholdDb = -40.0; cp.ratio = 6.0; cp.attackMs = 1.0; cp.releaseMs = 400.0; cp.kneeDb = 0.0;
    for (int b = 0; b < 3; ++b) { A.setBandParams (b, cp); Bm.setBandParams (b, cp); }
    std::vector<float> l (B), r (B);
    for (int k = 0; k < 30; ++k)
    {
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bm.process (jo, 2, B));
    }
    if (! drainToRest (A, Bm, B, 256)) ok (false, "precondition: the crossover tree reaches rest before the pause");
    A.setBandBypass (1, true); Bm.setBandBypass (1, true);        // A takes the falling edge WITH a bypassed band
    // A: stereo -> gap.  B: stereo -> the same audio time as REAL silence. The falling edge fires for A
    // on its first gap chunk, which is exactly where the duplicate call used to live.
    for (int off = 0; off < gap; )
    {
        const int n = std::min (B, gap - off);
        float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
        std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
        float* jo[2] = { z1.data(), z2.data() }; run (Bm.process (jo, 2, n));
        off += n;
    }
    std::vector<float> al (B), ar (B), bl (B), br (B);
    fillTone (al, 0, 300.0, 0.01f); ar = al; bl = al; br = al;
    float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
    run (A.process (ai, 2, B)); run (Bm.process (bi, 2, B));
    bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
    ok (same, "a gap through a multiband with a bypassed band equals the same time of silence");
}

//==============================================================================
// 6. THE TWO PLACES WHERE BEHAVIOUR AT `nch > 0` MOVED. Everything else in this change is bit-identical
//    to the base at a live width (verified separately, whole-tree, by hash); these two are not, and each
//    gets its own test rather than a sentence.

// (a) `MultibandProcessor` used to hand a BYPASSED band a row of NULL planes on a narrowing call. The
//     pointers are filled only in the band loop, which skips a bypassed band — so a band bypassed since
//     `prepare()` had never had them filled at all. Not a wrong number: a segfault, on untouched `main`.
static void multibandBypassedNarrowingDoesNotCrash()
{
    group ("law 11a — a bypassed band survives a narrowing call (it used to dereference null planes)");
    const int B = 96;
    multiband::MultibandCompressor<3> m;
    if (! m.prepare (kFs, B, 2)) { ok (false, "prepare"); return; }
    dynamics::CompressorParams cp;
    cp.thresholdDb = -30.0; cp.ratio = 5.0; cp.attackMs = 2.0; cp.releaseMs = 200.0;
    for (int b = 0; b < 3; ++b) m.setBandParams (b, cp);
    m.setBandBypass (1, true);                       // bypassed BEFORE the first call: never got planes
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int blk = 0; blk < 4; ++blk)
    {
        fillTone (l, blk * B, 300.0, 0.9f); r = l;
        float* io[2] = { l.data(), r.data() }; run (m.process (io, 2, B));
    }
    fillTone (l, 0, 300.0, 0.9f);
    float* io[1] = { l.data() };
    ok (m.process (io, 1, B), "narrowing 2 -> 1 with a band bypassed is accepted, not a crash");
    // ...and the surviving lane is finite and still processed
    bool finite = true, moved = false;
    for (int i = 0; i < B; ++i) { finite = finite && std::isfinite (l[(std::size_t) i]);
                                  moved = moved || std::fabs ((double) l[(std::size_t) i]) > 1.0e-6; }
    ok (finite, "the surviving lane is finite after the narrowing");
    ok (moved, "precondition: the surviving lane is actually carrying signal");
}

// (b) `NoiseGate` never dropped its per-lane sidechain high-pass for a lane that stopped — law 11a, which
//     every sibling has had since P18. The lane came back and its first sample measured `x - xPrev` against
//     audio from before the gap, which opens a closed gate on nothing.
static void gateDropsItsSidechainHighPass()
{
    group ("law 11a — NoiseGate drops a stopped lane's sidechain high-pass (89.99 dB before)");
    const int B = 128;
    dynamics::NoiseGate::Config cfg;
    cfg.floorDb = -90.0f; cfg.envAttackMs = 1.0f; cfg.envReleaseMs = 20.0f; cfg.holdMs = 5.0f; cfg.closeMs = 300.0f;
    // THE TONE RIDES THE LANE THAT STOPS and lane 0 stays digitally silent, so anything above zero on the
    // return came from the lane that left. `Bg` never narrows: it is the same audio time with lane 1 silent.
    dynamics::NoiseGate A, Bg;
    if (! (A.prepare (kFs, B, 2) && Bg.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    A.setConfig (cfg); Bg.setConfig (cfg); A.seedEnabled (true); Bg.seedEnabled (true);
    std::vector<float> z ((std::size_t) B, 0.0f), t ((std::size_t) B), a0 ((std::size_t) B), a1 ((std::size_t) B);
    for (int k = 0; k < 40; ++k)
    {
        fillTone (t, k * B, 500.0, 0.8f);
        a0 = z; a1 = t; float* io[2] = { a0.data(), a1.data() }; run (A.process (io, 2, B, true, -40.0f));
        a0 = z; a1 = t; float* jo[2] = { a0.data(), a1.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
    }
    for (int k = 0; k < 200; ++k)                    // A goes MONO; B keeps lane 1 and feeds it silence
    {
        a0 = z; float* io[1] = { a0.data() }; run (A.process (io, 1, B, true, -40.0f));
        a0 = z; a1 = z; float* jo[2] = { a0.data(), a1.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
    }
    // The RETURN probe is QUIET — below the open threshold — so only a gate that opened on a ghost shows.
    std::vector<float> pa ((std::size_t) B, 0.002f), pa1 ((std::size_t) B, 0.0f);
    std::vector<float> pb ((std::size_t) B, 0.002f), pb1 ((std::size_t) B, 0.0f);
    float* ai[2] = { pa.data(), pa1.data() }; float* bi[2] = { pb.data(), pb1.data() };
    run (A.process (ai, 2, B, true, -40.0f)); run (Bg.process (bi, 2, B, true, -40.0f));
    bool same = true; double worst = 0.0;
    for (int i = 0; i < B; ++i) { same = same && bitsEqual (pa[(std::size_t) i], pb[(std::size_t) i]);
                                  worst = std::fmax (worst, std::fabs ((double) pb[(std::size_t) i])); }
    ok (same, "the lane that left and came back does not re-open the gate on a ghost");
    ok (worst < 1.0e-5, "precondition: the gate really is CLOSED on the return, so an opening would show");
}

//==============================================================================
// 7. ROUND TWO. Everything below exists because a mutation of the real code SURVIVED round one. The
//    stand's score by round is in the report; this is the half of it that turned into tests.

// (a) AN EXTERNAL KEY IS STILL CONSUMED THROUGH A PAUSE. The programme stopped; the key did not. The
//     invariant says a gap equals the same call at a live width carrying silence, and THAT call runs the
//     detector on the key — so a pause that ignored it would not be the silence it is defined to equal.
static void compressorKeepsItsExternalKey()
{
    group ("law 11c — Compressor: a gap still consumes an external key");
    for (int gap : { 7, 63, 1000, 4801 })
    {
        const int B = 128;
        dynamics::CompressorParams p;
        p.thresholdDb = -40.0; p.ratio = 4.0; p.kneeDb = 0.0; p.attackMs = 1.0; p.releaseMs = 250.0;
        dynamics::Compressor A, Bc, Z;
        if (! (A.prepare (kFs, B, 2) && Bc.prepare (kFs, B, 2) && Z.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
        A.setParams (p); Bc.setParams (p); Z.setParams (p);
        std::vector<float> l (B), r (B), k1 (B), k2 (B);
        for (int blk = 0; blk < 8; ++blk)
        {
            fillTone (l, blk * B, 300.0, 0.9f); r = l; fillTone (k1, blk * B, 300.0, 0.9f); k2 = k1;
            float* io[2] = { l.data(), r.data() }; const float* key[2] = { k1.data(), k2.data() };
            run (A.process (io, 2, B, key, 2)); 
            fillTone (l, blk * B, 300.0, 0.9f); r = l;
            float* jo[2] = { l.data(), r.data() }; run (Bc.process (jo, 2, B, key, 2));
            fillTone (l, blk * B, 300.0, 0.9f); r = l;
            float* zo[2] = { l.data(), r.data() }; run (Z.process (zo, 2, B, key, 2));
        }
        // A: zero-width WITH the key still arriving.  Bc: a silent PROGRAMME at width 2, same key.
        // Z: zero-width with NO key — the control that proves the key is doing something.
        for (int off = 0; off < gap; )
        {
            const int n = std::min (B, gap - off);
            for (int i = 0; i < n; ++i) k1[(std::size_t) i] = 0.9f * (float) std::sin (2.0 * core::kPi * 300.0 * (off + i) / kFs);
            k2 = k1;
            const float* key[2] = { k1.data(), k2.data() };
            float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n, key, 2));
            std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
            float* jo[2] = { z1.data(), z2.data() }; run (Bc.process (jo, 2, n, key, 2));
            float* zo[2] = { nullptr, nullptr }; run (Z.process (zo, 0, n, nullptr, 0));
            off += n;
        }
        if (! bitsEqual ((float) A.gainReductionDb(), (float) Bc.gainReductionDb()))
            ok (false, "a keyed gap equals a keyed silence, gap " + std::to_string (gap));
        ok (! bitsEqual ((float) A.gainReductionDb(), (float) Z.gainReductionDb()),
            "precondition: the key is doing something — a keyless gap answers differently");
    }
    ok (true, "Compressor: a gap consumes its external key exactly as a silent block does");
}

// (b) THE CONTROL COUNTER'S CLOSED FORM, at a length no honest loop would ever be asked to walk. Two
//     things are under test that a short gap cannot reach: the widened arithmetic (`ksamp_ + rem` is an
//     int, and `rem` here is two billion), and the coefficient write that has to land in the remainder.
static void controlCounterAtHugeLengths()
{
    group ("law 11c — DynamicEqBand: the control counter closes correctly over an INT_MAX-sample pause");
    for (int K : { 3, 7, 16 })
        for (int pre : { 0, 1, 2, 5 })
        {
            const int B = 64;
            dynamiceq::DynamicEqBandParams p;
            p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
            p.attackMs = 1.0; p.releaseMs = 40.0; p.coeffUpdatePeriod = K;
            dynamiceq::DynamicEqBand A, Bd;
            if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
            A.setParams (p); Bd.setParams (p);
            std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
            for (int k = 0; k < 24; ++k)
            {
                fillTone (l, k * B, 3000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                fillTone (l, k * B, 3000.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
            }
            if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest");
            if (pre > 0)   // leave the counter off phase zero by a call of `pre` samples
            {
                std::vector<float> s1 ((std::size_t) pre, 0.0f), s2 ((std::size_t) pre, 0.0f);
                float* io[2] = { s1.data(), s2.data() }; run (A.process (io, 2, pre));
                std::vector<float> t1 ((std::size_t) pre, 0.0f), t2 ((std::size_t) pre, 0.0f);
                float* jo[2] = { t1.data(), t2.data() }; run (Bd.process (jo, 2, pre));
            }
            // A takes the whole INT_MAX in one zero-width call. B walks the SAME number of samples of
            // real silence, in chunks — and the two must land on the same counter phase, which the
            // returning audio is what makes visible. Two billion is chosen so `ksamp_ + rem` overflows an
            // int: the loop's own fixed-point exit is what makes it cheap enough to run in a test.
            // INT_MAX, not "a big number". `ksamp_ + rem` in `int` overflows only when the sum passes
            // 2147483647, and two billion plus a counter under sixteen does not: a mutation that drops
            // the widening survives every gap shorter than this one.
            { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 2147483647)); }
            {
                // The reference cannot literally walk two billion samples in a test, and it does not have
                // to: past the settling horizon the only thing still moving is the counter, so a length
                // CONGRUENT to 2e9 modulo K — and far past the horizon — lands on the same state. The
                // congruence is the whole point; an "about the same" length would silently change the
                // phase, which is exactly what this test is for.
                const long long M = 4000000LL + ((2147483647LL - 4000000LL) % (long long) K);
                std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f);
                for (long long done = 0; done < M; )
                {
                    const int n = (int) std::min (M - done, (long long) B);
                    std::fill (z1.begin(), z1.end(), 0.0f); std::fill (z2.begin(), z2.end(), 0.0f);
                    float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, n)); done += n;
                }
            }
            std::vector<float> al ((std::size_t) B), ar ((std::size_t) B), bl ((std::size_t) B), br ((std::size_t) B);
            fillTone (al, 0, 3000.0, 0.9f); ar = al; bl = al; br = al;
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
            bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
            if (! same) ok (false, "an INT_MAX-sample pause lands on the same counter phase, K " + std::to_string (K) + " pre " + std::to_string (pre));
        }
    ok (true, "DynamicEqBand: the counter's closed form matches the loop at a length that overflows an int");
}

// (c) THE COEFFICIENT WRITE IN THE REMAINDER. Once the continuous state is fixed the loop stops, and the
//     remaining samples are closed in form — including the `ksamp_ == 0` write, which is idempotent only
//     because `smooth` no longer moves. Dropping it leaves the audio filter designed for a delta the band
//     no longer holds. The gap here is long enough to reach the fixed point and the counter is placed so
//     the write MUST land inside the remainder.
static void coefficientWriteLandsInTheRemainder()
{
    group ("law 11c — DynamicEqBand: the coefficient write inside the skipped remainder is not lost");
    for (int K : { 4, 16, 64 })
    {
        const int B = 64;
        dynamiceq::DynamicEqBandParams p;
        p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
        p.attackMs = 1.0; p.releaseMs = 5.0;                    // FAST, so the fixed point arrives early
        p.coeffUpdatePeriod = K;
        dynamiceq::DynamicEqBand A, Bd;
        if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
        A.setParams (p); Bd.setParams (p);
        std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
        for (int k = 0; k < 24; ++k)
        {
            fillTone (l, k * B, 3000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
            fillTone (l, k * B, 3000.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
        }
        if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest");
        const double before = A.dynamicDeltaDb();
        ok (std::fabs (before) > 0.5, "precondition: the band is still holding a delta when the pause starts");
        // ONE call against ONE call. The law is stated at the SAME CALL BOUNDARIES, and here that is
        // load-bearing rather than pedantic: the once-per-call denormal flush is what turns the parked
        // subnormal into a real zero, so a single 300 000-sample gap against 4 688 blocked silent calls
        // would differ by exactly that flush — measured, `curGainDb_` -1.6815581571897805e-43 against 0 —
        // for a reason that is law 8's cadence and not law 11c's arithmetic.
        { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 300000)); }
        { std::vector<float> z1 (300000, 0.0f), z2 (300000, 0.0f); float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, 300000)); }
        ok (bitsEqual (A.dynamicDeltaDb(), Bd.dynamicDeltaDb()), "the delta after the skipped remainder matches silence, K " + std::to_string (K));
        ok (! bitsEqual (A.dynamicDeltaDb(), before), "precondition: the delta actually arrived during the pause");
        // ...and the FILTER agrees, which is what the coefficient write is for.
        std::vector<float> al ((std::size_t) B), ar ((std::size_t) B), bl ((std::size_t) B), br ((std::size_t) B);
        fillTone (al, 0, 3000.0, 0.9f); ar = al; bl = al; br = al;
        float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
        run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
        bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
        ok (same, "the audio filter is designed for the delta the band actually holds, K " + std::to_string (K));
    }
}

// (d) THE PAUSE OWES THE SAME LAW-8 FLUSH THE AUDIO PATH OWES. The silent recurrence parks on a SUBNORMAL
//     rather than on zero, so the once-per-call flush is what turns it into a real zero — and a pause that
//     skipped it would leave a subnormal in a feedback state and disagree with silence in the last bits.
static void thePauseFlushes()
{
    group ("law 11c — the pause carries the same once-per-call denormal flush as the audio path");
    for (int gap : { 20000, 60000, 200000 })
    {
        const int B = 64;
        dynamiceq::DynamicEqBandParams p;
        p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
        p.attackMs = 1.0; p.releaseMs = 30.0; p.coeffUpdatePeriod = 16;
        dynamiceq::DynamicEqBand A, Bd;
        if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
        A.setParams (p); Bd.setParams (p);
        std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
        for (int k = 0; k < 24; ++k)
        {
            fillTone (l, k * B, 3000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
            fillTone (l, k * B, 3000.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B));
        }
        if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest");
        // ONE zero-width call for the whole gap, against the SAME gap cut into blocks. The flush cadence
        // differs between them by construction — one flush against many — so this is not a bit-equality
        // claim, it is the claim that BOTH end at a state that behaves identically on the return, which
        // is what the flush exists to guarantee.
        { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, gap)); }
        {
            std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f);
            for (int off = 0; off < gap; off += B) { std::fill (z1.begin(), z1.end(), 0.0f); std::fill (z2.begin(), z2.end(), 0.0f);
                                                     float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, std::min (B, gap - off))); }
        }
        std::vector<float> al ((std::size_t) B), ar ((std::size_t) B), bl ((std::size_t) B), br ((std::size_t) B);
        fillTone (al, 0, 3000.0, 0.9f); ar = al; bl = al; br = al;
        float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
        run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
        bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
        ok (same, "one long gap and the same gap in blocks reach the same flushed state, gap " + std::to_string (gap));
    }
}

// (e) A POISONED STATE ENTERING A PAUSE. The silent loops break on two non-finite states in a row rather
//     than spinning for the whole gap; the stage must still come out finite, because law 8's flush is
//     what clears poison and the pause owes it.
// The pause owes law 8's flush, and the flush is OBSERVABLE: the silent recurrence parks on a subnormal,
// so a pause that skipped its flush leaves one in a feedback state. Past the horizon the delta must be
// EXACTLY zero — a subnormal there is the missing flush, and it is the one thing a numeric test can see
// of a law the architecture doc says has no gate at all.
static void thePauseLeavesNoSubnormal()
{
    group ("law 8 through law 11c — a long pause leaves the shared state at EXACTLY zero, not a subnormal");
    const int B = 64;
    dynamiceq::DynamicEqBandParams p;
    p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
    p.attackMs = 1.0; p.releaseMs = 5.0; p.coeffUpdatePeriod = 1;   // K = 1: every sample writes the seam
    dynamiceq::DynamicEqBand A;
    if (! A.prepare (kFs, 2)) { ok (false, "prepare"); return; }
    A.setParams (p);
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 24; ++k) { fillTone (l, k * B, 3000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B)); }
    ok (std::fabs (A.dynamicDeltaDb()) > 0.5, "precondition: the band is holding a delta before the pause");
    { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 300000)); }   // past the horizon
    { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 1)); }        // one more, to publish the flushed state
    ok (bitsEqual (A.dynamicDeltaDb(), 0.0), "the delta is exactly 0.0 after the pause, bit for bit");
}

// The same law-8 observable one storey up, on the two stages that publish a number the flush reaches.
// A pause that skipped its once-per-call flush leaves the follower parked on a SUBNORMAL instead of on
// zero; past the horizon these meters must read exactly +0.0. `TransientShaper` and `NoiseGate` have no
// equivalent window — their parked state is the detector's, which neither publishes — and that is stated
// here rather than left as an untested claim.
static void theGapFlushesWhereItCanBeSeen()
{
    group ("law 8 through law 11c — the meters read exactly +0.0 after a pause, not a subnormal");
    const int B = 128;
    {
        dynamics::CompressorParams p;
        p.thresholdDb = -40.0; p.ratio = 4.0; p.kneeDb = 0.0; p.attackMs = 1.0; p.releaseMs = 5.0;
        dynamics::Compressor c;
        if (! c.prepare (kFs, B, 2)) { ok (false, "prepare"); return; }
        c.setParams (p);
        std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
        for (int k = 0; k < 8; ++k) { fillTone (l, k * B, 300.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (c.process (io, 2, B)); }
        ok (std::fabs (c.gainReductionDb()) > 1.0, "precondition: the compressor is holding gain reduction");
        { float* io[2] = { nullptr, nullptr }; run (c.process (io, 0, 200000)); }
        ok (bitsEqual ((float) c.gainReductionDb(), 0.0f), "Compressor: the meter is exactly +0.0 after the pause");
    }
    {
        deesser::DeEsserParams p;
        p.mode = deesser::DeEsserMode::SplitBand; p.fc = 7000.0; p.thresholdDb = -40.0;
        p.ratio = 6.0; p.rangeDb = 12.0; p.attackMs = 1.0; p.releaseMs = 5.0;
        deesser::DeEsser d;
        if (! d.prepare (kFs, B, 2)) { ok (false, "prepare"); return; }
        d.setParams (p);
        std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
        for (int k = 0; k < 8; ++k) { fillTone (l, k * B, 7000.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (d.process (io, 2, B)); }
        ok (std::fabs (d.gainReductionDb()) > 0.3, "precondition: the de-esser is holding gain reduction");
        // TWO pauses, and the second is the one that reads. This stage publishes its meter from INSIDE
        // the loop, before the once-per-call flush, so after a single pause it holds the parked subnormal
        // — and so would a silent block, which is the whole claim. What the flush buys is that the NEXT
        // call starts from a real zero: a second pause writes `0 + c*(0 - 0)` and reads exactly +0.0,
        // where an unflushed follower would carry its subnormal straight through.
        { float* io[2] = { nullptr, nullptr }; run (d.process (io, 0, 200000)); }
        { float* io[2] = { nullptr, nullptr }; run (d.process (io, 0, 1)); }
        ok (bitsEqual ((float) d.gainReductionDb(), 0.0f), "DeEsser: the meter is exactly +0.0 after the pause");
    }
}

static void aPoisonedStateSurvivesAPause()
{
    group ("law 11c — a NaN reaching a stage before a pause does not survive it, and does not hang it");
    const int B = 64;
    dynamics::TransientShaperParams p;
    p.attackDb = 12.0; p.sustainDb = -9.0; p.threshold = 0.05; p.gainSmoothMs = 1.0;
    dynamics::TransientShaper s;
    if (! s.prepare (kFs, B, 2)) { ok (false, "prepare"); return; }
    s.setParams (p);
    std::vector<float> l ((std::size_t) B, 0.5f), r ((std::size_t) B, 0.5f);
    l[7] = std::numeric_limits<float>::quiet_NaN();
    { float* io[2] = { l.data(), r.data() }; run (s.process (io, 2, B)); }
    { float* io[2] = { nullptr, nullptr }; run (s.process (io, 0, 2000000000)); }   // must not spin
    std::vector<float> a ((std::size_t) B, 0.25f), b ((std::size_t) B, 0.25f);
    { float* io[2] = { a.data(), b.data() }; run (s.process (io, 2, B)); }
    bool finite = true; for (int i = 0; i < B; ++i) finite = finite && std::isfinite (a[(std::size_t) i]);
    ok (finite, "the stage is finite after a poisoned block and a two-billion-sample pause");
}

// (f) THE COMPOSITE'S DOUBLE CLOCK, at the width where it lived. The falling-edge call for a bypassed band
//     and the zero-width branch both reach for the same band; without the width guard the first chunk of
//     every gap advances that band by 2n.
static void multibandDoubleClockAtZeroWidth()
{
    group ("law 11c — MultibandProcessor clocks a BYPASSED band once per gap, not once per chunk");
    const int B = 64, gap = 48000;
    // A takes the gap in 750 chunks. `Bone` takes ONE chunk of the same length as A's first chunk and
    // then stops. The bypassed band is skipped on every chunk but the falling edge, so the two must
    // agree — and they would not if the `! fallingEdge` guard were missing, because A would then have
    // clocked the bypassed band by all 48 000 samples instead of by 64.
    multiband::MultibandCompressor<3> A, Bone;
    if (! (A.prepare (kFs, B, 2) && Bone.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    dynamics::CompressorParams cp;
    cp.thresholdDb = -45.0; cp.ratio = 8.0; cp.attackMs = 0.5; cp.releaseMs = 400.0; cp.kneeDb = 0.0;
    for (int b = 0; b < 3; ++b) { A.setBandParams (b, cp); Bone.setBandParams (b, cp); }
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 30; ++k)
    {
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bone.process (jo, 2, B));
    }
    if (! drainToRest (A, Bone, B, 256)) ok (false, "precondition: the crossover tree reaches rest before the pause");
    A.setBandBypass (1, true); Bone.setBandBypass (1, true);
    const double charged = A.bandGainReductionDb (1);
    ok (std::fabs (charged) > 1.0, "precondition: the bypassed band is holding real gain reduction");
    for (int off = 0; off < gap; off += B) { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, std::min (B, gap - off))); }
    { float* io[2] = { nullptr, nullptr }; run (Bone.process (io, 0, B)); }
    ok (bitsEqual ((float) A.bandGainReductionDb (1), (float) Bone.bandGainReductionDb (1)),
        "a 48 000-sample gap advances a bypassed band by ONE chunk, not by the whole gap");
    ok (! bitsEqual ((float) A.bandGainReductionDb (1), (float) charged),
        "precondition: that one chunk is visible — the band did move");
}

// ...and the residual this cannot remove: a bypassed band gets ONE chunk of the gap where live-width
// silence gives it none, because law 11 has no call that carries an edge without time. The test states
// the size of that residual so it cannot grow unnoticed.
static void bypassedBandResidualIsOneChunk()
{
    group ("law 11c — the bypassed-band residual is bounded to one chunk, and its size is stated");
    const int B = 64;
    multiband::MultibandCompressor<3> A, Bm;
    if (! (A.prepare (kFs, B, 2) && Bm.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    dynamics::CompressorParams cp;
    cp.thresholdDb = -45.0; cp.ratio = 8.0; cp.attackMs = 0.5; cp.releaseMs = 400.0; cp.kneeDb = 0.0;
    for (int b = 0; b < 3; ++b) { A.setBandParams (b, cp); Bm.setBandParams (b, cp); }
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 30; ++k)
    {
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        fillTone (l, k * B, 300.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bm.process (jo, 2, B));
    }
    if (! drainToRest (A, Bm, B, 256)) ok (false, "precondition: the crossover tree reaches rest before the pause");
    A.setBandBypass (1, true); Bm.setBandBypass (1, true);
    const double before = A.bandGainReductionDb (1);
    for (int off = 0; off < 48000; off += B) { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, std::min (B, 48000 - off))); }
    { std::vector<float> z1 (48000, 0.0f), z2 (48000, 0.0f);
      for (int off = 0; off < 48000; off += B) { std::fill (z1.begin(), z1.end(), 0.0f); std::fill (z2.begin(), z2.end(), 0.0f);
                                                 float* jo[2] = { z1.data(), z2.data() }; run (Bm.process (jo, 2, std::min (B, 48000 - off))); } }
    ok (bitsEqual ((float) Bm.bandGainReductionDb (1), (float) before),
        "live-width silence leaves a bypassed band FROZEN — that is what bypass means today");
    // THE ORACLE IS THE RECURRENCE, not a percentage window. A 5 % window accepts TWO chunks as easily
    // as one — at a 400 ms release two chunks move the reduction by 0.66 % — so it would pass against a
    // composite that clocked the bypassed band for the whole gap. `B` float steps of the follower's own
    // release, computed here from the parameters, is what one chunk actually is.
    const float c = (float) std::exp (-1.0 / (0.4 * kFs));       // releaseMs 400, the follower's coeff
    float expect = (float) before;
    for (int i = 0; i < B; ++i) expect = 0.0f + c * (expect - 0.0f);
    ok (bitsEqual ((float) A.bandGainReductionDb (1), expect),
        "a gap moves a bypassed band by exactly ONE chunk of release, bit for bit");
}

// (g) THE COUNTER'S WIDENING, at the only length that can see it. `ksamp_ + rem` in `int` overflows only
//     when the sum passes INT_MAX, and `rem` is `n` MINUS the samples the honest loop walked — so a gap
//     of INT_MAX on a band whose ballistics are still moving leaves `rem` a hundred thousand short of the
//     edge and the un-widened arithmetic survives. It is visible only when the loop breaks IMMEDIATELY,
//     i.e. when the shared state is already at its fixed point when the pause begins. That is not an
//     exotic state: it is every pause that follows another pause.
static void counterWideningOnASettledBand()
{
    group ("law 11c — DynamicEqBand: the counter widens correctly when the pause starts already settled");
    for (int K : { 3, 5, 7, 16 })
        for (int phase = 0; phase < K && phase < 4; ++phase)
        {
            const int B = 64;
            dynamiceq::DynamicEqBandParams p;
            p.freq = 3000.0; p.Q = 2.0; p.thresholdDb = -40.0; p.ratio = 6.0; p.rangeDb = 18.0;
            p.attackMs = 1.0; p.releaseMs = 0.4; p.coeffUpdatePeriod = K;
            dynamiceq::DynamicEqBand A, Bd;
            if (! (A.prepare (kFs, 2) && Bd.prepare (kFs, 2))) { ok (false, "prepare"); return; }
            A.setParams (p); Bd.setParams (p);
            std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
            for (int k = 0; k < 8; ++k) { fillTone (l, k * B, 3000.0, 0.9f); r = l;
                                          float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
                                          fillTone (l, k * B, 3000.0, 0.9f); r = l;
                                          float* jo[2] = { l.data(), r.data() }; run (Bd.process (jo, 2, B)); }
            if (! drainToRest (A, Bd, B)) ok (false, "precondition: the per-channel path reaches rest");
            // SETTLE the shared state first, so the pause below breaks on its first step and `rem` is the
            // whole of INT_MAX. Then place the counter at `phase`, which is what decides the overflow.
            { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 200000)); }
            { std::vector<float> z1 (200000, 0.0f), z2 (200000, 0.0f); float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, 200000)); }
            if (phase > 0)
            {
                float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, phase));
                std::vector<float> z1 ((std::size_t) phase, 0.0f), z2 ((std::size_t) phase, 0.0f);
                float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, phase));
            }
            { float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, 2147483647)); }
            {
                // CONGRUENT TO THE WHOLE OF A's REMAINING TIME, not to the gap alone: B has already
                // walked the same 200 000 settling samples, so the shortcut has to subtract them or the
                // two runs land on different counter phases. Written wrong the first time, and K = 5 and
                // K = 16 passed anyway because 200 000 is a multiple of both — the third blind form,
                // caught by the rows that are not.
                const long long M = 200000LL + ((2147483647LL - 200000LL) % (long long) K);
                std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f);
                for (long long done = 0; done < M; )
                { const int n = (int) std::min (M - done, (long long) B);
                  std::fill (z1.begin(), z1.end(), 0.0f); std::fill (z2.begin(), z2.end(), 0.0f);
                  float* jo[2] = { z1.data(), z2.data() }; run (Bd.process (jo, 2, n)); done += n; }
            }
            std::vector<float> al ((std::size_t) B), ar ((std::size_t) B), bl ((std::size_t) B), br ((std::size_t) B);
            fillTone (al, 0, 3000.0, 0.9f); ar = al; bl = al; br = al;
            float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
            run (A.process (ai, 2, B)); run (Bd.process (bi, 2, B));
            bool same = true; for (int i = 0; i < B; ++i) same = same && bitsEqual (al[(std::size_t) i], bl[(std::size_t) i]);
            if (! same) ok (false, "an INT_MAX pause from a settled state keeps the counter phase, K "
                                   + std::to_string (K) + " phase " + std::to_string (phase));
        }
    ok (true, "DynamicEqBand: the counter's widening holds at the length that overflows an int");
}

// (h) A CHILD'S REFUSAL ON THE NARROWING PATH REACHES THE CALLER. The verdict used to be `(void)`-cast
//     away, and the band loop skips a bypassed band, so a bypassed child that refused reported nothing
//     at all and the composite still returned true. Law 11's whole point is that a refusal is RETURNED.
static void aRefusedBypassedChildIsReported()
{
    group ("law 11 — a bypassed band that REFUSES the narrowing call is not silently accepted");
    const int B = 64;
    multiband::MultibandProcessor<dynamics::Compressor, 3> m;
    if (! m.prepare (kFs, B, 2, 0, [&] (dynamics::Compressor& c) { return c.prepare (kFs, B, 2, 50.0); }))
    { ok (false, "prepare"); return; }
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    fillTone (l, 0, 300.0, 0.5f); r = l;
    { float* io[2] = { l.data(), r.data() }; run (m.process (io, 2, B)); }
    m.setBandBypass (1, true);
    // Disarm band 1 on its own: `prepare` with a width of zero is refused, and law 11b says a refused
    // prepare leaves the object UNUSABLE. Every later call to it must therefore be refused too.
    ok (! m.band (1).prepare (kFs, B, 0, 50.0), "precondition: the band's own prepare(0) is refused");
    fillTone (l, 0, 300.0, 0.5f);
    float* io[1] = { l.data() };
    ok (! m.process (io, 1, B), "the composite reports the refusal instead of returning true");
}

// (i) THE FOLLOWER'S BRANCH IS PER SAMPLE, AND A CONSTANT TARGET DOES NOT FREEZE IT. `advanceConstant`
//     hands the same target to `process()` every sample precisely so the attack/release choice —
//     `|target| > |current|` — is made afresh each time; hoisting it out of the loop is the obvious
//     "optimisation" and it is wrong. It needs a target whose SIGN differs from the current reduction,
//     which happens when the mode changes under a live follower: every mode on its own is single-signed,
//     so no sweep over modes can reach it, and the whole suite passed against the hoisted version.
//     The trajectory below is chosen so the branch flips ON the zero crossing: release 0.5 takes -2 to
//     -0.5 (still |target| < |current|), then attack takes it to +1. A hoisted coefficient gives +0.25.
static void theFollowerBranchIsPerSample()
{
    group ("law 11c — a constant target does not freeze the follower's attack/release branch");
    // A release coefficient of EXACTLY 0.5f, so the arithmetic below is exact and the expected numbers
    // are the recurrence's own rather than a tolerance.
    const double relMs = 1000.0 / (kFs * std::log (2.0));
    for (double range2 : { 1.0, 4.0, 12.0 })
    {
        dynamics::GainReductionParams p;
        p.detector = dynamics::Detector::Peak; p.mode = dynamics::Mode::DownExpand;
        p.thresholdDb = -30.0; p.ratio = 4.0; p.kneeDb = 0.0; p.rangeDb = 2.0;
        p.attackMs = 0.0; p.releaseMs = relMs;                  // instant attack, half-per-sample release
        dynamics::GainReductionPath fast, honest;
        fast.prepare (kFs); honest.prepare (kFs);
        fast.setParams (p); honest.setParams (p);
        (void) fast.processSample (0.0f); (void) honest.processSample (0.0f);
        ok (fast.valueDb() < -1.0f, "precondition: the follower is holding a NEGATIVE reduction");
        // ...and now the target flips sign under it.
        p.mode = dynamics::Mode::UpCompress; p.rangeDb = range2;
        fast.setParams (p); honest.setParams (p);
        fast.advanceSilence (2);
        (void) honest.processSample (0.0f); (void) honest.processSample (0.0f);
        ok (bitsEqual (fast.valueDb(), honest.valueDb()),
            "the collapse crosses zero exactly as the honest loop does, range " + std::to_string (range2));
        ok (fast.valueDb() > 0.0f, "precondition: the trajectory really did cross zero — it is a boost now");
        // ...and further, in case the flip is one sample later at another range.
        for (int more : { 1, 3, 9, 40 })
        {
            fast.advanceSilence (more);
            for (int i = 0; i < more; ++i) (void) honest.processSample (0.0f);
            if (! bitsEqual (fast.valueDb(), honest.valueDb()))
                ok (false, "still equal after " + std::to_string (more) + " more, range " + std::to_string (range2));
        }
    }
    ok (true, "the follower's branch survives a sign change under a constant target");
}

// (j) THE GATE'S EXIT PREDICATE HAS FIVE TERMS, AND THREE OF THEM WERE NEVER EXERCISED. The silent loop
//     stops when the envelope, the VCA gain, the enable ramp, the HOLD COUNTER and the open flag all
//     stand still; drop any one and the loop can exit while that one is still moving. Two of the three
//     untested drops cost 90 dB under ordinary configuration, and neither is reachable from the settings
//     the group above uses.
static void theGatesExitPredicateNeedsAllOfIt()
{
    group ("law 11c — NoiseGate: the pause does not exit while the HOLD counter is still running");
    // A HOLD LONGER THAN THE DETECTOR'S PARK. At holdMs 1500 against a 20 ms release the envelope
    // reaches its fixed point around 31 000 samples while `hold` still has ~40 000 to count — so a
    // predicate that ignores `hold` sees a settled state and exits on the FIRST step of every call,
    // which makes the hold expire once per call instead of once per sample. Measured on that mutation:
    // the gap leaves the gate wide open at 1.0 where silence closes it to 3.16e-5, i.e. 90.00 dB.
    const int B = 128, gap = 120000;
    dynamics::NoiseGate::Config cfg;
    cfg.floorDb = -90.0f; cfg.envAttackMs = 1.0f; cfg.envReleaseMs = 20.0f;
    cfg.holdMs = 1500.0f; cfg.closeMs = 100.0f; cfg.sidechainHpHz = 4000.0f;
    dynamics::NoiseGate A, Bg;
    if (! (A.prepare (kFs, B, 2) && Bg.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    A.setConfig (cfg); Bg.setConfig (cfg); A.seedEnabled (true); Bg.seedEnabled (true);
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 24; ++k)
    {
        fillTone (l, k * B, 500.0, 0.8f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B, true, -40.0f));
        fillTone (l, k * B, 500.0, 0.8f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
    }
    for (int k = 0; k < 2; ++k)   // lanes to rest (4 kHz corner — see the group above)
    {
        std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f);
        float* io[2] = { z1.data(), z2.data() }; run (A.process (io, 2, B, true, -40.0f));
        std::vector<float> w1 ((std::size_t) B, 0.0f), w2 ((std::size_t) B, 0.0f);
        float* jo[2] = { w1.data(), w2.data() }; run (Bg.process (jo, 2, B, true, -40.0f));
    }
    ok (A.currentCoreGain() > 0.5f, "precondition: the gate is open, and its hold outlives the detector's park");
    for (int off = 0; off < gap; off += B)
    {
        const int n = std::min (B, gap - off);
        float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n, true, -40.0f));
        std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
        float* jo[2] = { z1.data(), z2.data() }; run (Bg.process (jo, 2, n, true, -40.0f));
    }
    ok (bitsEqual (A.currentCoreGain(), Bg.currentCoreGain()), "the hold counts the gap's samples, not its calls");
    ok (A.currentCoreGain() < 0.5f, "precondition: the hold DID expire inside the pause — the state moved");

    group ("law 11c — NoiseGate: the enable ramp advances per sample through a pause, not per call");
    // `on` TOGGLED DURING THE PAUSE, from a state where everything else has already parked. The enable
    // crossfade is then the only thing still moving, so a predicate that ignores it exits immediately
    // and the ramp advances one step per CALL. Measured on that mutation: 0.0159 against 1.0.
    dynamics::NoiseGate C, D;
    if (! (C.prepare (kFs, B, 2) && D.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    dynamics::NoiseGate::Config c2 = cfg; c2.holdMs = 5.0f; c2.enableMs = 25.0f;
    C.setConfig (c2); D.setConfig (c2); C.seedEnabled (true); D.seedEnabled (true);
    for (int k = 0; k < 24; ++k)
    {
        fillTone (l, k * B, 500.0, 0.8f); r = l; float* io[2] = { l.data(), r.data() }; run (C.process (io, 2, B, true, -40.0f));
        fillTone (l, k * B, 500.0, 0.8f); r = l; float* jo[2] = { l.data(), r.data() }; run (D.process (jo, 2, B, true, -40.0f));
    }
    for (int off = 0; off < 60000; off += B)   // settle everything with the gate ON
    {
        const int n = std::min (B, 60000 - off);
        float* io[2] = { nullptr, nullptr }; run (C.process (io, 0, n, true, -40.0f));
        std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
        float* jo[2] = { z1.data(), z2.data() }; run (D.process (jo, 2, n, true, -40.0f));
    }
    ok (C.currentGain() < 0.5f, "precondition: the gate is closed and everything else has parked");
    for (int off = 0; off < 4800; off += B)    // ...and NOW switch it off, mid-pause
    {
        const int n = std::min (B, 4800 - off);
        float* io[2] = { nullptr, nullptr }; run (C.process (io, 0, n, false, -40.0f));
        std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
        float* jo[2] = { z1.data(), z2.data() }; run (D.process (jo, 2, n, false, -40.0f));
    }
    ok (bitsEqual (C.currentGain(), D.currentGain()), "the enable ramp counts the gap's samples, not its calls");
    ok (C.currentGain() > 0.5f, "precondition: the ramp DID travel inside the pause — the state moved");
}

// (k) THE GAIN-REDUCTION TAP AT WIDTH ZERO. On the base a zero-width call left the tap untouched; it is
//     filled now, because the tap is a per-SAMPLE observable and a silent block of any width fills it
//     sample by sample. That is a contract change, so it gets a test rather than a sentence — and the
//     KEY is linked here rather than duplicated across channels, so the link is visible too.
static void theTapAndTheKeyAtWidthZero()
{
    group ("law 11c — Compressor: the gain-reduction tap is filled at width zero, and the key is LINKED");
    const int B = 128;
    dynamics::CompressorParams p;
    p.thresholdDb = -40.0; p.ratio = 4.0; p.kneeDb = 0.0; p.attackMs = 1.0; p.releaseMs = 250.0;
    p.link = dynamics::LinkMode::Max;
    dynamics::Compressor A, Bc;
    if (! (A.prepare (kFs, B, 2) && Bc.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    A.setParams (p); Bc.setParams (p);
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B), k1 ((std::size_t) B), k2 ((std::size_t) B);
    for (int blk = 0; blk < 8; ++blk)
    {
        fillTone (l, blk * B, 300.0, 0.9f); r = l;
        float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        fillTone (l, blk * B, 300.0, 0.9f); r = l;
        float* jo[2] = { l.data(), r.data() }; run (Bc.process (jo, 2, B));
    }
    // THE TWO KEY CHANNELS DIFFER, and the loud one is channel 1. A `Max` link must see it; a width-0
    // path that read only `key[0]` would follow the quiet channel and never notice.
    std::fill (k1.begin(), k1.end(), 0.0f);
    fillTone (k2, 0, 300.0, 0.9f);
    const float* key[2] = { k1.data(), k2.data() };
    std::vector<float> tapA ((std::size_t) B, 12345.0f), tapB ((std::size_t) B, 12345.0f);
    { float* io[2] = { nullptr, nullptr };
      run (A.process (io, 0, B, key, 2, dynamics::GainReductionTap { tapA.data(), B })); }
    { std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f); float* jo[2] = { z1.data(), z2.data() };
      run (Bc.process (jo, 2, B, key, 2, dynamics::GainReductionTap { tapB.data(), B })); }
    bool same = true, written = true, varies = false;
    for (int i = 0; i < B; ++i)
    {
        same    = same && bitsEqual (tapA[(std::size_t) i], tapB[(std::size_t) i]);
        written = written && ! bitsEqual (tapA[(std::size_t) i], 12345.0f);
        varies  = varies || ! bitsEqual (tapA[(std::size_t) i], tapA[0]);
    }
    ok (written, "the tap is WRITTEN at width zero — every slot, not just the first");
    ok (same, "and it holds the same samples a silent block of the same length would write");
    ok (varies, "precondition: the tap carries a per-sample trajectory, not one repeated number");
    ok (bitsEqual ((float) A.gainReductionDb(), (float) Bc.gainReductionDb()), "and the meter agrees");
}

// (l) WHAT the bypassed band was clocked ON. The narrowing fix has two halves — the planes exist (no
//     segfault) and they carry DIGITAL SILENCE — and only the first had a test: deleting the `fill_n`
//     leaves 444 checks green while the band is clocked on the previous chunk's split audio instead,
//     worth 6.9 dB out of digital silence on un-bypass.
static void theBypassedBandIsClockedOnSilence()
{
    group ("law 11c — a bypassed band's narrowing call carries SILENCE, not the last chunk's audio");
    const int B = 96;
    multiband::MultibandCompressor<3> A, Bm;
    if (! (A.prepare (kFs, B, 2) && Bm.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    dynamics::CompressorParams cp;
    cp.thresholdDb = -30.0; cp.ratio = 5.0; cp.attackMs = 2.0; cp.releaseMs = 2.0; cp.kneeDb = 0.0;
    for (int b = 0; b < 3; ++b) { A.setBandParams (b, cp); Bm.setBandParams (b, cp); }
    // Band 1 stays LIVE through the warm-up, so its planes are valid and the difference cannot be the
    // null-plane crash: this is the semantic half on its own.
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 6; ++k)
    {
        fillTone (l, k * B, 900.0, 0.9f); r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        fillTone (l, k * B, 900.0, 0.9f); r = l; float* jo[2] = { l.data(), r.data() }; run (Bm.process (jo, 2, B));
    }
    const double loud = A.bandGainReductionDb (1);
    ok (loud < -1.0, "precondition: the band is holding real gain reduction from the loud material");
    A.setBandBypass (1, true); Bm.setBandBypass (1, true);
    // A narrows (the falling edge fires and hands band 1 its silent block). B stays wide and feeds real
    // silence, which is what that block is supposed to BE. A fast release makes the two answers diverge
    // within one block: silence releases the band, the previous chunk's audio does not.
    fillTone (l, 0, 900.0, 0.9f);
    { float* io[1] = { l.data() }; run (A.process (io, 1, B)); }
    { std::vector<float> z1 ((std::size_t) B, 0.0f), z2 ((std::size_t) B, 0.0f); float* jo[2] = { z1.data(), z2.data() }; run (Bm.process (jo, 2, B)); }
    // WHAT THIS CAN AND CANNOT SAY, because the obvious oracle is not one. `Bm` is NOT a silence
    // reference: at a live width the band loop SKIPS a bypassed band, so `Bm`'s band 1 is frozen where
    // the loud material left it (measured, -18.204 against A's -6.697). That asymmetry is law 11c's one
    // inexact corner, documented in `MultibandProcessor.h` — it is the thing being compared against, not
    // a bug. So the assertion here is the direction: the edge RELEASED the band, which silence does and
    // the previous chunk's loud audio does not.
    // It is a weaker check than the rest of this file and that is stated rather than dressed up: a
    // mutation that deletes the `fill_n` and hands the band the previous chunk's split audio survives it
    // whenever that chunk's content in THIS band happens to be quiet. Closing it needs an observable for
    // "what the band was clocked on", which the composite does not expose today.
    ok (A.bandGainReductionDb (1) > loud + 0.5,
        "the edge RELEASED the bypassed band — silence does that, the previous chunk's audio does not ("
        + std::to_string (A.bandGainReductionDb (1)) + " dB from " + std::to_string (loud) + ")");
}

// (m) THE DE-ZIPPER'S TERM IN THE SHAPER'S EXIT PREDICATE — the exact sibling of the gate's HOLD term.
//     Give the gain smoother a longer time constant than either envelope and it is the last thing still
//     moving, so a predicate that ignores it exits early and freezes the gain mid-ramp.
static void theShapersDeZipperIsInThePredicate()
{
    group ("law 11c — TransientShaper: the pause does not exit while the de-zipper is still travelling");
    const int B = 128, gap = 20000;
    dynamics::TransientShaperParams p;
    // The envelopes park in a handful of samples; the de-zipper needs 200 ms. That ordering is the test.
    p.attackDb = 12.0; p.sustainDb = -9.0; p.threshold = 0.05;
    p.fastAttackMs = 0.05; p.fastReleaseMs = 0.05; p.slowAttackMs = 0.06; p.slowReleaseMs = 0.06;
    p.gainSmoothMs = 200.0;
    dynamics::TransientShaper A, Bs;
    if (! (A.prepare (kFs, B, 2) && Bs.prepare (kFs, B, 2))) { ok (false, "prepare"); return; }
    A.setParams (p); Bs.setParams (p);
    std::vector<float> l ((std::size_t) B), r ((std::size_t) B);
    for (int k = 0; k < 24; ++k)
    {
        for (int i = 0; i < B; ++i) { const int t = k * B + i; const float e = (t % 600 < 20) ? 1.0f : 0.02f;
                                      l[(std::size_t) i] = e * (float) std::sin (2.0 * core::kPi * 300.0 * t / kFs); }
        r = l; float* io[2] = { l.data(), r.data() }; run (A.process (io, 2, B));
        for (int i = 0; i < B; ++i) { const int t = k * B + i; const float e = (t % 600 < 20) ? 1.0f : 0.02f;
                                      l[(std::size_t) i] = e * (float) std::sin (2.0 * core::kPi * 300.0 * t / kFs); }
        r = l; float* jo[2] = { l.data(), r.data() }; run (Bs.process (jo, 2, B));
    }
    for (int off = 0; off < gap; off += B)
    {
        const int n = std::min (B, gap - off);
        float* io[2] = { nullptr, nullptr }; run (A.process (io, 0, n));
        std::vector<float> z1 ((std::size_t) n, 0.0f), z2 ((std::size_t) n, 0.0f);
        float* jo[2] = { z1.data(), z2.data() }; run (Bs.process (jo, 2, n));
    }
    std::vector<float> al ((std::size_t) B, 1.0f), ar ((std::size_t) B, 1.0f), bl ((std::size_t) B, 1.0f), br ((std::size_t) B, 1.0f);
    float* ai[2] = { al.data(), ar.data() }; float* bi[2] = { bl.data(), br.data() };
    run (A.process (ai, 2, B)); run (Bs.process (bi, 2, B));
    ok (bitsEqual (al[0], bl[0]), "the de-zipper arrived at the same place silence takes it");
    ok (! bitsEqual (bl[0], 1.0f), "precondition: the de-zipper was NOT at unity — it had somewhere to travel");
}

int main()
{
    std::printf ("law 11c — a pause is silence\n");
    floorPremise();
    collapseAgainstHonestLoop();
    fixedPointIsIdempotent();
    compressorInvariant();
    shaperInvariant();
    gateInvariant();
    dynamicEqBandInvariant();
    controlCounterPhase();
    deEsserInvariant();
    powerAmpInvariant();
    laneDynamicsInvariant();
    multibandDoesNotDoubleClock();
    multibandBypassedNarrowingDoesNotCrash();
    gateDropsItsSidechainHighPass();
    compressorKeepsItsExternalKey();
    controlCounterAtHugeLengths();
    coefficientWriteLandsInTheRemainder();
    thePauseFlushes();
    thePauseLeavesNoSubnormal();
    theGapFlushesWhereItCanBeSeen();
    aPoisonedStateSurvivesAPause();
    multibandDoubleClockAtZeroWidth();
    bypassedBandResidualIsOneChunk();
    counterWideningOnASettledBand();
    aRefusedBypassedChildIsReported();
    theFollowerBranchIsPerSample();
    theGatesExitPredicateNeedsAllOfIt();
    theTapAndTheKeyAtWidthZero();
    theBypassedBandIsClockedOnSilence();
    theShapersDeZipperIsInThePredicate();
    return felitronics::test::report();
}
