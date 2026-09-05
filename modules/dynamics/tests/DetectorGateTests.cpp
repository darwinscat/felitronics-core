// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P14, part two — PREVENTION. The four modules that fed a recursive sidechain from a raw signal now
// gate it, per sample, BEFORE the filter.
//
// WHY BEFORE AND NOT AFTER, since the difference is not obvious. The sidechain filter is recursive: a
// poisoned sample reaching it is permanent, and a gate placed on its OUTPUT reads the resulting NaN as
// silence — which turns "noise for ever" into "silence for ever" and fixes nothing. Measured on a
// 7 kHz band-pass at Q 2 feeding an Rms follower: gate-after diverges from an explicitly sanitised
// reference by 1222 to 5533 samples and 0.2 to 171 dB, the amount depending on the caller's block
// size; gate-before diverges by ZERO, at every block size. Part one's per-block recovery cannot
// substitute either — it restarts the filter from silence, which is a `reset()`, not a repair.
//
// WHAT AN INPUT GATE IS ENOUGH FOR, measured rather than assumed. The gate bounds the filter input to
// ±1e6, and the largest output any configuration produces from that — swept over every filter type ×
// 5 frequencies × Q up to 1000 × ±60 dB of gain, driven by a square wave at the filter's own centre
// frequency — is 1.98e9. Its square is 3.9e18 against a float maximum of 3.4e38, so the squaring that
// `MeanPower` and the Rms follower do cannot overflow, and a second gate on the filter's output would
// buy nothing. The output would have to reach 1.8e19 first.
//
// THE ACCEPTANCE IS THE GAIN TRAJECTORY, NOT THE AUDIO, and that is forced rather than chosen. These
// modules are SELF-KEYED: the poisoned sample is a sample of the PROGRAMME, which is also the key, and
// the audio path is deliberately not sanitised (a compressor is not a sanitiser). So at the poisoned
// sample the two outputs must differ, and a bit-exact null of the audio is not merely hard, it is
// wrong to ask for. What must be identical is what the DETECTOR decided: the gain reduction, the band
// delta, or an innocent channel's output, which carries the shared gain and nothing else.
// `LaneDynamics` is the exception — it takes its sidechain as a separate buffer, so poisoning only
// that leaves the audio untouched and the audio null is honest there.

#include <felitronics_test.h>

#include <felitronics/dynamics/TransientShaper.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/deesser/DeEsser.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;

static constexpr double kFs = 48000.0;
static const float kInf = std::numeric_limits<float>::infinity();
static const float kNan = std::numeric_limits<float>::quiet_NaN();

// The poison values and what an explicitly sanitised key would carry instead. The last two are the
// ones that matter for the CLAMP rather than the finiteness check: 1e20 squares to +Inf in float, and
// 2e6 is simply past the ceiling. A fixture without them tests half the gate.
struct Poison { const char* name; float bad; float sane; };
static const Poison kPoisons[] = {
    { "+Inf", kInf,     0.0f     },
    { "-Inf", -kInf,    0.0f     },
    { "NaN",  kNan,     0.0f     },
    { "1e20", 1.0e20f,  1.0e6f   },
    { "2e6",  2.0e6f,   1.0e6f   },
    { "-2e6", -2.0e6f, -1.0e6f   },
};

static float toneAt (int i, double f = 900.0, double amp = 0.35)
{
    return (float) (amp * std::sin (2.0 * core::kPi * f * i / kFs));
}

//==============================================================================
// TWO THINGS THIS FILE GOT WRONG FIRST, both worth keeping, because both produced convincing evidence
// for a defect that did not exist.
//
// ONE: THE BUFFER LENGTH MUST BE A MULTIPLE OF THE BLOCK SIZE. It was not — 24000 samples in blocks of
// 256 — so the last block ran 64 samples past the end of every channel. What that produced was not a
// crash but a STORY: two runs of `TransientShaper` over identical input appeared to differ, by up to
// 0.0106 on a 0.35 signal, on the wasm tier but not the desktop, at `-O0` as well as `-O2`, and on
// `origin/main` as well as here. Every one of those observations was real, and together they made a
// tidy case for a tier-specific non-determinism in the module. There is none: with the length fixed,
// all four modules repeat exactly, on both tiers. The reading past the end was picking up whatever
// happened to be adjacent, which differed between runs and between platforms. `AddressSanitizer`
// named it in one line after the measurements had been believed for an hour.
//
// TWO: TWO RUNS COMPARED BIT-FOR-BIT SHOULD SHARE THEIR STORAGE. Once the overrun was gone this
// stopped being load-bearing — the modules are deterministic and separate buffers agree — but it is
// kept because it removes the buffer address from the set of things a failure could mean, and that
// set was the expensive part.
struct SharedRig
{
    explicit SharedRig (int lanes, int samples) : n_ (samples), buf_ ((std::size_t) lanes, std::vector<float> ((std::size_t) samples)) {}
    float* lane (int i) noexcept { return buf_[(std::size_t) i].data(); }
    void fill (int i, const std::vector<float>& src) noexcept { std::memcpy (lane (i), src.data(), (std::size_t) n_ * sizeof (float)); }
    std::vector<float> snapshot (int i) const { return buf_[(std::size_t) i]; }
private:
    int n_;
    std::vector<std::vector<float>> buf_;
};

//==============================================================================
// TransientShaper — no gain accessor, so the innocent CHANNEL is the trajectory: its own input is
// untouched, so its output is bit-identical exactly when the shared gain is.
static void transientShaperSubstitutesItsKey()
{
    test::group ("TransientShaper — a poisoned channel does not move the gain the others get");

    const int n = 24064, blk = 256, at = 700;   // n IS A MULTIPLE OF blk — see the note at the top
    for (const auto& p : kPoisons)
     for (auto link : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
     {
         std::vector<float> src0 ((std::size_t) n), src1 ((std::size_t) n);
         for (int i = 0; i < n; ++i) { src0[(std::size_t) i] = toneAt (i, 900.0); src1[(std::size_t) i] = toneAt (i, 310.0); }

         SharedRig rig (2, n);
         auto run = [&] (float poisonValue)
         {
             rig.fill (0, src0); rig.fill (1, src1);
             rig.lane (0)[at] = poisonValue;
             dynamics::TransientShaper t; t.prepare (kFs, blk, 2);
             dynamics::TransientShaperParams tp;
             tp.attackDb = 8.0; tp.sustainDb = -4.0; tp.mix = 1.0; tp.link = link;
             t.setParams (tp);
             for (int off = 0; off < n; off += blk)
             {
                 float* io[2] { rig.lane (0) + off, rig.lane (1) + off };
                 t.process (io, 2, blk);
             }
             return rig.snapshot (1);                 // the innocent channel
         };
         // The clean case is run TWICE first, and the two are required to agree exactly. That is not
         // ceremony: an apparent failure of the check below is otherwise ambiguous between "the gate
         // is wrong" and "the module does not repeat", and this file spent an hour on the wrong branch
         // of exactly that ambiguity. Asserting determinism first makes the next assertion mean one
         // thing.
         const auto sanitised = run (p.sane), sanitisedAgain = run (p.sane), poisoned = run (p.bad);
         const std::string what = std::string ("TransientShaper ")
                                + (link == dynamics::LinkMode::Max ? "Max" : "MeanPower") + ", " + p.name;
         test::ok (std::memcmp (sanitised.data(), sanitisedAgain.data(), (std::size_t) n * sizeof (float)) == 0,
                   what + ": the module repeats exactly, so the next check means what it says");
         test::ok (std::memcmp (sanitised.data(), poisoned.data(), (std::size_t) n * sizeof (float)) == 0,
                   what + ": the innocent channel is BIT-IDENTICAL to the sanitised run");
     }
}

//==============================================================================
// DeEsser and DynamicEqBand expose the decision directly. Block size 1 so the reported value is a
// per-sample trajectory rather than an end-of-block sample of one.
static void theGainTrajectoriesAreSubstituted()
{
    test::group ("DeEsser / DynamicEqBand — the decision is bit-identical to the sanitised key");

    // THE BLOCK SIZE IS PART OF THE FIXTURE. At block 1 the per-block recovery from part one fires
    // after every sample, so a filter gated on its OUTPUT recovers almost instantly and the wrong
    // design passes — measured, that mutation survived until these larger sizes were added. Prevention
    // has to be shown where recovery is weak, which is the block size a host actually uses.
    const int n = 6000, at = 400;
    for (const auto& p : kPoisons)
    {
        for (auto mode : { deesser::DeEsserMode::DynamicEq, deesser::DeEsserMode::SplitBand })
         for (int blk : { 1, 64, 512 })
        {
            std::vector<double> traceA, traceB;
            for (int which = 0; which < 2; ++which)
            {
                // THE LEVEL IS PART OF THE FIXTURE TOO. At 0.4 the de-esser sits pinned at its full
                // 8 dB range for the entire run, and a saturated detector reports the same number
                // whatever it is fed — the wrong design passed until this dropped to a level where
                // the gain reduction is mid-range and moving. Bursts for the same reason.
                std::vector<float> x ((std::size_t) n);
                for (int i = 0; i < n; ++i)
                {
                    const double burst = ((i / 700) % 3 == 0) ? 2.5 : 1.0;
                    x[(std::size_t) i] = (float) (burst * (double) toneAt (i, 6500.0, 0.045));
                }
                x[(std::size_t) at] = which == 0 ? p.bad : p.sane;
                deesser::DeEsser d; d.prepare (kFs, blk, 1);
                deesser::DeEsserParams dp; dp.mode = mode; d.setParams (dp);
                auto& trace = which == 0 ? traceA : traceB;
                for (int off = 0; off + blk <= n; off += blk)
                { float* io[1] { x.data() + off }; d.process (io, 1, blk); trace.push_back (d.gainReductionDb()); }
            }
            int bad = 0;
            for (std::size_t i = 0; i < traceA.size(); ++i)
                if (std::memcmp (&traceA[i], &traceB[i], sizeof (double)) != 0) ++bad;
            test::ok (bad == 0, std::string ("DeEsser ")
                      + (mode == deesser::DeEsserMode::DynamicEq ? "DynamicEq" : "SplitBand")
                      + " blk " + std::to_string (blk) + ", " + p.name
                      + ": gain-reduction trajectory bit-identical (" + std::to_string (bad) + " of "
                      + std::to_string (traceA.size()) + " differ)");
        }

        for (auto link : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
         for (int blk : { 1, 64, 512 })
        {
            std::vector<double> traceA, traceB;
            for (int which = 0; which < 2; ++which)
            {
                std::vector<float> l ((std::size_t) n), r ((std::size_t) n);
                for (int i = 0; i < n; ++i) { l[(std::size_t) i] = toneAt (i, 400.0, 0.5); r[(std::size_t) i] = toneAt (i, 1500.0, 0.3); }
                l[(std::size_t) at] = which == 0 ? p.bad : p.sane;
                dynamiceq::DynamicEqBand e; e.prepare (kFs, 2);
                dynamiceq::DynamicEqBandParams ep; ep.link = link; ep.coeffUpdatePeriod = 1;
                e.setParams (ep);
                auto& trace = which == 0 ? traceA : traceB;
                for (int off = 0; off + blk <= n; off += blk)
                { float* io[2] { l.data() + off, r.data() + off }; e.process (io, 2, blk); trace.push_back (e.dynamicDeltaDb()); }
            }
            int bad = 0;
            for (std::size_t i = 0; i < traceA.size(); ++i)
                if (std::memcmp (&traceA[i], &traceB[i], sizeof (double)) != 0) ++bad;
            test::ok (bad == 0, std::string ("DynamicEqBand ")
                      + (link == dynamics::LinkMode::Max ? "Max" : "MeanPower")
                      + " blk " + std::to_string (blk) + ", " + p.name
                      + ": band-delta trajectory bit-identical (" + std::to_string (bad) + " differ)");
        }
    }
}

//==============================================================================
// LaneDynamics takes its sidechain separately, so here the AUDIO null is the honest one — and the
// mid/side lanes are the reason the gate sits before `laneSignal`'s arithmetic rather than after it:
// `0.5f * (Inf - Inf)` is a NaN indistinguishable from a legitimate zero.
static void laneDynamicsGatesBeforeTheMidSideArithmetic()
{
    test::group ("LaneDynamics — gated before the M/S subtraction, audio nulled against the sanitised key");

    const int n = 4096, blk = 128, at = 300;
    for (const auto& p : kPoisons)
     for (int poisonRight = 0; poisonRight < 2; ++poisonRight)
     {
         std::vector<float> src[2], sc[2];
         for (int c = 0; c < 2; ++c)
         {
             src[c].resize ((std::size_t) n); sc[c].resize ((std::size_t) n);
             for (int i = 0; i < n; ++i)
             {
                 // Bursts, not a steady tone: a stationary signal BECOMES the programme and the band
                 // correctly idles, which would make every configuration look alike for the wrong
                 // reason. The module's own suite says so and it is right.
                 const double burst = ((i / 700) % 3 == 0) ? 3.0 : 1.0;
                 const float v = (float) (burst * (double) toneAt (i, c == 0 ? 800.0 : 2100.0, 0.25));
                 src[c][(std::size_t) i] = v; sc[c][(std::size_t) i] = v;
             }
         }
         const int ch = poisonRight ? 1 : 0;

         // Shared storage, for the reason at the top of this file.
         SharedRig rig (4, n);
         auto run = [&] (float poisonValue)
         {
             for (int c = 0; c < 2; ++c) { rig.fill (c, src[c]); rig.fill (2 + c, sc[c]); }
             rig.lane (2 + ch)[at] = poisonValue;
             dynamiceq::LaneDynamics ld; ld.prepare (kFs, 2);
             // THE BAND ITSELF MUST BE ON, not only its lanes and its dynamics. `BandParams::on`
             // defaults to false and `processBand` then returns before touching the sidechain — the
             // first version of this test exercised nothing at all, and every mutation of the gate
             // survived it. Recipe taken from the module's own suite rather than guessed.
             eq::BandParams bp;
             bp.on = true; bp.type = eq::FilterType::Bell;
             for (auto lane : { eq::Lane::Stereo, eq::Lane::Mid, eq::Lane::Side })
             { auto& lp = bp.lane (lane); lp.on = true; lp.freq = 800.0; lp.Q = 1.0; lp.gainDb = 0.0; }
             bp.dyn.on = true; bp.dyn.rangeDb = -12.0;
             ld.setParams (bp);
             eq::EqBand band; band.prepare (kFs, 2); band.setParams (bp);
             for (int off = 0; off < n; off += blk)
             {
                 float* a[2] { rig.lane (0) + off, rig.lane (1) + off };
                 const float* s[2] { rig.lane (2) + off, rig.lane (3) + off };
                 ld.processBand (a, s, 2, blk, band);
             }
             std::vector<float> both;
             for (int c = 0; c < 2; ++c) { const auto v = rig.snapshot (c); both.insert (both.end(), v.begin(), v.end()); }
             return both;
         };
         const auto poisoned = run (p.bad), sanitised = run (p.sane);
         const int bad = std::memcmp (poisoned.data(), sanitised.data(), poisoned.size() * sizeof (float)) == 0 ? 0 : 1;
         test::ok (bad == 0, std::string ("LaneDynamics, ") + p.name + " in channel "
                   + std::to_string (ch) + ": the audio is bit-identical to the sanitised run");
     }
}

//==============================================================================
// ISOLATION, and where it actually lives. The first version of this checked that poisoning one
// channel left another channel's AUDIO untouched — and it failed, correctly: these modules apply ONE
// linked gain to every channel, so one channel's key moves them all by construction. That is what
// linking is for; a compressor that gained channels independently would shift the stereo image on
// every transient. So the invariant at THIS level is substitution, which the groups above check, and
// isolation belongs to the primitive — where part one already tests it directly on `eq::Svf`
// (poison channel 0, channel 1's filter output bit-identical).
//
// What is worth pinning here instead is that the gain really is SHARED, because "fixing" the failure
// above by making it per-channel would be a silent regression of exactly the property the linking
// exists to provide.
static void theGainIsSharedAndThatIsWhyIsolationIsAPrimitiveProperty()
{
    test::group ("the linked gain is shared — the reason a channel cannot be isolated at this level");

    const int n = 4000;
    std::vector<float> l ((std::size_t) n), r ((std::size_t) n), lRef ((std::size_t) n), rRef ((std::size_t) n);
    for (int i = 0; i < n; ++i)
    {
        // Channel 0 carries a loud burst the band will duck; channel 1 carries a steady quiet tone.
        const double burst = (i > n / 2 && i < n / 2 + 500) ? 0.9 : 0.02;
        l[(std::size_t) i] = lRef[(std::size_t) i] = (float) (burst * std::sin (2.0 * core::kPi * 1000.0 * i / kFs));
        r[(std::size_t) i] = rRef[(std::size_t) i] = toneAt (i, 1000.0, 0.2);
    }
    dynamiceq::DynamicEqBand e; e.prepare (kFs, 2);
    dynamiceq::DynamicEqBandParams ep; ep.link = dynamics::LinkMode::Max; ep.coeffUpdatePeriod = 1;
    ep.thresholdDb = -20.0; ep.ratio = 4.0; ep.rangeDb = 18.0;
    e.setParams (ep);
    for (int i = 0; i < n; ++i) { float* io[2] { l.data() + i, r.data() + i }; e.process (io, 2, 1); }

    double worst = 0.0;
    for (int i = n / 2 + 100; i < n / 2 + 400; ++i)
        worst = std::max (worst, std::fabs ((double) r[(std::size_t) i] - (double) rRef[(std::size_t) i]));
    test::ok (worst > 1.0e-4,
              "a burst in channel 0 audibly moves channel 1 (" + std::to_string (worst)
              + ") — the gain is LINKED, so per-channel isolation is not a property this level can have");
}

//==============================================================================
int main()
{
    std::printf ("P14 part two — the detector inputs are gated, before the filter\n");
    transientShaperSubstitutesItsKey();
    theGainTrajectoriesAreSubstituted();
    laneDynamicsGatesBeforeTheMidSideArithmetic();
    theGainIsSharedAndThatIsWhyIsolationIsAPrimitiveProperty();
    return felitronics::test::report();
}
