// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P14, part one — what the shared RECURSIVE PRIMITIVES do with a poisoned state.
//
// `state = in + c*(state - in)` never returns a NaN to a finite value, and `core::flushDenormal`
// steps straight over it: `fabs(NaN) < 1e-15f` is FALSE. So a denormal-only guard on a recursive
// kernel leaves exactly the state it would need to clear, and one bad sample is permanent. Measured
// on the real filters before this change — one +Inf into a band-pass followed by a healthy 7 kHz
// tone, with the per-block flush running the whole time: 479000 of the next 480000 outputs
// non-finite. NaN is identical. A finite 1e20 is NOT: the filters carry it without poisoning at all,
// which matters because it means the case for hardening them rests on Inf/NaN and nothing else.
//
// WHAT THIS FILE DOES AND DOES NOT ASSERT. Clearing poisoned state is a BOUND on the damage, not an
// undo. The flush is clocked by the owner, so a poisoned filter emits rubbish to the end of the
// caller's block and then restarts from silence — and restarting is precisely `reset()`, which is
// asserted here exactly. It follows that the recovered stream is NOT the stream a sanitised input
// would have produced, and that the damage window DEPENDS on the caller's block size. Both are
// pinned below, the second one deliberately, so that no later change can quietly claim a
// partition-independence this mechanism does not have.
//
// The audio path is where this is the ONLY remedy: the house rule is that it is not sanitised (a
// compressor is not a sanitiser — bounding the signal belongs to the limiter), so nothing may gate an
// audio filter's input. A DETECTOR is the opposite case and must be gated upstream instead; that is
// part two of P14 and not this file.
//
// `core::flushPoison`'s own contract — totality, idempotence, and agreement with `flushDenormal` on
// every finite input — is swept exhaustively by class in `core/tests/FlushPoisonTheoryTests.cpp`.
// What is left for here is the APPLICATION of it to four states, and the one place where the obvious
// application would have been wrong.

#include <felitronics_test.h>

#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

using namespace felitronics;

static constexpr double kFs = 48000.0;
static const float kInf = std::numeric_limits<float>::infinity();
static const float kNan = std::numeric_limits<float>::quiet_NaN();

static float toneAt (int i, double f = 7000.0, double amp = 0.5)
{
    return (float) (amp * std::sin (2.0 * core::kPi * f * i / kFs));
}

// The OLD guard, written out so the healthy-stream claim is checked against something rather than
// asserted. Only `Biquad` can be checked this way in-suite — its state is public; `Svf`'s is not, and
// for it the same claim is carried by the out-of-tree regression against origin/main (270
// configurations, 2.23 M samples, byte-identical), which is said here rather than left implied.
static void oldFlush (eq::Biquad& f) noexcept
{
    if (std::fabs (f.z1) < 1e-15f) f.z1 = 0.0f;
    if (std::fabs (f.z2) < 1e-15f) f.z2 = 0.0f;
}

static eq::BiquadCoeffs lowpass (double fc, double q)
{
    const double w = 2.0 * core::kPi * fc / kFs, a = std::sin (w) / (2.0 * q), cw = std::cos (w);
    const double a0 = 1.0 + a;
    eq::BiquadCoeffs c;
    c.b0 = ((1.0 - cw) * 0.5) / a0; c.b1 = (1.0 - cw) / a0; c.b2 = c.b0;
    c.a1 = (-2.0 * cw) / a0;        c.a2 = (1.0 - a) / a0;
    return c;
}

//==============================================================================
static void poisonIsPermanentWithoutTheFlush()
{
    test::group ("the defect: an owner that never flushes never recovers");

    // The flush is the owner's call, so "not calling it" is the honest way to show what the state
    // does on its own — and it is also the shape of every consumer that forgets to.
    for (int kind = 0; kind < 2; ++kind)
    {
        const float bad = (kind == 0) ? kInf : kNan;
        const char* nm = (kind == 0) ? "+Inf" : "NaN";
        {
            eq::Svf f; f.prepare (kFs, 2); f.setParams (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
            int nonFinite = 0;
            for (int i = 0; i < 4000; ++i)
            {
                const float y = f.processSample (0, (i == 100) ? bad : toneAt (i));
                if (i > 100 && ! std::isfinite (y)) ++nonFinite;
            }
            test::ok (nonFinite == 3899, std::string ("Svf: one ") + nm
                      + " with no flush leaves every later output non-finite (" + std::to_string (nonFinite) + "/3899)");
        }
        {
            eq::Biquad f; f.setCoeffs (lowpass (2000.0, 0.7));
            int nonFinite = 0;
            for (int i = 0; i < 4000; ++i)
            {
                const float y = f.processSample ((i == 100) ? bad : toneAt (i));
                if (i > 100 && ! std::isfinite (y)) ++nonFinite;
            }
            test::ok (nonFinite == 3899, std::string ("Biquad: same, ") + nm
                      + " (" + std::to_string (nonFinite) + "/3899)");
        }
    }

    // A FINITE 1e20 does not poison either filter, and pinning that matters: it is the reason the
    // case for hardening these primitives rests on Inf/NaN, and the reason a detector's 1e20 problem
    // is a problem of the SQUARE downstream rather than of the filter.
    {
        eq::Svf f; f.prepare (kFs, 2); f.setParams (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
        eq::Biquad b; b.setCoeffs (lowpass (2000.0, 0.7));
        int bad = 0;
        for (int i = 0; i < 4000; ++i)
        {
            const float x = (i == 100) ? 1.0e20f : toneAt (i);
            if (! std::isfinite (f.processSample (0, x))) ++bad;
            if (! std::isfinite (b.processSample (x))) ++bad;
        }
        test::ok (bad == 0, "a finite 1e20 poisons neither filter — not one non-finite output in 8000");
    }
}

//==============================================================================
static void theFlushRecoversAndIsExactlyAReset()
{
    test::group ("recovery: the flush clears the poison, and what it leaves is reset()");

    for (int kind = 0; kind < 2; ++kind)
    {
        const float bad = (kind == 0) ? kInf : kNan;
        const char* nm = (kind == 0) ? "+Inf" : "NaN";

        // Svf: charge, poison, flush at the block edge, then a healthy block must be finite AND alive.
        {
            eq::Svf f; f.prepare (kFs, 2); f.setParams (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
            for (int i = 0; i < 256; ++i) f.processSample (0, toneAt (i));
            f.processSample (0, bad);
            f.flushDenormals();
            double energy = 0.0; int nonFinite = 0;
            for (int i = 0; i < 2000; ++i)
            {
                const float y = f.processSample (0, toneAt (i));
                if (! std::isfinite (y)) ++nonFinite; else energy += std::fabs ((double) y);
            }
            test::ok (nonFinite == 0, std::string ("Svf recovers from ") + nm + ": nothing non-finite after the flush");
            test::ok (energy > 1.0, "and it is alive, not stuck at zero (a gate on the output would give zero)");
        }
        // ...and the recovered filter is bit-for-bit a reset one. That is the whole semantics: the
        // history is gone, so this is a restart and must not be described as a repair.
        {
            eq::Svf healed, fresh;
            healed.prepare (kFs, 2); fresh.prepare (kFs, 2);
            healed.setParams (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
            fresh.setParams  (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
            for (int i = 0; i < 256; ++i) healed.processSample (0, toneAt (i));
            healed.processSample (0, bad);
            healed.flushDenormals();
            fresh.reset();
            int diff = 0;
            for (int i = 0; i < 2000; ++i)
            {
                const float a = healed.processSample (0, toneAt (i));
                const float b = fresh.processSample  (0, toneAt (i));
                if (std::memcmp (&a, &b, 4) != 0) ++diff;
            }
            test::ok (diff == 0, std::string ("Svf after ") + nm + " + flush is bit-for-bit a reset() filter");
        }
        {
            eq::Biquad healed, fresh;
            healed.setCoeffs (lowpass (2000.0, 0.7)); fresh.setCoeffs (lowpass (2000.0, 0.7));
            for (int i = 0; i < 256; ++i) healed.processSample (toneAt (i));
            healed.processSample (bad);
            healed.flushDenormals();
            int diff = 0;
            for (int i = 0; i < 2000; ++i)
            {
                const float a = healed.processSample (toneAt (i));
                const float b = fresh.processSample  (toneAt (i));
                if (std::memcmp (&a, &b, 4) != 0) ++diff;
            }
            test::ok (diff == 0, std::string ("Biquad after ") + nm + " + flush is bit-for-bit a reset() filter");
        }
    }

    // One channel's poison must not touch another's history — the state is per channel and the flush
    // must not become a global reset by accident.
    {
        eq::Svf a, b;
        a.prepare (kFs, 2); b.prepare (kFs, 2);
        a.setParams (eq::FilterType::LowPass, 3000.0, 0.707, 0.0);
        b.setParams (eq::FilterType::LowPass, 3000.0, 0.707, 0.0);
        int diff = 0;
        for (int i = 0; i < 3000; ++i)
        {
            a.processSample (0, (i == 500) ? kInf : toneAt (i));       // poisoned channel
            b.processSample (0, toneAt (i));
            const float ya = a.processSample (1, toneAt (i, 300.0));   // the innocent neighbour
            const float yb = b.processSample (1, toneAt (i, 300.0));
            if (std::memcmp (&ya, &yb, 4) != 0) ++diff;
            if ((i % 128) == 127) { a.flushDenormals(); b.flushDenormals(); }
        }
        test::ok (diff == 0, "poisoning channel 0 leaves channel 1 bit-for-bit untouched, flush included");
    }
}

//==============================================================================
static void theDamageWindowIsTheCallersBlockAndSaysSo()
{
    test::group ("the bound is the OWNER's block — asserted, so it cannot later be called partition-independent");

    int windows[4] = {}; int k = 0;
    for (int blk : { 32, 128, 512, 4096 })
    {
        eq::Svf f; f.prepare (kFs, 2); f.setParams (eq::FilterType::BandPass, 7000.0, 2.0, 0.0);
        int nonFinite = 0;
        const int at = 8;                                   // early in the first block
        for (int i = 0; i < 3 * 4096; ++i)
        {
            const float y = f.processSample (0, (i == at) ? kInf : toneAt (i));
            if (i > at && ! std::isfinite (y)) ++nonFinite;
            if ((i % blk) == blk - 1) f.flushDenormals();
        }
        windows[k++] = nonFinite;
        test::ok (nonFinite <= blk, "block " + std::to_string (blk) + ": the damage is at most one block ("
                  + std::to_string (nonFinite) + " samples)");
        test::ok (nonFinite > 0, "and it is not zero — the flush bounds the damage, it does not prevent it");
    }
    test::ok (windows[0] < windows[3],
              "the window GROWS with the caller's block (" + std::to_string (windows[0]) + " against "
              + std::to_string (windows[3]) + ") — this mechanism is partition-DEPENDENT by construction");
}

//==============================================================================
static void healthyStreamsAreUntouched()
{
    test::group ("healthy streams: the new guard is a strict superset of the old one");

    // Checked against the old rule itself, on the streams whose state actually visits the region the
    // threshold cares about — a decaying tail, exact zeros of both signs, a level below the threshold,
    // a subnormal, and the threshold itself.
    struct Stream { const char* name; float (*gen) (int); };
    static const Stream streams[] = {
        { "decaying tail",   [] (int i) { return (float) (0.7 * std::exp (-0.004 * i) * std::sin (0.11 * i)); } },
        { "exact -0.0",      [] (int i) { return (i & 1) ? -0.0f : 0.0f; } },
        { "1e-16 tone",      [] (int i) { return (float) (1.0e-16 * std::sin (0.07 * i)); } },
        { "subnormal 1e-40", [] (int i) { return (i % 3 == 0) ? 1.0e-40f : 0.0f; } },
        { "alternating 1e-15", [] (int i) { return (i & 1) ? 1.0e-15f : -1.0e-15f; } },
        { "ordinary tone",   [] (int i) { return toneAt (i); } },
    };
    for (const auto& s : streams)
     for (double fc : { 80.0, 2000.0, 15000.0 })
      for (double q : { 0.5, 0.707, 8.0 })
      {
          eq::Biquad now, then;
          now.setCoeffs (lowpass (fc, q)); then.setCoeffs (lowpass (fc, q));
          int diff = 0;
          for (int i = 0; i < 6000; ++i)
          {
              const float x = s.gen (i);
              const float a = now.processSample (x);
              const float b = then.processSample (x);
              if (std::memcmp (&a, &b, 4) != 0) ++diff;
              if ((i % 64) == 63) { now.flushDenormals(); oldFlush (then); }
          }
          test::ok (diff == 0, std::string ("Biquad, ") + s.name + " at " + std::to_string ((int) fc)
                    + " Hz: bit-identical to the denormal-only guard");
      }
}

//==============================================================================
static void theFollowersKeepTheirOwnThresholds()
{
    test::group ("EnvelopeFollower — poison cleared WITHOUT undoing the power-domain threshold");

    // THE TRAP. `flushPoison` carries the house 1e-15, and in Rms mode the state is a POWER: 1e-15 of
    // power is an amplitude of 3.2e-8, i.e. -150 dBFS. Swapping it in wholesale would have silently
    // restored the exact defect P2 measured and fixed — -7.5 dB of gain reduction appearing or not
    // appearing depending on how the caller cut the stream. The non-finite half is what had to be
    // added; the 1e-30 threshold is what had to survive, and this is the check that says so.
    {
        const float quiet = (float) core::dbToGain (-160.0);          // well under the old 1e-15 of power
        dynamics::EnvelopeFollower f;
        f.prepare (kFs); f.setDetector (dynamics::Detector::Rms); f.setTimes (5.0, 5.0);
        float env = 0.0f;
        for (int blk = 0; blk < 40; ++blk)
        {
            for (int i = 0; i < 256; ++i) env = f.process (quiet);
            f.flushDenormals();
        }
        test::ok (env > 0.0f, "a -160 dBFS input still produces a NON-ZERO Rms envelope after 40 flushes");
        test::approx ((double) env, (double) quiet, (double) quiet * 0.05,
                      "and it tracks the input, rather than being zapped to silence by a 1e-15 power floor");
    }
    // ...and the poison half does work, in both detector modes.
    for (auto det : { dynamics::Detector::Peak, dynamics::Detector::Rms })
     for (float bad : { kInf, kNan })
     {
         dynamics::EnvelopeFollower f;
         f.prepare (kFs); f.setDetector (det); f.setTimes (5.0, 50.0);
         for (int i = 0; i < 100; ++i) f.process (0.4f);
         f.process (bad);
         test::ok (! std::isfinite (f.envelope()), "the poison does reach the state (otherwise this proves nothing)");
         f.flushDenormals();
         float last = 0.0f;
         for (int i = 0; i < 2000; ++i) last = f.process (0.4f);
         test::ok (std::isfinite (last) && last > 0.1f,
                   std::string (det == dynamics::Detector::Peak ? "Peak" : "Rms")
                   + ": the follower recovers and tracks again after the flush");
     }

    test::group ("GainReductionFollower — the same, for the ballistics on the gain");
    for (float bad : { kInf, kNan })
    {
        dynamics::GainReductionFollower g;
        g.prepare (kFs); g.setTimes (5.0, 50.0);
        for (int i = 0; i < 100; ++i) g.process (-6.0f);
        g.process (bad);
        test::ok (! std::isfinite (g.valueDb()), "poisoned");
        g.flushDenormals();
        float last = 0.0f;
        for (int i = 0; i < 2000; ++i) last = g.process (-6.0f);
        test::approx ((double) last, -6.0, 0.01, "recovers and converges on the target again");
    }
    // The denormal half must still behave exactly as before: a release to zero ends at exact zero.
    {
        dynamics::GainReductionFollower g;
        g.prepare (kFs); g.setTimes (1.0, 1.0);
        for (int i = 0; i < 200; ++i) g.process (-12.0f);
        for (int i = 0; i < 4000; ++i) g.process (0.0f);
        g.flushDenormals();
        test::ok (g.valueDb() == 0.0f, "a decayed follower still flushes to EXACT zero, as the denormal guard did");
    }
}

//==============================================================================
int main()
{
    std::printf ("P14 part one — poison recovery in the shared recursive primitives\n");
    poisonIsPermanentWithoutTheFlush();
    theFlushRecoversAndIsExactlyAReset();
    theDamageWindowIsTheCallersBlockAndSaysSo();
    healthyStreamsAreUntouched();
    theFollowersKeepTheirOwnThresholds();
    return felitronics::test::report();
}
