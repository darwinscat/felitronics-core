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
#include <felitronics/dynamics/TransientShaper.h>
#include <felitronics/dynamics/AutoLeveler.h>
#include <felitronics/deesser/DeEsser.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
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

// THE POISON SET, and -Inf is in it deliberately. A guard written as "NaN or +Inf" reads as complete
// and is not; a suite that feeds only those two blesses it. `std::isfinite` covers all three, so the
// cost of checking is nothing and the cost of NOT checking is a whole class of input.
struct Poison { const char* name; float value; };
static const Poison kPoisons[] = { { "+Inf", kInf }, { "-Inf", -kInf }, { "NaN", kNan } };

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
    for (const auto& p : kPoisons)
    {
        const float bad = p.value; const char* nm = p.name;
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

    for (const auto& p : kPoisons)
    {
        const float bad = p.value; const char* nm = p.name;

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
            a.processSample (0, (i == 500) ? -kInf : toneAt (i));      // poisoned channel (-Inf: the
                                                                      // value a "NaN or +Inf" guard
                                                                      // would step straight over)
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
    // The band between the two thresholds was unpinned: a mutation moving 1e-30 to 1e-20 survived the
    // whole tree, because nothing tracked a signal quiet enough to tell them apart. -230 dBFS is a
    // power of 1e-23, inside that band and still above `core::gainToDb`'s own floor.
    {
        const float veryQuiet = (float) core::dbToGain (-230.0);
        dynamics::EnvelopeFollower f;
        f.prepare (kFs); f.setDetector (dynamics::Detector::Rms); f.setTimes (5.0, 5.0);
        float env = 0.0f;
        for (int blk = 0; blk < 40; ++blk) { for (int i = 0; i < 256; ++i) env = f.process (veryQuiet); f.flushDenormals(); }
        test::approx ((double) env, (double) veryQuiet, (double) veryQuiet * 0.05,
                      "a -230 dBFS input is still tracked after 40 flushes — the threshold is 1e-30 of "
                      "POWER, and anything coarser silently swallows this band");
    }

    // ...and the poison half does work, in both detector modes.
    for (auto det : { dynamics::Detector::Peak, dynamics::Detector::Rms })
     for (const auto& p : kPoisons)
     {
         const float bad = p.value;
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
    for (const auto& p : kPoisons)
    {
        const float bad = p.value;
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
// PARTIAL poisoning: the case a direct Inf or NaN cannot reach, because it takes both state words at
// once. A finite input that OVERFLOWS on the way leaves one word non-finite and the other a large
// finite number, and a per-word flush then "heals" the filter into something that is not a reset
// filter at all — it carries the survivor forward and re-poisons itself from it. Measured on the
// biquad before the fix: two samples of SILENCE after the flush it emitted -3.27e38, then -Inf.
static void partialPoisoningIsHealedAtomically()
{
    test::group ("partial poisoning — the two state words are ONE state");

    // Reachable, and shown to be reachable rather than assumed: a finite 1e37 sine at Q 40.
    {
        eq::Biquad f; f.setCoeffs (lowpass (2000.0, 40.0));
        int partialAt = -1;
        for (int i = 0; i < 4000 && partialAt < 0; ++i)
        {
            f.processSample ((float) (1.0e37 * std::sin (2.0 * core::kPi * 2000.0 * i / kFs)));
            if (std::isfinite (f.z1) != std::isfinite (f.z2)) partialAt = i;
        }
        test::ok (partialAt > 0, "a FINITE input does leave one word non-finite and the other finite (sample "
                  + std::to_string (partialAt) + ")");

        f.flushDenormals();
        eq::Biquad fresh; fresh.setCoeffs (lowpass (2000.0, 40.0));
        int diff = 0;
        for (int k = 0; k < 64; ++k)
        {
            const float a = f.processSample (0.0f), b = fresh.processSample (0.0f);
            if (std::memcmp (&a, &b, 4) != 0) ++diff;
        }
        test::ok (diff == 0, "and after the flush it is bit-for-bit a reset filter over the next 64 samples of silence");
    }

    // The same shape for the Svf, asked of the output because the state is private — and asserting
    // BOTH halves, because the atomic branch has to fire when it should and not when it should not.
    // The same 1e37 drive overflows at Q 40 and does not at Q 8 (measured: 3852 non-finite outputs
    // against none), so the two Q values pin the two directions. A test that demanded "reset" of both
    // would be demanding that a healthy filter throw its history away.
    for (double q : { 8.0, 40.0 })
    {
        eq::Svf f, fresh;
        f.prepare (kFs, 2); fresh.prepare (kFs, 2);
        f.setParams (eq::FilterType::BandPass, 7000.0, q, 0.0);
        fresh.setParams (eq::FilterType::BandPass, 7000.0, q, 0.0);
        bool wentBad = false;
        for (int i = 0; i < 4000; ++i)
            if (! std::isfinite (f.processSample (0, (float) (1.0e37 * std::sin (2.0 * core::kPi * 7000.0 * i / kFs)))))
                wentBad = true;
        f.flushDenormals();
        int diff = 0;
        for (int k = 0; k < 64; ++k)
        {
            const float a = f.processSample (0, 0.0f), b = fresh.processSample (0, 0.0f);
            if (std::memcmp (&a, &b, 4) != 0) ++diff;
        }
        const std::string at = "Svf at Q " + std::to_string ((int) q);
        if (wentBad)
        {
            test::ok (q > 20.0, at + ": a finite 1e37 drive DOES overflow it here");
            test::ok (diff == 0, at + ": overflowed, then flushed, is bit-for-bit a reset filter");
        }
        else
        {
            test::ok (q < 20.0, at + ": the same drive leaves it finite here");
            test::ok (diff > 0, at + ": a HEALTHY filter keeps its history — the flush must not reset it");
        }
    }

    // A CONSTANT drive is what actually splits the Svf's pair — a swept sine does not, which is how a
    // 25920-point sweep of sines and squares concluded, wrongly, that the case was unreachable here.
    // fc 100 Hz at Q 2, 162 samples of 3e38: with the pair cleared word by word the next output on
    // SILENCE is 1.05e36; cleared atomically it is zero.
    {
        eq::Svf f, fresh;
        f.prepare (kFs, 2); fresh.prepare (kFs, 2);
        f.setParams (eq::FilterType::LowPass, 100.0, 2.0, 0.0);
        fresh.setParams (eq::FilterType::LowPass, 100.0, 2.0, 0.0);
        for (int i = 0; i < 162; ++i) f.processSample (0, 3.0e38f);
        f.flushDenormals();
        int diff = 0;
        for (int k = 0; k < 64; ++k)
        {
            const float a = f.processSample (0, 0.0f), b = fresh.processSample (0, 0.0f);
            if (std::memcmp (&a, &b, 4) != 0) ++diff;
        }
        test::ok (diff == 0, "Svf, constant 3e38 at fc 100 / Q 2: the pair is cleared atomically, so the "
                             "flushed filter is bit-for-bit a reset one");
    }

    // EVERY channel recovers, and it is checked for EVERY channel. A guard applied to channel 0 alone
    // — or a loop stopping one short — survives any test that poisons only channel 0 and treats the
    // rest as innocent bystanders, which is what the neighbour check above does by design.
    for (int c = 0; c < core::kMaxChannels; ++c)
    {
        eq::Svf f, fresh;
        f.prepare (kFs, core::kMaxChannels); fresh.prepare (kFs, core::kMaxChannels);
        f.setParams (eq::FilterType::LowPass, 3000.0, 0.707, 0.0);
        fresh.setParams (eq::FilterType::LowPass, 3000.0, 0.707, 0.0);
        for (int i = 0; i < 256; ++i) f.processSample (c, toneAt (i));
        f.processSample (c, kInf);
        f.flushDenormals();
        int diff = 0;
        for (int i = 0; i < 200; ++i)
        {
            const float a = f.processSample (c, toneAt (i));
            const float b = fresh.processSample (c, toneAt (i));
            if (std::memcmp (&a, &b, 4) != 0) ++diff;
        }
        test::ok (diff == 0, "channel " + std::to_string (c) + " recovers to a reset filter");
    }
}

//==============================================================================
// A finite input can poison a follower without going anywhere near 1e38, because the Rms state is a
// SQUARE: 1e20 * 1e20 overflows in float. Worth its own check, because the natural way to describe
// the bit-exactness claim — "only inputs above 1e38 differ" — is wrong for exactly this reason.
static void aFiniteInputCanPoisonThroughTheSquare()
{
    test::group ("the square is a poisoning route of its own — 1e20, not 1e38");

    dynamics::EnvelopeFollower f;
    f.prepare (kFs); f.setDetector (dynamics::Detector::Rms); f.setTimes (5.0, 5.0);
    for (int i = 0; i < 100; ++i) f.process (0.3f);
    f.process (1.0e20f);                                     // finite, positive, not alternating
    test::ok (! std::isfinite (f.envelope()),
              "a FINITE 1e20 leaves the Rms state non-finite — the square overflowed, not the input");
    f.flushDenormals();
    float last = 0.0f;
    for (int i = 0; i < 2000; ++i) last = f.process (0.3f);
    test::approx ((double) last, 0.3, 0.01, "and the flush brings it back to tracking the signal");

    // Peak mode does NOT square, so the same input is harmless there — which is the other half of
    // saying where the route actually is.
    dynamics::EnvelopeFollower p;
    p.prepare (kFs); p.setDetector (dynamics::Detector::Peak); p.setTimes (5.0, 5.0);
    for (int i = 0; i < 100; ++i) p.process (0.3f);
    p.process (1.0e20f);
    test::ok (std::isfinite (p.envelope()), "Peak mode does not square, so 1e20 leaves it finite");
}


//==============================================================================
// THE CONSUMERS, which is where the promise is actually kept or broken. Hardening four primitives
// does nothing for a product unless the module that owns them calls the flush — and until these
// checks existed, every one of those calls could be DELETED without a single suite going red.
// Mutations that survived the whole 18-binary tree: DynamicEqBand dropping `audio_.flushDenormals()`,
// dropping `side_` and `env_`, DeEsser dropping `xover_.flushDenormals()`, TransientShaper dropping
// both of its follower flushes. The contract is per module, so the test has to be per module.
static void theConsumersRecoverWithinOneBlock()
{
    test::group ("consumers — one bad sample, and the NEXT block is clean again");

    const double fs = kFs;
    const int n = 8192, blk = 512;
    auto tone = [] (int i) { return (float) (0.4 * std::sin (2.0 * core::kPi * 700.0 * i / kFs)); };

    // TransientShaper. Mono, because the stereo Max link would hide the poison — see the note below.
    for (auto link : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
    {
        std::vector<float> x ((std::size_t) n);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = tone (i);
        x[100] = kNan;
        dynamics::TransientShaper t; t.prepare (fs, blk, 1);
        dynamics::TransientShaperParams p; p.attackDb = 6.0; p.sustainDb = -3.0; p.mix = 1.0; p.link = link;
        t.setParams (p);
        for (int off = 0; off < n; off += blk) { float* io[1] { x.data() + off }; t.process (io, 1, blk); }
        int afterFirstBlock = 0; double energy = 0.0;
        for (int i = blk; i < n; ++i)
        {
            if (! std::isfinite (x[(std::size_t) i])) ++afterFirstBlock;
            else energy += std::fabs ((double) x[(std::size_t) i]);
        }
        const std::string nm = std::string ("TransientShaper (")
                             + (link == dynamics::LinkMode::Max ? "Max" : "MeanPower") + ")";
        test::ok (afterFirstBlock == 0, nm + ": everything after the poisoned block is finite ("
                  + std::to_string (afterFirstBlock) + ")");
        // ...AND THE SIGNAL IS BACK, measured against an unpoisoned run rather than against zero.
        // Clearing the de-zipper to zero instead of unity does not silence the stream for ever — the
        // one-pole climbs back over its own smoothing time — so a total-energy check cannot see it.
        // What it does is punch a hole: measured, the worst sample right after the flush is 0.030 of
        // the reference against 0.777 with unity. That is 30 dB, and it is what this pins.
        std::vector<float> ref ((std::size_t) n);
        for (int i = 0; i < n; ++i) ref[(std::size_t) i] = tone (i);
        {
            dynamics::TransientShaper t2; t2.prepare (fs, blk, 1); t2.setParams (p);
            for (int off = 0; off < n; off += blk) { float* io[1] { ref.data() + off }; t2.process (io, 1, blk); }
        }
        double worst = 1.0;
        for (int i = blk; i < blk + 2048 && i < n; ++i)
        {
            const double r = std::fabs ((double) ref[(std::size_t) i]);
            if (r > 0.05) worst = std::min (worst, std::fabs ((double) x[(std::size_t) i]) / r);
        }
        test::ok (worst > 0.5, nm + ": and no 30 dB hole where the de-zipper restarted (worst |y|/|ref| "
                  + std::to_string (worst) + "; clearing it to zero instead of unity gives 0.030)");
        test::ok (energy > 0.2 * 0.4 * (n - blk), nm + ": and the audio is there at all");
    }

    // ...and the reason the mono fixture is not an accident: a stereo Max link DROPS a NaN by
    // comparison (`NaN > mx` is false), so a stereo fixture would pass while proving nothing. Pinned,
    // so that nobody later "simplifies" the fixture into a blind one.
    {
        std::vector<float> a ((std::size_t) n), b ((std::size_t) n);
        for (int i = 0; i < n; ++i) a[(std::size_t) i] = b[(std::size_t) i] = tone (i);
        a[100] = kNan;
        dynamics::TransientShaper t; t.prepare (fs, blk, 2);
        dynamics::TransientShaperParams p; p.attackDb = 6.0; p.sustainDb = -3.0; p.link = dynamics::LinkMode::Max;
        t.setParams (p);
        for (int off = 0; off < n; off += blk) { float* io[2] { a.data() + off, b.data() + off }; t.process (io, 2, blk); }
        int bad = 0; for (int i = 0; i < n; ++i) if (! std::isfinite (b[(std::size_t) i])) ++bad;
        test::ok (bad == 0, "a stereo Max link drops a one-channel NaN by comparison — which is why the "
                            "fixture above is mono, and this is pinned so it stays that way");
    }

    // DeEsser — BOTH modes, because they are different signal paths and the default one returns early.
    // `DynamicEq` delegates to a dynamic-EQ band and never touches `xover_`, so a fixture that only
    // ran the default would bless a build with the crossover flush deleted. `SplitBand` is where that
    // crossover lives, in the AUDIO path, where a gate is not allowed and the flush is the only thing
    // standing between one sample and a permanently dead band.
    for (auto mode : { deesser::DeEsserMode::DynamicEq, deesser::DeEsserMode::SplitBand })
    {
        std::vector<float> x ((std::size_t) n);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 6500.0 * i / kFs));
        x[100] = kInf;
        deesser::DeEsser d; d.prepare (fs, blk, 1);
        deesser::DeEsserParams dp; dp.mode = mode; d.setParams (dp);
        for (int off = 0; off < n; off += blk) { float* io[1] { x.data() + off }; d.process (io, 1, blk); }
        int after = 0; double energy = 0.0;
        for (int i = blk; i < n; ++i)
        {
            if (! std::isfinite (x[(std::size_t) i])) ++after;
            else energy += std::fabs ((double) x[(std::size_t) i]);
        }
        const std::string nm = std::string ("DeEsser ")
                             + (mode == deesser::DeEsserMode::DynamicEq ? "DynamicEq" : "SplitBand");
        test::ok (after == 0, nm + ": everything after the poisoned block is finite (" + std::to_string (after) + ")");
        test::ok (energy > 0.2 * 0.4 * (n - blk), nm + ": and the audio came back");
        test::ok (std::isfinite (d.gainReductionDb()), nm + ": and the reported gain reduction is a number again");
    }

    // DynamicEqBand — two filters, one in the sidechain and one in the audio path, and the failure
    // here is INVISIBLE to a finiteness check. A poisoned detector reaches `std::max(NaN, 1e-9f)`,
    // which is NaN, and `core::gainToDb` turns that into its floor of -240 dB — a perfectly finite
    // number. The audio stays clean and the band simply stops responding for ever. So the check is
    // that the delta still MOVES with the signal, measured against a run that was never poisoned.
    {
        auto run = [&] (bool poison)
        {
            std::vector<float> x ((std::size_t) n);
            for (int i = 0; i < n; ++i)
            {
                const double amp = (i > n / 2) ? 0.7 : 0.02;          // a quiet half, then a loud one
                x[(std::size_t) i] = (float) (amp * std::sin (2.0 * core::kPi * 400.0 * i / kFs));
            }
            if (poison) x[100] = kInf;
            dynamiceq::DynamicEqBand e; e.prepare (fs, 1);
            int after = 0;
            for (int off = 0; off < n; off += blk)
            {
                float* io[1] { x.data() + off }; e.process (io, 1, blk);
                if (off >= blk) for (int i = 0; i < blk; ++i) if (! std::isfinite (x[(std::size_t) (off + i)])) ++after;
            }
            return std::pair<int, double> { after, e.dynamicDeltaDb() };
        };
        const auto clean = run (false), dirty = run (true);
        test::ok (dirty.first == 0, "DynamicEqBand: everything after the poisoned block is finite ("
                  + std::to_string (dirty.first) + ")");
        test::approx (dirty.second, clean.second, 0.5,
                      "and the band is still RESPONDING — its delta matches an unpoisoned run ("
                      + std::to_string (dirty.second) + " against " + std::to_string (clean.second)
                      + " dB). A finiteness check cannot see this: a dead detector reads -240 dB, which is finite.");
    }

    // AutoLeveler fails SILENTLY, so the check has to be about movement, not about finiteness: the
    // audio stays perfectly good while the makeup stops tracking for ever. Measured before the fix:
    // frozen at 0.000 dB across 4000 healthy blocks that were asking for +14.
    {
        dynamics::AutoLeveler al; al.prepare (fs);
        auto drive = [&] (double dry, double mix, int blocks)
        { for (int k = 0; k < blocks; ++k) { al.processBlock (dry, mix, true, blk); for (int i = 0; i < blk; ++i) al.getNextGain(); } };
        drive (0.01, 0.01, 200);
        const double before = al.currentGainDb();
        al.processBlock (std::numeric_limits<double>::quiet_NaN(), 0.01, true, blk);
        drive (0.25, 0.01, 4000);
        test::ok (std::fabs (al.currentGainDb() - before) > 10.0,
                  "AutoLeveler still TRACKS after a poisoned block (" + std::to_string (al.currentGainDb())
                  + " dB, was " + std::to_string (before) + ") — its failure mode is a frozen target, not a NaN");
        test::approx (al.currentGainDb(), 13.98, 0.5, "and it reaches the level the energies ask for");
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
    partialPoisoningIsHealedAtomically();
    aFiniteInputCanPoisonThroughTheSquare();
    theConsumersRecoverWithinOneBlock();
    return felitronics::test::report();
}
