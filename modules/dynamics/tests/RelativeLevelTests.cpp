// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::dynamics::RelativeLevel. Two carry most of the weight:
// LEVEL TRANSLATION INVARIANCE (the "one knob, any signal" promise) and NO RECTIFICATION BIAS (an
// asymmetric follower measured +3.9 dB high on a +-12 dB swing — this suite fails such a build).

#include <felitronics_test.h>
#include <felitronics/dynamics/RelativeLevel.h>
#include <felitronics/core/Math.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

static double feed (dynamics::RelativeLevel& rl, double fs, double amp, double ms, int K)
{
    const int n = (int) std::lround (ms * 1.0e-3 * fs);
    for (int i = 0; i < n; ++i)
    {
        rl.accumulate ((float) amp);
        if ((i % K) == K - 1) rl.update (K);
    }
    return rl.programDb();
}

// A programme swinging +-swingDb about trueDb — the adversarial case for any asymmetric follower.
static double swungEstimate (double fs, double trueDb, double swingDb, double modHz, double sec, int K)
{
    dynamics::RelativeLevel rl; rl.prepare (fs);
    const int n = (int) (fs * sec);
    for (int i = 0; i < n; ++i)
    {
        const double db = trueDb + swingDb * std::sin (2.0 * core::kPi * modHz * (double) i / fs);
        rl.accumulate ((float) core::dbToGain (db));
        if ((i % K) == K - 1) rl.update (K);
    }
    return rl.programDb();
}

int main()
{
    std::printf ("felitronics::dynamics::RelativeLevel tests\n");
    const double fs = 48000.0;
    const double e1 = std::exp (1.0);

    test::group ("level-translation invariance (the one-knob promise)");
    {
        double rel[3] = { 0, 0, 0 };
        int idx = 0;
        for (double refDb : { -70.0, -45.0, -20.0 })
        {
            dynamics::RelativeLevel rl; rl.prepare (fs);
            const double est = feed (rl, fs, core::dbToGain (refDb), 12000.0, 16);
            test::approx (est, refDb, 0.5, "estimate converges to the programme level");
            rel[idx++] = rl.relativeDb (refDb + 9.0);
        }
        test::approx (rel[1], rel[0], 0.05, "relative reading invariant across 25 dB");
        test::approx (rel[2], rel[0], 0.05, "relative reading invariant across 50 dB");
        test::approx (rel[0], 9.0, 0.5, "and it equals the actual excess");
    }

    // An asymmetric rise/fall estimator settles ABOVE the truth when the observation swings about it.
    // Measured +3.9 dB at +-12 dB / 8 Hz before this was made symmetric.
    test::group ("no rectification bias on a swinging programme");
    {
        for (double swing : { 2.0, 6.0, 12.0 })
            for (double hz : { 2.0, 8.0 })
            {
                const double est = swungEstimate (fs, -20.0, swing, hz, 40.0, 16);
                test::ok (std::fabs (est + 20.0) < 0.5,
                          "estimate stays within 0.5 dB of the true mean under a symmetric swing");
            }
    }

    // A gap is not a quieter programme: the estimate must barely move through a short one, yet a
    // genuinely quieter section must still arrive rather than latching forever.
    test::group ("gap slowdown: short pauses barely move it, long quiet sections still land");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        feed (rl, fs, core::dbToGain (-20.0), 12000.0, 16);
        const double before = rl.programDb();
        feed (rl, fs, core::dbToGain (-50.0), 400.0, 16);                 // 400 ms inter-phrase gap
        test::ok (before - rl.programDb() < 1.0, "a 400 ms gap moves the estimate less than 1 dB");

        feed (rl, fs, core::dbToGain (-50.0), 60000.0, 16);               // a genuinely quiet minute
        test::ok (rl.programDb() < -40.0, "a sustained quieter section is eventually learned");
    }

    // Every other test feeds a CONSTANT envelope, where mean-square and mean-amplitude agree. An
    // alternating envelope separates them: rms(0,1) = 1/sqrt(2) = -3.0103 dB, mean amplitude = -6.02.
    // This is what falsifies a build that averages amplitude and calls it rms.
    test::group ("the observable really is mean SQUARE, not mean amplitude");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        for (int i = 0; i < 48000 * 12; ++i)
        {
            rl.accumulate ((i & 1) ? 1.0f : 0.0f);
            if ((i % 16) == 15) rl.update (16);
        }
        test::approx (rl.programDb(), -3.0103, 0.05, "alternating 0/1 envelope reads as rms, not mean");
    }

    test::group ("control-rate update is interval-correct (K does not change the time constant)");
    {
        double at[3]; int idx = 0;
        for (int K : { 1, 16, 64 })
        {
            dynamics::RelativeLevel rl; rl.prepare (fs);
            feed (rl, fs, core::dbToGain (-60.0), 12000.0, K);
            at[idx++] = feed (rl, fs, core::dbToGain (-40.0), 2000.0, K);
        }
        test::approx (at[1], at[0], 0.25, "K=16 tracks like K=1");
        test::approx (at[2], at[0], 0.25, "K=64 tracks like K=1");
    }

    // Every other test runs at 48k; a hardcoded sample rate would hide here.
    test::group ("sample-rate invariance");
    {
        double at[3]; int idx = 0;
        for (double rate : { 44100.0, 48000.0, 96000.0 })
        {
            dynamics::RelativeLevel rl; rl.prepare (rate);
            feed (rl, rate, core::dbToGain (-60.0), 12000.0, 16);
            at[idx++] = feed (rl, rate, core::dbToGain (-40.0), 2000.0, 16);
        }
        test::approx (at[0], at[1], 0.25, "44.1k tracks like 48k in WALL CLOCK");
        test::approx (at[2], at[1], 0.25, "96k tracks like 48k in WALL CLOCK");
    }

    test::group ("seeding: the first valid observation sets the estimate outright");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        test::ok (! rl.seeded(), "starts unseeded");
        test::approx (rl.relativeDb (-30.0), 0.0, 1.0e-12, "pre-seed relative reading is transparent");
        for (int i = 0; i < 16; ++i) rl.accumulate ((float) core::dbToGain (-30.0));
        rl.update (16);
        test::ok (rl.seeded(), "seeded on first update");
        test::approx (rl.programDb(), -30.0, 0.1, "seeded AT the observation, not crawling from zero");
    }

    // The first observation after reset() is a detector still climbing its own attack ramp, so it
    // reads low; reset() arms the fast window so that error does not linger for a slew-capped second.
    test::group ("a low first observation self-corrects quickly (reset arms the fast window)");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        for (int i = 0; i < 16; ++i) rl.accumulate ((float) core::dbToGain (-60.0));   // ramp artefact
        rl.update (16);
        test::approx (rl.programDb(), -60.0, 0.1, "seeded low, as a real detector ramp would");
        feed (rl, fs, core::dbToGain (-25.0), 400.0, 16);
        test::ok (rl.programDb() > -38.0, "recovers most of a 35 dB seeding error within 400 ms");

        // Without the armed window the same error would still be crawling at the normal 2 s constant.
        dynamics::RelativeLevel slow; slow.prepare (fs);
        dynamics::RelativeLevel::Params noWindow; noWindow.fastWindowMs = 0.0;
        slow.setParams (noWindow); slow.reset();
        for (int i = 0; i < 16; ++i) slow.accumulate ((float) core::dbToGain (-60.0));
        slow.update (16);
        feed (slow, fs, core::dbToGain (-25.0), 400.0, 16);
        test::ok (slow.programDb() < rl.programDb() - 5.0, "and the window is what made that happen");
    }

    test::group ("slew cap binds in BOTH directions");
    {
        dynamics::RelativeLevel::Params p; p.slewDbPerSec = 6.0; p.timeMs = 1.0; p.fastWindowMs = 0.0;
        dynamics::RelativeLevel rl; rl.prepare (fs); rl.setParams (p); rl.reset();
        feed (rl, fs, core::dbToGain (-50.0), 3000.0, 16);
        const double before = rl.programDb();
        feed (rl, fs, core::dbToGain (-10.0), 100.0, 16);
        test::ok (rl.programDb() - before <= 6.0 * 0.1 + 0.2, "upward motion capped at slewDbPerSec");
        const double high = rl.programDb();
        feed (rl, fs, core::dbToGain (-60.0), 100.0, 16);
        test::ok (high - rl.programDb() <= 6.0 * 0.1 + 0.2, "downward motion capped too");
    }

    test::group ("silence HOLDS the estimate, and the activity flag has hysteresis");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        feed (rl, fs, core::dbToGain (-30.0), 12000.0, 16);
        const double held = rl.programDb();
        test::ok (rl.activity(), "active while the signal runs");
        feed (rl, fs, 0.0, 5000.0, 16);
        test::approx (rl.programDb(), held, 1.0e-9, "estimate bit-held through the pause");
        test::ok (! rl.activity(), "reports inactive so the consumer can gate its computer input");

        // Just above the floor but below floor+hysteresis: must NOT re-activate and chatter.
        const double p = rl.params().floorDb + 0.5 * rl.params().activityHystDb;
        feed (rl, fs, core::dbToGain (p), 200.0, 16);
        test::ok (! rl.activity(), "hysteresis keeps it inactive inside the band");
        feed (rl, fs, core::dbToGain (rl.params().floorDb + 12.0), 200.0, 16);
        test::ok (rl.activity(), "re-activates once clearly above the floor");
    }

    // A programme under the estimate's own clamp cannot be represented, so it must HOLD rather than
    // pin the estimate and silently break the invariance above.
    test::group ("no grey zone between the floor and the estimate clamp");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        test::ok (rl.params().floorDb >= rl.params().clampLoDb - 1.0e-12,
                  "floor is never below the clamp (sanitised)");
        dynamics::RelativeLevel::Params bad; bad.floorDb = -140.0; bad.clampLoDb = -80.0;
        rl.setParams (bad);
        test::ok (rl.params().floorDb >= rl.params().clampLoDb - 1.0e-12, "a hostile floor is raised");
        rl.reset();
        feed (rl, fs, core::dbToGain (-100.0), 3000.0, 16);
        test::ok (! rl.activity(), "a programme below the clamp reads inactive, not pinned-and-active");
    }

    test::group ("retune opens a bounded fast window that actually EXPIRES");
    {
        const double amp = core::dbToGain (-25.0);
        auto after = [&] (bool retune, double ms)
        {
            dynamics::RelativeLevel rl; rl.prepare (fs);
            feed (rl, fs, core::dbToGain (-60.0), 12000.0, 16);      // settle, window long spent
            if (retune) rl.retuned();
            return feed (rl, fs, amp, ms, 16);
        };
        test::ok (after (true, 200.0) > after (false, 200.0) + 1.0, "retuned estimator adapts faster");
        test::ok (after (true, 200.0) < -25.0 - 1.0e-9, "but does NOT jump to the observation");

        // Expiry, tested directly: a 500 ms window and a 5 s window must reach DIFFERENT places after
        // 2 s. A fastLeft_ that never decrements makes both windows effectively infinite and equal,
        // so this is what falsifies that mutation.
        auto withWindow = [&] (double windowMs)
        {
            dynamics::RelativeLevel rl; rl.prepare (fs);
            dynamics::RelativeLevel::Params p; p.fastWindowMs = windowMs;
            rl.setParams (p); rl.reset();
            feed (rl, fs, core::dbToGain (-60.0), 12000.0, 16);
            rl.retuned();
            return feed (rl, fs, amp, 2000.0, 16);
        };
        test::ok (withWindow (5000.0) > withWindow (500.0) + 1.0,
                  "a longer window gets further — so the short one really did expire");

        dynamics::RelativeLevel drag; drag.prepare (fs);
        feed (drag, fs, core::dbToGain (-60.0), 12000.0, 16);
        for (int i = 0; i < 200; ++i) { drag.retuned(); feed (drag, fs, amp, 1.0, 16); }
        test::ok (drag.programDb() < -25.0 - 1.0e-9, "spamming retune still cannot make it instantaneous");
    }

    test::group ("contract abuse cannot corrupt state");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        rl.update (16);                                   // nothing accumulated
        test::ok (! rl.seeded() && std::isfinite (rl.programDb()), "update with no samples is a no-op");
        for (int i = 0; i < 16; ++i) rl.accumulate ((float) core::dbToGain (-30.0));
        rl.update (16);
        rl.update (16);                                   // second update, accumulator already drained
        test::approx (rl.programDb(), -30.0, 0.1, "a repeated update does not move the estimate");
        rl.update (0);
        rl.update (-5);
        test::ok (std::isfinite (rl.programDb()), "zero/negative sample counts are ignored");
    }

    test::group ("poison in, finite out");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        feed (rl, fs, core::dbToGain (-30.0), 6000.0, 16);
        const double good = rl.programDb();
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const float inf = std::numeric_limits<float>::infinity();
        for (int i = 0; i < 64; ++i) { rl.accumulate (nan); rl.accumulate (inf); }
        rl.update (128);
        rl.flushDenormals();
        test::ok (std::isfinite (rl.programDb()), "an all-poison block cannot poison the estimate");
        test::approx (rl.programDb(), good, 1.0e-9, "it is rejected, not averaged in");

        // Mixed block: the mean must come from the VALID samples while dt still counts the whole span.
        for (int i = 0; i < 64; ++i) { rl.accumulate ((float) core::dbToGain (-30.0)); rl.accumulate (nan); }
        rl.update (128);
        test::approx (rl.programDb(), good, 0.5, "a half-poisoned block reads like the valid half");
    }

    test::group ("RT: no allocation in accumulate/update");
    {
        dynamics::RelativeLevel rl; rl.prepare (fs);
        feed (rl, fs, core::dbToGain (-25.0), 100.0, 16);
        const long before = g_allocs.load();
        feed (rl, fs, core::dbToGain (-25.0), 2000.0, 16);
        rl.retuned();
        feed (rl, fs, 0.0, 500.0, 16);
        rl.flushDenormals();
        test::okNoAlloc (g_allocs.load() == before, "no allocations on the audio path");
    }

    return test::report();
}
