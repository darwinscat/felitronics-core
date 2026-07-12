// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Adversarial falsification pass for AutoLeveler: explicit external constants plus edge schedules
// around the time-based transition window.

#include <felitronics_test.h>
#include <felitronics/dynamics/AutoLeveler.h>

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::dynamics::AutoLeveler;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kSr = 48000.0;
constexpr float  kMaxGain = 63.10f;
constexpr float  kMinGain = 0.0631f;

float gainDb (float g)
{
    return g > 0.0f ? std::max (-120.0f, 20.0f * std::log10 (g)) : -120.0f;
}

void runFor (AutoLeveler& a, double dry, double mix, bool enabled, int block, double seconds)
{
    const int blocks = (int) std::llround (seconds * kSr / block);
    for (int b = 0; b < blocks; ++b)
    {
        a.processBlock (dry, mix, enabled, block);
        for (int i = 0; i < block; ++i) a.getNextGain();
    }
}
} // namespace

int main()
{
    std::printf ("felitronics::dynamics AutoLeveler FALSIFY tests\n");

    group ("DeepSeek 5.1 materialized: dry==wet at -18 dBFS converges to unity");
    {
        const double ms = std::pow (10.0, -18.0 / 10.0);
        AutoLeveler a; a.prepare (kSr);
        runFor (a, ms, ms, true, 64, 2.0);
        approx (a.currentGain(), 1.0, 1.0e-3, "dry==wet at -18 dBFS mean-square holds unity gain");
    }

    group ("DeepSeek 5.2 materialized: wet 50 dB below dry clamps at kMatchMaxGain, never above");
    {
        AutoLeveler a; a.prepare (kSr);
        float peak = 0.0f;
        for (int b = 0; b < (int) std::llround (7.0 * kSr / 64.0); ++b)
        {
            a.processBlock (1.0, 1.0e-5, true, 64);     // sqrt ratio = +50 dB, beyond +36 dB clamp
            for (int i = 0; i < 64; ++i) peak = std::max (peak, a.getNextGain());
        }
        approx (a.currentGain(), kMaxGain, 0.02, "settled gain lands on the +36 dB clamp");
        ok (peak <= kMaxGain * 1.0001f, "applied gain never overshoots kMatchMaxGain");
    }

    group ("DeepSeek 5.3 adapted: a 28 dB target jump cannot exceed 4.5 dB after 0.5 s");
    {
        AutoLeveler a; a.prepare (kSr);
        runFor (a, 1.0, 1.0, true, 64, 2.0);
        runFor (a, 1.0, std::pow (10.0, -28.0 / 10.0), true, 64, 0.5);
        const float ceiling = std::pow (10.0f, 4.5f / 20.0f);
        ok (a.currentGain() <= ceiling * 1.002f, "0.5 s normal slew stays <= +4.5 dB");
    }

    group ("DeepSeek 5.4/5.5 adapted: transition windows are time-based, not error-threshold based");
    {
        AutoLeveler a; a.prepare (kSr);
        a.openTransitionWindow();
        int fastSamples = 0;
        while (a.snapGliding())
        {
            a.getNextGain();
            ++fastSamples;
        }
        ok (fastSamples == (int) std::lround (0.35 * kSr), "openTransitionWindow lasts exactly 0.35 seconds of samples");

        AutoLeveler b; b.prepare (kSr);
        b.openTransitionWindow();
        b.processBlock (1.0, 1.0e-3, true, (int) kSr);   // one hostile 1-second host block
        float prev = b.currentGain();
        bool slopeOk = true;
        for (int i = 0; i < (int) kSr; ++i)
        {
            const bool fast = b.snapGliding();
            const float next = b.getNextGain();
            const double stepDb = std::fabs (20.0 * std::log10 ((double) next / (double) prev));
            const double limit = (fast ? 40.0 : 9.0) / kSr + 1.0e-5;
            slopeOk = slopeOk && (stepDb <= limit);
            prev = next;
        }
        ok (slopeOk, "single huge transition block still obeys fast-then-normal per-sample ceilings");
    }

    group ("DeepSeek 5.6 adapted: dry silence with wet continuing holds after the dry floor, never rails");
    {
        AutoLeveler a; a.prepare (kSr);
        runFor (a, 1.0, 0.25, true, 64, 2.0);            // settle near +6 dB
        approx (a.currentGainDb(), 6.0, 0.1, "precondition: settled at +6 dB");

        runFor (a, 0.0, 0.25, true, 64, 3.0);            // dry follower decays to floor; target then holds
        const float heldDb = a.currentGainDb();
        runFor (a, 0.0, 0.25, true, 64, 3.0);
        approx (a.currentGainDb(), heldDb, 0.15, "after dry falls below floor, wet-only input holds the makeup");
        ok (a.currentGain() > kMinGain * 1.5f, "wet-only after dry silence does not slam to the -24 dB rail");
        ok (a.currentGain() < kMaxGain * 0.25f, "wet-only after dry silence does not slam upward to +36 dB");
    }

    group ("own attack: prepare() sample-rate fallback is deterministic below 1 kHz");
    {
        AutoLeveler a, b;
        a.prepare (999.0);
        b.prepare (48000.0);
        std::vector<float> ga, gb;
        for (int blk = 0; blk < 500; ++blk)
        {
            a.processBlock (1.0, 0.001, true, 64);
            b.processBlock (1.0, 0.001, true, 64);
            for (int i = 0; i < 64; ++i) { ga.push_back (a.getNextGain()); gb.push_back (b.getNextGain()); }
        }
        ok (ga == gb, "prepare(sr<=1000) follows the documented 48 kHz fallback bit-for-bit");
    }

    return felitronics::test::report();
}
