// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::basscomp::BassComp — the bass amp's COMP (moved from OrbitAmp).
//
//   BLEND 0, glided there: nothing computed, the block passes bit for bit; a steady note at the reference level keeps
//   its level at any BLEND in every character (the make-up is what the squeeze took); a threshold no note reaches
//   leaves the note alone; the output does not depend on how the host cuts its blocks (the cadence is counted from the
//   prepare); standing by, settled, it is fed and not heard; a prepare's snapAmount plays at once, a recall up from 0
//   warms before it is heard; the indicator reads the squeeze heard; the picture is there only while pictured and heard,
//   and says when its trace breaks; and nothing allocates in process().

#include <felitronics_test.h>
#include <alloc_counter.h>
#include <felitronics/basscomp/BassComp.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace felitronics;
using basscomp::BassComp;

namespace
{
    constexpr double fs = 48000.0;
    constexpr int    block = 512;

    const core::BypassFade::Span heard { 1.0f, 1.0f, false, 0 };
    const core::BypassFade::Span standingBy { 0.0f, 0.0f, false, 0 };

    /** A steady bass note: A1 and its two harmonics, at `rmsDb` dBFS RMS. */
    std::vector<float> note (double rmsDb, int n, double hz = 55.0)
    {
        std::vector<float> x ((size_t) n);
        const double norm = std::sqrt ((1.0 + 0.25 + 0.0625) / 2.0);
        const double a = std::pow (10.0, rmsDb / 20.0) / norm;
        for (int i = 0; i < n; ++i)
        {
            const double t = (double) i / fs;
            x[(size_t) i] = (float) (a * (std::sin (2.0 * core::kPi * hz * t) + 0.5 * std::sin (2.0 * core::kPi * 2.0 * hz * t)
                                          + 0.25 * std::sin (2.0 * core::kPi * 3.0 * hz * t)));
        }
        return x;
    }

    double rmsDb (const std::vector<float>& x, size_t from, size_t to)
    {
        double s = 0.0;
        for (size_t i = from; i < to; ++i)
            s += (double) x[i] * x[i];
        return 10.0 * std::log10 (s / (double) (to - from));
    }

    /** `x` through a fresh COMP in host blocks of `blk`, at a BLEND and settings that stand. */
    std::vector<float> run (std::vector<float> x, int blk, double blend, BassComp::Character c,
                            double thr = BassComp::defaultThresholdDb, double ratio = BassComp::defaultRatio,
                            const core::BypassFade::Span& link = heard)
    {
        BassComp comp;
        (void) comp.prepare (fs, block);
        comp.snapAmount (blend);
        for (size_t i = 0; i < x.size(); i += (size_t) blk)
        {
            const int n = (int) std::min<size_t> ((size_t) blk, x.size() - i);
            comp.process (x.data() + i, n, blend, link, BassComp::Attack::med, BassComp::Release::automatic, c, thr, ratio);
        }
        return x;
    }
}

int main()
{
    std::printf ("felitronics::basscomp BassComp tests\n");

    const auto x = note (BassComp::referenceDb, (int) (4.0 * fs));

    test::group ("BLEND 0: nothing computed, bit for bit");
    {
        const auto y = run (x, block, 0.0, BassComp::Character::opto);
        test::ok (std::memcmp (x.data(), y.data(), sizeof (float) * x.size()) == 0, "the note passes unchanged");
    }

    test::group ("the make-up is level-neutral: the reference note keeps its level at any BLEND");
    {
        const size_t from = (size_t) (2.0 * fs), to = x.size();
        for (const auto c : { BassComp::Character::smooth, BassComp::Character::punchy, BassComp::Character::opto,
                              BassComp::Character::fet, BassComp::Character::squash })
            for (const double blend : { 0.5, 1.0 })
            {
                const auto y = run (x, block, blend, c);
                test::approx (rmsDb (y, from, to), rmsDb (x, from, to), 1.2,
                              "character " + std::to_string ((int) c) + ", BLEND " + std::to_string (blend) + ": the note's level");
            }
        const auto deep = run (x, block, 1.0, BassComp::Character::smooth, -40.0, 8.0);
        test::approx (rmsDb (deep, from, to), rmsDb (x, from, to), 1.2, "deep (-40 dB, 8:1), SMOOTH: still its level");
    }

    test::group ("a threshold no note reaches: the note alone");
    {
        const auto quiet = note (-40.0, (int) (2.0 * fs));
        const auto y = run (quiet, block, 1.0, BassComp::Character::smooth, -12.0, 4.0);
        double worst = 0.0;
        for (size_t i = 0; i < y.size(); ++i)
            worst = std::max (worst, std::abs ((double) y[i] - (double) quiet[i]));
        test::ok (worst < 1.0e-5, "a -40 dB note under a -12 dB threshold: within 1e-5 of itself (" + std::to_string (worst) + ")");
    }

    test::group ("its time is its own: the host's blocks do not matter");
    {
        for (const auto c : { BassComp::Character::smooth, BassComp::Character::opto, BassComp::Character::fet })
        {
            const auto a = run (x, 512, 0.7, c, -30.0, 4.0);
            const auto b = run (x, 64, 0.7, c, -30.0, 4.0);
            const auto d = run (x, 333, 0.7, c, -30.0, 4.0);
            test::ok (std::memcmp (a.data(), b.data(), sizeof (float) * a.size()) == 0
                          && std::memcmp (a.data(), d.data(), sizeof (float) * a.size()) == 0,
                      "character " + std::to_string ((int) c) + ": blocks of 512, 64 and 333 play the same samples");
        }
    }

    test::group ("standing by, settled: fed, not heard");
    {
        const auto y = run (x, block, 1.0, BassComp::Character::smooth, -30.0, 4.0, standingBy);
        test::ok (std::memcmp (x.data(), y.data(), sizeof (float) * x.size()) == 0, "the dry, bit for bit");
    }

    test::group ("a prepare plays at once; a recall up from 0 warms first");
    {
        BassComp snapped, recalled;
        (void) snapped.prepare (fs, block);
        (void) recalled.prepare (fs, block);
        snapped.snapAmount (1.0);
        auto a = x, b = x;
        snapped.process (a.data(), block, 1.0, heard, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::smooth, -30.0, 4.0);
        recalled.process (b.data(), block, 1.0, heard, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::smooth, -30.0, 4.0);
        bool recalledDry = true;
        for (int i = 0; i < block; ++i)
            recalledDry = recalledDry && core::exactlyEqual (b[(size_t) i], x[(size_t) i]);
        test::ok (recalledDry, "recalled up from 0: its first block is fed, not heard (the warm)");
        bool snappedWorks = false;
        for (int i = 0; i < block; ++i)
            snappedWorks = snappedWorks || ! core::exactlyEqual (a[(size_t) i], x[(size_t) i]);
        test::ok (snappedWorks, "snapped at a prepare: heard from its first block");
    }

    test::group ("the indicator and the picture");
    {
        BassComp comp;
        (void) comp.prepare (fs, block);
        comp.snapAmount (1.0);
        auto y = x;
        comp.pictured = false;
        comp.process (y.data(), block, 1.0, heard, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::opto, -30.0, 4.0);
        test::ok (comp.picture() == nullptr, "not pictured: no picture");
        comp.pictured = true;
        comp.process (y.data() + block, block, 1.0, heard, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::opto, -30.0, 4.0);
        const auto* p = comp.picture();
        test::ok (p != nullptr && p->n == block && p->breaksTrace && p->rms, "pictured: the block's picture, its trace breaking (OPTO hears the RMS)");
        comp.process (y.data() + 2 * block, block, 1.0, heard, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::opto, -30.0, 4.0);
        p = comp.picture();
        test::ok (p != nullptr && ! p->breaksTrace, "the next one goes on the same trace");
        bool squeezes = p != nullptr;
        for (int i = 0; p != nullptr && i < p->n; ++i)
            squeezes = squeezes && p->grDb[i] <= 0.0f && core::exactlyEqual (p->in[i], x[(size_t) (2 * block + i)]);
        test::ok (squeezes, "IN is the block handed in, and the squeeze is never a lift");
        test::ok (comp.grDb.load() < -0.5f, "the indicator reads the squeeze heard (" + std::to_string (comp.grDb.load()) + " dB)");
        comp.process (y.data() + 3 * block, block, 1.0, standingBy, BassComp::Attack::med, BassComp::Release::automatic, BassComp::Character::opto, -30.0, 4.0);
        test::ok (comp.picture() == nullptr && core::exactlyEqual (comp.grDb.load(), 0.0f), "standing by: no picture, no squeeze shown");
    }

    test::group ("no allocation in process()");
    {
        BassComp comp;
        (void) comp.prepare (fs, block);
        comp.snapAmount (1.0);
        auto y = x;
        const long long before = alloc::count.load();
        for (int b = 0; b < 40; ++b)
            comp.process (y.data() + (size_t) b * block, block, 0.8, heard, BassComp::Attack::fast, BassComp::Release::automatic,
                          b < 20 ? BassComp::Character::smooth : BassComp::Character::squash, -30.0, 4.0);
        test::okNoAlloc (alloc::count.load() == before, "40 blocks, a character switch among them");
    }

    return test::report();
}
