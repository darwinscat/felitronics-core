// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::dynamics::CharacterComp — the bass amp's own compressor characters (moved from OrbitAmp).
//
//   a steady note settles where the static curve says, feed-forward and fed back, on the RMS and on the peaks; a
//   level under the knee is untouched; a two-stage release lets a short hit go at once and a long squeeze slowly;
//   GRAIN adds the third harmonic and leaves the note's level where the squeeze put it; DARK lets the top down; the
//   squeeze it reports is the one applied; seed() takes a squeeze over without a step; and it is the same per sample
//   however its calls are sliced.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/dynamics/CharacterComp.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace felitronics;
using dynamics::CharacterComp;

namespace
{
    constexpr double fs = 48000.0;

    std::vector<float> sine (double hz, double peak, int n)
    {
        std::vector<float> x ((size_t) n);
        for (int i = 0; i < n; ++i)
            x[(size_t) i] = (float) (peak * std::sin (2.0 * core::kPi * hz * (double) i / fs));
        return x;
    }

    double rmsDb (const std::vector<float>& x, size_t from, size_t to)
    {
        double s = 0.0;
        for (size_t i = from; i < to; ++i)
            s += (double) x[i] * x[i];
        return 10.0 * std::log10 (s / (double) (to - from));
    }

    /** The squeeze a steady note settles at, dB (0 or below), its last reported sample. */
    double settledGrDb (CharacterComp::Voice v, double thr, double ratio, double peak)
    {
        CharacterComp c;
        c.prepare (fs);
        c.set (v, thr, ratio, 8.0);
        auto x = sine (220.0, peak, (int) (3.0 * fs));
        const auto key = x;
        std::vector<float> gr (x.size());
        c.process (x.data(), key.data(), (int) x.size(), gr.data());
        double deepest = 0.0;
        for (size_t i = x.size() - 480; i < x.size(); ++i)   // the last 10 ms: two cycles of the ripple
            deepest = std::min (deepest, (double) gr[i]);
        return deepest;
    }
}

int main()
{
    std::printf ("felitronics::dynamics CharacterComp tests\n");

    test::group ("the static curve: where a steady note settles");
    {
        CharacterComp::Voice rms;
        rms.rms = true;
        rms.attackMs = 5.0;
        rms.releaseMs = 200.0;
        // A sine at peak 0.25: RMS -15.05 dBFS; at -30 dB, 4:1, the feed-forward curve takes 15.05 * 3/4 dB.
        const double want = CharacterComp::staticDeltaDb (false, -15.05, -30.0, 4.0, 8.0);
        test::approx (settledGrDb (rms, -30.0, 4.0, 0.25), want, 0.4, "feed-forward on the RMS");
        rms.feedback = true;
        const double fb = CharacterComp::staticDeltaDb (true, -15.05, -30.0, 4.0, 8.0);
        test::approx (settledGrDb (rms, -30.0, 4.0, 0.25), fb, 0.5, "fed back on the RMS");
        test::approx (fb, want, 0.2, "far above the knee the two curves keep the same ratio");
        // Fed back, the knee is R times wider in the input's terms: 2 dB under the threshold it already squeezes, by
        // what the quadratic solved here says ((sqrt (2.5) - 1) / 0.375 - 4 + 2 = -0.4503 dB), and feed-forward less.
        test::approx (CharacterComp::staticDeltaDb (true, -26.0, -24.0, 4.0, 8.0), (std::sqrt (2.5) - 1.0) / 0.375 - 4.0 + 2.0, 1.0e-9,
                      "fed back, in its wider knee: the quadratic's root");
        test::approx (CharacterComp::staticDeltaDb (false, -26.0, -24.0, 4.0, 8.0), -0.1875, 1.0e-12, "feed-forward there: (d + h)^2 / 4h * (1 - 1/R)");
        test::approx (settledGrDb (rms, -12.0, 4.0, 0.05), 0.0, 1.0e-9, "a level under the knee: no squeeze at all");
    }

    test::group ("the two-stage release: a short hit lets go at once, a long squeeze slowly");
    {
        CharacterComp::Voice opto;
        opto.rms = true;
        opto.feedback = true;
        opto.attackMs = 10.0;
        opto.releaseMs = 60.0;
        opto.memoryMs = 2000.0;
        opto.chargeMs = 800.0;
        const auto afterLetGo = [&] (double heldSeconds)
        {
            CharacterComp c;
            c.prepare (fs);
            c.set (opto, -30.0, 4.0, 8.0);
            auto x = sine (220.0, 0.25, (int) (heldSeconds * fs));
            const auto key = x;
            c.process (x.data(), key.data(), (int) x.size(), nullptr);
            std::vector<float> quiet ((size_t) (0.3 * fs), 0.0f);
            const auto k2 = quiet;
            c.process (quiet.data(), k2.data(), (int) quiet.size(), nullptr);
            return c.reductionDb();
        };
        const double shortHit = afterLetGo (0.05), longSqueeze = afterLetGo (3.0);
        test::ok (longSqueeze > 2.0 * shortHit + 0.5, "300 ms after: the long squeeze still holds more (" + std::to_string (longSqueeze)
                                                       + " dB) than the short hit (" + std::to_string (shortHit) + " dB)");
    }

    test::group ("GRAIN: the third harmonic, the note's level kept");
    {
        CharacterComp::Voice plain, fet;
        plain.attackMs = fet.attackMs = 1.0;
        plain.releaseMs = fet.releaseMs = 50.0;
        fet.grain = 0.16;
        const auto run = [&] (const CharacterComp::Voice& v)
        {
            CharacterComp c;
            c.prepare (fs);
            c.set (v, -30.0, 4.0, 8.0);
            auto x = sine (200.0, 0.25, (int) (2.0 * fs));
            const auto key = x;
            c.process (x.data(), key.data(), (int) x.size(), nullptr);
            return x;
        };
        const auto a = run (plain), b = run (fet);
        const size_t from = (size_t) fs, to = (size_t) (2.0 * fs) - 64;
        test::approx (rmsDb (b, from, to), rmsDb (a, from, to), 0.15, "the grain leaves the level where the squeeze put it");
        double diff = 0.0;
        for (size_t i = from; i < to; ++i)
            diff += std::pow ((double) b[i] - (double) a[i], 2.0);
        test::ok (diff > 0.0, "and it is there: the grained note is not the plain one");
    }

    test::group ("DARK: the top let down as it squeezes");
    {
        CharacterComp::Voice plain, dark;
        plain.attackMs = dark.attackMs = 3.0;
        plain.releaseMs = dark.releaseMs = 250.0;
        dark.darkHz = 3000.0;
        const auto top = [&] (const CharacterComp::Voice& v)
        {
            CharacterComp c;
            c.prepare (fs);
            c.set (v, -40.0, 8.0, 8.0);
            auto x = sine (8000.0, 0.5, (int) fs);
            const auto key = x;
            c.process (x.data(), key.data(), (int) x.size(), nullptr);
            return rmsDb (x, x.size() / 2, x.size());
        };
        test::ok (top (dark) < top (plain) - 3.0, "an 8 kHz tone squeezed hard comes out darker");
    }

    test::group ("the squeeze reported is the one applied; seed takes one over; slicing does not matter");
    {
        CharacterComp::Voice v;
        v.attackMs = 5.0;
        v.releaseMs = 100.0;
        CharacterComp a, b;
        a.prepare (fs);
        b.prepare (fs);
        a.set (v, -30.0, 4.0, 8.0);
        b.set (v, -30.0, 4.0, 8.0);
        auto x = sine (110.0, 0.3, 4800);
        auto y = x;
        const auto key = x;
        std::vector<float> ga (x.size()), gb (x.size());
        a.process (x.data(), key.data(), (int) x.size(), ga.data());
        for (size_t i = 0; i < y.size(); i += 37)
        {
            const int n = (int) std::min<size_t> (37, y.size() - i);
            b.process (y.data() + i, key.data() + i, n, gb.data() + i);
        }
        test::ok (std::memcmp (x.data(), y.data(), sizeof (float) * x.size()) == 0 && std::memcmp (ga.data(), gb.data(), sizeof (float) * ga.size()) == 0,
                  "one call and calls of 37: the same samples and the same squeeze");
        test::approx (-(double) ga.back(), a.reductionDb(), 1.0e-6, "reductionDb() is the last sample's squeeze");

        CharacterComp c;
        c.prepare (fs);
        c.set (v, -30.0, 4.0, 8.0);
        c.seed (6.0);
        test::approx (c.reductionDb(), 6.0, 1.0e-12, "seed(6): it holds 6 dB from its first sample");
        c.reset();
        test::approx (c.reductionDb(), 0.0, 1.0e-12, "reset(): nothing held");
    }

    return test::report();
}
