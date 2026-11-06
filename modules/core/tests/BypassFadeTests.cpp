// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::core::BypassFade — the one crossfade of a link that replaces the signal (moved from OrbitAmp).
//
//   settled, it reports no motion and costs nothing; a switch takes 15 ms of audio whatever the host's blocks — its
//   ramp counted within the block, then held; an empty block in a fade keeps the fade a fade; blend walks dry -> wet
//   for the ramp and holds after it; runs() keeps a link alive for the whole of its fade-out; a prepare's snap stands
//   at an end; and the fade decides as JUCE's approximatelyEqual did, to the sample.

#include <felitronics_test.h>
#include <felitronics/core/BypassFade.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

namespace
{
    /** Samples until a fade from on to off settles, advanced in blocks of `block`. */
    int fadeLength (double rate, int block)
    {
        core::BypassFade f;
        f.prepare (rate);
        f.snapTo (true);
        int samples = 0;
        for (int guard = 0; guard < 100000; ++guard)
        {
            const auto s = f.advance (block, false);
            if (! s.moving)
                return samples;
            samples += s.ramp;
        }
        return -1;
    }
}

int main()
{
    std::printf ("felitronics::core BypassFade tests\n");

    test::group ("settled: no motion, nothing to blend");
    {
        core::BypassFade f;
        f.prepare (48000.0);
        f.snapTo (true);
        const auto s = f.advance (512, true);
        test::ok (! s.moving && core::exactlyEqual (s.from, 1.0f) && core::exactlyEqual (s.to, 1.0f) && s.ramp == 0,
                  "on and asked on: 1 -> 1, not moving");
        f.snapTo (false);
        const auto o = f.advance (512, false);
        test::ok (! o.moving && core::exactlyEqual (o.to, 0.0f), "off and asked off: 0, not moving");
        test::ok (core::BypassFade::runs (s, true) && ! core::BypassFade::runs (o, false), "runs(): in, or not");
    }

    test::group ("a switch is 15 ms of audio, whatever the blocks");
    {
        for (const double rate : { 44100.0, 48000.0, 96000.0 })
        {
            const int want = (int) std::ceil (rate * 0.015);
            for (const int block : { 1, 32, 64, 333, 512, 4096 })
            {
                const int got = fadeLength (rate, block);
                test::ok (std::abs (got - want) <= 1, "rate " + std::to_string ((int) rate) + ", block " + std::to_string (block)
                                                      + ": " + std::to_string (got) + " samples, want ~" + std::to_string (want));
            }
        }
    }

    test::group ("the ramp is counted within the block, then held");
    {
        core::BypassFade f;
        f.prepare (48000.0);
        f.snapTo (false);
        const auto s = f.advance (4096, true);
        test::ok (s.moving && s.ramp >= 719 && s.ramp <= 721 && core::exactlyEqual (s.to, 1.0f), "one big block: the whole fade, 720 samples of it");

        std::vector<float> wet (4096, 1.0f), dry (4096, 0.0f);
        float* w[] { wet.data() };
        const float* d[] { dry.data() };
        core::BypassFade::blend (w, d, 1, 4096, s);
        test::ok (core::exactlyEqual (wet[0], 0.0f), "the first sample is the dry");
        test::ok (wet[(size_t) s.ramp / 2] > 0.45f && wet[(size_t) s.ramp / 2] < 0.55f, "half way through the ramp, half way");
        bool held = true;
        for (size_t i = (size_t) s.ramp; i < wet.size(); ++i)
            held = held && core::exactlyEqual (wet[i], 1.0f);
        test::ok (held, "after the ramp it holds at the wet, for the rest of the block");
        test::ok (! f.advance (4096, true).moving, "and the next block is settled");
    }

    test::group ("an empty block in a fade: it holds, still moving");
    {
        core::BypassFade f;
        f.prepare (48000.0);
        f.snapTo (true);
        const auto first = f.advance (64, false);
        const auto empty = f.advance (0, false);
        test::ok (first.moving && empty.moving && empty.ramp == 0
                      && core::exactlyEqual (empty.from, first.to) && core::exactlyEqual (empty.to, first.to),
                  "under way, no samples: the same place, still moving");
        test::ok (core::BypassFade::runs (empty, false), "a link fading out runs through it");
    }

    test::group ("the fade's own sequence of numbers — the plug-in's, to the bit");
    {
        // Settled within approximatelyEqual's reach is settled: a fade lands exactly on its end (min/max), and the
        // next advance reports no motion; reversed mid-way, it turns around from where it stands.
        core::BypassFade f;
        f.prepare (44100.0);
        f.snapTo (false);
        const auto a = f.advance (100, true);
        const auto b = f.advance (100, false);
        test::ok (core::exactlyEqual (b.from, a.to) && b.to < b.from && b.moving, "a turn mid-way starts where it stood");
        int guard = 0;
        while (f.advance (64, true).moving && ++guard < 1000) {}
        test::ok (! f.advance (64, true).moving && guard < 1000, "and lands, then rests");
    }

    return test::report();
}
