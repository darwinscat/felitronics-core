// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::eq::BassCuts and CutBiquad — the bass amp's own cuts (moved from OrbitAmp).
//
//   CutBiquad is a 2nd-order Butterworth: -3.01 dB at its corner, 12 dB/oct beyond it, unity in its pass band;
//   BassCuts settled and off is a wire, bit for bit; on, it is the biquad; a switch is a 15 ms crossfade and a
//   frequency move a glide that lands exactly; the output does not depend on how the host cuts its blocks (on its
//   32-sample cadence, bit for bit; off it, within an ulp where a fade has just landed); a NaN
//   that comes in leaves no NaN in the recursion; and nothing allocates.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/eq/BassCuts.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

using namespace felitronics;

namespace
{
    constexpr double fs = 48000.0;

    /** The steady gain, dB, of a cut on a sine at `hz` — after a second to settle, over the next one. */
    double gainDb (bool high, double cornerHz, double hz)
    {
        eq::CutBiquad b;
        b.design (high, cornerHz, fs);
        double in = 0.0, out = 0.0;
        for (int i = 0; i < (int) (2.0 * fs); ++i)
        {
            const float x = (float) std::sin (2.0 * core::kPi * hz * (double) i / fs);
            const float y = b.process (x);
            if (i >= (int) fs)
            {
                in  += (double) x * x;
                out += (double) y * y;
            }
        }
        return 10.0 * std::log10 (out / in);
    }

    std::vector<float> phrase (int n)
    {
        std::vector<float> x ((size_t) n);
        for (int i = 0; i < n; ++i)
        {
            const double t = (double) i / fs;
            x[(size_t) i] = (float) (0.3 * std::sin (2.0 * core::kPi * 41.2 * t) + 0.2 * std::sin (2.0 * core::kPi * 3000.0 * t)
                                     + 0.1 * std::sin (2.0 * core::kPi * 9.0 * t));
        }
        return x;
    }

    /** `x` through a stage, in host blocks of `block`, with `settings(i)` the settings of the block starting at i. */
    template <typename Settings>
    std::vector<float> run (std::vector<float> x, int block, Settings settings)
    {
        eq::BassCuts c;
        c.prepare (fs);
        for (int i = 0; i < (int) x.size(); i += block)
        {
            float* ch[] { x.data() + i };
            c.process (ch, 1, std::min (block, (int) x.size() - i), settings (i));
        }
        return x;
    }
}

int main()
{
    std::printf ("felitronics::eq BassCuts tests\n");

    test::group ("CutBiquad: a 2nd-order Butterworth");
    {
        test::approx (gainDb (true, 100.0, 100.0), -3.0103, 0.02, "high pass: -3 dB at its corner");
        test::approx (gainDb (true, 100.0, 50.0), -12.30, 0.1, "high pass: an octave under it, -12.3 dB");
        test::approx (gainDb (true, 100.0, 2000.0), 0.0, 0.01, "high pass: unity in its pass band");
        test::approx (gainDb (false, 4500.0, 4500.0), -3.0103, 0.05, "low pass: -3 dB at its corner");
        test::approx (gainDb (false, 4500.0, 200.0), 0.0, 0.01, "low pass: unity in its pass band");
    }

    const auto x = phrase ((int) fs);

    test::group ("settled and off: a wire, bit for bit");
    {
        const auto y = run (x, 512, [] (int) { return eq::BassCuts::Settings {}; });
        test::ok (std::memcmp (x.data(), y.data(), sizeof (float) * x.size()) == 0, "both cuts off: the input, unchanged");
    }

    test::group ("settled and on: the biquads, plainly");
    {
        eq::BassCuts::Settings s;
        s.hpfOn = true;
        s.hpfHz = 32.0;
        s.lpfOn = true;
        s.lpfHz = 4500.0;
        const auto y = run (x, 512, [&] (int) { return s; });
        eq::CutBiquad h, l;
        h.design (true, 32.0, fs);
        l.design (false, 4500.0, fs);
        bool same = true;
        for (size_t i = 0; i < x.size(); ++i)
            same = same && core::exactlyEqual (y[i], l.process (h.process (x[i])));
        test::ok (same, "the HPF then the LPF, sample for sample");
    }

    test::group ("the host's blocks do not matter");
    {
        const auto moves = [] (int i)
        {
            eq::BassCuts::Settings s;
            s.hpfOn = i >= 4800;
            s.hpfHz = i < 20000 ? 22.0 : 32.0;
            s.lpfOn = i < 30000;
            s.lpfHz = i < 12000 ? 4500.0 : 9500.0;
            return s;
        };
        // The settings change on boundaries every block size here shares (multiples of 96).
        const auto a = run (x, 96, [&] (int i) { return moves (i - i % 96); });
        const auto b = run (x, 32, [&] (int i) { return moves (i - i % 96); });
        const auto c = run (x, 1, [&] (int i) { return moves (i - i % 96); });
        test::ok (std::memcmp (a.data(), b.data(), sizeof (float) * a.size()) == 0, "blocks of 96 and 32: the same samples");
        // Sliced off the 32-sample cadence, the samples are the same but for the rest of a cadence step after a fade
        // lands: there a block that began the step fading keeps writing d + (y - d) * 1, which may round an ulp from the y a
        // fresh step writes. The plug-in's own arithmetic, moved unchanged — and held here to what it is.
        double worst = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            worst = std::max (worst, std::abs ((double) a[i] - (double) c[i]));
        test::ok (worst <= 1.0e-7, "blocks of 96 and 1: within an ulp of the samples (" + std::to_string (worst) + ")");
    }

    test::group ("a switch fades over 15 ms; a frequency glides and lands");
    {
        eq::BassCuts c;
        c.prepare (fs);
        eq::BassCuts::Settings s;
        std::vector<float> buf (64, 0.0f);
        float* ch[] { buf.data() };
        c.process (ch, 1, 64, s);
        s.hpfOn = true;
        int samples = 0;
        while (c.gainOf (0) < 1.0f && samples < 48000)
        {
            c.process (ch, 1, 64, s);
            samples += 64;
        }
        test::ok (samples >= 704 && samples <= 768, "on: 15 ms (720 samples) of crossfade, " + std::to_string (samples));

        s.hpfHz = 40.0;
        int glide = 0;
        while (c.hzOf (0) < 40.0 && glide < 48000)
        {
            c.process (ch, 1, 64, s);
            glide += 64;
        }
        test::ok (core::exactlyEqual (c.hzOf (0), 40.0), "the glide lands on the frequency asked, exactly");
        test::ok (glide > 1000 && glide < 3000, "...in ~30 ms to a percent and a little more to the end, " + std::to_string (glide));
    }

    test::group ("a NaN in, no NaN kept");
    {
        eq::BassCuts c;
        c.prepare (fs);
        eq::BassCuts::Settings s;
        s.hpfOn = true;
        s.lpfOn = true;
        std::vector<float> buf (64, 0.1f);
        buf[10] = std::numeric_limits<float>::quiet_NaN();
        float* ch[] { buf.data() };
        c.process (ch, 1, 64, s);
        std::vector<float> after (64, 0.1f);
        float* ch2[] { after.data() };
        c.process (ch2, 1, 64, s);
        bool finite = true;
        for (const float v : after)
            finite = finite && std::isfinite (v);
        test::ok (finite, "the block after a NaN is finite");
    }

    return test::report();
}
