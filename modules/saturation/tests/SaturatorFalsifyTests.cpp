// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Compact adversarial pins for Saturator hardening: sanitize-then-clamp order, stateful chunking,
// and failed-prepare no-op behavior.

#include <felitronics_test.h>
#include <felitronics/saturation/Saturator.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using Sat = felitronics::saturation::Saturator;
using Shape = felitronics::saturation::WaveShaper::Shape;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
bool bitEq (float a, float b)
{
    return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b);
}

bool equalBits (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (! bitEq (a[i], b[i])) return false;
    return true;
}
} // namespace

int main()
{
    std::printf ("felitronics::saturation Saturator FALSIFY tests\n");

    group ("DeepSeek 6.4 materialized: sanitize-then-clamp gate order, including -0 preservation");
    {
        Sat::Params p;
        p.mix = 0.0f;                                      // output the sanitized dry path
        p.outputDb = 0.0f;
        std::vector<float> x {
            std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity(),
            2.0e6f,
            -3.0e6f,
            5.0e5f,
            0.0f,
            -0.0f
        };
        const std::vector<float> want {
            0.0f, 0.0f, 0.0f, 1.0e6f, -1.0e6f, 5.0e5f, 0.0f, -0.0f
        };
        Sat s;
        ok (s.prepare (48000.0, 16, 1, 1), "prepare os=1");
        s.setParams (p);
        float* io[1] { x.data() };
        s.process (io, 1, (int) x.size());
        ok (equalBits (x, want), "gate maps [NaN,Inf,-Inf,2e6,-3e6,5e5,0,-0] by sanitize-then-clamp");
    }

    group ("DeepSeek 6.5 materialized: stateful chunking matches single-call bit-for-bit");
    {
        Sat::Params p;
        p.shape = Shape::Asym;
        p.driveDb = 9.0f;
        p.bias = 0.35f;
        p.mix = 0.65f;
        p.autoComp = 0.8f;
        p.dcBlockHz = 17.0f;

        std::vector<float> one (4096), split (4096);
        for (int i = 0; i < 4096; ++i)
            one[(std::size_t) i] = split[(std::size_t) i] =
                0.35f * (float) std::sin (0.011 * i) + 0.11f * (float) std::sin (0.071 * i);

        Sat a, b;
        ok (a.prepare (48000.0, 257, 1, 4) && b.prepare (48000.0, 257, 1, 4), "prepare stateful saturators");
        a.setParams (p); b.setParams (p);
        { float* io[1] { one.data() }; a.process (io, 1, (int) one.size()); }
        int pos = 0;
        const int chunks[] { 1, 64, 3, 257, 19, 128, 2, 511 };
        int ci = 0;
        while (pos < (int) split.size())
        {
            if ((ci % 5) == 0)
            {
                float* z[1] { split.data() + pos };
                b.process (z, 1, 0);
            }
            const int n = std::min (chunks[ci++ % 8], (int) split.size() - pos);
            float* io[1] { split.data() + pos };
            b.process (io, 1, n);
            pos += n;
        }
        ok (equalBits (one, split), "Asym/DC/dry-delay state is associative across hostile chunks");
    }

    group ("DeepSeek 6.6 adapted: in-range large finite values remain finite and gate-transparent");
    {
        Sat::Params p;
        p.mix = 0.0f;                                      // inspect the gate directly through dry delay 0
        std::vector<float> x { -999999.0f, -1.0f, -0.0f, 0.0f, 1.0f, 999999.0f };
        const auto want = x;
        Sat s;
        ok (s.prepare (48000.0, 16, 1, 1), "prepare os=1");
        s.setParams (p);
        float* io[1] { x.data() };
        s.process (io, 1, (int) x.size());
        ok (equalBits (x, want), "finite in-range values pass through the gate bit-identically");
    }

    group ("own attack: failed prepare leaves process as a no-op");
    {
        Sat s;
        ok (! s.prepare (-48000.0, 0, 2, 4), "invalid prepare fails");
        std::vector<float> a { 1.0f, -2.0f, 3.0f };
        const auto before = a;
        float* io[1] { a.data() };
        s.process (io, 1, (int) a.size());
        ok (equalBits (a, before), "process after failed prepare leaves caller buffer untouched");
    }

    return felitronics::test::report();
}
