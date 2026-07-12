// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Adversarial falsification pass for LevelProbe / ReferenceUnity. These are compact pins for the
// external case list plus implementation-boundary attacks not covered by the broader theory suite.

#include <felitronics_test.h>
#include <felitronics/measurement/ReferenceUnity.h>

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::measurement::referenceUnityGain;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
float gainOfMono (std::initializer_list<float> taps, double sampleRate = 48000.0)
{
    std::vector<float> h (taps);
    return referenceUnityGain (h.data(), (int) h.size(), sampleRate);
}
} // namespace

int main()
{
    std::printf ("felitronics::measurement LevelProbe/ReferenceUnity FALSIFY tests\n");

    group ("DeepSeek 2.1-2.5 materialized: delta scale, clamps, silence, pure delay");
    {
        approx (gainOfMono ({ 2.0f }), 0.5, 1.0e-5, "IR=[2] returns gain 0.5");
        approx (gainOfMono ({ 0.02f }), 31.6227766, 1.0e-4, "IR=[0.02] clamps at +30 dB");
        approx (gainOfMono ({ 50.0f }), 0.0316227766, 1.0e-6, "IR=[50] clamps at -30 dB");
        ok (gainOfMono ({ 0.0f, 0.0f, 0.0f, 0.0f }) == 1.0f, "all-zero IR returns exactly 1.0");
        approx (gainOfMono ({ 0.0f, 0.0f, 1.0f }), 1.0, 1.0e-4, "delay-only IR preserves unity");
    }

    group ("DeepSeek 2.6 adapted: tiny non-zero IR is below the documented -60 dB power floor");
    {
        ok (gainOfMono ({ 1.0e-12f }) == 1.0f, "IR=[1e-12] returns exactly 1.0 through the near-silent floor");
        ok (gainOfMono ({ 1.0e-6f }) == 1.0f, "IR at the strict floor boundary is treated as floor/hold");
        approx (gainOfMono ({ 1.01e-3f }), 31.6227766, 1.0e-4, "just above floor still clamps to +30 dB");
    }

    group ("own attack: very high-rate first-second cap must not overrun the max FFT workspace");
    {
        const int n = (1 << 21) + 257;                    // just beyond the implementation's max FFT length
        std::vector<float> h ((std::size_t) n, 0.0f);
        h[0] = 1.0f;
        const float g = referenceUnityGain (h.data(), n, 3.0e6);
        ok (std::isfinite (g), "high-rate long IR returns a finite gain");
        ok (g >= 0.0316227f && g <= 31.62278f, "high-rate long IR remains inside the documented clamp");
    }

    group ("own attack: stereo with one ultra-quiet channel still follows mean-power semantics");
    {
        std::vector<float> loud (512, 0.0f), quiet (512, 0.0f);
        loud[0] = 1.0f;
        quiet[0] = 1.0e-4f;                                // above zero, far below loud; mean ~= 0.5
        const float* ch[2] { loud.data(), quiet.data() };
        approx (referenceUnityGain (ch, 2, 512, 48000.0), std::sqrt (2.0), 3.0e-4,
                "one near-silent stereo channel yields the mean-power sqrt(2) gain");
    }

    return felitronics::test::report();
}
