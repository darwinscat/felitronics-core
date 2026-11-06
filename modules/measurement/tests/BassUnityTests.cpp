// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::measurement::bassReferenceUnityGain — the bass's unity (moved from OrbitAmp's cabinet).
//
//   an impulse passes at unity; a level is undone exactly; an IR that keeps the bass band and rolls the top off is
//   unity on the bass's reference where the guitar's calls it quieter; the family's clamp and its near-silent floor;
//   degenerate input leaves the IR alone.

#include <felitronics_test.h>
#include <felitronics/measurement/ReferenceUnity.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

int main()
{
    std::printf ("felitronics::measurement bass unity tests\n");
    constexpr double fs = 48000.0;

    test::group ("an impulse, and a level");
    {
        std::vector<float> ir (4096, 0.0f);
        ir[0] = 1.0f;
        test::approx (measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs), 1.0, 1.0e-6, "an impulse passes at unity");
        ir[0] = 0.5f;
        test::approx (measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs), 2.0, 1.0e-5, "half of one: twice the gain");
    }

    test::group ("the bass band decides, not the top");
    {
        // A one-pole low pass at 800 Hz: the bass band whole, the top rolled off.
        std::vector<float> ir (48000, 0.0f);
        const double a = std::exp (-2.0 * 3.141592653589793 * 800.0 / fs);
        double y = 0.0;
        for (size_t i = 0; i < ir.size(); ++i)
        {
            y = (i == 0 ? 1.0 - a : 0.0) + a * y;
            ir[i] = (float) y;
        }
        const float bass = measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs);
        const float guitar = measurement::referenceUnityGain (ir.data(), (int) ir.size(), fs);
        test::approx (20.0 * std::log10 (bass), 0.07, 0.1, "the bass's reference: within a tenth of a dB of unity");
        test::ok (guitar > bass * 1.05f, "the guitar's reference calls it quieter and lifts it more");
    }

    test::group ("the clamp, the floor, degenerate input");
    {
        std::vector<float> ir (4096, 0.0f);
        ir[0] = 0.01f;   // -40 dB: asks for +40, clamped to +30
        test::approx (measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs), std::pow (10.0, 30.0 / 20.0), 1.0e-3, "clamped at +30 dB");
        ir[0] = 1.0e-4f; // -80 dB: under the floor, left alone
        test::approx (measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs), 1.0, 0.0, "near-silent: 1");
        test::approx (measurement::bassReferenceUnityGain (nullptr, 10, fs), 1.0, 0.0, "no taps: 1");
        test::approx (measurement::bassReferenceUnityGain (ir.data(), 0, fs), 1.0, 0.0, "no length: 1");
        test::approx (measurement::bassReferenceUnityGain (ir.data(), 10, 0.0), 1.0, 0.0, "no rate: 1");
    }

    return test::report();
}
