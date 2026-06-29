// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the polyphase oversampler: upsample interpolation vs the analytic sine
// (correct + alias-free), inter-sample-peak revelation (the limiter's reason to exist), and the
// up→down round-trip == a delayed identity.

#include <felitronics_test.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

int main()
{
    std::printf ("felitronics::oversampling tests\n");
    const double sr = 48000.0;
    const int    tpp = 32;

    // --- 4x upsample interpolation matches the analytic sine (passband, no aliasing) ---
    test::group ("Oversampler 4x upsample == analytic sine");
    {
        const int L = 4, n = 512; const double f = 2000.0, A = 0.8;
        oversampling::PolyphaseOversampler os; test::ok (os.prepare (L, 1, tpp), "prepare 4x");
        std::vector<float> x (n); for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> y ((std::size_t) n * L);
        const float* xi[1] { x.data() }; float* yo[1] { y.data() };
        os.upsample (xi, 1, n, yo);

        const double gdOS = os.filterLatencyOversampled();        // (N-1)/2 OS samples
        double maxErr = 0.0;
        for (int i = 200; i < n * L - 50; ++i)
        {
            const double a = A * std::sin (2.0 * core::kPi * f * (i - gdOS) / (sr * L));
            maxErr = std::max (maxErr, std::fabs (a - y[(std::size_t) i]));
        }
        test::ok (maxErr < 0.02, "OS samples lie on the analytic sine (interp correct, images rejected)");
    }

    // --- reveals an inter-sample peak the baseband samples miss ---
    test::group ("Oversampler reveals inter-sample peak");
    {
        const int L = 4, n = 256; const double f = sr * 0.25, A = 1.0;
        oversampling::PolyphaseOversampler os; os.prepare (L, 1, tpp);
        std::vector<float> x (n); double sampMax = 0.0;
        for (int i = 0; i < n; ++i) { x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr + 0.7)); sampMax = std::max (sampMax, (double) std::fabs (x[i])); }
        std::vector<float> y ((std::size_t) n * L);
        const float* xi[1] { x.data() }; float* yo[1] { y.data() };
        os.upsample (xi, 1, n, yo);
        double osMax = 0.0; for (int i = 100; i < n * L - 100; ++i) osMax = std::max (osMax, (double) std::fabs (y[(std::size_t) i]));
        test::ok (osMax > sampMax + 0.01, "upsampled peak exceeds the sample peak (ISP found)");
        test::ok (osMax <= A * 1.02, "and ~ the true amplitude (no overshoot)");
    }

    // --- up → down round-trip == input delayed by latency (band-limited) ---
    test::group ("Oversampler up->down round-trip == delayed identity");
    {
        const int L = 4, n = 512; const double f = 500.0, A = 0.6;
        oversampling::PolyphaseOversampler os; os.prepare (L, 1, tpp);
        std::vector<float> x (n); for (int i = 0; i < n; ++i) x[i] = (float) (A * std::sin (2.0 * core::kPi * f * i / sr));
        std::vector<float> osb ((std::size_t) n * L), z (n);
        const float* xi[1] { x.data() }; float* ob[1] { osb.data() };
        os.upsample (xi, 1, n, ob);
        const float* obc[1] { osb.data() }; float* zo[1] { z.data() };
        os.downsample (obc, 1, n, zo);

        const int lat = os.latencySamples();
        test::ok (lat > 0, "reports a positive round-trip latency");
        double maxErr = 0.0;
        for (int i = lat + 60; i < n - 10; ++i) maxErr = std::max (maxErr, (double) std::fabs (z[(std::size_t) i] - x[(std::size_t) (i - lat)]));
        test::ok (maxErr < 0.03, "round-trip == input delayed by latency()");
    }

    return test::report();
}
