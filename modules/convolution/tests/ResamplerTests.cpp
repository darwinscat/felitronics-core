// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the offline Kaiser IR resampler: unity DC gain, output length, passband
// amplitude preservation, and stopband rejection (anti-aliasing on downsample) >= 55 dB.

#include <felitronics_test.h>
#include <felitronics/convolution/IrResampler.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace felitronics;

static double rmsRange (const std::vector<float>& v, int from, int to)
{
    double s = 0.0; int c = 0;
    for (int i = from; i < to && i < (int) v.size(); ++i) { s += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

int main()
{
    std::printf ("felitronics::convolution resampler tests\n");

    // --- unity DC gain: a constant resamples to the same constant ---
    test::group ("Resampler unity DC gain");
    {
        std::vector<float> in (1000, 0.5f);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        test::ok (! out.empty(), "produced output");
        double maxErr = 0.0;
        for (int i = 50; i < (int) out.size() - 50; ++i) maxErr = std::max (maxErr, (double) std::fabs (out[(std::size_t) i] - 0.5f));
        test::ok (maxErr < 5e-3, "constant preserved (DC gain == 1)");
    }

    // --- output length ~ inLen * ratio ---
    test::group ("Resampler output length");
    {
        std::vector<float> in (4800, 0.0f);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        const int expect = (int) std::llround (4800.0 * 44100.0 / 48000.0);
        test::ok (std::abs ((int) out.size() - expect) <= 1, "outLen ~ inLen*ratio");
    }

    // --- passband amplitude preserved (1 kHz through 48k -> 44.1k) ---
    test::group ("Resampler passband amplitude");
    {
        const int n = 4800; const double f = 1000.0, sr = 48000.0;
        std::vector<float> in (n);
        for (int i = 0; i < n; ++i) in[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * f * i / sr);
        auto out = convolution::resampleIr (in, 48000.0, 44100.0);
        const double amp = rmsRange (out, 100, (int) out.size() - 100) * std::sqrt (2.0);
        test::approx (amp, 1.0, 0.05, "1 kHz amplitude preserved across the SR change");
    }

    // --- anti-alias: a 20 kHz tone is rejected when downsampling 48k -> 24k (Nyquist 12 kHz) ---
    test::group ("Resampler anti-alias stopband (>= 55 dB)");
    {
        const int n = 9600; const double sr = 48000.0, fhi = 20000.0;
        std::vector<float> in (n);
        for (int i = 0; i < n; ++i) in[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * fhi * i / sr);
        auto out = convolution::resampleIr (in, 48000.0, 24000.0);
        const double inAmp = 1.0 / std::sqrt (2.0);                       // input RMS
        const double outR  = rmsRange (out, 100, (int) out.size() - 100);
        const double rejDb = 20.0 * std::log10 ((outR > 1e-12 ? outR : 1e-12) / inAmp);
        test::ok (rejDb < -55.0, "20 kHz tone rejected >= 55 dB on 48k->24k");
    }

    return test::report();
}
