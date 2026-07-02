// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/EqTypes.h>

#include <cmath>

namespace felitronics::eq
{

//==============================================================================
// teq::Svf — Cytomic / Zavalishin TPT state-variable filter (Andrew Simper, "Solving the
// continuous SVF equations using trapezoidal integration"). Zero-delay feedback → stays clean
// under fast fc / Q modulation, which is exactly what the search-mode sweep needs. One filter,
// per-channel integrator state. Coefficients are BLT-prewarped (tan), so it's accurate at the
// cutoff but — unlike the matched biquad — not Nyquist-perfect; that's the right trade for a
// transient swept band, while static treatment bands use the matched design.
class Svf
{
public:
    static constexpr int kMaxChannels = felitronics::core::kMaxChannels;   // single source of truth: teq/EqTypes.h

    void prepare (double sampleRate, int numChannels) noexcept
    {
        fs = sampleRate;
        ch = numChannels < 1 ? 1 : (numChannels > kMaxChannels ? kMaxChannels : numChannels);
        reset();
    }

    void reset() noexcept { for (int c = 0; c < kMaxChannels; ++c) { ic1[c] = 0.0f; ic2[c] = 0.0f; } }

    // Precondition: caller passes a finite freq; Q / gainDb are sanitised here defensively.
    void setParams (FilterType type, double freq, double Q, double gainDb) noexcept
    {
        if (Q < 1e-3) Q = 1e-3;                                   // guard 1/Q
        if (gainDb >  60.0) gainDb =  60.0;
        if (gainDb < -60.0) gainDb = -60.0;

        const double A = std::pow (10.0, gainDb / 40.0);
        double f = freq;
        if (f < 1.0) f = 1.0;
        if (f > 0.49 * fs) f = 0.49 * fs;
        double g = std::tan (kPi * f / fs);
        double k = 1.0 / Q;
        constexpr double kButterworth = 1.4142135623730951;      // sqrt(2): 2-pole Butterworth damping (shelves ignore Q)

        switch (type)
        {
            case FilterType::Bell:      k = 1.0 / (Q * A);                    m0 = 1.0;   m1 = k * (A * A - 1.0); m2 = 0.0;           break;
            case FilterType::LowShelf:  k = kButterworth; g /= std::sqrt (A); m0 = 1.0;   m1 = k * (A - 1.0);     m2 = (A * A - 1.0); break;
            case FilterType::HighShelf: k = kButterworth; g *= std::sqrt (A); m0 = A * A; m1 = k * (1.0 - A) * A; m2 = (1.0 - A * A); break;
            case FilterType::LowPass:                                         m0 = 0.0;   m1 = 0.0;               m2 = 1.0;           break;
            case FilterType::HighPass:                                        m0 = 1.0;   m1 = -k;                m2 = -1.0;          break;
            case FilterType::BandPass:                                        m0 = 0.0;   m1 = k;                 m2 = 0.0;           break;  // unity gain at centre
            case FilterType::AllPass:                                        m0 = 1.0;   m1 = -2.0 * k;          m2 = 0.0;           break;  // 2nd-order allpass: v0 - 2k·v1 (flat |H|=1)
            case FilterType::Notch:                                          m0 = 1.0;   m1 = -k;                m2 = 0.0;           break;  // notch = low + high = v0 - k·v1: exact null at the prewarped fc, unity at DC/Nyquist
            case FilterType::Tilt:                                           m0 = 1.0;   m1 = 0.0;               m2 = 0.0;           break;  // NOT sweepable (EqBand always runs the matched two-shelf tilt — a one-SVF tilt cannot hold a unity pivot past ~6 dB); pass-through kept defensively
        }

        a1 = 1.0 / (1.0 + g * (g + k));
        a2 = g * a1;
        a3 = g * a2;
    }

    inline float processSample (int c, float in) noexcept
    {
        const double v0 = in;
        const double v3 = v0 - ic2[c];
        const double v1 = a1 * ic1[c] + a2 * v3;
        const double v2 = ic2[c] + a2 * ic1[c] + a3 * v3;
        ic1[c] = (float) (2.0 * v1 - ic1[c]);
        ic2[c] = (float) (2.0 * v2 - ic2[c]);
        return (float) (m0 * v0 + m1 * v1 + m2 * v2);
    }

    // Per-block denormal guard: zap tiny integrator state to exact zero so decaying tails
    // don't sustain subnormals (CPU spikes). Inaudible (< ~-300 dB). Host should also set FTZ/DAZ.
    void flushDenormals() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c)
        {
            if (std::fabs (ic1[c]) < 1e-15f) ic1[c] = 0.0f;
            if (std::fabs (ic2[c]) < 1e-15f) ic2[c] = 0.0f;
        }
    }

private:
    double fs = 44100.0;
    int    ch = 2;
    double a1 = 0.0, a2 = 0.0, a3 = 0.0;
    double m0 = 1.0, m1 = 0.0, m2 = 0.0;
    float  ic1[kMaxChannels] {}, ic2[kMaxChannels] {};
};

} // namespace felitronics::eq
