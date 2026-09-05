// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>
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

    // Per-block guard on the integrator state: zap tiny values to exact zero so decaying tails don't
    // sustain subnormals (CPU spikes), AND zap non-finite ones, which is the part that is not an
    // optimisation. Inaudible (< ~-300 dB). A host may also set FTZ/DAZ.
    //
    // WHY POISON AND NOT ONLY DENORMALS. `ic1 = 2*v1 - ic1` is recursive, and a NaN in it reproduces
    // itself for ever: `fabs(NaN) < 1e-15f` is FALSE, so a denormal-only guard steps over exactly the
    // state it needs to clear. Measured on this filter — one +Inf sample into a band-pass, then a
    // healthy 7 kHz tone: 479000 of the next 480000 outputs non-finite, with the per-block flush
    // running throughout. NaN behaves identically. That is the "one bad sample forever" defect, and
    // for an AUDIO-path filter it is the only remedy available: the house rule is that the audio path
    // is not sanitised (a compressor is not a sanitiser), so nothing may gate this filter's input.
    // A DETECTOR path is different and must be gated upstream instead — see the note below.
    //
    // WHAT THIS BUYS AND WHAT IT DOES NOT. It bounds the damage; it does not undo it. Recovery happens
    // at the flush, which the OWNER clocks, so a poisoned filter emits rubbish until the end of the
    // caller's block — 0.3 ms at 32 samples, 84 ms at 4096, 170 ms at 8192 — and then restarts from
    // silence. Restarting is exactly `reset()`: the filter's history is gone, so the recovered stream
    // is NOT the stream a sanitised input would have produced, and this is deliberately not claimed.
    // A detector cannot live with either property, which is why detectors gate their input instead of
    // relying on this.
    //
    // BIT-IDENTICAL ON HEALTHY STREAMS, and that is provable rather than hoped: `core::flushPoison` is
    // a strict superset of `core::flushDenormal` — they agree on every FINITE input, exhaustively
    // swept by class in `core/tests/FlushPoisonTheoryTests.cpp`. The only finite input that reaches a
    // different outcome here is one that overflows the float state itself (sign-alternating |x| ≳
    // 1e38), where the new behaviour is the better one.
    void flushDenormals() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c)
        {
            core::flushPoison (ic1[c]);
            core::flushPoison (ic2[c]);
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
