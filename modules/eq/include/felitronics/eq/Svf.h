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

    // Precondition: the caller passes FINITE parameters. Q and gainDb are BOUNDED here, which is not
    // the same as sanitised and the difference matters: `Q < 1e-3` and both gain comparisons are false
    // for a NaN, so a NaN Q walks straight through into `pow`, `k`, `a1..a3` and `m0..m2`, and every
    // output is non-finite from then on — 4096 of 4096, measured. THAT poison is in the COEFFICIENTS,
    // not the state, so no amount of flushing reaches it; `reset()` does not either. It is the one
    // exception to "a bad sample is no longer permanent", and it is named here rather than left to be
    // discovered. Every consumer in this repository sanitises before calling (EqBand, DynamicEqBand,
    // DeEsser, LaneDynamics, MultibandSplitter, MonoBass, PowerAmpStage all do), so there is no live
    // defect — but the primitive's own API does not, and the comment used to say otherwise.
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
    // WHAT DOES AND DOES NOT REACH THIS STATE FROM A FINITE INPUT, stated carefully because the
    // convenient version of it is false. A single 1e20 does not poison any configuration tried (0 of
    // 675: 9 filter types × 5 frequencies × 5 Q values × 3 gains). A single 3e38 DOES, in 195 of those
    // same 675 — the poison threshold is 3.400e38 for a 7 kHz band-pass at Q 2 against an FLT_MAX of
    // 3.403e38, and it falls to 2.05e38 at 20 kHz, because `ic1 = 2*v1 - ic1` doubles the magnitude of
    // the state. Sign alternation is not required. So "only inputs above 1e38 differ" would have been
    // a statement about one fixture sitting 0.1 % under an edge, not about the filter.
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
    // caller's block and then restarts from silence. The window is therefore anywhere from 0 to one
    // block short of one, DEPENDING ON WHERE IN THE BLOCK THE SAMPLE LANDS — swept over every offset
    // at 48 kHz: 0..0.65 ms at 32 samples, 0..10.65 at 512, 0..85.31 at 4096, 0..170.65 at 8192. The
    // worst case is the one to design against, and quoting a single figure for it without saying
    // where the poison fell would be a real measurement of an unstated fixture.
    // Restarting is exactly `reset()`: the filter's history is gone, so the recovered stream is NOT
    // the stream a sanitised input would have produced, and this is deliberately not claimed. A
    // detector cannot live with either property, which is why detectors gate their input instead of
    // relying on this.
    //
    // BIT-IDENTICAL WHENEVER THE STATE IS FINITE, which is the claim that is actually provable: the
    // branch above reduces to the old per-word denormal flush the moment both words are finite, and
    // `core::flushPoison`'s agreement with `core::flushDenormal` on every finite value is swept
    // exhaustively by class in `core/tests/FlushPoisonTheoryTests.cpp`. Note the claim is about the
    // STATE and not about the input: a perfectly finite input can reach non-finite state by
    // overflowing on the way (a sign-alternating 1e37 here; a mere 1e20 in an Rms envelope, where the
    // state is a square). Saying "only inputs above 1e38 differ" would have been wrong for that
    // reason, and is not said.
    // ATOMIC PER CHANNEL, and that is not tidiness. The two integrators are ONE state joined by the
    // recursion, and clearing them word by word "heals" a filter into something that is not a reset
    // filter: it carries the surviving word forward. Measured on the sibling biquad, where the state
    // is reachable — a finite 1e37 sine at Q 40 puts `z1` at +Inf while `z2` is still a finite
    // -3.27e38; a per-word flush zeroes `z1`, and two samples of SILENCE later the filter emits
    // -3.27e38 and then -Inf. It re-poisons itself from the half that was left. So: if either word of
    // a channel is non-finite, BOTH go, which is the only thing that makes `reset()` an honest
    // description of the result.
    // FOR THIS FILTER SPECIFICALLY THAT IS A STRUCTURAL GUARANTEE RATHER THAN A DEMONSTRATED FIX, and
    // the difference is worth stating. Both integrators are read by both update lines, so whichever
    // one goes first takes the other with it on the next sample at the latest, and a block-edge flush
    // never sees exactly one bad word. Searched for a counterexample — 6 filter types × 5 frequencies
    // × 6 Q values × 6 amplitudes from 1e35 to 3e38 × 3 waveforms, this build against a word-by-word
    // one — and found none in 25920 outputs. The shape is kept anyway: the sibling biquad needs it
    // for real, and one contract with one shape beats two.
    void flushDenormals() noexcept
    {
        for (int c = 0; c < kMaxChannels; ++c)
        {
            if (! std::isfinite (ic1[c]) || ! std::isfinite (ic2[c])) { ic1[c] = 0.0f; ic2[c] = 0.0f; }
            else { core::flushDenormal (ic1[c]); core::flushDenormal (ic2[c]); }
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
