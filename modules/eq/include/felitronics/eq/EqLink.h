// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>
#include <felitronics/eq/MatchedBiquad.h>

#include <algorithm>

namespace felitronics::eq
{

/** One EQ link of the chain — the console grammar, spoken in matched biquads. Moved here from OrbitAmp
    (src/core/EqLink.h), its arithmetic unchanged: the amps' consoles and the bass amp's tone, the plug-in's and the
    web's.

    HPF and LPF with a real slope choice (6/12/18/24/48 dB/oct, Butterworth cascades), two
    shelves with free corners, two tone bells, and a third NARROW bell that switches in — the
    surgical slot search→treat will land in. No output level: the block's level used to ride
    here, and with no hand on any face left to reach it, it went on playing whatever an old
    session had left in it. No presence: with free shelf corners a second fixed high shelf is a
    duplicate. No linear phase either: it buys latency and pre-ring for a benefit a serial chain
    into a nonlinearity cannot hear — analog EQs are minimum-phase, and so is this.

    JUCE-free, allocation-free after prepare; magnitudeDb reads the same coefficients process
    runs, so the drawing can never drift from the sound. filterMagnitudeDb answers for the cut
    filters alone — the face draws them as their own walls. */
class EqLink
{
public:
    static constexpr int maxChannels = 2;

    struct Settings
    {
        bool   hpfOn    = false;
        double hpfHz    = 80.0;
        int    hpfSlope = 12;       // dB/oct: 6, 12, 18, 24 or 48

        double loDb = 0.0, loHz = 100.0;

        double b1Db = 0.0, b1Hz = 400.0,  b1Q = 1.0;
        double b2Db = 0.0, b2Hz = 2500.0, b2Q = 1.0;

        bool   b3On = false;
        double b3Db = 0.0, b3Hz = 3000.0, b3Q = 8.0;

        double hiDb = 0.0, hiHz = 8000.0;

        bool   lpfOn    = false;
        double lpfHz    = 10000.0;
        int    lpfSlope = 12;

        /** Member for member, as a defaulted == compares them (a NaN equals nothing) — spelt out with core::exactlyEqual so
            a consumer's -Wfloat-equal stays quiet. */
        bool operator== (const Settings& o) const noexcept
        {
            using core::exactlyEqual;
            return hpfOn == o.hpfOn && exactlyEqual (hpfHz, o.hpfHz) && hpfSlope == o.hpfSlope
                && exactlyEqual (loDb, o.loDb) && exactlyEqual (loHz, o.loHz)
                && exactlyEqual (b1Db, o.b1Db) && exactlyEqual (b1Hz, o.b1Hz) && exactlyEqual (b1Q, o.b1Q)
                && exactlyEqual (b2Db, o.b2Db) && exactlyEqual (b2Hz, o.b2Hz) && exactlyEqual (b2Q, o.b2Q)
                && b3On == o.b3On && exactlyEqual (b3Db, o.b3Db) && exactlyEqual (b3Hz, o.b3Hz) && exactlyEqual (b3Q, o.b3Q)
                && exactlyEqual (hiDb, o.hiDb) && exactlyEqual (hiHz, o.hiHz)
                && lpfOn == o.lpfOn && exactlyEqual (lpfHz, o.lpfHz) && lpfSlope == o.lpfSlope;
        }
    };

    void prepare (double sampleRate, int numChannels) noexcept
    {
        fs = sampleRate > 0.0 ? sampleRate : 48000.0;
        ch = std::clamp (numChannels, 1, maxChannels);
        design (current, true);
        reset();
    }

    void reset() noexcept
    {
        for (auto& band : state)
            for (auto& biquad : band)
                biquad.reset();
    }

    void setSettings (const Settings& s) noexcept
    {
        if (! (s == current))
            design (s, false);
    }

    const Settings& getSettings() const noexcept { return current; }

    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int n = std::min (numChannels, ch);

        for (int b = 0; b < numSections; ++b)
        {
            if (! active[b])
                continue;

            for (int c = 0; c < n; ++c)
            {
                auto& biquad = state[b][(size_t) c];
                auto* data   = channels[c];

                for (int i = 0; i < numSamples; ++i)
                    data[i] = biquad.processSample (data[i]);

                biquad.flushDenormals();
            }
        }
    }

    /** The full response in dB at a frequency — filters, shelves and bells. */
    double magnitudeDb (double freqHz) const noexcept
    {
        return sumDb (freqHz, 0, numSections);
    }

    /** The cut filters alone — what the face draws as the red walls. */
    double filterMagnitudeDb (double freqHz) const noexcept
    {
        return sumDb (freqHz, 0, filterSections);
    }

    /** The tone sections alone — shelves and bells, the bright curve's own voice. */
    double toneMagnitudeDb (double freqHz) const noexcept
    {
        return sumDb (freqHz, filterSections, numSections);
    }

private:
    // Sections: HPF (up to 4), LPF (up to 4), then LO / B1 / B2 / B3 / HI — cut filters first so
    // one index split serves both magnitude views.
    static constexpr int hpfSlots      = 4;
    static constexpr int lpfSlots      = 4;
    static constexpr int filterSections = hpfSlots + lpfSlots;
    static constexpr int numSections   = filterSections + 5;

    double sumDb (double freqHz, int from, int to) const noexcept
    {
        const double w = 2.0 * 3.14159265358979323846 * std::clamp (freqHz, 1.0, 0.499 * fs) / fs;

        double db = 0.0;
        for (int b = from; b < to; ++b)
            if (active[b])
                db += coeffs[b].magnitudeDb (w);

        return db;
    }

    /** Fills a run of slots with a Butterworth cascade of the asked slope. The Q ladder is the
        textbook one; the odd slopes lead with a first-order section. */
    void designCascade (int firstSlot, bool on, bool isHighpass, double hz, int slopeDb) noexcept
    {
        namespace m = matched;

        const double nyquist = 0.49 * fs;
        const double f = std::clamp (hz, 10.0, nyquist);

        BiquadCoeffs sections[4];
        int count = 0;

        const auto fo = [&] { return isHighpass ? m::highpass1 (f, fs) : m::lowpass1 (f, fs); };
        const auto bq = [&] (double q) { return isHighpass ? m::highpass (f, fs, q) : m::lowpass (f, fs, q); };

        switch (slopeDb)
        {
            case 6:  sections[count++] = fo(); break;
            default:
            case 12: sections[count++] = bq (0.70710678); break;
            case 18: sections[count++] = fo();
                     sections[count++] = bq (1.0); break;
            case 24: sections[count++] = bq (0.54119610);
                     sections[count++] = bq (1.30656296); break;
            case 48: sections[count++] = bq (0.50979558);
                     sections[count++] = bq (0.60134489);
                     sections[count++] = bq (0.89997622);
                     sections[count++] = bq (2.56291545); break;
        }

        for (int i = 0; i < 4; ++i)
            set (firstSlot + i, on && i < count,
                 i < count ? sections[i] : BiquadCoeffs {});
    }

    void design (const Settings& s, bool force) noexcept
    {
        namespace m = matched;

        const double nyquist = 0.49 * fs;
        const auto   clampHz = [nyquist] (double f) { return std::clamp (f, 10.0, nyquist); };

        designCascade (0,        s.hpfOn, true,  s.hpfHz, s.hpfSlope);
        designCascade (hpfSlots, s.lpfOn, false, s.lpfHz, s.lpfSlope);

        // A band at exactly 0 dB is a bypass — skipping it keeps a flat link bit-transparent.
        int b = filterSections;
        const auto on = [] (double db) { return ! core::exactlyEqual (db, 0.0); };
        set (b++, on (s.loDb),            m::lowShelfDb  (clampHz (s.loHz), fs, s.loDb));
        set (b++, on (s.b1Db),            m::peakingDb   (clampHz (s.b1Hz), fs, s.b1Q, s.b1Db));
        set (b++, on (s.b2Db),            m::peakingDb   (clampHz (s.b2Hz), fs, s.b2Q, s.b2Db));
        set (b++, s.b3On && on (s.b3Db),  m::peakingDb   (clampHz (s.b3Hz), fs, s.b3Q, s.b3Db));
        set (b++, on (s.hiDb),            m::highShelfDb (clampHz (s.hiHz), fs, s.hiDb));

        current = s;

        if (force)
            reset();
    }

    void set (int band, bool isActive, const BiquadCoeffs& c) noexcept
    {
        active[band] = isActive;
        coeffs[band] = c;

        for (auto& biquad : state[band])
            biquad.setCoeffs (c);
    }

    bool                          active[numSections] {};
    BiquadCoeffs coeffs[numSections] {};
    Biquad       state[numSections][maxChannels] {};

    Settings current;
    double   fs = 48000.0;
    int      ch = 2;
};

} // namespace felitronics::eq
