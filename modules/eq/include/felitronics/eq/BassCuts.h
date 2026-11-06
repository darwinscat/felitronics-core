// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace felitronics::eq
{

/** One 2nd-order Butterworth section, high- or low-pass (RBJ), transposed direct form II in double — the bass's
    cuts and its COMP's key. Moved here from OrbitAmp (src/core/BassCuts.h) with BassCuts, its arithmetic unchanged.
    No denormal flush of its own: in the plug-in it runs under the chain's ScopedNoDenormals, and it gives the same
    samples however a host cuts its blocks. */
struct CutBiquad
{
    double b0 = 1.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;
    double z1 = 0.0, z2 = 0.0;

    void design (bool high, double hz, double sampleRate) noexcept
    {
        constexpr double pi = 3.14159265358979323846;
        const double f  = std::clamp (hz, 1.0, 0.45 * sampleRate);
        const double w0 = 2.0 * pi * f / sampleRate;
        const double cw = std::cos (w0), alpha = std::sin (w0) / (2.0 * 0.7071067811865476);
        const double a0 = 1.0 + alpha;
        const double k  = high ? (1.0 + cw) * 0.5 : (1.0 - cw) * 0.5;
        b0 = k / a0;
        b1 = (high ? -2.0 * k : 2.0 * k) / a0;
        b2 = k / a0;
        a1 = -2.0 * cw / a0;
        a2 = (1.0 - alpha) / a0;
    }

    float process (float x) noexcept
    {
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return (float) y;
    }

    void reset() noexcept { z1 = z2 = 0.0; }
};

/** THE BASS AMP'S OWN CUTS (moved here from OrbitAmp, its arithmetic unchanged): an HPF and an LPF, 12 dB/oct each,
    every change of them a transition rather than a step — the stage after DI/CAB, and the gentle inner cuts across
    the chain.

      - ON / OFF is a crossfade of the cut against its input over BypassFade's 15 ms; once a cut has faded out
        its filter's memory is cleared, so turning it on again never replays a note from before.
      - A FREQUENCY moves by gliding — log-evenly, ~30 ms to within a percent, the coefficients redrawn every
        32 samples — so 5 <-> 4, FINGERS <-> SLAP or a pick from the menu is a sweep, not a jump.
      - Settled and on, the cut runs plainly; settled and off, it is not run at all: the stage is a wire.
      - A NaN that comes in heals at the end of its cadence step: the filter's memory is cleared, never kept.

    The console's EqLink is not this. Up to two channels; real-time safe, no allocation. */
class BassCuts
{
public:
    static constexpr int   maxChannels = 2;
    static constexpr int   cadence     = 32;     // samples between redrawn coefficients while a frequency glides
    static constexpr float fadeMs      = 15.0f;  // core::BypassFade's
    static constexpr double glideMs    = 30.0;   // to within about a percent

    struct Settings
    {
        bool   hpfOn = false;
        double hpfHz = 22.0;
        bool   lpfOn = false;
        double lpfHz = 4500.0;
    };

    void prepare (double sampleRate) noexcept
    {
        fs        = sampleRate > 0.0 ? sampleRate : 48000.0;
        fadeStep  = (float) (1000.0 / (fs * (double) fadeMs));
        glideKeep = std::exp (-(double) cadence / (glideMs * 0.001 * fs / 4.6));   // e^-4.6 ~ 1 % left after glideMs
        settled   = false;
    }

    /** Stand where `s` says at once, every memory cleared — a prepare's first block, not a move. */
    void snap (const Settings& s) noexcept
    {
        cuts[0].snap (true,  s.hpfOn, s.hpfHz, fs);
        cuts[1].snap (false, s.lpfOn, s.lpfHz, fs);
        settled = true;
    }

    void process (float* const* channels, int numChannels, int numSamples, const Settings& s) noexcept
    {
        if (! settled)
            snap (s);

        const int nc = std::min (numChannels, maxChannels);
        cuts[0].run (channels, nc, numSamples, s.hpfOn, s.hpfHz, fs, fadeStep, glideKeep);
        cuts[1].run (channels, nc, numSamples, s.lpfOn, s.lpfHz, fs, fadeStep, glideKeep);
    }

    /** Where a cut stands now — for a gate: its fade (0 off .. 1 on) and its frequency. */
    float  gainOf (int cut) const noexcept { return cuts[(size_t) cut].g; }
    double hzOf (int cut) const noexcept   { return cuts[(size_t) cut].hzNow; }

private:
    struct Cut
    {
        bool   high = true;
        float  g = 0.0f;           // the crossfade: 0 the input, 1 the cut
        double hzNow = 100.0;
        std::array<CutBiquad, maxChannels> bq;

        void snap (bool isHigh, bool on, double hz, double rate) noexcept
        {
            high  = isHigh;
            g     = on ? 1.0f : 0.0f;
            hzNow = hz;
            phase = 0;
            for (auto& b : bq)
            {
                b.design (high, hzNow, rate);
                b.reset();
            }
        }

        int phase = 0;   // samples into the glide's cadence — carried across calls, so the sweep keeps audio time

        void redesign (double rate) noexcept
        {
            for (auto& b : bq)
                b.design (high, hzNow, rate);
        }

        void run (float* const* ch, int nc, int n, bool on, double hzWant, double rate, float fade, double keep) noexcept
        {
            const float want = on ? 1.0f : 0.0f;

            // Off and faded out: not run, its memory already cleared when it got here. It stands at the frequency asked —
            // coefficients and all — so it comes back where its word says; and the cadence keeps time while it waits.
            if (g <= 0.0f && ! on)
            {
                if (hzNow < hzWant || hzNow > hzWant)
                {
                    hzNow = hzWant;
                    redesign (rate);
                }
                phase = (int) (((long long) phase + n) % cadence);
                return;
            }

            for (int start = 0; start < n;)
            {
                // A step of the glide on the cadence, counted in samples from the snap whatever the host's blocks:
                // log-evenly, snapped once within 0.1 %.
                if (phase == 0)
                {
                    if (std::abs (std::log (hzWant / hzNow)) > 0.001)
                    {
                        hzNow = std::exp (std::log (hzWant) + (std::log (hzNow) - std::log (hzWant)) * keep);
                        redesign (rate);
                    }
                    else if (hzNow < hzWant || hzNow > hzWant)   // an exact landing, once
                    {
                        hzNow = hzWant;
                        redesign (rate);
                    }
                }

                const int  len    = std::min (n - start, cadence - phase);
                const bool fading = std::abs (want - g) > 0.0f;
                for (int c = 0; c < nc; ++c)
                {
                    float* d  = ch[c] + start;
                    auto&  b  = bq[(size_t) c];
                    float  gg = g;
                    for (int i = 0; i < len; ++i)
                    {
                        const float y = b.process (d[i]);
                        if (fading)
                        {
                            gg = on ? std::min (1.0f, gg + fade) : std::max (0.0f, gg - fade);
                            d[i] = d[i] + (y - d[i]) * gg;
                        }
                        else
                            d[i] = y;
                    }
                    // A sample that was not a number has passed (the door downstream silences it), but it must not
                    // stay in the recursion: a NaN in z1/z2 feeds itself and the cut would be NaN from then on.
                    if (! std::isfinite (b.z1) || ! std::isfinite (b.z2))
                        b.reset();
                    if (c == nc - 1)
                        g = gg;
                }

                start += len;
                phase  = (phase + len) % cadence;

                // Faded out: the filter's memory goes with it; the rest of the block is the input, as it stands.
                if (g <= 0.0f && ! on)
                {
                    for (auto& b : bq)
                        b.reset();
                    phase = (int) (((long long) phase + (n - start)) % cadence);
                    return;
                }
            }
        }
    };

    std::array<Cut, 2> cuts;
    double fs = 48000.0, glideKeep = 0.9;
    float  fadeStep = 0.001f;
    bool   settled = false;
};

} // namespace felitronics::eq
