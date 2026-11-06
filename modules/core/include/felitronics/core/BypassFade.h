// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>   // approximatelyEqual — JUCE's, so the fade decides as it did in the plug-in

#include <algorithm>
#include <cmath>

namespace felitronics::core
{

/** THE ONE FADE, for every link that REPLACES the signal rather than adding to it. Moved here from OrbitAmp
    (src/core/BypassFade.h), JUCE-free, its arithmetic unchanged: the amps' chain and the web's play the same fade.

    An additive link — the delay, the room — needs none of this: stop feeding it and its own tail is the fade. A link
    that replaces has no tail to ride out; drop it out of the path and the waveform takes a step, and a step is a click.
    So it is crossfaded against the signal it was handed: dry at one end, what the block made at the other.

    It costs NOTHING while nothing is happening. A fade that is finished reports `moving == false`, and the call site
    then either runs the block plainly or skips it plainly — no copy of the signal, no per-sample multiply. Only the
    fifteen milliseconds after a switch pay for anything.

    The owner's thread sets nothing here: the audio thread reads its switches itself and this only follows what it
    read. Real-time safe; no allocation. */
class BypassFade
{
public:
    /** Milliseconds. Long enough that no step survives it, short enough that a switch still feels like a switch. */
    static constexpr float lengthMs = 15.0f;

    void prepare (double sampleRate) noexcept
    {
        perSample = sampleRate > 0.0
                      ? (float) (1000.0 / (sampleRate * (double) lengthMs))
                      : 1.0f;
    }

    /** No fade on the first block after a prepare: a chain that arrives switched off should be silent from its first
        sample, not fade in from wherever the last session left this. */
    void snapTo (bool on) noexcept { g = on ? 1.0f : 0.0f; }

    /** Where the blend stands at the start of this block and where it will stand at its end. `moving` false means the
        whole block sits at one end: nothing to blend, nothing to copy. */
    struct Span
    {
        float from = 1.0f, to = 1.0f;
        bool  moving = false;

        /** How many samples OF THIS BLOCK the fade actually occupies — a fifteen-millisecond fade stays fifteen
            milliseconds in a block of four thousand samples. */
        int   ramp = 0;
    };

    Span advance (int numSamples, bool target) noexcept
    {
        const float want = target ? 1.0f : 0.0f;

        if (approximatelyEqual (g, want))
            return { g, g, false };

        // A block with no samples moves nothing — but a fade under way is still under way: reporting it stopped would
        // tell the chain an outgoing link is gone, and the link would be cleared with its tail still sounding.
        if (numSamples <= 0)
            return { g, g, true, 0 };

        const float from   = g;
        const float needed = std::abs (want - g) / std::max (1.0e-9f, perSample);
        const int   ramp   = std::clamp ((int) std::ceil (needed), 1, numSamples);
        const float step   = perSample * (float) ramp;

        g = target ? std::min (want, g + step) : std::max (want, g - step);
        return { from, g, true, ramp };
    }

    /** Is the block's output wanted at all this block — because it is in, or because it is still on its way out and
        the tail of the crossfade needs something to fade FROM. */
    static bool runs (const Span& s, bool target) noexcept { return target || s.moving; }

    /** dry + (wet - dry) * g, with g walking from `from` to `to` across the block. `wet` is the buffer the block just
        wrote in place; `dry` is what it was handed. */
    static void blend (float* const* wet, const float* const* dry, int numChannels, int numSamples,
                       const Span& s) noexcept
    {
        if (numSamples <= 0)
            return;

        const int   ramp = std::clamp (s.ramp, 1, numSamples);
        const float step = (s.to - s.from) / (float) ramp;

        for (int ch = 0; ch < numChannels; ++ch)
        {
            float*       w = wet[ch];
            const float* d = dry[ch];

            // Walk for as long as the fade actually lasts, then HOLD: a fade that finishes a fifth of the way into a
            // big block must not go on sliding for the other four.
            for (int i = 0; i < numSamples; ++i)
            {
                const float gi = i < ramp ? s.from + step * (float) i : s.to;
                w[i] = d[i] + (w[i] - d[i]) * gi;
            }
        }
    }

private:
    float g         = 1.0f;
    float perSample = 1.0f;
};

} // namespace felitronics::core
