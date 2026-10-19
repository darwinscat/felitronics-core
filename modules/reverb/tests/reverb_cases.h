// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The reverb's render scenarios, shared by the golden tests here and by the out-of-tree parity harness that
// produced the goldens from the app's original class. JUCE-free; the golden inputs (impulse, LCG noise) call
// no libm, so every row regenerates the same INPUT bit for bit.
#pragma once

#include <felitronics/reverb/AmpReverb.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace felitronics::test::reverbcases
{

struct Lcg
{
    std::uint32_t s;
    float next() { s = s * 1664525u + 1013904223u; return (float) (std::int32_t) s * 4.6566128730773926e-10f; }
};

enum Input { kImpulse = 0, kNoiseBurst = 1, kSweep = 2 };

struct Knobs
{
    int   character  = 3;        // Plate, the shipping room
    float mix        = 0.3f;
    float decay      = 1.0f;
    float predelayMs = 0.0f;
    float hpfHz      = 120.0f;
};

struct Event { double atSec; Knobs knobs; bool feed; bool reset; };

struct Case
{
    const char* name;
    double sampleRate;
    int    block;
    int    nch;
    Input  input;
    double burstSec;
    double seconds;
    Knobs  knobs;
    std::vector<Event> events;
};

template <class SinFn, class ExpFn>
std::vector<std::vector<float>> makeInput (const Case& c, SinFn sinFn, ExpFn expFn)
{
    const int len = (int) (c.seconds * c.sampleRate);
    std::vector<std::vector<float>> io ((std::size_t) c.nch, std::vector<float> ((std::size_t) len, 0.0f));
    const int burst = std::min (len, (int) (c.burstSec * c.sampleRate));

    for (int ch = 0; ch < c.nch; ++ch)
    {
        auto& x = io[(std::size_t) ch];
        if (c.input == kImpulse)
        {
            if (ch == 0) x[0] = 1.0f;
            else if (len > 5) x[5] = -0.7f;
        }
        else if (c.input == kNoiseBurst)
        {
            Lcg r { 0x89ABCDEu + 104729u * (std::uint32_t) ch };
            for (int i = 0; i < burst; ++i)
                x[(std::size_t) i] = 0.5f * r.next();
        }
        else
        {
            const double f0 = 20.0, T = (double) burst / c.sampleRate, k = 6.907755278982137;   // ln(20000/20)
            for (int i = 0; i < burst; ++i)
            {
                const double t = (double) i / c.sampleRate;
                const double ph = 2.0 * 3.14159265358979323846 * f0 * T / k * (expFn (t / T * k) - 1.0);
                x[(std::size_t) i] = (float) (ch == 0 ? 0.5 * sinFn (ph) : 0.4 * sinFn (ph + 1.5707963267948966));
            }
        }
    }
    return io;
}

/** One render; the knobs are applied at the top of EVERY block, as a plugin's processBlock does. */
template <class Adapter>
void render (const Case& c, Adapter& ad, std::vector<std::vector<float>>& io)
{
    const int len = io.empty() ? 0 : (int) io[0].size();
    ad.prepare (c.sampleRate, c.block, c.nch);

    Knobs k = c.knobs;
    bool feed = true;
    std::size_t ev = 0;

    for (int pos = 0; pos < len; pos += c.block)
    {
        const int n = std::min (c.block, len - pos);
        while (ev < c.events.size() && c.events[ev].atSec * c.sampleRate <= (double) pos)
        {
            k = c.events[ev].knobs;
            feed = c.events[ev].feed;
            if (c.events[ev].reset)
                ad.reset();
            ++ev;
        }

        ad.apply (k);
        float* planes[2] = { io[0].data() + pos, c.nch > 1 ? io[1].data() + pos : nullptr };
        ad.process (planes, c.nch, n, feed);
    }
}

struct CoreAdapter
{
    felitronics::reverb::AmpReverb r;
    bool ok = true;

    void prepare (double sr, int block, int nch) { ok = r.prepare (sr, block, nch) && ok; }
    void reset() { r.reset(); }
    void apply (const Knobs& k)
    {
        felitronics::reverb::AmpReverbParams p;
        p.character  = (felitronics::reverb::Character) k.character;
        p.mix        = k.mix;
        p.decay      = k.decay;
        p.predelayMs = k.predelayMs;
        p.hpfHz      = k.hpfHz;
        r.setParams (p);
    }
    void process (float* const* io, int nch, int n, bool feed) { ok = r.process (io, nch, n, feed) && ok; }
};

inline std::vector<Case> goldenCases()
{
    std::vector<Case> v;
    static const char* names[] = { "ambience", "room", "hall", "plate", "spring", "modulated" };

    for (int c = 0; c < 6; ++c)
    {
        Knobs k; k.character = c; k.mix = 1.0f; k.decay = 1.3f; k.predelayMs = 12.0f; k.hpfHz = 150.0f;
        v.push_back ({ names[c], 48000.0, 64, 2, kNoiseBurst, 0.05, 0.5, k, {} });
    }

    // Mono at 44.1 kHz on an odd block, through a character change, a mix move, a stand-by and a restart.
    Knobs m; m.character = 2; m.mix = 0.8f; m.decay = 2.0f; m.predelayMs = 0.0f; m.hpfHz = 40.0f;
    Knobs m2 = m; m2.character = 5; m2.mix = 0.4f;
    Knobs m3 = m2; m3.predelayMs = 100.0f; m3.hpfHz = 500.0f; m3.decay = 0.5f;
    v.push_back ({ "mono-moves-333", 44100.0, 333, 1, kNoiseBurst, 0.08, 0.7, m,
                   { { 0.15, m2, true, false }, { 0.3, m3, true, false }, { 0.45, m3, false, false },
                     { 0.6, m3, true, true } } });

    // The impulse at 96 kHz, where the tunings are scaled furthest from the 44.1 kHz they were written for.
    Knobs h; h.character = 2; h.mix = 1.0f; h.decay = 0.7f; h.predelayMs = 37.0f; h.hpfHz = 80.0f;
    v.push_back ({ "hall-impulse-96k", 96000.0, 512, 2, kImpulse, 0.0, 0.4, h, {} });
    return v;
}

} // namespace felitronics::test::reverbcases
