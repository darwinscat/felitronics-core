// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The delay's render scenarios, shared by the golden tests here and by the out-of-tree parity harness that
// produced the goldens from the app's original class. JUCE-free, and libm-free on the golden path (impulse
// and LCG noise only), so every row regenerates the same INPUT bit for bit.
#pragma once

#include <felitronics/delay/TapeDelay.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace felitronics::test::delaycases
{

struct Lcg
{
    std::uint32_t s;
    float next() { s = s * 1664525u + 1013904223u; return (float) (std::int32_t) s * 4.6566128730773926e-10f; }
};

enum Input { kImpulse = 0, kNoiseBurst = 1, kSweep = 2 };

/** The knobs as the app's face holds them: DARK in percent, TIME free or synced. */
struct Knobs
{
    float timeMs   = 350.0f;
    bool  sync     = false;
    int   division = 5;
    double bpm     = 120.0;
    float repeats  = 0.35f;
    float darkPct  = 45.0f;
    float offsetMs = 0.0f;
    float mix      = 0.25f;
};

struct Event { double atSec; Knobs knobs; bool feed; bool reset; };

struct Case
{
    const char* name;
    double sampleRate;
    int    block;
    int    nch;
    Input  input;
    double burstSec;     // noise burst / sweep length; the rest is silence (the tail)
    double seconds;      // render length
    Knobs  knobs;
    std::vector<Event> events;
};

/** The input planes, deterministic. The sweep (harness only) calls libm; impulse and noise do not. */
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
            Lcg r { 0x1234567u + 7919u * (std::uint32_t) ch };
            for (int i = 0; i < burst; ++i)
                x[(std::size_t) i] = 0.5f * r.next();
        }
        else
        {
            // log sweep 20 Hz .. 20 kHz over the burst, the right channel in quadrature at 0.8
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

/** One render through any adapter that offers prepare / apply / process / reset. The knobs are applied at the
    top of EVERY block, as a plugin's processBlock does; an event lands on the first block starting at or
    after its time. */
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

/** The core stage, driven the way the app drives its own. */
struct CoreAdapter
{
    felitronics::delay::TapeDelay d;
    bool ok = true;

    void prepare (double sr, int block, int nch) { ok = d.prepare (sr, block, nch) && ok; }
    void reset() { d.reset(); }
    void apply (const Knobs& k)
    {
        felitronics::delay::TapeDelayParams p;
        p.timeMs   = felitronics::delay::timeMsFor (k.sync, k.division, k.bpm, k.timeMs);
        p.repeats  = k.repeats;
        p.darkHz   = felitronics::delay::darkHzForPercent (k.darkPct);
        p.offsetMs = k.offsetMs;
        p.mix      = k.mix;
        d.setParams (p);
    }
    void process (float* const* io, int nch, int n, bool feed) { ok = d.process (io, nch, n, feed) && ok; }
};

/** The golden scenarios: small, but each one through a different corner of the stage. */
inline std::vector<Case> goldenCases()
{
    std::vector<Case> v;

    Knobs a; a.mix = 1.0f;   // the default voice, wet at unity so the echoes are the observable
    v.push_back ({ "default-impulse", 48000.0, 64, 2, kImpulse, 0.0, 0.8, a, {} });

    Knobs b; b.timeMs = 120.0f; b.repeats = 0.95f; b.darkPct = 80.0f; b.offsetMs = 23.0f; b.mix = 1.0f;
    v.push_back ({ "dark-bloom-offset-noise", 48000.0, 64, 2, kNoiseBurst, 0.05, 0.8, b, {} });

    Knobs c; c.sync = true; c.division = 9; c.bpm = 140.0; c.repeats = 0.6f; c.darkPct = 10.0f; c.offsetMs = -12.0f; c.mix = 0.5f;
    Knobs c2 = c; c2.division = 5;                               // 1/8T -> 1/4: the head glides, the repeats bend
    Knobs c3 = c2; c3.offsetMs = 25.0f; c3.repeats = 1.0f;       // offset swings across, the loop at full
    v.push_back ({ "synced-glide-mono-333", 44100.0, 333, 1, kNoiseBurst, 0.08, 0.9, c,
                   { { 0.25, c2, true, false }, { 0.45, c3, true, false }, { 0.6, c3, false, false } } });

    Knobs d; d.timeMs = 40.0f; d.repeats = 0.8f; d.darkPct = 0.0f; d.offsetMs = -30.0f; d.mix = 0.7f;
    Knobs d2 = d; d2.timeMs = 300.0f; d2.offsetMs = 0.0f;
    v.push_back ({ "stereo-reset-96k-512", 96000.0, 512, 2, kNoiseBurst, 0.03, 0.5, d,
                   { { 0.2, d2, true, false }, { 0.35, d2, true, true } } });
    return v;
}

} // namespace felitronics::test::delaycases
