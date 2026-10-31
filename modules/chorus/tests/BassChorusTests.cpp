// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// chorus::BassChorus and its one-slider law. The effect: amount 0 is the dry bit for bit and computes nothing,
// from the very sample the mix lands; the voice's band is LR4's high half, so it and the band it leaves sum
// flat; the lows stay centred — the Side a tone under the crossover makes is bounded by the voice's own LR4
// attenuation there; the measured pitch deviation is the set one, the full depth on the lowest centre too; the
// two voices in quadrature; no knob clicks, the slider included. The stage
// contract: the law-11 properties at the module's own width ceiling, prepare binds and a refused one disarms,
// law 8a re-slicing, reset is a restart, a stopped channel drops its voice (11a), a pause is silence (11c),
// poison, law 8, no allocation. At 44.1, 48 and 96 kHz.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter
#include <law11_call_contract.h>
#include <law11c_pause.h>

#include <felitronics/chorus/BassChorus.h>
#include <felitronics/chorus/ChorusLaw.h>
#include <felitronics/eq/Crossover2.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using chorus::BassChorus;
using chorus::BassChorusParams;

namespace
{
constexpr double kPi = 3.14159265358979323846;
const double kRates[] = { 44100.0, 48000.0, 96000.0 };

using Buf = std::vector<float>;

std::string num (double v, int prec = 4)
{
    char b[64];
    std::snprintf (b, sizeof b, "%.*g", prec, v);
    return b;
}

std::string khz (double fs) { return num (fs / 1000.0, 3) + " kHz"; }

Buf tone (int n, double hz, double fs, float amp, double ph0 = 0.0)
{
    Buf v ((std::size_t) n);
    for (int i = 0; i < n; ++i)
        v[(std::size_t) i] = amp * (float) std::sin (2.0 * kPi * hz * (double) i / fs + ph0);
    return v;
}

Buf noise (int n, std::uint32_t seed, float amp)
{
    Buf v ((std::size_t) n);
    for (auto& s : v)
    {
        seed = seed * 1664525u + 1013904223u;
        s = amp * ((float) (seed >> 8) / 8388608.0f - 1.0f);
    }
    return v;
}

Buf plus (Buf a, const Buf& b)
{
    for (std::size_t i = 0; i < a.size(); ++i) a[i] += b[i];
    return a;
}

/** L (and R when given) through the stage in blocks of `block`. */
void render (BassChorus& s, Buf& l, Buf* r, int block)
{
    const int n = (int) l.size();
    for (int pos = 0; pos < n; pos += block)
    {
        const int m = std::min (block, n - pos);
        float* io[2] = { l.data() + pos, r != nullptr ? r->data() + pos : nullptr };
        test::run (s.process (io, r != nullptr ? 2 : 1, m));
    }
}

bool bitsEqual (const Buf& a, const Buf& b)
{
    return a.size() == b.size() && std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0;
}

/** Complex amplitude at `hz` over [a, b), Hann-windowed. */
std::complex<double> phasor (const Buf& x, int a, int b, double hz, double fs)
{
    std::complex<double> acc {};
    double wsum = 0.0;
    const int n = b - a;
    for (int i = 0; i < n; ++i)
    {
        const double w  = 0.5 - 0.5 * std::cos (2.0 * kPi * ((double) i + 0.5) / (double) n);
        const double ph = 2.0 * kPi * hz * (double) (a + i) / fs;
        acc  += w * (double) x[(std::size_t) (a + i)] * std::complex<double> (std::cos (ph), -std::sin (ph));
        wsum += w;
    }
    return 2.0 * acc / wsum;
}

/** |DTFT| of an impulse response at `hz`, dB. */
double magDb (const std::vector<double>& h, double hz, double fs)
{
    const std::complex<double> rot = std::polar (1.0, -2.0 * kPi * hz / fs);
    std::complex<double> e (1.0, 0.0), acc {};
    for (std::size_t i = 0; i < h.size(); ++i)
    {
        acc += h[i] * e;
        e *= rot;
        if ((i & 1023u) == 1023u)
            e /= std::abs (e);
    }
    return 20.0 * std::log10 (std::abs (acc));
}

//==============================================================================
// The law-11 properties, with the module's OWN width ceiling (2): the census macro pins it at the house maximum.
struct A_BassChorus
{
    using T = BassChorus;
    static const char* name() { return "felitronics::chorus::BassChorus"; }
    static constexpr int    maxCh           = 2;
    static constexpr bool   narrowOk        = true;
    static constexpr bool   chunks          = true;
    static constexpr bool   inPlace         = true;
    static constexpr bool   needsPrepare    = true;
    static constexpr bool   gapVisible      = true;
    static constexpr double silenceTol      = 1.0e-6;
    static constexpr bool   hasPrepareWidth = true;
    static constexpr int    prepareCeiling  = BassChorus::kMaxChannels;
    static bool prepWidth (T& s, int w) { return s.prepare (test::law11::kFs, test::law11::kMaxBlock, w); }
    static double blockObservable (T&) { return 0.0; }
    static bool prep (T& s)
    {
        // The lowest centre (5 ms) at full mix, so the short fixtures hear the voice.
        s.setParams ({ 1.0f, 1.0f, 3.0f, 150.0f, 2.0f });
        return s.prepare (test::law11::kFs, test::law11::kMaxBlock, 2);
    }
    static bool run (T& s, float* const* io, int nch, int n) { return s.process (io, nch, n); }
    static double observable (T&) { return 0.0; }
};

void contract()
{
    group ("contract: prepare binds and a refused one disarms; n == 0 is a no-op; the law-11 properties");

    BassChorus s;
    ok (! s.prepare (7999.0, 64, 2) && ! s.prepare (768001.0, 64, 2) && ! s.prepare (48000.0, 0, 2)
        && ! s.prepare (48000.0, 64, 0) && ! s.prepare (48000.0, 64, 3)
        && ! s.prepare (std::numeric_limits<double>::quiet_NaN(), 64, 2), "prepare refuses a rate, block or width it cannot honour");
    ok (s.prepare (8000.0, 1, 1) && s.prepare (768000.0, 64, 2) && s.prepare (48000.0, 64, 2), "...and takes the edges of its ranges");

    Buf l = noise (256, 7u, 0.5f), r = noise (256, 9u, 0.5f);
    const Buf l0 = l, r0 = r;
    float* io[2] = { l.data(), r.data() };
    ok (! s.prepare (48000.0, 64, 5), "precondition: a bad re-prepare is refused");
    ok (! s.process (io, 2, 256) && ! s.process (io, 1, 256) && ! s.process (io, 0, 256), "a refused re-prepare DISARMS: process refuses");
    ok (bitsEqual (l, l0) && bitsEqual (r, r0), "...and touches nothing");

    test::run (s.prepare (48000.0, 64, 2));
    ok (s.process (io, 2, 0) && bitsEqual (l, l0), "n == 0 is accepted and moves nothing");
    ok (! s.process (nullptr, 2, 16), "a null buffer with a width is refused");

    test::law11::allProperties<A_BassChorus>();
}

//==============================================================================
void theLaw()
{
    group ("the one-slider law: amount drives depth and mix together; the menu's knobs are left alone");
    const BassChorusParams menu { 9.0f, 0.9f, 2.5f, 333.0f, 12.0f };
    const auto p0 = chorus::paramsForAmount (0.0f, menu), ph = chorus::paramsForAmount (0.5f, menu),
               p1 = chorus::paramsForAmount (1.0f, menu);
    ok (! (p0.mix > 0.0f) && ! (p0.depthMs > 0.0f), "amount 0: mix 0 and depth 0");
    approx (ph.mix, 0.5, 1e-7, "amount 0.5: mix 0.5");
    approx (ph.depthMs, 1.5, 1e-6, "amount 0.5: depth 1.5 ms");
    approx (p1.mix, 1.0, 0.0, "amount 1: mix 1 (the voice at the dry's level)");
    approx (p1.depthMs, 3.0, 0.0, "amount 1: depth 3 ms");
    ok (core::sameBits (ph.rateHz, 2.5f) && core::sameBits (ph.crossoverHz, 333.0f) && core::sameBits (ph.delayMs, 12.0f),
        "the menu's rate, crossover and centre delay pass through untouched");
    const auto pn = chorus::paramsForAmount (std::numeric_limits<float>::quiet_NaN(), menu), pb = chorus::paramsForAmount (7.0f, menu),
               pm = chorus::paramsForAmount (-std::numeric_limits<float>::infinity(), menu);
    ok (core::sameBits (pn.mix, menu.mix) && core::sameBits (pn.depthMs, menu.depthMs) && core::sameBits (pm.mix, menu.mix),
        "a non-finite amount is ignored: the menu (the stage's current params) comes back untouched");
    ok (core::sameBits (pb.mix, 1.0f) && core::sameBits (pb.depthMs, 3.0f), "past 1 is 1");

    // The read-out: a pure function of the knobs, the LARGER (flat) excursion.
    BassChorusParams d = chorus::paramsForAmount (1.0f);
    d.rateHz = 0.7f;
    const double x = 2.0 * kPi * 0.7 * 0.003;
    approx (chorus::pitchDeviationCents (d), -1200.0 * std::log2 (1.0 - x), 1e-3,
            "amount 1 at the default 0.7 Hz: the flat excursion, 23.0 cents (the sharp one is 22.7)");
    const BassChorusParams wild { BassChorus::kMaxDepthMs, 1.0f, BassChorus::kMaxRateHz, 200.0f, BassChorus::kMaxDelayMs };
    approx (chorus::pitchDeviationCents (wild), 1713.0, 1.0, "10 Hz x 10 ms: 1713 cents flat, the larger of the two (+844 sharp)");
    const BassChorusParams capped { BassChorus::kMaxDepthMs, 1.0f, 1.0f, 200.0f, BassChorus::kMinDelayMs };
    const double held = (double) (BassChorus::kMinDelayMs - BassChorus::kHeadroomMs) * 0.001;
    approx (chorus::pitchDeviationCents (capped), -1200.0 * std::log2 (1.0 - 2.0 * kPi * 1.0 * held), 1e-3,
            "depth 10 ms on the lowest centre: the swing as the stage holds it, the centre less the headroom");
}

//==============================================================================
void amountZero()
{
    group ("amount 0: the dry bit for bit, nothing computed — fresh, and after the slider comes down");
    for (double fs : kRates)
    {
        BassChorus s;
        s.setParams (chorus::paramsForAmount (0.0f, s.params()));
        test::run (s.prepare (fs, 512, 2));
        Buf l = plus (noise ((int) fs, 11u, 0.7f), tone ((int) fs, 61.0, fs, 0.3f)), r = noise ((int) fs, 13u, 0.7f);
        l[100] = std::numeric_limits<float>::quiet_NaN();
        r[200] = -0.0f;
        const Buf l0 = l, r0 = r;
        int pos = 0;
        for (int b : { 1, 7, 64, 333, 512, 4096 })
            for (int k = 0; k < 3 && pos + b <= (int) l.size(); ++k, pos += b)
            {
                float* io[2] = { l.data() + pos, r.data() + pos };
                test::run (s.process (io, 2, b));
            }
        float* io[2] = { l.data() + pos, r.data() + pos };
        test::run (s.process (io, 2, (int) l.size() - pos));
        ok (bitsEqual (l, l0) && bitsEqual (r, r0) && s.isIdle(),
            "fresh at amount 0, " + khz (fs) + ": every bit of the input comes back (a NaN and a -0 included), and the stage idles");

        // From a full chorus down to 0: the mix lands on an exact 0, the stage idles, the dry is bit-exact.
        BassChorus c;
        c.setParams (chorus::paramsForAmount (0.8f, c.params()));
        test::run (c.prepare (fs, 512, 2));
        Buf a = tone ((int) (0.5 * fs), 997.0, fs, 0.5f), b = a;
        render (c, a, &b, 512);
        ok (! c.isIdle(), "precondition: at 0.8 the stage runs");
        c.setParams (chorus::paramsForAmount (0.0f, c.params()));
        int landed = -1;
        for (int k = 0; k < (int) (2.0 * fs) / 512 && landed < 0; ++k)
        {
            Buf x = tone (512, 997.0, fs, 0.5f), y = x;
            float* q[2] = { x.data(), y.data() };
            test::run (c.process (q, 2, 512));
            if (c.isIdle())
                landed = (k + 1) * 512;
        }
        ok (landed > 0 && landed <= (int) (0.025 * fs) + 512, "...down to 0 it idles once the counted mix ramp lands, after "
                                                       + num (1000.0 * landed / fs) + " ms (the 20 ms ramp and the block it ends in)");
        Buf x = plus (tone (4096, 997.0, fs, 0.5f), noise (4096, 5u, 0.2f)), y = noise (4096, 6u, 0.2f);
        const Buf x0 = x, y0 = y;
        render (c, x, &y, 300);
        ok (bitsEqual (x, x0) && bitsEqual (y, y0), "...and from then on the dry comes back bit for bit, " + khz (fs));

        // A wake is a restart: two stages with different pasts, idled at different grid phases, wake into the
        // same programme with the same bits.
        BassChorus wa, wb;
        for (auto* st : { &wa, &wb })
        {
            st->setParams (chorus::paramsForAmount (0.7f, st->params()));
            test::run (st->prepare (fs, 512, 2));
        }
        Buf h1 = noise (7777, 41u, 0.8f), h2 = noise (7777, 42u, 0.8f), g1 = tone (12345, 300.0, fs, 0.9f), g2 = g1;
        render (wa, h1, &h2, 512);
        render (wb, g1, &g2, 100);
        for (auto* st : { &wa, &wb })
            st->setParams (chorus::paramsForAmount (0.0f, st->params()));
        Buf q1 (3333, 0.1f), q2 (3333, 0.1f), q3 (5000, 0.1f), q4 (5000, 0.1f);
        render (wa, q1, &q2, 512);
        render (wb, q3, &q4, 77);
        ok (wa.isIdle() && wb.isIdle(), "precondition: both idle, " + khz (fs));
        Buf w1 = plus (tone (20000, 997.0, fs, 0.5f), noise (20000, 43u, 0.1f)), w2 = noise (20000, 44u, 0.3f), w3 = w1, w4 = w2;
        for (auto* st : { &wa, &wb })
            st->setParams (chorus::paramsForAmount (0.6f, st->params()));
        render (wa, w1, &w2, 512);
        render (wb, w3, &w4, 512);
        ok (bitsEqual (w1, w3) && bitsEqual (w2, w4), "a wake from idle is a restart: different pasts, the same bits after, " + khz (fs));
    }
}

//==============================================================================
void splitAndSumFlat()
{
    group ("depth 0: the voice and the band it leaves sum flat, ±0.1 dB 20 Hz – 20 kHz; the voice IS LR4's high half");
    struct Case { double fs; float delayMs; float xover; };
    for (const Case c : { Case { 44100.0, 10.0f, 200.0f }, Case { 48000.0, 8.0f, 200.0f }, Case { 96000.0, 8.0f, 200.0f },
                          Case { 48000.0, 8.0f, 100.0f },   Case { 48000.0, 8.0f, 500.0f }, Case { 48000.0, 8.0f, 1000.0f } })
    {
        BassChorus s;
        s.setParams ({ 0.0f, 1.0f, 0.7f, c.xover, c.delayMs });
        test::run (s.prepare (c.fs, 512, 1));
        Buf warm ((std::size_t) c.fs, 0.0f);   // a second of silence: the door is open after it
        render (s, warm, nullptr, 512);

        const int N = 1 << 15;
        Buf x ((std::size_t) N, 0.0f);
        x[0] = 1.0f;
        Buf y = x;
        render (s, y, nullptr, 512);

        // The voice, and LR4's LOW half of the same impulse from eq's own crossover, behind the same head.
        const int head = (int) std::lround ((double) c.delayMs * c.fs / 1000.0);
        eq::DeterministicCrossover2 xo;
        xo.prepare (c.fs, 1);
        xo.setFrequency (c.xover);
        std::vector<double> voice ((std::size_t) N), sum ((std::size_t) N, 0.0);
        for (int i = 0; i < N; ++i)
            voice[(std::size_t) i] = (double) y[(std::size_t) i] - (double) x[(std::size_t) i];
        for (int i = 0; i + head < N; ++i)
        {
            float lo = 0.0f, hi = 0.0f;
            xo.processSample (0, x[(std::size_t) i], lo, hi);
            sum[(std::size_t) (i + head)] = (double) lo;
        }
        for (int i = 0; i < N; ++i)
            sum[(std::size_t) i] += voice[(std::size_t) i];

        double worst = 0.0;
        for (int k = 0; k <= 240; ++k)
        {
            const double f = 20.0 * std::pow (1000.0, (double) k / 240.0);
            worst = std::max (worst, std::fabs (magDb (sum, f, c.fs)));
        }
        const std::string at = khz (c.fs) + ", crossover " + num (c.xover) + " Hz";
        ok (worst <= 0.1, "voice + LR4 low is flat: worst " + num (worst, 3) + " dB, " + at);
        // ...and the voice is the high HALF, not just some complement: -6.02 dB at the corner, -24.6 an octave
        // below (x - LP4 would also sum flat with LP4, and sits at +3.5 dB at the corner).
        approx (magDb (voice, c.xover, c.fs), -6.0206, 0.1, "the voice at the crossover, " + at);
        approx (magDb (voice, 0.5 * c.xover, c.fs), 20.0 * std::log10 (1.0 / 17.0), 0.2, "the voice an octave below, " + at);
    }

    group ("the read: a 4-point cubic Hermite, at its worst fraction (0.5)");
    {
        const double fs = 48000.0;
        auto voiceIr = [fs] (float delayMs)
        {
            BassChorus s;
            s.setParams ({ 0.0f, 1.0f, 0.7f, BassChorus::kMinCrossoverHz, delayMs });
            test::run (s.prepare (fs, 512, 1));
            Buf warm ((std::size_t) fs, 0.0f);
            render (s, warm, nullptr, 512);
            Buf x (8192, 0.0f);
            x[0] = 1.0f;
            Buf y = x;
            render (s, y, nullptr, 512);
            std::vector<double> v (y.size());
            for (std::size_t i = 0; i < y.size(); ++i)
                v[i] = (double) y[i] - (double) x[i];
            return v;
        };
        const auto whole = voiceIr (8.0f), half = voiceIr ((float) (384.5 / 48.0));   // 384 and 384.5 samples
        const double loss5 = magDb (half, 5000.0, fs) - magDb (whole, 5000.0, fs);
        const double loss10 = magDb (half, 10000.0, fs) - magDb (whole, 10000.0, fs);
        ok (loss5 > -0.1 && loss5 < 0.01, "at 5 kHz the half-sample read loses " + num (loss5, 3) + " dB (> -0.1; a linear read loses 0.47)");
        ok (loss10 > -0.7 && loss10 < 0.01, "at 10 kHz it loses " + num (loss10, 3) + " dB (> -0.7; a linear read loses 2.0)");
    }
}

//==============================================================================
void lowsMonoInPhase()
{
    group ("the lows stay with the dry: two and three octaves under the crossover L ~ R ~ dry (stereo bus, amount 1)");
    for (double fs : kRates)
    {
        BassChorus s;
        BassChorusParams menu;
        menu.crossoverHz = 200.0f;
        s.setParams (chorus::paramsForAmount (1.0f, menu));
        test::run (s.prepare (fs, 512, 2));
        const int N = (int) (3.0 * fs);
        const Buf dry = plus (plus (tone (N, 25.0, fs, 0.4f), tone (N, 50.0, fs, 0.4f, 1.0)), tone (N, 1000.0, fs, 0.3f, 2.0));
        Buf l = dry, r = dry;
        render (s, l, &r, 512);
        const int a = (int) fs, b = N;

        struct Probe { double hz; double tol; const char* what; };
        for (const Probe p : { Probe { 50.0, 0.01, "two octaves below (fc/4): -40 dB" }, Probe { 25.0, 0.001, "three octaves below (fc/8): -60 dB" } })
        {
            const auto X = phasor (dry, a, b, p.hz, fs), L = phasor (l, a, b, p.hz, fs), R = phasor (r, a, b, p.hz, fs);
            const double side = std::abs (L - R) / std::abs (X);
            const double off  = std::max (std::abs (L / X - 1.0), std::abs (R / X - 1.0));
            ok (side <= p.tol, std::string ("|L - R| re the dry at ") + num (p.hz) + " Hz, " + p.what + ": "
                               + num (20.0 * std::log10 (side), 3) + " dB, " + khz (fs));
            ok (off <= p.tol, std::string ("each side IS the dry at ") + num (p.hz) + " Hz: |L/dry - 1| = " + num (off, 3)
                              + " (" + num (20.0 * std::log10 (1.0 + off), 2) + " dB, " + num (std::asin (std::min (off, 1.0)) * 180.0 / kPi, 2)
                              + " deg at most), " + khz (fs));
        }

        // The precondition that makes the above mean something: above the crossover the voices DO spread.
        double sideHi = 0.0, mid = 0.0;
        for (int i = a; i < b; ++i)
        {
            const double d = (double) l[(std::size_t) i] - (double) r[(std::size_t) i];
            sideHi += d * d;
            mid += (double) dry[(std::size_t) i] * (double) dry[(std::size_t) i];
        }
        ok (sideHi / mid > 0.01, "precondition: above the crossover L and R differ (side/dry energy " + num (sideHi / mid, 3) + "), " + khz (fs));
    }
}

//==============================================================================
/** The voice's instantaneous frequency ratio minus 1, one value per cycle, from upward zero crossings. */
std::vector<double> deviation (const Buf& v, int start, double fs, double f0)
{
    std::vector<double> t, out;
    for (std::size_t i = (std::size_t) start + 1; i < v.size(); ++i)
        if (v[i - 1] < 0.0f && v[i] >= 0.0f)
            t.push_back ((double) (i - 1) + (double) v[i - 1] / ((double) v[i - 1] - (double) v[i]));
    for (std::size_t k = 1; k < t.size(); ++k)
        out.push_back (fs / (t[k] - t[k - 1]) / f0 - 1.0);
    return out;
}

void pitchDeviation()
{
    group ("the voice's pitch deviation is the set one, measured on a 1 kHz sine; the two voices in quadrature");
    struct Case { double fs; float rate; float depth; };
    for (const Case c : { Case { 44100.0, 2.0f, 3.0f }, Case { 48000.0, 0.7f, 1.5f }, Case { 96000.0, 2.0f, 3.0f }, Case { 48000.0, 5.0f, 0.5f } })
    {
        const double want = 2.0 * kPi * (double) c.rate * (double) c.depth * 0.001;   // peak ratio - 1
        BassChorus s, m;
        const BassChorusParams p { c.depth, 1.0f, c.rate, 200.0f, 8.0f };
        s.setParams (p);
        m.setParams (p);
        test::run (s.prepare (c.fs, 512, 2));
        test::run (m.prepare (c.fs, 512, 1));
        const int N = (int) (c.fs * (0.3 + 2.0 / (double) c.rate));
        const Buf x = tone (N, 1000.0, c.fs, 0.5f);
        Buf l = x, r = x, mono = x;
        render (s, l, &r, 512);
        render (m, mono, nullptr, 512);
        ok (bitsEqual (mono, l), "a mono stage is the stereo stage's left voice, bit for bit, " + khz (c.fs));

        Buf vl ((std::size_t) N), vr ((std::size_t) N);
        for (std::size_t i = 0; i < (std::size_t) N; ++i)
        {
            vl[i] = l[i] - x[i];
            vr[i] = r[i] - x[i];
        }
        const int start = (int) (0.3 * c.fs);
        const auto dl = deviation (vl, start, c.fs, 1000.0), dr = deviation (vr, start, c.fs, 1000.0);
        const auto [lmin, lmax] = std::minmax_element (dl.begin(), dl.end());
        const auto [rmin, rmax] = std::minmax_element (dr.begin(), dr.end());
        const std::string at = num (c.rate) + " Hz x " + num (c.depth) + " ms, " + khz (c.fs);
        const std::string got = " (set " + num (100.0 * want, 4) + " %, measured +" + num (100.0 * *lmax, 4) + " / -"
                                + num (-100.0 * *lmin, 4) + " left, +" + num (100.0 * *rmax, 4) + " / -" + num (-100.0 * *rmin, 4) + " right)";
        approx (*lmax, want, 0.02 * want, "left voice: peak sharp = 2pi·rate·depth, " + at + got);
        approx (-*lmin, want, 0.02 * want, "left voice: peak flat, " + at);
        approx (*rmax, want, 0.02 * want, "right voice: peak sharp, " + at);
        approx (-*rmin, want, 0.02 * want, "right voice: peak flat, " + at);

        double sxy = 0.0, sxx = 0.0, syy = 0.0;
        for (std::size_t k = 0; k < std::min (dl.size(), dr.size()); ++k)
        {
            sxy += dl[k] * dr[k];
            sxx += dl[k] * dl[k];
            syy += dr[k] * dr[k];
        }
        const double corr = sxy / std::sqrt (sxx * syy);
        // A quarter of the LFO's period later the right voice's deviation is the left one's: the voices are in
        // QUADRATURE, so their mono sum still swings instead of folding into a standing comb.
        const auto quarter = (std::size_t) std::lround (1000.0 / (4.0 * (double) c.rate));   // in cycles of the tone
        double best = 0.0;
        for (int sign : { -1, 1 })
        {
            double qxy = 0.0, qxx = 0.0, qyy = 0.0;
            for (std::size_t k = quarter; k + quarter < std::min (dl.size(), dr.size()); ++k)
            {
                const double a = dl[k], b = dr[sign > 0 ? k + quarter : k - quarter];
                qxy += a * b;
                qxx += a * a;
                qyy += b * b;
            }
            best = std::max (best, qxy / std::sqrt (qxx * qyy));
        }
        ok (std::fabs (corr) < 0.1 && best > 0.98, "the two voices swing in quadrature (correlation " + num (corr, 3)
                                                     + ", " + num (best, 4) + " a quarter period apart), " + at);
        const double cents = std::max (1200.0 * std::log2 (1.0 + *lmax), -1200.0 * std::log2 (1.0 + *lmin));
        approx (chorus::pitchDeviationCents (p), cents, 0.02 * cents,
                "the published deviation is the measured flat excursion (" + num (cents) + " cents), " + at);
    }
}

//==============================================================================
/** The largest third difference of a channel: a smooth sum of sines keeps it under A·(2 sin(ω/2))³, a step of
    size s puts 2s into it. Far more selective than the first difference, which the tone itself fills. */
double maxThirdDiff (const Buf& y, int from, std::size_t* at = nullptr)
{
    double worst = 0.0;
    for (std::size_t i = (std::size_t) std::max (from, 3); i < y.size(); ++i)
    {
        const double d = std::fabs ((double) y[i] - 3.0 * (double) y[i - 1] + 3.0 * (double) y[i - 2] - (double) y[i - 3]);
        if (d > worst)
        {
            worst = d;
            if (at != nullptr)
                *at = i;
        }
    }
    return worst;
}

void noClicks()
{
    group ("no click: the slider (from idle and back to it), rate, crossover, centre and depth all glide");
    for (double fs : kRates)
    {
        const double f = 997.0;
        const float  A = 0.5f;
        BassChorus s;
        s.setParams (chorus::paramsForAmount (0.0f, s.params()));
        test::run (s.prepare (fs, 64, 2));

        struct Ev { double at; int what; float v; };   // what: 0 amount, 1 rate, 2 crossover, 3 centre, 4 depth
        const Ev evs[] = { { 0.20, 0, 0.8f }, { 0.50, 0, 0.2f }, { 0.80, 1, 6.0f }, { 1.10, 2, 800.0f }, { 1.40, 3, 14.0f },
                           { 1.70, 4, 4.0f }, { 2.00, 0, 1.0f }, { 2.01, 1, 0.3f }, { 2.02, 2, 60.0f }, { 2.30, 0, 0.0f },
                           { 3.20, 0, 0.6f }, { 3.21, 3, 5.0f } };
        const int N = (int) (3.5 * fs);
        Buf l = tone (N, f, fs, A), r = l;
        std::size_t e = 0;
        bool idled = false, wokeAgain = false;
        for (int pos = 0; pos < N; pos += 64)
        {
            while (e < std::size (evs) && (double) pos >= evs[e].at * fs)
            {
                const Ev& ev = evs[e++];
                switch (ev.what)
                {
                    case 0: s.setParams (chorus::paramsForAmount (ev.v, s.params())); break;
                    case 1: s.setRateHz (ev.v); break;
                    case 2: s.setCrossoverHz (ev.v); break;
                    case 3: s.setDelayMs (ev.v); break;
                    default: s.setDepthMs (ev.v); break;
                }
            }
            float* io[2] = { l.data() + pos, r.data() + pos };
            test::run (s.process (io, 2, std::min (64, N - pos)));
            if ((double) pos > 2.3 * fs && s.isIdle())
                idled = true;
            if (idled && ! s.isIdle())
                wokeAgain = true;
        }
        // The largest pitch excursion the schedule makes: 6 Hz x 4 ms on the LFO, plus a 6 ms centre move and a
        // 3.4 ms depth move through their glides. A generous 0.3.
        const double w  = 2.0 * kPi * f / fs;
        const double d3 = (double) A * (std::pow (2.0 * std::sin (0.5 * w), 3.0) + std::pow (2.0 * std::sin (0.5 * w * 1.3), 3.0));
        const double bound = 1.25 * d3;
        std::size_t atL = 0, atR = 0;
        const double gl = maxThirdDiff (l, 0, &atL), gr = maxThirdDiff (r, 0, &atR);
        const double got = std::max (gl, gr);
        ok (got <= bound, "largest third difference " + num (got, 3) + " (at " + num ((double) (gl >= gr ? atL : atR) / fs, 5)
                          + " s) <= the smooth bound " + num (bound, 3) + " (a step of " + num (bound / 2.0, 2)
                          + " would break it), " + khz (fs));
        ok (idled && wokeAgain, "precondition: the stream went idle after 2.3 s and woke at 3.2 s, " + khz (fs));
    }
}

//==============================================================================
void reslicing()
{
    group ("law 8a: any slicing of one stream is bit-identical (glides, the door and the crossover in flight)");
    for (double fs : kRates)
    {
        const int N = (int) fs;
        const Buf srcL = plus (tone (N, 997.0, fs, 0.4f), noise (N, 21u, 0.2f)), srcR = plus (tone (N, 61.0, fs, 0.4f), noise (N, 22u, 0.2f));
        std::vector<Buf> outs;
        for (int block : { 1, 3, 64, 100, 517, N })
        {
            BassChorus s;
            s.setParams ({ 2.0f, 0.7f, 1.3f, 120.0f, 9.0f });
            test::run (s.prepare (fs, 512, 2));
            s.setCrossoverHz (700.0f);   // in flight from the first sample
            s.setDepthMs (4.0f);
            s.setDelayMs (15.0f);
            Buf l = srcL, r = srcR;
            render (s, l, &r, block);
            l.insert (l.end(), r.begin(), r.end());
            outs.push_back (l);
        }
        bool same = true;
        for (const auto& o : outs)
            same = same && bitsEqual (o, outs[0]);
        ok (same, "blocks of 1, 3, 64, 100, 517 and one whole call give the same bits, " + khz (fs));
    }
}

//==============================================================================
void resetAndStops()
{
    group ("reset() is a restart; a stopped channel drops its voice (11a); a pause is silence (11c)");
    const double fs = 48000.0;
    {
        BassChorus a, b;
        test::run (a.prepare (fs, 512, 2));
        test::run (b.prepare (fs, 512, 2));
        Buf x1 = noise (20000, 1u, 0.8f), y1 = noise (20000, 2u, 0.8f), x2 = tone (12345, 300.0, fs, 0.9f), y2 = x2;
        render (a, x1, &y1, 512);
        render (b, x2, &y2, 100);
        a.reset();
        b.reset();
        b.reset();   // idempotent
        Buf pa = plus (tone (30000, 997.0, fs, 0.5f), noise (30000, 3u, 0.1f)), qa = noise (30000, 4u, 0.3f), pb = pa, qb = qa;
        render (a, pa, &qa, 512);
        render (b, pb, &qb, 512);
        ok (bitsEqual (pa, pb) && bitsEqual (qa, qb), "two stages with different pasts answer the same programme with the same bits after reset()");
    }
    {
        BassChorus s;
        test::run (s.prepare (fs, 256, 2));
        Buf l = tone (24000, 997.0, fs, 0.5f), r = l;
        render (s, l, &r, 256);
        Buf m = tone (2560, 500.0, fs, 0.5f);
        render (s, m, nullptr, 256);             // the right channel stopped for these calls
        Buf zl (4800, 0.0f), zr (4800, 0.0f);
        render (s, zl, &zr, 256);
        double leak = 0.0;
        for (float v : zr)
            leak = std::max (leak, (double) std::fabs (v));
        ok (leak <= 0.0, "after mono calls, the stopped right channel's silence stays silent (leak " + num (leak) + ")");
    }
    namespace pz = test::law11c;
    int pausesEqual = 0, pausesMoved = 0, cases = 0;
    for (int gap : pz::kGaps)
    {
        // A one-sample pause would compare A with itself (B's silence would be empty): nothing to learn.
        if (gap < 2)
            continue;
        BassChorus A, B, C;
        for (auto* s : { &A, &B, &C })
        {
            s->setParams (chorus::paramsForAmount (0.6f, s->params()));
            test::run (s->prepare (fs, 512, 2));
            Buf l = tone (9600, 997.0, fs, 0.5f), r = l;
            render (*s, l, &r, 512);
        }
        // Drained in lockstep (all three get the same silence), to exact rest.
        bool rest = true;
        for (int k = 0; k < 100; ++k)
            for (auto* s : { &A, &B, &C })
            {
                Buf zl (512, 0.0f), zr (512, 0.0f);
                float* io[2] = { zl.data(), zr.data() };
                test::run (s->process (io, 2, 512));
                if (k == 99)
                    for (std::size_t i = 0; i < 512; ++i)
                        rest = rest && pz::bitsEqual (zl[i], 0.0f) && pz::bitsEqual (zr[i], 0.0f);
            }
        ok (rest, "precondition: drained to exact rest");
        for (auto* s : { &A, &B, &C })
        {
            s->setParams (chorus::paramsForAmount (0.9f, s->params()));   // the glides are in flight through the gap
            s->setRateHz (3.0f);
            s->setCrossoverHz (400.0f);
        }
        // A pauses for `gap`; B hears gap - 1 samples of silence and then pauses for one, so its lanes stop as A's
        // do (each lane's line and door are its own, dropped at the stop, law 11a) and only the SHARED clocks can
        // tell the two apart — which a pause must move exactly as silence does (law 11c).
        test::run (A.process (nullptr, 0, gap));
        Buf zl ((std::size_t) gap, 0.0f), zr = zl;
        float* zio[2] = { zl.data(), zr.data() };
        test::run (B.process (zio, 2, gap - 1));
        test::run (B.process (nullptr, 0, 1));
        std::vector<Buf> outs;
        for (auto* s : { &A, &B, &C })
        {
            Buf l = tone (4800, 997.0, fs, 0.5f), r = tone (4800, 1500.0, fs, 0.5f);
            render (*s, l, &r, 512);
            l.insert (l.end(), r.begin(), r.end());
            outs.push_back (l);
        }
        ++cases;
        pausesEqual += bitsEqual (outs[0], outs[1]) ? 1 : 0;
        pausesMoved += bitsEqual (outs[0], outs[2]) ? 0 : 1;
    }
    ok (pausesEqual == cases, "a pause of every length in kGaps leaves the stage where that much silence leaves it ("
                              + std::to_string (pausesEqual) + "/" + std::to_string (cases) + ")");
    ok (pausesMoved == cases, "...and it moved it: the tone after it differs from no gap at all ("
                              + std::to_string (pausesMoved) + "/" + std::to_string (cases) + ")");
}

//==============================================================================
void poisonAndLaw8()
{
    group ("poison: a bad input sample stays the dry's and is never voiced; an overflow heals by the next call");
    const double fs = 48000.0;
    {
        BassChorus a, b;
        a.setParams (chorus::paramsForAmount (1.0f, a.params()));
        b.setParams (chorus::paramsForAmount (1.0f, b.params()));
        test::run (a.prepare (fs, 512, 2));
        test::run (b.prepare (fs, 512, 2));
        Buf l = tone (24000, 997.0, fs, 0.5f), r = l;
        const int bad[] = { 5000, 5001, 9000 };
        const float val[] = { std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() };
        Buf lc = l, rc = r;
        for (int k = 0; k < 3; ++k)
        {
            l[(std::size_t) bad[k]] = val[k];
            lc[(std::size_t) bad[k]] = 0.0f;
        }
        render (a, l, &r, 512);
        render (b, lc, &rc, 512);
        bool sameElsewhere = bitsEqual (r, rc), badKept = true;
        for (std::size_t i = 0; i < l.size(); ++i)
        {
            const bool isBad = i == 5000 || i == 5001 || i == 9000;
            if (isBad)
                badKept = badKept && ! std::isfinite (l[i]);
            else
                sameElsewhere = sameElsewhere && core::sameBits (l[i], lc[i]);
        }
        ok (badKept, "the bad samples come out as they went in: the dry keeps them");
        ok (sameElsewhere, "every other sample is bit-identical to the same stream with those samples at 0: nothing was voiced");
    }
    {
        BassChorus s;
        s.setParams (chorus::paramsForAmount (1.0f, s.params()));
        test::run (s.prepare (fs, 512, 1));
        Buf l = tone (4096, 997.0, fs, 0.5f);
        render (s, l, nullptr, 512);
        Buf big (512, 0.0f);
        for (std::size_t i = 0; i < big.size(); ++i)
            big[i] = (i & 1u) ? -3.0e38f : 3.0e38f;
        render (s, big, nullptr, 512);
        Buf after = tone (48000, 997.0, fs, 0.5f);
        const Buf afterDry = after;
        render (s, after, nullptr, 512);
        bool noNaN = true;
        for (std::size_t i = 0; i < after.size(); ++i)
            noNaN = noNaN && ! (std::isfinite (afterDry[i]) && std::isnan (after[i]));
        ok (noNaN, "while the line replays the +-3e38 burst, no finite dry sample comes out as NaN (overflow is +-inf at worst)");
        // The line still replays what the burst wrote until the head has passed it; after that, nothing is left.
        const auto from = (std::size_t) (BassChorus::tailSeconds() * fs) + 512u;
        bool finite = true;
        for (std::size_t i = from; i < after.size(); ++i)
            finite = finite && std::isfinite (after[i]);
        ok (finite, "after a block of +-3e38 (the filters overflow), every sample past the line's reach is finite: no state was poisoned");
    }

    group ("law 8: a minute of silence after a chorus leaves no subnormal in the lines or the glides, and the silence is exact");
    for (double fs8 : kRates)
    {
        BassChorus s;
        s.setParams (chorus::paramsForAmount (1.0f, s.params()));
        test::run (s.prepare (fs8, 512, 2));
        Buf l = plus (tone ((int) (0.5 * fs8), 997.0, fs8, 0.5f), tone ((int) (0.5 * fs8), 41.2, fs8, 0.5f)), r = l;
        render (s, l, &r, 512);
        bool sub = false, silent = true;
        const int blocks = (int) (60.0 * fs8) / 512;
        for (int k = 0; k < blocks; ++k)
        {
            Buf zl (512, 0.0f), zr (512, 0.0f);
            float* io[2] = { zl.data(), zr.data() };
            test::run (s.process (io, 2, 512));
            sub = sub || s.anySubnormalState();
            if (k > blocks - (int) fs8 / 512)
                for (int i = 0; i < 512; ++i)
                    silent = silent && core::sameBits (zl[(std::size_t) i], 0.0f) && core::sameBits (zr[(std::size_t) i], 0.0f);
        }
        ok (! sub, "no subnormal in the lines or the glides at any block of the minute, " + khz (fs8));
        ok (silent, "the last second is exact zeros, " + khz (fs8));
    }
}


//==============================================================================
// Review round 1.

void depthCapCrossing()
{
    group ("the depth cap: a centre gliding into and out of the swing's limit bends the head with no pitch step");
    for (double fs : { 48000.0, 96000.0 })
    {
        const double f = 997.0;
        const float  A = 0.5f;
        const double w = 2.0 * kPi * f / fs;
        // Dry plus a voice up to 30 % sharp: 0.5 Hz x 10 ms on the LFO, plus a centre gliding ~21 ms at up to
        // ~0.1 sample per sample, doubled while the cap makes the swing follow the centre.
        const double bound = 1.25 * (double) A * (std::pow (2.0 * std::sin (0.5 * w), 3.0) + std::pow (2.0 * std::sin (0.5 * w * 1.3), 3.0));
        double worst = 0.0;
        for (int trial = 0; trial < 8; ++trial)   // the crossings land at eight different LFO phases
        {
            BassChorus s;
            s.setParams ({ BassChorus::kMaxDepthMs, 1.0f, 0.5f, 200.0f, BassChorus::kMaxDelayMs });
            test::run (s.prepare (fs, 64, 1));
            const int down = ((int) ((0.3 + 0.13 * trial) * fs) / 64) * 64, up = down + ((int) (0.8 * fs) / 64) * 64;
            const int N = up + (int) (0.8 * fs);
            Buf y = tone (N, f, fs, A);
            for (int pos = 0; pos < N; pos += 64)
            {
                if (pos == down) s.setDelayMs (BassChorus::kMinDelayMs);   // into the cap
                if (pos == up)   s.setDelayMs (BassChorus::kMaxDelayMs);   // and out of it
                float* io[1] = { y.data() + pos };
                test::run (s.process (io, 1, std::min (64, N - pos)));
            }
            worst = std::max (worst, maxThirdDiff (y, down));
        }
        ok (worst <= bound, "centre 25 -> min -> 25 ms with depth 10 ms: largest third difference " + num (worst, 3)
                            + " <= the smooth bound " + num (bound, 3) + ", " + khz (fs));
    }
}

void returningLane()
{
    group ("a lane that comes back after a narrow stretch or a gap fades its voice in: no step 10 ms later");
    const double fs = 48000.0, f = 1000.0;
    const double steady = 0.5 * 2.0 * std::sin (kPi * f / fs);   // the steady voice's own largest step
    for (int how = 0; how < 2; ++how)
    {
        BassChorus s;
        s.setParams ({ 0.0f, 1.0f, 0.7f, 200.0f, 10.0f });
        test::run (s.prepare (fs, 64, 2));
        Buf l = tone (14400, f, fs, 0.5f), r = l;
        render (s, l, &r, 64);
        if (how == 0)
        {
            Buf m = tone (9600, f, fs, 0.5f);
            render (s, m, nullptr, 64);                  // the right lane stops for 0.2 s
        }
        else
            test::run (s.process (nullptr, 0, 9600));    // both stop: a 0.2 s gap

        // The right lane comes back on a cosine crest (its first sample is 0.5). After mono calls the left one
        // never stopped and carries on its tone without a seam (24000 samples are 500 whole cycles); after a gap
        // it comes back on the crest too.
        Buf rr = tone (4800, f, fs, 0.5f, 0.5 * kPi), rl = how == 0 ? tone (4800, f, fs, 0.5f) : rr;
        const Buf dryL = rl, dryR = rr;
        render (s, rl, &rr, 64);
        double worstL = 0.0, worstR = 0.0;
        for (std::size_t i = 1; i < dryR.size(); ++i)
        {
            const auto step = [i] (const Buf& y, const Buf& d)
            { return std::fabs (((double) y[i] - d[i]) - ((double) y[i - 1] - d[i - 1])); };
            worstL = std::max (worstL, step (rl, dryL));
            worstR = std::max (worstR, step (rr, dryR));
        }
        const std::string what = how == 0 ? "after mono calls" : "after a zero-width gap";
        ok (worstR <= 1.05 * steady, "the returning right voice's largest step " + num (worstR, 3) + " <= its steady "
                                     + num (steady, 3) + " x 1.05, " + what);
        ok (worstL <= 1.05 * steady, "the left voice's largest step " + num (worstL, 3) + ", " + what);
    }
}

void tailBound()
{
    group ("tailSeconds(): the head's reach plus the high pass's ring, down to -100 dB of the input's peak");
    for (double fs : kRates)
        for (int glide = 0; glide < 2; ++glide)
        {
            const float lo = BassChorus::kMinCrossoverHz;
            BassChorus s;
            s.setParams ({ BassChorus::kMaxDepthMs, 1.0f, 0.5f, glide != 0 ? 1000.0f : lo, BassChorus::kMaxDelayMs });
            test::run (s.prepare (fs, 512, 2));
            const int N = (int) fs;
            Buf l = plus (tone (N, (double) lo, fs, 0.5f), noise (N, 31u, 0.4f)), r = plus (tone (N, (double) lo, fs, 0.5f, 1.0), noise (N, 32u, 0.4f));
            double peak = 0.0;
            for (std::size_t i = 0; i < l.size(); ++i)
                peak = std::max ({ peak, (double) std::fabs (l[i]), (double) std::fabs (r[i]) });
            render (s, l, &r, 512);
            if (glide != 0)
                s.setCrossoverHz (lo);   // the corner glides down as the input stops
            const int T = (int) std::ceil ((double) BassChorus::tailSeconds() * fs);
            Buf zl ((std::size_t) (T + 8192), 0.0f), zr = zl;
            render (s, zl, &zr, 512);
            double after = 0.0;
            for (std::size_t i = (std::size_t) T; i < zl.size(); ++i)
                after = std::max ({ after, (double) std::fabs (zl[i]), (double) std::fabs (zr[i]) });
            ok (after <= 1.0e-5 * peak, "after " + num (1000.0 * BassChorus::tailSeconds(), 3) + " ms of silence the voice is "
                                        + num (20.0 * std::log10 (std::max (after, 1e-30) / peak), 3) + " dB re the input's peak (<= -100), "
                                        + (glide != 0 ? "crossover gliding 1000 -> " : "crossover ") + num (lo) + " Hz, " + khz (fs));
        }
}

double lr4HighDb (double f, double fc)
{
    const double r4 = std::pow (f / fc, 4.0);
    return 20.0 * std::log10 (r4 / (1.0 + r4));
}

/** Side = (L - R)/2, as energy re the dry's, over seconds 1 … 4 of a lone tone at `f` on both sides at amount 1:
    every sideband the swing makes is in it, because nothing else is. */
double sideDb (double fs, double f, float xover, float rate, float centre)
{
    BassChorus s;
    BassChorusParams menu;
    menu.crossoverHz = xover;
    menu.rateHz = rate;
    menu.delayMs = centre;
    s.setParams (chorus::paramsForAmount (1.0f, menu));
    test::run (s.prepare (fs, 512, 2));
    const int N = (int) (4.0 * fs);
    const Buf dry = tone (N, f, fs, 0.5f);
    Buf l = dry, r = dry;
    render (s, l, &r, 512);
    double es = 0.0, ed = 0.0;
    for (int i = (int) fs; i < N; ++i)
    {
        const double side = 0.5 * ((double) l[(std::size_t) i] - (double) r[(std::size_t) i]);
        es += side * side;
        ed += (double) dry[(std::size_t) i] * (double) dry[(std::size_t) i];
    }
    return 10.0 * std::log10 (std::max (es, 1e-300) / ed);
}

void lowBandSide()
{
    group ("Side = (L - R)/2 of a tone under the crossover: at most the voice's own LR4 attenuation there, any setting");
    const double fs = 48000.0;
    struct Set { float rate, centre; };
    int cases = 0, held = 0;
    double worstMargin = -1e9;
    for (const float xover : { BassChorus::kMinCrossoverHz, 200.0f, BassChorus::kMaxCrossoverHz })
        for (const Set m : { Set { BassChorus::kMaxRateHz, BassChorus::kMaxDelayMs }, Set { BassChorus::kMaxRateHz, BassChorus::kMinDelayMs },
                             Set { BassChorus::kMinRateHz, BassChorus::kMaxDelayMs } })
            for (const double f : { 0.499 * xover, 0.25 * xover })   // the band's upper edge, and an octave under it
            {
                const double got = sideDb (fs, f, xover, m.rate, m.centre), bound = lr4HighDb (f, xover);
                ++cases;
                held += got <= bound + 0.2 ? 1 : 0;
                worstMargin = std::max (worstMargin, got - bound);
            }
    ok (held == cases, "Side <= the LR4 high half's level at the tone (+0.2 dB) in " + std::to_string (held) + "/" + std::to_string (cases)
                       + " cases: crossover 100 / 200 / 1000 Hz, rate and centre at their extremes, tones at 0.499 and 0.25 of the crossover"
                       + " (closest: " + num (worstMargin, 3) + " dB from the bound)");
    const double edge = sideDb (fs, 499.0, BassChorus::kMaxCrossoverHz, BassChorus::kMaxRateHz, BassChorus::kMaxDelayMs);
    ok (edge <= lr4HighDb (499.0, 1000.0) + 0.2, "the band's upper edge at the top crossover — 499 Hz under 1 kHz, 10 Hz, 25 ms: "
                                                 + num (edge, 3) + " dB (bound " + num (lr4HighDb (499.0, 1000.0), 3) + ")");

    // What the guarantee is FOR: the bass's open-string fundamentals, at the default crossover and at its floor.
    for (const float xover : { 200.0f, BassChorus::kMinCrossoverHz })
        for (const double f : { 30.87, 41.20, 55.00, 73.42, 98.00 })
        {
            const double got = sideDb (fs, f, xover, BassChorus::kMaxRateHz, BassChorus::kMaxDelayMs), bound = lr4HighDb (f, (double) xover);
            ok (got <= bound + 0.2, "open string at " + num (f) + " Hz, crossover " + num (xover) + " Hz: Side " + num (got, 3)
                                    + " dB re the dry (bound " + num (bound, 3) + ")");
        }
}

void fullDepthAtFloor()
{
    group ("the centre's floor leaves the slider's full depth whole: 3 ms at amount 1 on the lowest centre");
    BassChorusParams lo = chorus::paramsForAmount (1.0f), mid = lo;
    lo.rateHz = mid.rateHz = BassChorus::kMaxRateHz;
    lo.delayMs = BassChorus::kMinDelayMs;
    mid.delayMs = 12.0f;
    approx (chorus::pitchDeviationCents (lo), chorus::pitchDeviationCents (mid), 1e-3,
            "the read-out at 10 Hz on the lowest centre is the one on a 12 ms centre (" + num (chorus::pitchDeviationCents (mid)) + " cents)");

    const double fs = 48000.0, want = 2.0 * kPi * 2.0 * 0.003;
    BassChorusParams p = chorus::paramsForAmount (1.0f);
    p.rateHz = 2.0f;
    p.delayMs = BassChorus::kMinDelayMs;
    BassChorus s;
    s.setParams (p);
    test::run (s.prepare (fs, 512, 1));
    const int N = (int) (1.4 * fs);
    const Buf x = tone (N, 1000.0, fs, 0.5f);
    Buf y = x;
    render (s, y, nullptr, 512);
    Buf v ((std::size_t) N);
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = y[i] - x[i];
    const auto d = deviation (v, (int) (0.3 * fs), fs, 1000.0);
    const double peak = *std::max_element (d.begin(), d.end());
    approx (peak, want, 0.02 * want, "measured on the stage at 2 Hz on the lowest centre: +" + num (100.0 * peak, 4) + " % (set "
                                     + num (100.0 * want, 4) + " %, the full 3 ms)");
}

void nanAmountIgnored()
{
    group ("a non-finite amount changes nothing when the caller passes the stage's own params() as the menu");
    BassChorus a, b;
    for (auto* st : { &a, &b })
    {
        st->setParams (chorus::paramsForAmount (0.9f, st->params()));
        test::run (st->prepare (48000.0, 512, 2));
    }
    a.setParams (chorus::paramsForAmount (std::numeric_limits<float>::quiet_NaN(), a.params()));
    a.setParams (chorus::paramsForAmount (std::numeric_limits<float>::infinity(), a.params()));
    const auto pa = a.params(), pb = b.params();
    ok (core::sameBits (pa.depthMs, pb.depthMs) && core::sameBits (pa.mix, pb.mix) && core::sameBits (pa.rateHz, pb.rateHz)
            && core::sameBits (pa.crossoverHz, pb.crossoverHz) && core::sameBits (pa.delayMs, pb.delayMs) && core::sameBits (pa.mix, 0.9f),
        "after amount 0.9, NaN and +inf leave every knob where it was (mix " + num (pa.mix) + ", depth " + num (pa.depthMs) + " ms)");
    Buf l1 = plus (tone (20000, 997.0, 48000.0, 0.5f), noise (20000, 51u, 0.1f)), r1 = noise (20000, 52u, 0.3f), l2 = l1, r2 = r1;
    render (a, l1, &r1, 512);
    render (b, l2, &r2, 512);
    ok (bitsEqual (l1, l2) && bitsEqual (r1, r2), "...and the sound is the same bits as a stage never handed them");
}

void idleStopsAtTheLanding()
{
    group ("amount 0 in mid-call: the voice stops on the sample the mix ramp lands, the rest of the call is untouched");
    for (double fs : kRates)
    {
        BassChorus s;
        s.setParams (chorus::paramsForAmount (0.8f, s.params()));
        test::run (s.prepare (fs, 512, 2));
        Buf l = tone ((int) (0.2 * fs), 997.0, fs, 0.5f), r = l;
        render (s, l, &r, 512);
        const auto before = s.voicedSamples();
        s.setParams (chorus::paramsForAmount (0.0f, s.params()));
        const int N = (int) (10.0 * fs), ramp = (int) (BassChorus::kMixRampSeconds * fs + 0.5);
        Buf x = plus (tone (N, 997.0, fs, 0.5f), noise (N, 61u, 0.1f)), y = noise (N, 62u, 0.3f);
        const Buf x0 = x, y0 = y;
        float* io[2] = { x.data(), y.data() };
        test::run (s.process (io, 2, N));    // one call, ten seconds
        const auto ran = s.voicedSamples() - before;
        ok (ran < (std::uint64_t) ramp, "voiced " + std::to_string (ran) + " of " + std::to_string (N) + " samples (the ramp is "
                                        + std::to_string (ramp) + "), " + khz (fs));
        bool untouched = s.isIdle();
        for (std::size_t i = (std::size_t) ramp; i < (std::size_t) N; ++i)
            untouched = untouched && core::sameBits (x[i], x0[i]) && core::sameBits (y[i], y0[i]);
        ok (untouched, "...idle, and every sample from the landing on is the input's own, " + khz (fs));
    }
}

//==============================================================================
void noAllocation()
{
    group ("RT: process, reset and the setters never allocate");
    BassChorus s;
    test::run (s.prepare (96000.0, 256, 2));
    Buf l (1024, 0.25f), r (1024, -0.25f);
    float* lr[2] = { l.data(), r.data() };
    const long long before = alloc::count.load();
    s.setParams (chorus::paramsForAmount (0.7f, s.params()));
    s.setParams ({ 3.0f, 0.9f, 2.0f, 300.0f, 12.0f });
    test::run (s.process (lr, 2, 1024));
    test::run (s.process (lr, 1, 1024));
    test::run (s.process (nullptr, 0, 1024));
    s.setParams (chorus::paramsForAmount (0.0f, s.params()));
    for (int k = 0; k < 200; ++k)
        test::run (s.process (lr, 2, 1024));   // down to idle
    s.setParams (chorus::paramsForAmount (0.5f, s.params()));
    test::run (s.process (lr, 2, 1024));       // and the wake
    s.reset();
    s.reset();
    test::run (s.process (lr, 2, 1024));
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation: setters, process (stereo, mono, pause, idle, wake), reset");
}
} // namespace

int main()
{
    std::printf ("felitronics::chorus::BassChorus\n");
    contract();
    theLaw();
    amountZero();
    splitAndSumFlat();
    lowsMonoInPhase();
    pitchDeviation();
    noClicks();
    reslicing();
    resetAndStops();
    poisonAndLaw8();
    depthCapCrossing();
    returningLane();
    tailBound();
    lowBandSide();
    fullDepthAtFloor();
    nanAmountIgnored();
    idleStopsAtTheLanding();
    noAllocation();
    return felitronics::test::report();
}
