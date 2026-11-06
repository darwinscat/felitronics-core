// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::cabinet::CabinetIr — an amp's cabinet stage (moved from OrbitAmp, its DSP half).
//
//   a shot loads and lands: after the convolver's 50 ms crossfade the output is the IR's convolution at the level the
//   reference asks — the guitar's or the bass's — at sample zero; the post (a trimmed tail, a cut, the sign) is baked
//   into the taps, the level taken after it; a shot at another rate is resampled and keeps its level; a NaN in leaves
//   every later sample a number; out of the path it drains, and comes back with no ghost of what it played before;
//   a shot `takes` refuses changes nothing; and process() / idle() allocate nothing.

#include <felitronics_test.h>
#include <alloc_counter.h>
#include <felitronics/cabinet/CabinetIr.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using namespace felitronics;
using cabinet::CabinetIr;

namespace
{
    constexpr int bs = 256;

    /** A small cabinet-like IR at `rate`: a decaying, band-limited resonance with a zero-mean wiggle. */
    std::vector<float> shot (double rate, double seconds = 0.08)
    {
        std::vector<float> ir ((size_t) std::llround (seconds * rate));
        for (size_t i = 0; i < ir.size(); ++i)
        {
            const double t = (double) i / rate;
            ir[i] = (float) (std::exp (-t * 60.0) * (0.6 * std::cos (2.0 * core::kPi * 110.0 * t) + 0.3 * std::sin (2.0 * core::kPi * 1800.0 * t)));
        }
        ir[0] += 0.5f;
        return ir;
    }

    /** Blocks of a 0.1-peak 220 Hz tone through the stage, `count` of them; the last block's left channel kept. */
    std::vector<float> play (CabinetIr& cab, int count, double rate, long long& phase, bool on = true)
    {
        std::vector<float> l (bs), r (bs);
        for (int b = 0; b < count; ++b)
        {
            for (int i = 0; i < bs; ++i, ++phase)
                l[(size_t) i] = r[(size_t) i] = (float) (0.1 * std::sin (2.0 * core::kPi * 220.0 * (double) phase / rate));
            float* io[] { l.data(), r.data() };
            if (on)
                cab.process (io, 2, bs, true);
            else
                cab.idle (bs);
            cab.flushPending();
        }
        return l;
    }

    double rms (const std::vector<float>& x)
    {
        double s = 0.0;
        for (const float v : x)
            s += (double) v * v;
        return std::sqrt (s / (double) x.size());
    }
}

int main()
{
    std::printf ("felitronics::cabinet CabinetIr tests\n");

    test::group ("a shot lands, at the level its reference asks, at sample zero");
    {
        constexpr double fs = 48000.0;
        auto ir = shot (fs);
        for (const auto ref : { CabinetIr::Reference::guitar, CabinetIr::Reference::bass })
        {
            CabinetIr cab;
            cab.prepare (fs, bs, 2);
            cab.setReference (ref);
            test::ok (cab.load (ir.data(), (int) ir.size(), fs), "the shot is taken");
            long long phase = 0;
            const auto y = play (cab, 40, fs, phase);   // past the 50 ms crossfade and the IR's length

            const float g = ref == CabinetIr::Reference::bass
                              ? measurement::bassReferenceUnityGain (ir.data(), (int) ir.size(), fs)
                              : measurement::referenceUnityGain (ir.data(), (int) ir.size(), fs);
            // The direct convolution of the same tone with the IR times that gain, at the same samples.
            double worst = 0.0;
            const long long end = phase;
            for (int i = 0; i < bs; ++i)
            {
                const long long n = end - bs + i;
                double acc = 0.0;
                for (size_t k = 0; k < ir.size() && (long long) k <= n; ++k)
                    acc += (double) ir[k] * 0.1 * std::sin (2.0 * core::kPi * 220.0 * (double) (n - (long long) k) / fs);
                worst = std::max (worst, std::abs ((double) y[(size_t) i] - acc * (double) g));
            }
            test::ok (worst < 2.0e-4, std::string (ref == CabinetIr::Reference::bass ? "bass" : "guitar")
                                      + "'s reference: the convolution times its gain (" + std::to_string (worst) + " apart)");
            test::approx (cab.tailSeconds(), 0.08, 1.0e-3, "the tail it reports is the shot's");
        }
    }

    test::group ("the post is baked in; the level is taken after it");
    {
        constexpr double fs = 48000.0;
        auto ir = shot (fs, 0.4);
        CabinetIr cab;
        cab.prepare (fs, bs, 2);
        CabinetIr::Post trimmed;
        trimmed.trimOn = true;
        trimmed.trimFraction = 0.25f;
        cab.setPost (trimmed);
        test::ok (cab.load (ir.data(), (int) ir.size(), fs), "taken");
        test::approx (cab.tailSeconds(), 0.1, 1.0e-3, "a quarter of the tail kept");
        test::ok (cab.currentPost() == trimmed, "the post it stands at");
        CabinetIr::Post flipped;
        flipped.phase = true;
        cab.setPost (flipped);
        long long phase = 0;
        const auto a = play (cab, 40, fs, phase);
        CabinetIr plain;
        plain.prepare (fs, bs, 2);
        (void) plain.load (ir.data(), (int) ir.size(), fs);
        long long phase2 = 0;
        const auto b = play (plain, 40, fs, phase2);
        double worst = 0.0;
        for (int i = 0; i < bs; ++i)
            worst = std::max (worst, std::abs ((double) a[(size_t) i] + (double) b[(size_t) i]));
        test::ok (worst < 1.0e-6, "the sign flipped: the same cabinet, negated (" + std::to_string (worst) + ")");
    }

    test::group ("a shot at another rate keeps its level");
    {
        const auto ir = shot (48000.0);
        for (const double host : { 44100.0, 96000.0 })
        {
            CabinetIr at48, atHost;
            at48.prepare (48000.0, bs, 2);
            atHost.prepare (host, bs, 2);
            at48.setReference (CabinetIr::Reference::bass);
            atHost.setReference (CabinetIr::Reference::bass);
            (void) at48.load (ir.data(), (int) ir.size(), 48000.0);
            (void) atHost.load (ir.data(), (int) ir.size(), 48000.0);
            long long p1 = 0, p2 = 0;
            (void) play (at48, 40, 48000.0, p1);
            (void) play (atHost, 40, host, p2);
            double a = 0.0, b = 0.0;
            for (int k = 0; k < 80; ++k)
            {
                a += std::pow (rms (play (at48, 1, 48000.0, p1)), 2.0);
                b += std::pow (rms (play (atHost, 1, host, p2)), 2.0);
            }
            a = std::sqrt (a);
            b = std::sqrt (b);
            test::approx (20.0 * std::log10 (b / a), 0.0, 0.1, "48 kHz shot, " + std::to_string ((int) host) + " Hz session: the same level");
        }
    }

    test::group ("a NaN in: every later sample a number");
    {
        constexpr double fs = 48000.0;
        const auto ir = shot (fs);
        CabinetIr cab;
        cab.prepare (fs, bs, 2);
        (void) cab.load (ir.data(), (int) ir.size(), fs);
        long long phase = 0;
        (void) play (cab, 20, fs, phase);
        std::vector<float> l (bs, 0.05f), r (bs, 0.05f);
        l[17] = r[17] = std::numeric_limits<float>::quiet_NaN();
        float* io[] { l.data(), r.data() };
        cab.process (io, 2, bs, true);
        bool finite = true;
        for (int b = 0; b < 40; ++b)
        {
            const auto y = play (cab, 1, fs, phase);
            for (const float v : y)
                finite = finite && std::isfinite (v);
        }
        test::ok (finite, "after the NaN: numbers, block for block");
        test::ok (rms (play (cab, 1, fs, phase)) > 1.0e-3, "and the cabinet still sounds");
    }

    test::group ("out of the path it drains; back, no ghost");
    {
        constexpr double fs = 48000.0;
        const auto ir = shot (fs);
        CabinetIr cab;
        cab.prepare (fs, bs, 2);
        (void) cab.load (ir.data(), (int) ir.size(), fs);
        long long phase = 0;
        (void) play (cab, 40, fs, phase);
        const int drain = (int) std::ceil ((CabinetIr::maxSeconds + 1.0) * fs / bs);
        (void) play (cab, drain, fs, phase, false);   // out of the path: drained
        std::vector<float> l (bs, 0.0f), r (bs, 0.0f);
        float* io[] { l.data(), r.data() };
        cab.process (io, 2, bs, true);
        double peak = 0.0;
        for (const float v : l)
            peak = std::max (peak, (double) std::abs (v));
        test::ok (peak < 1.0e-9, "silence in after the drain: silence out (" + std::to_string (peak) + ")");
    }

    test::group ("a shot it does not take changes nothing");
    {
        test::ok (! CabinetIr::takes (0, 48000.0) && ! CabinetIr::takes (100, 4000.0) && ! CabinetIr::takes (100, 1.0e6)
                      && ! CabinetIr::takes ((long long) (6.0 * 48000.0), 48000.0) && CabinetIr::takes (24000, 48000.0),
                  "takes(): a length, a rate in 8 kHz … 768 kHz, at most 5 s");
        CabinetIr cab;
        cab.prepare (48000.0, bs, 2);
        const auto ir = shot (48000.0);
        (void) cab.load (ir.data(), (int) ir.size(), 48000.0);
        std::vector<float> long6 ((size_t) (6.0 * 48000.0), 0.1f);
        test::ok (! cab.load (long6.data(), (int) long6.size(), 48000.0) && ! cab.load (nullptr, 10, 48000.0), "refused");
        test::approx (cab.tailSeconds(), 0.08, 1.0e-3, "the cabinet it had stays");
    }

    test::group ("no allocation in process() and idle()");
    {
        constexpr double fs = 48000.0;
        const auto ir = shot (fs);
        CabinetIr cab;
        cab.prepare (fs, bs, 2);
        (void) cab.load (ir.data(), (int) ir.size(), fs);
        long long phase = 0;
        (void) play (cab, 10, fs, phase);
        std::vector<float> l (bs, 0.1f), r (bs, 0.1f);
        float* io[] { l.data(), r.data() };
        float* mono[] { l.data() };
        const long long before = alloc::count.load();
        for (int b = 0; b < 50; ++b)
        {
            cab.process (b % 2 == 0 ? io : mono, b % 2 == 0 ? 2 : 1, bs, true);
            cab.idle (bs);
        }
        test::okNoAlloc (alloc::count.load() == before, "stereo, mono and idle calls");
    }

    return test::report();
}
