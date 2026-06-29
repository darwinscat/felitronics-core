// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko <oleh@darwinscat.com> & Alisa <alisa@darwinscat.com>. Part of TabbyEQ — see LICENSE.

// Verifies the runtime engine by actually pushing audio through it: measured steady-state gain
// of a sine must agree with the analytic magnitudeDb() readout, smoothing converges, the swept
// SVF band matches the matched band, a 24 dB/oct HP rolls off ~24 dB/oct, and process() is clean
// (silence stays silent, no NaN). JUCE-free.

#include "TestUtil.h"

#include <teq/EqBand.h>
#include <teq/EqEngine.h>
#include <teq/Smoother.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace teq;

namespace
{
    // Steady-state gain (dB) of a mono in-place processing callback at frequency f.
    template <class Apply>
    double sineGainDb (Apply&& apply, double f, double fs, int warm = 600)
    {
        const int block = 256;                 // warm=600 -> ~3.2 s >> 30 ms smoothing + filter ring
        const double amp = 0.25, dp = 2.0 * kPi * f / fs;
        std::vector<float> buf ((size_t) block);
        double phase = 0.0;
        auto fill = [&] { for (int n = 0; n < block; ++n) { buf[(size_t) n] = (float) (amp * std::sin (phase));
                                                            phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi; } };

        for (int b = 0; b < warm; ++b) { fill(); float* p = buf.data(); apply (&p, 1, block); }

        fill();
        std::vector<float> in (buf.begin(), buf.end());
        { float* p = buf.data(); apply (&p, 1, block); }

        double ir = 0.0, orr = 0.0;
        for (int n = 0; n < block; ++n)
        {
            ir  += (double) in[(size_t) n]  * in[(size_t) n];
            orr += (double) buf[(size_t) n] * buf[(size_t) n];
        }
        return 20.0 * std::log10 (std::sqrt (orr / std::max (1e-30, ir)));
    }

    // Per-channel steady-state gain (dB) for an N-channel in-place callback. Every channel gets the
    // SAME sine but a DISTINCT amplitude, so any cross-channel state bleed shows up as a wrong gain.
    template <class Apply>
    std::vector<double> multiSineGainDb (Apply&& apply, int C, double f, double fs, int warm = 600)
    {
        const int block = 256;
        const double dp = 2.0 * kPi * f / fs;
        std::vector<std::vector<float>> bufs ((size_t) C, std::vector<float> ((size_t) block));
        std::vector<float*> ptr ((size_t) C);
        for (int c = 0; c < C; ++c) ptr[(size_t) c] = bufs[(size_t) c].data();

        double phase = 0.0;
        auto fill = [&]
        {
            for (int n = 0; n < block; ++n)
            {
                const double s = std::sin (phase);
                phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi;
                for (int c = 0; c < C; ++c) bufs[(size_t) c][(size_t) n] = (float) ((0.1 + 0.04 * c) * s);
            }
        };

        for (int b = 0; b < warm; ++b) { fill(); apply (ptr.data(), C, block); }

        fill();
        const std::vector<std::vector<float>> in = bufs;     // snapshot the inputs before in-place process
        apply (ptr.data(), C, block);

        std::vector<double> g ((size_t) C);
        for (int c = 0; c < C; ++c)
        {
            double ir = 0.0, orr = 0.0;
            for (int n = 0; n < block; ++n)
            {
                ir  += (double) in[(size_t) c][(size_t) n]   * in[(size_t) c][(size_t) n];
                orr += (double) bufs[(size_t) c][(size_t) n] * bufs[(size_t) c][(size_t) n];
            }
            g[(size_t) c] = 20.0 * std::log10 (std::sqrt (orr / std::max (1e-30, ir)));
        }
        return g;
    }

    bool anyNaN (const float* d, int n) { for (int i = 0; i < n; ++i) if (! std::isfinite (d[i])) return true; return false; }
}

void runEqEngineTests()
{
    using namespace teqtest;
    const double fs = 48000.0;

    group ("Smoother: converges to target, snap settles");
    {
        Smoother s; s.prepare (fs, 10.0); s.snap (0.0); s.setTarget (1.0);
        s.advance (48000);
        expectNear (s.value(), 1.0, 1e-3, "advance(1s) converges");
        s.snap (5.0);
        expectTrue (s.settled(), "snap settles");
    }

    group ("EqBand bell: measured gain == response() at centre, ~unity far off");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 6.0;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        const double resp = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 1000.0 / fs)));
        std::printf ("      bell @1k: measured=%.3f dB  response=%.3f dB\n", meas, resp);
        expectNear (meas, 6.0, 0.4,  "measured +6 dB");
        expectNear (resp, 6.0, 0.05, "response +6 dB");
        expectNear (meas, resp, 0.4, "measured ~ response");

        EqBand far; far.prepare (fs, 1);
        BandParams q; q.on = true; q.type = FilterType::Bell; q.freq = 1000.0; q.Q = 4.0; q.gainDb = 9.0;
        far.setParams (q);
        const double low = sineGainDb ([&] (float* const* ch, int nc, int n) { far.processBlock (ch, nc, n); }, 60.0, fs);
        expectNear (low, 0.0, 0.3, "unity an octave+ away");
    }

    group ("EqBand swept (SVF) bell tracks the same +6 dB");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 2000.0; p.Q = 3.0; p.gainDb = 6.0; p.swept = true;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 2000.0, fs);
        std::printf ("      swept bell @2k: measured=%.3f dB\n", meas);
        expectNear (meas, 6.0, 0.5, "SVF bell +6 dB");
    }

    group ("EqBand HP 24 dB/oct rolls off ~24 dB/octave in the stopband");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::HighPass; p.freq = 1000.0; p.Q = 0.707; p.slope = 24;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double g125 = sineGainDb (ap, 125.0, fs);
        const double g250 = sineGainDb (ap, 250.0, fs);
        std::printf ("      HP24: g(125)=%.2f  g(250)=%.2f  slope=%.1f dB/oct\n", g125, g250, g250 - g125);
        expectTrue ((g250 - g125) > 20.0 && (g250 - g125) < 28.0, "~24 dB/oct");
    }

    group ("EqEngine: two bells cascade, measured == magnitudeDb");
    {
        EqEngine eng; eng.prepare (fs, 512, 1);
        BandParams a; a.on = true; a.type = FilterType::Bell; a.freq = 300.0;  a.Q = 1.0; a.gainDb =  4.0;
        BandParams b; b.on = true; b.type = FilterType::Bell; b.freq = 3000.0; b.Q = 1.0; b.gainDb = -3.0;
        eng.setBand (0, a); eng.setBand (1, b);
        const double m = sineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, 300.0, fs);
        const double r = eng.magnitudeDb (300.0);
        std::printf ("      @300: measured=%.3f dB  magnitudeDb=%.3f dB\n", m, r);
        expectNear (m, r, 0.4, "measured ~ magnitudeDb");
        expectNear (r, 4.0, 0.6, "~+4 dB at band-0 centre");
    }

    group ("swept SVF BandPass has unity gain at centre");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::BandPass; p.freq = 1000.0; p.Q = 4.0; p.swept = true;
        band.setParams (p);
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        std::printf ("      swept BP @1k Q=4: measured=%.3f dB\n", meas);
        expectNear (meas, 0.0, 0.5, "unity (0 dB) at centre");
    }

    group ("swept SVF high shelf is Butterworth (plateau == gain, Q-independent)");
    {
        auto plateau = [&] (double Q)
        {
            EqBand band; band.prepare (fs, 1);
            BandParams p; p.on = true; p.type = FilterType::HighShelf; p.freq = 3000.0; p.Q = Q; p.gainDb = 12.0; p.swept = true;
            band.setParams (p);
            return sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 14000.0, fs);
        };
        const double g1 = plateau (1.0), g4 = plateau (4.0);
        std::printf ("      swept HS plateau @14k: Q=1 -> %.3f dB, Q=4 -> %.3f dB\n", g1, g4);
        expectNear (g1, 12.0, 0.6, "plateau ~ +12 dB");
        expectNear (g1, g4,   0.2, "Q ignored for shelves");
    }

    group ("parameter safety: bad Q / huge gain / NaN / Inf do not produce NaN");
    {
        const double inf = std::numeric_limits<double>::infinity();
        EqEngine eng; eng.prepare (fs, 256, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell;
        p.freq = std::nan (""); p.Q = 0.0; p.gainDb = inf;       // NaN freq, Q=0, +Inf gain
        eng.setBand (0, p);
        const int n = 256;
        std::vector<float> buf ((size_t) n);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.2f * (float) std::sin (2.0 * kPi * 1000.0 * i / fs);
        float* ch[1] = { buf.data() };
        for (int b = 0; b < 20; ++b) eng.process (ch, 1, n);
        expectTrue (! anyNaN (buf.data(), n),                                       "sanitised params -> finite output");
        expectTrue (std::isfinite (eng.magnitudeDb (1000.0)),                       "magnitudeDb finite");
        expectTrue (std::isfinite (EqEngine::magnitudeDbFor (&p, 1, 1000.0, fs)),   "magnitudeDbFor finite on bad params");
    }

    group ("denormal flush: tail decays to exact zero after silence");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 6.0;
        band.setParams (p);
        const int n = 64;
        std::vector<float> buf ((size_t) n, 0.0f);
        buf[0] = 1.0f;                                           // impulse, then silence
        float* ch[1] = { buf.data() };
        band.processBlock (ch, 1, n);
        for (int b = 0; b < 200; ++b) { for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n); }
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "tail flushed to exact zero (no denormal residue)");
    }

    group ("first setParams snaps (no ramp from defaults on load)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 6.0;
        band.setParams (p);
        // measured after only 2 blocks (< 30 ms smoothing): snapped -> already +6 dB, not mid-ramp.
        const double meas = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs, 2);
        std::printf ("      gain after 2 blocks: %.3f dB\n", meas);
        expectNear (meas, 6.0, 0.5, "target reached immediately (no ramp from defaults)");
    }

    group ("stateless magnitudeDbFor matches the live magnitudeDb");
    {
        EqEngine eng; eng.prepare (fs, 512, 1);
        BandParams arr[EqEngine::kMaxBands] {};
        arr[0].on = true; arr[0].type = FilterType::Bell;      arr[0].freq = 300.0;  arr[0].Q = 1.0; arr[0].gainDb =  4.0;
        arr[1].on = true; arr[1].type = FilterType::HighShelf; arr[1].freq = 8000.0;                 arr[1].gainDb = -5.0;
        eng.setBand (0, arr[0]); eng.setBand (1, arr[1]);
        const int n = 64; std::vector<float> z ((size_t) n, 0.0f); float* ch[1] = { z.data() };
        for (int b = 0; b < 4; ++b) eng.process (ch, 1, n);     // let the live coeffs update
        for (double f : { 100.0, 300.0, 1000.0, 8000.0, 15000.0 })
        {
            const double live      = eng.magnitudeDb (f);
            const double stateless = EqEngine::magnitudeDbFor (arr, EqEngine::kMaxBands, f, fs);
            expectNear (stateless, live, 1e-6, "stateless == live @" + std::to_string ((int) f));
        }
    }

    group ("swept HP slope=24 is single-stage: audio slope == response slope (~12 dB/oct)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.type = FilterType::HighPass; p.freq = 1000.0; p.Q = 0.707; p.slope = 24; p.swept = true;
        band.setParams (p);
        auto ap = [&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); };
        const double a125 = sineGainDb (ap, 125.0, fs), a250 = sineGainDb (ap, 250.0, fs);
        const double audioSlope = a250 - a125;
        const double r125 = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 125.0 / fs)));
        const double r250 = 20.0 * std::log10 (std::abs (band.response (2.0 * kPi * 250.0 / fs)));
        const double respSlope = r250 - r125;
        std::printf ("      swept HP24: audioSlope=%.1f  responseSlope=%.1f dB/oct\n", audioSlope, respSlope);
        expectTrue (audioSlope > 9.0 && audioSlope < 15.0, "audio ~12 dB/oct (single SVF, not 24)");
        expectNear (respSlope, audioSlope, 2.5,            "response agrees with audio (no 24-vs-12 mismatch)");
    }

    group ("variable HP slopes: analytic rolloff ~ N dB/oct");
    {
        for (int slope : { 6, 12, 24, 36, 48, 72, 96 })
        {
            BandParams p; p.on = true; p.type = FilterType::HighPass; p.freq = 1000.0; p.slope = slope;
            const double r1 = 20.0 * std::log10 (std::abs (bandResponse (p, fs, 2.0 * kPi * 250.0 / fs)));
            const double r2 = 20.0 * std::log10 (std::abs (bandResponse (p, fs, 2.0 * kPi * 500.0 / fs)));
            const double oct = r2 - r1;   // one octave (250->500) in the stopband
            std::printf ("      HP %2d dB/oct: octave rolloff = %.1f dB\n", slope, oct);
            expectTrue (oct > slope * 0.85 && oct < slope * 1.18, "HP slope ~ " + std::to_string (slope) + " dB/oct");
        }
    }

    group ("bypass: a bypassed band is transparent");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams p; p.on = true; p.bypass = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 12.0;
        band.setParams (p);
        const double g = sineGainDb ([&] (float* const* ch, int nc, int n) { band.processBlock (ch, nc, n); }, 1000.0, fs);
        expectNear (g, 0.0, 0.05, "bypassed band passes audio at 0 dB");
    }

    group ("turning a band off resets state (no stale tail on re-enable)");
    {
        EqBand band; band.prepare (fs, 1);
        BandParams on; on.on = true; on.type = FilterType::Bell; on.freq = 1000.0; on.Q = 2.0; on.gainDb = 6.0;
        band.setParams (on);
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); buf[0] = 1.0f; float* ch[1] = { buf.data() };
        band.processBlock (ch, 1, n);                              // build state from an impulse
        BandParams off = on; off.on = false; band.setParams (off);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f;
        band.processBlock (ch, 1, n);                              // off block -> resets state
        band.setParams (on);                                      // re-enable
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f;
        band.processBlock (ch, 1, n);                              // silent in -> must be silent out
        double tail = 0.0; for (int i = 0; i < n; ++i) tail += std::fabs (buf[(size_t) i]);
        expectTrue (tail == 0.0, "re-enabled band starts from a clean state");
    }

    group ("topology toggle (swept<->static, 24<->12) clears stale state");
    {
        const int n = 64; std::vector<float> buf ((size_t) n, 0.0f); float* ch[1] = { buf.data() };

        EqBand band; band.prepare (fs, 1);
        BandParams sw; sw.on = true; sw.type = FilterType::Bell; sw.freq = 1000.0; sw.Q = 4.0; sw.gainDb = 6.0; sw.swept = true;
        band.setParams (sw);
        buf[0] = 1.0f; band.processBlock (ch, 1, n);                                          // SVF builds state
        BandParams st = sw; st.swept = false; band.setParams (st);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);     // -> static (resets)
        band.setParams (sw);                                                                 // -> swept again
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; band.processBlock (ch, 1, n);
        double t1 = 0.0; for (int i = 0; i < n; ++i) t1 += std::fabs (buf[(size_t) i]);
        expectTrue (t1 == 0.0, "no stale SVF tail after swept<->static");

        EqBand hp; hp.prepare (fs, 1);
        BandParams p24; p24.on = true; p24.type = FilterType::HighPass; p24.freq = 1000.0; p24.Q = 0.707; p24.slope = 24;
        hp.setParams (p24);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; buf[0] = 1.0f; hp.processBlock (ch, 1, n);   // bq0+bq1 state
        BandParams p12 = p24; p12.slope = 12; hp.setParams (p12);
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; hp.processBlock (ch, 1, n);                  // -> single (resets)
        hp.setParams (p24);                                                                             // -> 24 again
        for (int i = 0; i < n; ++i) buf[(size_t) i] = 0.0f; hp.processBlock (ch, 1, n);
        double t2 = 0.0; for (int i = 0; i < n; ++i) t2 += std::fabs (buf[(size_t) i]);
        expectTrue (t2 == 0.0, "no stale bq1 tail after 24<->12 dB/oct");
    }

    group ("EqEngine: silence stays silent, no NaN (stereo)");
    {
        EqEngine eng; eng.prepare (fs, 256, 2);
        BandParams p; p.on = true; p.type = FilterType::HighShelf; p.freq = 8000.0; p.gainDb = 12.0;
        eng.setBand (0, p);
        const int n = 256;
        std::vector<float> L ((size_t) n, 0.0f), R ((size_t) n, 0.0f);
        float* ch[2] = { L.data(), R.data() };
        for (int b = 0; b < 10; ++b) eng.process (ch, 2, n);
        double energy = 0.0;
        for (int i = 0; i < n; ++i) energy += std::fabs (L[(size_t) i]) + std::fabs (R[(size_t) i]);
        expectTrue (energy == 0.0, "silence stays silent");
        expectTrue (! anyNaN (L.data(), n) && ! anyNaN (R.data(), n), "no NaN");
    }

    group ("EqEngine: 6-channel (5.1) static bell — every channel +6 dB, no cross-talk");
    {
        const int C = 6;
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 6.0;
        eng.setBand (0, p);
        const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 1000.0, fs);
        for (int c = 0; c < C; ++c)
            expectNear (g[(size_t) c], 6.0, 0.4, "ch " + std::to_string (c) + " static bell +6 dB");
    }

    group ("EqEngine: 8-channel (7.1) swept SVF bell — per-channel SVF state, every channel +6 dB");
    {
        const int C = 8;                                   // exercises Svf ic1[c]/ic2[c] for c = 0..7
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 2000.0; p.Q = 3.0; p.gainDb = 6.0; p.swept = true;
        eng.setBand (0, p);
        const auto g = multiSineGainDb ([&] (float* const* ch, int nc, int n) { eng.process (ch, nc, n); }, C, 2000.0, fs);
        for (int c = 0; c < C; ++c)
            expectNear (g[(size_t) c], 6.0, 0.5, "ch " + std::to_string (c) + " swept bell +6 dB");
    }

    group ("EqEngine: 16-channel (max cap) silence stays silent, no NaN");
    {
        const int C = teq::kMaxChannels;                   // 16: top of the supported range
        EqEngine eng; eng.prepare (fs, 256, C);
        BandParams p; p.on = true; p.type = FilterType::HighShelf; p.freq = 6000.0; p.gainDb = 9.0;
        eng.setBand (0, p);
        const int n = 256;
        std::vector<std::vector<float>> bufs ((size_t) C, std::vector<float> ((size_t) n, 0.0f));
        std::vector<float*> ptr ((size_t) C);
        for (int c = 0; c < C; ++c) ptr[(size_t) c] = bufs[(size_t) c].data();
        for (int b = 0; b < 10; ++b) eng.process (ptr.data(), C, n);
        double energy = 0.0; bool nan = false;
        for (int c = 0; c < C; ++c) for (int i = 0; i < n; ++i)
        { energy += std::fabs (bufs[(size_t) c][(size_t) i]); nan = nan || ! std::isfinite (bufs[(size_t) c][(size_t) i]); }
        expectTrue (energy == 0.0, "16ch silence stays silent");
        expectTrue (! nan,         "16ch no NaN");
    }

    group ("EqEngine: 16-channel independence — every channel bit-identical to its own mono engine");
    {
        const int C = teq::kMaxChannels;     // 16, each fed a DISTINCT (uncorrelated) frequency
        const int n = 256;
        auto signalAt = [&] (int c, int i) { return (float) (0.2 * std::sin (2.0 * kPi * (110.0 + 37.0 * c) * i / fs)); };

        // Run C independent single-channel engines, then the one 16-channel engine on the same inputs.
        // Per-channel state must be fully separate -> bit-for-bit equal. Any cross-talk / wrong index
        // / channel-count-dependent math (incl. on the top channel 15) makes maxErr != 0.
        auto checkIndependence = [&] (bool swept, const std::string& label)
        {
            BandParams p; p.on = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 3.0; p.gainDb = 6.0; p.swept = swept;

            std::vector<std::vector<float>> ref ((size_t) C, std::vector<float> ((size_t) n));
            for (int c = 0; c < C; ++c)
            {
                EqEngine mono; mono.prepare (fs, n, 1); mono.setBand (0, p);
                for (int i = 0; i < n; ++i) ref[(size_t) c][(size_t) i] = signalAt (c, i);
                float* m[1] = { ref[(size_t) c].data() };
                mono.process (m, 1, n);
            }

            EqEngine eng; eng.prepare (fs, n, C); eng.setBand (0, p);
            std::vector<std::vector<float>> buf ((size_t) C, std::vector<float> ((size_t) n));
            std::vector<float*> ptr ((size_t) C);
            for (int c = 0; c < C; ++c)
            {
                for (int i = 0; i < n; ++i) buf[(size_t) c][(size_t) i] = signalAt (c, i);
                ptr[(size_t) c] = buf[(size_t) c].data();
            }
            eng.process (ptr.data(), C, n);

            double maxErr = 0.0;
            for (int c = 0; c < C; ++c) for (int i = 0; i < n; ++i)
                maxErr = std::max (maxErr, (double) std::fabs (buf[(size_t) c][(size_t) i] - ref[(size_t) c][(size_t) i]));
            expectTrue (maxErr == 0.0, label + ": all 16 channels bit-identical to an independent engine");
        };

        checkIndependence (false, "static biquad");
        checkIndependence (true,  "swept SVF");
    }

    group ("EqBand M/S dual-mode: independent Mid & Side lanes; the idle axis is preserved");
    {
        const int block = 256;
        // steady-state gain (dB) in the Mid (side=false) or Side (side=true) domain at frequency f.
        // Drive a pure-Mid (L=R) or pure-Side (L=-R) sine so the other domain is exactly zero.
        auto msProbe = [&] (EqBand& band, bool side, double f)
        {
            const double amp = 0.25, dp = 2.0 * kPi * f / fs;
            std::vector<float> L ((size_t) block), R ((size_t) block);
            double phase = 0.0;
            auto fill = [&] { for (int n = 0; n < block; ++n) {
                                  const float v = (float) (amp * std::sin (phase));
                                  phase += dp; if (phase > 2.0 * kPi) phase -= 2.0 * kPi;
                                  L[(size_t) n] = v; R[(size_t) n] = side ? -v : v; } };
            for (int b = 0; b < 400; ++b) { fill(); float* ch[2] = { L.data(), R.data() }; band.processBlock (ch, 2, block); }
            fill();
            std::vector<double> inDom ((size_t) block);
            for (int n = 0; n < block; ++n) inDom[(size_t) n] = side ? 0.5 * (L[(size_t) n] - R[(size_t) n]) : 0.5 * (L[(size_t) n] + R[(size_t) n]);
            { float* ch[2] = { L.data(), R.data() }; band.processBlock (ch, 2, block); }
            double ir = 0.0, orr = 0.0;
            for (int n = 0; n < block; ++n) { const double dm = side ? 0.5 * (L[(size_t) n] - R[(size_t) n]) : 0.5 * (L[(size_t) n] + R[(size_t) n]);
                                              ir += inDom[(size_t) n] * inDom[(size_t) n]; orr += dm * dm; }
            return 10.0 * std::log10 (orr / std::max (1e-30, ir));
        };

        { EqBand b; b.prepare (fs, 2);                                  // Mid +12 @1k, Side flat
          BandParams p; p.on = true; p.ms = true;
          p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 12.0;
          p.sOn = true; p.sType = FilterType::Bell; p.sFreq = 1000.0; p.sQ = 2.0; p.sGainDb = 0.0;
          b.setParams (p);
          expectNear (msProbe (b, false, 1000.0), 12.0, 0.5, "Mid lane +12 dB at 1k");
          expectNear (msProbe (b, true,  1000.0),  0.0, 0.1, "Side lane flat at 1k"); }

        { EqBand b; b.prepare (fs, 2);                                  // Side -9 @3k, Mid flat
          BandParams p; p.on = true; p.ms = true;
          p.type = FilterType::Bell; p.freq = 3000.0; p.Q = 2.0; p.gainDb = 0.0;
          p.sOn = true; p.sType = FilterType::Bell; p.sFreq = 3000.0; p.sQ = 2.0; p.sGainDb = -9.0;
          b.setParams (p);
          expectNear (msProbe (b, true,  3000.0), -9.0, 0.5, "Side lane -9 dB at 3k");
          expectNear (msProbe (b, false, 3000.0),  0.0, 0.1, "Mid lane flat at 3k"); }

        { EqBand b; b.prepare (fs, 2);                                  // independent freq: Mid +6 @500, Side -6 @5k
          BandParams p; p.on = true; p.ms = true;
          p.type = FilterType::Bell; p.freq = 500.0;  p.Q = 2.0; p.gainDb =  6.0;
          p.sOn = true; p.sType = FilterType::Bell; p.sFreq = 5000.0; p.sQ = 2.0; p.sGainDb = -6.0;
          b.setParams (p);
          expectNear (msProbe (b, false, 500.0),   6.0, 0.5, "Mid +6 at its own freq");
          expectNear (msProbe (b, true,  5000.0), -6.0, 0.5, "Side -6 at its own freq");
          expectNear (msProbe (b, false, 5000.0),  0.0, 0.5, "Mid flat at the Side freq");
          expectNear (msProbe (b, true,  500.0),   0.0, 0.5, "Side flat at the Mid freq"); }

        { EqBand b; b.prepare (fs, 2);                                  // Side disabled -> Side axis bit-exact
          BandParams p; p.on = true; p.ms = true; p.sOn = false;
          p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 12.0;
          b.setParams (p);
          const int n = 200; std::vector<float> L ((size_t) n), R ((size_t) n), S0 ((size_t) n);
          for (int i = 0; i < n; ++i) { L[(size_t) i] = (float) (0.20 * std::sin (2.0 * kPi * 700.0 * i / fs));
                                        R[(size_t) i] = (float) (0.15 * std::sin (2.0 * kPi * 700.0 * i / fs + 0.7));
                                        S0[(size_t) i] = 0.5f * (L[(size_t) i] - R[(size_t) i]); }
          float* ch[2] = { L.data(), R.data() }; b.processBlock (ch, 2, n);
          double sErr = 0.0; for (int i = 0; i < n; ++i) sErr = std::max (sErr, (double) std::fabs (0.5f * (L[(size_t) i] - R[(size_t) i]) - S0[(size_t) i]));
          expectNear (sErr, 0.0, 1e-5, "Side disabled: Side axis (L-R) preserved to float precision"); }

        { EqBand b; b.prepare (fs, 2);                                  // Mid bypassed + Side active
          BandParams p; p.on = true; p.ms = true; p.bypass = true;
          p.sOn = true; p.sType = FilterType::Bell; p.sFreq = 2000.0; p.sQ = 2.0; p.sGainDb = 12.0;
          b.setParams (p);
          expectNear (msProbe (b, true,  2000.0), 12.0, 0.6, "Side processes while Mid is bypassed");
          expectNear (msProbe (b, false, 2000.0),  0.0, 0.1, "Mid axis preserved while bypassed"); }

        { EqBand b; b.prepare (fs, 1);                                  // mono + M/S -> Mid lane on the mono channel
          BandParams p; p.on = true; p.ms = true; p.type = FilterType::Bell; p.freq = 1000.0; p.Q = 2.0; p.gainDb = 6.0;
          b.setParams (p);
          const double g = sineGainDb ([&] (float* const* ch, int nc, int n) { b.processBlock (ch, nc, n); }, 1000.0, fs);
          expectNear (g, 6.0, 0.5, "mono + M/S: Mid lane applies to the mono channel"); }
    }
}
