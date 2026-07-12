// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Adversarial falsification pass for PowerAmpStage. The external case list is adapted to the real
// header contract: the DC blocker is 10 Hz, voicing is caller-supplied data, and the engine contract
// processes at most two channels.

#include <felitronics_test.h>
#include <felitronics/poweramp/PowerAmpStage.h>

#include <algorithm>
#include <cmath>
#include <vector>

using felitronics::poweramp::Params;
using felitronics::poweramp::PowerAmpStage;
using felitronics::poweramp::Voicing;
using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kSr = 48000.0;

std::vector<float> sine (int warm, int n, double hz, double amp)
{
    std::vector<float> x ((std::size_t) (warm + n + 128), 0.0f);
    for (int i = 0; i < (int) x.size(); ++i)
        x[(std::size_t) i] = (float) (amp * std::sin (2.0 * kPi * hz * (double) i / kSr));
    return x;
}

double magAt (const std::vector<float>& x, int start, int n, double hz)
{
    double re = 0.0, im = 0.0;
    const double w = 2.0 * kPi * hz / kSr;
    for (int i = 0; i < n; ++i)
    {
        const double s = x[(std::size_t) (start + i)];
        re += s * std::cos (w * i);
        im -= s * std::sin (w * i);
    }
    return std::sqrt (re * re + im * im) / (double) n;
}

double rms (const std::vector<float>& x, int start, int n)
{
    double e = 0.0;
    for (int i = 0; i < n; ++i) e += (double) x[(std::size_t) (start + i)] * x[(std::size_t) (start + i)];
    return std::sqrt (e / (double) n);
}

std::vector<float> runMono (Params p, Voicing v, std::vector<float> x, int maxBlock = 512)
{
    PowerAmpStage d;
    d.prepare (kSr, maxBlock, 4);
    d.setParams (p, v);
    float* io[1] { x.data() };
    d.process (io, 1, (int) x.size());
    return x;
}

double toneMag (double hz, Params p, Voicing v)
{
    const int warm = 48000, n = 48000;
    auto y = runMono (p, v, sine (warm, n, hz, 1.0e-4));
    return magAt (y, warm + 31, n, hz);
}

bool allFiniteAndBounded (const std::vector<float>& x)
{
    for (float v : x)
        if (! std::isfinite (v) || std::fabs (v) > 1.0e6f) return false;
    return true;
}
} // namespace

int main()
{
    std::printf ("felitronics::poweramp PowerAmpStage FALSIFY tests\n");

    group ("DeepSeek 4.1/4.2 adapted: real DC blocker corner is 10 Hz, not 7 Hz");
    {
        Params p; p.driveDb = 0.0f; p.autoComp = 1.0f;
        Voicing v;
        const double a10 = toneMag (10.0, p, v);
        const double a1k = toneMag (1000.0, p, v);
        const double cornerDb = 20.0 * std::log10 (a10 / a1k);
        std::printf ("      10 Hz / 1 kHz magnitude = %.3f dB\n", cornerDb);
        approx (cornerDb, -3.0, 0.75, "10 Hz sine is near the high-pass -3 dB corner");

        auto step = runMono (p, v, std::vector<float> (48000, 0.5f));
        ok (allFiniteAndBounded (step), "DC step output remains finite and clamp-bounded");
        ok (std::fabs (step.back()) < 1.0e-3f, "sustained DC step settles back near zero");
    }

    group ("DeepSeek 4.3 adapted: huge finite input is sanitized and rail-clamped finite");
    {
        Params p; p.driveDb = 36.0f; p.sag = 1.0f; p.autoComp = 1.0f;
        Voicing v; v.sagMaxDroop = 0.9f; v.sagFastMs = 1.0f; v.sagRecoveryMs = 20.0f;
        auto y = runMono (p, v, std::vector<float> (4096, 1.0e30f));
        ok (allFiniteAndBounded (y), "sustained huge finite input produces only finite samples inside +/-1e6");
    }

    group ("DeepSeek 4.4 materialized: autoComp=1 keeps 0 dB and +20 dB tiny-signal RMS close");
    {
        Voicing v;
        Params p0; p0.driveDb = 0.0f;  p0.autoComp = 1.0f;
        Params p1; p1.driveDb = 20.0f; p1.autoComp = 1.0f;
        auto x = sine (24000, 24000, 1000.0, 1.0e-4);
        const auto y0 = runMono (p0, v, x);
        const auto y1 = runMono (p1, v, x);
        const double deltaDb = 20.0 * std::log10 (rms (y1, 24000 + 31, 20000) / rms (y0, 24000 + 31, 20000));
        approx (deltaDb, 0.0, 1.5, "autoComp=1 keeps +20 dB drive within 1.5 dB of 0 dB drive");
    }

    group ("DeepSeek 4.5 adapted: presence max boosts 10 kHz / 100 Hz ratio by more than 3 dB");
    {
        Voicing v; v.presenceHz = 5000.0f; v.presenceMaxDb = 6.0f;
        Params off; off.presence = 0.0f; off.autoComp = 1.0f;
        Params on;  on.presence  = 1.0f; on.autoComp  = 1.0f;
        const double ratioOff = toneMag (10000.0, off, v) / toneMag (100.0, off, v);
        const double ratioOn  = toneMag (10000.0, on,  v) / toneMag (100.0, on,  v);
        const double boostDb = 20.0 * std::log10 (ratioOn / ratioOff);
        ok (boostDb > 3.0, "presence=1 boosts 10 kHz/100 Hz ratio by > 3 dB (got " + std::to_string (boostDb) + " dB)");
    }

    group ("own attack: process(n > maxBlock) is bit-identical to explicit maxBlock chunking");
    {
        Params p; p.driveDb = 12.0f; p.autoComp = 1.0f;
        Voicing v;
        std::vector<float> one (1500), split (1500);
        for (int i = 0; i < 1500; ++i)
            one[(std::size_t) i] = split[(std::size_t) i] = 0.2f * (float) std::sin (0.013 * i) + 0.01f * (float) (i & 7);

        PowerAmpStage a, b;
        a.prepare (kSr, 128, 4); b.prepare (kSr, 128, 4);
        a.setParams (p, v); b.setParams (p, v);
        { float* io[1] { one.data() }; a.process (io, 1, (int) one.size()); }
        for (int off = 0; off < (int) split.size(); off += 128)
        {
            b.setParams (p, v);
            float* io[1] { split.data() + off };
            b.process (io, 1, std::min (128, (int) split.size() - off));
        }
        ok (one == split, "oversized call matches explicit 128-sample chunking bit-for-bit");
    }

    group ("own attack: channel 2 is outside the stereo engine contract and remains untouched");
    {
        Params p; p.driveDb = 6.0f; p.autoComp = 1.0f;
        Voicing v;
        std::vector<float> a (512, 0.1f), b (512, -0.1f), sentinel (512, 123.0f);
        PowerAmpStage d; d.prepare (kSr, 512, 4); d.setParams (p, v);
        float* io[3] { a.data(), b.data(), sentinel.data() };
        d.process (io, 3, 512);
        bool untouched = true;
        for (float x : sentinel) untouched = untouched && (x == 123.0f);
        ok (untouched, "third channel is ignored without touching its buffer");
        ok (allFiniteAndBounded (a) && allFiniteAndBounded (b), "processed stereo channels remain finite");
    }

    return felitronics::test::report();
}
