// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// dehumbench — CPU of dehum::HumRemover<64> per 64-frame block, stereo, 24 harmonics, by path: playing (hum +
// plucks, the playing width), quiet (hum only: learn width, FLL, detector), learning (a LEARN running),
// listening (Auto with no hum: the detector alone, output untouched), digital silence, and switching (the
// mode flips 50 <-> 60 every 0.25 s: ramps and re-arms). Target: <= 0.5 %
// of a core at 48 kHz. Dev tool (top-level builds only): build Release, run manually; wall-clock numbers,
// best of five passes over 20 s of audio — read them as orders of magnitude between runs.
#include <felitronics/dehum/HumRemover.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace dh = felitronics::dehum;

enum class Path { Playing, Quiet, Learning, Listening, Silence, Switching };

static void fill (std::vector<float>& l, std::vector<float>& r, double fs, Path path)
{
    const double pi = 3.14159265358979323846;
    std::uint32_t s = 12345u;
    for (std::size_t n = 0; n < l.size(); ++n)
    {
        const double t = (double) n / fs;
        double hum = 0.0;
        if (path != Path::Listening && path != Path::Silence)
            for (int k = 1; k <= 24; ++k) hum += 0.0316 / k * std::sin (2.0 * pi * 50.0 * k * t + k);
        double note = 0.0;
        if (path == Path::Playing || path == Path::Listening)
        {
            const double tp = std::fmod (t, 0.5);
            note = 0.1 * std::exp (-tp) * std::sin (2.0 * pi * 329.63 * t);
        }
        s = s * 1664525u + 1013904223u;
        const double noise = path == Path::Silence ? 0.0 : 1e-4 * ((double) s / 4294967296.0 - 0.5);
        l[n] = (float) (hum + note + noise);
        r[n] = (float) (0.7 * hum + note + noise);
    }
}

static void bench (double fs, Path path, const char* name)
{
    const int N = (int) (20.0 * fs), B = 64;
    std::vector<float> l0 ((std::size_t) N), r0 ((std::size_t) N);
    fill (l0, r0, fs, path);
    double best = 1e9;
    for (int pass = 0; pass < 5; ++pass)
    {
        std::vector<float> l = l0, r = r0;
        dh::HumRemover<64> h;
        dh::HumRemoverParams p; p.mains = path == Path::Listening ? dh::Mains::Auto : dh::Mains::Hz50;
        h.setParams (p);
        if (! h.prepare (fs, B, 2)) { std::printf ("prepare failed\n"); return; }
        if (path == Path::Learning) h.beginLearn();
        const auto t0 = std::chrono::steady_clock::now();
        for (int o = 0; o < N; o += B)
        {
            if (path == Path::Switching && o % 12032 == 0) { p.mains = (o / 12032) % 2 ? dh::Mains::Hz60 : dh::Mains::Hz50; h.setParams (p); }
            float* io[2] { l.data() + o, r.data() + o };
            if (! h.process (io, 2, std::min (B, N - o))) { std::printf ("process refused\n"); return; }
        }
        const double sec = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
        best = std::min (best, sec);
        if (pass == 0 && l[(std::size_t) N - 1] > 1e30f) std::printf ("?");   // keep the output observable
    }
    const double blocks = (double) N / B;
    std::printf ("  %6.0f Hz  %-10s %7.3f us per 64-frame block   %6.3f %% of a core\n",
                 fs, name, 1e6 * best / blocks, 100.0 * best / (N / fs));
}

int main()
{
    std::printf ("dehumbench — HumRemover<64>, stereo, 24 harmonics (best of 5 x 20 s)\n");
    for (double fs : { 48000.0, 96000.0, 192000.0 })
    {
        bench (fs, Path::Playing,   "playing");
        bench (fs, Path::Quiet,     "quiet");
        bench (fs, Path::Learning,  "learning");
        bench (fs, Path::Listening, "listening");
        bench (fs, Path::Silence,   "silence");
        bench (fs, Path::Switching, "switching");
    }
    return 0;
}
