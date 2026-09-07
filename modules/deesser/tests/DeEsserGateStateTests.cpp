// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The eighth place in this repo where a channel that leaves and returns replays frozen recursions.
//
// SplitBand mode runs TWO per-channel recursions over `c < nc` with no ledger: the sidechain band-pass
// and the Linkwitz-Riley crossover. Measured before this suite, 1 s of stereo noise -> 2 s of mono
// silence -> stereo DIGITAL SILENCE: the returning channel emitted 6.25e-02 (-24.1 dBFS) at sample 0
// against 0.000 for the channel that stayed. DynamicEq mode was already safe, because it delegates to
// DynamicEqBand, which carries its own ledger.
//
// The existing DeEsser suite is green and always will be: every one of its cases processes a single
// channel, so a per-channel ledger is something it cannot have an opinion about.

#include <felitronics_test.h>
#include <felitronics/deesser/DeEsser.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <vector>

static std::atomic<int> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

static constexpr double kFs = 48000.0;

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static bool bitEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! (a[i] == b[i])) return false;
    return true;
}

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0005f; }
}

static deesser::DeEsserParams params (deesser::DeEsserMode mode)
{
    deesser::DeEsserParams p;
    p.mode = mode; p.fc = 6500.0; p.thresholdDb = -40.0; p.rangeDb = 10.0;
    return p;
}

int main()
{
    std::printf ("felitronics::deesser — execution-gate state\n");

    const int N = 128;

    group ("channel gate: a returning channel emits nothing from digital silence");
    for (int m = 0; m < 2; ++m)
    {
        const auto mode = (m == 0) ? deesser::DeEsserMode::SplitBand : deesser::DeEsserMode::DynamicEq;
        const char* name = (m == 0) ? "SplitBand" : "DynamicEq";

        deesser::DeEsser d; felitronics::test::run (d.prepare (kFs, N, 2));
        d.setParams (params (mode));

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* io[2] { L.data(), R.data() };

        double charged = 0.0;
        for (int k = 0; k < (int) (kFs * 1.0 / N); ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f);
            fillNoise (R, 7u + (unsigned) k);                 // only the right channel carries programme
            felitronics::test::run (d.process (io, 2, N));
            charged = std::fmax (charged, peakOf (R));
        }
        ok (charged > 0.05, std::string ("precondition ") + name + ": the right channel really was driven");

        for (int k = 0; k < (int) (kFs * 2.0 / N); ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (d.process (io, 1, N)); }

        double worst = 0.0;
        for (int k = 0; k < 20; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (d.process (io, 2, N));
            worst = std::fmax (worst, std::fmax (peakOf (L), peakOf (R)));
        }
        ok (worst == 0.0, std::string (name) + ": silence in, exact zero out (SplitBand was 6.25e-02)");
    }

    group ("isolation: the channel that never left keeps its history bit-exact");
    {
        // The detector is LINKED, so the shared half legitimately follows whichever channels are present.
        // Keep the leaving channel silent: the shared path then sees the same programme either way, and the
        // only thing that can differ is the per-channel state this fix touches.
        deesser::DeEsser dut, ref;
        felitronics::test::run (dut.prepare (kFs, N, 2)); felitronics::test::run (ref.prepare (kFs, N, 2));
        dut.setParams (params (deesser::DeEsserMode::SplitBand));
        ref.setParams (params (deesser::DeEsserMode::SplitBand));

        std::vector<float> d0 ((std::size_t) N), d1 ((std::size_t) N), r0 ((std::size_t) N), r1 ((std::size_t) N);
        float* dio[2] { d0.data(), d1.data() };
        float* rio[2] { r0.data(), r1.data() };

        bool equal = true; double energy = 0.0;
        for (int k = 0; k < 120; ++k)
        {
            fillNoise (d0, 91u + (unsigned) k);
            std::fill (d1.begin(), d1.end(), 0.0f);
            r0 = d0; r1 = d1;
            felitronics::test::run (dut.process (dio, (k >= 40 && k < 80) ? 1 : 2, N));
            felitronics::test::run (ref.process (rio, 2, N));
            equal = equal && bitEqual (d0, r0);
            energy += peakOf (d0);
        }
        ok (energy > 1.0, "precondition: the surviving channel carried signal");
        ok (equal, "channel 0 is bit-equal to the run where nothing ever left");
    }

    group ("RT-safety");
    {
        deesser::DeEsser d; felitronics::test::run (d.prepare (kFs, N, 4));
        d.setParams (params (deesser::DeEsserMode::SplitBand));
        std::vector<float> v[4];
        float* io[4] {};
        for (int c = 0; c < 4; ++c) { v[c].assign ((std::size_t) N, 0.05f); io[c] = v[c].data(); }
        felitronics::test::run (d.process (io, 4, N));
        const int before = g_allocs.load();
        for (int k = 0; k < 40; ++k) felitronics::test::run (d.process (io, (k % 3) + 2, N));
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across 40 blocks of changing width");
    }

    return felitronics::test::report();
}
