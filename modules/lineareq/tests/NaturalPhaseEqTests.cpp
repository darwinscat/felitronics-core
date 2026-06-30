// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for NaturalPhaseEq — the mixed-phase ("Natural") rendering. Decisive properties:
//   (1) a FLAT EQ renders a unit impulse at `bulkDelay` (unity pass-through at the reported latency);
//   (2) the realised magnitude matches the bank's target (mixed phase preserves |H|);
//   (3) latency = (1−k)·L/2 — strictly LESS than linear's L/2;
//   (4) the FIR is NON-symmetric and front-loaded (mixed phase, not linear), with less pre- than post-energy;
//   (5) process() never allocates (RT-safe) and runs mono + stereo.

#include <felitronics_test.h>
#include <felitronics/lineareq/NaturalPhaseEq.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using NPE = lineareq::NaturalPhaseEq;

static double firMagDb (const std::vector<float>& fir, double f, double sr)
{
    const double w = 2.0 * core::kPi * f / sr;
    double re = 0.0, im = 0.0;
    for (int k = 0; k < (int) fir.size(); ++k) { re += fir[(std::size_t) k] * std::cos (w * k); im -= fir[(std::size_t) k] * std::sin (w * k); }
    return 20.0 * std::log10 (std::max (1e-12, std::sqrt (re * re + im * im)));
}

int main()
{
    std::printf ("felitronics::lineareq NaturalPhaseEq tests\n");
    const double sr = 48000.0;
    const int Q = 1;                                                  // quality 1 → L = 4096

    NPE np; np.prepare (sr, 512, 2, Q, 0.5f);
    const int L = np.firSize();
    const int delay = np.latencySamples();                           // = (1−k)·L/2 = L/4 for k=0.5
    std::vector<float> firMid ((std::size_t) L), firSide ((std::size_t) L);
    eq::BandParams flat[1];

    // --- (1) flat EQ → a unit impulse at the bulk delay (unity pass-through at the latency) ---
    test::group ("NaturalPhaseEq flat EQ == unit impulse at the bulk delay");
    {
        np.buildFir (flat, 0, false, firMid.data());
        double maxOther = 0.0;
        for (int i = 0; i < L; ++i) if (i != delay) maxOther = std::max (maxOther, (double) std::fabs (firMid[(std::size_t) i]));
        test::approx (firMid[(std::size_t) delay], 1.0, 1e-3, "peak tap == 1.0 at bulkDelay (unity gain)");
        test::ok (maxOther < 1e-3, "every other tap ≈ 0 (a clean delayed delta)");
    }

    // --- (2) realised magnitude matches the bank's target (mixed phase preserves |H|) ---
    test::group ("NaturalPhaseEq magnitude matches the EQ target");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0;
        np.buildFir (b, 1, false, firMid.data());
        test::approx (firMagDb (firMid, 1000.0, sr),  6.0, 0.8, "+6 dB bell → ~+6 dB at 1 kHz");
        test::approx (firMagDb (firMid,  100.0, sr),  0.0, 0.6, "≈ 0 dB well below the bell");
        eq::BandParams c[1]; c[0].on = true; c[0].type = eq::FilterType::Bell; c[0].freq = 3000.0; c[0].Q = 3.0; c[0].gainDb = -8.0;
        np.buildFir (c, 1, false, firMid.data());
        test::approx (firMagDb (firMid, 3000.0, sr), -8.0, 0.8, "−8 dB cut → ~−8 dB at 3 kHz");
    }

    // --- (3) latency is (1−k)·L/2 and strictly lighter than linear's L/2 ---
    test::group ("NaturalPhaseEq latency lighter than linear");
    {
        test::ok (delay == L / 4, "k=0.5 → bulk delay == L/4");
        test::ok (delay < L / 2,  "lighter than a linear-phase FIR of the same length (L/2)");
    }

    // --- (4) the FIR is mixed phase: NON-symmetric + front-loaded (less pre- than post-energy) ---
    test::group ("NaturalPhaseEq FIR is mixed phase (not linear/symmetric)");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1200.0; b[0].Q = 2.0; b[0].gainDb = 8.0;
        np.buildFir (b, 1, false, firMid.data());
        double asym = 0.0; for (int i = 0; i < L; ++i) asym = std::max (asym, (double) std::fabs (firMid[(std::size_t) i] - firMid[(std::size_t) (L - 1 - i)]));
        test::ok (asym > 1e-3, "NOT symmetric (mixed phase, unlike linear)");
        double pre = 0.0, post = 0.0;
        for (int i = 0; i < delay; ++i)     pre  += std::fabs (firMid[(std::size_t) i]);
        for (int i = delay + 1; i < L; ++i) post += std::fabs (firMid[(std::size_t) i]);
        test::ok (pre < post, "less pre-ring than post-ring (energy front-loaded toward minimum phase)");
    }

    // --- (5) process() is RT-safe (no alloc) + runs mono and stereo ---
    test::group ("NaturalPhaseEq process no-alloc + mono/stereo");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 1.5; b[0].gainDb = 4.0;
        np.setBands (b, 1);                                          // build (message thread, allocates) BEFORE the snapshot
        std::vector<float> lch (512, 0.1f), rch (512, -0.1f);
        float* io[2] { lch.data(), rch.data() };
        np.process (io, 2, 512);                                     // consume the initial fade-in
        const long before = g_allocs.load();
        np.process (io, 2, 512);
        float* mono[1] { lch.data() };
        np.process (mono, 1, 512);                                   // mono path
        const long after = g_allocs.load();
        test::ok (after == before, "process() performed zero heap allocations (stereo + mono)");
        bool finite = true; for (float v : lch) finite = finite && std::isfinite (v);
        test::ok (finite, "output finite");
    }

    // --- (6) STEEP filter accuracy: the kept L taps must still track a steep HP's magnitude (council:
    //     codex — bulkDelay/L truncation is approximate; this bounds the error on a hard curve) ---
    test::group ("NaturalPhaseEq tracks a steep high-pass");
    {
        NPE hp; hp.prepare (sr, 512, 2, 2, 0.5f);                    // quality 2 → L = 8192 (headroom for ringing)
        std::vector<float> fir ((std::size_t) hp.firSize());
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::HighPass; b[0].freq = 1000.0; b[0].Q = 0.707; b[0].slope = 24;
        hp.buildFir (b, 1, false, fir.data());
        test::approx (firMagDb (fir, 4000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 4000.0, sr), 1.0, "pass-band (4 kHz) ≈ target");
        test::approx (firMagDb (fir, 2000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 2000.0, sr), 1.2, "near cutoff (2 kHz) ≈ target");
        test::ok (firMagDb (fir, 250.0, sr) < -14.0, "stop-band (250 Hz) strongly attenuated");
    }

    // --- (7) a deep narrow cut keeps its depth (steep/narrow feature survives the truncation) ---
    test::group ("NaturalPhaseEq keeps a deep narrow cut");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 2000.0; b[0].Q = 6.0; b[0].gainDb = -18.0;
        np.buildFir (b, 1, false, firMid.data());                   // np: L = 4096
        test::approx (firMagDb (firMid, 2000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 2000.0, sr), 2.0, "−18 dB narrow cut depth ≈ target");
        test::approx (firMagDb (firMid,  500.0, sr), 0.0, 0.6, "≈ 0 dB away from the cut");
    }

    return test::report();
}
