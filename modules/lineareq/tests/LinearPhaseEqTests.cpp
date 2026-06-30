// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the linear-phase EQ. The decisive properties (DSP council):
//   (1) a FLAT EQ renders an EXACT unit impulse at the centre tap N/2 — proves the build needs no run-time
//       gain hack (the 1/N inverse + the centre-1.0 Blackman-Harris suffice); (2) the FIR is symmetric →
//       exactly linear phase; (3) the realised magnitude matches the Eq bank's target (within window ripple);
//       (4) M/S builds independent Mid & Side IRs; (5) end-to-end the convolution delays by N/2 at unit gain
//       and applies the bell's gain to a tone; (6) process() never allocates; latency = N/2.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/lineareq/LinearPhaseEq.h>

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
using LPE = lineareq::LinearPhaseEq;

// |H(f)| in dB of an (N+1)-tap FIR via a direct DFT at one frequency.
static double firMagDb (const std::vector<float>& fir, double f, double sr)
{
    const double w = 2.0 * core::kPi * f / sr;
    double re = 0.0, im = 0.0;
    for (int k = 0; k < (int) fir.size(); ++k) { re += fir[(std::size_t) k] * std::cos (w * k); im -= fir[(std::size_t) k] * std::sin (w * k); }
    return 20.0 * std::log10 (std::max (1e-12, std::sqrt (re * re + im * im)));
}

int main()
{
    std::printf ("felitronics::lineareq tests\n");
    const double sr = 48000.0;
    const int Q = 0;                                                 // quality 0 → N = 4096 (fast; higher = same code)

    LPE lp; lp.prepare (sr, 512, 2, Q);
    const int N = lp.firSize();
    std::vector<float> firMid ((std::size_t) (N + 1)), firSide ((std::size_t) (N + 1));
    eq::BandParams flat[1];                                           // numBands=0 → composite magnitude ≡ 1

    // --- (1) THE key property: a flat EQ → an exact unit impulse at N/2 (no run-time gain hack) ---
    test::group ("LinearPhaseEq flat EQ == unit impulse at N/2");
    {
        lp.buildFir (flat, 0, false, firMid.data());
        double maxOther = 0.0; for (int i = 0; i <= N; ++i) if (i != N / 2) maxOther = std::max (maxOther, (double) std::fabs (firMid[(std::size_t) i]));
        test::approx (firMid[(std::size_t) (N / 2)], 1.0, 1e-3, "centre tap == 1.0 (unit gain, no firScale needed)");
        test::ok (maxOther < 1e-3, "every other tap ≈ 0 (a clean delta)");
    }

    // --- (2) the FIR is symmetric → EXACTLY linear phase ---
    test::group ("LinearPhaseEq FIR is symmetric (linear phase)");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1200.0; b[0].Q = 2.0; b[0].gainDb = 8.0;
        lp.buildFir (b, 1, false, firMid.data());
        double asym = 0; for (int i = 0; i <= N; ++i) asym = std::max (asym, (double) std::fabs (firMid[(std::size_t) i] - firMid[(std::size_t) (N - i)]));
        test::ok (asym < 1e-6, "fir[N/2+d] == fir[N/2−d] → symmetric → linear phase");
    }

    // --- (3) the realised magnitude matches the bank's target (within window ripple) ---
    test::group ("LinearPhaseEq magnitude matches the EQ target");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0;
        lp.buildFir (b, 1, false, firMid.data());
        test::approx (firMagDb (firMid, 1000.0, sr), 6.0, 0.7, "+6 dB bell → ~+6 dB at 1 kHz");
        test::approx (firMagDb (firMid,  100.0, sr), 0.0, 0.3, "flat away from the bell (~0 dB at 100 Hz)");
        const double target = eq::EqEngine::magnitudeDbFor (b, 1, 1000.0, sr);   // == the bank's own curve
        test::approx (firMagDb (firMid, 1000.0, sr), target, 0.7, "FIR magnitude tracks EqEngine::magnitudeDbFor");
    }

    // --- (4) M/S: independent Mid & Side IRs from one snapshot ---
    test::group ("LinearPhaseEq M/S builds independent Mid+Side IRs");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].ms = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0;
        b[0].sOn = true; b[0].sType = eq::FilterType::Bell; b[0].sFreq = 1000.0; b[0].sQ = 2.0; b[0].sGainDb = -6.0;
        lp.buildFir (b, 1, false, firMid.data());
        lp.buildFir (b, 1, true,  firSide.data());
        test::approx (firMagDb (firMid,  1000.0, sr),  6.0, 0.7, "Mid IR = +6 dB at 1 kHz");
        test::approx (firMagDb (firSide, 1000.0, sr), -6.0, 0.7, "Side IR = −6 dB at 1 kHz");
    }

    // --- (5a) end-to-end: a flat EQ delays an impulse by N/2 at unit gain ---
    test::group ("LinearPhaseEq end-to-end flat → delayed impulse");
    {
        LPE e; e.prepare (sr, 512, 2, Q); e.setBands (flat, 0);       // unit-impulse IR (crossfades in over ~40 ms)
        const int M = 12000, imp = 6000;                              // impulse well past the crossfade
        std::vector<float> L (M, 0.0f), R (M, 0.0f); L[(std::size_t) imp] = 1.0f; R[(std::size_t) imp] = 1.0f;
        for (int o = 0; o < M; o += 512) { float* io[2] { L.data() + o, R.data() + o }; e.process (io, 2, std::min (512, M - o)); }
        int peakIdx = 0; double peak = 0; for (int i = 0; i < M; ++i) if (std::fabs (L[(std::size_t) i]) > peak) { peak = std::fabs (L[(std::size_t) i]); peakIdx = i; }
        test::ok (peakIdx == imp + N / 2, "output impulse lands at input + N/2 (linear-phase latency)");
        test::approx (peak, 1.0, 0.02, "unit gain (flat EQ passes through, just delayed)");
    }

    // --- (5b) end-to-end: a +6 dB bell lifts a 1 kHz tone ~6 dB ---
    test::group ("LinearPhaseEq end-to-end tone gain");
    {
        LPE e; e.prepare (sr, 512, 2, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0;
        e.setBands (b, 1);
        const int M = 16000;
        std::vector<float> L (M), R (M); for (int i = 0; i < M; ++i) { L[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); R[(std::size_t) i] = L[(std::size_t) i]; }
        for (int o = 0; o < M; o += 512) { float* io[2] { L.data() + o, R.data() + o }; e.process (io, 2, std::min (512, M - o)); }
        double inSq = 0, outSq = 0; int from = M - 4000;              // tail: past crossfade + latency
        for (int i = from; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) L[(std::size_t) i] * L[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq / inSq), 6.0, 0.8, "1 kHz tone out ≈ +6 dB through the +6 dB bell");
    }

    // --- (6) process() never allocates (setBands is the message-thread build; process is RT-safe) ---
    test::group ("LinearPhaseEq process no-alloc");
    {
        LPE e; e.prepare (sr, 512, 2, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 2000.0; b[0].gainDb = -4.0; e.setBands (b, 1);
        std::vector<float> L (512, 0.3f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        e.process (io, 2, 512);
        const long before = g_allocs.load();
        e.process (io, 2, 512); e.process (io, 2, 512);
        test::ok (g_allocs.load() == before, "process() did not allocate");
        test::ok (e.latencySamples() == N / 2, "latencySamples() == N/2");
    }

    // --- (7) MONO path: 1-channel session → the Mid IR only, no crash, correct gain (refutes a review worry) ---
    test::group ("LinearPhaseEq mono path");
    {
        LPE e; e.prepare (sr, 512, 1, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0; e.setBands (b, 1);
        const int M = 16000; std::vector<float> x (M); for (int i = 0; i < M; ++i) x[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int o = 0; o < M; o += 512) { float* io[1] { x.data() + o }; e.process (io, 1, std::min (512, M - o)); }
        double inSq = 0, outSq = 0; for (int i = M - 4000; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) x[(std::size_t) i] * x[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq / inSq), 6.0, 0.8, "mono 1 kHz tone +6 dB through the Mid IR (no crash)");
    }

    // --- (8) an IR swap mid-stream is click-free AND response-correct after it settles (crossfade=N) ---
    test::group ("LinearPhaseEq IR swap settles clean");
    {
        LPE e; e.prepare (sr, 512, 2, Q); e.setBands (flat, 0);       // start flat (unit IR)
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].freq = 1000.0; b[0].Q = 2.0; b[0].gainDb = 6.0;
        const int M = 24000; std::vector<float> L (M), R (M); for (int i = 0; i < M; ++i) { L[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); R[(std::size_t) i] = L[(std::size_t) i]; }
        bool swapped = false;
        for (int o = 0; o < M; o += 512)
        {
            if (! swapped && o >= 8000 && e.setBands (b, 1)) swapped = true;   // swap flat→bell mid-stream
            float* io[2] { L.data() + o, R.data() + o }; e.process (io, 2, std::min (512, M - o));
        }
        double mx = 0; for (int i = 0; i < M; ++i) mx = std::max (mx, (double) std::fabs (L[(std::size_t) i]));
        double inSq = 0, outSq = 0; for (int i = M - 3000; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) L[(std::size_t) i] * L[(std::size_t) i]; }
        test::ok (mx < 1.1, "no full-scale spike across the swap (click-free)");
        test::approx (10.0 * std::log10 (outSq / inSq), 6.0, 0.8, "settled output = +6 dB (the new IR is response-correct, not under-primed)");
    }

    // --- (9) quality table + the unit-impulse property holds at a larger N too ---
    test::group ("LinearPhaseEq quality sizes");
    {
        const int expect[5] { 4096, 8192, 16384, 32768, 131072 };
        bool tab = true; for (int q = 0; q < 5; ++q) tab &= (LPE::firSizeForQuality (q) == expect[q]);
        test::ok (tab && LPE::firSizeForQuality (-1) == 4096 && LPE::firSizeForQuality (99) == 131072, "quality 0..4 → {4096,8192,16384,32768,131072}, clamped");
        LPE e; e.prepare (sr, 512, 2, 2);                            // q=2 → N=16384
        test::ok (e.firSize() == 16384 && e.latencySamples() == 8192, "q=2 → N=16384, latency 8192");
        std::vector<float> fir (16385); e.buildFir (flat, 0, false, fir.data());
        double mo = 0; for (int i = 0; i <= 16384; ++i) if (i != 8192) mo = std::max (mo, (double) std::fabs (fir[(std::size_t) i]));
        test::ok (std::fabs (fir[8192] - 1.0f) < 1e-3 && mo < 1e-3, "flat → unit impulse at N/2 holds at N=16384 too");
    }

    // --- (10) a steep linear-phase high-pass: flat pass-band, strongly (finitely) attenuated stop-band ---
    test::group ("LinearPhaseEq steep high-pass");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::HighPass; b[0].freq = 500.0; b[0].Q = 0.707; b[0].slope = 24;
        lp.buildFir (b, 1, false, firMid.data());
        test::ok (firMagDb (firMid, 2000.0, sr) > -1.0,  "pass-band (2 kHz, 2 oct up) ≈ 0 dB");
        test::ok (firMagDb (firMid,  125.0, sr) < -18.0, "stop-band (125 Hz, 2 oct down) strongly attenuated (finite, per the window)");
    }

    return test::report();
}
