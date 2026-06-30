// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the True-Peak meter. The decisive one is the canonical BS.1770 inter-sample case:
// a full-scale fs/4 sine sampled so the grid lands at ±0.707 (sample peak −3.01 dBFS) while the real peak is
// 0 dBFS — a true-peak meter must recover ~0 dBTP. Plus: TP ≥ sample peak always; a constant reads exactly
// (unity-DC phases); flat pass-band; factor by sample rate; stereo→max; ballistics; no-alloc; finite guard.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/analysis/TruePeakMeter.h>

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
using TPM = analysis::TruePeakMeter;

// run a mono tone through a fresh meter, return its true-peak dBTP
static double tonePeakDb (double f, double amp, double phase, int N, double sr)
{
    TPM m; m.prepare (sr, 1024, 1);
    std::vector<float> x (N); for (int i = 0; i < N; ++i) x[i] = (float) (amp * std::sin (2.0 * core::kPi * f * i / sr + phase));
    for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
    return m.truePeakDb();
}

// same, but raised-cosine fade-in so there is NO onset step (→ measures the STEADY tone, not a step overshoot)
static double tonePeakDbFaded (double f, double amp, double phase, int N, double sr)
{
    TPM m; m.prepare (sr, 1024, 1);
    const int fade = 2000;
    std::vector<float> x (N);
    for (int i = 0; i < N; ++i)
    {
        const double env = i < fade ? 0.5 * (1.0 - std::cos (core::kPi * i / fade)) : 1.0;
        x[i] = (float) (env * amp * std::sin (2.0 * core::kPi * f * i / sr + phase));
    }
    for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
    return m.truePeakDb();
}

int main()
{
    std::printf ("felitronics::analysis TruePeakMeter tests\n");
    const double sr = 48000.0;

    // --- THE canonical test: a 0 dBFS fs/4 sine whose samples straddle the peak reads ~0 dBTP, not −3 ---
    test::group ("TruePeak recovers the fs/4 inter-sample peak");
    {
        TPM m; m.prepare (sr, 1024, 1);
        const int N = 8000; std::vector<float> x (N);
        for (int i = 0; i < N; ++i) x[i] = (float) std::sin (2.0 * core::kPi * (sr / 4.0) * i / sr + core::kPi / 4.0);   // samples = ±0.707
        for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
        test::approx (m.samplePeakDb(), -3.0103, 0.05, "raw sample peak is −3.01 dBFS (grid misses the crest)");
        test::ok (m.truePeakDb() > -0.5 && m.truePeakDb() < 0.3, "true peak recovers ~0 dBTP (the real crest between samples)");
        test::ok (m.truePeakDb() - m.samplePeakDb() > 2.5, "true peak exceeds sample peak by ~3 dB (inter-sample recovery)");
    }

    // --- TP is never LESS than the sample peak (the interpolator passes through the grid) ---
    test::group ("TruePeak >= sample peak");
    {
        TPM m; m.prepare (sr, 1024, 1);
        unsigned long s = 5; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 6000; std::vector<float> x (N); for (int i = 0; i < N; ++i) x[i] = 0.6f * rng();
        for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
        test::ok (m.truePeakDb() >= m.samplePeakDb() - 1e-4, "true peak ≥ sample peak (max over a superset of the grid)");
    }

    // --- a STEADY constant reads EXACTLY (each polyphase branch has unity DC gain). The onset 0→0.5 is a
    //     real bandlimited step → it legitimately overshoots, so we measure a settled block, not the onset. ---
    test::group ("TruePeak of a steady constant is exact");
    {
        TPM m; m.prepare (sr, 1024, 1);
        std::vector<float> x (2048, 0.5f); const float* io[1] { x.data() };
        m.process (io, 1, 2048);                                       // fill the ring (absorb the onset step)
        m.process (io, 1, 2048);                                       // steady block — no step inside it
        test::approx (m.truePeakDbBlock(), core::gainToDb (0.5), 0.02, "steady constant 0.5 → exactly −6.02 dBTP (unity-DC phases)");
    }

    // --- pass-band is flat: full-scale sines from low to near-Nyquist all read ≈ 0 dBTP ---
    test::group ("TruePeak pass-band is flat (no droop)");
    {
        test::ok (std::fabs (tonePeakDb (1000.0,  1.0, 0.3, 8000, sr)) < 0.1, "1 kHz full-scale → ~0 dBTP (±0.1 dB)");
        test::ok (std::fabs (tonePeakDb (10000.0, 1.0, 0.3, 8000, sr)) < 0.2, "10 kHz full-scale → ~0 dBTP");
        test::ok (tonePeakDb (21000.0, 1.0, 0.3, 8000, sr) > -1.0,           "21 kHz (near Nyquist) → within ~1 dB (no large droop)");
    }

    // --- amplitude maps to dB correctly ---
    test::group ("TruePeak dB scaling");
    {
        test::approx (tonePeakDb (1000.0, 0.5,  0.0, 8000, sr), -6.0206, 0.1, "0.5 amplitude → −6.02 dBTP");
        test::approx (tonePeakDb (1000.0, 0.25, 0.0, 8000, sr), -12.041, 0.1, "0.25 amplitude → −12.04 dBTP");
    }

    // --- oversampling factor follows the source rate (BS.1770 analysis-rate targets) ---
    test::group ("TruePeak oversampling factor by rate");
    {
        TPM a, b, c, d;
        a.prepare (44100.0, 512, 1); b.prepare (48000.0, 512, 1); c.prepare (96000.0, 512, 1); d.prepare (192000.0, 512, 1);
        test::ok (a.oversampleFactor() == 4 && b.oversampleFactor() == 4, "44.1/48 kHz → 4×");
        test::ok (c.oversampleFactor() == 2, "96 kHz → 2×");
        test::ok (d.oversampleFactor() == 1, "192 kHz → 1× (grid already resolves the peak)");
    }

    // --- stereo reports the max across channels ---
    test::group ("TruePeak stereo = max across channels");
    {
        TPM m; m.prepare (sr, 1024, 2);
        const int N = 4000; std::vector<float> l (N), r (N);
        for (int i = 0; i < N; ++i) { l[i] = (float) (0.1 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); r[i] = (float) (0.8 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); }
        for (int o = 0; o < N; o += 1024) { const float* io[2] { l.data() + o, r.data() + o }; m.process (io, 2, std::min (1024, N - o)); }
        test::approx (m.truePeakDb(), core::gainToDb (0.8), 0.2, "loud R channel sets the reading (max across channels)");
    }

    // --- the display ballistic holds then decays; the authoritative max does NOT decay ---
    test::group ("TruePeak ballistics (display only)");
    {
        TPM m; analysis::TruePeakMeterParams p; p.holdMs = 50.0; p.decayDbPerSec = 100.0; m.prepare (sr, 1024, 1); m.setParams (p);
        const int burst = 1000, tail = 40000;
        std::vector<float> y (burst + tail, 0.0f); for (int i = 0; i < burst; ++i) y[i] = 0.9f;   // loud burst then silence
        for (int o = 0; o < (int) y.size(); o += 1024) { const float* io[1] { y.data() + o }; m.process (io, 1, std::min (1024, (int) y.size() - o)); }
        test::ok (m.truePeakDb() > core::gainToDb (0.9) - 0.5, "authoritative true peak captured the burst and HOLDS it (never decays)");
        test::ok (m.displayTruePeakDb() < m.truePeakDb() - 3.0, "display ballistic has decayed well below the captured peak");
    }

    // --- reset clears the maxima ---
    test::group ("TruePeak reset");
    {
        TPM m; m.prepare (sr, 1024, 1);
        std::vector<float> x (1024, 0.7f); const float* io[1] { x.data() }; m.process (io, 1, 1024);
        m.reset();
        test::ok (m.truePeakDb() < -100.0 && m.samplePeakDb() < -100.0, "after reset → −∞-ish (no peak held)");
    }

    // --- no allocation in process(); a non-finite input cannot poison the reading ---
    test::group ("TruePeak no-alloc + finite guard");
    {
        TPM m; m.prepare (sr, 1024, 2);
        std::vector<float> l (512, 0.3f), r (512, -0.2f); const float* io[2] { l.data(), r.data() };
        m.process (io, 2, 512);
        const long before = g_allocs.load();
        m.process (io, 2, 512); m.process (io, 2, 512);
        const bool noAlloc = (g_allocs.load() == before);
        l[100] = std::nanf (""); r[200] = INFINITY; m.process (io, 2, 512);
        test::okNoAlloc (noAlloc, "process() did not allocate");
        test::ok (std::isfinite (m.truePeakDb()), "NaN/inf input → finite reading (guarded)");
    }

    // --- frequency × phase sweep: quantify over-read (filter ripple/images) and under-read (4× grid) vs 0 dBTP ---
    test::group ("TruePeak over/under-read sweep (BS.1770 envelope)");
    {
        double worstUnder = 0.0, worstOver = -100.0, passbandFloor = 0.0;
        for (double fr = 0.02; fr < 0.49; fr += 0.01)
        {
            const double f = fr * sr;
            double bestPhase = -100.0, worstPhase = 100.0;
            for (int ph = 0; ph < 16; ++ph)
            {
                const double db = tonePeakDbFaded (f, 1.0, 2.0 * core::kPi * ph / 16.0, 6000, sr);
                bestPhase = std::max (bestPhase, db); worstPhase = std::min (worstPhase, db);
            }
            worstOver  = std::max (worstOver, bestPhase);            // highest a unit tone ever reads = over-read
            worstUnder = std::min (worstUnder, worstPhase);          // lowest reading = under-read
            if (fr <= 0.45) passbandFloor = std::min (passbandFloor, bestPhase);   // best-phase droop in the pass-band
        }
        std::printf ("    [sweep] over-read=%+.3f dB  under-read=%+.3f dB  passband best-phase floor=%+.3f dB\n", worstOver, worstUnder, passbandFloor);
        test::ok (worstOver  <  0.3,  "over-read ≤ +0.3 dB (no excessive ripple / image leakage for tones)");
        test::ok (worstUnder > -0.8,  "under-read ≥ −0.8 dB (within the inherent 4× grid limit, ~−0.69 dB at Nyquist)");
        test::ok (passbandFloor > -0.15, "pass-band flat: best-phase reading ≥ −0.15 dB out to 0.45·fs (no droop)");
    }

    // --- the 2× path (96 kHz) recovers the inter-sample peak too ---
    test::group ("TruePeak 2x path at 96 kHz");
    {
        const double sr2 = 96000.0;
        TPM m; m.prepare (sr2, 1024, 1);
        const int N = 8000; std::vector<float> x (N);
        for (int i = 0; i < N; ++i) x[i] = (float) std::sin (2.0 * core::kPi * (sr2 / 4.0) * i / sr2 + core::kPi / 4.0);
        for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
        std::printf ("    [96k] factor=%d truePeak=%+.3f samplePeak=%+.3f\n", m.oversampleFactor(), m.truePeakDb(), m.samplePeakDb());
        test::ok (m.oversampleFactor() == 2, "96 kHz selects 2×");
        test::ok (m.truePeakDb() - m.samplePeakDb() > 1.5, "2× recovers a chunk of the fs/4 inter-sample crest (above sample peak)");
    }

    // --- the 1× path (192 kHz) is sample-peak identical (the grid already resolves the crest) ---
    test::group ("TruePeak 1x path at 192 kHz");
    {
        TPM m; m.prepare (192000.0, 1024, 1);
        unsigned long s = 8; auto rng = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
        const int N = 4000; std::vector<float> x (N); for (int i = 0; i < N; ++i) x[i] = 0.7f * rng();
        for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
        test::ok (m.oversampleFactor() == 1, "192 kHz selects 1× (no oversampling)");
        test::ok (m.truePeakDb() == m.samplePeakDb(), "1× → true peak == sample peak (no interpolation)");
    }

    // --- a full-scale square wave (high inter-sample content) reads a finite, plausible dBTP ---
    test::group ("TruePeak of a full-scale square");
    {
        TPM m; m.prepare (sr, 1024, 1);
        const int N = 8000; std::vector<float> x (N);
        for (int i = 0; i < N; ++i) x[i] = ((i / 8) % 2) ? 1.0f : -1.0f;     // ~3 kHz square, ±1 full-scale
        for (int o = 0; o < N; o += 1024) { const float* io[1] { x.data() + o }; m.process (io, 1, std::min (1024, N - o)); }
        test::ok (std::isfinite (m.truePeakDb()) && m.truePeakDb() > 0.2 && m.truePeakDb() < 3.0, "square edges overshoot to a plausible +0.2..+3 dBTP (finite, no overflow)");
    }

    return test::report();
}
