// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the linear-phase EQ. The decisive properties:
//   (1) a FLAT EQ renders an EXACT unit impulse at the centre tap N/2 — proves the build needs no run-time
//       gain hack (the 1/N inverse + the centre-1.0 Blackman-Harris suffice); (2) the FIR is symmetric →
//       exactly linear phase; (3) the realised magnitude matches the Eq bank's target (within window ripple);
//       (4) M/S builds independent Mid & Side IRs; (5) end-to-end the convolution delays by N/2 at unit gain
//       and applies the bell's gain to a tone; (6) process() never allocates; latency = N/2.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/core/Math.h>
#include <felitronics/lineareq/LinearPhaseEq.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace felitronics;
using LPE = lineareq::LinearPhaseEq;

// P88. Level of a 1 kHz sine through `eq`, channel by channel, after the latency and the FIR have settled.
// `leftOnly` feeds the sine to channel 0 and silence to the rest (an M/S check). Returns dB per channel.
template <class Eq>
static std::vector<double> sineGainsDb (Eq& eq, int nch, int firLen, bool leftOnly)
{
    const double fs = 48000.0, f = 1000.0;                        // 48 samples a period
    const int settle = eq.latencySamples() + firLen;
    const int n = settle + 48 * 200;
    std::vector<std::vector<float>> io ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f));
    for (int c = 0; c < nch; ++c)
        if (c == 0 || ! leftOnly)
            for (int i = 0; i < n; ++i) io[(std::size_t) c][(std::size_t) i] = (float) (0.25 * std::sin (2.0 * core::kPi * f * i / fs));
    std::vector<float*> ptr ((std::size_t) nch);
    for (int done = 0; done < n; )
    {
        const int m = std::min (512, n - done);
        for (int c = 0; c < nch; ++c) ptr[(std::size_t) c] = io[(std::size_t) c].data() + done;
        felitronics::test::run (eq.process (ptr.data(), nch, m));
        done += m;
    }
    const double rmsIn = 0.25 / std::sqrt (2.0);
    std::vector<double> g;
    for (int c = 0; c < nch; ++c)
    {
        double e = 0.0;
        for (int i = settle; i < n; ++i) e += (double) io[(std::size_t) c][(std::size_t) i] * io[(std::size_t) c][(std::size_t) i];
        g.push_back (20.0 * std::log10 (std::sqrt (e / (double) (n - settle)) / rmsIn + 1e-30));
    }
    return g;
}

// P88. Noise through `eq` in 512-sample blocks; the output is kept only when `keep` is given.
template <class Eq>
static void noiseThrough (Eq& eq, int nch, int n, float amp, unsigned seed, std::vector<float>* keep)
{
    std::vector<std::vector<float>> io ((std::size_t) nch, std::vector<float> (512, 0.0f));
    std::vector<float*> ptr ((std::size_t) nch);
    for (int done = 0; done < n; )
    {
        const int m = std::min (512, n - done);
        for (int c = 0; c < nch; ++c)
            for (int i = 0; i < m; ++i) { seed = seed * 1664525u + 1013904223u; io[(std::size_t) c][(std::size_t) i] = amp * ((float) (seed >> 9) / 4194304.0f - 1.0f); }
        for (int c = 0; c < nch; ++c) ptr[(std::size_t) c] = io[(std::size_t) c].data();
        felitronics::test::run (eq.process (ptr.data(), nch, m));
        if (keep != nullptr)
            for (int c = 0; c < nch; ++c) keep->insert (keep->end(), io[(std::size_t) c].begin(), io[(std::size_t) c].begin() + m);
        done += m;
    }
}

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
        lp.buildFir (flat, 0, eq::Axis::Mid, firMid.data());
        double maxOther = 0.0; for (int i = 0; i <= N; ++i) if (i != N / 2) maxOther = std::max (maxOther, (double) std::fabs (firMid[(std::size_t) i]));
        test::approx (firMid[(std::size_t) (N / 2)], 1.0, 1e-3, "centre tap == 1.0 (unit gain, no firScale needed)");
        test::ok (maxOther < 1e-3, "every other tap ≈ 0 (a clean delta)");
    }

    // --- (2) the FIR is symmetric → EXACTLY linear phase ---
    test::group ("LinearPhaseEq FIR is symmetric (linear phase)");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell;
        b[0].lane (eq::Lane::Stereo).freq = 1200.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 8.0;
        lp.buildFir (b, 1, eq::Axis::Mid, firMid.data());
        double asym = 0; for (int i = 0; i <= N; ++i) asym = std::max (asym, (double) std::fabs (firMid[(std::size_t) i] - firMid[(std::size_t) (N - i)]));
        test::ok (asym < 1e-6, "fir[N/2+d] == fir[N/2−d] → symmetric → linear phase");
    }

    // --- (3) the realised magnitude matches the bank's target (within window ripple) ---
    test::group ("LinearPhaseEq magnitude matches the EQ target");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell;
        b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        lp.buildFir (b, 1, eq::Axis::Mid, firMid.data());
        test::approx (firMagDb (firMid, 1000.0, sr), 6.0, 0.7, "+6 dB bell → ~+6 dB at 1 kHz");
        test::approx (firMagDb (firMid,  100.0, sr), 0.0, 0.3, "flat away from the bell (~0 dB at 100 Hz)");
        const double target = eq::EqEngine::magnitudeDbFor (b, 1, 1000.0, sr);   // == the bank's own curve
        test::approx (firMagDb (firMid, 1000.0, sr), target, 0.7, "FIR magnitude tracks EqEngine::magnitudeDbFor");
    }

    // --- (4) M/S: independent Mid & Side IRs from one snapshot ---
    test::group ("LinearPhaseEq M/S builds independent Mid+Side IRs");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell;   // {m,s} split point (shared type)
        b[0].lane (eq::Lane::Stereo).on = false;
        b[0].lane (eq::Lane::Mid).on  = true; b[0].lane (eq::Lane::Mid).freq  = 1000.0; b[0].lane (eq::Lane::Mid).Q  = 2.0; b[0].lane (eq::Lane::Mid).gainDb  =  6.0;
        b[0].lane (eq::Lane::Side).on = true; b[0].lane (eq::Lane::Side).freq = 1000.0; b[0].lane (eq::Lane::Side).Q = 2.0; b[0].lane (eq::Lane::Side).gainDb = -6.0;
        lp.buildFir (b, 1, eq::Axis::Mid,  firMid.data());
        lp.buildFir (b, 1, eq::Axis::Side, firSide.data());
        test::approx (firMagDb (firMid,  1000.0, sr),  6.0, 0.7, "Mid IR = +6 dB at 1 kHz");
        test::approx (firMagDb (firSide, 1000.0, sr), -6.0, 0.7, "Side IR = −6 dB at 1 kHz");
    }

    // --- (5a) end-to-end: a flat EQ delays an impulse by N/2 at unit gain ---
    test::group ("LinearPhaseEq end-to-end flat → delayed impulse");
    {
        LPE e; e.prepare (sr, 512, 2, Q); e.setBands (flat, 0);       // unit-impulse IR (crossfades in over ~40 ms)
        const int M = 12000, imp = 6000;                              // impulse well past the crossfade
        std::vector<float> L (M, 0.0f), R (M, 0.0f); L[(std::size_t) imp] = 1.0f; R[(std::size_t) imp] = 1.0f;
        for (int o = 0; o < M; o += 512) { float* io[2] { L.data() + o, R.data() + o }; felitronics::test::run (e.process (io, 2, std::min (512, M - o))); }
        int peakIdx = 0; double peak = 0; for (int i = 0; i < M; ++i) if (std::fabs (L[(std::size_t) i]) > peak) { peak = std::fabs (L[(std::size_t) i]); peakIdx = i; }
        test::ok (peakIdx == imp + N / 2, "output impulse lands at input + N/2 (linear-phase latency)");
        test::approx (peak, 1.0, 0.02, "unit gain (flat EQ passes through, just delayed)");
    }

    // --- (5b) end-to-end: a +6 dB bell lifts a 1 kHz tone ~6 dB ---
    test::group ("LinearPhaseEq end-to-end tone gain");
    {
        LPE e; e.prepare (sr, 512, 2, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell;
        b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        e.setBands (b, 1);
        const int M = 16000;
        std::vector<float> L (M), R (M); for (int i = 0; i < M; ++i) { L[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); R[(std::size_t) i] = L[(std::size_t) i]; }
        for (int o = 0; o < M; o += 512) { float* io[2] { L.data() + o, R.data() + o }; felitronics::test::run (e.process (io, 2, std::min (512, M - o))); }
        double inSq = 0, outSq = 0; int from = M - 4000;              // tail: past crossfade + latency
        for (int i = from; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) L[(std::size_t) i] * L[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq / inSq), 6.0, 0.8, "1 kHz tone out ≈ +6 dB through the +6 dB bell");
    }

    // --- (6) process() never allocates (setBands is the message-thread build; process is RT-safe) ---
    test::group ("LinearPhaseEq process no-alloc");
    {
        LPE e; e.prepare (sr, 512, 2, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 2000.0; b[0].lane (eq::Lane::Stereo).gainDb = -4.0; e.setBands (b, 1);
        std::vector<float> L (512, 0.3f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        felitronics::test::run (e.process (io, 2, 512));
        const long long before = alloc::count.load();
        felitronics::test::run (e.process (io, 2, 512)); felitronics::test::run (e.process (io, 2, 512));
        test::okNoAlloc (alloc::count.load() == before, "process() did not allocate");
        test::ok (e.latencySamples() == N / 2, "latencySamples() == N/2");
    }

    // --- (7) MONO path: 1-channel session → the Mid IR only, no crash, correct gain ---
    test::group ("LinearPhaseEq mono path");
    {
        LPE e; e.prepare (sr, 512, 1, Q);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0; e.setBands (b, 1);
        const int M = 16000; std::vector<float> x (M); for (int i = 0; i < M; ++i) x[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int o = 0; o < M; o += 512) { float* io[1] { x.data() + o }; felitronics::test::run (e.process (io, 1, std::min (512, M - o))); }
        double inSq = 0, outSq = 0; for (int i = M - 4000; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) x[(std::size_t) i] * x[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq / inSq), 6.0, 0.8, "mono 1 kHz tone +6 dB through the Mid IR (no crash)");

        // v2 lanes rule: a non-stereo bus runs the ST lane ONLY — a {m}-only point must be TRANSPARENT
        // on mono (the IIR engine is; bank 0 is built from the ST-only composite so the FIR matches it).
        LPE t; t.prepare (sr, 512, 1, Q);
        eq::BandParams mb[1]; mb[0].on = true; mb[0].type = eq::FilterType::Bell; mb[0].lane (eq::Lane::Stereo).on = false;
        mb[0].lane (eq::Lane::Mid).on = true; mb[0].lane (eq::Lane::Mid).freq = 1000.0; mb[0].lane (eq::Lane::Mid).Q = 2.0; mb[0].lane (eq::Lane::Mid).gainDb = 12.0;
        t.setBands (mb, 1);
        std::vector<float> y (M); for (int i = 0; i < M; ++i) y[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int o = 0; o < M; o += 512) { float* io[1] { y.data() + o }; felitronics::test::run (t.process (io, 1, std::min (512, M - o))); }
        double inSq2 = 0, outSq2 = 0; for (int i = M - 4000; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq2 += s * s; outSq2 += (double) y[(std::size_t) i] * y[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq2 / inSq2), 0.0, 0.3, "mono {m}-only point is transparent (ST-only bank 0, matches the IIR engine)");
    }

    // --- (8) an IR swap mid-stream is click-free AND response-correct after it settles (crossfade=N) ---
    test::group ("LinearPhaseEq IR swap settles clean");
    {
        LPE e; e.prepare (sr, 512, 2, Q); e.setBands (flat, 0);       // start flat (unit IR)
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell;
        b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        const int M = 24000; std::vector<float> L (M), R (M); for (int i = 0; i < M; ++i) { L[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); R[(std::size_t) i] = L[(std::size_t) i]; }
        bool swapped = false;
        for (int o = 0; o < M; o += 512)
        {
            if (! swapped && o >= 8000 && e.setBands (b, 1)) swapped = true;   // swap flat→bell mid-stream
            float* io[2] { L.data() + o, R.data() + o }; felitronics::test::run (e.process (io, 2, std::min (512, M - o)));
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
        std::vector<float> fir (16385); e.buildFir (flat, 0, eq::Axis::Mid, fir.data());
        double mo = 0; for (int i = 0; i <= 16384; ++i) if (i != 8192) mo = std::max (mo, (double) std::fabs (fir[(std::size_t) i]));
        test::ok (std::fabs (fir[8192] - 1.0f) < 1e-3 && mo < 1e-3, "flat → unit impulse at N/2 holds at N=16384 too");
    }

    // --- (10) a steep linear-phase high-pass: flat pass-band, strongly (finitely) attenuated stop-band ---
    test::group ("LinearPhaseEq steep high-pass");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::HighPass;
        b[0].lane (eq::Lane::Stereo).freq = 500.0; b[0].lane (eq::Lane::Stereo).Q = 0.707; b[0].lane (eq::Lane::Stereo).slope = 24;
        lp.buildFir (b, 1, eq::Axis::Mid, firMid.data());
        test::ok (firMagDb (firMid, 2000.0, sr) > -1.0,  "pass-band (2 kHz, 2 oct up) ≈ 0 dB");
        test::ok (firMagDb (firMid,  125.0, sr) < -18.0, "stop-band (125 Hz, 2 oct down) strongly attenuated (finite, per the window)");
    }

    // --- lifecycle/misuse: setBands()/process() before prepare() must not touch empty buffers ---
    // Default N_ = 16384, but firMid_/spec_/time_/magBuf_/window_ are empty until prepare(); setBands() wrote
    // N+1 taps into nothing (OOB — silent on ARM64, heap corruption on x86-64 / bare-metal). Guarded now.
    test::group ("LinearPhaseEq: reject setBands/process before prepare");
    {
        lineareq::LinearPhaseEq eqm;                                  // NOT prepared
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        test::ok (! eqm.setBands (b, 1), "setBands() before prepare() returns false (no write into empty FIR)");
        float l[16] {}, r[16] {}; float* io[2] { l, r };
        test::ok (! eqm.process (io, 2, 16), "process() before prepare() is REFUSED (law 11)");
        eqm.reset();                                                  // safe on an unprepared engine
        test::ok (eqm.prepare (48000.0, 16, 2, 0), "prepare() after the rejected calls");
        test::ok (eqm.setBands (b, 1), "setBands() works once prepared");
        felitronics::test::run (eqm.process (io, 2, 16));
    }

    // --- 🔴 P88: a band move published a moment before a host restart is what plays after it --- The product
    // path of law 11e. `reset()` here forwards to MatrixConvolverNupc::reset() (or, above two channels, to one
    // per channel, staged then published), which used to DROP a publication the audio thread had not picked up:
    // measured through this class at 48 kHz, a +12 dB bell published over a settled flat curve read +0.00 dB
    // after the restart, and when it was the FIRST curve ever published the EQ answered exactly zero output
    // (the fixture's -600 dB floor) until the next band move. No convolver test can see a composite that
    // forgets to forward, or forwards to the wrong verb; this one can.
    test::group ("LinearPhaseEq: a restart keeps the curve published a moment before it (P88)");
    {
        eq::BandParams flatCurve[1];
        eq::BandParams boost[1]; boost[0].on = true; boost[0].type = eq::FilterType::Bell;
        boost[0].lane (eq::Lane::Stereo).freq = 1000.0; boost[0].lane (eq::Lane::Stereo).Q = 1.0; boost[0].lane (eq::Lane::Stereo).gainDb = 12.0;
        for (int nch : { 1, 2, 4 })                                   // 4: the per-channel stage-then-publish path
            for (int firstEver = 0; firstEver <= 1; ++firstEver)
            {
                LPE e;
                test::ok (e.prepare (48000.0, 512, nch, 0), "prepare");
                if (firstEver == 0)
                {
                    test::ok (e.setBands (flatCurve, 1), "the flat curve is accepted");
                    (void) sineGainsDb (e, nch, e.firSize(), false);       // …and fades in and settles
                }
                test::ok (e.setBands (boost, 1), "the +12 dB move is accepted");
                e.reset();                                                  // the host restarts before a block picks it up
                const std::vector<double> g = sineGainsDb (e, nch, e.firSize(), false);
                for (int c = 0; c < nch; ++c)
                    test::approx (g[(std::size_t) c], 12.0, 0.2, std::string (firstEver ? "first curve ever" : "move over a settled curve")
                                                                 + ", restarted: the bell still plays [" + std::to_string (nch)
                                                                 + " ch, channel " + std::to_string (c) + "]");

                // …and from the FIRST sample, not after a fade: law 11a independence. Two EQs fed different
                // audio, one restarted with the move still in flight and one after it settled, answer the
                // next programme with the same bits. A steady-state gain cannot tell reset() from a
                // composite that forwards to clearAudioState() — that verb leaves the fade running and
                // reaches the same gain 20 ms later — and the mutation stand showed exactly that mutant green.
                const auto restarted = [&] (bool settleFirst, float preAmp, unsigned preSeed)
                {
                    LPE x;
                    felitronics::test::run (x.prepare (48000.0, 512, nch, 0));
                    if (firstEver == 0)
                    {
                        test::ok (x.setBands (flatCurve, 1), "flat accepted");
                        noiseThrough (x, nch, 2 * x.firSize(), preAmp, preSeed, nullptr);
                    }
                    test::ok (x.setBands (boost, 1), "move accepted");
                    if (settleFirst) noiseThrough (x, nch, 4096, preAmp, preSeed + 7u, nullptr);   // the 20 ms fade is over
                    x.reset();
                    std::vector<float> out;
                    noiseThrough (x, nch, 8192, 0.3f, 424242u, &out);
                    return out;
                };
                const std::vector<float> inFlight = restarted (false, 0.5f, 11u), settled = restarted (true, 0.2f, 99u);
                bool same = inFlight.size() == settled.size();
                for (std::size_t i = 0; same && i < inFlight.size(); ++i) same = core::sameBits (inFlight[i], settled[i]);
                test::ok (same, std::string (firstEver ? "first curve ever" : "move over a settled curve")
                                + ", restarted in flight == restarted after it settled, every bit [" + std::to_string (nch) + " ch]");
            }
    }
    // …and through the M/S (MSDiag) operator, with a sine on the LEFT only so Mid and Side both carry
    // signal: Mid +6 dB and Side -6 dB decode to 0.5·(gM + gS) on the left and 0.5·(gM − gS) on the right.
    {
        eq::BandParams ms[1]; ms[0].on = true; ms[0].type = eq::FilterType::Bell;
        ms[0].lane (eq::Lane::Stereo).on = false;
        ms[0].lane (eq::Lane::Mid).on  = true; ms[0].lane (eq::Lane::Mid).freq  = 1000.0; ms[0].lane (eq::Lane::Mid).Q  = 1.0; ms[0].lane (eq::Lane::Mid).gainDb  =  6.0;
        ms[0].lane (eq::Lane::Side).on = true; ms[0].lane (eq::Lane::Side).freq = 1000.0; ms[0].lane (eq::Lane::Side).Q = 1.0; ms[0].lane (eq::Lane::Side).gainDb = -6.0;
        const double gM = std::pow (10.0, 6.0 / 20.0), gS = std::pow (10.0, -6.0 / 20.0);
        const double wantL = 20.0 * std::log10 (0.5 * (gM + gS)), wantR = 20.0 * std::log10 (0.5 * (gM - gS));
        eq::BandParams flatCurve[1];
        LPE e;
        test::ok (e.prepare (48000.0, 512, 2, 0), "prepare");
        test::ok (e.setBands (flatCurve, 1), "the flat curve is accepted");
        (void) sineGainsDb (e, 2, e.firSize(), true);
        test::ok (e.setBands (ms, 1), "the M/S move is accepted");
        e.reset();
        const std::vector<double> g = sineGainsDb (e, 2, e.firSize(), true);
        test::approx (g[0], wantL, 0.2, "M/S move, restarted: the left channel reads the decoded Mid+Side");
        test::approx (g[1], wantR, 0.2, "M/S move, restarted: the right channel reads the decoded Mid-Side");
        // …and the same independence as above, through the MSDiag operator with different L and R
        const auto restarted = [&] (bool settleFirst, float preAmp, unsigned preSeed)
        {
            LPE x;
            felitronics::test::run (x.prepare (48000.0, 512, 2, 0));
            test::ok (x.setBands (flatCurve, 1), "flat accepted");
            noiseThrough (x, 2, 2 * x.firSize(), preAmp, preSeed, nullptr);
            test::ok (x.setBands (ms, 1), "M/S move accepted");
            if (settleFirst) noiseThrough (x, 2, 4096, preAmp, preSeed + 7u, nullptr);
            x.reset();
            std::vector<float> out;
            noiseThrough (x, 2, 8192, 0.3f, 424242u, &out);
            return out;
        };
        const std::vector<float> inFlight = restarted (false, 0.5f, 11u), settled = restarted (true, 0.2f, 99u);
        bool same = inFlight.size() == settled.size();
        for (std::size_t i = 0; same && i < inFlight.size(); ++i) same = core::sameBits (inFlight[i], settled[i]);
        test::ok (same, "M/S move, restarted in flight == restarted after it settled, every bit");
    }

    return test::report();
}
