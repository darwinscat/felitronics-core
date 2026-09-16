// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for NaturalPhaseEq — the mixed-phase ("Natural") rendering. Decisive properties:
//   (1) a FLAT EQ renders a unit impulse at `bulkDelay` (unity pass-through at the reported latency);
//   (2) the realised magnitude matches the bank's target (mixed phase preserves |H|);
//   (3) latency = a FIXED L/4 — strictly LESS than linear's L/2, and the same for every k;
//   (4) the FIR is NON-symmetric and front-loaded (mixed phase, not linear), with less pre- than post-energy;
//   (5) process() never allocates (RT-safe) and runs mono + stereo.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/core/Math.h>
#include <felitronics/lineareq/NaturalPhaseEq.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace felitronics;
using NPE = lineareq::NaturalPhaseEq;

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
    const int delay = np.latencySamples();                           // = L/4, fixed for all k
    std::vector<float> firMid ((std::size_t) L), firSide ((std::size_t) L);
    eq::BandParams flat[1];

    // --- (1) flat EQ → a unit impulse at the bulk delay (unity pass-through at the latency) ---
    test::group ("NaturalPhaseEq flat EQ == unit impulse at the bulk delay");
    {
        np.buildFir (flat, 0, eq::Axis::Mid, firMid.data());
        double maxOther = 0.0;
        for (int i = 0; i < L; ++i) if (i != delay) maxOther = std::max (maxOther, (double) std::fabs (firMid[(std::size_t) i]));
        test::approx (firMid[(std::size_t) delay], 1.0, 1e-3, "peak tap == 1.0 at bulkDelay (unity gain)");
        test::ok (maxOther < 1e-3, "every other tap ≈ 0 (a clean delayed delta)");
    }

    // --- (2) realised magnitude matches the bank's target (mixed phase preserves |H|) ---
    test::group ("NaturalPhaseEq magnitude matches the EQ target");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        np.buildFir (b, 1, eq::Axis::Mid, firMid.data());
        test::approx (firMagDb (firMid, 1000.0, sr),  6.0, 0.8, "+6 dB bell → ~+6 dB at 1 kHz");
        test::approx (firMagDb (firMid,  100.0, sr),  0.0, 0.6, "≈ 0 dB well below the bell");
        eq::BandParams c[1]; c[0].on = true; c[0].type = eq::FilterType::Bell; c[0].lane (eq::Lane::Stereo).freq = 3000.0; c[0].lane (eq::Lane::Stereo).Q = 3.0; c[0].lane (eq::Lane::Stereo).gainDb = -8.0;
        np.buildFir (c, 1, eq::Axis::Mid, firMid.data());
        test::approx (firMagDb (firMid, 3000.0, sr), -8.0, 0.8, "−8 dB cut → ~−8 dB at 3 kHz");
    }

    // --- (3) latency is a FIXED L/4 and strictly lighter than linear's L/2 ---
    test::group ("NaturalPhaseEq latency lighter than linear");
    {
        test::ok (delay == L / 4, "k=0.5 → bulk delay == L/4");
        test::ok (delay < L / 2,  "lighter than a linear-phase FIR of the same length (L/2)");
    }

    // --- (4) the FIR is mixed phase: NON-symmetric + front-loaded (less pre- than post-energy) ---
    test::group ("NaturalPhaseEq FIR is mixed phase (not linear/symmetric)");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1200.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 8.0;
        np.buildFir (b, 1, eq::Axis::Mid, firMid.data());
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
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 1.5; b[0].lane (eq::Lane::Stereo).gainDb = 4.0;
        np.setBands (b, 1);                                          // build (message thread, allocates) BEFORE the snapshot
        std::vector<float> lch (512, 0.1f), rch (512, -0.1f);
        float* io[2] { lch.data(), rch.data() };
        felitronics::test::run (np.process (io, 2, 512));                                     // consume the initial fade-in
        const long long before = alloc::count.load();
        felitronics::test::run (np.process (io, 2, 512));
        float* mono[1] { lch.data() };
        // LAW 11(c): a stereo-prepared natural-phase EQ convolves a 2x2 matrix, so a 1-plane call cannot
        // be honoured. It used to be accepted and do NOTHING (the matrix convolver dropped it, void), and
        // this line called that "the mono path". It is a refusal now, and it says so.
        const bool monoRefused = ! np.process (mono, 1, 512);   // recorded here, ASSERTED after the snapshot:
        const long long after = alloc::count.load();   // test::ok builds a std::string and allocates
        test::ok (monoRefused, "a mono call on a STEREO-prepared engine is refused (law 11c)");
        test::okNoAlloc (after == before, "process() performed zero heap allocations (stereo + mono)");
        bool finite = true; for (float v : lch) finite = finite && std::isfinite (v);
        test::ok (finite, "output finite");

        // v2 lanes rule: a mono-prepared instance builds bank 0 from the ST-ONLY composite, so a
        // {m}-only point is TRANSPARENT on mono — matching the IIR engine's non-stereo behaviour.
        NPE t; t.prepare (sr, 512, 1, Q, 0.5f);
        eq::BandParams mb[1]; mb[0].on = true; mb[0].type = eq::FilterType::Bell; mb[0].lane (eq::Lane::Stereo).on = false;
        mb[0].lane (eq::Lane::Mid).on = true; mb[0].lane (eq::Lane::Mid).freq = 1000.0; mb[0].lane (eq::Lane::Mid).Q = 2.0; mb[0].lane (eq::Lane::Mid).gainDb = 12.0;
        t.setBands (mb, 1);
        const int M = 16000; std::vector<float> y ((std::size_t) M);
        for (int i = 0; i < M; ++i) y[(std::size_t) i] = (float) (0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int o = 0; o < M; o += 512) { float* io1[1] { y.data() + o }; felitronics::test::run (t.process (io1, 1, std::min (512, M - o))); }
        double inSq = 0, outSq = 0;
        for (int i = M - 4000; i < M; ++i) { const double s = 0.4 * std::sin (2.0 * core::kPi * 1000.0 * i / sr); inSq += s * s; outSq += (double) y[(std::size_t) i] * y[(std::size_t) i]; }
        test::approx (10.0 * std::log10 (outSq / inSq), 0.0, 0.3, "mono {m}-only point is transparent (ST-only bank 0, matches the IIR engine)");
    }

    // --- (6) STEEP filter accuracy: the kept L taps must still track a steep HP's magnitude
    //     (bulkDelay/L truncation is approximate; this bounds the error on a hard curve) ---
    test::group ("NaturalPhaseEq tracks a steep high-pass");
    {
        NPE hp; hp.prepare (sr, 512, 2, 2, 0.5f);                    // quality 2 → L = 8192 (headroom for ringing)
        std::vector<float> fir ((std::size_t) hp.firSize());
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::HighPass; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 0.707; b[0].lane (eq::Lane::Stereo).slope = 24;
        hp.buildFir (b, 1, eq::Axis::Mid, fir.data());
        test::approx (firMagDb (fir, 4000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 4000.0, sr), 1.0, "pass-band (4 kHz) ≈ target");
        test::approx (firMagDb (fir, 2000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 2000.0, sr), 1.2, "near cutoff (2 kHz) ≈ target");
        test::ok (firMagDb (fir, 250.0, sr) < -14.0, "stop-band (250 Hz) strongly attenuated");
    }

    // --- (7) a deep narrow cut keeps its depth (steep/narrow feature survives the truncation) ---
    test::group ("NaturalPhaseEq keeps a deep narrow cut");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 2000.0; b[0].lane (eq::Lane::Stereo).Q = 6.0; b[0].lane (eq::Lane::Stereo).gainDb = -18.0;
        np.buildFir (b, 1, eq::Axis::Mid, firMid.data());                   // np: L = 4096
        test::approx (firMagDb (firMid, 2000.0, sr), eq::EqEngine::magnitudeDbFor (b, 1, 2000.0, sr), 2.0, "−18 dB narrow cut depth ≈ target");
        test::approx (firMagDb (firMid,  500.0, sr), 0.0, 0.6, "≈ 0 dB away from the cut");
    }

    // --- (8) setBlend changes the phase LIVE but keeps the latency fixed (for a smooth knob, no re-prepare) ---
    test::group ("NaturalPhaseEq setBlend: phase moves, latency fixed");
    {
        NPE n2; n2.prepare (sr, 512, 2, 1, 0.5f);
        const int Ld = n2.firSize(), lat0 = n2.latencySamples();
        std::vector<float> fHalf ((std::size_t) Ld), fMin ((std::size_t) Ld);
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1500.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        n2.buildFir (b, 1, eq::Axis::Mid, fHalf.data());     // k = 0.5
        n2.setBlend (0.95f);
        n2.buildFir (b, 1, eq::Axis::Mid, fMin.data());      // k = 0.95 (near minimum phase)
        test::ok (n2.latencySamples() == lat0, "latency unchanged by the blend (fixed bulk delay → no re-prepare)");
        double diff = 0.0; for (int i = 0; i < Ld; ++i) diff = std::max (diff, (double) std::fabs (fHalf[(std::size_t) i] - fMin[(std::size_t) i]));
        test::ok (diff > 1e-3, "the FIR actually changed (phase blended)");
        double preH = 0.0, preM = 0.0; for (int i = 0; i < lat0; ++i) { preH += std::fabs (fHalf[(std::size_t) i]); preM += std::fabs (fMin[(std::size_t) i]); }
        test::ok (preM < preH, "k=0.95 has less pre-ring than k=0.5 (toward minimum phase)");
    }

    // --- (9) magnitude holds across the WHOLE blend range k∈[0,1], not just k=0.5. A panel review worried
    //     the fixed-L/4 bulk shift (asymmetric: L/4 pre-ring vs 3L/4 post) would clip pre-ring and wreck |H|
    //     at low k. Measured: it does NOT — the re-injected |H| keeps the realised magnitude on target for
    //     every k across the audible band (the only deviations sit in a steep HP's deep, inaudible sub-bass
    //     stop-band, and shrink with quality). This locks that in. ---
    test::group ("NaturalPhaseEq magnitude holds across all k (refutes the low-k truncation worry)");
    {
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1200.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        NPE m;
        for (float kk : { 0.0f, 0.25f, 0.5f, 0.75f, 1.0f })
        {
            m.prepare (sr, 512, 2, 1, kk);
            std::vector<float> fir ((std::size_t) m.firSize());
            m.buildFir (b, 1, eq::Axis::Mid, fir.data());
            double worst = 0.0;
            for (double f = 100.0; f < 16000.0; f *= 1.1)
                worst = std::max (worst, std::fabs (firMagDb (fir, f, sr) - eq::EqEngine::magnitudeDbFor (b, 1, f, sr)));
            char msg[112]; std::snprintf (msg, sizeof msg, "k=%.2f: realised |H| within 0.5 dB of target across the audible band", (double) kk);
            test::ok (worst < 0.5, msg);
        }
    }

    // --- lifecycle/misuse: setBands()/process() before prepare() must not touch empty buffers ---
    // Same shape as LinearPhaseEq, plus setBands() drives an UNPREPARED MixedPhaseFir (mp_.build) → empty
    // spec_/ceps_/h_. All silent on ARM64, OOB on x86-64 / bare-metal. Guarded now.
    test::group ("NaturalPhaseEq: reject setBands/process before prepare");
    {
        lineareq::NaturalPhaseEq eqm;                                 // NOT prepared
        eq::BandParams b[1]; b[0].on = true; b[0].type = eq::FilterType::Bell; b[0].lane (eq::Lane::Stereo).freq = 1000.0; b[0].lane (eq::Lane::Stereo).Q = 2.0; b[0].lane (eq::Lane::Stereo).gainDb = 6.0;
        test::ok (! eqm.setBands (b, 1), "setBands() before prepare() returns false (no write into empty FIR / unprepared MixedPhaseFir)");
        float l[16] {}, r[16] {}; float* io[2] { l, r };
        test::ok (! eqm.process (io, 2, 16), "process() before prepare() is REFUSED (law 11)");
        eqm.reset();
        test::ok (eqm.prepare (48000.0, 16, 2, 0, 0.5f), "prepare() after the rejected calls");
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
    test::group ("NaturalPhaseEq: a restart keeps the curve published a moment before it (P88)");
    {
        eq::BandParams flatCurve[1];
        eq::BandParams boost[1]; boost[0].on = true; boost[0].type = eq::FilterType::Bell;
        boost[0].lane (eq::Lane::Stereo).freq = 1000.0; boost[0].lane (eq::Lane::Stereo).Q = 1.0; boost[0].lane (eq::Lane::Stereo).gainDb = 12.0;
        for (int nch : { 1, 2, 4 })                                   // 4: the per-channel stage-then-publish path
            for (int firstEver = 0; firstEver <= 1; ++firstEver)
            {
                NPE e;
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
                    NPE x;
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

    return test::report();
}
