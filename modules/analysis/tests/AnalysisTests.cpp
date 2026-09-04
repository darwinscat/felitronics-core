// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::analysis::SpectrumTap (the lock-free SPSC handshake).

#include <felitronics_test.h>
#include <felitronics/analysis/SpectrumTap.h>
#include <felitronics/analysis/RollingSpectrumTap.h>
#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/CorrelationMeter.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
#include <chrono>

using namespace felitronics;

// Different orders are distinct, correctly-sized types (Law 5: configurable size without a fork).
static_assert (analysis::SpectrumTapT<10>::kSize == 1024, "order 10 -> 1024");
static_assert (analysis::SpectrumTap::kSize == 2048,       "default -> 2048");

// The rolling tap sizes to the MAX window; it serves any smaller order from the one ring.
static_assert (analysis::RollingSpectrumTap::kMaxSize == 16384,  "default MaxOrder 14 -> 16384");
static_assert (analysis::RollingSpectrumTapT<11>::kMaxSize == 2048, "MaxOrder 11 -> 2048");

volatile double g_perfSink = 0.0;   // keeps the timed loop from being optimised away

int main()
{
    std::printf ("felitronics::analysis tests\n");
    const int N = analysis::kSpectrumFftSize;
    static float dst[analysis::kSpectrumFftSize];   // static: keep the 8 KB off the stack

    // --- SPSC handshake: a frame becomes ready only on the wrap push, and tryPull re-arms ---
    test::group ("SpectrumTap SPSC handshake");
    {
        analysis::SpectrumTap tap;
        test::ok (! tap.tryPull (dst), "no frame before a full window");

        for (int i = 0; i < N; ++i) tap.push ((float) i);   // fills fifo[0..N-1], idx = N
        test::ok (! tap.tryPull (dst), "not ready until the wrap push");

        tap.push (-1.0f);                                    // wrap → snapshots [0..N-1]
        test::ok (tap.tryPull (dst), "frame ready after the wrap push");
        test::approx (dst[0],     0.0,            1e-9, "data[0] == 0");
        test::approx (dst[N - 1], (double) (N - 1), 1e-9, "data[N-1] == N-1");
        test::ok (! tap.tryPull (dst), "re-armed (consumed) after pull");
    }

    // --- no overwrite while a frame is unread: the reader sees the OLD frame, not the new fill ---
    test::group ("SpectrumTap preserves an unread frame");
    {
        analysis::SpectrumTap t2;
        for (int i = 0; i < N; ++i) t2.push (1.0f);
        t2.push (1.0f);                                      // ready, data == all 1.0
        for (int i = 0; i < N; ++i) t2.push (2.0f);          // another window, but unread → snapshot skipped
        test::ok (t2.tryPull (dst), "frame still ready");
        test::approx (dst[0], 1.0, 1e-9, "unread frame preserved (1.0, not overwritten by 2.0)");
    }

    // --- reset() clears the cursor + ready flag ---
    test::group ("SpectrumTap reset");
    {
        analysis::SpectrumTap t3;
        for (int i = 0; i < N + 1; ++i) t3.push (3.0f);      // make a frame ready
        t3.reset();
        test::ok (! t3.tryPull (dst), "reset clears ready");
    }

    //==========================================================================
    // RollingSpectrumTap — hop decoupled from window size; one 8192 ring serves any order.
    static float rdst[analysis::RollingSpectrumTap::kMaxSize];   // pull buffer (>= kMaxSize), static: off-stack
    int ro = -1;                                                 // reported order out-param
    // Push `count` samples of a ramp starting at value `start` (chronology check: data[i] == start+i).
    auto ramp = [] (auto& tap, long long count, long long start) noexcept {
        for (long long i = 0; i < count; ++i) tap.push ((float) (start + i));
    };

    // --- a frame carries the requested window N, the most-recent N samples in chronological order,
    //     and the order it was captured at; the first request is force-published ---
    test::group ("RollingSpectrumTap variable-N snapshot");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 1024, 0);                                    // wpos = 1024
        tap.publishIfDue (10, 512);                             // order 10 (N=1024), first → forced
        test::ok (tap.tryPull (rdst, ro), "first request is force-published");
        test::ok (ro == 10, "frame reports the captured order");
        test::approx (rdst[0],    0.0,    1e-9, "oldest of the last 1024 == 0");
        test::approx (rdst[1023], 1023.0, 1e-9, "newest == 1023");
        test::ok (! tap.tryPull (rdst, ro), "re-armed after pull");
    }

    // --- no frame until the ring holds a FULL window of history (no short/garbage first frame).
    //     Exercised at the MAX order (14 = 16384) — the largest window the ring must fill. ------------
    test::group ("RollingSpectrumTap waits for a full window");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 16000, 0);                                   // < 16384
        tap.publishIfDue (14, 1);                               // N=16384 > wpos → no-op
        test::ok (! tap.tryPull (rdst, ro), "no frame until >= N samples exist");
        ramp (tap, 400, 16000);                                 // wpos = 16400 >= 16384
        tap.publishIfDue (14, 1);
        test::ok (tap.tryPull (rdst, ro) && ro == 14, "frame once a full 16384 window exists");
    }

    // --- hop gates the cadence: no new frame until `hop` samples have accrued since the last publish ---
    test::group ("RollingSpectrumTap hop cadence");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 1024, 0);
        tap.publishIfDue (10, 512); test::ok (tap.tryPull (rdst, ro), "first frame");   // lastPublish = 1024
        ramp (tap, 256, 1024);                                  // wpos = 1280, delta 256 < 512
        tap.publishIfDue (10, 512);
        test::ok (! tap.tryPull (rdst, ro), "no frame before the hop elapses");
        ramp (tap, 256, 1280);                                  // wpos = 1536, delta 512 >= 512
        tap.publishIfDue (10, 512);
        test::ok (tap.tryPull (rdst, ro), "frame once the hop elapses");
    }

    // --- the frame reports the hop that HAPPENED: block rounding and a missed tick stretch it ---
    test::group ("RollingSpectrumTap reports the hop that happened");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        int hop = -1;
        ramp (tap, 1024, 0);
        tap.publishIfDue (10, 512);
        test::ok (tap.tryPull (rdst, ro, hop) && hop == 1024, "first frame: everything since reset (1024)");
        ramp (tap, 600, 1024);                                  // a 600-sample block crosses the 512 hop
        tap.publishIfDue (10, 512);
        test::ok (tap.tryPull (rdst, ro, hop) && hop == 600, "the hop reported is the delta that occurred (600), not the request (512)");
        ramp (tap, 512, 1624); tap.publishIfDue (10, 512);      // published…
        ramp (tap, 512, 2136); tap.publishIfDue (10, 512);      // …the reader is late, so this one is skipped
        test::ok (tap.tryPull (rdst, ro, hop) && hop == 512, "the frame the late reader takes reports its own hop (512)");
        ramp (tap, 100, 2648); tap.publishIfDue (10, 512);      // 512 skipped + 100 since the frame that was taken
        test::ok (tap.tryPull (rdst, ro, hop) && hop == 612, "the missed tick shows up in the NEXT frame's hop (612)");
        int ro2 = -1;
        ramp (tap, 512, 2748); tap.publishIfDue (10, 512);
        test::ok (tap.tryPull (rdst, ro2) && ro2 == 10, "the two-argument tryPull still works");
    }

    // --- an order change force-publishes immediately, bypassing the hop (a live resolution switch) ---
    test::group ("RollingSpectrumTap order change forces a frame");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 2048, 0);
        tap.publishIfDue (10, 1000000);                         // order 10, huge hop
        test::ok (tap.tryPull (rdst, ro) && ro == 10, "order 10 frame");
        tap.publishIfDue (11, 1000000);                         // order changed → forced despite the huge hop
        test::ok (tap.tryPull (rdst, ro) && ro == 11, "order change forces an immediate new-size frame");
    }

    // --- the wrap copy is chronologically correct once wpos exceeds the ring capacity. An order-13
    //     (8192) window over a wpos past the 16384 ring end straddles the wrap → the two-part copy. ----
    test::group ("RollingSpectrumTap wrap correctness");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        const long long total = (long long) analysis::RollingSpectrumTap::kMaxSize + 500;   // wpos past ring end
        ramp (tap, total, 0);
        tap.publishIfDue (13, 1);                               // last 8192 samples straddle the ring wrap
        test::ok (tap.tryPull (rdst, ro), "frame after wrap");
        test::approx (rdst[0],    (double) (total - 8192), 1e-9, "wrap: oldest of last 8192 correct");
        test::approx (rdst[8191], (double) (total - 1),    1e-9, "wrap: newest correct");
    }

    // --- an unread frame is preserved: a busy slot skips the snapshot, reader still sees the OLD frame ---
    test::group ("RollingSpectrumTap preserves an unread frame");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 1024, 0);
        tap.publishIfDue (10, 1);                               // frame A (data == 0..1023), ready
        ramp (tap, 1024, 1024);                                 // wpos = 2048, but slot unread
        tap.publishIfDue (10, 1);                               // ready → snapshot skipped
        test::ok (tap.tryPull (rdst, ro), "frame still ready");
        test::approx (rdst[0], 0.0, 1e-9, "unread frame A preserved (not overwritten by the newer fill)");
    }

    // --- reset restarts the PRODUCER without revoking the reader: an already-published frame stays
    //     drainable (no torn-frame race with a mid-pull reader), and the restarted ring re-arms the
    //     forced first frame + needs a fresh full window again ---
    test::group ("RollingSpectrumTap reset restarts the producer, keeps the mailbox");
    {
        analysis::RollingSpectrumTap tap; tap.reset();
        ramp (tap, 2048, 0);
        tap.publishIfDue (11, 1);                              // frame ready
        tap.reset();                                           // restart producer; mailbox left to the reader
        test::ok (tap.tryPull (rdst, ro), "reset does NOT revoke an in-flight frame (race-free with a mid-pull reader)");
        tap.publishIfDue (11, 1);                             // wpos == 0 after reset -> warmup gate, no frame yet
        test::ok (! tap.tryPull (rdst, ro), "reset restarted the ring (a full window is needed again)");
        ramp (tap, 2048, 0);
        tap.publishIfDue (11, 1000000);                       // forced despite the huge hop (lastOrder re-armed to -1)
        test::ok (tap.tryPull (rdst, ro), "reset re-armed the forced first frame");
    }

    // --- BS.1770 K-weighting: canonical 48 kHz coefficients + recompute at 44.1 k ---
    test::group ("KWeightingFilter canonical 48 kHz coefficients");
    {
        analysis::KWeightingFilter kw; kw.prepare (48000.0, 2);
        const auto s = kw.shelfCoeffs();
        test::approx (s.b0,  1.5351248595869702, 1e-9, "shelf b0");
        test::approx (s.b1, -2.6916961894063807, 1e-9, "shelf b1");
        test::approx (s.b2,  1.1983928108528501, 1e-9, "shelf b2");
        test::approx (s.a1, -1.6906592931824103, 1e-9, "shelf a1");
        test::approx (s.a2,  0.7324807742158501, 1e-9, "shelf a2");
        const auto h = kw.highpassCoeffs();
        test::approx (h.a1, -1.9900474548339797, 1e-9, "hp a1");
        test::approx (h.a2,  0.9900722503662099, 1e-9, "hp a2");
        analysis::KWeightingFilter kw2; kw2.prepare (44100.0, 2);
        test::ok (std::fabs (kw2.shelfCoeffs().b0 - s.b0) > 1e-4, "coeffs recomputed at 44.1 k (not hardcoded)");
    }

    // --- LAW 8: the state must reach EXACTLY zero on silence, not stick in the subnormals ---
    // The defect this guards is not theoretical and not small: with no flush the two TDF-II biquads decay
    // into the subnormal range after ~2.96 s of digital silence and then STICK — (z1,z2) = (-251u, +249u),
    // u = 2^-1074, maps to itself EXACTLY at 48 kHz — so on a CPU with no hardware FTZ (which is what a
    // browser gives; `wasm-audio` has no FP control register) silence costs 54x more than music.
    // Measured on x86-64 gcc 14.2 -O2, 2 s of signal: 0.299 ms for a tone against 16.077 ms for stuck
    // silence, and 0.384 ms once this flush is in. Apple Silicon runs subnormals at full speed and shows
    // nothing, which is exactly why this survived so long.
    test::group ("KWeightingFilter: law 8 — silence flushes to exact zero");
    {
        const double fs = 48000.0;
        auto exciteThenSilence = [fs] (analysis::KWeightingFilter& f, bool flush, double seconds)
        {
            for (int i = 0; i < (int) fs; ++i)                     // 1 s of tone: a zeroed state fed zeros
                f.process (0, std::sin (6.283185307179586 * 997.0 * i / fs));   // stays zero, so EXCITE first
            const int blocks = (int) (fs * seconds) / 512;
            for (int b = 0; b < blocks; ++b)
            {
                for (int i = 0; i < 512; ++i) f.process (0, 0.0);
                if (flush) f.flushDenormals();                    // exactly what LoudnessMeter does per block
            }
            return f.process (0, 0.0);
        };

        analysis::KWeightingFilter dirty; dirty.prepare (fs, 1);
        const double stuck = exciteThenSilence (dirty, false, 5.0);
        // NEGATIVE CONTROL. Without the flush the output after 5 s of silence is a nonzero SUBNORMAL — if
        // this ever reads 0 the positive test below has stopped proving anything.
        test::ok (stuck != 0.0 && std::fabs (stuck) < 1e-300,
                  "without the flush the state is a stuck subnormal, not zero");

        analysis::KWeightingFilter clean; clean.prepare (fs, 1);
        const double flushed = exciteThenSilence (clean, true, 5.0);
        test::ok (flushed == 0.0, "with the per-block flush the state is EXACTLY zero after 5 s of silence");

        // And it is reached long before the subnormal range is: the flush threshold is 1e-15, ~293 decimal
        // orders above the smallest normal double, so the sticky fixed point is never approached at all.
        analysis::KWeightingFilter quick; quick.prepare (fs, 1);
        test::ok (exciteThenSilence (quick, true, 1.0) == 0.0, "already exactly zero after 1 s");
    }

    // --- the flush must not make the answer depend on the host's block size ---
    // This is the test that killed the first design. Flushing at the end of LoudnessMeter::process() looks
    // equivalent and is not: it puts a numerical event on a boundary the CALLER chooses, so the same stream
    // split differently gives different numbers — and this repo tests the opposite (ProbeTests.cpp:212,
    // bit-exact across call sizes 1 … 100 003). Flushing per 10 ms sub-hop is caller-independent.
    test::group ("LoudnessMeter: the law-8 flush keeps chunk invariance bit-exact");
    {
        const double fs = 48000.0;
        const int n = (int) (fs * 3.0);
        std::vector<float> a ((size_t) n), b ((size_t) n);
        for (int i = 0; i < n; ++i)                      // 1 s tone, 1 s digital silence, 1 s tone: the
        {                                                // silence is what drives the state through 1e-15
            const double t = (double) i / fs;
            const double v = (t < 1.0 || t >= 2.0) ? 0.5 * std::sin (6.283185307179586 * 1000.0 * i / fs) : 0.0;
            a[(size_t) i] = b[(size_t) i] = (float) v;
        }
        const float* chans[2] { a.data(), b.data() };

        auto runChunked = [&] (int chunk)
        {
            analysis::LoudnessMeter m; m.prepare (fs, 2, 20.0);
            for (int off = 0; off < n; off += chunk)
            {
                const int take = (off + chunk <= n) ? chunk : (n - off);
                const float* part[2] { a.data() + off, b.data() + off };
                m.process (part, 2, take);
            }
            return m.integratedLufs();
        };
        (void) chans;
        const double ref = runChunked (n);
        for (const int chunk : { 1, 7, 64, 512, 8192, 100003 })
            test::ok (runChunked (chunk) == ref,
                      "chunk " + std::to_string (chunk) + " is bit-identical to one call");
    }

    // --- the flush cadence must beat the SHELF, at every rate the filter is valid for ---
    // The number that binds is not the RLB stage's 2.85 s but the shelf's 90 ms: pole 0.856 against 0.995,
    // so from the 1e-15 threshold the shelf reaches the subnormal floor in 3 983 samples at 44.1 kHz and
    // 4 324 at 48 kHz — 90 ms at every rate, since both scale together. A per-host-block flush at
    // kMaxBlockSize (8192 = 170 ms at 48 kHz) would therefore be too slow by 2x; the meter's 10 ms sub-hop
    // is 480 samples at 48 kHz, a 9x margin. This test drives the PRODUCTION path at four rates.
    test::group ("LoudnessMeter: the sub-hop flush beats the shelf stage at every rate");
    {
        for (const double fs : { 44100.0, 48000.0, 96000.0, 192000.0 })
        {
            analysis::KWeightingFilter f; f.prepare (fs, 1);
            const int sub = (int) std::lround (0.01 * fs);        // the meter's 10 ms sub-hop
            for (int i = 0; i < (int) fs; ++i)                    // 1 s of tone — a zeroed state fed zeros
                f.process (0, std::sin (6.283185307179586 * 997.0 * i / fs));   // stays zero, so excite first
            for (int k = 0; k < 500; ++k)                         // 5 s of silence, flushed per sub-hop
            { for (int i = 0; i < sub; ++i) f.process (0, 0.0); f.flushDenormals(); }
            test::ok (f.process (0, 0.0) == 0.0,
                      "exactly zero after 5 s of silence at " + std::to_string ((int) fs) + " Hz");
        }
    }

    // --- LAW 8, the part only a CLOCK can see: silence must not cost more than sound ---
    // This is the check that actually matters on the tier the law exists for, and it is written as a RATIO
    // on purpose: an absolute millisecond budget would be a lottery on a shared CI runner, while silence
    // and music are measured back to back so load hits both. On arm64 it passes by construction (Apple
    // Silicon runs subnormals at full speed) and proves nothing there — it earns its keep on the x86 rows,
    // where removing the flush takes the ratio from ~1.3 to ~54 (measured, i9-13900H, gcc 14.2 -O2).
    // Bound 10x: a factor of 4 below the real defect and 7 above the healthy state.
    test::group ("KWeightingFilter: silence is not more expensive than sound (law 8, x86 rows)");
    {
        const double fs = 48000.0;
        const int block = 512, seconds = 1;
        auto timeIt = [fs, block, seconds] (bool silence)
        {
            analysis::KWeightingFilter f; f.prepare (fs, 2);
            for (int i = 0; i < (int) fs; ++i) f.process (0, std::sin (6.283185307179586 * 997.0 * i / fs));
            for (int i = 0; i < (int) (fs * 4); ++i) f.process (0, 0.0);   // reach the sticky region
            std::vector<double> buf ((size_t) block);
            for (int i = 0; i < block; ++i)
                buf[(size_t) i] = silence ? 0.0 : 0.25 * std::sin (6.283185307179586 * 440.0 * i / fs);
            double best = 1e300;
            for (int rep = 0; rep < 3; ++rep)
            {
                const auto t0 = std::chrono::steady_clock::now();
                double acc = 0.0;
                for (int off = 0; off < (int) (fs * seconds); off += block)
                {
                    for (int i = 0; i < block; ++i) acc += f.process (0, buf[(size_t) i]);
                    f.flushDenormals();
                }
                const auto t1 = std::chrono::steady_clock::now();
                g_perfSink += acc;                                    // escape, or the loop is dead code
                best = std::min (best, std::chrono::duration<double, std::milli> (t1 - t0).count());
            }
            return best;
        };
        const double music = timeIt (false), silence = timeIt (true);
        const double ratio = music > 0.0 ? silence / music : 1.0;
        std::printf ("      [music %.3f ms, silence %.3f ms, ratio %.2fx]\n", music, silence, ratio);
        test::ok (ratio < 10.0, "silence costs less than 10x sound (ratio " + std::to_string (ratio) + ")");
    }

    // --- the flush must not move a reading: it only ever zeroes what is already far below the gate ---
    test::group ("KWeightingFilter: the flush changes no measured loudness");
    {
        const double fs = 48000.0;
        const int n = (int) fs;                                   // 1 s
        std::vector<float> l ((size_t) n), r ((size_t) n);
        for (int i = 0; i < n; ++i)
            l[(size_t) i] = r[(size_t) i] = (float) (0.5 * std::sin (6.283185307179586 * 1000.0 * i / fs));
        const float* chans[2] { l.data(), r.data() };

        analysis::LoudnessMeter m; m.prepare (fs, 2, 10.0);
        m.process (chans, 2, n);
        const double lufs = m.integratedLufs();
        // -3.01 LUFS is the analytic answer for a 1 kHz sine at 0.5 in both channels: a 1 kHz tone sits at
        // K-weighting unity, so LUFS = 10*log10(2 * 0.5^2/2) - 0.691 + ... — the conformance suite pins the
        // exact value; here the point is only that a flushed meter still produces a sane finite reading.
        test::ok (std::isfinite (lufs) && lufs < 0.0 && lufs > -10.0,
                  "a flushed meter still reads a sane integrated loudness (" + std::to_string (lufs) + " LUFS)");
    }

    // --- correlation: identical → +1, inverted → -1, quadrature → ~0 ---
    test::group ("CorrelationMeter");
    {
        const double sr = 48000.0, f = 500.0;
        analysis::CorrelationMeter cm; cm.prepare (sr, 50.0);
        for (int i = 0; i < (int) sr; ++i) { const float v = (float) std::sin (2.0 * core::kPi * f * i / sr); cm.process (v, v); }
        test::approx (cm.correlation(), 1.0, 0.01, "identical → +1");
        cm.reset();
        for (int i = 0; i < (int) sr; ++i) { const float v = (float) std::sin (2.0 * core::kPi * f * i / sr); cm.process (v, -v); }
        test::approx (cm.correlation(), -1.0, 0.01, "inverted → -1");
        cm.reset();
        for (int i = 0; i < (int) sr; ++i) { const float a = (float) std::sin (2.0 * core::kPi * f * i / sr), b = (float) std::cos (2.0 * core::kPi * f * i / sr); cm.process (a, b); }
        test::ok (std::fabs (cm.correlation()) < 0.1, "quadrature → ~0");
    }

    // --- loudness: a steady tone has momentary ≈ short-term ≈ integrated (gating doesn't skew it) ---
    test::group ("LoudnessMeter steady-tone consistency");
    {
        const double sr = 48000.0; const int n = (int) (5.0 * sr);
        std::vector<float> L (n), R (n);
        for (int i = 0; i < n; ++i) { const float v = (float) (0.5 * std::sin (2.0 * core::kPi * 1000.0 * i / sr)); L[i] = v; R[i] = v; }
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 10.0);
        const float* ch[2] { L.data(), R.data() };
        lm.process (ch, 2, n);
        const double m = lm.momentaryLufs(), st = lm.shortTermLufs(), ig = lm.integratedLufs();
        test::ok (std::isfinite (m) && m > -60.0, "momentary finite + sensible");
        test::approx (ig, m,  0.3, "integrated ≈ momentary for a steady tone");
        test::approx (st, m,  0.3, "short-term ≈ momentary for a steady tone");
        test::ok (lm.loudnessRangeLu() < 1.0, "LRA ≈ 0 for a steady tone (no loudness variation)");
    }

    // --- LRA reflects the spread of a two-level program (EBU Tech 3342: P95 − P10 of gated short-term) ---
    test::group ("LoudnessMeter LRA two-level program");
    {
        const double sr = 48000.0; const int seg = (int) (12.0 * sr);              // 12 s loud + 12 s quiet
        std::vector<float> x ((std::size_t) (2 * seg));
        for (int i = 0; i < seg; ++i) x[(std::size_t) i]         = (float) (0.5f  * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int i = 0; i < seg; ++i) x[(std::size_t) (seg + i)] = (float) (0.15f * std::sin (2.0 * core::kPi * 1000.0 * i / sr));   // ~10.5 dB down
        analysis::LoudnessMeter lm; lm.prepare (sr, 1, 30.0);
        const float* ch[1] { x.data() }; lm.process (ch, 1, 2 * seg);
        const double lra = lm.loudnessRangeLu();
        test::ok (lra > 7.0 && lra < 13.0, "loud↔quiet (~10.5 dB apart) → LRA ≈ the spread");
    }

    // --- a sub-3 s program has no short-term block → LRA is 0 by definition ---
    test::group ("LoudnessMeter LRA needs >= 3 s");
    {
        const double sr = 48000.0; const int n = (int) (2.0 * sr);
        std::vector<float> x ((std::size_t) n, 0.3f);
        analysis::LoudnessMeter lm; lm.prepare (sr, 1, 10.0);
        const float* ch[1] { x.data() }; lm.process (ch, 1, n);
        test::ok (lm.loudnessRangeLu() == 0.0, "< 3 s of audio → LRA = 0 (no short-term sample yet)");
    }

    // a two-level program builder (loud `seg` s at 0.5, then quiet `seg` s at `quietAmp`)
    auto twoLevel = [] (double sr, int seg, float quietAmp) {
        std::vector<float> x ((std::size_t) (2 * seg));
        for (int i = 0; i < seg; ++i) x[(std::size_t) i]         = (float) (0.5f     * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        for (int i = 0; i < seg; ++i) x[(std::size_t) (seg + i)] = (float) (quietAmp * std::sin (2.0 * core::kPi * 1000.0 * i / sr));
        return x;
    };

    // --- the relative gate is −20 LU (NOT −10): an 18 dB-down tail SURVIVES → LRA reflects it.
    //     (had the gate been −10 LU, the quiet tail would be trimmed and LRA would collapse toward 0.) ---
    test::group ("LoudnessMeter LRA −20 LU gate keeps wide dynamics");
    {
        const double sr = 48000.0; const int seg = (int) (12.0 * sr);
        auto x = twoLevel (sr, seg, 0.0630f);                                       // ≈ −18 dB quiet section
        analysis::LoudnessMeter lm; lm.prepare (sr, 1, 30.0);
        const float* ch[1] { x.data() }; lm.process (ch, 1, 2 * seg);
        const double lra = lm.loudnessRangeLu();
        test::ok (lra > 13.0 && lra < 22.0, "18 dB spread → LRA ≈ 18 (the quiet tail passes the −20 LU gate)");
    }

    // --- reset is deterministic; mono and dual-mono stereo give the SAME LRA (the +3 LU offset cancels) ---
    test::group ("LoudnessMeter LRA reset determinism + mono==stereo");
    {
        const double sr = 48000.0; const int seg = (int) (12.0 * sr);
        auto x = twoLevel (sr, seg, 0.15f);
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 30.0);
        const float* mono[1] { x.data() }; lm.process (mono, 1, 2 * seg);
        const double lraA = lm.loudnessRangeLu();
        lm.reset();
        lm.process (mono, 1, 2 * seg);
        const double lraB = lm.loudnessRangeLu();
        test::ok (lraB == lraA, "reset → identical LRA on the same input (deterministic)");
        lm.reset();
        const float* stereo[2] { x.data(), x.data() }; lm.process (stereo, 2, 2 * seg);   // dual-mono = +3 LU everywhere
        test::ok (std::fabs (lm.loudnessRangeLu() - lraA) < 0.2, "mono == dual-mono stereo LRA (a constant offset cancels in P95−P10)");
    }

    // --- lifecycle/misuse: process() before prepare() (empty hopRing) + fs<=0 (hopSamples 0 → /0 in finishHop) ---
    test::group ("LoudnessMeter: reject process before prepare; survive fs<=0");
    {
        analysis::LoudnessMeter lm;                                  // NOT prepared (hopRing empty; hopSamples defaults 4800)
        std::vector<float> sig (6000, 0.1f); const float* io[1] { sig.data() };
        lm.process (io, 1, 6000);                                    // > hopSamples → would finishHop into empty hopRing; must no-op
        test::ok (true, "process before prepare did not OOB the hop ring (ASan check)");
        lm.prepare (0.0, 2);                                         // fs<=0 → clamped, hopSamples>=1 (no /0 in finishHop)
        lm.process (io, 1, 6000);
        test::ok (true, "process after fs<=0 prepare did not divide by zero (UBSan check)");
    }

    // --- FALSIFICATION: one non-finite sample must not latch the correlation meter forever ---
    test::group ("CorrelationMeter recovers from a non-finite sample");
    {
        analysis::CorrelationMeter cm; cm.prepare (48000.0, 50.0);
        cm.process (std::numeric_limits<float>::quiet_NaN(), 0.5f);
        for (int i = 0; i < 48000; ++i)
        {
            const float v = (float) std::sin (2.0 * core::kPi * 997.0 * i / 48000.0);
            cm.process (0.5f * v, -0.5f * v);                          // anti-phase → ρ must read −1
        }
        test::approx (cm.correlation(), -1.0, 0.05, "one NaN sample must not latch the meter (anti-phase reads −1)");
    }

    // --- FALSIFICATION: a channel dropped mid-hop must not leak its partial energy into a later hop ---
    test::group ("LoudnessMeter: a dropped channel's partial hop can't leak later");
    {
        const double srr = 48000.0;
        analysis::LoudnessMeter lm; lm.prepare (srr, 2);
        const int half = 2400, hop = 4800;
        std::vector<float> zero ((std::size_t) hop, 0.0f), zhalf ((std::size_t) half, 0.0f), loud ((std::size_t) half);
        for (int i = 0; i < half; ++i) loud[(std::size_t) i] = (float) (0.45 * std::sin (2.0 * core::kPi * 1000.0 * i / srr));
        const float* st2[2] { zhalf.data(), loud.data() };
        lm.process (st2, 2, half);                                     // half a hop of loud R…
        const float* st0[2] { zero.data(), zero.data() };
        lm.process (st0, 2, 2300);                                     // …zeros IN STEREO so the K-filter ring decays
                                                                       // (the parked hop energy stays parked)…
        const float* mono[1] { zero.data() };
        for (int h = 0; h < 3; ++h) lm.process (mono, 1, hop);         // …then the channel count DROPS mid-hop
        for (int h = 0; h < 4; ++h) lm.process (st0, 2, hop);          // back to stereo: 400 ms of true silence
        // A leaked half-hop of −7 dBFS tone reads ≈ −20 LUFS here; the honest K-filter ring tail is ≈ −80.
        test::ok (lm.momentaryLufs() < -60.0, "momentary after 400 ms of silence ≈ silence (no stale channel energy)");
    }

    // --- BS.1770 absolute reference: a 997 Hz 0 dBFS mono sine reads −3.01 LKFS (both base rates) ---
    test::group ("LoudnessMeter 997 Hz 0 dBFS == −3.01 LUFS (BS.1770 reference)");
    {
        for (double srr : { 48000.0, 44100.0 })
        {
            analysis::LoudnessMeter lm; lm.prepare (srr, 1);
            const int n = (int) (5.0 * srr);
            std::vector<float> x ((std::size_t) n);
            for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) std::sin (2.0 * core::kPi * 997.0 * i / srr);
            const float* ch[1] { x.data() };
            lm.process (ch, 1, n);
            test::approx (lm.integratedLufs(), -3.01, 0.1, srr == 48000.0 ? "997 Hz @48k → −3.01 LKFS" : "997 Hz @44.1k → −3.01 LKFS");
        }
    }

    // --- steady DC is fully rejected by the RLB high-pass (K-weighting), not metered as loudness ---
    test::group ("LoudnessMeter steady DC is K-weighted away");
    {
        const double srr = 48000.0;
        analysis::LoudnessMeter lm; lm.prepare (srr, 1);
        std::vector<float> x ((std::size_t) (3.0 * srr), 0.5f);
        const float* ch[1] { x.data() };
        lm.process (ch, 1, (int) x.size());
        test::ok (lm.momentaryLufs() < -100.0, "steady 0.5 DC → momentary at the silence floor");
    }

    return test::report();
}
