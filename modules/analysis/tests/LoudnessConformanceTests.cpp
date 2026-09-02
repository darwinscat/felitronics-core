// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// EBU Tech 3341 conformance for felitronics::analysis::LoudnessMeter — the meter against the THEORY of
// ITU-R BS.1770-4, never against the implementation. Tech 3341 (2023) Table 1 defines its test signals
// mathematically — 1 kHz sines at set per-channel peak levels and durations, expected readings ±0.1 LU — so
// every fixture here is synthesized from the spec, at both base rates the family runs at:
//   * cases 1–2: calibration and linearity, for I, M and S alike;
//   * cases 3–5: the −10 LU relative gate, the −70 LUFS absolute gate, power-averaging across levels;
//   * case 6: the 5.0 channel weights (Ls/Rs at 1.41);
//   * cases 9 and 12: short-term and momentary settle to a constant over a periodic program.
// Then the parts a flat RMS meter would fail: the K in K-weighting (a high tone reads hot, a low one reads
// under), the −120 answer for what cannot be measured (silence, sub-gate signal, less than one 400 ms
// block), and chunk invariance (the same stream in odd-sized pieces reads identically).
//
// Grown from the Looper Cat conformance suite that gated the product's move onto this meter.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/analysis/LoudnessMeter.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace felitronics;

namespace
{
    // Tech 3341 Table 1 specifies 1 kHz. (BS.1770's own calibration is stated at 997 Hz, where the
    // K-filter's gain is +0.691 dB exactly; at 1 kHz it is +0.698 dB, a residual well inside ±0.1 LU.)
    constexpr double kToneHz = 1000.0;
    constexpr double kNone   = -120.0;     // the meter's "nothing survived the gates" answer

    using Steps = std::vector<std::pair<double, double>>;   // (dBFS per-channel peak, seconds)

    // Appends `seconds` of a dual-mono sine at `dbfs` peak amplitude, phase-continuous across calls through
    // the caller's running sample counter (the level steps of cases 3–5 must add no click energy), fed in
    // `chunk`-frame slices so the meter sees the block boundaries a real stream produces. `after` hears the
    // seconds fed so far after every slice — the hook the ballistics cases read M and S through.
    template <typename After>
    void feedSine (analysis::LoudnessMeter& lm, double sr, double hz, double dbfs, double seconds,
                   long long& sampleIndex, int chunk, After&& after)
    {
        const double amp = std::pow (10.0, dbfs / 20.0);
        const double w   = 2.0 * core::kPi * hz / sr;
        long long frames = std::llround (seconds * sr);
        std::vector<float> buf ((std::size_t) chunk);
        while (frames > 0)
        {
            const int n = (int) std::min<long long> (frames, chunk);
            for (int i = 0; i < n; ++i) buf[(std::size_t) i] = (float) (amp * std::sin (w * (double) sampleIndex++));
            const float* ch[2] { buf.data(), buf.data() };
            lm.process (ch, 2, n);
            frames -= n;
            after ((double) sampleIndex / sr);
        }
    }

    void feedSine (analysis::LoudnessMeter& lm, double sr, double hz, double dbfs, double seconds,
                   long long& sampleIndex, int chunk = 8192)
    {
        feedSine (lm, sr, hz, dbfs, seconds, sampleIndex, chunk, [] (double) {});
    }

    // One fresh stereo meter over a whole Tech 3341-style program of (dbfs, seconds) steps; `after` as above.
    template <typename After>
    analysis::LoudnessMeter run (double sr, const Steps& steps, int chunk, After&& after)
    {
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 200.0);
        long long idx = 0;
        for (const auto& [dbfs, seconds] : steps) feedSine (lm, sr, kToneHz, dbfs, seconds, idx, chunk, after);
        return lm;
    }

    double measure (double sr, const Steps& steps)
    {
        return run (sr, steps, 8192, [] (double) {}).integratedLufs();
    }

    std::string at (double sr, const char* what) { return std::string (what) + (sr == 48000.0 ? " @48k" : " @44.1k"); }
}

int main()
{
    std::printf ("felitronics::analysis LoudnessMeter — EBU Tech 3341 conformance\n");

    for (const double sr : { 48000.0, 44100.0 })
    {
        // --- Table 1: the integrated-loudness minimum requirements, ±0.1 LU ---
        test::group (at (sr, "Tech 3341 Table 1"));

        // Case 1: 1 kHz at −23 dBFS per channel, in phase on both channels, 20 s → M, S, I = −23.0 LUFS.
        // The arithmetic behind the round number: a sine's mean square is −3.01 dB under its peak, the
        // in-phase stereo sum adds +3.01 dB back, and the K-filter's +0.69 dB at the tone is what the −0.691
        // calibration term was chosen to cancel.
        {
            const auto lm = run (sr, { { -23.0, 20.0 } }, 8192, [] (double) {});
            test::approx (lm.integratedLufs(), -23.0, 0.1, at (sr, "case 1: I = −23.0 LUFS"));
            test::approx (lm.momentaryLufs(),  -23.0, 0.1, at (sr, "case 1: M = −23.0 LUFS"));
            test::approx (lm.shortTermLufs(),  -23.0, 0.1, at (sr, "case 1: S = −23.0 LUFS"));
        }

        // Case 2: the same tone at −33 dBFS → −33.0 (linearity), again for all three.
        {
            const auto lm = run (sr, { { -33.0, 20.0 } }, 8192, [] (double) {});
            test::approx (lm.integratedLufs(), -33.0, 0.1, at (sr, "case 2: I = −33.0 LUFS"));
            test::approx (lm.momentaryLufs(),  -33.0, 0.1, at (sr, "case 2: M = −33.0 LUFS"));
            test::approx (lm.shortTermLufs(),  -33.0, 0.1, at (sr, "case 2: S = −33.0 LUFS"));
        }

        // Case 3: quiet −36 dBFS shoulders around a −23 dBFS body. The relative gate (−10 LU under the
        // absolutely-gated mean) drops the shoulders; ungated averaging would read ≈ −24.2.
        test::approx (measure (sr, { { -36.0, 10.0 }, { -23.0, 60.0 }, { -36.0, 10.0 } }), -23.0, 0.1,
                      at (sr, "case 3: the relative gate drops the −36 dBFS shoulders"));

        // Case 4: case 3 wrapped in −72 dBFS. The spec's own robustness case: −72 sits under the absolute
        // gate AND far under the relative one, so this reads −23.0 whichever gate does the dropping — the
        // absolute gate is isolated further down, where only it can make the difference.
        test::approx (measure (sr, { { -72.0, 10.0 }, { -36.0, 10.0 }, { -23.0, 60.0 }, { -36.0, 10.0 }, { -72.0, 10.0 } }),
                      -23.0, 0.1, at (sr, "case 4: the −72 dBFS tails are dropped"));

        // Case 5: −26 / −20 / −26 dBFS at 20 / 20.1 / 20 s (the spec's 20.1 is deliberate) — everything gated
        // IN (all blocks within 10 LU of the mean), power-averaging to −23.0. A meter averaging amplitude or
        // dB instead of power reads −23.5 or worse.
        test::approx (measure (sr, { { -26.0, 20.0 }, { -20.0, 20.1 }, { -26.0, 20.0 } }), -23.0, 0.1,
                      at (sr, "case 5: three levels power-average to −23.0"));
    }

    const double sr = 48000.0;

    // --- Table 1 case 6: 5.0 channels at −28 (L, R), −24 (C), −30 (Ls, Rs) dBFS with the BS.1770 weights
    //     (1.0 front, 1.41 surround — a power weight, +1.5 dB) sum to −23.0. The one case that exercises
    //     setChannelWeight; the role→channel mapping is the host's. ---
    test::group ("Tech 3341 Table 1 case 6: 5.0 channel weights");
    {
        analysis::LoudnessMeter lm; lm.prepare (sr, 5, 30.0);
        lm.setChannelWeight (3, 1.41); lm.setChannelWeight (4, 1.41);
        const auto tone = [sr] (double dbfs) {
            std::vector<float> v ((std::size_t) std::llround (20.0 * sr));
            const double amp = std::pow (10.0, dbfs / 20.0);
            for (std::size_t i = 0; i < v.size(); ++i) v[i] = (float) (amp * std::sin (2.0 * core::kPi * kToneHz * (double) i / sr));
            return v;
        };
        const auto front = tone (-28.0), centre = tone (-24.0), surround = tone (-30.0);
        for (std::size_t o = 0; o < front.size(); o += 8192)
        {
            const int n = (int) std::min<std::size_t> (8192, front.size() - o);
            const float* ch[5] { front.data() + o, front.data() + o, centre.data() + o, surround.data() + o, surround.data() + o };
            lm.process (ch, 5, n);
        }
        test::approx (lm.integratedLufs(), -23.0, 0.1, "case 6: L/R −28, C −24, Ls/Rs −30 with weights → −23.0");
    }

    // --- Table 1 cases 9 and 12: over a program periodic in exactly one window, the windowed measure is a
    //     constant once the window is full — S over (1.34 s at −20, 1.66 s at −30) × 5 after 3 s, M over
    //     (0.18 s at −20, 0.22 s at −30) × 25 after 1 s — and that constant is the power mean, −23.0. ---
    test::group ("Tech 3341 Table 1 cases 9 and 12: S and M settle on a periodic program");
    {
        Steps nine, twelve;
        for (int r = 0; r < 5; ++r)  { nine.push_back ({ -20.0, 1.34 });  nine.push_back ({ -30.0, 1.66 }); }
        for (int r = 0; r < 25; ++r) { twelve.push_back ({ -20.0, 0.18 }); twelve.push_back ({ -30.0, 0.22 }); }
        const int hop = (int) std::lround (0.1 * sr);   // read the meter once per hop, as a display would
        double lo = 0.0, hi = -200.0;
        {
            analysis::LoudnessMeter lm; lm.prepare (sr, 2, 60.0);
            long long idx = 0;
            for (const auto& [dbfs, seconds] : nine)
                feedSine (lm, sr, kToneHz, dbfs, seconds, idx, hop, [&] (double t) { if (t >= 3.0) { lo = std::min (lo, lm.shortTermLufs()); hi = std::max (hi, lm.shortTermLufs()); } });
            test::approx (lo, -23.0, 0.1, "case 9: S after 3 s never reads under −23.0 ± 0.1");
            test::approx (hi, -23.0, 0.1, "case 9: S after 3 s never reads over −23.0 ± 0.1");
        }
        {
            analysis::LoudnessMeter lm; lm.prepare (sr, 2, 60.0); lo = 0.0; hi = -200.0;
            long long idx = 0;
            for (const auto& [dbfs, seconds] : twelve)
                feedSine (lm, sr, kToneHz, dbfs, seconds, idx, hop, [&] (double t) { if (t >= 1.0) { lo = std::min (lo, lm.momentaryLufs()); hi = std::max (hi, lm.momentaryLufs()); } });
            test::approx (lo, -23.0, 0.1, "case 12: M after 1 s never reads under −23.0 ± 0.1");
            test::approx (hi, -23.0, 0.1, "case 12: M after 1 s never reads over −23.0 ± 0.1");
        }
    }

    // --- the relative gate must be computed from ALL absolutely-gated blocks, not on the fly. Loud body FIRST
    //     is the control: a running threshold set by the loud head would drop the quiet tail too. Loud body LAST
    //     is the discriminating case: a running threshold would have admitted the whole quiet head (20 s at
    //     −36 dBFS), and the ungated power mean of that program is −24.2 LUFS — outside ±0.1. ---
    test::group ("relative gate is a two-pass over the whole program");
    {
        test::approx (measure (sr, { { -23.0, 60.0 }, { -36.0, 20.0 } }), -23.0, 0.1, "loud body first → −23.0 (the quiet tail is gated out)");
        test::approx (measure (sr, { { -36.0, 20.0 }, { -23.0, 60.0 } }), -23.0, 0.1, "loud body last → −23.0 (a running threshold would read −24.2)");
    }

    // --- what cannot be measured answers −120, never a made-up number ---
    test::group ("unmeasurable input answers the −120 sentinel");
    {
        // Digital silence: blocks exist, every one of them is under the absolute gate.
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 10.0);
        std::vector<float> zeros ((std::size_t) sr, 0.0f);
        const float* ch[2] { zeros.data(), zeros.data() };
        lm.process (ch, 2, (int) sr);
        test::ok (lm.integratedLufs() == kNone, "one second of digital silence → −120");

        // A tone wholly under the absolute gate: −80 dBFS reads −80 LUFS, below −70 — gated to nothing.
        test::ok (measure (sr, { { -80.0, 2.0 } }) == kNone, "a −80 dBFS tone (under the −70 LUFS gate) → −120");

        // 399 ms is one hop short of the first 400 ms gating block; exactly 400 ms is exactly one block.
        test::ok (measure (sr, { { -23.0, 0.399 } }) == kNone, "399 ms → no block yet → −120");
        test::ok (measure (sr, { { -23.0, 0.400 } }) > -70.0,  "400 ms → one block → a reading");

        // Any real gated mean averages blocks that each passed the −70 LUFS gate, so it lies above the gate by
        // construction: the gate is the line a consumer can draw between a reading and the sentinel.
        test::ok (measure (sr, { { -69.0, 5.0 } }) > -70.0, "a −69 dBFS tone (just over the gate) reads above −70");
    }

    // --- the K in K-weighting: shelf up high, RLB high-pass down low ---
    test::group ("K-weighting shapes the reading");
    {
        // 8 kHz sits near the top of the shelf's knee: the prototype's gain is +4 dB asymptotically, so the
        // reading lands a few dB hot of −23.
        analysis::LoudnessMeter hi; hi.prepare (sr, 2, 20.0);
        long long n1 = 0; feedSine (hi, sr, 8000.0, -23.0, 10.0, n1);
        test::ok (hi.integratedLufs() > -23.0 + 2.5 && hi.integratedLufs() < -23.0 + 4.5, "8 kHz at −23 dBFS reads +2.5..+4.5 LU hot (the shelf)");

        // 60 Hz is on the RLB high-pass slope (f0 ≈ 38 Hz): |H| ≈ −3 dB there.
        analysis::LoudnessMeter lo; lo.prepare (sr, 2, 20.0);
        long long n2 = 0; feedSine (lo, sr, 60.0, -23.0, 10.0, n2);
        test::ok (lo.integratedLufs() < -23.0 - 2.0 && lo.integratedLufs() > -23.0 - 4.5, "60 Hz at −23 dBFS reads 2..4.5 LU under (the RLB high-pass)");
    }

    // --- chunking must not matter: one call vs odd-sized pieces straddling every hop ---
    test::group ("chunk invariance");
    {
        analysis::LoudnessMeter one;  one.prepare (sr, 2, 20.0);
        analysis::LoudnessMeter many; many.prepare (sr, 2, 20.0);
        long long a = 0, b = 0;
        feedSine (one,  sr, kToneHz, -23.0, 5.0, a, 240000);   // one call
        feedSine (many, sr, kToneHz, -23.0, 5.0, b, 4097);     // odd chunks, straddling hops
        test::approx (many.integratedLufs(), one.integratedLufs(), 1e-9, "the same stream in 4097-frame pieces reads identically");
        test::approx (many.momentaryLufs(),  one.momentaryLufs(),  1e-9, "momentary too");
    }

    // --- the stated capacity is honoured to its last block: a program exactly as long as prepare() was told
    //     records every block (the +4 headroom covers the hop-count rounding), and reads as an uncapped one ---
    test::group ("maxDurationSec holds a program of exactly that length");
    {
        analysis::LoudnessMeter tight; tight.prepare (sr, 2, 3.0);
        analysis::LoudnessMeter roomy; roomy.prepare (sr, 2, 60.0);
        long long a = 0, b = 0;
        feedSine (tight, sr, kToneHz, -20.0, 3.0, a);
        feedSine (roomy, sr, kToneHz, -20.0, 3.0, b);
        test::approx (tight.integratedLufs(), roomy.integratedLufs(), 1e-9, "3 s of capacity holds 3 s of program");
    }

    return test::report();
}
