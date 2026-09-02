// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// EBU Tech 3341 conformance for felitronics::analysis::LoudnessMeter — the integrated measure against the
// THEORY of BS.1770-4, never against the implementation. Tech 3341 Table 1 defines its test signals
// mathematically (997 Hz sines at set levels and durations, expected integrated loudness ±0.1 LU), so every
// fixture here is synthesized from the spec: the calibration (case 1), linearity (case 2), the −10 LU relative
// gate (case 3), the −70 LUFS absolute gate (case 4) and power-averaging across levels (case 5) — at both base
// rates the family runs at. Then the parts a flat RMS meter would fail: the K in K-weighting (a high tone reads
// hot, a low one reads under), the −120 answer for what cannot be measured (silence, sub-gate signal, less than
// one 400 ms block), and chunk invariance (the same stream in odd-sized pieces reads identically).
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
    constexpr double kToneHz = 997.0;      // Tech 3341's reference tone
    constexpr double kNone   = -120.0;     // the meter's "nothing survived the gates" answer

    // Appends `seconds` of a dual-mono sine at `dbfs` peak amplitude, phase-continuous across calls through
    // the caller's running sample counter (the level steps of cases 3–5 must add no click energy), fed in
    // `chunk`-frame slices so the meter sees the block boundaries a real stream produces.
    void feedSine (analysis::LoudnessMeter& lm, double sr, double hz, double dbfs, double seconds,
                   long long& sampleIndex, int chunk = 8192)
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
        }
    }

    // One fresh stereo meter over a whole Tech 3341-style program of (dbfs, seconds) steps.
    double measure (double sr, const std::vector<std::pair<double, double>>& steps)
    {
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 200.0);
        long long idx = 0;
        for (const auto& [dbfs, seconds] : steps) feedSine (lm, sr, kToneHz, dbfs, seconds, idx);
        return lm.integratedLufs();
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

        // Case 1: 997 Hz at −23 dBFS, 20 s → −23.0 LUFS. The −0.691 calibration term exactly cancels the
        // K-filter's gain at 997 Hz — dual-mono stereo sums to the stated figure.
        test::approx (measure (sr, { { -23.0, 20.0 } }), -23.0, 0.1, at (sr, "case 1: −23 dBFS → −23.0 LUFS"));

        // Case 2: the same tone at −33 dBFS → −33.0 LUFS (linearity).
        test::approx (measure (sr, { { -33.0, 20.0 } }), -33.0, 0.1, at (sr, "case 2: −33 dBFS → −33.0 LUFS"));

        // Case 3: quiet −36 dBFS shoulders around a −23 dBFS body. The relative gate (−10 LU under the
        // absolutely-gated mean) drops the shoulders; ungated averaging would read ≈ −24.2.
        test::approx (measure (sr, { { -36.0, 10.0 }, { -23.0, 60.0 }, { -36.0, 10.0 } }), -23.0, 0.1,
                      at (sr, "case 3: the relative gate drops the −36 dBFS shoulders"));

        // Case 4: case 3 wrapped in −72 dBFS — under the −70 LUFS absolute gate, so those blocks never even
        // join the mean the relative gate is computed from.
        test::approx (measure (sr, { { -72.0, 10.0 }, { -36.0, 10.0 }, { -23.0, 60.0 }, { -36.0, 10.0 }, { -72.0, 10.0 } }),
                      -23.0, 0.1, at (sr, "case 4: the absolute gate drops the −72 dBFS tails"));

        // Case 5: −26 / −20 / −26 dBFS at 20 / 20.1 / 20 s — everything gated IN (all blocks within 10 LU of
        // the mean), power-averaging to −23.0.
        test::approx (measure (sr, { { -26.0, 20.0 }, { -20.0, 20.1 }, { -26.0, 20.0 } }), -23.0, 0.1,
                      at (sr, "case 5: three levels power-average to −23.0"));
    }

    const double sr = 48000.0;

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

        // 60 Hz is on the RLB high-pass slope (f0 ≈ 38 Hz, overdamped): |H| ≈ −3 dB there.
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
