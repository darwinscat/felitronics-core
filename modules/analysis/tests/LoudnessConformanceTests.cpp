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
// Then what Table 1 alone would let a wrong meter get away with, each pinned by a signal built so that only
// the property under test can change the reading: the absolute gate isolated from the relative one and at its
// boundary; the relative gate at −10 LU; every channel counted, in phase or not; the K in K-weighting at both
// rates against the spec's published coefficients; a burst that tells 75 % block overlap from none; the −120
// answer for what cannot be measured; chunk invariance; and process() allocating nothing.
//
// Grown from the Looper Cat conformance suite that gated the product's move onto this meter.

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/analysis/LoudnessMeter.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

// RT-safety witness: every allocation in this binary is counted (the TruePeakMeter suite's pattern).
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

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

    // How a tone is laid across the channels — the fixtures that tell a per-channel power sum from every
    // cheaper thing a meter could do instead (drop a channel, sum the samples, count one channel twice).
    enum class Layout { leftOnly, rightOnly, antiPhase, mono };

    double measureLayout (double sr, double dbfs, double seconds, Layout layout)
    {
        const int channels = layout == Layout::mono ? 1 : 2;
        analysis::LoudnessMeter lm; lm.prepare (sr, channels, 60.0);
        const double amp = std::pow (10.0, dbfs / 20.0);
        const double w   = 2.0 * core::kPi * kToneHz / sr;
        const long long frames = std::llround (seconds * sr);
        std::vector<float> a (8192), b (8192, 0.0f);
        for (long long o = 0; o < frames; o += 8192)
        {
            const int n = (int) std::min<long long> (8192, frames - o);
            for (int i = 0; i < n; ++i)
            {
                const float v = (float) (amp * std::sin (w * (double) (o + i)));
                switch (layout)
                {
                    case Layout::leftOnly:  a[(std::size_t) i] = v;  b[(std::size_t) i] = 0.0f; break;
                    case Layout::rightOnly: a[(std::size_t) i] = 0.0f; b[(std::size_t) i] = v;  break;
                    case Layout::antiPhase: a[(std::size_t) i] = v;  b[(std::size_t) i] = -v;   break;
                    case Layout::mono:      a[(std::size_t) i] = v;  break;
                }
            }
            const float* ch[2] { a.data(), b.data() };
            lm.process (ch, channels, n);
        }
        return lm.integratedLufs();
    }
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

        // 399 ms is one hop short of the first 400 ms gating block; exactly 400 ms is exactly one block — and
        // that block reads the tone, not merely something above the sentinel.
        test::ok (measure (sr, { { -23.0, 0.399 } }) == kNone, "399 ms → no block yet → −120");
        test::approx (measure (sr, { { -23.0, 0.400 } }), -23.0, 0.2, "400 ms → one block → −23.0");

        // Any real gated mean averages blocks that each passed the −70 LUFS gate, so it lies above the gate by
        // construction: the gate is the line a consumer can draw between a reading and the sentinel.
        test::approx (measure (sr, { { -69.0, 5.0 } }), -69.0, 0.1, "a −69 dBFS tone (just over the gate) reads −69.0");
    }

    // --- the absolute gate, isolated: Table 1's case 4 drops its −72 dBFS tails by the relative gate as much as
    //     by the absolute one, so a meter with no absolute gate at all passes it. Here the sub-gate part sits
    //     INSIDE the relative gate — only the absolute gate can drop it. ---
    test::group ("absolute gate at −70 LUFS, isolated from the relative gate");
    {
        // −65 dBFS body, −72 dBFS tail, equal lengths: the abs-gated mean is −65.0 and the relative gate
        // −75.0, so −72 would be kept by the relative gate — and a meter keeping it reads their mean, −67.2.
        test::approx (measure (sr, { { -65.0, 10.0 }, { -72.0, 10.0 } }), -65.0, 0.1, "−65 body + −72 tail → −65.0 (the tail is dropped by the absolute gate alone)");
        // The boundary itself: half a dB over the gate is a reading, half a dB under is the sentinel.
        test::approx (measure (sr, { { -69.5, 5.0 } }), -69.5, 0.1, "−69.5 dBFS (just over the gate) reads −69.5");
        test::ok (measure (sr, { { -70.5, 5.0 } }) == kNone, "−70.5 dBFS (just under the gate) → −120");
    }

    // --- the relative gate at −10 LU, not −8 (ATSC A/85's first edition) and not −12 ---
    test::group ("relative gate at −10 LU");
    {
        // 60 s at −23 then 10 s at −32: the abs-gated mean is −23.6 and the relative gate −33.6, so the −32
        // segment sits inside it and COUNTS — the reading is the mean of both, −23.58. At −34 the segment
        // falls outside and the reading is the body alone, −23.0. A −8 LU gate would drop both (reading −23.0
        // twice); a −12 LU gate would keep both (−23.58, then −23.61).
        test::approx (measure (sr, { { -23.0, 60.0 }, { -32.0, 10.0 } }), -23.58, 0.1, "a segment 9 LU under the body is counted → −23.58");
        test::approx (measure (sr, { { -23.0, 60.0 }, { -34.0, 10.0 } }), -23.00, 0.1, "a segment 11 LU under the body is dropped → −23.0");
    }

    // --- every channel counts, as power, in phase or not. Table 1 is dual-mono throughout, which a meter that
    //     reads one channel and adds 3 dB, or sums the samples before squaring, passes just the same. ---
    test::group ("channel handling: per-channel power, summed");
    {
        // One channel alone is the mono figure, −3.01 under the dual-mono one: −26.01. Either channel.
        test::approx (measureLayout (sr, -23.0, 10.0, Layout::leftOnly),  -26.01, 0.1, "left only → −26.01");
        test::approx (measureLayout (sr, -23.0, 10.0, Layout::rightOnly), -26.01, 0.1, "right only → −26.01");
        // Anti-phase stereo carries the same power as in-phase: a meter summing samples before squaring
        // would read silence.
        test::approx (measureLayout (sr, -23.0, 10.0, Layout::antiPhase), -23.0, 0.1, "L = −R → −23.0 (power, not a sample sum)");
        // A mono meter is the single-channel figure.
        test::approx (measureLayout (sr, -23.0, 10.0, Layout::mono), -26.01, 0.1, "mono → −26.01");
    }

    // --- the K in K-weighting, at both rates. Table 1 sits at the one frequency where the K-filter's gain and
    //     the −0.691 calibration cancel, so a flat meter passes all of it; away from 1 kHz the shape shows.
    //     Expected readings come from BS.1770-4's PUBLISHED 48 kHz coefficient table, evaluated analytically
    //     at the tone: |K| = −2.899 dB at 60 Hz (the RLB high-pass), +4.039 dB at 8 kHz (the shelf); a
    //     −23 dBFS dual-mono tone therefore reads −23 + |K| − 0.691. The 44.1 kHz design differs from the
    //     48 kHz table by under 0.005 dB at these tones. ---
    for (const double rate : { 48000.0, 44100.0 })
    {
        test::group (at (rate, "K-weighting shapes the reading"));
        analysis::LoudnessMeter hi; hi.prepare (rate, 2, 20.0);
        long long n1 = 0; feedSine (hi, rate, 8000.0, -23.0, 10.0, n1);
        test::approx (hi.integratedLufs(), -19.65, 0.15, at (rate, "8 kHz at −23 dBFS reads −19.65 (the shelf's +4.04 dB)"));

        analysis::LoudnessMeter lo; lo.prepare (rate, 2, 20.0);
        long long n2 = 0; feedSine (lo, rate, 60.0, -23.0, 10.0, n2);
        test::approx (lo.integratedLufs(), -26.59, 0.15, at (rate, "60 Hz at −23 dBFS reads −26.59 (the high-pass's −2.90 dB)"));
    }

    // --- 400 ms blocks at a 100 ms hop: steady tones cannot tell 75 % overlap from none, a burst can ---
    test::group ("block overlap: a burst tells 75 % from none");
    {
        // 200 ms of −20 dBFS tone from 300 ms to 500 ms in one second of silence, hop-aligned. Five
        // overlapping 400 ms blocks touch it, holding 1/4, 1/2, 1/2, 1/2 and 1/4 of a full block's power; the
        // silent blocks fall under the absolute gate and the relative gate (−10 LU under the mean) keeps all
        // five, so the reading is the tone −10·log10(0.4) = −3.98 LU: −23.98. Non-overlapping blocks would
        // hold 1/4 and 1/4 and read −6.02 LU under: −26.02.
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 10.0);
        std::vector<float> sig ((std::size_t) sr, 0.0f);
        const double amp = std::pow (10.0, -20.0 / 20.0);
        for (std::size_t i = (std::size_t) (0.3 * sr); i < (std::size_t) (0.5 * sr); ++i)
            sig[i] = (float) (amp * std::sin (2.0 * core::kPi * kToneHz * (double) (i - (std::size_t) (0.3 * sr)) / sr));
        for (std::size_t o = 0; o < sig.size(); o += 4097)
        {
            const int n = (int) std::min<std::size_t> (4097, sig.size() - o);
            const float* ch[2] { sig.data() + o, sig.data() + o };
            lm.process (ch, 2, n);
        }
        test::approx (lm.integratedLufs(), -20.0 - 3.98, 0.3, "a 200 ms burst reads the tone −3.98 LU (five overlapping blocks)");
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

    // --- the headline RT claim: process() allocates nothing, through hops, blocks and short-term samples ---
    test::group ("process() does not allocate");
    {
        analysis::LoudnessMeter lm; lm.prepare (sr, 2, 10.0);
        std::vector<float> l (4800, 0.3f), r (4800, -0.2f);
        const float* io[2] { l.data(), r.data() };
        lm.process (io, 2, 4800);                                   // warm: the first hop
        const long before = g_allocs.load();
        for (int i = 0; i < 40; ++i) lm.process (io, 2, 4800);     // 4 s: hops, blocks and short-term samples all fire
        const bool noAlloc = (g_allocs.load() == before);
        test::okNoAlloc (noAlloc, "40 hops of process() did not allocate");
        test::ok (std::isfinite (lm.integratedLufs()), "and the meter still reads");
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
