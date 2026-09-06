// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// A channel that stops being dithered and starts again. Two things freeze: the noise-shaper error
// history, which is sample memory, and the auto-blank counter, which is a CLOCK — it records how long
// this channel has been digitally silent, and that keeps being true while nobody feeds it.
//
// Measured before this suite: the returning channel emitted 3.58e-07 into digital silence — about six
// LSB of 24 bit — where a channel that never left emitted EXACTLY 0.0. Small, but this is the stage
// P6 F2 caught amplifying 1e-15 into three LSB of the export, so its counters matter.

#include <felitronics_test.h>
#include <felitronics/dither/Dither.h>

#include <cmath>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0005f; }
}

int main()
{
    std::printf ("felitronics::dither — execution-gate state\n");

    const int N = 128;

    group ("a returning channel blanks exactly like one that never left");
    {
        dither::Dither dut, ref;
        dut.prepare (48000.0, N, 2); ref.prepare (48000.0, N, 2);
        dither::DitherParams p; dut.setParams (p); ref.setParams (p);

        std::vector<float> a (N), b (N), ra (N), rb (N);
        float* dio[2] { a.data(), b.data() };
        float* rio[2] { ra.data(), rb.data() };

        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            std::fill (a.begin(), a.end(), 0.0f); fillNoise (b, 9u + (unsigned) k);
            ra = a; rb = b;
            dut.process (dio, 2, N); ref.process (rio, 2, N);
            charged = std::fmax (charged, peakOf (b));
        }
        ok (charged > 1.0e-4, "precondition: the channel under test really carried signal");

        for (int k = 0; k < 40; ++k)                       // the DUT loses it; the reference keeps it, on silence
        {
            std::fill (a.begin(), a.end(), 0.0f);
            std::fill (ra.begin(), ra.end(), 0.0f); std::fill (rb.begin(), rb.end(), 0.0f);
            dut.process (dio, 1, N); ref.process (rio, 2, N);
        }

        double wd = 0.0, wr = 0.0;
        for (int k = 0; k < 20; ++k)
        {
            std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
            std::fill (ra.begin(), ra.end(), 0.0f); std::fill (rb.begin(), rb.end(), 0.0f);
            dut.process (dio, 2, N); ref.process (rio, 2, N);
            wd = std::fmax (wd, peakOf (b)); wr = std::fmax (wr, peakOf (rb));
        }
        ok (wr == 0.0, "precondition: the reference blanked, so 'exactly zero' is the right bar");
        ok (wd == 0.0, "the returning channel is blanked too (was 3.58e-07, ~6 LSB of 24 bit)");
    }

    group ("the counter is a clock, not a latch");
    {
        // Advancing it once on the edge credits a single block and sends the returning channel back to the
        // start of the 4096-sample window. It has to keep counting for every channel nobody is feeding, so
        // a SHORT absence must leave it still un-blanked and a long one blanked — the counter is real time.
        dither::Dither d; d.prepare (48000.0, N, 2);
        dither::DitherParams p; p.autoBlankSamples = 4096; d.setParams (p);

        std::vector<float> a (N), b (N);
        float* io[2] { a.data(), b.data() };
        for (int k = 0; k < 40; ++k) { std::fill (a.begin(), a.end(), 0.0f); fillNoise (b, 3u + (unsigned) k); d.process (io, 2, N); }

        for (int k = 0; k < 4; ++k) { std::fill (a.begin(), a.end(), 0.0f); d.process (io, 1, N); }   // 512 samples away
        std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
        d.process (io, 2, N);
        ok (peakOf (b) > 0.0, "a short absence is not yet a blank — the window is real time, not a flag");

        for (int k = 0; k < 40; ++k) { std::fill (a.begin(), a.end(), 0.0f); d.process (io, 1, N); }  // well past it
        double worst = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::fill (a.begin(), a.end(), 0.0f); std::fill (b.begin(), b.end(), 0.0f);
            d.process (io, 2, N); worst = std::fmax (worst, peakOf (b));
        }
        ok (worst == 0.0, "and a long one is");
    }

    group ("the shaper history is dropped too, and here is how that can be seen at all");
    {
        // This one needs an oracle built on purpose. The noise source is a PCG advanced ONCE PER PROCESSED
        // SAMPLE, so two runs with different processing histories can never be compared bit for bit — which
        // is why an "exact zero out of silence" check cannot see the shaper at all: by the time silence has
        // blanked the channel, the blank path has cleared the history anyway.
        //
        // So: two instances that process EXACTLY the same number of samples on channel 1 — their PCG
        // streams stay in lockstep, because a draw happens per sample regardless of its value once
        // autoBlank is off — differing only in whether that channel's shaper was CHARGED before it left.
        // Both then lose the channel and get it back. If the drop clears the history, the two must agree
        // bit for bit; if it does not, the charged one replays it. Measured without the drop: 249 of 512
        // samples differ, worst 6.10e-05, which is two LSB of 16 bit.
        const int N = 64;
        dither::Dither hot, cold;
        hot.prepare (48000.0, N, 2); cold.prepare (48000.0, N, 2);
        dither::DitherParams p;
        p.bits = 16;                                    // a big LSB, so the shaper's feedback is visible
        p.shaping = dither::NoiseShaping::Weighted;
        p.autoBlank = false;                            // the blank path would otherwise clear the history
        hot.setParams (p); cold.setParams (p);

        std::vector<float> ha (N), hb (N), ca (N), cb (N);
        float* hio[2] { ha.data(), hb.data() };
        float* cio[2] { ca.data(), cb.data() };

        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            std::fill (ha.begin(), ha.end(), 0.0f); fillNoise (hb, 17u + (unsigned) k);
            std::fill (ca.begin(), ca.end(), 0.0f); std::fill (cb.begin(), cb.end(), 0.0f);
            hot.process (hio, 2, N); cold.process (cio, 2, N);
            charged = std::fmax (charged, peakOf (hb));
        }
        ok (charged > 1.0e-4, "precondition: one instance really charged its shaper and the other did not");

        for (int k = 0; k < 4; ++k)                     // channel 1 leaves on BOTH, so the PCGs stay aligned
        {
            std::fill (ha.begin(), ha.end(), 0.0f); std::fill (ca.begin(), ca.end(), 0.0f);
            hot.process (hio, 1, N); cold.process (cio, 1, N);
        }

        int differing = 0;
        for (int k = 0; k < 8; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                const float v = (float) (0.10 + 0.001 * std::sin (0.05 * (k * N + i)));
                ha[(std::size_t) i] = ca[(std::size_t) i] = 0.0f;
                hb[(std::size_t) i] = cb[(std::size_t) i] = v;
            }
            hot.process (hio, 2, N); cold.process (cio, 2, N);
            for (int i = 0; i < N; ++i) if (! (hb[(std::size_t) i] == cb[(std::size_t) i])) ++differing;
        }
        ok (differing == 0, "the returning channel carries no error history from before it left");
    }

    return felitronics::test::report();
}
