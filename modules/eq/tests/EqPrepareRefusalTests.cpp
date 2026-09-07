// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// prepare() refuses a sample rate it cannot honour, and the flush visits only the channels it prepared.
//
// The first is the second half of P6 F12: `Saturator::prepare` was closed then with "spelled positively
// so NaN fails", and the EQ, which validated NOTHING, was left. Only the mastering chain was protected,
// because it validates at its own входе; a direct consumer of `eq` got NaN in silence. Measured on the
// old code: prepare(0) and prepare(NaN) each put 63 of 64 output samples non-finite on a 0.25 input,
// and prepare(1.0) does the same to a HighPass or a Notch.
//
// The second is a cost fix with a contract behind it: the flush ran over kMaxChannels regardless of what
// was prepared, and since #126 its body is two isfinite tests and two stores per column rather than a
// flush — so a mono filter paid for sixteen, once per block, in the hottest primitive here.

#include <felitronics_test.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/eq/Svf.h>

#include <cmath>
#include <limits>
#include <vector>

using namespace felitronics;
using namespace felitronics::eq;
using felitronics::test::ok;
using felitronics::test::group;

static BandParams bell (double freq = 1000.0, double gainDb = 6.0, FilterType t = FilterType::Bell)
{
    BandParams p;
    p.on = true; p.type = t;
    LaneParams& st = p.lane (Lane::Stereo);
    st.on = true; st.freq = freq; st.Q = 1.0; st.gainDb = gainDb; st.slope = 24;
    return p;
}

static int nonFinite (const std::vector<float>& v)
{
    int n = 0; for (float x : v) if (! std::isfinite (x)) ++n; return n;
}

int main()
{
    std::printf ("felitronics::eq — prepare() refuses, and the flush is bounded\n");

    const int N = 64;
    const double bad[] { 0.0, -48000.0, std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(), 1.0e9 };
    const char* badName[] { "0", "negative", "NaN", "infinity", "1 GHz" };

    group ("EqEngine::prepare refuses a rate it cannot honour, and then does nothing");
    for (int i = 0; i < 5; ++i)
    {
        EqEngine e;
        ok (! e.prepare (bad[i], N, 2), std::string ("prepare(") + badName[i] + ") is refused");
        e.setBand (0, bell());
        std::vector<float> L ((std::size_t) N, 0.25f), R ((std::size_t) N, 0.25f);
        float* ch[2] { L.data(), R.data() };
        ok (! e.process (ch, 2, N), "...and the refusal is RETURNED, not silent (law 11)");
        ok (nonFinite (L) == 0 && nonFinite (R) == 0,
            std::string ("...and a refused engine emits nothing non-finite (") + badName[i] + ")");
        ok (L[0] == 0.25f && R[0] == 0.25f, "...and does not touch the buffer at all");
    }

    group ("EqBand::prepare does the same, because it is where the coefficients are built");
    for (int i = 0; i < 5; ++i)
    {
        // The band is the public entry LaneDynamics and the dynamic-EQ layer drive directly, so the
        // engine's guard is not the only one that has to hold.
        EqBand b;
        ok (! b.prepare (bad[i], 2), std::string ("EqBand::prepare(") + badName[i] + ") is refused");
        b.setParams (bell (1000.0, 6.0, FilterType::HighPass));
        std::vector<float> L ((std::size_t) N, 0.25f), R ((std::size_t) N, 0.25f);
        float* ch[2] { L.data(), R.data() };
        ok (! b.processBlock (ch, 2, N), "...and the refusal is RETURNED, not silent (law 11)");
        ok (nonFinite (L) == 0 && L[0] == 0.25f, std::string ("...and stays inert (") + badName[i] + ")");
    }

    group ("a rate it CAN honour is still accepted, and still filters");
    {
        EqEngine e;
        ok (e.prepare (48000.0, N, 2), "48 kHz is accepted");
        e.setBand (0, bell (1000.0, 12.0));
        std::vector<float> L ((std::size_t) N), R ((std::size_t) N);
        double peak = 0.0;
        for (int k = 0; k < 200; ++k)
        {
            for (int i = 0; i < N; ++i)
                L[(std::size_t) i] = R[(std::size_t) i] = (float) (0.25 * std::sin (2.0 * core::kPi * 1000.0 * (k * N + i) / 48000.0));
            float* ch[2] { L.data(), R.data() };
            felitronics::test::run (e.process (ch, 2, N));
            for (float x : L) peak = std::fmax (peak, (double) std::fabs (x));
        }
        ok (peak > 0.5, "a +12 dB bell at 1 kHz really lifts a 0.25 tone — the guard did not disable the EQ");
    }

    group ("the accepted domain is exactly the one the design can express");
    {
        // The decisive test here is a PROPERTY OVER THE WHOLE ACCEPTED DOMAIN, not a handful of bad values,
        // because this failure does not announce itself. Every band clamps its frequency to [10 Hz,
        // 0.49*fs], and `std::clamp` with lo above hi is a violated precondition: measured on the previous
        // guard, fs = 20 produced a perfectly plausible 0.377 and fs = 1 produced 7.2e+28. "The output is
        // finite" is therefore not a test for it, and neither is UBSan — a library precondition is not
        // language UB. So: for EVERY rate prepare() accepts, 0.49*fs must reach 10 Hz.
        int accepted = 0, refused = 0;
        bool domainHolds = true, agree = true, finiteOut = true;
        for (int i = 0; i <= 4000; ++i)
        {
            // dense across 20.408163…, the exact point where the domain becomes non-empty, then decades out
            const double fs = (i <= 2000) ? (19.0 + 0.001 * i)
                                          : std::pow (10.0, -1.0 + 8.0 * (double) (i - 2000) / 2000.0);
            EqEngine e; EqBand b;
            const bool okE = e.prepare (fs, 64, 2);
            const bool okB = b.prepare (fs, 2);
            agree = agree && (okE == okB);
            if (okE) { ++accepted; domainHolds = domainHolds && (0.49 * fs >= 10.0); }
            else       ++refused;
            if (okE)
            {
                e.setBand (0, bell (std::min (1000.0, 0.4 * fs), 6.0));
                std::vector<float> L ((std::size_t) 32, 0.25f), R ((std::size_t) 32, 0.25f);
                float* ch[2] { L.data(), R.data() };
                felitronics::test::run (e.process (ch, 2, 32));
                finiteOut = finiteOut && (nonFinite (L) == 0);
            }
        }
        ok (accepted > 100 && refused > 100, "precondition: the sweep straddles the boundary in both directions");
        ok (domainHolds, "every accepted rate leaves the design domain [10 Hz, 0.49*fs] non-empty");
        ok (agree, "the engine and the band accept and refuse exactly the same rates");
        ok (finiteOut, "and every accepted rate produces finite audio");
    }

    group ("Svf::flushDenormals visits the prepared channels only, and that changes no output");
    {
        // The columns past `ch` are never written, so flushing them was always a no-op. Assert exactly
        // that, which is the whole claim: a column outside the prepared count behaves identically whether
        // the flush visited it or not. The two-build differential against origin/main covers the rest.
        Svf used, fresh;
        used.prepare (48000.0, 1); fresh.prepare (48000.0, 1);
        used.setParams (FilterType::Bell, 900.0, 2.0, 9.0);
        fresh.setParams (FilterType::Bell, 900.0, 2.0, 9.0);
        for (int n = 0; n < 500; ++n)
        {
            (void) used.processSample (0, (float) std::sin (0.07 * n));   // only column 0 is ever driven
            used.flushDenormals();
        }
        bool same = true;
        for (int n = 0; n < 64; ++n)
        {
            const float x = (float) std::cos (0.02 * n);
            same = same && (used.processSample (1, x) == fresh.processSample (1, x));
        }
        ok (same, "a column outside the prepared count is untouched either way — flushing it was a no-op");

        Svf wide; wide.prepare (48000.0, 4);
        wide.setParams (FilterType::Bell, 900.0, 2.0, 9.0);
        bool finite = true;
        for (int n = 0; n < 500; ++n)
        {
            for (int c = 0; c < 4; ++c) finite = finite && std::isfinite (wide.processSample (c, (float) std::sin (0.05 * n + c)));
            wide.flushDenormals();
        }
        ok (finite, "and every PREPARED column is still flushed, at four channels as at one");
    }

    return felitronics::test::report();
}
