// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of OrbitCapture NAM — see LICENSE.

// Headless test for MagnitudeCurve.h — a magnitude table turned into a filter. The gate that
// matters is the one the old 1-pole fit could not pass: a Big Muff-shaped tone control (bass cut,
// treble boost, ~30 dB of tilt) has to come back out of the FIR within a fraction of a dB.

#include <felitronics_test.h>
#include <felitronics/lineareq/MagnitudeCurve.h>

#include <cmath>
#include <cstdio>
#include <vector>

using felitronics::test::ok;
using felitronics::test::group;
using namespace felitronics::lineareq;

namespace {

std::vector<double> grid(std::size_t n = 96) {
    std::vector<double> f(n);
    for (std::size_t i = 0; i < n; ++i) f[i] = 20.0 * std::pow(1000.0, (double) i / (double) (n - 1));
    return f;
}

// The FIR's magnitude response at `hz`, by direct evaluation of its DTFT — no FFT, so the test is
// checking the taps themselves rather than agreeing with the same maths that built them.
double firDbAt(const std::vector<float>& fir, double hz, double sr) {
    double re = 0.0, im = 0.0;
    const double w = 2.0 * 3.14159265358979 * hz / sr;
    for (std::size_t i = 0; i < fir.size(); ++i) {
        re += fir[i] * std::cos(w * (double) i);
        im -= fir[i] * std::sin(w * (double) i);
    }
    return 10.0 * std::log10(std::max(re * re + im * im, 1e-30));
}

} // namespace

int main() {
    std::printf("MagnitudeCurve tests\n");
    const auto f = grid();
    const double sr = 48000.0;

    group("a tilt no 1-pole could fit");
    {
        // The Big Muff tone control as measured: −14 dB at 100 Hz climbing to +16 dB at 6.4 kHz.
        std::vector<double> db(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
            const double t = std::log(f[i] / 100.0) / std::log(6400.0 / 100.0);
            db[i] = -14.0 + 30.0 * std::clamp(t, 0.0, 1.0);
        }
        const auto fir = magnitudeCurveToFir(db, f, sr);
        ok(! fir.empty(), "a curve that does something builds a filter");
        double worst = 0.0;
        for (double hz : { 120.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 6000.0 })
            worst = std::max(worst, std::abs(firDbAt(fir, hz, sr) - curveDbAt(db, f, hz)));
        std::printf("    worst error across the band: %.2f dB\n", worst);
        ok(worst < 1.0, "the FIR reproduces the measured tilt to within a dB");
    }

    // Outside the band where the measurement repeated, the curve is held at its edge value — the pack
    // ships the whole thing, the player declines to believe the part that did not reproduce.
    group("held outside the trusted band");
    {
        std::vector<double> db(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) db[i] = f[i] > 9000.0 ? 24.0 : 3.0;   // junk up top
        const auto held = heldOutsideBand(db, f, 80.0, 9000.0);
        ok(held.back() == 3.0, "the untrusted top octave is held at the last value that was trusted");
        ok(held.front() == 3.0 && held[f.size() / 2] == 3.0, "and the trusted part is untouched");

        const auto wild = magnitudeCurveToFir(db, f, sr);
        const auto sane = magnitudeCurveToFir(db, f, sr, 1024, 8192, 80.0, 9000.0);
        ok(firDbAt(wild, 16000.0, sr) > firDbAt(sane, 16000.0, sr) + 10.0,
           "…so the filter stops applying 24 dB of something nobody measured twice");
    }

    // The strongest warning the analysis can produce must not read as the widest permission. An empty
    // band means the curve reproduced at NO frequency; a guard that only recognised a proper band
    // returned the curve untouched, so "tested and failed everywhere" and "no band stated" behaved
    // identically — and the producer could emit exactly that, because the band's growth seeded at
    // 1 kHz without ever testing 1 kHz.
    group("an empty trusted band means apply nothing");
    {
        std::vector<double> db(f.size(), 6.0);
        const auto none = heldOutsideBand(db, f, 1000.0, 1000.0);
        ok(! none.empty() && none.front() == 0.0 && none.back() == 0.0,
           "a band with no width flattens the curve instead of licensing all of it");
        const auto untested = heldOutsideBand(db, f, 0.0, 0.0);
        ok(untested == db, "…while NO band stated is a different fact, and leaves the curve alone");
        ok(magnitudeCurveToFir(db, f, sr, 1024, 8192, 1000.0, 1000.0).empty(),
           "and nothing is built from it, so the control is silent rather than wrong");
    }

    group("a flat curve is not worth a convolution");
    {
        const std::vector<double> flat(f.size(), 0.0);
        ok(magnitudeCurveToFir(flat, f, sr).empty(),
           "the reference position is flat by construction — it must build nothing");
        ok(magnitudeCurveToFir({}, {}, sr).empty(), "and no curve builds nothing");
    }

    group("reading a curve between its grid points");
    {
        const std::vector<double> freq { 100.0, 1000.0, 10000.0 };
        const std::vector<double> db   {   0.0,  -12.0,     0.0 };
        ok(std::abs(curveDbAt(db, freq, 1000.0) + 12.0) < 1e-9, "an exact grid point reads back exactly");
        // Interpolated in LOG frequency: halfway from 100 to 1000 in octaves is ~316 Hz, not 550.
        ok(std::abs(curveDbAt(db, freq, 316.23) + 6.0) < 0.05, "the middle in log-f is the middle in dB");
        ok(curveDbAt(db, freq, 20.0) == 0.0 && curveDbAt(db, freq, 20000.0) == 0.0,
           "past either end the curve is held, never extrapolated into invention");
    }

    // A measured knob is CONTINUOUS — that is the whole reason for measuring instead of capturing. A
    // player has to be able to stand between two swept positions.
    group("between two swept positions");
    {
        const std::vector<std::vector<double>> positions {
            { 0.0, 0.0, 0.0 }, { -6.0, -12.0, -18.0 }
        };
        const std::vector<double> evenNorms { 0.0, 1.0 };
        const auto mid = curveAtPosition(positions, evenNorms, 0.5);
        ok(mid.size() == 3 && std::abs(mid[1] + 6.0) < 1e-9, "halfway between is halfway in dB");
        ok(curveAtPosition(positions, evenNorms, 0.0) == positions[0], "a whole position is itself");
        ok(curveAtPosition(positions, evenNorms, 7.0) == positions[1],
           "past the end HOLDS the last position - never extrapolated, which would invent decibels");

        // The axis matters. Sweep denser where the pot acts fastest and the swept positions are no
        // longer evenly spaced: at rotation 0.5 the pack's own norms land exactly on the middle curve,
        // while treating the array as evenly spaced would return a blend of its neighbours.
        const std::vector<std::vector<double>> three { { 0.0 }, { -10.0 }, { -30.0 } };
        const std::vector<double> uneven { 0.0, 0.5, 1.0 }, skewed { 0.0, 0.1, 1.0 };
        ok(std::abs(curveAtPosition(three, uneven, 0.5)[0] + 10.0) < 1e-9,
           "evenly swept: rotation 0.5 is the middle curve");
        ok(std::abs(curveAtPosition(three, skewed, 0.5)[0] + 18.9) < 0.05,
           "unevenly swept: the SAME rotation lands 8.9 dB away, because the pack says where it lands");
    }

    return felitronics::test::report();
}
