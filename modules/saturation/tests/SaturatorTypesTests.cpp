// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The Tube and Transistor curves (WaveShaper.h). What each group pins:
//
//   * shapeAt<S>(coeffs(), x) IS processSample(x), bit for bit — the parameter glide runs the first, the settled
//     stage the second, and a Shape that fell through to another formula in shapeAt would only sound in a glide.
//   * slopeAtZero() is the curve's own derivative at 0 (a central difference, several steps, driveDb 0..36).
//   * The peak normalisation: y(0) == 0, max |y| over [-1, 1] is 1 to an ulp, and Tube's peak is the negative side.
//   * Tube's precision: within 1e-6 of a double evaluation of tanh(u+0.3)-tanh(0.3) at every drive, and a 1-ulp
//     change of tanh's result moves it by less than 1e-5 (the proxy for native/wasm parity: libm differs).
//   * Transistor's clamp: s(u) <= 1 over every float in [0, 64] (a separate ctest entry — it is exhaustive), and a
//     huge input at a huge drive keeps its sign instead of collapsing to 0.
//   * Through the whole Saturator: the harmonic signature (Tube's H2 over H3, Transistor's H3 and no even
//     harmonics), Tube's DC removed by the blocker, and Tube's autoComp overshoot.
//
// Run with the argument `exhaustive` for the full float sweep of Transistor's curve; without it, the rest.

#include <felitronics_test.h>
#include <felitronics/saturation/Saturator.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using Sat = saturation::Saturator;
using WS = saturation::WaveShaper;
using Shape = WS::Shape;

namespace
{
constexpr double kFs = 48000.0;

const char* name (Shape s)
{
    switch (s)
    {
        case Shape::Tanh:       return "Tanh";
        case Shape::Atan:       return "Atan";
        case Shape::Cubic:      return "Cubic";
        case Shape::Asym:       return "Asym";
        case Shape::Tube:       return "Tube";
        case Shape::Transistor: return "Transistor";
    }
    return "?";
}

std::string num (double v, const char* fmt = "%.4g") { char b[64]; std::snprintf (b, sizeof b, fmt, v); return b; }

// ok() that also prints what it measured — the report cites these lines.
void okp (bool cond, const std::string& msg) { std::printf ("       %s\n", msg.c_str()); ok (cond, msg); }

// A running maximum a NaN cannot hide from: std::max (w, NaN) keeps w; here a NaN wins and STAYS (a later finite
// value would pass `! (v <= NaN)` and overwrite it), and every gate below compares with <=, which a NaN fails.
template <class T> void worse (T& w, T v) { if (! std::isnan (w) && ! (v <= w)) w = v; }

bool bitEq (float a, float b) { return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b); }

// The drive the Saturator designs from a driveDb (Saturator::design()).
float kOf (float driveDb) { return (float) (core::dbToGain ((double) driveDb) - 1.0); }

WS shaper (Shape s, float k, float bias = 0.0f)
{
    WS w; w.setShape (s); w.setBias (bias); w.setDrive (k);
    return w;
}

// The raw curve itself, at drive 1 and norm 1 — shapeAt() multiplies the drive and the norm by exactly 1.
float raw (Shape s, float u)
{
    const WS::Coeffs c { 1.0f, 0.0f, 0.0f, 1.0f };
    return s == Shape::Tube ? WS::shapeAt<Shape::Tube> (c, u) : WS::shapeAt<Shape::Transistor> (c, u);
}

const float kDrivesDb[] = { 0.0f, 0.1f, 0.5f, 1.0f, 2.0f, 3.0f, 4.5f, 6.0f, 9.0f, 12.0f, 18.0f, 24.0f, 30.0f, 36.0f };

//==============================================================================
// 2) shapeAt<S>(coeffs(), x) == processSample(x), bitwise.
template <Shape S>
void shapeAtIsProcessSample()
{
    long long n = 0, bad = 0;
    for (float db : { 0.0f, 0.1f, 1.0f, 3.0f, 6.0f, 12.0f, 24.0f, 36.0f, 80.0f })
        for (float bias : { 0.0f, 0.4f })                       // Tube and Transistor ignore it — shapeAt must too
        {
            const WS w = shaper (S, kOf (db), bias);
            const WS::Coeffs c = w.coeffs();
            auto one = [&] (float x) { ++n; if (! bitEq (WS::shapeAt<S> (c, x), w.processSample (x))) ++bad; };
            for (int i = -40000; i <= 40000; ++i) one ((float) i / 20000.0f);
            for (float x : { 0.0f, -0.0f, 1.0e-30f, -1.0e-30f, 1.0e-4f, 1.0e6f, -1.0e6f, 64.0f / c.drive, -64.0f / c.drive })
                one (x);
        }
    okp (bad == 0, std::string (name (S)) + ": shapeAt<S>(coeffs(), x) == processSample(x) bit for bit over "
                  + std::to_string (n) + " points (" + std::to_string (bad) + " differ)");
}

//==============================================================================
// 3) slopeAtZero against a central difference.
void slopeAtZero (Shape s)
{
    double worst = 0.0; bool allOk = true;
    for (float db : kDrivesDb)
    {
        const WS w = shaper (s, kOf (db));
        const double slope = (double) w.slopeAtZero();
        for (double hu : { 1.0e-2, 1.0e-3, 1.0e-4 })            // the step in u = kx
        {
            const float h = (float) (hu / (double) w.drive());
            const double fd = ((double) w.processSample (h) - (double) w.processSample (-h)) / (2.0 * (double) h);
            const double rel = std::fabs (fd / slope - 1.0);
            // truncation: the u³ term, |r'''(0)/(6 r'(0))| <= 0.25 for Tube, 0 for Transistor (no cubic term);
            // rounding: a few float ulps of each output.
            const double tol = hu * hu + 2.0e-6;
            worse (worst, rel);
            if (! (rel <= tol))
            {
                allOk = false;
                ok (false, std::string (name (s)) + " driveDb " + num (db) + " step " + num (hu) + ": slopeAtZero "
                           + num (slope, "%.9g") + " vs central difference " + num (fd, "%.9g"));
            }
        }
    }
    if (s == Shape::Tube)
    {
        float minDb = 0.0f, crossDb = 0.0f; double minSlope = 2.0;
        for (int i = 1; i <= 1200; ++i)
        {
            const float db = 0.01f * (float) i;
            const double sl = (double) shaper (s, kOf (db)).slopeAtZero();
            if (sl < minSlope) { minSlope = sl; minDb = db; }
            if (sl < 1.0) crossDb = db;
        }
        okp (minSlope < 1.0 && minSlope > 0.9, "Tube: normalised slope at 0 is below 1 for driveDb in (0, " + num (crossDb)
                                                + "], minimum " + num (minSlope, "%.4f") + " at driveDb " + num (minDb));
    }
    okp (allOk, std::string (name (s)) + ": slopeAtZero == central difference to h_u² + 2e-6 at steps 1e-2..1e-4, driveDb 0..36 "
               "(worst relative " + num (worst) + ")");
}

//==============================================================================
// 4) The peak normalisation.
void peakNormalisation (Shape s)
{
    const float up = 1.0f + 0x1p-23f, down = 1.0f - 0x1p-23f;    // one ulp of 1 upward, two downward
    std::vector<float> xs;
    for (int i = -65536; i <= 65536; ++i) xs.push_back ((float) i / 65536.0f);
    for (float x = 1.0f; xs.size() < 131073 + 4096; x = std::nextafter (x, 0.0f)) { xs.push_back (x); xs.push_back (-x); }
    bool zero = true, peak = true, negSide = true;
    double hiAll = 0.0, loAll = 2.0;
    for (float db : kDrivesDb)
    {
        const WS w = shaper (s, kOf (db));
        zero = zero && w.processSample (0.0f) == 0.0f && w.processSample (-0.0f) == 0.0f;
        float pos = 0.0f, neg = 0.0f;
        for (float x : xs)
        {
            const float y = w.processSample (x);
            if (x > 0.0f) worse (pos, y);
            if (x < 0.0f) worse (neg, -y);
        }
        float m = pos; worse (m, neg);
        worse (hiAll, (double) m); if (! ((double) m >= loAll)) loAll = (double) m;
        peak = peak && std::isfinite (pos) && std::isfinite (neg) && m <= up && m >= down;
        if (s == Shape::Tube) negSide = negSide && neg > pos && pos < 1.0f;
    }
    ok (zero, std::string (name (s)) + ": y(0) == 0 exactly at every drive");
    okp (peak, std::string (name (s)) + ": max |y| over [-1, 1] is 1 within [1 - 2^-23, 1 + 2^-23] at every drive (seen ["
              + num (loAll, "%.9g") + ", " + num (hiAll, "%.9g") + "])");
    if (s == Shape::Tube)
    {
        ok (negSide, "Tube: the negative side is the peak at every drive of the sweep");
        const WS w = shaper (s, kOf (6.0f));
        const float a = (float) std::pow (10.0, -3.0 / 20.0);
        const double ratio = 20.0 * std::log10 ((double) w.processSample (a) / -(double) w.processSample (-a));
        okp (ratio < 0.0, "Tube, driveDb 6, the curve at +/-(-3 dBFS): the positive peak sits " + num (ratio, "%.2f")
                          + " dB against the negative");
    }
}

//==============================================================================
// 5) Tube's precision.
void tubePrecision()
{
    const double c = std::tanh (0.3), q = 1.0 - c * c;
    auto ulpUp = [] (float v) { return (double) std::nextafter (v, 2.0f) - (double) v; };
    okp (std::fabs ((double) WS::kTubeC - c) <= ulpUp (WS::kTubeC) && std::fabs ((double) WS::kTubeQ - q) <= ulpUp (WS::kTubeQ),
        "Tube's literals: c == tanh(0.3) and q == 1 - c² in double, each within one float ulp (c off by "
        + num (std::fabs ((double) WS::kTubeC - c)) + ", q by " + num (std::fabs ((double) WS::kTubeQ - q)) + ")");

    for (float db : { 0.0f, 0.1f, 1.0f, 6.0f, 12.0f, 36.0f })
    {
        const WS w = shaper (Shape::Tube, kOf (db));
        const double k = (double) w.drive();
        auto rd = [&] (double x) { return std::tanh (k * x + 0.3) - std::tanh (0.3); };   // the curve, another way
        const double norm = 1.0 / std::fabs (rd (-1.0));
        // What the header does NOT do, printed for the record: the difference itself, in float.
        const float kf = w.drive();
        auto diff = [&] (float x) { return std::tanh (kf * x + 0.3f) - std::tanh (0.3f); };
        const float diffNorm = 1.0f / std::fabs (diff (-1.0f));
        double worst = 0.0, worstDiff = 0.0;
        for (int i = -100000; i <= 100000; ++i)
        {
            const float x = (float) i / 100000.0f;
            worse (worst, std::fabs ((double) w.processSample (x) - rd ((double) x) * norm));
            worse (worstDiff, std::fabs ((double) (diff (x) * diffNorm) - rd ((double) x) * norm));
        }
        okp (worst <= 1.0e-6, "Tube driveDb " + num (db) + ": float within 1e-6 of tanh(kx+0.3)-tanh(0.3) in double over [-1, 1] ("
                             + num (worst) + "; the difference in float would be off by " + num (worstDiff)
                             + (db == 0.0f ? ", and gives " + num ((double) diff (1.0e-4f)) + " at x = 1e-4" : std::string()) + ")");
    }

    // The native/wasm proxy. The replica below is Tube's arithmetic with tanh's result handed in; it is first
    // pinned to the header bit for bit, then fed tanh's neighbours.
    auto rep = [] (float t) { const float ct = WS::kTubeC * t; const float qt = WS::kTubeQ * t; return qt / (1.0f + ct); };
    long long pinned = 0, n = 0;
    double worst = 0.0;
    for (float db : kDrivesDb)
    {
        const WS w = shaper (Shape::Tube, kOf (db));
        const float k = w.drive();
        for (int i = -2000; i <= 2000; ++i)
        {
            const float x = (float) i / 2000.0f, t = std::tanh (k * x), tn = std::tanh (-k);
            ++n; if (! bitEq (rep (t) * w.coeffs().norm, w.processSample (x))) ++pinned;
            const float y = rep (t) / std::fabs (rep (tn));
            for (float t1 : { std::nextafter (t, 2.0f), std::nextafter (t, -2.0f) })
                for (float tn1 : { tn, std::nextafter (tn, 2.0f), std::nextafter (tn, -2.0f) })
                    worse (worst, (double) std::fabs (rep (t1) / std::fabs (rep (tn1)) - y));
        }
    }
    ok (pinned == 0, "PRECONDITION: the replica is Tube's arithmetic, bit for bit (" + std::to_string (pinned) + "/"
                     + std::to_string (n) + " differ)");
    okp (worst < 1.0e-5, "Tube: a 1-ulp change of tanh's result moves the output by < 1e-5 at every drive (worst " + num (worst) + ")");
}

//==============================================================================
// 6) Transistor: the clamp, the sign, the symmetry.
void transistorExhaustive()
{
    const std::uint32_t top = std::bit_cast<std::uint32_t> (64.0f);
    long long over = 0, asym = 0, reversals = 0;
    std::uint32_t worstUlps = 0; float worstAt = 0.0f;
    float prev = 0.0f;
    // The curve is not monotone in float (reported below), so the peak over x in [0, 1] at drive k is the largest
    // s(u) for u <= k, not s(k): that prefix maximum, times the norm, bounds every output of the stage's curve.
    std::vector<float> ks;
    for (float db : kDrivesDb) if (kOf (db) < 64.0f) ks.push_back (std::max (kOf (db), 1.0e-4f));
    std::sort (ks.begin(), ks.end());
    std::size_t next = 0; std::uint32_t worstPeak = 0; float worstPeakK = 0.0f;
    auto peakAt = [&] (float k, float m)
    {
        const float y = m * shaper (Shape::Transistor, k).coeffs().norm;
        const std::uint32_t d = y > 1.0f ? std::bit_cast<std::uint32_t> (y) - std::bit_cast<std::uint32_t> (1.0f) : 0u;
        if (d >= worstPeak) { worstPeak = d; worstPeakK = k; }
    };
    for (std::uint32_t b = 0; b <= top; ++b)
    {
        const float u = std::bit_cast<float> (b), s = raw (Shape::Transistor, u);
        while (next < ks.size() && u > ks[next]) peakAt (ks[next++], prev);
        if (! (s <= 1.0f)) ++over;
        if (! bitEq (raw (Shape::Transistor, -u), -s)) ++asym;
        if (s < prev)
        {
            ++reversals;
            const std::uint32_t d = std::bit_cast<std::uint32_t> (prev) - std::bit_cast<std::uint32_t> (s);
            if (d > worstUlps) { worstUlps = d; worstAt = u; }
        }
        prev = std::max (prev, s);
    }
    while (next < ks.size()) peakAt (ks[next++], prev);
    okp (worstPeak <= 1u, "the normalised peak over EVERY x in [0, 1] is at most 1 + 2^-23 at every drive of the sweep below the "
                          "clamp (worst " + std::to_string (worstPeak) + " ulp above 1, at k = " + num (worstPeakK) + ")");
    std::printf ("       Transistor s(u) over %u floats in [0, 64]: %lld reversals of monotonicity, worst %u ulp (at u = %.9g)\n",
                 top + 1u, reversals, worstUlps, (double) worstAt);
    ok (over == 0, "s(u) <= 1 exactly for every float in [0, 64] (" + std::to_string (over) + " above)");
    ok (raw (Shape::Transistor, 64.0f) == 1.0f && raw (Shape::Transistor, 1.0e30f) == 1.0f, "s(64) == 1 and the clamp holds it there");
    ok (asym == 0, "s(-u) == -s(u) bit for bit for every float in [0, 64]");
}

void transistorHuge()
{
    bool good = true; std::string seen;
    for (float db : { 36.0f, 80.0f })
    {
        const WS w = shaper (Shape::Transistor, kOf (db));
        const float yp = w.processSample (1.0e6f), yn = w.processSample (-1.0e6f);
        good = good && yp > 0.0f && yn < 0.0f;
        seen += " driveDb " + num (db) + ": " + num (yp) + " / " + num (yn) + ";";
    }
    okp (good, "Transistor: x = +/-1e6 keeps its sign and is never 0 at driveDb 36 and 80 (" + seen + ")");

    long long bad = 0;
    for (float db : kDrivesDb)
    {
        const WS w = shaper (Shape::Transistor, kOf (db));
        for (int i = 0; i <= 50000; ++i)
        {
            const float x = (float) i / 25000.0f;
            if (! (w.processSample (-x) == -w.processSample (x))) ++bad;
        }
    }
    ok (bad == 0, "Transistor: y(-x) == -y(x) exactly over a sweep (" + std::to_string (bad) + " differ)");
}

// Transistor against a double evaluation of u/(1+u⁴)^(1/4), normalised the same way. |u| <= k <= 62 here, so the
// clamp at 64 never acts.
void transistorPrecision()
{
    for (float db : { 0.0f, 0.1f, 1.0f, 6.0f, 12.0f, 36.0f })
    {
        const WS w = shaper (Shape::Transistor, kOf (db));
        const double k = (double) w.drive();
        auto sd = [] (double u) { return u / std::pow (1.0 + u * u * u * u, 0.25); };
        const double norm = 1.0 / sd (k);
        double worst = 0.0;
        for (int i = -100000; i <= 100000; ++i)
        {
            const float x = (float) i / 100000.0f;
            worse (worst, std::fabs ((double) w.processSample (x) - sd (k * (double) x) * norm));
        }
        okp (worst <= 1.0e-6, "Transistor driveDb " + num (db) + ": float within 1e-6 of u/(1+u^4)^(1/4) in double over [-1, 1] ("
                              + num (worst) + ")");
    }
}

//==============================================================================
// Through the whole Saturator: a bin-exact tone at os 4, mix 1, read over a settled window.
constexpr int kW = 16384, kWarm = 32768, kBin = 341;   // 341 · 48000 / 16384 = 999 Hz

std::vector<float> render (Shape s, float driveDb, double amp, float autoComp = 0.5f, float dcBlockHz = 10.0f)
{
    std::vector<float> x ((std::size_t) (kW + kWarm));
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] = (float) (amp * std::sin (2.0 * core::kPi * (double) kBin / kW * (double) i + 0.1));
    Sat st;
    Sat::Params p; p.shape = s; p.driveDb = driveDb; p.mix = 1.0f; p.autoComp = autoComp; p.dcBlockHz = dcBlockHz;
    st.setParams (p);
    ok (st.prepare (kFs, 512, 1, 4), "PRECONDITION: prepare");
    for (std::size_t o = 0; o < x.size(); o += 512)
    {
        float* io[1] { x.data() + o };
        felitronics::test::run (st.process (io, 1, (int) std::min<std::size_t> (512, x.size() - o)));
    }
    return std::vector<float> (x.begin() + kWarm, x.end());
}

// |X| at harmonic h of the tone, by a direct DFT over the window (the tone and its harmonics sit on bins).
double harmonic (const std::vector<float>& y, int h)
{
    double re = 0.0, im = 0.0;
    const double w = 2.0 * core::kPi * (double) (h * kBin) / kW;
    for (int i = 0; i < kW; ++i) { re += (double) y[(std::size_t) i] * std::cos (w * i); im -= (double) y[(std::size_t) i] * std::sin (w * i); }
    return 2.0 * std::hypot (re, im) / kW;
}

double dbc (const std::vector<float>& y, int h) { return 20.0 * std::log10 (harmonic (y, h) / harmonic (y, 1) + 1e-300); }

void harmonics()
{
    std::printf ("       H2..H7 in dBc (os 4, mix 1, autoComp 0.5, 999 Hz at 48 kHz):\n");
    for (Shape s : { Shape::Tanh, Shape::Tube, Shape::Transistor })
        for (float db : { 3.0f, 6.0f, 12.0f })
            for (double lvl : { -6.0, -3.0, 0.0 })
            {
                const std::vector<float> y = render (s, db, std::pow (10.0, lvl / 20.0));
                std::printf ("       %-10s drive %4.1f dB, %+4.1f dBFS:", name (s), (double) db, lvl);
                for (int h = 2; h <= 7; ++h) std::printf (" %8.2f", dbc (y, h));
                std::printf ("\n");
            }

    const std::vector<float> tube = render (Shape::Tube, 6.0f, std::pow (10.0, -3.0 / 20.0));
    const double h2 = dbc (tube, 2), h3 = dbc (tube, 3);
    okp (h2 - h3 >= 8.0, "Tube, driveDb 6, -3 dBFS: H2 leads H3 by >= 8 dB (H2 " + num (h2) + ", H3 " + num (h3) + " dBc)");

    const std::vector<float> tr = render (Shape::Transistor, 6.0f, std::pow (10.0, -3.0 / 20.0));
    const double t3 = dbc (tr, 3);
    okp (t3 >= -37.0 && t3 <= -33.0, "Transistor, driveDb 6, -3 dBFS: H3 in [-37, -33] dBc (" + num (t3) + ")");
    for (double lvl : { -3.0, 0.0 })
    {
        float cross = -1.0f;
        for (int i = 0; i <= 40 && cross < 0.0f; ++i)
        {
            const float db = 3.0f + 0.25f * (float) i;
            const double amp = std::pow (10.0, lvl / 20.0);
            if (dbc (render (Shape::Transistor, db, amp), 3) > dbc (render (Shape::Tanh, db, amp), 3)) cross = db;
        }
        okp (cross > 6.0f, "Transistor's H3 stays under Tanh's up to driveDb " + num (cross - 0.25f) + " at " + num (lvl)
                           + " dBFS (first above at " + num (cross) + ")");
    }
    double evenWorst = -1e9;
    for (int h = 2; h * kBin < kW / 2; h += 2) worse (evenWorst, dbc (tr, h));
    okp (evenWorst < -120.0, "Transistor: every even harmonic below Nyquist under -120 dBc (worst " + num (evenWorst) + ")");
}

//==============================================================================
// 8) Tube's DC goes to the blocker. The window holds whole periods of every harmonic, so its mean is the DC.
void tubeDc()
{
    const std::vector<float> y = render (Shape::Tube, 6.0f, 1.0);
    double sum = 0.0;
    for (float v : y) sum += (double) v;
    const double mean = sum / (double) y.size();
    okp (std::fabs (mean) < 1.0e-3, "Tube, full-scale sine, driveDb 6, blocker on: |mean| < 1e-3 over the settled tail (" + num (mean) + ")");
}

//==============================================================================
// 10) Tube's autoComp overshoot.
// The blocker takes the curve's DC out, which moves the settled waveform toward the smaller positive side; with it
// off (dcBlockHz 0) the output is the curve's own negative peak, which is also what a blocker that has not settled
// yet — a transient — lets through.
void tubeOvershoot()
{
    for (float hz : { 10.0f, 0.0f })
    {
        double worst = 0.0; float at = 0.0f;
        for (int i = 0; i <= 24; ++i)
        {
            const float db = 0.5f * (float) i;
            const std::vector<float> y = render (Shape::Tube, db, 1.0, 1.0f, hz);
            double m = 0.0;
            for (float v : y) worse (m, (double) std::fabs (v));
            if (! std::isnan (worst) && ! (m <= worst)) { worst = m; at = db; }
        }
        okp (worst <= 1.07, std::string ("Tube, autoComp 1, full-scale sine, driveDb 0..12, blocker ") + (hz > 0.0f ? "on" : "off")
                            + ": max |out| " + num (worst, "%.5f") + " (" + num (20.0 * std::log10 (worst), "%+.3f")
                            + " dB) at driveDb " + num (at) + ", pinned <= 1.07");
    }
}
} // namespace

int main (int argc, char** argv)
{
    if (argc > 1 && std::strcmp (argv[1], "exhaustive") == 0)
    {
        std::printf ("felitronics::saturation — Transistor's curve over every float in [0, 64]\n");
        group ("Transistor: s(u) <= 1, s(64) == 1, odd, monotonicity reported");
        transistorExhaustive();
        return felitronics::test::report();
    }
    std::printf ("felitronics::saturation — the Tube and Transistor curves\n");
    group ("shapeAt<S>(coeffs(), x) is processSample(x), bit for bit");
    shapeAtIsProcessSample<Shape::Tube>();
    shapeAtIsProcessSample<Shape::Transistor>();
    group ("slopeAtZero is the curve's derivative at 0");
    slopeAtZero (Shape::Tube);
    slopeAtZero (Shape::Transistor);
    group ("peak normalisation");
    peakNormalisation (Shape::Tube);
    peakNormalisation (Shape::Transistor);
    group ("Tube's precision and its sensitivity to tanh");
    tubePrecision();
    group ("Transistor at huge inputs, and its symmetry");
    transistorHuge();
    group ("Transistor's precision");
    transistorPrecision();
    group ("the harmonic signature through the Saturator");
    harmonics();
    group ("Tube's DC is removed");
    tubeDc();
    group ("Tube's autoComp overshoot");
    tubeOvershoot();
    return felitronics::test::report();
}
