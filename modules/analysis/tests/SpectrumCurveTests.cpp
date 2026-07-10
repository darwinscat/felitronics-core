// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis::offline — oracle + adversarial tests for logMagnitudeCurve / interferenceDb.
// Oracles are reference NULLs (a delta's spectrum is flat; a windowed tone peaks at its frequency; the
// interference formula is recomputed by hand). Adversarial cases feed every degenerate input the
// architecture consilium (codex + deepseek + gemini) flagged.

#include <felitronics_test.h>
#include <felitronics/analysis/offline/SpectrumCurve.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
namespace off = felitronics::analysis::offline;

namespace
{
const double kInf = std::numeric_limits<double>::infinity();
const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool allFinite (const std::vector<double>& x)
{
    for (double v : x) if (! std::isfinite (v)) return false;
    return true;
}
// Grid frequency of point i (must match the header's f = fLo·(fHi/fLo)^(i/(points-1))).
double gridF (int i, int points, double fLo, double fHi)
{
    return fLo * std::pow (fHi / fLo, (double) i / (double) (points - 1));
}
}

int main()
{
    std::printf ("felitronics::analysis::offline spectrum-curve tests\n");

    // ---- ORACLE 1: a delta's magnitude spectrum is flat → the curve is flat. ----
    group ("oracle: delta → flat curve");
    {
        std::vector<float> ir (16384, 0.0f); ir[0] = 1.0f;
        off::LogCurveSpec s; s.normalize = false;
        const auto c = off::logMagnitudeCurve (std::span<const float> (ir), 48000.0, s);
        ok (c.size() == 256, "256 points");
        double lo = 1e9, hi = -1e9;
        for (double v : c) { lo = std::min (lo, v); hi = std::max (hi, v); }
        ok (std::fabs (hi - 0.0) < 1e-6 && std::fabs (lo - 0.0) < 1e-6, "|FFT(δ)|=1 → 0 dB flat (unnormalized)");
    }

    // ---- ORACLE 2: a scaled delta is flat at 20·log10(A); normalize removes the scale. ----
    group ("oracle: scaled delta level + normalize");
    {
        std::vector<float> ir (8192, 0.0f); ir[0] = 2.0f;            // |FFT| = 2 → 20log10(2) = 6.0206 dB
        off::LogCurveSpec s; s.normalize = false;
        const auto c = off::logMagnitudeCurve (std::span<const float> (ir), 44100.0, s);
        approx (c[128], 20.0 * std::log10 (2.0), 1e-6, "flat at 20log10(2) dB");
        off::LogCurveSpec sn; sn.normalize = true;
        const auto cn = off::logMagnitudeCurve (std::span<const float> (ir), 44100.0, sn);
        approx (cn[128], 0.0, 1e-6, "normalize → peak 0 dB");
    }

    // ---- ORACLE 3: a windowed tone peaks at its frequency (log-grid mapping). ----
    group ("oracle: tone peaks at its frequency");
    {
        const double sr = 48000.0, f0 = 2000.0;
        const std::size_t N = 32768;
        std::vector<float> ir (N);
        for (std::size_t n = 0; n < N; ++n)
        {
            const double w = 0.5 * (1.0 - std::cos (2.0 * core::kPi * (double) n / (double) (N - 1)));  // Hann
            ir[n] = (float) (w * std::cos (2.0 * core::kPi * f0 * (double) n / sr));
        }
        const auto c = off::logMagnitudeCurve (std::span<const float> (ir), sr);
        int pk = 0; for (int i = 1; i < (int) c.size(); ++i) if (c[(std::size_t) i] > c[(std::size_t) pk]) pk = i;
        const double pf = gridF (pk, 256, 20.0, 20000.0);
        ok (std::fabs (pf - f0) / f0 < 0.05, "peak within 5% of 2 kHz");
    }

    // ---- ORACLE 4: interferenceDb recomputed by hand. ----
    group ("oracle: interferenceDb formula");
    {
        std::vector<double> a (8, 0.0), b (8, 0.0);                  // two mics at 0 dB → power 1 each
        std::vector<off::MicCurveView> mics { { std::span<const double> (a), true }, { std::span<const double> (b), true } };
        const double incoh = 10.0 * std::log10 (2.0);               // 3.0103 dB incoherent sum

        std::vector<double> cohReinforce (8, 20.0 * std::log10 (2.0));  // in-phase → +6.02 dB
        const auto r = off::interferenceDb (std::span<const double> (cohReinforce), std::span<const off::MicCurveView> (mics));
        approx (r[4], 20.0 * std::log10 (2.0) - incoh, 1e-9, "reinforcement = +6.02 − 3.01 dB");

        std::vector<double> cohCancel (8, -100.0);                   // deep cancellation
        const auto k = off::interferenceDb (std::span<const double> (cohCancel), std::span<const off::MicCurveView> (mics));
        approx (k[4], -100.0 - incoh, 1e-9, "cancellation = coherent − incoherent");
    }

    // ---- ORACLE 5: audibility mask + no-audible → 0. ----
    group ("oracle: audibility mask");
    {
        std::vector<double> a (4, 0.0), b (4, 0.0);
        std::vector<double> coh (4, 5.0);
        std::vector<off::MicCurveView> oneMuted { { std::span<const double> (a), true }, { std::span<const double> (b), false } };
        const auto r = off::interferenceDb (std::span<const double> (coh), std::span<const off::MicCurveView> (oneMuted));
        approx (r[0], 5.0 - 0.0, 1e-9, "muted mic drops out → incoherent = 1 mic (0 dB)");

        std::vector<off::MicCurveView> allMuted { { std::span<const double> (a), false } };
        const auto z = off::interferenceDb (std::span<const double> (coh), std::span<const off::MicCurveView> (allMuted));
        approx (z[0], 0.0, 0.0, "no audible mic → 0 dB (no signal, no interference)");
    }

    // ---- ADVERSARIAL: every degenerate input the consilium flagged. ----
    group ("adversarial: degenerate inputs");
    {
        std::vector<float> big (16384, 0.001f); for (std::size_t i = 0; i < big.size(); ++i) big[i] = (float) std::sin (0.05 * (double) i);

        ok (off::logMagnitudeCurve (std::span<const float> {}, 48000.0).empty(), "empty IR → empty");
        std::vector<float> tiny (4, 1.0f);
        ok (off::logMagnitudeCurve (std::span<const float> (tiny), 48000.0).empty(), "IR<8 → empty");
        ok (off::logMagnitudeCurve (std::span<const float> (big), 0.0).empty(), "sr=0 → empty");
        ok (off::logMagnitudeCurve (std::span<const float> (big), -48000.0).empty(), "sr<0 → empty");
        ok (off::logMagnitudeCurve (std::span<const float> (big), kNaN).empty(), "sr=NaN → empty");
        ok (off::logMagnitudeCurve (std::span<const float> (big), kInf).empty(), "sr=Inf → empty");

        off::LogCurveSpec bad; bad.fLo = -5.0;
        ok (off::logMagnitudeCurve (std::span<const float> (big), 48000.0, bad).empty(), "fLo<=0 → empty");
        off::LogCurveSpec inv; inv.fLo = 5000.0; inv.fHi = 1000.0;
        ok (off::logMagnitudeCurve (std::span<const float> (big), 48000.0, inv).empty(), "fHi<=fLo → empty");

        off::LogCurveSpec p1; p1.points = 1;                        // was: 1/(points-1) division by zero
        const auto c1 = off::logMagnitudeCurve (std::span<const float> (big), 48000.0, p1);
        ok (c1.size() == 2 && allFinite (c1), "points=1 clamped to 2, finite");

        off::LogCurveSpec hi; hi.fHi = 1.0e9;                       // > Nyquist → clamped, no OOB
        const auto ch = off::logMagnitudeCurve (std::span<const float> (big), 48000.0, hi);
        ok (ch.size() == 256 && allFinite (ch), "fHi>Nyquist clamped, finite");

        off::LogCurveSpec oz; oz.octaveDivisions = 0.0;             // was: pow(2, 1/0) = inf band
        ok (allFinite (off::logMagnitudeCurve (std::span<const float> (big), 48000.0, oz)), "octaveDivisions=0 healed, finite");

        off::LogCurveSpec ot; ot.octaveDivisions = 1e-300;         // review C1: tiny+ → halfBand=inf → (int)ceil(inf) was UB
        ok (allFinite (off::logMagnitudeCurve (std::span<const float> (big), 48000.0, ot)), "octaveDivisions=1e-300 clamped, no inf→int UB");
        off::LogCurveSpec oi; oi.octaveDivisions = kInf;
        ok (allFinite (off::logMagnitudeCurve (std::span<const float> (big), 48000.0, oi)), "octaveDivisions=Inf clamped, finite");

        // review C2: a non-pow2 minNfft must give the SAME curve as its pow2 round-up (binHz↔bins consistent).
        off::LogCurveSpec mnp; mnp.minNfft = 5000; mnp.normalize = false;   // rounds up to 8192
        off::LogCurveSpec mp;  mp.minNfft  = 8192; mp.normalize  = false;
        const auto cnp = off::logMagnitudeCurve (std::span<const float> (big), 48000.0, mnp);
        const auto cp  = off::logMagnitudeCurve (std::span<const float> (big), 48000.0, mp);
        bool binMatch = cnp.size() == cp.size();
        for (std::size_t i = 0; binMatch && i < cnp.size(); ++i) binMatch = std::fabs (cnp[i] - cp[i]) < 1e-12;
        ok (binMatch, "non-pow2 minNfft(5000) == pow2 minNfft(8192): no bin↔Hz skew");

        // review D2: a huge-but-finite sample must not overflow ΣM² to inf → NaN in normalize. Uses the
        // DOUBLE overload — float caps at ~3.4e38 (unreachable), but double 1e155 makes |FFT|²≈1e318=inf.
        std::vector<double> huge (16384, 0.0);
        for (std::size_t i = 0; i < huge.size(); ++i) huge[i] = 1.0e155 * std::sin (0.04 * (double) i);
        ok (allFinite (off::logMagnitudeCurve (std::span<const double> (huge), 48000.0)), "huge-finite input → finite curve (no ΣM² overflow NaN)");

        std::vector<float> nan (16384, 0.0f);                       // NaN/Inf samples neutralized to 0
        for (std::size_t i = 0; i < nan.size(); ++i) nan[i] = (float) std::sin (0.03 * (double) i);
        nan[100] = (float) kNaN; nan[200] = (float) kInf; nan[300] = -(float) kInf;
        ok (allFinite (off::logMagnitudeCurve (std::span<const float> (nan), 48000.0)), "NaN/Inf input → finite curve");

        std::vector<float> silent (16384, 0.0f);                    // silent + normalize must NOT flat-line at 0 dB
        const auto cs = off::logMagnitudeCurve (std::span<const float> (silent), 48000.0);
        ok (allFinite (cs) && cs[0] < -100.0, "silent IR stays at the floor, not 0 dB");

        // interferenceDb: a shorter mic curve than the coherent curve must not read OOB.
        std::vector<double> a (3, 0.0), coh (10, 0.0);
        std::vector<off::MicCurveView> mm { { std::span<const double> (a), true } };
        const auto ri = off::interferenceDb (std::span<const double> (coh), std::span<const off::MicCurveView> (mm));
        ok (ri.size() == 10 && allFinite (ri), "mic curve shorter than coherent → no OOB, finite");
        approx (ri[5], 0.0, 0.0, "past the mic curve → 0 (no audible mic at p)");

        // double overload compiles + matches the float overload for the same data.
        std::vector<double> bigd (big.begin(), big.end());
        const auto cf = off::logMagnitudeCurve (std::span<const float> (big), 48000.0);
        const auto cd = off::logMagnitudeCurve (std::span<const double> (bigd), 48000.0);
        bool same = cf.size() == cd.size();
        for (std::size_t i = 0; same && i < cf.size(); ++i) same = std::fabs (cf[i] - cd[i]) < 1e-9;
        ok (same, "float and double overloads agree");
    }

    return felitronics::test::report();
}
