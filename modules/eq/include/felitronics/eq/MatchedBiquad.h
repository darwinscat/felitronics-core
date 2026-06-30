// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>

#include <algorithm>
#include <cmath>
#include <complex>

//==============================================================================
// teq::matched — Vicanek matched-filter coefficient design. The matched design fits
// the analog prototype's MAGNITUDE near Nyquist far better than the RBJ cookbook
// (whose bilinear transform "cramps" high bells, steepens high shelves, narrows
// resonances). Same biquad runtime cost — just better coefficients.
//
// References (formulas cited inline):
//   Martin Vicanek, "Matched Second Order Digital Filters" (2016)
//       — bell / lowpass / highpass / bandpass.
//   Martin Vicanek, "Matched Two-Pole Digital Shelving Filters" (2024-2025)
//       — low / high shelf (2-pole Butterworth, no Q).
//
// JUCE-free (pure std) on purpose: this is a framework-agnostic module — drop the
// `teq/` folder into any plugin and `#include <felitronics/eq/MatchedBiquad.h>`.
//==============================================================================
namespace felitronics::eq
{

//==============================================================================
// One biquad. a0 normalised to 1; H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2).
struct BiquadCoeffs
{
    double b0 = 1.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;

    // |H(e^{jw})| at digital angular frequency w (rad/sample, 0..pi).
    double magnitude (double w) const noexcept
    {
        const std::complex<double> z1 = std::polar (1.0, -w);   // e^{-jw}
        const std::complex<double> z2 = z1 * z1;
        const std::complex<double> num = b0 + b1 * z1 + b2 * z2;
        const std::complex<double> den = 1.0 + a1 * z1 + a2 * z2;
        return std::abs (num / den);
    }

    double magnitudeDb (double w) const noexcept { return 20.0 * std::log10 (std::max (1e-12, magnitude (w))); }

    // Stable iff |a1| < 2 and |a1| - 1 < a2 < 1 (Vicanek eq. 3).
    bool isStable() const noexcept { return std::abs (a1) < 2.0 && (std::abs (a1) - 1.0) < a2 && a2 < 1.0; }
};

namespace detail
{
    inline double safeSqrt (double x) noexcept { return std::sqrt (std::max (0.0, x)); }

    // Common matched poles from impulse invariance (BiquadFits eq. 12). q = 1/(2Q).
    inline void matchedPoles (double w0, double Q, double& a1, double& a2) noexcept
    {
        const double q   = 1.0 / (2.0 * Q);
        const double eqw = std::exp (-q * w0);
        a1 = (q <= 1.0) ? -2.0 * eqw * std::cos  (w0 * std::sqrt (1.0 - q * q))
                        : -2.0 * eqw * std::cosh (w0 * std::sqrt (q * q - 1.0));
        a2 = std::exp (-2.0 * q * w0);
    }

    // Minimum-phase numerator from (B0,B1,B2) — BiquadFits eq. (29).
    inline void solveNumerator (double B0, double B1, double B2, BiquadCoeffs& c) noexcept
    {
        const double sB0 = safeSqrt (B0), sB1 = safeSqrt (B1);
        const double W   = 0.5 * (sB0 + sB1);
        c.b0 = 0.5 * (W + safeSqrt (W * W + B2));
        c.b1 = 0.5 * (sB0 - sB1);
        c.b2 = (c.b0 != 0.0) ? -B2 / (4.0 * c.b0) : 0.0;
    }

    // Shared solver for the 2-pole Butterworth shelf (2poleShelvingFits appendix), given the
    // Nyquist-normalised centre fc (= 2 f0/fs) and the working gain g. Returns the building
    // blocks v,w,a0,aa2,bb2; high/low shelf differ only in `g` and the b-scaling.
    struct ShelfTmp { double v, w, a0, aa2, bb2; };

    inline ShelfTmp shelfSolve (double fc, double g) noexcept
    {
        const double piHalf = kPi * 0.5;
        const double invg = 1.0 / g;
        const double fc2 = fc * fc, fc4 = fc2 * fc2;
        const double hny = (fc4 + g) / (fc4 + invg);                       // |H|^2 at Nyquist (eq. 10)

        const double f1 = fc / std::sqrt (0.160 + 1.543 * fc2);            // matching point 1 (eq. 11)
        const double f14 = f1 * f1 * f1 * f1;
        const double h1 = (fc4 + f14 * g) / (fc4 + f14 * invg);
        const double s1 = std::sin (piHalf * f1); const double phi1 = s1 * s1;

        const double f2 = fc / std::sqrt (0.947 + 3.806 * fc2);            // matching point 2 (eq. 11)
        const double f24 = f2 * f2 * f2 * f2;
        const double h2 = (fc4 + f24 * g) / (fc4 + f24 * invg);
        const double s2 = std::sin (piHalf * f2); const double phi2 = s2 * s2;

        const double d1 = (h1 - 1.0) * (1.0 - phi1);                       // linear system (eq. 12-13)
        const double c11 = -phi1 * d1, c12 = phi1 * phi1 * (hny - h1);
        const double d2 = (h2 - 1.0) * (1.0 - phi2);
        const double c21 = -phi2 * d2, c22 = phi2 * phi2 * (hny - h2);

        const double alfa1 = (c22 * d1 - c12 * d2) / (c11 * c22 - c12 * c21);   // (eq. 15)
        const double aa1 = (d1 - c11 * alfa1) / c12;
        const double bb1 = hny * aa1;

        ShelfTmp t;
        t.aa2 = 0.25 * (alfa1 - aa1);
        t.bb2 = 0.25 * (alfa1 - bb1);
        t.v = 0.5 * (1.0 + safeSqrt (aa1));
        t.w = 0.5 * (1.0 + safeSqrt (bb1));
        t.a0 = 0.5 * (t.v + safeSqrt (t.v * t.v + t.aa2));
        return t;
    }
}

namespace matched
{
    //==========================================================================
    // Peaking EQ (bell). gainLin = linear magnitude AT the centre (|H(w0)| = gainLin).
    // BiquadFits eq. (26),(27),(42)-(45),(29).
    inline BiquadCoeffs peaking (double f0, double fs, double Q, double gainLin) noexcept
    {
        // Symmetric (reciprocal) cut: design the mirror BOOST and invert 1/H, so -G and +G dB at the
        // same Q are exact vertical mirror images (the modern parametric norm). A matched bell has
        // gain-independent poles, so a directly-designed cut would be far narrower than its boost;
        // inverting fixes that. The boost numerator is minimum-phase (zeros inside the unit circle),
        // so the inverted cut is stable, and its magnitude is the exact reciprocal of a Nyquist-
        // matched response (still honest near Nyquist).
        if (gainLin < 1.0 && gainLin > 0.0)
        {
            const BiquadCoeffs b = peaking (f0, fs, Q, 1.0 / gainLin);
            const double inv = 1.0 / b.b0;
            return { inv, b.a1 * inv, b.a2 * inv, b.b1 * inv, b.b2 * inv };
        }

        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);

        const double s = std::sin (w0 * 0.5);
        const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
        const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
        const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

        const double G2 = gainLin * gainLin;
        const double B0 = A0;
        const double R1 = (A0 * phi0 + A1 * phi1 + A2 * phi2) * G2;
        const double R2 = (-A0 + A1 + 4.0 * (phi0 - phi1) * A2) * G2;
        const double B2 = (R1 - R2 * phi1 - B0) / (4.0 * phi1 * phi1);
        const double B1 = R2 + B0 + 4.0 * (phi1 - phi0) * B2;

        detail::solveNumerator (B0, B1, B2, c);
        return c;
    }

    inline BiquadCoeffs peakingDb (double f0, double fs, double Q, double gainDb) noexcept
    {
        return peaking (f0, fs, Q, std::pow (10.0, gainDb / 20.0));
    }

    //==========================================================================
    // Lowpass (resonant). Matches the analog |H(iw0)| = Q. Single zero (b2=0). eq. (30)-(34).
    inline BiquadCoeffs lowpass (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);

        const double s = std::sin (w0 * 0.5);
        const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
        const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
        const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

        const double B0 = A0;
        const double R1 = (A0 * phi0 + A1 * phi1 + A2 * phi2) * Q * Q;
        const double B1 = (R1 - B0 * phi0) / phi1;
        c.b0 = 0.5 * (t0 + detail::safeSqrt (B1));    // sqrt(B0) = 1+a1+a2 (eq. 34)
        c.b1 = t0 - c.b0;
        c.b2 = 0.0;
        return c;
    }

    //==========================================================================
    // Highpass. Double zero at z=1 (12 dB/oct). Matches analog |H(iw0)| = Q. eq. (35),(36).
    inline BiquadCoeffs highpass (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);

        const double s = std::sin (w0 * 0.5);
        const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
        const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
        const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

        c.b0 = Q * detail::safeSqrt (A0 * phi0 + A1 * phi1 + A2 * phi2) / (4.0 * phi1);
        c.b1 = -2.0 * c.b0;
        c.b2 = c.b0;
        return c;
    }

    //==========================================================================
    // Bandpass (unity gain at centre). Single zero at z=1. eq. (37)-(41).
    inline BiquadCoeffs bandpass (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);

        const double s = std::sin (w0 * 0.5);
        const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
        const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
        const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

        const double R1 = A0 * phi0 + A1 * phi1 + A2 * phi2;
        const double R2 = -A0 + A1 + 4.0 * (phi0 - phi1) * A2;
        const double B2 = (R1 - R2 * phi1) / (4.0 * phi1 * phi1);
        const double B1 = R2 + 4.0 * (phi1 - phi0) * B2;
        c.b1 = -0.5 * detail::safeSqrt (B1);
        c.b0 = 0.5 * (detail::safeSqrt (B2 + c.b1 * c.b1) - c.b1);
        c.b2 = -c.b0 - c.b1;
        return c;
    }

    //==========================================================================
    // Notch (band-stop): an infinite null at f0, unity at DC/Nyquist, width set by Q. Zeros sit on
    // the unit circle at +-w0; poles come from the matched-pole placement (Q-consistent with the
    // rest of the family), then the numerator is scaled to pass DC at unity.
    inline BiquadCoeffs notch (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);
        const double cw = std::cos (w0);
        const double k  = (1.0 + c.a1 + c.a2) / (2.0 - 2.0 * cw);
        c.b0 = k; c.b1 = -2.0 * k * cw; c.b2 = k;
        return c;
    }

    //==========================================================================
    // All-pass: |H| = 1 at every frequency (flat magnitude); the phase rotates 360° through f0 with
    // sharpness set by Q. Numerator is the reversed denominator, which forces unit magnitude.
    inline BiquadCoeffs allpass (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);
        c.b0 = c.a2; c.b1 = c.a1; c.b2 = 1.0;
        return c;
    }

    //==========================================================================
    // First-order low/high pass (6 dB/oct), bilinear. The matched 2nd-order machinery isn't needed
    // at this gentle slope (Nyquist cramping is negligible) — used as the odd section of a cascade.
    inline BiquadCoeffs lowpass1 (double f0, double fs) noexcept
    {
        BiquadCoeffs c;
        const double K = std::tan (kPi * f0 / fs), nrm = 1.0 / (1.0 + K);
        c.b0 = K * nrm; c.b1 = K * nrm; c.b2 = 0.0; c.a1 = (K - 1.0) * nrm; c.a2 = 0.0;
        return c;
    }
    inline BiquadCoeffs highpass1 (double f0, double fs) noexcept
    {
        BiquadCoeffs c;
        const double K = std::tan (kPi * f0 / fs), nrm = 1.0 / (1.0 + K);
        c.b0 = nrm; c.b1 = -nrm; c.b2 = 0.0; c.a1 = (K - 1.0) * nrm; c.a2 = 0.0;
        return c;
    }

    //==========================================================================
    // High shelf — matched 2-pole Butterworth. gainLin = the high-frequency plateau
    // (linear; |H| -> gainLin as f -> Nyquist+, 1.0 at DC). 2poleShelvingFits appendix A.1.
    inline BiquadCoeffs highShelf (double f0, double fs, double gainLin) noexcept
    {
        BiquadCoeffs c;
        const double fc = 2.0 * f0 / fs;                                   // normalised to Nyquist
        const double g  = (std::abs (1.0 - gainLin) < 1e-6) ? 1.00001 : gainLin;

        const auto t = detail::shelfSolve (fc, g);
        const double inva0 = 1.0 / t.a0;
        c.a1 = (1.0 - t.v) * inva0;
        c.a2 = -0.25 * t.aa2 * inva0 * inva0;
        c.b0 = (0.5 * (t.w + detail::safeSqrt (t.w * t.w + t.bb2))) * inva0;
        c.b1 = (1.0 - t.w) * inva0;
        c.b2 = (-0.25 * t.bb2 / c.b0) * inva0 * inva0;
        return c;
    }

    //==========================================================================
    // Low shelf — matched 2-pole Butterworth. gainLin = the low-frequency plateau
    // (linear; |H| -> gainLin at DC, 1.0 at high frequencies). 2poleShelvingFits appendix A.2.
    inline BiquadCoeffs lowShelf (double f0, double fs, double gainLin) noexcept
    {
        BiquadCoeffs c;
        const double fc = 2.0 * f0 / fs;
        const double g  = (std::abs (1.0 - gainLin) < 1e-6) ? 1.00001 : 1.0 / gainLin;
        const double invg = 1.0 / g;                                       // == gainLin

        const auto t = detail::shelfSolve (fc, g);
        const double inva0 = 1.0 / t.a0;
        const double ginva0 = invg * inva0;
        c.a1 = (1.0 - t.v) * inva0;
        c.a2 = -0.25 * t.aa2 * inva0 * inva0;
        const double b0raw = 0.5 * (t.w + detail::safeSqrt (t.w * t.w + t.bb2));
        c.b1 = (1.0 - t.w) * ginva0;
        c.b2 = (-0.25 * t.bb2 / b0raw) * ginva0;
        c.b0 = b0raw * ginva0;
        return c;
    }

    inline BiquadCoeffs highShelfDb (double f0, double fs, double gainDb) noexcept { return highShelf (f0, fs, std::pow (10.0, gainDb / 20.0)); }
    inline BiquadCoeffs lowShelfDb  (double f0, double fs, double gainDb) noexcept { return lowShelf  (f0, fs, std::pow (10.0, gainDb / 20.0)); }

    //==========================================================================
    // Resonant low/high shelf with a quality factor Q (the "expensive EQ" shelf). The poles are the
    // matched (w0, Q) pair shared with the bell / lowpass / highpass family; the numerator is fit so
    // the magnitude is EXACT at DC, Nyquist and the corner. The analog prototype has poles at (w0, Q)
    // and zeros at (w0*sqrt(A), Q), giving |H(w0)|^2 = Q^2 (A-1)^2 + A: a low Q is a gentle
    // (Butterworth-ish) shelf, a high Q overshoots (boost) / dips (cut) at the transition, like
    // Pro-Q / Neutron. Q ~ 0.707 reproduces a clean gentle shelf. Matched at Nyquist => high shelves
    // don't cramp. `high`=false -> low shelf (gain at DC); true -> high shelf (gain at Nyquist).
    inline BiquadCoeffs shelf (double f0, double fs, double gainLin, double Q, bool high) noexcept
    {
        // Cut = reciprocal of the mirror boost -> an exact mirror image (the resonance dip mirrors
        // the bump), the same trick the matched bell uses. If that inversion is unstable (a few
        // high-shelf cuts whose mirror-boost numerator goes non-minimum-phase right at Nyquist), fall
        // through to the direct design below, which is ALWAYS pole-stable (matched poles) and gives a
        // gentle dip at the very edge — where edge-of-band resonance is inaudible anyway.
        if (gainLin < 1.0 && gainLin > 0.0)
        {
            const BiquadCoeffs b = shelf (f0, fs, 1.0 / gainLin, Q, high);
            const double inv = 1.0 / b.b0;
            const BiquadCoeffs cut { inv, b.a1 * inv, b.a2 * inv, b.b1 * inv, b.b2 * inv };
            if (cut.isStable()) return cut;
        }

        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs;
        detail::matchedPoles (w0, Q, c.a1, c.a2);                         // matched poles -> always pole-stable

        const double s = std::sin (w0 * 0.5);
        const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
        const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
        const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

        const double A   = gainLin, A2g = A * A;
        const double B0  = high ? A0       : A2g * A0;                    // |H(DC)|^2  = (high ? 1 : A)^2
        const double B1  = high ? A2g * A1 : A1;                          // |H(Nyq)|^2 = (high ? A : 1)^2
        const double Hc2 = Q * Q * (A - 1.0) * (A - 1.0) + A;             // analog |H(w0)|^2
        const double denW0 = A0 * phi0 + A1 * phi1 + A2 * phi2;
        const double B2  = (Hc2 * denW0 - B0 * phi0 - B1 * phi1) / phi2;

        detail::solveNumerator (B0, B1, B2, c);
        return c;
    }

    inline BiquadCoeffs lowShelfQ    (double f0, double fs, double gainLin, double Q) noexcept { return shelf (f0, fs, gainLin, Q, false); }
    inline BiquadCoeffs highShelfQ   (double f0, double fs, double gainLin, double Q) noexcept { return shelf (f0, fs, gainLin, Q, true); }
    inline BiquadCoeffs lowShelfQDb  (double f0, double fs, double gainDb,  double Q) noexcept { return lowShelfQ  (f0, fs, std::pow (10.0, gainDb / 20.0), Q); }
    inline BiquadCoeffs highShelfQDb (double f0, double fs, double gainDb,  double Q) noexcept { return highShelfQ (f0, fs, std::pow (10.0, gainDb / 20.0), Q); }
}

//==============================================================================
// RBJ "Audio EQ Cookbook" designs — kept ONLY as the cramping baseline the unit tests
// measure matched against. Not used in the live signal path.
namespace rbj
{
    inline BiquadCoeffs peaking (double f0, double fs, double Q, double gainLin) noexcept
    {
        BiquadCoeffs c;
        const double A = std::sqrt (gainLin);
        const double w0 = 2.0 * kPi * f0 / fs, cw = std::cos (w0), alpha = std::sin (w0) / (2.0 * Q);
        const double a0 = 1.0 + alpha / A;
        c.b0 = (1.0 + alpha * A) / a0; c.b1 = (-2.0 * cw) / a0; c.b2 = (1.0 - alpha * A) / a0;
        c.a1 = (-2.0 * cw) / a0;       c.a2 = (1.0 - alpha / A) / a0;
        return c;
    }

    inline BiquadCoeffs lowpass (double f0, double fs, double Q) noexcept
    {
        BiquadCoeffs c;
        const double w0 = 2.0 * kPi * f0 / fs, cw = std::cos (w0), alpha = std::sin (w0) / (2.0 * Q);
        const double a0 = 1.0 + alpha;
        c.b0 = ((1.0 - cw) * 0.5) / a0; c.b1 = (1.0 - cw) / a0; c.b2 = ((1.0 - cw) * 0.5) / a0;
        c.a1 = (-2.0 * cw) / a0;        c.a2 = (1.0 - alpha) / a0;
        return c;
    }

    inline BiquadCoeffs highShelf (double f0, double fs, double gainLin) noexcept   // slope S = 1
    {
        BiquadCoeffs c;
        const double A = std::sqrt (gainLin);
        const double w0 = 2.0 * kPi * f0 / fs, cw = std::cos (w0), sw = std::sin (w0);
        const double alpha = (sw * 0.5) * std::sqrt ((A + 1.0 / A) * (1.0 / 1.0 - 1.0) + 2.0);
        const double tsa = 2.0 * std::sqrt (A) * alpha;
        const double a0 = (A + 1.0) - (A - 1.0) * cw + tsa;
        c.b0 =  A * ((A + 1.0) + (A - 1.0) * cw + tsa) / a0;
        c.b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cw) / a0;
        c.b2 =  A * ((A + 1.0) + (A - 1.0) * cw - tsa) / a0;
        c.a1 =  2.0 * ((A - 1.0) - (A + 1.0) * cw) / a0;
        c.a2 = ((A + 1.0) - (A - 1.0) * cw - tsa) / a0;
        return c;
    }
}

//==============================================================================
// Direct-form-II transposed runtime processor (float state). For the live path.
struct Biquad
{
    BiquadCoeffs c;
    float z1 = 0.0f, z2 = 0.0f;

    void reset() noexcept { z1 = z2 = 0.0f; }
    void setCoeffs (const BiquadCoeffs& nc) noexcept { c = nc; }

    inline float processSample (float x) noexcept
    {
        const double y = c.b0 * x + z1;
        z1 = (float) (c.b1 * x - c.a1 * y + z2);
        z2 = (float) (c.b2 * x - c.a2 * y);
        return (float) y;
    }

    // Per-block denormal guard — see Svf::flushDenormals.
    void flushDenormals() noexcept
    {
        if (std::fabs (z1) < 1e-15f) z1 = 0.0f;
        if (std::fabs (z2) < 1e-15f) z2 = 0.0f;
    }
};

} // namespace felitronics::eq
