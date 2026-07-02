// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqTypes.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

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

    // Feasibility of the 3-point (DC, corner, Nyquist) numerator fit: real min-phase b0,b2 exist iff
    // W²+B2 ≥ 0 (W = ½(√B0+√B1)). When it fails, the fit is infeasible — forcing it puts the numerator
    // zeros on/outside the unit circle (a spurious passband notch in a boost; a near-unit-circle pole,
    // hence a huge peak, in a reciprocal cut). Callers fall back to the non-resonant shelf there.
    inline bool shelfFitFeasible (double B0, double B1, double B2) noexcept
    {
        const double W = 0.5 * (safeSqrt (B0) + safeSqrt (B1));
        return W * W + B2 > 0.0;
    }

    // Exact min & max of the quadratic ratio (n2·t²+n1·t+n0)/(d2·t²+d1·t+d0) over t∈[lo,hi] (hi may be
    // +∞): endpoints + the ≤2 stationary points (roots of (n2d1−n1d2)t²+2(n2d0−n0d2)t+(n1d0−n0d1)=0).
    // Denominator assumed > 0 on the interval. Used to bound a biquad's |H|² (t=cos w∈[−1,1]) and the
    // analog resonant-shelf |H|² (t=(f/f0)²∈[0,∞)) in closed form — no frequency sampling.
    struct RatioBound { double lo, hi; };
    inline RatioBound ratioExtent (double n2, double n1, double n0, double d2, double d1, double d0,
                                   double lo, double hi) noexcept
    {
        auto val = [&] (double t) noexcept { return (n2 * t * t + n1 * t + n0) / (d2 * t * t + d1 * t + d0); };
        double mn = 1e300, mx = -1e300;
        auto upd = [&] (double t) noexcept { if (t >= lo && t <= hi) { const double v = val (t); mn = std::min (mn, v); mx = std::max (mx, v); } };
        upd (lo);
        if (std::isfinite (hi)) upd (hi);
        else if (d2 != 0.0) { const double v = n2 / d2; mn = std::min (mn, v); mx = std::max (mx, v); }   // t→∞
        const double a = n2 * d1 - n1 * d2, b = 2.0 * (n2 * d0 - n0 * d2), c = n1 * d0 - n0 * d1;
        if (std::abs (a) < 1e-300) { if (std::abs (b) > 1e-300) upd (-c / b); }
        else { const double disc = b * b - 4.0 * a * c; if (disc >= 0.0) { const double s = std::sqrt (disc); upd ((-b + s) / (2.0 * a)); upd ((-b - s) / (2.0 * a)); } }
        return { mn, mx };
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

        double B1, B2;
        if (w0 < 0.02)
        {
            // Tiny-w0 branch. The classic evaluation dies of catastrophic cancellation here: R1 and
            // R2 cancel to O(q²w0⁴)/O(q²w0²), then R1−R2·φ1 cancels again to O(q²w0⁶) — below double
            // precision for f0 ≲ 50 Hz at high Q, driving B1 negative → safeSqrt clamps b1 to −0 and
            // the CENTRE GAIN blows up (+28.6 dB at fs=192k, f0=10, Q=40; the formulas are exact,
            // the double arithmetic isn't). Evaluate B1,B2 by their Maclaurin series in w0 instead
            // (polynomials in q=1/(2Q) — valid for both the cos and cosh pole branches). Verified:
            // centre unity within 1.9e-4 dB over fs{44.1..192}k × Q[0.05,40] × f0[10,200 Hz], seam
            // jump ≤ 1.5e-5 dB at the gate; above the gate the classic path below is byte-identical.
            const double q = 1.0 / (2.0 * Q);
            const double q2 = q * q, q3 = q2 * q, q4 = q2 * q2, q5 = q4 * q, q6 = q4 * q2, q7 = q6 * q, q8 = q4 * q4;
            const double x = w0, x2 = x * x, x3 = x2 * x, x4 = x2 * x2, x5 = x4 * x, x6 = x4 * x2;
            B2 = x2 * ( (4.0 /  3.0) * q2
               + x  * (-8.0 /  3.0) * q3
               + x2 * ( 2.0 / 45.0  * q2 + 128.0 / 45.0 * q4)
               + x3 * (-4.0 / 45.0  * q3 -  32.0 / 15.0 * q5)
               + x4 * ( 1.0 / 630.0 * q2 +  94.0 / 945.0 * q4 + 44.0 / 35.0 * q6)
               + x5 * (-1.0 / 315.0 * q3 -  76.0 / 945.0 * q5 - 584.0 / 945.0 * q7)
               + x6 * ( 1.0 / 18900.0 * q2 + 209.0 / 56700.0 * q4 + 734.0 / 14175.0 * q6 + 3728.0 / 14175.0 * q8));
            B1 = x2 * ( (32.0 /  3.0) * q2
               + x  * (-64.0 /  3.0) * q3
               + x2 * (-128.0 / 45.0 * q2 + 1168.0 / 45.0 * q4)
               + x3 * ( 256.0 / 45.0 * q3 -  352.0 / 15.0 * q5)
               + x4 * (  82.0 / 315.0 * q2 - 6088.0 / 945.0 * q4 + 1072.0 / 63.0 * q6)
               + x5 * (-164.0 / 315.0 * q3 + 5008.0 / 945.0 * q5 - 1952.0 / 189.0 * q7)
               + x6 * ( -46.0 / 4725.0 * q2 + 8131.0 / 14175.0 * q4 - 49376.0 / 14175.0 * q6 + 76528.0 / 14175.0 * q8));
        }
        else
        {
            const double s = std::sin (w0 * 0.5);
            const double phi1 = s * s, phi0 = 1.0 - phi1, phi2 = 4.0 * phi0 * phi1;
            const double t0 = 1.0 + c.a1 + c.a2, t1 = 1.0 - c.a1 + c.a2;
            const double A0 = t0 * t0, A1 = t1 * t1, A2 = -4.0 * c.a2;

            const double R1 = A0 * phi0 + A1 * phi1 + A2 * phi2;
            const double R2 = -A0 + A1 + 4.0 * (phi0 - phi1) * A2;
            B2 = (R1 - R2 * phi1) / (4.0 * phi1 * phi1);
            B1 = R2 + 4.0 * (phi1 - phi0) * B2;
        }
        c.b1 = -0.5 * detail::safeSqrt (B1);
        c.b0 = 0.5 * (detail::safeSqrt (B2 + c.b1 * c.b1) - c.b1);
        c.b2 = -c.b0 - c.b1;
        return c;
    }

    //==========================================================================
    // Notch (band-stop): an infinite null at f0, unity at DC, width set by Q. Zeros sit on the
    // unit circle at +-w0; poles come from the matched-pole placement (Q-consistent with the rest
    // of the family), then the numerator is scaled to pass DC at unity. NOT unity at Nyquist for a
    // wide and/or high notch: the analog prototype itself keeps residual attenuation at fs/2
    // (|H|² = 1/(1+x²), x = BW·Ω/(Ω0²−Ω²)), and the heavily-damped matched poles under-shoot even
    // that, roughly DOUBLING the analog dB (f0=20 kHz, Q=0.5, fs=48 kHz: −25.4 dB at Nyquist vs
    // the analog −14.9; a narrow mid-band notch is back to ~0 dB well before fs/2). FROZEN legacy
    // contract — sessions rely on these exact bytes; the analog-faithful near-Nyquist behaviour
    // lives in notchCascade (sections ≥ 2), which alias-corrects its stagger poles.
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
    // Variable-order matched band-stop (Butterworth). `sections` (1..8) is the Butterworth LP
    // prototype order AND the biquad count: the analog low-pass→band-stop transform
    //     s_lp → BW·s / (s² + Ω0²),   Ω0 = 2π f0,  BW = Ω0 / Q
    // maps each prototype pole p to a quadratic  s² − (BW/p)·s + Ω0² = 0  whose two roots are a
    // geometric (log-symmetric) stagger pair about f0 (product of roots = Ω0²). Every section pins
    // its ZEROS at ±w0 (numerator normalised to unity at DC, exactly as matched::notch), so the
    // cascade keeps an order-fold INFINITE null that never drifts, a maximally-flat passband, and
    // skirts that steepen with `sections`; only the POLES stagger. Q stays the overall −3 dB
    // BANDWIDTH (Q = f0/BW), independent of order.
    //
    // Poles are placed by the matched-Z map z = exp(sT) — algebraically identical to
    // detail::matchedPoles — and because the LP→BS transform keeps every pole in the left
    // half-plane (roots share the sign of Re, and their sum has Re<0 ⇒ both in the LHP) each z
    // sits strictly inside the unit circle ⇒ every section is stable by construction. BUT exp()
    // WRAPS once Im(s)·T > π: a wide and/or high notch throws upper stagger poles ABOVE Nyquist
    // (only the upper root of a pair can — |r1||r2| = Ω0² pins the lower one below Ω0 < π/T),
    // and the alias lands mid-band where it stops providing the analog upper-side recovery. With
    // DC-normalised numerators the whole top end then sags to ≈ DOUBLE the analog dB (f0=15 kHz,
    // Q=0.5, fs=48 kHz, 4 sections: −49.5 dB at Nyquist vs the analog −25.0); heavily-damped
    // sections just below the wrap add a few dB of the same sag. Alias-gated repair, in order:
    //   1. REFIT every wrapped section: zeros stay pinned and DC stays unity, while the pole
    //      pair (2 DOF) is solved in closed form so the section interpolates its own analog
    //      magnitude — times the computable matched-Z error of the untouched sections, split
    //      evenly — EXACTLY at Nyquist and at wm = √(w0·π). The cascade then equals the analog
    //      band-stop at both points and tracks it between (the case above: exact at fs/2,
    //      ≤ 1.3 dB broadband). The Nyquist condition is linear in (u,v) = (−a1/2, a2), the wm
    //      condition reduces to one quadratic in s = 1+v; of two valid roots prefer the pole
    //      product nearest the matched-Z radius. No stable root ⇒ a double REAL pole solved for
    //      the Nyquist point alone (|H| is then a ratio of linears in cos w on each side of the
    //      null ⇒ provably monotone ⇒ can never bulge the passband).
    //   2. if nothing wraps yet the cascade misses the analog Nyquist point by > ±0.05 dB (f0
    //      parked near fs/2, poles hot but unwrapped), the top-angle section is refit the same way.
    //   3. a band-stop must never exceed unity (maximal flatness — the analog reference can't).
    //      The result is sampled on a fixed log grid; if it pokes > +0.005 dB anywhere the refit
    //      pole pairs are blended back toward matched-Z (7-step bisection on one global λ — the
    //      stability triangle is convex so every blend is stable; λ=0 is the legacy design, which
    //      never boosts). Only extreme corners blend (f0 ≳ 0.4·fs and wide).
    // When no gate fires the output is BIT-IDENTICAL to the legacy design (verified across a
    // 44.1/48/96 kHz × Q × order × f0 grid), and `sections`≤1 falls back to notch(f0,fs,Q)
    // BIT-FOR-BIT (the legacy single notch / swept anchor). Fixed arrays, noexcept, O(sections)
    // work — RT-safe at design rate. Fills out[0..n) and returns n (== clamped sections).
    // Preconditions (as across matched::): 0 < f0 < fs/2 and Q > 0 — designBand() guarantees them.
    inline int notchCascade (double f0, double fs, double Q, int sections, BiquadCoeffs* out) noexcept
    {
        const int m = sections < 1 ? 1 : (sections > 8 ? 8 : sections);
        if (m == 1) { out[0] = notch (f0, fs, Q); return 1; }

        const double w0  = 2.0 * kPi * f0 / fs;          // digital notch angle (shared zeros at ±w0)
        const double cw  = std::cos (w0);
        const double Om0 = 2.0 * kPi * f0;               // analog centre (rad/s)
        const double BW  = Om0 / Q;                       // analog −3 dB bandwidth (rad/s)
        const double T   = 1.0 / fs;

        // One matched notch biquad from a pair of analog-pole exponentials z1,z2 (a conjugate pair,
        // or two reals for an over-damped centre section): stable poles, zeros pinned at ±w0, gain
        // normalised to unity at DC — the same numerator the single matched::notch uses.
        auto emit = [&] (int idx, const std::complex<double>& z1, const std::complex<double>& z2) noexcept
        {
            const double a1 = -(z1 + z2).real();
            const double a2 =  (z1 * z2).real();
            const double k  = (1.0 + a1 + a2) / (2.0 - 2.0 * cw);
            out[idx] = { k, -2.0 * k * cw, k, a1, a2 };
        };
        // The same section from a pole pair given as (u, v) = (half the pole sum, pole product) —
        // the refit / blend parametrisation (a1 = −2u, a2 = v).
        auto emitUV = [&] (int idx, double u, double v) noexcept
        {
            const double a1 = -2.0 * u, a2 = v;
            const double k  = (1.0 + a1 + a2) / (2.0 - 2.0 * cw);
            out[idx] = { k, -2.0 * k * cw, k, a1, a2 };
        };

        // Pass 1 — legacy matched-Z placement (bit-identical ops), recording each section's analog
        // pole pair and the unwrapped digital angle Im(r)·T that decides aliasing.
        std::complex<double> ra[8], rb[8];
        double ang[8];
        int n = 0;
        for (int i = 1; i <= m; ++i)
        {
            // i-th Butterworth LP prototype pole on the unit circle in the LHP (angles in (π/2, 3π/2)).
            const double theta = kPi * (2.0 * i - 1.0) / (2.0 * m) + kPi * 0.5;
            const std::complex<double> p { std::cos (theta), std::sin (theta) };
            if (p.imag() < -1e-12) continue;             // one representative per conjugate prototype pair

            // LP→BS quadratic  s² + b·s + Ω0² = 0  with b = −BW/p. Take the larger-magnitude root
            // directly (its two terms reinforce — no cancellation), the smaller by Vieta (r1·r2 = Ω0²),
            // so a wide (low-Q) notch keeps full precision in its low stagger pole.
            const std::complex<double> b    = -BW / p;
            const std::complex<double> disc = std::sqrt (b * b - 4.0 * Om0 * Om0);
            const std::complex<double> r1   = (std::abs (-b + disc) >= std::abs (-b - disc))
                                                ? 0.5 * (-b + disc) : 0.5 * (-b - disc);
            const std::complex<double> r2   = (Om0 * Om0) / r1;

            if (std::abs (p.imag()) < 1e-12)
            {
                ra[n] = r1; rb[n] = r2; ang[n] = std::abs (r1.imag()) * T;        // centre pair: Im ≤ Ω0 < π/T,
                emit (n++, std::exp (r1 * T), std::exp (r2 * T));                 // never wraps (real proto pole)
            }
            else
            {
                const std::complex<double> zr1 = std::exp (r1 * T);              // complex proto pole → 2 stagger
                const std::complex<double> zr2 = std::exp (r2 * T);              // sections (each root ⊗ its conjugate)
                ra[n] = r1; rb[n] = std::conj (r1); ang[n] = std::abs (r1.imag()) * T;
                emit (n++, zr1, std::conj (zr1));
                ra[n] = r2; rb[n] = std::conj (r2); ang[n] = std::abs (r2.imag()) * T;
                emit (n++, zr2, std::conj (zr2));
            }
        }

        // Pass 2 — alias gate: sections whose analog pole wrapped under exp(sT) must be refit.
        bool refit[8] = {};
        int nR = 0;
        for (int i = 0; i < n; ++i)
            if (ang[i] > kPi) { refit[i] = true; ++nR; }

        const double Omn = kPi / T;                                              // Nyquist (rad/s)
        auto analogSec = [&] (int i, double Om) noexcept                         // section i's analog magnitude
        {                                                                        // (zeros ±iΩ0, DC-normalised)
            const std::complex<double> iw (0.0, Om);
            const double num = std::abs (Om0 * Om0 - Om * Om);
            const double den = std::abs (iw - ra[i]) * std::abs (iw - rb[i]);
            return std::abs (ra[i] * rb[i]) / (Om0 * Om0) * num / den;
        };
        auto analogCascade = [&] (double Om) noexcept                            // the closed-form reference
        {
            const double x = BW * Om / (Om0 * Om0 - Om * Om);
            return 1.0 / std::sqrt (1.0 + std::pow (std::abs (x), 2.0 * (double) m));
        };

        if (nR == 0)                                                             // nothing wrapped: promote the
        {                                                                        // top-angle section only if the
            double err = analogCascade (Omn);                                    // Nyquist point is off > ±0.05 dB
            for (int i = 0; i < n; ++i) err /= out[i].magnitude (kPi);
            if (err > 1.0058 || err < 0.9943)
            {
                int top = 0;
                for (int i = 1; i < n; ++i) if (ang[i] > ang[top]) top = i;
                refit[top] = true; nR = 1;
            }
            if (nR == 0) return n;                                               // BIT-IDENTICAL legacy output
        }

        // Matched-Z error of the sections we KEEP, at digital frequency w — the refit sections
        // absorb it (split evenly in magnitude) so the CASCADE, not just each section, lands on
        // the analog reference at the match points.
        auto residKept = [&] (double w) noexcept
        {
            double r = 1.0;
            for (int i = 0; i < n; ++i)
                if (! refit[i]) r *= analogSec (i, w / T) / out[i].magnitude (w);
            return r;
        };

        const double wm     = std::sqrt (w0 * kPi);                              // 2nd match point (log-mid of null→Nyquist)
        const double tm     = std::cos (wm);
        const double residN = std::pow (residKept (kPi), 1.0 / (double) nR);
        const double residM = std::pow (residKept (wm),  1.0 / (double) nR);

        double u0[8], v0[8], u1[8], v1[8];                                       // pole pairs: matched-Z … refit
        for (int i = 0; i < n; ++i)
        {
            if (! refit[i]) continue;
            u0[i] = -0.5 * out[i].a1; v0[i] = out[i].a2;

            const double Mn    = analogSec (i, Omn)    * residN;                 // Nyquist target
            const double Mm    = analogSec (i, wm / T) * residM;                 // wm target
            const double alpha = Mn * (1.0 - cw) / (1.0 + cw);

            // |H(π)| = Mn is LINEAR in (u,v): (1−2u+v)(1+c) = Mn(1−c)(1+2u+v) ⇒ u = g(1+v) with
            // g below; |H(wm)| = Mm then reduces to qa·s² + qb·s + qc = 0 in s = 1+v (numerator
            // power N(t) = 4(t−c)², denominator power D(t) = (1−v)² + 4u² − 4u(1+v)t + 4vt²).
            const double g  = (1.0 - alpha) / (2.0 * (1.0 + alpha));
            const double L  = (1.0 - 2.0 * g) * (1.0 - 2.0 * g) * 4.0 * (tm - cw) * (tm - cw)
                              / ((2.0 - 2.0 * cw) * (2.0 - 2.0 * cw));
            const double A_ = 1.0 + 4.0 * g * g - 4.0 * g * tm;
            const double B_ = 4.0 * (tm * tm - 1.0);
            const double qa = L - Mm * Mm * A_, qb = -Mm * Mm * B_, qc = Mm * Mm * B_;

            double s = -1.0; bool ok = false;
            if (std::abs (qa) > 1e-300)
            {
                const double dsc = qb * qb - 4.0 * qa * qc;
                if (dsc >= 0.0)
                {
                    const double sd = std::sqrt (dsc);
                    const double s1 = (-qb + sd) / (2.0 * qa), s2 = (-qb - sd) / (2.0 * qa);
                    const double Rz   = std::exp (ra[i].real() * T);             // matched-Z radius as the
                    const double want = 1.0 + Rz * Rz;                           // root-selection prior
                    const bool ok1 = s1 > 0.0 && s1 < 2.0, ok2 = s2 > 0.0 && s2 < 2.0;
                    if (ok1 && ok2) { s = std::abs (s1 - want) < std::abs (s2 - want) ? s1 : s2; ok = true; }
                    else if (ok1)   { s = s1; ok = true; }
                    else if (ok2)   { s = s2; ok = true; }
                }
            }
            else if (std::abs (qb) > 1e-300) { s = -qc / qb; ok = (s > 0.0 && s < 2.0); }

            bool solved = false;
            if (ok)
            {
                const double v = s - 1.0, u = g * s;
                const BiquadCoeffs c { 1.0, 0.0, 1.0, -2.0 * u, v };             // denominator-only stability probe
                if (c.isStable()) { u1[i] = u; v1[i] = v; solved = true; }
            }
            if (! solved)                                                        // 1-pt monotone fallback: double
            {                                                                    // real pole, Nyquist target exact
                const double sa = std::sqrt (std::max (0.0, alpha));
                const double Rp = (alpha >= 1.0) ? (sa - 1.0) / (sa + 1.0) : -(1.0 - sa) / (1.0 + sa);
                u1[i] = -Rp; v1[i] = Rp * Rp;
            }
        }

        // Pass 3 — maximal-flatness guard. Blend λ: 0 = legacy matched-Z pole pairs, 1 = full
        // refit. Sample the cascade on a fixed log grid; on any poke above unity bisect λ down.
        auto apply = [&] (double lam) noexcept
        {
            for (int i = 0; i < n; ++i)
                if (refit[i]) emitUV (i, u0[i] + lam * (u1[i] - u0[i]), v0[i] + lam * (v1[i] - v0[i]));
        };
        auto overUnity = [&] () noexcept
        {
            for (int j = 0; j <= 47; ++j)
            {
                const double w = kPi * std::pow (0.002, 1.0 - j / 47.0);
                double mag = 1.0;
                for (int i = 0; i < n; ++i) mag *= out[i].magnitude (w);
                if (mag > 1.0006) return true;                                   // > ~+0.005 dB
            }
            return false;
        };

        apply (1.0);
        if (overUnity())
        {
            double lo = 0.0, hi = 1.0;                                           // λ=0 == legacy: never boosts
            for (int it = 0; it < 7; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                apply (mid);
                if (overUnity()) hi = mid; else lo = mid;
            }
            apply (lo);
        }
        return n;
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

        // Feasibility guard: near Nyquist a low-Q, high-gain shelf makes the 3-point (DC/corner/Nyquist)
        // fit infeasible (W²+B2 ≤ 0). Forcing it puts the numerator zeros on/outside the unit circle → a
        // spurious deep passband notch in a boost, and via the reciprocal-cut path a near-unit-circle
        // pole → a huge peak. The infeasible corner is always LOW Q (no real resonance to lose), so fall
        // back to the smooth non-resonant matched Butterworth shelf: exact plateaus, monotone, no zeros
        // near the unit circle. Feasible designs (incl. every legitimate high-Q resonance) are untouched.
        if (! detail::shelfFitFeasible (B0, B1, B2))
            return high ? highShelf (f0, fs, gainLin) : lowShelf (f0, fs, gainLin);

        detail::solveNumerator (B0, B1, B2, c);

        // Post-design envelope guard for the FEASIBLE-but-pathological corner. A feasible fit can still
        // throw a spurious NOTCH into a boost (a +18 dB shelf dipping −23 dB mid-band → its reciprocal
        // CUT peaks +23 dB) or a spurious PEAK (a +6 dB low shelf reading +15 dB). Bound the realised
        // |H| by the INTENDED analog resonant-shelf envelope (poles (w0,Q), zeros (w0√A,Q)): a genuine
        // high-Q resonance — deep zero-driven dip AND high pole peak — stays INSIDE that envelope, while
        // a spurious excursion escapes it. Both computed in closed form (detail::ratioExtent).
        const double kq = 1.0 / (Q * Q) - 2.0;                          // analog |H|² = N(y)/D(y), y=(f/f0)²  (A already = gainLin)
        const auto an = detail::ratioExtent (high ? A * A : 1.0, A * kq, high ? 1.0 : A * A,   // N(y)
                                             1.0, kq, 1.0, 0.0, std::numeric_limits<double>::infinity());   // D(y), y∈[0,∞)
        const auto dg = detail::ratioExtent (4.0 * c.b0 * c.b2, 2.0 * c.b1 * (c.b0 + c.b2),    // realised |H|²
                                             (c.b0 - c.b2) * (c.b0 - c.b2) + c.b1 * c.b1,
                                             4.0 * c.a2, 2.0 * c.a1 * (1.0 + c.a2),
                                             (1.0 - c.a2) * (1.0 - c.a2) + c.a1 * c.a1, -1.0, 1.0);          // t=cos w∈[−1,1]
        const double m2 = 2.0;                                          // (10^(3/20))² ≈ +3 dB in power
        if (dg.hi > an.hi * m2 || dg.lo < an.lo / m2)                    // escapes the analog envelope by > 3 dB
            return high ? highShelf (f0, fs, gainLin) : lowShelf (f0, fs, gainLin);

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
