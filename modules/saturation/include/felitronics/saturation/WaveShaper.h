// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>

namespace felitronics::saturation
{

//==============================================================================
// felitronics::saturation::WaveShaper — a stateless soft-saturation transfer curve, PEAK-NORMALISED (all but
// Transformer, which is slope-normalised — see its note) so |x| <= 1 maps to |y| <= 1: the curve rounds the
// top + adds harmonics WITHOUT changing the full-scale level (the right shape for a mastering "glue"
// saturator that sits before the make-loud gain).
//
// `drive` (k) sets how much curve: k → 0 is ~linear (no effect), larger k = more harmonics. Curves:
// Tanh (odd, smooth/dark — the safe default), Atan (odd, a touch brighter/harder
// knee), Cubic (mostly 3rd, a cheap soft clipper), Asym (tube/triode — `bias` adds EVEN harmonics),
// Tube (a fixed tube-like bias: mostly 2nd, no knob of its own), Transistor (odd, cleaner than Tanh below
// the knee, harder past it), Transformer (odd; the curve below is only its static core — see the note):
//   Tanh        y = tanh(kx)/tanh(k)
//   Atan        y = atan(kx)/atan(k)
//   Cubic       y = h(kx)/h(k),  h(u)=1.5u-0.5u³ (|u|<1), sign(u) else
//   Asym        r(x)=tanh(k(x+b))-tanh(kb),  y = r(x)/max(|r(1)|,|r(-1)|)
//   Tube        r(u)=q·t/(1+c·t), t=tanh(u), u=kx, c=tanh(0.3), q=1-c²,  y = r(kx)/|r(-k)|
//   Transistor  s(u)=u/(1+u⁴)^(1/4),  y = s(kx)/s(k)
//   Transformer y = tanh(kx)/k                        (SLOPE-normalised, not peak-normalised)
// Pure function of the input (no per-sample state) → this is the kernel an oversampled Saturator wraps
// (run it at N× so the new harmonics stay below the base Nyquist). Asym's bias makes y(0)=0 but NOT a
// zero-mean output for music, so the Saturator follows it with a DC blocker; Tube gets the same blocker.
//
// TUBE is tanh(u+0.3)-tanh(0.3) — the addition theorem tanh(a+c) = (ta+tc)/(1+ta·tc) turns the difference
// into q·t/(1+c·t) — evaluated WITHOUT the difference: in float the difference cancels to exactly 0 at the
// drive floor (k = 1e-4, x = 1e-4) and is off by 5.1e-4 at driveDb 0, where this form stays within 1e-6 of
// a double evaluation of the same curve at every drive. `bias` does not reach it (Params::bias stays
// Asym's). The curve rises faster than it falls, so its peak is always the NEGATIVE side: at driveDb 6 and
// a -3 dBFS sine the positive peak sits 3.1 dB below the negative one, and H2 leads H3 by 9.3 dB. Its
// normalised small-signal slope k·q/|r(-k)| is BELOW 1 for driveDb in (0, ~5.7) — 0.936 at its minimum
// near 3.2 dB — so autoComp = 1 lifts the negative peak by up to +0.58 dB (a full-scale sine with the
// Saturator's DC blocker off, or before it settles; settled, the blocker recentres the wave and the peak
// is +0.17 dB). Asym with bias > 1/3 already does the same.
//
// TRANSISTOR has no cubic term (s = u - u⁵/4 + ...), so at low drive it is MORE transparent than Tanh —
// H3 -35.3 dBc at driveDb 6 and -3 dBFS against Tanh's -28.7; its H3 stays under Tanh's up to ~9 dB of
// drive at -3 dBFS (~7 dB at full scale) — and past its knee it is harder (more high odd harmonics). |u|
// is clamped to 64 before u⁴ — s(64) is exactly 1 in float, and without the clamp u⁴ overflows near
// |u| = 4.3e9 and a full-scale sample maps to 0 (driveDb is not clamped). In float s is not monotone (it
// steps back by up to 3 ulp), so its normalised peak over [-1, 1] is 1 to within one ulp.
//
// TRANSFORMER is the one shape that is not a pure function of the input: the Saturator runs it as a model with a
// per-channel flux state (see Saturator.h), and what this class holds is that model's static core sat(u) =
// tanh(ku)/k. It is normalised by its SLOPE, not its peak — norm = 1/k, so slopeAtZero() is 1 at every drive
// (the literal 1.0f, which makes the Saturator's drive-compensation exactly 1.0f) and a quiet signal passes at
// unity gain. |y| <= |x| still holds, so |x| <= 1 still maps to |y| <= 1, but the full-scale output is not 1:
// y(1) = tanh(k)/k < 1. processSample() here is that static curve alone, without the flux.
class WaveShaper
{
public:
    // Append only: a consumer pins these integers in an ABI. Never renumber.
    enum class Shape { Tanh, Atan, Cubic, Asym, Tube, Transistor, Transformer };

    // Tube's fixed bias c = tanh(0.3) and q = 1 - c², as the nearest floats (std::tanh is not constexpr).
    static constexpr float kTubeC = 0.2913126124515909f;
    static constexpr float kTubeQ = 0.9151369618266292f;

    void setShape (Shape s) noexcept { shape_ = s; updateNorm(); }
    void setDrive (float k) noexcept { drive_ = k > 1.0e-4f ? k : 1.0e-4f; updateNorm(); }   // k > 0
    void setBias  (float b) noexcept { bias_  = std::clamp (b, -0.95f, 0.95f); updateNorm(); } // Asym only

    float drive() const noexcept { return drive_; }

    // THE CURVE'S FOUR COEFFICIENTS, and the curve evaluated at explicit ones. `processSample()` reads the
    // members; `shapeAt<S>()` is the SAME arithmetic, operand for operand, with the coefficients handed in —
    // for a caller that moves them per sample (the Saturator's parameter glide interpolates them between two
    // designed sets). At the coefficients `coeffs()` returns it is `processSample()` exactly. Additive: nothing
    // above changes, and a consumer that never calls these gets the curve it always had.
    struct Coeffs { float drive = 1.0f, bias = 0.0f, biasTanh = 0.0f, norm = 1.0f; };
    Coeffs coeffs() const noexcept { return { drive_, bias_, biasTanh_, norm_ }; }

    template <Shape S>
    static float shapeAt (const Coeffs& c, float x) noexcept
    {
        if constexpr (S == Shape::Tanh)       return std::tanh (c.drive * x) * c.norm;
        else if constexpr (S == Shape::Atan)  return std::atan (c.drive * x) * c.norm;
        else if constexpr (S == Shape::Cubic) return cubicClip (c.drive * x) * c.norm;
        else if constexpr (S == Shape::Asym)  return (std::tanh (c.drive * (x + c.bias)) - c.biasTanh) * c.norm;
        else if constexpr (S == Shape::Tube)  return tubeRaw (c.drive * x) * c.norm;
        else if constexpr (S == Shape::Transistor) return transistorRaw (c.drive * x) * c.norm;
        else if constexpr (S == Shape::Transformer) return std::tanh (c.drive * x) * c.norm;
        else
        {
            // A new Shape gets its own branch above — a bare `else` would play it as the previous formula.
            static_assert (kUnhandled<S>, "WaveShaper::shapeAt: a Shape without its own branch");
            return x;
        }
    }

    // d y / d x at x = 0 — the small-signal gain. The Saturator uses it for drive-compensation
    // (auto-gain = slopeAtZero^(-amount)) so turning up drive doesn't change the low-level loudness.
    float slopeAtZero() const noexcept
    {
        switch (shape_)
        {
            case Shape::Tanh:  return norm_ * drive_;                                  // raw'(0) = 1
            case Shape::Atan:  return norm_ * drive_;                                  // raw'(0) = 1
            case Shape::Cubic: return norm_ * 1.5f * drive_;                           // raw'(0) = 1.5
            case Shape::Asym:  return norm_ * drive_ * (1.0f - biasTanh_ * biasTanh_); // sech²(kb)
            case Shape::Tube:  return norm_ * drive_ * kTubeQ;                         // raw'(0) = q
            case Shape::Transistor: return norm_ * drive_;                             // raw'(0) = 1
            case Shape::Transformer: return 1.0f;       // slope-normalised: exactly 1, not norm_ * drive_ rounded
        }
        return 1.0f;
    }

    float processSample (float x) const noexcept
    {
        switch (shape_)
        {
            case Shape::Tanh:  return std::tanh (drive_ * x) * norm_;
            case Shape::Atan:  return std::atan (drive_ * x) * norm_;
            case Shape::Cubic: return cubicClip (drive_ * x) * norm_;
            case Shape::Asym:  return (std::tanh (drive_ * (x + bias_)) - biasTanh_) * norm_;
            case Shape::Tube:  return tubeRaw (drive_ * x) * norm_;
            case Shape::Transistor: return transistorRaw (drive_ * x) * norm_;
            case Shape::Transformer: return std::tanh (drive_ * x) * norm_;
        }
        return x;
    }

private:
    // Soft cubic clipper: slope 1.5 at 0, reaches ±1 with zero slope at u=±1, flat beyond.
    static float cubicClip (float u) noexcept
    {
        if (u >=  1.0f) return  1.0f;
        if (u <= -1.0f) return -1.0f;
        return 1.5f * u - 0.5f * u * u * u;
    }

    // Tube's raw curve q·t/(1+c·t), t = tanh(u). The product and the sum in SEPARATE statements, so
    // `-ffp-contract=on` fuses neither and the value is the same with and without an FMA unit.
    static float tubeRaw (float u) noexcept
    {
        const float t  = std::tanh (u);
        const float ct = kTubeC * t;
        const float qt = kTubeQ * t;
        return qt / (1.0f + ct);
    }

    // Transistor's raw curve u/(1+u⁴)^(1/4) on |u| clamped to 64, the sign restored. Separate statements as
    // above: 1 + a⁴ is never an FMA.
    static float transistorRaw (float u) noexcept
    {
        const float a  = std::min (std::fabs (u), 64.0f);
        const float a2 = a * a;
        const float a4 = a2 * a2;
        const float d  = 1.0f + a4;
        return std::copysign (a / std::sqrt (std::sqrt (d)), u);
    }

    template <Shape> static constexpr bool kUnhandled = false;

    // Precompute the peak normaliser (raw curve value at the extreme input) so processSample stays cheap.
    void updateNorm() noexcept
    {
        biasTanh_ = std::tanh (drive_ * bias_);
        float raw = 1.0f;
        switch (shape_)
        {
            case Shape::Tanh:  raw = std::tanh (drive_); break;
            case Shape::Atan:  raw = std::atan (drive_); break;
            case Shape::Cubic: raw = cubicClip (drive_); break;
            case Shape::Asym:
            {
                const float rPos = std::tanh (drive_ * ( 1.0f + bias_)) - biasTanh_;
                const float rNeg = std::tanh (drive_ * (-1.0f + bias_)) - biasTanh_;
                raw = std::max (std::fabs (rPos), std::fabs (rNeg));
                break;
            }
            // The curve at x = -1, by the very arithmetic processSample() runs: its peak for any k > 0.
            case Shape::Tube:       raw = std::fabs (tubeRaw (-drive_)); break;
            case Shape::Transistor: raw = transistorRaw (drive_); break;
            // Not a peak: the slope of tanh(k·u) at 0, so the normalised curve's slope there is 1.
            case Shape::Transformer: raw = drive_; break;
        }
        norm_ = (raw > 1.0e-12f) ? 1.0f / raw : 1.0f;
    }

    Shape shape_   = Shape::Tanh;
    float drive_   = 1.0f;
    float bias_    = 0.0f;
    float biasTanh_ = 0.0f;
    float norm_    = 1.0f / 0.7615941559557649f;   // 1/tanh(1)
};

} // namespace felitronics::saturation
