// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>

namespace felitronics::saturation
{

//==============================================================================
// felitronics::saturation::WaveShaper — a stateless soft-saturation transfer curve, PEAK-NORMALISED so
// |x| <= 1 maps to |y| <= 1: the curve rounds the top + adds harmonics WITHOUT changing the full-scale
// level (the right shape for a mastering "glue" saturator that sits before the make-loud gain).
//
// `drive` (k) sets how much curve: k → 0 is ~linear (no effect), larger k = more harmonics. Curves
// (DSP-council pick): Tanh (odd, smooth/dark — the safe default), Atan (odd, a touch brighter/harder
// knee), Cubic (mostly 3rd, a cheap soft clipper), Asym (tube/triode — `bias` adds EVEN harmonics):
//   Tanh   y = tanh(kx)/tanh(k)
//   Atan   y = atan(kx)/atan(k)
//   Cubic  y = h(kx)/h(k),  h(u)=1.5u-0.5u³ (|u|<1), sign(u) else
//   Asym   r(x)=tanh(k(x+b))-tanh(kb),  y = r(x)/max(|r(1)|,|r(-1)|)
// Pure function of the input (no per-sample state) → this is the kernel an oversampled Saturator wraps
// (run it at N× so the new harmonics stay below the base Nyquist). Asym's bias makes y(0)=0 but NOT a
// zero-mean output for music, so the Saturator follows it with a DC blocker.
class WaveShaper
{
public:
    enum class Shape { Tanh, Atan, Cubic, Asym };

    void setShape (Shape s) noexcept { shape_ = s; updateNorm(); }
    void setDrive (float k) noexcept { drive_ = k > 1.0e-4f ? k : 1.0e-4f; updateNorm(); }   // k > 0
    void setBias  (float b) noexcept { bias_  = std::clamp (b, -0.95f, 0.95f); updateNorm(); } // Asym only

    float drive() const noexcept { return drive_; }

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
