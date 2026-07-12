// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/saturation/WaveShaper.h>

//==============================================================================
// felitronics::poweramp::TubeStage — the PURE tube transfer math, factored out of the trunk
// (PowerAmpStage) so the static curve can be unit-tested DIRECTLY (exact PP odd-symmetry, exact
// drive-comp unity, slope, boundedness) with no oversampling, no DFT. Header-only, std +
// felitronics::saturation only (embedded-safe). The shipping stage drives this same code, so the
// tests exercise the real kernel, not a parallel reimplementation.
//
//   SE (single-ended class A): y = g(u)            → asymmetric, even-harmonic rich
//   PP (push-pull class AB):   y = g(u+Vb) − g(−u+Vb)  → ODD by construction (even cancel)
// g = felitronics::saturation::WaveShaper(Asym) = (tanh(k(t+b)) − tanh(kb)) peak-normalised.
// `topoMix` blends PP(0) ↔ SE(1). The per-half bias mismatch `evenLeak` is the only thing that
// breaks PP's exact even-cancellation — by construction, so tests can assert it.
//
// Lifted verbatim from OrbitCab's tube power amp (cab::poweramp::TubeStage, TubeKernel.h). The
// per-tube voicing DATA (OrbitCab's kTubeVoicings presets) is product data and stays with the
// product — this kernel only consumes the plain numbers configure() is given.
//==============================================================================
namespace felitronics::poweramp
{

// Stateless composite tube transfer. configure() once per block from (smoothed) coeffs; at()
// per (oversampled) sample. slopeAtZero() is the composite small-signal slope incl. the pre-gain.
class TubeStage
{
public:
    void configure (float k, float biasSE, float vbPP, float evenLeak, float topoMix) noexcept
    {
        using Shape = felitronics::saturation::WaveShaper::Shape;
        vb_ = vbPP; topo_ = topoMix;
        se_.setShape  (Shape::Asym); se_.setDrive  (k); se_.setBias  (biasSE);
        pos_.setShape (Shape::Asym); pos_.setDrive (k); pos_.setBias ( evenLeak);
        neg_.setShape (Shape::Asym); neg_.setDrive (k); neg_.setBias (-evenLeak);
    }

    float at (float u) const noexcept
    {
        const float se = se_.processSample (u);
        const float pp = pos_.processSample (u + vb_) - neg_.processSample (-u + vb_);
        return (1.0f - topo_) * pp + topo_ * se;
    }

    // Per-sample PP-bias overload (dynamic bias-shift): `vbDelta` drifts the push-pull operating
    // point sample-by-sample (toward class-B under sag → crossover bloom), with NO per-block reconfigure
    // (which would break block-size determinism). SE (class A) has no crossover, so it stays unmodulated.
    // vbDelta == 0 is bit-identical to at(u). The shaper is memoryless, so re-evaluating it is free of state.
    float at (float u, float vbDelta) const noexcept
    {
        const float se = se_.processSample (u);
        const float vb = vb_ + vbDelta;
        const float pp = pos_.processSample (u + vb) - neg_.processSample (-u + vb);
        return (1.0f - topo_) * pp + topo_ * se;
    }

    // d/dx at(G·x) at x=0 — numeric central difference (topology-agnostic: covers PP's 2-eval
    // differential and the pre-gain G in one go, which felitronics' per-curve slopeAtZero cannot).
    float slopeAtZero (float G) const noexcept
    {
        const float d = 1.0e-4f;
        return (at (G * d) - at (-G * d)) / (2.0f * d);
    }

private:
    felitronics::saturation::WaveShaper se_, pos_, neg_;
    float vb_ = 0.30f, topo_ = 0.0f;
};

} // namespace felitronics::poweramp
