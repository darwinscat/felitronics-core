// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/FlushToZero.h>

#include <algorithm>
#include <cmath>

namespace felitronics::dynamics
{

//==============================================================================
// felitronics::dynamics::TransientShaper — attack/sustain shaping (SPL-Transient-Designer style). A FAST
// and a SLOW envelope on the linked sidechain, BOTH Peak (same detector type — different crest factors,
// e.g. Peak vs Rms, would read a steady sine as norm≈0.17 and falsely shape it). Their contrast tells the
// ATTACK phase (fast above slow, the onset) from the SUSTAIN phase (slow above fast, the decay):
//     norm = (fast − slow) / (fast + slow + ε)   ∈ [−1, 1]
//     shaped = deadzone(norm, threshold)         (a steady tone's micro-crest stays inside the deadzone)
//     gainDb = attackDb·max(shaped,0) + sustainDb·max(−shaped,0)     (two independent ± controls)
// A boost/cut on the attack sharpens/softens hits; on the sustain it lengthens/tightens tails. The gain is
// one-pole smoothed (no zipper). ZERO latency: the gain is computed from, and applied to, the SAME sample —
// the detector reacts within its 0.3 ms attack. (A proper anticipatory lookahead — gain pre-ramped over a
// delay window — is a future option; the naive "detect-now, apply-to-delayed" leads the audio and mis-times
// the boost, so it's intentionally not shipped.)
//
// This is GAIN MODULATION, not waveshaping → no oversampling. ONE linked gain drives all channels
// (image-preserving). Broadband; a per-band version is multiband::MultibandProcessor<TransientShaper>.
// RT-safe: prepare() allocates nothing on the hot path; process() is alloc/lock/throw-free, in place.
struct TransientShaperParams
{
    LinkMode link = LinkMode::Max;
    double attackDb     = 0.0;     // −24..24 — attack-phase boost(+) / cut(−)
    double sustainDb    = 0.0;     // −24..24 — sustain-phase boost(+) / cut(−)
    double fastAttackMs = 0.3,  fastReleaseMs = 20.0;
    double slowAttackMs = 15.0, slowReleaseMs = 150.0;
    double threshold    = 0.2;     // 0..0.9 — deadzone on |norm| so a steady tone's micro-crest isn't shaped
    double gainSmoothMs = 1.0;     // gain de-zipper time
    double mix          = 1.0;     // 0 dry … 1 fully shaped
};

class TransientShaper
{
public:
    void prepare (double sampleRate, int /*maxBlock*/ = 0, int /*maxChannels*/ = 2) noexcept
    {
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        fast_.prepare (fs_); fast_.setDetector (Detector::Peak);
        slow_.prepare (fs_); slow_.setDetector (Detector::Peak);   // SAME type as fast → steady tone = norm 0
        apply (params_);
        reset();
    }

    void reset() noexcept { fast_.reset(); slow_.reset(); gainSm_ = 1.0f; }

    void setParams (const TransientShaperParams& p) noexcept { params_ = p; apply (p); }
    static constexpr int latencySamples() noexcept { return 0; }

    void process (float* const* channels, int numChannels, int n) noexcept
    {
        const int nc = numChannels < core::kMaxChannels ? numChannels : core::kMaxChannels;
        if (nc <= 0) return;
        for (int i = 0; i < n; ++i)
        {
            const float linked = linkAmplitude (params_.link, channels, nc, i);
            const float fe = fast_.process (linked);
            const float se = slow_.process (linked);
            const float norm = (fe - se) / (fe + se + 1.0e-9f);                 // ∈ [−1, 1]
            // deadzone: a steady tone has a small fast/slow ripple — only a real transient clears `threshold`.
            const float a = std::fabs (norm);
            const float t = a <= threshold_ ? 0.0f : (a - threshold_) / (1.0f - threshold_);   // rescale [thr,1]→[0,1]
            const float shaped = norm < 0.0f ? -t : t;
            float gdb = attackDb_ * std::max (shaped, 0.0f) + sustainDb_ * std::max (-shaped, 0.0f);
            gdb = std::clamp (gdb, -24.0f, 24.0f);
            const float g = (float) core::dbToGain ((double) gdb);
            gainSm_ = g + smoothCoeff_ * (gainSm_ - g);                         // one-pole de-zipper

            const float m = mix_ * (gainSm_ - 1.0f);                            // 0 dry … applied at mix 1
            for (int c = 0; c < nc; ++c) channels[c][i] += channels[c][i] * m;  // = x*(1 + mix*(gain-1))
        }
        fast_.flushDenormals(); slow_.flushDenormals();
    }

private:
    static double finite (double v, double fallback) noexcept { return std::isfinite (v) ? v : fallback; }

    void apply (const TransientShaperParams& p) noexcept
    {
        fast_.setTimes (finite (p.fastAttackMs, 0.3), finite (p.fastReleaseMs, 20.0));
        slow_.setTimes (finite (p.slowAttackMs, 15.0), finite (p.slowReleaseMs, 150.0));
        attackDb_  = (float) std::clamp (finite (p.attackDb,  0.0), -36.0, 36.0);
        sustainDb_ = (float) std::clamp (finite (p.sustainDb, 0.0), -36.0, 36.0);
        threshold_ = (float) std::clamp (finite (p.threshold, 0.2), 0.0, 0.9);
        mix_       = (float) std::clamp (finite (p.mix, 1.0), 0.0, 1.0);
        const double tms = finite (p.gainSmoothMs, 1.0) * 0.001;
        smoothCoeff_ = (tms <= 0.0 || fs_ <= 0.0) ? 0.0f : (float) std::exp (-1.0 / (tms * fs_));
    }

    double fs_ = 48000.0;
    TransientShaperParams params_;
    EnvelopeFollower fast_, slow_;
    float attackDb_ = 0.0f, sustainDb_ = 0.0f, threshold_ = 0.2f, mix_ = 1.0f, smoothCoeff_ = 0.0f, gainSm_ = 1.0f;
};

} // namespace felitronics::dynamics
