// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/FlushToZero.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/EnvelopeFollower.h>

#include <cmath>

namespace felitronics::dynamics
{

// The DETECTOR half of a dynamics processor's parameters, on its own so that an OFFLINE analysis pass
// can be handed the SAME object the running processor holds. `CompressorParams` derives from this, so
// `detector.setParams (compressorParams)` is a base-class binding — there is no hand-written field
// mapping between the two, and therefore nothing that can drift when a field is added.
struct DetectorParams
{
    Detector detector    = Detector::Peak;   // Peak = instant |x| · Rms = a one-pole average of the power
    LinkMode link        = LinkMode::Max;    // how the key's channels collapse to one level
    double   rmsWindowMs = 5.0;              // only used when detector == Rms
};

//==============================================================================
// felitronics::dynamics::LinkedDetector — the whole detector path of a compressor as ONE object:
//
//     key frame → per-sample gate → channel link → envelope → linked level (linear amplitude)
//
// It exists because "the offline analysis ran the same detector as the compressor" has to be a
// STRUCTURAL fact rather than a promise kept by hand. Everything downstream of this (the static curve,
// the ballistics on the gain reduction) is already stateless or already a shared primitive, so this was
// the last piece that a second implementation could get subtly wrong: the compressor used to hand-roll
// its own one-pole on the power, and reproducing it offline meant copying five lines and hoping.
// "The same envelope" is a WITHIN-A-BUILD identity, and that is the useful reading: same object, same
// parameters, same sample rate, same flush cadence, same binary. Across toolchains the desktop tier
// declares `-ffp-contract=on` and libm is not bit-portable, so an offline pass compiled elsewhere is
// equal to the running one to rounding, not to the last bit.
//
// TIMES ARE SYMMETRIC, and that is a decision, not an omission. An RMS detector is an AVERAGE; an
// asymmetric attack/release on the power biases the tracked level of a steady tone (EnvelopeFollower's
// own note), so the window is applied to both. `Peak` is INSTANT for the same reason: in this topology
// the ballistics live on the gain reduction, downstream of the static curve, precisely so that the knee
// and the ratio cannot warp the attack and release.
//
// PRODUCT-NEUTRAL: no EQ. A sidechain filter is the product's (ADR §4 — "the core gives primitives;
// products compose them"), and it belongs UPSTREAM of this object: filter the key, then feed it here,
// and the very same filtered key goes to the offline pass.
//
// RT-safe: no allocation, lock, IO or throw; the only state is one float.
class LinkedDetector
{
public:
    void prepare (double sampleRate) noexcept
    {
        fs_ = (std::isfinite (sampleRate) && sampleRate > 0.0) ? sampleRate : 48000.0;
        env_.prepare (fs_);
        apply();
        reset();
    }

    void setParams (const DetectorParams& p) noexcept { p_ = p; apply(); }

    const DetectorParams& params()     const noexcept { return p_; }
    double                sampleRate() const noexcept { return fs_; }

    void reset() noexcept { env_.reset(); }

    // One frame of the KEY at `sampleIndex` → the linked detector level (linear amplitude).
    // `key[0 .. numKeyChannels-1]` are read; `numKeyChannels` is unrelated to the number of channels
    // the audio has, so a mono key on a stereo programme is one channel here and two there.
    // PRECONDITION: `key` is non-null and `numKeyChannels >= 1`. This is the primitive, so it does not
    // pay for a branch on either — `Compressor::process` is where "no key" is resolved into a source.
    inline float process (const float* const* key, int numKeyChannels, int sampleIndex) noexcept
    {
        return env_.process (linkAmplitudeGated (p_.link, key, numKeyChannels, sampleIndex));
    }

    // Same, for a caller that already holds one mono key sample. Gated identically, so the two entry
    // points cannot disagree about what the detector sees.
    inline float processSample (float keySample) noexcept { return env_.process (std::fabs (detectorGate (keySample))); }

    // The current level (linear amplitude) without advancing anything — for metering.
    float level() const noexcept { return env_.envelope(); }

    // Law 8: call once per block. The state is a one-pole that decays toward zero on silence.
    void flushDenormals() noexcept { env_.flushDenormals(); }

private:
    void apply() noexcept
    {
        env_.setDetector (p_.detector);                                     // converts the state, never reinterprets it
        const double t = (p_.detector == Detector::Rms) ? p_.rmsWindowMs : 0.0;
        env_.setTimes (t, t);                                               // symmetric — see the note above
    }

    EnvelopeFollower env_;
    DetectorParams   p_ {};
    double           fs_ = 48000.0;
};

} // namespace felitronics::dynamics
