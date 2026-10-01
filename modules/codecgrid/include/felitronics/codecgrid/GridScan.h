// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/codecgrid/GridCurve.h>
#include <felitronics/codecgrid/Mdct.h>
#include <felitronics/codecgrid/Mp3Hybrid.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

//==============================================================================
// felitronics::codecgrid::GridScan — ONE codec transform run at EVERY frame offset of one stretch of stereo PCM,
// in bounded steps: the curve of that transform over that stretch (GridCurve.h says what a curve is).
//
// THE FOUR TRANSFORMS, and what an offset means in each:
//   Mp3      the MPEG-1 Layer III hybrid filterbank, long blocks; 576 offsets, 576 lines read in 8 groups of 72.
//   AacSine  an MDCT of 2048 samples, hop 1024, the sine window;                     1024 offsets, 1024 lines.
//   AacKbd   the same with the Kaiser-Bessel-derived window, alpha 4;                1024 offsets, 1024 lines.
//   Celt     an MDCT of 1920 samples, hop 960, CELT's low-overlap window, over the stretch pre-emphasised as the
//            encoder does it (x[n] - 0.85 x[n-1]); 960 offsets, the first 800 lines read — CELT codes to 20 kHz.
// The stretch is taken AT THE CODEC'S RATE: a caller that suspects a grid at another rate converts first
// (BackResampler). Long blocks only everywhere: a frame coded with short blocks does not line up and reads as a
// frame without zeros.
//
// THE REFERENCE of every coefficient — its own rms over the stretch — is taken once, at offset 0, and every other
// offset is read against it: the level of a coefficient does not depend on where the frames are cut, only the
// zeros do.
//
// THE NUMBER OF FRAMES is the same at every offset of the MDCT transforms, floor ((n - frame) / hop), so that no
// offset reads a frame another cannot; for Mp3 it follows the polyphase phase and the granule start and differs
// by one between offsets, and the mean does not care.
//
// STEPS. begin() names the stretch and the curve to fill; step (k) does at most k offsets and says whether the
// curve is complete. Nothing is allocated after prepare(). The stretch must stay where it is until the last
// step, and until the last zeroProfile() call made for it.
//==============================================================================
namespace felitronics::codecgrid
{

enum class Transform : std::uint8_t { Mp3, AacSine, AacKbd, Celt };

inline constexpr int kTransforms = 4;
inline constexpr double kCeltPreEmphasis = 0.85000610;       // CELT's pre-emphasis coefficient at 48 kHz
inline constexpr int kAacFrame = 2048, kCeltFrame = 1920, kCeltOverlap = 120, kCeltLines = 800;

// The frame hop of a transform: the number of offsets it has, and the modulus of its grid phase.
constexpr int hopOf (Transform t) noexcept
{
    return t == Transform::Mp3 ? Mp3Hybrid::kLines : t == Transform::Celt ? kCeltFrame / 2 : kAacFrame / 2;
}

// The lines read per frame.
constexpr int linesOf (Transform t) noexcept
{
    return t == Transform::Mp3 ? Mp3Hybrid::kLines : t == Transform::Celt ? kCeltLines : kAacFrame / 2;
}

// The shortest stretch a transform can read at every one of its offsets with at least `frames` frames.
constexpr int minSamplesOf (Transform t, int frames) noexcept
{
    if (t == Transform::Mp3)
        return 31 + Mp3Hybrid::kWindow + 32 * (17 + 18 * (frames + 1) - 1);       // phase 31, granule start 17
    const int frame = t == Transform::Celt ? kCeltFrame : kAacFrame;
    return frame + frames * (frame / 2);
}

template <ComplexFftBackend Fft = MixedRadixFft>
class GridScan
{
public:
    static constexpr int kMaxLines = kAacFrame / 2;

    // What prepare() asks the heap for, for stretches of up to `maxSamples` samples.
    static std::uint64_t bytesFor (int maxSamples) noexcept
    {
        if (maxSamples <= 0) return 0;
        const std::uint64_t n = (std::uint64_t) maxSamples;
        const std::uint64_t sub = 2u * (n / 32u + 1u) * 32u;                         // the subband samples of both channels
        return sizeof (float) * (sub + 2u * n + 2u * (std::uint64_t) kMaxLines + (std::uint64_t) LevelAccumulator::floatsFor (kMaxLines))
             + sizeof (double) * ((std::uint64_t) LevelAccumulator::doublesFor (kMaxLines) + (std::uint64_t) kAacFrame + (std::uint64_t) kAacFrame / 2u + 1u)
             + 2u * (std::uint64_t) Mdct<Fft>::bytesFor (kAacFrame) + (std::uint64_t) Mdct<Fft>::bytesFor (kCeltFrame);
    }

    [[nodiscard]] bool prepare (int maxSamples)
    {
        prepared_ = false;
        active_ = false;
        if (maxSamples <= 0) return false;
        std::vector<double> w ((std::size_t) kAacFrame), scratch ((std::size_t) kAacFrame / 2 + 1);
        window::sine (w.data(), kAacFrame);
        if (! sine_.prepare (kAacFrame, w.data())) return false;
        window::kaiserBesselDerived (w.data(), kAacFrame, 4.0, scratch.data());
        if (! kbd_.prepare (kAacFrame, w.data())) return false;
        window::celt (w.data(), kCeltFrame, kCeltOverlap);
        if (! celt_.prepare (kCeltFrame, w.data())) return false;

        const std::size_t n = (std::size_t) maxSamples;
        sub_[0].assign ((n / 32u + 1u) * 32u, 0.0f);
        sub_[1].assign ((n / 32u + 1u) * 32u, 0.0f);
        pre_[0].assign (n, 0.0f);
        pre_[1].assign (n, 0.0f);
        lines_[0].assign ((std::size_t) kMaxLines, 0.0f);
        lines_[1].assign ((std::size_t) kMaxLines, 0.0f);
        floats_.assign (LevelAccumulator::floatsFor (kMaxLines), 0.0f);
        doubles_.assign (LevelAccumulator::doublesFor (kMaxLines), 0.0);
        maxSamples_ = maxSamples;
        prepared_ = true;
        return true;
    }

    bool prepared() const noexcept { return prepared_; }

    // Begin the scan of left[0 .. n), right[0 .. n) with one transform. `curve` takes hopOf (t) * kCells floats
    // and is complete when step() returns true. False — and nothing begun — when the object is not prepared, a
    // pointer is null, the stretch is longer than prepare() was told, or it is too short to give every offset at
    // least two frames.
    [[nodiscard]] bool begin (Transform t, const float* left, const float* right, int n, float* curve) noexcept
    {
        active_ = false;
        if (! prepared_ || left == nullptr || right == nullptr || curve == nullptr) return false;
        if (n > maxSamples_ || n < minSamplesOf (t, 2)) return false;
        transform_ = t;
        n_ = n;
        curve_ = curve;
        next_ = 0;
        phaseReady_ = -1;
        x_[0] = left;
        x_[1] = right;
        if (t == Transform::Celt)
        {
            for (int c = 0; c < 2; ++c)
            {
                const float* x = c == 0 ? left : right;
                float* p = pre_[c].data();
                p[0] = x[0];
                for (int i = 1; i < n; ++i) p[i] = x[i] - (float) kCeltPreEmphasis * x[i - 1];
                x_[c] = p;
            }
        }
        acc_.attach (floats_.data(), doubles_.data(), linesOf (t));
        frames_ = t == Transform::Mp3 ? 0 : (n - frameOf (t)) / hopOf (t);
        active_ = true;
        return true;
    }

    int offsets() const noexcept { return active_ ? hopOf (transform_) : 0; }
    int offsetsDone() const noexcept { return next_; }

    // At most `maxOffsets` offsets more. True when every offset is done (and when nothing was begun).
    bool step (int maxOffsets) noexcept
    {
        if (! active_) return true;
        const int total = hopOf (transform_);
        for (int i = 0; i < maxOffsets && next_ < total; ++i, ++next_)
        {
            // Mp3 walks its offsets polyphase phase by phase, so that the subband samples of a phase are computed
            // once for its 18 granule starts; the others walk them in order. Offset 0 is first in both.
            const int offset = transform_ == Transform::Mp3 ? (next_ % 18) * 32 + next_ / 18 : next_;
            if (next_ == 0)
            {
                acc_.beginReference();
                forEachFrame (offset, [this] (const float* l, const float* r) noexcept { acc_.addReferenceFrame (l, r); });
                if (! acc_.finishReference()) { active_ = false; return true; }
            }
            acc_.beginOffset();
            forEachFrame (offset, [this] (const float* l, const float* r) noexcept { acc_.addFrame (l, r); });
            acc_.finishOffset (curve_ + (std::size_t) offset * kCells);
        }
        return next_ >= total;
    }

    // The zero share at one offset of the stretch last scanned to completion: per band (32), the largest of the
    // four signals' shares of coefficients 40 dB under their own level. `map` (optional) takes up to `maxFrames`
    // rows of 32 bytes — the same share per frame, 0..255 — and the number of rows written is returned.
    int zeroProfile (int offset, float* shareByBand, std::uint8_t* map = nullptr, int maxFrames = 0) noexcept
    {
        for (int b = 0; b < kZeroBands; ++b) shareByBand[b] = 0.0f;
        if (! active_ || next_ < hopOf (transform_) || offset < 0 || offset >= hopOf (transform_)) return 0;
        // the reference is the scan's own: it is still in the accumulator
        std::uint32_t counts[kSignals * kZeroBands] {};
        int frames = 0, rows = 0;
        if (transform_ == Transform::Mp3) phaseReady_ = -1;
        forEachFrame (offset, [&] (const float* l, const float* r) noexcept
        {
            float share[kZeroBands];
            acc_.addZeroFrame (l, r, counts, share);
            if (map != nullptr && rows < maxFrames)
            {
                for (int b = 0; b < kZeroBands; ++b) map[(std::size_t) rows * kZeroBands + (std::size_t) b] = (std::uint8_t) (share[b] * 255.0f + 0.5f);
                ++rows;
            }
            ++frames;
        });
        if (frames <= 0) return 0;
        const int perBand = linesOf (transform_) / kZeroBands;
        for (int b = 0; b < kZeroBands; ++b)
        {
            std::uint32_t most = 0;
            for (int s = 0; s < kSignals; ++s) most = std::max (most, counts[s * kZeroBands + b]);
            shareByBand[b] = (float) most / (float) (perBand * frames);
        }
        return rows;
    }

private:
    static constexpr int frameOf (Transform t) noexcept { return t == Transform::Celt ? kCeltFrame : kAacFrame; }

    template <class F>
    void forEachFrame (int offset, F&& frame) noexcept
    {
        if (transform_ == Transform::Mp3)
        {
            const int phase = offset % 32, start = offset / 32;
            const int count = Mp3Hybrid::subbandSamples (n_, phase);
            if (phaseReady_ != phase)
            {
                mp3_.subbands (x_[0], n_, phase, sub_[0].data());
                mp3_.subbands (x_[1], n_, phase, sub_[1].data());
                phaseReady_ = phase;
            }
            const int granules = Mp3Hybrid::granules (count, start);
            for (int q = 0; q < granules; ++q)
            {
                const std::size_t first = (std::size_t) (start + Mp3Hybrid::kBlock * q) * Mp3Hybrid::kSubbands;
                mp3_.granule (sub_[0].data() + first, lines_[0].data());
                mp3_.granule (sub_[1].data() + first, lines_[1].data());
                frame (lines_[0].data(), lines_[1].data());
            }
            return;
        }
        Mdct<Fft>& mdct = transform_ == Transform::AacSine ? sine_ : transform_ == Transform::AacKbd ? kbd_ : celt_;
        const int hop = mdct.hop();
        for (int f = 0; f < frames_; ++f)
        {
            mdct.transform (x_[0] + offset + f * hop, lines_[0].data());
            mdct.transform (x_[1] + offset + f * hop, lines_[1].data());
            frame (lines_[0].data(), lines_[1].data());
        }
    }

    bool prepared_ = false, active_ = false;
    int maxSamples_ = 0, n_ = 0, next_ = 0, frames_ = 0, phaseReady_ = -1;
    Transform transform_ = Transform::Mp3;
    const float* x_[2] { nullptr, nullptr };
    float* curve_ = nullptr;
    Mp3Hybrid mp3_;
    Mdct<Fft> sine_, kbd_, celt_;
    LevelAccumulator acc_;
    std::vector<float> sub_[2], pre_[2], lines_[2], floats_;
    std::vector<double> doubles_;
};

} // namespace felitronics::codecgrid
