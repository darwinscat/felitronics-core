// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace felitronics::analysis
{

// How the channels of one frame become the one value the envelope is built from. The site's `mix` option,
// code for code: 'avr' · 'L' · 'R' · 'max'.
enum class PeakMix : int
{
    Average = 0,    // (sum of ALL channels) / channel count — signed; the absolute value is taken per box
    Left    = 1,    // channel 0
    Right   = 2,    // the LAST channel, not channel 1: on a 6-channel file 'R' reads channel 5
    Max     = 3     // the loudest |channel| of the frame; the box mean of those, no second abs
};

//==============================================================================
// felitronics::analysis::WaveformPeaks — the "bars" a waveform player draws: a whole file reduced to `buckets`
// absolute max-abs values, after a box-average decimation to about 8 kHz.
//
// THESE ARE NOT PEAKS IN ANY METERING SENSE, and the name says "waveform" for that reason. A bucket is the largest
// |mean| of ~125 µs boxes, so it sits at or below the sample peak of its stretch — up to the rounding of the box
// mean, which can lift a mean of equal samples a few ulps above them (2^30 frames of 1+3·2^-23 at an absurd rate
// read 4.3e-14 over) — and further below the true peak.
// Neither of the core's two true-peak meters is involved, and those two must not be confused with each other either:
// `fcore::Probe` (the 128-tap reference, what `fc_probe_tp_linear` reports) and `analysis::TruePeakMeter` (the spec's
// 48-tap filter, what the mastering chain runs) are DIFFERENT filters and read different numbers. How far apart is
// measured by P62, not stated here. A picture that must show a true peak next to these bars names which one it shows.
//
// THIS IS A PORT, AND THE SPEC IS EXECUTABLE: the site's `computePeaksFromBuffer` and `peaksFromWav`
// (audio-peaks.js). Before this header the same picture had four definitions — a Java sidecar generator
// decoding to 8 kHz s16, a python one, the page's JavaScript and a costume on top of the sidecars — and a demo
// file and an uploaded one could be drawn by different ones in the same interface. This is the one definition both
// roads are meant to take — the native CLI a server runs and the browser module; moving the site's generators and
// page onto it is the site's work and had not happened when this was written — and it is the JavaScript's to the bit: a cross-language NULL runs the same float32 PCM
// through both. Anything below that reads like an improvement is a change to the picture, not to this file.
//
// THE DEFINITION, step by step (the spec's variable names in brackets):
//   decim      = max(1, round(sampleRate / 8000))             JavaScript's Math.round — ties toward +infinity
//   decLen     = max(buckets, floor(frames / decim))
//   perBucket  = decLen / buckets                              [decPerBucket], a double
//   each frame → one value v by `PeakMix`; every `decim` values → one box, d = |mean| (Max: the mean itself)
//   the running max of d is emitted when the box counter reaches the bucket's boundary (bucket + 1) · perBucket
// Everything is binary64 in the spec's order — a sample is promoted to double before it is touched — so the
// arithmetic here IS the arithmetic there, and there is no product in it for FP contraction to fuse.
//
// TWO OUTPUT TYPES, BECAUSE THE SPEC HAS TWO. `computePeaksFromBuffer` stores the maxima in a plain Array, i.e.
// as doubles; `peaksFromWav` stores them in a Float32Array, i.e. rounded to float32. Same arithmetic, different
// last step, so the conversion is explicit and named: `peaks()` is the first, `peakAsFloat32()` the second. (What
// this header cannot reproduce is `peaksFromWav` on a 32-bit integer or float64 WAV: it reads those samples as
// doubles, and a float32 plane has already rounded them. That is a compatibility decision about the input, not a
// property of the reduction, and it stays outside this class.)
//
// THE SPEC'S KNOWN PROPERTIES ARE KEPT, not fixed — each is a change to the picture and each has a pinned witness:
//   * THE LAST BUCKET CAN BE LOST. The boundary is (bucket + 1) · perBucket, a product of an already rounded
//     quotient; at decLen = 2007, buckets = 1000 the last one is 2007.0000000000002, the counter tops out at
//     2007, and bucket 999 is never emitted — it reads 0.
//   * A tail shorter than one box is dropped: its frames never form a box.
//   * When floor(frames / decim) < buckets the boundaries outrun the boxes and the trailing buckets stay 0.
//   * The effective envelope rate depends on the file's rate: 44.1 kHz → decim 6 → 7350 Hz, 48 kHz → 8000 Hz.
//   * A NaN value never wins `d > max`, so a box containing one is skipped (its mean is NaN); +Inf wins. In 'max'
//     mode a NaN CHANNEL never wins `|x| > fv` either, so it reads as that frame's other channels.
//   * 'avr' divides per FRAME, channel 0 first: three channels at 16 kHz on [3f35e00d bf089722 3eb88b0a] and
//     [3ef22e5b bf11b04d bf347cf7] read 0x3fa6828cd5555556, where one division of all six samples at the end reads ...555.
//
// NOT COVERED BY PARITY, because it is runtime state (law 10): FTZ/DAZ. A subnormal sample peaks at 1.4e-45 in
// the spec and at 0 on a host thread that flushes denormals. Run with them off when the bits must match.
//
// STREAMING. `prepare()` takes the file's TOTAL length, because every boundary depends on it; a stream of
// unknown length is not supported and cannot be. `process()` may be called with any split of the file and gives
// the bits one call would — nothing here depends on where a call ends. The width is EXACT (law 11c): 'avr'
// divides by the channel count and 'R' reads the last channel, so a narrower call would be a different picture,
// and is refused. A call that would run past the prepared length is refused as a whole, before anything moves.
//
// RT-safe: prepare() is the only allocation; process() does no alloc / lock / throw.
class WaveformPeaks
{
public:
    static constexpr double kEnvelopeHz = 8000.0;    // the spec's ENVELOPE_HZ
    static constexpr int    kDefaultBuckets = 1000;   // the spec's PEAK_BUCKETS

    // CAPACITY GUARDS — not part of the picture. Buckets: 8 bytes each, and a browser worker cannot survive the
    // allocation of an absurd count (law 11d), so it is refused instead. Frames: every frame and box counter is
    // compared against a double boundary, which is exact only below 2^53.
    static constexpr int           kMaxBuckets = 1 << 20;
    static constexpr std::uint64_t kMaxFrames  = std::uint64_t (1) << 53;

    // JavaScript's Math.round for a finite value: the nearest integer, a tie toward +infinity. Not std::round
    // (ties away from zero) and not floor(x + 0.5), which rounds 0.49999999999999994 up.
    static double jsRound (double x) noexcept
    {
        const double f = std::floor (x);
        return (x - f >= 0.5) ? f + 1.0 : f;
    }

    // The decimation factor for a rate, exactly as the spec derives it. Refuses what is not a rate.
    [[nodiscard]] static bool decimationFor (double sampleRate, int& out) noexcept
    {
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return false;
        const double d = std::max (1.0, jsRound (sampleRate / kEnvelopeHz));
        if (! (d <= 2147483647.0)) return false;
        out = (int) d;
        return true;
    }

    [[nodiscard]] bool prepare (double sampleRate, int numChannels, std::uint64_t totalFrames,
                                int buckets = kDefaultBuckets, PeakMix mix = PeakMix::Average)
    {
        prepared_ = false;
        int decim = 0;
        if (! decimationFor (sampleRate, decim)) return false;
        if (numChannels < 1) return false;
        if (totalFrames < 1 || totalFrames > kMaxFrames) return false;
        if (buckets < 1 || buckets > kMaxBuckets) return false;
        if (mix != PeakMix::Average && mix != PeakMix::Left && mix != PeakMix::Right && mix != PeakMix::Max)
            return false;

        nch_    = numChannels;
        frames_ = totalFrames;
        decim_  = decim;
        mix_    = mix;
        const double decLen = std::max ((double) buckets, std::floor ((double) totalFrames / (double) decim));
        perBucket_ = decLen / (double) buckets;
        out_.assign ((std::size_t) buckets, 0.0);
        prepared_ = true;
        reset();
        return true;
    }

    // Back to the head of the same file: every bucket 0, nothing seen.
    void reset() noexcept
    {
        std::fill (out_.begin(), out_.end(), 0.0);
        bucket_ = 0; boundary_ = perBucket_; max_ = 0.0; acc_ = 0.0; accN_ = 0; boxes_ = 0; seen_ = 0;
    }

    [[nodiscard]] bool process (const float* const* channels, int numChannels, int n) noexcept
    {
        if (! prepared_ || channels == nullptr || numChannels != nch_ || n < 0) return false;
        if ((std::uint64_t) n > frames_ - seen_) return false;
        for (int c = 0; c < nch_; ++c) if (n > 0 && channels[c] == nullptr) return false;

        const int buckets = (int) out_.size();
        const bool absMode = (mix_ == PeakMix::Max);
        for (int j = 0; j < n; ++j)
        {
            double fv;
            if (mix_ == PeakMix::Left)       fv = (double) channels[0][j];
            else if (mix_ == PeakMix::Right) fv = (double) channels[nch_ - 1][j];
            else if (absMode)
            {
                fv = 0.0;
                for (int c = 0; c < nch_; ++c) { const double x = std::fabs ((double) channels[c][j]); if (x > fv) fv = x; }
            }
            else
            {
                double s = 0.0;
                for (int c = 0; c < nch_; ++c) s += (double) channels[c][j];
                fv = s / (double) nch_;
            }
            acc_ += fv;
            ++accN_;
            if (accN_ >= decim_)
            {
                const double mean = acc_ / (double) accN_;
                const double d = absMode ? mean : std::fabs (mean);
                if (d > max_) max_ = d;
                acc_ = 0.0; accN_ = 0; ++boxes_;
                if ((double) boxes_ >= boundary_ && bucket_ < buckets)
                {
                    out_[(std::size_t) bucket_++] = max_;
                    max_ = 0.0;
                    boundary_ = (double) (bucket_ + 1) * perBucket_;
                }
            }
        }
        seen_ += (std::uint64_t) n;
        return true;
    }

    [[nodiscard]] bool prepared() const noexcept { return prepared_; }
    // The whole prepared length has been seen: the buckets are final. Before that they are a prefix in progress.
    [[nodiscard]] bool complete() const noexcept { return prepared_ && seen_ == frames_; }

    std::span<const double> peaks() const noexcept { return { out_.data(), out_.size() }; }
    // The float32 form — `peaksFromWav`'s Float32Array store. Round to nearest, ties to even, as ToFloat32 is.
    float peakAsFloat32 (int i) const noexcept { return (float) out_[(std::size_t) i]; }

    int           buckets()      const noexcept { return (int) out_.size(); }
    int           decimation()   const noexcept { return decim_; }
    int           numChannels()  const noexcept { return nch_; }
    PeakMix       mix()          const noexcept { return mix_; }
    std::uint64_t totalFrames()  const noexcept { return frames_; }
    std::uint64_t framesSeen()   const noexcept { return seen_; }
    // How many buckets were actually emitted. Less than buckets() on a complete file is the spec's own
    // behaviour (see "THE LAST BUCKET CAN BE LOST"), and the only way a caller can tell a 0 from a silence.
    int           bucketsEmitted() const noexcept { return bucket_; }

private:
    std::vector<double> out_;
    double        perBucket_ = 0.0, boundary_ = 0.0, max_ = 0.0, acc_ = 0.0;
    std::uint64_t frames_ = 0, seen_ = 0, boxes_ = 0;
    int           nch_ = 0, decim_ = 1, accN_ = 0, bucket_ = 0;
    PeakMix       mix_ = PeakMix::Average;
    bool          prepared_ = false;
};

} // namespace felitronics::analysis
