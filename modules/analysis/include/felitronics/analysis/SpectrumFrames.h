// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/OfflineFft.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::SpectrumFrames — the shared FRAME PRODUCER the offline analyzers measure on.
// It is deliberately NOT a "long-term spectrum": it produces one calibrated power spectrum per WINDOW,
// on an absolute schedule, and the consumer decides what to do with each frame — average them all
// (SourceForensics), keep only the ones inside a quiet stretch (HumDetector), or fold a band of bins
// into semitones (LowEnd). An averaging primitive would fit the first and fit the other two not at all.
//
// WHY THIS IS SHARED AND THE AVERAGING IS NOT. Every known way a windowed analyzer silently stops being
// re-blocking-invariant lives in this layer, not above it: one FFT per process() call over "whatever
// arrived"; a frame trigger anchored to the frame START (`pos % hop == 0` fires at pos == hop, on a
// frame that does not exist yet, and never agrees with the end-anchored rule when N % hop != 0); a
// "samples since the last frame" counter that some path resets; zero-padding the tail inside the last
// process(), where "last" is unknowable; reducing the two physical halves of a wrapped ring separately;
// windowing a reused buffer twice. One implementation, tested once, spends that risk once.
//
// THE CONTRACT (law 8a — bit-identical under ARBITRARY re-slicing, which is stronger than law 11(a)'s
// same-boundaries promise):
//   * One integer clock. `totalSamples_` counts samples CONSUMED since reset(). push() writes one
//     channel's sample at the current time; tick() closes the sample and advances the clock once.
//   * Frame f covers exactly [f*hop, f*hop + N) and is emitted the instant sample f*hop + N - 1 has
//     been consumed — i.e. when `totalSamples_ >= N && (totalSamples_ - N) % hop == 0`. The trigger
//     reads the clock and nothing else, so no call boundary can move it.
//   * maxBlock sizes NOTHING here. A report must be identical for prepare(sr, 64, nch) and
//     prepare(sr, 8192, nch).
//   * THE TAIL IS NOT TRANSFORMED. finish() emits no extra frame; it publishes
//     tailUncoveredSamples() instead. A zero-padded frame is a DIFFERENT measurement (truncated window,
//     power scaled by (N-m)/N, plus leakage off the discontinuity), and a right-aligned tail frame
//     overlaps its predecessor and over-weights the end of the programme in any average. A programme
//     shorter than N produces no frame at all — that is `valid = false, ShorterThanWindow` upstairs,
//     never a zero that reads as "nothing found".
//   * A frame containing ANY hole (a non-finite sample, or a channel absent from that call) is emitted
//     with frameFinite(c) == false and is NOT transformed: holes are counted, never guessed.
//
// PRECISION. The transform is core::offline::fftInplace — double, in place, noexcept, and it does NOT
// allocate (the "ALLOCATES" banner on that header is about its convolve/magSpectrum wrappers). The
// scratch is owned by prepare(). Double, not the float core::fft seam, for three reasons: the float
// seam's backend concept has no allocation query, so an honest Storage/law-11d budget cannot include
// it; it carries a template parameter five sessions would have to agree on; and its pffft ON/OFF CI
// rows transform to different bits. A 2^17 frame costs ~11 Mflop, i.e. a second or two for a
// five-minute programme — this runs once per file.
//
// CALIBRATION. power(c)[k] = |X_k|^2 / (N * sum(w^2)), one-sided is NOT folded (bin k is bin k), so a
// band share is dimensionless and comparable with a time-domain energy ratio. Both window gains are
// published FROM THE WINDOW AS BUILT rather than from a textbook constant: windowPowerGain() = sum(w^2)/N
// (Parseval) and windowCoherentGain() = sum(w)/N (= 0.5 for Hann; a sinusoid's amplitude).
//
// RT: this is an OFFLINE instrument (message thread). prepare() allocates; push/tick/finish do not.
struct SpectrumFramesParams
{
    int fftOrder = 17;      // N = 1 << fftOrder. 17 at 48 kHz = 131072 samples, 2.73 s, 0.366 Hz bins —
                            // what it takes to RESOLVE 49.0 from 50.0 Hz (two Hann peaks need ~2 bins).
    int hop = 0;            // 0 means N/2. Must satisfy 0 < hop <= N.
};

class SpectrumFrames
{
public:
    // WHAT prepare() ASKS THE HEAP FOR (law 11d), from the function prepare() sizes itself with, so a
    // caller budgeting memory reads the numbers the object is actually built from.
    struct Storage
    {
        bool ok = false;
        std::size_t ringSamples = 0;     // float, per channel * channels
        std::size_t holeFlags   = 0;     // uint8, per channel * channels
        std::size_t windowD     = 0;     // double, N
        std::size_t fftComplex  = 0;     // complex<double>, N — ONE scratch, reused per channel
        std::size_t powerD      = 0;     // double, (N/2 + 1) * channels
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) sizeof (float) * ringSamples
                 + (std::uint64_t) holeFlags
                 + (std::uint64_t) sizeof (double) * (windowD + powerD)
                 + (std::uint64_t) sizeof (std::complex<double>) * fftComplex;
        }
    };

    static Storage storageFor (double sampleRate, int maxChannels, const SpectrumFramesParams& p) noexcept
    {
        Storage s;
        if (! (sampleRate > 0.0) || ! std::isfinite (sampleRate)) return s;
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return s;
        if (p.fftOrder < 4 || p.fftOrder > 22) return s;               // 16 .. 4194304 points
        const std::size_t n = (std::size_t) 1u << p.fftOrder;
        const std::size_t hop = p.hop == 0 ? n / 2u : (std::size_t) (p.hop < 0 ? 0 : p.hop);
        if (hop == 0 || hop > n) return s;
        const std::size_t ch = (std::size_t) maxChannels;
        s.ok = true;
        s.ringSamples = n * ch;
        s.holeFlags   = n * ch;
        s.windowD     = n;
        s.fftComplex  = n;
        s.powerD      = (n / 2u + 1u) * ch;
        return s;
    }

    void setParams (const SpectrumFramesParams& p) noexcept { params_ = p; }   // takes effect at the next prepare()

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        prepared_ = false;                                             // law 11b: disarm, validate, write
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! st.ok) return false;
        sampleRate_ = sampleRate;
        channels_   = maxChannels;
        n_   = (std::int64_t) 1 << params_.fftOrder;
        hop_ = params_.hop == 0 ? n_ / 2 : (std::int64_t) params_.hop;
        bins_ = (int) (n_ / 2 + 1);
        ring_.assign (st.ringSamples, 0.0f);
        hole_.assign (st.holeFlags, (std::uint8_t) 1);                 // an unwritten slot is a hole
        window_.assign (st.windowD, 0.0);
        scratch_.assign (st.fftComplex, std::complex<double> {});
        power_.assign (st.powerD, 0.0);
        buildWindow();
        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        totalSamples_ = 0;
        frameIndex_   = -1;
        frameStart_   = 0;
        finished_     = false;
        lastFrameEnd_ = 0;
        for (auto& v : ring_) v = 0.0f;
        for (auto& v : hole_) v = (std::uint8_t) 1;
        for (auto& v : power_) v = 0.0;
        for (int c = 0; c < core::kMaxChannels; ++c) { holesInWindow_[c] = 0; frameFinite_[c] = false; }
        if (prepared_) for (int c = 0; c < channels_; ++c) holesInWindow_[c] = n_;   // the ring starts all holes
    }

    // One channel's sample at the CURRENT time. `fed == false` marks a hole (the channel was absent from
    // this call, or the sample is non-finite): it is written as zero, counted, and poisons any frame that
    // contains it. Out-of-range channels are ignored rather than trapped: this is a sink.
    void push (int c, float x, bool fed) noexcept
    {
        if (! prepared_ || finished_ || c < 0 || c >= channels_) return;
        const bool hole = ! fed || ! std::isfinite (x);
        const std::size_t idx = (std::size_t) (totalSamples_ % n_) + (std::size_t) c * (std::size_t) n_;
        holesInWindow_[c] -= (std::int64_t) hole_[idx];                 // the slot being overwritten leaves the window
        hole_[idx] = (std::uint8_t) (hole ? 1 : 0);
        holesInWindow_[c] += (std::int64_t) (hole ? 1 : 0);
        ring_[idx] = hole ? 0.0f : x;
    }

    // Close the sample: advance the clock ONCE, after every channel of this sample has been pushed.
    // Returns true when a frame completed on this very sample — its power() is then readable until the
    // next tick() that completes another one.
    [[nodiscard]] bool tick() noexcept
    {
        if (! prepared_ || finished_) return false;
        ++totalSamples_;
        if (totalSamples_ < n_) return false;
        if ((totalSamples_ - n_) % hop_ != 0) return false;
        frameStart_   = totalSamples_ - n_;
        frameIndex_   = frameStart_ / hop_;
        lastFrameEnd_ = totalSamples_;
        for (int c = 0; c < channels_; ++c) transform (c);
        return true;
    }

    // End of stream. Emits NOTHING (see the tail contract above) and freezes: push/tick refuse until
    // reset(). Idempotent.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        finished_ = true;
    }

    // --- the frame just closed ---
    std::int64_t frameIndex() const noexcept { return frameIndex_; }                 // -1 before the first
    std::int64_t frameStart() const noexcept { return frameStart_; }
    std::int64_t frameEnd()   const noexcept { return frameStart_ + n_; }
    bool frameFinite (int c) const noexcept { return c >= 0 && c < channels_ && frameFinite_[c]; }
    // (n/2 + 1) doubles. All zero for a frame that held a hole — read frameFinite() before believing it.
    const double* power (int c) const noexcept
    {
        if (c < 0 || c >= channels_ || power_.empty()) return nullptr;
        return power_.data() + (std::size_t) c * (std::size_t) bins_;
    }

    // --- geometry, published so a consumer states its own resolution honestly ---
    int bins()   const noexcept { return bins_; }
    std::int64_t windowSamples() const noexcept { return n_; }
    std::int64_t hopSamples()    const noexcept { return hop_; }
    double binHz() const noexcept { return n_ > 0 ? sampleRate_ / (double) n_ : 0.0; }
    double sampleRate() const noexcept { return sampleRate_; }
    double windowPowerGain()    const noexcept { return sumW2_ / (double) n_; }       // Parseval
    double windowCoherentGain() const noexcept { return sumW_  / (double) n_; }       // 0.5 for Hann
    std::int64_t totalSamples() const noexcept { return totalSamples_; }
    // How much of the stream's end no frame covered. NOT an error — a number the consumer weighs.
    std::int64_t tailUncoveredSamples() const noexcept { return totalSamples_ - lastFrameEnd_; }
    std::int64_t frameCount() const noexcept { return frameIndex_ + 1; }
    bool isFinished() const noexcept { return finished_; }
    bool isPrepared() const noexcept { return prepared_; }
    static constexpr int latencySamples() noexcept { return 0; }                      // a read-only sink

private:
    void buildWindow() noexcept
    {
        // Hann, periodic (N, not N-1), built ONCE in double.
        //
        // TWO THINGS HERE ARE DELIBERATE AND BOTH WERE MEASURED.
        //
        // 1. core::det::cos, NOT std::cos. The system libm is not the same function on every row: the
        //    coefficients of this very window differ in 502 of 16384 places (order 14) and 4032 of
        //    131072 (order 17) between Apple's libm and musl's, and every power bin is multiplied by
        //    them. det::cos is one implementation compiled into every build, so the window is the same
        //    on the developer's Mac, on the CI row and in the browser. It is within 1 ulp of the
        //    correctly rounded value on the worst 64 of four million adversarial arguments — the same
        //    bound the system libm holds — so nothing is given up for it.
        //
        // 2. ONE QUADRANT, THREE REFLECTIONS. Computing every index from its own argument does NOT give
        //    a symmetric window: 2*pi*i/N and 2*pi*(N-i)/N are different doubles (5332 of 16383 differ
        //    from the exact reflection), so w[i] != w[N-i] in 10314 of 16383 places — on det::cos AND on
        //    the system libm alike, which is why no comparison against a libm or against a
        //    high-precision oracle could ever have shown it: both sides share the defect. It took an
        //    identity that needs no oracle at all to see it. Deriving the other three quadrants by
        //    index makes the symmetry EXACT by construction, and costs a quarter of the calls
        //    (32769 instead of 131072 at order 17). It also keeps every argument inside [0, pi/2],
        //    away from the neighbourhood of a multiple of pi/2 where argument reduction is hardest.
        const std::int64_t half = n_ / 2, quarter = n_ / 4;
        for (std::int64_t i = 0; i <= quarter; ++i)
        {
            const double c = core::det::cos (2.0 * core::kPi * (double) i / (double) n_);
            const double lo = 0.5 - 0.5 * c;           // w[i] and w[N-i]
            const double hi = 0.5 + 0.5 * c;           // w[N/2-i] and w[N/2+i]
            window_[(std::size_t) i] = lo;
            if (i > 0)             window_[(std::size_t) (n_ - i)]    = lo;
            window_[(std::size_t) (half - i)] = hi;
            if (half + i < n_)     window_[(std::size_t) (half + i)]  = hi;
        }
        // The sums are accumulated in INDEX order, separately from the fill, so their rounding order is
        // a property of the window and not of the order the quadrants happened to be written in.
        sumW_ = 0.0; sumW2_ = 0.0;
        for (std::int64_t i = 0; i < n_; ++i)
        {
            const double w = window_[(std::size_t) i];
            sumW_  += w;
            sumW2_ += w * w;
        }
    }

    void transform (int c) noexcept
    {
        const std::size_t base = (std::size_t) c * (std::size_t) n_;
        if (holesInWindow_[c] > 0)                                       // a holed frame is not guessed
        {
            frameFinite_[c] = false;
            double* p = power_.data() + (std::size_t) c * (std::size_t) bins_;
            for (int k = 0; k < bins_; ++k) p[k] = 0.0;
            return;
        }
        frameFinite_[c] = true;
        // Unwrap CHRONOLOGICALLY into the scratch — the oldest sample of the window sits at
        // ring[totalSamples_ % n], since exactly n samples have been written since it.
        const std::size_t start = (std::size_t) (totalSamples_ % n_);
        for (std::int64_t i = 0; i < n_; ++i)
        {
            const std::size_t r = base + (start + (std::size_t) i) % (std::size_t) n_;
            scratch_[(std::size_t) i] = std::complex<double> ((double) ring_[r] * window_[(std::size_t) i], 0.0);
        }
        core::offline::fftInplace (scratch_, -1);
        const double norm = 1.0 / ((double) n_ * sumW2_);
        double* p = power_.data() + (std::size_t) c * (std::size_t) bins_;
        for (int k = 0; k < bins_; ++k)
        {
            const std::complex<double>& z = scratch_[(std::size_t) k];
            p[k] = (z.real() * z.real() + z.imag() * z.imag()) * norm;
        }
    }

    SpectrumFramesParams params_ {};
    bool prepared_ = false, finished_ = false;
    double sampleRate_ = 48000.0;
    int channels_ = 0, bins_ = 0;
    std::int64_t n_ = 0, hop_ = 0;
    std::int64_t totalSamples_ = 0, frameIndex_ = -1, frameStart_ = 0, lastFrameEnd_ = 0;
    double sumW_ = 0.0, sumW2_ = 0.0;
    std::int64_t holesInWindow_[core::kMaxChannels] {};
    bool frameFinite_[core::kMaxChannels] {};
    std::vector<float> ring_;
    std::vector<std::uint8_t> hole_;
    std::vector<double> window_;
    std::vector<std::complex<double>> scratch_;
    std::vector<double> power_;
};

} // namespace felitronics::analysis
