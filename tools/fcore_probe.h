// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// fcore::Probe — the measurement body shared VERBATIM by the native reference CLI (tools/fcore_measure.cpp)
// and the wasm shim (tools/wasm/fc_probe.cpp), so that the P0 spike's claim — "the same code, the same
// doubles" — is a fact about one translation unit rather than a hope about two.
//
// Sharing the source removes ONE variable. It does not by itself make the two sides agree: different
// compilers (Apple clang vs the emscripten upstream), different libm, and different FP-contraction defaults
// still sit between this source and the doubles it produces. What closes the gap is the build contract:
//
//   * the native side MUST be compiled -ffp-contract=off. Baseline wasm has no scalar FMA instruction, so
//     emscripten cannot contract; a contracting native build differs. Measured on real audio: the true-peak
//     answer moves 6.75e-7 dB (exactly 1 float ulp in the polyphase accumulator) between contract=on and off.
//   * the wasm side MUST NOT be built -mrelaxed-simd. f64x2.relaxed_madd is implementation-defined (fused on
//     hosts with FMA, unfused elsewhere), which breaks determinism between MACHINES, not merely between
//     tiers. Plain -msimd128 lowers to f64x2.mul/add and is safe.
//
// WHAT TO COMPARE. Not the integrated LUFS: it is DISCONTINUOUS in its inputs, because the BS.1770 gates are
// strict comparisons — a 400 ms block sitting within ~1e-12 of a gate flips its inclusion between two builds
// and moves the reading by ~0.01 dB, seven orders above any sane tolerance. Compare gatingBlockEnergies():
// the pre-gate block energies are continuous in the input samples and can carry a bit-exactness claim. The
// true-peak LINEAR maximum is likewise continuous (a 1-ulp input moves a max by at most 1 ulp) and is
// compared directly; its dB form is derived once at the end and is for humans.
//
// TRUE-PEAK CONFIG IS PART OF THE CONTRACT. The core holds two different true-peak filters, and they disagree
// by 0.0009–0.0028 dB on real audio: this path (PolyphaseOversampler at 4x, 32 taps/phase = a 128-tap
// prototype, cutoff 0.90x Nyquist, Kaiser beta 9) and analysis::TruePeakMeter (48 taps, full base Nyquist,
// beta 8, factor chosen by rate). A shim that reached for the "obvious" TruePeakMeter would fail the spike's
// acceptance by ~3e-3 dB and the failure would read as a wasm bug. Both sides go through THIS class.

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/core/Config.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace fcore
{

class Probe
{
public:
    // The streaming step. Callers may hand process() any length: it walks the input in kChunk-frame steps
    // internally, so a whole-file buffer costs the same bounded scratch as a stream (the true-peak scratch is
    // kChunk*4 floats — 128 KB — not 4*frames, which would be 247 MB for a 5-minute stereo track on wasm32).
    // Chunking cannot change the arithmetic: LoudnessMeter::process() does identical per-sample work for any
    // n, and the oversampler's ring history makes upsample() a pure function of the samples seen so far —
    // verified bit-identical for chunk sizes 1 … 100003.
    static constexpr int kChunk         = 8192;
    static constexpr int kOsFactor      = 4;    // \ the true-peak filter this tool is the REFERENCE for;
    static constexpr int kOsTapsPerPhase = 32;  // / see the header comment before changing either

    // Audio sample rates, bounded to a range the downstream arithmetic survives. "Positive and finite" is NOT
    // enough: LoudnessMeter sizes its hop with `lround(0.01*fs)` into an int and its block store with
    // `ceil(maxDurationSec*fs)`, so an absurd-but-finite rate (1e300, or Number.MIN_VALUE from a page)
    // reaches an out-of-range lround — undefined behaviour — or asks for an impossible allocation. The bounds
    // are generous: every rate any product here will meet lies inside them by orders of magnitude.
    static constexpr double kMinSampleRate = 1000.0;
    static constexpr double kMaxSampleRate = 768000.0;

    // maxDurationSec sizes the meter's gating-block store. The default holds 4 h at any rate (1.27 MB) —
    // droppedBlocks() reads 0 for anything shorter, and a caller checks it rather than assuming.
    bool prepare (double sampleRate, int channels, double maxDurationSec = 4.0 * 3600.0)
    {
        prepared_ = false;
        finished_ = false;
        // The compound test rejects 0, negatives, NaN and +inf, and then the absurd-but-finite rates too.
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return false;
        if (channels < 1 || channels > felitronics::core::kMaxChannels) return false;
        if (! (maxDurationSec > 0.0) || ! std::isfinite (maxDurationSec)) return false;

        nc_ = channels;
        lm_.prepare (sampleRate, nc_, maxDurationSec);
        os_.assign ((std::size_t) nc_, {});
        for (auto& o : os_) if (! o.prepare (kOsFactor, 1, kOsTapsPerPhase)) return false;
        osBuf_.assign ((std::size_t) kChunk * (std::size_t) kOsFactor, 0.0f);
        maxTp_ = 0.0;
        samplePeak_ = 0.0;
        prepared_ = true;
        return true;
    }

    bool prepared() const noexcept { return prepared_; }

    // planar[c] holds n frames for channel c. Channels beyond the prepared count are ignored; fewer than
    // prepared is honoured as-is (the meter weights only what it is given).
    void process (const float* const* planar, int channels, long long n) noexcept
    {
        if (! prepared_ || finished_ || n <= 0) return;
        const int useCh = channels < nc_ ? channels : nc_;
        // A non-positive channel count is not "measure nothing": LoudnessMeter would still advance its hop
        // clock and record SILENT gating blocks, quietly lengthening the program with material that was never
        // submitted. Refuse instead.
        if (useCh <= 0) return;

        const float* view[felitronics::core::kMaxChannels] {};
        for (long long off = 0; off < n; off += kChunk)
        {
            const int m = (int) std::min<long long> (kChunk, n - off);
            for (int c = 0; c < useCh; ++c) view[c] = planar[c] + off;

            lm_.process (view, useCh, m);

            for (int c = 0; c < useCh; ++c)
            {
                for (int i = 0; i < m; ++i)
                    samplePeak_ = std::max (samplePeak_, (double) std::fabs (view[c][i]));

                const float* in[1] { view[c] };
                float*       out[1] { osBuf_.data() };
                os_[(std::size_t) c].upsample (in, 1, m, out);
                const int upN = m * kOsFactor;
                for (int k = 0; k < upN; ++k)
                    maxTp_ = std::max (maxTp_, (double) std::fabs (osBuf_[(std::size_t) k]));
            }
        }
    }

    // Ends the measurement: drains the polyphase FIR so a peak in the final samples is actually seen.
    //
    // WHY THIS IS NOT OPTIONAL. The oversampler is causal with a group delay of (N-1)/2 = 63.5 oversampled
    // samples, so the reconstruction of the last ~16 baseband samples never leaves the filter while input is
    // still being fed. Without draining, a file that ENDS on a transient is catastrophically under-read: a
    // unit impulse in the final sample measures 0.000071 instead of 0.881, and a buffer whose last ten
    // samples sit at 0.95 measures 0.0152 (−36 dBTP) when the true peak is 1.0625 — i.e. ABOVE full scale.
    // A true-peak tool that reports −36 dBTP for a clipping ending is worse than no tool.
    //
    // Idempotent, and it ends the stream: process() must not be called again afterwards, because the drain
    // has pushed silence through the filter history. Only the oversampler is drained — feeding the loudness
    // meter would append spurious silence to the program.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        finished_ = true;
        float zeros[kOsTapsPerPhase] {};                  // fixed: prepare() does all allocation
        for (int c = 0; c < nc_; ++c)
        {
            const float* in[1] { zeros };
            float*       out[1] { osBuf_.data() };
            os_[(std::size_t) c].upsample (in, 1, kOsTapsPerPhase, out);
            const int upN = kOsTapsPerPhase * kOsFactor;
            for (int k = 0; k < upN; ++k)
                maxTp_ = std::max (maxTp_, (double) std::fabs (osBuf_[(std::size_t) k]));
        }
    }

    // --- what a cross-toolchain check compares (continuous in the input samples) ---
    int           gatingBlockCount()    const noexcept { return lm_.gatingBlockCount(); }
    std::span<const double> gatingBlockEnergies() const noexcept { return lm_.gatingBlockEnergies(); }

    // The true peak is never below the SAMPLE peak: the reconstructed signal passes through the samples by
    // construction. Enforcing that as a floor costs nothing and makes a whole class of filter-side mistake
    // (a mis-sized drain, a wrong prototype, a bad phase) impossible to hide.
    double truePeakLinear() const noexcept { return std::max (maxTp_, samplePeak_); }
    double samplePeakLinear() const noexcept { return samplePeak_; }

    // --- what a human reads (derived; routed through log10, so not the bit-exactness surface) ---
    double integratedLufs() const noexcept { return lm_.integratedLufs(); }
    // The dB of truePeakLinear(), NOT of the oversampler maximum: the two differ whenever the sample peak
    // is the larger of the two (a lone impulse, anything ending abruptly), and a tool whose dB and linear
    // answers disagree is worse than one that is merely wrong.
    double truePeakDb()     const noexcept
    {
        const double tp = truePeakLinear();
        return 20.0 * std::log10 (tp > 1e-9 ? tp : 1e-9);
    }
    int    droppedBlocks()  const noexcept { return lm_.droppedBlocks(); }

private:
    felitronics::analysis::LoudnessMeter                        lm_;
    std::vector<felitronics::oversampling::PolyphaseOversampler> os_;
    std::vector<float>                                          osBuf_;
    double                                                      maxTp_ = 0.0;
    double                                                      samplePeak_ = 0.0;
    int                                                         nc_ = 0;
    bool                                                        prepared_ = false;
    bool                                                        finished_ = false;
};

} // namespace fcore
