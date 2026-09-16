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
// TRUE-PEAK CONFIG IS PART OF THE CONTRACT, AND IT LIVES IN ONE CLASS. The core holds two different true-peak
// filters: analysis::ReferenceTruePeakMeter (PolyphaseOversampler at 4x, 32 taps/phase = a 128-tap prototype,
// cutoff 0.90x Nyquist, Kaiser beta 9, at every rate) and analysis::TruePeakMeter (the spec's 12 taps/phase,
// full base Nyquist, beta 8, factor chosen by rate). They read different numbers — how far apart is measured
// and pinned by felitronics_truepeak_instrument_gap_tests, not written here. This probe is the REFERENCE: it
// measures through ReferenceTruePeakMeter, which is also what TargetLoudnessSolver aims a delivered ceiling
// with, so the file a solve delivers and the number this tool prints for it come from the same arithmetic. A
// shim that reached for the "obvious" TruePeakMeter would disagree with both, and the failure would read as a
// wasm bug. Both sides go through THIS class.

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>
#include <felitronics/analysis/StereoColumns.h>
#include <felitronics/analysis/WaveformPeaks.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace fcore
{

class Probe
{
public:
    // The streaming step. Callers may hand process() any length: it walks the input in kChunk-frame steps
    // internally, so a whole-file buffer costs the same bounded scratch as a stream (the true-peak meter walks
    // its own fixed scratch — never 4*frames, which would be 247 MB for a 5-minute stereo track on wasm32).
    // Chunking cannot change the arithmetic: LoudnessMeter::process() does identical per-sample work for any
    // n, and the oversampler's ring history makes upsample() a pure function of the samples seen so far —
    // verified bit-identical for chunk sizes 1 … 100003.
    static constexpr int kChunk         = 8192;
    // The reference filter's topology, named once in analysis::ReferenceTruePeakMeter and only re-exported here
    // (fc_probe_os_factor/_taps, fcore_measure's banner) — a literal here would be a second copy that could drift.
    static constexpr int kOsFactor       = felitronics::analysis::ReferenceTruePeakMeter::kFactor;
    static constexpr int kOsTapsPerPhase = felitronics::analysis::ReferenceTruePeakMeter::kTapsPerPhase;

    // Audio sample rates, bounded to a range an audio tool can mean. "Positive and finite" is NOT enough: an
    // absurd-but-finite rate (1e300, or Number.MIN_VALUE from a page) is no audio rate, and what each stage
    // downstream makes of one is not this tool's to find out — the meter, for one, now REFUSES a rate whose hop
    // overflows an int, where it used to reach an out-of-range lround (undefined behaviour). The probe refuses
    // such a rate itself, with its own answer.
    // THE FLOOR IS THE CORE'S (P51, core::kMinSampleRate = 8000), and it is not generous: it used to be 1000 here,
    // which sat INSIDE the band where the K-weighting shelf is past Nyquist — this probe read +3048.86 LUFS for the
    // CI fixture at 3300 Hz. Every analyzer behind the probe ABI reads the same constant. The ceiling is the
    // offline measurers' own.
    static constexpr double kMinSampleRate = felitronics::core::kMinSampleRate;
    static constexpr double kMaxSampleRate = 768000.0;

    // maxDurationSec sizes the meter's gating-block store. The default holds 4 h at any rate (1.27 MB) —
    // droppedBlocks() reads 0 for anything shorter, and a caller checks it rather than assuming. A duration whose
    // store the meter cannot represent — a block count past its int index, 3e8 s — is refused, because the meter
    // refuses it: a probe that ignored that answer would report prepared and measure nothing.
    // A REFUSAL LEAVES THE PROBE AS A FRESH ONE (law 11b), getters included: every refusal below drops the two
    // meters back to default-constructed ones, so a caller who ignored the `false` reads what a probe that was never
    // prepared reads (-120 LUFS, a zero peak, no blocks) and not the previous file. P51 made this matter: 1000..7999 Hz
    // used to RE-prepare the meters, and is a refusal now. (fc_probe's `haveResult` and fcore_measure's fresh probe
    // already kept their own callers safe; this is for everyone else.)
    bool prepare (double sampleRate, int channels, double maxDurationSec = 4.0 * 3600.0)
    {
        prepared_ = false;
        finished_ = false;
        // The compound test rejects 0, negatives, NaN and +inf, and then the absurd-but-finite rates too.
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)
            || channels < 1 || channels > felitronics::core::kMaxChannels
            || ! (maxDurationSec > 0.0) || ! std::isfinite (maxDurationSec)) return refuse();

        nc_ = channels;
        if (! lm_.prepare (sampleRate, nc_, maxDurationSec)) return refuse();
        if (! tp_.prepare (sampleRate, kChunk, nc_)) return refuse();     // prepare() also resets the running maxima
        prepared_ = true;
        return true;
    }

    bool prepared() const noexcept { return prepared_; }

    // planar[c] holds n frames for channel c. Channels beyond the prepared count are ignored; fewer than
    // prepared is honoured as-is (the meter weights only what it is given). A call NARROWER than the one before
    // it stops the channels it leaves out, and the true-peak meter drains them at that moment (law 11a, and why a
    // maximum drains rather than drops — ReferenceTruePeakMeter.h): a peak still inside their filter is measured
    // then, exactly as finish() would have measured it, and a channel that comes back starts from silence. Before
    // P62 the probe kept the history, so the reading of a narrowing stream that ends is unchanged, and one whose
    // channel RETURNS no longer replays audio from before its gap. Neither caller narrows: fcore_measure and
    // fc_probe hand every prepared channel to every call.
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

            (void) lm_.process (view, useCh, m);
            (void) tp_.process (view, useCh, m);      // prepared and 1 <= useCh <= nc_: it refuses only a null plane
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
    // has pushed silence through the filter history. Only the true-peak meter is drained — feeding the loudness
    // meter would append spurious silence to the program.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        finished_ = true;
        tp_.drain();
    }

    // --- what a cross-toolchain check compares (continuous in the input samples) ---
    int           gatingBlockCount()    const noexcept { return lm_.gatingBlockCount(); }
    std::span<const double> gatingBlockEnergies() const noexcept { return lm_.gatingBlockEnergies(); }

    // The true peak is never below the SAMPLE peak: the reconstructed signal passes through the samples by
    // construction. Enforcing that as a floor costs nothing and makes a whole class of filter-side mistake
    // (a mis-sized drain, a wrong prototype, a bad phase) impossible to hide.
    double truePeakLinear() const noexcept { return tp_.truePeakLinear(); }
    double samplePeakLinear() const noexcept { return tp_.samplePeakLinear(); }

    // --- what a human reads (derived; routed through log10, so not the bit-exactness surface) ---
    double integratedLufs() const noexcept { return lm_.integratedLufs(); }
    // The dB of truePeakLinear(), NOT of the oversampler maximum: the two differ whenever the sample peak
    // is the larger of the two (a lone impulse, anything ending abruptly), and a tool whose dB and linear
    // answers disagree is worse than one that is merely wrong.
    double truePeakDb()     const noexcept
    {
        // `det::log10`. BE ACCURATE ABOUT WHY, because the first version of this comment was not: this
        // value is NOT on the byte-diffed surface. parity.mjs prints `dbtp` only under --debug and only to
        // stderr; what CI diffs is stdout, which carries the true peak as a LINEAR bit pattern. The reason
        // to convert it is simpler and still good — it is a number a human reads from the CLI and from the
        // browser, and std::log10 differs between Apple's libm, glibc's and musl's, so the two would print
        // different digits for one file. It closes a discrepancy a user could see, not a parity hole.
        // THE 1e-9 FLOOR STAYS THIS FUNCTION'S OWN. It is NOT core::kGainToDbFloor (1e-12), and swapping
        // it for the shared one would change the number a silent file prints from -180 to -240 — a change
        // to what the tool says, smuggled in under a change to how it rounds. One thing at a time.
        const double tp = truePeakLinear();
        return 20.0 * felitronics::core::det::log10 (tp > 1e-9 ? tp : 1e-9);
    }
    int    droppedBlocks()  const noexcept { return lm_.droppedBlocks(); }
    // Forwarded so a caller can tell a MEASUREMENT from a best-effort number: non-zero means a
    // non-finite sample reached the loudness path and its 10 ms was recorded as silence.
    std::uint64_t nonFiniteSubHops() const noexcept { return lm_.nonFiniteSubHops(); }

private:
    bool refuse()
    {
        lm_ = felitronics::analysis::LoudnessMeter {};
        tp_ = felitronics::analysis::ReferenceTruePeakMeter {};
        nc_ = 0;
        return false;
    }

    felitronics::analysis::LoudnessMeter          lm_;
    felitronics::analysis::ReferenceTruePeakMeter tp_;
    int                                           nc_ = 0;
    bool                                          prepared_ = false;
    bool                                          finished_ = false;
};

// fcore::ShapeProbe — the waveform peaks and the stereo band of one file, the body shared VERBATIM by
// `fcore_measure waveform|stereo|needle|correlation` and by fc_probe_shapes_run in the wasm shim. Same reason as
// Probe above: "native and browser draw the same picture" is then a claim about one translation unit.
//
// A SEPARATE CLASS AND A SEPARATE ENTRY POINT, NOT AN OPTION ON Probe. The shapes need the file's total length
// before the first sample (every bucket and column boundary depends on it), which Probe's prepare() never asked
// for; and fc_probe_run's contract — "every call is a complete loudness measurement from scratch" — would
// otherwise have grown configuration state that outlives a call. The loudness path is untouched.
//
// THE BUILD CONTRACT above is not what makes THESE bits agree: both reductions are plain binary64 in the spec's
// order, and the one place contraction could reach — the stereo band's products — is pinned by a volatile store
// in the header itself, which survives any contraction mode.
class ShapeProbe
{
public:
    static constexpr int kChunk = Probe::kChunk;

    // THE RATE RANGE IS THE PROBE'S (the floor since P51, the ceiling since P104), checked here and not in
    // WaveformPeaks, which mirrors a page's JavaScript and only derives a decimation from the rate. Nothing outside it
    // is unstable in these two reductions; it is refused so that one ABI gives one answer about a rate —
    // `fc_probe_shapes_run` was the only run entry that took 44.1, and the only one that drew 1 MHz.
    // And a refusal RESETS THE PARTS (law 11b): `peaks()` and `stereo()` are public, and a refused call must not leave
    // them on the previous file — not their flags and not their data. Asking a part for zero channels would only clear
    // its flag (the diff-pass round: the old peaks, frame count, columns and RMS stayed readable), so each part is
    // replaced by a fresh one, as Probe replaces its meters; a refused ShapeProbe reads like a new one.
    bool prepare (double sampleRate, int channels, std::uint64_t frames, int buckets,
                  felitronics::analysis::PeakMix mix, int columns)
    {
        prepared_ = false;
        if (! (sampleRate >= Probe::kMinSampleRate && sampleRate <= Probe::kMaxSampleRate)   // NaN fails too
            || channels < 1 || channels > felitronics::core::kMaxChannels
            || ! peaks_.prepare (sampleRate, channels, frames, buckets, mix)
            || ! stereo_.prepare (channels, frames, columns))
        {
            peaks_  = felitronics::analysis::WaveformPeaks {};
            stereo_ = felitronics::analysis::StereoColumns {};
            nc_ = 0;
            return false;
        }
        nc_ = channels;
        prepared_ = true;
        return true;
    }

    // Walks the input in kChunk steps, as Probe does, so the native reader's 8192-frame reads and the shim's
    // whole-buffer call cut the stream at the same places — not that either reduction can tell.
    //
    // A CALL THAT CANNOT BE HONOURED IS REFUSED BEFORE ANYTHING MOVES (law 11): a null table or plane, or more frames
    // than the prepared length has left. Both are checked here rather than left to the two classes, because by the
    // time a later chunk reached their own checks the earlier chunks would already have been consumed.
    bool process (const float* const* planar, int channels, long long n) noexcept
    {
        if (! prepared_ || channels != nc_ || n < 0) return false;
        if (n == 0) return true;
        if (planar == nullptr) return false;
        for (int c = 0; c < nc_; ++c) if (planar[c] == nullptr) return false;
        if ((std::uint64_t) n > peaks_.totalFrames() - peaks_.framesSeen()) return false;
        const float* view[felitronics::core::kMaxChannels] {};
        for (long long off = 0; off < n; off += kChunk)
        {
            const int m = (int) std::min<long long> (kChunk, n - off);
            for (int c = 0; c < nc_; ++c) view[c] = planar[c] + off;
            if (! peaks_.process (view, nc_, m))  return false;
            if (! stereo_.process (view, nc_, m)) return false;
        }
        return true;
    }

    bool complete() const noexcept { return prepared_ && peaks_.complete() && stereo_.complete(); }
    const felitronics::analysis::WaveformPeaks& peaks()  const noexcept { return peaks_; }
    const felitronics::analysis::StereoColumns& stereo() const noexcept { return stereo_; }

private:
    felitronics::analysis::WaveformPeaks peaks_;
    felitronics::analysis::StereoColumns stereo_;
    int  nc_ = 0;
    bool prepared_ = false;
};

} // namespace fcore
