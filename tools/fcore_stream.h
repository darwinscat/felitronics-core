// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// fcore::StreamProbe — loudness and clipped runs of a stream that arrives piece by piece, read BETWEEN the
// pieces. The body shared verbatim by the `fc_stream_*` handles of the wasm shim (tools/wasm/fc_probe.cpp)
// and by `fcore_measure stream`, as fcore::ClipProbe is shared by the two `clips` roads. Nothing is measured
// here: analysis::DeterministicLoudnessMeter and analysis::ClipDetector do the measuring and are not touched.
//
// WHY THE DETERMINISTIC METER AND NOT LoudnessMeter. The system-libm meter does not give the same block
// energies natively and in wasm at every rate — `fcore_measure blocks` against parity.mjs on a 4 s
// make-fixture.mjs programme: equal at 22050 / 44100 / 48000 / 96000, different in 30 / 32 / 30 of 37 blocks at
// 8000 / 88200 / 192000. DetMath does not differ, and it is what ProgrammeReport has measured with since P80.
// fc_probe_run still sits on the system meter, so its numbers are NOT a reference for these.
//
// WHY NOT fcore::ClipProbe. Its whole contract is a FILE: it wants the length at prepare() and answers
// nothing until finish(). A stream has no length, and the point of this class is to be read before the end.
// So the detector is driven directly, and what ClipProbe's traps taught is kept where it still applies:
//   · A REFUSAL THAT CARRIED SAMPLES POISONS. Samples the caller handed over were not consumed, so every later
//     reading would describe less audio than the caller believes. Once poisoned, valid() is false for good.
//     A zero-length call carries nothing and destroys nothing.
//   · finish() IS IDEMPOTENT, and after it a call that carries samples is refused — and poisons, because the
//     frozen report no longer covers the stream as the caller sees it.
//   · RUNS READ BEFORE finish() ARE THE RUNS DECIDED SO FAR, not the runs of the audio so far: a run is decided
//     decisionDelaySamples() (20 ms) after it ends, and the last ones exist only after finish(). The run list
//     is append-only — a run, once readable, keeps its index and its fields — which is what lets a reader
//     poll it by index.
//
// THE BUILD CONTRACT of fcore_probe.h applies: ClipDetector's smooth-crest bound decides WHETHER a run exists,
// and a fused multiply-add moves it. Both roads build -ffp-contract=off; anything compared with them must too.

#include <felitronics/analysis/ClipDetector.h>
#include <felitronics/analysis/LoudnessMeter.h>

#include <cstdint>

namespace fcore
{

class StreamProbe
{
public:
    // THE METER'S STORE, one hour of gating blocks — the meter's own default and ProgrammeReportParams'
    // maxDurationSec. Blocks past it are counted in droppedBlocks() and not kept, so a stream longer than an
    // hour reads an integrated loudness of its first hour, and says so there.
    static constexpr double kMaxDurationSec = 3600.0;

    // Both instruments must accept the geometry. They do not refuse the same set on their own: the meter reads a
    // rate <= 0 or NaN as "not given" and measures at 48 kHz, the detector refuses it — so a stream opened with
    // rate 0 is refused here rather than half-measured.
    [[nodiscard]] bool prepare (double sampleRate, int channels)
    {
        prepared_ = finished_ = refused_ = false;
        channels_ = 0;
        if (channels < 1 || channels > felitronics::core::kMaxChannels) return false;
        det_.setParams (felitronics::analysis::ClipDetectorParams {});       // the default capacity, 65536 runs
        if (! meter_.prepare (sampleRate, channels, kMaxDurationSec)) return false;
        if (! det_.prepare (sampleRate, 0 /* nothing is sized by the block */, channels)) return false;
        channels_ = channels;
        prepared_ = true;
        return true;
    }

    // One piece of the stream, `channels()` planes of `n` samples. Both instruments or neither: the width and the
    // lifecycle are checked here BEFORE either moves, so neither can have consumed a piece the other refused.
    [[nodiscard]] bool process (const float* const* planar, int n) noexcept
    {
        if (! prepared_ || refused_) return false;
        if (n == 0) return ! finished_;
        if (n < 0 || finished_ || planar == nullptr) return poison();
        for (int c = 0; c < channels_; ++c) if (planar[c] == nullptr) return poison();
        // Defence in depth, unreachable today: every refusal either instrument has left is excluded above.
        if (! meter_.process (planar, channels_, n)) return poison();
        if (! det_.process (planar, channels_, n)) return poison();
        return true;
    }

    // The caller's own refusal of a piece it could not hand over — the shim's buffer checks. A piece that
    // carried samples and was not consumed is the same fault wherever it was caught.
    bool poison() noexcept { refused_ = true; return false; }

    // Decides the runs still waiting for their neighbourhood. The meter has no tail to flush: a sub-hop that
    // never closed is not a measurement, and the whole-buffer meter drops it the same way.
    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_ || refused_) return false;
        det_.finish();                                                          // idempotent itself
        finished_ = true;
        return true;
    }

    bool valid()    const noexcept { return prepared_ && ! refused_; }
    bool finished() const noexcept { return valid() && finished_; }
    int  channels() const noexcept { return valid() ? channels_ : 0; }

    // Samples per channel both instruments have consumed, read back from the detector's own clock.
    std::int64_t samples() const noexcept { return valid() ? det_.samplesProcessed() : 0; }

    // The instruments themselves, for the readings. Ungated: a caller reads them only while valid() holds.
    const felitronics::analysis::DeterministicLoudnessMeter& meter()    const noexcept { return meter_; }
    const felitronics::analysis::ClipDetector&               detector() const noexcept { return det_; }

private:
    felitronics::analysis::DeterministicLoudnessMeter meter_;
    felitronics::analysis::ClipDetector               det_;
    int  channels_ = 0;
    bool prepared_ = false, finished_ = false, refused_ = false;
};

} // namespace fcore
