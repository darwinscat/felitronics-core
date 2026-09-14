// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// fcore::ClipProbe — the clip report of one file, the body shared VERBATIM by `fcore_measure clips` and by
// fc_probe_clips_run in the wasm shim, exactly as fcore::Probe and fcore::ShapeProbe are (fcore_probe.h).
// Nothing is measured here: analysis::ClipDetector does the measuring and is not touched. What this class
// owns is the LIFECYCLE and the VALIDITY of a result, which is where all four of this mode's traps live.
//
// TRAP 1 — A RUN IS DECIDED LATE, SO A REPORT READ BEFORE finish() IS A LIE OF OMISSION.
// ClipDetector::decisionDelaySamples() is 20 ms of audio: a flat top is a candidate until every sample of its
// neighbourhood has been seen, so the last runs of a stream exist only after finish(). The detector's own
// accessors are deliberately TOTAL and answer partially while streaming (ClipDetector.h:245) — that is right
// for a meter being watched live, and wrong for a file report, where a short list reads as "clean". So the
// report path here REFUSES until finish() has succeeded: readClips() returns false and clears its output,
// and every accessor below answers zero.
//
// TRAP 2 — isFinished() IS NOT A VALIDITY FLAG, AND MEASURING SAYS SO. ClipDetector::prepare() disarms with
// `prepared_ = false` and then returns early on a refused argument (ClipDetector.h:172) — it does NOT clear
// finished_, the run list or the counters. Measured on this tree: a good run, finish(), then a refused
// prepare(0.0, …) leaves isFinished() == true, runCount() == 1 and samplePeak(0) == 0.625 — the PREVIOUS
// file's answer, behind a flag that says the measurement is done. A caller who ignores prepare()'s bool
// therefore gets a stale report that looks finished. Guarding on isFinished() alone would inherit that, so
// validity is this class's own bit: cleared at the head of every entry point, set only where a whole
// measurement has just succeeded. Same discipline, same reason, as `haveResult` / `haveShapes` in
// tools/wasm/fc_probe.cpp:104 — and IN ADDITION TO the shim's own `haveClips`, not instead of it: a call the
// ABI rejects in planarSpan() returns before prepare() is ever reached, so this flag would still be holding
// the previous file's result. Two flags, two different questions — "did this measurement complete" and "is
// there a current result to serve".
//
// TRAP 3 — A REFUSED process() STILL FINISHES, AND AN EMPTY REPORT READS AS A CLEAN FILE. ClipDetector's
// process() is [[nodiscard]] and refuses on a bad width or after finish() (ClipDetector.h:201); finish()
// then runs anyway and sets finished_. A road that dropped that bool would publish `runs 0 / complete 1 /
// peak 0` — a CERTIFICATE OF CLEANLINESS for a file nothing ever looked at, which is worse than a crash and
// worse than a short list. So a refusal here POISONS the measurement: once any call has been refused there is
// no result, whatever the caller does next. The same reason the two clocks are compared in finish(): an
// instrument must not be able to say "clean" about audio it did not see.
//
// TRAP 4 — TWO DIFFERENT TRUNCATIONS, WHICH MUST NOT BE CONFLATED. runsComplete() is about the DETECTOR's
// capacity (maxRuns): a file with more runs than that keeps counting and stores a prefix. The C ABI's output
// buffer has a capacity of its own, and a short copy there is the CALLER's business, not the file's. So
// `complete` reports only the first, the run copier returns how many runs it wrote, and the formatter refuses
// to print unless it got every stored run.
//
// WHY A SHARED CLASS AND NOT A SHARED REPORT READER. Every trap above is in setParams → prepare → process →
// finish, not in the field copying: maxRuns takes effect only at the next prepare() (ClipDetector.h:168) and
// the last runs exist only after finish(). Sharing only the reader would leave each road free to forget
// finish(), to set the capacity after prepare(), or to map the planes differently — while still calling the
// result "shared".
//
// THE BUILD CONTRACT of fcore_probe.h applies here too and is not decorative. The one expression in this
// measurement that a fused multiply-add can reach is the smooth-crest bound V = (tau + q)·stretch² + q
// (ClipDetector.h:510), and it decides WHETHER a run exists, not merely its last bit. Both roads are built
// -ffp-contract=off (tools/CMakeLists.txt for the CLI, tools/wasm/build.sh for the module); anything that
// compiles this header and compares its answers to theirs must be built the same way.

#include <felitronics/analysis/ClipDetector.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace fcore
{

// One finished clip report, as `fcore_measure clips` prints it and as the C ABI hands it out. The three
// configuration fields are carried rather than read back: the detector publishes neither its sample rate nor
// its configured capacity, and decisionDelaySamples() cannot be inverted to a rate (W = floor(sr·20/1000) is
// many-to-one).
struct ClipsReport
{
    bool         ok             = false;
    double       sampleRate     = 0.0;
    int          channels       = 0;
    std::int64_t frames         = 0;    // samples per channel that were fed
    std::int64_t maxRuns        = 0;    // the capacity the detector was prepared with
    std::int64_t decisionDelay  = 0;    // ClipDetector::decisionDelaySamples()
    std::int64_t runCount       = 0;    // every run found, stored or not
    std::int64_t storedRunCount = 0;    // the ones the list holds
    bool         complete       = false;// runsComplete(): whether the list is the whole list
    std::vector<double>                        peak;   // sample peak, one per channel
    std::vector<felitronics::analysis::ClipRun> runs;  // storedRunCount of them, in the detector's order

    void clear()
    {
        *this = ClipsReport {};
    }
};

class ClipProbe
{
public:
    // The adapter's own step. It is NOT a parameter of the measurement: law 8a says the report cannot see
    // where a stream was cut, and this class is one of the two places that could break that. It exists so the
    // native reader's 8192-frame reads and the shim's one whole-buffer call walk the detector identically.
    static constexpr int kChunk = 8192;

    // THE CAPACITY BOUND, ONE NUMBER FOR BOTH ROADS. ClipDetector::kMaxRunsLimit is 1<<24, and prepare() would
    // honour it: sizeof(ClipRun) is 40, so a 16-channel 768 kHz preparation at that capacity allocates 665 MiB
    // (measured: storageFor(768000, 16, 1<<24).bytes() = 697 641 088). This bound covers the RUN LIST only —
    // the deques and the pending queue are sized by the rate and the width, and at the same extremes they are
    // 25.3 MiB on their own, at any capacity including zero. That half is bounded by the core's own
    // kMaxSampleRate and kMaxChannels and is not narrowed here; 25 MiB is survivable, 665 is not.
    // The wasm module is built -fno-exceptions, where a failed operator new ABORTS — the page dies rather than
    // refusing — so a page's JavaScript must not be able to ask for it. 1<<20 runs is 40 MB and already far
    // past any real programme (the clipped test fixture makes about 500 runs per channel-second). It is
    // checked HERE rather than at each entry point so that the CLI and the module refuse the SAME set of
    // arguments: a capacity that is legal in one road and fatal in the other would make the parity diff
    // report a numeric failure for what is really a disagreement about the command line.
    static constexpr std::int64_t kMaxRuns = 1 << 20;

    // maxRuns is taken as int64 and checked here rather than narrowed at the call site: the C ABI hands it in
    // as a uint32_t, ClipDetectorParams::maxRuns is an `int` bounded by kMaxRunsLimit (1<<24), and an
    // unchecked cast of a value above INT_MAX is implementation-defined. 0 is legal and means "count, store
    // nothing" — prepare() accepts it (ClipDetector.h:158). Note that a capacity of 0 does not make a report
    // incomplete by itself: `complete` is `count <= capacity`, so a file with no runs at all is complete at
    // capacity 0, which is the right answer and not the obvious one.
    //
    // totalFrames IS REQUIRED, and it is the one thing this class asks for that the detector does not need.
    // ClipDetector is a pure streamer with no length in its contract, and without a length an adapter cannot
    // tell a legal call from one that runs off the end of the caller's buffer — `n` is simply believed, and a
    // wrong one is read past the allocation. Both roads know the length anyway (the CLI sizes the file before
    // reading it, the ABI is handed `frames`), and stating it here turns two duplicated "did the file deliver
    // what it was sized for" checks into one, in the class that can act on it. Same argument, same shape, as
    // fcore::ShapeProbe, which refuses a call longer than the prepared length has left (fcore_probe.h:220).
    // A 0-frame preparation is legal: an empty stream has an empty report.
    [[nodiscard]] bool prepare (double sampleRate, int channels, std::int64_t maxRuns, std::int64_t totalFrames) noexcept
    {
        valid_ = false;
        prepared_ = false;
        finished_ = false;
        refused_ = false;
        frames_ = 0;
        total_ = 0;
        if (channels < 1 || channels > felitronics::core::kMaxChannels) return false;
        if (maxRuns < 0 || maxRuns > kMaxRuns) return false;
        if (totalFrames < 0) return false;
        felitronics::analysis::ClipDetectorParams p;
        p.maxRuns = (int) maxRuns;
        det_.setParams (p);                                               // takes effect at the prepare() below
        if (! det_.prepare (sampleRate, kChunk, channels)) return false;
        sampleRate_ = sampleRate;
        channels_ = channels;
        maxRuns_ = maxRuns;
        total_ = totalFrames;
        prepared_ = true;
        return true;
    }

    // Walks the input in kChunk steps. A call that cannot be honoured is refused BEFORE anything moves
    // (law 11), because by the time a later chunk reached the detector's own checks the earlier chunks would
    // already have been consumed — and a half-eaten call is exactly the state law 8a forbids.
    //
    // THE PLANE STRIDE IS THE WHOLE POINT OF THIS LOOP. `planar[c] + off`, never `planar[c]`: the second
    // form replays the first chunk for every chunk and is invisible to any test that hands the detector
    // correct plane pointers itself.
    //
    // A NARROWER CALL IS LEGAL, and the adapter must not take that away. Law 11a: a channel the call does not
    // carry is a HOLE for those samples, which is how a channel that disappears mid-stream is expressed, and
    // the detector implements it (ClipDetector.h:209 feeds 0.0f with `fin = false` for c >= numChannels).
    // fcore::ShapeProbe requires `channels == nc_` exactly, and copying that here would have narrowed a
    // documented core behaviour into an error for no reason: `channels == 0` — every channel a hole — is
    // legal too, and then `planar` itself need not be a valid pointer. Only the first `channels` plane
    // pointers are read, which is the caller's whole obligation.
    //
    // WHAT THE POISON IS FOR, AND WHAT IT IS NOT FOR. It exists because a refused call means samples the
    // caller handed over were never consumed, and a report that then goes out describes less audio than the
    // caller believes it does. So a refusal poisons EXACTLY WHEN THE CALL CARRIED SAMPLES, and then it clears
    // validity on the spot rather than waiting for the next finish() — an accessor read in between would
    // otherwise still be serving the now-incomplete report. A zero-length call carries nothing: a malformed
    // one is answered `false`, as the detector answers it, and destroys nothing. That includes a zero-length
    // call after finish(), which is how a drain loop that has run dry ends.
    //
    // A call after finish() that DOES carry samples poisons, and that is not the same as law 12's "process()
    // refuses until reset()". The detector is right to refuse it; the adapter goes further and drops the
    // frozen report, because those samples were part of the stream as far as the caller is concerned and the
    // report no longer covers the stream. The whole premise of this header is that a caller can forget a bool.
    [[nodiscard]] bool process (const float* const* planar, int channels, long long n) noexcept
    {
        if (refused_ || ! prepared_) return false;
        auto poison = [this] { refused_ = true; valid_ = false; return false; };
        if (finished_ || channels < 0 || channels > channels_ || n < 0) return n != 0 ? poison() : false;
        if (n == 0) return true;
        // MORE FRAMES THAN THE STREAM HAS LEFT IS A REFUSAL, NOT A READ. This is the only check that can
        // catch an `n` that does not match the caller's buffer, and it has to come before the first sample
        // moves: once a chunk has been consumed the refusal would leave a half-eaten call behind, which is
        // the state law 8a forbids. It also puts the sample clock out of reach of overflow.
        if (n > total_ - frames_) return poison();
        if (channels > 0 && planar == nullptr) return poison();
        for (int c = 0; c < channels; ++c) if (planar[c] == nullptr) return poison();
        const float* view[felitronics::core::kMaxChannels] {};
        // `off += m`, not `off += kChunk`: m is min(kChunk, n - off), so off climbs to exactly n and the
        // addition cannot overflow for any n the type can hold. (The kChunk form is only safe where the
        // caller's n is bounded first, which is why the check above comes before the loop rather than after.)
        for (long long off = 0; off < n; )
        {
            const int m = (int) std::min<long long> (kChunk, n - off);
            for (int c = 0; c < channels; ++c) view[c] = planar[c] + off;
            // Defence in depth, and today it is unreachable: every guard ClipDetector::process has left —
            // an unprepared or finished detector, a negative or over-wide channel count, a negative n — is
            // excluded above. It stays because if it ever DID fire, the earlier chunks would already be
            // consumed, and a half-eaten call is the one state that must not become a report.
            if (! det_.process (view, channels, m)) return poison();
            off += m;
        }
        frames_ += n;
        return true;
    }

    // Decides the runs the neighbourhood had not yet reached, and only then is there a report to read.
    //
    // Idempotent, as the detector's own finish() is (law 12): calling it again answers the same verdict and
    // leaves the frozen report alone. The first draft cleared validity at the head of this function, which
    // made a second call — the natural thing for a caller that cannot easily tell whether it already
    // finished — silently destroy the answer. The suite caught it.
    //
    // THREE CLOCKS MUST AGREE OR THERE IS NO RESULT. total_ is the length the caller declared at prepare();
    // frames_ is what it actually handed over; samplesProcessed() is what the detector consumed, with only
    // the chunking loop above in between. A file that stopped short, a reader that dropped a block and a loop
    // that repeated a chunk are three different faults and each of them breaks one of these equalities; none
    // of them has any other symptom.
    // A REFUSAL ANYWHERE MEANS NO REPORT, even when the three clocks below happen to agree. Feeding the whole
    // stream, then making one malformed call that consumed nothing, then finishing, leaves a measurement that
    // is in fact complete — and it is still refused. The adapter does not adjudicate which refusals were
    // harmless: the caller saw a `false` and went on anyway, and that is the one situation where its own idea
    // of what it fed is known to be wrong. (The clock comparison would let this particular sequence through,
    // which is exactly why the flag is not redundant with it.)
    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_ || refused_) { valid_ = false; return false; }
        if (finished_) return valid_;      // idempotent (law 12) — and it does not discard the report it froze
        det_.finish();
        finished_ = true;
        valid_ = det_.samplesProcessed() == frames_ && frames_ == total_;
        return valid_;
    }

    // --- the report. Every one of these answers zero until finish() has succeeded on THIS preparation. ---
    bool valid() const noexcept { return valid_; }
    double       sampleRate()     const noexcept { return valid_ ? sampleRate_ : 0.0; }
    int          channels()       const noexcept { return valid_ ? channels_ : 0; }
    std::int64_t frames()         const noexcept { return valid_ ? frames_ : 0; }
    std::int64_t maxRuns()        const noexcept { return valid_ ? maxRuns_ : 0; }
    std::int64_t decisionDelay()  const noexcept { return valid_ ? det_.decisionDelaySamples() : 0; }
    std::int64_t runCount()       const noexcept { return valid_ ? det_.runCount() : 0; }
    std::int64_t storedRunCount() const noexcept { return valid_ ? det_.storedRunCount() : 0; }
    bool         complete()       const noexcept { return valid_ && det_.runsComplete(); }
    double       peak (int c)     const noexcept { return valid_ ? det_.samplePeak (c) : 0.0; }
    felitronics::analysis::ClipRun run (std::int64_t i) const noexcept
    {
        return valid_ ? det_.run (i) : felitronics::analysis::ClipRun {};
    }

private:
    felitronics::analysis::ClipDetector det_;
    double       sampleRate_ = 0.0;
    int          channels_   = 0;
    std::int64_t frames_     = 0;
    std::int64_t total_      = 0;
    std::int64_t maxRuns_    = 0;
    bool prepared_ = false, finished_ = false, valid_ = false, refused_ = false;
};

// Copies a finished probe into a report value. Returns false — and CLEARS `out`, so a reused report cannot
// carry the previous file's answer past a refusal — when the probe has no result.
inline bool readClips (const ClipProbe& p, ClipsReport& out)
{
    out.clear();
    if (! p.valid()) return false;
    out.sampleRate = p.sampleRate();
    out.channels = p.channels();
    out.frames = p.frames();
    out.maxRuns = p.maxRuns();
    out.decisionDelay = p.decisionDelay();
    out.runCount = p.runCount();
    out.storedRunCount = p.storedRunCount();
    out.complete = p.complete();
    out.peak.reserve ((std::size_t) out.channels);
    for (int c = 0; c < out.channels; ++c) out.peak.push_back (p.peak (c));
    out.runs.reserve ((std::size_t) out.storedRunCount);
    for (std::int64_t i = 0; i < out.storedRunCount; ++i) out.runs.push_back (p.run (i));
    out.ok = true;
    return true;
}

} // namespace fcore
