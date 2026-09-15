// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fc_probe — the P0 spike's C export: the thinnest possible skin over fcore::Probe, which is the SAME
// translation unit the native reference CLI runs (tools/fcore_probe.h). Nothing is measured here; this file
// only unpacks an ABI.
//
// ABI SHAPE. Planar float32, one pointer, channel c at `planar + c*frames`. That satisfies the ratified
// "PCM not files, planar float32 in heap" boundary while avoiding `const float* const*` across the wasm
// boundary — a pointer-to-pointer would mean building a table of i32 offsets in the heap and exporting
// HEAPU32 to write it, for no gain. It also maps 1:1 onto what JS already has: AudioBuffer.getChannelData(c)
// is planar, so the page does one HEAPF32.set() per channel and no de-interleave loop.
//
// `frames` is uint32_t, not int: wasm32 is a 32-bit target and a signed frame count invites an overflow that
// simply cannot happen on the 64-bit machine this core was written and tested on.
//
// WHAT THE PARITY CHECK READS. fc_probe_lufs / fc_probe_dbtp are the two numbers P0 asks for, but neither is
// the comparison surface: the integrated LUFS is discontinuous at the BS.1770 gates and the dB form of the
// true peak is routed through log10. The surface is fc_probe_block_energies() (pre-gate, continuous) plus
// fc_probe_tp_linear(). See tools/fcore_probe.h for why.

#include "fcore_clips.h"
#include "fcore_probe.h"

#include <felitronics/analysis/BandBursts.h>
#include <felitronics/analysis/HumDetector.h>
#include <felitronics/analysis/LowEnd.h>
#include <felitronics/analysis/SourceForensics.h>
#include <felitronics/analysis/ProgrammeReport.h>

#include <cstdint>
#include <cstring>
#include <vector>

#if defined(__EMSCRIPTEN__)
  #include <emscripten/emscripten.h>
  #include <emscripten/heap.h>              // emscripten_get_heap_size — NOT declared by emscripten.h
  #define FC_EXPORT extern "C" EMSCRIPTEN_KEEPALIVE
#else
  #define FC_EXPORT extern "C"
#endif

namespace
{
    // Function-local static rather than a file-scope global: no static-initialisation order to reason about,
    // and the spike is explicitly single-threaded (law 1 — the whole point of criterion 3).
    fcore::Probe& probe()
    {
        static fcore::Probe p;
        return p;
    }

    // Does [p, p+bytes) lie inside the wasm linear memory? A pointer can be aligned, and its length can fit a
    // 32-bit address space, and the span can still run off the end of the heap — an aligned pointer four
    // bytes below the top with frames=2 passes every other check here and then traps with "memory access out
    // of bounds". This cannot prove the caller actually owns the span (no ABI of this shape can), but it does
    // turn "the module dies" into "the call is refused".
    bool inHeap (const void* p, std::uint64_t bytes)
    {
#if defined(__EMSCRIPTEN__)
        const std::uint64_t base = (std::uint64_t) reinterpret_cast<std::uintptr_t> (p);
        if (base == 0) return false;
        const std::uint64_t end = base + bytes;
        if (end < base) return false;                                     // wrapped
        return end <= (std::uint64_t) emscripten_get_heap_size();
#else
        (void) p; (void) bytes;
        return true;                                                      // native: no linear memory to bound
#endif
    }

    // The planar input span: non-null, non-empty, a width the core has, 4-byte aligned, and inside the heap.
    bool planarSpan (const float* planar, std::uint32_t frames, std::uint32_t channels)
    {
        if (planar == nullptr || frames == 0) return false;
        if (channels < 1 || channels > (std::uint32_t) felitronics::core::kMaxChannels) return false;
        if ((reinterpret_cast<std::uintptr_t> (planar) & 0x3u) != 0) return false;   // a misaligned float* reads
                                                                                     // garbage in a release build
                                                                                     // and only traps under SAFE_HEAP
        // frames*channels must address real memory: on wasm32 the product is what a caller malloc'd, so a
        // wrapped one would hand us a window onto someone else's heap.
        const std::uint64_t bytes = (std::uint64_t) frames * (std::uint64_t) channels * sizeof (float);
        if (bytes > (std::uint64_t) 0xFFFFFFFFu) return false;
        return inHeap (planar, bytes);
    }

    // Ptr/size validation the core cannot do for us: everything below arrives from JS.
    bool viable (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
    {
        if (! planarSpan (planar, frames, channels)) return false;
        // rejects 0 / negative / NaN / +inf and the absurd-but-finite rates (see Probe::kMinSampleRate)
        return probe().prepare (sampleRate, (int) channels);
    }

    // An output span of `count` elements of `align` bytes each: non-null, aligned, inside the heap.
    bool outSpan (const void* out, std::uint32_t count, std::uint32_t align)
    {
        if (out == nullptr) return false;
        if ((reinterpret_cast<std::uintptr_t> (out) & (align - 1u)) != 0) return false;
        return inHeap (out, (std::uint64_t) count * align);
    }

    // Whether the getters have a result to report. Without this the contract would be an accident of where
    // Probe::prepare() happens to return: a bad sample rate is caught before the meter is touched, so the
    // getters would go on serving the PREVIOUS file's numbers to a caller who ignored the return value. An
    // ABI fed by a page's JavaScript should not have footguns that subtle.
    bool haveResult = false;

    bool run (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
    {
        haveResult = false;
        if (! viable (planar, frames, channels, sampleRate)) return false;
        const float* view[felitronics::core::kMaxChannels] {};
        for (std::uint32_t c = 0; c < channels; ++c) view[c] = planar + (std::size_t) c * (std::size_t) frames;
        probe().process (view, (int) channels, (long long) frames);
        probe().finish();          // drain the FIR, or a peak in the final samples goes unmeasured
        haveResult = true;
        return true;
    }
}

// Runs the measurement and leaves the result readable by the getters below. Returns 1 on success, 0 if the
// arguments were rejected.
//
// A REJECTED CALL CLEARS THE PREVIOUS RESULT: the getters read zero after a failure, never the last good run.
// That is enforced by a flag rather than left to fall out of where Probe::prepare() returns — a bad sample
// rate is caught before the meter is touched, so without the flag a caller who ignored the return value would
// be served the previous file's numbers. Every exported entry point below is a complete measurement from
// scratch: calling fc_probe_lufs() after fc_probe_run() re-runs the whole thing, and the getters then
// describe THAT run.
FC_EXPORT int fc_probe_run (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
{
    return run (planar, frames, channels, sampleRate) ? 1 : 0;
}

// The two numbers P0 names. Each is a full run; -120.0 / the dB floor come back if the arguments are refused.
FC_EXPORT double fc_probe_lufs (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
{
    if (! run (planar, frames, channels, sampleRate)) return -120.0;
    return probe().integratedLufs();
}

FC_EXPORT double fc_probe_dbtp (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
{
    if (! run (planar, frames, channels, sampleRate)) return 20.0 * std::log10 (1e-9);
    return probe().truePeakDb();
}

// --- the parity surface, read from the last successful run ---

FC_EXPORT double        fc_probe_tp_linear   (void) { return haveResult ? probe().truePeakLinear() : 0.0; }
FC_EXPORT double        fc_probe_sample_peak (void) { return haveResult ? probe().samplePeakLinear() : 0.0; }
FC_EXPORT std::uint32_t fc_probe_block_count  (void) { return haveResult ? (std::uint32_t) probe().gatingBlockCount() : 0u; }
FC_EXPORT int           fc_probe_dropped      (void) { return haveResult ? probe().droppedBlocks() : 0; }

// Copies min(count, cap) pre-gate block energies into `out` and returns how many were written. Doubles, so
// the page reads them through HEAPF64.
FC_EXPORT std::uint32_t fc_probe_block_energies (double* out, std::uint32_t cap)
{
    if (out == nullptr || ! haveResult) return 0;
    // ALIGNMENT, on the OUTPUT as well as the input. This file has refused a misaligned `const float*`
    // on the way in since P0 and checked only the ADDRESS RANGE on the way out — the same rule applied
    // in one direction. What the asymmetry costs depends on how the page reads the span back, and both
    // ways are bad: this repo's own consumers use `HEAPF64.subarray(ptr >>> 3, …)`
    // (`tools/wasm/parity.mjs`, `probe.html`), where a misaligned `ptr` TRUNCATES to a different index
    // and the caller silently reads someone else's doubles, while the byte-offset form
    // `new Float64Array(HEAPF64.buffer, ptr, n)` throws instead. A refusal replaces both with a zero the
    // caller can test. (Found while fixing the same class in fc_master.cpp, which had it on every
    // scalar out-parameter. NB our two consumers allocate with `_malloc` and are therefore aligned
    // today — this closes the hole rather than a live failure.)
    if ((reinterpret_cast<std::uintptr_t> (out) & 0x7u) != 0) return 0;
    const std::uint32_t n = (std::uint32_t) probe().gatingBlockCount();
    const std::uint32_t m = n < cap ? n : cap;
    if (! inHeap (out, (std::uint64_t) m * sizeof (double))) return 0;   // the output span must fit too
    if (m > 0) std::memcpy (out, probe().gatingBlockEnergies().data(), (std::size_t) m * sizeof (double));
    return m;
}

// --- the waveform peaks and the stereo band (P59a) ---
//
// NAMES THAT CANNOT BE READ AS ANOTHER NUMBER OF THIS ABI. `fc_probe_waveform_*` and not `fc_probe_peaks`: this file
// already reports two peaks, `fc_probe_sample_peak` and `fc_probe_tp_linear` (the true peak of the 128-tap REFERENCE
// filter, analysis::ReferenceTruePeakMeter — the one TargetLoudnessSolver aims a delivered ceiling with — and not
// analysis::TruePeakMeter's 48-tap spec filter: different filters, different numbers, the gap pinned by
// felitronics_truepeak_instrument_gap_tests), and a waveform bucket is neither: a
// box-averaged max-abs, not above the sample peak except by the rounding of a box mean. And
// `fc_probe_stereo_rms` and not `_loud`, the spec's name: it is an RMS, and it sits in the same ABI as
// `fc_probe_lufs`.
//
// ONE DEFINITION, TWO ROADS. The page draws these with its own JavaScript while the sidecars come from a server-side
// generator with a different definition; this is the road both are meant to take instead — fcore::ShapeProbe, the class
// `fcore_measure waveform|stereo|needle` runs natively. (Moving the site onto it is the site's work.) What is computed and why each output has the type it has:
// modules/analysis/include/felitronics/analysis/WaveformPeaks.h and StereoColumns.h.
//
// A SEPARATE RUN, NOT fc_probe_run WITH OPTIONS. The shapes need the file's length and three parameters before
// the first sample, and fc_probe_run's promise is a complete loudness measurement per call with no state carried
// between calls; a setter feeding it would be exactly that state. The loudness result is not touched by a shapes
// run, and a shapes result is not touched by a loudness run. Same rule as haveResult: a REJECTED shapes run
// clears the previous shapes result, so the getters below read zero rather than the last good file.
namespace
{
    fcore::ShapeProbe& shapes()
    {
        static fcore::ShapeProbe s;
        return s;
    }
    bool haveShapes = false;

    template <typename T, typename Src>
    std::uint32_t copyOut (T* out, std::uint32_t cap, std::uint32_t n, Src&& at)
    {
        const std::uint32_t m = n < cap ? n : cap;
        if (m == 0 || ! outSpan (out, m, (std::uint32_t) sizeof (T))) return 0;
        for (std::uint32_t i = 0; i < m; ++i) out[i] = at (i);
        return m;
    }
}

// `mix` is the PeakMix code: 0 'avr' · 1 'L' · 2 'R' · 3 'max'. Any other code is refused, not clamped. Returns 1
// when the whole file was reduced, 0 when the arguments were refused.
FC_EXPORT int fc_probe_shapes_run (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate,
                                   std::uint32_t buckets, std::int32_t mix, std::uint32_t columns)
{
    haveShapes = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    if (mix < 0 || mix > 3) return 0;
    if (buckets > 0x7FFFFFFFu || columns > 0x7FFFFFFFu) return 0;   // the core's own bounds are far below; this is
                                                                     // only the narrowing to its `int`
    auto& s = shapes();
    if (! s.prepare (sampleRate, (int) channels, frames, (int) buckets,
                     (felitronics::analysis::PeakMix) mix, (int) columns)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t c = 0; c < channels; ++c) view[c] = planar + (std::size_t) c * (std::size_t) frames;
    if (! s.process (view, (int) channels, (long long) frames) || ! s.complete()) return 0;
    haveShapes = true;
    return 1;
}

FC_EXPORT std::uint32_t fc_probe_waveform_count      (void) { return haveShapes ? (std::uint32_t) shapes().peaks().buckets() : 0u; }
FC_EXPORT std::uint32_t fc_probe_waveform_emitted    (void) { return haveShapes ? (std::uint32_t) shapes().peaks().bucketsEmitted() : 0u; }
FC_EXPORT std::uint32_t fc_probe_waveform_decimation (void) { return haveShapes ? (std::uint32_t) shapes().peaks().decimation() : 0u; }

// Copies min(count, cap) peaks and returns how many were written — 0 for a null, misaligned or out-of-heap
// buffer. `fc_probe_waveform_peaks` is `computePeaksFromBuffer`'s double form, read through HEAPF64;
// `fc_probe_waveform_peaks_f32` is `peaksFromWav`'s Float32Array form, read through HEAPF32.
FC_EXPORT std::uint32_t fc_probe_waveform_peaks (double* out, std::uint32_t cap)
{
    if (! haveShapes) return 0;
    const auto p = shapes().peaks().peaks();
    return copyOut (out, cap, (std::uint32_t) p.size(), [&] (std::uint32_t i) { return p[i]; });
}

FC_EXPORT std::uint32_t fc_probe_waveform_peaks_f32 (float* out, std::uint32_t cap)
{
    if (! haveShapes) return 0;
    const auto& w = shapes().peaks();
    return copyOut (out, cap, (std::uint32_t) w.buckets(), [&] (std::uint32_t i) { return w.peakAsFloat32 ((int) i); });
}

FC_EXPORT std::uint32_t fc_probe_stereo_cols     (void) { return haveShapes ? (std::uint32_t) shapes().stereo().columns() : 0u; }
FC_EXPORT int           fc_probe_stereo_is_mono  (void) { return haveShapes && shapes().stereo().isMono() ? 1 : 0; }
FC_EXPORT double        fc_probe_stereo_max_rms (void) { return haveShapes ? shapes().stereo().maxRms() : 0.0; }

FC_EXPORT std::uint32_t fc_probe_stereo_width (float* out, std::uint32_t cap)
{
    if (! haveShapes) return 0;
    const auto v = shapes().stereo().width();
    return copyOut (out, cap, (std::uint32_t) v.size(), [&] (std::uint32_t i) { return v[i]; });
}

FC_EXPORT std::uint32_t fc_probe_stereo_corr (float* out, std::uint32_t cap)
{
    if (! haveShapes) return 0;
    const auto v = shapes().stereo().correlation();
    return copyOut (out, cap, (std::uint32_t) v.size(), [&] (std::uint32_t i) { return v[i]; });
}

FC_EXPORT std::uint32_t fc_probe_stereo_rms (float* out, std::uint32_t cap)
{
    if (! haveShapes) return 0;
    const auto v = shapes().stereo().rms();
    return copyOut (out, cap, (std::uint32_t) v.size(), [&] (std::uint32_t i) { return v[i]; });
}

// The playhead needle over [from, to) of a planar buffer — `correlationOf` and `widthOf`, plus the stretch's
// RMS — written as three doubles into out3[0..2]. Stateless: it neither reads nor clears a shapes result.
// Channel 0 is L; channel 1 is R, or channel 0 again for a mono buffer. Returns 1, or 0 with nothing written.
FC_EXPORT int fc_probe_needle (const float* planar, std::uint32_t frames, std::uint32_t channels,
                               std::uint32_t from, std::uint32_t to, double* out3)
{
    if (! planarSpan (planar, frames, channels)) return 0;
    if (! outSpan (out3, 3u, 8u)) return 0;
    const float* L = planar;
    const float* R = channels > 1 ? planar + (std::size_t) frames : planar;
    felitronics::analysis::StereoColumns::Needle nd;
    if (! felitronics::analysis::StereoColumns::needle (L, R, frames, from, to, nd)) return 0;
    out3[0] = nd.correlation;
    out3[1] = nd.width;
    out3[2] = nd.rms;
    return 1;
}

// --- the clipped runs (P71) ---
//
// analysis::ClipDetector, through fcore::ClipProbe — the same adapter `fcore_measure clips` drives, so the
// two roads share a lifecycle and not merely a report reader. See tools/fcore_clips.h for the four traps that
// class exists to close; the two this file is responsible for are the last of them.
//
// ONE RUN, THEN GETTERS, like fc_probe_shapes_run and for the same reason: the capacity of the run list has
// to be chosen before the first sample (maxRuns takes effect at prepare()), and a setter feeding fc_probe_run
// would be exactly the state that entry point promises not to carry between calls. A REJECTED RUN CLEARS THE
// PREVIOUS RESULT — `haveClips` is cleared before any validation, the `haveResult` / `haveShapes` discipline
// of this file — so the getters answer zero rather than the last good file's runs.
//
// EVERY HEADER FIELD IS READ BACK FROM THE MEASUREMENT, NOT RECOMPUTED BY THE CALLER. fc_probe_clips_samples
// and _channels and _rate look redundant next to the arguments the caller just passed in; they are the
// opposite. If this shim ever fed the detector half the buffer, or the wrong width, a report whose header
// came from the harness's own arithmetic would still diff clean against the native tool — the numbers would
// agree because neither side asked the instrument. _delay has a second reason: it keeps JavaScript from
// re-deriving floor(sr*20/1000), which would be a second copy of a definition that already exists in C++.
//
// THE CAPACITY BOUND IS NOT OPTIONAL HERE. ClipDetector::kMaxRunsLimit (1<<24) would allocate 665 MB at
// 16 channels, and this module is built -fno-exceptions, where a failed allocation aborts the page instead of
// refusing. fcore::ClipProbe::kMaxRuns (1<<20) is the bound both roads apply.
namespace
{
    fcore::ClipProbe& clips()
    {
        static fcore::ClipProbe c;
        return c;
    }
    bool haveClips = false;

    // start, length, level, channel, sign, evidence — the fields of analysis::ClipRun, in the order the text
    // format prints them. Doubles throughout: `start` and `length` are int64 in the core, but planarSpan()
    // caps frames*channels*4 at 4 GiB, so a position in this ABI is below 2^30 and exact in a binary64 — while
    // an i64 return would need -sWASM_BIGINT (which this module does not build with) and would not match the
    // export whitelist's grep in tools/wasm/build.sh.
    constexpr std::uint32_t kClipRunStride = 6;
}

// Measures one planar buffer end to end — prepare, the whole file, finish — and leaves the result readable by
// the getters below. Returns 1, or 0 with every getter of this result cleared (fc_probe_clips_stride excepted:
// it is a property of the FORMAT, not of a measurement, and answers 6 always). `maxRuns` is the run list's
// capacity: 0 is legal and means "count them, store none" — and note that a file with no runs is then still
// COMPLETE, because completeness is `count <= capacity`.
//
// The capacity bound is not re-checked here. ClipProbe::prepare() applies it, `maxRuns` widens from uint32 to
// int64 without loss on the way in, and one rule in one place is the whole reason that constant lives in the
// shared header — the CLI and this entry point must refuse the same set or the parity diff reports a
// measurement failure for a disagreement about a command line.
FC_EXPORT int fc_probe_clips_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                  double sampleRate, std::uint32_t maxRuns)
{
    haveClips = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& c = clips();
    if (! c.prepare (sampleRate, (int) channels, (std::int64_t) maxRuns, (std::int64_t) frames)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! c.process (view, (int) channels, (long long) frames) || ! c.finish()) return 0;
    haveClips = true;
    return 1;
}

FC_EXPORT double        fc_probe_clips_rate     (void) { return haveClips ? clips().sampleRate() : 0.0; }
FC_EXPORT std::uint32_t fc_probe_clips_channels (void) { return haveClips ? (std::uint32_t) clips().channels() : 0u; }
FC_EXPORT std::uint32_t fc_probe_clips_samples  (void) { return haveClips ? (std::uint32_t) clips().frames() : 0u; }
FC_EXPORT std::uint32_t fc_probe_clips_max_runs (void) { return haveClips ? (std::uint32_t) clips().maxRuns() : 0u; }
FC_EXPORT std::uint32_t fc_probe_clips_delay    (void) { return haveClips ? (std::uint32_t) clips().decisionDelay() : 0u; }
FC_EXPORT std::uint32_t fc_probe_clips_count    (void) { return haveClips ? (std::uint32_t) clips().runCount() : 0u; }
FC_EXPORT std::uint32_t fc_probe_clips_stored   (void) { return haveClips ? (std::uint32_t) clips().storedRunCount() : 0u; }
FC_EXPORT int           fc_probe_clips_complete (void) { return haveClips && clips().complete() ? 1 : 0; }
FC_EXPORT std::uint32_t fc_probe_clips_stride   (void) { return kClipRunStride; }

// The per-channel sample peak, min(channels, cap) of them, read through HEAPF64.
FC_EXPORT std::uint32_t fc_probe_clips_peaks (double* out, std::uint32_t cap)
{
    if (! haveClips) return 0;
    const auto& c = clips();
    return copyOut (out, cap, (std::uint32_t) c.channels(), [&] (std::uint32_t i) { return c.peak ((int) i); });
}

// The runs, as many whole ones as fit, kClipRunStride doubles each.
//
// `cap` IS IN DOUBLES, LIKE EVERY OTHER COPIER IN THIS FILE, AND THE RETURN IS IN RUNS. That asymmetry is
// deliberate and it is the safe way round. Seven exports here share the shape `(T* out, uint32_t cap)` —
// fc_probe_block_energies, the two waveform copiers, the three stereo copiers, fc_probe_clips_peaks — and in
// every one of them `cap` is the number of ELEMENTS the buffer holds. A caller who reads this file, follows
// that convention and writes `p = _malloc(stored * 8); fc_probe_clips_runs(p, stored)` must not be handed a
// six-fold heap overwrite — and neither outSpan() nor inHeap() could see it, because they bound the linear
// memory and not the allocation. With `cap` in doubles that caller gets floor(stored/6) runs: too few, which
// is visible in its own output, instead of memory corruption that is not. (The crew's review round found this
// as a live hazard in the first shape of this entry point, where `cap` was in runs.)
//
// The RETURN is in runs because that is the number a reader needs, and because the two truncations must stay
// apart: a short buffer is the caller's business, which it can see and fix, while fc_probe_clips_complete()
// is the DETECTOR's capacity running out, which is a property of the file. A reader that conflated them would
// call a file incomplete because it passed a small buffer.
FC_EXPORT std::uint32_t fc_probe_clips_runs (double* out, std::uint32_t cap)
{
    if (! haveClips) return 0;
    const auto& c = clips();
    const std::uint32_t stored = (std::uint32_t) c.storedRunCount();
    const std::uint32_t fits = cap / kClipRunStride;               // whole runs only; a partial one is not written
    const std::uint32_t m = fits < stored ? fits : stored;
    if (m == 0) return 0;
    const std::uint32_t n = m * kClipRunStride;                   // m <= kMaxRuns (2^20), so this cannot wrap
    const std::uint32_t wrote = copyOut (out, n, n, [&] (std::uint32_t i)
    {
        const felitronics::analysis::ClipRun r = c.run ((std::int64_t) (i / kClipRunStride));
        switch (i % kClipRunStride)
        {
            case 0:  return (double) r.start;
            case 1:  return (double) r.length;
            case 2:  return r.level;
            case 3:  return (double) r.channel;
            case 4:  return (double) r.sign;
            default: return (double) (int) r.evidence;
        }
    });
    return wrote == n ? m : 0u;                                   // a refused span writes nothing and says so
}

// Build identity, so a mismatched artifact is obvious in a report rather than a mystery.
FC_EXPORT std::uint32_t fc_probe_os_factor       (void) { return (std::uint32_t) fcore::Probe::kOsFactor; }
FC_EXPORT std::uint32_t fc_probe_os_taps         (void) { return (std::uint32_t) fcore::Probe::kOsTapsPerPhase; }
FC_EXPORT std::uint32_t fc_probe_chunk           (void) { return (std::uint32_t) fcore::Probe::kChunk; }
FC_EXPORT std::uint32_t fc_probe_sizeof_longdouble (void) { return (std::uint32_t) sizeof (long double); }

//==============================================================================
// analysis::ProgrammeReport (P72) through the ABI.
//
// WHY THIS ONE IS SHAPED DIFFERENTLY FROM `clips`. ClipDetector publishes eight scalars and one list, so
// one export per field was honest. ProgrammeReport publishes on the order of a hundred NAMED values, each
// carrying a validity and a reason, and the report's visitor is deliberately THE ONE ENUMERATION of them
// ("a field added here is automatically compared and automatically printed; a hand-written list in either
// place would silently stop covering the new field, which is how a gate quietly stops being a gate").
// One export per field would rebuild exactly that hand-written list on this side, in C++ and again in
// JavaScript. So the rows and their NAMES both come out of the same walk: add a field to the visitor and
// it appears in the CLI, in the module and in the parity diff, with nothing to renumber.
//
// The names are a NUL-separated blob in visitor order — counts first, then values — and the rows are
// doubles in the same order. A row carries no name of its own: its name is its position. That is only
// safe because ONE walk produces both, which is the whole reason it is shaped this way.
namespace
{
    felitronics::analysis::ProgrammeReport& programme()
    {
        static felitronics::analysis::ProgrammeReport p;
        return p;
    }
    bool haveReport = false;

    constexpr std::uint32_t kReportCountStride = 2;   // channel, value
    constexpr std::uint32_t kReportValueStride = 4;   // channel, valid, reason, value
}

FC_EXPORT int fc_probe_report_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    haveReport = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& p = programme();
    if (! p.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! p.process (view, (int) channels, (int) frames)) return 0;
    p.finish();
    haveReport = true;
    return 1;
}

FC_EXPORT double fc_probe_report_samples (void)
{
    return haveReport ? (double) programme().report().totalSamples : 0.0;
}

FC_EXPORT std::uint32_t fc_probe_report_count_rows (void)
{
    if (! haveReport) return 0u;
    std::uint32_t n = 0;
    programme().report().visitCounts ([&n] (const char*, int, std::int64_t) { ++n; });
    return n;
}

FC_EXPORT std::uint32_t fc_probe_report_value_rows (void)
{
    if (! haveReport) return 0u;
    std::uint32_t n = 0;
    programme().report().visitValues ([&n] (const char*, int, const felitronics::analysis::ProgrammeValue&) { ++n; });
    return n;
}

FC_EXPORT std::uint32_t fc_probe_report_stride_counts (void) { return kReportCountStride; }
FC_EXPORT std::uint32_t fc_probe_report_stride_values (void) { return kReportValueStride; }

// The field names, NUL-separated, counts first then values, in the visitor's own order. Returns the bytes
// written, or the bytes REQUIRED when `cap` is 0 — so a caller sizes its buffer from the module rather
// than from a number it made up.
FC_EXPORT std::uint32_t fc_probe_report_names (char* out, std::uint32_t cap)
{
    if (! haveReport) return 0u;
    std::uint32_t need = 0;
    const auto& R = programme().report();
    auto measure = [&need] (const char* nm) { std::uint32_t k = 0; while (nm[k] != '\0') ++k; need += k + 1u; };
    R.visitCounts ([&] (const char* nm, int, std::int64_t) { measure (nm); });
    R.visitValues ([&] (const char* nm, int, const felitronics::analysis::ProgrammeValue&) { measure (nm); });
    if (cap == 0u) return need;
    if (out == nullptr || ! outSpan (out, cap, 1) || cap < need) return 0u;
    std::uint32_t at = 0;
    auto emit = [&] (const char* nm) { std::uint32_t k = 0; while (nm[k] != '\0') out[at++] = nm[k++]; out[at++] = '\0'; };
    R.visitCounts ([&] (const char* nm, int, std::int64_t) { emit (nm); });
    R.visitValues ([&] (const char* nm, int, const felitronics::analysis::ProgrammeValue&) { emit (nm); });
    return at;
}

// Rows into a caller-owned buffer, capacity in DOUBLES and the return in ROWS — the fc_probe_clips_runs
// rule. A short capacity writes the prefix that fits and says how many rows that was; it never
// half-writes a row, and a truncated read is the CALLER's business, never confused with the report's own
// completeness.
FC_EXPORT std::uint32_t fc_probe_report_counts (double* out, std::uint32_t cap)
{
    if (! haveReport || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const std::uint32_t room = cap / kReportCountStride;
    std::uint32_t at = 0;
    programme().report().visitCounts ([&] (const char*, int ch, std::int64_t v)
    {
        if (at >= room) return;
        out[at * kReportCountStride + 0] = (double) ch;
        out[at * kReportCountStride + 1] = (double) v;
        ++at;
    });
    return at;
}

FC_EXPORT std::uint32_t fc_probe_report_values (double* out, std::uint32_t cap)
{
    if (! haveReport || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const std::uint32_t room = cap / kReportValueStride;
    std::uint32_t at = 0;
    programme().report().visitValues ([&] (const char*, int ch, const felitronics::analysis::ProgrammeValue& v)
    {
        if (at >= room) return;
        out[at * kReportValueStride + 0] = (double) ch;
        out[at * kReportValueStride + 1] = v.valid ? 1.0 : 0.0;
        out[at * kReportValueStride + 2] = (double) (int) v.reason;
        out[at * kReportValueStride + 3] = v.value;      // the bits, untouched — invalid is a canonical +0.0
        ++at;
    });
    return at;
}

//==============================================================================
// analysis::BandBursts (P76) through the ABI.
//
// Unlike `report`, this mode's text is POSITIONAL — the CLI prints fixed columns, not named fields — so
// there is no name table to carry and the ABI is the plain shape: one scalar block in the order the CLI
// prints it, then the three variable-length lists (per channel, events, and the two sparse histograms),
// each into a caller-owned buffer whose capacity is mandatory.
//
// `enterDb` and `exitDb` are read from a default-constructed BandBurstsParams here because the detector
// does not publish them and the CLI does exactly the same. That is the one number on this road not read
// back out of the measurement; it is the documented default on both sides, so the diff still covers it,
// but it is named rather than hidden.
namespace
{
    felitronics::analysis::BandBursts& bursts()
    {
        static felitronics::analysis::BandBursts d;
        return d;
    }
    bool haveBursts = false;

    constexpr std::uint32_t kBurstsScalars    = 31;
    constexpr std::uint32_t kBurstsChanStride = 4;
    constexpr std::uint32_t kBurstsEvtStride  = 12;
    constexpr std::uint32_t kBurstsBinStride  = 2;
}

FC_EXPORT int fc_probe_bursts_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    haveBursts = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = bursts();
    d.setParams (felitronics::analysis::BandBurstsParams {});
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! d.process (view, (int) channels, (int) frames)) return 0;
    d.finish();
    haveBursts = true;
    return 1;
}

FC_EXPORT std::uint32_t fc_probe_bursts_scalars_len (void) { return kBurstsScalars; }
FC_EXPORT std::uint32_t fc_probe_bursts_chan_stride (void) { return kBurstsChanStride; }
FC_EXPORT std::uint32_t fc_probe_bursts_evt_stride  (void) { return kBurstsEvtStride; }
FC_EXPORT std::uint32_t fc_probe_bursts_bin_stride  (void) { return kBurstsBinStride; }

// The scalars, IN THE ORDER THE CLI PRINTS THEM. The order is the contract; the JavaScript half reads
// them by index and must not be edited without editing this.
FC_EXPORT std::uint32_t fc_probe_bursts_scalars (double* out, std::uint32_t cap)
{
    if (! haveBursts || out == nullptr || cap < kBurstsScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = bursts();
    const felitronics::analysis::BandBurstsParams bp {};
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();
    out[i++] = (double) d.channels();
    out[i++] = (double) d.hopSamples();
    out[i++] = (double) d.baselineHops();
    out[i++] = d.bandLowHz();
    out[i++] = d.bandHighHz();
    out[i++] = bp.enterDb;
    out[i++] = bp.exitDb;
    out[i++] = (double) fcore::Probe::kChunk;
    out[i++] = (double) d.samplesProcessed();
    out[i++] = (double) d.hopCount();
    out[i++] = (double) d.eligibleHops();
    out[i++] = (double) d.zeroBaselineHops();
    out[i++] = (double) d.burstHops();
    out[i++] = (double) d.damagedHops();
    out[i++] = (double) d.overflowSamples();
    out[i++] = (double) d.firstNonFiniteAt();
    out[i++] = (double) d.tailPartialSamples();
    out[i++] = d.tailPartialEnergy();
    out[i++] = d.eventsValid() ? 1.0 : 0.0;
    out[i++] = (double) (int) d.eventsInvalidReason();
    out[i++] = d.programmeEnergyValid() ? 1.0 : 0.0;
    out[i++] = (double) (int) d.programmeEnergyInvalidReason();
    out[i++] = (double) d.eventCount();
    out[i++] = (double) d.storedEventCount();
    out[i++] = d.eventsComplete() ? 1.0 : 0.0;
    out[i++] = (double) d.onsetCount();
    out[i++] = (double) d.intervalCount();
    out[i++] = (double) d.intervalOverflow();
    out[i++] = (double) d.modalIntervalHops();
    out[i++] = (double) d.modalIntervalMass();
    return i;
}

FC_EXPORT std::uint32_t fc_probe_bursts_chan (double* out, std::uint32_t cap)
{
    if (! haveBursts || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = bursts();
    const std::uint32_t room = cap / kBurstsChanStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels() && at < room; ++c, ++at)
    {
        out[at * kBurstsChanStride + 0] = (double) c;
        out[at * kBurstsChanStride + 1] = d.bandEnergy (c);
        out[at * kBurstsChanStride + 2] = (double) d.nonFiniteSamples (c);
        out[at * kBurstsChanStride + 3] = (double) d.absentSamples (c);
    }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_bursts_events (double* out, std::uint32_t cap)
{
    if (! haveBursts || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = bursts();
    const std::uint32_t room = cap / kBurstsEvtStride;
    std::uint32_t at = 0;
    for (std::int64_t i = 0; i < d.storedEventCount() && at < room; ++i, ++at)
    {
        const felitronics::analysis::BandBurst e = d.event (i);
        double* r = out + (std::size_t) at * kBurstsEvtStride;
        r[0] = (double) e.start;      r[1] = (double) e.length;  r[2] = (double) e.peakAt;
        r[3] = (double) e.hops;       r[4] = e.peakPower;        r[5] = e.peakBaseline;
        r[6] = e.peakExcessDb;        r[7] = e.peakWidePower;    r[8] = e.energy;
        r[9]  = e.touchedNonFinite ? 1.0 : 0.0;
        r[10] = e.baselineTouchedNonFinite ? 1.0 : 0.0;
        r[11] = e.closedByFinish ? 1.0 : 0.0;
    }
    return at;
}

// The two histograms are SPARSE on the CLI — it prints only non-zero bins — so the pairs (bin, mass) are
// what crosses, not a dense array. A dense one would make the JavaScript half decide which bins to print,
// i.e. re-implement a rule that lives in the C++ half.
FC_EXPORT std::uint32_t fc_probe_bursts_ioi (double* out, std::uint32_t cap)
{
    if (! haveBursts || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = bursts();
    const std::uint32_t room = cap / kBurstsBinStride;
    std::uint32_t at = 0;
    for (int b = 1; b <= felitronics::analysis::BandBursts::kIoiBins && at < room; ++b)
        if (d.intervalBin (b) != 0)
        { out[at * kBurstsBinStride + 0] = (double) b; out[at * kBurstsBinStride + 1] = (double) d.intervalBin (b); ++at; }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_bursts_lag (double* out, std::uint32_t cap)
{
    if (! haveBursts || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = bursts();
    const std::uint32_t room = cap / kBurstsBinStride;
    std::uint32_t at = 0;
    for (int b = 1; b <= felitronics::analysis::BandBursts::kMaxLag && at < room; ++b)
        if (d.lagBin (b) != 0)
        { out[at * kBurstsBinStride + 0] = (double) b; out[at * kBurstsBinStride + 1] = (double) d.lagBin (b); ++at; }
    return at;
}

//==============================================================================
// analysis::HumDetector (P74) through the ABI.
//
// This is the mode whose rows are not (field, channel) pairs but (channel, candidate, harmonic) triples,
// which is why the lists below carry their indices in the row rather than implying them from position:
// the harmonic table is channels x candidates x maxHarmonic and the CLI prints it nested, so a reader
// that inferred the coordinates from the row number would be re-deriving a loop rather than reading a
// measurement. `candidates` and `maxHarmonic` cross in the scalar block for the same reason — the JS half
// must not assume the geometry it is about to iterate.
namespace
{
    felitronics::analysis::HumDetector& hum()
    {
        static felitronics::analysis::HumDetector d;
        return d;
    }
    bool haveHum = false;

    constexpr std::uint32_t kHumScalars      = 8;
    constexpr std::uint32_t kHumChanStride   = 21;
    constexpr std::uint32_t kHumCandStride   = 19;
    constexpr std::uint32_t kHumHarmStride   = 8;
    constexpr std::uint32_t kHumStretchStride = 5;
}

FC_EXPORT int fc_probe_hum_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                double sampleRate)
{
    haveHum = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = hum();
    d.setParams (felitronics::analysis::HumDetectorParams {});
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! d.process (view, (int) channels, (int) frames)) return 0;
    d.finish();
    haveHum = true;
    return 1;
}

FC_EXPORT std::uint32_t fc_probe_hum_scalars_len   (void) { return kHumScalars; }
FC_EXPORT std::uint32_t fc_probe_hum_chan_stride   (void) { return kHumChanStride; }
FC_EXPORT std::uint32_t fc_probe_hum_cand_stride   (void) { return kHumCandStride; }
FC_EXPORT std::uint32_t fc_probe_hum_harm_stride   (void) { return kHumHarmStride; }
FC_EXPORT std::uint32_t fc_probe_hum_stretch_stride (void) { return kHumStretchStride; }

FC_EXPORT std::uint32_t fc_probe_hum_scalars (double* out, std::uint32_t cap)
{
    if (! haveHum || out == nullptr || cap < kHumScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = hum();
    const felitronics::analysis::HumDetectorParams hp {};
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();
    out[i++] = (double) d.channels();
    out[i++] = (double) d.geometry().order;
    out[i++] = (double) d.windowSamples();
    out[i++] = (double) d.hopSamples();
    out[i++] = d.binHz();
    out[i++] = (double) felitronics::analysis::HumDetector::kCandidates;
    out[i++] = (double) hp.maxHarmonic;
    return i;
}

FC_EXPORT std::uint32_t fc_probe_hum_chan (double* out, std::uint32_t cap)
{
    if (! haveHum || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = hum();
    const std::uint32_t room = cap / kHumChanStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels() && at < room; ++c, ++at)
    {
        const felitronics::analysis::HumReport r = d.report (c);
        double* w = out + (std::size_t) at * kHumChanStride;
        w[0] = (double) c;              w[1] = r.valid ? 1.0 : 0.0;   w[2] = (double) (int) r.reason;
        w[3] = (double) (int) r.mains;  w[4] = (double) r.baseHarmonic;
        w[5] = r.fundamentalObserved ? 1.0 : 0.0;  w[6] = r.fundamentalDerived ? 1.0 : 0.0;
        w[7] = r.fundamentalHz;         w[8]  = r.line.hz;            w[9]  = r.line.tonePower;
        w[10] = r.line.peakBinPower;    w[11] = r.line.floorPower;    w[12] = r.line.prominenceDb;
        w[13] = (double) r.frames;      w[14] = (double) r.finiteFrames;  w[15] = (double) r.holedFrames;
        w[16] = (double) r.quietFrames; w[17] = (double) r.quietStretches; w[18] = (double) r.storedStretches;
        w[19] = r.stretchesComplete ? 1.0 : 0.0;   w[20] = (double) r.tailUncoveredSamples;
    }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_hum_cand (double* out, std::uint32_t cap)
{
    if (! haveHum || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = hum();
    const std::uint32_t room = cap / kHumCandStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels(); ++c)
        for (int k = 0; k < felitronics::analysis::HumDetector::kCandidates && at < room; ++k, ++at)
        {
            const felitronics::analysis::HumCandidate q = d.candidate (c, k);
            double* w = out + (std::size_t) at * kHumCandStride;
            w[0] = (double) c;   w[1] = (double) k;   w[2] = q.nominalHz;
            w[3] = q.baseFound ? 1.0 : 0.0;           w[4] = (double) q.baseHarmonic;
            w[5] = q.fundamentalHz;
            w[6] = (double) q.stretchObservations;    w[7] = (double) q.stretchOffTolerance;
            w[8] = (double) q.frameObservations;
            w[9] = q.stretchSpreadHz;  w[10] = q.frameSpreadHz;  w[11] = q.maxIntraStretchSpreadHz;
            w[12] = q.stationary ? 1.0 : 0.0;  w[13] = q.passed ? 1.0 : 0.0;
            w[14] = (double) q.harmonicsObserved;     w[15] = (double) q.lowestHarmonicObserved;
            w[16] = q.windowPeak.found ? 1.0 : 0.0;   w[17] = q.windowPeak.hz;  w[18] = q.windowPeak.prominenceDb;
        }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_hum_harm (double* out, std::uint32_t cap)
{
    if (! haveHum || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = hum();
    const felitronics::analysis::HumDetectorParams hp {};
    const std::uint32_t room = cap / kHumHarmStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels(); ++c)
        for (int k = 0; k < felitronics::analysis::HumDetector::kCandidates; ++k)
            for (int h = 1; h <= hp.maxHarmonic && at < room; ++h, ++at)
            {
                const felitronics::analysis::HumHarmonic hh = d.harmonic (c, k, h);
                double* w = out + (std::size_t) at * kHumHarmStride;
                w[0] = (double) c;  w[1] = (double) k;  w[2] = (double) h;
                w[3] = hh.inBand ? 1.0 : 0.0;          w[4] = hh.peak.accepted ? 1.0 : 0.0;
                w[5] = hh.peak.hz;  w[6] = hh.peak.tonePower;  w[7] = hh.peak.prominenceDb;
            }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_hum_stretch (double* out, std::uint32_t cap)
{
    if (! haveHum || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = hum();
    const std::uint32_t room = cap / kHumStretchStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels(); ++c)
        for (std::int64_t i = 0; i < d.storedStretchCount (c) && at < room; ++i, ++at)
        {
            const felitronics::analysis::HumStretch st = d.stretch (c, i);
            double* w = out + (std::size_t) at * kHumStretchStride;
            w[0] = (double) c;  w[1] = (double) st.index;  w[2] = (double) st.startSample;
            w[3] = (double) st.endSample;  w[4] = (double) st.frames;
        }
    return at;
}

//==============================================================================
// analysis::SourceForensics (P73) through the ABI.
//
// The wall table has ONE MORE ROW THAN THERE ARE CHANNELS: index nc is the file's own aggregate wall, not
// a channel. That is the CLI's shape (`for c in 0..nc`) and it crosses unchanged, because collapsing it
// would make the JavaScript half decide what the aggregate is.
namespace
{
    felitronics::analysis::SourceForensics& forensics()
    {
        static felitronics::analysis::SourceForensics d;
        return d;
    }
    bool haveForensics = false;

    constexpr std::uint32_t kFxScalars    = 29;
    constexpr std::uint32_t kFxWallStride = 39;
    constexpr std::uint32_t kFxGridStride = 22;
}

FC_EXPORT int fc_probe_forensics_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                      double sampleRate)
{
    haveForensics = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = forensics();
    d.setParams (felitronics::analysis::SourceForensicsParams {});
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! d.process (view, (int) channels, (int) frames)) return 0;
    d.finish();
    haveForensics = true;
    return 1;
}

FC_EXPORT std::uint32_t fc_probe_forensics_scalars_len (void) { return kFxScalars; }
FC_EXPORT std::uint32_t fc_probe_forensics_wall_stride (void) { return kFxWallStride; }
FC_EXPORT std::uint32_t fc_probe_forensics_grid_stride (void) { return kFxGridStride; }
FC_EXPORT std::uint32_t fc_probe_forensics_khist_buckets (void)
{ return (std::uint32_t) felitronics::analysis::SourceForensics::gridExponentBuckets(); }

FC_EXPORT std::uint32_t fc_probe_forensics_scalars (double* out, std::uint32_t cap)
{
    if (! haveForensics || out == nullptr || cap < kFxScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = forensics();
    const felitronics::analysis::SourceForensicsParams fp {};
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();               out[i++] = (double) d.channels();
    out[i++] = (double) fp.fftOrder;         out[i++] = (double) d.hopSamples();
    out[i++] = (double) d.bins();            out[i++] = (double) d.binsPerCell();
    out[i++] = (double) d.exemptCells();     out[i++] = (double) d.distinctLimit();
    out[i++] = (double) d.plateauSpanCells();out[i++] = (double) d.floorSpanCells();
    out[i++] = d.cellHz();                   out[i++] = d.binHz();
    out[i++] = d.searchFromHz();             out[i++] = d.searchToHz();
    out[i++] = fp.cellWidthHz;               out[i++] = fp.searchFromHz;
    out[i++] = fp.plateauSpanHz;             out[i++] = fp.floorSpanHz;
    out[i++] = fp.transitionStartDb;         out[i++] = fp.transitionEndDb;
    out[i++] = fp.minDropDb;                 out[i++] = fp.maxTransitionHz;
    out[i++] = fp.nearNyquistFraction;       out[i++] = fp.emptyDb;
    out[i++] = fp.emptyMinHz;                out[i++] = fp.gridOutlierFraction;
    out[i++] = (double) d.samplesProcessed();out[i++] = (double) d.tailUncoveredSamples();
    out[i++] = (double) d.frames().frameCount();
    return i;
}

FC_EXPORT std::uint32_t fc_probe_forensics_wall (double* out, std::uint32_t cap)
{
    if (! haveForensics || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = forensics();
    const int nc = d.channels();
    const std::uint32_t room = cap / kFxWallStride;
    std::uint32_t at = 0;
    for (int c = 0; c <= nc && at < room; ++c, ++at)
    {
        const felitronics::analysis::SpectralWall w = (c < nc) ? d.wall (c) : d.wall();
        double* r = out + (std::size_t) at * kFxWallStride;
        r[0] = (double) c;
        r[1] = (double) (int) w.valid;                 r[2] = (double) (int) w.reason;
        r[3] = (double) (int) w.sharp;                 r[4] = (double) (int) w.nearNyquist;
        r[5] = (double) (int) w.transitionClipped;     r[6] = (double) (int) w.truncatedAtNyquist;
        r[7] = (double) w.exemptedCells;
        r[8] = (double) (int) w.secondValid;           r[9]  = (double) (int) w.secondSharp;
        r[10] = (double) (int) w.secondTransitionClipped; r[11] = (double) (int) w.secondTruncatedAtNyquist;
        r[12] = (double) (int) w.secondReason;
        r[13] = (double) (int) w.emptyAboveValid;      r[14] = (double) (int) w.emptyAboveReason;
        r[15] = (double) w.framesUsed;                 r[16] = (double) w.framesHoled;
        r[17] = w.cutoffHz;              r[18] = w.cutoffFractionOfNyquist; r[19] = w.steepestHz;
        r[20] = w.transitionEndHz;       r[21] = w.transitionHz;            r[22] = w.plateauPower;
        r[23] = w.floorLocalPower;       r[24] = w.maxAbovePower;           r[25] = w.sufMaxPower;
        r[26] = w.dropDb;                r[27] = w.strictDropDb;            r[28] = w.localDropDb;
        r[29] = w.recoveryDb;            r[30] = w.plateauSpreadDb;         r[31] = w.steepnessDbPerOctave;
        r[32] = w.secondCutoffHz;        r[33] = w.secondDropDb;            r[34] = w.secondTransitionHz;
        r[35] = w.emptyAboveHz;          r[36] = w.emptyAboveFractionOfNyquist;
        r[37] = w.emptyThresholdPower;   r[38] = w.peakCellPower;
    }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_forensics_grid (double* out, std::uint32_t cap)
{
    if (! haveForensics || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = forensics();
    const std::uint32_t room = cap / kFxGridStride;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels() && at < room; ++c, ++at)
    {
        const felitronics::analysis::SampleGrid g = d.sampleGrid (c);
        double* r = out + (std::size_t) at * kFxGridStride;
        r[0] = (double) c;
        r[1] = (double) (int) g.valid;            r[2] = (double) (int) g.reason;
        r[3] = (double) g.gridExponent;           r[4] = (double) (int) g.pcmCompatible;
        r[5] = (double) (int) g.outsidePcmRange;  r[6] = (double) g.minExactPcmBits;
        r[7] = (double) g.robustGridExponent;     r[8] = (double) g.robustPcmBits;
        r[9] = (double) g.alwaysZeroLowBits (24);
        r[10] = g.absPeak;    r[11] = g.sampleMin;   r[12] = g.sampleMax;
        r[13] = (double) g.nonZeroSamples;   r[14] = (double) g.zeroSamples;
        r[15] = (double) g.nonFiniteSamples; r[16] = (double) g.absentSamples;
        r[17] = (double) g.offGridSamples;   r[18] = (double) g.firstOffGridSample;
        r[19] = (double) g.firstMaxGridSample; r[20] = (double) g.distinctValues;
        r[21] = (double) (int) g.distinctComplete;
        }
    return at;
}

// One row per channel, `gridExponentBuckets()` wide, dense — the CLI prints every bucket.
FC_EXPORT std::uint32_t fc_probe_forensics_khist (double* out, std::uint32_t cap)
{
    if (! haveForensics || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = forensics();
    const std::uint32_t w = (std::uint32_t) felitronics::analysis::SourceForensics::gridExponentBuckets();
    const std::uint32_t room = cap / w;
    std::uint32_t at = 0;
    for (int c = 0; c < d.channels() && at < room; ++c, ++at)
    {
        const std::int64_t* h = d.gridExponentHistogram (c);
        for (std::uint32_t k = 0; k < w; ++k) out[(std::size_t) at * w + k] = (double) h[k];
    }
    return at;
}

//==============================================================================
// analysis::LowEnd (P75) through the ABI.
//
// One thing here is not a number: the note's NAME ("C#3"). It is a pure function of the peak MIDI value,
// so JavaScript could rebuild it from a pitch-class table — and that table would be a second copy of one
// that already exists in C++, which is the drift this whole ABI is shaped to avoid. So the name crosses
// from the module, formatted by the same function the CLI calls.
namespace
{
    felitronics::analysis::LowEnd& lowEnd()
    {
        static felitronics::analysis::LowEnd d;
        return d;
    }
    bool haveLowEnd = false;

    constexpr std::uint32_t kLeScalars     = 60;
    constexpr std::uint32_t kLeSeriesStride = 6;
    constexpr std::uint32_t kLeBandStride   = 11;
}

FC_EXPORT int fc_probe_lowend_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    haveLowEnd = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = lowEnd();
    d.setParams (felitronics::analysis::LowEndParams {});
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (! d.process (view, (int) channels, (int) frames)) return 0;
    if (! d.finish()) return 0;
    haveLowEnd = true;
    return 1;
}

FC_EXPORT std::uint32_t fc_probe_lowend_scalars_len  (void) { return kLeScalars; }
FC_EXPORT std::uint32_t fc_probe_lowend_series_stride (void) { return kLeSeriesStride; }
FC_EXPORT std::uint32_t fc_probe_lowend_band_stride   (void) { return kLeBandStride; }
FC_EXPORT std::uint32_t fc_probe_lowend_hist_bins     (void)
{ return (std::uint32_t) felitronics::analysis::LowEnd::kHistogramBins; }

FC_EXPORT std::uint32_t fc_probe_lowend_scalars (double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || cap < kLeScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = lowEnd();
    const felitronics::analysis::LowEndParams lp {};
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();            out[i++] = (double) d.channels();
    out[i++] = d.crossoverHz();           out[i++] = (double) lp.fftOrder;
    out[i++] = (double) d.hopSamples();   out[i++] = (double) d.blockSamples();
    out[i++] = (double) d.bandCount();    out[i++] = (double) fcore::Probe::kChunk;
    out[i++] = (double) (int) d.widthReason();  out[i++] = (double) (int) d.noteReason();
    out[i++] = (double) d.samplesProcessed();   out[i++] = (double) d.finiteSamples();
    out[i++] = (double) d.holeSamples();        out[i++] = (double) d.nonFiniteSamples();
    out[i++] = (double) d.absentSamples();      out[i++] = (double) d.filterNonFiniteSamples();
    out[i++] = d.lowMidEnergy();   out[i++] = d.lowSideEnergy();  out[i++] = d.highMidEnergy();
    out[i++] = d.highSideEnergy(); out[i++] = d.rawMidEnergy();   out[i++] = d.rawSideEnergy();
    out[i++] = d.lowSideFraction();
    out[i++] = d.highSideFraction();            out[i++] = d.rawSideFraction();
    out[i++] = (double) d.blockCount();         out[i++] = (double) d.storedBlockCount();
    out[i++] = d.blocksComplete() ? 1.0 : 0.0;  out[i++] = (double) d.histogramSamples();
    out[i++] = (double) d.worstFractionBlock(); out[i++] = d.worstFraction();  out[i++] = d.worstFractionEnergy();
    out[i++] = (double) d.peakEnergyBlock();    out[i++] = d.peakBlockEnergy(); out[i++] = d.peakEnergyBlockFraction();
    out[i++] = (double) d.peakSideEnergyBlock();out[i++] = d.peakBlockSideEnergy();
    out[i++] = d.peakLowSideAmplitude();        out[i++] = (double) d.peakLowSideAmplitudeAt();
    out[i++] = (double) d.usedFrames();         out[i++] = (double) d.holedFrames();
    out[i++] = (double) d.tailUncoveredSamples(); out[i++] = (double) d.windowSamples();
    out[i++] = (double) d.underResolvedBands();
    out[i++] = (double) d.peakBand();           out[i++] = (double) d.peakMidi();
    out[i++] = (double) d.peakDensityBand();    out[i++] = (double) d.secondBand();
    out[i++] = d.noteValid() ? 1.0 : 0.0;
    out[i++] = d.peakNoteHz();     out[i++] = d.peakCentroidHz();
    out[i++] = d.peakCentsOffset();out[i++] = d.peakBandSideFraction();
    out[i++] = d.frameEnergy();    out[i++] = d.bandRangeShare();
    out[i++] = d.backgroundDensity(); out[i++] = d.peakBandEnergy(); out[i++] = d.peakBandWidthHz();
    out[i++] = d.peakShare();      out[i++] = d.totalBandEnergy();
    return i;
}

// The note's name, formatted by the same functions the CLI calls. Returns the bytes written, or the bytes
// required when cap is 0. Empty when the note is not valid: the CLI prints a different line entirely then,
// and a name for a note that was not found would be a number that reads as a finding.
FC_EXPORT std::uint32_t fc_probe_lowend_note_name (char* out, std::uint32_t cap)
{
    if (! haveLowEnd) return 0u;
    const auto& d = lowEnd();
    if (! d.noteValid()) return 0u;
    const char* pc = felitronics::analysis::LowEnd::pitchClassName (d.peakMidi());
    const int oct = felitronics::analysis::LowEnd::noteOctave (d.peakMidi());
    char buf[16] {}; std::uint32_t n = 0;
    while (pc[n] != '\0' && n < 8) { buf[n] = pc[n]; ++n; }
    if (oct < 0) { buf[n++] = '-'; }
    const int a = oct < 0 ? -oct : oct;
    if (a >= 10) buf[n++] = (char) ('0' + (a / 10));
    buf[n++] = (char) ('0' + (a % 10));
    if (cap == 0u) return n;
    if (out == nullptr || ! outSpan (out, cap, 1) || cap < n) return 0u;
    for (std::uint32_t k = 0; k < n; ++k) out[k] = buf[k];
    return n;
}

FC_EXPORT std::uint32_t fc_probe_lowend_hist (double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = lowEnd();
    const std::uint32_t n = (std::uint32_t) felitronics::analysis::LowEnd::kHistogramBins;
    if (cap < n) return 0u;
    for (std::uint32_t i = 0; i < n; ++i) out[i] = (double) d.histogram ((int) i);
    return n;
}

FC_EXPORT std::uint32_t fc_probe_lowend_series (double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = lowEnd();
    const std::uint32_t room = cap / kLeSeriesStride;
    std::uint32_t at = 0;
    for (std::int64_t i = 0; i < d.storedBlockCount() && at < room; ++i, ++at)
    {
        const felitronics::analysis::LowEndBlock r = d.block (i);
        double* w = out + (std::size_t) at * kLeSeriesStride;
        w[0] = (double) r.index;  w[1] = (double) r.samples;  w[2] = (double) r.finiteSamples;
        w[3] = (double) r.holes;  w[4] = r.midEnergy;         w[5] = r.sideEnergy;
    }
    return at;
}

FC_EXPORT std::uint32_t fc_probe_lowend_bands (double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = lowEnd();
    const std::uint32_t room = cap / kLeBandStride;
    std::uint32_t at = 0;
    for (int b = 0; b < d.bandCount() && at < room; ++b, ++at)
    {
        const felitronics::analysis::LowEndBand r = d.band (b);
        double* w = out + (std::size_t) at * kLeBandStride;
        w[0] = (double) b;        w[1] = (double) r.midi;   w[2] = r.centreHz;
        w[3] = r.widthHz;         w[4] = r.binsPerBand;     w[5] = r.midEnergy;
        w[6] = r.sideEnergy;      w[7] = r.energy;          w[8] = r.density;
        w[9] = r.centroidHz;      w[10] = r.centsOffset;
    }
    return at;
}
