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
#include "fcore_stream.h"

#include <felitronics/analysis/BandBursts.h>
#include <felitronics/analysis/BandCrest.h>
#include <felitronics/analysis/HumDetector.h>
#include <felitronics/analysis/LowEnd.h>
#include <felitronics/analysis/PeakExcursions.h>
#include <felitronics/analysis/SourceForensics.h>
#include <felitronics/analysis/ProgrammeReport.h>

#include <cstdint>
#include <cstring>
#include <new>
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

    // THE ONE CHANNEL PREDICATE IN THIS FILE, read by the input checks below and by the five
    // fc_probe_<mode>_storage_bytes queries at the foot of it (P81). A second copy would be a seam the two
    // roads could drift apart along, which is the whole class this ABI keeps closing. The range is tested
    // BEFORE the narrowing to the analyzers' `int`: the conversion of a uint32 above INT_MAX is well defined
    // in C++20 and would land outside [1, kMaxChannels] anyway, but a width is refused here because it IS
    // out of range, not because a cast happened to carry it back out of range.
    bool geometry (std::uint32_t channels)
    {
        return channels >= 1 && channels <= (std::uint32_t) felitronics::core::kMaxChannels;
    }

    // The planar input span: non-null, non-empty, a width the core has, 4-byte aligned, and inside the heap.
    bool planarSpan (const float* planar, std::uint32_t frames, std::uint32_t channels)
    {
        if (planar == nullptr || frames == 0) return false;
        if (! geometry (channels)) return false;
        if ((reinterpret_cast<std::uintptr_t> (planar) & 0x3u) != 0) return false;   // a misaligned float* reads
                                                                                     // garbage in a release build
                                                                                     // and only traps under SAFE_HEAP
        // frames*channels must address real memory: on wasm32 the product is what a caller malloc'd, so a
        // wrapped one would hand us a window onto someone else's heap.
        const std::uint64_t bytes = (std::uint64_t) frames * (std::uint64_t) channels * sizeof (float);
        if (bytes > (std::uint64_t) 0xFFFFFFFFu) return false;
        return inHeap (planar, bytes);
    }

    // FOUR of the five offline analyzers (P77) accept an EMPTY programme and report on it — lowend does
    // NOT, because `fcore_measure lowend` refuses it too and the two roads must refuse the same set. — `fcore_measure report`
    // on /dev/null prints 2109 bytes of a perfectly good empty report. planarSpan() refuses frames == 0,
    // and rightly so for `clips`, whose emptiness is a FILE that was probably truncated; but here the
    // caller hands over a buffer, and "no audio" is a measurement, not a truncation. Refusing it in the
    // module while the CLI answers it is a byte-parity break, so these runs take this check instead.
    bool planarSpanOrEmpty (const float* planar, std::uint32_t frames, std::uint32_t channels)
    {
        if (! geometry (channels)) return false;
        if (frames == 0) return true;                                  // nothing to address, nothing to bound
        return planarSpan (planar, frames, channels);
    }

    // Ptr/size validation the core cannot do for us: everything below arrives from JS.
    bool viable (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
    {
        if (! planarSpan (planar, frames, channels)) return false;
        // rejects 0 / negative / NaN / +inf, every rate below the core's 8 kHz floor and the absurd-but-finite
        // ones above 768 kHz (see Probe::kMinSampleRate)
        return probe().prepare (sampleRate, (int) channels);
    }

    // An output span of `count` elements of `align` bytes each: non-null, aligned, inside the heap.
    bool outSpan (const void* out, std::uint32_t count, std::uint32_t align)
    {
        if (out == nullptr) return false;
        if ((reinterpret_cast<std::uintptr_t> (out) & (align - 1u)) != 0) return false;
        return inHeap (out, (std::uint64_t) count * align);
    }

    // `st.ok`, ALWAYS — never a bare `st.bytes()`. A REFUSED Storage IS NOT AN EMPTY ONE:
    // ProgrammeReport::Storage carries a DeterministicLoudnessMeter::Storage, whose bytes() has a constant
    // ring term (LoudnessMeter.h, `kSubRing`), so a default-constructed one reports 2400 bytes. A query that
    // forwarded that would quote a price for a measurement that cannot happen — and only for `report`, so a
    // test that probed the other four would not see it. (Measured on a9816e2; it is also why
    // ProgrammeReport.h's own "all zeros where prepare() refuses" note is true of the members and not of the
    // total.)
    template <typename Storage>
    double demand (const Storage& st) { return st.ok ? (double) st.bytes() : 0.0; }

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
// be served the previous file's numbers. Each of fc_probe_run, fc_probe_lufs and fc_probe_dbtp is a complete
// measurement from scratch: calling fc_probe_lufs() after fc_probe_run() re-runs the whole thing, and the
// getters then describe THAT run. (The fc_stream_* handles further down are the one road here that carries a
// measurement across calls, and they do it on purpose.)
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
    // The failure sentinel, spelled as the NUMBER it is rather than as a libm call that computes it.
    // `20.0 * std::log10 (1e-9)` is exactly -180.0 on Apple, glibc and musl, folded and at runtime
    // (measured, all six) — so this is bit-identical and buys two things. It stops a printed, diffed value
    // depending on a libm at all; and it closes the constant-fold hazard, where a compiler evaluates a
    // constant argument with its OWN high-precision arithmetic while the runtime call uses the row's libm,
    // and the same expression yields two different doubles in one binary. The neighbouring LUFS sentinel
    // at fc_probe_lufs() has always been a plain -120.0, which is the same decision made earlier.
    if (! run (planar, frames, channels, sampleRate)) return -180.0;
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
    // caps frames*channels*4 at 4 GiB, so a position in this ABI is below 2^30 and exact in a binary64.
    //
    // THE REASON USED TO BE WRITTEN HERE AS "an i64 would need -sWASM_BIGINT, which this module does not
    // build with". THAT IS FALSE on the toolchain this repo pins (emscripten 6.0.9, .github/workflows/ci.yml):
    // `WASM_BIGINT` defaults to TRUE in its src/settings.js, so an i64 DOES cross this boundary — measured,
    // as a JavaScript BigInt. The real reason for the double is what a BigInt then is on the page: a
    // different NUMERIC TYPE from every other value this ABI returns. `bigint + number` throws
    // TypeError ("Cannot mix BigInt and other types"), and JSON.stringify throws on it outright; only the
    // comparisons happen to work. A number a caller has to add to its own byte counts, or put in a report,
    // must not be one. (The same measurement corrected the sibling claim in tools/wasm/build.sh.)
    constexpr std::uint32_t kClipRunStride = 6;

    // Field k of a run in that order — ONE packing, read by fc_probe_clips_runs and by fc_stream_clips.
    double clipRunField (const felitronics::analysis::ClipRun& r, std::uint32_t k) noexcept
    {
        switch (k)
        {
            case 0:  return (double) r.start;
            case 1:  return (double) r.length;
            case 2:  return r.level;
            case 3:  return (double) r.channel;
            case 4:  return (double) r.sign;
            default: return (double) (int) r.evidence;
        }
    }
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
        return clipRunField (c.run ((std::int64_t) (i / kClipRunStride)), i % kClipRunStride);
    });
    return wrote == n ? m : 0u;                                   // a refused span writes nothing and says so
}

//==============================================================================
// THE STREAMING SURFACE (P120): loudness and clipped runs of a stream handed over piece by piece and read
// BETWEEN the pieces — fcore::StreamProbe (tools/fcore_stream.h), the class `fcore_measure stream` drives too.
// Read that header for why the meter is analysis::DeterministicLoudnessMeter and not the one fc_probe_run uses.
//
//   h = fc_stream_create (sampleRate, channels)          0 when refused
//   fc_stream_process (h, planar, frames)                1, or 0 — and a refused piece that carried samples
//                                                        poisons the stream: every reader below answers 0
//   fc_stream_loudness (h, out)                          fc_stream_loudness_fields() doubles: momentary,
//                                                        short-term, integrated (LUFS), samples consumed,
//                                                        gating blocks dropped past the one-hour store
//   fc_stream_clips_count (h)                            runs decided so far, stored or not
//   fc_stream_clips (h, from, out, cap)                  stored runs from index `from`, as fc_probe_clips_runs
//   fc_stream_finish (h)                                 decides the last runs; the stream takes no more audio
//   fc_stream_destroy (h)
//
// HANDLES, NOT A SINGLETON, because two streams must be measurable at once. A HANDLE IS A NUMBER, NEVER AN
// ADDRESS: a page hands it back, and an address handed back would be dereferenced whatever it was — stale after
// destroy, forged, or the 0 of a create that failed. So it is a serial looked up in a fixed table and never
// reused, and every entry point refuses one it does not know, which is planarSpan()'s "the call is refused, the
// module does not die" applied to the handle.
//
// AT MOST kMaxStreams AT ONCE. The page needs two. The bound is what a page can make the module hold: one stream
// asks for 3.1 MB at 48 kHz stereo and 29.5 MB at 768 kHz and 16 channels (the detector's deques, pending queue
// and 65536-run list, plus the meter's hour — measured natively), and this module is -fno-exceptions, where an
// allocation that fails aborts instead of refusing. Sixteen at the extreme stay far inside the 2 GiB heap.
namespace
{
    constexpr int kMaxStreams = 16;
    constexpr std::uint32_t kStreamLoudnessFields = 5;

    struct StreamSlot
    {
        fcore::StreamProbe* probe  = nullptr;
        std::uint32_t       handle = 0;
    };
    StreamSlot streamSlots[kMaxStreams] {};
    std::uint32_t lastStreamHandle = 0;

    fcore::StreamProbe* streamOf (std::uint32_t h) noexcept
    {
        if (h == 0) return nullptr;
        for (const auto& s : streamSlots) if (s.handle == h) return s.probe;
        return nullptr;
    }
}

FC_EXPORT std::uint32_t fc_stream_create (double sampleRate, std::uint32_t channels)
{
    if (! geometry (channels)) return 0u;
    StreamSlot* slot = nullptr;
    for (auto& s : streamSlots) if (s.probe == nullptr) { slot = &s; break; }
    if (slot == nullptr) return 0u;
    auto* p = new (std::nothrow) fcore::StreamProbe;
    if (p == nullptr) return 0u;
    if (! p->prepare (sampleRate, (int) channels)) { delete p; return 0u; }
    // Never 0 and never a live handle — a wrap past 2^32 creates would otherwise hand out a number in use.
    do { ++lastStreamHandle; } while (lastStreamHandle == 0u || streamOf (lastStreamHandle) != nullptr);
    slot->probe = p;
    slot->handle = lastStreamHandle;
    return slot->handle;
}

// `frames` samples of each of the stream's channels, planar as everywhere in this file. frames == 0 carries
// nothing and destroys nothing.
FC_EXPORT int fc_stream_process (std::uint32_t h, const float* planar, std::uint32_t frames)
{
    fcore::StreamProbe* s = streamOf (h);
    if (s == nullptr) return 0;
    if (frames == 0) return s->process (nullptr, 0) ? 1 : 0;
    const std::uint32_t channels = (std::uint32_t) s->channels();
    if (! planarSpan (planar, frames, channels)) { s->poison(); return 0; }
    const float* view[felitronics::core::kMaxChannels] {};
    for (std::uint32_t c = 0; c < channels; ++c) view[c] = planar + (std::size_t) c * (std::size_t) frames;
    // planarSpan() bounds frames*channels*4 by 4 GiB, so frames is below 2^30 and fits the instruments' int.
    return s->process (view, (int) frames) ? 1 : 0;
}

FC_EXPORT std::uint32_t fc_stream_loudness_fields (void) { return kStreamLoudnessFields; }

// Writes fc_stream_loudness_fields() doubles into `out` and returns 1, or writes nothing and returns 0.
FC_EXPORT int fc_stream_loudness (std::uint32_t h, double* out)
{
    const fcore::StreamProbe* s = streamOf (h);
    if (s == nullptr || ! s->valid() || ! outSpan (out, kStreamLoudnessFields, 8u)) return 0;
    const auto& m = s->meter();
    out[0] = m.momentaryLufs();
    out[1] = m.shortTermLufs();
    out[2] = m.integratedLufs();
    out[3] = (double) s->samples();
    out[4] = (double) m.droppedBlocks();
    return 1;
}

// Every run decided so far, including the ones past the list's capacity. A double, like every count here that
// a stream of unbounded length can grow.
FC_EXPORT double fc_stream_clips_count (std::uint32_t h)
{
    const fcore::StreamProbe* s = streamOf (h);
    return s != nullptr && s->valid() ? (double) s->detector().runCount() : 0.0;
}

FC_EXPORT std::uint32_t fc_stream_clips_stride (void) { return kClipRunStride; }

// The stored runs from index `from` on, as many whole ones as fit — `cap` IN DOUBLES and the return IN RUNS,
// the fc_probe_clips_runs rule and for its reason. A reader polls with `from` = the runs it already holds; 0
// means nothing new (or a refusal), never a short run.
FC_EXPORT std::uint32_t fc_stream_clips (std::uint32_t h, std::uint32_t from, double* out, std::uint32_t cap)
{
    const fcore::StreamProbe* s = streamOf (h);
    if (s == nullptr || ! s->valid()) return 0u;
    const auto& d = s->detector();
    const std::int64_t stored = d.storedRunCount();
    if ((std::int64_t) from >= stored) return 0u;
    const std::uint32_t left = (std::uint32_t) (stored - (std::int64_t) from);   // stored <= 65536
    const std::uint32_t fits = cap / kClipRunStride;
    const std::uint32_t m = fits < left ? fits : left;
    if (m == 0) return 0u;
    const std::uint32_t n = m * kClipRunStride;
    const std::uint32_t wrote = copyOut (out, n, n, [&] (std::uint32_t i)
    {
        return clipRunField (d.run ((std::int64_t) from + (std::int64_t) (i / kClipRunStride)), i % kClipRunStride);
    });
    return wrote == n ? m : 0u;
}

FC_EXPORT int fc_stream_finish (std::uint32_t h)
{
    fcore::StreamProbe* s = streamOf (h);
    return s != nullptr && s->finish() ? 1 : 0;
}

FC_EXPORT int fc_stream_destroy (std::uint32_t h)
{
    if (h == 0) return 0;
    for (auto& s : streamSlots)
        if (s.handle == h)
        {
            delete s.probe;
            s = StreamSlot {};
            return 1;
        }
    return 0;
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

    // ONE definition of the parameters this ABI measures with, read by BOTH roads — the run below and
    // fc_probe_report_storage_bytes at the foot of this file. Two `Params {}` literals would be two
    // definitions of the same thing, and the published price would go on matching the run only by
    // coincidence from the moment either one grew an argument. Same constant in all five modes.
    constexpr felitronics::analysis::ProgrammeReportParams kReportParams {};

    constexpr std::uint32_t kReportCountStride = 2;   // channel, value
    constexpr std::uint32_t kReportValueStride = 4;   // channel, valid, reason, value
}

FC_EXPORT int fc_probe_report_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    haveReport = false;
    if (! planarSpanOrEmpty (planar, frames, channels)) return 0;
    auto& p = programme();
    // Installed rather than left to the instance's own defaults, which are the same values: it makes this
    // road READ the constant the price query reads, so the two cannot drift apart silently. The other four
    // modes already called setParams() for their own reasons.
    p.setParams (kReportParams);
    if (! p.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    // Guarded on `frames`, not merely skipped later: `planar + k * frames` is undefined behaviour when
    // planar is null EVEN IF the offset is zero, and an empty programme is allowed to arrive with a null
    // pointer. JavaScript's _malloc(0) happens to return something non-null, so the manifestation is
    // theoretical from that road — and a direct C ABI caller is not obliged to be as lucky.
    if (frames != 0)
        for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (frames != 0 && ! p.process (view, (int) channels, (int) frames)) return 0;
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
// `enterDb` and `exitDb` are read from the parameters this ABI measures with — kBurstsParams, since P81 the
// same constant `fc_probe_bursts_run` installs and `fc_probe_bursts_storage_bytes` prices — and not from
// the finished measurement, because the detector does not report the thresholds it used back out of the
// report (it does publish its PENDING params through `params()`, which is a different thing: what the next
// prepare() would take, not what this one did). The CLI reads the same defaults the same way. That is the
// one number on this road not read back out of the measurement; it is the documented default on both
// sides, so the diff still covers it, but it is named rather than hidden.
namespace
{
    felitronics::analysis::BandBursts& bursts()
    {
        static felitronics::analysis::BandBursts d;
        return d;
    }
    bool haveBursts = false;

    // One definition, two roads — see kReportParams above.
    constexpr felitronics::analysis::BandBurstsParams kBurstsParams {};

    // WHAT THE LAST RUN ACTUALLY INSTALLED. Until K3 the scalars read `kBurstsParams` for `enterDb` and
    // `exitDb`, which was true while the only road installed that constant and became a LIE the moment a
    // second road could install something else: the two thresholds would have described the defaults while
    // the events came from the caller's. The detector does not report its own thresholds back out of a
    // finished measurement, so the honest fix is to remember what was installed and read THAT — one
    // variable written by every road into the detector, rather than a constant beside them.
    felitronics::analysis::BandBurstsParams installedBursts = kBurstsParams;

    constexpr std::uint32_t kBurstsScalars    = 32;
    constexpr std::uint32_t kBurstsChanStride = 4;
    constexpr std::uint32_t kBurstsEvtStride  = 12;
    constexpr std::uint32_t kBurstsBinStride  = 2;
}

namespace
{
    // THE ONE ROAD INTO THE DETECTOR. Both entry points below are this function with a different parameter
    // set, so the ordering that matters — disarm, validate, install, prepare, process, finish — has one
    // spelling and a caller cannot reach a state one road guards and the other does not.
    int runBursts (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate,
                   const felitronics::analysis::BandBurstsParams& p)
    {
        haveBursts = false;
        if (! planarSpanOrEmpty (planar, frames, channels)) return 0;
        auto& d = bursts();
        d.setParams (p);
        if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
        // INSTALLED ONLY ONCE THE PREPARE HELD, so that what is remembered is what actually ran. A refused
        // run is silent anyway — `haveBursts` is cleared on the first line, the established contract for
        // every road here — so no reader can see these; recording them before the prepare would simply
        // leave a false answer waiting for the NEXT successful run to be read beside.
        installedBursts = p;
        const float* view[felitronics::core::kMaxChannels] {};
        if (frames != 0)
            for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
        if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
        d.finish();
        haveBursts = true;
        return 1;
    }
}

FC_EXPORT int fc_probe_bursts_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    return runBursts (planar, frames, channels, sampleRate, kBurstsParams);
}

// K3 — THE SAME DETECTOR, THE CALLER'S BAND. The default 5-9 kHz band answers "sibilance"; the same
// machinery at 80 Hz - 8 kHz with a shorter baseline answers "how dense are the transients", which is a
// different question and not a different analyzer. `maxEvents` is NOT exposed: it is the event list's
// capacity, the counters keep counting past it, and `fc_probe_bursts_events_complete` already says whether
// the list is whole — a caller that could shrink it could only make the list lie about itself.
//
// EVERY VALUE IS THE CORE'S TO JUDGE. `prepare` refuses a band it cannot build (0.49*fs must clear the top
// corner), a non-positive hop, a baseline shorter than a hop; this returns 0 for all of it, exactly as the
// default road does, and a refused run leaves the getters SILENT — the contract every road here already
// had, and the one state a reader can act on without knowing which road refused.
FC_EXPORT int fc_probe_bursts_run_with (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                        double sampleRate, double bandLowHz, double bandHighHz,
                                        double hopMs, double baselineMs, double enterDb, double exitDb)
{
    felitronics::analysis::BandBurstsParams p;      // maxEvents keeps its default: see above
    p.bandLowHz = bandLowHz; p.bandHighHz = bandHighHz;
    p.hopMs = hopMs; p.baselineMs = baselineMs;
    p.enterDb = enterDb; p.exitDb = exitDb;
    return runBursts (planar, frames, channels, sampleRate, p);
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
    const auto& bp = installedBursts;    // WHAT THIS MEASUREMENT RAN AT — see the note on the variable
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
    out[i++] = d.onsetsPerSecond();      // K3 — over the JUDGED programme; the analyzer owns the denominator
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

    // One definition, two roads — see kReportParams above.
    constexpr felitronics::analysis::HumDetectorParams kHumParams {};

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
    if (! planarSpanOrEmpty (planar, frames, channels)) return 0;
    auto& d = hum();
    d.setParams (kHumParams);
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    // Guarded on `frames`, not merely skipped later: `planar + k * frames` is undefined behaviour when
    // planar is null EVEN IF the offset is zero, and an empty programme is allowed to arrive with a null
    // pointer. JavaScript's _malloc(0) happens to return something non-null, so the manifestation is
    // theoretical from that road — and a direct C ABI caller is not obliged to be as lucky.
    if (frames != 0)
        for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
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
    const auto& hp = kHumParams;         // the constant the run and the price read — not a fresh default
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
    const auto& hp = kHumParams;         // the constant the run and the price read — not a fresh default
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

    // One definition, two roads — see kReportParams above.
    constexpr felitronics::analysis::SourceForensicsParams kFxParams {};

    constexpr std::uint32_t kFxScalars    = 29;
    constexpr std::uint32_t kFxWallStride = 39;
    constexpr std::uint32_t kFxGridStride = 22;
}

FC_EXPORT int fc_probe_forensics_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                      double sampleRate)
{
    haveForensics = false;
    if (! planarSpanOrEmpty (planar, frames, channels)) return 0;
    auto& d = forensics();
    d.setParams (kFxParams);
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    // Guarded on `frames`, not merely skipped later: `planar + k * frames` is undefined behaviour when
    // planar is null EVEN IF the offset is zero, and an empty programme is allowed to arrive with a null
    // pointer. JavaScript's _malloc(0) happens to return something non-null, so the manifestation is
    // theoretical from that road — and a direct C ABI caller is not obliged to be as lucky.
    if (frames != 0)
        for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
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
    const auto& fp = kFxParams;          // the constant the run and the price read — not a fresh default
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

    // One definition, two roads — see kReportParams above.
    constexpr felitronics::analysis::LowEndParams kLeParams {};

    constexpr std::uint32_t kLeScalars     = 68;
    constexpr std::uint32_t kLeSeriesStride = 6;
    constexpr std::uint32_t kLeBandStride   = 13;
    constexpr std::uint32_t kLeLowestFields = 7;

    // WHAT THE LAST SUCCESSFUL RUN INSTALLED, not the compile-time default. The scalars publish the
    // geometry a caller must divide by, and the moment `_run_with` exists a constant there is a lie:
    // fftOrder, the note range and the duty threshold would keep reporting the defaults while the
    // measurement used something else. Same arrangement as the crest section above.
    felitronics::analysis::LowEndParams installedLe = kLeParams;
}

static int lowEndRunWith (const float* planar, std::uint32_t frames, std::uint32_t channels,
                          double sampleRate, const felitronics::analysis::LowEndParams& lp)
{
    haveLowEnd = false;
    // planarSpan, NOT planarSpanOrEmpty: `fcore_measure lowend` REFUSES an empty programme where the
    // other four report on one, and the two roads have to refuse the same set or the diff is comparing
    // different contracts. (That native refusal is LowEnd's own: an empty programme has no low end to
    // measure, and it says so by exiting rather than by publishing an invalid report.)
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = lowEnd();
    d.setParams (lp);
    if (! d.prepare (sampleRate, (int) fcore::Probe::kChunk, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    // Guarded on `frames`, not merely skipped later: `planar + k * frames` is undefined behaviour when
    // planar is null EVEN IF the offset is zero, and an empty programme is allowed to arrive with a null
    // pointer. JavaScript's _malloc(0) happens to return something non-null, so the manifestation is
    // theoretical from that road — and a direct C ABI caller is not obliged to be as lucky.
    if (frames != 0)
        for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
    if (! d.finish()) return 0;
    installedLe = lp;                    // only on success: a refused run leaves the last good geometry
    haveLowEnd = true;
    return 1;
}

FC_EXPORT int fc_probe_lowend_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                   double sampleRate)
{
    return lowEndRunWith (planar, frames, channels, sampleRate, kLeParams);
}

// The same measurement with the geometry a caller chose. Every argument is VALIDATED BY THE CORE, not
// here: `storageFor` refuses exactly what `prepare` refuses, and a second opinion in this file would be a
// second definition waiting to drift. A refused run answers 0 and changes nothing a reader can see.
FC_EXPORT int fc_probe_lowend_run_with (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                        double sampleRate, double crossoverHz, double lowNoteHz,
                                        double highNoteHz, std::int32_t fftOrder, double dutyThresholdDb,
                                        std::int32_t skipBlocks)
{
    felitronics::analysis::LowEndParams lp = kLeParams;
    lp.crossoverHz     = crossoverHz;
    lp.lowNoteHz       = lowNoteHz;
    lp.highNoteHz      = highNoteHz;
    lp.fftOrder        = (int) fftOrder;
    lp.dutyThresholdDb = dutyThresholdDb;
    lp.skipBlocks      = (int) skipBlocks;
    return lowEndRunWith (planar, frames, channels, sampleRate, lp);
}

FC_EXPORT std::uint32_t fc_probe_lowend_scalars_len  (void) { return kLeScalars; }
FC_EXPORT std::uint32_t fc_probe_lowend_lowest_fields (void) { return kLeLowestFields; }
FC_EXPORT std::uint32_t fc_probe_lowend_series_stride (void) { return kLeSeriesStride; }
FC_EXPORT std::uint32_t fc_probe_lowend_band_stride   (void) { return kLeBandStride; }
FC_EXPORT std::uint32_t fc_probe_lowend_hist_bins     (void)
{ return (std::uint32_t) felitronics::analysis::LowEnd::kHistogramBins; }

FC_EXPORT std::uint32_t fc_probe_lowend_scalars (double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || cap < kLeScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = lowEnd();
    const auto& lp = installedLe;        // what the RUN installed — see the note beside it
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
    // K9. The resolution boundary, so `underResolvedBands` above is actionable rather than a count; the
    // duty population, so a consumer never has to define it; the installed threshold, because a run may
    // not have used the default; and what the histogram did not count.
    out[i++] = (double) d.firstResolvedBand();  out[i++] = d.resolvedAboveHz();
    out[i++] = (double) felitronics::analysis::LowEnd::lobeBins();
    out[i++] = (double) d.dutyFrames();         out[i++] = d.dutyThresholdDb();
    out[i++] = (double) d.skippedBlocks();      out[i++] = (double) lp.skipBlocks;
    // The low band's share of the whole programme, named rather than left to be recomposed — see the
    // core's note: it is LR4-weighted, and content at the crossover counts at a quarter.
    out[i++] = d.infraLowShare();
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
        // K9: how OFTEN the band was there, and how loud it was when it was. The duty itself is
        // count/dutyFrames — one division both roads do identically, from two numbers both published.
        w[11] = (double) d.dutyCount (b);       w[12] = d.levelWhenOnDb (b);
    }
    return at;
}

//==============================================================================
// THE PRICE OF A MEASUREMENT, ASKED BEFORE IT IS PAID (P81).
//
// Every one of the five analyzers above publishes what it is about to need through a public
// `storageFor()` and only then allocates — law 11d, whose whole point is that the budget and the
// allocation are ONE function and therefore cannot drift apart. None of the five `fc_probe_<mode>_run`
// entry points passed that number out: JavaScript called `_run`, which went straight into an allocating
// `prepare()`, and there was nothing to ask the price with.
//
// WHY THAT MATTERS HERE AND NOT NATIVELY. This module is built `-fno-exceptions` (tools/wasm/build.sh),
// and there a failed allocation has no way to say no: `std::vector::assign` calls `operator new`, the
// throw is compiled out, and libc++ calls abort(). The C entry point then NEVER RETURNS — a caller
// reading `_fc_probe_hum_run(...) === 1` is not handed a 0, it is handed a thrown exception from a
// function it was reading a status out of.
//
// MEASURED, and measured more narrowly than the first telling of it. With emscripten's linear memory
// squeezed to 1 932 984 320 bytes of its 2 147 483 648 ceiling and a 64-byte input,
// `fc_probe_hum_run (ptr, 1, 16, 768000)` produced, on the checked module, exactly this:
// "Cannot enlarge memory, requested 2 252 435 248 bytes, but the limit is 2 147 483 648" ->
// "bad_alloc was thrown in -fno-exceptions mode" -> "Aborted(native code called abort())". What is NOT
// true, and was assumed before it was checked: the page does not die and the module is not poisoned. The
// WebAssembly.RuntimeError IS catchable in JavaScript, and after catching it the module went on
// answering — the getters of the aborted mode read 0, as the `have*` discipline promises, `_free` worked,
// and a later affordable run at 48 kHz stereo returned 1. That recovery is an observation and not a
// guarantee: abort() does not unwind, so whatever the allocation was half-way through is left where it
// stood. The defect this section closes is therefore the narrow, certain one — A MEASUREMENT THAT CANNOT
// BE PRICED AND CANNOT REFUSE — and not a crash story.
//
// And the prices are not small: `HumDetector::storageFor()` at 768 kHz and 16 channels asks for
// 352 688 184 bytes, the five together for 377 423 824, against a 2 GiB ceiling that also has to hold the
// caller's decoded input.
//
// WHAT THESE FIVE ANSWER, EXACTLY
//
//  * A DEMAND, NOT A FOOTPRINT, AND THE PAYLOAD ONLY. It is the sum of the byte counts the NEXT prepare()
//    at this geometry will pass to `operator new`, on a fresh object. It excludes, deliberately and
//    without apology: allocator metadata and alignment padding, fragmentation, the growth granularity of
//    the wasm heap, THE CALLER'S OWN INPUT BUFFER (frames * channels * 4, which only the caller knows),
//    and whatever the other four analyzers are still holding. It is also not what THIS analyzer is
//    holding: `std::vector::assign` never gives capacity back, so an instance that has already run wider
//    owns more than this number says (BandBursts::Storage spells the same caveat for itself), and one
//    that has already run at THIS geometry asks for nothing at all — except SourceForensics, whose
//    second prepare() at the same geometry asks for 688 144 of its 2 080 152 again. Measured, all of it.
//  * A POSITIVE ANSWER IS NOT A PROMISE. This entry point cannot reserve anything, and it is a number to
//    make a policy with, not an allocation that succeeded. That is also why the number is PUBLISHED here
//    rather than turned into a ceiling: a browser's memory budget is a product decision belonging to the
//    page, not a constant this file gets to invent. What the page is NOT given here is the scale to judge
//    the number against: emscripten's own glue keeps the heap maximum to itself, and while it can be read
//    out of the wasm bytes — tools/wasm/check-no-threads.mjs already parses that memory limit — that is a
//    build-time reading, not something a worker does at run time. If it is ever wanted it is one more fact
//    of module identity beside fc_probe_sizeof_longdouble, and still not a ceiling.
//  * ZERO MEANS THE GEOMETRY IS REFUSED, and it cannot be read as a free measurement — but the reason is
//    STRUCTURAL and not a measured minimum. An accepted geometry always buys something: four of the five
//    size at least one record per INPUT channel, and LowEnd — which does not, its spectrum being sized for
//    two fixed mid/side axes — still buys a band table and a block store. None of the five can be accepted
//    and cost nothing. A grid is only
//    a sample and says nothing about the rest of the domain: the smallest demand a fine sweep of the
//    accepted domain FOUND is 416 792 bytes (report, 8050 Hz, one channel; 158 240 — hum at 1517 Hz — before
//    P51 moved the rate floor to 8 kHz, and the same sweep reproduces that number on the old tree), which is
//    an observation. What
//    the gate pins is the equivalence, not a constant — felitronics_analysis_abi_tests walks a grid of
//    widths and rates, including each mode's own admission floor, and asserts that the query is positive
//    exactly where `_run` is accepted. The zero is a canonical +0.0, which is what lets a caller write
//    `need > 0` and a test write `! signbit`.
//  * NO AUDIO POINTER AND NO FRAME COUNT, on purpose, and NOT the argument list of `_run`. The price has
//    to be askable BEFORE the input buffer exists — that is the whole use — so these entry points neither
//    take nor dereference one, and `_run` keeps its own buffer checks (planarSpan / planarSpanOrEmpty,
//    and lowend's refusal of an empty programme, which is an input contract and not a geometry). What the
//    two roads must agree about is the GEOMETRY, and that is what the gate compares. Nothing is lost by
//    dropping the length: no preparation here is sized by the programme. The one analyzer sized by a
//    duration, ProgrammeReport, is sized by `ProgrammeReportParams::maxDurationSec` — a fixed hour — and a
//    longer programme overflows bounded stores and says so, it does not grow them.
//  * THE SAME PARAMETERS AND THE SAME BLOCK AS THE RUN. Each mode's defaults are one `constexpr` constant
//    read by both roads (kReportParams and its four siblings). `report` is the only one whose storageFor()
//    takes a maxBlock, and it is given `fcore::Probe::kChunk` — the value fc_probe_report_run() hands to
//    prepare(), not the caller's frame count.
//
// `double` AND NOT `std::uint64_t`. Not for the reason this file used to give — see the corrected note
// above fc_probe_clips_runs: i64 does cross this boundary on the pinned toolchain, as a BigInt, and a
// BigInt is the wrong shape for a byte count a page has to add to its own input size and put in a report:
// `bigint + number` throws TypeError and JSON.stringify refuses it, so every use site would have to
// convert it first — which is a conversion, and a place for a conversion to be wrong. Every demand in the
// accepted domain is far below 2^53, so a binary64
// carries it exactly; a `std::uint32_t` would arrive in JavaScript through a signed i32 and read negative
// for anything past 2^31, which the 2 GiB ceiling makes reachable in principle.
//
// THE NAME. `_storage_bytes` and not `_need`, which is what the sibling ABI calls the same idea
// (fc_master_need, DSP-ARCHITECTURE.md §law 11d): this side of the tree speaks `Storage` / `storageFor`,
// and the export is named after the function whose number it is carrying.

FC_EXPORT double fc_probe_report_storage_bytes (std::uint32_t channels, double sampleRate)
{
    if (! geometry (channels)) return 0.0;
    return demand (felitronics::analysis::ProgrammeReport::storageFor (
                       sampleRate, (int) fcore::Probe::kChunk, (int) channels, kReportParams));
}

//==================================================================================================
// K1 — BandCrest: the crest of a programme per 400 ms block and per band, and the PAIRED loss between a source
// and its master.
//
// TWO SLOTS, NOT ONE, and that is the whole reason this surface looks different from its neighbours. Every
// other analyzer here holds one static result because every other question is about one programme. This one
// is about the DIFFERENCE between two, and the loss is arithmetic that must not be reimplemented on the other
// side of the ABI: formed from linear cells it is one logarithm and invariant to the gain a master has over
// its source, and formed from two dB numbers in JavaScript it would be neither. So both runs live here and
// `fc_probe_crest_loss` does the comparing.
//
// SLOT 0 IS THE SOURCE AND SLOT 1 IS THE MASTER, and the order is load-bearing rather than a convention: the
// activity mask comes from the SOURCE alone, because a mask taken from the output can drop exactly the blocks
// the chain damaged most.
namespace
{
    felitronics::analysis::BandCrest& crest (int slot)
    {
        static felitronics::analysis::BandCrest a, b;
        return slot == 0 ? a : b;
    }
    bool haveCrest[2] { false, false };
    felitronics::analysis::BandCrestParams installedCrest[2] {};

    constexpr std::uint32_t kCrestScalars     = 22;   // 16 + one active-block count per band + the level
    constexpr std::uint32_t kCrestBlockStride = 15;   // 5 bands x { peak, meanSq, active }
    constexpr std::uint32_t kCrestLossFields  = 17;   // + outSilent and the coarse block lag

    int runCrest (int slot, const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate,
                  const felitronics::analysis::BandCrestParams& p)
    {
        if (slot < 0 || slot > 1) return 0;
        haveCrest[slot] = false;
        if (! planarSpanOrEmpty (planar, frames, channels)) return 0;
        auto& d = crest (slot);
        d.setParams (p);
        if (! d.prepare (sampleRate, (int) channels, (long long) frames)) return 0;
        installedCrest[slot] = p;
        const float* view[felitronics::core::kMaxChannels] {};
        if (frames != 0)
            for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
        if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
        d.finish();
        haveCrest[slot] = true;
        return 1;
    }
}

FC_EXPORT int fc_probe_crest_run (std::int32_t slot, const float* planar, std::uint32_t frames,
                                  std::uint32_t channels, double sampleRate)
{
    return runCrest ((int) slot, planar, frames, channels, sampleRate, felitronics::analysis::BandCrestParams {});
}

FC_EXPORT int fc_probe_crest_run_with (std::int32_t slot, const float* planar, std::uint32_t frames,
                                       std::uint32_t channels, double sampleRate,
                                       double edge0Hz, double edge1Hz, double edge2Hz,
                                       double hopMs, std::int32_t blockHops,
                                       double programmeFloorDb, double bandShareFloorDb)
{
    felitronics::analysis::BandCrestParams p;
    p.bandEdgeHz[0] = edge0Hz; p.bandEdgeHz[1] = edge1Hz; p.bandEdgeHz[2] = edge2Hz;
    p.hopMs = hopMs; p.blockHops = (int) blockHops;
    p.programmeFloorDb = programmeFloorDb; p.bandShareFloorDb = bandShareFloorDb;
    return runCrest ((int) slot, planar, frames, channels, sampleRate, p);
}

FC_EXPORT std::uint32_t fc_probe_crest_scalars_len  (void) { return kCrestScalars; }
FC_EXPORT std::uint32_t fc_probe_crest_block_stride (void) { return kCrestBlockStride; }
FC_EXPORT std::uint32_t fc_probe_crest_loss_len     (void) { return kCrestLossFields; }

// The scalars, IN THE ORDER THE CLI PRINTS THEM. The order is the contract.
FC_EXPORT std::uint32_t fc_probe_crest_scalars (std::int32_t slot, double* out, std::uint32_t cap)
{
    if (slot < 0 || slot > 1 || ! haveCrest[slot]) return 0u;
    if (out == nullptr || cap < kCrestScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = crest ((int) slot);
    const auto& p = installedCrest[slot];
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();
    out[i++] = (double) d.channels();
    out[i++] = (double) d.hopSamples();
    out[i++] = (double) d.blockHops();
    out[i++] = p.bandEdgeHz[0];
    out[i++] = p.bandEdgeHz[1];
    out[i++] = p.bandEdgeHz[2];
    out[i++] = p.programmeFloorDb;
    out[i++] = p.bandShareFloorDb;
    out[i++] = (double) d.samplesProcessed();
    out[i++] = (double) d.hopCount();
    // THE BASE-RATE CLOCK'S OWN COUNT, published rather than derived: it must equal the one above, and the two
    // are counted by different code over the same time. They disagreed once, by a whole hop.
    out[i++] = (double) d.basePeakHops();
    out[i++] = (double) d.blockCount();
    out[i++] = (double) d.nonFiniteSamples();
    out[i++] = (double) (int) d.invalidReason();
    out[i++] = d.fullBandPeakLin();     // the number the certifying true-peak meter also answers
    // THE POPULATION, per band, published rather than left to be counted from the block table. A consumer
    // deciding on a statistic needs to know how many blocks it rests on without walking the rows, and a count
    // it derives itself is a second definition of "active" waiting to drift from this one.
    for (int b = 0; b < felitronics::analysis::BandCrest::kBands; ++b)
        out[i++] = (double) d.activeBlocks (b);
    // The programme's own level in the gate's units, so a caller can set `programmeFloorDb` RELATIVE to it
    // without converting between two quantities that are not the same one. See the header.
    out[i++] = d.programmeMeanSquareDb();
    return i;
}

// One row per block: for each of the five bands, the linear peak, the mean square and the ACTIVE flag. Linear,
// not dB — the loss is a ratio of these and a dB store would put a logarithm in the middle of it.
FC_EXPORT std::uint32_t fc_probe_crest_blocks (std::int32_t slot, double* out, std::uint32_t cap)
{
    if (slot < 0 || slot > 1 || ! haveCrest[slot]) return 0u;
    if (out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = crest ((int) slot);
    const auto rows = (std::uint64_t) d.blockCount();
    if (rows * kCrestBlockStride > (std::uint64_t) cap) return 0u;        // all of it or none of it
    std::uint32_t i = 0;
    for (long long j = 0; j < d.blockCount(); ++j)
        for (int b = 0; b < felitronics::analysis::BandCrest::kBands; ++b)
        {
            out[i++] = d.blockPeakLin (j, b);
            out[i++] = d.blockMeanSq  (j, b);
            out[i++] = d.blockActive  (j, b) ? 1.0 : 0.0;
        }
    return (std::uint32_t) rows;
}

// THE COMPARISON, for one band: slot 0 as the source, slot 1 as the master. The mask is the source's.
FC_EXPORT std::uint32_t fc_probe_crest_loss (std::int32_t band, double* out, std::uint32_t cap)
{
    if (! haveCrest[0] || ! haveCrest[1]) return 0u;
    if (band < 0 || band >= felitronics::analysis::BandCrest::kBands) return 0u;
    if (out == nullptr || cap < kCrestLossFields || ! outSpan (out, cap, 8)) return 0u;
    static std::vector<double> scratch;
    const auto L = felitronics::analysis::bandCrestLoss (crest (0), crest (1), (int) band, scratch);
    std::uint32_t i = 0;
    out[i++] = (double) L.blocks;
    out[i++] = (double) L.inActive;
    out[i++] = (double) L.usable;
    out[i++] = L.p50Db;
    out[i++] = L.p95Db;
    out[i++] = L.cvar95Db;
    out[i++] = L.meanDb;
    out[i++] = L.maxDb;
    out[i++] = L.p5Db;
    out[i++] = (double) L.over1Db;
    out[i++] = (double) L.over3Db;
    out[i++] = (double) L.over6Db;
    out[i++] = L.peakShiftDb;
    out[i++] = L.levelShiftDb;
    out[i++] = (double) L.outSilent;
    out[i++] = (double) L.lagBlocks;
    out[i++] = L.valid ? 1.0 : 0.0;
    return i;
}

FC_EXPORT double fc_probe_crest_storage_bytes (std::uint32_t channels, double sampleRate, std::uint32_t frames)
{
    if (! geometry (channels)) return 0.0;
    const auto st = felitronics::analysis::BandCrest::storageFor (sampleRate, (int) channels, (long long) frames,
                                                                  felitronics::analysis::BandCrestParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}

FC_EXPORT double fc_probe_crest_storage_bytes_with (std::uint32_t channels, double sampleRate, std::uint32_t frames,
                                                    double edge0Hz, double edge1Hz, double edge2Hz,
                                                    double hopMs, std::int32_t blockHops,
                                                    double programmeFloorDb, double bandShareFloorDb)
{
    if (! geometry (channels)) return 0.0;
    felitronics::analysis::BandCrestParams p;
    p.bandEdgeHz[0] = edge0Hz; p.bandEdgeHz[1] = edge1Hz; p.bandEdgeHz[2] = edge2Hz;
    p.hopMs = hopMs; p.blockHops = (int) blockHops;
    p.programmeFloorDb = programmeFloorDb; p.bandShareFloorDb = bandShareFloorDb;
    const auto st = felitronics::analysis::BandCrest::storageFor (sampleRate, (int) channels, (long long) frames, p);
    return st.ok ? (double) st.bytes() : 0.0;
}

FC_EXPORT double fc_probe_bursts_storage_bytes (std::uint32_t channels, double sampleRate)
{
    if (! geometry (channels)) return 0.0;
    return demand (felitronics::analysis::BandBursts::storageFor (
                       sampleRate, (int) channels, kBurstsParams));
}

// K3's price, and it exists because the default one CANNOT answer for a parameterised run: the baseline
// ring is `round(baselineMs / hopMs)` hops of `hopMs` each, so both of those move the allocation and a page
// that sized itself by the default figure would be short exactly where it asked for a longer memory. It
// takes the same six values as `fc_probe_bursts_run_with` so the pair cannot drift apart.
FC_EXPORT double fc_probe_bursts_storage_bytes_with (std::uint32_t channels, double sampleRate,
                                                     double bandLowHz, double bandHighHz,
                                                     double hopMs, double baselineMs,
                                                     double enterDb, double exitDb)
{
    if (! geometry (channels)) return 0.0;
    felitronics::analysis::BandBurstsParams p;
    p.bandLowHz = bandLowHz; p.bandHighHz = bandHighHz;
    p.hopMs = hopMs; p.baselineMs = baselineMs;
    p.enterDb = enterDb; p.exitDb = exitDb;
    return demand (felitronics::analysis::BandBursts::storageFor (sampleRate, (int) channels, p));
}

FC_EXPORT double fc_probe_hum_storage_bytes (std::uint32_t channels, double sampleRate)
{
    if (! geometry (channels)) return 0.0;
    return demand (felitronics::analysis::HumDetector::storageFor (
                       sampleRate, (int) channels, kHumParams));
}

FC_EXPORT double fc_probe_forensics_storage_bytes (std::uint32_t channels, double sampleRate)
{
    if (! geometry (channels)) return 0.0;
    return demand (felitronics::analysis::SourceForensics::storageFor (
                       sampleRate, (int) channels, kFxParams));
}

//======================================================================================================
// K10 — excursions over a delivery ceiling, through analysis::PeakExcursions. The caller hands a render
// made WITHOUT the limiter and the ceiling it means to deliver at.
//
// EVERY AGGREGATE SURVIVES THE RUN LIST'S EXHAUSTION, so a page that asks only for scalars gets exact
// answers however long the programme is; the list is for coordinates, and `runs_complete` says whether it
// holds all of them.
namespace
{
    felitronics::analysis::PeakExcursions& excursions()
    {
        static felitronics::analysis::PeakExcursions e;
        return e;
    }
    bool haveExcursions = false;

    constexpr std::uint32_t kExScalars   = 22;
    constexpr std::uint32_t kExRunStride = 6;
    // What the last SUCCESSFUL run installed — never the compile-time default. The scalars publish the
    // ceiling and the window a caller must divide by, and a constant there becomes a lie the moment
    // `_run_with` exists: the crest section above learned that the expensive way.
    felitronics::analysis::PeakExcursions::Params installedEx {};
}

static int excursionsRunWith (const float* planar, std::uint32_t frames, std::uint32_t channels,
                              double sampleRate, const felitronics::analysis::PeakExcursions::Params& p)
{
    haveExcursions = false;
    if (! planarSpan (planar, frames, channels)) return 0;
    auto& d = excursions();
    d.setParams (p);
    if (! d.prepare (sampleRate, (int) channels)) return 0;
    const float* view[felitronics::core::kMaxChannels] {};
    // Guarded on `frames`: `planar + k * frames` is undefined behaviour when planar is null EVEN at a
    // zero offset, and an empty programme is allowed to arrive with one.
    if (frames != 0)
        for (std::uint32_t k = 0; k < channels; ++k) view[k] = planar + (std::size_t) k * (std::size_t) frames;
    if (frames != 0 && ! d.process (view, (int) channels, (int) frames)) return 0;
    if (! d.finish()) return 0;
    installedEx = p;
    haveExcursions = true;
    return 1;
}

FC_EXPORT int fc_probe_excursions_run (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                       double sampleRate, double thresholdDbtp)
{
    felitronics::analysis::PeakExcursions::Params p {};
    p.thresholdDbtp = thresholdDbtp;
    return excursionsRunWith (planar, frames, channels, sampleRate, p);
}

// The same measurement with the geometry a caller chose. Validated BY THE CORE, which is the only place
// that decides: a second opinion here would be a second definition. A refused run answers 0 and leaves
// the previous geometry rather than half of a new one.
FC_EXPORT int fc_probe_excursions_run_with (const float* planar, std::uint32_t frames, std::uint32_t channels,
                                            double sampleRate, double thresholdDbtp, double mergeMs,
                                            double e0, double e1, double e2, double e3, std::int32_t maxRuns)
{
    felitronics::analysis::PeakExcursions::Params p {};
    p.thresholdDbtp = thresholdDbtp;
    p.mergeMs       = mergeMs;
    p.classEdgesMs[0] = e0; p.classEdgesMs[1] = e1; p.classEdgesMs[2] = e2; p.classEdgesMs[3] = e3;
    p.maxRuns       = (int) maxRuns;
    return excursionsRunWith (planar, frames, channels, sampleRate, p);
}

FC_EXPORT std::uint32_t fc_probe_excursions_scalars_len (void) { return kExScalars; }
FC_EXPORT std::uint32_t fc_probe_excursions_run_stride  (void) { return kExRunStride; }
FC_EXPORT std::uint32_t fc_probe_excursions_classes     (void)
{ return (std::uint32_t) felitronics::analysis::PeakExcursions::kClasses; }
FC_EXPORT std::uint32_t fc_probe_excursions_crest_bins  (void)
{ return (std::uint32_t) felitronics::analysis::PeakExcursions::kCrestBins; }

FC_EXPORT std::uint32_t fc_probe_excursions_scalars (double* out, std::uint32_t cap)
{
    if (! haveExcursions || out == nullptr || cap < kExScalars || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = excursions();
    const auto& p = installedEx;
    std::uint32_t i = 0;
    out[i++] = d.sampleRate();              out[i++] = (double) d.channels();
    out[i++] = p.thresholdDbtp;             out[i++] = d.thresholdLinear();
    out[i++] = p.mergeMs;                   out[i++] = (double) d.samplesProcessed();
    out[i++] = (double) d.measuredOs();     out[i++] = (double) (int) d.reason();
    out[i++] = d.valid() ? 1.0 : 0.0;
    // The two peaks separately: the certificate is the larger of them, and an instrument publishing only
    // one would disagree with it without saying where.
    out[i++] = d.reconstructedPeak();       out[i++] = d.samplePeakLinear();
    out[i++] = d.truePeakLinear();
    out[i++] = (double) d.runCount();       out[i++] = (double) d.storedRunCount();
    out[i++] = d.runsComplete() ? 1.0 : 0.0;
    // Merge-free, both: a caller sweeping mergeMs keeps these two still while the count moves.
    out[i++] = (double) d.aboveOs();        out[i++] = d.occupancy();
    out[i++] = d.totalDose();               out[i++] = d.maxExcess();
    out[i++] = d.p90Ms();                   out[i++] = d.p90Saturated() ? 1.0 : 0.0;
    out[i++] = d.runsPerMinute();
    return i;
}

// One row per run: startOs, lengthOs, aboveOs, peak, dose, crestHz. Coordinates are OVERSAMPLED samples
// on the programme's own grid — the interpolator's 63.5-sample lag is already off, so `startOs / 4` is
// the input-sample position and the residual eighth of a sample is the FIR's half, not a rounding.
FC_EXPORT std::uint32_t fc_probe_excursions_runs (double* out, std::uint32_t cap)
{
    if (! haveExcursions || out == nullptr || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = excursions();
    const std::uint32_t room = cap / kExRunStride;
    std::uint32_t at = 0;
    for (std::int64_t r = 0; r < d.storedRunCount() && at < room; ++r, ++at)
    {
        const auto row = d.run (r);
        double* w = out + (std::size_t) at * kExRunStride;
        w[0] = (double) row.startOs;  w[1] = (double) row.lengthOs;  w[2] = (double) row.aboveOs;
        w[3] = row.peak;              w[4] = row.dose;
        w[5] = row.crestHz (d.sampleRate(), d.thresholdLinear());
    }
    return at;
}

// Per duration class: count then dose, in pairs. The classes are the installed edges plus "longer".
FC_EXPORT std::uint32_t fc_probe_excursions_classes_out (double* out, std::uint32_t cap)
{
    const std::uint32_t need = 2u * (std::uint32_t) felitronics::analysis::PeakExcursions::kClasses;
    if (! haveExcursions || out == nullptr || cap < need || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = excursions();
    std::uint32_t i = 0;
    for (int k = 0; k < felitronics::analysis::PeakExcursions::kClasses; ++k)
    {
        out[i++] = (double) d.classCount (k);
        out[i++] = d.classDose (k);
    }
    return i;
}

// The crest-frequency histogram, third-octaves from 20 Hz: low edge, count, dose, in triples. This is
// what replaces the `lowShare` that could not be measured — a consumer sweeps "is this bass" on its own
// corpus instead of receiving a boundary baked in at one frequency.
FC_EXPORT std::uint32_t fc_probe_excursions_crest (double* out, std::uint32_t cap)
{
    const std::uint32_t need = 3u * (std::uint32_t) felitronics::analysis::PeakExcursions::kCrestBins;
    if (! haveExcursions || out == nullptr || cap < need || ! outSpan (out, cap, 8)) return 0u;
    const auto& d = excursions();
    std::uint32_t i = 0;
    for (int b = 0; b < felitronics::analysis::PeakExcursions::kCrestBins; ++b)
    {
        out[i++] = felitronics::analysis::PeakExcursions::crestBinLowHz (b);
        out[i++] = (double) d.crestBinCount (b);
        out[i++] = d.crestBinDose (b);
    }
    return i;
}

// How tightly the peaks sit under the ceiling — the tell for a source that was true-peak LIMITED rather
// than clipped, whose maxima cluster with no flat tops for a clipping test to find. `minusDb` is the
// floor under the denominator; pass a large number for "every local maximum".
FC_EXPORT double fc_probe_excursions_ceiling_density (double minusDb, double withinDb)
{
    return haveExcursions ? excursions().ceilingDensityAbove (minusDb, withinDb) : 0.0;
}
FC_EXPORT double fc_probe_excursions_ceiling_maxima (void)
{
    return haveExcursions ? (double) excursions().ceilingMaxima() : 0.0;
}

FC_EXPORT double fc_probe_excursions_storage_bytes (std::uint32_t channels, double sampleRate,
                                                    double thresholdDbtp, double mergeMs, std::int32_t maxRuns)
{
    if (! geometry (channels)) return 0.0;
    felitronics::analysis::PeakExcursions::Params p {};
    p.thresholdDbtp = thresholdDbtp; p.mergeMs = mergeMs; p.maxRuns = (int) maxRuns;
    const auto st = felitronics::analysis::PeakExcursions::storageFor (sampleRate, (int) channels, p);
    return st.ok ? (double) st.bytes() : 0.0;
}

FC_EXPORT double fc_probe_lowend_storage_bytes (std::uint32_t channels, double sampleRate)
{
    if (! geometry (channels)) return 0.0;
    return demand (felitronics::analysis::LowEnd::storageFor (
                       sampleRate, (int) channels, kLeParams));
}

// The price of the geometry `_run_with` would take. Law 11d: the budget is the allocation, so a caller
// sizing a heap for a non-default order must be able to ask about THAT order, not about the default.
FC_EXPORT double fc_probe_lowend_storage_bytes_with (std::uint32_t channels, double sampleRate,
                                                     double crossoverHz, double lowNoteHz, double highNoteHz,
                                                     std::int32_t fftOrder, double dutyThresholdDb,
                                                     std::int32_t skipBlocks)
{
    if (! geometry (channels)) return 0.0;
    felitronics::analysis::LowEndParams lp = kLeParams;
    lp.crossoverHz = crossoverHz; lp.lowNoteHz = lowNoteHz; lp.highNoteHz = highNoteHz;
    lp.fftOrder = (int) fftOrder; lp.dutyThresholdDb = dutyThresholdDb; lp.skipBlocks = (int) skipBlocks;
    return demand (felitronics::analysis::LowEnd::storageFor (sampleRate, (int) channels, lp));
}

// HOW MANY 10 ms BLOCKS AN LR4 AT `crossoverHz` NEEDS TO FALL `dB` BELOW ITS OWN PEAK — the number
// `skipBlocks` wants, derived from the filter rather than written down on either side of the ABI. Exported
// because a consumer asked not to re-implement `u = 10.233` in its own language: a constant copied across
// a boundary is a second definition, and this one comes from a bisection on u*exp(1-u) = 10^(-dB/20).
// Pure: no instance, no state, nothing to prepare. 0 for arguments the class itself would refuse.
FC_EXPORT int fc_probe_lowend_settling_blocks (double sampleRate, double crossoverHz, double dB)
{
    return felitronics::analysis::LowEnd::settlingBlocks (sampleRate, crossoverHz, dB);
}

// The side fraction below a candidate crossover, from the band table — the SWEEP. One call per point, so
// a caller draws a curve. -1.0 — NOT 0.0 — for a frequency the crossover itself would refuse: 0.0 is a
// legitimate reading (a perfectly mono low end) and a consumer read one as the other. Read the core's note
// before using it as a number: outside [lowNoteHz, highNoteHz] the table is blind and the real filter is
// not, and on anti-phase content under 20 Hz the two give opposite answers.
FC_EXPORT double fc_probe_lowend_side_fraction_below (double hz)
{
    return haveLowEnd ? lowEnd().sideFractionBelow (hz) : 0.0;
}

// The lowest band present in at least `dutyMin` of the counted frames: band, midi, centreHz, count, duty,
// levelWhenOnDb, marginWhenOnDb. `band` is -1 when none qualifies, which is the core's own convention.
FC_EXPORT std::uint32_t fc_probe_lowend_lowest_occupied (double dutyMin, double* out, std::uint32_t cap)
{
    if (! haveLowEnd || out == nullptr || cap < kLeLowestFields || ! outSpan (out, cap, 8)) return 0u;
    const auto r = lowEnd().lowestOccupiedBand (dutyMin);
    std::uint32_t i = 0;
    out[i++] = (double) r.band;   out[i++] = (double) r.midi;  out[i++] = r.centreHz;
    out[i++] = (double) r.count;  out[i++] = r.duty;           out[i++] = r.levelWhenOnDb;
    out[i++] = r.marginWhenOnDb;
    return i;
}
