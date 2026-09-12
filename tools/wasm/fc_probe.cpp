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

#include "fcore_probe.h"

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
// already reports two peaks, `fc_probe_sample_peak` and `fc_probe_tp_linear` (the true peak of fcore::Probe's 128-tap
// REFERENCE filter — not analysis::TruePeakMeter's 48-tap spec filter, which the mastering chain runs: different
// filters, different numbers, the size of the gap is P62's to measure), and a waveform bucket is neither: a
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

// Build identity, so a mismatched artifact is obvious in a report rather than a mystery.
FC_EXPORT std::uint32_t fc_probe_os_factor       (void) { return (std::uint32_t) fcore::Probe::kOsFactor; }
FC_EXPORT std::uint32_t fc_probe_os_taps         (void) { return (std::uint32_t) fcore::Probe::kOsTapsPerPhase; }
FC_EXPORT std::uint32_t fc_probe_chunk           (void) { return (std::uint32_t) fcore::Probe::kChunk; }
FC_EXPORT std::uint32_t fc_probe_sizeof_longdouble (void) { return (std::uint32_t) sizeof (long double); }
