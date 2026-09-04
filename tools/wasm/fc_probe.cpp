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

    // Ptr/size validation the core cannot do for us: everything below arrives from JS.
    bool viable (const float* planar, std::uint32_t frames, std::uint32_t channels, double sampleRate)
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
        if (! inHeap (planar, bytes)) return false;
        // rejects 0 / negative / NaN / +inf and the absurd-but-finite rates (see Probe::kMinSampleRate)
        return probe().prepare (sampleRate, (int) channels);
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
    const std::uint32_t n = (std::uint32_t) probe().gatingBlockCount();
    const std::uint32_t m = n < cap ? n : cap;
    if (! inHeap (out, (std::uint64_t) m * sizeof (double))) return 0;   // the output span must fit too
    if (m > 0) std::memcpy (out, probe().gatingBlockEnergies().data(), (std::size_t) m * sizeof (double));
    return m;
}

// Build identity, so a mismatched artifact is obvious in a report rather than a mystery.
FC_EXPORT std::uint32_t fc_probe_os_factor       (void) { return (std::uint32_t) fcore::Probe::kOsFactor; }
FC_EXPORT std::uint32_t fc_probe_os_taps         (void) { return (std::uint32_t) fcore::Probe::kOsTapsPerPhase; }
FC_EXPORT std::uint32_t fc_probe_chunk           (void) { return (std::uint32_t) fcore::Probe::kChunk; }
FC_EXPORT std::uint32_t fc_probe_sizeof_longdouble (void) { return (std::uint32_t) sizeof (long double); }
