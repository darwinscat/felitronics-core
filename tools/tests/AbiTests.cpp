// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The C ABI in tools/wasm/fc_probe.cpp is the one surface that faces untrusted input: everything it receives
// was computed by JavaScript in a browser. Its guards therefore deserve tests, and they can have them —
// fc_probe.cpp compiles natively (the EMSCRIPTEN_KEEPALIVE macro degrades to a plain extern "C"), so the
// whole validation and addressing layer runs under ctest, ASan and UBSan like anything else.
//
// What is NOT covered here, because it needs a wasm heap: the detached-view discipline. That is exercised by
// tools/wasm/parity.mjs against real files.

#include <felitronics_test.h>

#include "fcore_probe.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

extern "C"
{
    int           fc_probe_run            (const float*, std::uint32_t, std::uint32_t, double);
    double        fc_probe_lufs           (const float*, std::uint32_t, std::uint32_t, double);
    double        fc_probe_dbtp           (const float*, std::uint32_t, std::uint32_t, double);
    double        fc_probe_tp_linear      (void);
    std::uint32_t fc_probe_block_count    (void);
    int           fc_probe_dropped        (void);
    std::uint32_t fc_probe_block_energies (double*, std::uint32_t);
    std::uint32_t fc_probe_os_factor      (void);
    std::uint32_t fc_probe_os_taps        (void);
    std::uint32_t fc_probe_chunk          (void);
    std::uint32_t fc_probe_sizeof_longdouble (void);
}

using namespace felitronics;

namespace
{
    // Planar-contiguous, the ABI's own layout: channel c occupies [c*frames, (c+1)*frames).
    std::vector<float> makePlanar (std::uint32_t frames, std::uint32_t channels, double sr, bool loudSecondChannel)
    {
        std::vector<float> v ((std::size_t) frames * channels, 0.0f);
        for (std::uint32_t c = 0; c < channels; ++c)
            for (std::uint32_t i = 0; i < frames; ++i)
            {
                const double amp = (c == 1 && loudSecondChannel) ? 0.5 : (c == 0 && loudSecondChannel ? 0.0 : 0.3);
                v[(std::size_t) c * frames + i] = (float) (amp * std::sin (2.0 * 3.14159265358979 * 440.0 * i / sr));
            }
        return v;
    }
}

int main()
{
    std::printf ("fc_probe — the wasm C ABI, exercised natively\n");
    const double sr = 48000.0;
    const std::uint32_t frames = 48000 * 3;

    test::group ("build identity is reported");
    {
        test::ok (fc_probe_os_factor() == 4 && fc_probe_os_taps() == 32, "the reference true-peak config");
        test::ok (fc_probe_chunk() == (std::uint32_t) fcore::Probe::kChunk, "the streaming step");
        test::ok (fc_probe_sizeof_longdouble() == sizeof (long double), "sizeof(long double) for this target");
    }

    // --- The ABI must produce EXACTLY what the class produces. If the planar addressing or the channel loop
    //     were wrong, the number would merely be a bit off — which is the hardest kind of wrong to notice. ---
    test::group ("ABI null: a run through the C entry point equals a direct fcore::Probe run");
    {
        const auto buf = makePlanar (frames, 2, sr, false);
        test::ok (fc_probe_run (buf.data(), frames, 2, sr) == 1, "the run was accepted");

        fcore::Probe direct;
        direct.prepare (sr, 2);
        const float* view[2] { buf.data(), buf.data() + frames };
        direct.process (view, 2, (long long) frames);

        test::ok (fc_probe_block_count() == (std::uint32_t) direct.gatingBlockCount(), "same block count");
        test::ok (fc_probe_tp_linear() == direct.truePeakLinear(), "same true-peak maximum, bit-for-bit");

        std::vector<double> got (fc_probe_block_count());
        const std::uint32_t n = fc_probe_block_energies (got.data(), (std::uint32_t) got.size());
        test::ok (n == got.size(), "every energy was copied out");
        bool same = true;
        for (std::uint32_t i = 0; i < n; ++i) if (got[i] != direct.gatingBlockEnergies()[i]) same = false;
        test::ok (same, "every pre-gate block energy is bit-identical to the class's own");
    }

    // --- Channel c really must be read from planar + c*frames. A shim that read channel 0 twice would give a
    //     plausible number on ordinary stereo material and be silently wrong. ---
    test::group ("planar addressing: the second channel is actually read");
    {
        const auto silentL = makePlanar (frames, 2, sr, true);      // ch0 silent, ch1 loud
        test::ok (fc_probe_run (silentL.data(), frames, 2, sr) == 1, "accepted");
        const double bothCh = fc_probe_lufs (silentL.data(), frames, 2, sr);

        // The same buffer read as MONO sees only channel 0 — which is silence.
        const double onlyCh0 = fc_probe_lufs (silentL.data(), frames, 1, sr);
        test::ok (onlyCh0 <= -100.0, "channel 0 alone is silence, as constructed");
        test::ok (bothCh > -60.0, "but the stereo run hears channel 1 — so c*frames addressing is right");
    }

    // --- Everything below arrives from JavaScript. ---
    test::group ("the ABI rejects what a page can hand it");
    {
        const auto buf = makePlanar (1024, 2, sr, false);
        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();

        test::ok (fc_probe_run (nullptr, 1024, 2, sr) == 0, "null pointer");
        test::ok (fc_probe_run (buf.data(), 0, 2, sr) == 0, "zero frames");
        test::ok (fc_probe_run (buf.data(), 1024, 0, sr) == 0, "zero channels");
        test::ok (fc_probe_run (buf.data(), 1024, 99, sr) == 0, "channels past kMaxChannels");
        test::ok (fc_probe_run (buf.data(), 1024, 2, 0.0) == 0, "sample rate 0");
        test::ok (fc_probe_run (buf.data(), 1024, 2, -48000.0) == 0, "negative sample rate");
        test::ok (fc_probe_run (buf.data(), 1024, 2, nan) == 0, "NaN sample rate");
        test::ok (fc_probe_run (buf.data(), 1024, 2, inf) == 0, "+inf sample rate (a naive `> 0` check lets it through)");
        test::ok (fc_probe_run (buf.data(), 1024, 2, 1e300) == 0, "an absurd but finite rate");
        test::ok (fc_probe_run (buf.data(), 1024, 2, 5e-324) == 0, "the smallest positive subnormal double");
        test::ok (fc_probe_run (buf.data(), 1024, 2, fcore::Probe::kMinSampleRate - 1.0) == 0, "below the accepted range");
        test::ok (fc_probe_run (buf.data(), 1024, 2, fcore::Probe::kMaxSampleRate + 1.0) == 0, "above the accepted range");

        // A misaligned float* reads garbage in a release wasm build and only traps under -sSAFE_HEAP.
        const char* raw = reinterpret_cast<const char*> (buf.data());
        const float* skewed = reinterpret_cast<const float*> (raw + 1);
        test::ok (fc_probe_run (skewed, 16, 1, sr) == 0, "a pointer that is not 4-byte aligned");

        // frames*channels*sizeof(float) must fit a 32-bit address space — on wasm32 the product IS the
        // allocation, so a wrapped one would hand the core a window onto unrelated heap.
        test::ok (fc_probe_run (buf.data(), 0xFFFFFFFFu, 2, sr) == 0, "a frames x channels product that overflows 32 bits");
    }

    test::group ("a rejected call CLEARS the previous result (the documented contract)");
    {
        const auto buf = makePlanar (frames, 2, sr, false);
        test::ok (fc_probe_run (buf.data(), frames, 2, sr) == 1, "a good run first");
        test::ok (fc_probe_block_count() > 0, "which produced blocks");
        test::ok (fc_probe_run (buf.data(), frames, 2, -1.0) == 0, "then a rejected one");
        test::ok (fc_probe_block_count() == 0, "the getters now read zero, NOT the previous run");
        test::ok (fc_probe_tp_linear() == 0.0, "true peak cleared too");
    }

    test::group ("fc_probe_block_energies respects its output capacity");
    {
        const auto buf = makePlanar (frames, 1, sr, false);
        test::ok (fc_probe_run (buf.data(), frames, 1, sr) == 1, "a run to read from");
        const std::uint32_t n = fc_probe_block_count();
        test::ok (n > 4, "enough blocks to truncate");

        std::vector<double> small (3, -1.0);
        test::ok (fc_probe_block_energies (small.data(), 3) == 3, "a small buffer is filled, not overrun");
        test::ok (small[0] != -1.0 && small[2] != -1.0, "and actually written");

        std::vector<double> big (n + 10, -1.0);
        test::ok (fc_probe_block_energies (big.data(), (std::uint32_t) big.size()) == n, "a large buffer gets exactly n");
        test::ok (big[n] == -1.0, "and nothing past n is touched");

        test::ok (fc_probe_block_energies (nullptr, 10) == 0, "a null output buffer writes nothing");
        test::ok (fc_probe_block_energies (big.data(), 0) == 0, "a zero capacity writes nothing");
    }

    test::group ("the scalar entry points agree with the surface they are derived from");
    {
        const auto buf = makePlanar (frames, 2, sr, false);
        const double lufs = fc_probe_lufs (buf.data(), frames, 2, sr);
        const double tpDb = fc_probe_dbtp (buf.data(), frames, 2, sr);
        test::ok (std::isfinite (lufs) && lufs < 0.0, "a sane loudness");
        test::approx (tpDb, 20.0 * std::log10 (fc_probe_tp_linear()), 1e-12, "dBTP is the dB form of the linear max");
        test::ok (fc_probe_dropped() == 0, "nothing dropped at this length");
    }

    test::group ("a rejected scalar call returns the documented sentinel");
    {
        test::ok (fc_probe_lufs (nullptr, 10, 2, 48000.0) == -120.0, "lufs sentinel");
        test::approx (fc_probe_dbtp (nullptr, 10, 2, 48000.0), 20.0 * std::log10 (1e-9), 1e-12, "dBTP floor");
    }

    return test::report();
}
