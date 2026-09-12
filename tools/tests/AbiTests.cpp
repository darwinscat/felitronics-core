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
#include <cstring>
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

    int           fc_probe_shapes_run          (const float*, std::uint32_t, std::uint32_t, double, std::uint32_t, std::int32_t, std::uint32_t);
    std::uint32_t fc_probe_waveform_count      (void);
    std::uint32_t fc_probe_waveform_emitted    (void);
    std::uint32_t fc_probe_waveform_decimation (void);
    std::uint32_t fc_probe_waveform_peaks      (double*, std::uint32_t);
    std::uint32_t fc_probe_waveform_peaks_f32  (float*, std::uint32_t);
    std::uint32_t fc_probe_stereo_cols         (void);
    int           fc_probe_stereo_is_mono      (void);
    double        fc_probe_stereo_max_rms      (void);
    std::uint32_t fc_probe_stereo_width        (float*, std::uint32_t);
    std::uint32_t fc_probe_stereo_corr         (float*, std::uint32_t);
    std::uint32_t fc_probe_stereo_rms          (float*, std::uint32_t);
    int           fc_probe_needle              (const float*, std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t, double*);
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

        // The OUTPUT is alignment-checked too, and it used not to be. The input side of this same file
        // has refused a misaligned pointer since P0; the output side checked only the address range,
        // which is the same rule applied in one direction. A page reads this span back as
        // `new Float64Array(HEAPF64.buffer, ptr, n)` and that THROWS on a ptr that is not a multiple of
        // 8, so the refusal replaces an exception in somebody's worker with a zero the caller can test.
        const char* rawOut = reinterpret_cast<const char*> (big.data());
        double* skewedOut = reinterpret_cast<double*> (const_cast<char*> (rawOut) + 1);
        test::ok (fc_probe_block_energies (skewedOut, 4) == 0, "a misaligned output buffer writes nothing");
        test::ok (big[0] != -1.0, "PRECONDITION: the aligned call above really did write, so this is a live check");
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

    // --- P59a: the waveform peaks and the stereo band. The same two questions as the loudness path: does the C
    //     entry point produce EXACTLY the class's bits, and does it refuse what a page can hand it. ---
    test::group ("shapes ABI null: fc_probe_shapes_run equals a direct fcore::ShapeProbe, bit for bit");
    {
        const std::uint32_t n = 100003;
        auto buf = makePlanar (n, 2, sr, false);
        for (std::uint32_t i = 0; i < n; i += 97) buf[(std::size_t) i + n] = -buf[i];   // some anti-phase
        test::ok (fc_probe_shapes_run (buf.data(), n, 2, 44100.0, 1100, 3, 1100) == 1, "the run was accepted");

        fcore::ShapeProbe direct;
        const float* view[2] { buf.data(), buf.data() + n };
        test::ok (direct.prepare (44100.0, 2, n, 1100, analysis::PeakMix::Max, 1100) && direct.process (view, 2, n) && direct.complete(),
                  "the direct run completed");

        std::vector<double> pk (1100); std::vector<float> pk32 (1100), w (1100), c (1100), l (1100);
        test::ok (fc_probe_waveform_count() == 1100 && fc_probe_waveform_peaks (pk.data(), 1100) == 1100 && fc_probe_waveform_peaks_f32 (pk32.data(), 1100) == 1100,
                  "1100 peaks, both forms");
        test::ok (std::memcmp (pk.data(), direct.peaks().peaks().data(), 1100 * sizeof (double)) == 0, "the double peaks are the class's bits");
        bool f32same = true;
        for (int i = 0; i < 1100; ++i) { const float want = direct.peaks().peakAsFloat32 (i); f32same = f32same && std::memcmp (&pk32[(std::size_t) i], &want, 4) == 0; }
        test::ok (f32same, "the float32 peaks are the class's float32 form");
        test::ok (fc_probe_waveform_emitted() == (std::uint32_t) direct.peaks().bucketsEmitted()
               && fc_probe_waveform_decimation() == (std::uint32_t) direct.peaks().decimation(), "emitted count and decimation forwarded");

        test::ok (fc_probe_stereo_cols() == 1100 && fc_probe_stereo_is_mono() == 0, "1100 columns, not mono");
        test::ok (fc_probe_stereo_width (w.data(), 1100) == 1100 && fc_probe_stereo_corr (c.data(), 1100) == 1100
               && fc_probe_stereo_rms (l.data(), 1100) == 1100, "three float32 arrays copied");
        test::ok (std::memcmp (w.data(), direct.stereo().width().data(), 1100 * 4) == 0
               && std::memcmp (c.data(), direct.stereo().correlation().data(), 1100 * 4) == 0
               && std::memcmp (l.data(), direct.stereo().rms().data(), 1100 * 4) == 0, "width / corr / loud are the class's bits");
        const double ml = fc_probe_stereo_max_rms(), dml = direct.stereo().maxRms();
        test::ok (std::memcmp (&ml, &dml, 8) == 0, "maxLoud is the class's double");

        double out3[3] {};
        analysis::StereoColumns::Needle nd;
        test::ok (fc_probe_needle (buf.data(), n, 2, 1000, 50000, out3) == 1
               && analysis::StereoColumns::needle (buf.data(), buf.data() + n, n, 1000, 50000, nd)
               && std::memcmp (&out3[0], &nd.correlation, 8) == 0 && std::memcmp (&out3[1], &nd.width, 8) == 0
               && std::memcmp (&out3[2], &nd.rms, 8) == 0, "the needle through the ABI is the class's three doubles");
        test::ok (fc_probe_needle (buf.data(), n, 1, 0, n, out3) == 1 && out3[0] == 1.0 && out3[1] == 0.0,
                  "a mono needle reads channel 0 against itself: corr +1, width 0");
    }

    test::group ("fcore::ShapeProbe refuses a call it cannot honour before anything moves (law 11)");
    {
        // The review round found both: a null table was dereferenced before either class could refuse it, and a call
        // one frame too long consumed its first 8192-frame chunk before the second chunk was refused — after which
        // even the correct call failed.
        const std::uint32_t n = 10000;
        const auto buf = makePlanar (n, 1, sr, false);
        const float* view[1] { buf.data() };
        fcore::ShapeProbe sp;
        test::ok (sp.prepare (8000.0, 1, n, 1, analysis::PeakMix::Left, 1), "prepared for 10000 frames");
        test::ok (! sp.process (nullptr, 1, 1), "a null plane table is refused");
        const float* nullPlane[1] { nullptr };
        test::ok (! sp.process (nullPlane, 1, 1), "a null plane is refused");
        test::ok (! sp.process (view, 1, (long long) n + 1), "10001 frames into 10000 is refused");
        test::ok (sp.peaks().framesSeen() == 0 && sp.stereo().framesSeen() == 0, "... and nothing moved: both counters still 0");
        test::ok (sp.process (view, 1, (long long) n) && sp.complete(), "the correct call then succeeds");
    }

    test::group ("the shapes ABI refuses what a page can hand it, and a refusal clears the previous result");
    {
        const std::uint32_t n = 4800;
        auto buf = makePlanar (n, 2, sr, false);
        test::ok (fc_probe_shapes_run (buf.data(), n, 2, sr, 1000, 0, 1100) == 1, "a good run first");
        test::ok (fc_probe_shapes_run (nullptr, n, 2, sr, 1000, 0, 1100) == 0, "null planes");
        test::ok (fc_probe_waveform_count() == 0 && fc_probe_stereo_cols() == 0 && fc_probe_stereo_max_rms() == 0.0,
                  "... and the previous result is gone, not served");
        test::ok (fc_probe_shapes_run (buf.data(), 0, 2, sr, 1000, 0, 1100) == 0, "zero frames");
        test::ok (fc_probe_shapes_run (buf.data(), n, 0, sr, 1000, 0, 1100) == 0
               && fc_probe_shapes_run (buf.data(), n, (std::uint32_t) core::kMaxChannels + 1, sr, 1000, 0, 1100) == 0, "a width the core does not have");
        test::ok (fc_probe_shapes_run (buf.data() + 1, n - 1, 2, sr, 1000, 0, 1100) == 0 || (reinterpret_cast<std::uintptr_t> (buf.data() + 1) & 3u) == 0,
                  "a misaligned float* (when the offset really is misaligned)");
        test::ok (fc_probe_shapes_run (buf.data(), n, 2, sr, 1000, -1, 1100) == 0 && fc_probe_shapes_run (buf.data(), n, 2, sr, 1000, 4, 1100) == 0,
                  "a mix code that names nothing is refused, not clamped");
        test::ok (fc_probe_shapes_run (buf.data(), n, 2, sr, 0, 0, 1100) == 0 && fc_probe_shapes_run (buf.data(), n, 2, sr, 1000, 0, 0) == 0
               && fc_probe_shapes_run (buf.data(), n, 2, sr, 0x80000000u, 0, 1100) == 0, "0 buckets, 0 columns, a count past int");
        test::ok (fc_probe_shapes_run (buf.data(), n, 2, std::numeric_limits<double>::quiet_NaN(), 1000, 0, 1100) == 0
               && fc_probe_shapes_run (buf.data(), n, 2, 0.0, 1000, 0, 1100) == 0, "a rate that is not a rate");

        test::ok (fc_probe_shapes_run (buf.data(), n, 2, sr, 1000, 0, 1100) == 1, "a good run again");
        std::vector<double> pk (1000, -7.0);
        test::ok (fc_probe_waveform_peaks (pk.data(), 10) == 10 && pk[10] == -7.0, "the capacity binds: 10 asked, 10 written, the 11th untouched");
        test::ok (fc_probe_waveform_peaks (nullptr, 10) == 0 && fc_probe_waveform_peaks (pk.data(), 0) == 0, "null or zero capacity writes nothing");
        std::vector<double> big (3);
        auto* odd = reinterpret_cast<double*> (reinterpret_cast<char*> (big.data()) + 1);
        test::ok (fc_probe_waveform_peaks (odd, 1) == 0, "a misaligned double* output is refused");
        double out3[3] {};
        test::ok (fc_probe_needle (buf.data(), n, 2, 10, 9, out3) == 0 && fc_probe_needle (buf.data(), n, 2, 0, n + 1, out3) == 0
               && fc_probe_needle (buf.data(), n, 2, 0, n, nullptr) == 0, "the needle refuses from > to, a stretch past the planes, a null output");
        test::ok (fc_probe_waveform_count() == 1000, "a needle call neither reads nor clears the shapes result");
    }

    return test::report();
}
