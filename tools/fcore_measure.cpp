// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_measure — a tiny streaming measurement CLI for validating felitronics-core's analysis against a
// reference (ffmpeg), and the NATIVE SIDE of the P0 wasm parity check. Reads interleaved 32-bit-float
// little-endian PCM (as ffmpeg emits with `-f f32le`).
//
// The measurement itself lives in fcore_probe.h, shared verbatim with the wasm shim — read that header for
// the build contract (-ffp-contract=off here, no -mrelaxed-simd there) and for why the true-peak filter
// config is part of that contract rather than an implementation detail.
//
//   lufs        → integrated loudness (LUFS)            ↔ ffmpeg ebur128 "I:"
//   truepeak    → max true peak (dBTP, 4× oversampled)  ↔ ffmpeg ebur128 "Peak:" (True Peak)
//   correlation → whole-file Pearson L/R correlation    (synthetic-validated; no clean ffmpeg single number)
//   blocks      → the CROSS-TOOLCHAIN SURFACE: every pre-gate 400 ms gating-block energy plus the true-peak
//                 linear maximum, as raw IEEE-754 bit patterns. Diffing two `blocks` outputs IS the parity
//                 test — see the note at the mode itself for why the gated scalars cannot be that test.
//
// Usage: fcore_measure <lufs|truepeak|correlation|blocks> <sampleRate> <channels> <raw.f32le> [--precise]
//
// The scalar modes print %.2f by default (tools/validate_ffmpeg.sh compares against ffmpeg's own two
// decimals); `--precise` switches them to %.17g plus the exact %a form. `blocks` is always exact.

#include "fcore_probe.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;

namespace
{
    constexpr int kChunk = fcore::Probe::kChunk;

    // Reads the file in kChunk-frame steps, de-interleaving into planar scratch and handing each step to
    // `sink`. The de-interleave is a pure float permutation — exact, and identical to what the wasm side's
    // caller does with an already-planar heap buffer.
    template <typename Sink>
    bool streamPlanar (std::FILE* f, int nc, Sink&& sink)
    {
        std::vector<float>              inter ((std::size_t) kChunk * (std::size_t) nc);
        std::vector<std::vector<float>> ch ((std::size_t) nc, std::vector<float> ((std::size_t) kChunk));
        std::vector<const float*>       cp ((std::size_t) nc);

        std::size_t got;
        while ((got = std::fread (inter.data(), sizeof (float), (std::size_t) kChunk * (std::size_t) nc, f)) > 0)
        {
            const int frames = (int) (got / (std::size_t) nc);
            if (frames <= 0) break;
            for (int i = 0; i < frames; ++i)
                for (int c = 0; c < nc; ++c)
                    ch[(std::size_t) c][(std::size_t) i] = inter[(std::size_t) (i * nc + c)];
            for (int c = 0; c < nc; ++c) cp[(std::size_t) c] = ch[(std::size_t) c].data();
            sink (cp.data(), frames);
        }
        return true;
    }

    std::uint64_t bits (double d) noexcept
    {
        std::uint64_t u;
        std::memcpy (&u, &d, sizeof u);
        return u;
    }

    void printScalar (double v, bool precise, const char* unit)
    {
        if (precise) std::printf ("%.17g  %a  %s\n", v, v, unit);
        else         std::printf ("%.2f\n", v);
    }
}

int main (int argc, char** argv)
{
    if (argc < 5)
    {
        std::fprintf (stderr,
            "usage: %s <lufs|truepeak|correlation|blocks> <sampleRate> <channels> <raw.f32le> [--precise]\n",
            argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    const double fs = std::atof (argv[2]);
    const int    nc = std::atoi (argv[3]);
    const bool precise = (argc > 5 && std::strcmp (argv[5], "--precise") == 0) || std::getenv ("FCORE_PRECISE") != nullptr;

    if (nc < 1 || nc > core::kMaxChannels || ! (fs > 0.0) || ! std::isfinite (fs))
    {
        std::fprintf (stderr, "bad sampleRate/channels (channels 1..%d, sampleRate finite and > 0)\n", core::kMaxChannels);
        return 2;
    }

    std::FILE* f = std::fopen (argv[4], "rb");
    if (! f) { std::perror ("open"); return 2; }

    if (mode == "correlation")
    {
        // NB: accumulates in `long double`, which is 8 bytes on arm64 macOS, 16 (IEEE quad) on wasm32 and
        // 80-bit x87 on x86-64 Linux — three tiers, three answers. Fine for this native-only sanity number;
        // it is NOT a cross-tier comparison surface, and a facade must not inherit the type.
        long double sumLR = 0.0L, sumLL = 0.0L, sumRR = 0.0L;
        streamPlanar (f, nc, [&] (const float* const* p, int n)
        {
            for (int i = 0; i < n; ++i)
            {
                const long double l = p[0][i];
                const long double r = nc > 1 ? p[1][i] : l;
                sumLR += l * r; sumLL += l * l; sumRR += r * r;
            }
        });
        std::fclose (f);
        const long double d = std::sqrt (sumLL * sumRR);
        const double corr = d > 1e-12L ? (double) std::clamp (sumLR / d, -1.0L, 1.0L) : 1.0;
        if (precise) std::printf ("%.17g  %a\n", corr, corr);
        else         std::printf ("%.3f\n", corr);
        return 0;
    }

    fcore::Probe probe;
    if (! probe.prepare (fs, nc))
    {
        std::fprintf (stderr, "probe.prepare failed\n");
        std::fclose (f);
        return 2;
    }
    streamPlanar (f, nc, [&] (const float* const* p, int n) { probe.process (p, nc, n); });
    std::fclose (f);
    probe.finish();     // drains the polyphase FIR — without it a peak in the final samples is not measured

    if (probe.droppedBlocks() != 0)
        std::fprintf (stderr, "warning: %d gating blocks dropped — program longer than the meter's capacity\n",
                      probe.droppedBlocks());

    if (mode == "lufs")          printScalar (probe.integratedLufs(), precise, "LUFS");
    else if (mode == "truepeak") printScalar (probe.truePeakDb(),     precise, "dBTP");
    else if (mode == "blocks")
    {
        // The parity surface. Everything here is CONTINUOUS in the input samples, so a bit-for-bit diff
        // between two toolchains is a meaningful equivalence test; the gated LUFS printed by `lufs` is not,
        // because a block within ~1e-12 of a gate flips inclusion and moves it by ~0.01 dB.
        // Emitted as raw IEEE-754 bit patterns, not %a and not decimal: the other side of this comparison is
        // JavaScript, which has no hex-float printing and whose decimal formatting is not C's. A 16-hex-digit
        // pattern is the one representation both sides produce identically, so `diff` IS the parity test.
        const auto e = probe.gatingBlockEnergies();
        const int    n = (int) e.size();
        // The rate goes out as a bit pattern too, not as text: C's %g and JavaScript's Number-to-string do
        // not agree on a fractional rate (48000.123456 prints as 48000.12346 at ten significant digits), so a
        // decimal header would break the whole-file diff while every measured bit matched.
        std::printf ("# fcore blocks v2 sr=%016llx ch=%d os=%dx%d chunk=%d\n",
                     (unsigned long long) bits (fs), nc, fcore::Probe::kOsFactor, fcore::Probe::kOsTapsPerPhase, kChunk);
        std::printf ("tp %016llx\n", (unsigned long long) bits (probe.truePeakLinear()));
        std::printf ("sp %016llx\n", (unsigned long long) bits (probe.samplePeakLinear()));
        std::printf ("blocks %d\n", n);
        for (int j = 0; j < n; ++j) std::printf ("%016llx\n", (unsigned long long) bits (e[(std::size_t) j]));
    }
    else { std::fprintf (stderr, "unknown mode '%s'\n", mode.c_str()); return 2; }

    return 0;
}
