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
//   correlation → whole-file normalised L/R correlation ΣLR/√(ΣLL·ΣRR) — the stereo band's own formula
//                 (analysis::StereoSums) over the whole file. Not Pearson's: nothing is centred.
//   blocks      → the CROSS-TOOLCHAIN SURFACE: every pre-gate 400 ms gating-block energy plus the true-peak
//                 linear maximum, as raw IEEE-754 bit patterns. Diffing two `blocks` outputs IS the parity
//                 test — see the note at the mode itself for why the gated scalars cannot be that test.
//   waveform    → the waveform peaks (analysis::WaveformPeaks): every bucket as a double AND as its float32
//                 form, bit patterns. A box-averaged max-abs, not above the sample peak except by rounding — NOT `truepeak`, which
//                 is fcore::Probe's reference true peak. [--buckets N] [--mix avr|L|R|max]
//   stereo      → the stereo band (analysis::StereoColumns): per column width / correlation / RMS as float32
//                 bit patterns, plus the maximum RMS as a double. RMS, not `lufs`. [--columns N]
//   needle      → correlation / width / RMS over [from, to) as double bit patterns. --from A --to B
//
// Usage: fcore_measure <mode> <sampleRate> <channels> <raw.f32le> [--precise] [mode options]
//
// The scalar modes print %.2f by default (tools/validate_ffmpeg.sh compares against ffmpeg's own two
// decimals); `--precise` switches them to %.17g plus the exact %a form. `blocks`, `waveform`, `stereo` and
// `needle` are always exact, and the node side (tools/wasm/shapes-parity.mjs) prints the same bytes.
//
// DECODE TO FLOAT, NEVER TO s16. This tool reads f32le: `ffmpeg -i x -f f32le out.f32`. The waveform, stereo and needle
// modes size the file before reading it (every boundary depends on the length), so they need a seekable file, not a pipe. An integer decode
// (`-f s16le`, as the old sidecar generator did) clamps a lossy file's samples above 0 dBFS and quantises the
// rest, and the peaks of that are not the peaks of the file.

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

    std::uint32_t bits32 (float x) noexcept
    {
        std::uint32_t u;
        std::memcpy (&u, &x, sizeof u);
        return u;
    }

    void printScalar (double v, bool precise, const char* unit)
    {
        if (precise) std::printf ("%.17g  %a  %s\n", v, v, unit);
        else         std::printf ("%.2f\n", v);
    }

    // The shape modes' options, parsed ONCE and strictly: every argument after the file is `--precise` or a known
    // `--name value` pair, each name at most once. An unknown or misspelt flag (`--bucket 2`) is refused, not ignored,
    // and a repeated one is refused rather than resolved — the node side of the parity check could not be relied on
    // to resolve a repeat the same way. (The older modes keep their old, lenient reading; nothing here changes them.)
    bool shapeOptions (int argc, char** argv, std::string& buckets, std::string& columns, std::string& mix,
                       std::string& from, std::string& to)
    {
        buckets = "1000"; columns = "1200"; mix = "avr"; from = "0"; to = "";
        bool seen[5] {};
        for (int i = 5; i < argc; ++i)
        {
            if (std::strcmp (argv[i], "--precise") == 0) continue;
            static constexpr const char* names[5] { "--buckets", "--columns", "--mix", "--from", "--to" };
            int k = 0;
            while (k < 5 && std::strcmp (argv[i], names[k]) != 0) ++k;
            if (k == 5 || seen[k] || i + 1 >= argc) return false;
            seen[k] = true;
            std::string& dst = k == 0 ? buckets : k == 1 ? columns : k == 2 ? mix : k == 3 ? from : to;
            dst = argv[++i];
        }
        return true;
    }

    // A finite, positive rate written as a plain decimal number and nothing else — `atof("8000Hz")` reads 8000.
    bool parseRate (const char* s, double& out)
    {
        char* end = nullptr;
        const double v = std::strtod (s, &end);
        if (end == s || *end != '\0' || ! (v > 0.0) || ! std::isfinite (v)) return false;
        out = v;
        return true;
    }

    // A whole non-negative decimal integer, nothing else — `atoi("12x")` would read 12.
    bool parseCount (const std::string& s, std::uint64_t& out)
    {
        if (s.empty() || s.size() > 18) return false;
        std::uint64_t v = 0;
        for (char ch : s) { if (ch < '0' || ch > '9') return false; v = v * 10 + (std::uint64_t) (ch - '0'); }
        out = v;
        return true;
    }

    // Frames in the file: its size over one interleaved frame. The size is read through the 64-bit file position
    // (`long` is 32 bits on Windows, and a stereo file past 2 GiB is not exotic), and a size that is not a whole
    // number of frames is REFUSED: a trailing partial frame is audio the measurement would silently not see.
    bool fileFrames (std::FILE* f, int nc, std::uint64_t& out)
    {
#if defined(_WIN32)
        if (_fseeki64 (f, 0, SEEK_END) != 0) return false;
        const long long size = _ftelli64 (f);
        if (size < 0 || _fseeki64 (f, 0, SEEK_SET) != 0) return false;
#else
        if (fseeko (f, 0, SEEK_END) != 0) return false;
        const long long size = (long long) ftello (f);
        if (size < 0 || fseeko (f, 0, SEEK_SET) != 0) return false;
#endif
        const std::uint64_t frameBytes = (std::uint64_t) nc * sizeof (float);
        if ((std::uint64_t) size % frameBytes != 0) return false;
        out = (std::uint64_t) size / frameBytes;
        return true;
    }
}

int main (int argc, char** argv)
{
    if (argc < 5)
    {
        std::fprintf (stderr,
            "usage: %s <lufs|truepeak|correlation|blocks|waveform|stereo|needle> <sampleRate> <channels> <raw.f32le>\n"
            "          [--precise] [--buckets N] [--mix avr|L|R|max] [--columns N] [--from A --to B]\n",
            argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    const double fs = std::atof (argv[2]);
    const int    nc = std::atoi (argv[3]);
    bool precise = std::getenv ("FCORE_PRECISE") != nullptr;
    for (int i = 5; i < argc; ++i) if (std::strcmp (argv[i], "--precise") == 0) precise = true;

    if (nc < 1 || nc > core::kMaxChannels || ! (fs > 0.0) || ! std::isfinite (fs))
    {
        std::fprintf (stderr, "bad sampleRate/channels (channels 1..%d, sampleRate finite and > 0)\n", core::kMaxChannels);
        return 2;
    }

    std::FILE* f = std::fopen (argv[4], "rb");
    if (! f) { std::perror ("open"); return 2; }

    if (mode == "correlation")
    {
        // The stereo band's formula over the whole file, streamed — the same StereoSums a column and the needle
        // accumulate, in binary64. (This mode used to accumulate in `long double`, a third definition of the
        // number and law 9's one sanctioned exception; neither is left.)
        analysis::StereoSums sums;
        streamPlanar (f, nc, [&] (const float* const* p, int n)
        {
            const float* L = p[0];
            const float* R = nc > 1 ? p[1] : p[0];
            for (int i = 0; i < n; ++i) sums.add (L[i], R[i]);
        });
        std::fclose (f);
        const double corr = sums.correlation();
        if (precise) std::printf ("%.17g  %a\n", corr, corr);
        else         std::printf ("%.3f\n", corr);
        return 0;
    }

    if (mode == "waveform" || mode == "stereo" || mode == "needle")
    {
        // Strict where the older modes are lenient: the rate and the width are checked as whole numbers here, since
        // `atoi("4294967297")` narrows to 1 channel and `atof("8000Hz")` reads 8000, and both would measure silently.
        double rate = 0.0; std::uint64_t width = 0;
        if (! parseRate (argv[2], rate) || ! parseCount (argv[3], width) || width < 1 || width > (std::uint64_t) core::kMaxChannels)
        {
            std::fprintf (stderr, "bad sampleRate/channels\n");
            std::fclose (f);
            return 2;
        }
        std::uint64_t frames = 0;
        if (! fileFrames (f, nc, frames) || frames == 0)
        {
            std::fprintf (stderr, "cannot size the file, it is not a whole number of %d-channel float32 frames, or it is empty\n", nc);
            std::fclose (f);
            return 2;
        }
        std::string sBuckets, sMix, sColumns, sFrom, sTo;
        std::uint64_t buckets = 0, columns = 0, from = 0, to = 0;
        if (! shapeOptions (argc, argv, sBuckets, sColumns, sMix, sFrom, sTo)
         || ! parseCount (sBuckets, buckets) || ! parseCount (sColumns, columns) || ! parseCount (sFrom, from)
         || (! sTo.empty() && ! parseCount (sTo, to)))
        {
            std::fprintf (stderr, "bad, unknown or repeated option\n");
            std::fclose (f);
            return 2;
        }
        const analysis::PeakMix mix = sMix == "avr" ? analysis::PeakMix::Average
                                    : sMix == "L"   ? analysis::PeakMix::Left
                                    : sMix == "R"   ? analysis::PeakMix::Right
                                    : sMix == "max" ? analysis::PeakMix::Max
                                    : (analysis::PeakMix) -1;
        if (sTo.empty()) to = frames;

        if (mode == "needle")
        {
            // The needle is over a stretch of planes already in memory; read the file whole.
            std::vector<float> inter ((std::size_t) (frames * (std::uint64_t) nc));
            const bool ok = std::fread (inter.data(), sizeof (float), inter.size(), f) == inter.size();
            std::fclose (f);
            std::vector<float> L ((std::size_t) frames), R ((std::size_t) frames);
            for (std::uint64_t i = 0; i < frames; ++i)
            {
                L[(std::size_t) i] = inter[(std::size_t) (i * (std::uint64_t) nc)];
                R[(std::size_t) i] = inter[(std::size_t) (i * (std::uint64_t) nc + (nc > 1 ? 1 : 0))];
            }
            analysis::StereoColumns::Needle nd;
            if (! ok || ! analysis::StereoColumns::needle (L.data(), R.data(), frames, from, to, nd))
            {
                std::fprintf (stderr, "needle refused: a stretch [from, to) inside the file is required\n");
                return 2;
            }
            std::printf ("# fcore needle v1 ch=%d frames=%llu from=%llu to=%llu\n", nc,
                         (unsigned long long) frames, (unsigned long long) from, (unsigned long long) to);
            std::printf ("corr %016llx\n",  (unsigned long long) bits (nd.correlation));
            std::printf ("width %016llx\n", (unsigned long long) bits (nd.width));
            std::printf ("rms %016llx\n",   (unsigned long long) bits (nd.rms));
            return 0;
        }

        fcore::ShapeProbe shapes;
        if (buckets > 0x7FFFFFFFu || columns > 0x7FFFFFFFu
         || ! shapes.prepare (fs, nc, frames, (int) buckets, mix, (int) columns))
        {
            std::fprintf (stderr, "shapes.prepare refused (buckets/columns 1..%d, mix avr|L|R|max)\n",
                          analysis::WaveformPeaks::kMaxBuckets);
            std::fclose (f);
            return 2;
        }
        bool ok = true;
        streamPlanar (f, nc, [&] (const float* const* p, int n) { ok = ok && shapes.process (p, nc, n); });
        std::fclose (f);
        if (! ok || ! shapes.complete())
        {
            std::fprintf (stderr, "the file did not deliver the frames it was sized for\n");
            return 2;
        }

        if (mode == "waveform")
        {
            const auto& w = shapes.peaks();
            std::printf ("# fcore waveform v1 sr=%016llx ch=%d frames=%llu buckets=%d mix=%s decim=%d emitted=%d\n",
                         (unsigned long long) bits (fs), nc, (unsigned long long) frames, w.buckets(), sMix.c_str(),
                         w.decimation(), w.bucketsEmitted());
            for (int i = 0; i < w.buckets(); ++i)
                std::printf ("%016llx %08lx\n", (unsigned long long) bits (w.peaks()[(std::size_t) i]),
                             (unsigned long) bits32 (w.peakAsFloat32 (i)));
        }
        else
        {
            const auto& s = shapes.stereo();
            std::printf ("# fcore stereo v1 ch=%d frames=%llu cols=%d mono=%d\n", nc, (unsigned long long) frames,
                         s.columns(), s.isMono() ? 1 : 0);
            std::printf ("maxrms %016llx\n", (unsigned long long) bits (s.maxRms()));
            for (int i = 0; i < s.columns(); ++i)
                std::printf ("%08lx %08lx %08lx\n", (unsigned long) bits32 (s.width()[(std::size_t) i]),
                             (unsigned long) bits32 (s.correlation()[(std::size_t) i]),
                             (unsigned long) bits32 (s.rms()[(std::size_t) i]));
        }
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
