// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_measure — a tiny streaming measurement CLI for validating felitronics-core's analysis against a
// reference (ffmpeg), and the NATIVE SIDE of the P0 wasm parity check. Reads interleaved 32-bit-float
// little-endian PCM (as ffmpeg emits with `-f f32le`).
//
// The measurement itself lives in fcore_probe.h (and, for `clips`, fcore_clips.h), shared verbatim with the
// wasm shim — read those headers for the build contract (-ffp-contract=off here, no -mrelaxed-simd there) and
// for why the true-peak filter config is part of that contract rather than an implementation detail.
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
//   clips       → the clipped runs (analysis::ClipDetector): how many were found, whether the list is whole,
//                 the sample peak of each channel, and every stored run — start, length, level, channel,
//                 polarity, evidence. Levels and peaks as bit patterns; the format lives in tools/fcore_clips_format.h
//                 and the wasm module prints it too, so a diff IS the parity test. [--max-runs N] [--chunk N]
//   stream      → the streaming surface (fcore::StreamProbe, tools/fcore_stream.h): the file fed in pieces of
//                 --chunk N frames (default 4096) and, after EACH piece, the deterministic meter's momentary,
//                 short-term and integrated LUFS as bit patterns, the samples consumed, the blocks dropped and
//                 the runs decided so far; then finish() and the runs only it decides. The wasm module prints
//                 the same bytes through its fc_stream_* handles (tools/wasm/stream-parity.mjs). [--chunk N]
//   report      → the whole-programme report (analysis::ProgrammeReport): DC, silence, tail, infra-low,
//                 stereo, PLR / LRA / short-term percentiles. Every scalar as a raw bit pattern with its
//                 validity and reason; every count as a decimal integer. Printed through the report's own
//                 field visitor, so the struct and this output cannot drift apart.
//   hum         → mains hum (analysis::HumDetector): per channel the validity reason, the mains nominal, the
//                 line's interpolated position / level / prominence, the comb, and the quiet stretches in
//                 sample coordinates — every float as a raw bit pattern, like `blocks`, so a future wasm
//                 comparison catches a flipped bit that %.17g would round away. `valid 0` is never "clean":
//                 read the reason. [--quiet-db X] [--order N]
//   lowend      → the vinyl low end (analysis::LowEnd): the integral Mid/Side energies of the LR4 low and high
//                 bands and of the unfiltered programme, every stored 10 ms block, the side-fraction histogram,
//                 the three extremum coordinates, and the full semitone band table with the dominant note —
//                 as raw IEEE-754 bit patterns, so `diff` between two toolchains IS the parity test. Every
//                 field the report publishes is here EXCEPT the law-8a trace, which exists only for the suite.
//                 An invalid note prints as `note INVALID reason N`, never as a note name. NOT `correlation`,
//                 which is a whole-file phase number with no band split.
//
//   forensics   → what the file WAS (analysis::SourceForensics): the spectral wall per channel and for the
// and forensics modes size the file before reading it (the first three because every boundary depends on the length,
// forensics so that a short read cannot come out as a successful measurement of a shorter programme), so they need a
//
// Usage: fcore_measure <mode> <sampleRate> <channels> <raw.f32le> [--precise] [mode options]
//
// The scalar modes print %.2f by default (tools/validate_ffmpeg.sh compares against ffmpeg's own two
// decimals); `--precise` switches them to %.17g plus the exact %a form. `blocks`, `waveform`, `stereo`,
// `needle` and `clips` are always exact, and the node side (tools/wasm/shapes-parity.mjs for the shapes,
// clips-parity.mjs for the runs) prints the same bytes.
//
// DECODE TO FLOAT, NEVER TO s16. This tool reads f32le: `ffmpeg -i x -f f32le out.f32`. The waveform, stereo, needle
// and clips modes size the file before reading it (the first three because every boundary depends on the length, clips
// because a short read must be a refusal rather than a clean report of a truncated file), so they need a seekable file,
// not a pipe. An integer decode
// (`-f s16le`, as the old sidecar generator did) clamps a lossy file's samples above 0 dBFS and quantises the
// rest, and the peaks of that are not the peaks of the file.

#include "fcore_clips_format.h"
#include "fcore_probe.h"
#include "fcore_stream.h"

#include <felitronics/analysis/ProgrammeReport.h>
#include <felitronics/analysis/SourceForensics.h>
#include <felitronics/analysis/HumDetector.h>
#include <felitronics/analysis/LowEnd.h>
#include <felitronics/analysis/BandBursts.h>
#include <felitronics/analysis/BandCrest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;

namespace
{
    constexpr int kChunk = fcore::Probe::kChunk;

    // fcore::ClipProbe spells its own step as a literal rather than deriving it from fcore::Probe: it must not
    // include fcore_probe.h, which would pull the loudness meter and the oversampler into the wasm module's
    // clips path for a constant. This is the one translation unit that sees both, so this is where the two are
    // tied — the `clips` mode's --chunk 8191/8192/8193 rows and the comments in fcore_clips.h about "the
    // native reader's 8192-frame reads" all go quietly false if they ever drift apart.
    static_assert (fcore::ClipProbe::kChunk == fcore::Probe::kChunk,
                   "the clips adapter's step and the file reader's step must be the same number");

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

    // A FINITE double, sign and zero allowed — `parseRate` insists on positive, which is right for a rate and
    // wrong for a threshold in dB, where 0 and a negative value are both ordinary. Same strictness otherwise:
    // the whole string must be the number, so "6dB" is refused rather than read as 6.
    bool parseFinite (const char* s, double& out)
    {
        char* end = nullptr;
        const double v = std::strtod (s, &end);
        if (end == s || *end != '\0' || ! std::isfinite (v)) return false;
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
            "usage: %s <lufs|truepeak|correlation|blocks|waveform|stereo|needle|clips|stream|report|hum|lowend|bursts|crest|forensics> <sampleRate> <channels> <raw.f32le>\n"
            "          [--precise] [--buckets N] [--mix avr|L|R|max] [--columns N] [--from A --to B]\n"
            "          [--max-runs N] [--chunk N]\n"
            "          [--quiet-db X] [--order N]\n",
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

    if (mode == "hum")
    {
        // THE FILE IS SIZED BEFORE IT IS READ, as clips/forensics/lowend already do it. streamPlanar()
        // floors a trailing partial frame and returns success regardless, so a file holding three floats
        // and declared stereo was measured to its break and reported as a whole programme — an instrument
        // certifying audio it never saw, which is the exact failure P71 says it closed. Three of the six
        // modes had inherited the lax path; this is the third of them.
        std::uint64_t declaredFrames = 0;
        if (! fileFrames (f, nc, declaredFrames))
        {
            std::fprintf (stderr, "cannot size the file, or it is not a whole number of %d-channel float32 frames\n", nc);
            std::fclose (f);
            return 2;
        }
        // STRICT ARGUMENTS, as clips/forensics/lowend already are. The shared parse at the top of main()
        // uses atoi/atof, which read "1.5" as 1 channel and "48000Hz" as 48000 — and then measure, happily,
        // the wrong thing. These three modes are new in this release, so tightening them breaks nothing;
        // the wasm harness refuses the same strings, and a refusal that does not match on both roads is a
        // parity break that a diff of two successful runs would never show.
        {
            double sRate = 0.0; std::uint64_t sWidth = 0;
            if (! parseRate (argv[2], sRate) || ! parseCount (argv[3], sWidth)
                || sWidth < 1 || sWidth > (std::uint64_t) core::kMaxChannels)
            {
                std::fprintf (stderr, "bad sampleRate/channels\n");
                std::fclose (f);
                return 2;
            }
        }
        // The whole file through analysis::HumDetector, streamed in kChunk steps — the answer is identical at
        // any slicing (law 8a), so the chunking is a convenience here and not part of the measurement.
        analysis::HumDetectorParams hp;
        for (int i = 5; i < argc; ++i)
        {
            if (std::strcmp (argv[i], "--quiet-db") == 0 && i + 1 < argc) hp.quietThresholdDb = std::atof (argv[++i]);
            else if (std::strcmp (argv[i], "--order") == 0 && i + 1 < argc) hp.fftOrder = std::atoi (argv[++i]);
        }
        analysis::HumDetector hd;
        hd.setParams (hp);
        if (! hd.prepare (fs, kChunk, nc))
        {
            std::fprintf (stderr, "hum: prepare refused these arguments (rate %g, %d channels)\n", fs, nc);
            std::fclose (f);
            return 2;
        }
        streamPlanar (f, nc, [&] (const float* const* p, int n) { (void) hd.process (p, nc, n); });
        std::fclose (f);
        hd.finish();
        if ((std::uint64_t) hd.samplesProcessed() != declaredFrames)
        {
            std::fprintf (stderr, "read %lld of %llu frames — refusing to report a partial measurement\n",
                          (long long) hd.samplesProcessed(), (unsigned long long) declaredFrames);
            return 2;
        }
        std::printf ("# fcore hum v1 sr=%016llx ch=%d order=%d n=%lld hop=%lld bin=%016llx\n",
                     (unsigned long long) bits (fs), nc, hd.geometry().order,
                     (long long) hd.windowSamples(), (long long) hd.hopSamples(),
                     (unsigned long long) bits (hd.binHz()));
        for (int c = 0; c < nc; ++c)
        {
            const analysis::HumReport r = hd.report (c);
            std::printf ("ch %d valid %d reason %d mains %d base %d fobs %d fderived %d\n",
                         c, r.valid ? 1 : 0, (int) r.reason, (int) r.mains, r.baseHarmonic,
                         r.fundamentalObserved ? 1 : 0, r.fundamentalDerived ? 1 : 0);
            std::printf ("ch %d f0 %016llx hz %016llx tone %016llx peakbin %016llx floor %016llx prom %016llx\n",
                         c, (unsigned long long) bits (r.fundamentalHz), (unsigned long long) bits (r.line.hz),
                         (unsigned long long) bits (r.line.tonePower), (unsigned long long) bits (r.line.peakBinPower),
                         (unsigned long long) bits (r.line.floorPower), (unsigned long long) bits (r.line.prominenceDb));
            std::printf ("ch %d frames %lld finite %lld holed %lld quiet %lld stretches %lld stored %lld complete %d tail %lld\n",
                         c, (long long) r.frames, (long long) r.finiteFrames, (long long) r.holedFrames,
                         (long long) r.quietFrames, (long long) r.quietStretches, (long long) r.storedStretches,
                         r.stretchesComplete ? 1 : 0, (long long) r.tailUncoveredSamples);
            for (int cand = 0; cand < analysis::HumDetector::kCandidates; ++cand)
            {
                const analysis::HumCandidate k = hd.candidate (c, cand);
                std::printf ("ch %d cand %016llx found %d base %d f0 %016llx sobs %lld soff %lld fobs %lld sspread %016llx fspread %016llx intra %016llx stat %d pass %d harm %d low %d\n",
                             c, (unsigned long long) bits (k.nominalHz), k.baseFound ? 1 : 0, k.baseHarmonic,
                             (unsigned long long) bits (k.fundamentalHz),
                             (long long) k.stretchObservations, (long long) k.stretchOffTolerance,
                             (long long) k.frameObservations,
                             (unsigned long long) bits (k.stretchSpreadHz), (unsigned long long) bits (k.frameSpreadHz),
                             (unsigned long long) bits (k.maxIntraStretchSpreadHz),
                             k.stationary ? 1 : 0, k.passed ? 1 : 0, k.harmonicsObserved, k.lowestHarmonicObserved);
                std::printf ("ch %d cand %d window %d hz %016llx prom %016llx\n", c, cand,
                             k.windowPeak.found ? 1 : 0, (unsigned long long) bits (k.windowPeak.hz),
                             (unsigned long long) bits (k.windowPeak.prominenceDb));
                for (int h = 1; h <= hp.maxHarmonic; ++h)
                {
                    const analysis::HumHarmonic hh = hd.harmonic (c, cand, h);
                    std::printf ("ch %d cand %d h %d inband %d acc %d hz %016llx tone %016llx prom %016llx\n",
                                 c, cand, h, hh.inBand ? 1 : 0, hh.peak.accepted ? 1 : 0,
                                 (unsigned long long) bits (hh.peak.hz), (unsigned long long) bits (hh.peak.tonePower),
                                 (unsigned long long) bits (hh.peak.prominenceDb));
                }
            }
            for (std::int64_t i = 0; i < hd.storedStretchCount (c); ++i)
            {
                const analysis::HumStretch st = hd.stretch (c, i);
                std::printf ("ch %d stretch %lld %lld %lld frames %lld\n", c, (long long) st.index,
                             (long long) st.startSample, (long long) st.endSample, (long long) st.frames);
            }
        }
        return 0;
    }

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

    if (mode == "clips")
    {
        // THE FILE IS SIZED FIRST, AND A SHORT READ IS A REFUSAL. streamPlanar() floors a trailing partial frame,
        // ignores what the sink answered and never looks at ferror(), so a truncated or unreadable file would
        // otherwise be measured to its break and reported as a clean whole — an instrument certifying audio it
        // never saw. fileFrames() refuses a size that is not a whole number of frames, the sink's verdict is kept,
        // and the frames the probe actually consumed are compared with the frames the file holds.
        //
        // AN EMPTY FILE IS REFUSED, as it is in waveform|stereo|needle and in the wasm module's own input guard
        // (fc_probe.cpp planarSpan, which rejects frames == 0). A zero-length programme has a perfectly good
        // clip report — no runs, zero peaks — and printing it was the other candidate here; refusing wins
        // because the most likely way to arrive at a zero-length file is truncation, and an instrument that
        // answers "clean" to a file that lost its contents is the failure this mode is built to avoid. The
        // zero-length stream itself is still covered, in felitronics_clips_exposure_tests, where it is a
        // measurement and not a file.
        double rate = 0.0; std::uint64_t width = 0;
        if (! parseRate (argv[2], rate) || ! parseCount (argv[3], width) || width < 1
            || width > (std::uint64_t) core::kMaxChannels)
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
        // Parsed with the shape modes' strictness and NOT with their parser: an unknown, misspelt or repeated
        // option is refused rather than ignored or silently resolved, because the node side of the parity check
        // could not be relied on to resolve it the same way.
        //
        // --chunk is the law-8a handle, and it is a real option rather than a test hook: the report must be the
        // same bytes however the stream is cut into process() calls, and without a way to ask for a different
        // cut from outside, that claim can only ever be checked from inside a test binary. 0 means "one call per
        // read block", which is what the tool does by itself; any N >= 1 sub-slices those blocks, so the reachable
        // calls run from a single sample up to the reader's own 8192-frame step. Calls LARGER than a read block
        // are reached the only place they can be — felitronics_clips_exposure_tests, which drives fcore::ClipProbe
        // (the same adapter both roads use) directly.
        std::uint64_t maxRuns = (std::uint64_t) analysis::ClipDetectorParams {}.maxRuns, chunk = 0;
        bool seen[2] {};
        for (int i = 5; i < argc; ++i)
        {
            if (std::strcmp (argv[i], "--precise") == 0) continue;
            const int k = std::strcmp (argv[i], "--max-runs") == 0 ? 0 : std::strcmp (argv[i], "--chunk") == 0 ? 1 : 2;
            if (k == 2 || seen[k] || i + 1 >= argc || ! parseCount (argv[i + 1], k == 0 ? maxRuns : chunk))
            {
                std::fprintf (stderr, "bad, unknown or repeated option\n");
                std::fclose (f);
                return 2;
            }
            seen[k] = true;
            ++i;
        }
        if (chunk > 0x7FFFFFFFu) { std::fprintf (stderr, "--chunk out of range\n"); std::fclose (f); return 2; }

        fcore::ClipProbe probe;
        if (maxRuns > (std::uint64_t) fcore::ClipProbe::kMaxRuns
            || frames > (std::uint64_t) std::numeric_limits<std::int64_t>::max()
            || ! probe.prepare (rate, nc, (std::int64_t) maxRuns, (std::int64_t) frames))
        {
            // The width is already checked twice above, so it cannot be the cause here and is not offered as one.
            std::fprintf (stderr, "clips.prepare refused (sampleRate %g..%g, --max-runs 0..%lld)\n",
                          analysis::ClipDetector::kMinSampleRate, analysis::ClipDetector::kMaxSampleRate,
                          (long long) fcore::ClipProbe::kMaxRuns);
            std::fclose (f);
            return 2;
        }
        bool ok = true;
        const long long step = chunk == 0 ? 0 : (long long) chunk;
        streamPlanar (f, nc, [&] (const float* const* p, int n)
        {
            if (step == 0) { ok = ok && probe.process (p, nc, n); return; }
            const float* view[core::kMaxChannels] {};
            for (long long off = 0; off < n && ok; off += step)
            {
                const long long m = std::min<long long> (step, (long long) n - off);
                for (int c = 0; c < nc; ++c) view[(std::size_t) c] = p[c] + off;
                ok = ok && probe.process (view, nc, m);
            }
        });
        std::fclose (f);
        // probe.finish() is the check now: it refuses unless the frames the file was sized for, the frames the
        // reader handed over and the samples the detector consumed are all the same number.
        fcore::ClipsReport rep;
        if (! ok || ! probe.finish() || ! readClips (probe, rep))
        {
            std::fprintf (stderr, "the file did not deliver the %llu frames it was sized for\n",
                          (unsigned long long) frames);
            return 2;
        }
        const std::string text = fcore::formatClips (rep);
        std::fwrite (text.data(), 1, text.size(), stdout);
        return 0;
    }

    if (mode == "stream")
    {
        // Sized first and read WHOLE, then cut into pieces of exactly --chunk frames with a shorter last one —
        // the schedule stream-parity.mjs cuts on its side, so the two roads read between the same samples. (The
        // 8192-frame reader of the other modes would restart the schedule at every read block.) An empty file and
        // a partial frame are refused, as `clips` refuses them. Nothing is printed until the whole run succeeded:
        // a refusal is an empty stdout, never a plausible prefix.
        double rate = 0.0; std::uint64_t width = 0;
        if (! parseRate (argv[2], rate) || ! parseCount (argv[3], width) || width < 1
            || width > (std::uint64_t) core::kMaxChannels)
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
        std::uint64_t chunk = 4096;
        bool seen = false;
        for (int i = 5; i < argc; ++i)
        {
            if (std::strcmp (argv[i], "--precise") == 0) continue;
            if (std::strcmp (argv[i], "--chunk") != 0 || seen || i + 1 >= argc || ! parseCount (argv[i + 1], chunk))
            {
                std::fprintf (stderr, "bad, unknown or repeated option\n");
                std::fclose (f);
                return 2;
            }
            seen = true;
            ++i;
        }
        if (chunk < 1 || chunk > 0x7FFFFFFFu) { std::fprintf (stderr, "--chunk out of range\n"); std::fclose (f); return 2; }

        std::vector<float> inter ((std::size_t) (frames * (std::uint64_t) nc));
        const bool whole = std::fread (inter.data(), sizeof (float), inter.size(), f) == inter.size();
        std::fclose (f);
        if (! whole) { std::fprintf (stderr, "the file did not deliver the %llu frames it was sized for\n", (unsigned long long) frames); return 2; }
        std::vector<std::vector<float>> planes ((std::size_t) nc, std::vector<float> ((std::size_t) frames));
        for (std::uint64_t i = 0; i < frames; ++i)
            for (int c = 0; c < nc; ++c)
                planes[(std::size_t) c][(std::size_t) i] = inter[(std::size_t) (i * (std::uint64_t) nc + (std::uint64_t) c)];

        fcore::StreamProbe sp;
        if (! sp.prepare (rate, nc))
        {
            std::fprintf (stderr, "stream.prepare refused (sampleRate %g..%g)\n",
                          analysis::ClipDetector::kMinSampleRate, analysis::ClipDetector::kMaxSampleRate);
            return 2;
        }
        std::string out;
        char line[256];
        std::snprintf (line, sizeof line, "# fcore stream v1 sr=%016llx ch=%d chunk=%llu\n",
                       (unsigned long long) bits (rate), nc, (unsigned long long) chunk);
        out += line;
        std::int64_t read = 0;
        auto drain = [&]
        {
            const auto& d = sp.detector();
            for (; read < d.storedRunCount(); ++read)
            {
                const analysis::ClipRun u = d.run (read);
                std::snprintf (line, sizeof line, "run %lld %lld %016llx %d %d %d\n", (long long) u.start,
                               (long long) u.length, (unsigned long long) bits (u.level), u.channel, u.sign,
                               (int) u.evidence);
                out += line;
            }
        };
        const float* view[core::kMaxChannels] {};
        for (std::uint64_t at = 0; at < frames; at += chunk)
        {
            const std::uint64_t m = std::min (chunk, frames - at);
            for (int c = 0; c < nc; ++c) view[c] = planes[(std::size_t) c].data() + at;
            if (! sp.process (view, (int) m)) { std::fprintf (stderr, "stream.process refused at frame %llu\n", (unsigned long long) at); return 2; }
            const auto& lm = sp.meter();
            std::snprintf (line, sizeof line, "at %lld m %016llx s %016llx i %016llx dropped %d runs %lld\n",
                           (long long) sp.samples(), (unsigned long long) bits (lm.momentaryLufs()),
                           (unsigned long long) bits (lm.shortTermLufs()), (unsigned long long) bits (lm.integratedLufs()),
                           lm.droppedBlocks(), (long long) sp.detector().runCount());
            out += line;
            drain();
        }
        if (! sp.finish()) { std::fprintf (stderr, "stream.finish refused\n"); return 2; }
        std::snprintf (line, sizeof line, "finish runs %lld\n", (long long) sp.detector().runCount());
        out += line;
        drain();
        std::fwrite (out.data(), 1, out.size(), stdout);
        return 0;
    }

    if (mode == "report")
    {
        // THE FILE IS SIZED BEFORE IT IS READ, as clips/forensics/lowend already do it. streamPlanar()
        // floors a trailing partial frame and returns success regardless, so a file holding three floats
        // and declared stereo was measured to its break and reported as a whole programme — an instrument
        // certifying audio it never saw, which is the exact failure P71 says it closed. Three of the six
        // modes had inherited the lax path; this is the third of them.
        std::uint64_t declaredFrames = 0;
        if (! fileFrames (f, nc, declaredFrames))
        {
            std::fprintf (stderr, "cannot size the file, or it is not a whole number of %d-channel float32 frames\n", nc);
            std::fclose (f);
            return 2;
        }
        // STRICT ARGUMENTS, as clips/forensics/lowend already are. The shared parse at the top of main()
        // uses atoi/atof, which read "1.5" as 1 channel and "48000Hz" as 48000 — and then measure, happily,
        // the wrong thing. These three modes are new in this release, so tightening them breaks nothing;
        // the wasm harness refuses the same strings, and a refusal that does not match on both roads is a
        // parity break that a diff of two successful runs would never show.
        {
            double sRate = 0.0; std::uint64_t sWidth = 0;
            if (! parseRate (argv[2], sRate) || ! parseCount (argv[3], sWidth)
                || sWidth < 1 || sWidth > (std::uint64_t) core::kMaxChannels)
            {
                std::fprintf (stderr, "bad sampleRate/channels\n");
                std::fclose (f);
                return 2;
            }
        }
        // THE WHOLE-PROGRAMME REPORT. Every floating-point number goes out as a raw IEEE-754 bit pattern,
        // exactly as `blocks` does and for the same reason: the other side of this comparison is
        // JavaScript, which has no hex-float printing and whose decimal formatting is not C's, so a
        // 16-hex-digit pattern is the one representation both sides produce identically and `diff` IS the
        // parity test. `%.17g` would not do either — it round-trips, but a flipped low bit can print the
        // same decimal on two libcs.
        //
        // Printed through the report's OWN field visitor rather than a list written out here, so a field
        // added to the struct appears in this output without this block being touched — the same
        // enumeration the law-8a gate compares through.
        //
        // maxBlock is kChunk because that is what streamPlanar hands over; it sizes the scratch and
        // nothing else, so the numbers do not depend on it.
        analysis::ProgrammeReport pr;
        if (! pr.prepare (fs, kChunk, nc))
        {
            std::fprintf (stderr, "report.prepare refused (rate, channels or a parameter out of range)\n");
            std::fclose (f);
            return 2;
        }
        bool accepted = true;
        streamPlanar (f, nc, [&] (const float* const* p, int n) { accepted = accepted && pr.process (p, nc, n); });
        std::fclose (f);
        if (! accepted)
        {
            std::fprintf (stderr, "the report refused a call\n");
            return 2;
        }
        pr.finish();
        if ((std::uint64_t) pr.samplesProcessed() != declaredFrames)
        {
            std::fprintf (stderr, "read %lld of %llu frames — refusing to report a partial measurement\n",
                          (long long) pr.samplesProcessed(), (unsigned long long) declaredFrames);
            return 2;
        }
        const auto& R = pr.report();
        std::printf ("# fcore report v1 sr=%016llx ch=%d samples=%lld\n",
                     (unsigned long long) bits (fs), nc, (long long) R.totalSamples);
        R.visitCounts ([] (const char* name, int ch, std::int64_t v)
                       { std::printf ("C %s %d %lld\n", name, ch, (long long) v); });
        R.visitValues ([] (const char* name, int ch, const analysis::ProgrammeValue& v)
                       { std::printf ("V %s %d %d %d %016llx\n", name, ch, v.valid ? 1 : 0,
                                      (int) v.reason, (unsigned long long) bits (v.value)); });
        return 0;
    }

    if (mode == "lowend")
    {
        // The vinyl low end. Streamed in kChunk steps, which is also the point: the report is bit-identical
        // under ANY slicing (law 8a), so the chunk size is not part of the measurement.
        // STRICT on its arguments, like the shape modes and unlike the older scalar ones: atoi("4294967297")
        // narrows to 1 channel and atof("8000Hz") reads 8000, and either would measure something silently.
        double rate = 0.0; std::uint64_t width = 0;
        if (! parseRate (argv[2], rate) || ! parseCount (argv[3], width)
            || width < 1 || width > (std::uint64_t) core::kMaxChannels)
        {
            std::fprintf (stderr, "bad sampleRate/channels\n");
            std::fclose (f);
            return 2;
        }
        // …and strict on its INPUT: a file that is not a whole number of frames, or that cannot be read
        // to the end, must not print a report and exit zero. (A directory opens successfully on macOS and
        // then fails every read, which used to come out as an empty report and a success.)
        std::uint64_t frames = 0;
        if (! fileFrames (f, nc, frames) || frames == 0)
        {
            std::fprintf (stderr, "cannot size the file, it is not a whole number of %d-channel float32 frames, or it is empty\n", nc);
            std::fclose (f);
            return 2;
        }
        analysis::LowEndParams lp;
        analysis::LowEnd le;
        le.setParams (lp);
        if (! le.prepare (rate, kChunk, nc))
        {
            std::fprintf (stderr, "LowEnd refused this geometry (rate, channels or note range)\n");
            std::fclose (f);
            return 2;
        }
        bool okAll = true;
        const bool read = streamPlanar (f, nc, [&] (const float* const* pp, int n) { okAll = le.process (pp, nc, n) && okAll; });
        std::fclose (f);
        if (! read || ! okAll || ! le.finish()) { std::fprintf (stderr, "LowEnd refused a chunk, or the file could not be read\n"); return 2; }
        if ((std::uint64_t) le.samplesProcessed() != frames)
        {
            std::fprintf (stderr, "read %lld of %llu frames — refusing to report a partial measurement\n",
                          (long long) le.samplesProcessed(), (unsigned long long) frames);
            return 2;
        }

        // Raw bit patterns, not %g: the other side of this comparison is JavaScript, whose decimal
        // formatting is not C's, so a 16-hex-digit pattern is the one representation both sides produce
        // identically. A decimal header would break a whole-file diff while every measured bit matched.
        std::printf ("# fcore lowend v2 sr=%016llx ch=%d xover=%016llx order=%d hop=%lld block=%lld bands=%d chunk=%d\n",
                     (unsigned long long) bits (fs), nc, (unsigned long long) bits (le.crossoverHz()),
                     lp.fftOrder, (long long) le.hopSamples(), (long long) le.blockSamples(), le.bandCount(), kChunk);
        std::printf ("reason %d %d\n", (int) le.widthReason(), (int) le.noteReason());
        std::printf ("samples %lld finite %lld holes %lld nonfinite %lld absent %lld overflow %lld\n",
                     (long long) le.samplesProcessed(), (long long) le.finiteSamples(), (long long) le.holeSamples(),
                     (long long) le.nonFiniteSamples(), (long long) le.absentSamples(),
                     (long long) le.filterNonFiniteSamples());
        const double energies[7] = { le.lowMidEnergy(), le.lowSideEnergy(), le.highMidEnergy(), le.highSideEnergy(),
                                     le.rawMidEnergy(), le.rawSideEnergy(), le.lowSideFraction() };
        static const char* const enames[7] = { "lowmid", "lowside", "highmid", "highside", "rawmid", "rawside", "lowfrac" };
        for (int i = 0; i < 7; ++i)
            std::printf ("%s %016llx\n", enames[i], (unsigned long long) bits (energies[i]));
        std::printf ("highfrac %016llx rawfrac %016llx\n",
                     (unsigned long long) bits (le.highSideFraction()), (unsigned long long) bits (le.rawSideFraction()));
        std::printf ("blocks %lld stored %lld complete %d histsamples %lld\n",
                     (long long) le.blockCount(), (long long) le.storedBlockCount(), le.blocksComplete() ? 1 : 0,
                     (long long) le.histogramSamples());
        for (int i = 0; i < analysis::LowEnd::kHistogramBins; ++i)
            std::printf ("h%02d %lld\n", i, (long long) le.histogram (i));
        std::printf ("worst %lld %016llx %016llx\n", (long long) le.worstFractionBlock(),
                     (unsigned long long) bits (le.worstFraction()), (unsigned long long) bits (le.worstFractionEnergy()));
        std::printf ("peakenergy %lld %016llx %016llx\n", (long long) le.peakEnergyBlock(),
                     (unsigned long long) bits (le.peakBlockEnergy()), (unsigned long long) bits (le.peakEnergyBlockFraction()));
        std::printf ("peakside %lld %016llx amp %016llx at %lld\n", (long long) le.peakSideEnergyBlock(),
                     (unsigned long long) bits (le.peakBlockSideEnergy()),
                     (unsigned long long) bits (le.peakLowSideAmplitude()), (long long) le.peakLowSideAmplitudeAt());
        std::printf ("frames used %lld holed %lld tail %lld window %lld underresolved %d\n",
                     (long long) le.usedFrames(), (long long) le.holedFrames(), (long long) le.tailUncoveredSamples(),
                     (long long) le.windowSamples(), le.underResolvedBands());
        // the 10 ms SERIES, every stored block: without it a diff cannot compare the quantity the
        // instrument publishes per block, which is where the wide-bass answer actually lives
        std::printf ("series index samples finite holes midEnergy sideEnergy\n");
        for (std::int64_t i = 0; i < le.storedBlockCount(); ++i)
        {
            const analysis::LowEndBlock r = le.block (i);
            std::printf ("s %lld %lld %lld %lld %016llx %016llx\n", (long long) r.index, (long long) r.samples,
                         (long long) r.finiteSamples, (long long) r.holes,
                         (unsigned long long) bits (r.midEnergy), (unsigned long long) bits (r.sideEnergy));
        }
        std::printf ("band midi centreHz widthHz binsPerBand midEnergy sideEnergy energy density centroidHz centsOffset dutyCount levelWhenOnDb\n");
        for (int b = 0; b < le.bandCount(); ++b)
        {
            const analysis::LowEndBand r = le.band (b);
            std::printf ("b %d %d %016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx %lld %016llx\n", b, r.midi,
                         (unsigned long long) bits (r.centreHz), (unsigned long long) bits (r.widthHz),
                         (unsigned long long) bits (r.binsPerBand), (unsigned long long) bits (r.midEnergy),
                         (unsigned long long) bits (r.sideEnergy), (unsigned long long) bits (r.energy),
                         (unsigned long long) bits (r.density), (unsigned long long) bits (r.centroidHz),
                         (unsigned long long) bits (r.centsOffset),
                         (long long) le.dutyCount (b), (unsigned long long) bits (le.levelWhenOnDb (b)));
        }
        std::printf ("peak %d %d density %d second %d\n", le.peakBand(), le.peakMidi(),
                     le.peakDensityBand(), le.secondBand());
        // Only when there IS a note. peakMidi() is canonically 0 for an invalid report, and feeding that
        // through the naming functions printed "note C-1" — an invalid answer wearing a real note's name,
        // which is precisely the number-that-reads-as-a-finding this instrument exists not to print.
        if (le.noteValid())
            std::printf ("note %s%d nominal %016llx centroid %016llx cents %016llx sidefrac %016llx\n",
                         analysis::LowEnd::pitchClassName (le.peakMidi()), analysis::LowEnd::noteOctave (le.peakMidi()),
                         (unsigned long long) bits (le.peakNoteHz()), (unsigned long long) bits (le.peakCentroidHz()),
                         (unsigned long long) bits (le.peakCentsOffset()), (unsigned long long) bits (le.peakBandSideFraction()));
        else
            std::printf ("note INVALID reason %d\n", (int) le.noteReason());
        // the dominance ratio is NOT printed as one number: its denominator is exactly zero for a tone in
        // digital silence. The three numbers it is made of are printed instead.
        std::printf ("frame %016llx rangeshare %016llx\n",
                     (unsigned long long) bits (le.frameEnergy()), (unsigned long long) bits (le.bandRangeShare()));
        std::printf ("background %016llx peakenergy %016llx peakwidth %016llx share %016llx total %016llx\n",
                     (unsigned long long) bits (le.backgroundDensity()), (unsigned long long) bits (le.peakBandEnergy()),
                     (unsigned long long) bits (le.peakBandWidthHz()), (unsigned long long) bits (le.peakShare()),
                     (unsigned long long) bits (le.totalBandEnergy()));
        // v2 — K9. The resolution boundary makes underResolvedBands actionable; the duty population is
        // published so a consumer never defines it twice; `skipped`/`asked` say what the histogram lost.
        std::printf ("resolved first %d above %016llx lobebins %d\n",
                     le.firstResolvedBand(), (unsigned long long) bits (le.resolvedAboveHz()),
                     analysis::LowEnd::lobeBins());
        std::printf ("duty frames %lld thresholddb %016llx skipped %lld asked %d\n",
                     (long long) le.dutyFrames(), (unsigned long long) bits (le.dutyThresholdDb()),
                     (long long) le.skippedBlocks(), lp.skipBlocks);
        return 0;
    }

    // K1 — BandCrest: the crest per 400 ms block and per band, and (with `--against`) the PAIRED loss between
    // this file as the SOURCE and that one as the master. Every float is a raw IEEE-754 bit pattern, as
    // `bursts` does it, so a later wasm comparison catches a flipped bit that decimal printing would round
    // away. The loss is computed HERE, in the core's own function, and not reassembled from the printed
    // per-run numbers: formed from linear cells it is one logarithm and invariant to the gain a master has
    // over its source, and a consumer re-deriving it from two dB columns would lose both properties.
    if (mode == "crest")
    {
        std::uint64_t declaredFrames = 0;
        if (! fileFrames (f, nc, declaredFrames))
        {
            std::fprintf (stderr, "cannot size the file, or it is not a whole number of %d-channel float32 frames\n", nc);
            std::fclose (f);
            return 2;
        }
        {
            double sRate = 0.0; std::uint64_t sWidth = 0;
            if (! parseRate (argv[2], sRate) || ! parseCount (argv[3], sWidth)
                || sWidth < 1 || sWidth > (std::uint64_t) core::kMaxChannels)
            {
                std::fprintf (stderr, "bad sampleRate/channels\n");
                std::fclose (f);
                return 2;
            }
        }
        analysis::BandCrestParams bp;
        std::string against;
        {
            struct Flag { const char* name; double* into; };
            const Flag flags[] = {
                { "--edge0",      &bp.bandEdgeHz[0] }, { "--edge1", &bp.bandEdgeHz[1] },
                { "--edge2",      &bp.bandEdgeHz[2] }, { "--hop-ms", &bp.hopMs },
                { "--floor-db",   &bp.programmeFloorDb }, { "--share-db", &bp.bandShareFloorDb },
            };
            for (int i = 5; i < argc; ++i)
            {
                const Flag* hit = nullptr;
                for (const Flag& fl : flags) if (std::strcmp (argv[i], fl.name) == 0) { hit = &fl; break; }
                if (hit != nullptr)
                {
                    if (i + 1 >= argc || ! parseFinite (argv[i + 1], *hit->into))
                    { std::fprintf (stderr, "crest: %s needs a finite number\n", hit->name); std::fclose (f); return 2; }
                    ++i; continue;
                }
                if (std::strcmp (argv[i], "--block-hops") == 0)
                {
                    std::uint64_t v = 0;
                    if (i + 1 >= argc || ! parseCount (argv[i + 1], v) || v < 1 || v > 64)
                    { std::fprintf (stderr, "crest: --block-hops needs 1..64\n"); std::fclose (f); return 2; }
                    bp.blockHops = (int) v; ++i; continue;
                }
                if (std::strcmp (argv[i], "--against") == 0)
                {
                    if (i + 1 >= argc) { std::fprintf (stderr, "crest: --against needs a file\n"); std::fclose (f); return 2; }
                    against = argv[i + 1]; ++i; continue;
                }
                // AN OPTION NOBODY KNOWS IS A REFUSAL — the rule the bursts mode learned the hard way.
                if (std::strncmp (argv[i], "--", 2) == 0 && std::strcmp (argv[i], "--precise") != 0)
                {
                    std::fprintf (stderr, "crest: unknown option %s (want --edge0 --edge1 --edge2 --hop-ms "
                                          "--block-hops --floor-db --share-db --against)\n", argv[i]);
                    std::fclose (f);
                    return 2;
                }
            }
        }

        auto runOne = [&] (std::FILE* fh, std::uint64_t frames, analysis::BandCrest& bc) -> bool
        {
            bc.setParams (bp);
            if (! bc.prepare (fs, nc, (long long) frames)) return false;
            bool okAll = true;
            streamPlanar (fh, nc, [&] (const float* const* p, int n) { okAll = bc.process (p, nc, n) && okAll; });
            bc.finish();
            return okAll;
        };

        analysis::BandCrest src;
        if (! runOne (f, declaredFrames, src))
        {
            std::fprintf (stderr, "crest: refused this configuration — edges %g/%g/%g Hz must rise and clear "
                                  "0.49 of the oversampled rate (%g Hz), hop %g ms, blockHops %d, at %g Hz\n",
                          bp.bandEdgeHz[0], bp.bandEdgeHz[1], bp.bandEdgeHz[2],
                          0.49 * fs * analysis::BandCrest::kFactor, bp.hopMs, bp.blockHops, fs);
            std::fclose (f);
            return 2;
        }
        std::fclose (f);
        if ((std::uint64_t) src.samplesProcessed() != declaredFrames)
        {
            std::fprintf (stderr, "read %lld of %llu frames — refusing to report a partial measurement\n",
                          (long long) src.samplesProcessed(), (unsigned long long) declaredFrames);
            return 2;
        }

        std::printf ("# fcore crest v1 sr=%016llx ch=%d hop=%d blockHops=%d e0=%016llx e1=%016llx e2=%016llx "
                     "floor=%016llx share=%016llx chunk=%d\n",
                     (unsigned long long) bits (fs), nc, src.hopSamples(), src.blockHops(),
                     (unsigned long long) bits (bp.bandEdgeHz[0]), (unsigned long long) bits (bp.bandEdgeHz[1]),
                     (unsigned long long) bits (bp.bandEdgeHz[2]), (unsigned long long) bits (bp.programmeFloorDb),
                     (unsigned long long) bits (bp.bandShareFloorDb), kChunk);
        std::printf ("samples %lld\n", (long long) src.samplesProcessed());
        std::printf ("hops %lld %lld %lld %lld %d\n", (long long) src.hopCount(), (long long) src.basePeakHops(),
                     (long long) src.blockCount(), (long long) src.nonFiniteSamples(),
                     (int) src.invalidReason());
        std::printf ("peak %016llx\n", (unsigned long long) bits (src.fullBandPeakLin()));
        std::printf ("active");
        for (int b = 0; b < analysis::BandCrest::kBands; ++b)
            std::printf (" %lld", (long long) src.activeBlocks (b));
        std::printf (" %016llx\n", (unsigned long long) bits (src.programmeMeanSquareDb()));
        for (long long j = 0; j < src.blockCount(); ++j)
        {
            std::printf ("b %lld", j);
            for (int b = 0; b < analysis::BandCrest::kBands; ++b)
                std::printf (" %016llx %016llx %d", (unsigned long long) bits (src.blockPeakLin (j, b)),
                             (unsigned long long) bits (src.blockMeanSq (j, b)), src.blockActive (j, b) ? 1 : 0);
            std::printf ("\n");
        }

        if (! against.empty())
        {
            std::FILE* g = std::fopen (against.c_str(), "rb");
            if (! g) { std::perror ("open --against"); return 2; }
            std::uint64_t gFrames = 0;
            if (! fileFrames (g, nc, gFrames))
            { std::fprintf (stderr, "--against: not a whole number of %d-channel float32 frames\n", nc); std::fclose (g); return 2; }
            analysis::BandCrest dst;
            if (! runOne (g, gFrames, dst)) { std::fprintf (stderr, "crest: --against refused\n"); std::fclose (g); return 2; }
            std::fclose (g);
            std::vector<double> scratch;
            for (int b = 0; b < analysis::BandCrest::kBands; ++b)
            {
                const auto L = analysis::bandCrestLoss (src, dst, b, scratch);
                std::printf ("loss %d %lld %lld %lld %016llx %016llx %016llx %016llx %016llx %016llx "
                             "%lld %lld %lld %016llx %016llx %lld %lld %d\n",
                             b, (long long) L.blocks, (long long) L.inActive, (long long) L.usable,
                             (unsigned long long) bits (L.p50Db), (unsigned long long) bits (L.p95Db),
                             (unsigned long long) bits (L.cvar95Db), (unsigned long long) bits (L.meanDb),
                             (unsigned long long) bits (L.maxDb), (unsigned long long) bits (L.p5Db),
                             (long long) L.over1Db, (long long) L.over3Db, (long long) L.over6Db,
                             (unsigned long long) bits (L.peakShiftDb), (unsigned long long) bits (L.levelShiftDb),
                             (long long) L.outSilent, (long long) L.lagBlocks, L.valid ? 1 : 0);
            }
        }
        return 0;
    }

    if (mode == "bursts")
    {
        // THE FILE IS SIZED BEFORE IT IS READ, as clips/forensics/lowend already do it. streamPlanar()
        // floors a trailing partial frame and returns success regardless, so a file holding three floats
        // and declared stereo was measured to its break and reported as a whole programme — an instrument
        // certifying audio it never saw, which is the exact failure P71 says it closed. Three of the six
        // modes had inherited the lax path; this is the third of them.
        std::uint64_t declaredFrames = 0;
        if (! fileFrames (f, nc, declaredFrames))
        {
            std::fprintf (stderr, "cannot size the file, or it is not a whole number of %d-channel float32 frames\n", nc);
            std::fclose (f);
            return 2;
        }
        // STRICT ARGUMENTS, as clips/forensics/lowend already are. The shared parse at the top of main()
        // uses atoi/atof, which read "1.5" as 1 channel and "48000Hz" as 48000 — and then measure, happily,
        // the wrong thing. These three modes are new in this release, so tightening them breaks nothing;
        // the wasm harness refuses the same strings, and a refusal that does not match on both roads is a
        // parity break that a diff of two successful runs would never show.
        {
            double sRate = 0.0; std::uint64_t sWidth = 0;
            if (! parseRate (argv[2], sRate) || ! parseCount (argv[3], sWidth)
                || sWidth < 1 || sWidth > (std::uint64_t) core::kMaxChannels)
            {
                std::fprintf (stderr, "bad sampleRate/channels\n");
                std::fclose (f);
                return 2;
            }
        }
        // analysis::BandBursts at its defaults: bursts in 5-9 kHz against the MEDIAN of a 2 s trailing
        // ring of 10 ms hops. Everything float is a RAW IEEE-754 BIT PATTERN, as `blocks` does it, so a
        // later wasm comparison catches a flipped bit that decimal printing would round away. The
        // histograms print only their NON-ZERO bins, which is a complete description of an integer
        // histogram and keeps a quiet file's output short.
        // K3 — THE BAND IS AN ARGUMENT NOW, and it has to be one HERE as well as in the wasm ABI. A road the
        // parity harness cannot drive is a road with no gate: `fc_probe_bursts_run_with` would have shipped
        // with its only evidence being that it compiles. Every flag is optional and defaults to the documented
        // value, so a command written before this still produces the same bytes.
        analysis::BandBursts det;
        analysis::BandBurstsParams bp;                    // the documented defaults, unless a flag moves one
        {
            struct Flag { const char* name; double* into; };
            const Flag flags[] = {
                { "--band-low",    &bp.bandLowHz },  { "--band-high",  &bp.bandHighHz },
                { "--hop-ms",      &bp.hopMs },      { "--baseline-ms", &bp.baselineMs },
                { "--enter-db",    &bp.enterDb },    { "--exit-db",    &bp.exitDb },
            };
            for (int i = 5; i < argc; ++i)
            {
                const Flag* hit = nullptr;
                for (const Flag& fl : flags) if (std::strcmp (argv[i], fl.name) == 0) { hit = &fl; break; }
                if (hit != nullptr)
                {
                    // A FLAG WITHOUT ITS VALUE IS A REFUSAL, not a silently kept default: `--enter-db` at the
                    // end of a command line would otherwise measure at 6 dB while the operator believes it
                    // set something. The same reasoning as the strict parse above.
                    if (i + 1 >= argc || ! parseFinite (argv[i + 1], *hit->into))
                    {
                        std::fprintf (stderr, "bursts: %s needs a finite number\n", hit->name);
                        std::fclose (f);
                        return 2;
                    }
                    ++i;            // the value is consumed; it must not be read as a name on the next turn
                    continue;
                }
                // AND A NAME NOBODY KNOWS IS A REFUSAL TOO. Without this the first version of the loop simply
                // did not match `--band-lo` and measured the DEFAULT 5-9 kHz band at exit 0, while the operator
                // read the command line and believed it had asked for 80 Hz. A typo that measures the wrong
                // thing silently is worse than one that measures nothing, and this is the same argument that
                // made a flag without its value a refusal — it was just applied to half the cases.
                if (std::strncmp (argv[i], "--", 2) == 0 && std::strcmp (argv[i], "--precise") != 0)
                {
                    std::fprintf (stderr, "bursts: unknown option %s (want --band-low --band-high --hop-ms "
                                          "--baseline-ms --enter-db --exit-db)\n", argv[i]);
                    std::fclose (f);
                    return 2;
                }
            }
        }
        det.setParams (bp);
        if (! det.prepare (fs, kChunk, nc))
        {
            // EVERY VALUE THE CORE JUDGED, not a guess at which one it disliked. The first version of this
            // named the band and the baseline, so `--enter-db 0` — refused because a threshold must be
            // positive and at least the exit one — reported a band problem that did not exist. The core
            // does not say which bound it hit; the honest thing is to print what it was handed.
            std::fprintf (stderr, "bursts: prepare refused this configuration — band %g-%g Hz (0.49 fs must "
                                  "clear the top corner), hop %g ms, baseline %g ms (>= one hop), "
                                  "enter %g dB (> 0), exit %g dB (0 <= exit <= enter), at %g Hz\n",
                          bp.bandLowHz, bp.bandHighHz, bp.hopMs, bp.baselineMs, bp.enterDb, bp.exitDb, fs);
            std::fclose (f);
            return 2;
        }
        streamPlanar (f, nc, [&] (const float* const* p, int n) { (void) det.process (p, nc, n); });
        std::fclose (f);
        det.finish();
        if ((std::uint64_t) det.samplesProcessed() != declaredFrames)
        {
            std::fprintf (stderr, "read %lld of %llu frames — refusing to report a partial measurement\n",
                          (long long) det.samplesProcessed(), (unsigned long long) declaredFrames);
            return 2;
        }

        std::printf ("# fcore bursts v1 sr=%016llx ch=%d hop=%d base=%d lo=%016llx hi=%016llx "
                     "enter=%016llx exit=%016llx chunk=%d\n",
                     (unsigned long long) bits (fs), nc, det.hopSamples(), det.baselineHops(),
                     (unsigned long long) bits (det.bandLowHz()), (unsigned long long) bits (det.bandHighHz()),
                     (unsigned long long) bits (bp.enterDb), (unsigned long long) bits (bp.exitDb), kChunk);
        std::printf ("samples %lld\n", (long long) det.samplesProcessed());
        std::printf ("hops %lld %lld %lld %lld\n", (long long) det.hopCount(), (long long) det.eligibleHops(),
                     (long long) det.zeroBaselineHops(), (long long) det.burstHops());
        std::printf ("damage %lld %lld %lld\n", (long long) det.damagedHops(),
                     (long long) det.overflowSamples(), (long long) det.firstNonFiniteAt());
        std::printf ("tail %lld %016llx\n", (long long) det.tailPartialSamples(),
                     (unsigned long long) bits (det.tailPartialEnergy()));
        std::printf ("valid %d %d %d %d\n", det.eventsValid() ? 1 : 0, (int) det.eventsInvalidReason(),
                     det.programmeEnergyValid() ? 1 : 0, (int) det.programmeEnergyInvalidReason());
        for (int c = 0; c < det.channels(); ++c)
            std::printf ("chan %d %016llx %lld %lld\n", c, (unsigned long long) bits (det.bandEnergy (c)),
                         (long long) det.nonFiniteSamples (c), (long long) det.absentSamples (c));
        std::printf ("events %lld %lld %d\n", (long long) det.eventCount(),
                     (long long) det.storedEventCount(), det.eventsComplete() ? 1 : 0);
        for (std::int64_t i = 0; i < det.storedEventCount(); ++i)
        {
            const analysis::BandBurst e = det.event (i);
            std::printf ("e %lld %lld %lld %lld %016llx %016llx %016llx %016llx %016llx %d%d%d\n",
                         (long long) e.start, (long long) e.length, (long long) e.peakAt, (long long) e.hops,
                         (unsigned long long) bits (e.peakPower), (unsigned long long) bits (e.peakBaseline),
                         (unsigned long long) bits (e.peakExcessDb), (unsigned long long) bits (e.peakWidePower),
                         (unsigned long long) bits (e.energy),
                         e.touchedNonFinite ? 1 : 0, e.baselineTouchedNonFinite ? 1 : 0,
                         e.closedByFinish ? 1 : 0);
        }
        std::printf ("onsets %lld %lld %lld %d %lld %016llx\n", (long long) det.onsetCount(),
                     (long long) det.intervalCount(), (long long) det.intervalOverflow(),
                     det.modalIntervalHops(), (long long) det.modalIntervalMass(),
                     (unsigned long long) bits (det.onsetsPerSecond()));
        for (int b = 1; b <= analysis::BandBursts::kIoiBins; ++b)
            if (det.intervalBin (b) != 0) std::printf ("ioi %d %lld\n", b, (long long) det.intervalBin (b));
        for (int b = 1; b <= analysis::BandBursts::kMaxLag; ++b)
            if (det.lagBin (b) != 0) std::printf ("lag %d %lld\n", b, (long long) det.lagBin (b));
        return 0;
    }

    if (mode == "forensics")
    {
        // The whole report, as bit patterns. The analyzer's own defaults are used and PRINTED, so a diff
        // between two toolchains compares the same instrument and not two configurations of it.
        // Strict where the top-level parse is lenient, like the shape modes: atoi("4294967297") narrows to
        // one channel and atof("48000Hz") reads 48000, and both would measure silently.
        double frate = 0.0;
        std::uint64_t fwidth = 0;
        if (! parseRate (argv[2], frate) || ! parseCount (argv[3], fwidth)
            || fwidth < 1 || fwidth > (std::uint64_t) core::kMaxChannels)
        {
            std::fprintf (stderr, "bad sampleRate/channels\n");
            std::fclose (f);
            return 2;
        }
        // The file is SIZED before it is read, so an input that is not a whole number of frames, or a read
        // that fails halfway, cannot come out as a successful measurement of a shorter programme.
        std::uint64_t fframes = 0;
        if (! fileFrames (f, nc, fframes))
        {
            std::fprintf (stderr, "cannot size the file, or it is not a whole number of %d-channel float32 frames\n", nc);
            std::fclose (f);
            return 2;
        }
        analysis::SourceForensics fx;
        const analysis::SourceForensicsParams fp;
        if (! fx.prepare (frate, kChunk, nc))
        {
            std::fprintf (stderr, "forensics.prepare refused (sample rate %g..%g)\n",
                          analysis::SourceForensics::kMinSampleRate, analysis::SourceForensics::kMaxSampleRate);
            std::fclose (f);
            return 2;
        }
        bool okAll = true;
        streamPlanar (f, nc, [&] (const float* const* p, int n) { okAll = fx.process (p, nc, n) && okAll; });
        const bool readError = std::ferror (f) != 0;
        std::fclose (f);
        if (! okAll) { std::fprintf (stderr, "forensics: a block was refused\n"); return 2; }
        if (readError) { std::fprintf (stderr, "forensics: the file could not be read to its end\n"); return 2; }
        fx.finish();
        if ((std::uint64_t) fx.samplesProcessed() != fframes)
        {
            std::fprintf (stderr, "forensics: the file delivered %lld of the %llu frames it was sized for\n",
                          (long long) fx.samplesProcessed(), (unsigned long long) fframes);
            return 2;
        }
        // EVERY parameter, so two builds that differ only in a threshold cannot print the same report.
        std::printf ("# fcore forensics v1 sr=%016llx ch=%d order=%d hop=%lld bins=%d percell=%d\n",
                     (unsigned long long) bits (frate), nc, fp.fftOrder, (long long) fx.hopSamples(),
                     fx.bins(), fx.binsPerCell());
        std::printf ("params exempt=%d distinctlimit=%d plateaucells=%d floorcells=%d", fx.exemptCells(),
                     fx.distinctLimit(), fx.plateauSpanCells(), fx.floorSpanCells());
        const double ps[] = { fx.cellHz(), fx.binHz(), fx.searchFromHz(), fx.searchToHz(), fp.cellWidthHz,
                              fp.searchFromHz, fp.plateauSpanHz, fp.floorSpanHz, fp.transitionStartDb,
                              fp.transitionEndDb, fp.minDropDb, fp.maxTransitionHz, fp.nearNyquistFraction,
                              fp.emptyDb, fp.emptyMinHz, fp.gridOutlierFraction };
        for (double d : ps) std::printf (" %016llx", (unsigned long long) bits (d));
        std::printf ("\n");
        std::printf ("samples %lld tail %lld frames %lld\n", (long long) fx.samplesProcessed(),
                     (long long) fx.tailUncoveredSamples(), (long long) fx.frames().frameCount());
        for (int c = 0; c <= nc; ++c)                     // per channel, then the file's own aggregate
        {
            const analysis::SpectralWall w = c < nc ? fx.wall (c) : fx.wall();
            std::printf ("wall %d valid=%d reason=%d sharp=%d nearnyq=%d clipped=%d trunc=%d exempted=%d"
                         " second=%d/%d/%d/%d secondreason=%d empty=%d emptyreason=%d frames=%lld/%lld\n",
                         c, (int) w.valid, (int) w.reason, (int) w.sharp, (int) w.nearNyquist,
                         (int) w.transitionClipped, (int) w.truncatedAtNyquist, w.exemptedCells,
                         (int) w.secondValid, (int) w.secondSharp, (int) w.secondTransitionClipped,
                         (int) w.secondTruncatedAtNyquist, (int) w.secondReason,
                         (int) w.emptyAboveValid, (int) w.emptyAboveReason,
                         (long long) w.framesUsed, (long long) w.framesHoled);
            const double ds[] = { w.cutoffHz, w.cutoffFractionOfNyquist, w.steepestHz, w.transitionEndHz,
                                  w.transitionHz, w.plateauPower, w.floorLocalPower, w.maxAbovePower,
                                  w.sufMaxPower, w.dropDb, w.strictDropDb, w.localDropDb, w.recoveryDb,
                                  w.plateauSpreadDb, w.steepnessDbPerOctave, w.secondCutoffHz, w.secondDropDb,
                                  w.secondTransitionHz, w.emptyAboveHz, w.emptyAboveFractionOfNyquist,
                                  w.emptyThresholdPower, w.peakCellPower };
            std::printf ("wall %d values", c);
            for (double d : ds) std::printf (" %016llx", (unsigned long long) bits (d));
            std::printf ("\n");
        }
        for (int c = 0; c < nc; ++c)
        {
            const analysis::SampleGrid g = fx.sampleGrid (c);
            std::printf ("grid %d valid=%d reason=%d k=%d pcm=%d outofrange=%d bits=%d robustk=%d"
                         " robustbits=%d zerolow24=%d peak=%016llx min=%016llx max=%016llx\n",
                         c, (int) g.valid, (int) g.reason, g.gridExponent, (int) g.pcmCompatible,
                         (int) g.outsidePcmRange, g.minExactPcmBits, g.robustGridExponent, g.robustPcmBits,
                         g.alwaysZeroLowBits (24), (unsigned long long) bits (g.absPeak),
                         (unsigned long long) bits (g.sampleMin), (unsigned long long) bits (g.sampleMax));
            std::printf ("grid %d nonzero=%lld zero=%lld nonfinite=%lld absent=%lld offgrid=%lld"
                         " firstoffgrid=%lld firstmaxk=%lld distinct=%lld complete=%d\n",
                         c, (long long) g.nonZeroSamples, (long long) g.zeroSamples,
                         (long long) g.nonFiniteSamples, (long long) g.absentSamples,
                         (long long) g.offGridSamples, (long long) g.firstOffGridSample,
                         (long long) g.firstMaxGridSample, (long long) g.distinctValues,
                         (int) g.distinctComplete);
            std::printf ("grid %d khist", c);
            const std::int64_t* h = fx.gridExponentHistogram (c);
            for (int k = 0; k < analysis::SourceForensics::gridExponentBuckets(); ++k)
                std::printf (" %lld", (long long) h[(std::size_t) k]);
            std::printf ("\n");
        }
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
            std::fprintf (stderr, "shapes.prepare refused (sample rate %g..%g, buckets/columns 1..%d, mix avr|L|R|max)\n",
                          fcore::Probe::kMinSampleRate, fcore::Probe::kMaxSampleRate, analysis::WaveformPeaks::kMaxBuckets);
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
        std::fprintf (stderr, "probe.prepare refused (sample rate %g..%g)\n",
                      fcore::Probe::kMinSampleRate, fcore::Probe::kMaxSampleRate);
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
