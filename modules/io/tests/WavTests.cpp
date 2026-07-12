// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::io — WAV write→read roundtrips at every supported depth,
// deinterleave, clamp/NaN handling, ragged-channel guard, and loud rejection of malformed / unsupported
// input (no silent zeros). Files go to the ctest working dir and are removed on the way out.

#include <felitronics_test.h>
#include <felitronics/io/Wav.h>

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace
{
std::vector<double> sine (std::size_t n, double f, double sr, double amp)
{
    std::vector<double> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = amp * std::sin (2.0 * M_PI * f * (double) i / sr);
    return x;
}

double maxErr (const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return 1e9;
    double e = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) e = std::max (e, std::fabs (a[i] - b[i]));
    return e;
}
} // namespace

int main()
{
    std::printf ("felitronics::io tests\n");
    const std::string tmp = "felitronics_io_tmp.wav";

    group ("f32 mono roundtrip");
    {
        const auto x = sine (4801, 997.0, 48000.0, 0.8);
        ok (io::writeWavMonoF32 (tmp, x, 48000.0), "writeWavMonoF32 succeeds");
        const auto w = io::readWav (tmp);
        ok (w.ok, "read ok");
        ok (w.ch.size() == 1 && w.frames() == x.size(), "1 ch, frame count preserved");
        approx (w.sr, 48000.0, 0.0, "sample rate preserved");
        ok (w.bits == 32 && w.is_float, "32-bit float flags");
        ok (maxErr (w.ch[0], x) < 1e-7, "float roundtrip within float precision");
    }

    group ("24-bit PCM stereo roundtrip + deinterleave");
    {
        const auto l = sine (2400, 440.0, 48000.0, 0.9);
        const auto r = sine (2400, 1234.0, 48000.0, 0.25);
        ok (io::writeWav (tmp, { l, r }, 48000.0, 24, false), "write 24-bit stereo");
        const auto w = io::readWav (tmp);
        ok (w.ok && w.ch.size() == 2 && w.frames() == 2400, "read ok, 2 ch");
        ok (w.bits == 24 && ! w.is_float, "24-bit PCM flags");
        ok (maxErr (w.ch[0], l) < 3e-7 && maxErr (w.ch[1], r) < 3e-7, "both channels within 24-bit LSB");
        ok (maxErr (w.ch[0], r) > 1e-2, "channels are NOT swapped/mixed");
    }

    group ("16-bit PCM roundtrip");
    {
        const auto x = sine (1000, 440.0, 44100.0, 0.7);
        ok (io::writeWav (tmp, { x }, 44100.0, 16, false), "write 16-bit");
        const auto w = io::readWav (tmp);
        ok (w.ok && w.bits == 16, "read ok");
        approx (w.sr, 44100.0, 0.0, "44.1k preserved");
        ok (maxErr (w.ch[0], x) < 4e-5, "within 16-bit LSB");
    }

    group ("24-bit sign extension");
    {
        const std::vector<double> x = { -0.5, -1.0, -1.0 / 8388608.0, 0.5 };
        ok (io::writeWav (tmp, { x }, 48000.0, 24, false), "write negatives");
        const auto w = io::readWav (tmp);
        ok (w.ok, "read ok");
        approx (w.ch[0][0], -0.5, 3e-7, "-0.5 survives");
        ok (w.ch[0][1] < -0.999, "-1.0 survives (no wrap to +)");
        ok (w.ch[0][2] < 0.0, "tiny negative stays negative (sign-extended)");
    }

    group ("clamp + non-finite on write");
    {
        const std::vector<double> x = { 2.0, -3.0, std::nan (""), INFINITY, 0.0 };
        ok (io::writeWav (tmp, { x }, 48000.0, 24, false), "write extremes");
        const auto w = io::readWav (tmp);
        ok (w.ok, "read ok");
        ok (w.ch[0][0] > 0.999 && w.ch[0][1] < -0.999, "±clamped to full scale");
        approx (w.ch[0][2], 0.0, 1e-12, "NaN → 0");
        approx (w.ch[0][3], 0.0, 1e-12, "Inf → 0");
    }

    group ("ragged channels + degenerate writes");
    {
        const std::vector<double> a (10, 0.1), b (5, 0.2);
        ok (io::writeWav (tmp, { a, b }, 48000.0, 16, false), "ragged write succeeds");
        const auto w = io::readWav (tmp);
        ok (w.ok && w.frames() == 5, "trimmed to shortest channel (no OOB)");
        ok (! io::writeWav (tmp, {}, 48000.0, 16, false), "no channels → false");
        ok (! io::writeWav (tmp, { std::vector<double>{} }, 48000.0, 16, false), "zero frames → false");
        ok (! io::writeWav (tmp, { a }, 48000.0, 8, false), "8-bit write unsupported → false");
        ok (! io::writeWav (tmp, { a }, 48000.0, 64, true), "f64 write unsupported → false");
    }

    group ("loud rejects (no silent zeros)");
    {
        ok (! io::readWav ("felitronics_io_definitely_missing.wav").ok, "missing file rejected");
        const std::uint8_t garbage[64] = { 'n', 'o', 't', 'a', 'w', 'a', 'v' };
        ok (! io::readWavMemory (garbage, sizeof (garbage)).ok, "garbage rejected");
        ok (! io::readWavMemory (garbage, 10).ok, "truncated (<44 bytes) rejected");

        // Hand-built canonical header claiming 8-bit PCM (unsupported) with a 4-byte data chunk.
        std::vector<std::uint8_t> u8;
        auto tag = [&] (const char* t) { u8.insert (u8.end(), t, t + 4); };
        auto u32 = [&] (std::uint32_t v) { for (int i = 0; i < 4; ++i) u8.push_back ((v >> (8 * i)) & 0xFF); };
        auto u16 = [&] (std::uint16_t v) { for (int i = 0; i < 2; ++i) u8.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); };
        tag ("RIFF"); u32 (36 + 4); tag ("WAVE");
        tag ("fmt "); u32 (16); u16 (1); u16 (1); u32 (48000); u32 (48000); u16 (1); u16 (8);
        tag ("data"); u32 (4); u8.insert (u8.end(), { 0x80, 0x80, 0x80, 0x80 });
        const auto w = io::readWavMemory (u8.data(), u8.size());
        ok (! w.ok && ! w.error.empty(), "8-bit PCM read rejected loudly");
    }


    group ("refutations from theory (crew round 1)");
    {
        auto tag = [] (std::vector<std::uint8_t>& o, const char* t) { o.insert (o.end(), t, t + 4); };
        auto u32 = [] (std::vector<std::uint8_t>& o, std::uint32_t v) { for (int i = 0; i < 4; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); };
        auto u16 = [] (std::vector<std::uint8_t>& o, std::uint16_t v) { for (int i = 0; i < 2; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); };

        // SPEC: a data chunk claiming more bytes than exist is corruption -> loud reject, no clamp.
        std::vector<std::uint8_t> trunc;
        tag (trunc, "RIFF"); u32 (trunc, 100); tag (trunc, "WAVE");
        tag (trunc, "fmt "); u32 (trunc, 16); u16 (trunc, 1); u16 (trunc, 1); u32 (trunc, 48000); u32 (trunc, 96000); u16 (trunc, 2); u16 (trunc, 16);
        tag (trunc, "data"); u32 (trunc, 4); trunc.push_back (0); trunc.push_back (0);
        ok (! io::readWavMemory (trunc.data(), trunc.size()).ok, "truncated data chunk rejected (never clamped)");

        // SPEC: EXTENSIBLE must carry its 40-byte body — a 16-byte fmt with tag 0xFFFE would read
        // the subtype out of the NEXT chunk's bytes.
        std::vector<std::uint8_t> ext;
        tag (ext, "RIFF"); u32 (ext, 100); tag (ext, "WAVE");
        tag (ext, "fmt "); u32 (ext, 16); u16 (ext, 0xFFFE); u16 (ext, 1); u32 (ext, 48000); u32 (ext, 96000); u16 (ext, 2); u16 (ext, 16);
        tag (ext, "data"); u32 (ext, 4); u16 (ext, 1); u16 (ext, 0);
        ok (! io::readWavMemory (ext.data(), ext.size()).ok, "undersized EXTENSIBLE fmt rejected");

        // SPEC: the writer refuses what the RIFF header cannot represent.
        std::vector<std::vector<double>> many (70000, std::vector<double> (1, 0.0));
        ok (! io::writeWav (tmp, many, 48000.0, 16, false), "u16 channel-count overflow refused");
        ok (! io::writeWav (tmp, { { 0.1, 0.2 } }, std::nan (""), 16, false), "non-finite sample rate refused");
        ok (! io::writeWav (tmp, { { 0.1 } }, -1.0, 16, false), "negative sample rate refused");
    }

    std::remove (tmp.c_str());
    return felitronics::test::report();
}
