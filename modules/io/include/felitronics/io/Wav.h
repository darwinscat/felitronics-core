// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::io — minimal, self-contained WAV read/write (no libsndfile, no JUCE). Moved verbatim
// from OrbitCapture's oc/wav.hpp so every Darwin's Cat product shares ONE implementation (OrbitCapture
// keeps an oc:: passthrough shim; WavData's field names stay as they were so that shim is a pure alias).
//
// Reads canonical PCM/float WAVE (16/24/32-bit int, 32/64-bit float, incl. WAVE_FORMAT_EXTENSIBLE),
// deinterleaves to per-channel doubles in [-1,1]. Rejects unsupported formats loudly (ok=false + error
// string) — never silent zeros. 512 MB DoS guard. Writes 16/24-bit PCM or 32-bit float; clamps,
// zeroes NaN/Inf, guards ragged channels. Little-endian (the WAV convention).
// OFFLINE, MESSAGE-THREAD-ONLY (allocates, touches the filesystem).
//==============================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace felitronics::io
{

struct WavData
{
    std::vector<std::vector<double>> ch;  // ch[c][n], normalized to [-1,1]
    double      sr       = 0.0;
    int         bits     = 0;
    bool        is_float = false;
    bool        ok       = false;
    std::string error;
    std::size_t frames() const { return ch.empty() ? 0 : ch[0].size(); }
};

namespace detail
{
    inline std::uint32_t rdU32 (const std::uint8_t* p)
    {
        return (std::uint32_t) p[0] | ((std::uint32_t) p[1] << 8) |
               ((std::uint32_t) p[2] << 16) | ((std::uint32_t) p[3] << 24);
    }
    inline std::uint16_t rdU16 (const std::uint8_t* p) { return (std::uint16_t) (p[0] | (p[1] << 8)); }

    inline void wrU32 (std::vector<std::uint8_t>& o, std::uint32_t v)
    {
        o.push_back (v & 0xFF); o.push_back ((v >> 8) & 0xFF);
        o.push_back ((v >> 16) & 0xFF); o.push_back ((v >> 24) & 0xFF);
    }
    inline void wrU16 (std::vector<std::uint8_t>& o, std::uint16_t v)
    {
        o.push_back (v & 0xFF); o.push_back ((v >> 8) & 0xFF);
    }
    inline void wrTag (std::vector<std::uint8_t>& o, const char* t) { o.insert (o.end(), t, t + 4); }
} // namespace detail

// Parse a whole WAVE image from memory (the file-reading wrapper is below). Exposed separately so a
// consumer with an embedded asset (BinaryData) decodes without a filesystem round-trip.
inline WavData readWavMemory (const std::uint8_t* bytes, std::size_t size)
{
    using namespace detail;
    WavData w;
    if (bytes == nullptr || size < 44) { w.error = "too small"; return w; }
    if (std::memcmp (bytes, "RIFF", 4) != 0 || std::memcmp (bytes + 8, "WAVE", 4) != 0)
    {
        w.error = "not a WAVE file"; return w;
    }
    std::size_t pos = 12;
    std::uint16_t fmt = 0, nch = 0, bits = 0; std::uint32_t rate = 0;
    const std::uint8_t* data = nullptr; std::size_t datalen = 0;
    while (pos + 8 <= size)
    {
        const std::uint8_t* id = bytes + pos;
        const std::uint32_t clen = rdU32 (bytes + pos + 4);
        const std::size_t body = pos + 8;
        // A chunk claiming more bytes than the file holds is corruption — reject LOUDLY (never
        // clamp: a clamped data chunk silently truncates audio). The checked advance also kills
        // the `body + clen` wrap that could loop forever on hostile input on 32-bit size_t.
        if (clen > size - body) { w.error = "corrupt WAV (chunk overruns the file)"; return w; }
        if (std::memcmp (id, "fmt ", 4) == 0 && clen >= 16)
        {
            fmt  = rdU16 (bytes + body);
            nch  = rdU16 (bytes + body + 2);
            rate = rdU32 (bytes + body + 4);
            bits = rdU16 (bytes + body + 14);
            if (fmt == 0xFFFE)
            {
                // WAVE_FORMAT_EXTENSIBLE: the subtype lives INSIDE this chunk (>= 40-byte body,
                // cbSize >= 22) and the SubFormat GUID tail must be the standard media GUID —
                // anything else is a crafted/corrupt header, not a format to trust.
                static const std::uint8_t kGuidTail[14] =
                    { 0x00,0x00, 0x00,0x00, 0x10,0x00, 0x80,0x00, 0x00,0xAA,0x00,0x38,0x9B,0x71 };
                if (clen < 40 || rdU16 (bytes + body + 16) < 22
                    || std::memcmp (bytes + body + 26, kGuidTail, sizeof (kGuidTail)) != 0)
                {
                    w.error = "corrupt WAV (bad WAVE_FORMAT_EXTENSIBLE header)"; return w;
                }
                fmt = rdU16 (bytes + body + 24);
            }
        }
        else if (std::memcmp (id, "data", 4) == 0)
        {
            data = bytes + body;
            datalen = clen;
        }
        const std::size_t advance = (std::size_t) clen + (clen & 1u);
        if (advance > size - body) break;                 // odd-length final chunk: parsed, done
        pos = body + advance;
    }
    if (! data || nch == 0 || bits == 0) { w.error = "missing fmt/data"; return w; }
    // Reject unsupported formats loudly — do NOT return silent zeros with ok=true.
    const bool supported = (fmt == 3 && (bits == 32 || bits == 64)) ||
                           (fmt == 1 && (bits == 16 || bits == 24 || bits == 32));
    if (! supported)
    {
        w.error = "unsupported WAV format (fmt=" + std::to_string (fmt) +
                  ", bits=" + std::to_string (bits) + ")";
        return w;
    }
    w.sr = rate; w.bits = bits; w.is_float = (fmt == 3);
    const std::size_t bytesPer = (std::size_t) bits / 8;
    const std::size_t frame = bytesPer * nch;
    if (frame == 0 || datalen % frame != 0) { w.error = "corrupt WAV (data not frame-aligned)"; return w; }
    const std::size_t nf = datalen / frame;
    w.ch.assign (nch, std::vector<double> (nf, 0.0));
    for (std::size_t i = 0; i < nf; ++i)
        for (std::uint16_t c = 0; c < nch; ++c)
        {
            const std::uint8_t* s = data + i * frame + c * bytesPer;
            double v = 0.0;
            if (fmt == 3 && bits == 32)      { float  fv; std::memcpy (&fv, s, 4); v = fv; }
            else if (fmt == 3 && bits == 64) { double dv; std::memcpy (&dv, s, 8); v = dv; }
            else if (bits == 16) { std::int16_t iv = (std::int16_t) rdU16 (s); v = iv / 32768.0; }
            else if (bits == 24)
            {
                std::int32_t iv = (std::int32_t) (s[0] | (s[1] << 8) | (s[2] << 16));
                if (iv & 0x800000) iv -= 0x1000000;   // sign-extend (no impl-defined cast)
                v = iv / 8388608.0;
            }
            else if (bits == 32) { std::int32_t iv = (std::int32_t) rdU32 (s); v = iv / 2147483648.0; }
            w.ch[c][i] = v;
        }
    w.ok = true; return w;
}

inline WavData readWav (const std::string& path)
{
    WavData w;
    std::FILE* f = std::fopen (path.c_str(), "rb");
    if (! f) { w.error = "cannot open " + path; return w; }
    std::fseek (f, 0, SEEK_END); long sz = std::ftell (f); std::fseek (f, 0, SEEK_SET);
    if (sz < 44) { std::fclose (f); w.error = "too small"; return w; }
    if (sz > (512L << 20)) { std::fclose (f); w.error = "file too large (>512 MB)"; return w; }  // DoS guard
    std::vector<std::uint8_t> b ((std::size_t) sz);
    if (std::fread (b.data(), 1, (std::size_t) sz, f) != (std::size_t) sz)
    {
        std::fclose (f); w.error = "read failed"; return w;
    }
    std::fclose (f);
    return readWavMemory (b.data(), b.size());
}

inline bool writeWav (const std::string& path, const std::vector<std::vector<double>>& ch,
                      double sr, int bits, bool is_float)
{
    using namespace detail;
    if (ch.empty()) return false;
    const bool supported = (is_float && bits == 32) || (! is_float && (bits == 16 || bits == 24));
    if (! supported) return false;
    std::size_t nf = ch[0].size();
    for (const auto& c : ch) nf = std::min (nf, c.size());   // guard ragged channels (no OOB)
    if (nf == 0) return false;
    // Refuse anything the RIFF header cannot represent — a silently truncated header field would
    // write a lying file (e.g. 70000 channels wrapping to a u16, NaN sample rates, >4 GB data).
    if (! std::isfinite (sr) || sr <= 0.0 || sr > 4294967295.0) return false;
    if (ch.size() > 0xFFFFu) return false;
    const std::uint64_t block64 = (std::uint64_t) ch.size() * (std::uint64_t) (bits / 8);
    const std::uint64_t data64  = (std::uint64_t) nf * block64;
    if (block64 > 0xFFFFu || data64 > 0xFFFFFFFFull - 44u) return false;
    const std::uint16_t nch = (std::uint16_t) ch.size();
    const std::uint16_t fmt = is_float ? 3 : 1;
    const std::uint32_t rate = (std::uint32_t) std::llround (sr);
    if ((std::uint64_t) rate * block64 > 0xFFFFFFFFull) return false;   // byteRate is a u32 too
    const std::uint16_t block = (std::uint16_t) block64;
    const std::uint32_t datalen = (std::uint32_t) data64;
    std::vector<std::uint8_t> o; o.reserve (44 + datalen);
    wrTag (o, "RIFF"); wrU32 (o, 36 + datalen); wrTag (o, "WAVE");
    wrTag (o, "fmt "); wrU32 (o, 16); wrU16 (o, fmt); wrU16 (o, nch);
    wrU32 (o, rate); wrU32 (o, rate * block); wrU16 (o, block); wrU16 (o, (std::uint16_t) bits);
    wrTag (o, "data"); wrU32 (o, datalen);
    auto clampd = [] (double v)
    {
        if (! std::isfinite (v)) return 0.0;                  // NaN/Inf → 0 (no llround UB)
        return v < -1.0 ? -1.0 : (v > 1.0 ? 1.0 : v);
    };
    for (std::size_t i = 0; i < nf; ++i)
        for (std::uint16_t c = 0; c < nch; ++c)
        {
            const double v = ch[c][i];
            if (is_float && bits == 32) { float fv = std::isfinite (v) ? (float) v : 0.0f; std::uint8_t t[4]; std::memcpy (t, &fv, 4); o.insert (o.end(), t, t + 4); }
            else if (bits == 16) { std::int16_t iv = (std::int16_t) std::llround (clampd (v) * 32767.0); wrU16 (o, (std::uint16_t) iv); }
            else if (bits == 24)
            {
                std::int32_t iv = (std::int32_t) std::llround (clampd (v) * 8388607.0);
                o.push_back (iv & 0xFF); o.push_back ((iv >> 8) & 0xFF); o.push_back ((iv >> 16) & 0xFF);
            }
        }
    std::FILE* f = std::fopen (path.c_str(), "wb");
    if (! f) return false;
    const bool okw = std::fwrite (o.data(), 1, o.size(), f) == o.size();
    std::fclose (f);
    return okw;
}

inline bool writeWavMonoF32 (const std::string& path, const std::vector<double>& x, double sr)
{
    return writeWav (path, { x }, sr, 32, true);
}

} // namespace felitronics::io
