// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::io's PCM GRID, held against the one felitronics::dither quantizes onto. Dither.h defines a
// `bits`-bit code as floor(v·2^(bits-1) + 0.5), clamped to [−2^(bits-1), 2^(bits-1) − 1], and readWav divides
// by 2^(bits-1). The writer used to scale by 2^(bits-1) − 1 under llround instead, so a master Dither had
// already put on the grid was quantized a SECOND time, onto a grid one LSB narrower and with no dither.
// Measured before the fix: 32767 of the 65536 16-bit codes came back moved by one LSB (every |k| > 16384 —
// 20000 was written as 19999, and -32768 as -32767), 8388607 of the 16777216 24-bit codes, and 30001 of the
// 48000 samples of a dithered 0.9 sine at 16 bit.
//
//   (1) write∘read is the identity on EVERY 16-bit code, and on 24-bit codes sampled across the range and
//       dense wherever the two grids part; read∘write is the identity on the bytes of a canonical image;
//   (2) Dither's own output, written and read back, is exactly the codes Dither chose — 16 and 24 bit, every
//       shaping, a programme that reaches both ends of the range;
//   (3) off the grid the writer rounds by Dither's rule: ties go UP on both signs, the clamp is in the double
//       domain (+1.0 lands on the top code, a huge finite value is not undefined behaviour), NaN/Inf → 0;
//   (4) THE PLANTED FAILURES: the old writer re-encoding the same payloads fails (1) and (2) by exactly the
//       counts the grid arithmetic predicts, and a writer on the right grid with llround's rounding fails (3)
//       on every negative tie — so a green run is the instruments seeing, not an instrument that is blind.
//
// Codes are read from the BYTES of the image, never only through readWav: the divisor there is the other half
// of the claim, and an instrument that went through it could not tell the two halves apart.

#include <felitronics_test.h>
#include <felitronics/dither/Dither.h>
#include <felitronics/io/Wav.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

namespace
{
constexpr double kPi = 3.14159265358979323846;

double full (int bits) { return bits == 16 ? 32768.0 : 8388608.0; }   // 2^(bits-1): the grid, both ways

// The codes an image carries, read from its bytes. Every image here is one channel behind the canonical
// 44-byte header, which is what the writer emits: the data chunk's size from its header, then RIFF's zero pad
// byte after it when that size is odd (a 24-bit image of an odd length). Anything else — including the pad
// missing — reads as empty and fails the count.
std::vector<long long> codesOf (const std::vector<std::uint8_t>& img, int bits)
{
    std::vector<long long> k;
    const std::size_t bp = (std::size_t) bits / 8;
    if (img.size() < 44 || std::memcmp (img.data() + 36, "data", 4) != 0) return k;
    const std::size_t len = (std::size_t) img[40] | ((std::size_t) img[41] << 8) | ((std::size_t) img[42] << 16) | ((std::size_t) img[43] << 24);
    if (len % bp != 0 || img.size() != 44 + len + (len & 1u) || ((len & 1u) && img.back() != 0)) return k;
    for (std::size_t off = 44; off < 44 + len; off += bp)
    {
        if (bits == 16) { k.push_back ((std::int16_t) (std::uint16_t) (img[off] | (img[off + 1] << 8))); continue; }
        std::int32_t v = img[off] | (img[off + 1] << 8) | (img[off + 2] << 16);
        if (v & 0x800000) v -= 0x1000000;
        k.push_back (v);
    }
    return k;
}

// A mono image holding exactly these codes, built without the writer — the reader's input for read∘write.
// Canonical, so an odd data chunk carries RIFF's pad byte, and the RIFF size counts it.
std::vector<std::uint8_t> imageOf (const std::vector<long long>& k, int bits)
{
    std::vector<std::uint8_t> o;
    auto u32 = [&] (std::uint32_t v) { for (int i = 0; i < 4; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); };
    auto u16 = [&] (std::uint16_t v) { for (int i = 0; i < 2; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); };
    auto tag = [&] (const char* t) { o.insert (o.end(), t, t + 4); };
    const std::uint32_t bp = (std::uint32_t) bits / 8, len = (std::uint32_t) k.size() * bp, pad = len & 1u;
    tag ("RIFF"); u32 (36 + len + pad); tag ("WAVE");
    tag ("fmt "); u32 (16); u16 (1); u16 (1); u32 (48000); u32 (48000 * bp); u16 ((std::uint16_t) bp); u16 ((std::uint16_t) bits);
    tag ("data"); u32 (len);
    for (long long c : k)
        for (std::uint32_t b = 0; b < bp; ++b) o.push_back ((std::uint8_t) (((unsigned long long) c >> (8 * b)) & 0xFF));
    if (pad) o.push_back (0);
    return o;
}

// A real image with its payload re-encoded by another quantizer, header untouched, so a planted writer is
// read through exactly the path the real one is.
template <typename Quantize>
std::vector<std::uint8_t> reencode (std::vector<std::uint8_t> img, const std::vector<double>& x, int bits, Quantize q)
{
    const std::size_t bp = (std::size_t) bits / 8;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const long long c = q (x[i]);
        for (std::size_t b = 0; b < bp; ++b) img[44 + i * bp + b] = (std::uint8_t) (((unsigned long long) c >> (8 * b)) & 0xFF);
    }
    return img;
}

// THE PLANTED FAILURE: the writer as it was — llround on the 2^(bits-1) − 1 grid.
std::vector<std::uint8_t> oldWriter (const std::vector<std::uint8_t>& img, const std::vector<double>& x, int bits)
{
    const double s = full (bits) - 1.0;
    return reencode (img, x, bits, [s] (double v)
    {
        const double c = std::isfinite (v) ? std::clamp (v, -1.0, 1.0) : 0.0;
        return (long long) std::llround (c * s);
    });
}

// The second one: the RIGHT grid, llround's rounding (half away from zero) instead of Dither's (half up).
std::vector<std::uint8_t> halfAwayWriter (const std::vector<std::uint8_t>& img, const std::vector<double>& x, int bits)
{
    const double f = full (bits);
    return reencode (img, x, bits, [f] (double v)
    {
        const double c = std::isfinite (v) ? std::clamp (v * f, -f, f - 1.0) : 0.0;
        return (long long) std::llround (c);
    });
}

// Dither.h's rule, restated from its header: the reference the writer is measured against off the grid.
long long ditherRule (double v, int bits)
{
    if (! std::isfinite (v)) return 0;
    const double f = full (bits), q = std::floor (v * f + 0.5);
    return q <= -f ? (long long) -f : q >= f - 1.0 ? (long long) (f - 1.0) : (long long) q;
}

// How many codes differ from what was meant (−1: the image did not even have the right length).
long long moved (const std::vector<long long>& got, const std::vector<long long>& want)
{
    if (got.size() != want.size()) return -1;
    long long n = 0;
    for (std::size_t i = 0; i < got.size(); ++i) n += got[i] != want[i];
    return n;
}

// How many samples readWav returns as anything but exactly want[i] / 2^(bits-1).
long long readMisses (const std::vector<std::uint8_t>& img, const std::vector<long long>& want, int bits)
{
    const auto w = io::readWavMemory (img.data(), img.size());
    if (! w.ok || w.ch.size() != 1 || w.frames() != want.size()) return -1;
    long long n = 0;
    for (std::size_t i = 0; i < want.size(); ++i) n += w.ch[0][i] != (double) want[i] / full (bits);
    return n;
}

// Where the old grid parts from this one: every |k| > 2^(bits-2) moves one LSB toward zero under it, and
// nothing else does (at |k| = 2^(bits-2) exactly, k ∓ 0.5 is a tie llround takes away from zero, back to k).
long long partedCount (const std::vector<long long>& k, int bits)
{
    const long long edge = 1LL << (bits - 2);
    return (long long) std::count_if (k.begin(), k.end(), [edge] (long long c) { return c > edge || c < -edge; });
}

// Every 16-bit code, in order.
std::vector<long long> codes16()
{
    std::vector<long long> k;
    for (long long c = -32768; c <= 32767; ++c) k.push_back (c);
    return k;
}

// 24-bit codes: a prime stride across the whole range, plus every code within 64 of 0, of ±2^21, of ±2^22
// (where the old grid began to part) and of both ends.
std::vector<long long> codes24()
{
    std::vector<long long> k;
    for (long long c = -8388608; c <= 8388607; c += 4099) k.push_back (c);
    for (long long centre : { 0LL, 2097152LL, -2097152LL, 4194304LL, -4194304LL })
        for (long long d = -64; d <= 64; ++d) k.push_back (centre + d);
    for (long long d = 0; d <= 64; ++d) { k.push_back (-8388608 + d); k.push_back (8388607 - d); }
    return k;
}

std::vector<double> valuesOf (const std::vector<long long>& k, int bits)
{
    std::vector<double> x;
    for (long long c : k) x.push_back ((double) c / full (bits));
    return x;
}

std::size_t indexOf (const std::vector<long long>& k, long long c)
{
    return (std::size_t) (std::find (k.begin(), k.end(), c) - k.begin());
}

// codes[indexOf (k, c)], or a value no code can be when the image did not decode — an instrument that read
// out of bounds would crash the suite instead of failing the check.
long long codeAt (const std::vector<long long>& codes, const std::vector<long long>& k, long long c)
{
    const std::size_t i = indexOf (k, c);
    return i < codes.size() ? codes[i] : std::numeric_limits<long long>::min();
}
} // namespace

int main()
{
    std::printf ("felitronics::io PCM grid tests\n");

    for (int bits : { 16, 24 })
    {
        const auto k = bits == 16 ? codes16() : codes24();
        const auto x = valuesOf (k, bits);
        const std::string at = " (" + std::to_string (bits) + "-bit, " + std::to_string (k.size()) + " codes)";

        group ("write∘read is the identity on every code" + at);
        {
            const auto img = io::writeWavMemory ({ x }, 48000.0, bits, false);
            ok (moved (codesOf (img, bits), k) == 0, "every code k/2^(bits-1) is written as k" + at);
            ok (readMisses (img, k, bits) == 0, "…and reads back as exactly k/2^(bits-1)" + at);

            // THE PLANTED FAILURE, through the same two instruments.
            const auto bad = oldWriter (img, x, bits);
            const auto badK = codesOf (bad, bits);
            const long long parted = partedCount (k, bits);
            ok (moved (badK, k) == parted && parted > 0,
                "planted: the old 2^(bits-1)-1 writer moves exactly the codes with |k| > 2^(bits-2)" + at);
            ok (readMisses (bad, k, bits) == parted, "planted: …and the reader sees every one of them" + at);
            if (bits == 16)
            {
                ok (parted == 32767, "planted: at 16 bit that is 32767 of the 65536 codes");
                ok (codeAt (badK, k, 20000) == 19999 && codeAt (badK, k, -32768) == -32767
                        && codeAt (badK, k, 16385) == 16384 && codeAt (badK, k, 16384) == 16384,
                    "planted: 20000 → 19999, -32768 → -32767, 16385 → 16384, and 16384 stays");
            }
            else
                ok (codeAt (badK, k, 4194305) == 4194304 && codeAt (badK, k, -8388608) == -8388607
                        && codeAt (badK, k, 4194304) == 4194304,
                    "planted: 4194305 → 4194304, -8388608 → -8388607, and 4194304 stays");
        }

        group ("read∘write is the identity on the bytes of a canonical image" + at);
        {
            const auto img = imageOf (k, bits);
            const auto w = io::readWavMemory (img.data(), img.size());
            ok (w.ok && w.bits == bits && w.frames() == k.size(), "the hand-built image reads" + at);
            ok (io::writeWavMemory (w.ch, w.sr, bits, false) == img, "…and writes back byte for byte" + at);
            ok (oldWriter (img, w.ch[0], bits) != img, "planted: the old writer does not give the bytes back" + at);
        }
    }

    // Dither has no amplitude control, and needs none for this: whatever noise it adds, its output is
    // code · 2^-(bits-1), and those codes are what the file has to carry. The rounding rule itself is pinned
    // off the grid by the next group.
    group ("Dither's output is written as the codes Dither chose");
    {
        // A hot 997 Hz sine, a ramp over the whole range, then plateaus at ±1.0 and past it: +0.09 dBFS is
        // 327 LSB over at 16 bit, more than the Psychoacoustic shaper's feedback can pull back (at most
        // Σ|h|·8 ≈ 171 LSB), so Dither reaches its clamp at both ends whatever its noise does.
        std::vector<float> programme;
        for (int i = 0; i < 48000; ++i) programme.push_back ((float) (0.9995 * std::sin (2.0 * kPi * 997.0 * i / 48000.0)));
        for (int i = 0; i < 48000; ++i) programme.push_back ((float) (-1.0 + 2.0 * i / 47999.0));
        for (float plateau : { -1.01f, 1.01f, -1.0f, 1.0f }) programme.insert (programme.end(), 512, plateau);
        const int n = (int) programme.size();
        const char* shapingName[] = { "None", "Weighted", "Psychoacoustic" };
        for (int bits : { 16, 24 })
            for (auto shaping : { dither::NoiseShaping::None, dither::NoiseShaping::Weighted, dither::NoiseShaping::Psychoacoustic })
            {
                const std::string at = " (" + std::to_string (bits) + "-bit, " + shapingName[(int) shaping] + ")";
                dither::Dither d;
                felitronics::test::run (d.prepare (48000.0, 4096, 1));
                dither::DitherParams p; p.bits = bits; p.shaping = shaping; d.setParams (p);
                auto y = programme;
                float* io[1] = { y.data() };
                felitronics::test::run (d.process (io, 1, n));

                std::vector<long long> k; long long offGrid = 0;
                for (float v : y)
                {
                    const double c = (double) v * full (bits);
                    offGrid += c != std::floor (c);
                    k.push_back ((long long) c);
                }
                ok (offGrid == 0, "precondition: Dither's output is on its grid" + at);
                ok (*std::min_element (k.begin(), k.end()) == (long long) -full (bits)
                        && *std::max_element (k.begin(), k.end()) == (long long) full (bits) - 1,
                    "precondition: the programme reaches both end codes" + at);

                const std::vector<double> x (y.begin(), y.end());
                const auto img = io::writeWavMemory ({ x }, 48000.0, bits, false);
                ok (moved (codesOf (img, bits), k) == 0, "every sample is written as Dither's code" + at);
                ok (readMisses (img, k, bits) == 0, "…and reads back as Dither's float, exactly" + at);

                const long long parted = partedCount (k, bits);
                ok (moved (codesOf (oldWriter (img, x, bits), bits), k) == parted && parted > n / 4,
                    "planted: the old writer re-quantizes every sample with |code| > 2^(bits-2), undithered" + at);
            }
    }

    group ("off the grid: Dither's rounding, Dither's clamp");
    {
        for (int bits : { 16, 24 })
        {
            const std::string at = " (" + std::to_string (bits) + "-bit)";
            // Every tie k + ½ at 16 bit; the 24-bit sample's. Half UP means k + 1 on BOTH signs — llround would
            // take a negative tie away from zero, to k.
            std::vector<long long> want; std::vector<double> ties;
            for (long long c : bits == 16 ? codes16() : codes24())
            {
                ties.push_back (((double) c + 0.5) / full (bits));
                want.push_back (std::min (c + 1, (long long) full (bits) - 1));
            }
            const auto img = io::writeWavMemory ({ ties }, 48000.0, bits, false);
            ok (moved (codesOf (img, bits), want) == 0, "every tie k + 1/2 is written as k + 1, both signs" + at);
            const long long negTies = (long long) std::count_if (want.begin(), want.end(), [] (long long c) { return c <= 0; });
            ok (moved (codesOf (halfAwayWriter (img, ties, bits), bits), want) == negTies && negTies > 0,
                "planted: llround's half-away rounding on the right grid misses exactly the negative ties" + at);
            ok (moved (codesOf (oldWriter (img, ties, bits), bits), want) > 0, "planted: the old writer misses ties too" + at);

            // Off-grid values in and past the range, from a fixed LCG, against the rule restated from Dither.h.
            std::vector<double> xs; std::vector<long long> ref;
            unsigned long long s = 0x9E3779B97F4A7C15ull;
            for (int i = 0; i < 20000; ++i)
            {
                s = s * 6364136223846793005ULL + 1442695040888963407ULL;
                const double v = ((double) (s >> 11) * 0x1.0p-53 * 2.0 - 1.0) * 1.05;
                xs.push_back (v); ref.push_back (ditherRule (v, bits));
            }
            const auto rimg = io::writeWavMemory ({ xs }, 48000.0, bits, false);
            ok (moved (codesOf (rimg, bits), ref) == 0, "20000 off-grid values in ±1.05 follow Dither's rule" + at);
            ok (moved (codesOf (oldWriter (rimg, xs, bits), bits), ref) > 0, "planted: the old writer does not" + at);

            // The ends. +1.0 is one LSB past the top code — two's complement has no +2^(bits-1) — and clamps to
            // it; −1.0 is the bottom code itself. A huge finite value clamps in the double domain (no UB in the
            // conversion); a non-finite one is 0, as it always was.
            const double f = full (bits);
            const long long top = (long long) f - 1, bottom = (long long) -f;
            const std::vector<double> ends = { 1.0, -1.0, 1e300, -1e300, DBL_MAX, -DBL_MAX,
                                               std::numeric_limits<double>::quiet_NaN(),
                                               std::numeric_limits<double>::infinity(),
                                               -std::numeric_limits<double>::infinity(),
                                               -0.0, std::numeric_limits<double>::denorm_min(),
                                               -std::numeric_limits<double>::denorm_min(),
                                               (f - 0.5) / f, (-f - 0.5) / f, (-f - 0.5000001) / f,
                                               (f - 1.5) / f, 0.5 / f, -0.5 / f };
            const std::vector<long long> endCodes = { top, bottom, top, bottom, top, bottom, 0, 0, 0,
                                                      0, 0, 0,
                                                      top, bottom, bottom,
                                                      top, 1, 0 };
            const auto eimg = io::writeWavMemory ({ ends }, 48000.0, bits, false);
            ok (moved (codesOf (eimg, bits), endCodes) == 0, "±1, ±huge, ±DBL_MAX, NaN, ±Inf, ±0, subnormals, the end ties" + at);
        }
    }

    return felitronics::test::report();
}
