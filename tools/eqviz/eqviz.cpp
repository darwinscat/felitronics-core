// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// eqviz — renders PNG "quality pictures" of every felitronics::eq filter type, so a reviewer can SEE
// filter quality instead of trusting their ears (the 2026-07 falsification campaign found 7 real bugs
// living in parameter-space corners nobody ever heard). Offline dev tool — no RT constraints.
//
// Per filter type it renders:
//   * a QUALITY HEATMAP: x = f0 (log), y = Q / order / gain (type-specific), colour = the worst
//     |digital − analog-reference| dB over the measurement region, MAXed with a per-type INVARIANT
//     penalty (centre gain, plateau, null depth, |H(fc)|=−3 dB, allpass flatness + group delay > 0)
//     so an excluded null region can never hide a broken design. Discrete bands: green < 0.1 dB,
//     yellow < 1, orange < 3, red ≥ 3, BLACK = unstable/non-finite design.
//   * a CURVE SHEET: representative configs — digital |H| vs the analog prototype vs (where honest)
//     the RBJ cookbook, an error track, and the MEASURED float audio path (adaptive-length impulse
//     through EqBand::processBlock → FFT) — catches "coefficients right, path wrong" bugs.
//
// Measurement region (consilium-agreed): f ∈ [20 Hz, 0.95·fs/2] (the last 5 % before Nyquist is the
// matched design's INTENDED deviation from analog — excluded by default so red always means "look"),
// analog reference > −50 dB, ±6 % around band-stop nulls excluded (the invariant covers the null).
// `--full` extends to 0.4999·fs for the honest edge picture.
//
// Usage: eqviz <outDir> [--fs 44100,96000] [--full]
// Writes PNGs + index.md + summary.csv into <outDir>. Exit code 1 on any unstable/non-finite design
// encountered (so CI can hard-fail on black pixels), 2 on usage errors.

#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqEngine.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics::eq;

//==============================================================================================
// Tiny PNG writer — 8-bit RGB, zlib STORED blocks (no compression; dev artifacts, determinism and
// zero third-party code beat file size here). PNG = sig, IHDR, IDAT (zlib: 0x78 0x01, stored-deflate
// blocks, Adler32), IEND; every chunk CRC32'd.
namespace png
{
    inline std::uint32_t crcTable (int i)
    {
        std::uint32_t c = (std::uint32_t) i;
        for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        return c;
    }
    inline std::uint32_t crc32 (const std::uint8_t* d, size_t n, std::uint32_t c = 0xFFFFFFFFu)
    {
        static std::uint32_t t[256]; static bool init = false;
        if (! init) { for (int i = 0; i < 256; ++i) t[i] = crcTable (i); init = true; }
        for (size_t i = 0; i < n; ++i) c = t[(c ^ d[i]) & 0xFF] ^ (c >> 8);
        return c;
    }
    inline void be32 (std::vector<std::uint8_t>& v, std::uint32_t x)
    {
        v.push_back ((std::uint8_t) (x >> 24)); v.push_back ((std::uint8_t) (x >> 16));
        v.push_back ((std::uint8_t) (x >> 8));  v.push_back ((std::uint8_t) x);
    }
    inline void chunk (std::vector<std::uint8_t>& out, const char* tag, const std::vector<std::uint8_t>& data)
    {
        be32 (out, (std::uint32_t) data.size());
        std::vector<std::uint8_t> body ((const std::uint8_t*) tag, (const std::uint8_t*) tag + 4);
        body.insert (body.end(), data.begin(), data.end());
        out.insert (out.end(), body.begin(), body.end());
        be32 (out, crc32 (body.data(), body.size()) ^ 0xFFFFFFFFu);
    }
    // rgb = w*h*3 bytes, row-major.
    inline bool write (const std::string& path, int w, int h, const std::vector<std::uint8_t>& rgb)
    {
        std::vector<std::uint8_t> raw;                                   // filter byte 0 + scanline
        raw.reserve ((size_t) h * ((size_t) w * 3 + 1));
        for (int y = 0; y < h; ++y)
        {
            raw.push_back (0);
            raw.insert (raw.end(), rgb.begin() + (size_t) y * w * 3, rgb.begin() + (size_t) (y + 1) * w * 3);
        }
        std::vector<std::uint8_t> idat { 0x78, 0x01 };                   // zlib header
        std::uint32_t s1 = 1, s2 = 0;                                    // Adler32
        for (std::uint8_t b : raw) { s1 = (s1 + b) % 65521u; s2 = (s2 + s1) % 65521u; }
        for (size_t off = 0; off < raw.size(); off += 65535)
        {
            const std::uint16_t len = (std::uint16_t) std::min<size_t> (65535, raw.size() - off);
            idat.push_back (off + len >= raw.size() ? 1 : 0);            // BFINAL, BTYPE=00 (stored)
            idat.push_back ((std::uint8_t) (len & 0xFF)); idat.push_back ((std::uint8_t) (len >> 8));
            const std::uint16_t nlen = (std::uint16_t) ~len;
            idat.push_back ((std::uint8_t) (nlen & 0xFF)); idat.push_back ((std::uint8_t) (nlen >> 8));
            idat.insert (idat.end(), raw.begin() + (long) off, raw.begin() + (long) (off + len));
        }
        be32 (idat, (s2 << 16) | s1);

        std::vector<std::uint8_t> out { 0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A };
        std::vector<std::uint8_t> ihdr;
        be32 (ihdr, (std::uint32_t) w); be32 (ihdr, (std::uint32_t) h);
        ihdr.push_back (8); ihdr.push_back (2);                          // 8-bit, truecolour RGB
        ihdr.push_back (0); ihdr.push_back (0); ihdr.push_back (0);
        chunk (out, "IHDR", ihdr);
        chunk (out, "IDAT", idat);
        chunk (out, "IEND", {});

        std::FILE* f = std::fopen (path.c_str(), "wb");
        if (! f) return false;
        const bool ok = std::fwrite (out.data(), 1, out.size(), f) == out.size();
        std::fclose (f);
        return ok;
    }
}

//==============================================================================================
// 5×7 bitmap font (caps, digits, the punctuation the plots use). Each glyph: 7 rows × 5 bits.
namespace font
{
    struct Glyph { char c; std::uint8_t r[7]; };
    static const Glyph kGlyphs[] = {
        {'0',{0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}},{'1',{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}},
        {'2',{0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}},{'3',{0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}},
        {'4',{0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}},{'5',{0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}},
        {'6',{0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}},{'7',{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}},
        {'8',{0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}},{'9',{0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}},
        {'A',{0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}},{'B',{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}},
        {'C',{0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}},{'D',{0x1C,0x12,0x11,0x11,0x11,0x12,0x1C}},
        {'E',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}},{'F',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}},
        {'G',{0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}},{'H',{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}},
        {'I',{0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}},{'J',{0x07,0x02,0x02,0x02,0x02,0x12,0x0C}},
        {'K',{0x11,0x12,0x14,0x18,0x14,0x12,0x11}},{'L',{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}},
        {'M',{0x11,0x1B,0x15,0x15,0x11,0x11,0x11}},{'N',{0x11,0x19,0x15,0x13,0x11,0x11,0x11}},
        {'O',{0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}},{'P',{0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}},
        {'Q',{0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}},{'R',{0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}},
        {'S',{0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}},{'T',{0x1F,0x04,0x04,0x04,0x04,0x04,0x04}},
        {'U',{0x11,0x11,0x11,0x11,0x11,0x11,0x0E}},{'V',{0x11,0x11,0x11,0x11,0x11,0x0A,0x04}},
        {'W',{0x11,0x11,0x11,0x15,0x15,0x1B,0x11}},{'X',{0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}},
        {'Y',{0x11,0x11,0x0A,0x04,0x04,0x04,0x04}},{'Z',{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}},
        {'-',{0x00,0x00,0x00,0x1F,0x00,0x00,0x00}},{'+',{0x00,0x04,0x04,0x1F,0x04,0x04,0x00}},
        {'.',{0x00,0x00,0x00,0x00,0x00,0x0C,0x0C}},{'/',{0x01,0x01,0x02,0x04,0x08,0x10,0x10}},
        {'(',{0x02,0x04,0x08,0x08,0x08,0x04,0x02}},{')',{0x08,0x04,0x02,0x02,0x02,0x04,0x08}},
        {'=',{0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}},{':',{0x00,0x0C,0x0C,0x00,0x0C,0x0C,0x00}},
        {'%',{0x19,0x19,0x02,0x04,0x08,0x13,0x13}},{',',{0x00,0x00,0x00,0x00,0x0C,0x04,0x08}},
        {'<',{0x02,0x04,0x08,0x10,0x08,0x04,0x02}},{'>',{0x08,0x04,0x02,0x01,0x02,0x04,0x08}},
        {' ',{0,0,0,0,0,0,0}},
    };
    inline const std::uint8_t* rows (char c)
    {
        if (c >= 'a' && c <= 'z') c = (char) (c - 'a' + 'A');
        for (const auto& g : kGlyphs) if (g.c == c) return g.r;
        return kGlyphs[sizeof (kGlyphs) / sizeof (Glyph) - 1].r;         // unknown -> space
    }
}

//==============================================================================================
// Canvas + drawing primitives.
struct RGB { std::uint8_t r, g, b; };
static const RGB kWhite {255,255,255}, kBlack {15,15,15}, kGrey {200,200,205}, kDarkGrey {110,110,115},
                 kGreen {58,168,82},  kYellow {247,201,49}, kOrange {247,126,33}, kRed {214,48,44},
                 kBlue {44,98,214},   kAnalog {30,160,70},  kRbj {150,150,155},   kMeas {214,48,120};

struct Canvas
{
    int w = 0, h = 0;
    std::vector<std::uint8_t> px;
    Canvas (int W, int H) : w (W), h (H), px ((size_t) W * H * 3, 255) {}
    void set (int x, int y, RGB c)
    {
        if (x < 0 || y < 0 || x >= w || y >= h) return;
        auto* p = &px[((size_t) y * w + x) * 3]; p[0] = c.r; p[1] = c.g; p[2] = c.b;
    }
    void fill (int x0, int y0, int x1, int y1, RGB c)
    { for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) set (x, y, c); }
    void line (int x0, int y0, int x1, int y1, RGB c)                     // Bresenham
    {
        int dx = std::abs (x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -std::abs (y1 - y0), sy = y0 < y1 ? 1 : -1, err = dx + dy;
        for (;;)
        {
            set (x0, y0, c);
            if (x0 == x1 && y0 == y1) break;
            const int e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }
    void text (int x, int y, const std::string& s, RGB c, int scale = 1)
    {
        for (char ch : s)
        {
            const std::uint8_t* g = font::rows (ch);
            for (int r = 0; r < 7; ++r) for (int b = 0; b < 5; ++b)
                if (g[r] & (1 << (4 - b)))
                    fill (x + b * scale, y + r * scale, x + b * scale + scale - 1, y + r * scale + scale - 1, c);
            x += 6 * scale;
        }
    }
};

//==============================================================================================
// Analog references (independent of the library code — the ground truth side of the pictures).
namespace ref
{
    // 2nd-order resonant low/high pass: |H(w0)| = Q.
    inline double lowpass (double f, double f0, double Q)
    {
        const double y = (f / f0) * (f / f0);
        return 1.0 / std::sqrt ((1.0 - y) * (1.0 - y) + y / (Q * Q));
    }
    inline double highpass (double f, double f0, double Q) { return lowpass (f0 * f0 / f, f0, Q) ; }
    // N-th order Butterworth LP/HP.
    inline double butterLP (double f, double f0, int N) { return 1.0 / std::sqrt (1.0 + std::pow (f / f0, 2.0 * N)); }
    inline double butterHP (double f, double f0, int N) { return 1.0 / std::sqrt (1.0 + std::pow (f0 / f, 2.0 * N)); }
    // Constant-Q peaking |H| of H(s) = (s² + (A/Q)ω0 s + ω0²)/(s² + (1/Q)ω0 s + ω0²), s = jΩ
    // (boost A ≥ 1; a cut is the exact reciprocal of the mirror boost, same as the matched design).
    inline double bellMag (double f, double f0, double Q, double A)
    {
        if (A < 1.0 && A > 0.0) return 1.0 / bellMag (f, f0, Q, 1.0 / A);
        const double r = f / f0, d = (1.0 - r * r);
        const double num = d * d + (A * r / Q) * (A * r / Q);
        const double den = d * d + (r / Q) * (r / Q);
        return std::sqrt (num / den);
    }
    // Resonant shelf prototype: poles (ω0, Q), zeros (ω0·√A, Q). k = 1/Q² − 2, y = (f/f0)².
    inline double shelfQ (double f, double f0, double Q, double A, bool high)
    {
        if (A < 1.0 && A > 0.0) return 1.0 / shelfQ (f, f0, Q, 1.0 / A, high);
        const double y = (f / f0) * (f / f0), k = 1.0 / (Q * Q) - 2.0;
        const double num = high ? (A * A * y * y + A * k * y + 1.0) : (y * y + A * k * y + A * A);
        const double den = y * y + k * y + 1.0;
        return std::sqrt (std::max (0.0, num) / std::max (1e-300, den));
    }
    // Non-resonant 2-pole Butterworth high shelf (Vicanek eq. 1): DC = 1, top = G.
    inline double butterHS (double f, double f0, double G)
    {
        const double sg = std::sqrt (G), r2 = (f / f0) * (f / f0);
        const double num = (1.0 - sg * r2) * (1.0 - sg * r2) + 2.0 * sg * r2;
        const double den = (1.0 - r2 / sg) * (1.0 - r2 / sg) + 2.0 * r2 / sg;
        return std::sqrt (num / den);
    }
    // Tilt about f0: lows −g, highs +g (product of the two non-resonant shelves; validated to 1e-13 dB).
    inline double tilt (double f, double f0, double G)                    // G = 10^(g/20)
    {
        const double fc4 = 1.0, x4 = std::pow (f / f0, 4.0);
        return (1.0 / G) * (fc4 + x4 * G) / (fc4 + x4 / G);
    }
    inline double bandpass (double f, double f0, double Q)
    {
        const double Om0 = f0, BW = f0 / Q, Om = f;
        const double d = Om0 * Om0 - Om * Om;
        return (BW * Om) / std::sqrt (d * d + BW * Om * BW * Om);
    }
    // Analog order-m Butterworth band-pass: |H(iΩ)|² = 1 / (1 + ((Ω0²−Ω²)/(BW·Ω))^{2m}) — the
    // power-complement of the band-stop, the reference the variable-order bandpassCascade nulls against.
    inline double bandpassM (double f, double f0, double Q, int m)
    {
        const double Om0 = f0, BW = f0 / Q, Om = f;
        if (Om <= 0.0) return 0.0;
        const double invx = (Om0 * Om0 - Om * Om) / (BW * Om);
        return 1.0 / std::sqrt (1.0 + std::pow (std::fabs (invx), 2.0 * (double) m));
    }
    inline double bandstop (double f, double f0, double Q, int m)
    {
        const double Om0 = f0, BW = f0 / Q, Om = f;
        const double x = BW * Om / (Om0 * Om0 - Om * Om);
        return 1.0 / std::sqrt (1.0 + std::pow (std::fabs (x), 2.0 * m));
    }
}

//==============================================================================================
// Digital side: designBand product response + stability.
static double digitalMag (const BandDesign& d, double w)
{
    std::complex<double> h { 1.0, 0.0 };
    for (int s = 0; s < d.n; ++s) h *= evalCoeffs (d.sec[s], w);
    return std::abs (h);
}
static bool designOk (const BandDesign& d)
{
    for (int s = 0; s < d.n; ++s)
    {
        const auto& c = d.sec[s];
        if (! (std::isfinite (c.b0) && std::isfinite (c.b1) && std::isfinite (c.b2)
            && std::isfinite (c.a1) && std::isfinite (c.a2))) return false;
        if (! c.isStable()) return false;
    }
    return true;
}
static double dB (double x) { return 20.0 * std::log10 (std::max (1e-300, x)); }

//==============================================================================================
// Heatmap machinery.
struct CellResult { double err = 0.0; bool broken = false; };
struct MapSpec
{
    std::string name;                                  // file/name stem, e.g. "bell_p12"
    std::string title;                                 // drawn on the image
    std::string yLabel;                                // "Q" / "ORDER" / "GAIN DB"
    int nx = 120, ny = 80;
    double yLo = 0.05, yHi = 40.0; bool yLog = true;   // y axis (Q default)
    // cell evaluator: returns worst error (dB) incl. invariant penalties; sets broken on bad designs.
    // f0 in Hz, y = the y-axis value (Q, order, or gain dB).
    CellResult (*eval) (double f0, double y, double fs, double fCap);
};

static bool gAnyBroken = false;
static bool gWriteFailed = false;

// Shared error scan: |digital − ref| dB over the region + feature points near f0. Both sides are
// FLOOR-CLAMPED at −50 dB first: a −60 dB skirt landing at −85 dB is silence-vs-silence, not a 25 dB
// "error" — but digital energy where the reference says silence (or vice versa) still shows fully.
static double kFloorDb = -50.0;
template <class RefFn>
static double scanError (const BandDesign& d, double fs, double f0, double fCap, bool notchNull, RefFn refMag)
{
    double worst = 0.0;
    auto probe = [&] (double f)
    {
        if (f < 20.0 || f > fCap) return;
        if (notchNull && std::fabs (f - f0) < 0.06 * f0) return;          // null neighbourhood -> invariant covers it
        const double rDb = std::max (kFloorDb, dB (refMag (f)));
        const double gDb = std::max (kFloorDb, dB (digitalMag (d, 2.0 * kPi * f / fs)));
        worst = std::max (worst, std::fabs (gDb - rDb));
    };
    for (int i = 0; i < 300; ++i)                                         // global log grid
        probe (20.0 * std::pow (fCap / 20.0, i / 299.0));
    for (int i = -12; i <= 12; ++i)                                       // feature cluster: ±1/3 oct around f0
        probe (f0 * std::pow (2.0, i / 36.0));
    for (double m : { 0.94, 0.9401, 1.0601, 1.062 }) probe (f0 * m);      // null-exclusion edges
    return worst;
}

//==============================================================================================
// Per-type cell evaluators (invariant penalties baked in).
static double gGain = 12.0;                                               // gain slice for bell/shelf maps

static CellResult evalBell (double f0, double Q, double fs, double fCap)
{
    BandParams p; p.on = true; p.type = FilterType::Bell; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).gainDb = gGain;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    const double A = std::pow (10.0, gGain / 20.0);
    double e = scanError (d, fs, f0, fCap, false, [&] (double f) { return ref::bellMag (f, f0, Q, A); });
    e = std::max (e, std::fabs (dB (digitalMag (d, 2.0 * kPi * f0 / fs)) - gGain));   // centre == gain
    return { e, false };
}
static CellResult evalShelf (double f0, double Q, double fs, double fCap, bool high)
{
    BandParams p; p.on = true; p.type = high ? FilterType::HighShelf : FilterType::LowShelf;
    p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).gainDb = gGain;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    const double A = std::pow (10.0, gGain / 20.0);

    // The resonant shelf's shipped CONTRACT (post feasibility/envelope guards) is: stay inside the
    // analog prototype's ENVELOPE, plateaus exact — the design intentionally degenerates to the
    // smooth Butterworth shelf where the fit turns infeasible, and it intentionally PINS the full
    // plateau at Nyquist (the matched no-cramping feature) where the analog is still mid-transition.
    // So the map measures ESCAPE from the analog envelope over the FULL analog curve — asymptotic
    // plateaus included — not pointwise tracking (the curve sheets show pointwise). A balloon /
    // spurious peak / spurious notch escapes the envelope -> red. Envelope and design are sampled
    // on the SAME frequency set (+ a dense cluster at f0 for the narrow high-Q corner peak).
    std::vector<double> freqs;
    for (int i = 0; i < 300; ++i) freqs.push_back (20.0 * std::pow (fCap / 20.0, i / 299.0));
    for (int i = -16; i <= 16; ++i) freqs.push_back (f0 * std::pow (2.0, i / 48.0));
    double refLo = std::min (0.0, gGain), refHi = std::max (0.0, gGain);   // both plateau asymptotes
    for (double f : freqs)
    {
        const double r = dB (ref::shelfQ (f, f0, Q, A, high));             // full analog curve (any f,
        refLo = std::min (refLo, r); refHi = std::max (refHi, r);          // even beyond the window)
    }
    refLo = std::max (refLo, kFloorDb);
    double e = 0.0;
    for (double f : freqs)
    {
        if (f < 20.0 || f > fCap) continue;
        const double g = std::max (kFloorDb, dB (digitalMag (d, 2.0 * kPi * f / fs)));
        // Overshoot above the envelope is STRICT (a band-shelf must never boost past its contract);
        // under-reach gets a 6 dB allowance — the guarded fallback legitimately under-delivers near
        // the band edge, but a shelf that silently stops shelving altogether still flags.
        e = std::max (e, std::max (g - refHi, (refLo - g) - 6.0));
    }
    // invariant: the FIXED plateau end must be exact (high shelf: DC = 0 dB; low shelf: Nyquist = 0 dB).
    e = std::max (e, std::fabs (dB (digitalMag (d, high ? 1e-4 : kPi)) - 0.0));
    return { std::max (e, 0.0), false };
}
static CellResult evalHiShelf (double f0, double Q, double fs, double fCap) { return evalShelf (f0, Q, fs, fCap, true); }
static CellResult evalLoShelf (double f0, double Q, double fs, double fCap) { return evalShelf (f0, Q, fs, fCap, false); }

static CellResult evalCutSlope (double f0, double orderD, double fs, double fCap, bool hp)
{
    const int order = std::clamp ((int) std::lround (orderD), 1, 16);
    BandParams p; p.on = true; p.type = hp ? FilterType::HighPass : FilterType::LowPass;
    p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = 0.707; p.lane (Lane::Stereo).slope = order * 6;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    auto rf = [&] (double f) { return hp ? ref::butterHP (f, f0, order) : ref::butterLP (f, f0, order); };
    double e = scanError (d, fs, f0, fCap, false, rf);
    e = std::max (e, std::fabs (dB (digitalMag (d, 2.0 * kPi * f0 / fs)) + 3.0103));  // -3.01 dB at fc, ANY order
    return { e, false };
}
static CellResult evalHP (double f0, double y, double fs, double fCap) { return evalCutSlope (f0, y, fs, fCap, true); }
static CellResult evalLP (double f0, double y, double fs, double fCap) { return evalCutSlope (f0, y, fs, fCap, false); }

static CellResult evalBandpass (double f0, double Q, double fs, double fCap)
{
    BandParams p; p.on = true; p.type = FilterType::BandPass; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    double e = scanError (d, fs, f0, fCap, false, [&] (double f) { return ref::bandpass (f, f0, Q); });
    e = std::max (e, std::fabs (dB (digitalMag (d, 2.0 * kPi * f0 / fs))));           // unity at centre, exact
    return { e, false };
}
static CellResult evalBandpassM (double f0, double Q, double fs, double fCap, int slope)
{
    BandParams p; p.on = true; p.type = FilterType::BandPass; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).slope = slope;
    const int m = std::clamp ((std::clamp (slope / 6, 1, 16) + 1) / 2, 1, 8);
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    double e = scanError (d, fs, f0, fCap, false, [&] (double f) { return ref::bandpassM (f, f0, Q, m); });
    e = std::max (e, std::fabs (dB (digitalMag (d, 2.0 * kPi * f0 / fs))));           // unity at centre, exact
    return { e, false };
}
static CellResult evalBandpassOrder (double f0, double orderD, double fs, double fCap)
{ return evalBandpassM (f0, 0.707, fs, fCap, std::clamp ((int) std::lround (orderD), 1, 16) * 6); }
static CellResult evalNotchM (double f0, double Q, double fs, double fCap, int slope)
{
    BandParams p; p.on = true; p.type = FilterType::Notch; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).slope = slope;
    const int m = std::clamp ((std::clamp (slope / 6, 1, 16) + 1) / 2, 1, 8);
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    double e = scanError (d, fs, f0, fCap, true, [&] (double f) { return ref::bandstop (f, f0, Q, m); });
    const double nullDb = dB (digitalMag (d, 2.0 * kPi * f0 / fs));
    if (nullDb > -60.0) e = std::max (e, nullDb + 60.0);                              // null must stay >= 60 dB deep
    e = std::max (e, std::fabs (dB (digitalMag (d, 1e-4))));                          // DC unity
    return { e, false };
}
static CellResult evalNotch1 (double f0, double Q, double fs, double fCap)  { return evalNotchM (f0, Q, fs, fCap, 12); }
static CellResult evalNotchOrder (double f0, double orderD, double fs, double fCap)
{ return evalNotchM (f0, 0.707, fs, fCap, std::clamp ((int) std::lround (orderD), 1, 16) * 6); }

static CellResult evalAllpass (double f0, double Q, double fs, double fCap)
{
    BandParams p; p.on = true; p.type = FilterType::AllPass; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    // |H| must be identically 1; group delay must stay positive (a wrong-phase allpass hides in |H|).
    double e = 0.0, prevPh = 0.0; bool first = true; double minGd = 1e9;
    for (int i = 0; i < 800; ++i)                                        // dense: keeps the unwrap step
    {                                                                    // ~1.4 rad << pi even at Q=40
        const double f = 20.0 * std::pow (fCap / 20.0, i / 799.0), w = 2.0 * kPi * f / fs;
        std::complex<double> h { 1.0, 0.0 };
        for (int s = 0; s < d.n; ++s) h *= evalCoeffs (d.sec[s], w);
        e = std::max (e, std::fabs (dB (std::abs (h))));
        double ph = std::arg (h);
        if (! first)
        {
            double dph = ph - prevPh;
            while (dph > kPi) dph -= 2.0 * kPi;
            while (dph < -kPi) dph += 2.0 * kPi;
            const double dw = w * (1.0 - std::pow (fCap / 20.0, -1.0 / 799.0));
            if (dw > 0) minGd = std::min (minGd, -dph / dw);
        }
        prevPh = ph; first = false;
    }
    if (minGd < -1e-6) e = std::max (e, 10.0);                                        // negative GD -> red
    return { e, false };
}
static CellResult evalTilt (double f0, double gainDb, double fs, double fCap)
{
    BandParams p; p.on = true; p.type = FilterType::Tilt; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).gainDb = gainDb; p.lane (Lane::Stereo).Q = 1.0;
    const BandDesign d = designBand (p, fs);
    if (! designOk (d)) { gAnyBroken = true; return { 0.0, true }; }
    const double G = std::pow (10.0, gainDb / 20.0);
    const double e = scanError (d, fs, f0, fCap, false, [&] (double f) { return ref::tilt (f, f0, G); });
    return { e, false };
}

//==============================================================================================
// Heatmap renderer.
static RGB bandColour (double errDb)
{
    if (errDb < 0.1) return kGreen;
    if (errDb < 1.0) return kYellow;
    if (errDb < 3.0) return kOrange;
    return kRed;
}
static void renderHeatmap (const MapSpec& spec, double fs, double fCap, const std::string& outDir,
                           std::FILE* csv, std::vector<std::string>& statsLines, std::vector<std::string>& imageBlocks)
{
    const int cell = 5, mL = 64, mT = 34, mB = 40, mR = 120;
    const int W = mL + spec.nx * cell + mR, H = mT + spec.ny * cell + mB;
    Canvas cv (W, H);
    const double f0Lo = 10.0, f0Hi = 0.49 * fs;

    int nGreen = 0, nBroken = 0; double worst = 0.0; std::string worstAt;
    for (int iy = 0; iy < spec.ny; ++iy)
        for (int ix = 0; ix < spec.nx; ++ix)
        {
            const double f0 = f0Lo * std::pow (f0Hi / f0Lo, ix / (double) (spec.nx - 1));
            const double y  = spec.yLog ? spec.yLo * std::pow (spec.yHi / spec.yLo, iy / (double) (spec.ny - 1))
                                        : spec.yLo + (spec.yHi - spec.yLo) * iy / (double) (spec.ny - 1);
            const CellResult r = spec.eval (f0, y, fs, fCap);
            RGB c = r.broken ? kBlack : bandColour (r.err);
            if (r.broken) ++nBroken; else if (r.err < 0.1) ++nGreen;
            if (! r.broken && r.err > worst)
            {
                worst = r.err;
                char buf[96]; std::snprintf (buf, sizeof (buf), "f0=%.0f %s=%.3g", f0, spec.yLabel.c_str(), y);
                worstAt = buf;
            }
            const int px = mL + ix * cell, py = mT + (spec.ny - 1 - iy) * cell;
            cv.fill (px, py, px + cell - 1, py + cell - 1, c);
        }

    // frame + axes
    cv.line (mL - 1, mT - 1, mL + spec.nx * cell, mT - 1, kDarkGrey);
    cv.line (mL - 1, mT + spec.ny * cell, mL + spec.nx * cell, mT + spec.ny * cell, kDarkGrey);
    cv.line (mL - 1, mT - 1, mL - 1, mT + spec.ny * cell, kDarkGrey);
    cv.line (mL + spec.nx * cell, mT - 1, mL + spec.nx * cell, mT + spec.ny * cell, kDarkGrey);
    for (double f : { 10.0, 100.0, 1000.0, 10000.0 })
    {
        if (f < f0Lo || f > f0Hi) continue;
        const int x = mL + (int) std::lround (std::log (f / f0Lo) / std::log (f0Hi / f0Lo) * (spec.nx * cell - 1));
        cv.line (x, mT + spec.ny * cell, x, mT + spec.ny * cell + 4, kBlack);
        const std::string lbl = f >= 1000.0 ? std::to_string ((int) (f / 1000)) + "K" : std::to_string ((int) f);
        cv.text (x - 3 * (int) lbl.size(), mT + spec.ny * cell + 8, lbl, kBlack);
    }
    auto yTicks = spec.yLog ? std::vector<double> { 0.1, 1.0, 10.0 }
                            : std::vector<double> { spec.yLo, 0.0, spec.yHi };
    for (double yv : yTicks)
    {
        if (yv < spec.yLo || yv > spec.yHi) continue;
        const double t = spec.yLog ? std::log (yv / spec.yLo) / std::log (spec.yHi / spec.yLo)
                                   : (yv - spec.yLo) / (spec.yHi - spec.yLo);
        const int y = mT + spec.ny * cell - 1 - (int) std::lround (t * (spec.ny * cell - 1));
        cv.line (mL - 5, y, mL - 1, y, kBlack);
        char buf[16]; std::snprintf (buf, sizeof (buf), "%g", yv);
        cv.text (4, y - 3, buf, kBlack);
    }
    cv.text (mL, 8, spec.title, kBlack, 2);
    cv.text (mL + spec.nx * cell / 2 - 30, H - 14, "F0 HZ", kBlack);
    cv.text (4, mT - 12, spec.yLabel, kBlack);

    // legend
    const int lx = mL + spec.nx * cell + 12; int ly = mT + 4;
    auto leg = [&] (RGB c, const char* t) { cv.fill (lx, ly, lx + 10, ly + 10, c); cv.text (lx + 14, ly + 2, t, kBlack); ly += 16; };
    leg (kGreen,  "<0.1 DB"); leg (kYellow, "<1 DB"); leg (kOrange, "<3 DB"); leg (kRed, ">3 DB"); leg (kBlack, "BROKEN");

    const std::string file = spec.name + "_fs" + std::to_string ((int) fs) + ".png";
    if (! png::write (outDir + "/" + file, W, H, cv.px)) gWriteFailed = true;
    const double pctGreen = 100.0 * nGreen / (double) (spec.nx * spec.ny);
    std::fprintf (csv, "%s,%d,%.1f,%.3f,\"%s\",%d\n", spec.name.c_str(), (int) fs, pctGreen, worst, worstAt.c_str(), nBroken);
    char line[256];
    std::snprintf (line, sizeof (line), "| %s | %d | %.1f%% | %.3f dB (%s) | %d |",
                   spec.title.c_str(), (int) fs, pctGreen, worst, worstAt.c_str(), nBroken);
    statsLines.push_back (line);
    imageBlocks.push_back ("### " + spec.title + " — " + std::to_string ((int) fs) + " HZ\n\n![" + file + "](" + file + ")\n");
    std::printf ("  %-28s fs=%-6d green=%5.1f%%  worst=%7.3f dB  broken=%d\n",
                 spec.name.c_str(), (int) fs, pctGreen, worst, nBroken);
}

//==============================================================================================
// Measured float audio path: adaptive-length impulse through EqBand -> own radix-2 FFT -> |H|.
static void fftRadix2 (std::vector<std::complex<double>>& a, bool inverse)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i)
    {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap (a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1)
    {
        const double ang = 2.0 * kPi / (double) len * (inverse ? 1.0 : -1.0);
        const std::complex<double> wl (std::cos (ang), std::sin (ang));
        for (size_t i = 0; i < n; i += len)
        {
            std::complex<double> w (1.0, 0.0);
            for (size_t k = 0; k < len / 2; ++k)
            {
                const auto u = a[i + k], v = a[i + k + len / 2] * w;
                a[i + k] = u + v; a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
}
// Returns |H| sampled at the FFT bin grid (size n/2+1) or empty if the IR never decayed (flagged).
static std::vector<double> measuredMag (const BandParams& p, double fs, size_t& nOut, bool& converged)
{
    size_t n = 1 << 16;
    std::vector<float> h;
    converged = false;
    for (; n <= (1u << 22); n <<= 2)
    {
        EqBand band; band.prepare (fs, 1);
        band.setParams (p);                                               // first set snaps
        h.assign (n, 0.0f); h[0] = 1.0f;
        for (size_t off = 0; off < n; off += 4096)
        {
            float* ch[1] = { h.data() + off };
            band.processBlock (ch, 1, (int) std::min<size_t> (4096, n - off));
        }
        float peak = 0.0f, tail = 0.0f;
        for (size_t i = 0; i < n; ++i) peak = std::max (peak, std::fabs (h[i]));
        for (size_t i = n - 4096; i < n; ++i) tail = std::max (tail, std::fabs (h[i]));
        if (tail < 1e-5f * std::max (1e-30f, peak)) { converged = true; break; }
    }
    nOut = std::min<size_t> (h.size(), 1u << 22);
    std::vector<std::complex<double>> a (nOut);
    for (size_t i = 0; i < nOut; ++i) a[i] = (double) h[i];
    fftRadix2 (a, false);
    std::vector<double> mag (nOut / 2 + 1);
    for (size_t i = 0; i < mag.size(); ++i) mag[i] = std::abs (a[i]);
    return mag;
}

//==============================================================================================
// Curve sheets: digital vs analog vs RBJ + measured float path + error track.
struct CurveCfg
{
    std::string label;
    BandParams p;
    double (*refMag) (double f, const BandParams& p);                    // analog reference
    bool rbj = false;                                                     // draw the RBJ baseline (bell only here)
    bool notchNull = false;
};

static void renderCurveSheet (const std::string& name, const std::string& title,
                              const std::vector<CurveCfg>& cfgs, double fs,
                              const std::string& outDir, std::vector<std::string>& imageBlocks)
{
    const int cols = 2, rows = ((int) cfgs.size() + 1) / 2;
    const int PW = 420, PH = 240, mL = 40, gap = 18, top = 30;
    const int W = mL + cols * (PW + gap), H = top + rows * (PH + 58 + gap);
    Canvas cv (W, H);
    cv.text (mL, 8, title, kBlack, 2);

    for (size_t ci = 0; ci < cfgs.size(); ++ci)
    {
        const auto& cfg = cfgs[ci];
        const int cx = mL + (int) (ci % cols) * (PW + gap);
        const int cy = top + (int) (ci / cols) * (PH + 58 + gap);
        const double fLo = 10.0, fHi = 0.4999 * fs;
        const double dbLo = -60.0, dbHi = 36.0;

        cv.fill (cx, cy, cx + PW, cy + PH, kWhite);
        for (double f : { 100.0, 1000.0, 10000.0 })
        {
            const int x = cx + (int) (std::log (f / fLo) / std::log (fHi / fLo) * PW);
            cv.line (x, cy, x, cy + PH, kGrey);
        }
        for (double db : { -48.0, -24.0, -12.0, 0.0, 12.0, 24.0 })
        {
            const int y = cy + (int) ((dbHi - db) / (dbHi - dbLo) * PH);
            cv.line (cx, y, cx + PW, y, db == 0.0 ? kDarkGrey : kGrey);
        }

        const BandDesign d = designBand (cfg.p, fs);
        auto X = [&] (double f)  { return cx + (int) (std::log (f / fLo) / std::log (fHi / fLo) * PW); };
        auto Y = [&] (double db) { return cy + (int) ((dbHi - std::clamp (db, dbLo, dbHi)) / (dbHi - dbLo) * PH); };

        // measured float path (drawn first, thick-ish, behind the analytic line)
        size_t nFft = 0; bool conv = false;
        const std::vector<double> meas = measuredMag (cfg.p, fs, nFft, conv);
        int lx = -1, ly = -1;
        for (int i = 0; i < 480; ++i)
        {
            const double f = fLo * std::pow (fHi / fLo, i / 479.0);
            const size_t bin = std::min (meas.size() - 1, (size_t) std::llround (f / fs * (double) nFft));
            const int x = X (f), y = Y (dB (meas[bin]));
            if (lx >= 0) { cv.line (lx, ly, x, y, kMeas); cv.line (lx, ly + 1, x, y + 1, kMeas); }
            lx = x; ly = y;
        }
        // analog reference (dashed)
        lx = -1; ly = -1;
        for (int i = 0; i < 480; ++i)
        {
            const double f = fLo * std::pow (fHi / fLo, i / 479.0);
            const int x = X (f), y = Y (dB (cfg.refMag (f, cfg.p)));
            if (lx >= 0 && (i / 4) % 2 == 0) cv.line (lx, ly, x, y, kAnalog);
            lx = x; ly = y;
        }
        // RBJ baseline (bell): shows the cookbook cramping the matched design avoids
        if (cfg.rbj)
        {
            const auto r = rbj::peaking (cfg.p.lane (Lane::Stereo).freq, fs, cfg.p.lane (Lane::Stereo).Q, std::pow (10.0, cfg.p.lane (Lane::Stereo).gainDb / 20.0));
            lx = -1; ly = -1;
            for (int i = 0; i < 480; ++i)
            {
                const double f = fLo * std::pow (fHi / fLo, i / 479.0);
                const int x = X (f), y = Y (r.magnitudeDb (2.0 * kPi * f / fs));
                if (lx >= 0) cv.line (lx, ly, x, y, kRbj);
                lx = x; ly = y;
            }
        }
        // digital analytic (the display/GUI truth)
        lx = -1; ly = -1;
        for (int i = 0; i < 480; ++i)
        {
            const double f = fLo * std::pow (fHi / fLo, i / 479.0);
            const int x = X (f), y = Y (dB (digitalMag (d, 2.0 * kPi * f / fs)));
            if (lx >= 0) cv.line (lx, ly, x, y, kBlue);
            lx = x; ly = y;
        }

        // error track: |analytic − ref| (orange) and |measured − analytic| (magenta)
        const int ey0 = cy + PH + 6, EH = 34;
        cv.fill (cx, ey0, cx + PW, ey0 + EH, {245,245,247});
        cv.line (cx, ey0 + EH, cx + PW, ey0 + EH, kDarkGrey);
        auto EY = [&] (double e) { return ey0 + EH - (int) (std::clamp (e, 0.0, 3.0) / 3.0 * EH); };
        int l1x=-1,l1y=-1,l2x=-1,l2y=-1;
        for (int i = 0; i < 480; ++i)
        {
            const double f = fLo * std::pow (fHi / fLo, i / 479.0);
            const bool skip = cfg.notchNull && std::fabs (f - cfg.p.lane (Lane::Stereo).freq) < 0.06 * cfg.p.lane (Lane::Stereo).freq;
            const double rm = cfg.refMag (f, cfg.p);
            const double am = digitalMag (d, 2.0 * kPi * f / fs);
            const size_t bin = std::min (meas.size() - 1, (size_t) std::llround (f / fs * (double) nFft));
            const int x = X (f);
            if (! skip && rm > 3.16e-3)
            {
                const int y = EY (std::fabs (dB (am) - dB (rm)));
                if (l1x >= 0) cv.line (l1x, l1y, x, y, kOrange);
                l1x = x; l1y = y;
            } else l1x = -1;
            if (! skip && am > 3.16e-3)
            {
                const int y = EY (std::fabs (dB (meas[bin]) - dB (am)));
                if (l2x >= 0) cv.line (l2x, l2y, x, y, kMeas);
                l2x = x; l2y = y;
            } else l2x = -1;
        }
        cv.text (cx, ey0 + EH + 4, cfg.label + (conv ? "" : "  (IR NOT DECAYED)"), kBlack);
    }
    // legend
    const int lgy = H - 12;
    cv.fill (mL, lgy, mL + 16, lgy + 2, kBlue);        cv.text (mL + 20, lgy - 3, "DIGITAL", kBlack);
    cv.fill (mL + 90, lgy, mL + 106, lgy + 2, kAnalog); cv.text (mL + 110, lgy - 3, "ANALOG", kBlack);
    cv.fill (mL + 180, lgy, mL + 196, lgy + 2, kMeas);  cv.text (mL + 200, lgy - 3, "MEASURED FLOAT", kBlack);
    cv.fill (mL + 310, lgy, mL + 326, lgy + 2, kRbj);   cv.text (mL + 330, lgy - 3, "RBJ", kBlack);

    const std::string file = name + "_fs" + std::to_string ((int) fs) + ".png";
    if (! png::write (outDir + "/" + file, W, H, cv.px)) gWriteFailed = true;
    imageBlocks.push_back ("### CURVES: " + title + " — " + std::to_string ((int) fs) + " HZ\n\n![" + file + "](" + file + ")\n");
    std::printf ("  %-28s fs=%-6d curve sheet\n", name.c_str(), (int) fs);
}

//==============================================================================================
// Curve-sheet reference adapters.
static double refBellP    (double f, const BandParams& p) { return ref::bellMag  (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q, std::pow (10.0, p.lane (Lane::Stereo).gainDb / 20.0)); }
static double refLoShelfP (double f, const BandParams& p) { return ref::shelfQ   (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q, std::pow (10.0, p.lane (Lane::Stereo).gainDb / 20.0), false); }
static double refHiShelfP (double f, const BandParams& p) { return ref::shelfQ   (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q, std::pow (10.0, p.lane (Lane::Stereo).gainDb / 20.0), true); }
static double refHPp      (double f, const BandParams& p) { return ref::butterHP (f, p.lane (Lane::Stereo).freq, std::clamp (p.lane (Lane::Stereo).slope / 6, 1, 16)); }
static double refLPp      (double f, const BandParams& p) { return ref::butterLP (f, p.lane (Lane::Stereo).freq, std::clamp (p.lane (Lane::Stereo).slope / 6, 1, 16)); }
static double refBPp      (double f, const BandParams& p) { return ref::bandpass (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q); }
static double refBPcasP    (double f, const BandParams& p)
{ return ref::bandpassM (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q, std::clamp ((std::clamp (p.lane (Lane::Stereo).slope / 6, 1, 16) + 1) / 2, 1, 8)); }
static double refNotchP   (double f, const BandParams& p)
{ return ref::bandstop (f, p.lane (Lane::Stereo).freq, p.lane (Lane::Stereo).Q, std::clamp ((std::clamp (p.lane (Lane::Stereo).slope / 6, 1, 16) + 1) / 2, 1, 8)); }
static double refAPp      (double f, const BandParams&)   { (void) f; return 1.0; }
static double refTiltP    (double f, const BandParams& p) { return ref::tilt     (f, p.lane (Lane::Stereo).freq, std::pow (10.0, p.lane (Lane::Stereo).gainDb / 20.0)); }

static BandParams mk (FilterType t, double f0, double Q, double g = 0.0, int slope = 12, bool swept = false)
{
    BandParams p; p.on = true; p.type = t; p.lane (Lane::Stereo).freq = f0; p.lane (Lane::Stereo).Q = Q; p.lane (Lane::Stereo).gainDb = g; p.lane (Lane::Stereo).slope = slope; p.swept = swept;
    return p;
}

//==============================================================================================
int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf (stderr, "usage: %s <outDir> [--fs 44100,96000] [--full]\n", argv[0]);
        return 2;
    }
    const std::string outDir = argv[1];
    std::vector<double> rates { 44100.0, 96000.0 };
    bool full = false;
    for (int i = 2; i < argc; ++i)
    {
        if (! std::strcmp (argv[i], "--full")) full = true;
        else if (! std::strcmp (argv[i], "--fs") && i + 1 < argc)
        {
            rates.clear();
            for (char* tok = std::strtok (argv[++i], ","); tok; tok = std::strtok (nullptr, ","))
                rates.push_back (std::atof (tok));
        }
    }

    std::FILE* csv = std::fopen ((outDir + "/summary.csv").c_str(), "w");
    if (! csv) { std::fprintf (stderr, "cannot write to %s (does the directory exist?)\n", outDir.c_str()); return 2; }
    std::fprintf (csv, "map,fs,pct_green,worst_db,worst_at,broken_cells\n");
    std::vector<std::string> statsLines;
    std::vector<std::string> imageBlocks;

    std::vector<MapSpec> maps;
    auto add = [&] (const char* n, const char* t, const char* yl, double lo, double hi, bool ylog,
                    CellResult (*fn) (double, double, double, double))
    { maps.push_back ({ n, t, yl, 120, 80, lo, hi, ylog, fn }); };

    for (double fs : rates)
    {
        const double fCap = (full ? 0.4999 : 0.95 * 0.5) * fs;
        std::printf ("== fs = %.0f  (error region up to %.0f Hz)%s ==\n", fs, fCap, full ? "  [FULL]" : "");

        // gain-sliced maps (bell/shelves): +12 and -12 dB
        for (double g : { 12.0, -12.0 })
        {
            gGain = g;
            const std::string sfx = g > 0 ? "_p12" : "_m12";
            const std::string gl  = g > 0 ? " +12DB" : " -12DB";
            maps.clear();
            add (("bell"    + sfx).c_str(), ("BELL"       + gl).c_str(), "Q", 0.05, 40.0, true, evalBell);
            add (("loshelf" + sfx).c_str(), ("LOW SHELF"  + gl).c_str(), "Q", 0.05, 40.0, true, evalLoShelf);
            add (("hishelf" + sfx).c_str(), ("HIGH SHELF" + gl).c_str(), "Q", 0.05, 40.0, true, evalHiShelf);
            for (const auto& m : maps) renderHeatmap (m, fs, fCap, outDir, csv, statsLines, imageBlocks);
        }
        // Q-independent / other-axis maps
        maps.clear();
        add ("bandpass", "BANDPASS",             "Q",       0.05, 40.0, true,  evalBandpass);
        add ("notch1",   "NOTCH SINGLE",         "Q",       0.05, 40.0, true,  evalNotch1);
        add ("allpass",  "ALLPASS (MAG+GD)",     "Q",       0.05, 40.0, true,  evalAllpass);
        add ("hp",       "HIGHPASS BUTTERWORTH", "ORDER",   1.0,  16.0, false, evalHP);
        add ("lp",       "LOWPASS BUTTERWORTH",  "ORDER",   1.0,  16.0, false, evalLP);
        add ("notchcas", "NOTCH CASCADE Q=.707", "ORDER",   1.0,  16.0, false, evalNotchOrder);
        add ("bpcascade","BANDPASS CASCADE Q=.707","ORDER", 1.0,  16.0, false, evalBandpassOrder);
        add ("tilt",     "TILT",                 "GAIN DB", -30.0, 30.0, false, evalTilt);
        for (const auto& m : maps) renderHeatmap (m, fs, fCap, outDir, csv, statsLines, imageBlocks);

        // curve sheets with the measured float path (sentinel + vanilla configs per type)
        renderCurveSheet ("curves_bell_shelf", "BELL + SHELVES", {
            { "BELL 1K Q2 +12",        mk (FilterType::Bell,      1000, 2.0,  12), refBellP,    true  },
            { "BELL 16K Q4 +6 (AIR)",  mk (FilterType::Bell,      16000, 4.0,  6), refBellP,    true  },
            { "HISHELF 18K Q0.7 +6",   mk (FilterType::HighShelf, 18000, 0.7,  6), refHiShelfP, false },
            { "LOSHELF 0.45FS Q0.05 +6", mk (FilterType::LowShelf, 0.45 * fs, 0.05, 6), refLoShelfP, false },
        }, fs, outDir, imageBlocks);
        renderCurveSheet ("curves_cuts_bp", "CUTS + BANDPASS", {
            { "HP 1K ORDER 3 (18DB/OCT)", mk (FilterType::HighPass, 1000, 0.707, 0, 18), refHPp, false },
            { "LP 4K ORDER 8",            mk (FilterType::LowPass,  4000, 0.707, 0, 48), refLPp, false },
            { "BP 10.6HZ Q40 (LOW CORNER)", mk (FilterType::BandPass, 10.61, 40.0, 0),   refBPp, false },
            { "BP CASCADE 1K Q2 M4 (48DB/OCT)", mk (FilterType::BandPass, 1000, 2.0, 0, 48), refBPcasP, false },
        }, fs, outDir, imageBlocks);
        renderCurveSheet ("curves_notch_ap_tilt", "NOTCH + ALLPASS + TILT", {
            { "NOTCH 15K Q0.5 M4 (ALIAS)", mk (FilterType::Notch, 15000, 0.5, 0, 48), refNotchP, false, true },
            { "NOTCH 1K Q3 M8",            mk (FilterType::Notch, 1000, 3.0, 0, 96),  refNotchP, false, true },
            { "ALLPASS 2K Q2",             mk (FilterType::AllPass, 2000, 2.0, 0),    refAPp,    false },
            { "TILT 1K +12 (SWEPT FLAG)",  mk (FilterType::Tilt, 1000, 1.0, 12, 12, true), refTiltP, false },
        }, fs, outDir, imageBlocks);
    }

    // index.md
    if (std::FILE* idx = std::fopen ((outDir + "/index.md").c_str(), "w"))
    {
        std::fprintf (idx, "# EQ Filter Quality Atlas\n\nEvery `felitronics::eq` filter type measured against an INDEPENDENT analog reference and rendered as quality maps. Generated by `tools/eqviz`.\n\nColour bands: green <0.1 dB · yellow <1 dB · orange <3 dB · red >=3 dB · black = unstable/non-finite.\n");
        std::fprintf (idx, "Error = worst |digital - analog reference| over f in [20 Hz, %s], both sides floor-clamped at -50 dB, MAXed with per-type invariants (centre gain, plateaus, null depth >= 60 dB, |H(fc)| = -3.01 dB, allpass |H|=1 + group delay > 0). Resonant shelves are scored against the analog ENVELOPE (their shipped contract), not pointwise.\n\n", full ? "0.4999 fs" : "0.95 fs/2");
        std::fprintf (idx, "READ THE MAPS AS DIFFS: some yellow/orange near the top of the band and at the extreme right edge is the matched design's INTENDED near-Nyquist trade (no cramping; first-order sections are bilinear; the single notch's frozen legacy sag; the BANDPASS CASCADE is a bilinear Butterworth band-pass, so its high-f0 columns roll off to the z=-1 Nyquist zero a touch faster than the analog skirt — deep-stopband, inaudible). A regression is NEW red/black, or orange spreading into the mid-band, relative to main's maps.\n\n");
        std::fprintf (idx, "## Summary\n\n| map | fs | green | worst | broken |\n|---|---|---|---|---|\n");
        for (const auto& l : statsLines) std::fprintf (idx, "%s\n", l.c_str());
        std::fprintf (idx, "\n## Maps\n\n");
        for (const auto& b : imageBlocks) std::fprintf (idx, "%s\n", b.c_str());
        std::fclose (idx);
    }
    std::fclose (csv);
    if (gWriteFailed) { std::fprintf (stderr, "PNG write failed.\n"); return 1; }
    if (gAnyBroken) { std::printf ("BROKEN designs encountered — see black cells.\n"); return 1; }
    std::printf ("done.\n");
    return 0;
}
