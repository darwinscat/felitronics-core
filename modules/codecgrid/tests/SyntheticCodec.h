// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// TEST FIXTURES: programmes that carry a codec's frame grid, made here, with no codec and no audio file.
//
// A lossy codec does three things that matter to this module: it transforms a frame, it sets the small
// coefficients to zero, and a decoder transforms back. These fixtures do exactly that and nothing else — no
// psychoacoustics, no bit allocation, no bitstream: the smallest `zeroShare` of each frame's coefficients become
// zero, the rest pass unquantised. That is a coarser thing than a codec and it is enough: the zeros are what the
// scan reads.
//
//   mdctCoded   an MDCT codec with any window satisfying Princen-Bradley (AAC's two windows, CELT's), frames
//               starting at `offset`; with a pre-emphasis coefficient it codes the pre-emphasised signal and
//               de-emphasises the result, as CELT does.
//   mp3Coded    the Layer III hybrid: the analysis is the module's own Mp3Hybrid (nulled against the standard in
//               Mp3HybridTests), the SYNTHESIS is written here from the standard's decoder — inverse alias
//               butterflies, the 18 -> 36 IMDCT with overlap-add, the frequency inversion, the polyphase
//               synthesis with the D window — and is the half of the format the module does not contain.
//
// Each returns the coded signal ALIGNED with its input (sample i of the output is sample i of the input), so the
// grid of a programme coded at `offset` is at `offset`; and each reports the residual of its own round trip with
// nothing zeroed, which the suites assert before they believe a fixture.

#include <felitronics/codecgrid/GridScan.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace synthetic
{
constexpr double kPi = 3.14159265358979323846;

struct Stereo
{
    std::vector<float> left, right;
};

// Programme material: three decorrelated noises of different colour and a few partials, left and right sharing
// most of it and each carrying something of its own.
inline Stereo programme (int n, std::uint32_t seed, double level = 0.25)
{
    Stereo s;
    s.left.resize ((std::size_t) n);
    s.right.resize ((std::size_t) n);
    std::uint32_t r = seed * 2654435761u + 12345u;
    auto white = [&r]() noexcept { r = r * 1664525u + 1013904223u; return (double) (r >> 8) / 8388608.0 - 1.0; };
    double lpC = 0.0, lpL = 0.0, lpR = 0.0, lp2 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double wc = white(), wl = white(), wr = white();
        lpC += 0.02 * (wc - lpC);            // dark
        lp2 += 0.3 * (wc - lp2);             // mid
        lpL += 0.1 * (wl - lpL);
        lpR += 0.1 * (wr - lpR);
        const double tone = 0.5 * std::sin (2.0 * kPi * 0.0113 * i) + 0.3 * std::sin (2.0 * kPi * 0.0931 * i + 0.7) + 0.2 * std::sin (2.0 * kPi * 0.2377 * i + 2.1);
        const double centre = 4.0 * lpC + 0.8 * lp2 + 0.15 * wc + 0.4 * tone;
        s.left[(std::size_t) i]  = (float) (level * (centre + 1.2 * lpL + 0.05 * wl));
        s.right[(std::size_t) i] = (float) (level * (centre + 1.2 * lpR + 0.05 * wr));
    }
    return s;
}

// Zero the smallest `share` of v[0 .. n) by magnitude.
inline void zeroSmallest (float* v, int n, double share, std::vector<float>& scratch)
{
    const int k = (int) std::floor (share * (double) n);
    if (k <= 0) return;
    scratch.resize ((std::size_t) n);
    for (int i = 0; i < n; ++i) scratch[(std::size_t) i] = std::fabs (v[i]);
    std::nth_element (scratch.begin(), scratch.begin() + (k - 1), scratch.end());
    const float threshold = scratch[(std::size_t) (k - 1)];
    for (int i = 0; i < n; ++i) if (std::fabs (v[i]) <= threshold) v[i] = 0.0f;
}

//==============================================================================
// One channel through an MDCT codec.
inline std::vector<float> mdctCodedChannel (const std::vector<float>& x, int frame, const std::vector<double>& window, int offset,
                                            double zeroShare, double preEmphasis)
{
    using felitronics::codecgrid::Mdct;
    const int n = (int) x.size(), m = frame / 2, q = m / 2;
    std::vector<float> in (x);
    if (preEmphasis > 0.0)
        for (int i = n - 1; i >= 1; --i) in[(std::size_t) i] = x[(std::size_t) i] - (float) preEmphasis * x[(std::size_t) (i - 1)];

    Mdct<> analysis, dct4;
    std::vector<double> ones ((std::size_t) frame, 1.0);
    if (! analysis.prepare (frame, window.data()) || ! dct4.prepare (frame, ones.data())) return {};

    std::vector<double> acc ((std::size_t) n, 0.0);
    std::vector<float> buf ((std::size_t) frame), coeff ((std::size_t) m), u ((std::size_t) m), scratch;
    for (int start = offset - 2 * m; start < n; start += m)
    {
        for (int i = 0; i < frame; ++i)
        {
            const int at = start + i;
            buf[(std::size_t) i] = at >= 0 && at < n ? in[(std::size_t) at] : 0.0f;
        }
        analysis.transform (buf.data(), coeff.data());
        zeroSmallest (coeff.data(), m, zeroShare, scratch);

        // the DCT-IV of the coefficients, through the same engine with a rectangular window: a frame
        // (a, b, c, d) = (v2, 0, 0, -v1) folds to (v1, v2)
        for (int i = 0; i < q; ++i)
        {
            buf[(std::size_t) i] = coeff[(std::size_t) (q + i)];
            buf[(std::size_t) (q + i)] = 0.0f;
            buf[(std::size_t) (m + i)] = 0.0f;
            buf[(std::size_t) (m + q + i)] = -coeff[(std::size_t) i];
        }
        dct4.transform (buf.data(), u.data());
        // unfold (the transpose of the fold), window, overlap-add: (u2, -u2 reversed, -u1 reversed, -u1) * 2 / M
        for (int i = 0; i < q; ++i)
        {
            const double a = u[(std::size_t) (q + i)], b = -u[(std::size_t) (m - 1 - i)], c = -u[(std::size_t) (q - 1 - i)], d = -u[(std::size_t) i];
            const double v[4] { a, b, c, d };
            for (int part = 0; part < 4; ++part)
            {
                const int at = start + part * q + i;
                if (at >= 0 && at < n) acc[(std::size_t) at] += 2.0 / m * v[part] * window[(std::size_t) (part * q + i)];
            }
        }
    }
    std::vector<float> y ((std::size_t) n);
    if (preEmphasis > 0.0)
    {
        double prev = 0.0;
        for (int i = 0; i < n; ++i) { prev = acc[(std::size_t) i] + (i > 0 ? preEmphasis * prev : 0.0); y[(std::size_t) i] = (float) prev; }
    }
    else
        for (int i = 0; i < n; ++i) y[(std::size_t) i] = (float) acc[(std::size_t) i];
    return y;
}

inline Stereo mdctCoded (const Stereo& x, felitronics::codecgrid::Transform t, int offset, double zeroShare)
{
    using namespace felitronics::codecgrid;
    const int frame = t == Transform::Celt ? kCeltFrame : kAacFrame;
    std::vector<double> w ((std::size_t) frame), scratch ((std::size_t) frame / 2 + 1);
    if (t == Transform::AacSine) window::sine (w.data(), frame);
    else if (t == Transform::AacKbd) window::kaiserBesselDerived (w.data(), frame, 4.0, scratch.data());
    else window::celt (w.data(), frame, kCeltOverlap);
    const double pre = t == Transform::Celt ? kCeltPreEmphasis : 0.0;
    return { mdctCodedChannel (x.left, frame, w, offset, zeroShare, pre), mdctCodedChannel (x.right, frame, w, offset, zeroShare, pre) };
}

//==============================================================================
// One channel through the Layer III hybrid filterbank and back.
inline std::vector<float> mp3CodedChannel (const std::vector<float>& x, int offset, double zeroShare)
{
    using felitronics::codecgrid::Mp3Hybrid;
    using felitronics::codecgrid::kMp3SynthesisWindow;
    const int n = (int) x.size(), pad = 4 * 576, total = n + 2 * pad;
    const int phase = offset % 32, start = offset / 32;
    std::vector<float> xe ((std::size_t) total, 0.0f);
    std::copy (x.begin(), x.end(), xe.begin() + pad);

    const Mp3Hybrid bank;
    const int count = Mp3Hybrid::subbandSamples (total, phase);
    std::vector<float> sub ((std::size_t) count * 32);
    bank.subbands (xe.data(), total, phase, sub.data());

    // the decoder's half, from the standard: per granule the inverse butterflies and the IMDCT, overlap-added
    // into subband samples; then the frequency inversion undone
    std::vector<double> rec ((std::size_t) count * 32, 0.0);
    static const double ci[8] { -0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037 };
    std::vector<double> cosTable ((std::size_t) 36 * 18), win (36);
    for (int i = 0; i < 36; ++i)
    {
        win[(std::size_t) i] = std::sin (kPi / 36.0 * (i + 0.5));
        for (int k = 0; k < 18; ++k) cosTable[(std::size_t) (i * 18 + k)] = std::cos (kPi / 72.0 * (2.0 * i + 1.0 + 18.0) * (2.0 * k + 1.0));
    }
    std::vector<float> scratch;
    const int granules = Mp3Hybrid::granules (count, start);
    for (int g = 0; g < granules; ++g)
    {
        const int first = start + 18 * g;
        float lines[576];
        bank.granule (sub.data() + (std::size_t) first * 32, lines);
        zeroSmallest (lines, 576, zeroShare, scratch);
        double l[576];
        for (int i = 0; i < 576; ++i) l[i] = lines[i];
        for (int sb = 1; sb < 32; ++sb)
            for (int i = 0; i < 8; ++i)
            {
                const double cs = 1.0 / std::sqrt (1.0 + ci[i] * ci[i]), ca = ci[i] / std::sqrt (1.0 + ci[i] * ci[i]);
                const double lo = l[18 * sb - 1 - i], hi = l[18 * sb + i];
                l[18 * sb - 1 - i] = lo * cs - hi * ca;                     // the transpose of the encoder's rotation
                l[18 * sb + i] = hi * cs + lo * ca;
            }
        for (int sb = 0; sb < 32; ++sb)
            for (int i = 0; i < 36; ++i)
            {
                double acc = 0.0;
                for (int k = 0; k < 18; ++k) acc += l[sb * 18 + k] * cosTable[(std::size_t) (i * 18 + k)];
                rec[(std::size_t) ((first + i) * 32 + sb)] += acc * win[(std::size_t) i] / 9.0;
            }
    }
    for (int t = start; t < count; ++t)
        if ((((t - start) % 18) & 1) != 0)
            for (int sb = 1; sb < 32; sb += 2) rec[(std::size_t) (t * 32 + sb)] = -rec[(std::size_t) (t * 32 + sb)];

    // the polyphase synthesis of the standard's decoder (figure 3-A.2)
    std::vector<double> matrix ((std::size_t) 64 * 32);
    for (int i = 0; i < 64; ++i)
        for (int k = 0; k < 32; ++k) matrix[(std::size_t) (i * 32 + k)] = std::cos ((16.0 + i) * (2.0 * k + 1.0) * kPi / 64.0);
    std::vector<double> v (1024, 0.0), ye ((std::size_t) total, 0.0);
    for (int t = 0; t < count; ++t)
    {
        for (int i = 1023; i >= 64; --i) v[(std::size_t) i] = v[(std::size_t) (i - 64)];
        for (int i = 0; i < 64; ++i)
        {
            double acc = 0.0;
            for (int k = 0; k < 32; ++k) acc += matrix[(std::size_t) (i * 32 + k)] * rec[(std::size_t) (t * 32 + k)];
            v[(std::size_t) i] = acc;
        }
        for (int j = 0; j < 32; ++j)
        {
            double acc = 0.0;
            for (int i = 0; i < 8; ++i)
            {
                acc += v[(std::size_t) (i * 128 + j)] * kMp3SynthesisWindow[(std::size_t) (i * 64 + j)];
                acc += v[(std::size_t) (i * 128 + 96 + j)] * kMp3SynthesisWindow[(std::size_t) (i * 64 + 32 + j)];
            }
            // analysis and synthesis together delay by 481 samples: the 32 samples emitted with subband sample t
            // are the input's samples phase + 32 t - 1 .. phase + 32 t + 30
            const int at = phase + 32 * t + j - 1;
            if (at >= 0 && at < total) ye[(std::size_t) at] = acc;
        }
    }
    std::vector<float> y ((std::size_t) n);
    for (int i = 0; i < n; ++i) y[(std::size_t) i] = (float) ye[(std::size_t) (pad + i)];
    return y;
}

inline Stereo mp3Coded (const Stereo& x, int offset, double zeroShare)
{
    return { mp3CodedChannel (x.left, offset, zeroShare), mp3CodedChannel (x.right, offset, zeroShare) };
}

inline Stereo coded (const Stereo& x, felitronics::codecgrid::Transform t, int offset, double zeroShare)
{
    return t == felitronics::codecgrid::Transform::Mp3 ? mp3Coded (x, offset, zeroShare) : mdctCoded (x, t, offset, zeroShare);
}

// The level of a - b against the level of a, in dB (how deep a round trip nulls). -300 when identical.
inline double residualDb (const std::vector<float>& a, const std::vector<float>& b, int skip)
{
    double num = 0.0, den = 0.0;
    for (std::size_t i = (std::size_t) skip; i + (std::size_t) skip < a.size(); ++i)
    {
        const double d = (double) a[i] - (double) b[i];
        num += d * d;
        den += (double) a[i] * (double) a[i];
    }
    return den > 0.0 && num > 0.0 ? 10.0 * std::log10 (num / den) : -300.0;
}
}
