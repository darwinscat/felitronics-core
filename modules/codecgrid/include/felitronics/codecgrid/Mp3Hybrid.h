// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/codecgrid/Mp3Window.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>

#include <array>
#include <cmath>
#include <cstddef>

//==============================================================================
// felitronics::codecgrid::Mp3Hybrid — the ANALYSIS half of the MPEG-1 Layer III hybrid filterbank, as the
// standard's encoder runs it (ISO/IEC 11172-3, 2.4.3.4 and annex C): the 32-band polyphase analysis, then per
// subband a 36-point MDCT of two consecutive blocks of 18 subband samples (long blocks, the sine window), the
// frequency inversion of odd subbands, and the alias-reduction butterflies. 576 lines per granule.
//
// WHY AN ENCODER'S FILTERBANK LIVES IN A MEASUREMENT MODULE. A decoder's output, analysed again with the
// encoder's own transform FROM THE SAME SAMPLE, gives back the dequantised lines — and the lines the quantiser
// set to zero come back as zeros. One sample off, they are filled by their neighbours. So this class exists to
// be run at every one of the 576 possible granule offsets of a stretch of PCM (codecgrid::GridScan does that);
// it is not a codec and nothing here encodes anything.
//
// WHAT IS AND IS NOT MODELLED. Long blocks only (block type 0): short and mixed blocks use other windows, and a
// granule coded with them does not line up — it reads as a granule without zeros, which costs evidence and
// produces none. No scalefactors, no quantiser, no bitstream.
//
// THE OFFSET. A granule offset o in 0..575 is a polyphase phase o % 32 and a granule start o / 32 in subband
// samples: subband sample t of phase p is computed from x[p + 32 t .. p + 32 t + 511], the newest sample last,
// and granule q of start g takes subband samples g + 18 q .. g + 18 q + 35.
//
// ARITHMETIC. Single precision in the transform, tables designed in double through core::det and rounded once:
// the same tables on every platform. Measured on the research prototype before this was written: the scan gives
// the same verdicts in float32 as in double (score within 0.93-1.02 on identical samples).
//
// The MDCT is computed as a fold of the 36 windowed samples to 18 followed by an 18 x 18 DCT-IV — the same
// numbers as the 36 x 18 matrix the standard prints, at half the multiplications. Mp3HybridTests nulls it
// against that matrix written out directly.
//==============================================================================
namespace felitronics::codecgrid
{

class Mp3Hybrid
{
public:
    static constexpr int kSubbands = 32;
    static constexpr int kBlock = 18;                      // subband samples per granule, and lines per subband
    static constexpr int kLines = kSubbands * kBlock;      // 576
    static constexpr int kOffsets = kLines;                // every granule offset, in samples
    static constexpr int kWindow = kMp3WindowLength;       // 512: the polyphase prototype's length
    static constexpr int kAliasPairs = 8;

    Mp3Hybrid() noexcept
    {
        for (int i = 0; i < kWindow; ++i) c_[(std::size_t) i] = (float) (kMp3SynthesisWindow[(std::size_t) i] / 32.0);
        // the matrixing of the analysis: cos ((2 sb + 1) (j - 16) pi / 64)
        for (int sb = 0; sb < kSubbands; ++sb)
            for (int j = 0; j < 64; ++j)
                poly_[(std::size_t) sb][(std::size_t) j] =
                    (float) core::det::cos ((double) ((2 * sb + 1) * (j - 16)) * core::kPi / 64.0);
        for (int n = 0; n < 2 * kBlock; ++n)
            win_[(std::size_t) n] = (float) core::det::sin (core::kPi / 36.0 * ((double) n + 0.5));
        // the DCT-IV the fold leaves: cos (pi / 18 (n + 1/2) (k + 1/2)), written with integers so that the
        // argument is one rounding of an exact product
        for (int n = 0; n < kBlock; ++n)
            for (int k = 0; k < kBlock; ++k)
                dct_[(std::size_t) n][(std::size_t) k] =
                    (float) core::det::cos ((double) ((2 * n + 1) * (2 * k + 1)) * core::kPi / 72.0);
        // the alias-reduction coefficients, table 3-B.9: cs = 1 / sqrt (1 + c^2), ca = c / sqrt (1 + c^2)
        constexpr double ci[kAliasPairs] { -0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037 };
        for (int i = 0; i < kAliasPairs; ++i)
        {
            const double r = std::sqrt (1.0 + ci[i] * ci[i]);
            cs_[(std::size_t) i] = (float) (1.0 / r);
            ca_[(std::size_t) i] = (float) (ci[i] / r);
        }
    }

    // How many subband samples a stretch of `n` samples yields at polyphase phase `phase` (0..31).
    static constexpr int subbandSamples (int n, int phase) noexcept
    {
        return n - phase < kWindow ? 0 : (n - phase - kWindow) / kSubbands + 1;
    }

    // How many granules `subbandSamples` subband samples yield from granule start `g` (0..17): a granule needs
    // its own block and the one after it.
    static constexpr int granules (int subbandSamples, int g) noexcept
    {
        const int blocks = (subbandSamples - g) / kBlock;
        return blocks < 2 ? 0 : blocks - 1;
    }

    // The polyphase analysis of one channel at one phase. `out` takes subbandSamples (n, phase) * 32 floats:
    // out[t * 32 + sb].
    void subbands (const float* x, int n, int phase, float* out) const noexcept
    {
        const int count = subbandSamples (n, phase);
        for (int t = 0; t < count; ++t)
        {
            // z[i] = C[i] * (the sample i steps back from the newest), folded to 64: y[j] = sum_m z[j + 64 m]
            const float* newest = x + phase + kSubbands * t + (kWindow - 1);
            float y[64];
            for (int j = 0; j < 64; ++j)
            {
                float acc = 0.0f;
                for (int m = 0; m < 8; ++m) acc += newest[-(j + 64 * m)] * c_[(std::size_t) (j + 64 * m)];
                y[j] = acc;
            }
            float* o = out + (std::size_t) t * kSubbands;
            for (int sb = 0; sb < kSubbands; ++sb)
            {
                const float* row = poly_[(std::size_t) sb].data();
                float acc = 0.0f;
                for (int j = 0; j < 64; ++j) acc += y[j] * row[j];
                o[sb] = acc;
            }
        }
    }

    // One granule: `sub` points at subband sample (granule start + 18 q) of the array subbands() filled, and 36
    // subband samples are read from there. `lines` takes 576 floats, subband-major: lines[sb * 18 + k].
    void granule (const float* sub, float* lines) const noexcept
    {
        for (int sb = 0; sb < kSubbands; ++sb)
        {
            // the 36 inputs of this subband, windowed, with the frequency inversion: the odd sample of an odd
            // subband changes sign, counted inside each block of 18
            float w[2 * kBlock];
            for (int n = 0; n < 2 * kBlock; ++n)
            {
                const float v = sub[(std::size_t) n * kSubbands + (std::size_t) sb] * win_[(std::size_t) n];
                w[n] = ((sb & 1) != 0 && ((n % kBlock) & 1) != 0) ? -v : v;
            }
            // the fold: quarters a b c d of nine samples each -> u = (-c reversed - d, a - b reversed)
            float u[kBlock];
            for (int n = 0; n < kBlock / 2; ++n)
            {
                u[n] = -w[kBlock + kBlock / 2 - 1 - n] - w[kBlock + kBlock / 2 + n];
                u[kBlock / 2 + n] = w[n] - w[kBlock - 1 - n];
            }
            float* o = lines + (std::size_t) sb * kBlock;
            for (int k = 0; k < kBlock; ++k) o[k] = 0.0f;
            for (int n = 0; n < kBlock; ++n)
            {
                const float un = u[n];
                const float* row = dct_[(std::size_t) n].data();
                for (int k = 0; k < kBlock; ++k) o[k] += un * row[k];
            }
        }
        // the alias reduction between neighbouring subbands, in the encoder's direction
        for (int sb = 1; sb < kSubbands; ++sb)
            for (int i = 0; i < kAliasPairs; ++i)
            {
                float& lo = lines[(std::size_t) (kBlock * sb - 1 - i)];
                float& hi = lines[(std::size_t) (kBlock * sb + i)];
                const float l = lo, h = hi;
                lo = l * cs_[(std::size_t) i] + h * ca_[(std::size_t) i];
                hi = h * cs_[(std::size_t) i] - l * ca_[(std::size_t) i];
            }
    }

private:
    std::array<float, kWindow> c_ {};
    std::array<std::array<float, 64>, kSubbands> poly_ {};
    std::array<float, 2 * kBlock> win_ {};
    std::array<std::array<float, kBlock>, kBlock> dct_ {};
    std::array<float, kAliasPairs> cs_ {}, ca_ {};
};

} // namespace felitronics::codecgrid
