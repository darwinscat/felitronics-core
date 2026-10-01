// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/codecgrid/Windows.h>
#include <felitronics/core/DetMath.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <vector>

//==============================================================================
// felitronics::codecgrid — the MDCT of the transform codecs (AAC-LC, CELT), its windows, and the complex FFT it
// runs on.
//
// THE TRANSFORM. X[k] = sum_{n < 2M} w[n] x[n] cos (pi / M (n + 1/2 + M/2) (k + 1/2)), k < M — 2M samples in,
// M coefficients out, hop M. It is computed the usual way: the 2M windowed samples fold to M, the M fold to a
// complex sequence of M/2 with a pre-twiddle, ONE complex FFT of M/2 points, a post-twiddle. So an AAC frame
// (2048 samples) costs a 512-point complex FFT and a CELT frame (1920 samples) a 480-point one.
//
// THE SEAM, AND WHY IT IS NOT core::fft::RealFftBackend. That concept is the convolvers': a real-input transform
// of a power-of-two size whose spectrum layout is the backend's own business, with a spectral multiply-add. An
// MDCT wants none of that and two things it does not offer — a COMPLEX transform, and 480 points, which is
// 2^5 x 3 x 5: CELT's frame is 20 ms at 48 kHz, not a power of two. So this module states the little it needs
// (ComplexFftBackend below) and ships a scalar mixed-radix reference for it. A compiled SIMD backend plugs in
// through the same concept.
//
// ARITHMETIC. Single precision throughout the transform; twiddles and windows designed in double through
// core::det and rounded once, so the tables are the same on every platform. The transform's own rounding is
// NOT promised bit-identical across platforms or backends — contraction and SIMD differ — and nothing downstream
// needs it to be: what is read off these coefficients is floored 80 dB under their own level.
//==============================================================================
namespace felitronics::codecgrid
{

// An unnormalised forward complex FFT of a size fixed at prepare(): out[k] = sum_n in[n] exp (-2 pi i k n / N).
// `in` and `out` are distinct arrays of N. prepare() is the only allocation; forward() is no alloc/lock/throw.
template <class B>
concept ComplexFftBackend = requires (B b, const B cb, const std::complex<float>* in, std::complex<float>* out, int n) {
    { B::supported (n) } noexcept -> std::same_as<bool>;
    { b.prepare (n) } -> std::same_as<bool>;
    { cb.size() } noexcept -> std::same_as<int>;
    { cb.forward (in, out) } noexcept -> std::same_as<void>;
};

//==============================================================================
// MixedRadixFft — the scalar reference: any N = 2^a 3^b 5^c, decimation in time. The factors are taken 5, 3,
// then one 2 when the power of two is odd, then 4s — so the leaves of the recursion are 4-point transforms
// wherever N has a factor of four, and the butterflies with the most arithmetic run the fewest times. Each radix
// has its own butterfly (the ones kissfft made familiar).
//
// The arithmetic is spelled out on real and imaginary parts. std::complex's operator* must honour infinities
// and NaNs the way Annex G of C says, which costs a branch per product and keeps the loops from vectorising;
// nothing here needs it — a non-finite input gives a non-finite output either way.
//==============================================================================
class MixedRadixFft
{
public:
    using Complex = std::complex<float>;
    static constexpr int kMaxSize = 1 << 20;

    static constexpr bool supported (int n) noexcept
    {
        if (n < 1 || n > kMaxSize) return false;
        for (int p : { 2, 3, 5 }) while (n % p == 0) n /= p;
        return n == 1;
    }

    [[nodiscard]] bool prepare (int n)
    {
        n_ = 0;
        stages_ = 0;
        if (! supported (n)) return false;
        int rest = n, twos = 0;
        for (int t = n; t % 2 == 0; t /= 2) ++twos;
        auto take = [this, &rest] (int p) noexcept { rest /= p; radix_[stages_] = p; remain_[stages_] = rest; ++stages_; };
        while (rest % 5 == 0) take (5);
        while (rest % 3 == 0) take (3);
        if ((twos & 1) != 0) take (2);
        while (rest % 4 == 0) take (4);
        twiddle_.assign (2u * (std::size_t) n, 0.0f);
        for (int k = 0; k < n; ++k)
        {
            const double a = 2.0 * core::kPi * (double) k / (double) n;
            twiddle_[2u * (std::size_t) k] = (float) core::det::cos (a);
            twiddle_[2u * (std::size_t) k + 1u] = (float) -core::det::sin (a);
        }
        n_ = n;
        return true;
    }

    int size() const noexcept { return n_; }

    // std::complex<float> is, by the standard, two floats: real then imaginary.
    void forward (const Complex* in, Complex* out) const noexcept
    {
        if (n_ <= 0) return;
        if (n_ == 1) { out[0] = in[0]; return; }
        work (reinterpret_cast<float*> (out), reinterpret_cast<const float*> (in), 1, 0);
    }

private:
    static constexpr int kMaxStages = 24;

    void work (float* out, const float* in, int stride, int stage) const noexcept
    {
        const int p = radix_[stage], m = remain_[stage];
        if (m == 1) { leaf (out, in, stride, p); return; }
        for (int i = 0; i < p; ++i)
            work (out + 2 * i * m, in + 2 * i * stride, stride * p, stage + 1);
        const float* tw = twiddle_.data();
        const int s = stride;
        float* f0 = out;
        float* f1 = out + 2 * m;
        if (p == 2)
        {
            for (int k = 0; k < m; ++k)
            {
                const float wr = tw[2 * k * s], wi = tw[2 * k * s + 1];
                const float tr = f1[2 * k] * wr - f1[2 * k + 1] * wi, ti = f1[2 * k] * wi + f1[2 * k + 1] * wr;
                f1[2 * k] = f0[2 * k] - tr;      f1[2 * k + 1] = f0[2 * k + 1] - ti;
                f0[2 * k] += tr;                 f0[2 * k + 1] += ti;
            }
        }
        else if (p == 4)
        {
            float* f2 = out + 4 * m;
            float* f3 = out + 6 * m;
            for (int k = 0; k < m; ++k)
            {
                const float w1r = tw[2 * k * s], w1i = tw[2 * k * s + 1];
                const float w2r = tw[4 * k * s], w2i = tw[4 * k * s + 1];
                const float w3r = tw[6 * k * s], w3i = tw[6 * k * s + 1];
                const float ar = f1[2 * k] * w1r - f1[2 * k + 1] * w1i, ai = f1[2 * k] * w1i + f1[2 * k + 1] * w1r;
                const float br = f2[2 * k] * w2r - f2[2 * k + 1] * w2i, bi = f2[2 * k] * w2i + f2[2 * k + 1] * w2r;
                const float cr = f3[2 * k] * w3r - f3[2 * k + 1] * w3i, ci = f3[2 * k] * w3i + f3[2 * k + 1] * w3r;
                const float dr = f0[2 * k] - br, di = f0[2 * k + 1] - bi;
                const float sr = f0[2 * k] + br, si = f0[2 * k + 1] + bi;
                const float er = ar + cr, ei = ai + ci, gr = ar - cr, gi = ai - ci;
                f0[2 * k] = sr + er;             f0[2 * k + 1] = si + ei;
                f2[2 * k] = sr - er;             f2[2 * k + 1] = si - ei;
                f1[2 * k] = dr + gi;             f1[2 * k + 1] = di - gr;
                f3[2 * k] = dr - gi;             f3[2 * k + 1] = di + gr;
            }
        }
        else if (p == 3)
        {
            float* f2 = out + 4 * m;
            const float epi = tw[2 * m * s + 1];                 // the imaginary part of exp (-2 pi i / 3)
            for (int k = 0; k < m; ++k)
            {
                const float w1r = tw[2 * k * s], w1i = tw[2 * k * s + 1];
                const float w2r = tw[4 * k * s], w2i = tw[4 * k * s + 1];
                const float ar = f1[2 * k] * w1r - f1[2 * k + 1] * w1i, ai = f1[2 * k] * w1i + f1[2 * k + 1] * w1r;
                const float br = f2[2 * k] * w2r - f2[2 * k + 1] * w2i, bi = f2[2 * k] * w2i + f2[2 * k + 1] * w2r;
                const float sr = ar + br, si = ai + bi;
                const float dr = (ar - br) * epi, di = (ai - bi) * epi;
                const float hr = f0[2 * k] - 0.5f * sr, hi = f0[2 * k + 1] - 0.5f * si;
                f0[2 * k] += sr;                 f0[2 * k + 1] += si;
                f2[2 * k] = hr + di;             f2[2 * k + 1] = hi - dr;
                f1[2 * k] = hr - di;             f1[2 * k + 1] = hi + dr;
            }
        }
        else                                                                    // 5
        {
            float* f2 = out + 4 * m;
            float* f3 = out + 6 * m;
            float* f4 = out + 8 * m;
            const float yar = tw[2 * m * s], yai = tw[2 * m * s + 1];      // exp (-2 pi i / 5)
            const float ybr = tw[4 * m * s], ybi = tw[4 * m * s + 1];      // exp (-4 pi i / 5)
            for (int k = 0; k < m; ++k)
            {
                const float w1r = tw[2 * k * s], w1i = tw[2 * k * s + 1];
                const float w2r = tw[4 * k * s], w2i = tw[4 * k * s + 1];
                const float w3r = tw[6 * k * s], w3i = tw[6 * k * s + 1];
                const float w4r = tw[8 * k * s], w4i = tw[8 * k * s + 1];
                const float s0r = f0[2 * k], s0i = f0[2 * k + 1];
                const float s1r = f1[2 * k] * w1r - f1[2 * k + 1] * w1i, s1i = f1[2 * k] * w1i + f1[2 * k + 1] * w1r;
                const float s2r = f2[2 * k] * w2r - f2[2 * k + 1] * w2i, s2i = f2[2 * k] * w2i + f2[2 * k + 1] * w2r;
                const float s3r = f3[2 * k] * w3r - f3[2 * k + 1] * w3i, s3i = f3[2 * k] * w3i + f3[2 * k + 1] * w3r;
                const float s4r = f4[2 * k] * w4r - f4[2 * k + 1] * w4i, s4i = f4[2 * k] * w4i + f4[2 * k + 1] * w4r;
                const float s7r = s1r + s4r, s7i = s1i + s4i, s10r = s1r - s4r, s10i = s1i - s4i;
                const float s8r = s2r + s3r, s8i = s2i + s3i, s9r = s2r - s3r, s9i = s2i - s3i;
                f0[2 * k] = s0r + s7r + s8r;     f0[2 * k + 1] = s0i + s7i + s8i;
                const float s5r = s0r + s7r * yar + s8r * ybr, s5i = s0i + s7i * yar + s8i * ybr;
                const float s6r = s10i * yai + s9i * ybi, s6i = -s10r * yai - s9r * ybi;
                f1[2 * k] = s5r - s6r;           f1[2 * k + 1] = s5i - s6i;
                f4[2 * k] = s5r + s6r;           f4[2 * k + 1] = s5i + s6i;
                const float s11r = s0r + s7r * ybr + s8r * yar, s11i = s0i + s7i * ybr + s8i * yar;
                const float s12r = -s10i * ybi + s9i * yai, s12i = s10r * ybi - s9r * yai;
                f2[2 * k] = s11r + s12r;         f2[2 * k + 1] = s11i + s12i;
                f3[2 * k] = s11r - s12r;         f3[2 * k + 1] = s11i - s12i;
            }
        }
    }

    // A transform of p points taken `stride` apart: the leaf of the recursion.
    void leaf (float* out, const float* in, int stride, int p) const noexcept
    {
        const int s = 2 * stride;
        if (p == 4)
        {
            const float a0r = in[0], a0i = in[1], a1r = in[s], a1i = in[s + 1], a2r = in[2 * s], a2i = in[2 * s + 1], a3r = in[3 * s], a3i = in[3 * s + 1];
            const float sr = a0r + a2r, si = a0i + a2i, dr = a0r - a2r, di = a0i - a2i;
            const float er = a1r + a3r, ei = a1i + a3i, gr = a1r - a3r, gi = a1i - a3i;
            out[0] = sr + er;    out[1] = si + ei;
            out[4] = sr - er;    out[5] = si - ei;
            out[2] = dr + gi;    out[3] = di - gr;
            out[6] = dr - gi;    out[7] = di + gr;
            return;
        }
        if (p == 2)
        {
            const float a0r = in[0], a0i = in[1], a1r = in[s], a1i = in[s + 1];
            out[0] = a0r + a1r;  out[1] = a0i + a1i;
            out[2] = a0r - a1r;  out[3] = a0i - a1i;
            return;
        }
        const float* tw = twiddle_.data();
        const int step = n_ / p;
        for (int k = 0; k < p; ++k)
        {
            float accR = in[0], accI = in[1];
            for (int q = 1; q < p; ++q)
            {
                const int t = 2 * (((k * q) % p) * step);
                const float xr = in[q * s], xi = in[q * s + 1];
                accR += xr * tw[t] - xi * tw[t + 1];
                accI += xr * tw[t + 1] + xi * tw[t];
            }
            out[2 * k] = accR;
            out[2 * k + 1] = accI;
        }
    }

    int n_ = 0, stages_ = 0;
    int radix_[kMaxStages] {}, remain_[kMaxStages] {};
    std::vector<float> twiddle_;            // cos, -sin, interleaved
};

static_assert (ComplexFftBackend<MixedRadixFft>, "MixedRadixFft must satisfy the seam");

//==============================================================================
// Mdct — one frame of 2M samples to M coefficients, on a ComplexFftBackend of M/2 points.
//==============================================================================
template <ComplexFftBackend Fft = MixedRadixFft>
class Mdct
{
public:
    using Complex = std::complex<float>;

    // The frame lengths this class takes: a multiple of four whose quarter the backend transforms.
    static bool supported (int frame) noexcept { return frame >= 8 && frame % 4 == 0 && Fft::supported (frame / 4); }

    // window: `frame` doubles. Allocates; false on a frame the backend cannot transform.
    [[nodiscard]] bool prepare (int frame, const double* window)
    {
        hop_ = 0;
        if (! supported (frame) || window == nullptr) return false;
        const int hop = frame / 2, quarter = hop / 2;
        if (! fft_.prepare (quarter)) return false;
        window_.assign ((std::size_t) frame, 0.0f);
        for (int n = 0; n < frame; ++n) window_[(std::size_t) n] = (float) window[n];
        pre_.assign ((std::size_t) quarter, Complex {});
        post_.assign ((std::size_t) quarter, Complex {});
        for (int n = 0; n < quarter; ++n)
        {
            // exp (-i pi (4n + 1) / (4M)) and exp (-i pi n / M)
            const double a = core::kPi * (double) (4 * n + 1) / (double) (4 * hop);
            const double b = core::kPi * (double) n / (double) hop;
            pre_[(std::size_t) n]  = Complex ((float) core::det::cos (a), (float) -core::det::sin (a));
            post_[(std::size_t) n] = Complex ((float) core::det::cos (b), (float) -core::det::sin (b));
        }
        fold_.assign ((std::size_t) hop, 0.0f);
        in_.assign ((std::size_t) quarter, Complex {});
        out_.assign ((std::size_t) quarter, Complex {});
        hop_ = hop;
        return true;
    }

    int hop() const noexcept { return hop_; }           // M: the coefficients per frame, and the hop
    int frame() const noexcept { return 2 * hop_; }

    // x: 2M samples. coefficients: M floats. Not const: the transform's scratch is the object's.
    void transform (const float* x, float* coefficients) noexcept
    {
        const int m = hop_, q = m / 2;
        const float* w = window_.data();
        float* u = fold_.data();
        for (int i = 0; i < q; ++i)
        {
            u[i]     = -x[m + q - 1 - i] * w[m + q - 1 - i] - x[m + q + i] * w[m + q + i];
            u[q + i] =  x[i] * w[i] - x[m - 1 - i] * w[m - 1 - i];
        }
        float* in = reinterpret_cast<float*> (in_.data());
        const float* pre = reinterpret_cast<const float*> (pre_.data());
        for (int n = 0; n < q; ++n)
        {
            const float ar = u[2 * n], ai = u[m - 1 - 2 * n], wr = pre[2 * n], wi = pre[2 * n + 1];
            in[2 * n] = ar * wr - ai * wi;
            in[2 * n + 1] = ar * wi + ai * wr;
        }
        fft_.forward (in_.data(), out_.data());
        const float* out = reinterpret_cast<const float*> (out_.data());
        const float* post = reinterpret_cast<const float*> (post_.data());
        for (int k = 0; k < q; ++k)
        {
            const float yr = out[2 * k], yi = out[2 * k + 1], wr = post[2 * k], wi = post[2 * k + 1];
            coefficients[2 * k] = yr * wr - yi * wi;
            coefficients[m - 1 - 2 * k] = -(yr * wi + yi * wr);
        }
    }

    static std::size_t bytesFor (int frame) noexcept
    {
        const std::size_t hop = (std::size_t) frame / 2u, quarter = hop / 2u;
        return sizeof (float) * ((std::size_t) frame + hop) + sizeof (Complex) * 5u * quarter;   // window, fold; pre, post, in, out, the backend's twiddles
    }

private:
    Fft fft_;
    int hop_ = 0;
    std::vector<float> window_, fold_;
    std::vector<Complex> pre_, post_, in_, out_;
};

} // namespace felitronics::codecgrid
