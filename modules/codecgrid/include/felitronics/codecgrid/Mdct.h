// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

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
// MixedRadixFft — the scalar reference: any N = 2^a 3^b 5^c, decimation in time by the factors 4, 2, 3, 5 in
// that order (the decomposition kissfft uses). Radix 2 and 4 have their own butterflies; 3 and 5 share a
// generic one.
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
        int rest = n;
        for (int p : { 4, 2, 3, 5 })
            while (rest % p == 0)
            {
                rest /= p;
                radix_[stages_] = p;
                remain_[stages_] = rest;
                ++stages_;
            }
        twiddle_.assign ((std::size_t) n, Complex {});
        for (int k = 0; k < n; ++k)
        {
            const double a = 2.0 * core::kPi * (double) k / (double) n;
            twiddle_[(std::size_t) k] = Complex ((float) core::det::cos (a), (float) -core::det::sin (a));
        }
        n_ = n;
        return true;
    }

    int size() const noexcept { return n_; }

    void forward (const Complex* in, Complex* out) const noexcept
    {
        if (n_ <= 0) return;
        if (n_ == 1) { out[0] = in[0]; return; }
        work (out, in, 1, 0);
    }

private:
    static constexpr int kMaxStages = 24;

    void work (Complex* out, const Complex* in, int stride, int stage) const noexcept
    {
        const int p = radix_[stage], m = remain_[stage];
        if (m == 1)
            for (int i = 0; i < p; ++i) out[i] = in[(std::size_t) i * (std::size_t) stride];
        else
            for (int i = 0; i < p; ++i) work (out + (std::size_t) i * (std::size_t) m, in + (std::size_t) i * (std::size_t) stride, stride * p, stage + 1);

        const Complex* tw = twiddle_.data();
        if (p == 2)
        {
            for (int k = 0; k < m; ++k)
            {
                const Complex t = out[m + k] * tw[(std::size_t) k * (std::size_t) stride];
                out[m + k] = out[k] - t;
                out[k] += t;
            }
        }
        else if (p == 4)
        {
            for (int k = 0; k < m; ++k)
            {
                const std::size_t s = (std::size_t) k * (std::size_t) stride;
                const Complex a = out[m + k] * tw[s], b = out[2 * m + k] * tw[2 * s], c = out[3 * m + k] * tw[3 * s];
                const Complex d = out[k] - b;
                out[k] += b;
                const Complex e = a + c, f = a - c;
                out[2 * m + k] = out[k] - e;
                out[k] += e;
                out[m + k]     = Complex (d.real() + f.imag(), d.imag() - f.real());
                out[3 * m + k] = Complex (d.real() - f.imag(), d.imag() + f.real());
            }
        }
        else
        {
            Complex scratch[5];
            for (int u = 0; u < m; ++u)
            {
                for (int q = 0; q < p; ++q) scratch[q] = out[u + q * m];
                for (int q1 = 0; q1 < p; ++q1)
                {
                    const int k = u + q1 * m;
                    int idx = 0;
                    Complex acc = scratch[0];
                    for (int q = 1; q < p; ++q)
                    {
                        idx += stride * k;
                        if (idx >= n_) idx %= n_;
                        acc += scratch[q] * tw[(std::size_t) idx];
                    }
                    out[k] = acc;
                }
            }
        }
    }

    int n_ = 0, stages_ = 0;
    int radix_[kMaxStages] {}, remain_[kMaxStages] {};
    std::vector<Complex> twiddle_;
};

static_assert (ComplexFftBackend<MixedRadixFft>, "MixedRadixFft must satisfy the seam");

//==============================================================================
// The windows. Each fills w[0 .. 2M - 1] for a frame of 2M samples.
//==============================================================================
namespace window
{
    // sin (pi / (2M) (n + 1/2)) — AAC's sine window, and MP3's at 2M = 36.
    inline void sine (double* w, int frame) noexcept
    {
        for (int n = 0; n < frame; ++n) w[n] = core::det::sin (core::kPi / (double) frame * ((double) n + 0.5));
    }

    namespace detail
    {
        // The modified Bessel function I0 by its series; 60 terms settle it far past double for the arguments used.
        inline double bessel0 (double x) noexcept
        {
            double sum = 1.0, term = 1.0;
            for (int k = 1; k < 60; ++k)
            {
                const double h = x / (2.0 * (double) k);
                term *= h * h;
                sum += term;
            }
            return sum;
        }
    }

    // The Kaiser-Bessel-derived window: the square root of the running sum of a Kaiser window of M + 1 points
    // with beta = pi alpha, mirrored. AAC's long window uses alpha = 4. `scratch` takes M + 1 doubles.
    inline void kaiserBesselDerived (double* w, int frame, double alpha, double* scratch) noexcept
    {
        const int half = frame / 2;
        const double beta = core::kPi * alpha, denom = detail::bessel0 (beta);
        double acc = 0.0;
        for (int i = 0; i <= half; ++i)
        {
            const double r = 2.0 * (double) i / (double) half - 1.0;
            acc += detail::bessel0 (beta * std::sqrt (1.0 - r * r)) / denom;
            scratch[i] = acc;
        }
        for (int i = 0; i < half; ++i)
        {
            w[i] = std::sqrt (scratch[i] / acc);
            w[frame - 1 - i] = w[i];
        }
    }

    // CELT's long-block window (RFC 6716, 4.3.7): of the 2M samples only the middle carries anything — zeros,
    // a rise of `overlap` samples, ones, the fall, zeros. The rise is sin (pi/2 sin^2 (pi/2 (i + 1/2) / overlap)),
    // which is power-complementary with its own mirror.
    inline void celt (double* w, int frame, int overlap) noexcept
    {
        const int hop = frame / 2, edge = (hop - overlap) / 2;
        for (int n = 0; n < frame; ++n) w[n] = 0.0;
        for (int i = 0; i < overlap; ++i)
        {
            const double s = core::det::sin (0.5 * core::kPi * ((double) i + 0.5) / (double) overlap);
            const double v = core::det::sin (0.5 * core::kPi * s * s);
            w[edge + i] = v;
            w[frame - 1 - edge - i] = v;
        }
        for (int n = edge + overlap; n < frame - edge - overlap; ++n) w[n] = 1.0;
    }
}

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
        for (int n = 0; n < q; ++n) in_[(std::size_t) n] = Complex (u[2 * n], u[m - 1 - 2 * n]) * pre_[(std::size_t) n];
        fft_.forward (in_.data(), out_.data());
        for (int k = 0; k < q; ++k)
        {
            const Complex y = out_[(std::size_t) k] * post_[(std::size_t) k];
            coefficients[2 * k] = y.real();
            coefficients[m - 1 - 2 * k] = -y.imag();
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
