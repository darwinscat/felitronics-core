// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <new>
#include <type_traits>
#include <vector>

//==============================================================================
// felitronics::core::fft — the FFT SEAM (the architecture keystone). A real-input FFT reached through a
// COMPILE-TIME backend (template/concept, never a vtable in the convolver hot path). Each backend owns
// its OWN spectrum layout (so a SIMD backend — pffft/vDSP — never pays an O(N) repack to a forced
// lowest-common-denominator); a consumer (the partitioned convolver) only ever calls forward / inverse /
// spectralMultiplyAdd and never indexes spectrum bins. A plan is built in prepare() (all allocation
// there); forward()/inverse()/spectralMultiplyAdd() are RT-safe (no alloc/lock/throw).
//
// Ships a scalar radix-2 reference backend so the core self-tests JUCE-free with no heavy dep. Real
// backends (pffft/kissfft = BSD on desktop+wasm; juce::dsp::FFT in the JUCE adapter; CMSIS embedded)
// plug in later as compiled targets satisfying RealFftBackend — pffft's native zconvolve_accumulate maps
// straight onto spectralMultiplyAdd.
//==============================================================================
namespace felitronics::core::fft
{

inline constexpr bool isPow2 (int n) noexcept { return n > 0 && (n & (n - 1)) == 0; }

// A real-FFT backend for size N (pow2). `spectrumFloats(N)` is the backend's spectrum buffer length in
// floats (layout is the backend's business). `inverse` is 1/N-normalized so a round-trip is identity.
// `spectralMultiplyAdd(a,b,acc)` does acc += a (.*) b in the backend's own layout (the convolver's MAC).
// LAYOUT LINEARITY — the layout must be a fixed LINEAR map of the transform: an elementwise float add/scale of
// two spectra must equal the spectrum of the same add/scale in the time domain. MatrixConvolver relies on this
// to build its on-the-fly M/S views (½(X_L ± X_R), per float). Holds for the scalar packed-Hermitian layout and
// pffft's z-order (both are permutations of the linear DFT); a polar magnitude/phase layout would NOT qualify.
template <class B>
concept RealFftBackend = requires (B b, const float* r, float* w, const float* s, float* acc, int n) {
    { B::spectrumFloats (n) } noexcept -> std::convertible_to<int>;
    { b.prepare (n) }         noexcept -> std::same_as<bool>;
    { b.forward (r, w) }      noexcept -> std::same_as<void>;   // real[N]      -> spectrum[spectrumFloats(N)]
    { b.inverse (s, w) }      noexcept -> std::same_as<void>;   // spectrum     -> real[N]  (1/N normalized)
    { b.spectralMultiplyAdd (s, s, acc) } noexcept -> std::same_as<void>;
};

// A backend whose spectrum is the packed-Hermitian layout the lineareq FIR designers hand-write directly
// (s[0]=DC, s[1]=Nyquist, s[2k]/s[2k+1]=Re[k]/Im[k]). ONLY such a backend may drive the DESIGN-time FFTs; a
// SIMD backend with its own opaque layout (pffft's z-order) must NOT (it would silently corrupt the FIR).
// The design FFTs are pinned to ScalarRadix2Real and static_assert this, so a wrong instantiation fails to
// COMPILE rather than mis-designing at runtime.
template <class F>
concept PackedHermitianSpectrum = RealFftBackend<F>
    && requires { requires std::same_as<std::remove_cvref_t<decltype (F::kPackedHermitianSpectrum)>, bool>;   // a bool member (not int/enum/fn that is merely truthy)
                  requires F::kPackedHermitianSpectrum; };                                                    // ...and it is true

//==============================================================================
// SEAM STORAGE — the frequency-domain buffers a SIMD backend convolves IN PLACE (pffft's
// zconvolve_accumulate, vDSP) must be SIMD-aligned; std::vector's default allocator does not guarantee it,
// and a backend-side padding trick is impossible because the convolver memcpys whole spectra between FDL
// rows (a byte copy preserves offset-in-buffer, not absolute alignment). So the alignment is engine-side:
// the convolvers hold their seam-crossing buffers in AlignedVector. kSeamAlignment=64 is a cacheline —
// a safe superset of every real-FFT SIMD width (SSE/NEON 16, AVX 32). ROW-STRIDE PRECONDITION: FDL / IR
// partition rows are packed at j*spectrumFloats(N), so a row is aligned only if spectrumFloats(N)*sizeof(float)
// is a multiple of kSeamAlignment. Holds for ScalarRadix2Real and pffft (spectrumFloats==N; pow2 N>=32 → a
// 64-multiple). A backend whose natural spectrum length is NOT a multiple (an N+2-float unpacked-Hermitian
// layout — juce::dsp::FFT, kissfft-real) MUST round its reported spectrumFloats up (harmless: rows are
// zero-init, forward writes only the natural length, the backend never reads the slack).
// Backend-agnostic: the scalar reference ignores the stronger alignment (unobservable in its output). All
// allocation is in prepare() (message thread), so the throwing aligned operator new is RT-fine.
//
// SEAM CONTRACT (relied on across the two backend instances a convolver holds — the message-thread build
// FFT and the audio FFT): for a given backend TYPE and size N, spectra are interchangeable — the layout is
// a pure function of (type, N), independent of the instance. True for ScalarRadix2Real and pffft; a plan/
// wisdom-dependent backend (e.g. FFTW with FFTW_MEASURE) would violate it and must not be used here.
inline constexpr std::size_t kSeamAlignment = 64;

template <class T, std::size_t Alignment = kSeamAlignment>
struct SeamAllocator
{
    static_assert ((Alignment & (Alignment - 1)) == 0, "kSeamAlignment must be a power of two");
    static_assert (Alignment >= alignof (T),            "kSeamAlignment must be >= alignof(T)");

    using value_type = T;
    template <class U> struct rebind { using other = SeamAllocator<U, Alignment>; };

    SeamAllocator() noexcept = default;
    template <class U> SeamAllocator (const SeamAllocator<U, Alignment>&) noexcept {}

    T* allocate (std::size_t n)
    {
        if (n > static_cast<std::size_t> (-1) / sizeof (T)) throw std::bad_array_new_length();
        return static_cast<T*> (::operator new (n * sizeof (T), std::align_val_t (Alignment)));
    }
    void deallocate (T* p, std::size_t) noexcept { ::operator delete (p, std::align_val_t (Alignment)); }

    template <class U> bool operator== (const SeamAllocator<U, Alignment>&) const noexcept { return true; }
    template <class U> bool operator!= (const SeamAllocator<U, Alignment>&) const noexcept { return false; }
};

// A std::vector whose storage is kSeamAlignment-aligned — for every convolver buffer that crosses the FFT
// seam (FDL rows, IR partition spectra, the input/view/accumulator spectra, and the forward/inverse real
// scratch). Time-domain-only buffers (head taps, cached tails) stay plain std::vector.
template <class T> using AlignedVector = std::vector<T, SeamAllocator<T>>;

//==============================================================================
// ScalarRadix2Real — the JUCE-free reference + spike default. Correctness-first (a full complex radix-2
// transform of the real input; a real backend is ~2× faster but this is the analytic reference the tests
// trust). Packed spectrum (N floats): s[0]=Re[0] (DC), s[1]=Re[N/2] (Nyquist), then s[2k]=Re[k],
// s[2k+1]=Im[k] for 1<=k<N/2.
class ScalarRadix2Real
{
public:
    static constexpr int  spectrumFloats (int n) noexcept { return n; }
    static constexpr bool kPackedHermitianSpectrum = true;   // s[0]=DC, s[1]=Nyq, s[2k]/s[2k+1]=Re/Im — the layout lineareq hand-writes

    bool prepare (int n) noexcept
    {
        if (! isPow2 (n) || n < 4) { n_ = 0; return false; }   // reject → UNPREPARED (house rule): a failed
        n_ = n;                                                //   RE-prepare must not keep the stale old plan
        scratch_.assign ((std::size_t) n, std::complex<float> {});   // ALLOC here (prepare) only
        return true;
    }

    int size() const noexcept { return n_; }

    void forward (const float* real, float* spec) noexcept
    {
        if (n_ <= 0) return;                                // unprepared — scratch_ empty (per-transform guard, not per-sample)
        for (int i = 0; i < n_; ++i) scratch_[(std::size_t) i] = std::complex<float> (real[i], 0.0f);
        transform (scratch_.data(), n_, false);
        spec[0] = scratch_[0].real();                       // DC (real)
        spec[1] = scratch_[(std::size_t) (n_ / 2)].real();  // Nyquist (real)
        for (int k = 1; k < n_ / 2; ++k)
        {
            spec[2 * k]     = scratch_[(std::size_t) k].real();
            spec[2 * k + 1] = scratch_[(std::size_t) k].imag();
        }
    }

    void inverse (const float* spec, float* real) noexcept
    {
        if (n_ <= 0) return;                                // unprepared — scratch_ empty
        scratch_[0]                            = std::complex<float> (spec[0], 0.0f);
        scratch_[(std::size_t) (n_ / 2)]       = std::complex<float> (spec[1], 0.0f);
        for (int k = 1; k < n_ / 2; ++k)       // Hermitian symmetry rebuilds the upper half
        {
            const std::complex<float> c (spec[2 * k], spec[2 * k + 1]);
            scratch_[(std::size_t) k]          = c;
            scratch_[(std::size_t) (n_ - k)]   = std::conj (c);
        }
        transform (scratch_.data(), n_, true);
        const float inv = 1.0f / (float) n_;
        for (int i = 0; i < n_; ++i) real[i] = scratch_[(std::size_t) i].real() * inv;
    }

    // acc[] += a[] (complex .*) b[] in the packed layout. DC + Nyquist are real-only.
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) const noexcept
    {
        acc[0] += a[0] * b[0];
        acc[1] += a[1] * b[1];
        for (int k = 1; k < n_ / 2; ++k)
        {
            const float ar = a[2 * k], ai = a[2 * k + 1], br = b[2 * k], bi = b[2 * k + 1];
            acc[2 * k]     += ar * br - ai * bi;
            acc[2 * k + 1] += ar * bi + ai * br;
        }
    }

private:
    // Iterative Cooley-Tukey radix-2, in place. Twiddles accumulated in double for reference accuracy;
    // butterflies in float (the hot-path type). Caller normalizes the inverse.
    static void transform (std::complex<float>* a, int n, bool inverse) noexcept
    {
        for (int i = 1, j = 0; i < n; ++i)              // bit-reversal permutation
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap (a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1)
        {
            const double ang = 2.0 * kPi / len * (inverse ? 1.0 : -1.0);
            const std::complex<double> wlen (std::cos (ang), std::sin (ang));
            for (int i = 0; i < n; i += len)
            {
                std::complex<double> w (1.0, 0.0);
                for (int k = 0; k < len / 2; ++k)
                {
                    const std::complex<float> u = a[i + k];
                    const std::complex<float> v = a[i + k + len / 2]
                                                * std::complex<float> ((float) w.real(), (float) w.imag());
                    a[i + k]            = u + v;
                    a[i + k + len / 2]  = u - v;
                    w *= wlen;
                }
            }
        }
    }

    int n_ = 0;
    std::vector<std::complex<float>> scratch_;
};

static_assert (RealFftBackend<ScalarRadix2Real>, "ScalarRadix2Real must satisfy the seam");
static_assert (PackedHermitianSpectrum<ScalarRadix2Real>,
               "ScalarRadix2Real is the packed-Hermitian design reference the lineareq FIR designers hand-pack");

// The default real-FFT backend (swap per tier later via the template seam).
using DefaultRealFft = ScalarRadix2Real;

} // namespace felitronics::core::fft
