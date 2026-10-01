// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid's complex FFT, windows and MDCT against their definitions written out directly in double: the DFT as
// an N x N sum, the MDCT as a 2M x M sum, the windows from their formulas and from the property each was built
// for (Princen-Bradley: w[n]^2 + w[n + M]^2 = 1).

#include <felitronics/codecgrid/Mdct.h>

#include "felitronics_test.h"

#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

using namespace felitronics::codecgrid;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
constexpr double kPi = 3.14159265358979323846;

struct Lcg
{
    std::uint32_t s;
    double next() noexcept { s = s * 1664525u + 1013904223u; return (double) (s >> 8) / 8388608.0 - 1.0; }   // [-1, 1)
};
}

int main()
{
    std::printf ("felitronics::codecgrid MDCT tests\n");

    group ("which sizes the scalar FFT takes");
    {
        ok (MixedRadixFft::supported (480) && MixedRadixFft::supported (512) && MixedRadixFft::supported (1) && MixedRadixFft::supported (15),
            "2^a 3^b 5^c: 480, 512, 15 and 1");
        ok (! MixedRadixFft::supported (7) && ! MixedRadixFft::supported (14) && ! MixedRadixFft::supported (0) && ! MixedRadixFft::supported (-4),
            "a factor of 7, zero and a negative size are refused");
        ok (! MixedRadixFft::supported (MixedRadixFft::kMaxSize * 2), "a size past the ceiling is refused");
        MixedRadixFft f;
        ok (! f.prepare (7) && f.size() == 0, "prepare() refuses what supported() refuses, and stays unprepared");
        ok (f.prepare (480) && f.size() == 480, "prepare (480) is accepted");
        ok (! f.prepare (0) && f.size() == 0, "a refused prepare() after a good one leaves nothing prepared");
    }

    group ("the FFT nulls against the DFT sum");
    {
        for (int n : { 1, 2, 3, 4, 5, 6, 8, 15, 16, 30, 60, 64, 120, 125, 243, 480, 512, 960, 1024 })
        {
            MixedRadixFft fft;
            if (! felitronics::test::run (fft.prepare (n))) continue;
            Lcg rng { 31u + (std::uint32_t) n };
            std::vector<std::complex<float>> in ((std::size_t) n), out ((std::size_t) n);
            for (auto& v : in) v = std::complex<float> ((float) rng.next(), (float) rng.next());
            fft.forward (in.data(), out.data());
            double worst = 0.0, peak = 0.0;
            for (int k = 0; k < n; ++k)
            {
                std::complex<double> acc (0.0, 0.0);
                for (int j = 0; j < n; ++j)
                {
                    const double a = -2.0 * kPi * (double) ((long long) k * j % n) / (double) n;
                    acc += std::complex<double> (in[(std::size_t) j].real(), in[(std::size_t) j].imag()) * std::complex<double> (std::cos (a), std::sin (a));
                }
                worst = std::fmax (worst, std::abs (acc - std::complex<double> (out[(std::size_t) k].real(), out[(std::size_t) k].imag())));
                peak = std::fmax (peak, std::abs (acc));
            }
            ok (std::isfinite (worst) && worst <= 4.0e-6 * std::fmax (peak, 1.0),
                "N = " + std::to_string (n) + ": worst " + std::to_string (worst) + " against a peak of " + std::to_string (peak));
        }
    }

    group ("the windows are the formulas, and each satisfies Princen-Bradley");
    {
        for (int frame : { 36, 1920, 2048 })
        {
            std::vector<double> w ((std::size_t) frame), k ((std::size_t) frame), scratch ((std::size_t) frame / 2 + 1);
            window::sine (w.data(), frame);
            window::kaiserBesselDerived (k.data(), frame, 4.0, scratch.data());
            double worstSine = 0.0, worstKbd = 0.0, formula = 0.0;
            for (int n = 0; n < frame / 2; ++n)
            {
                worstSine = std::fmax (worstSine, std::fabs (w[(std::size_t) n] * w[(std::size_t) n] + w[(std::size_t) (n + frame / 2)] * w[(std::size_t) (n + frame / 2)] - 1.0));
                worstKbd  = std::fmax (worstKbd,  std::fabs (k[(std::size_t) n] * k[(std::size_t) n] + k[(std::size_t) (n + frame / 2)] * k[(std::size_t) (n + frame / 2)] - 1.0));
            }
            for (int n = 0; n < frame; ++n) formula = std::fmax (formula, std::fabs (w[(std::size_t) n] - std::sin (kPi / frame * (n + 0.5))));
            const std::string tag = "frame " + std::to_string (frame) + ": ";
            ok (formula < 1.0e-15, tag + "the sine window is sin (pi / N (n + 1/2)) to the last bit or two");
            ok (worstSine < 1.0e-14 && worstKbd < 1.0e-13, tag + "sine and Kaiser-Bessel-derived windows are power-complementary");
            bool symmetric = true, rising = true;
            for (int n = 0; n < frame / 2; ++n)
            {
                symmetric = symmetric && std::fabs (k[(std::size_t) n] - k[(std::size_t) (frame - 1 - n)]) < 1.0e-16;
                if (n > 0) rising = rising && k[(std::size_t) n] > k[(std::size_t) (n - 1)];
            }
            // how near zero it starts depends on the length: a 2048-point window starts under 1e-3, a 36-point one cannot
            ok (symmetric && rising && k[0] > 0.0 && k[0] < (frame >= 1920 ? 1.0e-3 : 0.1), tag + "the KBD window is symmetric, rises monotonically, and starts near zero");
        }
        // the KBD window differs from the sine window where it should: narrower skirts, so smaller at the ends
        std::vector<double> s (2048), k (2048), scratch (1025);
        window::sine (s.data(), 2048);
        window::kaiserBesselDerived (k.data(), 2048, 4.0, scratch.data());
        ok (k[100] < 0.5 * s[100] && std::fabs (k[1023] - 1.0) < std::fabs (s[1023] - 1.0), "KBD (alpha 4) is below the sine window at the edge and closer to one at the centre");

        std::vector<double> c (1920);
        window::celt (c.data(), 1920, 120);
        bool zeros = true, ones = true;
        for (int n = 0; n < 420; ++n) zeros = zeros && std::fpclassify (c[(std::size_t) n]) == FP_ZERO && std::fpclassify (c[(std::size_t) (1919 - n)]) == FP_ZERO;
        for (int n = 540; n < 1380; ++n) ones = ones && std::fabs (c[(std::size_t) n] - 1.0) < 1.0e-300;
        ok (zeros && ones, "CELT's long window: 420 zeros at each end, ones through the middle");
        double worst = 0.0;
        for (int n = 0; n < 960; ++n) worst = std::fmax (worst, std::fabs (c[(std::size_t) n] * c[(std::size_t) n] + c[(std::size_t) (n + 960)] * c[(std::size_t) (n + 960)] - 1.0));
        ok (worst < 1.0e-14, "...and it is power-complementary over its 120-sample overlap");
        ok (std::fabs (c[420] - std::sin (0.5 * kPi * std::pow (std::sin (0.5 * kPi * 0.5 / 120.0), 2.0))) < 1.0e-15, "its first rising sample is the RFC's formula");
    }

    group ("the MDCT nulls against its definition, for the three frames the module uses");
    {
        struct Case { int frame; const char* name; };
        for (const Case cs : { Case { 2048, "AAC sine" }, Case { 2048, "AAC KBD" }, Case { 1920, "CELT" }, Case { 16, "tiny" } })
        {
            const int frame = cs.frame, m = frame / 2;
            std::vector<double> w ((std::size_t) frame), scratch ((std::size_t) m + 1);
            if (std::string (cs.name) == "AAC KBD") window::kaiserBesselDerived (w.data(), frame, 4.0, scratch.data());
            else if (std::string (cs.name) == "CELT") window::celt (w.data(), frame, 120);
            else window::sine (w.data(), frame);

            Mdct<> mdct;
            if (! felitronics::test::run (mdct.prepare (frame, w.data()))) continue;
            ok (mdct.hop() == m && mdct.frame() == frame, std::string (cs.name) + ": hop and frame are published");
            Lcg rng { 5u + (std::uint32_t) frame };
            std::vector<float> x ((std::size_t) frame), got ((std::size_t) m);
            for (auto& v : x) v = (float) (0.5 * rng.next());
            mdct.transform (x.data(), got.data());
            double worst = 0.0, peak = 0.0;
            for (int k = 0; k < m; ++k)
            {
                double acc = 0.0;
                for (int n = 0; n < frame; ++n)
                    acc += (double) (float) w[(std::size_t) n] * (double) x[(std::size_t) n] * std::cos (kPi / m * (n + 0.5 + 0.5 * m) * (k + 0.5));
                worst = std::fmax (worst, std::fabs (acc - (double) got[(std::size_t) k]));
                peak = std::fmax (peak, std::fabs (acc));
            }
            std::printf ("    %-9s frame %4d: peak %.3f, worst difference %.3e (%.1f dB under the peak)\n", cs.name, frame, peak, worst,
                         20.0 * std::log10 (peak / std::fmax (worst, 1e-300)));
            ok (std::isfinite (worst) && peak > 0.1 && worst < 1.0e-5 * peak, std::string (cs.name) + ": single-precision agreement with the 2M x M sum");
        }
    }

    group ("what the MDCT refuses");
    {
        Mdct<> mdct;
        std::vector<double> w (2048, 1.0);
        ok (! Mdct<>::supported (2046) && ! Mdct<>::supported (4) && ! Mdct<>::supported (0) && ! Mdct<>::supported (28),
            "a frame that is not a multiple of four, too short, or whose quarter has a factor of 7");
        ok (! mdct.prepare (28, w.data()) && mdct.hop() == 0, "prepare() refuses the same, and stays unprepared");
        ok (! mdct.prepare (2048, nullptr) && mdct.hop() == 0, "a null window is refused");
        ok (mdct.prepare (2048, w.data()) && mdct.hop() == 1024, "a good prepare()");
        ok (! mdct.prepare (30, w.data()) && mdct.hop() == 0, "a refused prepare() after a good one leaves nothing prepared");
    }

    return felitronics::test::report();
}
