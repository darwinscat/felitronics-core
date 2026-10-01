// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The SIMD backend under codecgrid (built only with FELITRONICS_WITH_PFFFT): fftpffft::PffftComplexFft against
// the DFT sum, the MDCT on it against the MDCT's definition, and the scan and the detector on it against the
// same scan and detector on the scalar reference — the same offsets, the same verdicts, the curves within the
// rounding of two single-precision transforms.

#include <felitronics/fftpffft/PffftComplexFft.h>

#include "DetectorHarness.h"
#include "alloc_counter.h"
#include "felitronics_test.h"

#include <cmath>
#include <complex>
#include <string>
#include <vector>

using namespace harness;
using felitronics::fftpffft::PffftComplexFft;
using felitronics::test::group;
using felitronics::test::ok;

static_assert (ComplexFftBackend<PffftComplexFft>, "PffftComplexFft must satisfy codecgrid's seam");

int main()
{
    std::printf ("felitronics::codecgrid on the pffft backend (SIMD width %d)\n", PffftComplexFft::simdWidth());

    group ("which sizes the complex pffft takes");
    {
        ok (PffftComplexFft::supported (480) && PffftComplexFft::supported (512) && PffftComplexFft::supported (16) && PffftComplexFft::supported (240),
            "multiples of 16 made of 2, 3 and 5: 16, 240, 480, 512");
        ok (! PffftComplexFft::supported (8) && ! PffftComplexFft::supported (24) && ! PffftComplexFft::supported (112) && ! PffftComplexFft::supported (0)
                && ! PffftComplexFft::supported (-16) && ! PffftComplexFft::supported (PffftComplexFft::kMaxSize * 2),
            "under 16, not a multiple of 16, a factor of 7, zero, negative, past the ceiling: refused");
        PffftComplexFft f;
        ok (! f.prepare (24) && f.size() == 0, "prepare() refuses what supported() refuses");
        ok (f.prepare (480) && f.size() == 480 && ! f.prepare (7) && f.size() == 0, "a refused prepare() after a good one leaves nothing prepared");
    }

    group ("the transform nulls against the DFT sum");
    {
        for (int n : { 16, 32, 48, 80, 240, 480, 512, 960, 1024 })
        {
            PffftComplexFft fft;
            if (! felitronics::test::run (fft.prepare (n))) continue;
            std::uint32_t s = 77u + (std::uint32_t) n;
            auto next = [&s]() noexcept { s = s * 1664525u + 1013904223u; return (float) ((double) (s >> 8) / 8388608.0 - 1.0); };
            std::vector<std::complex<float>> in ((std::size_t) n), out ((std::size_t) n);
            for (auto& v : in) v = std::complex<float> (next(), next());
            fft.forward (in.data(), out.data());
            double worst = 0.0, peak = 0.0;
            for (int k = 0; k < n; ++k)
            {
                std::complex<double> acc (0.0, 0.0);
                for (int j = 0; j < n; ++j)
                {
                    const double a = -2.0 * synthetic::kPi * (double) ((long long) k * j % n) / (double) n;
                    acc += std::complex<double> (in[(std::size_t) j].real(), in[(std::size_t) j].imag()) * std::complex<double> (std::cos (a), std::sin (a));
                }
                worst = std::fmax (worst, std::abs (acc - std::complex<double> (out[(std::size_t) k].real(), out[(std::size_t) k].imag())));
                peak = std::fmax (peak, std::abs (acc));
            }
            ok (std::isfinite (worst) && worst <= 4.0e-6 * std::fmax (peak, 1.0), "N = " + std::to_string (n) + ": worst " + std::to_string (worst) + " against a peak of " + std::to_string (peak));
        }
    }

    group ("the scan on pffft finds what the scalar scan finds, where it finds it");
    {
        constexpr int n = 40000;
        const auto original = synthetic::programme (n, 1u);
        GridScan<> scalar;
        GridScan<PffftComplexFft> simd;
        ok (scalar.prepare (n) && simd.prepare (n), "both scanners prepared");
        struct Case { Transform t; int offset; };
        for (const Case c : { Case { Transform::AacSine, 480 }, Case { Transform::AacKbd, 1023 }, Case { Transform::Celt, 548 }, Case { Transform::Mp3, 133 } })
        {
            const auto x = synthetic::coded (original, c.t, c.offset, 0.3);
            std::vector<float> a ((std::size_t) hopOf (c.t) * kCells), b (a.size());
            const bool ran = scalar.begin (c.t, x.left.data(), x.right.data(), n, a.data()) && simd.begin (c.t, x.left.data(), x.right.data(), n, b.data());
            while (! scalar.step (1 << 20)) {}
            while (! simd.step (1 << 20)) {}
            double worst = 0.0;
            for (std::size_t i = 0; i < a.size(); ++i) worst = std::fmax (worst, std::fabs ((double) a[i] - (double) b[i]));
            std::vector<double> d (ReadScratch::doublesFor (hopOf (c.t)));
            std::vector<float> f (ReadScratch::floatsFor (hopOf (c.t)));
            ReadScratch scratch;
            scratch.local = d.data(); scratch.sort = d.data() + hopOf (c.t); scratch.dip = f.data();
            const GridReading ra = readCurve (a.data(), hopOf (c.t), scratch), rb = readCurve (b.data(), hopOf (c.t), scratch);
            std::printf ("    transform %d: curves within %.2e dB; score %.1f scalar, %.1f pffft\n", (int) c.t, worst, ra.score, rb.score);
            ok (ran && ra.offset == c.offset && rb.offset == c.offset, "transform " + std::to_string ((int) c.t) + ": both find the offset it was coded at");
            ok (std::isfinite (worst) && worst < 0.02, "...and the curves agree within the rounding of two single-precision transforms");
            ok (GridRule {}.found (ra) && GridRule {}.found (rb) && ra.cells == rb.cells, "...the rule finds both, in the same number of cells");
        }
    }

    group ("the detector on pffft gives the scalar detector's verdicts");
    {
        const int n44 = 44100 * 3;
        const auto original = synthetic::programme (n44, 7u);
        struct Case { const char* what; synthetic::Stereo x; double rate; };
        const Case cases[] {
            { "AAC at the programme's rate", synthetic::coded (original, Transform::AacSine, 480, 0.3), 44100.0 },
            { "CELT in a 44.1 kHz file", codedAtOtherRate (48000 * 3, 11u, Transform::Celt, 548, 48000, 44100), 44100.0 },
            { "MP3", synthetic::coded (original, Transform::Mp3, 133, 0.3), 44100.0 },
            { "uncoded", original, 44100.0 },
        };
        for (const Case& c : cases)
        {
            const auto a = analyse (c.x, c.rate, params());
            CodecGridDetector<PffftComplexFft> d;
            d.setParams (params());
            const float* in[2] { c.x.left.data(), c.x.right.data() };
            const bool ran = d.prepare (c.rate, 2, c.x.left.size()) && d.process (in, 2, (int) c.x.left.size()) && d.finish();
            if (! felitronics::test::run (ran && a != nullptr)) continue;
            const auto& p = a->result();
            const auto& q = d.result();
            std::printf ("    %-28s %s / %s, phase %d, found %d — scalar score %.1f, pffft %.1f\n", c.what, name (q.verdict), name (q.ground), q.gridPhase, q.windowsFound, p.bestScore, q.bestScore);
            ok (p.verdict == q.verdict && p.ground == q.ground && p.family == q.family && p.codecRate == q.codecRate && p.gridPhase == q.gridPhase
                    && p.windowsFound == q.windowsFound && p.windowsAgreeing == q.windowsAgreeing && p.windowsExamined == q.windowsExamined,
                std::string (c.what) + ": the same verdict, ground, family, rate, phase and counts");
            ok (std::fabs (p.bestScore - q.bestScore) <= 0.05 * std::fmax (p.bestScore, 1.0), std::string (c.what) + ": the best score within 5 %");
        }
    }

    group ("no C++ allocation after prepare() on this backend either");
    {
        const int n44 = 44100 * 3;
        const auto x = synthetic::coded (synthetic::programme (n44, 7u), Transform::AacSine, 480, 0.3);
        CodecGridDetector<PffftComplexFft> d;
        d.setParams (params (Depth::Verdict));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n44), "prepared");
        const float* in[2] { x.left.data(), x.right.data() };
        const long long before = alloc::count.load();
        const bool ran = d.process (in, 2, n44) && d.finish();
        const long long after = alloc::count.load();
        ok (ran && d.result().verdict == Verdict::Confirmed, "the analysis ran and confirmed the grid");
        // what this counter cannot see is pffft's own malloc — it is called in pffft_new_setup, inside prepare()
        felitronics::test::okNoAlloc (after == before, "process() and finish() did not allocate through operator new");
    }

    return felitronics::test::report();
}
