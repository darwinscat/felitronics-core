// Cross-language NULL dumper for felitronics::analysis::offline. Dumps C++ curve/interference results
// as raw little-endian f64 for an INDEPENDENT numpy recompute (check.py). Dev-only, NOT in ctest.
#include <felitronics/analysis/offline/SpectrumCurve.h>
#include <cstdio>
#include <cmath>
#include <span>
#include <string>
#include <vector>

static void dump (const std::string& p, const std::vector<double>& v)
{
    std::FILE* f = std::fopen (p.c_str(), "wb");
    std::fwrite (v.data(), sizeof (double), v.size(), f);
    std::fclose (f);
}

int main (int argc, char** argv)
{
    namespace off = felitronics::analysis::offline;
    const std::string dir = (argc > 1) ? argv[1] : ".";
    auto p = [&] (const char* n) { return dir + "/" + n; };

    // A structured IR: two decaying resonances → a non-trivial curve.
    const double sr = 48000.0;
    const std::size_t N = 8000;
    std::vector<float> ir (N);
    for (std::size_t n = 0; n < N; ++n)
    {
        const double t = (double) n;
        ir[n] = (float) (std::exp (-t / 1500.0) * std::sin (2.0 * M_PI * 1500.0 * t / sr)
                       + 0.5 * std::exp (-t / 400.0) * std::sin (2.0 * M_PI * 5000.0 * t / sr));
    }
    std::vector<double> ird (ir.begin(), ir.end());
    dump (p ("ir.f64"), ird);

    off::LogCurveSpec raw;  raw.normalize = false;
    off::LogCurveSpec norm; norm.normalize = true;
    dump (p ("curve_raw.f64"),  off::logMagnitudeCurve (std::span<const float> (ir), sr, raw));
    dump (p ("curve_norm.f64"), off::logMagnitudeCurve (std::span<const float> (ir), sr, norm));

    std::FILE* fs = std::fopen (p ("spec.txt").c_str(), "w");
    std::fprintf (fs, "sr %.17g\nN %zu\nfLo %.17g\nfHi %.17g\npoints %d\noctaveDivisions %.17g\nminNfft %zu\n",
                  sr, N, raw.fLo, raw.fHi, raw.points, raw.octaveDivisions, (std::size_t) raw.minNfft);
    std::fclose (fs);

    // interferenceDb: two mic curves + a coherent curve (arbitrary dB), dumped for a numpy recompute.
    std::vector<double> micA (64), micB (64), coh (64);
    for (std::size_t i = 0; i < 64; ++i)
    {
        micA[i] = -6.0 + 4.0 * std::sin (0.20 * (double) i);
        micB[i] = -9.0 + 3.0 * std::cos (0.15 * (double) i);
        coh[i]  = -3.0 + 5.0 * std::sin (0.10 * (double) i + 0.5);
    }
    std::vector<off::MicCurveView> mics { { std::span<const double> (micA), true },
                                          { std::span<const double> (micB), true } };
    dump (p ("micA.f64"), micA);
    dump (p ("micB.f64"), micB);
    dump (p ("coh.f64"),  coh);
    dump (p ("interf.f64"), off::interferenceDb (std::span<const double> (coh), std::span<const off::MicCurveView> (mics)));

    std::printf ("dumped to %s\n", dir.c_str());
    return 0;
}
