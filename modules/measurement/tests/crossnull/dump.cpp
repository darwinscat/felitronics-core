// Cross-language NULL dumper for felitronics::measurement. Dumps C++ results as raw
// little-endian f64 for an INDEPENDENT numpy recompute to null against. NOT part of ctest
// (keeps core Python-free) — a one-off dev verification per Oleh's "ноль-тестами со scipy/numpy".
#include <felitronics/measurement/Convolve.h>
#include <felitronics/measurement/Sweep.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

static void dump (const std::string& path, const std::vector<double>& v)
{
    std::FILE* f = std::fopen (path.c_str(), "wb");
    std::fwrite (v.data(), sizeof (double), v.size(), f);
    std::fclose (f);
}

int main (int argc, char** argv)
{
    using namespace felitronics::measurement;
    const std::string dir = (argc > 1) ? argv[1] : ".";
    auto p = [&] (const char* n) { return dir + "/" + n; };

    // --- (1) convolve(): FFT-based. NULL vs numpy.convolve (direct time-domain). Two sizes. ---
    for (int c = 0; c < 2; ++c)
    {
        const std::size_t nx = (c == 0) ? 500 : 2000;
        const std::size_t nh = (c == 0) ? 300 : 2000;
        std::vector<double> x (nx), h (nh);
        for (std::size_t i = 0; i < nx; ++i)
            x[i] = std::sin (0.017 * (double) i) + 0.3 * std::cos (0.11 * (double) i + 0.7);
        for (std::size_t i = 0; i < nh; ++i)
            h[i] = std::exp (-0.004 * (double) i) * std::sin (0.05 * (double) i + 0.2);
        const auto y = convolve (x, h);
        dump (p (c == 0 ? "conv0_x.f64" : "conv1_x.f64"), x);
        dump (p (c == 0 ? "conv0_h.f64" : "conv1_h.f64"), h);
        dump (p (c == 0 ? "conv0_y.f64" : "conv1_y.f64"), y);
    }

    // --- (2) magSpectrum(): radix-2 FFT magnitude. NULL vs numpy.fft. ---
    {
        const std::size_t nx = 1000, nfft = 2048;
        std::vector<double> x (nx);
        for (std::size_t i = 0; i < nx; ++i)
            x[i] = std::sin (0.09 * (double) i) + 0.5 * std::sin (0.31 * (double) i + 1.1);
        const auto m = magSpectrum (x, nfft);
        dump (p ("mag_x.f64"), x);
        dump (p ("mag.f64"), m);
    }

    // --- (3) makeSweep(): analytic ESS + fade. NULL vs an independent numpy recompute. ---
    {
        SweepSpec s; s.f1 = 100.0; s.f2 = 4000.0; s.durationSeconds = 0.2;
        s.sampleRate = 16000.0; s.amplitude = 0.5; s.fadeSeconds = 0.01; s.tailSeconds = 0.05;
        const auto sw = makeSweep (s);
        dump (p ("sweep.f64"), sw.signal);
        dump (p ("inverse.f64"), sw.inverse);
        // Emit the SANITIZED spec so numpy recomputes from the SAME truth.
        std::FILE* f = std::fopen (p ("sweep_spec.txt").c_str(), "w");
        std::fprintf (f, "f1 %.17g\nf2 %.17g\ndur %.17g\nsr %.17g\namp %.17g\nfade %.17g\ntail %.17g\nsweepLen %zu\nharmonicL %.17g\n",
                      sw.spec.f1, sw.spec.f2, sw.spec.durationSeconds, sw.spec.sampleRate,
                      sw.spec.amplitude, sw.spec.fadeSeconds, sw.spec.tailSeconds, sw.sweepLen, sw.harmonicL);
        std::fclose (f);
    }
    std::printf ("dumped to %s\n", dir.c_str());
    return 0;
}
