// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

// fcore_measure — a tiny streaming measurement CLI for validating felitronics-core's analysis against a
// reference (ffmpeg). Reads interleaved 32-bit-float little-endian PCM (as ffmpeg emits with `-f f32le`)
// and prints ONE number:
//   lufs        → integrated loudness (LUFS)            ↔ ffmpeg ebur128 "I:"
//   truepeak    → max true peak (dBTP, 4× oversampled)  ↔ ffmpeg ebur128 "Peak:" (True Peak)
//   correlation → whole-file Pearson L/R correlation    (synthetic-validated; no clean ffmpeg single number)
//
// Usage: fcore_measure <lufs|truepeak|correlation> <sampleRate> <channels> <raw.f32le>

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;

int main (int argc, char** argv)
{
    if (argc < 5)
    {
        std::fprintf (stderr, "usage: %s <lufs|truepeak|correlation> <sampleRate> <channels> <raw.f32le>\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    const double fs = std::atof (argv[2]);
    const int    nc = std::atoi (argv[3]);
    if (nc < 1 || fs <= 0.0) { std::fprintf (stderr, "bad sampleRate/channels\n"); return 2; }

    std::FILE* f = std::fopen (argv[4], "rb");
    if (! f) { std::perror ("open"); return 2; }

    const int CHUNK = 8192;                          // frames per read
    std::vector<float> inter ((std::size_t) CHUNK * nc);
    std::vector<std::vector<float>> ch ((std::size_t) nc, std::vector<float> ((std::size_t) CHUNK));
    std::vector<const float*> cp ((std::size_t) nc);

    analysis::LoudnessMeter lm;
    std::vector<oversampling::PolyphaseOversampler> os;
    std::vector<float> osbuf;
    double maxTp = 0.0;
    long double sumLR = 0.0L, sumLL = 0.0L, sumRR = 0.0L;

    if (mode == "lufs")      lm.prepare (fs, nc, 4.0 * 3600.0);
    else if (mode == "truepeak") { os.resize ((std::size_t) nc); for (auto& o : os) o.prepare (4, 1, 32); osbuf.assign ((std::size_t) CHUNK * 4, 0.0f); }

    std::size_t got;
    while ((got = std::fread (inter.data(), sizeof (float), (std::size_t) CHUNK * nc, f)) > 0)
    {
        const int frames = (int) (got / (std::size_t) nc);
        for (int i = 0; i < frames; ++i)
            for (int c = 0; c < nc; ++c) ch[(std::size_t) c][(std::size_t) i] = inter[(std::size_t) (i * nc + c)];
        for (int c = 0; c < nc; ++c) cp[(std::size_t) c] = ch[(std::size_t) c].data();

        if (mode == "lufs")
        {
            lm.process (cp.data(), nc, frames);
        }
        else if (mode == "truepeak")
        {
            for (int c = 0; c < nc; ++c)
            {
                const float* ib[1] { ch[(std::size_t) c].data() };
                float* ob[1] { osbuf.data() };
                os[(std::size_t) c].upsample (ib, 1, frames, ob);
                for (int k = 0; k < frames * 4; ++k) maxTp = std::max (maxTp, (double) std::fabs (osbuf[(std::size_t) k]));
            }
        }
        else // correlation
        {
            for (int i = 0; i < frames; ++i)
            {
                const long double l = ch[0][(std::size_t) i];
                const long double r = nc > 1 ? ch[1][(std::size_t) i] : l;
                sumLR += l * r; sumLL += l * l; sumRR += r * r;
            }
        }
    }
    std::fclose (f);

    if (mode == "lufs")            std::printf ("%.2f\n", lm.integratedLufs());
    else if (mode == "truepeak")   std::printf ("%.2f\n", 20.0 * std::log10 (maxTp > 1e-9 ? maxTp : 1e-9));
    else                           { const long double d = std::sqrt (sumLL * sumRR); std::printf ("%.3f\n", d > 1e-12L ? (double) std::clamp (sumLR / d, -1.0L, 1.0L) : 1.0); }

    return 0;
}
