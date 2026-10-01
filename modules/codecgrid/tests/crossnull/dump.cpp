// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Dev-only: writes what check.py nulls against numpy. For each of the four transforms a programme coded at a
// known offset (SyntheticCodec.h), as interleaved float32, with the curve GridScan reads from it and the reading
// of that curve; and a stretch through BackResampler in both directions. See README.md.

#include <felitronics/codecgrid/BackResampler.h>
#include <felitronics/codecgrid/GridScan.h>

#include "../SyntheticCodec.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace felitronics::codecgrid;

namespace
{
template <class T>
void write (const std::string& path, const std::vector<T>& v)
{
    if (FILE* f = std::fopen (path.c_str(), "wb")) { std::fwrite (v.data(), sizeof (T), v.size(), f); std::fclose (f); }
}
}

int main (int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    constexpr int n = 60000;
    const auto original = synthetic::programme (n, 2026u);
    GridScan<> scanner;
    if (! scanner.prepare (n)) return 1;

    FILE* readings = std::fopen ((dir + "/readings.txt").c_str(), "w");
    if (readings == nullptr) return 1;
    struct Case { Transform t; const char* name; int offset; };
    for (const Case c : { Case { Transform::Mp3, "mp3", 133 }, Case { Transform::AacSine, "aac_sine", 480 }, Case { Transform::AacKbd, "aac_kbd", 1001 },
                          Case { Transform::Celt, "celt", 548 } })
    {
        const auto x = synthetic::coded (original, c.t, c.offset, 0.3);
        std::vector<float> pcm ((std::size_t) 2 * n), curve ((std::size_t) hopOf (c.t) * kCells);
        for (int i = 0; i < n; ++i) { pcm[(std::size_t) (2 * i)] = x.left[(std::size_t) i]; pcm[(std::size_t) (2 * i + 1)] = x.right[(std::size_t) i]; }
        if (! scanner.begin (c.t, x.left.data(), x.right.data(), n, curve.data())) return 1;
        while (! scanner.step (64)) {}
        std::vector<double> d (ReadScratch::doublesFor (hopOf (c.t)));
        std::vector<float> f (ReadScratch::floatsFor (hopOf (c.t)));
        ReadScratch scratch;
        scratch.local = d.data(); scratch.sort = d.data() + hopOf (c.t); scratch.dip = f.data();
        const GridReading r = readCurve (curve.data(), hopOf (c.t), scratch);
        std::vector<float> share (kZeroBands);
        (void) scanner.zeroProfile (r.offset, share.data());
        write (dir + "/" + c.name + ".pcm.f32", pcm);
        write (dir + "/" + c.name + ".curve.f32", curve);
        write (dir + "/" + c.name + ".zeros.f32", share);
        std::fprintf (readings, "%s %d %d %.9g %.9g %.9g %d\n", c.name, c.offset, r.offset, r.score, r.second, r.localDipSum, r.cells);
    }
    std::fclose (readings);

    for (const auto ratio : { BackResampler::ratioFor (44100, 48000), BackResampler::ratioFor (48000, 44100) })
    {
        BackResampler rs;
        if (! rs.prepare (ratio.up, ratio.down)) return 1;
        const int count = (int) rs.outputLength (n);
        std::vector<float> out ((std::size_t) count);
        rs.resample (original.left.data(), n, 0, count, out.data());
        write (dir + "/resample_" + std::to_string (ratio.up) + "_" + std::to_string (ratio.down) + ".f32", out);
    }
    write (dir + "/original_left.f32", original.left);
    std::printf ("dumped to %s\n", dir.c_str());
    return 0;
}
