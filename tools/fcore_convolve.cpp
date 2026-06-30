// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_convolve — mono convolution CLI for cross-checking against ffmpeg `afir`. Reads a mono f32le
// input + a mono f32le IR, convolves with the zero-latency PartitionedConvolver, writes a mono f32le
// output the same length as the input.  (Note: the convolver is also unit-tested == direct convolution,
// which is an EXACT reference — afir is itself a partitioned approximation, so this is a sanity cross-check.)
//
// Usage: fcore_convolve <in.f32le> <ir.f32le> <out.f32le>

#include <felitronics/convolution/PartitionedConvolver.h>

#include <cstdio>
#include <vector>

static std::vector<float> readAll (const char* path)
{
    std::vector<float> v;
    std::FILE* f = std::fopen (path, "rb");
    if (! f) return v;
    float b[4096]; std::size_t r;
    while ((r = std::fread (b, sizeof (float), 4096, f)) > 0) v.insert (v.end(), b, b + r);
    std::fclose (f);
    return v;
}

int main (int argc, char** argv)
{
    if (argc < 4) { std::fprintf (stderr, "usage: %s <in.f32le> <ir.f32le> <out.f32le>\n", argv[0]); return 2; }

    std::vector<float> in = readAll (argv[1]);
    std::vector<float> ir = readAll (argv[2]);
    if (in.empty() || ir.empty()) { std::fprintf (stderr, "empty input or IR\n"); return 2; }

    felitronics::convolution::PartitionedConvolver<> conv;
    conv.prepare (128, (int) ir.size());
    conv.setIr (ir.data(), (int) ir.size());

    std::vector<float> out (in.size(), 0.0f);
    conv.process (in.data(), out.data(), (int) in.size());

    std::FILE* f = std::fopen (argv[3], "wb");
    if (! f) { std::perror ("open out"); return 2; }
    std::fwrite (out.data(), sizeof (float), out.size(), f);
    std::fclose (f);
    return 0;
}
