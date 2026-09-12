// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for analysis::WaveformPeaks, analysis::StereoColumns and StereoSums — ports whose spec is
// executable (the site's audio-peaks.js and stereo-meter.js). The cross-language NULL — the same float32 PCM through
// node and through these classes, on real material — runs out of tree; what is pinned HERE is:
//   * the spec's own numbers on witnesses built from exactly representable values, each computed by the site's
//     JavaScript in node and copied as a literal (an oracle repeats the number; the mutation stand proves the
//     repetition is not idle);
//   * a reference of a DIFFERENT construction (materialise every box first, then bucket; walk each column's
//     stretch directly) against the streaming classes, bit for bit, on seeded programmes;
//   * split invariance, the law-11 call contract, no allocation in process().
//
// THE CONTRACTION WITNESSES AND WHAT THEY CAN SEE — measured, not assumed (P59a mutation stand + probes). The stereo
// band's products are pinned by volatile stores; W5 and W6f read a width that differs by an ulp when `mid += m*m`
// is fused. Whether a compiler fuses depends on the flag, the shape of the source and the call site:
//                                   pinned    pin removed, one expression    pin removed, two statements
//   Apple clang 21 arm64  on O2/O3    ok      FUSED (witnesses fail)          not fused (witnesses cannot see it)
//                         fast O2/O3  ok      not fused                       not fused
//   gcc 14.2 x86-64       on/fast O2  ok      not fused (no FMA without -march)
//                         -march=native: on  ok  FUSED                        not fused;  fast: FUSED
//   gcc 14 arm64 (docker) on O2/O3    ok      FUSED                           not fused;  fast: FUSED
// So this target keeps the tree's `on`: it catches the natural removal (one expression) on the arm64 rows, and the
// pin holds in every cell. The two-statement removal is visible only where gcc-style cross-statement fusion runs with
// FMA available — no in-tree target, CI row or tier builds that (fcore_measure and the wasm modules are -ffp-contract=off);
// the out-of-tree NULL's C++ side, built -ffp-contract=fast -march=native, is what catches it. Witness inputs go through
// a volatile (`opaque`) as a precaution against a compiler evaluating the reduction at compile time.
// The reference below stores its own products, so it stays the spec's arithmetic under any flag.

#include <felitronics_test.h>
#include <felitronics/analysis/StereoColumns.h>
#include <felitronics/analysis/WaveformPeaks.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <new>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new      (std::size_t s, std::align_val_t a)
{
    g_allocs.fetch_add (1, std::memory_order_relaxed);
    return std::aligned_alloc ((std::size_t) a, ((s ? s : 1) + (std::size_t) a - 1) / (std::size_t) a * (std::size_t) a);
}
void* operator new[]    (std::size_t s, std::align_val_t a) { return operator new (s, a); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete   (void* p, std::align_val_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::align_val_t) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t) noexcept { std::free (p); }

using namespace felitronics;
using analysis::PeakMix;
using analysis::StereoColumns;
using analysis::StereoSums;
using analysis::WaveformPeaks;

namespace
{
using Planes = std::vector<std::vector<float>>;

std::uint64_t bits64 (double d)  { std::uint64_t u; std::memcpy (&u, &d, 8); return u; }
std::uint32_t bits32 (float f)   { std::uint32_t u; std::memcpy (&u, &f, 4); return u; }
float         fromBits (std::uint32_t u) { float f; std::memcpy (&f, &u, 4); return f; }

// Witness inputs read back through a volatile, so the compiler cannot evaluate the reduction at compile time. Without
// this the contraction witnesses are dead: clang at -O3 constant-folds needle() over a `const float[2]` with exact,
// unfused arithmetic, and a header with its pin removed passed them (measured on the mutation stand).
std::vector<float> opaque (std::initializer_list<float> xs)
{
    std::vector<float> v;
    for (float x : xs) { volatile float y = x; v.push_back ((float) y); }
    return v;
}

std::vector<const float*> ptrs (const Planes& p, std::size_t off = 0)
{
    std::vector<const float*> v;
    for (const auto& c : p) v.push_back (c.data() + off);
    return v;
}

struct PeaksOut { std::vector<double> peaks; int emitted = 0; bool ok = false; };

// Feeds the planes in the given split (a list of call lengths, repeated until the file is done).
PeaksOut runPeaks (const Planes& p, double sr, int buckets, PeakMix mix, const std::vector<int>& split)
{
    PeaksOut r;
    WaveformPeaks w;
    const auto frames = (std::uint64_t) p[0].size();
    if (! w.prepare (sr, (int) p.size(), frames, buckets, mix)) return r;
    std::uint64_t at = 0; std::size_t k = 0;
    while (at < frames)
    {
        const int n = (int) std::min<std::uint64_t> ((std::uint64_t) split[k++ % split.size()], frames - at);
        const auto v = ptrs (p, (std::size_t) at);
        if (! w.process (v.data(), (int) p.size(), n)) return r;
        at += (std::uint64_t) n;
    }
    r.ok = w.complete();
    r.peaks.assign (w.peaks().begin(), w.peaks().end());
    r.emitted = w.bucketsEmitted();
    return r;
}

struct StereoOut { std::vector<float> width, corr, loud; double maxLoud = 0.0; int cols = 0; bool mono = false, ok = false; };

StereoOut runStereo (const Planes& p, int columns, const std::vector<int>& split)
{
    StereoOut r;
    StereoColumns s;
    const auto frames = (std::uint64_t) p[0].size();
    if (! s.prepare ((int) p.size(), frames, columns)) return r;
    std::uint64_t at = 0; std::size_t k = 0;
    while (at < frames)
    {
        const int n = (int) std::min<std::uint64_t> ((std::uint64_t) split[k++ % split.size()], frames - at);
        const auto v = ptrs (p, (std::size_t) at);
        if (! s.process (v.data(), (int) p.size(), n)) return r;
        at += (std::uint64_t) n;
    }
    r.ok = s.complete();
    r.width.assign (s.width().begin(), s.width().end());
    r.corr.assign (s.correlation().begin(), s.correlation().end());
    r.loud.assign (s.rms().begin(), s.rms().end());
    r.maxLoud = s.maxRms(); r.cols = s.columns(); r.mono = s.isMono();
    return r;
}

bool samePeaks (const PeaksOut& a, const PeaksOut& b)
{
    if (a.ok != b.ok || a.emitted != b.emitted || a.peaks.size() != b.peaks.size()) return false;
    for (std::size_t i = 0; i < a.peaks.size(); ++i) if (bits64 (a.peaks[i]) != bits64 (b.peaks[i])) return false;
    return true;
}

bool sameF (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (bits32 (a[i]) != bits32 (b[i])) return false;
    return true;
}

bool sameStereo (const StereoOut& a, const StereoOut& b)
{
    return a.ok == b.ok && a.cols == b.cols && a.mono == b.mono && bits64 (a.maxLoud) == bits64 (b.maxLoud)
        && sameF (a.width, b.width) && sameF (a.corr, b.corr) && sameF (a.loud, b.loud);
}

//------------------------------------------------------------------------------------------------------------
// THE REFERENCE, OF A DIFFERENT CONSTRUCTION. Peaks: every box mean is materialised into an array first, and the
// buckets are then cut from that array — no running state shared with the frame loop. Stereo: each column's
// stretch is walked directly from its two ends, the way the spec itself walks it, rather than streamed.
// Both are the SPEC'S arithmetic (it is the thing being ported), written as a second, independent program.
PeaksOut refPeaks (const Planes& p, double sr, int buckets, PeakMix mix)
{
    PeaksOut r;
    const std::size_t nch = p.size(), len = p[0].size();
    double decimD = std::floor (sr / 8000.0);
    if (sr / 8000.0 - decimD >= 0.5) decimD += 1.0;
    const std::size_t decim = (std::size_t) std::max (1.0, decimD);

    std::vector<double> boxes;
    for (std::size_t b = 0; (b + 1) * decim <= len; ++b)
    {
        double acc = 0.0;
        for (std::size_t j = b * decim; j < (b + 1) * decim; ++j)
        {
            double fv = 0.0;
            if (mix == PeakMix::Left)       fv = p[0][j];
            else if (mix == PeakMix::Right) fv = p[nch - 1][j];
            else if (mix == PeakMix::Max)   { for (std::size_t c = 0; c < nch; ++c) { const double x = std::fabs ((double) p[c][j]); if (x > fv) fv = x; } }
            else                            { double s = 0.0; for (std::size_t c = 0; c < nch; ++c) s += (double) p[c][j]; fv = s / (double) nch; }
            acc += fv;
        }
        const double mean = acc / (double) decim;
        boxes.push_back (mix == PeakMix::Max ? mean : std::fabs (mean));
    }

    const double decLen = std::max ((double) buckets, std::floor ((double) len / (double) decim));
    const double per = decLen / (double) buckets;
    r.peaks.assign ((std::size_t) buckets, 0.0);
    std::size_t bucket = 0; double mx = 0.0;
    for (std::size_t b = 0; b < boxes.size(); ++b)
    {
        if (boxes[b] > mx) mx = boxes[b];
        if ((double) (b + 1) >= (double) (bucket + 1) * per && bucket < (std::size_t) buckets)
        {
            r.peaks[bucket++] = mx;
            mx = 0.0;
        }
    }
    r.emitted = (int) bucket;
    r.ok = true;
    return r;
}

StereoOut refStereo (const Planes& p, int columns)
{
    StereoOut r;
    const std::size_t len = p[0].size();
    const std::vector<float>& L = p[0];
    const std::vector<float>& R = p.size() > 1 ? p[1] : p[0];
    r.cols = (int) std::max<std::size_t> (1, std::min<std::size_t> ((std::size_t) columns, len));
    r.mono = p.size() < 2;
    const double per = (double) len / (double) r.cols;
    for (int i = 0; i < r.cols; ++i)
    {
        const auto start = (std::size_t) std::floor ((double) i * per);
        const auto end   = std::min (len, (std::size_t) std::floor ((double) (i + 1) * per));
        double sLL = 0, sRR = 0, sLR = 0, sM = 0, sS = 0;
        for (std::size_t j = start; j < end; ++j)
        {
            const double l = L[j], r2 = R[j];
            volatile double q0 = l * l, q1 = r2 * r2, q2 = l * r2;
            const double m = (l + r2) * 0.5, s = (l - r2) * 0.5;
            volatile double q3 = m * m, q4 = s * s;
            sLL += q0; sRR += q1; sLR += q2; sM += q3; sS += q4;
        }
        const double n = (double) (end > start ? end - start : 1);
        const double mr = std::sqrt (sM / n), sr2 = std::sqrt (sS / n);
        r.width.push_back ((float) ((mr + sr2) > 1e-9 ? sr2 / (mr + sr2) : 0.0));
        const double den = std::sqrt (sLL * sRR);
        double c = 1.0;
        if (den > 1e-12) { c = sLR / den; c = c > 1.0 ? 1.0 : c; c = c < -1.0 ? -1.0 : c; }
        r.corr.push_back ((float) c);
        const double rms = std::sqrt ((sLL + sRR) / (2.0 * n));
        r.loud.push_back ((float) rms);
        if (rms > r.maxLoud) r.maxLoud = rms;
    }
    r.ok = true;
    return r;
}

// A seeded programme. Integer PRNG, values drawn from the set a decoder actually produces (16-bit and 24-bit grids,
// raw float32), with silence stretches, full-scale samples, a DC run and an anti-phase run — nothing transcendental.
Planes programme (std::uint32_t seed, int nch, std::size_t len)
{
    std::uint64_t st = 0x9E3779B97F4A7C15ull ^ seed;
    auto next = [&] { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; };
    Planes p ((std::size_t) nch, std::vector<float> (len));
    for (std::size_t j = 0; j < len; ++j)
    {
        const std::uint64_t kind = (j / 997) % 5;
        for (int c = 0; c < nch; ++c)
        {
            const std::uint64_t u = next();
            float v;
            switch (u % 3)
            {
                case 0:  v = (float) ((std::int16_t) (u >> 16)) / 32768.0f; break;
                case 1:  v = (float) ((std::int32_t) ((u >> 8) & 0xFFFFFF) - 0x800000) / 8388608.0f; break;
                default: { std::uint32_t w = (std::uint32_t) (u >> 32); v = (float) ((double) w / 4294967296.0 * 2.0 - 1.0); }
            }
            if (kind == 1) v = 0.0f;                                  // silence: the guards' territory
            if (kind == 2 && (u >> 40) % 11 == 0) v = (c % 2) ? -1.0f : 1.0f;
            if (kind == 3) v = 0.25f;                                 // DC
            if (kind == 4 && c == 1) v = -p[0][j];                    // anti-phase
            p[(std::size_t) c][j] = v;
        }
    }
    return p;
}

Planes mono (std::vector<float> x) { for (auto& v : x) { volatile float y = v; v = y; } return Planes { std::move (x) }; }
} // namespace

int main()
{
    std::printf ("felitronics::analysis waveform shape tests (WaveformPeaks, StereoColumns)\n");

    //==========================================================================================================
    test::group ("spec numbers: the site's JavaScript, computed in node, copied as literals");
    {
        // W1 — THE LOST LAST BUCKET. 2007 frames at 8 kHz (decim 1), 1000 buckets: the last boundary is
        // (2007/1000)·1000 = 2007.0000000000002, the box counter tops out at 2007, bucket 999 is never emitted —
        // so the one non-zero sample, in frame 2006, is nowhere in the output. node: 0 non-zero buckets.
        Planes x = mono (std::vector<float> (2007, 0.0f)); x[0][2006] = 1.0f;
        const PeaksOut w = runPeaks (x, 8000.0, 1000, PeakMix::Left, { 2007 });
        bool allZero = true; for (double v : w.peaks) allZero = allZero && bits64 (v) == 0;
        test::ok (w.ok && allZero && w.emitted == 999, "W1: 2007 frames, 1000 buckets -> the last bucket is never emitted (999 emitted, all 0)");

        // W2 — THE COLUMN BOUNDARY IS floor(i · (len/cols)). 1206 frames, 1200 columns, one impulse at frame 200:
        // the spec's column 200 starts at 200 (the exact partition says 201). node: loud[199] = 0,
        // loud[200] = 0.7071067690849304 (0x3f3504f3), maxLoud = 0x3fe6a09e667f3bcd.
        Planes y = mono (std::vector<float> (1206, 0.0f)); y[0][200] = 1.0f;
        const StereoOut s = runStereo (y, 1200, { 1206 });
        test::ok (s.ok && s.cols == 1200 && bits32 (s.loud[199]) == 0u && bits32 (s.loud[200]) == 0x3f3504f3u,
                  "W2: 1206 / 1200 — column 200 holds the impulse at frame 200 (loud 0x3f3504f3), column 199 is silent");
        test::ok (bits64 (s.maxLoud) == 0x3fe6a09e667f3bcdull, "W2: maxLoud is the unrounded double 0x3fe6a09e667f3bcd");

        // W3 — JavaScript's Math.round: 20 kHz / 8000 = 2.5 -> 3 (a banker's rint would say 2 and read 0).
        // 12 kHz / 8000 = 1.5 -> 2. node: 0x3fd5555555555555 and 0.5.
        test::ok (bits64 (runPeaks (mono ({ 0, 0, 1 }), 20000.0, 1, PeakMix::Left, { 3 }).peaks[0]) == 0x3fd5555555555555ull,
                  "W3: 20 kHz -> decim 3, [0,0,1] reads 1/3 (0x3fd5555555555555)");
        test::ok (runPeaks (mono ({ 0, 1, 0, 0 }), 12000.0, 1, PeakMix::Left, { 4 }).peaks[0] == 0.5,
                  "W3: 12 kHz -> decim 2, [0,1,0,0] reads 0.5");
        int d = 0;
        test::ok (WaveformPeaks::decimationFor (44100.0, d) && d == 6, "decim(44.1 kHz) = 6 (an envelope at 7350 Hz)");
        test::ok (WaveformPeaks::decimationFor (1000.0, d) && d == 1, "decim(1 kHz) = max(1, round(0.125)) = 1");

        // W4 — THE BOX ACCUMULATOR IS SEQUENTIAL ACROSS CALLS. [1, 2^-53, 2^-53] at 24 kHz (decim 3): each addition
        // rounds back to 1, mean 1/3 = 0x3fd5555555555555. Summing a call's frames separately reads ...557.
        for (const auto& split : std::vector<std::vector<int>> { { 3 }, { 1, 2 }, { 2, 1 }, { 1 } })
            test::ok (bits64 (runPeaks (mono ({ 1.0f, 0x1p-53f, 0x1p-53f }), 24000.0, 1, PeakMix::Left, split).peaks[0]) == 0x3fd5555555555555ull,
                      "W4: [1, 2^-53, 2^-53] reads 0x3fd5555555555555 in a split of " + std::to_string (split[0]));

        // W5 — THE CONTRACTION WITNESS. L = [1+2^-23, 2^-30], R = [2^-25, 1+2^-23]: widthOf reads 0x3fdffffff7c00012;
        // `mid += m*m` fused into one FMA reads ...0010. (A product of two float32 samples is exact in double, so
        // l*l cannot tell; m*m can.) node: widthOf 0.49999999231658976.
        const auto L5 = opaque ({ 1.0f + 0x1p-23f, 0x1p-30f }), R5 = opaque ({ 0x1p-25f, 1.0f + 0x1p-23f });
        StereoColumns::Needle nd;
        test::ok (StereoColumns::needle (L5.data(), R5.data(), 2, 0, 2, nd) && bits64 (nd.width) == 0x3fdffffff7c00012ull,
                  "W5: the needle width is 0x3fdffffff7c00012 — the products are not fused");
        // W6f — the NULL's own catch: six frames of a seeded 3-channel programme at 8 kHz (item syn-8000-3-6), L and R
        // the first two channels. node: width 0x3fe34195a140630a, correlation 0xbfdfe94ee90f9bef, RMS 0x3fd5bb4cf415e0bb;
        // fused, the width reads ...6309.
        const auto L6 = opaque ({ fromBits (0x33c264cdu), fromBits (0xbf59de04u), fromBits (0xbebc9000u), fromBits (0x3f027c00u), fromBits (0x3d26d180u), fromBits (0x34097a6eu) });
        const auto R6 = opaque ({ fromBits (0x35442d9eu), fromBits (0x34d15606u), fromBits (0x3e95b12cu), fromBits (0xbead613cu), fromBits (0x3e78e800u), fromBits (0x3d5fc000u) });
        test::ok (StereoColumns::needle (L6.data(), R6.data(), 6, 0, 6, nd) && bits64 (nd.width) == 0x3fe34195a140630aull
               && bits64 (nd.correlation) == 0xbfdfe94ee90f9befull && bits64 (nd.rms) == 0x3fd5bb4cf415e0bbull,
                  "W6f: syn-8000-3-6 needle — width 0x3fe34195a140630a, corr 0xbfdfe94ee90f9bef, RMS 0x3fd5bb4cf415e0bb");

        // W6 — THE MIX MODES on three channels [-1,1] [0.25,0.25] [0.5,0.5] at 16 kHz, one bucket.
        // node: L 0, R 0.5 (the LAST channel), avr 0.25, max 1 (abs before the channel max). Stereo with five
        // columns asked: 2 columns (len), width [0.625, 0.375], corr [-1, 1] — channels 0 and 1.
        Planes six3 { { -1.0f, 1.0f }, { 0.25f, 0.25f }, { 0.5f, 0.5f } };
        test::ok (runPeaks (six3, 16000.0, 1, PeakMix::Left,    { 2 }).peaks[0] == 0.0,  "W6: L = 0");
        test::ok (runPeaks (six3, 16000.0, 1, PeakMix::Right,   { 2 }).peaks[0] == 0.5,  "W6: R = 0.5 (the last channel, not channel 1)");
        test::ok (runPeaks (six3, 16000.0, 1, PeakMix::Average, { 2 }).peaks[0] == 0.25, "W6: avr = 0.25 (all channels)");
        test::ok (runPeaks (six3, 16000.0, 1, PeakMix::Max,     { 2 }).peaks[0] == 1.0,  "W6: max = 1 (|channel| before the max)");
        const StereoOut s6 = runStereo (six3, 5, { 2 });
        test::ok (s6.cols == 2 && ! s6.mono && s6.width[0] == 0.625f && s6.width[1] == 0.375f && s6.corr[0] == -1.0f && s6.corr[1] == 1.0f,
                  "W6: stereo reads channels 0 and 1 — 2 columns, width [0.625, 0.375], corr [-1, 1]");
        Planes six6; for (int c = 0; c < 6; ++c) six6.push_back ({ (float) c / 8.0f, (float) c / 8.0f });
        test::ok (runPeaks (six6, 8000.0, 1, PeakMix::Right, { 2 }).peaks[0] == 0.625 && runPeaks (six6, 8000.0, 1, PeakMix::Average, { 2 }).peaks[0] == 0.3125,
                  "W16: six channels c/8 — R reads channel 5 (0.625), avr the mean of all six (0.3125)");

        // W7 — NON-FINITE. A NaN box never wins `d > max` (node: [NaN,0,0.5,0.5] at 16 kHz reads 0.5); +Inf wins
        // (reads Infinity). Stereo L = +Inf, R = 1: corr NaN (0x7fc00000), width NaN, loud +Inf, maxLoud +Inf.
        test::ok (runPeaks (mono ({ NAN, 0.0f, 0.5f, 0.5f }), 16000.0, 1, PeakMix::Left, { 4 }).peaks[0] == 0.5, "W7: a NaN box is skipped");
        test::ok (std::isinf (runPeaks (mono ({ INFINITY, 0.0f, 0.5f, 0.5f }), 16000.0, 1, PeakMix::Left, { 4 }).peaks[0]), "W7: +Inf wins");
        // ... and a NaN box LAST in its bucket: [0.5,0.5,NaN,0] still reads 0.5 (node). A max written as std::max(d, max)
        // returns the NaN here — with the NaN box first, as above, it would not have shown.
        test::ok (runPeaks (mono ({ 0.5f, 0.5f, NAN, 0.0f }), 16000.0, 1, PeakMix::Left, { 4 }).peaks[0] == 0.5, "W7b: a NaN box after the max does not replace it");
        const StereoOut s7 = runStereo (Planes { { INFINITY }, { 1.0f } }, 1, { 1 });
        test::ok (std::isnan (s7.corr[0]) && std::isnan (s7.width[0]) && std::isinf (s7.loud[0]) && std::isinf (s7.maxLoud),
                  "W7: stereo L=+Inf R=1 -> corr NaN, width NaN, loud +Inf, maxLoud +Inf (the clamp lets NaN through)");
        const StereoOut s7n = runStereo (Planes { { NAN }, { 1.0f } }, 1, { 1 });
        test::ok (s7n.width[0] == 0.0f && s7n.corr[0] == 1.0f && std::isnan (s7n.loud[0]) && s7n.maxLoud == 0.0,
                  "W7: stereo L=NaN -> width 0, corr +1, loud NaN, and NaN never becomes maxLoud");

        // W8 — SHORT FILES. [0.1 0.3 0.5 0.7 0.2 0.4] at 16 kHz, 5 buckets: 3 boxes < 5 buckets, so the trailing
        // buckets stay 0. node: [0.20000000670552254, 0.5999999940395355, 0.30000000447034836, 0, 0].
        const PeaksOut w8 = runPeaks (mono ({ 0.1f, 0.3f, 0.5f, 0.7f, 0.2f, 0.4f }), 16000.0, 5, PeakMix::Left, { 6 });
        test::ok (w8.peaks[0] == 0.20000000670552254 && w8.peaks[1] == 0.5999999940395355 && w8.peaks[2] == 0.30000000447034836
               && w8.peaks[3] == 0.0 && w8.peaks[4] == 0.0 && w8.emitted == 3, "W8: 3 boxes into 5 buckets -> the last two stay 0");
        const StereoOut s13 = runStereo (mono ({ 1.0f, 0.5f, 0.25f }), 1100, { 3 });
        test::ok (s13.cols == 3 && s13.mono && s13.width == std::vector<float> (3, 0.0f) && s13.corr == std::vector<float> (3, 1.0f),
                  "W13: 3 frames, 1100 columns asked -> 3 columns, mono: width 0, corr +1");

        // W9/W10 — PROMOTION AND THE TWO OUTPUT TYPES. [1, 1+2^-23] at 16 kHz: the double peak is
        // 1.0000000596046448 (0x3ff0000010000000) and its float32 form is 1 (a tie, to even); the column loudness
        // is float32 0x3f800001 while maxLoud is the unrounded 0x3ff0000010000008.
        WaveformPeaks w10;
        const auto x10 = opaque ({ 1.0f, 1.0f + 0x1p-23f });
        const float* p10[1] { x10.data() };
        test::run (w10.prepare (16000.0, 1, 2, 1, PeakMix::Left));
        test::run (w10.process (p10, 1, 2));
        test::ok (bits64 (w10.peaks()[0]) == 0x3ff0000010000000ull && bits32 (w10.peakAsFloat32 (0)) == 0x3f800000u,
                  "W10: peak 0x3ff0000010000000 as a double, 0x3f800000 as float32 (the tie rounds to even)");
        const StereoOut s10 = runStereo (mono ({ 1.0f, 1.0f + 0x1p-23f }), 1, { 2 });
        test::ok (bits32 (s10.loud[0]) == 0x3f800001u && bits64 (s10.maxLoud) == 0x3ff0000010000008ull,
                  "W10: loud[0] 0x3f800001, maxLoud 0x3ff0000010000008 — not recomputed from the rounded array");
        const StereoOut s9 = runStereo (mono ({ 1.0f + 0x1p-23f }), 1, { 1 });
        test::ok (bits32 (s9.loud[0]) == 0x3f800001u, "W9: loud of [1+2^-23] is 0x3f800001 — the sample is promoted before l*l");

        // W12 — SUBNORMALS AND THE SIGN OF ZERO. L = [2, 0], R = [-2^-149, 2]: the needle correlation is -2^-150
        // (0xb690000000000000), the float32 column rounds it to -0 (0x80000000).
        const auto L12 = opaque ({ 2.0f, 0.0f }), R12 = opaque ({ -0x1p-149f, 2.0f });
        test::ok (StereoColumns::needle (L12.data(), R12.data(), 2, 0, 2, nd) && bits64 (nd.correlation) == 0xb690000000000000ull,
                  "W12: needle correlation -2^-150 (0xb690000000000000)");
        const StereoOut s12 = runStereo (Planes { { 2.0f, 0.0f }, { -0x1p-149f, 2.0f } }, 1, { 2 });
        test::ok (bits32 (s12.corr[0]) == 0x80000000u, "W12: the float32 column is -0 (0x80000000)");

        // W14 — NOT PEARSON: L = (1, 2), R = (2, 1) reads 0.8 (Pearson's coefficient would be -1).
        const auto L14 = opaque ({ 1.0f, 2.0f }), R14 = opaque ({ 2.0f, 1.0f });
        test::ok (StereoColumns::needle (L14.data(), R14.data(), 2, 0, 2, nd) && nd.correlation == 0.8, "W14: uncentred correlation of (1,2),(2,1) is 0.8");

        // W15 — 44.1 kHz decim 6 on [1,-1,0,0,0,0, 0.5,0,0,0,0,0]: box 1 mean 0, box 2 mean 1/12. node: 0x3fb5555555555555.
        test::ok (bits64 (runPeaks (mono ({ 1, -1, 0, 0, 0, 0, 0.5f, 0, 0, 0, 0, 0 }), 44100.0, 1, PeakMix::Left, { 12 }).peaks[0]) == 0x3fb5555555555555ull,
                  "W15: 44.1 kHz boxes of 6 — 1/12 (0x3fb5555555555555)");

        // The crew's witnesses (P59a pre-start round), each re-computed by the site's JS in node.
        // W17 — ONE DIVISION PER FRAME, IN CHANNEL ORDER. Three channels at 16 kHz (decim 2), frames
        // [3f35e00d bf089722 3eb88b0a] and [3ef22e5b bf11b04d bf347cf7]: 'avr' reads 0x3fa6828cd5555556. Summing
        // all six and dividing once reads ...555 (node); two channels divide by 2 exactly and could never tell.
        Planes three { { fromBits (0x3f35e00du), fromBits (0x3ef22e5bu) }, { fromBits (0xbf089722u), fromBits (0xbf11b04du) },
                       { fromBits (0x3eb88b0au), fromBits (0xbf347cf7u) } };
        test::ok (bits64 (runPeaks (three, 16000.0, 1, PeakMix::Average, { 2 }).peaks[0]) == 0x3fa6828cd5555556ull,
                  "W17: three-channel avr divides per frame — 0x3fa6828cd5555556");

        // W18/W19 — BOUNDARIES ARE MULTIPLIED, NEVER ACCUMULATED. 44.1 kHz, 74070 frames (decLen 12345), box 7407
        // (0-based) at 0.9, the rest 0.1: 600·12.345 = 7407 exactly, so bucket 599 is emitted after box 7406 and box
        // 7407 opens bucket 600 — 0.1 / 0.9. A running sum (7407.000000000039) emits bucket 599 one box later, after
        // box 7407, and swaps them: 0.9 / 0.1 (node, both). Box 7408 would read the same under both, which is what this
        // witness used to use — the review round found it could not tell. Stereo: 44106 frames, 1200 columns, one
        // impulse at 14702: 400·per = 14702.000000000002, so column 400 starts at 14702 and holds it, with 37 frames
        // where a running sum's column would have 36 (loud 0x3e2aaaab).
        std::vector<float> x18 (74070, 0.1f);
        for (std::size_t j = 7407 * 6; j < 7408 * 6; ++j) x18[j] = 0.9f;
        const PeaksOut w18 = runPeaks (mono (x18), 44100.0, 1000, PeakMix::Left, { 74070 });
        test::ok (bits64 (w18.peaks[599]) == 0x3fb99999a0000000ull && bits64 (w18.peaks[600]) == 0x3fecccccc0000000ull,
                  "W18: bucket 599 = 0.1, bucket 600 = 0.9 (the boundary is 600 x perBucket, not a running sum)");
        std::vector<float> x19 (44106, 0.0f); x19[14702] = 1.0f;
        const StereoOut s19 = runStereo (mono (x19), 1200, { 44106 });
        test::ok (bits32 (s19.loud[399]) == 0u && bits32 (s19.loud[400]) == 0x3e2aaaabu,
                  "W19: 44106 frames, the impulse at 14702 is in column 400 (0x3e2aaaab), not 399");

        // W20 — THE LAST FRAME CAN BELONG TO NO COLUMN. 1206 frames / 1200 columns and 1120 / 1100: floor(cols·per)
        // is len - 1, so an impulse in the last frame is read by nothing — maxLoud 0.
        std::vector<float> x20 (1206, 0.0f); x20[1205] = 1.0f;
        std::vector<float> x20b (1120, 0.0f); x20b[1119] = 1.0f;
        test::ok (runStereo (mono (x20), 1200, { 1206 }).maxLoud == 0.0 && runStereo (mono (x20b), 1100, { 1120 }).maxLoud == 0.0,
                  "W20: 1206/1200 and 1120/1100 — a last-frame impulse is in no column (maxLoud 0)");

        // W21 — SUBNORMALS ARE MEASURED, NOT FLUSHED: L = R = [1e-40f] reads maxLoud 0x37a16c2000000000 (loud
        // 0x000116c2); a 1.4e-45 sample peaks at 0x36a0000000000000. A host running this under FTZ/DAZ reads 0.
        const StereoOut s21 = runStereo (Planes { { 1e-40f }, { 1e-40f } }, 1, { 1 });
        test::ok (bits64 (s21.maxLoud) == 0x37a16c2000000000ull && bits32 (s21.loud[0]) == 0x000116c2u, "W21: a subnormal pair's loudness is 0x37a16c2000000000");
        test::ok (bits64 (runPeaks (mono ({ 1.4e-45f }), 8000.0, 1, PeakMix::Average, { 1 }).peaks[0]) == 0x36a0000000000000ull, "W21: a subnormal peak is 0x36a0000000000000");

        // W23 — A SECOND CONTRACTION WITNESS, from random-looking float32: L = [bdcdafb0 3f5b9584], R = [bdd5e88f
        // bcd01fac]: widthOf 0x3fe03d22ac2492b1 (fused: ...b0). Its float32 column does not show it — the double
        // needle is where contraction is visible.
        const auto L23 = opaque ({ fromBits (0xbdcdafb0u), fromBits (0x3f5b9584u) }), R23 = opaque ({ fromBits (0xbdd5e88fu), fromBits (0xbcd01facu) });
        test::ok (StereoColumns::needle (L23.data(), R23.data(), 2, 0, 2, nd) && bits64 (nd.width) == 0x3fe03d22ac2492b1ull, "W23: needle width 0x3fe03d22ac2492b1");

        // W25 — THE WIDTH GUARD IS 1e-9, NOT 0: L = [1e-10, 1e-10], R = -L has no mid and a side RMS of 1e-10, under
        // the guard, so widthOf reads 0 (node) — a guard of 0 would read 1.
        const auto L25 = opaque ({ 1e-10f, 1e-10f }), R25 = opaque ({ -1e-10f, -1e-10f });
        test::ok (StereoColumns::needle (L25.data(), R25.data(), 2, 0, 2, nd) && nd.width == 0.0, "W25: a side RMS of 1e-10 is under the 1e-9 guard — width 0");

        // W24 — THE TAIL IS NOT FLUSHED. 44.1 kHz, [0.1 x6, 0.9], 2 buckets: the 0.9 never completes a box, and
        // there is no finish() to push it: [0x3fb99999a0000000, 0].
        const PeaksOut w24 = runPeaks (mono ({ 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.9f }), 44100.0, 2, PeakMix::Left, { 1 });
        test::ok (bits64 (w24.peaks[0]) == 0x3fb99999a0000000ull && w24.peaks[1] == 0.0, "W24: a partial last box is dropped, not flushed");
    }

    //==========================================================================================================
    test::group ("reference of a different construction, bit for bit, on seeded programmes");
    {
        const double rates[] { 8000.0, 12000.0, 16000.0, 20000.0, 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0, 1000.0 };
        const std::size_t lens[] { 1, 2, 5, 6, 7, 999, 1000, 1001, 2007, 6000, 6006, 12047, 48000, 100003 };
        const int bucketsList[] { 1, 3, 1000, 1100 };
        const int columnsList[] { 1, 7, 1100, 1200 };
        const int chans[] { 1, 2, 3, 6 };
        int peakCases = 0, peakBad = 0, stereoCases = 0, stereoBad = 0;
        std::uint32_t seed = 1;
        for (std::size_t len : lens)
            for (int nch : chans)
            {
                const Planes p = programme (seed++, nch, len);
                for (double sr : rates)
                    for (int b : bucketsList)
                        for (PeakMix mix : { PeakMix::Average, PeakMix::Left, PeakMix::Right, PeakMix::Max })
                        {
                            ++peakCases;
                            if (! samePeaks (runPeaks (p, sr, b, mix, { (int) len }), refPeaks (p, sr, b, mix))) ++peakBad;
                        }
                for (int cols : columnsList)
                {
                    ++stereoCases;
                    if (! sameStereo (runStereo (p, cols, { (int) len }), refStereo (p, cols))) ++stereoBad;
                }
            }
        test::ok (peakBad == 0,   "peaks: " + std::to_string (peakCases) + " cases, " + std::to_string (peakBad) + " differ from the reference");
        test::ok (stereoBad == 0, "stereo: " + std::to_string (stereoCases) + " cases, " + std::to_string (stereoBad) + " differ from the reference");

        // The needle IS a column: over a column's own stretch its double numbers round to that column's float32.
        const Planes p = programme (777, 2, 12047);
        StereoColumns sc;
        test::run (sc.prepare (2, 12047, 1100));
        const auto v = ptrs (p);
        test::run (sc.process (v.data(), 2, 12047));
        int needleBad = 0;
        for (int i = 0; i < sc.columns(); ++i)
        {
            StereoColumns::Needle nd;
            if (! StereoColumns::needle (p[0].data(), p[1].data(), 12047, sc.columnStart (i), sc.columnEnd (i), nd)
                || bits32 ((float) nd.width) != bits32 (sc.width()[(std::size_t) i])
                || bits32 ((float) nd.correlation) != bits32 (sc.correlation()[(std::size_t) i])
                || bits32 ((float) nd.rms) != bits32 (sc.rms()[(std::size_t) i])) ++needleBad;
        }
        test::ok (needleBad == 0, "the needle over each column's stretch rounds to that column, all 1100 columns");
    }

    //==========================================================================================================
    test::group ("a split changes nothing — the Probe's set, and splits cut exactly on bucket and column edges");
    {
        const std::size_t len = 250007;
        const Planes p = programme (4242, 2, len);
        const PeaksOut  wholeP = runPeaks (p, 44100.0, 1100, PeakMix::Average, { (int) len });
        const StereoOut wholeS = runStereo (p, 1100, { (int) len });
        for (int chunk : { 1, 7, 999, 4096, 8192, 8193, 100003 })
        {
            test::ok (samePeaks (runPeaks (p, 44100.0, 1100, PeakMix::Average, { chunk }), wholeP), "peaks: chunk " + std::to_string (chunk));
            test::ok (sameStereo (runStereo (p, 1100, { chunk }), wholeS), "stereo: chunk " + std::to_string (chunk));
        }
        // Edges ON the features: a split whose every call ends exactly where a bucket is emitted (a box boundary
        // times the bucket's box count) and one ending exactly at every column end, and both one frame off.
        WaveformPeaks probe;
        test::run (probe.prepare (44100.0, 2, len, 1100, PeakMix::Average));
        const int decim = probe.decimation();
        StereoColumns cprobe;
        test::run (cprobe.prepare (2, len, 1100));
        std::vector<std::size_t> bucketEdges, columnEdges;
        {
            const double decLen = std::max (1100.0, std::floor ((double) len / decim));
            const double per = decLen / 1100.0;
            std::size_t bucket = 0;
            for (std::size_t box = 1; box * (std::size_t) decim <= len; ++box)
                if ((double) box >= (double) (bucket + 1) * per) { bucketEdges.push_back (box * (std::size_t) decim); ++bucket; }
            for (int i = 0; i < cprobe.columns(); ++i) columnEdges.push_back ((std::size_t) cprobe.columnEnd (i));
        }
        auto splitAt = [&] (const std::vector<std::size_t>& edges, long shift)
        {
            std::vector<int> s; std::size_t prev = 0;
            for (std::size_t e : edges)
            {
                const long at = std::clamp ((long) e + shift, (long) prev, (long) len);
                if ((std::size_t) at > prev) { s.push_back ((int) ((std::size_t) at - prev)); prev = (std::size_t) at; }
            }
            if (prev < len) s.push_back ((int) (len - prev));
            return s;
        };
        for (long shift : { 0L, -1L, 1L })
        {
            const auto sb = splitAt (bucketEdges, shift);
            const auto scs = splitAt (columnEdges, shift);
            // a split list is consumed cyclically by runPeaks; these cover the file exactly once
            test::ok (samePeaks (runPeaks (p, 44100.0, 1100, PeakMix::Average, sb), wholeP),
                      "peaks: calls ending on every bucket edge, shifted " + std::to_string (shift));
            test::ok (sameStereo (runStereo (p, 1100, scs), wholeS),
                      "stereo: calls ending on every column edge, shifted " + std::to_string (shift));
        }
    }

    //==========================================================================================================
    test::group ("the call contract (law 11): prepare refuses, width is exact, a call past the file is refused whole");
    {
        WaveformPeaks w;
        test::ok (! w.prepare (0.0, 1, 10) && ! w.prepare (NAN, 1, 10) && ! w.prepare (INFINITY, 1, 10) && ! w.prepare (-8000.0, 1, 10),
                  "peaks: a rate that is not a rate is refused");
        test::ok (! w.prepare (48000.0, 0, 10) && ! w.prepare (48000.0, 1, 0) && ! w.prepare (48000.0, 1, WaveformPeaks::kMaxFrames + 1),
                  "peaks: no channel, no frame, or past 2^53 frames is refused");
        test::ok (! w.prepare (48000.0, 1, 10, 0) && ! w.prepare (48000.0, 1, 10, WaveformPeaks::kMaxBuckets + 1)
               && ! w.prepare (48000.0, 1, 10, 10, (PeakMix) 4) && ! w.prepare (48000.0, 1, 10, 10, (PeakMix) -1),
                  "peaks: 0 buckets, too many buckets, or a mix code that names nothing is refused");
        test::ok (! w.prepared(), "peaks: a refused prepare leaves the object unprepared");

        const Planes p = programme (9, 2, 5000);
        const PeaksOut ref = runPeaks (p, 48000.0, 100, PeakMix::Average, { 5000 });
        test::run (w.prepare (48000.0, 2, 5000, 100, PeakMix::Average));
        auto v = ptrs (p);
        test::ok (! w.process (v.data(), 1, 100) && ! w.process (v.data(), 3, 100), "peaks: a narrower or wider call is refused");
        test::ok (! w.process (nullptr, 2, 100) && ! w.process (v.data(), 2, -1), "peaks: null planes or a negative length is refused");
        test::ok (w.framesSeen() == 0, "peaks: ... and nothing moved");
        test::ok (w.process (v.data(), 2, 0) && w.framesSeen() == 0, "peaks: n = 0 is a legal no-op");
        test::run (w.process (v.data(), 2, 4000));
        v = ptrs (p, 4000);
        test::ok (! w.process (v.data(), 2, 1001) && w.framesSeen() == 4000, "peaks: a call past the prepared length is refused whole");
        test::run (w.process (v.data(), 2, 1000));
        test::ok (w.complete() && std::memcmp (w.peaks().data(), ref.peaks.data(), 100 * sizeof (double)) == 0,
                  "peaks: after the refusals the file still reduces to the one-call bits");
        w.reset();
        v = ptrs (p);
        test::run (w.process (v.data(), 2, 5000));
        test::ok (std::memcmp (w.peaks().data(), ref.peaks.data(), 100 * sizeof (double)) == 0, "peaks: reset() and a replay give the same bits");

        StereoColumns s;
        test::ok (! s.prepare (0, 10) && ! s.prepare (2, 0) && ! s.prepare (2, 10, 0) && ! s.prepare (2, 10, StereoColumns::kMaxColumns + 1),
                  "stereo: no channel, no frame, 0 or too many columns is refused");
        const StereoOut sref = runStereo (p, 64, { 5000 });
        test::run (s.prepare (2, 5000, 64));
        v = ptrs (p);
        test::ok (! s.process (v.data(), 1, 10) && ! s.process (v.data(), 2, 5001) && s.framesSeen() == 0,
                  "stereo: a narrower call or one past the file is refused, nothing moved");
        test::run (s.process (v.data(), 2, 5000));
        test::ok (s.complete() && sameStereo (StereoOut { { s.width().begin(), s.width().end() }, { s.correlation().begin(), s.correlation().end() },
                                                          { s.rms().begin(), s.rms().end() }, s.maxRms(), s.columns(), s.isMono(), true }, sref),
                  "stereo: the file still reduces to the one-call bits");

        StereoColumns::Needle nd;
        test::ok (! StereoColumns::needle (p[0].data(), p[1].data(), 5000, 10, 9, nd) && ! StereoColumns::needle (p[0].data(), p[1].data(), 5000, 0, 5001, nd)
               && ! StereoColumns::needle (nullptr, p[1].data(), 5000, 0, 10, nd), "needle: from > to, past the planes, or null is refused");
        test::ok (StereoColumns::needle (p[0].data(), p[1].data(), 5000, 7, 7, nd) && nd.correlation == 1.0 && nd.width == 0.0 && nd.rms == 0.0,
                  "needle: an empty stretch reads corr +1, width 0, loud 0 (n floored at 1)");
    }

    //==========================================================================================================
    test::group ("process() does not allocate");
    {
        const Planes p = programme (31, 2, 70000);
        WaveformPeaks w; StereoColumns s;
        test::run (w.prepare (48000.0, 2, 70000, 1100, PeakMix::Max));
        test::run (s.prepare (2, 70000, 1100));
        const auto v = ptrs (p);
        const long before = g_allocs.load();
        const bool a = w.process (v.data(), 2, 70000);
        const bool b = s.process (v.data(), 2, 70000);
        StereoColumns::Needle nd;
        const bool c = StereoColumns::needle (p[0].data(), p[1].data(), 70000, 0, 70000, nd);
        test::okNoAlloc (g_allocs.load() == before && a && b && c, "WaveformPeaks::process, StereoColumns::process and needle() allocate nothing");
    }

    return test::report();
}
