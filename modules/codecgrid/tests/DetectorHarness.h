// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// What the detector suites share: feeding a programme in a given split, and making the programmes — coded at the
// programme's rate, coded at the other rate and converted, coded in part.

#include <felitronics/codecgrid/CodecGridDetector.h>

#include "SyntheticCodec.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace harness
{
using namespace felitronics::codecgrid;
using Detector = CodecGridDetector<>;

inline const char* name (Verdict v)
{
    return v == Verdict::Confirmed ? "Confirmed" : v == Verdict::InPlaces ? "InPlaces" : v == Verdict::SeveralGrids ? "SeveralGrids" : v == Verdict::NoGrid ? "NoGrid" : "NotExamined";
}
inline const char* name (Ground g)
{
    return g == Ground::Both ? "Both" : g == Ground::TwoWindows ? "TwoWindows" : g == Ground::PhaseAgreement ? "PhaseAgreement" : "None";
}

// The stretches of the suites are a quarter of a second: every property under test is one of the detector's
// logic, and an eighth of the arithmetic proves it as well. The rule's numbers were measured at 2 s; the fixtures
// zero 30 % of every frame and stand far above them at any length.
inline constexpr double kStretchSeconds = 0.25;

inline CodecGridParams params (Depth depth = Depth::Follow, double windowSeconds = kStretchSeconds)
{
    CodecGridParams p;
    p.depth = depth;
    p.windowSeconds = windowSeconds;
    return p;
}

// Feed x in calls of `block` frames (0: one call). Returns false on a refused call.
inline bool feed (Detector& det, const synthetic::Stereo& x, int channels, int block = 0, std::size_t frames = SIZE_MAX)
{
    const std::size_t n = std::min (frames, x.left.size());
    std::size_t at = 0;
    while (at < n)
    {
        const std::size_t take = block > 0 ? std::min ((std::size_t) block, n - at) : n - at;
        const float* in[2] { x.left.data() + at, x.right.data() + at };
        if (! det.process (in, channels, (int) take)) return false;
        at += take;
    }
    return true;
}

inline std::unique_ptr<Detector> analyse (const synthetic::Stereo& x, double rate, const CodecGridParams& p, int channels = 2, int block = 0)
{
    auto det = std::make_unique<Detector>();
    det->setParams (p);
    if (! det->prepare (rate, channels, x.left.size())) return nullptr;
    if (! feed (*det, x, channels, block)) return nullptr;
    if (! det->finish()) return nullptr;
    return det;
}

inline int scans (const Detector& det)
{
    int n = 0;
    for (int w = 0; w < det.windows(); ++w)
        for (int h = 0; h < kHypotheses; ++h) n += det.window (w).hypotheses[h].scanned ? 1 : 0;
    return n;
}

// A programme coded at `codecRate` and converted, zero phase, to `fileRate`.
inline synthetic::Stereo codedAtOtherRate (int frames, std::uint32_t seed, Transform t, int offset, int codecRate, int fileRate)
{
    const auto original = synthetic::programme (frames, seed);
    const auto coded = synthetic::coded (original, t, offset, 0.3);
    const auto ratio = BackResampler::ratioFor (codecRate, fileRate);
    BackResampler rs;
    synthetic::Stereo out;
    if (! ratio.ok || ! rs.prepare (ratio.up, ratio.down)) return out;
    const int count = (int) rs.outputLength (frames);
    out.left.resize ((std::size_t) count);
    out.right.resize ((std::size_t) count);
    rs.resample (coded.left.data(), frames, 0, count, out.left.data());
    rs.resample (coded.right.data(), frames, 0, count, out.right.data());
    return out;
}

// `x` with the frames [from, to) taken from `y`.
inline synthetic::Stereo spliced (synthetic::Stereo x, const synthetic::Stereo& y, std::size_t from, std::size_t to)
{
    for (std::size_t i = from; i < to && i < x.left.size(); ++i)
    {
        x.left[i] = y.left[i];
        x.right[i] = y.right[i];
    }
    return x;
}

// Everything a report says, compared field by field (the struct has padding; its bytes are not its value).
inline bool sameReport (const Detector& a, const Detector& b)
{
    const auto& p = a.result();
    const auto& q = b.result();
    bool same = p.verdict == q.verdict && p.reason == q.reason && p.ground == q.ground && p.family == q.family && p.transform == q.transform
             && p.codecRate == q.codecRate && p.hop == q.hop && p.gridPhase == q.gridPhase && p.windows == q.windows
             && p.windowsExamined == q.windowsExamined && p.windowsFound == q.windowsFound && p.windowsAgreeing == q.windowsAgreeing
             && p.bestWindow == q.bestWindow && p.bestHypothesis == q.bestHypothesis && p.bestOffset == q.bestOffset
             && std::memcmp (&p.bestScore, &q.bestScore, sizeof (double)) == 0 && std::memcmp (p.zeroByBand, q.zeroByBand, sizeof (p.zeroByBand)) == 0
             && std::memcmp (&p.zeroShare, &q.zeroShare, sizeof (float)) == 0 && p.nonFiniteSamples == q.nonFiniteSamples
             && a.windows() == b.windows() && a.curveOffsets() == b.curveOffsets() && a.zeroMapRows() == b.zeroMapRows();
    if (! same) return false;
    for (int w = 0; w < a.windows(); ++w)
    {
        const auto& u = a.window (w);
        const auto& v = b.window (w);
        same = same && u.startFrame == v.startFrame && u.complete == v.complete && u.silent == v.silent && u.best == v.best
            && std::memcmp (&u.rmsDb, &v.rmsDb, sizeof (double)) == 0;
        for (int h = 0; h < kHypotheses; ++h)
        {
            const auto& x = u.hypotheses[h];
            const auto& y = v.hypotheses[h];
            same = same && x.scanned == y.scanned && x.found == y.found && x.several == y.several && x.gridPhase == y.gridPhase
                && x.reading.offset == y.reading.offset && x.reading.cells == y.reading.cells
                && std::memcmp (&x.reading.score, &y.reading.score, sizeof (double)) == 0 && std::memcmp (&x.reading.second, &y.reading.second, sizeof (double)) == 0
                && std::memcmp (&x.reading.localDipSum, &y.reading.localDipSum, sizeof (double)) == 0;
        }
    }
    if (a.curveOffsets() > 0) same = same && std::memcmp (a.curve(), b.curve(), (std::size_t) a.curveOffsets() * kCells * sizeof (float)) == 0;
    if (a.zeroMapRows() > 0) same = same && std::memcmp (a.zeroMap(), b.zeroMap(), (std::size_t) a.zeroMapRows() * kZeroBands) == 0;
    return same;
}
}
