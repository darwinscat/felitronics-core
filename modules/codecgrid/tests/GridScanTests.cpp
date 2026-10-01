// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid::GridScan on programmes that carry a known grid (SyntheticCodec.h): the scan must find the offset the
// programme was coded at, with the transform it was coded with — and must find nothing in the same programme
// uncoded, nor with another codec's transform. Plus the form: steps of any size give the same curve, a stretch
// leaves nothing behind for the next, and no call after prepare() allocates.

#include <felitronics/codecgrid/BackResampler.h>
#include <felitronics/codecgrid/GridScan.h>

#include "SyntheticCodec.h"
#include "alloc_counter.h"
#include "felitronics_test.h"

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics::codecgrid;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
const char* nameOf (Transform t) { return t == Transform::Mp3 ? "MP3" : t == Transform::AacSine ? "AAC sine" : t == Transform::AacKbd ? "AAC KBD" : "CELT"; }

struct Scanned
{
    std::vector<float> curve;
    GridReading reading;
    bool ran = false;
};

Scanned scan (GridScan<>& scanner, Transform t, const synthetic::Stereo& x, int stepSize = 1 << 20)
{
    Scanned s;
    s.curve.assign ((std::size_t) hopOf (t) * kCells, 0.0f);
    if (! scanner.begin (t, x.left.data(), x.right.data(), (int) x.left.size(), s.curve.data())) return s;
    while (! scanner.step (stepSize)) {}
    std::vector<double> d (ReadScratch::doublesFor (hopOf (t)));
    std::vector<float> f (ReadScratch::floatsFor (hopOf (t)));
    ReadScratch scratch;
    scratch.local = d.data();
    scratch.sort = d.data() + hopOf (t);
    scratch.dip = f.data();
    s.reading = readCurve (s.curve.data(), hopOf (t), scratch);
    s.ran = true;
    return s;
}
}

int main()
{
    std::printf ("felitronics::codecgrid scan tests\n");
    constexpr int n = 40000;                       // 0.9 s at 44.1 kHz: 37 AAC frames, 66 granules
    constexpr double zeroShare = 0.3;
    const synthetic::Stereo original = synthetic::programme (n, 1u);
    const GridRule rule;

    GridScan<> scanner;
    ok (scanner.prepare (n), "prepare() for stretches of 40000 samples");

    group ("the fixtures are what they say: with nothing zeroed each round trip returns its input");
    {
        for (Transform t : { Transform::Mp3, Transform::AacSine, Transform::AacKbd, Transform::Celt })
        {
            const auto back = synthetic::coded (original, t, 7, 0.0);
            const double left = synthetic::residualDb (original.left, back.left, 4096), right = synthetic::residualDb (original.right, back.right, 4096);
            std::printf ("    %-9s round trip residual: %.1f dB left, %.1f dB right\n", nameOf (t), left, right);
            // the MDCT fixtures reconstruct to single precision; the Layer III filterbank is near-perfect by
            // design, not perfect — its own ripple and aliasing sit about 90 dB down
            ok (back.left.size() == original.left.size() && std::fmax (left, right) < (t == Transform::Mp3 ? -70.0 : -100.0),
                std::string (nameOf (t)) + ": the fixture reconstructs its input, in place");
        }
    }

    group ("a programme coded at an offset is found at that offset, by its own transform");
    {
        struct Case { Transform t; int offset; };
        for (const Case c : { Case { Transform::Mp3, 0 }, Case { Transform::Mp3, 133 }, Case { Transform::Mp3, 575 },
                              Case { Transform::AacSine, 0 }, Case { Transform::AacSine, 480 }, Case { Transform::AacKbd, 1023 },
                              Case { Transform::Celt, 0 }, Case { Transform::Celt, 548 } })
        {
            const auto x = synthetic::coded (original, c.t, c.offset, zeroShare);
            const Scanned s = scan (scanner, c.t, x);
            const std::string tag = std::string (nameOf (c.t)) + " coded at " + std::to_string (c.offset) + ": ";
            ok (s.ran, tag + "the scan ran");
            std::printf ("    %-9s at %4d: found %4d, score %7.1f, runner-up %5.1f, dip %6.1f dB, cells %2d\n", nameOf (c.t), c.offset,
                         s.reading.offset, s.reading.score, s.reading.second, s.reading.localDipSum, s.reading.cells);
            ok (s.reading.offset == c.offset, tag + "the offset is the one it was coded at");
            const auto repeat = scanner.frameRepeat (s.reading.offset);
            ok (repeat.frames > 0 && repeat.period == 0, tag + "its frames do not repeat: the score stands");
            ok (rule.found (s.reading) && s.reading.cells >= 24, tag + "the rule finds it, in at least 24 of the 32 cells");

            float share[kZeroBands];
            std::vector<std::uint8_t> map ((std::size_t) 200 * kZeroBands, 0);
            const int rows = scanner.zeroProfile (s.reading.offset, share, map.data(), 200);
            double mean = 0.0;
            for (float v : share) mean += v;
            mean /= kZeroBands;
            ok (rows > 20 && rows <= 200, tag + "the zero map has a row per frame");
            // 30 % of every frame's coefficients were zeroed per channel; the share published is the largest of
            // the four signals, so it is at least that of one channel and cannot be far above it
            ok (mean > 0.2 && mean < 0.5, tag + "the zero share at the found offset is about the share that was zeroed (" + std::to_string (mean) + ")");
            float off[kZeroBands];
            (void) scanner.zeroProfile ((s.reading.offset + hopOf (c.t) / 2) % hopOf (c.t), off);
            double away = 0.0;
            for (float v : off) away += v;
            ok (away / kZeroBands < 0.25 * mean, tag + "half a frame away the zeros are gone");
        }
    }

    group ("what must not be found: the programme uncoded, and a grid read with another codec's transform");
    {
        for (Transform t : { Transform::Mp3, Transform::AacSine, Transform::AacKbd, Transform::Celt })
        {
            const Scanned s = scan (scanner, t, original);
            std::printf ("    uncoded, %-9s score %5.2f, dip %5.2f dB, cells %d\n", nameOf (t), s.reading.score, s.reading.localDipSum, s.reading.cells);
            ok (s.ran && ! rule.found (s.reading) && ! rule.several (s.reading), std::string ("the uncoded programme carries no ") + nameOf (t) + " grid");
        }
        const auto aac = synthetic::coded (original, Transform::AacSine, 480, zeroShare);
        const auto mp3 = synthetic::coded (original, Transform::Mp3, 133, zeroShare);
        const auto celt = synthetic::coded (original, Transform::Celt, 548, zeroShare);
        struct Cross { const synthetic::Stereo* x; Transform t; const char* what; };
        for (const Cross c : { Cross { &aac, Transform::Mp3, "an AAC grid read as MP3" }, Cross { &aac, Transform::Celt, "an AAC grid read as CELT" },
                               Cross { &mp3, Transform::AacSine, "an MP3 grid read as AAC" }, Cross { &mp3, Transform::Celt, "an MP3 grid read as CELT" },
                               Cross { &celt, Transform::AacSine, "a CELT grid read as AAC" }, Cross { &celt, Transform::Mp3, "a CELT grid read as MP3" } })
        {
            const Scanned s = scan (scanner, c.t, *c.x);
            std::printf ("    %-26s score %5.2f, dip %5.2f dB, cells %d\n", c.what, s.reading.score, s.reading.localDipSum, s.reading.cells);
            ok (s.ran && ! rule.found (s.reading), std::string (c.what) + " is not found");
        }
    }

    group ("the curve does not depend on how the scan was stepped, nor on what was scanned before");
    {
        const auto aac = synthetic::coded (original, Transform::AacKbd, 300, zeroShare);
        const auto mp3 = synthetic::coded (original, Transform::Mp3, 133, zeroShare);
        for (Transform t : { Transform::Mp3, Transform::AacKbd, Transform::Celt })
        {
            const synthetic::Stereo& x = t == Transform::Mp3 ? mp3 : aac;
            const Scanned whole = scan (scanner, t, x);
            const Scanned single = scan (scanner, t, x, 1);
            const Scanned odd = scan (scanner, t, x, 37);
            (void) scan (scanner, t == Transform::Mp3 ? Transform::Celt : Transform::Mp3, original);        // something else in between
            const Scanned again = scan (scanner, t, x);
            const std::size_t bytes = whole.curve.size() * sizeof (float);
            ok (whole.ran && std::memcmp (whole.curve.data(), single.curve.data(), bytes) == 0 && std::memcmp (whole.curve.data(), odd.curve.data(), bytes) == 0,
                std::string (nameOf (t)) + ": one offset per step, 37 per step and all at once give the same bits");
            ok (std::memcmp (whole.curve.data(), again.curve.data(), bytes) == 0, std::string (nameOf (t)) + ": the same stretch after another scan gives the same bits");
        }
    }

    group ("what begin() refuses, and that a refusal begins nothing");
    {
        std::vector<float> curve ((std::size_t) 1024 * kCells, 123.0f);
        GridScan<> fresh;
        ok (! fresh.begin (Transform::AacSine, original.left.data(), original.right.data(), n, curve.data()), "an unprepared scanner refuses");
        ok (fresh.step (10) && fresh.offsets() == 0, "...and step() on it is a completed nothing");
        ok (! scanner.begin (Transform::AacSine, nullptr, original.right.data(), n, curve.data()), "a null left plane");
        ok (! scanner.begin (Transform::AacSine, original.left.data(), nullptr, n, curve.data()), "a null right plane");
        ok (! scanner.begin (Transform::AacSine, original.left.data(), original.right.data(), n, nullptr), "a null curve");
        ok (! scanner.begin (Transform::AacSine, original.left.data(), original.right.data(), n + 1, curve.data()), "a stretch longer than prepare() was told");
        for (Transform t : { Transform::Mp3, Transform::AacSine, Transform::AacKbd, Transform::Celt })
        {
            const int least = minSamplesOf (t, 2);
            ok (! scanner.begin (t, original.left.data(), original.right.data(), least - 1, curve.data()), std::string (nameOf (t)) + ": one sample short of two frames at every offset is refused");
            ok (scanner.begin (t, original.left.data(), original.right.data(), least, curve.data()), std::string (nameOf (t)) + ": exactly two frames at every offset is accepted");
            while (! scanner.step (64)) {}
            bool finite = true;
            for (int i = 0; i < hopOf (t) * kCells; ++i) finite = finite && std::isfinite (curve[(std::size_t) i]);
            ok (finite, std::string (nameOf (t)) + ": ...and its curve is finite at every offset");
        }
        ok (! scanner.begin (Transform::AacSine, original.left.data(), original.right.data(), 0, curve.data()) && scanner.offsets() == 0, "an empty stretch begins nothing");
        float share[kZeroBands];
        ok (scanner.zeroProfile (0, share) == 0, "a zero profile with no completed scan is empty");
        GridScan<> none;
        ok (! none.prepare (0) && ! none.prepare (-5) && ! none.prepared(), "prepare() refuses a capacity of nothing");
    }

    group ("silence, a constant and non-finite input do not break the scan");
    {
        synthetic::Stereo z { std::vector<float> ((std::size_t) n, 0.0f), std::vector<float> ((std::size_t) n, 0.0f) };
        for (Transform t : { Transform::Mp3, Transform::AacSine, Transform::Celt })
        {
            const Scanned s = scan (scanner, t, z);
            bool finite = true;
            for (float v : s.curve) finite = finite && std::isfinite (v);
            ok (s.ran && finite && ! rule.found (s.reading), std::string (nameOf (t)) + ": digital silence gives a finite curve and no grid");
        }
        synthetic::Stereo bad = original;
        bad.left[20000] = std::nanf ("");
        bad.right[20001] = INFINITY;
        const Scanned s = scan (scanner, Transform::AacSine, bad);
        bool finite = true;
        for (float v : s.curve) finite = finite && std::isfinite (v);
        ok (s.ran && finite, "a NaN and an infinity in the stretch leave every cell finite");
        ok (std::isfinite (s.reading.score) && ! rule.found (s.reading), "...and are not read as a grid");
    }

    group ("no call after prepare() allocates");
    {
        const auto aac = synthetic::coded (original, Transform::AacSine, 480, zeroShare);
        std::vector<float> curve ((std::size_t) 1024 * kCells);
        float share[kZeroBands];
        std::vector<std::uint8_t> map ((std::size_t) 64 * kZeroBands);
        bool accepted = true;
        const auto before = alloc::count.load();
        for (Transform t : { Transform::Mp3, Transform::AacSine, Transform::AacKbd, Transform::Celt })
        {
            accepted = accepted && scanner.begin (t, aac.left.data(), aac.right.data(), n, curve.data());
            while (! scanner.step (100)) {}
            (void) scanner.zeroProfile (5, share, map.data(), 64);
        }
        const auto after = alloc::count.load();
        ok (accepted, "the four scans were accepted");
        felitronics::test::okNoAlloc (after == before, "begin(), step() and zeroProfile() did not allocate");
    }

    group ("a rhythm on the frame's period makes a dip the rule accepts, and the frames say it is not a grid");
    {
        // One noise burst, the same samples every time, every 5292 samples of 44.1 kHz: a sixteenth note at 125 bpm,
        // 120 ms, six frames of CELT. Nothing in it was ever coded. Two seconds, taken back to 48 kHz as the
        // detector does for that hypothesis.
        const int n44 = 90000, n48 = 96000, period = 5292, length = 2646;
        BackResampler rs;
        GridScan<> wide;
        int fooled = 0;
        if (felitronics::test::run (rs.prepare (160, 147) && wide.prepare (n48)))
            for (std::uint32_t seed : { 1u, 4u, 22u, 38u, 54u })      // 22, 38 and 54 are the three of the first sixty that fool the rule
            {
                std::uint32_t state = seed * 2654435761u;
                const auto noise = [&state]() noexcept
                {
                    float v = 0.0f;
                    for (int i = 0; i < 4; ++i)
                    {
                        state = state * 1664525u + 1013904223u;
                        v += (float) (state >> 8) * (1.0f / 16777216.0f);
                    }
                    return (v - 2.0f) * 1.7320508f;
                };
                std::vector<float> burst ((std::size_t) length);
                float previous = noise();
                for (int i = 0; i < length; ++i)
                {
                    const float v = noise();
                    burst[(std::size_t) i] = (v - previous) * (float) std::exp (-(double) i / 352.8);
                    previous = v;
                }
                std::vector<float> l44 ((std::size_t) n44, 0.0f), r44 ((std::size_t) n44, 0.0f);
                for (int at = 0; at + length < n44; at += period)
                    for (int i = 0; i < length; ++i)
                    {
                        l44[(std::size_t) (at + i)] += 0.072f * burst[(std::size_t) i];
                        r44[(std::size_t) (at + i)] += 0.072f * burst[(std::size_t) i];
                    }
                for (int i = 0; i < n44; ++i)
                {
                    l44[(std::size_t) i] += 0.0006f * noise();
                    r44[(std::size_t) i] += 0.0006f * noise();
                }
                synthetic::Stereo x;
                x.left.assign ((std::size_t) n48, 0.0f);
                x.right.assign ((std::size_t) n48, 0.0f);
                rs.resample (l44.data(), n44, 200, n48, x.left.data());
                rs.resample (r44.data(), n44, 200, n48, x.right.data());
                const Scanned s = scan (wide, Transform::Celt, x);
                const auto repeat = wide.frameRepeat (s.reading.offset);
                const double factor = repeatFactor ((double) repeat.repeat, repeat.period, repeat.frames);
                GridReading worth = s.reading;
                worth.score *= factor;
                worth.second *= factor;
                const bool found = rule.found (s.reading);
                if (found) ++fooled;
                std::printf ("    burst %u: score %5.1f, runner-up %4.1f, dip %5.1f dB, cells %2d%s; %d frames repeat every %d by %.3f: the score is worth %.1f\n",
                             (unsigned) seed, s.reading.score, s.reading.second, s.reading.localDipSum, s.reading.cells, found ? ", FOUND by the rule" : "",
                             repeat.frames, repeat.period, (double) repeat.repeat, worth.score);
                const std::string tag = "burst " + std::to_string (seed) + ": ";
                ok (s.ran && repeat.frames > 90 && repeat.period == 6 && repeat.repeat > 0.75f, tag + "the frames repeat every six");
                ok (! rule.found (worth) && worth.score < 0.75 * rule.minScore, tag + "and for what it is worth among six different frames the score is far under the rule's");
            }
        ok (fooled >= 1, "at face value the rule finds a grid in at least one of them");
        ok (wide.frameRepeat (-1).frames == 0 && wide.frameRepeat (hopOf (Transform::Celt)).frames == 0, "an offset that is not there has no reading");
        GridScan<> fresh;
        ok (fresh.prepare (4096) && fresh.frameRepeat (0).frames == 0, "nor has a scan that was never run");
    }

    group ("the back resampler: zero phase, unity gain, and the published length");
    {
        ok (BackResampler::ratioFor (44100, 48000).ok && BackResampler::ratioFor (44100, 48000).up == 160 && BackResampler::ratioFor (44100, 48000).down == 147, "44.1 -> 48 kHz is 160 : 147");
        ok (BackResampler::ratioFor (48000, 44100).up == 147 && BackResampler::ratioFor (48000, 44100).down == 160, "48 -> 44.1 kHz is 147 : 160");
        ok (! BackResampler::ratioFor (44100, 44100).ok && ! BackResampler::ratioFor (96000, 48000).ok, "any other pair is not offered");
        BackResampler bad;
        ok (! bad.prepare (0, 147) && ! bad.prepare (160, 0) && ! bad.prepare (5000, 1) && ! bad.prepared(), "prepare() refuses a ratio of nothing or of thousands");
        for (const auto ratio : { BackResampler::ratioFor (44100, 48000), BackResampler::ratioFor (48000, 44100) })
        {
            BackResampler rs;
            if (! felitronics::test::run (rs.prepare (ratio.up, ratio.down))) continue;
            const std::string tag = std::to_string (ratio.up) + ":" + std::to_string (ratio.down) + " ";
            const int in = 147 * 160 * 2;                                // 47040: both rates divide it
            ok (rs.outputLength (in) == (std::int64_t) in * ratio.up / ratio.down && rs.outputLength (0) == 0
                    && rs.outputLength (1) == (ratio.up > ratio.down ? 2 : 1) && rs.outputLength (ratio.down) == ratio.up && rs.outputLength (ratio.down + 1) == ratio.up + (ratio.up > ratio.down ? 2 : 1),
                tag + "the output length is ceil (n up / down)");
            // a constant stays the constant (unity gain at DC), away from the ends
            std::vector<float> one ((std::size_t) in, 1.0f), out (4000);
            rs.resample (one.data(), in, 10000, 4000, out.data());
            double worst = 0.0;
            for (float v : out) worst = std::fmax (worst, std::fabs ((double) v - 1.0));
            ok (worst < 2.0e-5, tag + "a constant comes out as the same constant (worst " + std::to_string (worst) + ")");
            // zero phase: a slow sine sampled at the output's own instants j down / up
            const double f = 0.01;                                        // cycles per input sample, far inside the passband
            std::vector<float> sine ((std::size_t) in);
            for (int i = 0; i < in; ++i) sine[(std::size_t) i] = (float) std::sin (2.0 * synthetic::kPi * f * i);
            rs.resample (sine.data(), in, 10000, 4000, out.data());
            worst = 0.0;
            for (int c = 0; c < 4000; ++c)
                worst = std::fmax (worst, std::fabs ((double) out[(std::size_t) c] - std::sin (2.0 * synthetic::kPi * f * (double) (10000 + c) * ratio.down / ratio.up)));
            ok (worst < 1.0e-3, tag + "output sample j sits at input time j down / up (worst " + std::to_string (worst) + ")");
            // an impulse at input sample k lands on output sample k up / down when that is whole: no latency
            std::vector<float> impulse ((std::size_t) in, 0.0f);
            const int k = ratio.down * 100;
            impulse[(std::size_t) k] = 1.0f;
            rs.resample (impulse.data(), in, (std::int64_t) k * ratio.up / ratio.down - 50, 101, out.data());
            int at = 0;
            for (int c = 1; c < 101; ++c) if (std::fabs (out[(std::size_t) c]) > std::fabs (out[(std::size_t) at])) at = c;
            bool symmetric = true;
            for (int c = 1; c <= 50; ++c) symmetric = symmetric && std::fabs ((double) out[(std::size_t) (50 + c)] - (double) out[(std::size_t) (50 - c)]) < 1.0e-6;
            ok (at == 50 && symmetric, tag + "an impulse comes out centred on its own instant, symmetric about it");
            // the ends: input taken as zero outside, no read past the array
            rs.resample (one.data(), in, -20, 40, out.data());
            rs.resample (one.data(), in, rs.outputLength (in) - 20, 40, out.data());
            ok (std::isfinite (out[0]) && std::isfinite (out[39]), tag + "outputs before the first and after the last input sample are computed from zeros");
        }
    }

    return felitronics::test::report();
}
