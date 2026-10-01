// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid::CodecGridDetector's FORM — what the call contract of this library (DSP-ARCHITECTURE.md, law 11)
// asks of an offline measurer: prepare() disarms, validates and only then writes; what it asks the heap for is
// published first and is exactly what it takes; process() refuses a request it cannot honour in full and
// consumes nothing of it; the report does not depend on how the programme was split into calls; nothing after
// prepare() allocates; finish is idempotent; a refused prepare() leaves nothing readable.

#include "DetectorHarness.h"
#include "alloc_counter.h"
#include "felitronics_test.h"

#include <cmath>
#include <limits>
#include <string>

using namespace harness;
using felitronics::test::group;
using felitronics::test::ok;

int main()
{
    std::printf ("felitronics::codecgrid detector law tests\n");
    // Three seconds and stretches of a quarter of a second: the form does not depend on the length of a stretch,
    // and the fixture's grid stands far above the rule at any length.
    constexpr int n = 44100 * 3;
    constexpr double kStretch = kStretchSeconds;
    const auto lawParams = [] (Depth depth) { return params (depth, kStretch); };
    const auto original = synthetic::programme (n, 7u);
    const auto coded = synthetic::coded (original, Transform::AacSine, 480, 0.3);
    const float* planes[2] { coded.left.data(), coded.right.data() };

    group ("prepare(): what is refused, and that a refusal leaves nothing readable");
    {
        Detector d;
        d.setParams (lawParams (Depth::Follow));
        const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
        ok (! d.prepare (nan, 2, 1000) && ! d.prepare (inf, 2, 1000) && ! d.prepare (-44100.0, 2, 1000) && ! d.prepare (0.0, 2, 1000), "a rate that is NaN, infinite, negative or zero");
        ok (! d.prepare (7999.0, 2, 1000) && ! d.prepare (768001.0, 2, 1000), "a rate under the library's floor or over the ceiling");
        ok (d.prepare (8000.0, 2, 1000) && d.prepare (768000.0, 2, 1000), "...and the floor and the ceiling themselves are accepted");
        ok (! d.prepare (44100.0, 0, 1000) && ! d.prepare (44100.0, -1, 1000) && ! d.prepare (44100.0, felitronics::core::kMaxChannels + 1, 1000), "no channel, a negative count, one more than the library's width");
        ok (d.prepare (44100.0, felitronics::core::kMaxChannels, 1000) && d.prepare (44100.0, 1, 1000), "the widest and the narrowest widths are accepted");
        ok (! d.prepare (44100.0, 2, Detector::kMaxFrames + 1) && d.prepare (44100.0, 2, 0), "a length past the ceiling is refused; an empty programme is accepted");
        for (double seconds : { nan, inf, 0.0, -1.0, 0.2499, 8.0001 })
        {
            auto p = params();
            p.windowSeconds = seconds;
            d.setParams (p);
            ok (! d.prepare (44100.0, 2, (std::uint64_t) n), "a stretch length of " + std::to_string (seconds) + " s is refused");
            ok (! Detector::storageFor (44100.0, 2, (std::uint64_t) n, p).ok && Detector::storageFor (44100.0, 2, (std::uint64_t) n, p).bytes() == 0, "...and its storage says so and asks for nothing");
        }
        auto bad = params();
        bad.phaseTolerance = -1;
        d.setParams (bad);
        ok (! d.prepare (44100.0, 2, (std::uint64_t) n), "a negative phase tolerance is refused");
        bad = params();
        bad.phaseAgreeWindows = 1;
        d.setParams (bad);
        ok (! d.prepare (44100.0, 2, (std::uint64_t) n), "a phase ground of one stretch is refused: one stretch agrees with itself");

        // a finished analysis, then a refused prepare(): the old report must be gone
        d.setParams (lawParams (Depth::Verdict));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n) && d.process (planes, 2, n) && d.finish(), "a finished analysis");
        ok (d.result().verdict == Verdict::Confirmed && d.finished() && d.curveOffsets() > 0, "...with a report");
        ok (! d.prepare (nan, 2, (std::uint64_t) n), "then a refused prepare()");
        ok (d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::NotFinished && d.result().family == Family::None
                && d.result().bestWindow < 0 && d.windows() == 0 && d.curveOffsets() == 0 && d.zeroMapRows() == 0 && ! d.finished(),
            "leaves no verdict, no stretch, no curve: nothing of the old report is readable");
        ok (! d.process (planes, 2, 10) && ! d.finish() && ! d.finishStep() && d.progress() <= 0.0, "...and the object is unprepared: process() and finish() refuse");
        ok (d.window (-1).best < 0 && d.window (99).best < 0 && d.hypothesis (-1).codecRate == 0 && d.hypothesis (99).codecRate == 0, "reading a stretch or a hypothesis that is not there returns an empty one");
    }

    group ("process(): a request that cannot be honoured in full is refused whole and consumes nothing");
    {
        Detector d;
        d.setParams (lawParams (Depth::Verdict));
        ok (! d.process (planes, 2, 100), "before prepare()");
        ok (d.prepare (44100.0, 2, (std::uint64_t) n), "prepare()");
        ok (! d.process (planes, 2, -1) && ! d.process (planes, -2, 100), "a negative length or width");
        ok (! d.process (nullptr, 2, 100), "no planes at all");
        ok (! d.process (planes, 1, 100) && ! d.process (planes, 3, 100), "a width that is not the prepared one, narrower or wider");
        const float* holed[2] { coded.left.data(), nullptr };
        ok (! d.process (holed, 2, 100), "a missing plane");
        ok (d.process (planes, 2, 0) && d.process (nullptr, 0, 0) == false, "an empty call is accepted at the prepared width, and refused at another");
        ok (! d.process (planes, 2, n + 1), "one frame more than the prepared length, in one call");
        // every refusal above consumed nothing: the whole programme still fits, and gives the report of a clean run
        ok (d.process (planes, 2, n), "after all those refusals the whole programme is still accepted");
        ok (! d.process (planes, 2, 1), "and then not one frame more");
        ok (d.finish(), "finish()");
        const auto clean = analyse (coded, 44100.0, lawParams (Depth::Verdict));
        ok (clean != nullptr && sameReport (d, *clean), "the report is that of a run with no refused call in it");
        ok (! d.process (planes, 2, 0) && ! d.process (planes, 2, 10), "after finish() every process() is refused, the empty one too");
        ok (d.finish() && d.finishStep() && sameReport (d, *clean), "finish() and finishStep() again change nothing");
        d.reset();
        ok (d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::NotFinished && ! d.finished() && d.windows() == 8, "reset(): no report, the stretches still placed");
        ok (d.process (planes, 2, n) && d.finish() && sameReport (d, *clean), "the same programme again after reset() gives the same report");
        // finishing begins with the first finishStep(): process() is refused from then on
        d.reset();
        ok (d.process (planes, 2, n / 2), "half the programme");
        ok (! d.finishStep(), "one finish step: work remains");
        ok (! d.process (planes, 2, 10), "...and process() is refused once finishing has begun");
    }

    group ("the report does not depend on how the programme was split into calls");
    {
        const auto whole = analyse (coded, 44100.0, lawParams (Depth::Follow));
        if (felitronics::test::run (whole != nullptr))
        {
            for (int block : { 1, 7, 64, 4096, 44099, 100000 })
            {
                const auto split = analyse (coded, 44100.0, lawParams (Depth::Follow), 2, block);
                ok (split != nullptr && sameReport (*whole, *split), "calls of " + std::to_string (block) + " frames: the same report, the same curve, the same zero map");
            }
            ok (whole->result().verdict == Verdict::Confirmed && whole->result().windowsFound == 8, "...and that report is a confirmed grid read in all eight stretches, so every stretch's copy was compared");
            // an uneven split whose boundaries fall inside stretches and margins
            Detector d;
            d.setParams (lawParams (Depth::Follow));
            bool fed = d.prepare (44100.0, 2, (std::uint64_t) n);
            std::size_t at = 0;
            std::uint32_t r = 12345u;
            while (fed && at < (std::size_t) n)
            {
                r = r * 1664525u + 1013904223u;
                const std::size_t take = std::min ((std::size_t) n - at, (std::size_t) (1u + (r >> 16) % 30000u));
                const float* in[2] { coded.left.data() + at, coded.right.data() + at };
                fed = d.process (in, 2, (int) take);
                at += take;
            }
            ok (fed && d.finish() && sameReport (*whole, d), "calls of random lengths: the same report");
        }
    }

    group ("the bounded steps: the analysis by finishStep() is the analysis by finish()");
    {
        const auto whole = analyse (coded, 44100.0, lawParams (Depth::Follow));
        Detector d;
        d.setParams (lawParams (Depth::Follow));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n) && d.process (planes, 2, n), "prepared and fed");
        ok (d.progress() <= 0.0 && ! d.finished(), "before the first step: no progress, not finished");
        int steps = 0;
        bool monotone = true, readable = true;
        double last = 0.0;
        while (! d.finishStep())
        {
            ++steps;
            const double p = d.progress();
            monotone = monotone && p >= 0.0 && p <= 1.0;
            last = p;
            // a report read in the middle is the unfinished one, never a half-written verdict
            readable = readable && d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::NotFinished;
        }
        std::printf ("    %d finish steps; progress just before the last one %.3f\n", steps, last);
        ok (steps > 100, "the analysis is cut into many steps");
        ok (monotone, "progress stays within 0 .. 1");
        ok (readable, "between steps the report says NotFinished and nothing else");
        ok (d.finished() && d.progress() >= 1.0, "then finished, progress 1");
        ok (whole != nullptr && sameReport (*whole, d), "and the report is the one finish() gives");
    }

    group ("a programme that ended early is examined in the stretches that arrived whole");
    {
        Detector d;
        d.setParams (lawParams (Depth::Exhaustive));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n), "prepared for the whole programme");
        ok (feed (d, coded, 2, 4096, (std::size_t) n / 2), "fed half of it");
        ok (d.finish(), "finish()");
        int complete = 0;
        bool order = true;
        for (int w = 0; w < d.windows(); ++w)
        {
            if (d.window (w).complete) ++complete;
            if (w > 0) order = order && ! (d.window (w).complete && ! d.window (w - 1).complete);
            bool scanned = false;
            for (int h = 0; h < kHypotheses; ++h) scanned = scanned || d.window (w).hypotheses[h].scanned;
            order = order && scanned == d.window (w).complete;
        }
        std::printf ("    %d of %d stretches arrived whole\n", complete, d.windows());
        ok (complete >= 2 && complete <= 4 && order, "the stretches of the first half are complete and scanned, the later ones neither");
        ok (d.result().verdict == Verdict::Confirmed && d.result().windowsExamined == complete && d.result().windowsFound == complete, "and the grid is confirmed from them");
        // nothing fed at all
        d.reset();
        ok (d.finish() && d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::TooShort && d.result().windowsExamined == 0, "nothing fed: NotExamined, TooShort");
    }

    group ("non-finite samples are read as zero, counted, and do not poison the report");
    {
        synthetic::Stereo bad = coded;
        const auto g = Detector::geometryFor (44100.0, kStretch);
        const std::uint64_t s0 = Detector::windowStart (g, (std::uint64_t) n, 0), s7 = Detector::windowStart (g, (std::uint64_t) n, 7);
        // inside a stretch and clear of every other stretch's margin, so each is met exactly once: the head of the
        // first stretch (the second one's margin starts 0.05 s in) and the tail of the last (the seventh's ends 0.2 s in)
        bad.left[(std::size_t) s0 + 100u] = std::numeric_limits<float>::quiet_NaN();
        bad.right[(std::size_t) s0 + 1000u] = std::numeric_limits<float>::infinity();
        bad.left[(std::size_t) s7 + 9500u] = -std::numeric_limits<float>::infinity();
        bad.left[10] = std::numeric_limits<float>::quiet_NaN();             // outside every stretch: not this instrument's to count
        const auto d = analyse (bad, 44100.0, lawParams (Depth::Exhaustive));
        if (felitronics::test::run (d != nullptr))
        {
            ok (d->result().nonFiniteSamples == 3, "the three non-finite samples inside the stretches are counted (" + std::to_string (d->result().nonFiniteSamples) + ")");
            ok (d->result().verdict == Verdict::Confirmed && d->result().gridPhase == 480 && d->result().windowsFound == 8, "the grid is still confirmed in all eight stretches");
            bool finite = std::isfinite (d->result().bestScore);
            for (int w = 0; w < 8; ++w)
                for (int h = 0; h < kHypotheses; ++h)
                    finite = finite && std::isfinite (d->window (w).hypotheses[h].reading.score) && std::isfinite (d->window (w).hypotheses[h].reading.second);
            for (int i = 0; i < d->curveOffsets() * kCells; ++i) finite = finite && std::isfinite (d->curve()[(std::size_t) i]);
            ok (finite, "every score and every cell of the curve is finite");
        }
    }

    group ("what prepare() takes is what storageFor() said, and nothing after prepare() allocates");
    {
        for (double rate : { 44100.0, 48000.0 })
            for (double seconds : { 0.5, 2.0 })
            {
                const auto p = params (Depth::Follow, seconds);
                const std::uint64_t total = (std::uint64_t) (rate * 25.0);
                const auto st = Detector::storageFor (rate, 2, total, p);
                auto d = std::make_unique<Detector>();
                d->setParams (p);
                const long long before = alloc::bytes.load();
                const bool prepared = d->prepare (rate, 2, total);
                const long long took = alloc::bytes.load() - before;
                std::printf ("    %.1f kHz, stretches of %.1f s: storageFor says %llu bytes, prepare() asked for %lld\n", rate / 1000.0, seconds, (unsigned long long) st.bytes(), took);
                ok (prepared && st.ok && st.examined && st.windows == 8, "prepared; the storage is that of eight stretches");
               #if defined(_LIBCPP_VERSION)
                ok ((unsigned long long) took == st.bytes(), "a fresh object's prepare() asks the heap for exactly storageFor().bytes()");
               #else
                ok (took > 0, "prepare() allocated  [byte-exact accounting is asserted on libc++]");
               #endif
            }
        // a programme that cannot be examined asks for nothing
        {
            auto d = std::make_unique<Detector>();
            d->setParams (params());
            const long long before = alloc::count.load();
            const bool a = d->prepare (96000.0, 2, 1000000), b = d->prepare (44100.0, 2, 100);
            felitronics::test::okNoAlloc (alloc::count.load() == before, "an unsupported rate and a too-short programme allocate nothing");
            ok (a && b && Detector::storageFor (96000.0, 2, 1000000, params()).bytes() == 0, "...and their storage is zero bytes");
        }
        Detector d;
        d.setParams (lawParams (Depth::Follow));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n), "prepared");
        const long long before = alloc::count.load();
        bool accepted = feed (d, coded, 2, 4096);
        accepted = accepted && d.finish();
        const bool confirmed = d.result().verdict == Verdict::Confirmed;
        (void) d.progress();
        d.reset();
        accepted = accepted && feed (d, original, 2, 0) && d.finish();
        const long long after = alloc::count.load();
        ok (accepted && confirmed && d.result().verdict == Verdict::NoGrid, "two whole analyses ran: a coded programme, reset(), an uncoded one");
        felitronics::test::okNoAlloc (after == before, "process(), finish(), progress() and reset() did not allocate");
    }

    group ("setParams() takes effect at the next prepare(), not before");
    {
        Detector d;
        d.setParams (lawParams (Depth::Verdict));
        ok (d.prepare (44100.0, 2, (std::uint64_t) n), "prepared with depth Verdict");
        d.setParams (params (Depth::Exhaustive, 1.0));
        ok (d.params().depth == Depth::Verdict && std::fabs (d.params().windowSeconds - kStretch) < 1.0e-12, "params() still says what prepare() installed");
        ok (d.process (planes, 2, n) && d.finish() && scans (d) == kHypotheses + 1, "and the analysis ran with it: eight scans");
        ok (d.prepare (44100.0, 2, (std::uint64_t) n) && d.params().depth == Depth::Exhaustive && std::fabs (d.params().windowSeconds - 1.0) < 1.0e-12, "the next prepare() installs the new ones");
    }

    return felitronics::test::report();
}
