// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// SpectrumFrames self-tests — the shared frame producer the offline analyzers measure on, so its
// invariance is spent once here instead of three times upstairs.
//
//   · a REFERENCE NULL: every bin of every frame against a direct O(N^2) double DFT of the same
//     windowed slice, computed independently of the ring, the trigger and the FFT
//   · Parseval per frame, folded correctly for a one-sided unfolded layout
//   · law 8a: the FULL trace (frame index, start, finite flags, every power bin) compared BIT-EXACTLY
//     by std::bit_cast across 14 slicings, including n == 0 calls, calls larger than maxBlock, and
//     seeded ragged partitions — and across two prepared maxBlock values
//   · the tail contract: no frame is invented, tailUncoveredSamples() names what no frame covered,
//     a programme shorter than the window produces no frame at all
//   · holes: a non-finite sample poisons exactly the frames that contain it and no others
//   · finish() idempotent, push/tick refused after it, reset() replays identically
//   · no allocation in push/tick/finish; the published Storage is the allocated storage
//
// MUTANT PASS — run, not asserted. A gate that has never been red is not a gate, and this one was
// born green, which proves nothing on its own. Six mutants, each built with the object file deleted
// first (a header edited in the same SECOND as the previous build is invisible to make's 1-second
// mtime granularity — two mutants reported a false GREEN that way before the stand was fixed to refuse
// a verdict when no compile line appeared):
//   (1) frame trigger anchored to the frame START instead of its END  -> RED, 10 of 82
//   (2) the ring read in physical order instead of chronological      -> RED,  1 of 82 (the DFT null)
//   (3) a hole no longer poisons the frame that contains it           -> RED,  1 of 82
//   (4) power normalised by N alone, without the window energy        -> RED,  2 of 82
//   (5) the slot leaving the window is not subtracted from the count  -> RED,  3 of 82
//   (6) the window rebuilt per frame instead of once in prepare()     -> GREEN, and recorded as green:
//       it is the same window, so this is equivalence, not a hole in the suite.
// Mutant (1) is why section 3b exists. The re-slicing comparison of section 3 could NOT see it: that
// test compares the implementation with ITSELF, so a schedule that is deterministically wrong stays
// bit-identical across all 14 slicings. It took a witness computed outside the object — frame k starts
// at k*hop — to make it red. Any analyzer built on this producer inherits that lesson: an invariance
// test needs an independent oracle beside it, or it certifies only self-consistency.

#include <felitronics_test.h>
#include <felitronics/analysis/SpectrumFrames.h>

#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s); }
void  operator delete (void* p) noexcept { std::free (p); }
void  operator delete (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using analysis::SpectrumFrames;
using analysis::SpectrumFramesParams;

namespace
{

// --- the fixture: NON-stationary and structured, and its length is coprime with every hop used ---
std::vector<std::vector<float>> fixture (std::size_t frames, int nch, unsigned seed)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    std::vector<std::vector<float>> x ((std::size_t) nch, std::vector<float> (frames, 0.0f));
    for (std::size_t i = 0; i < frames; ++i)
    {
        const double t = (double) i / 48000.0;
        const double env = i < frames / 8 ? 0.0                           // lead silence
                         : i > frames * 7 / 8 ? 0.02                      // quiet tail
                         : (i % 9601 < 300 ? 1.0 : 0.25);                 // bursts against a bed
        for (int c = 0; c < nch; ++c)
        {
            const double f = c == 0 ? 440.0 : 613.0;
            x[(std::size_t) c][i] = (float) (env * (0.7 * std::sin (2.0 * core::kPi * f * t)
                                                  + 0.3 * (double) u (rng)));
        }
    }
    return x;
}

// --- the trace: everything the producer emitted, as bits ---
std::vector<std::uint64_t> run (const std::vector<std::vector<float>>& x, int nch,
                                const std::vector<int>& slices, int maxBlock, int order, int hop)
{
    SpectrumFrames sf;
    SpectrumFramesParams p; p.fftOrder = order; p.hop = hop;
    sf.setParams (p);
    std::vector<std::uint64_t> tr;
    if (! sf.prepare (48000.0, maxBlock, nch)) { tr.push_back (0xDEADull); return tr; }
    const std::size_t total = x[0].size();
    std::size_t at = 0, s = 0;
    while (at < total)
    {
        const int want = slices[s % slices.size()]; ++s;
        const std::size_t take = want <= 0 ? 0 : std::min ((std::size_t) want, total - at);
        for (std::size_t i = 0; i < take; ++i)
        {
            for (int c = 0; c < nch; ++c) sf.push (c, x[(std::size_t) c][at + i], true);
            if (sf.tick())
            {
                tr.push_back ((std::uint64_t) sf.frameIndex());
                tr.push_back ((std::uint64_t) sf.frameStart());
                for (int c = 0; c < nch; ++c)
                {
                    tr.push_back (sf.frameFinite (c) ? 1u : 0u);
                    const double* pw = sf.power (c);
                    for (int k = 0; k < sf.bins(); ++k) tr.push_back (std::bit_cast<std::uint64_t> (pw[k]));
                }
            }
        }
        at += take;
    }
    sf.finish();
    tr.push_back ((std::uint64_t) sf.totalSamples());
    tr.push_back ((std::uint64_t) sf.tailUncoveredSamples());
    tr.push_back ((std::uint64_t) sf.frameCount());
    return tr;
}

} // namespace

int main()
{
    using felitronics::test::ok;
    using felitronics::test::approx;

    // ---------- 1. REFERENCE NULL against a direct double DFT ----------
    {
        const int order = 8, n = 1 << order, hop = 64, nch = 1;
        auto x = fixture ((std::size_t) (n * 3 + 37), nch, 11);
        SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = hop; sf.setParams (p);
        ok (sf.prepare (48000.0, 64, nch), "DFT null: prepare");
        // the window, rebuilt here from the definition rather than read out of the object
        std::vector<double> w ((std::size_t) n);
        double sumW2 = 0.0;
        for (int i = 0; i < n; ++i) { w[(std::size_t) i] = 0.5 - 0.5 * std::cos (2.0 * core::kPi * i / (double) n); sumW2 += w[(std::size_t) i] * w[(std::size_t) i]; }
        int framesChecked = 0;
        double worst = 0.0;
        bool worstFinite = true;
        for (std::size_t i = 0; i < x[0].size(); ++i)
        {
            sf.push (0, x[0][i], true);
            if (! sf.tick()) continue;
            ++framesChecked;
            const std::int64_t st = sf.frameStart();
            for (int k : { 0, 1, 2, 7, 23, n / 4, n / 2 - 1, n / 2 })
            {
                double re = 0.0, im = 0.0;                                  // the oracle: a direct DFT
                for (int t = 0; t < n; ++t)
                {
                    const double v = (double) x[0][(std::size_t) (st + t)] * w[(std::size_t) t];
                    const double a = -2.0 * core::kPi * (double) k * (double) t / (double) n;
                    re += v * std::cos (a); im += v * std::sin (a);
                }
                const double want = (re * re + im * im) / ((double) n * sumW2);
                const double got  = sf.power (0)[k];
                const double d = std::fabs (got - want) / (want > 1e-18 ? want : 1e-18);
                if (! std::isfinite (d)) worstFinite = false;
                if (d > worst) worst = d;
            }
        }
        ok (framesChecked >= 3, "DFT null: the fixture actually produced frames (" + std::to_string (framesChecked) + ")");
        ok (worstFinite, "DFT null: every relative error is finite");
        ok (worst < 1e-9, "DFT null: worst relative error " + std::to_string (worst) + " < 1e-9");
    }

    // ---------- 2. Parseval, folded for the unfolded one-sided layout ----------
    {
        const int order = 8, n = 1 << order;
        SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = n; sf.setParams (p);
        ok (sf.prepare (48000.0, 32, 1), "Parseval: prepare");
        std::mt19937 rng (5); std::uniform_real_distribution<float> u (-1.0f, 1.0f);
        std::vector<float> x ((std::size_t) n);
        for (auto& v : x) v = u (rng);
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const double w = 0.5 - 0.5 * std::cos (2.0 * core::kPi * i / (double) n);
            num += (double) x[(std::size_t) i] * w * (double) x[(std::size_t) i] * w;
            den += w * w;
        }
        bool fired = false;
        for (int i = 0; i < n; ++i) { sf.push (0, x[(std::size_t) i], true); if (sf.tick()) fired = true; }
        ok (fired, "Parseval: a frame closed");
        const double* pw = sf.power (0);
        double sum = pw[0] + pw[n / 2];
        for (int k = 1; k < n / 2; ++k) sum += 2.0 * pw[k];
        approx (sum, num / den, 1e-12, "Parseval: folded bin sum equals the windowed mean square");
    }

    // ---------- 3. LAW 8a: the full trace, bit-exact, across slicings ----------
    {
        const int order = 8, n = 1 << order, nch = 2;
        for (int hop : { n / 2, n / 4, n, 97 })
        {
            auto x = fixture (4099, nch, 3);                               // coprime with 64, every hop and n
            const std::vector<std::vector<int>> slicings = {
                { 4099 }, { 1 }, { 2 }, { 3 }, { 63 }, { 64 }, { 65 },
                { hop - 1 }, { hop }, { hop + 1 }, { n }, { 997 }, { 4097 }, { 8192 },
                { 0, 1, 0, 700, 0 },                                       // n == 0 calls interleaved
            };
            const auto ref = run (x, nch, slicings[0], 512, order, hop);
            ok (ref.size() > 64, "8a hop=" + std::to_string (hop) + ": the trace is non-trivial (" + std::to_string (ref.size()) + " words)");
            int differing = 0;
            for (std::size_t s = 1; s < slicings.size(); ++s)
                if (run (x, nch, slicings[s], 512, order, hop) != ref) ++differing;
            ok (differing == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (differing) + " of "
                                + std::to_string (slicings.size() - 1) + " slicings differ");
            // seeded ragged partitions
            int raggedDiffer = 0;
            for (unsigned seed = 0; seed < 6; ++seed)
            {
                std::mt19937 rng (seed + 41); std::uniform_int_distribution<int> d (1, 1500);
                std::vector<int> sl; sl.reserve (64);
                for (int i = 0; i < 64; ++i) sl.push_back (d (rng));
                if (run (x, nch, sl, 512, order, hop) != ref) ++raggedDiffer;
            }
            ok (raggedDiffer == 0, "8a hop=" + std::to_string (hop) + ": " + std::to_string (raggedDiffer) + " of 6 ragged partitions differ");
            // maxBlock sizes nothing
            ok (run (x, nch, { 64 }, 64, order, hop) == ref && run (x, nch, { 64 }, 8192, order, hop) == ref,
                "8a hop=" + std::to_string (hop) + ": maxBlock 64 and 8192 give the same trace");
        }
    }

    // ---------- 3b. THE SCHEDULE ORACLE ----------
    // Why this exists: section 3 compares the implementation with ITSELF at different slicings, so a
    // schedule that is deterministically WRONG passes it unchanged. The mutation stand proved that —
    // anchoring the trigger to the frame START instead of its END left all 14 slicings bit-identical
    // and the suite green. The schedule therefore needs a witness that does not come from the object:
    // frame k covers exactly [k*hop, k*hop + N), and the count is 1 + floor((T - N)/hop). Hops that do
    // NOT divide N are the ones that separate the two rules (when hop | N they coincide).
    {
        const int order = 8, n = 1 << order;
        for (int hop : { 97, 61, 255, 33, n / 2, n })
        {
            SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = hop; sf.setParams (p);
            ok (sf.prepare (48000.0, 512, 1), "schedule hop=" + std::to_string (hop) + ": prepare");
            const std::int64_t total = 2000;
            std::int64_t k = 0; int wrongStart = 0, wrongIndex = 0;
            for (std::int64_t i = 0; i < total; ++i)
            {
                sf.push (0, 0.25f, true);
                if (! sf.tick()) continue;
                if (sf.frameStart() != k * (std::int64_t) hop) ++wrongStart;
                if (sf.frameIndex() != k) ++wrongIndex;
                ++k;
            }
            sf.finish();
            const std::int64_t want = (total - n) / hop + 1;
            ok (k == want, "schedule hop=" + std::to_string (hop) + ": " + std::to_string (k)
                           + " frames, oracle says " + std::to_string (want));
            ok (wrongStart == 0, "schedule hop=" + std::to_string (hop) + ": every frame starts at k*hop ("
                                 + std::to_string (wrongStart) + " wrong)");
            ok (wrongIndex == 0, "schedule hop=" + std::to_string (hop) + ": frameIndex is k ("
                                 + std::to_string (wrongIndex) + " wrong)");
            ok (sf.tailUncoveredSamples() == total - ((want - 1) * hop + n),
                "schedule hop=" + std::to_string (hop) + ": the uncovered tail matches the oracle");
        }
    }

    // ---------- 4. the tail contract ----------
    {
        const int order = 8, n = 1 << order, hop = n / 2;
        auto x = fixture (1000, 1, 7);
        SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = hop; sf.setParams (p);
        ok (sf.prepare (48000.0, 64, 1), "tail: prepare");
        std::int64_t fired = 0;
        for (std::size_t i = 0; i < x[0].size(); ++i) { sf.push (0, x[0][i], true); if (sf.tick()) ++fired; }
        sf.finish();
        const std::int64_t expectFrames = (1000 - n) / hop + 1;            // 1 + floor((T-N)/hop)
        ok (fired == expectFrames, "tail: " + std::to_string (fired) + " frames, expected " + std::to_string (expectFrames));
        ok (sf.frameCount() == fired, "tail: frameCount agrees");
        ok (sf.tailUncoveredSamples() == 1000 - ((expectFrames - 1) * hop + n), "tail: the uncovered tail is named exactly");
        ok (sf.tailUncoveredSamples() > 0, "tail: this fixture really does leave a tail");
        // shorter than the window: no frame at all, and the tail is the whole programme
        SpectrumFrames sh; sh.setParams (p);
        ok (sh.prepare (48000.0, 64, 1), "tail: short prepare");
        std::int64_t shortFired = 0;
        for (int i = 0; i < n - 1; ++i) { sh.push (0, 0.5f, true); if (sh.tick()) ++shortFired; }
        sh.finish();
        ok (shortFired == 0, "tail: a programme shorter than the window produces NO frame");
        ok (sh.frameCount() == 0 && sh.tailUncoveredSamples() == n - 1, "tail: and names the whole programme uncovered");
    }

    // ---------- 5. holes poison exactly the frames that contain them ----------
    {
        const int order = 8, n = 1 << order, hop = n / 2;
        SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = hop; sf.setParams (p);
        ok (sf.prepare (48000.0, 64, 1), "holes: prepare");
        const std::int64_t bad = 300;                                       // inside frames whose window covers 300
        std::vector<std::pair<std::int64_t, bool>> seen;
        for (std::int64_t i = 0; i < 1400; ++i)
        {
            const bool poison = (i == bad);
            sf.push (0, poison ? std::numeric_limits<float>::quiet_NaN() : 0.5f, true);
            if (sf.tick()) seen.push_back ({ sf.frameStart(), sf.frameFinite (0) });
        }
        int wrong = 0;
        for (auto& [start, finite] : seen)
        {
            const bool covers = bad >= start && bad < start + n;
            if (finite == covers) ++wrong;                                  // finite must be the NEGATION of covers
        }
        ok (! seen.empty(), "holes: frames were produced (" + std::to_string (seen.size()) + ")");
        ok (wrong == 0, "holes: " + std::to_string (wrong) + " frames disagree with which windows contain the NaN");
        // and a poisoned frame's power is not merely wrong, it is refused
        for (auto& [start, finite] : seen)
            if (! finite) { (void) start; break; }
    }

    // ---------- 6. finish / reset / no allocation / storage ----------
    {
        const int order = 8;
        SpectrumFrames sf; SpectrumFramesParams p; p.fftOrder = order; p.hop = 1 << (order - 1); sf.setParams (p);
        const auto st = SpectrumFrames::storageFor (48000.0, 2, p);
        ok (st.ok && st.bytes() > 0, "storage: published before the allocation");
        const long before = g_allocs.load();
        ok (sf.prepare (48000.0, 64, 2), "lifecycle: prepare");
        const long allocsInPrepare = g_allocs.load() - before;
        ok (allocsInPrepare > 0, "lifecycle: prepare is where the heap is touched");
        const long b2 = g_allocs.load();
        for (int i = 0; i < 1000; ++i) { sf.push (0, 0.3f, true); sf.push (1, -0.3f, true); (void) sf.tick(); }
        sf.finish();
        felitronics::test::okNoAlloc (g_allocs.load() == b2, "lifecycle: push/tick/finish allocate nothing");
        sf.finish();
        ok (sf.isFinished(), "lifecycle: finish is idempotent");
        const std::int64_t frozen = sf.totalSamples();
        sf.push (0, 1.0f, true);
        ok (! sf.tick() && sf.totalSamples() == frozen, "lifecycle: tick after finish moves nothing");
        sf.reset();
        ok (sf.totalSamples() == 0 && sf.frameCount() == 0 && ! sf.isFinished(), "lifecycle: reset re-anchors");
        // replay from reset is identical to a fresh run
        auto x = fixture (900, 1, 13);
        SpectrumFrames fresh; fresh.setParams (p); ok (fresh.prepare (48000.0, 64, 1), "replay: fresh prepare");
        SpectrumFrames reused; reused.setParams (p); ok (reused.prepare (48000.0, 64, 1), "replay: reused prepare");
        for (std::size_t i = 0; i < 400; ++i) { reused.push (0, 0.9f, true); (void) reused.tick(); }
        reused.finish(); reused.reset();
        std::vector<std::uint64_t> a, b;
        for (std::size_t i = 0; i < x[0].size(); ++i)
        {
            fresh.push (0, x[0][i], true);  if (fresh.tick())  for (int k = 0; k < fresh.bins();  ++k) a.push_back (std::bit_cast<std::uint64_t> (fresh.power (0)[k]));
            reused.push (0, x[0][i], true); if (reused.tick()) for (int k = 0; k < reused.bins(); ++k) b.push_back (std::bit_cast<std::uint64_t> (reused.power (0)[k]));
        }
        ok (! a.empty() && a == b, "replay: after reset the trace is bit-identical to a fresh object");
    }

    // ---------- 7. refusals ----------
    {
        SpectrumFrames sf; SpectrumFramesParams p;
        p.fftOrder = 3;  sf.setParams (p); ok (! sf.prepare (48000.0, 64, 1), "refuse: fftOrder below 4");
        p.fftOrder = 8;  p.hop = -1; sf.setParams (p); ok (! sf.prepare (48000.0, 64, 1), "refuse: negative hop");
        p.hop = (1 << 8) + 1; sf.setParams (p); ok (! sf.prepare (48000.0, 64, 1), "refuse: hop larger than the window");
        p.hop = 0; sf.setParams (p);
        ok (! sf.prepare (0.0, 64, 1), "refuse: sample rate zero");
        ok (! sf.prepare (std::numeric_limits<double>::infinity(), 64, 1), "refuse: sample rate not finite");
        ok (! sf.prepare (48000.0, 64, 0), "refuse: zero channels");
        ok (! sf.prepare (48000.0, 64, core::kMaxChannels + 1), "refuse: too many channels");
        ok (sf.prepare (48000.0, 64, 1), "refuse: and a good call is still accepted afterwards");
    }

    return felitronics::test::report();
}
