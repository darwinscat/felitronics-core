// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Theory-first FALSIFICATION suite for core::DryAligner. Where DryAlignerTests.cpp is example-based
// (hand-picked taps/blocks), this suite derives expectations from the DOCUMENTED CONTRACT ONLY and
// attacks them with an INDEPENDENT reference model + randomized fuzzing, then checks the algebraic
// properties the contract implies (linearity, channel independence, reset==fresh, prepare-clears).
//
// The contract (from DryAligner.h, verbatim intent):
//   * advance() stages scratch = `io` delayed by the CURRENT block's tap D, reading history from a
//     per-channel ring fed EVERY block. The delayed sample is a PURE COPY of the input from exactly D
//     samples earlier — "no arithmetic, no interpolation"; D == 0 is an exact identity copy.
//   * D is clamped to [0, capacity-1]; numChannels to [0, prepared channels]; numSamples to
//     [0, prepared maxBlock]. A kept contract makes every clamp a no-op; a violated one processes the
//     prefix that fits and the cursor advances by the CLAMPED count.
//   * prepare() clears history (warm-up region reads zeros); reset() restores freshly-prepared state.
//   * delayed(ch) is a real copy that survives the caller overwriting `io` in place.
//
// INDEPENDENT REFERENCE (no shared code with DryAligner): per channel, an ever-growing history vector
// indexed by an absolute global sample counter. Output at global time g (current block's tap D) is
// history[g - D], or 0.0f for g - D < 0 (pre-stream silence). This model has NO ring, NO wraparound,
// NO shared cursor — a ring off-by-one / bad wrap / stale-tap / cursor-desync in DryAligner diverges
// from it. It is valid EXACTLY for the contract's canonical usage: every prepared channel fed every
// block (numChannels held == prepared). Count-clamp corners are attacked separately with hand-derived
// expectations, and the variable-numChannels case (which the contract UNDER-SPECIFIES) is documented
// and pinned below.

#include <felitronics_test.h>
#include <felitronics/core/DryAligner.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;

namespace
{
    //==========================================================================
    // Deterministic, portable RNG (own LCG — no <random> distribution, so the exact schedule that
    // trips a failure reproduces bit-for-bit on any toolchain/host).
    struct Rng
    {
        std::uint32_t s;
        explicit Rng (std::uint32_t seed) : s (seed ? seed : 0x9E3779B9u) {}
        std::uint32_t u32() { s = s * 1664525u + 1013904223u; return s; }
        // Uniform-ish int in [lo, hi] inclusive (small ranges → negligible modulo bias, determinism is
        // what matters — both the reference and the DUT see the SAME drawn schedule).
        int range (int lo, int hi) { if (hi <= lo) return lo; return lo + (int) (u32() % (std::uint32_t) (hi - lo + 1)); }
        // A finite float in [-scale, scale). Never inf / nan (bit-exact float== is the oracle here).
        float finiteF (float scale) { return ((float) ((double) u32() / 4294967296.0) - 0.5f) * 2.0f * scale; }
    };

    // Bit-exact float equality — compares the raw 32-bit pattern (so -0.0 != +0.0 and a stray bit in a
    // "pure copy" is caught). The contract promises bit-exactness, so the oracle must be bit-level.
    inline bool bitEq (float a, float b) noexcept
    {
        std::uint32_t x, y; std::memcpy (&x, &a, 4); std::memcpy (&y, &b, 4); return x == y;
    }
    inline std::uint32_t bits (float a) noexcept { std::uint32_t x; std::memcpy (&x, &a, 4); return x; }

    //==========================================================================
    // Independent reference model — infinite per-channel history + explicit global indexing.
    // Models the CANONICAL usage (every prepared channel fed every block). advance() mirrors the
    // DOCUMENTED clamps (D, numChannels, numSamples) and the "cursor advances by the CLAMPED count"
    // rule, but computes the delayed sample by direct history lookup, not a ring.
    struct RefAligner
    {
        int channels = 1, cap = 2, maxBlock = 1;
        long long processed = 0;                       // absolute samples committed per channel (== cursor)
        std::vector<std::vector<float>> hist;          // hist[ch][globalIndex]

        void prepare (int nc, int mb, int c)
        {
            channels = std::max (1, nc);
            maxBlock = std::max (1, mb);
            cap      = std::max (2, c);
            hist.assign ((std::size_t) channels, {});
            processed = 0;
        }

        void reset()
        {
            for (auto& h : hist) h.clear();
            processed = 0;
        }

        // io[ch][0..numSamples) is the raw block. Returns the effective (clamped) sample count; writes
        // the staged delayed copy for every prepared channel into out[ch][0..ns).  REQUIRES numChannels
        // >= channels (canonical "feed all" usage) so the infinite-history model stays exact.
        int advance (const std::vector<std::vector<float>>& io, int numChannels, int numSamples, int delaySamples,
                     std::vector<std::vector<float>>& out)
        {
            const int D  = std::clamp (delaySamples, 0, cap - 1);
            const int nc = std::clamp (numChannels,  0, channels);
            const int ns = std::clamp (numSamples,   0, maxBlock);

            for (int ch = 0; ch < nc; ++ch)
                for (int i = 0; i < ns; ++i)
                    hist[(std::size_t) ch].push_back (io[(std::size_t) ch][(std::size_t) i]);   // commit fed samples

            out.assign ((std::size_t) channels, std::vector<float> ((std::size_t) std::max (0, ns), 0.0f));
            for (int ch = 0; ch < nc; ++ch)
                for (int i = 0; i < ns; ++i)
                {
                    const long long g   = processed + i;
                    const long long idx = g - D;
                    out[(std::size_t) ch][(std::size_t) i] = (idx < 0) ? 0.0f : hist[(std::size_t) ch][(std::size_t) idx];
                }
            processed += ns;   // cursor advances by the CLAMPED count
            return ns;
        }
    };

    //==========================================================================
    // One randomized fuzz schedule: prepare(channels, maxBlock, cap); run a sequence of blocks with
    // per-block ragged sizes and a tap pattern; assert DryAligner's staged scratch == the reference,
    // bit-for-bit, every channel every sample. Returns "" on success, else a human-readable divergence.
    enum class TapMode { Constant, Ramp, RandomJumps };

    std::string runSchedule (std::uint32_t seed)
    {
        Rng rng (seed);

        const int channels = rng.range (1, 8);
        const int cap      = rng.range (2, 512);
        const int maxBlock = rng.range (1, 256);
        const int nBlocks  = rng.range (1, 30);
        const TapMode mode = (TapMode) rng.range (0, 2);
        static const float scales[] = { 1.0e-3f, 1.0f, 1.0e3f, 1.0e6f };
        const float  scale = scales[rng.range (0, 3)];

        core::DryAligner dut; dut.prepare (channels, maxBlock, cap);
        RefAligner       ref; ref.prepare (channels, maxBlock, cap);

        if (dut.capacity() != std::max (2, cap) || dut.numChannels() != std::max (1, channels))
            return "prepare() mis-reported capacity/numChannels";

        const int constD = rng.range (-8, cap + 8);                 // Constant mode: one tap for the schedule
        int rampD = rng.range (0, cap - 1);                         // Ramp mode: swept tap
        int rampStep = rng.range (1, 7) * (rng.range (0, 1) ? 1 : -1);

        for (int b = 0; b < nBlocks; ++b)
        {
            // --- tap for this block (deliberately reaches OUT of range sometimes → exercise the clamp) ---
            int Draw;
            switch (mode)
            {
                case TapMode::Constant:    Draw = constD; break;
                case TapMode::Ramp:        rampD += rampStep; if (rampD < 0 || rampD > cap - 1) { rampStep = -rampStep; rampD += 2 * rampStep; } Draw = rampD; break;
                default:                   { const int pick = rng.range (0, 9);                   // random jumps incl. 0, cap-1, OOB
                                             Draw = pick == 0 ? 0 : pick == 1 ? cap - 1
                                                  : pick == 2 ? rng.range (-8, -1)
                                                  : pick == 3 ? rng.range (cap, cap + 8)
                                                  : rng.range (0, cap - 1); } break;
            }

            // --- block size: mostly valid, sometimes 0 / negative / over-maxBlock (violation) ---
            int nsRaw;
            { int pick = rng.range (0, 11);
              nsRaw = pick == 0 ? 0
                    : pick == 1 ? rng.range (-4, -1)                       // negative → clamps to 0
                    : pick == 2 ? rng.range (maxBlock + 1, maxBlock + 40)  // over-max → clamps to maxBlock
                    : rng.range (1, maxBlock); }
            const int bufLen = std::max (0, nsRaw);

            // --- channels passed: always >= prepared (extras must be ignored → nc clamp) ---
            const int passC = channels + rng.range (0, 3);

            // --- distinct random finite input; feed the SAME buffers to DUT and reference ---
            std::vector<std::vector<float>> buf ((std::size_t) passC);
            for (int ch = 0; ch < passC; ++ch)
            {
                buf[(std::size_t) ch].resize ((std::size_t) bufLen);
                for (int i = 0; i < bufLen; ++i) buf[(std::size_t) ch][(std::size_t) i] = rng.finiteF (scale);
            }
            std::vector<const float*> io ((std::size_t) passC);
            for (int ch = 0; ch < passC; ++ch) io[(std::size_t) ch] = buf[(std::size_t) ch].data();

            dut.advance (io.data(), passC, nsRaw, Draw);

            std::vector<std::vector<float>> refOut;
            const int ns = ref.advance (buf, passC, nsRaw, Draw, refOut);

            for (int ch = 0; ch < channels; ++ch)
                for (int i = 0; i < ns; ++i)
                    if (! bitEq (dut.delayed (ch)[i], refOut[(std::size_t) ch][(std::size_t) i]))
                    {
                        char m[512];
                        std::snprintf (m, sizeof m,
                            "seed=%u ch=%d cap=%d maxBlock=%d block=%d/%d D=%d ns=%d i=%d  dut=%08x ref=%08x",
                            seed, channels, cap, maxBlock, b, nBlocks, Draw, ns, i,
                            bits (dut.delayed (ch)[i]), bits (refOut[(std::size_t) ch][(std::size_t) i]));
                        return std::string (m);
                    }
        }
        return {};
    }
}

int main()
{
    std::printf ("felitronics::core DryAligner THEORY (falsification) tests\n");

    //==========================================================================
    test::group ("fuzz-equivalence vs independent infinite-history reference (thousands of schedules)");
    {
        const int N = 4000;
        std::string firstFail;
        int passed = 0;
        for (int k = 0; k < N; ++k)
        {
            const auto fail = runSchedule (0xA5F00D01u + (std::uint32_t) k * 2654435761u);
            if (fail.empty()) ++passed;
            else { firstFail = fail; break; }
        }
        test::ok (firstFail.empty(),
                  firstFail.empty() ? (std::to_string (passed) + " randomized schedules staged bit-exactly")
                                    : ("staged scratch diverged from reference: " + firstFail));
    }

    //==========================================================================
    // Corner schedules the random draw might under-sample — force the extremes.
    test::group ("forced corners: cap=2, maxBlock=1, blocks >> capacity, cap-1 tap, D=0 identity");
    {
        const std::uint32_t seeds[] = { 1u, 2u, 3u, 7u, 42u, 255u, 4096u, 0xDEADu, 0xBEEFu, 0xFFFFu };
        // These just run more schedules through the same oracle with fixed low-entropy seeds; combined
        // with the 4000 random ones above they blanket cap=2..512 / maxBlock=1..256.
        std::string firstFail;
        for (auto sd : seeds) { auto f = runSchedule (sd); if (! f.empty()) { firstFail = f; break; } }
        test::ok (firstFail.empty(), firstFail.empty() ? "all forced-seed schedules bit-exact" : firstFail);

        // Explicit cap=2 (minimum ring): a D=1 tap on 512-sample blocks (blocks 256x the ring) must still
        // emit x[n-1] for every n>=1 — the tightest wrap the ring can do.
        core::DryAligner a; a.prepare (1, 512, 2);
        RefAligner r; r.prepare (1, 512, 2);
        Rng rng (0xC0FFEEu);
        bool exact = true;
        for (int b = 0; b < 20 && exact; ++b)
        {
            std::vector<std::vector<float>> buf (1); buf[0].resize (512);
            for (auto& v : buf[0]) v = rng.finiteF (1.0f);
            const float* io[1] = { buf[0].data() };
            a.advance (io, 1, 512, 1);
            std::vector<std::vector<float>> ro; r.advance (buf, 1, 512, 1, ro);
            for (int i = 0; i < 512 && exact; ++i) exact = bitEq (a.delayed (0)[i], ro[0][(std::size_t) i]);
        }
        test::ok (exact, "cap=2 ring, 512-sample blocks, D=1: bit-exact against reference");
    }

    //==========================================================================
    // Superposition — the aligner is LINEAR time-varying: under one fixed schedule, staging a·x + b·y
    // equals a·stage(x) + b·stage(y). Value-equality (==) rather than bit-equality, because the warm-up
    // silence is a literal +0.0f on the direct path but a·0 + b·0 on the recombined path (signed-zero
    // only) — the linearity claim is about VALUES, and the arithmetic in the interior is the identical
    // expression a·x[sel] + b·y[sel] evaluated the same way, so it matches exactly there.
    test::group ("superposition: stage(a*x + b*y) == a*stage(x) + b*stage(y)");
    {
        bool ok = true;
        for (int trial = 0; trial < 200 && ok; ++trial)
        {
            Rng rng (0x5150u + (std::uint32_t) trial);
            const int channels = rng.range (1, 4);
            const int cap      = rng.range (2, 256);
            const int maxBlock = rng.range (1, 128);
            const int nBlocks  = rng.range (1, 16);
            const float a = rng.finiteF (2.0f), b = rng.finiteF (2.0f);

            // Pre-draw the whole schedule + inputs so all three runs are IDENTICAL time-variation.
            std::vector<int> nsSeq (nBlocks), dSeq (nBlocks);
            std::vector<std::vector<std::vector<float>>> X (nBlocks), Y (nBlocks);   // [block][ch][i]
            for (int bkt = 0; bkt < nBlocks; ++bkt)
            {
                nsSeq[bkt] = rng.range (0, maxBlock);
                dSeq[bkt]  = rng.range (0, cap - 1);
                X[bkt].assign ((std::size_t) channels, {});
                Y[bkt].assign ((std::size_t) channels, {});
                for (int ch = 0; ch < channels; ++ch)
                {
                    X[bkt][(std::size_t) ch].resize ((std::size_t) nsSeq[bkt]);
                    Y[bkt][(std::size_t) ch].resize ((std::size_t) nsSeq[bkt]);
                    for (int i = 0; i < nsSeq[bkt]; ++i)
                    {
                        X[bkt][(std::size_t) ch][(std::size_t) i] = rng.finiteF (1.0f);
                        Y[bkt][(std::size_t) ch][(std::size_t) i] = rng.finiteF (1.0f);
                    }
                }
            }

            core::DryAligner ax; ax.prepare (channels, maxBlock, cap);   // runs on the combined input
            core::DryAligner sx; sx.prepare (channels, maxBlock, cap);   // runs on x
            core::DryAligner sy; sy.prepare (channels, maxBlock, cap);   // runs on y

            for (int bkt = 0; bkt < nBlocks && ok; ++bkt)
            {
                const int ns = nsSeq[bkt], D = dSeq[bkt];
                std::vector<std::vector<float>> comb ((std::size_t) channels);
                std::vector<const float*> ioc (channels), iox (channels), ioy (channels);
                for (int ch = 0; ch < channels; ++ch)
                {
                    comb[(std::size_t) ch].resize ((std::size_t) ns);
                    for (int i = 0; i < ns; ++i)
                        comb[(std::size_t) ch][(std::size_t) i] = a * X[bkt][(std::size_t) ch][(std::size_t) i]
                                                                + b * Y[bkt][(std::size_t) ch][(std::size_t) i];
                    ioc[(std::size_t) ch] = comb[(std::size_t) ch].data();
                    iox[(std::size_t) ch] = X[bkt][(std::size_t) ch].data();
                    ioy[(std::size_t) ch] = Y[bkt][(std::size_t) ch].data();
                }
                ax.advance (ioc.data(), channels, ns, D);
                sx.advance (iox.data(), channels, ns, D);
                sy.advance (ioy.data(), channels, ns, D);

                for (int ch = 0; ch < channels && ok; ++ch)
                    for (int i = 0; i < ns && ok; ++i)
                    {
                        const float lhs = ax.delayed (ch)[i];
                        const float rhs = a * sx.delayed (ch)[i] + b * sy.delayed (ch)[i];
                        ok = (lhs == rhs);
                    }
            }
        }
        test::ok (ok, "linear time-varying: recombination matches the combined staging");
    }

    //==========================================================================
    // Channel independence — poisoning one channel's input must not perturb another channel's output
    // (separate ring rows; the shared cursor advances identically regardless of sample values).
    test::group ("channel independence: poisoning one channel never changes another's output");
    {
        bool ok = true;
        for (int trial = 0; trial < 100 && ok; ++trial)
        {
            Rng rng (0x0C0Du + (std::uint32_t) trial);
            const int channels = rng.range (2, 6);
            const int cap      = rng.range (2, 256);
            const int maxBlock = rng.range (1, 128);
            const int nBlocks  = rng.range (1, 12);
            const int victim   = rng.range (0, channels - 1);
            const int poison   = (victim + 1 + rng.range (0, channels - 2)) % channels;

            std::vector<int> nsSeq (nBlocks), dSeq (nBlocks);
            std::vector<std::vector<std::vector<float>>> in (nBlocks);
            for (int bkt = 0; bkt < nBlocks; ++bkt)
            {
                nsSeq[bkt] = rng.range (0, maxBlock);
                dSeq[bkt]  = rng.range (-4, cap + 4);
                in[bkt].assign ((std::size_t) channels, {});
                for (int ch = 0; ch < channels; ++ch)
                {
                    in[bkt][(std::size_t) ch].resize ((std::size_t) nsSeq[bkt]);
                    for (int i = 0; i < nsSeq[bkt]; ++i) in[bkt][(std::size_t) ch][(std::size_t) i] = rng.finiteF (1.0f);
                }
            }

            auto runCollectVictim = [&] (bool poisonRun)
            {
                core::DryAligner a; a.prepare (channels, maxBlock, cap);
                std::vector<float> collected;
                for (int bkt = 0; bkt < nBlocks; ++bkt)
                {
                    const int ns = nsSeq[bkt];
                    std::vector<std::vector<float>> buf = in[bkt];
                    if (poisonRun)
                        for (int i = 0; i < ns; ++i) buf[(std::size_t) poison][(std::size_t) i] = rng.finiteF (9.0e5f) + 1.0f;
                    std::vector<const float*> io (channels);
                    for (int ch = 0; ch < channels; ++ch) io[(std::size_t) ch] = buf[(std::size_t) ch].data();
                    a.advance (io.data(), channels, ns, dSeq[bkt]);
                    for (int i = 0; i < ns; ++i) collected.push_back (a.delayed (victim)[i]);
                }
                return collected;
            };

            const auto clean    = runCollectVictim (false);
            const auto poisoned  = runCollectVictim (true);
            if (clean.size() != poisoned.size()) { ok = false; }
            else for (std::size_t i = 0; i < clean.size() && ok; ++i) ok = bitEq (clean[i], poisoned[i]);
        }
        test::ok (ok, "victim channel output is bit-identical whether or not another channel is poisoned");
    }

    //==========================================================================
    // prepare() clears history — the warm-up region reads EXACT zeros, and re-prepare after warming
    // clears again (documented: "prepare clears").
    test::group ("prepare() clears: warm-up region is exact zeros, re-prepare re-clears");
    {
        const int cap = 200, D = 137, blk = 64;
        core::DryAligner a; a.prepare (2, blk, cap);
        Rng rng (0xABCDu);
        // Warm it thoroughly.
        for (int b = 0; b < 10; ++b)
        {
            std::vector<std::vector<float>> buf (2, std::vector<float> (blk));
            for (auto& row : buf) for (auto& v : row) v = rng.finiteF (1.0f);
            const float* io[2] = { buf[0].data(), buf[1].data() };
            a.advance (io, 2, blk, D);
        }
        // Re-prepare must wipe history: first block's first D samples are pre-stream silence == 0.
        a.prepare (2, blk, cap);
        std::vector<std::vector<float>> buf (2, std::vector<float> (blk));
        for (auto& row : buf) for (auto& v : row) v = rng.finiteF (1.0f);
        const float* io[2] = { buf[0].data(), buf[1].data() };
        a.advance (io, 2, blk, D);
        bool zeros = true;
        for (int ch = 0; ch < 2; ++ch)
            for (int i = 0; i < D && i < blk && zeros; ++i)
                zeros = bitEq (a.delayed (ch)[i], 0.0f);
        test::ok (zeros, "post-prepare warm-up staged as exact +0.0f (no stale history survived)");
    }

    //==========================================================================
    // reset() == freshly-prepared — an identical schedule after reset() reproduces the from-fresh run
    // bit-for-bit, and matches a separate instance that was only ever prepared (never warmed).
    test::group ("reset() == freshly-prepared state (bit-identical on a repeated schedule)");
    {
        auto captureRun = [] (core::DryAligner& a, std::uint32_t seed, int channels, int maxBlock, int cap, int nBlocks)
        {
            Rng rng (seed);
            std::vector<float> out;
            for (int b = 0; b < nBlocks; ++b)
            {
                const int ns = rng.range (0, maxBlock);
                const int D  = rng.range (-4, cap + 4);
                std::vector<std::vector<float>> buf ((std::size_t) channels, std::vector<float> ((std::size_t) ns));
                for (auto& row : buf) for (auto& v : row) v = rng.finiteF (1.0f);
                std::vector<const float*> io (channels);
                for (int ch = 0; ch < channels; ++ch) io[(std::size_t) ch] = buf[(std::size_t) ch].data();
                a.advance (io.data(), channels, ns, D);
                for (int ch = 0; ch < channels; ++ch) for (int i = 0; i < ns; ++i) out.push_back (a.delayed (ch)[i]);
            }
            return out;
        };

        const int channels = 3, maxBlock = 100, cap = 300, nBlocks = 25;
        const std::uint32_t sched = 0x7E57u;

        core::DryAligner a; a.prepare (channels, maxBlock, cap);
        const auto firstRun = captureRun (a, sched, channels, maxBlock, cap, nBlocks);
        // Perturb with an UNRELATED schedule, then reset() and replay the ORIGINAL schedule.
        (void) captureRun (a, 0x1234u, channels, maxBlock, cap, nBlocks);
        a.reset();
        const auto afterReset = captureRun (a, sched, channels, maxBlock, cap, nBlocks);

        core::DryAligner fresh; fresh.prepare (channels, maxBlock, cap);   // never warmed
        const auto freshRun = captureRun (fresh, sched, channels, maxBlock, cap, nBlocks);

        bool eqReset = firstRun.size() == afterReset.size();
        for (std::size_t i = 0; i < firstRun.size() && eqReset; ++i) eqReset = bitEq (firstRun[i], afterReset[i]);
        bool eqFresh = firstRun.size() == freshRun.size();
        for (std::size_t i = 0; i < firstRun.size() && eqFresh; ++i) eqFresh = bitEq (firstRun[i], freshRun[i]);

        test::ok (eqReset, "reset() then replay == from-fresh replay (bit-identical)");
        test::ok (eqFresh, "reset()-based state == a never-warmed prepared instance (bit-identical)");
    }

    //==========================================================================
    // delayed() is a real copy — it survives the caller overwriting `io` in place (the whole reason
    // scratch exists). Also proves "no arithmetic": non-finite / signed-zero bit patterns pass through
    // UNTOUCHED (a pure copy preserves the exact bits; any arithmetic would canonicalise a NaN).
    test::group ("staged copy survives buffer overwrite; non-finite values pass through bit-exact");
    {
        const int blk = 16;
        core::DryAligner a; a.prepare (1, blk, 32);
        std::vector<float> live (blk);
        // Adversarial bit patterns: +0, -0, +inf, -inf, a signalling-ish NaN, a quiet NaN, FLT_MAX, tiny.
        const std::uint32_t pat[] = { 0x00000000u, 0x80000000u, 0x7F800000u, 0xFF800000u,
                                      0x7F800001u, 0x7FC00000u, 0x7F7FFFFFu, 0x00000001u,
                                      0x3F800000u, 0xBF800000u, 0x40490FDBu, 0xC0490FDBu,
                                      0x00800000u, 0x807FFFFFu, 0x4B7FFFFFu, 0xCB000000u };
        for (int i = 0; i < blk; ++i) std::memcpy (&live[(std::size_t) i], &pat[i], 4);
        const float* io[1] = { live.data() };
        a.advance (io, 1, blk, 0);                                   // D=0 identity
        std::fill (live.begin(), live.end(), -99.0f);                // caller trashes the live buffer
        bool intact = true;
        for (int i = 0; i < blk && intact; ++i) intact = (bits (a.delayed (0)[i]) == pat[i]);
        test::ok (intact, "scratch preserves the exact pre-overwrite bits — incl. NaN/inf/-0 (no arithmetic)");
    }

    //==========================================================================
    // Contract-violation clamps, hand-derived (NOT via the fuzz oracle) so the documented behavior is
    // legible: numSamples > maxBlock stages the prefix that fits AND commits the cursor by the clamped
    // count; numChannels > prepared processes the channel prefix and ignores the extras.
    test::group ("clamp semantics: prefix staged, cursor commits by the CLAMPED count");
    {
        // (1) numSamples > maxBlock. prepare maxBlock=8, cap=64. Feed a 30-sample block at D=0 → only the
        //     first 8 are staged. Then a second block at D=8 must read back exactly that 8-sample prefix,
        //     PROVING the cursor advanced by 8 (the clamp), not by 30.
        core::DryAligner a; a.prepare (1, 8, 64);
        std::vector<float> big (30), nxt (8);
        for (int i = 0; i < 30; ++i) big[(std::size_t) i] = (float) (i + 1);       // 1..30, distinct
        for (int i = 0; i < 8;  ++i) nxt[(std::size_t) i] = (float) (-i - 1);       // distinct, disjoint
        const float* io1[1] = { big.data() };
        a.advance (io1, 1, 30, 0);
        bool prefix = true;
        for (int i = 0; i < 8 && prefix; ++i) prefix = bitEq (a.delayed (0)[i], (float) (i + 1));
        test::ok (prefix, "numSamples>maxBlock: first maxBlock samples staged bit-exactly (identity)");

        const float* io2[1] = { nxt.data() };
        a.advance (io2, 1, 8, 8);                                                   // D=8 reaches back a full block
        bool commit = true;
        for (int i = 0; i < 8 && commit; ++i) commit = bitEq (a.delayed (0)[i], (float) (i + 1));
        test::ok (commit, "cursor advanced by the CLAMPED count (8): D=8 reads back the staged prefix");

        // (2) numChannels > prepared. prepare 2 ch; pass 5 → extras ignored, both real channels staged.
        core::DryAligner c; c.prepare (2, 8, 32);
        std::vector<float> ch0 (8), ch1 (8), extra (8, 0.0f);
        for (int i = 0; i < 8; ++i) { ch0[(std::size_t) i] = (float) (i + 100); ch1[(std::size_t) i] = (float) (i + 200); }
        const float* io3[5] = { ch0.data(), ch1.data(), extra.data(), extra.data(), extra.data() };
        c.advance (io3, 5, 8, 0);
        bool twoOk = true;
        for (int i = 0; i < 8 && twoOk; ++i) twoOk = bitEq (c.delayed (0)[i], ch0[(std::size_t) i]) && bitEq (c.delayed (1)[i], ch1[(std::size_t) i]);
        test::ok (twoOk, "numChannels>prepared: the prepared-channel prefix is staged, extras ignored");
    }

    //==========================================================================
    // FINDING — contract UNDER-SPECIFICATION, pinned. The doc says the ring is "fed EVERY block (so it
    // stays warm)" yet also permits "numChannels <= prepared". If a caller actually passes FEWER channels
    // on some block, the SHARED cursor still advances, so the skipped channel's ring row develops a hole:
    // its later reads land on positions that were never written for it (pre-warm zeros) rather than on its
    // own delayed samples. i.e. per-channel alignment is only well-defined when numChannels is held
    // CONSTANT across advance() calls — the "<= prepared" wording under-specifies the variable case.
    // This is a DOC-PRECISION finding, not a valid-path defect (constant-numChannels usage is bit-exact,
    // proven above). We PIN the actual behavior so it is regression-locked and any change is deliberate.
    test::group ("[finding/pinned] variable numChannels desyncs the skipped channel (contract under-specifies)");
    {
        core::DryAligner a; a.prepare (2, 4, 8);
        // Block 1: feed ONLY channel 0 (numChannels=1). Channel 1 is NOT fed; cursor still advances by 4.
        std::vector<float> c0 (4), c1 (4);
        for (int i = 0; i < 4; ++i) { c0[(std::size_t) i] = (float) (i + 1); c1[(std::size_t) i] = (float) (i + 11); }
        const float* io1[2] = { c0.data(), c1.data() };
        a.advance (io1, 1, 4, 0);                       // nc clamps to... 1 (only ch0 written)
        // Block 2: feed BOTH channels at D=4 — a full-block reach-back.
        std::vector<float> d0 (4), d1 (4);
        for (int i = 0; i < 4; ++i) { d0[(std::size_t) i] = (float) (i + 21); d1[(std::size_t) i] = (float) (i + 31); }
        const float* io2[2] = { d0.data(), d1.data() };
        a.advance (io2, 2, 4, 4);
        // Channel 0 WAS fed in block 1, so its D=4 read returns block-1's samples (1..4) — aligned.
        bool ch0Aligned = true;
        for (int i = 0; i < 4 && ch0Aligned; ++i) ch0Aligned = bitEq (a.delayed (0)[i], (float) (i + 1));
        // Channel 1 was NOT fed in block 1 → its ring row is still zero there → D=4 reads +0.0f, NOT its
        // "input from 4 samples ago" (which never existed). This is the desync we are documenting.
        bool ch1ReadsZeros = true;
        for (int i = 0; i < 4 && ch1ReadsZeros; ++i) ch1ReadsZeros = bitEq (a.delayed (1)[i], 0.0f);
        test::ok (ch0Aligned,   "fed channel stays aligned across the variable-numChannels boundary");
        test::ok (ch1ReadsZeros, "[pinned] skipped channel reads pre-warm zeros — the under-specified case");
    }

    return test::report();
}
