// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore::Probe is the measurement body shared VERBATIM between the native reference CLI and the wasm shim,
// so the whole P0 parity claim rests on it. What that claim needs, and what this suite pins:
//
//   * the TU is genuinely built without FP contraction — the acceptance condition is a compiler flag, and a
//     flag nobody tests is a flag that silently disappears;
//   * the refactor into Probe changed no arithmetic — nulled bit-for-bit against the hand-rolled path the
//     tool used before;
//   * chunking cannot move a bit, at any chunk size, because the wasm side chunks differently from the CLI;
//   * the true-peak CONFIG is pinned, because the core holds a second, different true-peak filter and
//     reaching for it would fail the spike by ~3e-3 dB while looking like a wasm bug;
//   * the adversarial surface (rates, channel counts, degenerate lengths, non-finite samples) behaves the way
//     it is DOCUMENTED to behave, not the way one hopes.
//
// Bit-exact throughout where bit-exactness is the actual claim: test::approx would hide precisely the class
// of defect this suite exists to catch.

#include <felitronics_test.h>

#include "fcore_probe.h"

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/core/Math.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// RT-safety witness: every allocation in this binary is counted (the LoudnessMeter suite's pattern).
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

namespace
{
    // Bit identity, not near-equality: two builds either computed the same double or they did not.
    bool sameBits (double a, double b) noexcept
    {
        std::uint64_t x, y;
        std::memcpy (&x, &a, sizeof x);
        std::memcpy (&y, &b, sizeof y);
        return x == y;
    }

    using Planar = std::vector<std::vector<float>>;

    std::vector<const float*> ptrs (const Planar& p)
    {
        std::vector<const float*> v (p.size());
        for (std::size_t c = 0; c < p.size(); ++c) v[c] = p[c].data();
        return v;
    }

    // A deterministic, reproducible program with transients (so the true-peak path has something to find),
    // a sustained tone (so gating blocks differ) and a quiet stretch (so the relative gate has work).
    Planar makeProgram (double sr, int nc, double seconds)
    {
        const long long n = (long long) std::llround (seconds * sr);
        Planar p ((std::size_t) nc, std::vector<float> ((std::size_t) n, 0.0f));
        std::uint32_t rng = 0x9E3779B9u;                                   // xorshift; no <random>, no locale, no surprises
        for (long long i = 0; i < n; ++i)
        {
            rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
            const double t     = (double) i / sr;
            const double quiet = (t > seconds * 0.6) ? 0.02 : 1.0;         // a stretch the relative gate drops
            const double tone  = 0.25 * std::sin (2.0 * core::kPi * 440.0 * t);
            const double noise = 0.02 * ((double) (rng >> 8) / 8388608.0 - 1.0);
            const double click = ((i % (long long) std::llround (sr * 0.5)) == 0) ? 0.7 : 0.0;   // inter-sample peaks
            for (int c = 0; c < nc; ++c)
                p[(std::size_t) c][(std::size_t) i] =
                    (float) (quiet * (tone * (c == 0 ? 1.0 : 0.93) + noise) + click);
        }
        return p;
    }

    // The comparison surface: everything a cross-toolchain diff looks at.
    struct Surface
    {
        std::vector<double> blocks;
        double              tp = 0.0;
        bool                operator== (const Surface& o) const
        {
            if (blocks.size() != o.blocks.size() || ! sameBits (tp, o.tp)) return false;
            for (std::size_t i = 0; i < blocks.size(); ++i) if (! sameBits (blocks[i], o.blocks[i])) return false;
            return true;
        }
    };

    Surface surfaceOf (const fcore::Probe& p)
    {
        Surface s;
        s.tp = p.truePeakLinear();
        const auto e = p.gatingBlockEnergies();
        s.blocks.assign (e.begin(), e.end());
        return s;
    }

    // Runs the program through a fresh Probe, handing it `chunk` frames at a time (0 = the whole buffer).
    // `drain` mirrors how the real tools run — they always finish(); the refactor-null test turns it off so it
    // can compare against the pre-refactor arithmetic, which had no drain.
    Surface runChunked (const Planar& prog, double sr, int nc, long long chunk, bool drain = true)
    {
        fcore::Probe p;
        p.prepare (sr, nc);
        const auto  base = ptrs (prog);
        const long long n = (long long) prog[0].size();
        if (chunk <= 0) { p.process (base.data(), nc, n); }
        else
        {
            std::vector<const float*> view ((std::size_t) nc);
            for (long long off = 0; off < n; off += chunk)
            {
                const long long m = std::min (chunk, n - off);
                for (int c = 0; c < nc; ++c) view[(std::size_t) c] = prog[(std::size_t) c].data() + off;
                p.process (view.data(), nc, m);
            }
        }
        if (drain) p.finish();
        return surfaceOf (p);
    }
}

int main()
{
    std::printf ("fcore::Probe — the shared P0 measurement body\n");

    // --- The acceptance condition is a COMPILER FLAG. Assert it took effect, or it will vanish in a refactor
    //     and the parity check will start failing for a reason nobody connects to the build.
    //
    //     WHERE THIS ACTUALLY GUARDS. It asserts the weaker, true thing — "no contraction happened here" —
    //     because that is all a program can observe. On arm64, `fmadd` is in the baseline ISA, so dropping
    //     the flag makes this fail immediately: the dev machines are covered. On BASELINE x86-64 there is no
    //     FMA instruction at all, so the compiler cannot contract even at gcc's -ffp-contract=fast default
    //     and this group passes with or without the flag — verified. Add -march=native (or anything else
    //     enabling FMA) on the same x86 machine and it fails correctly, also verified.
    //     So: a default-flags Linux CI job does NOT verify the flag. That is not a hole in the test, it is a
    //     hole in what an x86 CI job can see, and it is why the flag is set explicitly in CMake rather than
    //     left to be caught by tests. ---
    test::group ("build contract: this TU does not contract a*b+c into an FMA");
    {
        // (1+2⁻²⁷)² = 1 + 2⁻²⁶ + 2⁻⁵⁴. Rounded to double that is 1 + 2⁻²⁶ (2⁻⁵⁴ is below half an ulp), so
        // mul-then-add against −(1+2⁻²⁶) is EXACTLY zero — while a single fused rounding keeps the 2⁻⁵⁴.
        volatile double a = 1.0 + std::ldexp (1.0, -27);
        volatile double b = 1.0 + std::ldexp (1.0, -27);
        volatile double c = -(1.0 + std::ldexp (1.0, -26));
        const double separate = a * b + c;
        const double fused    = std::fma ((double) a, (double) b, (double) c);
        test::ok (fused == std::ldexp (1.0, -54), "the constants really do separate a fused result from a rounded one");
        test::ok (separate == 0.0, "a*b+c was rounded twice — no contraction happened in this TU");
        test::ok (! sameBits (separate, fused),
                  "so wasm parity holds; where the target CAN fuse (arm64, or x86 with FMA enabled) this is "
                  "exactly the assertion that dies when -ffp-contract=off is dropped");
    }

    // --- Probe replaced a hand-rolled loop in fcore_measure.cpp. It must have changed NOTHING. ---
    test::group ("refactor null: Probe reproduces the hand-rolled path bit-for-bit");
    {
        const double sr = 48000.0; const int nc = 2;
        const Planar prog = makeProgram (sr, nc, 6.0);
        const long long n = (long long) prog[0].size();

        // the arithmetic the tool ran before fcore_probe.h existed
        analysis::LoudnessMeter lm; lm.prepare (sr, nc, 4.0 * 3600.0);
        std::vector<oversampling::PolyphaseOversampler> os ((std::size_t) nc);
        for (auto& o : os) o.prepare (4, 1, 32);
        std::vector<float> osbuf ((std::size_t) fcore::Probe::kChunk * 4);
        double maxTp = 0.0;
        for (long long off = 0; off < n; off += fcore::Probe::kChunk)
        {
            const int m = (int) std::min<long long> (fcore::Probe::kChunk, n - off);
            const float* view[2] { prog[0].data() + off, prog[1].data() + off };
            lm.process (view, nc, m);
            for (int c = 0; c < nc; ++c)
            {
                const float* in[1] { view[c] }; float* out[1] { osbuf.data() };
                os[(std::size_t) c].upsample (in, 1, m, out);
                for (int k = 0; k < m * 4; ++k) maxTp = std::max (maxTp, (double) std::fabs (osbuf[(std::size_t) k]));
            }
        }
        // The pre-refactor tool had neither the drain nor the sample-peak floor — those are deliberate
        // additions, tested in their own groups. Apply the floor here so this null compares like with like:
        // what it must prove is that the block energies and the oversampler's own maximum did not move.
        double legacySamplePeak = 0.0;
        for (int c = 0; c < nc; ++c)
            for (float v : prog[(std::size_t) c]) legacySamplePeak = std::max (legacySamplePeak, (double) std::fabs (v));

        Surface legacy; legacy.tp = std::max (maxTp, legacySamplePeak);
        const auto le = lm.gatingBlockEnergies();
        legacy.blocks.assign (le.begin(), le.end());

        const Surface viaProbe = runChunked (prog, sr, nc, 0, /*drain*/ false);
        test::ok (legacy.blocks.size() > 40, "the fixture produced a meaningful number of gating blocks");
        test::ok (viaProbe == legacy, "every pre-gate block energy and the true-peak maximum are bit-identical");
    }

    // --- The CLI streams in 8192-frame reads; the wasm shim hands over a whole heap buffer. If chunking moved
    //     a single bit, the parity test would be measuring the harness instead of the toolchains. ---
    test::group ("chunk invariance is BIT-exact, not approximate");
    {
        const double sr = 48000.0; const int nc = 2;
        const Planar prog = makeProgram (sr, nc, 5.0);
        const Surface whole = runChunked (prog, sr, nc, 0);
        for (const long long chunk : { 1LL, 7LL, 999LL, 4096LL, 8192LL, 8193LL, 100003LL })
            test::ok (runChunked (prog, sr, nc, chunk) == whole,
                      "chunk " + std::to_string (chunk) + " gives a bit-identical surface");
    }

    // --- Same input, same answer, twice. A stray static or an uninitialised accumulator dies here. ---
    test::group ("determinism: the same program measured twice is bit-identical");
    {
        const Planar prog = makeProgram (44100.0, 1, 4.0);
        test::ok (runChunked (prog, 44100.0, 1, 0) == runChunked (prog, 44100.0, 1, 0), "run twice, same bits");
    }

    // --- prepare() is the guard rail an ABI leans on: the wasm shim will hand it whatever JS produced. ---
    test::group ("prepare() rejects every unusable configuration");
    {
        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        struct Case { double sr; int ch; bool want; const char* what; };
        const Case cases[] {
            { 48000.0,  2, true,  "48 kHz stereo" },
            { 44100.0,  1, true,  "44.1 kHz mono" },
            { 192000.0, 2, true,  "192 kHz stereo" },
            { 48000.0,  core::kMaxChannels, true, "the maximum channel count" },
            { 0.0,      2, false, "sampleRate 0" },
            { -48000.0, 2, false, "negative sampleRate" },
            { nan,      2, false, "NaN sampleRate" },
            { inf,      2, false, "+inf sampleRate (passes a naive `sr > 0` check)" },
            { 48000.0,  0, false, "channels 0" },
            { 48000.0, -1, false, "negative channels" },
            { 48000.0,  core::kMaxChannels + 1, false, "channels past kMaxChannels" },
            // "Positive and finite" is not enough. LoudnessMeter sizes its hop with lround(0.01*fs) into an
            // int and its block store with ceil(maxDurationSec*fs); an absurd-but-finite rate reaches an
            // out-of-range lround (undefined behaviour) or demands an impossible allocation. A page can hand
            // us Number.MIN_VALUE as easily as 48000.
            { fcore::Probe::kMinSampleRate,       2, true,  "the lowest accepted rate" },
            { fcore::Probe::kMaxSampleRate,       2, true,  "the highest accepted rate" },
            { fcore::Probe::kMinSampleRate - 1.0, 2, false, "just below the lowest" },
            { fcore::Probe::kMaxSampleRate + 1.0, 2, false, "just above the highest" },
            { 1e300,    2, false, "an absurd but finite rate (lround would be out of range)" },
            { 5e-324,   2, false, "the smallest positive subnormal double" },
        };
        for (const auto& k : cases)
        {
            fcore::Probe p;
            test::ok (p.prepare (k.sr, k.ch) == k.want, std::string ("prepare: ") + k.what);
            test::ok (p.prepared() == k.want, std::string ("prepared() agrees: ") + k.what);
        }
    }

    test::group ("an unprepared Probe is inert, not a crash");
    {
        fcore::Probe p;
        std::vector<float> buf (256, 0.5f);
        const float* io[1] { buf.data() };
        p.process (io, 1, 256);                                   // must be a silent no-op
        test::ok (p.gatingBlockCount() == 0, "no blocks");
        test::ok (p.truePeakLinear() == 0.0, "no true peak");
    }

    // --- Degenerate lengths: the wasm shim will be handed frames=0 and frames=1 by a fuzzing page. ---
    test::group ("degenerate lengths and silence");
    {
        fcore::Probe p; p.prepare (48000.0, 2);
        std::vector<float> z (1, 0.0f);
        const float* io[2] { z.data(), z.data() };
        p.process (io, 2, 0);
        test::ok (p.gatingBlockCount() == 0 && p.truePeakLinear() == 0.0, "0 frames measures nothing");
        p.process (io, 2, 1);
        test::ok (p.gatingBlockCount() == 0, "1 frame is not yet a gating block");
        test::ok (p.truePeakLinear() == 0.0, "and silence has no true peak");
        test::ok (std::isfinite (p.integratedLufs()), "integratedLufs() stays finite on silence");
        test::approx (p.integratedLufs(), -120.0, 1e-12, "silence reads the −120 sentinel");
        test::ok (p.truePeakDb() == 20.0 * std::log10 (1e-9), "true-peak dB clamps rather than returning −inf");
    }

    test::group ("negative lengths are ignored rather than read out of bounds");
    {
        fcore::Probe p; p.prepare (48000.0, 1);
        std::vector<float> z (16, 0.5f);
        const float* io[1] { z.data() };
        p.process (io, 1, -1);
        p.process (io, 1, -100000);
        test::ok (p.gatingBlockCount() == 0 && p.truePeakLinear() == 0.0, "a negative frame count measures nothing");
    }

    test::group ("a non-positive channel count does not silently lengthen the program");
    {
        // Passing 0 channels is not "measure nothing" by default: LoudnessMeter would still advance its hop
        // clock and record SILENT gating blocks, appending material that was never submitted.
        fcore::Probe p; p.prepare (48000.0, 2);
        std::vector<float> s (48000, 0.5f);
        const float* io[2] { s.data(), s.data() };
        p.process (io, 0, (long long) s.size());
        p.process (io, -3, (long long) s.size());
        test::ok (p.gatingBlockCount() == 0, "no gating blocks were invented");
        test::ok (p.truePeakLinear() == 0.0 && p.samplePeakLinear() == 0.0, "and nothing was measured");
        p.process (io, 2, (long long) s.size());
        test::ok (p.gatingBlockCount() > 0, "a real call afterwards still works");
    }

    test::group ("maxDurationSec is validated too");
    {
        fcore::Probe p;
        test::ok (! p.prepare (48000.0, 2, 0.0), "zero capacity");
        test::ok (! p.prepare (48000.0, 2, -1.0), "negative capacity");
        test::ok (! p.prepare (48000.0, 2, std::numeric_limits<double>::quiet_NaN()), "NaN capacity");
        test::ok (! p.prepare (48000.0, 2, std::numeric_limits<double>::infinity()), "infinite capacity");
        test::ok (p.prepare (48000.0, 2, 10.0), "a sane one is accepted");
    }

    // --- Full scale and beyond: true peak must exceed sample peak on a signal built to have inter-sample
    //     overs, or the oversampler is not actually doing anything. ---
    test::group ("true peak exceeds sample peak on an inter-sample over");
    {
        const double sr = 48000.0;
        // ±1 alternating at Nyquist/2 has its reconstruction peaks BETWEEN the samples.
        std::vector<float> s ((std::size_t) 48000);
        for (std::size_t i = 0; i < s.size(); ++i) s[i] = (i % 4 < 2) ? 0.9f : -0.9f;
        fcore::Probe p; p.prepare (sr, 1);
        const float* io[1] { s.data() };
        p.process (io, 1, (long long) s.size());
        test::ok (p.truePeakLinear() > 0.9, "the reconstructed peak sits above the 0.9 sample peak");
        test::ok (p.truePeakLinear() < 2.0, "and not absurdly above it");
        test::ok (std::isfinite (p.truePeakDb()), "the dB form is finite");
    }

    // --- finish() drains the FIR. Without it a peak in the final samples is not merely imprecise, it is
    //     missing: the oversampler is causal with a 63.5-oversampled-sample group delay, so the last ~16
    //     baseband samples never leave the filter while input is still arriving. ---
    test::group ("finish(): a transient at the very end of the buffer is measured, not lost");
    {
        // A file that stops dead on a loud passage. Its true peak is ABOVE full scale; without draining, the
        // tool reports about −36 dBTP.
        std::vector<float> s (48000, 0.0f);
        for (std::size_t i = s.size() - 10; i < s.size(); ++i) s[i] = 0.95f;
        const float* io[1] { s.data() };

        fcore::Probe undrained; undrained.prepare (48000.0, 1);
        undrained.process (io, 1, (long long) s.size());
        const double withoutFinish = undrained.truePeakLinear();

        fcore::Probe drained; drained.prepare (48000.0, 1);
        drained.process (io, 1, (long long) s.size());
        drained.finish();
        const double withFinish = drained.truePeakLinear();

        // The reference: the same signal with room to drain naturally.
        std::vector<float> padded = s; padded.resize (s.size() + 256, 0.0f);
        fcore::Probe ref; ref.prepare (48000.0, 1);
        const float* io2[1] { padded.data() };
        ref.process (io2, 1, (long long) padded.size());
        ref.finish();

        test::ok (withFinish > 1.0, "the true peak of an abrupt 0.95 ending is ABOVE full scale");
        test::approx (withFinish, ref.truePeakLinear(), 1e-12, "finish() reaches what natural drain reaches");

        // The two mitigations do different amounts of work, and it is worth pinning both. The oversampler
        // alone sees almost nothing here — 0.0152, i.e. −36 dBTP for a signal that is over full scale. The
        // sample-peak floor drags that up to the sample peak on its own, which is most of the rescue; finish()
        // then supplies the inter-sample part the floor cannot know about.
        test::ok (sameBits (withoutFinish, (double) (float) 0.95f),
                  "undrained, the reading falls back to the sample-peak floor — the FIR contributed nothing");
        test::ok (withFinish > withoutFinish * 1.1,
                  "and the drain still adds about a dB of inter-sample peak the floor cannot see");

        drained.finish();
        test::approx (drained.truePeakLinear(), withFinish, 0.0, "finish() is idempotent");
        drained.process (io, 1, 16);
        test::approx (drained.truePeakLinear(), withFinish, 0.0, "and process() after finish() is a no-op");
    }

    test::group ("the sample peak is a hard floor under the true peak");
    {
        for (const double amp : { 0.1, 0.5, 0.999, 1.0 })
        {
            std::vector<float> s (4096, 0.0f);
            s[2048] = (float) amp;                              // a lone impulse, mid-buffer
            fcore::Probe p; p.prepare (48000.0, 1);
            const float* io[1] { s.data() };
            p.process (io, 1, (long long) s.size());
            p.finish();
            test::ok (p.samplePeakLinear() == (double) (float) amp, "the sample peak is exactly the loudest sample");
            test::ok (p.truePeakLinear() >= p.samplePeakLinear(),
                      "true peak >= sample peak at amplitude " + std::to_string (amp));
        }
    }

    // --- The dB form and the linear form must be the same number. They were not: truePeakDb() was computed
    //     from the oversampler maximum while truePeakLinear() applied the sample-peak floor, so the two
    //     disagreed on exactly the signals the floor exists for. Caught in review; pinned here. ---
    test::group ("truePeakDb() is the dB of truePeakLinear(), including on the signals where the floor bites");
    {
        for (const bool endOnTransient : { false, true })
        {
            std::vector<float> s (4096, 0.0f);
            s[endOnTransient ? s.size() - 1 : 2048] = 1.0f;      // a lone impulse: the floor is the answer
            fcore::Probe p; p.prepare (48000.0, 1);
            const float* io[1] { s.data() };
            p.process (io, 1, (long long) s.size());
            p.finish();
            test::approx (p.truePeakDb(), 20.0 * std::log10 (p.truePeakLinear()), 1e-12,
                          endOnTransient ? "impulse at the very end" : "impulse mid-buffer");
        }
    }

    // --- The drain length must actually be enough to empty the ring, not merely "enough in practice". ---
    test::group ("finish() drains the whole ring: one more zero would add nothing");
    {
        std::vector<float> s (2048, 0.0f);
        s[s.size() - 1] = 1.0f;
        const float* io[1] { s.data() };

        fcore::Probe drained; drained.prepare (48000.0, 1);
        drained.process (io, 1, (long long) s.size());
        drained.finish();

        // The same signal followed by MORE silence than finish() pushes must reach the same answer: if 32
        // zeros left anything in the filter, the longer tail would find it.
        std::vector<float> padded = s; padded.resize (s.size() + 4 * fcore::Probe::kOsTapsPerPhase, 0.0f);
        fcore::Probe longer; longer.prepare (48000.0, 1);
        const float* io2[1] { padded.data() };
        longer.process (io2, 1, (long long) padded.size());
        longer.finish();

        test::approx (drained.truePeakLinear(), longer.truePeakLinear(), 0.0,
                      "32 zeros flush the 32-sample ring exactly — four times as many find nothing more");
    }

    // --- The wasm shim keeps ONE Probe and re-prepares it per call, so a stale byte between runs would be a
    //     silent, cross-file corruption. Prove re-prepare is a real reset at the Probe level. ---
    test::group ("re-prepare is a clean reset: A → B → A reproduces A exactly");
    {
        const Planar a = makeProgram (48000.0, 2, 3.0);
        const Planar b = makeProgram (44100.0, 1, 2.0);

        fcore::Probe p;
        auto runOn = [&p] (const Planar& prog, double sr, int nc)
        {
            p.prepare (sr, nc);
            const auto base = ptrs (prog);
            p.process (base.data(), nc, (long long) prog[0].size());
            p.finish();
            return surfaceOf (p);
        };

        const Surface first  = runOn (a, 48000.0, 2);
        (void)                 runOn (b, 44100.0, 1);
        const Surface again  = runOn (a, 48000.0, 2);
        test::ok (first == again, "a different rate and channel count in between left nothing behind");
    }

    // --- Non-finite SAMPLES. Documenting what actually happens, because the honest answer is unpleasant and
    //     a future reader must not mistake silence for safety. ---
    test::group ("non-finite samples: the documented (unpleasant) behaviour");
    {
        fcore::Probe p; p.prepare (48000.0, 1);
        std::vector<float> good (8192, 0.5f);
        const float* io[1] { good.data() };
        p.process (io, 1, 8192);
        const double tpBefore = p.truePeakLinear();
        test::ok (tpBefore > 0.4, "a real reading first");

        const double lufsBefore = p.integratedLufs();

        std::vector<float> bad (8192, 0.5f);
        bad[100] = std::numeric_limits<float>::quiet_NaN();
        const float* io2[1] { bad.data() };
        p.process (io2, 1, 8192);

        // The two paths behave OPPOSITELY, and the difference is the whole point:
        //
        //  * the true peak RECOVERS. The oversampler's history is a 32-sample ring, so the NaN shifts out of
        //    it; `std::max(x, NaN)` returns `x` meanwhile, so the maximum is merely blind for those samples
        //    and then works again. Feed it something LOUDER afterwards and it duly rises — which is exactly
        //    why the naive test (feeding the same level after the NaN) cannot tell recovery from a freeze.
        //  * the loudness does NOT recover. K-weighting is an IIR: once its state is NaN it is NaN forever,
        //    so every later block energy is NaN. And the reading never goes NaN either, which would at least
        //    be visible — `NaN > absT` is false, so the absolute gate silently DROPS every poisoned block and
        //    the meter keeps reporting a healthy number computed over the fraction of the program that
        //    predates the NaN. droppedBlocks() does not mention it; it counts capacity overflow only.
        //
        // Recorded for the fix session, pinned here so it cannot change unnoticed. NOT fixed in the spike.
        std::vector<float> louder (8192, 0.95f);
        const float* io3[1] { louder.data() };
        p.process (io3, 1, 8192);
        test::ok (p.truePeakLinear() > tpBefore + 0.1,
                  "the true peak RECOVERS once the NaN shifts out of the 32-sample FIR ring");

        test::ok (! std::isnan (p.integratedLufs()), "the loudness reading does not go NaN, which would be visible");
        test::ok (sameBits (p.integratedLufs(), lufsBefore),
                  "it is frozen instead: every post-NaN block is silently dropped by the absolute gate");
        int nanBlocks = 0;
        for (int j = 0; j < p.gatingBlockCount(); ++j)
            if (std::isnan (p.gatingBlockEnergies()[j])) ++nanBlocks;
        test::ok (nanBlocks > 0, "the dropped blocks are there in the vector, as NaN — invisible to the reading");
        test::ok (p.droppedBlocks() == 0, "and droppedBlocks() does NOT report them (it counts capacity only)");
    }

    // --- +inf is not NaN and does not behave like it. Worth its own group, because the loudness answer it
    //     produces is the most misleading one this code can emit: it reads as SILENCE. ---
    test::group ("+inf samples: the failure mode that reads as silence");
    {
        const double inf = std::numeric_limits<float>::infinity();
        {
            // +inf before the first gating block closes: every block that ever forms is NaN (the K-weighting
            // IIR is poisoned from that sample on), all are dropped by the absolute gate, and the meter has
            // nothing left to average — so it answers −120, the sentinel that means "silence".
            fcore::Probe p; p.prepare (48000.0, 1);
            std::vector<float> s (48000, 0.5f);
            s[1000] = (float) inf;
            const float* io[1] { s.data() };
            p.process (io, 1, (long long) s.size());
            p.finish();
            test::approx (p.integratedLufs(), -120.0, 1e-12,
                          "a loud signal containing one +inf reads as SILENCE, not as an error");
            test::ok (std::isinf (p.truePeakDb()) || p.truePeakDb() > 100.0,
                      "while the true peak reports the infinity honestly");
        }
        {
            // +inf after real program has been measured: the earlier blocks survive, the later ones are
            // dropped, and the reading freezes at whatever the pre-inf material said.
            fcore::Probe p; p.prepare (48000.0, 1);
            std::vector<float> good (48000, 0.5f);
            const float* g[1] { good.data() };
            p.process (g, 1, (long long) good.size());
            const double before = p.integratedLufs();
            std::vector<float> bad (48000, 0.5f);
            bad[100] = (float) inf;
            const float* b[1] { bad.data() };
            p.process (b, 1, (long long) bad.size());
            p.process (g, 1, (long long) good.size());
            test::ok (std::isfinite (p.integratedLufs()), "the reading stays finite and plausible");
            test::ok (sameBits (p.integratedLufs(), before), "and frozen at the pre-inf value");
        }
    }

    // --- The gating-block cadence, at every rate the family plausibly meets. Also the evidence table for the
    //     escalated cross-tier determinism question, which is about rates. ---
    test::group ("gating-block cadence holds at every rate");
    {
        for (const double sr : { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
        {
            const Planar prog = makeProgram (sr, 2, 3.0);
            fcore::Probe p; p.prepare (sr, 2);
            const auto base = ptrs (prog);
            p.process (base.data(), 2, (long long) prog[0].size());
            test::ok (p.gatingBlockCount() == 27,
                      "3 s at " + std::to_string ((int) sr) + " Hz → 30 hops less the 3 before the first block");
            test::ok (p.droppedBlocks() == 0, "and nothing was dropped");
            test::ok (std::isfinite (p.integratedLufs()) && p.integratedLufs() < 0.0, "the reading is sane");
        }
    }

    // --- The true-peak CONFIG is part of the parity contract. The core holds a SECOND true-peak filter with a
    //     different length, cutoff and window; swapping to it would fail the spike by ~3e-3 dB while looking
    //     like a wasm bug. Pin both the constants and the fact that the two really do disagree. ---
    test::group ("true-peak config is pinned, and the core's other true-peak filter really does differ");
    {
        test::ok (fcore::Probe::kOsFactor == 4, "4× oversampling, as the reference tool has always used");
        test::ok (fcore::Probe::kOsTapsPerPhase == 32, "32 taps/phase — a 128-tap prototype");

        // A 12 kHz sine phased to miss the sample crests: the classic inter-sample-peak fixture, and the one
        // place the two filter DESIGNS actually diverge. On a broadband impulse they now agree exactly, since
        // both land on the sample-peak floor; on ordinary music they agree to 0.0009–0.0028 dB.
        const double sr = 48000.0;
        const long long n = (long long) (sr * 2);
        Planar prog (2, std::vector<float> ((std::size_t) n));
        for (long long i = 0; i < n; ++i)
        {
            const float v = (float) (0.9 * std::sin (2.0 * core::kPi * 12000.0 * (double) i / sr + 0.7));
            prog[0][(std::size_t) i] = prog[1][(std::size_t) i] = v;
        }

        fcore::Probe p; p.prepare (sr, 2);
        const auto base = ptrs (prog);
        p.process (base.data(), 2, n);
        p.finish();

        analysis::TruePeakMeter tpm; tpm.prepare (sr, fcore::Probe::kChunk, 2);
        tpm.process (base.data(), 2, (int) n);

        const double a = p.truePeakDb(), b = tpm.truePeakDb();
        test::ok (std::isfinite (a) && std::isfinite (b), "both filters produced a reading");
        test::ok (a > 20.0 * std::log10 (p.samplePeakLinear()) + 1.0,
                  "the inter-sample peak really is well above the sample peak here");
        test::ok (! sameBits (a, b),
                  "the two designs are NOT bit-identical — a shim reaching for TruePeakMeter fails parity");
        test::ok (std::fabs (a - b) > 0.01 && std::fabs (a - b) < 0.5,
                  "they differ by tens of a dB on high-frequency content, not by a rounding step");
        // NB: both under-read ffmpeg's ebur128 here by ~0.5 dB (it reports −0.3, we report −0.81 / −0.88), so
        // ffmpeg does NOT arbitrate between them — that is a separate question about oversampling factor, not
        // about which of these two prototypes is right. Recorded, not resolved.
    }

    // --- The RT claim the core lives by, on the shared body too. ---
    test::group ("process() does not allocate");
    {
        const Planar prog = makeProgram (48000.0, 2, 2.0);
        fcore::Probe p; p.prepare (48000.0, 2);
        const auto base = ptrs (prog);
        p.process (base.data(), 2, 4800);                       // warm: first hop, first block
        const long before = g_allocs.load();
        p.process (base.data(), 2, (long long) prog[0].size()); // every internal chunk step, hops and blocks
        test::okNoAlloc (g_allocs.load() == before, "a full pass through process() allocated nothing");
        test::ok (std::isfinite (p.integratedLufs()), "and it still reads");
    }

    // --- Capacity: past the prepared duration blocks are COUNTED, not silently lost, and the CLI warns. ---
    test::group ("capacity overflow is reported, not silent");
    {
        fcore::Probe p;
        test::ok (p.prepare (48000.0, 1, 1.0), "prepare for 1 s of capacity");
        const Planar prog = makeProgram (48000.0, 1, 3.0);
        const auto base = ptrs (prog);
        p.process (base.data(), 1, (long long) prog[0].size());
        test::ok (p.droppedBlocks() > 0, "3 s into 1 s of capacity drops blocks");
        test::ok (std::isfinite (p.integratedLufs()), "and the reading still stands over what was kept");
    }

    return test::report();
}
