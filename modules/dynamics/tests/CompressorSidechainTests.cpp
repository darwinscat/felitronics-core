// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the Compressor's EXTERNAL DETECTOR (sidechain / key) input, and for the
// hostile state sequences that surface around an argument which is present on some blocks and absent on
// others. That shape is exactly what the last two defects in this repo lived in, and both times the test
// that was supposed to pin the behaviour asserted the SHAPE OF THE CALL instead of a property under an
// unfriendly ORDER of states — so nearly everything here is a sequence, not a single call.
//
// The organising idea of the suite is a REFERENCE NULL. `Oracle` below rebuilds the whole compressor
// out of the module's PUBLIC primitives — LinkedDetector, GainComputer, GainReductionFollower,
// core::DelayLine — configured by handing it the very same CompressorParams object. If that assembly
// nulls the real thing bit for bit, then two things are true at once: the behaviour is pinned against
// any future drift, and "an offline pass can run the same detector the compressor runs" is a fact about
// the code rather than a claim in a comment. The next task (P3) stands on that.

#include <felitronics_test.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/LinkedDetector.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

// global allocation counter (no-alloc-in-process proof)
static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

// Compile-time contract detectors for the call shapes the API deliberately does and does not offer.
namespace detail
{
    template <class C, class = void> struct HasFourArgProcess : std::false_type {};
    template <class C> struct HasFourArgProcess<C, std::void_t<decltype (std::declval<C&>().process (
        std::declval<float* const*>(), 0, 0, std::declval<const float* const*>()))>> : std::true_type {};

    template <class C, class = void> struct HasFiveArgProcess : std::false_type {};
    template <class C> struct HasFiveArgProcess<C, std::void_t<decltype (std::declval<C&>().process (
        std::declval<float* const*>(), 0, 0, std::declval<const float* const*>(), 0))>> : std::true_type {};
}

static constexpr double kFs = 48000.0;
static const float kInf = std::numeric_limits<float>::infinity();
static const float kNan = std::numeric_limits<float>::quiet_NaN();

//==============================================================================
// The ORACLE: the documented topology, reassembled from the public primitives only.
struct Oracle
{
    dynamics::LinkedDetector        det;
    dynamics::GainComputer          gc;
    dynamics::GainReductionFollower gr;
    std::vector<core::DelayLine>    delays;
    double makeupDb = 0.0;
    int    look = 0, lastNc = 0;

    void prepare (double fs, int maxChannels, double maxLookaheadMs, const dynamics::CompressorParams& p)
    {
        det.prepare (fs);
        gr.prepare (fs);
        const int maxLook = (int) std::ceil (maxLookaheadMs * 0.001 * fs);
        delays.assign ((std::size_t) maxChannels, core::DelayLine {});
        for (auto& d : delays) d.prepare (maxLook);
        setParams (fs, maxLook, p);
        for (auto& d : delays) d.reset();
        det.reset(); gr.reset(); lastNc = 0;
    }

    void setParams (double fs, int maxLook, const dynamics::CompressorParams& p)
    {
        gc.setMode (p.mode); gc.setThresholdDb (p.thresholdDb); gc.setRatio (p.ratio);
        gc.setKneeDb (p.kneeDb); gc.setRangeDb (p.rangeDb);
        gr.setTimes (p.attackMs, p.releaseMs);
        det.setParams (p);                                   // the SAME object, sliced to its detector half
        const double ms = std::isfinite (p.lookaheadMs) ? std::clamp (p.lookaheadMs, 0.0, 250.0) : 0.0;
        const int newLook = std::clamp ((int) std::lround (ms * 0.001 * fs), 0, maxLook);
        const bool changed = (newLook != look);
        look = newLook;
        for (auto& d : delays) d.setDelay (look);
        if (changed) for (auto& d : delays) d.reset();
        makeupDb = (std::isfinite (p.makeupDb) ? p.makeupDb : 0.0) + (p.autoMakeup ? dynamics::autoMakeupDb (gc) : 0.0);
    }

    void process (float* const* io, int nc, int n, const float* const* key = nullptr, int keyNc = 0)
    {
        if (lastNc > 0 && nc > lastNc) for (int c = lastNc; c < nc; ++c) delays[(std::size_t) c].reset();
        lastNc = nc;
        const float* const* dc = key; int dn = keyNc;
        if (key == nullptr || keyNc <= 0) { dc = io; dn = nc; }
        for (int i = 0; i < n; ++i)
        {
            const float  level    = det.process (dc, dn, i);
            const float  targetDb = (float) gc.deltaDb (core::gainToDb (level));
            const float  grDb     = gr.process (targetDb);
            const float  gain     = (float) core::dbToGain ((double) grDb + makeupDb);
            for (int c = 0; c < nc; ++c) { const float x = io[c][i]; io[c][i] = delays[(std::size_t) c].process (x) * gain; }
        }
        det.flushDenormals(); gr.flushDenormals();
    }
};

//==============================================================================
static std::vector<float> tone (int n, double f, double amp, double phase = 0.0)
{
    std::vector<float> v ((std::size_t) n);
    for (int i = 0; i < n; ++i) v[(std::size_t) i] = (float) (amp * std::sin (2.0 * core::kPi * f * i / kFs + phase));
    return v;
}

static std::size_t countDiff (const std::vector<float>& a, const std::vector<float>& b)
{
    std::size_t d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) if (std::memcmp (&a[i], &b[i], sizeof (float)) != 0) ++d;
    return d;
}

static bool allFinite (const std::vector<float>& v)
{
    for (float x : v) if (! std::isfinite (x)) return false;
    return true;
}

static double rmsOf (const std::vector<float>& v, int from, int to)
{
    double s = 0.0; int c = 0;
    for (int i = from; i < to && i < (int) v.size(); ++i) { s += (double) v[(std::size_t) i] * v[(std::size_t) i]; ++c; }
    return c ? std::sqrt (s / c) : 0.0;
}

// A wide-band, wide-dynamics probe signal: an amplitude-swept low tone, a bright layer, periodic bursts.
static std::vector<float> probeSignal (int n, double detune)
{
    std::vector<float> v ((std::size_t) n);
    for (int i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        const double env = 0.05 + 0.9 * std::fabs (std::sin (2.0 * core::kPi * 1.3 * t));
        double x = env * (0.6 * std::sin (2.0 * core::kPi * (55.0 + detune) * t) + 0.25 * std::sin (2.0 * core::kPi * 2700.0 * t));
        if ((i % 977) < 24) x += 0.8 * std::sin (2.0 * core::kPi * 3000.0 * t);
        v[(std::size_t) i] = (float) x;
    }
    return v;
}

static dynamics::CompressorParams baseParams()
{
    dynamics::CompressorParams p;
    p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 5.0; p.link = dynamics::LinkMode::Max;
    p.thresholdDb = -24.0; p.ratio = 4.0; p.kneeDb = 0.0; p.rangeDb = 40.0;
    p.attackMs = 3.0; p.releaseMs = 60.0; p.makeupDb = 0.0; p.lookaheadMs = 0.0;
    return p;
}

//==============================================================================
int main()
{
    std::printf ("felitronics::dynamics Compressor external-key (sidechain) tests\n");

    //==========================================================================
    // A. THE CONTRACT
    //==========================================================================

    // Acceptance 1, in its exact wording: feeding `io` as the key must be TODAY'S output, bit for bit.
    // It is not a coincidence that it holds — at every sample index the detector reads before anything
    // writes — so this pins the aliasing order as much as the routing.
    test::group ("key == io is the self-keyed path, bit for bit");
    {
        std::size_t worst = 0; int cases = 0;
        for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
          for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
            for (double look : { 0.0, 2.5 })
              for (int nch = 1; nch <= 2; ++nch)
                for (int blk : { 1, 37, 512 })
                {
                    const int n = 4096;
                    auto s0 = probeSignal (n, 0.0), s1 = probeSignal (n, 11.0);
                    auto a0 = s0, a1 = s1, b0 = s0, b1 = s1;
                    dynamics::Compressor ca, cb;
                    test::ok (ca.prepare (kFs, blk, nch, 20.0) && cb.prepare (kFs, blk, nch, 20.0), "prepare");
                    auto p = baseParams(); p.detector = d; p.link = l; p.lookaheadMs = look;
                    ca.setParams (p); cb.setParams (p);
                    for (int off = 0; off < n; off += blk)
                    {
                        const int m = std::min (blk, n - off);
                        float* ia[2] { a0.data() + off, a1.data() + off };
                        float* ib[2] { b0.data() + off, b1.data() + off };
                        const float* kb[2] { b0.data() + off, b1.data() + off };   // the io buffer itself
                        ca.process (ia, nch, m);
                        cb.process (ib, nch, m, kb, nch);
                        }
                    worst += countDiff (a0, b0) + (nch == 2 ? countDiff (a1, b1) : 0);
                    ++cases;
                }
        test::ok (worst == 0, "self-keyed == io-as-key over " + std::to_string (cases) + " configurations, bit for bit");
    }

    // Acceptance 4. The detector path is assemblable from OUTSIDE out of the same types with the same
    // params — demonstrated in the strongest available form: the whole compressor rebuilt from public
    // primitives nulls the real one exactly, keyed and self-keyed alike.
    test::group ("an external assembly from the public primitives nulls the compressor exactly");
    {
        std::size_t diff = 0; int cases = 0;
        for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
          for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
            for (auto m : { dynamics::Mode::DownCompress, dynamics::Mode::UpCompress, dynamics::Mode::DownExpand })
              for (double look : { 0.0, 1.0, 6.25 })
                for (double win : { 0.5, 12.0 })
                  for (bool keyed : { false, true })
                  {
                      const int n = 3000, nch = 2;
                      auto s0 = probeSignal (n, 0.0), s1 = probeSignal (n, 7.0);
                      auto key0 = tone (n, 700.0, 0.35), key1 = tone (n, 90.0, 0.6);
                      auto a0 = s0, a1 = s1, b0 = s0, b1 = s1;

                      auto p = baseParams();
                      p.detector = d; p.link = l; p.mode = m; p.lookaheadMs = look; p.rmsWindowMs = win;
                      p.makeupDb = 1.25; p.autoMakeup = true; p.kneeDb = 4.0;

                      dynamics::Compressor c; test::ok (c.prepare (kFs, 512, nch, 20.0), "prepare");
                      c.setParams (p);
                      Oracle o; o.prepare (kFs, nch, 20.0, p);

                      float* ia[2] { a0.data(), a1.data() };
                      float* ib[2] { b0.data(), b1.data() };
                      const float* k[2] { key0.data(), key1.data() };
                      if (keyed) { c.process (ia, nch, n, k, nch); o.process (ib, nch, n, k, nch); }
                      else       { c.process (ia, nch, n);          o.process (ib, nch, n); }

                      diff += countDiff (a0, b0) + countDiff (a1, b1);
                      ++cases;
                  }
        test::ok (diff == 0, "oracle nulls the compressor over " + std::to_string (cases) + " configurations, bit for bit");
    }

    // The detector must read the KEY and nothing else, and the gain must land on `io` and nothing else.
    // Both directions, because an implementation that filtered `io` instead would pass only one of them.
    test::group ("isolation: the key drives, the programme is driven");
    {
        const int n = 4096;
        auto loud = tone (n, 400.0, 0.7);
        std::vector<float> silence ((std::size_t) n, 0.0f);
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 8.0; p.attackMs = 1.0;

        {   // silent programme, loud key: the programme stays EXACTLY silent, but gain reduction happens
            auto a = silence; float* io[1] { a.data() }; const float* k[1] { loud.data() };
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            c.process (io, 1, n, k, 1);
            bool zero = true; for (float v : a) zero &= (v == 0.0f);
            test::ok (zero, "silent programme + loud key: output is exactly silent (no key leakage into io)");
            test::ok (c.gainReductionDb() < -10.0, "silent programme + loud key: the key still drives the gain reduction");
        }
        {   // loud programme, silent key: nothing happens AT ALL — gain is exactly 1.0
            auto a = loud; float* io[1] { a.data() }; const float* k[1] { silence.data() };
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            c.process (io, 1, n, k, 1);
            test::ok (countDiff (a, loud) == 0, "loud programme + silent key: output is the input, bit for bit");
        }
    }

    // A mono key on a stereo programme is a first-class case, not a workaround, and it must agree with
    // the same key duplicated into two channels. MeanPower is the interesting one: it divides by the
    // KEY's channel count, so mean(k^2, k^2) == k^2 exactly — the doubling and the halving cancel.
    test::group ("a mono key on a stereo programme == the same key duplicated");
    {
        for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
        {
            const int n = 3000;
            auto s0 = probeSignal (n, 0.0), s1 = probeSignal (n, 5.0);
            auto key = tone (n, 120.0, 0.55);
            auto a0 = s0, a1 = s1, b0 = s0, b1 = s1;
            auto p = baseParams(); p.link = l;

            dynamics::Compressor ca, cb;
            test::ok (ca.prepare (kFs, 512, 2) && cb.prepare (kFs, 512, 2), "prepare");
            ca.setParams (p); cb.setParams (p);
            float* ia[2] { a0.data(), a1.data() }; const float* k1[1] { key.data() };
            float* ib[2] { b0.data(), b1.data() }; const float* k2[2] { key.data(), key.data() };
            ca.process (ia, 2, n, k1, 1);
            cb.process (ib, 2, n, k2, 2);
            test::ok (countDiff (a0, b0) + countDiff (a1, b1) == 0,
                      std::string ("mono key == duplicated stereo key (") + (l == dynamics::LinkMode::Max ? "Max" : "MeanPower") + ")");
        }
    }

    // The key must be READ-ONLY. It is `const float* const*` in the signature, but a product may well be
    // metering or displaying the same buffer, so pin that nothing writes through it.
    test::group ("the key buffer is not written");
    {
        const int n = 2048;
        auto key = probeSignal (n, 3.0);
        const auto keyCopy = key;
        auto a = probeSignal (n, 0.0);
        float* io[1] { a.data() }; const float* k[1] { key.data() };
        auto p = baseParams(); p.thresholdDb = -40.0; p.ratio = 10.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1, 10.0), "prepare");
        p.lookaheadMs = 3.0; c.setParams (p);
        c.process (io, 1, n, k, 1);
        test::ok (countDiff (key, keyCopy) == 0, "the compressor did not modify the key buffer");
    }

    // Block-length invariance with a key: the per-sample loop must carry the key offset correctly.
    // A key pointer that is not advanced per block passes a single-block test and fails this one.
    test::group ("block-length invariance with a key");
    {
        const int n = 5000;
        auto s = probeSignal (n, 0.0);
        auto key = probeSignal (n, 400.0);
        std::vector<float> ref;
        for (int blk : { 1, 3, 64, 512, 5000 })
        {
            auto a = s; float* io[1] { a.data() };
            dynamics::Compressor c; test::ok (c.prepare (kFs, blk, 1, 10.0), "prepare");
            auto p = baseParams(); p.lookaheadMs = 2.0; c.setParams (p);
            for (int off = 0; off < n; off += blk)
            {
                const int m = std::min (blk, n - off);
                float* sub[1] { a.data() + off };
                const float* k[1] { key.data() + off };
                c.process (sub, 1, m, k, 1);
            }
            if (ref.empty()) ref = a;
            else test::ok (countDiff (a, ref) == 0, "block " + std::to_string (blk) + " == block 1, bit for bit");
        }
    }

    // A block LARGER than the maxBlock passed to prepare() must behave — the compressor sizes no scratch
    // from it, so an offline caller handing it a whole file is a normal thing to do.
    test::group ("a block larger than maxBlock is processed, not dropped");
    {
        const int n = 40000;
        auto s = probeSignal (n, 0.0); auto key = probeSignal (n, 0.0);
        auto a = s, b = s;
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 6.0;
        dynamics::Compressor c1, c2;
        test::ok (c1.prepare (kFs, 64, 1) && c2.prepare (kFs, n, 1), "prepare");
        c1.setParams (p); c2.setParams (p);
        float* io1[1] { a.data() }; float* io2[1] { b.data() };
        const float* k[1] { key.data() };
        c1.process (io1, 1, n, k, 1);
        c2.process (io2, 1, n, k, 1);
        test::ok (countDiff (a, b) == 0, "maxBlock 64 fed 40000 samples == maxBlock 40000, bit for bit");
        test::ok (rmsOf (a, 20000, n) < rmsOf (s, 20000, n) * 0.9, "and it actually compressed (not a silent bypass)");
    }

    // Time alignment. An impulse in the key at t must move the gain at t, and the programme must come
    // out `lookahead` samples later — a one-sample slip in either direction shows up here and nowhere
    // else. The oracle test above cannot catch it: the oracle would slip identically.
    test::group ("key/programme alignment under lookahead");
    {
        for (double lookMs : { 0.0, 1.0 })
        {
            const int n = 2048, t = 500;
            std::vector<float> a ((std::size_t) n, 0.0f), key ((std::size_t) n, 0.0f);
            a[(std::size_t) t] = 0.9f; key[(std::size_t) t] = 0.9f;
            const auto in = a;
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1, 10.0), "prepare");
            auto p = baseParams();
            p.detector = dynamics::Detector::Peak; p.thresholdDb = -12.0; p.ratio = 4.0;
            // A 10 s release makes the decay across the lookahead negligible (< 1e-4 dB), so what this
            // measures is WHERE the sample lands, not how far the follower has crept by the time it does.
            p.attackMs = 0.0; p.releaseMs = 10000.0; p.lookaheadMs = lookMs;
            c.setParams (p);
            const int look = c.latencySamples();
            float* io[1] { a.data() }; const float* k[1] { key.data() };
            c.process (io, 1, n, k, 1);
            // static target for 0.9 at -12 dBFS / 4:1  ->  -(20log10(0.9)+12)*0.75 = -8.31 dB
            const double lvlDb = core::gainToDb ((double) in[(std::size_t) t]);
            const double want  = (double) in[(std::size_t) t] * core::dbToGain (-(lvlDb + 12.0) * 0.75);
            test::approx ((double) a[(std::size_t) (t + look)], want, 5e-4,
                          "the impulse lands at t+lookahead with the key's own gain reduction (look " + std::to_string (look) + ")");
            double before = 0.0;
            for (int i = 0; i < t + look; ++i) before = std::max (before, (double) std::fabs (a[(std::size_t) i]));
            test::ok (before < 1e-9, "nothing before t+lookahead (no off-by-one on the key)");
        }
    }

    //==========================================================================
    // B. THE ACCEPTANCE SIGNAL — 50 Hz sine + 3 kHz click, high-passed key
    //==========================================================================
    // Fixture discipline (CLAUDE.md): the tone is TAPERED at both ends, because a tone that starts at
    // full scale is itself a broadband click IN THE KEY and would drive the gain at t=0; and the truth
    // of the fixture — what the filter actually does at 50 Hz and at 3 kHz — is MEASURED here rather
    // than assumed from the design equations.
    test::group ("a high-passed key: the click drives the gain reduction, the 50 Hz sine does not");
    {
        const int N = 48000, n0 = 24000, taper = (int) (0.010 * kFs);
        std::vector<float> io ((std::size_t) N, 0.0f);
        for (int i = 0; i < N; ++i)
        {
            double w = 1.0;
            if (i < taper)      w = 0.5 * (1.0 - std::cos (core::kPi * i / taper));
            if (i >= N - taper) w = 0.5 * (1.0 - std::cos (core::kPi * (N - 1 - i) / taper));
            io[(std::size_t) i] = (float) (w * 0.5 * std::sin (2.0 * core::kPi * 50.0 * i / kFs));
        }
        for (int k = 0; k < 48; ++k)                                    // 3 cycles of 3 kHz, zero at both ends
            io[(std::size_t) (n0 + k)] += (float) (0.9 * std::sin (core::kPi * k / 8.0));

        // The key filter lives HERE, in the caller — which is the whole architectural point of P2.
        const auto co = eq::matched::highpass (300.0, kFs, 0.70710678118654752);
        eq::Biquad hp; hp.setCoeffs (co);
        std::vector<float> key ((std::size_t) N);
        for (int i = 0; i < N; ++i) key[(std::size_t) i] = hp.processSample (io[(std::size_t) i]);

        auto magAt = [&] (double f) {
            eq::Biquad b; b.setCoeffs (co); double mx = 0.0;
            const int M = (int) (kFs / f * 40);
            for (int i = 0; i < M; ++i) { const float y = b.processSample ((float) std::sin (2.0 * core::kPi * f * i / kFs)); if (i > M / 2) mx = std::max (mx, (double) std::fabs (y)); }
            return mx;
        };
        test::approx (20.0 * std::log10 (magAt (50.0)),   -31.1, 1.0, "fixture: the key filter is -31 dB at 50 Hz");
        test::approx (20.0 * std::log10 (magAt (3000.0)),   0.0, 0.5, "fixture: the key filter is flat at 3 kHz");

        auto p = baseParams();
        p.detector = dynamics::Detector::Peak; p.link = dynamics::LinkMode::Max;
        p.thresholdDb = -12.0; p.ratio = 4.0; p.kneeDb = 0.0;
        p.attackMs = 0.1; p.releaseMs = 50.0; p.lookaheadMs = 0.0; p.makeupDb = 0.0;

        auto run = [&] (bool external, std::vector<double>& gr) {
            auto a = io;
            dynamics::Compressor c; test::ok (c.prepare (kFs, N, 1, 5.0), "prepare"); c.setParams (p);
            gr.assign ((std::size_t) N, 0.0);
            for (int i = 0; i < N; ++i)
            {
                float* s[1] { a.data() + i };
                const float* k[1] { key.data() + i };
                if (external) c.process (s, 1, 1, k, 1); else c.process (s, 1, 1);
                gr[(std::size_t) i] = c.gainReductionDb();
            }
            return a;
        };

        std::vector<double> grExt, grSelf;
        const auto yExt  = run (true,  grExt);
        const auto ySelf = run (false, grSelf);

        // Control: WITHOUT the high-pass, the 50 Hz sine is what drives the compressor. This is the
        // defect P2 exists to remove, stated as a number.
        double selfDeepest = 0.0;
        for (int i = taper; i < n0 - 1000; ++i) selfDeepest = std::min (selfDeepest, grSelf[(std::size_t) i]);
        test::ok (selfDeepest < -3.0, "control: self-keyed, the sine alone pulls the gain down (measured " + std::to_string (selfDeepest) + " dB)");

        // With the high-passed key: not "small", EXACTLY nothing. The key sits 25 dB under a hard knee,
        // so the static curve returns +-0 and the applied gain is 1.0f exactly — which makes the whole
        // pre-click stretch a bit-exact null. Any leak of the programme into the detector breaks it.
        int nonZero = 0;
        for (int i = 0; i < n0; ++i) if (! (grExt[(std::size_t) i] == 0.0)) ++nonZero;
        test::ok (nonZero == 0, "keyed: gain reduction is exactly 0 dB for the whole 500 ms of sine before the click");
        std::size_t moved = 0;
        for (int i = 0; i < n0; ++i) if (std::memcmp (&yExt[(std::size_t) i], &io[(std::size_t) i], sizeof (float)) != 0) ++moved;
        test::ok (moved == 0, "keyed: and the output is the input over that stretch, bit for bit");

        // The click, on the other hand, is seen: the static target for 0.9 through the filter is
        // -8.83 dB, and a 0.1 ms attack reaches most of it inside the 48-sample burst.
        double deepest = 0.0; int argmin = -1;
        for (int i = n0; i < n0 + 480; ++i) if (grExt[(std::size_t) i] < deepest) { deepest = grExt[(std::size_t) i]; argmin = i; }
        test::ok (deepest < -6.5 && deepest > -9.0, "keyed: the click pulls the gain down to " + std::to_string (deepest) + " dB (static target -8.83)");
        test::ok (argmin >= n0 && argmin <= n0 + 64, "keyed: the deepest reduction sits inside the click, not elsewhere");
        test::ok (std::fabs (grExt[(std::size_t) (n0 + 5 * 2400)]) < 0.15, "keyed: released to within 0.15 dB after five release constants");

        // The headline: one signal, one compressor, the key the only difference. Self-keyed, the 50 Hz
        // sine alone earns real reduction and the click is buried in it; keyed, the sine earns exactly
        // none and the click earns more than the sine ever did.
        test::ok (deepest < selfDeepest,
                  "the key moved the compressor's attention from the sine (" + std::to_string (selfDeepest)
                  + " dB) to the click (" + std::to_string (deepest) + " dB)");
    }

    //==========================================================================
    // C. HOSTILE SEQUENCES
    //==========================================================================

    // The key appears in the middle of a stream and later goes away. Neither is an error and neither may
    // discontinue anything by itself: the detector is a running average and simply starts averaging a
    // different signal. Pinned against the oracle, which models exactly that.
    test::group ("the key appears mid-stream, then disappears");
    {
        const int blk = 256, blocks = 40;
        auto s = probeSignal (blk * blocks, 0.0);
        auto key = tone (blk * blocks, 800.0, 0.5);
        auto a = s, b = s;
        auto p = baseParams(); p.thresholdDb = -26.0; p.ratio = 4.0; p.lookaheadMs = 1.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, blk, 1, 10.0), "prepare"); c.setParams (p);
        Oracle o; o.prepare (kFs, 1, 10.0, p);
        double biggestJump = 0.0, prevGr = 0.0;
        for (int bi = 0; bi < blocks; ++bi)
        {
            const int off = bi * blk;
            float* ia[1] { a.data() + off }; float* ib[1] { b.data() + off };
            const float* k[1] { key.data() + off };
            const bool keyed = (bi >= 10 && bi < 25);                  // on at block 10, off at block 25
            if (keyed) { c.process (ia, 1, blk, k, 1); o.process (ib, 1, blk, k, 1); }
            else       { c.process (ia, 1, blk);        o.process (ib, 1, blk); }
            if (bi == 10 || bi == 25) biggestJump = std::max (biggestJump, std::fabs (c.gainReductionDb() - prevGr));
            prevGr = c.gainReductionDb();
            if (c.latencySamples() != 48) { test::ok (false, "latency moved at block " + std::to_string (bi)); break; }
        }
        test::ok (c.latencySamples() == 48, "latencySamples() never moved across the routing changes");
        test::ok (countDiff (a, b) == 0, "keyed/unkeyed alternation nulls the oracle, bit for bit");
        test::ok (allFinite (a), "and stays finite");
        test::ok (biggestJump < 40.0, "no state was silently dropped at the routing change (jump " + std::to_string (biggestJump) + " dB)");
    }

    // The key's own channel count changes while the stream runs. Only `key[0..keyNc-1]` may ever be
    // read; under ASan a regression here is an out-of-bounds read rather than a wrong number.
    test::group ("the key's channel count changes mid-stream");
    {
        const int blk = 128, blocks = 24;
        auto s0 = probeSignal (blk * blocks, 0.0), s1 = probeSignal (blk * blocks, 9.0);
        auto k0 = tone (blk * blocks, 300.0, 0.5), k1 = tone (blk * blocks, 900.0, 0.4);
        auto a0 = s0, a1 = s1, b0 = s0, b1 = s1;
        auto p = baseParams(); p.link = dynamics::LinkMode::MeanPower; p.thresholdDb = -28.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, blk, 2, 10.0), "prepare"); c.setParams (p);
        Oracle o; o.prepare (kFs, 2, 10.0, p);
        for (int bi = 0; bi < blocks; ++bi)
        {
            const int off = bi * blk;
            const int keyNc = (bi % 3 == 0) ? 1 : 2;                   // 1, 2, 2, 1, 2, 2, ...
            float* ia[2] { a0.data() + off, a1.data() + off };
            float* ib[2] { b0.data() + off, b1.data() + off };
            const float* k[2] { k0.data() + off, k1.data() + off };
            c.process (ia, 2, blk, k, keyNc);
            o.process (ib, 2, blk, k, keyNc);
            if (c.latencySamples() != 0) { test::ok (false, "latency moved at block " + std::to_string (bi)); break; }
        }
        test::ok (c.latencySamples() == 0, "latencySamples() never moved across the key-count changes");
        test::ok (countDiff (a0, b0) + countDiff (a1, b1) == 0, "a key count that changes every block nulls the oracle");
        test::ok (allFinite (a0) && allFinite (a1), "and stays finite");
    }

    // reset() has to clear EVERYTHING a keyed stream can leave behind, or a reused instance is not a
    // fresh one — and every host reuses instances.
    test::group ("reset() after a keyed stream leaves a fresh instance");
    {
        const int n = 2000;
        auto s = probeSignal (n, 0.0); auto key = tone (n, 250.0, 0.8);
        auto p = baseParams(); p.thresholdDb = -34.0; p.ratio = 8.0; p.lookaheadMs = 4.0;

        dynamics::Compressor used, fresh;
        test::ok (used.prepare (kFs, 512, 1, 10.0) && fresh.prepare (kFs, 512, 1, 10.0), "prepare");
        used.setParams (p); fresh.setParams (p);
        {   // dirty the instance: keyed blocks, unkeyed blocks, a loud key
            auto d = s; float* io[1] { d.data() }; const float* k[1] { key.data() };
            used.process (io, 1, n, k, 1);
            used.process (io, 1, n);
            used.reset();
        }
        auto a = s, b = s;
        float* ia[1] { a.data() }; float* ib[1] { b.data() };
        const float* k[1] { key.data() };
        used.process (ia, 1, n, k, 1);
        fresh.process (ib, 1, n, k, 1);
        test::ok (countDiff (a, b) == 0, "reset instance == fresh instance, bit for bit");
    }

    // A single poisoned KEY sample must not be able to end the stream. On the unguarded code this was
    // permanent in every mode: a bypass under DownCompress, silence under DownExpand, and +rangeDb of
    // BOOST under UpCompress. A finite 1e20 is in the list because it squares to +Inf in float.
    test::group ("one poisoned key sample does not end the stream");
    {
        for (float poison : { kInf, -kInf, kNan, 1.0e20f, -1.0e20f })
          for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
            for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
              for (int nch = 1; nch <= 2; ++nch)
              {
                  const int n = 24000;
                  auto p = baseParams();
                  p.detector = d; p.link = l; p.thresholdDb = -24.0; p.ratio = 4.0;
                  p.attackMs = 1.0; p.releaseMs = 30.0; p.rmsWindowMs = 5.0;

                  // a key 12 dB over threshold: the healthy answer is about -9 dB of reduction
                  auto key0 = tone (n, 500.0, core::dbToGain (-12.0) * std::sqrt (2.0));
                  auto key1 = key0;
                  key0[(std::size_t) 100] = poison;                    // the single bad sample
                  auto a0 = tone (n, 500.0, 0.5), a1 = tone (n, 500.0, 0.5);

                  dynamics::Compressor c; test::ok (c.prepare (kFs, n, nch), "prepare"); c.setParams (p);
                  float* io[2] { a0.data(), a1.data() };
                  const float* k[2] { key0.data(), key1.data() };
                  c.process (io, nch, n, k, nch);
                  test::ok (allFinite (a0) && (nch == 1 || allFinite (a1)), "output stays finite");
                  test::ok (c.gainReductionDb() < -6.0,
                            "the detector recovered and is still compressing (GR " + std::to_string (c.gainReductionDb()) + " dB)");
              }
    }

    // The same poison, but under UpCompress — the mode where the failure was not "nothing happens" but
    // "+60 dB, permanently".
    test::group ("a poisoned key under UpCompress does not become a permanent boost");
    {
        for (float poison : { kInf, kNan, 1.0e20f })
        {
            const int n = 24000;
            auto p = baseParams();
            p.detector = dynamics::Detector::Rms; p.mode = dynamics::Mode::UpCompress;
            p.thresholdDb = -20.0; p.ratio = 3.0; p.rangeDb = 60.0; p.attackMs = 5.0; p.releaseMs = 40.0;
            auto key = tone (n, 500.0, 0.5); key[(std::size_t) 100] = poison;
            auto a = tone (n, 500.0, 0.5);
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            float* io[1] { a.data() }; const float* k[1] { key.data() };
            c.process (io, 1, n, k, 1);
            test::ok (allFinite (a), "output stays finite");
            test::ok (c.gainReductionDb() < 1.0, "no runaway boost (GR " + std::to_string (c.gainReductionDb()) + " dB, cap was +60)");
            test::ok (rmsOf (a, n - 4000, n) < 2.0, "and the output level is sane");
        }
    }

    // A ratio of exactly 1 is the honest "1:1, off" setting AND the fallback for any invalid ratio. It
    // makes the curve's slope exactly zero, and an infinite level then computes Inf*0 = NaN, which the
    // range clamp's two ordered comparisons both let through.
    test::group ("ratio 1:1 with an infinite key sample stays finite");
    {
        for (double ratio : { 1.0, 0.5, std::numeric_limits<double>::quiet_NaN() })
        {
            const int n = 8192;
            auto p = baseParams(); p.detector = dynamics::Detector::Peak; p.ratio = ratio;
            auto key = tone (n, 500.0, 0.5); key[(std::size_t) 50] = kInf;
            auto a = tone (n, 500.0, 0.5);
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            float* io[1] { a.data() }; const float* k[1] { key.data() };
            c.process (io, 1, n, k, 1);
            test::ok (allFinite (a), "output finite at ratio " + std::to_string (ratio));
            test::ok (countDiff (a, tone (n, 500.0, 0.5)) == 0, "and 1:1 is transparent, bit for bit");
        }
    }

    // A non-finite PROGRAMME sample is a different matter: the compressor is not a sanitiser and passes
    // it on (there is no recursive state on the audio path). What must NOT happen is the detector state
    // being poisoned by it in the self-keyed case.
    test::group ("a non-finite programme sample does not poison the detector");
    {
        for (float poison : { kInf, kNan })
        {
            const int n = 24000;
            auto p = baseParams(); p.detector = dynamics::Detector::Rms; p.thresholdDb = -24.0; p.ratio = 4.0;
            auto a = tone (n, 500.0, core::dbToGain (-12.0) * std::sqrt (2.0));
            a[(std::size_t) 100] = poison;
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            float* io[1] { a.data() };
            c.process (io, 1, n);
            test::ok (c.gainReductionDb() < -6.0, "self-keyed detector recovered from a poisoned programme sample");
            bool tailFinite = true; for (int i = 1000; i < n; ++i) tailFinite &= (bool) std::isfinite (a[(std::size_t) i]);
            test::ok (tailFinite, "and everything after it is finite");
        }
    }

    // THE #119 SHAPE, transplanted. A channel that sits out blocks comes back holding audio from before
    // it left, and the lookahead delay is where it hides. Measured on the unguarded code: 0.75 out of
    // digital silence, 417 ms later.
    test::group ("stereo -> mono -> stereo does not re-emit the returning channel's old audio");
    {
        const int look = 480;
        auto p = baseParams(); p.thresholdDb = 24.0; p.lookaheadMs = 1000.0 * look / kFs;   // no compression: isolate the delay
        dynamics::Compressor c; test::ok (c.prepare (kFs, 4096, 2, 50.0), "prepare"); c.setParams (p);
        test::ok (c.latencySamples() == look, "lookahead as expected");

        std::vector<float> L (2048, 0.0f), R (2048, 0.0f);
        const int n1 = 1000;
        for (int i = n1 - look; i < n1; ++i) R[(std::size_t) i] = 0.75f;    // a marker only R will hold
        float* ch2[2] { L.data(), R.data() };
        c.process (ch2, 2, n1);

        std::vector<float> mono (20000, 0.0f);
        float* ch1[1] { mono.data() };
        c.process (ch1, 1, 20000);                                          // R sits out for 417 ms
        test::ok (c.latencySamples() == look, "latency unchanged by the channel-count change");

        std::vector<float> L3 (2048, 0.0f), R3 (2048, 0.0f);                // silence in
        float* ch3[2] { L3.data(), R3.data() };
        c.process (ch3, 2, 1024);
        double worst = 0.0;
        for (int i = 0; i < 1024; ++i) worst = std::max (worst, (double) std::fabs (R3[(std::size_t) i]));
        test::ok (worst == 0.0, "silence in, silence out on the returning channel (was 0.75)");

        // And the live channel must NOT have been reset along with it: the shared gain keeps running.
        std::vector<float> L4 (2048, 0.5f), R4 (2048, 0.5f);
        float* ch4[2] { L4.data(), R4.data() };
        c.process (ch4, 2, 1024);
        test::ok (countDiff (L4, R4) == 0, "both channels still receive the same gain (no image shift)");
    }

    // Switching the detector between Peak and Rms mid-stream. `env` is |x| in one mode and x^2 in the
    // other, so flipping the enum without converting reinterprets a level as a power — measured at
    // +20 dB on the primitive, and -8.5 dB of unearned reduction through the compressor.
    test::group ("switching the detector mid-stream carries no stale level");
    {
        const int n = (int) kFs * 2;
        auto p = baseParams();
        p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 10.0; p.thresholdDb = -40.0; p.ratio = 4.0;
        p.attackMs = 5.0; p.releaseMs = 40.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);

        auto loud = tone (n, 400.0, 0.7);
        { auto a = loud; float* io[1] { a.data() }; c.process (io, 1, n); }
        test::ok (c.gainReductionDb() < -10.0, "loud RMS pass compressed");

        p.detector = dynamics::Detector::Peak; c.setParams (p);
        auto quiet = tone (n * 2, 400.0, core::dbToGain (-80.0));
        { auto a = quiet; float* io[1] { a.data() }; c.process (io, 1, n * 2); }

        p.detector = dynamics::Detector::Rms; c.setParams (p);
        { auto a = quiet; float* io[1] { a.data() }; c.process (io, 1, n * 2); }
        test::ok (std::fabs (c.gainReductionDb()) < 0.5,
                  "-80 dBFS material earns no reduction after Rms -> Peak -> Rms (was -8.5 dB)");
    }

    // Changing the lookahead moves the read pointer of a live ring: raising it re-emits audio already
    // delivered, lowering it skips over some. Measured on the unguarded code: 48 samples of 0.9 into
    // silence. The reported latency moves too, which is what makes this a resync rather than automation.
    test::group ("changing the lookahead mid-stream does not re-emit delivered audio");
    {
        const int n = 512;
        auto p = baseParams(); p.thresholdDb = 24.0; p.lookaheadMs = 0.0;   // no compression: isolate the delay
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1, 10.0), "prepare"); c.setParams (p);
        std::vector<float> loud ((std::size_t) n, 0.9f);
        { auto a = loud; float* io[1] { a.data() }; c.process (io, 1, n); }
        test::ok (c.latencySamples() == 0, "lookahead 0 reported");

        p.lookaheadMs = 1.0; c.setParams (p);                               // 48 samples
        test::ok (c.latencySamples() == 48, "the new lookahead is reported");
        std::vector<float> silence ((std::size_t) n, 0.0f);
        { float* io[1] { silence.data() }; c.process (io, 1, n); }
        double worst = 0.0; for (float v : silence) worst = std::max (worst, (double) std::fabs (v));
        test::ok (worst == 0.0, "silence in, silence out after the lookahead change (was 0.9)");
    }

    // The strongest form of "the gate works": a poisoned key must produce EXACTLY the stream a key with
    // the sanitised value produces. `isfinite` on the output only says the damage was survivable; this
    // says there was no damage — the gate is a substitution, not a repair.
    test::group ("a poisoned key sample equals the explicitly sanitised one, bit for bit");
    {
        const std::pair<float, float> subs[] {
            { kInf, 0.0f }, { -kInf, 0.0f }, { kNan, 0.0f },
            { 1.0e20f, 1.0e6f }, { -1.0e20f, -1.0e6f }, { 2.0e6f, 1.0e6f }
        };
        for (auto [bad, good] : subs)
          for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
            for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
            {
                const int n = 6000, nch = 2;
                auto p = baseParams(); p.detector = d; p.link = l; p.thresholdDb = -30.0; p.ratio = 4.0;
                auto k0 = tone (n, 500.0, 0.4), k1 = tone (n, 310.0, 0.3);
                auto k0b = k0; k0b[(std::size_t) 77] = bad;
                auto k0g = k0; k0g[(std::size_t) 77] = good;
                auto a0 = probeSignal (n, 0.0), a1 = probeSignal (n, 4.0);
                auto b0 = a0, b1 = a1;

                dynamics::Compressor ca, cb;
                test::ok (ca.prepare (kFs, n, nch) && cb.prepare (kFs, n, nch), "prepare");
                ca.setParams (p); cb.setParams (p);
                float* ia[2] { a0.data(), a1.data() }; const float* ka[2] { k0b.data(), k1.data() };
                float* ib[2] { b0.data(), b1.data() }; const float* kb[2] { k0g.data(), k1.data() };
                ca.process (ia, nch, n, ka, nch);
                cb.process (ib, nch, n, kb, nch);
                test::ok (countDiff (a0, b0) + countDiff (a1, b1) == 0,
                          "poisoned key == sanitised key, bit for bit");
            }
    }

    // A non-null key POINTER with a count of zero or less is self-keyed, exactly as `nullptr` is. A
    // caller that computes the count can reach this; nothing else in the suite passes that pair.
    test::group ("a non-null key with no channels is the self-keyed path");
    {
        for (int keyNc : { 0, -1, -7 })
        {
            const int n = 3000;
            auto p = baseParams(); p.thresholdDb = -26.0; p.ratio = 4.0; p.lookaheadMs = 1.0;
            auto s0 = probeSignal (n, 0.0);
            auto loudKey = tone (n, 700.0, 0.95);            // would give VERY different gain if read
            auto a = s0, b = s0;
            dynamics::Compressor ca, cb;
            test::ok (ca.prepare (kFs, n, 1, 10.0) && cb.prepare (kFs, n, 1, 10.0), "prepare");
            ca.setParams (p); cb.setParams (p);
            float* ia[1] { a.data() }; float* ib[1] { b.data() };
            const float* k[1] { loudKey.data() };
            ca.process (ia, 1, n);                            // three-argument form
            cb.process (ib, 1, n, k, keyNc);                  // non-null key, count <= 0
            test::ok (countDiff (a, b) == 0,
                      "key with count " + std::to_string (keyNc) + " == self-keyed, bit for bit");
        }
    }

    // The key's channel count is the KEY's, with no relation to the programme's and no silent clamp to
    // kMaxChannels. A 17-channel key whose only content is in channel 17 proves both at once: the count
    // is honoured past the house channel limit, and every plane is read.
    test::group ("a key with more channels than kMaxChannels is read in full");
    {
        const int n = 4096, keyNc = core::kMaxChannels + 1;
        std::vector<std::vector<float>> planes ((std::size_t) keyNc, std::vector<float> ((std::size_t) n, 0.0f));
        planes[(std::size_t) (keyNc - 1)] = tone (n, 500.0, 0.9);      // only the LAST plane carries signal
        std::vector<const float*> ptrs ((std::size_t) keyNc);
        for (int c = 0; c < keyNc; ++c) ptrs[(std::size_t) c] = planes[(std::size_t) c].data();

        auto p = baseParams(); p.detector = dynamics::Detector::Peak; p.link = dynamics::LinkMode::Max;
        p.thresholdDb = -20.0; p.ratio = 4.0; p.attackMs = 1.0;
        auto a = tone (n, 400.0, 0.2); float* io[1] { a.data() };
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
        c.process (io, 1, n, ptrs.data(), keyNc);
        test::ok (c.gainReductionDb() < -5.0, "the 17th key channel drove the gain (count is not clamped)");
    }

    // Detector behaviour must not depend on how the caller cut the stream into blocks — the property an
    // offline analysis pass and a real-time pass have to share, and the reason P3 can trust either.
    // The flush is what breaks it: it fires once per call, and in RMS mode `env` is a POWER, so the
    // house 1e-15 threshold reaches up to an amplitude of -150 dBFS.
    test::group ("the detector does not depend on the caller's block partition");
    {
        for (double amp : { 1.0e-8, 1.0e-7, 3.0e-8, 1.0e-6 })
        {
            const int n = 10000;
            auto p = baseParams();
            p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 5.0;
            p.thresholdDb = -170.0; p.ratio = 4.0; p.kneeDb = 0.0; p.rangeDb = 60.0;
            p.attackMs = 0.0; p.releaseMs = 0.0;
            double gr[2];
            for (int mode = 0; mode < 2; ++mode)
            {
                std::vector<float> a ((std::size_t) n, (float) amp);
                dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
                if (mode == 0) { float* io[1] { a.data() }; c.process (io, 1, n); }
                else for (int i = 0; i < n; ++i) { float* io[1] { a.data() + i }; c.process (io, 1, 1); }
                gr[mode] = c.gainReductionDb();
            }
            test::approx (gr[0], gr[1], 1e-9,
                          "one call == n one-sample calls at " + std::to_string (20.0 * std::log10 (amp)) + " dBFS");
        }
    }

    // Finite parameters, permanent damage: a huge but finite range let a finite threshold produce a
    // delta of ~1e300, which is -Inf once cast to float; the GR follower then computed `-Inf + c*(0 -
    // -Inf)` = NaN and never came back, so ordinary parameters restored afterwards still produced NaN.
    test::group ("finite but extreme parameters cannot poison the gain permanently");
    {
        const int n = 4096;
        auto sane = baseParams(); sane.thresholdDb = -24.0; sane.ratio = 4.0; sane.rangeDb = 40.0;
        auto wild = sane; wild.thresholdDb = -1.0e300; wild.rangeDb = 1.0e300;

        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare");
        auto a = tone (n, 400.0, 0.5); float* io[1] { a.data() };
        c.setParams (wild); c.process (io, 1, n);
        test::ok (allFinite (a), "the extreme-parameter block itself stays finite");
        auto b = tone (n, 400.0, 0.5); float* io2[1] { b.data() };
        c.setParams (sane); c.process (io2, 1, n);
        test::ok (allFinite (b), "and ordinary parameters afterwards produce finite audio");
        test::ok (std::isfinite (c.gainReductionDb()), "the gain-reduction state recovered");

        // A makeup gain of +1000 dB is 1e50 — +Inf in float — and `0.0f * Inf` is NaN, so a silent
        // passage came out poisoned rather than silent.
        auto huge = sane; huge.makeupDb = 1000.0;
        dynamics::Compressor c2; test::ok (c2.prepare (kFs, n, 1), "prepare");
        c2.setParams (huge);
        std::vector<float> quiet ((std::size_t) n, 0.0f); float* io3[1] { quiet.data() };
        c2.process (io3, 1, n);
        test::ok (allFinite (quiet), "an absurd makeup does not turn silence into NaN");

        // THE SUM, which bounding each term on its own does not bound: an UpCompress curve parked at
        // +rangeDb plus a large makeup is +800 dB = 1e40, +Inf in float, and `0.0f * Inf` is NaN.
        auto both = sane;
        both.detector = dynamics::Detector::Peak; both.mode = dynamics::Mode::UpCompress;
        both.thresholdDb = 1000.0; both.ratio = 2.0; both.kneeDb = 0.0; both.rangeDb = 400.0;
        both.makeupDb = 400.0; both.attackMs = 0.0; both.releaseMs = 0.0;
        dynamics::Compressor c3; test::ok (c3.prepare (kFs, n, 1), "prepare");
        c3.setParams (both);
        std::vector<float> sil ((std::size_t) n, 0.0f), key ((std::size_t) n, 0.5f);
        float* io4[1] { sil.data() }; const float* k4[1] { key.data() };
        c3.process (io4, 1, n, k4, 1);
        test::ok (allFinite (sil), "a +rangeDb boost PLUS a +400 dB makeup still produces finite audio");
        test::ok (std::isfinite (c3.gainReductionDb()), "and finite gain-reduction state");
    }

    // A refused call must be a GAP, not a silent partial pass: nothing advanced, and the stream picks up
    // exactly where it left off when a well-formed call arrives.
    test::group ("a refused call leaves the stream exactly where it was");
    {
        const int n = 256;
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 6.0; p.lookaheadMs = 1.0;
        auto s = probeSignal (n * 4, 0.0);

        dynamics::Compressor a, b;
        test::ok (a.prepare (kFs, n, 1, 10.0) && b.prepare (kFs, n, 1, 10.0), "prepare");
        a.setParams (p); b.setParams (p);
        auto ya = s, yb = s;
        std::vector<float> spare ((std::size_t) n, 0.7f);
        for (int bi = 0; bi < 4; ++bi)
        {
            float* ia[1] { ya.data() + bi * n }; float* ib[1] { yb.data() + bi * n };
            if (bi == 2)                                     // one refused call in the middle, for `a` only
            {
                float* two[2] { spare.data(), spare.data() };
                a.process (two, 2, n);                       // 2 channels on a 1-channel instance: refused
            }
            a.process (ia, 1, n);
            b.process (ib, 1, n);
        }
        test::ok (countDiff (ya, yb) == 0, "a refused call changed nothing about the stream that followed");
    }

    // An absurdly long time constant must be very slow, not DEAD. `exp(-1/(t*fs))` rounds to exactly
    // 1.0f around 1e6 ms at 48 kHz, and a coefficient of 1 makes `x = x` — the follower stops for good.
    test::group ("an absurd time constant is slow, not frozen");
    {
        const int n = 8192;
        for (double relMs : { 1.0e6, 1.0e9, std::numeric_limits<double>::infinity() })
        {
            auto p = baseParams();
            p.detector = dynamics::Detector::Peak; p.thresholdDb = -30.0; p.ratio = 4.0;
            p.attackMs = 0.5; p.releaseMs = relMs;
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            auto loud = tone (n, 400.0, 0.8); float* io[1] { loud.data() };
            c.process (io, 1, n);
            test::ok (c.gainReductionDb() < -5.0, "a loud pass still earns reduction with releaseMs = " + std::to_string (relMs));
            test::ok (std::isfinite (c.gainReductionDb()), "and the state stays finite");
        }
    }

    // The partial reset has to be PARTIAL. A blunt reset() on the channel-count change would also
    // satisfy "the returning channel is silent" — and would additionally throw away the surviving
    // channel's in-flight lookahead audio and snap the shared gain to 0 dB. Compare the survivor
    // against a run that never changed shape: it must be untouched.
    test::group ("a channel-count change does not disturb the channels that stayed");
    {
        const int look = 480, blk = 1024;
        auto p = baseParams();
        p.detector = dynamics::Detector::Peak; p.thresholdDb = -20.0; p.ratio = 4.0;
        p.attackMs = 2.0; p.releaseMs = 80.0; p.lookaheadMs = 1000.0 * look / kFs;

        auto L = probeSignal (blk * 4, 0.0), R = probeSignal (blk * 4, 6.0);

        // reference: always stereo
        auto rL = L, rR = R;
        dynamics::Compressor cr; test::ok (cr.prepare (kFs, blk, 2, 50.0), "prepare"); cr.setParams (p);
        for (int bi = 0; bi < 4; ++bi) { float* io[2] { rL.data() + bi * blk, rR.data() + bi * blk }; cr.process (io, 2, blk); }

        // under test: stereo, stereo, MONO (R sits out), stereo
        auto tL = L, tR = R;
        dynamics::Compressor ct; test::ok (ct.prepare (kFs, blk, 2, 50.0), "prepare"); ct.setParams (p);
        double grBefore = 0.0, grAfter = 0.0;
        for (int bi = 0; bi < 4; ++bi)
        {
            float* io[2] { tL.data() + bi * blk, tR.data() + bi * blk };
            if (bi == 2) { grBefore = ct.gainReductionDb(); ct.process (io, 1, blk); grAfter = ct.gainReductionDb(); }
            else ct.process (io, 2, blk);
        }
        // The mono block itself sees only L, so L can differ from the reference from there on; what must
        // NOT happen is the survivor losing its in-flight audio at the moment R leaves.
        test::ok (countDiff (std::vector<float> (tL.begin(), tL.begin() + 2 * blk),
                             std::vector<float> (rL.begin(), rL.begin() + 2 * blk)) == 0,
                  "the surviving channel is bit-identical up to the change");
        test::ok (std::fabs (grAfter - grBefore) < 12.0,
                  "the shared gain kept running across the change (jump " + std::to_string (grAfter - grBefore) + " dB, a full reset would snap to 0)");
        double energy = 0.0;
        for (int i = 2 * blk; i < 2 * blk + look; ++i) energy += std::fabs ((double) tL[(std::size_t) i]);
        test::ok (energy > 1e-3, "the survivor's lookahead pipeline was not emptied (no hole of silence)");
    }

    // A coefficient of exactly 1.0f is not "very slow", it is DEAD: `x = x` forever. In the detector
    // that means the envelope never leaves zero and NOTHING is ever compressed; in the gain-reduction
    // follower it means the reduction never leaves zero on attack, and never returns on release.
    test::group ("an absurd time constant creeps; it does not stop");
    {
        {   // detector: an RMS window of 1e6 ms rounds exp() to 1.0f at 48 kHz
            const int n = 8192;
            auto p = baseParams();
            p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 1.0e6;
            p.thresholdDb = -100.0; p.ratio = 4.0; p.kneeDb = 0.0; p.attackMs = 0.0; p.releaseMs = 0.0;
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            auto loud = tone (n, 400.0, 0.8); float* io[1] { loud.data() };
            c.process (io, 1, n);
            test::ok (c.gainReductionDb() < -0.5,
                      "the envelope still climbed out of zero with a 1e6 ms window (GR " + std::to_string (c.gainReductionDb()) + " dB)");
        }
        {   // GR follower on ATTACK: an attack of 1e6 ms must still move the gain
            const int n = 8192;
            auto p = baseParams();
            p.detector = dynamics::Detector::Peak; p.thresholdDb = -30.0; p.ratio = 4.0;
            p.attackMs = 1.0e6; p.releaseMs = 50.0;
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            auto loud = tone (n, 400.0, 0.8); float* io[1] { loud.data() };
            c.process (io, 1, n);
            test::ok (c.gainReductionDb() < 0.0, "the gain reduction left zero with a 1e6 ms attack");
        }
        {   // GR follower on RELEASE: an INFINITE release must fall back to instant, not freeze
            const int n = 8192;
            auto p = baseParams();
            p.detector = dynamics::Detector::Peak; p.thresholdDb = -30.0; p.ratio = 4.0;
            p.attackMs = 0.5; p.releaseMs = std::numeric_limits<double>::infinity();
            dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
            auto loud = tone (n, 400.0, 0.8); float* io[1] { loud.data() };
            c.process (io, 1, n);
            test::ok (c.gainReductionDb() < -5.0, "a loud pass earned reduction");
            std::vector<float> quiet ((std::size_t) n, 0.0f); float* io2[1] { quiet.data() };
            c.process (io2, 1, n);
            test::ok (c.gainReductionDb() > -0.5,
                      "and an infinite release still released (GR " + std::to_string (c.gainReductionDb()) + " dB, a frozen coefficient would hold it)");
        }
    }

    // GainComputer's own promise, tested on the primitive rather than through the compressor — the
    // gate makes this route unreachable from in there, but the curve is public and P3 drives it
    // directly with detector levels.
    test::group ("GainComputer returns a finite delta for every level");
    {
        for (double ratio : { 1.0, 0.5, 2.0, 1000.0 })
          for (auto m : { dynamics::Mode::DownCompress, dynamics::Mode::UpCompress, dynamics::Mode::DownExpand })
            for (double knee : { 0.0, 6.0 })
            {
                dynamics::GainComputer gc;
                gc.setMode (m); gc.setRatio (ratio); gc.setKneeDb (knee);
                gc.setThresholdDb (-18.0); gc.setRangeDb (60.0);
                for (double lvl : { -1.0e308, -240.0, -18.0, 0.0, 120.0, 1.0e308,
                                    std::numeric_limits<double>::infinity(),
                                    -std::numeric_limits<double>::infinity() })
                    test::ok (std::isfinite (gc.deltaDb (lvl)), "finite delta at an extreme level");
            }
        // and an absurd range cannot produce a delta that overflows a float
        dynamics::GainComputer gc;
        gc.setMode (dynamics::Mode::DownCompress); gc.setRatio (4.0); gc.setKneeDb (0.0);
        gc.setThresholdDb (-1.0e300); gc.setRangeDb (1.0e300);
        test::ok (std::isfinite ((float) gc.deltaDb (0.0)), "an extreme threshold/range pair survives the float cast");
    }

    // Law 8, as a PROPERTY and not as "the method exists". The detector is a one-pole that decays
    // toward zero and never reaches it in finite time on its own, so on a long silence its state walks
    // into the subnormal range and stays there — 10-100x the cost per sample on any CPU without
    // hardware FTZ. What pins the flush is that the state is EXACTLY zero afterwards. (This repo has
    // found three "the method is there, nobody calls it" defects; only a state assertion catches them.)
    test::group ("law 8: the detector state reaches exact zero on silence");
    {
        const int n = 200000;
        auto p = baseParams();
        p.detector = dynamics::Detector::Rms; p.rmsWindowMs = 50.0; p.thresholdDb = -60.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1), "prepare"); c.setParams (p);
        auto loud = tone (4096, 400.0, 0.8); float* io[1] { loud.data() };
        c.process (io, 1, 4096);
        test::ok (c.detectorLevel() > 0.1f, "the detector is charged");
        std::vector<float> silence ((std::size_t) n, 0.0f); float* io2[1] { silence.data() };
        c.process (io2, 1, n);
        test::ok (c.detectorLevel() == 0.0f, "and is exactly zero after silence, not a subnormal tail");
        test::ok (c.gainReductionDb() == 0.0, "the gain-reduction state landed on exact zero too");
    }

    //==========================================================================
    // D. THE PREPARED CONTRACT
    //==========================================================================

    test::group ("prepare() refuses what it cannot honour");
    {
        dynamics::Compressor c;
        test::ok (! c.prepare (std::numeric_limits<double>::quiet_NaN(), 512, 2), "NaN sample rate refused");
        test::ok (! c.prepare (0.0, 512, 2), "zero sample rate refused");
        test::ok (! c.prepare (-48000.0, 512, 2), "negative sample rate refused");
        test::ok (! c.prepare (1.0e9, 512, 2), "absurd sample rate refused");
        test::ok (! c.prepare (kFs, 0, 2), "zero maxBlock refused");
        test::ok (! c.prepare (kFs, 512, 0), "zero channels refused");
        test::ok (! c.prepare (kFs, 512, core::kMaxChannels + 1), "more than kMaxChannels refused, not clamped");
        test::ok (! c.prepare (kFs, 512, 2, std::numeric_limits<double>::quiet_NaN()), "NaN maxLookahead refused");
        test::ok (! c.prepare (kFs, 512, 2, -1.0), "negative maxLookahead refused");
        test::ok (! c.isPrepared(), "and none of that left it looking prepared");
        test::ok (c.latencySamples() == 0, "an unprepared compressor reports no latency");

        // An unprepared instance must do NOTHING, not something.
        std::vector<float> a (256, 0.5f); const auto in = a;
        float* io[1] { a.data() };
        c.process (io, 1, 256);
        test::ok (countDiff (a, in) == 0, "an unprepared process() leaves the buffer untouched");

        test::ok (c.prepare (kFs, 512, 2), "a sane configuration is accepted");
        test::ok (c.isPrepared(), "and reports itself prepared");
    }

    // More channels than prepared: REFUSED, and the reason is the module's central promise. Processing
    // a prefix leaves the surplus channels neither gained nor delayed — a level mismatch plus a
    // lookahead-sized skew, which is the image shift the channel link exists to prevent.
    test::group ("more channels than prepared are refused, not half-processed");
    {
        const int n = 512;
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 8.0; p.lookaheadMs = 1.0;
        dynamics::Compressor c; test::ok (c.prepare (kFs, n, 1, 10.0), "prepare for mono"); c.setParams (p);
        std::vector<float> L ((std::size_t) n, 0.9f), R ((std::size_t) n, 0.9f);
        const auto Lin = L, Rin = R;
        float* io[2] { L.data(), R.data() };
        c.process (io, 2, n);
        test::ok (countDiff (L, Lin) == 0 && countDiff (R, Rin) == 0,
                  "a stereo buffer on a mono-prepared instance is returned untouched");
    }

    test::group ("no allocation in the keyed process()");
    {
        const int n = 512;
        std::vector<float> a ((std::size_t) n, 0.2f), b ((std::size_t) n, 0.2f), k ((std::size_t) n, 0.4f);
        float* io[2] { a.data(), b.data() };
        const float* key[1] { k.data() };
        dynamics::Compressor c;
        test::ok (c.prepare (kFs, n, 2, 10.0), "prepare");
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 4.0; p.lookaheadMs = 2.0;
        c.setParams (p);
        c.process (io, 2, n, key, 1);                       // warm every branch before counting
        const long before = g_allocs.load();
        c.process (io, 2, n, key, 1);
        c.process (io, 2, n);
        c.process (io, 2, n, key, 1);
        const long after = g_allocs.load();
        test::okNoAlloc (after == before, "keyed and unkeyed process() performed zero heap allocations");
    }

    // THE MISSING PROPERTY, and the one products actually depend on: `setParams()` with UNCHANGED
    // parameters must be a no-op. Every consumer here pushes params every block, so anything in
    // apply() that is not idempotent runs hundreds of times a second. Two ways it was not: converting
    // the detector state on every call (in Rms that squares a power per block, `env -> env^(2^k)`,
    // which collapses to zero and silently stops the compressor), and clearing the delay rings on
    // every call (a hole of `lookahead` zeros per block).
    test::group ("setParams() with unchanged parameters is a no-op");
    {
        for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
          for (double look : { 0.0, 2.0 })
          {
              const int blk = 64, blocks = 60;
              auto p = baseParams();
              p.detector = d; p.rmsWindowMs = 5.0; p.thresholdDb = -26.0; p.ratio = 4.0; p.lookaheadMs = look;
              auto s0 = probeSignal (blk * blocks, 0.0);
              auto a = s0, b = s0;

              dynamics::Compressor ca, cb;
              test::ok (ca.prepare (kFs, blk, 1, 10.0) && cb.prepare (kFs, blk, 1, 10.0), "prepare");
              ca.setParams (p); cb.setParams (p);
              for (int bi = 0; bi < blocks; ++bi)
              {
                  float* ia[1] { a.data() + bi * blk }; float* ib[1] { b.data() + bi * blk };
                  ca.process (ia, 1, blk);
                  cb.setParams (p);                                  // ...every single block
                  cb.process (ib, 1, blk);
              }
              test::ok (countDiff (a, b) == 0, "params pushed every block == pushed once, bit for bit");
              test::ok (rmsOf (b, blk * blocks / 2, blk * blocks) > 1e-3, "and it is still compressing, not collapsed to silence");
          }
    }

    // The survivor must be intact THROUGH AND AFTER the return block, not merely up to the departure.
    // A delay-only blunt reset (clear every line on the increase, rather than only the idle ones)
    // satisfies every other assertion here and zeroes the survivor's in-flight lookahead audio.
    test::group ("the surviving channel keeps its in-flight audio across the return");
    {
        const int look = 64, blk = 512;
        auto p = baseParams(); p.thresholdDb = 24.0; p.lookaheadMs = 1000.0 * look / kFs;   // gain exactly 1
        dynamics::Compressor c; test::ok (c.prepare (kFs, blk, 2, 10.0), "prepare"); c.setParams (p);
        test::ok (c.latencySamples() == look, "lookahead as expected");

        auto L = probeSignal (blk * 4, 0.0), R = probeSignal (blk * 4, 6.0);
        const auto Lin = L;
        for (int bi = 0; bi < 4; ++bi)
        {
            float* io[2] { L.data() + bi * blk, R.data() + bi * blk };
            c.process (io, bi == 2 ? 1 : 2, blk);                    // block 2 is mono; block 3 is the RETURN
        }
        // With gain == 1 the compressor is a pure delay, so L must be its own input shifted by `look`
        // at EVERY index past the first `look` samples — including across the return block.
        std::size_t wrong = 0;
        for (int i = look; i < blk * 4; ++i)
            if (std::memcmp (&L[(std::size_t) i], &Lin[(std::size_t) (i - look)], sizeof (float)) != 0) ++wrong;
        test::ok (wrong == 0, "the survivor is its own input delayed by the lookahead, through the return block");
    }

    // A refused re-prepare must leave the instance UNPREPARED on the stale configuration, not quietly
    // prepared on it. Testing refusal only on a fresh instance cannot see the difference.
    test::group ("a refused re-prepare unprepares a working instance");
    {
        const int n = 512;
        dynamics::Compressor c;
        test::ok (c.prepare (kFs, n, 2, 10.0), "prepare succeeds");
        auto p = baseParams(); p.thresholdDb = -30.0; p.ratio = 6.0; p.lookaheadMs = 2.0;
        c.setParams (p);
        auto warm = probeSignal (n, 0.0); float* io[2] { warm.data(), warm.data() };
        {   // charge the state on a real stream
            auto a = probeSignal (n, 0.0), b = probeSignal (n, 3.0);
            float* w[2] { a.data(), b.data() };
            c.process (w, 2, n);
        }
        test::ok (c.isPrepared() && c.latencySamples() == 96, "prepared and reporting its latency");

        test::ok (! c.prepare (std::numeric_limits<double>::quiet_NaN(), n, 2), "the re-prepare is refused");
        test::ok (! c.isPrepared(), "and the instance is now unprepared");
        test::ok (c.latencySamples() == 0, "reporting no latency rather than a stale one");
        auto before = probeSignal (n, 5.0); auto after = before;
        float* q[2] { after.data(), after.data() };
        c.process (q, 2, n);
        test::ok (countDiff (after, before) == 0, "and process() is a byte-for-byte no-op");
        (void) io;
    }

    // `processSample` is a public entry point with its own gate; nothing exercised it.
    test::group ("LinkedDetector::processSample matches the one-channel frame path");
    {
        for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
        {
            dynamics::DetectorParams dp; dp.detector = d; dp.link = dynamics::LinkMode::Max; dp.rmsWindowMs = 4.0;
            dynamics::LinkedDetector a, b;
            a.prepare (kFs); b.prepare (kFs); a.setParams (dp); b.setParams (dp);
            auto x = probeSignal (4000, 0.0);
            x[10] = kInf; x[500] = kNan; x[1200] = 1.0e20f; x[2000] = -3.0e6f;   // through the gate
            const float* ch[1] { x.data() };
            std::size_t diff = 0;
            for (int i = 0; i < 4000; ++i)
            {
                const float u = a.processSample (x[(std::size_t) i]);
                const float v = b.process (ch, 1, i);
                if (std::memcmp (&u, &v, sizeof (float)) != 0) ++diff;
            }
            test::ok (diff == 0, "processSample == process(key,1,i), bit for bit, poison included");
        }
    }

    // A count of zero must not dereference plane 0 — the spelling an offline caller reaches for.
    test::group ("LinkedDetector with no key channels reads nothing");
    {
        dynamics::DetectorParams dp; dp.detector = dynamics::Detector::Peak;
        dynamics::LinkedDetector d; d.prepare (kFs); d.setParams (dp);
        test::ok (dynamics::linkAmplitudeGated (dynamics::LinkMode::Max, nullptr, 0, 0) == 0.0f,
                  "a zero count returns silence without touching the pointer");
        test::ok (dynamics::linkAmplitude (dynamics::LinkMode::MeanPower, nullptr, -3, 0) == 0.0f,
                  "and so does a negative one");
        test::ok (d.process (nullptr, 0, 0) == 0.0f, "the detector agrees");
    }

    // Independent, ANALYTIC checks — the Oracle above shares every primitive with the implementation,
    // so a bug inside one of them would null on both sides. These do not: each is a closed-form value.
    test::group ("analytic checks that share no code with the implementation");
    {
        // (a) MeanPower over a frame is sqrt(mean of squares) over the KEY's own count.
        for (int nc : { 1, 2, 3, 5 })
        {
            std::vector<std::vector<float>> planes ((std::size_t) nc, std::vector<float> (1, 0.0f));
            planes[(std::size_t) (nc - 1)][0] = 1.0f;                       // one hot lane
            std::vector<const float*> ptrs ((std::size_t) nc);
            for (int c = 0; c < nc; ++c) ptrs[(std::size_t) c] = planes[(std::size_t) c].data();
            const double want = std::sqrt (1.0 / nc);
            test::approx ((double) dynamics::linkAmplitudeGated (dynamics::LinkMode::MeanPower, ptrs.data(), nc, 0),
                          want, 1e-6, "MeanPower over " + std::to_string (nc) + " key channels == sqrt(1/n)");
            test::approx ((double) dynamics::linkAmplitudeGated (dynamics::LinkMode::Max, ptrs.data(), nc, 0),
                          1.0, 1e-6, "Max over the same frame == 1");
        }

        // (b) The RMS envelope's step response, in closed form, at THREE sample rates — the coefficient
        // is the only place the rate enters, and every other fixture in this file runs at 48 kHz.
        for (double fs : { 44100.0, 48000.0, 96000.0 })
          for (double winMs : { 2.0, 25.0 })
          {
              dynamics::DetectorParams dp;
              dp.detector = dynamics::Detector::Rms; dp.link = dynamics::LinkMode::Max; dp.rmsWindowMs = winMs;
              dynamics::LinkedDetector det; det.prepare (fs); det.setParams (dp);
              const float A = 0.6f;
              std::vector<float> step (60000, A);
              const float* ch[1] { step.data() };
              const int N = (int) (winMs * 0.001 * fs * 3.0);              // three time constants
              float lvl = 0.0f;
              for (int i = 0; i < N; ++i) lvl = det.process (ch, 1, i);
              // env is a one-pole on the POWER: env[n] = A^2 * (1 - c^n), level = sqrt(env)
              const double c = std::exp (-1.0 / (winMs * 0.001 * fs));
              const double want = (double) A * std::sqrt (1.0 - std::pow (c, (double) N));
              test::approx ((double) lvl, want, 1e-4,
                            "RMS step at " + std::to_string ((int) fs) + " Hz / " + std::to_string (winMs) + " ms == A*sqrt(1-c^n)");
          }

        // (c) The gain-reduction follower's step response, in closed form.
        {
            const double atkMs = 4.0, target = -9.0;
            dynamics::GainReductionFollower gr; gr.prepare (kFs); gr.setTimes (atkMs, 500.0);
            const int N = 300;
            float v = 0.0f;
            for (int i = 0; i < N; ++i) v = gr.process ((float) target);
            const double c = std::exp (-1.0 / (atkMs * 0.001 * kFs));
            test::approx ((double) v, target * (1.0 - std::pow (c, (double) N)), 1e-4, "GR step == t*(1-c^n)");
        }

        // (d) The static curve, from the textbook piecewise formula, hard knee.
        {
            dynamics::GainComputer gc;
            gc.setMode (dynamics::Mode::DownCompress); gc.setThresholdDb (-20.0);
            gc.setRatio (4.0); gc.setKneeDb (0.0); gc.setRangeDb (60.0);
            for (double lvl : { -40.0, -20.0, -12.0, 0.0, 20.0 })
            {
                const double over = lvl - (-20.0);
                const double want = -(over > 0.0 ? over : 0.0) * (1.0 - 1.0 / 4.0);
                test::approx (gc.deltaDb (lvl), want, 1e-12, "curve at " + std::to_string (lvl) + " dB");
            }
        }
    }

    // The four-argument call must NOT exist: a bare key pointer with no count is the shape that reads
    // past a mono key array on a stereo bus. A default argument on `numKeyChannels` would silently
    // bring it back, and no runtime test can see that.
    test::group ("a key pointer without a count does not compile");
    {
        static_assert (! detail::HasFourArgProcess<dynamics::Compressor>::value,
                       "process(io, nch, n, key) must not be callable — the count has to travel with the pointer");
        static_assert (detail::HasFiveArgProcess<dynamics::Compressor>::value,
                       "process(io, nch, n, key, keyCh) must be callable");
        test::ok (! detail::HasFourArgProcess<dynamics::Compressor>::value, "no four-argument process() exists");
        test::ok (detail::HasFiveArgProcess<dynamics::Compressor>::value, "the five-argument form does");
    }

    //==========================================================================
    // E. THE DETECTOR AS A STANDALONE OBJECT (what P3 will hold)
    //==========================================================================

    // The gated and ungated link must be the SAME function on healthy audio. They are generated from
    // one implementation so they cannot drift, and this is the assertion that keeps it that way.
    test::group ("the gated link is bit-transparent on finite audio");
    {
        const int n = 20000;
        auto c0 = probeSignal (n, 0.0), c1 = probeSignal (n, 13.0);
        for (int i = 0; i < n; i += 997) { c0[(std::size_t) i] *= 1.0e5f; c1[(std::size_t) i] *= -9.9e5f; }  // large, still inside +-1e6
        const float* ch[2] { c0.data(), c1.data() };
        std::size_t diff = 0;
        for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
            for (int nc = 1; nc <= 2; ++nc)
                for (int i = 0; i < n; ++i)
                {
                    const float a = dynamics::linkAmplitude (l, ch, nc, i);
                    const float b = dynamics::linkAmplitudeGated (l, ch, nc, i);
                    if (std::memcmp (&a, &b, sizeof (float)) != 0) ++diff;
                }
        test::ok (diff == 0, "gated == ungated for every finite sample within +-1e6");

        // and it is NOT transparent where it must not be
        std::vector<float> bad { kInf, kNan, 1.0e20f, -1.0e20f };
        const float* bch[1] { bad.data() };
        for (int i = 0; i < 4; ++i)
            test::ok (std::isfinite (dynamics::linkAmplitudeGated (dynamics::LinkMode::Max, bch, 1, i))
                      && dynamics::linkAmplitudeGated (dynamics::LinkMode::Max, bch, 1, i) <= 1.0e6f,
                      "the gate bounds a poisoned sample");
    }

    // The envelope an offline pass computes must be the envelope the compressor ran — not close to it,
    // the same. This is the P3 precondition stated directly on the detector object.
    test::group ("an offline LinkedDetector reproduces the compressor's own envelope");
    {
        for (auto d : { dynamics::Detector::Peak, dynamics::Detector::Rms })
          for (auto l : { dynamics::LinkMode::Max, dynamics::LinkMode::MeanPower })
          {
              const int n = 4000;
              auto k0 = probeSignal (n, 0.0), k1 = probeSignal (n, 6.0);
              auto a0 = tone (n, 400.0, 0.3), a1 = tone (n, 400.0, 0.3);
              auto p = baseParams(); p.detector = d; p.link = l; p.rmsWindowMs = 8.0;

              dynamics::Compressor c; test::ok (c.prepare (kFs, n, 2), "prepare"); c.setParams (p);
              dynamics::LinkedDetector off;
              off.prepare (kFs);
              off.setParams (p);                                   // the SAME params object, sliced

              const float* k[2] { k0.data(), k1.data() };
              std::size_t diff = 0;
              for (int i = 0; i < n; ++i)
              {
                  float* io[2] { a0.data() + i, a1.data() + i };
                  const float* ks[2] { k0.data() + i, k1.data() + i };
                  c.process (io, 2, 1, ks, 1 + 1);                 // one sample at a time, so the meter is per-sample
                  const float mine = off.process (k, 2, i);
                  const float theirs = c.detectorLevel();
                  if (std::memcmp (&mine, &theirs, sizeof (float)) != 0) ++diff;
              }
              test::ok (diff == 0, "offline detector == the compressor's detector, sample for sample, bit for bit");
          }
    }

    // The primitive that made the mode switch safe, tested on its own: the reported AMPLITUDE is
    // continuous across a Peak <-> Rms change instead of being reinterpreted.
    test::group ("EnvelopeFollower carries its level across a detector change");
    {
        dynamics::EnvelopeFollower f;
        f.prepare (kFs); f.setDetector (dynamics::Detector::Peak); f.setTimes (0.0, 0.0);
        f.process (0.25f);
        const float peakLevel = f.envelope();
        f.setDetector (dynamics::Detector::Rms);
        test::approx ((double) f.envelope(), (double) peakLevel, 1e-6, "Peak -> Rms keeps the amplitude (was a +20 dB step)");
        f.setDetector (dynamics::Detector::Peak);
        test::approx ((double) f.envelope(), (double) peakLevel, 1e-6, "and back again");

        // an infinite averaging window must not FREEZE the envelope
        dynamics::EnvelopeFollower g;
        g.prepare (kFs); g.setDetector (dynamics::Detector::Rms);
        g.setTimes (std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
        g.process (0.5f);
        test::approx ((double) g.envelope(), 0.5, 1e-6, "an infinite window falls back to instant, not to a frozen state");
    }

    return test::report();
}
