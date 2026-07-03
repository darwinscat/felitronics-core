// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the stereo module. The MonoBass battery is deliberately adversarial — each
// test tries to PROVE the filter does NOT work: measured side magnitude must equal the BLT-warped
// analytic LR4 curve (w + r⁴)/(1+r⁴), r = tan(πf/fs)/tan(πfc/fs), across the band AND at a high fc where
// the warp actually bites; a sample-exact reference NULL against a hand-rolled M/S → SVF-LR4 → blend →
// decode chain built from the same primitives; the mono-sum invariant on random stereo; bit-exact mono
// passthrough; correlation → +1 below fc; click-free bypass-boundary crossfades with a stale-state proof
// (re-entry from bypass must NOT replay old filter tails); fc/param abuse (near-Nyquist clamp, non-finite
// rejection); denormal flush to exact zero; strict numChannels==2 gating; and no allocation in process().
// Plus the Mid/Side matrix round-trip and the StereoWidth suite (mono-fold invariant et al.).

#include <felitronics_test.h>
#include <felitronics/core/Math.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/stereo/MidSide.h>
#include <felitronics/stereo/MonoBass.h>
#include <felitronics/stereo/StereoWidth.h>

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;

static double rmsTail (const std::vector<float>& v, int from)
{
    double e = 0.0; int c = 0;
    for (int i = from; i < (int) v.size(); ++i) { e += (double) v[i] * v[i]; ++c; }
    return c ? std::sqrt (e / (double) c) : 0.0;
}

// Deterministic uncorrelated stereo noise (shared by the MonoBass + StereoWidth batteries).
static void rngPair (int N, unsigned long long seed, std::vector<float>& L, std::vector<float>& R)
{
    unsigned long long s = seed; auto u = [&] { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (float) ((s >> 40) & 0xffff) / 32768.0f - 1.0f; };
    L.resize (N); R.resize (N); for (int i = 0; i < N; ++i) { L[i] = 0.6f * u(); R[i] = 0.6f * u(); }
}

// Steady-state amplitude of the f-Hz tone in v over [from, end) by quadrature projection. Callers keep
// an integer number of cycles in the window (f integral, window 1 s / 0.5 s) so there is no leakage bias.
static double toneAmp (const std::vector<float>& v, int from, double f, double sr)
{
    double a = 0.0, b = 0.0; const int nWin = (int) v.size() - from;
    for (int i = 0; i < nWin; ++i)
    {
        const double ph = 2.0 * core::kPi * f * (double) (from + i) / sr;
        a += v[from + i] * std::sin (ph); b += v[from + i] * std::cos (ph);
    }
    return 2.0 * std::sqrt (a * a + b * b) / (double) nWin;
}

// The analytic blended-side magnitude of the LR4 elliptical topology, at the BLT-warped frequency ratio
// r = tan(πf/fs)/tan(πfc/fs):  |S'/S| = (w + r⁴)/(1 + r⁴).  w=0 is the pure HP4 rolloff; w=1 is unity.
// (The TPT SVF is exactly the bilinear transform of the analog prototype, so this — not f/fc — is the
// truth the measurement must match; f/fc is only its low-frequency approximation.)
static double analyticSide (double f, double fc, double sr, double w)
{
    const double r  = std::tan (core::kPi * f / sr) / std::tan (core::kPi * fc / sr);
    const double r4 = r * r * r * r;
    return (w + r4) / (1.0 + r4);
}

// Drive a settled MonoBass with an anti-phase tone (pure Side, M = 0) and measure |S_out/S_in| steady-state.
static double measuredSide (double f, double fc, float w, double sr)
{
    stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency ((float) fc); mb.setLowWidth (w); mb.reset();
    const int warm = 24000, win = 48000, N = warm + win;
    std::vector<float> L (N), R (N);
    for (int i = 0; i < N; ++i) { const float v = 0.5f * (float) std::sin (2.0 * core::kPi * f * i / sr); L[i] = v; R[i] = -v; }
    float* io[2] { L.data(), R.data() };
    mb.process (io, 2, N);
    std::vector<float> side (N); for (int i = 0; i < N; ++i) side[i] = 0.5f * (L[i] - R[i]);
    return toneAmp (side, warm, f, sr) / 0.5;
}

int main()
{
    std::printf ("felitronics::stereo tests\n");
    const double sr = 48000.0, pi = core::kPi;

    // --- Mid/Side round-trip ---
    test::group ("MidSide identity");
    {
        float m, s, l, r;
        stereo::MidSide::encode (0.3f, -0.7f, m, s);
        stereo::MidSide::decode (m, s, l, r);
        test::approx (l,  0.3, 1e-6, "decode(encode) recovers L");
        test::approx (r, -0.7, 1e-6, "decode(encode) recovers R");
        test::approx (m, -0.2, 1e-6, "M = 1/2(L+R)");
        test::approx (s,  0.5, 1e-6, "S = 1/2(L-R)");
    }

    //==========================================================================
    // MonoBass — the adversarial battery.

    // --- the headline property: IN-PHASE content is NOT filtered, bit-exact (vs a full crossover, which would) ---
    test::group ("MonoBass: mono content passes UNFILTERED (bit-exact)");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.0f); mb.reset();          // default fc = 120
        test::ok (mb.frequency() == 120.0f, "default crossover is 120 Hz");
        const int N = 6000; std::vector<float> L (N), R (N), ref (N);
        for (int i = 0; i < N; ++i)
        {
            const float v = 0.5f * (float) std::sin (2.0 * pi * 50.0 * i / sr)
                          + 0.3f * (float) std::sin (2.0 * pi * 3000.0 * i / sr);
            L[i] = v; R[i] = v; ref[i] = v;
        }
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - ref[i])); md = std::max (md, (double) std::fabs (R[i] - ref[i])); }
        test::ok (md == 0.0, "L==R in -> out == in BIT-EXACT (S=0; float 1/2(L+L)=L is exact) — kept bass never filtered");

        std::vector<float> mono, dummy; rngPair (2048, 11, mono, dummy);                  // aliased dual-mono buffer
        auto buf = mono; float* alias[2] { buf.data(), buf.data() };
        mb.reset(); mb.process (alias, 2, 2048);
        md = 0; for (int i = 0; i < 2048; ++i) md = std::max (md, (double) std::fabs (buf[i] - mono[i]));
        test::ok (md == 0.0, "aliased io[0]==io[1] (dual-mono) -> untouched bit-exact");
    }

    // --- measured side magnitude == the BLT-warped analytic LR4 curve (wrong Q / order / warp all fail this) ---
    test::group ("MonoBass: measured == analytic side rolloff (BLT-warped LR4, lowWidth=0)");
    {
        for (double f : { 24.0, 40.0, 60.0, 120.0, 240.0, 480.0, 1200.0, 4000.0 })
        {
            const double got = measuredSide (f, 120.0, 0.0f, sr), want = analyticSide (f, 120.0, sr, 0.0);
            test::approx (got, want, 2e-5 + 0.005 * want, "side |H| at " + std::to_string ((int) f) + " Hz (fc=120)");
        }
        test::approx (measuredSide (120.0, 120.0, 0.0f, sr), 0.5, 1e-3, "exactly -6 dB at fc (the LR point)");
        test::ok (measuredSide (1200.0, 120.0, 0.0f, sr) > 0.99885, "transparent at 10x fc (within 0.01 dB of unity)");
        // High fc — where the BLT warp actually bites (r != f/fc by ~24% here): proves the SVF realises
        // the warped curve and that the analytic reference (not the naive f/fc one) is the right truth.
        for (double f : { 4000.0, 12000.0 })
        {
            const double got = measuredSide (f, 8000.0, 0.0f, sr), want = analyticSide (f, 8000.0, sr, 0.0);
            test::approx (got, want, 2e-5 + 0.005 * want, "warp-sensitive side |H| at " + std::to_string ((int) f) + " Hz (fc=8k)");
        }
    }

    // --- the lowWidth blend follows (w + r^4)/(1 + r^4) exactly — bump-free, monotone, exact endpoints ---
    test::group ("MonoBass: blend law (w + r^4)/(1 + r^4)");
    {
        for (double f : { 40.0, 120.0, 480.0 })
            for (float w : { 0.25f, 0.5f, 0.75f })
            {
                const double got = measuredSide (f, 120.0, w, sr), want = analyticSide (f, 120.0, sr, (double) w);
                test::approx (got, want, 2e-5 + 0.005 * want,
                              "blend |H| at " + std::to_string ((int) f) + " Hz, w=" + std::to_string (w));
            }
        const double a0 = measuredSide (40.0, 120.0, 0.0f, sr),  a25 = measuredSide (40.0, 120.0, 0.25f, sr),
                     a5 = measuredSide (40.0, 120.0, 0.5f, sr),  a75 = measuredSide (40.0, 120.0, 0.75f, sr);
        test::ok (a0 < a25 && a25 < a5 && a5 < a75 && a75 < 1.0, "below-fc side level strictly monotone in lowWidth");
    }

    // --- reference NULL: MonoBass == a hand-rolled M/S -> 4xSVF LR4 -> blend -> decode chain, sample-exact ---
    test::group ("MonoBass: reference NULL vs hand-rolled chain");
    {
        const float fcF = 97.3f, wF = 0.37f;
        const int N = 12000; std::vector<float> L0, R0; rngPair (N, 41, L0, R0);

        auto L = L0, R = R0;
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (fcF); mb.setLowWidth (wF); mb.reset();
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);

        constexpr double Q = 0.7071067811865476;              // the same Butterworth literal Crossover2 uses
        eq::Svf lp1, lp2, hp1, hp2;
        lp1.prepare (sr, 1); lp2.prepare (sr, 1); hp1.prepare (sr, 1); hp2.prepare (sr, 1);
        lp1.setParams (eq::FilterType::LowPass,  fcF, Q, 0.0); lp2.setParams (eq::FilterType::LowPass,  fcF, Q, 0.0);
        hp1.setParams (eq::FilterType::HighPass, fcF, Q, 0.0); hp2.setParams (eq::FilterType::HighPass, fcF, Q, 0.0);
        double md = 0;
        for (int i = 0; i < N; ++i)
        {
            const float xf = 0.0f;                            // settled, not full-wide -> pure wet (mirrors the impl)
            float m, s; stereo::MidSide::encode (L0[i], R0[i], m, s);
            const float lp = lp2.processSample (0, lp1.processSample (0, s));
            const float hp = hp2.processSample (0, hp1.processSample (0, s));
            const float wet  = wF * lp + hp;
            const float sOut = xf * s + (1.0f - xf) * wet;
            float l, r; stereo::MidSide::decode (m, sOut, l, r);
            md = std::max (md, (double) std::fabs (l - L[i])); md = std::max (md, (double) std::fabs (r - R[i]));
        }
        test::ok (md == 0.0, "null vs the same-primitive reference chain is SAMPLE-EXACT (Crossover2 dedup changed nothing)");
    }

    // --- the mono-sum invariant: 1/2(L'+R') == 1/2(L+R) — the Mid is never touched, even mid-fade ---
    test::group ("MonoBass: mono-sum invariant (Mid untouched, incl. through fades)");
    {
        const int N = 6000; std::vector<float> L0, R0; rngPair (N, 7, L0, R0); auto L = L0, R = R0;
        stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.5f); mb.reset();
        float* a[2] { L.data(), R.data() };                   mb.process (a, 2, 2000);   // settled w=0.5
        mb.setLowWidth (1.0f);
        float* b[2] { L.data() + 2000, R.data() + 2000 };     mb.process (b, 2, 500);    // mid-fade toward full-wide
        mb.setLowWidth (0.2f);
        float* c[2] { L.data() + 2500, R.data() + 2500 };     mb.process (c, 2, N - 2500);   // retargeted mid-fade
        double worst = 0.0;
        for (int i = 0; i < N; ++i) worst = std::max (worst, (double) std::fabs (0.5f * (L[i] + R[i]) - 0.5f * (L0[i] + R0[i])));
        test::ok (worst < 1e-6, "mono fold preserved through settled + fading states (bound: decode halves round separately)");
    }

    // --- interchannel correlation -> +1 below fc (mixed mono bass + anti-phase low stereo) ---
    test::group ("MonoBass: correlation -> +1 below fc");
    {
        const int N = 96000; std::vector<float> L (N), R (N);
        for (int i = 0; i < N; ++i)
        {
            const float bass = 0.5f * (float) std::sin (2.0 * pi * 60.0 * i / sr);
            const float nse  = 0.3f * (float) std::sin (2.0 * pi * 45.0 * i / sr);
            L[i] = bass + nse; R[i] = bass - nse;
        }
        auto corrTail = [&] (const std::vector<float>& A, const std::vector<float>& B)
        {
            double ab = 0, aa = 0, bb = 0;
            for (int i = N / 2; i < N; ++i) { ab += (double) A[i] * B[i]; aa += (double) A[i] * A[i]; bb += (double) B[i] * B[i]; }
            return ab / std::sqrt (aa * bb);
        };
        const double corrIn = corrTail (L, R);
        stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.0f); mb.reset();
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
        const double corrOut = corrTail (L, R);
        test::ok (corrIn < 0.6, "input decorrelated below fc (corr ~ 0.47)");
        test::ok (corrOut > 0.995, "output correlation -> +1 below fc (low side collapsed)");
    }

    // --- bypass paths are bit-exact; gating is STRICTLY stereo ---
    test::group ("MonoBass: bypass + strict numChannels==2 gating");
    {
        const int N = 512; std::vector<float> L0, R0; rngPair (N, 99, L0, R0);
        auto md2 = [&] (const std::vector<float>& A, const std::vector<float>& B, const std::vector<float>& A0, const std::vector<float>& B0)
        {
            double m = 0; for (int i = 0; i < N; ++i) { m = std::max (m, (double) std::fabs (A[i] - A0[i])); m = std::max (m, (double) std::fabs (B[i] - B0[i])); }
            return m;
        };
        { auto L = L0, R = R0; stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (1.0f); mb.reset();
          float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
          test::ok (md2 (L, R, L0, R0) == 0.0, "settled lowWidth=1 -> exact passthrough"); }
        { auto L = L0, R = R0; stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.0f); mb.reset(); mb.setEnabled (false);
          float* io[2] { L.data(), R.data() }; mb.process (io, 2, N);
          test::ok (md2 (L, R, L0, R0) == 0.0, "disabled -> exact passthrough"); }
        { auto L = L0; stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.0f); mb.reset();
          float* io[1] { L.data() }; mb.process (io, 1, N);
          double m = 0; for (int i = 0; i < N; ++i) m = std::max (m, (double) std::fabs (L[i] - L0[i]));
          test::ok (m == 0.0, "mono call (numChannels=1) -> untouched"); }
        {   // 5.1-style call: ALL channels must pass through — the front pair is NOT silently treated
            std::vector<std::vector<float>> ch (6), ref (6);
            for (int c = 0; c < 6; ++c) { std::vector<float> d; rngPair (N, 200 + (unsigned long long) c, ch[c], d); ref[c] = ch[c]; }
            float* io[6]; for (int c = 0; c < 6; ++c) io[c] = ch[c].data();
            stereo::MonoBass mb; mb.prepare (sr); mb.setLowWidth (0.0f); mb.reset();
            mb.process (io, 6, N);
            double m = 0; for (int c = 0; c < 6; ++c) for (int i = 0; i < N; ++i) m = std::max (m, (double) std::fabs (ch[c][i] - ref[c][i]));
            test::ok (m == 0.0, "surround call (numChannels=6) -> ALL channels untouched (strict stereo gate)"); }
    }

    // --- the bypass boundary: fading to/from full-wide is click-free, and re-entry has NO stale tails ---
    test::group ("MonoBass: bypass boundary click-free + no stale tails");
    {
        // Worst case on purpose: pure Side tone AT fc, where the wet allpass is -1 vs the raw +1 —
        // an unsmoothed switch would step by up to 2A; the crossfade must keep steps at carrier-slew
        // scale. fc is deliberately NOT an integer divisor of the segment (121.7 Hz + a phase offset):
        // 120 Hz over 4800-sample calls puts every call boundary on a zero crossing, where a hard
        // switch would be invisible — the topology change happens at call boundaries, so they must
        // land at non-zero carrier phase for this test to have teeth.
        const int seg = 4800, N = 6 * seg; const float A = 0.5f; const double fTone = 121.7;
        std::vector<float> L (N), R (N);
        for (int i = 0; i < N; ++i) { const float v = A * (float) std::sin (2.0 * pi * fTone * i / sr + 0.7); L[i] = v; R[i] = -v; }
        double inSlew = 0; for (int i = 1; i < N; ++i) inSlew = std::max (inSlew, (double) std::fabs (0.5f * (L[i] - R[i]) - 0.5f * (L[i - 1] - R[i - 1])));

        std::vector<float> Li = L, Ri = R;                                    // keep the input for the bypassed-segment check
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency ((float) fTone); mb.setLowWidth (0.0f); mb.reset();
        auto blk = [&] (int k) { float* io[2] { L.data() + k * seg, R.data() + k * seg }; mb.process (io, 2, seg); };
        blk (0);                                                              // settled mono-making
        mb.setLowWidth (1.0f); blk (1);                                       // 20 ms fade out, settles inside
        blk (2);                                                              // fully bypassed segment
        mb.setLowWidth (0.0f); blk (3); blk (4); blk (5);                     // fade back in from a CLEAN state
        double worstStep = 0; bool finite = true;
        for (int i = 1; i < N; ++i)
        {
            const double sNow = 0.5f * (L[i] - R[i]), sPrev = 0.5f * (L[i - 1] - R[i - 1]);
            worstStep = std::max (worstStep, std::fabs (sNow - sPrev));
            finite = finite && std::isfinite (L[i]) && std::isfinite (R[i]);
        }
        test::ok (finite, "finite throughout the fade out / bypass / fade in cycle");
        test::ok (worstStep < 3.0 * inSlew, "side steps stay at carrier-slew scale (no 2A hard-switch click at the boundary)");
        double mdByp = 0; for (int i = 2 * seg; i < 3 * seg; ++i) { mdByp = std::max (mdByp, (double) std::fabs (L[i] - Li[i])); mdByp = std::max (mdByp, (double) std::fabs (R[i] - Ri[i])); }
        test::ok (mdByp == 0.0, "the settled full-wide segment is a bit-exact passthrough");

        // Stale-tail proof: hammer a LOUD low tone into bypass, then come back on SILENCE — any old filter
        // state would ring out audibly; the reset-on-entering-bypass must make the return exactly silent.
        stereo::MonoBass mb2; mb2.prepare (sr); mb2.setFrequency (120.0f); mb2.setLowWidth (0.0f); mb2.reset();
        std::vector<float> l (seg), r (seg); float* io2[2] { l.data(), r.data() };
        auto loud = [&] { for (int i = 0; i < seg; ++i) { const float v = 0.9f * (float) std::sin (2.0 * pi * 30.0 * i / sr); l[i] = v; r[i] = -v; } };
        loud(); mb2.process (io2, 2, seg);                                    // warm, loud low side
        mb2.setLowWidth (1.0f); loud(); mb2.process (io2, 2, seg);            // fade out on loud material
        loud(); mb2.process (io2, 2, seg);                                    // settled bypass -> state reset here
        mb2.setLowWidth (0.0f);
        std::fill (l.begin(), l.end(), 0.0f); std::fill (r.begin(), r.end(), 0.0f);
        mb2.process (io2, 2, seg);                                            // re-enter on silence
        double tail = 0; for (int i = 0; i < seg; ++i) { tail = std::max (tail, (double) std::fabs (l[i])); tail = std::max (tail, (double) std::fabs (r[i])); }
        test::ok (tail == 0.0, "re-entry from bypass on silence is EXACTLY silent (no stale filter tails)");
    }

    // --- a live lowWidth step inside (0,1) is smoothed (no click), and lands on the analytic blend ---
    test::group ("MonoBass: lowWidth step is click-free and lands on the blend law");
    {
        const int N = 96000; std::vector<float> L (N), R (N);
        for (int i = 0; i < N; ++i) { const float v = 0.5f * (float) std::sin (2.0 * pi * 60.0 * i / sr); L[i] = v; R[i] = -v; }
        double inSlew = 0; for (int i = 1; i < N; ++i) inSlew = std::max (inSlew, (double) std::fabs (0.5f * (L[i] - R[i]) - 0.5f * (L[i - 1] - R[i - 1])));
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (120.0f); mb.setLowWidth (0.0f); mb.reset();
        float* a[2] { L.data(), R.data() };                   mb.process (a, 2, N / 2);
        mb.setLowWidth (0.8f);                                                // STEP — must ramp, not jump
        float* b[2] { L.data() + N / 2, R.data() + N / 2 };   mb.process (b, 2, N / 2);
        double worstStep = 0;
        for (int i = 1; i < N; ++i) worstStep = std::max (worstStep, std::fabs ((double) (0.5f * (L[i] - R[i])) - (double) (0.5f * (L[i - 1] - R[i - 1]))));
        test::ok (worstStep < 3.0 * inSlew, "side steps stay at carrier-slew scale through the w ramp (unsmoothed jump would be ~0.3)");
        std::vector<float> side (N); for (int i = 0; i < N; ++i) side[i] = 0.5f * (L[i] - R[i]);
        const double got = toneAmp (side, N - 24000, 60.0, sr) / 0.5, want = analyticSide (60.0, 120.0, sr, 0.8);
        test::approx (got, want, 2e-5 + 0.005 * want, "post-ramp steady state == the w=0.8 blend law");
    }

    // --- parameter abuse: near-Nyquist clamp, low clamp, prepare re-clamp, non-finite rejection ---
    test::group ("MonoBass: fc / param abuse stays sane");
    {
        stereo::MonoBass mb; mb.prepare (sr);
        mb.setFrequency (1.0e9f);  test::ok (mb.frequency() == (float) (0.45 * sr), "fc clamps to 0.45*fs");
        mb.setFrequency (1.0f);    test::ok (mb.frequency() == 20.0f, "fc clamps up to 20 Hz");
        mb.setFrequency (std::nanf ("")); mb.setFrequency (INFINITY); mb.setFrequency (-INFINITY);
        test::ok (mb.frequency() == 20.0f, "setFrequency(NaN/inf) ignored — last good value kept");
        mb.setLowWidth (0.6f); mb.setLowWidth (std::nanf ("")); mb.setLowWidth (-INFINITY);
        test::ok (mb.lowWidth() == 0.6f, "setLowWidth(NaN/-inf) ignored — last good value kept");

        mb.setFrequency (5000.0f); mb.prepare (8000.0);
        test::ok (mb.frequency() == 3600.0f, "prepare() at a lower fs re-clamps fc to the new 0.45*fs");

        // prepare() abuse: non-finite / absurdly low rates must not reach the clamp with lo > hi (UB)
        // or the smoother ramp with an inf sample count.
        { stereo::MonoBass m4; m4.prepare (INFINITY);
          test::ok (m4.frequency() == 120.0f, "prepare(inf) -> 48 kHz fallback, fc clamp sane"); }
        { stereo::MonoBass m4; m4.prepare (std::nan ("")); m4.setFrequency (200.0f);
          test::ok (m4.frequency() == 200.0f, "prepare(NaN) -> 48 kHz fallback, setters live"); }
        { stereo::MonoBass m4; m4.prepare (30.0);                            // 0.45*fs < kMinFreq: hi must not sink below lo
          test::ok (m4.frequency() == 20.0f, "prepare(30 Hz) -> fc pinned at the 20 Hz floor (no lo>hi clamp UB)");
          std::vector<float> l (256, 0.5f), r (256, -0.5f); float* io4[2] { l.data(), r.data() };
          m4.process (io4, 2, 256);
          bool fin = true; for (int i = 0; i < 256; ++i) fin = fin && std::isfinite (l[i]) && std::isfinite (r[i]);
          test::ok (fin, "processing at an absurd fs stays finite"); }

        // Hammer noise through the extreme corners — the output must stay finite everywhere.
        bool finite = true;
        for (float fc : { 20.0f, (float) (0.45 * sr) })
        {
            stereo::MonoBass m2; m2.prepare (sr); m2.setFrequency (fc); m2.setLowWidth (0.0f); m2.reset();
            std::vector<float> L, R; rngPair (48000, 3 + (unsigned long long) fc, L, R);
            float* io[2] { L.data(), R.data() };
            for (int k = 0; k < 48000 / 512; ++k) { float* b[2] { L.data() + k * 512, R.data() + k * 512 }; m2.process (b, 2, 512); }
            (void) io;
            for (int i = 0; i < 48000; ++i) finite = finite && std::isfinite (L[i]) && std::isfinite (R[i]);
        }
        stereo::MonoBass m3; m3.prepare (0.0);                               // absurd fs -> 48 kHz fallback
        std::vector<float> L, R; rngPair (2048, 5, L, R); float* io[2] { L.data(), R.data() };
        m3.process (io, 2, 2048);
        for (int i = 0; i < 2048; ++i) finite = finite && std::isfinite (L[i]) && std::isfinite (R[i]);
        test::ok (finite, "fc at both clamp rails + prepare(0) fallback: output finite under noise");
    }

    // --- loud -> silence: the state must flush to EXACT zero (no denormal-sustained tails) ---
    test::group ("MonoBass: denormal flush -> exact zero after silence");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (120.0f); mb.setLowWidth (0.0f); mb.reset();
        std::vector<float> L (4800), R (4800);
        for (int i = 0; i < 4800; ++i) { const float v = 0.9f * (float) std::sin (2.0 * pi * 120.0 * i / sr); L[i] = v; R[i] = -v; }
        float* io[2] { L.data(), R.data() }; mb.process (io, 2, 4800);
        std::vector<float> zl (512), zr (512); float* zio[2] { zl.data(), zr.data() };
        bool finite = true;
        for (int k = 0; k < 94; ++k)                                          // ~1 s of true silence, re-zeroed per block
        {
            std::fill (zl.begin(), zl.end(), 0.0f); std::fill (zr.begin(), zr.end(), 0.0f);
            mb.process (zio, 2, 512);
            for (int i = 0; i < 512; ++i) finite = finite && std::isfinite (zl[i]) && std::isfinite (zr[i]);
        }
        double lastMax = 0; for (int i = 0; i < 512; ++i) { lastMax = std::max (lastMax, (double) std::fabs (zl[i])); lastMax = std::max (lastMax, (double) std::fabs (zr[i])); }
        test::ok (finite, "decay stays finite");
        test::ok (lastMax == 0.0, "state flushed -> output EXACTLY zero after 1 s of silence (no denormal hum)");
    }

    // --- no allocation in process(), including through param moves and a bypass fade ---
    test::group ("MonoBass: no-alloc");
    {
        stereo::MonoBass mb; mb.prepare (sr); mb.setFrequency (120.0f); mb.setLowWidth (0.0f); mb.reset();
        std::vector<float> L (512, 0.2f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        mb.process (io, 2, 512);
        const long before = g_allocs.load();
        mb.process (io, 2, 512);
        mb.setLowWidth (1.0f);  mb.process (io, 2, 512);                     // fade + settle into bypass
        mb.process (io, 2, 512);                                             // bypassed (state reset path)
        mb.setLowWidth (0.3f); mb.setFrequency (90.0f); mb.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process()/setters did not allocate (incl. bypass transitions)");
    }

    //==========================================================================
    // StereoWidth: a Mid/Side side-gain. The headline safety property is the mono-fold invariant.

    // --- THE mono-fold invariant: ½(L'+R') == ½(L+R) for ANY width (we touch only the Side) ---
    test::group ("StereoWidth: mono-fold invariant");
    {
        const int N = 4000; std::vector<float> L0, R0; rngPair (N, 7, L0, R0);
        double worst = 0.0;
        for (float w : { 0.0f, 0.5f, 1.3f, 2.0f })
        {
            std::vector<float> L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (w);
            float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
            for (int i = 0; i < N; ++i) worst = std::max (worst, (double) std::fabs (0.5f * (L[i] + R[i]) - 0.5f * (L0[i] + R0[i])));
        }
        test::ok (worst < 1e-6, "mono sum unchanged at width ∈ {0, .5, 1.3, 2} (widening can't weaken the fold)");
    }

    // --- neutral (width=1, gain=1) is BIT-EXACT passthrough (the early-return bypass) ---
    test::group ("StereoWidth: neutral is bit-exact");
    {
        const int N = 512; std::vector<float> L0, R0; rngPair (N, 3, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr);                              // defaults: width 1, gain 1
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "width=1, gain=1 → exact passthrough");
    }

    // --- width=0 collapses to the mono sum (L'==R'==M), exactly ---
    test::group ("StereoWidth: width=0 → mono");
    {
        const int N = 1000; std::vector<float> L0, R0; rngPair (N, 9, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (0.0f); sw.reset();   // snap past the smoothing ramp
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double diff = 0, toM = 0; for (int i = 0; i < N; ++i) { diff = std::max (diff, (double) std::fabs (L[i] - R[i])); toM = std::max (toM, (double) std::fabs (L[i] - 0.5f * (L0[i] + R0[i]))); }
        test::ok (diff == 0.0 && toM < 1e-6, "width=0 → L==R==½(L+R) (pure mono)");
    }

    // --- width=2 doubles the side: (L'-R') == 2·(L-R) ---
    test::group ("StereoWidth: width=2 doubles the side");
    {
        const int N = 1000; std::vector<float> L0, R0; rngPair (N, 21, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (2.0f); sw.reset();   // snap past the smoothing ramp
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs ((L[i] - R[i]) - 2.0f * (L0[i] - R0[i])));
        test::ok (md < 1e-5, "side difference scaled ×2 (the S axis tracks width)");
    }

    // --- a MONO source stays bit-exact for any width (you cannot synthesise stereo from S=0) ---
    test::group ("StereoWidth: mono source is invariant");
    {
        const int N = 600; std::vector<float> in; std::vector<float> dummy; rngPair (N, 5, in, dummy);
        auto L = in, R = in;                                                  // L==R → S=0
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (2.0f);
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - in[i])); md = std::max (md, (double) std::fabs (R[i] - in[i])); }
        test::ok (md == 0.0, "L==R in → unchanged out at width=2 (width scales zero side)");
    }

    // --- outputGain is a clean linear trim on the whole signal ---
    test::group ("StereoWidth: outputGain trims level");
    {
        const int N = 600; std::vector<float> L0, R0; rngPair (N, 13, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.0f); sw.setOutputGain (2.0f); sw.reset();
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - 2.0f * L0[i])); md = std::max (md, (double) std::fabs (R[i] - 2.0f * R0[i])); }
        test::ok (md < 1e-6, "width=1, gain=2 → output = 2× input (M and S both scale)");
    }

    // --- side energy tracks width monotonically (and ≈ proportional) ---
    test::group ("StereoWidth: side energy ∝ width");
    {
        const int N = 4000; std::vector<float> L0, R0; rngPair (N, 31, L0, R0);
        auto sideRms = [&] (float w) {
            std::vector<float> L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (w); sw.reset();
            float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
            std::vector<float> sde (N); for (int i = 0; i < N; ++i) sde[i] = 0.5f * (L[i] - R[i]);
            return rmsTail (sde, 0);
        };
        const double s05 = sideRms (0.5f), s10 = sideRms (1.0f), s15 = sideRms (1.5f);
        test::ok (s05 < s10 && s10 < s15 && std::fabs (s15 / s10 - 1.5) < 0.02, "narrower→less / wider→more side, ratio ≈ width");
    }

    // --- bypass paths: < 2 channels and disabled are untouched ---
    test::group ("StereoWidth: bypass");
    {
        const int N = 256; std::vector<float> L0, R0; rngPair (N, 99, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.8f);
        float* mono[1] { L.data() }; sw.process (mono, 1, N);                 // < 2 ch → passthrough
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (L[i] - L0[i]));
        test::ok (md == 0.0, "mono call (numChannels<2) → untouched");
        sw.setEnabled (false); float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        md = 0; for (int i = 0; i < N; ++i) { md = std::max (md, (double) std::fabs (L[i] - L0[i])); md = std::max (md, (double) std::fabs (R[i] - R0[i])); }
        test::ok (md == 0.0, "disabled → untouched");
    }

    // --- no allocation in process() ---
    test::group ("StereoWidth: no-alloc");
    {
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.7f); sw.setOutputGain (1.1f);
        std::vector<float> L (512, 0.3f), R (512, -0.2f); float* io[2] { L.data(), R.data() };
        sw.process (io, 2, 512);
        const long before = g_allocs.load();
        sw.process (io, 2, 512); sw.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate");
    }

    // --- a live width step is SMOOTHED (no click) — and the mono fold stays invariant through the ramp ---
    test::group ("StereoWidth: width step is click-free");
    {
        const int N = 4000; std::vector<float> L0 (N), R0 (N);
        for (int i = 0; i < N; ++i) { L0[i] = 0.5f * (float) std::sin (2.0 * pi * 500.0 * i / sr); R0[i] = 0.5f * (float) std::sin (2.0 * pi * 500.0 * i / sr + 0.6); }
        double inSlew = 0; for (int i = 1; i < N; ++i) inSlew = std::max (inSlew, (double) std::fabs (0.5f * (L0[i] - R0[i]) - 0.5f * (L0[i - 1] - R0[i - 1])));
        auto L = L0, R = R0; stereo::StereoWidth sw; sw.prepare (sr); sw.reset();        // settled at width 1
        float* a[2] { L.data(), R.data() }; sw.process (a, 2, N / 2);                     // width 1 (near-passthrough)
        sw.setWidth (2.0f);                                                              // STEP — must ramp, not jump
        float* b[2] { L.data() + N / 2, R.data() + N / 2 }; sw.process (b, 2, N / 2);
        double worstStep = 0, worstMono = 0;
        for (int i = 1; i < N; ++i)
        {
            worstStep = std::max (worstStep, (double) std::fabs (0.5f * (L[i] - R[i]) - 0.5f * (L[i - 1] - R[i - 1])));
            worstMono = std::max (worstMono, (double) std::fabs (0.5f * (L[i] + R[i]) - 0.5f * (L0[i] + R0[i])));
        }
        test::ok (worstStep < 4.0 * inSlew, "width 1→2 ramps smoothly (side step ≪ the instant 2× jump an unsmoothed step gives)");
        test::ok (worstMono < 1e-6, "mono fold stays invariant THROUGH the width ramp (Mid untouched per sample)");
    }

    // --- a stray non-finite parameter is rejected, never poisons the stream ---
    test::group ("StereoWidth: non-finite params rejected");
    {
        const int N = 256; std::vector<float> L0, R0; rngPair (N, 55, L0, R0); auto L = L0, R = R0;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.5f); sw.setOutputGain (1.2f); sw.reset();
        sw.setWidth (std::nanf ("")); sw.setOutputGain (INFINITY); sw.setWidth (-INFINITY);   // all must be ignored
        test::ok (sw.width() == 1.5f && sw.outputGain() == 1.2f, "setWidth(NaN)/setOutputGain(±inf) ignored — last good value kept");
        float* io[2] { L.data(), R.data() }; sw.process (io, 2, N);
        bool fin = true; for (int i = 0; i < N; ++i) if (! std::isfinite (L[i]) || ! std::isfinite (R[i])) fin = false;
        test::ok (fin, "output stays finite after a NaN/inf parameter poke");
    }

    // --- aliased channels (io[0] == io[1], dual-mono) apply gain ONCE, not squared ---
    test::group ("StereoWidth: aliased L==R applies gain once");
    {
        const int N = 200; std::vector<float> in, dummy; rngPair (N, 77, in, dummy); auto buf = in;
        stereo::StereoWidth sw; sw.prepare (sr); sw.setWidth (1.7f); sw.setOutputGain (2.0f); sw.reset();
        float* io[2] { buf.data(), buf.data() };                                  // L and R alias one buffer
        sw.process (io, 2, N);
        double md = 0; for (int i = 0; i < N; ++i) md = std::max (md, (double) std::fabs (buf[i] - 2.0f * in[i]));
        test::ok (md < 1e-5, "io[0]==io[1] → output = gain·in (gain applied once; width inert on S=0)");
    }

    return test::report();
}
