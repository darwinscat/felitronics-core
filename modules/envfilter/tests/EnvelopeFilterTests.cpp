// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// envfilter::EnvelopeFilter and its one-slider law. The effect: the cutoff tracks the envelope (measured from the
// audio with a probe tone, at two levels, both directions, three rates); stability at maximum resonance under the
// fastest modulation and on a full-scale slap; no zipper in a swept output; the clean lows hold the low end. The
// stage: the law-11 properties, amount 0 = bit-exact dry, a slider move glides without a step, re-slicing (8a),
// a stopped channel and a pause (11a/11c), reset is a restart, poison, law 8, mono = stereo, no allocation.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter
#include <law11_call_contract.h>
#include <law11c_pause.h>

#include <felitronics/envfilter/EnvelopeFilter.h>
#include <felitronics/envfilter/EnvelopeFilterLaw.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace core = felitronics::core;   // ADAPT expands core::kMaxChannels here
namespace ef   = felitronics::envfilter;
namespace eq   = felitronics::eq;
using felitronics::test::ok;
using felitronics::test::group;
using felitronics::test::run;

namespace
{
using Stage  = ef::EnvelopeFilter;
using Params = ef::EnvelopeFilterParams;
using Buf    = std::vector<std::vector<float>>;

constexpr double kPi = 3.14159265358979323846;
const double kRates[] = { 44100.0, 48000.0, 96000.0 };

std::string fmt (double v, int prec = 2)
{
    char b[64];
    std::snprintf (b, sizeof (b), "%.*f", prec, v);
    return b;
}

/** A check that also prints its measurement: the numbers are the point of these tests, so they are on the record
    whether or not the bound holds. */
void okv (bool cond, const std::string& msg)
{
    ok (cond, msg);
    std::printf ("      %s\n", msg.c_str());
}

double dbOf (double a) { return 20.0 * std::log10 (std::max (a, 1e-300)); }

/** Process a whole buffer in `block`-sized calls; false if any call was refused. */
bool render (Stage& s, Buf& io, int block)
{
    const int nch = (int) io.size();
    const int n   = (int) io[0].size();
    bool all = true;
    for (int at = 0; at < n; at += block)
    {
        float* p[2] = { io[0].data() + at, nch > 1 ? io[1].data() + at : nullptr };
        all = run (s.process (p, nch, std::min (block, n - at))) && all;
    }
    return all;
}

/** A tone with Hann-shaped gates of `edge` seconds at both ends (the fixture rule: taper both ends). */
void addBurst (std::vector<float>& x, double fs, double t0, double len, double edge, double hz, double amp)
{
    const int i0 = (int) (t0 * fs), n = (int) (len * fs), e = std::max (1, (int) (edge * fs));
    for (int i = 0; i < n && i0 + i < (int) x.size(); ++i)
    {
        double g = 1.0;
        if (i < e)          g = 0.5 - 0.5 * std::cos (kPi * i / e);
        else if (i > n - e) g = 0.5 - 0.5 * std::cos (kPi * (n - i) / e);
        x[(std::size_t) (i0 + i)] += (float) (amp * g * std::sin (2.0 * kPi * hz * i / fs));
    }
}

/** The amplitude of the `hz` component in x[a, b), Hann-windowed (sidelobes fall 18 dB/oct). */
double toneAmp (const std::vector<float>& x, int a, int b, double hz, double fs)
{
    double re = 0.0, im = 0.0, wsum = 0.0;
    const int n = b - a;
    for (int i = 0; i < n; ++i)
    {
        const double w = 0.5 - 0.5 * std::cos (2.0 * kPi * (i + 0.5) / n);
        const double v = w * (double) x[(std::size_t) (a + i)];
        re += v * std::cos (2.0 * kPi * hz * (a + i) / fs);
        im -= v * std::sin (2.0 * kPi * hz * (a + i) / fs);
        wsum += w;
    }
    return 2.0 * std::sqrt (re * re + im * im) / wsum;
}

/** The TPT low pass's exact digital magnitude at `hz` for a cutoff `fc` (the bilinear map of 1/(s² + ks + 1)). */
double lpMag (double hz, double fc, double q, double fs)
{
    const double w = std::tan (kPi * hz / fs) / std::tan (kPi * fc / fs);
    const double re = 1.0 - w * w, im = w / q;
    return 1.0 / std::sqrt (re * re + im * im);
}

double maxStep (const std::vector<float>& y, int a, int b)
{
    double m = 0.0;
    for (int i = std::max (1, a); i < b && i < (int) y.size(); ++i)
        m = std::max (m, (double) std::fabs (y[(std::size_t) i] - y[(std::size_t) (i - 1)]));
    return m;
}

Params voice (ef::Mode m = ef::Mode::lowPass, ef::Direction d = ef::Direction::up, float q = 4.0f)
{
    Params p;
    p.mode = m; p.direction = d; p.resonance = q;
    return p;
}

//==============================================================================
void sliderLaw()
{
    group ("the one-slider law: amount -> mix and sensitivity, the menu passes through");
    const auto p0 = ef::paramsForAmount (0.0f), pq = ef::paramsForAmount (0.25f), ph = ef::paramsForAmount (0.5f),
               p1 = ef::paramsForAmount (1.0f);
    ok (p0.mix <= 0.0f && pq.mix == 0.5f && ph.mix == 1.0f && p1.mix == 1.0f, "mix 0, 0.5, 1, 1 at 0, 1/4, 1/2, 1");
    ok (p0.sensitivityDb == -6.0f && pq.sensitivityDb == 0.0f && ph.sensitivityDb == 6.0f && p1.sensitivityDb == 18.0f,
        "sensitivity -6, 0, +6, +18 dB at 0, 1/4, 1/2, 1");
    bool mono = true;
    for (int i = 1; i <= 100; ++i)
    {
        const auto a = ef::paramsForAmount ((float) (i - 1) / 100.0f), b = ef::paramsForAmount ((float) i / 100.0f);
        mono = mono && b.mix >= a.mix && b.sensitivityDb > a.sensitivityDb;
    }
    ok (mono, "monotone over the travel: more slider is never less wah");
    Params menu = voice (ef::Mode::bandPass, ef::Direction::down, 7.0f);
    menu.keepLows = 0.3f; menu.sweepLoHz = 200.0f; menu.releaseMs = 300.0f;
    const auto pm = ef::paramsForAmount (0.6f, menu);
    ok (pm.mode == ef::Mode::bandPass && pm.direction == ef::Direction::down && pm.resonance == 7.0f
        && pm.keepLows == 0.3f && pm.sweepLoHz == 200.0f && pm.releaseMs == 300.0f, "the voice is the menu's, untouched");
    ok (ef::paramsForAmount (std::numeric_limits<float>::quiet_NaN()).mix <= 0.0f && ef::paramsForAmount (-1.0f).mix <= 0.0f
        && ef::paramsForAmount (2.0f).sensitivityDb == 18.0f, "a NaN or below-range amount is 0, above-range is 1");
}

//==============================================================================
using felitronics::test::law11::kFs;
using felitronics::test::law11::kMaxBlock;
using felitronics::test::law11::kPrepCh;

ADAPT (A_EnvelopeFilter, felitronics::envfilter::EnvelopeFilter, kPrepCh, true, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh);
         auto p = felitronics::envfilter::paramsForAmount (1.0f);           // a moving, resonant band pass, so
         p.mode = felitronics::envfilter::Mode::bandPass; p.resonance = 8.0f; // the short fixture is far from dry
         p.keepLows = 0.5f; p.attackMs = 0.5f; p.releaseMs = 20.0f;
         s.setParams (p); return okp; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

void contract()
{
    group ("law 11: the call contract, swept (P1-P8), and the stage's own prepare refusals");
    felitronics::test::law11::allProperties<A_EnvelopeFilter>();

    Stage s;
    ok (! s.prepare (7999.0, 64, 2) && ! s.prepare (768001.0, 64, 2) && ! s.prepare (std::nan (""), 64, 2)
        && ! s.prepare (48000.0, 0, 2) && ! s.prepare (48000.0, 64, 0) && ! s.prepare (48000.0, 64, 3),
        "prepare refuses a rate outside 8 k ... 768 k, a NaN rate, maxBlock 0, and 0 or 3 channels");
    ok (s.prepare (8000.0, 1, 1) && s.prepare (768000.0, 8192, 2), "...and accepts both ends");
    ok (! s.prepare (48000.0, 64, 3), "precondition: a refused re-prepare");
    std::vector<float> a (64, 0.5f);
    float* p[1] = { a.data() };
    ok (! s.process (p, 1, 64) && a[0] == 0.5f, "after it the stage refuses and touches nothing (law 11b)");
    ok (Stage::latencySamples() == 0, "zero latency");
}

//==============================================================================
void dryAtZero()
{
    group ("amount 0: bit-exact dry, NaN payloads included, and nothing computed — fresh, and after a move to 0");
    for (double fs : kRates)
        for (int nch = 1; nch <= 2; ++nch)
        {
            Stage s;
            run (s.prepare (fs, 256, 2));
            s.setParams (ef::paramsForAmount (0.0f, voice (ef::Mode::bandPass)));
            Buf io (2, std::vector<float> ((std::size_t) fs, 0.0f));
            felitronics::test::law11::Rng r { 0x5EEDu };
            for (auto& c : io) for (auto& v : c) v = r.next();
            std::uint32_t nanBits = 0x7fc12345u;
            std::memcpy (&io[0][100], &nanBits, 4);
            io[0][200] = std::numeric_limits<float>::infinity();
            io.resize ((std::size_t) nch);
            const Buf in = io;
            render (s, io, 256);
            okv (std::memcmp (io[0].data(), in[0].data(), io[0].size() * 4) == 0
                && (nch < 2 || std::memcmp (io[1].data(), in[1].data(), io[1].size() * 4) == 0)
                && s.isAsleep(), "fresh at amount 0: every bit of the input comes back, asleep — "
                                 + fmt (fs, 0) + " Hz, " + std::to_string (nch) + " ch");
        }

    // A move to 0 glides out and then sleeps: from the landing on, the output is the input bit for bit.
    for (double fs : kRates)
    {
        Stage s;
        run (s.prepare (fs, 128, 1));
        s.setParams (ef::paramsForAmount (0.8f, voice (ef::Mode::bandPass, ef::Direction::up, 8.0f)));
        const int n = (int) (2.0 * fs);
        Buf io (1, std::vector<float> ((std::size_t) n, 0.0f));
        addBurst (io[0], fs, 0.0, 2.0, 0.01, 98.0, 0.3);
        const Buf in = io;
        const int at = (int) (0.5 * fs) / 128 * 128;
        bool all = true;
        for (int i = 0; i < n; i += 128)
        {
            if (i == at) s.setParams (ef::paramsForAmount (0.0f, voice (ef::Mode::bandPass, ef::Direction::up, 8.0f)));
            float* p[1] = { io[0].data() + i };
            all = run (s.process (p, 1, std::min (128, n - i))) && all;
        }
        int lastDiff = -1;
        for (int i = 0; i < n; ++i)
            if (! core::sameBits (io[0][(std::size_t) i], in[0][(std::size_t) i])) lastDiff = i;
        const double landMs = 1000.0 * (lastDiff - at) / fs;
        okv (all && s.isAsleep() && lastDiff > at && landMs < 500.0,
            "after a move to 0 the glide lands and the stage sleeps: dry bit for bit " + fmt (landMs, 0)
            + " ms after the move (bound 500) — " + fmt (fs, 0) + " Hz");
    }
}

//==============================================================================
// A slider move never steps: the largest sample-to-sample step of the output around each move stays within the
// largest step of the steady output just before or after it (the tone's own slope, through whatever the filter
// makes of it). The jumps are deliberately large: off -> 0.6, -> 1, -> 0.2, -> 0.
void sliderGlides()
{
    group ("a slider move glides: no step larger than the steady output's own");
    const double fs = 48000.0;
    const int    n  = (int) (2.6 * fs);
    const double hz = 82.41, amp = 0.5;
    struct Move { double t; float amount; };
    const Move moves[] = { { 0.3, 0.6f }, { 0.8, 1.0f }, { 1.3, 0.2f }, { 1.8, 0.0f } };

    for (auto mode : { ef::Mode::lowPass, ef::Mode::bandPass })
    {
        Stage s;
        run (s.prepare (fs, 64, 1));
        const Params menu = voice (mode, ef::Direction::up, 6.0f);
        s.setParams (ef::paramsForAmount (0.0f, menu));
        std::vector<float> y ((std::size_t) n, 0.0f);
        addBurst (y, fs, 0.0, 2.6, 0.05, hz, amp);
        std::size_t k = 0;
        for (int i = 0; i < n; i += 64)
        {
            if (k < std::size (moves) && i >= (int) (moves[k].t * fs)) s.setParams (ef::paramsForAmount (moves[k++].amount, menu));
            float* p[1] = { y.data() + i };
            run (s.process (p, 1, std::min (64, n - i)));
        }
        double worst = 0.0;
        for (const auto& m : moves)
        {
            const int t = (int) (m.t * fs) / 64 * 64;
            const double around = maxStep (y, t, t + (int) (0.15 * fs));
            const double steady = std::max (maxStep (y, t - (int) (0.1 * fs), t),
                                            maxStep (y, t + (int) (0.35 * fs), t + (int) (0.45 * fs)));
            worst = std::max (worst, around / steady);
        }
        okv (worst <= 1.25, std::string (mode == ef::Mode::lowPass ? "low pass" : "band pass")
                           + ": the largest step around a move is " + fmt (worst, 3) + "x the steady one (bound 1.25)");
    }

    // The menu's settings ride the same glides: a mode flip, a direction flip, the resonance, the range, keepLows.
    {
        Stage s;
        run (s.prepare (fs, 64, 1));
        Params p = ef::paramsForAmount (0.7f, voice (ef::Mode::lowPass, ef::Direction::up, 2.0f));
        s.setParams (p);
        std::vector<float> y ((std::size_t) n, 0.0f);
        addBurst (y, fs, 0.0, 2.6, 0.05, hz, amp);
        const double ts[] = { 0.3, 0.7, 1.1, 1.5, 1.9 };
        std::size_t k = 0;
        for (int i = 0; i < n; i += 64)
        {
            if (k < std::size (ts) && i >= (int) (ts[k] * fs))
            {
                switch (k)
                {
                    case 0:  p.mode = ef::Mode::bandPass; break;
                    case 1:  p.direction = ef::Direction::down; break;
                    case 2:  p.resonance = 12.0f; break;
                    case 3:  p.sweepLoHz = 400.0f; p.sweepHiHz = 1200.0f; break;
                    default: p.keepLows = 0.0f; break;
                }
                s.setParams (p);
                ++k;
            }
            float* q[1] = { y.data() + i };
            run (s.process (q, 1, std::min (64, n - i)));
        }
        double worst = 0.0;
        for (double tt : ts)
        {
            const int t = (int) (tt * fs) / 64 * 64;
            const double around = maxStep (y, t, t + (int) (0.15 * fs));
            const double steady = std::max (maxStep (y, t - (int) (0.1 * fs), t),
                                            maxStep (y, t + (int) (0.25 * fs), t + (int) (0.35 * fs)));
            worst = std::max (worst, around / steady);
        }
        okv (worst <= 1.25, "menu moves (mode, direction, Q 2 -> 12, range, keepLows): the largest step is "
                           + fmt (worst, 3) + "x the steady one (bound 1.25)");
    }
}

//==============================================================================
// THE CUTOFF TRACKS THE ENVELOPE. Two Hann-gated bursts of a 110 Hz tone, -24 and -12 dBFS, carry a probe at
// 1100 Hz and -70 dBFS (too quiet to move the detector: +0.04 dB on the softer burst; an even multiple of 110, so
// none of the ripple's sidebands — the odd multiples — lands on it). The low pass at Q 0.707 makes the probe's
// gain a measurement of the cutoff, and it is compared with the analytic response at the cutoff the stage
// publishes: the published sweep is the one in the audio.
void tracking()
{
    group ("the cutoff tracks the envelope: the louder hit opens further — measured in the audio");
    const double levelsDb[2] = { -24.0, -12.0 };
    const double sensDb = 0.0, q = 0.70710678;

    for (double fs : kRates)
        for (auto dir : { ef::Direction::up, ef::Direction::down })
        {
            Stage s;
            run (s.prepare (fs, 16, 1));
            Params p = voice (ef::Mode::lowPass, dir, (float) q);
            p.sensitivityDb = (float) sensDb; p.keepLows = 0.0f; p.mix = 1.0f;
            s.setParams (p);

            const double probeHz = 1100.0, probeAmp = std::pow (10.0, -70.0 / 20.0);
            const double t0[2] = { 0.2, 1.2 }, len = 0.4;
            const int n = (int) (2.4 * fs);
            std::vector<float> x ((std::size_t) n, 0.0f);
            for (int b = 0; b < 2; ++b)
            {
                addBurst (x, fs, t0[b], len, 0.015, 110.0, std::pow (10.0, levelsDb[b] / 20.0));
                addBurst (x, fs, t0[b], len, 0.015, probeHz, probeAmp);
            }
            std::vector<float> y = x;
            std::vector<double> cut ((std::size_t) (n / 16 + 1), 0.0);
            for (int i = 0; i < n; i += 16)
            {
                float* io[1] = { y.data() + i };
                run (s.process (io, 1, std::min (16, n - i)));
                cut[(std::size_t) (i / 16)] = (double) s.shownCutoffHz();
            }

            const std::string tag = fmt (fs, 0) + " Hz " + (dir == ef::Direction::up ? "up" : "down");
            double lfMean[2] {}, gainDb[2] {}, predDb[2] {};
            const double lfLo = std::log2 (150.0), lfHi = std::log2 (3000.0);
            bool lawOk = true;
            for (int b = 0; b < 2; ++b)
            {
                const int a = (int) ((t0[b] + 0.15) * fs), e = (int) ((t0[b] + 0.35) * fs);
                double acc = 0.0; int cnt = 0;
                for (int i = a / 16; i < e / 16; ++i) { acc += std::log2 (cut[(std::size_t) i]); ++cnt; }
                lfMean[b] = acc / cnt;
                const double pos = 1.0 + (levelsDb[b] + sensDb) / Stage::kSpanDb;
                const double want = dir == ef::Direction::up ? lfLo + pos * (lfHi - lfLo) : lfHi - pos * (lfHi - lfLo);
                lawOk = lawOk && std::fabs (lfMean[b] - want) <= 0.1;
                gainDb[b] = dbOf (toneAmp (y, a, e, probeHz, fs) / toneAmp (x, a, e, probeHz, fs));
                predDb[b] = dbOf (lpMag (probeHz, std::exp2 (lfMean[b]), q, fs));
            }
            const double sweep = std::fabs (lfMean[1] - lfMean[0]);
            const double want  = 12.0 / Stage::kSpanDb * (lfHi - lfLo);
            okv (dir == ef::Direction::up ? lfMean[1] > lfMean[0] : lfMean[1] < lfMean[0],
                tag + ": the louder burst moves the cutoff further from rest ("
                + fmt (std::exp2 (lfMean[0]), 0) + " -> " + fmt (std::exp2 (lfMean[1]), 0) + " Hz)");
            okv (std::fabs (sweep - want) <= 0.1, tag + ": 12 dB more moves it " + fmt (sweep, 3) + " oct (the law: "
                                                 + fmt (want, 3) + ", tol 0.1)");
            ok (lawOk, tag + ": both bursts sit where the law puts them, within 0.1 oct");
            okv (std::fabs (gainDb[0] - predDb[0]) <= 0.5 && std::fabs (gainDb[1] - predDb[1]) <= 0.5,
                tag + ": the probe's measured gain matches the published cutoff's response ("
                + fmt (gainDb[0]) + " vs " + fmt (predDb[0]) + ", " + fmt (gainDb[1]) + " vs " + fmt (predDb[1]) + " dB, tol 0.5)");
            const double open = dir == ef::Direction::up ? gainDb[1] - gainDb[0] : gainDb[0] - gainDb[1];
            okv (open >= 10.0, tag + ": measured, the probe passes " + fmt (open, 1) + " dB more "
                              + (dir == ef::Direction::up ? "under the louder" : "under the softer") + " burst (>= 10)");
            // ...and the release brings it home: at rest within 0.6 s of the second burst's end.
            const double rest = dir == ef::Direction::up ? 150.0 : 3000.0;
            const double end  = cut[(std::size_t) ((int) ((t0[1] + len + 0.6) * fs) / 16)];
            okv (std::fabs (end / rest - 1.0) < 0.01, tag + ": back at rest 0.6 s after the note (" + fmt (end, 1) + " Hz)");
        }
}

//==============================================================================
void stability()
{
    group ("stable and finite at Q 12 under the fastest modulation, and on a full-scale slap");
    // (a) THE KERNEL, driven harder than the stage ever can: the cutoff flips between its two extremes (20 Hz and
    // 0.45 fs) EVERY sample, and on a random schedule, at Q 12, on full-scale noise for ten seconds.
    const double fs = 48000.0;
    for (int pattern = 0; pattern < 2; ++pattern)
    {
        ef::detail::SweptSvf f;
        const double gLo = std::tan (kPi * 20.0 / fs), gHi = std::tan (kPi * 0.45), k = 1.0 / 12.0;
        felitronics::test::law11::Rng r { 0xC0FFEEu + (std::uint32_t) pattern };
        double peak = 0.0; bool finite = true;
        const int n = (int) (10.0 * fs), tail = (int) (2.0 * fs);
        double tailPeak = 0.0;
        for (int i = 0; i < n + tail; ++i)
        {
            const bool hi = pattern == 0 ? (i & 1) != 0 : r.next() > 0.0f;
            f.design (hi ? gHi : gLo, k);
            double lp = 0.0, bp = 0.0;
            f.tick (0, i < n ? (double) r.next() : 0.0, lp, bp);
            if ((i & 63) == 63) f.flush (1);
            finite = finite && std::isfinite (lp) && std::isfinite (bp);
            const double m = std::max (std::fabs (lp), std::fabs (bp));
            if (i < n) peak = std::max (peak, m); else if (i > n + tail / 2) tailPeak = std::max (tailPeak, m);
        }
        okv (finite && peak <= 32.0, std::string (pattern == 0 ? "every-sample flips" : "random flips")
                                    + ": finite, and the largest tap is " + fmt (peak, 2) + " on +-1 noise (bound 32)");
        ok (tailPeak == 0.0, "...and once the input stops it decays to exact zero within a second ("
                             + std::to_string (tailPeak) + ")");
    }

    // (b) THE STAGE at its fastest: attack 0.1 ms, release 5 ms, the widest range, Q 12, +36 dB, both modes and
    // directions, on a full-scale slap — a 41 Hz square at 0 dBFS cut by full-scale clicks of both signs.
    for (auto mode : { ef::Mode::lowPass, ef::Mode::bandPass })
        for (auto dir : { ef::Direction::up, ef::Direction::down })
        {
            Stage s;
            run (s.prepare (fs, 256, 2));
            Params p = voice (mode, dir, 12.0f);
            p.attackMs = 0.1f; p.releaseMs = 5.0f; p.sweepLoHz = 20.0f; p.sweepHiHz = 12000.0f;
            p.sensitivityDb = 36.0f; p.keepLows = 0.0f; p.mix = 1.0f;
            s.setParams (p);
            const int n = (int) (10.0 * fs);
            Buf io (2, std::vector<float> ((std::size_t) n, 0.0f));
            for (int i = 0; i < (int) (2.0 * fs); ++i)
            {
                const float sq = std::fmod (i * 41.0 / fs, 1.0) < 0.5 ? 1.0f : -1.0f;
                const bool click = i % 2400 < 3;
                io[0][(std::size_t) i] = click ? (i % 4800 < 3 ? 1.0f : -1.0f) : ((i / 12000) % 2 ? sq : 0.0f);
                io[1][(std::size_t) i] = -io[0][(std::size_t) i];
            }
            render (s, io, 256);
            double peak = 0.0, tailPeak = 0.0; bool finite = true;
            for (int c = 0; c < 2; ++c)
                for (int i = 0; i < n; ++i)
                {
                    const double v = io[(std::size_t) c][(std::size_t) i];
                    finite = finite && std::isfinite (v);
                    if (i < (int) (2.0 * fs)) peak = std::max (peak, std::fabs (v));
                    else if (i > (int) (9.0 * fs)) tailPeak = std::max (tailPeak, std::fabs (v));
                }
            const std::string tag = std::string (mode == ef::Mode::lowPass ? "low pass " : "band pass ")
                                  + (dir == ef::Direction::up ? "up" : "down");
            okv (finite && peak <= 2.0, tag + ": finite, peak " + fmt (dbOf (peak)) + " dBFS on a 0 dBFS slap (bound +6)");
            ok (tailPeak == 0.0, tag + ": the ring dies to exact zero (checked 7 s after the slap; a 20 Hz pole at Q 12 "
                                       "rings for ~6.6 s down to the 1e-15 flush)");
        }
}

//==============================================================================
// NO ZIPPER. A 110 Hz tone in Hann-gated bursts sweeps the cutoff open and shut; the swept output may only hold
// the tone and its slow modulation. What it must NOT hold is energy far above both — a cutoff that moved in steps
// would put its steps' broadband edges there. Measured through an 8-pole high pass at 6 kHz.
void noZipper()
{
    group ("no zipper: a swept output holds nothing far above the tone");
    for (double fs : kRates)
    {
        Stage s;
        run (s.prepare (fs, 64, 1));
        Params p = voice (ef::Mode::lowPass, ef::Direction::up, 4.0f);
        p.sensitivityDb = 6.0f; p.keepLows = 0.0f; p.mix = 1.0f;
        s.setParams (p);
        const int n = (int) (2.0 * fs);
        std::vector<float> x ((std::size_t) n, 0.0f);
        for (int b = 0; b < 4; ++b) addBurst (x, fs, 0.05 + 0.5 * b, 0.35, 0.05, 110.0, b % 2 ? 0.5 : 0.2);
        std::vector<float> y = x;
        float* io[1] = { y.data() };
        run (s.process (io, 1, n));

        double lo = 1e9, hi = -1e9;
        Stage t;   // the same stage again, read per block, so the test knows the sweep really went wide
        run (t.prepare (fs, 64, 1)); t.setParams (p);
        std::vector<float> z = x;
        for (int i = 0; i < n; i += 64)
        {
            float* q[1] = { z.data() + i };
            run (t.process (q, 1, std::min (64, n - i)));
            lo = std::min (lo, (double) t.shownCutoffHz()); hi = std::max (hi, (double) t.shownCutoffHz());
        }

        const auto hfRatioDb = [&] (const std::vector<float>& v)
        {
            eq::DeterministicSvf hp[4];
            for (auto& h : hp) { h.prepare (fs, 1); h.setParams (eq::FilterType::HighPass, 6000.0, 0.70710678, 0.0); }
            double eh = 0.0, et = 0.0;
            for (float sample : v)
            {
                float u = sample;
                for (auto& h : hp) u = h.processSample (0, u);
                eh += (double) u * u; et += (double) sample * sample;
            }
            return 10.0 * std::log10 (std::max (eh, 1e-300) / et);
        };
        const double out = hfRatioDb (y), in = hfRatioDb (x);
        okv (hi / lo > 10.0, fmt (fs, 0) + " Hz: precondition — the sweep went wide (" + fmt (lo, 0) + " ... " + fmt (hi, 0) + " Hz)");
        okv (out <= -120.0, fmt (fs, 0) + " Hz: energy above 6 kHz is " + fmt (out, 1) + " dB under the whole (bound -120; the input's own: "
                           + fmt (in, 1) + ")");
    }
}

//==============================================================================
// THE CLEAN LOWS. A band pass driven wide open is the case that loses the note; keepLows 1 must hold the low end.
void keepsLows()
{
    group ("the clean lows: the low end stays within a stated dB, and without them it does not");
    const double hzs[] = { 30.87, 41.2, 55.0, 73.42 };
    for (double fs : kRates)
        for (auto mode : { ef::Mode::bandPass, ef::Mode::lowPass })
        {
            double worstKeep = 0.0, leastLoss = 1e9;
            for (double hz : hzs)
                for (float keep : { 1.0f, 0.0f })
                {
                    Stage s;
                    run (s.prepare (fs, 256, 1));
                    Params p = ef::paramsForAmount (1.0f, voice (mode, ef::Direction::up, 12.0f));
                    p.keepLows = keep;
                    s.setParams (p);
                    const int n = (int) (1.5 * fs);
                    std::vector<float> x ((std::size_t) n, 0.0f);
                    addBurst (x, fs, 0.0, 1.5, 0.02, hz, 0.25);
                    std::vector<float> y = x;
                    float* io[1] = { y.data() };
                    run (s.process (io, 1, n));
                    const int a = (int) (0.8 * fs), b = (int) (1.4 * fs);
                    const double g = dbOf (toneAmp (y, a, b, hz, fs) / toneAmp (x, a, b, hz, fs));
                    if (keep > 0.5f) worstKeep = std::max (worstKeep, std::fabs (g));
                    else             leastLoss = std::min (leastLoss, -g);
                }
            const std::string tag = fmt (fs, 0) + " Hz " + (mode == ef::Mode::bandPass ? "band pass" : "low pass Q 12");
            okv (worstKeep <= 1.5, tag + ": keepLows 1 holds 31 ... 73 Hz within " + fmt (worstKeep) + " dB of the input (bound 1.5)");
            okv (leastLoss >= 9.0, tag + ": ...where keepLows 0 loses at least " + fmt (leastLoss, 1) + " dB (>= 9)");
        }
}

//==============================================================================
void reslicing()
{
    group ("law 8a: the output does not depend on how the stream is sliced");
    for (double fs : kRates)
    {
        const int n = (int) (1.5 * fs);
        Buf x (2, std::vector<float> ((std::size_t) n, 0.0f));
        for (int b = 0; b < 3; ++b) { addBurst (x[0], fs, 0.1 + 0.45 * b, 0.3, 0.005, 55.0 * (b + 1), 0.4); addBurst (x[1], fs, 0.1 + 0.45 * b, 0.3, 0.005, 82.0, 0.2); }
        const Params p = ef::paramsForAmount (0.8f, voice (ef::Mode::bandPass, ef::Direction::up, 9.0f));
        Stage a, b;
        run (a.prepare (fs, n, 2)); run (b.prepare (fs, 512, 2));
        a.setParams (p); b.setParams (p);
        Buf ya = x, yb = x;
        float* pa[2] = { ya[0].data(), ya[1].data() };
        run (a.process (pa, 2, n));
        felitronics::test::law11::Rng r { 0x51CEu };
        for (int at = 0; at < n;)
        {
            const int len = std::min (n - at, 1 + (int) ((r.next() * 0.5f + 0.5f) * 300.0f));
            float* pb[2] = { yb[0].data() + at, yb[1].data() + at };
            run (b.process (pb, 2, len));
            at += len;
        }
        okv (std::memcmp (ya[0].data(), yb[0].data(), (std::size_t) n * 4) == 0
            && std::memcmp (ya[1].data(), yb[1].data(), (std::size_t) n * 4) == 0,
            fmt (fs, 0) + " Hz: one whole-file call == 1 ... 300-sample calls, bit for bit");
    }
}

//==============================================================================
void stoppedChannelPauseReset()
{
    group ("law 11a / 11c / reset: a stopped channel returns silent, a pause is silence, a reset is a restart");
    const double fs = 48000.0;
    const Params p = ef::paramsForAmount (1.0f, voice (ef::Mode::bandPass, ef::Direction::up, 12.0f));
    {
        Stage s;
        run (s.prepare (fs, 256, 2)); s.setParams (p);
        Buf io (2, std::vector<float> (4800, 0.0f));
        addBurst (io[0], fs, 0.0, 0.1, 0.002, 200.0, 0.9); addBurst (io[1], fs, 0.0, 0.1, 0.002, 200.0, 0.9);
        render (s, io, 256);
        std::vector<float> mono (4800, 0.3f);
        float* m[1] = { mono.data() };
        run (s.process (m, 1, 4800));          // the right channel stops here, ringing
        Buf back (2, std::vector<float> (2048, 0.0f));
        back[0].assign (2048, 0.3f);
        float* b[2] = { back[0].data(), back[1].data() };
        run (s.process (b, 2, 2048));
        double leak = 0.0;
        for (float v : back[1]) leak = std::max (leak, (double) std::fabs (v));
        ok (leak == 0.0, "a channel that stopped while ringing comes back into silence with exact zero (law 11a)");
    }

    // 11c: two stages at rest, a glide put in flight (mix, resonance, direction, mode, range), then A pauses and
    // B hears silence for the same length; the shared state, and the next note, must be the same bits.
    using felitronics::test::law11c::kGaps;
    int cases = 0, same = 0, moved = 0;
    for (int gap : kGaps)
    {
        Stage A, B, C;
        for (Stage* s : { &A, &B, &C })
        {
            run (s->prepare (fs, 256, 2)); s->setParams (p);
            Buf z (2, std::vector<float> (256, 0.0f));
            render (*s, z, 256);
            Params q = p; q.mix = 0.4f; q.resonance = 2.0f; q.direction = ef::Direction::down; q.mode = ef::Mode::lowPass;
            q.sweepHiHz = 1500.0f;
            s->setParams (q);
        }
        run (A.process (nullptr, 0, gap));
        Buf z (2, std::vector<float> ((std::size_t) gap, 0.0f));
        render (B, z, 256);
        const bool sharedSame = core::sameBits (A.shownCutoffHz(), B.shownCutoffHz())
                             && A.currentMix() == B.currentMix() && A.currentResonance() == B.currentResonance();
        const bool sharedMoved = A.currentMix() != C.currentMix();
        Buf ta (2, std::vector<float> (4096, 0.0f));
        addBurst (ta[0], fs, 0.0, 4096 / fs, 0.001, 150.0, 0.7); addBurst (ta[1], fs, 0.0, 4096 / fs, 0.001, 97.0, 0.5);
        Buf tb = ta, tc = ta;
        render (A, ta, 256); render (B, tb, 256); render (C, tc, 256);
        const bool eqAB = std::memcmp (ta[0].data(), tb[0].data(), 4096 * 4) == 0 && std::memcmp (ta[1].data(), tb[1].data(), 4096 * 4) == 0;
        const bool neC  = std::memcmp (ta[0].data(), tc[0].data(), 4096 * 4) != 0;
        ++cases; if (eqAB && sharedSame) ++same; if (neC && sharedMoved) ++moved;
        if (! (eqAB && sharedSame)) ok (false, "pause != silence at gap " + std::to_string (gap));
    }
    ok (same == cases, "a pause is silence: the shared state and the next note are the same bits, all "
                       + std::to_string (cases) + " gaps (law 11c)");
    ok (moved == cases, "...and every gap MOVED the shared state (it is not a blind comparison)");

    // reset is a restart: two stages fed different programmes answer the same programme alike after reset().
    {
        Stage a, b;
        run (a.prepare (fs, 256, 2)); run (b.prepare (fs, 256, 2));
        a.setParams (p); b.setParams (p);
        Buf xa (2, std::vector<float> (9600, 0.0f)), xb = xa;
        addBurst (xa[0], fs, 0.0, 0.2, 0.002, 60.0, 0.9); addBurst (xb[1], fs, 0.0, 0.2, 0.002, 300.0, 0.2);
        render (a, xa, 256); render (b, xb, 256);
        a.reset(); b.reset(); b.reset();
        Buf ya (2, std::vector<float> (4800, 0.0f));
        addBurst (ya[0], fs, 0.0, 0.1, 0.002, 98.0, 0.6); ya[1] = ya[0];
        Buf yb = ya;
        render (a, ya, 256); render (b, yb, 256);
        ok (std::memcmp (ya[0].data(), yb[0].data(), 4800 * 4) == 0 && std::memcmp (ya[1].data(), yb[1].data(), 4800 * 4) == 0,
            "after reset() two stages with different histories answer the same bits");
    }
}

//==============================================================================
void poisonAndLaw8()
{
    group ("poison: a non-finite input is the dry's and nobody else's; law 8: no subnormal state, ever");
    const double fs = 48000.0;
    const Params p = ef::paramsForAmount (1.0f, voice (ef::Mode::bandPass, ef::Direction::up, 12.0f));
    {
        const int n = 24000;
        Buf x (2, std::vector<float> ((std::size_t) n, 0.0f));
        addBurst (x[0], fs, 0.0, 0.5, 0.002, 73.4, 0.5); addBurst (x[1], fs, 0.0, 0.5, 0.002, 110.0, 0.5);
        Buf poisoned = x;
        const int at[] = { 300, 5000, 5001, 12345 };
        const float bad[] = { std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(),
                              -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN() };
        for (int k = 0; k < 4; ++k) { poisoned[(std::size_t) (k & 1)][(std::size_t) at[k]] = bad[k]; x[(std::size_t) (k & 1)][(std::size_t) at[k]] = 0.0f; }
        Stage a, b;
        run (a.prepare (fs, 256, 2)); run (b.prepare (fs, 256, 2)); a.setParams (p); b.setParams (p);
        render (a, poisoned, 256); render (b, x, 256);
        int differ = 0, kept = 0;
        for (int c = 0; c < 2; ++c)
            for (int i = 0; i < n; ++i)
            {
                const float u = poisoned[(std::size_t) c][(std::size_t) i], v = x[(std::size_t) c][(std::size_t) i];
                if (! core::sameBits (u, v)) { ++differ; if (! std::isfinite (u)) ++kept; }
            }
        ok (differ == 4 && kept == 4, "NaN / +inf / -inf in: those 4 samples go out as they came, and every other sample is the "
                                      "bits of the same stream with zeros there (" + std::to_string (differ) + " differ)");
    }
    {
        Stage s;
        run (s.prepare (fs, 256, 1)); s.setParams (p);
        std::vector<float> big (512);
        for (int i = 0; i < 512; ++i) big[(std::size_t) i] = i % 2 ? 3.0e38f : -3.0e38f;
        float* q[1] = { big.data() };
        run (s.process (q, 1, 512));
        std::vector<float> after (4800, 0.0f);
        addBurst (after, fs, 0.0, 0.1, 0.002, 98.0, 0.5);
        float* r[1] = { after.data() };
        run (s.process (r, 1, 4800));
        bool finite = true;
        for (float v : after) finite = finite && std::isfinite (v);
        ok (finite, "a finite +-3e38 that overflows the filter is healed at the end of its call: the next call is finite");
    }
    // Law 8: a resonant ring, then a minute of silence — in blocks, and as ONE call (the 8a case: a whole-file
    // call must flush inside itself).
    for (int whole = 0; whole < 2; ++whole)
    {
        Stage s;
        run (s.prepare (fs, whole ? (int) (61.0 * fs) : 512, 2)); s.setParams (p);
        Buf io (2, std::vector<float> ((std::size_t) (61.0 * fs), 0.0f));
        addBurst (io[0], fs, 0.0, 0.5, 0.001, 140.0, 0.9); addBurst (io[1], fs, 0.0, 0.5, 0.001, 40.0, 0.9);
        bool sub = false, quiet = true;
        if (whole)
        {
            float* q[2] = { io[0].data(), io[1].data() };
            run (s.process (q, 2, (int) io[0].size()));
            sub = s.anySubnormalState();
        }
        else
            for (std::size_t at = 0; at < io[0].size(); at += 512)
            {
                float* q[2] = { io[0].data() + at, io[1].data() + at };
                run (s.process (q, 2, (int) std::min<std::size_t> (512, io[0].size() - at)));
                sub = sub || s.anySubnormalState();
            }
        for (int c = 0; c < 2; ++c)
            for (std::size_t i = (std::size_t) (2.0 * fs); i < io[0].size(); ++i)
                quiet = quiet && io[(std::size_t) c][i] == 0.0f;
        ok (! sub, std::string (whole ? "one 61 s call" : "512-sample calls") + ": no state reaches the subnormal band");
        ok (quiet, std::string (whole ? "one 61 s call" : "512-sample calls") + ": from 1.5 s after the ring the output is exact zero");
    }
}

//==============================================================================
void monoStereoAndLanding()
{
    group ("mono and stereo; the glides land");
    const double fs = 48000.0;
    const Params p = ef::paramsForAmount (0.7f, voice (ef::Mode::lowPass, ef::Direction::up, 5.0f));
    {
        const int n = 19200;
        std::vector<float> m ((std::size_t) n, 0.0f);
        addBurst (m, fs, 0.0, 0.4, 0.003, 61.7, 0.6);
        Buf st (2, m);
        Stage a, b;
        run (a.prepare (fs, 256, 1)); run (b.prepare (fs, 256, 2)); a.setParams (p); b.setParams (p);
        Buf mo (1, m);
        render (a, mo, 256); render (b, st, 256);
        ok (std::memcmp (mo[0].data(), st[0].data(), (std::size_t) n * 4) == 0 && std::memcmp (mo[0].data(), st[1].data(), (std::size_t) n * 4) == 0,
            "stereo with L == R is the mono stage on each side, bit for bit");
    }
    {
        // Linked: the louder channel moves both sides' cutoff; the quieter side alone would not open it so far.
        Stage a, b;
        run (a.prepare (fs, 256, 2)); run (b.prepare (fs, 256, 1)); a.setParams (p); b.setParams (p);
        Buf st (2, std::vector<float> (9600, 0.0f));
        addBurst (st[0], fs, 0.0, 0.2, 0.003, 98.0, 0.1); addBurst (st[1], fs, 0.0, 0.2, 0.003, 98.0, 0.5);
        Buf loud (1, st[1]);
        render (a, st, 256); render (b, loud, 256);
        ok (core::sameBits (a.shownCutoffHz(), b.shownCutoffHz()), "linked: one envelope, the louder side's, moves both channels");
    }
    {
        Stage s;
        run (s.prepare (fs, 256, 1)); s.setParams (p);
        std::vector<float> x (48000, 0.0f);
        float* q[1] = { x.data() };
        run (s.process (q, 1, 256));
        Params r = p; r.mix = 0.3f; r.resonance = 9.0f;
        s.setParams (r);
        run (s.process (q, 1, 48000));
        ok (s.currentMix() == (double) 0.3f && s.currentResonance() == 9.0,
            "a second after a move the glides stand exactly on their targets");
        ok (! s.anySubnormalState(), "...with nothing parked in the subnormal band");
    }
}

//==============================================================================
void noAllocation()
{
    group ("RT: process, reset and the setters never allocate");
    Stage s;
    run (s.prepare (96000.0, 256, 2));
    std::vector<float> l (1024, 0.25f), r (1024, -0.25f);
    float* lr[2] = { l.data(), r.data() };
    const long long before = alloc::count.load();
    s.setParams (ef::paramsForAmount (0.9f, voice (ef::Mode::bandPass)));
    run (s.process (lr, 2, 1024));
    run (s.process (lr, 1, 1024));
    run (s.process (nullptr, 0, 1024));
    s.setParams (ef::paramsForAmount (0.0f));
    run (s.process (lr, 2, 1024));
    for (int k = 0; k < 40; ++k) run (s.process (lr, 2, 1024));   // glides out and sleeps
    s.setParams (ef::paramsForAmount (0.5f));
    run (s.process (lr, 2, 1024));                                 // wakes
    s.reset(); s.reset();
    run (s.process (lr, 2, 1024));
    felitronics::test::okNoAlloc (alloc::count.load() == before,
                                  "no allocation: setParams, process (stereo, mono, pause, sleep, wake), reset");
}
} // namespace

int main()
{
    std::printf ("felitronics::envfilter::EnvelopeFilter\n");
    sliderLaw();
    contract();
    dryAtZero();
    sliderGlides();
    tracking();
    stability();
    noZipper();
    keepsLows();
    reslicing();
    stoppedChannelPauseReset();
    poisonAndLaw8();
    monoStereoAndLanding();
    noAllocation();
    return felitronics::test::report();
}
