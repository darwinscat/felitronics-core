// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// ClipDetector self-tests. The truth is BY CONSTRUCTION: a signal is clamped here, at a known ceiling, and the clamp
// records which samples it changed — so "was it found" is a comparison with a fact, not with another detector.
//
//   · clamped material (float / 24 / 16 bit, dithered, turned down by up to 40 dB, 44.1–192 kHz) — found run by run
//   · clean material that is loud, flat or both (full-scale sines at every phase, the 16-bit rail, squares, PWM,
//     quantised sub-bass, sample-and-hold, 8 bit, noise, silence, gates, DC) — never reported
//   · a whole-file REFERENCE written straight from the rule, independent of the streaming machinery (deques, pending
//     queue, holes, finish) — the engine must null it run for run on randomised material
//   · slicing, polarity, channel order and power-of-two level invariance; non-finite input; the law-11 contract;
//     no allocation in process()/finish(); the published storage is the allocated storage.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included
#include <felitronics/analysis/ClipDetector.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace felitronics;
using CD = analysis::ClipDetector;
using Planes = std::vector<std::vector<float>>;
static constexpr double kTau2Pi = 6.283185307179586;

// ============================================================================================== signals + truth

struct Truth { int ch; long start, len; };

// music-like: bass with harmonics and a slow envelope, a mid tone, filtered noise bursts (hats) — seeded
static std::vector<double> programme (long N, double sr, unsigned seed)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<double> U (-1.0, 1.0);
    std::vector<double> x ((std::size_t) N);
    const double fb = 38.0 + 30.0 * (0.5 + 0.5 * U (rng)), fm = 300.0 + 400.0 * (0.5 + 0.5 * U (rng));
    double lp = 0.0;
    const double a = std::exp (-kTau2Pi * 6000.0 / sr);
    for (long n = 0; n < N; ++n)
    {
        const double t = (double) n / sr;
        const double beat = std::fmod (t, 0.5);
        const double env = std::exp (-beat * 4.0);
        double v = env * (std::sin (kTau2Pi * fb * t) + 0.5 * std::sin (kTau2Pi * 2 * fb * t + 0.4) + 0.25 * std::sin (kTau2Pi * 3 * fb * t + 1.1));
        v += 0.3 * std::sin (kTau2Pi * fm * t) * (0.6 + 0.4 * std::sin (kTau2Pi * 0.3 * t));
        lp = a * lp + (1.0 - a) * U (rng);
        v += (std::fmod (t, 0.25) < 0.05 ? 0.25 : 0.03) * (U (rng) - lp);
        x[(std::size_t) n] = v;
    }
    return x;
}

static void normalisePeak (std::vector<double>& x, double peak)
{
    double m = 0.0; for (double v : x) m = std::max (m, std::fabs (v));
    for (double& v : x) v *= peak / m;
}

// clamp at [lo, hi]; the truth is exactly the samples the clamp changed
static std::vector<Truth> clampTruth (std::vector<double>& x, int ch, double lo, double hi)
{
    std::vector<Truth> t;
    long start = -1;
    for (long n = 0; n <= (long) x.size(); ++n)
    {
        bool c = false;
        if (n < (long) x.size())
        {
            double& v = x[(std::size_t) n];
            if (v > hi) { v = hi; c = true; } else if (v < lo) { v = lo; c = true; }
        }
        if (c && start < 0) start = n;
        if (! c && start >= 0) { t.push_back ({ ch, start, n - start }); start = -1; }
    }
    return t;
}

// deliver: 0 = float32, 16 / 24 = integer PCM read back as code / 2^(bits-1); tpdf = ±1 LSB triangular dither
static std::vector<float> deliver (const std::vector<double>& x, int bits, bool tpdf = false, unsigned seed = 1)
{
    std::vector<float> y (x.size());
    std::mt19937 rng (seed);
    std::uniform_real_distribution<double> U (-0.5, 0.5);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (bits == 0) { y[i] = (float) x[i]; continue; }
        const double S = std::ldexp (1.0, bits - 1);
        double k = std::nearbyint (x[i] * S + (tpdf ? U (rng) + U (rng) : 0.0));
        k = std::clamp (k, -S, S - 1.0);
        y[i] = (float) (k / S);
    }
    return y;
}

// ============================================================================================== running the engine

struct Report
{
    bool clipped = false, complete = true;
    std::vector<std::tuple<int, long, long, int, double, int>> runs;   // channel, start, length, sign, level, evidence (engine order)
    std::vector<double> peak, dc;
    std::vector<long long> perChRuns, clippedSamples, longest, nonFinite;
};

static Report runEngine (const Planes& x, double sr, int chunk /* 0 = whole, -1 = random */, unsigned seed = 7, int maxRuns = 1 << 20)
{
    CD d; d.setParams ({ maxRuns });
    const int nch = (int) x.size();
    test::run (d.prepare (sr, 512, nch));
    const long N = (long) x[0].size();
    std::mt19937 rng (seed);
    std::vector<const float*> p ((std::size_t) nch);
    for (long pos = 0; pos < N;)
    {
        long len = chunk == 0 ? N : chunk > 0 ? chunk : 1 + (long) (rng() % 5000u);
        len = std::min (len, N - pos);
        for (int c = 0; c < nch; ++c) p[(std::size_t) c] = x[(std::size_t) c].data() + pos;
        test::run (d.process (p.data(), nch, (int) len));
        pos += len;
    }
    d.finish();
    Report r;
    r.clipped = d.clipped(); r.complete = d.runsComplete();
    for (std::int64_t i = 0; i < d.storedRunCount(); ++i)
    {
        const auto& u = d.run (i);
        r.runs.emplace_back (u.channel, (long) u.start, (long) u.length, u.sign, u.level, (int) u.evidence);
    }
    for (int c = 0; c < nch; ++c)
    {
        r.peak.push_back (d.samplePeak (c)); r.dc.push_back (d.dcOffset (c));
        r.perChRuns.push_back (d.runCount (c)); r.clippedSamples.push_back (d.clippedSamples (c));
        r.longest.push_back (d.longestRun (c)); r.nonFinite.push_back (d.nonFiniteSamples (c));
    }
    return r;
}

static bool sameReport (const Report& a, const Report& b)
{
    if (a.clipped != b.clipped || a.complete != b.complete || a.runs.size() != b.runs.size()) return false;
    for (std::size_t i = 0; i < a.runs.size(); ++i) if (a.runs[i] != b.runs[i]) return false;       // order included
    for (std::size_t c = 0; c < a.peak.size(); ++c)
        if (a.peak[c] != b.peak[c] || a.dc[c] != b.dc[c] || a.perChRuns[c] != b.perChRuns[c] || a.clippedSamples[c] != b.clippedSamples[c]
            || a.longest[c] != b.longest[c] || a.nonFinite[c] != b.nonFinite[c]) return false;
    return true;
}

// ============================================================================================== the whole-file reference
// The rule written straight from its definition, over whole arrays: q(t) as a table, bands by direct scan, the window by
// direct scan, sides by index. Nothing is shared with the engine but the constants.
struct RefRun { int ch; long start, len; };

static std::vector<RefRun> reference (const Planes& planes, double sr)
{
    std::vector<RefRun> out;
    const long W = (long) (sr * CD::kWindowMs / 1000.0);
    for (std::size_t c = 0; c < planes.size(); ++c)
    {
        const std::vector<float>& y = planes[c];
        const long N = (long) y.size();
        auto fin = [&] (long i) { return i >= 0 && i < N && std::isfinite (y[(std::size_t) i]); };
        auto Y = [&] (long i) { return (double) y[(std::size_t) i]; };
        std::vector<double> qAt ((std::size_t) N);
        double step = std::numeric_limits<double>::infinity();
        int grid = 0;                                                     // by repeated doubling, not by bits
        for (long i = 0; i < N; ++i)
        {
            if (fin (i) && fin (i - 1)) { const double d = std::fabs (Y (i) - Y (i - 1)); if (d > 0) step = std::min (step, d); }
            if (fin (i) && Y (i) != 0.0)
            {
                double v = std::fabs (Y (i)); int k = 0;
                while (k <= 24 && v != std::floor (v)) { v *= 2.0; ++k; }
                grid = std::max (grid, k);
            }
            qAt[(std::size_t) i] = std::min (step, grid <= 24 ? std::ldexp (1.0, -grid) : std::numeric_limits<double>::infinity());
        }
        double ceil[2] = { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
        long k = 0;
        while (k < N)
        {
            if (! fin (k)) { ++k; continue; }
            double lo = Y (k), hi = lo, sum = lo;
            long e = k + 1;
            while (e < N && fin (e))
            {
                const double nlo = std::min (lo, Y (e)), nhi = std::max (hi, Y (e));
                if (nhi - nlo > CD::kBandQuanta * qAt[(std::size_t) e]) break;
                lo = nlo; hi = nhi; sum += Y (e); ++e;
            }
            const long L = e - k;
            const double qq = qAt[(std::size_t) std::min (N - 1, e + 2)], tau = CD::kBandQuanta * qq;
            if (L >= 2 && std::isfinite (qq))
            {
                const double m = sum / (double) L;
                const bool kin = fin (k - 1), kout = fin (e);
                int type = 0;
                if (kin || kout)
                {
                    if ((! kin || Y (k - 1) < lo) && (! kout || Y (e) < lo)) type = 1;
                    else if ((! kin || Y (k - 1) > hi) && (! kout || Y (e) > hi)) type = -1;
                }
                bool windowOk = type != 0;
                for (long i = std::max (0L, k - W); windowOk && i < std::min (N, e + W); ++i)
                    if ((i < k || i >= e) && fin (i) && (type > 0 ? Y (i) - hi > tau : lo - Y (i) > tau)) windowOk = false;
                if (windowOk)
                {
                    auto atLevel = [&] (long i) { return std::fabs (Y (i) - m) <= tau; };
                    auto isStep = [&] (long a, long b, long c2) {
                        if (! fin (a) || ! fin (b) || ! fin (c2)) return true;
                        return (std::fabs (Y (a) - Y (b)) <= tau && ! atLevel (b)) || (std::fabs (Y (b) - Y (c2)) <= tau && ! atLevel (b));
                    };
                    double act = 0; long na = 0;
                    for (long i = k - 1; i > k - 1 - CD::kFlankSteps && fin (i) && fin (i - 1); --i)
                        if (! atLevel (i) && ! atLevel (i - 1)) { act += std::fabs (Y (i) - Y (i - 1)); ++na; }
                    for (long i = e; i < e + CD::kFlankSteps && fin (i) && fin (i + 1); ++i)
                        if (! atLevel (i) && ! atLevel (i + 1)) { act += std::fabs (Y (i + 1) - Y (i)); ++na; }
                    const double sigma = na ? act / (double) na : 0.0;
                    double chance = 1.0;                                          // ((r + q)/sigma)^(L-1), stopped at the bound
                    const double ratio = sigma > 0 ? (hi - lo + qq) / sigma : 1.0;
                    for (long i = 0; ratio < 1.0 && i < L - 1 && chance > CD::kChance; ++i) chance *= ratio;
                    const bool unlikely = ratio < 1.0 && chance <= CD::kChance;
                    const double stretch = (double) (L + 1) / (double) (L - 1);
                    const double V = (tau + qq) * (stretch * stretch) + qq;
                    const double s = type;
                    bool ramp = L >= CD::kMinRampLength && unlikely && kin && kout;
                    if (ramp) ramp = ! isStep (k - 1, k - 2, k - 3) && s * (m - Y (k - 1)) >= V && ! isStep (e, e + 1, e + 2) && s * (m - Y (e)) >= V;
                    auto turns = [&] (long a, long b, long c2) {                // the growth of the outward steps
                        const double d1 = s * (m - Y (a)), d2 = s * (m - Y (b)), d3 = s * (m - Y (c2));
                        return d2 - d1 > 0 ? ((d3 - d2) - (d2 - d1)) / (d2 - d1) : 1e300;
                    };
                    if (ramp && hi != lo) ramp = turns (k - 1, k - 2, k - 3) <= CD::kRoughTurnRatio && turns (e, e + 1, e + 2) <= CD::kRoughTurnRatio;
                    const int ci = type > 0 ? 0 : 1;
                    const bool atCeil = std::isfinite (ceil[ci]) && std::fabs (m - ceil[ci]) <= tau;
                    if (ramp || atCeil) out.push_back ({ (int) c, k, L });
                    if (ramp) ceil[ci] = m;
                }
            }
            k = e;
        }
    }
    return out;
}

static std::vector<RefRun> sortedRuns (const Report& r)
{
    std::vector<RefRun> v;
    for (const auto& u : r.runs) v.push_back ({ std::get<0> (u), std::get<1> (u), std::get<2> (u) });
    std::sort (v.begin(), v.end(), [] (const RefRun& a, const RefRun& b) { return std::tie (a.ch, a.start) < std::tie (b.ch, b.start); });
    return v;
}

static bool sameRuns (const std::vector<RefRun>& a, const std::vector<RefRun>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (a[i].ch != b[i].ch || a[i].start != b[i].start || a[i].len != b[i].len) return false;
    return true;
}

// fraction of truth runs of length >= minLen overlapped by a reported run on their channel; and reported runs that
// overlap no truth run at all
static void recall (const std::vector<Truth>& truth, const Report& r, long minLen, double& frac, long& outside)
{
    long n = 0, hit = 0;
    std::vector<bool> used (r.runs.size(), false);
    for (const Truth& t : truth)
    {
        bool h = false;
        for (std::size_t i = 0; i < r.runs.size(); ++i)
        {
            const auto& u = r.runs[i];
            if (std::get<0> (u) == t.ch && std::get<1> (u) < t.start + t.len && t.start < std::get<1> (u) + std::get<2> (u)) { h = true; used[i] = true; }
        }
        if (t.len >= minLen) { ++n; if (h) ++hit; }
    }
    frac = n ? (double) hit / (double) n : 1.0;
    outside = (long) std::count (used.begin(), used.end(), false);
}

// ============================================================================================== tests

int main()
{
    std::printf ("felitronics::analysis ClipDetector tests\n");

    // ---------------------------------------------------------------------------------------------- found
    test::group ("a clamp is found run by run — float / 24 / 16 bit / dithered, turned down, at every rate");
    {
        struct Case { double sr; int bits; bool tpdf; double driveDb, attenDb; long minLen; double wantRecall; };
        // wantRecall is what this construction measures, pinned just under it: a regression moves it, a mutant falls through it.
        // Two rule changes moved them, and the pins moved with the numbers. The rough-band turn test: float at 0 dB 0.9951 ->
        // 0.9902 (a neighbour within 2q joins the band, which is then no longer exactly flat), dithered at -6 dB 1.0000 ->
        // 0.9750. The chance bound 1e-5 -> 1e-6 (real clean recordings, see the header): float at 0 dB 0.9902 -> 0.9853,
        // 96 kHz at -12 dB 0.9957 -> 0.9830, and the shallow 1 dB clamp 1.0000 -> 0.7500 — a shallow clamp is short, and
        // a short run is where a coincidence bound bites.
        const Case cases[] = {
            { 48000, 0,  false, 6,  0,   2, 0.98 },  { 48000, 24, false, 6, -1,   2, 0.995 }, { 44100, 16, false, 6,  0,  2, 0.99 },
            { 96000, 16, false, 3, -12,  3, 0.98 },  { 192000, 24, false, 12, -6, 2, 0.995 }, { 48000, 16, true, 6, -6,  3, 0.97 },
            { 44100, 0,  false, 1,  -3,  3, 0.70 },  { 48000, 16, false, 6, -40, 6, 0.95 },
        };
        unsigned seed = 11;
        for (const Case& cs : cases)
        {
            const long N = (long) (cs.sr * 3.0);
            Planes planes; std::vector<Truth> truth;
            for (int ch = 0; ch < 2; ++ch)
            {
                auto x = programme (N, cs.sr, seed++);
                normalisePeak (x, 1.0);
                for (double& v : x) v *= std::pow (10.0, cs.driveDb / 20.0);
                auto t = clampTruth (x, ch, -1.0, 1.0);
                truth.insert (truth.end(), t.begin(), t.end());
                for (double& v : x) v *= std::pow (10.0, cs.attenDb / 20.0);
                planes.push_back (deliver (x, cs.bits, cs.tpdf, seed));
            }
            const Report r = runEngine (planes, cs.sr, 4096);
            double frac; long outside;
            recall (truth, r, cs.minLen, frac, outside);
            char msg[256];
            std::snprintf (msg, sizeof msg, "sr %.0f bits %d%s drive %.0f dB atten %.0f dB: run recall (len >= %ld) %.4f >= %.3f, %ld outside truth, %zu truth runs",
                           cs.sr, cs.bits, cs.tpdf ? " tpdf" : "", cs.driveDb, cs.attenDb, cs.minLen, frac, cs.wantRecall, outside, truth.size());
            test::ok (r.clipped && frac >= cs.wantRecall, msg);
            // a crest that reaches the ceiling's code without passing it cannot be told from one that did: at -40 dB in 16 bit
            // the ceiling is 328 codes and tau is 2 of them. Those are the runs "outside", and they stay under 1 % of the truth.
            test::ok (outside * 100 <= (long) truth.size(), std::string ("runs outside the clamp are under 1 % of the truth runs — ") + msg);
        }
    }

    test::group ("a float clamp is found where it is: reported runs ARE clamp runs, start and length");
    {
        const double sr = 48000;
        auto x = programme (144000, sr, 5);
        normalisePeak (x, 1.0);
        for (double& v : x) v *= 2.0;
        const auto truth = clampTruth (x, 0, -1.0, 1.0);
        for (double& v : x) v *= 0.3;                                    // turned down by a non-power-of-two in float
        const Report r = runEngine ({ deliver (x, 0) }, sr, 1000);
        long exact = 0;
        for (const auto& u : r.runs)
            for (const Truth& t : truth) if (std::get<1> (u) == t.start && std::get<2> (u) == t.len) { ++exact; break; }
        double frac; long outside;
        recall (truth, r, 2, frac, outside);
        // The runs it misses here are the first few of each polarity: short runs between dips, before any ramp has proven
        // the ceiling (the ceiling only looks forward — measured 5 of 376, all starting before sample 900).
        char msg[200]; std::snprintf (msg, sizeof msg, "%ld of %zu reported runs match a clamp run exactly; recall (len >= 2) %.4f >= 0.985; %ld outside",
                                      exact, r.runs.size(), frac, outside);
        test::ok (r.runs.size() > 100 && exact == (long) r.runs.size() && frac >= 0.985 && outside == 0, msg);
    }

    // ---------------------------------------------------------------------------------------------- silent
    test::group ("clean material is never reported: loud, flat, stepped, quantised, noisy or silent");
    {
        int reported = 0, files = 0;
        auto expectClean = [&] (const Planes& p, double sr, const std::string& what) {
            ++files;
            const Report r = runEngine (p, sr, 0);
            if (r.clipped) { ++reported; test::ok (false, what + ": reported " + std::to_string (r.runs.size()) + " runs"); }
        };
        const double srs[] = { 44100, 48000, 96000, 192000 };
        for (double sr : srs)
            for (double f : { 20.0, 31.0, 100.0, 997.0, 1000.0, 3000.0, sr / 8, sr / 6, sr / 4, 15000.0 })
                for (int ph = 0; ph < 3; ++ph)
                    for (int bits : { 0, 16, 24 })
                    {
                        // phase 0: random; 1: crest exactly between two samples (an exactly equal pair); 2: fs/4 at 45°
                        const long N = std::max ((long) (sr / 4), (long) (6 * sr / f));
                        const double phase = ph == 0 ? 1.234 : ph == 1 ? kTau2Pi / 4 - kTau2Pi * f * 0.5 / sr : kTau2Pi / 8;
                        std::vector<double> x ((std::size_t) N);
                        for (long n = 0; n < N; ++n) x[(std::size_t) n] = std::sin (kTau2Pi * f * n / sr + phase);
                        normalisePeak (x, bits == 16 ? 32767.0 / 32768.0 : bits == 24 ? 8388607.0 / 8388608.0 : 1.0);
                        expectClean ({ deliver (x, bits) }, sr, "sine " + std::to_string (f) + " Hz at the rail, phase " + std::to_string (ph) + ", bits " + std::to_string (bits));
                    }
        for (double sr : { 44100.0, 48000.0, 96000.0 })
            for (double f : { 50.0, 110.0, 1000.0, 5512.5 })
                for (double duty : { 0.5, 0.1, 0.33, 0.9 })
                    for (double amp : { 1.0, 0.5, 0.03 })
                        for (int bits : { 0, 16 })
                        {
                            const long N = (long) sr / 2;
                            std::vector<double> x ((std::size_t) N);
                            for (long n = 0; n < N; ++n) x[(std::size_t) n] = std::fmod (f * n / sr, 1.0) < duty ? amp : -amp;
                            if (bits == 16) for (double& v : x) v = std::min (v, 32767.0 / 32768.0);
                            expectClean ({ deliver (x, bits) }, sr, "pulse " + std::to_string (f) + " Hz duty " + std::to_string (duty));
                        }
        {   // PWM with a swept duty, fractional period
            const double sr = 48000; const long N = 96000;
            std::vector<double> x ((std::size_t) N); double ph = 0;
            for (long n = 0; n < N; ++n) { ph = std::fmod (ph + 220.7 / sr, 1.0); x[(std::size_t) n] = ph < 0.5 + 0.35 * std::sin (kTau2Pi * 0.5 * n / sr) ? 0.5 : -0.5; }
            expectClean ({ deliver (x, 16) }, sr, "PWM, swept duty");
        }
        for (double sr : { 96000.0, 192000.0, 48000.0 })                   // sub-bass: long identical 16-bit codes at every crest
            for (double f : { 20.0, 30.0, 41.2, 55.0 })
                for (double level : { 1.0, 0.5, 0.1, 0.01 })
                    for (int bits : { 16, 24 })
                    {
                        const long N = (long) sr;
                        std::vector<double> x ((std::size_t) N);
                        for (long n = 0; n < N; ++n) { const double t = n / sr; x[(std::size_t) n] = std::sin (kTau2Pi * f * t) * (0.8 + 0.2 * std::cos (kTau2Pi * 0.7 * t)) + 0.15 * std::sin (kTau2Pi * 2 * f * t + 0.3); }
                        normalisePeak (x, level * 32767.0 / 32768.0);
                        expectClean ({ deliver (x, bits) }, sr, "sub-bass " + std::to_string (f) + " Hz level " + std::to_string (level));
                    }
        struct Shape { std::vector<int> pre, band, post; };
        const Shape shapes[] = {
            // the shape measured on a false run of real 16-bit bass: outward distances 7, 15, 31, 55, 89, 134 quanta (a turn of
            // 1.0), seven samples within two codes
            { { 134, 89, 55, 31, 15, 7 }, { -1, 1, 1, 1, 1, 1, -1 }, { 8, 18, 37, 67, 111, 170 } },
            // a gentler turn, 0.75 and 0.7 — the middle of what the false runs showed — on flanks steep enough to pass the chance bound
            { { 420, 300, 200, 120, 60, 29, 15, 7 }, { -1, 1, 1, 1, 1, 1, -1 }, { 8, 18, 35, 70, 130, 210, 310, 430 } },
        };
        for (const Shape& sh : shapes)
            for (double peakCodes : { 12000.0, 3000.0 })
            {   // a sustained note whose crest is flatter than a parabola: it TURNS at its top, entering faster than the
                // parabolic bound allows. Bound (1) passes it; the turn test does not.
                std::vector<double> half;
                for (int d : sh.pre) half.push_back (-d);
                for (int b : sh.band) half.push_back (b);
                for (int d : sh.post) half.push_back (-d);
                const double from = -sh.post.back(), to = -2.0 * peakCodes + sh.post.back();
                for (int i = 1; i <= 220; ++i) half.push_back (from + (to - from) * (0.5 - 0.5 * std::cos (kTau2Pi * 0.5 * i / 220.0)));
                std::vector<double> period = half;
                for (double v : half) period.push_back (-2.0 * peakCodes - v);
                std::vector<float> y (96000);
                for (std::size_t n = 0; n < y.size(); ++n) y[n] = (float) (std::nearbyint (peakCodes + period[n % period.size()]) / 32768.0);
                expectClean ({ y }, 48000, "a note whose crest turns flatter than a parabola (shape " + std::to_string (&sh - shapes) + "), peak "
                             + std::to_string (peakCodes) + " codes");
            }
        {   // a crest flat to the fourth order: 2/3 cos − 1/6 cos 2θ repeats one 16-bit code for ~119 samples at 96 k
            const double sr = 96000; const long N = 96000;
            std::vector<double> x ((std::size_t) N);
            for (long n = 0; n < N; ++n) { const double th = kTau2Pi * 30.0 * n / sr; x[(std::size_t) n] = 2.0 / 3.0 * std::cos (th) - std::cos (2 * th) / 6.0; }
            expectClean ({ deliver (x, 16) }, sr, "fourth-order-flat crest");
        }
        {   // programme: peak-normalised to the rail (float and 16 bit), 8 bit, sample-and-hold, gated, on DC, as noise
            const double sr = 48000; const long N = 240000;
            for (unsigned seed : { 1u, 2u, 3u })
            {
                auto x = programme (N, sr, seed);
                normalisePeak (x, 32767.0 / 32768.0);
                expectClean ({ deliver (x, 16) }, sr, "programme at the 16-bit rail");
                expectClean ({ deliver (x, 0) }, sr, "programme at full scale, float");
                expectClean ({ deliver (x, 24) }, sr, "programme at the 24-bit rail");
                std::vector<double> h = x;
                for (std::size_t i = 0; i < h.size(); ++i) h[i] = x[i - i % 4];
                expectClean ({ deliver (h, 16) }, sr, "sample-and-hold x4");
                std::vector<double> b = x;
                for (double& v : b) v = std::nearbyint (v * 127.0) / 128.0;
                expectClean ({ deliver (b, 0) }, sr, "8-bit programme");
                std::vector<double> g = x;
                for (std::size_t i = 0; i < g.size(); ++i) if ((i / 12000) % 2) g[i] = 0.0;
                expectClean ({ deliver (g, 16) }, sr, "hard-gated programme");
                std::vector<double> d = x;
                for (double& v : d) v = 0.1 + 0.5 * v;
                expectClean ({ deliver (d, 16) }, sr, "programme on +0.1 DC");
                expectClean ({ deliver (x, 16, true, seed) }, sr, "programme dithered at -0.0 dB (never clamped: the dither stays inside)");
            }
            // noise-like drums (snare and hats are white noise) at 16 bit, loud and very quiet: the material where samples
            // land inside a band by coincidence. A chance bound of 1e-2 reports runs here; 1e-5 reports none.
            for (unsigned dseed : { 2u, 4u })
                for (double dpeak : { 0.5, 0.01 })
                {
                    const long ND = (long) (sr * 20);
                    std::mt19937 dr (dseed); std::uniform_real_distribution<double> DU (-1.0, 1.0);
                    std::vector<double> dx ((std::size_t) ND);
                    for (long n = 0; n < ND; ++n)
                    {
                        const double t = n / sr, bt = std::fmod (t, 0.5), stt = std::fmod (t + 0.25, 0.5), ht = std::fmod (t, 0.25);
                        double v = 0.9 * std::exp (-bt * 8) * std::sin (kTau2Pi * (45 * bt + 90 * (1 - std::exp (-bt * 28)) / 28));
                        v += 0.45 * std::exp (-stt * 22) * (DU (dr) + 0.6 * std::sin (kTau2Pi * 185 * stt));
                        v += 0.15 * std::exp (-ht * 90) * DU (dr);
                        v += 0.05 * (std::sin (kTau2Pi * 220 * t) + std::sin (kTau2Pi * 277.18 * t) + std::sin (kTau2Pi * 329.63 * t));
                        dx[(std::size_t) n] = v;
                    }
                    normalisePeak (dx, dpeak);
                    expectClean ({ deliver (dx, 16) }, sr, "noise-like drums at peak " + std::to_string (dpeak) + ", seed " + std::to_string (dseed));
                }
            std::mt19937 rng (9); std::uniform_real_distribution<double> U (-0.999, 0.999);
            std::vector<double> w ((std::size_t) N); for (double& v : w) v = U (rng);
            expectClean ({ deliver (w, 16) }, sr, "white noise");
            expectClean ({ std::vector<float> ((std::size_t) N, 0.0f) }, sr, "digital silence");
            std::vector<double> imp ((std::size_t) N, 0.0); imp[1000] = 1.0; imp[50000] = -1.0;
            expectClean ({ deliver (imp, 0) }, sr, "single impulses");
            for (int bits : { 0, 16 })
            {   // slow music first (so the smallest step is small), then digital silence carrying one-sample clicks
                std::vector<double> mc ((std::size_t) N, 0.0);
                for (long n = 0; n < N / 2; ++n) mc[(std::size_t) n] = 0.5 * std::sin (kTau2Pi * 55.0 * n / sr) * std::exp (-3.0 * n / sr);
                for (long n = N / 2; n < N; n += 700) mc[(std::size_t) n] = (n / 700) % 2 ? 0.6 : -0.6;
                expectClean ({ deliver (mc, bits) }, sr, "one-sample clicks in digital silence after slow music, bits " + std::to_string (bits));
                std::vector<double> cut ((std::size_t) N, 0.0);                 // loud music cut hard to silence, then clicks
                auto pg = programme (N / 2, sr, 17u + (unsigned) bits); normalisePeak (pg, 0.9);
                for (long n = 0; n < N / 2; ++n) cut[(std::size_t) n] = pg[(std::size_t) n];
                for (long n = N / 2 + 350; n < N; n += 700) cut[(std::size_t) n] = (n / 700) % 2 ? 0.5 : -0.5;
                expectClean ({ deliver (cut, bits) }, sr, "loud music cut hard to digital silence, then one-sample clicks, bits " + std::to_string (bits));
            }
        }
        char msg[96]; std::snprintf (msg, sizeof msg, "%d of %d clean signals reported", reported, files);
        test::ok (reported == 0 && files > 700, msg);
    }

    // ---------------------------------------------------------------------------------------------- the reference null
    test::group ("the streaming engine nulls the whole-file reference, run for run, on randomised material");
    {
        std::mt19937 rng (58);
        std::uniform_real_distribution<double> U (0.0, 1.0);
        int files = 0, equal = 0; long runs = 0;
        for (int it = 0; it < 60; ++it)
        {
            const double sr = std::array<double, 5> { 8000, 44100, 48000, 96000, 192000 }[(std::size_t) it % 5];
            const int nch = 1 + it % 3;
            const long N = (long) (sr * (0.3 + 1.2 * U (rng)));
            const int bits = std::array<int, 3> { 0, 16, 24 }[(std::size_t) (it / 3) % 3];
            Planes planes;
            for (int c = 0; c < nch; ++c)
            {
                std::vector<double> x;
                const int kind = (it + c) % 6;
                if (kind == 5)
                {   // pairs at fs/4 over a clamped programme — the most candidates the pending queue can be asked to hold
                    x.assign ((std::size_t) N, 0.0);
                    for (long n = 0; n < N; ++n) x[(std::size_t) n] = ((n / 2) % 2 ? -0.4 : 0.4);
                    auto y = programme (N, sr, (unsigned) it);
                    for (long n = N / 2; n < N; ++n) x[(std::size_t) n] = std::clamp (y[(std::size_t) n], -0.8, 0.8);
                }
                else
                {
                    x = programme (N, sr, (unsigned) (it * 7 + c));
                    normalisePeak (x, 1.0);
                    const double drive = kind == 0 ? 1.0 : std::pow (10.0, (0.5 + 12.0 * U (rng)) / 20.0);
                    for (double& v : x) v *= drive;
                    clampTruth (x, c, kind == 3 ? -0.7 : -1.0, kind == 4 ? 0.9 : 1.0);
                    const double att = std::pow (10.0, -30.0 * U (rng) / 20.0);
                    for (double& v : x) v *= att;
                }
                planes.push_back (deliver (x, bits, it % 4 == 1, (unsigned) it));
                if (it % 5 == 2)                                          // holes: NaN / ±Inf, one inside a flat top if there is one
                    for (int h = 0; h < 6; ++h)
                    {
                        const std::size_t at = (std::size_t) (U (rng) * (double) (N - 1));
                        planes.back()[at] = h % 3 == 0 ? std::numeric_limits<float>::quiet_NaN() : h % 3 == 1 ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity();
                    }
            }
            const auto ref = reference (planes, sr);
            const Report r = runEngine (planes, sr, -1, (unsigned) it);
            ++files; runs += (long) ref.size();
            if (sameRuns (ref, sortedRuns (r))) ++equal;
            else std::printf ("    reference null differs on item %d: reference %zu runs, engine %zu\n", it, ref.size(), r.runs.size());
        }
        char msg[128]; std::snprintf (msg, sizeof msg, "%d of %d items identical (%ld reference runs)", equal, files, runs);
        test::ok (equal == files && runs > 1000, msg);
    }

    test::group ("the reference null on tiny streams (W = 20 samples): window edges, stream edges, holes, dense ceilings");
    {
        // Short streams at 1 kHz put every run within reach of a window edge, a stream edge or a hole, and a dense
        // pattern of two-sample flats at a proven ceiling followed by an overshoot fills the pending queue to its bound.
        std::mt19937 rng (1158);
        std::uniform_real_distribution<double> U (0.0, 1.0);
        int items = 0, equal = 0; long refRuns = 0;
        const float nan = std::numeric_limits<float>::quiet_NaN(), inf = std::numeric_limits<float>::infinity();
        for (int it = 0; it < 4000; ++it)
        {
            const int nch = 1 + (int) (U (rng) * 2.0);
            const long N = 20 + (long) (U (rng) * 400.0);
            Planes planes;
            for (int c = 0; c < nch; ++c)
            {
                std::vector<double> x;
                double cur = U (rng) * 2.0 - 1.0;
                const double P = 0.3 + 0.6 * U (rng);
                while ((long) x.size() < N)
                {
                    const int kind = (int) (U (rng) * 8.0);
                    const int len = 1 + (int) (U (rng) * 10.0);
                    if (kind == 0) for (int i = 0; i < len; ++i) x.push_back (cur);                                 // flat where it is
                    else if (kind == 1) { cur = (U (rng) < 0.5 ? P : -P); for (int i = 0; i < len; ++i) x.push_back (cur); }   // flat at a ceiling
                    else if (kind == 2) { const double to = U (rng) * 2.0 - 1.0; for (int i = 1; i <= len; ++i) x.push_back (cur + (to - cur) * i / len); cur = to; }
                    else if (kind == 3) for (int i = 0; i < len; ++i) x.push_back (cur = U (rng) * 2.0 - 1.0);       // noise
                    else if (kind == 4) { x.push_back (cur - 0.1 * U (rng)); x.push_back (cur); }                     // a dip
                    else if (kind == 5) for (int i = 0; i < len; ++i) { x.push_back (P); x.push_back (P); x.push_back (P - 0.2 * U (rng)); }   // dense flats
                    else if (kind == 6) x.push_back (P + 0.05 + 0.2 * U (rng));                                       // an overshoot
                    else { const double to = U (rng) < 0.5 ? P : -P; for (int i = 1; i <= len; ++i) x.push_back (cur + (to - cur) * std::sin (1.5707963 * i / len)); cur = to; }
                }
                x.resize ((std::size_t) N);
                const int grid = (int) (U (rng) * 3.0);                                                          // float, 1/64 grid, 1/4096 grid
                if (grid) for (double& v : x) v = std::nearbyint (v * (grid == 1 ? 64.0 : 4096.0)) / (grid == 1 ? 64.0 : 4096.0);
                const double gain = U (rng) < 0.3 ? 0.37 : std::ldexp (1.0, -(int) (U (rng) * 4.0));
                std::vector<float> f ((std::size_t) N);
                for (long n = 0; n < N; ++n) f[(std::size_t) n] = (float) (x[(std::size_t) n] * gain);
                if (U (rng) < 0.3) { const std::size_t at = (std::size_t) (U (rng) * (double) N); f[at] = U (rng) < 0.5 ? nan : (U (rng) < 0.5 ? inf : -inf); }
                planes.push_back (f);
            }
            const auto ref = reference (planes, 1000.0);
            const Report r = runEngine (planes, 1000.0, it % 3 == 0 ? 1 : -1, (unsigned) it);
            ++items; refRuns += (long) ref.size();
            if (sameRuns (ref, sortedRuns (r))) ++equal;
            else if (items - equal <= 3) std::printf ("    tiny-stream null differs on item %d: reference %zu runs, engine %zu\n", it, ref.size(), r.runs.size());
        }
        char msg[128]; std::snprintf (msg, sizeof msg, "%d of %d tiny streams identical (%ld reference runs)", equal, items, refRuns);
        test::ok (equal == items && refRuns > 2000, msg);
    }

    test::group ("a neighbour equal to the band's lowest sample is flat with it, not below it");
    {
        // On a 1/64 grid (q = 1/64, tau = 2/64) at 1 kHz: a ramp of 12 proves a ceiling at 40, then the greedy scan ends a band
        // {36, 38} at 38 because 39 would widen it past tau, and starts {39, 38, 40, 40} — whose lowest sample IS its left
        // neighbour. That band sits at the ceiling, but it is no top: its left side is level with it.
        const double g = 1.0 / 64.0;
        std::vector<float> x;
        auto put = [&] (std::initializer_list<int> v) { for (int u : v) x.push_back ((float) (u * g)); };
        put ({ 0, 1, 0, -1, -10, 0, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 20, 10, 0, -10, 0 });   // the grid shown, a ramp run at 40
        for (int i = 0; i < 30; ++i) put ({ 0, -3 });
        put ({ 30, 36, 38, 39, 38, 40, 40, 30, 20, 10, 0, -10, 0 });
        for (int i = 0; i < 30; ++i) put ({ -3, 0 });
        const Report r = runEngine ({ x }, 1000.0, 1);
        const auto ref = reference ({ x }, 1000.0);
        test::ok (r.runs.size() == 1 && std::get<1> (r.runs[0]) == 9 && std::get<2> (r.runs[0]) == 12 && sameRuns (ref, sortedRuns (r)),
                  "only the proven run is reported (" + std::to_string (r.runs.size()) + " runs)");
    }

    test::group ("the pending queue holds every candidate a window can end: a full queue never decides early");
    {
        // W = 20 samples at 1 kHz, so W/2 + 1 two-sample bands can end inside one window. After a ramp proves a ceiling at
        // 40, pairs at 40 and 30 alternate — every pair at 40 sits at the ceiling — and then one sample at 50 passes it, within
        // the window of all of them. Eight bands pending is more than a queue of W/4 + 2 would hold.
        // Decided on their full window the pairs are rejected; decided early, before the 50 arrives, they would not be.
        const double g = 1.0 / 64.0;
        std::vector<float> x;
        auto put = [&] (std::initializer_list<int> v) { for (int u : v) x.push_back ((float) (u * g)); };
        put ({ 0, 1, 0, -1, -10, 0, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 20, 10, 0, -10, 0 });
        for (int i = 0; i < 15; ++i) put ({ 0, -3 });
        for (int i = 0; i < 4; ++i) put ({ 40, 40, 30, 30 });                                        // eight bands end in 16 samples
        put ({ 20, 50, 10, 0 });                                                                     // passed within the window of every one
        for (int i = 0; i < 30; ++i) put ({ -3, 0 });
        const Report r = runEngine ({ x }, 1000.0, 1);
        test::ok (r.runs.size() == 1 && std::get<1> (r.runs[0]) == 9 && sameRuns (reference ({ x }, 1000.0), sortedRuns (r)),
                  "only the proven run is reported (" + std::to_string (r.runs.size()) + " runs)");
    }

    test::group ("a stream edge proves nothing: a signal cut into digital silence at the start or the end is not a clamp");
    {
        for (int bits : { 0, 16 })
        {
            const double sr = 48000; const long N = 48000;
            std::vector<double> x ((std::size_t) N, 0.0);
            for (long n = 0; n < N / 2; ++n) x[(std::size_t) n] = 0.8 * std::sin (kTau2Pi * 30.0 * n / sr + 0.3);   // ends mid-swing
            test::ok (! runEngine ({ deliver (x, bits) }, sr, 0).clipped, "cut into trailing silence, bits " + std::to_string (bits));
            std::reverse (x.begin(), x.end());
            test::ok (! runEngine ({ deliver (x, bits) }, sr, 0).clipped, "leading silence cut into a signal, bits " + std::to_string (bits));
        }
    }

    // ---------------------------------------------------------------------------------------------- invariances
    test::group ("invariances: slicing, polarity, channel order, power-of-two level");
    {
        const double sr = 44100;
        Planes st;
        std::vector<Truth> truth;
        for (int c = 0; c < 2; ++c)
        {
            auto x = programme (132300, sr, 40u + (unsigned) c);
            normalisePeak (x, 1.0);
            for (double& v : x) v *= 2.2;
            clampTruth (x, c, -1.0, 1.0);
            for (double& v : x) v *= 0.25;
            st.push_back (deliver (x, 16));
        }
        const Report whole = runEngine (st, sr, 0);
        test::ok (whole.clipped && whole.runs.size() > 100, "the fixture is clipped");
        test::ok (sameReport (whole, runEngine (st, sr, 1)), "one-sample calls: the identical report (run order, peak and DC bits)");
        test::ok (sameReport (whole, runEngine (st, sr, -1, 3)), "random call lengths: the identical report");
        test::ok (sameReport (whole, runEngine (st, sr, 64)), "64-sample calls: the identical report");

        Planes inv = st; for (auto& p : inv) for (float& v : p) v = -v;
        const Report ri = runEngine (inv, sr, 0);
        bool polOk = ri.runs.size() == whole.runs.size();
        for (std::size_t i = 0; polOk && i < ri.runs.size(); ++i)
            polOk = std::get<1> (ri.runs[i]) == std::get<1> (whole.runs[i]) && std::get<2> (ri.runs[i]) == std::get<2> (whole.runs[i])
                 && std::get<3> (ri.runs[i]) == -std::get<3> (whole.runs[i]) && std::get<4> (ri.runs[i]) == -std::get<4> (whole.runs[i]);
        test::ok (polOk, "polarity inverted: the same runs, sign and level negated");

        Planes sw { st[1], st[0] };
        const Report rs = runEngine (sw, sr, 0);
        auto a = sortedRuns (whole), b = sortedRuns (rs);
        for (auto& u : b) u.ch = 1 - u.ch;
        std::sort (b.begin(), b.end(), [] (const RefRun& p, const RefRun& q) { return std::tie (p.ch, p.start) < std::tie (q.ch, q.start); });
        test::ok (sameRuns (a, b), "channels swapped: the same runs on the swapped channels");

        for (int k : { 1, 3, 6 })
        {
            Planes sc = st; for (auto& p : sc) for (float& v : p) v = std::ldexp (v, -k);
            const Report rk = runEngine (sc, sr, 0);
            test::ok (sameRuns (sortedRuns (whole), sortedRuns (rk)), "turned down by 2^-" + std::to_string (k) + " (exact): the same runs");
        }
    }

    // ---------------------------------------------------------------------------------------------- non-finite
    test::group ("non-finite samples neither crash nor move a verdict, in either direction");
    {
        const double sr = 48000; const long N = 96000;
        auto x = programme (N, sr, 77);
        normalisePeak (x, 0.9);
        Planes clean { deliver (x, 16) };
        auto y = x; for (double& v : y) v *= 2.5;
        const auto truth = clampTruth (y, 0, -0.9, 0.9);
        Planes clipped { deliver (y, 16) };
        long inRun = -1;
        // a run from the second half: its halves have one unknown side each, so they are found as CEILING runs, which needs
        // the ceiling proven earlier in the stream
        for (const Truth& t : truth) if (t.len >= 9 && t.start > N / 2) { inRun = t.start + t.len / 2; break; }
        test::ok (inRun > 0, "the fixture has a run of >= 9 in its second half to put a NaN into");
        const float nan = std::numeric_limits<float>::quiet_NaN(), inf = std::numeric_limits<float>::infinity();
        const long marks[] = { 0, 5, 1000, 20000, 40000, N - 1 };
        for (auto* p : { &clean, &clipped })
        {
            int m = 0;
            for (long at : marks) (*p)[0][(std::size_t) at] = m++ % 3 == 0 ? nan : (m % 3 == 1 ? inf : -inf);
        }
        clipped[0][(std::size_t) inRun] = nan;
        const Report rc = runEngine (clean, sr, 333), rk = runEngine (clipped, sr, 333);
        test::ok (! rc.clipped, "clean + NaN / ±Inf: still clean");
        test::ok (rk.clipped, "clipped + NaN / ±Inf (one inside a run): still clipped");
        test::ok (rc.nonFinite[0] == 6 && rk.nonFinite[0] == 7, "non-finite samples are counted");
        test::ok (std::isfinite (rc.peak[0]) && rc.peak[0] <= 0.9 + 1e-6 && std::isfinite (rc.dc[0]), "peak and DC ignore them");
        bool bothHalves = false, spans = false;
        for (const auto& u : rk.runs)
        {
            if (std::get<1> (u) + std::get<2> (u) == inRun) for (const auto& w : rk.runs) if (std::get<1> (w) == inRun + 1) bothHalves = true;
            if (std::get<1> (u) <= inRun && inRun < std::get<1> (u) + std::get<2> (u)) spans = true;
        }
        test::ok (bothHalves && ! spans, "a NaN inside a flat top splits it: both halves are reported, no run covers the hole");
        Planes allNan { std::vector<float> ((std::size_t) N, nan) };
        const Report rn = runEngine (allNan, sr, 0);
        test::ok (! rn.clipped && rn.nonFinite[0] == N && rn.peak[0] == 0.0, "an all-NaN stream: nothing, and every sample counted");
    }

    // ---------------------------------------------------------------------------------------------- report fields
    test::group ("the report: peak, DC, per-channel counts, capacity");
    {
        const double sr = 48000; const long N = 48000;
        std::vector<double> x ((std::size_t) N);
        for (long n = 0; n < N; ++n) x[(std::size_t) n] = 0.125 + 0.5 * std::sin (kTau2Pi * 100.0 * n / sr);   // 480 samples a period: whole periods
        auto y = x; for (double& v : y) v *= 1.6;
        const auto truth = clampTruth (y, 1, -0.5, 0.5);
        const Report r = runEngine ({ deliver (x, 0), deliver (y, 0) }, sr, 0);
        test::approx (r.peak[0], 0.625, 1e-6, "sample peak of channel 0");
        test::approx (r.dc[0], 0.125, 1e-6, "DC of channel 0 (whole periods)");
        test::ok (r.perChRuns[0] == 0 && r.perChRuns[1] > 0 && ! r.runs.empty(), "only channel 1 carries runs");
        long sumLen = 0, longest = 0;
        for (const auto& u : r.runs) { sumLen += std::get<2> (u); longest = std::max (longest, std::get<2> (u)); test::ok (std::get<0> (u) == 1, "run channel"); }
        test::ok (r.clippedSamples[1] == sumLen && r.longest[1] == longest, "clippedSamples and longestRun agree with the run list");
        bool levels = true;
        // the level is the band's mean: a slow sine's neighbours within 2q of the ceiling belong to the band (measured 0.4999612)
        for (const auto& u : r.runs) levels = levels && std::fabs (std::get<4> (u) - (std::get<3> (u) > 0 ? 0.5 : -0.5)) < 1e-4;
        test::ok (levels, "each run's level is its ceiling, with its sign");

        bool rampFirst[2] = { false, false }, seen[2] = { false, false }, ceilingsAtARamp = true; int nCeil = 0;
        std::vector<double> rampLevels[2];
        for (const auto& u : r.runs)
        {
            const int pol = std::get<3> (u) > 0 ? 0 : 1;
            const bool isRamp = std::get<5> (u) == (int) analysis::ClipEvidence::Ramp;
            if (! seen[pol]) { seen[pol] = true; rampFirst[pol] = isRamp; }
            if (isRamp) rampLevels[pol].push_back (std::get<4> (u));
            else
            {
                ++nCeil;
                bool near = false;
                for (double lv : rampLevels[pol]) near = near || std::fabs (lv - std::get<4> (u)) < 1e-3;
                ceilingsAtARamp = ceilingsAtARamp && near;
            }
        }
        {   // a programme's clamp: short runs between dips come in as Ceiling evidence
            auto pg = programme (96000, sr, 12); normalisePeak (pg, 1.0); for (double& v : pg) v *= 2.0; clampTruth (pg, 0, -1, 1);
            const Report rp = runEngine ({ deliver (pg, 0) }, sr, 0);
            int ramps = 0, ceils = 0;
            for (const auto& u : rp.runs) (std::get<5> (u) == (int) analysis::ClipEvidence::Ramp ? ramps : ceils)++;
            test::ok (ramps > 0 && ceils > 0, "a clamped programme carries both kinds of evidence (" + std::to_string (ramps) + " ramp, " + std::to_string (ceils) + " ceiling)");
        }
        // (only the tops are found on this fixture: a pure tone locked to the grid steps onto its floor by 0.002 every period, which
        // is then the smallest step the stream ever shows — q is that loose, and V(L) with it. Music has slow places; a tone does not.)
        test::ok (seen[0] && (! seen[0] || rampFirst[0]) && (! seen[1] || rampFirst[1]), "the first run of each polarity is Ramp evidence — a ceiling is never assumed");
        test::ok (ceilingsAtARamp, "every Ceiling run sits at the level of an earlier Ramp run of its polarity (" + std::to_string (nCeil) + " of them)");

        {   // DC is summed with compensation: {1e20, 1, -1e20} averages 1/3, which a plain double sum loses entirely
            std::vector<float> big (30000);
            for (std::size_t i = 0; i < big.size(); ++i) big[i] = i % 3 == 0 ? 1.0e20f : i % 3 == 1 ? 1.0f : -1.0e20f;
            test::approx (runEngine ({ big }, sr, 1000).dc[0], 1.0 / 3.0, 1e-9, "DC with compensated summation");
            std::vector<float> early ((std::size_t) N, 0.1f); early[10] = -0.75f;
            const Report re = runEngine ({ early }, sr, 777);
            test::ok (re.peak[0] == 0.75, "a single early peak is kept for the whole stream");
        }

        const Report capped = runEngine ({ deliver (x, 0), deliver (y, 0) }, sr, 0, 7, 3);
        test::ok (capped.runs.size() == 3 && ! capped.complete && capped.perChRuns[1] == r.perChRuns[1], "capacity 3: three stored, the count keeps counting, runsComplete() is false");
        const Report exact = runEngine ({ deliver (x, 0), deliver (y, 0) }, sr, 0, 7, (int) r.runs.size());
        test::ok (exact.complete && exact.runs.size() == r.runs.size(), "capacity exactly the run count: complete");
    }

    // ---------------------------------------------------------------------------------------------- law 11 contract
    test::group ("every report accessor is total: out-of-range indices and an unprepared object answer empty");
    {
        CD u;
        test::ok (u.samplePeak (0) == 0.0 && u.dcOffset (0) == 0.0 && u.runCount (0) == 0 && u.longestRun (3) == 0
                  && u.run (0).length == 0 && u.samplePeak() == 0.0, "before prepare(): zeros, no read past an empty vector");
        const double sr = 48000; const long N = 48000;
        std::vector<double> x ((std::size_t) N);
        for (long n = 0; n < N; ++n) x[(std::size_t) n] = std::clamp (1.6 * std::sin (kTau2Pi * 100.0 * n / sr + 0.3), -0.5, 0.5);
        CD d; d.setParams ({ 2 });
        test::run (d.prepare (sr, 512, 2));
        const auto y = deliver (x, 16);
        const float* io[2] { y.data(), y.data() };
        test::run (d.process (io, 2, (int) N));
        d.finish();
        test::ok (d.runCount() > 2 && d.storedRunCount() == 2, "the fixture overflows a two-run list");
        test::ok (d.run (2).length == 0 && d.run (-1).length == 0 && d.run (d.runCount()).length == 0 && d.run (1).length > 0,
                  "run(i) past the stored list, or negative: an empty run");
        test::ok (d.samplePeak (2) == 0.0 && d.samplePeak (-1) == 0.0 && d.nonFiniteSamples (16) == 0 && d.clippedSamples (-5) == 0,
                  "a channel outside the prepared width: zeros");
    }

    test::group ("law 11: refusals, clock-only calls, stopped channels, finish, reset, prepare");
    {
        CD d;
        const float z[4] {}; const float* io[2] { z, z };
        test::ok (! d.process (io, 1, 4), "unprepared: refused");
        test::ok (! d.prepare (0.0, 64, 2) && ! d.prepare (std::numeric_limits<double>::quiet_NaN(), 64, 2) && ! d.prepare (999.0, 64, 2)
                  && ! d.prepare (768001.0, 64, 2) && ! d.prepare (48000, 64, 0) && ! d.prepare (48000, 64, core::kMaxChannels + 1), "prepare refuses what it cannot honour");
        test::run (d.prepare (48000, 64, 2));
        d.setParams ({ -1 });
        test::ok (! d.prepare (48000, 64, 2) && ! d.process (io, 1, 4), "a refused prepare leaves the object unusable, not on its previous build");
        d.setParams ({});
        test::run (d.prepare (48000, 64, 2));
        test::ok (! d.process (io, -1, 4) && ! d.process (io, 1, -1), "malformed: refused");
        test::ok (! d.process (io, 3, 4), "wider than prepared: refused");
        test::ok (d.samplesProcessed() == 0, "a refused call moves no clock");
        test::ok (d.process (io, 2, 0) && d.samplesProcessed() == 0, "n == 0: accepted, no time");
        test::ok (d.process (nullptr, 0, 480) && d.samplesProcessed() == 480, "nch == 0: a clock-only call advances time");
        test::ok (d.nonFiniteSamples (0) == 0, "a stopped channel's missing samples are holes, not non-finite input");
        d.finish();
        test::ok (d.isFinished() && ! d.process (io, 2, 4), "after finish(): refused");
        d.finish();
        test::ok (d.isFinished(), "finish() is idempotent");
        d.reset();
        test::ok (! d.isFinished() && d.samplesProcessed() == 0 && d.process (io, 2, 4), "reset() re-arms");

        // a channel that stops mid-run and comes back: its gap is a hole, the other channel is untouched
        const double sr = 48000; const long N = 96000;
        auto x = programme (N, sr, 3); normalisePeak (x, 1.0); for (double& v : x) v *= 2.0;
        const auto truthGap = clampTruth (x, 0, -1, 1);   // no-op (already clamped): the runs, to cut the gap through one
        long gapAt = 30000;
        for (const Truth& tr : truthGap) if (tr.start > 25000 && tr.len >= 8) { gapAt = tr.start + tr.len / 2; break; }
        const auto a = deliver (x, 0);
        const Report solo = runEngine ({ a }, sr, 0);
        CD e; test::run (e.prepare (sr, 512, 2));
        const float* pp[2];
        for (long pos = 0; pos < N;)
        {   // calls cut at the gap's edges; the gap starts inside a flat top of channel 1
            const long next = pos < gapAt ? gapAt : pos < gapAt + 10000 ? gapAt + 10000 : N;
            const int len = (int) std::min (1000L, next - pos);
            pp[0] = a.data() + pos; pp[1] = a.data() + pos;
            const bool gap = pos >= gapAt && pos < gapAt + 10000;
            test::run (e.process (pp, gap ? 1 : 2, len));
            pos += len;
        }
        e.finish();
        std::vector<RefRun> ch0, ch1;
        for (std::int64_t i = 0; i < e.storedRunCount(); ++i) (e.run (i).channel == 0 ? ch0 : ch1).push_back ({ 0, (long) e.run (i).start, (long) e.run (i).length });
        test::ok (sameRuns (ch0, sortedRuns (solo)), "channel 0, fed throughout: exactly the solo result");
        test::ok (gapAt != 30000, "the gap starts inside a flat top");
        Planes gapped { a }; for (long n = gapAt; n < gapAt + 10000; ++n) gapped[0][(std::size_t) n] = std::numeric_limits<float>::quiet_NaN();
        auto ref = reference (gapped, sr);
        test::ok (sameRuns (ch1, ref), "channel 1, stopped for 10000 samples: exactly the reference with a hole there");
    }

    // ---------------------------------------------------------------------------------------------- RT + storage
    test::group ("RT: process(), finish() and reset() allocate nothing; prepare() allocates exactly storageFor()");
    {
        const double sr = 96000; const long N = 192000;
        auto x = programme (N, sr, 21); normalisePeak (x, 1.0); for (double& v : x) v *= 3.0; clampTruth (x, 0, -1, 1);
        const auto a = deliver (x, 16);
        CD d; d.setParams ({ 5000 });
        const long long b0 = alloc::bytes.load();
        const long long al0 = alloc::count.load();
        test::run (d.prepare (sr, 512, 3));
        // both deltas are read into locals BEFORE the check: the message is a std::string that allocates, and gcc
        // builds a call's arguments in its own order — reading the counter inside the call counts the message too
        const long long used = alloc::bytes.load() - b0;
        const long long blocks = alloc::count.load() - al0;
        const auto st = CD::storageFor (sr, 3, 5000);
        test::ok (st.ok && (long long) st.bytes() == used && blocks == 4, "prepare() requested exactly storageFor().bytes() in four blocks ("
                  + std::to_string (used) + " bytes, " + std::to_string (blocks) + " blocks)");
        test::ok (! CD::storageFor (500.0, 3, 10).ok && ! CD::storageFor (sr, 0, 10).ok && ! CD::storageFor (sr, 3, -1).ok && CD::storageFor (sr, 3, 0).ok,
                  "storageFor() refuses exactly what prepare() refuses");
        const float* io[3] { a.data(), a.data(), a.data() };
        const long long before = alloc::count.load();
        for (long pos = 0; pos + 777 <= N; pos += 777) { io[0] = io[1] = io[2] = a.data() + pos; (void) d.process (io, 3, 777); }
        (void) d.process (nullptr, 0, 5000);
        d.finish();
        d.reset();
        const bool noAlloc = alloc::count.load() == before;   // read before the message string exists
        test::okNoAlloc (noAlloc, "no allocation in process / clock-only / finish / reset");
    }

    // ---------------------------------------------------------------------------------------------- what it does not see
    test::group ("what it does not see: a clamp disturbed afterwards by more than the band is no longer flat");
    {
        // NB a short smoother is NOT such a disturbance — [1/4 1/2 1/4] leaves the inside of every long flat top exactly at
        // the ceiling, and that is found. Noise wider than tau is; so, measured, are AAC, MP3 and resampling.
        const double sr = 48000; const long N = 144000;
        auto x = programme (N, sr, 8); normalisePeak (x, 1.0); for (double& v : x) v *= 2.0;
        clampTruth (x, 0, -1, 1);
        for (double& v : x) v *= 0.5;
        auto withNoise = [&] (double lsb) {
            std::mt19937 rng (4); std::uniform_real_distribution<double> U (-lsb / 32768.0, lsb / 32768.0);
            std::vector<double> y = x; for (double& v : y) v += U (rng);
            return deliver (y, 16);
        };
        // measured on this fixture, runs of >= 3 found: ±0 LSB 0.991, ±1 0.964, ±2 0.495, ±4 and more 0.000
        test::ok (runEngine ({ withNoise (1.0) }, sr, 0).clipped, "the clamp with ±1 LSB of noise after it (inside the band): found");
        test::ok (! runEngine ({ withNoise (32.0) }, sr, 0).clipped, "the same clamp with ±32 LSB of noise after it: not found (the limit is published)");
    }

    test::group ("one PCM, two histories: a one-sample pulse train on a flat level, entered by a moving signal");
    {
        // After slow music, a flat level at -0.7 carrying one-sample pulses to +0.7 is bit for bit a signal clamped at -0.7
        // whose excursions inside the range last one sample. The first flat stretch is entered by the music moving and
        // left by a one-sample dip — the evidence a clamp leaves. Reported, and pinned so a change to it is a decision.
        // (The same pulses with nothing moving into the level — from rest, or from silence — are not: see the clean group.)
        const double sr = 48000; const long N = 48000;
        std::vector<double> pl ((std::size_t) N);
        for (long n = 0; n < N; ++n) pl[(std::size_t) n] = n < N / 2 ? 0.3 * std::sin (kTau2Pi * 40.0 * n / sr) : (n % 9 == 0 ? 0.7 : -0.7);
        test::ok (runEngine ({ deliver (pl, 16) }, sr, 0).clipped, "reported as clipped");
        std::vector<double> rest ((std::size_t) N);
        for (long n = 0; n < N; ++n) rest[(std::size_t) n] = n % 9 == 0 ? 0.7 : -0.7;
        test::ok (! runEngine ({ deliver (rest, 16) }, sr, 0).clipped, "the same pulse train with no moving entry: not reported");
    }

    return test::report();
}
