// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The Transformer shape (Saturator.h, WaveShaper.h) — the one shape with a model state. What each group pins:
//
//   * REFERENCE NULL at os 1: the stage against the model's three lines written here in double, on an LF-rich
//     full-scale signal, at four drives and two rates.
//   * Its behaviour at os 4: THD against the same double model run at 4·fs, THD falling ~4x per octave, and the top
//     octaves clean where Tanh at the same drive is not. The same THD at 44.1, 48 and 96 kHz.
//   * Transparency: a -60 dBFS multitone comes out as Tanh at driveDb 0 lets it through; mix 0.5 does not comb.
//   * The peak law |y| <= max|x| at os 1, and at os 4 an overshoot no worse than Tanh's.
//   * slopeAtZero() is the literal 1.0f, so drive-compensation is exactly 1.0f.
//   * The model state's lifecycle: NaN, the denormal flush, chunking, reset, a channel that leaves and returns, a
//     switch away from the shape and back, its bytes, and the glide landing on the settled floats.

#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter
#include <felitronics/saturation/Saturator.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
using Sat = saturation::Saturator;
using WS = saturation::WaveShaper;
using Shape = WS::Shape;
using Buf = std::vector<std::vector<float>>;

namespace
{
constexpr double kPi = 3.14159265358979323846;

std::string num (double v, const char* fmt = "%.4g") { char b[64]; std::snprintf (b, sizeof b, fmt, v); return b; }

void okp (bool cond, const std::string& msg) { std::printf ("       %s\n", msg.c_str()); ok (cond, msg); }

// A running maximum a NaN cannot hide from: a NaN wins and STAYS, and every gate compares with <=, which it fails.
template <class T> void worse (T& w, T v) { if (! std::isnan (w) && ! (v <= w)) w = v; }

bool bitEq (float a, float b) { return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b); }

long long diffs (const std::vector<float>& a, const std::vector<float>& b, std::size_t from, std::size_t to)
{
    long long d = 0;
    for (std::size_t i = from; i < to; ++i) d += bitEq (a[i], b[i]) ? 0 : 1;
    return d;
}

Sat::Params P (Shape s, float driveDb, float mix = 1.0f, float autoComp = 1.0f, float outDb = 0.0f)
{
    Sat::Params p; p.shape = s; p.driveDb = driveDb; p.mix = mix; p.autoComp = autoComp; p.outputDb = outDb;
    return p;
}

bool prep (Sat& s, double fs, int nch, int os, int maxBlock = 512)
{
    const bool okP = s.prepare (fs, maxBlock, nch, os);
    if (! okP) ok (false, "PRECONDITION: prepare");
    return okP;
}

// `y` through a prepared stage, in calls of `block` samples on the first `nch` channels.
void through (Sat& s, Buf& y, int block, std::size_t from = 0, std::size_t to = 0, int nch = -1)
{
    if (to == 0) to = y[0].size();
    if (nch < 0) nch = (int) y.size();
    for (std::size_t o = from; o < to; )
    {
        const int m = (int) std::min<std::size_t> ((std::size_t) block, to - o);
        float* io[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) io[c] = y[(std::size_t) c].data() + o;
        felitronics::test::run (s.process (io, nch, m));
        o += (std::size_t) m;
    }
}

std::vector<float> render (const Sat::Params& p, double fs, int os, std::vector<float> x, int block = 512)
{
    Sat s; s.setParams (p);
    if (! prep (s, fs, 1, os)) return {};
    Buf y { std::move (x) };
    through (s, y, block);
    return y[0];
}

//==============================================================================
// Signals.
std::vector<float> sine (double fs, double hz, double amp, std::size_t n)
{
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = (float) (amp * std::sin (2.0 * kPi * hz * (double) i / fs));
    return x;
}

// 20..200 Hz tones, a 30 Hz square and a DC step, scaled to a peak of 1. `v` varies the phases.
std::vector<float> lfRich (double fs, std::size_t n, int v = 0)
{
    std::vector<double> d (n);
    double peak = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / fs, ph = 0.7 * v;
        const double sq = std::sin (2.0 * kPi * 30.0 * t + ph) >= 0.0 ? 1.0 : -1.0;
        d[i] = 0.22 * std::sin (2.0 * kPi * 20.0 * t + ph) + 0.18 * std::sin (2.0 * kPi * 47.0 * t + 0.3 + ph)
             + 0.15 * std::sin (2.0 * kPi * 110.0 * t + 1.1) + 0.12 * std::sin (2.0 * kPi * 200.0 * t + 2.0 + ph)
             + 0.2 * sq + (i >= n * 2 / 5 ? 0.13 : -0.13);
        peak = std::max (peak, std::fabs (d[i]));
    }
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = (float) (d[i] / peak);
    return x;
}

// Tones from `lo` to `hi` Hz — `lo`·1.6^j below `hi`, and `hi` itself — scaled to a peak of `peak`.
std::vector<float> multitone (double fs, std::size_t n, double lo, double hi, double peak)
{
    std::vector<double> fr;
    for (double f = lo; f < hi; f *= 1.6) fr.push_back (f);
    fr.push_back (hi);
    std::vector<double> d (n, 0.0);
    for (std::size_t j = 0; j < fr.size(); ++j)
        for (std::size_t i = 0; i < n; ++i) d[i] += std::sin (2.0 * kPi * fr[j] * (double) i / fs + 0.9 * (double) j);
    double m = 0.0;
    for (double v : d) m = std::max (m, std::fabs (v));
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = (float) (peak * d[i] / m);
    return x;
}

//==============================================================================
// THE ORACLE: the model written from its definition, in double, touching no production code.
//   L += a·(x − L),  a = 1 − exp(−2π·40/fs);   t = tanh(k·L);   y = t/k + (1 − t²)·(x − L)
// k = 10^(driveDb/20) − 1 with the WaveShaper's documented drive floor 1e-4.
std::vector<double> oracle (const std::vector<double>& x, double fs, double driveDb)
{
    const double a = 1.0 - std::exp (-2.0 * kPi * 40.0 / fs);
    const double k = std::max (1.0e-4, std::pow (10.0, driveDb / 20.0) - 1.0);
    double L = 0.0;
    std::vector<double> y (x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        L += a * (x[i] - L);
        const double t = std::tanh (k * L);
        y[i] = t / k + (1.0 - t * t) * (x[i] - L);
    }
    return y;
}

// Amplitude of the component at `cycles` cycles per window y[0, n) (Goertzel, double).
template <class T>
double amp (const T* y, std::size_t n, double cycles)
{
    const double w = 2.0 * kPi * cycles / (double) n, c2 = 2.0 * std::cos (w);
    double s1 = 0.0, s2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) { const double s0 = (double) y[i] + c2 * s1 - s2; s2 = s1; s1 = s0; }
    const double pw = s1 * s1 + s2 * s2 - c2 * s1 * s2;
    return 2.0 * std::sqrt (std::max (0.0, pw)) / (double) n;
}

// THD of a bin-exact tone at `cycles` cycles per window: every harmonic below `nyqCycles`, over the fundamental.
template <class T>
double thd (const T* y, std::size_t n, double cycles, double nyqCycles)
{
    double h = 0.0;
    for (int k = 2; (double) k * cycles < nyqCycles; ++k) { const double a = amp (y, n, (double) k * cycles); h += a * a; }
    return std::sqrt (h) / amp (y, n, cycles);
}

// A full-scale sine at `hz` through the stage at os 4 (1 s window after 0.5 s): its THD below the base Nyquist.
double stageThd (Shape s, double fs, double hz, float driveDb)
{
    const std::size_t warm = (std::size_t) (fs / 2.0), n = (std::size_t) fs;
    const std::vector<float> y = render (P (s, driveDb), fs, 4, sine (fs, hz, 1.0, warm + n));
    return thd (y.data() + warm, n, hz, fs / 2.0);
}

// The same sine through the double model at 4·fs, with the same window in seconds and the same harmonics.
double oracleThd (double fs, double hz, double driveDb)
{
    const double f4 = 4.0 * fs;
    const std::size_t warm = (std::size_t) (f4 / 2.0), n = (std::size_t) f4;
    std::vector<double> x (warm + n);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::sin (2.0 * kPi * hz * (double) i / f4);
    const std::vector<double> y = oracle (x, f4, driveDb);
    return thd (y.data() + warm, n, hz, fs / 2.0);
}

//==============================================================================
void referenceNull()
{
    double worstAll = 0.0;
    for (double fs : { 44100.0, 96000.0 })
    {
        const std::vector<float> x = lfRich (fs, (std::size_t) (1.5 * fs));
        const std::vector<double> xd (x.begin(), x.end());
        for (float db : { 0.0f, 6.0f, 12.0f, 24.0f })
        {
            const std::vector<float> y = render (P (Shape::Transformer, db), fs, 1, x);
            const std::vector<double> r = oracle (xd, fs, (double) db);
            double worst = y.size() == r.size() ? 0.0 : std::numeric_limits<double>::quiet_NaN();
            for (std::size_t i = 0; i < y.size() && i < r.size(); ++i) worse (worst, std::fabs ((double) y[i] - r[i]));
            std::printf ("       %5.1f kHz, driveDb %4.1f: max |stage - double model| %.3g\n", fs / 1000.0, (double) db, worst);
            worse (worstAll, worst);
        }
    }
    okp (worstAll <= 1.0e-5, "os 1, mix 1, autoComp 1: the stage nulls against the double model to <= 1e-5 at driveDb 0/6/12/24, "
                             "44.1 and 96 kHz (worst " + num (worstAll, "%.3g") + ")");
}

void behaviourAtOs4()
{
    const double fs = 48000.0;
    std::printf ("       full-scale sine, driveDb 6, os 4 at 48 kHz: THD (%%) stage | double model at 4*fs\n");
    bool within = true;
    double t320 = 0.0, t640 = 0.0;
    for (double hz : { 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 5000.0 })
    {
        const double st = stageThd (Shape::Transformer, fs, hz, 6.0f), orc = oracleThd (fs, hz, 6.0);
        std::printf ("       %6.0f Hz: %9.4f | %9.4f\n", hz, 100.0 * st, 100.0 * orc);
        if (hz == 40.0 || hz == 160.0 || hz == 640.0) within = within && std::fabs (st / orc - 1.0) <= 0.1;
        if (hz == 320.0) t320 = st;
        if (hz == 640.0) t640 = st;
    }
    okp (within, "THD at 40, 160 and 640 Hz within 10 % of the double model's own");
    const double ratio = t320 / t640;
    okp (ratio >= 3.4 && ratio <= 4.6, "THD(320) / THD(640) in [3.4, 4.6] (" + num (ratio) + ")");
    const double tr5 = stageThd (Shape::Transformer, fs, 5000.0, 6.0f), ta5 = stageThd (Shape::Tanh, fs, 5000.0, 6.0f);
    okp (tr5 < 1.0e-4 && ta5 > 1.0e-2, "5 kHz: Transformer THD " + num (100.0 * tr5) + " % < 0.01 %, Tanh's at the same drive "
                                        + num (100.0 * ta5) + " % > 1 %");
}

void rateIndependence()
{
    double lo = 1e9, hi = -1e9;
    for (double fs : { 44100.0, 48000.0, 96000.0 })
    {
        const double db = 20.0 * std::log10 (stageThd (Shape::Transformer, fs, 80.0, 6.0f));
        std::printf ("       %5.1f kHz: THD of a full-scale 80 Hz sine, driveDb 6, os 4: %.3f dB\n", fs / 1000.0, db);
        if (std::isnan (db) || db < lo) lo = db;
        worse (hi, db);
    }
    okp (hi - lo <= 0.1, "80 Hz THD at 44.1 / 48 / 96 kHz within 0.1 dB (spread " + num (hi - lo, "%.4f") + " dB)");
}

void smallSignal()
{
    const double fs = 48000.0;
    const std::vector<float> x = multitone (fs, 48000, 20.0, 20000.0, 1.0e-3);
    const std::vector<float> a = render (P (Shape::Transformer, 6.0f), fs, 4, x);
    const std::vector<float> b = render (P (Shape::Tanh, 0.0f), fs, 4, x);
    double worst = a.size() == b.size() && ! a.empty() ? 0.0 : std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) worse (worst, std::fabs ((double) a[i] - (double) b[i]));
    okp (worst <= 1.0e-6, "-60 dBFS multitone 20 Hz..20 kHz, os 4: Transformer at driveDb 6 against Tanh at driveDb 0, max |diff| "
                          + num (worst, "%.3g") + " <= 1e-6");
}

void mixAndLatency()
{
    bool same = true;
    for (int os : { 1, 2, 4, 8 })
    {
        Sat t, h;
        t.setParams (P (Shape::Transformer, 6.0f)); h.setParams (P (Shape::Tanh, 6.0f));
        same = same && prep (t, 48000.0, 2, os) && prep (h, 48000.0, 2, os) && t.latencySamples() == h.latencySamples();
    }
    ok (same, "latencySamples() equals Tanh's at os 1, 2, 4 and 8");
    // In the oversampler's passband, where its round trip is a pure delay to far below 1e-6 at this level.
    const double fs = 48000.0;
    const std::vector<float> x = multitone (fs, 48000, 20.0, 12000.0, 1.0e-3);
    Sat s; s.setParams (P (Shape::Transformer, 6.0f, 0.5f));
    if (! prep (s, fs, 1, 4)) return;
    const std::size_t lat = (std::size_t) s.latencySamples();
    Buf y { x };
    through (s, y, 512);
    double worst = lat > 0 ? 0.0 : std::numeric_limits<double>::quiet_NaN(), comb = 0.0;
    for (std::size_t i = lat + 2048; i < x.size(); ++i)
    {
        worse (worst, std::fabs ((double) y[0][i] - (double) x[i - lat]));
        worse (comb, std::fabs ((double) y[0][i] - (double) x[i]));
    }
    okp (worst <= 1.0e-6, "mix 0.5, -60 dBFS multitone 20 Hz..12 kHz, os 4: nulls against the input delayed by the latency ("
                          + std::to_string (lat) + " samples) to " + num (worst, "%.3g") + " (against the undelayed input: "
                          + num (comb, "%.3g") + ")");
}

//==============================================================================
void peakLaw()
{
    const double fs = 48000.0;
    const std::size_t n = (std::size_t) (1.5 * fs);
    struct M { const char* name; std::vector<float> x; };
    std::vector<M> ms;
    {
        std::vector<float> a (n), b (n), c (n), d (n);
        const float steps[] { 1.0f, -1.0f, 0.5f, -1.0f, 1.0f, 0.0f };
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double) i / fs, tone = 0.8 * std::sin (2.0 * kPi * 25.0 * t);
            a[i] = std::sin (2.0 * kPi * 20.0 * t) >= 0.0 ? 1.0f : -1.0f;
            b[i] = std::sin (2.0 * kPi * 55.0 * t) >= 0.0 ? 0.9f : -0.9f;
            c[i] = steps[std::min<std::size_t> (5, i / (n / 6))];
            d[i] = (float) ((i % 1920) < 96 ? (tone >= 0.0 ? -1.0 : 1.0) : tone);   // 2 ms bursts against the tone
        }
        ms = { { "20 Hz square", a }, { "55 Hz square 0.9", b }, { "DC steps", c }, { "25 Hz tone + opposite bursts", d } };
    }
    double worstRatio = 0.0;
    for (const M& m : ms)
    {
        double in = 0.0;
        for (float v : m.x) worse (in, (double) std::fabs (v));
        for (float db : { 6.0f, 12.0f, 18.0f, 24.0f, 30.0f, 36.0f })
        {
            const std::vector<float> y = render (P (Shape::Transformer, db), fs, 1, m.x);
            double out = y.empty() ? std::numeric_limits<double>::quiet_NaN() : 0.0;
            for (float v : y) worse (out, (double) std::fabs (v));
            worse (worstRatio, out / in);
        }
    }
    okp (worstRatio <= 1.0 + 1.0e-6, "os 1: max |out| <= max |in| * (1 + 1e-6) on squares, DC steps and bursts, driveDb 6..36 "
                                     "(worst ratio " + num (worstRatio, "%.9g") + ")");

    std::printf ("       os 4, autoComp 0: overshoot 20*log10(max|out| / max|in|) in dB, Transformer | Tanh\n");
    bool noWorse = true;
    for (const M& m : ms)
    {
        double in = 0.0;
        for (float v : m.x) worse (in, (double) std::fabs (v));
        std::printf ("       %-30s", m.name);
        for (float db : { 6.0f, 12.0f, 24.0f, 36.0f })
        {
            double o[2] {};
            int j = 0;
            for (Shape s : { Shape::Transformer, Shape::Tanh })
            {
                const std::vector<float> y = render (P (s, db, 1.0f, 0.0f), fs, 4, m.x);
                double out = y.empty() ? std::numeric_limits<double>::quiet_NaN() : 0.0;
                for (float v : y) worse (out, (double) std::fabs (v));
                o[j++] = 20.0 * std::log10 (out / in);
            }
            std::printf (" | %4.0f dB: %+.3f %+.3f", (double) db, o[0], o[1]);
            noWorse = noWorse && o[0] <= o[1] + 0.1;
        }
        std::printf ("\n");
    }
    okp (noWorse, "os 4: Transformer's overshoot <= Tanh's + 0.1 dB on every row");
}

void unitSlope()
{
    bool one = true;
    for (int i = 0; i <= 8000; ++i)
    {
        WS w; w.setShape (Shape::Transformer); w.setDrive ((float) (std::pow (10.0, 0.01 * i / 20.0) - 1.0));
        one = one && bitEq (w.slopeAtZero(), 1.0f);
    }
    ok (one, "slopeAtZero() == 1.0f exactly at driveDb 0..80 in 0.01 dB steps");
    // comp = slopeAtZero^(-autoComp) multiplies the wet; exactly 1.0f means autoComp changes no bit of the output.
    // At the drives where (1/k)·k rounds away from 1 in float — where a slope computed as norm·k would not be 1.
    std::vector<float> drives;
    for (int i = 1; i <= 4000 && drives.size() < 4; ++i)
    {
        const float db = 0.01f * (float) i, k = (float) (std::pow (10.0, (double) db / 20.0) - 1.0);
        const float nrm = 1.0f / k;
        if (! bitEq (nrm * k, 1.0f)) drives.push_back (db);
    }
    drives.push_back (6.0f);
    long long bad = drives.size() == 5 ? 0 : -1;
    const std::vector<float> x = lfRich (48000.0, 24000);
    for (float db : drives)
    {
        const std::vector<float> r = render (P (Shape::Transformer, db, 1.0f, 0.0f), 48000.0, 1, x);
        for (float ac : { 0.5f, 1.0f })
            bad += diffs (r, render (P (Shape::Transformer, db, 1.0f, ac), 48000.0, 1, x), 0, r.size());
    }
    std::string at;
    for (float db : drives) at += (at.empty() ? "" : ", ") + num (db);
    okp (bad == 0, "the stage's drive-compensation is exactly 1.0f: autoComp 0.5 and 1 render autoComp 0's bits at driveDb " + at
                   + " (all but 6 are drives where (1/k)*k != 1 in float; " + std::to_string (bad) + " differ)");
}

//==============================================================================
void lifecycle()
{
    const double fs = 48000.0;
    // NaN, then silence: the gate keeps the model finite, and the flush takes L to exact zero.
    for (int os : { 1, 4 })
        for (float db : { 0.0f, 12.0f })
        {
            std::vector<float> x = lfRich (fs, 24000);
            x.insert (x.end(), 1000, std::numeric_limits<float>::quiet_NaN());
            x.insert (x.end(), (std::size_t) fs, 0.0f);
            const std::vector<float> y = render (P (Shape::Transformer, db), fs, os, x);
            bool finite = ! y.empty(), zero = ! y.empty();
            for (float v : y) finite = finite && std::isfinite (v);
            for (std::size_t i = y.size() > 4800 ? y.size() - 4800 : 0; i < y.size(); ++i) zero = zero && y[i] == 0.0f;
            ok (finite && zero, "os " + std::to_string (os) + ", driveDb " + num (db) + ": a NaN burst leaves the output finite, "
                                "and after 1 s of silence it is exact 0.0f");
        }

    // Chunking, with a drive glide in the middle of the stream: n = 1, 7, 64, 4096 render the same bits.
    {
        const std::size_t n = 86016, at = 28672;              // both multiples of 7 · 4096
        Buf x { lfRich (fs, n, 1), lfRich (fs, n, 2) };
        Buf ref;
        long long bad = 0;
        for (int blk : { 4096, 1, 7, 64 })
        {
            Sat s; s.setParams (P (Shape::Transformer, 6.0f));
            if (! prep (s, fs, 2, 4)) return;
            Buf y = x;
            through (s, y, blk, 0, at);
            s.setParams (P (Shape::Transformer, 15.0f, 0.8f));
            through (s, y, blk, at, n);
            if (ref.empty()) ref = y;
            else for (int c = 0; c < 2; ++c) bad += diffs (ref[(std::size_t) c], y[(std::size_t) c], 0, n);
        }
        ok (bad == 0, "calls of 1, 7, 64 and 4096 samples, a drive glide included, render the same bits (" + std::to_string (bad) + " differ)");
    }

    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        const std::size_t n = 14400;
        // reset() is a fresh prepare(), bit for bit — after a loud LF passage that left the flux far from 0.
        {
            Sat a, b;
            a.setParams (P (Shape::Transformer, 12.0f)); b.setParams (P (Shape::Transformer, 12.0f));
            if (! prep (a, fs, 1, os) || ! prep (b, fs, 1, os)) return;
            Buf warm { lfRich (fs, n, 3) }, ya { lfRich (fs, n, 4) }, yb = ya;
            through (a, warm, 512);
            a.reset();
            through (a, ya, 512);
            through (b, yb, 512);
            const long long d = diffs (ya[0], yb[0], 0, n);
            ok (d == 0, "reset() renders what a fresh prepare() renders" + tag + " (" + std::to_string (d) + " differ)");
        }
        // A channel that leaves and returns starts from L = 0: from the return on, it is a fresh stage's channel.
        {
            Sat a, b;
            a.setParams (P (Shape::Transformer, 12.0f, 0.7f)); b.setParams (P (Shape::Transformer, 12.0f, 0.7f));
            if (! prep (a, fs, 2, os) || ! prep (b, fs, 2, os)) return;
            Buf ya { lfRich (fs, 3 * n, 5), lfRich (fs, 3 * n, 6) };
            through (a, ya, 512, 0, n);
            through (a, ya, 512, n, 2 * n, 1);                          // channel 1 is away
            through (a, ya, 512, 2 * n, 3 * n);
            Buf yb { std::vector<float> (ya[0].size()), lfRich (fs, 3 * n, 6) };
            through (b, yb, 512, 2 * n, 3 * n);
            const long long d = diffs (ya[1], yb[1], 2 * n, 3 * n);
            ok (d == 0, "a channel that leaves and returns renders what a fresh stage renders from the return on" + tag
                        + " (" + std::to_string (d) + " differ)");
        }
        // Transformer -> Tanh -> Transformer starts from L = 0: the same as a stage that never ran the Transformer.
        {
            Sat a, b;
            a.setParams (P (Shape::Transformer, 12.0f)); b.setParams (P (Shape::Tanh, 12.0f));
            if (! prep (a, fs, 1, os) || ! prep (b, fs, 1, os)) return;
            Buf ya { lfRich (fs, 3 * n, 7) }, yb = ya;
            through (a, ya, 512, 0, n);
            through (b, yb, 512, 0, n);
            a.setParams (P (Shape::Tanh, 12.0f));
            through (a, ya, 512, n, 2 * n);
            through (b, yb, 512, n, 2 * n);
            a.setParams (P (Shape::Transformer, 12.0f)); b.setParams (P (Shape::Transformer, 12.0f));
            through (a, ya, 512, 2 * n, 3 * n);
            through (b, yb, 512, 2 * n, 3 * n);
            const long long d = diffs (ya[0], yb[0], 2 * n, 3 * n);
            ok (d == 0, "Transformer -> Tanh -> Transformer renders what a stage that never ran the Transformer renders" + tag
                        + " (" + std::to_string (d) + " differ)");
        }
    }

    // The model state's bytes: Storage says them, prepare() asks for exactly that, and they are
    // kModelFloats floats per channel on top of everything else.
    {
        bool okAll = Sat::kModelFloats == 2;   // slot 0 is the flux; the second is Tape's
        const int taps = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase;
        struct G { double fs; int mb, ch, os; };
        for (const G& g : { G { 48000.0, 512, 2, 4 }, G { 44100.0, 256, 1, 1 }, G { 96000.0, 300, 6, 2 } })
        {
            Sat::Storage st;
            okAll = okAll && Sat::storageFor (g.fs, g.mb, g.ch, g.os, taps, st);
            Sat::Storage without = st; without.model = 0;
            Sat s; s.setParams (P (Shape::Transformer, 6.0f));
            const long long before = alloc::bytes.load();
            const bool okP = s.prepare (g.fs, g.mb, g.ch, g.os, taps);
            const long long got = alloc::bytes.load() - before;
            const bool row = okP && got == (long long) st.bytes() && st.model == (std::size_t) (g.ch * Sat::kModelFloats)
                          && st.bytes() - without.bytes() == sizeof (float) * (std::uint64_t) (g.ch * Sat::kModelFloats);
            std::printf ("       fs %.0f, block %d, %d ch, os %d: Storage::bytes() %llu, prepare() asked for %lld\n",
                         g.fs, g.mb, g.ch, g.os, (unsigned long long) st.bytes(), got);
            okAll = okAll && row;
        }
        ok (okAll, "Storage::bytes() == what prepare() allocated, and the model state is kModelFloats (2) floats per channel of it");
    }
}

// A drive glide lands on the settled floats: past its last period and the downsampler's memory, bit-identical to the
// stage that always had the target (the flux L depends on the input alone, so it is the same in both).
void glideLands()
{
    const double fs = 48000.0;
    const std::size_t n = 16000, at = 3000;
    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        const Sat::Params a = P (Shape::Transformer, 3.0f), b = P (Shape::Transformer, 15.0f, 0.7f, 0.8f, -2.0f);
        Sat m, t;
        m.setParams (a); t.setParams (b);
        if (! prep (m, fs, 2, os) || ! prep (t, fs, 2, os)) return;
        const Buf x { lfRich (fs, n, 8), lfRich (fs, n, 9) };
        Buf ym = x, yt = x;
        through (m, ym, 256, 0, at);
        m.setParams (b);
        through (m, ym, 256, at, n);
        through (t, yt, 256);
        const std::size_t start = (at / 64 + 1) * 64;
        const std::size_t landed = start + (std::size_t) m.glideTicks() * 64 + (std::size_t) m.latencySamples() + 1;
        long long after = 0, during = 0;
        for (int c = 0; c < 2; ++c)
        {
            after  += diffs (ym[(std::size_t) c], yt[(std::size_t) c], landed, n);
            during += diffs (ym[(std::size_t) c], yt[(std::size_t) c], start, landed - 64);
        }
        ok (m.glideTicks() > 0 && after == 0 && during > 0, "a drive glide lands bit-identical to the target stage" + tag + " ("
                                                            + std::to_string (after) + " differ after the landing, "
                                                            + std::to_string (during) + " during the glide)");
    }
}
} // namespace

int main()
{
    std::printf ("felitronics::saturation — the Transformer\n");
    group ("reference NULL at os 1 against the model in double");
    referenceNull();
    group ("behaviour at os 4: THD against the double model at 4*fs, and the top octaves clean");
    behaviourAtOs4();
    group ("rate independence");
    rateIndependence();
    group ("small signal: transparent");
    smallSignal();
    group ("mix and latency");
    mixAndLatency();
    group ("the peak law");
    peakLaw();
    group ("unit slope: drive-compensation is exactly 1.0f");
    unitSlope();
    group ("the model state's lifecycle");
    lifecycle();
    group ("the glide lands on the settled floats");
    glideLands();
    return felitronics::test::report();
}
