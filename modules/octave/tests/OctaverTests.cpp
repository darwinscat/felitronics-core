// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// octave::Octaver and its slider law. The tracker's own claims, each MEASURED on the stage's output: the sub's
// fundamental is half (and a quarter) of the note's at 41.2 and 30.9 Hz, on 44.1/48/96 kHz; it locks within two
// cycles, and a quieter note after a louder one locks too; it stays locked through a 3-second, 40 dB decay; it does
// not double-trigger on a sawtooth or a plucked string; a muted note's sub ends within its own periods, and a lone
// click makes none; it is an exact zero between notes and does not sputter on a decay through the gate; the dry is
// untouched sample for sample. Then the stage contract: the law-11 properties, amount 0 is bit-exact dry with
// nothing computed, a slider move glides on a straight ramp, TONE and SUB HP retune without a click, a pause is
// silence (gate, glides and all), reset restarts, law 8, poison, mono and stereo, prepare binds, no allocation.
//
// THE SUB ALONE is read exactly from a stereo call whose right channel is silent: the stage tracks the mean of its
// channels, so {2x, 0} tracks x, and the right channel comes back as 0 + sub — the sub, bit for bit.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter
#include <law11_call_contract.h>

#include <felitronics/octave/OctaveLaw.h>
#include <felitronics/octave/Octaver.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// The law-11 properties, swept by the shared harness: −1 at unity and −2 at half, so the sub joins the buffer once
// the harness's 220 Hz programme locks.
namespace felitronics::test::law11
{
ADAPT (A_Octaver, felitronics::octave::Octaver, kPrepCh, true, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh);
         felitronics::octave::OctaverParams op;
         op.sub1 = 1.0f; op.sub2 = 0.5f;
         s.setParams (op); return okp; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });
} // namespace felitronics::test::law11

namespace
{
using felitronics::test::ok;
using felitronics::test::group;
namespace oc = felitronics::octave;
using Buf = std::vector<float>;

constexpr double kPi  = 3.14159265358979323846;
constexpr double kLowE = 41.2034, kLowB = 30.8677;

//------------------------------------------------------------------------------ fixtures
std::size_t samples (double fs, double sec) { return (std::size_t) std::lround (sec * fs); }

void fadeIn (Buf& x, double fs, double sec)
{
    const std::size_t n = std::min (x.size(), samples (fs, sec));
    for (std::size_t i = 0; i < n; ++i) x[i] *= (float) ((double) i / (double) n);
}

void fadeOut (Buf& x, double fs, double sec)
{
    const std::size_t n = std::min (x.size(), samples (fs, sec));
    for (std::size_t i = 0; i < n; ++i) x[x.size() - 1 - i] *= (float) ((double) i / (double) n);
}

/** A sine of `amp` peak, `decayDb` down by its end, faded in over `fadeSec`. */
Buf sine (double fs, double hz, double amp, double sec, double fadeSec = 0.01, double decayDb = 0.0)
{
    Buf x (samples (fs, sec));
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double t = (double) i / fs;
        x[i] = (float) (amp * std::pow (10.0, -decayDb * t / sec / 20.0) * std::sin (2.0 * kPi * hz * t));
    }
    fadeIn (x, fs, fadeSec);
    return x;
}

/** A naive sawtooth — every harmonic at 1/n, and the reset's full-height jump. */
Buf saw (double fs, double hz, double amp, double sec)
{
    Buf x (samples (fs, sec));
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double ph = hz * (double) i / fs;
        x[i] = (float) (amp * (2.0 * (ph - std::floor (ph)) - 1.0));
    }
    fadeIn (x, fs, 0.002);
    return x;
}

/** A plucked string, modal: plucked at `beta` of its length, heard by a velocity pickup at `gamma`, the upper modes
    damped faster (`bright` per mode) and a little sharp (stiffness). A pickup near the bridge puts the fundamental
    UNDER the second and third harmonics — the real case a bass octaver has to survive. Peak-normalised to `amp`. */
Buf pluck (double fs, double hz, double sec, double beta, double gamma, double bright, double amp)
{
    struct Mode { double a, f, d; };
    std::vector<Mode> modes;
    for (int n = 1; n <= 80; ++n)
    {
        const double fn = hz * n * (1.0 + 1.0e-4 * n * n);
        if (fn > 0.45 * fs) break;
        modes.push_back ({ std::sin (n * kPi * beta) * std::sin (n * kPi * gamma) / n * std::pow (bright, n - 1), fn,
                           0.6 + 0.08 * n * n * hz / 40.0 });
    }
    Buf x (samples (fs, sec));
    double peak = 0.0;
    std::vector<double> v (x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double t = (double) i / fs;
        double s = 0.0;
        for (const auto& m : modes) s += m.a * std::exp (-m.d * t) * std::sin (2.0 * kPi * m.f * t);
        v[i] = s;
        peak = std::max (peak, std::fabs (s));
    }
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = (float) (amp * v[i] / peak);
    return x;
}

/** The pluck's true fundamental (its modes are slightly sharp). */
double pluckF1 (double hz) { return hz * (1.0 + 1.0e-4); }

/** Karplus–Strong with a pluck-shaped (triangle) excitation; its period is N + 1/2 samples. */
Buf karplus (double fs, int N, double sec, double beta, double amp)
{
    std::vector<double> line ((std::size_t) N);
    const int apex = std::max (1, (int) (beta * N));
    for (int i = 0; i < N; ++i)
        line[(std::size_t) i] = i < apex ? (double) i / apex : (double) (N - i) / (N - apex);
    double mean = 0.0;
    for (double v : line) mean += v;
    mean /= N;
    for (double& v : line) v -= mean;
    Buf x (samples (fs, sec));
    int at = 0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const int nx = at + 1 == N ? 0 : at + 1;
        const double y = line[(std::size_t) at];
        line[(std::size_t) at] = 0.997 * 0.5 * (line[(std::size_t) at] + line[(std::size_t) nx]);
        x[i] = (float) (amp * y);
        at = nx;
    }
    return x;
}

/** Uniform white noise of `amp` peak (a 32-bit LCG). */
Buf noise (std::size_t n, double amp, std::uint32_t seed)
{
    Buf x (n);
    for (auto& v : x)
    {
        seed = seed * 1664525u + 1013904223u;
        v = (float) (amp * ((double) (std::int32_t) seed / 2147483648.0));
    }
    return x;
}

void append (Buf& a, const Buf& b) { a.insert (a.end(), b.begin(), b.end()); }
void addInto (Buf& a, const Buf& b) { for (std::size_t i = 0; i < std::min (a.size(), b.size()); ++i) a[i] += b[i]; }

//------------------------------------------------------------------------------ rendering
oc::OctaverParams levels (float sub1, float sub2)
{
    oc::OctaverParams p;
    p.sub1 = sub1; p.sub2 = sub2;
    return p;
}

void runBlocks (oc::Octaver& o, std::vector<float*> planes, std::size_t n, std::size_t from = 0, int block = 64)
{
    for (std::size_t at = from; at < n; at += (std::size_t) block)
    {
        const int len = (int) std::min ((std::size_t) block, n - at);
        std::vector<float*> p;
        for (auto* q : planes) p.push_back (q + at);
        felitronics::test::run (o.process (p.data(), (int) p.size(), len));
    }
}

/** Mono through a fresh stage at `fs`: the input plus the sub. */
Buf mix (double fs, Buf x, oc::OctaverParams p, int block = 64)
{
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, block, 1));
    o.setParams (p);
    runBlocks (o, { x.data() }, x.size(), 0, block);
    return x;
}

/** The sub ALONE, bit for bit: {2x, 0} through a stereo stage tracks exactly x, and the silent channel comes back
    as the sub. */
Buf subOf (double fs, const Buf& x, oc::OctaverParams p, int block = 64)
{
    Buf l (x.size()), r (x.size(), 0.0f);
    for (std::size_t i = 0; i < x.size(); ++i) l[i] = 2.0f * x[i];
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, block, 2));
    o.setParams (p);
    runBlocks (o, { l.data(), r.data() }, x.size(), 0, block);
    return r;
}

//------------------------------------------------------------------------------ measuring
/** The local peak around each sample: block maxima over `w`-sample blocks, then the max of the five blocks around. */
std::vector<double> localPeak (const Buf& y, std::size_t w)
{
    const std::size_t nb = (y.size() + w - 1) / w;
    std::vector<double> bm (nb, 0.0), lp (y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) bm[i / w] = std::max (bm[i / w], (double) std::fabs (y[i]));
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        const std::size_t b = i / w;
        double m = 0.0;
        for (std::size_t k = b >= 2 ? b - 2 : 0; k <= std::min (nb - 1, b + 2); ++k) m = std::max (m, bm[k]);
        lp[i] = m;
    }
    return lp;
}

/** The sub's EDGES: crossings of ±half its local peak with hysteresis (the high pass turns each step of the square
    into a pulse and a smaller lobe of the other sign; a plain zero crossing would count the lobe). `notePeriod` sizes
    the local window. Interpolated times; `risingOnly` keeps the upward ones. */
std::vector<double> edges (const Buf& y, double fs, double fromSec, double toSec, double notePeriod, bool risingOnly)
{
    const auto lp = localPeak (y, std::max<std::size_t> (1, (std::size_t) (notePeriod / 2.0)));
    std::vector<double> c;
    const std::size_t a = std::max<std::size_t> (1, samples (fs, fromSec)), b = std::min (y.size(), samples (fs, toSec));
    int side = 0;   // -1 below -h, +1 above +h
    for (std::size_t i = a; i < b; ++i)
    {
        const double h = 0.5 * lp[i];
        if (h <= 0.0) { side = 0; continue; }
        const double p = y[i - 1], q = y[i];
        if (q > h && side < 0)
        {
            c.push_back ((double) (i - 1) + (h - p) / (q - p));
            side = 1;
        }
        else if (q < -h && side > 0)
        {
            if (! risingOnly) c.push_back ((double) (i - 1) + (-h - p) / (q - p));
            side = -1;
        }
        else if (side == 0)
            side = q > h ? 1 : (q < -h ? -1 : 0);
    }
    return c;
}

/** The fundamental from the rising edges: (count − 1) periods over the time they span. */
double edgeHz (const Buf& y, double fs, double fromSec, double toSec, double notePeriod)
{
    const auto c = edges (y, fs, fromSec, toSec, notePeriod, true);
    if (c.size() < 3) return 0.0;
    return (double) (c.size() - 1) * fs / (c.back() - c.front());
}

/** Every interval between successive edges must be one note period (±tol) — a lost edge makes one twice as long, a
    double trigger two short ones — AND there must be as many edges as the window holds note periods (±2): a silent
    output has no bad interval and is still a failure. Returns the bad intervals; `count` reports the edges. */
int badHalfPeriods (const Buf& y, double fs, double fromSec, double toSec, double periodSamples, double tol, int* count = nullptr)
{
    const auto c = edges (y, fs, fromSec, toSec, periodSamples, false);
    int bad = 0;
    for (std::size_t i = 1; i < c.size(); ++i)
        if (std::fabs ((c[i] - c[i - 1]) / periodSamples - 1.0) > tol) ++bad;
    const double expected = (toSec - fromSec) * fs / periodSamples;
    if (std::fabs ((double) c.size() - expected) > 2.0) bad += 1000;   // missing (or extra) edges: never clean
    if (count) *count = (int) c.size();
    return bad;
}

/** Level in dB of `hz` over [from, to), Hann-windowed (Goertzel in double). */
double levelDb (const Buf& y, double fs, double fromSec, double toSec, double hz)
{
    const std::size_t a = samples (fs, fromSec), b = std::min (y.size(), samples (fs, toSec));
    const double w = 2.0 * kPi * hz / fs;
    double re = 0.0, im = 0.0;
    const double N = (double) (b - a);
    for (std::size_t i = a; i < b; ++i)
    {
        const double win = 0.5 - 0.5 * std::cos (2.0 * kPi * (double) (i - a) / N);
        re += win * y[i] * std::cos (w * (double) i);
        im -= win * y[i] * std::sin (w * (double) i);
    }
    return 20.0 * std::log10 (std::sqrt (re * re + im * im) / N + 1e-300);
}

double peakOf (const Buf& y, double fs, double fromSec, double toSec)
{
    double p = 0.0;
    for (std::size_t i = samples (fs, fromSec); i < std::min (y.size(), samples (fs, toSec)); ++i)
        p = std::max (p, (double) std::fabs (y[i]));
    return p;
}

/** The first sample at or after `from` from which the output never again exceeds `level` (in absolute value). */
std::size_t quietFrom (const Buf& y, std::size_t from, double level)
{
    std::size_t q = from;
    for (std::size_t i = from; i < y.size(); ++i) if (std::fabs (y[i]) > level) q = i + 1;
    return q;
}

std::string fmt (const char* f, double a, double b = 0.0, double c = 0.0)
{
    char buf[256];
    std::snprintf (buf, sizeof buf, f, a, b, c);
    return buf;
}

//------------------------------------------------------------------------------ the tracker's claims
void pitch()
{
    group ("the sub is exactly half (and a quarter) of the note — measured on the output, 44.1 / 48 / 96 kHz");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (double f0 : { kLowE, kLowB })
        {
            const Buf x  = sine (fs, f0, 0.5, 4.0);
            const Buf y1 = subOf (fs, x, levels (1.0f, 0.0f));
            const Buf y2 = subOf (fs, x, levels (0.0f, 1.0f));

            const double h1 = edgeHz (y1, fs, 1.0, 4.0, fs / f0), h2 = edgeHz (y2, fs, 1.0, 4.0, 2.0 * fs / f0);
            const double e1 = levelDb (y1, fs, 1.0, 4.0, f0 / 2.0) - levelDb (y1, fs, 1.0, 4.0, f0);
            const double e2 = levelDb (y2, fs, 1.0, 4.0, f0 / 4.0) - levelDb (y2, fs, 1.0, 4.0, f0 / 2.0);
            std::printf ("    %5.1f kHz, %.4f Hz: -1 at %.5f Hz (x%.6f), -2 at %.5f Hz (x%.6f); "
                         "the even partial under the sub: %.1f / %.1f dB\n",
                         fs / 1000.0, f0, h1, h1 / f0, h2, h2 / f0, e1, e2);
            const std::string at = fmt ("%.0f Hz, %.2f Hz", fs, f0);
            ok (std::fabs (h1 / (f0 / 2.0) - 1.0) < 1.0e-4, at + ": -1 octave is f/2 within 1e-4");
            ok (std::fabs (h2 / (f0 / 4.0) - 1.0) < 1.0e-4, at + ": -2 octaves is f/4 within 1e-4");
            // A divider that skipped or doubled an edge would break the square's symmetry and put energy on the
            // NOTE's own frequency (for -1) — a 50 % square has none, and neither does its high-passed self.
            ok (e1 > 40.0, at + ": -1 carries nothing at f — 40 dB under its fundamental");
            ok (e2 > 40.0, at + ": -2 carries nothing at f/2 — 40 dB under its fundamental");
        }
}

/** A note from silence, swept over 16 onset phases x both polarities: when the fade opens, when the sub first
    reaches half its own steady peak, and whether the lock ever lets go in the first 300 ms after it opens. */
void lockFromSilence()
{
    group ("lock: a note from silence, every onset phase and polarity — opens within two cycles, half level soon after, never lets go");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (const auto nb : { std::pair<double, double> { kLowB, 1.75 }, std::pair<double, double> { kLowE, 1.85 },
                               std::pair<double, double> { 97.9989, 2.25 } })
        {
            const double f0 = nb.first;
            double worstOpen = 0.0, worstHalf = 0.0, worstDrop = 0.0, sumHalf = 0.0;
            int stray = 0;
            for (int c = 0; c < 32; ++c)
            {
                Buf x (samples (fs, 0.2), 0.0f);
                Buf n = sine (fs, f0, c < 16 ? 0.5 : -0.5, 1.0, 0.0);
                for (std::size_t i = 0; i < n.size(); ++i)
                    n[i] = (float) ((c < 16 ? 0.5 : -0.5) * std::sin (2.0 * kPi * f0 * (double) i / fs + (c % 16) * kPi / 8.0));
                fadeIn (n, fs, 0.002);
                append (x, n);

                Buf l (x.size()), r (x.size(), 0.0f);
                for (std::size_t i = 0; i < x.size(); ++i) l[i] = 2.0f * x[i];
                oc::Octaver o;
                felitronics::test::run (o.prepare (fs, 64, 2));
                o.setParams (levels (1.0f, 0.0f));
                std::vector<char> locked (x.size(), 0);
                for (std::size_t at = 0; at < x.size(); at += 16)
                {
                    const int len = (int) std::min<std::size_t> (16, x.size() - at);
                    float* io[2] = { l.data() + at, r.data() + at };
                    felitronics::test::run (o.process (io, 2, len));
                    for (int k = 0; k < len; ++k) locked[at + (std::size_t) k] = o.isLocked() ? 1 : 0;
                }
                const std::size_t onset = samples (fs, 0.2);
                const double steady = peakOf (r, fs, 0.9, 1.2);
                std::size_t open = onset, half = onset;
                while (open < r.size() && ! locked[open]) ++open;
                while (half < r.size() && std::fabs (r[half]) < 0.5 * steady) ++half;
                std::size_t drop = 0;
                for (std::size_t i = open; i < std::min (r.size(), open + samples (fs, 0.3)); ++i) drop += locked[i] ? 0 : 1;
                worstOpen = std::max (worstOpen, (double) (open - onset) * f0 / fs);
                worstHalf = std::max (worstHalf, (double) (half - onset) * f0 / fs);
                sumHalf  += (double) (half - onset) * f0 / fs;
                worstDrop = std::max (worstDrop, (double) drop / fs * 1000.0);
                stray += badHalfPeriods (r, fs, (double) half / fs + 2.0 / f0, 1.2, fs / f0, 0.03) == 0 ? 0 : 1;
            }
            std::printf ("    %5.1f kHz, %5.2f Hz, 32 onsets: the fade opens by %.2f cycles; half level by %.2f (mean %.2f) = %.1f ms; "
                         "unlocked after opening %.1f ms; onsets with a stray half-period after %d\n",
                         fs / 1000.0, f0, worstOpen, worstHalf, sumHalf / 32.0, worstHalf / f0 * 1000.0, worstDrop, stray);
            const std::string w = fmt ("%.0f Hz, %.2f Hz", fs, f0);
            ok (worstOpen <= 2.0, w + ": the fade opens within two cycles at every onset phase and polarity");
            ok (worstHalf <= nb.second, w + fmt (": half level within %.2f cycles at every onset", nb.second));
            ok (worstDrop == 0.0, w + ": once open, the lock never lets go in the first 300 ms");
            ok (stray == 0, w + ": no stray half-period from the second cycle after half level, at any onset");
        }
}

/** A louder note, then — legato — a new one at the same level or 12 / 20 dB under, three pitches, eight onset phases:
    the WHOLE transient after the louder note's sub has faded must stay under the new note's own steady level. */
void relock()
{
    group ("relock: after a louder note the new note's sub reopens at its own level — every pitch, phase and drop");
    const double fs = 48000.0;
    double worst = 0.0;
    std::string where;
    for (double hz : { kLowE, 55.0, 97.9989 })
        for (double dB : { 0.0, -12.0, -20.0 })
            for (int ph = 0; ph < 8; ++ph)
            {
                Buf x = sine (fs, kLowE, 0.5, 0.5, 0.002);
                const std::size_t s0 = x.size();
                Buf n (samples (fs, 0.7));
                for (std::size_t i = 0; i < n.size(); ++i)
                    n[i] = (float) (0.5 * std::pow (10.0, dB / 20.0) * std::sin (2.0 * kPi * hz * (double) i / fs + ph * kPi / 4.0));
                fadeIn (n, fs, 0.0005);
                append (x, n);
                const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
                const double steady = peakOf (y, fs, 0.9, 1.2);
                const double from = (double) s0 / fs + 1.25 / kLowE + oc::Octaver::kLockRampSeconds + 0.001;   // the old sub's fade is over
                const double r = peakOf (y, fs, from, 0.9) / steady;
                if (r > worst) { worst = r; where = fmt ("%.1f Hz %+.0f dB phase %.0f", hz, dB, ph); }
            }
    std::printf ("    the largest sub from the old note's fade on, re the new note's steady: %.2f (%s; 7188cfc: 6.29)\n", worst, where.c_str());
    ok (worst <= 1.6, "the relock transient stays under 1.6 of the new note's own steady level, everywhere");
}

void quieterNote()
{
    group ("a quieter note after a louder one locks too — played over it, or after a 60 ms gap");
    const double fs = 48000.0;
    for (double dB : { -12.0, -20.0 })
        for (double gapMs : { 0.0, 60.0 })
        {
            Buf x = sine (fs, kLowE, 0.5, 0.5, 0.002);
            const std::size_t s0 = x.size() + samples (fs, gapMs / 1000.0);
            x.resize (s0, 0.0f);
            append (x, sine (fs, 55.0, 0.5 * std::pow (10.0, dB / 20.0), 1.0, 0.002));
            const Buf y = subOf (fs, x, levels (1.0f, 0.0f));

            // Locked: the first edge of the new note from which every interval is its own period (±5 %).
            const double T = fs / 55.0;
            const auto c = edges (y, fs, (double) s0 / fs, (double) y.size() / fs, T, false);
            std::size_t from = c.size();
            for (std::size_t k = c.size(); k-- > 1; )
            {
                if (std::fabs ((c[k] - c[k - 1]) / T - 1.0) > 0.05) break;
                from = k - 1;
            }
            const double cycles = from < c.size() ? (c[from] - (double) s0) / T : 99.0;
            const double steady = peakOf (y, fs, (double) s0 / fs + 0.6, (double) y.size() / fs);
            // The louder note's sub may run on until its own next edge is a quarter period late, then fades.
            const double gone = 1.25 / kLowE + 0.015;
            const double early = peakOf (y, fs, (double) s0 / fs + 0.020, (double) s0 / fs + 0.040);
            const double under = peakOf (y, fs, (double) s0 / fs + gone, (double) s0 / fs + gone + 0.020);
            std::printf ("    %+.0f dB, gap %2.0f ms: steady from %.2f cycles of the new note; 20…40 ms in the sub peaks at %.2f of the new "
                         "note's steady, %.0f…%.0f ms in at %.2f\n", dB, gapMs, cycles, early / steady, gone * 1000.0, gone * 1000.0 + 20.0, under / steady);
            const std::string w = fmt ("%+.0f dB, gap %.0f ms", dB, gapMs);
            const double bound = gapMs > 0.0 ? 2.5 : 4.5;   // played over the louder note, the lock has to be lost first
            ok (steady > 0.01 && c.size() > 40, w + ": precondition — the new note has a sub");
            ok (cycles <= bound, w + fmt (": the new note's sub is steady within %.1f of its cycles", bound));
            if (gapMs <= 0.0)   // played over the louder note: its sub must not run on at its own level (3971592 held it
                                // there, 3.6 / 9.0 times the new note's steady sub 20…40 ms in)
                ok (early / steady < 1.0, w + ": 20…40 ms in, the louder note's sub is gone (under the new note's own level)");
        }
}

void decay()
{
    group ("a 3-second decay of 40 dB (-6 to -46 dBFS, the gate at -50): locked all the way, the sub following");
    const double fs = 48000.0;
    for (double f0 : { kLowE, kLowB })
    {
        const Buf x = sine (fs, f0, 0.5, 3.0, 0.002, 40.0);
        const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
        int count = 0;
        const int bad = badHalfPeriods (y, fs, 0.1, 2.98, fs / f0, 0.03, &count);
        const double expected = (2.98 - 0.1) * f0;
        const double drop = 20.0 * std::log10 (peakOf (y, fs, 2.85, 2.98) / peakOf (y, fs, 0.1, 0.23));
        const double want = -40.0 * ((2.85 + 2.98) - (0.1 + 0.23)) / 2.0 / 3.0;
        std::printf ("    %.2f Hz: %d edges (%.1f expected), %d intervals off by more than 3 %%; the sub fell %.1f dB (the note %.1f)\n",
                     f0, count, expected, bad, drop, want);
        const std::string w = fmt ("%.2f Hz", f0);
        ok (bad == 0, w + ": every half-period is one note period through the decay, none lost");
        ok (std::fabs (drop - want) < 3.0, w + ": the sub falls with the note (within 3 dB)");
    }
}

void harmonics()
{
    group ("no double trigger on harmonic-rich notes: a sawtooth, a modal pluck (bridge pickup), Karplus-Strong");
    const double fs = 48000.0;
    for (double f0 : { kLowB, kLowE, 55.0, 73.4162, 97.9989 })
    {
        struct Case { const char* name; Buf x; double f1; };
        const int N = (int) std::lround (fs / f0);
        std::vector<Case> cases;
        cases.push_back ({ "saw",            saw (fs, f0, 0.5, 1.5),                          f0 });
        cases.push_back ({ "pluck bridge",   pluck (fs, f0, 1.5, 0.20, 0.08, 0.97, 0.5),       pluckF1 (f0) });
        cases.push_back ({ "pluck middle",   pluck (fs, f0, 1.5, 0.25, 0.20, 0.95, 0.5),       pluckF1 (f0) });
        cases.push_back ({ "Karplus-Strong", karplus (fs, N, 1.5, 0.2, 0.5),                   fs / (N + 0.5) });
        for (auto& c : cases)
        {
            const Buf y = subOf (fs, c.x, levels (1.0f, 0.0f));
            int count = 0;
            const int bad = badHalfPeriods (y, fs, 0.15, 1.45, fs / c.f1, 0.05, &count);
            const double h2 = levelDb (c.x, fs, 0.0, 0.5, 2.0 * c.f1) - levelDb (c.x, fs, 0.0, 0.5, c.f1);
            if (bad > 0)
                std::printf ("    %s at %.2f Hz: %d bad intervals, %d edges\n", c.name, f0, bad, count);
            ok (bad == 0, fmt ("%.2f Hz ", f0) + c.name + fmt (" (2nd harmonic %+.1f dB re the fundamental): every half-period is one note period, none lost", h2));
        }
    }

    // The limit, measured and PINNED, on the hardest shape there is: the fundamental and the 2nd harmonic ALONE (no
    // odd partial to mark the period), the 2nd in eight phases. Where does the comparator start following the 2nd
    // harmonic instead? (Odd partials only help.)
    const struct { double f0; int pin; } limits[] = { { kLowB, -10 }, { kLowE, -12 }, { 55.0, -16 }, { 97.9989, -20 } };
    for (const auto& L : limits)
    {
        const double f0 = L.f0;
        int lastClean = 1;
        for (int under = 0; under <= 30; ++under)
        {
            bool clean = true;
            for (int ph = 0; ph < 8 && clean; ++ph)
            {
                Buf x (samples (fs, 1.5));
                const double a1 = std::pow (10.0, -under / 20.0);
                for (std::size_t i = 0; i < x.size(); ++i)
                {
                    const double th = 2.0 * kPi * f0 * (double) i / fs;
                    x[i] = (float) (0.3 * (a1 * std::sin (th) + std::sin (2.0 * th + ph * kPi / 4.0)));
                }
                fadeIn (x, fs, 0.002);
                const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
                clean = badHalfPeriods (y, fs, 0.15, 1.45, fs / f0, 0.05) == 0;
            }
            if (! clean) break;
            lastClean = -under;
        }
        std::printf ("    %.2f Hz, the fundamental and the 2nd harmonic alone: tracks in every phase down to %d dB under the 2nd\n",
                     f0, lastClean);
        ok (lastClean <= L.pin, fmt ("%.2f Hz: tracks a fundamental %.0f dB under its 2nd harmonic, in every phase", f0, -L.pin));
    }
}

void mutedTail()
{
    group ("a muted note: the sub lets go within its own periods of its last edge — no DC plateau, no thud (at SUB HP 10 Hz)");
    const double fs = 48000.0;
    for (double f0 : { kLowB, kLowE, 55.0, 73.4162, 97.9989 })
    {
        Buf x = sine (fs, f0, 0.5, 1.0, 0.002);
        const std::size_t m = samples (fs, 0.6);
        for (std::size_t i = m; i < x.size(); ++i) x[i] *= (float) std::exp (-(double) (i - m) / (0.002 * fs));   // a palm mute

        // One sample per call, so the tracker's last edge is known to the sample (edgeCount() steps on it).
        Buf l (x.size()), y (x.size(), 0.0f);
        for (std::size_t i = 0; i < x.size(); ++i) l[i] = 2.0f * x[i];
        // At the menu's lowest SUB HP, which lets a held square through as the DC plateau it is; at the default
        // 25 Hz the high pass alone already turns it into a pulse that decays in about 15 ms.
        oc::OctaverParams p = levels (1.0f, 0.0f);
        p.subHpHz = oc::Octaver::kMinSubHpHz;
        oc::Octaver o;
        felitronics::test::run (o.prepare (fs, 64, 2));
        o.setParams (p);
        std::int64_t edgesBefore = 0;
        std::size_t last = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            float* io[2] = { l.data() + i, y.data() + i };
            felitronics::test::run (o.process (io, 2, 1));
            if (o.edgeCount() != edgesBefore) { edgesBefore = o.edgeCount(); if (i <= m) last = i; }   // later ones are
                                                                                                        // the residual's:
                                                                                                        // one edge, no lock
        }
        const double P = peakOf (y, fs, 0.3, 0.6), T = fs / f0;
        const std::size_t down20 = quietFrom (y, m, 0.1 * P), zero = quietFrom (y, m, 0.0);
        // Dying at 5/4 of a period (capped by the 1/22 s loss), then the 4 ms fade to an exact zero.
        const double dying = std::min (1.25 * T, std::ceil (fs / oc::Octaver::kMinTrackHz)) + 1.0;
        const double levBound  = dying + oc::Octaver::kLockRampSeconds * fs;
        const double zeroBound = dying + (oc::Octaver::kLockRampSeconds + 0.001) * fs;
        std::printf ("    %.2f Hz: last edge %.1f ms before the mute; from it, under -20 dB re steady after %.1f ms (bound %.1f), "
                     "an exact zero after %.1f ms (bound %.1f)\n", f0, ((double) m - (double) last) / fs * 1000.0,
                     ((double) down20 - (double) last) / fs * 1000.0, levBound / fs * 1000.0,
                     ((double) zero - (double) last) / fs * 1000.0, zeroBound / fs * 1000.0);
        const std::string w = fmt ("%.2f Hz", f0);
        ok (P > 0.2 && last > samples (fs, 0.5) && last <= m, w + ": precondition — a steady sub, and its last edge just before the mute");
        ok ((double) down20 - (double) last <= levBound, w + ": under -20 dB re steady by 5/4 of a period + the 4 ms fade after the last edge");
        ok ((double) zero - (double) last <= zeroBound, w + ": an exact zero by 5/4 of a period + the fade (+1 ms)");
    }
}

void loneClick()
{
    group ("a lone click makes no sub: one 80 Hz cycle, sharp or smooth (2 ms ramps), every phase and polarity, over silence or hum");
    const double fs = 48000.0;
    for (double humDb : { -300.0, -80.0, -65.0 })
    {
        double worst = 0.0, worstNote = 1.0;
        for (double pk : { 0.3, 0.05 })
            for (int c = 0; c < 32 + 1; ++c)   // 32 smooth phases/polarities, and the sharp negative-first cycle
            {
                Buf x (samples (fs, 0.4), 0.0f);
                const std::size_t c0 = samples (fs, 0.1), L = samples (fs, 0.0125), R = samples (fs, 0.002);
                for (std::size_t i = 0; i < L; ++i)
                {
                    if (c == 32)
                        x[c0 + i] = (float) (-pk * std::sin (2.0 * kPi * 80.0 * (double) i / fs));
                    else
                    {
                        const double g = std::min (1.0, std::min ((double) i / (double) R, (double) (L - i) / (double) R));
                        x[c0 + i] = (float) ((c < 16 ? pk : -pk) * g * std::sin (2.0 * kPi * 80.0 * (double) i / fs + (c % 16) * kPi / 16.0));
                    }
                }
                append (x, sine (fs, kLowE, 0.5, 0.5, 0.002));   // ...then a real note, which must still lock
                for (std::size_t i = 0; i < x.size(); ++i) x[i] += (float) (std::pow (10.0, humDb / 20.0) * std::sin (2.0 * kPi * 60.0 * (double) i / fs + 0.3));
                const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
                worst = std::max (worst, peakOf (y, fs, 0.0, 0.4));
                worstNote = std::min (worstNote, peakOf (y, fs, 0.6, 0.9));
            }
        std::printf ("    hum %s: the sub's largest sample around 66 clicks %g; the weakest note after one %.3f\n",
                     humDb < -200.0 ? "none   " : fmt ("%.0f dBFS", humDb).c_str(), worst, worstNote);
        ok (worst == 0.0, (humDb < -200.0 ? std::string ("no hum") : fmt ("hum %.0f dBFS", humDb)) + ": no click makes any sub");
        ok (worstNote > 0.2, (humDb < -200.0 ? std::string ("no hum") : fmt ("hum %.0f dBFS", humDb)) + ": ...and the note after each one locks");
    }

    // ...and one SMOOTH cycle at every bass pitch (B0 … G2), 16 phases x both polarities, plain or with 0.5 / 2 ms end
    // ramps: the band filters ring on after the source stops, and their ringing must not finish the proof (0faf723
    // opened on it: 0.30 … 0.37 of sub from one cycle at 30.9 / 41.2 / 55 Hz).
    double worstLow = 0.0;
    std::string whereLow;
    for (double f : { kLowB, kLowE, 55.0, 73.4162, 97.9989 })
        for (double fade : { 0.0, 0.0005, 0.002 })
            for (int c = 0; c < 32; ++c)
            {
                Buf x (samples (fs, 0.4), 0.0f);
                const std::size_t start = samples (fs, 0.1), len = (std::size_t) (fs / f);
                for (std::size_t i = 0; i < len; ++i)
                {
                    const double g = fade > 0.0 ? std::min ({ 1.0, (double) i / (fade * fs), (double) (len - 1 - i) / (fade * fs) }) : 1.0;
                    x[start + i] = (float) ((c < 16 ? 0.3 : -0.3) * g * std::sin (2.0 * kPi * f * (double) i / fs + (c % 16) * kPi / 8.0));
                }
                const double p = peakOf (subOf (fs, x, levels (1.0f, 0.0f)), fs, 0.0, 0.4);
                if (p > worstLow) { worstLow = p; whereLow = fmt ("%.2f Hz, %.1f ms ramps, case %.0f", f, fade * 1000.0, c); }
            }
    std::printf ("    one smooth cycle at B0 … G2, 3 ramps x 32 phases/polarities: the sub's largest sample %g %s\n", worstLow, whereLow.c_str());
    ok (worstLow == 0.0, "one cycle at any bass pitch, phase, polarity or ramp makes no sub");

    // A burst of TWO cycles is a note: it is heard (the rule is two rises and two falls, not a length).
    Buf b2 (samples (fs, 0.3), 0.0f);
    for (std::size_t i = 0; i < samples (fs, 0.025); ++i) b2[samples (fs, 0.1) + i] = (float) (0.3 * std::sin (2.0 * kPi * 80.0 * (double) i / fs));
    ok (peakOf (subOf (fs, b2, levels (1.0f, 0.0f)), fs, 0.1, 0.3) > 0.05, "...while a two-cycle burst is a note, and is heard");
}

/** Hum under the gate (a DI's 60 Hz at -80 / -65 dBFS) supplies nothing: a low note from cold locks on its own
    edges and never lets go, and after a note the hum alone leaves the output bit-exact. */
void humUnderTheGate()
{
    group ("hum under the gate proves nothing: a low note over it locks cleanly; after the note, the hum alone is untouched");
    const double fs = 48000.0;
    for (double humDb : { -80.0, -65.0 })
    {
        const auto hum = [&] (Buf& x) { for (std::size_t i = 0; i < x.size(); ++i) x[i] += (float) (std::pow (10.0, humDb / 20.0) * std::sin (2.0 * kPi * 60.0 * (double) i / fs + 0.3)); };
        double worstDrop = 0.0, worstOpen = 0.0;
        for (double f0 : { kLowB, kLowE })
            for (int ph = 0; ph < 4; ++ph)
            {
                Buf x (samples (fs, 1.1), 0.0f);
                for (std::size_t i = samples (fs, 0.3); i < x.size(); ++i)
                    x[i] = (float) (0.5 * std::sin (2.0 * kPi * f0 * (double) (i - samples (fs, 0.3)) / fs + ph * kPi / 2.0));
                hum (x);
                Buf l (x.size()), r (x.size(), 0.0f);
                for (std::size_t i = 0; i < x.size(); ++i) l[i] = 2.0f * x[i];
                oc::Octaver o;
                felitronics::test::run (o.prepare (fs, 64, 2));
                o.setParams (levels (1.0f, 0.0f));
                std::vector<char> locked (x.size(), 0);
                for (std::size_t at = 0; at < x.size(); ++at)
                {
                    float* io[2] = { l.data() + at, r.data() + at };
                    felitronics::test::run (o.process (io, 2, 1));
                    locked[at] = o.isLocked() ? 1 : 0;
                }
                const std::size_t on = samples (fs, 0.3);
                std::size_t open = on, drop = 0;
                while (open < x.size() && ! locked[open]) ++open;
                for (std::size_t i = open; i < on + samples (fs, 0.3); ++i) drop += locked[i] ? 0 : 1;
                worstOpen = std::max (worstOpen, (double) (open - on) * f0 / fs);
                worstDrop = std::max (worstDrop, (double) drop / fs * 1000.0);
            }
        // A note, then the hum alone.
        Buf t = sine (fs, kLowE, 0.5, 0.5, 0.002);
        fadeOut (t, fs, 0.01);
        append (t, Buf (samples (fs, 1.0), 0.0f));
        hum (t);
        const Buf full = mix (fs, t, levels (1.0f, 0.0f));
        bool untouched = true;
        for (std::size_t i = samples (fs, 0.65); i < t.size(); ++i) untouched = untouched && std::memcmp (&full[i], &t[i], sizeof (float)) == 0;
        std::printf ("    hum %.0f dBFS: low B / E from cold open by %.2f cycles, unlocked after opening %.1f ms (7188cfc: locked on the hum 1.1 ms in, "
                     "then dropped ~10 ms); after the note the hum alone is %s\n", humDb, worstOpen, worstDrop, untouched ? "bit-exact" : "TOUCHED");
        const std::string w = fmt ("hum %.0f dBFS", humDb);
        ok (worstOpen <= 2.0 && worstDrop == 0.0, w + ": low B and E lock on their own edges within two cycles and never let go");
        ok (untouched, w + ": from 150 ms after a note, the hum alone comes out bit for bit");
    }

    // A loud G2 stopping dead over hum at -58 dBFS: the old note's follower decays through the gate's opening level
    // while the hum supplies crossings — the hum must not be acquired on that tail (0faf723 acquired it 153 ms after
    // the note, sub to +250 ms at 0.0024).
    {
        Buf x (samples (fs, 2.0), 0.0f);
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = (float) (std::pow (10.0, -58.0 / 20.0) * std::sin (2.0 * kPi * 60.0 * (double) i / fs + 0.3));
        const std::size_t on = samples (fs, 0.3), off = samples (fs, 1.1);
        for (std::size_t i = on; i < off; ++i)
            x[i] += (float) (0.5 * std::min (1.0, (double) (i - on) / 96.0) * std::sin (2.0 * kPi * 97.9989 * (double) (i - on) / fs));
        const Buf y = subOf (fs, x, oc::paramsForAmount (1.0f));
        std::size_t last = off;
        for (std::size_t i = off; i < y.size(); ++i) if (y[i] != 0.0f) last = i;
        std::printf ("    G2 stopping dead over -58 dBFS hum: the sub's last sample %.1f ms after the note\n", (double) (last - off) / fs * 1000.0);
        ok (peakOf (y, fs, 0.5, 1.0) > 0.2 && (double) (last - off) / fs < 0.1, "a loud G2 stopping over -58 dBFS hum: no sub from 100 ms after it — the hum is not acquired");
    }
}

/** A note much quieter than the one before it — played at once, its phase restarting — and a quiet note after a DC
    pulse: the source's own peak must forget the old signal (13d1d3f held it for ever: 0.004 after a 0.5 low B, or
    0.01 after a +0.5 DC pulse, never acquired — the gate open, the sub exactly zero for ten seconds). */
void quietAfterLoud()
{
    group ("a much quieter note after a loud one (and after a DC pulse) is acquired — the source's peak forgets");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
    {
        for (double q : { 0.004, 0.05, 0.005 })   // codex's case, 20 dB and 40 dB under
        {
            Buf x (samples (fs, 3.0), 0.0f);
            const std::size_t s0 = samples (fs, 0.5);
            for (std::size_t i = 0; i < s0; ++i) x[i] = (float) (0.5 * std::sin (2.0 * kPi * kLowB * (double) i / fs));
            fadeIn (x, fs, 0.002);
            for (std::size_t i = s0; i < x.size(); ++i) x[i] = (float) (q * std::sin (2.0 * kPi * kLowE * (double) (i - s0) / fs));
            const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
            const double steady = peakOf (y, fs, 2.0, 3.0);
            std::size_t back = samples (fs, 0.5 + 1.25 / kLowB + 0.006);   // past the loud note's own fade
            while (back < y.size() && y[back] == 0.0f) ++back;
            const double ms = (double) (back - s0) / fs * 1000.0;
            ok (steady > 0.5 * q && ms < 150.0 && badHalfPeriods (y, fs, 1.0, 3.0, fs / kLowE, 0.05) == 0,
                fmt ("%.0f Hz: a %.3f low E right after a 0.5 low B ", fs, q) + fmt ("has its sub %.1f ms after the switch, and tracks", ms));
        }
        Buf d (samples (fs, 2.0), 0.0f);
        for (std::size_t i = samples (fs, 0.1); i < samples (fs, 0.12); ++i) d[i] = 0.5f;
        for (std::size_t i = samples (fs, 0.3); i < d.size(); ++i) d[i] = (float) (0.01 * std::sin (2.0 * kPi * kLowE * (double) (i - samples (fs, 0.3)) / fs));
        const Buf y = subOf (fs, d, levels (1.0f, 0.0f));
        std::size_t back = samples (fs, 0.3);
        while (back < y.size() && y[back] == 0.0f) ++back;
        ok ((double) (back - samples (fs, 0.3)) / fs < 0.1 && badHalfPeriods (y, fs, 0.6, 2.0, fs / kLowE, 0.05) == 0,
            fmt ("%.0f Hz: a 0.01 low E after a +0.5 DC pulse has its sub ", fs) + fmt ("%.1f ms in, and tracks", (double) (back - samples (fs, 0.3)) / fs * 1000.0));
    }
}

/** The supported range starts at 25 Hz: smooth, harmonic notes (four partials, a band-limited saw) at 25, 27.5 (drop
    A) and 30.9 Hz (low B) acquire in BOTH polarities within two and a half cycles. (Below it the 1/22 s deadline is too
    close to the period: codex's 22.2 Hz saw acquires, but only after 0.2 … 0.9 s, depending on its polarity.) */
void supportedRange()
{
    group ("the supported range starts at 25 Hz: smooth notes at 25 / 27.5 / 30.9 Hz acquire in both polarities");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (double f : { 25.0, 27.5, kLowB })
            for (double sign : { 1.0, -1.0 })
            {
                Buf x (samples (fs, 1.5));
                for (std::size_t i = 0; i < x.size(); ++i)
                {
                    double v = 0.0;
                    for (int k = 1; k <= 4; ++k) v += std::sin (2.0 * kPi * k * f * (double) i / fs) / k;
                    x[i] = (float) (sign * 0.6 / kPi * v);
                }
                fadeIn (x, fs, 0.002);
                const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
                std::size_t first = 0;
                while (first < y.size() && y[first] == 0.0f) ++first;
                const double cycles = (double) first * f / fs;
                ok (cycles <= 2.5 && badHalfPeriods (y, fs, 0.4, 1.5, fs / f, 0.05) == 0,
                    fmt ("%.0f Hz, %.1f Hz, ", fs, f) + fmt ("polarity %+.0f: the sub starts %.2f cycles in and tracks", sign, cycles));
            }
}

/** A sustained note right at the gate, whose follower ripples across the opening level (0.84 … 1.28 times it): it is
    acquired from its fresh level over a period, not at the instant of each crossing (0faf723 never acquired it: two
    seconds of exact zero, every rate, both polarities). */
void nearTheGate()
{
    group ("a sustained low B right at the gate — its follower rippling across the opening level — is acquired, every rate and polarity");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (double sign : { 1.0, -1.0 })
        {
            Buf x (samples (fs, 2.0));
            for (std::size_t i = 0; i < x.size(); ++i)
            {
                const double ph = 2.0 * kPi * kLowB * (double) i / fs;
                x[i] = (float) (sign * 0.0024 * (std::sin (ph) + 0.9 * std::sin (2.0 * ph + 0.8)));
            }
            fadeIn (x, fs, 0.002);
            const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
            std::size_t first = 0;
            while (first < y.size() && y[first] == 0.0f) ++first;
            const double at = (double) first / fs;
            ok (at < 0.15 && badHalfPeriods (y, fs, 0.4, 2.0, fs / kLowB, 0.05) == 0,
                fmt ("%.0f Hz, polarity %+.0f: the sub starts ", fs, sign) + fmt ("%.1f ms in and tracks for the whole note", at * 1000.0));
        }
}

void silenceBetweenNotes()
{
    group ("between notes the sub is an exact zero, and a decay through the gate does not sputter");
    const double fs = 48000.0;

    // A note, a gap with a -70 dBFS noise floor, the next note.
    Buf x = sine (fs, kLowE, 0.5, 1.0, 0.002);
    fadeOut (x, fs, 0.01);
    append (x, Buf (samples (fs, 1.0), 0.0f));
    append (x, sine (fs, 55.0, 0.5, 0.6, 0.002));
    addInto (x, noise (x.size(), std::pow (10.0, -70.0 / 20.0), 0x0C7Au));
    const Buf y = subOf (fs, x, levels (1.0f, 0.0f));
    const Buf full = mix (fs, x, levels (1.0f, 0.0f));

    std::size_t last = samples (fs, 1.0);
    for (std::size_t i = samples (fs, 1.0); i < samples (fs, 2.0); ++i) if (y[i] != 0.0f) last = i;
    const double after = (double) (last - samples (fs, 1.0)) / fs * 1000.0;
    bool dryInGap = true;
    for (std::size_t i = samples (fs, 1.15); i < samples (fs, 2.0); ++i) dryInGap = dryInGap && std::memcmp (&full[i], &x[i], sizeof (float)) == 0;
    const double nextPeak = peakOf (y, fs, 2.0, 2.06);
    std::printf ("    the sub's last non-zero sample: %.1f ms after the note ended; the next note's sub within 3 cycles: %.3f\n", after, nextPeak);
    ok (dryInGap && peakOf (y, fs, 1.15, 2.0) == 0.0, "the gap (from 150 ms after the note, -70 dBFS noise in it): output = input, bit for bit");
    ok (nextPeak > 0.2, "...and the next note's sub is back within 3 of its cycles");

    // ...and a note that never opens the gate (-53 dBFS: the comparator follows it, the gate at -50 stays shut) is
    // neither heard nor counted.
    {
        const Buf n = sine (fs, kLowE, std::pow (10.0, -53.0 / 20.0), 1.0, 0.002);
        oc::Octaver o;
        felitronics::test::run (o.prepare (fs, 64, 1));
        Buf m = n;
        runBlocks (o, { m.data() }, m.size());
        ok (std::memcmp (m.data(), n.data(), m.size() * sizeof (float)) == 0 && o.edgeCount() == 0 && ! o.isLocked(),
            "a -53 dBFS note under a -50 gate: output = input, no edge counted, not locked");
    }

    // A note decaying 60 dB through the gate (-50) over a -80 dBFS floor: the sub stops ONCE.
    Buf d = sine (fs, kLowE, 0.1, 3.0, 0.002, 60.0);
    addInto (d, noise (d.size(), std::pow (10.0, -80.0 / 20.0), 0xD3CAu));
    const Buf z = subOf (fs, d, levels (1.0f, 0.0f));
    int runs = 0;
    std::size_t quiet = 0, end = 0;
    bool inRun = false;
    const std::size_t gap = samples (fs, 0.05);
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (z[i] != 0.0f) { if (! inRun) ++runs; inRun = true; quiet = 0; end = i; }
        else if (inRun && ++quiet >= gap) inRun = false;
    }
    std::printf ("    a decay through the gate: %d stretch(es) of sub, the last sample at %.3f s (the note crosses -50 dBFS at %.2f s)\n",
                 runs, (double) end / fs, 3.0 * 30.0 / 60.0);
    ok (runs == 1, "a decay through the gate stops the sub once — no sputter");
}

void latency()
{
    group ("no added latency: each output sample is its input sample plus the sub, bit for bit");
    const double fs = 48000.0;
    Buf x = pluck (fs, kLowE, 1.0, 0.2, 0.1, 0.96, 0.6);
    append (x, pluck (fs, 55.0, 1.0, 0.15, 0.1, 0.97, 0.7));
    const Buf out = mix (fs, x, levels (1.0f, 0.5f));
    const Buf sub = subOf (fs, x, levels (1.0f, 0.5f));
    std::size_t same = 0;
    for (std::size_t i = 0; i < x.size(); ++i) if (out[i] == x[i] + sub[i]) ++same;
    ok (same == x.size(), fmt ("every output sample is the input sample plus the sub (%.0f of %.0f)", (double) same, (double) x.size()));
    ok (oc::Octaver::latencySamples() == 0, "latencySamples() == 0");
    ok (peakOf (sub, fs, 0.2, 2.0) > 0.1, "precondition: the sub is there to be added");
}

void subHighPass()
{
    group ("SUB HP: the sub's 2-pole high pass is where the Butterworth formula puts it");
    const double fs = 48000.0;
    const Buf x = sine (fs, kLowB, 0.5, 4.0);
    auto at = [&] (float hp) { oc::OctaverParams p = levels (1.0f, 0.0f); p.subHpHz = hp; return levelDb (subOf (fs, x, p), fs, 1.0, 4.0, kLowB / 2.0); };
    const double f = kLowB / 2.0, lo = at (10.0f), mid = at (25.0f), hi = at (60.0f);
    const auto bw = [f] (double fc) { return -10.0 * std::log10 (1.0 + std::pow (fc / f, 4.0)); };
    std::printf ("    the -1 of low B (%.2f Hz) at SUB HP 10 / 25 / 60 Hz: %.2f / %.2f / %.2f dB re 10 Hz (formula %.2f / %.2f)\n",
                 f, 0.0, mid - lo, hi - lo, bw (25.0) - bw (10.0), bw (60.0) - bw (10.0));
    ok (std::fabs ((mid - lo) - (bw (25.0) - bw (10.0))) < 0.5 && std::fabs ((hi - lo) - (bw (60.0) - bw (10.0))) < 0.5,
        "the measured cut at 25 and 60 Hz matches the formula within 0.5 dB");
}

//------------------------------------------------------------------------------ the slider and the contract
void sliderLaw()
{
    group ("the one-slider law: a master over the menu's tops; 0 is an exact zero; current params fed back are stable");
    oc::OctaverParams menu;
    menu.sub1Top = 2.0f; menu.sub2Top = 1.5f; menu.toneHz = 333.0f; menu.subHpHz = 40.0f; menu.gateDb = -61.0f;
    const auto z = oc::paramsForAmount (0.0f, menu);
    ok (felitronics::core::sameBits (z.sub1, 0.0f) && felitronics::core::sameBits (z.sub2, 0.0f),
        "amount 0 -> both sub levels +0.0 exactly, whatever the tops");
    bool mono = true, voice = true;
    float prev = -1.0f;
    for (int k = 0; k <= 1000; ++k)
    {
        const float a = (float) k / 1000.0f;
        const auto p = oc::paramsForAmount (a, menu);
        mono = mono && p.sub1 > prev && p.sub2 == 1.5f * oc::sliderTaper (a) + 0.0f;
        prev = p.sub1;
        voice = voice && p.sub1Top == 2.0f && p.sub2Top == 1.5f && p.toneHz == 333.0f && p.subHpHz == 40.0f && p.gateDb == -61.0f;
    }
    ok (mono, "sub1 rises strictly with the slider; sub2 is its top times the same taper");
    ok (voice, "the tops, tone, sub HP and gate pass through from the menu at every slider position");
    for (float a : { 0.25f, 0.5f, 0.75f, 1.0f })
        std::printf ("    amount %.2f -> sub1 %.4f (%+.1f dB) with the default tops\n", a, oc::paramsForAmount (a).sub1, 20.0 * std::log10 (oc::paramsForAmount (a).sub1));
    ok (oc::paramsForAmount (1.0f).sub1 == oc::OctaverParams {}.sub1Top && oc::paramsForAmount (1.0f).sub2 == 0.0f,
        "amount 1 with the default menu -> -1 at its top (1.25), -2 off");
    ok (oc::paramsForAmount (-1.0f, menu).sub1 == 0.0f && oc::paramsForAmount (2.0f, menu).sub1 == 2.0f, "out of range is clamped");

    // The family's convention: the stage's CURRENT params are the menu. Fed back, the law gives the same levels at
    // the same amount (3971592's tops-in-the-level-fields halved them each call: 0.75 -> 0.28), and comes back from 0.
    const auto once  = oc::paramsForAmount (0.5f, menu);
    const auto twice = oc::paramsForAmount (0.5f, once);
    const auto off   = oc::paramsForAmount (0.0f, twice);
    const auto back  = oc::paramsForAmount (0.5f, off);
    ok (twice.sub1 == once.sub1 && twice.sub2 == once.sub2 && back.sub1 == once.sub1 && back.sub2 == once.sub2,
        fmt ("current params fed back: %.4f -> %.4f -> 0 -> %.4f (stable, and back from zero)", once.sub1, twice.sub1, back.sub1));
    const auto nanP = oc::paramsForAmount (std::numeric_limits<float>::quiet_NaN(), once);
    ok (std::memcmp (&nanP, &once, sizeof (once)) == 0, "a NaN amount returns the menu (the current params) untouched");
}

void amountZero()
{
    group ("amount 0 is bit-exact dry and computes nothing — from the start, and again after a move and back");
    const double fs = 48000.0;
    Buf l = pluck (fs, kLowE, 1.0, 0.2, 0.1, 0.96, 0.6), r = l;
    addInto (l, noise (l.size(), 0.01, 7u));
    addInto (r, noise (r.size(), 0.02, 8u));
    const Buf l0 = l, r0 = r;

    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (oc::paramsForAmount (0.0f));
    bool idleAll = true;
    for (std::size_t at = 0; at < l.size(); at += 64)
    {
        const int n = (int) std::min<std::size_t> (64, l.size() - at);
        float* io[2] = { l.data() + at, r.data() + at };
        felitronics::test::run (o.process (io, 2, n));
        idleAll = idleAll && o.isIdle();
    }
    ok (std::memcmp (l.data(), l0.data(), l.size() * sizeof (float)) == 0
        && std::memcmp (r.data(), r0.data(), r.size() * sizeof (float)) == 0, "a fresh stage at amount 0: both channels bit-exact");
    ok (idleAll && o.edgeCount() == 0, "...idle from the first call: the tracker never ran");

    // Up to 0.7 for half a second, back to 0: once the glide lands and the next grid point passes, exact again.
    Buf m = sine (fs, kLowE, 0.5, 1.5);
    const Buf m0 = m;
    const std::size_t up = samples (fs, 0.2), down = samples (fs, 0.7);
    std::size_t exactFrom = 0;
    bool movedInBetween = false;
    for (std::size_t at = 0; at < m.size(); at += 64)
    {
        if (at == up / 64 * 64)   o.setParams (oc::paramsForAmount (0.7f));
        if (at == down / 64 * 64) o.setParams (oc::paramsForAmount (0.0f));
        float* io[1] = { m.data() + at };
        felitronics::test::run (o.process (io, 1, 64));
        if (at > down && o.isIdle() && exactFrom == 0) exactFrom = at + 64;
    }
    for (std::size_t i = up; i < down; ++i) movedInBetween = movedInBetween || m[i] != m0[i];
    // The wake was mid-note: a restart through the door, and from the third cycle on the sub's edges are the note's.
    Buf woke (m.size());
    for (std::size_t i = 0; i < m.size(); ++i) woke[i] = m[i] - m0[i];
    int wakeEdges = 0;
    const int wakeBad = badHalfPeriods (woke, fs, (double) up / fs + 3.0 / kLowE, (double) down / fs, fs / kLowE, 0.03, &wakeEdges);
    ok (wakeBad == 0, fmt ("woken mid-note: every half-period from the third cycle is the note's period (%.0f edges)", wakeEdges));
    bool exactAfter = exactFrom > 0;
    for (std::size_t i = exactFrom; exactAfter && i < m.size(); ++i) exactAfter = std::memcmp (&m[i], &m0[i], sizeof (float)) == 0;
    std::printf ("    idle again %.1f ms after the slider returned to 0\n", (double) (exactFrom - down / 64 * 64) / fs * 1000.0);
    ok (movedInBetween, "precondition: at 0.7 the sub was there");
    ok (exactAfter && (double) (exactFrom - down / 64 * 64) / fs <= oc::Octaver::kLevelRampSeconds + 2.0 * 64.0 / fs,
        "after the move back to 0: bit-exact from the first grid point after the glide lands");
}

void glide()
{
    group ("a slider move glides: a straight counted ramp, landed in 20 ms, no step larger than the steady sound's own");
    const double fs = 48000.0;
    const Buf x = sine (fs, kLowE, 0.5, 3.0);
    const Buf steady = mix (fs, x, oc::paramsForAmount (1.0f));
    double stepSteady = 0.0, stepDry = 0.0;
    for (std::size_t i = samples (fs, 0.5); i < x.size(); ++i)
    {
        stepSteady = std::max (stepSteady, (double) std::fabs (steady[i] - steady[i - 1]));
        stepDry    = std::max (stepDry, (double) std::fabs (x[i] - x[i - 1]));
    }

    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 1));
    o.setParams (oc::paramsForAmount (0.2f));
    Buf y = x;
    double step = 0.0;
    int offTrack = 0, trackPoints = 0;
    const std::size_t jump = samples (fs, 1.0) / 64 * 64;
    const float from = oc::paramsForAmount (0.2f).sub1, to = oc::paramsForAmount (1.0f).sub1;
    const int len = (int) std::lround (oc::Octaver::kLevelRampSeconds * fs);
    for (std::size_t at = 0; at < y.size(); at += 64)
    {
        if (at == jump) o.setParams (oc::paramsForAmount (1.0f));
        if (at == samples (fs, 2.0) / 64 * 64) o.setParams (oc::paramsForAmount (0.0f));
        float* io[1] = { y.data() + at };
        felitronics::test::run (o.process (io, 1, 64));
        if (at >= jump && at < jump + (std::size_t) len + 128)
        {
            // The ramp's own arithmetic, recomputed: from + (to - from) * (pos / len), pos = samples since the jump.
            const int pos = std::min (len, (int) (at + 64 - jump));
            const float want = pos >= len ? to : from + (to - from) * ((float) pos * (1.0f / (float) len));
            ++trackPoints;
            if (! felitronics::core::sameBits (o.sub1Level(), want)) ++offTrack;
        }
    }
    for (std::size_t i = 1; i < y.size(); ++i) step = std::max (step, (double) std::fabs (y[i] - y[i - 1]));
    std::printf ("    largest step: %.5f through 0.2 -> 1 -> 0 jumps; %.5f steady at 1; %.5f dry alone\n", step, stepSteady, stepDry);
    ok (offTrack == 0 && trackPoints > 15, fmt ("the level follows the straight ramp bit for bit and lands on its target (%.0f block points)", trackPoints));
    ok (step <= 1.02 * std::max (stepSteady, stepDry), "the jumps make no step larger than the steady sound's own (2 %)");
}

void retune()
{
    group ("TONE and SUB HP retune without a click: a per-sample glide over the whole range, and no added curvature, 44.1 / 48 / 96 kHz");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
    {
        // The WHOLE glide, both ways, over the whole range: the corner moves every sample, by no more than 0.5 %,
        // and lands exactly. (7188cfc re-designed every 16 samples: 7.25 % at a time.)
        double worst = 1.0;
        bool landed = true;
        for (int dir = 0; dir < 2; ++dir)
        {
            oc::Octaver o;
            felitronics::test::run (o.prepare (fs, 64, 1));
            oc::OctaverParams p = levels (1.0f, 0.0f);
            p.toneHz  = dir == 0 ? oc::Octaver::kMinToneHz : oc::Octaver::kMaxToneHz;
            p.subHpHz = dir == 0 ? oc::Octaver::kMinSubHpHz : oc::Octaver::kMaxSubHpHz;
            o.setParams (p);
            float z = 0.0f;
            float* io[1] = { &z };
            felitronics::test::run (o.process (io, 1, 1));
            o.setToneHz (dir == 0 ? oc::Octaver::kMaxToneHz : oc::Octaver::kMinToneHz);
            o.setSubHpHz (dir == 0 ? oc::Octaver::kMaxSubHpHz : oc::Octaver::kMinSubHpHz);
            double t0 = o.toneHzNow(), h0 = o.subHpHzNow();
            for (int i = 0; i < (int) (0.03 * fs); ++i)
            {
                z = 0.0f;
                felitronics::test::run (o.process (io, 1, 1));
                const double t1 = o.toneHzNow(), h1 = o.subHpHzNow();
                worst = std::max ({ worst, t1 / t0, t0 / t1, h1 / h0, h0 / h1 });
                t0 = t1; h0 = h1;
            }
            landed = landed && std::fabs (t0 / (dir == 0 ? oc::Octaver::kMaxToneHz : oc::Octaver::kMinToneHz) - 1.0) < 1.0e-5
                            && std::fabs (h0 / (dir == 0 ? oc::Octaver::kMaxSubHpHz : oc::Octaver::kMinSubHpHz) - 1.0) < 1.0e-5;
        }
        std::printf ("    %5.1f kHz: TONE 30 <-> 2000 Hz and SUB HP 10 <-> 80 Hz glides — the largest step between two samples x%.4f\n", fs / 1000.0, worst);
        ok (worst <= 1.005 && landed, fmt ("%.0f Hz: the corners move every sample by at most 0.5 %% over the whole range and land", fs));

        const Buf x = sine (fs, kLowE, 0.5, 1.4);
        struct Move { float fromTone, toTone, fromHp, toHp; const char* name; };
        for (const Move mv : { Move { 120.0f, 2000.0f, 25.0f, 25.0f, "tone up" }, Move { 2000.0f, 120.0f, 25.0f, 25.0f, "tone down" },
                               Move { 120.0f, 120.0f, 25.0f, 80.0f, "sub HP up" } })
        {
            oc::OctaverParams held = levels (1.0f, 0.0f);
            held.toneHz = mv.fromTone; held.subHpHz = mv.fromHp;
            const Buf b = subOf (fs, x, held);
            // The held output's own curvature scale: its largest second difference over a whole sub period.
            double scale = 0.0;
            for (std::size_t i = samples (fs, 1.0); i < samples (fs, 1.0) + (std::size_t) (2.0 * fs / kLowE); ++i)
                scale = std::max (scale, (double) std::fabs (b[i] - 2.0f * b[i - 1] + b[i - 2]));
            std::size_t e = samples (fs, 1.0);
            while (! (b[e - 1] < 0.0f && b[e] >= 0.0f)) ++e;   // a transition of the sub
            double worst = 0.0;
            for (double offMs : { -3.0, -1.0, 0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0 })
            {
                const std::size_t n0 = (std::size_t) ((double) e + offMs * fs / 1000.0);
                Buf l (x.size()), r (x.size(), 0.0f);
                for (std::size_t i = 0; i < x.size(); ++i) l[i] = 2.0f * x[i];
                oc::Octaver o;
                felitronics::test::run (o.prepare (fs, 64, 2));
                o.setParams (held);
                runBlocks (o, { l.data(), r.data() }, n0);
                oc::OctaverParams moved = held;
                moved.toneHz = mv.toTone; moved.subHpHz = mv.toHp;
                o.setParams (moved);
                runBlocks (o, { l.data(), r.data() }, x.size(), n0);
                for (std::size_t i = n0; i < n0 + samples (fs, 0.002); ++i)
                    worst = std::max (worst, std::fabs ((double) (r[i] - 2.0f * r[i - 1] + r[i - 2]) - (double) (b[i] - 2.0f * b[i - 1] + b[i - 2])) / scale);
            }
            std::printf ("    %5.1f kHz %-9s: the curvature it adds in the first 2 ms, worst over 9 phases: %.2f of the held sound's largest\n",
                         fs / 1000.0, mv.name, worst);
            ok (worst <= 1.0, fmt ("%.0f Hz ", fs) + mv.name + ": the retune adds no more curvature than the held sound's own largest");
        }
    }
}

void pauseIsSilence()
{
    group ("law 11c: a pause IS silence — the gate holds, the glides travel, and what comes back is bit-identical");
    const double fs = 48000.0;
    int same = 0, cases = 0, heard = 0;
    for (int gapLen : { 1, 2, 7, 64, 65, 480, 4801, 48000 })
    {
        oc::Octaver a, b;
        felitronics::test::run (a.prepare (fs, 64, 2) && b.prepare (fs, 64, 2));
        a.setParams (oc::paramsForAmount (0.5f));
        b.setParams (oc::paramsForAmount (0.5f));
        Buf t = sine (fs, kLowE, 0.5, 0.5);
        Buf u = t, v = t, w = t;
        runBlocks (a, { t.data(), u.data() }, t.size());
        runBlocks (b, { v.data(), w.data() }, v.size());
        a.setParams (oc::paramsForAmount (1.0f));   // a level glide in flight...
        b.setParams (oc::paramsForAmount (1.0f));
        a.setToneHz (300.0f);                        // ...and a corner glide
        b.setToneHz (300.0f);
        Buf z1 ((std::size_t) gapLen, 0.0f), z2 = z1;
        float* iz[2] = { z1.data(), z2.data() };
        felitronics::test::run (a.process (nullptr, 0, gapLen));   // a pause...
        felitronics::test::run (b.process (iz, 2, gapLen));        // ...and the same length of silence
        // A quiet note on its return (-44 dBFS: over the gate's opening level, so it is acquired again whatever the
        // gate held through the gap — the BITS must still be the same).
        Buf q1 = sine (fs, kLowE, std::pow (10.0, -44.0 / 20.0), 1.0, 0.0), q2 = q1, q3 = q1, q4 = q1;
        runBlocks (a, { q1.data(), q2.data() }, q1.size());
        runBlocks (b, { q3.data(), q4.data() }, q3.size());
        ++cases;
        if (std::memcmp (q1.data(), q3.data(), q1.size() * sizeof (float)) == 0
            && std::memcmp (q2.data(), q4.data(), q2.size() * sizeof (float)) == 0) ++same;
        const Buf qa = sine (fs, kLowE, std::pow (10.0, -44.0 / 20.0), 1.0, 0.0);
        double d = 0.0;
        for (std::size_t i = samples (fs, 0.2); i < q3.size(); ++i) d = std::max (d, (double) std::fabs (q3[i] - qa[i]));
        if (d > 1.0e-3) ++heard;
        if (gapLen == 1 || gapLen == 48000)
            std::printf ("    gap %5d: the quiet note's sub after it peaks at %.2e (through silence) and the two returns are %s\n",
                         gapLen, d, std::memcmp (q1.data(), q3.data(), q1.size() * sizeof (float)) == 0 ? "bit-identical" : "DIFFERENT");
    }
    ok (same == cases, fmt ("after a pause, every sample equals the same stream after silence (%.0f of %.0f gap lengths)", same, cases));
    ok (heard == cases, "precondition: every return has its sub");

    // A long pause parks: ten minutes in one call costs no more than its first second and a half.
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (levels (1.0f, 1.0f));
    Buf t = sine (fs, kLowE, 0.5, 0.5), u = t;
    runBlocks (o, { t.data(), u.data() }, t.size());
    const auto t0 = std::chrono::steady_clock::now();
    felitronics::test::run (o.process (nullptr, 0, (int) (600.0 * fs)));
    const double ms = std::chrono::duration<double, std::milli> (std::chrono::steady_clock::now() - t0).count();
    std::printf ("    a ten-minute pause in one call: %.2f ms\n", ms);
    ok (ms < 200.0, "a ten-minute pause in one call parks on the silent fixed point (well under a second of work)");

    // Idle through a pause stays idle.
    oc::Octaver idle;
    felitronics::test::run (idle.prepare (fs, 64, 1));
    idle.setParams (oc::paramsForAmount (0.0f));
    felitronics::test::run (idle.process (nullptr, 0, 4800));
    ok (idle.isIdle() && idle.edgeCount() == 0, "a pause at amount 0 leaves the stage idle");
}

void resetRestarts()
{
    group ("reset() is a restart: two stages fed different notes answer the next programme with the same bits");
    const double fs = 48000.0;
    oc::Octaver a, b;
    felitronics::test::run (a.prepare (fs, 64, 1) && b.prepare (fs, 64, 1));
    a.setParams (levels (1.0f, 0.7f));
    b.setParams (levels (1.0f, 0.7f));
    Buf p = sine (fs, kLowE, 0.5, 0.4), q = sine (fs, 55.0, 0.3, 0.37);
    runBlocks (a, { p.data() }, p.size());
    runBlocks (b, { q.data() }, q.size());
    b.setParams (oc::paramsForAmount (0.3f));   // glides in flight on one of them
    b.setToneHz (900.0f);
    b.setParams (levels (1.0f, 0.7f));
    b.setToneHz (120.0f);
    a.reset(); a.reset();                       // idempotent
    b.reset();
    Buf x = pluck (fs, kLowB, 1.0, 0.2, 0.1, 0.96, 0.6), y = x;
    runBlocks (a, { x.data() }, x.size());
    runBlocks (b, { y.data() }, y.size());
    ok (std::memcmp (x.data(), y.data(), x.size() * sizeof (float)) == 0, "the same programme after reset: bit-identical");
}

void law8()
{
    group ("law 8: no subnormal state through a note's tail, long silence, or a subnormal-level input");
    const double fs = 48000.0;
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 1));
    o.setParams (levels (1.0f, 1.0f));
    Buf x = sine (fs, kLowE, 0.5, 0.5);
    append (x, Buf (samples (fs, 3.0), 0.0f));
    Buf tiny = sine (fs, kLowE, 1.0e-39, 0.5);   // the input itself subnormal
    append (x, tiny);
    int subn = 0;
    for (std::size_t at = 0; at < x.size(); at += 64)
    {
        float* io[1] = { x.data() + at };
        felitronics::test::run (o.process (io, 1, (int) std::min<std::size_t> (64, x.size() - at)));
        if (o.anySubnormalState()) ++subn;
    }
    ok (subn == 0, fmt ("no subnormal word at any grid point (%.0f blocks)", (double) (x.size() / 64)));
}

void poison()
{
    group ("poison: NaN / Inf / 3e38 in — the dry keeps it, the tracker reads it as nothing, the sub tracks on");
    const double fs = 48000.0;
    Buf l = sine (fs, kLowE, 0.5, 4.0), r = l;
    const float nan = std::numeric_limits<float>::quiet_NaN(), inf = std::numeric_limits<float>::infinity();
    l[48000] = nan; r[60000] = inf; l[72000] = -inf; r[72000] = -inf; l[84000] = 3.0e38f;
    const Buf l0 = l, r0 = r;

    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (levels (1.0f, 0.0f));
    runBlocks (o, { l.data(), r.data() }, l.size());
    int strays = 0, kept = 0;
    for (std::size_t i = 0; i < l.size(); ++i)
    {
        for (int ch = 0; ch < 2; ++ch)
        {
            const float in = ch == 0 ? l0[i] : r0[i], out = ch == 0 ? l[i] : r[i];
            if (std::isfinite (in)) strays += std::isfinite (out) ? 0 : 1;
            else kept += (std::isnan (in) ? std::isnan (out) : out == in) ? 1 : 0;
        }
    }
    ok (strays == 0, "every finite input sample comes out finite");
    ok (kept == 4, "the four non-finite samples come out as they went in (the dry is not ours)");

    const Buf y = subOf (fs, l0, levels (1.0f, 0.0f));
    const double hz = edgeHz (y, fs, 2.5, 4.0, fs / kLowE);
    std::printf ("    after the events: the sub at %.5f Hz; its peak anywhere %.3f\n", hz, peakOf (y, fs, 0.0, 4.0));
    ok (std::fabs (hz / (kLowE / 2.0) - 1.0) < 1.0e-3, "...and the sub is f/2 again after them");
    ok (peakOf (y, fs, 0.0, 4.0) < 2.0 * oc::Octaver::kInputLimit, "...and stayed bounded by the tracker's input limit");
}

void monoStereo()
{
    group ("mono and stereo: the sum is tracked, the same sub goes to both, a dual-mono pair plays as mono");
    const double fs = 48000.0;
    const Buf x = pluck (fs, kLowE, 1.0, 0.2, 0.1, 0.96, 0.6);
    const Buf mono = mix (fs, x, levels (1.0f, 0.5f));

    Buf l = x, r = x;
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (levels (1.0f, 0.5f));
    runBlocks (o, { l.data(), r.data() }, l.size());
    ok (std::memcmp (l.data(), mono.data(), l.size() * sizeof (float)) == 0
        && std::memcmp (r.data(), mono.data(), r.size() * sizeof (float)) == 0, "L = R: each channel is bit-identical to mono");

    // L != R: the sub of the pair is the sub of their sum (read alone from {L + R, 0}), added to each channel.
    Buf a = x, b = noise (x.size(), 0.2, 99u);
    addInto (b, pluck (fs, 55.0, 1.0, 0.2, 0.1, 0.96, 0.4));
    const Buf a0 = a, b0 = b;
    Buf sum (a.size()), zero (a.size(), 0.0f);
    for (std::size_t i = 0; i < a.size(); ++i) sum[i] = a[i] + b[i];
    oc::Octaver s, t;
    felitronics::test::run (s.prepare (fs, 64, 2) && t.prepare (fs, 64, 2));
    s.setParams (levels (1.0f, 0.5f));
    t.setParams (levels (1.0f, 0.5f));
    runBlocks (s, { a.data(), b.data() }, a.size());
    runBlocks (t, { sum.data(), zero.data() }, sum.size());
    std::size_t same = 0;
    for (std::size_t i = 0; i < a.size(); ++i) if (a[i] == a0[i] + zero[i] && b[i] == b0[i] + zero[i]) ++same;
    ok (same == a.size() && peakOf (zero, fs, 0.2, 1.0) > 0.05, "L != R: each channel is its input plus the one sub of the pair, bit for bit");

    // A narrow call on a stereo stage is the mono stage.
    Buf n1 = x;
    oc::Octaver w;
    felitronics::test::run (w.prepare (fs, 64, 2));
    w.setParams (levels (1.0f, 0.5f));
    runBlocks (w, { n1.data() }, n1.size());
    ok (std::memcmp (n1.data(), mono.data(), n1.size() * sizeof (float)) == 0, "a mono call on a stereo-prepared stage: bit-identical to mono");
}

void prepareBinds()
{
    group ("prepare binds: rates outside 8 kHz … 768 kHz, a block under 1, widths outside 1 … 2 are refused and disarm");
    int refused = 0, cases = 0;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    struct P { double fs; int block, ch; };
    for (const P p : { P { 7999.0, 64, 2 }, P { 768001.0, 64, 2 }, P { nan, 64, 2 }, P { 48000.0, 0, 2 },
                       P { 48000.0, 64, 0 }, P { 48000.0, 64, 3 }, P { 48000.0, 64, -1 } })
    {
        oc::Octaver o;
        felitronics::test::run (o.prepare (48000.0, 64, 2));
        ++cases;
        if (! o.prepare (p.fs, p.block, p.ch))
        {
            Buf z (64, 0.25f);
            float* io[1] = { z.data() };
            if (! o.process (io, 1, 64) && ! o.isPrepared()) ++refused;
        }
    }
    ok (refused == cases, fmt ("all %.0f refused, and the refused stage refuses process()", cases));
    for (double fs : { 8000.0, 44100.0, 768000.0 })
    {
        oc::Octaver o;
        ok (o.prepare (fs, 1, 1), fmt ("%.0f Hz, block 1, mono: accepted", fs));
    }
}

void noAllocation()
{
    group ("RT: process, reset and the setters never allocate");
    oc::Octaver o;
    felitronics::test::run (o.prepare (48000.0, 512, 2));
    Buf l = sine (48000.0, kLowE, 0.5, 0.2), r = l;
    const long long before = alloc::count.load();
    o.setParams (oc::paramsForAmount (0.8f));
    o.setToneHz (200.0f);
    o.setSubHpHz (40.0f);
    o.setGateDb (-60.0f);
    for (std::size_t at = 0; at + 512 <= l.size(); at += 512)
    {
        float* io[2] = { l.data() + at, r.data() + at };
        felitronics::test::run (o.process (io, 2, 512));
        felitronics::test::run (o.process (io, 1, 100));
        felitronics::test::run (o.process (nullptr, 0, 77));
    }
    o.reset();
    felitronics::test::okNoAlloc (alloc::count.load() == before, "no allocation: setters, process (stereo, mono, pause), reset");
}
} // namespace

int main()
{
    std::printf ("felitronics octave — Octaver\n");
    felitronics::test::law11::allProperties<felitronics::test::law11::A_Octaver>();
    sliderLaw();
    pitch();
    lockFromSilence();
    relock();
    quieterNote();
    decay();
    harmonics();
    mutedTail();
    loneClick();
    humUnderTheGate();
    nearTheGate();
    quietAfterLoud();
    supportedRange();
    silenceBetweenNotes();
    latency();
    subHighPass();
    amountZero();
    glide();
    retune();
    pauseIsSilence();
    resetRestarts();
    law8();
    poison();
    monoStereo();
    prepareBinds();
    noAllocation();
    return felitronics::test::report();
}
