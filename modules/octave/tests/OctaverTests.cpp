// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// octave::Octaver and its slider law. The tracker's own claims, each MEASURED on the stage's output: the sub's
// fundamental is half (and a quarter) of the note's at 41.2 and 30.9 Hz, on 44.1/48/96 kHz; it locks within two
// cycles; it stays locked through a 3-second, 40 dB decay; it does not double-trigger on a sawtooth or a plucked
// string; it is an exact zero between notes and does not sputter on a decay through the gate; the dry is untouched
// sample for sample. Then the stage contract: the law-11 properties, amount 0 is bit-exact dry with nothing
// computed, a slider move glides, a pause is silence for the glides and drops the tracker, reset restarts, law 8,
// poison, mono and stereo, prepare binds, no allocation.
#include <felitronics_test.h>
#include <alloc_counter.h>   // installs the allocation counter
#include <law11_call_contract.h>

#include <felitronics/octave/OctaveLaw.h>
#include <felitronics/octave/Octaver.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// The law-11 properties, swept by the shared harness: −1 at unity, −2 at half and the dry at half, so the buffer
// moves from the first sample and the sub joins once the harness's 220 Hz programme locks.
namespace felitronics::test::law11
{
ADAPT (A_Octaver, felitronics::octave::Octaver, kPrepCh, true, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh);
         felitronics::octave::OctaverParams op;
         op.sub1 = 1.0f; op.sub2 = 0.5f; op.dry = 0.5f;
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
oc::OctaverParams levels (float sub1, float sub2, float dry)
{
    oc::OctaverParams p;
    p.sub1 = sub1; p.sub2 = sub2; p.dry = dry;
    return p;
}

void renderInto (oc::Octaver& o, Buf& x, int block = 64)
{
    for (std::size_t at = 0; at < x.size(); at += (std::size_t) block)
    {
        const int n = (int) std::min ((std::size_t) block, x.size() - at);
        float* io[1] = { x.data() + at };
        felitronics::test::run (o.process (io, 1, n));
    }
}

/** Mono through a fresh stage at `fs` with the given levels. */
Buf render (double fs, Buf x, oc::OctaverParams p, int block = 64)
{
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, block, 1));
    o.setParams (p);
    renderInto (o, x, block);
    return x;
}

//------------------------------------------------------------------------------ measuring
/** Strict sign changes in [from, to), interpolated; `rising` keeps only the upward ones. */
std::vector<double> crossings (const Buf& y, double fs, double fromSec, double toSec, bool risingOnly)
{
    std::vector<double> c;
    const std::size_t a = std::max<std::size_t> (1, samples (fs, fromSec)), b = std::min (y.size(), samples (fs, toSec));
    for (std::size_t i = a; i < b; ++i)
    {
        const double p = y[i - 1], q = y[i];
        const bool up = p < 0.0 && q > 0.0, down = p > 0.0 && q < 0.0;
        if (up || (down && ! risingOnly))
            c.push_back ((double) (i - 1) + p / (p - q));
    }
    return c;
}

/** The fundamental from rising zero crossings: (count − 1) periods over the time they span. */
double crossingHz (const Buf& y, double fs, double fromSec, double toSec)
{
    const auto c = crossings (y, fs, fromSec, toSec, true);
    if (c.size() < 3) return 0.0;
    return (double) (c.size() - 1) * fs / (c.back() - c.front());
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

/** Every interval between successive sign changes must be one input period (±tol): a missed edge makes one
    interval twice as long, a double trigger makes two short ones. Returns the count of bad intervals. */
int badHalfPeriods (const Buf& y, double fs, double fromSec, double toSec, double periodSamples, double tol, int* count = nullptr)
{
    const auto c = crossings (y, fs, fromSec, toSec, false);
    int bad = 0;
    for (std::size_t i = 1; i < c.size(); ++i)
        if (std::fabs ((c[i] - c[i - 1]) / periodSamples - 1.0) > tol) ++bad;
    if (count) *count = (int) c.size();
    return bad;
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
            const Buf y1 = render (fs, x, levels (1.0f, 0.0f, 0.0f));
            const Buf y2 = render (fs, x, levels (0.0f, 1.0f, 0.0f));

            const double h1 = crossingHz (y1, fs, 1.0, 4.0), h2 = crossingHz (y2, fs, 1.0, 4.0);
            const double e1 = levelDb (y1, fs, 1.0, 4.0, f0 / 2.0) - levelDb (y1, fs, 1.0, 4.0, f0);
            const double e2 = levelDb (y2, fs, 1.0, 4.0, f0 / 4.0) - levelDb (y2, fs, 1.0, 4.0, f0 / 2.0);
            std::printf ("    %5.1f kHz, %.4f Hz: -1 at %.5f Hz (x%.6f), -2 at %.5f Hz (x%.6f); "
                         "the even partial under the sub: %.1f / %.1f dB\n",
                         fs / 1000.0, f0, h1, h1 / f0, h2, h2 / f0, e1, e2);
            const std::string at = fmt ("%.0f Hz, %.2f Hz", fs, f0);
            ok (std::fabs (h1 / (f0 / 2.0) - 1.0) < 1.0e-4, at + ": -1 octave is f/2 within 1e-4");
            ok (std::fabs (h2 / (f0 / 4.0) - 1.0) < 1.0e-4, at + ": -2 octaves is f/4 within 1e-4");
            // A divider that skipped or doubled an edge would break the square's symmetry and put energy on the
            // NOTE's own frequency (for -1) — a 50 % square has none.
            ok (e1 > 40.0, at + ": -1 carries nothing at f — 40 dB under its fundamental");
            ok (e2 > 40.0, at + ": -2 carries nothing at f/2 — 40 dB under its fundamental");
        }
}

void lock()
{
    group ("lock: a note from silence locks within one cycle; no stray edge from then on, and steady from the third cycle");
    for (double fs : { 44100.0, 48000.0, 96000.0 })
        for (double f0 : { kLowE, kLowB })
        {
            Buf x (samples (fs, 0.2), 0.0f);
            append (x, sine (fs, f0, 0.5, 2.0, 0.002));
            const Buf y = render (fs, x, levels (1.0f, 0.0f, 0.0f));

            const double steady = peakOf (y, fs, 1.5, 2.2);
            const std::size_t onset = samples (fs, 0.2);
            std::size_t at = onset;
            while (at < y.size() && std::fabs (y[at]) < 0.5 * steady) ++at;
            const double cycles = (double) (at - onset) * f0 / fs;
            // A stray edge makes a half-period of about 1/2, a lost one about 2: ±10 % separates both from the
            // band filters still settling over the first cycles (the second half-period runs 5.6 % long).
            const int stray  = badHalfPeriods (y, fs, (double) at / fs, 2.2, fs / f0, 0.10);
            const int wobble = badHalfPeriods (y, fs, (double) at / fs + 3.0 / f0, 2.2, fs / f0, 0.01);
            std::printf ("    %5.1f kHz, %.2f Hz: half the steady level %.2f cycles after the onset; stray edges %d, "
                         "half-periods off by > 1 %% from the third cycle %d\n", fs / 1000.0, f0, cycles, stray, wobble);
            const std::string w = fmt ("%.0f Hz, %.2f Hz", fs, f0);
            ok (cycles <= 1.0, w + ": locks within one cycle");
            ok (stray == 0, w + ": no stray or lost edge from the lock on");
            ok (wobble == 0, w + ": every half-period within 1 % of the note's from the third cycle");
        }
}

void decay()
{
    group ("a 3-second decay of 40 dB (-6 to -46 dBFS, the gate at -50): locked all the way, the sub following");
    const double fs = 48000.0;
    for (double f0 : { kLowE, kLowB })
    {
        const Buf x = sine (fs, f0, 0.5, 3.0, 0.002, 40.0);
        const Buf y = render (fs, x, levels (1.0f, 0.0f, 0.0f));
        int count = 0;
        const int bad = badHalfPeriods (y, fs, 0.1, 2.98, fs / f0, 0.03, &count);
        const double expected = (2.98 - 0.1) * f0;
        const double drop = 20.0 * std::log10 (peakOf (y, fs, 2.85, 2.98) / peakOf (y, fs, 0.1, 0.23));
        const double want = -40.0 * ((2.85 + 2.98) - (0.1 + 0.23)) / 2.0 / 3.0;
        std::printf ("    %.2f Hz: %d half-periods (%.1f expected), %d off by more than 3 %%; the sub fell %.1f dB (the note %.1f)\n",
                     f0, count, expected, bad, drop, want);
        const std::string w = fmt ("%.2f Hz", f0);
        ok (bad == 0, w + ": every half-period is one note period through the decay");
        ok (std::fabs ((double) count - expected) <= 2.0, w + ": no edge lost or added through the decay");
        ok (std::fabs (drop - want) < 3.0, w + ": the sub falls with the note (within 3 dB)");
    }
}

void harmonics()
{
    group ("no double trigger on harmonic-rich notes: a sawtooth, a modal pluck (bridge pickup), Karplus-Strong");
    const double fs = 48000.0;
    int worst = 0;
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
            const Buf y = render (fs, c.x, levels (1.0f, 0.0f, 0.0f));
            int count = 0;
            const int bad = badHalfPeriods (y, fs, 0.15, 1.45, fs / c.f1, 0.05, &count);
            const double expected = (1.45 - 0.15) * c.f1;   // one sign change per note period
            const double h2 = levelDb (c.x, fs, 0.0, 0.5, 2.0 * c.f1) - levelDb (c.x, fs, 0.0, 0.5, c.f1);
            worst = std::max (worst, bad);
            if (bad > 0 || std::fabs ((double) count - expected) > 2.0)
                std::printf ("    %s at %.2f Hz: %d bad half-periods, %d of %.0f\n", c.name, f0, bad, count, expected);
            ok (bad == 0, fmt ("%.2f Hz ", f0) + c.name + fmt (" (2nd harmonic %+.1f dB re the fundamental): every half-period is one note period", h2));
            ok (std::fabs ((double) count - expected) <= 2.0, fmt ("%.2f Hz ", f0) + c.name + ": no edge lost");
        }
    }

    // The limit, measured rather than claimed, on the hardest shape there is: the fundamental and the 2nd harmonic
    // ALONE (no odd partial to mark the period), the 2nd in eight phases. Where does the comparator start following
    // the 2nd harmonic instead? Odd partials only help — with 3-6 present it tracked to the end of the sweep.
    for (double f0 : { kLowB, kLowE, 55.0, 97.9989 })
    {
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
                const Buf y = render (fs, x, levels (1.0f, 0.0f, 0.0f));
                clean = badHalfPeriods (y, fs, 0.15, 1.45, fs / f0, 0.05) == 0;
            }
            if (! clean) break;
            lastClean = -under;
        }
        std::printf ("    %.2f Hz, the fundamental and the 2nd harmonic alone: tracks in every phase down to %d dB under the 2nd\n",
                     f0, lastClean);
        ok (lastClean <= -6, fmt ("%.2f Hz: tracks a fundamental 6 dB under its 2nd harmonic, in every phase", f0));
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
    const Buf y = render (fs, x, levels (1.0f, 0.0f, 0.0f));

    std::size_t last = samples (fs, 1.0);
    for (std::size_t i = samples (fs, 1.0); i < samples (fs, 2.0); ++i) if (y[i] != 0.0f) last = i;
    const double after = (double) (last - samples (fs, 1.0)) / fs * 1000.0;
    double gapPeak = 0.0;
    for (std::size_t i = samples (fs, 1.15); i < samples (fs, 2.0); ++i) gapPeak = std::max (gapPeak, (double) std::fabs (y[i]));
    const double nextPeak = peakOf (y, fs, 2.0, 2.06);
    std::printf ("    the sub's last non-zero sample: %.1f ms after the note ended; peak in the gap from 150 ms: %g; "
                 "the next note's sub within 3 cycles: %.3f\n", after, gapPeak, nextPeak);
    ok (gapPeak == 0.0, "the gap (from 150 ms after the note, -70 dBFS noise in it) is an exact zero");
    ok (nextPeak > 0.2, "...and the next note's sub is back within 3 of its cycles");

    // A note decaying 60 dB through the gate (-50) over a -80 dBFS floor: the sub stops ONCE.
    Buf d = sine (fs, kLowE, 0.1, 3.0, 0.002, 60.0);
    addInto (d, noise (d.size(), std::pow (10.0, -80.0 / 20.0), 0xD3CAu));
    const Buf z = render (fs, d, levels (1.0f, 0.0f, 0.0f));
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
    group ("no added latency: the dry is untouched sample for sample, the sub is added on the same sample");
    const double fs = 48000.0;
    Buf x = pluck (fs, kLowE, 1.0, 0.2, 0.1, 0.96, 0.6);
    append (x, pluck (fs, 55.0, 1.0, 0.15, 0.1, 0.97, 0.7));
    const Buf wet = render (fs, x, levels (1.0f, 0.5f, 1.0f));
    const Buf sub = render (fs, x, levels (1.0f, 0.5f, 0.0f));
    std::size_t same = 0;
    for (std::size_t i = 0; i < x.size(); ++i) if (wet[i] == x[i] + sub[i]) ++same;
    ok (same == x.size(), fmt ("every output sample is the input sample plus the sub (%.0f of %.0f)", (double) same, (double) x.size()));
    ok (oc::Octaver::latencySamples() == 0, "latencySamples() == 0");
    ok (peakOf (sub, fs, 0.2, 2.0) > 0.1, "precondition: the sub is there to be added");
}

//------------------------------------------------------------------------------ the slider and the contract
void sliderLaw()
{
    group ("the one-slider law: -1 rises under a unity dry, -2 stays out, 0 is an exact zero");
    ok (felitronics::core::sameBits (oc::paramsForAmount (0.0f).sub1, 0.0f), "amount 0 -> sub1 is +0.0 exactly");
    bool mono = true, rest = true;
    float prev = -1.0f;
    for (int k = 0; k <= 1000; ++k)
    {
        const auto p = oc::paramsForAmount ((float) k / 1000.0f);
        mono = mono && p.sub1 > prev;
        prev = p.sub1;
        rest = rest && p.sub2 == 0.0f && p.dry == 1.0f;
    }
    ok (mono, "sub1 strictly rises with the slider");
    ok (rest, "sub2 is 0 and dry is 1 at every slider position");
    for (float a : { 0.25f, 0.5f, 0.75f, 1.0f })
        std::printf ("    amount %.2f -> sub1 %.4f (%+.1f dB)\n", a, oc::paramsForAmount (a).sub1, 20.0 * std::log10 (oc::paramsForAmount (a).sub1));
    ok (oc::paramsForAmount (1.0f).sub1 == oc::kSliderTopSub1, "amount 1 -> the top");
    ok (oc::paramsForAmount (std::numeric_limits<float>::quiet_NaN()).sub1 == 0.0f
        && oc::paramsForAmount (-1.0f).sub1 == 0.0f && oc::paramsForAmount (2.0f).sub1 == oc::kSliderTopSub1,
        "NaN reads as 0, out of range is clamped");
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
    bool restAll = true;
    for (std::size_t at = 0; at < l.size(); at += 64)
    {
        const int n = (int) std::min<std::size_t> (64, l.size() - at);
        float* io[2] = { l.data() + at, r.data() + at };
        felitronics::test::run (o.process (io, 2, n));
        restAll = restAll && o.isAtRest();
    }
    ok (std::memcmp (l.data(), l0.data(), l.size() * sizeof (float)) == 0
        && std::memcmp (r.data(), r0.data(), r.size() * sizeof (float)) == 0, "a fresh stage at amount 0: both channels bit-exact");
    ok (restAll && o.edgeCount() == 0, "...at rest from the first call: the tracker never ran");

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
        if (at > down && o.isAtRest() && exactFrom == 0) exactFrom = at + 64;
    }
    for (std::size_t i = up; i < down; ++i) movedInBetween = movedInBetween || m[i] != m0[i];
    bool exactAfter = exactFrom > 0;
    for (std::size_t i = exactFrom; exactAfter && i < m.size(); ++i) exactAfter = std::memcmp (&m[i], &m0[i], sizeof (float)) == 0;
    std::printf ("    back at rest %.1f ms after the slider returned to 0\n", (double) (exactFrom - down / 64 * 64) / fs * 1000.0);
    ok (movedInBetween, "precondition: at 0.7 the sub was there");
    ok (exactAfter && (double) (exactFrom - down / 64 * 64) / fs <= oc::Octaver::kLevelRampSeconds + 2.0 * 64.0 / fs,
        "after the move back to 0: bit-exact from the first grid point after the glide lands");
}

void glide()
{
    group ("a slider move glides: no step larger than the steady sound's own");
    const double fs = 48000.0;
    const Buf x = sine (fs, kLowE, 0.5, 3.0);
    const Buf steady = render (fs, x, oc::paramsForAmount (1.0f));
    double stepSteady = 0.0, stepDry = 0.0;
    for (std::size_t i = samples (fs, 0.5); i < x.size(); ++i)
    {
        stepSteady = std::max (stepSteady, (double) std::fabs (steady[i] - steady[i - 1]));
        stepDry    = std::max (stepDry, (double) std::fabs (x[i] - x[i - 1]));
    }

    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 1));
    o.setParams (oc::paramsForAmount (0.0f));
    Buf y = x;
    double step = 0.0;
    bool landed = false;
    for (std::size_t at = 0; at < y.size(); at += 64)
    {
        if (at == samples (fs, 1.0) / 64 * 64) o.setParams (oc::paramsForAmount (1.0f));
        if (at == samples (fs, 2.0) / 64 * 64) o.setParams (oc::paramsForAmount (0.2f));
        if (at == (samples (fs, 1.0) + samples (fs, 0.03)) / 64 * 64) landed = o.sub1Level() == oc::kSliderTopSub1;
        float* io[1] = { y.data() + at };
        felitronics::test::run (o.process (io, 1, 64));
    }
    for (std::size_t i = 1; i < y.size(); ++i) step = std::max (step, (double) std::fabs (y[i] - y[i - 1]));
    std::printf ("    largest step: %.5f through 0 -> 1 -> 0.2 jumps; %.5f steady at 1; %.5f dry alone\n", step, stepSteady, stepDry);
    ok (step <= 1.02 * std::max (stepSteady, stepDry), "the jumps make no step larger than the steady sound's own (2 %)");
    ok (landed, "...and the level has landed on its target 30 ms after the move");
}

void pauseIsSilence()
{
    group ("law 11c: a pause moves the level glides exactly as silence would, and drops the tracker (law 11a)");
    const double fs = 48000.0;
    int same = 0, cases = 0;
    for (int gapLen : { 1, 2, 7, 15, 63, 64, 65, 127, 480, 1000, 4801 })
    {
        oc::Octaver a, b;
        felitronics::test::run (a.prepare (fs, 64, 2) && b.prepare (fs, 64, 2));
        a.setParams (oc::paramsForAmount (0.5f));
        b.setParams (oc::paramsForAmount (0.5f));
        Buf t = sine (fs, kLowE, 0.5, 0.3);
        Buf u = t, v = t, w = t;
        float* ia[2] = { t.data(), u.data() };
        float* ib[2] = { v.data(), w.data() };
        felitronics::test::run (a.process (ia, 2, (int) t.size()) && b.process (ib, 2, (int) t.size()));
        a.setParams (oc::paramsForAmount (1.0f));   // a glide in flight...
        b.setParams (oc::paramsForAmount (1.0f));
        a.setDry (0.3f);
        b.setDry (0.3f);
        Buf z1 ((std::size_t) gapLen, 0.0f), z2 = z1;
        float* iz[2] = { z1.data(), z2.data() };
        felitronics::test::run (a.process (nullptr, 0, gapLen));   // ...through a pause
        felitronics::test::run (b.process (iz, 2, gapLen));        // ...and through silence
        ++cases;
        if (felitronics::core::sameBits (a.sub1Level(), b.sub1Level()) && felitronics::core::sameBits (a.dryLevel(), b.dryLevel())) ++same;
    }
    ok (same == cases, fmt ("the glides stand on the same bits after a pause as after the same silence (%.0f of %.0f gap lengths)", same, cases));

    // ...and the sub cannot ring into the return: silence after a pause is an exact zero at once.
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (levels (1.0f, 1.0f, 1.0f));
    Buf t = sine (fs, kLowE, 0.5, 0.5), u = t;
    float* io[2] = { t.data(), u.data() };
    felitronics::test::run (o.process (io, 2, (int) t.size()));
    const bool wasLocked = o.isLocked();
    felitronics::test::run (o.process (nullptr, 0, 1));
    Buf z1 (64, 0.0f), z2 = z1;
    float* iz[2] = { z1.data(), z2.data() };
    felitronics::test::run (o.process (iz, 2, 64));
    ok (wasLocked && peakOf (z1, fs, 0.0, 1.0) == 0.0 && peakOf (z2, fs, 0.0, 1.0) == 0.0,
        "a locked sub, a one-sample pause, then silence: an exact zero from the first sample");
}

void resetRestarts()
{
    group ("reset() is a restart: two stages fed different notes answer the next programme with the same bits");
    const double fs = 48000.0;
    oc::Octaver a, b;
    felitronics::test::run (a.prepare (fs, 64, 1) && b.prepare (fs, 64, 1));
    a.setParams (levels (1.0f, 0.7f, 1.0f));
    b.setParams (levels (1.0f, 0.7f, 1.0f));
    Buf p = sine (fs, kLowE, 0.5, 0.4), q = sine (fs, 55.0, 0.3, 0.37);
    renderInto (a, p);
    renderInto (b, q);
    b.setParams (oc::paramsForAmount (0.3f));   // a glide in flight on one of them
    b.setParams (levels (1.0f, 0.7f, 1.0f));
    a.reset(); a.reset();                       // idempotent
    b.reset();
    Buf x = pluck (fs, kLowB, 1.0, 0.2, 0.1, 0.96, 0.6), y = x;
    renderInto (a, x);
    renderInto (b, y);
    ok (std::memcmp (x.data(), y.data(), x.size() * sizeof (float)) == 0, "the same programme after reset: bit-identical");
}

void law8()
{
    group ("law 8: no subnormal state through a note's tail, long silence, or a subnormal-level input");
    const double fs = 48000.0;
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 1));
    o.setParams (levels (1.0f, 1.0f, 1.0f));
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
    o.setParams (levels (1.0f, 0.0f, 1.0f));
    for (std::size_t at = 0; at < l.size(); at += 64)
    {
        float* io[2] = { l.data() + at, r.data() + at };
        felitronics::test::run (o.process (io, 2, 64));
    }
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

    const Buf s = l0;   // the sub alone, mono through the same events
    const Buf y = render (fs, s, levels (1.0f, 0.0f, 0.0f));
    const double hz = crossingHz (y, fs, 2.5, 4.0);
    std::printf ("    after the events: the sub at %.5f Hz; its peak anywhere %.3f\n", hz, peakOf (y, fs, 0.0, 4.0));
    ok (std::fabs (hz / (kLowE / 2.0) - 1.0) < 1.0e-3, "...and the sub is f/2 again after them");
    ok (peakOf (y, fs, 0.0, 4.0) < 2.0 * oc::Octaver::kInputLimit, "...and stayed bounded by the tracker's input limit");
}

void monoStereo()
{
    group ("mono and stereo: the sum is tracked, the same sub goes to both, a dual-mono pair plays as mono");
    const double fs = 48000.0;
    Buf x = pluck (fs, kLowE, 1.0, 0.2, 0.1, 0.96, 0.6);
    const Buf mono = render (fs, x, levels (1.0f, 0.5f, 1.0f));

    Buf l = x, r = x;
    oc::Octaver o;
    felitronics::test::run (o.prepare (fs, 64, 2));
    o.setParams (levels (1.0f, 0.5f, 1.0f));
    for (std::size_t at = 0; at < l.size(); at += 64)
    {
        const int n = (int) std::min<std::size_t> (64, l.size() - at);
        float* io[2] = { l.data() + at, r.data() + at };
        felitronics::test::run (o.process (io, 2, n));
    }
    ok (std::memcmp (l.data(), mono.data(), l.size() * sizeof (float)) == 0
        && std::memcmp (r.data(), mono.data(), r.size() * sizeof (float)) == 0, "L = R: each channel is bit-identical to mono");

    Buf a = x, b = noise (x.size(), 0.2, 99u);
    addInto (b, pluck (fs, 55.0, 1.0, 0.2, 0.1, 0.96, 0.4));
    oc::Octaver s;
    felitronics::test::run (s.prepare (fs, 64, 2));
    s.setParams (levels (1.0f, 0.5f, 0.0f));
    for (std::size_t at = 0; at < a.size(); at += 64)
    {
        const int n = (int) std::min<std::size_t> (64, a.size() - at);
        float* io[2] = { a.data() + at, b.data() + at };
        felitronics::test::run (s.process (io, 2, n));
    }
    ok (std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0 && peakOf (a, fs, 0.2, 1.0) > 0.05,
        "L != R: the same sub on both channels");

    // A narrow call on a stereo stage is the mono stage.
    Buf n1 = x;
    oc::Octaver w;
    felitronics::test::run (w.prepare (fs, 64, 2));
    w.setParams (levels (1.0f, 0.5f, 1.0f));
    renderInto (w, n1);
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
    lock();
    decay();
    harmonics();
    silenceBetweenNotes();
    latency();
    amountZero();
    glide();
    pauseIsSilence();
    resetRestarts();
    law8();
    poison();
    monoStereo();
    prepareBinds();
    noAllocation();
    return felitronics::test::report();
}
