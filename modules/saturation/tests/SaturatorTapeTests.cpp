// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The Tape shape (Saturator.h, WaveShaper.h) — Tanh's core between a pre-emphasis E and its exact inverse D. What each
// group pins:
//
//   * REFERENCE NULL at os 1: the stage against E, the core and D written here in double from the analog formula, on
//     full-scale tones, a square and a click train, at four drives and two rates.
//   * Linearity: a -60 dBFS multitone comes out as Tanh at the same drive lets it through; at the drive floor a
//     full-scale one does too (D·E = 1).
//   * The top saturates first, measured on the fundamental: the 0.5 dB compression point of 10 kHz sits below 1 kHz's
//     by E's own analog gain difference, and Tanh's two sit together. The same 10 kHz point at 44.1, 48 and 96 kHz.
//   * Peaks against Tanh's; latency and a comb-free mix; the bypass below 14.1 kHz at os 1 renders Tanh's bits.
//   * The model state's lifecycle: NaN, the denormal flush, chunking, reset, a channel that leaves and returns, a
//     switch away and back, a DIRECT switch to and from the Transformer, its bytes, and the glide.

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
using Shape = saturation::WaveShaper::Shape;
using Buf = std::vector<std::vector<float>>;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string num (double v, const char* fmt = "%.4g") { char b[64]; std::snprintf (b, sizeof b, fmt, v); return b; }

void okp (bool cond, const std::string& msg) { std::printf ("       %s\n", msg.c_str()); ok (cond, msg); }

// A running maximum a NaN cannot hide from: a NaN wins and STAYS, and every gate compares with <=, which it fails.
template <class T> void worse (T& w, T v) { if (! std::isnan (w) && ! (v <= w)) w = v; }

bool bitEq (float a, float b) { return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b); }

long long diffs (const std::vector<float>& a, const std::vector<float>& b, std::size_t from, std::size_t to)
{
    if (a.size() < to || b.size() < to) return -1;
    long long d = 0;
    for (std::size_t i = from; i < to; ++i) d += bitEq (a[i], b[i]) ? 0 : 1;
    return d;
}

double maxAbsDiff (const std::vector<float>& a, const std::vector<float>& b, std::size_t from = 0)
{
    double w = a.size() == b.size() && a.size() > from ? 0.0 : kNaN;
    for (std::size_t i = from; i < a.size() && i < b.size(); ++i) worse (w, std::fabs ((double) a[i] - (double) b[i]));
    return w;
}

double peakOf (const std::vector<float>& y)
{
    double p = y.empty() ? kNaN : 0.0;
    for (float v : y) worse (p, (double) std::fabs (v));
    return p;
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

// A naive (not band-limited) square at `hz`, ±`amp`.
std::vector<float> square (double fs, double hz, double amp, std::size_t n)
{
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = std::sin (2.0 * kPi * hz * (double) i / fs + 0.1) >= 0.0 ? (float) amp : (float) -amp;
    return x;
}

// A band-limited click train: 100 Hz harmonics in cosine phase up to 18 kHz, scaled to a peak of 1.
std::vector<float> clicks (double fs, std::size_t n)
{
    std::vector<double> d (n, 0.0);
    for (int h = 1; 100.0 * h <= 18000.0; ++h)
        for (std::size_t i = 0; i < n; ++i) d[i] += std::cos (2.0 * kPi * 100.0 * h * (double) i / fs);
    double m = 0.0;
    for (double v : d) m = std::max (m, std::fabs (v));
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = (float) (d[i] / m);
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

// Full-scale material for the NULL: tones 200 Hz..15 kHz (peak 1), then a 1 kHz square (±1), then single-sample ±1
// clicks on silence, a third each.
std::vector<float> fullScale (double fs, std::size_t n)
{
    const std::size_t t = n / 3;
    std::vector<float> x = multitone (fs, t, 200.0, 15000.0, 1.0);
    const std::vector<float> sq = square (fs, 1000.0, 1.0, t);
    x.insert (x.end(), sq.begin(), sq.end());
    for (std::size_t i = 0; x.size() < n; ++i) x.push_back (i % 480 == 0 ? ((i / 480) % 2 ? -1.0f : 1.0f) : 0.0f);
    return x;
}

// Varied programme for the lifecycle: LF and HF tones, a burst of a 3 kHz square, `v` varies the phases.
std::vector<float> programme (double fs, std::size_t n, int v)
{
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / fs, ph = 0.7 * v;
        const double sq = (i / 2000) % 3 == 1 ? (std::sin (2.0 * kPi * 3000.0 * t + ph) >= 0.0 ? 0.3 : -0.3) : 0.0;
        x[i] = (float) (0.35 * std::sin (2.0 * kPi * 47.0 * t + ph) + 0.2 * std::sin (2.0 * kPi * 1300.0 * t + 1.1)
                        + 0.2 * std::sin (2.0 * kPi * 9100.0 * t + 2.0 + ph) + sq);
    }
    return x;
}

// Amplitude of the component at `cycles` cycles per window y[0, n) (Goertzel, double).
double amp (const float* y, std::size_t n, double cycles)
{
    const double w = 2.0 * kPi * cycles / (double) n, c2 = 2.0 * std::cos (w);
    double s1 = 0.0, s2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) { const double s0 = (double) y[i] + c2 * s1 - s2; s2 = s1; s1 = s0; }
    const double pw = s1 * s1 + s2 * s2 - c2 * s1 * s2;
    return 2.0 * std::sqrt (std::max (0.0, pw)) / (double) n;
}

//==============================================================================
// THE ORACLE, from the analog formula and touching no production code:
//   E(s) = (1 + s·τ) / (1 + s·τ/g),   τ = 50 µs,  g = 10^(6/20)   (the corners ω1 = 1/τ and ω2 = g/τ)
//   D(s) = 1 / E(s)
// Each first-order factor (1 + s/ωc) goes through its own bilinear transform s = Kc·(1 − z⁻¹)/(1 + z⁻¹), with Kc chosen
// so that ωc maps onto itself, Kc = ωc / tan(ωc / (2·fs)); the two (1 + z⁻¹) denominators cancel. Run in direct form I.
// The core y = tanh(k·e)/tanh(k), k = max(1e-4, 10^(driveDb/20) − 1), and autoComp 1 multiplies by 1/slope = tanh(k)/k.
std::vector<double> oracle (const std::vector<double>& x, double fs, double driveDb)
{
    const double tau = 50.0e-6, g = std::pow (10.0, 6.0 / 20.0);
    const double w1 = 1.0 / tau, w2 = g / tau;
    const double K1 = w1 / std::tan (w1 / (2.0 * fs)), K2 = w2 / std::tan (w2 / (2.0 * fs));
    // (1 + s/ωc) → [(1 + z⁻¹) + (Kc/ωc)(1 − z⁻¹)] / (1 + z⁻¹)
    const double n0 = 1.0 + K1 / w1, n1 = 1.0 - K1 / w1;     // E's numerator
    const double d0 = 1.0 + K2 / w2, d1 = 1.0 - K2 / w2;     // E's denominator
    const double k = std::max (1.0e-4, std::pow (10.0, driveDb / 20.0) - 1.0);
    const double comp = std::tanh (k) / k;
    std::vector<double> y (x.size());
    double xp = 0.0, ep = 0.0, wp = 0.0, yp = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double e = (n0 * x[i] + n1 * xp - d1 * ep) / d0;
        const double w = std::tanh (k * e) / std::tanh (k);
        const double v = (d0 * w + d1 * wp - n1 * yp) / n0;  // D: numerator and denominator swapped
        xp = x[i]; ep = e; wp = w; yp = v;
        y[i] = comp * v;
    }
    return y;
}

// |E(j·2π·f)| from the analog formula.
double analogE (double f)
{
    const double tau = 50.0e-6, g = std::pow (10.0, 6.0 / 20.0), w = 2.0 * kPi * f;
    return std::sqrt (1.0 + (w * tau) * (w * tau)) / std::sqrt (1.0 + (w * tau / g) * (w * tau / g));
}

//==============================================================================
void referenceNull()
{
    double worstAll = 0.0;
    for (double fs : { 44100.0, 96000.0 })
    {
        const std::vector<float> x = fullScale (fs, (std::size_t) (1.5 * fs));
        const std::vector<double> xd (x.begin(), x.end());
        for (float db : { 0.0f, 6.0f, 12.0f, 24.0f })
        {
            const std::vector<float> y = render (P (Shape::Tape, db), fs, 1, x);
            const std::vector<double> r = oracle (xd, fs, (double) db);
            double worst = y.size() == r.size() ? 0.0 : kNaN;
            for (std::size_t i = 0; i < y.size() && i < r.size(); ++i) worse (worst, std::fabs ((double) y[i] - r[i]));
            std::printf ("       %5.1f kHz, driveDb %4.1f: max |stage - double model| %.3g\n", fs / 1000.0, (double) db, worst);
            worse (worstAll, worst);
        }
    }
    okp (worstAll <= 1.0e-5, "os 1, mix 1, autoComp 1: the stage nulls against the double model to <= 1e-5 at driveDb 0/6/12/24, "
                             "44.1 and 96 kHz (worst " + num (worstAll, "%.3g") + ")");
}

void linearity()
{
    const double fs = 48000.0;
    const std::vector<float> quiet = multitone (fs, 48000, 20.0, 20000.0, 1.0e-3);
    for (float db : { 6.0f, 12.0f })
    {
        const double d = maxAbsDiff (render (P (Shape::Tape, db), fs, 4, quiet), render (P (Shape::Tanh, db), fs, 4, quiet));
        okp (d <= 1.0e-6, "-60 dBFS multitone 20 Hz..20 kHz, os 4, driveDb " + num (db) + ": Tape against Tanh, max |diff| "
                          + num (d, "%.3g") + " <= 1e-6");
    }
    const std::vector<float> loud = multitone (fs, 48000, 20.0, 20000.0, 1.0);
    const double d = maxAbsDiff (render (P (Shape::Tape, 0.0f), fs, 4, loud), render (P (Shape::Tanh, 0.0f), fs, 4, loud));
    okp (d <= 1.0e-5, "full-scale multitone, os 4, driveDb 0 (k = 1e-4): Tape against Tanh, max |diff| " + num (d, "%.3g")
                      + " <= 1e-5 (D*E = 1)");
}

//==============================================================================
// The fundamental's gain through the stage (os 4, driveDb 6, autoComp 0) for a sine at `levelDb` dBFS: 0.1 s measured
// after 0.1 s, bin-exact at every rate used here.
double fundGain (Shape s, double fs, double hz, double levelDb)
{
    const std::size_t n = (std::size_t) (fs / 10.0);
    const double a = std::pow (10.0, levelDb / 20.0);
    const std::vector<float> y = render (P (s, 6.0f, 1.0f, 0.0f), fs, 4, sine (fs, hz, a, 2 * n));
    return y.size() == 2 * n ? amp (y.data() + n, n, hz * (double) n / fs) / a : kNaN;
}

// The input level (dBFS) at which the fundamental is compressed by 0.5 dB against its gain at -60 dBFS. NaN if any
// measurement was.
double compressionPoint (Shape s, double fs, double hz)
{
    const double g0 = fundGain (s, fs, hz, -60.0);
    double lo = -40.0, hi = 12.0;
    bool nan = std::isnan (g0);
    for (int it = 0; it < 30; ++it)
    {
        const double mid = 0.5 * (lo + hi), c = 20.0 * std::log10 (fundGain (s, fs, hz, mid) / g0);
        nan = nan || std::isnan (c);
        (c > -0.5 ? lo : hi) = mid;
    }
    return nan ? kNaN : 0.5 * (lo + hi);
}

void hfFirst()
{
    const double fs = 48000.0;
    const double want = 20.0 * std::log10 (analogE (10000.0) / analogE (1000.0));
    const double t1 = compressionPoint (Shape::Tape, fs, 1000.0), t10 = compressionPoint (Shape::Tape, fs, 10000.0);
    const double h1 = compressionPoint (Shape::Tanh, fs, 1000.0), h10 = compressionPoint (Shape::Tanh, fs, 10000.0);
    std::printf ("       48 kHz, os 4, driveDb 6, autoComp 0: 0.5 dB compression of the fundamental at (dBFS)\n");
    std::printf ("       Tape 1 kHz %.3f, 10 kHz %.3f | Tanh 1 kHz %.3f, 10 kHz %.3f\n", t1, t10, h1, h10);
    okp (std::fabs ((t1 - t10) - want) <= 0.2, "Tape: 10 kHz compresses " + num (t1 - t10, "%.3f") + " dB before 1 kHz, E's analog "
                                               "|E(10k)|/|E(1k)| is " + num (want, "%.3f") + " dB (within 0.2)");
    okp (std::fabs (h1 - h10) <= 0.1, "Tanh (control): 10 kHz and 1 kHz compress at the same level (" + num (h1 - h10, "%.3f") + " dB, within 0.1)");
}

void rateIndependence()
{
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    bool nan = false;
    for (double fs : { 44100.0, 48000.0, 96000.0 })
    {
        const double c = compressionPoint (Shape::Tape, fs, 10000.0);
        std::printf ("       %5.1f kHz, os 4: Tape's 10 kHz 0.5 dB compression point %.3f dBFS\n", fs / 1000.0, c);
        nan = nan || std::isnan (c);
        lo = std::min (lo, c); hi = std::max (hi, c);
    }
    const double spread = nan ? kNaN : hi - lo;
    okp (spread <= 0.1, "44.1 / 48 / 96 kHz within 0.1 dB (spread " + num (spread, "%.4f") + " dB)");
}

//==============================================================================
void peaks()
{
    const double fs = 48000.0;
    const std::size_t n = (std::size_t) (0.5 * fs);
    struct M { const char* name; std::vector<float> x; };
    const std::vector<M> ms { { "1 kHz sine", sine (fs, 1000.0, 1.0, n) }, { "5 kHz sine", sine (fs, 5000.0, 1.0, n) },
                              { "10 kHz sine", sine (fs, 10000.0, 1.0, n) }, { "1 kHz square", square (fs, 1000.0, 1.0, n) },
                              { "band-limited clicks", clicks (fs, n) } };
    std::printf ("       os 4, autoComp 1: 20*log10(max|out| / max|in|) in dB, Tape | Tanh\n");
    bool within = true;
    for (const M& m : ms)
    {
        const double in = peakOf (m.x);
        std::printf ("       %-20s", m.name);
        for (float db : { 0.0f, 6.0f, 12.0f })
        {
            const double tp = 20.0 * std::log10 (peakOf (render (P (Shape::Tape, db), fs, 4, m.x)) / in);
            const double th = 20.0 * std::log10 (peakOf (render (P (Shape::Tanh, db), fs, 4, m.x)) / in);
            std::printf (" | %2.0f dB: %+.3f %+.3f", (double) db, tp, th);
            within = within && tp <= th + 0.5;
        }
        std::printf ("\n");
    }
    okp (within, "Tape's peak <= Tanh's + 0.5 dB on every row");
}

void mixAndLatency()
{
    bool same = true;
    for (int os : { 1, 2, 4, 8 })
    {
        Sat t, h;
        t.setParams (P (Shape::Tape, 6.0f)); h.setParams (P (Shape::Tanh, 6.0f));
        same = same && prep (t, 48000.0, 2, os) && prep (h, 48000.0, 2, os) && t.latencySamples() == h.latencySamples();
    }
    ok (same, "latencySamples() equals Tanh's at os 1, 2, 4 and 8");
    const double fs = 48000.0;
    const std::vector<float> x = multitone (fs, 48000, 20.0, 12000.0, 1.0e-3);
    Sat s; s.setParams (P (Shape::Tape, 6.0f, 0.5f));
    if (! prep (s, fs, 1, 4)) return;
    const std::size_t lat = (std::size_t) s.latencySamples();
    Buf y { x };
    through (s, y, 512);
    double worst = lat > 0 ? 0.0 : kNaN, comb = 0.0;
    for (std::size_t i = lat + 2048; i < x.size(); ++i)
    {
        worse (worst, std::fabs ((double) y[0][i] - (double) x[i - lat]));
        worse (comb, std::fabs ((double) y[0][i] - (double) x[i]));
    }
    okp (worst <= 1.0e-6, "mix 0.5, -60 dBFS multitone 20 Hz..12 kHz, os 4: nulls against the input delayed by the latency ("
                          + std::to_string (lat) + " samples) to " + num (worst, "%.3g") + " (against the undelayed input: "
                          + num (comb, "%.3g") + ")");
}

// Below the bypass rate (f2 = 6351 Hz >= 0.45·fsOs, os 1 under ~14.1 kHz) Tape is Tanh bit for bit, settled and
// gliding; at 16 kHz os 1 the pair runs and the two differ.
void bypass()
{
    const double f2 = Sat::kTapeCornerHz * std::pow (10.0, Sat::kTapeShelfDb / 20.0);
    std::printf ("       f1 %.1f Hz, f2 %.1f Hz: the pair is bypassed at os 1 below %.1f Hz\n", Sat::kTapeCornerHz, f2,
                 f2 / Sat::kTapeMaxCorner);
    for (double fs : { 8000.0, 16000.0 })
    {
        const std::size_t n = (std::size_t) fs;
        std::vector<float> ya[2], yb[2];
        for (int j = 0; j < 2; ++j)
        {
            Sat s; s.setParams (P (j == 0 ? Shape::Tape : Shape::Tanh, 6.0f));
            if (! prep (s, fs, 1, 1)) return;
            Buf y { programme (fs, n, 3) };
            through (s, y, 256, 0, n / 2);
            s.setParams (P (j == 0 ? Shape::Tape : Shape::Tanh, 14.0f, 0.9f));
            through (s, y, 256, n / 2, n);
            (j == 0 ? ya : yb)[0] = y[0];
        }
        const long long d = diffs (ya[0], yb[0], 0, n);
        if (fs == 8000.0) okp (d == 0, "8 kHz, os 1: Tape renders Tanh's bits, settled and through a drive glide (" + std::to_string (d) + " differ)");
        else              okp (d > 0, "control, 16 kHz, os 1: the pair runs and Tape differs from Tanh (" + std::to_string (d) + " samples)");
    }
}

//==============================================================================
void lifecycle()
{
    const double fs = 48000.0;
    // NaN, then silence: the gate keeps the model finite, and the flush takes both states to exact zero.
    for (int os : { 1, 4 })
        for (float db : { 0.0f, 12.0f })
        {
            std::vector<float> x = programme (fs, 24000, 1);
            x.insert (x.end(), 1000, std::numeric_limits<float>::quiet_NaN());
            x.insert (x.end(), (std::size_t) fs, 0.0f);
            const std::vector<float> y = render (P (Shape::Tape, db), fs, os, x);
            bool finite = ! y.empty(), zero = ! y.empty();
            for (float v : y) finite = finite && std::isfinite (v);
            for (std::size_t i = y.size() > 4800 ? y.size() - 4800 : 0; i < y.size(); ++i) zero = zero && y[i] == 0.0f;
            ok (finite && zero, "os " + std::to_string (os) + ", driveDb " + num (db) + ": a NaN burst leaves the output finite, "
                                "and after 1 s of silence it is exact 0.0f");
        }

    // Chunking, with a drive glide in the middle of the stream: n = 1, 7, 64, 4096 render the same bits.
    {
        const std::size_t n = 86016, at = 28672;              // both multiples of 7 · 4096
        Buf x { programme (fs, n, 1), programme (fs, n, 2) };
        Buf ref;
        long long bad = 0;
        for (int blk : { 4096, 1, 7, 64 })
        {
            Sat s; s.setParams (P (Shape::Tape, 6.0f));
            if (! prep (s, fs, 2, 4)) return;
            Buf y = x;
            through (s, y, blk, 0, at);
            s.setParams (P (Shape::Tape, 15.0f, 0.8f));
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
        // reset() is a fresh prepare(), bit for bit — after a loud passage that left both states away from 0.
        {
            Sat a, b;
            a.setParams (P (Shape::Tape, 12.0f)); b.setParams (P (Shape::Tape, 12.0f));
            if (! prep (a, fs, 1, os) || ! prep (b, fs, 1, os)) return;
            Buf warm { programme (fs, n, 3) }, ya { programme (fs, n, 4) }, yb = ya;
            through (a, warm, 512);
            a.reset();
            through (a, ya, 512);
            through (b, yb, 512);
            const long long d = diffs (ya[0], yb[0], 0, n);
            ok (d == 0, "reset() renders what a fresh prepare() renders" + tag + " (" + std::to_string (d) + " differ)");
        }
        // A channel that leaves and returns starts from zero: from the return on, it is a fresh stage's channel.
        {
            Sat a, b;
            a.setParams (P (Shape::Tape, 12.0f, 0.7f)); b.setParams (P (Shape::Tape, 12.0f, 0.7f));
            if (! prep (a, fs, 2, os) || ! prep (b, fs, 2, os)) return;
            Buf ya { programme (fs, 3 * n, 5), programme (fs, 3 * n, 6) };
            through (a, ya, 512, 0, n);
            through (a, ya, 512, n, 2 * n, 1);                          // channel 1 is away
            through (a, ya, 512, 2 * n, 3 * n);
            Buf yb { std::vector<float> (ya[0].size()), programme (fs, 3 * n, 6) };
            through (b, yb, 512, 2 * n, 3 * n);
            const long long d = diffs (ya[1], yb[1], 2 * n, 3 * n);
            ok (d == 0, "a channel that leaves and returns renders what a fresh stage renders from the return on" + tag
                        + " (" + std::to_string (d) + " differ)");
        }
        // Tape -> Tanh -> Tape starts from zero: the same as a stage that never ran Tape.
        {
            Sat a, b;
            a.setParams (P (Shape::Tape, 12.0f)); b.setParams (P (Shape::Tanh, 12.0f));
            if (! prep (a, fs, 1, os) || ! prep (b, fs, 1, os)) return;
            Buf ya { programme (fs, 3 * n, 7) }, yb = ya;
            through (a, ya, 512, 0, n);
            through (b, yb, 512, 0, n);
            a.setParams (P (Shape::Tanh, 12.0f));
            through (a, ya, 512, n, 2 * n);
            through (b, yb, 512, n, 2 * n);
            a.setParams (P (Shape::Tape, 12.0f)); b.setParams (P (Shape::Tape, 12.0f));
            through (a, ya, 512, 2 * n, 3 * n);
            through (b, yb, 512, 2 * n, 3 * n);
            const long long d = diffs (ya[0], yb[0], 2 * n, 3 * n);
            ok (d == 0, "Tape -> Tanh -> Tape renders what a stage that never ran Tape renders" + tag
                        + " (" + std::to_string (d) + " differ)");
        }
    }

    // A DIRECT switch between the two model shapes keeps the model running, so it must hand the new model zeroed slots.
    // At os 1 (no oversampler memory, no dry delay) the switched stage then IS a fresh stage of the new shape fed only
    // what came after the switch.
    for (const Shape from : { Shape::Transformer, Shape::Tape })
    {
        const Shape to = from == Shape::Tape ? Shape::Transformer : Shape::Tape;
        const std::size_t n = 14400;
        Sat a, b;
        a.setParams (P (from, 12.0f)); b.setParams (P (to, 12.0f));
        if (! prep (a, fs, 2, 1) || ! prep (b, fs, 2, 1)) return;
        Buf ya { programme (fs, 2 * n, 8), programme (fs, 2 * n, 9) }, yb = ya;
        through (a, ya, 512, 0, n);
        a.setParams (P (to, 12.0f));
        through (a, ya, 512, n, 2 * n);
        through (b, yb, 512, n, 2 * n);
        long long d = 0;
        for (int c = 0; c < 2; ++c) d += diffs (ya[(std::size_t) c], yb[(std::size_t) c], n, 2 * n);
        ok (d == 0, std::string (from == Shape::Tape ? "Tape -> Transformer" : "Transformer -> Tape")
                    + ", directly, os 1: from the switch on, a fresh stage of the new shape (" + std::to_string (d) + " differ)");
    }

    // The model state's bytes: Storage says them, prepare() asks for exactly that, and they are kModelFloats (2) floats
    // per channel — one more per channel than the Transformer alone needed.
    {
        bool okAll = Sat::kModelFloats == 2;
        const int taps = oversampling::PolyphaseOversampler::kDefaultTapsPerPhase;
        struct G { double fs; int mb, ch, os; };
        for (const G& g : { G { 48000.0, 512, 2, 4 }, G { 44100.0, 256, 1, 1 }, G { 96000.0, 300, 6, 2 } })
        {
            Sat::Storage st;
            okAll = okAll && Sat::storageFor (g.fs, g.mb, g.ch, g.os, taps, st);
            Sat::Storage without = st; without.model = 0;
            Sat s; s.setParams (P (Shape::Tape, 6.0f));
            const long long before = alloc::bytes.load();
            const bool okP = s.prepare (g.fs, g.mb, g.ch, g.os, taps);
            const long long got = alloc::bytes.load() - before;
            const bool row = okP && got == (long long) st.bytes() && st.model == (std::size_t) (g.ch * 2)
                          && st.bytes() - without.bytes() == sizeof (float) * (std::uint64_t) (g.ch * 2);
            std::printf ("       fs %.0f, block %d, %d ch, os %d: Storage::bytes() %llu, prepare() asked for %lld\n",
                         g.fs, g.mb, g.ch, g.os, (unsigned long long) st.bytes(), got);
            okAll = okAll && row;
        }
        ok (okAll, "Storage::bytes() == what prepare() allocated, and the model state is 2 floats per channel of it");
    }
}

//==============================================================================
// A drive glide lands on the settled floats. The pre-emphasis depends on the input alone, but the de-emphasis state
// carries what the gliding core put into it, so the two stages meet once that decays: gated bit-identical from the
// landing plus `settle` samples, and the first sample from which they agree is printed.
void glideLands()
{
    const double fs = 48000.0;
    const std::size_t n = 16000, at = 3000, settle = 256;
    for (int os : { 1, 4 })
    {
        const std::string tag = " (os " + std::to_string (os) + ")";
        const Sat::Params a = P (Shape::Tape, 3.0f), b = P (Shape::Tape, 15.0f, 0.7f, 0.8f, -2.0f);
        Sat m, t;
        m.setParams (a); t.setParams (b);
        if (! prep (m, fs, 2, os) || ! prep (t, fs, 2, os)) return;
        const Buf x { programme (fs, n, 8), programme (fs, n, 9) };
        Buf ym = x, yt = x;
        through (m, ym, 256, 0, at);
        m.setParams (b);
        through (m, ym, 256, at, n);
        through (t, yt, 256);
        const std::size_t start = (at / 64 + 1) * 64;
        const std::size_t landed = start + (std::size_t) m.glideTicks() * 64 + (std::size_t) m.latencySamples() + 1;
        long long after = 0, during = 0;
        std::size_t agree = landed;
        for (int c = 0; c < 2; ++c)
        {
            after  += diffs (ym[(std::size_t) c], yt[(std::size_t) c], landed + settle, n);
            during += diffs (ym[(std::size_t) c], yt[(std::size_t) c], start, landed - 64);
            for (std::size_t i = landed; i < n; ++i)
                if (! bitEq (ym[(std::size_t) c][i], yt[(std::size_t) c][i])) agree = std::max (agree, i + 1);
        }
        okp (m.glideTicks() > 0 && after == 0 && during > 0, "a drive glide lands bit-identical to the target stage" + tag
                                                             + ": from " + std::to_string (agree - landed) + " samples after the landing ("
                                                             + std::to_string (after) + " differ past landing + " + std::to_string (settle)
                                                             + ", " + std::to_string (during) + " during the glide)");
    }
}

// The glide path runs the same model: a 6 -> 6.001 dB glide on full-scale material stays next to the settled Tape at
// 6 dB throughout, where Tanh at the same drive is far from it.
void glideIsTape()
{
    const double fs = 48000.0;
    const std::size_t n = 12000, at = 2000;
    const std::vector<float> x = fullScale (fs, n);
    double glided = 0.0, gap = std::numeric_limits<double>::infinity();
    for (int os : { 1, 4 })
    {
        Sat g; g.setParams (P (Shape::Tape, 6.0f));
        if (! prep (g, fs, 1, os)) return;
        Buf y { x };
        through (g, y, 256, 0, at);
        g.setParams (P (Shape::Tape, 6.001f));
        const bool gliding = g.isGliding();
        through (g, y, 256, at, n);
        const std::vector<float> settled = render (P (Shape::Tape, 6.0f), fs, os, x, 256);
        const std::vector<float> th = render (P (Shape::Tanh, 6.0f), fs, os, x, 256);
        const double dg = gliding ? maxAbsDiff (y[0], settled, at) : kNaN, dt = maxAbsDiff (th, settled, at);
        worse (glided, dg);
        gap = std::isnan (gap) || std::isnan (dt) ? kNaN : std::min (gap, dt);   // the control's smallest gap
    }
    okp (glided <= 1.0e-3 && gap >= 1.0e-1, "os 1 and 4: a 6 -> 6.001 dB glide stays within " + num (glided, "%.3g")
                                             + " of the settled Tape at 6 dB (<= 1e-3); Tanh at 6 dB is " + num (gap, "%.3g")
                                             + " from it (>= 0.1)");
}
} // namespace

int main()
{
    std::printf ("felitronics::saturation — Tape\n");
    group ("reference NULL at os 1 against E, the core and D in double");
    referenceNull();
    group ("linearity: Tanh's small-signal path, and D*E = 1");
    linearity();
    group ("the top saturates first, on the fundamental");
    hfFirst();
    group ("rate independence");
    rateIndependence();
    group ("peaks against Tanh's");
    peaks();
    group ("mix and latency");
    mixAndLatency();
    group ("the bypass below 14.1 kHz at os 1");
    bypass();
    group ("the model state's lifecycle");
    lifecycle();
    group ("the glide");
    glideLands();
    glideIsTape();
    return felitronics::test::report();
}
