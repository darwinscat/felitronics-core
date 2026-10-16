// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// dehum::HumRemover knowing the mains: the FLL on detuned and drifting mains, Auto 50/60 on every hum shape
// (full, missing fundamental, odd-only, even-only, opposite-phase stereo), stickiness while playing, the
// false-learning fixtures (held notes beside harmonics, vibrato, a bend, a chord — and the tracker must not
// move under them), plucks beside harmonics, steps at mode switches, and the ghost a note leaves.
//
// THE CLEAN PASS. "Within x dB of the clean pass" compares the stage, fed hum + note, with the same stage
// fed the note alone in its LTI configuration (dehumtest::lti: the playing width, no tracking) — what a fixed
// 0.5 Hz notch comb does to the note. The difference is what the hum, the learner and the tracker did.
#include <felitronics_test.h>
#include "dehum_fixtures.h"

#include <cstdio>
#include <cstring>
#include <string>

using namespace dehumtest;
namespace dh = felitronics::dehum;
using felitronics::test::ok;
using felitronics::test::group;
using HR = dh::HumRemover<>;

//------------------------------------------------------------------------------------------------------
struct HumSpec
{
    double f0 = 50.0, a1 = 0.01, slope = 0.0, rampFrom = 0.0, rampTo = 0.0, sign = 1.0, until = 1e9;
    int    kmin = 1, kmax = 24, parity = 0;   // parity 1 = odd harmonics only, 2 = even only
    std::uint32_t seed = 1;
};
static double humCycles (const HumSpec& s, double t)
{
    double c = s.f0 * t;
    if (s.slope != 0.0 && t > s.rampFrom)
    {
        const double te = std::min (t, s.rampTo) - s.rampFrom;
        c += 0.5 * s.slope * te * te;
        if (t > s.rampTo) c += s.slope * (s.rampTo - s.rampFrom) * (t - s.rampTo);
    }
    return c;
}
static double humFreq (const HumSpec& s, double t)
{
    return t <= s.rampFrom ? s.f0 : s.f0 + s.slope * (std::min (t, s.rampTo) - s.rampFrom);
}
static void addHum (std::vector<float>& x, const HumSpec& s, double fs, std::int64_t from = 0)
{
    Lcg r { s.seed }; double ph[64];
    for (double& p : ph) p = 2.0 * kPi * r.uniform();
    for (std::int64_t n = from; n < (std::int64_t) x.size(); ++n)
    {
        const double t = (double) n / fs;
        if (t >= s.until) break;
        const double c = humCycles (s, t);
        double v = 0.0;
        for (int k = s.kmin; k <= s.kmax; ++k)
        {
            if ((s.parity == 1 && k % 2 == 0) || (s.parity == 2 && k % 2 == 1)) continue;
            const double kc = k * c;
            v += s.a1 / k * std::sin (2.0 * kPi * (kc - std::floor (kc)) + ph[k - 1]);
        }
        x[(std::size_t) n] += (float) (s.sign * v);
    }
}

struct Note
{
    double f0 = 0.0, amp = 0.1, B = 1e-4, start = 0.0, dur = 1.0, tau = 0.0, attack = 0.010;
    double vibRate = 0.0, vibDepth = 0.0, bendRatio = 1.0, bendSec = 0.0, release = 0.010;
    int    partials = 6;
    double partialHz (int p) const { return p * f0 * std::sqrt (1.0 + B * p * p); }
};
static void addNote (std::vector<float>& x, const Note& nt, double fs)
{
    const std::int64_t n0 = (std::int64_t) std::llround (nt.start * fs), n1 = std::min<std::int64_t> ((std::int64_t) x.size(), (std::int64_t) std::llround ((nt.start + nt.dur) * fs));
    double cyc[16] {};
    for (std::int64_t n = n0; n < n1; ++n)
    {
        const double t = (double) (n - n0) / fs, left = (double) (n1 - n) / fs;
        double ratio = 1.0;
        if (nt.vibRate > 0.0) ratio *= 1.0 + nt.vibDepth * std::sin (2.0 * kPi * nt.vibRate * t);
        if (nt.bendSec > 0.0) ratio *= t < nt.bendSec ? 1.0 + (nt.bendRatio - 1.0) * t / nt.bendSec : nt.bendRatio;
        double env = t < nt.attack ? 0.5 - 0.5 * std::cos (kPi * t / nt.attack) : 1.0;
        if (left < nt.release) env *= 0.5 - 0.5 * std::cos (kPi * left / nt.release);
        if (nt.tau > 0.0) env *= std::exp (-t / nt.tau);
        double v = 0.0;
        for (int p = 1; p <= nt.partials; ++p)
        {
            cyc[p] += nt.partialHz (p) * ratio / fs;
            cyc[p] -= std::floor (cyc[p]);
            v += nt.amp / p * std::sin (2.0 * kPi * cyc[p] + 0.3 * p);
        }
        x[(std::size_t) n] += (float) (env * v);
    }
}

// Telemetry per block of the run.
struct Tel { std::vector<double> f; std::vector<int> det, quiet, run; int block = 256; };
static void render (HR& h, std::vector<float>& l, std::vector<float>* r, Tel* tel, int block = 256)
{
    if (tel) tel->block = block;
    for (std::size_t o = 0; o < l.size(); o += (std::size_t) block)
    {
        const int n = (int) std::min<std::size_t> ((std::size_t) block, l.size() - o);
        float* io[2] { l.data() + o, r ? r->data() + o : nullptr };
        felitronics::test::run (h.process (io, r ? 2 : 1, n));
        if (tel)
        {
            tel->f.push_back (h.mainsHz());
            tel->det.push_back (h.detected() == dh::Mains::Hz50 ? 50 : (h.detected() == dh::Mains::Hz60 ? 60 : 0));
            tel->quiet.push_back (h.quiet() ? 1 : 0);
            tel->run.push_back (h.running() ? 1 : 0);
        }
    }
}
static double tOf (const Tel& t, std::size_t i, double fs) { return (double) ((i + 1) * (std::size_t) t.block) / fs; }

static std::vector<float> processed (const std::vector<float>& in, const dh::HumRemoverParams& p, double fs, Tel* tel = nullptr)
{
    HR h; h.setParams (p); felitronics::test::run (h.prepare (fs, 512, 1));
    std::vector<float> x = in;
    render (h, x, nullptr, tel);
    return x;
}

static double rmsDiffDb (const std::vector<float>& a, const std::vector<float>& b, std::int64_t from, std::int64_t to)
{
    double e = 0.0, r = 0.0;
    for (std::int64_t i = from; i < to; ++i)
    {
        const double d = (double) a[(std::size_t) i] - (double) b[(std::size_t) i];
        e += d * d;
        r += (double) b[(std::size_t) i] * (double) b[(std::size_t) i];
    }
    return 10.0 * std::log10 (std::max (e, 1e-300) / std::max (r, 1e-300));
}

static double worstNull (const std::vector<float>& y, const std::vector<float>& x, std::int64_t from, std::int64_t len, double f0, double fs)
{
    std::vector<double> fr;
    for (int k = 1; k <= 24; ++k) fr.push_back (f0 * k);
    const auto ay = fitAmps (y, from, len, fr, fs), ax = fitAmps (x, from, len, fr, fs);
    double w = -400.0;
    for (std::size_t k = 0; k < fr.size(); ++k) w = std::max (w, dbOf (ay[k] / ax[k]));
    return w;
}

static std::string fmt (double v, int prec = 3) { char b[32]; std::snprintf (b, sizeof b, "%.*f", prec, v); return b; }

//------------------------------------------------------------------------------------------------------
// DETUNED MAINS: Auto, hum at 49.6 / 50.4 / 59.5 / 60.4 Hz (-30 dBFS, 24 harmonics). |f̂ − f| < 0.01 Hz by
// 10 s, and k = 24 nulled by >= 40 dB. The coarse estimate alone must bring it within 0.05 Hz by 3 s.
static void testDetuned()
{
    group ("FLL: detuned mains -> |f - f_true| < 0.01 Hz by 10 s, k24 <= -40 dB");
    struct C { double fs, f; };
    for (C c : { C { 48000.0, 49.6 }, C { 48000.0, 50.4 }, C { 48000.0, 59.5 }, C { 48000.0, 60.4 }, C { 192000.0, 49.6 } })
    {
        HumSpec hs; hs.f0 = c.f; hs.a1 = 0.0316; hs.seed = 3u;
        std::vector<float> x ((std::size_t) (12.0 * c.fs), 0.0f);
        addHum (x, hs, c.fs);
        Tel tel;
        const auto y = processed (x, params (dh::Mains::Auto, 0.5, 3.0), c.fs, &tel);
        double worstLate = 0.0, at3 = 1e9, lockT = 0.0;
        for (std::size_t i = 0; i < tel.f.size(); ++i)
        {
            const double t = tOf (tel, i, c.fs), e = std::fabs (tel.f[i] - c.f);
            if (t >= 10.0) worstLate = std::max (worstLate, e);
            if (t >= 3.0 && at3 > 1e8) at3 = e;
            if (e >= 0.01) lockT = t;
        }
        std::vector<double> fr;
        for (int k = 1; k <= 24; ++k) fr.push_back (k * c.f);
        const std::int64_t S = (std::int64_t) (11.0 * c.fs), W = (std::int64_t) c.fs - 1000;
        const auto ao = fitAmps (y, S, W, fr, c.fs), ai = fitAmps (x, S, W, fr, c.fs);
        double worstNull = -400.0;
        for (std::size_t k = 0; k < fr.size(); ++k) worstNull = std::max (worstNull, dbOf (ao[k] / ai[k]));
        const double k24 = dbOf (ao[23] / ai[23]);
        const std::string tag = std::to_string ((int) c.fs) + " Hz, mains " + fmt (c.f, 1);
        std::printf ("      %s: locked (< 0.01 Hz for good) at %.2f s; |err| at 3 s %.4f, after 10 s <= %.5f Hz; k24 %.1f dB, worst harmonic %.1f dB\n",
                     tag.c_str(), lockT, at3, worstLate, k24, worstNull);
        ok (worstLate < 0.01, tag + ": |f - f_true| < 0.01 Hz from 10 s on (" + fmt (worstLate, 5) + ")");
        ok (at3 < 0.05, tag + ": the coarse estimate has it within 0.05 Hz by 3 s (" + fmt (at3, 4) + ")");
        ok (k24 <= -40.0, tag + ": k24 nulled to " + fmt (k24, 1) + " dB");
    }
}

// HISS: hum at 50.3 Hz, -40 dBFS, under white noise at -60 dBFS RMS. The estimate's jitter once locked is
// what the FLL's damping buys; printed, and bounded by what keeps k24 within a tenth of the learn width.
static void testHiss()
{
    group ("FLL under hiss: the estimate's jitter once locked");
    const double fs = 48000.0, f = 50.3;
    HumSpec hs; hs.f0 = f; hs.a1 = 0.01; hs.seed = 5u;
    std::vector<float> x ((std::size_t) (15.0 * fs), 0.0f);
    addHum (x, hs, fs);
    Lcg rng { 77u };
    for (float& v : x) v += (float) (0.001 * std::sqrt (3.0)) * rng.bipolar();
    Tel tel;
    processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
    double e2 = 0.0, worst = 0.0; int n = 0;
    for (std::size_t i = 0; i < tel.f.size(); ++i)
        if (tOf (tel, i, fs) >= 5.0) { const double e = tel.f[i] - f; e2 += e * e; worst = std::max (worst, std::fabs (e)); ++n; }
    const double rmsErr = std::sqrt (e2 / std::max (n, 1));
    std::printf ("      50.3 Hz at -40 dBFS + hiss at -60 dBFS: |f - f_true| rms %.5f Hz, worst %.5f Hz (from 5 s on)\n", rmsErr, worst);
    ok (worst < 0.0125, "under hiss the estimate stays within 0.0125 Hz (k24 within 0.3 Hz, a tenth of the learn width): " + fmt (worst, 5));
}

// SLOW MAINS DRIFT: 50 Hz, then +-0.05 Hz over 20 s, then steady — tracked within 0.01 Hz throughout.
static void testRamp()
{
    group ("FLL: a slow mains ramp (+-0.05 Hz over 20 s) is tracked within 0.01 Hz");
    const double fs = 48000.0;
    for (double slope : { 0.0025, -0.0025 })
    {
        HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.slope = slope; hs.rampFrom = 4.0; hs.rampTo = 24.0; hs.seed = 9u;
        std::vector<float> x ((std::size_t) (26.0 * fs), 0.0f);
        addHum (x, hs, fs);
        Tel tel;
        processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        double worst = 0.0;
        for (std::size_t i = 0; i < tel.f.size(); ++i)
        {
            const double t = tOf (tel, i, fs);
            if (t >= 3.0) worst = std::max (worst, std::fabs (tel.f[i] - humFreq (hs, t)));
        }
        std::printf ("      ramp %+.3f Hz over 20 s: worst |f - f_true| from 3 s on %.5f Hz\n", slope * 20.0, worst);
        ok (worst < 0.01, "ramp " + fmt (slope * 20.0, 3) + " Hz: tracked within 0.01 Hz (" + fmt (worst, 5) + ")");
    }
}

//------------------------------------------------------------------------------------------------------
// AUTO on every hum shape: decided within 1.5 s, never the other family.
static void testAutoFamilies()
{
    group ("Auto: hum-only decided in <= 1.5 s — full, missing fundamental, odd-only, even-only, opposite-phase stereo");
    const double fs = 48000.0;
    struct C { const char* name; double f0; int kmin, parity; bool opposite; };
    const C cases[] = { { "full", 50, 1, 0, false }, { "full", 60, 1, 0, false },
                        { "missing fundamental", 50, 2, 0, false }, { "missing fundamental", 60, 2, 0, false },
                        { "odd only", 50, 1, 1, false }, { "odd only", 60, 1, 1, false },
                        { "even only", 50, 1, 2, false }, { "even only", 60, 1, 2, false },
                        { "opposite-phase stereo", 50, 1, 0, true }, { "opposite-phase stereo", 60, 1, 0, true } };
    for (const C& c : cases)
    {
        HumSpec hs; hs.f0 = c.f0; hs.a1 = 0.01; hs.kmin = c.kmin; hs.parity = c.parity; hs.seed = 21u;
        std::vector<float> l ((std::size_t) (3.0 * fs), 0.0f), r = l;
        addHum (l, hs, fs);
        hs.sign = c.opposite ? -1.0 : 1.0;
        addHum (r, hs, fs);
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 2));
        Tel tel;
        render (h, l, &r, &tel);
        const int want = (int) c.f0;
        double tDecide = -1.0; bool wrong = false;
        for (std::size_t i = 0; i < tel.det.size(); ++i)
        {
            if (tel.det[i] != 0 && tel.det[i] != want) wrong = true;
            if (tDecide < 0.0 && tel.det[i] == want) tDecide = tOf (tel, i, fs);
        }
        const std::string tag = std::string (c.name) + " " + std::to_string (want) + " Hz";
        std::printf ("      %-32s decided at %.3f s, running at %.3f Hz\n", tag.c_str(), tDecide, h.mainsHz());
        ok (tDecide > 0.0 && tDecide <= 1.5 && ! wrong, tag + ": decided " + std::to_string (want) + " in " + fmt (tDecide) + " s, never the other");
        ok (h.running() && std::fabs (h.mainsHz() - c.f0) < 0.05, tag + ": running at the mains frequency");
    }
}

// A FIXED MODE ON THE WRONG FAMILY (the user's error) learns nothing it should not: its coarse estimate is
// taken only where its family wins, so it stays at its own nominal frequency.
static void testWrongFamily()
{
    group ("a fixed mode on the other family's hum stays at its nominal frequency");
    const double fs = 48000.0;
    for (double f : { 50.0, 60.0 })
    {
        HumSpec hs; hs.f0 = f; hs.a1 = 0.01; hs.seed = 8u;
        std::vector<float> x ((std::size_t) (5.0 * fs), 0.0f);
        addHum (x, hs, fs);
        const dh::Mains m = f > 55.0 ? dh::Mains::Hz50 : dh::Mains::Hz60;
        const double nominal = f > 55.0 ? 50.0 : 60.0;
        Tel tel;
        processed (x, params (m, 0.5, 3.0), fs, &tel);
        double dev = 0.0;
        for (double v : tel.f) dev = std::max (dev, std::fabs (v - nominal));
        ok (dev <= 0.05, "mode " + fmt (nominal, 0) + " on " + fmt (f, 0) + " Hz hum: the estimate stays within 0.05 Hz of "
                         + fmt (nominal, 0) + " (" + fmt (dev, 4) + ")");
    }
}

//------------------------------------------------------------------------------------------------------
// B1 (61.7 Hz) at -20 dBFS beside 60 Hz, on 50 Hz mains. Once decided it stays 50 and the tracker does not
// move under the note; a note there from the start never makes Auto say 60.
static void testB1()
{
    group ("Auto: 50 Hz hum + a B1 61.7 Hz note at -20 dBFS -> still 50, the tracker not pulled");
    const double fs = 48000.0;
    HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.seed = 4u;
    {
        std::vector<float> x ((std::size_t) (12.0 * fs), 0.0f);
        addHum (x, hs, fs);
        Note b1; b1.f0 = 61.735; b1.amp = 0.1; b1.start = 2.5; b1.dur = 8.0;
        addNote (x, b1, fs);
        Tel tel;
        processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        double fAt = 0.0, dev = 0.0; bool stay = true;
        for (std::size_t i = 0; i < tel.f.size(); ++i)
        {
            const double t = tOf (tel, i, fs);
            if (t < 2.5) fAt = tel.f[i];
            else { dev = std::max (dev, std::fabs (tel.f[i] - fAt)); stay = stay && tel.det[i] == 50 && tel.run[i] == 1; }
        }
        std::printf ("      decided before the note; through 8 s of B1 the estimate moved %.5f Hz\n", dev);
        ok (stay, "decided 50 before the note, still 50 through it and after");
        ok (dev <= 0.05, "the tracker is not pulled by B1 (moved " + fmt (dev, 5) + " Hz, bound 0.05)");
    }
    {
        // hum at -50 dBFS: B1's leakage into the 60/120 bins then beats the 50 family by more than 6 dB, so
        // only the quiet gate (broadband power within 12 dB of the bins) keeps Auto from saying 60
        HumSpec quietHum = hs; quietHum.a1 = 0.00316;
        std::vector<float> x ((std::size_t) (9.0 * fs), 0.0f);
        addHum (x, quietHum, fs);
        Note b1; b1.f0 = 61.735; b1.amp = 0.1; b1.start = 0.2; b1.dur = 6.0;
        addNote (x, b1, fs);
        Tel tel;
        processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        bool never60 = true;
        for (int d : tel.det) never60 = never60 && d != 60;
        ok (never60, "B1 sounding from the start never makes Auto say 60");
        ok (tel.det.back() == 50, "and once the note is gone it decides 50");
    }
}

// NEVER FLIPS WHILE PLAYING: decided 50, then 10 s of loud notes that put their energy in the 60/120 bins.
static void testNeverFlips()
{
    group ("Auto: never flips while playing (10 s of B1 / A#1 / B2 / A#2 over 50 Hz hum)");
    const double fs = 48000.0;
    HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.seed = 6u;
    std::vector<float> x ((std::size_t) (13.0 * fs), 0.0f);
    addHum (x, hs, fs);
    const double fr[] = { 61.735, 58.270, 123.471, 116.541 };
    for (int i = 0; i < 20; ++i)
    {
        Note n; n.f0 = fr[i % 4]; n.amp = 0.2; n.start = 2.5 + 0.5 * i; n.dur = 1.5; n.tau = 1.0; n.attack = 0.001;
        addNote (x, n, fs);
    }
    Tel tel;
    processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
    bool stay = true;
    for (std::size_t i = 0; i < tel.det.size(); ++i)
        if (tOf (tel, i, fs) >= 2.5) stay = stay && tel.det[i] == 50 && tel.run[i] == 1;
    ok (stay, "decided 50 before the playing, 50 throughout it");
}

//------------------------------------------------------------------------------------------------------
// FALSE LEARNING: hum locked for 3 s, then a note held 8 s beside a harmonic. The tracker must not move
// (bound printed); the output must equal the clean pass (rms error, and per partial for steady notes).
static void testFalseLearning()
{
    group ("false learning: held notes beside harmonics, vibrato, a bend, a chord vs the clean pass");
    const double fs = 48000.0;
    struct C { const char* name; double mains; std::vector<Note> notes; bool steady; };
    auto note = [] (double f0, double amp) { Note n; n.f0 = f0; n.amp = amp; n.start = 3.0; n.dur = 8.0; return n; };
    Note vib = note (659.255, 0.1); vib.vibRate = 5.0; vib.vibDepth = 0.004;
    Note bend = note (640.0, 0.1); bend.bendRatio = 665.0 / 640.0; bend.bendSec = 1.0;
    std::vector<Note> chord;
    for (double f : { 82.407, 123.471, 164.814, 207.652, 246.942, 329.628 }) chord.push_back (note (f, 0.05));
    const C cases[] = { { "F4 349.2 Hz held 8 s (50 Hz mains, k7 350)", 50.0, { note (349.228, 0.1) }, true },
                        { "E5 659.3 Hz held 8 s (60 Hz mains, k11 660)", 60.0, { note (659.255, 0.1) }, true },
                        { "A2 110 Hz, 5th partial 550.7 (50 Hz mains, k11 550)", 50.0, { note (110.0, 0.1) }, true },
                        { "E5 with 5 Hz vibrato +-0.4 % (60 Hz mains)", 60.0, { vib }, false },
                        { "a bend 640 -> 665 Hz through 660 (60 Hz mains)", 60.0, { bend }, false },
                        { "an E major chord, six notes (60 Hz mains)", 60.0, chord, false } };
    for (const C& c : cases)
    {
        std::vector<float> clean ((std::size_t) (12.0 * fs), 0.0f);
        for (const Note& n : c.notes) addNote (clean, n, fs);
        std::vector<float> x = clean;
        HumSpec hs; hs.f0 = c.mains; hs.a1 = 0.01; hs.seed = 12u;
        addHum (x, hs, fs);
        Tel tel;
        const auto y = processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        const auto yc = processed (clean, lti (c.mains > 55.0 ? dh::Mains::Hz60 : dh::Mains::Hz50, 0.5), fs);
        double fOn = 0.0, dev = 0.0;
        for (std::size_t i = 0; i < tel.f.size(); ++i)
        {
            const double t = tOf (tel, i, fs);
            if (t < 3.0) fOn = tel.f[i];
            else if (t < 11.0) dev = std::max (dev, std::fabs (tel.f[i] - fOn));
        }
        const double err = rmsDiffDb (y, yc, (std::int64_t) (3.0 * fs), (std::int64_t) (11.0 * fs));
        double worstP = 0.0;
        if (c.steady)
        {
            std::vector<double> fr;
            for (int p = 1; p <= 6; ++p) fr.push_back (c.notes[0].partialHz (p));
            const auto a = fitAmps (y, (std::int64_t) (7.0 * fs), (std::int64_t) (3.5 * fs), fr, fs);
            const auto b = fitAmps (yc, (std::int64_t) (7.0 * fs), (std::int64_t) (3.5 * fs), fr, fs);
            for (std::size_t p = 0; p < fr.size(); ++p) worstP = std::max (worstP, std::fabs (dbOf (a[p] / b[p])));
        }
        std::printf ("      %-52s tracker moved %.5f Hz; out - clean %.1f dB re clean%s\n", c.name, dev, err,
                     c.steady ? (", worst partial " + fmt (worstP, 3) + " dB").c_str() : "");
        ok (tel.det.back() == (int) c.mains && fOn > 0.0, std::string (c.name) + ": precondition — locked before the note");
        ok (dev <= 0.005, std::string (c.name) + ": the tracker held through the note (moved " + fmt (dev, 5) + " Hz, bound 0.005)");
        ok (err <= -40.0, std::string (c.name) + ": output - clean pass " + fmt (err, 1) + " dB re the clean pass");
        if (c.steady) ok (worstP <= 1.0, std::string (c.name) + ": every partial within 1 dB of the clean pass (" + fmt (worstP, 3) + ")");
    }
}

//------------------------------------------------------------------------------------------------------
// PLUCKS beside harmonics (G3 / B3 / E4 / E5, tau 1 s, an onset every 0.5 s) after lock: the output equals
// the clean pass to -40 dB, and in the 50 ms after each onset nothing beyond it reaches -60 dB re the pluck.
static void testPlucks()
{
    group ("plucks beside harmonics: output - clean <= -40 dB, no onset artefact > -60 dB re the pluck");
    const double fs = 48000.0;
    for (int pass = 0; pass < 4; ++pass)
    {
        const double mains = pass == 1 ? 50.0 : 60.0;
        const int partials = pass >= 2 ? 1 : 6;          // passes 2, 3: pure tones (harmonics), the smallest crest factor
        const double attack = pass == 3 ? 0.040 : 0.001; // pass 3: volume swells — the residual rises slowly
        std::vector<float> clean ((std::size_t) (10.0 * fs), 0.0f);
        const double fr[] = { 196.0, 246.942, 329.628, 659.255 };
        std::vector<double> onsets;
        for (int i = 0; i < 12; ++i)
        {
            Note n; n.f0 = fr[i % 4]; n.amp = 0.1; n.start = 3.0 + 0.5 * i; n.tau = 1.0; n.attack = attack; n.release = 0.0; n.partials = partials;
            n.dur = 10.0 - n.start;
            addNote (clean, n, fs);
            onsets.push_back (n.start);
        }
        // a -200 dBFS floor under the clean pass: its input must not be digital black before the first pluck,
        // or its correction fades in from black there and the reference itself moves (the stage, fed hum, never is)
        { Lcg r { 3u }; for (float& v : clean) v += 1e-10f * r.bipolar(); }
        std::vector<float> x = clean;
        HumSpec hs; hs.f0 = mains; hs.a1 = 0.01; hs.seed = 14u;
        addHum (x, hs, fs);
        const auto y  = processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs);
        const auto yc = processed (clean, lti (mains > 55.0 ? dh::Mains::Hz60 : dh::Mains::Hz50, 0.5), fs);
        const double err = rmsDiffDb (y, yc, (std::int64_t) (3.0 * fs), (std::int64_t) (9.5 * fs));
        double worstOnset = -400.0;
        for (double t0 : onsets)
        {
            const std::int64_t a = std::llround (t0 * fs), b = a + (std::int64_t) (0.05 * fs);
            double d = 0.0, pk = 0.0;
            for (std::int64_t i = a; i < b; ++i)
            {
                d  = std::max (d, std::fabs ((double) y[(std::size_t) i] - (double) yc[(std::size_t) i]));
                pk = std::max (pk, std::fabs ((double) clean[(std::size_t) i]));
            }
            worstOnset = std::max (worstOnset, dbOf (d / pk));
        }
        const std::string tag = std::to_string ((int) mains) + " Hz mains" + (partials == 1 ? (attack > 0.01 ? ", pure-tone swells" : ", pure tones") : "");
        std::printf ("      %s: out - clean %.1f dB re clean; worst onset window %.1f dB re the pluck peak\n", tag.c_str(), err, worstOnset);
        ok (err <= -40.0, tag + ": output - clean pass " + fmt (err, 1) + " dB");
        ok (worstOnset <= -60.0, tag + ": the 50 ms after every onset stays within -60 dB of the pluck (" + fmt (worstOnset, 1) + ")");
    }
}

//------------------------------------------------------------------------------------------------------
// STEPS at mode switches at arbitrary sample offsets: the correction (in - out) has no step. Measured as
// its second difference, which a step of size s shows as about s; the smooth -30 dBFS hum's own is ≈ -68 dBFS.
static double maxSecondDiff (const std::vector<float>& in, const std::vector<float>& out, std::int64_t from, std::int64_t to)
{
    double m = 0.0;
    for (std::int64_t i = std::max<std::int64_t> (from, 2); i < to; ++i)
    {
        auto c = [&] (std::int64_t j) { return (double) in[(std::size_t) j] - (double) out[(std::size_t) j]; };
        m = std::max (m, std::fabs (c (i) - 2.0 * c (i - 1) + c (i - 2)));
    }
    return m;
}
static void testSwitchSteps()
{
    group ("mode switches at arbitrary offsets: no step > -60 dBFS in the correction");
    const double fs = 48000.0;
    struct Ev { double at; dh::Mains to; };
    struct C { const char* name; double humF0, humF1, humSwitch; dh::Mains start; std::vector<Ev> ev; };
    const C cases[] = { { "Off -> 50", 50, 50, 99, dh::Mains::Off, { { 2.0 + 33.0 / fs, dh::Mains::Hz50 } } },
                        { "50 -> 60", 50, 50, 99, dh::Mains::Hz50, { { 3.0 + 17.0 / fs, dh::Mains::Hz60 } } },
                        { "60 -> Off", 60, 60, 99, dh::Mains::Hz60, { { 3.0 + 5.0 / fs, dh::Mains::Off } } },
                        { "Auto re-decision 50 -> 60 (the hum changes at 4.0003 s)", 50, 60, 4.0003, dh::Mains::Auto, {} } };
    for (const C& c : cases)
    {
        std::vector<float> x ((std::size_t) (8.0 * fs), 0.0f);
        HumSpec a; a.f0 = c.humF0; a.a1 = 0.0316; a.seed = 2u; a.until = c.humSwitch;
        addHum (x, a, fs);
        if (c.humSwitch < 90.0) { HumSpec b = a; b.f0 = c.humF1; b.until = 1e9; addHum (x, b, fs, std::llround (c.humSwitch * fs)); }
        HR h; h.setParams (params (c.start, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::vector<float> y = x;
        std::int64_t o = 0;
        for (const Ev& e : c.ev)
        {
            const std::int64_t at = std::llround (e.at * fs);
            float* io[1] { y.data() + o }; felitronics::test::run (h.process (io, 1, (int) (at - o)));
            h.setParams (params (e.to, 0.5, 3.0));
            o = at;
        }
        { float* io[1] { y.data() + o }; felitronics::test::run (h.process (io, 1, (int) ((std::int64_t) y.size() - o))); }
        const double m = maxSecondDiff (x, y, 0, (std::int64_t) y.size());
        std::printf ("      %-56s max |second difference of the correction| %.1f dBFS\n", c.name, dbOf (m));
        ok (dbOf (m) <= -60.0, std::string (c.name) + ": no step > -60 dBFS (" + fmt (dbOf (m), 1) + ")");
        if (c.ev.empty()) ok (h.detected() == dh::Mains::Hz60 && h.running() && std::fabs (h.mainsHz() - 60.0) < 0.05, "  and Auto re-decided 60");
    }
}

//------------------------------------------------------------------------------------------------------
// THE GHOST: a 330 Hz tone (20 Hz below k7 = 350 on 50 Hz mains) held 2 s over locked hum, then stopped.
// The notch's state rings out at 350 Hz. Bound: it starts at the LTI value A·(BW/2)/sqrt(Δ² + (BW/2)²)
// (−38.1 dB re the tone) and falls below -60 dB re the tone at the playing width's pace (2.5·τ = 1.6 s) —
// quiet does not return while the residual carries it. No step.
static void ghostRun (const dh::HumRemoverParams& p, double& startDb, double& t60, double& stepDb, double& tQuiet)
{
    const double fs = 48000.0, stop = 5.0;
    std::vector<float> x ((std::size_t) (9.0 * fs), 0.0f);
    HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.seed = 17u;
    addHum (x, hs, fs);
    Note n; n.f0 = 330.0; n.partials = 1; n.amp = 0.1; n.start = 3.0; n.dur = 2.0;
    addNote (x, n, fs);
    Tel tel;
    const auto y = processed (x, p, fs, &tel);
    tQuiet = -1.0;
    for (std::size_t i = 0; i < tel.quiet.size(); ++i)
        if (tOf (tel, i, fs) > stop && tel.quiet[i] && tQuiet < 0.0) tQuiet = tOf (tel, i, fs) - stop;
    const std::int64_t W = 1920;   // 40 ms, two periods
    startDb = -400.0; t60 = -1.0;
    for (std::int64_t a = std::llround (stop * fs); a + W < (std::int64_t) y.size(); a += 480)
    {
        const double g = dbOf (fitAmp (y, a, W, 350.0, fs) / n.amp);
        if (a == std::llround (stop * fs)) startDb = g;
        if (g > -60.0) t60 = (double) (a + W - std::llround (stop * fs)) / fs;
    }
    stepDb = dbOf (maxSecondDiff (x, y, std::llround (stop * fs), std::llround ((stop + 2.0) * fs)));
}
static void testGhost()
{
    group ("the ghost a note leaves: 330 Hz held 2 s, 20 Hz from k7, then stopped");
    const double theory = dbOf (0.25 / std::sqrt (20.0 * 20.0 + 0.25 * 0.25));
    double s1, t1, st1, q1, s2, t2, st2, q2;
    ghostRun (params (dh::Mains::Auto, 0.5, 3.0), s1, t1, st1, q1);
    ghostRun (params (dh::Mains::Auto, 0.5, 0.5), s2, t2, st2, q2);
    std::printf ("      LTI start %.1f dB re the tone. Learn width resuming (quiet again %.2f s after the stop): starts %.1f dB,"
                 " below -60 dB after %.2f s. Playing width only: starts %.1f dB, below -60 dB after %.2f s."
                 " Correction second difference over the 2 s after the stop: %.1f / %.1f dBFS\n",
                 theory, q1, s1, t1, s2, t2, st1, st2);
    ok (s1 <= theory + 1.0, "the ghost starts at the LTI value (" + fmt (s1, 1) + " dB, theory " + fmt (theory, 1) + ")");
    // NOTE SAFETY WINS: quiet now asks the residual to sit at its floor, and the ghost IS residual, so the
    // learn width does not come back while it rings — it decays at the playing width, like the LTI notch
    // (2.53·τ = 1.61 s to fall 22 dB). It used to fall below -60 dB in 0.70 s, when quiet admitted it.
    ok (t1 > 0.0 && t1 <= 2.0, "it falls below -60 dB re the note at the playing width's pace, within 2.0 s (" + fmt (t1, 2) + " s)");
    ok (std::fabs (t1 - t2) <= 0.2, "the same as with the playing width alone (" + fmt (t2, 2) + " s): the learn width does not chase a ghost");
    ok (st1 <= -60.0, "the switch back to the learn width steps nothing (" + fmt (st1, 1) + " dBFS)");
}

//------------------------------------------------------------------------------------------------------
// LEARN. 64-sample calls on the grid throughout, so a swap that ends on a boundary ends a call.
static void block64 (HR& h, std::vector<float>& x, std::int64_t& o, std::int64_t until)
{
    while (o < until)
    {
        const int n = (int) std::min<std::int64_t> (64, until - o);
        float* io[1] { x.data() + o };
        felitronics::test::run (h.process (io, 1, n));
        o += n;
    }
}

static void testLearnHumOnly()
{
    group ("LEARN on hum only: mains, hz within 0.02 Hz, floorAfterDb <= hum - 40 dB; cold and warm");
    const double fs = 48000.0;
    for (bool warm : { false, true })
    {
        HumSpec hs; hs.f0 = 59.8; hs.a1 = 0.0316; hs.seed = 33u;
        const double pre = warm ? 3.0 : 0.0, len = 2.5;
        std::vector<float> x ((std::size_t) ((pre + len + 2.0) * fs), 0.0f);
        addHum (x, hs, fs);
        const std::vector<float> in = x;
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::int64_t o = 0;
        block64 (h, x, o, (std::int64_t) (pre * fs));
        h.beginLearn();
        block64 (h, x, o, (std::int64_t) ((pre + len) * fs));
        const auto st = h.learnStatus();
        const auto r = h.endLearn (true);
        const std::string tag = warm ? "warm (3 s of hum first)" : "cold (learn from the first sample)";
        std::printf ("      %s: mains %d, %.4f Hz, hum %.1f dB, k1 %.1f dB, floor %.1f dBFS from %.0f ms, locked %d\n", tag.c_str(),
                     r.mains == dh::Mains::Hz60 ? 60 : (r.mains == dh::Mains::Hz50 ? 50 : 0), r.hz, r.humDb, r.harmonicDb[0],
                     r.floorAfterDb, r.settledAtMs, (int) h.locked());
        ok (st == HR::LearnStatus::Running && r.verdict == HR::LearnVerdict::Learned, tag + ": Learned");
        ok (r.mains == dh::Mains::Hz60 && std::fabs (r.hz - 59.8) < 0.02, tag + ": 60 Hz family at " + fmt (r.hz, 4) + " Hz");
        ok (r.floorAfterDb <= -30.0 - 40.0, tag + ": floorAfterDb " + fmt (r.floorAfterDb, 1) + " dBFS <= hum (-30) - 40");
        ok (std::fabs (r.harmonicDb[0] + 30.0) < 0.5 && std::fabs (r.harmonicDb[9] + 30.0 + 20.0) < 0.5,
            tag + ": the harmonic levels are the hum's (k1 " + fmt (r.harmonicDb[0], 2) + ", k10 " + fmt (r.harmonicDb[9], 2) + " dB)");
        if (warm)
        {
            const double st = dbOf (maxSecondDiff (in, x, (std::int64_t) (pre * fs) - 4800, (std::int64_t) (pre * fs) + 9600));
            ok (st <= -60.0, tag + ": beginLearn ramps the old model out — no step (" + fmt (st, 1) + " dBFS)");
        }
        block64 (h, x, o, (std::int64_t) x.size());
        const double nul = worstNull (x, in, (std::int64_t) ((pre + len + 1.0) * fs), (std::int64_t) fs - 2000, 59.8, fs);
        ok (h.running() && h.detected() == dh::Mains::Hz60 && nul <= -50.0 && h.locked(),
            tag + ": committed — it keeps cancelling (" + fmt (nul, 1) + " dB) and reports locked");
    }
}

static void testLearnDisturbedAndRestore()
{
    group ("LEARN disturbed by a pluck at 1.2 s: Disturbed, frozen in that period, nothing committed, restored bit for bit");
    const double fs = 48000.0;
    for (bool commit : { false, true })
    {
        HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.0316; hs.seed = 35u;
        std::vector<float> x ((std::size_t) (7.0 * fs), 0.0f);
        addHum (x, hs, fs);
        Note n; n.f0 = 329.628; n.amp = 0.1; n.start = 3.0 + 1.2 + 13.0 / fs; n.dur = 1.0; n.tau = 1.0; n.attack = 0.001;
        addNote (x, n, fs);
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::int64_t o = 0;
        block64 (h, x, o, (std::int64_t) (3.0 * fs));
        HR::Weight snap[25]; for (int k = 1; k <= 24; ++k) snap[k] = h.weight (0, k);
        const double fSnap = h.mainsHz(); const auto dSnap = h.detected();
        h.beginLearn();
        const std::int64_t onset = std::llround (n.start * fs), onsetBlockEnd = (onset / 64 + 1) * 64;
        block64 (h, x, o, onsetBlockEnd - 64);
        const bool runningBefore = h.learnStatus() == HR::LearnStatus::Running;
        block64 (h, x, o, onsetBlockEnd);
        const bool disturbedNow = h.learnStatus() == HR::LearnStatus::Disturbed;
        HR::Weight w1[25]; for (int k = 1; k <= 24; ++k) w1[k] = h.weight (0, k);
        block64 (h, x, o, (std::int64_t) (5.5 * fs));
        bool frozen = true;
        for (int k = 1; k <= 24; ++k) { const auto w = h.weight (0, k); frozen = frozen && std::memcmp (&w, &w1[k], sizeof w) == 0; }
        const auto r = h.endLearn (commit);
        int guard = 0;
        while (h.swapping() && guard++ < 100) block64 (h, x, o, o + 64);
        bool restored = ! h.swapping();
        for (int k = 1; k <= 24; ++k) { const auto w = h.weight (0, k); restored = restored && std::memcmp (&w, &snap[k], sizeof w) == 0; }
        const double fNow = h.mainsHz();
        const std::string tag = commit ? "endLearn(true)" : "endLearn(false)";
        ok (runningBefore && disturbedNow, tag + ": Running up to the pluck, Disturbed in the period that holds its onset");
        ok (frozen, tag + ": adaptation frozen from the trigger on (weights bit-identical to the end of the learn)");
        ok (r.verdict == HR::LearnVerdict::Disturbed, tag + ": the report says Disturbed");
        ok (restored && std::memcmp (&fNow, &fSnap, sizeof fNow) == 0 && h.detected() == dSnap,
            tag + ": nothing committed — weights, frequency and decision restored bit for bit (" + std::to_string (guard) + " calls of ramp)");
    }
}

static void testLearnNothingHeard()
{
    group ("LEARN hears nothing: digital silence and -95 dBFS hum");
    const double fs = 48000.0;
    for (int which = 0; which < 2; ++which)
    {
        std::vector<float> x ((std::size_t) (3.0 * fs), 0.0f);
        if (which == 1) { HumSpec hs; hs.f0 = 50.0; hs.a1 = 1.778e-5; hs.seed = 36u; addHum (x, hs, fs); }
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::int64_t o = 0;
        h.beginLearn();
        block64 (h, x, o, (std::int64_t) (1.5 * fs));
        const auto st = h.learnStatus();
        block64 (h, x, o, (std::int64_t) (2.5 * fs));
        const auto r = h.endLearn (true);
        const std::string tag = which == 0 ? "digital silence" : "-95 dBFS hum";
        std::printf ("      %s: mains %d, hum %.1f dB, verdict %d\n", tag.c_str(), r.mains == dh::Mains::Hz50 ? 50 : (r.mains == dh::Mains::Hz60 ? 60 : 0), r.humDb, (int) r.verdict);
        ok (st == HR::LearnStatus::NothingHeard && r.verdict == HR::LearnVerdict::NothingHeard, tag + ": NothingHeard, during and in the report");
        ok (h.detected() == dh::Mains::Off, tag + ": and nothing was committed");
    }
}

//------------------------------------------------------------------------------------------------------
// G2 — WHAT QUIET ADMITS. A note is not quiet, however stationary; a real hum change still is.
static void testQuietAdmitsNoNote()
{
    group ("quiet admits no note: C1 a soft note under the hum, C2 a bass note at a bin, C3 a 60 Hz tone over 50 Hz hum");
    const double fs = 48000.0;
    {   // C1: hum -30 dBFS, F4 349.2 Hz (0.77 Hz from k7) held 8 s at -40 dBFS — under the hum
        std::vector<float> clean ((std::size_t) (12.0 * fs), 0.0f);
        Note n; n.f0 = 349.228; n.amp = 0.01; n.start = 3.0; n.dur = 8.0; n.partials = 1;
        addNote (clean, n, fs);
        std::vector<float> x = clean; HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.0316; hs.seed = 51u; addHum (x, hs, fs);
        const auto y = processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs), yc = processed (clean, lti (dh::Mains::Hz50, 0.5), fs);
        const std::int64_t S = (std::int64_t) (7.0 * fs), W = (std::int64_t) (3.5 * fs);
        const double d = dbOf (fitAmp (y, S, W, n.f0, fs) / fitAmp (yc, S, W, n.f0, fs)), raw = dbOf (fitAmp (y, S, W, n.f0, fs) / fitAmp (clean, S, W, n.f0, fs));
        std::printf ("      C1: F4 at -40 dBFS under -30 dBFS hum: %.2f dB vs the clean pass, %.2f dB vs the untouched note\n", d, raw);
        ok (std::fabs (d) <= 0.5, "C1: the soft note is treated as the 0.5 Hz notch treats it (" + fmt (d, 2) + " dB vs the clean pass)");
    }
    {   // C2: G1 49.0 Hz held 8 s at -20 dBFS on 50 Hz mains (hum -40 dBFS)
        std::vector<float> clean ((std::size_t) (12.0 * fs), 0.0f);
        Note n; n.f0 = 48.999; n.amp = 0.1; n.start = 3.0; n.dur = 8.0;
        addNote (clean, n, fs);
        std::vector<float> x = clean; HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.seed = 52u; addHum (x, hs, fs);
        Tel tel;
        const auto y = processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel), yc = processed (clean, lti (dh::Mains::Hz50, 0.5), fs);
        double fOn = 0.0, dev = 0.0;
        for (std::size_t i = 0; i < tel.f.size(); ++i) { const double t = tOf (tel, i, fs); if (t < 3.0) fOn = tel.f[i]; else dev = std::max (dev, std::fabs (tel.f[i] - fOn)); }
        const std::int64_t S = (std::int64_t) (7.0 * fs), W = (std::int64_t) (3.5 * fs);
        const double d = dbOf (fitAmp (y, S, W, n.partialHz (1), fs) / fitAmp (yc, S, W, n.partialHz (1), fs));
        std::printf ("      C2: G1 49.0 Hz held 8 s: tracker moved %.4f Hz; fundamental %.2f dB vs the clean pass\n", dev, d);
        ok (dev <= 0.005, "C2: the tracker does not walk to the note (" + fmt (dev, 4) + " Hz)");
        ok (std::fabs (d) <= 0.5, "C2: the note is not learned at the fast width (" + fmt (d, 2) + " dB)");
    }
    {   // C3: locked 50 Hz hum, then an exact 60 + 120 Hz tone for 3 s
        const double fsx = fs;
        std::vector<float> x ((std::size_t) (8.0 * fsx), 0.0f);
        HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.01; hs.seed = 53u; addHum (x, hs, fsx);
        Note t60; t60.f0 = 60.0; t60.B = 0.0; t60.amp = 0.1; t60.partials = 2; t60.start = 3.0; t60.dur = 3.0;
        addNote (x, t60, fsx);
        Tel tel;
        processed (x, params (dh::Mains::Auto, 0.5, 3.0), fsx, &tel);
        bool stay = true;
        for (std::size_t i = 0; i < tel.det.size(); ++i) if (tOf (tel, i, fsx) >= 2.0) stay = stay && tel.det[i] == 50 && tel.run[i] == 1;
        ok (stay, "C3: a sustained exact 60/120 Hz tone over the 50 Hz hum never flips Auto");
    }
}

static void testHumChanges()
{
    group ("a real hum change still re-qualifies: +10 dB louder, a 50 -> 60 swap");
    const double fs = 48000.0;
    {
        std::vector<float> x ((std::size_t) (8.0 * fs), 0.0f);
        HumSpec hs; hs.f0 = 50.1; hs.a1 = 0.01; hs.seed = 54u; addHum (x, hs, fs);
        const std::int64_t c = std::llround (4.0003 * fs);
        for (std::size_t i = (std::size_t) c; i < x.size(); ++i) x[i] *= 3.1623f;
        const std::vector<float> in = x;
        Tel tel;
        const auto y = processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        double back = 0.0;   // the last moment, within 3 s of the change, at which it was not quiet
        for (std::size_t i = 0; i < tel.quiet.size(); ++i)
        {
            const double t = tOf (tel, i, fs);
            if (t > 4.0003 && t < 7.0 && ! tel.quiet[i]) back = t - 4.0003;
        }
        const double nul = worstNull (y, in, (std::int64_t) (6.5 * fs), (std::int64_t) fs, 50.1, fs);
        std::printf ("      +10 dB louder hum: quiet again %.2f s after the change; null 2.5 s after %.1f dB\n", back, nul);
        ok (back <= 1.0, "+10 dB: quiet re-qualifies within 1.0 s (" + fmt (back, 2) + " s)");
        ok (nul <= -50.0, "+10 dB: and the new hum is nulled (" + fmt (nul, 1) + " dB)");
    }
    {
        std::vector<float> x ((std::size_t) (8.0 * fs), 0.0f);
        HumSpec a; a.f0 = 50.0; a.a1 = 0.01; a.seed = 55u; a.until = 4.0003; addHum (x, a, fs);
        HumSpec b = a; b.f0 = 60.0; b.until = 1e9; addHum (x, b, fs, std::llround (4.0003 * fs));
        Tel tel;
        processed (x, params (dh::Mains::Auto, 0.5, 3.0), fs, &tel);
        double t60 = -1.0;
        for (std::size_t i = 0; i < tel.det.size(); ++i) if (tel.det[i] == 60 && t60 < 0.0) t60 = tOf (tel, i, fs) - 4.0003;
        std::printf ("      50 -> 60 hum swap: Auto re-decided %.2f s after it\n", t60);
        ok (t60 > 0.0 && t60 <= 1.5, "swap: Auto re-decides 60 within 1.5 s (" + fmt (t60, 2) + " s)");
    }
}

static void testLearnAfterChanges()
{
    group ("LEARN: C4 pressed while the hum fades back in, after the rig changed, and a slow swell into the floor window");
    const double fs = 48000.0;
    auto noise = [] (std::vector<float>& x, std::int64_t from, std::int64_t to, double rms, std::uint32_t seed)
    {
        Lcg r { seed };
        for (std::int64_t i = from; i < to; ++i) x[(std::size_t) i] += (float) (rms * std::sqrt (3.0)) * r.bipolar();
    };
    for (int sc = 0; sc < 3; ++sc)
    {
        const char* name[] = { "C4: pressed 0.2 s after the hum faded back in (floor still low)", "after the rig changed (noise +15 dB, hum +6 dB)",
                               "a slow swell (0.3 s attack, -50 dBFS) into the floor window" };
        std::vector<float> x ((std::size_t) (10.0 * fs), 0.0f);
        HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.0316; hs.seed = 56u; addHum (x, hs, fs);
        double learnAt = 3.0;
        if (sc == 0)
        {
            learnAt = 6.7;   // the fade-in itself is an onset; LEARN is pressed once it is done, the floor still old
            for (std::size_t i = 0; i < x.size(); ++i)
            {
                const double t = (double) i / fs;
                const double g = t < 3.0 ? 1.0 : (t < 6.0 ? 0.0 : std::min (1.0, (t - 6.0) / 0.5));
                x[i] *= (float) g;
            }
            noise (x, 0, (std::int64_t) x.size(), 3.16e-5, 7u);
        }
        else if (sc == 1)
        {
            noise (x, 0, (std::int64_t) (3.0 * fs), 3.16e-5, 8u);
            for (std::size_t i = (std::size_t) (3.0 * fs); i < x.size(); ++i) x[i] *= 2.0f;
            noise (x, (std::int64_t) (3.0 * fs), (std::int64_t) x.size(), 1.78e-4, 9u);
            learnAt = 3.3;
        }
        else
        {
            Note sw; sw.f0 = 329.628; sw.partials = 1; sw.amp = 0.00316; sw.attack = 0.3; sw.start = 3.0 + 0.9; sw.dur = 3.0;
            addNote (x, sw, fs);
        }
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::int64_t o = 0;
        block64 (h, x, o, (std::int64_t) (learnAt * fs));
        h.beginLearn();
        block64 (h, x, o, (std::int64_t) ((learnAt + 0.98) * fs));
        const auto stEarly = h.learnStatus();
        block64 (h, x, o, (std::int64_t) ((learnAt + 2.5) * fs));
        const auto r = h.endLearn (true);
        std::printf ("      %-58s status at 0.98 s %d, verdict %d, %.4f Hz, floor %.1f dBFS\n", name[sc], (int) stEarly, (int) r.verdict, r.hz, r.floorAfterDb);
        if (sc < 2) ok (r.verdict == HR::LearnVerdict::Learned && std::fabs (r.hz - 50.0) < 0.02, std::string (name[sc]) + ": Learned, not Disturbed");
        else
        {
            ok (stEarly == HR::LearnStatus::Running, std::string (name[sc]) + ": no per-sample trigger caught it (Running at 0.98 s)");
            ok (r.verdict == HR::LearnVerdict::Disturbed, std::string (name[sc]) + ": Disturbed by the floor-window rule, not committed");
        }
    }
}

//------------------------------------------------------------------------------------------------------
// G4 — THE LEARN CONTRACT. E1 the report says whether the model was kept; E2 a floor window of a moment is
// TooShort; E3 a learn begun while a cancelled one still ramps out keeps the ORIGINAL snapshot; E4 honest
// verdicts off a learn, and a pluck early in a cold learn is Disturbed; I the FLL is not diluted by a
// stopped channel's held weights.
static HR warmOn (std::vector<float>& x, std::int64_t& o, dh::Mains m, double until, double fs)
{
    HR h; h.setParams (params (m, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
    block64 (h, x, o, (std::int64_t) (until * fs));
    return h;
}

static void testLearnContract()
{
    group ("LEARN contract: kept, TooShort, begin during a restore, NotLearning, a pluck early in a cold learn");
    const double fs = 48000.0;
    auto humOf = [&] (double f0, double secs, std::uint32_t seed)
    {
        std::vector<float> x ((std::size_t) (secs * fs), 0.0f);
        HumSpec hs; hs.f0 = f0; hs.a1 = 0.0316; hs.seed = seed; addHum (x, hs, fs);
        return x;
    };
    {   // E1: kept only where the mode keeps the learned family
        struct C { dh::Mains m; bool kept; const char* name; };
        for (C c : { C { dh::Mains::Auto, true, "Auto" }, C { dh::Mains::Hz60, true, "60" }, C { dh::Mains::Hz50, false, "50 (hum is 60)" }, C { dh::Mains::Off, false, "Off" } })
        {
            auto x = humOf (60.0, 6.0, 71u); std::int64_t o = 0;
            HR h = warmOn (x, o, c.m, 2.0, fs);
            h.beginLearn(); block64 (h, x, o, (std::int64_t) (4.5 * fs));
            const auto r = h.endLearn (true);
            ok (r.verdict == HR::LearnVerdict::Learned && r.kept == c.kept, std::string ("E1: commit in mode ") + c.name + ": Learned, kept = " + (r.kept ? "true" : "false"));
        }
    }
    {   // E2: ended 10 ms into the floor window -> TooShort; 0.3 s into it -> Learned
        for (double extra : { 0.01, 0.30 })
        {
            auto x = humOf (50.0, 6.0, 72u); std::int64_t o = 0;
            HR h = warmOn (x, o, dh::Mains::Auto, 3.0, fs);
            h.beginLearn(); block64 (h, x, o, (std::int64_t) ((3.0 + 1.0 + extra) * fs));
            const auto r = h.endLearn (true);
            const bool want = extra > 0.1;
            ok (want ? r.verdict == HR::LearnVerdict::Learned : r.verdict == HR::LearnVerdict::TooShort,
                "E2: ended " + fmt (extra, 2) + " s into the floor window: " + (want ? "Learned" : "TooShort") + " (verdict " + std::to_string ((int) r.verdict) + ")");
        }
    }
    {   // E3: begin -> end(false) -> begin at once (the restore still ramping out) -> Disturbed -> end(false)
        auto x = humOf (50.0, 9.0, 73u); std::int64_t o = 0;
        Note n; n.f0 = 329.628; n.amp = 0.1; n.start = 3.5 + 1.5; n.dur = 1.0; n.tau = 1.0; n.attack = 0.001; addNote (x, n, fs);
        HR h = warmOn (x, o, dh::Mains::Auto, 3.0, fs);
        HR::Weight snap[25]; for (int k = 1; k <= 24; ++k) snap[k] = h.weight (0, k);
        const double fSnap = h.mainsHz();
        h.beginLearn(); block64 (h, x, o, (std::int64_t) (3.5 * fs));
        (void) h.endLearn (false);
        const bool swappingAtRebegin = h.swapping();
        h.beginLearn(); block64 (h, x, o, (std::int64_t) (6.0 * fs));
        const auto r = h.endLearn (false);
        int guard = 0; while (h.swapping() && guard++ < 100) block64 (h, x, o, o + 64);
        bool same = ! h.swapping() && r.verdict == HR::LearnVerdict::Disturbed;
        for (int k = 1; k <= 24; ++k) { const auto w = h.weight (0, k); same = same && std::memcmp (&w, &snap[k], sizeof w) == 0; }
        const double fNow = h.mainsHz();
        ok (swappingAtRebegin, "E3: precondition — the second LEARN begins while the first one's restore still ramps");
        ok (same && std::memcmp (&fNow, &fSnap, sizeof fNow) == 0, "E3: cancelling the second LEARN restores the ORIGINAL model, bit for bit");
    }
    {   // E4: honest verdicts off a learn
        auto x = humOf (50.0, 3.0, 74u); std::int64_t o = 0;
        HR h = warmOn (x, o, dh::Mains::Auto, 2.0, fs);
        ok (h.endLearn (true).verdict == HR::LearnVerdict::NotLearning, "E4: endLearn with no learn running says NotLearning");
        h.beginLearn(); block64 (h, x, o, (std::int64_t) (2.5 * fs));
        felitronics::test::run (h.prepare (fs, 512, 1));
        ok (! h.learning() && h.learnStatus() == HR::LearnStatus::Idle && h.endLearn (true).verdict == HR::LearnVerdict::NotLearning,
            "E4: prepare() during a learn cancels it (a restart): Idle, and endLearn says NotLearning");
    }
    {   // E4: a pluck 0.1 s into a COLD learn (no floor yet) is Disturbed
        auto x = humOf (50.0, 3.0, 75u);
        Note n; n.f0 = 329.628; n.amp = 0.1; n.start = 0.1; n.dur = 1.0; n.tau = 1.0; n.attack = 0.001; addNote (x, n, fs);
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 1));
        std::int64_t o = 0;
        h.beginLearn(); block64 (h, x, o, (std::int64_t) (2.5 * fs));
        const auto r = h.endLearn (true);
        ok (r.verdict == HR::LearnVerdict::Disturbed, "E4: a pluck 0.1 s into a cold learn makes it Disturbed (verdict " + std::to_string ((int) r.verdict) + ")");
    }
}

// I: the FLL with a stopped channel holding its (large) model tracks a mains ramp as a mono stage does.
static void testFllLiveChannels()
{
    group ("FLL: a stopped channel's held weights do not dilute the tracker (ramp tracked as in mono)");
    const double fs = 48000.0;
    HumSpec hs; hs.f0 = 50.0; hs.a1 = 0.0316; hs.slope = 0.0025; hs.rampFrom = 3.0; hs.rampTo = 23.0; hs.seed = 76u;
    std::vector<float> x ((std::size_t) (23.0 * fs), 0.0f);
    addHum (x, hs, fs);
    double worst[2] = { 0.0, 0.0 };
    for (int stereo = 0; stereo < 2; ++stereo)
    {
        HR h; h.setParams (params (dh::Mains::Auto, 0.5, 3.0)); felitronics::test::run (h.prepare (fs, 512, 2));
        std::vector<float> l = x, r = x;
        for (std::size_t o = 0; o < l.size(); o += 256)
        {
            const int n = (int) std::min<std::size_t> (256, l.size() - o);
            const bool both = stereo == 1 && (double) o / fs < 3.0;
            float* io[2] { l.data() + o, r.data() + o };
            felitronics::test::run (h.process (io, both ? 2 : 1, n));
            const double t = (double) (o + (std::size_t) n) / fs;
            if (t > 4.0) worst[stereo] = std::max (worst[stereo], std::fabs (h.mainsHz() - humFreq (hs, t)));
        }
    }
    std::printf ("      ramp lag: mono %.5f Hz, stereo with channel 1 stopped after 3 s %.5f Hz\n", worst[0], worst[1]);
    ok (worst[1] <= 1.1 * worst[0] + 1e-5, "the stopped channel does not slow the tracker (" + fmt (worst[1], 5) + " vs " + fmt (worst[0], 5) + " Hz)");
}

int main()
{
    std::printf ("felitronics::dehum — HumRemover knowing the mains\n");
    testDetuned();
    testRamp();
    testHiss();
    testAutoFamilies();
    testWrongFamily();
    testB1();
    testNeverFlips();
    testFalseLearning();
    testPlucks();
    testSwitchSteps();
    testGhost();
    testLearnHumOnly();
    testLearnDisturbedAndRestore();
    testLearnNothingHeard();
    testQuietAdmitsNoNote();
    testHumChanges();
    testLearnAfterChanges();
    testLearnContract();
    testFllLiveChannels();
    return felitronics::test::report();
}
