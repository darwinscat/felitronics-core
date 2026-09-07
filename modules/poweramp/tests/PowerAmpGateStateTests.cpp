// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// What a channel of this stage does when it stops being fed and is fed again later.
//
// Six recursions hide behind one `ch < nCh` here: the oversampler's polyphase histories, the DC blocker's
// x1/y1, the output transformer's two poles, and five per-channel SVF columns (presence, depth, mid, and
// the two virtual-load filters). None of them decays while the channel is absent. Measured before this
// suite: 0.649 (-3.8 dBFS) out of DIGITAL SILENCE after a stereo -> mono -> stereo excursion.
//
// `sag` is one supply shared by both channels and is deliberately left alone — it is SUPPOSED to follow
// whichever channels are actually present, so the second test asserts that it still does.

#include <felitronics_test.h>
#include <felitronics/poweramp/PowerAmpStage.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <vector>

static std::atomic<int> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

static constexpr double kFs = 48000.0;

static double peakOf (const std::vector<float>& v)
{
    double m = 0.0; for (float x : v) m = std::fmax (m, (double) std::fabs (x)); return m;
}

static bool bitEqual (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! (a[i] == b[i])) return false;
    return true;
}

static void fillNoise (std::vector<float>& v, unsigned seed)
{
    unsigned s = seed;
    for (auto& x : v) { s = s * 1664525u + 1013904223u; x = (float) ((int) (s >> 8) % 2000 - 1000) * 0.0004f; }
}

// Every per-channel filter turned ON, so the excursion has all six recursions to freeze. BOTH halves are
// needed: the knobs live in Params, but every filter they drive is gated by a VOICING depth that ships at
// zero (presenceMaxDb, depthMaxDb, midDb, loadResDb, loadRiseDb, otHfHz — each documented as "0 => inert").
// With a default Voicing the whole feel section is bit-for-bit absent no matter what the knobs say, which
// is exactly how a first version of this fixture measured a stage whose SVF columns never ran at all.
static poweramp::Params loudParams()
{
    poweramp::Params p;
    p.driveDb = 12.0f; p.presence = 0.8f; p.depth = 0.8f; p.load = 0.8f; p.iron = 0.6f; p.sag = 0.5f;
    return p;
}

static poweramp::Voicing liveVoicing()
{
    poweramp::Voicing v;
    v.presenceMaxDb = 6.0f; v.depthMaxDb = 5.0f; v.midDb = -4.0f;
    v.loadResDb = 5.0f; v.loadRiseDb = 4.0f; v.otHfHz = 6000.0f;
    v.sagMaxDroop = 0.25f; v.sagBiasDepth = 0.3f; v.nfbOpen = 0.4f;
    return v;
}

// Drive a stage to the end of a charge and report the last block's peak — used to prove the fixture is
// LIVE before anything is concluded from it.
static double lastPeakWith (const poweramp::Params& p, const poweramp::Voicing& v)
{
    const int N = 128;
    poweramp::PowerAmpStage a; a.prepare (kFs, N, 4);
    a.setParams (p, v);
    std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
    float* io[2] { L.data(), R.data() };
    double last = 0.0;
    for (int k = 0; k < 40; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            L[(std::size_t) i] = 0.0f;
            R[(std::size_t) i] = (float) (0.9 * std::sin (0.05 * (k * N + i)));
        }
        felitronics::test::run (a.process (io, 2, N));
        last = peakOf (R);
    }
    return last;
}

int main()
{
    std::printf ("felitronics::poweramp — execution-gate state\n");

    const int N = 128;

    group ("the fixture is live before anything is concluded from it");
    {
        poweramp::Params off = loudParams();
        off.presence = off.depth = off.load = off.iron = off.sag = 0.0f;
        const double withFeels    = lastPeakWith (loudParams(), liveVoicing());
        const double withoutFeels = lastPeakWith (off,          liveVoicing());
        ok (std::fabs (withFeels - withoutFeels) > 1.0e-4,
            "precondition: the feel filters actually change the output, so their state is real");
    }

    group ("channel gate: a returning channel emits nothing from digital silence");
    {
        poweramp::PowerAmpStage a; a.prepare (kFs, N, 4);
        a.setParams (loudParams(), liveVoicing());

        std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
        float* io[2] { L.data(), R.data() };
        double charged = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                L[(std::size_t) i] = 0.0f;
                R[(std::size_t) i] = (float) (0.9 * std::sin (0.05 * (k * N + i)));
            }
            felitronics::test::run (a.process (io, 2, N));
            charged = std::fmax (charged, peakOf (R));
        }
        ok (charged > 0.05, "precondition: the right channel really was driven");

        for (int k = 0; k < 20; ++k) { std::fill (L.begin(), L.end(), 0.0f); felitronics::test::run (a.process (io, 1, N)); }

        double worst = 0.0;
        for (int k = 0; k < 10; ++k)
        {
            std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f);
            felitronics::test::run (a.process (io, 2, N));
            worst = std::fmax (worst, std::fmax (peakOf (L), peakOf (R)));
        }
        ok (worst == 0.0, "silence in, exact zero out after 2 -> 1 -> 2 (was 0.649 = -3.8 dBFS)");
    }

    group ("the stage's OWN gates: a knob glided to zero stops a real recursion");
    {
        // Same defect, at a CONSTANT channel count. The knob glides (~25 ms at block rate), so the gate
        // closes about 84 blocks after the write — with the filter's state holding live signal. Everything
        // else in the chain is an FIR or a one-pole that a long silent stretch drains EXACTLY, which is
        // what makes "zero" the right oracle and is asserted as a precondition before the knob comes back.
        struct Gate { const char* name; float poweramp::Params::* knob; double wasDbfs; };
        const Gate gates[4] {
            { "presence", &poweramp::Params::presence, -44.9 },
            { "depth",    &poweramp::Params::depth,    -59.7 },
            { "load",     &poweramp::Params::load,     -35.5 },
            { "iron",     &poweramp::Params::iron,     -40.2 },
        };

        for (const Gate& g : gates)
        {
            poweramp::PowerAmpStage a; a.prepare (kFs, N, 4);
            poweramp::Params on = loudParams();
            on.presence = on.depth = on.load = on.iron = 0.0f;   // exactly ONE gate under test at a time
            on.*(g.knob) = 0.8f;
            a.setParams (on, liveVoicing());

            std::vector<float> L ((std::size_t) N, 0.0f), R ((std::size_t) N, 0.0f);
            float* io[2] { L.data(), R.data() };
            auto tone = [&] (int k) { for (int i = 0; i < N; ++i) L[(std::size_t) i] = R[(std::size_t) i] = (float) (0.8 * std::sin (0.07 * (k * N + i))); };
            auto hush = [&] { std::fill (L.begin(), L.end(), 0.0f); std::fill (R.begin(), R.end(), 0.0f); };

            double charged = 0.0;
            for (int k = 0; k < 40; ++k) { tone (k); felitronics::test::run (a.process (io, 2, N)); charged = std::fmax (charged, peakOf (L)); }
            ok (charged > 0.01, std::string ("precondition ") + g.name + ": the gate really was open and driven");

            poweramp::Params off = on; off.*(g.knob) = 0.0f;      // the knob goes to zero...
            a.setParams (off, liveVoicing());
            for (int k = 40; k < 200; ++k) { tone (k); felitronics::test::run (a.process (io, 2, N)); }   // ...and the gate closes mid-signal

            for (int k = 0; k < 1500; ++k) { hush(); felitronics::test::run (a.process (io, 2, N)); }     // drain everything that still runs
            ok (peakOf (L) == 0.0 && peakOf (R) == 0.0,
                std::string ("precondition ") + g.name + ": the rest of the chain drained to EXACT zero");

            a.setParams (on, liveVoicing());                                     // the knob comes back
            double worst = 0.0;
            for (int k = 0; k < 200; ++k) { hush(); felitronics::test::run (a.process (io, 2, N)); worst = std::fmax (worst, std::fmax (peakOf (L), peakOf (R))); }
            char msg[160];
            std::snprintf (msg, sizeof msg, "%s gate off -> on: exact zero out of silence (was %.1f dBFS)", g.name, g.wasDbfs);
            ok (worst == 0.0, msg);
        }
    }

    group ("isolation: the channel that never left keeps its history bit-exact");
    {
        poweramp::PowerAmpStage dut, ref;
        dut.prepare (kFs, N, 4); ref.prepare (kFs, N, 4);
        dut.setParams (loudParams(), liveVoicing());
        ref.setParams (loudParams(), liveVoicing());

        std::vector<float> d0 ((std::size_t) N), d1 ((std::size_t) N), r0 ((std::size_t) N), r1 ((std::size_t) N);
        float* dio[2] { d0.data(), d1.data() };
        float* rio[2] { r0.data(), r1.data() };

        // The supply (`sag`) is shared, so channel 0 legitimately moves when channel 1 leaves. Keep the
        // excursion channel SILENT: then the shared supply sees the same programme either way and the
        // only thing that can differ is the per-channel state this fix touches.
        bool equal = true; double energy = 0.0;
        for (int k = 0; k < 80; ++k)
        {
            fillNoise (d0, (unsigned) (11u + (unsigned) k));
            std::fill (d1.begin(), d1.end(), 0.0f);
            r0 = d0; r1 = d1;
            felitronics::test::run (dut.process (dio, (k >= 25 && k < 50) ? 1 : 2, N));
            felitronics::test::run (ref.process (rio, 2, N));
            equal = equal && bitEqual (d0, r0);
            energy += peakOf (d0);
        }
        ok (energy > 1.0, "precondition: the surviving channel carried signal");
        ok (equal, "channel 0 is bit-equal to the run where nothing ever left");
    }

    group ("RT-safety");
    {
        poweramp::PowerAmpStage a; a.prepare (kFs, N, 4);
        a.setParams (loudParams(), liveVoicing());
        std::vector<float> L ((std::size_t) N, 0.05f), R ((std::size_t) N, 0.05f);
        float* io[2] { L.data(), R.data() };
        felitronics::test::run (a.process (io, 2, N));
        const int before = g_allocs.load();
        for (int k = 0; k < 40; ++k) felitronics::test::run (a.process (io, (k % 2) + 1, N));
        felitronics::test::okNoAlloc (g_allocs.load() == before, "no allocation across 40 blocks of changing width");
    }

    return felitronics::test::report();
}
