// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::analysis::ReferenceTruePeakMeter — the reference true peak (P62).
//
// The decisive group is the first: the class IS `oversampling::PolyphaseOversampler` at 4x / 32 taps per phase,
// one per channel, maximum of |x| over the 4x stream, floored at the sample peak, drained with 32 zeros — the
// arithmetic `fcore::Probe` ran by hand before P62 — recomputed here by hand and compared BIT FOR BIT. A
// reference that is merely close to itself is not a reference. How far it reads from the OTHER meter is
// felitronics_truepeak_instrument_gap_tests' question; law 11 is felitronics_call_contract_tests'.

#include <felitronics_test.h>
#include <felitronics/analysis/ReferenceTruePeakMeter.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

// The byte counter counts what the CONTAINER asked for, which is the quantity a budget states. MSVC's STL on x86/x64
// asks operator new for sizeof(void*) + 31 more on a block of 4096 bytes or more (its hand alignment, <xmemory>) —
// the reference meter's scratch is such a block — so the counter takes that back off, exactly as
// LoudnessConformanceTests.cpp does and for the reason written there; the first law-11d check proves the correction
// is this STL's. Release only: iterator debugging adds proxies no counter can tell from storage.
#if defined(_MSVC_STL_VERSION) && (defined(_M_IX86) || defined(_M_X64))
#  if defined(_DEBUG)
static constexpr std::size_t kStlBigPad = 2 * sizeof (void*) + 31;
#  else
static constexpr std::size_t kStlBigPad = sizeof (void*) + 31;
#  endif
#else
static constexpr std::size_t kStlBigPad = 0;
#endif
static constexpr std::size_t kStlBigBlock = 4096;
static std::atomic<long>        g_allocs { 0 };
static std::atomic<std::size_t> g_bytes  { 0 };
static std::size_t containerBytes (std::size_t s) noexcept
{
    return kStlBigPad != 0 && s >= kStlBigBlock + kStlBigPad ? s - kStlBigPad : s;
}
void* operator new      (std::size_t s) { g_allocs.fetch_add (1); g_bytes.fetch_add (containerBytes (s)); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1); g_bytes.fetch_add (containerBytes (s)); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using RTP = analysis::ReferenceTruePeakMeter;

namespace
{
struct Lcg { std::uint32_t s; float next() { s = s * 1664525u + 1013904223u; return (float) (std::int32_t) s * 4.6566128730773926e-10f; } };

std::vector<std::vector<float>> programme (int nch, int n, std::uint32_t seed)
{
    std::vector<std::vector<float>> x ((std::size_t) nch, std::vector<float> ((std::size_t) n));
    Lcg r { seed };
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < n; ++i)
        {
            // Bright enough that the reconstruction exceeds the grid: noise plus a near-Nyquist tone.
            const double t = (double) i;
            x[(std::size_t) c][(std::size_t) i] = (float) (0.45 * r.next() + 0.4 * std::sin (2.0 * 3.141592653589793 * 0.37 * t + c));
        }
    return x;
}

// THE HAND-ROLLED PATH — the loop fcore::Probe carried before P62, written out again rather than shared.
struct HandRolled { double tp = 0.0, sp = 0.0; };
HandRolled handRolled (const std::vector<std::vector<float>>& x)
{
    HandRolled h;
    const int n = (int) x[0].size();
    for (const auto& ch : x)
    {
        oversampling::PolyphaseOversampler os;
        test::ok (os.prepare (4, 1, 32), "the hand-rolled oversampler prepares");
        std::vector<float> buf ((std::size_t) (n + 32) * 4u);
        const float* in[1] { ch.data() };
        float*       out[1] { buf.data() };
        os.upsample (in, 1, n, out);
        const float zeros[32] {};
        const float* zin[1] { zeros };
        float*       zout[1] { buf.data() + (std::size_t) n * 4u };
        os.upsample (zin, 1, 32, zout);
        for (float v : ch)  h.sp = std::max (h.sp, (double) std::fabs (v));
        for (float v : buf) h.tp = std::max (h.tp, (double) std::fabs (v));
    }
    h.tp = std::max (h.tp, h.sp);
    return h;
}

double feed (RTP& m, const std::vector<std::vector<float>>& x, int step)
{
    const int nch = (int) x.size(), n = (int) x[0].size();
    for (int off = 0; off < n; off += step)
    {
        const int k = std::min (step, n - off);
        const float* p[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) p[c] = x[(std::size_t) c].data() + off;
        test::run (m.process (p, nch, k));
    }
    m.drain();
    return m.truePeakLinear();
}
} // namespace

int main()
{
    std::printf ("felitronics::analysis ReferenceTruePeakMeter tests\n");
    const double fs = 48000.0;

    test::group ("the reference IS the 4x / 32-tap oversampler, floored and drained — bit for bit");
    {
        const auto x = programme (2, 30011, 7u);
        const HandRolled h = handRolled (x);
        RTP m; test::run (m.prepare (fs, 4096, 2));
        const double got = feed (m, x, 4096);
        test::ok (got == h.tp, "the reading equals the hand-rolled maximum exactly");
        test::ok (m.samplePeakLinear() == h.sp, "the sample peak equals the hand-rolled one exactly");
        test::ok (h.tp > h.sp * 1.001, "precondition — the programme's reconstruction exceeds its grid");
    }

    test::group ("chunking is BIT-invisible, across the internal walk's boundary and past it");
    {
        const auto x = programme (2, 100003, 99u);
        RTP ref; test::run (ref.prepare (fs, 100003, 2));
        const double whole = feed (ref, x, 100003);
        for (int step : { 1, 7, 1023, 1024, 1025, 4096, 65537 })
        {
            RTP m; test::run (m.prepare (fs, step, 2));
            test::ok (feed (m, x, step) == whole, "calls of " + std::to_string (step) + " frames read the whole-buffer number");
        }
    }

    test::group ("the rate does not shape the filter: the same samples read the same at every rate");
    {
        const auto x = programme (1, 20000, 3u);
        RTP a; test::run (a.prepare (44100.0, 20000, 1));
        const double base = feed (a, x, 20000);
        for (double r : { 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
        {
            RTP m; test::run (m.prepare (r, 20000, 1));
            test::ok (feed (m, x, 20000) == base, "identical at " + std::to_string ((int) r) + " Hz");
        }
    }

    test::group ("drain(): a peak in the last samples is measured, and a second drain adds nothing");
    {
        std::vector<std::vector<float>> x (1, std::vector<float> (4000, 0.0f));
        for (int i = 3990; i < 4000; ++i) x[0][(std::size_t) i] = 0.95f;
        RTP m; test::run (m.prepare (fs, 4000, 1));
        const float* p[1] { x[0].data() };
        test::run (m.process (p, 1, 4000));
        const double undrainedOs = m.truePeakLinearBlock();
        m.drain();
        const double drained = m.truePeakLinear();
        std::printf ("    block before the drain %.6f, drained %.6f\n", undrainedOs, drained);
        test::ok (drained > 1.0, "the ending's overshoot above full scale is seen only once drained");
        test::ok (undrainedOs < drained, "and was not visible before the drain");
        m.drain();
        test::ok (m.truePeakLinear() == drained, "one more drain moves nothing: the ring is already silence");
    }

    // The example TargetLoudnessSolver's measuring-rig comment quotes: a programme ending `..., 0, 1, 1`.
    test::group ("a programme ending on [0, 1, 1] reads +0 dBTP undrained and +1.833993 dBTP drained");
    {
        std::vector<std::vector<float>> x (1, std::vector<float> (64, 0.0f));
        x[0][62] = 1.0f; x[0][63] = 1.0f;
        RTP m; test::run (m.prepare (fs, 64, 1));
        const float* p[1] { x[0].data() };
        test::run (m.process (p, 1, 64));
        test::ok (m.truePeakDb() == 0.0, "undrained: the sample-peak floor, exactly 0 dBTP");
        m.drain();
        std::printf ("    drained %.6f dBTP\n", m.truePeakDb());
        test::approx (m.truePeakDb(), 1.833993, 5.0e-6, "drained: +1.833993 dBTP");
    }

    test::group ("the sample peak is a hard floor, and dB is the dB of the linear reading");
    {
        std::vector<std::vector<float>> x (1, std::vector<float> (512, 0.0f));
        x[0][300] = 0.5f;                                                  // a lone impulse: the FIR reads under it
        RTP m; test::run (m.prepare (fs, 512, 1));
        feed (m, x, 512);
        test::ok (m.truePeakLinear() >= 0.5 && m.samplePeakLinear() == 0.5, "the impulse reads at least its own sample");
        // `gainToDbDet`, which is what truePeakDb() now calls. Written as `core::gainToDb` this assertion
        // still PASSED — det and the system libm agree at this fixture's peak — so it was a sentence that
        // had become false about the code while staying green, which is worse than a red test.
        test::ok (m.truePeakDb() == core::gainToDbDet (m.truePeakLinear()), "truePeakDb() is gainToDbDet of the linear reading");
    }

    // THE FALLING EDGE DRAINS. The review round's sequence: the right channel ends on a peak still inside the FIR,
    // then one mono block. Dropping that history (what a delay line does under law 11a) read 0.95 where the
    // retained-and-drained reading — the hand-rolled path, which never forgets a channel — is 1.0625: an over lost.
    test::group ("a channel that stops is DRAINED: its peak still inside the filter is measured, and it returns silent");
    {
        std::vector<float> left (256, 0.0f), right (256, 0.0f), quiet (256, 0.0f);
        for (int i = 246; i < 256; ++i) right[(std::size_t) i] = 0.95f;   // an ending the grid does not show
        RTP m; test::run (m.prepare (fs, 256, 2));
        const float* both[2] { left.data(), right.data() };
        test::run (m.process (both, 2, 256));
        const float* mono[1] { quiet.data() };
        test::run (m.process (mono, 1, 256));                              // channel 1 stops here
        const double atStop = m.truePeakLinear();
        const HandRolled h = handRolled ({ left, right });                 // every channel kept, drained at the end
        std::printf ("    reading after the stop %.9f, hand-rolled %.9f\n", atStop, h.tp);
        test::ok (atStop == h.tp, "the stop reads exactly what keeping the channel and draining it reads");
        test::ok (h.tp > 1.0, "precondition — that reading is an over the grid does not show");
        m.drain();
        test::ok (m.truePeakLinear() == atStop, "and a drain afterwards adds nothing: the stopped channel is already silent");
        const float* silent[2] { quiet.data(), quiet.data() };
        test::run (m.process (silent, 2, 256));
        test::ok (m.truePeakLinearBlock() == 0.0, "stereo silence after the stop reads exactly zero: nothing replays");
    }

    test::group ("a zero-width call is a pause: every channel is drained, and the block reports that tail");
    {
        std::vector<float> x (256, 0.0f);
        for (int i = 246; i < 256; ++i) x[(std::size_t) i] = 0.95f;
        RTP m; test::run (m.prepare (fs, 256, 1));
        const float* p[1] { x.data() };
        test::run (m.process (p, 1, 256));
        test::run (m.process (nullptr, 0, 64));
        test::ok (m.truePeakLinearBlock() > 1.0 && m.truePeakLinear() == m.truePeakLinearBlock(),
                  "the pause's block reading is the drained tail, and it reached the running maximum");
    }

    test::group ("prepare() refuses what it cannot honour, and a refusal leaves the meter unusable");
    {
        RTP m;
        test::ok (m.prepare (fs, 64, 2), "a sane configuration prepares");
        for (double bad : { 0.0, -48000.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() })
            test::ok (! m.prepare (bad, 64, 2), "refuses the rate " + std::to_string (bad));
        test::ok (! m.prepare (fs, 64, 0) && ! m.prepare (fs, 64, core::kMaxChannels + 1), "refuses 0 and kMaxChannels + 1 channels");
        std::vector<float> s (64, 0.1f);
        const float* p[1] { s.data() };
        test::ok (! m.process (p, 1, 64) && ! m.isPrepared(), "after a refused prepare(), process() refuses too");
        test::ok (! m.prepare (fs, 0, 2) && ! m.prepare (fs, -1, 2), "refuses a block of no samples");
        test::ok (! RTP::storageFor (fs, 64, 0).ok && RTP::storageFor (fs, 64, 0).bytes() == 0 && RTP::storageFor (fs, 0, 2).bytes() == 0,
                  "and the budget for a refusal is zero");
    }

    test::group ("a null plane is refused before anything moves");
    {
        RTP m; test::run (m.prepare (fs, 64, 2));
        std::vector<float> s (64, 0.5f);
        const float* p[2] { s.data(), nullptr };
        test::ok (! m.process (p, 2, 64), "refused");
        test::ok (m.truePeakLinear() == 0.0 && m.samplePeakLinear() == 0.0, "and nothing was measured");
    }

    test::group ("law 11d: prepare() asks the heap for exactly storageFor(); process() and drain() ask for nothing");
    {
        const std::size_t before = g_bytes.load();
        {
            std::vector<float> v;
            v.assign (4096, 0.0f);                        // 16 384 B: a padded block on MSVC's STL
            volatile float* sink = v.data();              // observed, so the allocation cannot be elided
            sink[0] = 1.0f;
        }
        const std::size_t counted = g_bytes.load() - before;
        test::ok (counted == 4096u * sizeof (float), "the byte counter counts a big vector as its container asked (" + std::to_string (counted) + ")");
    }
    for (int nch : { 1, 2, 6, core::kMaxChannels })
    {
        const auto x = programme (nch, 5000, 11u);
        const std::size_t b0 = g_bytes.load();
        RTP m;
        test::run (m.prepare (96000.0, 5000, nch));
        const std::size_t asked = g_bytes.load() - b0;
        test::ok (asked == RTP::storageFor (96000.0, 5000, nch).bytes(),
                  std::to_string (nch) + " ch: " + std::to_string (asked) + " bytes allocated, budget " + std::to_string (RTP::storageFor (96000.0, 5000, nch).bytes()));
        const long a0 = g_allocs.load();
        const float* p[core::kMaxChannels] {};
        for (int c = 0; c < nch; ++c) p[c] = x[(std::size_t) c].data();
        const bool accepted = m.process (p, nch, 5000);
        m.drain();
        const long allocated = g_allocs.load() - a0;
        test::run (accepted);
        test::okNoAlloc (allocated == 0, std::to_string (nch) + " ch: process() and drain() allocate nothing");
    }

    return test::report();
}
