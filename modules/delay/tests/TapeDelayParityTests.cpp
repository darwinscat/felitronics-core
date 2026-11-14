// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// delay::TapeDelay against its FROZEN self (TapeDelayFrozen.h, the class as it stood before it had more than one
// head): at the default head layout the shipping class must be today's sound TO THE BIT. Not within a tolerance —
// the goldens allow 1e-6 and cannot prove that — but sample for sample, compared as bytes, so a signed zero or a
// last-ulp difference counts. Two ways in:
//   · the named scenarios: every golden case and the re-slicing case, rendered through both;
//   · random SESSIONS: a seeded script of prepares (refused ones too), knob setters (absurd values, re-sends and
//     one-ulp nudges included), process calls of every width the contract names (stereo, a narrower mono call,
//     the pause nch == 0, refused widths and extents), stand-by flips, resets and poison, sent to both classes
//     alike. After EVERY call the verdict, the audio and every public tap must agree bit for bit. The script is
//     counted as it runs, and the counts are gated, so a generator that never reaches a case fails here rather
//     than passing blind.
#include <felitronics_test.h>

#include "TapeDelayFrozen.h"
#include "delay_cases.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;
namespace dc  = felitronics::test::delaycases;
namespace frz = felitronics::delay::frozen;

namespace
{
template <class T>
bool sameBits (T a, T b)
{
    return std::memcmp (&a, &b, sizeof (T)) == 0;
}

bool sameBits (const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() && (a.empty() || std::memcmp (a.data(), b.data(), a.size() * sizeof (float)) == 0);
}

/** A 64-bit LCG read from its HIGH bits: the low bits of a power-of-two LCG cycle with short periods (the last
    one alternates), which is how a generator ends up never reaching a case. */
struct Rng
{
    std::uint64_t s;
    explicit Rng (std::uint64_t seed) : s (seed * 0x9E3779B97F4A7C15ull + 0xD1B54A32D192ED03ull) { u32(); u32(); }
    std::uint32_t u32() { s = s * 6364136223846793005ull + 1442695040888963407ull; return (std::uint32_t) (s >> 32); }
    int  below (int n) { return (int) (((std::uint64_t) u32() * (std::uint64_t) n) >> 32); }
    int  range (int lo, int hi) { return lo + below (hi - lo + 1); }
    bool chance (int percent) { return below (100) < percent; }
    float unit() { return (float) (u32() >> 8) * (1.0f / 16777216.0f); }   // [0, 1)
    float signedUnit() { return 2.0f * unit() - 1.0f; }
};

/** The frozen class driven through the render scenarios, as delay_cases.h's CoreAdapter drives the shipping one. */
struct FrozenAdapter
{
    frz::TapeDelay d;
    bool ok = true;

    void prepare (double sr, int block, int nch) { ok = d.prepare (sr, block, nch) && ok; }
    void reset() { d.reset(); }
    void apply (const dc::Knobs& k)
    {
        frz::TapeDelayParams p;
        p.timeMs   = delay::timeMsFor (k.sync, k.division, k.bpm, k.timeMs);
        p.repeats  = k.repeats;
        p.darkHz   = delay::darkHzForPercent (k.darkPct);
        p.offsetMs = k.offsetMs;
        p.mix      = k.mix;
        d.setParams (p);
    }
    void process (float* const* io, int nch, int n, bool feed) { ok = d.process (io, nch, n, feed) && ok; }
};

bool sameRender (const dc::Case& c)
{
    double (*const noSin) (double) = [] (double) { return 0.0; };
    double (*const noExp) (double) = [] (double) { return 0.0; };
    auto a = dc::makeInput (c, noSin, noExp);
    auto b = a;
    dc::CoreAdapter now;
    FrozenAdapter   then;
    dc::render (c, now, a);
    dc::render (c, then, b);
    bool same = now.ok && then.ok && a.size() == b.size();
    for (std::size_t ch = 0; same && ch < a.size(); ++ch)
        same = sameBits (a[ch], b[ch]);
    return same;
}

void namedScenarios()
{
    group ("the named scenarios: every golden case and the re-slicing case, byte for byte against the frozen class");
    for (const auto& c : dc::goldenCases())
        ok (sameRender (c), std::string (c.name) + ": the same bytes as the frozen class");

    dc::Knobs k; k.timeMs = 23.0f; k.repeats = 0.97f; k.darkPct = 60.0f; k.offsetMs = -11.0f; k.mix = 0.8f;
    bool same = true;
    for (int block : { 4096, 1, 64, 333, 517 })
        same = same && sameRender ({ "slice", 48000.0, block, 2, dc::kNoiseBurst, 0.02, 0.25, k, {} });
    ok (same, "the re-slicing case at blocks 4096, 1, 64, 333 and 517: the same bytes as the frozen class");
}

//==============================================================================
// THE SESSIONS
//==============================================================================

/** What a session's script reached — gated at the end, so a case the generator never produced fails loudly. */
struct Census
{
    long long calls = 0, refusedCalls = 0, stereo = 0, narrower = 0, monoPrepared = 0, pauses = 0, zeroLength = 0;
    long long longCalls = 0, poisonCalls = 0, subnormalInput = 0, feedFlips = 0, standByCalls = 0, resets = 0;
    long long prepares = 0, refusedPrepares = 0, setParams = 0, singleSetters = 0, absurdKnobs = 0, reSends = 0;
    long long ulpNudges = 0, samples = 0, mismatches = 0;
};

/** The last value the script sent to each knob, so a re-send and a one-ulp nudge have something to repeat. */
struct Sent
{
    float timeMs = 350.0f, repeats = 0.35f, darkHz = 4000.0f, offsetMs = 0.0f, mix = 0.25f;
};

float knobValue (Rng& r, float lo, float hi, bool logScale, float last, Census& n)
{
    const float big = std::numeric_limits<float>::max();
    switch (r.below (24))
    {
        case 0:  ++n.absurdKnobs; return std::numeric_limits<float>::quiet_NaN();
        case 1:  ++n.absurdKnobs; return r.chance (50) ? std::numeric_limits<float>::infinity()
                                                       : -std::numeric_limits<float>::infinity();
        case 2:  ++n.absurdKnobs; return r.chance (50) ? big : -big;
        case 3:  ++n.absurdKnobs; return r.chance (50) ? 1.0e-30f : -1.0e-30f;
        case 4:  ++n.absurdKnobs; return r.chance (50) ? 0.0f : -0.0f;
        case 5:  ++n.absurdKnobs; return std::numeric_limits<float>::denorm_min();
        case 6:
        case 7:  ++n.reSends; return last;
        case 8:
        case 9:  ++n.ulpNudges; return std::nextafter (last, r.chance (50) ? big : -big);
        case 10: return lo;
        case 11: return hi;
        default:
            if (logScale)
                return lo * (float) std::pow ((double) hi / (double) lo, (double) r.unit());
            return lo + (hi - lo) * r.unit();
    }
}

/** One input plane's worth of signal, of a kind drawn per call. */
void fillInput (Rng& r, std::vector<float>& x, int kind, Census& n)
{
    const float amp = r.chance (20) ? 1.0f : 0.05f + 0.9f * r.unit();
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        switch (kind)
        {
            case 0:  x[i] = 0.0f; break;                                                        // silence
            case 1:  x[i] = amp * r.signedUnit(); break;                                         // noise
            case 2:  x[i] = amp * (float) std::sin (0.0491 * (double) i + 0.3); break;           // a tone
            case 3:  x[i] = r.chance (2) ? amp * r.signedUnit() : 0.0f; break;                   // sparse impulses
            case 4:  x[i] = (float) (r.chance (50) ? 1 : -1) * std::numeric_limits<float>::denorm_min()
                           * (float) r.range (1, 1000); break;                                   // subnormal input
            default: x[i] = amp * r.signedUnit(); break;
        }
    }
    if (kind == 5 && ! x.empty())   // poison: a NaN or an infinity somewhere in the call
    {
        x[(std::size_t) r.below ((int) x.size())] = std::numeric_limits<float>::quiet_NaN();
        x[(std::size_t) r.below ((int) x.size())] = r.chance (50) ? std::numeric_limits<float>::infinity()
                                                                  : -std::numeric_limits<float>::infinity();
        ++n.poisonCalls;
    }
    if (kind == 4)
        ++n.subnormalInput;
}

/** Every public tap of the two, bit for bit. The subnormal witness walks both whole lines, so it is asked only
    when `deep` (every 16th call and at the end of a session) — asked after every call it was the suite's cost. */
template <class Now>
bool sameTaps (const Now& a, const frz::TapeDelay& b, bool deep)
{
    return a.isPrepared() == b.isPrepared() && a.maxChannels() == b.maxChannels() && a.maxBlock() == b.maxBlock()
        && sameBits (a.tailSeconds(), b.tailSeconds()) && sameBits (a.shownTimeMs(), b.shownTimeMs())
        && sameBits (a.shownOffsetMs(), b.shownOffsetMs()) && sameBits (a.inputPeak(), b.inputPeak())
        && sameBits (a.headSamples(), b.headSamples()) && sameBits (a.offsetHeadSamples (0), b.offsetHeadSamples (0))
        && sameBits (a.offsetHeadSamples (1), b.offsetHeadSamples (1))
        && (! deep || a.anySubnormalState() == b.anySubnormalState());
}

/** The same knob through the same setter on both. */
template <class Now>
void setOne (Now& a, frz::TapeDelay& b, int which, float v, Sent& sent)
{
    switch (which)
    {
        case 0:  a.setTimeMs (v);   b.setTimeMs (v);   sent.timeMs = v;   break;
        case 1:  a.setRepeats (v);  b.setRepeats (v);  sent.repeats = v;  break;
        case 2:  a.setDarkHz (v);   b.setDarkHz (v);   sent.darkHz = v;   break;
        case 3:  a.setOffsetMs (v); b.setOffsetMs (v); sent.offsetMs = v; break;
        default: a.setMix (v);      b.setMix (v);      sent.mix = v;      break;
    }
}

float drawKnob (Rng& r, int which, const Sent& s, Census& n)
{
    switch (which)
    {
        case 0:  return knobValue (r, 0.5f, 6200.0f, true, s.timeMs, n);
        case 1:  return knobValue (r, -0.1f, 1.1f, false, s.repeats, n);
        case 2:  return knobValue (r, 0.5f, 120000.0f, true, s.darkHz, n);
        case 3:  return knobValue (r, -35.0f, 35.0f, false, s.offsetMs, n);
        default: return knobValue (r, -0.1f, 1.1f, false, s.mix, n);
    }
}

/** One session: `ops` calls drawn from `seed`, both classes driven alike, compared after every call. Returns the
    index of the first call after which they differ, or -1. */
template <class Now>
int session (std::uint64_t seed, int ops, Census& n)
{
    Rng r { seed };
    Now a;
    frz::TapeDelay b;
    Sent sent;

    static constexpr double kRates[] = { 1000.0, 8000.0, 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 192000.0 };
    int maxCh = 0;
    bool feed = true;
    std::vector<float> ia[2], ib[2], keep[2];

    if (r.chance (50))   // knobs before the first prepare, as a plugin's constructor hands its state over
    {
        const int which = r.below (5);
        setOne (a, b, which, drawKnob (r, which, sent, n), sent);
    }

    for (int op = 0; op < ops; ++op)
    {
        // A disarmed instance is prepared again soon: a script that sat refused for long stretches would spend
        // its calls proving one rule over and over.
        const int pick = op == 0 || (maxCh == 0 && r.chance (40)) ? 0 : r.below (100);
        bool verdictsAgree = true, audioAgrees = true;

        if (pick < 3)   // prepare — now and then a refused one, which disarms
        {
            double sr = kRates[r.below ((int) std::size (kRates))];
            if (r.chance (4)) sr = 768000.0;
            int block = std::array<int, 4> { 1, 64, 512, 4096 }[(std::size_t) r.below (4)];
            int ch = r.chance (30) ? 1 : 2;
            const bool refuse = op > 0 && r.chance (8);
            if (refuse)
            {
                switch (r.below (4))
                {
                    case 0:  sr = 999.0; break;
                    case 1:  sr = std::numeric_limits<double>::quiet_NaN(); break;
                    case 2:  ch = 3; break;
                    default: block = 0; break;
                }
            }
            const bool va = a.prepare (sr, block, ch), vb = b.prepare (sr, block, ch);
            verdictsAgree = va == vb && va == ! refuse;
            maxCh = va ? ch : 0;
            ++n.prepares;
            if (refuse) ++n.refusedPrepares;
            if (va && ch == 1) ++n.monoPrepared;
        }
        else if (pick < 13)   // setParams, every field drawn
        {
            frz::TapeDelayParams pb;
            pb.timeMs   = drawKnob (r, 0, sent, n);
            pb.repeats  = drawKnob (r, 1, sent, n);
            pb.darkHz   = drawKnob (r, 2, sent, n);
            pb.offsetMs = drawKnob (r, 3, sent, n);
            pb.mix      = drawKnob (r, 4, sent, n);
            delay::TapeDelayParams pa;
            pa.timeMs = pb.timeMs; pa.repeats = pb.repeats; pa.darkHz = pb.darkHz; pa.offsetMs = pb.offsetMs; pa.mix = pb.mix;
            a.setParams (pa);
            b.setParams (pb);
            sent = { pb.timeMs, pb.repeats, pb.darkHz, pb.offsetMs, pb.mix };
            ++n.setParams;
        }
        else if (pick < 25)   // one setter
        {
            const int which = r.below (5);
            setOne (a, b, which, drawKnob (r, which, sent, n), sent);
            ++n.singleSetters;
        }
        else if (pick < 30)   // reset — twice now and then: the second must be the cheap no-op
        {
            a.reset(); b.reset();
            if (r.chance (30)) { a.reset(); b.reset(); }
            ++n.resets;
        }
        else   // process
        {
            if (r.chance (12)) { feed = ! feed; ++n.feedFlips; }

            int nch;
            const int w = r.below (100);
            if (w < 12)       nch = 0;
            else if (w < 40 || (maxCh == 1 && w < 94)) nch = 1;   // a mono preparation mostly hears mono
            else if (w < 96)  nch = 2;
            else              nch = r.chance (50) ? -1 : 3;   // refused: a negative width, or wider than any prepare

            int len;
            const int e = r.below (100);
            if (e < 3)        len = 0;
            else if (e < 8)   len = 1;
            else if (e < 35)  len = r.range (2, 64);
            else if (e < 82)  len = r.range (65, 1024);
            else if (e < 98)  len = r.range (1025, 4096);
            else              len = -r.range (1, 64);         // refused: a negative extent

            const int planes = nch > 0 ? std::min (nch, 2) : 0;
            const int kind = r.below (6);
            for (int ch = 0; ch < 2; ++ch)
            {
                ia[ch].assign ((std::size_t) std::max (len, 0), 0.0f);
                if (ch < planes) fillInput (r, ia[ch], kind, n);
                ib[ch] = ia[ch];
                keep[ch] = ia[ch];
            }
            const bool nullPlanes = nch > 0 && r.chance (3);   // refused: planes that are not there
            float* pa[2] = { ia[0].data(), ia[1].data() };
            float* pb[2] = { ib[0].data(), ib[1].data() };
            const bool va = a.process (nullPlanes || nch == 0 ? nullptr : pa, nch, len, feed);
            const bool vb = b.process (nullPlanes || nch == 0 ? nullptr : pb, nch, len, feed);
            verdictsAgree = va == vb;
            for (int ch = 0; ch < 2; ++ch)
                audioAgrees = audioAgrees && sameBits (ia[ch], ib[ch]) && (va || sameBits (ia[ch], keep[ch]));

            ++n.calls;
            if (! va) ++n.refusedCalls;
            else if (len > 0)
            {
                n.samples += len;
                if (nch == 0) ++n.pauses;
                else if (nch < maxCh) ++n.narrower;
                else if (nch == 2) ++n.stereo;
                if (len > 1024) ++n.longCalls;
                if (! feed) ++n.standByCalls;
            }
            else if (len == 0) ++n.zeroLength;
        }

        if (! (verdictsAgree && audioAgrees && sameTaps (a, b, op % 16 == 15 || op + 1 == ops)))
        {
            ++n.mismatches;
            return op;
        }
    }
    return -1;
}

template <class Now>
void sessions (int count, int ops)
{
    group ("random sessions: every call of every kind, verdicts, audio and taps byte for byte against the frozen class");
    Census n;
    int failed = 0;
    for (int s = 0; s < count; ++s)
    {
        const int at = session<Now> ((std::uint64_t) s + 1u, ops, n);
        if (at >= 0)
        {
            if (failed < 5)
                std::printf ("      [session %d: the classes part after call %d]\n", s + 1, at);
            ++failed;
        }
    }
    ok (failed == 0, std::to_string (count) + " sessions of " + std::to_string (ops) + " calls: "
                     + std::to_string (failed) + " parted from the frozen class");

    std::printf ("      [census: %lld calls (%lld refused), %lld stereo, %lld narrower, %lld pauses, %lld zero-length, "
                 "%lld long, %lld poison, %lld subnormal-input, %lld stand-by, %lld feed flips, %lld resets, "
                 "%lld prepares (%lld refused, %lld mono), %lld setParams, %lld single setters, %lld absurd knobs, "
                 "%lld re-sends, %lld one-ulp nudges, %lld samples]\n",
                 n.calls, n.refusedCalls, n.stereo, n.narrower, n.pauses, n.zeroLength, n.longCalls, n.poisonCalls,
                 n.subnormalInput, n.standByCalls, n.feedFlips, n.resets, n.prepares, n.refusedPrepares, n.monoPrepared,
                 n.setParams, n.singleSetters, n.absurdKnobs, n.reSends, n.ulpNudges, n.samples);

    // The script reached every case it names — a floor on each, not a hope.
    ok (n.stereo >= 1000 && n.narrower >= 200 && n.pauses >= 200 && n.zeroLength >= 50 && n.longCalls >= 200,
        "the census: stereo, narrower mono calls, pauses, zero-length and long calls all reached");
    ok (n.refusedCalls >= 100 && n.refusedPrepares >= 10 && n.monoPrepared >= 20 && n.prepares >= 100,
        "the census: refused calls, refused prepares and mono preparations all reached");
    ok (n.poisonCalls >= 50 && n.subnormalInput >= 50 && n.standByCalls >= 200 && n.feedFlips >= 100 && n.resets >= 100,
        "the census: poison, subnormal input, stand-by, feed flips and resets all reached");
    ok (n.setParams >= 200 && n.singleSetters >= 200 && n.absurdKnobs >= 200 && n.reSends >= 100 && n.ulpNudges >= 100,
        "the census: whole and single knob sets, absurd values, re-sends and one-ulp nudges all reached");
}
} // namespace

int main()
{
    std::printf ("felitronics::delay::TapeDelay against its frozen self\n");
    namedScenarios();
    sessions<delay::TapeDelay> (60, 400);
    return felitronics::test::report();
}
