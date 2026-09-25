// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

//==================================================================================================
// LAW 11 (DSP-ARCHITECTURE.md §2) — the property harness of CallContractTests: the fixture helpers,
// the ADAPT macro that describes a stage, and the seven properties P1–P7 swept by allProperties<A>().
// Moved out of modules/mastering/tests/CallContractTests.cpp verbatim, so felitronics-guitar-core
// sweeps its stages with the SAME properties instead of a copy that could drift. The census of
// stages and the named cases stay with each repository's own suite.
//
// Include in ONE translation unit per executable, then `using namespace felitronics::test::law11;`.
// ADAPT expands `core::kMaxChannels` at the use site, so that TU also says
// `namespace core = felitronics::core;`.
//==================================================================================================
#pragma once

#include <felitronics_test.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace felitronics::test::law11 {

namespace core = felitronics::core;
using felitronics::test::ok;
using felitronics::test::group;

constexpr double kFs       = 48000.0;
constexpr int    kMaxBlock = 64;          // deliberately SMALL, so "past maxBlock" is easy to reach
constexpr int    kPrepCh   = 2;
constexpr double kPi       = 3.14159265358979323846;   // MSVC has no M_PI without _USE_MATH_DEFINES

// A deterministic, non-degenerate programme: every stage below must actually DO something to it, or
// the properties would hold on a stage that had stopped working.
struct Rng { std::uint32_t s; float next() { s = s * 1664525u + 1013904223u; return (float) (std::int32_t) s * 4.6566128730773926e-10f; } };

inline void fill (std::vector<std::vector<float>>& b, int nch, int n, std::uint32_t seed)
{
    Rng r { seed };
    b.assign ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f));
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < n; ++i)
            b[(std::size_t) c][(std::size_t) i] = 0.5f * std::sin (2.0 * kPi * (220.0 + 30.0 * c) * i / kFs) + 0.1f * r.next();
}

inline std::vector<float*> planes (std::vector<std::vector<float>>& b)
{
    std::vector<float*> p; p.reserve (b.size());
    for (auto& v : b) p.push_back (v.data());
    return p;
}

inline bool bitEqual (const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        if (a[c].size() != b[c].size()) return false;
        for (std::size_t i = 0; i < a[c].size(); ++i) if (a[c][i] != b[c][i]) return false;
    }
    return true;
}

// `OBS` is the stage's OWN observable: the number a chunk-equivalence comparison should use. For an
// in-place stage the buffer IS the observable and OBS is unused; for a METER the buffer is untouched by
// construction, so comparing it would be comparing two copies of the input — the blind fixture this
// project keeps finding. A meter therefore names its reading here instead.
#define ADAPT(NAME, TYPE, MAXCH, NARROW, CHUNKS, INPLACE, PREP, RUN, OBS, PREPW)                          \
    struct NAME {                                                                                  \
        using T = TYPE;                                                                            \
        static const char* name() { return #TYPE; }                                                \
        static constexpr int  maxCh    = (MAXCH);                                                  \
        static constexpr bool narrowOk = (NARROW);                                                 \
        static constexpr bool chunks   = (CHUNKS);                                                 \
        static constexpr bool inPlace  = (INPLACE);                                                \
        static constexpr bool needsPrepare = true;                                                 \
        static constexpr bool gapVisible   = true;                                                 \
        static constexpr double silenceTol = 1.0e-6;   /* what SILENT means for THIS stage */     \
        static bool prepWidth (T& s, int w) { PREPW }                                              \
        static constexpr bool hasPrepareWidth = true;                                              \
        static constexpr int  prepareCeiling  = core::kMaxChannels;  /* the MODULE's own limit */  \
        static double blockObservable (T& s) { (void) s; return 0.0; }                             \
        static bool prep (T& s) { PREP }                                                           \
        static bool run (T& s, float* const* io, int nch, int n) { RUN }                           \
        static double observable (T& s) { OBS }                                                    \
    }

//==================================================================================================
// P1 — the REJECTED domain, swept whole. Every width above the prepared one and every negative
// extent must be refused, must SAY it was refused, and must leave the caller's buffer untouched.
//==================================================================================================
template <class A>
void rejectedDomain()
{
    typename A::T s;
    ok (A::prep (s), std::string (A::name()) + ": the fixture prepares (precondition)");

    int cases = 0, refused = 0, untouched = 0;
    // widths above the prepared one, right up to two past the house maximum...
    for (int nch = A::maxCh + 1; nch <= core::kMaxChannels + 2; ++nch)
        for (int n : { 0, 1, 7, kMaxBlock, kMaxBlock + 1, 3 * kMaxBlock + 5 })
        {
            std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), std::max (n, 1), 0xB0B0u + (std::uint32_t) (nch * 97 + n));
            const auto before = b;
            auto p = planes (b);
            ++cases;
            if (! A::run (s, p.data(), nch, n)) ++refused;
            if (bitEqual (b, before)) ++untouched;
        }
    // ...and, where a narrow call is MEANINGLESS (law 11c), every width below the prepared one. Without
    // this the 11(c) clause had no coverage at all: a mutation that answered `true` to a narrow call and
    // processed nothing survived both this suite and MatrixConvolver's own.
    if constexpr (! A::narrowOk)
        for (int nch = 0; nch < A::maxCh; ++nch)
            for (int n : { 1, 7, kMaxBlock, kMaxBlock + 1 })
            {
                std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), n, 0x7A11u + (std::uint32_t) (nch * 53 + n));
                const auto before = b;
                auto p = planes (b);
                ++cases;
                if (! A::run (s, p.data(), nch, n)) ++refused;
                if (bitEqual (b, before)) ++untouched;
            }

    // ...and every malformed extent.
    for (int nch : { -1, -3, 0, 1, A::maxCh })
        for (int n : { -1, -2, -1000, 0, 64, 517 })   // a negative width with a VALID length too: sweeping
                                                      // it only alongside a negative length let the width
                                                      // check be deleted without the suite noticing
        {
            if (nch >= 0 && n >= 0) continue;                  // that combination is the ACCEPTED domain
            std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), std::max (n, 128), 0xC0C0u + (std::uint32_t) (nch * 13 - n));
            const auto before = b;
            auto p = planes (b);
            ++cases;
            if (! A::run (s, p.data(), nch, n)) ++refused;
            if (bitEqual (b, before)) ++untouched;
        }
    ok (cases > 0, std::string (A::name()) + ": the rejected sweep is non-empty (precondition)");
    ok (refused == cases,  std::string (A::name()) + ": all " + std::to_string (cases) + " rejected inputs are REFUSED and say so");
    ok (untouched == cases, std::string (A::name()) + ": ...and every one left the buffer bit-identical");
    if constexpr (! A::inPlace)
    {
        // A METER never writes the buffer, so "the buffer is untouched" is true of a stage that reset
        // its whole state on every refusal. Its OBSERVABLE is the reading, and that is what must not
        // move — a mutation that called reset() and then returned false survived without this.
        const double before = A::observable (s);
        std::vector<std::vector<float>> b; fill (b, A::maxCh, 128, 0xBEEFu);
        auto p = planes (b);
        ok (A::run (s, p.data(), A::maxCh, 128), "precondition: a valid call is accepted");
        const double moved = A::observable (s);
        std::vector<std::vector<float>> junk; fill (junk, A::maxCh + 1, 128, 0xF00Fu);
        auto jp = planes (junk);
        ok (! A::run (s, jp.data(), A::maxCh + 1, 128), "precondition: the injected call IS refused");
        ok (A::observable (s) == moved, std::string (A::name()) + ": a refused call did not move the READING either");
        (void) before;
    }
}

//==================================================================================================
// P2 — the ACCEPTED domain, swept whole. Includes lengths far past maxBlock (law 11a) and the
// degenerate widths and lengths (law 11d).
//==================================================================================================
template <class A>
void acceptedDomain()
{
    typename A::T s;
    ok (A::prep (s), std::string (A::name()) + ": the fixture prepares (precondition)");

    int cases = 0, accepted = 0;
    const int lo = A::narrowOk ? 0 : A::maxCh;
    for (int nch = lo; nch <= A::maxCh; ++nch)
        for (int n : { 0, 1, 2, 63, kMaxBlock, kMaxBlock + 1, 127, 128, 191, 1000, 4099 })
        {
            std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), std::max (n, 1), 0xD0D0u + (std::uint32_t) (nch * 31 + n));
            auto p = planes (b);
            ++cases;
            if (A::run (s, p.data(), nch, n)) ++accepted;
        }
    ok (cases > 0, std::string (A::name()) + ": the accepted sweep is non-empty (precondition)");
    ok (accepted == cases, std::string (A::name()) + ": all " + std::to_string (cases) + " accepted inputs are HONOURED");
}

//==================================================================================================
// P3 — a refused call is indistinguishable from one never made. The buffer test above cannot see
// the object's STATE; this one can, because the stream that follows would diverge if it had moved.
//==================================================================================================
template <class A>
void refusalIsInert()
{
    typename A::T a, b;
    ok (A::prep (a) && A::prep (b), std::string (A::name()) + ": both instances prepare (precondition)");

    const int N = 128, blocks = 12;
    std::vector<std::vector<float>> ya, yb;
    double energy = 0.0;
    for (int k = 0; k < blocks; ++k)
    {
        fill (ya, A::maxCh, N, 0xE0E0u + (std::uint32_t) k);
        yb = ya;
        if (k == 5)   // a REFUSED call, on `a` only: one width past the prepared one
        {
            std::vector<std::vector<float>> junk; fill (junk, A::maxCh + 1, N, 0x9999u);
            auto jp = planes (junk);
            ok (! A::run (a, jp.data(), A::maxCh + 1, N), std::string (A::name()) + ": the injected call IS refused (precondition)");
        }
        // The precondition is about the PROBE, so it is measured on the INPUT: a stage with a long FIR
        // (LinearPhaseEq's is N/2 = thousands of samples) is still priming here and its OUTPUT is zeros,
        // which would make a precondition read on the output false about a perfectly live fixture.
        for (int c = 0; c < A::maxCh; ++c) for (int i = 0; i < N; ++i) energy += std::fabs (ya[(std::size_t) c][(std::size_t) i]);
        auto pa = planes (ya); auto pb = planes (yb);
        felitronics::test::run (A::run (a, pa.data(), A::maxCh, N));
        felitronics::test::run (A::run (b, pb.data(), A::maxCh, N));
        if (! bitEqual (ya, yb)) { ok (false, std::string (A::name()) + ": a refused call moved the STATE"); return; }
        // A METER never writes the buffer, so the comparison above is between two copies of the input —
        // the blind fixture this file warns about, in this file. Compare its per-block READING too: a
        // mutation that zeroed the polyphase history and returned false survived five suites without it.
        if constexpr (! A::inPlace)
            if (A::blockObservable (a) != A::blockObservable (b))
            { ok (false, std::string (A::name()) + ": a refused call moved the READING"); return; }
    }
    ok (energy > 1.0, std::string (A::name()) + ": precondition — the probe stream carried signal (sum |x| = "
                      + std::to_string ((long long) energy) + ")");
    ok (true, std::string (A::name()) + ": a refused call leaves the following stream bit-identical (state untouched)");
}

//==================================================================================================
// P4 — chunking equivalence (law 11a). One long call == the caller's own maxBlock-sized calls.
// The lengths deliberately include a non-multiple, so the short final chunk is exercised.
//==================================================================================================
template <class A>
void chunkEquivalence()
{
    if constexpr (! A::chunks) return;
    else
    {
        int cases = 0, equal = 0; double spread = 0.0;
        for (int N : { kMaxBlock + 1, 2 * kMaxBlock, 3 * kMaxBlock + 17, 517 })
        {
            typename A::T one, many;
            if (! (A::prep (one) && A::prep (many))) { ok (false, std::string (A::name()) + ": prepare failed"); return; }

            std::vector<std::vector<float>> a; fill (a, A::maxCh, N, 0xF00Du + (std::uint32_t) N);
            std::vector<std::vector<float>> b = a;
            const auto in = a;

            auto pa = planes (a);
            felitronics::test::run (A::run (one, pa.data(), A::maxCh, N));

            for (int off = 0; off < N; )
            {
                const int m = std::min (N - off, kMaxBlock);
                std::vector<float*> pb;
                for (int c = 0; c < A::maxCh; ++c) pb.push_back (b[(std::size_t) c].data() + off);
                felitronics::test::run (A::run (many, pb.data(), A::maxCh, m));
                off += m;
            }
            ++cases;
            if constexpr (A::inPlace)
            {
                if (bitEqual (a, b)) ++equal;
                for (int c = 0; c < A::maxCh; ++c)
                    for (int i = 0; i < N; ++i)
                        spread = std::max (spread, (double) std::fabs (a[(std::size_t) c][(std::size_t) i] - in[(std::size_t) c][(std::size_t) i]));
            }
            else
            {
                const double r1 = A::observable (one), r2 = A::observable (many);
                if (r1 == r2) ++equal;
                spread = std::max (spread, std::fabs (r1));   // the READING is what moved, not the buffer
            }
        }
        ok (spread > 1e-6, std::string (A::name()) + ": precondition — the stage MOVED the signal (max |out-in| = "
                           + std::to_string (spread) + "), so the comparison is not between two copies of the input");
        ok (equal == cases, std::string (A::name()) + ": a long call is bit-identical to the caller's own "
                            + std::to_string (kMaxBlock) + "-sample calls (" + std::to_string (equal) + "/" + std::to_string (cases) + ")");
    }
}

//==================================================================================================
// P5 — an UNPREPARED instance refuses everything and touches nothing. "Sensible member defaults" are
// not a prepared object: dither::Dither's PCG streams are `state = inc = 0` before prepare(), which
// makes every draw exactly -0.5, so its "TPDF" is a constant -1 LSB — a DC offset, not dither.
//==================================================================================================
template <class A>
void unpreparedRefuses()
{
    typename A::T s;                                   // deliberately NOT prepared
    // TWO stages have no unprepared state to speak of: `stereo::MonoBass` and `stereo::StereoWidth` are
    // pure M/S arithmetic whose default members ARE a valid 48 kHz configuration (prepare() only changes
    // the rate). There is nothing for them to refuse, so the property they owe is the opposite one —
    // that the default state is usable — and asserting a refusal here would be asserting a wish.
    if constexpr (! A::needsPrepare)
    {
        int cases = 0, accepted = 0;
        for (int nch : { 0, 1, A::maxCh })
            for (int n : { 0, 1, 64, 1000 })
            {
                std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), std::max (n, 1), 0xA11Cu);
                auto p = planes (b);
                ++cases;
                if (A::run (s, p.data(), nch, n)) ++accepted;
            }
        ok (accepted == cases, std::string (A::name()) + ": its DEFAULT state is a valid configuration, so all "
                               + std::to_string (cases) + " calls are accepted without a prepare()");
        return;
    }
    int cases = 0, refused = 0, untouched = 0;
    for (int nch : { 0, 1, A::maxCh })
        for (int n : { 0, 1, 64, 1000 })
        {
            std::vector<std::vector<float>> b; fill (b, std::max (nch, 1), std::max (n, 1), 0xA11Cu + (std::uint32_t) (nch * 7 + n));
            const auto before = b;
            auto p = planes (b);
            ++cases;
            if (! A::run (s, p.data(), nch, n)) ++refused;
            if (bitEqual (b, before)) ++untouched;
        }
    ok (cases > 0, std::string (A::name()) + ": the unprepared sweep is non-empty (precondition)");
    ok (refused == cases,  std::string (A::name()) + ": an UNPREPARED instance refuses all "
                           + std::to_string (cases) + " calls and says so");
    ok (untouched == cases, std::string (A::name()) + ": ...and touches nothing");
}

//==================================================================================================
// P6 — A ZERO-WIDTH STRETCH IS A GAP, AND WHAT COMES BACK IS SILENCE. Law 11a as a property rather
// than as one named case. This is the clause the diff-pass consilium found six modules violating,
// every one of them a COMPOSITE that answered `true` at nch == 0 instead of passing the gap down:
// MultibandCompressor -10.2 dBFS out of digital silence, MultibandWidth -6.06, LinearPhaseEq -4.4,
// DeEsser -14.8, RigPlayer -12.7, TruePeakMeter -0.92 dBTP. One named case covered ONE of them.
//==================================================================================================
template <class A>
void zeroWidthGapIsSilent()
{
    if constexpr (! A::narrowOk || ! A::gapVisible) return;   // no gap, or nothing that can show one
    else
    {
        typename A::T s;
        if (! A::prep (s)) { ok (false, std::string (A::name()) + ": prepare failed"); return; }

        // TAIL runs past dither's 4096-sample auto-blank window on purpose: at 24 blocks the property
        // was empty for Dither — a mutation that swallowed the gap passed, because the tail never got
        // far enough for the blank to engage and the tolerance covered the dither noise itself.
        const int N = 64, TONE = 48, GAP = 48, TAIL = 80;
        std::vector<std::vector<float>> t;
        double drive = 0.0;
        for (int k = 0; k < TONE; ++k)                 // a real programme, so there is something to freeze
        {
            fill (t, A::maxCh, N, 0x60A7u + (std::uint32_t) k);
            auto p = planes (t);
            felitronics::test::run (A::run (s, p.data(), A::maxCh, N));
            for (int c = 0; c < A::maxCh; ++c) for (int i = 0; i < N; ++i) drive += std::fabs (t[(std::size_t) c][(std::size_t) i]);
        }
        ok (drive > 1.0, std::string (A::name()) + ": precondition — the programme carried signal (sum |x| = "
                         + std::to_string ((long long) drive) + ")");

        std::vector<std::vector<float>> z ((std::size_t) A::maxCh, std::vector<float> ((std::size_t) N, 0.0f));
        auto zp = planes (z);
        for (int k = 0; k < GAP; ++k) felitronics::test::run (A::run (s, zp.data(), 0, N));   // THE GAP

        double leak = 0.0;
        for (int k = 0; k < TAIL; ++k)                 // ...and now DIGITAL SILENCE goes in
        {
            std::vector<std::vector<float>> q ((std::size_t) A::maxCh, std::vector<float> ((std::size_t) N, 0.0f));
            auto qp = planes (q);
            felitronics::test::run (A::run (s, qp.data(), A::maxCh, N));
            if constexpr (A::inPlace)
                for (int c = 0; c < A::maxCh; ++c) for (int i = 0; i < N; ++i) leak = std::max (leak, (double) std::fabs (q[(std::size_t) c][(std::size_t) i]));
            else if (k == 0)
                // A METER writes nothing, so its leak is in the READING — and it has to be the PER-BLOCK
                // one. The cumulative figures (`truePeakDb()`, `shortTermLufs()`) are memory by design;
                // asking them about a gap measures what they exist for, not what leaked through it.
                leak = std::max (leak, std::pow (10.0, A::blockObservable (s) * 0.05));
        }
        ok (leak <= A::silenceTol, std::string (A::name()) + ": digital silence after a zero-width GAP stays silent (leak "
                                   + std::to_string (leak) + ", tolerance " + std::to_string (A::silenceTol) + ")");
    }
}

//==================================================================================================
// P7 — LAW 11(b) ON THE PREPARE SIDE. An observable refusal in process() is worth nothing if prepare()
// already agreed to a width it cannot honour: the caller is told nothing there and then makes a
// perfectly legal call whose surplus planes go out dry. Nothing tested this clause, and three modules
// were still clamping when the mutation stand asked.
//==================================================================================================
template <class A>
void prepareIsBinding()
{
    if constexpr (! A::needsPrepare || ! A::hasPrepareWidth) return;   // no prepare-side width to bind
    else
    {
        int cases = 0, refused = 0;
        // The sweep has to reach the widths between the MODULE's own ceiling and the house maximum —
        // not the fixture's chosen width, which is what `maxCh` is. `MultibandWidth::prepare(..., 4)`
        // sat in exactly that gap: accepted, after which process(io, 4, n) moved the buffer by 0.5358
        // and THEN refused. `prepareCeiling` is core::kMaxChannels unless the module says otherwise.
        std::vector<int> widths { core::kMaxChannels + 1, core::kMaxChannels + 4, 0, -1 };
        for (int w = A::prepareCeiling + 1; w <= core::kMaxChannels; ++w) widths.push_back (w);
        for (int w : widths)
        {
            typename A::T s;
            ++cases;
            if (! A::prepWidth (s, w)) ++refused;
        }
        ok (cases > 0, std::string (A::name()) + ": the prepare sweep is non-empty (precondition)");
        ok (refused == cases, std::string (A::name()) + ": prepare() REFUSES all " + std::to_string (cases)
                              + " widths it cannot honour, instead of clamping them");
    }
}

template <class A>
void allProperties()
{
    group (A::name());
    rejectedDomain<A>();
    acceptedDomain<A>();
    refusalIsInert<A>();
    chunkEquivalence<A>();
    unpreparedRefuses<A>();
    zeroWidthGapIsSilent<A>();
    prepareIsBinding<A>();
}

} // namespace felitronics::test::law11
