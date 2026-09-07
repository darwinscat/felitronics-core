// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

//==================================================================================================
// LAW 11 (DSP-ARCHITECTURE.md §2) — "the call is a request against a prepared capacity".
//
// One suite over EVERY block-level entry point in the core, written as PROPERTIES over the whole
// ACCEPTED domain and the whole REJECTED domain rather than as points. Four properties per stage:
//
//   P1  REJECTED DOMAIN IS REFUSED, AND THE REFUSAL IS OBSERVABLE. Every width above the prepared
//       one, and every negative extent, returns false — and leaves the caller's buffer bit-identical.
//   P2  ACCEPTED DOMAIN IS ACCEPTED. Every width in [0, maxChannels] and every non-negative length,
//       including lengths far past maxBlock, returns true.
//   P3  A REFUSED CALL IS INDISTINGUISHABLE FROM ONE NEVER MADE. The stream that follows a refused
//       call is bit-identical to the stream of an instance that never saw it — so the object's state
//       did not move either, which a buffer comparison alone cannot show.
//   P4  CHUNKING EQUIVALENCE. For a stage that chunks (law 11a), one call of length N is bit-identical
//       to the caller having made the same maxBlock-sized calls itself.
//
// The stage list is the census. Five entry points are driven by their OWN suites instead, and the
// reason is the same in every case — they need a fixture this one cannot build: `nam::NamStage` and
// `rigplayer::RigPlayer` need a loaded model (see NamStageTests / RigPlayerTests),
// `dynamiceq::LaneDynamics` needs an EqBand and a captured sidechain, `mastering::MasteringChain` has
// an EXACT width and its own block-invariance suite, and `neural::NeuralStage` is a template over a
// backend (NeuralTests). Saying "if a module has a process(), it is here" would be false, and the
// diff-pass consilium checked.
//==================================================================================================

#include <felitronics_test.h>

#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/convolution/CabConvolver.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/deesser/DeEsser.h>
#include <felitronics/dither/Dither.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/dynamiceq/LaneDynamics.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/NoiseGate.h>
#include <felitronics/dynamics/TransientShaper.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/lineareq/LinearPhaseEq.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/multiband/MultibandCompressor.h>
#include <felitronics/multiband/MultibandWidth.h>
#include <felitronics/mastering/OfflineRenderer.h>
#include <felitronics/lineareq/NaturalPhaseEq.h>
#include <felitronics/poweramp/PowerAmpStage.h>
#include <felitronics/saturation/Saturator.h>
#include <felitronics/stereo/MonoBass.h>
#include <felitronics/stereo/StereoWidth.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

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

void fill (std::vector<std::vector<float>>& b, int nch, int n, std::uint32_t seed)
{
    Rng r { seed };
    b.assign ((std::size_t) nch, std::vector<float> ((std::size_t) n, 0.0f));
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < n; ++i)
            b[(std::size_t) c][(std::size_t) i] = 0.5f * std::sin (2.0 * kPi * (220.0 + 30.0 * c) * i / kFs) + 0.1f * r.next();
}

std::vector<float*> planes (std::vector<std::vector<float>>& b)
{
    std::vector<float*> p; p.reserve (b.size());
    for (auto& v : b) p.push_back (v.data());
    return p;
}

bool bitEqual (const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        if (a[c].size() != b[c].size()) return false;
        for (std::size_t i = 0; i < a[c].size(); ++i) if (a[c][i] != b[c][i]) return false;
    }
    return true;
}

//==================================================================================================
// The stage adapters. Each says how to prepare it, how to drive it, and which of law 11's clauses
// it is claiming — `narrowOk` false means "a narrower call is meaningless here" (law 11c).
//==================================================================================================
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

ADAPT (A_Saturator, felitronics::saturation::Saturator, kPrepCh, true, true, true,
       { return s.prepare (kFs, kMaxBlock, kPrepCh, 4, 32); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 4, 32); });

ADAPT (A_Limiter, felitronics::limiter::TruePeakLimiter, kPrepCh, true, true, true,
       { return s.prepare (kFs, kMaxBlock, kPrepCh, { 1.0, 4, 32 }); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, { 1.0, 4, 32 }); });

ADAPT (A_Compressor, felitronics::dynamics::Compressor, kPrepCh, true, true, true,
       { return s.prepare (kFs, kMaxBlock, kPrepCh, 10.0); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 10.0); });

ADAPT (A_NoiseGate, felitronics::dynamics::NoiseGate, kPrepCh, true, true, true,
       { return s.prepare (kFs, kMaxBlock, kPrepCh); },
       // The threshold is ABOVE the fixture's level on purpose: an OPEN gate multiplies by 1 and a
       // silently truncated tail is then invisible — the mutation stand proved that with a -30 dB
       // threshold nothing here noticed the clamp coming back.
       { return s.process (io, nch, n, true, 0.0f); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_TransientShaper, felitronics::dynamics::TransientShaper, kPrepCh, true, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh);
         felitronics::dynamics::TransientShaperParams tp;   // the defaults are INERT — dial it in, or the
         tp.attackDb = 9.0; tp.sustainDb = -6.0;              // chunk comparison compares two copies of the input
         s.setParams (tp); return okp; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_DitherBase, felitronics::dither::Dither, kPrepCh, true, true, true,
       { felitronics::dither::DitherParams p; p.bits = 16; s.setParams (p); return s.prepare (kFs, kMaxBlock, kPrepCh); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_DeEsser, felitronics::deesser::DeEsser, kPrepCh, true, true, true,
       { return s.prepare (kFs, kMaxBlock, kPrepCh); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_DynamicEqBand, felitronics::dynamiceq::DynamicEqBand, kPrepCh, true, true, true,
       { return s.prepare (kFs, kPrepCh); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, w); });

ADAPT (A_EqBand, felitronics::eq::EqBand, kPrepCh, true, true, true,
       { felitronics::eq::BandParams p; p.on = true; p.type = felitronics::eq::FilterType::Bell;
         auto& l = p.lane (felitronics::eq::Lane::Stereo); l.on = true; l.freq = 1000.0; l.Q = 1.0; l.gainDb = 9.0;
         const bool okp = s.prepare (kFs, kPrepCh); s.setParams (p); return okp; },
       { return s.processBlock (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, w); });

ADAPT (A_EqEngine, felitronics::eq::EqEngine, kPrepCh, true, true, true,
       { felitronics::eq::BandParams p; p.on = true; p.type = felitronics::eq::FilterType::Bell;
         auto& l = p.lane (felitronics::eq::Lane::Stereo); l.on = true; l.freq = 1000.0; l.Q = 1.0; l.gainDb = 9.0;
         const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh); s.setBand (0, p); return okp; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_MonoBassBase, felitronics::stereo::MonoBass, 2, true, true, true,
       { felitronics::test::run (s.prepare (kFs)); s.setLowWidth (0.0f); s.reset(); return true; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, 0, w); });

ADAPT (A_StereoWidthBase, felitronics::stereo::StereoWidth, 2, true, true, true,
       { felitronics::test::run (s.prepare (kFs)); s.setWidth (1.6f); s.reset(); return true; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, 0, w); });

ADAPT (A_PowerAmpBase, felitronics::poweramp::PowerAmpStage, 2, true, true, true,
       { felitronics::poweramp::Params p; p.driveDb = 8.0f; p.autoComp = 1.0f;
         felitronics::poweramp::Voicing v; s.prepare (kFs, kMaxBlock, 4); s.setParams (p, v); return true; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { (void) w; (void) s; return true; });

ADAPT (A_TruePeakMeterBase, felitronics::analysis::TruePeakMeter, kPrepCh, true, true, false,
       { felitronics::test::run (s.prepare (kFs, kMaxBlock, kPrepCh)); return true; },
       { return s.process ((const float* const*) io, nch, n); },
       { return s.truePeakDb(); },
       { return s.prepare (kFs, kMaxBlock, w); });

ADAPT (A_LoudnessMeterBase, felitronics::analysis::LoudnessMeter, kPrepCh, true, true, false,
       { felitronics::test::run (s.prepare (kFs, kPrepCh)); return true; },
       { return s.process ((const float* const*) io, nch, n); },
       { return s.shortTermLufs(); },
       { return s.prepare (kFs, w); });

using MbWidth4 = felitronics::multiband::MultibandWidth<4>;
ADAPT (A_MultibandWidthBase, MbWidth4, 2, true, true, true,
       { return s.prepare (kFs, kMaxBlock, 2); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w); });

using NatPhaseEq = felitronics::lineareq::NaturalPhaseEq;
ADAPT (A_NaturalPhaseEq, NatPhaseEq, kPrepCh, false, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh, 0, 0.5f);
         felitronics::eq::BandParams b[1]; b[0].on = true; b[0].type = felitronics::eq::FilterType::Bell;
         auto& l = b[0].lane (felitronics::eq::Lane::Stereo); l.on = true; l.freq = 1000.0; l.Q = 1.5; l.gainDb = 6.0;
         return okp && s.setBands (b, 1); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 0, 0.5f); });

using MbComp4 = felitronics::multiband::MultibandCompressor<4>;
ADAPT (A_MultibandComp, MbComp4, kPrepCh, true, true, true,
       // The LOOKAHEAD is what freezes across a gap, and it defaults to 0 — with the default params the
       // zero-width property had nothing to find here, and the mutation that removed the fix survived.
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh, 10.0);
         felitronics::dynamics::CompressorParams cp; cp.lookaheadMs = 5.0; cp.thresholdDb = -30.0; cp.ratio = 4.0;
         for (int b = 0; b < 4; ++b) s.setBandParams (b, cp);
         // ...and NOT at mix 1. The parallel dry sum is `d + mix*(wet - d)`, which at mix == 1 is `wet`
         // EXACTLY in IEEE — so a frozen dry delay cancels itself and the gap property reads zero on a
         // stage that is leaking the whole line. My own probe sat on that value and read 0.
         s.setMix (0.5f);
         return okp; },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 10.0); });

ADAPT (A_LinearPhaseEq, felitronics::lineareq::LinearPhaseEq, kPrepCh, false, true, true,
       { const bool okp = s.prepare (kFs, kMaxBlock, kPrepCh, 0);
         felitronics::eq::BandParams b[1]; b[0].on = true; b[0].type = felitronics::eq::FilterType::Bell;
         auto& l = b[0].lane (felitronics::eq::Lane::Stereo); l.on = true; l.freq = 1000.0; l.Q = 1.5; l.gainDb = 6.0;
         return okp && s.setBands (b, 1); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 0); });

// The engine and the raw matrix convolver are driven directly, not only through CabConvolver: the
// mutation stand caught nothing when their width guards were reverted, because nothing reached them.
using ConvEng = felitronics::convolution::ConvolutionEngine<>;
ADAPT (A_ConvEngine, ConvEng, 2, true, true, true,
       { const bool okp = s.prepare (128, 512, 256, 2);
         std::vector<float> ir (256, 0.0f); ir[0] = 1.0f; ir[7] = 0.5f;
         return okp && s.setIr (ir.data(), (int) ir.size()); },
       { return s.process ((const float* const*) io, io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (128, 512, 256, w); });

using MatConv = felitronics::convolution::MatrixConvolver<>;
ADAPT (A_MatrixConvolverBase, MatConv, 2, false, true, true,
       { const bool okp = s.prepare (128, 512, 256, 2);
         std::vector<float> ir (256, 0.0f); ir[0] = 1.0f; ir[7] = 0.5f;
         return okp && s.setIr (ir.data(), (int) ir.size()); },
       { return s.process ((const float* const*) io, io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (128, 512, 256, w); });

using CabConv = felitronics::convolution::CabConvolver;
inline bool prepCab (CabConv& s)
{
    const bool okp = s.prepare (kFs, kMaxBlock, 2, 0.25, false);
    std::vector<float> ir (512, 0.0f); ir[0] = 1.0f; ir[9] = 0.4f;
    const float* irp[2];
    irp[0] = ir.data(); irp[1] = ir.data();
    s.loadIR (irp, 2, (int) ir.size(), kFs);
    for (int k = 0; k < 200 && s.isBusy(); ++k) (void) s.flushPending();   // let the warm swap land
    return okp;
}
ADAPT (A_CabConvolverBase, CabConv, 2, false, true, true,
       { return prepCab (s); },
       { return s.process (io, nch, n); },
       { (void) s; return 0.0; },
       { return s.prepare (kFs, kMaxBlock, w, 0.25, false); });

// `dither::Dither` is a NOISE SOURCE by construction: its job is to emit into digital silence, and it
// only mutes once autoBlank has counted 4096 silent samples. Its floor is therefore the dither itself
// — about 2 LSB of 16 bit — not zero. Naming that here keeps it inside the property rather than
// exempting it; a frozen-audio leak would be orders of magnitude above this.
struct A_Dither : A_DitherBase { static constexpr double silenceTol = 1.5e-4; };

// Modules whose OWN prepare-side ceiling is two, not core::kMaxChannels — a fixed stereo topology.
struct A_MultibandWidth : A_MultibandWidthBase  { static constexpr int prepareCeiling = 2; };
struct A_CabConvolver   : A_CabConvolverBase    { static constexpr int prepareCeiling = 2; };
struct A_MatrixConvolver: A_MatrixConvolverBase { static constexpr int prepareCeiling = 2; };

// (The exemption that used to sit here — "MonoBass and StereoWidth have no unprepared state" — was
// FALSE and is gone. A default-constructed MonoBass has a crossover with no coefficients, so its fold
// does not fold: measured on a pure side signal at lowWidth 0, the side band came out at -6.02 dB where
// a prepared object kills it to -54.22 — 48.199 dB at 30 Hz. Both stereo stages have a prepared_ gate
// now, and answer P5 like every other stage.)

// `poweramp::PowerAmpStage::prepare` takes no channel count at all — its width is the compile-time
// `kMaxCh`, so there is no prepare-side width to bind and P7 has nothing to ask it.
struct A_PowerAmp : A_PowerAmpBase { static constexpr bool hasPrepareWidth = false; };

// A METER'S GAP IS VISIBLE ONLY IN A PER-BLOCK READING. `truePeakDb()` is a running maximum and
// `shortTermLufs()` is a three-second window: both are SUPPOSED to remember the programme, so asking
// them whether a gap left something behind measures the memory they exist for. TruePeakMeter has a
// per-block figure and answers honestly; LoudnessMeter has none, and says so.
struct A_TruePeakMeter : A_TruePeakMeterBase
{ static double blockObservable (T& s) { return s.truePeakDbBlock(); } };
struct A_LoudnessMeter : A_LoudnessMeterBase
{ static constexpr bool gapVisible = false; };

// The two stages whose DEFAULT state is already a valid configuration — see unpreparedRefuses().
// Their prepare() DOES take a width, and now honours it — the exemption that used to sit here said
// otherwise and the suite rested on it. Their ceiling is two, like every other fixed stereo stage.
struct A_MonoBass    : A_MonoBassBase    { static constexpr int prepareCeiling = 2; };
struct A_StereoWidth : A_StereoWidthBase { static constexpr int prepareCeiling = 2; };

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

} // namespace

//==================================================================================================
// The named cases. Each of these is a mutation the property sweep above did NOT catch — the sweep
// asks "is the contract honoured", and these ask "is the DEFECT the contract was written for gone".
//==================================================================================================
void namedCases()
{
    group ("law 11a — a zero-width stretch is a GAP: the falling edge fires (Compressor lookahead)");
    {
        // The measurement this clause exists for: on untouched main, 5 ms of lookahead, a tone, 4800
        // samples of zero-width calls, then stereo DIGITAL SILENCE -> 0.280315 out of the silence
        // (-11.05 dBFS), last non-zero at sample 239 = the whole 240-sample line replaying.
        felitronics::dynamics::Compressor c;
        ok (c.prepare (kFs, 64, 2, 50.0), "precondition: prepared");
        felitronics::dynamics::CompressorParams cp; cp.lookaheadMs = 5.0; c.setParams (cp);
        const int TONE = 4800, GAP = 4800, TAIL = 512;   // every length here is a whole number of 64-sample
                                                        // blocks: 480 was not, and the last block ran 32
                                                        // floats past the end of the plane (the wasm tier
                                                        // caught it; ASan had not been re-run since)
        std::vector<std::vector<float>> t; fill (t, 2, TONE, 0x2468u);
        double inEnergy = 0.0;
        for (int i = 0; i < TONE; ++i) inEnergy += std::fabs (t[0][(std::size_t) i]);
        ok (inEnergy > 100.0, "precondition: the tone carried signal (sum |x| = " + std::to_string ((long long) inEnergy) + ")");
        for (int off = 0; off < TONE; off += 64)
        { float* sub[2] { t[0].data() + off, t[1].data() + off }; felitronics::test::run (c.process (sub, 2, 64)); }

        std::vector<float> z ((std::size_t) 64, 0.0f);
        float* zio[2] { z.data(), z.data() };
        for (int off = 0; off < GAP; off += 64) felitronics::test::run (c.process (zio, 0, 64));

        std::vector<std::vector<float>> sil ((std::size_t) 2, std::vector<float> ((std::size_t) TAIL, 0.0f));
        for (int off = 0; off < TAIL; off += 64)
        { float* sub[2] { sil[0].data() + off, sil[1].data() + off }; felitronics::test::run (c.process (sub, 2, 64)); }
        double pk = 0.0; for (int i = 0; i < TAIL; ++i) pk = std::max (pk, (double) std::fabs (sil[0][(std::size_t) i]));
        ok (pk == 0.0, "digital silence in, digital silence out: the gap dropped the lookahead lines (was 0.280315)");
    }

    group ("law 11d — EqEngine spends audio time on a zero-width call, exactly as EqBand does");
    {
        // Before: 10240 samples of (nch = 0) calls left the two 11.08 dB apart on the same glide.
        auto mk = [] (double f) {
            felitronics::eq::BandParams p; p.on = true; p.type = felitronics::eq::FilterType::Bell;
            auto& l = p.lane (felitronics::eq::Lane::Stereo); l.on = true; l.freq = f; l.Q = 1.0; l.gainDb = 12.0;
            return p; };
        felitronics::eq::EqBand   band; ok (band.prepare (kFs, 2), "precondition: band prepared");   // EqEngine::prepare
        // gives its bands EqBand::prepare's DEFAULT smoothing, so the standalone one must take it too — else
        // this measures two different glide rates instead of two different clocks.
        felitronics::eq::EqEngine eng;  ok (eng.prepare (kFs, 512, 2),   "precondition: engine prepared");
        band.setParams (mk (500.0)); eng.setBand (0, mk (500.0));
        std::vector<float> warm ((std::size_t) 64, 0.0f);
        float* wio[2] { warm.data(), warm.data() };
        for (int k = 0; k < 200; ++k) { felitronics::test::run (band.processBlock (wio, 2, 64)); felitronics::test::run (eng.process (wio, 2, 64)); }
        const double w = 2.0 * kPi * 4000.0 / kFs;
        const double b0 = 20.0 * std::log10 (std::max (1e-12, std::abs (band.response (w, felitronics::eq::Axis::Mid))));
        const double e0 = eng.magnitudeDb (4000.0, felitronics::eq::Axis::Mid);
        band.setParams (mk (4000.0)); eng.setBand (0, mk (4000.0));
        std::vector<float> z ((std::size_t) 512, 0.0f);
        float* zio[2] { z.data(), z.data() };
        for (int k = 0; k < 20; ++k) { felitronics::test::run (band.processBlock (zio, 0, 512)); felitronics::test::run (eng.process (zio, 0, 512)); }
        const double b1 = 20.0 * std::log10 (std::max (1e-12, std::abs (band.response (w, felitronics::eq::Axis::Mid))));
        const double e1 = eng.magnitudeDb (4000.0, felitronics::eq::Axis::Mid);
        ok (std::fabs (b1 - b0) > 1.0, "precondition: the compared quantity MOVED (band glided "
                                       + std::to_string (b1 - b0) + " dB over the zero-width stretch)");
        ok (std::fabs ((b1 - b0) - (e1 - e0)) < 1e-9,
            "engine and band spend the same audio time on a zero-width call (was 11.08 dB apart)");
    }

    group ("law 11d — MonoBass does not fire its bypass edge on an EMPTY call");
    {
        felitronics::stereo::MonoBass a, b;
        felitronics::test::run (a.prepare (kFs)); felitronics::test::run (b.prepare (kFs));
        a.setLowWidth (0.0f); b.setLowWidth (0.0f); a.reset(); b.reset();
        std::vector<std::vector<float>> xa, xb;
        double e = 0.0; bool equal = true;
        for (int k = 0; k < 8; ++k)
        {
            fill (xa, 2, 128, 0x5150u + (std::uint32_t) k); xb = xa;
            // The edge that used to fire lives in the BYPASS branch, so the empty call has to reach it:
            // a stereo-width call at lowWidth 0 is NOT bypassed and never gets there. `nch = 1, n = 0`.
            if (k == 3) felitronics::test::run (a.process (planes (xa).data(), 1, 0));   // the EMPTY call
            auto pa = planes (xa); auto pb = planes (xb);
            felitronics::test::run (a.process (pa.data(), 2, 128));
            felitronics::test::run (b.process (pb.data(), 2, 128));
            for (int i = 0; i < 128; ++i) e += std::fabs (xa[0][(std::size_t) i]);
            equal = equal && bitEqual (xa, xb);
        }
        ok (e > 1.0, "precondition: the stream carried signal");
        ok (equal, "an n == 0 call is a TRUE no-op: it did not clear the crossover");
    }

    group ("law 11a — NoiseGate: phase B may not reach past what phase A produced");
    {
        felitronics::dynamics::NoiseGate g;
        ok (g.prepare (kFs, 512, 2), "precondition: prepared for maxBlock 512");
        std::vector<float> quiet ((std::size_t) 512, 1.0e-4f), loud ((std::size_t) 512, 0.5f);
        const float* qk[2] { quiet.data(), quiet.data() };
        const float* lk[2] { loud.data(),  loud.data()  };
        for (int k = 0; k < 200; ++k) felitronics::test::run (g.analyse (qk, 2, 512, true, -30.0f));
        ok (g.analysedSamples() == 512, "precondition: 512 samples of curve exist");
        felitronics::test::run (g.analyse (lk, 2, 256, true, -30.0f));
        ok (g.analysedSamples() == 256, "a shorter analyse() shortens the valid curve");
        std::vector<float> unity ((std::size_t) 512, 1.0f);
        float* io[1] { unity.data() };
        ok (! g.applyGain (io, 1, 512), "applying 512 of a 256-sample curve is REFUSED (was: 90.0 dB on unanalysed material)");
        bool untouched = true; for (float v : unity) untouched = untouched && (v == 1.0f);
        ok (untouched, "...and the refused apply touched nothing");
        ok (! g.analyse (qk, 2, 513, true, -30.0f), "an over-long analyse() is refused (its result must outlive the call)");
        // ...and being refused, it is INERT — law 11's own invariant. The curve it did not produce does not
        // erase the one that is there. That is safe in the pairing this API exists for, because phase B is
        // called with the same n and is refused by the same bound.
        ok (g.analysedSamples() == 256, "...and a REFUSED analyse leaves the previous curve exactly as it was");
        ok (! g.applyGain (io, 1, 512), "...so the over-long apply is still refused by the 256-sample bound");
        ok (g.applyGain (io, 1, 256),   "...and the 256 that ARE analysed still apply");
    }

    group ("law 11b — CabConvolver::prepare refuses a width it cannot honour, and a refusal disarms it");
    {
        felitronics::convolution::CabConvolver cab;
        ok (! cab.prepare (kFs, 64, 4, 0.25, false), "prepare(numChannels = 4) is REFUSED (its ceiling is 2, not kMaxChannels)");
        ok (! cab.prepare (kFs, 64, 0, 0.25, false), "prepare(numChannels = 0) is REFUSED");
        ok (cab.prepare (kFs, 64, 2, 0.25, false),   "...and 2 is accepted");
        // ...and a LATER refusal must leave it unprepared rather than quietly keeping the old build:
        // process() has to read that flag, which it did not.
        ok (! cab.prepare (kFs, 64, 4, 0.25, false), "a re-prepare at an impossible width is refused");
        std::vector<std::vector<float>> b; fill (b, 2, 64, 0xCAB1u);
        auto p = planes (b);
        ok (! cab.process (p.data(), 2, 64), "...and the object is now UNPREPARED, so process() refuses too");
    }

    group ("law 11 — a refused band must not move LaneDynamics' control seams");
    {
        felitronics::dynamiceq::LaneDynamics ld;
        ok (ld.prepare (kFs, 2), "precondition: the lane layer is prepared for 2");
        felitronics::eq::EqBand band;
        ok (band.prepare (kFs, 1), "precondition: the BAND is prepared for 1 — the geometry disagrees");
        felitronics::eq::BandParams lp; lp.on = true; lp.type = felitronics::eq::FilterType::Bell;
        lp.dyn.on = true; lp.dyn.rangeDb = -6.0; lp.dyn.thrDb = -40.0; lp.dyn.thrAuto = false;
        auto& ll = lp.lane (felitronics::eq::Lane::Stereo);
        ll.on = true; ll.freq = 1000.0; ll.Q = 1.0; ll.gainDb = 0.0;
        ld.setParams (lp); band.setParams (lp);
        std::vector<std::vector<float>> a; fill (a, 2, 512, 0x1A2Bu);
        auto ap = planes (a);
        std::vector<const float*> sc { a[0].data(), a[1].data() };
        const double before = ld.deltaDb (felitronics::eq::Lane::Stereo);
        ok (! ld.processBand (ap.data(), sc.data(), 2, 512, band), "the call is REFUSED (the band cannot take 2)");
        ok (ld.deltaDb (felitronics::eq::Lane::Stereo) == before,
            "...and the control seam did not move either (was 0 -> -0.0345 dB on a refused call)");
    }

    group ("law 11 — EqEngine::captureSectionInput is the OTHER two-phase API, and owes the same");
    {
        felitronics::eq::EqEngine eng;
        ok (eng.prepare (kFs, 64, 2), "precondition: prepared for 2 channels, maxBlock 64");
        std::vector<std::vector<float>> b; fill (b, 3, 64, 0xEC01u);
        std::vector<const float*> cp { b[0].data(), b[1].data(), b[2].data() };
        ok (eng.captureSectionInput (cp.data(), 2, 32) != nullptr, "a capture inside the domain is honoured");
        ok (eng.sectionInputSamples() == 32, "...and records its length");
        ok (eng.captureSectionInput (cp.data(), 3, 32) == nullptr, "a capture WIDER than prepared is REFUSED, not narrowed");
        ok (eng.sectionInputSamples() == 32, "...and being refused, it left the previous capture intact");
        ok (eng.captureSectionInput (cp.data(), 2, 65) == nullptr, "a capture LONGER than maxBlock is refused");
        ok (eng.sectionInputSamples() == 32, "...and that refusal is inert too (was: a valid 32 became 0)");
    }

    group ("law 11 — NoiseGate: n == 0 is the ONE true no-op, even in phase A");
    {
        felitronics::dynamics::NoiseGate g;
        ok (g.prepare (kFs, 512, 2), "precondition: prepared");
        std::vector<float> key ((std::size_t) 512, 0.2f);
        const float* kp[2] { key.data(), key.data() };
        felitronics::test::run (g.analyse (kp, 2, 64, true, -60.0f));
        ok (g.analysedSamples() == 64, "precondition: a 64-sample curve exists");
        felitronics::test::run (g.analyse (kp, 2, 0, true, -60.0f));
        ok (g.analysedSamples() == 64, "an n == 0 analyse changes NOTHING — not even the curve's length");
        std::vector<float> unity ((std::size_t) 64, 1.0f);
        float* io[1] { unity.data() };
        ok (g.applyGain (io, 1, 64), "...so phase B still has its 64 samples to apply");
    }
}

//==================================================================================================
// The named cases the SECOND diff pass produced. Each one is a fix whose absence a consilium seat
// measured and which no property here could see.
//==================================================================================================
void secondPassCases()
{
    group ("law 11b — a refused prepare() adopts nothing and leaves the object UNUSABLE");
    {
        // The two halves are only compatible in one order — disarm, validate, write — and the
        // observable contract is these three lines, not "the old fields survive": storing the width
        // before validating the block size was a heap-buffer-overflow in render().
        felitronics::mastering::OfflineRenderer r;
        ok (r.prepare (1, 256), "precondition: a 1-channel, 256-sample renderer");
        ok (! r.prepare (2, 0), "prepare(2, 0) is refused — the block size is impossible");
        ok (r.maxChannels() != 2, "...and the REFUSED call's own width was not adopted");
        felitronics::mastering::MasteringChain chain;
        ok (chain.prepare (kFs, 1), "precondition: a 1-channel chain to render through");
        std::vector<float> buf ((std::size_t) 256, 0.1f);
        float* p[1] { buf.data() };
        ok (! r.render (chain, (const float* const*) p, p, 1, 256),
            "...and the renderer is UNUSABLE until a prepare() succeeds (it used to write past its scratch)");
        ok (r.prepare (1, 256) && r.render (chain, (const float* const*) p, p, 1, 256),
            "...and a successful prepare() brings it back");
    }

    group ("law 11b — a refused prepare() leaves the object UNPREPARED, not armed on the old build");
    {
        felitronics::multiband::MultibandCompressor<4> mc;
        ok (mc.prepare (kFs, kMaxBlock, 2, 10.0), "precondition: prepared for 2");
        ok (! mc.prepare (kFs, kMaxBlock, 0, 10.0), "prepare(0) is refused");
        std::vector<std::vector<float>> b; fill (b, 2, kMaxBlock, 0xAAB1u);
        auto p = planes (b);
        ok (! mc.process (p.data(), 2, kMaxBlock), "...and the composite is DISARMED, not still running the old build");
    }

    group ("law 11a — MultibandProcessor drops the DRY line too, and tells its bypassed bands");
    {
        // At mix 1 the dry sum cancels exactly; this runs at 0.5, where the whole 240-sample line used
        // to replay: 0.25 out of digital silence (-12.0 dBFS), last non-zero at sample 239.
        for (int bypassDuringGap = 0; bypassDuringGap < 2; ++bypassDuringGap)
        {
            felitronics::multiband::MultibandCompressor<4> mc;
            ok (mc.prepare (kFs, 64, 2, 10.0), "precondition: prepared");
            felitronics::dynamics::CompressorParams cp; cp.lookaheadMs = 5.0; cp.thresholdDb = -30.0; cp.ratio = 4.0;
            for (int b = 0; b < 4; ++b) mc.setBandParams (b, cp);
            mc.setMix (0.5f);
            std::vector<std::vector<float>> t;
            for (int k = 0; k < 75; ++k)
            {
                fill (t, 2, 64, 0x51D0u + (std::uint32_t) k);
                auto p = planes (t);
                felitronics::test::run (mc.process (p.data(), 2, 64));
            }
            if (bypassDuringGap) for (int b = 0; b < 4; ++b) mc.setBandBypass (b, true);
            std::vector<std::vector<float>> z ((std::size_t) 2, std::vector<float> (64, 0.0f));
            auto zp = planes (z);
            for (int k = 0; k < 75; ++k) felitronics::test::run (mc.process (zp.data(), 0, 64));
            if (bypassDuringGap) for (int b = 0; b < 4; ++b) mc.setBandBypass (b, false);
            double leak = 0.0;
            for (int k = 0; k < 8; ++k)
            {
                std::vector<std::vector<float>> q ((std::size_t) 2, std::vector<float> (64, 0.0f));
                auto qp = planes (q);
                felitronics::test::run (mc.process (qp.data(), 2, 64));
                for (auto& v : q) for (float x : v) leak = std::max (leak, (double) std::fabs (x));
            }
            ok (leak == 0.0, std::string ("digital silence after a gap stays silent at mix 0.5")
                             + (bypassDuringGap ? " (bands BYPASSED during the gap)" : "")
                             + " — leak " + std::to_string (leak));
        }
    }

    group ("law 11a — LoudnessMeter's K-weighting is per-channel memory, like TruePeakMeter's history");
    {
        felitronics::analysis::LoudnessMeter lm;
        ok (lm.prepare (kFs, 2), "precondition: prepared for 2");
        std::vector<std::vector<float>> t;
        for (int k = 0; k < 750; ++k)                       // one second of 60 Hz at 0.9
        {
            t.assign (2, std::vector<float> (64, 0.0f));
            for (int c = 0; c < 2; ++c) for (int i = 0; i < 64; ++i)
                t[(std::size_t) c][(std::size_t) i] = 0.9f * (float) std::sin (2.0 * kPi * 60.0 * (k * 64 + i) / kFs);
            auto p = planes (t);
            felitronics::test::run (lm.process ((const float* const*) p.data(), 2, 64));
        }
        std::vector<std::vector<float>> z ((std::size_t) 2, std::vector<float> (64, 0.0f));
        auto zp = planes (z);
        for (int k = 0; k < 750; ++k) felitronics::test::run (lm.process ((const float* const*) zp.data(), 0, 64));
        for (int k = 0; k < 320; ++k) felitronics::test::run (lm.process ((const float* const*) zp.data(), 2, 64));
        ok (lm.momentaryLufs() < -100.0,
            "digital silence after a gap reads silence, not the frozen filter spilling out (was -29.19 LUFS): "
            + std::to_string (lm.momentaryLufs()));
    }

    group ("law 11 — LaneDynamics passes the gap to its own lanes, not only to the band");
    {
        felitronics::dynamiceq::LaneDynamics ld;
        ok (ld.prepare (kFs, 2), "precondition: prepared");
        felitronics::eq::EqBand band;
        ok (band.prepare (kFs, 2), "precondition: the band matches");
        felitronics::eq::BandParams bp; bp.on = true; bp.type = felitronics::eq::FilterType::Bell;
        bp.dyn.on = true; bp.dyn.rangeDb = -12.0; bp.dyn.thrDb = -40.0; bp.dyn.thrAuto = false;
        auto& l = bp.lane (felitronics::eq::Lane::Stereo);
        l.on = true; l.freq = 1000.0; l.Q = 1.0; l.gainDb = 0.0;
        ld.setParams (bp); band.setParams (bp);
        std::vector<std::vector<float>> t;
        for (int k = 0; k < 375; ++k)
        {
            t.assign (2, std::vector<float> (128, 0.0f));
            for (int c = 0; c < 2; ++c) for (int i = 0; i < 128; ++i)
                t[(std::size_t) c][(std::size_t) i] = 0.7f * (float) std::sin (2.0 * kPi * 1000.0 * (k * 128 + i) / kFs);
            auto p = planes (t);
            std::vector<const float*> sc { t[0].data(), t[1].data() };
            felitronics::test::run (ld.processBand (p.data(), sc.data(), 2, 128, band));
        }
        const double earned = ld.deltaDb (felitronics::eq::Lane::Stereo);
        ok (earned < -1.0, "precondition: the lane earned real reduction (" + std::to_string (earned) + " dB)");
        std::vector<std::vector<float>> z ((std::size_t) 2, std::vector<float> (128, 0.0f));
        auto zp = planes (z);
        std::vector<const float*> zsc { z[0].data(), z[1].data() };
        for (int k = 0; k < 375; ++k) felitronics::test::run (ld.processBand (zp.data(), zsc.data(), 0, 128, band));
        ok (std::fabs (ld.deltaDb (felitronics::eq::Lane::Stereo)) < 0.01,
            "a second of GAP releases the lane, as a second of silence would (was held at -12.000 dB): "
            + std::to_string (ld.deltaDb (felitronics::eq::Lane::Stereo)));
    }

    group ("law 11d — captureSectionInput: n == 0 is the ONE true no-op there too");
    {
        felitronics::eq::EqEngine eng;
        ok (eng.prepare (kFs, 64, 2), "precondition: prepared");
        std::vector<std::vector<float>> b; fill (b, 2, 64, 0xC0DEu);
        std::vector<const float*> cp { b[0].data(), b[1].data() };
        ok (eng.captureSectionInput (cp.data(), 2, 32) != nullptr, "a capture of 32 is honoured");
        ok (eng.sectionInputSamples() == 32, "precondition: 32 recorded");
        ok (eng.captureSectionInput (cp.data(), 2, 0) == nullptr, "an EMPTY capture returns nothing to read");
        ok (eng.sectionInputSamples() == 32, "...and changes nothing — n == 0 is a no-op, not a wipe");
    }
}

//==================================================================================================
// The named cases the PR round produced — a THIRD wave, and again every one of them was a regression
// introduced while fixing the wave before. The shape repeats: a fix without a test is a fix that the
// next edit removes for free.
//==================================================================================================
void prRoundCases()
{
    group ("law 11d — an EMPTY call changes no measurement (LoudnessMeter)");
    {
        // `process(io, 0, 0)` is the ONE true no-op. It used to reach the falling edge, and the edge
        // used to clear the sub-hop accumulator as well as the filter — so an empty call inserted into
        // a tone moved the reading by -3.02 LU, one zero-width sample turned -9.11 LUFS into -120, and
        // the evidence of a non-finite input went from 1 to 0.
        felitronics::analysis::LoudnessMeter a, b;
        ok (a.prepare (kFs, 2) && b.prepare (kFs, 2), "precondition: both prepared");
        std::vector<std::vector<float>> t;
        for (int k = 0; k < 200; ++k)
        {
            t.assign (2, std::vector<float> (64, 0.0f));
            for (int c = 0; c < 2; ++c) for (int i = 0; i < 64; ++i)
                t[(std::size_t) c][(std::size_t) i] = 0.5f * (float) std::sin (2.0 * kPi * 997.0 * (k * 64 + i) / kFs);
            auto p = planes (t);
            if (k == 100) felitronics::test::run (a.process ((const float* const*) p.data(), 0, 0));   // THE EMPTY CALL
            felitronics::test::run (a.process ((const float* const*) p.data(), 2, 64));
            felitronics::test::run (b.process ((const float* const*) p.data(), 2, 64));
        }
        ok (a.momentaryLufs() > -60.0, "precondition: the tone registered (" + std::to_string (a.momentaryLufs()) + " LUFS)");
        ok (a.momentaryLufs() == b.momentaryLufs(),
            "an empty call changed the measurement by nothing at all (was -3.02 LU)");
    }

    group ("law 11a — a falling edge drops MEMORY, never measurement already counted");
    {
        // One zero-width SAMPLE in the middle of a sub-hop used to erase the energy counted so far.
        // Compared against ITSELF across the gap, not against a second instance: a one-sample zero-width
        // call spends one sample of audio time, so a reference that did not receive it is a sub-hop out
        // of phase and the difference would be about that instead.
        felitronics::analysis::LoudnessMeter a;
        ok (a.prepare (kFs, 2), "precondition: prepared");
        std::vector<std::vector<float>> t;
        double before = 0.0;
        for (int k = 0; k < 200; ++k)
        {
            t.assign (2, std::vector<float> (64, 0.7f));
            auto p = planes (t);
            if (k == 100) { before = a.momentaryLufs();
                            felitronics::test::run (a.process ((const float* const*) p.data(), 0, 1)); }
            felitronics::test::run (a.process ((const float* const*) p.data(), 2, 64));
        }
        ok (before > -60.0, "precondition: the tone was registering before the gap (" + std::to_string (before) + " LUFS)");
        ok (std::fabs (a.momentaryLufs() - before) < 3.0,
            "a one-sample gap did not erase the sub-hop's energy (it used to collapse to -120): "
            + std::to_string (before) + " -> " + std::to_string (a.momentaryLufs()));

        // The same loss, where it is BINARY rather than a fraction of a 400 ms window: a non-finite input
        // is recorded per sub-hop, and clearing the accumulator on a gap threw the evidence away (1 -> 0).
        felitronics::analysis::LoudnessMeter m;
        ok (m.prepare (kFs, 2), "precondition: prepared");
        std::vector<std::vector<float>> bad ((std::size_t) 2, std::vector<float> (64, 0.1f));
        bad[0][8] = std::numeric_limits<float>::infinity();
        auto bp2 = planes (bad);
        felitronics::test::run (m.process ((const float* const*) bp2.data(), 2, 64));
        felitronics::test::run (m.process ((const float* const*) bp2.data(), 0, 1));   // the gap, mid-sub-hop
        std::vector<std::vector<float>> good ((std::size_t) 2, std::vector<float> (64, 0.1f));
        auto gp = planes (good);
        for (int k = 0; k < 12; ++k) felitronics::test::run (m.process ((const float* const*) gp.data(), 2, 64));
        ok (m.nonFiniteSubHops() >= 1,
            "...and the evidence of a non-finite input survived the gap (it used to go 1 -> 0): "
            + std::to_string (m.nonFiniteSubHops()));
    }

    group ("law 11b — DISARM comes first: a refused re-prepare must not leave the old build answering");
    {
        {
            felitronics::dynamiceq::LaneDynamics ld;
            ok (ld.prepare (kFs, 2), "LaneDynamics: precondition prepared");
            ok (! ld.prepare (kFs, 0), "...a refused re-prepare");
            felitronics::eq::EqBand band; ok (band.prepare (kFs, 2), "precondition: a band to drive");
            std::vector<std::vector<float>> b; fill (b, 2, 64, 0x1D1Du);
            auto p = planes (b);
            std::vector<const float*> sc { b[0].data(), b[1].data() };
            ok (! ld.processBand (p.data(), sc.data(), 2, 64, band), "...disarms it");
        }
        {
            felitronics::multiband::MultibandWidth<4> mw;
            ok (mw.prepare (kFs, kMaxBlock, 2), "MultibandWidth: precondition prepared");
            ok (! mw.prepare (kFs, kMaxBlock, 4), "...a refused re-prepare (its ceiling is 2)");
            std::vector<std::vector<float>> b; fill (b, 2, kMaxBlock, 0x2E2Eu);
            auto p = planes (b);
            ok (! mw.process (p.data(), 2, kMaxBlock), "...disarms it, even though the refusal is its OWN");
        }
    }

    group ("law 11 — LaneDynamics takes the band's verdict BEFORE moving its own state");
    {
        felitronics::dynamiceq::LaneDynamics ld;
        ok (ld.prepare (kFs, 2), "precondition: prepared for 2");
        felitronics::eq::EqBand band;                  // deliberately NOT prepared: the only way to make
                                                       // the band refuse a call this layer will forward
        felitronics::eq::BandParams bp; bp.on = true; bp.type = felitronics::eq::FilterType::Bell;
        bp.dyn.on = true; bp.dyn.rangeDb = -12.0; bp.dyn.thrDb = -40.0; bp.dyn.thrAuto = false;
        auto& l = bp.lane (felitronics::eq::Lane::Stereo);
        l.on = true; l.freq = 1000.0; l.Q = 1.0; l.gainDb = 0.0;
        ld.setParams (bp);
        // EARN a real reduction first, through a band that works — otherwise deltaDb is 0 before and
        // after, and "unchanged" is true of a fixture that proves nothing.
        felitronics::eq::EqBand live; ok (live.prepare (kFs, 2), "precondition: a working band to earn with");
        live.setParams (bp);
        std::vector<std::vector<float>> t;
        for (int k = 0; k < 375; ++k)
        {
            t.assign (2, std::vector<float> (128, 0.0f));
            for (int c = 0; c < 2; ++c) for (int i = 0; i < 128; ++i)
                t[(std::size_t) c][(std::size_t) i] = 0.7f * (float) std::sin (2.0 * kPi * 1000.0 * (k * 128 + i) / kFs);
            auto p = planes (t);
            std::vector<const float*> sc { t[0].data(), t[1].data() };
            felitronics::test::run (ld.processBand (p.data(), sc.data(), 2, 128, live));
        }
        const double before = ld.deltaDb (felitronics::eq::Lane::Stereo);
        ok (before < -1.0, "precondition: the lane earned real reduction (" + std::to_string (before) + " dB)");
        std::vector<std::vector<float>> z ((std::size_t) 2, std::vector<float> (64, 0.0f));
        auto zp = planes (z);
        std::vector<const float*> zsc { z[0].data(), z[1].data() };
        ok (! ld.processBand (zp.data(), zsc.data(), 0, 64, band), "the zero-width call is refused (the band is unprepared)");
        ok (ld.deltaDb (felitronics::eq::Lane::Stereo) == before,
            "...and disengage() did NOT run ahead of that verdict (was -11.9999 -> 0)");
    }

    group ("law 11 — the fused NoiseGate and its two-phase halves answer alike");
    {
        felitronics::dynamics::NoiseGate g;
        ok (g.prepare (kFs, 512, 2), "precondition: prepared");
        std::vector<float> key ((std::size_t) 64, 0.2f);
        const float* kp[2] { key.data(), key.data() };
        felitronics::test::run (g.analyse (kp, 2, 64, true, -60.0f));
        ok (g.analysedSamples() == 64, "precondition: a 64-sample curve exists");
        float* io[1] { key.data() };
        felitronics::test::run (g.process (io, 0, 64, true, -60.0f));    // accepted, but no lanes ran
        ok (g.analysedSamples() == 0,
            "a fused call with no lanes drops the curve, exactly as analyse(key, 0, n) does");
    }

    group ("law 11a — the bands hear a NARROWING, not only a zero-width gap");
    {
        felitronics::multiband::MultibandCompressor<4> mc;
        ok (mc.prepare (kFs, 64, 2, 10.0), "precondition: prepared");
        felitronics::dynamics::CompressorParams cp; cp.lookaheadMs = 5.0; cp.thresholdDb = -30.0; cp.ratio = 4.0;
        for (int b = 0; b < 4; ++b) mc.setBandParams (b, cp);
        mc.setMix (0.5f);
        std::vector<std::vector<float>> t;
        for (int k = 0; k < 75; ++k)
        {
            fill (t, 2, 64, 0x8B8Bu + (std::uint32_t) k);
            auto p = planes (t);
            felitronics::test::run (mc.process (p.data(), 2, 64));
        }
        // The NARROWING is what law 11 owns. (Toggling bypass across the same stretch leaks too — 0.44
        // out of digital silence — but that is P18's "a bypassed cell is a stopped cell" class, older and
        // wider than this law, and it is recorded rather than folded in here.)
        std::vector<std::vector<float>> z ((std::size_t) 2, std::vector<float> (64, 0.0f));
        auto zp = planes (z);
        for (int k = 0; k < 75; ++k) felitronics::test::run (mc.process (zp.data(), 1, 64));   // a NARROWING, not a gap
        // The falling edge is about the plane that STOPPED and came back. Plane 0 never stopped, so a
        // decaying tail there is its own legitimate history — the P18 tests measure the returning plane
        // for exactly this reason.
        double leak = 0.0, stayed = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::vector<std::vector<float>> q ((std::size_t) 2, std::vector<float> (64, 0.0f));
            auto qp = planes (q);
            felitronics::test::run (mc.process (qp.data(), 2, 64));
            for (float x : q[1]) leak   = std::max (leak,   (double) std::fabs (x));
            for (float x : q[0]) stayed = std::max (stayed, (double) std::fabs (x));
        }
        ok (leak == 0.0, "the RETURNING plane brings nothing back after a narrowing (was -16.6378 dBFS): "
                         + std::to_string (leak));
        std::printf ("      [narrowing] the plane that STAYED carries %.6f — its own history, not a leak\n", stayed);

        // ...and the same narrowing with the bands BYPASSED across it: they receive nothing on the
        // ordinary path, so the composite has to hand them the falling edge itself.
        felitronics::multiband::MultibandCompressor<4> mb;
        ok (mb.prepare (kFs, 64, 2, 10.0), "precondition: a second instance");
        for (int b = 0; b < 4; ++b) mb.setBandParams (b, cp);
        mb.setMix (0.5f);
        for (int k = 0; k < 75; ++k)
        {
            fill (t, 2, 64, 0x9C9Cu + (std::uint32_t) k);
            auto p = planes (t);
            felitronics::test::run (mb.process (p.data(), 2, 64));
        }
        for (int b = 0; b < 4; ++b) mb.setBandBypass (b, true);
        for (int k = 0; k < 75; ++k) felitronics::test::run (mb.process (zp.data(), 1, 64));
        for (int b = 0; b < 4; ++b) mb.setBandBypass (b, false);
        double bypLeak = 0.0;
        for (int k = 0; k < 8; ++k)
        {
            std::vector<std::vector<float>> q ((std::size_t) 2, std::vector<float> (64, 0.0f));
            auto qp = planes (q);
            felitronics::test::run (mb.process (qp.data(), 2, 64));
            for (float x : q[1]) bypLeak = std::max (bypLeak, (double) std::fabs (x));
        }
        ok (bypLeak == 0.0, "a BYPASSED band is told about the narrowing too (was -16.6378 dBFS): "
                            + std::to_string (bypLeak));
    }
}

int main()
{
    std::printf ("felitronics law 11 — the caller's contract\n");

    allProperties<A_Saturator>();
    allProperties<A_Limiter>();
    allProperties<A_Compressor>();
    allProperties<A_NoiseGate>();
    allProperties<A_TransientShaper>();
    allProperties<A_Dither>();
    allProperties<A_DeEsser>();
    allProperties<A_DynamicEqBand>();
    allProperties<A_EqBand>();
    allProperties<A_EqEngine>();
    allProperties<A_MonoBass>();
    allProperties<A_StereoWidth>();
    allProperties<A_PowerAmp>();
    allProperties<A_TruePeakMeter>();
    allProperties<A_LoudnessMeter>();
    allProperties<A_MultibandComp>();
    allProperties<A_LinearPhaseEq>();
    allProperties<A_CabConvolver>();
    allProperties<A_ConvEngine>();
    allProperties<A_MatrixConvolver>();
    allProperties<A_MultibandWidth>();
    allProperties<A_NaturalPhaseEq>();

    namedCases();
    secondPassCases();
    prRoundCases();

    return felitronics::test::report();
}
