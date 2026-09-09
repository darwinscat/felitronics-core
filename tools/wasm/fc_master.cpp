// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fc_master — the implementation of the C ABI declared in tools/fc_master_abi.h. Read that header first:
// it carries the contract, and this file carries only the unpacking of it.
//
// NOTHING HERE COMPUTES ANYTHING. Every number handed back was produced by `felitronics::mastering` and
// read out of it. What this file does is four things and no fifth: validate what a page can hand it,
// translate enum codes, copy fields, and own handles. If a line of this file ever looks like DSP, it is
// a defect, and `fcore_master --selftest` is the gate that finds it — the same programme through this
// ABI and through a direct C++ call, one binary, one machine, compared bit for bit.
//
// It compiles natively as well as under emscripten: EMSCRIPTEN_KEEPALIVE degrades to a plain extern "C",
// so the whole validation and addressing layer runs under ctest, ASan and UBSan like anything else.
// That is the same arrangement fc_probe.cpp uses and for the same reason.

#include "fc_master_abi.h"

#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <memory>
#include <new>

#if defined(__EMSCRIPTEN__)
  #include <emscripten/emscripten.h>
  #include <emscripten/heap.h>                 // emscripten_get_heap_size — NOT declared by emscripten.h
  #define FC_EXPORT extern "C" EMSCRIPTEN_KEEPALIVE
#else
  #define FC_EXPORT extern "C"
#endif

using namespace felitronics;
using namespace felitronics::mastering;

namespace
{

//==============================================================================
// LAYOUT PINS
//
// The ABI's job is to survive a change to the C++ structs it mirrors, and "survive" means FAIL TO BUILD
// rather than compile into a silently different meaning. Three mechanisms, because no one of them is
// enough on its own:
//
//  1. `static_assert` on every enum code. A reordering of the C++ enumerators is a change of meaning
//     with no diagnostic anywhere else.
//  2. `static_assert` on the SIZE of each mirrored struct. Catches a field appended.
//  3. A STRUCTURED BINDING with the member names spelled out. This is the one that catches what the
//     others cannot: a `bool` dropped into existing padding does not move `sizeof`, so a size pin sails
//     past it, while a binding of the wrong arity is a hard error — "decomposes into 14 elements, but
//     only 13 names were provided". `CompressorParams` cannot be decomposed at all (it and its base
//     both have members), so that one is pinned by size and by its base's binding instead.
static_assert ((int) eq::FilterType::Bell      == FC_FILTER_BELL);
static_assert ((int) eq::FilterType::LowShelf  == FC_FILTER_LOW_SHELF);
static_assert ((int) eq::FilterType::HighShelf == FC_FILTER_HIGH_SHELF);
static_assert ((int) eq::FilterType::HighPass  == FC_FILTER_HIGH_PASS);
static_assert ((int) eq::FilterType::LowPass   == FC_FILTER_LOW_PASS);
static_assert ((int) eq::FilterType::BandPass  == FC_FILTER_BAND_PASS);
static_assert ((int) eq::FilterType::Notch     == FC_FILTER_NOTCH);
static_assert ((int) eq::FilterType::AllPass   == FC_FILTER_ALL_PASS);
static_assert ((int) eq::FilterType::Tilt      == FC_FILTER_TILT);

static_assert ((int) dynamics::Detector::Peak      == FC_DETECTOR_PEAK);
static_assert ((int) dynamics::Detector::Rms       == FC_DETECTOR_RMS);
static_assert ((int) dynamics::LinkMode::Max       == FC_LINK_MAX);
static_assert ((int) dynamics::LinkMode::MeanPower == FC_LINK_MEAN_POWER);
static_assert ((int) dynamics::Mode::DownCompress  == FC_COMP_DOWN_COMPRESS);
static_assert ((int) dynamics::Mode::UpCompress    == FC_COMP_UP_COMPRESS);
static_assert ((int) dynamics::Mode::DownExpand    == FC_COMP_DOWN_EXPAND);

static_assert ((int) saturation::WaveShaper::Shape::Tanh  == FC_SHAPE_TANH);
static_assert ((int) saturation::WaveShaper::Shape::Atan  == FC_SHAPE_ATAN);
static_assert ((int) saturation::WaveShaper::Shape::Cubic == FC_SHAPE_CUBIC);
static_assert ((int) saturation::WaveShaper::Shape::Asym  == FC_SHAPE_ASYM);

static_assert ((int) dither::NoiseShaping::None          == FC_SHAPING_NONE);
static_assert ((int) dither::NoiseShaping::Weighted      == FC_SHAPING_WEIGHTED);
static_assert ((int) dither::NoiseShaping::Psychoacoustic == FC_SHAPING_PSYCHO);

static_assert ((int) MasteringSolveStatus::Solved                  == FC_SOLVE_SOLVED);
static_assert ((int) MasteringSolveStatus::TargetUnreachable       == FC_SOLVE_TARGET_UNREACHABLE);
static_assert ((int) MasteringSolveStatus::UpstreamViolation       == FC_SOLVE_UPSTREAM_VIOLATION);
static_assert ((int) MasteringSolveStatus::TargetBetweenAchievable == FC_SOLVE_TARGET_BETWEEN);
static_assert ((int) MasteringSolveStatus::PassLimit               == FC_SOLVE_PASS_LIMIT);
static_assert ((int) MasteringSolveStatus::MeasurementInvalid      == FC_SOLVE_MEASUREMENT_INVALID);
static_assert ((int) MasteringSolveStatus::RenderFailed            == FC_SOLVE_RENDER_FAILED);
static_assert ((int) MasteringSolveStatus::NotPrepared             == FC_SOLVE_NOT_PREPARED);
static_assert ((int) MasteringSolveStatus::InvalidRequest          == FC_SOLVE_INVALID_REQUEST);

static_assert ((int) MasteringConstraint::None                    == FC_CONSTRAINT_NONE);
static_assert ((int) MasteringConstraint::TruePeakCeiling         == FC_CONSTRAINT_TRUE_PEAK);
static_assert ((int) MasteringConstraint::LimiterGainReduction    == FC_CONSTRAINT_LIMITER_GR);
static_assert ((int) MasteringConstraint::PeakToLoudness          == FC_CONSTRAINT_PLR);
static_assert ((int) MasteringConstraint::LoudnessRange           == FC_CONSTRAINT_LRA);
static_assert ((int) MasteringConstraint::GainRange               == FC_CONSTRAINT_GAIN_RANGE);
static_assert ((int) MasteringConstraint::CompressorGainReduction == FC_CONSTRAINT_COMPRESSOR_GR);

static_assert ((int) GrStatistic::Mean == FC_GR_MEAN);
static_assert ((int) GrStatistic::P95  == FC_GR_P95);
static_assert ((int) GrStatistic::Max  == FC_GR_MAX);

static_assert (eq::EqEngine::kMaxBands == FC_MAX_EQ_BANDS, "the ABI mirrors every band or it mirrors none");
static_assert (eq::kNumLanes           == FC_MAX_EQ_LANES, "same for lanes");

// The SIZE pins. The comment above promised these and the first draft of this file did not have them,
// which is the "comment that cannot be reproduced" class this repository keeps closing. They are also
// the ONLY pin available for `dynamics::CompressorParams`, which cannot be decomposed at all — it and
// its base both have members — so a field added to it or to `GainReductionParams` moves a number here
// and nothing else in the build.
static_assert (sizeof (eq::LaneParams)                 == 40);
static_assert (sizeof (eq::DynParams)                  == 48);
static_assert (sizeof (eq::BandParams)                 == 264);
static_assert (sizeof (stereo::MonoBassParams)         == 12);
static_assert (sizeof (dynamics::DetectorParams)       == 16);
static_assert (sizeof (dynamics::GainReductionParams)  == 72);
static_assert (sizeof (dynamics::CompressorParams)     == 96);
static_assert (sizeof (saturation::Saturator::Params)  == 28);
static_assert (sizeof (limiter::TruePeakLimiterParams) == 16);
static_assert (sizeof (dither::DitherParams)           == 24);
static_assert (sizeof (MasteringChainConfig)           == 48);
static_assert (sizeof (MasteringChainParams)           == 6544);
static_assert (sizeof (MasteringChainResolved)         == 72);

// The arity pins. Declared in a never-called function so they cost nothing and read as what they are.
[[maybe_unused]] void layoutPins()
{
    eq::LaneParams lane {};
    auto& [l_on, l_freq, l_q, l_gain, l_slope, l_bypass] = lane;
    (void) l_on; (void) l_freq; (void) l_q; (void) l_gain; (void) l_slope; (void) l_bypass;

    eq::DynParams dyn {};
    auto& [d_on, d_range, d_thr, d_auto, d_atk, d_rel] = dyn;
    (void) d_on; (void) d_range; (void) d_thr; (void) d_auto; (void) d_atk; (void) d_rel;

    eq::BandParams band {};
    auto& [b_on, b_type, b_swept, b_bypass, b_dyn, b_lanes] = band;
    (void) b_on; (void) b_type; (void) b_swept; (void) b_bypass; (void) b_dyn; (void) b_lanes;

    stereo::MonoBassParams mb {};
    auto& [mb_en, mb_f, mb_w] = mb;
    (void) mb_en; (void) mb_f; (void) mb_w;

    // `CompressorParams` inherits, so it cannot be decomposed; its BASE can, and the derived part is
    // pinned by size below. Between them, a field added at either level is a build error.
    dynamics::DetectorParams det {};
    auto& [det_d, det_link, det_rms] = det;
    (void) det_d; (void) det_link; (void) det_rms;

    saturation::Saturator::Params clip {};
    auto& [c_shape, c_drive, c_bias, c_mix, c_out, c_auto, c_dc] = clip;
    (void) c_shape; (void) c_drive; (void) c_bias; (void) c_mix; (void) c_out; (void) c_auto; (void) c_dc;

    limiter::TruePeakLimiterParams lim {};
    auto& [lim_ceil, lim_rel] = lim;
    (void) lim_ceil; (void) lim_rel;

    dither::DitherParams dit {};
    auto& [dit_bits, dit_shape, dit_seed, dit_blank, dit_blankn] = dit;
    (void) dit_bits; (void) dit_shape; (void) dit_seed; (void) dit_blank; (void) dit_blankn;

    MasteringChainConfig cfg {};
    auto& [k_block, k_eq, k_mb, k_comp, k_clip, k_lim, k_dith,
           k_clook, k_llook, k_os, k_taps, k_hpf] = cfg;
    (void) k_block; (void) k_eq; (void) k_mb; (void) k_comp; (void) k_clip; (void) k_lim; (void) k_dith;
    (void) k_clook; (void) k_llook; (void) k_os; (void) k_taps; (void) k_hpf;

    MasteringChainParams prm {};
    auto& [p_in, p_pre, p_eq, p_mb, p_comp, p_clip, p_lim, p_dith,
           p_bE, p_bM, p_bC, p_bK, p_bL, p_bD] = prm;
    (void) p_in; (void) p_pre; (void) p_eq; (void) p_mb; (void) p_comp; (void) p_clip; (void) p_lim;
    (void) p_dith; (void) p_bE; (void) p_bM; (void) p_bC; (void) p_bK; (void) p_bL; (void) p_bD;

    MasteringChainResolved res {};
    auto& [r_lat, r_blk, r_clook, r_clip, r_lim, r_llook, r_os, r_ctap, r_ltap,
           r_ceil, r_rel, r_mb] = res;
    (void) r_lat; (void) r_blk; (void) r_clook; (void) r_clip; (void) r_lim; (void) r_llook;
    (void) r_os; (void) r_ctap; (void) r_ltap; (void) r_ceil; (void) r_rel; (void) r_mb;
}

//==============================================================================
// MEMORY VALIDATION — everything below arrives from a page's JavaScript.

// Does [p, p+bytes) lie inside the wasm linear memory? A pointer can be aligned, and its length can fit
// a 32-bit address space, and the span can still run off the end of the heap. This cannot prove the
// caller OWNS the span — no ABI of this shape can — but it turns "the module dies" into "the call is
// refused". Same guard, same wording as fc_probe.cpp; deliberately not generalised into a shared helper,
// because the two files are read separately and a reader of either should see the rule.
bool inHeap (const void* p, std::uint64_t bytes) noexcept
{
#if defined(__EMSCRIPTEN__)
    const std::uint64_t base = (std::uint64_t) reinterpret_cast<std::uintptr_t> (p);
    if (base == 0) return false;
    const std::uint64_t end = base + bytes;
    if (end < base) return false;                                       // wrapped
    return end <= (std::uint64_t) emscripten_get_heap_size();
#else
    (void) p; (void) bytes;
    return true;                                                        // native: no linear memory to bound
#endif
}

bool aligned4 (const void* p) noexcept { return (reinterpret_cast<std::uintptr_t> (p) & 0x3u) == 0; }

// A SCALAR out-parameter is an address from JavaScript too, and it was the one class this file checked
// only for null. `fc_master_latency(h, (int32_t*) heapSize)` then traps the module where every other
// bad address is a refusal — the barrier the facade exists to provide, with a hole in it exactly where
// the value is small enough to look harmless. Natively `inHeap` cannot answer, so this is a wasm-tier
// guard and says so; the alignment half works everywhere.
template <typename T>
fc_status checkScalarOut (const T* p) noexcept
{
    if (p == nullptr) return FC_ERR_NULL;
    if ((reinterpret_cast<std::uintptr_t> (p) & (alignof (T) - 1)) != 0) return FC_ERR_ALIGNMENT;
    if (! inHeap (p, sizeof (T))) return FC_ERR_SPAN;
    return FC_OK;
}

// A planar audio span: `channels * frames` floats starting at `p`. The product is computed in 64 bits
// BEFORE the multiply can wrap, because on wasm32 that product IS the caller's allocation and a wrapped
// one would hand the core a window onto unrelated heap.
fc_status checkAudio (const void* p, std::uint32_t frames, int channels) noexcept
{
    if (p == nullptr) return FC_ERR_NULL;
    if (! aligned4 (p)) return FC_ERR_ALIGNMENT;
    const std::uint64_t bytes = (std::uint64_t) frames * (std::uint64_t) channels * sizeof (float);
    if (bytes > (std::uint64_t) 0xFFFFFFFFu) return FC_ERR_SPAN;
    if (! inHeap (p, bytes)) return FC_ERR_SPAN;
    return FC_OK;
}

// Does a scalar out-parameter sit inside an audio span this call is about to write? The overlap rule
// for `in`/`out` does not see this class at all, and the consequence is silent: `fc_master_flush(h,
// out, D, (uint32_t*) out)` writes the whole drain and THEN overwrites the first sample with the frame
// count. The caller gets audio whose first four bytes are a small integer, and nothing anywhere says so.
bool aliasesSpan (const void* scalar, std::size_t scalarBytes, const void* span, std::uint64_t spanBytes) noexcept
{
    const auto a = (std::uint64_t) reinterpret_cast<std::uintptr_t> (scalar);
    const auto b = (std::uint64_t) reinterpret_cast<std::uintptr_t> (span);
    return (a < b + spanBytes) && (b < a + (std::uint64_t) scalarBytes);
}

// Two planar spans either coincide exactly or do not touch. A PARTIAL overlap is refused rather than
// handled: `memmove` would preserve the audio but destroy part of a buffer the caller declared `const`,
// and a contract that silently eats its own input is worse than one that says no.
bool partiallyOverlaps (const void* a, const void* b, std::uint64_t bytes) noexcept
{
    if (a == b) return false;
    const auto x = (std::uint64_t) reinterpret_cast<std::uintptr_t> (a);
    const auto y = (std::uint64_t) reinterpret_cast<std::uintptr_t> (b);
    return (x < y + bytes) && (y < x + bytes);
}

template <typename T>
fc_status checkHeaderIn (const T* p) noexcept
{
    if (p == nullptr) return FC_ERR_NULL;
    if ((reinterpret_cast<std::uintptr_t> (p) & 0x7u) != 0) return FC_ERR_ALIGNMENT;
    // THE HEADER FIRST, THEN THE REST OF THE SPAN, and the two steps are separate on purpose. Bounding
    // the whole struct before reading the version means a caller who got the version wrong is told
    // FC_ERR_SPAN — the wrong diagnosis, on the one field whose job is to catch exactly that mistake.
    // Eight bytes is what the version costs to read, so eight bytes is what is bounded first.
    if (! inHeap (p, sizeof (fc_header))) return FC_ERR_SPAN;
    if (p->header.abiVersion != FC_MASTER_ABI_VERSION) return FC_ERR_ABI_VERSION;
    if (p->header.structSize != (std::uint32_t) sizeof (T)) return FC_ERR_STRUCT_SIZE;
    if (! inHeap (p, sizeof (T))) return FC_ERR_SPAN;
    return FC_OK;
}

// An OUT struct is checked the same way, because the caller states which layout it expects to be
// written — a facade that wrote its own layout into a buffer sized for another one is the same defect
// with the arrow reversed.
template <typename T>
fc_status checkHeaderOut (T* p) noexcept { return checkHeaderIn (const_cast<const T*> (p)); }

template <typename T>
void stampHeader (T* p) noexcept
{
    p->header.abiVersion = FC_MASTER_ABI_VERSION;
    p->header.structSize = (std::uint32_t) sizeof (T);
}

// A parameter field the core has no verdict for. Ranges are NOT checked here — the core clamps them by
// design and reports what it clamped to — but a non-finite value has no clamp to fall into: it becomes
// 0 dB silently. See the note in fc_master_abi.h.
bool fin (double v) noexcept { return std::isfinite (v); }
bool fin (float v)  noexcept { return std::isfinite (v); }

//==============================================================================
// ENUM TRANSLATION
//
// A code from JS names something or it names nothing; nothing is a refusal, never a clamp to the first
// enumerator. Written as an explicit table rather than a cast so that a C++ enumerator added upstream
// is caught by the static_asserts above rather than silently acquiring a code.
#define FC_MAP_ENUM(fn, CppType, ...)                                          \
    bool fn (std::int32_t code, CppType& out) noexcept                          \
    {                                                                           \
        switch (code) { __VA_ARGS__ default: return false; }                    \
        return true;                                                            \
    }
#define FC_CASE(code, val) case code: out = val; break;

FC_MAP_ENUM (mapFilterType, eq::FilterType,
    FC_CASE (FC_FILTER_BELL,       eq::FilterType::Bell)
    FC_CASE (FC_FILTER_LOW_SHELF,  eq::FilterType::LowShelf)
    FC_CASE (FC_FILTER_HIGH_SHELF, eq::FilterType::HighShelf)
    FC_CASE (FC_FILTER_HIGH_PASS,  eq::FilterType::HighPass)
    FC_CASE (FC_FILTER_LOW_PASS,   eq::FilterType::LowPass)
    FC_CASE (FC_FILTER_BAND_PASS,  eq::FilterType::BandPass)
    FC_CASE (FC_FILTER_NOTCH,      eq::FilterType::Notch)
    FC_CASE (FC_FILTER_ALL_PASS,   eq::FilterType::AllPass)
    FC_CASE (FC_FILTER_TILT,       eq::FilterType::Tilt))

FC_MAP_ENUM (mapDetector, dynamics::Detector,
    FC_CASE (FC_DETECTOR_PEAK, dynamics::Detector::Peak)
    FC_CASE (FC_DETECTOR_RMS,  dynamics::Detector::Rms))

FC_MAP_ENUM (mapLink, dynamics::LinkMode,
    FC_CASE (FC_LINK_MAX,        dynamics::LinkMode::Max)
    FC_CASE (FC_LINK_MEAN_POWER, dynamics::LinkMode::MeanPower))

FC_MAP_ENUM (mapCompMode, dynamics::Mode,
    FC_CASE (FC_COMP_DOWN_COMPRESS, dynamics::Mode::DownCompress)
    FC_CASE (FC_COMP_UP_COMPRESS,   dynamics::Mode::UpCompress)
    FC_CASE (FC_COMP_DOWN_EXPAND,   dynamics::Mode::DownExpand))

FC_MAP_ENUM (mapShape, saturation::WaveShaper::Shape,
    FC_CASE (FC_SHAPE_TANH,  saturation::WaveShaper::Shape::Tanh)
    FC_CASE (FC_SHAPE_ATAN,  saturation::WaveShaper::Shape::Atan)
    FC_CASE (FC_SHAPE_CUBIC, saturation::WaveShaper::Shape::Cubic)
    FC_CASE (FC_SHAPE_ASYM,  saturation::WaveShaper::Shape::Asym))

FC_MAP_ENUM (mapShaping, dither::NoiseShaping,
    FC_CASE (FC_SHAPING_NONE,     dither::NoiseShaping::None)
    FC_CASE (FC_SHAPING_WEIGHTED, dither::NoiseShaping::Weighted)
    FC_CASE (FC_SHAPING_PSYCHO,   dither::NoiseShaping::Psychoacoustic))

FC_MAP_ENUM (mapGrStat, GrStatistic,
    FC_CASE (FC_GR_MEAN, GrStatistic::Mean)
    FC_CASE (FC_GR_P95,  GrStatistic::P95)
    FC_CASE (FC_GR_MAX,  GrStatistic::Max))

#undef FC_CASE
#undef FC_MAP_ENUM

//==============================================================================
// FIELD MAPPING — the other half of what the thinness law permits. Copies, never conversions of meaning:
// a field that is `float` in the core crosses as `float`, because a `double` C-POD field feeding a
// `float` core field adds a narrowing the direct C++ path does not have, and that shows up as a
// bit-exactness failure with no bug behind it.

fc_status toCore (const fc_master_config& c, MasteringChainConfig& out) noexcept
{
    if (! fin (c.sampleRate) || ! fin (c.compressorLookaheadMs)
        || ! fin (c.limiterLookaheadMs) || ! fin (c.sidechainHpfHz)) return FC_ERR_NON_FINITE;
    out.internalBlock         = c.internalBlock;
    out.eq                    = c.eq != 0;
    out.monoBass              = c.monoBass != 0;
    out.compressor            = c.compressor != 0;
    out.clipper               = c.clipper != 0;
    out.limiter               = c.limiter != 0;
    out.dither                = c.dither != 0;
    out.compressorLookaheadMs = c.compressorLookaheadMs;
    out.limiterLookaheadMs    = c.limiterLookaheadMs;
    out.oversampleFactor      = c.oversampleFactor;
    out.tapsPerPhase          = c.tapsPerPhase;
    out.sidechainHpfHz        = c.sidechainHpfHz;
    return FC_OK;
}

fc_status toCore (const fc_master_params& p, MasteringChainParams& out) noexcept
{
    if (! fin (p.inputGainDb) || ! fin (p.preLimiterGainDb)) return FC_ERR_NON_FINITE;

    for (int b = 0; b < FC_MAX_EQ_BANDS; ++b)
    {
        const fc_eq_band& src = p.eqBands[b];
        eq::BandParams&   dst = out.eqBands[b];
        if (! mapFilterType (src.type, dst.type)) return FC_ERR_ENUM;
        dst.on     = src.on     != 0;
        dst.swept  = src.swept  != 0;
        dst.bypass = src.bypass != 0;

        if (! fin (src.dyn.rangeDb) || ! fin (src.dyn.thrDb)
            || ! fin (src.dyn.atk) || ! fin (src.dyn.rel)) return FC_ERR_NON_FINITE;
        dst.dyn.on      = src.dyn.on != 0;
        dst.dyn.rangeDb = src.dyn.rangeDb;
        dst.dyn.thrDb   = src.dyn.thrDb;
        dst.dyn.thrAuto = src.dyn.thrAuto != 0;
        dst.dyn.atk     = src.dyn.atk;
        dst.dyn.rel     = src.dyn.rel;

        for (int l = 0; l < FC_MAX_EQ_LANES; ++l)
        {
            const fc_eq_lane& sl = src.lanes[l];
            eq::LaneParams&   dl = dst.lanes[l];
            if (! fin (sl.freq) || ! fin (sl.q) || ! fin (sl.gainDb)) return FC_ERR_NON_FINITE;
            dl.on     = sl.on != 0;
            dl.freq   = sl.freq;
            dl.Q      = sl.q;
            dl.gainDb = sl.gainDb;
            dl.slope  = sl.slope;
            dl.bypass = sl.bypass != 0;
        }
    }

    if (! fin (p.monoBass.frequencyHz) || ! fin (p.monoBass.lowWidth)) return FC_ERR_NON_FINITE;
    out.monoBass.enabled     = p.monoBass.enabled != 0;
    out.monoBass.frequencyHz = p.monoBass.frequencyHz;
    out.monoBass.lowWidth    = p.monoBass.lowWidth;

    if (! mapDetector (p.compressor.detector, out.compressor.detector)) return FC_ERR_ENUM;
    if (! mapLink (p.compressor.link, out.compressor.link)) return FC_ERR_ENUM;
    if (! mapCompMode (p.compressor.mode, out.compressor.mode)) return FC_ERR_ENUM;
    if (! fin (p.compressor.rmsWindowMs) || ! fin (p.compressor.thresholdDb) || ! fin (p.compressor.ratio)
        || ! fin (p.compressor.kneeDb) || ! fin (p.compressor.rangeDb) || ! fin (p.compressor.attackMs)
        || ! fin (p.compressor.releaseMs) || ! fin (p.compressor.makeupDb)) return FC_ERR_NON_FINITE;
    out.compressor.rmsWindowMs = p.compressor.rmsWindowMs;
    out.compressor.thresholdDb = p.compressor.thresholdDb;
    out.compressor.ratio       = p.compressor.ratio;
    out.compressor.kneeDb      = p.compressor.kneeDb;
    out.compressor.rangeDb     = p.compressor.rangeDb;
    out.compressor.attackMs    = p.compressor.attackMs;
    out.compressor.releaseMs   = p.compressor.releaseMs;
    out.compressor.makeupDb    = p.compressor.makeupDb;
    out.compressor.autoMakeup  = p.compressor.autoMakeup != 0;
    // `lookaheadMs` is not mapped and has no ABI field: the chain overwrites it from the config on every
    // apply, so a field here would be a knob that does nothing.

    if (! mapShape (p.clipper.shape, out.clipper.shape)) return FC_ERR_ENUM;
    if (! fin (p.clipper.driveDb) || ! fin (p.clipper.bias) || ! fin (p.clipper.mix)
        || ! fin (p.clipper.outputDb) || ! fin (p.clipper.autoComp)
        || ! fin (p.clipper.dcBlockHz)) return FC_ERR_NON_FINITE;
    out.clipper.driveDb   = p.clipper.driveDb;
    out.clipper.bias      = p.clipper.bias;
    out.clipper.mix       = p.clipper.mix;
    out.clipper.outputDb  = p.clipper.outputDb;
    out.clipper.autoComp  = p.clipper.autoComp;
    out.clipper.dcBlockHz = p.clipper.dcBlockHz;

    if (! fin (p.limiter.ceilingDbTp) || ! fin (p.limiter.releaseMs)) return FC_ERR_NON_FINITE;
    out.limiter.ceilingDbTp = p.limiter.ceilingDbTp;
    out.limiter.releaseMs   = p.limiter.releaseMs;

    if (! mapShaping (p.dither.shaping, out.dither.shaping)) return FC_ERR_ENUM;
    out.dither.bits             = p.dither.bits;
    // Two halves, because a 64-bit seed does not survive a JS Number: 0x853c49e6748fea9b loses its low
    // bits through a double, and a dither stream seeded from a truncated value is a different render.
    out.dither.seed             = ((std::uint64_t) p.dither.seedHi << 32) | (std::uint64_t) p.dither.seedLo;
    out.dither.autoBlank        = p.dither.autoBlank != 0;
    out.dither.autoBlankSamples = p.dither.autoBlankSamples;

    out.inputGainDb      = p.inputGainDb;
    out.preLimiterGainDb = p.preLimiterGainDb;
    out.bypassEq         = p.bypassEq         != 0;
    out.bypassMonoBass   = p.bypassMonoBass   != 0;
    out.bypassCompressor = p.bypassCompressor != 0;
    out.bypassClipper    = p.bypassClipper    != 0;
    out.bypassLimiter    = p.bypassLimiter    != 0;
    out.bypassDither     = p.bypassDither     != 0;
    return FC_OK;
}

fc_status toCore (const fc_loudness_request& r, LoudnessRequest& out) noexcept
{
    if (! mapGrStat (r.limiterGr.statistic, out.limiterGr.statistic)) return FC_ERR_ENUM;
    if (! mapGrStat (r.compressorGr.statistic, out.compressorGr.statistic)) return FC_ERR_ENUM;
    // NOTE the asymmetry, and it is the request's own contract rather than an oversight: NaN and the
    // infinities are MEANINGFUL here. `limitDb == +inf` is "no limit", `minPlrDb == -inf` is "no limit",
    // `inputLoudnessRangeLu == NaN` is "not supplied", `initialGainDb == NaN` is "use the params' own".
    // A blanket isfinite() sweep over this struct would switch three constraints off by accident. The
    // solver validates each field on its own terms and refuses with `InvalidRequest`.
    out.targetLufs            = r.targetLufs;
    out.toleranceLu           = r.toleranceLu;
    out.maxTruePeakDbTp       = r.maxTruePeakDbTp;
    out.truePeakAimDb         = r.truePeakAimDb;
    out.limiterGr.limitDb     = r.limiterGr.limitDb;
    out.compressorGr.limitDb  = r.compressorGr.limitDb;
    out.minPlrDb              = r.minPlrDb;
    out.maxLraLossLu          = r.maxLraLossLu;
    out.inputLoudnessRangeLu  = r.inputLoudnessRangeLu;
    out.activityThresholdDb   = r.activityThresholdDb;
    out.maxPasses             = r.maxPasses;
    out.initialGainDb         = r.initialGainDb;
    return FC_OK;
}

void fromCore (const MasteringChainResolved& r, int tapOs, fc_master_resolved& out) noexcept
{
    stampHeader (&out);
    out.latencySamples      = r.latencySamples;
    out.internalBlock       = r.internalBlock;
    out.compressorLookahead = r.compressorLookahead;
    out.clipperLatency      = r.clipperLatency;
    out.limiterLatency      = r.limiterLatency;
    out.limiterLookahead    = r.limiterLookahead;
    out.oversampleFactor    = r.oversampleFactor;
    out.compressorTapOffset = r.compressorTapOffset;
    out.limiterTapOffset    = r.limiterTapOffset;
    out.limiterCeilingDbTp  = r.limiterCeilingDbTp;
    out.limiterReleaseMs    = r.limiterReleaseMs;
    out.monoBass.enabled     = r.monoBass.enabled ? 1 : 0;
    out.monoBass.frequencyHz = r.monoBass.frequencyHz;
    out.monoBass.lowWidth    = r.monoBass.lowWidth;
    out.tapOversampleFactor = tapOs;
}

void fromCore (const GainReductionStats& s, fc_gr_stats& out) noexcept
{
    out.meanDb = s.meanDb; out.p95Db = s.p95Db; out.maxDb = s.maxDb;
    out.activeFraction = s.activeFraction;
    out.frames = s.frames; out.nonFinite = s.nonFinite; out.aboveRange = s.aboveRange;
    out.valid = s.valid ? 1 : 0;
}

void fromCore (const MasterMeasurement& m, fc_measurement& out) noexcept
{
    stampHeader (&out);
    out.integratedLufs  = m.integratedLufs;
    out.truePeakDbTp    = m.truePeakDbTp;
    out.samplePeakDb    = m.samplePeakDb;
    out.loudnessRangeLu = m.loudnessRangeLu;
    out.plrDb           = m.plrDb;
    fromCore (m.compressor, out.compressor);
    fromCore (m.limiter,    out.limiter);
    out.limiterMaxReconstructedPeakDb = m.limiterMaxReconstructedPeakDb;
    out.latencySamples   = m.latencySamples;
    out.gatingBlocks     = m.gatingBlocks;
    out.droppedBlocks    = m.droppedBlocks;
    out.nonFiniteSubHops = m.nonFiniteSubHops;
    out.loudnessValid    = m.loudnessValid ? 1 : 0;
    out.lraValid         = m.lraValid ? 1 : 0;
}

//==============================================================================
// HANDLES
//
// An index and a generation, packed into 32 bits, plus the KIND of object the slot holds. Measured
// reason rather than caution: under emscripten, destroying a facade instance and creating another
// returned the SAME address 19 times out of 19 with emmalloc and with dlmalloc, while native macOS
// malloc reused it 0 times out of 20. So a raw-pointer handle is a use-after-free that the developer's
// own machine never reproduces and the shipping tier reproduces always — and a linear memory has no
// unmapped page to trap on.
//
// This is global mutable state, which law 6 forbids — to `modules/`. This file is `tools/`, is
// single-threaded by contract (law 1: the core is called synchronously; the worker owns the loop), and
// the alternative is worse.
constexpr int kMaxHandles = 8;

enum class Kind : std::uint8_t { Free = 0, Master = 1, Solution = 2 };

struct MasterInstance
{
    MasteringChainConfig cfg {};
    MasteringChain       chain;
    OfflineRenderer      renderer;
    TargetLoudnessSolver solver;
    // 64 bits, because 32 overflows on a legal stream: 4096-frame calls wrap `framesIn` after about
    // 24 h 51 min at 48 kHz and 6 h 13 min at 192 kHz, and a counter that can read zero after having
    // been non-zero is worse than no counter.
    std::uint64_t framesIn = 0, framesFlushed = 0;
    bool audioSeen = false;         // set by process(); configure() refuses once this is true

    // A SOLVE THAT RAN LEAVES THE CHAIN HOLDING THE SOLVER'S PARAMETERS, not the caller's, and standing
    // wherever its last pass ended. Measured: `process()` straight after a solve returns FC_OK and
    // differs from the delivered render in 558 691 of 576 000 samples, and from the CONFIGURED render in
    // 575 998 — a third render nobody chose, out of a handle whose stats read `framesIn = 0` as if it
    // were fresh. `reset()` does not fix it either: it clears the audio state and keeps the solver's
    // gain and ceiling. So the handle is marked, and `process`/`flush` refuse until `configure` puts a
    // known parameter set back. Law 11 at the handle level: a call that cannot be honoured as the caller
    // means it is refused rather than answered with something plausible.
    bool solverRan = false;
};

struct Slot
{
    Kind kind = Kind::Free;
    // 24 BITS OF GENERATION, and the width is the whole design rather than a spare-bits accident. An
    // 8-bit generation forces a choice between ABA (reuse the numbers) and RETIREMENT (spend the slot),
    // and retirement turns "8 live objects" into a LIFETIME BUDGET: 8 x 255 = 2040 create/destroy
    // cycles per page load, after which every correct create is refused for ever. That is reachable —
    // the reference CLI's own render loop creates a handle per programme, so a worker written from it
    // would die on its 2041st file. At 24 bits the budget is 8 x 16.7 million and the question stops
    // being one; a stale handle is still refused, which was the point.
    std::uint32_t gen = 1;          // starts at 1 so a zeroed handle is never valid
    std::unique_ptr<MasterInstance> master;
    std::unique_ptr<LoudnessSolution> solution;
};

Slot g_slots[kMaxHandles];

// handle = (gen << 8) | (index + 1). Zero is never valid, so a zeroed variable in JS is a refusal.
std::uint32_t packHandle (int idx, std::uint32_t gen) noexcept
{
    return ((gen & 0x00FFFFFFu) << 8) | (std::uint32_t) (idx + 1);
}

Slot* lookup (std::uint32_t h, Kind want) noexcept
{
    if (h == 0) return nullptr;
    const int idx = (int) (h & 0xFFu) - 1;
    if (idx < 0 || idx >= kMaxHandles) return nullptr;
    Slot& s = g_slots[idx];
    if (s.kind != want) return nullptr;
    // The generation IS the "never issued" test now that it fills the rest of the word: a fabricated
    // value differs from the live one in one of 24 bits rather than sharing the object with 65 535
    // other spellings, which is what an 8-bit generation left behind.
    if ((s.gen & 0x00FFFFFFu) != ((h >> 8) & 0x00FFFFFFu)) return nullptr;
    return &s;
}

int allocSlot (Kind kind, std::uint32_t& outHandle) noexcept
{
    for (int i = 0; i < kMaxHandles; ++i)
    {
        if (g_slots[i].kind != Kind::Free) continue;
        g_slots[i].kind = kind;
        outHandle = packHandle (i, g_slots[i].gen);
        return i;
    }
    return -1;
}

// A slot whose handle was NEVER HANDED OUT goes back untouched. `freeSlot` exists to make an issued
// handle stale, and that costs a generation; doing it for a create that failed spends the table down
// for nothing.
void abandonSlot (Slot& s) noexcept
{
    s.master.reset();
    s.solution.reset();
    s.kind = Kind::Free;
}

void freeSlot (Slot& s) noexcept
{
    s.master.reset();
    s.solution.reset();
    s.kind = Kind::Free;
    // Bump the generation so the handle just destroyed can never address the next object here. 24 bits
    // wide, and it wraps rather than retiring the slot: 16.7 million destroys of ONE slot before a
    // number repeats is not a budget anybody meets, while retirement was one that a per-file render
    // loop meets in an afternoon.
    s.gen = (s.gen + 1) & 0x00FFFFFFu;
    if (s.gen == 0u) s.gen = 1u;
}

// `MasteringChain::process` takes planar pointers; the ABI carries one pointer and a stride. Building
// the table on the stack is the whole of the translation.
void planes (float* base, std::uint32_t stride, int nch, float** out) noexcept
{
    for (int c = 0; c < nch; ++c) out[c] = base + (std::size_t) c * (std::size_t) stride;
}

}   // namespace

//==============================================================================
// ENTRY POINTS

FC_EXPORT fc_status fc_master_create (const fc_master_config* cfg, fc_master* out)
{
    // `*out` IS NOT TOUCHED UNTIL THE CALL SUCCEEDS. It used to be cleared on entry, which reads as the
    // careful thing and is not: a caller reusing a variable that still holds a LIVE handle would have it
    // wiped by a call that failed on the version field, and the object it named would then be
    // unreachable and undestroyable. A refused call is indistinguishable from one never made — that
    // includes its arguments.
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    if (const fc_status st = checkHeaderIn (cfg); st != FC_OK) return st;

    MasteringChainConfig cc {};
    if (const fc_status st = toCore (*cfg, cc); st != FC_OK) return st;
    if (cfg->channels < 1 || cfg->channels > core::kMaxChannels) return FC_ERR_REFUSED_BY_CORE;

    std::uint32_t h = 0;
    const int idx = allocSlot (Kind::Master, h);
    if (idx < 0) return FC_ERR_EXHAUSTED;

    // `new` here, not in process(): the RT-safety claim is about the audio path, and this is a worker
    // call. Under -fno-exceptions a failure aborts rather than returning — stated in fc_master_abi.h as
    // a non-promise rather than pretended away.
    g_slots[idx].master = std::make_unique<MasterInstance>();
    auto& m = *g_slots[idx].master;
    m.cfg = cc;

    // The renderer's block size has NO effect on the result — the chain re-blocks everything to its own
    // quantum — so it is chosen once, here, for the solver's tap sizing and never exposed.
    if (! m.renderer.prepare (cfg->channels, 4096)
        || ! m.chain.prepare (cfg->sampleRate, cfg->channels, cc))
    {
        // `abandonSlot`, not `freeSlot`: no handle ever left this function, so there is nothing for a
        // stale one to alias and no reason to spend a generation. Spending one here is not free —
        // generations are 8 bits and a slot RETIRES when they run out, so 2040 refused creates would
        // exhaust a table that had never issued a single handle.
        abandonSlot (g_slots[idx]);
        return FC_ERR_REFUSED_BY_CORE;
    }
    *out = h;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_configure (fc_master h, const fc_master_params* params,
                                         fc_master_resolved* resolved)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkHeaderIn (params); st != FC_OK) return st;
    if (const fc_status st = checkHeaderOut (resolved); st != FC_OK) return st;

    auto& m = *s->master;
    // Re-preparing would silently discard a stream in progress. Refused, and `fc_master_reset` is how a
    // caller gets back to the head of one — see the long note in fc_master_abi.h for why this call has
    // to re-prepare at all.
    if (m.audioSeen) return FC_ERR_STATE;

    MasteringChainParams cp {};
    if (const fc_status st = toCore (*params, cp); st != FC_OK) return st;

    // NOTHING HAS MOVED UNTIL HERE. A refused call above left the chain exactly as it was, which is what
    // makes "a refused call is indistinguishable from one never made" true at this boundary too.
    m.chain.setParams (cp);
    if (! m.chain.prepare (m.chain.sampleRate(), m.chain.numChannels(), m.cfg))
        return FC_ERR_REFUSED_BY_CORE;

    m.framesIn = m.framesFlushed = 0;
    m.solverRan = false;                    // a known parameter set is back in the chain
    fromCore (m.chain.resolved(), m.chain.tapOversampleFactor(), *resolved);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_resolved_get (fc_master h, fc_master_resolved* out)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkHeaderOut (out); st != FC_OK) return st;
    fromCore (s->master->chain.resolved(), s->master->chain.tapOversampleFactor(), *out);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_process (fc_master h, const float* in, float* out, std::uint32_t frames)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();

    if (frames == 0)
    {
        // `n == 0` is the one true no-op in law 11: no time, no falling edge, nothing. It is not an
        // error, and it must not set `audioSeen` either — no audio was seen.
        return FC_OK;
    }
    if (m.solverRan) return FC_ERR_STATE;   // the chain holds the SOLVER's parameters — see MasterInstance
    if (frames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;   // the core takes `int`

    if (const fc_status st = checkAudio (in, frames, nch); st != FC_OK) return st;
    if (const fc_status st = checkAudio (out, frames, nch); st != FC_OK) return st;
    const std::uint64_t bytes = (std::uint64_t) frames * (std::uint64_t) nch * sizeof (float);
    if (partiallyOverlaps (in, out, bytes)) return FC_ERR_SPAN;

    if (in != out) std::memcpy (out, in, (std::size_t) bytes);        // transport, not arithmetic

    float* pl[core::kMaxChannels] {};
    planes (out, frames, nch, pl);
    if (! m.chain.process (pl, nch, (int) frames)) return FC_ERR_REFUSED_BY_CORE;

    m.framesIn += frames;
    m.audioSeen = true;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_flush (fc_master h, float* out, std::uint32_t capacity, std::uint32_t* written)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (written); st != FC_OK) return st;
    *written = 0;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();

    if (m.solverRan) return FC_ERR_STATE;   // as process(): the chain's configuration is not the caller's
    if (capacity == 0) return FC_ERR_CAPACITY;
    if (capacity > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    // A capacity below the latency cannot drain the tail, and the core keeps no arrears, so a second
    // call would not continue this drain — it would start another one. Refused as a whole (law 11).
    if ((std::int64_t) capacity < (std::int64_t) m.chain.latencySamples()) return FC_ERR_CAPACITY;
    if (const fc_status st = checkAudio (out, capacity, nch); st != FC_OK) return st;
    if (aliasesSpan (written, sizeof (*written), out,
                     (std::uint64_t) capacity * (std::uint64_t) nch * sizeof (float))) return FC_ERR_SPAN;

    float* pl[core::kMaxChannels] {};
    planes (out, capacity, nch, pl);
    const int n = m.chain.flush (pl, nch, (int) capacity);
    if (n <= 0) return FC_ERR_REFUSED_BY_CORE;
    *written = (std::uint32_t) n;
    m.framesFlushed += (std::uint32_t) n;
    m.audioSeen = true;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_latency (fc_master h, std::int32_t* out)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    *out = s->master->chain.latencySamples();
    return FC_OK;
}

FC_EXPORT fc_status fc_master_get_stats (fc_master h, fc_master_stats* out)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkHeaderOut (out); st != FC_OK) return st;
    out->framesIn      = s->master->framesIn;
    out->framesFlushed = s->master->framesFlushed;
    out->nonFiniteIn   = s->master->chain.nonFiniteInputSamples();
    return FC_OK;
}

FC_EXPORT fc_status fc_master_reset (fc_master h)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    s->master->chain.reset();
    s->master->framesIn = s->master->framesFlushed = 0;
    s->master->audioSeen = false;
    // `solverRan` is NOT cleared here, and that is the point: `reset()` clears audio state and leaves
    // the solver's gain and ceiling standing, so a reset-then-process would still render a parameter
    // set the caller never asked for. Only `configure` puts a known one back.
    return FC_OK;
}

FC_EXPORT fc_status fc_master_destroy (fc_master h)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    freeSlot (*s);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_measure_lra (fc_master h, const float* in, std::uint32_t frames, double* out)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();
    if (frames == 0) return FC_ERR_REFUSED_BY_CORE;
    if (frames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    if (const fc_status st = checkAudio (in, frames, nch); st != FC_OK) return st;

    if (! m.solver.isPrepared()
        && ! m.solver.prepare (m.chain.sampleRate(), nch, m.renderer.blockSize(),
                               m.chain.internalBlock(), m.chain.tapOversampleFactor()))
        return FC_ERR_REFUSED_BY_CORE;

    const float* pl[core::kMaxChannels] {};
    for (int c = 0; c < nch; ++c) pl[c] = in + (std::size_t) c * (std::size_t) frames;
    double lra = 0.0;
    // A measurement has no gate to hide behind, so this one really is a refusal: the core now says false
    // when its own meter flagged the programme, and a facade that published the number anyway would be
    // publishing a measurement it had been told not to trust.
    if (! m.solver.measureInputLoudnessRange (pl, nch, (int) frames, lra)) return FC_ERR_REFUSED_BY_CORE;
    *out = lra;
    return FC_OK;
}

// BS.1770 CHANNEL WEIGHTS, forwarded to the solver's meters. Not decoration and not a convenience: the
// standard weights Ls/Rs at 1.41 and EXCLUDES LFE, and the core says in as many words that its default
// of 1.0 everywhere is correct for mono and stereo and wrong for surround. This ABI accepts up to
// `kMaxChannels`, so without this entry point a correct surround search could not be expressed through
// it at all — which would make the facade a NARROWER road than the C++ API, and the thinness law is
// about both directions. The host-layout-to-role mapping stays outside, exactly as the core says.
FC_EXPORT fc_status fc_master_set_channel_weight (fc_master h, std::int32_t channel, double weight)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (channel < 0 || channel >= core::kMaxChannels) return FC_ERR_RANGE;
    if (! std::isfinite (weight) || weight < 0.0) return FC_ERR_NON_FINITE;
    auto& m = *s->master;
    if (! m.solver.isPrepared()
        && ! m.solver.prepare (m.chain.sampleRate(), m.chain.numChannels(), m.renderer.blockSize(),
                               m.chain.internalBlock(), m.chain.tapOversampleFactor()))
        return FC_ERR_REFUSED_BY_CORE;
    m.solver.setChannelWeight (channel, weight);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_solve (fc_master h, const fc_master_params* params,
                                     const fc_loudness_request* req,
                                     const float* in, float* out, std::uint32_t frames,
                                     fc_solution* out_solution)
{
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (out_solution); st != FC_OK) return st;
    *out_solution = 0;
    if (const fc_status st = checkHeaderIn (params); st != FC_OK) return st;
    if (const fc_status st = checkHeaderIn (req); st != FC_OK) return st;

    auto& m = *s->master;
    const int nch = m.chain.numChannels();
    // `frames == 0` is NOT refused here. The solver has its own answer for it — `InvalidRequest`, a
    // VERDICT — and this file's contract says FC_OK means a verdict was obtained and that it never
    // guesses a reason the core has. Refusing it here made the ABI answer for the core on one input and
    // forward the core on every other, which is two policies for one question.
    if (frames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    if (const fc_status st = checkAudio (in, frames, nch); st != FC_OK) return st;
    if (const fc_status st = checkAudio (out, frames, nch); st != FC_OK) return st;
    // A SEARCH cannot render in place: every pass after the first would read the previous pass's master
    // as its input. The core refuses that per channel; here the whole planar block is one span, so any
    // touching at all is refused, equality included.
    const std::uint64_t bytes = (std::uint64_t) frames * (std::uint64_t) nch * sizeof (float);
    if (in == out || partiallyOverlaps (in, out, bytes)) return FC_ERR_SPAN;
    if (aliasesSpan (out_solution, sizeof (*out_solution), out, bytes)) return FC_ERR_SPAN;

    MasteringChainParams cp {};
    if (const fc_status st = toCore (*params, cp); st != FC_OK) return st;
    LoudnessRequest lr {};
    if (const fc_status st = toCore (*req, lr); st != FC_OK) return st;

    if (! m.solver.isPrepared()
        && ! m.solver.prepare (m.chain.sampleRate(), nch, m.renderer.blockSize(),
                               m.chain.internalBlock(), m.chain.tapOversampleFactor()))
        return FC_ERR_REFUSED_BY_CORE;

    std::uint32_t sh = 0;
    const int idx = allocSlot (Kind::Solution, sh);
    if (idx < 0) return FC_ERR_EXHAUSTED;

    const float* ip[core::kMaxChannels] {};
    float*       op[core::kMaxChannels] {};
    for (int c = 0; c < nch; ++c)
    {
        ip[c] = in  + (std::size_t) c * (std::size_t) frames;
        op[c] = out + (std::size_t) c * (std::size_t) frames;
    }

    g_slots[idx].solution = std::make_unique<LoudnessSolution> (
        m.solver.solve (m.chain, m.renderer, cp, ip, op, nch, (int) frames, lr));

    // The search drives the renderer, which RESETS the chain on every pass — so where it ran, the
    // handle's streaming state is gone and its counters would be lying if they survived.
    //
    // WHERE IT DID NOT RUN, NOTHING MAY MOVE, and this used to be unconditional. `NotPrepared` and
    // `InvalidRequest` are returned before `OfflineRenderer::render()` is reached, so the chain still
    // holds whatever was in its FIFO — and clearing `audioSeen` there disarmed the guard that stops
    // `configure` from re-preparing over a live stream. Measured: `process(64)` (configure correctly
    // refused with FC_ERR_STATE), then a solve with no target, and the very same configure was then
    // ACCEPTED and destroyed the 64 frames. A refused call has to be indistinguishable from one never
    // made even when what refused it was the core.
    const MasteringSolveStatus verdict = g_slots[idx].solution->status;
    if (verdict != MasteringSolveStatus::NotPrepared && verdict != MasteringSolveStatus::InvalidRequest)
    {
        m.framesIn = m.framesFlushed = 0;
        m.audioSeen = false;
        m.solverRan = true;         // and now process/flush refuse until configure — see MasterInstance
    }

    *out_solution = sh;
    return FC_OK;   // A VERDICT WAS OBTAINED. Whether it is `Solved` is the summary's business.
}

FC_EXPORT fc_status fc_solution_summary_get (fc_solution sh, fc_solution_summary* out)
{
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkHeaderOut (out); st != FC_OK) return st;
    const LoudnessSolution& v = *s->solution;
    out->status              = (std::int32_t) v.status;
    out->binding             = (std::int32_t) v.binding;
    out->alsoViolated        = v.alsoViolated;
    out->preLimiterGainDb    = v.preLimiterGainDb;
    out->ceilingDbTp         = v.ceilingDbTp;
    out->passes              = v.passes;
    out->logCount            = v.logCount;
    out->activityThresholdDb = v.activityThresholdDb;
    out->achievedBelowLufs   = v.achievedBelowLufs;
    out->achievedAboveLufs   = v.achievedAboveLufs;
    out->gainBelowDb         = v.gainBelowDb;
    out->gainAboveDb         = v.gainAboveDb;
    return FC_OK;
}

FC_EXPORT fc_status fc_solution_measurement (fc_solution sh, fc_measurement* out)
{
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkHeaderOut (out); st != FC_OK) return st;
    fromCore (s->solution->measured, *out);
    return FC_OK;
}

FC_EXPORT fc_status fc_solution_log (fc_solution sh, fc_solve_pass* out, std::uint32_t cap,
                                     std::uint32_t* written)
{
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (written); st != FC_OK) return st;
    *written = 0;
    if (cap == 0) return FC_OK;                       // asking for nothing is not an error
    if (out == nullptr) return FC_ERR_NULL;
    if ((reinterpret_cast<std::uintptr_t> (out) & 0x7u) != 0) return FC_ERR_ALIGNMENT;
    if (! inHeap (out, (std::uint64_t) cap * sizeof (fc_solve_pass))) return FC_ERR_SPAN;

    const LoudnessSolution& v = *s->solution;
    const std::uint32_t n = (std::uint32_t) v.logCount < cap ? (std::uint32_t) v.logCount : cap;
    for (std::uint32_t i = 0; i < n; ++i)
    {
        const SolvePassRecord& r = v.log[i];
        out[i].gainDb          = r.gainDb;
        out[i].ceilingDb       = r.ceilingDb;
        out[i].integratedLufs  = r.integratedLufs;
        out[i].truePeakDbTp    = r.truePeakDbTp;
        out[i].plrDb           = r.plrDb;
        out[i].limiterMaxGrDb  = r.limiterMaxGrDb;
        out[i].loudnessRangeLu = r.loudnessRangeLu;
        out[i].violated        = r.violated;
    }
    *written = n;
    return FC_OK;
}

FC_EXPORT fc_status fc_solution_destroy (fc_solution sh)
{
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    freeSlot (*s);
    return FC_OK;
}

//==============================================================================
// DEFAULTS — the core's own, written through the same mapping every other value crosses by.

FC_EXPORT void fc_master_config_default (fc_master_config* out)
{
    if (out == nullptr) return;
    const MasteringChainConfig d {};
    std::memset (out, 0, sizeof (*out));
    stampHeader (out);
    // LEFT AT ZERO, DELIBERATELY. There is no core default for either, and writing one here would be
    // this file choosing a geometry for every caller who forgot to — the same objection that keeps
    // `targetLufs` and `maxTruePeakDbTp` at NaN in the request. `fc_master_create` refuses both, so a
    // caller who forgets is told rather than silently given 48 kHz stereo.
    out->sampleRate            = 0.0;
    out->channels              = 0;
    out->internalBlock         = d.internalBlock;
    out->eq                    = d.eq ? 1 : 0;
    out->monoBass              = d.monoBass ? 1 : 0;
    out->compressor            = d.compressor ? 1 : 0;
    out->clipper               = d.clipper ? 1 : 0;
    out->limiter               = d.limiter ? 1 : 0;
    out->dither                = d.dither ? 1 : 0;
    out->compressorLookaheadMs = d.compressorLookaheadMs;
    out->limiterLookaheadMs    = d.limiterLookaheadMs;
    out->oversampleFactor      = d.oversampleFactor;
    out->tapsPerPhase          = d.tapsPerPhase;
    out->sidechainHpfHz        = d.sidechainHpfHz;
}

FC_EXPORT void fc_master_params_default (fc_master_params* out)
{
    if (out == nullptr) return;
    const MasteringChainParams d {};
    std::memset (out, 0, sizeof (*out));
    stampHeader (out);
    out->inputGainDb      = d.inputGainDb;
    out->preLimiterGainDb = d.preLimiterGainDb;

    for (int b = 0; b < FC_MAX_EQ_BANDS; ++b)
    {
        const eq::BandParams& sb = d.eqBands[b];
        fc_eq_band& db = out->eqBands[b];
        db.on     = sb.on ? 1 : 0;
        db.type   = (std::int32_t) sb.type;
        db.swept  = sb.swept ? 1 : 0;
        db.bypass = sb.bypass ? 1 : 0;
        db.dyn.on      = sb.dyn.on ? 1 : 0;
        db.dyn.rangeDb = sb.dyn.rangeDb;
        db.dyn.thrDb   = sb.dyn.thrDb;
        db.dyn.thrAuto = sb.dyn.thrAuto ? 1 : 0;
        db.dyn.atk     = sb.dyn.atk;
        db.dyn.rel     = sb.dyn.rel;
        for (int l = 0; l < FC_MAX_EQ_LANES; ++l)
        {
            const eq::LaneParams& sl = sb.lanes[l];
            fc_eq_lane& dl = db.lanes[l];
            dl.on     = sl.on ? 1 : 0;
            dl.freq   = sl.freq;
            dl.q      = sl.Q;
            dl.gainDb = sl.gainDb;
            dl.slope  = sl.slope;
            dl.bypass = sl.bypass ? 1 : 0;
        }
    }

    out->monoBass.enabled     = d.monoBass.enabled ? 1 : 0;
    out->monoBass.frequencyHz = d.monoBass.frequencyHz;
    out->monoBass.lowWidth    = d.monoBass.lowWidth;

    out->compressor.detector    = (std::int32_t) d.compressor.detector;
    out->compressor.link        = (std::int32_t) d.compressor.link;
    out->compressor.rmsWindowMs = d.compressor.rmsWindowMs;
    out->compressor.mode        = (std::int32_t) d.compressor.mode;
    out->compressor.thresholdDb = d.compressor.thresholdDb;
    out->compressor.ratio       = d.compressor.ratio;
    out->compressor.kneeDb      = d.compressor.kneeDb;
    out->compressor.rangeDb     = d.compressor.rangeDb;
    out->compressor.attackMs    = d.compressor.attackMs;
    out->compressor.releaseMs   = d.compressor.releaseMs;
    out->compressor.makeupDb    = d.compressor.makeupDb;
    out->compressor.autoMakeup  = d.compressor.autoMakeup ? 1 : 0;

    out->clipper.shape     = (std::int32_t) d.clipper.shape;
    out->clipper.driveDb   = d.clipper.driveDb;
    out->clipper.bias      = d.clipper.bias;
    out->clipper.mix       = d.clipper.mix;
    out->clipper.outputDb  = d.clipper.outputDb;
    out->clipper.autoComp  = d.clipper.autoComp;
    out->clipper.dcBlockHz = d.clipper.dcBlockHz;

    out->limiter.ceilingDbTp = d.limiter.ceilingDbTp;
    out->limiter.releaseMs   = d.limiter.releaseMs;

    out->dither.bits             = d.dither.bits;
    out->dither.shaping          = (std::int32_t) d.dither.shaping;
    out->dither.seedLo           = (std::uint32_t) (d.dither.seed & 0xFFFFFFFFull);
    out->dither.seedHi           = (std::uint32_t) (d.dither.seed >> 32);
    out->dither.autoBlank        = d.dither.autoBlank ? 1 : 0;
    out->dither.autoBlankSamples = d.dither.autoBlankSamples;

    out->bypassEq         = d.bypassEq ? 1 : 0;
    out->bypassMonoBass   = d.bypassMonoBass ? 1 : 0;
    out->bypassCompressor = d.bypassCompressor ? 1 : 0;
    out->bypassClipper    = d.bypassClipper ? 1 : 0;
    out->bypassLimiter    = d.bypassLimiter ? 1 : 0;
    out->bypassDither     = d.bypassDither ? 1 : 0;
}

FC_EXPORT void fc_loudness_request_default (fc_loudness_request* out)
{
    if (out == nullptr) return;
    const LoudnessRequest d {};
    std::memset (out, 0, sizeof (*out));
    stampHeader (out);
    // The two REQUIRED fields stay NaN: the core ships no default target because "-14 LUFS, -1 dBTP" is
    // a delivery policy and the core is product-neutral by rule. Carrying them here would be this file
    // choosing that policy for every caller who forgot to.
    out->targetLufs             = d.targetLufs;
    out->maxTruePeakDbTp        = d.maxTruePeakDbTp;
    out->toleranceLu            = d.toleranceLu;
    out->truePeakAimDb          = d.truePeakAimDb;
    out->limiterGr.limitDb      = d.limiterGr.limitDb;
    out->limiterGr.statistic    = (std::int32_t) d.limiterGr.statistic;
    out->compressorGr.limitDb   = d.compressorGr.limitDb;
    out->compressorGr.statistic = (std::int32_t) d.compressorGr.statistic;
    out->minPlrDb               = d.minPlrDb;
    out->maxLraLossLu           = d.maxLraLossLu;
    out->inputLoudnessRangeLu   = d.inputLoudnessRangeLu;
    out->activityThresholdDb    = d.activityThresholdDb;
    out->maxPasses              = d.maxPasses;
    out->initialGainDb          = d.initialGainDb;
}

//==============================================================================
FC_EXPORT std::uint32_t fc_master_abi_version  (void) { return FC_MASTER_ABI_VERSION; }
FC_EXPORT std::uint32_t fc_master_max_channels (void) { return (std::uint32_t) core::kMaxChannels; }
FC_EXPORT std::uint32_t fc_master_max_eq_bands (void) { return (std::uint32_t) FC_MAX_EQ_BANDS; }
FC_EXPORT std::uint32_t fc_master_sizeof_params(void) { return (std::uint32_t) sizeof (fc_master_params); }
FC_EXPORT std::uint32_t fc_master_sizeof_config(void) { return (std::uint32_t) sizeof (fc_master_config); }
