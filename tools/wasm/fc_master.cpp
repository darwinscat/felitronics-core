// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fc_master — the implementation of the C ABI declared in tools/fc_master_abi.h. Read that header first:
// it carries the contract, and this file carries only the unpacking of it.
//
// NOTHING HERE COMPUTES ANYTHING. Every number handed back was produced by `felitronics::mastering` and
// read out of it — save `fc_need.facadeBytes`, this file's own `sizeof`, forwarded beside the core's numbers
// and never added to them. What this file does is five things and no sixth: validate what a page can hand
// it, translate enum codes, copy fields, own handles, and refuse every call once one has not returned (the
// poison — see FC_GUARD). If a line of this file ever looks like DSP, it is
// a defect, and `fcore_master selftest` is the gate that finds it — the same programme through this
// ABI and through a direct C++ call, one binary, one machine, compared bit for bit.
//
// It compiles natively as well as under emscripten: EMSCRIPTEN_KEEPALIVE degrades to a plain extern "C",
// so the whole validation and addressing layer runs under ctest, ASan and UBSan like anything else.
// That is the same arrangement fc_probe.cpp uses and for the same reason.

#include "fc_master_abi.h"

#include <felitronics/mastering/DeliveredMastering.h>
#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
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
static_assert (sizeof (MasteringChainParams)           == 6552);   // + compressorMix (P60) — see below
static_assert (sizeof (MasteringChainResolved)         == 80);     // + compressorMix (P60) — see below

// THE PIN THAT WORKS THROUGH INHERITANCE. A structured binding cannot decompose a type whose base has
// members, and `sizeof` is blind to a field that lands in existing padding — so between them the two
// mechanisms above leave the whole `CompressorParams` chain unguarded. Aggregate initialisation is not
// blind to either: brace-initialising with EXACTLY the expected list of member types compiles only
// while that list is the whole of the type. The negative control is half the pin: without it, a
// too-short list would still compile through value-initialisation of the rest.
template <class T, class... A>
concept BraceInit = requires { T { A {}... }; };

static_assert (BraceInit<dynamics::GainReductionParams,
                         dynamics::DetectorParams, dynamics::Mode,
                         double, double, double, double, double, double>,
               "GainReductionParams grew or lost a member — update the ABI mapping with it");
static_assert (! BraceInit<dynamics::GainReductionParams,
                           dynamics::DetectorParams, dynamics::Mode,
                           double, double, double, double, double, double, double>,
               "negative control: the pin above must be exact, not a lower bound");

static_assert (BraceInit<dynamics::CompressorParams,
                         dynamics::GainReductionParams, double, bool, double>,
               "CompressorParams grew or lost a member — update the ABI mapping with it");
static_assert (! BraceInit<dynamics::CompressorParams,
                           dynamics::GainReductionParams, double, bool, double, double>,
               "negative control");

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

    // `CompressorParams` inherits, so it cannot be decomposed, and NEITHER CAN ITS BASE: the chain is
    // `CompressorParams : GainReductionParams : DetectorParams`, and only the great-grandparent
    // decomposes. An earlier version of this comment claimed the size pins covered the rest and that
    // was FALSE, measured: a `bool` dropped into the padding after `mode` leaves
    // `sizeof(GainReductionParams)` at 72, and one after `autoMakeup` leaves `sizeof(CompressorParams)`
    // at 96 — both pins blind on exactly the type the arity check exists for, and nine mapped fields
    // resting on them. The BraceInit pins below are what actually catch it.
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

    // `p_mix` and `r_mix` arrived with P60 as a build break here, and are mapped since ABI v3 — `toCore` and
    // `fromCore` below, at the end of `fc_master_params` and `fc_master_resolved`.
    MasteringChainParams prm {};
    auto& [p_in, p_pre, p_eq, p_mb, p_comp, p_clip, p_lim, p_dith,
           p_bE, p_bM, p_bC, p_bK, p_bL, p_bD, p_mix] = prm;
    (void) p_in; (void) p_pre; (void) p_eq; (void) p_mb; (void) p_comp; (void) p_clip; (void) p_lim;
    (void) p_dith; (void) p_bE; (void) p_bM; (void) p_bC; (void) p_bK; (void) p_bL; (void) p_bD; (void) p_mix;

    MasteringChainResolved res {};
    auto& [r_lat, r_blk, r_clook, r_clip, r_lim, r_llook, r_os, r_ctap, r_ltap,
           r_ceil, r_rel, r_mb, r_mix] = res;
    (void) r_lat; (void) r_blk; (void) r_clook; (void) r_clip; (void) r_lim; (void) r_llook;
    (void) r_os; (void) r_ctap; (void) r_ltap; (void) r_ceil; (void) r_rel; (void) r_mb; (void) r_mix;
}

//==============================================================================
// THE ABI'S OWN LAYOUT — the pins that make VERSIONING's rules mechanical rather than a promise.
//
// The size table is the single source (rule 5); `sizeFor` reads it, and everything that checks a caller's
// `structSize` goes through `sizeFor` with the caller's version. The pins below tie the table to this build's
// structs: the newest row of every struct is its `sizeof`, a struct's rows only ever grow, no struct with a header
// ends in implicit padding (rule 4), and every top-level field sits at the offset the version that added it published.
// What these offset and size pins make a build error is a change that MOVES something they read: a top-level field
// moved, a field inserted in front of others, the last field retyped (FC_ENDS_AT takes its size from the compiler), a
// retype that shifts the field after it. A change that moves nothing builds: a field dropped into padding (an `int32_t`
// after `fc_loudness_request::maxPasses`), or a retype that keeps every size and offset — `int32_t` to `uint32_t`, a
// `double` narrowed to `float` in front of another `double` — except in v4's fields, which are pinned by TYPE as well
// (beside the frozen sizes below); v1..v3's are not. Inside a frozen nested struct only the size is pinned, so two of its
// fields swapped build too; tools/wasm/layout-check.mjs sees that one under ctest, by the offset of every field
// `fcore_master layout` lists — it compares no types.
template <typename T> struct AbiId;   // declared, never defined: an unmapped struct fails to compile
template <> struct AbiId<fc_master_config>    { static constexpr int id = FC_STRUCT_CONFIG; };
template <> struct AbiId<fc_master_params>    { static constexpr int id = FC_STRUCT_PARAMS; };
template <> struct AbiId<fc_master_resolved>  { static constexpr int id = FC_STRUCT_RESOLVED; };
template <> struct AbiId<fc_master_stats>     { static constexpr int id = FC_STRUCT_STATS; };
template <> struct AbiId<fc_need>             { static constexpr int id = FC_STRUCT_NEED; };
template <> struct AbiId<fc_loudness_request> { static constexpr int id = FC_STRUCT_REQUEST; };
template <> struct AbiId<fc_measurement>      { static constexpr int id = FC_STRUCT_MEASUREMENT; };
template <> struct AbiId<fc_solution_summary> { static constexpr int id = FC_STRUCT_SUMMARY; };

constexpr std::uint32_t sizeFor (int id, std::uint32_t version) noexcept
{
    if (version < 1u || version > FC_MASTER_ABI_VERSION) return 0u;
    std::uint32_t since = 0u, bytes = 0u;
#define FC_ROW(sid, v, b) \
    if ((int) (sid) == id && (std::uint32_t) (v) <= version && (std::uint32_t) (v) > since) { since = (v); bytes = (b); }
    FC_MASTER_STRUCT_SIZES (FC_ROW)
#undef FC_ROW
    return bytes;
}

// Rows of one struct, in version order, strictly grow — and no two rows claim the same version. The ids walked are
// the table's own, not a hand-kept range, and a struct may be ABSENT (0) before the version that introduced it
// (rule 5) — but once it exists it never shrinks and never vanishes.
constexpr int tableMaxId() noexcept
{
    int m = -1;
#define FC_MAXID(sid, v, b) if ((int) (sid) > m) m = (int) (sid); (void) (v); (void) (b);
    FC_MASTER_STRUCT_SIZES (FC_MAXID)
#undef FC_MAXID
    return m;
}
constexpr bool tableGrows() noexcept
{
    for (int id = 0; id <= tableMaxId(); ++id)
    {
        std::uint32_t prev = 0u;
        for (std::uint32_t v = 1u; v <= FC_MASTER_ABI_VERSION; ++v)
        {
            const std::uint32_t s = sizeFor (id, v);
            if (prev != 0u && s < prev) return false;      // shrank, or vanished after it existed
            prev = s;
        }
        if (prev == 0u) return false;                      // an id with no size at the current version
        int rows = 0;
#define FC_COUNT(sid, v, b) if ((int) (sid) == id) { ++rows; (void) (v); (void) (b); }
        FC_MASTER_STRUCT_SIZES (FC_COUNT)
#undef FC_COUNT
        // A struct whose size did not change between two rows would be a row that says nothing — or a version
        // bump nobody can detect by size (rule 4). Counting the distinct sizes against the rows catches both.
        int distinct = 0; std::uint32_t last = 0u;
        for (std::uint32_t v = 1u; v <= FC_MASTER_ABI_VERSION; ++v)
            if (const std::uint32_t s = sizeFor (id, v); s != last) { ++distinct; last = s; }
        if (distinct != rows) return false;
    }
    return true;
}
static_assert (tableGrows(), "FC_MASTER_STRUCT_SIZES: a struct's rows must strictly grow, one row per size");

template <typename T>
constexpr bool newestRowIsSizeof() noexcept { return sizeFor (AbiId<T>::id, FC_MASTER_ABI_VERSION) == sizeof (T); }
static_assert (newestRowIsSizeof<fc_master_config>());
static_assert (newestRowIsSizeof<fc_master_params>());
static_assert (newestRowIsSizeof<fc_master_resolved>());
static_assert (newestRowIsSizeof<fc_master_stats>());
static_assert (newestRowIsSizeof<fc_need>());
static_assert (newestRowIsSizeof<fc_loudness_request>());
static_assert (newestRowIsSizeof<fc_measurement>());
static_assert (newestRowIsSizeof<fc_solution_summary>());

// Rule 4: no struct with a header ends in padding the layout does not name.
// The size of the last field is the compiler's (`sizeof (T::last)`), not a literal: a `double` retyped to `float`
// leaves the struct's size and the field's offset where they were and opens four bytes of padding — a literal 8 here
// would have agreed with both.
#define FC_ENDS_AT(T, last) \
    static_assert (sizeof (T) == offsetof (T, last) + sizeof (T::last), #T " ends in implicit padding — name it (rule 4)")
FC_ENDS_AT (fc_master_config,    deliveryRate);
FC_ENDS_AT (fc_master_params,    compressorMix);
FC_ENDS_AT (fc_master_resolved,  compressorMix);
FC_ENDS_AT (fc_master_stats,     nonFiniteIn);
FC_ENDS_AT (fc_need,             _pad0);
FC_ENDS_AT (fc_loudness_request, initialGainDb);
FC_ENDS_AT (fc_measurement,      limiterGrTraceValid);
FC_ENDS_AT (fc_solution_summary, gainAboveDb);
// and the types the table's sizes were computed from
static_assert (sizeof (fc_master_config::deliveryRate) == 8 && sizeof (fc_master_params::compressorMix) == 8
               && sizeof (fc_master_resolved::compressorMix) == 8);
#undef FC_ENDS_AT

// Rule 3: what never grows. A size moved here is a frozen struct that grew.
static_assert (sizeof (fc_header) == 8);
static_assert (sizeof (fc_eq_lane) == 40 && sizeof (fc_eq_dyn) == 48 && sizeof (fc_eq_band) == 264);
static_assert (sizeof (fc_mono_bass) == 12 && sizeof (fc_compressor) == 88 && sizeof (fc_clipper) == 28);
static_assert (sizeof (fc_limiter) == 16 && sizeof (fc_dither) == 24 && sizeof (fc_gr_limit) == 16);
static_assert (sizeof (fc_solve_pass) == 64 && sizeof (fc_gr_stats) == 64);
static_assert (sizeof (fc_gr_trace_bucket) == 24);                                                          // v4, frozen
// v4's fields by TYPE as well as by offset: a `double` retyped to `float` keeps this struct's size and every offset (the
// four bytes become padding), so neither pin above would see it, and JavaScript would read eight bytes where four were
// written (the code-review round, by mutation).
static_assert (std::is_same_v<decltype (fc_gr_trace_bucket::maxDb), double> && std::is_same_v<decltype (fc_gr_trace_bucket::meanDb), double>
               && std::is_same_v<decltype (fc_gr_trace_bucket::samples), uint32_t> && std::is_same_v<decltype (fc_gr_trace_bucket::nonFinite), uint32_t>);
static_assert (std::is_same_v<decltype (fc_measurement::compressorGrTraceBuckets), int32_t> && std::is_same_v<decltype (fc_measurement::limiterGrTraceBuckets), int32_t>
               && std::is_same_v<decltype (fc_measurement::compressorGrTraceValid), int32_t> && std::is_same_v<decltype (fc_measurement::limiterGrTraceValid), int32_t>);

// Every top-level field at the offset it was published at. A field inserted in front of others moves a number; one dropped
// into padding moves none (see the top of this section).
#define FC_AT(T, f, off) static_assert (offsetof (T, f) == (off), #T "::" #f " moved")
FC_AT (fc_master_config, header, 0);           FC_AT (fc_master_config, sampleRate, 8);
FC_AT (fc_master_config, channels, 16);        FC_AT (fc_master_config, internalBlock, 20);
FC_AT (fc_master_config, eq, 24);              FC_AT (fc_master_config, monoBass, 28);
FC_AT (fc_master_config, compressor, 32);      FC_AT (fc_master_config, clipper, 36);
FC_AT (fc_master_config, limiter, 40);         FC_AT (fc_master_config, dither, 44);
FC_AT (fc_master_config, compressorLookaheadMs, 48); FC_AT (fc_master_config, limiterLookaheadMs, 56);
FC_AT (fc_master_config, oversampleFactor, 64); FC_AT (fc_master_config, tapsPerPhase, 68);
FC_AT (fc_master_config, sidechainHpfHz, 72);  FC_AT (fc_master_config, deliveryRate, 80);        // v2

FC_AT (fc_master_params, header, 0);           FC_AT (fc_master_params, inputGainDb, 8);
FC_AT (fc_master_params, preLimiterGainDb, 16); FC_AT (fc_master_params, eqBands, 24);
FC_AT (fc_master_params, monoBass, 6360);      FC_AT (fc_master_params, compressor, 6376);
FC_AT (fc_master_params, clipper, 6464);       FC_AT (fc_master_params, limiter, 6496);
FC_AT (fc_master_params, dither, 6512);        FC_AT (fc_master_params, bypassEq, 6536);
FC_AT (fc_master_params, bypassMonoBass, 6540); FC_AT (fc_master_params, bypassCompressor, 6544);
FC_AT (fc_master_params, bypassClipper, 6548); FC_AT (fc_master_params, bypassLimiter, 6552);
FC_AT (fc_master_params, bypassDither, 6556); FC_AT (fc_master_params, compressorMix, 6560);             // v3

FC_AT (fc_master_resolved, header, 0);         FC_AT (fc_master_resolved, latencySamples, 8);
FC_AT (fc_master_resolved, internalBlock, 12); FC_AT (fc_master_resolved, compressorLookahead, 16);
FC_AT (fc_master_resolved, clipperLatency, 20); FC_AT (fc_master_resolved, limiterLatency, 24);
FC_AT (fc_master_resolved, limiterLookahead, 28); FC_AT (fc_master_resolved, oversampleFactor, 32);
FC_AT (fc_master_resolved, compressorTapOffset, 36); FC_AT (fc_master_resolved, limiterTapOffset, 40);
FC_AT (fc_master_resolved, limiterCeilingDbTp, 48); FC_AT (fc_master_resolved, limiterReleaseMs, 56);
FC_AT (fc_master_resolved, monoBass, 64);      FC_AT (fc_master_resolved, tapOversampleFactor, 76);
FC_AT (fc_master_resolved, compressorMix, 80);                                                              // v3

FC_AT (fc_master_stats, header, 0);            FC_AT (fc_master_stats, framesIn, 8);
FC_AT (fc_master_stats, framesFlushed, 16);    FC_AT (fc_master_stats, nonFiniteIn, 24);

FC_AT (fc_need, header, 0);                    FC_AT (fc_need, callBytes, 8);
FC_AT (fc_need, solverPrepareBytes, 16);       FC_AT (fc_need, facadeBytes, 24);
FC_AT (fc_need, solverPrepared, 32);           FC_AT (fc_need, _pad0, 36);

FC_AT (fc_loudness_request, header, 0);        FC_AT (fc_loudness_request, targetLufs, 8);
FC_AT (fc_loudness_request, toleranceLu, 16);  FC_AT (fc_loudness_request, maxTruePeakDbTp, 24);
FC_AT (fc_loudness_request, truePeakAimDb, 32); FC_AT (fc_loudness_request, limiterGr, 40);
FC_AT (fc_loudness_request, compressorGr, 56); FC_AT (fc_loudness_request, minPlrDb, 72);
FC_AT (fc_loudness_request, maxLraLossLu, 80); FC_AT (fc_loudness_request, inputLoudnessRangeLu, 88);
FC_AT (fc_loudness_request, activityThresholdDb, 96); FC_AT (fc_loudness_request, maxPasses, 104);
FC_AT (fc_loudness_request, initialGainDb, 112);

FC_AT (fc_measurement, header, 0);             FC_AT (fc_measurement, integratedLufs, 8);
FC_AT (fc_measurement, truePeakDbTp, 16);      FC_AT (fc_measurement, samplePeakDb, 24);
FC_AT (fc_measurement, loudnessRangeLu, 32);   FC_AT (fc_measurement, plrDb, 40);
FC_AT (fc_measurement, compressor, 48);        FC_AT (fc_measurement, limiter, 112);
FC_AT (fc_measurement, limiterMaxReconstructedPeakDb, 176); FC_AT (fc_measurement, latencySamples, 184);
FC_AT (fc_measurement, gatingBlocks, 188);     FC_AT (fc_measurement, droppedBlocks, 192);
FC_AT (fc_measurement, nonFiniteSubHops, 196); FC_AT (fc_measurement, loudnessValid, 200);
FC_AT (fc_measurement, lraValid, 204);
FC_AT (fc_measurement, compressorGrTraceBuckets, 208); FC_AT (fc_measurement, limiterGrTraceBuckets, 212);   // v4
FC_AT (fc_measurement, compressorGrTraceValid, 216);   FC_AT (fc_measurement, limiterGrTraceValid, 220);     // v4
FC_AT (fc_gr_trace_bucket, maxDb, 0);          FC_AT (fc_gr_trace_bucket, meanDb, 8);                        // v4
FC_AT (fc_gr_trace_bucket, samples, 16);       FC_AT (fc_gr_trace_bucket, nonFinite, 20);                    // v4

FC_AT (fc_solution_summary, header, 0);        FC_AT (fc_solution_summary, status, 8);
FC_AT (fc_solution_summary, binding, 12);      FC_AT (fc_solution_summary, alsoViolated, 16);
FC_AT (fc_solution_summary, preLimiterGainDb, 24); FC_AT (fc_solution_summary, ceilingDbTp, 32);
FC_AT (fc_solution_summary, passes, 40);       FC_AT (fc_solution_summary, logCount, 44);
FC_AT (fc_solution_summary, activityThresholdDb, 48); FC_AT (fc_solution_summary, achievedBelowLufs, 56);
FC_AT (fc_solution_summary, achievedAboveLufs, 64); FC_AT (fc_solution_summary, gainBelowDb, 72);
FC_AT (fc_solution_summary, gainAboveDb, 80);
#undef FC_AT

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

// Rules 6 and 7 of VERSIONING. `bytes` comes back as the size of the CALLER's version, and it is the number every
// later check on that struct uses — a bound or an alias test taken with this build's `sizeof` would refuse a
// legal call whose older, shorter struct sits right against the next thing in the heap.
template <typename T>
fc_status checkHeader (const T* p, std::uint32_t& bytes) noexcept
{
    if (p == nullptr) return FC_ERR_NULL;
    if ((reinterpret_cast<std::uintptr_t> (p) & 0x7u) != 0) return FC_ERR_ALIGNMENT;
    // THE HEADER FIRST, THEN THE REST OF THE SPAN, and the two steps are separate on purpose. Bounding
    // the whole struct before reading the version means a caller who got the version wrong is told
    // FC_ERR_SPAN — the wrong diagnosis, on the one field whose job is to catch exactly that mistake.
    // Eight bytes is what the version costs to read, so eight bytes is what is bounded first.
    if (! inHeap (p, sizeof (fc_header))) return FC_ERR_SPAN;
    const std::uint32_t size = sizeFor (AbiId<T>::id, p->header.abiVersion);
    if (size == 0u) return FC_ERR_ABI_VERSION;              // older than v1, or newer than this build
    if (p->header.structSize != size) return FC_ERR_STRUCT_SIZE;
    if (! inHeap (p, size)) return FC_ERR_SPAN;
    bytes = size;
    return FC_OK;
}

// The defaults writers, at this build's layout — defined with the entry points that publish them.
void writeDefaults (fc_master_config& out) noexcept;
void writeDefaults (fc_master_params& out) noexcept;
void writeDefaults (fc_loudness_request& out) noexcept;

// An IN struct as THIS build lays it out: its defaults first, then exactly the caller's bytes over them. A v1
// config read this way has `deliveryRate == 0` — the default, which is v1 — and never the bytes that follow the
// caller's 80, whatever they hold.
template <typename T>
T loadIn (const T* p, std::uint32_t bytes) noexcept
{
    T local {};
    writeDefaults (local);
    std::memcpy (&local, p, bytes);
    return local;
}

// The OUT half. The caller's header stays as it is — an echo — and exactly `bytes` are written, header excluded.
template <typename T>
void writeOut (T* p, const T& local, std::uint32_t bytes) noexcept
{
    std::memcpy (static_cast<void*> ((unsigned char*) p + sizeof (fc_header)),
                 (const unsigned char*) &local + sizeof (fc_header), bytes - sizeof (fc_header));
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
        || ! fin (c.limiterLookaheadMs) || ! fin (c.sidechainHpfHz) || ! fin (c.deliveryRate)) return FC_ERR_NON_FINITE;
    // `deliveryRate` has no core field to land in, and that is the design: it is not chain topology (a chain
    // does not emit a variable frame count) but the choice of which core object stands in front of the chain.
    // `fc_master_create` reads it from the struct.
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
    // v3. Finite is this file's check; the RANGE is the core's, which clamps to [0, 1] and reports the applied
    // value in `resolved` — and would map a NaN to 1 without saying so, which is why a NaN stops here.
    if (! fin (p.compressorMix)) return FC_ERR_NON_FINITE;
    out.compressorMix    = p.compressorMix;
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
    out.compressorMix       = r.compressorMix;          // v3
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
    // v2: a DELIVERING handle (`fc_master_config::deliveryRate != 0`). The chain above then runs at the delivery
    // rate and this stands in front of it; unprepared and empty on a handle that does not deliver.
    DeliveredMastering   delivered;
    bool                 delivering = false;
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

// The renderer's block size has NO effect on the result — the chain re-blocks everything to its own
// quantum — so it is chosen once, here, for the solver's tap sizing and never exposed. Named rather than
// spelled twice now that `fc_master_need_create` has to budget the very renderer `fc_master_create` builds.
constexpr int kRendererBlock = 4096;

// The verdict and the budget of a create — ONE core expression per kind of handle, chosen by the field that
// decides the kind. Shared by `fc_master_create` and `fc_master_need_create`, so the two cannot disagree.
std::uint64_t createBytesFor (const fc_master_config& c, const MasteringChainConfig& cc) noexcept
{
    return c.deliveryRate != 0.0
        ? DeliveredMastering::createBytes (c.sampleRate, c.deliveryRate, c.channels, cc, kRendererBlock)
        : mastering::createBytes (c.sampleRate, c.channels, cc, kRendererBlock);
}

// THE CALL THAT NEVER RETURNED. Under -fno-exceptions an exhausted heap ABORTS inside a core call, and the module is
// not stopped by it: emscripten lets the page call again, the heap is intact, and every object is wherever the abort
// left it. Measured on this facade in wasm32 (P41): an abort inside `fc_master_solve`, then `fc_master_process` on the
// same handle answered FC_OK and rendered at the solver's own pass-1 gain, because `solverRan` is written after the
// solve returns — +12 dB in that replay, which asked for that start, and MasterAbiTests repeats it. And
// kMaxHandles − 1 = 7 such solves left the handle table full for good (FC_ERR_EXHAUSTED with one live object, whether
// the memory came back or not). An instance in that state does not fail; it LIES.
//
// So it is POISONED instead. Every status-returning entry point marks a call in progress on entry and clears the mark
// only on a normal return. Finding the mark already set on entry means the previous call never returned — an abort, a
// trap, or natively an exception that escaped — and from then on every such call answers FC_ERR_POISONED and touches
// nothing. The way out is a new module instance. The mark survives an abort because nothing unwinds; natively it
// survives an escaping exception because the guard, destroyed by the unwinding, sees the exception and keeps it.
// Single-threaded by contract (law 1) and no entry point calls another, so "in progress on entry" means "abandoned" —
// or RE-ENTERED from inside an allocation (a native new_handler), which this module does not support and cannot tell
// apart: both are answered FC_ERR_POISONED. See "THE MODULE IS NOT RE-ENTRANT" in fc_master_abi.h.
// `volatile`, not for threads (there are none) but because the mark's whole job is to be seen by a LATER call after
// this one never returned: nothing obliges an optimizer to keep a store whose only reader follows an abort, so the
// mechanism is not left to what it can prove. (The code-review round; measured kept at -O2, -O3 and -flto either way.)
volatile bool g_inCall = false, g_poisoned = false;

struct CallGuard
{
    const int unwinding = std::uncaught_exceptions();
    ~CallGuard() { if (std::uncaught_exceptions() == unwinding) g_inCall = false; }
};

// The first statement of every status-returning entry point — ahead of the handle check, because an abandoned
// module has no handle this file can vouch for.
#define FC_GUARD                                                                    \
    if (g_poisoned || g_inCall) { g_poisoned = true; return FC_ERR_POISONED; }     \
    g_inCall = true;                                                                \
    const CallGuard fcGuard_ {}

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
    // wide, and it WRAPS rather than retiring the slot: 16.7 million destroys of one slot before a
    // number repeats is not a budget anybody meets, while the 8-bit retirement this replaced was one a
    // per-file render loop meets in an afternoon (8 x 255 = 2040 objects per page load, refused creates
    // included). Generation 0 is skipped so a handle never equals a small integer a caller might pass
    // by accident; that costs one comparison and cannot be observed by a test, which is why it is
    // written down here rather than pinned.
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
    FC_GUARD;
    // `*out` IS NOT TOUCHED UNTIL THE CALL SUCCEEDS. It used to be cleared on entry, which reads as the
    // careful thing and is not: a caller reusing a variable that still holds a LIVE handle would have it
    // wiped by a call that failed on the version field, and the object it named would then be
    // unreachable and undestroyable. A refused call is indistinguishable from one never made — that
    // includes its arguments.
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    std::uint32_t cfgBytes = 0;
    if (const fc_status st = checkHeader (cfg, cfgBytes); st != FC_OK) return st;
    const fc_master_config c = loadIn (cfg, cfgBytes);

    MasteringChainConfig cc {};
    if (const fc_status st = toCore (c, cc); st != FC_OK) return st;
    if (c.channels < 1 || c.channels > core::kMaxChannels) return FC_ERR_REFUSED_BY_CORE;
    const bool   delivering = c.deliveryRate != 0.0;
    const double chainRate  = delivering ? c.deliveryRate : c.sampleRate;

    std::uint32_t h = 0;
    const int idx = allocSlot (Kind::Master, h);
    if (idx < 0) return FC_ERR_EXHAUSTED;

    // ADMITTED BEFORE ANYTHING IS BUILT (law 11d). This call used to construct the instance, prepare the
    // renderer, and then let the chain allocate its way down to the first stage that refused — measured,
    // on the default geometry 392 408 bytes asked for and handed straight back on a 20 Hz rate, 394 456 on a
    // 300 ms compressor lookahead, 51 288 even when the chain refused on its own front door — and 1 668 312
    // for that same lookahead at sixteen channels and an 8192-sample quantum, since these numbers scale
    // with the geometry. On the wasm tier an
    // allocation that cannot be served is not a refusal at all (law 11d), so bytes asked for on the way to
    // saying "no" are bytes that can end the module instead.
    //
    // ONE CORE CALL DECIDES AND BUDGETS. `createBytes` returns 0 for exactly the geometries the core will
    // not build, so the verdict here and the number `fc_master_need_create` publishes cannot disagree —
    // they are the same expression. Arithmetic of this file's own is what that avoids. A delivering handle asks
    // the one core function that covers the chain at the delivery rate AND the converter.
    if (createBytesFor (c, cc) == 0u)
    {
        // `abandonSlot`, not `freeSlot`: no handle ever left this function, so there is nothing for a
        // stale one to alias and no reason to spend a generation. (With the generation at 24 bits this
        // is no longer load-bearing — it was, at 8, where 2040 refused creates exhausted a table that
        // had never issued a handle — but the distinction is the honest one and it costs nothing.)
        abandonSlot (g_slots[idx]);
        return FC_ERR_REFUSED_BY_CORE;
    }

    // `new` here, not in process(): the RT-safety claim is about the audio path, and this is a worker
    // call. Under -fno-exceptions a failure aborts rather than returning, and the instance is then POISONED
    // rather than left answering — see FC_GUARD above and "POISON COMES FIRST" in fc_master_abi.h.
    g_slots[idx].master = std::make_unique<MasterInstance>();
    auto& m = *g_slots[idx].master;
    m.cfg = cc;

    // Neither can refuse now — the geometry was admitted above, and that equivalence is pinned across the
    // whole rate x width x topology matrix rather than asserted. Kept as a belt: a stage that grows a new
    // refusal must fail loudly here rather than leave a half-built instance answering calls.
    if (! m.renderer.prepare (c.channels, kRendererBlock)
        || ! m.chain.prepare (chainRate, c.channels, cc)
        || (delivering && ! m.delivered.prepare (c.sampleRate, c.deliveryRate, c.channels, kRendererBlock)))
    {
        abandonSlot (g_slots[idx]);
        return FC_ERR_REFUSED_BY_CORE;
    }
    m.delivering = delivering;
    *out = h;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_configure (fc_master h, const fc_master_params* params,
                                         fc_master_resolved* resolved)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    // THE HANDLE'S STATE COMES SECOND, before the structs are looked at, because a call that is illegal
    // for THIS HANDLE is illegal whatever else it carries — and the contract states the order, which is
    // what makes one malformed call have one answer. It used to sit after the header checks, so a
    // configure during a stream that also had a stale ABI version answered FC_ERR_ABI_VERSION and sent
    // the caller to fix the wrong thing.
    //
    // Re-preparing would silently discard a stream in progress. Refused, and `fc_master_reset` is how a
    // caller gets back to the head of one — see the long note in fc_master_abi.h for why this call has
    // to re-prepare at all.
    if (m.audioSeen) return FC_ERR_STATE;
    std::uint32_t prmBytes = 0, resBytes = 0;
    if (const fc_status st = checkHeader (params, prmBytes); st != FC_OK) return st;
    if (const fc_status st = checkHeader (resolved, resBytes); st != FC_OK) return st;

    MasteringChainParams cp {};
    if (const fc_status st = toCore (loadIn (params, prmBytes), cp); st != FC_OK) return st;

    // NOTHING HAS MOVED UNTIL HERE. A refused call above left the chain exactly as it was, which is what
    // makes "a refused call is indistinguishable from one never made" true at this boundary too.
    m.chain.setParams (cp);
    if (! m.chain.prepare (m.chain.sampleRate(), m.chain.numChannels(), m.cfg))
        return FC_ERR_REFUSED_BY_CORE;

    m.framesIn = m.framesFlushed = 0;
    m.solverRan = false;                    // a known parameter set is back in the chain
    fc_master_resolved r {};
    fromCore (m.chain.resolved(), m.chain.tapOversampleFactor(), r);
    writeOut (resolved, r, resBytes);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_resolved_get (fc_master h, fc_master_resolved* out)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    fc_master_resolved r {};
    fromCore (s->master->chain.resolved(), s->master->chain.tapOversampleFactor(), r);
    writeOut (out, r, bytes);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_process (fc_master h, const float* in, float* out, std::uint32_t frames)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();

    // THE HANDLE'S STATE COMES BEFORE THE NO-OP. `n == 0` moves nothing, but the contract says this
    // handle cannot be processed at all until it is configured, and a call that answers FC_OK on a
    // handle the header says is refused makes the sentence false for one input.
    if (m.delivering) return FC_ERR_STATE;  // two lengths cannot share one stride — fc_master_render_delivered
    if (m.solverRan) return FC_ERR_STATE;   // the chain holds the SOLVER's parameters — see MasterInstance
    if (frames == 0)
    {
        // `n == 0` is the one true no-op in law 11: no time, no falling edge, nothing. It is not an
        // error, and it must not set `audioSeen` either — no audio was seen.
        return FC_OK;
    }
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
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();
    if (m.delivering) return FC_ERR_STATE;  // as process(): a delivering handle has no stream to drain
    if (m.solverRan) return FC_ERR_STATE;   // as process(): the chain's configuration is not the caller's
    if (const fc_status st = checkScalarOut (written); st != FC_OK) return st;
    if (capacity == 0) return FC_ERR_CAPACITY;
    if (capacity > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    // A capacity below the latency cannot drain the tail, and the core keeps no arrears, so a second
    // call would not continue this drain — it would start another one. Refused as a whole (law 11).
    if ((std::int64_t) capacity < (std::int64_t) m.chain.latencySamples()) return FC_ERR_CAPACITY;
    if (const fc_status st = checkAudio (out, capacity, nch); st != FC_OK) return st;
    if (aliasesSpan (written, sizeof (*written), out,
                     (std::uint64_t) capacity * (std::uint64_t) nch * sizeof (float))) return FC_ERR_SPAN;
    // CLEARED ONLY ONCE EVERY REFUSAL IS BEHIND US. Clearing on entry meant a call refused for aliasing
    // had already written a zero into the caller's audio — a refusal that moved something.
    *written = 0;

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
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    *out = s->master->chain.latencySamples();
    return FC_OK;
}

FC_EXPORT fc_status fc_master_get_stats (fc_master h, fc_master_stats* out)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    fc_master_stats v {};
    v.framesIn      = s->master->framesIn;
    v.framesFlushed = s->master->framesFlushed;
    // A DELIVERING handle's programmes are gated by the converter, ahead of the conversion, so the chain behind it
    // never sees the bad sample; the count that means "samples of the caller's programme that were replaced" is the
    // converter's. One or the other, chosen by the handle's kind — never a sum.
    v.nonFiniteIn   = s->master->delivering ? s->master->delivered.nonFiniteInputSamples()
                                            : s->master->chain.nonFiniteInputSamples();
    writeOut (out, v, bytes);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_reset (fc_master h)
{
    FC_GUARD;
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
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    freeSlot (*s);
    return FC_OK;
}

// Forwarding only: every number is the core's (TargetLoudnessSolver's budgets) or this file's own `sizeof`, and
// none is added to another here — the page sums what applies.
FC_EXPORT fc_status fc_master_need (fc_master h, std::int32_t op, std::uint32_t frames, fc_need* out)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    if (frames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;   // the core takes `int`
    auto& m = *s->master;
    const double fs  = m.chain.sampleRate();
    const int    nch = m.chain.numChannels();
    // A DELIVERING HANDLE budgets the calls it can make, for the `frames` the caller hands IN — see fc_need. A
    // delivered length the core cannot receive is the same refusal the call itself would give.
    const double srcFs = m.delivered.sourceRate(), dstFs = m.delivered.deliveryRate();
    if (m.delivering && op != FC_NEED_CONFIGURE)
    {
        const long long d = DeliveredMastering::deliveredFrames (srcFs, dstFs, (long long) frames);
        if (d < 0 || d > 0x7FFFFFFFLL) return FC_ERR_RANGE;
    }
    std::uint64_t call = 0, facade = 0;
    bool solverOp = true;                  // the solver's two fields are NEUTRAL for a configure
    switch (op)
    {
        case FC_NEED_SOLVE:       call = m.delivering ? DeliveredMastering::solveBytes (srcFs, dstFs, nch, (long long) frames)
                                                      : TargetLoudnessSolver::solveBytes (fs, nch, (int) frames);
                                  facade = sizeof (LoudnessSolution);                     break;
        case FC_NEED_MEASURE_LRA: call = m.delivering ? DeliveredMastering::measureRangeBytes (srcFs, dstFs, nch, (long long) frames)
                                                      : TargetLoudnessSolver::measureRangeBytes (fs, (int) frames);
                                  facade = 0;                                             break;
        // A re-preparation has no programme, so the count is not ignored — it is REQUIRED to be 0. An
        // argument a call reads as nothing is an argument a caller can be wrong about for ever.
        case FC_NEED_CONFIGURE:   if (frames != 0u) return FC_ERR_RANGE;
                                  // The core's own answer about the chain in front of it: `configure`
                                  // re-prepares at the handle's own rate, width and config, and a chain
                                  // that already holds that geometry asks for nothing.
                                  call = m.chain.reprepareBytes (fs, nch, m.cfg);
                                  facade = 0; solverOp = false;                           break;
        default:                  return FC_ERR_ENUM;
    }
    fc_need v {};
    v.callBytes          = call;
    // The arguments are the ones this file hands the solver's prepare() itself — see fc_master_solve.
    v.solverPrepareBytes = solverOp ? TargetLoudnessSolver::prepareBytes (m.renderer.blockSize(),
                                                                         m.chain.internalBlock(),
                                                                         m.chain.tapOversampleFactor())
                                    : 0u;
    v.facadeBytes        = facade;
    v.solverPrepared     = (solverOp && m.solver.isPrepared()) ? 1 : 0;
    writeOut (out, v, bytes);
    return FC_OK;
}

// A DRY RUN OF fc_master_create, and the only budget with no handle to ask — see fc_need and the note at
// the declaration. The order is the create's own: poison, then this call's own out-pointer and the config's
// header, then the mapping, then a free slot, then the geometry. Every one of those refusals is the status
// the create would return, so FC_OK here means the create can only fail on the heap itself.
FC_EXPORT fc_status fc_master_need_create (const fc_master_config* cfg, fc_need* out)
{
    FC_GUARD;
    std::uint32_t outBytes = 0, cfgBytes = 0;
    if (const fc_status st = checkHeader (out, outBytes); st != FC_OK) return st;
    if (const fc_status st = checkHeader (cfg, cfgBytes); st != FC_OK) return st;
    const fc_master_config c = loadIn (cfg, cfgBytes);

    MasteringChainConfig cc {};
    if (const fc_status st = toCore (c, cc); st != FC_OK) return st;
    if (c.channels < 1 || c.channels > core::kMaxChannels) return FC_ERR_REFUSED_BY_CORE;

    // THE TABLE IS PART OF THE ANSWER. Without this a budget was published for a create that returns
    // FC_ERR_EXHAUSTED with every slot taken — a number for a call that cannot be made.
    bool haveSlot = false;
    for (int i = 0; i < kMaxHandles && ! haveSlot; ++i) haveSlot = (g_slots[i].kind == Kind::Free);
    if (! haveSlot) return FC_ERR_EXHAUSTED;

    // The same expression `fc_master_create` decides by, so the two cannot disagree: 0 is exactly the
    // geometries the core will not build, and it never reaches `out`.
    const std::uint64_t call = createBytesFor (c, cc);
    if (call == 0u) return FC_ERR_REFUSED_BY_CORE;

    fc_need v {};
    v.callBytes          = call;
    v.solverPrepareBytes = 0u;                        // neutral: a create does not touch the search
    v.facadeBytes        = sizeof (MasterInstance);
    v.solverPrepared     = 0;
    writeOut (out, v, outBytes);
    return FC_OK;
}

FC_EXPORT fc_status fc_master_measure_lra (fc_master h, const float* in, std::uint32_t frames, double* out)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (out); st != FC_OK) return st;
    auto& m = *s->master;
    const int nch = m.chain.numChannels();
    // `frames == 0` is the CORE's answer to give (`measureInputLoudnessRange` returns false on it), the
    // same reasoning as in `fc_master_solve`. Removing the facade's own copy of that verdict leaves the
    // observable result identical and the policy singular.
    if (frames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    if (m.delivering)
    {
        const long long d = DeliveredMastering::deliveredFrames (m.delivered.sourceRate(), m.delivered.deliveryRate(),
                                                                 (long long) frames);
        if (d < 0 || d > 0x7FFFFFFFLL) return FC_ERR_RANGE;
    }
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
    // publishing a measurement it had been told not to trust. On a delivering handle the core converts
    // first and measures the delivered programme — the one the search will meter.
    const bool measured = m.delivering
        ? m.delivered.measureInputLoudnessRange (m.solver, pl, nch, (long long) frames, lra)
        : m.solver.measureInputLoudnessRange (pl, nch, (int) frames, lra);
    if (! measured) return FC_ERR_REFUSED_BY_CORE;
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
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (channel < 0 || channel >= core::kMaxChannels) return FC_ERR_RANGE;
    if (! std::isfinite (weight)) return FC_ERR_NON_FINITE;
    // A finite negative weight is out of RANGE, and calling it non-finite was a lie about a number the
    // caller can see. (This is the one place the facade does bound a value: the core's own setter is a
    // silent no-op on it, so there is no verdict to forward — the same reason non-finite fields are
    // refused here.)
    if (weight < 0.0) return FC_ERR_RANGE;
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
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    // THE HANDLE'S STATE FIRST, and a search is not exempt from it. `configure` refuses to re-prepare
    // over a stream in progress; `solve` used to reset the chain out from under exactly such a stream
    // and answer FC_OK — measured, `process(64)` then a valid solve and the 64 frames in the FIFO were
    // gone with no refusal anywhere. One policy: `fc_master_reset` is the way out of a stream, for both.
    if (m.delivering) return FC_ERR_STATE;  // fc_master_solve_delivered: one stride cannot carry two lengths
    if (m.audioSeen) return FC_ERR_STATE;
    if (const fc_status st = checkScalarOut (out_solution); st != FC_OK) return st;
    std::uint32_t prmBytes = 0, reqBytes = 0;
    if (const fc_status st = checkHeader (params, prmBytes); st != FC_OK) return st;
    if (const fc_status st = checkHeader (req, reqBytes); st != FC_OK) return st;

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
    // THE OUT-HANDLE MAY NOT POINT INTO ANY SPAN THIS CALL READS OR WRITES — not just `out`. It used to
    // be cleared on entry, so an `out_solution` inside `in` zeroed an input SAMPLE before the search saw
    // it (measured: the delivered audio then differed from an honest solve in 381 075 of 384 000
    // samples), and one inside `params` zeroed a parameter FIELD before the mapping read it, after which
    // the caller read the handle's bits back as that field's value. Nothing is written until FC_OK now,
    // and all three spans are checked.
    if (aliasesSpan (out_solution, sizeof (*out_solution), out, bytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), in, bytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), params, prmBytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), req, reqBytes)) return FC_ERR_SPAN;

    MasteringChainParams cp {};
    if (const fc_status st = toCore (loadIn (params, prmBytes), cp); st != FC_OK) return st;
    LoudnessRequest lr {};
    if (const fc_status st = toCore (loadIn (req, reqBytes), lr); st != FC_OK) return st;

    if (! m.solver.isPrepared()
        && ! m.solver.prepare (m.chain.sampleRate(), nch, m.renderer.blockSize(),
                               m.chain.internalBlock(), m.chain.tapOversampleFactor()))
        return FC_ERR_REFUSED_BY_CORE;

    std::uint32_t sh = 0;
    const int idx = allocSlot (Kind::Solution, sh);
    if (idx < 0) return FC_ERR_EXHAUSTED;
    // (`*out_solution` is still untouched here, and stays so until the bottom of this function.)

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
    FC_GUARD;
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    const LoudnessSolution& v = *s->solution;
    fc_solution_summary o {};
    o.status              = (std::int32_t) v.status;
    o.binding             = (std::int32_t) v.binding;
    o.alsoViolated        = v.alsoViolated;
    o.preLimiterGainDb    = v.preLimiterGainDb;
    o.ceilingDbTp         = v.ceilingDbTp;
    o.passes              = v.passes;
    o.logCount            = v.logCount;
    o.activityThresholdDb = v.activityThresholdDb;
    o.achievedBelowLufs   = v.achievedBelowLufs;
    o.achievedAboveLufs   = v.achievedAboveLufs;
    o.gainBelowDb         = v.gainBelowDb;
    o.gainAboveDb         = v.gainAboveDb;
    writeOut (out, o, bytes);
    return FC_OK;
}

FC_EXPORT fc_status fc_solution_measurement (fc_solution sh, fc_measurement* out)
{
    FC_GUARD;
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    fc_measurement o {};
    fromCore (s->solution->measured, o);
    // v4: the traces are the solution's, not the measurement's, so they are read from it here. A caller at v1..v3
    // gets none of these four — `writeOut` stops at its `structSize`.
    const LoudnessSolution& v = *s->solution;
    o.compressorGrTraceBuckets = v.compressorTrace.buckets;
    o.limiterGrTraceBuckets    = v.limiterTrace.buckets;
    o.compressorGrTraceValid   = v.compressorTrace.valid ? 1 : 0;
    o.limiterGrTraceValid      = v.limiterTrace.valid ? 1 : 0;
    writeOut (out, o, bytes);
    return FC_OK;
}

// `written` MAY NOT POINT INTO THE RECORDS, and is cleared only once every refusal is behind the call — the order
// `fc_master_flush` takes. It used to be cleared on entry, before `out` was checked at all: a `written` inside the
// buffer then took a zero into the caller's records on a call refused for alignment or span, and a successful call
// wrote the count over a record it had just copied, answering FC_OK.
FC_EXPORT fc_status fc_solution_log (fc_solution sh, fc_solve_pass* out, std::uint32_t cap,
                                     std::uint32_t* written)
{
    FC_GUARD;
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (written); st != FC_OK) return st;
    if (cap > 0)                                      // asking for nothing is not an error, and names no span
    {
        if (out == nullptr) return FC_ERR_NULL;
        if ((reinterpret_cast<std::uintptr_t> (out) & 0x7u) != 0) return FC_ERR_ALIGNMENT;
        if (! inHeap (out, (std::uint64_t) cap * sizeof (fc_solve_pass))) return FC_ERR_SPAN;
        if (aliasesSpan (written, sizeof (*written), out, (std::uint64_t) cap * sizeof (fc_solve_pass))) return FC_ERR_SPAN;
    }
    *written = 0;

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

// v4 — the order is the header's: poison, handle, `written`, then `out` only when there is something to write
// into (null, alignment, the span, and `written` not inside it), then `stage` — a field value, checked with `cap == 0`
// as well, so a stage code that names nothing is never answered FC_OK. `written` is cleared only once every refusal
// is behind us, as in `fc_master_flush`: a `written` that pointed into the buckets used to be zeroed by a call that
// then went on to write them — and a successful call wrote the count over the first bucket's `samples`.
FC_EXPORT fc_status fc_solution_gr_trace (fc_solution sh, std::int32_t stage, fc_gr_trace_bucket* out,
                                          std::uint32_t cap, std::uint32_t* written)
{
    FC_GUARD;
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    if (const fc_status st = checkScalarOut (written); st != FC_OK) return st;
    if (cap > 0)
    {
        if (out == nullptr) return FC_ERR_NULL;
        if ((reinterpret_cast<std::uintptr_t> (out) & 0x7u) != 0) return FC_ERR_ALIGNMENT;
        if (! inHeap (out, (std::uint64_t) cap * sizeof (fc_gr_trace_bucket))) return FC_ERR_SPAN;
        if (aliasesSpan (written, sizeof (*written), out, (std::uint64_t) cap * sizeof (fc_gr_trace_bucket))) return FC_ERR_SPAN;
    }
    const LoudnessSolution& v = *s->solution;
    const GainReductionTrace* t = nullptr;
    switch (stage)                             // on the CODE, not a cast to fc_gr_stage: an int outside the enum's
    {                                          // range is not a value of it, which is the case this must refuse
        case FC_GR_STAGE_COMPRESSOR: t = &v.compressorTrace; break;
        case FC_GR_STAGE_LIMITER:    t = &v.limiterTrace;    break;
        default:                     return FC_ERR_ENUM;
    }
    *written = 0;

    const std::uint32_t n = (std::uint32_t) t->buckets < cap ? (std::uint32_t) t->buckets : cap;
    for (std::uint32_t i = 0; i < n; ++i)
    {
        const GainReductionTraceBucket& b = t->bucket[i];
        out[i].maxDb     = b.maxDb;
        out[i].meanDb    = b.meanDb;
        out[i].samples   = b.samples;
        out[i].nonFinite = b.nonFinite;
    }
    *written = n;
    return FC_OK;
}

FC_EXPORT fc_status fc_solution_destroy (fc_solution sh)
{
    FC_GUARD;
    Slot* s = lookup (sh, Kind::Solution);
    if (s == nullptr) return FC_ERR_HANDLE;
    freeSlot (*s);
    return FC_OK;
}

//==============================================================================
// THE DELIVERING HANDLE (v2) — forwarding to `mastering::DeliveredMastering`, which owns the composition. This
// file checks what a page can hand it and copies; the length arithmetic, the conversion and the render are the
// core's, and `fcore_master selftest` calls the same class directly and compares.

namespace
{
// The two lengths, in the header's order: NARROWING first — either count, or the delivered length the core
// computes, past INT_MAX — and only then a caller's `outFrames` that is not that length. A comparison, not a
// derivation: the number `outFrames` is held against comes out of the core.
fc_status checkDeliveredLengths (const MasterInstance& m, std::uint32_t inFrames, std::uint32_t outFrames) noexcept
{
    if (inFrames > (std::uint32_t) 0x7FFFFFFFu || outFrames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    const long long d = DeliveredMastering::deliveredFrames (m.delivered.sourceRate(), m.delivered.deliveryRate(),
                                                             (long long) inFrames);
    if (d < 0 || d > 0x7FFFFFFFLL) return FC_ERR_RANGE;
    if ((long long) outFrames != d) return FC_ERR_CAPACITY;
    return FC_OK;
}
}   // namespace

FC_EXPORT fc_status fc_master_delivered_frames (fc_master h, std::uint32_t inFrames, std::uint32_t* outFrames)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    if (! m.delivering) return FC_ERR_STATE;
    if (const fc_status st = checkScalarOut (outFrames); st != FC_OK) return st;
    if (inFrames > (std::uint32_t) 0x7FFFFFFFu) return FC_ERR_RANGE;
    const long long d = DeliveredMastering::deliveredFrames (m.delivered.sourceRate(), m.delivered.deliveryRate(),
                                                             (long long) inFrames);
    if (d < 0 || d > 0x7FFFFFFFLL) return FC_ERR_RANGE;
    *outFrames = (std::uint32_t) d;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_render_delivered (fc_master h, const float* in, std::uint32_t inFrames,
                                                float* out, std::uint32_t outFrames)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    if (! m.delivering) return FC_ERR_STATE;
    if (m.solverRan) return FC_ERR_STATE;   // the chain holds the SOLVER's parameters — see MasterInstance
    if (const fc_status st = checkDeliveredLengths (m, inFrames, outFrames); st != FC_OK) return st;
    const int nch = m.chain.numChannels();
    if (const fc_status st = checkAudio (in, inFrames, nch); st != FC_OK) return st;
    if (const fc_status st = checkAudio (out, outFrames, nch); st != FC_OK) return st;
    // The conversion writes `out` while it reads `in`, and the two have different strides — no overlap of any
    // kind is meaningful, equality included.
    const std::uint64_t inBytes  = (std::uint64_t) inFrames  * (std::uint64_t) nch * sizeof (float);
    const std::uint64_t outBytes = (std::uint64_t) outFrames * (std::uint64_t) nch * sizeof (float);
    if ((in == out && inBytes + outBytes > 0u) || aliasesSpan (in, inBytes, out, outBytes)) return FC_ERR_SPAN;

    const float* ip[core::kMaxChannels] {};
    float*       op[core::kMaxChannels] {};
    for (int c = 0; c < nch; ++c)
    {
        ip[c] = in  + (std::size_t) c * (std::size_t) inFrames;
        op[c] = out + (std::size_t) c * (std::size_t) outFrames;
    }
    if (! m.delivered.render (m.chain, m.renderer, ip, nch, (long long) inFrames, op, (long long) outFrames))
        return FC_ERR_REFUSED_BY_CORE;
    return FC_OK;
}

FC_EXPORT fc_status fc_master_solve_delivered (fc_master h, const fc_master_params* params,
                                               const fc_loudness_request* req,
                                               const float* in, std::uint32_t inFrames,
                                               float* out, std::uint32_t outFrames,
                                               fc_solution* out_solution)
{
    FC_GUARD;
    Slot* s = lookup (h, Kind::Master);
    if (s == nullptr) return FC_ERR_HANDLE;
    auto& m = *s->master;
    if (! m.delivering) return FC_ERR_STATE;
    if (const fc_status st = checkScalarOut (out_solution); st != FC_OK) return st;
    std::uint32_t prmBytes = 0, reqBytes = 0;
    if (const fc_status st = checkHeader (params, prmBytes); st != FC_OK) return st;
    if (const fc_status st = checkHeader (req, reqBytes); st != FC_OK) return st;
    if (const fc_status st = checkDeliveredLengths (m, inFrames, outFrames); st != FC_OK) return st;

    const int nch = m.chain.numChannels();
    if (const fc_status st = checkAudio (in, inFrames, nch); st != FC_OK) return st;
    if (const fc_status st = checkAudio (out, outFrames, nch); st != FC_OK) return st;
    const std::uint64_t inBytes  = (std::uint64_t) inFrames  * (std::uint64_t) nch * sizeof (float);
    const std::uint64_t outBytes = (std::uint64_t) outFrames * (std::uint64_t) nch * sizeof (float);
    // As `fc_master_solve`: a search cannot render over what it reads, so any touching is refused, equality included.
    if (in == out || aliasesSpan (in, inBytes, out, outBytes)) return FC_ERR_SPAN;
    // THE AUDIO OUTPUT MAY NOT OVERWRITE THE CALLER'S CONST STRUCTS either: they are copied before the search, so the
    // call would succeed — and leave the caller's parameter set or request rewritten with audio (the diverse-testing
    // round). At the caller's size, like every other check on them.
    if (aliasesSpan (params, prmBytes, out, outBytes) || aliasesSpan (req, reqBytes, out, outBytes)) return FC_ERR_SPAN;
    if (aliasesSpan (out_solution, sizeof (*out_solution), out, outBytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), in, inBytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), params, prmBytes)
        || aliasesSpan (out_solution, sizeof (*out_solution), req, reqBytes)) return FC_ERR_SPAN;

    MasteringChainParams cp {};
    if (const fc_status st = toCore (loadIn (params, prmBytes), cp); st != FC_OK) return st;
    LoudnessRequest lr {};
    if (const fc_status st = toCore (loadIn (req, reqBytes), lr); st != FC_OK) return st;

    // THE SLOT BEFORE THE SOLVER'S PREPARATION, the other way round from `fc_master_solve`: a full table is a refusal
    // that must not have prepared anything on its way — on the wasm tier a preparation that cannot be served is not a
    // status but the end of the module (the diverse-testing round). A refused preparation hands the slot back.
    std::uint32_t sh = 0;
    const int idx = allocSlot (Kind::Solution, sh);
    if (idx < 0) return FC_ERR_EXHAUSTED;
    if (! m.solver.isPrepared()
        && ! m.solver.prepare (m.chain.sampleRate(), nch, m.renderer.blockSize(),
                               m.chain.internalBlock(), m.chain.tapOversampleFactor()))
    {
        abandonSlot (g_slots[idx]);
        return FC_ERR_REFUSED_BY_CORE;
    }

    const float* ip[core::kMaxChannels] {};
    float*       op[core::kMaxChannels] {};
    for (int c = 0; c < nch; ++c)
    {
        ip[c] = in  + (std::size_t) c * (std::size_t) inFrames;
        op[c] = out + (std::size_t) c * (std::size_t) outFrames;
    }
    g_slots[idx].solution = std::make_unique<LoudnessSolution> (
        m.delivered.solve (m.solver, m.chain, m.renderer, cp, ip, nch, (long long) inFrames,
                           op, (long long) outFrames, lr));

    // As `fc_master_solve`: where the search ran, the chain holds its parameters; where it did not, nothing moved.
    const MasteringSolveStatus verdict = g_slots[idx].solution->status;
    if (verdict != MasteringSolveStatus::NotPrepared && verdict != MasteringSolveStatus::InvalidRequest)
    {
        m.framesIn = m.framesFlushed = 0;
        m.audioSeen = false;
        m.solverRan = true;
    }
    *out_solution = sh;
    return FC_OK;
}

//==============================================================================
// DEFAULTS — the core's own, written through the same mapping every other value crosses by.
//
// `writeDefaults` fills THIS build's layout. The published writers are two views of it (rule 8 of VERSIONING):
// `fc_*_defaults` write the caller's stamped version, `fc_*_default` are frozen at v1.

namespace
{
template <typename T>
void writeFrozenV1 (T* out) noexcept
{
    if (out == nullptr) return;
    T local {};
    writeDefaults (local);
    const std::uint32_t v1 = sizeFor (AbiId<T>::id, 1u);
    local.header.abiVersion = 1u;
    local.header.structSize = v1;
    std::memcpy (static_cast<void*> (out), &local, v1);
}

template <typename T>
fc_status writeVersioned (T* out) noexcept
{
    std::uint32_t bytes = 0;
    if (const fc_status st = checkHeader (out, bytes); st != FC_OK) return st;
    T local {};
    writeDefaults (local);
    writeOut (out, local, bytes);
    return FC_OK;
}

void stamp (fc_header& h, std::uint32_t size) noexcept
{
    h.abiVersion = FC_MASTER_ABI_VERSION;
    h.structSize = size;
}

void writeDefaults (fc_master_config& o) noexcept
{
    const MasteringChainConfig d {};
    std::memset (static_cast<void*> (&o), 0, sizeof (o));
    stamp (o.header, (std::uint32_t) sizeof (o));
    // LEFT AT ZERO, DELIBERATELY. There is no core default for either, and writing one here would be
    // this file choosing a geometry for every caller who forgot to — the same objection that keeps
    // `targetLufs` and `maxTruePeakDbTp` at NaN in the request. `fc_master_create` refuses both, so a
    // caller who forgets is told rather than silently given 48 kHz stereo.
    o.sampleRate            = 0.0;
    o.channels              = 0;
    o.internalBlock         = d.internalBlock;
    o.eq                    = d.eq ? 1 : 0;
    o.monoBass              = d.monoBass ? 1 : 0;
    o.compressor            = d.compressor ? 1 : 0;
    o.clipper               = d.clipper ? 1 : 0;
    o.limiter               = d.limiter ? 1 : 0;
    o.dither                = d.dither ? 1 : 0;
    o.compressorLookaheadMs = d.compressorLookaheadMs;
    o.limiterLookaheadMs    = d.limiterLookaheadMs;
    o.oversampleFactor      = d.oversampleFactor;
    o.tapsPerPhase          = d.tapsPerPhase;
    o.sidechainHpfHz        = d.sidechainHpfHz;
    // v2. Not a core default — the core has no such field — but the ABI's encoding of "no conversion", and
    // the value under which a v2 config is exactly a v1 config (rule 2).
    o.deliveryRate          = 0.0;
}

void writeDefaults (fc_master_params& o) noexcept
{
    fc_master_params* out = &o;
    const MasteringChainParams d {};
    std::memset (static_cast<void*> (out), 0, sizeof (*out));
    stamp (out->header, (std::uint32_t) sizeof (*out));
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
    out->compressorMix    = d.compressorMix;            // v3: 1, the chain before the field existed
}

void writeDefaults (fc_loudness_request& o) noexcept
{
    fc_loudness_request* out = &o;
    const LoudnessRequest d {};
    std::memset (static_cast<void*> (out), 0, sizeof (*out));
    stamp (out->header, (std::uint32_t) sizeof (*out));
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
}   // namespace

FC_EXPORT fc_status fc_master_config_defaults (fc_master_config* out)       { FC_GUARD; return writeVersioned (out); }
FC_EXPORT fc_status fc_master_params_defaults (fc_master_params* out)       { FC_GUARD; return writeVersioned (out); }
FC_EXPORT fc_status fc_loudness_request_defaults (fc_loudness_request* out) { FC_GUARD; return writeVersioned (out); }

FC_EXPORT void fc_master_config_default (fc_master_config* out)             { writeFrozenV1 (out); }
FC_EXPORT void fc_master_params_default (fc_master_params* out)             { writeFrozenV1 (out); }
FC_EXPORT void fc_loudness_request_default (fc_loudness_request* out)       { writeFrozenV1 (out); }

//==============================================================================
FC_EXPORT std::uint32_t fc_master_abi_version  (void) { return FC_MASTER_ABI_VERSION; }
FC_EXPORT std::uint32_t fc_master_max_channels (void) { return (std::uint32_t) core::kMaxChannels; }
FC_EXPORT std::uint32_t fc_master_max_eq_bands (void) { return (std::uint32_t) FC_MAX_EQ_BANDS; }
FC_EXPORT std::uint32_t fc_master_sizeof (std::int32_t id, std::uint32_t version) { return sizeFor (id, version); }
// Frozen at v1 — see the declaration.
FC_EXPORT std::uint32_t fc_master_sizeof_params(void) { return sizeFor (FC_STRUCT_PARAMS, 1u); }
FC_EXPORT std::uint32_t fc_master_sizeof_config(void) { return sizeFor (FC_STRUCT_CONFIG, 1u); }
