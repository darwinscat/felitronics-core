// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#ifndef FC_MASTER_ABI_H
#define FC_MASTER_ABI_H

// fc_master — the C ABI over `felitronics::mastering`, for a browser worker (wasm32) and for the native
// reference CLI that proves the two agree. It is the DECLARATION only: the implementation is
// tools/wasm/fc_master.cpp, which compiles natively as well (EMSCRIPTEN_KEEPALIVE degrades to extern "C").
//
// ====================================================================================
// THE FACADE IS THIN, AND THAT IS A LAW RATHER THAN A STYLE
// ====================================================================================
// The desktop application links `felitronics::mastering` AS C++, past this file entirely. So this
// surface may never become the only road to a capability: everything it can do has to be reachable
// without it. What is duplicated here is exactly two things — the ENUM CODES and the FIELD MAPPING.
// No arithmetic. Not one clamp, not one default, not one derived number. Every value this header hands
// back was computed by the core and read out of it — save one, `fc_need.facadeBytes`, the facade's own
// `sizeof`, named as such and never added to a core number; where the core cannot answer, this file
// refuses rather than inventing.
//
// The mechanical proof of that is `fcore_master selftest`: the same programme rendered through this
// ABI and through a direct C++ call, in ONE binary on ONE machine, compared bit for bit. A facade that
// had grown arithmetic of its own would have to be wrong in exactly the same way twice to pass it.
//
// ====================================================================================
// SHAPE
// ====================================================================================
// * PLANAR float32, ONE pointer, channel c at `planar + c * frames`. Same choice and same reason as
//   fc_probe.cpp: a `const float* const*` across the wasm boundary would mean building a table of i32
//   offsets in the heap and exporting HEAPU32 to write it, for no gain — and AudioBuffer.getChannelData(c)
//   is already planar, so the page does one HEAPF32.set() per channel and no de-interleave loop.
// * NO CALLBACKS. The block loop lives in JS inside the worker; nothing here calls back into the host.
// * `uint32_t` frame counts. wasm32 is a 32-bit target and a signed frame count invites an overflow
//   that cannot happen on the 64-bit machine this core was written on.
// * EVERY entry point returns `fc_status`. Nothing returns a value in band with an error.
//
// ====================================================================================
// WHY A STATUS CODE HERE, WHEN THE CORE ANSWERS `bool`
// ====================================================================================
// Law 11 (docs/DSP-ARCHITECTURE.md §2) makes every block-level `process()` a `[[nodiscard]] bool` and
// says plainly that a second refusal idiom for the same concept is how a core ends up with four
// policies. This is not that second idiom, and the difference is the CALLER.
//
// The law's argument for a bare `bool` is that "the reason is always visible at the call site — the
// caller knows what it passed". That argument holds for C++ code that passed three arguments by hand.
// It does NOT hold here: the caller is JavaScript that marshalled a 6 KB struct of several hundred
// fields into linear memory, and "which field was refused" is not visible at its call site at all.
//
// So the codes below name THIS FILE'S OWN checks — the ones the core does not make and cannot make:
// a handle that was destroyed, a span that leaves the wasm heap, a misaligned pointer, a struct whose
// version this build does not know, an enum code that names nothing. The core's own refusal maps to
// exactly ONE code, `FC_REFUSED_BY_CORE`, and this file never guesses a reason for it, because it does
// not have one: `bool` is what the core said.
//
// `fc_master_last_error(h)` was considered and rejected: it cannot answer for a handle that is itself
// invalid, it lets a stale reason outlive the call that produced it, and every getter then has to
// decide whether it clears it.
//
// ====================================================================================
// WHAT THIS FILE DOES *NOT* CHECK, AND WHY REFUSING WOULD BE WORSE
// ====================================================================================
// IT NEVER RANGE-CHECKS A PARAMETER VALUE. The core clamps them BY DESIGN and reports what it clamped
// to — `TruePeakLimiter` folds `ceilingDbTp` into [-200, +60] and the mastering suite pins that a
// ceiling of 1e308 comes back as 60. A facade that refused "+5 dBTP" would diverge from the direct C++
// call, which renders it clamped; worse, it would need its own table of every stage's range, which is
// exactly the field mapping that falls out of step the first time a stage moves one. So a value out of
// range is the CORE's business and comes back through `fc_master_resolved`.
//
// What it does check about a value is that it is FINITE, because there the core has no verdict to
// forward: `MasteringChain` maps a NaN gain to 0 dB silently, so a page that wrote NaN into a slider
// would get a render with no gain, no refusal and no way to find out.
//
// IT NEVER REFUSES A CALL OVER A NON-FINITE AUDIO SAMPLE either — see fc_master_process.

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

//==============================================================================
// VERSIONING
//
// Every struct that crosses this boundary begins with the same two fields, and they are the first two
// things read: `abiVersion` says which contract the caller was compiled against, `structSize` says how
// big it believed the struct to be. A build that does not know the version refuses; a size that does
// not match the version's size refuses. Neither is a warning — a struct read at the wrong layout is a
// silently plausible parameter set, which is precisely the failure this pair exists to prevent.
//
// EXACT size match, not ">=", for v1 — and v1 said the first version to add a field would state its own
// rule. v2 is that version (`fc_master_config::deliveryRate` and the delivering entry points), so the rule
// follows, and it is written so that the NEXT bump is a row in a table rather than an event. v3 is the first
// bump made by it — `compressorMix` (P60), two rows in the size table and a field at the end of two structs.
//
// A NEW CODE IS NOT A NEW VERSION, and the rule for codes is written here rather than left to be inferred
// from the one for structs. A status or op code is only ever APPENDED — an existing code never changes
// meaning or value — so a caller built against an older header meets a code it does not know exactly
// where it already has to handle "not FC_OK", and nothing it relied on moved. (FC_ERR_POISONED and the
// FC_NEED_* ops came in this way.)
//
// ------------------------------------------------------------------------------------------------------
// THE COMPATIBILITY RULE (v2 onward)
// ------------------------------------------------------------------------------------------------------
// 1. WHAT MOVES THE VERSION. One number for the whole ABI, and it moves when the surface a caller may use
//    grows: a struct with a header gains a field, OR an entry point is added. The second half is not
//    pedantry — a page's only protection against calling an export the module does not have (a TypeError,
//    not a status) is the loader's "module version >= page version" gate, so an entry point that arrived
//    without a bump would pass that gate and still be missing.
//
// 2. WHAT A BUMP MAY DO. One logical addition, which may append fields to the END of any number of structs
//    that begin with `fc_header` (v3, P60's `compressorMix`, is one bump that appends to `fc_master_params` and
//    to `fc_master_resolved`). An IN field carries a default under which the call behaves exactly as the
//    previous version did, bit for bit — `deliveryRate = 0` is v1. An OUT field is written with what the
//    core computed.
//
// 3. WHAT NEVER GROWS. `fc_header`; every struct nested BY VALUE (`fc_eq_band`, `fc_compressor`, `fc_gr_limit`
//    …) — a field inside `fc_compressor` would move everything after it in `fc_master_params`; and every element
//    of a header-less array (`fc_solve_pass`, written with a stride of its `sizeof`). A new stage therefore
//    arrives as fields at the end of a top-level struct, or as a NEW nested struct appended whole at the end —
//    frozen from the version that introduced it.
//
// 4. NO IMPLICIT TAIL PADDING in a struct with a header: a struct whose last field ends short of its `sizeof`
//    names the gap (`_pad0`), and the facade pins `sizeof == end of the last field`. Two reasons, neither of
//    them about JavaScript (fc-master-layout.mjs aligns fields itself): a grown struct must be strictly BIGGER
//    than its previous version, or a caller that bumped the layout and forgot the version is read at the old
//    layout with no refusal; and an old caller's padding bytes are uninitialised memory that the facade copies
//    over its defaults — a new field placed in them would be read from garbage.
//
// 5. THE SIZE OF EVERY (struct, version) IS ONE ROW OF `FC_MASTER_STRUCT_SIZES` below, and nowhere else. A
//    struct that did not change in a version inherits its previous row; a struct introduced in version N has no
//    size before N. `fc_master_sizeof(id, version)` publishes the table, 0 for a pair it does not have.
//
// 6. READING AN IN STRUCT. The first 8 bytes are bounded; the version must be 1..FC_MASTER_ABI_VERSION, else
//    FC_ERR_ABI_VERSION; `structSize` must be that version's row, else FC_ERR_STRUCT_SIZE; then exactly
//    `structSize` bytes are bounded and copied OVER THE CURRENT LAYOUT FILLED WITH ITS DEFAULTS — the same writer
//    `fc_*_defaults` uses. Every later bound and alias check on that struct uses `structSize`, not this build's
//    `sizeof`. A version NEWER than the build is refused, and that is a decision, not a gap: the build cannot
//    know what the extra bytes mean, and reading a caller's intent as "whatever fits" is the silently plausible
//    parameter set this pair of fields exists to prevent.
//
// 7. WRITING AN OUT STRUCT. Same checks. The caller's header is left as it is (an echo), and exactly its
//    `structSize` bytes are written — a field past them is not written at all. Stamping this build's own
//    version and size into a buffer laid out for an older one would advertise room the caller never allocated,
//    and the next call through that buffer would write past it.
//
// 8. DEFAULTS. `void fc_*_default(out)` cannot know how big `out` is, so they are FROZEN at v1 — they stamp v1
//    and write exactly v1's bytes, for callers compiled against v1. Everything newer uses `fc_*_defaults(out)`,
//    which requires the caller's stamp (FC_INIT, or `Struct.init()` in JS) and writes that version's size.
//    ⚠ A caller that fills a struct from a FROZEN writer and then sets a newer field has written it past the
//    stamped size, where it is not read: `deliveryRate` set that way is ignored and the handle does not convert,
//    and `compressorMix` set that way renders at mix 1. That is what rule 6 must do with a v1 struct, and it is
//    why the JS accessor refuses a field past the stamp.
//
// WHAT ONE BUMP TOUCHES — each a line, and the list is the whole of it:
//   * this file: the field(s), a row per grown struct in FC_MASTER_STRUCT_SIZES, FC_MASTER_ABI_VERSION, and
//     the entry-point declarations if any;
//   * tools/wasm/fc_master.cpp: the layout pins (sizeof / offsetof / arity), `toCore`/`fromCore`, the defaults
//     writer, and the entry points — and, for a NEW struct with a header, its `AbiId` specialisation;
//   * tools/wasm/fc-master-layout.mjs: the fields in STRUCTS and FC_MASTER_ABI_VERSION — and the copy of that
//     file the site ships (see TRANSITION below);
//   * tools/fcore_master.cpp: the hand mirror in `directRender`, the key=value parser, the offset table behind
//     `fcore_master layout`, and the selftest fixture (every new field moved off its default);
//   * tools/tests/MasterAbiTests.cpp: the version matrix.
//
// TRANSITION. The rule makes v3 cheap for a page written against v2; it cannot reach back into a page already
// shipped against v1, whose loader requires `version === 1` and fails on a v2 module before its first call.
// The move from v1 to v2 on the site is therefore a coordinated release of the worker and the module together.
#define FC_MASTER_ABI_VERSION 3u

typedef struct fc_header
{
    uint32_t abiVersion;       // FC_MASTER_ABI_VERSION
    uint32_t structSize;       // sizeof(the struct this header begins)
} fc_header;

// Stamp a struct's header. AN OUT STRUCT NEEDS THIS TOO, and that is not an oversight: the caller is
// the one saying which layout it has room for, so a facade writing its own idea of the layout into a
// buffer sized for another one is the same defect with the arrow reversed. Two stores; a page does the
// same two through HEAPU32. Not a function, because it must work before any handle exists and from C.
#define FC_INIT(s) do { (s).header.abiVersion = FC_MASTER_ABI_VERSION;                       \
                        (s).header.structSize = (uint32_t) sizeof (s); } while (0)

//==============================================================================
// STATUS
//
// Order is not accidental: it is the order the checks run in, which law 11 makes part of the contract
// so that one malformed call has one answer.
//
//   POISON -> handle -> the handle's STATE -> out-parameters (null, alignment, in-heap) -> struct
//   headers (the first 8 bytes bounded, then version, then size, THEN the rest of the struct's span) ->
//   NARROWING (a count this ABI cannot hand the core's `int`) -> audio spans and their aliases -> field
//   values -> core
//
// Three details of that are load-bearing rather than incidental. THE HANDLE'S STATE comes second
// because a call that is illegal for this handle is illegal whatever else it carries. OUT-PARAMETERS
// come before the input structs because a call that cannot report its result must not perform it. And
// NARROWING comes before the spans because it is the only check that can still fire: a frame count past
// INT_MAX makes a byte span past 32 bits too, so a span check placed first would answer every such call
// with FC_ERR_SPAN and the specific diagnosis would be unreachable.
//
// POISON COMES FIRST, AND IT IS FOR EVER (law 11d, docs/DSP-ARCHITECTURE.md). Under -fno-exceptions an exhausted
// heap does not come back as a status: the module ABORTS inside the core, and the page receives a JavaScript
// RuntimeError instead of a return value. What the abort does NOT do is stop the module — emscripten lets the
// page call again, with every object wherever the abort left it — and such an instance does not fail, it LIES.
// Measured on v0.30.0 in wasm32: an abort inside `fc_master_solve`, then `fc_master_process` on the same handle
// answered FC_OK and rendered at the search's own pass-1 gain — +12 dB in that replay, which asked for that
// start (`initialGainDb`; MasterAbiTests repeats the scenario) — and kMaxHandles − 1 = 7 such aborts left the
// handle table full for good, because each reserves its solution slot before the search runs. So every entry
// point that returns `fc_status` marks a call in progress and clears the mark only on a normal return. Finding
// it set means an earlier call never returned — an abort, a trap, or natively an exception that escaped — and
// from then on every such call answers FC_ERR_POISONED, writes nothing and touches nothing. There is no way back
// inside the instance: the page discards it and instantiates a new one. The entry points that return no status —
// the frozen v1 `*_default` writers, `fc_master_sizeof` and the build-identity queries — read no instance state
// and stay callable. The versioned `*_defaults` writers return a status and are therefore guarded like every other
// such entry point: a poisoned module refuses them too.
//
// THE MODULE IS NOT RE-ENTRANT, and the poison is what says so. An entry point called while another is still
// running — from a new_handler, a signal handler, anything the runtime runs inside an allocation this file made —
// cannot be told apart from the first call after an abandoned one, and is answered FC_ERR_POISONED, for good. (The
// code-review round: a native new_handler that destroyed a spare handle during `fc_master_create` used to work and
// now poisons the module. A browser page cannot install one; a native host must not.)

typedef enum fc_status
{
    FC_OK                  =  0,

    FC_ERR_HANDLE          =  1,   // null, never issued, already destroyed, or from a previous generation
    FC_ERR_ABI_VERSION     =  2,   // `abiVersion` is not one this build knows
    FC_ERR_STRUCT_SIZE     =  3,   // `structSize` disagrees with this build's sizeof for that version
    FC_ERR_NULL            =  4,   // a required pointer was null
    FC_ERR_ALIGNMENT       =  5,   // a float*/double* that is not 4-/8-byte aligned
    FC_ERR_SPAN            =  6,   // the span leaves the heap, or frames*channels overflows 32 bits
    FC_ERR_ENUM            =  7,   // an enum field carries a code this ABI does not define
    FC_ERR_RANGE           =  8,   // a count this ABI must narrow cannot be narrowed — `frames` above
                                   // INT_MAX, which the core's own `int` cannot receive. It is NOT used
                                   // for parameter VALUES: see the note below
    FC_ERR_CAPACITY        =  9,   // a caller-owned output buffer cannot hold what the call would write
    FC_ERR_STATE           = 10,   // the call is legal but not HERE — see the entry point's own note
    FC_ERR_NON_FINITE      = 11,   // a NaN or an infinity where the contract admits neither
    FC_ERR_REFUSED_BY_CORE = 12,   // the core returned false. This ABI does not know why, and says so
    FC_ERR_EXHAUSTED       = 13,   // no free slot in the handle table
    FC_ERR_POISONED        = 14    // an earlier call into this module never returned: the instance is
                                   // abandoned, and nothing but a new one answers — see above
} fc_status;

//==============================================================================
// ENUM CODES
//
// These are the duplication the thinness law permits, and they are duplicated ON PURPOSE rather than
// cast from the C++ enums: a `static_cast<eq::FilterType>(code)` would silently follow any reordering
// of the C++ enumerators, which is a change of meaning with no diagnostic anywhere. The mapping is a
// switch with no default fall-through, so a C++ enumerator added without a code here fails to compile,
// and a code arriving from JS that names nothing is refused with FC_ERR_ENUM rather than clamped.
typedef enum fc_filter_type
{
    FC_FILTER_BELL = 0, FC_FILTER_LOW_SHELF = 1, FC_FILTER_HIGH_SHELF = 2,
    FC_FILTER_HIGH_PASS = 3, FC_FILTER_LOW_PASS = 4, FC_FILTER_BAND_PASS = 5,
    FC_FILTER_NOTCH = 6, FC_FILTER_ALL_PASS = 7, FC_FILTER_TILT = 8
} fc_filter_type;

typedef enum fc_detector { FC_DETECTOR_PEAK = 0, FC_DETECTOR_RMS = 1 } fc_detector;

// `MeanPower`, not "average": sqrt(mean(ch^2)). The name is the core's and is kept, because a code
// called FC_LINK_AVERAGE would be a second name for one thing, and the second name is where a mapping
// starts drifting.
typedef enum fc_link_mode { FC_LINK_MAX = 0, FC_LINK_MEAN_POWER = 1 } fc_link_mode;

// Three, not two: `DownExpand` exists and a compressor mode field that could not express it would be a
// knob this ABI silently removed.
typedef enum fc_comp_mode
{
    FC_COMP_DOWN_COMPRESS = 0, FC_COMP_UP_COMPRESS = 1, FC_COMP_DOWN_EXPAND = 2
} fc_comp_mode;

typedef enum fc_shape
{
    FC_SHAPE_TANH = 0, FC_SHAPE_ATAN = 1, FC_SHAPE_CUBIC = 2, FC_SHAPE_ASYM = 3
} fc_shape;
typedef enum fc_noise_shaping { FC_SHAPING_NONE = 0, FC_SHAPING_WEIGHTED = 1, FC_SHAPING_PSYCHO = 2 } fc_noise_shaping;

// Mirrors mastering::MasteringSolveStatus. FC_OK from fc_master_solve means A VERDICT WAS OBTAINED —
// it does not mean the target was met. That is this field.
typedef enum fc_solve_status
{
    FC_SOLVE_SOLVED = 0, FC_SOLVE_TARGET_UNREACHABLE = 1, FC_SOLVE_UPSTREAM_VIOLATION = 2,
    FC_SOLVE_TARGET_BETWEEN = 3, FC_SOLVE_PASS_LIMIT = 4, FC_SOLVE_MEASUREMENT_INVALID = 5,
    FC_SOLVE_RENDER_FAILED = 6, FC_SOLVE_NOT_PREPARED = 7, FC_SOLVE_INVALID_REQUEST = 8
} fc_solve_status;

// Mirrors mastering::MasteringConstraint.
typedef enum fc_constraint
{
    FC_CONSTRAINT_NONE = 0, FC_CONSTRAINT_TRUE_PEAK = 1, FC_CONSTRAINT_LIMITER_GR = 2,
    FC_CONSTRAINT_PLR = 3, FC_CONSTRAINT_LRA = 4, FC_CONSTRAINT_GAIN_RANGE = 5,
    FC_CONSTRAINT_COMPRESSOR_GR = 6
} fc_constraint;

// Mirrors mastering::GrStatistic. Which statistic a gain-reduction limit binds is part of the limit's
// TYPE and never a hidden convention.
typedef enum fc_gr_statistic { FC_GR_MEAN = 0, FC_GR_P95 = 1, FC_GR_MAX = 2 } fc_gr_statistic;

//==============================================================================
// TOPOLOGY — fixed for the life of a handle
//
// This mirrors `mastering::MasteringChainConfig`, and it is a SEPARATE struct from the per-block
// parameters for the reason the C++ side gives: every field here moves latencySamples(), and a moving
// latency is a host resynchronisation event rather than automation. A create() that took only a sample
// rate and a channel count — the shape this task's brief proposed — would leave mono-bass and the
// clipper absent for the life of the handle with no way to ask for them, because no per-block parameter
// can bring an absent stage into existence.
typedef struct fc_master_config
{
    fc_header header;

    double  sampleRate;
    int32_t channels;

    int32_t internalBlock;          // THE QUANTUM. See the block-invariance note at fc_master_process

    int32_t eq;                     // int32_t and not a C `bool`: `_Bool` has no ABI this file controls,
    int32_t monoBass;               // and JS writes these through HEAP32
    int32_t compressor;
    int32_t clipper;
    int32_t limiter;
    int32_t dither;

    double  compressorLookaheadMs;
    double  limiterLookaheadMs;
    int32_t oversampleFactor;
    int32_t tapsPerPhase;
    double  sidechainHpfHz;

    // ---- v2 ----
    // THE RATE THE PROGRAMME IS DELIVERED AT, and the one field that decides what kind of handle this is:
    //   * 0 — no conversion: the chain runs at `sampleRate`, and the handle is exactly a v1 handle.
    //   * anything else — a DELIVERING handle. The chain, the renderer and the solver run at `deliveryRate`,
    //     and `DeliveryConverter` (SRC first, over the whole programme) stands in front of them. EQUAL RATES
    //     INCLUDED: `deliveryRate == sampleRate` is a delivering handle whose converter copies bits, so a page
    //     has one path — the `*_delivered` calls — whatever rate the user picked.
    // A delivering handle has no streaming path: `fc_master_process`, `fc_master_flush` and `fc_master_solve`
    // answer FC_ERR_STATE on it, because two lengths cannot share one planar stride. A rate pair the resampler
    // does not support (anything but integer rates it can plan — see core::DeliveryResampler) is refused at
    // create with FC_ERR_REFUSED_BY_CORE; a non-finite one with FC_ERR_NON_FINITE.
    double  deliveryRate;
} fc_master_config;

//==============================================================================
// PER-BLOCK PARAMETERS
//
// A field-for-field mirror of `mastering::MasteringChainParams` and everything it contains. Nothing is
// summarised, nothing is omitted and nothing is renamed: an ABI that exposed a "useful subset" would be
// a road that ends, and the thinness law says this surface may not be the only road to a capability —
// which is only true if it is not a narrower one either.
#define FC_MAX_EQ_BANDS 24
#define FC_MAX_EQ_LANES 5

typedef struct fc_eq_lane
{
    int32_t on;
    double  freq;
    double  q;
    double  gainDb;
    int32_t slope;
    int32_t bypass;
} fc_eq_lane;

typedef struct fc_eq_dyn
{
    int32_t on;
    double  rangeDb;
    double  thrDb;
    int32_t thrAuto;
    double  atk;
    double  rel;
} fc_eq_dyn;

typedef struct fc_eq_band
{
    int32_t     on;
    int32_t     type;               // fc_filter_type
    int32_t     swept;
    int32_t     bypass;
    fc_eq_dyn   dyn;
    fc_eq_lane  lanes[FC_MAX_EQ_LANES];
} fc_eq_band;

typedef struct fc_mono_bass
{
    int32_t enabled;
    float   frequencyHz;
    float   lowWidth;
} fc_mono_bass;

typedef struct fc_compressor
{
    int32_t detector;               // fc_detector
    int32_t link;                   // fc_link_mode
    double  rmsWindowMs;
    int32_t mode;                   // fc_comp_mode
    double  thresholdDb;
    double  ratio;
    double  kneeDb;
    double  rangeDb;
    double  attackMs;
    double  releaseMs;
    double  makeupDb;
    int32_t autoMakeup;
    // `lookaheadMs` is deliberately ABSENT. It is topology — it moves latency — and the chain
    // overwrites it from the config on every apply. A field here would be a knob that does nothing,
    // which is the silent-no-op class this project keeps closing.
} fc_compressor;

typedef struct fc_clipper
{
    int32_t shape;                  // fc_shape
    float   driveDb;
    float   bias;
    float   mix;
    float   outputDb;
    float   autoComp;
    float   dcBlockHz;
} fc_clipper;

typedef struct fc_limiter
{
    double ceilingDbTp;
    double releaseMs;
} fc_limiter;

typedef struct fc_dither
{
    int32_t  bits;
    int32_t  shaping;               // fc_noise_shaping
    // The seed is 64 bits and JS has no 64-bit integer in a Number: 0x853c49e6748fea9b loses its low
    // bits through a double. Two uint32 halves cross the boundary exactly, and the page writes them as
    // two HEAPU32 stores. `seedLo` is the low half.
    uint32_t seedLo;
    uint32_t seedHi;
    int32_t  autoBlank;
    int32_t  autoBlankSamples;
} fc_dither;

typedef struct fc_master_params
{
    fc_header header;

    double inputGainDb;
    double preLimiterGainDb;

    fc_eq_band    eqBands[FC_MAX_EQ_BANDS];
    fc_mono_bass  monoBass;
    fc_compressor compressor;
    fc_clipper    clipper;
    fc_limiter    limiter;
    fc_dither     dither;

    int32_t bypassEq, bypassMonoBass, bypassCompressor;
    int32_t bypassClipper, bypassLimiter, bypassDither;

    // ---- v3 ----
    // PARALLEL COMPRESSION (P60): the compressor stage's output is `(1 - mix) * dry + mix * compressed`, the dry
    // path aligned to the compressor's lookahead. A field of the CHAIN, not of `fc_compressor` — rule 3 of
    // VERSIONING, and the core's own choice besides (`dynamics::Compressor` keeps dry/wet out on purpose). The
    // default 1 renders the chain as it was before the field existed, bit for bit, so a v1 or v2 parameter set
    // read over this build's defaults renders exactly what it did. Clamped to [0, 1] by the core and read back
    // through `fc_master_resolved::compressorMix`; a non-finite value is refused here (FC_ERR_NON_FINITE), because
    // the core would map it to 1 without a word.
    double compressorMix;
} fc_master_params;

//==============================================================================
// WHAT THE CHAIN ACTUALLY APPLIED
//
// A mirror of `mastering::MasteringChainResolved`. Every number in it is READ OUT of the prepared
// chain, never recomputed here — which is what makes it worth reading at all.
typedef struct fc_master_resolved
{
    fc_header header;

    int32_t latencySamples;
    int32_t internalBlock;
    int32_t compressorLookahead;
    int32_t clipperLatency;
    int32_t limiterLatency;
    int32_t limiterLookahead;
    int32_t oversampleFactor;
    int32_t compressorTapOffset;
    int32_t limiterTapOffset;
    double  limiterCeilingDbTp;
    double  limiterReleaseMs;
    fc_mono_bass monoBass;

    // The stride of the OVERSAMPLED taps, which is NOT `oversampleFactor` — that one answers "what
    // factor is this chain oversampling at" and reports the clipper's when there is no limiter, while
    // this answers "how long must my limiter tap buffer be", which with no limiter is one per frame.
    int32_t tapOversampleFactor;

    // ---- v3 ----
    // The mix the chain APPLIED — the clamped value the stage holds, not the request — and 0 on a topology with no
    // compressor, where there is nothing to mix.
    double compressorMix;
} fc_master_resolved;

//==============================================================================
// STATISTICS
//
// What the handle can say about the stream it has seen. Deliberately small: this is not a measurement
// API — `fc_probe_*` already is one, and a second one here would be a second definition of loudness.
typedef struct fc_master_stats
{
    fc_header header;

    // 64 bits, because 32 overflows on a legal stream: 4096-frame calls wrap after about 24 h 51 min at
    // 48 kHz and 6 h 13 min at 192 kHz, and a counter that can read zero after having been non-zero is
    // worse than no counter.
    uint64_t framesIn;              // frames handed to fc_master_process since the last reset
    uint64_t framesFlushed;         // frames drained by fc_master_flush since the last reset
    // Input samples the chain's own gate had to replace because they were not finite, read straight out
    // of `MasteringChain::nonFiniteInputSamples()`. Not recomputed here, and not derived: a second scan
    // of the audio would be a second definition of "non-finite", and the count of internal quanta —
    // which an earlier draft of this struct carried — would have been this file re-deriving the chain's
    // own quantum accounting. The chain owns both; this reports one and does not invent the other.
    // ON A DELIVERING HANDLE it is the core's count for the programme the last `*_delivered` call or converting
    // `fc_master_measure_lra` was handed, at the source rate: the converter gates every input sample ahead of the
    // conversion (one bad sample must not become a kernel's worth of zeroes), so the chain behind it sees none.
    uint64_t nonFiniteIn;
} fc_master_stats;

//==============================================================================
// WHAT A CALL WILL ASK THE HEAP FOR (law 11d, docs/DSP-ARCHITECTURE.md)
//
// On the wasm tier exhaustion is not a status (see POISON above), so a page that must not lose its worker
// budgets BEFORE the call. Every number here but `facadeBytes` is computed by the core with the very functions
// its prepare() sizes itself with — a budget cannot drift from its allocation — and `facadeBytes` is the facade's
// own `sizeof`; all are forwarded field by field, never summed: the page adds what applies. REQUESTED bytes, not a promise that a heap can serve them: allocator headers,
// alignment and fragmentation are the page's margin to keep.
//
// The op is a code rather than a field per call, which is how `configure` joined without moving this
// struct — and how `create`, which has no handle to ask, gets an entry point of its own below instead.
typedef enum fc_need_op
{
    FC_NEED_SOLVE       = 0,
    FC_NEED_MEASURE_LRA = 1,
    // The chain's own storage, re-prepared. `frames` must be 0 for this op — a re-preparation has no
    // programme length — and anything else is FC_ERR_RANGE rather than ignored.
    FC_NEED_CONFIGURE   = 2
} fc_need_op;

typedef struct fc_need
{
    fc_header header;

    // What the core's own requests occupy AT ONCE during one such call. What bounds it differs by op, and
    // saying so is law 11d's own instruction ("what each number bounds is written where it is defined"):
    //
    //   * SOLVE — one PASS. The search builds its meters per pass and frees them at the pass's end, so
    //     this is one pass, plus the true-peak drain the first solve keeps (afterwards an upper bound by
    //     exactly that drain).
    //   * MEASURE_LRA — one meter, and 0 for a programme too short to have a range, which the call
    //     refuses before building one.
    //   * CREATE — the SUM of what the call requests, which is what it holds: everything a create asks
    //     for, it keeps until the instance is destroyed. It exceeds the peak by the part of the seed each
    //     of the chain's three dry aligners hands back when its preparation re-sizes it — 12 bytes for
    //     an aligner re-sized whole, 4 for a mono compressor whose lookahead rounds to 0 samples — and, on a
    //     plain handle, by nothing else. A DELIVERING create also counts the converter's transient FIR
    //     prototypes, which its preparation designs the phase tables from and frees before it builds the
    //     histories: there the sum exceeds the peak by those too, and it is still exactly what the call requests.
    //   * ON A DELIVERING HANDLE (`deliveryRate != 0`), SOLVE and MEASURE_LRA budget the calls that handle can
    //     make — `fc_master_solve_delivered` and the converting `fc_master_measure_lra` — and `frames` is still
    //     the count the caller hands IN, at the source rate. Each adds the converted programme, which the call
    //     holds from start to end (0 at equal rates, where the input is read in place), to what the solver
    //     builds at the DELIVERED length — so a range too short to measure is judged on the delivered length,
    //     which is what the meter sees. `fc_master_render_delivered` has no op because it asks for nothing: it
    //     converts into the caller's output and renders there. A delivered length past INT_MAX is FC_ERR_RANGE,
    //     as the call itself would be.
    //     WHAT THAT MEANS FOR A PAGE: a delivered solve lives beside THREE programme-sized buffers — the page's
    //     input, the page's output at the delivery rate, and the converted input. Stereo 44.1 -> 192 kHz is
    //     3,424,800 bytes per second of input, so a 2 GiB heap (emscripten's growth limit) holds about ten
    //     minutes of it before anything else is counted; 44.1 -> 96 kHz, about nineteen.
    //   * CONFIGURE — 0 when the chain already holds the geometry it is being re-prepared at, which is
    //     every configure this ABI can make (the rate, the width and the config are the handle's own).
    //     The chain re-uses its EQ engine and assigns every buffer to the length it already has, so this
    //     is exactly zero rather than nearly zero: it used to be 345 224 bytes, a second 331 KiB EQ engine
    //     built before the first was destroyed plus three temporaries.
    uint64_t callBytes;
    uint64_t solverPrepareBytes;  // the search's one-time preparation — tap buffers and gain-reduction histograms —
                                  // done lazily by the first solve, measure_lra or channel weight on the handle.
                                  // NEUTRAL (0) for CREATE and CONFIGURE: neither call touches the search
    uint64_t facadeBytes;         // this file's own object for the call: a solve's solution record, an instance
                                  // record for CREATE; 0 for measure_lra and for CONFIGURE
    int32_t  solverPrepared;      // 1 once that preparation has happened: `solverPrepareBytes` is then already spent.
                                  // NEUTRAL (0) for CREATE and CONFIGURE
    int32_t  _pad0;               // the tail padding, NAMED — see rule 4 of VERSIONING. Always written 0.
} fc_need;

//==============================================================================
// THE LOUDNESS SEARCH
//
// A versioned C-POD request in, an OPAQUE HANDLE out, and the per-pass log copied into a buffer the
// CALLER owns. Marshalling `mastering::LoudnessSolution` whole would be ~2.3 KiB per call of C++ enums,
// `bool`, padding and a 32-entry log that almost every caller drops on the floor.
typedef struct fc_gr_limit
{
    double  limitDb;                // +infinity = no limit. NOT "any non-finite": -infinity is an
                                    // unsatisfiable limit and must stay one
    int32_t statistic;              // fc_gr_statistic
} fc_gr_limit;

typedef struct fc_loudness_request
{
    fc_header header;

    double targetLufs;              // REQUIRED — the core ships no default target, by decision
    double toleranceLu;
    double maxTruePeakDbTp;         // REQUIRED
    double truePeakAimDb;

    fc_gr_limit limiterGr;
    fc_gr_limit compressorGr;

    double  minPlrDb;
    double  maxLraLossLu;
    double  inputLoudnessRangeLu;   // NaN = not supplied, which switches the range constraint off
    double  activityThresholdDb;
    int32_t maxPasses;
    double  initialGainDb;          // NaN = use the params' own
} fc_loudness_request;

typedef struct fc_solve_pass
{
    double   gainDb, ceilingDb;
    double   integratedLufs, truePeakDbTp, plrDb;
    double   limiterMaxGrDb, loudnessRangeLu;
    uint32_t violated;              // bitmask over (1 << (fc_constraint - 1))
} fc_solve_pass;

typedef struct fc_gr_stats
{
    double   meanDb, p95Db, maxDb, activeFraction;
    uint64_t frames, nonFinite, aboveRange;
    int32_t  valid;
} fc_gr_stats;

typedef struct fc_measurement
{
    fc_header header;

    double integratedLufs, truePeakDbTp, samplePeakDb, loudnessRangeLu, plrDb;
    fc_gr_stats compressor, limiter;
    double  limiterMaxReconstructedPeakDb;
    int32_t latencySamples, gatingBlocks, droppedBlocks, nonFiniteSubHops;
    int32_t loudnessValid, lraValid;
} fc_measurement;

typedef struct fc_solution_summary
{
    fc_header header;

    int32_t  status;                // fc_solve_status
    int32_t  binding;               // fc_constraint
    uint32_t alsoViolated;
    double   preLimiterGainDb, ceilingDbTp;
    int32_t  passes, logCount;
    double   activityThresholdDb;
    double   achievedBelowLufs, achievedAboveLufs, gainBelowDb, gainAboveDb;
} fc_solution_summary;

//==============================================================================
// THE SIZE TABLE — rule 5 of VERSIONING. One row per (struct, first version at that size). A bump that grows a
// struct adds a row here and touches no other row; a bump that does not grow a struct adds nothing for it.
//
// The codes are part of the ABI (a page passes them to `fc_master_sizeof`) and are only ever appended.
typedef enum fc_struct_id
{
    FC_STRUCT_CONFIG = 0, FC_STRUCT_PARAMS = 1, FC_STRUCT_RESOLVED = 2, FC_STRUCT_STATS = 3,
    FC_STRUCT_NEED = 4, FC_STRUCT_REQUEST = 5, FC_STRUCT_MEASUREMENT = 6, FC_STRUCT_SUMMARY = 7
} fc_struct_id;

//      X(id,                     since, bytes)
#define FC_MASTER_STRUCT_SIZES(X)                 \
        X(FC_STRUCT_CONFIG,       1,       80)    \
        X(FC_STRUCT_CONFIG,       2,       88)    \
        X(FC_STRUCT_PARAMS,       1,     6560)    \
        X(FC_STRUCT_PARAMS,       3,     6568)    \
        X(FC_STRUCT_RESOLVED,     1,       80)    \
        X(FC_STRUCT_RESOLVED,     3,       88)    \
        X(FC_STRUCT_STATS,        1,       32)    \
        X(FC_STRUCT_NEED,         1,       40)    \
        X(FC_STRUCT_REQUEST,      1,      120)    \
        X(FC_STRUCT_MEASUREMENT,  1,      208)    \
        X(FC_STRUCT_SUMMARY,      1,       88)

//==============================================================================
// HANDLES
//
// A handle is an INDEX AND A GENERATION packed into 32 bits, not a pointer. A raw pointer cast to
// uint32_t would make `destroy` followed by any further call a use-after-free that wasm does not trap:
// there is no unmapped page in a linear memory, so the read succeeds and returns whatever the allocator
// has since put there. With a generation, a stale handle is a refusal (FC_ERR_HANDLE) instead.
// Zero is never a valid handle.
//
// THE GENERATION IS 24 BITS WIDE, and the width is the design rather than spare space. A narrow
// generation forces a choice between ABA (reuse the numbers) and RETIREMENT (spend the slot), and
// retirement turns a table of eight live objects into a LIFETIME BUDGET — at 8 bits, 2040 create/destroy
// cycles per page load, after which every correct create is refused for ever. That is reachable: the
// reference CLI creates a handle per programme, so a worker written from it would die on its 2041st
// file. At 24 bits the budget stops being a question and a fabricated handle still differs from the
// live one in one of 24 bits.
typedef uint32_t fc_master;
typedef uint32_t fc_solution;

//==============================================================================
// ENTRY POINTS

// Create a chain. `cfg` carries the sample rate, the channel count and the topology; the chain is
// prepared with the core's own default parameters, so latency and geometry are answerable at once.
//
// THE HANDLE IS WRITTEN ONLY ON FC_OK, and `out` is not touched at all otherwise — not even cleared.
// Clearing on entry reads as the careful thing and is not: a caller reusing a variable that still held
// a LIVE handle would have it wiped by a call that failed on the version field, and the object it named
// would be unreachable and undestroyable. A COUNT out-parameter (`written`) is the opposite and IS
// cleared first, because zero is the truthful answer for a call that wrote nothing.
fc_status fc_master_create (const fc_master_config* cfg, fc_master* out);

// Apply a parameter set and report back what the chain will actually run.
//
// WHY THIS RE-PREPARES, AND WHY IT IS REFUSED MID-STREAM. `MasteringChain::setParams()` defers to the
// next internal quantum — deliberately, so a parameter change lands at the same place in the stream
// however the caller cut it — and `resolved()` reads the stage state, i.e. the last APPLIED set. So a
// configure that wrote the parameters and read `resolved()` straight back would hand the caller the
// PREVIOUS set: measured on the full chain at 48 kHz, a **5.0000 dB** error on the limiter ceiling and
// **149.968 ms** on its release (a first set at -1 dBTP / 50 ms read back after a second at -6 / 200),
// with the mono-bass corner reading 120 Hz where 250 had been asked for — and on a freshly prepared
// chain it reports the core's own defaults rather than anything the caller asked for. The suite's own
// fixture uses a different pair of sets, so the numbers here are the measurement's, not that test's.
//
// Applying early instead is WORSE, and that was measured too: `dither::Dither::setParams()` reseeds its
// RNG whenever the seed changes, and `eq::EqBand::setParams()` SNAPS while uninitialised and glides
// after, so N configure calls with no audio between them do not equal one call with the last set —
// which would make a render depend on how many times a knob was moved before the button was pressed,
// and `OfflineRenderer` promises exactly the opposite.
//
// What is left is the order the core itself blesses: "configure, then prepare" — MasteringChain.h says
// in as many words that this is the order a C-ABI facade takes, and `prepare()` applies the pending set
// and then resets. So this call stores the parameters and re-prepares, which makes `resolved` exact and
// makes N calls identical to the last one alone. It COSTS NOTHING at the heap — it used to rebuild the
// 331 KiB EQ engine and pay 345 224 bytes to change nothing, and the chain now re-uses it and assigns every
// buffer to the length it already has (`FC_NEED_CONFIGURE` publishes the 0, and a suite pins it). It is
// still NOT real-time — it is a worker call between renders, and it says so.
//
// It is therefore REFUSED with FC_ERR_STATE once audio has been handed to this handle, because
// re-preparing would silently discard the stream. `fc_master_reset` is how a caller gets back to the
// head of one.
fc_status fc_master_configure (fc_master h, const fc_master_params* params, fc_master_resolved* resolved);

// Read back the resolved geometry without changing anything.
fc_status fc_master_resolved_get (fc_master h, fc_master_resolved* out);

// Process `frames` frames. `in` and `out` are planar with stride `frames`; they may be equal, and must
// not overlap in any other way.
//
// REFUSED WITH FC_ERR_STATE ON A DELIVERING HANDLE (`deliveryRate != 0`), and so are `fc_master_flush` and
// `fc_master_solve`: the output of a conversion is a different number of frames than its input, and one
// planar stride cannot describe both. A delivering handle renders through `fc_master_render_delivered`.
//
// REFUSED WITH FC_ERR_STATE ON A HANDLE THAT HAS SOLVED, until `fc_master_configure` runs. A search
// resets the chain on every pass and leaves its OWN gain and ceiling in it, standing wherever its last
// pass ended — so a process call here renders a parameter set the caller never chose. Measured: the
// output differs from the delivered render in 558 691 of 576 000 samples and from the CONFIGURED render
// in 575 998. `fc_master_reset` does NOT lift the refusal, because it clears audio state and keeps the
// solver's parameters; only a configure puts a known set back.
//
// BLOCK INVARIANCE SURVIVES THIS BOUNDARY because nothing here re-blocks anything: the call is handed
// to `MasteringChain::process()` in one piece and the chain accumulates into its own fixed quantum,
// which is the only reason "same input, same output, whatever the caller's block size" is a theorem
// rather than a hope. The FIFO lives in the chain and persists across ABI calls, so cutting a programme
// into 1-sample calls and into one whole-file call give the same bits.
//
// A NON-FINITE INPUT SAMPLE IS NOT A REFUSAL HERE, and that is a decision with a reason. The chain
// already gates every input sample (isfinite/clamp, ahead of every stage), so its behaviour on a bad
// sample is ONE rule that does not depend on which stages are on. Refusing the call instead would
// throw away every GOOD sample travelling with the bad one — and the damage would scale with the
// caller's block size, which is the one thing this whole design exists to make irrelevant. The count
// comes back in fc_master_get_stats instead, so the damage is visible without being amplified.
fc_status fc_master_process (fc_master h, const float* in, float* out, uint32_t frames);

// Drain the chain's latency. Writes min(latencySamples(), capacity) frames into `out` (planar, stride
// `capacity`) and reports how many through `written`.
//
// `written` MAY NOT POINT INTO `out`. The in/out overlap rule does not see that class at all, and the
// consequence is silent: the drain is written and then its first sample is overwritten by the frame
// count, so the caller gets audio whose first four bytes are a small integer. Refused with FC_ERR_SPAN.
// Also refused on a handle that has solved, for the reason at `fc_master_process`.
//
// CAPACITY IS AN ARGUMENT because the caller owns the buffer — the same rule as
// fc_probe_block_energies(out, cap). AND A CAPACITY BELOW THE LATENCY IS REFUSED, before anything moves,
// with FC_ERR_CAPACITY. That is not caution, it is the only honest reading of the word `flush`: the
// underlying `MasteringChain::flush()` keeps NO ARREARS — every call pushes min(D, capacity) fresh zeros
// through the chain and advances it — so a caller looping "until it returns zero" would never finish and
// would keep feeding the chain silence. A request that cannot be honoured in full is refused as a whole
// (law 11). `fc_master_latency` is how a caller sizes the buffer.
fc_status fc_master_flush (fc_master h, float* out, uint32_t capacity, uint32_t* written);

// The chain's latency in frames — what `flush` needs capacity for.
fc_status fc_master_latency (fc_master h, int32_t* out);

fc_status fc_master_get_stats (fc_master h, fc_master_stats* out);
fc_status fc_master_reset (fc_master h);
fc_status fc_master_destroy (fc_master h);

// The budget of the next such call on this handle for `frames` of programme — see fc_need. A `frames` past INT_MAX
// is refused with FC_ERR_RANGE (the core counts in `int`), and then an `op` that is not an fc_need_op with
// FC_ERR_ENUM — narrowing before field values, as everywhere. FC_NEED_CONFIGURE takes `frames == 0` and
// refuses anything else with FC_ERR_RANGE: a re-preparation has no programme length, and an ignored
// argument is an argument a caller can be wrong about for ever. Reads the handle and moves nothing.
//
// IT ANSWERS A NUMBER, NEVER A PERMISSION. The handle's own state is not consulted: a solve is budgeted
// while a stream is in progress even though `fc_master_solve` would refuse it with FC_ERR_STATE, and so is
// a configure. That is deliberate — the cost of a call does not depend on the moment it is made, and a page
// deciding whether to reset a stream and re-configure needs the number precisely when the call would be
// refused. Everything this entry point is asked about belongs to a handle whose geometry was admitted at
// `fc_master_create`, so there is nothing here that could be inadmissible. `fc_master_need_create` below is
// the opposite case and behaves the opposite way.
fc_status fc_master_need (fc_master h, int32_t op, uint32_t frames, fc_need* out);

// THE BUDGET OF A `fc_master_create` THAT HAS NOT HAPPENED — the one budget that cannot be asked through a
// handle, because the handle is what the call would produce.
//
// IT IS A DRY RUN, and that is the difference from `fc_master_need` above: every refusal `fc_master_create`
// can reach before its first allocation is returned HERE, with the same status — a malformed or
// wrong-version struct, a field the mapping rejects (FC_ERR_ENUM), a width or a geometry the core will not
// build (FC_ERR_REFUSED_BY_CORE), a handle table with no free slot (FC_ERR_EXHAUSTED). So FC_OK carries a
// budget that is ALWAYS non-zero, and a caller never has to read "0 bytes" as three different things.
//
// It is also the only way to ask whether a configuration is buildable WITHOUT paying for the attempt: a
// create refused by the core used to ask the heap for 392 408 bytes on its way to saying no on the default
// geometry, and 1 668 312 at sixteen channels and an 8192-sample quantum, on a
// tier where an allocation that fails is not a refusal but the end of the module (law 11d).
//
// `frames` does not appear: a create has no programme. The solver fields come back neutral. Allocates
// nothing itself and touches no handle.
fc_status fc_master_need_create (const fc_master_config* cfg, fc_need* out);

// The input's loudness range, for the LRA constraint — which is a DELTA and therefore needs both ends.
// Stateless by construction: it returns the number and the caller puts it into the request, so it
// cannot outlive the programme it describes.
//
// UNLIKE `process`, THIS REFUSES A POISONED PROGRAMME. A measurement has no gate to hide behind: one
// non-finite sample poisons the K-weighted state, its 10 ms is recorded as silence, and the meter goes
// on returning a plausible integrated number computed over part of the programme. The core makes that
// visible (`nonFiniteSubHops()`) and now refuses on it.
//
// BE EXACT ABOUT THE GRANULARITY, because the obvious claim is stronger than the truth: the counter is
// incremented when a 10 ms SUB-HOP COMPLETES, so a non-finite sample inside the final, incomplete
// sub-hop — under 480 frames at 48 kHz — is not counted and this call returns FC_OK. Closing that would
// mean scanning the input here, which is a second definition of "non-finite" and a second pass over the
// audio; the honest move is to say where the line is. Everything from one completed sub-hop onward is
// refused, which is every case that can move the number this call exists to produce.
//
// ON A DELIVERING HANDLE IT CONVERTS FIRST and measures the programme at the delivery rate — the programme the
// search will meter, so the range constraint compares a programme with itself. `frames` is the input count, and
// the call is refused (FC_ERR_REFUSED_BY_CORE, having converted nothing) when the DELIVERED length is too short
// to have a range. A delivered length past INT_MAX is FC_ERR_RANGE. Converting rather than refusing is the
// cheaper road for a page: the alternative is a second, non-delivering handle — two of the eight slots.
fc_status fc_master_measure_lra (fc_master h, const float* in, uint32_t frames, double* out);

// BS.1770 CHANNEL WEIGHTS for the solver's meters. The standard weights Ls/Rs at 1.41 and EXCLUDES LFE,
// and the core's default of 1.0 everywhere is correct for mono and stereo and wrong for surround. This
// ABI accepts up to `kMaxChannels`, so without this a correct surround search could not be expressed
// through it — and a facade that is a NARROWER road than the C++ API breaks the thinness law just as a
// wider one would. The host-layout-to-role mapping stays outside, as the core says.
fc_status fc_master_set_channel_weight (fc_master h, int32_t channel, double weight);

// Run the loudness search over the whole programme. `in` and `out` are planar with stride `frames` and
// MUST NOT be the same buffer: a search reads its input again on every pass, and rendering in place
// would make every pass after the first read the previous pass's master.
//
// FC_OK means A VERDICT WAS OBTAINED. Whether the target was met is `fc_solution_summary::status`.
// The solution handle is independent of the chain handle and outlives it; destroy it separately.
//
// THAT INCLUDES A DEGENERATE LENGTH: `frames == 0` is handed to the solver, which answers
// `InvalidRequest`, rather than being refused here. Answering for the core on one input while
// forwarding it on every other is two policies for one question, and this file's rule is that it never
// guesses a reason the core has.
//
// A SOLVE THAT RAN LEAVES THE HANDLE UNUSABLE FOR STREAMING until the next `fc_master_configure` — see
// `fc_master_process`. The delivered audio is already in `out`; a caller that wants to stream with the
// gain the search found configures with it.
// `params` travels with the call rather than being taken from the handle, and that is the C++
// signature's own choice: `TargetLoudnessSolver::solve()` takes the parameter set BY VALUE and its
// header says why — reading it back from the chain would give the last APPLIED set (a quantum stale)
// and the caller's UNCLAMPED request, and a search that read its own actuator through either would be
// measuring a number it did not apply.
fc_status fc_master_solve (fc_master h, const fc_master_params* params, const fc_loudness_request* req,
                           const float* in, float* out, uint32_t frames, fc_solution* out_solution);

//==============================================================================
// THE DELIVERING HANDLE (v2) — `fc_master_config::deliveryRate != 0`. SRC FIRST: the programme is converted to
// the delivery rate as one whole-programme operation (`mastering::DeliveryConverter`), and everything that
// measures or limits then sees the rate that is delivered, with the dither still the last thing to touch it.
// The composition lives in the core (`mastering::DeliveredMastering`), not here: this file forwards to it, and
// `fcore_master selftest` compares the two bit for bit.
//
// Each of these answers FC_ERR_STATE on a handle that does not deliver (`deliveryRate == 0`).
//
// THE CHECK ORDER FOR THE TWO LENGTHS is the header's own (narrowing before spans), made explicit because there
// are two: `inFrames`, `outFrames`, or the delivered length the core computes, past INT_MAX -> FC_ERR_RANGE;
// only then `outFrames` != the delivered length -> FC_ERR_CAPACITY; then the spans. `in` has stride `inFrames`,
// `out` has stride `outFrames`, and the two spans may not touch at all.
//
// A NON-FINITE INPUT SAMPLE is gated ahead of the conversion by the chain's own rule (NaN/inf -> 0, clamp +-1e6) and
// counted (`fc_master_stats::nonFiniteIn`): converting a NaN is bit-identical to converting a zero in its place, as
// processing one is on a plain handle. Without that a windowed sinc would spread one bad sample over its kernel.

// The exact delivered length for `inFrames` of programme — `DeliveryConverter::deliveredFrames`, read out, never
// recomputed: the natural `ceil(inFrames * deliveryRate / sampleRate)` in double is wrong on real lengths (147
// frames at 44.1 -> 48 kHz is 160, and the double says 161). `*outFrames` is written only on FC_OK.
fc_status fc_master_delivered_frames (fc_master h, uint32_t inFrames, uint32_t* outFrames);

// A render at the parameters of the last `fc_master_configure`, over the whole programme. Asks the heap for
// nothing: it converts into `out` and renders there. Refused with FC_ERR_STATE after a solve, until a configure,
// for the reason given at `fc_master_process`.
fc_status fc_master_render_delivered (fc_master h, const float* in, uint32_t inFrames, float* out, uint32_t outFrames);

// `fc_master_solve` over the delivered programme. Everything said there holds — FC_OK is a verdict, `params`
// travel with the call, the handle is left unusable for rendering until the next configure, a zero-length
// programme is forwarded to the solver for its verdict. The budget is FC_NEED_SOLVE on this handle.
fc_status fc_master_solve_delivered (fc_master h, const fc_master_params* params, const fc_loudness_request* req,
                                     const float* in, uint32_t inFrames, float* out, uint32_t outFrames,
                                     fc_solution* out_solution);

fc_status fc_solution_summary_get (fc_solution s, fc_solution_summary* out);
fc_status fc_solution_measurement (fc_solution s, fc_measurement* out);
// Copies min(logCount, cap) pass records into `out` and reports how many were written. Same ownership
// rule as everywhere else here: the buffer is the caller's and its capacity is binding.
fc_status fc_solution_log (fc_solution s, fc_solve_pass* out, uint32_t cap, uint32_t* written);
fc_status fc_solution_destroy (fc_solution s);

// THE CORE'S OWN DEFAULTS, written through the same mapping every other value crosses by.
//
// A zeroed struct is NOT a valid parameter set and is not close to one: `eq::BandParams` has its Stereo
// lane ON, `Dither` defaults to 24 bits with a non-zero seed and weighted shaping, and
// `MonoBassParams::enabled` defaults true while the CONFIG's `monoBass` defaults false. Measured, a
// memset(0) parameter set against `MasteringChainParams{}` differs in 287 998 of 288 000 samples, worst
// 0.874 full scale. So a caller starts from these and overwrites what it means to change; a caller that
// starts from zeroed memory is rendering something nobody chose.
//
// The FROZEN writers stamp the header themselves (at v1), so their result is immediately usable as an argument;
// the VERSIONED ones leave the caller's stamp as it is. NB both config writers leave `sampleRate` and `channels`
// at ZERO on purpose: the core has no
// default for either, so writing one would be this file choosing a geometry for every caller who forgot
// to state one. `fc_master_create` refuses both, which is how the caller finds out.
//
// VERSIONED: `fc_*_defaults` require `out` to carry a stamped header of any version this build knows (FC_INIT,
// or `Struct.init()` in JS), leave that header as it is, and write exactly that version's size — the only
// writers a caller newer than v1 should use. Same refusals as any OUT struct (FC_ERR_NULL, _ALIGNMENT, _SPAN,
// _ABI_VERSION, _STRUCT_SIZE) and guarded like every status-returning entry point.
fc_status fc_master_params_defaults (fc_master_params* out);
fc_status fc_master_config_defaults (fc_master_config* out);
fc_status fc_loudness_request_defaults (fc_loudness_request* out);

// FROZEN AT v1 (rule 8 of VERSIONING): these stamp v1 and write exactly v1's bytes, however large the struct the
// caller compiled against has since become. They exist for callers built against v1 and nothing else — a v2
// field set after one of these lies past the stamped size and is not read.
void fc_master_params_default (fc_master_params* out);
void fc_master_config_default (fc_master_config* out);
void fc_loudness_request_default (fc_loudness_request* out);

// Build identity, so a mismatched artifact is obvious in a report rather than a mystery. These take no
// handle and cannot fail.
uint32_t fc_master_abi_version (void);
uint32_t fc_master_max_channels (void);
uint32_t fc_master_max_eq_bands (void);
// The size table (rule 5): the bytes of struct `id` at `version`, 0 for a pair this build does not have — an
// unknown id, a version it does not know, or a struct that did not exist yet at that version.
uint32_t fc_master_sizeof (int32_t id, uint32_t version);
// FROZEN AT v1, like the `_default` writers: the sizes a v1 loader compares against. A v1 page then fails on
// the version check, which says what is wrong, instead of on a size, which would blame its own layout file.
uint32_t fc_master_sizeof_params (void);
uint32_t fc_master_sizeof_config (void);

#ifdef __cplusplus
}   // extern "C"
#endif

#endif   // FC_MASTER_ABI_H
