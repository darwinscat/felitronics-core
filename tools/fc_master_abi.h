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
// EXACT size match, not ">=", for v1. A tolerant read is a promise about how future fields will be
// laid out, and there is no future field yet to test that promise against; the first version to add
// one states its own rule then. Refusing early costs a caller a rebuild and costs nobody a wrong render.
//
// A NEW CODE IS NOT A NEW VERSION, and the rule for codes is written here rather than left to be inferred
// from the one for structs. A status or op code is only ever APPENDED — an existing code never changes
// meaning or value — so a caller built against an older header meets a code it does not know exactly
// where it already has to handle "not FC_OK", and nothing it relied on moved. The version moves only
// when a struct's layout does. (FC_ERR_POISONED and the FC_NEED_* ops came in this way.)
#define FC_MASTER_ABI_VERSION 1u

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
// the `*_default` writers and the build-identity queries — read no instance state and stay callable.
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
// The op is a code rather than a field per call so that the calls not yet budgeted here (create, configure —
// the chain's own storage) can join without moving this struct.
typedef enum fc_need_op { FC_NEED_SOLVE = 0, FC_NEED_MEASURE_LRA = 1 } fc_need_op;

typedef struct fc_need
{
    fc_header header;

    uint64_t callBytes;           // the core's PEAK request during ONE such call. A solve builds its meters per pass
                                  // and frees them at the pass's end, so this is one pass — plus the true-peak
                                  // drain, which the first solve keeps (so afterwards an upper bound by that drain)
    uint64_t solverPrepareBytes;  // the search's one-time preparation — tap buffers and gain-reduction histograms —
                                  // done lazily by the first solve, measure_lra or channel weight on the handle
    uint64_t facadeBytes;         // this file's own object for the call: a solve's solution record; 0 for measure_lra
    int32_t  solverPrepared;      // 1 once that preparation has happened: `solverPrepareBytes` is then already spent
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
// makes N calls identical to the last one alone. It costs a reallocation (the EQ engine is 331 KiB) and
// it is NOT real-time — it is a worker call between renders, and it says so.
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
// FC_ERR_ENUM — narrowing before field values, as everywhere. Reads the handle and moves nothing.
fc_status fc_master_need (fc_master h, int32_t op, uint32_t frames, fc_need* out);

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
// `header` is filled in too, so the result is immediately usable as an argument. NB
// `fc_master_config_default` leaves `sampleRate` and `channels` at ZERO on purpose: the core has no
// default for either, so writing one would be this file choosing a geometry for every caller who forgot
// to state one. `fc_master_create` refuses both, which is how the caller finds out.
void fc_master_params_default (fc_master_params* out);
void fc_master_config_default (fc_master_config* out);
void fc_loudness_request_default (fc_loudness_request* out);

// Build identity, so a mismatched artifact is obvious in a report rather than a mystery. These take no
// handle and cannot fail.
uint32_t fc_master_abi_version (void);
uint32_t fc_master_max_channels (void);
uint32_t fc_master_max_eq_bands (void);
uint32_t fc_master_sizeof_params (void);
uint32_t fc_master_sizeof_config (void);

#ifdef __cplusplus
}   // extern "C"
#endif

#endif   // FC_MASTER_ABI_H
