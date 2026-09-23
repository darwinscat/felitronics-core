// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// `FC_DOMAINS` of tools/wasm/fc-master-layout.mjs against the RUNNING ABI, one real call per bound:
//
//   felitronics_master_domains_tests <path to fc-master-layout.mjs>
//
// The table is the ARGUMENT, not a copy carried here, and every probe is DERIVED from its own `min` and `max`:
// the values tried are those bounds and one step either side. A bound that no longer matches the code therefore
// lands on the wrong side of the real boundary. The path is required; without it the suite refuses to run.
//
// WHAT EACH `edge` IS CHECKED WITH:
//   'refuse'  the bound is accepted; one step past it AND well past it are the named FC_ERR_*. Both, because a
//             bound moved inwards still refuses one step past itself — the real boundary is further out.
//   'clamp'   both are ACCEPTED, and where `resolved` names a read-back, the value applied AT the bound and
//             BEYOND it is the same number while one step INSIDE it is a different one.
//   'verdict' `fc_master_solve` answers FC_OK either way; the SUMMARY carries FC_SOLVE_INVALID_REQUEST outside.
//   'free'    a huge finite value of the field's own width is accepted, and the row carries no bound.
//   'any'     0, 1, -1 and a huge value are accepted, the field is an integer, and the row carries no bound.
// `nonFinite` is checked on every row, and 'none' is a claim like any other: it asserts the field is an integer.
//
// TWO SAMPLE RATES for every row whose bound is written in `sr`.
//
// WHAT IT DOES NOT CATCH: a row DEMOTED to a weaker claim — `clamp` rewritten as `free` or `any` with the
// bounds dropped — since an unbounded row asserts only that huge values are accepted, which a clamped field
// also does. The structural half of that class is closed (a floating field may not be `any`, an integer may not
// be `free`, neither may carry a bound). What is left cannot be closed by testing the claim, because a field
// the code does not bound and a field whose EFFECT saturates look the same from outside.
//
// AND WHAT NOTHING HERE CAN PIN is printed by name on every run: a clamp with no read-back can be shown to be
// a clamp rather than a refusal, and no more.

#include <felitronics_test.h>

#include "fc_master_abi.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using felitronics::test::ok;
using felitronics::test::group;

namespace
{

constexpr double        kFsA         = 48000.0;    // every row
constexpr double        kFsB         = 96000.0;    // ... and the rate-dependent rows again here
constexpr std::uint32_t kSolveFrames = 4800;       // 0.1 s: `admits()` runs before the first pass, so the
                                                   // programme only has to exist

//==============================================================================
// THE TABLE, READ OUT OF THE JAVASCRIPT

// A bound: a number, nothing at all, or `<a>*sr` / `<a>/sr` / `<a>*sr*os` — the four forms `domainBound()`
// evaluates, `sr` being the chain rate and `os` the oversample factor.
struct Bound
{
    bool   present = false;
    double a       = 0.0;
    char   op      = 0;            // 0 = a plain number · '*' = a*sr · '/' = a/sr · 'o' = a*sr*os
    double at (double sr, double os) const noexcept
    {
        return op == '*' ? a * sr : op == '/' ? a / sr : op == 'o' ? a * sr * os : a;
    }
    bool scaled() const noexcept { return op != 0; }
};

struct Row
{
    std::string field, unit, open, edge, err, nonFinite, resolved;
    Bound       min, max;
};

bool parseBound (std::string t, Bound& out)
{
    while (! t.empty() && (t.front() == ' ' || t.front() == '\t')) t.erase (t.begin());
    while (! t.empty() && (t.back()  == ' ' || t.back()  == '\t')) t.pop_back();
    if (t == "null") { out = Bound {}; return true; }
    if (t.size() >= 2 && t.front() == '\'' && t.back() == '\'')
    {
        static const std::regex re (R"(^(-?[0-9.eE+-]+)(\*sr\*os|\*sr|/sr)$)");
        std::smatch m;
        const std::string body = t.substr (1, t.size() - 2);
        if (! std::regex_match (body, m, re)) return false;
        const std::string tail = m[2].str();
        out.present = true;
        out.a       = std::strtod (m[1].str().c_str(), nullptr);
        out.op      = tail == "/sr" ? '/' : tail == "*sr" ? '*' : 'o';
        return true;
    }
    char* end = nullptr;
    const double v = std::strtod (t.c_str(), &end);
    if (end == t.c_str() || *end != '\0') return false;
    out.present = true; out.a = v; out.op = 0;
    return true;
}

// The block is REQUIRED to be found, and every `{ field:` inside it REQUIRED to parse: a partial read would
// leave rows unchecked and the run green.
bool readDomains (const char* path, std::vector<Row>& rows, std::string& why)
{
    std::ifstream f (path);
    if (! f) { why = "cannot open " + std::string (path); return false; }
    std::stringstream ss; ss << f.rdbuf();
    const std::string text = ss.str();

    const std::size_t from = text.find ("export const FC_DOMAINS = [");
    if (from == std::string::npos) { why = "no `export const FC_DOMAINS = [` in the file"; return false; }
    const std::size_t to = text.find ("\n];", from);
    if (to == std::string::npos) { why = "the FC_DOMAINS block is not terminated"; return false; }
    const std::string block = text.substr (from, to - from);

    static const std::regex re (
        R"(\{\s*field:\s*'([^']*)',\s*unit:\s*'([^']*)',\s*min:\s*([^,]+),\s*max:\s*([^,]+),)"
        R"(\s*open:\s*'([^']*)',\s*edge:\s*'([^']*)',\s*err:\s*'([^']*)',)"
        R"(\s*nonFinite:\s*'([^']*)',\s*resolved:\s*'([^']*)',\s*depends:)");

    std::size_t declared = 0;
    for (std::size_t i = block.find ("{ field:"); i != std::string::npos; i = block.find ("{ field:", i + 1)) ++declared;

    for (auto it = std::sregex_iterator (block.begin(), block.end(), re); it != std::sregex_iterator(); ++it)
    {
        Row r;
        r.field = (*it)[1]; r.unit = (*it)[2];
        if (! parseBound ((*it)[3], r.min)) { why = r.field + ": `min` is not a bound"; return false; }
        if (! parseBound ((*it)[4], r.max)) { why = r.field + ": `max` is not a bound"; return false; }
        r.open = (*it)[5]; r.edge = (*it)[6]; r.err = (*it)[7];
        r.nonFinite = (*it)[8]; r.resolved = (*it)[9];
        rows.push_back (r);
    }
    if (rows.size() != declared)
    {
        why = "parsed " + std::to_string (rows.size()) + " of " + std::to_string (declared) + " rows";
        return false;
    }
    return ! rows.empty();
}

// An `err` name back to the code it names. Unknown is a failure, never a skip.
bool statusNamed (const std::string& name, fc_status& out)
{
    static const struct { const char* n; fc_status s; } kMap[] = {
        { "FC_OK", FC_OK }, { "FC_ERR_HANDLE", FC_ERR_HANDLE }, { "FC_ERR_ABI_VERSION", FC_ERR_ABI_VERSION },
        { "FC_ERR_STRUCT_SIZE", FC_ERR_STRUCT_SIZE }, { "FC_ERR_NULL", FC_ERR_NULL },
        { "FC_ERR_ALIGNMENT", FC_ERR_ALIGNMENT }, { "FC_ERR_SPAN", FC_ERR_SPAN }, { "FC_ERR_ENUM", FC_ERR_ENUM },
        { "FC_ERR_RANGE", FC_ERR_RANGE }, { "FC_ERR_CAPACITY", FC_ERR_CAPACITY }, { "FC_ERR_STATE", FC_ERR_STATE },
        { "FC_ERR_NON_FINITE", FC_ERR_NON_FINITE }, { "FC_ERR_REFUSED_BY_CORE", FC_ERR_REFUSED_BY_CORE },
        { "FC_ERR_EXHAUSTED", FC_ERR_EXHAUSTED }, { "FC_ERR_POISONED", FC_ERR_POISONED },
        { "FC_ERR_CANCELLED", FC_ERR_CANCELLED },
    };
    for (const auto& e : kMap) if (name == e.n) { out = e.s; return true; }
    return false;
}

const char* statusName (fc_status s)
{
    static const char* n[] = { "FC_OK", "FC_ERR_HANDLE", "FC_ERR_ABI_VERSION", "FC_ERR_STRUCT_SIZE", "FC_ERR_NULL",
                               "FC_ERR_ALIGNMENT", "FC_ERR_SPAN", "FC_ERR_ENUM", "FC_ERR_RANGE", "FC_ERR_CAPACITY",
                               "FC_ERR_STATE", "FC_ERR_NON_FINITE", "FC_ERR_REFUSED_BY_CORE", "FC_ERR_EXHAUSTED",
                               "FC_ERR_POISONED", "FC_ERR_CANCELLED" };
    return (s >= 0 && s <= FC_ERR_CANCELLED) ? n[(int) s] : "FC_STATUS(?)";
}

//==============================================================================
// THE KNOBS — one writer per row, and the row's own path is the key

enum class Where { Config, Params, Request };
// F64W is a double that must carry a WHOLE number — `deliveryRate`, whose plan refuses a fraction before it
// looks at the range, so a probe off a bound by one per cent would be refused for the wrong reason.
enum class Ty    { F64, F64W, F32, I32, U32 };

struct World
{
    fc_master_config    cfg {};
    fc_master_params    prm {};
    fc_loudness_request req {};
};

// `stepRel`/`stepAbs` are the RESOLUTION OF THIS FIELD'S READ-BACK — the smallest change at a bound the
// read-back can still tell apart, and therefore how tightly the bound is pinned. 0 takes the default for the
// field's type and the row's `edge`: a status boundary is sharp, a read-back is not. `arm` is the parameter
// set under which the bound is observable at all, applied BEFORE the probed field is written; its argument
// says which bound is being approached, for the fields whose two ends need different signals.
struct Knob
{
    const char* field;
    Where       where;
    Ty          type;
    void      (*set) (World&, double);
    double      stepRel = 0.0;
    double      stepAbs = 0.0;
    void      (*arm) (World&, bool isMin) = nullptr;
};

// A double into an integer field WITHOUT the out-of-range conversion, which is undefined and is exactly what a
// probe at 1e300 would reach.
std::int32_t toI32 (double v) noexcept
{
    if (! (v > -2147483649.0)) return (std::int32_t) -2147483647 - 1;
    if (! (v <  2147483648.0)) return (std::int32_t)  2147483647;
    return (std::int32_t) v;
}
std::uint32_t toU32 (double v) noexcept
{
    if (! (v > -2147483649.0)) return 0x80000000u;
    if (! (v <  4294967296.0)) return 0xFFFFFFFFu;
    return v < 0.0 ? (std::uint32_t) (std::int64_t) v : (std::uint32_t) v;
}

#define CFG_D(path, member) { path, Where::Config,  Ty::F64, [] (World& w, double v) { w.cfg.member = v; } }
#define CFG_I(path, member) { path, Where::Config,  Ty::I32, [] (World& w, double v) { w.cfg.member = toI32 (v); } }
#define PRM_D(path, member) { path, Where::Params,  Ty::F64, [] (World& w, double v) { w.prm.member = v; } }
#define PRM_F(path, member) { path, Where::Params,  Ty::F32, [] (World& w, double v) { w.prm.member = (float) v; } }
#define PRM_I(path, member) { path, Where::Params,  Ty::I32, [] (World& w, double v) { w.prm.member = toI32 (v); } }
#define PRM_U(path, member) { path, Where::Params,  Ty::U32, [] (World& w, double v) { w.prm.member = toU32 (v); } }
#define REQ_D(path, member) { path, Where::Request, Ty::F64, [] (World& w, double v) { w.req.member = v; } }
#define REQ_I(path, member) { path, Where::Request, Ty::I32, [] (World& w, double v) { w.req.member = toI32 (v); } }

// The EQ fields are written into the FIRST and the LAST slot of every array they live in: `toCore` walks all
// bands and all lanes, and a probe confined to band 0 lane 0 would not see a walk that had stopped short.
#define EQB_D(path, member) { path, Where::Params, Ty::F64, [] (World& w, double v) \
    { w.prm.eqBands[0].member = v; w.prm.eqBands[FC_MAX_EQ_BANDS - 1].member = v; } }
#define EQB_I(path, member) { path, Where::Params, Ty::I32, [] (World& w, double v) \
    { w.prm.eqBands[0].member = toI32 (v); w.prm.eqBands[FC_MAX_EQ_BANDS - 1].member = toI32 (v); } }
#define EQL_D(path, member) { path, Where::Params, Ty::F64, [] (World& w, double v) \
    { w.prm.eqBands[0].lanes[0].member = v; w.prm.eqBands[FC_MAX_EQ_BANDS - 1].lanes[FC_MAX_EQ_LANES - 1].member = v; } }
#define EQL_I(path, member) { path, Where::Params, Ty::I32, [] (World& w, double v) \
    { w.prm.eqBands[0].lanes[0].member = toI32 (v); w.prm.eqBands[FC_MAX_EQ_BANDS - 1].lanes[FC_MAX_EQ_LANES - 1].member = toI32 (v); } }

// The same three, with a probe resolution and an arming.
#define PRM_DX(path, member, rel, ab, armfn) { path, Where::Params, Ty::F64, \
    [] (World& w, double v) { w.prm.member = v; }, rel, ab, armfn }
#define EQB_DX(path, member, rel, ab, armfn) { path, Where::Params, Ty::F64, \
    [] (World& w, double v) { w.prm.eqBands[0].member = v; w.prm.eqBands[FC_MAX_EQ_BANDS - 1].member = v; }, \
    rel, ab, armfn }
#define REQ_DX(path, member, rel, ab, armfn) { path, Where::Request, Ty::F64, \
    [] (World& w, double v) { w.req.member = v; }, rel, ab, armfn }

// THE ARMINGS. Each is a parameter set under which one field's bound is REACHABLE: an idle stage, a delta
// already hard against its cap, or a quantiser whose step is coarser than the probe all make a clamp
// invisible, and an invisible clamp cannot be pinned to the place the table puts it.
void armNoDither (World& w, bool) { w.cfg.dither = 0; }

// The compressor's range caps the DELTA, so it bites only where the delta would exceed it.
void armCompCap (World& w, bool)
{
    w.cfg.dither = 0;
    w.prm.compressor.thresholdDb = -2000.0;
}

// The search moves the gain after its first render, so what it reports is the gain it was GIVEN only while
// there is one pass to give it to.
void armOnePass (World& w, bool) { w.req.maxPasses = 1; }

// THE POINT'S DYNAMICS, ARMED AND ABSOLUTE. `rangeDb` caps the delta, so a clamp on it is visible only
// where the reduction is hard against that cap — which a threshold every sample is over is what buys.
// The same arming makes the delta travel from 0 to that cap at the start of the programme, which is
// what a clamp on `atk`/`rel` needs to be seen at all.
void armDyn (World& w, bool)
{
    w.cfg.dither = 0;
    w.prm.eqBands[0].dyn.on      = 1;
    w.prm.eqBands[0].dyn.thrAuto = 0;
    w.prm.eqBands[0].dyn.thrDb   = -120.0;
    w.prm.eqBands[0].dyn.rangeDb = -30.0;
}

// THE THRESHOLD is measured against the band probe's own level, so its two bounds are reachable only
// where the programme sits near them — and between them the computer must be OFF its cap, since a
// reduction pinned there is the same reduction at every threshold. The input gain is what moves the
// programme to meet each end.
void armDynThr (World& w, bool isMin)
{
    armDyn (w, isMin);
    w.prm.inputGainDb = isMin ? -60.0 : 60.0;
}

// A BALLISTICS KNOB CHANGES NOTHING UNLESS THE DELTA MOVES BOTH WAYS, and the arming above only ever
// moves it one: with the threshold 96 dB under this programme the computer sits on its cap and the
// reduction never releases (measured over the 18 quanta of the fixture: eighteen falls, not one rise).
// A threshold inside the programme's own band level runs both halves.
void armDynMove (World& w, bool isMin)
{
    armDyn (w, isMin);
    w.prm.eqBands[0].dyn.thrDb = -60.0;
}

const Knob kKnobs[] = {
    CFG_D ("fc_master_config.sampleRate", sampleRate),
    // THE MONO-BASS STAGE IS TURNED OFF BY THIS ONE WRITER, unconditionally: the stage admits a width of 2 and
    // nothing else, which is an interaction with its own check below rather than part of this row.
    { "fc_master_config.channels", Where::Config, Ty::I32,
      [] (World& w, double v) { w.cfg.channels = toI32 (v); w.cfg.monoBass = 0; } },
    CFG_I ("fc_master_config.internalBlock", internalBlock),
    CFG_I ("fc_master_config.eq", eq),
    CFG_I ("fc_master_config.monoBass", monoBass),
    CFG_I ("fc_master_config.compressor", compressor),
    CFG_I ("fc_master_config.clipper", clipper),
    CFG_I ("fc_master_config.limiter", limiter),
    CFG_I ("fc_master_config.dither", dither),
    CFG_D ("fc_master_config.compressorLookaheadMs", compressorLookaheadMs),
    CFG_D ("fc_master_config.limiterLookaheadMs", limiterLookaheadMs),
    CFG_I ("fc_master_config.oversampleFactor", oversampleFactor),
    CFG_I ("fc_master_config.tapsPerPhase", tapsPerPhase),
    CFG_D ("fc_master_config.sidechainHpfHz", sidechainHpfHz),
    { "fc_master_config.deliveryRate", Where::Config, Ty::F64W,
      [] (World& w, double v) { w.cfg.deliveryRate = std::round (v); } },   // whole hertz, or no route at all

    PRM_DX ("fc_master_params.inputGainDb", inputGainDb, 1.0e-6, 1.0e-6, armNoDither),
    PRM_DX ("fc_master_params.preLimiterGainDb", preLimiterGainDb, 1.0e-6, 1.0e-6, armNoDither),
    PRM_DX ("fc_master_params.compressorMix", compressorMix, 1.0e-6, 1.0e-6, nullptr),
    PRM_I ("fc_master_params.limiterDualRelease", limiterDualRelease),
    PRM_DX ("fc_master_params.limiterSlowReleaseMs", limiterSlowReleaseMs, 1.0e-4, 1.0e-4, nullptr),
    PRM_I  ("fc_master_params.peakClipper", peakClipper),
    PRM_DX ("fc_master_params.peakClipperOverCeilingDb", peakClipperOverCeilingDb, 1.0e-9, 1.0e-9, nullptr),
    PRM_DX ("fc_master_params.peakClipperKneeDb", peakClipperKneeDb, 1.0e-9, 1.0e-9, nullptr),
    PRM_I ("fc_master_params.bypassEq", bypassEq),
    PRM_I ("fc_master_params.bypassMonoBass", bypassMonoBass),
    PRM_I ("fc_master_params.bypassCompressor", bypassCompressor),
    PRM_I ("fc_master_params.bypassClipper", bypassClipper),
    PRM_I ("fc_master_params.bypassLimiter", bypassLimiter),
    PRM_I ("fc_master_params.bypassDither", bypassDither),

    EQB_I ("fc_master_params.eqBands[].on", on),
    EQB_I ("fc_master_params.eqBands[].type", type),
    EQB_I ("fc_master_params.eqBands[].swept", swept),
    EQB_I ("fc_master_params.eqBands[].bypass", bypass),
    EQB_I ("fc_master_params.eqBands[].dyn.on", dyn.on),
    // A RENDER IS THE READ-BACK for the four `dyn` clamps, and the probe resolution is stated because the
    // default billionth is below what a float coefficient can carry: the bell's gain, the computer's cap
    // and the follower's time constant all reach the samples through a narrowing.
    EQB_DX ("fc_master_params.eqBands[].dyn.rangeDb", dyn.rangeDb, 1.0e-4, 1.0e-4, armDyn),
    EQB_DX ("fc_master_params.eqBands[].dyn.thrDb", dyn.thrDb, 1.0e-4, 1.0e-4, armDynThr),
    EQB_I ("fc_master_params.eqBands[].dyn.thrAuto", dyn.thrAuto),
    EQB_DX ("fc_master_params.eqBands[].dyn.atk", dyn.atk, 1.0e-3, 1.0e-3, armDynMove),
    EQB_DX ("fc_master_params.eqBands[].dyn.rel", dyn.rel, 1.0e-3, 1.0e-3, armDynMove),
    EQL_I ("fc_master_params.eqBands[].lanes[].on", on),
    EQL_D ("fc_master_params.eqBands[].lanes[].freq", freq),
    EQL_D ("fc_master_params.eqBands[].lanes[].q", q),
    EQL_D ("fc_master_params.eqBands[].lanes[].gainDb", gainDb),
    EQL_I ("fc_master_params.eqBands[].lanes[].slope", slope),
    EQL_I ("fc_master_params.eqBands[].lanes[].bypass", bypass),

    PRM_I ("fc_master_params.monoBass.enabled", monoBass.enabled),
    PRM_F ("fc_master_params.monoBass.frequencyHz", monoBass.frequencyHz),
    PRM_F ("fc_master_params.monoBass.lowWidth", monoBass.lowWidth),

    PRM_I ("fc_master_params.compressor.detector", compressor.detector),
    PRM_I ("fc_master_params.compressor.link", compressor.link),
    PRM_DX ("fc_master_params.compressor.rmsWindowMs", compressor.rmsWindowMs, 1.0e-6, 1.0e-2, armNoDither),
    PRM_I ("fc_master_params.compressor.mode", compressor.mode),
    PRM_D ("fc_master_params.compressor.thresholdDb", compressor.thresholdDb),
    PRM_DX ("fc_master_params.compressor.ratio", compressor.ratio, 1.0e-5, 1.0e-5, armNoDither),
    PRM_DX ("fc_master_params.compressor.kneeDb", compressor.kneeDb, 1.0e-6, 1.0, armNoDither),
    PRM_DX ("fc_master_params.compressor.rangeDb", compressor.rangeDb, 1.0e-6, 1.0e-5, armCompCap),
    PRM_DX ("fc_master_params.compressor.attackMs", compressor.attackMs, 1.0e-6, 1.0e-2, armNoDither),
    PRM_DX ("fc_master_params.compressor.releaseMs", compressor.releaseMs, 1.0e-6, 1.0e-2, armNoDither),
    PRM_D ("fc_master_params.compressor.makeupDb", compressor.makeupDb),
    PRM_I ("fc_master_params.compressor.autoMakeup", compressor.autoMakeup),

    PRM_I ("fc_master_params.clipper.shape", clipper.shape),
    PRM_F ("fc_master_params.clipper.driveDb", clipper.driveDb),
    PRM_F ("fc_master_params.clipper.bias", clipper.bias),
    PRM_F ("fc_master_params.clipper.mix", clipper.mix),
    PRM_F ("fc_master_params.clipper.outputDb", clipper.outputDb),
    PRM_F ("fc_master_params.clipper.autoComp", clipper.autoComp),
    PRM_F ("fc_master_params.clipper.dcBlockHz", clipper.dcBlockHz),

    PRM_D ("fc_master_params.limiter.ceilingDbTp", limiter.ceilingDbTp),
    PRM_DX ("fc_master_params.limiter.releaseMs", limiter.releaseMs, 1.0e-4, 1.0e-4, nullptr),

    PRM_I ("fc_master_params.dither.bits", dither.bits),
    PRM_I ("fc_master_params.dither.shaping", dither.shaping),
    PRM_U ("fc_master_params.dither.seedLo", dither.seedLo),
    PRM_U ("fc_master_params.dither.seedHi", dither.seedHi),
    PRM_I ("fc_master_params.dither.autoBlank", dither.autoBlank),
    PRM_I ("fc_master_params.dither.autoBlankSamples", dither.autoBlankSamples),

    REQ_D ("fc_loudness_request.targetLufs", targetLufs),
    REQ_D ("fc_loudness_request.toleranceLu", toleranceLu),
    REQ_D ("fc_loudness_request.maxTruePeakDbTp", maxTruePeakDbTp),
    REQ_D ("fc_loudness_request.truePeakAimDb", truePeakAimDb),
    REQ_D ("fc_loudness_request.limiterGr.limitDb", limiterGr.limitDb),
    REQ_I ("fc_loudness_request.limiterGr.statistic", limiterGr.statistic),
    REQ_D ("fc_loudness_request.compressorGr.limitDb", compressorGr.limitDb),
    REQ_I ("fc_loudness_request.compressorGr.statistic", compressorGr.statistic),
    REQ_D ("fc_loudness_request.minPlrDb", minPlrDb),
    REQ_D ("fc_loudness_request.maxLraLossLu", maxLraLossLu),
    REQ_D ("fc_loudness_request.inputLoudnessRangeLu", inputLoudnessRangeLu),
    REQ_D ("fc_loudness_request.activityThresholdDb", activityThresholdDb),
    REQ_I ("fc_loudness_request.maxPasses", maxPasses),
    REQ_DX ("fc_loudness_request.initialGainDb", initialGainDb, 1.0e-9, 1.0e-9, armOnePass),
    REQ_I ("fc_loudness_request.grTraceBuckets", grTraceBuckets),
    REQ_D ("fc_loudness_request.limiterGrQuantile", limiterGrQuantile),
    REQ_D ("fc_loudness_request.compressorGrQuantile", compressorGrQuantile),
    // v10 / K11 — the gate on the limiter's input. It decides nothing the solver judges, so its row's `edge` is
    // `free`: there is no clamp and no refusal to probe, only the value crossing and being echoed back.
    REQ_D ("fc_loudness_request.limiterActiveInputDb", limiterActiveInputDb),
};

const Knob* knobFor (const std::string& field)
{
    for (const Knob& k : kKnobs) if (field == k.field) return &k;
    return nullptr;
}

//==============================================================================
// ONE PROBE — a real call, through the C entry points, with one field moved

// ONE GEOMETRY THAT LEAVES NOTHING IDLE — an idle stage makes a clamp on its own parameter invisible. Every
// stage is on; band 0 is audible and its dynamics armed; the compressor runs its RMS detector, the only setting
// under which `rmsWindowMs` is read; the clipper runs Asym, the only shape that reads `bias` and enables the DC
// blocker; the second limiter envelope is on, since its resolved slow release reads 0 while it is off.
World baseWorld (double sr)
{
    World w;
    FC_INIT (w.cfg); fc_master_config_defaults (&w.cfg);
    w.cfg.sampleRate = sr;
    w.cfg.channels   = 2;
    w.cfg.eq = w.cfg.monoBass = w.cfg.compressor = w.cfg.clipper = w.cfg.limiter = w.cfg.dither = 1;

    FC_INIT (w.prm); fc_master_params_defaults (&w.prm);
    w.prm.limiterDualRelease = 1;
    w.prm.eqBands[0].on = 1;
    w.prm.eqBands[0].type = 0;                      // FC_FILTER_BELL
    w.prm.eqBands[0].lanes[0].on     = 1;
    w.prm.eqBands[0].lanes[0].freq   = 1000.0;
    w.prm.eqBands[0].lanes[0].q      = 1.0;
    w.prm.eqBands[0].lanes[0].gainDb = 6.0;
    w.prm.eqBands[0].dyn.on      = 1;
    w.prm.eqBands[0].dyn.rangeDb = -6.0;
    w.prm.eqBands[0].dyn.thrDb   = -12.0;
    w.prm.compressor.detector = 1;                  // FC_DETECTOR_RMS
    w.prm.clipper.shape       = 3;                  // FC_SHAPE_ASYM
    w.prm.clipper.driveDb     = 6.0f;

    FC_INIT (w.req); fc_loudness_request_defaults (&w.req);
    w.req.targetLufs      = -14.0;
    w.req.maxTruePeakDbTp =  -1.0;
    w.req.maxPasses       =   2;
    return w;
}

// The magnitude of band 0 on the Stereo axis, on a log grid: the read-back for the EQ lane fields, which
// `fc_master_resolved` does not carry.
std::vector<double> eqCurve (const fc_master_params& p, double sr)
{
    constexpr std::uint32_t kN = 24;
    std::vector<double> f (kN), db (kN);
    const double lo = 10.0, hi = 0.49 * sr;
    for (std::uint32_t i = 0; i < kN; ++i) f[i] = lo * std::pow (hi / lo, (double) i / (double) (kN - 1));
    std::uint32_t written = 0;
    if (fc_master_eq_curve (&p, sr, 0, 0, f.data(), kN, db.data(), kN, &written) != FC_OK || written != kN)
        return {};
    return db;
}

bool summaryField (const fc_solution_summary& sum, const std::string& path, double& out)
{
    if (path == "preLimiterGainDb") { out = sum.preLimiterGainDb; return true; }
    if (path == "ceilingDbTp")      { out = sum.ceilingDbTp;      return true; }
    return false;
}

bool resolvedField (const fc_master_resolved& r, const std::string& path, double& out)
{
    if (path == "internalBlock")        { out = (double) r.internalBlock;    return true; }
    if (path == "oversampleFactor")     { out = (double) r.oversampleFactor; return true; }
    if (path == "limiterLookahead")     { out = (double) r.limiterLookahead; return true; }
    if (path == "limiterCeilingDbTp")   { out = r.limiterCeilingDbTp;        return true; }
    if (path == "limiterReleaseMs")     { out = r.limiterReleaseMs;          return true; }
    if (path == "limiterSlowReleaseMs") { out = r.limiterSlowReleaseMs;      return true; }
    if (path == "peakClipperThresholdDbTp") { out = r.peakClipperThresholdDbTp; return true; }
    if (path == "compressorMix")        { out = r.compressorMix;             return true; }
    if (path == "monoBass.enabled")     { out = (double) r.monoBass.enabled; return true; }
    if (path == "monoBass.frequencyHz") { out = (double) r.monoBass.frequencyHz; return true; }
    if (path == "monoBass.lowWidth")    { out = (double) r.monoBass.lowWidth;    return true; }
    return false;
}

// EXACTLY the same number, spelled without `==`. The comparison must be exact: a tolerance would pass a clamp
// that had moved by less than it.
bool sameValue (double a, double b) noexcept { return ! (a < b) && ! (b < a); }

bool sameCurve (const std::vector<double>& a, const std::vector<double>& b) noexcept
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (! sameValue (a[i], b[i])) return false;
    return true;
}

struct Result
{
    fc_status           status  = FC_OK;   // create's for a config row, configure's for a params row, solve's for a request
    int                 verdict = -1;      // the SOLUTION's status, or -1 where no solve was made
    bool                haveRead = false;
    std::vector<double> read;              // the applied value, or the EQ curve
    bool                readFailed = false;
};

// BROADBAND and deterministic, so every filter's passband and every detector's range are excited; one sequence
// per channel, so no fault in the planar addressing leaves the two planes agreeing. Built once; the caller
// supplies the destination, since a search may not render in place.
const std::vector<float>& programme()
{
    static const std::vector<float> in = [] {
        std::vector<float> v ((std::size_t) kSolveFrames * 2, 0.0f);
        std::uint32_t s0 = 0x2545F491u, s1 = 0x9E3779B9u;
        for (std::uint32_t i = 0; i < kSolveFrames; ++i)
        {
            s0 = s0 * 1664525u + 1013904223u;
            s1 = s1 * 1664525u + 1013904223u;
            v[i]                = 0.50f * ((float) (s0 >> 8) * (1.0f / 8388608.0f) - 1.0f);
            v[kSolveFrames + i] = 0.40f * ((float) (s1 >> 8) * (1.0f / 8388608.0f) - 1.0f);
        }
        return v;
    }();
    return in;
}

Result probe (const Knob& k, double value, double sr, const std::string& resolved, bool isMin = true)
{
    World w = baseWorld (sr);
    if (k.arm != nullptr) k.arm (w, isMin);
    k.set (w, value);

    Result out;
    fc_master h = 0;
    const fc_status cr = fc_master_create (&w.cfg, &h);
    if (k.where == Where::Config) out.status = cr;
    if (cr != FC_OK) { if (k.where != Where::Config) out.status = cr; return out; }

    fc_master_resolved r {}; FC_INIT (r);
    const fc_status cf = fc_master_configure (h, &w.prm, &r);
    if (k.where == Where::Params) out.status = cf;
    else if (k.where == Where::Config && cf != FC_OK) out.status = cf;   // an accepted geometry must configure

    if (cf == FC_OK && ! resolved.empty() && resolved.rfind ("summary.", 0) != 0)
    {
        if (resolved == "eqCurve") { out.read = eqCurve (w.prm, sr); out.haveRead = ! out.read.empty(); }
        else if (resolved == "render")
        {
            // WHAT THE CHAIN DOES WITH IT, for a clamp `fc_master_resolved` does not report. Same build, same
            // input, same parameters, same bits — so two renders that agree are two values the stage could
            // not tell apart.
            std::vector<float> dst ((std::size_t) kSolveFrames * 2, 0.0f);
            if (fc_master_process (h, programme().data(), dst.data(), kSolveFrames) == FC_OK)
            {
                out.read.assign (dst.begin(), dst.end());
                out.haveRead = true;
            }
        }
        else
        {
            double v = 0.0;
            if (resolvedField (r, resolved, v)) { out.read = { v }; out.haveRead = true; }
            else out.readFailed = true;
        }
    }

    if (cf == FC_OK && k.where == Where::Request)
    {
        std::vector<float> dst ((std::size_t) kSolveFrames * 2, 0.0f);
        fc_solution s = 0;
        out.status = fc_master_solve (h, &w.prm, &w.req, programme().data(), dst.data(), kSolveFrames, &s);
        if (out.status == FC_OK)
        {
            fc_solution_summary sum {}; FC_INIT (sum);
            if (fc_solution_summary_get (s, &sum) == FC_OK)
            {
                out.verdict = sum.status;
                if (resolved.rfind ("summary.", 0) == 0)
                {
                    double v = 0.0;
                    if (summaryField (sum, resolved.substr (8), v)) { out.read = { v }; out.haveRead = true; }
                    else out.readFailed = true;
                }
            }
            fc_solution_destroy (s);
        }
    }
    fc_master_destroy (h);
    return out;
}

//==============================================================================
// THE ROW CHECKS

constexpr int kInvalidRequest = 8;      // FC_SOLVE_INVALID_REQUEST

std::string at (const Row& row, double sr, const std::string& what, double v)
{
    char buf[320];
    std::snprintf (buf, sizeof buf, "%s @ %g Hz: %s = %.10g", row.field.c_str(), sr, what.c_str(), v);
    return buf;
}

// ONE STEP OUTSIDE A BOUND IS THE RESOLUTION OF WHAT ANSWERS, and the bound is pinned to within it. A whole
// number steps by 1. A status is sharp, so a refusal steps by a billionth. A read-back is not, and the field
// says how coarse its own is (`stepRel`/`stepAbs`); the default is a float ulp's worth for an `f32` field and
// a billionth for an `f64` one.
double stepFor (const Knob& k, const Row& row, double bound)
{
    if (k.type == Ty::I32 || k.type == Ty::U32 || k.type == Ty::F64W) return 1.0;
    const bool   sharp = (row.edge == "refuse" || row.edge == "verdict");
    const double tiny  = (! sharp && k.type == Ty::F32) ? 1.0e-6 : 1.0e-9;
    const double rel   = k.stepRel > 0.0 ? k.stepRel : tiny;
    const double abs   = k.stepAbs > 0.0 ? k.stepAbs : tiny;
    const double sc    = std::fabs (bound) * rel;
    return sc > abs ? sc : abs;
}

// Was this value ADMITTED? For a config or a params row that is the status; for a request row the solver's
// verdict is the answer and the status is FC_OK either way, which is the trap this column exists to name.
bool admitted (const Row& row, const Result& res)
{
    if (row.field.rfind ("fc_loudness_request.", 0) == 0 && res.verdict >= 0)
        return res.status == FC_OK && res.verdict != kInvalidRequest;
    return res.status == FC_OK;
}

void checkBound (const Row& row, const Knob& k, double sr, bool isMin)
{
    const Bound& b = isMin ? row.min : row.max;
    if (! b.present) return;
    const double bound  = b.at (sr, (double) baseWorld (sr).cfg.oversampleFactor);
    const double step   = stepFor (k, row, bound);
    const double out    = isMin ? bound - step : bound + step;      // one step past it
    const double in     = isMin ? bound + step : bound - step;      // one step inside it
    const bool   isOpen = (row.open == (isMin ? "min" : "max")) || row.open == "both";

    // With an OPEN bound the bound itself is already outside, and the first admitted value is one step in.
    const double outside = isOpen ? bound : out;
    const double inside  = isOpen ? in    : bound;

    if (row.edge == "refuse" || row.edge == "verdict")
    {
        // ONE STEP OUT IS NOT ENOUGH ON ITS OWN: a bound moved INWARDS still refuses one step past itself,
        // the real boundary being further out. This second probe sits WELL outside — half a positive bound,
        // ten steps past a zero or negative one — where a bound that had crept inwards would have to ACCEPT.
        const double far = isMin ? (bound > 0.0 ? bound * 0.5 : bound - 10.0 * step)
                                 : (bound > 0.0 ? bound * 2.0 : bound + 10.0 * step);
        const Result bad = probe (k, outside, sr, "", isMin);
        const Result away = probe (k, far, sr, "", isMin);
        const Result good = probe (k, inside, sr, row.resolved, isMin);
        if (row.edge == "refuse")
        {
            fc_status want = FC_OK;
            ok (statusNamed (row.err, want), row.field + ": `err` names no status (" + row.err + ")");
            ok (bad.status == want, at (row, sr, "past the bound answers", outside)
                                    + " — " + statusName (bad.status) + ", not " + row.err);
            ok (away.status == want, at (row, sr, "well past the bound answers", far)
                                     + " — " + statusName (away.status) + ", not " + row.err);
        }
        else
        {
            ok (row.err == "FC_SOLVE_INVALID_REQUEST", row.field + ": a `verdict` row must name FC_SOLVE_INVALID_REQUEST");
            ok (bad.status == FC_OK, at (row, sr, "past the bound keeps the STATUS at FC_OK", outside));
            ok (bad.verdict == kInvalidRequest, at (row, sr, "past the bound is InvalidRequest", outside)
                                                + " — verdict " + std::to_string (bad.verdict));
            ok (away.verdict == kInvalidRequest, at (row, sr, "well past the bound is InvalidRequest", far)
                                                 + " — verdict " + std::to_string (away.verdict));
        }
        ok (admitted (row, good), at (row, sr, "the bound itself is admitted", inside)
                                  + " — " + statusName (good.status) + "/verdict " + std::to_string (good.verdict));
        // A whole-number field with a published read-back must read back as asked.
        if (good.haveRead && good.read.size() == 1 && (k.type == Ty::I32 || k.type == Ty::U32))
            ok (sameValue (good.read[0], inside), at (row, sr, "reads back as asked through resolved." + row.resolved, inside)
                                        + " — got " + std::to_string (good.read[0]));
        ok (! good.readFailed, row.field + ": `resolved` names no field of fc_master_resolved (" + row.resolved + ")");
        return;
    }

    if (row.edge == "clamp")
    {
        const Result beyond = probe (k, outside, sr, row.resolved, isMin);
        const Result atB    = probe (k, inside,  sr, row.resolved, isMin);
        const Result within = probe (k, isMin ? inside + step : inside - step, sr, row.resolved, isMin);
        ok (admitted (row, beyond), at (row, sr, "past the bound is ACCEPTED, not refused", outside)
                                    + " — " + statusName (beyond.status));
        ok (admitted (row, atB), at (row, sr, "the bound itself is accepted", inside));
        ok (! atB.readFailed, row.field + ": `resolved` names no field of fc_master_resolved (" + row.resolved + ")");
        if (row.resolved.empty()) return;
        ok (beyond.haveRead && atB.haveRead && within.haveRead, row.field + ": the read-back could not be taken");
        if (! (beyond.haveRead && atB.haveRead && within.haveRead)) return;
        ok (sameCurve (beyond.read, atB.read), at (row, sr, "is CLAMPED at the bound (the same value applies beyond it)", outside));
        ok (! sameCurve (within.read, atB.read), at (row, sr, "the clamp is no further in than the table says", inside));
        return;
    }
}

bool isFloating (const Knob& k) noexcept { return k.type != Ty::I32 && k.type != Ty::U32; }

void checkRow (const Row& row, const Knob& k, double sr)
{
    const double big = k.type == Ty::F32 ? 1.0e30 : isFloating (k) ? 1.0e300 : 1073741824.0;

    if (row.edge == "refuse" || row.edge == "verdict" || row.edge == "clamp")
    {
        checkBound (row, k, sr, true);
        checkBound (row, k, sr, false);
    }
    else if (row.edge == "free")
    {
        ok (! row.min.present && ! row.max.present, row.field + ": a `free` row may not carry a bound");
        // AND IT HAS TO BE A FLOATING FIELD: an integer the code does not bound is `any`.
        ok (isFloating (k), row.field + ": `free` is for a floating field; an unbounded integer is `any`");
        for (const double v : { big, -big })
            ok (admitted (row, probe (k, v, sr, "")), at (row, sr, "a huge finite value is admitted", v));
    }
    else if (row.edge == "any")
    {
        // `any` is honest only for a field whose every INTEGER value is a value — a flag read as `!= 0`, an
        // opaque code — so a floating field may not claim it, and the row may carry no bound.
        ok (! isFloating (k), row.field + ": a floating field cannot be `any` — it has a clamp or a refusal");
        ok (! row.min.present && ! row.max.present, row.field + ": an `any` row may not carry a bound");
        for (const double v : { 0.0, 1.0, -1.0, big })
            ok (admitted (row, probe (k, v, sr, "")), at (row, sr, "every value is admitted", v));
    }
    else
    {
        ok (false, row.field + ": unknown edge `" + row.edge + "`");
    }

    // ── the non-finite column ─────────────────────────────────────────────────────────────────────
    // 'none' IS A CLAIM, not a way out: it says the field is an integer and has no non-finite value to answer
    // for, and is checked as one.
    ok (isFloating (k) != (row.nonFinite == "none"),
        row.field + (isFloating (k) ? ": a floating field may not answer `none` for a non-finite value"
                                    : ": an integer field has no non-finite value and must answer `none`"));

    const double qnan = std::numeric_limits<double>::quiet_NaN();
    const double posInf = std::numeric_limits<double>::infinity();
    if (row.nonFinite == "refuse")
    {
        for (const double v : { qnan, posInf, -posInf })
            ok (probe (k, v, sr, "").status == FC_ERR_NON_FINITE, at (row, sr, "a non-finite value is FC_ERR_NON_FINITE", v));
    }
    else if (row.nonFinite == "verdict")
    {
        for (const double v : { qnan, posInf, -posInf })
        {
            const Result res = probe (k, v, sr, "");
            ok (res.status == FC_OK && res.verdict == kInvalidRequest,
                at (row, sr, "a non-finite value is InvalidRequest", v));
        }
    }
    else if (row.nonFinite == "nan-verdict")
    {
        const Result bad = probe (k, qnan, sr, "");
        ok (bad.status == FC_OK && bad.verdict == kInvalidRequest, at (row, sr, "a NaN is InvalidRequest", qnan));
        for (const double v : { posInf, -posInf })
            ok (admitted (row, probe (k, v, sr, "")), at (row, sr, "an infinity is an ordinary value here", v));
    }
    else if (row.nonFinite == "off")
    {
        for (const double v : { qnan, posInf, -posInf })
            ok (admitted (row, probe (k, v, sr, "")), at (row, sr, "a non-finite value switches the field off", v));
    }
    else
    {
        ok (row.nonFinite == "none", row.field + ": unknown nonFinite `" + row.nonFinite + "`");
    }
}

//==============================================================================
// WHAT THE TABLE SAYS IN WORDS — the `depends` column, run
//
// Every row above is probed against ONE base geometry, so a bound that moves with another FIELD is pinned only
// in that geometry. These are the checks for the rest of it.

const Row* rowFor (const std::vector<Row>& rows, const char* field)
{
    for (const Row& r : rows) if (r.field == field) return &r;
    return nullptr;
}

fc_status createWith (const fc_master_config& c)
{
    fc_master h = 0;
    const fc_status st = fc_master_create (&c, &h);
    if (st == FC_OK) fc_master_destroy (h);
    return st;
}

// `fc_master_configure` at this geometry, with `resolved.limiterLookahead` (in SAMPLES) read back.
bool lookaheadSamples (double sr, double ms, int& out)
{
    World w = baseWorld (sr);
    w.cfg.limiterLookaheadMs = ms;
    fc_master h = 0;
    if (fc_master_create (&w.cfg, &h) != FC_OK) return false;
    fc_master_resolved r {}; FC_INIT (r);
    const bool good = fc_master_configure (h, &w.prm, &r) == FC_OK;
    out = r.limiterLookahead;
    fc_master_destroy (h);
    return good;
}

std::vector<double> renderOf (const World& in)
{
    World w = in;
    fc_master h = 0;
    if (fc_master_create (&w.cfg, &h) != FC_OK) return {};
    fc_master_resolved r {}; FC_INIT (r);
    std::vector<double> got;
    if (fc_master_configure (h, &w.prm, &r) == FC_OK)
    {
        std::vector<float> dst ((std::size_t) kSolveFrames * 2, 0.0f);
        if (fc_master_process (h, programme().data(), dst.data(), kSolveFrames) == FC_OK)
            got.assign (dst.begin(), dst.end());
    }
    fc_master_destroy (h);
    return got;
}

std::vector<double> curveWithSlope (double sr, double slope)
{
    World w = baseWorld (sr);
    w.prm.eqBands[0].type = 3;                          // FC_FILTER_HIGH_PASS — the slope is read by the cuts
    w.prm.eqBands[0].lanes[0].slope = toI32 (slope);
    return eqCurve (w.prm, sr);
}

void checkDependencies (const std::vector<Row>& rows)
{
    group ("the `depends` column, run");

    // channels × monoBass — the stage folds a stereo pair and nothing else.
    {
        World w = baseWorld (kFsA);
        w.cfg.monoBass = 0; w.cfg.channels = 16;
        ok (createWith (w.cfg) == FC_OK, "16 channels are admitted with the mono-bass stage off");
        w.cfg.channels = 17;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "17 are refused whatever the stage does");
        w.cfg.monoBass = 1; w.cfg.channels = 16;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "and with the stage ON a width of 16 is refused");
        w.cfg.channels = 1;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "as is a width of 1");
        w.cfg.channels = 2;
        ok (createWith (w.cfg) == FC_OK, "2 is the one width the stage admits");
    }

    // oversampleFactor / tapsPerPhase — the ceilings belong to the two stages that oversample.
    {
        World w = baseWorld (kFsA);
        w.cfg.limiter = 0; w.cfg.clipper = 1;
        w.cfg.oversampleFactor = 64;
        ok (createWith (w.cfg) == FC_OK, "the clipper alone takes an oversample factor of 64");
        w.cfg.oversampleFactor = 65;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "and refuses 65");
        w.cfg.clipper = 0;
        w.cfg.oversampleFactor = 1000000; w.cfg.tapsPerPhase = 100000;
        ok (createWith (w.cfg) == FC_OK, "with neither stage on, nothing above the floors is refused at all");
        w.cfg.oversampleFactor = 1;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "the floor of 2 stands even then");
    }

    // compressorLookaheadMs — the 250 ms ceiling is the compressor's own.
    {
        World w = baseWorld (kFsA);
        w.cfg.compressor = 0; w.cfg.compressorLookaheadMs = 300.0;
        ok (createWith (w.cfg) == FC_OK, "300 ms of compressor lookahead are admitted with the stage off");
        w.cfg.compressor = 1;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "and refused with it on");
    }

    // limiterLookaheadMs — a clamp whose applied value is readable only in SAMPLES.
    {
        const Row* row = rowFor (rows, "fc_master_config.limiterLookaheadMs");
        ok (row != nullptr, "the table carries a row for limiterLookaheadMs");
        if (row != nullptr)
            for (const double sr : { kFsA, kFsB })
            {
                const std::string where = " at " + std::to_string ((int) sr) + " Hz";
                const double lo = row->min.at (sr, 1.0), hi = row->max.at (sr, 1.0);
                int floorS = 0, belowS = 0, twiceS = 0, capS = 0, beyondS = 0, insideS = 0;
                const bool got = lookaheadSamples (sr, lo, floorS) && lookaheadSamples (sr, lo * 0.5, belowS)
                              && lookaheadSamples (sr, lo * 2.0, twiceS)
                              && lookaheadSamples (sr, hi, capS) && lookaheadSamples (sr, hi * 1000.0, beyondS)
                              && lookaheadSamples (sr, hi - 1000.0 / sr, insideS);
                ok (got, "the limiter lookahead reads back at every probe" + where);
                ok (floorS == 2, "the table's floor is exactly 2 baseband samples" + where);
                ok (belowS == floorS, "below it the lookahead is clamped, not refused" + where);
                // The floor is the SMALLEST time that rounds to 2 samples, so twice it must be 4 — a floor
                // stated too low still reads 2 samples, because the clamp has already absorbed it.
                ok (twiceS == 4, "and twice the floor is 4 samples, which a floor stated too low is not" + where);
                ok (capS == beyondS, "the table's ceiling is where the lookahead stops growing" + where);
                ok (insideS == capS - 1, "and one sample inside it the lookahead is exactly one shorter" + where);
            }
        World w = baseWorld (kFsA);
        w.cfg.limiterLookaheadMs = -1.0;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "a NEGATIVE lookahead is refused rather than clamped");
    }

    // deliveryRate — the interval is necessary and not sufficient.
    {
        World w = baseWorld (kFsA);
        w.cfg.deliveryRate = 0.0;
        ok (createWith (w.cfg) == FC_OK, "a delivery rate of 0 is `no conversion`");
        w.cfg.deliveryRate = 44100.0;
        ok (createWith (w.cfg) == FC_OK, "48000 -> 44100 has a route");
        w.cfg.deliveryRate = 44100.5;
        ok (createWith (w.cfg) == FC_ERR_REFUSED_BY_CORE, "and half a hertz has none: the rate must be whole");
    }

    // eqBands[].lanes[].slope — read as slope/6 POLES; the clamp shows in the curve and nowhere else.
    {
        const Row* row = rowFor (rows, "fc_master_params.eqBands[].lanes[].slope");
        ok (row != nullptr, "the table carries a row for the lane slope");
        if (row != nullptr)
        {
            // THE STEP IS SIX, the divisor the design reads the slope through: one pole. A probe closer than
            // that shares a pole with the bound and cannot tell a floor of 6 from one of -10000.
            constexpr double kPole = 6.0;
            const double lo = row->min.at (kFsA, 1.0), hi = row->max.at (kFsA, 1.0);
            const std::vector<double> atLo = curveWithSlope (kFsA, lo);
            const std::vector<double> atHi = curveWithSlope (kFsA, hi);
            ok (! atLo.empty() && ! atHi.empty(), "the slope curve is answered");
            ok (sameCurve (atLo, curveWithSlope (kFsA, lo - kPole)), "one pole below the floor is the floor");
            ok (sameCurve (atLo, curveWithSlope (kFsA, -1.0e6)), "a hugely negative slope included");
            ok (! sameCurve (atLo, curveWithSlope (kFsA, lo + kPole)),
                "and one pole ABOVE it is a different filter, which a floor stated too low is not");
            ok (sameCurve (atHi, curveWithSlope (kFsA, hi + kPole)), "one pole above the ceiling is the ceiling");
            ok (sameCurve (atHi, curveWithSlope (kFsA, 1.0e9)), "a hugely positive slope included");
            ok (! sameCurve (atHi, curveWithSlope (kFsA, hi - kPole)), "and one pole below it is a different filter");
            ok (! sameCurve (atLo, atHi), "the two ends are different filters");
        }
    }

    // THE `dyn` GROUP IS LIVE, which is what makes the four rows above pinnable by a render at all. The
    // chain's EQ stage drives `dynamiceq::LaneDynamics` into each point's delta seam, so an armed point
    // is a different render — and `dyn.on` is the gate: with it off, every other field in the group is
    // carried and changes nothing.
    {
        World w = baseWorld (kFsA);
        w.cfg.dither = 0;
        w.prm.eqBands[0].dyn.on = 1; w.prm.eqBands[0].dyn.thrAuto = 0;
        w.prm.eqBands[0].dyn.thrDb = -30.0; w.prm.eqBands[0].dyn.rangeDb = -30.0;
        const std::vector<double> armed = renderOf (w);
        World off = w; off.prm.eqBands[0].dyn.on = 0;
        const std::vector<double> unarmed = renderOf (off);
        ok (! armed.empty() && ! unarmed.empty(), "the dynamics probe renders");
        ok (! sameCurve (armed, unarmed),
            "the band dynamics MOVE the render: an armed point is not the same samples as dyn.on = 0");

        // ...AND THE GROUP IS GATED BY THAT ONE FLAG. Every other field moved to its far end, with
        // `dyn.on` off, is the same render — the claim `dyn.on`'s own row makes by being a flag.
        World other = off;
        other.prm.eqBands[0].dyn.rangeDb = 30.0;
        other.prm.eqBands[0].dyn.thrDb   = 24.0;
        other.prm.eqBands[0].dyn.thrAuto = 1;
        other.prm.eqBands[0].dyn.atk     = 0.0;
        other.prm.eqBands[0].dyn.rel     = 1.0;
        ok (sameCurve (unarmed, renderOf (other)),
            "with dyn.on = 0 the rest of the group changes no sample");

        // RANGE 0 IS NO DYNAMICS, whatever `dyn.on` says — the sentence its own row carries.
        World zero = w; zero.prm.eqBands[0].dyn.rangeDb = 0.0;
        ok (sameCurve (unarmed, renderOf (zero)), "an armed point with rangeDb = 0 renders as an unarmed one");
    }

    // dither.bits — a depth of 32 or more is a BYPASS, accepted like every other value.
    {
        World w = baseWorld (kFsA);
        fc_master h = 0;
        ok (fc_master_create (&w.cfg, &h) == FC_OK, "a handle for the dither depths");
        if (h != 0)
        {
            fc_master_resolved r {}; FC_INIT (r);
            w.prm.dither.bits = 32;
            ok (fc_master_configure (h, &w.prm, &r) == FC_OK, "a dither depth of 32 is accepted, not refused");
            fc_master_destroy (h);
        }
    }
}

}   // namespace

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf (stderr, "usage: %s <path to tools/wasm/fc-master-layout.mjs>\n", argv[0]);
        return 2;
    }

    std::vector<Row> rows;
    std::string why;
    if (! readDomains (argv[1], rows, why))
    {
        std::fprintf (stderr, "FC_DOMAINS: %s\n", why.c_str());
        return 2;
    }
    std::printf ("felitronics_master_domains_tests — %zu rows from %s\n", rows.size(), argv[1]);

    // ── both directions, before anything is probed ────────────────────────────────────────────────
    group ("every row is driven and every knob is described");
    for (const Row& r : rows)
        ok (knobFor (r.field) != nullptr, r.field + ": a row with nothing here that can write it");
    for (const Knob& k : kKnobs)
        ok (rowFor (rows, k.field) != nullptr, std::string (k.field) + ": written here and absent from FC_DOMAINS");

    group ("the domains, at 48 kHz");
    for (const Row& r : rows)
        if (const Knob* k = knobFor (r.field)) checkRow (r, *k, kFsA);

    // A bound written in `sr` is a formula, and a formula checked at one rate is a constant.
    group ("the rate-dependent domains again, at 96 kHz");
    int scaled = 0;
    for (const Row& r : rows)
        if (r.min.scaled() || r.max.scaled())
            if (const Knob* k = knobFor (r.field)) { ++scaled; checkRow (r, *k, kFsB); }
    ok (scaled >= 5, "the table still carries rate-dependent bounds to re-check (" + std::to_string (scaled) + ")");

    checkDependencies (rows);

    // THE HOLE, NAMED. A clamp with no read-back can be shown to be a clamp rather than a refusal, and no more:
    // where the bound SITS rests on the reading of the code. Two are pinned by a check of their own above and
    // are excluded; the rest are printed on every run.
    static const char* kPinnedByName[] = { "fc_master_config.limiterLookaheadMs",
                                           "fc_master_params.eqBands[].lanes[].slope" };
    for (const char* f : kPinnedByName)
    {
        const Row* r = rowFor (rows, f);
        ok (r != nullptr && r->edge == "clamp" && r->resolved.empty(),
            std::string (f) + ": pinned by name above, so it must still be the clamp with no read-back that needs it");
    }
    std::printf ("\n  clamps with no read-back and no check of their own — accepted, bound NOT pinned:\n");
    int blind = 0;
    for (const Row& r : rows)
    {
        bool named = false;
        for (const char* f : kPinnedByName) named = named || r.field == f;
        if (r.edge == "clamp" && r.resolved.empty() && ! named) { ++blind; std::printf ("      %s\n", r.field.c_str()); }
    }
    std::printf ("    %d of %zu rows\n", blind, rows.size());

    return felitronics::test::report();
}
