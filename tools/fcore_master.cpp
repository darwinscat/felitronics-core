// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_master — the native CLI that drives the C ABI (tools/fc_master_abi.h, implemented in
// tools/wasm/fc_master.cpp), and the reference that proves the ABI adds nothing.
//
//   fcore_master render   <sampleRate> <channels> <in.f32le> <out.f32le> [key=value ...]
//   fcore_master solve    <sampleRate> <channels> <in.f32le> <out.f32le> target=<LUFS> tp=<dBTP> [key=value ...]
//   fcore_master lra      <sampleRate> <channels> <in.f32le>
//   fcore_master selftest [sampleRate] [channels]
//
// I/O is interleaved 32-bit-float little-endian PCM, exactly what `ffmpeg -f f32le` emits and reads and
// what tools/fcore_measure.cpp already speaks, so a harness can put this binary in a pipeline next to
// ffmpeg without a converter in between.
//
// ====================================================================================
// WHY `selftest` IS THE ACCEPTANCE AND NOT A SMOKE TEST
// ====================================================================================
// It renders one programme twice — once through the C ABI, once by calling
// `mastering::OfflineRenderer` directly — and compares the two BIT FOR BIT, in this binary, on this
// machine. That scope is the point: one libm, one set of FP flags, one instruction set. Anything that
// differs is a marshalling, ownership or ordering fault in the facade, because nothing else is left.
// (Whether two TIERS agree is a different question with a different answer and its own task; it is not
// this one, and folding it in here would weaken the only check that can be exact.)
//
// THE TWO PATHS MUST DIFFER IN NOTHING BUT THE ABI, which turns out to be four separate disciplines,
// every one of them learned from a measurement rather than reasoned out:
//
//   * SAME FP FLAGS, and this target does NOT inherit `tools/`' `-ffp-contract=off`. Measured on one
//     source: `on` against `off` differs in 203 269 of 288 000 samples. Worse, two TUs of one binary
//     built with different flags measured as AGREEING — because the linker merged the header-only
//     instantiations and whichever TU came first won. Agreement by link order is exactly the
//     false-green this file exists to prevent, so the flag is stated once, for the whole target.
//   * SAME LIFECYCLE ORDER. `setParams` then `prepare` is not the same render as `prepare` then
//     `setParams`: measured, 59 259 of 80 000 samples and 0.0715 full scale apart, all of it
//     `stereo::MonoBass` (its `reset()` SNAPS the width, its setter RAMPS it over 20 ms). The ABI's
//     `configure` re-prepares, i.e. takes the first order, so the direct path here does the same.
//   * SAME DEFAULTS. A zeroed parameter struct is not `MasteringChainParams{}` — 287 998 of 288 000
//     samples apart — so both paths start from `fc_master_params_default()`.
//   * A RESET BEFORE EVERY RENDER. `OfflineRenderer` resets the chain itself; the ABI has no renderer,
//     so a second programme through the same handle without `fc_master_reset` differs from the first by
//     151 884 samples. The loop below resets, and the harness must too.

#include "fc_master_abi.h"

#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;
using namespace felitronics::mastering;

namespace
{

const char* statusName (fc_status s)
{
    switch (s)
    {
        case FC_OK:                  return "OK";
        case FC_ERR_HANDLE:          return "HANDLE";
        case FC_ERR_ABI_VERSION:     return "ABI_VERSION";
        case FC_ERR_STRUCT_SIZE:     return "STRUCT_SIZE";
        case FC_ERR_NULL:            return "NULL";
        case FC_ERR_ALIGNMENT:       return "ALIGNMENT";
        case FC_ERR_SPAN:            return "SPAN";
        case FC_ERR_ENUM:            return "ENUM";
        case FC_ERR_RANGE:           return "RANGE";
        case FC_ERR_CAPACITY:        return "CAPACITY";
        case FC_ERR_STATE:           return "STATE";
        case FC_ERR_NON_FINITE:      return "NON_FINITE";
        case FC_ERR_REFUSED_BY_CORE: return "REFUSED_BY_CORE";
        case FC_ERR_EXHAUSTED:       return "EXHAUSTED";
    }
    return "?";
}

//==============================================================================
// KEY=VALUE PARAMETERS
//
// Deliberately flat and deliberately STRICT: an unknown key is an error rather than a shrug, because a
// harness that mistypes `lim.ceiling` and gets the default silently produces a whole table of numbers
// describing settings nobody chose. Same reason the ABI refuses an enum code it does not know.
struct Args
{
    fc_master_config cfg {};
    fc_master_params prm {};
    fc_loudness_request req {};
    std::uint32_t block = 4096;
};

// STRICT number parsing. `atof`/`atoi` stop at the first character they do not understand and report
// nothing, so `inputGainDb=oops` becomes 0 dB and `comp.ratio=4oops` becomes 4 — a harness row computed
// from settings nobody chose, which is the whole failure mode this CLI's strictness exists to prevent.
bool num (const std::string& v, double& out)
{
    if (v.empty()) return false;
    char* end = nullptr;
    const double d = std::strtod (v.c_str(), &end);
    if (end == v.c_str() || *end != '\0') return false;
    out = d;
    return true;
}

bool inum (const std::string& v, long& out)
{
    if (v.empty()) return false;
    char* end = nullptr;
    const long i = std::strtol (v.c_str(), &end, 10);
    if (end == v.c_str() || *end != '\0') return false;
    out = i;
    return true;
}

bool parseBool (const std::string& v, std::int32_t& out)
{
    if (v == "1" || v == "on" || v == "true")  { out = 1; return true; }
    if (v == "0" || v == "off" || v == "false") { out = 0; return true; }
    return false;
}

bool parseEnumName (const std::string& v, const char* const* names, int count, std::int32_t& out)
{
    for (int i = 0; i < count; ++i) if (v == names[i]) { out = i; return true; }
    return false;
}

// `bandN.field` — lane 0 (Stereo) only. The five other lanes are reachable through the ABI and are not
// reachable through this CLI, which is a limit of the CLI and is said here rather than discovered: the
// harness's presets drive a high-pass and a handful of tone controls, all of them Stereo.
bool applyBandKey (Args& a, const std::string& key, const std::string& val)
{
    if (key.rfind ("band", 0) != 0) return false;
    const std::size_t dot = key.find ('.');
    if (dot == std::string::npos) return false;
    // The INDEX must be the whole of what sits between "band" and the dot. `atoi` would read `bandXYZ`
    // as band 0 and quietly apply the value to the wrong band.
    long idxL = 0;
    if (! inum (key.substr (4, dot - 4), idxL)) return false;
    const int idx = (int) idxL;
    if (idx < 0 || idx >= FC_MAX_EQ_BANDS) return false;
    const std::string f = key.substr (dot + 1);
    fc_eq_band& b = a.prm.eqBands[idx];
    static const char* kTypes[] { "bell", "lowshelf", "highshelf", "highpass", "lowpass",
                                  "bandpass", "notch", "allpass", "tilt" };
    if (f == "on")     return parseBool (val, b.on);
    if (f == "bypass") return parseBool (val, b.bypass);
    if (f == "swept")  return parseBool (val, b.swept);
    if (f == "type")   return parseEnumName (val, kTypes, 9, b.type);
    if (f == "laneon") return parseBool (val, b.lanes[0].on);
    double d = 0.0; long i = 0;
    if (f == "freq")  { if (! num (val, d)) return false; b.lanes[0].freq   = d; return true; }
    if (f == "q")     { if (! num (val, d)) return false; b.lanes[0].q      = d; return true; }
    if (f == "gain")  { if (! num (val, d)) return false; b.lanes[0].gainDb = d; return true; }
    if (f == "slope") { if (! inum (val, i)) return false; b.lanes[0].slope = (std::int32_t) i; return true; }
    return false;
}

bool applyKey (Args& a, const std::string& key, const std::string& val)
{
    if (applyBandKey (a, key, val)) return true;

    // Every numeric key below reads through these two, so a value that is not a whole number is a
    // refusal rather than a silent zero. `dNeeded`/`iNeeded` are false when the value did not parse;
    // each use tests them, so one malformed value cannot slip past as "the default".
    double d = 0.0; long i = 0;
    const bool dOk = num (val, d);
    const bool iOk = inum (val, i);
    (void) dOk;
    #define FC_D(expr) do { if (! dOk) return false; expr; return true; } while (0)
    #define FC_I(expr) do { if (! iOk) return false; expr; return true; } while (0)

    if (key == "block")          { FC_I (a.block = (std::uint32_t) (i > 0 ? i : 1)); }
    if (key == "internalBlock")  { FC_I (a.cfg.internalBlock = i); }
    if (key == "oversample")     { FC_I (a.cfg.oversampleFactor = i); }
    if (key == "taps")           { FC_I (a.cfg.tapsPerPhase = i); }
    if (key == "compLookaheadMs"){ FC_D (a.cfg.compressorLookaheadMs = d); }
    if (key == "limLookaheadMs") { FC_D (a.cfg.limiterLookaheadMs = d); }
    if (key == "sidechainHpfHz") { FC_D (a.cfg.sidechainHpfHz = d); }
    if (key == "eq")             return parseBool (val, a.cfg.eq);
    if (key == "monoBass")       return parseBool (val, a.cfg.monoBass);
    if (key == "compressor")     return parseBool (val, a.cfg.compressor);
    if (key == "clipper")        return parseBool (val, a.cfg.clipper);
    if (key == "limiter")        return parseBool (val, a.cfg.limiter);
    if (key == "dither")         return parseBool (val, a.cfg.dither);

    if (key == "inputGainDb")      { FC_D (a.prm.inputGainDb = d); }
    if (key == "preLimiterGainDb") { FC_D (a.prm.preLimiterGainDb = d); }

    if (key == "mb.on")    return parseBool (val, a.prm.monoBass.enabled);
    if (key == "mb.freq")  { FC_D (a.prm.monoBass.frequencyHz = (float) d); }
    if (key == "mb.width") { FC_D (a.prm.monoBass.lowWidth = (float) d); }

    static const char* kDet[]   { "peak", "rms" };
    static const char* kLink[]  { "max", "meanpower" };
    static const char* kMode[]  { "downcompress", "upcompress", "downexpand" };
    static const char* kShape[] { "tanh", "atan", "cubic", "asym" };
    static const char* kShap[]  { "none", "weighted", "psycho" };
    static const char* kGrSt[]  { "mean", "p95", "max" };

    if (key == "comp.detector")  return parseEnumName (val, kDet,  2, a.prm.compressor.detector);
    if (key == "comp.link")      return parseEnumName (val, kLink, 2, a.prm.compressor.link);
    if (key == "comp.mode")      return parseEnumName (val, kMode, 3, a.prm.compressor.mode);
    if (key == "comp.rmsMs")     { FC_D (a.prm.compressor.rmsWindowMs = d); }
    if (key == "comp.threshold") { FC_D (a.prm.compressor.thresholdDb = d); }
    if (key == "comp.ratio")     { FC_D (a.prm.compressor.ratio = d); }
    if (key == "comp.knee")      { FC_D (a.prm.compressor.kneeDb = d); }
    if (key == "comp.range")     { FC_D (a.prm.compressor.rangeDb = d); }
    if (key == "comp.attack")    { FC_D (a.prm.compressor.attackMs = d); }
    if (key == "comp.release")   { FC_D (a.prm.compressor.releaseMs = d); }
    if (key == "comp.makeup")    { FC_D (a.prm.compressor.makeupDb = d); }
    if (key == "comp.autoMakeup")return parseBool (val, a.prm.compressor.autoMakeup);

    if (key == "clip.shape")  return parseEnumName (val, kShape, 4, a.prm.clipper.shape);
    if (key == "clip.drive")  { FC_D (a.prm.clipper.driveDb  = (float) d); }
    if (key == "clip.bias")   { FC_D (a.prm.clipper.bias     = (float) d); }
    if (key == "clip.mix")    { FC_D (a.prm.clipper.mix      = (float) d); }
    if (key == "clip.output") { FC_D (a.prm.clipper.outputDb = (float) d); }
    if (key == "clip.autoComp") { FC_D (a.prm.clipper.autoComp = (float) d); }
    if (key == "clip.dcBlockHz"){ FC_D (a.prm.clipper.dcBlockHz = (float) d); }

    if (key == "lim.ceiling") { FC_D (a.prm.limiter.ceilingDbTp = d); }
    if (key == "lim.release") { FC_D (a.prm.limiter.releaseMs = d); }

    if (key == "dith.bits")    { FC_I (a.prm.dither.bits = i); }
    if (key == "dith.shaping") return parseEnumName (val, kShap, 3, a.prm.dither.shaping);
    if (key == "dith.seedLo")  { char* e = nullptr; const unsigned long u = std::strtoul (val.c_str(), &e, 0);
                                 if (e == val.c_str() || *e != '\0') return false;
                                 a.prm.dither.seedLo = (std::uint32_t) u; return true; }
    if (key == "dith.seedHi")  { char* e = nullptr; const unsigned long u = std::strtoul (val.c_str(), &e, 0);
                                 if (e == val.c_str() || *e != '\0') return false;
                                 a.prm.dither.seedHi = (std::uint32_t) u; return true; }
    if (key == "dith.autoBlank") return parseBool (val, a.prm.dither.autoBlank);
    if (key == "dith.autoBlankSamples") { FC_I (a.prm.dither.autoBlankSamples = i); }

    if (key == "bypass.eq")         return parseBool (val, a.prm.bypassEq);
    if (key == "bypass.monoBass")   return parseBool (val, a.prm.bypassMonoBass);
    if (key == "bypass.compressor") return parseBool (val, a.prm.bypassCompressor);
    if (key == "bypass.clipper")    return parseBool (val, a.prm.bypassClipper);
    if (key == "bypass.limiter")    return parseBool (val, a.prm.bypassLimiter);
    if (key == "bypass.dither")     return parseBool (val, a.prm.bypassDither);

    if (key == "target")     { FC_D (a.req.targetLufs = d); }
    if (key == "tp")         { FC_D (a.req.maxTruePeakDbTp = d); }
    if (key == "tolerance")  { FC_D (a.req.toleranceLu = d); }
    if (key == "tpAim")      { FC_D (a.req.truePeakAimDb = d); }
    if (key == "maxPasses")  { FC_I (a.req.maxPasses = i); }
    if (key == "initialGain"){ FC_D (a.req.initialGainDb = d); }
    if (key == "minPlr")     { FC_D (a.req.minPlrDb = d); }
    if (key == "maxLraLoss") { FC_D (a.req.maxLraLossLu = d); }
    if (key == "inputLra")   { FC_D (a.req.inputLoudnessRangeLu = d); }
    if (key == "limGrLimit") { FC_D (a.req.limiterGr.limitDb = d); }
    if (key == "limGrStat")  return parseEnumName (val, kGrSt, 3, a.req.limiterGr.statistic);
    if (key == "compGrLimit"){ FC_D (a.req.compressorGr.limitDb = d); }
    if (key == "compGrStat") return parseEnumName (val, kGrSt, 3, a.req.compressorGr.statistic);
    if (key == "activityDb") { FC_D (a.req.activityThresholdDb = d); }
    #undef FC_D
    #undef FC_I
    return false;
}

bool parseArgs (Args& a, int argc, char** argv, int from)
{
    fc_master_config_default (&a.cfg);
    fc_master_params_default (&a.prm);
    fc_loudness_request_default (&a.req);
    for (int i = from; i < argc; ++i)
    {
        const std::string s = argv[i];
        const std::size_t eq = s.find ('=');
        if (eq == std::string::npos) { std::fprintf (stderr, "not a key=value: %s\n", s.c_str()); return false; }
        if (! applyKey (a, s.substr (0, eq), s.substr (eq + 1)))
        { std::fprintf (stderr, "unknown or malformed key: %s\n", s.c_str()); return false; }
    }
    return true;
}

//==============================================================================
// FILE I/O — interleaved f32le in, planar in memory, interleaved f32le out.

bool readInterleaved (const char* path, int nc, std::vector<float>& planar, std::size_t& frames)
{
    std::FILE* f = std::fopen (path, "rb");
    if (f == nullptr) { std::perror ("open"); return false; }
    std::vector<float> inter;
    float buf[8192];
    std::size_t got;
    while ((got = std::fread (buf, sizeof (float), 8192, f)) > 0) inter.insert (inter.end(), buf, buf + got);
    // A READ ERROR IS NOT A SHORT FILE. Without this a truncated read becomes a successful render of a
    // PREFIX — the harness gets a row, the row gets a number, and the number is of a programme nobody
    // supplied. The same reason the ABI refuses rather than processing what it can.
    if (std::ferror (f) != 0) { std::fclose (f); std::fprintf (stderr, "read error\n"); return false; }
    std::fclose (f);
    frames = inter.size() / (std::size_t) nc;
    planar.assign (frames * (std::size_t) nc, 0.0f);
    for (std::size_t i = 0; i < frames; ++i)
        for (int c = 0; c < nc; ++c)
            planar[(std::size_t) c * frames + i] = inter[i * (std::size_t) nc + (std::size_t) c];
    return true;
}

bool writeInterleaved (const char* path, int nc, const std::vector<float>& planar, std::size_t frames)
{
    std::FILE* f = std::fopen (path, "wb");
    if (f == nullptr) { std::perror ("open"); return false; }
    std::vector<float> inter (frames * (std::size_t) nc);
    for (std::size_t i = 0; i < frames; ++i)
        for (int c = 0; c < nc; ++c)
            inter[i * (std::size_t) nc + (std::size_t) c] = planar[(std::size_t) c * frames + i];
    const bool ok = std::fwrite (inter.data(), sizeof (float), inter.size(), f) == inter.size();
    std::fclose (f);
    return ok;
}

//==============================================================================
// THE ABI RENDER — the block loop a JS worker will run, written once, here.
//
// `out[n] = y[n + D]`: the chain is a streaming object, so the first D output samples are its own
// priming and the last D leave it only after the input has ended. This is `OfflineRenderer`'s formula
// and it is re-stated rather than re-derived: the CLI cannot call that class through the C ABI, so the
// arithmetic exists in both places and the selftest is what keeps them equal.
// `inPlace = false` drives the OUT-OF-PLACE path — `fc_master_process(h, in, out, n)` with two distinct
// buffers, which is what a page does (an input view and an output view) and which the in-place calls
// never exercise: a crew round deleted the facade's `memcpy` outright and every check stayed green,
// because nothing in this file or in the suite ever passed two different pointers with content behind
// them. The copy is transport, but transport that is never run is transport that is never tested.
bool abiRender (const Args& a, const std::vector<float>& in, std::size_t frames, int nc,
                std::vector<float>& out, fc_master_resolved& res, bool inPlace = true)
{
    fc_master h = 0;
    if (const fc_status st = fc_master_create (&a.cfg, &h); st != FC_OK)
    { std::fprintf (stderr, "create: %s\n", statusName (st)); return false; }

    FC_INIT (res);
    if (const fc_status st = fc_master_configure (h, &a.prm, &res); st != FC_OK)
    { std::fprintf (stderr, "configure: %s\n", statusName (st)); fc_master_destroy (h); return false; }

    const std::size_t D = (std::size_t) res.latencySamples;
    std::vector<float> stream ((frames + D) * (std::size_t) nc, 0.0f);
    const std::size_t stride = frames + D;

    // Feed the programme, then D zeros, in `block`-frame slices. Planar with stride `stride`, so a slice
    // is not contiguous: it goes through a scratch buffer whose stride is the slice length, which is the
    // ABI's own layout rule applied honestly rather than by aliasing into the middle of a plane.
    std::vector<float> slice ((std::size_t) a.block * (std::size_t) nc);
    std::vector<float> dest  ((std::size_t) a.block * (std::size_t) nc, 0.0f);
    std::size_t written = 0;
    for (std::size_t off = 0; off < frames; )
    {
        const std::size_t m = std::min<std::size_t> (a.block, frames - off);
        for (int c = 0; c < nc; ++c)
            std::memcpy (slice.data() + (std::size_t) c * m, in.data() + (std::size_t) c * frames + off,
                         m * sizeof (float));
        float* dst = inPlace ? slice.data() : dest.data();
        if (const fc_status st = fc_master_process (h, slice.data(), dst, (std::uint32_t) m); st != FC_OK)
        { std::fprintf (stderr, "process: %s\n", statusName (st)); fc_master_destroy (h); return false; }
        for (int c = 0; c < nc; ++c)
            std::memcpy (stream.data() + (std::size_t) c * stride + written, dst + (std::size_t) c * m,
                         m * sizeof (float));
        written += m; off += m;
    }

    if (D > 0)
    {
        std::vector<float> tail (D * (std::size_t) nc, 0.0f);
        std::uint32_t got = 0;
        if (const fc_status st = fc_master_flush (h, tail.data(), (std::uint32_t) D, &got); st != FC_OK)
        { std::fprintf (stderr, "flush: %s\n", statusName (st)); fc_master_destroy (h); return false; }
        for (int c = 0; c < nc; ++c)
            std::memcpy (stream.data() + (std::size_t) c * stride + written,
                         tail.data() + (std::size_t) c * D, (std::size_t) got * sizeof (float));
        written += got;
    }

    out.assign (frames * (std::size_t) nc, 0.0f);
    for (int c = 0; c < nc; ++c)
        for (std::size_t n = 0; n < frames; ++n)
            out[(std::size_t) c * frames + n] = stream[(std::size_t) c * stride + n + D];

    fc_master_destroy (h);
    return true;
}

// The same render through the C++ API, with the SAME lifecycle order the ABI takes (parameters written
// BEFORE prepare, which is what `fc_master_configure` does when it re-prepares).
bool directRender (const Args& a, const std::vector<float>& in, std::size_t frames, int nc,
                   std::vector<float>& out)
{
    MasteringChainConfig cc {};
    cc.internalBlock         = a.cfg.internalBlock;
    cc.eq                    = a.cfg.eq != 0;
    cc.monoBass              = a.cfg.monoBass != 0;
    cc.compressor            = a.cfg.compressor != 0;
    cc.clipper               = a.cfg.clipper != 0;
    cc.limiter               = a.cfg.limiter != 0;
    cc.dither                = a.cfg.dither != 0;
    cc.compressorLookaheadMs = a.cfg.compressorLookaheadMs;
    cc.limiterLookaheadMs    = a.cfg.limiterLookaheadMs;
    cc.oversampleFactor      = a.cfg.oversampleFactor;
    cc.tapsPerPhase          = a.cfg.tapsPerPhase;
    cc.sidechainHpfHz        = a.cfg.sidechainHpfHz;

    // The parameter set comes from the SAME mapper the ABI uses — through the ABI's own defaults and a
    // throwaway handle would be circular, so it is built here by hand from the same `Args`. This is the
    // one place the CLI mirrors the mapping, and the selftest is what proves the mirror is faithful.
    MasteringChainParams cp {};
    cp.inputGainDb      = a.prm.inputGainDb;
    cp.preLimiterGainDb = a.prm.preLimiterGainDb;
    for (int b = 0; b < FC_MAX_EQ_BANDS; ++b)
    {
        const fc_eq_band& sb = a.prm.eqBands[b];
        eq::BandParams& db = cp.eqBands[b];
        db.on = sb.on != 0; db.type = (eq::FilterType) sb.type;
        db.swept = sb.swept != 0; db.bypass = sb.bypass != 0;
        db.dyn.on = sb.dyn.on != 0; db.dyn.rangeDb = sb.dyn.rangeDb; db.dyn.thrDb = sb.dyn.thrDb;
        db.dyn.thrAuto = sb.dyn.thrAuto != 0; db.dyn.atk = sb.dyn.atk; db.dyn.rel = sb.dyn.rel;
        for (int l = 0; l < FC_MAX_EQ_LANES; ++l)
        {
            const fc_eq_lane& sl = sb.lanes[l];
            eq::LaneParams& dl = db.lanes[l];
            dl.on = sl.on != 0; dl.freq = sl.freq; dl.Q = sl.q; dl.gainDb = sl.gainDb;
            dl.slope = sl.slope; dl.bypass = sl.bypass != 0;
        }
    }
    cp.monoBass.enabled = a.prm.monoBass.enabled != 0;
    cp.monoBass.frequencyHz = a.prm.monoBass.frequencyHz;
    cp.monoBass.lowWidth = a.prm.monoBass.lowWidth;
    cp.compressor.detector = (dynamics::Detector) a.prm.compressor.detector;
    cp.compressor.link = (dynamics::LinkMode) a.prm.compressor.link;
    cp.compressor.mode = (dynamics::Mode) a.prm.compressor.mode;
    cp.compressor.rmsWindowMs = a.prm.compressor.rmsWindowMs;
    cp.compressor.thresholdDb = a.prm.compressor.thresholdDb;
    cp.compressor.ratio = a.prm.compressor.ratio;
    cp.compressor.kneeDb = a.prm.compressor.kneeDb;
    cp.compressor.rangeDb = a.prm.compressor.rangeDb;
    cp.compressor.attackMs = a.prm.compressor.attackMs;
    cp.compressor.releaseMs = a.prm.compressor.releaseMs;
    cp.compressor.makeupDb = a.prm.compressor.makeupDb;
    cp.compressor.autoMakeup = a.prm.compressor.autoMakeup != 0;
    cp.clipper.shape = (saturation::WaveShaper::Shape) a.prm.clipper.shape;
    cp.clipper.driveDb = a.prm.clipper.driveDb;
    cp.clipper.bias = a.prm.clipper.bias;
    cp.clipper.mix = a.prm.clipper.mix;
    cp.clipper.outputDb = a.prm.clipper.outputDb;
    cp.clipper.autoComp = a.prm.clipper.autoComp;
    cp.clipper.dcBlockHz = a.prm.clipper.dcBlockHz;
    cp.limiter.ceilingDbTp = a.prm.limiter.ceilingDbTp;
    cp.limiter.releaseMs = a.prm.limiter.releaseMs;
    cp.dither.bits = a.prm.dither.bits;
    cp.dither.shaping = (dither::NoiseShaping) a.prm.dither.shaping;
    cp.dither.seed = ((std::uint64_t) a.prm.dither.seedHi << 32) | (std::uint64_t) a.prm.dither.seedLo;
    cp.dither.autoBlank = a.prm.dither.autoBlank != 0;
    cp.dither.autoBlankSamples = a.prm.dither.autoBlankSamples;
    cp.bypassEq = a.prm.bypassEq != 0;
    cp.bypassMonoBass = a.prm.bypassMonoBass != 0;
    cp.bypassCompressor = a.prm.bypassCompressor != 0;
    cp.bypassClipper = a.prm.bypassClipper != 0;
    cp.bypassLimiter = a.prm.bypassLimiter != 0;
    cp.bypassDither = a.prm.bypassDither != 0;

    MasteringChain chain;
    OfflineRenderer r;
    if (! r.prepare (nc, (int) a.block)) return false;
    chain.setParams (cp);                                        // SAME ORDER as fc_master_configure
    if (! chain.prepare (a.cfg.sampleRate, nc, cc)) return false;

    out.assign (frames * (std::size_t) nc, 0.0f);
    const float* ip[core::kMaxChannels] {};
    float*       op[core::kMaxChannels] {};
    for (int c = 0; c < nc; ++c)
    {
        ip[c] = in.data() + (std::size_t) c * frames;
        op[c] = out.data() + (std::size_t) c * frames;
    }
    return r.render (chain, ip, op, nc, (int) frames);
}

//==============================================================================
// SELFTEST

std::size_t bitDiff (const std::vector<float>& a, const std::vector<float>& b, double& worst)
{
    worst = 0.0;
    if (a.size() != b.size()) return a.size() + b.size();
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (std::memcmp (&a[i], &b[i], sizeof (float)) != 0)
        { ++n; worst = std::max (worst, (double) std::fabs (a[i] - b[i])); }
    return n;
}

// A programme with something for every stage to do: two tones, a slow envelope so the compressor moves,
// and periodic transients so the limiter and the clipper are not asleep. A fixture on which every stage
// is inert would make this whole comparison pass for the wrong reason.
std::vector<float> programme (double fs, int nc, std::size_t frames)
{
    std::vector<float> v (frames * (std::size_t) nc, 0.0f);
    for (int c = 0; c < nc; ++c)
        for (std::size_t i = 0; i < frames; ++i)
        {
            const double t = (double) i / fs;
            const double env = 0.25 + 0.55 * (0.5 + 0.5 * std::sin (2.0 * 3.14159265358979 * 0.7 * t));
            double x = env * (0.6 * std::sin (2.0 * 3.14159265358979 * (220.0 + 55.0 * c) * t)
                            + 0.3 * std::sin (2.0 * 3.14159265358979 * 3100.0 * t));
            if (i % 9600 < 24) x += 0.55;                      // transients, 5 per second
            v[(std::size_t) c * frames + i] = (float) x;
        }
    return v;
}

int selftest (double fs, int nc)
{
    int failures = 0;
    auto check = [&failures] (bool ok, const char* what, const char* detail = "")
    {
        std::printf ("  [%s] %s%s%s\n", ok ? "ok" : "FAIL", what, *detail ? " — " : "", detail);
        if (! ok) ++failures;
    };

    const std::size_t frames = (std::size_t) (fs * 4.0);
    const auto in = programme (fs, nc, frames);

    Args a;
    fc_master_config_default (&a.cfg);
    fc_master_params_default (&a.prm);
    fc_loudness_request_default (&a.req);
    a.cfg.sampleRate = fs;
    a.cfg.channels   = nc;
    a.cfg.monoBass   = (nc == 2) ? 1 : 0;                 // stereo-only stage; exercised when it can be
    a.cfg.clipper    = 1;                                 // every stage present, so nothing is untested

    // EVERY field moved off its default, and none of them representable in binary32. Two separate
    // failure modes are being closed here and both were measured: a value that survives a stray
    // double->float narrowing unchanged (0.5, -1.0, 60) lets that narrowing pass, and a field left at
    // its DEFAULT lets a MISSING mapping pass — because the direct path builds its parameters from the
    // same `Args`, so a dropped field lands on the same default on both sides and the bit-compare is
    // green. A crew round's own stand dropped nine fields this way and every one survived.
    a.cfg.compressorLookaheadMs = 1.7;
    a.cfg.limiterLookaheadMs    = 1.3;
    a.cfg.tapsPerPhase          = 48;
    a.cfg.oversampleFactor      = 2;
    a.cfg.sidechainHpfHz        = 47.0;
    a.cfg.internalBlock         = 128;
    a.prm.compressor.detector   = FC_DETECTOR_RMS;     // NOT the default: a wrong-but-VALID enum
    a.prm.compressor.link       = FC_LINK_MEAN_POWER;  // translation is invisible on default values
    a.prm.compressor.mode       = FC_COMP_DOWN_COMPRESS;
    a.prm.clipper.shape         = FC_SHAPE_ATAN;
    a.prm.dither.shaping        = FC_SHAPING_PSYCHO;
    a.prm.dither.autoBlank      = 0;
    a.prm.dither.autoBlankSamples = 3777;
    a.prm.compressor.autoMakeup = 1;
    a.prm.eqBands[2].on         = 1;
    a.prm.eqBands[2].type       = FC_FILTER_HIGH_SHELF;
    a.prm.eqBands[2].swept      = 1;
    a.prm.eqBands[2].lanes[3].on = 1;                  // a lane that is NOT lane 0
    a.prm.eqBands[2].lanes[3].freq = 7331.7;
    a.prm.eqBands[2].lanes[3].gainDb = 1.7;
    a.prm.eqBands[2].lanes[3].q = 0.77;
    a.prm.eqBands[2].lanes[4].on = 1;
    a.prm.eqBands[2].lanes[4].bypass = 1;              // `bypass` is not `on`, and both are mapped
    a.prm.eqBands[2].dyn.on     = 1;
    a.prm.eqBands[2].dyn.rangeDb = -3.7;
    a.prm.eqBands[2].dyn.thrDb  = -27.3;
    a.prm.eqBands[2].dyn.thrAuto = 0;
    a.prm.eqBands[2].dyn.atk    = 0.37;
    a.prm.eqBands[2].dyn.rel    = 0.63;
    a.prm.inputGainDb        = -2.7;
    a.prm.preLimiterGainDb   =  1.3;
    a.prm.eqBands[0].on      = 1;
    a.prm.eqBands[0].type    = FC_FILTER_HIGH_PASS;
    a.prm.eqBands[0].lanes[0].on = 1;
    a.prm.eqBands[0].lanes[0].freq = 31.7;
    a.prm.eqBands[0].lanes[0].slope = 12;
    a.prm.eqBands[1].on      = 1;
    a.prm.eqBands[1].type    = FC_FILTER_BELL;
    a.prm.eqBands[1].lanes[0].on = 1;
    a.prm.eqBands[1].lanes[0].freq = 2137.3;
    a.prm.eqBands[1].lanes[0].q = 1.7;
    a.prm.eqBands[1].lanes[0].gainDb = -3.1;
    a.prm.monoBass.enabled     = 1;
    a.prm.monoBass.frequencyHz = 123.7f;
    a.prm.monoBass.lowWidth    = 0.3f;
    a.prm.compressor.thresholdDb = -17.3;
    a.prm.compressor.ratio       = 2.7;
    a.prm.compressor.kneeDb      = 4.3;
    a.prm.compressor.attackMs    = 7.3;
    a.prm.compressor.releaseMs   = 137.0;
    a.prm.clipper.driveDb        = 3.7f;
    // NOT at their defaults, because a field left at its default lets a MISSING mapping pass this test:
    // the mutation stand dropped `clipper.mix` and a bypass flag and both survived a suite in which
    // every one of those fields still held the value the default writer had put there.
    a.prm.clipper.mix            = 0.83f;
    a.prm.clipper.bias           = 0.07f;
    a.prm.clipper.outputDb       = -0.7f;
    a.prm.clipper.autoComp       = 0.37f;
    a.prm.clipper.dcBlockHz      = 13.0f;
    // NO BYPASS FLAG IS SET IN THE MAIN SET, and that is a correction rather than an omission. An
    // earlier draft set `bypassMonoBass` here to cover the bypass mapping — and thereby DISABLED the
    // stage whose snap-versus-ramp is the only thing the lifecycle-order check can see, so the check
    // that had caught `prepare`-then-`setParams` stopped catching it. A fixture that covers one thing
    // by switching off another is not coverage. The bypass flags get their own comparison below.
    a.prm.compressor.rangeDb     = 37.0;
    a.prm.compressor.makeupDb    = 1.7;
    a.prm.compressor.rmsWindowMs = 7.7;
    a.prm.limiter.ceilingDbTp    = -1.3;
    a.prm.limiter.releaseMs      = 77.0;
    a.prm.dither.bits            = 24;
    a.prm.dither.seedLo          = 0x748fea9bu;
    a.prm.dither.seedHi          = 0x853c49e6u;

    std::printf ("fcore_master selftest — %g Hz, %d ch, %zu frames\n", fs, nc, frames);

    // --- 1. THE ACCEPTANCE: the ABI and a direct C++ call, one binary, one machine ------------------
    std::vector<float> viaAbi, viaCpp;
    fc_master_resolved res {};
    const bool ranAbi = abiRender (a, in, frames, nc, viaAbi, res);
    const bool ranCpp = directRender (a, in, frames, nc, viaCpp);
    check (ranAbi && ranCpp, "both paths rendered");
    if (ranAbi && ranCpp)
    {
        double worst = 0.0;
        const std::size_t d = bitDiff (viaAbi, viaCpp, worst);
        char msg[128];
        std::snprintf (msg, sizeof msg, "%zu of %zu samples differ, worst %.9g", d, viaAbi.size(), worst);
        check (d == 0, "ABI render is BIT-IDENTICAL to the direct C++ render", msg);

        // A fixture on which the chain does nothing would pass the line above for the wrong reason.
        double amp = 0.0, delta = 0.0;
        for (std::size_t i = 0; i < viaAbi.size(); ++i)
        {
            amp = std::max (amp, (double) std::fabs (viaAbi[i]));
            delta = std::max (delta, (double) std::fabs (viaAbi[i] - in[i]));
        }
        char pre[128];
        std::snprintf (pre, sizeof pre, "output peaks at %.4f and differs from the input by up to %.4f", amp, delta);
        check (amp > 0.05 && delta > 0.01, "PRECONDITION: the chain actually did something", pre);
    }

    // --- 1b. THE OUT-OF-PLACE TRANSPORT PATH --------------------------------------------------------
    // The one a page actually uses, and the one the in-place calls above cannot see. Deleting the
    // facade's `memcpy` left every other check in this binary green.
    {
        std::vector<float> viaCopy; fc_master_resolved r2 {};
        if (abiRender (a, in, frames, nc, viaCopy, r2, /*inPlace*/ false))
        {
            double worst = 0.0;
            const std::size_t d = bitDiff (viaAbi, viaCopy, worst);
            char msg[160];
            std::snprintf (msg, sizeof msg, "%zu differ, worst %.9g", d, worst);
            check (d == 0, "in != out gives the same render as in == out", msg);
        }
        else check (false, "the out-of-place render ran");
    }

    // --- 1c. THE BYPASS FLAGS, in their own comparison ----------------------------------------------
    // Separate from the main set on purpose: setting a bypass flag in the main fixture switches OFF the
    // stage whose behaviour the lifecycle-order check depends on, so covering the flags there quietly
    // removed a check instead of adding one.
    {
        Args b = a;
        b.prm.bypassEq = 1; b.prm.bypassMonoBass = 1; b.prm.bypassCompressor = 1;
        b.prm.bypassClipper = 1; b.prm.bypassLimiter = 1; b.prm.bypassDither = 1;
        std::vector<float> bAbi, bCpp; fc_master_resolved rb {};
        const bool ranB = abiRender (b, in, frames, nc, bAbi, rb) && directRender (b, in, frames, nc, bCpp);
        check (ranB, "the all-bypassed render ran through both paths");
        if (ranB)
        {
            double worst = 0.0;
            const std::size_t d = bitDiff (bAbi, bCpp, worst);
            char msg[160];
            std::snprintf (msg, sizeof msg, "%zu differ, worst %.9g", d, worst);
            check (d == 0, "every bypass flag maps identically through the ABI", msg);
            // And the flags must actually DO something, or this comparison passes for the wrong reason.
            double delta = 0.0;
            for (std::size_t i = 0; i < bAbi.size(); ++i)
                delta = std::max (delta, (double) std::fabs (bAbi[i] - viaAbi[i]));
            char pre[128];
            std::snprintf (pre, sizeof pre, "bypassed vs active differ by up to %.4f", delta);
            check (delta > 0.01, "PRECONDITION: the bypass flags changed the render", pre);
        }
    }

    // --- 2. BLOCK INDEPENDENCE THROUGH THE ABI ------------------------------------------------------
    // The chain's fixed internal quantum is what makes this a theorem rather than a hope; this checks
    // that the ABI does not route around it, which it would the moment it re-blocked anything itself.
    for (std::uint32_t blk : { 1u, 337u, 4096u, (std::uint32_t) frames })
    {
        Args b = a; b.block = blk;
        std::vector<float> other; fc_master_resolved r2 {};
        if (! abiRender (b, in, frames, nc, other, r2)) { check (false, "block-size render"); continue; }
        double worst = 0.0;
        const std::size_t d = bitDiff (viaAbi, other, worst);
        char msg[160];
        std::snprintf (msg, sizeof msg, "block %u: %zu differ, worst %.9g", blk, d, worst);
        check (d == 0, "block-independent through the ABI", msg);
    }

    // --- 3. THE RESET DISCIPLINE --------------------------------------------------------------------
    // A second programme through the same handle WITHOUT a reset is a different render; with one it is
    // the same. Both halves are asserted, because only the pair says the reset is what did it.
    {
        fc_master h = 0; fc_master_resolved r3 {}; FC_INIT (r3);
        std::vector<float> first, second, third;
        const fc_status stC = fc_master_create (&a.cfg, &h);
        const fc_status stK = (stC == FC_OK) ? fc_master_configure (h, &a.prm, &r3) : stC;
        const bool created = (stC == FC_OK && stK == FC_OK);
        char cm[96];
        std::snprintf (cm, sizeof cm, "create=%s configure=%s", statusName (stC), statusName (stK));
        check (created, "handle for the reset check", cm);
        if (created)
        {
            auto runOnce = [&] (std::vector<float>& outv)
            {
                const std::size_t D = (std::size_t) r3.latencySamples;
                std::vector<float> work (in);
                outv.assign (frames * (std::size_t) nc, 0.0f);
                fc_master_process (h, work.data(), work.data(), (std::uint32_t) frames);
                std::vector<float> tail (D * (std::size_t) nc, 0.0f);
                std::uint32_t got = 0;
                fc_master_flush (h, tail.data(), (std::uint32_t) D, &got);
                for (int c = 0; c < nc; ++c)
                    for (std::size_t n = 0; n < frames; ++n)
                        outv[(std::size_t) c * frames + n] = (n + D < frames)
                            ? work[(std::size_t) c * frames + n + D]
                            : tail[(std::size_t) c * D + (n + D - frames)];
            };
            runOnce (first);
            runOnce (second);                                   // no reset
            fc_master_reset (h);
            runOnce (third);                                    // after a reset
            double w1 = 0.0, w2 = 0.0;
            const std::size_t d1 = bitDiff (first, second, w1);
            const std::size_t d2 = bitDiff (first, third, w2);
            char m1[160], m2[160];
            std::snprintf (m1, sizeof m1, "%zu differ, worst %.9g", d1, w1);
            std::snprintf (m2, sizeof m2, "%zu differ, worst %.9g", d2, w2);
            check (d1 > 0, "a second programme WITHOUT reset really is a different render", m1);
            check (d2 == 0, "and identical after fc_master_reset", m2);
            fc_master_destroy (h);
        }
    }

    // --- 4. THE HANDLE IS NOT A POINTER --------------------------------------------------------------
    {
        fc_master h = 0;
        check (fc_master_create (&a.cfg, &h) == FC_OK && h != 0, "create issues a non-zero handle");
        check (fc_master_destroy (h) == FC_OK, "destroy accepts it once");
        check (fc_master_destroy (h) == FC_ERR_HANDLE, "and refuses it the second time");
        fc_master again = 0;
        check (fc_master_create (&a.cfg, &again) == FC_OK, "a new handle after the destroy");
        check (again != h, "which is NOT the old value, whatever the allocator did with the memory");
        std::int32_t lat = 0;
        check (fc_master_latency (h, &lat) == FC_ERR_HANDLE, "the stale handle is refused, not aliased");
        fc_master_destroy (again);
    }

    std::printf ("%s — %d failure(s)\n", failures == 0 ? "PASS" : "FAIL", failures);
    return failures == 0 ? 0 : 1;
}

}   // namespace

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf (stderr,
            "usage:\n"
            "  %s render   <sampleRate> <channels> <in.f32le> <out.f32le> [key=value ...]\n"
            "  %s solve    <sampleRate> <channels> <in.f32le> <out.f32le> target=<LUFS> tp=<dBTP> [key=value ...]\n"
            "  %s lra      <sampleRate> <channels> <in.f32le>\n"
            "  %s selftest [sampleRate] [channels]\n", argv[0], argv[0], argv[0], argv[0]);
        return 2;
    }
    const std::string mode = argv[1];

    if (mode == "selftest")
    {
        const double fs = argc > 2 ? std::atof (argv[2]) : 48000.0;
        const int    nc = argc > 3 ? std::atoi (argv[3]) : 2;
        return selftest (fs, nc);
    }

    if (mode == "lra")
    {
        if (argc < 5) { std::fprintf (stderr, "lra needs <sampleRate> <channels> <in.f32le>\n"); return 2; }
        Args a;
        if (! parseArgs (a, argc, argv, 5)) return 2;
        a.cfg.sampleRate = std::atof (argv[2]);
        a.cfg.channels   = std::atoi (argv[3]);
        // BEFORE the read, not after: `readInterleaved` divides by the channel count, so a zero here is
        // an integer division by zero — undefined behaviour where the contract promises a refusal.
        if (a.cfg.channels < 1 || a.cfg.channels > core::kMaxChannels)
        { std::fprintf (stderr, "bad channel count\n"); return 2; }
        std::vector<float> in; std::size_t frames = 0;
        if (! readInterleaved (argv[4], a.cfg.channels, in, frames)) return 2;
        if (frames > (std::size_t) 0x7FFFFFFFu)
        { std::fprintf (stderr, "input longer than the ABI's frame count\n"); return 2; }
        fc_master h = 0;
        if (const fc_status st = fc_master_create (&a.cfg, &h); st != FC_OK)
        { std::fprintf (stderr, "create: %s\n", statusName (st)); return 2; }
        double lra = 0.0;
        const fc_status st = fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &lra);
        fc_master_destroy (h);
        if (st != FC_OK) { std::fprintf (stderr, "lra: %s\n", statusName (st)); return 1; }
        std::printf ("%.17g\n", lra);
        return 0;
    }

    if (argc < 6) { std::fprintf (stderr, "need <sampleRate> <channels> <in.f32le> <out.f32le>\n"); return 2; }
    Args a;
    if (! parseArgs (a, argc, argv, 6)) return 2;
    a.cfg.sampleRate = std::atof (argv[2]);
    a.cfg.channels   = std::atoi (argv[3]);
    const int nc = a.cfg.channels;
    if (nc < 1 || nc > core::kMaxChannels) { std::fprintf (stderr, "bad channel count\n"); return 2; }

    std::vector<float> in; std::size_t frames = 0;
    if (! readInterleaved (argv[4], nc, in, frames)) return 2;
    if (frames == 0) { std::fprintf (stderr, "empty input\n"); return 2; }
    // The ABI counts frames in 32 bits. Narrowing silently would not merely lose length: the planar
    // STRIDE is the frame count, so channel 1's pointer would land inside channel 0's plane and the
    // render would be of a signal that does not exist.
    if (frames > (std::size_t) 0x7FFFFFFFu)
    { std::fprintf (stderr, "input longer than the ABI's frame count (%zu frames)\n", frames); return 2; }

    if (mode == "render")
    {
        std::vector<float> out; fc_master_resolved res {};
        if (! abiRender (a, in, frames, nc, out, res)) return 1;
        if (! writeInterleaved (argv[5], nc, out, frames)) return 1;
        // The resolved geometry on stderr, so a harness can read the numbers without parsing the audio.
        std::fprintf (stderr, "latency=%d internalBlock=%d ceiling=%.17g release=%.17g\n",
                      res.latencySamples, res.internalBlock, res.limiterCeilingDbTp, res.limiterReleaseMs);
        return 0;
    }

    if (mode == "solve")
    {
        fc_master h = 0;
        if (const fc_status st = fc_master_create (&a.cfg, &h); st != FC_OK)
        { std::fprintf (stderr, "create: %s\n", statusName (st)); return 2; }
        std::vector<float> out (in.size(), 0.0f);
        fc_solution sol = 0;
        const fc_status st = fc_master_solve (h, &a.prm, &a.req, in.data(), out.data(),
                                              (std::uint32_t) frames, &sol);
        if (st != FC_OK) { std::fprintf (stderr, "solve: %s\n", statusName (st)); fc_master_destroy (h); return 1; }

        fc_solution_summary sum {}; FC_INIT (sum);
        fc_measurement meas {};      FC_INIT (meas);
        fc_solution_summary_get (sol, &sum);
        fc_solution_measurement (sol, &meas);
        std::printf ("status=%d binding=%d gain=%.17g ceiling=%.17g passes=%d\n",
                     sum.status, sum.binding, sum.preLimiterGainDb, sum.ceilingDbTp, sum.passes);
        std::printf ("I=%.17g TP=%.17g LRA=%.17g PLR=%.17g loudnessValid=%d lraValid=%d\n",
                     meas.integratedLufs, meas.truePeakDbTp, meas.loudnessRangeLu, meas.plrDb,
                     meas.loudnessValid, meas.lraValid);
        // A VERDICT IS NOT A RENDER. `InvalidRequest` and `NotPrepared` are returned before the solver
        // touches the output, so `out` is still the zero buffer it was allocated as — writing it would
        // hand a harness a file of the right length, full of digital silence, with exit status 0 and an
        // existing output overwritten. The status line has already been printed; the file is not.
        const bool delivered = (sum.status != FC_SOLVE_INVALID_REQUEST
                             && sum.status != FC_SOLVE_NOT_PREPARED
                             && sum.status != FC_SOLVE_RENDER_FAILED);
        bool ok = true;
        if (delivered) ok = writeInterleaved (argv[5], nc, out, frames);
        else std::fprintf (stderr, "no render was delivered (status=%d) — the output file is NOT written\n",
                           sum.status);
        fc_solution_destroy (sol);
        fc_master_destroy (h);
        return (delivered && ok) ? 0 : 1;
    }

    std::fprintf (stderr, "unknown mode '%s'\n", mode.c_str());
    return 2;
}
