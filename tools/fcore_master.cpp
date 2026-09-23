// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_master — the native CLI that drives the C ABI (tools/fc_master_abi.h, implemented in
// tools/wasm/fc_master.cpp), and the reference that proves the ABI adds nothing.
//
//   fcore_master render   <sampleRate> <channels> <in.f32le> <out.f32le> [key=value ...]
//   fcore_master solve    <sampleRate> <channels> <in.f32le> <out.f32le> target=<LUFS> tp=<dBTP> [key=value ...]
//   fcore_master lra      <sampleRate> <channels> <in.f32le>
//   fcore_master selftest [sampleRate] [channels]
//   fcore_master layout                                  — the ABI's struct offsets as JSON, for layout-check.mjs
//
// `delivery=<rate>` on render / solve / lra makes a DELIVERING handle (ABI v2): SRC first, and the output file
// holds the DELIVERED length at <rate> — `fc_master_delivered_frames`, not the input's frame count.
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
//     samples apart — so both paths start from the ABI's own defaults writer (`fc_*_defaults`, the versioned one:
//     the frozen v1 `fc_*_default` would stamp v1 and hide every newer field from the facade).
//   * A RESET BEFORE EVERY RENDER. `OfflineRenderer` resets the chain itself; the ABI has no renderer,
//     so a second programme through the same handle without `fc_master_reset` is a different render —
//     the check below prints the count for the fixture it actually ran, rather than carrying a number
//     from a fixture that has since changed. The loop resets, and the harness must too.

#include "fc_master_abi.h"

#include <felitronics/mastering/DeliveredMastering.h>
#include <felitronics/mastering/LoudnessSolver.h>
#include <felitronics/mastering/MasteringChain.h>
#include <felitronics/mastering/OfflineRenderer.h>

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <string>
#include <system_error>
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
        case FC_ERR_POISONED:        return "POISONED";
        case FC_ERR_CANCELLED:       return "CANCELLED";
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
    struct Weight { std::int32_t channel; double value; };
    std::vector<Weight> weights;
};

// THE CURRENT LAYOUT'S DEFAULTS, through the VERSIONED writers. Not `fc_*_default`: those are frozen at v1 and
// stamp v1, so a `delivery=` key written after them would lie past the stamped size and never be read — the
// handle would not convert, and the file would come out at the source rate with exit status 0.
bool initArgs (Args& a)
{
    FC_INIT (a.cfg); FC_INIT (a.prm); FC_INIT (a.req);
    return fc_master_config_defaults (&a.cfg) == FC_OK
        && fc_master_params_defaults (&a.prm) == FC_OK
        && fc_loudness_request_defaults (&a.req) == FC_OK;
}

// STRICT number parsing. `atof`/`atoi` stop at the first character they do not understand and report
// nothing, so `inputGainDb=oops` becomes 0 dB and `comp.ratio=4oops` becomes 4 — a harness row computed
// from settings nobody chose, which is the whole failure mode this CLI's strictness exists to prevent.
bool num (const std::string& v, double& out)
{
    if (v.empty()) return false;
    char* end = nullptr;
    errno = 0;
    const double d = std::strtod (v.c_str(), &end);
    if (end == v.c_str() || *end != '\0' || errno == ERANGE) return false;
    out = d;
    return true;
}

// `errno` IS PART OF THE CHECK, and leaving it out made this guard TIER-DEPENDENT — which is exactly
// the class the wasm row exists to catch, and it caught this one. `long` is 64 bits here and 32 bits on
// wasm32, so `strtol("4294967296")` returns the value on the desktop row and SATURATES to LONG_MAX on
// wasm32 with ERANGE set. Without the errno test the saturated value then passed the range check that
// the guard is, and the block-of-zero hang the guard exists to prevent came back — on the tier the
// facade is written for, and only there.
bool inum (const std::string& v, long& out)
{
    if (v.empty()) return false;
    char* end = nullptr;
    errno = 0;
    const long i = std::strtol (v.c_str(), &end, 10);
    if (end == v.c_str() || *end != '\0' || errno == ERANGE) return false;
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
    // RANGE BEFORE THE CAST. `(int) 4294967296` is 0, so a band index far out of range would land on
    // band 0 and quietly apply the value to the wrong band — the same shape as the `block` hang above.
    if (idxL < 0 || idxL >= FC_MAX_EQ_BANDS) return false;
    const int idx = (int) idxL;
    const std::string f = key.substr (dot + 1);
    fc_eq_band& b = a.prm.eqBands[idx];
    static const char* kTypes[] { "bell", "lowshelf", "highshelf", "highpass", "lowpass",
                                  "bandpass", "notch", "allpass", "tilt" };
    if (f == "on")     return parseBool (val, b.on);
    if (f == "bypass") return parseBool (val, b.bypass);
    if (f == "swept")  return parseBool (val, b.swept);
    if (f == "type")   return parseEnumName (val, kTypes, 9, b.type);
    if (f == "laneon") return parseBool (val, b.lanes[0].on);
    // The point's dynamics. `dyn.range` carries the DIRECTION in its sign: negative cuts as the band
    // gets loud, positive boosts. 0 is no dynamics whatever `dyn.on` says.
    if (f == "dyn.on")      return parseBool (val, b.dyn.on);
    if (f == "dyn.thrAuto") return parseBool (val, b.dyn.thrAuto);
    double d = 0.0; long i = 0;
    if (f == "dyn.range") { if (! num (val, d)) return false; b.dyn.rangeDb = d; return true; }
    if (f == "dyn.thr")   { if (! num (val, d)) return false; b.dyn.thrDb   = d; return true; }
    if (f == "dyn.atk")   { if (! num (val, d)) return false; b.dyn.atk     = d; return true; }
    if (f == "dyn.rel")   { if (! num (val, d)) return false; b.dyn.rel     = d; return true; }
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

    // A `block` that does not fit a positive int is REFUSED. `(uint32_t) 4294967296` is 0, and a block
    // of zero makes the render loop stop advancing — a hang, out of a typo, with no message.
    if (key == "block")          { if (! iOk || i < 1 || i > 0x7FFFFFFF) return false;
                                   a.block = (std::uint32_t) i; return true; }
    if (key == "internalBlock")  { FC_I (a.cfg.internalBlock = i); }
    if (key == "oversample")     { FC_I (a.cfg.oversampleFactor = i); }
    if (key == "taps")           { FC_I (a.cfg.tapsPerPhase = i); }
    if (key == "compLookaheadMs"){ FC_D (a.cfg.compressorLookaheadMs = d); }
    if (key == "limLookaheadMs") { FC_D (a.cfg.limiterLookaheadMs = d); }
    if (key == "sidechainHpfHz") { FC_D (a.cfg.sidechainHpfHz = d); }
    if (key == "delivery")       { FC_D (a.cfg.deliveryRate = d); }
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
    static const char* kGrSt[]  { "mean", "p95", "max", "percentile" };

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
    if (key == "comp.mix")       { FC_D (a.prm.compressorMix = d); }

    if (key == "clip.shape")  return parseEnumName (val, kShape, 4, a.prm.clipper.shape);
    if (key == "clip.drive")  { FC_D (a.prm.clipper.driveDb  = (float) d); }
    if (key == "clip.bias")   { FC_D (a.prm.clipper.bias     = (float) d); }
    if (key == "clip.mix")    { FC_D (a.prm.clipper.mix      = (float) d); }
    if (key == "clip.output") { FC_D (a.prm.clipper.outputDb = (float) d); }
    if (key == "clip.autoComp") { FC_D (a.prm.clipper.autoComp = (float) d); }
    if (key == "clip.dcBlockHz"){ FC_D (a.prm.clipper.dcBlockHz = (float) d); }

    if (key == "lim.ceiling") { FC_D (a.prm.limiter.ceilingDbTp = d); }
    if (key == "lim.release") { FC_D (a.prm.limiter.releaseMs = d); }
    if (key == "lim.dual")    return parseBool (val, a.prm.limiterDualRelease);
    if (key == "lim.slowRelease") { FC_D (a.prm.limiterSlowReleaseMs = d); }
    // K13 — the peak clipper inside the limiter's island. `over` is dB ABOVE the ceiling, not a level.
    if (key == "lim.peakClip")     return parseBool (val, a.prm.peakClipper);
    if (key == "lim.clipOver")     { FC_D (a.prm.peakClipperOverCeilingDb = d); }
    if (key == "lim.clipKnee")     { FC_D (a.prm.peakClipperKneeDb = d); }

    if (key == "dith.bits")    { FC_I (a.prm.dither.bits = i); }
    if (key == "dith.shaping") return parseEnumName (val, kShap, 3, a.prm.dither.shaping);
    if (key == "dith.seedLo")  { char* e = nullptr; errno = 0; const unsigned long u = std::strtoul (val.c_str(), &e, 0);
                                 if (e == val.c_str() || *e != '\0' || errno == ERANGE) return false;
                                 a.prm.dither.seedLo = (std::uint32_t) u; return true; }
    if (key == "dith.seedHi")  { char* e = nullptr; errno = 0; const unsigned long u = std::strtoul (val.c_str(), &e, 0);
                                 if (e == val.c_str() || *e != '\0' || errno == ERANGE) return false;
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
    if (key == "limGrStat")  return parseEnumName (val, kGrSt, 4, a.req.limiterGr.statistic);
    if (key == "compGrLimit"){ FC_D (a.req.compressorGr.limitDb = d); }
    if (key == "compGrStat") return parseEnumName (val, kGrSt, 4, a.req.compressorGr.statistic);
    if (key == "limActiveDb"){ FC_D (a.req.limiterActiveInputDb = d); }   // v10 — K11's gate on the limiter's input
    if (key == "limGrQ")     { FC_D (a.req.limiterGrQuantile = d); }
    if (key == "compGrQ")    { FC_D (a.req.compressorGrQuantile = d); }
    if (key == "activityDb") { FC_D (a.req.activityThresholdDb = d); }
    if (key == "grTraceBuckets") { if (! iOk || i < -0x7FFFFFFFL - 1L || i > 0x7FFFFFFFL) return false;
                                   a.req.grTraceBuckets = (std::int32_t) i; return true; }

    // `weight<N>=<w>` — the BS.1770 channel weights. The ABI grew an entry point for them and nothing
    // called it, which made the capability reachable in principle and not in practice: a surround run
    // through this CLI was measured at 1.0 everywhere, where the standard says Ls/Rs 1.41 and LFE 0.
    if (key.rfind ("weight", 0) == 0)
    {
        long ch = 0;
        if (! inum (key.substr (6), ch) || ! dOk) return false;
        if (ch < 0 || ch >= core::kMaxChannels) return false;
        a.weights.push_back ({ (std::int32_t) ch, d });
        return true;
    }
    #undef FC_D
    #undef FC_I
    return false;
}

// The POSITIONAL arguments deserve the same strictness as the keys: `fcore_master render 48000oops 2 …`
// used to become 48000, and `2oops` became 2, so the whole run described a geometry nobody typed.
bool positional (const char* text, double& out) { return num (std::string (text), out); }
bool positional (const char* text, long& out)   { return inum (std::string (text), out); }

bool parseArgs (Args& a, int argc, char** argv, int from)
{
    if (! initArgs (a)) { std::fprintf (stderr, "the ABI refused its own defaults\n"); return false; }
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
    // A TRAILING PARTIAL FRAME IS NOT A SHORT FILE EITHER. Integer division swallowed it, so a file of
    // 8193 floats at two channels rendered its first 4096 frames and exited 0 — a harness row for a
    // programme that is not the one on disk. Same rule as everywhere else here: refuse, do not truncate.
    if (inter.size() % (std::size_t) nc != 0)
    { std::fprintf (stderr, "input is not a whole number of %d-channel frames\n", nc); return false; }
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

    for (const auto& w : a.weights)
        if (const fc_status st = fc_master_set_channel_weight (h, w.channel, w.value); st != FC_OK)
        { std::fprintf (stderr, "weight%d: %s\n", w.channel, statusName (st)); fc_master_destroy (h); return false; }

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

// THE DELIVERED ABI RENDER (v2) — create -> configure -> fc_master_delivered_frames -> fc_master_render_delivered.
// There is no block loop and no latency drop here, and that is the point of the entry point: the delivered
// length, the trim and the drain are the core's, and the output is `outFrames` long at the delivery rate.
bool abiRenderDelivered (const Args& a, const std::vector<float>& in, std::size_t frames, int nc,
                         std::vector<float>& out, std::size_t& outFrames, fc_master_resolved& res)
{
    fc_master h = 0;
    if (const fc_status st = fc_master_create (&a.cfg, &h); st != FC_OK)
    { std::fprintf (stderr, "create: %s\n", statusName (st)); return false; }
    auto fail = [h] (const char* what, fc_status st) { std::fprintf (stderr, "%s: %s\n", what, statusName (st));
                                                       fc_master_destroy (h); return false; };
    for (const auto& w : a.weights)
        if (const fc_status st = fc_master_set_channel_weight (h, w.channel, w.value); st != FC_OK) return fail ("weight", st);
    FC_INIT (res);
    if (const fc_status st = fc_master_configure (h, &a.prm, &res); st != FC_OK) return fail ("configure", st);
    std::uint32_t d = 0;
    if (const fc_status st = fc_master_delivered_frames (h, (std::uint32_t) frames, &d); st != FC_OK)
        return fail ("delivered_frames", st);
    outFrames = d;
    out.assign ((std::size_t) d * (std::size_t) nc, 0.0f);
    if (const fc_status st = fc_master_render_delivered (h, in.data(), (std::uint32_t) frames, out.data(), d); st != FC_OK)
        return fail ("render_delivered", st);
    fc_master_destroy (h);
    return true;
}

// The ABI's config and parameters, mirrored by hand into the core's — see the note inside.
void mirror (const Args& a, MasteringChainConfig& cc, MasteringChainParams& cp)
{
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
    cp.limiter.dualRelease = a.prm.limiterDualRelease != 0;
    cp.limiter.slowReleaseMs = a.prm.limiterSlowReleaseMs;
    cp.limiter.peakClip      = a.prm.peakClipper != 0;
    cp.limiter.overCeilingDb = a.prm.peakClipperOverCeilingDb;
    cp.limiter.kneeDb        = a.prm.peakClipperKneeDb;
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
    cp.compressorMix = a.prm.compressorMix;
}

// The same render through the C++ API, with the SAME lifecycle order the ABI takes (parameters written
// BEFORE prepare, which is what `fc_master_configure` does when it re-prepares).
bool directRender (const Args& a, const std::vector<float>& in, std::size_t frames, int nc,
                   std::vector<float>& out)
{
    MasteringChainConfig cc {};
    MasteringChainParams cp {};
    mirror (a, cc, cp);

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

// The delivered render through the C++ API: the chain at the DELIVERY rate, and `mastering::DeliveredMastering`
// — the very class the facade forwards to — in front of it. What the selftest compares is therefore the ABI's
// marshalling against a core call, not a hand-written conversion against another one. The renderer and the
// converter run at `a.block` here and at the facade's own 4096 there; both are free by contract, and a
// difference would be a real one.
bool directRenderDelivered (const Args& a, const std::vector<float>& in, std::size_t frames, int nc,
                            std::vector<float>& out, std::size_t& outFrames)
{
    MasteringChainConfig cc {};
    MasteringChainParams cp {};
    mirror (a, cc, cp);

    MasteringChain chain;
    OfflineRenderer r;
    DeliveredMastering dm;
    if (! r.prepare (nc, (int) a.block)) return false;
    chain.setParams (cp);
    if (! chain.prepare (a.cfg.deliveryRate, nc, cc)) return false;
    if (! dm.prepare (a.cfg.sampleRate, a.cfg.deliveryRate, nc, (int) a.block)) return false;
    const long long d = DeliveredMastering::deliveredFrames (a.cfg.sampleRate, a.cfg.deliveryRate, (long long) frames);
    if (d < 0) return false;
    outFrames = (std::size_t) d;
    out.assign (outFrames * (std::size_t) nc, 0.0f);
    const float* ip[core::kMaxChannels] {};
    float*       op[core::kMaxChannels] {};
    for (int c = 0; c < nc; ++c)
    {
        ip[c] = in.data() + (std::size_t) c * frames;
        op[c] = out.data() + (std::size_t) c * outFrames;
    }
    return dm.render (chain, r, ip, nc, (long long) frames, op, d);
}

//==============================================================================
// `fcore_master layout` — EVERY FIELD OF EVERY STRUCT the JS layout describes, at the offset THIS compiler put it,
// for tools/wasm/layout-check.mjs to hold fc-master-layout.mjs against. A total size alone cannot see two fields
// of one type swapped in the JS list, and that permutation writes a page's value into the wrong knob.
//
// The list is a transcription and says so; what keeps it honest is that `offsetof` does not compile on a field
// that is not there, and that the checker demands the SAME SET of names in both directions — a field the JS has
// and this list lacks is reported, and so is the reverse. Nested fields are listed under their own struct.
#define FC_LAYOUT_FIELDS(X)                                                                                         \
    X (fc_header, abiVersion) X (fc_header, structSize)                                                            \
    X (fc_master_config, header) X (fc_master_config, sampleRate) X (fc_master_config, channels)                   \
    X (fc_master_config, internalBlock) X (fc_master_config, eq) X (fc_master_config, monoBass)                    \
    X (fc_master_config, compressor) X (fc_master_config, clipper) X (fc_master_config, limiter)                   \
    X (fc_master_config, dither) X (fc_master_config, compressorLookaheadMs)                                       \
    X (fc_master_config, limiterLookaheadMs) X (fc_master_config, oversampleFactor)                                \
    X (fc_master_config, tapsPerPhase) X (fc_master_config, sidechainHpfHz) X (fc_master_config, deliveryRate)      \
    X (fc_eq_lane, on) X (fc_eq_lane, freq) X (fc_eq_lane, q) X (fc_eq_lane, gainDb) X (fc_eq_lane, slope)         \
    X (fc_eq_lane, bypass)                                                                                         \
    X (fc_eq_dyn, on) X (fc_eq_dyn, rangeDb) X (fc_eq_dyn, thrDb) X (fc_eq_dyn, thrAuto) X (fc_eq_dyn, atk)        \
    X (fc_eq_dyn, rel)                                                                                             \
    X (fc_eq_band, on) X (fc_eq_band, type) X (fc_eq_band, swept) X (fc_eq_band, bypass) X (fc_eq_band, dyn)       \
    X (fc_eq_band, lanes)                                                                                          \
    X (fc_mono_bass, enabled) X (fc_mono_bass, frequencyHz) X (fc_mono_bass, lowWidth)                             \
    X (fc_compressor, detector) X (fc_compressor, link) X (fc_compressor, rmsWindowMs) X (fc_compressor, mode)     \
    X (fc_compressor, thresholdDb) X (fc_compressor, ratio) X (fc_compressor, kneeDb) X (fc_compressor, rangeDb)   \
    X (fc_compressor, attackMs) X (fc_compressor, releaseMs) X (fc_compressor, makeupDb)                           \
    X (fc_compressor, autoMakeup)                                                                                  \
    X (fc_clipper, shape) X (fc_clipper, driveDb) X (fc_clipper, bias) X (fc_clipper, mix)                         \
    X (fc_clipper, outputDb) X (fc_clipper, autoComp) X (fc_clipper, dcBlockHz)                                    \
    X (fc_limiter, ceilingDbTp) X (fc_limiter, releaseMs)                                                          \
    X (fc_dither, bits) X (fc_dither, shaping) X (fc_dither, seedLo) X (fc_dither, seedHi)                         \
    X (fc_dither, autoBlank) X (fc_dither, autoBlankSamples)                                                       \
    X (fc_master_params, header) X (fc_master_params, inputGainDb) X (fc_master_params, preLimiterGainDb)          \
    X (fc_master_params, eqBands) X (fc_master_params, monoBass) X (fc_master_params, compressor)                  \
    X (fc_master_params, clipper) X (fc_master_params, limiter) X (fc_master_params, dither)                       \
    X (fc_master_params, bypassEq) X (fc_master_params, bypassMonoBass) X (fc_master_params, bypassCompressor)     \
    X (fc_master_params, bypassClipper) X (fc_master_params, bypassLimiter) X (fc_master_params, bypassDither)     \
    X (fc_master_params, compressorMix) X (fc_master_params, limiterDualRelease) X (fc_master_params, _pad0)       \
    X (fc_master_params, limiterSlowReleaseMs)                                                                     \
    X (fc_master_params, peakClipper) X (fc_master_params, _pad1)                                                  \
    X (fc_master_params, peakClipperOverCeilingDb) X (fc_master_params, peakClipperKneeDb)                         \
    X (fc_master_resolved, header) X (fc_master_resolved, latencySamples) X (fc_master_resolved, internalBlock)    \
    X (fc_master_resolved, compressorLookahead) X (fc_master_resolved, clipperLatency)                             \
    X (fc_master_resolved, limiterLatency) X (fc_master_resolved, limiterLookahead)                                \
    X (fc_master_resolved, oversampleFactor) X (fc_master_resolved, compressorTapOffset)                           \
    X (fc_master_resolved, limiterTapOffset) X (fc_master_resolved, limiterCeilingDbTp)                            \
    X (fc_master_resolved, limiterReleaseMs) X (fc_master_resolved, monoBass)                                      \
    X (fc_master_resolved, tapOversampleFactor) X (fc_master_resolved, compressorMix)                              \
    X (fc_master_resolved, limiterSlowReleaseMs) X (fc_master_resolved, peakClipperThresholdDbTp)                  \
    X (fc_master_stats, header) X (fc_master_stats, framesIn) X (fc_master_stats, framesFlushed)                   \
    X (fc_master_stats, nonFiniteIn)                                                                               \
    X (fc_need, header) X (fc_need, callBytes) X (fc_need, solverPrepareBytes) X (fc_need, facadeBytes)            \
    X (fc_need, solverPrepared) X (fc_need, _pad0)                                                                 \
    X (fc_gr_limit, limitDb) X (fc_gr_limit, statistic)                                                            \
    X (fc_loudness_request, header) X (fc_loudness_request, targetLufs) X (fc_loudness_request, toleranceLu)       \
    X (fc_loudness_request, maxTruePeakDbTp) X (fc_loudness_request, truePeakAimDb)                                \
    X (fc_loudness_request, limiterGr) X (fc_loudness_request, compressorGr) X (fc_loudness_request, minPlrDb)     \
    X (fc_loudness_request, maxLraLossLu) X (fc_loudness_request, inputLoudnessRangeLu)                            \
    X (fc_loudness_request, activityThresholdDb) X (fc_loudness_request, maxPasses)                                \
    X (fc_loudness_request, initialGainDb) X (fc_loudness_request, grTraceBuckets) X (fc_loudness_request, _pad0)  \
    X (fc_loudness_request, limiterGrQuantile) X (fc_loudness_request, compressorGrQuantile)                        \
    X (fc_loudness_request, limiterActiveInputDb)                                                                  \
    X (fc_gr_active_stats, header) X (fc_gr_active_stats, stats) X (fc_gr_active_stats, windows)                   \
    X (fc_gr_active_stats, activeWindows) X (fc_gr_active_stats, thresholdDb)                                      \
    X (fc_solve_pass, gainDb) X (fc_solve_pass, ceilingDb) X (fc_solve_pass, integratedLufs)                       \
    X (fc_solve_pass, truePeakDbTp) X (fc_solve_pass, plrDb) X (fc_solve_pass, limiterMaxGrDb)                     \
    X (fc_solve_pass, loudnessRangeLu) X (fc_solve_pass, violated)                                                 \
    X (fc_gr_stats, meanDb) X (fc_gr_stats, p95Db) X (fc_gr_stats, maxDb) X (fc_gr_stats, activeFraction)          \
    X (fc_gr_stats, frames) X (fc_gr_stats, nonFinite) X (fc_gr_stats, aboveRange) X (fc_gr_stats, valid)          \
    X (fc_measurement, header) X (fc_measurement, integratedLufs) X (fc_measurement, truePeakDbTp)                 \
    X (fc_measurement, samplePeakDb) X (fc_measurement, loudnessRangeLu) X (fc_measurement, plrDb)                 \
    X (fc_measurement, compressor) X (fc_measurement, limiter) X (fc_measurement, limiterMaxReconstructedPeakDb)    \
    X (fc_measurement, latencySamples) X (fc_measurement, gatingBlocks) X (fc_measurement, droppedBlocks)           \
    X (fc_measurement, nonFiniteSubHops) X (fc_measurement, loudnessValid) X (fc_measurement, lraValid)            \
    X (fc_measurement, compressorGrTraceBuckets) X (fc_measurement, limiterGrTraceBuckets)                         \
    X (fc_measurement, compressorGrTraceValid) X (fc_measurement, limiterGrTraceValid)                             \
    X (fc_measurement, peakClipReductionMaxDb) X (fc_measurement, peakClipReductionP95Db)                          \
    X (fc_measurement, peakClipOccupancy) X (fc_measurement, peakClipRuns)                                         \
    X (fc_measurement, peakClipRunSamplesTotal) X (fc_measurement, peakClipLongestRunSamples)                      \
    X (fc_gr_trace_bucket, maxDb) X (fc_gr_trace_bucket, meanDb) X (fc_gr_trace_bucket, samples)                   \
    X (fc_gr_trace_bucket, nonFinite)                                                                              \
    X (fc_gr_trace_bucket64, maxDb) X (fc_gr_trace_bucket64, meanDb) X (fc_gr_trace_bucket64, samples)             \
    X (fc_gr_trace_bucket64, nonFinite)                                                                            \
    X (fc_progress, stage) X (fc_progress, pass) X (fc_progress, maxPasses) X (fc_progress, hasRecord)             \
    X (fc_progress, fraction) X (fc_progress, record)                                                              \
    X (fc_solution_summary, header) X (fc_solution_summary, status) X (fc_solution_summary, binding)               \
    X (fc_solution_summary, alsoViolated) X (fc_solution_summary, preLimiterGainDb)                                \
    X (fc_solution_summary, ceilingDbTp) X (fc_solution_summary, passes) X (fc_solution_summary, logCount)         \
    X (fc_solution_summary, activityThresholdDb) X (fc_solution_summary, achievedBelowLufs)                        \
    X (fc_solution_summary, achievedAboveLufs) X (fc_solution_summary, gainBelowDb)                                \
    X (fc_solution_summary, gainAboveDb)

// One line per fact: `V <abi version>`, `S <struct> <sizeof>`, `F <struct> <field> <offset>`, and for every struct
// with a header `T <struct> <id> <fc_master_sizeof(id, current)>` — the published table, not this file's sizeof.
int printLayout()
{
    std::printf ("V %u\n", fc_master_abi_version());
    const char* last = "";
#define FC_PRINT(T, f)                                                                              \
    if (std::strcmp (last, #T) != 0) { std::printf ("S " #T " %zu\n", sizeof (T)); last = #T; }     \
    std::printf ("F " #T " " #f " %zu\n", offsetof (T, f));
    FC_LAYOUT_FIELDS (FC_PRINT)
#undef FC_PRINT
    const struct { const char* name; int id; } headered[] {
        { "fc_master_config", FC_STRUCT_CONFIG }, { "fc_master_params", FC_STRUCT_PARAMS },
        { "fc_master_resolved", FC_STRUCT_RESOLVED }, { "fc_master_stats", FC_STRUCT_STATS },
        { "fc_need", FC_STRUCT_NEED }, { "fc_loudness_request", FC_STRUCT_REQUEST },
        { "fc_measurement", FC_STRUCT_MEASUREMENT }, { "fc_solution_summary", FC_STRUCT_SUMMARY },
        { "fc_gr_active_stats", FC_STRUCT_GR_ACTIVE } };
    for (const auto& s : headered)
        std::printf ("T %s %d %u\n", s.name, s.id, fc_master_sizeof (s.id, fc_master_abi_version()));
    return 0;
}

//==============================================================================
// THE TWO FILE-WRITING COMMANDS. Functions rather than `main`'s body so the selftest can run them and read back
// what they wrote — the file's length is part of the contract, and a length nobody reads back is a claim.
//
// THE DELIVERED LENGTH IS THE FILE'S LENGTH. Every buffer, the de-planarisation and the write take `outFrames`
// on a delivering handle — a file written with the input's count would be SHORT on an upsample, and its second
// channel would be read out of the middle of the first plane.

int cmdRender (const Args& a, const std::vector<float>& in, std::size_t frames, int nc, const char* path)
{
    std::vector<float> out; fc_master_resolved res {};
    std::size_t outFrames = frames;
    if (a.cfg.deliveryRate != 0.0 ? ! abiRenderDelivered (a, in, frames, nc, out, outFrames, res)
                                  : ! abiRender (a, in, frames, nc, out, res)) return 1;
    if (! writeInterleaved (path, nc, out, outFrames)) return 1;
    // The resolved geometry on stderr, so a harness can read the numbers without parsing the audio.
    std::fprintf (stderr, "latency=%d internalBlock=%d ceiling=%.17g release=%.17g\n",
                  res.latencySamples, res.internalBlock, res.limiterCeilingDbTp, res.limiterReleaseMs);
    return 0;
}

int cmdSolve (const Args& a, const std::vector<float>& in, std::size_t frames, int nc, const char* path)
{
    const bool delivering = a.cfg.deliveryRate != 0.0;
    fc_master h = 0;
    if (const fc_status st = fc_master_create (&a.cfg, &h); st != FC_OK)
    { std::fprintf (stderr, "create: %s\n", statusName (st)); return 2; }
    for (const auto& w : a.weights)
        if (fc_master_set_channel_weight (h, w.channel, w.value) != FC_OK)
        { std::fprintf (stderr, "weight%d rejected\n", w.channel); fc_master_destroy (h); return 2; }
    std::uint32_t outFrames = (std::uint32_t) frames;
    if (delivering)
        if (const fc_status st = fc_master_delivered_frames (h, (std::uint32_t) frames, &outFrames); st != FC_OK)
        { std::fprintf (stderr, "delivered_frames: %s\n", statusName (st)); fc_master_destroy (h); return 1; }
    std::vector<float> out ((std::size_t) outFrames * (std::size_t) nc, 0.0f);
    fc_solution sol = 0;
    const fc_status st = delivering
        ? fc_master_solve_delivered (h, &a.prm, &a.req, in.data(), (std::uint32_t) frames,
                                     out.data(), outFrames, &sol)
        : fc_master_solve (h, &a.prm, &a.req, in.data(), out.data(), (std::uint32_t) frames, &sol);
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
    // v4 — each stage's trace, summarised for a harness: its bucket count, whether it is a measurement, the largest
    // bucket maximum and the first bucket that holds it. The buckets themselves are `fc_solution_gr_trace`.
    for (const auto& [label, code] : { std::pair<const char*, int> { "comp", FC_GR_STAGE_COMPRESSOR },
                                       std::pair<const char*, int> { "lim",  FC_GR_STAGE_LIMITER } })
    {
        const int buckets = code == FC_GR_STAGE_LIMITER ? meas.limiterGrTraceBuckets : meas.compressorGrTraceBuckets;
        std::vector<fc_gr_trace_bucket64> tb ((std::size_t) std::max (buckets, 0));
        std::uint32_t w = 0;
        (void) fc_solution_gr_trace64 (sol, code, tb.data(), (std::uint32_t) tb.size(), &w);
        double mx = 0.0; std::uint32_t at = 0;
        for (std::uint32_t i = 0; i < w; ++i) if (tb[i].maxDb > mx) { mx = tb[i].maxDb; at = i; }
        const int valid = code == FC_GR_STAGE_LIMITER ? meas.limiterGrTraceValid : meas.compressorGrTraceValid;
        std::printf ("%sTrace buckets=%u valid=%d max=%.17g at=%u\n", label, w, valid, mx, at);
    }
    // A VERDICT IS NOT A RENDER. `InvalidRequest` and `NotPrepared` are returned before the solver
    // touches the output, so `out` is still the zero buffer it was allocated as — writing it would
    // hand a harness a file of the right length, full of digital silence, with exit status 0 and an
    // existing output overwritten. The status line has already been printed; the file is not.
    const bool delivered = (sum.status != FC_SOLVE_INVALID_REQUEST
                         && sum.status != FC_SOLVE_NOT_PREPARED
                         && sum.status != FC_SOLVE_RENDER_FAILED);
    bool ok = true;
    if (delivered) ok = writeInterleaved (path, nc, out, outFrames);
    else std::fprintf (stderr, "no render was delivered (status=%d) — the output file is NOT written\n",
                       sum.status);
    fc_solution_destroy (sol);
    fc_master_destroy (h);
    return (delivered && ok) ? 0 : 1;
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
// periodic transients so the limiter and the clipper are not asleep, DIFFERENT CONTENT PER CHANNEL so a
// planar-addressing fault is visible, and a stretch of EXACT DIGITAL SILENCE at the end so the dither's
// auto-blank has something to blank — without which `autoBlankSamples` maps for free. A fixture on which
// a stage is inert makes the whole comparison pass for the wrong reason, one field at a time.
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
            // The last eighth is exact digital silence — the auto-blank's own trigger.
            if (i >= frames - frames / 8) x = 0.0;
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
    check (initArgs (a), "the versioned defaults writers accept a struct stamped at this build's version");
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
    // ASYM with a bias, not Atan: the DC blocker only has anything to remove when the shaper is
    // asymmetric, so `dcBlockHz` is unmappable-for-free on any symmetric shape.
    a.prm.clipper.shape         = FC_SHAPE_ASYM;
    a.prm.dither.shaping        = FC_SHAPING_PSYCHO;
    // `autoBlank` STAYS ON, and this is the second time the same mistake was caught in this fixture:
    // moving a field off its default switched off the mechanism the NEXT field controls, so
    // `autoBlankSamples` became unmappable for free — exactly as setting `bypassMonoBass` here had
    // earlier disabled the stage the lifecycle-order check depends on. Measured: with auto-blank off,
    // 3777 against the default 4096 renders identically; with it on and a silent tail, they differ in
    // 412 bytes. The OFF value is covered by the dither-only comparison instead, where it blinds nothing.
    a.prm.dither.autoBlank      = 1;
    a.prm.dither.autoBlankSamples = 3777;
    a.prm.compressor.autoMakeup = 1;
    a.prm.eqBands[2].on         = 1;
    a.prm.eqBands[2].type       = FC_FILTER_HIGH_SHELF;
    // `swept` is inert unless the band is a CUT/notch/band-pass with ONLY the Stereo lane enabled
    // (`EqBand.h:156`), so band 2 — a high shelf on lanes 3 and 4 — cannot show it and the flag mapped
    // for free. It goes on band 0, the high-pass, which is exactly that configuration.
    a.prm.eqBands[0].swept      = 1;
    a.prm.eqBands[2].lanes[3].on = 1;                  // a lane that is NOT lane 0
    a.prm.eqBands[2].lanes[3].freq = 7331.7;
    a.prm.eqBands[2].lanes[3].gainDb = 1.7;
    a.prm.eqBands[2].lanes[3].q = 0.77;
    // A BYPASSED LANE MUST BE DISTINGUISHABLE FROM AN ABSENT ONE, or `bypass` is unmapped for free: a
    // lane at its default 0 dB does nothing whether it is bypassed or not, so the flag needs a lane
    // that WOULD be audible.
    a.prm.eqBands[2].lanes[4].on = 1;
    a.prm.eqBands[2].lanes[4].freq = 313.7;
    a.prm.eqBands[2].lanes[4].gainDb = 5.3;
    a.prm.eqBands[2].lanes[4].q = 1.3;
    a.prm.eqBands[2].lanes[4].bypass = 1;              // `bypass` is not `on`, and both are mapped
    a.prm.eqBands[2].dyn.on     = 1;
    // Engaged, not merely enabled: an absolute threshold well inside the programme's own level, and a
    // range big enough to hear. Otherwise `thrAuto` (absolute vs relative) changes nothing measurable.
    a.prm.eqBands[2].dyn.rangeDb = -9.3;
    a.prm.eqBands[2].dyn.thrDb  = -20.3;
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
    // v3. Off its default like every other field, and not a binary32-exact value: at 1 a dropped mapping renders
    // the chain before the field existed and the compare stays green.
    a.prm.compressorMix          = 0.73;
    a.prm.limiter.ceilingDbTp    = -1.3;
    a.prm.limiter.releaseMs      = 77.0;
    // v6
    a.prm.limiterDualRelease     = 1;
    a.prm.limiterSlowReleaseMs   = 173.7;
    a.prm.peakClipper              = 1;
    a.prm.peakClipperOverCeilingDb = 1.7;
    a.prm.peakClipperKneeDb        = 0.4;
    a.req.grTraceBuckets         = 4099;
    // v8, off their defaults like every other field, and the statistics moved with them so the numbers are read
    // rather than merely carried.
    a.req.limiterGr.statistic    = FC_GR_PERCENTILE;
    a.req.compressorGr.statistic = FC_GR_PERCENTILE;
    a.req.limiterGrQuantile      = 0.877;
    a.req.compressorGrQuantile   = 0.611;
    // v10 — EVERY NEW FIELD MOVED OFF ITS DEFAULT, which is what this fixture is for: a field left at its
    // default crosses the ABI identically whether it is marshalled or forgotten.
    a.req.limiterActiveInputDb   = -47.5;
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

    // --- 1d. A MINIMAL TOPOLOGY -----------------------------------------------------------------
    // Dither alone. Two reasons, both of them measured rather than tidy: a field is only pinned in a
    // topology where it can be SEEN, and `dither.autoBlankSamples` cannot be seen through the full
    // chain — the limiter's tail is not exact digital zero, so the auto-blank never fires and 3777
    // against the default 4096 renders identically. Through a dither-only chain the same pair differs
    // in 433 bytes. It also exercises the corner where almost every stage is absent, which the full
    // fixture cannot.
    {
        // TWICE, with the auto-blank ON and OFF, because one run cannot cover both: with the blank OFF
        // its sample count is inert, and with it ON the count is what the fixture sees. Two comparisons
        // over one topology is what it takes for a flag AND the field it gates to be pinned at once.
        std::vector<float> blanked;
        for (int autoBlank = 1; autoBlank >= 0; --autoBlank)
        {
            Args b = a;
            b.cfg.eq = b.cfg.monoBass = b.cfg.compressor = b.cfg.clipper = b.cfg.limiter = 0;
            b.cfg.dither = 1;
            b.prm.dither.autoBlank = autoBlank;
            std::vector<float> mAbi, mCpp; fc_master_resolved rm {};
            const bool ranM = abiRender (b, in, frames, nc, mAbi, rm)
                           && directRender (b, in, frames, nc, mCpp);
            check (ranM, autoBlank ? "the dither-only render ran (auto-blank ON)"
                                   : "the dither-only render ran (auto-blank OFF)");
            if (! ranM) continue;
            double worst = 0.0;
            const std::size_t d = bitDiff (mAbi, mCpp, worst);
            char msg[160];
            std::snprintf (msg, sizeof msg, "auto-blank %s: %zu differ, worst %.9g",
                           autoBlank ? "on" : "off", d, worst);
            check (d == 0, "a topology with almost every stage ABSENT maps identically", msg);
            if (autoBlank) blanked = mAbi;
            else
            {
                // PRECONDITION for the pair: the flag has to CHANGE the render, or both runs are the
                // same experiment and neither pins anything about it.
                double delta = 0.0;
                for (std::size_t i = 0; i < mAbi.size() && i < blanked.size(); ++i)
                    delta = std::max (delta, (double) std::fabs (mAbi[i] - blanked[i]));
                char pre[128];
                std::snprintf (pre, sizeof pre, "on against off differ by up to %.3g", delta);
                check (delta > 0.0, "PRECONDITION: the auto-blank changes this render", pre);
            }
            char lat[96];
            std::snprintf (lat, sizeof lat, "latency %d against the full chain's %d",
                           rm.latencySamples, res.latencySamples);
            check (rm.latencySamples < res.latencySamples,
                   "PRECONDITION: the absent stages really are absent — the latency dropped", lat);
        }
    }

    // --- 1e. THE DUAL RELEASE (v6), on a 1 kHz tone 400 ms loud and 200 ms quiet: mapped identically, changing the
    // render, and not read from a v5-stamped parameter set.
    {
        std::vector<float> held (frames * (std::size_t) nc, 0.0f);
        for (int c = 0; c < nc; ++c)
            for (std::size_t i = 0; i < frames; ++i)
            {
                const double t = (double) i / fs;
                const double amp = std::fmod (t, 0.6) < 0.4 ? 0.9 : 0.3;
                held[(std::size_t) c * frames + i] = (float) (amp * std::sin (2.0 * 3.14159265358979 * (1000.0 + 7.0 * c) * t));
            }
        Args b = a;
        b.prm.preLimiterGainDb = 12.0;
        auto render = [&] (const Args& x, std::vector<float>& viaAbiX, std::vector<float>& viaCppX)
        {
            fc_master_resolved rx {};
            return abiRender (x, held, frames, nc, viaAbiX, rx) && directRender (x, held, frames, nc, viaCppX);
        };
        std::vector<float> dAbi, dCpp, offAbi, offCpp, slowAbi, slowCpp, v5Abi, v5Cpp;
        // EVERY field past the stamped size, not just v6's: a v5 caller cannot reach v11's clipper either,
        // so the baseline it is compared against must hold the clipper at ITS defaults too. Written out
        // rather than left at the fixture's values — the fixture deliberately turns the clipper ON.
        Args off = b;  off.prm.limiterDualRelease = 0;
        off.prm.peakClipper = 0;
        off.prm.peakClipperOverCeilingDb = felitronics::limiter::TruePeakLimiterParams {}.overCeilingDb;
        off.prm.peakClipperKneeDb        = felitronics::limiter::TruePeakLimiterParams {}.kneeDb;
        Args slow = b; slow.prm.limiterSlowReleaseMs = 200.0;
        Args v5 = b;   v5.prm.header.abiVersion = 5u; v5.prm.header.structSize = 6568u;
        const bool ran = render (b, dAbi, dCpp) && render (off, offAbi, offCpp) && render (slow, slowAbi, slowCpp)
                      && render (v5, v5Abi, v5Cpp);
        check (ran, "the dual-release renders ran through both paths");
        if (ran)
        {
            double worst = 0.0;
            char msg[160];
            std::size_t d = bitDiff (dAbi, dCpp, worst) + bitDiff (offAbi, offCpp, worst) + bitDiff (slowAbi, slowCpp, worst);
            std::snprintf (msg, sizeof msg, "%zu differ", d);
            check (d == 0, "the dual release maps identically through the ABI — on, off, and at another slow release", msg);
            const std::size_t onOff = bitDiff (dAbi, offAbi, worst), slowMoved = bitDiff (dAbi, slowAbi, worst);
            std::snprintf (msg, sizeof msg, "on against off: %zu differ; 173.7 against 200 ms: %zu differ", onOff, slowMoved);
            check (onOff > 0 && slowMoved > 0, "PRECONDITION: the flag and the slow release both change this render", msg);
            d = bitDiff (v5Abi, offAbi, worst);
            std::snprintf (msg, sizeof msg, "%zu differ from the dual release off", d);
            check (d == 0, "a v5-stamped parameter set does not read the two fields past its 6568 bytes", msg);
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


    // --- 5. THE CLI'S OWN ARGUMENT PARSER --------------------------------------------------------
    // It lives in this binary and had no test at all, so its strictness was a claim in a comment. A
    // mutation stand removed the full-parse check and nothing went red; another removed the channel
    // check that stands in front of a division and nothing went red either.
    {
        Args t;
        (void) initArgs (t);
        check (applyKey (t, "lim.ceiling", "-1.5"), "a good numeric key is accepted");
        check (applyKey (t, "delivery", "44100") && t.cfg.deliveryRate == 44100.0, "the delivery rate is a key");
        check (! applyKey (t, "delivery", "44.1k"), "and a malformed one is refused, not read as 0 — which is no conversion");
        check (applyKey (t, "comp.mix", "0.25") && t.prm.compressorMix == 0.25 && t.prm.clipper.mix != 0.25f,
               "comp.mix is the CHAIN's compressor mix, not the clipper's");
        check (! applyKey (t, "comp.mix", "x"), "and a malformed mix is refused");
        check (! applyKey (t, "lim.ceiling", "oops"), "a non-numeric value is REFUSED, not read as 0");
        check (! applyKey (t, "comp.ratio", "4oops"), "and so is a numeric prefix with a tail");
        check (! applyKey (t, "comp.ratio", ""), "and an empty value");
        check (! applyKey (t, "block", "4294967296"),
               "a block that cannot fit a positive int is refused — it used to narrow to 0 and HANG");
        check (! applyKey (t, "block", "0"), "and so is a zero block");
        check (applyKey (t, "block", "128"), "while a real one is accepted");
        check (! applyKey (t, "band4294967296.gain", "6"),
               "a band index past the array is refused — it used to narrow onto band 0");
        check (! applyKey (t, "bandXYZ.gain", "6"), "and so is a non-numeric one");
        check (applyKey (t, "band3.gain", "6"), "while a real band is accepted");
        check (! applyKey (t, "nosuchkey", "1"), "an unknown key is an error, not a shrug");
        check (! applyKey (t, "comp.detector", "sideways"), "and so is an unknown enum name");
        check (applyKey (t, "comp.detector", "rms"), "while a real one is accepted");
        double d = 0.0; long i = 0;
        check (positional ("48000", d) && d == 48000.0, "a positional rate parses");
        check (! positional ("48000oops", d), "and one with a tail does not");
        check (! positional ("2x", i), "nor does a channel count with a tail");
        check (positional ("2", i) && i == 2, "while a real one does");
    }

    // --- 6. ABI v2: A v1 CALLER IS A v1 CALLER ----------------------------------------------------------
    // `a` is a v2 config with `deliveryRate = 0` written explicitly, and section 1 has already shown it renders
    // what the direct C++ path renders. The same struct stamped v1 — and carrying 96 000 in the bytes past its
    // 80, where a facade reading this build's `sizeof` would find a delivery rate — must render the same bits and
    // must make a handle that does not deliver. (Against the v1 BUILD itself is the third corner of this check,
    // and it is run outside this binary: an in-binary comparison cannot see a change both paths share.)
    {
        // The oracle for an older caller is the CURRENT caller with every newer field at its previous-version
        // value, written explicitly — and the direct C++ path at those values, so the pair is not two new roads
        // agreeing with each other.
        Args prev = a;
        prev.cfg.deliveryRate = 0.0;
        prev.prm.compressorMix = 1.0;
        prev.prm.limiterDualRelease = 0;
        prev.prm.peakClipper = 0;                                              // v11
        prev.prm.peakClipperOverCeilingDb = felitronics::limiter::TruePeakLimiterParams {}.overCeilingDb;
        prev.prm.peakClipperKneeDb        = felitronics::limiter::TruePeakLimiterParams {}.kneeDb;
        std::vector<float> viaPrev, viaPrevCpp; fc_master_resolved rp {};
        const bool ranPrev = abiRender (prev, in, frames, nc, viaPrev, rp) && directRender (prev, in, frames, nc, viaPrevCpp);
        double wp = 0.0;
        const std::size_t dp = ranPrev ? bitDiff (viaPrev, viaPrevCpp, wp) : 1u;
        check (ranPrev && dp == 0, "the current caller at the previous versions' values renders the core's bits");
        double moved = 0.0;
        for (std::size_t i = 0; ranPrev && i < viaPrev.size(); ++i) moved = std::max (moved, (double) std::fabs (viaPrev[i] - viaAbi[i]));
        check (moved > 1e-3, "PRECONDITION: and compressorMix 0.73 really changed the main render against mix 1");

        for (const std::uint32_t ver : { 1u, 2u })
        {
            Args old = a;                                         // the newer fields hold NON-default values ...
            old.cfg.header.abiVersion = ver;  old.cfg.header.structSize = ver == 1u ? 80u : 88u;
            old.prm.header.abiVersion = ver;  old.prm.header.structSize = 6560u;
            if (ver == 1u) old.cfg.deliveryRate = 96000.0;        // ... past the stamped size: must never be read
            old.prm.compressorMix = 0.2;
            std::vector<float> viaOld; fc_master_resolved ro {};
            const bool ran = abiRender (old, in, frames, nc, viaOld, ro);
            char what[96];
            std::snprintf (what, sizeof what, "a v%u-stamped config and parameter set are accepted by the v%u build",
                           ver, (unsigned) FC_MASTER_ABI_VERSION);
            check (ran, what);
            if (! ran) continue;
            double worst = 0.0;
            const std::size_t d = bitDiff (viaPrev, viaOld, worst);
            char msg[128];
            std::snprintf (msg, sizeof msg, "v%u: %zu of %zu samples differ, worst %.9g", ver, d, viaOld.size(), worst);
            check (d == 0, "an older caller renders the current caller's bits at the previous values — nothing past its size was read", msg);
        }
        Args v1 = a;
        v1.cfg.header.abiVersion = 1u;  v1.cfg.header.structSize = 80u;
        v1.cfg.deliveryRate = 96000.0;
        fc_master h = 0; std::uint32_t df = 0;
        const bool made = fc_master_create (&v1.cfg, &h) == FC_OK;
        check (made && fc_master_delivered_frames (h, 1000u, &df) == FC_ERR_STATE,
               "and the delivery rate in the bytes past v1's 80 was not read: the handle does not deliver");
        if (made) fc_master_destroy (h);
    }

    // --- 7. THE DELIVERING HANDLE, through the ABI and through the core ---------------------------------
    // Down, up, and EQUAL rates. The direct path calls `mastering::DeliveredMastering` — the class the facade
    // forwards to — so a difference is marshalling. At equal rates there is a stronger oracle than that: the
    // converter copies bits, so the delivered render must be the plain render of section 1, bit for bit.
    // DOWN is the highest delivery rate below `fs` and UP the lowest above; at 44.1 and 192 kHz one of them does not
    // exist, and the output says so rather than running one direction twice.
    const double kDeliveryRates[] { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };
    double down = 0.0, up = 0.0;
    for (double r : kDeliveryRates) { if (r < fs) down = r; if (r > fs && up == 0.0) up = r; }
    if (down == 0.0 || up == 0.0)
        std::printf ("  (note: %g Hz has no %s delivery rate — that direction is not exercised)\n", fs, down == 0.0 ? "lower" : "higher");
    std::vector<double> pairs;
    for (double r : { down, up, fs }) if (r != 0.0) pairs.push_back (r);
    // A temp name no other run shares: two selftests at once (two build trees, one machine) must not read each
    // other's files. The address of a local in this process and the clock make it unique enough for that.
    const std::string runTag = std::to_string ((unsigned long long) (std::uintptr_t) &pairs) + "_"
                             + std::to_string ((long long) std::chrono::steady_clock::now().time_since_epoch().count());
    for (double dr : pairs)
    {
        Args b = a;
        b.cfg.deliveryRate = dr;
        std::vector<float> dAbi, dCpp; std::size_t nAbi = 0, nCpp = 0; fc_master_resolved rd {};
        const bool ran = abiRenderDelivered (b, in, frames, nc, dAbi, nAbi, rd)
                      && directRenderDelivered (b, in, frames, nc, dCpp, nCpp);
        char what[96];
        std::snprintf (what, sizeof what, "%g -> %g: the delivered render ran through both paths", fs, dr);
        check (ran, what);
        if (! ran) continue;
        const long long expect = DeliveredMastering::deliveredFrames (fs, dr, (long long) frames);
        char len[128];
        std::snprintf (len, sizeof len, "ABI %zu, direct %zu, the converter's formula %lld", nAbi, nCpp, expect);
        check ((long long) nAbi == expect && nCpp == nAbi, "the delivered length is the core's, on both paths", len);
        double worst = 0.0;
        const std::size_t d = bitDiff (dAbi, dCpp, worst);
        char msg[128];
        std::snprintf (msg, sizeof msg, "%g -> %g: %zu of %zu differ, worst %.9g", fs, dr, d, dAbi.size(), worst);
        check (d == 0, "the delivered render through the ABI is BIT-IDENTICAL to the core's", msg);
        double amp = 0.0;
        for (float v : dAbi) amp = std::max (amp, (double) std::fabs (v));
        char pre[96];
        std::snprintf (pre, sizeof pre, "peaks at %.4f", amp);
        check (amp > 0.05, "PRECONDITION: the delivered render is not silence", pre);
        if (dr == fs)
        {
            const std::size_t e = bitDiff (dAbi, viaAbi, worst);
            std::snprintf (msg, sizeof msg, "%zu of %zu differ, worst %.9g", e, dAbi.size(), worst);
            check (e == 0, "EQUAL RATES: the delivered render is the plain render, bit for bit", msg);
        }

        // The file a command writes is `outFrames` long and its channels are the delivered planes — read back
        // from disk, not inferred from the buffer that was handed to the writer.
        std::error_code ec;
        const auto dir = std::filesystem::temp_directory_path (ec);
        if (ec) { check (false, "a temp directory for the file check"); continue; }
        const std::string path = (dir / ("fcore_master_selftest_" + runTag + "_" + std::to_string ((long long) dr) + ".f32")).string();
        const int rc = cmdRender (b, in, frames, nc, path.c_str());
        std::vector<float> back; std::size_t backFrames = 0;
        const bool read = rc == 0 && readInterleaved (path.c_str(), nc, back, backFrames);
        std::remove (path.c_str());
        char fl[128];
        std::snprintf (fl, sizeof fl, "rc %d, %zu frames on disk against %zu delivered", rc, backFrames, nAbi);
        check (read && backFrames == nAbi, "the render command writes the DELIVERED length", fl);
        if (read && backFrames == nAbi)
        {
            const std::size_t f = bitDiff (back, dAbi, worst);
            std::snprintf (fl, sizeof fl, "%zu of %zu differ", f, back.size());
            check (f == 0, "and each channel of the file is its own delivered plane", fl);
        }
    }

    // --- 8. THE DELIVERED SEARCH AND THE DELIVERED RANGE ------------------------------------------------
    {
        Args b = a;
        b.cfg.deliveryRate = pairs.front() == fs ? 48000.0 : pairs.front();
        b.req.targetLufs = -14.0; b.req.maxTruePeakDbTp = -1.0;
        MasteringChainConfig cc {}; MasteringChainParams cp {};
        mirror (b, cc, cp);
        LoudnessRequest lr {};
        lr.targetLufs = -14.0; lr.maxTruePeakDbTp = -1.0;
        lr.grTraceBuckets = b.req.grTraceBuckets;
        lr.limiterGr.statistic    = GrStatistic::Percentile;    // v8
        lr.compressorGr.statistic = GrStatistic::Percentile;
        lr.limiterGr.quantile     = b.req.limiterGrQuantile;
        lr.compressorGr.quantile  = b.req.compressorGrQuantile;

        fc_master h = 0;
        const bool made = fc_master_create (&b.cfg, &h) == FC_OK;
        check (made, "a delivering handle for the search");
        if (made)
        {
            std::uint32_t d = 0;
            (void) fc_master_delivered_frames (h, (std::uint32_t) frames, &d);
            double lraAbi = -1.0;
            const fc_status stL = fc_master_measure_lra (h, in.data(), (std::uint32_t) frames, &lraAbi);
            std::vector<float> sAbi ((std::size_t) d * (std::size_t) nc, 0.0f);
            fc_solution sol = 0;
            const fc_status stS = fc_master_solve_delivered (h, &b.prm, &b.req, in.data(), (std::uint32_t) frames,
                                                             sAbi.data(), d, &sol);
            fc_solution_summary sum {}; FC_INIT (sum);
            const bool solved = stS == FC_OK && fc_solution_summary_get (sol, &sum) == FC_OK;

            // The direct path, with the facade's own geometry: its renderer block and the solver it prepares.
            MasteringChain chain; OfflineRenderer r; TargetLoudnessSolver solver; DeliveredMastering dm;
            const bool built = r.prepare (nc, 4096) && chain.prepare (b.cfg.deliveryRate, nc, cc)
                            && dm.prepare (fs, b.cfg.deliveryRate, nc, 4096)
                            && solver.prepare (b.cfg.deliveryRate, nc, r.blockSize(), chain.internalBlock(),
                                               chain.tapOversampleFactor());
            check (built, "the direct delivered search is built");
            double lraCpp = -2.0;
            const float* ip[core::kMaxChannels] {};
            for (int c = 0; c < nc; ++c) ip[c] = in.data() + (std::size_t) c * frames;
            const bool lraOk = built && dm.measureInputLoudnessRange (solver, ip, nc, (long long) frames, lraCpp);
            char lm[128];
            std::snprintf (lm, sizeof lm, "ABI %s %.17g, core %.17g", statusName (stL), lraAbi, lraCpp);
            check (stL == FC_OK && lraOk && lraAbi == lraCpp, "the delivered range through the ABI is the core's", lm);

            std::vector<float> sCpp ((std::size_t) d * (std::size_t) nc, 0.0f);
            float* op[core::kMaxChannels] {};
            for (int c = 0; c < nc; ++c) op[c] = sCpp.data() + (std::size_t) c * (std::size_t) d;
            LoudnessSolution direct;
            if (built) direct = dm.solve (solver, chain, r, cp, ip, nc, (long long) frames, op, (long long) d, lr);
            char sm[160];
            std::snprintf (sm, sizeof sm, "ABI %s status %d gain %.17g, core status %d gain %.17g", statusName (stS),
                           sum.status, sum.preLimiterGainDb, (int) direct.status, direct.preLimiterGainDb);
            check (solved && sum.status == (int) direct.status && sum.preLimiterGainDb == direct.preLimiterGainDb
                   && sum.ceilingDbTp == direct.ceilingDbTp, "the delivered search's verdict is the core's", sm);
            double worst = 0.0;
            const std::size_t dd = bitDiff (sAbi, sCpp, worst);
            std::snprintf (sm, sizeof sm, "%zu of %zu differ, worst %.9g", dd, sAbi.size(), worst);
            check (solved && dd == 0, "and its delivered audio is bit-identical", sm);

            // v4 — the gain-reduction traces through the ABI are the core's, bit for bit, both stages, and the
            // measurement's four trace fields say what the core's traces say.
            {
                fc_measurement meas {}; FC_INIT (meas);
                const bool mOk = solved && fc_solution_measurement (sol, &meas) == FC_OK;
                check (mOk && meas.compressorGrTraceBuckets == direct.compressorTrace.buckets
                       && meas.limiterGrTraceBuckets == direct.limiterTrace.buckets
                       && meas.compressorGrTraceValid == (direct.compressorTrace.valid ? 1 : 0)
                       && meas.limiterGrTraceValid == (direct.limiterTrace.valid ? 1 : 0)
                       && direct.limiterTrace.buckets == b.req.grTraceBuckets,
                       "the trace's bucket counts (the request's, v6) and validity through the ABI are the core's");
                std::size_t traceDiff = 0;
                for (const auto& [code, t] : { std::pair<int, const GainReductionTrace*> { FC_GR_STAGE_COMPRESSOR, &direct.compressorTrace },
                                               std::pair<int, const GainReductionTrace*> { FC_GR_STAGE_LIMITER,    &direct.limiterTrace } })
                {
                    std::vector<fc_gr_trace_bucket64> tb ((std::size_t) t->buckets);
                    std::vector<fc_gr_trace_bucket>   t4 ((std::size_t) t->buckets);
                    std::uint32_t w = 0, w4 = 0;
                    if (! solved || fc_solution_gr_trace64 (sol, code, tb.data(), (std::uint32_t) tb.size(), &w) != FC_OK
                        || fc_solution_gr_trace (sol, code, t4.data(), (std::uint32_t) t4.size(), &w4) != FC_OK
                        || (int) w != t->buckets || w4 != w) { ++traceDiff; continue; }
                    for (std::uint32_t i = 0; i < w; ++i)
                        if (std::memcmp (&tb[i].maxDb, &t->bucket[i].maxDb, 8) != 0 || std::memcmp (&tb[i].meanDb, &t->bucket[i].meanDb, 8) != 0
                            || tb[i].samples != t->bucket[i].samples || tb[i].nonFinite != t->bucket[i].nonFinite
                            || std::memcmp (&t4[i].maxDb, &t->bucket[i].maxDb, 8) != 0 || std::memcmp (&t4[i].meanDb, &t->bucket[i].meanDb, 8) != 0
                            || (std::uint64_t) t4[i].samples != t->bucket[i].samples || (std::uint64_t) t4[i].nonFinite != t->bucket[i].nonFinite) ++traceDiff;
                }
                std::snprintf (sm, sizeof sm, "%zu buckets differ", traceDiff);
                check (traceDiff == 0, "and both traces through the ABI are the core's, bit for bit", sm);
            }

            // v8 — the quantile read-back is the core's own, at every fraction and on both stages, and it refuses
            // where the core refuses. `fc_solution_gr_quantile` states no definition of its own, so this is the
            // whole of the acceptance for it.
            {
                std::size_t qDiff = 0;
                for (const auto& [code, st] : { std::pair<int, GrStage> { FC_GR_STAGE_COMPRESSOR, GrStage::Compressor },
                                                std::pair<int, GrStage> { FC_GR_STAGE_LIMITER,    GrStage::Limiter } })
                    for (const double q : { 0.05, 0.5, 0.611, 0.877, 0.95, 1.0 })
                    {
                        double want = 0.0, got = 0.0;
                        const bool coreOk = direct.grQuantile (st, q, want);
                        const fc_status qs = solved ? fc_solution_gr_quantile (sol, code, q, &got)
                                                    : FC_ERR_HANDLE;
                        if (coreOk != (qs == FC_OK) || (coreOk && std::memcmp (&want, &got, 8) != 0)) ++qDiff;
                    }
                std::snprintf (sm, sizeof sm, "%zu of 12 readings differ", qDiff);
                check (solved && qDiff == 0, "the quantile read-back through the ABI is the core's, bit for bit", sm);
                double sink = 12345.0;
                const fc_status qNan = solved ? fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER,
                                                                         std::numeric_limits<double>::quiet_NaN(), &sink)
                                              : FC_ERR_HANDLE;
                const fc_status qHi  = solved ? fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 1.5, &sink) : FC_ERR_HANDLE;
                const fc_status qLo  = solved ? fc_solution_gr_quantile (sol, FC_GR_STAGE_LIMITER, 0.0, &sink) : FC_ERR_HANDLE;
                const fc_status qEnum = solved ? fc_solution_gr_quantile (sol, 7, 0.5, &sink) : FC_ERR_HANDLE;
                check (qNan == FC_ERR_NON_FINITE && qHi == FC_ERR_RANGE && qLo == FC_ERR_RANGE
                       && qEnum == FC_ERR_ENUM && sink == 12345.0,
                       "and its refusals are named apart and write nothing");
            }
            if (solved) fc_solution_destroy (sol);
            fc_master_destroy (h);

            // The SOLVE command's file, read back like the render's: the delivered length, and each plane its own.
            std::error_code ec;
            const auto dir = std::filesystem::temp_directory_path (ec);
            const std::string path = (dir / ("fcore_master_selftest_" + runTag + "_solve.f32")).string();
            const int rc = ec ? -1 : cmdSolve (b, in, frames, nc, path.c_str());
            std::vector<float> back; std::size_t backFrames = 0;
            const bool read = rc == 0 && readInterleaved (path.c_str(), nc, back, backFrames);
            std::remove (path.c_str());
            char fl[128];
            std::snprintf (fl, sizeof fl, "rc %d, %zu frames on disk against %u delivered", rc, backFrames, d);
            check (read && backFrames == d && bitDiff (back, sAbi, worst) == 0,
                   "the solve command writes the DELIVERED length, and the delivered master in it", fl);
        }
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
            "  %s selftest [sampleRate] [channels]\n"
            "  %s layout\n"
            "  (delivery=<rate> on render/solve/lra: SRC first, the output at <rate>)\n",
            argv[0], argv[0], argv[0], argv[0], argv[0]);
        return 2;
    }
    const std::string mode = argv[1];

    if (mode == "layout") return printLayout();

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
        double sr = 0.0; long nc = 0;
        if (! positional (argv[2], sr) || ! positional (argv[3], nc))
        { std::fprintf (stderr, "sampleRate and channels must be numbers\n"); return 2; }
        a.cfg.sampleRate = sr;
        a.cfg.channels   = (std::int32_t) nc;
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
        for (const auto& w : a.weights)
            if (fc_master_set_channel_weight (h, w.channel, w.value) != FC_OK)
            { std::fprintf (stderr, "weight%d rejected\n", w.channel); fc_master_destroy (h); return 2; }
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
    double sr = 0.0; long ncL = 0;
    if (! positional (argv[2], sr) || ! positional (argv[3], ncL))
    { std::fprintf (stderr, "sampleRate and channels must be numbers\n"); return 2; }
    a.cfg.sampleRate = sr;
    if (ncL < 1 || ncL > core::kMaxChannels) { std::fprintf (stderr, "bad channel count\n"); return 2; }
    a.cfg.channels   = (std::int32_t) ncL;
    const int nc = a.cfg.channels;

    std::vector<float> in; std::size_t frames = 0;
    if (! readInterleaved (argv[4], nc, in, frames)) return 2;
    if (frames == 0) { std::fprintf (stderr, "empty input\n"); return 2; }
    // The ABI counts frames in 32 bits. Narrowing silently would not merely lose length: the planar
    // STRIDE is the frame count, so channel 1's pointer would land inside channel 0's plane and the
    // render would be of a signal that does not exist.
    if (frames > (std::size_t) 0x7FFFFFFFu)
    { std::fprintf (stderr, "input longer than the ABI's frame count (%zu frames)\n", frames); return 2; }

    if (mode == "render") return cmdRender (a, in, frames, nc, argv[5]);
    if (mode == "solve")  return cmdSolve (a, in, frames, nc, argv[5]);

    std::fprintf (stderr, "unknown mode '%s'\n", mode.c_str());
    return 2;
}
