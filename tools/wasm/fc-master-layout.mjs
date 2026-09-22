// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The JavaScript side of tools/fc_master_abi.h: the struct layouts, and a reader/writer over the wasm
// heap that uses them. It is imported by the node parity harness here and COPIED VERBATIM into the site
// that ships the module — one description of the layout, not two, for the reason the ABI header gives
// about its own mappings: a second copy is where a mapping starts drifting.
//
// ==================================================================================================
// WHY THIS COMPUTES OFFSETS INSTEAD OF LISTING THEM
// ==================================================================================================
// A list of numbers is a transcription, and a transcription of 200 offsets is wrong somewhere. What is
// written below is the thing that does not change — each struct's FIELDS, IN ORDER, WITH THEIR TYPES,
// copied from the header one line at a time — and the offsets fall out of the C layout algorithm, which
// for these five scalar types is four lines: align each field to its own alignment, the struct's
// alignment is the widest member's, the size is rounded up to it. (wasm32 aligns f64/u64 to 8, the same
// as the x86-64 and arm64 hosts the native reference runs on, so one computation serves every tier.)
//
// ==================================================================================================
// AND WHY A WRONG LAYOUT CANNOT RENDER ANYTHING
// ==================================================================================================
// `structSize` is not decoration. Every struct crossing this boundary starts with {abiVersion,
// structSize}, this file stamps them from the size it computed, and fc_master.cpp compares that against
// the size its table gives for this file's VERSION and answers FC_ERR_STRUCT_SIZE when they differ — before
// reading one further byte. So a field added, dropped or mistyped here is a REFUSAL on the first call, not a
// plausible parameter set. `fc_master_sizeof(id, version)` is checked against every struct at load besides,
// which turns the refusal into a message that says what is wrong rather than which call noticed.
//
// What a size cannot catch is a permutation that preserves the total — two f64 fields swapped. That one is
// caught by tools/wasm/layout-check.mjs, which holds every field's offset here against the offset the compiler
// gave it (`fcore_master layout`), in both directions, under ctest.
//
// ==================================================================================================
// VERSIONS (tools/fc_master_abi.h, VERSIONING — read the rule there)
// ==================================================================================================
// This file describes ONE version, FC_MASTER_ABI_VERSION below, and every struct it stamps is that version. A
// module at the same version or newer reads it (a newer module reads an older caller's struct over its own
// defaults); an older module is refused at load. `set()` and `get()` refuse any field of a struct whose header is
// not stamped at THIS file's version and size. That is stricter than "a field past the stamp", on purpose: a page
// written against one version has no business holding a struct of another, and the one way it comes to hold one is
// a FROZEN v1 `_fc_*_default` writer — which the strict rule reports at the first field written (`sampleRate`),
// instead of at the first newer one (`deliveryRate`, `compressorMix`), where the module would silently not read
// it. Stamp with `init()` and fill with `_fc_*_defaults`.

// ── the layout algorithm ──────────────────────────────────────────────────────────────────────────
const SCALAR = { i32: 4, u32: 4, f32: 4, f64: 8, u64: 8 };

// Each entry: [name, type, count?]. `type` is a scalar key above or another struct's name.
// Copied from tools/fc_master_abi.h, in order. Nothing is summarised and nothing is renamed — the ABI
// says its own surface may not be a narrower road than the C++ API, and a partial mirror here would
// make it one the moment a page wanted a field this file had skipped.
const STRUCTS = {
    fc_header: [['abiVersion', 'u32'], ['structSize', 'u32']],

    fc_master_config: [
        ['header', 'fc_header'],
        ['sampleRate', 'f64'], ['channels', 'i32'], ['internalBlock', 'i32'],
        ['eq', 'i32'], ['monoBass', 'i32'], ['compressor', 'i32'],
        ['clipper', 'i32'], ['limiter', 'i32'], ['dither', 'i32'],
        ['compressorLookaheadMs', 'f64'], ['limiterLookaheadMs', 'f64'],
        ['oversampleFactor', 'i32'], ['tapsPerPhase', 'i32'], ['sidechainHpfHz', 'f64'],
        ['deliveryRate', 'f64'],                                // v2
    ],

    fc_eq_lane: [
        ['on', 'i32'], ['freq', 'f64'], ['q', 'f64'], ['gainDb', 'f64'],
        ['slope', 'i32'], ['bypass', 'i32'],
    ],
    fc_eq_dyn: [
        ['on', 'i32'], ['rangeDb', 'f64'], ['thrDb', 'f64'], ['thrAuto', 'i32'],
        ['atk', 'f64'], ['rel', 'f64'],
    ],
    fc_eq_band: [
        ['on', 'i32'], ['type', 'i32'], ['swept', 'i32'], ['bypass', 'i32'],
        ['dyn', 'fc_eq_dyn'], ['lanes', 'fc_eq_lane', 5],       // FC_MAX_EQ_LANES
    ],
    fc_mono_bass:  [['enabled', 'i32'], ['frequencyHz', 'f32'], ['lowWidth', 'f32']],
    fc_compressor: [
        ['detector', 'i32'], ['link', 'i32'], ['rmsWindowMs', 'f64'], ['mode', 'i32'],
        ['thresholdDb', 'f64'], ['ratio', 'f64'], ['kneeDb', 'f64'], ['rangeDb', 'f64'],
        ['attackMs', 'f64'], ['releaseMs', 'f64'], ['makeupDb', 'f64'], ['autoMakeup', 'i32'],
    ],
    fc_clipper: [
        ['shape', 'i32'], ['driveDb', 'f32'], ['bias', 'f32'], ['mix', 'f32'],
        ['outputDb', 'f32'], ['autoComp', 'f32'], ['dcBlockHz', 'f32'],
    ],
    fc_limiter: [['ceilingDbTp', 'f64'], ['releaseMs', 'f64']],
    fc_dither: [
        ['bits', 'i32'], ['shaping', 'i32'],
        ['seedLo', 'u32'], ['seedHi', 'u32'],       // 64 bits, split: a Number cannot hold the low bits
        ['autoBlank', 'i32'], ['autoBlankSamples', 'i32'],
    ],

    fc_master_params: [
        ['header', 'fc_header'],
        ['inputGainDb', 'f64'], ['preLimiterGainDb', 'f64'],
        ['eqBands', 'fc_eq_band', 24],                          // FC_MAX_EQ_BANDS
        ['monoBass', 'fc_mono_bass'], ['compressor', 'fc_compressor'],
        ['clipper', 'fc_clipper'], ['limiter', 'fc_limiter'], ['dither', 'fc_dither'],
        ['bypassEq', 'i32'], ['bypassMonoBass', 'i32'], ['bypassCompressor', 'i32'],
        ['bypassClipper', 'i32'], ['bypassLimiter', 'i32'], ['bypassDither', 'i32'],
        ['compressorMix', 'f64'],                               // v3
        ['limiterDualRelease', 'i32'], ['_pad0', 'i32'],        // v6
        ['limiterSlowReleaseMs', 'f64'],                        // v6
    ],

    fc_master_resolved: [
        ['header', 'fc_header'],
        ['latencySamples', 'i32'], ['internalBlock', 'i32'], ['compressorLookahead', 'i32'],
        ['clipperLatency', 'i32'], ['limiterLatency', 'i32'], ['limiterLookahead', 'i32'],
        ['oversampleFactor', 'i32'], ['compressorTapOffset', 'i32'], ['limiterTapOffset', 'i32'],
        ['limiterCeilingDbTp', 'f64'], ['limiterReleaseMs', 'f64'],
        ['monoBass', 'fc_mono_bass'], ['tapOversampleFactor', 'i32'],
        ['compressorMix', 'f64'],                               // v3
        ['limiterSlowReleaseMs', 'f64'],                        // v6
    ],

    fc_master_stats: [
        ['header', 'fc_header'],
        ['framesIn', 'u64'], ['framesFlushed', 'u64'], ['nonFiniteIn', 'u64'],
    ],

    fc_need: [
        ['header', 'fc_header'],
        ['callBytes', 'u64'], ['solverPrepareBytes', 'u64'], ['facadeBytes', 'u64'],
        ['solverPrepared', 'i32'], ['_pad0', 'i32'],            // the tail padding, named (VERSIONING rule 4)
    ],

    fc_gr_limit: [['limitDb', 'f64'], ['statistic', 'i32']],

    fc_loudness_request: [
        ['header', 'fc_header'],
        ['targetLufs', 'f64'], ['toleranceLu', 'f64'], ['maxTruePeakDbTp', 'f64'],
        ['truePeakAimDb', 'f64'],
        ['limiterGr', 'fc_gr_limit'], ['compressorGr', 'fc_gr_limit'],
        ['minPlrDb', 'f64'], ['maxLraLossLu', 'f64'], ['inputLoudnessRangeLu', 'f64'],
        ['activityThresholdDb', 'f64'], ['maxPasses', 'i32'], ['initialGainDb', 'f64'],
        ['grTraceBuckets', 'i32'], ['_pad0', 'i32'],            // v6
        ['limiterGrQuantile', 'f64'], ['compressorGrQuantile', 'f64'],   // v8
        ['limiterActiveInputDb', 'f64'],                                 // v10 — K11's gate, dBFS
    ],

    fc_solve_pass: [
        ['gainDb', 'f64'], ['ceilingDb', 'f64'], ['integratedLufs', 'f64'],
        ['truePeakDbTp', 'f64'], ['plrDb', 'f64'], ['limiterMaxGrDb', 'f64'],
        ['loudnessRangeLu', 'f64'], ['violated', 'u32'],
    ],

    fc_gr_stats: [
        ['meanDb', 'f64'], ['p95Db', 'f64'], ['maxDb', 'f64'], ['activeFraction', 'f64'],
        ['frames', 'u64'], ['nonFinite', 'u64'], ['aboveRange', 'u64'], ['valid', 'i32'],
    ],

    fc_measurement: [
        ['header', 'fc_header'],
        ['integratedLufs', 'f64'], ['truePeakDbTp', 'f64'], ['samplePeakDb', 'f64'],
        ['loudnessRangeLu', 'f64'], ['plrDb', 'f64'],
        ['compressor', 'fc_gr_stats'], ['limiter', 'fc_gr_stats'],
        ['limiterMaxReconstructedPeakDb', 'f64'],
        ['latencySamples', 'i32'], ['gatingBlocks', 'i32'], ['droppedBlocks', 'i32'],
        ['nonFiniteSubHops', 'i32'], ['loudnessValid', 'i32'], ['lraValid', 'i32'],
        ['compressorGrTraceBuckets', 'i32'], ['limiterGrTraceBuckets', 'i32'],            // v4
        ['compressorGrTraceValid', 'i32'], ['limiterGrTraceValid', 'i32'],                // v4
    ],

    // v10 — `_fc_solution_gr_active_stats`: the limiter's statistics over the windows its INPUT reached the gate.
    // `stats` is the same frozen fc_gr_stats, every field of it over the accepted windows only.
    fc_gr_active_stats: [
        ['header', 'fc_header'],
        ['stats', 'fc_gr_stats'],
        ['windows', 'u64'], ['activeWindows', 'u64'], ['thresholdDb', 'f64'],
    ],

    // v4 — one bucket of `_fc_solution_gr_trace`, header-less (read with a stride of its size, like fc_solve_pass).
    fc_gr_trace_bucket: [
        ['maxDb', 'f64'], ['meanDb', 'f64'], ['samples', 'u32'], ['nonFinite', 'u32'],
    ],

    // v6 — one bucket of `_fc_solution_gr_trace64`, header-less.
    fc_gr_trace_bucket64: [
        ['maxDb', 'f64'], ['meanDb', 'f64'], ['samples', 'u64'], ['nonFinite', 'u64'],
    ],

    fc_progress: [
        ['stage', 'i32'], ['pass', 'i32'], ['maxPasses', 'i32'], ['hasRecord', 'i32'],
        ['fraction', 'f64'], ['record', 'fc_solve_pass'],
    ],

    fc_solution_summary: [
        ['header', 'fc_header'],
        ['status', 'i32'], ['binding', 'i32'], ['alsoViolated', 'u32'],
        ['preLimiterGainDb', 'f64'], ['ceilingDbTp', 'f64'],
        ['passes', 'i32'], ['logCount', 'i32'], ['activityThresholdDb', 'f64'],
        ['achievedBelowLufs', 'f64'], ['achievedAboveLufs', 'f64'],
        ['gainBelowDb', 'f64'], ['gainAboveDb', 'f64'],
    ],
};

// fc_struct_id — the codes `_fc_master_sizeof(id, version)` takes, for every struct that begins with a header.
export const STRUCT_IDS = {
    fc_master_config: 0, fc_master_params: 1, fc_master_resolved: 2, fc_master_stats: 3,
    fc_need: 4, fc_loudness_request: 5, fc_measurement: 6, fc_solution_summary: 7,
};

export const structNames = () => Object.keys(STRUCTS);

const layouts = new Map();

// Returns { size, align, fields: Map<name, {offset, type, count, stride}> } — computed once, cached.
export function layoutOf (name) {
    const hit = layouts.get(name);
    if (hit) return hit;
    const def = STRUCTS[name];
    if (!def) throw new Error(`fc-master-layout: no such struct ${name}`);

    // Placed in the map before the fields are walked so a cycle is a stack overflow at build time
    // rather than an infinite loop; there are none, and there is to be none.
    const out = { size: 0, align: 1, fields: new Map() };
    layouts.set(name, out);

    let off = 0;
    for (const [field, type, count] of def) {
        const el = SCALAR[type] !== undefined
            ? { size: SCALAR[type], align: SCALAR[type] }
            : layoutOf(type);
        off = align(off, el.align);
        out.fields.set(field, { offset: off, type, count: count || 0, stride: el.size });
        off += el.size * (count || 1);
        if (el.align > out.align) out.align = el.align;
    }
    out.size = align(off, out.align);
    return out;
}

const align = (n, a) => (n + a - 1) & ~(a - 1);

export const sizeOf = name => layoutOf(name).size;

// ── the accessor ──────────────────────────────────────────────────────────────────────────────────
//
// THE DETACHED-VIEW TRAP, the same one parity.mjs documents. With ALLOW_MEMORY_GROWTH a `memory.grow`
// replaces the underlying ArrayBuffer and every view onto the old one detaches: writes vanish, reads
// throw. `prepare()` and `fc_master_solve` both grow it. So this class holds NO view across a call —
// it re-reads the module's heap and builds a fresh DataView on every single access. A DataView costs
// nothing to construct next to a wasm render, and a rule that depends on nobody reordering the code is
// not a rule.
export class Struct {
    constructor (module, name, ptr) {
        this.m = module;
        this.name = name;
        this.ptr = ptr;
        this.layout = layoutOf(name);
    }

    // The exact address of a field, following dots and [i] through nested structs and arrays:
    //   s.addr('eqBands[3].lanes[0].gainDb')
    addr (path) {
        let base = this.ptr, layout = this.layout;
        for (const step of path.split('.')) {
            const m = /^([A-Za-z_]\w*)(?:\[(\d+)])?$/.exec(step);
            if (!m) throw new Error(`fc-master-layout: bad path step '${step}' in '${path}'`);
            const f = layout.fields.get(m[1]);
            if (!f) throw new Error(`fc-master-layout: ${layout === this.layout ? this.name : '(nested)'} has no field '${m[1]}'`);
            const idx = m[2] === undefined ? 0 : Number(m[2]);
            // An index past the array would address the NEXT field and write a plausible value into it.
            if (idx !== 0 && idx >= f.count) throw new Error(`fc-master-layout: '${m[1]}[${idx}]' is past its ${f.count}`);
            base += f.offset + idx * f.stride;
            if (SCALAR[f.type] === undefined) layout = layoutOf(f.type);
            else { layout = null; this._t = f.type; }
        }
        return base;
    }

    _view () { return new DataView(this.m.HEAPF32.buffer); }

    // The stamp check `get` and `set` share — see VERSIONS at the top of this file.
    _stamped (path, v) {
        if (STRUCT_IDS[this.name] === undefined || path.startsWith('header')) return;
        const ver = v.getUint32(this.ptr, true), size = v.getUint32(this.ptr + 4, true);
        if (ver !== FC_MASTER_ABI_VERSION || size !== this.layout.size)
            throw new Error(`fc-master-layout: ${this.name} is stamped v${ver}/${size} B, not this file's `
                          + `v${FC_MASTER_ABI_VERSION}/${this.layout.size} B — stamp it with init() (and fill an input `
                          + `struct with _${this.name}_defaults, not a frozen v1 _default writer) before '${path}'`);
    }

    get (path) {
        const at = this.addr(path), t = this._t, v = this._view();
        this._stamped(path, v);
        switch (t) {
            case 'i32': return v.getInt32(at, true);
            case 'u32': return v.getUint32(at, true);
            case 'f32': return v.getFloat32(at, true);
            case 'f64': return v.getFloat64(at, true);
            case 'u64': return v.getBigUint64(at, true);
            default: throw new Error(`fc-master-layout: '${path}' is a struct, not a value`);
        }
    }

    set (path, value) {
        const at = this.addr(path), t = this._t, v = this._view();
        this._stamped(path, v);
        switch (t) {
            case 'i32': v.setInt32(at, value, true); break;
            case 'u32': v.setUint32(at, value, true); break;
            case 'f32': v.setFloat32(at, value, true); break;
            case 'f64': v.setFloat64(at, value, true); break;
            case 'u64': v.setBigUint64(at, BigInt(value), true); break;
            default: throw new Error(`fc-master-layout: '${path}' is a struct, not a value`);
        }
        return this;
    }

    setAll (obj) { for (const k of Object.keys(obj)) this.set(k, obj[k]); return this; }

    // FC_INIT: the two stores every struct on this boundary begins with. An OUT struct needs them too —
    // the caller is the one saying which layout it has room for.
    init () {
        const v = this._view();
        v.setUint32(this.ptr + 0, FC_MASTER_ABI_VERSION, true);
        v.setUint32(this.ptr + 4, this.layout.size, true);
        return this;
    }

    // Everything in the struct as a plain object, arrays included. For reporting, never for a hot loop.
    toJSON (name = this.name, base = this.ptr) {
        const out = {};
        for (const [field, f] of layoutOf(name).fields) {
            const one = (at) => {
                if (SCALAR[f.type] !== undefined) {
                    const v = this._view();
                    switch (f.type) {
                        case 'i32': return v.getInt32(at, true);
                        case 'u32': return v.getUint32(at, true);
                        case 'f32': return v.getFloat32(at, true);
                        case 'f64': return v.getFloat64(at, true);
                        case 'u64': return v.getBigUint64(at, true);
                    }
                }
                return this.toJSON(f.type, at);
            };
            out[field] = f.count
                ? Array.from({ length: f.count }, (_, i) => one(base + f.offset + i * f.stride))
                : one(base + f.offset);
        }
        return out;
    }
}

export const FC_MASTER_ABI_VERSION = 10;

// The status codes, in the order fc_master_abi.h declares them — so a refusal reaches a human as a name.
export const FC_STATUS = [
    'FC_OK', 'FC_ERR_HANDLE', 'FC_ERR_ABI_VERSION', 'FC_ERR_STRUCT_SIZE', 'FC_ERR_NULL',
    'FC_ERR_ALIGNMENT', 'FC_ERR_SPAN', 'FC_ERR_ENUM', 'FC_ERR_RANGE', 'FC_ERR_CAPACITY',
    'FC_ERR_STATE', 'FC_ERR_NON_FINITE', 'FC_ERR_REFUSED_BY_CORE', 'FC_ERR_EXHAUSTED',
    'FC_ERR_POISONED', 'FC_ERR_CANCELLED',
];
export const statusName = s => FC_STATUS[s] ?? `FC_STATUS(${s})`;

export const FC_SOLVE_STATUS = [
    'Solved', 'TargetUnreachable', 'UpstreamViolation', 'TargetBetween', 'PassLimit',
    'MeasurementInvalid', 'RenderFailed', 'NotPrepared', 'InvalidRequest',
];
export const FC_CONSTRAINT = [
    'None', 'TruePeak', 'LimiterGr', 'Plr', 'Lra', 'GainRange', 'CompressorGr',
];

// `fc_solution_summary.alsoViolated` — the INDEX is the bit position and the value is the constraint that
// bit stands for, so bit `i` is `FC_CONSTRAINT[i + 1]` and `None` has no bit. `binding` is INCLUDED in the
// mask. Held against the core's `constraintBit()` by layout-check.mjs.
export const FC_CONSTRAINT_BITS = FC_CONSTRAINT.slice(1);
export const constraintsOf = mask => FC_CONSTRAINT_BITS.filter((_, bit) => (mask & (1 << bit)) !== 0);

// ── the enum codes ────────────────────────────────────────────────────────────────────────────────
//
// One list per enum a caller writes into a struct or passes to an entry point. THE INDEX IS THE CODE. A name is
// the header suffix with the underscores dropped and the words capitalised; `FC_ENUMS` below names the header
// enum and the prefix each list mirrors, and layout-check.mjs holds every list against it — order, value and
// letters, in both directions.
export const FC_EQ_AXIS       = ['Stereo', 'Left', 'Right', 'Mid', 'Side'];                    // v7, `lane` of eq_curve
export const FC_FILTER_TYPE   = ['Bell', 'LowShelf', 'HighShelf', 'HighPass', 'LowPass',
                                 'BandPass', 'Notch', 'AllPass', 'Tilt'];                      // fc_eq_band.type
export const FC_DETECTOR      = ['Peak', 'Rms'];                                               // fc_compressor.detector
export const FC_LINK_MODE     = ['Max', 'MeanPower'];                                          // fc_compressor.link
export const FC_COMP_MODE     = ['DownCompress', 'UpCompress', 'DownExpand'];                  // fc_compressor.mode
export const FC_SHAPE         = ['Tanh', 'Atan', 'Cubic', 'Asym'];                             // fc_clipper.shape
export const FC_NOISE_SHAPING = ['None', 'Weighted', 'Psycho'];                                // fc_dither.shaping
export const FC_GR_STATISTIC  = ['Mean', 'P95', 'Max', 'Percentile'];                          // fc_gr_limit.statistic
export const FC_GR_STAGE      = ['Compressor', 'Limiter'];                                     // `stage` of gr_trace
export const FC_PROGRESS_STAGE= ['Convert', 'Lra', 'Pass', 'Final', 'Render'];                 // fc_progress.stage

// Which header enum each list mirrors, and under which prefix — the gate reads this rather than a list of its
// own, so a list added here without a header enum to hold it against is a failure and not an omission.
export const FC_ENUMS = {
    FC_EQ_AXIS:        { enum: 'fc_eq_axis',        prefix: 'FC_EQ_AXIS_',   names: FC_EQ_AXIS },
    FC_FILTER_TYPE:    { enum: 'fc_filter_type',    prefix: 'FC_FILTER_',    names: FC_FILTER_TYPE },
    FC_DETECTOR:       { enum: 'fc_detector',       prefix: 'FC_DETECTOR_',  names: FC_DETECTOR },
    FC_LINK_MODE:      { enum: 'fc_link_mode',      prefix: 'FC_LINK_',      names: FC_LINK_MODE },
    FC_COMP_MODE:      { enum: 'fc_comp_mode',      prefix: 'FC_COMP_',      names: FC_COMP_MODE },
    FC_SHAPE:          { enum: 'fc_shape',          prefix: 'FC_SHAPE_',     names: FC_SHAPE },
    FC_NOISE_SHAPING:  { enum: 'fc_noise_shaping',  prefix: 'FC_SHAPING_',   names: FC_NOISE_SHAPING },
    FC_GR_STATISTIC:   { enum: 'fc_gr_statistic',   prefix: 'FC_GR_',        names: FC_GR_STATISTIC },
    FC_GR_STAGE:       { enum: 'fc_gr_stage',       prefix: 'FC_GR_STAGE_',  names: FC_GR_STAGE },
    FC_PROGRESS_STAGE: { enum: 'fc_progress_stage', prefix: 'FC_PROGRESS_',  names: FC_PROGRESS_STAGE },
    FC_SOLVE_STATUS:   { enum: 'fc_solve_status',   prefix: 'FC_SOLVE_',     names: FC_SOLVE_STATUS },
    FC_CONSTRAINT:     { enum: 'fc_constraint',     prefix: 'FC_CONSTRAINT_', names: FC_CONSTRAINT },
};

// ==================================================================================================
// THE FIELD DOMAINS
// ==================================================================================================
// What every input field of `fc_master_config`, `fc_master_params` and `fc_loudness_request` ADMITS, where the
// boundary is, and what a value on the far side of it does. Every number is the code's: the checks in `toCore`
// and `fc_master_create` (tools/wasm/fc_master.cpp), each module's own `prepare`/`setParams` clamps, and
// `TargetLoudnessSolver::admits`. A field the code does not bound is written down as unbounded.
//
// NO ROW HAS A LIST-VALUED DOMAIN: every one is an interval, an enum code range, or unbounded. Two that look
// list-valued are not — `oversampleFactor` admits every INTEGER in its interval (the mastering chain builds the
// Kaiser oversampler, which has no power-of-two restriction), and `dither.bits` refuses nothing.
//
// tools/tests/MasterDomainsTests.cpp holds this table against the running ABI with a real call per bound.
//
// THE COLUMNS
//   field      the path from the struct root; `[]` stands for EVERY element (all 24 bands, all 5 lanes)
//   unit       dB · dBTP · dB/oct · LUFS · LU · Hz · ms · x (a multiplier) · fraction · frames · samples ·
//              count · bits · code (an opaque 32-bit value) · flag (0 / non-0) · enum:<NAME> (an index into
//              that list above)
//   min, max   the bound, `null` where the code has none, or a formula in `sr` — the CHAIN sample rate (the
//              delivery rate on a delivering handle) — and `os`, `fc_master_config.oversampleFactor`, which
//              `domainBound()` below evaluates
//   open       which bound is EXCLUSIVE ('max' means the bound itself is already outside), '' when both are inclusive
//   edge       what a value outside the interval does:
//                'refuse'  the call is refused with `err`, nothing moves
//                'clamp'   ⚠ the value is SILENTLY pulled to the bound and the call succeeds — `resolved` says
//                          where the applied value can be read back, and '' means it cannot be read anywhere
//                'verdict' `fc_master_solve` still answers FC_OK; the SOLUTION carries `err` (an FC_SOLVE_* name),
//                          so this one has to be caught by reading the summary, not the status
//                'free'    the code bounds nothing — only the finiteness check in `nonFinite` stands
//                'any'     every value of the field type means something (a flag read as `!= 0`, an opaque code)
//   err        the refusal that names this boundary, '' where there is none
//   nonFinite  'refuse'      NaN and both infinities are FC_ERR_NON_FINITE
//              'verdict'     any non-finite value is FC_SOLVE_INVALID_REQUEST
//              'nan-verdict' a NaN is FC_SOLVE_INVALID_REQUEST; the infinities are ordinary values of this field
//              'off'         a non-finite value switches the field off / means "not supplied"
//              'none'        an integer field
//   resolved   where the APPLIED value can be read back: a field of `fc_master_resolved`, `summary.<field>`
//              (a field of the solution `fc_solution_summary_get` answers), `eqCurve` (the magnitude
//              `fc_master_eq_curve` answers), `render` (what the chain puts out), or '' — the applied value
//              is not observable through this ABI at all
//   depends    what else moves this domain — the sample rate, or another field
//
// 'refuse' AGAINST 'clamp' IS THE DISTINCTION AN INTERFACE NEEDS. A refusal is a status on the call that made
// it. A clamp is silent, and shows only in what was rendered: every 'clamp' row is a control that must bound
// itself, because the core will not complain.
export const FC_DOMAINS = [
    // ── fc_master_config — the topology, fixed for the life of a handle (`fc_master_create`) ──────
    { field: 'fc_master_config.sampleRate', unit: 'Hz', min: 8000, max: 3e6, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'refuse', resolved: '', depends: 'deliveryRate: when it is not 0 this is the SOURCE rate and the chain runs at the delivery rate, and the pair must have a resampling route' },
    { field: 'fc_master_config.channels', unit: 'count', min: 1, max: 16, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'none', resolved: '', depends: 'monoBass: with the mono-bass stage on the only admitted width is 2' },
    { field: 'fc_master_config.internalBlock', unit: 'frames', min: 8, max: 8192, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'none', resolved: 'internalBlock', depends: '' },
    { field: 'fc_master_config.eq', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_config.monoBass', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'none', resolved: '', depends: 'channels: on with a width other than 2 the create is refused, because the stage would leave the buffer untouched' },
    { field: 'fc_master_config.compressor', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_config.clipper', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_config.limiter', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_config.dither', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_config.compressorLookaheadMs', unit: 'ms', min: 0, max: 250, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'refuse', resolved: '', depends: 'compressor: the 250 ms ceiling is the compressor stage own and is only reached with the stage on; a negative value is refused either way' },
    { field: 'fc_master_config.limiterLookaheadMs', unit: 'ms', min: '2000/sr', max: 20, open: '', edge: 'clamp', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'refuse', resolved: '', depends: 'sampleRate: the floor is 2 baseband samples. A NEGATIVE value is refused; everything from 0 up is clamped into the interval, and the applied lookahead is readable only in SAMPLES as resolved.limiterLookahead' },
    { field: 'fc_master_config.oversampleFactor', unit: 'x', min: 2, max: 16, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'none', resolved: 'oversampleFactor', depends: 'limiter and clipper: 16 is the limiter ceiling, 64 the clipper ceiling with no limiter, and with both stages off nothing above 2 is refused at all - and the resolved value is then 0, since no stage oversamples' },
    { field: 'fc_master_config.tapsPerPhase', unit: 'count', min: 4, max: 1024, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'none', resolved: '', depends: 'limiter and clipper: the 1024 ceiling is the oversampler own, so with both stages off only the floor of 4 stands' },
    { field: 'fc_master_config.sidechainHpfHz', unit: 'Hz', min: 0, max: '0.5*sr', open: 'max', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'refuse', resolved: '', depends: 'sampleRate: the Nyquist itself is already refused' },
    { field: 'fc_master_config.deliveryRate', unit: 'Hz', min: 8000, max: 3e6, open: '', edge: 'refuse', err: 'FC_ERR_REFUSED_BY_CORE', nonFinite: 'refuse', resolved: '', depends: '0 means no conversion. Inside the interval the value must also be a WHOLE number of hertz and share a resampling route with sampleRate, so the interval is necessary and not sufficient' },

    // ── fc_master_params — per-configure, and every numeric field of it is CLAMPED rather than refused ──
    // The only refusals this struct can make are the two the ABI makes on its own: a non-finite number, and an
    // enum code this version does not define. Nothing in `fc_master_configure` below the mapping can say no —
    // the chain re-prepares at the geometry the handle already has — so a value out of range is a clamp, always.
    { field: 'fc_master_params.inputGainDb', unit: 'dB', min: -60, max: 60, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.preLimiterGainDb', unit: 'dB', min: -60, max: 60, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.compressorMix', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'compressorMix', depends: 'compressor: 0 is reported without the stage, whatever was asked for' },
    { field: 'fc_master_params.limiterDualRelease', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.limiterSlowReleaseMs', unit: 'ms', min: '8000/sr', max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'limiterSlowReleaseMs', depends: 'sampleRate: the floor is 8 baseband samples. limiterDualRelease: the resolved value reads 0 while the second envelope is off' },
    { field: 'fc_master_params.bypassEq', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.bypassMonoBass', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.bypassCompressor', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.bypassClipper', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.bypassLimiter', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.bypassDither', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },

    // ── fc_master_params.eqBands[] — every band, and `lanes[]` every lane of it ───────────────────
    { field: 'fc_master_params.eqBands[].on', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].type', unit: 'enum:FC_FILTER_TYPE', min: 0, max: 8, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].swept', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].bypass', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].dyn.on', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].dyn.rangeDb', unit: 'dB', min: -30, max: 30, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'dyn.on gates the whole group. The SIGN chooses the direction - negative cuts as the band gets loud, positive boosts - and the MAGNITUDE is what the interval bounds. 0 is no dynamics whatever dyn.on says' },
    { field: 'fc_master_params.eqBands[].dyn.thrDb', unit: 'dB', min: -120, max: 24, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'dyn.thrAuto: an ABSOLUTE dBFS threshold, read ONLY while the automatic threshold is off. What it is measured against is the band probe level, so inputGainDb moves what crossing it means' },
    { field: 'fc_master_params.eqBands[].dyn.thrAuto', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].dyn.atk', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'dyn.on, and a delta that MOVES: 0.5 is the automatic value derived from the lane fc/Q and the ends are deviations from it (1/4x .. 4x), so a band whose reduction never changes renders the same at every setting' },
    { field: 'fc_master_params.eqBands[].dyn.rel', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'as dyn.atk, and never shorter than the attack it follows' },
    { field: 'fc_master_params.eqBands[].lanes[].on', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.eqBands[].lanes[].freq', unit: 'Hz', min: 10, max: '0.49*sr', open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'eqCurve', depends: 'sampleRate' },
    { field: 'fc_master_params.eqBands[].lanes[].q', unit: 'x', min: 0.05, max: 40, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'eqCurve', depends: '' },
    { field: 'fc_master_params.eqBands[].lanes[].gainDb', unit: 'dB', min: -30, max: 30, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'eqCurve', depends: '' },
    { field: 'fc_master_params.eqBands[].lanes[].slope', unit: 'dB/oct', min: 6, max: 96, open: '', edge: 'clamp', err: '', nonFinite: 'none', resolved: '', depends: 'eqBands[].type: read only by the cut, notch and band-pass types, and read as slope/6 POLES clamped to [1, 16] — so only multiples of 6 are distinct and everything below 6 behaves as 6' },
    { field: 'fc_master_params.eqBands[].lanes[].bypass', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },

    // ── fc_master_params.monoBass ─────────────────────────────────────────────────────────────────
    { field: 'fc_master_params.monoBass.enabled', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: 'monoBass.enabled', depends: 'fc_master_config.monoBass: the resolved flag is the AND of the two, so it reads 0 whenever the stage is not in the chain' },
    { field: 'fc_master_params.monoBass.frequencyHz', unit: 'Hz', min: 20, max: '0.45*sr', open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'monoBass.frequencyHz', depends: 'sampleRate' },
    { field: 'fc_master_params.monoBass.lowWidth', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'monoBass.lowWidth', depends: '' },

    // ── fc_master_params.compressor ───────────────────────────────────────────────────────────────
    { field: 'fc_master_params.compressor.detector', unit: 'enum:FC_DETECTOR', min: 0, max: 1, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.compressor.link', unit: 'enum:FC_LINK_MODE', min: 0, max: 1, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.compressor.rmsWindowMs', unit: 'ms', min: 0, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'compressor.detector: read only by the Rms detector. At or below 0 the envelope is instant' },
    { field: 'fc_master_params.compressor.mode', unit: 'enum:FC_COMP_MODE', min: 0, max: 2, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.compressor.thresholdDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'refuse', resolved: '', depends: 'the detector level enters the curve through a floor of -240 dB, so a threshold below that makes digital silence an ACTIVE sample' },
    { field: 'fc_master_params.compressor.ratio', unit: 'x', min: 1, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.compressor.kneeDb', unit: 'dB', min: 0, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.compressor.rangeDb', unit: 'dB', min: 0, max: 400, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.compressor.attackMs', unit: 'ms', min: 0, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'at or below 0 the ballistics are instant, and a time long enough to round the coefficient to 1 is backed off so the envelope never freezes' },
    { field: 'fc_master_params.compressor.releaseMs', unit: 'ms', min: 0, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'as attackMs' },
    { field: 'fc_master_params.compressor.makeupDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'refuse', resolved: '', depends: 'what is bounded is the SUM of this and the gain reduction, inside the render, and not this field' },
    { field: 'fc_master_params.compressor.autoMakeup', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },

    // ── fc_master_params.clipper ──────────────────────────────────────────────────────────────────
    { field: 'fc_master_params.clipper.shape', unit: 'enum:FC_SHAPE', min: 0, max: 3, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.clipper.driveDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'refuse', resolved: '', depends: 'the drive the shaper runs is dbToGain(driveDb) - 1 floored at 1e-4, so every driveDb at or below 20*log10(1 + 1e-4) = 8.685455e-4 dB is the same linear stage. The figure this row used to give, 0.00087, is that threshold ROUNDED THE WRONG WAY: 8.7e-4 is past it, so the value the row named as identical to zero drive already renders differently. A threshold in prose rounds toward the safe side or not at all. NB the floored stage is still not a BYPASSED one — it runs the oversampler round trip; only bypassClipper skips that, and clipper.mix = 0 does not (see MasteringChain.h)' },
    { field: 'fc_master_params.clipper.bias', unit: 'fraction', min: -0.95, max: 0.95, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'clipper.shape: read only by Asym' },
    { field: 'fc_master_params.clipper.mix', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.clipper.outputDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'refuse', resolved: '', depends: '' },
    { field: 'fc_master_params.clipper.autoComp', unit: 'fraction', min: 0, max: 1, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: '' },
    { field: 'fc_master_params.clipper.dcBlockHz', unit: 'Hz', min: 0, max: '0.49*sr*os', open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'render', depends: 'sampleRate AND fc_master_config.oversampleFactor: the blocker runs in the OVERSAMPLED domain. clipper.shape gates it — only Asym enables it at all' },

    // ── fc_master_params.limiter ──────────────────────────────────────────────────────────────────
    { field: 'fc_master_params.limiter.ceilingDbTp', unit: 'dBTP', min: -200, max: 60, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'limiterCeilingDbTp', depends: '' },
    { field: 'fc_master_params.limiter.releaseMs', unit: 'ms', min: '8000/sr', max: null, open: '', edge: 'clamp', err: '', nonFinite: 'refuse', resolved: 'limiterReleaseMs', depends: 'sampleRate: the floor is 8 baseband samples. A time long enough to round the coefficient to 1 is backed off, so the recovery never stops' },

    // ── fc_master_params.dither ───────────────────────────────────────────────────────────────────
    { field: 'fc_master_params.dither.bits', unit: 'bits', min: 2, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'none', resolved: 'render', depends: 'nothing is refused. The quantiser scale is built from bits clamped to [2, 31], so everything below 2 is 2 - but the TOP is not a clamp: any value of 32 or more BYPASSES the stage entirely (a float export), which is why there is no max here' },
    { field: 'fc_master_params.dither.shaping', unit: 'enum:FC_NOISE_SHAPING', min: 0, max: 2, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.dither.seedLo', unit: 'code', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: 'the low half of the 64-bit seed; every 32-bit value is a seed' },
    { field: 'fc_master_params.dither.seedHi', unit: 'code', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: 'the high half of the 64-bit seed' },
    { field: 'fc_master_params.dither.autoBlank', unit: 'flag', min: null, max: null, open: '', edge: 'any', err: '', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_master_params.dither.autoBlankSamples', unit: 'samples', min: 1, max: null, open: '', edge: 'clamp', err: '', nonFinite: 'none', resolved: 'render', depends: '' },

    // ── fc_loudness_request — judged by the SOLVER, so its refusals arrive as a VERDICT ───────────
    // `fc_master_solve` answers FC_OK for every row below (the two enum codes excepted, which the mapping
    // refuses before the solver is reached) and the solution carries FC_SOLVE_INVALID_REQUEST. Read the summary.
    // The infinities are MEANINGFUL here and are not swept away: `limitDb` at +infinity is "no limit",
    // `minPlrDb` at -infinity is "no limit", a non-finite `inputLoudnessRangeLu` is "not supplied" and a
    // non-finite `initialGainDb` is "use the parameter set own gain".
    { field: 'fc_loudness_request.targetLufs', unit: 'LUFS', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'verdict', resolved: '', depends: 'REQUIRED: the defaults writer leaves it NaN on purpose, because a delivery target is a product policy' },
    { field: 'fc_loudness_request.toleranceLu', unit: 'LU', min: 0, max: null, open: '', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'verdict', resolved: '', depends: '' },
    { field: 'fc_loudness_request.maxTruePeakDbTp', unit: 'dBTP', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'verdict', resolved: '', depends: 'REQUIRED, and it is the DELIVERED measured peak, not the limiter setting, which the solver derives' },
    { field: 'fc_loudness_request.truePeakAimDb', unit: 'dB', min: 0, max: null, open: '', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'verdict', resolved: '', depends: '' },
    { field: 'fc_loudness_request.limiterGr.limitDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'nan-verdict', resolved: '', depends: '+infinity is the OFF value and the default; a NaN is malformed' },
    { field: 'fc_loudness_request.limiterGr.statistic', unit: 'enum:FC_GR_STATISTIC', min: 0, max: 3, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_loudness_request.compressorGr.limitDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'nan-verdict', resolved: '', depends: 'as limiterGr.limitDb' },
    { field: 'fc_loudness_request.compressorGr.statistic', unit: 'enum:FC_GR_STATISTIC', min: 0, max: 3, open: '', edge: 'refuse', err: 'FC_ERR_ENUM', nonFinite: 'none', resolved: '', depends: '' },
    { field: 'fc_loudness_request.minPlrDb', unit: 'dB', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'nan-verdict', resolved: '', depends: '-infinity is the OFF value and the default' },
    { field: 'fc_loudness_request.maxLraLossLu', unit: 'LU', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'nan-verdict', resolved: '', depends: '+infinity is the OFF value and the default. It is a LOSS against inputLoudnessRangeLu, not an absolute floor' },
    { field: 'fc_loudness_request.inputLoudnessRangeLu', unit: 'LU', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'off', resolved: '', depends: 'maxLraLossLu: it is the other end of that delta, and a non-finite value switches the pair off' },
    { field: 'fc_loudness_request.activityThresholdDb', unit: 'dB', min: 0, max: null, open: '', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'verdict', resolved: '', depends: '' },
    { field: 'fc_loudness_request.maxPasses', unit: 'count', min: 1, max: 32, open: '', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'none', resolved: '', depends: 'it bounds the SEARCH; the delivered solution can cost two renders more, and the summary passes count says so' },
    { field: 'fc_loudness_request.initialGainDb', unit: 'dB', min: -60, max: 60, open: '', edge: 'clamp', err: '', nonFinite: 'off', resolved: 'summary.preLimiterGainDb', depends: 'the search clamps the starting gain to the chain own +-60 dB before its first render. A non-finite value means start from the parameter set own preLimiterGainDb' },
    { field: 'fc_loudness_request.grTraceBuckets', unit: 'count', min: 1, max: 65536, open: '', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'none', resolved: '', depends: 'the trace actually built has min(this, programme frames) buckets' },
    { field: 'fc_loudness_request.limiterGrQuantile', unit: 'fraction', min: 0, max: 1, open: 'min', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'verdict', resolved: '', depends: 'admitted WHATEVER the statistic is, so a 0 here is refused even when the limit does not read it' },
    { field: 'fc_loudness_request.compressorGrQuantile', unit: 'fraction', min: 0, max: 1, open: 'min', edge: 'verdict', err: 'FC_SOLVE_INVALID_REQUEST', nonFinite: 'verdict', resolved: '', depends: 'as limiterGrQuantile' },
    { field: 'fc_loudness_request.limiterActiveInputDb', unit: 'dBFS', min: null, max: null, open: '', edge: 'free', err: '', nonFinite: 'off', resolved: '', depends: 'K11. The gate on the LIMITER INPUT for fc_solution_gr_active_stats, at the limiter node (after preLimiterGainDb). Decides nothing the solver judges — every limit still reads the ungated distribution. Not clamped and not refused: -inf accepts every window that carried any non-zero input, +inf accepts none, NaN reads as -inf. `off` is the nearest word this vocabulary has and its prose is looser than the truth: a non-finite value is ADMITTED and each one means something, rather than one of them switching the field off. Default -60' },
];

// A bound of FC_DOMAINS at a given CHAIN sample rate and oversample factor: a number passes through, `null`
// stays null, and the three formula forms `<a>*sr`, `<a>/sr` and `<a>*sr*os` are evaluated. There is no fourth
// form; anything else throws.
export function domainBound (bound, sampleRate, oversampleFactor = 1) {
    if (bound === null || typeof bound === 'number') return bound;
    const m = /^(-?[0-9.eE+-]+)(\*sr\*os|\*sr|\/sr)$/.exec(String(bound));
    if (!m) throw new Error(`fc-master-layout: '${bound}' is not a domain bound — use a number, null, `
                          + `<a>*sr, <a>/sr or <a>*sr*os`);
    if (m[2] === '/sr') return Number(m[1]) / sampleRate;
    return Number(m[1]) * sampleRate * (m[2] === '*sr*os' ? oversampleFactor : 1);
}

// The one check this file can make about itself before anything is rendered. Called by every loader; throws
// rather than returning, because a mismatch here means every later refusal would be FC_ERR_STRUCT_SIZE with no
// explanation attached.
//
// THE VERSION FIRST, AND "AT LEAST", NOT "EQUAL": a newer module reads this file's structs (VERSIONING rule 6) and
// has every entry point this file can call; an older one may lack both. Then EVERY struct with a header, at THIS
// file's version, from the module's own size table.
export function assertLayoutMatches (module) {
    const v = module._fc_master_abi_version();
    if (!(v >= FC_MASTER_ABI_VERSION))
        throw new Error(`fc-master-layout: ABI v${v} in the module is older than v${FC_MASTER_ABI_VERSION} here — `
                      + `the module cannot read this file's structs or has not got its entry points`);
    for (const [name, id] of Object.entries(STRUCT_IDS)) {
        const ours = sizeOf(name), theirs = module._fc_master_sizeof(id, FC_MASTER_ABI_VERSION);
        if (ours !== theirs)
            throw new Error(`fc-master-layout: ${name} is ${ours} bytes here and ${theirs} in the module at `
                          + `v${FC_MASTER_ABI_VERSION} — this file and tools/fc_master_abi.h have fallen out of step`);
    }
}
