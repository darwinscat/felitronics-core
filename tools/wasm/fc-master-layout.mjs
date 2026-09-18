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

export const FC_MASTER_ABI_VERSION = 7;

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
