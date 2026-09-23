// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The wasm side of P81's `fc_probe_<mode>_storage_bytes` — the price of a measurement, asked before it is
// paid.
//
//   node storage-probe.mjs build/fcprobe.node.js table    # the same 475 rows the native tier prints
//   node storage-probe.mjs build/fcprobe.node.js check    # the assertions that only this tier can make
//
// WHY A HARNESS OF ITS OWN, WHEN felitronics_analysis_abi_tests ALREADY HAMMERS THESE FUNCTIONS. That
// suite compiles fc_probe.cpp NATIVELY. It cannot see whether a name survived the link into the module,
// and `cmp` cannot detect a function nobody called — an export missing from the artifact would leave
// every native check green and the page holding `undefined`. So:
//
//   `table` prints the demands as the MODULE computes them, for
//           `felitronics_analysis_abi_tests --storage-table` to be compared against. That is a CROSS-TIER
//           oracle: wasm32 is a 32-bit target with its own std::size_t and its own struct layouts, so the
//           two tables agreeing row for row is a statement about the arithmetic and not about one build
//           being compared with itself. Only integer rates are listed, so both sides print "48000" and
//           the diff is of numbers rather than of two printf dialects.
//   `check` asserts what only this tier can: that all five names are callable on the Module at all, and
//           that in THIS tier too the price is positive exactly where `_run` is accepted.
//
// No HEAP view is held across a call into the module. An accepted `_run` may allocate, and `memory.grow`
// inside it DETACHES the ArrayBuffer every existing view was made on. A detached view is not a view onto
// freed memory: its length becomes 0, an indexed read gives `undefined`, an indexed write reaches nothing,
// and `.set()` throws TypeError. That is why every view here is taken AFTER the call and copied out with
// `.slice()` before the next one.

import { createRequire } from 'node:module';
import { resolve } from 'node:path';

const [, , modPath, cmdArg] = process.argv;
const cmd = cmdArg ?? 'table';
const refuse = m => { console.error(m); process.exit(2); };
if (!modPath || (cmd !== 'table' && cmd !== 'check'))
    refuse('usage: node storage-probe.mjs <module.js> [table|check]');

const require = createRequire(import.meta.url);
const M = await require(resolve(modPath))();

// `acceptsEmpty` is the one deliberate asymmetry between the price and the run: four of the five report
// on an empty programme, and `lowend` refuses one because `fcore_measure lowend` refuses it too.
const MODES = ['report', 'bursts', 'hum', 'forensics', 'lowend', 'stereobursts'];
const ACCEPTS_EMPTY = { report: true, bursts: true, hum: true, forensics: true, lowend: false,
                        stereobursts: true };

// THE FIRST ASSERTION, AND IT IS MADE IN BOTH MODES. A name that did not reach the artifact is `undefined`
// here, and calling it would throw a TypeError three lines later with nothing to say about why.
const query = {}, run = {};
for (const m of MODES) {
    const q = M[`_fc_probe_${m}_storage_bytes`], r = M[`_fc_probe_${m}_run`];
    if (typeof q !== 'function') refuse(`_fc_probe_${m}_storage_bytes is ${typeof q} on this module — it did not reach the artifact`);
    if (typeof r !== 'function') refuse(`_fc_probe_${m}_run is ${typeof r} on this module`);
    query[m] = q; run[m] = r;
}

// The table the native tier prints, row for row. Integer rates only — see the header.
const TABLE_WIDTHS = [1, 2, 6, 16, 17];
const TABLE_RATES = [999, 1000, 7999, 8000, 11025, 12000, 16000, 22050, 44100, 48000,
                     88200, 96000, 176400, 192000, 352800, 384000, 705600, 768000, 768001];

if (cmd === 'table') {
    let out = '';
    for (const m of MODES)
        for (const ch of TABLE_WIDTHS)
            for (const sr of TABLE_RATES)
                out += `${m} ${ch} ${sr} ${query[m](ch, sr)}\n`;
    process.stdout.write(out);
    process.exit(0);
}

// ---- check ----------------------------------------------------------------------------------------
let checks = 0, bad = 0;
const ok = (cond, what) => { ++checks; if (!cond) { ++bad; console.error(`FAIL: ${what}`); } };

// K1 — the crest price and its two entry points, on the artifact. Its geometry includes the programme LENGTH
// (the cell store is per hop), which is why it takes a third argument and cannot be a row of the table above.
{
    const q = M._fc_probe_crest_storage_bytes, qw = M._fc_probe_crest_storage_bytes_with;
    if (typeof q !== 'function')
        refuse(`_fc_probe_crest_storage_bytes is ${typeof q} on this module — it did not reach the artifact`);
    if (typeof qw !== 'function')
        refuse(`_fc_probe_crest_storage_bytes_with is ${typeof qw} on this module`);
    for (const n of ['_fc_probe_crest_run', '_fc_probe_crest_run_with', '_fc_probe_crest_scalars',
                     '_fc_probe_crest_blocks', '_fc_probe_crest_loss'])
        if (typeof M[n] !== 'function') refuse(`${n} is ${typeof M[n]} on this module`);
    const oneSec = q(2, 48000, 48000);
    ok(oneSec > 0, 'crest: one second of stereo at 48 kHz is priced above zero');
    ok(q(2, 48000, 480000) > oneSec, 'crest: ten times the programme costs more — the store is per hop');
    ok(qw(2, 48000, 48000, 120, 2000, 6000, 100, 4, -70, -40) === oneSec,
       'crest: the parameterised price at the documented defaults is the default price');
    ok(qw(2, 48000, 48000, 3000, 2000, 6000, 100, 4, -70, -40) === 0,
       'crest: edges that do not rise are priced at the canonical zero');
    ok(q(0, 48000, 48000) === 0 && q(99, 48000, 48000) === 0,
       'crest: a width the core does not have is priced at zero');
}

// K3 — the PARAMETERISED bursts price, which has a different arity and so cannot be a row of the table
// above. It still needs a gate on the artifact: the failure that matters is the name not reaching the
// module at all, and after that, the two claims the native suite makes about it. Same shape as the rest
// of this file, and it runs before the table's own checks because it touches nothing they read.
{
    const qw = M._fc_probe_bursts_storage_bytes_with;
    if (typeof qw !== 'function')
        refuse(`_fc_probe_bursts_storage_bytes_with is ${typeof qw} on this module — it did not reach the artifact`);
    if (typeof M._fc_probe_bursts_run_with !== 'function')
        refuse(`_fc_probe_bursts_run_with is ${typeof M._fc_probe_bursts_run_with} on this module`);
    const dflt = M._fc_probe_bursts_storage_bytes(2, 48000);
    ok(qw(2, 48000, 5000, 9000, 10, 2000, 6, 3) === dflt,
       'bursts: the parameterised price at the documented defaults is the default price');
    ok(qw(2, 48000, 5000, 9000, 10, 8000, 6, 3) > dflt,
       'bursts: a four-times-longer baseline costs more than the default');
    ok(qw(2, 48000, 20000, 30000, 10, 2000, 6, 3) === 0,
       'bursts: a band the run refuses is priced at the canonical zero');
}

// A real, non-empty programme: `lowend` refuses an empty one where the other four report on it, and the
// equivalence below is about GEOMETRY, so it must not be confounded with an input-contract refusal.
const FRAMES = 512, WIDEST = 16;
const ptr = M._malloc(FRAMES * WIDEST * 4);
if (!ptr) refuse('wasm OOM on the probe input');
{
    const v = new Float32Array(FRAMES * WIDEST);
    for (let i = 0; i < v.length; ++i) v[i] = Math.sin(i * 0.031) * 0.4;
    M.HEAPF32.set(v, ptr >>> 2);            // the view is taken and dropped before any call into the module
}

// asking the price must not disturb a result, in this tier as in the native one
ok(run.hum(ptr, FRAMES, 2, 48000) === 1, 'hum_run accepts the probe programme');
const scalarsBefore = (() => {
    const p = M._malloc(64 * 8);
    const n = M._fc_probe_hum_scalars(p, 64);
    const v = new Float64Array(M.HEAPF64.buffer, p, n).slice();
    M._free(p);
    return Array.from(v).join(',');
})();
for (const m of MODES) { query[m](16, 768000); query[m](0, 48000); query[m](2, -1); }
const scalarsAfter = (() => {
    const p = M._malloc(64 * 8);
    const n = M._fc_probe_hum_scalars(p, 64);
    const v = new Float64Array(M.HEAPF64.buffer, p, n).slice();
    M._free(p);
    return Array.from(v).join(',');
})();
ok(scalarsBefore === scalarsAfter && scalarsBefore.length > 0,
   'hum_scalars reads the same after a round of price queries');

// The refusal sets, in this tier — and over THREE programme lengths, the empty one included. With a
// single non-empty length this loop agreed with itself while a run that ignored what `prepare()` returned
// accepted a geometry priced at zero: `_run` skips `process()` when frames is 0, so an unprepared
// analyzer still reached `finish()`, still set its flag and still answered 1. Measured in the native
// suite and replayed here on a module built from the mutated source, where this check passed 10/10.
// The low bound is the core's 8000 Hz (P51): 8000 - 2**-40 is the double one ulp under it. 44.1 is the unit error
// the floor exists for; 1000 was the old floor (999 the row under it) and 2000 forensics' own — all refusals now.
const RATES = [0, -1, 44.1, 999, 1000, 2000, 7999, 8000 - 2 ** -40, 8000, 8130.5, 18367, 18368, 44100, 48000,
               96000, 192000, 768000, 768000.5, Infinity, -Infinity, NaN];
const WIDTHS = [0, 1, 2, 6, 16, 17, 0xFFFFFFFF];
const COUNTS = [0, 64, FRAMES];
let rows = 0, disagreed = 0, spoke = 0;
const silent = { report: '_fc_probe_report_counts', bursts: '_fc_probe_bursts_scalars',
                 hum: '_fc_probe_hum_scalars', forensics: '_fc_probe_forensics_scalars',
                 lowend: '_fc_probe_lowend_scalars',
                 stereobursts: '_fc_probe_stereobursts_scalars' };
const scratch = M._malloc(64 * 8);
if (!scratch) refuse('wasm OOM on the probe scratch');
for (const m of MODES)
    for (const sr of RATES)
        for (const ch of WIDTHS)
            for (const n of COUNTS) {
                const priceable = query[m](ch, sr) > 0;
                const expected = priceable && (n !== 0 || ACCEPTS_EMPTY[m]);
                const runnable = run[m](ptr, n, ch, sr) === 1;
                ++rows;
                if (runnable !== expected) {
                    ++disagreed;
                    console.error(`FAIL: ${m} priced ${priceable ? '>0' : '0'} but _run ${runnable ? 'accepted' : 'refused'} at ${ch} x ${sr} with ${n} frames`);
                }
                // and a refused run must leave nothing readable behind
                if (!runnable && M[silent[m]](scratch, 60) !== 0) ++spoke;
            }
ok(disagreed === 0, `the price and the run agree about all ${rows} rows in the wasm tier`);
ok(spoke === 0, 'and every refused run of those left its getters silent');
M._free(scratch);

// the canonical zero, and the headline number, IN THIS TIER — wasm32 sizes its own structs
for (const m of MODES) ok(Object.is(query[m](2, 0), 0), `${m}: a refused geometry quotes +0.0, not -0.0 or NaN`);
ok(query.hum(16, 768000) === 352688184, `hum(16, 768000) is 352688184 here too, not ${query.hum(16, 768000)}`);
ok(query.hum(1, 8130.5) === 1110776, `hum(1, 8130.5) is 1110776 here too — a fractional rate survives the boundary`);
for (const m of MODES)
    ok(query[m](2, 8000 - 2 ** -40) === 0 && query[m](2, 44.1) === 0,
       `${m}: one ulp under the 8000 Hz floor and 44.1 are refused here too`);

M._free(ptr);
console.log(`${checks} checks, ${bad} failures`);
process.exit(bad === 0 ? 0 : 1);
