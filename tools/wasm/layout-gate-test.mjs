// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// The two refusals fc-master-layout.mjs makes on a page's behalf, tested without a module:
//
//   node layout-gate-test.mjs
//
// * `assertLayoutMatches` — a module at this file's version or NEWER loads, an older one does not, and a module
//   whose size table disagrees with this file about any struct at this file's version does not either.
// * `Struct.set` — refuses a field on a struct whose header is not stamped at this file's version and size. That
//   is exactly a struct filled by a frozen v1 `_fc_*_default` writer and then given a newer field, which the
//   module would never read (VERSIONING rules 6 and 8 in tools/fc_master_abi.h).
//
// Everything here is JavaScript logic over a fake heap: whether the real module's table agrees with this file is
// layout-check.mjs's question (against the compiler) and master-parity.mjs's (against the running module).

import { Struct, sizeOf, assertLayoutMatches, STRUCT_IDS, FC_MASTER_ABI_VERSION } from './fc-master-layout.mjs';

let failures = 0, checks = 0;
const check = (pass, what) => { ++checks; if (!pass) { console.log(`  FAIL: ${what}`); ++failures; } };
const throws = (fn) => { try { fn(); return false; } catch { return true; } };

const V = FC_MASTER_ABI_VERSION;
const moduleAt = (version, sizeAt = (name) => sizeOf(name)) => ({
    _fc_master_abi_version: () => version,
    _fc_master_sizeof: (id, v) => {
        const name = Object.keys(STRUCT_IDS).find(n => STRUCT_IDS[n] === id);
        return (name && v >= 1 && v <= version) ? sizeAt(name, v) : 0;
    },
});

// ── the load gate ─────────────────────────────────────────────────────────────────────────────────
check(!throws(() => assertLayoutMatches(moduleAt(V))), 'a module at this file\'s version loads');
check(!throws(() => assertLayoutMatches(moduleAt(V + 1))), 'a NEWER module loads — it reads this file\'s structs');
check(throws(() => assertLayoutMatches(moduleAt(V - 1))), 'an OLDER module is refused');
for (const name of Object.keys(STRUCT_IDS))
    check(throws(() => assertLayoutMatches(moduleAt(V, (n) => sizeOf(n) + (n === name ? 8 : 0)))),
          `a module whose table disagrees about ${name} is refused`);
let askedAt = -1;
assertLayoutMatches({ _fc_master_abi_version: () => V + 5,
                      _fc_master_sizeof: (id, v) => { askedAt = v; const n = Object.keys(STRUCT_IDS).find(k => STRUCT_IDS[k] === id); return sizeOf(n); } });
check(askedAt === V, `the sizes are asked at THIS file's version (v${V}), not the module's (asked v${askedAt})`);

// ── the write gate ────────────────────────────────────────────────────────────────────────────────
const heap = { HEAPF32: new Float32Array(new ArrayBuffer(1 << 16)) };
const view = () => new DataView(heap.HEAPF32.buffer);
const ptr = 64;

const cfg = new Struct(heap, 'fc_master_config', ptr);
view().setUint32(ptr, 1, true); view().setUint32(ptr + 4, 80, true);          // what a frozen v1 writer stamps
check(throws(() => cfg.set('deliveryRate', 96000)), 'a v1-stamped config refuses `deliveryRate` — the module would not read it');
check(throws(() => cfg.set('sampleRate', 48000)), 'and refuses every other field too: the whole struct is at the wrong version');
check(!throws(() => cfg.set('header.structSize', 88)), 'while its header stays writable, so a caller can stamp it');
cfg.init();
check(!throws(() => cfg.set('deliveryRate', 96000)) && cfg.get('deliveryRate') === 96000,
      'stamped by init(), the field is written and reads back');

const prm = new Struct(heap, 'fc_master_params', ptr);
view().setUint32(ptr, V, true); view().setUint32(ptr + 4, sizeOf('fc_master_params') - 8, true);
check(throws(() => prm.set('compressorMix', 0.5)), 'this file\'s version with an older SIZE is refused as well');
prm.init();
check(!throws(() => prm.set('compressorMix', 0.5)) && prm.get('compressorMix') === 0.5, 'and accepted once stamped');

const res = new Struct(heap, 'fc_master_resolved', ptr);
view().setUint32(ptr, 1, true); view().setUint32(ptr + 4, 80, true);          // an OUT struct stamped v1
check(throws(() => res.get('compressorMix')), 'get() is gated too: a v1-stamped resolved does not read a field the module never wrote');
res.init();
check(!throws(() => res.get('compressorMix')), 'and reads once stamped at this file\'s version');

const pass = new Struct(heap, 'fc_solve_pass', ptr);                             // no header: nothing to check
check(!throws(() => pass.set('gainDb', 1.5)), 'a struct without a header is not gated');

console.log(`layout-gate-test: ${checks} checks, ${failures} failure(s)`);
process.exit(failures === 0 ? 0 : 1);
