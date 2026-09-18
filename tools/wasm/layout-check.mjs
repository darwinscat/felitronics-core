// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// fc-master-layout.mjs against the COMPILER, field by field:
//
//   node layout-check.mjs <fcore_master | fcore_master.js>
//
// `fcore_master layout` prints the offset of every field of every struct the ABI carries, as the compiler laid it
// out — natively, or as wasm32 laid it out when the argument is the emscripten build, which is the tier a page
// actually runs. This holds the JavaScript layout against it in BOTH directions: a field whose offset differs, a
// field the JS lists and the compiler has not got, and a field the compiler has and the JS does not list. A size
// check alone cannot see two fields of one type swapped, and that is a page writing its value into the wrong knob.
//
// It also holds the version and the size table: this file's FC_MASTER_ABI_VERSION must be the build's, and every
// struct with a header must be the size the build's own table publishes for it.
//
// AND IT HOLDS `FC_EQ_AXIS` AGAINST THE HEADER'S `fc_eq_axis`, which no offset can reach: a name list is not a
// layout, and its index IS the code a page passes as `lane`. The enumerators are read out of
// tools/fc_master_abi.h and compared entry by entry — declaration order, value and spelling — so a permutation
// there fails here instead of renaming every lane on the page.

import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { resolve, dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { layoutOf, structNames, STRUCT_IDS, FC_MASTER_ABI_VERSION, FC_EQ_AXIS } from './fc-master-layout.mjs';

const [, , bin] = process.argv;
if (!bin) { console.error('usage: node layout-check.mjs <fcore_master | fcore_master.js>'); process.exit(2); }

const r = bin.endsWith('.js')
    ? spawnSync(process.execPath, [resolve(bin), 'layout'], { encoding: 'utf8' })
    : spawnSync(resolve(bin), ['layout'], { encoding: 'utf8' });
if (r.status !== 0) { console.error(`fcore_master layout: exit ${r.status}\n${r.stderr}`); process.exit(1); }

let version = -1;
const sizes = new Map(), fields = new Map(), table = new Map();
for (const line of r.stdout.split('\n')) {
    const t = line.trim().split(/\s+/);
    if (t[0] === 'V') version = Number(t[1]);
    else if (t[0] === 'S') sizes.set(t[1], Number(t[2]));
    else if (t[0] === 'F') { if (!fields.has(t[1])) fields.set(t[1], new Map()); fields.get(t[1]).set(t[2], Number(t[3])); }
    else if (t[0] === 'T') table.set(t[1], { id: Number(t[2]), size: Number(t[3]) });
}

let failures = 0;
const check = (pass, what) => { if (!pass) { console.log(`  FAIL: ${what}`); ++failures; } };

check(version === FC_MASTER_ABI_VERSION, `ABI version: the build speaks v${version}, this file v${FC_MASTER_ABI_VERSION}`);

let compared = 0;
for (const name of structNames()) {
    const js = layoutOf(name), c = fields.get(name);
    check(c !== undefined, `${name}: described here and absent from the build's field list`);
    if (!c) continue;
    check(sizes.get(name) === js.size, `${name}: ${js.size} B here, ${sizes.get(name)} B in the build`);
    for (const [field, f] of js.fields) {
        ++compared;
        check(c.has(field), `${name}.${field}: listed here and not in the build`);
        if (c.has(field)) check(c.get(field) === f.offset, `${name}.${field}: offset ${f.offset} here, ${c.get(field)} in the build`);
    }
    for (const field of c.keys()) check(js.fields.has(field), `${name}.${field}: in the build and missing here`);
}
for (const name of fields.keys()) check(structNames().includes(name), `${name}: in the build and not described here`);

for (const [name, id] of Object.entries(STRUCT_IDS)) {
    const row = table.get(name);
    check(row !== undefined && row.id === id, `${name}: struct id ${id} here, ${row && row.id} in the build`);
    check(row !== undefined && row.size === layoutOf(name).size,
          `${name}: the build's size table says ${row && row.size} B at v${version}, this file ${layoutOf(name).size} B`);
}
for (const name of table.keys()) check(STRUCT_IDS[name] !== undefined, `${name}: has a header in the build and no id here`);

// ── FC_EQ_AXIS against the header's enum ──────────────────────────────────────────────────────────
// The block is REQUIRED to be found: a regex that matched nothing would compare an empty list and pass.
const header = readFileSync(join(dirname(fileURLToPath(import.meta.url)), '..', 'fc_master_abi.h'), 'utf8');
const block = /typedef\s+enum\s+fc_eq_axis\s*\{([^}]*)\}/.exec(header);
check(block !== null, 'fc_eq_axis: no such enum in tools/fc_master_abi.h');
if (block) {
    const codes = [...block[1].matchAll(/FC_EQ_AXIS_([A-Z0-9_]+)\s*=\s*(\d+)/g)]
        .map(m => ({ suffix: m[1], value: Number(m[2]) }));
    check(codes.length === FC_EQ_AXIS.length,
          `fc_eq_axis: ${codes.length} codes in the header, ${FC_EQ_AXIS.length} names here`);
    codes.forEach((c, i) => {
        check(c.value === i, `FC_EQ_AXIS_${c.suffix}: declared at index ${i} in the header and numbered ${c.value}`);
        check(FC_EQ_AXIS[i] !== undefined && FC_EQ_AXIS[i].toUpperCase() === c.suffix,
              `fc_eq_axis[${i}]: the header says ${c.suffix}, this file says ${FC_EQ_AXIS[i]}`);
    });
}

console.log(`layout-check: ${structNames().length} structs, ${compared} fields, `
          + `${FC_EQ_AXIS.length} axis codes, ${failures} failure(s)`);
process.exit(failures === 0 ? 0 : 1);
