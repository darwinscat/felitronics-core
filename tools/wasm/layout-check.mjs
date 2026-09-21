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
// AND IT HOLDS EVERY ENUM LIST OF `FC_ENUMS` AGAINST THE HEADER'S OWN ENUM, which no offset can reach: a name
// list is not a layout, and its index IS the code a page writes into a struct or passes to an entry point. The
// enumerators are read out of tools/fc_master_abi.h and compared entry by entry — declaration order, value and
// letters, in both directions. The word breaks are the JavaScript's own (`LowShelf` for `FC_FILTER_LOW_SHELF`),
// so the comparison drops the underscores; nothing else about the spelling is free.
//
// AND IT HOLDS `FC_DOMAINS` AGAINST THE LAYOUTS — not the domain VALUES, which are behaviour and belong to
// tools/tests/MasterDomainsTests.cpp, but that every row names a value field described here, that its bound
// spellings parse, that an `enum:` unit names a list of FC_ENUMS and that `resolved` names a field of
// `fc_master_resolved`.

import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { resolve, dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { layoutOf, structNames, STRUCT_IDS, FC_MASTER_ABI_VERSION,
         FC_ENUMS, FC_DOMAINS, domainBound } from './fc-master-layout.mjs';

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

// ── FC_ENUMS against the header's enums ───────────────────────────────────────────────────────────
// Each block is REQUIRED to be found: a regex that matched nothing would compare an empty list and pass.
const header = readFileSync(join(dirname(fileURLToPath(import.meta.url)), '..', 'fc_master_abi.h'), 'utf8');
let codesCompared = 0;
for (const [jsName, spec] of Object.entries(FC_ENUMS)) {
    const block = new RegExp(`typedef\\s+enum\\s+${spec.enum}\\s*\\{([^}]*)\\}`).exec(header);
    check(block !== null, `${spec.enum}: no such enum in tools/fc_master_abi.h`);
    if (!block) continue;
    const codes = [...block[1].matchAll(new RegExp(`${spec.prefix}([A-Z0-9_]+)\\s*=\\s*(\\d+)`, 'g'))]
        .map(m => ({ suffix: m[1], value: Number(m[2]) }));
    check(codes.length === spec.names.length,
          `${spec.enum}: ${codes.length} codes in the header, ${spec.names.length} names in ${jsName}`);
    codes.forEach((c, i) => {
        ++codesCompared;
        check(c.value === i, `${spec.prefix}${c.suffix}: declared at index ${i} in the header and numbered ${c.value}`);
        check(spec.names[i] !== undefined && spec.names[i].toUpperCase() === c.suffix.replace(/_/g, ''),
              `${jsName}[${i}]: the header says ${spec.prefix}${c.suffix}, this file says ${spec.names[i]}`);
    });
}

// ── FC_DOMAINS against the struct layouts ─────────────────────────────────────────────────────────
const fieldExists = (path) => {
    const parts = path.split('.');
    let layout;
    try { layout = layoutOf(parts[0]); } catch { return false; }
    for (let i = 1; i < parts.length; ++i) {
        const raw = parts[i], step = raw.replace(/\[\]$/, '');
        const f = layout.fields.get(step);
        if (!f) return false;
        if (raw.endsWith('[]') !== (f.count > 0)) return false;    // `[]` where there is no array, or the reverse
        let inner = null;
        try { inner = layoutOf(f.type); } catch { inner = null; }  // a scalar type has no layout to walk into
        if (inner === null) return i === parts.length - 1;         // a row must end on a VALUE, never on a struct
        layout = inner;
    }
    return false;
};
// Every scalar leaf of the three INPUT structs, as FC_DOMAINS spells it: `[]` for an array, the `header` and
// the named padding skipped. A field with no row is a field a page has no domain for.
const inputLeaves = (name) => {
    const out = [];
    const walk = (struct, prefix) => {
        for (const [field, f] of layoutOf(struct).fields) {
            if (field === 'header' || field.startsWith('_pad')) continue;
            const step = prefix + field + (f.count > 0 ? '[]' : '');
            let inner = null;
            try { inner = layoutOf(f.type); } catch { inner = null; }
            if (inner === null) out.push(step); else walk(f.type, step + '.');
        }
    };
    walk(name, name + '.');
    return out;
};

const described = new Set(FC_DOMAINS.map(r => r.field));
let leaves = 0;
for (const struct of ['fc_master_config', 'fc_master_params', 'fc_loudness_request'])
    for (const leaf of inputLeaves(struct)) {
        ++leaves;
        check(described.has(leaf), `FC_DOMAINS: ${leaf} is an input field with no domain row`);
    }

const seen = new Set();
for (const row of FC_DOMAINS) {
    check(!seen.has(row.field), `FC_DOMAINS: ${row.field} is listed twice`);
    seen.add(row.field);
    check(fieldExists(row.field), `FC_DOMAINS: ${row.field} is not a value field of any struct described here`);
    for (const which of ['min', 'max'])
        try { domainBound(row[which], 48000); } catch (e) { check(false, `FC_DOMAINS: ${row.field}.${which}: ${e.message}`); }
    if (row.unit.startsWith('enum:'))
        check(FC_ENUMS[row.unit.slice(5)] !== undefined, `FC_DOMAINS: ${row.field} names ${row.unit}, which is not in FC_ENUMS`);
    if (row.resolved !== '' && row.resolved !== 'eqCurve' && row.resolved !== 'render')
        check(fieldExists(`fc_master_resolved.${row.resolved}`),
              `FC_DOMAINS: ${row.field} reads back through fc_master_resolved.${row.resolved}, which does not exist`);
}

console.log(`layout-check: ${structNames().length} structs, ${compared} fields, `
          + `${codesCompared} enum codes, ${FC_DOMAINS.length} domain rows over ${leaves} input fields, `
          + `${failures} failure(s)`);
process.exit(failures === 0 ? 0 : 1);
