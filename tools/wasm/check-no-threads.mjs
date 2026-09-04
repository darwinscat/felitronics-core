// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// P0 acceptance criterion 3, proven from the ARTIFACT rather than from a page loading. "It works without
// COOP/COEP" is evidence; "the module declares no shared memory and imports nothing thread-shaped" is proof,
// and it needs no browser, no server and no wasm-objdump (which is not installed on most machines — the
// build script's earlier objdump step silently checked nothing).
//
// Usage: node check-no-threads.mjs <module.wasm> [glue.js]   — one module (what tools/wasm/build.sh does)
//        node check-no-threads.mjs --all <dir>               — every .wasm under dir, each paired with its
//                                                              sibling .js; one line per module + a summary
// Exits non-zero if anything thread-shaped is found.
//
// WHY THE IMPORT SECTION IS PARSED. The first version read only the memory section (id 5) and printed
// "(none declared — memory is imported or absent)" otherwise. Measured on emsdk 6.0.9: a `-pthread` build
// declares NO memory section at all — the shared memory is IMPORTED from JS — so the one check that was
// supposed to catch a threaded module returned "memories declared: 0" and passed. The memory test was
// vacuous in exactly the case it existed for; only the glue string scan caught it. Imports are parsed now,
// so the proof holds on the wasm alone, with or without glue.

import { readFileSync, readdirSync, statSync, lstatSync } from 'node:fs';
import { join, extname } from 'node:path';

// LEB128 unsigned, as every wasm length and index is encoded.
function uleb (b, i) {
    let v = 0, shift = 0, byte;
    do { byte = b[i++]; v += (byte & 0x7f) * 2 ** shift; shift += 7; } while (byte & 0x80);
    return [v, i];
}

// limits := flags min [max].  bit0 = has-max, bit1 = SHARED, bit2 = 64-bit index (memory64).
function limits (b, p) {
    const flags = b[p++];
    let min, max = null;
    [min, p] = uleb(b, p);
    if (flags & 0x1) [max, p] = uleb(b, p);
    return [{ flags, min, max, shared: Boolean(flags & 0x2) }, p];
}

function auditWasm (path) {
    const buf = readFileSync(path);
    if (buf.length < 8 || buf.readUInt32LE(0) !== 0x6d736100) return { error: 'not a wasm module' };
    // VALIDATE BEFORE PARSING. Without this the hand walk below can be fed a module that is not wasm at all
    // and report it CLEAN: `0061736d01000000050101` says "memory section, count 1" and then stops, so the
    // descriptor read runs off the end, `undefined & flags` invents 0, and the checker cheerfully announces
    // an unshared memory. A checker that green-lights garbage is worse than none. WebAssembly.validate is
    // already in node, rejects exactly that module, and also covers the version word and every truncation.
    if (!WebAssembly.validate(buf)) return { error: 'not a VALID wasm module (WebAssembly.validate refused it)' };

    const memories = [], threadImports = [];
    let unknownKind = null;
    let i = 8;
    while (i < buf.length) {
        const id = buf[i++];
        let size; [size, i] = uleb(buf, i);
        const end = i + size;
        if (id === 2) {                                   // import section: vec(mod, field, desc)
            let count, p; [count, p] = uleb(buf, i);
            for (let k = 0; k < count && p < end; ++k) {
                let n;
                [n, p] = uleb(buf, p); const mod = buf.toString('utf8', p, p + n); p += n;
                [n, p] = uleb(buf, p); const fld = buf.toString('utf8', p, p + n); p += n;
                const kind = buf[p++];
                if      (kind === 0x00) { let t; [t, p] = uleb(buf, p); }              // func: typeidx
                else if (kind === 0x01) { p++; let l; [l, p] = limits(buf, p); }       // table: reftype + limits
                else if (kind === 0x02) { let m; [m, p] = limits(buf, p);              // memory: limits
                                          memories.push({ ...m, from: `${mod}.${fld}` }); }
                else if (kind === 0x03) { p += 2; }                                    // global: valtype + mut
                else if (kind === 0x04) { p++; let t; [t, p] = uleb(buf, p); }         // tag: attribute + typeidx
                else {
                    // An unknown descriptor means the rest of this section cannot be decoded — and the old
                    // `break` here was a SILENT PASS: a tag import (0x04, what -fwasm-exceptions emits) hid
                    // every import after it, so a shared memory or a pthread import following one was never
                    // seen and the module was reported clean. Refuse instead, exactly as an absent memory is
                    // refused: this file's whole claim is proof from the artifact, and a parser that gives up
                    // quietly proves nothing.
                    unknownKind = `0x${kind?.toString(16) ?? '??'} at import ${k} (${mod}.${fld})`;
                    break;
                }
                if (p > end) { unknownKind = `import section overran its own size (${p} > ${end})`; break; }
                if (/pthread|_emscripten_(init|thread)_|thread_spawn/i.test(`${mod}.${fld}`)) threadImports.push(`${mod}.${fld}`);
            }
        } else if (id === 5) {                            // memory section: vec(limits)
            let count, p; [count, p] = uleb(buf, i);
            for (let m = 0; m < count && p < end; ++m) { let l; [l, p] = limits(buf, p); memories.push({ ...l, from: null }); }
        }
        i = end;
    }
    return { memories, threadImports, unknownKind };
}

// Tokens that actually discriminate. A bare 'pthread' does NOT: an -sASSERTIONS build emits
// assert(... "Module.pthreadMainPrefixURL option was removed" ...), a deprecation string with no thread
// support behind it, so every one of the 78 tier artifacts "failed" on it. Measured counts, threaded build
// vs a tier test vs the P0 web artifact: PThread 31/0/0, ENVIRONMENT_IS_PTHREAD 20/0/0, _emscripten_thread
// 22/0/0, Atomics. 5/0/0, pthread_create 3/0/0, SharedArrayBuffer 1/0/0, new Worker 1/0/0.
const GLUE_PATTERNS = ['SharedArrayBuffer', 'PThread', 'ENVIRONMENT_IS_PTHREAD', 'pthread_create',
                       '_emscripten_thread', 'new Worker', 'Atomics.'];

function auditGlue (path) {
    const glue = readFileSync(path, 'utf8');
    return GLUE_PATTERNS.map(pat => ({ pat, n: glue.split(pat).length - 1 }));
}

function reportOne (wasmPath, gluePath, verbose) {
    const { error, memories, threadImports, unknownKind } = auditWasm(wasmPath);
    if (error) { console.error(`  ${wasmPath}: ${error}`); return 1; }
    let failures = 0;
    if (unknownKind) { console.error(`    *** cannot decode the import section: ${unknownKind}`); failures++; }

    if (verbose) console.log(`  memories: ${memories.length}`);
    for (const m of memories) {
        const pages = n => n === null ? 'none' : `${n} pages (${(n * 64 / 1024).toFixed(1)} MiB)`;
        const where = m.from ? `imported ${m.from}` : 'declared';
        if (verbose) console.log(`    ${where}: min ${pages(m.min)}, max ${pages(m.max)}, shared=${m.shared}`);
        if (m.shared) { console.error(`    *** SHARED MEMORY (${where}) — this module is threaded`); failures++; }
    }
    // A module with neither a declared nor an imported memory is not "clean", it is unreadable to this
    // check — say so rather than passing it silently, which is the bug this rewrite exists to remove.
    if (memories.length === 0) { console.error(`    *** no memory declared OR imported in ${wasmPath} — cannot prove anything`); failures++; }

    for (const name of threadImports) { console.error(`    *** thread-shaped import: ${name}`); failures++; }

    if (gluePath) {
        for (const { pat, n } of auditGlue(gluePath)) {
            if (verbose) console.log(`  glue mentions ${pat.padEnd(18)} ${n} time(s)`);
            if (n > 0) { console.error(`    *** the glue references ${pat}`); failures++; }
        }
    }
    return failures;
}

const args = process.argv.slice(2);
if (args[0] === '--all') {
    const root = args[1];
    if (!root) { console.error('usage: node check-no-threads.mjs --all <dir>'); process.exit(2); }
    const found = [];
    (function walk (d) {
        for (const e of readdirSync(d)) {
            const p = join(d, e);
            if (lstatSync(p).isSymbolicLink()) continue;          // never follow a link out of the tree, or in a circle
            if (statSync(p).isDirectory()) walk(p);
            else if (extname(p) === '.wasm') found.push(p);
        }
    })(root);
    if (found.length === 0) { console.error(`no .wasm under ${root} — nothing was audited, which is not a pass`); process.exit(2); }
    let bad = 0;
    for (const w of found.sort()) {
        const glue = w.replace(/\.wasm$/, '.js');
        let g = null; try { statSync(glue); g = glue; } catch { /* built without glue */ }
        const f = reportOne(w, g, false);
        console.log(`  ${f === 0 ? 'ok  ' : 'FAIL'}  ${w}`);
        if (f) bad++;
    }
    console.log(`\n  ${found.length} module(s) audited, ${bad} thread-shaped`);
    process.exit(bad ? 1 : 0);
}

const [wasmPath, gluePath] = args;
if (!wasmPath) { console.error('usage: node check-no-threads.mjs <module.wasm> [glue.js]  |  --all <dir>'); process.exit(2); }
const failures = reportOne(wasmPath, gluePath, true);
if (failures) { console.error(`  FAILED: ${failures} thread-shaped finding(s)`); process.exit(1); }
console.log('  no threads: memory is unshared and the glue is clean');
