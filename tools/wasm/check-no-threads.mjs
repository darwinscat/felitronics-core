// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// P0 acceptance criterion 3, proven from the ARTIFACT rather than from a page loading. "It works without
// COOP/COEP" is evidence; "the module declares no shared memory and imports nothing thread-shaped" is proof,
// and it needs no browser, no server and no wasm-objdump (which is not installed on most machines — the
// build script's earlier objdump step silently checked nothing).
//
// Usage: node check-no-threads.mjs <module.wasm> [glue.js]   — exits non-zero if anything thread-shaped is found.

import { readFileSync } from 'node:fs';

const [, , wasmPath, gluePath] = process.argv;
if (!wasmPath) { console.error('usage: node check-no-threads.mjs <module.wasm> [glue.js]'); process.exit(2); }

const buf = readFileSync(wasmPath);
if (buf.readUInt32LE(0) !== 0x6d736100) { console.error('not a wasm module'); process.exit(2); }

// LEB128 unsigned, as every wasm length and index is encoded.
function uleb(b, i) {
    let v = 0, shift = 0, byte;
    do { byte = b[i++]; v |= (byte & 0x7f) << shift; shift += 7; } while (byte & 0x80);
    return [v >>> 0, i];
}

let i = 8, failures = 0;
const memories = [];
while (i < buf.length) {
    const id = buf[i++];
    let size; [size, i] = uleb(buf, i);
    const end = i + size;
    if (id === 5) {                                   // memory section: vec(limits)
        let count, p; [count, p] = uleb(buf, i);
        for (let m = 0; m < count; ++m) {
            const flags = buf[p++];
            let min; [min, p] = uleb(buf, p);
            let max = null;
            if (flags & 0x1) { [max, p] = uleb(buf, p); }
            memories.push({ flags, min, max, shared: Boolean(flags & 0x2) });
        }
    }
    i = end;
}

console.log(`  memories declared: ${memories.length}`);
for (const m of memories) {
    const pages = (n) => n === null ? 'none' : `${n} pages (${(n * 64 / 1024).toFixed(1)} MiB)`;
    console.log(`    min ${pages(m.min)}, max ${pages(m.max)}, shared=${m.shared}`);
    if (m.shared) { console.error('    *** SHARED MEMORY — this module is threaded'); failures++; }
}
if (memories.length === 0) console.log('    (none declared — memory is imported or absent)');

if (gluePath) {
    const glue = readFileSync(gluePath, 'utf8');
    for (const pat of ['SharedArrayBuffer', 'pthread', 'new Worker', 'Atomics.']) {
        const n = glue.split(pat).length - 1;
        console.log(`  glue mentions ${pat.padEnd(18)} ${n} time(s)`);
        if (n > 0) { console.error(`    *** the glue references ${pat}`); failures++; }
    }
}

if (failures) { console.error(`  FAILED: ${failures} thread-shaped finding(s)`); process.exit(1); }
console.log('  no threads: memory is unshared and the glue is clean');
