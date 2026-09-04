// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// DSP-ARCHITECTURE.md §2 law 9, enforced instead of merely written down: no `long double` in core code.
// The type has no fixed meaning — measured 8 bytes / 53 mantissa bits on arm64 macOS, 16 / 64 on x86-64,
// 16 / 113 on wasm32 — so one source file gives three different answers, and the two 16-byte answers are
// DIFFERENT FORMATS. The law names one sanctioned exception (`tools/fcore_measure.cpp`'s native-only
// `correlation` accumulator) and the ABI probe that measures the type on purpose; both live under tools/,
// which this lint does not scan, so neither needs an allowlist that could rot into a loophole.
//
// Usage: node tools/lint/check-no-long-double.mjs [--self-test] [paths...]
//        default paths: modules/*/include and modules/*/src
// Exits non-zero if the TYPE is found. Exit 0 and a one-line summary otherwise.
//
// WHY THIS IS NOT A GREP. `grep -rn 'long double'` cannot work here, and the reason is specific rather
// than fastidious: the words appear in PROSE constantly — law 9's own explanation says them repeatedly,
// and any header that documents why it uses `double` will too. The sibling lint in the same CI job
// (`#include <thread>`) is safe as a bare grep because that string cannot occur accidentally in a
// sentence; this one is not. So the file is lexed first — block comments, line comments, string, char
// and raw-string literals are blanked (keeping line numbers) — and only then matched.
//
// WHAT IT MATCHES, and what it deliberately does not:
//   long double / double long   both legal spellings of the type (`double long x;` compiles)
//   1.0L, 1e5L, 0x1.8p3L        an L-suffixed FLOATING literal IS a long double, and is how the type
//                               sneaks in without ever being named. Integer suffixes (100L, 0xFFUL,
//                               100LL) are NOT matched: the rule requires a '.' or a binary/decimal
//                               exponent, which no integer literal has.
//   NOT matched: a type reached through a macro, an alias, or token pasting, and anything a library
//   returns. That would need a Clang AST check filtered to project locations — stronger, and out of
//   proportion to a rule whose realistic violation is someone typing the words.

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, sep } from 'node:path';

const SOURCE_EXT = /\.(h|hpp|hh|inl|c|cc|cpp|cxx)$/;

// Blank out everything a compiler would not read as code, preserving newlines so line numbers survive.
export function stripNonCode (src)
{
    let out = '';
    const blank = (s) => s.replace(/[^\n]/g, ' ');
    for (let i = 0; i < src.length;)
    {
        const two = src.slice(i, i + 2);
        if (two === '/*') { const e = src.indexOf('*/', i + 2); const end = e < 0 ? src.length : e + 2; out += blank(src.slice(i, end)); i = end; continue; }
        if (two === '//') { const e = src.indexOf('\n', i); const end = e < 0 ? src.length : e; out += blank(src.slice(i, end)); i = end; continue; }
        // Raw string: R"delim( ... )delim"  — its body may legally contain anything, comments included.
        const raw = /^(?:u8|u|U|L)?R"([^()\\ \t\n]{0,16})\(/.exec(src.slice(i, i + 24));
        if (raw) { const close = ')' + raw[1] + '"'; const e = src.indexOf(close, i + raw[0].length); const end = e < 0 ? src.length : e + close.length; out += blank(src.slice(i, end)); i = end; continue; }
        if (src[i] === '"' || src[i] === '\'')
        {
            const q = src[i]; let j = i + 1;
            while (j < src.length && src[j] !== q && src[j] !== '\n') { if (src[j] === '\\') j++; j++; }
            const end = Math.min(j + 1, src.length); out += blank(src.slice(i, end)); i = end; continue;
        }
        out += src[i]; i++;
    }
    return out;
}

// Both spellings of the type, plus an L-suffixed floating literal (which IS one without naming it).
const TYPE_RE = /(?<![A-Za-z0-9_])(long\s+double|double\s+long)(?![A-Za-z0-9_])/;
const LSUFFIX_RE = /(?<![A-Za-z0-9_.])(?:\d+\.\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+|0[xX][0-9a-fA-F]*\.?[0-9a-fA-F]*[pP][+-]?\d+)[lL](?![A-Za-z0-9_])/;

export function findViolations (text)
{
    const hits = [];
    stripNonCode(text).split('\n').forEach((line, n) =>
    {
        const t = TYPE_RE.exec(line);      if (t) hits.push({ line: n + 1, what: t[1].replace(/\s+/g, ' ') });
        const l = LSUFFIX_RE.exec(line);   if (l) hits.push({ line: n + 1, what: l[0] + '  (long double literal)' });
    });
    return hits;
}

function walk (dir, acc)
{
    for (const e of readdirSync(dir))
    {
        const p = join(dir, e);
        if (statSync(p).isDirectory()) walk(p, acc);
        else if (SOURCE_EXT.test(e)) acc.push(p);
    }
    return acc;
}

function defaultRoots ()
{
    const roots = [];
    for (const m of readdirSync('modules'))
        for (const sub of ['include', 'src'])
        {
            const p = join('modules', m, sub);
            try { if (statSync(p).isDirectory()) roots.push(p); } catch { /* module has no such dir */ }
        }
    return roots;
}

// The matcher's own negative control: prose must pass, every spelling of the type must fail. A lint that
// silently stopped matching would pass its CI step forever; this is what notices.
function selfTest ()
{
    const cases = [
        ['// a comment that says long double must not trip it',            0],
        ['const char* s = "long double";',                                 0],
        ['auto x = R"(long double)";',                                     0],
        ['/* long double\n   across lines */ int ok = 1;',                 0],
        ['int n = 100L; unsigned u = 0xFFUL; long long b = 100LL;',        0],
        ['double d = 1.0; float f = 1.0f;',                                0],
        ['long double x = 0.0;',                                           1],
        ['double long x = 0.0;',                                           1],
        ['auto x = 1.0L;',                                                 1],
        ['auto x = 1e5L;',                                                 1],
        ['auto x = 0x1.8p3L;',                                             1],
        ['long   double   spaced;',                                        1],
    ];
    let bad = 0;
    for (const [src, want] of cases)
    {
        const got = findViolations(src).length ? 1 : 0;
        if (got !== want) { console.error(`  SELF-TEST FAIL: wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    if (bad) { console.error(`long-double lint self-test: ${bad}/${cases.length} cases wrong`); process.exit(1); }
    console.log(`long-double lint self-test: ${cases.length}/${cases.length} cases correct`);
}

const args = process.argv.slice(2);
if (args.includes('--self-test')) { selfTest(); process.exit(0); }

const roots = args.length ? args : defaultRoots();
const files = [];
for (const r of roots) { try { statSync(r).isDirectory() ? walk(r, files) : files.push(r); } catch { console.error(`no such path: ${r}`); process.exit(2); } }

let found = 0;
for (const f of files)
    for (const h of findViolations(readFileSync(f, 'utf8')))
    { console.error(`${f}:${h.line}: long double in core code — DSP-ARCHITECTURE.md §2 law 9 (use double, or Kahan/Neumaier summation)`); console.error(`    ${h.what}`); found++; }

if (found) { console.error(`\n^^ ${found} violation(s) of law 9 across ${files.length} files in ${roots.map(r => r.split(sep).join('/')).join(', ')}`); process.exit(1); }
console.log(`no long double in ${files.length} core source files (${roots.length} roots) — law 9 holds`);
