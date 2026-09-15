// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// P52, enforced instead of merely fixed: ONE allocation counter, in test_support/alloc_counter.h.
//
// WHAT WENT WRONG THAT THIS PREVENTS. A suite proves "process() does not allocate" by replacing the global
// allocation functions and reading a delta. The set of functions to replace is EIGHT, not two: C++17 routes
// any object whose `alignof` exceeds __STDCPP_DEFAULT_NEW_ALIGNMENT__ through
// `operator new(std::size_t, std::align_val_t)`, a different function, and the tree is full of such objects
// (`core::SeamAllocator<64>`, `alignof(eq::EqEngine) == 64`). On `a9816e2` the idiom was copied into 61
// translation units; 11 of the copies installed the over-aligned form and 50 did not. Nothing told them
// apart, because nothing was looking — the same sentence meant two different measurements depending on
// which copy a suite happened to start from.
//
// So the copies are gone and there is one header. This lint is what keeps it that way: a new suite that
// hand-rolls a counter again gets a red CI step naming the line, instead of a silent half-measurement that
// looks exactly like the real thing.
//
// The RUNTIME half of the gate lives in the header itself and does not depend on this file: before main()
// it asks the counter for one allocation through each of the eight forms and aborts the run, naming them,
// if any went uncounted. This one catches what a runtime probe structurally cannot — a counter that is
// complete today and quietly loses a form tomorrow, in a file that never includes the shared header at all.
//
// Usage: node tools/lint/check-alloc-counter.mjs [--self-test] [paths...]
//        default paths: modules, tools, test_support
// Exits non-zero if a replacement is defined outside the shared header.
//
// WHY THIS IS NOT A BARE GREP, same reason as its sibling: the words "operator new" are all over the prose
// in this tree — the header above says them a dozen times, and so does every comment explaining why a
// budget is counted. So the file is lexed first (comments, strings, char and raw-string literals blanked,
// line numbers preserved) and only then matched.
//
// WHAT IT MATCHES: a definition or declaration of `operator new` / `operator delete`, in any of their
// twenty forms, **at global scope** — which is the only place a replacement can live. Scope is decided by
// brace depth over the stripped text, not by indentation, and the signature may be split across lines
// (`void*` on one, `operator new(...)` on the next) because a line-at-a-time matcher misses exactly that.
// WHAT IT DELIBERATELY DOES NOT MATCH: a CLASS-scope `operator new` (inside a `struct`/`class` body, or
// spelled `T::operator new`), which is not a replacement and is nobody's business here; a member or local
// declaration at any depth; a CALL to `::operator new(n)`, which is how a test legitimately reaches the
// allocator; and a form inside a named namespace, which likewise replaces nothing.
// WHAT IT MATCHES THAT IS NOT STRICTLY A REPLACEMENT: a global PLACEMENT form such as
// `void* operator new (std::size_t, Arena&)`. That is deliberate rather than sloppy — this tree has none,
// and a test that wants one is defining a global allocation function of its own, which is the conversation
// this lint exists to force.

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, sep } from 'node:path';

const SOURCE_EXT = /\.(h|hpp|hh|inl|c|cc|cpp|cxx)$/;

// The one file allowed to define them — the shared counter itself.
const ALLOWED = /(^|[\\/])test_support[\\/]alloc_counter\.h$/;

// Blank out everything a compiler would not read as code, preserving newlines so line numbers survive.
// (Kept in step with tools/lint/check-no-long-double.mjs, which needs the same lexer for the same reason.)
export function stripNonCode (src)
{
    let out = '';
    const blank = (s) => s.replace(/[^\n]/g, ' ');
    for (let i = 0; i < src.length;)
    {
        const two = src.slice(i, i + 2);
        if (two === '/*') { const e = src.indexOf('*/', i + 2); const end = e < 0 ? src.length : e + 2; out += blank(src.slice(i, end)); i = end; continue; }
        if (two === '//') { const e = src.indexOf('\n', i); const end = e < 0 ? src.length : e; out += blank(src.slice(i, end)); i = end; continue; }
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

// ANCHORED ON THE OPERATOR, NOT ON THE RETURN TYPE. Matching `void* ... operator new` reads naturally and
// is wrong: two crew seats broke it in a minute with shapes that all compile —
//   void*\noperator new (std::size_t)          the return type on its own line (clang-format does this)
//   [[nodiscard]] void* operator new (...)      an attribute in front
//   __attribute__((malloc)) void* operator new  likewise
//   extern "C++" void* operator new (...)       a linkage specifier
//   auto operator new (std::size_t) -> void*    a trailing return type, with no `void` in front at all
// Every one of those is a working replacement that a return-type matcher lets through, which is the silent
// half-measurement the lint exists to prevent. So the anchor is `operator new` / `operator delete` itself,
// at global brace depth, and what is REJECTED is the short list of things that can precede it and make it
// something other than a declaration.
const OPERATOR_RE = /\boperator\s+(new|delete)\s*(\[\s*\])?\s*\(/g;

// What may precede a call or a qualified name rather than a declaration: `::operator new (n)` (the legal way
// a test reaches the allocator), `T::operator new`, `= operator new(8)`, `return operator new(8)`, and an
// operator name appearing as a function argument.
const NOT_A_DECLARATION = /(?:::|=|\(|,|\.|->|\breturn|\bnew|\bdelete)\s*$/;

export function findViolations (text)
{
    const code = stripNonCode(text);

    // Brace depth at every offset, so scope is decided by the code rather than by how it is indented. A
    // replacement is only a replacement at global scope; anything nested is a member, a local, or lives in a
    // namespace, and none of those replace the global allocation functions.
    const depth = new Int32Array(code.length);
    let d = 0;
    for (let i = 0; i < code.length; i++) { const c = code[i]; if (c === '}') d--; depth[i] = d; if (c === '{') d++; }

    const hits = [];
    OPERATOR_RE.lastIndex = 0;
    for (let m; (m = OPERATOR_RE.exec(code)) !== null;)
    {
        if (depth[m.index] !== 0) continue;                              // a member, a local, or in a namespace
        if (NOT_A_DECLARATION.test(code.slice(0, m.index))) continue;    // a call, or a qualified name
        const line = 1 + (code.slice(0, m.index).match(/\n/g) || []).length;
        hits.push({ line, what: `global operator ${m[1]}${m[2] ? '[]' : ''}` });
    }
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

// The matcher's own negative control: every historical shape of the hand-rolled counter must fail, and the
// things that merely look like one must not. A lint that silently stopped matching would pass its CI step
// for ever; this is what notices. The failing cases below are copied from the 61 copies P52 removed.
function selfTest ()
{
    const cases = [
        ['// a comment about operator new being counted',                                          0],
        ['const char* s = "void* operator new (std::size_t)";',                                    0],
        ['auto x = R"(void* operator new (std::size_t s))";',                                      0],
        ['/* void* operator new\n   across lines */ int ok = 1;',                                  0],
        ['void* raw = ::operator new (n);',                                                        0],
        ['    void* p = ::operator new (n, std::align_val_t (64));',                               0],
        ['::operator delete (p, std::align_val_t (64));',                                          0],
        ['static void* countedAlignedNew (std::size_t s, std::size_t a) { return nullptr; }',      0],
        ['struct T { static void* operator_new_like(); };',                                        0],
        // The nineteenth case and its neighbours, from the crew round: a class-specific allocator is NOT a
        // replacement (it was falsely flagged), and a signature broken after the return type IS one (it was
        // silently missed — a private counter formatted that way would have reintroduced P52 through CI).
        ['struct T { void* operator new (std::size_t n) { return std::malloc (n); } };',            0],
        ['struct Arena; struct T { void* operator new (std::size_t n, Arena&); };',                 0],
        ['namespace ns { void* operator new (std::size_t s) { return nullptr; } }',                 0],
        ['void f() { void* p = ::operator new (8); ::operator delete (p); }',                       0],
        ['void*\noperator new (std::size_t n) { return std::malloc (n ? n : 1); }',                 1],
        ['void\noperator delete (void* p) noexcept { std::free (p); }',                             1],
        ['inline\nvoid* operator new[] (std::size_t n) { return nullptr; }',                        1],
        // ...and the shapes the two crew seats got past the return-type matcher. Every one compiles and
        // every one is a working replacement.
        ['[[nodiscard]] void* operator new (std::size_t n) { return nullptr; }',                     1],
        ['__attribute__((malloc)) void* operator new (std::size_t n) { return nullptr; }',           1],
        ['extern "C++" void* operator new (std::size_t n) { return nullptr; }',                      1],
        ['auto operator new (std::size_t n) -> void* { return nullptr; }',                           1],
        ['void* operator new (std::size_t n, std::align_val_t a, const std::nothrow_t&) noexcept;',  1],
        ['void* g = operator new (8);',                                                              0],
        ['void* h() { return operator new (8); }',                                                   0],
        ['void* operator new (std::size_t s) { return std::malloc (s ? s : 1); }',                 1],
        ['void* operator new[] (std::size_t s) { return std::malloc (s); }',                       1],
        ['void* operator new (std::size_t s, std::align_val_t a) { return nullptr; }',             1],
        ['void* operator new (std::size_t s, const std::nothrow_t&) noexcept { return nullptr; }', 1],
        ['void  operator delete (void* p) noexcept { std::free (p); }',                            1],
        ['void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }',             1],
        ['void  operator delete (void* p, std::align_val_t) noexcept { std::free (p); }',          1],
        ['inline void* operator new (std::size_t s);',                                             1],
        ['void*operator new(std::size_t s){return 0;}',                                            1],
    ];
    let bad = 0;
    for (const [src, want] of cases)
    {
        const got = findViolations(src).length ? 1 : 0;
        if (got !== want) { console.error(`  SELF-TEST FAIL: wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    if (bad) { console.error(`alloc-counter lint self-test: ${bad}/${cases.length} cases wrong`); process.exit(1); }
    console.log(`alloc-counter lint self-test: ${cases.length}/${cases.length} cases correct`);
}

const args = process.argv.slice(2);
if (args.includes('--self-test')) { selfTest(); process.exit(0); }

const DEFAULT_ROOTS = ['modules', 'tools', 'test_support'];
const roots = args.length ? args : DEFAULT_ROOTS;
const files = [];
for (const r of roots) { try { statSync(r).isDirectory() ? walk(r, files) : files.push(r); } catch { console.error(`no such path: ${r}`); process.exit(2); } }

let found = 0, allowed = 0;
for (const f of files)
{
    const norm = f.split(sep).join('/');
    const hits = findViolations(readFileSync(f, 'utf8'));
    if (ALLOWED.test(norm)) { allowed += hits.length; continue; }
    for (const h of hits)
    {
        console.error(`${norm}:${h.line}: ${h.what} replaced outside test_support/alloc_counter.h — P52.`);
        console.error(`    A private counter replaces the two default-aligned forms and misses every over-aligned`);
        console.error(`    allocation, so "no heap allocation" silently stops meaning it. Include <alloc_counter.h>`);
        console.error(`    instead (exactly one TU per executable) and read alloc::count / alloc::bytes.`);
        found++;
    }
}

if (found) { console.error(`\n^^ ${found} private allocation-counter definition(s) across ${files.length} files in ${roots.join(', ')}`); process.exit(1); }
// The lint's own anchor check. A scan that never reaches the shared header would report a clean tree having
// proved nothing — a lint keyed on one thing goes blind on that one thing silently. Only meaningful for the
// default roots; an explicit path list is a human narrowing the scan on purpose.
if (allowed === 0 && roots === DEFAULT_ROOTS)
{
    console.error(`the shared counter defines nothing — test_support/alloc_counter.h was not scanned or no longer`);
    console.error(`defines the replacements, so a clean result here would prove nothing. Fix the scan, not this line.`);
    process.exit(1);
}
console.log(`one allocation counter: ${allowed} replacements, all in test_support/alloc_counter.h (${files.length} files scanned)`);
