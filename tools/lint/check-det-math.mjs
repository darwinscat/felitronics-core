// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// THE libm AUDIT, ENFORCED. DSP-ARCHITECTURE.md law 10 and the `core::det` header say the system libm is
// not the same function on every row — measured: over the ranges this library uses, Apple's and glibc's
// log10 disagree at 4325 of 200000 points, their tan at 7274 of 20918, their pow10 at 54342 of 132001
// (41 %). Every dB THRESHOLD is built from pow10, and a threshold that moves by an ulp moves a DECISION.
// P79/P80 routed the values that cross a platform boundary through `core::det`. This is what keeps them
// there, and what makes a new `std::cos` in the wrong place a red build instead of a discovery two
// releases later.
//
// WHY THIS IS NOT A LIST OF CALL SITES, WHICH IS THE OBVIOUS DESIGN AND THE WRONG ONE.
// The analyzers whose output is diffed byte-for-byte against the wasm module contain, today, ZERO direct
// system transcendental calls — all sixteen of their transcendental calls already spell `det::`. Their
// entire remaining exposure was INDIRECT, through four ordinary-looking functions:
//     core::gainToDb / core::dbToGain      (std::log10 / std::pow)
//     core::offline::fftInplace            (std::cos / std::sin twiddle seeds)
//     PolyphaseOversampler::designFilter   (std::sin FIR taps)
// A lint that greps for `std::` in those analyzers finds nothing and passes, which is exactly what
// happened for three releases. So this lint has TWO halves: a direct ban, and a CARRIER ban — a call to
// a function that is itself known to reach libm counts as reaching libm.
//
// THE FOUR RULES
//   1. ZONE     — inside the deterministic zone, a governed system call is red.
//   2. CARRIER  — inside the zone, a call to a declared carrier is red (that is the indirect half).
//   3. MANIFEST — every file OUTSIDE the zone that calls libm must appear in the manifest with the exact
//                 multiset of functions it calls. A new call, a removed call, a changed function: red.
//   4. CLOSURE  — every file the parity entry points can #include must be in the zone or the manifest.
//
// WHY THE ZONE IS A TYPED LIST AND THE CLOSURE IS ONLY A NET. The elegant design is to COMPUTE the zone
// as the #include closure of the three translation units whose outputs CI diffs, and let it maintain
// itself. It was tried here and it is wrong, for a reason worth writing down: `fcore_measure` links BOTH
// regimes into one binary — its `report` mode is the deterministic analyzer and its `master` mode is the
// shipped RT chain (DetMath.h says exactly this about it). So the closure sweeps in MatchedBiquad,
// Compressor, Smoother and EnvelopeFollower — 52 files, most of them TabbyEQ's and OrbitCab's own DSP,
// whose `std::tan` is their sound and must not be touched. A ban over the closure would have demanded
// the one change this task is forbidden to make.
//
// So the closure is kept, but demoted to what it is actually good at: DISCOVERY. Rule 4 says every file
// reachable from an entry point must be CLASSIFIED — in the zone (must be deterministic) or in the
// manifest (deliberately system, with a reason). A new `#include` that drags in unclassified libm is
// then red, which is the property the closure was wanted for, without the ban it cannot support.
// The zone itself is an explicit list because membership is a claim about where a number GOES, and no
// include graph knows that.
//
// WHAT THIS LINT CANNOT DO, said plainly rather than discovered later:
//   · It cannot tell a per-sample call from a coefficient call. That is a property of the CALLER's loop,
//     not of the callee's text, and no lexer resolves it. `cadence` in the manifest is a HUMAN
//     classification; this file checks that it is present and self-consistent, never that it is true.
//   · It cannot resolve C++ overloads. `std::abs` is integer abs, float abs, or COMPLEX abs (which is
//     `hypot`, a libm call) depending on its argument's type. It is governed only in files that mention
//     `std::complex`, and that is a heuristic, stated as one.
//   · It cannot see libm reached through a template parameter, a function pointer, or a macro. The
//     SystemMath/DetMath policy routing is pinned by static_assert in MathPolicyTests.cpp instead; that
//     is the right tool for it and this is not.
//
// Usage: node tools/lint/check-det-math.mjs [--self-test] [--report] [--propose]
//   --self-test  run the matcher's own negative controls and exit
//   --report     print the full inventory (file, line, scope, call) and exit 0 — for an audit, not a gate
//   --propose    print manifest lines for files that have none, marked UNCLASSIFIED. It never writes the
//                manifest and never marks anything allowed: a human types the reason or the build stays red.

import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';

//==============================================================================
// THE LEXER. Reused wholesale from check-no-long-double.mjs, and for the same reason: the words this
// lint matches occur constantly in PROSE. `core::det`'s own header says "std::cos" a dozen times while
// explaining why not to call it. Comments, strings, char and raw-string literals are blanked, newlines
// preserved so line numbers survive.
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

//==============================================================================
// WHAT IS GOVERNED. The scalar transcendentals, every spelling. `sqrt` is deliberately ABSENT: IEEE-754
// requires it correctly rounded, so it is the same function on every row and `core::DetMath` forwards it
// to std:: on purpose. `fabs`, `floor`, `ceil`, `round`, `fmod`, `copysign` likewise — exact operations,
// not approximations.
const SCALAR = ['cos','sin','tan','acos','asin','atan','atan2','cosh','sinh','tanh','acosh','asinh','atanh',
                'exp','exp2','exp10','expm1','log','log2','log10','log1p','pow','cbrt','hypot','erf','erfc',
                'tgamma','lgamma'];
// Reached through <complex>, and invisible to a by-name search for the list above. `std::polar(r,t)` is
// cos+sin; complex `std::exp/log/pow/sqrt` are transcendental; `std::arg` is atan2. Measured on this
// tree: `std::abs(std::complex<double>)` IS `std::hypot` on Apple, glibc and musl (same checksum as an
// explicit hypot call on each), and the three rows give three DIFFERENT checksums.
const COMPLEX_FNS = ['polar','arg','proj'];

const FN_ALT = SCALAR.join('|');
// std::fn( / ::fn( / bare fn( — but not a member call (x.fn / p->fn) and not some_other_ns::fn(.
const CALL_RE   = new RegExp(`(?<![A-Za-z0-9_])(?:(std::|::)|(?<![.>:]))(${FN_ALT})f?[ \\t]*\\(`, 'g');
const CPLX_RE   = new RegExp(`(?<![A-Za-z0-9_])std::(${COMPLEX_FNS.join('|')})[ \\t]*\\(`, 'g');
const CPLX_ABS  = /(?<![A-Za-z0-9_])std::abs[ \t]*\(/g;
const DET_RE    = new RegExp(`(?<![A-Za-z0-9_])det::(${FN_ALT})[ \\t]*\\(`, 'g');
// A DEFINITION is not a call. `inline double cos (double x) noexcept` names the function being defined;
// counting it made `core` read as 33 sites when it has 20. The tell is a return type right before the name.
const DEFN_RE = /(?:^|[;{}(,]|\b(?:static|inline|constexpr|friend|virtual|explicit))\s*(?:const\s+)?(?:unsigned\s+|signed\s+)?(?:double|float|int|long|short|char|bool|auto|void|size_t|std::size_t|T)\s*[&*]?\s*$/;

//==============================================================================
// THE CARRIERS — ordinary functions in this tree that REACH libm. A call to one of these from inside the
// zone is a libm call, and this is the half of the lint that the parity analyzers actually needed.
// Each entry names the file that defines it, so rule 4 can check the carrier still HAS a system call in
// it: a carrier that was converted and left in this list would silently forbid something harmless, and a
// carrier removed from the list while still reaching libm would silently allow the real thing.
const CARRIERS = [
    { name: 'gainToDb',  defined: 'modules/core/include/felitronics/core/Math.h',
      why: 'std::log10 — use core::gainToDbDet for a value that crosses a row boundary' },
    { name: 'dbToGain',  defined: 'modules/core/include/felitronics/core/Math.h',
      why: 'std::pow(10,x), and Apple clang lowers it to __exp10, a DIFFERENT function: measured, 27 of 20000 dB values differ from a genuine pow(10,x) in the same binary' },
];

// Files that define the deterministic floor itself, or that are the system policy BY DESIGN. These are
// the implementation, not consumers of it, so the ban cannot apply to them.
const IMPLEMENTATION = new Set([
    'modules/core/include/felitronics/core/DetMath.h',   // det:: itself, and SystemMath's deliberate std:: bodies
    'modules/core/include/felitronics/core/Math.h',      // defines both dB spellings; gainToDb's std::log10 is the point
]);

const ENTRY_POINTS = [
    'tools/wasm/fc_probe.cpp',     // the wasm probe — CI diffs its lufs/dbtp text against the native CLI
    'tools/wasm/fc_master.cpp',    // the mastering ABI
    'tools/fcore_measure.cpp',     // the native CLI, the other side of every one of those diffs
];

// THE DETERMINISTIC ZONE. Every file here produces, or is directly consumed by, a number that CI compares
// BYTE FOR BYTE between the native CLI and the wasm module (tools/wasm/*-parity.mjs, and the plain `diff`
// steps in .github/workflows/ci.yml). Membership is a claim about where a value GOES, which is why it is
// written rather than computed — see the note at the top of this file about the closure.
const ZONE = new Set([
    // the analyzers whose output is diffed, one per parity harness
    'modules/analysis/include/felitronics/analysis/ProgrammeReport.h',
    'modules/analysis/include/felitronics/analysis/SourceForensics.h',
    'modules/analysis/include/felitronics/analysis/HumDetector.h',
    'modules/analysis/include/felitronics/analysis/LowEnd.h',
    'modules/analysis/include/felitronics/analysis/BandBursts.h',
    'modules/analysis/include/felitronics/analysis/ClipDetector.h',
    'modules/analysis/include/felitronics/analysis/WaveformPeaks.h',
    'modules/analysis/include/felitronics/analysis/StereoColumns.h',
    'modules/analysis/include/felitronics/analysis/SpectrumFrames.h',
    'modules/analysis/include/felitronics/analysis/ReferenceTruePeakMeter.h',
    // the shared floor those analyzers stand on
    'modules/core/include/felitronics/core/OfflineFft.h',
    'modules/oversampling/include/felitronics/oversampling/PolyphaseOversampler.h',
    // and the tools that PRINT the diffed text. The long-double lint deliberately skips tools/; for this
    // lint that would be a hole exactly on the surface being defended — fcore_probe.h computes the dBTP
    // that `diff native.txt wasm.txt` compares, with its own floor and, until P80, its own std::log10.
    'tools/fcore_probe.h',
    'tools/fcore_clips.h',
]);

// IN-ZONE EXCEPTIONS, each of which had to be argued rather than waved through. A line in the zone may
// call libm only if it carries a `// libm-ok: <reason>` marker AND appears here; a marker without an
// entry is red, and an entry whose marker has gone is red. Two locks, because one of them is a comment.
const ZONE_EXCEPTIONS = [
    { file: 'modules/core/include/felitronics/core/OfflineFft.h', fn: 'abs',
      why: '`std::abs(std::complex<double>)` IS std::hypot — measured, identical checksums to an explicit hypot on Apple, glibc and musl, and three DIFFERENT checksums between those rows. There is no det::hypot to move it to, and sqrt(norm(z)) is not a rewrite, it is a different (less accurate, differently-overflowing) function. It is recorded rather than converted because magSpectrum feeds analysis/offline/SpectrumCurve and measurement/CaptureGate, neither of which is in a byte diff: the analyzers that ARE diffed take their magnitudes from SpectrumFrames, which uses re*re+im*im and calls no libm at all.' },
];
const INCLUDE_ROOTS = [];   // filled from modules/*/include below

//==============================================================================
function walk (dir, acc) { for (const e of readdirSync(dir)) { const p = join(dir, e); if (statSync(p).isDirectory()) walk(p, acc); else if (/\.(h|hpp|inl|c|cc|cpp|cxx)$/.test(e)) acc.push(p); } return acc; }

function sourceFiles ()
{
    const files = [];
    for (const m of readdirSync('modules'))
        for (const sub of ['include', 'src'])
        { const p = join('modules', m, sub); try { if (statSync(p).isDirectory()) walk(p, files); } catch { /* module has no such dir */ } }
    // tools/ IS scanned, deliberately and unlike the long-double lint. The parity SURFACE lives there:
    // tools/fcore_probe.h computes the dBTP that CI diffs, with its own floor and its own log10. A
    // modules-only audit would have declared that path clean while the printed number was unpinned.
    for (const p of ['tools']) { try { if (statSync(p).isDirectory()) walk(p, files); } catch {} }
    return files.filter (f => ! /\/(tests|bench)\//.test (f)).sort();
}

function resolveInclude (inc)
{
    for (const r of INCLUDE_ROOTS) { const p = join(r, inc); if (existsSync(p)) return p; }
    return null;
}

// The zone: everything reachable by #include from the entry points. Returns a Set of repo-relative paths.
function computeZone ()
{
    for (const m of readdirSync('modules'))
    { const p = join('modules', m, 'include'); try { if (statSync(p).isDirectory()) INCLUDE_ROOTS.push(p); } catch {} }
    const zone = new Set(), queue = [];
    for (const e of ENTRY_POINTS) if (existsSync(e)) { zone.add(e); queue.push(e); }
    while (queue.length)
    {
        const f = queue.shift();
        let text; try { text = readFileSync(f, 'utf8'); } catch { continue; }
        const code = stripNonCode(text);
        for (const m of code.matchAll(/#[ \t]*include[ \t]*[<"]([^">]+)[">]/g))
        {
            const inc = m[1];
            let p = null;
            if (inc.startsWith('felitronics/')) p = resolveInclude(inc);
            else { const sib = join(dirname(f), inc); if (existsSync(sib)) p = sib; }     // tools' own headers
            if (p && ! zone.has(p)) { zone.add(p); queue.push(p); }
        }
    }
    return zone;
}

//==============================================================================
// Nearest preceding line that looks like a function definition — for the DIAGNOSTIC and for the manifest's
// scope grouping. A heuristic, and only ever used to say WHERE something is, never to decide whether it is
// allowed.
function enclosing (lines, n)
{
    for (let i = n; i >= 0 && i > n - 400; i--)
    {
        const m = /^\s{0,8}(?:\[\[[^\]]*\]\]\s*)?(?:static\s+|inline\s+|constexpr\s+|virtual\s+|explicit\s+|friend\s+|template\s*<[^>]*>\s*)*(?:[A-Za-z_][\w:<>,\s*&]*?\s+)?([A-Za-z_]\w*)\s*\([^;]*\)\s*(?:const\s*)?(?:noexcept\s*)?(?:->[^{;]*)?\s*\{?\s*$/.exec(lines[i]);
        if (m && !['if','for','while','switch','return','else','catch','do'].includes(m[1])) return m[1];
    }
    return '?';
}

export function scanText (text, opts = {})
{
    const code = stripNonCode(text);
    const lines = code.split('\n');
    const hits = [];
    const hasComplex = /std::complex/.test(code);
    lines.forEach((line, i) =>
    {
        for (const m of line.matchAll(CALL_RE))
        {
            if (DEFN_RE.test(line.slice(0, m.index))) continue;                 // a definition, not a call
            // det::log2(...) inside DetMath's own log2 is recursion in the deterministic implementation,
            // not a libm call. Only the implementation file can contain these.
            hits.push({ line: i + 1, fn: m[2], ns: m[1] || '', kind: 'scalar', scope: enclosing(lines, i) });
        }
        for (const m of line.matchAll(CPLX_RE))
            hits.push({ line: i + 1, fn: m[1], ns: 'std::', kind: 'complex', scope: enclosing(lines, i) });
        if (hasComplex && opts.complexAbs !== false)
            for (const m of line.matchAll(CPLX_ABS))
                hits.push({ line: i + 1, fn: 'abs', ns: 'std::', kind: 'complex-abs', scope: enclosing(lines, i) });
    });
    return hits;
}

export function scanCarriers (text, names)
{
    const code = stripNonCode(text);
    const lines = code.split('\n');
    const hits = [];
    const alt = names.join('|');
    if (! alt) return hits;
    const re = new RegExp(`(?<![A-Za-z0-9_])(?:core::)?(${alt})[ \\t]*\\(`, 'g');
    lines.forEach((line, i) =>
    {
        for (const m of line.matchAll(re))
        {
            if (DEFN_RE.test(line.slice(0, m.index))) continue;
            hits.push({ line: i + 1, fn: m[1], scope: enclosing(lines, i) });
        }
    });
    return hits;
}

//==============================================================================
// THE MANIFEST. One line per file, not one per call — 52 lines a human re-reads, not 195 nobody does.
// Each line pins the exact MULTISET of functions that file calls, so a call added, removed, or swapped
// for another moves the line and forces a review; whitespace and line moves do not.
//     <path>  <fn*count fn*count ...>  <disposition>  <reason>
const MANIFEST_PATH = 'tools/lint/det-math-manifest.txt';
const DISPOSITIONS = new Set(['retain-rt', 'retain-local', 'retain-offline', 'owner-decision', 'UNCLASSIFIED']);

function parseManifest ()
{
    if (! existsSync(MANIFEST_PATH)) return { entries: new Map(), missing: true };
    const entries = new Map();
    for (const raw of readFileSync(MANIFEST_PATH, 'utf8').split('\n'))
    {
        const line = raw.trim();
        if (! line || line.startsWith('#')) continue;
        const m = /^(\S+)\s+\[([^\]]*)\]\s+(\S+)\s+(.+)$/.exec(line);
        if (! m) { entries.set('__PARSE_ERROR__' + line, null); continue; }
        entries.set(m[1], { multiset: m[2].trim(), disposition: m[3], reason: m[4].trim() });
    }
    return { entries, missing: false };
}

function multisetOf (hits)
{
    const c = {};
    for (const h of hits) c[h.fn] = (c[h.fn] || 0) + 1;
    return Object.keys(c).sort().map(k => `${k}*${c[k]}`).join(' ');
}

//==============================================================================
function selfTest ()
{
    const cases = [
        // [source, expected number of governed hits]
        ['// a comment that says std::cos and std::pow must not trip it',                 0],
        ['const char* s = "std::cos (x)";',                                               0],
        ['auto r = R"(std::sin (y))";',                                                   0],
        ['/* std::tan\n   across lines */ int ok = 1;',                                   0],
        ['double y = std::sqrt (x);',                                                     0],   // sqrt is exact — not governed
        ['double y = std::fabs (x) + std::floor (z);',                                    0],
        ['inline double cos (double x) noexcept { return 1.0; }',                         0],   // a DEFINITION
        ['double y = std::cos (x);',                                                      1],
        ['double y = std::pow (10.0, x / 20.0);',                                         1],
        ['double a = std::cos (x), b = std::sin (y);',                                    2],   // two on one line
        ['double y = ::log10 (x);',                                                       1],
        ['double y = cosf (x);',                                                          1],   // the float spelling
        ['double y = std::cos(x);',                                                       1],   // no space before (
        ['double y = obj.cos (x);',                                                       0],   // a member, not libm
        ['double y = p->sin (x);',                                                        0],
        ['double y = mine::cos (x);',                                                     0],   // another namespace
        ['auto z = std::polar (1.0, -w);',                                                1],   // hidden cos+sin
        ['std::complex<double> c; double m = std::abs (c);',                              1],   // hidden hypot
        ['double m = std::abs (x);',                                                      0],   // no <complex> in the TU
    ];
    let bad = 0;
    for (const [src, want] of cases)
    {
        const got = scanText(src).length;
        if (got !== want) { console.error(`  SELF-TEST FAIL: wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    // The carrier matcher has its own controls: it must see a call and ignore the definition.
    const carrierCases = [
        ['const double d = core::gainToDb (x);',   1],
        ['const double d = gainToDb (x);',         1],
        ['inline double gainToDb (double g) noexcept { return 0.0; }', 0],
        ['// gainToDb is mentioned here in prose',  0],
        ['const double d = core::gainToDbDet (x);', 0],   // the DETERMINISTIC spelling is not a carrier
    ];
    for (const [src, want] of carrierCases)
    {
        const got = scanCarriers(src, ['gainToDb', 'dbToGain']).length;
        if (got !== want) { console.error(`  SELF-TEST FAIL (carrier): wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    if (bad) { console.error(`det-math lint self-test: ${bad} of ${cases.length + carrierCases.length} cases wrong`); process.exit(1); }
    console.log(`det-math lint self-test: ${cases.length + carrierCases.length}/${cases.length + carrierCases.length} cases correct`);
}

//==============================================================================
const args = process.argv.slice(2);
if (args.includes('--self-test')) { selfTest(); process.exit(0); }

const closure = computeZone();          // rule 4's net — NOT the ban set; see the note at the top
const zone = ZONE;
const files = sourceFiles();
const carrierNames = CARRIERS.map(c => c.name);
const violations = [];
const inventory = [];
const perFile = new Map();
const usedExceptions = new Set();

for (const f of files)
{
    const text = readFileSync(f, 'utf8');
    const rawLines = text.split('\n');
    const hits = scanText(text);
    const carrierHits = IMPLEMENTATION.has(f) ? [] : scanCarriers(text, carrierNames);
    if (hits.length) perFile.set(f, hits);
    for (const h of hits) inventory.push({ f, ...h });

    if (IMPLEMENTATION.has(f)) continue;                 // det:: itself and the dB definitions

    if (zone.has(f))
    {
        for (const h of hits)
        {
            // The marker must be ON the calling line, and the exception must be recorded. Either alone
            // is not enough: a comment anyone can type is not an approval, and an approval nobody can see
            // at the call site is not a warning.
            const marked = /\/\/[^\n]*libm-ok:/.test(rawLines[h.line - 1] || '');
            const listed = ZONE_EXCEPTIONS.find(e => e.file === f && e.fn === h.fn);
            if (marked && listed) { usedExceptions.add(f + '::' + h.fn); continue; }
            if (marked && ! listed)
            { violations.push({ f, line: h.line, rule: 'ZONE',
                                msg: `${h.ns}${h.fn}() carries a "libm-ok" marker but there is no entry for it in ZONE_EXCEPTIONS. A marker is a note to a reader; the entry is where the argument has to be written down.` }); continue; }
            violations.push({ f, line: h.line, rule: 'ZONE',
                              msg: `${h.ns}${h.fn}() in the deterministic zone (scope ${h.scope}). This file's numbers are compared BYTE FOR BYTE between the native CLI and the wasm module, and ${h.fn} is not the same function on those rows. Use core::det::${h.fn} — but check its DOMAIN first, which is not the same as std's: det::tan returns NaN at |x| >= 2^24 and det::pow returns NaN at x <= 0, and only det::cos/det::sin are pinned against the mpmath oracle (the other six are checked by identities), so a new use needs its own argument-range argument. If it genuinely must stay system, mark the line "// libm-ok: <why>" and add an entry to ZONE_EXCEPTIONS.` });
        }
        for (const h of carrierHits)
        {
            const c = CARRIERS.find(x => x.name === h.fn);
            const marked = /\/\/[^\n]*libm-ok:/.test(rawLines[h.line - 1] || '');
            const listed = ZONE_EXCEPTIONS.find(e => e.file === f && e.fn === h.fn);
            if (marked && listed) { usedExceptions.add(f + '::' + h.fn); continue; }
            violations.push({ f, line: h.line, rule: 'CARRIER',
                              msg: `${h.fn}() in the deterministic zone (scope ${h.scope}) — ${c.why}. It reads as ordinary arithmetic and is a libm call: that is the whole reason this rule exists, because the analyzers here contain no direct std:: call at all and were exposed entirely through functions that look like this one.` });
        }
    }
}

// AND THE EXCEPTIONS MUST NOT ROT EITHER. One that no longer matches anything is an argument left
// standing for a call that is gone — exactly the stale allowance this lint exists to prevent elsewhere.
for (const e of ZONE_EXCEPTIONS)
    if (! usedExceptions.has(e.file + '::' + e.fn))
        violations.push({ f: e.file, line: 0, rule: 'ZONE-EXCEPTION-ROT',
                          msg: `ZONE_EXCEPTIONS allows ${e.fn}() here, but no marked call to it was found. Either the call went away (remove the entry) or its "// libm-ok:" marker did (put it back) — an unused allowance is how a list stops meaning anything.` });

// Rule 4 — THE CARRIER LIST MUST NOT ROT. A carrier that no longer reaches libm would forbid something
// harmless forever; one deleted from this list while still reaching libm would let the real thing through.
// So each declared carrier is checked against the file that defines it.
for (const c of CARRIERS)
{
    if (! existsSync(c.defined))
    { violations.push({ f: c.defined, line: 0, rule: 'CARRIER-ROT', msg: `carrier ${c.name} names a file that does not exist` }); continue; }
    const body = stripNonCode(readFileSync(c.defined, 'utf8'));
    const defLine = body.split('\n').findIndex(l => new RegExp(`\\b${c.name}\\b[ \\t]*\\(`).test(l) && DEFN_RE.test(l.slice(0, l.indexOf(c.name))));
    if (defLine < 0)
        violations.push({ f: c.defined, line: 0, rule: 'CARRIER-ROT',
                          msg: `carrier ${c.name} is declared in this lint but no definition of it was found in ${c.defined}. Either it moved (update CARRIERS) or it is gone (remove it) — a carrier list nobody checks is a list that stops being true.` });
}

// Rule 3 — the manifest, for everything outside the zone.
const { entries, missing } = parseManifest();
if (missing)
    violations.push({ f: MANIFEST_PATH, line: 0, rule: 'MANIFEST', msg: 'the manifest does not exist; run with --propose to generate a starting point (every line UNCLASSIFIED until a human writes a reason)' });
else
{
    for (const [k, v] of entries)
        if (v === null) violations.push({ f: MANIFEST_PATH, line: 0, rule: 'MANIFEST', msg: `unparseable line: ${k.replace('__PARSE_ERROR__', '')}` });

    for (const [f, hits] of perFile)
    {
        if (zone.has(f) || IMPLEMENTATION.has(f)) continue;
        const want = multisetOf(hits);
        const e = entries.get(f);
        if (! e)
        { violations.push({ f, line: hits[0].line, rule: closure.has(f) ? 'CLOSURE' : 'MANIFEST',
                            msg: `calls libm [${want}] and has no manifest entry.`
                                 + (closure.has(f) ? ` THIS FILE IS REACHABLE BY #include FROM A PARITY ENTRY POINT (${ENTRY_POINTS.join(', ')}), so it is compiled into the binaries whose outputs CI diffs — read it before classifying it.`
                                                   : ` It is not reachable from any parity entry point, so this is bookkeeping rather than a hazard — but it is still a decision somebody made.`)
                                 + ` Add a line to ${MANIFEST_PATH} (--propose prints a starting one, marked UNCLASSIFIED so it cannot pass by accident).` }); continue; }
        if (e.disposition === 'UNCLASSIFIED')
            violations.push({ f, line: hits[0].line, rule: 'MANIFEST', msg: `manifest entry is still UNCLASSIFIED — someone has to say what these calls are and why they stay` });
        else if (! DISPOSITIONS.has(e.disposition))
            violations.push({ f, line: 0, rule: 'MANIFEST', msg: `unknown disposition "${e.disposition}" (expected one of ${[...DISPOSITIONS].join(', ')})` });
        if (e.multiset !== want)
            violations.push({ f, line: hits[0].line, rule: 'MANIFEST',
                              msg: `the calls in this file CHANGED.\n      manifest: [${e.multiset}]\n      actual:   [${want}]\n      A libm call was added, removed or swapped for another. Re-read the file, then update the line — this is the review, not an obstacle to it.` });
    }
    for (const [f, e] of entries)
    {
        if (f.startsWith('__PARSE_ERROR__')) continue;
        if (! perFile.has(f) || zone.has(f))
            violations.push({ f, line: 0, rule: 'MANIFEST',
                              msg: ! existsSync(f) ? 'manifest names a file that no longer exists — remove the line'
                                 : zone.has(f)     ? 'this file is now inside the deterministic zone, where the ban applies and a manifest entry means nothing — remove the line and convert the calls'
                                                   : 'manifest entry for a file with no libm calls left — remove the line (a stale allowance is how a list stops meaning anything)' });
    }
}

//==============================================================================
if (args.includes('--report'))
{
    console.log(`# det-math inventory — ${inventory.length} governed calls in ${perFile.size} files`);
    console.log(`# deterministic zone: ${zone.size} files, the #include closure of ${ENTRY_POINTS.join(', ')}`);
    for (const h of inventory) console.log(`${h.f}:${h.line}\t${h.ns}${h.fn}\t${h.scope}\t${zone.has(h.f) ? 'ZONE' : 'outside'}`);
    process.exit(0);
}
if (args.includes('--propose'))
{
    console.log(`# proposed lines for files with no manifest entry — EVERY ONE IS UNCLASSIFIED ON PURPOSE.`);
    console.log(`# The lint stays red until a human replaces UNCLASSIFIED with a disposition and a real reason.`);
    for (const [f, hits] of perFile)
    {
        if (zone.has(f) || IMPLEMENTATION.has(f) || entries.has(f)) continue;
        console.log(`${f}  [${multisetOf(hits)}]  UNCLASSIFIED  <what these are, and why they stay system>`);
    }
    process.exit(0);
}

if (violations.length)
{
    for (const v of violations)
        console.error(`${v.f}${v.line ? ':' + v.line : ''}: [${v.rule}] ${v.msg}`);
    console.error(`\n^^ ${violations.length} violation(s). ${inventory.length} governed calls in ${perFile.size} files; deterministic zone is ${zone.size} files.`);
    process.exit(1);
}
console.log(`det-math: ${inventory.length} governed libm calls in ${perFile.size} files; ${zone.size}-file deterministic zone is clean (direct + carriers); manifest matches.`);
