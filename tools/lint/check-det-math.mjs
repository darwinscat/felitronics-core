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
// The analyzers whose output is diffed byte-for-byte against the wasm module (felitronics-mastering-core's,
// linted by this script with --satellite) contained, when this was written, ZERO direct
// system transcendental calls — all thirty of their transcendental calls already spell `det::` (11 pow10,
// 11 log10, 4 log2, 2 sin, 1 exp2, 1 cos). Their entire remaining exposure was INDIRECT, through four
// ordinary-looking functions:
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
// itself. It was tried here and it is wrong, for a reason worth writing down: `fcore_measure` (now in
// felitronics-mastering-core) links BOTH regimes into one binary — its `report` mode is the deterministic analyzer and its `master` mode is the
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
// Usage: node tools/lint/check-det-math.mjs [--self-test] [--report] [--propose] [--satellite]
//   --self-test  run the matcher's own negative controls and exit
//   --report     print the full inventory (file, line, scope, call) and exit 0 — for an audit, not a gate
//   --propose    print manifest lines for files that have none, marked UNCLASSIFIED. It never writes the
//                manifest and never marks anything allowed: a human types the reason or the build stays red.
//   --satellite  run from ANOTHER repository's root (felitronics-mastering-core, felitronics-guitar-core),
//                with this script taken from the felitronics-core checkout that repository builds against.
//                Two passes, both whole: core's tree against the lists in this file and core's manifest,
//                exactly as core's CI runs it; then the satellite's ./modules and ./tools against ITS lists
//                (tools/lint/det-math-zone.txt — see ZONE FILE below) and its own manifest, with core's
//                carriers. Nothing is skipped on either side; see TWO REPOSITORIES, ONE GATE below.

import { readFileSync, readdirSync, statSync, existsSync, realpathSync } from 'node:fs';
import { join, dirname, relative, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

// Everything below the matcher runs ONLY when this file is invoked as a program. Without this an
// `import` of it — which the self-test cases and any harness that wants scanText() must do — ran the
// whole gate as a side effect, printed its verdict and could exit the importing process.
const RUN_AS_PROGRAM = process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1];

// The marker check needs the opposite view from the matcher: comments KEPT, strings blanked. Reading the
// raw line instead let a string literal containing "// libm-ok:" authorise the call beside it, which an
// adversarial round used. Everything here is deliberately the mirror of stripNonCode().
export function stripStringsKeepComments (src)
{
    let out = '';
    const blank = (s) => s.replace(/[^\n]/g, ' ');
    for (let i = 0; i < src.length;)
    {
        const two = src.slice(i, i + 2);
        if (two === '//') { const e = src.indexOf('\n', i); const end = e < 0 ? src.length : e; out += src.slice(i, end); i = end; continue; }
        if (two === '/*') { const e = src.indexOf('*/', i + 2); const end = e < 0 ? src.length : e + 2; out += src.slice(i, end); i = end; continue; }
        if (src[i] === '\'' && /[0-9a-fA-F]/.test(src[i - 1] || '') && /[0-9a-fA-F]/.test(src[i + 1] || '')) { out += ' '; i++; continue; }
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
        // A ' BETWEEN DIGITS IS A SEPARATOR, NOT A CHARACTER LITERAL. `if (n < 1'000) return std::cos(x);`
        // otherwise opens a "literal" that swallows the rest of the line, and the call in it vanishes —
        // measured: that exact line scanned as zero hits before this case existed.
        if (src[i] === '\'' && /[0-9a-fA-F]/.test(src[i - 1] || '') && /[0-9a-fA-F]/.test(src[i + 1] || ''))
        { out += ' '; i++; continue; }
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
// `\s*` and not `[ \t]*`: a call may put its opening parenthesis on the NEXT LINE, and a line-by-line
// matcher misses it entirely. Measured — `return std::cos\n    (x);` scanned as zero hits, in the zone,
// with no marker needed because no call was recorded. Everything below therefore matches against the whole
// stripped TEXT and derives the line number from the match index.
// `\\)?\\s*\\(` after the name: `(std::cos)(x)` is a legal call that wraps the name in parentheses — it
// suppresses macro expansion and reads as a deliberate way round a name — and without this it matched
// nothing at all. Measured: an adversarial round used exactly that spelling to move a printed dBTP.
const CALL_RE   = new RegExp(`(?<![A-Za-z0-9_])(?:(std::|::)|(?<![.>:]))(${FN_ALT})(f?)\\s*\\)?\\s*\\(`, 'g');
const CPLX_RE   = new RegExp(`(?<![A-Za-z0-9_])std::(${COMPLEX_FNS.join('|')})\\s*\\(`, 'g');
// std::abs and std::sqrt on a COMPLEX are hypot and a complex square root; the scalar overloads are exact
// and deliberately ungoverned, so these are matched only where the file spells std::complex at all.
const CPLX_OVER = /(?<![A-Za-z0-9_])std::(abs|sqrt)\s*\(/g;
// `pow10` is deliberately absent from SCALAR — it is not a libm name, so there is no `std::pow10` to ban.
// But `det::pow10` exists and is the single most-used deterministic call in the analyzers, so the det side
// has to know it or every count this file prints is short by a third.
const DET_RE    = new RegExp(`(?<![A-Za-z0-9_])det::(${FN_ALT}|pow10)\\s*\\(`, 'g');
// A DEFINITION is not a call. `inline double cos (double x) noexcept` names the function being defined;
// counting it made `core` read as 33 sites when it has 20. The tell is a return type right before the name.
// A DEFINITION is not a call. The return type must begin the LINE (after keywords): anchored only at a
// `(`/`,`/`;` it also matched `g(T * std::cos(x))`, where `T *` is a multiplication — measured as zero
// hits. `T` as a bare type name is dropped for the same reason; template code spells it in a signature
// that starts the line anyway.
const DEFN_RE = /^\s*(?:\[\[[^\]]*\]\]\s*)?(?:(?:static|inline|constexpr|friend|virtual|explicit|template\s*<[^>]*>)\s+)*(?:const\s+)?(?:unsigned\s+|signed\s+)?(?:double|float|int|long|short|char|bool|auto|void|size_t|std::size_t)\s*[&*]?\s*$/;

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

// The parity entry points of THIS repository: none since the C ABIs moved to felitronics-mastering-core, whose own
// list names them (its tools/lint/det-math-zone.txt). A tree that grows a native-vs-wasm comparison again
// lists its entry point here, and rule 4's net comes back with it.
const ENTRY_POINTS = [];

// THE DETERMINISTIC ZONE of this repository. Every file here produces, or is directly consumed by, a number
// that is compared BYTE FOR BYTE between a native build and the wasm module. Those comparisons run in
// felitronics-mastering-core, whose analyzers and probe stand on these four files; its own zone lists the
// rest, and its CI runs this script with --satellite, which runs this tree's pass as well. Membership is a
// claim about where a value GOES, which is why it is written rather than computed — see the note at the top
// of this file about the closure.
const ZONE = new Set([
    // the reference true-peak instrument: what certifies a delivered file, and what the probe prints
    'modules/analysis/include/felitronics/analysis/ReferenceTruePeakMeter.h',
    // the shared floor the offline analyzers stand on
    'modules/core/include/felitronics/core/OfflineFft.h',
    'modules/oversampling/include/felitronics/oversampling/PolyphaseOversampler.h',
    // P31: not diffed by a parity harness today (no tool builds it), but it promises the same filter bits on
    // every row — its suite pins them by hash on each one — and it sits in the mastering ABI's include
    // closure through Saturator.h. A system sine in its design would break that promise silently; here it is red.
    'modules/oversampling/include/felitronics/oversampling/CascadeOversampler.h',
]);

// IN-ZONE EXCEPTIONS, each of which had to be argued rather than waved through. A line in the zone may
// call libm only if it carries a `// libm-ok: <reason>` marker AND appears here; a marker without an
// entry is red, and an entry whose marker has gone is red. Two locks, because one of them is a comment.
const ZONE_EXCEPTIONS = [
    { file: 'modules/core/include/felitronics/core/OfflineFft.h', fn: 'abs', count: 1,
      why: '`std::abs(std::complex<double>)` IS std::hypot — measured, identical checksums to an explicit hypot on Apple, glibc and musl, and three DIFFERENT checksums between those rows. There is no det::hypot to move it to, and sqrt(norm(z)) is not a rewrite, it is a different (less accurate, differently-overflowing) function. It is recorded rather than converted because magSpectrum feeds analysis/offline/SpectrumCurve and measurement/CaptureGate, neither of which is in a byte diff: the analyzers that ARE diffed (felitronics-mastering-core) take their magnitudes from SpectrumFrames, which uses re*re+im*im and calls no libm at all.' },
];
//==============================================================================
// TWO REPOSITORIES, ONE GATE. felitronics-core owns this script, the carriers and the deterministic floor;
// a satellite repository runs the SAME script from its own root with --satellite. Every file is named by the
// tree it lives in and its path inside that tree, never by a bare relative path: both repositories have a
// modules/analysis/..., and a list entry that meant one of them must not quietly cover the other.
//   core pass       core's ./modules and ./tools against the lists above and core's manifest — the whole
//                   gate, on the core checkout the satellite builds with, so a satellite run cannot be green
//                   on a core whose own lists have rotted.
//   satellite pass  the satellite's ./modules and ./tools against its lists and its manifest, with core's
//                   carriers. Its parity entry points' #include closure walks INTO core's headers; a core file
//                   reached that way is judged by the core pass, and labelled there as reachable from one.
// Each pass checks its own lists for rot (a zone entry, an entry point or an exception naming a file or a call
// that is not there any more), so neither side's allowances can outlive what they allowed.
const CORE_ROOT = realpathSync(join(dirname(fileURLToPath(import.meta.url)), '..', '..'));

// THE ZONE FILE of a satellite: tools/lint/det-math-zone.txt at its root. One entry per line, `#` comments:
//     entry      <path>                  <why>    a parity entry point: its #include closure is rule 4's net
//     zone       <path>                  <why>    a file whose numbers are compared byte for byte across rows
//     exception  <path>  <fn>*<count>    <why>    an in-zone libm call, argued for (the line marker is required too)
// A repository with no such file has no zone and no entry points: every libm call in it is a manifest line.
// That is a legitimate state (felitronics-guitar-core), so an ABSENT file is not an error — which is why the
// repository that does have a zone pins it with a planted-violation control in its own CI.
const ZONE_FILE = 'tools/lint/det-math-zone.txt';

export function parseZoneFile (text)
{
    const lists = { entryPoints: [], zone: [], exceptions: [], errors: [] };
    const seen = new Set();
    text.split('\n').forEach ((raw, i) =>
    {
        const line = raw.trim();
        if (! line || line.startsWith('#')) return;
        const err = (msg) => lists.errors.push({ line: i + 1, msg });
        let m;
        if ((m = /^(entry|zone)\s+(\S+)\s+(\S.*)$/.exec(line)))
        {
            if (seen.has(m[1] + ' ' + m[2])) return err(`${m[1]} ${m[2]} is listed twice`);
            seen.add(m[1] + ' ' + m[2]);
            (m[1] === 'entry' ? lists.entryPoints : lists.zone).push({ path: m[2], why: m[3].trim() });
        }
        else if ((m = /^exception\s+(\S+)\s+([A-Za-z0-9_]+)\*([1-9][0-9]*)\s+(\S.*)$/.exec(line)))
            lists.exceptions.push({ file: m[1], fn: m[2], count: Number(m[3]), why: m[4].trim() });
        else
            err(`unparseable line: ${line} — expected "entry <path> <why>", "zone <path> <why>" or "exception <path> <fn>*<count> <why>", and every one of them needs its reason`);
    });
    return lists;
}

//==============================================================================
function walk (dir, acc) { for (const e of readdirSync(dir)) { const p = join(dir, e); if (statSync(p).isDirectory()) walk(p, acc); else if (/\.(h|hpp|inl|c|cc|cpp|cxx)$/.test(e)) acc.push(p); } return acc; }

function relOf (root, p) { return relative(root, p).split(sep).join('/'); }

function sourceFiles (root)
{
    const files = [];
    const modules = join(root, 'modules');
    for (const m of (existsSync(modules) ? readdirSync(modules) : []))
        for (const sub of ['include', 'src'])
        { const p = join(modules, m, sub); try { if (statSync(p).isDirectory()) walk(p, files); } catch { /* module has no such dir */ } }
    // tools/ IS scanned, deliberately and unlike the long-double lint. The parity SURFACE lives there:
    // felitronics-mastering-core's tools/fcore_probe.h computes the dBTP that CI diffs, with its own floor
    // and, once, its own log10. A modules-only audit would have declared that path clean while the printed
    // number was unpinned.
    for (const p of [join(root, 'tools')]) { try { if (statSync(p).isDirectory()) walk(p, files); } catch {} }
    return files.map (f => relOf(root, f)).filter (f => ! /\/(tests|bench)\//.test (f)).sort();
}

function includeRootsOf (root)
{
    const roots = [];
    const modules = join(root, 'modules');
    for (const m of (existsSync(modules) ? readdirSync(modules) : []))
    { const p = join(modules, m, 'include'); try { if (statSync(p).isDirectory()) roots.push(p); } catch {} }
    return roots;
}

// Rule 4's net: everything the parity entry points can #include, across every tree on the include path (the
// running repository's first, then core's). Returns Map<tree, Set<path in that tree>>. An include of
// felitronics/... that resolves in NO tree, or in TWO, is a violation rather than a silent gap: a header
// this lint cannot find is one it cannot audit, and one that two repositories both provide is one where it
// may be auditing a different file from the one the compiler takes.
function computeClosure (trees, violations, display)
{
    const reached = new Map(trees.map (t => [t, new Set()]));
    const queue = [];
    const add = (t, rel) => { if (! reached.get(t).has(rel)) { reached.get(t).add(rel); queue.push([t, rel]); } };
    for (const t of trees) for (const e of t.entryPoints) if (existsSync(join(t.root, e))) add(t, e);
    while (queue.length)
    {
        const [t, rel] = queue.shift();
        const abs = join(t.root, rel);
        let text; try { text = readFileSync(abs, 'utf8'); } catch { continue; }
        // Comments only. stripNonCode() also blanks STRING literals, and `#include "fcore_probe.h"` IS a
        // string literal — so every quoted include vanished before this regex saw it, and the closure
        // silently lost the tools' own headers, which is where the parity surface lives.
        const code = text.replace(/\/\*[\s\S]*?\*\//g, m => m.replace(/[^\n]/g, ' '))
                         .replace(/\/\/[^\n]*/g, m => m.replace(/[^\n]/g, ' '));
        for (const m of code.matchAll(/#[ \t]*include[ \t]*[<"]([^">]+)[">]/g))
        {
            const inc = m[1];
            if (inc.startsWith('felitronics/'))
            {
                const hits = [];
                for (const tt of trees)
                    for (const r of tt.includeRoots)
                    { const p = join(r, inc); if (existsSync(p)) { hits.push([tt, relOf(tt.root, p)]); break; } }
                const line = lineAt(code, m.index);
                if (hits.length === 0)
                    violations.push({ f: display(t, rel), line, rule: 'CLOSURE',
                                      msg: `#include <${inc}> resolves in no include root this lint knows. A header it cannot find is a header it cannot audit.` });
                else if (hits.length > 1)
                    violations.push({ f: display(t, rel), line, rule: 'CLOSURE',
                                      msg: `#include <${inc}> resolves in two repositories: ${hits.map (([tt, r]) => display(tt, r)).join(' and ')}. One header, one owner — otherwise this lint may audit a different file from the one the compiler takes.` });
                if (hits.length) add(hits[0][0], hits[0][1]);
            }
            else
            {
                // beside the file, then on the tools include path of its own tree: tools/wasm/fc_probe.cpp
                // includes "fcore_probe.h", which lives in tools/ and arrives through -I, not as a sibling.
                const sib = join(dirname(abs), inc), viaTools = join(t.root, 'tools', inc);
                const p = existsSync(sib) ? sib : (existsSync(viaTools) ? viaTools : null);
                if (p) add(t, relOf(t.root, p));
            }
        }
    }
    return reached;
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

// Line number from a match index, without splitting the text into lines (which is what made the matcher
// blind to a call whose parenthesis sits on the next line).
function lineAt (text, idx) { let n = 1; for (let i = 0; i < idx; i++) if (text[i] === '\n') n++; return n; }

// The enclosing function, for the DIAGNOSTIC and for reading the manifest — never for deciding whether a
// call is allowed.
function enclosingAt (text, idx)
{
    const upto = text.slice(0, idx).split('\n');
    return enclosing(upto, upto.length - 1);
}

export function scanText (text, opts = {})
{
    const code = stripNonCode(text);
    const hits = [];
    const hasComplex = /std::complex/.test(code);
    const push = (m, fn, ns, kind, suffix) =>
        hits.push({ line: lineAt(code, m.index), fn: fn + (suffix || ''), ns, kind, scope: enclosingAt(code, m.index) });

    CALL_RE.lastIndex = 0;
    for (const m of code.matchAll(CALL_RE))
    {
        const lineStart = code.lastIndexOf('\n', m.index) + 1;
        if (DEFN_RE.test(code.slice(lineStart, m.index))) continue;          // a definition, not a call
        // The `f` suffix is KEPT in the name. Normalising `sinf` to `sin` would let one be swapped for the
        // other without moving a manifest line, and they are different functions with different rounding.
        push (m, m[2], m[1] || '', 'scalar', m[3]);
    }
    CPLX_RE.lastIndex = 0;
    for (const m of code.matchAll(CPLX_RE)) push (m, m[1], 'std::', 'complex');
    if (hasComplex && opts.complexAbs !== false)
    {
        CPLX_OVER.lastIndex = 0;
        for (const m of code.matchAll(CPLX_OVER)) push (m, m[1], 'std::', 'complex-overload');
    }
    return hits;
}

// A CALL to a carrier — `core::gainToDb(x)` or the unqualified `gainToDb(x)`. NOT `o.gainToDb(x)`,
// `p->gainToDb(x)` or `other::gainToDb(x)`: all three matched before the lookbehind below, because the
// suffix alone was enough and the qualification was optional.
//
// AND NOT A FILE'S OWN FUNCTION OF THE SAME NAME. `dynamics/NoiseGate.h` defines a private
// `static float dbToGain (float)` and calls it unqualified. Counting that as `core::dbToGain` put four
// phantom carrier calls into its manifest line — a manifest that says a file reaches libm through a function it never calls is
// wrong in the direction that matters, because it reads as an audited fact. Their own `std::pow` is
// already counted as a direct call; this only stops it being counted twice under someone else's name.
export function scanCarriers (text, names)
{
    const code = stripNonCode(text);
    const hits = [];
    if (! names.length) return hits;
    // TWO alternatives, not one with an optional prefix: a single `(?:core::)?` followed by a lookbehind
    // that forbids `::` rejects `core::gainToDb` as well, which silently disables this whole rule. Either
    // the call is explicitly `core::NAME`, or it is a bare NAME with nothing qualifying it.
    const alt = names.join('|');
    // Names this file DEFINES itself: an unqualified call to one of them is a call to the local one.
    const own = new Set();
    for (const m of code.matchAll(new RegExp(`(?:^|\\n)[^\\n]*?\\b(?:double|float|int|auto)\\s+(${alt})\\s*\\(`, 'g')))
        own.add(m[1]);
    const re = new RegExp(`(?:(?<![A-Za-z0-9_])core::(${alt})|(?<![A-Za-z0-9_.>:])(${alt}))\\s*\\)?\\s*\\(`, 'g');
    for (const m of code.matchAll(re))
    {
        const lineStart = code.lastIndexOf('\n', m.index) + 1;
        if (DEFN_RE.test(code.slice(lineStart, m.index))) continue;
        const qualified = m[1] !== undefined;
        const name = m[1] || m[2];
        if (! qualified && own.has(name)) continue;          // the file's own function of that name
        hits.push({ line: lineAt(code, m.index), fn: name, scope: enclosingAt(code, m.index) });
    }
    return hits;
}

//==============================================================================
// THE MANIFEST. One line per file, not one per call — 52 lines a human re-reads, not 195 nobody does.
// Each line pins the exact MULTISET of functions that file calls, so a call added, removed, or swapped
// for another moves the line and forces a review; whitespace and line moves do not.
//     <path>  <fn*count fn*count ...>  <disposition>  <reason>
// Each repository has its own, at this path under its root.
const MANIFEST_PATH = 'tools/lint/det-math-manifest.txt';
const DISPOSITIONS = new Set(['retain-rt', 'retain-local', 'retain-offline', 'owner-decision', 'UNCLASSIFIED']);

function parseManifest (root)
{
    const path = join(root, MANIFEST_PATH);
    if (! existsSync(path)) return { entries: new Map(), missing: true };
    const entries = new Map();
    for (const raw of readFileSync(path, 'utf8').split('\n'))
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
        // --- the holes a review round found, each now a case so it cannot reopen ---
        ['double y = std::cos\n    (x);',                                                 1],   // ( on the NEXT line
        ['double y = g(T * std::cos (x));',                                               1],   // `T *` is a multiply, not a return type
        ["double y = (n < 1'000) ? std::cos (x) : 0.0;",                                  1],   // digit separator is not a char literal
        ['std::complex<double> z; auto r = std::sqrt (z);',                               1],   // complex sqrt IS transcendental
        ['double r = std::sqrt (x);',                                                     0],   // scalar sqrt is exact — never governed
        ['inline double cos (double x)\n{ return 1.0; }',                                 0],   // definition split over lines
        // --- the bypasses the ADVERSARIAL round executed, each now a case ---
        ['double y = (std::cos)(x);',                                                     1],   // the name in parentheses
        ['double y = (std::cos) (x);',                                                    1],
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
        ['const double d = o.gainToDb (x);',        0],   // a MEMBER of the same name is not this function
        ['const double d = p->gainToDb (x);',       0],
        ['const double d = other::gainToDb (x);',   0],   // another namespace's function of the same name
        ['const double d = core::gainToDb\n    (x);', 1],  // ( on the next line, the carrier half
        ['const double d = (core::gainToDb)(x);',   1],   // ...and with the name in parentheses
    ];
    for (const [src, want] of carrierCases)
    {
        const got = scanCarriers(src, ['gainToDb', 'dbToGain']).length;
        if (got !== want) { console.error(`  SELF-TEST FAIL (carrier): wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    // The MARKER lexer is the mirror of the matcher's: comments kept, strings blanked. A string that
    // contains "// libm-ok:" must not authorise the call beside it — an adversarial round used exactly
    // that, and the check was reading the raw line.
    const markerCases = [
        ['double y = std::cos (x);   // libm-ok: a real comment',  true ],
        ['const char* s = "// libm-ok: not a comment";',           false],
        ['const char* s = "/* libm-ok: nor this */";',             false],
    ];
    for (const [src, want] of markerCases)
    {
        const got = /\/\/[^\n]*libm-ok:/.test(stripStringsKeepComments(src));
        if (got !== want) { console.error(`  SELF-TEST FAIL (marker): wanted ${want}, got ${got} for: ${JSON.stringify(src)}`); bad++; }
    }
    // The ZONE FILE parser: a satellite's lists. A line it cannot read is an error, never an empty list —
    // a misspelt keyword that parsed as "nothing" would be a zone that covers nothing and reads as coverage.
    const zoneCases = [
        ['zone modules/a/include/x.h the analyzer whose report is diffed',            { zone: 1, errors: 0 }],
        ['entry tools/wasm/fc_probe.cpp the wasm probe',                              { entryPoints: 1, errors: 0 }],
        ['exception modules/a/include/x.h abs*1 std::abs of a complex is hypot',      { exceptions: 1, errors: 0 }],
        ['# a comment\n\n   ',                                                         { errors: 0 }],
        ['zones modules/a/include/x.h misspelt keyword',                              { zone: 0, errors: 1 }],
        ['zone modules/a/include/x.h',                                                { zone: 0, errors: 1 }],   // no reason
        ['exception modules/a/include/x.h abs the count is missing',                  { exceptions: 0, errors: 1 }],
        ['zone a.h why\nzone a.h why again',                                          { zone: 1, errors: 1 }],   // listed twice
    ];
    for (const [src, want] of zoneCases)
    {
        const got = parseZoneFile(src);
        for (const [k, n] of Object.entries(want))
            if (got[k].length !== n) { console.error(`  SELF-TEST FAIL (zone file): wanted ${n} ${k}, got ${got[k].length} for: ${JSON.stringify(src)}`); bad++; }
    }
    const total = cases.length + carrierCases.length + markerCases.length + zoneCases.length;
    if (bad) { console.error(`det-math lint self-test: ${bad} of ${total} cases wrong`); process.exit(1); }
    console.log(`det-math lint self-test: ${total}/${total} cases correct`);
}

//==============================================================================
// ONE PASS over one tree: rules 1-3 against that tree's own lists and manifest, plus the rot checks of those
// lists. `closure` is the set of this tree's files reached from any parity entry point (rule 4's net), and
// `carriers` the carrier NAMES to look for — core's, wherever the tree is.
function lintTree (t, closure, carriers, entryPointsShown, display, violations)
{
    const V = (f, line, rule, msg) => violations.push({ f: display(t, f), line, rule, msg });
    // EVERY FILE THE PARITY CLOSURE REACHES IS SCANNED, whatever it is called. The directory walk takes only
    // the usual source extensions, and a `#include "tables.inc"` from a zone file was read by the closure and
    // audited by nobody: a system call planted in it passed on both sides while the same call one level up
    // failed. The closure is the net, so what it catches is inventoried.
    const files = [...new Set(sourceFiles(t.root).concat([...closure].filter (f => ! /\/(tests|bench)\//.test (f))))].sort();
    const scanned = new Set(files);
    const carrierNames = carriers.map(c => c.name);
    const inventory = [];
    const perFile = new Map();
    const exceptionUses = new Map();

    // THE LISTS MUST NAME FILES THAT ARE THERE. A zone entry for a file that moved covers nothing, and the
    // count printed at the end kept reading as if it did — the move of the mastering modules out of this
    // repository is exactly the change that would have left twelve such entries behind, green.
    for (const f of t.zone)
        if (! scanned.has(f))
            V(f, 0, 'ZONE-ROT', `the deterministic zone names a file this lint does not scan (moved, renamed, or under tests/). A zone entry that covers nothing reads as coverage — fix the path or remove the entry.`);
    for (const f of t.entryPoints)
        if (! existsSync(join(t.root, f)))
            V(f, 0, 'ZONE-ROT', `a parity entry point that does not exist: its #include closure is empty, and rule 4 with it.`);
    for (const f of t.implementation)
        if (! scanned.has(f))
            V(f, 0, 'ZONE-ROT', `IMPLEMENTATION names a file this lint does not scan — the exemption now covers nothing.`);

    for (const f of files)
    {
        const text = readFileSync(join(t.root, f), 'utf8');
        const markerLines = stripStringsKeepComments(text).split('\n');   // comments kept, strings blanked
        const hits = scanText(text);
        const carrierHits = t.implementation.has(f) ? [] : scanCarriers(text, carrierNames);
        // CARRIER CALLS COUNT TOWARDS THE MANIFEST TOO, and this is not tidiness. Without it the rule only
        // looked inside the zone, and a file outside it could be reverted from `gainToDbDet` to `gainToDb`
        // with the gate still green — measured on LoudnessSolver::peakDb, which is the one place the
        // certificate and the report must agree bit for bit, i.e. exactly the regression this lint was
        // written after. They are spelled `name()` so they cannot collide with a scalar of the same name.
        const all = hits.concat (carrierHits.map (h => ({ ...h, fn: h.fn + '()', ns: '', kind: 'carrier' })));
        if (all.length) perFile.set(f, all);
        for (const h of all) inventory.push({ f, ...h });

        if (t.implementation.has(f)) continue;               // det:: itself and the dB definitions

        if (t.zone.has(f))
        {
            for (const h of hits)
            {
                // The marker must be ON the calling line, and the exception must be recorded. Either alone
                // is not enough: a comment anyone can type is not an approval, and an approval nobody can see
                // at the call site is not a warning.
                const marked = /\/\/[^\n]*libm-ok:/.test(markerLines[h.line - 1] || '');
                const listed = t.exceptions.find(e => e.file === f && e.fn === h.fn);
                if (marked && listed) { exceptionUses.set(f + '::' + h.fn, (exceptionUses.get(f + '::' + h.fn) || 0) + 1); continue; }
                if (marked && ! listed)
                { V(f, h.line, 'ZONE', `${h.ns}${h.fn}() carries a "libm-ok" marker but there is no entry for it in ZONE_EXCEPTIONS. A marker is a note to a reader; the entry is where the argument has to be written down.`); continue; }
                V(f, h.line, 'ZONE', `${h.ns}${h.fn}() in the deterministic zone (scope ${h.scope}). This file's numbers are compared BYTE FOR BYTE between the native CLI and the wasm module, and ${h.fn} is not the same function on those rows. Use core::det::${h.fn} — but check its DOMAIN first, which is not the same as std's: det::tan returns NaN at |x| >= 2^24 and det::pow returns NaN at x <= 0, and only det::cos/det::sin are pinned against the mpmath oracle (the other six are checked by identities), so a new use needs its own argument-range argument. If it genuinely must stay system, mark the line "// libm-ok: <why>" and add an entry to ZONE_EXCEPTIONS.`);
            }
            for (const h of carrierHits)
            {
                const c = carriers.find(x => x.name === h.fn);
                const marked = /\/\/[^\n]*libm-ok:/.test(markerLines[h.line - 1] || '');
                const listed = t.exceptions.find(e => e.file === f && e.fn === h.fn);
                if (marked && listed) { exceptionUses.set(f + '::' + h.fn, (exceptionUses.get(f + '::' + h.fn) || 0) + 1); continue; }
                V(f, h.line, 'CARRIER', `${h.fn}() in the deterministic zone (scope ${h.scope}) — ${c.why}. It reads as ordinary arithmetic and is a libm call: that is the whole reason this rule exists, because the analyzers here contain no direct std:: call at all and were exposed entirely through functions that look like this one.`);
            }
        }
    }

    // AND THE EXCEPTIONS MUST NOT ROT EITHER. One that no longer matches anything is an argument left
    // standing for a call that is gone — exactly the stale allowance this lint exists to prevent elsewhere.
    for (const e of t.exceptions)
    {
        const used = exceptionUses.get(e.file + '::' + e.fn) || 0;
        if (used === 0)
            V(e.file, 0, 'ZONE-EXCEPTION-ROT', `ZONE_EXCEPTIONS allows ${e.fn}() here, but no marked call to it was found. Either the call went away (remove the entry) or its "// libm-ok:" marker did (put it back) — an unused allowance is how a list stops meaning anything.`);
        else if (used !== e.count)
            V(e.file, 0, 'ZONE-EXCEPTION-ROT', `ZONE_EXCEPTIONS allows ${e.count} marked ${e.fn}() call(s) here and found ${used}. The written argument is about specific calls; another one needs its own, not a share of this one.`);
    }

    // Rule 4 — THE CARRIER LIST MUST NOT ROT. A carrier that no longer reaches libm would forbid something
    // harmless forever; one deleted from this list while still reaching libm would let the real thing through.
    // So each declared carrier is checked against the file that defines it — by the tree that defines it.
    for (const c of t.carriers)
    {
        const path = join(t.root, c.defined);
        if (! existsSync(path))
        { V(c.defined, 0, 'CARRIER-ROT', `carrier ${c.name} names a file that does not exist`); continue; }
        const body = stripNonCode(readFileSync(path, 'utf8'));
        const lines = body.split('\n');
        const defLine = lines.findIndex(l => new RegExp(`\\b${c.name}\\b\\s*\\(`).test(l) && DEFN_RE.test(l.slice(0, l.indexOf(c.name))));
        if (defLine < 0)
            V(c.defined, 0, 'CARRIER-ROT', `carrier ${c.name} is declared in this lint but no definition of it was found in ${c.defined}. Either it moved (update CARRIERS) or it is gone (remove it) — a carrier list nobody checks is a list that stops being true.`);
        else
        {
            // AND THE BODY MUST STILL REACH libm. Checking only that a definition EXISTS leaves the other
            // half of the rot: a carrier that was converted to det:: would go on forbidding something
            // harmless in the zone forever, and the comment above claimed this was checked when it was not.
            // The window is the definition line plus the few that can hold a one-expression body.
            const window = lines.slice(defLine, defLine + 6).join('\n');
            if (scanText(window).length === 0)
                V(c.defined, defLine + 1, 'CARRIER-ROT', `carrier ${c.name} no longer calls a system transcendental in its first lines. If it was converted, it is not a carrier any more — remove it from CARRIERS, or the zone keeps refusing a call that is now safe.`);
        }
    }

    // Rule 3 — the manifest, for everything outside the zone.
    const { entries, missing } = parseManifest(t.root);
    if (missing)
        V(MANIFEST_PATH, 0, 'MANIFEST', 'the manifest does not exist; run with --propose to generate a starting point (every line UNCLASSIFIED until a human writes a reason)');
    else
    {
        for (const [k, v] of entries)
            if (v === null) V(MANIFEST_PATH, 0, 'MANIFEST', `unparseable line: ${k.replace('__PARSE_ERROR__', '')}`);

        for (const [f, hits] of perFile)
        {
            if (t.zone.has(f) || t.implementation.has(f)) continue;
            const want = multisetOf(hits);
            const e = entries.get(f);
            if (! e)
            { V(f, hits[0].line, closure.has(f) ? 'CLOSURE' : 'MANIFEST',
                `calls libm [${want}] and has no manifest entry.`
                + (closure.has(f) ? ` THIS FILE IS REACHABLE BY #include FROM A PARITY ENTRY POINT (${entryPointsShown.join(', ')}), so it is compiled into the binaries whose outputs CI diffs — read it before classifying it.`
                                  : ` It is not reachable from any parity entry point, so this is bookkeeping rather than a hazard — but it is still a decision somebody made.`)
                + ` Add a line to ${display(t, MANIFEST_PATH)} (--propose prints a starting one, marked UNCLASSIFIED so it cannot pass by accident).`); continue; }
            // `--propose` emits a placeholder reason. Changing only the DISPOSITION in front of it and leaving
            // the placeholder made the gate green with nobody having written anything — measured.
            if (/^<.*>$/.test(e.reason) || e.reason.includes('why they stay system'))
                V(f, hits[0].line, 'MANIFEST', `the manifest entry still carries --propose's placeholder reason. The disposition is not the classification; the sentence after it is.`);
            if (e.disposition === 'UNCLASSIFIED')
                V(f, hits[0].line, 'MANIFEST', `manifest entry is still UNCLASSIFIED — someone has to say what these calls are and why they stay`);
            else if (! DISPOSITIONS.has(e.disposition))
                V(f, 0, 'MANIFEST', `unknown disposition "${e.disposition}" (expected one of ${[...DISPOSITIONS].join(', ')})`);
            if (e.multiset !== want)
                V(f, hits[0].line, 'MANIFEST', `the calls in this file CHANGED.\n      manifest: [${e.multiset}]\n      actual:   [${want}]\n      A libm call was added, removed or swapped for another. Re-read the file, then update the line — this is the review, not an obstacle to it.`);
        }
        for (const [f, e] of entries)
        {
            if (f.startsWith('__PARSE_ERROR__')) continue;
            if (! perFile.has(f) || t.zone.has(f))
                V(f, 0, 'MANIFEST', ! existsSync(join(t.root, f)) ? 'manifest names a file that no longer exists — remove the line'
                                  : t.zone.has(f)                   ? 'this file is now inside the deterministic zone, where the ban applies and a manifest entry means nothing — remove the line and convert the calls'
                                                                    : 'manifest entry for a file with no libm calls left — remove the line (a stale allowance is how a list stops meaning anything)');
        }
    }
    return { t, inventory, perFile, entries: entries || new Map() };
}

//==============================================================================
const args = process.argv.slice(2);
if (! RUN_AS_PROGRAM) { /* imported for its matcher; the gate below is not ours to run */ }
else {
if (args.includes('--self-test')) { selfTest(); process.exit(0); }

// WHERE THIS RUNS, stated rather than assumed. Without --satellite the tree is core itself, so the working
// directory must BE this script's checkout; with it, the working directory is another repository and core is
// the checkout this script came from. Either mix-up used to lint the wrong tree against the wrong lists.
const SATELLITE = args.includes('--satellite');
const CWD = realpathSync(process.cwd());
if (! existsSync(join(CORE_ROOT, 'modules/core/include/felitronics/core/DetMath.h')))
{ console.error(`check-det-math: ${CORE_ROOT} does not look like a felitronics-core checkout (no core/DetMath.h) — this script must be run from inside one`); process.exit(2); }
if (SATELLITE && CWD === CORE_ROOT)
{ console.error(`check-det-math: --satellite is for ANOTHER repository's root; this is felitronics-core's own (${CORE_ROOT}). Run it without --satellite.`); process.exit(2); }
if (! SATELLITE && CWD !== CORE_ROOT)
{ console.error(`check-det-math: run from felitronics-core's root (${CORE_ROOT}), or pass --satellite to lint the repository in ${CWD} against this core.`); process.exit(2); }

const violations = [];
const core = { name: 'felitronics-core', root: CORE_ROOT, zone: ZONE, entryPoints: ENTRY_POINTS, exceptions: ZONE_EXCEPTIONS,
               carriers: CARRIERS, implementation: IMPLEMENTATION, includeRoots: includeRootsOf(CORE_ROOT) };
let local = core;
if (SATELLITE)
{
    const zonePath = join(CWD, ZONE_FILE);
    const lists = existsSync(zonePath) ? parseZoneFile(readFileSync(zonePath, 'utf8')) : { entryPoints: [], zone: [], exceptions: [], errors: [] };
    for (const e of lists.errors) violations.push({ f: ZONE_FILE, line: e.line, rule: 'LISTS', msg: e.msg });
    local = { name: 'this repository', root: CWD, zone: new Set(lists.zone.map(z => z.path)), entryPoints: lists.entryPoints.map(e => e.path),
              exceptions: lists.exceptions, carriers: [], implementation: new Set(), includeRoots: includeRootsOf(CWD),
              zoneFile: existsSync(zonePath) };
}
const trees = SATELLITE ? [local, core] : [core];
const display = (t, rel) => t.root === CWD ? rel : relOf(CWD, join(t.root, rel));
const entryPointsShown = trees.flatMap (t => t.entryPoints.map (e => display(t, e)));

const reached = computeClosure(trees, violations, display);          // rule 4's net — NOT the ban set; see the note at the top
const results = trees.slice().reverse().map (t =>                     // core's pass first, then the satellite's
    lintTree(t, reached.get(t), CARRIERS, entryPointsShown, display, violations));

//==============================================================================
if (args.includes('--report'))
{
    for (const { t, inventory, perFile } of results)
    {
        console.log(`# det-math inventory${SATELLITE ? ' [' + t.name + ']' : ''} — ${inventory.length} governed calls in ${perFile.size} files`);
        console.log(`# deterministic zone: ${t.zone.size} files (a written list — see ZONE; the #include closure of`);
        console.log(`#   ${entryPointsShown.join(', ') || '(no parity entry points)'} is the DISCOVERY net, not the zone)`);
        for (const h of inventory) console.log(`${display(t, h.f)}:${h.line}\t${h.ns}${h.fn}\t${h.scope}\t${t.zone.has(h.f) ? 'ZONE' : 'outside'}`);
    }
    process.exit(0);
}
if (args.includes('--propose'))
{
    const { t, perFile, entries } = results[results.length - 1];      // the running repository's own tree
    console.log(`# proposed lines for files with no manifest entry — EVERY ONE IS UNCLASSIFIED ON PURPOSE.`);
    console.log(`# The lint stays red until a human replaces UNCLASSIFIED with a disposition and a real reason.`);
    for (const [f, hits] of perFile)
    {
        if (t.zone.has(f) || t.implementation.has(f) || entries.has(f)) continue;
        console.log(`${f}  [${multisetOf(hits)}]  UNCLASSIFIED  <what these are, and why they stay system>`);
    }
    process.exit(0);
}

const summary = ({ t, inventory, perFile }) =>
    `${inventory.length} governed libm calls in ${perFile.size} files; ${t.zone.size}-file deterministic zone`;
if (violations.length)
{
    for (const v of violations)
        console.error(`${v.f}${v.line ? ':' + v.line : ''}: [${v.rule}] ${v.msg}`);
    console.error(`\n^^ ${violations.length} violation(s). ` + results.map (r => (SATELLITE ? `[${r.t.name}] ` : '')
        + `${r.inventory.length} governed calls in ${r.perFile.size} files; deterministic zone is ${r.t.zone.size} files.`).join(' '));
    process.exit(1);
}
if (! SATELLITE)
    console.log(`det-math: ${summary(results[0])} is clean (direct + carriers); manifest matches.`);
else
    for (const r of results)
        console.log(`det-math [${r.t.name}${r.t === core ? ' at ' + (relOf(CWD, CORE_ROOT) || '.') : ''}]: ${summary(r)} is clean (direct + carriers); manifest matches.`
                    + (r.t === local ? ` Lists: ${local.zoneFile ? ZONE_FILE : 'none (no ' + ZONE_FILE + ')'}; ${local.entryPoints.length} parity entry point(s), closure ${[...reached.values()].reduce((n, s) => n + s.size, 0)} files, ${reached.get(core).size} of them in felitronics-core.` : ''));
}
