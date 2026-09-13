// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Folds `changelog.d/*.md` into CHANGELOG.md. See changelog.d/README.md for WHY the fragments exist;
// this file is only the folding.
//
//   node tools/changelog-collect.mjs --preview            what the next release would say, writes nothing
//   node tools/changelog-collect.mjs --release v0.32.0    folds, renames the heading, deletes the fragments
//
// Order is by the leading task id (p58 before p60, p9 before p10 — numeric, not lexicographic), then by
// name for anything that does not carry one. README.md is not a fragment.

import { readdirSync, readFileSync, writeFileSync, unlinkSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const dir  = join(root, 'changelog.d');
const book = join(root, 'CHANGELOG.md');

const args   = process.argv.slice(2);
const preview = args.includes('--preview');
const version = (args[args.indexOf('--release') + 1] || '').trim();

if (preview === (args.includes('--release')))
{
    console.error('usage: changelog-collect.mjs --preview | --release vX.Y.Z');
    process.exit(2);
}
if (! preview && ! /^v\d+\.\d+\.\d+$/.test(version))
{
    console.error(`--release needs a version like v0.32.0, got "${version}"`);
    process.exit(2);
}

// `p57b` is task 57's second half and belongs right after 57, so the id is a NUMBER and an optional
// letter — not `\b`-terminated, which is what made `p57b` and `p59b` unparseable and sorted them to the
// end, behind p64, in the first release this tool saw.
const rank = (name) => {
    const m = /^p(\d+)([a-z]*)/i.exec(name);
    return m ? [Number(m[1]), m[2].toLowerCase()] : [Number.MAX_SAFE_INTEGER, ''];
};

const files = readdirSync(dir)
    .filter(f => f.endsWith('.md') && f !== 'README.md')
    .sort((a, b) => {
        const [na, sa] = rank(a), [nb, sb] = rank(b);
        return (na - nb) || sa.localeCompare(sb) || a.localeCompare(b);
    });

if (files.length === 0) { console.error('changelog.d/ holds no fragments'); process.exit(preview ? 0 : 1); }

const folded = files
    .map(f => readFileSync(join(dir, f), 'utf8').replace(/\s+$/, ''))
    .join('\n\n');

if (preview)
{
    console.log(`# ${files.length} fragment(s), in release order\n`);
    for (const f of files) console.log(`  ${f}`);
    console.log('\n' + '-'.repeat(78) + '\n');
    console.log(folded);
    process.exit(0);
}

const date = new Date().toISOString().slice(0, 10);
let text = readFileSync(book, 'utf8');

// Either there is an `## Unreleased` section to close, or the release opens its own above the newest one.
if (/^## Unreleased\s*$/m.test(text))
{
    const start = text.search(/^## Unreleased\s*$/m);
    const rest  = text.slice(start);
    const next  = rest.slice(1).search(/^## /m);           // the heading of the previous release, if any
    const body  = next === -1 ? rest : rest.slice(0, next + 1);
    const tail  = next === -1 ? ''   : rest.slice(next + 1);
    const kept  = body.replace(/^## Unreleased\s*$/m, '').replace(/\s+$/, '');
    text = text.slice(0, start)
         + `## ${version} — ${date}\n`
         + (kept ? kept + '\n\n' : '\n')
         + folded + '\n\n'
         + tail;
}
else
{
    const first = text.search(/^## /m);
    const at    = first === -1 ? text.length : first;
    text = text.slice(0, at) + `## ${version} — ${date}\n\n` + folded + '\n\n' + text.slice(at);
}

writeFileSync(book, text);
for (const f of files) unlinkSync(join(dir, f));
console.log(`CHANGELOG.md: ${version} — ${date}, ${files.length} fragment(s) folded and removed.`);
