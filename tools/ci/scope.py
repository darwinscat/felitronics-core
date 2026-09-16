#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which tests does a pull request need? Scope by SUBJECT, derived from the tree.

    python3 tools/ci/scope.py <base> [<head>]   a pull request: the files that differ between base and head
    python3 tools/ci/scope.py --all             anything else: a push, a release branch, a schedule, a dispatch

Prints a summary and, under GitHub Actions, writes to $GITHUB_OUTPUT:

    mode    none      only prose changed (docs/, changelog.d/, *.md) -- nothing to test
            selected  only modules/ changed -- those modules, everything that includes them, and nothing else
            all       anything else: tools/, CMake, CI, the root, or not a pull request at all
    regex   the ctest -R pattern for `selected`, anchored; empty otherwise
    abi     true when the selection reaches a module that tools/ is built on, so the native-vs-wasm parity
            and the ABI tests must run; the plumbing itself (tools/wasm) changing is already `all`

Every mapping is DERIVED, none is listed by hand -- in this repository a hand list is the first thing to go
stale:
    test    -> module   the CMakeLists.txt under modules/<m>/ that declares it with add_test(NAME ...)
    module  -> users    #include <felitronics/<m>/...> anywhere in modules/<user>/ (headers, sources, tests)
    tools   -> modules  the same includes read from tools/, then closed forward through the modules' own

A selection that resolves to no test at all is not trusted: it falls back to `all`, because an empty filter
reads as a green run that tested nothing.
"""
import os
import re
import subprocess
import sys
from collections import defaultdict, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INCLUDE = re.compile(r'#\s*include\s*[<"]felitronics/([a-z0-9]+)/')
ADD_TEST = re.compile(r'add_test\s*\(\s*NAME\s+([A-Za-z0-9_]+)')
SOURCE_SUFFIXES = {'.h', '.hpp', '.cpp', '.cc', '.inl'}


def is_prose(path):
    return path.startswith(('docs/', 'changelog.d/')) or path.endswith('.md')


def included_modules(directory):
    found = set()
    if not directory.is_dir():
        return found
    for f in directory.rglob('*'):
        if f.suffix in SOURCE_SUFFIXES and f.is_file():
            found |= set(INCLUDE.findall(f.read_text(encoding='utf-8', errors='replace')))
    return found


def declared_tests(directory):
    names = set()
    for cml in directory.rglob('CMakeLists.txt'):
        names |= set(ADD_TEST.findall(cml.read_text(encoding='utf-8', errors='replace')))
    return names


def emit(mode, regex='', abi=False, note=''):
    print(f'scope: {mode}' + (f' -- {note}' if note else ''))
    out = os.environ.get('GITHUB_OUTPUT')
    if out:
        with open(out, 'a', encoding='utf-8') as fh:
            fh.write(f'mode={mode}\nregex={regex}\nabi={"true" if abi else "false"}\n')


def main(argv):
    if argv[:1] == ['--all']:
        emit('all', note='not a pull request')
        return 0
    if not argv:
        print(__doc__)
        return 2
    base, head = argv[0], (argv[1] if len(argv) > 1 else 'HEAD')
    changed = subprocess.run(['git', 'diff', '--name-only', base, head], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout.split()

    modules = sorted(p.name for p in (ROOT / 'modules').iterdir() if p.is_dir())
    touched, elsewhere = set(), []
    for path in changed:
        if is_prose(path):
            continue
        m = re.match(r'modules/([^/]+)/', path)
        if m and m.group(1) in modules:
            touched.add(m.group(1))
        else:
            elsewhere.append(path)

    if elsewhere:
        shown = ', '.join(elsewhere[:3]) + (' ...' if len(elsewhere) > 3 else '')
        emit('all', note=f'outside modules/: {shown}')
        return 0
    if not touched:
        emit('none', note=f'{len(changed)} file(s), all prose')
        return 0

    users, forward = defaultdict(set), {}
    for m in modules:
        uses = included_modules(ROOT / 'modules' / m) - {m}
        forward[m] = uses
        for dep in uses:
            users[dep].add(m)

    closure, queue = set(touched), deque(touched)
    while queue:
        for user in users[queue.popleft()]:
            if user not in closure:
                closure.add(user)
                queue.append(user)

    reached, queue = set(), deque(included_modules(ROOT / 'tools'))
    while queue:
        m = queue.popleft()
        if m in reached or m not in forward:
            continue
        reached.add(m)
        queue.extend(forward[m])
    abi = bool(closure & reached)

    tests = set()
    for m in closure:
        tests |= declared_tests(ROOT / 'modules' / m)
    if abi:
        tests |= declared_tests(ROOT / 'tools')
    if not tests:
        emit('all', note=f'{sorted(closure)} declare no test -- not trusting an empty filter')
        return 0

    regex = '^(' + '|'.join(sorted(tests)) + ')$'
    emit('selected', regex, abi,
         note=f'changed {sorted(touched)}; with users {len(closure)} module(s), {len(tests)} test(s)'
              f'{"; reaches the ABI" if abi else ""}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
