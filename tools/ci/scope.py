#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which tests does a pull request need? Scope by SUBJECT, derived from the tree.

    python3 tools/ci/scope.py <base> [<head>]   a pull request: the files that differ between base and head
    python3 tools/ci/scope.py --all             anything else: a push, a release branch, a schedule, a dispatch

Prints a summary and, under GitHub Actions, writes to $GITHUB_OUTPUT:

    mode    none      only prose changed (docs/, changelog.d/, *.md, and a release's version line) -- nothing to test
            selected  only modules/ changed -- those modules, everything that includes them, and nothing else
            all       anything else: tools/, CMake, CI, the root, or not a pull request at all
    regex   the ctest -R pattern for `selected`, anchored; empty otherwise

Every mapping is DERIVED, none is listed by hand -- in this repository a hand list is the first thing to go
stale:
    test    -> module   the CMakeLists.txt under modules/<m>/ that declares it with add_test(NAME ...); a name
                        built from a variable (a foreach over test files) matches any value of that variable
    module  -> users    #include <felitronics/<m>/...> anywhere in modules/<user>/ (headers, sources, tests)

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
ADD_TEST = re.compile(r'add_test\s*\(\s*NAME\s+([A-Za-z0-9_${}]+)')
VARIABLE = re.compile(r'\$\{[^}]*\}')
SOURCE_SUFFIXES = {'.h', '.hpp', '.cpp', '.cc', '.inl'}


def is_prose(path):
    return path.startswith(('docs/', 'changelog.d/')) or path.endswith('.md')


def version_bump_only(base, head, path):
    # A release pull request folds the changelog and moves ONE line: `project(felitronics_core VERSION x.y.z
    # LANGUAGES CXX)`. That is prose for testing purposes -- the code is the tree main already tested -- and the
    # whole matrix runs on the release's merge to main, which is the tree the tag goes on.
    if path != 'CMakeLists.txt':
        return False
    diff = subprocess.run(['git', 'diff', '-U0', base, head, '--', path], cwd=ROOT,
                          capture_output=True, text=True, check=True).stdout.splitlines()
    changed = [l for l in diff if l[:1] in '+-' and not l.startswith(('+++', '---'))]
    return bool(changed) and all(re.fullmatch(r'[+-]project\(felitronics_core VERSION [0-9]+\.[0-9]+\.[0-9]+ LANGUAGES CXX\)', l)
                                 for l in changed)


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
        for name in ADD_TEST.findall(cml.read_text(encoding='utf-8', errors='replace')):
            # A test named through a CMake variable (felitronics_<m>_${tl}_tests) read literally is its prefix
            # alone, which matches no test, so its suites would silently leave every selection. Stand for any
            # value of the variable instead -- selecting too much is visible, selecting too little is green.
            names.add(VARIABLE.sub('[A-Za-z0-9_]*', name))
    return names


def emit(mode, regex='', note=''):
    print(f'scope: {mode}' + (f' -- {note}' if note else ''))
    out = os.environ.get('GITHUB_OUTPUT')
    if out:
        with open(out, 'a', encoding='utf-8') as fh:
            fh.write(f'mode={mode}\nregex={regex}\n')


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
        if is_prose(path) or version_bump_only(base, head, path):
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

    users = defaultdict(set)
    for m in modules:
        uses = included_modules(ROOT / 'modules' / m) - {m}
        for dep in uses:
            users[dep].add(m)

    closure, queue = set(touched), deque(touched)
    while queue:
        for user in users[queue.popleft()]:
            if user not in closure:
                closure.add(user)
                queue.append(user)

    tests = set()
    for m in closure:
        tests |= declared_tests(ROOT / 'modules' / m)
    if not tests:
        emit('all', note=f'{sorted(closure)} declare no test -- not trusting an empty filter')
        return 0

    regex = '^(' + '|'.join(sorted(tests)) + ')$'
    emit('selected', regex,
         note=f'changed {sorted(touched)}; with users {len(closure)} module(s), {len(tests)} test(s)')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
