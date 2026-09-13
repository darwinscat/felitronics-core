<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `changelog.d/` — one file per task, so two branches never collide

Every branch that changes behaviour writes its release note **here, as a new file**, and leaves
`CHANGELOG.md` alone. Git cannot conflict on two branches adding two different files; it conflicts every
single time they append to the same section of the same file, which is what `## Unreleased` used to be.
That cost five rebase rounds in one day before this directory existed.

**Name the file for the task**: `p58-clip-detector.md`, `p60-compressor-mix.md` — the task id, a dash, a
couple of words. The id decides the order the notes appear in the release.

**Write exactly what used to go under `## Unreleased`** and nothing more: a `### module · module — the
headline` line, then the body. No version number, no date, no `## Unreleased` heading — the release adds
those.

**At release time** `node tools/changelog-collect.mjs --release vX.Y.Z` folds every fragment into
`CHANGELOG.md` under the new version heading and deletes them, in one commit on the release branch,
where there is nobody to conflict with. `--preview` prints what the next release would say without
touching anything — that is how you read the accumulated notes now that they live apart.
