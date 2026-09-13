### `tools` — release notes move to `changelog.d/`, one file per task

Every branch that changes behaviour now writes its note as a NEW FILE under `changelog.d/` instead of appending to
`## Unreleased`. Git cannot conflict on two branches adding two different files; it conflicted on that one section
five times in a single day, and each time cost a full rebase-build-push-wait round on work that was already green.

`node tools/changelog-collect.mjs --release vX.Y.Z` folds the fragments into `CHANGELOG.md` under the new heading and
deletes them — one commit on the release branch, where there is nobody to conflict with. `--preview` prints what the
next release would say, which is how the accumulated notes stay readable now that they live apart.
