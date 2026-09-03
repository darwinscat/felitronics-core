<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Contributing to felitronics-core

**Short version: issues yes, code no.**

## Bug reports, questions and ideas — welcome

Open an issue. "Your release ballistics misbehave below 8 kHz", a question about an API, a
suggestion for a module — all of these are genuinely useful and get acted on. If a fix or a
feature lands because you raised it, you are credited in `CHANGELOG.md`.

## Code — not accepted

Pull requests containing code are closed unmerged, automatically, regardless of how good the
change is. This is a standing policy, not a verdict on your patch.

**Why.** felitronics-core is AGPL-3.0-or-later and the copyright is held solely by its two
authors, Oleh Tsymaienko and Alisa Lafoks. A merged contribution leaves its author holding
copyright on those lines permanently. From that moment the project can no longer change its own
licensing terms — dual-licensing, granting an exception, shipping it inside something else —
without tracking down every past contributor and obtaining each one's agreement. Nobody can
foresee a decade of a project's life; declining code costs nothing today, so that is the side
we err on.

If your patch would have been genuinely good, **describe the change in an issue instead**.
Ideas are not copyrightable — only their expression is — so a well-described problem or approach
can be implemented here freely. You get the fix you wanted, with credit, and the project keeps
its options open.

## What this policy does NOT do

It takes away none of the freedoms the AGPL grants you. You may use this library, study it,
modify it, and distribute your own modified version under the same licence. Forking is expressly
fine. This is only about what gets merged *here*.

## Third-party code

The optional compiled backends vendor third-party code (pffft, NeuralAmpModelerCore, Eigen,
nlohmann/json, namz) under permissive licences, recorded in `THIRD_PARTY_NOTICES.md`. Those are
governed by their own terms and are unaffected by the policy above. New third-party dependencies
are not added without review — see the watch-list at the bottom of that file.
