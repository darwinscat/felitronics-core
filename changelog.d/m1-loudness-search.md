<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering — the loudness search stops at the limit it may not break

**`mastering::TargetLoudnessSolver`**, for the three limits that grow with drive — the limiter's gain reduction, the
peak-to-loudness ratio and the loudness-range loss:

- **A target past such a limit is answered at the limit.** When the next step asks for a drive at or past the
  smallest one a render broke a limit at, the search renders inside the bracket between that render and the loudest
  one that kept every limit instead — regula falsi on the limit's own statistic, held inside the bracket, the
  midpoint where a statistic is not measured — and delivers the loudest render that keeps them as
  `TargetUnreachable`, `binding` what the smallest breaking drive broke. It stops once the breaking render, at the ceiling its true peak asks for, is under the
  target's tolerance and within `toleranceLu` of that render, or at `maxPasses`. It used to step on loudness alone and
  deliver the feasible render it happened to have: after a first render at the starting gain, the source at its own
  level. Nothing is assumed about how a limit moves with drive beyond where to look: every render the bound chooses is
  measured, and a render that keeps a limit above one that broke it ends the bracket.
- **Whenever the search without the limit solves on a render that keeps it, the search with the limit is the same
  search**, render for render and bit for bit. Without a render that kept every limit there is no bracket, and the
  search is the one before.
- **The first render whose limiter works, after an idle one, is aimed 0.15 dB under the aim**, or by the largest
  overshoot measured on a working render at or under its drive. It used to be aimed at the promise itself, land over
  it by the between-sample overshoot and spend a render on the correction: on the suite's programme a -12 LUFS target
  now takes two renders where it took three.

No field and no ABI change.
