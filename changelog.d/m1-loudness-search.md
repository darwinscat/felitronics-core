<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### mastering — the loudness search stops at the limit it may not break

**`mastering::TargetLoudnessSolver`**, for the three limits that grow with drive — the limiter's gain reduction, the
peak-to-loudness ratio and the loudness-range loss:

- **A target that needs such a limit broken is answered at the limit.** Once a render breaks one, the next renders
  bracket the largest drive that keeps it — between the loudest render that kept it (or the drive the limiter starts
  working at) and the quietest that broke it, by regula falsi on the limit's own statistic — and the call delivers the
  loudest render that keeps it as `TargetUnreachable`, with `binding` the limit broken at the smallest breaking drive.
  The search stops when that render is within `toleranceLu` of the smallest breaking one at the ceiling its true peak
  asks for — unless that one is inside the target's tolerance, where a render that keeps the limit is still `Solved` —
  or at `maxPasses`. It used to step on loudness alone and deliver the feasible render it happened to have: after a
  first render at the starting gain, that render — the source at its own level.
- **A limit the target keeps changes nothing**: the renders and the audio are bit-identical to the same request
  without it. A limit the first render already breaks is still `UpstreamViolation`, and now also under a ceiling below
  the promise when the limiter is idle on that render — neither knob moves it there. A limit whose statistic is not
  monotone in drive can end the search short of it; the delivered render keeps every limit whenever a render that
  keeps them was seen, and `Solved` still means the delivered render keeps them all.
- **The first render whose limiter works, after an idle one, is aimed 0.15 dB under the aim**, or by the largest
  overshoot measured on a working render at or under its drive. It used to be aimed at the promise itself, land over
  it by the between-sample overshoot and spend a render on the correction: on the suite's programme a -12 LUFS target
  now takes two renders where it took three.

No field and no ABI change.
