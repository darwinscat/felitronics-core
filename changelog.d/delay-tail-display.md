### delay — `TapeDelay`'s tail and picture follow the knobs as they stand

Three small departures from the one-head class, none of which moves a sample (the parity test still holds the
default layout's audio to the frozen class byte for byte):

- **A REPEATS of exactly 0 always lands.** The setter's guard took a request for 0 for a leftover of ~1e-38 and kept
  that; the loop never heard it, but `tailSeconds()` counted echoes that were not there (0.54 s instead of 0.35 s
  at the default TIME).
- **Every stored OFFSET refreshes the tail.** The offset is stored on every call, but the tail was worked out again
  only when it moved past the same loose guard, so a run of one-ulp moves walked the offset away from its tail.
- **`reset()` shows the targets.** `shownTimeMs`, `shownOffsetMs` and `shownHeadMs` show where the next call puts
  the heads from the reset on, instead of where the stopped stream left them until the next call.

`targets()` returns the knobs as the setters left them; the tests rebuild a fresh stage from it to hold the tail to
the formula of the targets as they stand.
