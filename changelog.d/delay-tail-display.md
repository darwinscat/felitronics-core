### delay — `TapeDelay`'s tail and picture follow the knobs as they stand

Three small departures from the one-head class. They are about the taps, not the sound: the parity test still
holds the default layout's audio to the frozen class byte for byte through every run it plays — tiny REPEATS and
one-ulp OFFSET runs included. (A REPEATS that is kept differently can steer a later guard differently, so a
contrived walk of hundreds of millions of requests at the guard's edge can end on other samples.)

- **A REPEATS of exactly 0 always lands.** The setter's guard took a request for 0 for a leftover of ~1e-38 and kept
  that; the loop never heard it, but `tailSeconds()` counted echoes that were not there.
- **Every stored OFFSET refreshes the tail.** The offset is stored on every call, but the tail was worked out again
  only when it moved past the same loose guard, so a run of one-ulp moves walked the offset away from its tail.
- **`reset()` and `prepare()` show the targets.** `shownTimeMs`, `shownOffsetMs` and `shownHeadMs` show where the
  next call puts the heads, instead of where the stopped stream left them until the next call. A poison restart
  inside a call still shows the heads that call played; the next call shows the snap.

`targets()` returns the knobs as the setters left them; the tests rebuild a fresh stage from it to hold the tail to
the formula of the targets as they stand.
