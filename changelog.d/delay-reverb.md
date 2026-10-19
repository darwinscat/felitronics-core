### delay · reverb — two new modules: the amp plugins' echo and room, moved out of the product

`felitronics::delay` and `felitronics::reverb` — header-only, JUCE-free; `delay` depends on `core`, `reverb` on
`core` and `eq`. They are OrbitAmp's delay and reverb, moved here so every amp plugin and the browser demo run one
implementation; the signal paths are the app's, expression for expression, with every transcendental in
`core::det`, so every row computes the same bits up to FMA contraction (law 10).

- **`delay::TapeDelay`** — a tape head, not a sampler: TIME is a target the head glides to (~150 ms) and the glide
  bends the repeats' pitch; every pass through the loop goes through a one-pole low pass (DARK) and a fixed `tanh`
  saturator; OFFSET holds one channel's wet back by up to ±31 ms. The dry is never touched.
- **`delay/DelayLaw.h`** — the amps' knob law: the 13-step sync ladder (`kDivisionNames`, `kDivisionBeats`,
  `divisionMs`, the tempo floored at 40 BPM so no division outruns the 6.1 s line), the sync policy
  (`timeMsFor`), the DARK knob's log-even map from 16 kHz to 800 Hz (`darkHzForPercent`) and the knob ranges.
- **`reverb::AmpReverb`** — the six characters (Ambience, Room, Hall, Plate, Spring, Modulated with its slow
  chorus), DECAY ×0.5 … ×2 through the room size, PREDELAY 0 … 100 ms and a wet-only HPF; the tail is ADDED to an
  untouched dry. Built on **`reverb::Freeverb`** — Jezar's Freeverb (2000, public domain), written here from the
  published algorithm; the behaviour of the framework Freeverb the amps were voiced on was established by running
  it, not by reading it into this code.
- **Contract**: `prepare (rate, maxBlock, maxChannels)` allocates (1 or 2 channels, 1 … 768 kHz) and a refused
  one DISARMS; `[[nodiscard]] bool process (io, nch, n, feed)`; `reset()` is a stream restart (clocks re-anchored,
  glides landed); knobs are bounded before any arithmetic; both stages in the law-11 census and the law-11c pause
  suite. Law 8: flushed per sample at the house 1e-15 (+0.008 % RT); the delay's heads glide in double and land.
- **Law 11, a new census property for every stage**: P8, a refused re-prepare disarms. Every existing stage
  already did.
- **Distance to the app** (an out-of-tree harness renders the app's own classes beside these): the delay within
  3.0e-7 at fixed knobs (ulp-level: `det` against the system libm); the room within 1.2e-10 for the unmodulated characters and 5.4e-5 for MODULATED
  (both stages primed past the app's own first-10-ms leak); `Freeverb` against the framework's within 2.8e-14. Regression goldens
  of the core's own renders (25 KB) guard the sound in CI.
- **Intel builds of the plugins move by up to 1.04e-6 (-120 dBFS) from today's app**: the framework's Freeverb
  quantises its state on x86-64 and not on arm64; these modules quantise nowhere — the arm64 sound everywhere.
- **Other departures from the app, on purpose**: a fresh room starts at its parameters (the app leaked up to
  0.8 × mix of the dry into the first 10 ms after its first prepare); the delay's heads land on their targets (the
  app's float motor parked up to 2.3 ms short at 6 s); `reset()` is a true restart; absurd knob values can no
  longer reach a NaN read position; the tail estimates are right from construction; a call longer than `maxBlock`
  is processed whole.
