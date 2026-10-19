### delay · reverb — two new modules: the amp plugins' echo and room, moved out of the product

`felitronics::delay` and `felitronics::reverb` — header-only, JUCE-free; `delay` depends on `core`, `reverb` on
`core` and `eq`. They are OrbitAmp's delay and reverb, moved here so every amp plugin and the browser demo run one
implementation; the arithmetic is the app's, expression for expression.

- **`delay::TapeDelay`** — a tape head, not a sampler: TIME is a target the head glides to (~150 ms) and the glide
  bends the repeats' pitch; every pass through the loop goes through a one-pole low pass (DARK) and a fixed `tanh`
  saturator; OFFSET holds one channel's wet back by up to ±31 ms. The dry is never touched. With it, the amps'
  parameter law: the 13-step sync ladder (`kDivisionNames`, `kDivisionBeats`, `divisionMs`) and the DARK knob's
  log-even map from 16 kHz to 800 Hz (`darkHzForPercent`).
- **`reverb::Freeverb`** — Jezar's Freeverb (2000, public domain), written here from the published algorithm and
  matched by measurement to the behaviour of the Freeverb the amps were voiced on: the published tunings scaled in
  integer arithmetic, its gain, damping, room, wet, dry, width and freeze scaling, and its 10 ms linear parameter
  glides. Against that implementation directly: 144 of 144 renders bit-identical on arm64 (random parameters with
  freeze, mono, stereo and alternating, 22.05 … 192 kHz, a 1e-7 input).
- **`reverb::AmpReverb`** — the six characters (Ambience, Room, Hall, Plate, Spring, Modulated with its slow
  chorus), DECAY ×0.5 … ×2 through the room size, PREDELAY 0 … 100 ms and a wet-only HPF; the tail is ADDED to
  an untouched dry.
- **Parity** with the app's own classes, by an out-of-tree harness (impulse, noise, sweep and their silent tails;
  44.1/48/96 kHz; blocks 1, 64, 512, 333; mono and stereo; knob moves, stand-by, restarts, one-ulp jitter): delay
  1584 renders bit-identical above 1e-30, reverb 1824 renders within 4.9e-16 (1640 bit-identical); with the law-8
  flushes removed, every render of both is bit-identical. Excerpts of the app's renders are kept as goldens
  (25 KB), so CI guards the sound without the app.
- **Contract**: `prepare (rate, maxBlock, maxChannels)` allocates (1 or 2 channels; wider is refused);
  `[[nodiscard]] bool process (io, nch, n, feed)`; idempotent `reset()`; both stages registered in the law-11
  census. Law 8: the loops are flushed per sample at 1e-30 — not the house 1e-15, because these loops integrate
  what they flush (1.2e-10 on a HALL impulse at 1e-15); written up in `docs/LAW8-AUDIT.md`.
- **Departures from the app, on purpose**: the flushes; non-finite knob values are ignored; a call longer than
  `maxBlock` is processed whole (the app's reverb dropped the rest); a channel that stops is cleared and a
  zero-width call is a pause (law 11); the delay's `prepare` clears its loop filter. And on x86-64 the app's own
  Freeverb quantises its state (it does not on arm64): this module does not, so it is the arm64 sound everywhere,
  1.04e-6 from the x86-64 app.
