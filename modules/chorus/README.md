<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::chorus`

The bass amp's chorus — `BassChorus` — and the one-slider law its face drives it with. Header-only,
JUCE-free, depends on `core` and `eq` (the crossover's `DetMath` Svf).

```cmake
target_link_libraries(app PRIVATE felitronics::chorus)
```

## `BassChorus`

The dry plus a VOICE: a short copy of the signal read by a head that swings on a sine around a centre
delay. The bass rule is a rule about the voice: it is made from the input's Linkwitz–Riley high band (LR4's
high half, two Butterworth high passes at the crossover), so nothing below the crossover is ever delayed,
modulated or spread — the lows stay the dry's, mono and in phase. The dry itself is never split: it goes
out sample for sample as it came in, at every setting, so it carries no crossover phase and there is no
switch when the effect comes in. On a stereo call the two sides carry two voices whose LFOs swing in
opposite phase (one sharp while the other is flat); a mono call is one voice.

```cpp
felitronics::chorus::BassChorus c;
c.setParams ({ .depthMs = 1.5f, .mix = 0.5f, .rateHz = 0.7f, .crossoverHz = 200.0f, .delayMs = 8.0f });
c.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 … 768 kHz
c.setAmount (0.6f);                              // the face's slider: depth and mix from the law
c.setRateHz (1.2f); c.setCrossoverHz (250.0f);   // the right-click menu
c.process (io, nch, n);                          // [[nodiscard]] bool, in place
c.reset();                                       // stream restart, idempotent
c.isIdle(); c.pitchDeviationCents();             // amount 0 computes nothing; the set wobble, in cents
```

- **Voice band:** LR4's high half. Measured: -6.02 dB at the crossover, -24.6 dB an octave below; with the
  head held still, the voice plus LR4's low half of the same input is flat to 4e-6 dB from 20 Hz to 20 kHz.
- **The swing:** a sine LFO (0.05 … 10 Hz), depth 0 … 10 ms of head travel around a 2 … 25 ms centre, held
  1 ms clear of the write head. Peak pitch deviation `2π · rate · depth` — 22.7 cents at the law's full
  3 ms and the 0.7 Hz default. A sine, not a triangle: a triangle's corners are pitch steps.
- **The read:** 4-point cubic Hermite (Catmull–Rom). At its worst fraction it loses 0.04 dB at 5 kHz and
  0.53 dB at 10 kHz (48 kHz); a linear read loses 0.6 and 2.0 dB, and its loss swings with the LFO.
- **Mix:** the voice is ADDED; 1 is the voice at the dry's level (the classic 50/50 chorus).

## The one-slider law (`law::paramsFor`, `BassChorus::setAmount`)

`amount` 0..1 → `depthMs = 3 · amount`, `mix = amount`: the wobble and the voice's level rise together,
linearly, so the pitch deviation (∝ depth) and the voice's level both read evenly along the slider. Rate,
crossover and centre delay are the menu's and the slider never moves them — a menu RATE then means the same
thing at every slider position, and the wobble scales with it as a pedal's does (faster is more pitch).
**Amount 0 is mix 0, and at mix 0 the stage idles:** once the mix has glided onto an exact zero, `process`
returns with the buffer untouched and computes nothing. Pulled off zero, the voice restarts from silence (its
line and filters cleared, its LFO at phase 0) while the mix glides up from 0.

## Contract

- `prepare` allocates the lines (35 ms + reach per channel, a power of two); `process`/`reset`/setters never
  allocate, lock or throw. A refused `prepare` **disarms** until a prepare succeeds (law 11b).
- Law 11: a negative extent, an unprepared stage or a width above the prepared one is refused before anything
  moves; any `n` is processed in full. A channel at an index `>= nch` has stopped — its line and filter state
  are dropped once (11a). `nch == 0` is a pause: the LFO, the glides, the door and the grid move as through
  silence (11c).
- `reset()` is a stream restart: lines and filters cleared, the write clock, LFO and grid at zero, every glide
  on its target. Bit-identical under any re-slicing (law 8a).
- Glides: the head (depth, centre) and the mix glide critically damped — two one-poles in cascade, 40 / 80 /
  10 ms each — so the head's speed never steps (a one-pole's first sample is a pitch jump: 6 % for a 3 ms depth
  at 50 ms). Rate (50 ms) and crossover (30 ms, coefficients per sample) are one-poles. After `prepare`, `reset`
  or a wake the voice's input opens over 10 ms.
- Law 8: the Svf's flush on `core::StateGrid`, the line flushed per write, the glides land.
- Poison: a non-finite input sample is not voiced (the dry keeps it); an overflowing filter writes silence
  and is healed at the end of the call. Knobs: a non-finite value is ignored, a finite one clamped.
- Threads: one thread owns the stage. Math: `core::det` (the LFO's sin, the glides' exp, the published log2).

## Listening

`tools/ChorusRender.cpp` (target `felitronics_chorus_render`, built with the tests) runs a bass DI — the WAVs
it is given, or a synthesised Karplus–Strong line — through the stage at amounts 0.25 / 0.5 / 0.75 / 1 and at
0.75 with a 500 Hz crossover, and writes 24-bit stereo WAVs with the dry beside them.
