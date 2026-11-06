<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::chorus`

The bass amp's chorus — `BassChorus` — and the one-slider law its face drives it with (`ChorusLaw.h`).
Header-only, JUCE-free, depends on `core` and `eq` (the crossover's `DetMath` Svf).

```cmake
target_link_libraries(app PRIVATE felitronics::chorus)
```

## `BassChorus`

The dry plus a VOICE: a short copy of the signal read by a head that swings on a sine around a centre
delay. The bass rule is a rule about the voice: it is made from the input's Linkwitz–Riley high band (LR4's
high half, two Butterworth high passes at the crossover). The dry is never split: it goes out sample for
sample as it came in, at every setting, so it carries no crossover phase and there is no switch when the
effect comes in. On a stereo call the two sides carry two voices whose LFOs are in quadrature (a quarter
period apart); a mono call is one voice.

```cpp
felitronics::chorus::BassChorus c;
c.setParams ({ .depthMs = 1.5f, .mix = 0.5f, .rateHz = 0.7f, .crossoverHz = 200.0f, .delayMs = 8.0f });
c.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 … 768 kHz
c.setParams (felitronics::chorus::paramsForAmount (amount, c.params()));   // the face's slider
c.setRateHz (1.2f); c.setCrossoverHz (250.0f);   // the right-click menu
c.process (io, nch, n);                          // [[nodiscard]] bool, in place
c.reset();                                       // stream restart, idempotent
c.isIdle();                                      // amount 0 computes nothing
BassChorus::tailSeconds();                       // 71 ms: the voice's tail at any setting, to -100 dB
```

- **Voice band:** LR4's high half — -6.02 dB at the crossover, -24.6 dB an octave below, -48 two octaves
  below. With the head held still, the voice plus LR4's low half of the same input sums flat to 4e-6 dB from
  20 Hz to 20 kHz. **The lows are attenuated in the voice, not excluded:** the crossover's floor is 100 Hz,
  where the low B and E fundamentals enter the voice at least 31 dB down. The Side, `(L − R)/2`, of a tone at
  `f` is at most the voice's own level there, `(f/fc)⁴ / (1 + (f/fc)⁴)` as energy re the tone's, at any rate
  and centre (sidebands included, measured): the open strings' fundamentals read -76 dB (low B) to -28 dB (G2)
  at the default 200 Hz crossover, -52 to -9 dB at the 100 Hz floor.
- **The swing:** a sine LFO (0.05 … 10 Hz), depth 0 … 10 ms of head travel around a 5 … 25 ms centre. The
  head is held 1 ms clear of the write head by a smooth minimum (a 1 ms knee), so a centre gliding into that
  limit bends the head instead of kinking its speed. The centre's 5 ms floor is the law's full depth plus that
  1 ms plus the knee: the slider's full 3 ms is never trimmed. The frequency ratio swings by ±2π·rate·depth; the flat excursion is
  the larger in cents.
- **Quadrature, not opposite phases:** opposite voices spread wider, but their mono sum is a standing comb that
  breathes like a flanger; quadrature voices still swing when a PA or a DI sums the pair. An adapter that wants
  two voices from a mono source puts it in two DISTINCT buffers.
- **The read:** 4-point cubic Hermite (Catmull–Rom), in double. At its worst fraction it loses 0.04 dB at
  5 kHz and 0.53 dB at 10 kHz (48 kHz); a linear read loses 0.47 and 2.0 dB, and its loss swings with the LFO.
- **Level:** the voice is ADDED, so the voiced band rises — about +3 dB at amount 1 for uncorrelated content,
  up to +6 dB on coherent peaks — while the lows are unchanged. Budget that headroom before a limiter.

## The one-slider law (`ChorusLaw.h`)

`paramsForAmount (amount, menu)` owns two fields, `depthMs = 3 · amount` and `mix = amount`: the wobble and
the voice's level rise together, linearly. Rate, crossover and centre delay are the menu's and pass through
untouched — a menu RATE then means the same at every slider position, and the wobble scales with it as a
pedal's does. The menu must not expose depth or mix. A non-finite amount is ignored on the condition the call
keeps: the caller passes the stage's current `params()` as the menu, so `menu` comes back untouched and nothing
changes (a menu kept apart, or the default, would bring its own depth and mix). A finite one is clamped to 0..1.
`pitchDeviationCents (params)` is the face's read-out: a pure function, the flat excursion in cents (23.0 at
amount 1 and 0.7 Hz).

**Amount 0 is mix 0, and at mix 0 the stage idles:** the voice stops on the very sample the counted 20 ms mix
ramp lands on an exact zero — the rest of that call is the input's own — and from then on `process` returns
with the buffer untouched and computes nothing (`voicedSamples()` is the odometer that shows it). **A wake is a restart** — lines,
filters, LFO, write clock and law-8 grid back at zero, every knob on its target — and the voice opens through
its 10 ms door.

## Contract

- `prepare` allocates the lines (35 ms + reach per channel, a power of two); `process`/`reset`/setters never
  allocate, lock or throw. A refused `prepare` **disarms** until a prepare succeeds (law 11b).
- Law 11: a negative extent, an unprepared stage or a width above the prepared one is refused before anything
  moves; any `n` is processed in full. A channel at an index `>= nch` has stopped — its line, filter state and
  door are dropped once (11a), so it returns as a new stream that fades in, while the other lane and the LFO
  carry on. `nch == 0` is a pause: the LFO, the glides, the mix ramp and the grid move as through silence
  (11c).
- `reset()` is a stream restart. Bit-identical under any re-slicing (law 8a).
- Glides: the head's depth and centre glide critically damped — two one-poles in cascade, 75 and 80 ms each —
  so the head's speed never steps. Rate (50 ms) and crossover (30 ms, coefficients per sample) are one-poles;
  the mix is a counted linear ramp.
- Law 8: the Svf's flush on `core::StateGrid` (the eq module's contract), the line flushed per write, the
  glides land. The law-8 test sees the lines and the glides; the Svf's integrators stand on eq's own tests.
- Poison: a non-finite input sample is not voiced (the dry keeps it); an overflowing filter writes silence
  and is healed at the end of the call; the double read bounded to the float range turns a near-FLT_MAX burst
  into ±inf at worst, never NaN. Knobs: a non-finite value is ignored, a finite one clamped.
- Threads: one thread owns the stage; the face reads the law's pure functions. Math: `core::det` (the LFO's
  sin and cos, the glides' exp, the read-out's log2).

## Listening

`tools/ChorusRender.cpp` (target `felitronics_chorus_render`, built with the tests) runs a bass DI — the WAVs
it is given, or a synthesised Karplus–Strong line — through the stage at amounts 0.25 / 0.5 / 0.75 / 1, at
0.75 with a 500 Hz crossover, and amount 1 folded to mono, and writes 24-bit stereo WAVs beside the dry.
