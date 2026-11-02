<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::reverb`

The amp plugins' six-character room, `AmpReverb`, and the Freeverb it is built on. Header-only, JUCE-free,
depends on `core` and `eq` (the wet HPF). Moved here from OrbitAmp's `ReverbStage`, which sat on the
framework's Freeverb and so had to live in the product; now every amp plugin and the browser demo run one
implementation. **A host uses `AmpReverb`**; `Freeverb` is the primitive it owns.

```cmake
target_link_libraries(app PRIVATE felitronics::reverb)
```

## `AmpReverb`

Mix is the hero, the CHARACTER is the title (Ambience, Room, Hall, Plate, Spring, Modulated — each a room
size, damping and width; Modulated adds a slow 0.45 Hz chorus on the tail), and DECAY (×0.5 … ×2, through
the room size), PREDELAY (0 … 100 ms) and the wet-only HPF are refinements. The dry is never touched; the
whole wet chain runs on a scratch copy and is added.

```cpp
felitronics::reverb::AmpReverb r;
r.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 1 … 768 kHz
r.setParams ({ .character = Character::plate, .mix = 0.3f, .decay = 1, .predelayMs = 0, .hpfHz = 120 });
r.process (io, nch, n, feed);                    // [[nodiscard]] bool, in place; feed=false = stand-by
r.reset();                                       // stream restart, idempotent
r.tailSeconds(); r.addedWet (ch); r.addedLength();
```

`kCharacterNames` carries the six labels. The per-knob setters are public, guarded like the app's;
`setParams` calls them in the app's order.

## `Freeverb`

Jezar's Freeverb (Dreampoint, 2000, public domain): eight lowpass-feedback combs into four allpasses per
channel, the right channel spread 23 samples. Its constants: the published tunings scaled as
`(int) rate * tuning / 44100`; input gain 0.015; damping × 0.4; feedback `room × 0.28 + 0.7`; wet × 3 split
by width; dry × 2; freeze at ≥ 0.5; a linear 10 ms glide on damping, feedback, dry and both wet gains.
`processStereo` / `processMono` are `[[nodiscard]] bool` and refuse when unprepared.

## Contract

- `prepare` allocates; `process`/`reset`/setters never allocate, lock or throw. A refused `prepare`
  **disarms**. A prepared room starts AT its set parameters — nothing glides in.
- Law 11: refusals as in the delay; a call longer than the scratch is processed in chunks of
  `min(maxBlock, 8192)`, the same as the caller chunking (`addedWet` then holds the last chunk; a whole-file
  `maxBlock` costs only the capped scratch). One channel runs the room's mono path; a channel that stopped is
  cleared once; `nch == 0` is a pause during which the room's glides and the chorus LFO run on as through
  silence.
- `reset()` is a stream restart: every line cleared, every clock re-anchored (line positions, the chorus
  LFO's phase, the law-8 grid) and every glide landed on its target.
- `feed` fades, it does not switch: what enters the wet chain ramps linearly to nothing, or back, over
  `kFeedRampSeconds` (5 ms), so a flip on a loud note sends no click into the room. The room's own tail is
  never ramped, and a held `feed` is bit-identical to the old switch. A pause moves the ramp as silence would;
  after `prepare`/`reset` it starts settled on the next call's `feed` — a fresh room does not fade in.
- Law 8: Freeverb flushes its comb filters and every comb/allpass write below 1e-15 per sample
  (`core::flushDenormal` — its buffers are the state); the HPF rides `core::StateGrid` with the eq module's
  own flush. A long silence ends in exact zero, and stays there.
- Knobs: non-finite values are ignored, an unknown character index too; finite ones are clamped.
- Threads: prepare, reset, the setters and process belong to one thread; `tailSeconds()` may be read from
  any thread; `addedWet`/`addedLength` belong to the process thread.
- Math: `core::det` (`log2` for DECAY and the tail, `sin` for the chorus LFO — per sample, ~5 ns), so the
  same bits on every row up to FMA contraction (law 10).

## Provenance and parity

`AmpReverb` is OrbitAmp's own `ReverbStage` (feat/titan, by the same authors), expression for expression.
`Freeverb` was written here from the published, public-domain description of Jezar's algorithm and its
well-known constants. The framework's Freeverb (JUCE's `juce::Reverb`), on which the app's stage ran, was
only ever RUN, never copied: its behaviour — the integer tuning scale, the glide, the mono path leaving wet2
alone, no state quantisation on arm64 — was established by measuring its output against this one's.
`tests/reverb_golden.h` holds regression excerpts of this module's own renders, so CI guards the sound.

Measured on arm64 against the app's class (impulse, noise, sweep and their silent tails; 44.1/48/96 kHz;
blocks 1, 64, 512, 333; mono and stereo; every character at mix/decay/predelay/HPF corners and middles;
character walks, stand-by and restarts, a per-block DECAY sweep, one-ulp jitter). Both stages are primed
after `prepare` with floor(0.01 · rate) samples of silence — the app's first prepare glides in from
Freeverb's published initial parameters, and through silence that glide carries nothing — and compared from
sample 0: **892 of 1728 fixed-knob renders bit-identical; the five unmodulated characters within 1.2e-10**
(the flush and `det::log2`), **MODULATED within 5.4e-5** (its LFO in `det::sin` and its read position
pinned against contraction, where the app's is the system `sin` and a fused multiply-add). With the knobs
moving: the DECAY sweep and the one-ulp jitter bit-identical, the character walk 3.5e-5 and stand-by with
restarts 2.9e-5 (both through MODULATED). `Freeverb` against the framework's
directly (random parameters with freeze, dry and width; mono, stereo, alternating; parameters set before
the first prepare; 22.05 … 192 kHz; a 1e-7 input): 73 of 144 bit-identical, the rest within 2.8e-14 —
the law-8 flush.

**Intel builds of the plugins move.** Built for x86-64 (measured under Rosetta), the framework's Freeverb
quantises its comb and allpass state to 0.1f's ulp (≈ 7.5e-9); on arm64 it does not. This module does not
quantise anywhere — it is the arm64 sound on every platform (owner's decision) — so on Intel the plugins
move by up to 1.04e-6 (-120 dBFS) from today's app.

**Where it departs from the app, on purpose:** the Intel line above; `core::det` instead of the system libm;
the law-8 flushes; a fresh room no longer glides in from Freeverb's published initial parameters (the app's
first 10 ms after its first prepare added up to 0.8 × mix of the dry); `reset()` re-anchors the clocks and
lands the glides (the app's chorus kept its phase across a reset); `prepare` disarms on refusal; the tail
estimate is right from construction (the app reported 1.000 s, the formula says 1.176); a call longer than
`maxBlock` is processed whole (the app dropped the rest); and law 11's stopped-channel and pause rules.
