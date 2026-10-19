<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::reverb`

Freeverb, and the amp plugins' six-character room built on it. Header-only, JUCE-free, depends on `core`
and `eq` (the wet HPF). Moved here from OrbitAmp's `ReverbStage`, which sat on the framework's Freeverb and
so had to live in the product; now every amp plugin and the browser demo run one implementation.

```cmake
target_link_libraries(app PRIVATE felitronics::reverb)
```

## `Freeverb`

Jezar's Freeverb (Dreampoint, 2000, public domain) written here from the published algorithm: eight
lowpass-feedback combs into four allpasses per channel, the right channel spread 23 samples. It reproduces
the behaviour of the Freeverb the amps were voiced on (the one JUCE ships as `juce::Reverb`), established
by measuring that implementation's output, not by copying it: the published tunings scaled as
`(int) rate * tuning / 44100`; input gain 0.015; damping × 0.4; feedback `room × 0.28 + 0.7`; wet × 3 split
by width; dry × 2; freeze at ≥ 0.5; and a linear 10 ms glide on damping, feedback, dry and both wet gains.

```cpp
felitronics::reverb::Freeverb f;
f.prepare (sampleRate);                 // [[nodiscard]] bool; allocates; lands every glide
f.setParams ({ .roomSize = 0.5f, .damping = 0.5f, .wetLevel = 0.33f, .dryLevel = 0.4f, .width = 1, .freezeMode = 0 });
f.processStereo (left, right, n); f.processMono (x, n);
f.reset(); f.clearChannel (ch); f.advanceGlides (n);
```

## `AmpReverb`

Mix is the hero, the CHARACTER is the title (Ambience, Room, Hall, Plate, Spring, Modulated — each a room
size, damping and width; Modulated adds a slow 0.45 Hz chorus on the tail), and DECAY (×0.5 … ×2, through
the room size), PREDELAY (0 … 100 ms) and the wet-only HPF are refinements. The dry is never touched; the
whole wet chain runs on a scratch copy and is added.

```cpp
felitronics::reverb::AmpReverb r;
r.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels
r.setParams ({ .character = Character::plate, .mix = 0.3f, .decay = 1, .predelayMs = 0, .hpfHz = 120 });
r.process (io, nch, n, feed);                    // [[nodiscard]] bool, in place; feed=false = stand-by
r.reset();                                       // stream restart, idempotent
r.tailSeconds(); r.addedWet (ch); r.addedLength();
```

`kCharacterNames` carries the six labels. The per-knob setters are public, guarded like the app's;
`setParams` calls them in the app's order.

## Contract

- `prepare` allocates; `process`/`reset`/setters never allocate, lock or throw.
- Law 11: refusals as in the delay; a call longer than `maxBlock` is processed in `maxBlock` chunks, the
  same as the caller chunking (`addedWet` then holds the last chunk). One channel runs the room's mono
  path; a channel that stopped is cleared once; `nch == 0` is a pause during which the glides and the
  chorus LFO run on.
- Law 8: Freeverb flushes its comb filters and every comb/allpass write below 1e-30 **per sample** (its
  buffers are the state; see the header for why not 1e-15); the HPF rides `core::StateGrid` with the eq
  module's own flush. A long silence ends in exact zero.
- Poison: a non-finite input sample never enters the wet chain; a non-finite tail is dropped whole and the
  room restarts.
- libm: `std::sin` (the chorus LFO, per sample), `std::log2`, `std::log` are the system functions on
  purpose — they are the app's (`retain-rt` in the manifest).

## Parity with the app

An out-of-tree harness renders OrbitAmp's own `ReverbStage` (feat/titan, on `juce::Reverb` of JUCE 8.0.16) and
`AmpReverb` side by side — impulse, noise burst, log sweep into a silent tail; 44.1/48/96 kHz; blocks 1,
64, 512, 333; mono and stereo; every character at mix/decay/predelay/HPF corners and middles; a character
walk, stand-by and restarts, a per-block DECAY sweep (each move restarting the room's glide), one-ulp
jitter. 1824 renders on arm64: 1640 bit-identical, the rest within 4.9e-16 — the HPF's 1e-15 flush, and
nothing else: with the flushes removed all 1824 are bit-identical. `Freeverb` against `juce::Reverb`
directly (random parameters with freeze, dry and width; mono, stereo, alternating; parameters set before
the first prepare; 22.05 … 192 kHz; a 1e-7 input): 144 of 144 bit-identical. `tests/reverb_golden.h`
keeps excerpts of the app's renders so CI guards this without JUCE.

**x86-64 is different, in the app itself.** Built for x86-64 (measured under Rosetta), the framework's
Freeverb quantises its comb and allpass state to 0.1f's ulp (≈ 7.5e-9); on arm64 it does not. So the
shipped amps already differ between the two by up to 1.04e-6 (-120 dBFS) — the harness's x86-64 build
reads exactly that against this module, and 0 for the delay. This module does not snap anywhere: it is the arm64 sound on every
platform, and on x86-64 it differs from the x86-64 app by that same amount.

**Where it departs from the app, on purpose:** the law-8 flushes; non-finite knob values are ignored;
an unknown character index is ignored; a call past `maxBlock` is processed whole (the app dropped the
tail); and law 11's stopped-channel and pause rules.
