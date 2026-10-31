<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::envfilter`

The bass amp's envelope filter — `EnvelopeFilter`, the funk auto-wah of the Mu-Tron III family — and the
one-slider law a face drives it with (`EnvelopeFilterLaw.h`). Header-only, JUCE-free, real-time safe, no heap
at all; depends on `core`, `dynamics` (the follower) and `eq` (the splits).

```cmake
target_link_libraries(app PRIVATE felitronics::envfilter)
```

## `EnvelopeFilter`

A resonant filter whose cutoff follows how hard the string is hit. A peak follower (attack, release) on the
loudest channel — or on an external KEY, e.g. the signal before the amp's compressor and octaver — reads the
playing; its level in dB, pushed by SENSITIVITY, places the cutoff on a log-frequency line between the range's
ends — fully open when envelope + sensitivity reaches 0 dBFS, closed 30 dB below. UP rests at the low end and a
hit opens it; DOWN rests at the high end and a hit pushes it down. The playing reaches the cutoff through a 4 ms
smoother; the hand's edits (direction, range, sensitivity) glide over 20 ms first, so a direction flip is a
sweep, not a hit. The filter is redesigned every sample, so it never steps.

The filter is a Cytomic TPT state-variable filter (`eq::Svf`'s equations), low-pass or band-pass, a mode change
crossfading the two over 10 ms. Both are scaled by 1/√max(1, Q): at the cutoff they are equal at every Q and the
peak grows as √Q; the price is the low pass's passband, **−6.0 dB at Q 4 and −10.8 dB at Q 12** (pinned by the
tests). KEEP LOWS (0..1) takes the band under `keepLowsHz` out of the filter and adds it back clean, so the note
keeps its weight — and each mode has the split that cannot cancel it:

- **low pass:** the lows are `x − HP(x)`, an exact complement of the input, and the dry is `x`. The low pass's
  passband is in phase with `x`, so the two only add: over a grid of cutoff × Q × corner × mix the band under the
  corner never drops more than 0.07 dB (the complement lifts it at most +1.25 dB).
- **band pass:** an LR4 pair, whose halves are in phase at every frequency. A band pass never has a negative real
  part, so it can only add to the low half: the band under the corner keeps at least the LR4 low share (−0.5 dB
  at half the corner, −1.9 dB at 0.7, −6 dB at the corner). The dry is the same LR4 allpass, so a mix never
  notches. (The complement here let a band pass parked just above the corner take a 98 Hz note down 13.6 dB.)

The default corner, 150 Hz, sits a fifth above G2, so every open-string-to-G2 fundamental keeps ≥ −1.45 dB in
band pass with the filter wide open.

```cpp
felitronics::envfilter::EnvelopeFilter f;
f.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 … 768 kHz
f.setParams (felitronics::envfilter::paramsForAmount (slider, menu));
f.process (io, nch, n);                          // [[nodiscard]] bool, in place, self-keyed
f.process (io, nch, n, key, numKeyChannels);     // ...keyed (nullptr / 0 = self-keyed)
f.reset();                                       // stream restart, idempotent
f.shownCutoffHz();                               // relaxed atomic, any thread — for the picture
```

| Knob | Range | Default | |
|---|---|---|---|
| `sensitivityDb` | −24 … +36 dB | +6 | **law-owned** — how far a given hit opens the sweep |
| `attackMs` / `releaseMs` | 0.1 … 100 / 5 … 2000 ms | 3 / 120 | the follower; the release is the "wow" |
| `sweepLoHz` / `sweepHiHz` | 20 … 12000 Hz | 150 / 3000 | the range (either order; never above 0.45 fs) |
| `direction` | up / down | up | |
| `mode` | lowPass / bandPass | lowPass | |
| `resonance` | Q 0.5 … 12 | 4 | |
| `keepLows` / `keepLowsHz` | 0 … 1 / 30 … 300 Hz | 1 / 150 | the clean lows |
| `mix` | 0 … 1 | 1 | **law-owned** — 0 = dry, and the stage goes idle |

Everything glides (20 ms, one-pole, and lands) except the follower's attack/release, which retune at once; the
mode door and the engage door are counted 10 ms ramps.

**Amount 0 is not there.** Setting `mix` to 0 closes the engage door (10 ms, counted, lands exactly); from then
on the stage is idle (`isIdle()`): it returns the buffer untouched — bit-exact, NaN payloads included — and
computes nothing. A wake is a restart: every state from silence, the door opening from dry, the mix gliding in.

Cost, measured on deb (gcc 14, i9-13900H): 1.0 % of real time at 48 kHz stereo. Most of it is the per-sample
sweep's `core::det` log2/exp2/tan (the same bits on every row); with the system's three it is 0.53 %.

## `EnvelopeFilterLaw.h` — the one slider

`paramsForAmount (amount, menu)`: the slider's lower half raises MIX 0 → 1 (all filter from 0.5), and the
whole travel raises SENSITIVITY −6 → +18 dB. At a quarter the filter sits beside the dry note and only hard
notes speak (a −12 dBFS note opens 60 % of the way); at the top every attack hotter than −18 dBFS opens it fully
and the note's decay closes it. `mix` and `sensitivityDb` are the law's — a face must not put them in its menu;
the mode, resonance, direction, range, follower times and clean lows are the voice, passed through from `menu`.
RESONANCE is a menu knob, but its audible effect grows with the slider: the mix brings in the low pass's
1/√Q passband, so the slider's lower half is also a level change on the harmonics. A non-finite amount returns
`menu` unchanged; a finite one is clamped to 0..1.

## Tests

`tests/EnvelopeFilterTests.cpp` — the effect: the cutoff tracks the envelope, measured in the audio with a probe
tone (two burst levels, both directions, 44.1/48/96 kHz: 12 dB more moves it 1.729 oct, the law's number; the
probe's gain matches the published cutoff's response within 0.01 dB) and follows an external key; the kernel at
Q 12 with the cutoff flipping between 20 Hz and 0.45 fs every sample stays finite; a full-scale slap at the
fastest settings peaks at −1.9 dBFS; no zipper (−150 dB above 6 kHz); the clean lows over 360 cases of cutoff ×
Q × corner × mix × tone; the passband and mode levels. The stage: the law-11 properties, amount 0, slider and menu
moves without a step, the hand's edits as sweeps (≥ 30 ms), re-slicing, 11a/11c (keyed too), reset, poison (an
overflow healed at the end of its call), law 8, the 0.45 fs cap, allocation.

`tools/EnvelopeFilterRender.cpp` (`felitronics_envfilter_render <outDir> [inDir]`) renders a bass DI —
every `*.wav` in `inDir`, or a synthesised plucked-string line — at several slider positions to 24-bit WAVs.
