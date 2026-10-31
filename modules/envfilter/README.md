<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::envfilter`

The bass amp's envelope filter — `EnvelopeFilter`, the funk auto-wah of the Mu-Tron III family — and the
one-slider law a face drives it with (`EnvelopeFilterLaw.h`). Header-only, JUCE-free, real-time safe, no heap
at all; depends on `core`, `dynamics` (the follower) and `eq` (the split's `Svf`).

```cmake
target_link_libraries(app PRIVATE felitronics::envfilter)
```

## `EnvelopeFilter`

A resonant filter whose cutoff follows how hard the string is hit. A peak follower (attack, release) on the
loudest channel reads the playing; its level in dB, pushed by SENSITIVITY, places the cutoff on a log-frequency
line between the range's ends — fully open when envelope + sensitivity reaches 0 dBFS, closed 30 dB below.
UP rests at the low end and a hit opens it; DOWN rests at the high end and a hit pushes it down. The cutoff
goes through one 4 ms smoother and the filter is redesigned every sample, so it never steps.

The filter is a Cytomic TPT state-variable filter (`eq::Svf`'s equations) giving low-pass and band-pass from
one state: a mode change crossfades the two over 10 ms. Both are scaled by 1/√max(1, Q), so the peak grows as
√Q in either mode. KEEP LOWS (0..1) takes everything under `keepLowsHz` out of the filter and adds it back
clean, so the note keeps its weight; the lows are `x − HP(x)`, an exact complement of the input, which is what
lets the dry/wet mix blend without notching the bass (an LR4 pair would sum to an allpass that is −1 at its
crossover).

```cpp
felitronics::envfilter::EnvelopeFilter f;
f.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 … 768 kHz
f.setParams (felitronics::envfilter::paramsForAmount (slider, menu));
f.process (io, nch, n);                          // [[nodiscard]] bool, in place
f.reset();                                       // stream restart, idempotent
f.shownCutoffHz();                               // relaxed atomic, any thread — for the picture
```

| Knob | Range | Default | |
|---|---|---|---|
| `sensitivityDb` | −24 … +36 dB | +6 | how far a given hit opens the sweep |
| `attackMs` / `releaseMs` | 0.1 … 100 / 5 … 2000 ms | 3 / 120 | the follower; the release is the "wow" |
| `sweepLoHz` / `sweepHiHz` | 20 … 12000 Hz | 150 / 3000 | the range (either order; never above 0.45 fs) |
| `direction` | up / down | up | |
| `mode` | lowPass / bandPass | lowPass | |
| `resonance` | Q 0.5 … 12 | 4 | |
| `keepLows` / `keepLowsHz` | 0 … 1 / 30 … 300 Hz | 1 / 100 | the clean lows |
| `mix` | 0 … 1 | 1 | 0 = dry, and the stage sleeps |

mix, sensitivity, resonance and keepLows glide (20 ms, and land); the range and the direction move the cutoff
through its smoother; attack/release and keepLowsHz retune at once (menu settings).

**Amount 0 is not there.** When `mix` stands at 0 the stage returns the buffer untouched — bit-exact, NaN
payloads included — and computes nothing; a move to 0 glides out (~0.41 s to land) and then sleeps. Waking
starts fresh, the mix gliding in from 0.

## `EnvelopeFilterLaw.h` — the one slider

`paramsForAmount (amount, menu)`: the slider's lower half raises MIX 0 → 1 (all filter from 0.5), and the
whole travel raises SENSITIVITY −6 → +18 dB. At a quarter the filter sits beside the dry note and only hard
notes speak (a −12 dBFS note opens 60 % of the way); at the top every fingered note opens it fully and closes
it as it decays. The mode, resonance, direction, range, follower times and clean lows are the voice — the
menu's, passed through untouched.

## Tests

`tests/EnvelopeFilterTests.cpp`: the cutoff tracks the envelope, measured in the audio with a probe tone (two
burst levels, both directions, 44.1/48/96 kHz: 12 dB more moves it 1.729 oct, the law's number; the probe's gain
matches the published cutoff's analytic response within 0.01 dB); the kernel at Q 12 with the cutoff flipping
between 20 Hz and 0.45 fs every sample stays finite (largest tap 22.6 on ±1 noise); a full-scale slap at the
fastest settings peaks at −0.55 dBFS; no zipper (energy above 6 kHz −150 dB under the whole, the input's own
−152 … −154; a cutoff held for 32 samples reads −105); keepLows 1 holds 31 … 73 Hz within 1.23 dB where 0 loses
43 dB; plus the law-11 properties, amount 0, glides, re-slicing, 11a/11c, reset, poison, law 8 and allocation.

`tools/EnvelopeFilterRender.cpp` (`felitronics_envfilter_render <outDir> [inDir]`) renders a bass DI —
every `*.wav` in `inDir`, or a synthesised plucked-string line — at several slider positions to 24-bit WAVs.
