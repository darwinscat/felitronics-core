<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::delay`

The amp plugins' echo — `TapeDelay` — and the parameter law their faces drive it with. Header-only,
JUCE-free, depends on `core` only. Moved here from OrbitAmp's `DelayStage` so every amp plugin and the
browser demo run one implementation.

```cmake
target_link_libraries(app PRIVATE felitronics::delay)
```

## `TapeDelay`

A tape head, not a sampler: TIME is a target the read head glides to (~150 ms), and the glide bends the
repeats' pitch. Every pass through the loop — the first echo included — goes through a one-pole low pass
(DARK) and a fixed `tanh` saturator, so repeats darken and compress as they recirculate. OFFSET holds one
channel's wet back by up to ±31 ms (its own ~30 ms glide). The dry is never touched; the wet is added.

```cpp
felitronics::delay::TapeDelay d;
d.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels; allocates the lines
d.setParams ({ .timeMs = 350, .repeats = 0.35f, .darkHz = 4000, .offsetMs = 0, .mix = 0.25f });
d.process (io, nch, n, feed);                    // [[nodiscard]] bool, in place; feed=false = stand-by
d.reset();                                       // stream restart, idempotent (cheap to call every block)
d.tailSeconds(); d.shownTimeMs(); d.shownOffsetMs(); d.inputPeak();   // relaxed atomics, any thread
```

The setters (`setTimeMs`, `setRepeats`, `setDarkHz`, `setOffsetMs`, `setMix`) are public too, each guarded
against re-sending the value it already holds; `setParams` calls them in that order.

**The parameter law** of the shipped amps, so a host lands on the app's numbers: `kDivisionNames` /
`kDivisionBeats` (13 divisions, 1/1 … 1/16T), `divisionMs (division, bpm)`, `darkHzForPercent (pct)`
(16 kHz → 800 Hz, log-even), and the knob ranges (`kTimeMinMs` … `kDarkDefaultPct`).

## Contract

- `prepare` allocates (6.1 s + 31 ms of line per channel); `process`/`reset`/setters never allocate,
  lock or throw. A refused `prepare` leaves the instance as it was.
- Law 11: an unprepared stage, a width above the prepared one, or a negative extent is refused before
  anything moves; any `n` is processed in full. A channel at an index `>= nch` has stopped — its line and
  loop filter are dropped once, and its offset head snaps to its target on return. `nch == 0` is a pause:
  the shared motor glides on as through silence.
- Law 8: the loop filter's state and the offset glide are flushed below 1e-30 **per sample** (not the
  house 1e-15: the loop integrates what it flushes — see the header). Bit-identical under any re-slicing.
- Poison: a non-finite input sample is not recorded (the dry keeps it); a non-finite loop state restarts
  the line.
- libm: `std::tanh` (per sample), `std::exp`, `std::pow`, `std::log` are the system functions on purpose —
  they are the app's (`retain-rt` in `tools/lint/det-math-manifest.txt`).

## Parity with the app

An out-of-tree harness renders OrbitAmp's own `DelayStage` (feat/titan) and this class side by side —
impulse, noise burst, log sweep, each into a silent tail; 44.1/48/96 kHz; blocks 1, 64, 512, 333; mono and
stereo; six free settings, all 13 divisions at 120 BPM plus 1/1 at 40 and 1/16T at 240; knob moves,
stand-by, restarts, and one-ulp knob jitter. 1584 renders on arm64: every sample above 1e-30 is
bit-identical, and with the law-8 flush removed all 1584 are bit-identical. The picture taps agree exactly.
`tests/delay_golden.h` keeps excerpts of the app's renders so CI guards this without the app.

**Where it departs from the app, on purpose:** the law-8 flush above; a non-finite knob value is ignored
(the app's clamps passed a NaN through); `setDarkHz` is clamped to 1 Hz … 100 kHz (the app's mapping never
leaves 800 … 16000); `prepare` zeroes the loop filter (the app kept the previous session's value across a
re-prepare); and law 11's stopped-channel and pause rules (the app froze a channel it stopped feeding).
