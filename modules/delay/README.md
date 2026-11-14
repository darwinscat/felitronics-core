<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::delay`

The amp plugins' echo — `TapeDelay` — and the knob law their faces drive it with (`DelayLaw.h`).
Header-only, JUCE-free, depends on `core` only. Moved here from OrbitAmp's `DelayStage` so every amp
plugin and the browser demo run one implementation.

```cmake
target_link_libraries(app PRIVATE felitronics::delay)
```

## `TapeDelay`

A tape head, not a sampler: TIME is a target the motor glides to (a ~150 ms time constant), and the glide
bends the repeats' pitch. Every pass through the loop — the first echo included — goes through a one-pole low
pass (DARK) and a fixed `tanh` saturator, so repeats darken and compress as they recirculate. OFFSET holds one
channel's wet back by up to ±31 ms (its own ~30 ms glide). The dry is never touched; the wet is added.

**Up to `kMaxHeads` (4) playback heads** read the one tape. Each has a POSITION, a fraction of TIME in (0, 1]
— TIME stays the one motor, so turning it moves every head, and a head moved on its own glides under TIME's law
and bends that repeat — a LEVEL (0 … 1, 0 is off) and a PAN. The loop feeds back from 1x TIME, never from a
head, so the heads' figure repeats once per TIME and decays with REPEATS: heads at 1/3, 2/3 and 1 play a
triplet, at 1/4 … 1 four on the beat, and the third echo can be louder than the second. PAN is a balance on the
head's own lines at constant power, the centre at unity on both sides (+3 dB on its side at hard pan); on dual
mono — one guitar on both lines, as the amps feed it — heads panned apart are a ping-pong. A mono call aims
every head at its level, pan aside. Level and pan arrive over 20 ms ramps — linear in gain, so constant power is
where they land and a jump across the field dips up to 3 dB on the way; a silent head stands where it was sent. The default layout, one head
at 1x, level 1, centre, is the one-head stage byte for byte — every existing preset plays as it did.

```cpp
felitronics::delay::TapeDelay d;
d.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 1 … 768 kHz
felitronics::delay::TapeDelayParams p { .timeMs = 350, .repeats = 0.35f, .darkHz = 4000, .offsetMs = 0, .mix = 0.25f };
p.heads[0] = { .position = 1.0f / 3, .level = 0.4f, .pan = -1.0f };   // a ping-pong triplet
p.heads[1] = { .position = 2.0f / 3, .level = 0.7f, .pan = +1.0f };
p.heads[2] = { .position = 1.0f,     .level = 1.0f, .pan =  0.0f };
d.setParams (p);                                 // the whole state, heads included
d.setHead (1, { .position = 0.5f, .level = 0.7f, .pan = 1.0f });   // or one head (an index off the array is ignored)
d.process (io, nch, n, feed);                    // [[nodiscard]] bool, in place; feed=false = stand-by
d.reset();                                       // stream restart, idempotent (cheap to call every block)
d.tailSeconds(); d.shownTimeMs(); d.shownOffsetMs(); d.shownHeadMs (k); d.inputPeak();   // relaxed atomics, any thread
```

The setters (`setTimeMs`, `setRepeats`, `setDarkHz`, `setOffsetMs`, `setMix`, `setHead`, `setHeads`) are public
too, each guarded against re-sending the value it already holds; `setParams` calls them in that order. It is the
whole state: a `TapeDelayParams` built without heads carries the default layout, so a caller that sets heads with
`setHead` and hands over a fresh struct every block hands the default back every block — use one door.

## `DelayLaw.h` — the knob law

`kDivisionNames` / `kDivisionBeats` (13 divisions: the plain values 1/1 … 1/16, each with its dotted value
just before it and its triplet just after), `divisionMs (division, bpm)` (the tempo floored at 40 BPM, so the
slowest value is 6000 ms and fits the 6.1 s line), `timeMsFor (sync, division, bpm, freeMs)` (the sync
policy; whose tempo is the caller's choice), `darkHzForPercent (pct)` (16 kHz → 800 Hz, log-even), and the
knob ranges (`kTimeMinMs` … `kDarkDefaultPct`). All of it in `core::det`, so every row maps a knob to the
same number.

## Contract

- `prepare` allocates (6.1 s + 31 ms of line per channel); `process`/`reset`/setters never allocate, lock
  or throw. A refused `prepare` **disarms**: `process` refuses until a prepare succeeds.
- Law 11: an unprepared stage, a width above the prepared one, or a negative extent is refused before
  anything moves; any `n` is processed in full. A channel at an index `>= nch` has stopped — its line and
  loop filter are dropped once, and its offset head snaps to its target on return. `nch == 0` is a pause:
  the shared motor glides on as through silence.
- `reset()` is a stream restart: lines cleared, write clocks re-anchored at zero, heads re-snapped — what
  follows does not depend on what came before.
- `feed` fades, it does not switch: what enters the record head ramps linearly to nothing, or back, over
  `kFeedRampSeconds` (5 ms), so a flip on a loud note leaves no tick in the repeats. The recirculation is
  never ramped, and a held `feed` is bit-identical to the old switch. A pause moves the ramp as silence would;
  after `prepare`/`reset` it starts settled on the next call's `feed` — a fresh stage does not fade in.
- Law 8: the loop filter's state and the offset glide are flushed below 1e-15 per sample
  (`core::flushDenormal`). Bit-identical under any re-slicing while the state is finite (law 8a: a poison
  restart keeps the call's clock). The heads glide in double and land exactly on their targets.
- Knobs: a non-finite value is ignored, a finite one clamped (TIME 1 … 6100 ms, DARK 1 Hz … 100 kHz,
  OFFSET ±31 ms, REPEATS and MIX 0 … 1) — no knob can reach an index off the line. A head's three fields are
  bounded each on its own: POSITION to [smallest normal float, 1] (anything under one sample of delay plays at
  one sample, the motor's own floor), LEVEL to 0 … 1 (under 1e-15 it is 0), PAN to −1 … +1.
- The heads' controls are sample-clocked like the motor: positions glide per sample and land, the gains ride
  counted ramps that land on target, a pause steps them as silence would, and a silent head stands on its target.
  A mono call ramps each head's side 0 to its level (pan aside) and a stereo call back — no step at the width
  edge. After `prepare`/`reset` the positions snap with the motor, and the gains on the first call that carries a
  channel. Bit-identical under any re-slicing, live head moves included, while the state is finite.
- `tailSeconds()` counts the heads: stacked heads add and hard pan is +3 dB, so the first echo is taken at the sum
  of each head's level times its louder side; the default layout gives the one-head stage's estimate to the bit.
- Threads: prepare, reset, the setters and process belong to one thread; the taps above may be read from any
  thread (each its own atomic — the heads' taps are not one snapshot).
- Math: `core::det` for every transcendental — `exp2`/`log2`/`pow` for the coefficients and the tail, a
  det-built `tanh` per sample in the loop (~5 ns a call) — so the same bits on every row up to FMA
  contraction (law 10: arm64 fuses `a*b + c`, x86-64 and wasm cannot). The pan law is rational,
  gL = (1 − p)/√(1 + p²) and gR = (1 + p)/√(1 + p²), every operation correctly rounded, so it needs no `det`
  call; the gain ramps are pinned against contraction (`core::det::mulAdd`).

## Provenance and parity

The signal path is OrbitAmp's own `DelayStage` (feat/titan, by the same authors), expression for
expression. What changed on the way in is listed below; an out-of-tree harness renders the app's class and
this one side by side and reports the distance. `tests/delay_golden.h` holds regression excerpts of this
module's own renders, so CI guards the sound; `tests/TapeDelayFrozen.h` is the one-head class as it stood before
the heads, and `tests/TapeDelayParityTests.cpp` holds the default layout to it byte for byte.

Measured on arm64 against the app's class (impulse, noise, sweep and their silent tails; 44.1/48/96 kHz;
blocks 1, 64, 512, 333; mono and stereo; six free settings, all 13 divisions, 1/1 at 40 and 1/16T at 240):
**at fixed knobs, at most 3.0e-7 in 1512 renders** — `core::det`'s `tanh` and `exp` against the system
libm's, ulp-level. With the knobs moving the distance is the heads' landing: 2.1e-5 through stand-by and
restarts, 6.3e-5 under one-ulp knob jitter, and 0.55 when TIME travels to 1.5 s — not a level error but
the app's echo arriving up to 28 samples early, where its float motor parked.

**Where it departs from the app, on purpose:** `core::det` instead of the system libm; the law-8 flushes;
the heads glide in double and land (the app's float motor parks up to 1/(2·coeff) ulps short — 112 samples,
2.3 ms, at 6 s); `reset()` re-anchors the write clocks (the app's read fraction after a reset depended on how
long the stream had run: the float position near 2^18 has 1/32-sample steps); knobs are bounded before any
arithmetic (the app turned `setTimeMs(FLT_MAX)` into an infinite target and then a NaN read position);
`prepare` disarms on refusal and clears the loop filter (the app kept the old session's value); the tail
estimate is right from construction (the app reported 0.35 s at its default knobs, the formula says 2.653);
and law 11's stopped-channel and pause rules (the app froze a channel it stopped feeding).
