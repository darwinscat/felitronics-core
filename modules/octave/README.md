<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::octave`

The bass amp's sub-octave — `Octaver` — and the law its one slider drives it with (`OctaveLaw.h`). Header-only,
JUCE-free, depends on `core`, `eq` (`Svf`) and `dynamics` (`EnvelopeFollower`).

```cmake
target_link_libraries(app PRIVATE felitronics::octave)
```

## `Octaver`

The analog octaver of the OC-2 family, not a pitch shifter: it follows the string's zero crossings and divides
them, so nothing waits for a window and the dry is never delayed (`latencySamples() == 0`).

- **Tracker.** The playing channels summed (÷ nch), band-limited to the bass fundamental range (2-pole high pass
  20 Hz, 2-pole low pass 400 Hz), tilted by a one-pole at 50 Hz, squared up by a comparator whose thresholds are
  0.88 of each polarity's own peak (40 ms release). An edge sooner than 1/480 s after the last is ignored.
- **Divider.** Each accepted rising edge flips −1; every second flip flips −2. The step lands mid-sample where the
  crossing lay.
- **Voice.** Each square × the band's peak over the last cycle (up at once, down at an edge, 10 ms glide) × a lock
  gain, through its own 2-pole low pass (TONE), then × its level. At level 1 the sub's peak matches the band's.
- **Gate.** The sub plays while the band is above `gateDb` (6 dB hysteresis, 50 ms hold) and an edge came within
  1/22 s; otherwise it fades out over 4 ms. Between notes it is an exact zero.
- **Monophonic.** A chord's sum has no one fundamental; the comparator follows whatever the sum crosses — usually a
  loud low root, otherwise the beating, which comes out as a growl.

```cpp
felitronics::octave::Octaver o;
o.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 kHz … 768 kHz
o.setParams ({ .sub1 = 1, .sub2 = 0, .dry = 1, .toneHz = 120, .gateDb = -50 });
o.process (io, nch, n);                          // [[nodiscard]] bool, in place
o.reset();                                       // stream restart, idempotent
```

The setters (`setSub1`, `setSub2`, `setDry`, `setToneHz`, `setGateDb`) are public too. The three levels glide over
20 ms on a counted ramp that lands exactly; at both sub levels 0 and dry 1 the stage drops its tracker at the next
64-sample grid point and from then on does not touch the buffer (bit-exact dry, nothing computed).

## The one slider (`OctaveLaw.h`)

`paramsForAmount (amount)`: the −1 octave rises under a unity dry, `sub1 = 1.25 · (a/2 + a²/2)` — 0 → exact zero,
0.25 → −14.2 dB, 0.5 → −6.6 dB, 0.75 → −1.7 dB, 1 → +1.9 dB re the band's peak. The −2 octave stays out of the
slider: on a bass it lands at 7.7 Hz (low B) to 24.5 Hz (G2), where only its square's upper partials are heard
while its fundamental spends headroom and speaker excursion. It is in the menu (`sub2`).

## How the tracker was chosen

A prototype swept the comparator's input filter (2- and 4-pole low passes at 50 … 200 Hz, cascades of one-poles
at 25 … 70 Hz), the threshold fraction (0.3 … 0.94), separate or shared polarity peaks, and the threshold release
(20 … 150 ms) over sines, sawtooths, modal plucks (bridge, middle and neck pickups), noise-excited Karplus–Strong
and weak-fundamental tones at 30.9 … 392 Hz. Two findings decided it: a single symmetric threshold STICKS on a
sawtooth (its filtered wave swings twice as far one way — a 30.9 Hz saw made two edges in three seconds), so each
polarity has its own peak; and thresholds near the peaks (0.88) with a release shorter than the lowest period
reject a harmonic's second bump far better than any filter alone (the best variant with one symmetric threshold
mis-tracked 34 of 168 fixtures, with per-polarity thresholds at 0.3 … 0.7 of the peak 20 of 168; this one 5 of 240, all five on low B: four white-noise-excited Karplus–Strong strings and one
tone whose fundamental sits 8 dB under its 2nd harmonic in an unlucky phase — each followed the 2nd harmonic for the
whole note, so the sub came out at the note's own pitch).

Measured limit (`OctaverTests`, the fundamental and the 2nd harmonic alone, eight phases): it tracks while the
fundamental is no more than 11 dB under the 2nd harmonic at low B, 13 at low E, 17 at A1, 21 at G2. Odd partials
mark the period and only help.

## Listening

`fcore_octave_render <out-dir> [in.wav ...]` renders a bass DI (or, with none, a synthesised bass line) through the
slider at 0.25 / 0.5 / 0.75 / 1, the −2 octave under the dry, and each octave alone, as 24-bit WAVs.
