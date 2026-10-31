<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# `felitronics::octave`

The bass amp's sub-octave — `Octaver` — and the law its one slider drives it with (`OctaveLaw.h`). Header-only,
JUCE-free, depends on `core`, `eq` (`Svf`) and `dynamics` (`EnvelopeFollower`).

```cmake
target_link_libraries(app PRIVATE felitronics::octave)
```

## `Octaver`

The analog octaver of the OC-2 family, not a pitch shifter: it follows the string's zero crossings and divides
them, so nothing waits for a window. The dry is never touched: each output sample is its input sample plus the sub,
the same sub on every channel (`latencySamples() == 0`). There is no dry level — the amp has its own mix.

- **Tracker.** The playing channels summed (÷ nch), band-limited to the bass fundamental range (2-pole high pass
  20 Hz, 2-pole low pass 400 Hz), tilted by a one-pole at 50 Hz, squared up by a comparator whose thresholds are
  0.88 of each polarity's own peak (40 ms release). An edge sooner than 1/480 s after the last is ignored.
- **The note's life.** The period is the longest of the last three edge intervals. An edge a quarter period late:
  the note is dying, and the sub fades out over 4 ms after the filters (back in on the next edge). Two periods without
  one (at most 1/22 s): the note is lost, the thresholds drop, and the next note must show TWO edges a plausible
  period apart before the sub plays — a lone click in a gap makes no sub.
- **Divider.** Each accepted rising edge flips −1; every second flip flips −2. The step lands mid-sample where the
  crossing lay.
- **Voice.** Each square × the band's peak over the last cycle (up at once, down at an edge, 10 ms glide), through
  its own 2-pole low pass (TONE) and 2-pole high pass (SUB HP), then × its level × the lock fade. At level 1 the
  square swings to the band's peak before the filters.
- **Gate.** The sub plays while the band is above `gateDb` (6 dB hysteresis, 50 ms hold) and the note is locked and
  not dying. Between notes it is an exact zero. A gate down at −90 dBFS opens on a DI's own noise floor and the
  tracker follows the hiss; −50 is the default.
- **Monophonic.** A chord's sum has no one fundamental; the comparator follows whatever the sum crosses — usually a
  loud low root, otherwise the beating, which comes out as a growl.

```cpp
felitronics::octave::Octaver o;
o.prepare (sampleRate, maxBlock, maxChannels);   // [[nodiscard]] bool; 1 or 2 channels, 8 kHz … 768 kHz
o.setParams (felitronics::octave::paramsForAmount (amount, menu));
o.process (io, nch, n);                          // [[nodiscard]] bool, in place: io += sub
o.reset();                                       // stream restart, idempotent
o.isLocked();                                    // relaxed atomic, any thread — a tuner light
```

The setters (`setSub1`, `setSub2`, `setToneHz`, `setSubHpHz`, `setGateDb`) are public too. The two levels glide over
20 ms on a counted ramp that lands exactly; TONE and SUB HP glide to a new corner in log frequency over 20 ms,
re-designed every 16 samples of audio time (an instant swap under live state is a click: 1075 times the held filter's
curvature for 120 → 2000 Hz in mid-edge). At both levels 0 the stage drops its tracker at the next 64-sample grid
point and from then on does not touch the buffer (`isIdle()`, bit-exact dry, nothing computed); a wake is a restart
through a 10 ms input door. A pause (`nch == 0`) runs the whole instrument as that much silence, and parks on the
silent fixed point (a ten-minute pause in one call: about 1 ms).

## The one slider (`OctaveLaw.h`)

`paramsForAmount (amount, menu)`: a MASTER over the octaves the menu voiced. It owns `sub1` and `sub2`: each is the
menu's value of that field — the level at the slider's top — times `c(a) = a/2 + a²/2`. Tone, sub HP and gate pass
through. Amount 0 is an exact zero for both whatever the menu says (idle); a non-finite amount comes back as NaN in
the two owned fields, which the stage's setters ignore, so the levels stand. With the default menu (−1 at 1.25,
−2 off): 0.25 → −14.2 dB, 0.5 → −6.6 dB, 0.75 → −1.7 dB, 1 → +1.9 dB re the band's peak.

The −2 octave is off in the default menu: on a bass it lands at 7.7 Hz (low B) to 24.5 Hz (G2), where only its
square's upper partials are heard while its fundamental spends headroom and speaker excursion.

## What it sends down the chain

At the slider's top the sub's square swings to 1.25 × the band's peak before its filters, so the sum can peak well
over the dry (the default renders: −10.0 dBFS in, see `~/CLAUDE_FILES/bass-fx/octave/README.txt` for the measured
peaks). SUB HP (default 25 Hz, 2-pole Butterworth) keeps what no amp model or speaker should be asked for out of it:
the −1 octave is cut by 9.0 dB on the open B (15.4 Hz), 5.0 dB on the open E (20.6 Hz), 2.3 dB on the open A
(27.5 Hz), 0.9 dB on D (36.7 Hz), 0.3 dB on G (49.0 Hz); the −2 by 20.4 dB on B (7.7 Hz) and 15.5 dB on E (10.3 Hz).
The measured cut matches the formula within 0.5 dB (`OctaverTests`). At the default the high pass also turns a held
square into a pulse that decays in about 15 ms. OCTAVE runs before an envelope filter: its sub raises that
filter's detector level, so the filter opens further with OCTAVE up.

## How the tracker was chosen

A prototype swept the comparator's input filter (2- and 4-pole low passes at 50 … 200 Hz, cascades of one-poles
at 25 … 70 Hz), the threshold fraction (0.3 … 0.94), separate or shared polarity peaks, and the threshold release
(20 … 150 ms) over sines, sawtooths, modal plucks (bridge, middle and neck pickups), noise-excited Karplus–Strong
and weak-fundamental tones at 30.9 … 392 Hz. Two findings decided it: a single symmetric threshold STICKS on a
sawtooth (its filtered wave swings twice as far one way — a 30.9 Hz saw made two edges in three seconds), so each
polarity has its own peak; and thresholds near the peaks (0.88) with a release shorter than the lowest period
reject a harmonic's second bump far better than any filter alone (the best variant with one symmetric threshold
mis-tracked 34 of 168 fixtures, with per-polarity thresholds at 0.3 … 0.7 of the peak 20 of 168; this one 5 of 240,
all five on low B: four white-noise-excited Karplus–Strong strings and one tone whose fundamental sits 8 dB under
its 2nd harmonic in an unlucky phase — each followed the 2nd harmonic for the whole note, so the sub came out at the
note's own pitch).

Measured limit (`OctaverTests`, the fundamental and the 2nd harmonic alone, eight phases; pinned at 10 / 12 / 16 /
20): it tracks while the fundamental is no more than 11 dB under the 2nd harmonic at low B, 13 at low E, 17 at A1,
21 at G2. Odd partials mark the period and only help.

Two things were tried for the note's life and are NOT in it. Re-anchoring the thresholds when an edge is late fed
itself on a weak fundamental (one double trigger halved the period, the next true edge looked late, the lowered
thresholds took the 2nd harmonic's bump — the low-B limit fell from −11 to −3 dB). And releasing the level before the
filters handed the sub's high pass a falling step, which came out as one more pulse; the fade after the filters ends
in an exact zero instead.

## Listening

`fcore_octave_render <out-dir> [in.wav ...]` renders a bass DI (or, with none, a synthesised bass line) through the
slider at 0.25 / 0.5 / 0.75 / 1, the −2 octave under the dry, and each octave alone, as 24-bit WAVs.
