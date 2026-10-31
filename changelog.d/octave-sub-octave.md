### octave — new module: `Octaver`, the bass amp's sub-octave, and its one-slider law

`felitronics::octave::Octaver` is an analog-style octaver of the OC-2 family, not a pitch shifter.
- **Tracker.** The playing channels are summed, band-limited to the bass fundamental range (20 … 400 Hz) and tilted
  toward the lowest partial. A comparator squares that up, with thresholds at 0.88 of each polarity's own decaying
  peak and a floor 20 dB under the band's own level.
- **Divider.** Flip-flops divide its rising edges by 2 and by 4.
- **Voices.** Each divided square is scaled by the band's peak over the last cycle and goes through its own low pass
  (TONE) and high pass (SUB HP, 25 Hz). One sub is added to every channel; the dry is never touched
  (`latencySamples() == 0`, no dry level).

The note's life is measured in its own periods.
- Edges count only with the gate open.
- A lost note comes back after two rises and two falls, each fresh twice over: the band's peak over the last period
  stands at the gate's opening level, and the raw input itself has just swung that way. So hum under the gate proves
  nothing, a note whose follower ripples across the gate still acquires, and no one-cycle click — any bass pitch,
  phase, polarity or ramp — makes a sub.
- Deadlines come only from a period confirmed by two agreeing intervals. An edge a quarter period late fades the sub
  out over 4 ms, after the filters, to an exact zero; two periods without one lose the note.
- On reopening, the voices start primed at the new note's level.
- A gate (6 dB hysteresis, 50 ms hold) keeps the sub an exact zero between notes.

Knobs: `sub1`, `sub2` (the levels now), `sub1Top`, `sub2Top` (the slider's top, for the law), `toneHz`, `subHpHz`,
`gateDb`. Levels glide over 20 ms on a counted ramp; TONE and SUB HP glide in log frequency, re-designed every sample.

`OctaveLaw.h`'s `paramsForAmount (amount, menu)` is the face's one slider, a master over the menu's tops with taper
`a/2 + a²/2`: 0 → exact zero, 0.5 → −6.6 dB, 1 → +1.9 dB re the band's peak with the default menu. −2 is off by
default. It never writes the tops, so the current params can be fed back as the menu, and a non-finite amount returns
the menu untouched. With both levels at 0 the stage goes idle (`isIdle()`): the dry bit for bit, nothing computed. A
wake restarts it through a 10 ms input door. A pause (`nch == 0`) runs it exactly as that much silence and parks on
the silent fixed point.

Measured (`OctaverTests`, gcc 14 / MSVC / wasm):
- **Pitch:** the sub is f/2 and f/4 at 41.2 and 30.9 Hz on 44.1, 48 and 96 kHz (asserted within 1e-4; it reads
  0.500000 and 0.250000 to six places).
- **Onset, from silence, worst over 16 phases × both polarities:** the fade opens within 1.51 cycles on low B,
  1.61 on low E and 1.89 on G2. Half level comes within 1.59 / 1.71 / 2.14 cycles (51 / 42 / 22 ms).
- **Relock:** after a louder note, the reopened sub stays within 1.45× the new note's own steady level.
- **Decay:** lock holds through a 40 dB, 3-second decay to −46 dBFS.
- **Harmonic-rich notes:** no stray edge on sawtooths, modal plucks and Karplus–Strong from low B to G2.
- **Weak fundamental:** tracks while the fundamental is no more than 11 dB (low B) … 21 dB (G2) under the 2nd
  harmonic in its worst phase.
- **Muted note:** the sub is an exact zero 5/4 of a period + 4 ms after the note's last edge.
- **Hum:** with 60 Hz hum at −80 / −65 dBFS, low B and E lock on their own edges and the hum alone passes bit-exact;
  a loud G2 stopping dead over −58 dBFS hum leaves no sub 100 ms on.
- **Near the gate:** a sustained low B whose follower ripples 0.84 … 1.28× the opening level is acquired within 90 ms.
- **Retune:** TONE and SUB HP move at most 0.5 % per sample over their whole range.
