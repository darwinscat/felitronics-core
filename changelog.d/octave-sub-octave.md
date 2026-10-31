### octave — new module: `Octaver`, the bass amp's sub-octave, and its one-slider law

`felitronics::octave::Octaver` is an analog-style octaver of the OC-2 family, not a pitch shifter. The playing
channels are summed, band-limited to the bass fundamental range (20 … 400 Hz) and tilted toward the lowest partial.
A comparator squares that up, with thresholds at 0.88 of each polarity's own decaying peak, and flip-flops divide
its rising edges by 2 and by 4. Each divided square is scaled by the band's peak over the last cycle and goes
through its own low pass (TONE) and high pass (SUB HP, 25 Hz). One sub is then added to every channel; the dry is
never touched (`latencySamples() == 0`, no dry level).

The note's life is measured in its own periods. An edge a quarter period late fades the sub out over 4 ms, after the
filters, to an exact zero, so a muted note leaves no DC plateau. Two periods without an edge (at most 1/22 s) drop
the thresholds, and the next note needs two edges a plausible period apart, so a lone click makes no sub. A gate
(6 dB hysteresis, 50 ms hold) keeps it an exact zero between notes. Knobs: `sub1`, `sub2`, `toneHz`, `subHpHz`,
`gateDb`; levels glide over 20 ms on a counted ramp, TONE and SUB HP in log frequency.

`OctaveLaw.h`'s `paramsForAmount (amount, menu)` is the face's one slider, a master over the menu's −1 and −2
levels with taper `a/2 + a²/2`: 0 → exact zero, 0.5 → −6.6 dB, 1 → +1.9 dB re the band's peak with the default
menu. −2 is off by default. Tone, sub HP and gate pass through from the menu. With both levels at 0 the stage goes
idle (`isIdle()`): the dry bit for bit, nothing computed. A wake restarts it through a 10 ms input door. A pause
(`nch == 0`) runs it exactly as that much silence and parks on the silent fixed point.

Measured (`OctaverTests`, gcc 14 / MSVC / wasm):
- **Pitch:** the sub is f/2 and f/4 at 41.2 and 30.9 Hz on 44.1, 48 and 96 kHz (asserted within 1e-4; it reads
  0.500000 and 0.250000 to six places).
- **Lock:** half its level 1.05 … 1.17 cycles after a note starts.
- **Decay:** lock holds through a 40 dB, 3-second decay to −46 dBFS.
- **Harmonic-rich notes:** no stray edge on sawtooths, modal plucks and Karplus–Strong from low B to G2.
- **Weak fundamental:** tracks while the fundamental is no more than 11 dB (low B) … 21 dB (G2) under the 2nd
  harmonic in its worst phase.
- **Muted note:** the sub is an exact zero 5/4 of a period + 4 ms after the note's last edge.
- **Retune:** TONE and SUB HP retune adding at most 0.12 of the held sound's own largest curvature.
