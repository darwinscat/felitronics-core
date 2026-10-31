### octave — new module: `Octaver`, the bass amp's sub-octave, and its one-slider law

`felitronics::octave::Octaver` is an analog-style octaver of the OC-2 family, not a pitch shifter: the playing
channels are summed, band-limited to the bass fundamental range (20 … 400 Hz) and tilted toward the lowest partial,
a comparator squares that up with thresholds at 0.88 of each polarity's own decaying peak, and flip-flops divide its
rising edges by 2 and by 4. Each divided square is scaled by the band's peak over the last cycle, rounded by its own
low pass (TONE) and set by its level; one sub goes to every channel, and the dry is never delayed
(`latencySamples() == 0`). A gate (6 dB hysteresis, 50 ms hold) and a lock timeout (no edge for 1/22 s) fade the
sub out over 4 ms, so it is an exact zero between notes. Knobs: `sub1`, `sub2`, `dry`, `toneHz`, `gateDb`; the levels
glide over 20 ms on a counted ramp. `OctaveLaw.h`'s `paramsForAmount (amount)` is the face's one slider: the −1
octave rises under a unity dry (0 → exact zero, 0.5 → −6.6 dB, 1 → +1.9 dB re the band's peak); the −2 octave stays
in the menu. At both sub levels 0 and dry 1 the stage drops its tracker and does not touch the buffer.

Measured (`OctaverTests`): the sub is f/2 and f/4 at 41.2 and 30.9 Hz on 44.1, 48 and 96 kHz (asserted within 1e-4; it reads 0.500000 and
0.250000 to six places), with the
note's own frequency more than 100 dB under the −1 sub's fundamental; it reaches half its level 0.66 … 0.74 cycles
after a note starts; it holds lock through a 40 dB, 3-second decay to −46 dBFS (every half-period one note period;
the sub falls 36.9 dB with the note's 36.7); no stray edge on sawtooths, modal plucks and Karplus–Strong from low B
to G2; it tracks while the fundamental is no more than 11 dB (low B) … 21 dB (G2) under the 2nd harmonic in its
worst phase. The stage sweeps the shared law-11 properties; a pause moves the glides as silence would and drops the
tracker.
