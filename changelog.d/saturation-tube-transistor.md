### saturation — Tube and Transistor shapes

`WaveShaper::Shape` gains `Tube` and `Transistor`, appended as 4 and 5 so no existing value moves. Neither has a knob of
its own: drive, mix, output and autoComp work as for every shape, and `Params::bias` stays Asym's. **Tube** is
tanh(u+0.3) - tanh(0.3), evaluated as q·t/(1+c·t) so it keeps float precision at low drive (within 1e-6 of a double
evaluation at every drive; the plain difference is off by 5.1e-4 at driveDb 0). It adds mostly the 2nd harmonic — H2
leads H3 by 9.3 dB at driveDb 6 and -3 dBFS — its negative peak is the larger one, and the Saturator's DC blocker runs
for it as for Asym. With autoComp 1 its small-signal slope below 1 (0.936 at the least, near driveDb 3) lifts a
full-scale sine's peak by up to +0.58 dB before the blocker settles and +0.17 dB after. **Transistor** is
u/(1+u⁴)^(1/4): odd, with no cubic term, so it is cleaner than Tanh at low drive (H3 -35.3 dBc against Tanh's -28.7 at
driveDb 6 and -3 dBFS) and harder past its knee; |u| is clamped to 64 so a huge input at a huge drive keeps its sign.
Tanh, Atan, Cubic and Asym are unchanged bit for bit, settled and gliding.
