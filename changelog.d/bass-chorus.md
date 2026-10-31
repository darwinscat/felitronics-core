### chorus — new module: `BassChorus`, the bass amp's chorus, whose lows stay dry, mono and in phase

`felitronics::chorus::BassChorus` is the dry plus a voice: a copy of the input's Linkwitz–Riley high band (LR4's
high half, at an adjustable crossover, 40 … 1000 Hz) read by a head that swings on a sine around a centre delay.
Nothing below the crossover is delayed, modulated or spread, and the dry is never split — it leaves the stage
sample for sample as it came in, at every setting. On a stereo bus the two sides carry two voices in opposite LFO
phase; a mono call is one voice. The read is a 4-point cubic Hermite.

The one-slider law (`chorus::law::paramsFor`, `BassChorus::setAmount`): `amount` 0..1 drives depth (0 … 3 ms)
and mix (0 … 1) together; rate and crossover are the menu's. Amount 0 is the dry bit for bit and computes
nothing: the stage idles once its mix has glided onto an exact zero.

Measured (`modules/chorus/tests/BassChorusTests.cpp`, 44.1 / 48 / 96 kHz): the voice plus LR4's low half of the
same input sums flat to 4e-6 dB from 20 Hz to 20 kHz (the bar is 0.1 dB), and the voice reads -6.02 dB at the
crossover and -24.6 dB an octave below; at amount 1, two octaves under a 200 Hz crossover, L − R is -61 dB under
the dry and each side is the dry within 0.03 dB and 0.2°, three octaves under -91 dB; the pitch deviation
measured on a 1 kHz sine is the set `2π · rate · depth` within 0.3 %, the two voices in opposite phase
(correlation ≤ -0.997); a stream that idles, wakes and moves every knob keeps its largest third difference at
half the bound of a smooth signal. The head and the mix glide critically damped: with one-pole glides the
slider's drop to zero was a 6 % pitch jump of the voice, 1.1x to 3.7x past that bound. Plus the law-11
properties, re-slicing, restart, pause, stopped channel, poison, a minute of silence with no subnormal, and no
allocation.
