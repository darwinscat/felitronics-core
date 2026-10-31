### chorus — new module: `BassChorus`, the bass amp's chorus, whose lows stay dry and centred

`felitronics::chorus::BassChorus` is the dry plus a voice: a copy of the input's Linkwitz–Riley high band (LR4's
high half, at an adjustable crossover, 100 … 1000 Hz) read by a head that swings on a sine around a centre delay.
The dry is never split — it leaves the stage sample for sample as it came in, at every setting — and the low end is
kept out of the voice by LR4's 24 dB/oct: attenuated, not excluded. On a stereo bus the two sides carry two voices
in quadrature, so a mono sum still swings; a mono call is one voice. The read is a 4-point cubic Hermite in double.

The one-slider law is `ChorusLaw.h`: `paramsForAmount (amount, menu)` drives depth (0 … 3 ms) and mix (0 … 1)
together and passes rate, crossover and centre delay through from the menu; a non-finite amount is ignored.
`pitchDeviationCents (params)` is the read-out a face may print. Amount 0 is the dry bit for bit and computes
nothing: the stage idles on the sample its counted mix ramp lands on zero, and a wake is a restart.

Measured (`modules/chorus/tests/BassChorusTests.cpp`, 44.1 / 48 / 96 kHz): the voice plus LR4's low half of the
same input sums flat to 4e-6 dB from 20 Hz to 20 kHz, and the voice reads -6.02 dB at the crossover and -24.6 dB
an octave below; the Side, (L − R)/2, of a tone under the crossover stays under the voice's own LR4 level there
at every rate, centre and crossover extreme, sidebands included — the open strings' fundamentals read -76 dB (low
B) to -28 dB (G2) at the default 200 Hz crossover; the full 3 ms depth holds on the lowest centre (5 ms); the
voice stops on the sample the mix lands on zero; the pitch deviation measured on a 1 kHz sine is the set one within 0.3 %, the two voices a quarter
period apart; a stream that idles, wakes, moves every knob, and glides the centre into and out of the depth's
limit keeps its largest third difference under the bound of a smooth signal; a lane that returns after a narrow
stretch or a gap fades in, its largest step the steady voice's own; `tailSeconds()` (71 ms) leaves the voice
under -100 dB at the lowest crossover. The voiced band rises about +3 dB at amount 1 (up to +6 on coherent
peaks): budget that before a limiter. Plus the law-11 properties, re-slicing, restart, pause, stopped channel,
poison, a minute of silence with no subnormal in the lines or the glides, and no allocation.
