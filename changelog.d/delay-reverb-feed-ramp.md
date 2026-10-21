### delay · reverb — `feed` fades over 5 ms instead of switching, so a stand-by flip no longer ticks in the repeats or clicks into the room

`TapeDelay::process (…, feed)` and `AmpReverb::process (…, feed)` switched what enters the line or the room at
the call boundary. A flip on a loud note was a full-scale step at the record head, which the delay then played
back as a tick at every repeat, and which the reverb carried through its tail as a click. The input now goes
through a linear ramp of `kFeedRampSeconds` (5 ms, counted in samples at the prepared rate: 240 at 48 kHz) from
1 to 0 or from 0 to 1. A flip in mid-ramp reverses from where the ramp stands. The recirculation and the room's
own tail are not ramped.

Measured on a 0 dBFS 375 Hz tone flipped on a crest at a block boundary: the delay's largest sample-to-sample
step where the first and second repeats replay the flip falls from 14.8x, 13.3x and 15.6x the steady echo's own
step (off, on, and off-then-on within 128 samples) to 0.98x, 1.05x and 1.00x. In a PLATE, the room's largest
sample above 4 kHz in the 250 ms after the flip falls from -24.7, -24.7 and -23.5 dB re its steady answer to
-70.1, -70.2 and -63.5 dB.

With `feed` held, the output is bit-identical to before: at the ends of the ramp the input is the plain sample
or the plain zero, with no multiply. This was checked on 15.8 M output samples (the golden scenarios, 84 random
renders from 1 kHz to 192 kHz with held `feed`, and 24 stand-by renders), and no golden moved. After
`prepare()` or `reset()` the ramp starts settled on the next call's `feed`, so a fresh insert does not fade in
its first note. A pause (`nch == 0`) moves the ramp as the same length of silence would (law 11c). The ramp
lands on exactly 0 or 1 because it is counted, not glided.
