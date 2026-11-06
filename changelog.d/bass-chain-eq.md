### eq — `BassCuts` and `CutBiquad`, the bass amp's own cuts

Moved here from OrbitAmp, JUCE-free and their arithmetic unchanged, so the bass amps' plug-ins and their web port
play the same filters.

`felitronics::eq::CutBiquad` is a 2nd-order Butterworth high or low pass (RBJ, transposed direct form II in double).
`felitronics::eq::BassCuts` is the bass amp's HPF and LPF, 12 dB/oct each, up to two channels: settled and off it is
a wire, bit for bit; on, the biquads plainly; a switch is a 15 ms crossfade (and a cut faded out forgets its memory);
a frequency move glides log-evenly, ~30 ms to a percent, its coefficients redrawn on a 32-sample cadence counted from
the snap, and lands exactly; a NaN that came in is not kept in the recursion.

Measured (`modules/eq/tests/BassCutsTests.cpp`): -3.01 dB at the corner, -12.3 dB an octave under a high pass, unity
in the pass band; off a wire and on the biquads, sample for sample; blocks of 96 and 32 the same samples, and blocks
of 1 within an ulp (1.0e-7) — the rest of a cadence step after a fade lands writes `d + (y - d) · 1`, which can round
an ulp from `y`: the plug-in's arithmetic, held to what it is; the 15 ms switch and the glide's exact landing; finite
after a NaN.
