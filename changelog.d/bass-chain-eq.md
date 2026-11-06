### eq — `BassCuts` and `CutBiquad`, the bass amp's own cuts; `EqLink`, one EQ link of an amp's chain

Both moved here from OrbitAmp, JUCE-free and their arithmetic unchanged, so the amps' plug-ins and their web port
play the same filters.

`felitronics::eq::CutBiquad` is a 2nd-order Butterworth high or low pass (RBJ, transposed direct form II in double).
`felitronics::eq::BassCuts` is the bass amp's HPF and LPF, 12 dB/oct each, up to two channels: settled and off it is
a wire, bit for bit; on, the biquads plainly; a switch is a 15 ms crossfade (and a cut faded out forgets its memory);
a frequency move glides log-evenly, ~30 ms to a percent, its coefficients redrawn on a 32-sample cadence counted from
the snap, and lands exactly; a NaN that came in is not kept in the recursion.

`felitronics::eq::EqLink` is the console grammar in matched biquads: HPF and LPF with a real slope (6, 12, 18, 24 or
48 dB/oct), two shelves with free corners, two tone bells and a narrow third that switches in. A band at exactly 0 dB
is not run, so a flat link is bit-transparent; `magnitudeDb`, `filterMagnitudeDb` and `toneMagnitudeDb` read the
coefficients `process` runs. Its `Settings` compare member for member, spelt out with `core::exactlyEqual`.

Measured (`modules/eq/tests/BassCutsTests.cpp`, `EqLinkTests.cpp`, the latter moved with the class): -3.01 dB at the
corner, -12.3 dB an octave under a high pass, unity in the pass band; off a wire and on the biquads, sample for
sample; blocks of 96 and 32 the same samples, and blocks of 1 within an ulp (1.0e-7) — the rest of a cadence step
after a fade lands writes `d + (y - d) · 1`, which can round an ulp from `y`: the plug-in's arithmetic, held to what
it is; the 15 ms switch and the glide's exact landing; finite after a NaN. EqLink's gate as it was: the shelves'
plateaus, the bells' centres, the slope ladder two octaves under an HPF, and the curve equal to the sound.
