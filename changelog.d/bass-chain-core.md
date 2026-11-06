### core — `BypassFade`, the crossfade of a link that replaces the signal; `approximatelyEqual` and `roundToInt` with JUCE's semantics

`felitronics::core::BypassFade` moved here from OrbitAmp (`src/core/BypassFade.h`), JUCE-free and its arithmetic
unchanged: a link that replaces the signal (a captured amp, a cabinet, a compressor's blend) is crossfaded against
the signal it was handed over 15 ms of audio, whatever the host's blocks — `advance (n, target)` returns the block's
`Span` (from, to, moving, the ramp counted within the block), `blend` walks dry → wet for the ramp and holds after
it, `runs` keeps a link alive through its fade-out, `snapTo` stands a prepare at an end. Settled, it reports no
motion and costs nothing. An empty block in a fade keeps it a fade.

Two spellings came with it, in `Math.h`, so code moved out of a JUCE plug-in decides as it did there:
`approximatelyEqual (a, b)` is JUCE 8's at its default tolerance (a non-finite pair compares exactly; otherwise the
difference passes at or under the smallest normal value or one epsilon of the larger magnitude), and
`roundToInt (v)` is JUCE's 1.5 · 2^52 trick, ties to even — not `std::lround`, which rounds them away from zero.

Measured (`modules/core/tests/BypassFadeTests.cpp`, `CoreTests.cpp`): a switch settles in 15 ms of audio (±1
sample) at 44.1, 48 and 96 kHz in blocks of 1 … 4096; one 4096-sample block carries the whole 720-sample fade at
48 kHz and holds at the wet after it; the empty block, the turn mid-way and the landing; the two helpers on their
edges (one and two ulps, the smallest normal, a subnormal, infinities and NaN; ties at ±2.5, 3.5, 0.5f, 1.5f).
