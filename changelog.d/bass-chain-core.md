### core — `approximatelyEqual` and `roundToInt` with JUCE's semantics

Two spellings in `Math.h`, so code moved out of a JUCE plug-in decides as it did there: `approximatelyEqual (a, b)`
is JUCE 8's at its default tolerance (a non-finite pair compares exactly; otherwise the difference passes at or
under the smallest normal value or one epsilon of the larger magnitude), and `roundToInt (v)` is JUCE's 1.5 · 2^52
trick, ties to even — not `std::lround`, which rounds them away from zero.

Measured (`modules/core/tests/CoreTests.cpp`): one and two ulps, the smallest normal, a subnormal, infinities and
NaN; ties at ±2.5, 3.5, 0.5f, 1.5f; the 32-bit range's ends.
