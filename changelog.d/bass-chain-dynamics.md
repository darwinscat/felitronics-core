### dynamics — `CharacterComp`, the bass amp's own compressor characters (OPTO, FET, SQUASH)

`felitronics::dynamics::CharacterComp` moved here from OrbitAmp (`src/core/CompCharacters.h`), JUCE-free and its
arithmetic unchanged. One mono compressor whose character (`Voice`) sets the curve — feed-forward, or fed back, its
knee R times wider in the input's terms — the detector (the RMS over 10 ms, or the peaks), a release in one stage or
two (a memory that sinks in while a squeeze is held and lets go slowly), and the gain element's colour growing with
the squeeze: GRAIN (odd harmonics worked at twice the rate, added as a residual with no fundamental of its own) and
DARK (the top let down toward a one-pole low pass). `staticDeltaDb` draws the curve; `seed` takes a squeeze over
from another detector without a step. The detector's key is the caller's: no filter inside.

Measured (`modules/dynamics/tests/CharacterCompTests.cpp`): a steady note settles where `staticDeltaDb` says, fed
forward and fed back (within 0.5 dB); fed back, 2 dB under the threshold, the curve is the quadratic's root
(-0.4503 dB at 4:1, an 8 dB knee) where feed-forward takes -0.1875; a level under the knee is untouched; 300 ms after
it is let go a long squeeze still holds more than twice a short hit's; GRAIN keeps the note's level within 0.15 dB;
DARK takes an 8 kHz tone more than 3 dB lower; one call and calls of 37 samples give the same samples and the same
squeeze; `seed` and `reset`.
