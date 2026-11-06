### measurement — `bassReferenceUnityGain`, the bass's unity for a cabinet IR

`felitronics::measurement::bassReferenceUnityGain (taps, n, rate)` moved here from OrbitAmp's cabinet
(`src/core/CabinetIr.h`, `bassUnityGain`), its arithmetic unchanged: the gain that passes the bass's reference — white
noise through second-order Butterworth cuts at 40 and 120 Hz, the band a bass line's fundamentals fill — at unity
RMS. `referenceUnityGain`'s reference weighs the top a bass cabinet rolls off and left every bass cabinet 4 to 8 dB
hot on a real bass line; against this one the bass's boxes land within a dB of unity. The family's arithmetic
otherwise: the first second, the ±30 dB clamp, a near-silent IR (under -60 dB) left alone. Mono.

Measured (`modules/measurement/tests/BassUnityTests.cpp`): an impulse at unity and half of one at twice the gain; an
800 Hz one-pole low pass within a tenth of a dB of unity on the bass's reference while the guitar's lifts it more;
the clamp at +30 dB, the floor, and degenerate input.
