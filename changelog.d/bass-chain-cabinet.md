### cabinet — a new module: `CabinetIr`, an amp's cabinet stage

`felitronics::cabinet::CabinetIr` moved here from OrbitAmp (`src/core/CabinetIr.h`) — its DSP half, JUCE-free and
its arithmetic unchanged; header-only, depends on `core`, `eq`, `measurement` and `convolution`. One IR shot comes in
as taps at its own rate (`load (taps, n, rate)`, false for a shot `takes` refuses: no length, a rate outside 8 kHz …
768 kHz, or longer than 5 s) — decoding a file stays the caller's. What a player does to it — the tail trimmed on a
half-cosine, the 6 … 48 dB/oct cuts in matched biquads, the sign — is baked into the taps (`setPost`), the level taken
after it at the shot's own rate against the guitar's reference or the bass's (`setReference`,
`measurement::bassReferenceUnityGain`), and the taps resampled to the session's rate once per rate before
`convolution::CabConvolver` (zero latency, the 50 ms click-free swap) is handed them; silence under -80 dBFS at
either end is cut first. On the audio thread `process` keeps a non-finite sample out of the convolution and silences
a poisoned output for the impulse's length; `idle` drains the whole schedule once the cabinet stands out of the path,
so it comes back with no ghost, and runs a cabinet handed over meanwhile in on silence. `prepare`, `load`,
`setPost`, `setReference` and `flushPending` are the producer's side and take one lock; `process` and `idle` take none.

Measured (`modules/cabinet/tests/CabinetIrTests.cpp`): after the crossfade the output is the direct convolution
times the reference's gain, guitar and bass, at sample zero (within 2e-4); a quarter of the tail kept reports a
quarter of the length, and the sign flipped is the same cabinet negated; a 48 kHz shot in a 44.1 or 96 kHz session
keeps its level within 0.1 dB; after a NaN every sample is a number and the cabinet still sounds; drained, silence
in is silence out; a refused shot leaves the cabinet it had; no allocation in `process` (stereo and mono) or `idle`.
