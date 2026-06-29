<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. -->

# Convolution — core API is ready (for OrbitCab T7 de-JUCE)

The core convolution path is built + tested (JUCE-free). OrbitCab's T7 (de-JUCE the cab path) now has a
concrete target. The only adapter-side piece is a fast FFT backend (JUCE has one); everything else is core.

## What core provides (all green, `felitronics::convolution` / `felitronics::core::fft`)

- **`core::fft` seam** — a COMPILE-TIME real-FFT backend (concept `RealFftBackend`): `spectrumFloats(n)`,
  `prepare(n)`, `forward(real→spec)`, `inverse(spec→real)` (1/N), `spectralMultiplyAdd(a,b,acc)`. The
  spectrum layout is the backend's own (no forced repack). Default = the scalar reference (tests/wasm).
- **`PartitionedConvolver<Fft>`** — zero-latency (direct head + uniform-partitioned overlap-save tail),
  fixed partition + arbitrary host block, RT-safe, verified == direct convolution under hostile variable
  block splits + no-alloc.
- **`ConvolutionEngine<Fft>`** — the PRODUCTION wrapper: `setIr()` (message thread) builds the new IR
  into an inactive engine; `process()` (audio) CROSSFADES old→new (click-free), race-free, zero latency,
  `isBusy()` for coalescing. This replaces `juce::dsp::Convolution`'s internal swap.
- **`resampleIr()`** — offline Kaiser windowed-sinc IR resampler (≥ ~70 dB), for loading an IR at a
  different SR than the host (do NOT use Catmull-Rom for IRs).

## The adapter's only job: a JUCE FFT backend (desktop speed)

Write, in OrbitCab's JUCE adapter, a `RealFftBackend` wrapping `juce::dsp::FFT` so desktop gets SIMD speed:

```cpp
struct JuceFftBackend {
    static int  spectrumFloats (int n) noexcept { return n + 2; }   // JUCE real layout: (n/2+1) interleaved complex
    bool prepare (int n) noexcept;                                   // build juce::dsp::FFT(log2(n)) + scratch
    void forward (const float* real, float* spec) noexcept;         // performRealOnlyForwardTransform → JUCE's packed layout
    void inverse (const float* spec, float* real) noexcept;         // performRealOnlyInverseTransform, then /N
    void spectralMultiplyAdd (const float* a, const float* b, float* acc) noexcept;  // complex MAC on JUCE's layout
};
static_assert (felitronics::core::fft::RealFftBackend<JuceFftBackend>);
```
Then `felitronics::convolution::ConvolutionEngine<JuceFftBackend>`. Verify the backend against the core's
scalar reference with the SAME parity tests (round-trip, `spectralMultiplyAdd == circular conv`) — those
are the contract. (pffft is the alternative for a JUCE-free desktop/wasm build, BSD.)

## Wiring (replaces `cab::Convolver` + the JUCE bits of `cab::IRSlot`)

- `cab::Convolver` (`src/core/Convolver.h`, wraps `juce::dsp::Convolution`) → `ConvolutionEngine<JuceFftBackend>`.
  `loadIR(...)` → `resampleIr()` (if SR differs) then `engine.setIr(...)` (off-thread build + click-free swap).
- `cab::IRSlot::processWet` → `engine.process(...)` (the per-slot HPF/LPF SVF is the separate T3b swap to
  `felitronics::eq::Svf`); drop `juce::dsp::AudioBlock`/`ProcessContextReplacing`/`juce::AudioBuffer`.
- Latency stays ~0 (the convolver is zero-latency). Golden-test the cab output vs the `juce::dsp::Convolution`
  version (within tolerance — it's a versioned-behaviour swap).

After T7 the whole cab path is JUCE-free → OrbitCab matches TabbyEQ's posture (only the adapter touches JUCE).
