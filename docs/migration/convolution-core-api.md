<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# Convolution — core API ready + the verified JUCE spec (OrbitCab cab-path de-JUCE)

The core convolution path is built + tested (JUCE-free, all green). This doc is the **authoritative spec**
for replacing OrbitCab's `juce::dsp::Convolution` without changing the sound. The JUCE behaviours below were
**read out of the vendored source** (`build/_deps/juce-src/.../juce_Convolution.cpp`) and must be replicated
byte-for-byte by the adapter — they were verified line-by-line against that source.

## What core provides (`felitronics::convolution` / `felitronics::core::fft`, all green)

- **`core::fft` seam** — COMPILE-TIME real-FFT backend (concept `RealFftBackend`). Default = scalar reference.
- **`PartitionedConvolver<Fft>`** — zero-latency (direct head + uniform-partitioned overlap-save tail),
  verified == direct convolution under hostile variable-block splits + no-alloc.
- **`ConvolutionEngine<Fft, MaxChannels = 2>`** — PRODUCTION wrapper, **multi-channel with ONE lockstep
  crossfade**. `prepare(partitionSize, maxIrSamples, crossfadeSamples, numChannels)`; `setIr(ir,len)`
  (mono → broadcast to all channels) or `setIr(irPerCh,nch,len)` (true-stereo); `process(in[],out[],nch,n)`
  (planar) or `process(in,out,n)` (mono). The crossfade position is shared across channels and advanced
  once per sample → a stereo IR swap can never move the image (proven: L==R bit-identical through a swap).
  `isBusy()` for coalescing; `latencySamples()==0`.
- **`resampleIr(in,inLen,inSr,outSr,cfg)`** → `std::vector<float>` — offline Kaiser windowed-sinc (~80 dB),
  **DC gain normalized to 1**. Message-thread only (allocates, double math).

## JUCE behaviour the adapter MUST replicate (verified vs `juce_Convolution.cpp`)

| Behaviour | JUCE source | What the adapter does |
|---|---|---|
| **Normalise::no still applies a gain** | `makeEngine`, `applyGain(irSr/hostSr)` (`:792`) | after any resample, multiply IR by `irSr/hostSr` |
| **Normalise::yes = energy norm, −18 dB** | `calculateNormalisationFactor` (`:623`) | `g = (E<1e-8) ? 1 : 0.125/sqrt(E)`, `E = max_ch Σ ir[n]²`, then multiply |
| **IR swap crossfades over 50 ms** | `CrossoverMixer::prepare` `smoother.reset(sr, 0.05)` (`:951`) | `crossfadeSamples = round(0.05 * hostSr)` |
| **default convolution is zero-latency** | `Latency{0}` (`:1206`), mc-latency `0` (`:458`); OrbitCab reports only NAM latency (`PluginProcessor.cpp:810`) | keep `latencySamples()==0` → PDC + wet/dry alignment unchanged |
| **resample only when rates differ** | `resampleImpulseResponse` ratio 1 ≈ identity | **guard: call `resampleIr` ONLY if `irSr != hostSr`** (our Kaiser at ratio 1 would still LPF at fc≈0.475 and color a host-rate IR) |

⚠️ **Resampler is NOT identical.** JUCE uses `ResamplingAudioSource` (linear interp + 2nd-order LPF);
ours is Kaiser windowed-sinc (higher quality). For a **non-host-rate** IR the two won't null-test.
**Decision (desktop quality #1): accept ours** — it's the better resampler and cab IRs are mostly 48 k;
the 44.1 k→48 k difference is tiny and in our favour. (To go bug-for-bug instead, add a JUCE-matching
linear+LPF mode.) Record this as a deliberate, documented improvement.

## Adapter load path (replaces `cab::Convolver::loadIR`)

```
ir = decoded planar IR @ irSr
if (irSr != hostSr) ir[c] = resampleIr(ir[c], irSr, hostSr)        // skip at host rate (see guard above)
if (normalise == yes) g = (E<1e-8 ? 1 : 0.125/sqrt(E)),  E = max_c Σ ir[c][n]²   // byte-fallback path
else                  g = irSr / hostSr                                            // the normal cab path
multiply every ir[c] by g
engine.setIr(ir.data(), len)            // mono  → broadcast (juce Stereo::yes)
engine.setIr(irPtrs, nch, len)          // stereo IR → per channel
```
Off-thread build + `isBusy()` coalescing (mirror the existing reload-coalescing). `prepare(P, maxIr,
round(0.05*hostSr), nch)`. Sensible sizes: `partitionSize = 128`, `maxIrSamples = ceil(maxIrSeconds*hostSr)`
(OrbitCab caps decode at 20 s; a shorter cap is a separate product decision, not this swap).

## The adapter's FFT backend (desktop speed)

Wrap `juce::dsp::FFT` as a `RealFftBackend` (`spectrumFloats(n)=n+2`, `forward`/`inverse` via
`performRealOnly*Transform` + `/N`, `spectralMultiplyAdd` = complex MAC on JUCE's packed layout), then use
`ConvolutionEngine<JuceFftBackend, 2>`. Verify it against the scalar reference with the SAME parity tests
(round-trip, `spectralMultiplyAdd == circular conv`). pffft (BSD) is the JUCE-free desktop/wasm alternative.

## Golden test (prove sound-preserving)

Run `juce::dsp::Convolution` and the new engine side-by-side, **aligned to the first changed sample** (JUCE's
load is async). Tiers: (1) **steady-state** — impulse / −18 dBFS noise / real DI through mono + true-stereo
IRs at host rate, block-size sweep {1,17,64,128,maxBlock,random}, tolerance max-abs ≤ `2e-5`, null RMS ≤
−100 dBFS; (2) **preprocessing** — non-host-rate IR + `Normalise::yes` (expect the resampler delta — assert
the *gain/energy* matches, not bit-equality); (3) **swap** — IR A→B at a known boundary, 50 ms fade, output
matches steady-B after `crossfade+irLen`.

After this swap the whole cab path is JUCE-free → OrbitCab matches TabbyEQ's posture (only the adapter
touches JUCE). The per-slot HPF/LPF is the separate SVF swap (`felitronics::eq::Svf`); the parameter
smoothers are the separate `felitronics::core::LinearSmoother` swap.
