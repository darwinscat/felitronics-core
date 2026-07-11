<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# felitronics-core — what's inside (quick map)

JUCE-free, RT-safe, header-only C++20 DSP. **17 modules · 38 test suites · 5612 checks green**
(also green under ASan+UBSan). Full design: [`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md).

## Foundations (building blocks)

| Module | What | Key types |
|---|---|---|
| `core` | maths, smoothing, delay, RT FFT seam + offline double FFT, denormal flush | `Math`, `Smoother` + `LinearSmoother`, `DelayLine`, `Fft`, `offline::` (`convolve`/`magSpectrum`), `FlushToZero` |
| `eq` | filters + the EQ engine + the multiband split | `Svf` (Cytomic), `MatchedBiquad` (Vicanek), `EqEngine`, `Crossover2` (LR4), `MultibandSplitter` |
| `dynamics` | the detector/gain toolkit + compressor + gate + transient | `EnvelopeFollower`, `GainComputer`, `Compressor`, `NoiseGate`, `TransientShaper`, `ChannelLinker` |
| `oversampling` | polyphase windowed-sinc up/down (alias-free / true-peak) | `PolyphaseOversampler` |
| `convolution` | zero-latency partitioned IR convolver | `PartitionedConvolver`, `ConvolutionEngine` (click-free IR swap), `IrResampler` |
| `lineareq` | linear- & mixed-phase FIR EQ over partitioned convolution | `LinearPhaseEq` (5 quality steps), `NaturalPhaseEq` + `MixedPhaseFir` (φ=k·φ_min "Natural" blend) |
| `neural` | process-only inference seam (NAM lives in the adapter) | `Inference`, `NeuralStage` |

## Mastering chain

| Module | What | Key types |
|---|---|---|
| `saturation` | oversampled soft-saturation | `WaveShaper` (Tanh/Atan/Cubic/Asym), `Saturator` |
| `stereo` | mid/side image tools | `MidSide`, `MonoBass` (bass mono-maker / elliptical), `StereoWidth` (mono-fold-safe) |
| `dynamiceq` | level-driven EQ band (cut/boost when loud/quiet) | `DynamicEqBand` |
| `deesser` | sibilance control, 2 topologies | `DeEsser` (surgical dynamic-EQ · classic split-band) |
| `multiband` | split → per-band processor → recombine (LR4, allpass-flat) | `MultibandProcessor`, `MultibandCompressor`, `MultibandWidth` |
| `dither` | export bit-depth reduction | `Dither` (TPDF + noise shaping; 16/20/24-bit) |
| `limiter` | brick-wall ceiling | `TruePeakLimiter` (oversample → limit → down) |

## Meters · analysis · measurement

| Module | What | Key types |
|---|---|---|
| `analysis` | RT metering + the FFT tap; **offline display curves** (`::offline`) | `LoudnessMeter` (LUFS M/S/I + **LRA**), `TruePeakMeter` (dBTP, BS.1770-4), `CorrelationMeter`, `KWeightingFilter`, `SpectrumTap`, `offline::logMagnitudeCurve` (1/N-oct, log-f), `offline::interferenceDb` |
| `measurement` | **offline** IR capture: ESS/Farina sweep + deconv, IR post, capture gate, multi-mic align, fine time/polarity align by cross-correlation (message-thread, double) | `Sweep`, `Deconvolve`, `IrPost`, `CaptureGate`, `MicSetAlign`, `XcorrAlign`, `ModelGuess` |
| `blend` | **offline** multi-mic IR blend engine: per-mic gain/phase/shift/HPF/LPF + master, solo/mute — the canonical home of the blend defaults | `StripParams`/`MasterParams`, `Filter`, `blendIrs`, `processedMic`, `Overlay` (`makeOverlay` — the one-call mix-view facade) |

**Build & test:** `cmake -S . -B build -DFELITRONICS_BUILD_TESTS=ON && cmake --build build -j && ctest --test-dir build`

**A full mastering chain is now buildable in core:** saturation → dynamic-EQ → de-esser → multiband comp →
stereo width → transient → mono-bass → dither, metered by true-peak (dBTP) + LUFS/LRA.
