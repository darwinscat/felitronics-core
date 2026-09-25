<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# felitronics-core — what's inside (quick map)

JUCE-free, RT-safe, mostly header-only C++20 DSP. Optional compiled backends stay OFF by default,
and the suite is also green under ASan+UBSan. Full design: [`DSP-ARCHITECTURE.md`](DSP-ARCHITECTURE.md).

**Before adding or changing a kernel that carries state**, read
[`LAW8-AUDIT.md`](LAW8-AUDIT.md) — every recursive kernel in the repo has a verdict there, including
the ones that are clean and why, plus the six-question checklist for a new one. It exists because the
same audit had been started three times.

## Foundations (building blocks)

| Module | What | Key types |
|---|---|---|
| `core` | maths, smoothing, delay, RT FFT seam + offline double FFT, denormal flush | `Math`, `Smoother` + `LinearSmoother`, `DelayLine`, `Fft`, `offline::` (`convolve`/`magSpectrum`), `FlushToZero`, `RtStreams` (Audition/Conv/Rec — the RT swap discipline as types) |
| `eq` | filters + the EQ engine + the multiband split | `Svf` (Cytomic), `MatchedBiquad` (Vicanek), `EqEngine`, `Crossover2` (LR4), `MultibandSplitter` |
| `dynamics` | the detector/gain toolkit + compressor (with an external key) + gate + transient | `EnvelopeFollower`, `LinkedDetector`, `GainComputer`, `Compressor`, `NoiseGate`, `TransientShaper`, `ChannelLinker` |
| `oversampling` | polyphase windowed-sinc up/down (alias-free / true-peak); a strict cascade flat to 20 kHz; the `Topology` switch the nonlinear stages offer | `PolyphaseOversampler`, `CascadeOversampler`, `Oversampler` |
| `convolution` | zero-latency partitioned IR convolver + the shared cabinet loader (reference-unity normalization OR, un-normalized, the rate factor a resample costs; resample-on-load, click-free latest-wins swaps) | `PartitionedConvolver`, `MatrixConvolverNupc`, `ConvolutionEngine`, `IrResampler`, `CabConvolver` |
| `lineareq` | linear- & mixed-phase FIR EQ over partitioned convolution | `LinearPhaseEq` (5 quality steps), `NaturalPhaseEq` + `MixedPhaseFir` (φ=k·φ_min "Natural" blend) |
| `neural` | process-only inference seam + swap-safe model holder (backend-free, header-only) | `Inference`, `NeuralStage` |

The guitar-amp modules — `nam` (the NeuralAmpModelerCore backend), `rigplayer` (the `.orbitrig` pack player) and
`poweramp` (the tube power-amp trunk) — live in [felitronics-guitar-core](https://github.com/darwinscat/felitronics-guitar-core).

## Mastering stages

| Module | What | Key types |
|---|---|---|
| `saturation` | oversampled soft-saturation | `WaveShaper` (Tanh/Atan/Cubic/Asym), `Saturator` |
| `stereo` | mid/side image tools | `MidSide`, `MonoBass` (bass mono-maker / elliptical), `StereoWidth` (mono-fold-safe) |
| `dynamiceq` | level-driven EQ band (cut/boost when loud/quiet) | `DynamicEqBand`, `LaneDynamics` (drives an `eq::EqBand` point's per-lane delta seam — what makes the mastering chain's EQ points dynamic) |
| `deesser` | sibilance control, 2 topologies | `DeEsser` (surgical dynamic-EQ · classic split-band) |
| `multiband` | split → per-band processor → recombine (LR4, allpass-flat) | `MultibandProcessor`, `MultibandCompressor`, `MultibandWidth` |
| `dither` | export bit-depth reduction | `Dither` (TPDF + noise shaping; 16/20/24-bit) |
| `limiter` | brick-wall ceiling | `TruePeakLimiter` (oversample → limit → down) |

## Meters · analysis · measurement

| Module | What | Key types |
|---|---|---|
| `analysis` | RT metering, the analyser taps and the spectrum panes (JUCE-free display pipelines); **offline display curves** (`::offline`) | `LoudnessMeter` (LUFS M/S/I + **LRA**), `TruePeakMeter` (dBTP, BS.1770-4 — the spec's short filter, for live display), `ReferenceTruePeakMeter` (dBTP, the 4×/128-tap reference a delivered file is certified and a delivered ceiling aimed with), `CorrelationMeter`, `KWeightingFilter`, `SpectrumTap`, `RollingSpectrumTap` (hop ≠ window, reports the hop that happened), `PlotMap`, `SpectrumPaneT<Fft>` (single FFT, the classic look), `MultiResSpectrumPaneT<…, Fft>` (constant-Q from several FFT lengths — `docs/ANALYZER-MULTIRES.md`), `MultiResSpectrumPaneFastT<…, Fft>` (the same pane, ~1.8–2.3× cheaper on pffft, fill bit-identical — `docs/PERF-ANALYZER-MULTIRES.md`), `offline::logMagnitudeCurve` (1/N-oct, log-f), `offline::interferenceDb` |
| `measurement` | **offline** IR capture: ESS/Farina sweep + deconv, IR post, capture gate (+ standalone sweepless peak/flat-top clip scan), multi-mic align, fine time/polarity align by cross-correlation (message-thread, double) | `Sweep`, `Deconvolve`, `IrPost`, `CaptureGate`, `PeakClip`, `MicSetAlign`, `XcorrAlign`, `ModelGuess` |
| `blend` | **offline** multi-mic IR blend engine: per-mic gain/phase/shift/HPF/LPF + master, solo/mute — the canonical home of the blend defaults | `StripParams`/`MasterParams`, `Filter`, `blendIrs`, `processedMic`, `Overlay` (`makeOverlay` — the one-call mix-view facade) |
| `io` | **offline** file I/O: minimal self-contained WAV read/write (moved from OrbitCapture's `oc/wav.hpp`; zero-dep, loud rejects, memory + file readers) | `WavData`, `readWav`, `readWavMemory`, `writeWav`, `writeWavMonoF32` |

**Build & test:** `cmake -S . -B build -DFELITRONICS_BUILD_TESTS=ON && cmake --build build -j && ctest --test-dir build`.
Add `-DFELITRONICS_WITH_PFFFT=ON` for the optional compiled SIMD FFT backend.

**Every stage of a mastering chain is in core:** saturation → dynamic-EQ → de-esser → multiband comp →
stereo width → transient → mono-bass → dither, metered by true-peak (dBTP) + LUFS/LRA.

The chain BUILT from them — `mastering` (`MasteringChain`, `OfflineRenderer`, the target-loudness solver,
delivery at another rate) — and the offline programme analyzers over the meters above (`ProgrammeReport`,
`SourceForensics`, `HumDetector`, `LowEnd`, `BandBursts`, `ClipDetector`, `WaveformPeaks` / `StereoColumns`, …,
target `felitronics::analysis_offline`) live in
[felitronics-mastering-core](https://github.com/darwinscat/felitronics-mastering-core), with their C ABIs and
the wasm build.
