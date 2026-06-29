<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# felitronics-core

A shared, **framework-agnostic, JUCE-free, real-time-safe** C++20 DSP core for the Darwin's Cat /
Felitronics product family — EQ, dynamics, convolution, limiting, analysis, neural amp, … — built as
**independent modules** that any product composes: plugins (TabbyEQ, a future compressor/limiter),
the guitar amp plugin (OrbitCab), and future **WASM** and **embedded (hardware)** builds.

The core stays pure C++ so it compiles everywhere; host frameworks (JUCE, a WASM shell, firmware)
live only in thin per-platform *adapters*.

➡ **Start here: [`docs/DSP-ARCHITECTURE.md`](docs/DSP-ARCHITECTURE.md)** — the living architecture /
ADR. Companions: [`docs/dsp-architecture-third-review.md`](docs/dsp-architecture-third-review.md) (panel
review) and [`docs/dsp-shared-dsp-review.md`](docs/dsp-shared-dsp-review.md) (cross-plugin shared-DSP +
de-JUCE plan).

## Modules

| Module | Namespace | Target | What | Status |
|---|---|---|---|---|
| `core` | `felitronics::core` | `felitronics::core` | `Sample` alias, size config (SSOT), `Math`, `Smoother`, software denormal flush + `ScopedFlushToZero`, `DelayLine`, the **`Fft` seam** (+ scalar backend) | header-only |
| `analysis` | `felitronics::analysis` | `felitronics::analysis` | `SpectrumTap` (SPSC) + BS.1770 `LoudnessMeter` (LUFS, ffmpeg-validated) + `CorrelationMeter` | header-only |
| `dynamics` | `felitronics::dynamics` | `felitronics::dynamics` | `EnvelopeFollower` + `GainComputer` + `ChannelLinker` + `GainReductionFollower` + a broadband `Compressor` | header-only |
| `eq` | `felitronics::eq` | `felitronics::eq` (+ `teq::core` compat) | matched biquads (Vicanek) + Cytomic SVF + `EqBand` + `EqEngine` — migrated from TabbyEQ's `teq/` | header-only |
| `convolution` | `felitronics::convolution` | `felitronics::convolution` | zero-latency partitioned FFT convolver + click-free `ConvolutionEngine` (dual-engine IR swap) + offline Kaiser IR resampler | header-only |
| `oversampling` | `felitronics::oversampling` | `felitronics::oversampling` | polyphase windowed-sinc FIR up/down sampler (true-peak / alias-free) | header-only |
| `limiter` | `felitronics::limiter` | `felitronics::limiter` | brickwall **true-peak** limiter (guaranteed ceiling; oversample→limit→downsample) | header-only |
| `neural` | `felitronics::neural` | `felitronics::neural` | process-only `Inference` seam + swap-safe `NeuralStage` (the NAM/Eigen backend + model loading live in the adapter) | header-only |

Every module is JUCE-free, declares its target tiers, and ships its own JUCE-free self-tests
(measured audio == analytic curve — the `teq` discipline).

## Build & test

JUCE-free, so configuring does **not** fetch JUCE:

```sh
cmake -S . -B build -DFELITRONICS_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Consume (pinned by tag, like JUCE)

```cmake
include(FetchContent)
FetchContent_Declare(felitronics_core
    GIT_REPOSITORY <repo-url>
    GIT_TAG        <pinned-tag>)
FetchContent_MakeAvailable(felitronics_core)

target_link_libraries(app PRIVATE felitronics::eq felitronics::dynamics)   # link only what you use
# #include <felitronics/eq/EqEngine.h>     (new code)
# #include <teq/EqEngine.h>                (transitional compat — teq:: aliases felitronics::eq)
```

> Status: **bootstrapping.** Light modules (`core`, `analysis`, `dynamics`, `eq`) land first; `eq` is
> migrated out of TabbyEQ's `teq/` with a `teq::` compat shim so consumers repoint one at a time. The
> FFT-seam spike has **landed** (`core/Fft` + `convolution`, the architecture's keystone); the remaining
> heavy modules (`neural`, `limiter`, `oversampling`) follow.

License: **AGPL-3.0-or-later** (SPDX header on every file).
