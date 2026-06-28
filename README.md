<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# felitronics-core

A shared, **framework-agnostic, JUCE-free, real-time-safe** C++ DSP core for the Darwin's Cat /
Felitronics product family — EQ, dynamics, convolution, limiting, analysis, neural amp, … — built as
**independent modules** that any product composes: plugins (TabbyEQ, a future compressor/limiter),
the guitar amp plugin, and future **WASM** and **embedded (hardware)** builds.

The core stays pure C++ so it compiles everywhere; host frameworks (JUCE, a WASM shell, firmware)
live only in thin per-platform *adapters*.

➡ **Start here: [`docs/DSP-ARCHITECTURE.md`](docs/DSP-ARCHITECTURE.md)** — the living architecture /
ADR: module layout, portability laws (WASM + embedded are first-class), dependency seams, build/repo
strategy, migration plan, and per-product consumer notes.

> Status: bootstrapping. The first content is the architecture doc; modules land next
> (`dynamics` first, then `eq` is migrated out of TabbyEQ's `teq/`).

License: **AGPL-3.0-or-later**.
