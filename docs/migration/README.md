<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->
<!-- Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. -->

# Migration playbook — repointing the products onto felitronics-core

The generic, public map. **The executable, per-plugin agent prompts live in each plugin's own
`docs/migration/`** — that keeps any product-private detail out of this public core repo, and puts each
task where the agent actually works.

## The compat contract (why a repoint is ~zero source change)

`felitronics-core` ships, alongside the canonical `felitronics::eq`, a transitional compat layer:
- a **real CMake target `teq_core`** (+ `teq::core` alias) that carries `felitronics::eq`;
- `<teq/*.h>` **re-export shims** that alias the moved names back into `namespace teq`
  (`teq::EqEngine`, `teq::BandParams`, `teq::Svf`, `teq::Smoother`, `teq::SpectrumTap`, `teq::kMaxChannels`, …).

So a consumer that says `#include <teq/EqEngine.h>` / `teq::EqEngine` keeps compiling after it repoints
its `FetchContent` from the old `teq` source to `felitronics-core` — **only the CMake fetch block
changes.** felitronics-core is JUCE-free, so this fetch never pulls JUCE.

## Sequence (each step keeps the plugin releasable + verified)

1. **T1 — TabbyEQ repoint.** Stop owning `teq/` locally; `FetchContent` felitronics-core, link
   `teq::core`. Build VST3/AU/Standalone/CLAP + run the (now-upstream) tests + `auval`.
   → `tabby-eq/docs/migration/`.
2. **T2 — OrbitCab repoint.** Swap the `teq` fetch (was `darwinscat/tabby-eq`, `SOURCE_SUBDIR teq`) for
   felitronics-core; keep `teq::core`. Build VST3/AU/CLAP + cab tests + `auval`.
   → `orbitcab/docs/migration/`.
3. **T3 — OrbitCab dedups** (each behind a golden magnitude/phase or output test, because each is a
   *versioned-behaviour* swap): `cab::SpectrumTap` → `felitronics::analysis::SpectrumTap`
   (`FloatVectorOperations`→`tryPull`); the per-slot JUCE SVF → `felitronics::eq::Svf`;
   `juce::SmoothedValue`/`LinearSmoothedValue` → `felitronics::core::Smoother` (linear→exponential glide).
   → `orbitcab/docs/migration/`.
4. **T6 — FFT-seam + zero-latency partitioned-convolution spike** (the keystone, core-side). Must also
   prove an offline IR resampler (≳60 dB SNR), variable host block size, and RT-safe partition GC — not
   just a raw FFT. → `orbitcab/docs/migration/` (scope) + a core spike.
5. **T7 — OrbitCab cab-path de-JUCE** (`cab::Convolver`/`IRSlot`/`CabEngine` off `juce::dsp`), gated on T6.

## Acceptance (every task)

Green build of **all formats** the plugin ships + its test suite + `auval` (AU). DSP output is
versioned behaviour — any audible change is gated by a golden-vector test, and a deliberate change is
recorded as such. Never leave `main` un-releasable; do the work on a branch + PR.

## Why (the payoff)

Write-once / harden-once / test-once DSP across the product family, and the JUCE-free core that makes
desktop unit-testing, a future WASM demo, and (far future) embedded reachable — without paying their
tax now. Desktop is the gate.
