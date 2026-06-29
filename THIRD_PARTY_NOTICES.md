<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Third-party notices — felitronics-core

felitronics-core is **AGPL-3.0-or-later**. Every source file carries an SPDX header. New third-party
deps must be AGPL-compatible (BSD / MIT / Apache / MPL-2.0) and are recorded here, per module.

## Light modules (header-only, shipping today)

| Module | Third-party code | Notes |
|---|---|---|
| `core` | **none** | Original code. |
| `analysis` | **none** | Original code (`SpectrumTap`). FFT impl is plugged by the adapter via the seam — its licence is recorded by whoever links it (JUCE `juce::dsp::FFT`, or pffft / kissfft = BSD). |
| `dynamics` | **none** | Original code. |
| `eq` | **none copied** | The matched-filter **method** is from Martin Vicanek's papers — *"Matched Second Order Digital Filters"* (2016) and *"Matched Two-Pole Digital Shelving Filters"* (2024–2025) — **cited inline in `MatchedBiquad.h`**; no third-party code is copied. The Cytomic/Zavalishin TPT SVF (`Svf.h`) likewise implements a published method, not copied code. |

## Heavy modules (planned — deps land with the module)

These are **not built yet**; listed so the AGPL-compatibility is pre-cleared:

| Module | Planned dep | Licence | AGPL-compatible |
|---|---|---|---|
| `convolution` | pffft **or** kissfft (FFT backend) | BSD | yes |
| `convolution` | (a high-quality offline IR resampler — windowed-sinc, original) | — | n/a |
| `neural` | NeuralAmpModelerCore (NAM) | MIT | yes |
| `neural` | Eigen | MPL-2.0 | yes |
| `neural` | nlohmann/json | MIT | yes |

Watch-list (do **not** add without review): ONNX Runtime, anything GPL-incompatible or with a stricter
copyleft. Each heavy module must also be verified to build under `-fno-exceptions` (some deps use
throwing asserts) before it claims the `wasm-audio` / `embedded-fpu` tiers.
