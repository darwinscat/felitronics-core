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

## Optional compiled modules (opt-in via a CMake option)

| Module | Third-party code | Licence | AGPL-compatible |
|---|---|---|---|
| `fftpffft` | **pffft** (Julien Pommier, 2013), derived from FFTPACKv4 (Dr Paul Swarztrauber / NCAR, 1985) | BSD-style (FFTPACK5 / UCAR) | yes |
| `nam` | **NeuralAmpModelerCore** at `b5a68c3ebed5035a91d9207219346c81e8e3ce8e` | MIT | yes |
| `nam` | **Eigen**, supplied by NeuralAmpModelerCore's `Dependencies/eigen` submodule at that pin | MPL-2.0 | yes |
| `nam` | **nlohmann/json**, vendored inside NeuralAmpModelerCore at that pin | MIT | yes |
| `nam` | **namz** `v1.1.1` at `9be9ed8448ab2b72ebfd2605808a4b1f0b24c75a` (`https://github.com/darwinscat/namz`) | MIT | yes |

`felitronics::fftpffft` — the optional SIMD FFT backend, built only with `-DFELITRONICS_WITH_PFFFT=ON`
(default **OFF**, so the default consumer stays header-only with no external code). It vendors the 2-file
original pffft — `modules/fftpffft/pffft/pffft.c` and `pffft.h` — at upstream commit
`09796885cd5b9da5692242de2df0d81e5e1f3d21` (`https://bitbucket.org/jpommier/pffft`). Those two files are
vendored **pristine**: they retain their upstream BSD/FFTPACK header and, by deliberate exemption from the
"every file carries an SPDX header" house rule, are **not** given a felitronics SPDX header — so upstream
diffs stay byte-clean. The FFTPACK/UCAR licence requires the copyright notice be **reproduced in the
documentation of binary distributions**; that obligation flows to any product that ships this backend (e.g.
OrbitCab must carry the pffft notice if it enables the SIMD path).

`felitronics::nam` — the optional compiled NAM inference backend, built only with
`-DFELITRONICS_WITH_NAM=ON` (default **OFF**). NeuralAmpModelerCore and namz are fetched at the pins
above unless local source overrides are supplied; Eigen and nlohmann/json come from the pinned NAM
tree. The target carries NAM's architecture-registration whole-archive link contract transitively.

## Heavy modules (planned — deps land with the module)

These are **not built yet**; listed so the AGPL-compatibility is pre-cleared:

| Module | Planned dep | Licence | AGPL-compatible |
|---|---|---|---|
| `convolution` | FFT backend — **shipped** as the optional `fftpffft` module (pffft; see above). kissfft stays a possible alternative. | BSD | yes |
| `convolution` | (a high-quality offline IR resampler — windowed-sinc, original) | — | n/a |

Watch-list (do **not** add without review): ONNX Runtime, anything GPL-incompatible or with a stricter
copyleft. Each heavy module must also be verified to build under `-fno-exceptions` (some deps use
throwing asserts) before it claims the `wasm-audio` / `embedded-fpu` tiers.
