<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

## v0.4.0 — the #1 performance debt resolved

The scalar-FFT / `O(P)` direct-head convolution bottleneck — long-convolution cost that **exploded at host
block 2048+** — is fixed. Cost is now **block-INDEPENDENT**, zero-latency, and JUCE-free.

- **feat(fftpffft):** a new optional, compiled SIMD FFT backend `felitronics::fftpffft::PffftRealFft`
  (`-DFELITRONICS_WITH_PFFFT=ON`, default OFF) — vendored 2-file pffft, hidden-visibility, cross-backend NULL
  parity tested on x86-64/SSE + arm64/NEON incl. ASan/UBSan. (#25)
- **feat(lineareq):** the audio FFT backend is now a template parameter; the design-time FFTs stay pinned to
  the scalar packed-Hermitian layout and the split is **compile-enforced** — a SIMD backend cannot silently
  corrupt a designed FIR. (#24)
- **feat(lineareq):** the convolver partition is **decoupled from the host block** (fixed internal `P=128`).
  A 131072-tap linear-phase EQ is ~2.0 %RT with pffft at every host block (was ~39 %RT @ block 8192 on the
  old scalar path). (#26)
- **feat(convolution):** SIMD-aligned every FFT-seam buffer (`core::fft::SeamAllocator<64>`). (#23)
- **docs:** `PERF-SCALAR-FFT-BOTTLENECK.md` marked RESOLVED.

Known follow-up (tracked separately): the fixed-`P=128` convolver is block-independent but 2–13× more CPU than
`juce::dsp::Convolution` at host blocks ≥ 512; a non-uniform (Gardner) partitioned convolver is planned to
beat it. See `docs/PERF-CONVOLVER-JUCE-GAP.md`.

Prior versions: see the `v0.1.x` – `v0.3.0` git tags.
