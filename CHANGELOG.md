<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

## v0.5.0 — the non-uniform (Gardner) convolver: block-independent, cheaper at small buffers

The v0.4.0 follow-up is delivered. A **non-uniform partitioned (Gardner 1995) convolver** — a time-domain head +
geometrically growing overlap-save FFT stages — replaces the fixed-`P=128` path. It is **block-INDEPENDENT**,
**true sample-zero-latency**, and **cheaper than `juce::dsp::Convolution` across the small, low-latency buffers a
live rig runs** (JUCE only wins the mean at large power-of-two blocks — a theorem of true zero-latency, not a
shortfall). Flat mean **~0.6 %RT on an Apple M5 Pro / ~1.08 % on an Intel i9-13900H** at every DAW buffer from 16
to 4096, from one `prepare()`.

- **feat(convolution):** `NonUniformConvolver<Fft>` — the mono zero-latency NUPC primitive (head `P0` + capped
  octave-doubling to `B_max`, its own frequency-domain delay line per stage).
- **feat(convolution):** `MatrixConvolverNupc<Fft>` — the shipping 2×2 matrix convolver on one raw-L/R history:
  all four routings (mono / LRDiag / MSDiag / Full) + a **click-free 2-slot smoothstep crossfade** for live IR
  swaps. A NULL-verified drop-in for `MatrixConvolver`.
- **feat(lineareq):** the linear- & mixed-phase EQ now convolves on `MatrixConvolverNupc` (A/B transparent — the
  change is CPU + zero latency only).
- **fix(convolution):** removed the long cold-prime crossfade — a cold FDL already yields the exact causal
  convolution, so the ~2.7 s first-activation fade only attenuated correct output and never completed for renders
  shorter than it (offline/short renders came out ~10 dB down). Every swap now uses the short anti-click fade
  (in `NonUniformConvolver`, `MatrixConvolver`, `ConvolutionEngine`).
- **fix(convolution):** `setIr()` broadcasts a mono IR to both channels on a stereo instance (was rejected); a
  `static_assert` pins `state_` lock-free; `maxIrSamples` capped against a stage-offset overflow.
- **docs/tools:** `PERF-NUPC-VS-JUCE.md` (a two-machine 4→8192 fine log-ladder sweep — Apple M5 Pro + Intel i9-13900H,
  adaptive 3–10 warmed reps — rendered to an in-repo SVG chart), `PERF-CONVOLVER-JUCE-GAP.md` (the design ADR);
  `fftbench` head-to-head + `FCORE_FINE_SWEEP`; `tools/plot-convolver-sweep.py`.

Verified by an architecture consilium (per phase), a Popper falsification campaign (differential fuzzer vs
`PartitionedConvolver` + double precision, ASan/UBSan/TSan), and a release-candidate crew review. `ctest` green.

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
