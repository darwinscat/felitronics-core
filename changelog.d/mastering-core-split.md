### mastering · the offline analyzers · the C ABIs — moved to felitronics-mastering-core (BREAKING)

`felitronics::mastering` (the chain, `OfflineRenderer`, the target-loudness solver, delivery at another rate)
and `felitronics::analysis_offline` — the twelve offline programme analyzers: `ProgrammeReport`,
`SourceForensics`, `HumDetector`, `LowEnd`, `BandBursts`, `StereoBandBursts`, `BandCrest`, `PeakExcursions`,
`ClipDetector`, `WaveformPeaks`, `StereoColumns`, `SpectrumFrames` — now live in
[felitronics-mastering-core](https://github.com/darwinscat/felitronics-mastering-core), with their suites, the
two C ABIs (`fc_master`, `fc_probe`), their native CLIs (`fcore_master`, `fcore_measure`) and suites,
`tools/wasm/build.sh` with every native-vs-wasm comparison, and their det-math manifest lines. Target names,
namespaces and header spellings are unchanged (`<felitronics/analysis/ProgrammeReport.h>` included as before),
so a consumer changes its CMake and nothing else: make felitronics-core available first, then
felitronics-mastering-core. A consumer of the wasm modules builds them there.

Core keeps every mastering STAGE (`limiter`, `dither`, `deesser`, `multiband`, `stereo`, `dynamiceq`,
`saturation`), the base meters (`LoudnessMeter`, `TruePeakMeter`, `ReferenceTruePeakMeter`, `KWeightingFilter`,
`CorrelationMeter`, the spectrum panes), `DeliveryResampler`, `StreamResampler`, the `wasm-audio` tier and
`felitronics::test_support`. What changed here for the move:

- `modules/laws/` — NOT a module: the law-11 census (`CallContractTests`) and law 11c (`PauseIsSilenceTests`)
  moved out of `modules/mastering`, with explicit link lists. `tools/ci/scope.py` selects them whenever any
  module they drive changes, derived from their `#include` lines like every other suite. The census's
  OfflineRenderer / MasteringChain case went with the chain: 528 checks here became 522 + 6 there.
- `MathPolicyTests` keeps the SystemMath/DetMath half; the analyzers' six ownership static_asserts moved
  with them.
- `tools/lint/check-det-math.mjs --satellite` runs two whole passes — core's tree against core's lists and
  manifest, then the other repository's against its own (`tools/lint/det-math-zone.txt` there: zone, parity
  entry points, exceptions) — with nothing skipped on either side. Every list entry must name a file that
  exists (`ZONE-ROT`), and an include that resolves in both repositories is refused. Core's zone is now the
  four files the analyzers stand on (`ReferenceTruePeakMeter`, `OfflineFft`, `PolyphaseOversampler`,
  `CascadeOversampler`); core has no parity entry point left.
- `modules/laws/tests/PublishedNumbersTests.cpp` (new): the stages' own law-11d numbers — `storageFor()` against
  what `prepare()` allocates at each stage's clamps, `latencyFor()` against the prepared latency for the
  compressor, limiter and saturator, `EqEngine::objectBytes()`, the delay line's re-clamped tap, `MonoBass`
  `setParams()` — which the mastering chain's suite pinned and nothing else in core did. The meter's
  conformance suite pins its chunk invariance bit for bit on every pre-gate block energy.
- The det-math lint scans every file the parity closure reaches, whatever its extension (a `.inc` included
  from a zone file was read by the closure and audited by nobody).
- CI: the wasm job keeps the tier, the no-threads audit, the long-double gates and the three lints with
  their planted-violation controls (moved onto core's own zone; a new control plants a zone entry for a
  missing file); the wasm-spike build, the native-vs-wasm NULL steps and the storage-price step moved. The
  scope job's `abi` output is gone — nothing in this tree is an ABI any more.

134 -> 101 tests (pffft on): 36 moved — the 34 suites of those modules and tools, plus the two split out of the
law census and the math-policy gate on the way — carrying 10067 of the 30560 checks, and one was added. Every
suite that stayed prints the number it printed before except the census (528 -> 522, the 6 went with the chain)
and the loudness conformance suite (286 -> 294, the bit-exact chunking above).
