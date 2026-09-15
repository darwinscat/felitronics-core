### Added

- **The five offline analyzers are callable from JavaScript, and the wasm module's answer is the native
  tool's answer BIT FOR BIT.** `ProgrammeReport`, `SourceForensics`, `HumDetector`, `LowEnd` and
  `BandBursts` reach `tools/wasm/fc_probe.cpp` through one `_run` entry point and caller-owned row
  buffers with mandatory capacities, and each has a parity harness that reproduces
  `fcore_measure <mode>` exactly — a `diff` of the two IS the test. Measured on every mode at 48 and
  44.1 kHz, against the release module and the checked debug module: **zero differing bytes**, over 1302
  lines of raw IEEE-754 bit patterns per rate.

  This is the acceptance the task was written with. It had been weakened to a wasm-only comparison after
  a crew seat measured that byte-exact native-vs-wasm parity was **unattainable** for `hum`, `lowend` and
  `forensics` on the old numerics; `core::det` made the original criterion reachable, so it is the one
  being met.

- **A CI step that enforces it**, five analyzers × two rates × two modules, refusing to count a
  comparison whose output is empty.

### Fixed

- `tools/wasm/build.sh` had no include root for `modules/eq` or `modules/stereo`, which three of the five
  analyzers include from — the module could not compile at all. The same omission on the native side left
  `felitronics_abi_tests` and `felitronics_clips_exposure_tests` linking `felitronics::analysis` when they
  needed `felitronics::analysis_offline`.
- `tools/CMakeLists.txt` carried **four** `target_link_libraries(fcore_measure …)` lines, a residue of
  merging six branches that each added the one they needed. Collapsed to one.

### Notes

- `report` is the only mode whose text is NAMED rather than positional, and its ~100 field names cross
  the ABI **from the module**, produced by the same visitor walk that produced the rows. The visitor is
  deliberately the single enumeration of the report's fields; a name list rebuilt in C++ and again in
  JavaScript would be the second and third copies of it, and the first field added would stop being
  covered without anything failing. The same rule puts `lowend`'s note NAME on the module's side of the
  boundary rather than rebuilding a pitch-class table in JavaScript.
