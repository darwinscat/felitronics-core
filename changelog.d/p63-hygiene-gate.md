### `mastering` — `LoudnessSolver.h` and `DeliveredMastering.h` enter the strict header gate

Both headers were outside `felitronics_header_hygiene`: the solver did not compile under the downstream flag set,
and `DeliveredMastering.h` includes it, so every mastering header built on top would have stayed out too. gcc 14.2
reported seven diagnostics, all in `LoudnessSolver.h`; `DeliveredMastering.h` added none of its own.

- **Six `-Wfloat-equal`** — the "no limit" sentinels (`GainReductionLimit::off()`, the `minPlrDb = -inf` and
  `maxLraLossLu = +inf` tests in `worstExcess()` and `violatedMask()`) and the infeasible tie-break on equal excess
  — now go through `core::exactlyEqual`, like the rest of core. `!=` became `! exactlyEqual`, which is the same
  predicate for every IEEE value, NaN included.
- **One `-Wconversion`, `uint64_t -> int` for `MasterMeasurement::nonFiniteSubHops`** — explicit cast, not a
  defect. The meter in `measure()` is a local, zeroed by its own prepare, fed exactly one `process()` of `frames`
  samples, and its counter moves by at most one per sub-hop of at least one sample, so it cannot exceed `frames`,
  an `int`. The comment at the cast says so, and names what would break it.
- **Proof.** The gate builds warning-free on gcc 14.2 (Debian), Apple clang and emscripten 6.0.9; the full suite
  passes. Behaviour is unchanged by construction: a probe that renders, solves and measures LRA (plain, constrained,
  infeasible, delivered 44.1 -> 48 kHz) prints byte-identical results against `main`, and its `.text` section
  compiled against `main` and against this branch (gcc 14.2 `-O3`) is byte-identical.
