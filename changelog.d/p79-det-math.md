### Added

- **`felitronics::core::det` — a transcendental floor that computes the same bits on every row.**
  `cos`, `sin`, `tan`, `log2`, `log10`, `exp2` and `pow10`, each within 2 ulp of the correctly rounded
  value (measured against mpmath at 60 digits, not against a libm), each pinned against FP contraction so
  the answer does not depend on whether the row has an FMA instruction.

### Changed

- **The five offline analyzers and `analysis::SpectrumFrames` now measure on `core::det`, not on the
  system libm.** They were never reproducible across rows: over the ranges they actually use, `pow10`
  differed in 41 % of results between Apple's libm and both Linux ones, `tan` in 35 %, `log10` in 2 %,
  and the Hann window differed in 4032 of 131072 coefficients at order 17 — which every power bin is
  multiplied by. The odd row was not wasm but **Apple**: glibc and musl agree almost everywhere, so the
  developer's Mac computed something neither CI nor the browser did.

  `fcore_measure report|forensics|hum|lowend|bursts` on the same 10 s programme is now **byte-identical
  across Apple clang/arm64, gcc 14/glibc x86-64 and wasm32/musl** — 1299 lines, zero differences, where
  `lowend` alone differed in 28 lines before.

- **`SpectrumFrames` builds one quadrant of its Hann window and reflects the rest by index.** The window
  was never symmetric — `w[i] != w[N-i]` in 10314 of 16383 places, because `2*pi*i/N` and `2*pi*(N-i)/N`
  are different doubles — on the system libm and on `det` alike, which is why no comparison between them
  could ever have shown it. Reflecting makes the symmetry exact and costs a quarter of the calls.

### Notes

- The deterministic path is for COEFFICIENTS, never for a per-sample loop: a window built once in
  `prepare()`, log2(N) twiddle seeds per transform, a handful of thresholds per `setParams`, one `log10`
  per reported value. Measured cost is +0.5 ms on one order-17 window against a 20-30 ms analyzer run.
  The RT modules (`eq::Svf`, `analysis::KWeightingFilter`) are untouched by this entry.
