### core · analysis · oversampling · tools — the three libm derivations, closed, and the lint that keeps them closed

v0.33.0 named three places where a number that crosses a platform boundary was still derived through the
system libm, and left them open. All three are closed here, none of them moved a bit of any shipped audio,
and the audit that found them is a gate rather than a reading: `tools/lint/check-det-math.mjs` plus
`tools/lint/det-math-manifest.txt` enumerate every governed transcendental call in the tree — 331 of them
across 61 files — classify each one, and go red when a call is added, removed, or swapped for another.

**The oversampler did not need the frozen tap table the plan assumed.** `PolyphaseOversampler::designFilter`
moves to `det::sin` outright: the taps are narrowed to float, and that narrowing discards 29 of the bits the
libms can disagree about. The system and deterministic designs differ at 18 of 128 taps in double and at
zero of 128 in float, and the final float arrays are byte-identical across Apple clang/arm64, gcc 14/glibc,
emcc/musl and MSVC/UCRT. The margin is not luck: perturbing every `sin()` by a deliberate k ulp leaves all
128 taps unmoved up to k = 2^24, against a real spread of one to three. The tests pin all 128 taps against
the table the OLD path produced on those four rows.

**The FFT needed its seeds converted and nothing else.** `core::offline::fftInplace` takes its stage seeds
from `det::cos`/`det::sin`; the butterfly is left alone, because it was measured not to contract under
`-ffp-contract=on` on any row, with both controls in place (an FMA canary that fuses in the same build, and
a 1-ulp twiddle perturbation that moves the checksum). Pinning O(N log N) products against a flag no shipped
road uses would have been a real cost for nothing. Note the scope of what was OBSERVED, because it is not
the scope of what now holds: the before-and-after hash equality was taken on a 2^16 transform, whose 16 seed
angles are the ones that measurement exercised, while the consumers run at order 17 by default and
`HumDetector`'s AUTO order reaches 21 at 768 kHz. Cross-row agreement at those larger sizes is not observed,
it is by construction — `det` is one implementation compiled into every row — and that is the whole reason
for converting the seeds rather than measuring them again. `core::fft::ScalarRadix2Real` keeps its system seeds by
decision, listed in the manifest as `retain-rt`: that one is the partitioned convolver's, i.e. shipped audio,
and a deterministic route there is an additive class beside it, not a replacement.

**`gainToDb` was worse than described, and the shared function does not move at all.** It is called once per
sample on four paths, worst of them inside `TruePeakLimiter`'s oversampled loop, in the module that is most
of a render, and `det::log10` is 7.4x a system call. So the CONSUMERS move instead: `core::gainToDbDet` sits
beside `core::gainToDb`, sharing the one floor constant, and the programme report, the clip detector, the
reference true-peak meter and the loudness solver call it. `gainToDb` itself is bit-identical and the same
speed before and after.

**Three rounds against the work, and each found something the reading could not.** The gate did not catch
the regression it was written after — reverting `LoudnessSolver::peakDb` to `core::gainToDb` left the lint
green, because carrier calls were only checked inside the deterministic zone; carrier calls now join the
manifest's multiset. The matcher had four evasions and the manifest four phantom entries. The `det::` pin
was six points, which is not a pin: dropping one `volatile` from `exp2Frac` changed `det::pow10` across
100 000 arguments while all six pinned values still matched, so the pin is now a dense per-function checksum
captured on three rows. A NaN clamp written as `std::max (gain, floor)` returns NaN instead of the floor
with the whole suite green. A solver gate pinned only through digital silence let its threshold be raised
tenfold while a peak of 2e-10 made the report read -200 dB and the certificate -193.98.
