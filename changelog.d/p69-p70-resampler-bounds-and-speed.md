### convolution — the bounds the IR resampler did not have, and the fifteenfold it did not need to cost

**`convolution::resampleIr` is 12-22x faster, and not one output bit moved.** `besselI0` ran once per tap:
37.7 ms for a one-second IR at 44.1 kHz, 164.6 ms at 192 kHz, on the message thread, every cabinet change.
Measured split of that loop: the Bessel series 78% of it, `std::sin` 13%. What replaces it is not a table
and not an approximation — it is an EXACT memo on the fractional phase, so the kernel arithmetic is
untouched and the answer is the same double it always was:

* a tap's weight depends on the output position only through `xx = t - k`, and the window's indices are
  `c + j` for the same offsets `j` about `c = floor(t)`;
* `frac = t - c` is EXACT for every `c >= 0` (at zero it is `t`; above it, Sterbenz), so `t` IS `c + frac`
  and the exact real behind every `xx` is `frac - j` — the SAME real for two outputs that share a `frac`,
  whatever their `c`. IEEE subtraction is correctly rounded, so both produce the same double from it. `xx`
  itself need not be representable; only the sameness of the number being rounded matters.
* `c < 0` is excluded because `frac` is not exact there: at 8 -> 48 kHz outputs 0 and 6 carry the same
  `frac` to the bit and five of their 64 weights still differ in the last place. A fixture in the suite
  holds an input where that difference reaches the float32 result.

It pays because audio rates are small rationals: a one-second 48 -> 44.1 kHz resample has 913 distinct
phases for 44100 outputs and reuses 98% of its windows, and 48 -> 96 kHz has two. A rate with no period at
all (a corrupt file rate, 48000 -> 44101) reuses nothing, and the memo gives up after 4096 misses without a
hit rather than make that case slower — measured 1.0x. Measured after: 2.7 ms at 44.1 kHz, 8.0 ms at
192 kHz. The suite runs every case through a frozen copy of the un-memoized kernel and compares every
output float as bits.

**Bounds that were missing.**

* **`kMaxResampleSamples` (16.7M samples, 64 MiB of float per channel).** The old length gate only kept the
  arithmetic addressable, so a result one sample under `INT_MAX` was 8.6 GB asked of the heap in a single
  call, on the message thread — and a file rate does not have to be absurd to ask: 1.2 Hz against a 48 kHz
  host is a ratio of 40000, and a one-second IR then wants 7.7 GB. An output past the bound is REFUSED, the
  rule P67 already ratified for a result that cannot be addressed. The figure equals
  `MatrixConvolverNupc::kMaxIrSamples` on purpose: the resampler can produce anything the convolver can hold.
* **`IrResampleConfig` has ranges.** `beta` was gated by `isfinite` where a RANGE was meant, so `beta = 1e300`
  passed and the 64-term series returned a silently wrong number (it only overflows to `+inf`, and the window
  to NaN, at about 1.36e4 — long past where its answers stopped being answers). `kMaxBeta = 53.0` is where
  that series stops meeting its OWN convergence test, measured at 53.038057. `kMaxHalfTaps = 4096` bounds one
  output's tap loop, which `halfTaps = 1e9` did not. A value that is not a REQUEST (a NaN, a negative beta, a
  radius below one) is repaired as before; a request past a ceiling is refused, because answering it with a
  smaller kernel would hand back a filter the caller did not ask for.
* The vector overload checks a length before narrowing it to `int`.

**`CabConvolver::prepare` refuses what it cannot honour (law 11(b)), where it used to guess.** A rate that
was not given was answered with the factory 48 kHz and a convolver was then sized from the answer; an
infinite rate made two float-to-int conversions undefined (on arm64 that produced `LLONG_MAX` and a
one-sample "crossfade", i.e. a hard switch). The rate must now be in `(0, kMaxSampleRate = 3e6]` — the house
figure — and `maxIrSeconds` must be a non-negative number; the IR budget saturates in double BEFORE the cast,
through the public `maxIrSamplesFor`, which is total for every argument. A refused prepare now also clears
the pending-retry geometry it can no longer publish: that survived a refusal, and `isBusy()` then answered
true for the life of the object while `flushPending()` could never publish, because it returns on
`! prepared_`. The clear happens before every refusal, so a refused WIDTH (`prepare(..., 4)`) drops a
pending load too, which it did not before — three behaviour changes in this call, all three named.

**An out-of-bounds write in the reference-gain path.** `normalizationGain` took the first second of the IR
as its analysis window but sized its transform at `N`, which stops doubling at 1<<21; above a megasample of
window the two parted company and the copy ran off the end of the buffer — measured under AddressSanitizer
as a 12 MB heap-buffer-overflow WRITE for a 3 000 000-sample IR at a 4 MHz host, reachable through the public
`loadIR`. The window is now the first second OR the transform, whichever is shorter.
