### Added — `analysis::StereoBandBursts`: which axis carried the burst

The band detector run on **Mid and Side**, with the other axis captured *at each event's own peak hop* and
carried in the same row. Published apart, the two would be lists a caller aligns by index, and the first
off-by-one reads one burst's Side against another burst's Mid — a mistake that looks like a stereo finding
rather than like a bug.

**Absence is structural, not a threshold.** Side is absent when the programme has **one channel**, and only
then. Two bit-identical channels give an exactly zero Side, which the engine already reports as
`eventsValid() == false` with no eligible hop — "explicitly absent" without anyone inventing a dual-mono
predicate. A near-mono programme has a real Side carrying a real residue, and `exactZeroHops()` beside the
band energies lets a caller set its own rule instead of receiving one baked in here.

**The dome is published, because the share alone is not comparable.** The band is an LR4 pair, so for a
pure tone the band-over-wideband ratio peaks strictly below 1 — and not at the band's geometric mean but at
the geometric mean of the *prewarped* corners. With `a = tan(πf_low/fs)` and `b = tan(πf_high/fs)` the
maximum sits at `tan(πf/fs) = √(ab)` and equals `b⁸/(a²+b²)⁴`: **−3.9885 dB at 48 kHz**, −3.8559 at 44.1,
−4.5072 at 96. A caller thresholding "the band holds most of this hop" against a textbook number is wrong
by the dome, and wrong again at another rate. The closed form is checked *through the filter* — a tone is
run and the band/wideband ratio read off the hop trace — and agrees to better than **1e-4 dB**; swapping
the geometric mean for the arithmetic one moves it 0.072 dB and turns the test red.

**The oracle is the scaled-copy theorem, not a second detector.** If `R = k·L` then `M = (1+k)/2·L` and
`S = (1−k)/2·L`: both axes are scaled copies of the same signal, so a detector judging each hop against its
*own* baseline must fire identically on both, at the same hops, whatever `k` is — while the cross reading
must be exactly `((1+k)/(1−k))²` in power. At 0.18 dB of imbalance, small enough that a reader looking at
levels would call the programme centred, that is **+39.69 dB**, and one fixture pins the event logic and
the cross arithmetic at once.

Both roads: `fcore_measure stereobursts` and `fc_probe_stereobursts_*`, byte-identical, gated in CI across
two rates and both modules — with the fixture asserted to fire on **both** axes, since a centred fixture
would let two "nothing on Side" reports diff clean while the cross arithmetic never ran.

### Fixed — a published demand that differed between tiers

`StereoBandBursts::Storage::bytes()` added `2 * sizeof (BandBursts)`. The engines are **members**, not heap
allocations, so that both double-counts and — a `std::vector` holding pointers — publishes a different
number on wasm32 than on a 64-bit host: 3 697 472 against 3 697 328, exactly two engines' worth of three
pointers. A demand a page sizes a heap from cannot depend on which tier answered. Caught by the cross-tier
storage gate on the first run, before the class had shipped.

### Fixed — a published frequency that would differ between rows

`domeHz` was `sampleRate * atan(t) / π`, the obvious spelling. `atan` is a **system transcendental** and is
not the same function on every row, and this number crosses the ABI — so macOS, glibc and wasm could each
publish a different dome. The libm audit refused it. It is now inverted by **bisection on `det::tan`**,
which is strictly increasing on `(0, π/2)` with the band's own corners bracketing the answer: comparisons
and halves only, nothing reaching libm. The result differs from `atan`'s by **one ulp**, and the dome
*share* is bit-identical. `sqrt` stays — IEEE-754 requires it correctly rounded, so it already is the same
everywhere.

### Fixed — two refusal sets that did not match

The K3c parity harness refused an **empty programme** (`malloc(0)` returns 0 and reads exactly like an
allocation failure) and exited **1** rather than 2 on a missing file, while the native road accepted the
first and refused the second with 2. Both are the half of parity a diff of two *successful* runs can never
show. The module also used `planarSpan` where it needed `planarSpanOrEmpty` — the trap that function's own
comment describes, walked into by copying the excursions block, where refusing an empty file *is* the
native behaviour.
