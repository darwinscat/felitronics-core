### Added — K5: what a dynamic EQ band did, read by the (band, lane) pair

Three entry points — `fc_solution_band_gr_stats`, `fc_solution_band_gr_active_stats` and
`fc_solution_band_gr_trace` — and three refusal codes. **No struct grows**, so the version moves for the
entry points alone, as it did at v5, v7 and v9.

**Per lane, not per band**, and that is the core's shape rather than the facade's. The `dyn` block is shared
by a point's lanes, but each lane has its own probe, its own level and therefore its own delta: a band with
Mid and Side both enabled has two different answers at once, and "the band's GR" is not one number.

**The unit is not the compressor's.** This is the depth of the dynamic delta *at the bell's centre*, in dB —
not a change in loudness, and not one multiplier over the whole signal. The statistics carry `|delta|`; the
sign is the sign of `dyn.rangeDb`, which the caller already has, and the delta saturates at exactly that
value (measured: a request of 99 dB reaches 30.0000 and no further, because the core clamps to 30).

### Two halves, and neither is the other

Over the **whole programme**, `activeFraction` says how *often* the band worked. Over the **windows in which
the delta was not zero**, the statistics say how *deep* it went when it did. On this tree's own de-esser
fixture the two read **1.795 dB against 0.315** — a ratio of 5.7 — and the test pins the *ratio*, because the
absolute pair is a property of the fixture and the solver's gain rather than of the statistic.

The pair is **not circular only because both are published**: the active number is conditioned and the
whole-programme `activeFraction` states the condition. A consumer deciding "is a de-esser needed" from the
active p95 alone will always see a busy band.

**The gate is the delta itself**, echoed as `-inf` — which this ABI's own gate documentation defines as
"every window that carried any non-zero input at all". A planted control that gates on something always on
collapses the active half into the whole-programme one: p95 0.315 against 0.315, ratio **1.000000**, all
3751 windows admitted instead of 1916.

**`aboveRange` is structurally zero for a band**, and the doc says so rather than the field implying that
nothing exceeded anything: the distribution spans 0…30 dB because the core clamps `|delta|` to exactly that.
The field is kept because `fc_gr_stats` is one shape shared with the compressor and the limiter.

### Only armed pairs have an answer

"Armed" is the three refusals turned inside out — `dyn.on`, `rangeDb != 0`, the lane enabled — and the
reason is computed **once**, in the solver beside the arming itself, and carried in the solution. A facade
re-deriving it from the parameters would be a second definition of one predicate.

`FC_ERR_BAND_NOT_DYNAMIC`, `FC_ERR_BAND_INERT` (armed at range 0 — a different thing to show a user than
"off") and `FC_ERR_LANE_OFF`; an out-of-range index is `FC_ERR_RANGE`, as `fc_master_eq_dyn_times` reads the
same two arguments. That call deliberately does *not* refuse an off lane, because ballistics exist whatever
the switch says; a statistic of a render does not, and the difference is documented at both.

### The price, taken from the allocator

**80 016 bytes per armed pair** — two histograms over 0…30 dB at 0.01 (24 008 each) and one trace at the
default 1000 buckets (32 000) — measured on a live solve at 80 827. No armed pair costs nothing.

Counted through `operator new` and not `sizeof`, because `sizeof(QuantileHistogram)` is 120 bytes while its
bins are a separate allocation: an estimate built from `sizeof` came out **2500× too small** and nearly
shipped a full 24×5 grid that would have cost 73 MB a solve. And counted over **one pass**, because the
accumulators live inside a render and are freed at its end — a cumulative reading across a two-pass solve
says how much was allocated, not how much is held at once.
