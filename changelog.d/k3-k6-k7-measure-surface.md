### analysis · dynamiceq · ABI — the measurement surface the automatic mastering mode asks for: a band that is an argument, the milliseconds behind a deviation knob, and the alignment stated as a contract

Three items requested by darwinscat.com's `mastering-v2`, which decides the shape of a chain from measurements
of the input and then has to compare its output against that input. Each one closes a question the surface
could not answer, and two of them close a way it could answer WRONGLY.

**K3 — `fc_probe_bursts_run_with`: the same detector, pointed somewhere else.** `analysis::BandBursts` at its
documented 5–9 kHz band answers "sibilance". The same machinery at 80 Hz – 8 kHz with a shorter baseline
answers "how dense are the transients", which decides whether a clipper goes before a limiter or after — and
that is a different question, not a different analyzer. The new entry point takes `bandLowHz`, `bandHighHz`,
`hopMs`, `baselineMs`, `enterDb`, `exitDb`; `maxEvents` is deliberately not exposed, because the counters keep
counting past the event list and `eventsComplete` already says whether the list is whole, so a caller that
could shrink it could only make the list lie about itself. `fc_probe_bursts_storage_bytes_with` prices it: the
baseline ring is `round(baselineMs / hopMs)` hops of `hopMs` each, so both of those move the allocation and a
page sizing itself by the default figure would be short exactly where it asked for a longer memory.

`fcore_measure bursts` takes the same six values as flags, with the same names, defaults and refusals. That is
not a convenience: a road the parity harness cannot drive is a road with no gate, and `fc_probe_bursts_run_with`
would otherwise have shipped with its only evidence being that it compiles. Native and wasm are byte-identical
on the parameterised road, and their refusal sets agree on nine malformed or out-of-domain argument lists.

**A number that was true and was about to become a lie.** The scalar block reported `enterDb` and `exitDb` by
reading the compile-time default, which was correct while one road existed and would have described 6 / 3 dB
while the events came from the caller's numbers the moment a second one did. The module now remembers what the
last successful run installed and reports that.

**`BandBursts::onsetsPerSecond()`** — the whole of the question is the DENOMINATOR, which is why it lives in the
analyzer rather than in each consumer's arithmetic. It is the JUDGED programme, `eligibleHops()` long, not the
file: the first `baselineHops` have no surroundings to be measured against and are not judged, so dividing by
the file's length reports a density over a stretch where no onset could have been found, and reads low on
exactly the short programmes where it matters most.

**K7 — `fc_master_eq_dyn_times`: the milliseconds behind a deviation knob.** `fc_eq_dyn::atk` and `::rel` are
not times. They are deviations in `[0, 1]` around an automatic value the core derives from the band's own fc/Q
(0.5 is auto, 0 four times faster, 1 four times slower), and that automatic value is a field of no struct — so
a caller could not tell "50 ms" from "the slow rail" by looking at what it had written. The call answers what
the follower is actually set to.

It is **per lane, not per band**, which is the core's shape rather than a choice made here: the times come from
the sidechain probe, the probe sits on the lane, and a point's up-to-five lanes carry their own freq/Q while
`dyn` is shared. And it runs `dynamiceq::LaneDynamics::ballisticsFor`, the one expression the chain runs. That
function is new, and it exists because the rails are part of the answer: a lane's freq and Q reach the producers
RAW (`MasteringChain` hands them the caller's parameters, not the band's clamped copy) and are railed inside
`LaneDynamics` to `[10 Hz, 0.49·fs]` and `[0.05, 40]`. A readback that had asked `BandBallistics` with the raw
pair would have answered for 5 Hz where the probe sits at 10. Both audio call sites now go through it too, so
there is one place where the rails live.

A rate below `core::kMinSampleRate` is refused rather than substituted, unlike the audio path, which reads an
unusable rate as 48000 because a `prepare` has already refused one by then. Nothing has refused anything here,
and an answer computed at a rate the caller did not ask about is indistinguishable from one it did.

**K6 — the delivered render is aligned with its input, and the rest is time, not index.** No code changed: what
was missing was the claim, stated where a consumer reads it and pinned where it can fail. `OfflineRenderer`'s
contract is arithmetic — `out[n] = y[n + D]` — so the latency is already off, output sample n is input sample n
processed, and the chain's tail is in the output rather than cut off; a caller that subtracts `latencySamples`
is introducing an error rather than removing one. The other half is not latency at all: on a delivering handle
the two sides have different rates, different lengths and different sample grids, so windows belong in TIME and
an index-for-index comparison measures the resampling ratio and calls it drift.

Three instruments pin it, deliberately of different construction, because one instrument agreeing with itself is
what the first version of this check did while it was wrong: bit-for-bit identity with every stage bypassed,
correlation lag with the stages running, and where a lone impulse lands. The energy centroid was tried as the
second instrument and rejected on evidence — it disagreed by 3.2 ms with the stages running, and that was the
limiter telling the truth (it attenuates ahead of the peak through its lookahead and stays down through its
release, so more of a burst's tail is pulled down than its head). A centroid measures the envelope a dynamics
stage exists to reshape, so it can only answer an alignment question with every stage bypassed, which is where
it is used and where it is exact.

**ABI v9** — one entry point, no struct grew, and therefore no row in the size table, exactly as v7.
