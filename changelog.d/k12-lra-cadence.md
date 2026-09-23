### analysis · mastering — every short-term series is sampled at 10 Hz, as EBU Tech 3342 requires

`LoudnessMeter` took one 3 s short-term sample a second. Tech 3342 §3.1 asks for more: *"A minimum block
overlap of 2.9 s between consecutive analysis windows (i.e. ≥10 Hz sampling of the loudness level) is
required."* With a 3 s window a 2.9 s overlap **is** a 100 ms step, so the two phrasings are one requirement,
and 100 ms is this meter's hop (exactly, at every rate whose hundredth is a whole number — see the last
paragraph). The requirement entered in V3 (January 2016); the meter had kept libebur128's 1 Hz cadence, a
2 s overlap. It now takes one sample per hop.

**This is compliance with the clause's text, not the repair of a failed conformance test.** Tech 3342's own
Table 1 passes at either cadence — 10.000 / 5.000 / 20.000 / 15.000 LU, well inside the ±1 LU tolerance —
and that table is now in the suite as `tech3342MinimumRequirements()`, cases 1–4, with 5–6 named as not
covered rather than left to be inferred.

**What it changes.** Not "nothing on long programmes", which is what a first measurement on three fixtures
suggested and is wrong. On a square envelope whose states last exactly the window's 3 s, the range reads
20.0 LU at 1 Hz and 9.5 LU at 10 Hz, and the gap does not close with length: 1 Hz answers a flat 20.0 at 12,
20, 30, 60 and 120 s, and 10 Hz answers 9.5 at all of them but 20 s, which is 9.3 — the same figures at 44.1
and 48 kHz. The mechanism is the percentiles. At 1 Hz such a programme offers the window six phases, a sixth
of them sitting wholly inside a loud state and a sixth wholly inside a quiet one,
so P95 and P10 both land on an extreme and the answer is the envelope's full 20 dB swing. At 10 Hz there are
sixty phases, each aligned one is 1.7 % — under the 5 % P95 reaches for — and the fifty-odd mixed windows
that 1 Hz never looked at fill the distribution between them. On that fixture the two cadences agree exactly
away from the window's timescale — 0 LU at 1.5 s states, 20 LU at 5 s and 7.5 s — but that is the fixture, not
a regime: its states are 20 dB apart, so both clear the −20 LU relative gate. Drop the quiet state to −30 dB
and 5 s states read 4.7 LU against 7.7, 7.5 s read 4.7 against 6.9, 10 s read 4.7 against 5.7.

**Cost.** Ten times the short-term entries: for a ten-minute programme the store is 48 KB where it was 4.8.
Every allocation oracle that spells a meter moved with it and says why — the solver's per-pass figure, the
`measure_lra` and `solve` budgets in the C ABI, and the storage table in the conformance suite.

**New:** `LoudnessMeter::shortTermCount()` — the number of short-term samples a range was read from. The
cadence is part of what the standard specifies, and a count is the only way a caller, or a test, can see
which cadence a build runs.

**Fixed, a defect this change would have introduced.** The accepted programme length was bounded at
`INT_MAX - 4`, sized for the block store's `floor(hops) + 4`; the short-term store's `ceil(hops) + 8` then
reaches `INT_MAX + 4`, and the guard `stCount < (int) stE.size()` compares against a narrowed negative, so
every short-term sample is refused from the first one on — a range of 0.0 LU presented as a measurement,
with `droppedBlocks()` reading 0. The bound is now `INT_MAX - 8`, the tighter of the two stores, pinned by
an allocation-free check at the exact boundary. (Unreachable in practice — the store would be 17 GB — and
reachable on paper, which is where a bound lives.)

**New, and a defect that predates this change:** `LoudnessMeter::droppedShortTermSamples()`. The range is
read from its own store, which overflows on its own schedule, and nothing counted the overflow — a programme
fed past its declared length got a range computed over whatever survived, with no way to know. The margin
also narrowed with the cadence: at one sample a second it took an overfeed of about a hundred hops; at one a
hop, 38 are enough. `droppedBlocks()` documented itself as covering the range too, which is half right in a
way worth stating: the block store always fills FIRST, 30 hops earlier, so a caller that checks it reads 0
does have both answers whole — but a non-zero count does **not** mean the range was truncated, and for three
seconds of overfeed it climbs while the range is still perfect. The ordering is now pinned as an ordering —
blocks always first, at both rates, at whole hops AND at a fractional one, and stated in hops rather than in
seconds, because the first draft of that test asked for "3.8 s" and got a length one hop short of what it
meant.

Two things the change made visible rather than caused:

- The C ABI's `measure_lra` fixture swung 1.9 LU on paper but measured 0.5–0.6: its states were 1 s long, so
  the 3 s window never sat inside one and always averaged them. The precondition asking for more than 0.5 LU
  was passing by a single 0.1 LU histogram bin, and the cadence spent that bin. The states are now 4 s, the
  measured range is the programme's real 1.9 LU at **both** cadences. What needs that range is not the check
  it was written under — a bit-equality between the ABI's answer and the core's catches an added offset even
  on a range of zero — but the other defect the block guards against, one that reads channel 0 for every
  plane: channel 0 is the steady tone, so that reads 0.0 LU against a truth of 1.9.
- The LRA-constrained solver search closes this fixture in seven renders where it closed in six. The bracket
  is judged against the measured input range, and that range moved. But the SEARCH still takes six: the
  seventh render is the delivery re-render, bit-identical to the fourth. Before the change the sixth probe
  happened to keep the limit and was delivered as it stood; now it breaks it by a 0.05 LU quantum — the range
  reads on a 0.1 LU staircase and the 2.05 LU limit sits between the two readings that probe can land on. The
  bound is therefore not loosened to seven: the re-render is identified and subtracted, so the six the test
  used to carry (regula falsi, not bisection) stays a search count, and the stale-end halving that makes it
  regula falsi is now asserted directly from the log.

**`ProgrammeReport` moves with it — the same clause governs it.** That class computes EBU Tech 3342 from its
own short-term series rather than reading `LoudnessMeter::loudnessRangeLu()`, so that the gated observation
count is known exactly, and it sampled that series once a second. For the length of this work the core
therefore held **two numbers under one standard's name that disagreed by 10.5 LU** — on the 3 s-state
envelope, at 44.1 and 48 kHz, over 12, 30 and 120 s programmes, `lraLu` read 20.0 where the meter read 9.5.
`kObservationHops` is now 10, both series sample at 10 Hz, and they agree to 0.3 LU again — the tolerance
being summation order (this class oldest-first, the meter newest-first) and nothing else.

Six published fields move with that cadence: `lraLu`, `shortTermP10`, `shortTermP50`, `shortTermP95`,
`shortTermSpreadLu` and `shortTermObservations` — the last by a factor of ten, since it counts observations. A consumer calibrated
against them will see the change; that is the point, not a side effect.

**The guard that should have caught it, and why it did not.** A check asserting the two classes agree to
0.3 LU already existed and passed the whole time. Its fixture is long steady steps, where every window that
is not straddling a step reads the same value and the sampling cadence cannot matter — a guard naming the
right pair and blind by fixture. It is now held on a square envelope with 3 s states as well, which is the
shape that can tell two cadences apart, and the independent oracle beside it writes out BOTH numbers (300
sub-hops for the window, 10 for the step) rather than importing either from the header.

**One boundary moved with it, worth naming.** `ProgrammeReport` refuses with `LoudnessCapacityExceeded`
past `maxDurationSec` — but the store carries a fixed margin behind that line, and the margin is counted in
OBSERVATIONS, so ten times as many observations is a tenth as much slack in seconds. Measured against a 10 s
declaration at 48 kHz: the refusal used to turn at 20 s of programme and now turns at 13.8. The contract did
not change; the undeclared slack behind it shrank to match its own words, from about a hundred hops to 38. A
caller that fed a few seconds over its declaration and got away with it will now be refused. Pinned in hops
at both rates, which nothing did before — the check that stood there was a disjunction ("either fitted or
says it did not") and asserted nothing either way.

**Budgets, for anything that reads them.** `fc_master_need(FC_NEED_MEASURE_LRA)` at the ABI suite's geometry
is **3296 B** where it was 2936; `fc_master_need(FC_NEED_SOLVE)` for 1 s of stereo at the default 1000
buckets is **1 047 848 B** where it was 1 047 704. Struct layouts and the ABI version are untouched — these
are computed, not frozen. `ProgrammeReport`'s short-term store at the default one-hour ceiling is 288 KB
where it was 28.9.

**A compliance caveat the work surfaced, deliberately left alone.** A sub-hop is `lround(0.01·fs)` samples,
so the hop is 100 ms only where a hundredth of the rate is whole. At 48 and 44.1 kHz the cadence is 10 Hz
exactly; at 22050 Hz it is 9.9774 Hz — *under* §3.1's minimum, by 0.23 % — and at 8050 Hz 9.9383 Hz, with the window
rounding to 3.0068 s and 3.0186 s against a specified 3 s. The rounding cuts the other way too: 8049 Hz takes
80.49 down to 80, giving 10.0613 Hz and a 2.9817 s window — fast enough, and too short. That is the meter's
whole sub-hop grid, older than this change, and moving it would move every number the class produces at
those rates. Recorded as its own item rather than touched here: the rates that ship are 44.1, 48 and 96 kHz,
where the cadence is 10 Hz exactly.
