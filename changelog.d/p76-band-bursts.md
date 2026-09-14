<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### analysis — `BandBursts`: bursts in a band, measured against the band's own past

A new offline analyzer, `felitronics::analysis::BandBursts`, finds where one frequency band (5–9 kHz by
default, the corners are parameters) rises above ITS OWN SURROUNDINGS: the coordinates of each burst, how
far above it rose, and how regularly the bursts recur. It is the honest way to talk about sibilance and
harshness without a target curve — an excess over a programme's own past needs no reference the core does
not have, and it survives the two things that destroy a balance measurement, a different genre and a
different overall level.

**No FFT.** A band is a filter, so the band is `high` of an `eq::Crossover2` at the low corner fed into
`low` of a second at the high corner — two LR4 skirts, zero latency, no transform per hop and no second
source of cross-platform divergence. Group delay at the band edges (≈ 0.14 ms, 1.4 % of a 10 ms hop) is
named in the header and deliberately not compensated.

**The baseline is the MEDIAN of a trailing ring of hops, not the mean**, and that is the decision the
module turns on. A mean baseline switches the detector off exactly where bursts are densest: a steady
train of duty cycle `d` over a silent floor reads `1/d`, so it is invisible once `d ≥ 1/10^(enterDb/10)`
— 25.1 % at the default 6 dB, which sixteenth-note hats are past. Three further failures were measured
during the design round: a 50 %-duty pattern reads 3.01 dB and never fires; moving one onset 10 ms
earlier flips it from found to missed; and 151 zeroed hops make the programme resuming at its ordinary
level read as a 6.1 dB burst. A median is immune to all four, raises the duty bound to 50 %, makes that
bound threshold-independent, and — because it SELECTS one of the observed hop energies rather than
summing them — removes the summation-order and drift questions instead of answering them.

**Periodicity is reported, never judged.** Onsets feed two integer histograms: the spacing of adjacent
events, and the bounded-lag autocorrelation of the onset train. The second is not redundant — a hi-hat
with one hit in five missing turns a clean period `P` into `P, 2P` in the first while the second still
peaks at `P`. Both survive event-list exhaustion. The core publishes the counts and the modal spacing; it
does not publish "this is a hi-hat" or "this is a problem".

**Law 8a**, bit-identical under arbitrary re-slicing (same binary, same channel-presence timeline): one
integer clock, the sample loop outside and the channel loop inside, hop boundaries from a counted
`nextHopEnd_` rather than a modulo, thresholds turned into ratios once in `prepare()` so no logarithm
decides anything, and denormal maintenance clocked by `core::StateGrid` instead of by the `process()`
boundary. `maxBlock` sizes nothing. Law 11d: `storageFor()` publishes the whole demand before anything is
allocated, and `process()`/`finish()` allocate nothing.

Two gates that input sanitising cannot cover are named and closed: `eq::Svf` narrows its state to float,
so a **finite** input can overflow it (successive ±3e38 make the 5 kHz high-pass emit −inf on the second
sample), and a hop whose baseline is **exactly zero** — a silent lead-in — would otherwise open an event
at any level and publish +∞ dB. The filter output is gated as well as the input, and a zero-baseline hop
is not judged and is counted. Events also carry whether damage sat in the BASELINE they were measured
against, because a burst can be manufactured entirely by a hole before it.

Absolute powers in the report (`peakPower`, `bandEnergy`, …) are uncalibrated filter-output power: the LR4
pair's passband peak is −3.99 dB at 6791 Hz at 48 kHz and moves with the sample rate (−1.11 dB at
22.05 kHz, −4.66 at 384 kHz), so they are comparable within one prepared instrument and not across rates.
Every ratio is immune to a GAIN — exactly, over the whole range a delivered programme occupies — and only
partly immune to a change of SPECTRUM: when a burst and its baseline have different shapes each side is
weighted by a different point of that rate-dependent dome, which is worth about half a dB between 44.1 and
96 kHz. Both limits are measured and named in the header rather than claimed away.

`fcore_measure bursts <rate> <channels> <raw.f32le>` prints the whole report as raw IEEE-754 bit patterns,
the way `blocks` does, so a future wasm comparison catches a flipped bit that decimal printing would
round away.
