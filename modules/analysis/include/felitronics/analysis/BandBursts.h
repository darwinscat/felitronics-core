// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/DetMath.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/StateGrid.h>
#include <felitronics/eq/Crossover2.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::BandBursts — finds BURSTS of energy in one frequency band (5–9 kHz by default):
// the places where the band rises above ITS OWN SURROUNDINGS, with their coordinates, how far above they
// rose, and how REGULARLY they recur.
//
// WHY A BURST AND NOT A BALANCE. "Too much sibilance", "harshness" and "too much air" are usually asked of
// a target curve — a claim that this programme's 5–9 kHz should sit where some other programme's does.
// That claim needs a reference the core does not have and cannot honestly invent. An EXCESS OVER ITS OWN
// SURROUNDINGS needs none: it is a statement about this programme only, and it survives the two things
// that destroy a balance measurement — a different genre, and a different overall level. The instrument
// therefore measures a RATIO of the band against itself over time, never the band against full scale and
// never the band against a curve.
//
// A CONSEQUENCE, STATED FIRST BECAUSE IT IS THE POINT: the report does not change when the whole programme
// is turned up. Multiplying every sample by a power of two scales every hop energy and the baseline by the
// same exact factor, so every comparison, every coordinate and every dB in the report is bit-identical.
// If a change of level moved a number here, this would be a level meter wearing a burst detector's name.
//   THE ONE EXCEPTION, because an unconditional version of that sentence is false: `core::flushDenormal`
//   zeroes filter state below a FIXED 1e-15, which is not a homogeneous operation, so at amplitudes near
//   and below that floor the two renders diverge. Measured during the review round, on a fixture that is
//   NOT in the suite and is described here so the figure is not an unstated measurement: a 6 kHz
//   programme of 40 hops at amplitude 1e-15, three hops at 1e-14 and ten more at 1e-15, against a 20-hop
//   ring at 48 kHz, reports a peak excess of 20.64179744144726 dB, while the same programme doubled
//   sample by sample reports 20.526561047779559 dB. A second seat measured the underlying divergence on
//   27 of 256 grid alignments. Overflow to inf at the top and float subnormals at the bottom
//   break it for the same reason. Over the whole range a delivered programme actually occupies — roughly
//   1e-9 to 1 — the invariance is exact, and that is what the suite pins.
//
// NO FFT, AND THAT IS DELIBERATE. A band is a FILTER; a spectrum is a detour that costs a transform per hop
// and adds a second source of cross-platform divergence for a number that is then summed back into one
// scalar anyway. The band is `high` of an `eq::Crossover2` at `bandLowHz` fed into `low` of a second at
// `bandHighHz` — two 4th-order Linkwitz-Riley skirts, minimum phase, zero latency by construction
// (`eq::Crossover2::latencySamples() == 0`).
//   BOTH OUTPUTS OF BOTH CROSSOVERS ARE COMPUTED AND TWO ARE DISCARDED. `Crossover2::processSample()` has
//   no one-sided form (`Crossover2.h:52`), so all EIGHT SVF sections run, hold recursive state and are
//   flushed; only four of them are read. That is the price of not writing a fifth private LR4 in this
//   repository. It is paid once per sample per channel offline, and the discarded sections are live state
//   rather than idle — which is why the flush below covers both crossovers whole.
//   GROUP DELAY IS NAMED AND NOT COMPENSATED, and the figure is the DIGITAL one. The analog LR4 formula
//   2·√2/ω_c (90 µs at 5 kHz) understates the prewarped filter, so the cascade's phase was differentiated
//   numerically instead: across 4.5–9.6 kHz the worst delay is 0.160 ms at 44.1 kHz, 0.160 at 48, 0.158
//   at 96 and 0.157 at 192 kHz — nearly rate-independent in TIME, and ≤ 1.6 % of a default 10 ms hop.
//   In SAMPLES it is not small — that same worst case is 7.67 samples at 48 kHz and 30.2 at 192 — so
//   this is a hop-scale argument and not a sub-sample one.
//   AND IT BLOWS UP WHERE THE HIGH CORNER MEETS 0.49·fs. The BLT prewarp stretches without bound as a
//   corner approaches Nyquist: at 18368 Hz — the first INTEGER rate at which the default band is accepted
//   (the bound is 9000/0.49 = 18367.35, and every rate above it is accepted) — the delay at 9 kHz is
//   2.54 ms and the worst in-band value is 3.04 ms, a THIRD of a hop. So the
//   "fraction of a hop" claim holds from about 22 kHz upward and is stated with that bound rather than
//   as a property of the filter. It is still not compensated: the delay is a smooth function of
//   frequency, so there is no single number to subtract, and an event's coordinates are deliberately
//   reported on the INPUT timeline.
//
// THE BAND IS A DOME, NOT A PLATEAU, AND ITS GAIN MOVES WITH THE SAMPLE RATE. This is why nothing here is
// called "band energy in dBFS". The LR4 pair does not multiply the band by one: computed over the cascade
// HP4(5k)²·LP4(9k)², its passband PEAK is −3.99 dB at 6791 Hz at 48 kHz, and 5 kHz and 9 kHz each sit
// 2.59 dB below that peak — the corners are the −6 dB points of their own skirts, which is what an LR4
// crossover is for and is not a flat band. Worse for absolute use, the peak gain is PREWARP-dependent and
// therefore rate-dependent: −1.11 dB at 22.05 kHz, −3.86 at 44.1, −3.99 at 48, −4.51 at 96, −4.66 at
// 384 kHz. A product comparing the absolute band power of a 44.1 kHz master against a 96 kHz one would be
// wrong by 0.65 dB. So: every ABSOLUTE power published here (`peakPower`, `peakBaseline`, `bandEnergy`,
// `tailPartialEnergy`, the trace's `energy`) is UNCALIBRATED FILTER-OUTPUT POWER, comparable within one
// prepared instrument and not across sample rates.
//   A RATIO IS IMMUNE TO A GAIN, AND ONLY PARTLY IMMUNE TO A CHANGE OF SPECTRUM — the weaker claim is the
//   true one. When the burst and its baseline have the same spectral shape the dome cancels exactly and
//   the excess is a pure power ratio. When they do NOT, each side is weighted by a different point of a
//   dome whose shape moves with the rate, and the excess moves with it: a 5 kHz baseline under a 7 kHz
//   burst measures 22.66 dB at 44.1 kHz and 22.20 dB at 96 kHz. Half a dB, from the rate alone. Every
//   gain-invariance claim in this file is exact; no cross-RATE claim is made for a spectral change.
//
// THE RULE.
//   · The programme is cut into HOPS of `hopMs` (10 ms). A hop's ENERGY is the sum over channels and over
//     the hop's samples of the squared band signal — one programme-wide number, because a burst in a mix
//     is one event and not one per channel. Per-channel band energy TOTALS are published separately, so a
//     consumer can still see which channel carries the band.
//   · A hop's BASELINE is the MEDIAN of the `baselineMs / hopMs` hops BEFORE it (200 hops = 2 s), and the
//     hop itself is not in it: the thing being measured is not part of its own surroundings. The hop is
//     DECIDED against the ring and only then WRITTEN into it. What the other order costs is NOT the
//     7.7 dB an earlier draft quoted — that is mean arithmetic, and on a ring of 200 equal energies
//     including a 1000× hop leaves the median exactly where it was, a 0 dB change. What it really does is
//     shift the window by one hop (the write evicts hop k−N early), which moves the median on any
//     material that is not constant. The mutant is caught by the oracle, not by a dB figure.
//   · A hop is ABOVE when `energy > baseline · 10^(enterDb/10)`, and an open burst STAYS open while
//     `energy > baseline · 10^(exitDb/10)`. An EVENT is a maximal run of above-hops. Its start, length,
//     integrated energy, peak hop, the peak's excess in dB, and the WIDEBAND power at that peak are
//     published.
//   · Onsets feed two integer histograms — the spacing of ADJACENT events, and the bounded-lag
//     AUTOCORRELATION of the onset train — which is how a regular hi-hat is told from single flashes.
//     The core publishes the counts. It does NOT publish "this is a hi-hat" or "this is a problem".
//
// WHY THE MEDIAN AND NOT THE MEAN — the decision this module turns on. A mean baseline is a burst detector
// that switches itself off exactly where bursts are densest, and the bound is embarrassingly low: a steady
// train of bursts with duty cycle d over a silent floor has mean = d·peak, so the excess is 1/d and the
// train is INVISIBLE once d ≥ 1/10^(enterDb/10) = 25.1 % at the default 6 dB. Sixteenth-note hats are past
// that. The mean fails four more ways, each found by a different seat of the design round:
//   · 5 loud hops against 5 silent ones — a 50 % duty pattern, i.e. most drum programmes — measures
//     3.01 dB and never fires at all;
//   · moving one onset 10 ms earlier flips it from found to missed, because one burst hop already in the
//     ring lifts the mean by 1.5 % while a ×4 burst clears 6 dB by only 0.5 %;
//   · 151 zeroed hops (a dropout) leave a mean of 0.245·B, and the programme RESUMING at its ordinary
//     level then reads as a 6.1 dB burst — a false event manufactured by the hole before it;
//   · the FIRST hat of a track is measured against a hat-free mean and reads ≈ 20 dB while every later
//     identical hat, measured against a mean the hats themselves lifted, reads ≈ 13 dB. The "worst
//     excess" coordinate then depends on programme ORDER, which for an offline instrument is pure loss.
// The median is immune to all five: a value that is not in the middle of the ring cannot move it. The
// duty-cycle bound rises from 25.1 % to 50 %, and — the part worth the most — it stops depending on the
// threshold at all. Above 50 % duty the bursts ARE the surroundings and the instrument reports no events.
//   AND THERE IS NO WARNING AT THAT BOUNDARY. `burstHops()/eligibleHops()` is a useful duty figure BELOW
//   the bound, but it cannot signal the bound itself: at 50 % duty the median already sits at the high
//   level, so nothing is above threshold, so the detected duty reads ZERO — the same output a clean
//   programme gives. An earlier draft of this comment claimed that ratio was the diagnostic; it is not,
//   and there is no other. That is why the limit is stated here, in the header, rather than left to be
//   detected: a consumer who needs to distinguish "no bursts" from "bursts everywhere" must compare the
//   band's energy against a second instrument, or widen `baselineMs` past the pattern.
//   THE MEDIAN ALSO REMOVES AN ARITHMETIC QUESTION INSTEAD OF ANSWERING IT. It is a SELECTION, not an
//   accumulation: the baseline is always one of the observed hop energies, bit for bit, so there is no
//   summation order to pin, no compensation to argue about and no drift over a million hops. The window
//   sum it replaces would have needed all three — `LoudnessMeter.h:320` recomputes its window from the
//   ring every hop for exactly that reason, and note it walks newest→oldest, one more convention that
//   would have needed its own pinned test (and whose reversal no re-slicing test can see, because both
//   orders are slicing-invariant). None of it is needed here. `baselineHops/2` is the upper median for an
//   even count; there is no interpolation, because interpolating would put a number in the baseline that
//   no hop ever measured.
//
// A STEP CANNOT OUTLIVE THE BASELINE'S ADAPTATION, and that much is a counting theorem. Burst hops enter
// the ring like any others, and the median sits at sorted index `baselineHops/2`, so that index becomes a
// burst hop once ⌈baselineHops/2⌉ of them are in the ring — at which point the ratio falls to about 1 and
// the event ends. For a step to a CONSTANT level well above the exit ratio the length is exactly that:
// 100 hops (1.00 s) at the defaults, and measured 1, 2, 3 hops at `baselineHops` = 1, 3, 5, which is
// ⌈N/2⌉ and not `N/2`.
//   IT IS NOT A BOUND ON EVERY EVENT, and the earlier draft of this comment claimed it was. A signal that
//   keeps RISING holds an event open indefinitely, because the median chases it from behind and never
//   catches up: a 15 %-per-hop ramp over 40 hops against a 20-hop ring is ONE 40-hop event. So the
//   theorem is about a step, the cap is about a step, and a long event means either a step or a ramp —
//   which the event's own energy and peak coordinate distinguish.
//   AND THAT IS ALSO THE SIGNATURE OF A LEVEL STEP, which is the one artifact a local baseline cannot
//   avoid. A file that opens with a quiet lead-in (a fade-in, a dithered floor) and then plays produces
//   ONE event at the start, ⌈baselineHops/2⌉ hops long, with a huge excess; so does any hard section
//   change. A local-baseline instrument cannot tell a step from a burst at the moment of the step — but
//   the LENGTH tells them apart afterwards, and it is published: a sibilant is 5–15 hops, a step is
//   ⌈baselineHops/2⌉ hops and ends by adaptation. A consumer reads the length before the excess.
//
// WHAT IT CANNOT SEE, AND SAYS SO.
//   · A burst pattern at or above 50 % duty, and a step's excess past ⌈baselineHops/2⌉ hops (both above).
//   · A band whose CONTENT is narrower than about 1/hopMs. Resolution here is set by the band-time
//     product: a 4 kHz band over a 10 ms hop has ≈ 2·B·H = 80 degrees of freedom, so a stationary
//     BROADBAND programme's hop energy varies by only 0.6375 dB (1σ) and the default 6 dB threshold is
//     9.41σ away. The χ²₈₀ upper tail at 10^0.6 of the mean is 3.4e-30 — a COMPUTED figure, quoted as a
//     scale and not as a measurement: the suite's stationary-bright test asserts zero events over some
//     hundreds of judged hops, which is all a test of that length can assert. Narrow the
//     content to ≈ 20 Hz and the hop energy becomes nearly exponential, a few percent of hops exceed
//     6 dB, and the detector fires. It is not wrong to: a 20 Hz-wide band at 6 kHz genuinely fluctuates
//     by more than 6 dB several times a second. That is a tremolo, and the instrument reports the tremolo.
//   · The first `baselineHops` hops have no baseline and are not judged. With the defaults the first
//     eligible hop is hop 200, which ENDS at 2.01 s — so a programme needs `(baselineHops + 1)` whole
//     hops before anything is looked at, which is why the reason below is `NoEligibleHop` and not
//     "shorter than the baseline".
//   · An open hi-hat whose decay never falls back to within `exitDb` OF the median merges its hits into
//     one event — the stay test is `energy > baseline · 10^(exitDb/10)`, so an event runs while the band
//     is still ABOVE the median by more than `exitDb`, and equality closes it —
//     the hysteresis that stops one sibilant SPLITTING is the same mechanism that MERGES a fast roll.
//     The event's length and `hops` are the evidence that it did.
//   · Order. The adjacent-spacing histogram is a distribution, so a pattern and a shuffle of the same
//     spacings read alike. The autocorrelation recovers SOME of what that loses but not all of it, and
//     in particular it is REVERSAL-BLIND: onsets at {0,10,30,60} and {0,30,50,60} have the same multiset
//     of pair distances {10,20,30,30,50,60} and therefore identical adjacent histograms AND identical
//     lag histograms. The full onset sequence is in the event list's `start` fields for as long as the
//     list holds, and that is the only place order survives.
//   · Tempo drift, in `modalIntervalMass()`: the window is ± ONE HOP, so its tolerance is ±1/P of the
//     period — ±8.3 % at a 12-hop period, ±0.2 % at 512. A hand-played hat whose jitter exceeds one hop
//     spreads past the window and reads as less periodic than it is. The full histogram and the
//     autocorrelation are published because that scalar alone would lie.
//
// NO EVIDENCE IS NOT A ZERO. Two situations produce no number rather than a reassuring one:
//   · No hop was ever eligible — `eventsValid()` is false with reason `NoEligibleHop`. An empty event
//     list then means "not looked at", not "clean".
//   · A hop whose baseline is EXACTLY zero — more than half its ring is digital silence — is not judged.
//     An excess over nothing is not a number. Without this gate the first non-zero sample after a silent
//     lead-in opens an event at ANY level, and its peak then FREEZES on the silence boundary, because a
//     zero denominator makes the strict cross-multiplication below reject every later hop however much
//     worse it is. (The published dB would NOT be +∞ — `closeEvent`'s own `peakDen_ > 0.0` guard writes
//     0.0 — which is if anything worse: a fabricated event reporting 0.0 dB of excess.) Counted in
//     `zeroBaselineHops()`.
// An invalid quantity reads 0.0. NaN is never used to mean "no value", and no published double is ever
// non-finite.
//
// NON-FINITE INPUT, AND THE PART THAT IS NOT ABOUT THE INPUT. A non-finite sample is counted, its first
// coordinate recorded, and a canonical 0.0f substituted BEFORE the filters, so nothing non-finite reaches
// an IIR, a square or an accumulator. That alone is NOT enough here: `eq::Svf` narrows its state to float
// (`Svf.h:91`) and a FINITE input can overflow it — successive +3e38f, −3e38f drive the 5 kHz HP4
// non-finite on the SECOND sample (its first 2nd-order section emits −inf and the second turns that into
// NaN), which `Svf.h:100` documents for 3e38 and `Svf.h:134` for a sign-alternating 1e37. So the FILTER
// OUTPUT is gated too, and a hop that hit either gate is DAMAGED: counted, flagged on any event that
// contains one, and flagged on any event whose BASELINE ring held one (`baselineTouchedNonFinite` — a
// burst can be manufactured entirely by damage that happened before it, so damage to the REFERENCE has to
// be visible on the event and not only in a global counter). The filter state is healed at the next
// `core::StateGrid` boundary; the hop and ring values an inf would have poisoned are not healable, which
// is why they are gated rather than flushed.
//   AND THE FLUSH IS THE GRID'S ALONE — no `healPoison()` at the end of `process()`. Law 8a(a) says the
//   poison half of a flush runs at the end of every call as well, and for an AUDIO path that is right: a
//   NaN's only quality is how soon it goes, and law 8a explicitly stops promising invariance "once a
//   filter has overflowed" (`docs/DSP-ARCHITECTURE.md:162`). A READ-ONLY INSTRUMENT trades the other way.
//   Recovery timed by the call boundary would make the samples after an overflow a function of where the
//   caller cut, so the set of DAMAGED hops — a published quantity — would stop being a property of the
//   programme. Grid-only recovery is bounded by 64 samples, lands on the same absolute sample under every
//   slicing, and every sample it fails to save is gated and counted anyway. The deviation is deliberate,
//   and `eq::Crossover2::resetChannel()` (`Crossover2.h:36`) is likewise never called: a channel that
//   stops being fed is handed the silence it is receiving (law 11a), because clearing its four columns
//   would put the caller's channel-count timeline into the filter state.
//
// LAW 8a — THE REPORT IS BIT-IDENTICAL UNDER ARBITRARY RE-SLICING. Same input bits, same parameters, same
// binary, and the same per-sample map of which channels are present: every field of the report and every
// hop of the trace is bit-for-bit the same however the caller cuts the stream into `process()` calls.
// What buys it: one integer clock (`totalSamples_`), the sample loop OUTSIDE and the channel loop INSIDE
// (the other order sums the shared hop accumulator as c0i0 c0i1 c1i0 c1i1 on one call and c0i0 c1i0 c0i1
// c1i1 on two — different order, different bits), hop boundaries from a counted `nextHopEnd_` and never
// from `pos % hop`, a baseline that is a selection rather than a sum, thresholds turned into ratios once
// in `prepare()` so no logarithm takes part in a decision, and denormal maintenance clocked by
// `core::StateGrid` — never at the end of `process()`, which is wherever the CALLER chose to cut.
// `maxBlock` sizes nothing at all and may not select a window, a hop or a branch.
//   WHAT IS NOT CLAIMED, precisely. (1) Bit-identity ACROSS platforms, toolchains or contraction modes.
//   The tree builds with `-ffp-contract=on` and the wasm tier with `off`; inside `Svf::processSample`
//   `a1*ic1 + a2*v3` is a contractible form with a 53-bit `a1`, so a band sample can differ between rows
//   by a fused multiply-add. Every product this file forms is nevertheless STORED before it is added, per
//   law 10's blanket pin (`analysis::StereoSums::add`, `StereoColumns.h:35`) — a product of two float32s
//   is exact in double and could not have fused into a different number, but the rule reads "every
//   product" so that nobody has to re-derive which ones are safe. (2) A channel-count timeline is the
//   caller's, not the audio's (`StateGrid.h:34`): re-slicing must preserve which channels are present at
//   each absolute sample, and "omit channel 1" is a different programme from "feed channel 1 zeros".
//   (3) WAS a caveat about `enterDb`/`exitDb` becoming ratios through `std::pow`, which is not correctly
//   rounded on every libm, so a 1-ulp difference could move a hop sitting exactly on the threshold. Since
//   v0.33.0 they go through `core::det::pow10` instead, which IS the same function on every row — and it
//   mattered more than it looked: over the dB range these thresholds live in, 41 % of `std::pow(10, x)`
//   results differ between Apple's libm and both Linux ones (measured). A threshold that moves by an ulp
//   does not move a number by an ulp; it flips a decision.
//
// FORM. setParams / prepare / process / finish / reset, the shape of `analysis::ClipDetector`. Every
// parameter takes effect at the next `prepare()`. `prepare()` disarms, validates, sizes itself through
// the public `storageFor()` and only then allocates (law 11d); `process()` and `finish()` allocate
// nothing. `finish()` is idempotent, closes an open event at coordinate T, publishes the uncovered tail
// as a number, and freezes the report; `process()` refuses until `reset()`.
struct BandBurstsParams
{
    double bandLowHz   = 5000.0;    // the band, as a pair of Linkwitz-Riley corners
    double bandHighHz  = 9000.0;
    double hopMs       = 10.0;      // the observation
    double baselineMs  = 2000.0;    // the surroundings; baselineHops = round(baselineMs / hopMs)
    double enterDb     = 6.0;       // a hop opens a burst above this much over its baseline
    double exitDb      = 3.0;       // and the burst stays open above this much — see the note below
    int    maxEvents   = 1 << 14;   // capacity of the event list; the counters keep counting past it

    // WHY TWO THRESHOLDS. One dipping hop inside one sibilant would otherwise split it into two events,
    // and a split is not cosmetic: the spacing histograms read the GAPS, so an event cut in half invents
    // a short interval that never happened. 6 dB / 3 dB are chosen for the default band, not derived —
    // a sibilant or a cymbal in a mix sits 6–15 dB over its 2 s surroundings, while broadband programme
    // material modulates that band by well under 3 dB (see the band-time product above). They are POWER
    // ratios, hence 10^(dB/10) and not 10^(dB/20).
};

// Why a quantity has no number. `None` is the only value that means the quantity is present.
enum class BandBurstsInvalid : std::uint8_t
{
    None           = 0,
    // No hop was ever judged. TWO causes, and `zeroBaselineHops()` tells them apart: the programme never
    // reached hop `baselineHops`, OR it did and every baseline was exactly zero (digital silence in more
    // than half of every ring — an all-zero programme of any length reports this reason for ever).
    NoEligibleHop  = 1,
    NonFiniteInput = 2    // a non-finite sample, or an overflowed filter, occurred in the programme
};

// One burst. Coordinates and continuous evidence; no verdict.
struct BandBurst
{
    std::int64_t start        = 0;      // first sample of the first above-hop, in samples since reset()
    std::int64_t length       = 0;      // samples; whole hops, except an event closed by finish() (to T)
    std::int64_t peakAt       = 0;      // first sample of the hop with the greatest excess
    double       peakPower    = 0.0;    // that hop's mean square, raw and linear (uncalibrated — see above)
    double       peakBaseline = 0.0;    // the baseline it was measured against, same units
    double       peakExcessDb = 0.0;    // 10·log10 of the ratio, formed ONCE from the accumulators
    double       peakWidePower = 0.0;   // WIDEBAND mean square of that same hop, before the filters
    double       energy       = 0.0;    // sum of the event's hop energies, raw and linear — the burst's dose
    std::int64_t hops         = 0;
    // Damage anywhere inside this event's published extent [start, start+length): one of the hops it
    // counted, OR the partial tail an event closed by finish() reaches over. So it can be true while
    // every one of the `hops` full hops was clean — the extent is what it describes, not `hops`.
    bool         touchedNonFinite         = false;
    bool         baselineTouchedNonFinite = false;  // a damaged hop sat in the baseline this event was judged against
    bool         closedByFinish           = false;  // the stream ended while the event was open
};

// One closed hop, handed to an observer AT THE MOMENT IT CLOSES — never at the exit of process(). This is
// the band-energy curve as the instrument actually sees it, and it is also what makes law 8a testable:
// comparing only the final report is not enough. This project has its own measured case —
// `docs/LAW8-KWEIGHTING.md:78`, where the final LUFS and dBTP matched bit for bit while 5 of 97
// intermediate block energies had moved — and a maximum, a percentile, a gate and a histogram can all
// absorb a re-slicing defect. Everything here is RAW and LINEAR, before any log, ratio or threshold.
//
// WIDEBAND IS HERE FOR THE QUESTION THE RATIO CANNOT ANSWER. "The band rose above its own past" does not
// separate a hat that is loud only for a quiet passage from a sibilant in a full mix; the band's share of
// the broadband energy does. It is a second coordinate, not a second verdict — the core never compares
// the two.
struct BandBurstsHopTrace
{
    std::int64_t hopIndex    = 0;
    std::int64_t startSample = 0;
    std::int64_t endSample   = 0;       // exclusive
    double       energy      = 0.0;     // the hop's summed squared BAND signal, raw
    double       wideEnergy  = 0.0;     // the same over the unfiltered (finite-substituted) input
    double       baseline    = 0.0;     // the median of its ring, raw; 0 when it was not eligible
    bool         full        = false;   // false only for the single partial hop finish() closes
    bool         eligible    = false;   // the ring was full AND its median was above zero
    bool         damaged     = false;   // a non-finite sample, or an overflowed filter, in this hop
    bool         above       = false;   // the decision taken for this hop (enter or stay, as applicable)
    bool         inEvent     = false;   // event state AFTER this hop
    std::int64_t eventStart  = -1;      // the open event's start, or -1
};

// PRECONDITION: an observer MUST NOT call back into the detector it was called from. It is a read-only
// notification delivered from inside the sample loop, at a point where the hop it describes has not yet
// been committed — a `finish()` from the callback re-enters hop closing and publishes the just-finished
// full hop a second time as a partial one, and a `process()` from it would interleave two sample loops.
// Copy what you need and return.
using BandBurstsHopObserver = void (*) (void* user, const BandBurstsHopTrace& hop);

class BandBursts
{
public:
    static constexpr int    kIoiBins         = 512;     // adjacent-onset spacings, one hop per bin
    static constexpr int    kMaxLag          = 512;     // onset autocorrelation, lags 1..512 hops
    static constexpr double kMinSampleRate   = core::kMinSampleRate;   // P51: the core's floor, one number
    static constexpr double kMaxSampleRate   = 768000.0;
    static constexpr double kMinHopMs        = 1.0;
    static constexpr double kMaxHopMs        = 1000.0;
    static constexpr int    kMaxBaselineHops = 1 << 16;
    static constexpr int    kMaxEventsLimit  = 1 << 24;
    static constexpr double kMaxThresholdDb  = 200.0;   // keeps 10^(dB/10) finite and refuses +inf
    // `Svf::setParams` clamps a corner into [1, 0.49·fs] SILENTLY (`Svf.h:61`). Refusing the same range
    // here is what keeps the filter from measuring a band the caller never asked for: at 16 kHz the
    // default HIGH corner alone would clamp from 9000 Hz to 7840 (the 5 kHz corner is already under it),
    // leaving a 5 kHz–7.84 kHz band reported as if it were 5–9 kHz. The admission condition is
    // `bandHighHz <= 0.49·fs`, so the default band needs fs >= 9000/0.49 = 18367.35 Hz — and 18368 Hz,
    // the first integer above it, is accepted.
    static constexpr double kMaxCornerFraction = 0.49;

    // WHAT prepare() ALLOCATES, from the function prepare() sizes and validates itself with (law 11d).
    // `ok == false` on exactly the arguments prepare() refuses, and then every count is zero.
    struct Storage
    {
        bool          ok             = false;
        int           hopSamples     = 0;
        int           baselineHops   = 0;
        std::size_t   ringEntries    = 0;    // the baseline ring
        std::size_t   scratchEntries = 0;    // the copy nth_element selects the median in
        std::size_t   flagEntries    = 0;    // one damage flag per ring slot
        std::size_t   eventEntries   = 0;
        std::size_t   channels       = 0;
        // The demand of THIS prepare(). It is not a running footprint: `std::vector::assign` never
        // releases capacity, so after a second, smaller prepare() the object still owns the larger
        // buffers while `bytes()` reports the smaller demand. Law 11d is satisfied either way —
        // prepare() never allocates more than it published — but the number is a demand, not a census.
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) ringEntries    * sizeof (double)
                 + (std::uint64_t) scratchEntries * sizeof (double)
                 + (std::uint64_t) flagEntries    * sizeof (std::uint8_t)
                 + (std::uint64_t) eventEntries   * sizeof (BandBurst)
                 + (std::uint64_t) channels       * sizeof (Channel);
        }
    };

    static Storage storageFor (double sampleRate, int maxChannels, const BandBurstsParams& p) noexcept
    {
        Storage st;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return st;   // NaN fails too
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        if (p.maxEvents < 0 || p.maxEvents > kMaxEventsLimit) return st;
        if (! (p.hopMs >= kMinHopMs && p.hopMs <= kMaxHopMs)) return st;
        // The upper bound is not decoration: `llround` below is UNDEFINED outside int64, and
        // baselineMs = 1e30 (or DBL_MAX, or +inf) reaches it. Bound the input, then convert.
        if (! (p.baselineMs >= p.hopMs && p.baselineMs <= kMaxHopMs * (double) kMaxBaselineHops)) return st;
        // `enterDb > 0.0` is TRUE for +inf, and an infinite ratio disables detection for ever behind a
        // report that still reads valid. Both thresholds must be finite and bounded so that
        // 10^(dB/10) is a finite number.
        if (! (p.enterDb > 0.0 && p.enterDb <= kMaxThresholdDb)) return st;
        if (! (p.exitDb >= 0.0 && p.exitDb <= p.enterDb)) return st;
        if (! (p.bandLowHz >= 1.0) || ! (p.bandHighHz > p.bandLowHz)) return st;
        // Validated on the FLOAT-ROUNDED corners, because those are the ones the filter is given
        // (`Crossover2::setFrequency` takes a float). Checking the doubles would let a corner round up
        // past 0.49·fs into `Svf`'s silent clamp, and would accept 5000 / 5000.0001 as a band when both
        // round to the same float and the band is empty.
        const float loF = (float) p.bandLowHz, hiF = (float) p.bandHighHz;
        if (! ((double) loF >= 1.0) || ! ((double) hiF > (double) loF)) return st;
        if (! ((double) hiF <= kMaxCornerFraction * sampleRate)) return st;

        // Both of these are checked against 1 even though the bounds above already imply it. A zero hop
        // would make `nextHopEnd_` unreachable and close NO hop for the whole programme — a green,
        // valid-looking, all-zero report — and a zero baseline count would make `k % baselineHops` an
        // integer division by zero. Neither may depend on arithmetic elsewhere staying as it is.
        const std::int64_t hop = std::llround (sampleRate * p.hopMs / 1000.0);
        if (hop < 1) return st;
        const std::int64_t nb = std::llround (p.baselineMs / p.hopMs);
        if (nb < 1 || nb > kMaxBaselineHops) return st;

        st.ok             = true;
        st.hopSamples     = (int) hop;
        st.baselineHops   = (int) nb;
        st.ringEntries    = (std::size_t) nb;
        st.scratchEntries = (std::size_t) nb;
        st.flagEntries    = (std::size_t) nb;
        st.eventEntries   = (std::size_t) p.maxEvents;
        st.channels       = (std::size_t) maxChannels;
        return st;
    }

    void setParams (const BandBurstsParams& p) noexcept { params_ = p; }   // all of it: at the next prepare()
    const BandBurstsParams& params() const noexcept { return params_; }

    // The hop observer is not a structural parameter and sizes nothing; it may be set at any time and is
    // not touched by reset(). Passing nullptr detaches it.
    void setObserver (BandBurstsHopObserver fn, void* user) noexcept { observer_ = fn; observerUser_ = user; }

    // ARMEDNESS IS DERIVED FROM THE STORAGE, not from the flag alone. The implicitly generated MOVE
    // steals the vectors and COPIES `prepared_`, so a moved-from detector would otherwise stay armed and
    // index an empty `chans_` on the next `process()` — the defect a code-review round found here, and
    // the same shape as the moved-from chain deref already known elsewhere in this tree. Asking the
    // storage instead makes a moved-from object simply unprepared, which is what it is. `events_` is
    // deliberately not consulted: `maxEvents == 0` is a legal configuration with an empty event list.
    bool isPrepared() const noexcept { return prepared_ && ! chans_.empty() && ! ring_.empty(); }

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        // Law 11b: disarm, validate, write — AND DISARM MEANS THE REPORT TOO, the way SourceForensics
        // already spells it. A refused prepare() after a finished measurement used to leave isFinished()
        // answering true and the previous report still readable, so an instance reconfigured with bad
        // arguments kept certifying the programme before it. Found by the release round; P71 had already
        // named the same shape a defect.
        prepared_ = false;
        finished_ = false;
        channels_ = 0;
        const Storage st = storageFor (sampleRate, maxChannels, params_);
        if (! st.ok) return false;

        sampleRate_   = sampleRate;
        channels_     = maxChannels;
        hopSamples_   = st.hopSamples;
        baselineHops_ = st.baselineHops;

        // The thresholds become RATIOS here, once, so that no logarithm ever takes part in a decision.
        enterRatio_ = core::det::pow10 (params_.enterDb / 10.0);
        exitRatio_  = core::det::pow10 (params_.exitDb  / 10.0);

        xLow_.prepare  (sampleRate, maxChannels);
        xHigh_.prepare (sampleRate, maxChannels);
        xLow_.setFrequency  ((float) params_.bandLowHz);         // Crossover2 carries its corner as float;
        xHigh_.setFrequency ((float) params_.bandHighHz);        // the report publishes the EFFECTIVE value

        ring_.assign      (st.ringEntries,    0.0);
        scratch_.assign   (st.scratchEntries, 0.0);
        ringFlags_.assign (st.flagEntries,    (std::uint8_t) 0);
        events_.assign    (st.eventEntries,   BandBurst {});
        chans_.assign     (st.channels,       Channel {});

        prepared_ = true;
        reset();
        return true;
    }

    // Re-anchors EVERYTHING: the clock, the filters, the grid, the ring and its damage bookkeeping, the
    // event machine, both histograms, every counter and every reason.
    void reset() noexcept
    {
        for (auto& ch : chans_) ch = Channel {};
        std::fill (ring_.begin(), ring_.end(), 0.0);
        std::fill (ringFlags_.begin(), ringFlags_.end(), (std::uint8_t) 0);
        for (auto& b : onsetBits_) b = 0u;
        for (auto& v : ioi_) v = 0;
        for (auto& v : lag_) v = 0;

        xLow_.reset(); xHigh_.reset(); grid_.reset();

        totalSamples_ = 0;
        nextHopEnd_   = (std::int64_t) hopSamples_;               // hopSamples_ >= 1 whenever prepared_
        hopCount_ = 0; hopSumSq_ = 0.0; hopWideSumSq_ = 0.0; hopDamaged_ = false;
        eligibleHops_ = 0; zeroBaselineHops_ = 0; burstHops_ = 0; damagedHops_ = 0; damagedInRing_ = 0;
        overflowSamples_ = 0; firstNonFiniteAt_ = -1;
        eventCount_ = 0; intervalCount_ = 0; ioiOverflow_ = 0; onsetCount_ = 0;
        lastOnsetHop_ = 0; haveLastOnset_ = false;
        inEvent_ = false; evStart_ = 0; evHops_ = 0; evEnergy_ = 0.0;
        evTouched_ = false; evBaseTouched_ = false;
        peakNum_ = 0.0; peakDen_ = 0.0; peakWide_ = 0.0; peakAt_ = 0;
        tailPartialSamples_ = 0; tailPartialEnergy_ = 0.0;
        finished_ = false;
    }

    static constexpr int latencySamples() noexcept { return 0; }          // a read-only sink

    // READ-ONLY. Law 11: malformed → unprepared → finished → nch > maxChannels → n == 0 → run.
    // A legal call of any length is consumed IN FULL; a refused one consumes nothing. `nch == 0` is a
    // legal call in which every prepared channel is absent, and it still spends its samples of audio time.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! isPrepared() || finished_) return false;
        if (numChannels > channels_) return false;
        if (n == 0) return true;

        for (int i = 0; i < n; ++i)
        {
            // Sample OUTSIDE, channel INSIDE — law 8a's summation order for the shared hop accumulators.
            for (int c = 0; c < channels_; ++c)
            {
                Channel& ch = chans_[(std::size_t) c];
                float x = 0.0f;
                if (c < numChannels)
                {
                    const float raw = in[c][i];
                    if (std::isfinite (raw)) { x = raw; }
                    else
                    {
                        ++ch.nonFinite;
                        if (firstNonFiniteAt_ < 0) firstNonFiniteAt_ = totalSamples_;
                        hopDamaged_ = true;                       // canonical 0.0f goes on into the filters
                    }
                }
                else { ++ch.absent; }                             // law 11a: a hole, and it is counted

                // Law 10's pin: every product is STORED before it is added. A product of two float32s is
                // exact in double and could not fuse into a different number, but the rule reads "every
                // product" so nobody re-derives which ones are safe (`StereoColumns.h:35`).
                volatile double prod;
                prod = (double) x * (double) x; hopWideSumSq_ += prod;

                float low = 0.0f, high = 0.0f, band = 0.0f, discard = 0.0f;
                xLow_.processSample  (c, x,    low,  high);       // high = HP4 above bandLowHz
                xHigh_.processSample (c, high, band, discard);    // band = LP4 below bandHighHz

                // The second gate, and the one input sanitising cannot stand in for: Svf's float state
                // can overflow from a FINITE input (`Svf.h:100`). One inf here would make this hop, and
                // then the whole ring, NaN for ever, and no flush can undo a poisoned accumulator.
                if (std::isfinite (band))
                {
                    prod = (double) band * (double) band;
                    hopSumSq_ += prod;
                    addCompensated (ch.energy, ch.energyComp, prod);
                }
                else { ++overflowSamples_; hopDamaged_ = true; }
            }

            ++totalSamples_;                                      // exactly once, after every channel of t

            // Maintenance on AUDIO time. StateGrid::advance returns true when the sample just consumed
            // ENDED a period, so the flush lands on the same absolute sample however the caller sliced.
            if (grid_.advance (1)) { xLow_.flushDenormals(); xHigh_.flushDenormals(); }

            // Frame geometry without a modulo: hop k is exactly [k·H, k·H + H).
            while (totalSamples_ == nextHopEnd_)
            {
                closeHop ((std::int64_t) hopSamples_, true);
                nextHopEnd_ += (std::int64_t) hopSamples_;
            }
        }
        return true;
    }

    // End of stream. Closes the one partial hop by its ACTUAL length, closes an open event at T, and
    // freezes the report. Idempotent.
    void finish() noexcept
    {
        if (! isPrepared() || finished_) return;
        const std::int64_t covered = hopCount_ * (std::int64_t) hopSamples_;
        tailPartialSamples_ = totalSamples_ - covered;
        if (tailPartialSamples_ > 0)
        {
            // The partial hop is NAMED, not judged. It does not enter the ring (a short hop among
            // equal-duration ones is a silent re-weighting of the median) and it cannot open an event
            // (a fraction of the observation the threshold was set on is not that observation). Every
            // sample it holds is already inside the per-channel totals and every counter.
            tailPartialEnergy_ = hopSumSq_;
            // An event still open here will be closed at T, so the partial tail lies INSIDE its
            // published extent; its damage has to reach the event's flag, which the full-hop path
            // below can never do for a hop that is not full.
            if (hopDamaged_ && inEvent_) evTouched_ = true;
            closeHop (tailPartialSamples_, false);
        }
        if (inEvent_) closeEvent (totalSamples_, true);
        finished_ = true;
    }

    // --- the report (partial while streaming, final after finish()) ---
    // Every accessor is total: an index outside the report answers an empty record or zero rather than
    // reading past a vector.
    bool isFinished() const noexcept { return finished_; }
    std::int64_t samplesProcessed() const noexcept { return totalSamples_; }
    int channels() const noexcept { return channels_; }
    double sampleRate() const noexcept { return sampleRate_; }
    int hopSamples() const noexcept { return hopSamples_; }
    int baselineHops() const noexcept { return baselineHops_; }
    // The EFFECTIVE corners: Crossover2 carries its frequency as float, so this is the band that was
    // actually measured, not the double that was asked for. ONLY MEANINGFUL ONCE PREPARED: before
    // prepare() (and on a moved-from object, which `isPrepared()` reports as unprepared) these return
    // `Crossover2`'s own 1000 Hz default for BOTH corners — a low == high band that `storageFor` would
    // itself refuse. Ask `isPrepared()` before reading the geometry.
    double bandLowHz() const noexcept { return (double) xLow_.frequency(); }
    double bandHighHz() const noexcept { return (double) xHigh_.frequency(); }

    std::int64_t hopCount() const noexcept { return hopCount_; }
    std::int64_t eligibleHops() const noexcept { return eligibleHops_; }
    std::int64_t zeroBaselineHops() const noexcept { return zeroBaselineHops_; }
    std::int64_t burstHops() const noexcept { return burstHops_; }
    // Hops that held a non-finite input sample or an overflowed filter output. NOTE: the partial tail
    // counts here if it was damaged, and it is NOT a hop — so this can exceed `hopCount()`, and on a
    // programme shorter than one hop it can be 1 while `hopCount()` is 0. Do not form a ratio of the two
    // without bounding it.
    std::int64_t damagedHops() const noexcept { return damagedHops_; }
    std::int64_t tailPartialSamples() const noexcept { return tailPartialSamples_; }
    double tailPartialEnergy() const noexcept { return tailPartialEnergy_; }

    std::int64_t eventCount() const noexcept { return eventCount_; }
    std::int64_t storedEventCount() const noexcept
        { return std::min<std::int64_t> (eventCount_, (std::int64_t) events_.size()); }
    bool eventsComplete() const noexcept { return eventCount_ <= (std::int64_t) events_.size(); }
    BandBurst event (std::int64_t i) const noexcept
        { return i >= 0 && i < storedEventCount() ? events_[(std::size_t) i] : BandBurst {}; }

    bool eventsValid() const noexcept { return eligibleHops_ > 0; }
    BandBurstsInvalid eventsInvalidReason() const noexcept
        { return eligibleHops_ > 0 ? BandBurstsInvalid::None : BandBurstsInvalid::NoEligibleHop; }

    // Per-channel band energy over the whole programme — the one thing a programme-wide sum loses. A
    // compensated sum, because this one runs to tens of millions of terms. It is the total over the
    // samples that WERE finite, which is a defined quantity and is still published when the flag below
    // is false: `ClipDetector` sets the precedent by publishing peak and DC over the finite samples and
    // counting the rest (`ClipDetector.h:261`). The flag says the programme had damage; it does not
    // blank a number that is well defined over what was measurable.
    double bandEnergy (int c) const noexcept
        { return has (c) ? chans_[(std::size_t) c].energy + chans_[(std::size_t) c].energyComp : 0.0; }
    std::int64_t nonFiniteSamples (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].nonFinite : 0; }
    std::int64_t absentSamples (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].absent : 0; }
    // CHANNEL-samples, not samples: one frame that overflows every one of 16 prepared channels counts 16
    // here. (`nonFiniteSamples(c)` beside it is per channel, which is the asymmetry to keep in mind.)
    std::int64_t overflowSamples() const noexcept { return overflowSamples_; }
    std::int64_t firstNonFiniteAt() const noexcept { return firstNonFiniteAt_; }   // −1 = none

    bool programmeEnergyValid() const noexcept { return firstNonFiniteAt_ < 0 && overflowSamples_ == 0; }
    BandBurstsInvalid programmeEnergyInvalidReason() const noexcept
        { return programmeEnergyValid() ? BandBurstsInvalid::None : BandBurstsInvalid::NonFiniteInput; }

    // --- periodicity: integer counts, no arithmetic, exactly invariant, and no verdict ---
    // The ADJACENT spacing histogram answers "what is the typical gap"; the AUTOCORRELATION answers "is
    // there a period", which is not the same question. A hi-hat with one hit missed turns a clean P into
    // P, 2P in the first (at 30 % dropout, 70 % of the histogram's mass stays at P — the (1−q)² = 49 % often
// quoted for this is the surviving fraction of the ORIGINAL pair count, a different denominator) while
// barely perturbing
    // a comb in the second, because removing one onset costs O(n) of O(n²) pairs. A perfectly periodic
    // train at a NON-INTEGER period puts mass around every multiple, so the spacing is readable from the
    // tenth peak even where bin-1 resolution is 20 % of the period. The adjacent histogram is the m = 1
    // slice of the autocorrelation and is otherwise reconstructable from the event list's `start` fields
    // — its own unique value is that it survives `maxEvents` exhaustion, which the autocorrelation does
    // too. Neither is a periodicity SCORE: with two or three onsets a mode carries 100 % of the mass and
    // means nothing, so `intervalCount()` and `onsetCount()` sit beside them and the consumer decides
    // what is enough. The useful range is about 5 to 512 hops (50 ms – 5.1 s at the default hop): below
    // ~5 the 10 ms grid quantises a period into two adjacent bins, above 512 it is in the overflow count.
    std::int64_t onsetCount() const noexcept { return onsetCount_; }
    std::int64_t intervalCount() const noexcept { return intervalCount_; }
    std::int64_t intervalOverflow() const noexcept { return ioiOverflow_; }        // spacings past kIoiBins
    std::int64_t intervalBin (int hops) const noexcept                             // spacing of `hops` hops
        { return hops >= 1 && hops <= kIoiBins ? ioi_[(std::size_t) (hops - 1)] : 0; }
    std::int64_t lagBin (int hops) const noexcept                                  // onset pairs `hops` apart
        { return hops >= 1 && hops <= kMaxLag ? lag_[(std::size_t) (hops - 1)] : 0; }

    // The most common adjacent spacing; ties go to the SHORTER one. 0 when there is no interval at all.
    int modalIntervalHops() const noexcept
    {
        int best = 0; std::int64_t bestN = 0;
        for (int b = 1; b <= kIoiBins; ++b)
        {
            const std::int64_t v = ioi_[(std::size_t) (b - 1)];
            if (v > bestN) { bestN = v; best = b; }               // strict >: the first (shortest) wins
        }
        return best;
    }
    // How many spacings sit within one hop of the mode — the jitter a real player has. 0 when no mode.
    std::int64_t modalIntervalMass() const noexcept
    {
        const int m = modalIntervalHops();
        if (m == 0) return 0;
        return intervalBin (m - 1) + intervalBin (m) + intervalBin (m + 1);
    }

using CrossoverType = eq::DeterministicCrossover2;   // asserted by the math-policy suite

private:
    struct Channel
    {
        double       energy     = 0.0;    // compensated sum of the squared band signal
        double       energyComp = 0.0;
        std::int64_t nonFinite  = 0;
        std::int64_t absent     = 0;
    };

    bool has (int c) const noexcept { return c >= 0 && c < (int) chans_.size(); }

    // Neumaier compensation. Two constant scalars and one fixed recurrence — legal under law 8a, and
    // `ClipDetector.h:393` already carries its DC sum this way.
    static void addCompensated (double& sum, double& comp, double v) noexcept
    {
        const double t = sum + v;
        comp += (std::fabs (sum) >= std::fabs (v)) ? ((sum - t) + v) : ((v - t) + sum);
        sum = t;
    }

    // The baseline: the upper median of the ring, SELECTED and not computed. No allocation — the standard
    // libraries implement `nth_element` as an in-place partition-based quickselect over the scratch this
    // object already owns (not "introselect": libc++'s has no introspective depth fallback, which is a
    // complexity note and not a correctness one). What matters for the report is that the VALUE left at
    // the midpoint is the k-th smallest whatever pivots an implementation chose, so the baseline is
    // portable even though the permutation around it is not. The ring never holds a non-finite value —
    // both gates are upstream of it — so the ordering is total and the selection is unambiguous.
    double baselineNow() noexcept
    {
        std::copy (ring_.begin(), ring_.end(), scratch_.begin());
        const auto mid = scratch_.begin() + (std::ptrdiff_t) (baselineHops_ / 2);
        std::nth_element (scratch_.begin(), mid, scratch_.end());
        return *mid;
    }

    bool onsetBit (std::int64_t hop) const noexcept
    {
        const std::size_t idx = (std::size_t) (hop % (std::int64_t) kMaxLag);
        return ((onsetBits_[idx / 64u] >> (idx % 64u)) & 1ull) != 0ull;
    }
    void writeOnsetBit (std::int64_t hop, bool v) noexcept
    {
        const std::size_t idx = (std::size_t) (hop % (std::int64_t) kMaxLag);
        const std::uint64_t m = 1ull << (idx % 64u);
        if (v) onsetBits_[idx / 64u] |= m;
        else   onsetBits_[idx / 64u] &= ~m;
    }

    // An onset at hop k: one adjacent spacing, and one pair count for every earlier onset within kMaxLag.
    // Both survive event-list exhaustion, because neither reads the list — the previous onset is a scalar
    // and the pair history is a bit ring (law 11: exhaustion is data, and the counters go on).
    void noteOnset (std::int64_t k) noexcept
    {
        ++onsetCount_;
        if (haveLastOnset_)
        {
            const std::int64_t d = k - lastOnsetHop_;
            ++intervalCount_;
            if (d >= 1 && d <= (std::int64_t) kIoiBins) ++ioi_[(std::size_t) (d - 1)];
            else                                        ++ioiOverflow_;
        }
        lastOnsetHop_ = k; haveLastOnset_ = true;

        for (int L = 1; L <= kMaxLag; ++L)
        {
            const std::int64_t prev = k - (std::int64_t) L;
            if (prev < 0) break;
            if (onsetBit (prev)) ++lag_[(std::size_t) (L - 1)];
        }
    }

    void closeEvent (std::int64_t endSample, bool byFinish) noexcept
    {
        BandBurst e;
        e.start  = evStart_;
        e.length = endSample - evStart_;
        e.peakAt = peakAt_;
        e.hops   = evHops_;
        e.energy = evEnergy_;
        e.touchedNonFinite         = evTouched_;
        e.baselineTouchedNonFinite = evBaseTouched_;
        e.closedByFinish           = byFinish;
        // Raw linear evidence, then the dB — ONCE, and from the ACCUMULATORS rather than from the two
        // published scalars. `10·(log10 a − log10 b)` would have been the natural spelling and is a trap:
        // it is not invariant under a global gain change (10·(log10 20 − log10 3) and
        // 10·(log10 80 − log10 12) differ in the last bit), while a ratio of exactly-scaled operands is.
        // peakDen_ > 0 holds for every event: a hop with a zero baseline is never eligible.
        const double hop = (double) hopSamples_;
        e.peakPower     = peakNum_ / hop;
        e.peakBaseline  = peakDen_ / hop;
        e.peakWidePower = peakWide_ / hop;
        e.peakExcessDb  = peakDen_ > 0.0 ? 10.0 * core::det::log10 (peakNum_ / peakDen_) : 0.0;
        if (eventCount_ < (std::int64_t) events_.size()) events_[(std::size_t) eventCount_] = e;
        ++eventCount_;              // law 11: exhaustion is DATA — the count and both histograms go on
        inEvent_ = false;
    }

    void closeHop (std::int64_t lengthSamples, bool full) noexcept
    {
        const std::int64_t k      = hopCount_;
        const std::int64_t startS = k * (std::int64_t) hopSamples_;
        const std::int64_t endS   = startS + lengthSamples;
        const double       energy = hopSumSq_;
        const double       wide   = hopWideSumSq_;
        const bool         dmg    = hopDamaged_;

        double baseline = 0.0;
        bool   eligible = false, above = false, opened = false;

        if (full)
        {
            // DECIDE against the ring, and only then write into it. Slot k % baselineHops still holds
            // hop k − baselineHops at this point, so the window walked is exactly [k − N, k − 1].
            if (k >= (std::int64_t) baselineHops_)
            {
                baseline = baselineNow();
                if (baseline > 0.0) { eligible = true; ++eligibleHops_; }
                else                { ++zeroBaselineHops_; }
            }

            if (eligible)
            {
                above = energy > baseline * (inEvent_ ? exitRatio_ : enterRatio_);

                if (! inEvent_ && above)
                {
                    inEvent_ = true; opened = true;
                    evStart_ = startS; evHops_ = 0; evEnergy_ = 0.0;
                    evTouched_ = false; evBaseTouched_ = false;
                    // SEED the peak from the opening hop. Leaving it at 0/0 would make the strict
                    // cross-multiplication below reject every candidate, and the event would report a
                    // peak of zero at sample zero.
                    peakNum_ = energy; peakDen_ = baseline; peakWide_ = wide; peakAt_ = startS;
                    noteOnset (k);
                }

                if (inEvent_)
                {
                    if (above)
                    {
                        ++evHops_; ++burstHops_; evEnergy_ += energy;
                        if (dmg) evTouched_ = true;
                        if (damagedInRing_ > 0) evBaseTouched_ = true;
                        // Greatest EXCESS, compared without dividing. Strict >, so a tie keeps the
                        // EARLIER hop. Both products are bounded well inside double: the largest legal
                        // energy is hopSamples·kMaxChannels·FLT_MAX² = 768000·16·1.158e77 = 1.423e84, so
                        // the product is 2.025e168 against a DBL_MAX of 1.8e308; and the smallest
                        // non-zero energy is one float32's square, 1.96e-90, so the ratio the dB is
                        // taken of cannot overflow either. Law 10's pin applies here too.
                        volatile double lhs, rhs;
                        lhs = energy * peakDen_;
                        rhs = peakNum_ * baseline;
                        if (lhs > rhs)
                            { peakNum_ = energy; peakDen_ = baseline; peakWide_ = wide; peakAt_ = startS; }
                    }
                    else { closeEvent (startS, false); }
                }
            }
            else if (inEvent_)
            {
                // An event cannot continue through a hop that has no baseline to continue against.
                closeEvent (startS, false);
            }
        }

        if (observer_ != nullptr)
        {
            BandBurstsHopTrace t;
            t.hopIndex = k; t.startSample = startS; t.endSample = endS;
            t.energy = energy; t.wideEnergy = wide; t.baseline = baseline;
            t.full = full; t.eligible = eligible; t.damaged = dmg; t.above = above;
            t.inEvent = inEvent_; t.eventStart = inEvent_ ? evStart_ : (std::int64_t) -1;
            observer_ (observerUser_, t);
        }

        if (dmg) ++damagedHops_;

        if (full)
        {
            writeOnsetBit (k, opened);

            // Into the ring, in absolute hop order. A burst hop and a damaged hop both go in: the ring
            // is the programme's own history, and leaving either out would make the window's span a
            // function of the programme instead of a constant.
            const std::size_t slot = (std::size_t) (k % (std::int64_t) baselineHops_);
            if (ringFlags_[slot] != 0u) --damagedInRing_;      // the hop this evicts (0 before the ring fills)
            ring_[slot]      = energy;
            ringFlags_[slot] = dmg ? (std::uint8_t) 1u : (std::uint8_t) 0u;
            if (dmg) ++damagedInRing_;

            ++hopCount_;
            hopSumSq_ = 0.0; hopWideSumSq_ = 0.0; hopDamaged_ = false;
        }
    }

    BandBurstsParams params_ {};
    BandBurstsHopObserver observer_ = nullptr;
    void* observerUser_ = nullptr;

    CrossoverType xLow_, xHigh_;   // declared THROUGH the public alias, so the assertion cannot drift from the member   // deterministic coefficients (see core::DetMath)
    core::StateGrid  grid_;

    double sampleRate_   = 48000.0;
    int    channels_     = 0;
    int    hopSamples_   = 0;
    int    baselineHops_ = 0;
    double enterRatio_   = 0.0, exitRatio_ = 0.0;
    bool   prepared_     = false, finished_ = false;

    std::vector<double>       ring_, scratch_;
    std::vector<std::uint8_t> ringFlags_;
    std::vector<BandBurst>    events_;
    std::vector<Channel>      chans_;

    std::int64_t totalSamples_ = 0, nextHopEnd_ = 0, hopCount_ = 0;
    double       hopSumSq_ = 0.0, hopWideSumSq_ = 0.0;
    bool         hopDamaged_ = false;

    std::int64_t eligibleHops_ = 0, zeroBaselineHops_ = 0, burstHops_ = 0, damagedHops_ = 0;
    int          damagedInRing_ = 0;
    std::int64_t overflowSamples_ = 0, firstNonFiniteAt_ = -1;
    std::int64_t tailPartialSamples_ = 0;
    double       tailPartialEnergy_ = 0.0;

    bool         inEvent_ = false;
    std::int64_t evStart_ = 0, evHops_ = 0, peakAt_ = 0, eventCount_ = 0;
    double       evEnergy_ = 0.0, peakNum_ = 0.0, peakDen_ = 0.0, peakWide_ = 0.0;
    bool         evTouched_ = false, evBaseTouched_ = false;

    std::uint64_t onsetBits_[(std::size_t) kMaxLag / 64u] {};
    std::int64_t  ioi_[(std::size_t) kIoiBins] {};
    std::int64_t  lag_[(std::size_t) kMaxLag] {};
    std::int64_t  onsetCount_ = 0, intervalCount_ = 0, ioiOverflow_ = 0, lastOnsetHop_ = 0;
    bool          haveLastOnset_ = false;
};

} // namespace felitronics::analysis
