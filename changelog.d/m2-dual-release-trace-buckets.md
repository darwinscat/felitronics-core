<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### limiter · mastering · tools — the limiter's dual release, the caller's gain-reduction trace size; C ABI v6

**`limiter::TruePeakLimiter` — a dual release.** `TruePeakLimiterParams::dualRelease`, off by default, and
`slowReleaseMs` (200 ms). On, a second envelope — instant attack, exponential release over `slowReleaseMs` — follows the
least reduction the gain law required over the last 50 ms (`TruePeakLimiter::kSlowWindowMs`, `slowWindowSamplesFor`),
and the reduction applied is the larger of the two envelopes'; `releaseMs` is then the fast one's. Both envelopes read
only the required reduction, so the loudness search's scale law holds under it (5.754e-07 on its 3x3 grid) and the
reduction stays monotone in drive (a 20-step drive sweep gives back no reduction in any tap sample). Off, the output is
the single release bit for bit, whatever `slowReleaseMs` holds. Switching it on starts the slow envelope from rest;
switching it off continues the fast envelope from the reduction applied at the switch.
- The slow envelope binds only where the required reduction stays below 0 dB for the whole window and is followed by
  less than the fast envelope keeps. At the 1 ms lookahead a programme whose peak falls under the ceiling within any
  50 ms — bass zero crossings do — never engages it: its render is the fast release's, bit for bit.
- `effectiveSlowReleaseMs()` reads the floored slow release back, 0 while off; `MasteringChainResolved::limiterSlowReleaseMs`
  carries it. The window is prepared whether or not the dual release is on: `Storage::slowWindow`, 115 200 B at 48 kHz
  and 4x.
- The true-peak characterisation under it is in the header: the same worst rows as the single release on the same
  witnesses, and a lower excess where the slow envelope works (a held plateau at a 1 ms fast release, 2x: +0.0786 dB
  against +0.2206).

**`mastering::TargetLoudnessSolver` — the trace's size is the caller's.** `LoudnessRequest::grTraceBuckets`, 1000 by
default, 1..65536, anything else `InvalidRequest` before any pass. Each trace holds min(grTraceBuckets, frames) buckets,
on the heap (`GainReductionTrace::bucket` is a vector; `kDefaultBuckets` 1000, `kMaxBuckets` 65536), and a bucket's
`samples` and `nonFinite` are 64-bit: one bucket fed 2^32 + 5 samples counts them. `solveBytes (fs, nch, frames,
grTraceBuckets)` and `DeliveredMastering::solveBytes (…, grTraceBuckets)` count the two traces — 64 000 B at the default
— and a solution record no longer carries them by value. At the default the trace is the one before, bit for bit.

**C ABI v6.** `fc_master_params` gains `limiterDualRelease`, `_pad0` and `limiterSlowReleaseMs` (6584 B),
`fc_master_resolved` gains `limiterSlowReleaseMs` (96 B), `fc_loudness_request` gains `grTraceBuckets` and `_pad0`
(128 B); at their defaults a call is v5's. `fc_solution_gr_trace64` copies `fc_gr_trace_bucket64`, the 64-bit counts;
`fc_solution_gr_trace` refuses a bucket whose count does not fit 32 bits with `FC_ERR_RANGE`, writing nothing.
`fc_master_need_solve (h, req, frames, out)` budgets a solve for a request; `FC_NEED_SOLVE` through `fc_master_need`
budgets the default 1000 buckets, and now counts the traces too. `fc-master-layout.mjs` is v6; `fcore_master` takes
`lim.dual=`, `lim.slowRelease=` and `grTraceBuckets=`.
