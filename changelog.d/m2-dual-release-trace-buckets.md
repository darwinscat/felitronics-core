<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### limiter · mastering · tools — the limiter's dual release, the caller's gain-reduction trace size; C ABI v6

**`limiter::TruePeakLimiter` — dual release.** `TruePeakLimiterParams::dualRelease` (off) and `slowReleaseMs` (200 ms).
On, a second envelope with instant attack and release `slowReleaseMs` takes as input the least, over the last 70 ms
(`kSlowWindowMs`), of the largest reduction required within 20 ms (`kSlowBridgeMs`); the applied reduction is the larger
of the two envelopes', and `releaseMs` is the fast one's. Both envelopes read only the required reduction, so the
loudness search's scale law and the reduction's monotonicity in drive hold. Off, the output is the single release bit
for bit, whatever `slowReleaseMs` holds. Switching on starts the slow envelope and its windows from 0 dB; switching off
continues the fast envelope from the applied reduction. `effectiveSlowReleaseMs()` (0 while off) and
`MasteringChainResolved::limiterSlowReleaseMs`. Both windows are prepared whether or not it is on
(`Storage::slowWindow`, `slowBridge`).

**`mastering::TargetLoudnessSolver` — the trace's size is the caller's.** `LoudnessRequest::grTraceBuckets`: 1000 by
default, 1..65536, otherwise `InvalidRequest` before any pass. Each trace holds min(grTraceBuckets, frames) buckets on
the heap (`GainReductionTrace::bucket`, `kDefaultBuckets`, `kMaxBuckets`); a bucket's `samples` and `nonFinite` are
64-bit. `solveBytes (fs, nch, frames, grTraceBuckets)` and `DeliveredMastering::solveBytes (…, grTraceBuckets)` count
both traces. At 1000 buckets the trace is the one before, bit for bit.

**C ABI v6.** `fc_master_params` gains `limiterDualRelease`, `_pad0`, `limiterSlowReleaseMs` (6584 B);
`fc_master_resolved` gains `limiterSlowReleaseMs` (96 B); `fc_loudness_request` gains `grTraceBuckets`, `_pad0`
(128 B). `fc_solution_gr_trace64` copies `fc_gr_trace_bucket64`; `fc_solution_gr_trace` answers `FC_ERR_RANGE` for a
count past 32 bits. `fc_master_need_solve (h, req, frames, out)` budgets a solve for a request; `FC_NEED_SOLVE` budgets
1000 buckets. `fc-master-layout.mjs` is v6; `fcore_master` takes `lim.dual=`, `lim.slowRelease=`, `grTraceBuckets=`.
