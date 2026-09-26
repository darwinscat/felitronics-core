### limiter — a lower ceiling glides in over 2 ms; the bound holds against the ceiling in force

A lower `ceilingDbTp` acted on the next oversampled sample — the attack is instant by design — so a ceiling pulled
down while the limiter held a tone stepped the output: through the mastering chain on a -12 dBFS 227 Hz sine, 0 ->
-18 dBTP put max|Δ²y| at -29.9 dBFS against -73.1 for the steady tone. The ceiling the detector compares against now
falls to a lower target along a one-pole of `TruePeakLimiter::kCeilingGlideMs` = 2 ms per oversampled sample (the
product and the sum in separate statements, landing exactly once the residual is under 1e-9 dB): -64.4. 2 ms was
measured against 1 ms (-52.1) and 5 ms (-62.8). A HIGHER ceiling still lands at once — the release already makes that
move smooth, and gliding it too would only delay recovery — and every write before the first sample after
prepare()/reset() lands at once, so a limiter configured before its first sample renders what it always did.

THE BOUND STILL HOLDS, against the ceiling IN FORCE: every emitted oversampled sample is inside its own detector
window, so its gain is at most 10^((c_i - smaxDb_i)/20) with smax_i at least its own magnitude — |out_i| <= 10^(c_i/20)
whatever path c takes; the dual release takes the min of two reductions each at least that deep, and the peak clipper
acts before the window sees the sample and now rides the gliding ceiling too (`clipAt()`), so it never cuts at a level
the ceiling has not reached. During a downward glide c_i lies between the old ceiling and the new one: the limit is
lowered over a couple of milliseconds instead of in one sample, and never exceeded. `effectiveCeilingDbTp()` and
`clipThresholdDbTp()` keep reporting the target; `ceilingNowDbTp()` reads the one in force. New
`felitronics_limiter_ceiling_glide_tests`: the bound replayed per oversampled sample through the taps against the
recursion's c_i (clipper off and on, single and dual release), the landing, the instant rise, the snap, law 8a under
per-sample, 64 and ragged cuts, and the edge under -60 dBFS.

THE REVIEW ROUND (codex astra): a clock-only call froze the glide — a lower ceiling written during a gap stood at the old
one until audio returned; the gap now spends it as audio time, one step per oversampled sample, the same steps the audio
loop takes (a gap cut in pieces lands where one call does). And "never exceeded" carries the limiter's own float
rounding — the gain is a float of a float dB, so an emitted sample can sit an ulp or two above 10^(c_i/20) (1e-6 dB
measured, contraction on and off alike), the margin the static ceiling always had; the doc says so and the test's
tolerance is 1e-5 dB, a few ulps, where it was 1e-4.

