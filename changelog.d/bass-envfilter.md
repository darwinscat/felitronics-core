### envfilter — a new module: `EnvelopeFilter`, the bass amp's envelope filter (auto-wah), and its one-slider law

`felitronics::envfilter` — header-only, JUCE-free, no heap; depends on `core`, `dynamics` and `eq`. A peak
`dynamics::EnvelopeFollower` on the loudest channel — or on an external key, `process (io, nch, n, key, nk)`, the
dynamics house shape — places the cutoff of a resonant Cytomic TPT state-variable filter on a log-frequency sweep
(150 Hz … 3 kHz by default): fully open when envelope + sensitivity reaches 0 dBFS, closed 30 dB below; UP rests
low and a hit opens it, DOWN rests high and a hit pushes it down. The playing reaches the cutoff through a 4 ms
smoother, the hand's edits (direction, range, sensitivity) glide 20 ms first, and the filter is redesigned every
sample. Low pass or band pass (a 10 ms crossfade), both scaled by 1/√max(1, Q) (Q 0.5 … 12): the low pass's
passband sits at −6.0 dB at Q 4 and −10.8 dB at Q 12. `keepLows` takes the band under `keepLowsHz` (150 Hz) out
of the filter and adds it back clean, with a split per mode that cannot cancel it: an exact complement in low pass
(dry = the input), an in-phase LR4 pair in band pass (dry = its allpass). The doors never crossfade the two phase
families: the low pass has a 10 ms door from the input, the band pass a 70 ms ladder (a crossfade to the LR4
allpass parked at 0.45 fs, a 50 ms sweep of it down to the corner, a crossfade to the band pass), so a wake, a
disengage and a mode change keep a tone at or above the corner within 0.92 dB of its settled level. Everything
glides and lands; `mix` 0 closes the doors (counted), after which the stage is idle: the buffer is returned
untouched and nothing is computed; a wake is a restart. Mono or stereo (linked), 8 … 768 kHz, zero latency, law 8/8a on `core::StateGrid`,
law 11 (an external key is consumed through a gap). About 1.04 % of real time at 48 kHz stereo with the sweep moving and the glides at rest (deb, i9-13900H) — not a
worst case: both modes always run, and a corner edit or the ladder redesigns the splits or the allpass per sample.

`EnvelopeFilterLaw.h` — `paramsForAmount (amount, menu)`: mix rises 0 → 1 over the slider's lower half,
sensitivity −6 → +18 dB over the whole travel; `mix` and `sensitivityDb` are law-owned, everything else is the
menu's; a non-finite amount returns `menu` unchanged.

`core` — DetMath.h's cost note names the two stages that put `core::det` on a per-sample path deliberately
(`delay::TapeDelay`, `envfilter::EnvelopeFilter`), with the measured price.

Measured (`felitronics_envfilter_tests`, 44.1 / 48 / 96 kHz): 12 dB more input moves the cutoff 1.729 oct, the
law's number, and a probe tone's gain under each burst matches the published cutoff's analytic response within
0.01 dB; the kernel at Q 12 with its cutoff flipped between 20 Hz and 0.45 fs every sample stays finite (largest
tap 22.6 on ±1 noise); a swept output holds −150 dB of energy above 6 kHz (a cutoff held for 32 samples: −105);
over 360 cases of cutoff × Q × corner × mix × tone the band under the corner never drops below its floor (0 dB in
low pass, the LR4 low share in band pass) — the band pass parked at 104 Hz over a 100 Hz corner that took a 98 Hz
note down 13.6 dB now gives it +2.4 dB; a corner retuned under a sustained note, a direction flip and a range
edit move without a step (≤ 1.04× the steady step; the flip takes 51 ms to cover 90 % of its travel).
