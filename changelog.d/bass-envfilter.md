### envfilter — a new module: `EnvelopeFilter`, the bass amp's envelope filter (auto-wah), and its one-slider law

`felitronics::envfilter` — header-only, JUCE-free, no heap; depends on `core`, `dynamics` and `eq`. A peak
`dynamics::EnvelopeFollower` on the loudest channel (attack 3 ms, release 120 ms by default) places the cutoff
of a resonant Cytomic TPT state-variable filter on a log-frequency sweep (150 Hz … 3 kHz by default): fully
open when envelope + sensitivity reaches 0 dBFS, closed 30 dB below; UP rests low and a hit opens it, DOWN rests
high and a hit pushes it down. The cutoff passes one 4 ms smoother and the filter is redesigned every sample.
Low-pass and band-pass come from one state and crossfade over 10 ms; both are scaled by 1/√max(1, Q) (Q 0.5 … 12).
`keepLows` takes everything under `keepLowsHz` (100 Hz) out of the filter and adds it back clean; the lows are
the exact complement `x − HP(x)`, so the dry/wet mix never notches the bass. mix, sensitivity, resonance and
keepLows glide (20 ms) and land; `mix` 0 puts the stage to sleep: the buffer is returned untouched and nothing is
computed. Mono or stereo (linked), 8 … 768 kHz, zero latency, law 8/8a on `core::StateGrid`, law 11.

`EnvelopeFilterLaw.h` — `paramsForAmount (amount, menu)`: mix rises 0 → 1 over the slider's lower half,
sensitivity −6 → +18 dB over the whole travel; everything else is the menu's.

Measured (`felitronics_envfilter_tests`, 44.1 / 48 / 96 kHz): 12 dB more input moves the cutoff 1.729 oct, the
law's number, and a probe tone's gain under each burst matches the published cutoff's analytic response within
0.01 dB; the kernel at Q 12 with its cutoff flipped between 20 Hz and 0.45 fs every sample stays finite (largest
tap 22.6 on ±1 noise); a 0 dBFS slap at the fastest settings peaks at −0.55 dBFS; a swept output holds −150 dB
of energy above 6 kHz (a cutoff held for 32 samples: −105); keepLows 1 holds 31 … 73 Hz within 1.23 dB of the
input, where 0 loses 43 dB in band-pass.
