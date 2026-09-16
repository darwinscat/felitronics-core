<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### convolution — an IR whose header claims a rate under 8000 Hz loads as is, like any other broken rate

**What `convolution::CabConvolver::loadIR` now promises about the IR's rate:** a KNOWN rate is a finite rate of
8000 Hz (`core::kMinSampleRate`) or more, and only a known rate off the host's is resampled. Everything else —
NaN, zero, negative, ±inf, and now any rate under 8000 Hz — is UNKNOWN, and an unknown rate loads the taps as
they are, with no rate factor (P67's rule: the samples are fine, only the metadata is broken, and refusing
would play silence). Above, nothing changed.

**Why:** a header that says 44.1 is kilohertz written as hertz, and it was trusted: a 4096-tap cabinet on a
48 kHz host was resampled x1088 into 4 458 231 taps — measured 3.76 s of `loadIR`, now 0.2 ms — and on the
verbatim (reverb) path scaled by the rate factor as well. Now it is the NaN load, bit for bit — the staged taps, `irNormalizationGain()`,
`irNormalizationGainDb()` and what the convolver plays, on both paths.

**What that moves, and what it does not:** the rate factor on the verbatim path now starts at 8000 / 3e6 =
2.67e-3 (−51.48 dB, an 8 kHz IR on a 3 MHz host) instead of about 4.7e-10. A 4096-tap cabinet at any known
rate costs at most x375 of itself now (a 3 MHz host); on a 48 kHz host it is x6 and 25 ms at 8001 Hz. What a
known rate can still ask for is unchanged and is bounded by the resampler's 2^24-tap output, not by the rate:
an 8001 Hz IR of 2.8 million taps (a 350-second file) took 17.4 s of `loadIR` on a 48 kHz host, measured.
