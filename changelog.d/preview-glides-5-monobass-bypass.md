### stereo — MonoBass::setBypass: a host bypass of the whole island that rides its own fades

`setBypass (true)` fades the bass's crossfade to dry and the air's plateau to 0 dB (`kSmoothingMs` each, the fades
`enabled` uses since the corners-glide note), after which the island retires on the sample clock — a bit-exact
passthrough that returns before touching the buffer — and `setBypass (false)` fades both back in from filters
restarted at zero. It OVERRIDES the parameters without replacing them: `params()` and `air()` keep reporting the
caller's settings, and a write made while bypassed cannot switch either tool back on (the bypass is part of both fade
targets). Before the stream's first sample it lands at once, so a stage bypassed from the start is untouched input.
felitronics-mastering-core's MasteringChain uses it for `bypassMonoBass`, which used to skip the stage with a reset on
both edges: -21.5 dBFS max|Δ²y| into the bypass and -48.7 out of it on a wide 103.7 Hz tone, -74.3 and -75.5 now. Pinned
in `felitronics_monobass_glide_tests`.
