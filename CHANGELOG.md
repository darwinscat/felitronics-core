<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Changelog

Notable changes to felitronics-core. Releases are git tags (`vX.Y.Z`); the project VERSION lives in
`CMakeLists.txt`.

## Unreleased

- **BREAKING (behaviour), `dynamics`, `deesser`, `dynamiceq`, `poweramp`, `multiband`, `core`: LAW 11c —
  A PAUSE IS SILENCE.** A call with `nch == 0, n > 0` now advances a stage's SHARED, one-per-instance
  ballistics exactly as `n` samples of digital silence at a live width would, instead of freezing them.
  Law 11(d) already called such a call a GAP IN THE STREAM rather than a no-op — audio time PASSED — and
  a detector standing still through passing time contradicts the sentence that made the grid advance.
  Per-channel memory keeps being dropped exactly as before (11a/11d are untouched), and `reset()` is not
  the answer either: it claims a stream RESTART where the caller said a gap.
  - **What the freeze cost, measured.** `dynamics::Compressor` held **-25.311 dB of gain reduction
    through a full second of gap** where the same second of silence releases to -2.076, and dipped the
    return by **-22.11 dB** on material below its threshold (-24.09 dB through ten seconds);
    `dynamiceq::DynamicEqBand` held -21.774 against -2.985; `deesser::DeEsser` -8.000 against -1.104;
    `dynamics::TransientShaper` moved the return by 1.76-2.60 dB. The loudest is `dynamics::NoiseGate`,
    and it is not a shifted envelope but a state machine that never fired: a gate that has to CLOSE
    through a pause stayed wide open, and a -54 dBFS tone on the return came out at -54 where silence
    gates it to -144 — **89.99 dB, 100 % of the construction ceiling.** All of those are now 0.00 dB.
  - **Seven addresses, chosen by MECHANISM.** The five above, plus `dynamiceq::LaneDynamics` (which
    deliberately DISENGAGED — the "dynamics were switched off" verb, not the "time passed" one; at width
    zero its STEREO lane now runs the control loop on silence, while L/R/M/S take the same "this lane
    stopped" branch they take at width one, because `laneRuns()` gates them on `nc == 2`) and
    `poweramp::PowerAmpStage`, whose one shared sag supply
    and thirteen block-rate glides stopped dead on a gap. `multiband::MultibandProcessor` forwards a gap
    to every band exactly once — it used to forward it TWICE to a bypassed band, invisible under freeze
    and a double clock under this law.
  - **`limiter::TruePeakLimiter` is deliberately NOT included.** It does a full `reset()` on a gap, which
    is wrong under any answer (gain reduction -5.08 dB to 0.00, a 48-sample hole on the return), but its
    ballistics are a sliding lookahead window rather than an exponential, and its per-channel state is
    audio that has not been emitted yet. That is its own item (P29) with its own rule.
  - **Bit-exactness, with its area.** Every call with `nch > 0` is bit-identical to the base — verified by
    hash over all seven stages, every mode x detector x link x control period, with width sweeps,
    self-keyed and externally keyed, and the gain-reduction tap. Two exceptions, each with its own test:
    `NoiseGate` now drops a stopped lane's sidechain high-pass (law 11a, which it never had — worth
    **89.99 dB** on a narrowing as well as on a gap), and `MultibandProcessor` no longer hands a bypassed
    band a row of NULL planes on a narrowing call, **which was a segfault, not a wrong number** — and the
    same fix has a second, non-crashing half: a bypassed band whose planes were already valid used to be
    clocked on the PREVIOUS chunk's split audio at a narrowing edge and is now clocked on digital silence
    (3.12 dB on its meter, 1.6 dB out of silence on un-bypass). The segfault itself needs the band
    bypassed BEFORE its first call at that width — once it has run live, the planes are filled.
  - **Two contract changes at width zero**, both consequences of "a pause is the same call carrying
    silence" and both tested: a `GainReductionTap` handed to `Compressor::process` is now FILLED on a
    zero-width call where it used to be left untouched, and an external key handed to the same call is now
    DEREFERENCED where it previously read nothing — so a key passed at width zero must be valid for `n`
    samples, exactly as at any other width.
  - **Cost.** The silent recurrence is autonomous, so it reaches a bitwise fixed point and the rest of the
    pause is free; and once the detector level reaches `core::kGainToDbFloor` the curve's output is a
    constant, so the per-sample work collapses to one multiply-add. Three of the seven have no such floor
    and run their full body until they park: one zero-width call covering a MINUTE costs `Compressor`
    1.49 ms, `NoiseGate` 0.20, `TransientShaper` 3.65, `DeEsser` 4.08 and `DynamicEqBand` 6.43 — against a
    2.67 ms callback budget. The same minute delivered as 128-sample calls costs **0.00304 ms** in its
    worst call, so an RT caller is unaffected and an offline caller that hands a whole transport jump as
    one call pays it once. `pow(c, n)` is deliberately not used —
    it is a different number from `n` rounded multiplications. The horizon is a property of the TIME
    CONSTANT, not of the pause: 23 609 samples at a 5 ms release, 457 808 at 100 ms, 4 461 677 at 1 s, and
    not reached in 200 000 000 at the coefficient cap.
  - **New API.** `core::kGainToDbFloor` and `core::sameBits`; `EnvelopeFollower::stateWord`/
    `advanceSilence`, `GainReductionFollower::advanceConstant`, `GainReductionPath::advanceSilence`,
    `LinkedDetector::stateWord`/`advanceSilence`, `poweramp::PowerAmpStage::sagDroop`. Additive only.
  - `dynamiceq::LaneDynamics`' park counter is now a saturating `long long`: as a `long` it is 32 bits on
    the MSVC row, where `parked += n` at 48 kHz is signed overflow after 12.4 hours — unreachable while a
    gap disengaged a lane in one step, reachable the moment a pause became something a lane spends.
- **BREAKING (reported latency), `nam`:** **`NamStage::latencySamples()` was 2.16 samples too long at
  44.1 kHz and up to 3.3 at 88.2 kHz.** It reported `ceil(3·hostSR/modelRunSR) + 3` — a guess at "~3
  samples of lookahead per stage" — where the geometry is exact and one line away in the same
  repository: `StreamResampler.h` says the identity ratio "passes the signal with a clean 2-sample
  delay", and that is true at EVERY ratio, because `reset()` leaves 3 leading history zeros with
  `pos = 1.0`, so output *k* reads input position *k·inPerOut − 2*. The round trip is therefore
  `2 + 2·hostSR/modelRunSR` host samples, now reported rounded to nearest: **6 → 4 at 44.1 kHz, 9 → 6
  at 96 kHz, 9 → 6 at 88.2 kHz, 5 → 3 at 22.05 kHz**; unchanged (0) at the model's own rate, where the
  resampler is not in the path at all. Hosts using the reported number for delay compensation move by
  that much — and **twice that** in the two shipped hosts, which sum a preamp and a poweramp stage
  (4 samples at 44.1 kHz, 6 at 96). **Nothing inside `rigplayer` moves**: slot alignment runs on
  `AlignmentTable::delayOf()` → `blendDelay()`/`lagTail_`, and none of those reads `latencySamples()`
  at all — `RigPlayer` only republishes the max of the two slots outward, and both changed identically.
  **No audio sample changes inside this repository. Downstream, audio does move, and it moves into
  alignment:** OrbitCab delays its dry/bypass path by this same number
  (`src/poweramp/PowerAmpRouter.cpp`, `src/core/CabEngine.cpp`) and orbit-amp does the same at the dry
  end of its crossfade, so the wet path sat at the true 3.84 samples while the dry was held at the
  reported 6 — a 2.16-sample mismatch whose first comb notch fell at ~10.2 kHz during an on↔off
  crossfade (3.0 samples and ~16 kHz at 96 kHz). It is now 0.16 and 0.00. ⚠️ **Three OrbitCab tests pin a number that is now known to be WRONG**
  (`tests/PowerAmpRouterAlignTests.cpp`, three `expectEquals(L, ceil(3·sr/48000) + 3)`). They will fail
  on the next core bump, and the fix is to replace the pinned value with the geometry
  `2 + 2·hostSR/modelRunSR` — **not** to restore the old formula in core. The old tests pinned the FORMULA, which is why
  nothing caught it; the new one measures the delay from the carrier phase of the shipped round trip
  (3.8375 / 6.0000 / 5.6750 samples, matching the geometry to four decimals) and asserts the reported
  integer is the nearest one to it.
- **`core`, docs:** **`StreamResampler`'s header claimed transparency it does not have, and now carries
  the measurement instead.** The old justification — *"the driven nonlinear stage masks the
  interpolation images"* — had no number behind it, and the quantity that had since been measured was a
  different one. Measured (`docs/STREAM-RESAMPLER-COST.md`, new): a phase-dependent kernel is a linear
  periodically time-varying filter whose per-phase gain has Fourier coefficients `H(Ω+2πk)`, so the
  "amplitude modulation" and the "interpolation images" are **one mechanism**, not two. The NAM round
  trip at 44.1 ↔ 48 kHz costs **−4.17 dB coherent and −9.27 dB worst-phase at 17.64 kHz** (−2.59/−5.14
  at 15 kHz, −5.48/−14.79 at 20 kHz), from a composite period of exactly 147 output samples. **That is
  ONE round trip; OrbitCab runs two `NamStage`s IN SERIES** (preamp → EQ → poweramp, and its own
  `updateLatency()` comment says so), i.e. four of these stages — measured **−9.03/−13.16 at 17.64 kHz
  and −12.09/−17.85 at 20 kHz**, where doubling the decibels would say −8.35/−18.53. On that chain the
  BEST phase falls from −0.61 dB to −5.20, so the top octave is down at every phase rather than only at
  some. `rigplayer` runs its two stages in parallel and stays on the one-round-trip row. That set
  is complete for the shipped priming but is a LINE through the two stages' phase torus, not the full
  product: over all 160 integer alignments the worst phase barely moves (−9.29 against −9.27) while the
  coherent carrier spans −3.59…−6.83, because it is an interference term between the stages. **The decimating direction has no stopband at
  all**: at phase *t = 0* the weights are `(0,1,0,0)`, a bare sample pick, so a tone above the output
  Nyquist survives at −3 dB rms / 0.0 dB SAMPLE peak (a time-domain fact with a kernel reason: at
  `t = 0` the weights are `(0,1,0,0)` so |M(0)| = 1 at every frequency, and the phase grid has points
  within 1/147 of zero where |M| is still −0.002 dB; across the measured rows no spectral line exceeds
  −4.67 dB) and folds back as TWO strong components — `44100 − g` at
  about −5 dB and `g − 3900` at about −7 dB — i.e. across **18.15–22.05 kHz**, not one top slice. Against the model's own
  aliasing floor the OUTPUT leg alone sits 3–24 dB below it on a high-gain capture but **up to +9.8 dB
  above it on a clean one**, at every level from 17.5 kHz up (the whole rate-match: above in 22 of 30 tone × level cells, up to
  +12.3 dB), and **driving harder does not help** — across a 42 dB sweep the error-to-signal ratio is
  flat in 16–22 kHz where the artifacts live and grows +10…+14 dB in 0–4 kHz where they do not, because
  the nonlinearity **demodulates** the input leg's images into the audible range: a 20 kHz tone at
  −18 dBFS into a high-gain capture returns a **100 Hz line at −17.7 dBFS, 14.5 dB louder than its own
  carrier**, against −174.6 dBFS through an ideal round trip. **No kernel change here**: the header now states the cost, the
  new `felitronics_core_streamresampler_lptv_tests` (75 checks) pins the table, the period-147 closure,
  the 0 dB decimation peak and the criterion itself — the round trip adds **−8.84 dBc at 17.5 kHz**, and
  a `tanh` has to be driven to **`tanh(6.2x)`** (bisected) before its own folding reaches that, so below
  a near-square-wave drive the rate-match is the LOUDER artifact. The candidate comparison is in the
  document for the product decision — including the part the two tone axes get wrong on their own: a
  32-tap sinc flattens both axes and is still **5 dB worse in the bass** on real DI through a driven
  capture, because its band edge feeds the same demodulation from a different cause. A 64-tap one at a
  0.99 cutoff is better in every band at every drive, for **+2.5 %** of what the stage already spends on
  the model and **61.4 host samples** of delay against today's 3.84.

- **BREAKING (behaviour + latency), `oversampling`, `saturation`, `limiter`, `poweramp`:** **the shipped
  `tapsPerPhase` default rises from 32 to 64, and the reason is aliasing, not the pass band.**
  `PolyphaseOversampler`'s cutoff is FIXED at 0.90 × baseband Nyquist, so its transition band has to fit
  between 0.45 fs and the fold at 0.50 fs, and `tapsPerPhase` is the only thing that decides whether it
  does. The design DECLARES its own target one line from the taps — `beta = 9.0`, a ~90 dB Kaiser
  stopband — and **at 32 taps it delivered 27 dB.** Everything above 0.50 fs folds straight back into the
  audio band, so that is not a nicety. Measured end to end on a `Saturator` at 0.17 fs and 4×, TOTAL
  non-harmonic energy — every bin that is not a harmonic BELOW Nyquist, so a harmonic that folded to get
  there counts as the aliasing it is:

    drive     +0       +6      +12      +18      +24      +30      +36  dB
    32     −135.3    −60.1    −48.2    −44.3    −31.8    −25.3    −22.9  dBc
    64     −130.8   −132.3   −101.0    −50.3    −31.3    −24.8    −22.4  dBc

  **At the drive this stage is used at — `Saturator::Params` calls 1..6 dB the mastering range — the
  taps remove 50 to 72 dB of aliasing.** They stop helping above about +24 dB, where the harmonic series
  reaches past the OS Nyquist and folds INSIDE the oversampled domain, which no decimation filter can
  reach: that is the oversampling FACTOR's axis, and there the total is 0.5 dB worse, because a flatter
  pass band also delivers what had already folded. The component the taps own at every drive is the
  transition-band leakage — the 3rd harmonic of a 0.17 fs tone lands at 0.51 fs and folds to 0.49 fs,
  where 32 taps left it at −45.5 dBc and 64 puts it at −139.4, into the float noise. (That last pair is
  one tone landing in one of the kernel's nulls; the honest figure for what 64 taps give ANYWHERE in the
  fold region is the worst case below, −90.5 dB, not the null.) Worst rejection over
  the whole fold region, worst of factors 2/4/8: 32 → −26.9 dB, 48 → −51.1, 56 → −76.2, 57 → −82.6,
  **58 → −90.7**, 59 → −89.8, 60 → −90.7, **64 → −90.5**, 96 → −94.3. Read as a knee, not a step:
  below ~58 the transition is genuinely unfinished, at ~58 it reaches the Kaiser window's floor and then
  RIPPLES there by about a dB, so 59 and 61 fall a tenth of a dB short while 58, 60 and 64 clear it.
  "The first taps count meeting −90.0" is therefore 58, and that integer is an artefact of a hard bar on
  a rippling quantity. 64 is the round number past the knee; above it further taps buy pass-band width
  rather than rejection — which is why the answer is not 96.
  - **The pass band is the COROLLARY, and it is the half that was already written down.** Two oversampled
    stages in series — a clipper in front of a limiter, the real assembly — cost **−1.549 dB at 17.6 kHz
    and −6.033 at 18.5 kHz at 44.1 kHz** on the old default, i.e. every consumer that built that chain got
    a ~19 kHz lowpass silently. They now cost **+0.000 and −0.610**. `felitronics::mastering` has passed 64
    explicitly since it was written and is **bit-identical** across this change (verified over 942 912
    float32 values, 18 configurations; the same stand shows the DEFAULT paths differing, so it is not a
    blind null). A caller passing `tapsPerPhase` explicitly is likewise bit-identical — with one
    exception, which is the new guard below rather than the default: a topology that used to be accepted
    and is now refused (`factor` above 64, `tapsPerPhase` above 1024) does not "produce the same bits",
    it produces a `false` from `prepare()`.
  - **LATENCY MOVES, and it is host-visible.** Every affected stage reports `tapsPerPhase − 1`, so
    **31 → 63** samples; `TruePeakLimiter` at its default 1 ms lookahead and 48 kHz goes **79 → 111**.
    Read it from `latencySamples()`, which is what `mastering::MasteringChain` already does.
  - **CPU ROUGHLY DOUBLES in the FIR** — "no CPU cost" would be false here. Measured at 48 kHz, stereo,
    4×, block 512: `Saturator` 0.97 → 2.19 %RT, `TruePeakLimiter` 1.10 → 2.22 %RT, the pair in series
    **2.13 → 4.44 %RT** (2.09×).
  - **`poweramp::PowerAmpStage` gets the knob it never had** — `prepare (sampleRate, maxBlock,
    oversampleFactor, tapsPerPhase)` — and takes the same default. It used to hardcode 32 with no way for
    a caller to say otherwise, and its own aliasing gate certified that as adequate **because the gate's
    analysis window stopped at 10 kHz**, below where the transition-band leakage lands. Over the whole
    band the same gate reads **−68.7 dBc at 32 taps, which FAILS its own −70 dBc bar**, against −77.4 at
    64; a 3 kHz fundamental at +12 dB of drive goes −56.2 → −75.9 dBc. The window is widened to Nyquist
    and the suite prints the whole map at BOTH taps counts. The stage's CPU roughly doubles with everyone
    else's — **1.24 → 2.43 %RT** — and its round trip goes 31 → 63 samples (+0.67 ms at 48 kHz). What the
    taps do NOT fix there: the map's hot cells (a 3 kHz fundamental at +24 dB reads −40.9 dBc at either
    taps count) are the tube's harmonics folding inside the oversampled domain, which is the FACTOR's
    axis. OrbitCab keeps its own copy of the stage at an explicit 32 and is unaffected.
  - **`TruePeakLimiter`'s delivered-excess characterisation is RE-DERIVED**, as its own header demanded
    ("widen the pass band and this term has to be re-derived"). The worst tone the round trip passes flat
    is now **2fs/5 rather than fs/3**, so the grid-geometry term goes **+0.302 → +0.436 dB at 4×** and
    **+0.075 → +0.108 at 8×** (unchanged at 2×, where fs/3 still dominates at +1.250). The old figures
    were a property of the lowpass, not of the limiter — `mastering` has been running at +0.436 all along.
    2fs/5 is now a witness in the ceiling matrix, and the CHARACTERISATION battery — everything reached
    through `Setup`/`renderAt`, including the reference oracle it nulls against — runs on the SHIPPED
    topology instead of a hardcoded 32. Roughly two dozen guard and plumbing checks still spell 32 out on
    purpose; those are about refusal and resumption, not about the filter. **The one place the sharper filter costs** is the degenerate corner with
    lookahead AND release both at their floors, where re-band-limiting a step rings more: **+0.36 / +0.28 /
    +0.26 dB at 2× / 4× / 8×**, pinned. The ON-GRID bound, the only thing actually promised, is unchanged.
  - **Source-level API break beyond the defaults:** `PowerAmpStage::prepare` gains a fourth parameter,
    so its *type* changes. Ordinary calls still compile (the parameter is defaulted), but anything that
    names the function's type — `void (PowerAmpStage::*)(double, int, int)`, a `std::function` built
    from it, an explicit `&PowerAmpStage::prepare` cast — does not.
  - **`PolyphaseOversampler::prepare` now refuses `tapsPerPhase` above `kMaxTapsPerPhase` (1024) and
    `factor` above `kMaxFactor` (64)** —
    both factors of `N = factor * tapsPerPhase` had to be bounded because the PRODUCT is what overflows a
    signed int before allocating, and either argument alone can do it. `TruePeakLimiter` guarded both at
    its own gate already; `Saturator` passes an unbounded `oversampleFactor` straight through, and
    `prepare(INT_MAX, 1)` on the default taps is UBSan-confirmed overflow followed by a `length_error` —
    a terminate under the wasm tier's `-fno-exceptions`.
  - **`oversampling`'s own suite measured none of this**: it ran entirely at a hardcoded 32 on tones of
    500 Hz and 2 kHz at 48 kHz (0.010 and 0.042 fs), owning the default and never looking where the
    default decides anything. It now pins the stopband, the aliasing and the two-stage droop surface,
    two-sidedly, with its oracle's liveness asserted.

- **BREAKING (API + behaviour), every module with a block-level `process()`:** **the core had four
  different answers to "the caller passed something other than what was prepared", and every one of them
  was SILENT.** A census over all 24 modules found three answers for the width, four for the length, and a
  case that fits neither. The single rule is now **law 11** in `docs/DSP-ARCHITECTURE.md`, and every
  block-level entry point returns its verdict: **`[[nodiscard]] bool process(...)`**, `true` = "accepted
  and honoured in full". That is the API break — a call site that ignored the result now warns, and warns
  in exactly the place where the contract is decided.
  - **The LENGTH is a capacity, not a limit.** `maxBlock` sizes scratch; `process()` chunks, so any
    `n >= 0` is processed IN FULL and the chunked pass is bit-identical to the caller having chunked it
    itself at the same boundaries. Three stages truncated silently before: **`dynamics::NoiseGate`** let
    everything past `maxBlock` out UNGATED — 3840 of 4096 samples at **+89.99 dB** over the gated ones,
    which is 100 % of the construction ceiling (`-floorDb` = 90 dB); **`nam::NamStage`** let it bypass the
    amp model **bit-identical to its input** — 448 of 512 samples on a 64-sample prepare; and
    **`multiband::MultibandProcessor`** dropped the ENTIRE call, not one sample touched. (Dropping
    `NamStage`'s clamp without chunking would have been far worse than the defect: an unclamped `n` is a
    heap overflow in the backend's scratch and a resize — an allocation — inside NAM, on the audio thread.)
  - **The exception, and it is named in the law:** a two-phase API whose first phase RETURNS a buffer of
    `maxBlock` cannot chunk, because the result has to outlive the call. `NoiseGate::analyse`/`applyGain`
    and `EqEngine::captureSectionInput` therefore REFUSE an over-long call rather than truncate it.
    `NoiseGate` also now records how far the curve is valid: phase B used to be bounded by the buffer's
    CAPACITY rather than by what phase A wrote, so a short analysis followed by a long apply multiplied
    the tail by the PREVIOUS call's curve — measured, **90.0 dB of attenuation on material that was never
    analysed**. `NoiseGate::analysedSamples()` is the new accessor; `NoiseGate::prepare()` now returns
    `bool` and honours `maxChannels`, which it used to ignore.
  - **The WIDTH is a limit: `nch > maxChannels` refuses the WHOLE call, before anything moves.** It used
    to process a prefix in eleven modules. A prefix is not the safer half-measure it looks like — the
    surplus channels are unprocessed either way, and processing some of them only hides the fault while
    the processed ones acquire a latency and a gain the others do not: `limiter::TruePeakLimiter` prepared
    for 2 and called with 4 emitted the surplus **+7.02 dB over its ceiling** AND left the two it did
    process **79 samples late** relative to them, which combs at 304 Hz on any fold-down.
    `saturation::Saturator`'s surplus came out **bit-identical to the input**, with no saturation at all.
    **Affected:** `TruePeakLimiter`, `Saturator`, `EqBand`, `EqEngine`, `Dither` (which was bounded by
    `core::kMaxChannels` rather than by its own prepared width), `DeEsser`, `DynamicEqBand`,
    `TransientShaper` (which ignored `maxChannels` entirely), `PowerAmpStage`, `ConvolutionEngine`,
    `MultibandProcessor`, `LinearPhaseEq`, `NaturalPhaseEq`, `LoudnessMeter`, `TruePeakMeter`, `MonoBass`,
    `StereoWidth`, `RigPlayer`, `NamStage`.
  - **A NARROWER call stays legal** (it is the falling-edge mode P18 defined) **except where it is
    meaningless**, and there it is refused OBSERVABLY: `convolution::MatrixConvolver` and
    `MatrixConvolverNupc` need exactly `channels_` planes — a 2x2 matrix needs both inputs to compute
    either output. They used to drop such a call and write NOTHING, so a caller that pre-zeroed its output
    got digital silence and no way to find out (P18 F35). Their width is EXACT now in both directions: a
    WIDER call used to be accepted with the extra planes left dry. This propagates to `CabConvolver`,
    `LinearPhaseEq` and `NaturalPhaseEq`, whose "mono path" on a stereo-prepared engine had in fact been
    doing nothing at all.
  - **BREAKING (behaviour): a call with `n > 0` and NO channels is a GAP in the stream, not a no-op.** It
    spends audio time (the law-8a grid advances) AND every channel at index >= `nch` counts as stopped, so
    P18's falling edge fires — at `nch == 0`, for all of them. Saying "time passed" and "nobody stopped"
    in one breath reopens exactly the defect P18 closed: measured on untouched `main`,
    `dynamics::Compressor` with 5 ms of lookahead, a tone, 4800 samples of zero-width calls, then stereo
    DIGITAL SILENCE emitted **0.280315 out of the silence (-11.05 dBFS)**, last non-zero at sample 239 —
    the whole 240-sample lookahead line, note for note. `Compressor` also gains the NARROWING half of that
    edge, which it never had. `eq::EqEngine` used to drop a zero-width call outright while the same
    `EqBand` driven directly advanced its grid: on one 500 -> 4000 Hz glide and 10240 samples of such
    calls the two had diverged by **11.08 dB**. Negative `n` or `nch` is malformed and refused, never
    clamped into an index.
  - **`prepare()` is binding too, and refuses what it cannot honour** — an observable refusal in
    `process()` is worth nothing if `prepare()` already lied. `convolution::CabConvolver::prepare()`
    silently clamped `numChannels` to 2, after which `process(io, 4, n)` was a well-formed call that left
    planes 2-3 DRY; it returns `bool` now, as do `NoiseGate::prepare`, `Dither::prepare`,
    `DeEsser::prepare`, `DynamicEqBand::prepare` and `TransientShaper::prepare`.
    `MultibandProcessor` gained a `prepared_` flag (its `prepare()` could already fail and `process()`
    ran anyway), and `ConvolutionEngine::prepare()` now clears the participation ledger it left stale.
  - **`neural::Inference` changed:** the concept requires `process(io, nc, n) -> bool`. Any custom
    inference backend must return its verdict.
  - Two integer-overflow fixes that came with the "any `n`" promise: the chunk loops in
    `saturation::Saturator` and `poweramp::PowerAmpStage` advanced by `maxBlock` rather than by the length
    actually taken (signed overflow near `INT_MAX`, which `TruePeakLimiter` had already fixed), and
    `dither::Dither`'s blank counter added `n` to an `int` before clamping it.
  - **A refused `prepare()` writes NOTHING and leaves the object UNPREPARED.** Both halves cost a defect
    while this law was being applied and are now part of it: storing one argument before validating the
    next left a new WIDTH beside an old buffer (`mastering::OfflineRenderer` — a heap-buffer-overflow
    ASan caught), and refusing before reaching an inner `prepare()` left the object ARMED on its previous
    build (`multiband::MultibandCompressor` — `prepare(2)`, a refused `prepare(0)`, and `process()` still
    ran). **Also `[[nodiscard]] bool` now:** `rigplayer::RigPlayer::prepare` (it CLAMPED the width — the
    literal defect this law describes), `multiband::MultibandWidth::prepare` (its ceiling is 2, not
    `kMaxChannels`, because its band is a fixed stereo stage), `analysis::TruePeakMeter::prepare`,
    `analysis::LoudnessMeter::prepare`, `dynamiceq::LaneDynamics::prepare`,
    `mastering::OfflineRenderer::prepare`.
  - **`stereo::MonoBass` and `stereo::StereoWidth` now refuse `process()` before `prepare()`.** Their
    member defaults looked like a valid configuration and are not: MonoBass's crossover has no
    coefficients until `prepare()` runs, so a default-constructed object passes the side band it exists
    to fold at **-6.02 dB where a prepared one kills it to -54.22** — 48.2 dB at 30 Hz.
  - **The gap reaches per-channel state one layer further down than the first pass saw.** A composite
    must PASS a zero-width call to what it wraps, not answer for it: `multiband::MultibandProcessor` now
    also drops the parallel DRY delay (the whole 240-sample line replayed — 0.25 out of digital silence,
    -12.0 dBFS, invisible at mix 1 where the sum cancels exactly) and notifies BYPASSED bands (-11.83
    dBFS); `dynamiceq::LaneDynamics` releases its own lanes (a lane held -12.000 dB through a second of
    gap and pushed the return down by -5.25 dB); `analysis::LoudnessMeter` clears the K-weighting state
    of a stopped channel, like `TruePeakMeter`'s history (momentary read **-29.19 LUFS** into digital
    silence, against -120.00). New: `analysis::KWeightingFilter::resetChannel(c)`.
  - **`eq::EqEngine::captureSectionInput` is the law's other two-phase API and owes the same:** a width
    above the prepared one is refused rather than narrowed, and a refused OR empty call is inert — an
    `n == 0` capture used to wipe a valid 32-sample one.
  - **Also newly REFUSED where these used to clamp or run:** `prepare(..., maxChannels = 0)` on
    `eq::EqBand`, `eq::EqEngine`, `saturation::Saturator`, `lineareq::LinearPhaseEq`,
    `lineareq::NaturalPhaseEq`, `multiband::MultibandProcessor` (hence `MultibandCompressor`) — a JUCE
    host with a disabled bus passes exactly that; `prepare(..., maxBlock/blockSize <= 0)` on
    `dynamics::NoiseGate`, `multiband::MultibandProcessor`, `mastering::OfflineRenderer` and
    `rigplayer::RigPlayer`; and `process()` BEFORE `prepare()` on `dither::Dither`,
    `deesser::DeEsser`, `dynamics::TransientShaper`, `dynamiceq::DynamicEqBand`, `dynamiceq::LaneDynamics`
    and `dynamics::NoiseGate`, which used to run on their member defaults.
  - **`limiter::TruePeakLimiter` now treats a zero-width call as the channel-count change it is** and
    resets — measured, gain reduction from -5.08 dB to 0.00 and a 48-sample hole on the return. That is
    the discontinuity its header already accepts for a width change; it is called out here because it is
    a THIRD answer to law 11c's open question, arrived at by the width rule rather than chosen.
  - **`convolution::PartitionedConvolver::process` and `NonUniformConvolver::process` also return
    `bool`.** They are mono and take no channel count, so "every block-level entry point" would otherwise
    have been a claim with two exceptions.
  - **`nam::NamStage` has no falling edge**: a stereo model's second instance keeps its receptive field
    across a gap (measured, 0.5 out of digital silence at a 16-sample field). Law 11a's guarantee does
    not reach it, and the fix belongs with the model half — recorded, not silently claimed.
  - **`stereo::MonoBass::prepare` and `StereoWidth::prepare` return `[[nodiscard]] bool` and honour the
    `maxChannels` they take** — they used to accept it and ignore it, which is the thing law 11(b) calls
    lying about a contract. `mastering::MasteringChain::prepare` propagates the MonoBass verdict.
  - **Migration:** check the return value. **Two-phase callers first:** `NoiseGate::analyse` now REFUSES
    a block longer than `maxBlock` where it used to truncate, so a consumer that hands it a raw host
    block and ignores the verdict gets the WHOLE block ungated instead of only the tail. Use the fused
    `process()`, which chunks, or chunk at the call site. `if (! stage.process (io, nch, n)) { /* your geometry is wrong */ }`
    A consumer that today passes more channels than it prepared for, or a block longer than a
    capacity-bearing `analyse()`, was already getting broken audio — it was just not being told.

- **fix(core, eq, stereo):** **law 8's denormal flush ran once per `process()` call, so the caller's
  block size decided where a numerical event landed.** New `core::StateGrid` — a phase counter over
  AUDIO samples, period 64, re-anchored by `reset()` — and `eq::EqBand` (hence `EqEngine`) and
  `stereo::MonoBass` now do their periodic maintenance there. The rule is written up as **law 8a** in
  `docs/DSP-ARCHITECTURE.md`; two modules had already reached it independently (`analysis::LoudnessMeter`
  per 10 ms sub-hop, `saturation::Saturator` per sample).
  - **BREAKING (behaviour), `eq` and `stereo`:** output moves. On settled parameters the change is
    confined to the sub-audible: 8.3 % of samples over a 27 M-sample sweep, worst difference 1.95e-13,
    loudest differing sample −127 dBFS. With a parameter ramp in flight it is large and intended
    (worst 1.25 full scale) — the smoothers now advance by exactly one grid period per tick instead of
    by the whole call, which is what makes a ramping band slicing-invariant at all.
  - What it buys, measured: a whole-file render of a bell into silence went from **38 522 of 40 000 tail
    samples differing** from a one-sample-at-a-time render to **0**; from **37 678 subnormal tail
    samples** (the 10–100× stall on any CPU without hardware FTZ) to **0**; and one `+Inf` input sample,
    which used to poison **479 900 of the next 480 000 samples** of a whole-file render, is now healed
    within one grid period (28 samples). A host with a block SHORTER than a period keeps its immediate
    recovery: the `isfinite` half of the flush still runs at the end of every call
    (`eq::Biquad::healPoison()`, `eq::Svf::healPoison()`, `eq::Crossover2::healPoison()` — additive).
  - **BREAKING (behaviour), `dynamiceq::LaneDynamics`:** it drives `eq::EqBand` in 16-sample control
    chunks, so the band's STATIC freq/Q/gain glide used to advance every 16 samples and now advances
    every 64. Measured on a 500 → 5000 Hz +12 dB edit with 30 ms smoothing: designs per 100 ms
    300 → 75, largest single step 2.15 → 3.22 dB — still finer than any real host block gave before
    (a 512-sample host took 10 steps of 11.8 dB), but it moved, and a dynamic band is where edit
    smoothness is most visible. The gain DELTA itself is unaffected: it is an arrival, still consumed
    at the producer's 16-sample cadence.
  - **BREAKING (behaviour), `stereo::MonoBass`:** settling into the full-wide bypass is now decided per
    sample instead of at the top of the next call. The samples in between used to take the M/S round
    trip, which is not the identity in float — 1 LSB of 24 bit (5.96e-08) on 22 953 samples of the
    re-slicing sweep, and it made the output depend on where the caller cut.
- **BREAKING (behaviour), `eq::EqBand::reset()`:** it is now a real STREAM RESTART — it re-snaps the
  freq/Q/gain smoothers, redesigns, and clears `initialized` so the first parameter write after it snaps,
  exactly as the first write after `prepare()` does. It used to clear filter state and leave the
  smoothers mid-glide, so a second render of the same programme started from a different design: measured
  **0.51 full scale** against a freshly prepared chain. `mastering::MasteringChain` carried that as a
  `forceSnap_` workaround (writing every band with its lanes off and then writing them back); the
  workaround is deleted, and the chain is bit-identical without it across 20 scenarios × 1 152 000 samples.
- **NEW API, `eq::EqBand::clearAudioState()` / `eq::EqEngine::clearAudioState()`:** a STOP — clear what
  the previous audio left behind (filter memory, participation ledgers, dynamic seams, design-key caches)
  and leave the parameter epoch and the grid phase alone. This is what a consumer skipping the engine for
  a while actually wants at a bypass edge; `reset()` there would now also snap every ramp in flight and
  turn the next parameter write into a hard step. **`mastering::MasteringChain` was updated; an external
  consumer that calls `EqEngine::reset()` on an EQ-off edge (OrbitCab's `AmpEq::process`) should switch
  to `clearAudioState()` to keep today's behaviour.**
- **fix(core):** `core::Smoother::advance(n)` memoises `pow(coeff, n)` on `n`, dropped whenever `coeff`
  changes. Bit-identical by construction (`pow` is a pure function) — pinned by a test, and by the
  bit-exact null of every existing consumer — and it is what keeps a transcendental per parameter per
  grid period off the audio thread.

- **fix(eq, saturation, dynamiceq, deesser, poweramp, convolution, multiband, dither):** **a filter that
  is not called does not decay — it FREEZES**, and replays a signal from before the gap when it is called
  again. Eleven instances of one shape: a cell of per-channel or per-lane recursion behind a gate the
  surrounding object keeps running through. All eleven were present in every release up to this one.
  Measured out of DIGITAL SILENCE on the input: `ConvolutionEngine` 7.44e-01 (**−2.6 dBFS**),
  `MultibandProcessor` 1.55e-01 (−16.2), `DeEsser` in SplitBand 6.25e-02 (−24.1), `PowerAmpStage`
  1.67e-02 (−35.5) on its own knob gates and 0.649 (−3.8) on a channel change, `EqBand` **+8.12 dBFS**,
  `LaneDynamics` 11.97 dB of unearned gain reduction, `DynamicEqBand` 0.388 eleven samples in,
  `Saturator` 0.9337, `EqBand`'s swept branch 0.690, `Dither` six LSB of 24 bit.
  - The rule is **the falling edge of PARTICIPATION**, per cell: a cell that ran on the previous
    sample-bearing call and does not run on this one loses its sample memory, whatever stopped it. NOT
    "the channel count grew" — `EqBand`'s single-signal lanes gate on `nc == 2` exactly, so `2→3→2`
    retires them on the INCREASE and revives them on the DECREASE, and the compressor's `nc > lastNc_`
    rule never fires there. Only sample memory is dropped: gain seams, smoothers and applied-value caches
    are current CONTROL and are preserved.
- **BREAKING (behaviour):** the gates above stop and restart at a CONSTANT channel count too, so fixing
  them changes what those transitions sound like. Named, because each is now pinned by a test:
  - a whole-band idle in `eq::EqBand` clears signal history only — it used to call the full `reset()`, so
    an inert companion lane on 0 dB decided whether another lane's commanded gain survived. The full
    reset, seams included, now belongs to explicit `reset()` alone.
  - switching a lane off, flipping `swept`, or turning `dyn.on` off now clears that path's delta filter
    state; `Saturator`'s Asym↔symmetric toggle clears its DC blocker; `PowerAmpStage`'s presence, depth,
    load and iron gates clear their filters and its sag envelope.
  - a parameter ramp in `eq::EqBand` now advances while the band is idle. It always did for lanes that
    were not running; the fully-idle band was the one case that fell out of the rule, so the same edit
    landed at two different times depending on whether an unrelated lane happened to be on.
  - a parked lane's PROGRAMME ESTIMATE in `dynamiceq::LaneDynamics` is now duration-aware: kept below a
    quarter of the estimator's own averaging constant, fast-adapted below four times it, discarded beyond.
    Keeping it unconditionally cost 17.87 dB for 3.36 s against a lane that heard the change — in boost
    mode, 18 dB of unearned BOOST. The trade is named rather than hidden: after a long SILENT park the
    first loud material is no longer ducked.
- **BREAKING (eq, source):** `EqEngine::prepare` is now `[[nodiscard]] bool` and REFUSES a sample rate it
  cannot honour; `EqBand::prepare` returns `bool` and a refused band stays inert. It validated nothing
  before, and `fs` reaches every coefficient through `tan(pi*f/fs)`: measured, `prepare(0)` and
  `prepare(NaN)` each put 63 of 64 output samples non-finite on a 0.25 input, and `prepare(1.0)` does the
  same to a HighPass or a Notch. Only the mastering chain was protected, because it validates at its own
  entrance. This is the second half of the fix `Saturator::prepare` received in v0.26.0.
- **BREAKING (eq, behaviour):** `EqEngine::captureSectionInput` records the WIDTH it captured, not only
  the length, and hands back `nullptr` for the columns outside it — the same refusal it already made for
  an over-long block. They used to point at the PREVIOUS block's audio, so a consumer reading one column
  too far detected on a signal that was not that block's, silently and on plausible data.
  `sectionInputChannels()` reports the width.
- **feat(eq, oversampling):** `resetChannel(int)` on `Svf`, `Crossover2`, `MultibandSplitter` and
  `PolyphaseOversampler` — purely additive. Clearing one channel's state without restarting the channels
  that never left is what the isolation half of the fix above needs, and none of them could do it.
  `PolyphaseOversampler::resetChannel` clears the ring POSITIONS with the samples.
- **perf(eq):** `Svf::flushDenormals` visits the PREPARED channels instead of `kMaxChannels`. The columns
  past `ch` are never written and `reset()` zeroes them, so visiting them was a no-op that cost the same
  as real work — since v0.26.0 the body is two `isfinite` tests and two stores per column rather than a
  flush, so a mono `Svf` paid for sixteen, once per block, in the primitive every band, lane, crossover
  and probe calls.

## v0.29.0 — the loudness tag answers for the slot that is sounding, and the chain gives the same bits however you cut it (`rigplayer`, `mastering`)

- **feat(rigplayer):** `soundingLoudness()` — the tag of the model actually carrying the sound, with
  `tagged`, `blended` and `slot` beside the number. `modelLoudness()` / `modelHasLoudness()` always
  read slot 0, and slots are handed out by knot parity: on an odd capture the whole sound comes from
  slot 1 while slot 0 holds a neighbour at zero weight, so a host drawing "the tag of the model
  sounding" drew the silent one's — or called a tagged model untagged. Three reviewers found it
  independently.
  - The slot is chosen by the APPLIED weight, not the requested one: while a model is loading or
    warming, the old slot is what is audible, and a face has to name what is heard.
  - `blended` is set only when the weight is strictly between the ends AND the two slots hold
    different files — at rest both slots hold the same capture, and one number is then the whole
    truth. During a real crossfade there are two tags in the sound, and an API that returns one
    without saying so invites the same quiet wrongness this readout exists to end.
- **BREAKING:** `modelLoudness()` and `modelHasLoudness()` are REMOVED rather than kept as aliases.
  Their documented meaning is "slot 0" — the defect itself — so an alias would go on answering wrongly
  for anyone who did not recompile.

<!-- `mastering` shipped in this tag (PR #132, merged before the v0.29.0 release PR #134) but its
     entry sat under Unreleased until it was moved here, so this file claimed the module was not
     released for three versions. The text below is that entry, unchanged. -->

- **feat(mastering):** a new module, and the first one that is a **composition** rather than a stage:
  `MasteringChain` — `gain → EQ → [M/S mono-bass] → compressor (optional internal sidechain HPF) →
  [soft clipper] → gain → true-peak limiter → dither` — plus `OfflineRenderer` over it. No new DSP; every
  stage already shipped. What is new is the part composition kept getting wrong.
  - **A FIXED INTERNAL QUANTUM, not a promise.** The chain never hands a caller's block boundary to a
    stage: it buffers and calls every stage with exactly `internalBlock` samples. That is not tidiness,
    it is the only thing that makes "same input, same output, whatever the block size" TRUE here.
    Measured on this tree: `eq::EqBand` and `stereo::MonoBass` flush filter state once per *call*, so
    renders at block 4096 and block 1 differ in 4721 and 8908 samples; a band parked at 0 dB moves that
    divergence out of the tail and into the programme (2083 samples of a 50000-sample tone, from sample
    24); `dither::Dither`'s auto-blank then amplifies 1e-15 into **three LSB of a 24-bit master** by
    deciding whether the tail is exactly zero; and the compressor is not exempt either — a key of
    `[1.3e-15, 1.5e-12]` with the RMS window set so the follower coefficient is exactly 0.5 gives
    0.478396237 in one 2-sample call against 0.478396297 in two 1-sample calls. Costs `internalBlock`
    samples of declared latency and, measured, **no CPU at all**: 2.26–2.30 %RT flat for every quantum
    from 16 to 4096, so the size is chosen for latency and never for speed. It also closes a Law 8 hole
    nobody had looked at: a whole-file call left 37678 of a 38000-sample tail SUBNORMAL, because
    `FlushToZero`'s "a block is too short to re-traverse the gap" is false for a big block.
  - **Latency is READ BACK from the stages, never computed.** A chain that derives the number itself can
    agree with its own bypass aligners while disagreeing with the stage — `Compressor::prepare(..,
    maxLookaheadMs = 50)` caps the lookahead at 2400 samples, so a chain asking for 60 ms and computing
    `lround(2880)` would hold both at 2880 and pass a bypass null while sitting 480 samples out. The
    suite finds the delay by SEARCHING for the shift that nulls, and measures the active chain's group
    delay from a DFT bin's phase at 25 Hz — a frequency whose period exceeds twice the delay, so the
    answer is not ambiguous modulo a period.
  - **Bypass is decided per stage by measurement, not by uniformity.** The compressor gets an exactly
    transparent WARM bypass out of its own curve (`ratio = 1` → slope exactly 0 → gain exactly `1.0f`,
    sign of zero included, detector still tracking). The clipper and the limiter are skipped with a
    `core::DryAligner` holding their PDC — the clipper's own `mix = 0` is bit-exact for ordinary audio
    but normalises `-0.0f`, and the limiter has no bypass at all: at a ceiling of +60 dBTP it still costs
    the 0.90·Nyquist round trip. Latency never moves when a bypass is toggled.
  - **`tapsPerPhase` defaults to 64**, against the stages' own 32. The prototype loss is a round trip, so
    it doubles in dB, and clipper + limiter in series double it again: at 44.1 kHz that is **−1.549 dB at
    17.6 kHz**, −6.033 at 18.5 and −16.131 at 19.4. At 64 taps: +0.000 / −0.610 / −10.182, for +64
    samples. The table is pinned in the suite — the limiter's own tests stop at 0.357·fs and did not
    notice that its header's droop figures are one filter pass labelled as the round trip.
  - **The channel count is EXACT.** `process()` refuses any count but the prepared one, touching neither
    buffer nor state, so `[A, refused, B]` is bit-identical to `[A, B]`. The stages disagree among
    themselves — the compressor refuses a wider call, the limiter, saturator and EQ process a PREFIX,
    mono-bass ignores anything that is not exactly two — and reconciling that from outside is the
    field-mapping that falls out of step.
  - **One input gate**, in the shape `Saturator` and `TruePeakLimiter` already use, so a poisoned sample
    is bit-identical to the sanitised one it stands for, whichever stages are on. And the chain validates
    its own sample rate positively: `Saturator::prepare(NaN)` returns **true** and emits NaN, and
    `EqEngine::prepare` does not validate at all.
  - `OfflineRenderer` is defined by a formula rather than a description: `out[n] = y[n + D]`, where `y`
    is the chain's output for the input followed by `D` zeros. Same length as the input, aligned, and the
    last `D` frames present — the ffmpeg chain this replaces never emitted them, so a click 4 ms before
    the end of a file disappeared.
  - **Coverage: 193 checks across two binaries, and 27 of 27 mutations of the module are caught.** The
    second binary exists because the obvious latency test is unsound on its own; everything in it
    re-derives the answer from outside the chain. Two of the mutations found real defects while the suite
    was green — `reset()` resumed an interrupted EQ parameter ramp (0.51 of difference, full scale,
    against a fresh chain) and a bypass toggle read a cold aligner.
- **feat(stereo):** `MonoBassParams` + `MonoBass::setParams()` / `params()`. Additive; `setParams` calls
  the three existing setters, so every clamp and rejection is unchanged and the audio is bit-identical
  (pinned). `params()` reads back the RESOLVED values — a corner below 20 Hz comes back clamped.
- **docs(ADR §4):** the "no cross-module glue in the core" rule is amended to what the tree has actually
  been doing since `dynamiceq`: a composite belongs here when IT is the unit under test; the voicing
  stays in the product.

## v0.28.0 — the host states the whole number, and the player does the subtracting (`rigplayer`)

- **feat(rigplayer):** `setHostInputDb()` / `setHostOutputDb()` — a host with a fader of its own now
  states the WHOLE level it wants a device played at, and the player uses it in place of the pack's
  `chain[].input_db` / `chain[].output_db`. `std::nullopt` (the default) hands the level back to the
  pack, so a plugin that only plays packs is unaffected and needs no call.
  - v0.27.0 left the arithmetic with the host: it sent `hand − what the pack says`, applied outside
    the player. That cannot be made correct. The subtraction has to know which pack is loaded at the
    moment it happens, and a host reads that from its own document — which diverges from the loaded
    pack whenever somebody edits it, a rebuild is in flight, or a device is switched mid-build. Each
    of those left the level silently wrong, in the same class as the double application the two keys
    were added to end. The player is the only place that cannot disagree with itself about which pack
    it holds. `hostInputDb()` / `hostOutputDb()` read the hand back.
  - The hand survives a `load()`: it belongs to the bench, not to the pack.
- **fix(rigplayer):** a pack swap no longer publishes unity levels in the middle of itself. `load()`
  called `unload()`, which resets both gains, and republished the real ones only at the end — a
  callback landing in between heard neither pack's level. The stage's levels are now published as
  soon as the stage is known, and one function decides what is applied.
- **test(rigplayer):** a fixture that is not a pure scalar. Every model in these tests was a linear
  gain, and through a pure gain the two levels are indistinguishable — moving where the player applies
  them would have failed nothing. A Linear NAM WITH BIAS adds an offset that input scaling leaves
  alone and output scaling takes down, so the release's central claim (input is drive, output is
  volume) is now something a swap would break.

## v0.27.0 — a pack carries its own two levels, and normalizing is the default (`rigplayer`, `nam`)

- **feat(dynamics):** the gain-reduction path became **one object**: `GainReductionPath` carries
  detector → curve → ballistics together, `CompressorParams : GainReductionParams : DetectorParams`
  completes that layering, and `GainReductionTap` publishes the per-sample gain reduction. Alongside them
  `dynamics::offline::{QuantileHistogram, EnvelopeAnalyzer, ThresholdSolver}` — the detector-domain
  analysis that finds the threshold delivering a wanted gain reduction by MEASURING rather than by
  inverting a curve. None of this was declared when it shipped; the refactor itself was bit-identical to
  the previous compressor across 23 328 configurations × 210 million samples.
- **BREAKING (dynamics, source):** completing that layering moved `thresholdDb` and its neighbours into
  `GainReductionParams`, so **`decltype (&CompressorParams::thresholdDb)` is now
  `double GainReductionParams::*`, not `double CompressorParams::*`**. Exact-type reflection,
  serialisation tables and any `T C::*` template argument spelled with the derived class stop compiling.
  The header says so at `Compressor.h`; this file did not. See also the v0.26.0 note above, which is the
  first half of the same change.
- **feat(rigplayer):** a pack states **how hard it is fed and how loud it leaves** — namz 4.1.0's
  `chain[].input_db` and `chain[].output_db`, applied by the player and by nothing else. Packs are not
  balanced against each other (a Big Muff leaves some 12 dB louder than a clean preamp, and a boost is
  hotter still than the preamp it feeds) while inside a pack the models sit within ±0.4 dB of one
  another, so the level that differs belongs to the DEVICE. Until now the only place to put it was
  `files[].input_db`, spread across every entry, which made one key mean two things at once.
  - The input level goes **first, ahead of the dry copy**, and deliberately NOT into `chainGain_` where
    `extendDb` lives: a blend knob mixes one guitar with itself, and feeding the models less while the
    dry side is fed as ever turns the mix into two instruments at two volumes. `extendDb` stays where it
    is — it is a trick played inside the pack past the top capture, which the dry path leaving at the
    input jack never saw.
  - The output level goes **last, after the mix**. One number for the whole stage, and applied any
    earlier it would ride the wet side alone and move the blend the pack states for that position.
  - Neither is behind `setInputTrims()`, which keeps gating exactly what it always gated:
    `files[].input_db`, the trim of one alias against its neighbour. `stageInputDb()` /
    `stageOutputDb()` read back what the player is applying — **a host with a fader of its own must
    send only its deviation from these**, or the level lands twice.
- **change(rigplayer):** `setNormalize()` now defaults to **true**. A model's `metadata.loudness` tag is
  a contract, not a listener's option: with it off, every capture plays at whatever level the hardware
  happened to give, and no two packs can be compared at all. A host that wants the old behaviour must
  now ask for it.
- **build(nam):** namz pinned at **v4.1.0**, which is where the two keys are.

## v0.26.0 — a tier that stops being aspirational, an external key for the compressor, and six kernels that finally arrive at zero (`build`, `dynamics`, `analysis`, `saturation`, `tools`)

- **feat(build):** the **`wasm-audio` tier is a gate**, not a paragraph. DSP-ARCHITECTURE.md §2 had
  named it since it was written and always marked it aspirational. There is now a `wasm-audio` CMake
  preset and a CI job that builds every default module for wasm32 with exceptions and RTTI off and no
  pthreads, runs the whole suite in node, and audits every emitted artifact — green under a checked
  `SAFE_HEAP` configuration too. What it does and does not prove is written down in
  `docs/WASM-AUDIO-TIER.md`: a thread is not a build error, and law 8 is not covered by it.
- **fix(core):** the one `throw` that closed the exception-free tier. `core/Fft.h`'s
  `SeamAllocator::allocate` held the only `throw` in the shipped core, and `-fno-exceptions` makes a
  `throw` a hard PARSE error — so the tier the ADR promises could not be built at all. Guarded the way
  libc++ guards `__throw_bad_array_new_length`: throw where exceptions exist, `std::abort()` where they
  do not.
- **BREAKING (dynamics, source):** `CompressorParams` stopped being a flat struct and became
  `CompressorParams : GainReductionParams : DetectorParams`, so its own fields moved into base classes.
  Three things break, and this release did not say so:
  - **Designated initializers stop compiling.** `CompressorParams{ .ratio = 4.0 }` is now
    *"field designator 'ratio' does not refer to any field in type 'CompressorParams'"* — a designator
    may only name a DIRECT member, and `ratio` lives in `GainReductionParams` now.
  - **Positional brace-initialization stops compiling too**, loudly rather than silently: brace elision
    fills the base first, so the pre-existing `CompressorParams{ -12.0, 4.0 }` now tries to initialize
    `DetectorParams::Detector` and `DetectorParams::LinkMode` from doubles. Loud is the good outcome —
    the alternative would have been the same spelling quietly assigning to different fields.
  - **`std::is_standard_layout_v<CompressorParams>` is now false** (it was true), which matters to
    anything doing `offsetof`, C interop, or exact-layout serialisation. It remains an aggregate and
    trivially copyable.
  - Assigning to the fields by name (`p.ratio = 4.0;`) is unaffected, which is how most callers write it
    and why this went unnoticed.
- **feat(dynamics):** an **external key for the compressor**. On a mastering bus the kick and the bass
  decide the gain reduction of the whole mix — defect number one of the ffmpeg chain this core
  replaces, and it would have been reproduced exactly, because the ADR rightly forbids putting an EQ
  inside the compressor. So the EQ does not come in: the input goes out. Eleven ways the module could
  not be trusted with a key are closed with it.
- **fix(dynamics):** the limiter's promise matches its measurement — **twenty** contract defects, against
  the six the task named. The worst was on no list: only the channels passed to `process()` advanced
  their state, so a channel-count change mid-stream went **+19.8 dB** over the ceiling. Topology moved
  to `prepare()`, input gated, oversized blocks chunked, floors stated.
- **fix:** **law 8, finished** — six kernels that never actually arrived at zero, not the three the
  earlier pass left open, and the entry recorded as harmless was the most expensive of them. Stated
  properly the mechanism is not "denormals": `x <- t + r*(x-t)` never ARRIVES, because a residual of
  `k` ulps maps to itself for every `k <= 0.5/(1-r)`, so the decay stops dead while the state is still
  finite and keeps feeding whatever is downstream.
- **fix(analysis,multiband):** silence stops costing **54x** more than music. `KWeightingFilter` was the
  last kernel in core+analysis+oversampling with no software denormal flush; on zero input its two
  TDF-II biquads reach a subnormal pair that maps to itself exactly and freeze there — measured out to
  600 s. On a CPU with no hardware FTZ, which is exactly what a browser gives, that is what silence
  then costs.
- **fix(analysis):** one bad sample stops making the loudness meter **lie quietly**. K-weighting is an
  IIR, so a single NaN made its state NaN forever, every later 400 ms block energy NaN, and
  `NaN > absT` false — so the absolute gate silently dropped all of them and the meter averaged only
  the part that predated the NaN, reporting a healthy, plausible number for a programme it had stopped
  measuring.
- **feat(analysis):** the **pre-gate block energies**, for a comparison that cannot jump. Integrated
  LUFS is discontinuous in its own inputs — BS.1770's gates are strict comparisons, so a block within
  ~1e-12 of a threshold flips inclusion between two builds and moves the reading by ~0.01 dB. No
  scalar tolerance on a gated quantity is a robust equivalence measure; the energies before either
  gate are continuous in the input samples, and cross-toolchain work needs that surface.
- **test(analysis,tools):** true-peak gets its **first external judge**. Both paths had been tested only
  against their own documented behaviour; EBU Tech 3341 §2.6 defines an acceptance envelope nobody here
  wrote, and it is now applied — tests 15-23 of the official EBU Loudness Test Set v05 pass on both
  implementations.
- **fix(build):** **law 10** — say what FP contraction we want instead of inheriting three answers. The
  arm64 Linux CI row went red on its first run on untouched `main`, which is what it was added for:
  `a*b + c` may fuse into one FMA with a single rounding, and the toolchains disagree about when.
  `-ffp-contract=on` is stated; it costs the shipping tier nothing, because it is already clang's
  default.
- **test(saturation):** a noise floor the build **measures for itself**, and the DC blocker's pole and
  numerator, neither of which was tested anywhere. The poison-containment floor had been a constant
  fitted twice to whichever toolchain last went red; below it the assertion has no defined answer,
  because the state residual is a random walk with an absorbing zero. It is now a control run — every
  input sample moved one ulp — because a real build's own noise floor reads 20 ulp = 1.19e-06, above
  the constant the line used to carry. Mutation testing then showed the corner frequency could be
  hardcoded and the numerator's zero moved without a single one of 197 checks noticing; both are
  asserted directly now, two-sided.
- **feat(tools):** `fcore::Probe` — one measurement body for the native reference and for wasm — plus
  the C ABI, build recipe, parity harness and page behind it, and an artifact audit that proves
  no-threads from the binary and keeps the two JS sides from drifting.

## v0.25.0 — a knob that clicks states its filter, and NAM is pinned to a release (`rigplayer`, `nam`)

- **feat(rigplayer):** a switch shipping its bands **per position**. namz 4.0.0 gave the format
  `positions[].sections` — each position states the filter it IS, with no travel law, because a
  switch's positions are words with an order and no angle for a gain to travel on. The player had no
  path to it: such a knob went down the curve path, met an empty grid, matched `0 == 0` on the length
  test and bypassed the filter. Flat, with every check along the way passing. `sectionsAtValue` builds
  that position's bands by name and `bandsPerPosition` sends the knob to the band path.
- **fix(rigplayer):** `setSwitch` refuses a value no position declares. It used to store any word,
  report success and read it back on the knob while the sound was the reference — a typo of one letter
  was accepted and heard as nothing. A dial still takes any degree of its travel, swept or not.
- **build(nam):** NeuralAmpModelerCore is pinned to **v0.5.4**, a release rather than a commit.
  `b5a68c3` was the tip of upstream main when OrbitCab pinned it and stopped being so forty-five
  minutes later; it was inherited here with the NAM path and never moved. A real model through a
  second of deterministic signal is **bit-identical** across the two pins, so nothing already captured
  sounds new.
- **build:** namz is fetched at **v4.0.0**, the schema the player now speaks.

## v0.24.0 — the constant-Q analyzer for half the CPU, and silence stops costing more than sound (`analysis`)

- **feat(analysis):** `MultiResSpectrumPaneFast` — the same pane as `MultiResSpectrumPane`, computed
  differently. A `perf` profile said **43.98 % of a tick was inside libm**, and the band integration
  everyone assumed was the cost was 7.6 % while the FFT itself was 4.6 %. Three exact changes remove
  35 349 of the tick's 37 111 transcendental calls: the peak trace is kept in **power** instead of dB
  (which had cost a `log10` per bin to make the value the peak law compares against and then an `exp`
  per bin to undo it — a transform and its own inverse); each column's geometry (owning tier, seam,
  blend, fractional bin edges, display tilt) is derived **once into a fixed-size plan and shared by
  the fill and the peak**, which previously each computed all of it; and DC and Nyquist are peeled out
  of the magnitude loop. **1.82× on an M5 Pro and 2.33× on an i9-13900H with pffft** (1.20× / 1.34× on
  the scalar FFT — the same absolute saving, since none of it depends on the transform). On x86 the
  fast pane now costs less than a single classic 16384 pane. It is a **sibling**: `MultiResSpectrumPane`
  is untouched, and the two are held together by a paired NULL in which **the fill is bit-identical**
  and the peak is inside the sibling's own float-dB quantisation. Two divergences are pinned rather
  than hidden — a negative `peakFallDb` is clamped instead of reaching an infinity, and below −120 dB
  the band integral's own conditioning (a difference of two prefix sums scaled by the loudest bin in
  the tier) stops the sibling being a reference at all, which a test measures rather than asserts.
  New: `docs/PERF-ANALYZER-MULTIRES.md`, with the profile, the tables, the machines they were taken
  on, and the lossy ideas that were measured and rejected.
- **fix(analysis):** `MultiResSpectrumPane` no longer leaves its state stuck in the subnormals on
  digital silence. The fill smoother is `pw ← (1−c)·pw`, and in float32 that descent does not end at
  zero: at the smallest subnormal `c·pw` rounds away and `pw` stops moving, leaving all ~10 755 bins
  doing subnormal arithmetic on the message thread with nothing setting FTZ/DAZ. On x86 the pane
  therefore got **slower the longer the transport stayed stopped** — 505.4 µs on settled silence
  against 167.6 µs with FTZ forced, a **3.0× penalty**, while on ARM it was free and so invisible.
  The smoothed power is now flushed to a true zero below `1e-30` — deliberately **not**
  `core::flushDenormal`, whose `1e-15` is an amplitude threshold and would erase the −150…−200 dB bins
  this pane deliberately sums. No reading moves. New `tierBinPower` accessor: `tierBinDb` floors at
  −200 and could not tell a bin that reached zero from one stuck at `1e-45`. The classic `SpectrumPane`
  never had this — it smooths dB, which is already floored.
- **docs(analysis):** `ANALYZER-MULTIRES.md` drops two claims the code stopped honouring in v0.22.2 —
  a shelf repeating the last band above Nyquist (it reads the floor), and a −120 dB floor (it is −200).

## v0.23.0 — the loudness meter meets Tech 3341 (`analysis`)

- **fix(analysis):** `LoudnessMeter` momentary and short-term now accumulate on a 10 ms sub-hop
  (M = the last 40, S = the last 300) instead of the 100 ms gating hop, so either window lands
  within 10 ms of any event. EBU Tech 3341's file-based cases 10 and 13 slide a 3 s / 400 ms tone
  in 150 ms / 20 ms steps and expect the maximum to read the tone ±0.1 LU at every offset; a 400 ms
  burst 40 ms off the hop grid read 0.45 LU low, and case 13 failed at 16 of its 20 offsets. **M and
  S readings on transients change** by up to that much; the integrated measure is untouched in
  substance — every tenth sub-hop closes the same 400 ms block at the same 100 ms hop, LRA keeps its
  1 s cadence — and reads as before to the suite's tolerances (not bit-for-bit: 40 partial sums
  where there were 4). The sub-hop ring is 300 doubles; `process()` still allocates nothing.
- **fix(analysis):** `LoudnessMeter` sizes its block store by hops at the prepared rate rather than
  by seconds (a hop is 10 × lround (0.01·fs) samples — 100 ms only where fs is a multiple of 100),
  and blocks that arrive past `maxDurationSec` are counted in the new `droppedBlocks()` accessor
  instead of vanishing under a straight-faced reading. Additive API; nothing throws on the audio
  thread. A caller that must not lose a block sizes `prepare()` for its longest program and checks
  it reads 0.
- **test(analysis):** `felitronics_loudness_conformance_tests` — EBU Tech 3341 (2023) Table 1
  synthesized from the spec's text: cases 1–5 at 48 and 44.1 kHz (I, M and S on cases 1–2), case 6
  (the 5.0 channel weights), cases 9 and 12 (S and M settle on a periodic program), cases 10 and 13
  (S and M at every offset); then what Table 1 lets a wrong meter get away with, each pinned by a
  signal only the property under test can move — the absolute gate isolated from the relative one
  and at its boundary, the relative gate at −10 LU, every channel as power (in phase or not), the K
  shape at both rates against the published 48 kHz coefficients, a burst that tells 75 % block
  overlap from none, chunk invariance, and `process()` allocating nothing. 59 checks. Grown from
  the Looper Cat suite that gated the product's move onto this meter; crew-reviewed (Codex,
  DeepSeek, Antigravity, a fresh Opus mutating the meter).

## v0.22.2 — fill and peak keep the same bins (`analysis`)

- **fix(analysis):** `MultiResSpectrumPane` shares the classic pane's −200 dB floor. Its fill
  counted every positive bin while its peak-hold was floored at −120 dB, so a band of quiet bins
  read higher on the fill than on the peak trace and the two crossed near the plot's right edge. A
  −130 dB tone now reads −130 untilted — below any plot bottom, not flattened — and the peak trace
  never sits below the fill.

## v0.22.1 — the floor is silence, not a shelf (`analysis`)

- **fix(analysis):** both spectrum panes clamped a reading at −120 dB and then added the display
  tilt, so silence came out as a straight line rising at the tilt's slope — with +6 dB/oct it stood
  at −96 dB at 16 kHz, in plain view on a 120 dB range. `MultiResSpectrumPane` now applies the tilt
  to the power and floors afterwards (`readDb (f, fs, tilt, pivot)`); `SpectrumPane` drops its
  internal floor to −200 dB, deep below any plot bottom, so a tilted floor can never surface and a
  bin with real energy at −130 dB keeps it. Both panes freeze the tilt past Nyquist, where the
  columns repeat the last band that fits. Consumers that pinned the classic pane's −120 floor in
  their own tests move to `SpectrumPane::kFloorDb`.

## v0.22.0 — the analyzer reads constant-Q, and rides SIMD (`analysis`, `fftpffft`)

- **feat(analysis):** `MultiResSpectrumPane` — a constant-Q analyzer from several FFT lengths at
  once. One frame from the rolling tap feeds 16384 / 4096 / 1024-point Hann suffixes that share the
  frame's end (the short tier reports a transient first); a reading is the power in a 1/24-octave
  band, integrated over the tier's bins with fractional edges from double prefix sums and normalised
  by the window's measured ENBW — the one quantity two FFT lengths agree on for both a sine and
  noise, which is what makes a seam invisible. Tiers are used where their bin is at least two per
  band (seams 811 / 3246 Hz at 48 kHz), crossfaded in power over a third of an octave; below the
  longest tier's bin the lows are bin-limited and say so. A tier shorter than the frame hop
  Welch-averages as many half-overlapped windows as reach back over it, so a click between two
  frames is never missed; the per-bin smoothing runs on power, not dB, because a log-domain average
  carries a bias that depends on how the bin was fed. Design and physics: `docs/ANALYZER-MULTIRES.md`.
  Crew-reviewed twice (codex, deepseek, Fable — the latter with simulations); property tests, not
  golden files (`felitronics_multires_spectrum_tests`, 166 checks).
- **feat(analysis):** `RollingSpectrumTap::tryPull (dst, order, hop)` reports the samples that
  entered the ring since the previous publish — the hop that happened, which a block boundary or a
  missed UI tick stretches past the one requested. The two-argument form stays.
- **feat(analysis):** `SpectrumPane` becomes `SpectrumPaneT<Fft>` (constrained to the packed-Hermitian
  layout its bin loop reads) with `SpectrumPane` the scalar alias, so every consumer reads as before;
  it gains `reset()` (the next ingest seeds) for a consumer that returns to it after drawing another
  pane. `MultiResSpectrumPaneT<MaxOrder, MaxTiers, Fft>` likewise.
- **feat(fftpffft):** `PffftOrderedRealFft` — pffft's real transform in canonical order, which is
  exactly the packed-Hermitian layout the scalar reference writes (F(0) and F(N/2) in the first slot,
  then interleaved Re/Im, the e^{−jωn} sign), so it advertises `kPackedHermitianSpectrum` and is
  admissible wherever bins are read. Beside the z-order `PffftRealFft`, not instead: the convolvers
  keep their vectorised MAC; the analyzers get their SIMD. Nulled float for float against the scalar
  (DC / Nyquist exact, Re / Im ≤ 2.4e-7 of full scale) with basis vectors pinning every slot, and both
  panes on both backends within 0.02 dB. Per tick on an M-series Mac: a classic 16384 pane 267 → 73 µs,
  the multi-res pane 381 → 126 µs.
- **test:** `felitronics_spectrum_pane_perf_tests` — the panes' cost per UI tick on the scalar FFT,
  printed for the record with loose ceilings; the pffft suite gains the scalar-vs-SIMD comparison.

## v0.21.2 — the pack's input trims become a switch

- **feat(rigplayer):** the per-file `input_db` trims (a linked setting plays its neighbour softer)
  can be switched off — `setInputTrims(bool)`, ON by default. OFF feeds every capture at unity:
  for a library shot at one honest level the stated attenuations only push a capture's drive
  around, and into a nonlinear model a few dB less in is a lot less out. Read on the audio side,
  riding the existing slot ramps — the toggle lands on the next block, click-free.

## v0.21.1 — a stale landing dies with its pack; a failed load is refused

- **fix(rigplayer):** `unload()` right after `deliver()`, with no audio block between: the published
  landing no longer lands a stale model in the wiped law over an emptied stage — it dies before the
  forget is posted.
- **fix(nam):** a load that failed is refused (`BlendState::refused`), not asked for again on every
  block: a file that fails identically every time cost a fetch and a parse per service tick, for
  ever. The refusal lifts when the request names something else for the slot, or on any landing; the
  slot keeps its old capture and counts as at rest, so it may sleep.

## v0.21.0 — a slot at rest goes cold (`rigplayer`, `nam`)

- **feat(rigplayer):** a slot that has stood silent — weight exactly zero, request unchanged — for
  `RigPlayer::kColdAfterSeconds` (2 s) goes cold: its model is not run, not mixed, and stays loaded;
  its delay line is cleared as it falls asleep. The next change of request wakes it warm-up first, by
  the landing path, with the model's own field — no load, nothing unfed heard. A slot with any weight
  never sleeps: between two captures both models run, on one they do not. `setColdAfterSeconds()`
  (zero or less = never), `slotCold()`, `coldBlocks()`. Measured in OrbitAmp's block with the dial
  at rest: two passes for one sound, 4.5 % of a P-core and 14.3 % of an E-core — now one pass.
- **feat(nam):** BlendLaw owns the flag — `BlendPolicy::coldAfterSamples` (default never),
  `BlendState::cold` and `still`, `blendSameRequest()`. A cold slot counts as unfed; the wake is
  `blendLanded()` with the model and the need already held; the rest is counted as played; a slot
  asked for nothing is at rest too.
- **fix(rigplayer):** the models run on the planes they are given, not on the width the player was
  prepared for: a host prepared for stereo that plays one plane no longer pushes the silent second
  plane through both networks (8.8 % → 4.5 % of real time in OrbitAmp's block).

## v0.20.0 — the pack player comes home (`rigplayer`)

- **feat(rigplayer):** a new header-only module, `felitronics::rigplayer` — one device of a
  `.orbitrig` pack, playing: which captures sound for a panel (`namz::rig`'s policy joined to
  `pickBlend` along the gain dial), the crossfade between them (`nam`'s BlendLaw, once per block on
  the audio thread), the tone knobs as the pack describes them (sections → biquads, a curve → one
  minimum-phase FIR per side), a blend knob's dry path, the models' alignment from the pack's
  `lag_samples`. JUCE-free, no thread of its own: a load is a job the host runs anywhere but the audio
  thread (`takeLoadJob` / `run` / `deliver`). Moved, not copied, from OrbitCapture NAM with its tests,
  so the capture app and a plugin play through ONE player; `modules/rigplayer/README.md` is the
  host's contract. Self-gates on `FELITRONICS_WITH_NAM`.
- **build(nam):** the module's own namz pin rises from v1.1.1 to v3.1.0 (`namz_rig.h` with `tone`
  and `lag_samples`).
- **fix(nam):** BlendLaw's exact-zero swap test is spelled `<= 0.0` on a weight clamped to [0, 1] —
  the same zero, and GCC's `-Wfloat-equal` in the header-hygiene gate is satisfied.

## v0.19.0 — a model is built apart from its stage (`nam`)

- **feat(nam):** `NamStage::prepareModel(bytes, sampleRate, maxBlock)` is the heavy half of a load —
  the bytes parsed, both instances built, the rate-match and the prewarm prepared — done on any
  thread but the audio one and touching no stage; `install(prepared)` is the light half, a pointer
  swap on the message thread through the same pending machinery a load uses. A pack player crossing
  a capture paid up to 20 ms of its drawing thread per landing (measured in OrbitCapture NAM,
  2026-08-30); now it hands that work to a worker and installs the result. `loadModelFromMemory` is
  the two calls in a row — one path, two entry points — and behaves as before. A model prepared for
  other numbers than the stage runs at is prepared again inside `install`, at the old cost; the rate
  contract is judged at `install`, since only a stage can. `NamBackend` binds the normalize flag when
  it meets its stage (`bindNormalize`) instead of at construction.

## v0.18.0 — a WAV can be handed over as bytes (`io`)

- **feat(io):** `writeWavMemory` returns the encoded WAV image instead of writing it to a path. A
  library that keeps its audio in a database needs the image itself, and a temp file on the way there
  is a file to lose. This is the encoder that was already inside `writeWav`, lifted out unchanged:
  `writeWav` is now that call plus one `fwrite`, so every rule holds for both by construction — an
  empty result means the request was refused rather than a truncated file written, ragged channels
  are guarded against reading past an end, and a header RIFF cannot represent (more than 65535
  channels, a sample rate that is not finite, data past 4 GB) is refused outright rather than
  silently wrapped into a lying field. A test asserts the file and the memory image are byte for byte
  identical, so the two cannot drift apart later.

## v0.17.0 — one law owns which capture is audible (`nam`)

- **feat(nam):** `BlendLaw.h` — two model slots play the same input, a knob between two captures asks
  for a fraction of each, and this header is the single writer of that number. It comes from the
  OrbitCapture NAM player, where three places wrote it — the blend, the code following the nearest
  capture, and a warm-up gate bolted on later. Each was defensible alone; together they replaced
  models under a live gain and swung the weight across its whole range inside single 10 ms blocks
  (measured on one sweep of a nine-capture device: 62 loads, 401 parked retries, a full 1.000 swing).
  Seven attempts to fix that by ear in one day, three of which made it worse — hence a law with
  properties provable by construction: the weight moves at most one step per block whatever happens,
  a model is replaced only in a slot whose weight is exactly zero, and progress needs no timeout
  because some slot's goal is always a rail. What gives is instantaneous accuracy — the incoming
  model arrives ~132 ms late and the outgoing one carries alone until it does, which reads as a knob
  rather than as a fault.
- **feat(nam):** the mixing arithmetic lives with the law rather than in the host: `blendMix` and
  `blendDelay`, with the properties that make a slip arithmetic instead of a listening exercise —
  weight 0 is the first slot alone, 1 the second, halfway their linear average (two coherent halves
  stay at unity, not +3 dB), a ramp lands exactly on the weight its block ends with, and a delay line
  carries across a block boundary with nothing lost. The host had lost one of its two model calls, so
  the mix was raw DI against a model; that is audible only as loudness rippling with the knob, and a
  law about weights cannot see it.
- **feat(nam):** `NamStage` reports the receptive field a model must be fed before it means anything,
  and the law counts it BEFORE the block rather than after — one block of conservatism instead of a
  few hundred samples of a network that has not heard the last 132 ms.
- **fix(nam):** an unfed slot stays silent even when that means silence. The old "something must
  sound" rule named slot 0, and a `NamStage` with no model is a PASSTHROUGH — so a device change put
  the raw DI on the output, some ten decibels above any normalised capture. The step now carries its
  own gain, ramped by the same clamp as the weight.
- **refactor(nam):** `blendDelay` takes a pointer and a capacity instead of a template over the
  tail's array size, which no caller with a `std::array` or a runtime-sized buffer could satisfy
  without copying the function — which is exactly what the host had done, leaving 45 assertions
  guarding code that never ran.

## v0.16.0 — the plot's own maths comes home, and 6 dB/oct stops diving (`analysis`, `eq`)

- **feat(analysis):** `PlotMap.h` + `SpectrumPane.h` graduate from TabbyEQ's `eqview` incubator.
  `PlotMap` is the log-frequency / dB coordinate map of an EQ plot (freq↔x, dB↔y, both ways);
  `SpectrumPane` is the analyzer pipeline — Hann → real FFT → smoothed dB with peak-hold, resolved
  into liquid log-frequency columns. Header-only, JUCE-free, no plugin includes. A second consumer
  is what unlocked the move: OrbitAmp's EQ links need the same two, and TabbyEQ now reads them from
  here like everyone else.
- **fix(eq):** `matched::lowpass1` no longer dives at Nyquist. Bilinear puts a zero at z = −1, so
  every 6 dB/oct low pass fell to −inf in the top octaves instead of rolling off gently — and since
  the first-order section is the odd member of every variable-slope cascade, orders 6 / 18 / 30 / 42
  inherited it. The replacement is a magnitude-matched one-pole: exact at DC, at f0 (−3.01 dB) and
  at Nyquist (the analog value), its pole taken from the quadratic whose roots are reciprocal, so
  the small root is always the stable one. `highpass1` keeps bilinear — its zero belongs at DC,
  which is where the analog filter has one.
- **test(eq):** the first-order sections are pinned against the analog prototype
  `|H| = 1/sqrt(1 + (f/f0)²)` over 44.1 / 48 / 96 / 192 kHz × 20 Hz … 20 kHz: unity at DC, −3.01 dB
  at f0 (4e−08), the analog magnitude **at Nyquist** (7e−15 — a bilinear section reads −inf and
  fails outright), monotone descent with no ripple between the matched points, and ≤ 0.97 dB from
  the prototype in between.

## v0.15.0 — a magnitude table becomes a filter (`lineareq`)

- **feat(lineareq):** `MagnitudeCurve.h` — dB against a log-frequency grid turned into a runnable
  minimum-phase FIR. A measurement produces a table, a pack ships a table, and every reader of one
  needs exactly this; it lived in OrbitCapture NAM until a second consumer appeared, which is one
  copy earlier than the house rule allows. `logFreqGrid`, `curveDbAt`, `heldOutsideBand`,
  `magnitudeCurveToFir`, `curveAtPosition`.
- Minimum phase via `MixedPhaseFir` at k = 1: the magnitude is exact, there is no pre-ring and no
  bulk delay, so switching a control does not shift the audio in time. It is not a fitted 1-pole —
  measured on real hardware the best fit missed a Big Muff tone control by 15 to 65 dB across its
  30 dB of tilt, so the curve is the only honest description.
- `curveAtPosition` interpolates against the producer's stated `norm`, never the array index: a
  control measured at 0/30/150/300 on a 300-degree dial has norms 0, 0.1, 0.5, 1.0, and treating
  those as evenly spaced puts a knob halfway four to six decibels wrong.
- **No product policy travels with it.** Trust thresholds and what a pack may claim stay with the
  format that owns them; the trusted band arrives here as two numbers.

## v0.14.0 — the auto-first dynamic EQ across placement lanes (`dynamics`, `eq`, `dynamiceq`)

- **feat(dynamics):** `RelativeLevel` — a slow programme-level estimator with an offset, so a
  threshold can be expressed relative to what a band's own region normally does rather than in
  absolute dBFS (meaningless when the same setting is shared by lanes sitting 20 dB apart).
  `BandBallistics` derives attack/release from a band's own `fc` and `Q`: a band's envelope rises
  with its inverse bandwidth `Q/fc`, not with its period, so one 0..1 "deviation" pair means the
  same thing on a 60 Hz band and a 7 kHz one. Both JUCE-free, software-flushed, control-rate safe.
- **feat(eq):** `BandParams` gains point-level `DynParams { on, rangeDb, thrDb, thrAuto, atk, rel }`
  and `EqBand` gains the **delta seam** `setLaneDeltaDb(Lane, dB)` — the band accepts a NUMBER and
  applies it with its own `Svf` inside the lane, after the matched static biquad and before the M/S
  fold. `eq` therefore takes no dependency on `dynamics`. With `dyn.on` false a band is
  bit-identical to one built before dynamics existed.
- **feat(eq):** `EqEngine::captureSectionInput()` + `bandAt()`/`bandCount()`. A dynamics layer must
  detect on the section's own input: in a series chain a band's input is the previous bands' OUTPUT,
  so their moving deltas would modulate later detectors at overlapping frequencies and the chain
  would pump. The engine deliberately does not hard-code the interleaved loop — it hands out the
  bands so a composition layer can run "delta, then band" in chain order.
- **feat(dynamiceq):** `LaneDynamics` — the lanes-aware composition layer: per-lane sidechain probe,
  envelope follower, `RelativeLevel`, `GainComputer` (ratio and knee fixed internally), GR
  ballistics and the seam write, at a K=16 control rate. Reuses the primitives rather than
  `DynamicEqBand`, whose single-SVF topology would forfeit the matched static curve.
- **fix:** the hardening rounds behind the above — a per-sample control path (interval-correct
  coefficients; a block-size-invariance test pins it), gate parity, no live state carried across
  `reset()`, per-lane state drop on disengage, no seed-from-warmup duck at transport start, a
  hot-signal clamp, and `core::fastGainToDb` on the detector paths (max error 0.0001 dB) to keep the
  per-sample path affordable.
- **fix(dynamiceq):** `laneSignal`'s `switch` names the Stereo lane instead of leaning on `default:`
  — `-Wswitch-enum` is part of JUCE's recommended warning set, so the old form fired in every
  consumer that included the header.

## v0.13.1 — nlohmann rides with the exported namz headers (`felitronics::nam`)

- **fix(nam):** namz's rig headers (`namz_rig*.h`) publicly include `<nlohmann/json.hpp>`;
  consumers reaching them through `felitronics::nam` (OrbitCab's v0.13.0 migration) could not
  resolve the include. The module now exports namz + nlohmann include dirs as a usage
  requirement of its shipped headers — NAM/Eigen stay PRIVATE.

## v0.13.0 — CabConvolver + opt-in NAM backend (`felitronics::convolution`, `felitronics::nam`)

- **feat(convolution):** `CabConvolver` — the product-level cab IR wrapper moves in from OrbitCab,
  DSP byte-identical: reference-unity RMS normalization (2 kHz-shaped reference, ±30 dB clamp,
  −60 dB near-silence floor), fixed NUPC schedule (head 128, 50 ms click-free crossfade),
  mono-broadcast / LRDiag true-stereo publishing, latest-wins retry of a load rejected
  mid-crossfade, staged-taps accessor for offline blend analysis. The
  `FELITRONICS_WITH_PFFFT` backend selection now propagates uniformly through the module's
  INTERFACE (fftpffft registers before convolution; hygiene/test targets deduped).
- **feat(nam):** `felitronics::nam` — NEW opt-in compiled module (`FELITRONICS_WITH_NAM`, default
  OFF, CMake ≥ 3.24). `NamStage` is OrbitCab's `cab::AmpStage` ported verbatim onto
  `felitronics::neural::NeuralStage` + `felitronics::core::StreamResampler`: dual-instance true
  stereo, −18 dB loudness makeup with per-model trim, model-rate contract, bounded retire queue
  with a lossless last-wins pending intent. Loads raw `.nam` and packed `.namz` (namz used
  directly on std buffers; NAMZ_IMPLEMENTATION is compiled exactly ONCE here — consumers include
  `<namz.h>` without the define and link this archive). NeuralAmpModelerCore pinned by SHA and
  linked WHOLE_ARCHIVE (architecture self-registration) with PRIVATE usage requirements — no
  Eigen/nlohmann/NAM headers leak into consumer TUs; namz pinned by immutable commit; both
  overridable by fail-loud local-source cache vars. ctest: analytic Linear FIR, namz round-trip +
  v1 wire compat, loudness/trim makeup, rate-contract refusal against a live model, 70+
  frozen-audio swap stress (mirror-publication discipline, deferred clear), true-stereo
  independence, exact 96 kHz latency pin, resampled/truncated IR loads, and a dedicated
  `EIGEN_RUNTIME_NO_MALLOC` gate (Linear + WaveNet allocation-free; LSTM/ConvNet documented as
  upstream-allocating at the pin, load-compatible by design). Consumed by OrbitCab (cab::AmpStage
  / cab::Convolver become aliases) and OrbitCapture NAM's Queue audition player.
- **test(support):** `approx()` NaN-proofed (unordered comparisons now FAIL); CI adds NAM rows on
  all three OSes (first MSVC nam_core build) + the sanitizer job, and keeps default OFF/OFF rows.

## v0.12.0 — decoupled-hop rolling analyzer tap (`felitronics::analysis`)

- **feat(analysis):** `RollingSpectrumTap` — a lock-free SPSC analyzer tap whose snapshot cadence
  (hop) is decoupled from its analysis window size. One rolling ring (max order 14 = 16384) serves any
  FFT size ≤ MaxOrder from a single buffer: `publishIfDue(order, hop)` copies the most-recent
  `N = 1<<order` samples — chronologically, across the ring wrap — into a single-slot immutable mailbox
  and force-publishes when the order changes, so a consumer can offer a **selectable analyzer FFT size
  with a click-free live switch** at a steady UI-rate cadence (overlapping windows for large N, gapped
  for small N) — no per-order ring duplication, one write per sample. `tryPull(dst, outOrder)` reports
  the order the frame was captured at so a wrong-size frame is discarded across a switch. `reset()`
  restarts the producer only and never revokes a mid-pull reader — closing a torn-frame race that a
  reset-clears-ready design would have when `prepareToPlay` runs against a live GUI reader. Tear-free by
  the same acquire/release ownership handoff `SpectrumTap` uses, plus per-frame order/size metadata.
  Header-only, JUCE-free. ctest: variable-N snapshot, warmup gate, hop cadence, forced publish on order
  change, chronological wrap copy, race-free reset. Consumed by tabby-eq's analyzer-resolution feature.

## v0.11.0 — UTF-8 → ASCII romanization for filename slugs (`felitronics::text`)

- **feat(text):** `felitronics::text` — JUCE-free UTF-8 → ASCII romanization for filename slugs.
  `decodeUtf8` (malformed bytes skipped, never guessed), `romanize` (one code point → ASCII: German
  umlauts / ß / Œ, common Latin-1 accents, and full Russian Cyrillic incl. the multi-letter cases
  Ж→Zh, Х→Kh, Щ→Shch, Ю→Yu; hard/soft signs vanish), and `transliterate` (the whole-string fold).
  Extracted from OrbitCapture NAM's `Matrix.h` so every Darwin's Cat product shares ONE romanization
  table — the RAW text stays in metadata; this is for FILENAMES, where cross-platform sync (macOS NFD
  ↔ Windows/Linux NFC) and legacy ASCII parsers demand plain ASCII and a never-empty result.
  Header-only, zero deps, header-hygiene clean under the strict downstream warning set. 41
  falsification checks (every mapping asserted; intentional collisions like е/э→"e" pinned so callers
  detect them).

## v0.10.0 — file I/O module + more OrbitCapture/OrbitCab promotions

Each promotion behind the extraction bar (theory-first falsification tests + adversarial crew
codex/deepseek + NULL where possible), landed as OrbitCapture NAM (the second capture product) and
OrbitCab consume them.

- **feat(io):** new zero-dependency `felitronics::io` module — minimal self-contained WAV
  read/write (`readWav`/`readWavMemory`/`writeWav`/`writeWavMonoF32`), moved from OrbitCapture's
  `oc/wav.hpp`. Crew-hardened: corrupt/truncated chunks are rejected loudly (never clamped),
  `WAVE_FORMAT_EXTENSIBLE` requires its full body + SubFormat GUID, checked chunk advance (no
  32-bit wrap), writer refuses headers it cannot represent (u16/u32 overflow, non-finite rates),
  data must be frame-aligned. WavTests pin every case.
- **feat(measurement):** `PeakClip.h` — standalone `scanPeakClip` (peak dBFS + flat-top clip run)
  extracted from `gateRecording` for sweepless consumers (NAM reamp takes); a NaN breaks a run,
  `clipRunSamples` clamped to >= 1. `detectLeadingSilence` — OrbitCab's forward-scan head-trim onset.
- **feat(core,saturation):** NaN/Inf poison guards + chunking hardening (promoted from OrbitCab),
  with theory-first falsification suites.

## v0.9.0 — RT stream types (`core::RtStreams`) + model guess + the mix-view overlay facade

Three promotions from OrbitCapture, each behind the extraction bar (tests + adversarial crew
codex/deepseek + NULL where possible) — landed as a second capture product arrived to consume them.

- **feat(core):** `RtStreams` — the RT buffer-swap discipline AS TYPES (`AuditionStream`,
  `ConvStream`, `RecStream`): flag DOWN → mutate within fixed/reserved storage → flag UP; std +
  atomics only. Crew-hardened with four pinned fixes: a zero-length publish no longer arms an
  empty looping clip (reader OOB), `RecStream::push` release-publishes `len` (an ARM harvest could
  read an unsynchronized sample), the audition buffer is fixed-size (no reliance on
  assign-within-capacity), `reserve()` guards int lengths. Reader rules (ACQUIRE gates, per-block
  flag reload) and the known one-block quiescence gap (+ the deferred reader-ack epoch design) are
  documented in the header. Suite includes a discipline-respecting concurrent smoke; ASan/UBSan
  clean.
- **feat(measurement):** `ModelGuess` — gear-model detection in a free-form file name against a
  catalog. Conservative contract: exact fingerprint (the entry's last token) beats the 3+-digit
  bare-number fallback; ANY ambiguity → no guess (a wrong guess poisons imported metadata
  forever). Crew fix: locale-FREE tokenizer (`std::isalnum` could admit high-bit bytes under a
  single-byte locale).
- **feat(blend):** `Overlay`/`makeOverlay` — the one-call mix-view facade: per-mic curves (post
  filters + gain), the blend curve (post Master), and the interference column (PRE-Master basis so
  a Master rolloff never reads as phase cancellation, then faded by the Master's |H| — a display
  weight by deliberate product decision, documented). NULLed to the hand-composed primitives at
  machine epsilon; analytical pins (+3.01 dB in-phase twins, deep 180° notch, |H| only
  attenuates). Crew fixes: bit-parity gain expression, no partial overlays on degenerate input,
  non-finite params heal to defaults (the offline convention).

## v0.8.0 — mic-blend engine (`felitronics::blend`) + fine IR alignment (`measurement::XcorrAlign`)

The MIX side of the IR-capture family joins the capture side (v0.7.0) in core — both extracted from
OrbitCapture (written portable by design), so OrbitCab and other consumers share ONE numerical
fingerprint for blending and aligning multi-mic IRs.

- **feat(blend):** a new `felitronics::blend` module — the **offline multi-mic IR blend engine**:
  `StripParams`/`MasterParams` (per-mic gain / phase / fractional time shift / HPF / LPF + master,
  solo/mute audibility rules), per-section-Q Butterworth HPF/LPF, Hilbert-based phase rotation,
  windowed-sinc fractional shift (positive `shiftMs` = delay), and `blendIrs` (weighted sum +
  master chain). Canonical defaults live HERE (80 Hz/24 dB · 8 kHz/12 dB) — consumers must not
  re-declare them (a drifted default silently re-voices saved mixes). Extraction was gated by a
  **byte-NULL** against the app's previous in-tree engine; OrbitCapture's `ocap::` blend names are
  now `using`-shims over this module.
- **feat(measurement):** `XcorrAlign` — **fine time/polarity alignment** of an IR against a
  reference by normalized cross-correlation (`xcorrAlign` / `xcorrAlignSet`, ±maxLag samples,
  fractional result): per-lag normalization with a Cauchy–Schwarz corr ≤ 1 bound, an onset-delta
  guard that REFUSES (corr = 0) when the two onsets sit further apart than the search range,
  polarity from the best normalized lag, subnormal-safe denominators. corr = 0 is the "no
  confident suggestion" contract — callers must leave such a channel untouched.
- **robustness (crew-hardened):** the adversarial consilium (codex + deepseek) hit XcorrAlign
  before merge and found real bugs: a fixed-denominator normalization that could rank a wrong lag
  above the true one AND report corr > 1; a confident wrong lag + false invert when the true lag
  lies beyond `maxLag` (repro'd at corr 0.91); `inf` via `√(eR·eC)` underflow on subnormal
  energies; a window off-by-one (+ the analysis-window floor raised 16 → 64). Each fix carries a
  pinned counterexample test, and the whole search is NULLed against a brute-force oracle over
  randomized signals.

## v0.7.0 — offline IR measurement + display curves (`felitronics::measurement`, `felitronics::analysis::offline`)

The **offline** (message-thread, allocating, double-precision) half of an IR-capture pipeline, extracted
from OrbitCapture (written portable by design) so the capture math lives in core instead of the app. This is
deliberately **NOT** the RT path — real-time consumers still use the float `felitronics::core::fft` seam;
these transform whole ~6 s captures where the numerical floor must sit far below the analog chain's, so they
run in `double`. Every function clamps/heals non-finite params; correctness is oracle- + **numpy-cross-NULL**-
anchored (an independent `numpy`/direct-time-domain recompute nulls the C++ to machine epsilon) and
ASan/UBSan-clean.

- **feat(measurement):** a new `felitronics::measurement` module — an **exponential sine sweep (ESS / Farina)**
  generator + matched inverse (`Sweep`), **Farina deconvolution** with a latency-absorbing onset search past
  the harmonic region (`Deconvolve`), IR post (onset / trim / peak-normalize, `IrPost`), a pre-deconv
  **capture-quality gate** (clip / non-finite / sweep-presence / SNR, `CaptureGate`), and **multi-mic
  common-onset alignment** that preserves the inter-mic comb (`MicSetAlign`). NULL-verified: sweep⊛inverse≈δ,
  known-answer in-band magnitude, and the direct-time-domain convolution vs `numpy.convolve`.
- **feat(analysis):** `felitronics::analysis::offline` display curves — `logMagnitudeCurve` (a 1/N-octave
  RMS-power-smoothed magnitude on a log-f grid, in dB; energy-preserving `10·log10(mean|X|²) = 20·log10(rms)`)
  and `interferenceDb` (where a multi-mic blend cancels or reinforces vs an incoherent power sum). Namespaced
  `::offline` to keep the RT-metering contract of the rest of `felitronics::analysis` intact.
- **refactor(core):** the shared offline double FFT (`nextPow2` / `detail::fftInplace` / `convolve` /
  `magSpectrum`) is **promoted** `measurement` → `felitronics::core::offline` (`core/OfflineFft.h`) now that a
  second offline-FFT consumer (the display curves) exists. `measurement/Convolve.h` re-exports it — the
  measurement API and all its tests are unchanged.
- **robustness (crew-hardened):** an adversarial "break-it" consilium (deepseek + antigravity + Fable) found
  9 edge-case bugs in `measurement` and 4 in the display curves (a tiny/inf smoothing band → `(int)ceil(inf)`
  UB; a non-pow2 `minNfft` → binHz↔bins skew; a huge-finite sample → `ΣM²` overflow → NaN; a top-bin band
  inversion) — each fixed **with a regression test**. Verified false alarms were rejected against the code.

## v0.6.0 — noise gate (`felitronics::dynamics::NoiseGate`)

A new `felitronics::dynamics` primitive: a dual-detection (ISP Decimator "G-String" style) **noise gate** —
architecturally distinct from the continuous `Compressor` (a bistable **Schmitt trigger** with hysteresis +
hold, a **LINEAR-fast open / EXP-slow close** VCA, a closed floor, and an on/off **enable crossfade**). It
**composes the module kit** (`ChannelLinker` linked key + a Peak `EnvelopeFollower`) and adds the gate-specific
state machine on top. Extracted from OrbitCab's in-amp gate (written portable by design) so plugins don't
reinvent it.

- **feat(dynamics):** `NoiseGate` — a two-phase **keyed** API (`analyse()` fills a per-sample gain curve from
  the clean KEY; `applyGain()` attenuates a possibly-different downstream buffer, so any latency between the two
  is free lookahead) plus a self-keyed `process()` convenience. A `Config` voicing struct (defaults = OrbitCab's
  shipped tuning), `seedEnabled()` for a restored on-state, `currentGain()`/`currentCoreGain()` for a GR meter;
  zero latency.
- **RT / robustness:** NaN/Inf-safe (detector-input clamp + non-finite heal) and denormal-safe in software
  (Law 8); no-allocation-in-`process()` proven; block-split NULL + sample-rate-invariance + adversarial
  fail-open / transient / low-note-chatter tests.

## v0.5.0 — the non-uniform (Gardner) convolver: block-independent, cheaper at small buffers

The v0.4.0 follow-up is delivered. A **non-uniform partitioned (Gardner 1995) convolver** — a time-domain head +
geometrically growing overlap-save FFT stages — replaces the fixed-`P=128` path. It is **block-INDEPENDENT**,
**true sample-zero-latency**, and **cheaper than `juce::dsp::Convolution` across the small, low-latency buffers a
live rig runs** (JUCE only wins the mean at large power-of-two blocks — a theorem of true zero-latency, not a
shortfall). Flat mean **~0.6 %RT on an Apple M5 Pro / ~1.08 % on an Intel i9-13900H** at every DAW buffer from 16
to 4096, from one `prepare()`.

- **feat(convolution):** `NonUniformConvolver<Fft>` — the mono zero-latency NUPC primitive (head `P0` + capped
  octave-doubling to `B_max`, its own frequency-domain delay line per stage).
- **feat(convolution):** `MatrixConvolverNupc<Fft>` — the shipping 2×2 matrix convolver on one raw-L/R history:
  all four routings (mono / LRDiag / MSDiag / Full) + a **click-free 2-slot smoothstep crossfade** for live IR
  swaps. A NULL-verified drop-in for `MatrixConvolver`.
- **feat(lineareq):** the linear- & mixed-phase EQ now convolves on `MatrixConvolverNupc` (A/B transparent — the
  change is CPU + zero latency only).
- **fix(convolution):** removed the long cold-prime crossfade — a cold FDL already yields the exact causal
  convolution, so the ~2.7 s first-activation fade only attenuated correct output and never completed for renders
  shorter than it (offline/short renders came out ~10 dB down). Every swap now uses the short anti-click fade
  (in `NonUniformConvolver`, `MatrixConvolver`, `ConvolutionEngine`).
- **fix(convolution):** `setIr()` broadcasts a mono IR to both channels on a stereo instance (was rejected); a
  `static_assert` pins `state_` lock-free; `maxIrSamples` capped against a stage-offset overflow.
- **docs/tools:** `PERF-NUPC-VS-JUCE.md` (a two-machine 4→8192 fine log-ladder sweep — Apple M5 Pro + Intel i9-13900H,
  adaptive 3–10 warmed reps — rendered to an in-repo SVG chart), `PERF-CONVOLVER-JUCE-GAP.md` (the design ADR);
  `fftbench` head-to-head + `FCORE_FINE_SWEEP`; `tools/plot-convolver-sweep.py`.

Verified by an architecture consilium (per phase), a Popper falsification campaign (differential fuzzer vs
`PartitionedConvolver` + double precision, ASan/UBSan/TSan), and a release-candidate crew review. `ctest` green.

## v0.4.0 — the #1 performance debt resolved

The scalar-FFT / `O(P)` direct-head convolution bottleneck — long-convolution cost that **exploded at host
block 2048+** — is fixed. Cost is now **block-INDEPENDENT**, zero-latency, and JUCE-free.

- **feat(fftpffft):** a new optional, compiled SIMD FFT backend `felitronics::fftpffft::PffftRealFft`
  (`-DFELITRONICS_WITH_PFFFT=ON`, default OFF) — vendored 2-file pffft, hidden-visibility, cross-backend NULL
  parity tested on x86-64/SSE + arm64/NEON incl. ASan/UBSan. (#25)
- **feat(lineareq):** the audio FFT backend is now a template parameter; the design-time FFTs stay pinned to
  the scalar packed-Hermitian layout and the split is **compile-enforced** — a SIMD backend cannot silently
  corrupt a designed FIR. (#24)
- **feat(lineareq):** the convolver partition is **decoupled from the host block** (fixed internal `P=128`).
  A 131072-tap linear-phase EQ is ~2.0 %RT with pffft at every host block (was ~39 %RT @ block 8192 on the
  old scalar path). (#26)
- **feat(convolution):** SIMD-aligned every FFT-seam buffer (`core::fft::SeamAllocator<64>`). (#23)
- **docs:** `PERF-SCALAR-FFT-BOTTLENECK.md` marked RESOLVED.

Known follow-up (tracked separately): the fixed-`P=128` convolver is block-independent but 2–13× more CPU than
`juce::dsp::Convolution` at host blocks ≥ 512; a non-uniform (Gardner) partitioned convolver is planned to
beat it. See `docs/PERF-CONVOLVER-JUCE-GAP.md`.

Prior versions: see the `v0.1.x` – `v0.3.0` git tags.
