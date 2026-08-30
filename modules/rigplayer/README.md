# felitronics::rigplayer — the pack player, and what a host owes it

One device of a `.orbitrig` pack, playing. Hand it the pack's structures (`namz::rig::Rig`, as
`namz::rig::loadRigManifest` returns them) and a way to get a file's bytes, tell it where the knobs
stand, feed it audio. Header-only, JUCE-free, no thread of its own. Moved here from OrbitCapture NAM
so the capture app and a plugin play through ONE player; the app is the reference host.

Signal order, per block:

    in → dry copy → pre tone (bands, curve) → chain trim → ┬→ trim A → model A → align ─┐
                                                           └→ trim B → model B → align ─┴→ mix
       → post tone (bands, curve) → dry/wet blend → out

- WHICH files sound is `namz::rig`'s policy (a turned control is law) joined to `pickBlend` along
  the gain dial — `RigSelection.h`, pure.
- WHEN a model may be replaced, and how loudly each is heard, is `felitronics::nam::BlendLaw` — run
  once per block on the audio thread, the only writer of the weight.
- A tone knob plays as the pack describes it: `sections` become biquads, `positions` (a curve)
  become one minimum-phase FIR per side. A blend knob mixes the DI through the dry path's response.
- The models' offsets come from the pack (`lag_samples`); a pack without them can be measured by the
  host (`AlignmentTable`) and handed in.

## Threads — the whole contract

| call | thread |
|---|---|
| `prepare(rate, maxBlock, channels ≤ 2)` | message, never while `process()` runs |
| `load(rig, source)` / `unload()` | message |
| `setDial(name, degrees)` / `setSwitch(name, value)` / `setToneOverride` / `setBlendShape` / `setNormalize` / `setAlignment` / `setColdAfterSeconds` | message |
| `service()` then `takeLoadJob()` | message, from a timer — a few times a second at least |
| `RigPlayer::run(job)` | ANY thread but the audio one — a pool, a worker, or right here |
| `deliver(loaded)` | message — strictly; it installs into the stages |
| `process(io, channels, n)` | audio — allocates nothing, locks nothing, touches no file |
| every read-out (`liveMix`, `heldFileId`, `selection`, `knobValue`, the drawing readers…) | message |

**Loading is a job the host runs.** `service()` never loads. When the law wants a model,
`takeLoadJob()` hands out ONE job at a time — slot, `files[].id`, the bytes if already fetched, the
`ModelSource`, the numbers the stages were prepared with. `run(job)` is a pure function of it: the
bytes through the source when the job has none, then the heavy half of a load
(`NamStage::prepareModel` — parsing, two instances, the prewarm; some twenty milliseconds for a
WaveNet). `deliver(loaded)` is the light half. **A job taken MUST come back through `deliver()`** —
with a null model if it failed — or the law waits for it forever. A job that returns after `unload()`
is dropped by the player; the host need not track it. `serviceHere()` runs the same path on one
thread, for a host without a worker and for a test.

**`ModelSource` is called from wherever the host runs `run()`.** Make it safe to call off the message
thread: it is called at most once at a time (one job in flight), but not from the thread that opened
the pack. The reference host reads from a `juce::ZipFile` kept open, one entry per call.

**A slot at rest goes cold.** A dial parked on a capture keeps the neighbouring capture in the other
slot at exactly zero, warm and waiting — and its network used to run every block for nothing (measured
in OrbitAmp's block, prepared for two planes and fed one: 4.5 % of a P-core, 14 % of an E-core). After
`RigPlayer::kColdAfterSeconds` (2 s) at exactly zero under an unchanged request the slot goes COLD: its
model is not run — nor mixed — and it stays loaded: nothing is fetched, nothing freed; its delay
line is cleared as it falls asleep, as a landing clears it. The law owns the flag (`BlendState::cold`)
and wakes the slot on the first block of the next change of request, warm-up first — the same
re-landing a load ends with, with the model's own field — so nothing unfed is ever heard, and the
first turn after a rest trails the hand by one warm-up (some 130 ms for a WaveNet, then the ordinary
slew; a model that declares no field is heard at once). A slot with any weight is never cold: between two captures both models run, on one
they do not. `setColdAfterSeconds(s)` sets the rest (zero or less = never); `slotCold(i)` and
`coldBlocks(i)` (since `clearCounters()`) read it out beside `warmBlocks()`, for a dump or a badge.

## What a host must do

1. `prepare()` in the host's own prepare; re-prepare on a rate or block change (the filters are
   designed for a rate; a model built for another rate is prepared again on install).
2. `load()` once per pack/device; the player opens on the pack's defaults (`namz::rig::defaultSettings`)
   and, if that combination was never captured, on the closest one that was.
3. A timer: `service()`, then `while (auto job = takeLoadJob()) …` — run it where you like, bring it
   back with `deliver()` on the message thread.
4. Report `latencySamples()` to the host after every `deliver()` and `prepare()` — it is the models'
   rate-matching; the alignment delays are relative and the reference model carries none.
5. Knobs by NAME, the pack's names: `setDial(name, degrees)` for anything with a sweep (a captured
   dial, a tone knob, a blend knob), `setSwitch(name, value)` for a token. `knobValue(name)` reads
   any of them back as the pack spells it; `dialDegrees()` is the crossfade dial's angle,
   `settings()` the captured combination.

## Drawing — the one honest source of the curve

The tone as it PLAYS, not as it was measured: `commonGrid()` (1/12-octave, 20 Hz – 20 kHz) with
`curveDb(side)` — the summed curve-form knobs of a side, empty when none — and `curveActive(side)`
(a FIR, or a wire); `bands(side)` — the biquads of the section-form knobs at their current
positions; `blendGains()` — the dry and wet gains for the blend knob's position. A host that draws
from these draws what is sounding; recomputing from the pack's ladders draws something else.

## Not the player's business

The source (a loop, the jack), a cabinet IR, output routing, a pre-model trim the bench owns, level
metering, presets — the chain AROUND the player. In OrbitCapture NAM that is `AuditionEngine`
(the JUCE audio side) and `AuditionController` (the message-thread wiring, the pool that runs jobs).
