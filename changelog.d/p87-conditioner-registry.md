<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

### nam — the ledger answers for the WHOLE model, conditioner and heads included

`felitronics::nam::detail` is one registry with three questions — how far a capture's memory reaches,
whether anything in it is recurrent, and whether anything in it is charged NAM's partitioned-FFT ring —
and **two readers spend its answers**: law 11a's drain (a lane the host stops handing over is fed silence
for that long) and P47's stream restart. It walked `config.submodels` and stopped there, so a capture
whose **CONDITIONER is a whole model of its own** (`config.condition_dsp`, which NAM builds with `get_dsp`
like any other model) hid that model's memory from both. The two were short by exactly the same amount,
which is what said it was one defect in one ledger rather than two: **0.905147969723** left after a full
drain against **0.905148267746** after a restart, on the same capture.

- **The composition is a SUM, and that is the number rather than the wording.** The conditioner's output
  is the network's conditioning input, so the two memories are in SERIES — which is also how NAM's own
  arithmetic composes them. Measured from outside the registry, on a loaded two-layer stack with a
  2500-sample Linear conditioner: the impulse's last non-zero sample is **5000**, where the same stack
  without a conditioner reaches 2501 and a worst-of would have answered 2502.
  What the registry promises is an **upper bound**, and the slack is named: the condition enters a layer
  AFTER that layer's own convolution, so the true reach is `Mₒ + H + max(0, M_c − L₀)` and the answer is
  `Mₒ + H + M_c + 1` — over by `L₀ + 1`, never short.
- **All three questions walk the branch, not one of them.** A `Linear` conditioner is charged the
  2048-sample ring (its instance owns the same engine, and its ring feeds the layer arrays through the
  mixin); an `LSTM` conditioner makes the whole capture recurrent, because the cell's state enters every
  layer through a memoryless 1×1 and no finite silence empties it. A mutation stand confirmed each of the
  four one-function-only variants goes red, as do "walk it but discard the field" and "worst-of instead
  of the sum".
- **And two more branches of the same config were not being read**, both of which the sum needs in order
  to be an upper bound at all. A layer array ends in a causal head rechannel whose kernel is the layer's
  own `head.kernel_size`, worth `kernel − 1` samples — `example_models/A2.nam` spells it 16 — and a
  post-stack `config.head` is worth `Σ(kᵢ − 1)`. A `ConvNet` keeps its whole stack in a TOP-LEVEL
  `dilations` array and answered **zero**. Each is normally hidden by NAM's own answer, which `NamStage`
  raises this number with; none of them is hidden behind a **slimmable** WaveNet, which answers zero for
  everything. Measured on ones that load: a slimmable WaveNet with a 16-tap head reads 2 against an
  impulse reaching 16, and one with a ConvNet conditioner reads 2 against an impulse reaching 15.
- **The area of the change, bit-for-bit.** Rendered through both readers at 44.1 / 48 / 96 kHz, NAM's nine
  shipped example captures are **byte-identical** except the two that carry a conditioner, whose drain
  grows by **one sample** (`wavenet_a2_max` 31 → 32, `wavenet_condition_dsp` 45 → 46).
- **The answer is monotone, and it took a review round to make that true.** The first cut of this change
  added the ConvNet reader as a link in a first-non-zero CHAIN, ahead of the declared-field reader — and
  a chain can take a number away. A `Linear` capture carrying a stray `dilations` array loads (its parser
  reads neither key) and the chain answered **2** for a 4999-sample impulse response that answers 4999
  without the stray key: a lane draining 2050 where it needs 4999, i.e. a regression introduced by the
  fix. The three sources are now the stack, then the WORST of the other two. In the same round the
  container branch stopped RETURNING: NAM dispatches on the `architecture` string and never on shape, so
  a `"WaveNet"` carrying a stray `submodels` array loads with its whole stack behind that return
  (measured: 1 answered for a model reaching 4200), and the return also made `isRecurrent` charge a
  conditioner the field was ignoring — the very divergence this change exists to close. Both are gated.
  With those two, every term is added and none is replaced, so no capture can drain shorter than before.
  The one named exception is a number that is not representable at all: a field spelled `1e300` is now
  REFUSED rather than cast, because the cast is undefined behaviour.
- **A dead key must never silence a live reader, and a diverse-testing round found two more of those.**
  NAM's ConvNet parser does not read `layers`, so a ConvNet carrying `"layers": []` had its whole stack
  suppressed — behind a slimmable capture's conditioner, where NAM answers zero, that drained 102 for a
  model reaching 131 and handed back 29 samples. And a boolean is a number to `get<int>()` but not to
  `is_number()`: NAM loads `"dilations":[true,true]` as `[1,1]` and answers 3 for it, where this file
  answered nothing. Both gated.
- **A value that does not fit an `int` is now REFUSED rather than clamped**, because the model NAM built
  does not contain it: `"dilations":[4294967396]` builds a network whose dilation is 100, and clamping to
  INT_MAX made `reset()` spend 2 147 483 646 samples — a measured **23.7 s of synchronous audio-thread
  work** — for a capture that remembers a hundred. Refused reads as absent, and NAM's own answer covers
  what it did build.
- **The gates were checked by mutation, twice, in an isolated copy of the tree:** 21 variants of this
  file, 20 red. The one survivor is an `is_object()` guard proved equivalent (nlohmann's `contains()` is
  false for every non-object, without throwing). A sweep of 3200 generated loadable configs found no
  answer short of the model's measured impulse reach outside law 11a's named recurrent exception.
- **Not fixed here, registered:** the hybrid slimmable-wrapper shape, where NAM dispatches on a top-level
  `config.layers[i].slimmable` marker and the real config then hides under `config.model`. It loads, NAM
  answers 0 and the registry answers 0 however long its stack is (measured 314 samples of leak on one).
  Reading it means restating NAM's dispatch heuristic, which is its own decision and its own number.
  *(Closed by P92 without that restatement — see `p92-unplaced-shape-ceiling.md`.)*
- **And one that is not this module's at all, registered with its reproducer:** a `condition_dsp` that is
  a `SlimmableContainer` with a WaveNet submodel, or a slimmable WaveNet, **crashes the host on load** —
  a null write inside NAM's own `WaveNet::_set_condition_array` during `DSP::prewarm()`, because neither
  class overrides `SetMaxBufferSize` and the conditioner path is the one place that is the only call they
  get. Reproduces on pure `nam::get_dsp` + `DSP::Reset` with none of this code in the path, and
  `prepareModel`'s catch-all cannot catch a SIGSEGV.
