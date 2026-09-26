### stereo — MonoBass: the corners glide, `enabled` fades, the first write snaps, and setAir re-tunes on a live move

Measured through the mastering chain on a 103.7 Hz tone with L and R 90° apart (a large Side), max|Δ²y| where the
change reaches the output against -84.2 dBFS for the steady tone: the crossover corner 60 -> 250 Hz stepped at
-45.3 dBFS (its four SVFs were redesigned at once) and 60 -> 70 Hz at -69.3; `enabled` false -> true read -48.7 and
true -> false -21.5 (a hard switch with a crossover reset). The air shelf's corner and `enabled` switched as hard,
which the Δ² metric cannot see on a 5 kHz tone: through a Q-2 notch at the tone, which it can, corner 6000 -> 3000 Hz
read -25.9 and `enabled` true -> false -26.3, where the plateau's own 20 ms glide reads -66.3.

* THE CROSSOVER CORNER rides a one-pole `core::Smoother` (`kFreqSmoothMs` = 30 ms, the EQ's time constant), advanced
  by exactly 64 samples at every `core::StateGrid` boundary with the SVFs redesigned there when it moved — EqBand's
  idiom. A stretch the island skips (bypassed, full-wide, a mono or a clock-only call) still crosses its boundaries:
  the corner walks through it as through audio. Now -73.0 and -83.0.
* `enabled` IS THE xf CROSSFADE: off fades the Side to dry as a full-wide lowWidth does (20 ms), and the bit-exact
  passthrough follows once it settles; on fades back in from a crossover restarted at zero. Now -75.5 and -74.3. It
  used to be a hard toggle for parity with StereoWidth; the mastering chain's live preview is the product that
  disagrees, and a host that wants a hard bypass has one — not calling the stage.
* THE AIR SHELF'S CORNER rides a Smoother of its own on the grid, and its `enabled` glides the plateau to and from
  0 dB (where the shelf is skipped). Notch reading now -56.9 for the corner and -66.4 / -66.3 for the enable edges.
* EVERY WRITE BEFORE THE FIRST SAMPLE after prepare()/reset() SNAPS — lowWidth and the air plateau included, which
  used to glide from the prepared values when written between prepare() and the first sample. So "set, then prepare"
  and "prepare, then set" rendered differently (through the mastering chain, 71 529 and 59 030 samples of a
  one-second render for a lowWidth and a plateau write); they are one behaviour now, the set-then-prepare one, whose
  bits are unchanged. `reset()` lands both corner glides and designs the crossover at its corner.
* setAir RE-TUNES ITS WIDTH BAND ON A LIVE CORNER MOVE. `wasFreq` was captured after the assignment, so the band was
  re-tuned only when the clamp moved the corner — contrary to the note at `retuneWidthBand()` — and a live 6000 ->
  8000 Hz move kept summing across both bands: `airJudgedSamples()` read 2000 where the new band had judged 1000.

THE REVIEW ROUND (codex astra) found five defects in the first spelling of this, all fixed: the linear ramps (width,
crossfade, plateau) froze through a skipped stretch while the corners walked on — a disable followed by a pause ran
its whole fade on return — and now spend it sample by sample; with the bass already idle an air fade-out that landed
mid-call never retired the island, so the rest of the call stayed in the M/S round trip and depended on the cut
(7 471 floats between whole-block and one-sample renders) — the island now retires on the air's settle too, on the
sample clock; that exit did not reset the shelf, which then replayed a 1.42e-5 tail on the next enable; the width
ramp froze while a disabled bass let the air keep the island running; and the air's band measurement stopped the
moment `enabled` went false although the shelf was still audibly fading — it covers the fade now (an enabled 0 dB
plateau is measured as before). Its test critique was right too: a noisy programme's own Δ² hid a hard step, and
nothing observed the skipped time, so the click checks run on pure tones against an absolute -60 dBFS bound with a
spliced hard switch as the precondition, and new accessors (`crossoverDesignHz()`, `airDesignHz()`, `crossfade()`,
`lowWidthNow()`) let a test see where the glides stand after a pause, a mono stretch and an idle island.

New suite `felitronics_monobass_glide_tests`: the two configuration orders bit-identical for seven parameters, the
snap after a reset() mid-glide, a timeline of writes (both corners with a retarget, both enables, lowWidth, the
plateau) bit-identical under per-sample, 64, 100, 4096 and ragged-with-empty-calls cuts with a clock-only pause and a
mono stretch in it, a corner glide started while the island is idle, the fades' passthrough and edges, the setAir
re-tune, the skipped-time positions, the air retiring the island, no allocation. Mutations caught: `wasFreq` after the
assignment, an extra corner tick per call, lowWidth gliding after a reset, no corner ticks in a skip, no linear-ramp
steps in a skip, an instantaneous crossover, no air retire, a width frozen while the bass is idle.
