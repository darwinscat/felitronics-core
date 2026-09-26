### dynamics · dynamiceq — the Compressor's makeup glides; LaneDynamics can release a switched-off point instead of snapping it

**Compressor makeup.** `makeupDb`, and `autoMakeup` with every curve change that moves it, used to be added as is, so a
live change stepped the output: through the mastering chain on a -12 dBFS sine at 4:1 (steady -76.8 dBFS max|Δ²y|),
makeup 0 -> +3 read -30.5 dBFS, autoMakeup switched on -21.4, and a -18 -> -30 dB threshold move with autoMakeup on
-19.5 (the curve's own change rides the GR ballistics; its auto-makeup did not). The makeup now glides linearly in dB
over `Compressor::kMakeupGlideMs` = 30 ms, one step per processed sample (double accumulator, landing on the target
itself), and a clock-only call spends the glide sample by sample as the silence it stands for would: -73.8, -69.4 and
-67.7. Every write before the first sample after prepare()/reset() snaps, so a compressor configured before its first
sample renders what it always did. New `felitronics_dynamics_makeup_glide_tests`: both configuration orders agree, the
glide replayed sample for sample below the threshold (where the output IS the delayed input times the makeup), no step
from makeup or autoMakeup (a threshold move is held to the same move's edge without autoMakeup, scaled by level), law
8a under per-sample, 64 and ragged cuts with a clock-only pause mid-glide, a glide spanned by a pause landed on return,
repeated writes free, no allocation. Mutations caught: no glide, a pause that does not spend it.

**LaneDynamics — RELEASE ON DISENGAGE, opt-in** (`setReleaseOnDisengage`, off by default, so TabbyEQ and every other
consumer is bit-identical). A point whose dynamics are switched off mid-duck — dyn.on false, rangeDb 0, or no sidechain
— used to disengage on the edge: every lane's delta to 0 at once, a step in the band's gain (through the mastering
chain, a -9 dB range ducking a -12 dBFS tone: -16.1 dBFS against -82.1 for the steady tone). Opted in, each lane's
gain-reduction follower is driven to 0 dB through its own release ballistics on the same control grid, and the point
disengages exactly as before once every |delta| is under `kReleaseFloorDb` (1e-6 dB); a clock-only pause spends the
release like silence, and switching back on mid-release carries on from where it stands. EqBand ignores deltas while
its own `dyn.on` is false, so for as long as it releases this layer writes the band back with `dyn.on` true and
restores the caller's value when it is done — the delta bell at 0 dB is transparent, so from that sample the band
renders bit for bit what a band that never had dynamics renders (pinned). Now -55.3 dBFS: 39 dB down, and what is left
is the band's 16-sample control grid, the same zipper its own attack and release have while engaged. It is opt-in
because it writes the band: TabbyEQ re-writes every band every block (so the flag would flip back and forth, a
redesign each time) and stops calling processBand() for a point once no point is dynamic (so the release would never
finish) — both would have to change before it can opt in. New `felitronics_dynamiceq_release_tests`.
