### dehum — a new module: `HumRemover` removes mains hum (50/60 Hz and its harmonics) from an input

`felitronics::dehum::HumRemover<MaxHarmonics = 64>` — header-only, depends on `core` only, zero latency,
no heap. A phase-coherent synchronous canceller: one complex weight per harmonic per channel, adapted by
LMS against references from one uint32 mains phase. With a fixed step it is exactly a notch comb, and the
step is set from a width in hertz, so every harmonic gets the same absolute width: 0.5 Hz while playing
(a partial 0.5 Hz from a harmonic loses 0.97 dB, 3 Hz away 0.03 dB), 3 Hz while the input is quiet.

- **Modes** Off / 50 / 60 / Auto. Auto decides from Goertzel bins {50, 100} vs {60, 120} on 0.5 s windows
  (summed over channels as powers): a family must win by 6 dB in three consecutive quiet windows; ambiguous
  → it stays dry; decided only in quiet, so it is sticky while playing. Hum-only input decides in ≈ 1.0 s,
  missing-fundamental, odd-only, even-only and opposite-phase stereo hum included.
- **Frequency**: a coarse estimate from the bins' phase advance, then an FLL from the weights' own rotation
  (gain K·τ = ½, slewed ≤ 0.05 Hz/s, bounded by `lockRangeHz`, default ±1 Hz; it sums only the channels that
  adapted, and moves the increment only on the grid). It runs only while the window is hum and the coarse
  estimate agrees with it:
  detuned mains (49.6 … 60.4 Hz) lock within 0.01 Hz in ≈ 1.2 s; a held note beside a harmonic does not
  move it. While playing the frequency holds — the header gives the arithmetic for a drift without a quiet gap.
- **Quiet / active**: quiet needs the window to be hum (broadband within 12 dB of the hum bins), the
  canceller's residual within 12 dB of its own floor (a minimum over 4 s that rises ≤ 2 dB/s — a note is
  not quiet however stationary, a louder hum is told apart and re-qualifies in ≈ 0.5 s), and no trigger in
  the last 50 ms (the residual's per-sample trigger catches a pluck's first samples). Note safety wins over
  a fast clean-up: the ghost a note leaves decays at the playing width (−60 dB re the note in ≈ 1.6 s).
  Auto moves off a decided family only when that family has lost 10 dB. `active()` is the 64-sample power
  12 dB over the floor for 20 ms.
- **LEARN**: `beginLearn()` / `learnStatus()` / `endLearn(commit)` → `LearnReport` { mains, hz, humDb,
  harmonicDb[], floorAfterDb (peak |out| over the settled window), settledAtMs, verdict, kept }. Disturbed
  (a trigger, or a floor window that rises) freezes adaptation from that sample; a floor window under
  0.25 s is TooShort; `kept` says whether the stage now runs the learned model; anything but a committed,
  Learned learn restores the snapshot — weights, frequency, decision and its own phase frame — bit for bit
  through a 20 ms ramp. `endLearn` with no learn running says NotLearning.
- **Telemetry**: `detected()`, `mainsHz()`, `oscillatorHz()`, `locked()`, `humDb()`, `harmonicDb(k)`,
  `quiet()`, `active()`.
- **Laws**: bit-identical under re-slicing (8a), the weights flushed and healed on the grid (8), a pause is
  silence for every shared clock (11c), a stopped channel keeps its model (11a), `prepare()` disarms and
  refuses (11b), and a restart is a fresh instance bit for bit; the weight update cannot contract under any
  `-ffp-contract` mode (10); digital black in is bit-exact black out, and the correction fades out over the
  64 zeros and back in after them, so region ends and starts do not step. Harmonics stop where the top of
  the lock range would cross maxHz or 0.9·Nyquist. setParams / beginLearn / endLearn belong on the process
  thread or strictly between process calls.
- **Bench**: `tools/dehumbench.cpp` — ≈ 3.4 µs per 64-frame block, 0.26 % of a core at 48 kHz stereo,
  24 harmonics.
- **Cannot**: remove hum made after its input (e.g. inside a physical amplifier downstream).
