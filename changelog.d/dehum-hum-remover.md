### dehum — a new module: `HumRemover` removes mains hum (50/60 Hz and its harmonics) from an input

`felitronics::dehum::HumRemover<MaxHarmonics = 64>` — header-only, depends on `core` only, zero latency,
no heap. A phase-coherent synchronous canceller: one complex weight per harmonic per channel, adapted by
LMS against references from one uint32 mains phase. With a fixed step it is exactly a notch comb, and the
step is set from a width in hertz, so every harmonic gets the same absolute width: 0.5 Hz while playing
(a partial 0.5 Hz from a harmonic loses 0.97 dB, 3 Hz away 0.03 dB), 3 Hz while the input is quiet.

- **Modes** Off / 50 / 60 / Auto. Auto decides from Goertzel bins {50, 100} vs {60, 120} on 0.5 s windows
  (summed over channels as powers): a family must win by 6 dB in three consecutive quiet windows; ambiguous
  → it stays dry; decided only in quiet, so it is sticky while playing. Hum-only input decides in ≈ 1.0 s,
  missing-fundamental, odd-only, even-only and opposite-phase stereo hum included. Moving off a decided family
  needs that family to have lost 10 dB; a real swap re-decides in ≈ 1.3 s.
- **Frequency**: a coarse estimate from the bins' phase advance, the median of up to four in a row, acquired only
  from a clean set (two or more within 0.05 Hz), so a knock or a palm on the strings does not throw it off. Then an
  FLL from the weights' own rotation (K·τ = ½, slewed ≤ 0.05 Hz/s, bounded by `lockRangeHz`, default ±1 Hz). It runs
  only in quiet, while the coarse estimate agrees and the canceller explains the input — broadband, or, where hiss
  makes that impossible (a guitar DI), in the weights' own narrowband terms. Detuned mains (49.6 … 60.4 Hz) lock
  within 0.01 Hz in ≈ 1.2 s; a held note beside a harmonic does not move it; a single-coil DI locks ≈ 1–2.5 s after
  Auto decides. While playing the frequency holds, and the first quiet gap re-locks a drift (0.8–1.05 s for 0.06 Hz).
- **Quiet / active**: quiet needs (a) the window to be hum — broadband within 12 dB of the hum bins, or the
  fundamental a clean line at the mains bin (a humbucker's hum is only ≈ 4 dB over its hiss), which serves the
  tracker and lends the learn width only near the broadband rule; (b) the residual within 12 dB of its own floor
  (a 4 s minimum rising ≤ 2 dB/s — a note is not quiet however stationary; a louder hum is told apart and
  re-qualifies in ≈ 0.5 s); and (c) no per-sample trigger in the last 50 ms (the residual catches a pluck's first
  samples). `active()` is the 64-sample power 12 dB over the floor for 20 ms.
- **The ghost a note leaves is returned**: a note beside a harmonic leaves its beat in that harmonic's weight, which
  would ring out at the playing width (1.56 s to −60 dB) — behind high gain an after-ring after every chug. An event
  opens on a trigger with a snapshot of the weights before the attack; when the input is back near its floor, a note
  of up to 3 s returns the weights to the snapshot over 40 ms. Every return is checked over the next 16 periods
  against the residual the weights it replaced would have made on the same samples, and taken back if those did
  better. No return when the hum dropped under the note (hands on the strings) or the frequency moved. Measured:
  the 330 Hz ghost under −60 dB re the tone 0.07 s after the stop; palm-muted chugs over a single coil leave
  ≤ −62 dB re the chug 50–90 ms after each end (it was −37).
- **LEARN**: `beginLearn()` / `learnStatus()` / `endLearn(commit)` → `LearnReport` { mains, hz, humDb,
  harmonicDb[], floorAfterDb (peak |out| over the settled window), settledAtMs, verdict, kept }. Disturbed
  (a trigger, or a floor window that rises) freezes adaptation from that sample; a floor window under
  0.25 s is TooShort; `kept` says whether the stage now runs the learned model; anything but a committed,
  Learned learn restores the snapshot — weights, frequency, decision and its own phase frame — bit for bit
  through a 20 ms ramp. `endLearn` with no learn running says NotLearning.
- **Telemetry**: `detected()`, `mainsHz()`, `oscillatorHz()`, `locked()` (judged on a whole hop, kept while the
  FLL is held), `humDb()`, `harmonicDb(k)`, `quiet()`, `active()`.
- **Laws**: bit-identical under re-slicing (8a), the weights flushed and healed on the grid (8), a pause is
  silence for every shared clock (11c), a stopped channel keeps its model (11a), `prepare()` disarms and
  refuses (11b), and a restart is a fresh instance bit for bit; the weight update cannot contract under any
  `-ffp-contract` mode (10); digital black in is bit-exact black out, and the correction fades out over the
  64 zeros and back in after them, so region ends and starts do not step. Harmonics stop where the top of
  the lock range would cross maxHz or 0.9·Nyquist. setParams / beginLearn / endLearn belong on the process
  thread or strictly between process calls.
- **Bench**: `tools/dehumbench.cpp` — ≈ 4.0 µs per 64-frame block, 0.30 % of a core at 48 kHz stereo,
  24 harmonics (Apple M5 Pro). Footprint 39 KiB at 24 harmonics, 104 KiB at 64, plus ≈ 9 KiB.
- **Limits, stated in the header**: a slow attack (≈ 10 ms) within ≈ 100 ms of the previous note may trip no
  trigger, and that note keeps the old ghost; a note held longer than 3 s gets no return, so a partial of it
  inside a harmonic's notch rings out at the playing width — the A string is the plain case (5·110 = 11·50,
  6·110 = 11·60).
- **Cannot**: remove hum made after its input (e.g. inside a physical amplifier downstream).
