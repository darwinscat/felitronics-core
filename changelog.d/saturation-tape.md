### saturation — Tape shape

`WaveShaper::Shape` gains `Tape`, appended as 7 so no existing value moves. The Saturator runs it as a model: per channel
and per oversampled sample, a pre-emphasis E(s) = (1 + s/ω1)/(1 + s/ω2), Tanh's core tanh(k·e)/tanh(k), then the exact
inverse D = 1/E, with f1 = 1/(2π·50 µs) = 3183.1 Hz (the NAB 15 ips corner) and f2 = f1·10^(6/20) = 6351.1 Hz, a +6 dB
high shelf. The top reaches the core louder and saturates first: at driveDb 6 (48 kHz, os 4, autoComp 0) the
fundamental of a 10 kHz sine is 0.5 dB compressed at -11.104 dBFS and that of a 1 kHz sine at -6.429 — 4.675 dB apart,
where the analog |E(10k)|/|E(1k)| is 4.645 dB — while Tanh compresses both at -6.126. The 10 kHz point stays within
0.0285 dB across 44.1, 48 and 96 kHz. Since D·E = 1, a -60 dBFS multitone comes out within 6.4e-10 of Tanh at the same
drive, a full-scale one within 4.17e-7 at driveDb 0, and drive-compensation is Tanh's. Each section is a first-order
bilinear transform with its own corner prewarped (both exact at every rate), in transposed direct form II. On
full-scale 1/5/10 kHz sines, a 1 kHz square and band-limited clicks at driveDb 0/6/12 (os 4, autoComp 1), Tape's peak
is Tanh's or lower to the printed 0.001 dB. At os 1 the stage nulls against a double-precision evaluation of E, the
core and D to 2.43e-7. Where f2 >= 0.45·fsOs (os 1 below 14113.6 Hz) the pair is bypassed and Tape renders Tanh's bits.
Tape has no knob of its own (`Params::bias` is not read), no head bump, no HF roll-off, no hysteresis, no DC blocker
and no added latency; a standalone `WaveShaper` set to Tape is Tanh's curve, operand for operand.

`Saturator::kModelFloats` rises from 1 to 2 (the Transformer keeps slot 0, Tape uses slots 0 and 1), so `Storage` and
`Storage::bytes()` grow by one float per channel. The state is zeroed by prepare() and reset(), flushed per sample and
gated like the Transformer's. A direct switch between the Transformer and Tape starts the new model from zeroed state:
at os 1 the switched stage renders what a fresh stage of the new shape renders from the switch on. A drive glide lands
on the settled stage's bits 20 samples after the landing at os 1 and 9 at os 4, once the de-emphasis state has
forgotten the glide. Tanh, Atan, Cubic, Asym, Tube, Transistor and Transformer are unchanged bit for bit, settled and
gliding.
