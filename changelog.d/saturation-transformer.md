### saturation — Transformer shape

`WaveShaper::Shape` gains `Transformer`, appended as 6 so no existing value moves. It is the one shape with memory: per
channel and per oversampled sample the Saturator runs a flux proxy L (a leaky integrator at 40 Hz, unity DC gain) and
outputs sat(L) + sat'(L)·(x − L) with sat(u) = tanh(k·u)/k — the exact derivative of the saturated flux, with no
differentiator filter. The low end saturates and the top passes: at driveDb 6 on a full-scale sine (48 kHz, os 4) THD
is 8.95 % at 20 Hz, 1.40 % at 160 Hz, 0.096 % at 640 Hz and 0.0016 % at 5 kHz, where Tanh has 6.6 %; the same 80 Hz
THD at 44.1, 48 and 96 kHz to 0.002 dB. It has no knob of its own (drive, mix, output and autoComp; `Params::bias` is
not read), no hysteresis, no DC blocker and no added latency. The curve is SLOPE-normalised (norm = 1/k, so
`slopeAtZero()` is exactly 1.0f and drive-compensation is exactly 1.0f): a -60 dBFS multitone comes out within 1.2e-10
of Tanh at driveDb 0, and the model's output peak never exceeds its input's. At os 1 the stage nulls against a
double-precision evaluation of the model to 9.3e-7. The model state is one float per channel
(`Saturator::kModelFloats`), in `Storage` and `Storage::bytes()` — which grow by exactly that — zeroed by prepare() and
reset(), flushed per sample, and gated like the DC blocker, so a channel that leaves and returns, or a switch to another
shape and back, starts it from zero. Tanh, Atan, Cubic, Asym, Tube and Transistor are unchanged bit for bit, settled
and gliding.
