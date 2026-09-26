### io — the WAV writer quantizes on `dither`'s grid, so a dithered master is written as the codes it was given

`writeWav`/`writeWavMemory` scaled 16- and 24-bit PCM by 2^(bits-1) − 1 (32767, 8388607) under `llround`, while
`readWav` divides by 2^(bits-1) and `dither::Dither` defines a code as `floor(v·2^(bits-1) + 0.5)` clamped to
[−2^(bits-1), 2^(bits-1) − 1]. A master Dither had already put on the grid was therefore quantized a second time,
onto a grid one LSB narrower, with no dither. Measured before the fix, writing every code k as k/2^(bits-1): 32767
of the 65536 16-bit codes came back moved by one LSB toward zero (every |k| > 16384 — 20000 was written as 19999,
−32768 as −32767), 8388607 of the 16777216 24-bit codes, and 30001 of 48000 samples of a Weighted-dithered 0.9
sine at 16 bit (30002 at 24 bit).

The writer now uses Dither's rule exactly — the same scale, mid-tread round-half-up on both signs, the clamp in
the double domain before the conversion, NaN/Inf → 0 as before. The reader was already on the grid and is
unchanged. So write∘read is the identity on every code, read∘write gives a canonical image's bytes back, and a
Dither output is written as the codes Dither chose. What a caller sees change: −1.0 is now the bottom code
(−32768, was −32767); +1.0 is still the top code; an off-grid sample past half scale can land one LSB away from
where it used to; a sample exactly on a half-code rounds up on both signs, where `llround` took a negative one
away from zero. Images written before are not byte-identical to images written now; the 32-bit float writer is
untouched.

New suite `felitronics_io_grid_tests` (59 checks, links `dither` to test the agreement; the io library stays
zero-dep): every 16-bit code and a 24-bit sample dense where the two grids part, both directions; Dither's output
at 16 and 24 bit under all three shapings; ties, the ends, ±DBL_MAX, NaN, ±Inf and subnormals against Dither's rule.
Planted failures run through the same instruments: the old writer moves exactly the codes with |k| > 2^(bits-2),
and a writer on the right grid with `llround`'s rounding misses exactly the negative ties.
