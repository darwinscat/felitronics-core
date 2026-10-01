<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# codecgrid — cross-language NULL (numpy / scipy)

A **dev-only** cross-implementation check, deliberately **not** wired into `ctest` — core stays Python-free, and
CI correctness is covered by the in-tree suites, which null every transform against its definition written out
directly. This harness is the other half: the same numbers recomputed by a library that shares no code with the
module.

`dump.cpp` makes four programmes that carry a known frame grid (`../SyntheticCodec.h`: a transform, the smallest
30 % of every frame's coefficients set to zero, the inverse — no codec and no audio file), scans each with
`GridScan`, and writes the PCM, the curve, the reading and the zero share; it also runs a stretch through
`BackResampler` both ways. `check.py` recomputes, from the PCM alone:

- the **MP3** hybrid lines from the standard's formulas in numpy (the window table is read out of
  `Mp3Window.h`, the only thing shared);
- the **AAC** and **CELT** MDCT through `scipy.fft.dct` (type IV) of the folded frame, with the windows from
  their formulas (`numpy.kaiser` for the Kaiser-Bessel-derived one);
- the **curve** in double precision with a logarithm per coefficient (the module multiplies in single
  precision and takes one per sixteen frames);
- the **reading** with `scipy.ndimage.median_filter (mode='mirror')` and `numpy.median`;
- the **resampler** with `scipy.signal.resample_poly (window=('kaiser', 14.0))`.

## Run

```sh
python3 -m venv venv && ./venv/bin/pip install numpy scipy
c++ -std=c++20 -O2 -I ../../../core/include -I ../../include dump.cpp -o dump
mkdir -p out && ./dump out
./venv/bin/python check.py out ../../include/felitronics/codecgrid/Mp3Window.h
```

## Last result (2026-10-01, numpy 2.5.3, scipy 1.18.1, Apple clang, arm64)

```
PASS mp3 curve: max |difference| 6.89e-05 dB over 576 x 32 cells
PASS mp3 reading: offset 133 (coded at 133), score 849.700 vs 849.700, second 3.059 vs 3.059, dip 469.714 vs 469.714, cells 32 vs 32
PASS mp3 zero share: max |difference| 1.46e-08 over 32 bands, mean share 0.303
PASS aac_sine curve: max |difference| 4.03e-05 dB over 1024 x 32 cells
PASS aac_sine reading: offset 480 (coded at 480), score 922.649 vs 922.649, second 3.111 vs 3.111, dip 475.633 vs 475.633, cells 32 vs 32
PASS aac_sine zero share: max |difference| 1.28e-08 over 32 bands, mean share 0.304
PASS aac_kbd curve: max |difference| 3.88e-05 dB over 1024 x 32 cells
PASS aac_kbd reading: offset 1001 (coded at 1001), score 922.084 vs 922.084, second 3.809 vs 3.809, dip 473.149 vs 473.149, cells 32 vs 32
PASS aac_kbd zero share: max |difference| 1.28e-08 over 32 bands, mean share 0.304
PASS celt curve: max |difference| 2.25e-05 dB over 960 x 32 cells
PASS celt reading: offset 548 (coded at 548), score 816.901 vs 816.901, second 4.603 vs 4.603, dip 493.967 vs 493.967, cells 32 vs 32
PASS celt zero share: max |difference| 1.76e-08 over 32 bands, mean share 0.318
PASS resample 160:147: length 65307 vs 65307, max |difference| 9.53e-08 (peak 0.532)
PASS resample 147:160: length 55125 vs 55125, max |difference| 1.08e-07 (peak 0.491)
ALL NULLS PASS
```

The curves agree to the single precision of the module's transform; the readings agree to the digits printed.
`out/`, `dump` and `venv/` are not tracked.
