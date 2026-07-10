# measurement — cross-language NULL (numpy)

A **dev-only** cross-implementation check: the C++ measurement math dumped to raw f64, then
NULLed against an **independent numpy recompute**. This is deliberately **NOT** wired into `ctest`
— core stays Python-free, and CI correctness is covered by the in-tree oracle + adversarial C++
suites. This harness exists so the "prove it against a standard library" claim is reproducible.

Each reference is computed from scratch, not a re-run of the same code:
- `convolve()` (our radix-2 **FFT** convolution) vs `numpy.convolve` — a **direct time-domain**
  O(n²) convolution that shares no code path with our FFT.
- `magSpectrum()` (our radix-2 FFT) vs `numpy.fft`.
- `makeSweep()` signal vs the analytic ESS formula recomputed in numpy from the *sanitized* spec.
- `sweep ⊛ inverse ≈ unit-peak δ` via numpy's independent MAC of the C++-built inverse.

## Run

```sh
python3 -m venv venv && ./venv/bin/pip install numpy
c++ -std=c++20 -O2 \
    -I ../../../core/include -I ../../include \
    dump.cpp -o dump
./dump .
./venv/bin/python check.py .
```

## Last result (2026-07-10, numpy 2.5.1)

```
PASS convolve[0] (n=500x300):   max|Δ|=6.004e-13  rel=2.052e-14
PASS convolve[1] (n=2000x2000): max|Δ|=1.631e-12  rel=5.574e-14
PASS magSpectrum (nfft=2048):   max|Δ|=4.280e-12  rel=8.877e-15
PASS makeSweep signal (N=3200): max|Δ|=0.000e+00  rel=0.000e+00   (bit-identical)
PASS makeSweep tail silence:    max|Δ|=0.000e+00
PASS sweep⊛inverse peak=1.0000, floor/peak=1.67e-06 @ ±2000 from δ
ALL NULLS PASS
```

Machine-epsilon agreement against a wholly independent implementation; the sweep matches the
analytic recompute bit-for-bit. `*.f64` dumps and `venv/` are gitignored.
