# analysis::offline — cross-language NULL (numpy)

Dev-only cross-implementation check (NOT in `ctest` — core stays Python-free). numpy recomputes the
1/12-octave RMS-power smoothing and the interference formula from scratch, then NULLs against the C++
`logMagnitudeCurve` / `interferenceDb` dumps.

## Run
```sh
python3 -m venv venv && ./venv/bin/pip install numpy
c++ -std=c++20 -O2 -I ../../../core/include -I ../../include dump.cpp -o dump
./dump . && ./venv/bin/python check.py .
```

## Last result (2026-07-10, numpy 2.5.1)
```
PASS logMagnitudeCurve raw:  max|Δ|=1.186e-11
PASS logMagnitudeCurve norm: max|Δ|=1.180e-11
PASS interferenceDb:         max|Δ|=8.882e-16
ALL NULLS PASS
```
