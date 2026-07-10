#!/usr/bin/env python3
# Independent numpy NULL of felitronics::analysis::offline C++ dumps. The 1/12-oct RMS-power smoothing
# and the interference formula are recomputed from scratch (numpy FFT + hand arithmetic), then nulled
# against the C++ curve/interference dumps.
import sys, numpy as np

DIR = sys.argv[1] if len(sys.argv) > 1 else "."
def rd(n): return np.fromfile(f"{DIR}/{n}", dtype="<f8")

spec = {}
for line in open(f"{DIR}/spec.txt"):
    k, v = line.split(); spec[k] = float(v)
sr = spec["sr"]; N = int(spec["N"]); fLo = spec["fLo"]; fHi = spec["fHi"]
points = int(spec["points"]); octDiv = spec["octaveDivisions"]; minNfft = int(spec["minNfft"])

fails = 0
def null(name, got, want, atol):
    global fails
    got = np.asarray(got); want = np.asarray(want)
    if got.shape != want.shape:
        print(f"FAIL {name}: shape {got.shape} != {want.shape}"); fails += 1; return
    err = float(np.max(np.abs(got - want))) if got.size else 0.0
    print(f"{'PASS' if err <= atol else 'FAIL'} {name}: max|Δ|={err:.3e} (atol={atol:.0e}, n={got.size})")
    if err > atol: fails += 1

def next_pow2(n):
    p = 1
    while p < n: p <<= 1
    return p

def curve(ir, normalize):
    nfft = max(minNfft, next_pow2(len(ir)))
    M = np.abs(np.fft.fft(ir, nfft))[:nfft // 2]        # matches magSpectrum: nfft/2 bins, no Nyquist
    binHz = sr / nfft
    nyq = 0.5 * sr
    fhi = min(fHi, nyq)
    half = 2.0 ** (1.0 / (2.0 * octDiv))
    out = np.empty(points)
    lastBin = len(M) - 1
    for i in range(points):
        f = fLo * (fhi / fLo) ** (i / (points - 1))
        blo = int(np.floor(f / half / binHz)); bhi = int(np.ceil(f * half / binHz))
        blo = max(1, blo); bhi = min(lastBin, max(bhi, blo))
        seg = M[blo:bhi + 1]
        out[i] = 10.0 * np.log10(np.mean(seg ** 2) + 1e-20) if seg.size else -200.0
    if normalize:
        pk = out.max()
        if pk > -199.0: out -= pk
    return out

ir = rd("ir.f64")
null("logMagnitudeCurve raw",  rd("curve_raw.f64"),  curve(ir, False), 1e-9)
null("logMagnitudeCurve norm", rd("curve_norm.f64"), curve(ir, True),  1e-9)

micA, micB, coh = rd("micA.f64"), rd("micB.f64"), rd("coh.f64")
pw = 10.0 ** (micA / 10.0) + 10.0 ** (micB / 10.0)
ref_interf = coh - 10.0 * np.log10(pw + 1e-20)          # both mics audible
null("interferenceDb", rd("interf.f64"), ref_interf, 1e-9)

print("\n" + ("ALL NULLS PASS" if fails == 0 else f"{fails} FAILURES"))
sys.exit(1 if fails else 0)
