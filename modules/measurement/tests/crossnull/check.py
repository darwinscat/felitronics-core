#!/usr/bin/env python3
# Independent numpy NULL of felitronics::measurement C++ dumps. Each reference is recomputed
# from scratch (numpy.convolve = DIRECT time-domain, numpy.fft, analytic ESS) — a genuine
# cross-implementation null, not a re-run of the same code.
import sys, numpy as np

DIR = sys.argv[1] if len(sys.argv) > 1 else "."
def rd(n): return np.fromfile(f"{DIR}/{n}", dtype="<f8")

fails = 0
def null(name, got, want, atol):
    global fails
    got = np.asarray(got); want = np.asarray(want)
    if got.shape != want.shape:
        print(f"FAIL {name}: shape {got.shape} != {want.shape}"); fails += 1; return
    err = float(np.max(np.abs(got - want))) if got.size else 0.0
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    rel = err / scale if scale > 0 else err
    tag = "PASS" if err <= atol else "FAIL"
    if tag == "FAIL": fails += 1
    print(f"{tag} {name}: max|Δ|={err:.3e}  rel={rel:.3e}  (atol={atol:.0e}, n={got.size})")

# (1) convolve vs numpy.convolve (direct time-domain — independent of any FFT)
for c in (0, 1):
    x, h, y = rd(f"conv{c}_x.f64"), rd(f"conv{c}_h.f64"), rd(f"conv{c}_y.f64")
    ref = np.convolve(x, h)                       # full linear conv, O(n^2), no FFT
    null(f"convolve[{c}] (n={x.size}x{h.size})", y, ref, 1e-9)

# (2) magSpectrum vs numpy.fft
x = rd("mag_x.f64"); m = rd("mag.f64")
nfft = 2048
xp = np.zeros(nfft); xp[:x.size] = x
ref = np.abs(np.fft.fft(xp))[:nfft // 2]
null(f"magSpectrum (nfft={nfft})", m, ref, 1e-9)

# (3) makeSweep signal vs analytic ESS recompute
spec = {}
for line in open(f"{DIR}/sweep_spec.txt"):
    k, v = line.split(); spec[k] = float(v)
f1, f2, T, sr = spec["f1"], spec["f2"], spec["dur"], spec["sr"]
amp, fade = spec["amp"], spec["fade"]
N = int(spec["sweepLen"])
w1, w2 = 2*np.pi*f1, 2*np.pi*f2
R = np.log(w2 / w1)
n = np.arange(N)
t = n / sr
phi = (w1 * T / R) * (np.exp((t / T) * R) - 1.0)
env = np.ones(N)
Nf = int(round(fade * sr))
if Nf > 0:
    i0 = n < Nf
    env[i0] = 0.5 * (1 - np.cos(np.pi * n[i0] / Nf))
    i1 = n >= N - Nf
    env[i1] = 0.5 * (1 - np.cos(np.pi * (N - 1 - n[i1]) / Nf))
ref_sweep = amp * np.sin(phi) * env
sweep = rd("sweep.f64")
null(f"makeSweep signal (N={N})", sweep[:N], ref_sweep, 1e-12)
# tail must be exact silence
null("makeSweep tail silence", sweep[N:], np.zeros(sweep.size - N), 0.0)

# (3b) sweep_proper ⊛ inverse ≈ unit-peak δ (independent numpy MAC of C++ inverse)
inv = rd("inverse.f64")
d = np.convolve(sweep[:N], inv)                  # numpy direct conv of the C++-built inverse
pk = int(np.argmax(np.abs(d)))
peak = abs(d[pk])
far = np.concatenate([d[:max(0, pk-2000)], d[pk+2000:]])
floor = float(np.max(np.abs(far))) if far.size else 0.0
print(f"{'PASS' if abs(peak-1.0) < 0.05 else 'FAIL'} sweep⊛inverse peak={peak:.4f} (want ~1.0), "
      f"floor/peak={floor/peak:.2e} @ ±2000 from δ")
if abs(peak - 1.0) >= 0.05: fails += 1

print("\n" + ("ALL NULLS PASS" if fails == 0 else f"{fails} FAILURES"))
sys.exit(1 if fails else 0)
