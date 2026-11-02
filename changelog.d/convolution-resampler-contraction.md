### convolution — `resampleIr` pins its arithmetic against FMA: a GCC x86-64-v3 build resamples to the bits of a baseline one

GitHub's `ubuntu-latest` moved to Ubuntu 26.04 in October 2026. That image is the amd64v3 build Canonical's cloud
images carry, and its GCC (14.3 and 15.2) is configured for `-march=x86-64-v3`, which has an FMA. Under the
policy's `-ffp-contract=on` (law 10), GCC fused `acc += in[k] * w` in both of `resampleIr`'s accumulations. It did
not fuse the test's frozen unmemoized copy, whose `in[k]` is a call to `std::vector::operator[]`. The P70
differential case with a rectangular window then differed at output 2: its normalised weights are 1/8, 3/8, 3/8,
1/8 to an ulp, so an exact answer can sit on a float32 midpoint. `main` was red on the three ubuntu-latest rows.

The products in the kernel that meet an add now go through `core::det::mulAdd` / `core::det::mul`, a `volatile`
pin that holds under every contraction mode. Those are the accumulation, the window's `1 - r²` and the weight
`wsum` sums. The test's copy is pinned the same way.

At x86-64-v3 no FMA is left in either copy, including under `-ffp-contract=fast`. `besselI0`'s series stays
unpinned: it is two statements, both copies share it, and no FMA was emitted there. A pin there cost 8% for no
bit that moved.

On a mismatch the P70 case now prints both floats in hex. The CI Linux rows print whether the compiler's default
target has FMA.

Measured over 4.0 M output floats. The corpus is every standard rate pair from 22.05 to 192 kHz and four kernel
configs, on inputs built from integer bits:
- **Baseline x86-64:** bit-identical to before.
- **x86-64-v3:** now bit-identical to baseline, at `on` and at `fast`. Before, the two differed in 4337 floats: the
  same input resampled differently depending on the x86 target a consumer compiled for.
- **arm64 (Apple clang):** 4338 floats moved by the unfusing. One output float now differs from glibc
  x86-64's; `sin` is the one libm call in the kernel.

Cost on a 4 s IR at 48 → 44.1 kHz: 11.1 → 12.7 ms on baseline x86-64 (i9-13900H, GCC 15.2), 8.0 → 7.1 ms on
arm64 macOS (best of 7, three runs each).
