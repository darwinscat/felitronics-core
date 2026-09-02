<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Multi-resolution spectrum pane (`analysis`) — design note

*Status: design under crew review (2026-09-02). Companion to `SpectrumPane.h` / `RollingSpectrumTap.h`.*

## The problem a single FFT cannot solve

A log-frequency display wants **constant Q**: the same *fraction of an octave* per reading everywhere.
One FFT gives constant **Hz** per bin. At 48 kHz:

| N | bin | window | 40 Hz is… | 10 kHz column (1/24 oct = 289 Hz) holds… |
|---|---|---|---|---|
| 2048 | 23.4 Hz | 43 ms | half an octave per bin | 12 bins |
| 16384 | 2.9 Hz | 341 ms | 1/14 oct per bin | 99 bins, **341 ms behind** |

Short window: the lows are unreadable. Long window: the lows read, the highs smear in time and the
whole display goes sluggish. TabbyEQ today hands that trade-off to the user as a *Resolution* menu.
Real analyzers do not choose — they run **several windows at once**, long for the low decades, short
for the high ones, stitched by frequency (Smaart's MTW is the well-known example). That is what this
pane does.

## What already exists

- `RollingSpectrumTap` publishes the **most recent `1<<order` samples** at a UI-rate hop, any order up
  to its ring (16384). A short window is a **suffix** of a long one, so **one frame carries every
  tier**: tier *k* transforms the last `N_k` samples of the frame. All tiers share the frame's *end*
  time — the short tier reacts to a transient first, the long tier confirms it later. That is the
  intended behaviour ("the top does not smear"), not an artefact.
- `SpectrumPane` owns the per-bin pipeline (Hann → real FFT → dB → one-pole smoothing + peak-hold)
  and the "liquid" column rule. It stays as it is — OrbitAmp consumes it, TabbyEQ keeps it as the
  classic fixed-size mode. The new pane is a sibling, not a replacement.

## The physics that dictates the level convention

Stitching two FFT lengths is only invisible if both tiers **read the same number for the same
signal**. For a sine that is easy (Hann single-bin compensation: the peak bin reads the amplitude
at every N). For broadband material it is not: noise power **per bin** is proportional to bin width
(`6σ²/N` with Hann amplitude normalisation), so a 4× shorter window reads its noise floor **6 dB
higher**, bin for bin. No per-bin rule survives a seam:

- *max over the column's bins* (today's rule): at the seam the long tier's column holds 4 bins each
  at ¼ power, `E[max of 4] ≈ 0.52×` the short tier's single bin → **−2.8 dB step**, and the ratio
  converges to −6 dB as columns widen. Not a placement problem — no seam position fixes it.
- *mean per bin*: exactly the 6 dB.

Only a **bandwidth-integrated** quantity is tier-invariant: band power = Σ(bin power) over a band
of width B is `∝ B`, independent of N (Parseval). So the pane reads **power in a fractional-octave
band** — the constant-Q / RTA convention — and that is the *definition* of a constant-Q display, not
merely a fix.

### The reading, precisely

For frequency *f* and band width `B(f) = f·(2^(o/2) − 2^(−o/2))`, `o = bandOctaves` (default 1/24):

1. **Fractional-bin integration.** Each bin's power is spread uniformly over its width; the band
   integrates that piecewise-constant density with *fractional* edge overlap, from **double**
   prefix sums (`O(bins)` per tick, `O(1)` per column). Whole-bin membership would alternate 1↔2
   bins across neighbouring columns — a 3 dB sawtooth on noise. Fractional edges are continuous.
2. **ENBW normalisation.** Divide by the window's equivalent noise bandwidth in bins
   (`N·Σw²/(Σw)²` = 1.5 for Hann, computed from the actual window). A full-scale sine whose main
   lobe lies inside the band then reads **0 dBFS**; band power is Parseval-exact (the bands sum to
   the signal power). Nothing reads above 0 dBFS.
3. **Bin-limited lows.** Where the band is narrower than a bin (`B < B_bin`), the band cannot
   resolve; the reading is the local **bin's** power (density × `B_bin`, i.e. the integral scaled
   by `B_bin/B`), interpolated across the bin pair the band touches. Continuous at `B = B_bin`. A
   sine there reads the peak bin alone: **−1.76 dB** (the rest of its lobe sits in the neighbours).
   Noise reads flat there and rises **+3 dB/oct** (white) / flat (pink) once bands are wider than a
   bin. That knee is the true resolution limit of the longest tier, made visible instead of hidden.

**Consequence for the display tilt:** per-bin analyzers add ~+4.5 dB/oct to make music look flat;
under constant-Q pink noise is already flat, so the natural tilt is ~0…+1.5 dB/oct. The host's
tilt control keeps its meaning (dB/oct added); its *default* wants to differ per mode.

## Tier selection and seams

`tier(f)` = the **shortest** tier whose bin is fine enough: `B_bin_k ≤ B(f) / binsPerBand`
(`binsPerBand` default 1). Seam of tier *k*: `f_k = binsPerBand · B_bin_k / (2^(o/2) − 2^(−o/2))`.

Default tiers `{14, 12, 10}` = 16384 / 4096 / 1024, 1/24 oct, 48 kHz: bin-limited below **101 Hz**,
16384 to **405 Hz**, 4096 to **1.62 kHz**, 1024 above (21 ms window). Seams scale with fs (96 kHz:
203 / 810 / 3240 Hz). 4 tiers (`…, 8`) would give 256-sample windows above ~13 kHz — 5 ms of a
33 ms tick — no visible gain; three is the default, the set is configurable.

**Seam blend.** Above each seam the reading crossfades **in dB** from the longer tier to the shorter
over `blendOctaves` (default 1/3). Band power makes both tiers agree exactly on noise; a sine near a
1-bin-wide band can still disagree by up to ~1.8 dB (lobe partly outside, scallop), and the blend
turns any residue into a slope, never a step.

## Per-tier state and cost

Each tier keeps its own per-bin **dB** smoothing (one-pole per tick, `smoothCoeff`) and peak-hold
(`peakFallDb`) — SpectrumPane's law, so the "liquid" feel stays; the short tiers move faster simply
because their window is shorter. `starve()` = hold, then fade after ~0.5 s, as before. Per tick:
K FFTs (16384 + 4096 + 1024 ≈ 1.3× the long one), an O(bins) dB→power pass + prefix sums, then
O(columns) reads. Non-finite bins are treated as the floor (a NaN must not poison the smoothing).

## Interface

```cpp
template <int MaxOrder = RollingSpectrumTap::kMaxOrder, int MaxTiers = 4>
struct MultiResSpectrumPaneT {
    float  peakFallDb = 0.8f, smoothCoeff = 0.25f;   // as SpectrumPane
    double bandOctaves = 1.0/24, blendOctaves = 1.0/3, binsPerBand = 1.0;
    void   setTiers (const int* orders, int count);   // sorted desc, dedup, clamped to [8, MaxOrder]; FFT plans prepared here (message thread — the one allocating call)
    int    frameOrder() const;                        // what to request from the tap = the longest tier
    float* frameInput(); void ingest (int order); void starve();
    template <class Emit> void buildColumns (const PlotMap&, double fs, double tiltDbPerOct, double tiltPivotHz, Emit&&) const;   // same emit as SpectrumPane
    int    tierAt (double f, double fs) const; double readDb (double f, double fs) const; double readPeakDb (double f, double fs) const;   // introspection / tests
};
using MultiResSpectrumPane = MultiResSpectrumPaneT<>;
```

Storage is fixed and flat-packed across tiers (Σ bins over distinct orders ≤ `kMaxSize + MaxTiers`);
`ingest` / `starve` / `buildColumns` never allocate. A frame shorter than a tier leaves that tier
holding (starved); `ingest` clamps the order. Message-thread only; the audio thread is untouched.

## Tests (property, no golden files)

1. **Sine oracle** — full-scale sine in every tier's region and in the bin-limited lows: reading in
   [−3, 0] dB everywhere; within ±0.3 dB of 0 where the band holds the lobe (≥ 3 bins/band).
2. **Noise seam invariance** — white noise, converged: the reading just below vs just above each
   seam (and through the blend) differs by ≤ 1 dB; the local slope is +3 dB/oct ± tolerance. Pink
   noise reads flat ± tolerance across all seams.
3. **Time resolution** — a burst in the last 512 samples of the frame: the 1024 tier reads within a
   few dB of steady state after one tick, the 16384 tier ~15 dB down (512/16384) — the tiers
   provably read different time windows.
4. **Tier map** — `tierAt` non-increasing in f; seams match the formula; they double at 2× fs.
5. **Dynamics parity** — smoothing convergence, 0.8 dB/tick peak-hold decay, hold-then-fade starve:
   the same numbers `SpectrumPane` pins.
6. **Adversarial** — unsorted / duplicate / out-of-range tiers; a frame shorter than the longest
   tier; order out of range; f outside (0, fs/2]; width 0 map; NaN/Inf samples; ingest without fill.
7. **No-alloc** — `ingest` + `buildColumns` under the global operator-new counter.
8. **Column geometry** — N+1 ascending points, i = 0 at x = 0, tilt about the pivot.
