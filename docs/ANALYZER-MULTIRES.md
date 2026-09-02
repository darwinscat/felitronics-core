<!-- SPDX-License-Identifier: AGPL-3.0-or-later -->

# Multi-resolution spectrum pane (`analysis`) — design note

*Status: crew-reviewed design (codex + deepseek + Fable, 2026-09-02), implemented as
`MultiResSpectrumPane.h`. Companion to `SpectrumPane.h` / `RollingSpectrumTap.h`.*

## The problem a single FFT cannot solve

A log-frequency display wants **constant Q**: the same *fraction of an octave* per reading everywhere.
One FFT gives constant **Hz** per bin. At 48 kHz:

| N | bin | window | 40 Hz is… | 10 kHz band (1/24 oct = 289 Hz) holds… |
|---|---|---|---|---|
| 2048 | 23.4 Hz | 43 ms | half an octave per bin | 12 bins |
| 16384 | 2.9 Hz | 341 ms | 1/10 oct per bin | 99 bins, **341 ms behind** |

Short window: the lows are unreadable. Long window: the lows read, the highs smear in time and the
whole display goes sluggish. TabbyEQ today hands that trade-off to the user as a *Resolution* menu.
Real analyzers do not choose — they run **several windows at once**, long for the low decades, short
for the high ones, stitched by frequency (Smaart's MTW is the well-known example). That is what this
pane does.

## What already exists

- `RollingSpectrumTap` publishes the **most recent `1<<order` samples** at a UI-rate hop, any order up
  to its ring (16384). A short window is a **suffix** of a long one, so **one frame carries every
  tier**. All tiers share the frame's *end*: the short tier reacts to a transient first, the long tier
  confirms it later. Centring the short tiers instead would put their highs 160 ms behind — wrong for
  a live EQ. The frame is read-only to the pane; every tier copies its suffix into scratch and
  windows the copy (windowing in place, as `SpectrumPane` may, would corrupt the next tier's input).
- `SpectrumPane` owns the classic per-bin pipeline and the "liquid" column rule. It stays as it is —
  OrbitAmp consumes it, TabbyEQ keeps it as the fixed-size mode. The new pane is a sibling, not a
  replacement: its level convention, transient response and natural tilt differ.

## The hop gap (found in review)

At 30 fps the consumer's hop is ~1600 samples. A 1024-sample suffix sees 21 ms of every 33 ms: a
click can land in the other 12 ms and never enter the short tier — the region it owns shows nothing,
peak-hold included. Shortening the hop does not help (the tap's mailbox is single-slot; the effective
hop is the UI period). So the pane is told the hop (`coverSamples`) and **a tier shorter than the hop
analyses as many 50 %-overlapped suffix windows as it takes to reach back over it, averaging their
power (Welch)**: at 48 kHz the 1024 tier runs three (the earliest starts 2048 samples back — the
cost of one 4096 FFT), the longer tiers one. No blind gap, and the 1-bin variance drops ~2× for free.

**The hop that happened, not the one requested** (both code reviews caught this): the tap publishes
at a block boundary once the requested hop has elapsed, and not at all while the reader still owns
the slot — with 2048-sample blocks the real hop is 2048, a missed UI tick doubles it, a 100 ms stall
makes it 4800. So `RollingSpectrumTap` now publishes the samples that entered since its previous
publish next to the order (`tryPull (dst, order, hop)`), and the consumer hands *that* to the pane
per frame: a stall simply means more sub-windows (9 of 1024, 2 of 4096 for 4800 samples), never a
gap. The requested hop is only a lower bound.

## The physics that dictates the level convention

Stitching two FFT lengths is only invisible if both tiers **read the same number for the same
signal**. For a sine that is easy (coherent-gain compensation: the peak bin reads the amplitude at
every N). For broadband material it is not: noise power **per bin** is proportional to bin width
(`6σ²/N` with the amplitude normalisation), so a 4× shorter window reads its noise floor **6 dB
higher**, bin for bin. No fixed per-bin scalar survives a seam:

- *max over the column's bins* (today's rule): at the seam the long tier's column holds 4 bins each
  at ¼ power; for independent bins `E[max of 4] = μ·H₄`, `10·log10(H₄/4) = −2.8 dB`; with the Hann's
  adjacent-bin correlation the measured step is **−3.3 dB**, converging to −6 dB as columns widen.
  Not a placement problem — no seam position fixes it.
- *mean per bin*: exactly the 6 dB.

PSD normalisation is tier-invariant for noise, coherent amplitude for tones; no single per-bin number
is both. A **bandwidth-integrated** quantity is: band power = Σ(bin power) over a band of width B is
`∝ B`, independent of N (Parseval). So the pane reads **power in a fractional-octave band** — the
constant-Q / RTA convention, which is the *definition* of a constant-Q display.

### The reading, precisely

For frequency *f* and band `[f·2^(−o/2), f·2^(o/2)]`, `o = bandOctaves` (default 1/24, width
`B(f) = 0.028882·f`):

1. **Per-bin power** `P_k = 4·|X_k|² / (N·Σw²)` — the coherent-gain amplitude `|X_k|/(Σw/2)`
   squared, divided by the window's equivalent noise bandwidth in bins (`N·Σw²/(Σw)²`, which is
   `1.5·N/(N−1)` for the symmetric Hann; both are measured from the window as built). DC and Nyquist
   are **half-width cells** — their one-sided weight ½ falls out of the geometry.
2. **Fractional-bin integration.** Each bin's power is spread uniformly over its cell; the band
   integrates that piecewise-constant density with *fractional* edge overlap, from **double** prefix
   sums (`O(bins)` per tick, `O(1)` per band). Whole-bin membership would alternate 1↔2 bins across
   neighbouring columns — a 3 dB sawtooth on noise. Bins at the −120 dB floor count as zero, so
   silence reads the floor whatever the band width.
3. **Calibration.** A full-scale sine whose main lobe lies inside the band reads **≈ 0 dBFS**
   (measured 0.00 / −0.09 dB at ≥ 3 bins per band). White noise of RMS σ reads
   `10·log10(4σ²B/fs)` — the level of the sine that would carry the band's power, i.e. 2× the
   band's mean square. Disjoint bands add up to the windowed frame's power. A signal *can* read above
   0 dBFS (a full-scale square's fundamental is +2.1 dB), as it can today.
4. **Bin-limited lows.** Where the band is narrower than a bin (`B < B_bin`), the band cannot
   resolve; the reading is the local **cell's** power (the density × one cell), interpolated across
   the two bins the band touches, continuous at `B = B_bin`. A bin-centred sine there reads the peak
   bin alone, **−1.76 dB**; half-way between bins **−3.2 dB** (Hann scallop). Noise reads flat there
   (white) / falls 3 dB/oct (pink), and once bands are wider than a bin rises **+3 dB/oct** (white) /
   flat (pink). That knee is the true resolution limit of the longest tier, shown instead of hidden.

**Consequence for the display tilt:** per-bin analyzers add ~+4.5 dB/oct to make music look flat;
under constant-Q pink noise is already flat, so the natural tilt is ~0…+1.5 dB/oct. The host's tilt
control keeps its meaning (dB/oct added); its *default* wants to differ per mode.

## Tier selection and seams

`tier(f)` = the **shortest** tier whose bin is fine enough: `B_bin_k ≤ B(f) / binsPerBand`. Seam of
tier *k*: `f_k = binsPerBand · B_bin_k / (2^(o/2) − 2^(−o/2))`.

`binsPerBand` defaults to **2** (crew consensus): a 1-bin band is an exponential estimate (5.6 dB
std per tick) and reads a half-bin sine at −3.2 dB right where the boost-sweep search tone lives;
at 2 bins the worst sine is −0.8 dB and the variance halves, for seams one octave higher. Default
tiers `{14, 12, 10}` = 16384 / 4096 / 1024, 1/24 oct, 48 kHz: bin-limited below **101 Hz**, 16384 to
**811 Hz**, 4096 to **3.25 kHz**, 1024 above (21 ms windows, Welch-covered). Seams scale with fs
(96 kHz: 203 / 1623 / 6492 Hz). A 4th tier of 256 would serve above 13 kHz with 5 ms windows — it
works (13 Welch sub-windows over the hop) but buys nothing visible; the set is configurable.

**Seam blend.** Above each seam the reading crossfades — **in power**, where band contributions
add — from the longer tier to the shorter over `blendOctaves` (default 1/3), never wider than the
spacing to the next seam. Band power makes both tiers agree on noise in expectation; a sine's lobe
partly outside a narrow band still disagrees by up to ~1 dB, and the blend turns any residue into a
slope, never a step. The invariance is a statement about expected **power**: the dB the display shows
is the log of a smoothed estimate, and two tiers with different variances carry slightly different
log biases (Jensen) — measured ≤ 0.3 dB at the Fast preset on a rolling stream, bounded by test 6.

## Per-tier state — and why it smooths power, not dB

Each tier keeps its own per-bin **power** one-pole (`smoothCoeff` per tick) and a **peak-hold in dB**
falling `peakFallDb` per tick. `SpectrumPane` smooths in dB; that cannot be kept here. A log-domain
average is a geometric mean whose bias (−2.5 dB on an exponential bin) is uniform only if every bin
is fed the same way — and the Welch tiers are fed by several windows, the long tier's frames overlap
90 % while the short tier's do not. Fable measured the leak: **1.07 dB** step at the 4096 seam on
stationary white noise with dB smoothing, moving with the hop; **0.2 dB** with power smoothing,
unbiased against the analytic band power. Power smoothing also means silence fades as a release
(−1.25 dB/tick at 0.25) rather than collapsing 30 dB in a tick, and the attack settles in ~6 ticks
instead of ~17 — "Medium" is not the same motion as the classic pane's. The peak trace integrates the
per-bin holds — a persistence envelope over the band, ~0.5 dB above the smoothed fill on noise, not
the power of any one instant. `starve()` = hold, then fade at −3 dB/tick (faster than the peak falls,
so the peak never sinks under the fill). A non-finite *sample* is dropped before the transform (one
NaN must not silence a whole tier); every positive bin adds to a band however small (eight −123 dB
components are a −114 dB band); only the final reading is floored.

Per tick: K FFTs (16384 + 4096 + 3×1024 ≈ 1.5× the long one), an `O(bins)` normalise + prefix pass,
then `O(columns)` reads.

## Interface

```cpp
template <int MaxOrder = RollingSpectrumTap::kMaxOrder, int MaxTiers = 4>
struct MultiResSpectrumPaneT {
    float  peakFallDb = 0.8f, smoothCoeff = 0.25f;                      // as SpectrumPane's numbers; the law is on power
    double bandOctaves = 1.0/24, blendOctaves = 1.0/3, binsPerBand = 2.0;
    int    coverSamples = 0;                                            // the consumer's hop → Welch sub-windows
    int    setTiers (const int* orders, int count);                     // clamp [8, MaxOrder], dedup, longest-first; the one allocating call; returns the count
    void   reset();                                                     // silence; the next frame seeds (no fade-in)
    int    frameOrder() const;                                          // request THIS from the tap = the longest tier
    float* frameInput(); void ingest (int order); void starve();
    template <class Emit> void buildColumns (const PlotMap&, double fs, double tiltDbPerOct, double tiltPivotHz, Emit&&) const;   // same emit as SpectrumPane
    // introspection: tierCount / tierOrder / tierBins / seamHz (k, fs) / tierAt (f, fs) / enbwBins (k) / subWindows (k, order)
    //                tierBinDb / tierBinPeak (k, i) · tierBandDb / tierBandPeakDb (k, f, fs) · readDb / readPeakDb (f, fs)
};
using MultiResSpectrumPane = MultiResSpectrumPaneT<>;
```

Storage is fixed and flat-packed across tiers (Σ bins over distinct orders ≤ `kMaxSize + MaxTiers`);
`ingest` / `starve` / `buildColumns` / reads never allocate. The object is ~0.8 MB at order 14 —
hold it by `unique_ptr`. A frame shorter than a tier leaves that tier holding, and reads fall back
to the tiers that have seen a frame; `ingest` clamps the order; f above Nyquist repeats the last band
that fits (tier choice and blend included). Message-thread only.

**Consumer notes (TabbyEQ):** request order `frameOrder()` from the tap and pass the hop the tap
*reported* for each frame as `coverSamples`; keep the classic `specResolution` for the fixed-size
mode and store the multi flag separately; `buildSpectrumPaths` becomes a template over the pane
type; the tilt default drops per mode; the analyzer speed presets map 1:1 onto `smoothCoeff` /
`peakFallDb` (a different motion, see above); reset the pane pair you switch TO (`SpectrumPane`
gained `reset()` for this) so it seeds instead of resuming a stale spectrum.

## Tests (`felitronics_multires_spectrum_tests`, property, no golden files)

1. **Direct-DFT oracle** — each tier's bins equal a hand-computed DFT of the exact suffix within
   0.02 dB (pins the offset, the window, `4|X|²/(N·Σw²)`, the packed layout).
2. **Sine oracle** — a full-scale sine swept 30 Hz–20 kHz: |reading| ≤ 0.3 dB wherever the band holds
   the lobe (≥ 4 bins), never below −3.5 dB or above +0.2 dB; the bin-limited −1.76 / −3.2 dB pins; a
   sweep through every seam + blend steps ≤ 0.5 dB per 1/60 oct.
3. **DC half cell** — a unit constant: raw DC bin `4/ENBW`, a bin-limited read at DC half of it.
4. **Welch cover** — sub-window counts at hop 0 / 1024 / 1025 / 1600 / 4800 / INT_MAX; a click inside
   the hop gap is invisible at cover 0 and registers at cover 1600; each tier's bins equal the **mean of
   direct DFTs of every sub-window at its exact offset** (a 75 %-overlap mutant passed the suite
   without this one). `RollingSpectrumTap` reports the hop that happened (block rounding, a missed tick).
5. **Time resolution** — a burst in the last 512 samples: 1024 tier ≈ −3 dB (the trailing half of its
   Hann), 16384 tier below −40 dB (its Hann tail); a burst centred under the 1024 window: ≈ 0 dB vs
   25+ dB down.
6. **Noise ensemble** (fixed PRNG, 800 independent frames, α = 0.05) — at every seam the longer and
   shorter tier agree within 0.6 dB over the blend zone; +3 dB/oct across all tiers; the analytic
   `10·log10(4σ²B/fs)` at 1.2 kHz within 0.5 dB; the seam residue at hops 800 / 1600 / 3200 (the
   reading must not depend on the UI rate); pink noise flat within 0.7 dB/oct; and a **rolling
   stream** as the tap delivers it (90 %-overlapped frames) at the Fast preset, residue averaged over
   300 ticks: fill bias ≤ 0.5 dB, peak-trace bias ≤ 0.8 dB at every hop.
7. **Parseval** — contiguous bands tiling 400 Hz–20 kHz sum to the bin power over the same interval
   within 0.02 dB (fractional edges, no double counting).
8. **Fractional integration** — 60 random bands (1/200–1/4 oct, 20 Hz–12 kHz, all tiers) vs an
   `O(N)` reference including the bin-limited rule and the half cells, within 0.01 dB.
9. **Smoothing + peak law** — the one-pole on power and the 0.8 dB/tick hold, tick by tick; starve
   holds 15 ticks then fades; `reset` is silence.
10. **Silence** — every band at the floor after reset and on a zero frame; a −100 dB tone leaves the
    far bands at the floor; eight −123 dB components read as one −114 dB band, not the floor.
11. **Adversarial** — NaN/Inf samples (dropped: the sine next to them survives within 0.3 dB); f ≤ 0,
    NaN, Inf, above Nyquist; fs ≤ 0 / NaN; a frame shorter than the longest tier (it holds, the others
    update, reads fall back to them); order 99; blend width 0; a non-positive band width; a bin-limited
    read at Nyquist is half the Nyquist bin; tiers {14,13,12,11} with a 1.5-oct blend sweep without a
    step; unsorted / duplicate / out-of-range / empty / null tier lists.
12. **Column geometry + no-alloc** — N+1 ascending points from x = 0, every column is
    `specDbToY (readDb + tilt)`, the 257-column minimum; ingest / starve / buildColumns / reads under
    the global operator-new counter.
