// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>
#include "RollingSpectrumTap.h"                        // kMaxOrder — the frame the tap can publish
#include <felitronics/analysis/PlotMap.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

//==============================================================================
// MultiResSpectrumPaneFast — the SAME analyzer as MultiResSpectrumPane, computed differently.
//
// Every number this pane draws is defined by its sibling and by docs/ANALYZER-MULTIRES.md: one frame
// feeding several tiers, a reading that is the POWER IN A FRACTIONAL-OCTAVE BAND integrated from
// double prefix sums with fractional edges, DC/Nyquist half cells, power-domain smoothing, a power
// crossfade above each seam, a peak trace that falls peakFallDb per tick. NOTHING about the physics
// changes here — read the sibling header for it. This one only stops paying for the same answer
// twice. It is a SIBLING, not a replacement: MultiResSpectrumPane keeps its behaviour and its tests.
//
// Three exact changes, none of them an approximation:
//
//   • THE PEAK TRACE LIVES IN POWER. The sibling stores the peak in dB, which costs a log10 for every
//     bin every tick to produce the dB the peak law compares against, and then an exp for every bin
//     to convert the hold BACK to power for its prefix sum — 10 755 of each at the default ladder, a
//     transform and its own inverse. In power the law "fall peakFallDb per tick" is a multiply by the
//     constant 10^(−peakFallDb/10), so both vanish. max() commutes with a monotone map, so this is
//     the same trace, not an approximation of it. The per-bin dB array goes with them; tierBinDb
//     computes its one value on demand.
//
//   • ONE COLUMN PLAN, SHARED BY BOTH TRACES. In the sibling, readDb and readPeakDb differ by a
//     single pointer — the prefix array — yet each re-derives the whole geometry of the column: which
//     tier owns this frequency, the seam, the blend width and weight, the band's fractional bin
//     edges, the display tilt. All of it depends on the plot map, the sample rate and the tuning, and
//     none of it on which trace is being read or on the frame that just arrived. So it is computed
//     ONCE per column into a fixed-size plan and applied to both prefix arrays: a band becomes four
//     prefix loads and a multiply-add. The plan is rebuilt only when its inputs change (a resize, a
//     new sample rate, a tuning edit, a tier becoming valid) — never per tick.
//
//   • DC AND NYQUIST ARE PEELED OUT of the magnitude loop, leaving a contiguous interior the compiler
//     can see the bounds of. Bin for bin the arithmetic is unchanged.
//
// The reads (readDb / readPeakDb / tierBandDb / tierAt / seamHz …) are deliberately left as the
// sibling writes them: they are the reference the plan is NULL-tested against, and the diagnostic
// path tests use. Only buildColumns takes the fast route.
//
// Message-thread only. prepare/setTiers is still the one allocating call; ingest, starve,
// buildColumns and the plan rebuild allocate nothing.
namespace felitronics::analysis
{

// The FFT backend is a template parameter constrained to the packed-Hermitian layout the bin loop reads:
// the scalar reference by default, or fftpffft::PffftOrderedRealFft for SIMD (the tiers' 16384 + 4096 +
// 3×1024 points per tick are the one place the analyzer's cost lives).
template <int MaxOrder = RollingSpectrumTap::kMaxOrder, int MaxTiers = 4,
          felitronics::core::fft::PackedHermitianSpectrum Fft = felitronics::core::fft::DefaultRealFft>
struct MultiResSpectrumPaneFastT
{
    using FftType = Fft;

    static constexpr int kMaxOrder = MaxOrder;
    static constexpr int kMinOrder = 8;                                  // 256 — below that a tier is a smear, not a window
    static constexpr int kMaxSize  = 1 << MaxOrder;
    static constexpr int kMaxBins  = kMaxSize / 2 + 1;
    static constexpr int kMaxTiers = MaxTiers;
    static_assert (MaxOrder >= kMinOrder && MaxOrder <= 24, "MultiResSpectrumPaneFast: MaxOrder out of range");
    static_assert (MaxTiers >= 1 && MaxTiers <= MaxOrder - kMinOrder + 1, "MultiResSpectrumPaneFast: more tiers than distinct orders");

    // Capacity of the flat-packed per-bin storage: the sum of (N/2 + 1) over any MaxTiers DISTINCT orders
    // ≤ MaxOrder is < 2^MaxOrder + MaxTiers (a geometric tail plus one per tier). Windows: Σ N < 2·kMaxSize.
    static constexpr int kBinCapacity = kMaxSize + MaxTiers;
    static constexpr int kSumCapacity = kBinCapacity + MaxTiers;          // one extra prefix entry per tier
    static constexpr int kWinCapacity = 2 * kMaxSize;

    static constexpr float  kFloorDb    = -200.0f;                       // SpectrumPane's floor too: deep below every plot bottom
    static constexpr double kFloorPower = 1.0e-20;                       // 10^(kFloorDb/10): a READING below it is the floor
    static constexpr double kMaxPower   = 1.0e12;                        // +120 dB: a bin power is clamped here before it narrows to float
    // Denormal guard for the smoothed POWER state. On an all-zero frame the fill smoother is
    // pw ← (1−c)·pw, which in float32 descends into the subnormals and STICKS there: at the smallest
    // subnormal c·pw rounds to zero, so pw stops changing and every bin is left doing subnormal
    // arithmetic for the rest of the session. On x86 that is measurable — a tick on silence cost 3.0×
    // a tick on signal (505 vs 168 us, i9-13900H, FTZ off, the message thread never sets it) — and it
    // gets WORSE the longer the transport stays stopped, which is exactly backwards. So flush the
    // state to a true zero once it is far below anything a reading can use.
    // NOT core::flushDenormal: its 1e-15 is written for AMPLITUDE state, and here the array holds
    // POWER, where 1e-15 is −150 dB — well inside the range this pane deliberately sums (v0.22.2:
    // "every positive bin counts, however small"). At 1e-30 all kBinCapacity bins together carry
    // < 1e-25 (−250 dB), below the −200 floor even under the steepest display tilt, so no reading
    // moves; and 1e-30 is still ~8 binades above the smallest NORMAL float (1.18e-38).
    static constexpr float  kFlushPower = 1.0e-30f;
    static constexpr int    kMaxColumns = 901;                           // buildColumns emits i = 0..N, N ≤ 900
    // The floor as the peak array actually stores it. The sibling's peak rests at kFloorDb and its
    // prefix sum drops anything AT the floor by testing `pk > kFloorDb`; here the equivalent test is
    // against this sentinel and not against kFloorPower, so "a hold resting on the floor counts as
    // nothing" holds BY CONSTRUCTION. Comparing the stored float to the double kFloorPower would work
    // today only because 1e-20 happens to round DOWN in float (9.9999997e-21); a floor whose constant
    // rounded the other way would silently start contributing ~1e-20 from every bin, and a freshly
    // reset pane would read −161 dB instead of the floor.
    static constexpr float  kFloorPowerF = (float) kFloorPower;

    //==============================================================================
    // Tuning (message thread; take effect on the next tick / read).
    float  peakFallDb   = 0.8f;          // peak-hold decay per tick (~24 dB/s at 30 Hz)
    float  smoothCoeff  = 0.25f;         // per-tick one-pole on POWER toward the new frame (the analyzer "speed")
    double bandOctaves  = 1.0 / 24.0;    // width of the constant-Q band a reading integrates
    double blendOctaves = 1.0 / 3.0;     // crossfade width above each seam (0 = hard seam)
    double binsPerBand  = 2.0;           // a tier is used where its bin is at least this many times finer than the band
    int    coverSamples = 0;             // the hop the tap reported for the frame; a tier shorter than it Welch-averages sub-windows over it

    MultiResSpectrumPaneFastT()
    {
        const int def[3] = { MaxOrder, MaxOrder - 2, MaxOrder - 4 };     // 16384 / 4096 / 1024 at order 14
        setTiers (def, 3);
    }

    //==============================================================================
    // Setup — the ONE allocating call (FFT plans). Orders are clamped to [kMinOrder, kMaxOrder],
    // de-duplicated and kept longest-first; more distinct orders than MaxTiers keeps the longest ones; an
    // empty list means the single longest window. Resets every tier's state (the storage layout changes);
    // the next frame seeds directly. Returns the number of tiers in effect.
    int setTiers (const int* orders, int count)
    {
        unsigned mask = 0;
        for (int i = 0; i < count && orders != nullptr; ++i)
            mask |= 1u << (unsigned) std::clamp (orders[i], kMinOrder, kMaxOrder);
        if (mask == 0) mask = 1u << (unsigned) kMaxOrder;

        numTiers = 0;
        int binOff = 0, winOff = 0, sumOff = 0;
        for (int o = kMaxOrder; o >= kMinOrder && numTiers < MaxTiers; --o)
        {
            if ((mask & (1u << (unsigned) o)) == 0) continue;
            Tier& t = tiers[(std::size_t) numTiers];
            t.order = o; t.size = 1 << o; t.bins = t.size / 2 + 1;
            t.binOffset = binOff; t.winOffset = winOff; t.sumOffset = sumOff;
            binOff += t.bins; winOff += t.size; sumOff += t.bins + 1;

            fft[(std::size_t) numTiers].prepare (t.size);              // alloc here (setup), never per tick

            // Symmetric Hann; the amplitude compensation and the ENBW come from the window as built, not
            // from the textbook constants, so the calibration is exact for THIS window.
            double sw = 0.0, sw2 = 0.0;
            constexpr double pi = 3.14159265358979323846;
            for (int i = 0; i < t.size; ++i)
            {
                const double w = 0.5 - 0.5 * std::cos (2.0 * pi * (double) i / (double) (t.size - 1));
                window[(std::size_t) (t.winOffset + i)] = (float) w;
                sw += w; sw2 += w * w;
            }
            t.enbwBins  = (double) t.size * sw2 / (sw * sw);             // 1.5·N/(N−1) for this Hann
            t.powerNorm = 4.0 / ((double) t.size * sw2);                 // P_k = 4|X_k|² / (N·Σw²) = (|X_k|/(Σw/2))² / ENBW
            ++numTiers;
        }
        ++layoutEpoch;                                                   // the prefix offsets moved
        reset();
        return numTiers;
    }

    // Forget every frame: all bins to silence, the next frame seeds directly (no fade-in). For transport
    // jumps / stream discontinuities. Keeps the tiers and the tuning.
    void reset() noexcept
    {
        frameArmed  = false;
        starveTicks = 0;
        for (int k = 0; k < numTiers; ++k)
        {
            Tier& t = tiers[(std::size_t) k];
            t.valid = false;
            for (int i = 0; i < t.bins; ++i)
            {
                specPow    [(std::size_t) (t.binOffset + i)] = 0.0f;
                specPeakPow[(std::size_t) (t.binOffset + i)] = kFloorPowerF;   // == kFloorDb
            }
            rebuildSums (t);
        }
        planKey = PlanKey {};                                            // the layout moved: the plan is stale
    }

    int tierCount()         const noexcept { return numTiers; }
    int tierOrder (int k)   const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].order : 0; }
    int tierBins  (int k)   const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].bins  : 0; }
    int frameOrder()        const noexcept { return numTiers > 0 ? tiers[0].order : kMaxOrder; }   // request THIS from the tap
    int frameSize()         const noexcept { return 1 << frameOrder(); }
    double enbwBins (int k) const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].enbwBins : 1.0; }

    // Per-bin introspection (tests / diagnostics): the smoothed power's dB and the peak-hold of tier k, bin i.
    float tierBinDb   (int k, int i) const noexcept { return inRange (k, i) ? (float) toDb ((double) specPow    [(std::size_t) (tiers[(std::size_t) k].binOffset + i)]) : kFloorDb; }
    float tierBinPeak (int k, int i) const noexcept { return inRange (k, i) ? (float) toDb ((double) specPeakPow[(std::size_t) (tiers[(std::size_t) k].binOffset + i)]) : kFloorDb; }
    // The RAW smoothed power behind tierBinDb. tierBinDb floors at −200 dB, so it cannot tell a bin
    // that reached zero from one stuck in the subnormals — this is what the denormal-flush test reads.
    float tierBinPower (int k, int i) const noexcept { return inRange (k, i) ? specPow [(std::size_t) (tiers[(std::size_t) k].binOffset + i)] : 0.0f; }

    // How many sub-windows tier k analyses per frame for the current coverSamples and a frame of `order`:
    // the fewest half-overlapped windows whose earliest one starts at or before n − coverSamples (one
    // window already spans t.size), never more than the frame holds.
    int subWindows (int k, int order) const noexcept
    {
        if (k < 0 || k >= numTiers) return 0;
        const Tier& t = tiers[(std::size_t) k];
        const int n = 1 << std::clamp (order, 0, kMaxOrder);
        if (t.size > n) return 0;
        const int hop   = t.size / 2;
        const int cover = std::clamp (coverSamples, 0, kMaxSize);
        const int need  = 1 + (std::max (0, cover - t.size) + hop - 1) / hop;
        const int room  = 1 + (n - t.size) / hop;
        return std::clamp (need, 1, room);
    }

    //==============================================================================
    // Frames. Fill frameInput() with (1<<order) samples, then ingest (order) — exactly one ingest per fill
    // (the debug guard trips otherwise). The buffer holds any order ≤ kMaxOrder. A tier longer than the
    // frame cannot be computed from it and simply holds (it is starved for that tick); the normal case is
    // order == frameOrder(), which serves every tier. The frame is read-only to the pane.
    float* frameInput() noexcept { frameArmed = true; return frame.data(); }

    void ingest (int order) noexcept
    {
        fallMul = fallMultiplier();
        assert (frameArmed && "MultiResSpectrumPaneFast: fill frameInput() before every ingest()");
        frameArmed  = false;
        starveTicks = 0;
        order = std::clamp (order, 0, kMaxOrder);
        const int n = 1 << order;
        for (int k = 0; k < numTiers; ++k)
        {
            Tier& t = tiers[(std::size_t) k];
            if (t.size > n) continue;                                    // frame too short for this tier: hold
            transformTier (k, n, subWindows (k, order), ! t.valid);      // first frame after reset seeds (no fade-in)
            t.valid = true;
        }
    }

    // No new frame this tick is NORMAL (a window arrives at up to ~30 fps vs the UI timer). Hold the last
    // spectrum; only fade once genuinely starved (audio stopped ~0.5 s). Bounded counter.
    void starve() noexcept
    {
        fallMul = fallMultiplier();
        if (starveTicks < 16) ++starveTicks;
        for (int k = 0; k < numTiers; ++k)
        {
            const Tier& t = tiers[(std::size_t) k];
            float* pw = specPow.data()     + t.binOffset;
            float* pk = specPeakPow.data() + t.binOffset;
            for (int i = 0; i < t.bins; ++i)
            {
                if (starveTicks > 15)
                {
                    pw[i] *= 0.5f;                                       // −3 dB/tick (60 dB in 0.7 s, as the classic pane's
                    if (pw[i] < kFlushPower) pw[i] = 0.0f;               // fade) — and faster than the peak falls, so the
                }                                                        // peak trace never sinks under the fill.
                pk[i] = clampPeak (pk[i] * fallMul);                     // ×10^(−fall/10) IS −fall dB, resting at the floor
            }
            rebuildSums (t);
        }
    }

    //==============================================================================
    // Reads (message thread). f in Hz, fs the sample rate the frame was captured at.

    // Seam of tier k (k ≥ 1): the frequency from which tier k's bin is fine enough for the band. Tier 0
    // (the longest) has no seam — it is the fallback below every other tier's.
    double seamHz (int k, double fs) const noexcept
    {
        if (k <= 0 || k >= numTiers || ! (fs > 0.0)) return 0.0;
        return std::max (1e-3, binsPerBand) * (fs / (double) tiers[(std::size_t) k].size) / bandWidthFactor();
    }

    // The tier a reading at f is taken from (before the seam blend): the shortest whose seam f has passed.
    // A tier that has not seen a frame yet (a frame shorter than it, before the first full one) is
    // skipped, so the display falls back to what IS there; before any frame at all, plain geometry.
    int tierAt (double f, double fs) const noexcept
    {
        bool anyValid = false;
        for (int j = 0; j < numTiers; ++j) anyValid = anyValid || tiers[(std::size_t) j].valid;
        int k = -1;
        for (int j = 0; j < numTiers; ++j)
        {
            if (anyValid && ! tiers[(std::size_t) j].valid) continue;
            if (k < 0 || f >= seamHz (j, fs)) k = j; else break;
        }
        return k < 0 ? 0 : k;
    }

    // The blend above seam k is at most blendOctaves wide and never reaches the next seam.
    double blendWidthOctaves (int k, double fs) const noexcept
    {
        double w = blendOctaves;
        if (k + 1 < numTiers && seamHz (k, fs) > 0.0) w = std::min (w, std::log2 (seamHz (k + 1, fs) / seamHz (k, fs)));
        return w;
    }

    // One tier's band-power reading at f (no blend) — for tests and diagnostics.
    double tierBandDb     (int k, double f, double fs) const noexcept { return (k >= 0 && k < numTiers) ? toDb (bandPower (tiers[(std::size_t) k], sumDb.data(),   f, fs)) : (double) kFloorDb; }
    double tierBandPeakDb (int k, double f, double fs) const noexcept { return (k >= 0 && k < numTiers) ? toDb (bandPower (tiers[(std::size_t) k], sumPeak.data(), f, fs)) : (double) kFloorDb; }

    // The stitched readings — what the display shows before the tilt. Floored at kFloorDb.
    double readDb     (double f, double fs) const noexcept { return toDb (readPower (sumDb.data(),   f, fs)); }
    double readPeakDb (double f, double fs) const noexcept { return toDb (readPower (sumPeak.data(), f, fs)); }

    // The stitched readings WITH the display tilt (dB/oct about tiltPivotHz). The tilt is a gain on the
    // signal, so it is applied to the POWER and the floor comes after: silence stays at the floor whatever
    // the tilt, and a band just above the floor is lifted like any other. (Adding the tilt to a floored
    // dB value drew silence as a straight line rising at the tilt's slope — the floor "lying on the plot
    // and tilted", seen at once on a 120 dB range.) Past Nyquist there is nothing to read and nothing is
    // read: the floor. (Holding the last band flat to the end of the axis drew a shelf there.)
    double readDb     (double f, double fs, double tiltDbPerOct, double tiltPivotHz) const noexcept { return toDb (readPower (sumDb.data(),   f, fs) * tiltGain (f, fs, tiltDbPerOct, tiltPivotHz)); }
    double readPeakDb (double f, double fs, double tiltDbPerOct, double tiltPivotHz) const noexcept { return toDb (readPower (sumPeak.data(), f, fs) * tiltGain (f, fs, tiltDbPerOct, tiltPivotHz)); }

    // Sample the stitched spectrum into N+1 log-frequency columns and emit plot points through the map —
    // the same contract as SpectrumPane::buildColumns: emit (int i, float x, float yFill, float yPeak),
    // i = 0..N, x ascending from 0. Every column is specDbToY (readDb (f, fs, tilt, pivot)).
    template <class Emit>
    void buildColumns (const PlotMap& pm, double fs, double tiltDbPerOct, double tiltPivotHz, Emit&& emit) const
    {
        const int N = std::clamp ((int) pm.width, 256, kMaxColumns - 1);
        ensurePlan (pm, fs, tiltDbPerOct, tiltPivotHz, N);
        const double* S = sumDb.data();
        const double* P = sumPeak.data();
        for (int i = 0; i <= N; ++i)
        {
            const Column& c = plan[(std::size_t) i];
            emit (i, c.x, pm.specDbToY (toDb (planPower (c, S) * c.tilt)),
                          pm.specDbToY (toDb (planPower (c, P) * c.tilt)));
        }
    }

private:
    struct Tier
    {
        int    order = 0, size = 0, bins = 0;
        int    binOffset = 0, winOffset = 0, sumOffset = 0;
        double powerNorm = 0.0, enbwBins = 1.0;
        bool   valid = false;                                            // ingested at least once since reset
    };

    bool inRange (int k, int i) const noexcept { return k >= 0 && k < numTiers && i >= 0 && i < tiers[(std::size_t) k].bins; }

    // "Fall peakFallDb per tick" as a power multiplier. One pow per ingest/starve, not per bin. A
    // non-finite setting must not poison the hold (the sibling lets max(NaN, …) through).
    float fallMultiplier() const noexcept
    {
        return std::isfinite (peakFallDb) ? (float) std::pow (10.0, -(double) peakFallDb / 10.0) : 1.0f;
    }

    // The hold's resting place is the floor, and its ceiling is the same one a bin power gets. The upper
    // clamp also bounds a NEGATIVE peakFallDb, which in the dB domain merely grew without limit but in
    // power would reach an infinity in a few hundred ticks. NaN lands on the floor.
    static float clampPeak (float v) noexcept
    {
        if (! (v >= kFloorPowerF)) return kFloorPowerF;
        return v > (float) kMaxPower ? (float) kMaxPower : v;
    }

    //==============================================================================
    // THE COLUMN PLAN. Everything buildColumns needs that does not depend on the frame: which tier owns
    // the column, its band's fractional bin edges resolved to ABSOLUTE prefix indices, the blend weight
    // and the display tilt. A reading is then four prefix loads per contributing tier. Rebuilt only when
    // one of its inputs changes — a resize, a new sample rate, a tuning edit, a tier becoming valid, a
    // setTiers — so a steady display pays for it once, not every tick.
    struct Edge                                                          // one tier's band as a fixed form
    {
        int    bLo = 0, bHi = 0;                                         // absolute indices into sumDb / sumPeak
        double frLo = 0.0, frHi = 0.0;                                   // fractional part of each edge
        double span = 0.0;                                               // the band's width in bins; 0 = contributes nothing
        double cell = 0.0;                                               // >0 only when bin-limited (span < 1), and then it is the cell
        // span and cell are kept APART, not folded into one factor, so the bin-limited branch evaluates
        // p / span * cell in the sibling's own order — folding it to p * (cell/span) rounds differently
        // and the fill trace would stop being bit-identical for no gain.
    };

    struct Column
    {
        Edge   hi {}, lo {};                                             // hi = the tier tierAt picked; lo = the longer one
        double w     = 1.0;                                              // weight on hi (blend only)
        double tilt  = 1.0;                                              // display tilt as a power gain
        float  x     = 0.0f;
        bool   blend = false;
    };

    struct PlanKey
    {
        int      n = -1, nTiers = -1;
        unsigned validMask = 0, epoch = 0;
        float    width = 0.0f;
        double   fMin = 0.0, fMax = 0.0, fs = 0.0, tilt = 0.0, pivot = 0.0;
        double   bandOct = 0.0, blendOct = 0.0, binsPer = 0.0;
        bool operator== (const PlanKey&) const = default;
    };

    // bandPower()'s edges, resolved once. Mirrors it line for line — including the bin-limited rule and
    // the Nyquist clamp — so the plan is the same arithmetic, only hoisted.
    Edge bandEdge (const Tier& t, double f, double fs) const noexcept
    {
        Edge e;
        if (! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs) || f > 0.5 * fs) return e;
        const double binHz = fs / (double) t.size;
        const double h     = bandHalfFactor();
        const double uMax  = (double) t.bins - 0.5;
        const double uLo   = std::clamp ((f / h) / binHz + 0.5, 0.5, uMax);
        const double uHi   = std::clamp ((f * h) / binHz + 0.5, 0.5, uMax);
        const double span  = uHi - uLo;
        if (! (span > 0.0)) return e;                                    // span stays 0
        const int bLo = (int) uLo, bHi = (int) uHi;                      // uLo,uHi ∈ [0.5, bins−0.5] → b ∈ [0, bins−1]
        e.bLo = t.sumOffset + bLo; e.frLo = uLo - (double) bLo;
        e.bHi = t.sumOffset + bHi; e.frHi = uHi - (double) bHi;
        e.span = span;
        if (span < 1.0)
        {
            const int centre = (int) (f / binHz + 0.5);
            e.cell = (centre <= 0 || centre >= t.bins - 1) ? 0.5 : 1.0;
        }
        return e;
    }

    static double edgePower (const Edge& e, const double* S) noexcept
    {
        if (! (e.span > 0.0)) return 0.0;
        const double a = S[e.bLo] + e.frLo * (S[e.bLo + 1] - S[e.bLo]);
        const double b = S[e.bHi] + e.frHi * (S[e.bHi + 1] - S[e.bHi]);
        double p = b - a;
        if (e.cell > 0.0) p = p / e.span * e.cell;                       // the bin-limited rule, in the sibling's order
        return p > 0.0 ? p : 0.0;
    }

    static double planPower (const Column& c, const double* S) noexcept
    {
        const double p = edgePower (c.hi, S);
        return c.blend ? (1.0 - c.w) * edgePower (c.lo, S) + c.w * p : p;
    }

    void ensurePlan (const PlotMap& pm, double fs, double tiltDbPerOct, double tiltPivotHz, int n) const noexcept
    {
        PlanKey key;
        key.n = n; key.width = pm.width; key.fMin = pm.freqMin; key.fMax = pm.freqMax;
        key.fs = fs; key.tilt = tiltDbPerOct; key.pivot = tiltPivotHz;
        key.bandOct = bandOctaves; key.blendOct = blendOctaves; key.binsPer = binsPerBand;
        key.nTiers = numTiers; key.epoch = layoutEpoch;
        for (int k = 0; k < numTiers; ++k) if (tiers[(std::size_t) k].valid) key.validMask |= 1u << (unsigned) k;
        if (key == planKey) return;                                      // the common case: nothing moved
        planKey = key;

        for (int i = 0; i <= n; ++i)
        {
            Column& c = plan[(std::size_t) i];
            c = Column {};
            c.x = (float) i / (float) n * pm.width;
            const double f = pm.xToFreq (c.x);
            c.tilt = tiltGain (f, fs, tiltDbPerOct, tiltPivotHz);
            // readPower()'s guards: past Nyquist, or nothing to read, the column IS the floor — both
            // edges keep span 0, so planPower returns 0 and toDb floors it.
            if (numTiers == 0 || ! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs) || f > 0.5 * fs)
                continue;
            const int k = tierAt (f, fs);
            c.hi = bandEdge (tiers[(std::size_t) k], f, fs);
            if (k > 0 && blendOctaves > 0.0 && tiers[(std::size_t) (k - 1)].valid)
            {
                const double width = blendWidthOctaves (k, fs);
                const double w     = width > 0.0 ? std::log2 (f / seamHz (k, fs)) / width : 1.0;
                if (w < 1.0) { c.blend = true; c.w = w; c.lo = bandEdge (tiers[(std::size_t) (k - 1)], f, fs); }
            }
        }
    }

    static double toDb (double p) noexcept
    {
        return (p > kFloorPower) ? 10.0 * std::log10 (p) : (double) kFloorDb;
    }

    // The display tilt as a power gain; a non-positive / non-finite pivot means no tilt.
    static double tiltGain (double f, double fs, double tiltDbPerOct, double tiltPivotHz) noexcept
    {
        if (! (tiltPivotHz > 0.0) || ! std::isfinite (tiltPivotHz) || ! std::isfinite (tiltDbPerOct) || ! (f > 0.0) || ! (fs > 0.0)) return 1.0;
        return std::pow (10.0, tiltDbPerOct * std::log2 (f / tiltPivotHz) / 10.0);
    }

    double bandHalfFactor()  const noexcept { return std::exp2 (0.5 * std::max (1e-4, bandOctaves)); }   // 2^(o/2)
    double bandWidthFactor() const noexcept { const double h = bandHalfFactor(); return h - 1.0 / h; } // B(f) = f · (2^(o/2) − 2^(−o/2))

    // Tier k on the frame of n samples: `subs` Hann-windowed suffixes at 50 % hop, ending at the frame's end
    // and stepping back (the newest first), power averaged per bin, then the per-bin smoothing + peak-hold.
    void transformTier (int k, int n, int subs, bool seed) noexcept
    {
        const Tier&  t = tiers[(std::size_t) k];
        const float* w = window.data() + t.winOffset;
        const int    half = t.size / 2;
        double* a = acc.data();
        for (int i = 0; i < t.bins; ++i) a[i] = 0.0;

        for (int s = 0; s < subs; ++s)
        {
            const float* src = frame.data() + (n - t.size - s * half);
            for (int i = 0; i < t.size; ++i)                             // a NaN/Inf SAMPLE would poison every bin of the
            {                                                            // transform — it is dropped here, the rest survive
                const float v = src[i];
                work[(std::size_t) i] = std::isfinite (v) ? v * w[i] : 0.0f;
            }
            fft[(std::size_t) k].forward (work.data(), spec.data());   // real[N] -> [DC, Nyq, re1, im1, …]
            // DC and Nyquist are the only real-only cells; peeling them leaves a contiguous interior
            // whose bounds the compiler can see. Bin for bin the arithmetic is what the sibling does.
            {
                const double dc = spec[0], nq = spec[1];
                a[0]    += dc * dc;
                a[half] += nq * nq;
            }
            for (int i = 1; i < half; ++i)
            {
                const double re = spec[(std::size_t) (2 * i)], im = spec[(std::size_t) (2 * i + 1)];
                a[i] += re * re + im * im;
            }
        }

        const double norm = t.powerNorm / (double) subs;
        float* pw = specPow.data()     + t.binOffset;
        float* pk = specPeakPow.data() + t.binOffset;
        for (int i = 0; i < t.bins; ++i)
        {
            double p = a[i] * norm;
            if (! (p >= 0.0)) p = 0.0;                                   // NaN → silence, never poison
            if (p > kMaxPower) p = kMaxPower;                            // …and nothing narrows to a float infinity
            if (seed) pw[i] = (float) p;
            else      pw[i] += smoothCoeff * ((float) p - pw[i]);
            if (std::fabs (pw[i]) < kFlushPower) pw[i] = 0.0f;            // never leave the state in the subnormals
            // The sibling compares the hold against toDb(pw), which floors at −200; here the same floor
            // is kFloorPower, and the comparison happens in power. Same trace, no transcendental.
            const float floored = pw[i] > kFloorPowerF ? pw[i] : kFloorPowerF;
            pk[i] = seed ? floored : clampPeak (std::max (pk[i] * fallMul, floored));
        }
        rebuildSums (t);
    }

    // Prefix sums of per-bin POWER in double: S[0] = 0, S[i+1] = S[i] + P_i. A band's power is then two
    // lookups whatever its width, and a fractional bin edge is a linear fraction of P_i. Every positive
    // bin counts, however small — eight components at −123 dB ARE a −114 dB band — and only the final
    // reading is floored. The peak trace integrates the per-bin holds (a persistence envelope); a hold
    // resting on the floor counts as nothing.
    void rebuildSums (const Tier& t) noexcept
    {
        const float* pw = specPow.data()     + t.binOffset;
        const float* pk = specPeakPow.data() + t.binOffset;
        double* S = sumDb.data()   + t.sumOffset;
        double* P = sumPeak.data() + t.sumOffset;
        S[0] = 0.0; P[0] = 0.0;
        for (int i = 0; i < t.bins; ++i)
        {
            const double p = (double) pw[i];
            S[i + 1] = S[i] + (p > 0.0 ? p : 0.0);
            P[i + 1] = P[i] + (pk[i] > kFloorPowerF ? (double) pk[i] : 0.0);
        }
    }

    static double prefixAt (const double* S, int bins, double u) noexcept   // ∫ density from 0 to u (bin units)
    {
        if (u <= 0.0) return 0.0;
        if (u >= (double) bins) return S[bins];
        const int    b  = (int) u;
        const double fr = u - (double) b;
        return S[b] + fr * (S[b + 1] - S[b]);
    }

    // Band power at f from one tier's prefix sums (sums = the flat array; the tier's offset is applied
    // here). Linear power in "full-scale sine = 1" units; 0 when the band is empty or invalid.
    // Bin i's cell is [(i − ½)·binHz, (i + ½)·binHz]; in bin units u = f/binHz + ½ a cell spans [i, i+1).
    // The one-sided axis runs u ∈ [½, bins − ½]: DC and Nyquist are half cells, which IS their one-sided
    // weight. A band reaching past Nyquist is clipped to it (the reading is what lies below Nyquist);
    // f itself above Nyquist is nothing — read() returns the floor there.
    double bandPower (const Tier& t, const double* sums, double f, double fs) const noexcept
    {
        if (! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs) || f > 0.5 * fs) return 0.0;
        const double* S     = sums + t.sumOffset;
        const double  binHz = fs / (double) t.size;
        const double  h     = bandHalfFactor();
        const double  uMax  = (double) t.bins - 0.5;
        const double uLo  = std::clamp ((f / h) / binHz + 0.5, 0.5, uMax);
        const double uHi  = std::clamp ((f * h) / binHz + 0.5, 0.5, uMax);
        const double span = uHi - uLo;                                   // the band's width in bins (edge-clamped)
        if (! (span > 0.0)) return 0.0;
        double p = prefixAt (S, t.bins, uHi) - prefixAt (S, t.bins, uLo);
        if (span < 1.0)                                                  // bin-limited: the local density × one cell
        {
            const int centre = (int) (f / binHz + 0.5);
            const double cell = (centre <= 0 || centre >= t.bins - 1) ? 0.5 : 1.0;   // DC / Nyquist cells are half-width
            p = p / span * cell;
        }
        return (p > 0.0) ? p : 0.0;
    }

    // The stitched reading as linear power (0 = nothing); toDb() floors it.
    double readPower (const double* sums, double f, double fs) const noexcept
    {
        if (numTiers == 0 || ! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs)) return 0.0;
        if (f > 0.5 * fs) return 0.0;                                    // above Nyquist there is nothing: the floor
        const int k = tierAt (f, fs);
        double p = bandPower (tiers[(std::size_t) k], sums, f, fs);
        if (k > 0 && blendOctaves > 0.0 && tiers[(std::size_t) (k - 1)].valid)
        {
            const double width = blendWidthOctaves (k, fs);
            const double w = width > 0.0 ? std::log2 (f / seamHz (k, fs)) / width : 1.0;   // 0 at the seam → 1 a blend width above it
            if (w < 1.0)                                                                    // power crossfade: contributions add in power
                p = (1.0 - w) * bandPower (tiers[(std::size_t) (k - 1)], sums, f, fs) + w * p;
        }
        return p;
    }

    std::array<Tier, (std::size_t) MaxTiers>          tiers {};
    int                                               numTiers = 0;
    Fft                                               fft[(std::size_t) MaxTiers];

    std::array<float, (std::size_t) kMaxSize>         frame {};         // the consumer's frame (frameInput); read-only to the pane
    std::array<float, (std::size_t) kMaxSize>         work {};          // one sub-window, windowed
    std::array<float, (std::size_t) kMaxSize>         spec {};          // packed real-FFT output
    std::array<double, (std::size_t) kMaxBins>        acc {};           // one tier's power, summed over its sub-windows
    std::array<float, (std::size_t) kWinCapacity>     window {};        // Hann per tier, flat-packed
    std::array<float, (std::size_t) kBinCapacity>     specPow {};       // per-bin smoothed POWER, flat-packed by tier
    std::array<float, (std::size_t) kBinCapacity>     specPeakPow {};   // per-bin peak-hold, in POWER (the sibling keeps dB)
    std::array<double, (std::size_t) kSumCapacity>    sumDb {};         // prefix sums of bin power (bins+1 per tier)
    std::array<double, (std::size_t) kSumCapacity>    sumPeak {};       // prefix sums of the peak-hold power
    int      starveTicks = 0;
    bool     frameArmed  = false;                                        // debug protocol guard (see frameInput)
    float    fallMul     = 1.0f;                                         // 10^(−peakFallDb/10), refreshed per ingest/starve
    unsigned layoutEpoch = 0;                                            // bumped by setTiers: the prefix offsets moved

    // The plan is display state, not analysis state: buildColumns is const, so it is mutable. Fixed size,
    // so a rebuild allocates nothing (~80 KB beside the pane's own ~0.8 MB).
    mutable std::array<Column, (std::size_t) kMaxColumns> plan {};
    mutable PlanKey                                       planKey {};
};

// Drop-in default: tiers from the tap's 16384-sample frame (16384 / 4096 / 1024).
using MultiResSpectrumPaneFast = MultiResSpectrumPaneFastT<>;

} // namespace felitronics::analysis
