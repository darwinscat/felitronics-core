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
// MultiResSpectrumPane — a constant-Q analyzer built from several FFT lengths at once. The design and
// the physics behind its level convention are in docs/ANALYZER-MULTIRES.md; the short version:
//
//   • ONE frame (the most recent 1<<frameOrder() samples, as RollingSpectrumTap publishes them) feeds
//     every TIER: tier k Hann-windows and transforms suffixes of it. A short window is a suffix of the
//     long one, so all tiers share the frame's END — the short tier reports a transient first, the
//     long tier confirms it later. The lows come from the long window (fine bins), the highs from the
//     short one (a short time window), stitched by frequency. The frame is never modified: each tier
//     copies its suffix into scratch and windows the copy.
//   • A short tier covers only part of the interval between two frames (1024 samples of a ~1600-sample
//     hop at 30 fps), so a transient could fall in the gap. Tell the pane the hop (coverSamples) and a
//     tier shorter than it analyses several 50 %-overlapping sub-windows reaching back over the hop
//     and averages their power (Welch) — no blind gap, and less variance for free. The longest tier
//     is one window, as before.
//   • A log axis wants constant Q, and a seam between two FFT lengths is only invisible if both tiers
//     read the SAME number for the same signal. No fixed per-bin scalar can (noise power per bin is
//     proportional to bin width — a 4× shorter window sits 6 dB higher, bin for bin; a max over the
//     column's bins is not invariant either). A bandwidth-integrated quantity is. So a reading is the
//     POWER IN A FRACTIONAL-OCTAVE BAND (bandOctaves, default 1/24), integrated over the tier's bins
//     with fractional edge overlap (double prefix sums — whole-bin membership would sawtooth 1↔2 bins
//     between neighbouring columns). Per interior bin the power is 4·|X_k|²/(N·Σw²): coherent-gain
//     amplitude divided by the window's equivalent noise bandwidth in bins (1.5·N/(N−1) for this Hann,
//     measured from the window as built). DC and Nyquist are half-width cells (one-sided weight ½).
//     Calibration: a full-scale sine whose main lobe lies inside the band reads ≈ 0 dBFS; white noise
//     of RMS σ in a band B reads 10·log10(4σ²B/fs) — the level of the sine that would carry the band's
//     power, i.e. 2× the band's mean square; disjoint bands add up to the windowed frame's power
//     (Parseval). A signal can read above 0 dBFS (a full-scale square's fundamental is +2.1 dB).
//   • Where the band is narrower than a bin (the lows of the longest tier) the band cannot resolve:
//     the reading is the LOCAL BIN's power (the density the band sees × one cell), continuous at
//     B = B_bin. A bin-centred sine there reads its peak bin alone, −1.76 dB (−3.2 dB half-way between
//     bins, Hann scallop); noise reads flat there and rises +3 dB/oct (white) / flat (pink) once bands
//     are wider than a bin. That knee is the true resolution limit of the longest window, shown.
//   • The tier used at f is the SHORTEST whose bin is fine enough: B_bin ≤ B(f)/binsPerBand. Above
//     each seam the reading crossfades — in POWER, the domain where band contributions add — from the
//     longer tier to the shorter over blendOctaves, so a sine near a narrow band (its lobe partly
//     outside; scallop) never steps.
//
// Per bin, each tier smooths POWER with a one-pole per tick (smoothCoeff) and holds a PEAK in dB that
// falls peakFallDb per tick (SpectrumPane's peak law). Smoothing power, not dB, is what keeps the seam
// honest: a log-domain average is a geometric mean whose bias depends on how many windows fed it, and
// the Welch tiers are fed by several. It also means silence fades as a release (−1.25 dB/tick at 0.25)
// rather than collapsing. The peak trace integrates per-bin holds, so on a band it is a persistence
// envelope, not the power of any one instant. starve() = hold, then fade after ~0.5 s of no frames.
//
// Message-thread only; the audio thread is never touched (frames arrive through a lock-free tap).
// Storage is fixed and flat-packed across tiers; ingest()/starve()/buildColumns() never allocate.
// setTiers() is the one allocating call (it prepares the FFT plans). The object is large (~0.7 MB at
// order 14) — hold it by unique_ptr or as a member of a heap-allocated owner, not on the stack.
//
// SpectrumPane stays as it is: this is a sibling for a constant-Q display, not a replacement.
namespace felitronics::analysis
{

template <int MaxOrder = RollingSpectrumTap::kMaxOrder, int MaxTiers = 4>
struct MultiResSpectrumPaneT
{
    using Fft = felitronics::core::fft::DefaultRealFft;
    static_assert (Fft::kPackedHermitianSpectrum, "MultiResSpectrumPane reads the packed [DC, Nyq, re, im…] layout; another backend needs an adapter");

    static constexpr int kMaxOrder = MaxOrder;
    static constexpr int kMinOrder = 8;                                  // 256 — below that a tier is a smear, not a window
    static constexpr int kMaxSize  = 1 << MaxOrder;
    static constexpr int kMaxBins  = kMaxSize / 2 + 1;
    static constexpr int kMaxTiers = MaxTiers;
    static_assert (MaxOrder >= kMinOrder && MaxOrder <= 24, "MultiResSpectrumPane: MaxOrder out of range");
    static_assert (MaxTiers >= 1 && MaxTiers <= MaxOrder - kMinOrder + 1, "MultiResSpectrumPane: more tiers than distinct orders");

    // Capacity of the flat-packed per-bin storage: the sum of (N/2 + 1) over any MaxTiers DISTINCT orders
    // ≤ MaxOrder is < 2^MaxOrder + MaxTiers (a geometric tail plus one per tier). Windows: Σ N < 2·kMaxSize.
    static constexpr int kBinCapacity = kMaxSize + MaxTiers;
    static constexpr int kSumCapacity = kBinCapacity + MaxTiers;          // one extra prefix entry per tier
    static constexpr int kWinCapacity = 2 * kMaxSize;

    static constexpr float  kFloorDb    = -120.0f;
    static constexpr double kFloorPower = 1.0e-12;                       // 10^(kFloorDb/10): below it a bin is silence

    //==============================================================================
    // Tuning (message thread; take effect on the next tick / read).
    float  peakFallDb   = 0.8f;          // peak-hold decay per tick (~24 dB/s at 30 Hz)
    float  smoothCoeff  = 0.25f;         // per-tick one-pole on POWER toward the new frame (the analyzer "speed")
    double bandOctaves  = 1.0 / 24.0;    // width of the constant-Q band a reading integrates
    double blendOctaves = 1.0 / 3.0;     // crossfade width above each seam (0 = hard seam)
    double binsPerBand  = 2.0;           // a tier is used where its bin is at least this many times finer than the band
    int    coverSamples = 0;             // the consumer's frame hop; a tier shorter than it Welch-averages sub-windows over it

    MultiResSpectrumPaneT()
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
                specPow [(std::size_t) (t.binOffset + i)] = 0.0f;
                specDb  [(std::size_t) (t.binOffset + i)] = kFloorDb;
                specPeak[(std::size_t) (t.binOffset + i)] = kFloorDb;
            }
            rebuildSums (t);
        }
    }

    int tierCount()         const noexcept { return numTiers; }
    int tierOrder (int k)   const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].order : 0; }
    int tierBins  (int k)   const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].bins  : 0; }
    int frameOrder()        const noexcept { return numTiers > 0 ? tiers[0].order : kMaxOrder; }   // request THIS from the tap
    int frameSize()         const noexcept { return 1 << frameOrder(); }
    double enbwBins (int k) const noexcept { return (k >= 0 && k < numTiers) ? tiers[(std::size_t) k].enbwBins : 1.0; }

    // Per-bin introspection (tests / diagnostics): the smoothed power's dB and the peak-hold of tier k, bin i.
    float tierBinDb   (int k, int i) const noexcept { return inRange (k, i) ? specDb  [(std::size_t) (tiers[(std::size_t) k].binOffset + i)] : kFloorDb; }
    float tierBinPeak (int k, int i) const noexcept { return inRange (k, i) ? specPeak[(std::size_t) (tiers[(std::size_t) k].binOffset + i)] : kFloorDb; }

    // How many sub-windows tier k analyses per frame for the current coverSamples and a frame of `order`.
    int subWindows (int k, int order) const noexcept
    {
        if (k < 0 || k >= numTiers) return 0;
        const Tier& t = tiers[(std::size_t) k];
        const int n = 1 << std::clamp (order, 0, kMaxOrder);
        if (t.size > n) return 0;
        const int hop  = t.size / 2;
        const int need = 1 + (std::max (0, coverSamples - hop) + hop - 1) / hop;   // windows to reach back over the hop
        const int room = 1 + (n - t.size) / hop;                                    // windows the frame can hold
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
        assert (frameArmed && "MultiResSpectrumPane: fill frameInput() before every ingest()");
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
        if (starveTicks < 16) ++starveTicks;
        for (int k = 0; k < numTiers; ++k)
        {
            const Tier& t = tiers[(std::size_t) k];
            float* pw = specPow.data()  + t.binOffset;
            float* db = specDb.data()   + t.binOffset;
            float* pk = specPeak.data() + t.binOffset;
            for (int i = 0; i < t.bins; ++i)
            {
                if (starveTicks > 15)
                {
                    pw[i] *= 0.9f;                                       // −0.46 dB/tick, ~14 dB/s: a fade, not a cut
                    db[i]  = (float) toDb ((double) pw[i]);
                }
                pk[i] = std::max (kFloorDb, pk[i] - peakFallDb);
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
    int tierAt (double f, double fs) const noexcept
    {
        int k = 0;
        for (int j = 1; j < numTiers; ++j)
            if (f >= seamHz (j, fs)) k = j; else break;
        return k;
    }

    // One tier's band-power reading at f (no blend) — for tests and diagnostics.
    double tierBandDb     (int k, double f, double fs) const noexcept { return (k >= 0 && k < numTiers) ? toDb (bandPower (tiers[(std::size_t) k], sumDb.data(),   f, fs)) : (double) kFloorDb; }
    double tierBandPeakDb (int k, double f, double fs) const noexcept { return (k >= 0 && k < numTiers) ? toDb (bandPower (tiers[(std::size_t) k], sumPeak.data(), f, fs)) : (double) kFloorDb; }

    // The stitched readings — what the display shows (before the tilt).
    double readDb     (double f, double fs) const noexcept { return read (sumDb.data(),   f, fs); }
    double readPeakDb (double f, double fs) const noexcept { return read (sumPeak.data(), f, fs); }

    // Sample the stitched spectrum into N+1 log-frequency columns and emit plot points through the map —
    // the same contract as SpectrumPane::buildColumns: emit (int i, float x, float yFill, float yPeak),
    // i = 0..N, x ascending from 0. The pink-noise display tilt is added around tiltPivotHz.
    template <class Emit>
    void buildColumns (const PlotMap& pm, double fs, double tiltDbPerOct, double tiltPivotHz, Emit&& emit) const
    {
        const int N = std::clamp ((int) pm.width, 256, 900);
        for (int i = 0; i <= N; ++i)
        {
            const float  x = (float) i / (float) N * pm.width;
            const double f = pm.xToFreq (x);
            const double tilt = tiltDbPerOct * std::log2 (f / tiltPivotHz);
            emit (i, x, pm.specDbToY (readDb (f, fs) + tilt), pm.specDbToY (readPeakDb (f, fs) + tilt));
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

    static double toDb (double p) noexcept
    {
        return (p > kFloorPower) ? 10.0 * std::log10 (p) : (double) kFloorDb;
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
            for (int i = 0; i < t.size; ++i) work[(std::size_t) i] = src[i] * w[i];
            fft[(std::size_t) k].forward (work.data(), spec.data());   // real[N] -> [DC, Nyq, re1, im1, …]
            for (int i = 0; i < t.bins; ++i)
            {
                const double re = (i == 0) ? spec[0] : (i == half) ? spec[1] : spec[(std::size_t) (2 * i)];
                const double im = (i == 0 || i == half) ? 0.0 : spec[(std::size_t) (2 * i + 1)];
                a[i] += re * re + im * im;
            }
        }

        const double norm = t.powerNorm / (double) subs;
        float* pw = specPow.data()  + t.binOffset;
        float* db = specDb.data()   + t.binOffset;
        float* pk = specPeak.data() + t.binOffset;
        for (int i = 0; i < t.bins; ++i)
        {
            double p = a[i] * norm;
            if (! (p >= 0.0) || ! std::isfinite (p)) p = 0.0;             // a NaN/Inf bin is silence, never poison
            if (seed) pw[i] = (float) p;
            else      pw[i] += smoothCoeff * ((float) p - pw[i]);
            db[i] = (float) toDb ((double) pw[i]);
            pk[i] = seed ? db[i] : std::max (pk[i] - peakFallDb, db[i]);
        }
        rebuildSums (t);
    }

    // Prefix sums of per-bin POWER in double: S[0] = 0, S[i+1] = S[i] + P_i. A band's power is then two
    // lookups whatever its width, and a fractional bin edge is a linear fraction of P_i. Bins at the
    // silence floor count as zero, so silence reads the floor whatever the band width. The peak trace
    // integrates the per-bin holds (a persistence envelope).
    void rebuildSums (const Tier& t) noexcept
    {
        constexpr double k = 0.23025850929940458;                        // ln(10)/10
        const float* pw = specPow.data()  + t.binOffset;
        const float* pk = specPeak.data() + t.binOffset;
        double* S = sumDb.data()   + t.sumOffset;
        double* P = sumPeak.data() + t.sumOffset;
        S[0] = 0.0; P[0] = 0.0;
        for (int i = 0; i < t.bins; ++i)
        {
            const double p = (double) pw[i];
            S[i + 1] = S[i] + (p > kFloorPower ? p : 0.0);
            P[i + 1] = P[i] + (pk[i] > kFloorDb ? std::exp ((double) pk[i] * k) : 0.0);
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
    // weight. A band reaching past Nyquist is clipped; f itself is held at the last band that fits, so
    // columns above Nyquist repeat the top reading (as SpectrumPane's do).
    double bandPower (const Tier& t, const double* sums, double f, double fs) const noexcept
    {
        if (! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs)) return 0.0;
        const double* S     = sums + t.sumOffset;
        const double  binHz = fs / (double) t.size;
        const double  h     = bandHalfFactor();
        const double  uMax  = (double) t.bins - 0.5;
        f = std::min (f, (0.5 * fs) / h);
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

    double read (const double* sums, double f, double fs) const noexcept
    {
        if (numTiers == 0 || ! (f > 0.0) || ! (fs > 0.0) || ! std::isfinite (f) || ! std::isfinite (fs)) return (double) kFloorDb;
        const int k = tierAt (f, fs);
        double p = bandPower (tiers[(std::size_t) k], sums, f, fs);
        if (k > 0 && blendOctaves > 0.0)
        {
            const double w = std::log2 (f / seamHz (k, fs)) / blendOctaves;   // 0 at the seam → 1 a blend width above it
            if (w < 1.0)                                                       // power crossfade: contributions add in power
                p = (1.0 - w) * bandPower (tiers[(std::size_t) (k - 1)], sums, f, fs) + w * p;
        }
        return toDb (p);
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
    std::array<float, (std::size_t) kBinCapacity>     specDb {};        // its dB (what the peak law compares against)
    std::array<float, (std::size_t) kBinCapacity>     specPeak {};      // per-bin peak-hold, dB
    std::array<double, (std::size_t) kSumCapacity>    sumDb {};         // prefix sums of bin power (bins+1 per tier)
    std::array<double, (std::size_t) kSumCapacity>    sumPeak {};       // prefix sums of the peak-hold power
    int  starveTicks = 0;
    bool frameArmed  = false;                                            // debug protocol guard (see frameInput)
};

// Drop-in default: tiers from the tap's 16384-sample frame (16384 / 4096 / 1024).
using MultiResSpectrumPane = MultiResSpectrumPaneT<>;

} // namespace felitronics::analysis
