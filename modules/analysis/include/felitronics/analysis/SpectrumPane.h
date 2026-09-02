// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Fft.h>
#include "RollingSpectrumTap.h"                        // kMaxOrder / kMaxSize — the max analysis window
#include <felitronics/analysis/PlotMap.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

//==============================================================================
// SpectrumPane — the analyzer's JUCE-free pipeline + state. Incubated in TabbyEQ, moved home
// once OrbitAmp became the second consumer.
// Pull-fed frames go Hann window → real FFT → per-bin dB with smoothing + a slow-decay peak-hold;
// buildColumns() then samples the bins into log-frequency plot columns — the "liquid" rule: where a
// column spans LESS than one bin, interpolate between neighbours (smooth lows), where it spans
// several, take the bin PEAK (narrow HF resonances stay visible) — applies the pink-noise display
// tilt and maps through PlotMap. The owner draws the emitted points (path building stays JUCE-side)
// and owns the frame SOURCE: fill frameInput() and call ingest(order), or starve() when no frame
// arrived this tick. Message-thread only; the audio thread is never touched (frames arrive via a
// lock-free tap).
//
// RESOLUTION: the analysis size N is chosen at RUNTIME (1024 … 16384 = FFT order 10..kMaxOrder, which
// rides RollingSpectrumTap's). All FFT plans are prebuilt in the ctor, so switching resolution costs no allocation — a
// live switch just starts FFTing at the new N. On the frame where the order changes, the dB state is
// SEEDED directly from the new bins (not smoothed up from the −120 floor) so the display cuts cleanly
// to the new resolution instead of fading in. Storage is sized to the maximum window; a smaller N
// uses the leading portion.
namespace felitronics::analysis
{

// The FFT backend is a template parameter constrained to the packed-Hermitian layout the bin loop reads
// ([DC, Nyq, re1, im1, …]): the scalar reference by default, or fftpffft::PffftOrderedRealFft for SIMD.
// `SpectrumPane` (the alias below) is the scalar one and reads exactly as before.
template <felitronics::core::fft::PackedHermitianSpectrum Fft = felitronics::core::fft::DefaultRealFft>
struct SpectrumPaneT
{
    using FftType = Fft;
    static constexpr int kMaxOrder = felitronics::analysis::RollingSpectrumTap::kMaxOrder;   // 14
    static constexpr int kMaxSize  = felitronics::analysis::RollingSpectrumTap::kMaxSize;    // 16384
    static constexpr int kMaxBins  = kMaxSize / 2 + 1;
    static constexpr int kMinOrder = 10;                                                     // 1024
    // The internal floor. It sits far below any plot bottom on purpose: the display tilt is ADDED to a
    // bin's dB, so a bin clamped at the floor and then tilted would surface as a straight line rising at
    // the tilt's slope (with a −120 floor and +6 dB/oct that line stood at −94 dB at 20 kHz — silence
    // "lying on the plot and tilted"). At −200 the tilted floor stays below every range; a bin that has
    // real energy at −130 keeps it and is lifted like any other.
    static constexpr float kFloorDb  = -200.0f;

    float peakFallDb  = 0.8f;         // peak-hold decay per tick (~24 dB/s at 30 Hz)
    float smoothCoeff = 0.25f;        // per-tick smoothing toward the new frame (the analyzer "speed":
                                      // smaller = slower attack/release, larger = snappier)

    SpectrumPaneT()
    {
        for (int o = kMinOrder; o <= kMaxOrder; ++o)
            fft[(std::size_t) (o - kMinOrder)].prepare (1 << o);   // plan/alloc once, on the owner's (UI) thread
        buildWindow (order);          // Hann for the default order
        specDb.fill (kFloorDb);
        specPeak.fill (kFloorDb);
    }

    // The active analysis size (message-thread reads for buildColumns geometry).
    int  activeOrder()   const noexcept { return order; }
    int  activeFftSize() const noexcept { return fftSize; }

    // Fill with up to (1<<order) samples (the tap frame), then call ingest(order). PROTOCOL: exactly
    // one ingest() per fill — the window is applied IN PLACE, so ingesting the same buffer twice would
    // Hann² the data. The debug guard trips on the first violation. The buffer is max-sized, so it
    // holds any order's frame; ingest(order) says how many of those samples are live.
    float* frameInput() noexcept { frameArmed = true; return fftBuf.data(); }

    // Transform the frame at FFT `newOrder` (10..kMaxOrder). When the order changes vs the last ingest, the
    // window is rebuilt and the dB/peak state is seeded straight from the new bins (a clean cut, no
    // fade-in); otherwise the usual per-tick smoothing + peak-hold runs.
    void ingest (int newOrder) noexcept
    {
        assert (frameArmed && "SpectrumPane: fill frameInput() before every ingest()");
        frameArmed = false;
        starveTicks = 0;

        if (newOrder < kMinOrder) newOrder = kMinOrder;
        if (newOrder > kMaxOrder) newOrder = kMaxOrder;
        const bool orderChanged = (newOrder != order) || seedNext;   // reset() asks for a seed as an order change would
        seedNext = false;
        if (orderChanged)
        {
            order   = newOrder;
            fftSize = 1 << newOrder;
            numBins = fftSize / 2 + 1;
            buildWindow (newOrder);
        }

        for (int i = 0; i < fftSize; ++i) fftBuf[(std::size_t) i] *= window[(std::size_t) i];
        fft[(std::size_t) (order - kMinOrder)].forward (fftBuf.data(), spec.data());   // real[N] -> [DC, Nyq, re1,im1, …]

        const double norm = (double) fftSize * 0.25;              // Hann single-bin compensation
        for (int i = 0; i < numBins; ++i)
        {
            const float re  = (i == 0) ? spec[0] : (i == fftSize / 2) ? spec[1] : spec[(size_t) (2 * i)];
            const float im  = (i == 0 || i == fftSize / 2) ? 0.0f : spec[(size_t) (2 * i + 1)];
            const float mag = std::sqrt (re * re + im * im);      // |bin| of the unnormalized forward
            const double g  = (double) mag / norm;                // == juce::Decibels::gainToDecibels (g, kFloorDb)
            const double db = g > 0.0 ? std::max ((double) kFloorDb, 20.0 * std::log10 (g)) : (double) kFloorDb;
            if (orderChanged)                                     // seed the new-resolution bins directly (no fade-in)
            {
                specDb[(size_t) i]   = (float) db;
                specPeak[(size_t) i] = (float) db;
            }
            else
            {
                specDb[(size_t) i]  += smoothCoeff * ((float) db - specDb[(size_t) i]);          // smooth toward target
                specPeak[(size_t) i] = std::max (specPeak[(size_t) i] - peakFallDb, specDb[(size_t) i]);
            }
        }
    }

    // Forget every frame: bins to the floor, and the next ingest SEEDS from its bins (a clean cut, no
    // fade-in from −120) — for a consumer that returns to this pane after drawing another, or a stream
    // discontinuity. Keeps the order, the plans and the tuning.
    void reset() noexcept
    {
        specDb.fill (kFloorDb);
        specPeak.fill (kFloorDb);
        starveTicks = 0;
        frameArmed  = false;
        seedNext    = true;
    }

    // No new frame this tick is NORMAL (a window arrives at up to ~30 fps vs the 30 Hz timer). Hold the
    // last spectrum; only fade once genuinely starved (audio stopped ~0.5 s). Bounded to the active bins.
    void starve() noexcept
    {
        if (starveTicks < 16) ++starveTicks;                      // bounded — never overflows over long silence
        for (int i = 0; i < numBins; ++i)
        {
            if (starveTicks > 15) specDb[(size_t) i] += 0.05f * (kFloorDb - specDb[(size_t) i]);
            specPeak[(size_t) i] = std::max (kFloorDb, specPeak[(size_t) i] - peakFallDb);
        }
    }

    // Sample the bins into N+1 log-frequency columns and emit plot points through the map.
    // emit (int i, float x, float yFill, float yPeak) is called with i = 0..N, x ascending.
    template <class Emit>
    void buildColumns (const PlotMap& pm, double fs, double tiltDbPerOct, double tiltPivotHz, Emit&& emit) const
    {
        const double binPerHz = (double) fftSize / fs;
        const int    N = std::clamp ((int) pm.width, 256, 900);

        // The FIRST column is a zero-width point at freqMin — seed prevBin from ITS bin, not from
        // DC: a prevBin = 0 seed made column 0 peak-scan (0, curBin] and leak DC/bin-1 energy into
        // the left edge whenever freqMin sits above bin 1 (low sample rates, zoomed-in maps).
        // Found by the crew's theory suite (codex #7); on the classic 20 Hz axis at ≥44.1 kHz both
        // seeds are bin 0, so the fix is bit-identical there.
        int prevBin = std::clamp ((int) std::floor (pm.xToFreq (0.0f) * binPerHz), 0, numBins - 1);
        const double nyquist = 0.5 * fs;
        for (int i = 0; i <= N; ++i)
        {
            const float  x = (float) i / (float) N * pm.width;
            const double f = pm.xToFreq (x);
            const int    curBin = std::clamp ((int) std::floor (f * binPerHz), 0, numBins - 1);
            // Above Nyquist there is nothing: the floor (the columns used to repeat the Nyquist bin to the
            // end of the axis — a shelf). The FIRST column past it still runs the peak scan, so the bins
            // between the last sub-Nyquist column and the Nyquist bin are drawn once, never skipped.
            if (f > nyquist && prevBin >= numBins - 1)
            {
                emit (i, x, pm.specDbToY ((double) kFloorDb), pm.specDbToY ((double) kFloorDb));
                continue;
            }
            const float  tilt = (float) (tiltDbPerOct * std::log2 (f / tiltPivotHz));   // pink-noise comp
            emit (i, x, pm.specDbToY (column (specDb,   f, binPerHz, prevBin, curBin) + tilt),
                        pm.specDbToY (column (specPeak, f, binPerHz, prevBin, curBin) + tilt));
            prevBin = curBin;
        }
    }

private:
    template <class Arr>
    float column (const Arr& raw, double f, double binPerHz, int loBin, int hiBin) const noexcept
    {
        if (hiBin <= loBin)                                       // < 1 bin this column -> interpolate
        {
            const double bf = std::clamp (f * binPerHz, 0.0, (double) (numBins - 2));
            const int b0 = (int) bf; const float t = (float) (bf - (double) b0);
            return raw[(size_t) b0] + t * (raw[(size_t) (b0 + 1)] - raw[(size_t) b0]);
        }
        float m = -200.0f;                                        // >= 1 bin -> peak (max) = detail
        for (int b = std::max (0, loBin + 1); b <= std::min (numBins - 1, hiBin); ++b)
            m = std::max (m, raw[(size_t) b]);
        return m;
    }

    void buildWindow (int ord) noexcept
    {
        const int    n  = 1 << ord;
        constexpr float pi = 3.14159265358979323846f;
        for (int i = 0; i < n; ++i)
            window[(size_t) i] = 0.5f - 0.5f * std::cos (2.0f * pi * (float) i / (float) (n - 1));
    }

    // Prebuilt real-FFT plans for orders kMinOrder..kMaxOrder (index = order − kMinOrder); switching
    // resolution picks a plan — never allocates.
    Fft fft[(std::size_t) (kMaxOrder - kMinOrder + 1)];

    int order   = 11;                                             // active FFT order (default 2048)
    int fftSize = 1 << 11;
    int numBins = (1 << 11) / 2 + 1;

    std::array<float, kMaxSize>     window {};                    // Hann for the active order (leading fftSize used)
    std::array<float, kMaxSize * 2> fftBuf {};                    // frame input + FFT workspace
    std::array<float, kMaxSize>     spec {};                      // packed real-FFT output
    std::array<float, kMaxBins>     specDb {};
    std::array<float, kMaxBins>     specPeak {};                  // slow-decay peak-hold
    int starveTicks = 0;                                          // consecutive ticks with no new frame
    bool frameArmed = false;                                      // debug protocol guard (see frameInput)
    bool seedNext   = false;                                      // set by reset(): the next ingest seeds
};

// The pane every consumer has been using: the scalar reference FFT.
using SpectrumPane = SpectrumPaneT<>;

} // namespace felitronics::analysis
