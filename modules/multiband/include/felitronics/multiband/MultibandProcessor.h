// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/MultibandSplitter.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <type_traits>
#include <array>
#include <vector>

namespace felitronics::multiband
{

//==============================================================================
// felitronics::multiband::MultibandProcessor<BandProcessor, MaxBands> — split a signal into N LR4 bands
// (eq::MultibandSplitter), run ONE independent BandProcessor per band, and recombine. The generic engine
// behind a multiband compressor / saturator / width.
//
// BandProcessor concept: process(float* const* io, int nch, int n) [in place] · int latencySamples() · reset().
// prepare() is given a callable that prepares each band instance (their prepare signatures differ).
//
// Per-band bypass + solo. A global parallel dry/wet whose DRY is the band-sum ALLPASS RECONSTRUCTION (delayed
// to match the wet's latency) — NOT the raw input, which would comb against the allpass-phase wet. Bands
// are aligned to the max band latency via per-band
// DelayLines so a band with lookahead doesn't desync the sum. RT-safe: prepare() allocates, process() is
// alloc/lock/throw-free. n ≤ maxBlock.
template <class BandProcessor, int MaxBands = 4>
class MultibandProcessor
{
    static_assert (MaxBands >= 2 && MaxBands <= 8, "MultibandProcessor supports 2..8 bands");

public:
    // `maxAlignSamples` bounds the per-band latency-alignment delay (≥ the largest band lookahead in samples).
    template <class PrepareBand>
    bool prepare (double sampleRate, int maxBlock, int maxChannels, int maxAlignSamples, PrepareBand&& prepareBand)
    {
        prepared_ = false;                                   // any early return below leaves it unprepared
        if (sampleRate <= 0.0 || maxBlock < 1) return false;
        // Law 11(b): prepare() is BINDING. A silently clamped width makes process()'s refusal a fiction —
        // the caller is told nothing here and then legitimately hands over more planes than exist.
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return false;
        fs_ = sampleRate; maxBlock_ = maxBlock;
        channels_ = maxChannels;
        if (maxAlignSamples < 0) maxAlignSamples = 0;

        splitter_.prepare (sampleRate, channels_);
        bandBuf_.assign ((std::size_t) MaxBands * (std::size_t) channels_ * (std::size_t) maxBlock_, 0.0f);
        dryBuf_.assign  ((std::size_t) channels_ * (std::size_t) maxBlock_, 0.0f);
        // A band processor whose own prepare() can FAIL (dynamics::Compressor refuses a configuration
        // it cannot honour rather than half-honouring it) must be able to say so through here, or the
        // caller learns about it as "the multiband does nothing" much later. Every band is still
        // prepared — returning early would leave the rest of them uninitialised — and the verdict is
        // ANDed. `if constexpr` keeps void-returning band processors (stereo::StereoWidth) working.
        bool bandsOk = true;
        for (int b = 0; b < MaxBands; ++b)
        {
            bandPtrs_[(std::size_t) b].assign ((std::size_t) channels_, nullptr);
            if constexpr (std::is_void_v<decltype (prepareBand (proc_[0]))>) prepareBand (proc_[(std::size_t) b]);
            else bandsOk = prepareBand (proc_[(std::size_t) b]) && bandsOk;
        }

        align_.assign ((std::size_t) MaxBands * (std::size_t) channels_, core::DelayLine {});
        for (auto& d : align_) d.prepare (maxAlignSamples);
        dryDelay_.assign ((std::size_t) channels_, core::DelayLine {});
        for (auto& d : dryDelay_) d.prepare (maxAlignSamples);

        refreshLatency(); reset();
        prepared_ = bandsOk;                                 // a band that refused leaves the composite unusable
        return bandsOk;
    }

    bool isPrepared() const noexcept { return prepared_; }

    void reset() noexcept
    {
        splitter_.reset();
        for (auto& p : proc_) p.reset();
        for (auto& d : align_) d.reset();
        ranNc_ = 0;   // nothing has run, so nothing can be stopping
        for (auto& d : dryDelay_) d.reset();
    }

    bool setNumBands (int n) noexcept { if (! splitter_.setNumBands (n)) return false; refreshLatency(); reset(); return true; }
    void setCrossovers (const float* hz, int count) noexcept { splitter_.setCrossovers (hz, count); }
    int  numBands() const noexcept { return splitter_.numBands(); }

    BandProcessor&       band (int b)       noexcept { return proc_[(std::size_t) std::clamp (b, 0, MaxBands - 1)]; }
    const BandProcessor& band (int b) const noexcept { return proc_[(std::size_t) std::clamp (b, 0, MaxBands - 1)]; }

    // Set a band's params + re-derive latency (a lookahead change shifts alignment). OOB index → ignored.
    template <class Params>
    void setBandParams (int b, const Params& p) noexcept { if ((unsigned) b < (unsigned) MaxBands) { proc_[(std::size_t) b].setParams (p); refreshLatency(); } }

    void setBandBypass (int b, bool x) noexcept { if ((unsigned) b < (unsigned) MaxBands) { bypass_[(std::size_t) b] = x; refreshLatency(); } }
    void setBandSolo   (int b, bool x) noexcept { if ((unsigned) b < (unsigned) MaxBands) solo_[(std::size_t) b] = x; }
    void setMix        (float m) noexcept       { mix_ = std::clamp (m, 0.0f, 1.0f); }
    void refreshLatency() noexcept { recomputeLatency(); }

    int  latencySamples() const noexcept { return latency_; }

    // Law 11. Geometry is checked BEFORE anything moves, so a refused call is indistinguishable from one
    // never made; then the LENGTH is chunked (it used to drop the whole call at n > maxBlock — measured,
    // max |out - in| over a 1024-sample call at maxBlock 256 was exactly 0: not one sample was touched).
    // Each chunk is one full pass of the four steps below, so the result is bit-identical to the caller
    // having made the same maxBlock-sized calls itself — every piece of state (the crossover, the band
    // processors, the alignment and dry delay lines) lives in members and carries across the seam.
    [[nodiscard]] bool process (float* const* io, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_) return false;
        if (numChannels > channels_) return false;          // width is a LIMIT — law 11(b)
        if (n == 0) return true;                            // no samples: no time, no edge, nothing
        bool ok = true;
        for (int off = 0; off < n; )
        {
            const int m = std::min (n - off, maxBlock_);
            float* sub[core::kMaxChannels] {};
            for (int c = 0; c < numChannels; ++c) sub[(std::size_t) c] = io[c] + off;
            ok = runChunk (sub, numChannels, m) && ok;
            off += m;                                       // `off += maxBlock_` could step past INT_MAX
        }
        return ok;
    }

private:
    // ONE chunk: nc is already validated and m <= maxBlock_.
    bool runChunk (float* const* io, int numChannels, int n) noexcept
    {
        const int nc = numChannels;
        const int nb = splitter_.numBands();

        // A channel that stops being split keeps the whole crossover tree for its column — every Svf in
        // every crossover and every allpass compensator — plus its per-band alignment delay lines. None of
        // it decays while the channel is away, and it replays on return: measured 1.55e-01 (-16.2 dBFS)
        // out of DIGITAL SILENCE. Per channel, never wholesale: the channels that stayed owe it nothing.
        // The band processors are not touched here — each carries its own ledger, or does not need one.
        // A BYPASSED band is stopped whatever the width does — it has been receiving nothing all along.
        // Telling it only at nc == 0 left the narrowing case leaking: stereo -> bypass -> mono -> bands
        // back on -> stereo silence emitted -16.6378 dBFS on the right, last non-zero at sample 239,
        // the whole lookahead line.
        for (int b = 0; b < nb; ++b)
            if (bypass_[(std::size_t) b])
                for (int c = nc; c < ranNc_; ++c)
                {
                    if constexpr (std::is_void_v<decltype (proc_[0].process (bandPtrs_[0].data(), nc, n))>)
                        proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), nc, n);
                    else
                        (void) proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), nc, n);
                    break;                          // one call per band carries the whole falling edge
                }
        for (int c = nc; c < ranNc_; ++c)
        {
            splitter_.resetChannel (c);
            for (int b = 0; b < MaxBands; ++b) align_[(std::size_t) (b * channels_ + c)].reset();
            // ...AND the parallel dry line, which is per-channel memory exactly like the rest. It was
            // missed because the obvious fixture cannot see it: at mix == 1 the sum is
            // `d + 1.0f * (wet - d)`, which is `wet` EXACTLY in IEEE, so the frozen dry cancels itself.
            // At mix 0.5 the whole 240-sample line replays: 0.25 out of DIGITAL SILENCE (-12.0 dBFS),
            // last non-zero at sample 239. My own zero-width probe ran at mix 1 and read 0.
            dryDelay_[(std::size_t) c].reset();
        }
        ranNc_ = nc;
        if (nc == 0)
        {
            // LAW 11(a)+(d): the gap has to reach the BANDS. Returning here left every band processor
            // holding its own frozen per-channel state — measured through MultibandCompressor with 5 ms
            // of lookahead: 0.308 out of DIGITAL SILENCE (-10.2 dBFS), last non-zero at sample 239, the
            // whole delay line. A composite that swallows the gap reproduces, one storey up, exactly the
            // defect its own bands were fixed for.
            // A BYPASSED band is stopped too — it has been receiving nothing all along, and the gap is
            // the moment its own per-channel state has to go, or it replays on un-bypass. Measured with
            // the bands bypassed during the gap: -11.83 dBFS out of digital silence.
            bool zeroOk = true;
            for (int b = 0; b < nb; ++b)
            {
                if constexpr (std::is_void_v<decltype (proc_[0].process (bandPtrs_[0].data(), nc, n))>)
                    proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), 0, n);
                else
                    zeroOk = proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), 0, n) && zeroOk;
            }
            return zeroOk;
        }

        // 1) split into per-band planar buffers; capture the allpass-reconstructed dry for the parallel mix
        float tmp[(std::size_t) MaxBands] {};
        for (int c = 0; c < nc; ++c)
        {
            float* dry = dryData (c);
            for (int i = 0; i < n; ++i)
            {
                splitter_.splitSample (c, io[c][i], tmp);
                dry[i] = splitter_.sumSample (tmp);
                for (int b = 0; b < nb; ++b) bandData (b, c)[i] = tmp[b];
            }
        }
        // LAW 8. The crossover's SVFs are feedback state and nothing was flushing them: the per-band
        // Compressor flushes itself (Compressor.h), and MultibandSplitter has offered flushDenormals()
        // (MultibandSplitter.h) since it was written, but no caller in this module ever invoked it — so on
        // silence the crossover state decayed into the subnormals and stayed there, at 10-100x the cost on
        // any CPU without hardware FTZ. Once per block, after the split loop, exactly like every other
        // kernel in the core. Found while fixing the same law's violation in analysis::KWeightingFilter.
        splitter_.flushDenormals();

        // 2) process each band in place (bypassed bands keep their split signal)
        // A band's own verdict is ANDed in. It cannot be checked before the split (the bands run on the
        // SPLIT signal, which does not exist yet), so a band that refuses leaves this chunk partially
        // processed — the composite's "refused ⇒ nothing moved" guarantee covers the geometry THIS class
        // validates, which is the part a caller can get wrong. A band refusing at a width and length this
        // class already validated means someone re-prepared that band alone, through band(i).
        // `if constexpr` keeps void-returning band processors (stereo::StereoWidth) working, exactly as
        // prepare() already does for prepareBand.
        bool bandsOk = true;
        for (int b = 0; b < nb; ++b)
        {
            if (bypass_[(std::size_t) b]) continue;
            for (int c = 0; c < nc; ++c) bandPtrs_[(std::size_t) b][(std::size_t) c] = bandData (b, c);
            if constexpr (std::is_void_v<decltype (proc_[0].process (bandPtrs_[0].data(), nc, n))>)
                proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), nc, n);
            else
                bandsOk = proc_[(std::size_t) b].process (bandPtrs_[(std::size_t) b].data(), nc, n) && bandsOk;
        }

        // 3) align every band to the max band latency
        for (int b = 0; b < nb; ++b)
            for (int c = 0; c < nc; ++c)
            {
                core::DelayLine& d = align_[(std::size_t) (b * channels_ + c)];
                float* x = bandData (b, c);
                for (int i = 0; i < n; ++i) x[i] = d.process (x[i]);
            }

        // 4) sum (honouring solo) + global dry/wet against the latency-matched allpass dry
        bool anySolo = false; for (int b = 0; b < nb; ++b) anySolo |= solo_[(std::size_t) b];
        for (int c = 0; c < nc; ++c)
        {
            float* dry = dryData (c);
            for (int i = 0; i < n; ++i)
            {
                for (int b = 0; b < nb; ++b) tmp[b] = (! anySolo || solo_[(std::size_t) b]) ? bandData (b, c)[i] : 0.0f;
                const float wet = splitter_.sumSample (tmp);
                const float d   = dryDelay_[(std::size_t) c].process (dry[i]);   // always advance the delay
                // Solo is a monitor: hear ONLY the soloed bands (full), bypassing the parallel dry — else the
                // full-band dry would leak the non-soloed bands back in at mix < 1.
                io[c][i] = anySolo ? wet : (d + mix_ * (wet - d));
            }
        }
        return bandsOk;
    }

    void recomputeLatency() noexcept
    {
        const int nb = splitter_.numBands();
        latency_ = 0;
        for (int b = 0; b < nb; ++b) latency_ = std::max (latency_, bypass_[(std::size_t) b] ? 0 : proc_[(std::size_t) b].latencySamples());
        for (int b = 0; b < nb; ++b)
        {
            const int bl = bypass_[(std::size_t) b] ? 0 : proc_[(std::size_t) b].latencySamples();
            const int extra = std::max (0, latency_ - bl);
            for (int c = 0; c < channels_; ++c) align_[(std::size_t) (b * channels_ + c)].setDelay (extra);
        }
        for (auto& d : dryDelay_) d.setDelay (latency_);
    }

    float* bandData (int b, int c) noexcept { return bandBuf_.data() + ((std::size_t) b * (std::size_t) channels_ + (std::size_t) c) * (std::size_t) maxBlock_; }
    float* dryData  (int c)        noexcept { return dryBuf_.data()  + (std::size_t) c * (std::size_t) maxBlock_; }

    int    ranNc_ = 0;                    // channels that advanced state on the previous call
    double fs_ = 48000.0;
    int maxBlock_ = 0, channels_ = 0, latency_ = 0;
    bool prepared_ = false;               // true only after a prepare() in which every band succeeded
    float mix_ = 1.0f;

    eq::MultibandSplitter<MaxBands> splitter_;
    std::array<BandProcessor, (std::size_t) MaxBands> proc_ {};
    std::array<bool, (std::size_t) MaxBands> bypass_ {}, solo_ {};
    std::array<std::vector<float*>, (std::size_t) MaxBands> bandPtrs_;
    std::vector<float> bandBuf_, dryBuf_;
    std::vector<core::DelayLine> align_, dryDelay_;
};

} // namespace felitronics::multiband
