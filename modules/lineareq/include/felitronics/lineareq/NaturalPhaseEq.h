// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqEngine.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/lineareq/MixedPhaseFir.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <memory>
#include <atomic>
#include <cmath>
#include <complex>
#include <vector>

namespace felitronics::lineareq
{

//==============================================================================
// felitronics::lineareq::NaturalPhaseEq — the "Natural phase" rendering of an eq::EqEngine bank: the
// mastering middle ground between Zero-Latency (minimum-phase IIR) and Linear (zero-phase FIR).
//
// It builds MIXED-PHASE FIRs (MixedPhaseFir, cepstral blend φ = k·φ_min) and convolves with the click-free
// convolution::MatrixConvolver, picking a topology per snapshot (basis detection):
//   • MSDiag / LRDiag (diagonal): blend the two AXIS composites directly (Mid+Side, or Left+Right) — the
//     min-phase blend distributes over the products in a diagonal entry, so this is today's exact path and
//     cost (2 cepstra), no matrix machinery.
//   • Full (mixed L/R AND M/S lanes): the blend does NOT distribute over the SUMS inside an H_MS entry
//     (0.5(H_M±H_S)) — min phase is multiplicative, not additive. So it is applied PER SCALAR LANE response
//     first (MixedPhaseFir::buildSpectrum), THEN the complex 2×2 entries are composed on the D-grid
//     (eq::accumulateBand), and each entry is inverse-transformed to its own causal FIR. Costs up to
//     one cepstrum PER ACTIVE (band,lane) — only engaged when domains are genuinely mixed.
//
// At k≈0.5 the phase stays close to flat (most of linear's benefit) while pre-ringing drops sharply — and
// the reported PDC latency is a FIXED L/4 (¼ the FIR, the same for every k AND every topology), FAR below
// linear's L/2. Every entry FIR shares that one bulk delay.
//
//  • Design FFT D = 8·L masks cepstral time-aliasing on steep filters.
//  • Each impulse is shifted causal by a FIXED `bulkDelay = L/4` (that shift IS the latency) and truncated
//    to L taps with a TAIL taper (NO centred window — that would gut the front-loaded impulse).
//  • A flat EQ renders a unit impulse at `bulkDelay` ⇒ exact unity-gain pass-through at the reported latency.
//
// Same threading contract as LinearPhaseEq: setBands()/buildFir()/prepare() are RT-UNSAFE (message thread,
// host-serialised); process() is RT-safe. To change quality, re-prepare(); the blend k changes LIVE via
// setBlend() (the bulk delay is fixed, so PDC never moves — no re-prepare).
class NaturalPhaseEq
{
public:
    static constexpr int kNumQuality = 4;                       // kept-IR length (L) per quality — lighter than Linear
    static int firSizeForQuality (int q) noexcept
    {
        static constexpr int sizes[kNumQuality] = { 2048, 4096, 8192, 16384 };
        return sizes[std::clamp (q, 0, kNumQuality - 1)];
    }

    // quality picks L (the convolution IR length). k ∈ [0,1] blends phase (0 linear … 1 minimum). RT-UNSAFE.
    bool prepare (double sampleRate, int maxBlock, int numChannels, int quality, float k = 0.5f) noexcept
    {
        prepared_ = false;                                      // any early return below leaves it unprepared
        fs_       = sampleRate > 0.0 ? sampleRate : 48000.0;
        maxBlock_ = std::max (1, maxBlock);
        channels_ = std::clamp (numChannels, 1, core::kMaxChannels);
        k_.store (std::clamp (k, 0.0f, 1.0f), std::memory_order_relaxed);
        L_        = firSizeForQuality (quality);
        D_        = 8 * L_;                                      // cepstral design size (≥ 8× the kept IR)
        bulkDelay_ = L_ / 4;                                     // FIXED causal shift → a STABLE reported latency

        if (! mp_.prepare (D_)) return false;
        if (! dFft_.prepare (D_)) return false;                 // size-D IFFT for the composed matrix entries (Full)
        magBuf_.assign ((std::size_t) (D_ / 2 + 1), 0.0f);
        for (auto& f : fir_) f.assign ((std::size_t) L_, 0.0f); // up to 4 entry IRs
        for (auto& g : laneRe_) g.assign ((std::size_t) (D_ / 2 + 1), 0.0);   // per-lane blended spectra (Full)
        for (auto& g : laneIm_) g.assign ((std::size_t) (D_ / 2 + 1), 0.0);
        for (auto& g : accRe_)  g.assign ((std::size_t) (D_ / 2 + 1), 0.0);   // composed 2×2 entry spectra (Full)
        for (auto& g : accIm_)  g.assign ((std::size_t) (D_ / 2 + 1), 0.0);
        packSpec_.assign  ((std::size_t) D_, 0.0f);
        entryTime_.assign ((std::size_t) D_, 0.0f);
        computeTaper();

        int part = 64; while (part < maxBlock_) part <<= 1;     // partition ≥ maxBlock, pow2
        const int warmXfade = std::max (2 * part, (int) std::lround (0.02 * fs_));   // short anti-click fade
        chConv_.clear();
        if (channels_ > 2)                                      // non-stereo: per-channel mono ST-only operators
        {
            for (int c = 0; c < channels_; ++c)
            {
                chConv_.push_back (std::make_unique<Conv>());
                if (! chConv_.back()->prepare (part, L_, warmXfade, 1)) return false;
            }
        }
        else if (! conv_.prepare (part, L_, warmXfade, channels_ == 1 ? 1 : 2)) return false;
        prepared_ = true;                                       // fully built — setBands()/buildFir() may now run
        return true;
    }

    void reset() noexcept { conv_.reset(); for (auto& c : chConv_) c->reset(); }

    int  firSize() const noexcept { return L_; }
    int  latencySamples() const noexcept { return bulkDelay_; }   // the causal bulk shift (convolver adds 0)
    bool isBusy() const noexcept
    {
        if (conv_.isBusy()) return true;
        for (auto& c : chConv_) if (c->isBusy()) return true;
        return false;
    }
    float blend() const noexcept { return k_.load (std::memory_order_relaxed); }

    // Change the phase blend LIVE — no re-prepare (the bulk delay / reported latency is fixed). The next
    // setBands() rebuilds the FIRs with the new k via the host's off-thread builder (a click-free swap).
    void setBlend (float k) noexcept { k_.store (std::clamp (k, 0.0f, 1.0f), std::memory_order_relaxed); }

    // RT-UNSAFE (message thread): rebuild the mixed-phase entry FIRs from a snapshot, pick the topology
    // (basis detection), hand the operator to the convolver for a click-free swap. False if a swap is still
    // crossfading (host coalesces with the latest snapshot).
    bool setBands (const eq::BandParams* bands, int numBands) noexcept
    {
        if (! prepared_) return false;
        if (isBusy()) return false;                            // coalesce BEFORE any expensive rebuild

        // NON-STEREO bus (channels_ > 2): the v2 engine runs ONLY the ST lane, per channel — so does
        // the FIR: one ST-only-composite IR convolved independently on every channel.
        if (channels_ > 2)
        {
            buildFir (bands, numBands, eq::Axis::Stereo, fir_[0].data());
            const float* one[1] { fir_[0].data() };
            bool ok = true;                                     // stage ALL channels first (expensive builds)...
            for (auto& c : chConv_) ok = c->stageOperator (Conv::Topology::LRDiag, one, 1, L_) && ok;
            if (! ok) return false;
            for (auto& c : chConv_) c->publishStaged();         // ...then publish back-to-back: lockstep fades
            return true;
        }
                         // unprepared — FIR buffers + MixedPhaseFir are empty

        // MONO (channels_ == 1): the v2 engine runs ONLY the ST lane on a non-stereo bus → bank 0 is the
        // ST-only composite (Axis::Stereo), and the FIR matches the IIR (LANES.md §FIR, mono).
        if (channels_ == 1)
        {
            buildFir (bands, numBands, eq::Axis::Stereo, fir_[0].data());
            return conv_.setIr (fir_[0].data(), L_);
        }

        using Topo = Conv::Topology;
        const eq::MatrixBasis basis = eq::detectBasis (bands, numBands);
        Topo topo = Topo::MSDiag;
        int  numBanks = 2;
        switch (basis)
        {
            case eq::MatrixBasis::MSDiag:                 // diagonal in M/S — blend axis composites (today's path)
                buildFir (bands, numBands, eq::Axis::Mid,  fir_[0].data());
                buildFir (bands, numBands, eq::Axis::Side, fir_[1].data());
                topo = Topo::MSDiag; numBanks = 2;
                break;
            case eq::MatrixBasis::LRDiag:                 // diag(H_LL,H_RR) — blend Left/Right composites
                buildFir (bands, numBands, eq::Axis::Left,  fir_[0].data());
                buildFir (bands, numBands, eq::Axis::Right, fir_[1].data());
                topo = Topo::LRDiag; numBanks = 2;
                break;
            case eq::MatrixBasis::Full:                   // mixed — blend per lane, compose complex entries
                buildFullBanks (bands, numBands);
                topo = Topo::Full; numBanks = 4;
                break;
        }
        const float* banks[4] { fir_[0].data(), fir_[1].data(), fir_[2].data(), fir_[3].data() };
        return conv_.setOperator (topo, banks, numBanks, L_);
    }

    // RT-UNSAFE: build ONE axis's mixed-phase FIR (from its composite MAGNITUDE) into out[0..L). Exposed for
    // tests / a host that drives the convolver itself. magnitude → MixedPhaseFir (phase k·φ_min) → causal
    // shift by bulkDelay + taper.
    void buildFir (const eq::BandParams* bands, int numBands, eq::Axis axis, float* out) noexcept
    {
        if (! prepared_) return;                               // unprepared — magBuf_/taper_ empty + MixedPhaseFir not built
        eq::EqEngine::magnitudeGridFor (bands, numBands, fs_, magBuf_.data(), D_ / 2 + 1, axis);
        const float* h = mp_.build (magBuf_.data(), k_.load (std::memory_order_relaxed));   // D-point mixed-phase impulse
        shiftTaper (h, out);
    }

    // RT-safe (audio thread), in place. Stereo (≥2 ch) → the picked topology; mono → bank 0. `n` ≤ maxBlock.
    // CONTRACT: the call's channel topology must match prepare() — bank 0 is built for the PREPARED count;
    // re-prepare on a bus change, as the adapter does.
    void process (float* const* io, int numChannels, int n) noexcept
    {
        if (n <= 0) return;
        if (channels_ > 2)                                     // non-stereo: ST-only IR per channel
        {
            const int nc = numChannels < channels_ ? numChannels : channels_;
            for (int c = 0; c < nc; ++c) { float* one[1] { io[c] }; chConv_[(std::size_t) c]->process (one, one, 1, n); }
            return;
        }
        const int nc = numChannels >= 2 ? 2 : 1;
        conv_.process (io, io, nc, n);                         // raw L/R in → MatrixConvolver routes + decodes in place
    }

private:
    using Conv = convolution::MatrixConvolver<core::fft::DefaultRealFft>;

    // Extract L causal taps from a D-point mixed-phase impulse: peak at h[0], pre-ring wraps to h[D-1..];
    // shift right by bulkDelay so out[bulkDelay] = h[0]. Tail (+ light head) taper suppresses truncation ripple.
    void shiftTaper (const float* h, float* out) noexcept
    {
        for (int i = 0; i < L_; ++i)
        {
            int idx = i - bulkDelay_;
            if (idx < 0) idx += D_;
            out[(std::size_t) i] = h[(std::size_t) idx] * taper_[(std::size_t) i];
        }
    }

    // Full topology: blend each ACTIVE scalar lane response per band (MixedPhaseFir::buildSpectrum), compose
    // the complex 2×2 on the D-grid (eq::accumulateBand — the SAME composition the analytic matrix uses),
    // then inverse-transform each entry to its own causal FIR. Idle lanes contribute exact identity (skipped).
    void buildFullBanks (const eq::BandParams* bands, int numBands) noexcept
    {
        const int H = D_ / 2;
        const float k = k_.load (std::memory_order_relaxed);

        for (int bk = 0; bk <= H; ++bk)                        // accumulator = identity 2×2
        {
            accRe_[0][(std::size_t) bk] = 1.0; accIm_[0][(std::size_t) bk] = 0.0;   // hLL
            accRe_[1][(std::size_t) bk] = 0.0; accIm_[1][(std::size_t) bk] = 0.0;   // hLR
            accRe_[2][(std::size_t) bk] = 0.0; accIm_[2][(std::size_t) bk] = 0.0;   // hRL
            accRe_[3][(std::size_t) bk] = 1.0; accIm_[3][(std::size_t) bk] = 0.0;   // hRR
        }

        for (int i = 0; i < numBands; ++i)
        {
            const eq::BandParams& b = bands[i];
            bool active[eq::kNumLanes];
            for (int li = 0; li < eq::kNumLanes; ++li)
            {
                const eq::Lane l = static_cast<eq::Lane> (li);
                active[li] = eq::laneActive (b, l);
                if (active[li])
                {
                    fillLaneMag (b, l, magBuf_.data());        // |H_lane| on the D-grid (analytic)
                    mp_.buildSpectrum (magBuf_.data(), k, laneRe_[(std::size_t) li].data(), laneIm_[(std::size_t) li].data());
                }
            }
            for (int bk = 0; bk <= H; ++bk)                    // fold this band into acc, per bin
            {
                auto lane = [&] (int li) noexcept -> std::complex<double>
                {
                    return active[li] ? std::complex<double> (laneRe_[(std::size_t) li][(std::size_t) bk],
                                                              laneIm_[(std::size_t) li][(std::size_t) bk])
                                      : std::complex<double> (1.0, 0.0);
                };
                eq::ResponseMatrix acc {
                    { accRe_[0][(std::size_t) bk], accIm_[0][(std::size_t) bk] },
                    { accRe_[1][(std::size_t) bk], accIm_[1][(std::size_t) bk] },
                    { accRe_[2][(std::size_t) bk], accIm_[2][(std::size_t) bk] },
                    { accRe_[3][(std::size_t) bk], accIm_[3][(std::size_t) bk] } };
                eq::accumulateBand (acc, lane (0), lane (1), lane (2), lane (3), lane (4));   // ST, L, R, M, S
                accRe_[0][(std::size_t) bk] = acc.hLL.real(); accIm_[0][(std::size_t) bk] = acc.hLL.imag();
                accRe_[1][(std::size_t) bk] = acc.hLR.real(); accIm_[1][(std::size_t) bk] = acc.hLR.imag();
                accRe_[2][(std::size_t) bk] = acc.hRL.real(); accIm_[2][(std::size_t) bk] = acc.hRL.imag();
                accRe_[3][(std::size_t) bk] = acc.hRR.real(); accIm_[3][(std::size_t) bk] = acc.hRR.imag();
            }
        }

        for (int e = 0; e < 4; ++e)                            // each entry: pack Hermitian spectrum → IFFT → shift+taper
        {
            packSpec_[0] = (float) accRe_[(std::size_t) e][0];               // DC (real)
            packSpec_[1] = (float) accRe_[(std::size_t) e][(std::size_t) H]; // Nyquist (real)
            for (int bk = 1; bk < H; ++bk)
            {
                packSpec_[(std::size_t) (2 * bk)]     = (float) accRe_[(std::size_t) e][(std::size_t) bk];
                packSpec_[(std::size_t) (2 * bk + 1)] = (float) accIm_[(std::size_t) e][(std::size_t) bk];
            }
            dFft_.inverse (packSpec_.data(), entryTime_.data());   // D-point mixed-phase impulse (circular)
            shiftTaper (entryTime_.data(), fir_[(std::size_t) e].data());
        }
    }

    // |H_lane| of ONE lane of ONE band on the D-grid (DC…Nyquist), analytic (no FFT) — same grid as
    // magnitudeGridFor: bin k ↔ w = 2πk/D.
    void fillLaneMag (const eq::BandParams& b, eq::Lane l, float* out) noexcept
    {
        const int H = D_ / 2;
        const eq::BandParams v = eq::laneView (b, l);
        for (int k = 0; k <= H; ++k)
        {
            const double w = 2.0 * core::kPi * (double) k / (double) D_;
            out[(std::size_t) k] = (float) std::abs (eq::bandResponse (v, fs_, w));
        }
    }

    // Tukey-ish taper: a short rising head (suppress any clipped pre-ring) + a longer falling tail
    // (suppress post-ring truncation). NOT a centred window — the bulk of the impulse stays untouched.
    void computeTaper() noexcept
    {
        taper_.assign ((std::size_t) L_, 1.0f);
        const int head = std::min (bulkDelay_ / 2, L_ / 16);                 // ≤ half the pre-ring region
        const int tail = std::max (1, L_ / 8);                              // last ~12.5%
        for (int i = 0; i < head; ++i)
            taper_[(std::size_t) i] = (float) (0.5 * (1.0 - std::cos (core::kPi * (double) i / (double) head)));
        for (int i = 0; i < tail; ++i)
        {
            const int idx = L_ - 1 - i;
            taper_[(std::size_t) idx] = (float) (0.5 * (1.0 - std::cos (core::kPi * (double) i / (double) tail)));
        }
    }

    double fs_ = 48000.0;
    int maxBlock_ = 512, channels_ = 2;
    int L_ = 4096, D_ = 32768, bulkDelay_ = 1024;
    bool prepared_ = false;                                    // true only after a fully-successful prepare()
    std::atomic<float> k_ { 0.5f };   // phase blend — live-settable without a re-prepare (the builder reads it)

    MixedPhaseFir<core::fft::DefaultRealFft> mp_;                            // cepstral mixed-phase FIR designer
    core::fft::DefaultRealFft dFft_;                                         // size-D IFFT for composed matrix entries
    Conv conv_;
    std::vector<std::unique_ptr<Conv>> chConv_;                 // channels_ > 2: one mono ST-only operator per channel                                                              // matrix operator convolver, click-free swap

    std::vector<float> magBuf_, taper_, packSpec_, entryTime_;
    std::vector<float>  fir_[4];                                             // up to 4 entry IRs (diag topologies use [0..1])
    std::vector<double> laneRe_[eq::kNumLanes], laneIm_[eq::kNumLanes];      // per-lane blended spectra (Full)
    std::vector<double> accRe_[4], accIm_[4];                                // composed 2×2 entry spectra (Full)
};

} // namespace felitronics::lineareq
