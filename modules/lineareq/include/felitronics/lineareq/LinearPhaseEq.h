// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/eq/EqEngine.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <memory>
#include <cmath>
#include <vector>

namespace felitronics::lineareq
{

//==============================================================================
// felitronics::lineareq::LinearPhaseEq — a linear-phase rendering of an eq::EqEngine bank. It samples the
// bank's exact 2×2 transfer matrix on the FFT grid with every SCALAR LANE response replaced by its
// zero-phase MAGNITUDE (eq::matrixResponseZeroPhase), so all matrix entries are REAL (but SIGNED) → each
// entry inverse-transforms to a SYMMETRIC FIR: the EQ's amplitude shape, ZERO phase shift, constant group
// delay = N/2.
//
// Built to FIX the plugin prototype it's extracted from (which stuttered rebuilding a huge FIR on every
// node drag):
//   • CORE owns NO thread. setBands()/buildFir()/prepare() are RT-UNSAFE (FFT/IFFT/window) — the HOST calls
//     them off the audio thread and DEBOUNCES. process() is RT-safe. setBands() returns false while a swap
//     is mid-crossfade (isBusy()) → the host coalesces with the latest snapshot.
//   • Each entry FIR is N+1 symmetric taps (type-I linear phase): a real (signed) grid → core::fft inverse
//     → fftshift → 4-term Blackman-Harris (peaks EXACTLY 1.0 at the centre tap N/2). So a flat EQ renders an
//     exact unit impulse at N/2 — no run-time gain hack.
//   • Convolution is convolution::MatrixConvolver on ONE canonical raw-L/R history, picking a topology per
//     snapshot (basis detection): MSDiag (2 IRs, the classic M/S path), LRDiag (2 IRs, direct L/R) or Full
//     (4 IRs + cross sums). Every topology shares the same warm history, so even ENABLING the first L lane
//     mid-playback (an MSDiag→Full topology change) crossfades click-free.
//
// Latency = N/2 (the FIR group delay; the convolver adds 0) — 2048…65536 samples, so this is an OFFLINE /
// MASTERING tool, not for live monitoring. CONTRACTS: setBands()/buildFir()/prepare() are message-thread
// (host-serialised — a single producer); reset() is audio-thread (or externally synced) and cancels any
// pending swap. To change quality (N) re-prepare(). RT-safe: process() never allocates/locks/throws.
template <core::fft::RealFftBackend AudioFft = core::fft::DefaultRealFft>
class BasicLinearPhaseEq
{
public:
    static constexpr int kNumQuality = 5;                       // Low / Medium / High / Very High / Maximum
    static int firSizeForQuality (int q) noexcept
    {
        // Power-of-two FIR lengths (our FFT is radix-2). Latency = N/2 ≈ FabFilter's linear-phase ladder
        // — 2048 / 4096 / 8192 / 16384 / 65536 vs FabFilter's 3072 / 5120 / 9216 / 17408 / 66560 @ 44.1k:
        // comparable, converging to equal at Maximum.
        static constexpr int sizes[kNumQuality] = { 4096, 8192, 16384, 32768, 131072 };
        return sizes[std::clamp (q, 0, kNumQuality - 1)];
    }

    // RT-UNSAFE (message thread): allocates the FIRs/FFT/convolver. quality ∈ [0,4] picks the FIR length.
    bool prepare (double sampleRate, int maxBlock, int numChannels, int quality) noexcept
    {
        prepared_ = false;                                     // any early return below leaves it unprepared
        fs_ = sampleRate > 0.0 ? sampleRate : 48000.0;
        maxBlock_ = std::max (1, maxBlock);
        channels_ = std::clamp (numChannels, 1, core::kMaxChannels);
        N_ = firSizeForQuality (quality);

        if (! buildFft_.prepare (N_)) return false;
        magBuf_.assign ((std::size_t) (N_ / 2 + 1), 0.0f);
        spec_.assign   ((std::size_t) DesignFft::spectrumFloats (N_), 0.0f);   // packed-Hermitian scratch (== N_ for the pinned scalar)
        time_.assign   ((std::size_t) N_, 0.0f);
        for (auto& g : grid_) g.assign ((std::size_t) (N_ / 2 + 1), 0.0f);   // 4 signed-real entry grids (Full)
        for (auto& f : fir_)  f.assign ((std::size_t) (N_ + 1), 0.0f);       // up to 4 entry IRs
        computeWindow();

        int part = 64; while (part < maxBlock_) part <<= 1;     // partition ≥ maxBlock, pow2
        // SHORT crossfade for interactive swaps (MatrixConvolver keeps ONE warm L/R history shared by two
        // operator slots, so a swapped operator is response-correct immediately — only a short anti-click
        // fade is needed, ≥ the partition P). The engine derives its OWN long fade for the COLD first
        // activation, so an EQ band drag lands in ~tens of ms, not a full FIR. isBusy() is then true only
        // for that short window → the host's coalescing rebuild stays responsive.
        const int warmXfade = std::max (2 * part, (int) std::lround (0.02 * fs_));   // ≥ 2P and ≥ ~20 ms
        chConv_.clear();
        if (channels_ > 2)                                      // non-stereo: per-channel mono ST-only operators
        {
            for (int c = 0; c < channels_; ++c)
            {
                chConv_.push_back (std::make_unique<Conv>());
                if (! chConv_.back()->prepare (part, N_ + 1, warmXfade, 1)) return false;
            }
        }
        else if (! conv_.prepare (part, N_ + 1, warmXfade, channels_ == 1 ? 1 : 2)) return false;
        prepared_ = true;                                      // fully built — setBands()/buildFir() may now run
        return true;
    }

    void reset() noexcept { conv_.reset(); for (auto& c : chConv_) c->reset(); }

    int  firSize() const noexcept { return N_; }
    int  latencySamples() const noexcept { return N_ / 2; }     // symmetric FIR group delay (convolver adds 0)
    bool isBusy() const noexcept
    {
        if (conv_.isBusy()) return true;
        for (auto& c : chConv_) if (c->isBusy()) return true;
        return false;
    }     // a swap is mid-crossfade → coalesce setBands()

    // RT-UNSAFE (message thread, off the audio thread): rebuild the entry FIRs from a caller-owned
    // BandParams snapshot, pick the convolver topology (basis detection), and hand the operator over for a
    // click-free swap. Returns false if a swap is still crossfading (isBusy()) — the host retries with the
    // LATEST snapshot once it's free.
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
            for (auto& c : chConv_) ok = c->stageOperator (Conv::Topology::LRDiag, one, 1, N_ + 1) && ok;
            if (! ok) return false;
            for (auto& c : chConv_) c->publishStaged();         // ...then publish back-to-back: lockstep fades
            return true;
        }
                            // unprepared — the FIR/scratch buffers are empty

        // MONO (channels_ == 1): the v2 engine runs ONLY the ST lane on a non-stereo bus, so bank 0 is the
        // ST-only composite (Axis::Stereo) and the FIR matches the IIR (LANES.md §FIR, mono).
        if (channels_ == 1)
        {
            buildFir (bands, numBands, eq::Axis::Stereo, fir_[0].data());
            return conv_.setIr (fir_[0].data(), N_ + 1);
        }

        using Topo = Conv::Topology;
        const eq::MatrixBasis basis = eq::detectBasis (bands, numBands);
        Topo topo = Topo::MSDiag;
        int  numBanks = 2;
        switch (basis)
        {
            case eq::MatrixBasis::MSDiag:                 // no active L/R lane ⇒ diagonal in M/S (today's 2-IR path)
                buildFir (bands, numBands, eq::Axis::Mid,  fir_[0].data());
                buildFir (bands, numBands, eq::Axis::Side, fir_[1].data());
                topo = Topo::MSDiag; numBanks = 2;
                break;
            case eq::MatrixBasis::LRDiag:                 // no active M/S lane ⇒ diag(H_LL,H_RR) — direct L/R
                buildFir (bands, numBands, eq::Axis::Left,  fir_[0].data());
                buildFir (bands, numBands, eq::Axis::Right, fir_[1].data());
                topo = Topo::LRDiag; numBanks = 2;
                break;
            case eq::MatrixBasis::Full:                   // mixed ⇒ the 4-entry matrix (signed real → symmetric)
                eq::EqEngine::matrixGridZeroPhase (bands, numBands, fs_, grid_[0].data(), grid_[1].data(),
                                                   grid_[2].data(), grid_[3].data(), N_ / 2 + 1);
                for (int k = 0; k < 4; ++k) gridToSymmetricFir (grid_[k].data(), fir_[k].data());
                topo = Topo::Full; numBanks = 4;
                break;
        }
        const float* banks[4] { fir_[0].data(), fir_[1].data(), fir_[2].data(), fir_[3].data() };
        return conv_.setOperator (topo, banks, numBanks, N_ + 1);
    }

    // RT-UNSAFE: build ONE axis's symmetric zero-phase FIR (from its composite MAGNITUDE) into out[0..N]
    // (N+1 taps). Exposed for tests / a host that drives the convolver itself.
    void buildFir (const eq::BandParams* bands, int numBands, eq::Axis axis, float* out) noexcept
    {
        if (! prepared_) return;                               // unprepared — magBuf_/spec_/time_/window_ are empty
        eq::EqEngine::magnitudeGridFor (bands, numBands, fs_, magBuf_.data(), N_ / 2 + 1, axis);
        gridToSymmetricFir (magBuf_.data(), out);
    }

    // RT-safe (audio thread), in place. Stereo (≥2 ch) → the picked topology; mono → bank 0. `n` ≤ maxBlock.
    // CONTRACT: the call's channel topology must match prepare() — bank 0 is built for the PREPARED count
    // (mono → ST-only composite, stereo → the matrix); re-prepare on a bus change, as the adapter does.
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
    using DesignFft = core::fft::ScalarRadix2Real;   // FIR design hand-packs the scalar packed-Hermitian layout → pinned
    static_assert (core::fft::PackedHermitianSpectrum<DesignFft>, "design FFT must be packed-Hermitian (gridToSymmetricFir)");
    using Conv      = convolution::MatrixConvolver<AudioFft>;   // audio path — swappable to a SIMD backend

    // A real (signed) grid of N/2+1 points → a symmetric zero-phase FIR of N+1 taps. Shared by the axis
    // path (magnitude grid, always ≥ 0) and the Full entry path (signed matrix entries — an off-diagonal
    // can go negative; a signed real spectrum still inverse-transforms to a SYMMETRIC IR, keeping the sign).
    void gridToSymmetricFir (const float* grid, float* out) noexcept
    {
        const int N = N_;
        spec_[0] = grid[0];                                     // DC (real)
        spec_[1] = grid[(std::size_t) (N / 2)];                 // Nyquist (real)
        for (int k = 1; k < N / 2; ++k)                         // zero phase: imag = 0
        {
            spec_[(std::size_t) (2 * k)]     = grid[(std::size_t) k];
            spec_[(std::size_t) (2 * k + 1)] = 0.0f;
        }
        buildFft_.inverse (spec_.data(), time_.data());         // time_[0..N-1], 1/N-normalised, even-symmetric about 0
        for (int i = 0; i <= N; ++i)                            // fftshift (peak → N/2) + window → N+1 symmetric taps
            out[(std::size_t) i] = time_[(std::size_t) ((i + N / 2) % N)] * window_[(std::size_t) i];
    }

    void computeWindow() noexcept                              // 4-term Blackman-Harris over N+1 taps (peak 1.0 at N/2)
    {
        window_.assign ((std::size_t) (N_ + 1), 0.0f);
        constexpr double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
        const double M1 = (double) N_;                          // (taps - 1) = N
        for (int i = 0; i <= N_; ++i)
        {
            const double x = (double) i / M1;
            window_[(std::size_t) i] = (float) (a0 - a1 * std::cos (2.0 * core::kPi * x)
                                                   + a2 * std::cos (4.0 * core::kPi * x)
                                                   - a3 * std::cos (6.0 * core::kPi * x));
        }
    }

    double fs_ = 48000.0;
    int maxBlock_ = 512, channels_ = 2, N_ = 16384;
    bool prepared_ = false;                                    // true only after a fully-successful prepare()

    DesignFft buildFft_;                                       // size-N IFFT for the FIR design (pinned scalar)
    Conv conv_;
    std::vector<std::unique_ptr<Conv>> chConv_;                 // channels_ > 2: one mono ST-only operator per channel                                                 // matrix operator convolver (MSDiag / LRDiag / Full), click-free swap

    std::vector<float> magBuf_, spec_, time_, window_;
    std::vector<float> grid_[4];                                // 4 signed-real entry grids (Full only)
    std::vector<float> fir_[4];                                 // up to 4 entry IRs (MSDiag/LRDiag use [0..1])
};

// Keep the bare `lineareq::LinearPhaseEq` spelling working (default = the scalar audio backend); a consumer
// wanting a SIMD audio path spells BasicLinearPhaseEq<SomeSimdFft> for any core::fft::RealFftBackend (e.g. the
// option-gated pffft adapter) — the FIR-design FFTs stay scalar regardless.
using LinearPhaseEq = BasicLinearPhaseEq<>;

} // namespace felitronics::lineareq
