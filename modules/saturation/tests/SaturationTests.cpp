// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for the saturation module: peak-normalised WaveShaper curves (odd vs even
// harmonics), and the oversampled Saturator (clean dry/wet, ~linear at zero drive, peak-safe under
// full-scale drive, no-alloc in process(), and the DC blocker neutralising the asymmetric curve's offset).

#include <felitronics_test.h>
#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/saturation/Saturator.h>

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

static std::atomic<long> g_allocs { 0 };
void* operator new      (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void* operator new[]    (std::size_t s) { g_allocs.fetch_add (1, std::memory_order_relaxed); return std::malloc (s ? s : 1); }
void  operator delete   (void* p) noexcept { std::free (p); }
void  operator delete[] (void* p) noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t) noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t) noexcept { std::free (p); }

using namespace felitronics;
using Shape = saturation::WaveShaper::Shape;

//==============================================================================
// FROZEN pre-hardening Saturator — a verbatim copy of the class BEFORE the poison-guard/chunking
// hardening (WaveShaper/PolyphaseOversampler/DelayLine underneath are the live, untouched headers).
// The NULL test below streams identical finite scenarios through this and the live class and
// bit-compares: same TU, same compiler, same flags — any output difference is BY the guards, and
// there must be none. Fidelity of this copy was verified against byte goldens captured from a build
// of the pre-change tree (six scenarios, byte-identical).
namespace prechange
{
    namespace core = felitronics::core;
    namespace oversampling = felitronics::oversampling;
    using felitronics::saturation::WaveShaper;

    class Saturator
    {
    public:
        struct Params
        {
            WaveShaper::Shape shape = WaveShaper::Shape::Tanh;
            float driveDb   = 3.0f;
            float bias      = 0.0f;
            float mix       = 1.0f;
            float outputDb  = 0.0f;
            float autoComp  = 0.5f;
            float dcBlockHz = 10.0f;
        };

        bool prepare (double sampleRate, int maxBlock, int maxChannels, int oversampleFactor = 4, int tapsPerPhase = 32)
        {
            prepared_ = false;
            if (sampleRate <= 0.0 || maxBlock < 1) return false;
            fs_       = sampleRate;
            maxBlock_ = maxBlock;
            channels_ = std::clamp (maxChannels, 1, core::kMaxChannels);
            os_       = (oversampleFactor >= 2) ? oversampleFactor : 1;
            if (os_ > 1 && ! ovs_.prepare (os_, channels_, tapsPerPhase)) return false;

            osBuf_.assign  ((std::size_t) channels_ * (std::size_t) (maxBlock_ * os_), 0.0f);
            wetBuf_.assign ((std::size_t) channels_ * (std::size_t) maxBlock_,         0.0f);
            osPtrs_.assign ((std::size_t) channels_, nullptr);
            wetPtrs_.assign((std::size_t) channels_, nullptr);
            dcX1_.assign   ((std::size_t) channels_, 0.0f);
            dcY1_.assign   ((std::size_t) channels_, 0.0f);
            dryDelay_.assign ((std::size_t) channels_, core::DelayLine {});
            const int lat = (os_ > 1) ? ovs_.latencySamples() : 0;
            for (auto& d : dryDelay_) { d.prepare (lat); d.setDelay (lat); }
            applyParams();
            prepared_ = true;
            return true;
        }

        void reset() noexcept
        {
            if (os_ > 1) ovs_.reset();
            std::fill (dcX1_.begin(), dcX1_.end(), 0.0f);
            std::fill (dcY1_.begin(), dcY1_.end(), 0.0f);
            for (auto& d : dryDelay_) d.reset();
        }

        int  latencySamples() const noexcept { return os_ > 1 ? ovs_.latencySamples() : 0; }

        void setParams (const Params& p) noexcept { params_ = p; applyParams(); }

        void process (float* const* io, int numChannels, int n) noexcept
        {
            const int nc = std::min (numChannels, channels_);
            if (! prepared_ || nc <= 0 || n <= 0 || n > maxBlock_) return;
            const int osN = n * os_;
            for (int c = 0; c < nc; ++c)
            {
                osPtrs_[(std::size_t) c]  = &osBuf_[(std::size_t) c * (std::size_t) (maxBlock_ * os_)];
                wetPtrs_[(std::size_t) c] = &wetBuf_[(std::size_t) c * (std::size_t) maxBlock_];
            }

            if (os_ > 1) ovs_.upsample (io, nc, n, osPtrs_.data());
            else for (int c = 0; c < nc; ++c) std::copy (io[c], io[c] + n, osPtrs_[(std::size_t) c]);

            for (int c = 0; c < nc; ++c)
            {
                float* b = osPtrs_[(std::size_t) c];
                if (dcEnabled_)
                {
                    float x1 = dcX1_[(std::size_t) c], y1 = dcY1_[(std::size_t) c];
                    for (int i = 0; i < osN; ++i)
                    {
                        const float w  = shaper_.processSample (b[i]);
                        const float dc = w - x1 + dcR_ * y1;
                        x1 = w; y1 = dc;
                        b[i] = dc;
                    }
                    dcX1_[(std::size_t) c] = x1; dcY1_[(std::size_t) c] = y1;
                }
                else
                {
                    for (int i = 0; i < osN; ++i) b[i] = shaper_.processSample (b[i]);
                }
            }

            if (os_ > 1) ovs_.downsample (osPtrs_.data(), nc, n, wetPtrs_.data());
            else for (int c = 0; c < nc; ++c) std::copy (osPtrs_[(std::size_t) c], osPtrs_[(std::size_t) c] + n, wetPtrs_[(std::size_t) c]);

            for (int c = 0; c < nc; ++c)
            {
                core::DelayLine& dl = dryDelay_[(std::size_t) c];
                for (int i = 0; i < n; ++i)
                {
                    const float dry = dl.process (io[c][i]);
                    const float wet = comp_ * wetPtrs_[(std::size_t) c][i];
                    io[c][i] = outGain_ * ((1.0f - mix_) * dry + mix_ * wet);
                }
            }
        }

    private:
        static float finite (float v, float fallback) noexcept { return std::isfinite (v) ? v : fallback; }

        void applyParams() noexcept
        {
            const float driveDb  = finite (params_.driveDb,   3.0f);
            const float bias     = finite (params_.bias,      0.0f);
            const float autoComp = finite (params_.autoComp,  0.5f);
            const float dcHz     = finite (params_.dcBlockHz, 10.0f);
            shaper_.setShape (params_.shape);
            shaper_.setBias  (bias);
            shaper_.setDrive ((float) (core::dbToGain (driveDb) - 1.0));
            comp_    = (float) std::pow ((double) std::max (1.0e-6f, shaper_.slopeAtZero()),
                                         (double) -std::clamp (autoComp, 0.0f, 1.0f));
            mix_     = std::clamp (finite (params_.mix, 1.0f), 0.0f, 1.0f);
            outGain_ = (float) core::dbToGain (finite (params_.outputDb, 0.0f));
            const double fsOs = fs_ * (double) os_;
            const double fc   = std::clamp ((double) dcHz, 0.0, 0.49 * fsOs);
            dcR_ = (fc <= 0.0) ? 0.0f : (float) std::exp (-2.0 * core::kPi * fc / fsOs);
            dcEnabled_ = (dcHz > 0.0f) && (params_.shape == WaveShaper::Shape::Asym);
        }

        Params      params_ {};
        WaveShaper  shaper_ {};
        oversampling::PolyphaseOversampler ovs_;

        double fs_ = 0.0;
        int    maxBlock_ = 0, channels_ = 0, os_ = 1;
        float  comp_ = 1.0f, dcR_ = 0.0f, mix_ = 1.0f, outGain_ = 1.0f;
        bool   dcEnabled_ = false;
        bool   prepared_  = false;

        std::vector<float>  osBuf_, wetBuf_;
        std::vector<float*> osPtrs_, wetPtrs_;
        std::vector<float>  dcX1_, dcY1_;
        std::vector<core::DelayLine> dryDelay_;
    };
} // namespace prechange

//==============================================================================
// Deterministic scenario suite for the hardening NULL: finite, in-range signals covering os=1/2/4,
// every curve family, partial mix, edge buffers (silence / DC / ±FS / 1e-20 tiny / ramps) and a finite
// mid-stream param sweep. The SAME code generated the pre-change byte goldens (scratch build of the
// unmodified tree) — keep it in sync with that generator if scenarios ever change.
namespace nulltest
{
    struct Rng
    {
        explicit Rng (uint32_t seed) : s (seed) {}
        float next() { s = s * 1664525u + 1013904223u; return (float) (int32_t) s * 4.6566128730773926e-10f; }
        uint32_t s;
    };

    template <class Sat>
    std::vector<float> runScenario (int id)
    {
        std::vector<float> out, ch0, ch1;
        auto push = [&] (float* const* io, int nc, int n)
        {
            for (int c = 0; c < nc; ++c) out.insert (out.end(), io[c], io[c] + n);
        };

        Sat s;
        typename Sat::Params p;

        switch (id)
        {
            case 1:   // os=1, 2ch, Tanh on random noise, partial mix + trim
            {
                s.prepare (48000.0, 512, 2, 1);
                p.shape = Shape::Tanh; p.driveDb = 18.0f; p.mix = 0.7f; p.outputDb = -1.0f; p.autoComp = 0.5f;
                s.setParams (p);
                Rng r (0xC0FFEE01u);
                ch0.resize (512); ch1.resize (512);
                float* io[2] { ch0.data(), ch1.data() };
                for (int b = 0; b < 6; ++b)
                {
                    for (int i = 0; i < 512; ++i) { ch0[i] = r.next(); ch1[i] = r.next(); }
                    s.process (io, 2, 512);
                    push (io, 2, 512);
                }
                break;
            }
            case 2:   // os=4, 2ch, Tanh on two sines, full wet + autoComp 1
            {
                s.prepare (48000.0, 256, 2, 4);
                p.shape = Shape::Tanh; p.driveDb = 12.0f; p.mix = 1.0f; p.autoComp = 1.0f;
                s.setParams (p);
                ch0.resize (256); ch1.resize (256);
                float* io[2] { ch0.data(), ch1.data() };
                for (int b = 0; b < 8; ++b)
                {
                    for (int i = 0; i < 256; ++i)
                    {
                        const int t = b * 256 + i;
                        ch0[i] = 0.9f * (float) std::sin (2.0 * 3.141592653589793 *  997.0 * t / 48000.0);
                        ch1[i] = 0.9f * (float) std::sin (2.0 * 3.141592653589793 * 1499.0 * t / 48000.0 + 0.5);
                    }
                    s.process (io, 2, 256);
                    push (io, 2, 256);
                }
                break;
            }
            case 3:   // os=4, 1ch, Asym (DC blocker active) on sine+noise, half mix
            {
                s.prepare (48000.0, 512, 1, 4);
                p.shape = Shape::Asym; p.driveDb = 9.0f; p.bias = 0.4f; p.mix = 0.5f; p.outputDb = 2.0f; p.autoComp = 0.25f;
                s.setParams (p);
                Rng r (0xBADF00D5u);
                ch0.resize (512);
                float* io[1] { ch0.data() };
                for (int b = 0; b < 8; ++b)
                {
                    for (int i = 0; i < 512; ++i)
                    {
                        const int t = b * 512 + i;
                        ch0[i] = 0.7f * (float) std::sin (2.0 * 3.141592653589793 * 187.5 * t / 48000.0) + 0.2f * r.next();
                    }
                    s.process (io, 1, 512);
                    push (io, 1, 512);
                }
                break;
            }
            case 4:   // os=4, 2ch, Asym on EDGE buffers: silence, DC, ±FS alternation, tiny 1e-20, ramps
            {
                s.prepare (48000.0, 512, 2, 4);
                p.shape = Shape::Asym; p.driveDb = 24.0f; p.bias = -0.3f; p.mix = 1.0f; p.autoComp = 0.5f;
                s.setParams (p);
                Rng r (0x5EED0004u);
                ch0.resize (512); ch1.resize (512);
                float* io[2] { ch0.data(), ch1.data() };
                for (int b = 0; b < 5; ++b)
                {
                    for (int i = 0; i < 512; ++i)
                    {
                        switch (b)
                        {
                            case 0: ch0[i] = 0.0f;                          ch1[i] = 0.0f;                          break;
                            case 1: ch0[i] = 0.9f;                          ch1[i] = -0.9f;                         break;
                            case 2: ch0[i] = (i & 1) ? -1.0f : 1.0f;        ch1[i] = (i & 1) ? 1.0f : -1.0f;        break;
                            case 3: ch0[i] = 1.0e-20f * r.next();           ch1[i] = 1.0e-20f * r.next();           break;
                            default: { const float t = -1.0f + 2.0f * (float) i / 511.0f; ch0[i] = t; ch1[i] = -t; } break;
                        }
                    }
                    s.process (io, 2, 512);
                    push (io, 2, 512);
                }
                break;
            }
            case 5:   // os=4, 2ch, finite param sweep between blocks (shape switch mid-stream), odd block size
            {
                s.prepare (44100.0, 333, 2, 4);
                Rng r (0xAB12CD34u);
                ch0.resize (333); ch1.resize (333);
                float* io[2] { ch0.data(), ch1.data() };
                const float mixes[6] { 1.0f, 0.8f, 0.6f, 0.4f, 0.2f, 0.0f };
                for (int b = 0; b < 6; ++b)
                {
                    typename Sat::Params q;
                    q.shape    = b < 3 ? Shape::Tanh : Shape::Asym;
                    q.driveDb  = 3.0f + 3.0f * (float) b;
                    q.bias     = 0.1f * (float) b - 0.2f;
                    q.mix      = mixes[b];
                    q.outputDb = (float) b - 3.0f;
                    q.autoComp = 0.2f * (float) b;
                    s.setParams (q);
                    for (int i = 0; i < 333; ++i) { ch0[i] = r.next(); ch1[i] = r.next(); }
                    s.process (io, 2, 333);
                    push (io, 2, 333);
                }
                break;
            }
            case 6:   // os=2, 1ch, Cubic on random, custom tapsPerPhase, 96 kHz
            {
                s.prepare (96000.0, 512, 1, 2, 16);
                p.shape = Shape::Cubic; p.driveDb = 15.0f; p.mix = 0.9f; p.outputDb = -3.0f; p.autoComp = 1.0f;
                s.setParams (p);
                Rng r (0x0F0F0F0Fu);
                ch0.resize (512);
                float* io[1] { ch0.data() };
                for (int b = 0; b < 4; ++b)
                {
                    for (int i = 0; i < 512; ++i) ch0[i] = r.next();
                    s.process (io, 1, 512);
                    push (io, 1, 512);
                }
                break;
            }
        }
        return out;
    }
} // namespace nulltest

int main()
{
    std::printf ("felitronics::saturation tests\n");
    const double pi = core::kPi;

    // --- WaveShaper transfer curves ---
    test::group ("WaveShaper peak-normalised curves");
    {
        saturation::WaveShaper ws; ws.setShape (Shape::Tanh); ws.setDrive (2.0f);
        test::approx (ws.processSample ( 1.0f),  1.0, 1e-5, "tanh: peak preserved at +1");
        test::approx (ws.processSample (-1.0f), -1.0, 1e-5, "tanh: peak preserved at -1");
        test::approx (ws.processSample ( 0.0f),  0.0, 1e-9, "tanh: through origin");
        test::ok (ws.processSample (0.5f) > 0.5f, "tanh: small-signal boost (norm slope > 1)");
        test::approx (ws.processSample (0.3f) + ws.processSample (-0.3f), 0.0, 1e-6, "tanh: odd (no even harmonics)");
        test::ok (ws.slopeAtZero() > 1.0f, "tanh: slopeAtZero > 1 for drive > 0");

        saturation::WaveShaper a; a.setShape (Shape::Asym); a.setDrive (2.0f); a.setBias (0.4f);
        test::approx (a.processSample (0.0f), 0.0, 1e-6, "asym: through origin");
        test::ok (std::fabs (a.processSample (0.3f) + a.processSample (-0.3f)) > 1e-3, "asym: NOT odd -> even harmonics");
        test::ok (std::fabs (a.processSample ( 1.0f)) <= 1.0001f, "asym: peak bounded (+)");
        test::ok (std::fabs (a.processSample (-1.0f)) <= 1.0001f, "asym: peak bounded (-)");

        saturation::WaveShaper c; c.setShape (Shape::Cubic); c.setDrive (1.5f);
        test::approx (c.processSample (1.0f), 1.0, 1e-5, "cubic: peak preserved at +1");
        test::ok (c.processSample (5.0f) <= 1.0001f, "cubic: hard-bounded beyond the knee");
    }

    // --- Saturator: clean dry/wet + ~linear at zero drive (os = 1) ---
    test::group ("Saturator dry/wet + linearity (os=1)");
    {
        saturation::Saturator s; test::ok (s.prepare (48000.0, 512, 2, 1), "prepare os=1");
        test::ok (s.latencySamples() == 0, "os=1 -> zero latency");

        saturation::Saturator::Params p; p.driveDb = 18.0f; p.mix = 0.0f; p.outputDb = 0.0f; s.setParams (p);
        float a[6] { 0.1f, -0.2f, 0.3f, -0.4f, 0.5f, -0.6f }, ref[6];
        for (int i = 0; i < 6; ++i) ref[i] = a[i];
        float b6[6] {}; float* io[2] { a, b6 };
        s.process (io, 2, 6);
        double md = 0; for (int i = 0; i < 6; ++i) md = std::max (md, (double) std::fabs (a[i] - ref[i]));
        test::ok (md < 1e-6, "mix=0 -> exact dry passthrough");

        saturation::Saturator s2; s2.prepare (48000.0, 512, 1, 1);
        saturation::Saturator::Params lp; lp.driveDb = 0.0f; lp.mix = 1.0f; lp.autoComp = 0.0f; s2.setParams (lp);
        std::vector<float> x (512), y (512);
        for (int i = 0; i < 512; ++i) { x[i] = 0.5f * (float) std::sin (2.0 * pi * 1000.0 * i / 48000.0); y[i] = x[i]; }
        float* io2[1] { y.data() };
        s2.process (io2, 1, 512);
        double dl = 0; for (int i = 64; i < 512; ++i) dl = std::max (dl, (double) std::fabs (y[i] - x[i]));
        test::ok (dl < 5e-3, "drive 0 -> ~linear passthrough (os=1)");
    }

    // --- Saturator oversampled (os = 4): latency, peak-safety, no-alloc ---
    test::group ("Saturator oversampled (os=4)");
    {
        saturation::Saturator s; s.prepare (48000.0, 512, 2, 4);              // tpp = 32 default
        test::ok (s.latencySamples() == 31, "os=4 round-trip latency == tpp-1 == 31");

        saturation::Saturator::Params p; p.shape = Shape::Tanh; p.driveDb = 18.0f; p.mix = 1.0f; p.autoComp = 0.0f; s.setParams (p);
        std::vector<float> L (512), R (512); float* io[2] { L.data(), R.data() };
        double peak = 0;
        for (int blk = 0; blk < 8; ++blk)
        {
            for (int i = 0; i < 512; ++i) { const int n = blk * 512 + i; const float v = (float) std::sin (2.0 * pi * 1000.0 * n / 48000.0); L[i] = v; R[i] = v; }
            s.process (io, 2, 512);
            if (blk >= 4) for (int i = 0; i < 512; ++i) peak = std::max (peak, (double) std::fabs (L[i]));
        }
        test::ok (peak < 1.10, "peak-safe: full-scale sine stays ~bounded (peak-normalised curve)");

        for (int i = 0; i < 512; ++i) { L[i] = 0.3f; R[i] = -0.3f; }
        const long before = g_allocs.load();
        s.process (io, 2, 512); s.process (io, 2, 512);
        test::okNoAlloc (g_allocs.load() == before, "process() did not allocate (os=4)");
    }

    // --- Asymmetric curve -> the DC blocker removes the offset (with a contrast where it's off) ---
    test::group ("Saturator Asym -> DC-blocked");
    {
        const double f = 187.5;   // exactly 4 periods / 1024-sample block @48k -> a block mean = DC only
        saturation::Saturator::Params p; p.shape = Shape::Asym; p.driveDb = 12.0f; p.bias = 0.4f; p.mix = 1.0f; p.autoComp = 0.0f;

        saturation::Saturator on; on.prepare (48000.0, 1024, 1, 4); on.setParams (p);
        saturation::Saturator off; off.prepare (48000.0, 1024, 1, 4);
        saturation::Saturator::Params p0 = p; p0.dcBlockHz = 0.0f; off.setParams (p0);

        std::vector<float> xo (1024), xf (1024); double meanOn = 0, meanOff = 0;
        for (int blk = 0; blk < 8; ++blk)
        {
            for (int i = 0; i < 1024; ++i) { const int n = blk * 1024 + i; const float v = 0.7f * (float) std::sin (2.0 * pi * f * n / 48000.0); xo[i] = v; xf[i] = v; }
            float* ioOn[1] { xo.data() }; float* ioOff[1] { xf.data() };
            on.process (ioOn, 1, 1024); off.process (ioOff, 1, 1024);
            if (blk == 7)
            {
                double sOn = 0, sOff = 0; for (int i = 0; i < 1024; ++i) { sOn += xo[i]; sOff += xf[i]; }
                meanOn = sOn / 1024.0; meanOff = sOff / 1024.0;
            }
        }
        test::ok (std::fabs (meanOn)  < 1e-3, "DC blocker ON  -> output block-mean ~ 0");
        test::ok (std::fabs (meanOff) > 1e-3, "DC blocker OFF -> asym curve leaves a measurable DC offset");
    }

    // --- lifecycle/misuse: process() after a FAILED prepare (partial init) or with n>maxBlock must not OOB ---
    // prepare() sets channels_ (clamped >=1) BEFORE ovs_.prepare() can fail — a failed prepare left channels_>=1
    // with empty osBuf_/osPtrs_, and process() (nc>=1) indexed them. Guarded now; an oversized block is CHUNKED
    // to maxBlock internally (fully processed, never overruns osBuf_) — see the chunking-equivalence group below.
    test::group ("Saturator: reject process before / after failed prepare; oversized blocks chunk safely");
    {
        saturation::Saturator sat;                                   // NOT prepared (channels_ == 0)
        float a[32] {}, b[32] {}; float* io[2] { a, b };
        sat.process (io, 2, 16);                                     // unprepared → safe no-op
        test::ok (! sat.prepare (48000.0, 16, 2, 4, 2), "prepare(tapsPerPhase=2) fails (oversampler rejects <4)");
        sat.process (io, 2, 16);                                     // FAILED prepare → no-op, must not index empty osBuf_
        test::ok (sat.prepare (48000.0, 16, 2, 4, 32), "prepare valid");
        sat.process (io, 2, 16);                                     // works
        sat.process (io, 2, 32);                                     // n=32 > maxBlock=16 → chunked 16+16, must not overrun osBuf_
        test::ok (true, "no OOB across failed-prepare / oversized-block process (ASan/UBSan is the real check)");
    }

    // --- FALSIFICATION: at os>1 the dry/wet mix must be latency-aligned (undelayed dry combs the wet) ---
    test::group ("Saturator os=4 dry/wet mix is latency-aligned (no comb)");
    {
        const int n = 4096;
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (0.5 * std::sin (2.0 * pi * 1000.0 * i / 48000.0));
        std::vector<float> y = x; float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p;
        p.driveDb = 0.0f; p.autoComp = 0.0f; p.mix = 0.5f; p.outputDb = 0.0f;   // ≈linear curve → wet ≈ delayed dry
        sat.setParams (p);
        sat.process (ch, 1, n);
        auto rmsHalf = [] (const std::vector<float>& v) {
            double s2 = 0.0; const int from = (int) v.size() / 2;
            for (int i = from; i < (int) v.size(); ++i) s2 += (double) v[i] * v[i];
            return std::sqrt (s2 / ((int) v.size() - from)); };
        // aligned: 0.5·x + 0.5·x == x; unaligned by 31 samples @1 kHz: |0.5 + 0.5·e^{-jθ}| ≈ 0.44
        test::approx (rmsHalf (y) / rmsHalf (x), 1.0, 0.05, "mix=0.5 at ~zero drive is level-neutral (dry delayed to match wet)");
    }

    // --- autoComp=1 → unity small-signal gain regardless of drive ---
    test::group ("Saturator autoComp=1 small-signal unity");
    {
        const int n = 4096;
        std::vector<float> x (n);
        for (int i = 0; i < n; ++i) x[i] = (float) (0.005 * std::sin (2.0 * pi * 500.0 * i / 48000.0));
        std::vector<float> y = x; float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p; p.driveDb = 12.0f; p.autoComp = 1.0f; p.mix = 1.0f;
        sat.setParams (p);
        sat.process (ch, 1, n);
        double sx = 0.0, sy = 0.0;
        for (int i = n / 2; i < n; ++i) { sx += (double) x[i] * x[i]; sy += (double) y[i] * y[i]; }
        test::approx (std::sqrt (sy / sx), 1.0, 0.02, "12 dB drive + autoComp 1 → tiny signal passes at unity");
    }

    // --- non-finite params must not poison the stream (house rule: clamp + reject non-finite) ---
    test::group ("Saturator non-finite params rejected");
    {
        const int n = 512;
        std::vector<float> y (n);
        for (int i = 0; i < n; ++i) y[i] = (float) (0.5 * std::sin (2.0 * pi * 997.0 * i / 48000.0));
        float* ch[1] { y.data() };
        saturation::Saturator sat; sat.prepare (48000.0, n, 1, 4);
        saturation::Saturator::Params p;
        p.shape = Shape::Asym;
        const float qnan = std::numeric_limits<float>::quiet_NaN();
        p.bias = qnan; p.mix = qnan; p.outputDb = qnan; p.autoComp = qnan; p.dcBlockHz = qnan; p.driveDb = qnan;
        sat.setParams (p);
        sat.process (ch, 1, n);
        bool finite = true; for (float v : y) finite &= (bool) std::isfinite (v);
        test::ok (finite, "all-NaN params → finite output (fallbacks applied)");
    }

    // --- hardening NULL: the poison guards are bit-transparent on finite, in-range input. Six
    //     deterministic scenarios run through the LIVE class and the frozen PRE-CHANGE copy above;
    //     outputs must be bit-identical (memcmp) — the guards may only ever touch poisoned streams. ---
    test::group ("Saturator hardening NULL: bit-identical to the frozen pre-change engine");
    {
        for (int id = 1; id <= 6; ++id)
        {
            const std::vector<float> now = nulltest::runScenario<saturation::Saturator> (id);
            const std::vector<float> pre = nulltest::runScenario<prechange::Saturator> (id);
            test::ok (now.size() == pre.size() && ! now.empty(),
                      "scenario " + std::to_string (id) + ": output sizes match");
            test::ok (now.size() == pre.size()
                        && std::memcmp (now.data(), pre.data(), now.size() * sizeof (float)) == 0,
                      "scenario " + std::to_string (id) + ": bit-identical output");
        }
    }

    // --- poison recovery: injected NaN/Inf samples must not propagate or lodge in state ---
    test::group ("Saturator poison recovery (NaN/Inf injection)");
    {
        const float qnan = std::numeric_limits<float>::quiet_NaN();
        const float pinf = std::numeric_limits<float>::infinity();

        // Tanh path: all state is FIR/delay memory (finite history) → after one clean block the
        // poisoned twin must reconverge to a never-poisoned twin EXACTLY (bit-identical blocks).
        {
            saturation::Saturator clean, hit;
            clean.prepare (48000.0, 512, 2, 4); hit.prepare (48000.0, 512, 2, 4);
            saturation::Saturator::Params p; p.shape = Shape::Tanh; p.driveDb = 15.0f; p.mix = 0.8f; p.autoComp = 0.5f;
            clean.setParams (p); hit.setParams (p);
            std::vector<float> a0 (512), a1 (512), b0 (512), b1 (512);
            float* ioA[2] { a0.data(), a1.data() }; float* ioB[2] { b0.data(), b1.data() };
            bool allFinite = true, lateEqual = true;
            for (int blk = 0; blk < 8; ++blk)
            {
                for (int i = 0; i < 512; ++i)
                {
                    const int t = blk * 512 + i;
                    const float v = 0.8f * (float) std::sin (2.0 * pi * 331.0 * t / 48000.0);
                    a0[i] = v; a1[i] = -v; b0[i] = v; b1[i] = -v;
                }
                if (blk == 2) { b0[7] = qnan; b0[100] = pinf; b0[255] = -pinf; b1[300] = qnan; }
                clean.process (ioA, 2, 512); hit.process (ioB, 2, 512);
                for (int i = 0; i < 512; ++i) allFinite &= std::isfinite (b0[i]) && std::isfinite (b1[i]);
                if (blk >= 4)
                    lateEqual &= std::memcmp (a0.data(), b0.data(), sizeof (float) * 512) == 0
                              && std::memcmp (a1.data(), b1.data(), sizeof (float) * 512) == 0;
            }
            test::ok (allFinite, "Tanh: output stays finite through and after the poison block");
            test::ok (lateEqual, "Tanh: reconverges BIT-EXACTLY to a never-poisoned twin (FIR/delay state flushed)");
        }

        // Asym path adds the recursive DC blocker (IIR) → poison must still not propagate, and the
        // state difference vs the clean twin must decay away (exponential, ~0.51x per 512 block).
        {
            saturation::Saturator clean, hit;
            clean.prepare (48000.0, 512, 1, 4); hit.prepare (48000.0, 512, 1, 4);
            saturation::Saturator::Params p; p.shape = Shape::Asym; p.driveDb = 12.0f; p.bias = 0.35f; p.mix = 1.0f;
            clean.setParams (p); hit.setParams (p);
            std::vector<float> a (512), b (512);
            float* ioA[1] { a.data() }; float* ioB[1] { b.data() };
            bool allFinite = true; double lateDiff = 0.0;
            for (int blk = 0; blk < 24; ++blk)
            {
                for (int i = 0; i < 512; ++i)
                {
                    const int t = blk * 512 + i;
                    const float v = 0.7f * (float) std::sin (2.0 * pi * 441.0 * t / 48000.0);
                    a[i] = v; b[i] = v;
                }
                if (blk == 2) { b[0] = qnan; b[128] = pinf; b[400] = -pinf; }
                clean.process (ioA, 1, 512); hit.process (ioB, 1, 512);
                for (int i = 0; i < 512; ++i) allFinite &= (bool) std::isfinite (b[i]);
                if (blk >= 20)
                    for (int i = 0; i < 512; ++i) lateDiff = std::max (lateDiff, (double) std::fabs (a[i] - b[i]));
            }
            test::ok (allFinite, "Asym+DC: output stays finite through and after the poison block");
            test::ok (lateDiff < 1.0e-3, "Asym+DC: state difference vs the clean twin decays out (IIR recovers)");
        }
    }

    // --- long-buffer chunking equivalence: one 8192-sample call == manual 512-chunks == the internal
    //     driver chunking an 8192 call against maxBlock=500 (uneven tail chunk) — all bit-identical ---
    test::group ("Saturator chunking equivalence (8192 == chunked, bit-identical)");
    {
        const int N = 8192;
        std::vector<float> x (N);
        nulltest::Rng r (0x00DDBA11u);
        for (int i = 0; i < N; ++i) x[i] = 0.9f * r.next();
        saturation::Saturator::Params p; p.shape = Shape::Tanh; p.driveDb = 14.0f; p.mix = 0.6f; p.autoComp = 0.5f;

        std::vector<float> yA = x, yB = x, yC = x;

        saturation::Saturator sA; sA.prepare (48000.0, N, 1, 4); sA.setParams (p);
        { float* io[1] { yA.data() }; sA.process (io, 1, N); }                       // one monolithic call

        saturation::Saturator sB; sB.prepare (48000.0, N, 1, 4); sB.setParams (p);
        for (int off = 0; off < N; off += 512)                                       // manual 16 x 512 calls
        { float* io[1] { yB.data() + off }; sB.process (io, 1, 512); }

        saturation::Saturator sC; sC.prepare (48000.0, 500, 1, 4); sC.setParams (p); // maxBlock 500 → the driver
        { float* io[1] { yC.data() }; sC.process (io, 1, N); }                       //   chunks 16x500 + 192

        test::ok (std::memcmp (yA.data(), yB.data(), sizeof (float) * N) == 0,
                  "one 8192 call == 16x512 manual calls (state streams across calls)");
        test::ok (std::memcmp (yA.data(), yC.data(), sizeof (float) * N) == 0,
                  "one 8192 call == internal maxBlock=500 chunking (uneven tail fully processed)");
        // falsification vs the OLD behavior (n > maxBlock was a silent no-op): the buffer must differ
        // from the raw input — i.e. the oversized call really processed, it didn't just return.
        test::ok (std::memcmp (yC.data(), x.data(), sizeof (float) * N) != 0,
                  "the oversized call transformed the buffer (old code silently no-opped it)");
    }

    return test::report();
}
