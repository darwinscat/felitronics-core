// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// THE PROOF SUITE for felitronics::limiter::TruePeakLimiter's ceiling. `LimiterTests.cpp` stays the
// fast unit/lifecycle suite; this binary is the numerically heavy one that answers a single question:
//
//     does the delivered true peak stay at or below the ceiling on material that is KNOWN to break a
//     lookahead sample-peak limiter, and is that answer measured by something independent?
//
// WHY IT EXISTS. A product migration rests on this limiter's bound being provable where the
// predecessor's (ffmpeg `alimiter`) is heuristic: measured over a 36-row matrix, that predecessor
// ships true peak above its own ceiling in 17 rows and above 0 dBFS in 3, worst case +2.7 dB on a
// click train under ~+88 dB of makeup — while reporting success. That promise was resting on 19
// assertions.
//
// WHAT THIS SUITE FOUND, STATED UP FRONT BECAUSE IT CHANGES A PRODUCT DECISION AND NOT ONLY A TEST:
//
//   1. The bound is EXACT where the limiter enforces it — on its own oversampled grid. Verified
//      directly (`osGridPeakDb()` below) across every witness and factor: never above the ceiling.
//      That is the algebra, and it is what separates this limiter from a gain-ramp heuristic.
//   2. The DELIVERED bound is not the same statement, and it is not flat. What the limiter enforces on
//      an F x fs grid, a listener (and a conformance meter) hears on the continuum between grid points.
//      TWO mechanisms live in that gap and they behave oppositely:
//        (a) GRID GEOMETRY, closed form. A tone at fs*p/q exceeds the ceiling by exactly
//            -20*log10(cos(pi/M')), M' = the distinct magnitude phases it visits on the grid. Measured
//            +1.250 / +0.302 / +0.076 dB at 2x / 4x / 8x. This is the half that oversampling fixes.
//        (b) GAIN MODULATION, which oversampling does NOT fix. The attack is instantaneous, so the
//            limited product is not band-limited and the downsampler re-band-limits it on the way out.
//            On dense material at a 1 ms release: +1.13 / +0.96 / +0.91 dB at 2x / 4x / 8x — 8x costs
//            twice the work, cuts (a) by four, and does not move this at all. It answers to RELEASE
//            SPEED and to density instead. A single hard edge is nearly free (+0.04 dB), so it is the
//            continuous re-modulation that costs, not any one transition.
//   3. So the shipping claim has to be "provable bound on the internal grid, plus a derate", the derate
//      is a design input rather than an accident, and it is NOT bought by oversampling harder:
//      ~0.50 dB at a release of 50 ms or slower, ~1.0 dB below that, at both 4x and 8x. The cheap
//      lever is a release floor. The module header's "GUARANTEES the output true peak <= ceiling
//      (within the down-sampler's tiny ripple)" overstates what the code proves, and the ripple is not
//      what the gap is made of.
//
// This is still a categorically better bound than the predecessor's: (a) is a closed form that shrinks
// with the factor, (b) is governed by one parameter we choose, and neither depends on the material's
// crest factor or on the makeup gain (the algebra is scale-invariant — verified at +64, +76 and
// +88 dB). The predecessor's +2.7 dB excess is none of those things. The comparison drawn here is
// "we do not carry their defect", never "we match them"; nothing in this file is nulled against
// `alimiter`.
//
// EVERY TRUE-PEAK NUMBER NAMES ITS PATH. The core holds two true-peak designs that disagree by up to
// 0.13 dB on nothing but the sample-grid offset of one band-limited signal, so a bare "dBTP" is not a
// measurement. The primary oracle here is `test::tp::truePeakDbFft` — spectral zero-padding, i.e. the
// definition of band-limited reconstruction rather than an implementation of it — cross-checked by an
// independently designed windowed-sinc reconstruction, and both are validated in-suite against the
// five analytically-known EBU Tech 3341 Table 1 signals before anything else runs.

#include <felitronics_test.h>
#include <truepeak_oracle.h>
#include <truepeak_witnesses.h>

#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace felitronics;
namespace tp  = felitronics::test::tp;
namespace tpw = felitronics::test::tpw;

//==============================================================================
// A knobbed transcription of the limiter's loop. Two jobs, and the first one is what makes the second
// trustworthy: in its FAITHFUL configuration it must reproduce the real class BIT FOR BIT (asserted
// below), which is what stops this reference from drifting into a different limiter that proves
// nothing; and with a knob thrown it must FAIL the suite, which is what stops the suite from being a
// set of assertions that no implementation error could violate. That second property is not
// hypothetical insurance — this sprint already produced a test that passed both before and after the
// behaviour it was written to pin was changed.
//
// It also exposes the one number the real class cannot: the peak ON the internal oversampled grid,
// where the algebraic bound actually lives.
struct Weakening
{
    bool releaseTakesMax    = false;   // `min` -> `max` in the gain law: release outruns the requirement
    bool basebandDomain     = false;   // limit SAMPLE peaks at the base rate — the predecessor's mode
    bool gainOnUndelayed    = false;   // gain applied to the current sample instead of the delayed one
    bool windowOfOne        = false;   // sliding max sees only the newest sample, so lookahead is a delay
    bool insertBeforeExpire = false;   // the deque ordering bug: a decreasing run overwrites the max
    bool perChannelGain     = false;   // unlinked per-channel gain — moves the stereo image
};

// The deque with the ordering bug, kept here rather than in the module so nothing shippable carries a
// sabotage switch. Identical to the real one except that expiry runs AFTER the insert.
class BrokenSlidingMax
{
public:
    void prepare (int maxWindow)
    {
        cap = maxWindow < 1 ? 1 : maxWindow;
        v.assign ((std::size_t) cap, 0.0f); ix.assign ((std::size_t) cap, 0);
        W = cap; reset();
    }
    void reset() noexcept { head = tail = count = 0; n = 0; }
    void setWindow (int w) noexcept { W = w < 1 ? 1 : (w > cap ? cap : w); }
    float push (float x) noexcept
    {
        while (count > 0) { const int b = (tail - 1 + cap) % cap; if (v[(std::size_t) b] <= x) { tail = b; --count; } else break; }
        v[(std::size_t) tail] = x; ix[(std::size_t) tail] = n; tail = (tail + 1) % cap; ++count;
        while (count > 0 && ix[(std::size_t) head] <= n - (std::int64_t) W) { head = (head + 1) % cap; --count; }
        ++n;
        return v[(std::size_t) head];
    }
private:
    int cap = 1, W = 1; std::vector<float> v; std::vector<std::int64_t> ix;
    int head = 0, tail = 0, count = 0; std::int64_t n = 0;
};

class ReferenceLimiter
{
public:
    void prepare (double sampleRate, int maxBlock, int maxChannels, int factor, const Weakening& w)
    {
        wk = w; fs = sampleRate; maxCh = maxChannels; maxBlock_ = maxBlock;
        F  = factor < 2 ? 2 : factor;                                    // the real class clamps identically
        (void) os.prepare (F, maxCh, 32); (void) down.prepare (F, maxCh, 32);
        osBuf.assign ((std::size_t) maxCh, std::vector<float> ((std::size_t) maxBlock * (std::size_t) F, 0.0f));
        osPtrs.assign ((std::size_t) maxCh, nullptr);
        const int maxLookOS = (int) std::ceil (20.0 * 0.001 * fs) * F;
        osDelays.assign ((std::size_t) maxCh, core::DelayLine {});
        for (auto& d : osDelays) d.prepare (maxLookOS);
        (void) slide.prepare (maxLookOS + 1); (void) broken.prepare (maxLookOS + 1);
        perCh.assign ((std::size_t) maxCh, detail_slide {});
        for (auto& s : perCh) s.s.prepare (maxLookOS + 1);
        gridPeak = 0.0; grDb = 0.0f;
    }

    void setParams (double ceilingDbTp, double releaseMs, double lookaheadMs)
    {
        ceilingDb    = ceilingDbTp;
        lookBaseband = (int) std::lround (lookaheadMs * 0.001 * fs);
        // In the baseband weakening there is no oversampled rate, so the lookahead is in BASE samples:
        // scaling it by F there would quietly give the predecessor's structure four times the lookahead
        // it was asked for, and the row that compares the two would be comparing two different settings.
        const int lookOS = wk.basebandDomain ? lookBaseband : lookBaseband * F;
        for (auto& d : osDelays) d.setDelay (lookOS);
        const int W = wk.windowOfOne ? 1 : lookOS + 1;
        slide.setWindow (W); broken.setWindow (W);
        for (auto& s : perCh) s.s.setWindow (W);
        const double t = releaseMs * 0.001 * fs * (double) (wk.basebandDomain ? 1 : F);
        relCoef = (t <= 0.0) ? 0.0f : (float) std::exp (-1.0 / t);
    }

    void process (float* const* channels, int numChannels, int numSamples) noexcept
    {
        const int nc = numChannels < maxCh ? numChannels : maxCh;
        if (nc <= 0 || numSamples <= 0 || numSamples > maxBlock_) return;

        if (wk.basebandDomain)                                            // no oversampling at all
        {
            for (int i = 0; i < numSamples; ++i)
            {
                float linked = 0.0f;
                for (int c = 0; c < nc; ++c) linked = std::max (linked, std::fabs (channels[c][i]));
                step (linked, nc);
                for (int c = 0; c < nc; ++c)
                {
                    const float x = channels[c][i];
                    channels[c][i] = (wk.gainOnUndelayed ? x : osDelays[(std::size_t) c].process (x)) * lastGain;
                }
            }
            return;
        }

        const int osN = numSamples * F;
        for (int c = 0; c < nc; ++c) osPtrs[(std::size_t) c] = osBuf[(std::size_t) c].data();
        os.upsample (channels, nc, numSamples, osPtrs.data());

        for (int i = 0; i < osN; ++i)
        {
            float linked = 0.0f;
            for (int c = 0; c < nc; ++c) linked = std::max (linked, std::fabs (osBuf[(std::size_t) c][(std::size_t) i]));
            if (! wk.perChannelGain) step (linked, nc);

            for (int c = 0; c < nc; ++c)
            {
                float g = lastGain;
                if (wk.perChannelGain)
                {
                    const float a = std::fabs (osBuf[(std::size_t) c][(std::size_t) i]);
                    auto& st = perCh[(std::size_t) c];
                    const float smax = st.s.push (a);
                    double raw = ceilingDb - core::gainToDb ((double) smax);
                    if (raw > 0.0) raw = 0.0;
                    st.gr = std::min ((float) raw, st.gr * relCoef);
                    g = (float) core::dbToGain ((double) st.gr);
                }
                const float x = osBuf[(std::size_t) c][(std::size_t) i];
                const float y = (wk.gainOnUndelayed ? x : osDelays[(std::size_t) c].process (x)) * g;
                osBuf[(std::size_t) c][(std::size_t) i] = y;
                gridPeak = std::max (gridPeak, (double) std::fabs (y));
            }
        }
        down.downsample ((const float* const*) osPtrs.data(), nc, numSamples, channels);
        core::flushDenormal (grDb);
    }

    // Peak ON the internal oversampled grid — where the algebra applies — since prepare().
    double osGridPeakDb() const { return 20.0 * std::log10 (gridPeak > 1e-15 ? gridPeak : 1e-15); }
    double gainReductionDb() const { return grDb; }

private:
    struct detail_slide { limiter::detail::SlidingMax s; float gr = 0.0f; };

    void step (float linked, int) noexcept
    {
        const float smax = wk.insertBeforeExpire ? broken.push (linked) : slide.push (linked);
        double raw = ceilingDb - core::gainToDb ((double) smax);
        if (raw > 0.0) raw = 0.0;
        grDb = wk.releaseTakesMax ? std::max ((float) raw, grDb * relCoef)
                                  : std::min ((float) raw, grDb * relCoef);
        lastGain = (float) core::dbToGain ((double) grDb);
    }

    Weakening wk;
    oversampling::PolyphaseOversampler os, down;
    std::vector<std::vector<float>> osBuf;
    std::vector<float*>             osPtrs;
    std::vector<core::DelayLine>    osDelays;
    limiter::detail::SlidingMax     slide;
    BrokenSlidingMax                broken;
    std::vector<detail_slide>       perCh;
    double fs = 48000.0, ceilingDb = -1.0, gridPeak = 0.0;
    int maxCh = 0, maxBlock_ = 0, F = 4, lookBaseband = 0;
    float relCoef = 0.0f, grDb = 0.0f, lastGain = 1.0f;
};

//==============================================================================
struct Setup { double ceilingDb = -1.0, releaseMs = 50.0, lookaheadMs = 1.0; int factor = 4; };

// Render through the REAL class. The drain matters and is not padding: the limiter delays by
// latencySamples(), so without it the tail of the witness never comes out, the buffer stops
// mid-waveform, and both oracles — which interpolate the periodic extension — would honestly report
// the overshoot of a discontinuity the signal never contained. Blocked at `block` (0 = one shot) so the
// same helper serves the block-invariance NULL. Topology goes to prepare(), parameters after it.
static std::vector<std::vector<float>> renderAt (const std::vector<std::vector<float>>& in, double sr,
                                                 const Setup& s, int block, double* grOut = nullptr)
{
    const int nch = (int) in.size(), n = (int) in[0].size();
    limiter::TruePeakLimiter lim;
    const int maxBlock = block > 0 ? block : n;
    (void) lim.prepare (sr, maxBlock, nch, { s.lookaheadMs, s.factor, 32 });   // topology is prepare-time only now
    limiter::TruePeakLimiterParams p;
    p.ceilingDbTp = s.ceilingDb; p.releaseMs = s.releaseMs;
    lim.setParams (p);
    const int drain = lim.latencySamples() + 64;
    std::vector<std::vector<float>> out ((std::size_t) nch, std::vector<float> ((std::size_t) (n + drain), 0.0f));
    for (int c = 0; c < nch; ++c) std::copy (in[(std::size_t) c].begin(), in[(std::size_t) c].end(), out[(std::size_t) c].begin());
    std::vector<float*> ptrs ((std::size_t) nch);
    for (int i = 0; i < n + drain; i += maxBlock)
    {
        const int k = std::min (maxBlock, n + drain - i);
        for (int c = 0; c < nch; ++c) ptrs[(std::size_t) c] = out[(std::size_t) c].data() + i;
        lim.process (ptrs.data(), nch, k);
    }
    if (grOut) *grOut = lim.gainReductionDb();
    return out;
}

static std::vector<std::vector<float>> renderRef (const std::vector<std::vector<float>>& in, double sr,
                                                  const Setup& s, const Weakening& w, double* gridDb = nullptr)
{
    const int nch = (int) in.size(), n = (int) in[0].size();
    ReferenceLimiter ref;
    (void) ref.prepare (sr, n, nch, s.factor, w);
    ref.setParams (s.ceilingDb, s.releaseMs, s.lookaheadMs);
    const int drain = (32 - 1) + (int) std::lround (s.lookaheadMs * 0.001 * sr) + 64;
    std::vector<std::vector<float>> out ((std::size_t) nch, std::vector<float> ((std::size_t) (n + drain), 0.0f));
    for (int c = 0; c < nch; ++c) std::copy (in[(std::size_t) c].begin(), in[(std::size_t) c].end(), out[(std::size_t) c].begin());
    // The reference is prepared for maxBlock == n, so drive it in n-sized chunks.
    std::vector<float*> ptrs ((std::size_t) nch);
    for (int i = 0; i < n + drain; i += n)
    {
        const int k = std::min (n, n + drain - i);
        for (int c = 0; c < nch; ++c) ptrs[(std::size_t) c] = out[(std::size_t) c].data() + i;
        ref.process (ptrs.data(), nch, k);
    }
    if (gridDb) *gridDb = ref.osGridPeakDb();
    return out;
}

// First index at which two renders differ, or -1 for bit-identical. The comparison is EXACT on
// purpose — every caller here is a NULL test, where "close enough" would let a state bug through, so
// the repo's downstream -Wfloat-equal is being contradicted deliberately rather than by accident.
static int firstDifference (const std::vector<float>& a, const std::vector<float>& b)
{
    const std::size_t n = std::min (a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) if (a[i] != b[i]) return (int) i;
    return a.size() == b.size() ? -1 : (int) n;
}

static std::string dbs (double v) { char b[32]; std::snprintf (b, sizeof b, "%+.4f", v); return b; }

//==============================================================================
int main()
{
    std::printf ("felitronics::limiter true-peak CEILING PROOF\n");
    const double sr = 48000.0;

    // ---------------------------------------------------------------- the oracle, before anything else
    // An oracle that has drifted moves every budget in this file silently, so it is checked against
    // external, analytically-known values first. EBU Tech 3341 Table 1's signals are continuous sines:
    // their reconstruction maximum IS the amplitude, so the truth here is arithmetic, not a golden.
    test::group ("Oracles validated against EBU Tech 3341 Table 1 (analytic truth)");
    for (const auto& c : test::ebu3341::kTruePeakCases)
    {
        std::vector<float> s;
        test::ebu3341::synthesize (c, sr, 0.10, s);
        const double truth = test::ebu3341::oracleDb (c);
        test::approx (tp::truePeakDbFft (s, 16),  truth, 0.001, "FFT zero-pad oracle == analytic, EBU #" + std::to_string (c.number));
        test::approx (tp::truePeakDbSinc (s),     truth, 0.001, "windowed-sinc cross-check == analytic, EBU #" + std::to_string (c.number));
    }
    // The five Table 1 phases (0, 45, 60, 67.5 degrees) all land ON the oracle's own interpolation grid,
    // so the agreement above validates its scaling and not its resolution. One tone whose denominator is
    // coprime to the pad does the other half — and its residual has a derivable bound, the oracle's own
    // grid term, rather than a picked number.
    {
        const int n = 4096;
        std::vector<float> t ((std::size_t) n);
        for (int i = 0; i < n; ++i)
            t[(std::size_t) i] = (float) (tpw::taper (i, n, sr) * 0.5 * std::sin (2.0 * core::kPi * (double) i / 7.0 + 0.3));
        const double bound = -20.0 * std::log10 (std::cos (core::kPi / (double) tpw::gridPhaseCount (1, 7, 32)));
        const double err = 20.0 * std::log10 (0.5) - tp::truePeakDbFft (t);
        test::ok (err >= -0.001 && err <= bound + 0.002,
                  "on fs/7, whose phases do NOT lie on the oracle's grid, it under-reads by " + dbs (err)
                  + " dB — within its own derived grid bound of " + dbs (bound));
    }

    // And the status-quo reading, which this suite deliberately does NOT use as its judge: a fresh 8x
    // instance of the limiter's OWN oversampler class. It under-reads, i.e. it errs toward PASS.
    {
        std::vector<float> s;
        test::ebu3341::synthesize (test::ebu3341::kTruePeakCases[1], sr, 0.10, s);   // #16, fs/4 phase 45
        oversampling::PolyphaseOversampler m; (void) m.prepare (8, 1, 32);
        std::vector<float> osb (s.size() * 8);
        const float* xi[1] { s.data() }; float* oo[1] { osb.data() };
        m.upsample (xi, 1, (int) s.size(), oo);
        double mx = 0.0; for (float v : osb) mx = std::max (mx, (double) std::fabs (v));
        const double polyDb = 20.0 * std::log10 (mx);
        test::ok (polyDb < test::ebu3341::oracleDb (test::ebu3341::kTruePeakCases[1]) - 0.01,
                  "[TP path: PolyphaseOversampler 8x/256-tap/0.90-Nyquist] under-reads the truth by "
                  + dbs (test::ebu3341::oracleDb (test::ebu3341::kTruePeakCases[1]) - polyDb) + " dB — why it is not the judge here");
    }

    // ---------------------------------------------------------------- the witnesses are what they claim
    test::group ("Witness self-checks (a malformed envelope survives every assertion about a maximum)");
    {
        const auto click = tpw::clickTrain (sr, 0.20, 3.0, -6.0, 0.5);
        test::ok (tpw::endsAreSilent (click), "click train begins and ends in silence (both oracles require it)");
        test::ok (tpw::crestFactorDb (click) > 18.0, "click train crest factor > 18 dB (got " + dbs (tpw::crestFactorDb (click)) + ")");
        const double ispExcess = tp::truePeakDbFft (click, 16) - tp::samplePeakDb (click);
        test::ok (ispExcess > 3.0, "click train's inter-sample peaks sit > 3 dB above its sample peaks (got "
                                   + dbs (ispExcess) + ") — without this a sample-peak limiter would pass");
        const auto onGrid = tpw::clickTrain (sr, 0.20, 3.0, -6.0, 0.0);
        test::ok (tp::truePeakDbFft (onGrid, 16) - tp::samplePeakDb (onGrid) < 0.01,
                  "at offset 0 the grid lands ON the crest and the same fixture discriminates nothing — the control");
        const auto tone = tpw::gridTone (sr, 0.05, 1, 3, -6.0, 0.125);
        test::ok (tpw::endsAreSilent (tone), "grid tone begins and ends in silence");
        const auto noise = tpw::denseNoise (sr, 0.10, -12.0);
        test::ok (tpw::endsAreSilent (noise), "dense noise begins and ends in silence");
        test::approx (test::ebu3341::fadeAsymmetry ((long long) (sr * 0.20), (long long) (sr * 0.010)), 0.0, 0.0,
                      "the shared taper's two ramps mirror exactly");
        // The two oracles are independent constructions; on the witness that matters they must agree,
        // or the number this suite spends its budget on is one construction's opinion.
        test::approx (tp::truePeakDbSinc (click), tp::truePeakDbFft (click, 16), 0.01,
                      "the two independent oracles agree on the click train to 0.01 dB");
    }

    // ---------------------------------------------------------------- the reference must BE the limiter
    test::group ("Reference transcription == the real class, bit for bit (else it proves nothing)");
    for (int F : { 2, 4, 8 })
    {
        const auto click = tpw::clickTrain (sr, 0.10, 3.0, 26.0, 0.5);
        std::vector<std::vector<float>> in { click };
        const Setup s { -1.0, 50.0, 1.0, F };
        const auto real = renderAt (in, sr, s, 0);
        const auto ref  = renderRef (in, sr, s, Weakening {});
        const int d = firstDifference (real[0], ref[0]);
        test::ok (d < 0, "F=" + std::to_string (F) + ": faithful reference is bit-identical to TruePeakLimiter"
                         + (d < 0 ? "" : " (first difference at sample " + std::to_string (d) + ")"));
    }

    // ---------------------------------------------------------------- THE ALGEBRA: exact on the grid
    // This is the part that is a theorem. smax over the lookahead window necessarily includes the very
    // sample being emitted, and the release branch can only ever KEEP more reduction than required, so
    // gain <= ceiling/smax and the emitted grid sample cannot exceed the ceiling. Asserted on the
    // reference, which the NULL above ties to the shipping class.
    test::group ("The bound is EXACT on the limiter's own oversampled grid (the theorem)");
    for (int F : { 2, 4, 8 })
    {
        const double C = -2.2;                       // the predecessor's worst row, for continuity of numbers
        for (double makeup : { 64.0, 88.0 })
        {
            const auto click = tpw::clickTrain (sr, 0.10, 3.0, -60.0 + makeup, 0.5);
            std::vector<std::vector<float>> in { click };
            double gridDb = 0.0;
            renderRef (in, sr, Setup { C, 50.0, 1.0, F }, Weakening {}, &gridDb);
            test::ok (gridDb <= C + 0.001, "F=" + std::to_string (F) + ", makeup " + dbs (makeup)
                      + " dB: grid peak " + dbs (gridDb) + " <= ceiling " + dbs (C));
        }
        const auto noise = tpw::denseNoise (sr, 0.10, -7.0);
        std::vector<std::vector<float>> in { noise };
        double gridDb = 0.0;
        renderRef (in, sr, Setup { -1.0, 10.0, 1.0, F }, Weakening {}, &gridDb);
        test::ok (gridDb <= -1.0 + 0.001, "F=" + std::to_string (F) + ", dense noise: grid peak "
                                          + dbs (gridDb) + " <= ceiling -1.0000");
    }

    // ---------------------------------------------------------------- THE DELIVERED CEILING: closed form
    // A two-sided NULL, deliberately. The upper side is the ceiling claim; the LOWER side is what keeps
    // the test honest — a limiter that ducked 3 dB below its ceiling would satisfy any one-sided
    // assertion while being useless, and this suite exists because 19 one-sided assertions were not a
    // safety net.
    test::group ("Delivered true peak == ceiling + the grid's closed form (two-sided, +-0.015 dB)");
    {
        struct Tone { int p, q; const char* name; };
        const Tone tones[] { { 1, 3, "fs/3" }, { 1, 4, "fs/4" }, { 4, 11, "4fs/11" } };
        for (const auto& t : tones)
            for (int F : { 2, 4, 8 })
            {
                const double C = -1.0;
                const int M = tpw::gridPhaseCount (t.p, t.q, F);
                const double want = tpw::gridBreachDb (t.p, t.q, F);
                // The worst phase is FOUND, not derived: the internal grid sits half an oversampled
                // sample off the input grid (the prototype's group delay (N-1)/2 is a half-integer), so
                // p/q alone does not say where the crest lands. A full circle in half-step increments,
                // and the oracle left at its default density — measuring these rows through a pad of 8
                // reported 0.003 dB where the truth is 0.089, because at that density the ORACLE's
                // phases coincide with the limiter's (see truepeak_oracle.h). Both grids have to be
                // placed against the feature, not just the one being tested.
                double worst = -1e9;
                for (int k = 0; k < 2 * M; ++k)
                {
                    const auto x = tpw::gridTone (sr, 0.05, t.p, t.q, 12.0, (double) k / (double) (2 * M));
                    std::vector<std::vector<float>> in { x };
                    const auto y = renderAt (in, sr, Setup { C, 50.0, 1.0, F }, 0);
                    worst = std::max (worst, tp::truePeakDbFft (y[0]));
                }
                test::approx (worst - C, want, 0.015,
                              std::string (t.name) + " at F=" + std::to_string (F) + ": delivered excess == -20log10(cos(pi/"
                              + std::to_string (M) + "))");
            }
    }

    // ---------------------------------------------------------------- the transient witness, all factors
    // The material that breaks the predecessor, at the predecessor's own operating point. The budget is
    // DERIVED — the worst in-band grid excess for the factor, plus one documented modulation allowance —
    // not a number chosen to make the row pass.
    test::group ("Click train under +64/+76/+88 dB makeup: delivered TP within the derived budget");
    {
        const double C = -2.2;
        for (int F : { 2, 4, 8 })
        {
            const double budget = tpw::deliveredBudgetDb (F);
            for (double makeup : { 64.0, 76.0, 88.0 })
            {
                const auto click = tpw::clickTrain (sr, 0.15, 3.0, -60.0 + makeup, 0.5);
                std::vector<std::vector<float>> in { click };
                double gr = 0.0;
                const auto y = renderAt (in, sr, Setup { C, 50.0, 1.0, F }, 0, &gr);
                const double got = tp::truePeakDbFft (y[0], 16);
                test::ok (got <= C + budget, "F=" + std::to_string (F) + " makeup " + dbs (makeup) + " dB: delivered "
                          + dbs (got) + " <= " + dbs (C + budget) + " (GR " + dbs (gr) + " dB)");
                test::ok (got > C - 1.5, "F=" + std::to_string (F) + " makeup " + dbs (makeup)
                          + " dB: and it is not ducking far BELOW the ceiling (" + dbs (got) + ")");
            }
        }
    }

    // ---------------------------------------------------------------- the oracles agree on the OUTPUT too
    // Agreeing on the witnesses is not the same as agreeing on what the limiter did to them: the limited
    // waveform is flatter, its crests sit at nearly equal heights, and that is exactly the shape where
    // the cross-check's top-K candidate restriction is least obviously safe. So it is checked there.
    test::group ("Both oracles agree on LIMITED output (where the cross-check is least comfortable)");
    {
        for (int F : { 2, 4, 8 })
        {
            const auto click = tpw::clickTrain (sr, 0.10, 3.0, 26.0, 0.5);
            std::vector<std::vector<float>> in { click };
            const auto y = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, F }, 0);
            test::approx (tp::truePeakDbSinc (y[0]), tp::truePeakDbFft (y[0]), 0.02,
                          "F=" + std::to_string (F) + ": windowed-sinc cross-check == FFT oracle on limited material");
        }
    }

    // ---------------------------------------------------------------- the OTHER mechanism, measured
    // Dense broadband material is where the modulation term lives, and it is the term that does NOT
    // shrink with the oversampling factor — which is why the budget carries it separately from the
    // grid's closed form. Asserted two-sided: the upper side is the ceiling claim, the lower side
    // records that the mechanism is really there, so a future change that quietly removes the excess
    // has to come here and say so.
    test::group ("Dense material: the modulation term, which oversampling does not buy off");
    {
        const double C = -1.0;
        double deepAt4 = 0.0, deepAt8 = 0.0;
        // At DEPTH, and at a slow release. An earlier draft ran these at about 1 dB of gain reduction
        // and at that depth the term is still climbing, which is how a release-dependent budget got
        // asserted for a mechanism that turns out to saturate instead.
        for (double rel : { 100.0, 50.0, 1.0 })
            for (int F : { 2, 4, 8 })
            {
                const auto noise = tpw::denseNoise (sr, 0.20, 10.0);            // ~11 dB of reduction
                std::vector<std::vector<float>> in { noise };
                const auto y = renderAt (in, sr, Setup { C, rel, 1.0, F }, 0);
                const double over = tp::truePeakDbFft (y[0]) - C;
                test::ok (over <= tpw::deliveredBudgetDb (F),
                          "F=" + std::to_string (F) + " release " + dbs (rel) + " ms: dense excess " + dbs (over)
                          + " <= budget " + dbs (tpw::deliveredBudgetDb (F)));
                if (rel >= 50.0 && F == 4) deepAt4 = std::max (deepAt4, over);
                if (rel >= 50.0 && F == 8) deepAt8 = std::max (deepAt8, over);
            }
        // Saturation: a 100 ms release does not rescue it. Asserted because the opposite was the
        // intuitive conclusion and it was wrong.
        test::ok (deepAt4 > 0.7, "at depth the excess survives a 50-100 ms release (" + dbs (deepAt4)
                                 + " dB) — a release floor is not the lever it looks like");
        // And the one thing that DOES buy it off: spectral tilt. This row is the realistic-master datum,
        // and the reason the worst case above is not the number a normal mix would see.
        {
            const auto pink = tpw::denseNoise (sr, 0.20, 10.0, 0x9E3779B97F4A7C15ULL, 0.9);
            std::vector<std::vector<float>> in { pink };
            const auto y = renderAt (in, sr, Setup { C, 50.0, 1.0, 4 }, 0);
            const double over = tp::truePeakDbFft (y[0]) - C;
            test::ok (over < 0.25, "the same depth on roughly music-like material costs only " + dbs (over)
                                   + " dB, against " + dbs (deepAt4) + " on flat noise");
        }
        // Click train at the fastest release a caller might legitimately pick.
        for (int F : { 2, 4, 8 })
        {
            const auto click = tpw::clickTrain (sr, 0.15, 3.0, 26.0, 0.5);
            std::vector<std::vector<float>> in { click };
            const auto y = renderAt (in, sr, Setup { C, 0.1, 1.0, F }, 0);
            const double over = tp::truePeakDbFft (y[0]) - C;
            test::ok (over <= tpw::deliveredBudgetDb (F), "F=" + std::to_string (F)
                      + ", click at a 0.1 ms release: " + dbs (over) + " <= " + dbs (tpw::deliveredBudgetDb (F)));
        }
        // THE claim, asserted rather than described: 8x costs twice the work of 4x and cuts the grid
        // term by a factor of four, and it does not move this excess at all. So a delivery promise
        // cannot be bought with oversampling — it has to be bought with a release floor or a derate.
        test::approx (deepAt8, deepAt4, 0.15, "dense excess at 8x == at 4x (" + dbs (deepAt8) + " vs " + dbs (deepAt4)
                                              + ") while the grid term falls from " + dbs (tpw::gridBreachDb (1, 3, 4))
                                              + " to " + dbs (tpw::gridBreachDb (1, 3, 8)) + " dB");
        test::ok (deepAt4 > 4.0 * tpw::gridBreachDb (1, 3, 8),
                  "and at 8x it is " + dbs (deepAt8 / tpw::gridBreachDb (1, 3, 8))
                  + "x the entire grid term — this mechanism, not geometry, sets the derate");
        // One hard transition, by contrast, is nearly free — so the cost above is the CONTINUOUS
        // re-modulation of dense material, not any single edge. Worth pinning: it is the difference
        // between "add a release floor" and "redesign the gain path".
        for (double step : { 12.0, 60.0 })
        {
            const auto edge = tpw::hardEdge (sr, 0.10, 997.0, -1.0, step, 0.05);
            std::vector<std::vector<float>> in { edge };
            const auto y = renderAt (in, sr, Setup { -1.0, 1.0, 1.0, 4 }, 0);
            const double over = tp::truePeakDbFft (y[0]) - (-1.0);
            test::ok (over < 0.10, "a single " + dbs (step) + " dB edge at a 1 ms release costs only "
                                   + dbs (over) + " dB — density is what costs, not transitions");
        }
    }

    // ---------------------------------------------------------------- lookahead 0 is a legal parameter
    // And it removes the bound: with no delay there is nothing for the sliding maximum to look ahead
    // OF, so the structure degenerates into a clipper at the oversampled rate. `lookaheadMs` is not
    // clamped to a minimum, so a caller can reach this by asking politely. Pinned as behaviour, with
    // the number, because the composite must not ship it.
    // WAS a pinned defect: lookaheadMs = 0 was legal and left the bound behind (+2.05 dB over) — with
    // no delay there is nothing to look ahead OF, and the structure degenerated into a clipper at the
    // oversampled rate. Anything rounding to zero baseband samples did it, so the floor is not a
    // millisecond number: it is TWO BASEBAND SAMPLES, the smallest value at which every witness lands
    // back inside the delivered envelope at both conformance rates and every factor (worst margin
    // 0.35 dB). At 48 kHz that is 42 microseconds — invisible at any musical setting.
    test::group ("lookaheadMs below the measured floor is clamped, visibly");
    {
        const double C = -1.0;
        const auto click = tpw::clickTrain (sr, 0.10, 3.0, 26.0, 0.5);
        std::vector<std::vector<float>> in { click };
        for (int F : { 2, 4, 8 })
        {
            const auto none = renderAt (in, sr, Setup { C, 50.0, 0.0, F }, 0);
            const double over = tp::truePeakDbFft (none[0]) - C;
            test::ok (over <= tpw::deliveredBudgetDb (F), "F=" + std::to_string (F)
                      + ": a requested lookahead of 0 delivers " + dbs (over) + " dB over, inside the budget"
                      " (it used to be +2.05)");
        }
        limiter::TruePeakLimiter lim;
        (void) lim.prepare (sr, 512, 1, { 0.0, 4, 32 });
        test::ok (lim.lookaheadSamples() == 2, "and the clamp is not silent: lookaheadSamples() reports "
                                               + std::to_string (lim.lookaheadSamples()));
        limiter::TruePeakLimiter big;
        (void) big.prepare (44100.0, 512, 1, { 0.0, 4, 32 });
        test::ok (big.lookaheadSamples() == 2, "the floor is in SAMPLES, so it is the same at 44.1 kHz");
        // The other end, and it is not symmetric bookkeeping: an absurd FINITE lookahead used to reach
        // std::lround() before the clamp, and the out-of-range result then clamped UPWARD to the
        // minimum — a caller asking for the largest lookahead silently got the smallest.
        limiter::TruePeakLimiter huge;
        (void) huge.prepare (sr, 512, 1, { 1.0e300, 4, 32 });
        test::ok (huge.lookaheadSamples() == (int) std::lround (20.0 * 0.001 * sr),
                  "an absurd finite lookahead clamps to the 20 ms maximum, not to the minimum (got "
                  + std::to_string (huge.lookaheadSamples()) + ")");
    }

    // ---------------------------------------------------------------- releaseMs = 0 is legal too
    // The same shape of hazard as lookahead 0, found by an adversarial witness rather than by reading:
    // at releaseMs = 0 the release coefficient is 0, so the gain is free to jump on every oversampled
    // sample, and the downsampler's signed step response then rides on top of a ceiling-flat plateau
    // instead of a decaying tail. Measured on the plateau witness, edge offset swept: +1.48 dB over at
    // 4x. One tenth of a millisecond of release already halves it. The budget above is stated for
    // release > 0 and this is why.
    test::group ("releaseMs below the measured floor is clamped, visibly");
    {
        // WAS a pinned defect of the same family, found by an adversarial witness rather than by reading:
        // at releaseMs = 0 the release coefficient is 0, the gain is free to jump on every oversampled
        // sample, and the downsampler's step response then rides on a ceiling-flat plateau (+1.48 dB
        // over). The floor is EIGHT baseband samples — the smallest at which every witness, including the
        // click train at that release, lands back inside the envelope (worst margin 0.23 dB). 167 us at
        // 48 kHz, far below any musical release.
        const double C = -1.0;
        for (int F : { 2, 4, 8 })
        {
            double worst = -1e9;
            for (int k = 0; k < 8; ++k)
            {
                const auto x = tpw::plateau (sr, 0.17, 60.0, 18.0, 0.042, 0.104, (double) k / 8.0);
                std::vector<std::vector<float>> in { x };
                worst = std::max (worst, tp::truePeakDbFft (renderAt (in, sr, Setup { C, 0.0, 1.0, F }, 0)[0]));
            }
            test::ok (worst - C <= tpw::deliveredBudgetDb (F), "F=" + std::to_string (F)
                      + ": a requested release of 0 delivers " + dbs (worst - C) + " dB over, inside the budget"
                      " (it used to be +1.48)");
            const auto click = tpw::clickTrain (sr, 0.12, 3.0, 26.0, 0.5);
            std::vector<std::vector<float>> ci { click };
            const double c2 = tp::truePeakDbFft (renderAt (ci, sr, Setup { C, 0.0, 1.0, F }, 0)[0]) - C;
            test::ok (c2 <= tpw::deliveredBudgetDb (F), "F=" + std::to_string (F)
                      + ": and the click train at that release too (" + dbs (c2) + ")");
        }
        limiter::TruePeakLimiter lim;
        (void) lim.prepare (sr, 512, 1, { 1.0, 4, 32 });
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = 0.0;
        lim.setParams (p);
        // Tolerance, not exactness, and deliberately: effectiveReleaseMs() is derived BACK from the
        // float coefficient the limiter actually uses, so it reports what the audio path does rather
        // than what the caller asked for. The round trip through float costs a few parts per million.
        test::approx (lim.effectiveReleaseMs(), 8000.0 / sr, 1.0e-4,
                      "and the clamp is visible: effectiveReleaseMs() reports " + dbs (lim.effectiveReleaseMs()) + " ms");
    }

    // ---------------------------------------------------------------- scale invariance of the excess
    // The predecessor collapses specifically under enormous makeup. This limiter's algebra is
    // scale-invariant, and that is a claim worth pinning rather than assuming: the same excess must
    // come out at +64 and at +88 dB.
    test::group ("The excess does not grow with makeup (the predecessor's failure mode is absent)");
    for (int F : { 2, 4, 8 })
    {
        const double C = -2.2;
        double lo = 0.0, hi = 0.0;
        for (int j = 0; j < 2; ++j)
        {
            const auto click = tpw::clickTrain (sr, 0.15, 3.0, -60.0 + (j == 0 ? 64.0 : 88.0), 0.5);
            std::vector<std::vector<float>> in { click };
            const auto y = renderAt (in, sr, Setup { C, 50.0, 1.0, F }, 0);
            (j == 0 ? lo : hi) = tp::truePeakDbFft (y[0], 16);
        }
        test::approx (hi, lo, 0.02, "F=" + std::to_string (F) + ": +88 dB of makeup delivers the same peak as +64 dB");
    }

    // ---------------------------------------------------------------- "we do not carry their defect"
    // Not a null against ffmpeg — the reference's baseband mode is the STRUCTURE the predecessor uses
    // (bound the samples, not the reconstruction), and the point is the size of the gap between that
    // structure and this one on the same witness.
    test::group ("Sample-peak limiting (the predecessor's structure) breaches where ours holds");
    {
        const double C = -2.2;
        const auto click = tpw::clickTrain (sr, 0.15, 3.0, 26.0, 0.5);
        std::vector<std::vector<float>> in { click };
        Weakening base; base.basebandDomain = true;
        const auto bad = renderRef (in, sr, Setup { C, 50.0, 1.0, 4 }, base);
        const double badTp = tp::truePeakDbFft (bad[0], 16);
        const auto good = renderAt (in, sr, Setup { C, 50.0, 1.0, 4 }, 0);
        const double goodTp = tp::truePeakDbFft (good[0], 16);
        test::ok (badTp > C + 2.5, "a sample-peak limiter ships " + dbs (badTp - C) + " dB over the same ceiling");
        test::ok (goodTp < badTp - 2.0, "the true-peak limiter is at least 2 dB tighter on the identical witness ("
                                        + dbs (goodTp - C) + " vs " + dbs (badTp - C) + " dB over)");
    }

    // ---------------------------------------------------------------- ROLLBACK: each knob must break it
    // Run at a FAST release on purpose. At release 100 ms a one-sample sliding window nearly passes
    // (the release cannot climb far in one oversampled step, so the missing window is almost invisible);
    // at 1 ms it is unmistakable. A negative control that only bites at one setting is worth knowing
    // about, so the setting is named here rather than discovered later.
    test::group ("Falsification: an artificially weakened limiter FAILS this suite");
    {
        const double C = -2.2;
        const auto click = tpw::clickTrain (sr, 0.15, 3.0, 26.0, 0.5);
        std::vector<std::vector<float>> mono { click };
        const Setup fast { C, 1.0, 1.0, 4 };
        const double budget = tpw::deliveredBudgetDb (4);

        struct Knob { const char* name; Weakening w; };
        const Knob knobs[] {
            { "release takes max (gain outruns the requirement)", [] { Weakening w; w.releaseTakesMax = true; return w; }() },
            { "limiting in the sample domain, not the OS domain", [] { Weakening w; w.basebandDomain = true; return w; }() },
            { "gain applied to the undelayed sample",             [] { Weakening w; w.gainOnUndelayed = true; return w; }() },
            { "sliding-max window of one",                        [] { Weakening w; w.windowOfOne = true; return w; }() },
        };
        for (const auto& k : knobs)
        {
            const auto y = renderRef (mono, sr, fast, k.w);
            const double got = tp::truePeakDbFft (y[0], 16) - C;
            test::ok (got > budget, std::string ("caught: ") + k.name + " ships " + dbs (got)
                                    + " dB over, past the " + dbs (budget) + " dB budget");
        }

        // The deque ordering bug needs a strictly-decreasing run LONGER than the deque, which only
        // happens when the window is wide. At the 1 ms lookahead every audio witness in this file is
        // blind to it — recorded, because a negative control that silently does not bite is exactly the
        // failure this group exists to prevent.
        {
            const int n = (int) (sr * 0.06);
            std::vector<float> ramp ((std::size_t) n, 0.0f);
            const int len = (int) (sr * 0.030);                        // 30 ms decreasing > the 20 ms window
            for (int i = 0; i < len; ++i) ramp[(std::size_t) i] = (float) (4.0 * (1.0 - 0.7 * (double) i / (double) len));
            for (int i = 0; i < n; ++i) ramp[(std::size_t) i] = (float) (ramp[(std::size_t) i] * tpw::taper (i, n, sr));
            std::vector<std::vector<float>> in { ramp };
            const Setup wide { C, 0.5, 20.0, 4 };
            Weakening w; w.insertBeforeExpire = true;
            const auto bad  = renderRef (in, sr, wide, w);
            const auto good = renderRef (in, sr, wide, Weakening {});
            const double badTp  = tp::truePeakDbFft (bad[0], 16);
            const double goodTp = tp::truePeakDbFft (good[0], 16);
            test::ok (badTp > C + budget, "caught: insert-before-expire ships " + dbs (badTp - C)
                                          + " dB over — but ONLY on a decreasing run at a 20 ms lookahead");
            test::ok (goodTp <= C + budget, "and the faithful deque holds on that same ramp (" + dbs (goodTp - C) + " dB)");
        }

        // Unlinked gain is invisible to every mono witness; it needs a deliberately asymmetric pair.
        {
            auto l = tpw::clickTrain (sr, 0.10, 3.0, 26.0, 0.5);
            std::vector<float> r ((std::size_t) l.size());
            for (std::size_t i = 0; i < l.size(); ++i) r[i] = l[i] * 0.01f;      // R sits exactly 40 dB below L
            std::vector<std::vector<float>> in { l, r };
            Weakening w; w.perChannelGain = true;
            const auto bad  = renderRef (in, sr, Setup { C, 50.0, 1.0, 4 }, w);
            const auto good = renderRef (in, sr, Setup { C, 50.0, 1.0, 4 }, Weakening {});
            const double badRatio  = tp::truePeakDbFft (bad[1], 16)  - tp::truePeakDbFft (bad[0], 16);
            const double goodRatio = tp::truePeakDbFft (good[1], 16) - tp::truePeakDbFft (good[0], 16);
            test::approx (goodRatio, -40.0, 0.5, "linked gain preserves the 40 dB L/R relationship (image intact)");
            test::ok (badRatio > goodRatio + 10.0, "caught: unlinked gain moves the image by "
                                                   + dbs (badRatio - goodRatio) + " dB");
        }
    }

    // ---------------------------------------------------------------- linking, on the SHIPPING class
    // The negative control above compares two configurations of the reference, which cannot notice a
    // detector that stopped looking at channel 1. This one drives the real limiter with a deliberately
    // asymmetric pair: the loud channel must duck the quiet one by the same gain, so the 40 dB
    // relationship between them survives. A per-channel detector would hold the quiet channel at unity
    // and move the image by tens of dB.
    test::group ("Channel linking on the real class (the image does not move)");
    {
        auto loud = tpw::clickTrain (sr, 0.10, 3.0, 26.0, 0.5);
        std::vector<float> quiet ((std::size_t) loud.size());
        for (std::size_t i = 0; i < loud.size(); ++i) quiet[i] = loud[i] * 0.01f;
        // BOTH orientations, and that is not symmetry for its own sake: with the loud channel first, a
        // detector that only ever reads channel 0 still sees the loud one and produces the identical
        // output. Measured — the one-orientation version of this group passed that mutation.
        for (int which = 0; which < 2; ++which)
        {
            std::vector<std::vector<float>> in = which == 0 ? std::vector<std::vector<float>> { loud, quiet }
                                                            : std::vector<std::vector<float>> { quiet, loud };
            const char* name = which == 0 ? "loud on the left" : "loud on the RIGHT (channel 0 is quiet)";
            for (int F : { 2, 4, 8 })
            {
                const auto y = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, F }, 0);
                const double lo = tp::truePeakDbFft (y[which == 0 ? 1 : 0]);
                const double hi = tp::truePeakDbFft (y[which == 0 ? 0 : 1]);
                test::approx (lo - hi, -40.0, 0.5, std::string (name) + ", F=" + std::to_string (F)
                                                   + ": the 40 dB relationship survives limiting (" + dbs (lo - hi) + ")");
                test::ok (hi <= -1.0 + tpw::deliveredBudgetDb (F),
                          std::string (name) + ", F=" + std::to_string (F) + ": and the loud channel IS limited ("
                          + dbs (hi) + " dBFS)");
            }
        }
        // And a stereo-prepared limiter handed a single channel must limit it, not walk off the end of
        // the pointer array — a host that drops to mono mid-session is not exotic.
        limiter::TruePeakLimiter lim;
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0;
        (void) lim.prepare (sr, (int) loud.size(), 2, { 1.0, 4, 32 }); lim.setParams (p);
        std::vector<float> mono = loud; float* one[1] { mono.data() };
        lim.process (one, 1, (int) loud.size());
        test::ok (tp::samplePeakDb (mono) < -0.5, "a stereo-prepared limiter still limits a mono call ("
                                                  + dbs (tp::samplePeakDb (mono)) + " dBFS)");
    }

    // ---------------------------------------------------------------- parameters actually take effect
    // Everything above that measures the shipping class does it at one operating point: 48 kHz, 1 ms
    // lookahead, default taps. That is enough to prove a ceiling and not enough to prove the parameters
    // mean anything — a limiter that silently hard-coded any of the three would sail through. These are
    // the cheap assertions that close it.
    test::group ("Parameters are not decorative: rate, lookahead and taps all reach the state");
    {
        struct Row { double rate; double lookMs; int taps; };
        const Row rows[] { { 48000.0, 1.0, 32 }, { 44100.0, 1.0, 32 }, { 48000.0, 20.0, 32 },
                           { 48000.0, 1.0, 16 }, { 48000.0, 1.0, 64 } };
        for (const auto& r : rows)
        {
            limiter::TruePeakLimiter lim;
            limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0;
            (void) lim.prepare (r.rate, 512, 1, { r.lookMs, 4, r.taps }); lim.setParams (p);
            const int want = (r.taps - 1) + (int) std::lround (r.lookMs * 0.001 * r.rate);
            test::ok (lim.latencySamples() == want,
                      std::to_string ((int) r.rate) + " Hz, " + dbs (r.lookMs) + " ms, " + std::to_string (r.taps)
                      + " taps: latency " + std::to_string (lim.latencySamples()) + " == " + std::to_string (want));
        }
        // Release is a time, not a flag: the same excursion must recover more slowly at 500 ms than at
        // 50 ms, measured at a fixed distance after the event.
        const int n = (int) (sr * 0.25);
        std::vector<float> x ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (0.2 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));
        for (int i = 1000; i < 1200; ++i) x[(std::size_t) i] *= 40.0f;               // a loud excursion
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (x[(std::size_t) i] * tpw::taper (i, n, sr));
        std::vector<std::vector<float>> in { x };
        auto levelAt = [] (const std::vector<float>& v, int a, int b)
        { double s2 = 0.0; for (int i = a; i < b; ++i) s2 += (double) v[(std::size_t) i] * v[(std::size_t) i];
          return 20.0 * std::log10 (std::sqrt (s2 / (double) (b - a)) + 1e-30); };
        // A 20 ms lookahead must not merely change the reported latency — the sliding window has to be
        // that wide in the audio path too. Nothing else in this suite runs the SHIPPING class at a
        // lookahead above 1 ms, so a deque sized one slot short at the 20 ms cap would be invisible.
        {
            const int M = tpw::gridPhaseCount (1, 3, 4);
            double worst = -1e9;
            for (int k = 0; k < 2 * M; ++k)
            {
                const auto t = tpw::gridTone (sr, 0.08, 1, 3, 12.0, (double) k / (double) (2 * M));
                std::vector<std::vector<float>> ti { t };
                worst = std::max (worst, tp::truePeakDbFft (renderAt (ti, sr, Setup { -1.0, 1.0, 20.0, 4 }, 0)[0]));
            }
            test::approx (worst + 1.0, tpw::gridBreachDb (1, 3, 4), 0.02,
                          "at a 20 ms lookahead and a 1 ms release the bound still holds exactly (" + dbs (worst + 1.0) + ")");
        }
        // More channels than were prepared must be refused, not walked past.
        {
            limiter::TruePeakLimiter narrow;
            limiter::TruePeakLimiterParams np; np.ceilingDbTp = -1.0;
            (void) narrow.prepare (sr, 256, 1, { 1.0, 4, 32 }); narrow.setParams (np);
            std::vector<float> ca (256, 0.9f), cb (256, 0.9f);
            float* two[2] { ca.data(), cb.data() };
            narrow.process (two, 2, 256);                                 // prepared for 1, handed 2
            test::ok (cb[128] == 0.9f, "a mono-prepared limiter leaves the extra channel untouched rather than writing it");
        }
        // setParams AFTER prepare has to reach the state as well; every render helper here sets before
        // preparing, which would let a setParams that stored its argument and forgot to apply it pass.
        {
            limiter::TruePeakLimiter late;
            limiter::TruePeakLimiterParams lp; lp.ceilingDbTp = -1.0;
            (void) late.prepare (sr, 1024, 1, { 1.0, 4, 32 });
            lp.ceilingDbTp = -12.0; lp.releaseMs = 50.0;
            late.setParams (lp);                                          // only NOW is the ceiling set
            std::vector<float> hot (1024);
            for (int i = 0; i < 1024; ++i) hot[(std::size_t) i] = (float) (0.9 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));
            float* hp[1] { hot.data() };
            late.process (hp, 1, 1024);
            test::ok (tp::samplePeakDb (hot) < -11.0, "a ceiling set after prepare() is the one that applies ("
                                                      + dbs (tp::samplePeakDb (hot)) + " dBFS)");
        }
        const auto fastR = renderAt (in, sr, Setup { -1.0, 50.0,  1.0, 4 }, 0);
        const auto slowR = renderAt (in, sr, Setup { -1.0, 500.0, 1.0, 4 }, 0);
        const int a = 3000, b = 5000;                                                // ~40-100 ms after the burst
        test::ok (levelAt (slowR[0], a, b) < levelAt (fastR[0], a, b) - 1.0,
                  "release 500 ms is still " + dbs (levelAt (fastR[0], a, b) - levelAt (slowR[0], a, b))
                  + " dB further down than 50 ms, 40 ms after the burst");
    }

    // ---------------------------------------------------------------- coverage the module simply lacked
    test::group ("Block-size invariance is BIT-EXACT (streaming state carries no block dependence)");
    {
        const auto click = tpw::clickTrain (sr, 0.10, 3.0, 20.0, 0.5);
        std::vector<std::vector<float>> in { click };
        const Setup s { -1.0, 50.0, 1.0, 4 };
        const auto whole = renderAt (in, sr, s, 0);
        for (int b : { 1, 7, 64, 512 })
        {
            const auto part = renderAt (in, sr, s, b);
            const int d = firstDifference (whole[0], part[0]);
            test::ok (d < 0, "block " + std::to_string (b) + " == one-shot, sample for sample"
                             + (d < 0 ? "" : " (differs at " + std::to_string (d) + ")"));
        }
    }

    test::group ("Requested 1x oversampling is the clamped 2x, bit for bit (there is no 1x path)");
    {
        const auto click = tpw::clickTrain (sr, 0.05, 3.0, 20.0, 0.5);
        std::vector<std::vector<float>> in { click };
        const auto one = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, 1 }, 0);
        const auto two = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, 2 }, 0);
        test::ok (firstDifference (one[0], two[0]) < 0,
                  "factor 1 renders identically to factor 2 — the predecessor's os1 mode cannot be reproduced here");
    }

    test::group ("Gain never rises above unity, and deep reduction recovers");
    {
        // Not `gainReductionDb() <= 0`, which is true by construction and therefore proves nothing:
        // the property worth asserting is that the limiter never makes anything LOUDER than it arrived,
        // on a witness that spends most of its length in recovery from a deep excursion.
        const auto click = tpw::clickTrain (sr, 0.30, 5.0, 28.0, 0.5);
        std::vector<std::vector<float>> in { click };
        double gr = 0.0;
        const auto y = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, 4 }, 0, &gr);
        test::ok (tp::truePeakDbFft (y[0]) <= tp::truePeakDbFft (click) + 0.001,
                  "output true peak never exceeds the input's — no upward gain anywhere in the recovery");
        const auto quiet = tpw::gridTone (sr, 0.30, 1, 48, -20.0, 0.0);
        std::vector<std::vector<float>> qin { quiet };
        double gr2 = 0.0;
        const auto q = renderAt (qin, sr, Setup { -1.0, 50.0, 1.0, 4 }, 0, &gr2);
        test::approx (tp::truePeakDbFft (q[0], 16), -20.0, 0.1, "a signal 19 dB below the ceiling passes untouched");
        test::ok (gr2 > -0.05, "and draws no gain reduction");
    }

    // FINDING, recorded here and not fixed (it is outside this task): reset() clears osBuf, the
    // lookahead delays, the sliding max and grDb — but never calls os.reset(), and
    // PolyphaseOversampler holds per-channel up/down ring histories. Measured: after processing a
    // 0.9 FS tone, reset(), then feeding SILENCE, the output is NOT silent — it reaches the ceiling
    // (-1.00 dBFS) and takes 109 samples to decay, where a freshly-prepared instance gives silence. A
    // host calling reset() on transport stop therefore gets a burst of the previous stream. The fix is
    // one line in reset(); it belongs to whoever owns that decision, not to a test.
    //
    // What is asserted is the pair of properties that hold BOTH before and after that fix, so this
    // group keeps its value either way: the residue cannot outlive the filter's own history, and —
    // the property the product actually rests on — it is still bounded by the ceiling.
    test::group ("reset() leaves no residue beyond the filter history, and none above the ceiling");
    {
        const int n = 6000;          // NOT 4096: at that one length the residue happens to land just under
                                     // the ceiling, and the row passed for the wrong reason. See below.
        std::vector<float> loud ((std::size_t) n), silence ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i) loud[(std::size_t) i] = (float) (0.9 * std::sin (2.0 * core::kPi * 1000.0 * (double) i / sr));
        limiter::TruePeakLimiter a, b;
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0;
        (void) a.prepare (sr, n, 1, { 1.0, 4, 32 }); a.setParams (p);
        (void) b.prepare (sr, n, 1, { 1.0, 4, 32 }); b.setParams (p);
        std::vector<float> warm = loud; float* w[1] { warm.data() };
        a.process (w, 1, n);
        a.reset();
        std::vector<float> ya = silence, yb = silence;
        float* pa[1] { ya.data() }; float* pb[1] { yb.data() };
        a.process (pa, 1, n);
        b.process (pb, 1, n);
        int last = -1;
        for (int i = 0; i < n; ++i) if (ya[(std::size_t) i] != yb[(std::size_t) i]) last = i;   // exact, as above
        // The residue cannot outlive the reported latency plus the FIR tail. Written from
        // latencySamples() rather than from a constant 32: the lookahead is part of it, so a fixed
        // number is a hidden 48 kHz assumption — the same construction runs to sample 157 at 96 kHz and
        // 253 at 192 kHz.
        const int historyBaseband = 32;
        // WAS bounded-but-present: reset() cleared the delays, the deque and the gain but never the
        // oversampler's FIR histories, so SILENCE after a reset came out at the ceiling for ~110 samples
        // of the previous stream, peaking 0.036 dB ABOVE it. The old form of this row asserted only that
        // the residue died within the filter history, which was true before and after — it could not
        // close the defect. Equality with a fresh instance can.
        test::ok (last < 0, "silence after reset() is bit-identical to a freshly prepared instance"
                            + std::string (last < 0 ? "" : " (differs to sample " + std::to_string (last) + ")"));
        // And it is bounded — but NOT by the ceiling, which is the sharper half of this finding. reset()
        // leaves the downsampler holding ceiling-level oversampled samples and then feeds it zeros; the
        // decimation FIR reconstructs that step with its own overshoot, so the residue peaks slightly
        // ABOVE the ceiling. Measured +0.036 dB at most buffer lengths (and, by luck, +0.001 at 4096 —
        // which is why the first version of this row asserted "at or below the ceiling" and passed).
        test::ok (tp::samplePeakDb (ya) <= -299.0, "and it is silence, not the previous stream at the ceiling ("
                                                  + dbs (tp::samplePeakDb (ya)) + " dBFS)");
        // The gain state specifically must be cleared, and silence cannot show it. Two details matter and
        // both were found by trying the mutation: the reset has to follow the LOUD pass directly (a
        // silence pass in between lets the release walk grDb back toward 0 on its own), and the level has
        // to be read in a SHORT window right after the reset — measuring a peak over the whole buffer
        // finds the end, where the release has recovered anyway, and the mutation survives.
        // The two defects have to be separated or neither can be measured. A 64-sample silent pass
        // BEFORE the reset flushes the oversampler's 32-sample history (so the residue above cannot
        // drive the gain down and be mistaken for a retained one), while a 500 ms release leaves grDb
        // essentially where the loud pass left it (64 samples move it by 0.3 %). Then reset() is the
        // only thing that can bring the gain back to unity, and the measurement can be taken at once.
        limiter::TruePeakLimiter c;
        limiter::TruePeakLimiterParams cp; cp.ceilingDbTp = -1.0; cp.releaseMs = 500.0;
        (void) c.prepare (sr, n, 1, { 1.0, 4, 32 }); c.setParams (cp);
        std::vector<float> hot = loud; float* ph[1] { hot.data() };
        for (int i = 0; i < n; ++i) hot[(std::size_t) i] = (float) (hot[(std::size_t) i] * 30.0);   // ~29 dB of reduction
        c.process (ph, 1, n);
        test::ok (c.gainReductionDb() < -15.0, "the loud pass really did drive the gain deep ("
                                                + dbs (c.gainReductionDb()) + " dB)");
        std::vector<float> flush (64, 0.0f); float* pf[1] { flush.data() };
        c.process (pf, 1, 64);                                       // flush the filter history, keep grDb
        c.reset();
        std::vector<float> soft ((std::size_t) n);
        for (int i = 0; i < n; ++i) soft[(std::size_t) i] = (float) (0.05 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));
        std::vector<float> ys = soft; float* ps[1] { ys.data() };
        c.process (ps, 1, n);
        const std::vector<float> justAfter (ys.begin() + 2 * historyBaseband, ys.begin() + 2 * historyBaseband + 256);
        test::approx (tp::samplePeakDb (justAfter), -26.02, 0.2,
                      "immediately after reset, quiet material passes at unity — the gain state really is cleared");
    }

    // ---------------------------------------------------------------- what the round trip costs up top
    // Not a ceiling property, but a mastering-relevant one that no test pinned: the limiter's own
    // 0.90 x Nyquist prototype is in the signal path even when nothing is being limited, so the top of
    // the band is attenuated on the way through. Pinned deliberately, because a redesign that widens
    // the pass band should have to change this number on purpose.
    test::group ("Pass-band cost of the round trip (deliberate, pinned, mastering-relevant)");
    {
        struct Row { int p, q; double maxDroopDb; const char* where; };
        const Row rows[] { { 1, 6, 0.02, "fs/6  (8.0 kHz at 48k)" },
                           { 1, 4, 0.02, "fs/4  (12.0 kHz at 48k)" },
                           { 5, 14, 0.10, "5fs/14 (17.1 kHz at 48k)" } };
        for (const auto& r : rows)
        {
            const auto x = tpw::gridTone (sr, 0.10, r.p, r.q, -20.0, 0.0);
            std::vector<std::vector<float>> in { x };
            const auto y = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, 4 }, 0);
            const double droop = tp::truePeakDbFft (x, 16) - tp::truePeakDbFft (y[0], 16);
            test::ok (droop <= r.maxDroopDb, std::string ("round-trip droop at ") + r.where + " is "
                                             + dbs (droop) + " dB (budget " + dbs (r.maxDroopDb) + ")");
        }
    }

    // ...AND THE THREE FREQUENCIES THE HEADER'S OWN FIGURES NAME, which nothing measured until now.
    // That gap is why the header stood for so long saying the ROUND TRIP costs -0.40 / -4.0 / -6.0 dB
    // at 0.40 / 0.44 / 0.45 fs: those are ONE filter pass. The prototype is applied twice — once
    // interpolating, once decimating — so the dB double, and the real figures are the ones pinned here.
    // Two-sided on purpose: an upper budget alone would have accepted the wrong number just as happily.
    test::group ("Round-trip droop where the header quotes it (two-sided, the doc-defect witness)");
    {
        struct Row { int p, q; double expectDb; const char* where; };
        const Row rows[] { { 2,  5, 0.775,  "0.40 fs (19.2 kHz at 48k)" },
                           { 11, 25, 8.065, "0.44 fs (21.1 kHz at 48k)" },
                           { 9, 20, 12.041, "0.45 fs (21.6 kHz at 48k)" } };
        for (const auto& r : rows)
        {
            const auto x = tpw::gridTone (sr, 0.10, r.p, r.q, -20.0, 0.0);
            std::vector<std::vector<float>> in { x };
            const auto y = renderAt (in, sr, Setup { -1.0, 50.0, 1.0, 4 }, 0);
            const double droop = tp::truePeakDbFft (x, 16) - tp::truePeakDbFft (y[0], 16);
            test::approx (droop, r.expectDb, 0.25, std::string ("round-trip droop at ") + r.where);
        }
    }

    // ---------------------------------------------------------------- the hole, pinned as a hole
    // A block larger than maxBlock is REJECTED, and rejection means the audio is returned untouched —
    // so an oversized block is delivered UNLIMITED, silently. The existing suite asserts only that this
    // path does not read out of bounds. It is pinned here as behaviour because the composite that will
    // sit on top of this limiter has to chunk its input, and "it no-ops" and "it ships your peaks" are
    // the same sentence read twice.
    // WAS a pinned defect: a block above maxBlock came back untouched — unlimited, and silently, which
    // for a module whose only promise is a ceiling is the worst possible failure. It is chunked now, and
    // the assertion is EQUIVALENCE rather than "the peak looks limited": a peak-only check would pass a
    // chunker that carried its streaming state wrongly across the boundary.
    test::group ("A block above maxBlock is chunked, and equals the caller having chunked it");
    {
        const int n = 2048;
        std::vector<float> x ((std::size_t) n, 0.0f);
        for (int i = 0; i < n; ++i) x[(std::size_t) i] = (float) (10.0 * std::sin (2.0 * core::kPi * 1000.0 * (double) i / sr));
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = -1.0; p.releaseMs = 50.0;
        for (int mb : { 512, 100 })                      // 100 leaves a remainder chunk, the interesting one
        {
            limiter::TruePeakLimiter one, many;
            test::ok (one.prepare (sr, mb, 1, { 1.0, 4, 32 }) && many.prepare (sr, mb, 1, { 1.0, 4, 32 }),
                      "maxBlock " + std::to_string (mb) + ": both prepared");
            one.setParams (p); many.setParams (p);
            std::vector<float> a2 = x, b2 = x;
            float* pa[1] { a2.data() };
            one.process (pa, 1, n);                                          // one oversized call
            for (int off = 0; off < n; off += mb)
            { float* pb[1] { b2.data() + off }; many.process (pb, 1, std::min (mb, n - off)); }
            test::ok (firstDifference (a2, b2) < 0, "maxBlock " + std::to_string (mb)
                      + ": one call of " + std::to_string (n) + " == the caller chunking it, sample for sample");
            test::ok (tp::samplePeakDb (a2) < -0.9, "maxBlock " + std::to_string (mb)
                      + ": and the whole buffer really is limited (" + dbs (tp::samplePeakDb (a2)) + " dBFS)");
        }
    }

    // ---------------------------------------------------------------- two holes found by the review round
    // Both WERE holes pinned by the previous task and are closed now; the assertions below are the
    // positive form, and each of them fails if its fix is undone.
    test::group ("Non-finite input, and a lookahead change that can no longer be expressed");
    {
        const double C = -1.0;
        // (a) WAS two pinned defects. One +Inf drove grDb to -inf, where -inf * relCoef stays -inf, and
        // the rest of the stream became digital silence; a NaN additionally lodged in the sliding-max
        // deque, which pops on `v[b] <= x` and therefore never pops past a NaN. A gate at the input makes
        // both unreachable — the oversampler copies raw input into its history, so anything downstream of
        // upsample() would already be too late. The assertion is the MAPPING, not merely the absence of
        // NaN: the poisoned buffer must render bit-identically to the buffer with the substitution
        // already applied. A finite 3e38 is in the same table because the clamp, not the isfinite test,
        // is what catches it.
        const int n = 4096;
        for (int kind = 0; kind < 3; ++kind)
        {
            const float poison = kind == 0 ? std::numeric_limits<float>::infinity()
                               : kind == 1 ? std::numeric_limits<float>::quiet_NaN()
                                           : 3.0e38f;
            const float substitute = kind == 2 ? 1.0e6f : 0.0f;
            const char* what = kind == 0 ? "+Inf" : kind == 1 ? "NaN" : "a finite 3e38";
            std::vector<float> bad ((std::size_t) n), good ((std::size_t) n);
            for (int i = 0; i < n; ++i)
            {
                const float v = (float) (4.0 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));   // hot: it must LIMIT
                bad[(std::size_t) i] = good[(std::size_t) i] = v;
            }
            bad[1000] = poison; good[1000] = substitute;
            std::vector<std::vector<float>> bi { bad }, gi { good };
            const auto yb = renderAt (bi, sr, Setup { C, 50.0, 1.0, 4 }, 0);
            const auto yg = renderAt (gi, sr, Setup { C, 50.0, 1.0, 4 }, 0);
            test::ok (firstDifference (yb[0], yg[0]) < 0,
                      std::string ("one ") + what + " renders bit-identically to the value the gate substitutes");
            bool finite = true;
            for (float v : yb[0]) finite &= (bool) std::isfinite (v);
            test::ok (finite, std::string ("one ") + what + " leaves no non-finite sample in the output");
            double tail = 0.0;
            for (int i = 2000; i + 200 < (int) yb[0].size(); ++i) tail = std::max (tail, (double) std::fabs (yb[0][(std::size_t) i]));
            // The substitution level decides how fast the stream comes back, so the two cases get
            // different thresholds rather than one loose one. A gate value of 0 costs nothing; the 1e6
            // clamp is +120 dBFS, so it legitimately ducks and climbs back over a release time — what
            // matters is that it CLIMBS, where +Inf used to leave digital silence (-600 dBFS) for good.
            const double tailDb = 20.0 * std::log10 (tail > 1e-30 ? tail : 1e-30);
            test::ok (tailDb > -100.0, std::string ("the stream after ") + what + " is alive ("
                      + dbs (tailDb) + " dBFS) — +Inf used to mute it permanently");
            if (kind != 2)
                test::ok (tailDb > -3.0, std::string ("and after ") + what
                          + " it is back at the ceiling within the buffer (" + dbs (tailDb) + " dBFS)");
        }
        // A FINITE but absurd ceiling was the same kill through another door: (float)(-1e308 - smaxDb)
        // overflows to -inf. Raised by the review round, not by the task list, and verified here.
        for (double absurd : { -1.0e30, -1.0e308, 1.0e308 })
        {
            std::vector<float> x ((std::size_t) 2048);
            for (int i = 0; i < 2048; ++i) x[(std::size_t) i] = (float) (0.1 * std::sin (0.1 * (double) i));
            limiter::TruePeakLimiter lim; (void) lim.prepare (sr, 2048, 1, { 1.0, 4, 32 });
            limiter::TruePeakLimiterParams pp; pp.ceilingDbTp = absurd; pp.releaseMs = 50.0;
            lim.setParams (pp);
            float* ch[1] { x.data() };
            lim.process (ch, 1, 2048); lim.process (ch, 1, 2048);
            test::ok (std::isfinite (lim.gainReductionDb()), "an absurd finite ceiling (" + dbs (absurd)
                      + " dBTP) leaves the gain state finite, not stuck at -inf");
            test::ok (lim.effectiveCeilingDbTp() >= -200.0 && lim.effectiveCeilingDbTp() <= 60.0,
                      "and the clamp is visible: effectiveCeilingDbTp() reports " + dbs (lim.effectiveCeilingDbTp()));
        }

        // (b) Raising lookaheadMs mid-stream used to RE-EMIT audio that had already gone out, at a gain
        // computed for the quiet material after it: +28.8 dB over the ceiling. setDelay moves only a read
        // offset, and the widened window cannot recover detector entries the narrow one evicted. The
        // parameter now lives in the prepare-time config, so the call that did it does not COMPILE — the
        // assertion is therefore about the contract rather than a number: latency is fixed for the life
        // of a prepared stream, and per-block automation cannot move it.
        {
            limiter::TruePeakLimiter a, b2;
            (void) a.prepare (sr, 64, 1, { 0.5, 4, 32 });
            (void) b2.prepare (sr, 64, 1, { 5.0, 4, 32 });
            test::ok (a.latencySamples() != b2.latencySamples(),
                      "lookahead still changes latency (" + std::to_string (a.latencySamples()) + " vs "
                      + std::to_string (b2.latencySamples()) + ") — which is why it cannot be a per-block parameter");
            limiter::TruePeakLimiterParams pp; pp.ceilingDbTp = C; pp.releaseMs = 1.0;
            std::vector<float> burst ((std::size_t) 8192, 0.0f);
            for (int i = 470; i < 500; ++i)
                burst[(std::size_t) i] = (float) (std::pow (10.0, 30.0 / 20.0) * std::sin (2.0 * core::kPi * 0.25 * (double) i + 0.7));
            const int before = a.latencySamples();
            float* ch[1];
            for (int i = 0; i < 8192; i += 64)
            {
                pp.releaseMs = (i == 640) ? 2.0 : 1.0;                    // per-block automation, every block
                a.setParams (pp);
                ch[0] = burst.data() + i; a.process (ch, 1, 64);
            }
            test::ok (a.latencySamples() == before, "and per-block setParams cannot move it");
            test::ok (tp::samplePeakDb (burst) - C <= tpw::deliveredBudgetDb (4),
                      "the burst is delivered " + dbs (tp::samplePeakDb (burst) - C)
                      + " dB over the ceiling, inside the budget — it used to be +28.4");
        }
    }

    // ---------------------------------------------------------------- what the review round found
    // Every one of these closes a defect the first version of this task shipped past. They are grouped
    // because they share a shape: state that is PER CHANNEL meeting state that is SHARED, and validation
    // that checked one end of a range.
    test::group ("Channel count, validation ranges, and a release long enough to stop being one");
    {
        const double C = -1.0;
        // A channel that sits out a few blocks used to come back re-emitting audio whose detector
        // entries had expired globally while it was away: +19.76 dB over the ceiling, at every factor.
        for (int F : { 2, 4, 8 })
        {
            limiter::TruePeakLimiter lim;
            (void) lim.prepare (sr, 4096, 2, { 1.0, F, 32 });
            limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = 50.0;
            lim.setParams (p);
            std::vector<float> L (40, 0.0f), R (40, 0.0f);
            R[0] = 10.0f;                                                 // +20 dBFS, RIGHT channel only
            { float* io[2] { L.data(), R.data() }; lim.process (io, 2, 40); }
            std::vector<float> mono (24000, 0.0f);
            { float* io[1] { mono.data() }; lim.process (io, 1, 24000); } // the right channel freezes here
            std::vector<float> L2 (4096, 0.0f), R2 (4096, 0.0f);
            { float* io[2] { L2.data(), R2.data() }; lim.process (io, 2, 4096); }
            test::ok (tp::samplePeakDb (R2) <= C, "F=" + std::to_string (F)
                      + ": stereo -> mono -> stereo leaves nothing above the ceiling on the resumed channel ("
                      + dbs (tp::samplePeakDb (R2)) + " dBFS; it used to be +18.76)");
        }
        // Validation covered one end of two ranges only.
        {
            limiter::TruePeakLimiter lim;
            test::ok (! lim.prepare (sr, 64, 1, { 1.0, 4, std::numeric_limits<int>::max() }),
                      "tapsPerPhase = INT_MAX is refused (it overflowed factor*taps inside the oversampler)");
            test::ok (! lim.prepare (1.0, 64, 1, { 1.0, 4, 32 }),
                      "a 1 Hz rate is refused (20 ms could not hold the minimum lookahead, and the clamp "
                      "became std::clamp(x, 2, 1) — lo > hi is undefined)");
            test::ok (lim.prepare (8000.0, 64, 1, { 1.0, 4, 32 }), "but a legitimately low rate still prepares");
            test::ok (lim.isPrepared() && lim.oversampleFactor() == 4,
                      "and the getters report the live topology (isPrepared, oversampleFactor "
                      + std::to_string (lim.oversampleFactor()) + ")");
        }
        // A release long enough that exp(-1/t) rounds to 1.0f stopped being a release at all.
        for (double rel : { 1000.0, 175000.0, 1.0e9 })
        {
            limiter::TruePeakLimiter lim; (void) lim.prepare (sr, 512, 1, { 1.0, 4, 32 });
            limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = rel;
            lim.setParams (p);
            std::vector<float> hot (512, 4.0f); float* io[1] { hot.data() };
            lim.process (io, 1, 512);
            const double engaged = lim.gainReductionDb();
            std::vector<float> quiet (480000, 0.0f); float* q[1] { quiet.data() };
            lim.process (q, 1, 480000);                                   // ten seconds of silence
            test::ok (engaged < -1.0, "release " + dbs (rel) + " ms: the limiter engaged (" + dbs (engaged) + " dB)");
            test::ok (lim.gainReductionDb() > engaged, "release " + dbs (rel)
                      + " ms: and the gain still recovers, however slowly (" + dbs (engaged) + " -> "
                      + dbs (lim.gainReductionDb()) + " dB) — at 175 s it used to be stuck for good");
        }
        // The gate was only ever tested with POSITIVE poison.
        for (int kind = 0; kind < 2; ++kind)
        {
            const float poison = kind == 0 ? -std::numeric_limits<float>::infinity() : -3.0e38f;
            const float substitute = kind == 0 ? 0.0f : -1.0e6f;
            std::vector<float> bad (4096), good (4096);
            for (int i = 0; i < 4096; ++i)
            {
                const float v = (float) (4.0 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));
                bad[(std::size_t) i] = good[(std::size_t) i] = v;
            }
            bad[1000] = poison; good[1000] = substitute;
            std::vector<std::vector<float>> bi { bad }, gi { good };
            test::ok (firstDifference (renderAt (bi, sr, Setup { C, 50.0, 1.0, 4 }, 0)[0],
                                       renderAt (gi, sr, Setup { C, 50.0, 1.0, 4 }, 0)[0]) < 0,
                      std::string (kind == 0 ? "-Inf" : "a finite -3e38")
                      + " maps to its substitute too — the gate is not one-sided");
        }
        // Chunking and reset were only ever nulled in mono, so a per-channel offset or a per-channel
        // reset could have been wrong in the second channel and nothing would have said so.
        {
            const int n = 2048;
            std::vector<float> a1 (n), a2 (n), b1 (n), b2 (n);
            for (int i = 0; i < n; ++i)
            {
                a1[(std::size_t) i] = a2[(std::size_t) i] = (float) (10.0 * std::sin (2.0 * core::kPi * 1000.0 * (double) i / sr));
                b1[(std::size_t) i] = b2[(std::size_t) i] = (float) (7.0 * std::sin (2.0 * core::kPi * 1700.0 * (double) i / sr + 0.9));
            }
            limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = 50.0;
            limiter::TruePeakLimiter one, many;
            (void) one.prepare (sr, 512, 2, { 1.0, 4, 32 }); (void) many.prepare (sr, 512, 2, { 1.0, 4, 32 });
            one.setParams (p); many.setParams (p);
            { float* io[2] { a1.data(), b1.data() }; one.process (io, 2, n); }
            for (int off = 0; off < n; off += 512)
            { float* io[2] { a2.data() + off, b2.data() + off }; many.process (io, 2, 512); }
            test::ok (firstDifference (a1, a2) < 0 && firstDifference (b1, b2) < 0,
                      "STEREO: one oversized call == the caller chunking it, on both channels");
            limiter::TruePeakLimiter warm, fresh;
            (void) warm.prepare (sr, n, 2, { 1.0, 4, 32 }); (void) fresh.prepare (sr, n, 2, { 1.0, 4, 32 });
            warm.setParams (p); fresh.setParams (p);
            std::vector<float> w1 = a1, w2 = b1;
            { float* io[2] { w1.data(), w2.data() }; warm.process (io, 2, n); }
            warm.reset();
            std::vector<float> x1 (n, 0.0f), x2 (n, 0.0f), y1 (n, 0.0f), y2 (n, 0.0f);
            { float* io[2] { x1.data(), x2.data() }; warm.process (io, 2, n); }
            { float* io[2] { y1.data(), y2.data() }; fresh.process (io, 2, n); }
            test::ok (firstDifference (x1, y1) < 0 && firstDifference (x2, y2) < 0,
                      "STEREO: reset() equals a fresh instance on both channels");
        }
    }

    // ---------------------------------------------------------------- the contract surface itself
    // The review round mutated the header line by line and found that most of what lives OUTSIDE the
    // sample loop was unpinned: validation ranges, fallbacks, getters, and the state a failed re-prepare
    // leaves behind. The reference NULL cannot see any of it — the transcription has no guards at all —
    // so it needs its own rows.
    test::group ("Validation ranges, fallbacks and getters are pinned, not just the sample loop");
    {
        const double C = -1.0;
        limiter::TruePeakLimiter lim;
        test::ok (! lim.prepare (1.0e300, 512, 1),                "an absurd sample rate is refused");
        test::ok (! lim.prepare (sr, 512, 1, { -1.0, 4, 32 }),    "a negative lookahead is refused");
        test::ok (! lim.prepare (sr, 512, 1, { 1.0, 32, 32 }),    "a factor above the cap is REFUSED, not silently reduced");
        test::ok (! lim.prepare (sr, 512, core::kMaxChannels + 1),"more channels than supported is REFUSED — clamping would have "
                                                                 "let the surplus pass through unlimited");
        test::ok (lim.prepare (sr, 512, 2, { 1.0, 8, 32 }), "a valid configuration prepares");
        test::ok (lim.isPrepared() && lim.oversampleFactor() == 8 && lim.lookaheadSamples() == 48,
                  "and the getters report the live topology (factor " + std::to_string (lim.oversampleFactor())
                  + ", lookahead " + std::to_string (lim.lookaheadSamples()) + " samples)");
        const int wasLatency = lim.latencySamples();
        test::ok (! lim.prepare (sr, 512, 1, { 1.0, 4, 2 }), "then a bad re-prepare fails");
        test::ok (! lim.isPrepared() && lim.latencySamples() == 0 && lim.lookaheadSamples() == 0
                  && lim.oversampleFactor() == 0,
                  "and every accessor reports nothing-prepared rather than the previous topology (was latency "
                  + std::to_string (wasLatency) + ")");

        // maxBlock is a scratch-buffer hint now, so a big one must PREPARE rather than fail — failing
        // would leave process() returning the buffer untouched, the very defect this task closed.
        limiter::TruePeakLimiter big;
        test::ok (big.prepare (sr, 1 << 23, 1, { 1.0, 4, 32 }), "an oversized maxBlock is clamped, not refused");
        {
            limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = 50.0; big.setParams (p);
            std::vector<float> x (4096, 4.0f); float* io[1] { x.data() };
            big.process (io, 1, 4096);
            test::ok (tp::samplePeakDb (x) < -0.9, "and it still limits (" + dbs (tp::samplePeakDb (x)) + " dBFS)");
        }

        // Non-finite parameters fall back to the documented DEFAULTS, which "the output stayed finite"
        // does not pin — the fallback value could be anything and that assertion would still pass.
        limiter::TruePeakLimiter fb;
        test::ok (fb.prepare (sr, 512, 1, { 1.0, 4, 32 }), "prepared for the fallback checks");
        limiter::TruePeakLimiterParams np;
        np.ceilingDbTp = std::numeric_limits<double>::quiet_NaN();
        np.releaseMs   = std::numeric_limits<double>::quiet_NaN();
        fb.setParams (np);
        test::approx (fb.effectiveCeilingDbTp(), -1.0, 1.0e-12, "a NaN ceiling falls back to the -1 dBTP default");
        test::approx (fb.effectiveReleaseMs(), 100.0, 0.05, "a NaN release falls back to the 100 ms default"); // 0.05: the float coefficient round trip
        limiter::TruePeakLimiterParams cp; cp.ceilingDbTp = -1.0e9; cp.releaseMs = 50.0;
        fb.setParams (cp);
        test::approx (fb.effectiveCeilingDbTp(), -200.0, 1.0e-12, "and an absurd finite ceiling lands on the documented floor");
    }

    test::group ("Poison and reset, on the channel and at the moment the mono rows could not reach");
    {
        const double C = -1.0;
        // The gate loops over channels; a loop that only ever read channel 0 passed every mono row.
        std::vector<float> l (4096), r0 (4096), r1 (4096);
        for (int i = 0; i < 4096; ++i)
        {
            const float v = (float) (4.0 * std::sin (2.0 * core::kPi * 997.0 * (double) i / sr));
            l[(std::size_t) i] = v; r0[(std::size_t) i] = r1[(std::size_t) i] = v * 0.5f;
        }
        r0[1000] = std::numeric_limits<float>::quiet_NaN();
        r1[1000] = 0.0f;
        std::vector<std::vector<float>> bi { l, r0 }, gi { l, r1 };
        const auto yb = renderAt (bi, sr, Setup { C, 50.0, 1.0, 4 }, 0);
        const auto yg = renderAt (gi, sr, Setup { C, 50.0, 1.0, 4 }, 0);
        test::ok (firstDifference (yb[1], yg[1]) < 0 && firstDifference (yb[0], yg[0]) < 0,
                  "a NaN on channel 1 is gated exactly like one on channel 0 — both channels stay identical "
                  "to the substituted render");

        // reset() must clear the DETECTOR, not only the gain: without slide.reset() the first sample after
        // a reset still saw the old window and pulled the gain 27 dB down.
        limiter::TruePeakLimiter lim;
        test::ok (lim.prepare (sr, 4096, 1, { 1.0, 4, 32 }), "prepared for the reset checks");
        limiter::TruePeakLimiterParams p; p.ceilingDbTp = C; p.releaseMs = 50.0; lim.setParams (p);
        std::vector<float> hot (4096, 4.0f); float* ph[1] { hot.data() };
        lim.process (ph, 1, 4096);
        lim.reset();
        std::vector<float> one (1, 0.0f); float* po[1] { one.data() };
        lim.process (po, 1, 1);
        test::approx (lim.gainReductionDb(), 0.0, 0.0,
                      "one silent sample after reset() draws exactly no gain reduction (the detector was cleared too)");

        // And the gain state must reach EXACT zero after limiting, not a denormal tail — which is the one
        // observable difference the poison/denormal flush makes.
        lim.process (ph, 1, 4096);
        std::vector<float> quiet (120000, 0.0f); float* pq[1] { quiet.data() };
        lim.process (pq, 1, 120000);                                  // 2.5 s of silence
        test::approx (lim.gainReductionDb(), 0.0, 0.0,
                      "and after limiting plus 2.5 s of silence the gain state is exactly 0, not a denormal");
    }

    // ---------------------------------------------------------------- the corner the floors do NOT cover
    // Each floor was measured with the OTHER parameter at a musical value. Both at once is a different
    // point, and it sits OUTSIDE the envelope — bringing it inside would need floors at 24 baseband
    // samples (0.5 ms), which is a musical setting and would change the sound of a legitimate one. So it
    // is characterised here with its number instead of being clamped away, and the header says so.
    test::group ("Both floors at once: outside the envelope, on purpose, with the number pinned");
    {
        const double C = -1.0;
        for (int F : { 2, 4, 8 })
        {
            double worst = -1e9;
            std::vector<std::vector<float>> ws { tpw::clickTrain (sr, 0.12, 3.0, 26.0, 0.5),
                                                 tpw::denseNoise (sr, 0.12, 10.0),
                                                 tpw::denseNoise (sr, 0.12, 20.0) };
            for (auto& w : ws)
            {
                std::vector<std::vector<float>> in { w };
                worst = std::max (worst, tp::truePeakDbFft (renderAt (in, sr, Setup { C, 0.0, 0.0, F }, 0)[0]) - C);
            }
            test::ok (worst <= tpw::deliveredBudgetDb (F) + 0.5,
                      "F=" + std::to_string (F) + ": at BOTH floors the worst witness is " + dbs (worst)
                      + " dB over — outside the " + dbs (tpw::deliveredBudgetDb (F))
                      + " dB envelope, and bounded well inside the unfloored +2.8");
            if (F == 8)
                test::ok (worst > tpw::deliveredBudgetDb (F),
                          "and at 8x it really is outside it (" + dbs (worst) + " > " + dbs (tpw::deliveredBudgetDb (F))
                          + ") — this row exists so that fact cannot be forgotten");
        }
    }

    // ---------------------------------------------------------------- 44.1 kHz, and the rates above it
    test::group ("Rate independence at the two rates the conformance claim covers (44.1 / 48 kHz)");
    for (double rate : { 44100.0, 48000.0 })
    {
        const double C = -1.0;
        const int F = 4;
        const int M = tpw::gridPhaseCount (1, 3, F);
        double worst = -1e9;
        for (int k = 0; k < 2 * M; ++k)
        {
            const auto x = tpw::gridTone (rate, 0.05, 1, 3, 12.0, (double) k / (double) (2 * M));
            std::vector<std::vector<float>> in { x };
            const auto y = renderAt (in, rate, Setup { C, 50.0, 1.0, F }, 0);
            worst = std::max (worst, tp::truePeakDbFft (y[0]));
        }
        test::approx (worst - C, tpw::gridBreachDb (1, 3, F), 0.015,
                      "fs/3 excess at " + std::to_string ((int) rate) + " Hz matches the closed form (rate-independent)");
    }

    return test::report();
}
