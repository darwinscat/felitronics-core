// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Config.h>
#include <felitronics/core/Math.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace felitronics::analysis
{

//==============================================================================
// felitronics::analysis::ClipDetector — finds SAMPLE CLIPPING in delivered audio: the runs of samples whose peaks hit
// a ceiling and became FLAT, with their positions per channel, plus the sample peak and DC offset of each channel.
//
// WHAT CLIPPING IS, AND THE TWO THINGS THAT FOLLOW. A clamp replaces the top of every excursion past its ceiling with
// the ceiling itself, so the file carries runs of equal samples where a smooth signal would have kept moving. Hence:
//   · FLATNESS, NOT LOUDNESS. Consecutive loud samples are not evidence — a full-scale sine and a master normalised to
//     0 dBFS have them and are clean. Nothing here compares a sample with full scale.
//   · LEVEL-FREE. Clipped and then turned down is still clipped: the flat tops moved down. The ceiling is found where
//     the flat tops are, at whatever level, per channel and per polarity.
//
// THE RULE. Per channel, the stream is cut into BANDS — maximal runs of finite samples whose range is <= tau = 2q.
// q is the file's quantum as far as the stream has shown it: the smaller of the smallest non-zero step between adjacent
// samples and the coarsest 2^-k grid (k <= 24) every sample lies on — integer PCM sits on its grid exactly. Neither can
// fall below the true quantum, so tau never under-estimates the rounding a flat top carries (2q admits TPDF dither of
// ±1 LSB), and both only shrink as the stream goes on, so an early decision is the cautious one. A band of
// L >= 2 samples whose known neighbours both lie below it is a TOP (both above: a FLOOR); a side that borders a
// non-finite sample or the stream edge is unknown and is not used. A top/floor that no finite sample within ±20 ms
// passes by more than tau is a candidate, and it is a clipped run on either of two kinds of evidence:
//
//   RAMP — the flat top itself proves the clamp, so both its sides must be known (a side at a hole or at the stream edge
//     proves nothing: a signal cut into trailing digital silence is not a clamp) and all of these must hold:
//     (1) THE SMOOTH-CREST BOUND. A parabolic crest that stays inside a band of width tau for L samples cannot ENTER
//         it faster than V(L) = (tau + q)·((L+1)/(L−1))² + q: held within tau over L samples its curvature is
//         <= 8·tau/(L−1)², so the sample outside sits at most that curvature × (L+1)²/8 lower, rounding (q) included on
//         both ends. The bound is per side and the centred crest is its worst case, so one fast side would already be
//         inconsistent with a parabola; both are asked for, as margin against crests that are not parabolas (asking
//         one side only finds 5 more grid-locked pure tones of 426 clamped synthetic signals and, at the final bounds,
//         no false run — an option, not taken).
//     (2) THE CHANCE BOUND. Noise-like material (hats, snare, dither) is not a smooth crest: its samples can land inside
//         a band by coincidence. With flank activity sigma — the mean |step| over up to 8 steps before the entry and 8
//         after the exit, SKIPPING every step that touches this band's own level, because a step onto or along the
//         ceiling is the clamp, not the signal's activity — the chance that L samples fall inside the band's actual
//         range r is ~ ((r + q)/sigma)^(L−1); it must be <= 1e-6 (evaluated by multiplication, not logarithms, so no libm
//         decides a verdict and every tier agrees). The bound is set on REAL recordings, where it is
//         tight: over 72 clean files made from 18 full-length tracks and stems (16 bit at 0, −12 and −24 dB from the
//         peak, 24 bit at −6 dB) the false runs number 434 / 135 / 19 / 1 / 0 / 0 at 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 /
//         1e-7, while 52 minutes of synthetic programme were clean already at 1e-3. 1e-6 is one decade past the last
//         false run (one run, in a stem whose history is not known); 1e-7 would buy a second decade for 0.2 % of the
//         runs of length 2, 5 % of a dithered clamp's runs, and three of the six files driven less than 1 dB over the
//         ceiling that 1e-6 finds. The skip is what keeps a click in digital
//         silence and a one-sample pulse on a flat baseline out: counted, their single spike made a flat stretch look
//         like a clamp; skipped, they have no activity at all.
//     (3) NOT A STEP FROM REST. A square wave, a PWM or a sample-and-hold staircase arrives at its flat top by a jump
//         from ANOTHER flat top; a clamp is approached by a signal that is moving. A flat pair among the two samples
//         before the entry (after the exit) at a level other than this band's means a step, and a step is no ramp.
//         This is the whole answer to "why is a square wave not clipped": by itself its flat tops are never evidence.
//         Its limit is honest, too — a clamp driven so hard that it swings from one ceiling to the other in <= 2
//         samples IS that square wave, and no local feature can separate the two (one PCM, two histories).
//     (4) A BAND THAT IS NOT EXACTLY FLAT MUST BE MET BY A LINE, NOT BY A TURN. Bound (1) is for parabolas, and real
//         crests are often flatter — bass through an amp, a compressor, a saturator — and enter a 1–2-quantum band
//         faster than V(L). What separates them is the approach: a crest TURNS at its top, so going outward its steps
//         grow — by 2/(2n+1) for a parabola whose apex is n samples in, by more for a flatter one; a clamp is reached by
//         a signal still moving, whose steps barely change. On each side, over the three samples outside,
//         (step2 − step1)/step1 must be <= 0.4. Measured on the real recordings above: the false runs that passed
//         (1)–(3) turned by 0.58–1.38, the ramps of clamps by −0.15 at the median and <= 0.48 for 90 % of them; at the
//         final chance bound this test is what keeps six runs in one clean bass stem out (at 1e-5 it removed 18 of 19),
//         and every threshold from 0.3 to 0.75 gives the same result. An exactly flat band (every sample equal) is what
//         an undithered clamp leaves and is exempt: shallow clamps turn like crests do, and asking it of them lost 3 of
//         84 clamped files. A dithered clamp is rough; the ramps it fails are recovered by CEILING evidence.
//   CEILING — once a RAMP run has fixed a channel's ceiling of that polarity, any candidate of L >= 2 at that level
//     (within tau) is clipped too: a clamp repeats its level exactly. This is what finds the short runs, the step-shaped
//     runs of hard-driven material and the halves of a run a non-finite sample split, and it cannot fire in a file that
//     never produced a ramp. It looks forward only: the few runs of a polarity that come before its first ramp are missed.
//
// WHAT IT DOES NOT SEE, AND SAYS SO. A single clamped sample has no flat top (a run of 1 is found only by accident).
// A clamp that a lossy codec, a resampler, a post-clip filter or added noise has since reshaped is no longer flat, and
// flatness is the evidence; the fraction lost is measured and published with the module (CHANGELOG), not hidden: AAC,
// MP3 and resampling after a clamp leave 0 % of its runs findable, noise after a 16-bit clamp 96 % at ±1 LSB, 50 % at
// ±2 and none from ±4. A pure tone locked to the sample grid can hide its clamp: the stream never shows it a slow step,
// so q stays loose (music has slow places). Soft saturation
// is not a clamp. True peak is a different quantity (analysis::TruePeakMeter): a file can be over 0 dBTP with no
// clipped sample and clipped with every sample under −6 dBFS. And measurement::PeakClip — a capture-gating probe with
// its own threshold — is unrelated.
//
// FORM. prepare / process / reset / finish. process() is READ-ONLY and allocates nothing; its result is a function of
// the samples alone, bit-identical under any slicing of the stream into calls (a per-sample state machine, frame by
// frame). A decision is made 20 ms of audio after a run ends (the neighbourhood it must not exceed), so finish() is
// called once at the end of the stream to decide the last runs; after it the report is final and process() refuses
// until reset(). Non-finite samples are holes: they break a run, count in nonFiniteSamples(), and are never part of a
// peak, a DC or a decision. A channel not fed by a call that carries samples (index >= nch, law 11a) is a hole for
// those samples too. The run list holds up to ClipDetectorParams::maxRuns entries (set before prepare()); every count
// keeps counting past it and runsComplete() says whether the list is whole.
struct ClipDetectorParams
{
    int maxRuns = 1 << 16;          // capacity of the run list; takes effect at the next prepare()
};

enum class ClipEvidence : std::uint8_t
{
    Ramp    = 1,                    // the run itself passed the smooth-crest, chance and step tests
    Ceiling = 2                     // the run sits at a ceiling a Ramp run of this channel and polarity already fixed
};

struct ClipRun
{
    std::int64_t start  = 0;        // first sample of the run, in samples since reset()
    std::int64_t length = 0;        // samples
    double       level  = 0.0;      // mean value of the run's samples — the ceiling as it sits in the file
    int          channel = 0;
    int          sign   = 0;        // +1 a flat top, −1 a flat floor
    ClipEvidence evidence = ClipEvidence::Ramp;
};

class ClipDetector
{
public:
    static constexpr double kWindowMs        = 20.0;    // the neighbourhood a flat top may not be passed in
    static constexpr double kBandQuanta      = 2.0;     // tau = 2q
    static constexpr int    kFlankSteps      = 8;       // sigma's reach on each side
    static constexpr double kChance          = 1.0e-6;  // the chance bound
    static constexpr double kRoughTurnRatio  = 0.4;     // a band that is not exactly flat: outward steps may grow by <= 40 %
    static constexpr int    kMinRampLength   = 3;
    static constexpr int    kMinCeilingLength = 2;
    static constexpr double kMinSampleRate   = 1000.0;
    static constexpr double kMaxSampleRate   = 768000.0;
    static constexpr int    kMaxRunsLimit    = 1 << 24;

    // WHAT prepare() ALLOCATES, from the function prepare() sizes and validates itself with (law 11d). ok == false on
    // exactly the arguments prepare() refuses, and then every count is zero.
    struct Storage
    {
        bool        ok = false;
        std::int64_t windowSamples = 0;
        std::size_t dequeEntries   = 0;     // two monotone deques per channel
        std::size_t pendingEntries = 0;     // candidates waiting for their neighbourhood, all channels
        std::size_t runEntries     = 0;
        std::size_t channels       = 0;
        std::uint64_t bytes() const noexcept
        {
            return (std::uint64_t) dequeEntries * sizeof (WindowEntry) + (std::uint64_t) pendingEntries * sizeof (Pending)
                 + (std::uint64_t) runEntries * sizeof (ClipRun) + (std::uint64_t) channels * sizeof (Channel);
        }
    };
    static Storage storageFor (double sampleRate, int maxChannels, int maxRuns) noexcept
    {
        Storage st;
        if (! (sampleRate >= kMinSampleRate && sampleRate <= kMaxSampleRate)) return st;   // NaN fails too
        if (maxChannels < 1 || maxChannels > core::kMaxChannels) return st;
        if (maxRuns < 0 || maxRuns > kMaxRunsLimit) return st;
        st.ok = true;
        st.windowSamples = windowFor (sampleRate);
        st.channels = (std::size_t) maxChannels;
        st.dequeEntries = st.channels * 2u * dequeCapacityFor (st.windowSamples);
        st.pendingEntries = st.channels * pendingCapacityFor (st.windowSamples);
        st.runEntries = (std::size_t) maxRuns;
        return st;
    }

    void setParams (const ClipDetectorParams& p) noexcept { params_ = p; }   // maxRuns: at the next prepare()

    [[nodiscard]] bool prepare (double sampleRate, int /*maxBlock: nothing is sized by it*/, int maxChannels) noexcept
    {
        prepared_ = false;                                               // law 11b: disarm, validate, write
        const Storage st = storageFor (sampleRate, maxChannels, params_.maxRuns);
        if (! st.ok) return false;
        sampleRate_ = sampleRate;
        channels_ = maxChannels;
        W_ = st.windowSamples;
        dequeCap_ = dequeCapacityFor (W_);
        pendingCap_ = pendingCapacityFor (W_);
        deques_.assign (st.dequeEntries, WindowEntry {});
        pending_.assign (st.pendingEntries, Pending {});
        runs_.assign (st.runEntries, ClipRun {});
        chans_.assign (st.channels, Channel {});
        prepared_ = true;
        reset();
        return true;
    }

    void reset() noexcept
    {
        for (auto& ch : chans_) ch = Channel {};
        t_ = 0;
        runCount_ = 0;
        finished_ = false;
    }

    static constexpr int latencySamples() noexcept { return 0; }          // a read-only sink
    std::int64_t decisionDelaySamples() const noexcept { return W_; }      // how far behind the stream a run is decided

    // READ-ONLY. Law 11: malformed → unprepared → finished → nch > maxChannels → n == 0 → run (nch == 0: all holes).
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (! prepared_ || finished_) return false;
        if (numChannels > channels_) return false;
        if (n == 0) return true;
        for (int i = 0; i < n; ++i)
        {
            for (int c = 0; c < channels_; ++c)
            {
                if (c < numChannels) step (c, in[c][i], true);
                else                 step (c, 0.0f, false);
            }
            ++t_;
        }
        return true;
    }

    // End of stream: the last runs are decided with the neighbourhood the stream had. Idempotent.
    void finish() noexcept
    {
        if (! prepared_ || finished_) return;
        for (int c = 0; c < channels_; ++c)
        {
            Channel& ch = chans_[(std::size_t) c];
            if (ch.inBand) closeBand (c, t_);
            for (int k = 0; k < ch.pendingCount; ++k) pendingAt (c, k).finalQ = ch.q;   // no e+2 sample will come
        }
        for (;;)                                                          // decide in (end, channel) order, as the stream would
        {
            int best = -1; std::int64_t bestEnd = 0;
            for (int c = 0; c < channels_; ++c)
            {
                const Channel& ch = chans_[(std::size_t) c];
                if (ch.pendingCount == 0) continue;
                const std::int64_t e = pendingAt (c, 0).end;
                if (best < 0 || e < bestEnd) { best = c; bestEnd = e; }
            }
            if (best < 0) break;
            decideFront (best);
        }
        finished_ = true;
    }

    // --- the report (partial while streaming, final after finish()) ---
    // Every accessor is total: an index outside the report (a run past storedRunCount(), a channel outside the prepared
    // width, anything before prepare()) answers an empty run or zero rather than reading past a vector (P58 crew round).
    bool isFinished() const noexcept { return finished_; }
    bool clipped() const noexcept { return runCount_ > 0; }
    std::int64_t runCount() const noexcept { return runCount_; }
    std::int64_t storedRunCount() const noexcept { return std::min<std::int64_t> (runCount_, (std::int64_t) runs_.size()); }
    bool runsComplete() const noexcept { return runCount_ <= (std::int64_t) runs_.size(); }
    ClipRun run (std::int64_t i) const noexcept { return i >= 0 && i < storedRunCount() ? runs_[(std::size_t) i] : ClipRun {}; }
    std::int64_t samplesProcessed() const noexcept { return t_; }
    int channels() const noexcept { return channels_; }

    double samplePeak (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].peak : 0.0; }
    double samplePeakDb (int c) const noexcept { return toDb (samplePeak (c)); }
    double samplePeak() const noexcept { double m = 0; for (const auto& ch : chans_) m = std::max (m, ch.peak); return m; }
    double samplePeakDb() const noexcept { return toDb (samplePeak()); }
    double dcOffset (int c) const noexcept                                 // mean of the finite samples
    {
        if (! has (c)) return 0.0;
        const Channel& ch = chans_[(std::size_t) c];
        return ch.finite > 0 ? (ch.dcSum + ch.dcComp) / (double) ch.finite : 0.0;
    }
    std::int64_t finiteSamples (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].finite : 0; }
    std::int64_t nonFiniteSamples (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].nonFinite : 0; }
    std::int64_t runCount (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].runs : 0; }
    std::int64_t clippedSamples (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].clippedSamples : 0; }
    std::int64_t longestRun (int c) const noexcept { return has (c) ? chans_[(std::size_t) c].longest : 0; }

private:
    static constexpr int kPre = kFlankSteps + 1;        // samples kept before a band: y[k-1] … y[k-9]
    static constexpr int kPost = kFlankSteps + 1;       // samples collected after it: y[e] … y[e+8]

    struct WindowEntry { std::int64_t t = 0; float v = 0.0f; };

    struct Pending
    {
        std::int64_t start = 0, end = 0;                 // band [start, end)
        double lo = 0, hi = 0, sum = 0;
        double preMax = 0, preMin = 0;                   // over the finite samples of [start - W, start)
        double q = 0;                                    // q as of sample end + 2
        double finalQ = 0;                               // q at finish(), for a run whose end + 2 never arrived
        bool qTaken = false;
        int postCount = 0;
        float pre[kPre] {};
        float post[kPost] {};
    };

    struct Channel
    {
        double q = std::numeric_limits<double>::infinity();        // min (minStep, lattice)
        double minStep = std::numeric_limits<double>::infinity();
        int latticeK = 0;                                          // the coarsest 2^-k grid so far; > 24 = none
        double latticeQ = 1.0;                                     // 2^-latticeK, or infinity past 24
        float last = std::numeric_limits<float>::quiet_NaN();
        float hist[kPre] {};                             // hist[0] newest; NaN = hole or before the stream
        bool histInit = false;
        bool inBand = false;
        std::int64_t bStart = 0;
        double bLo = 0, bHi = 0, bSum = 0, bPreMax = 0, bPreMin = 0;
        float bPre[kPre] {};
        int dqMaxHead = 0, dqMaxCount = 0, dqMinHead = 0, dqMinCount = 0;
        int pendingHead = 0, pendingCount = 0;
        double ceiling[2] {};                            // [top, floor]
        bool ceilingSet[2] { false, false };
        double peak = 0, dcSum = 0, dcComp = 0;
        std::int64_t finite = 0, nonFinite = 0, runs = 0, clippedSamples = 0, longest = 0;
    };

    static std::int64_t windowFor (double sampleRate) noexcept { return (std::int64_t) (sampleRate * kWindowMs / 1000.0); }
    static std::size_t dequeCapacityFor (std::int64_t W) noexcept { return (std::size_t) W + 2u; }
    // Bands tile the stream and a candidate is >= 2 samples long, so at most W/2 + 1 of them end inside one window.
    static std::size_t pendingCapacityFor (std::int64_t W) noexcept { return (std::size_t) W / 2u + 2u; }
    // ratio^n <= kChance, by plain multiplication rather than log10/pow: no libm, so every tier reaches the same verdict.
    // The product only falls (ratio < 1) and stops at the bound, so it neither underflows nor runs past the first answer.
    static bool chanceBelow (double ratio, std::int64_t n) noexcept
    {
        if (! (ratio < 1.0)) return false;
        double prod = 1.0;
        for (std::int64_t i = 0; i < n; ++i)
            if ((prod *= ratio) <= kChance) return true;
        return false;
    }
    static constexpr int kMaxLatticeK = 24;
    // the smallest k >= 0 with x·2^k an integer (a finite non-zero float); a subnormal answers more than kMaxLatticeK
    static int gridExponent (float x) noexcept
    {
        const std::uint32_t bits = std::bit_cast<std::uint32_t> (x);
        const int ex = (int) ((bits >> 23) & 0xFFu);
        if (ex == 0) return kMaxLatticeK + 1;
        const std::uint32_t mant = (bits & 0x7FFFFFu) | 0x800000u;
        return std::max (0, 150 - ex - std::countr_zero (mant));
    }
    static double toDb (double lin) noexcept { return lin > 1.0e-10 ? core::gainToDb (lin) : -200.0; }
    bool has (int c) const noexcept { return c >= 0 && c < (int) chans_.size(); }
    static constexpr float kHole = std::numeric_limits<float>::quiet_NaN();

    // --- monotone deques over the last W finite samples (max and min), per channel ---
    WindowEntry& dq (int c, int which, int i) noexcept
    {
        return deques_[(std::size_t) c * 2u * dequeCap_ + (std::size_t) which * dequeCap_ + (std::size_t) i];
    }
    void dqPush (int c, std::int64_t t, float v) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        const int cap = (int) dequeCap_;
        auto at = [cap] (int head, int i) { const int j = head + i; return j >= cap ? j - cap : j; };   // head, i < cap
        // Expire from the front first: without it a deque is trimmed only by a query, and a long flat band issues none,
        // so a monotone stretch of W samples before it plus the band's own entries overran the ring by one (ASan).
        // Every query asks for times >= t - W, so nothing a query can see is dropped, and the ring holds <= W + 1 entries.
        while (ch.dqMaxCount > 0 && dq (c, 0, ch.dqMaxHead).t < t - W_) { ch.dqMaxHead = at (ch.dqMaxHead, 1); --ch.dqMaxCount; }
        while (ch.dqMinCount > 0 && dq (c, 1, ch.dqMinHead).t < t - W_) { ch.dqMinHead = at (ch.dqMinHead, 1); --ch.dqMinCount; }
        {   // max: drop from the back everything not larger than v
            while (ch.dqMaxCount > 0 && dq (c, 0, at (ch.dqMaxHead, ch.dqMaxCount - 1)).v <= v) --ch.dqMaxCount;
            dq (c, 0, at (ch.dqMaxHead, ch.dqMaxCount)) = { t, v }; ++ch.dqMaxCount;
        }
        {
            while (ch.dqMinCount > 0 && dq (c, 1, at (ch.dqMinHead, ch.dqMinCount - 1)).v >= v) --ch.dqMinCount;
            dq (c, 1, at (ch.dqMinHead, ch.dqMinCount)) = { t, v }; ++ch.dqMinCount;
        }
    }
    // extremes of the finite samples at times >= from (the deques never hold anything older than the window)
    void dqQuery (int c, std::int64_t from, double& mx, double& mn) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        const int cap = (int) dequeCap_;
        while (ch.dqMaxCount > 0 && dq (c, 0, ch.dqMaxHead).t < from) { ch.dqMaxHead = (ch.dqMaxHead + 1) % cap; --ch.dqMaxCount; }
        while (ch.dqMinCount > 0 && dq (c, 1, ch.dqMinHead).t < from) { ch.dqMinHead = (ch.dqMinHead + 1) % cap; --ch.dqMinCount; }
        mx = ch.dqMaxCount > 0 ? (double) dq (c, 0, ch.dqMaxHead).v : -std::numeric_limits<double>::infinity();
        mn = ch.dqMinCount > 0 ? (double) dq (c, 1, ch.dqMinHead).v : std::numeric_limits<double>::infinity();
    }

    Pending& pendingAt (int c, int k) noexcept
    {
        const Channel& ch = chans_[(std::size_t) c];
        return pending_[(std::size_t) c * pendingCap_ + (std::size_t) ((ch.pendingHead + k) % (int) pendingCap_)];
    }

    void step (int c, float x, bool fed) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        if (! ch.histInit) { std::fill (std::begin (ch.hist), std::end (ch.hist), kHole); ch.histInit = true; }
        const bool fin = fed && std::isfinite (x);
        if (fed && ! fin) ++ch.nonFinite;
        const float v = fin ? x : kHole;
        if (fin)
        {
            const double a = std::fabs ((double) x);
            if (a > ch.peak) ch.peak = a;
            const double s = ch.dcSum + (double) x;                       // Neumaier
            ch.dcComp += std::fabs (ch.dcSum) >= std::fabs ((double) x) ? (ch.dcSum - s) + (double) x : ((double) x - s) + ch.dcSum;
            ch.dcSum = s;
            ++ch.finite;
            if (std::isfinite (ch.last))
            {
                const double d = std::fabs ((double) x - (double) ch.last);
                if (d > 0.0 && d < ch.minStep) ch.minStep = d;
            }
            if (! core::exactlyEqual (x, 0.0f))
            {
                const int k = gridExponent (x);
                if (k > ch.latticeK)
                {
                    ch.latticeK = k;
                    ch.latticeQ = k <= kMaxLatticeK ? std::ldexp (1.0, -k) : std::numeric_limits<double>::infinity();
                }
            }
            ch.q = std::min (ch.minStep, ch.latticeQ);
        }

        // bands
        if (ch.inBand)
        {
            if (fin)
            {
                const double nlo = std::min (ch.bLo, (double) x), nhi = std::max (ch.bHi, (double) x);
                if (nhi - nlo <= kBandQuanta * ch.q) { ch.bLo = nlo; ch.bHi = nhi; ch.bSum += (double) x; }
                else { closeBand (c, t_); startBand (c, x); }
            }
            else closeBand (c, t_);
        }
        else if (fin) startBand (c, x);

        if (fin) dqPush (c, t_, x);
        for (int i = kPre - 1; i > 0; --i) ch.hist[i] = ch.hist[i - 1];
        ch.hist[0] = v;
        ch.last = v;

        // the samples after each band's end, and q as of end + 2
        for (int k = ch.pendingCount - 1; k >= 0; --k)
        {
            Pending& p = pendingAt (c, k);
            if (p.postCount >= kPost) break;                             // older ones are complete too
            p.post[p.postCount++] = v;
            if (p.postCount == 3) { p.q = ch.q; p.qTaken = true; }
        }
        while (ch.pendingCount > 0 && pendingAt (c, 0).end + W_ - 1 <= t_) decideFront (c);
    }

    void startBand (int c, float x) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        ch.inBand = true;
        ch.bStart = t_;
        ch.bLo = ch.bHi = ch.bSum = (double) x;
        std::copy (std::begin (ch.hist), std::end (ch.hist), std::begin (ch.bPre));   // y[t-1] … y[t-9], before x is pushed
        dqQuery (c, t_ - W_, ch.bPreMax, ch.bPreMin);
    }

    void closeBand (int c, std::int64_t end) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        ch.inBand = false;
        if (end - ch.bStart < kMinCeilingLength) return;                  // a single sample is never flat
        if (ch.pendingCount == (int) pendingCap_) decideFront (c);        // unreachable by the capacity bound; kept safe
        Pending& p = pendingAt (c, ch.pendingCount);
        ++ch.pendingCount;
        p.start = ch.bStart; p.end = end;
        p.lo = ch.bLo; p.hi = ch.bHi; p.sum = ch.bSum;
        p.preMax = ch.bPreMax; p.preMin = ch.bPreMin;
        p.q = p.finalQ = 0.0; p.qTaken = false;
        p.postCount = 0;
        std::copy (std::begin (ch.bPre), std::end (ch.bPre), std::begin (p.pre));
    }

    void decideFront (int c) noexcept
    {
        Channel& ch = chans_[(std::size_t) c];
        const Pending p = pendingAt (c, 0);
        ch.pendingHead = (ch.pendingHead + 1) % (int) pendingCap_;
        --ch.pendingCount;

        const std::int64_t L = p.end - p.start;
        const double q = p.qTaken ? p.q : p.finalQ;
        if (! std::isfinite (q)) return;
        const double tau = kBandQuanta * q;
        const double m = p.sum / (double) L;
        auto post = [&p] (int i) { return i < p.postCount ? (double) p.post[i] : std::numeric_limits<double>::quiet_NaN(); };
        auto pre = [&p] (int i) { return (double) p.pre[i]; };
        const bool knownIn = std::isfinite (pre (0)), knownOut = std::isfinite (post (0));
        if (! knownIn && ! knownOut) return;
        int type = 0;
        {
            const bool upIn = ! knownIn || pre (0) < p.lo, upOut = ! knownOut || post (0) < p.lo;
            const bool dnIn = ! knownIn || pre (0) > p.hi, dnOut = ! knownOut || post (0) > p.hi;
            if (upIn && upOut) type = 1; else if (dnIn && dnOut) type = -1;
        }
        if (type == 0) return;

        double postMax, postMin;
        dqQuery (c, p.end, postMax, postMin);                             // [end, end + W) — all that has arrived of it
        if (type > 0 ? std::max (p.preMax, postMax) - p.hi > tau : p.lo - std::min (p.preMin, postMin) > tau) return;

        const double s = (double) type;
        auto atLevel = [&] (double y) { return std::fabs (y - m) <= tau; };
        auto stepSide = [&] (double a, double b, double c2) {             // a unknown sample counts as a step
            if (! std::isfinite (a) || ! std::isfinite (b) || ! std::isfinite (c2)) return true;
            return (std::fabs (a - b) <= tau && ! atLevel (b)) || (std::fabs (b - c2) <= tau && ! atLevel (b));
        };
        double act = 0.0; int na = 0;
        for (int i = 0; i < kFlankSteps && std::isfinite (pre (i)) && std::isfinite (pre (i + 1)); ++i)
            if (! atLevel (pre (i)) && ! atLevel (pre (i + 1))) { act += std::fabs (pre (i) - pre (i + 1)); ++na; }
        for (int i = 0; i < kFlankSteps && std::isfinite (post (i)) && std::isfinite (post (i + 1)); ++i)
            if (! atLevel (post (i)) && ! atLevel (post (i + 1))) { act += std::fabs (post (i + 1) - post (i)); ++na; }
        const double sigma = na > 0 ? act / (double) na : 0.0;
        const double stretch = (double) (L + 1) / (double) (L - 1);
        const double V = (tau + q) * (stretch * stretch) + q;

        // how much the outward steps grow, from the band's level: (step2 - step1) / step1 over the three samples of a side
        auto turn = [&] (double a, double b, double c2) {
            const double d1 = s * (m - a), d2 = s * (m - b), d3 = s * (m - c2);
            const double s1 = d2 - d1, s2 = d3 - d2;
            return s1 > 0.0 ? (s2 - s1) / s1 : std::numeric_limits<double>::infinity();
        };
        bool ramp = L >= kMinRampLength && chanceBelow (sigma > 0.0 ? (p.hi - p.lo + q) / sigma : 1.0, L - 1) && knownIn && knownOut
                 && ! stepSide (pre (0), pre (1), pre (2)) && s * (m - pre (0)) >= V
                 && ! stepSide (post (0), post (1), post (2)) && s * (m - post (0)) >= V;
        if (ramp && p.hi > p.lo)
            ramp = turn (pre (0), pre (1), pre (2)) <= kRoughTurnRatio && turn (post (0), post (1), post (2)) <= kRoughTurnRatio;
        const int ci = type > 0 ? 0 : 1;
        const bool atCeiling = ch.ceilingSet[ci] && std::fabs (m - ch.ceiling[ci]) <= tau;
        if (! ramp && ! atCeiling) return;
        if (ramp) { ch.ceiling[ci] = m; ch.ceilingSet[ci] = true; }

        ++ch.runs;
        ch.clippedSamples += L;
        ch.longest = std::max (ch.longest, L);
        if (runCount_ < (std::int64_t) runs_.size())
        {
            ClipRun& r = runs_[(std::size_t) runCount_];
            r.start = p.start; r.length = L; r.level = m; r.channel = c; r.sign = type;
            r.evidence = ramp ? ClipEvidence::Ramp : ClipEvidence::Ceiling;
        }
        ++runCount_;
    }

    double sampleRate_ = 48000.0;
    int channels_ = 0;
    std::int64_t W_ = 0;
    std::size_t dequeCap_ = 0, pendingCap_ = 0;
    bool prepared_ = false, finished_ = false;
    std::int64_t t_ = 0, runCount_ = 0;
    ClipDetectorParams params_;
    std::vector<WindowEntry> deques_;
    std::vector<Pending> pending_;
    std::vector<ClipRun> runs_;
    std::vector<Channel> chans_;
};

} // namespace felitronics::analysis
