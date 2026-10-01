// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/codecgrid/BackResampler.h>
#include <felitronics/codecgrid/GridCurve.h>
#include <felitronics/codecgrid/GridScan.h>
#include <felitronics/core/Config.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

//==============================================================================
// felitronics::codecgrid::CodecGridDetector — does this programme carry the frame grid of a lossy codec?
//
// WHAT IT MEASURES. A decoder's output, analysed again with the codec's own transform from the same sample, gives
// the quantiser's zeros back; one sample off, it does not. Eight stretches of the programme are analysed with
// seven hypotheses — the MPEG-1 Layer III filterbank, the AAC MDCT with each of its two windows, each at the
// programme's rate and at the other of 44.1 / 48 kHz, and CELT's MDCT at 48 kHz — at every frame offset, and each
// (stretch, hypothesis) is reduced to a READING (GridCurve.h). The report says which grid stands out, how far, in
// how many stretches, and whether its PHASE — the offset counted from sample 0 of the programme, modulo the
// codec's hop — is the same in all of them.
//
// WHAT IT DOES NOT SAY. It does not say "this is an MP3". It reports a grid; a programme that went through the
// same transform with coefficients zeroed by something that is not a codec carries the same grid and gets the
// same report. And `NoGrid` is "no grid of the transforms tried, on whole-sample offsets, in the stretches
// examined" — never "lossless". Measured holes: HE-AAC reads at the threshold (its core runs at half the rate), a
// coded sample under a dense mix is not seen, material cut up after decoding has no single phase, and a
// programme converted to another rate AND cropped has its grid a fraction of a sample off.
//
// THE TWO GROUNDS FOR `Confirmed`, either of which is enough:
//   TwoWindows      two stretches in which the rule finds a grid (GridRule: score, depth, breadth, uniqueness),
//                   of one codec family at one rate, on EXACTLY the same phase;
//   PhaseAgreement  the best offset of one hypothesis lands on one phase (within a sample) in at least four of
//                   the eight stretches — whatever the scores. On a programme without a grid the best offsets
//                   are scattered, and four of eight landing together by chance is about 1e-5 per hypothesis;
//                   measured, 146 lossless programmes never had more than three. This ground is what sees a weak
//                   grid everywhere (HE-AAC from a streaming service: every stretch under the rule's threshold,
//                   all eight on one phase) and a grid half a sample off after a rate conversion and a crop.
// One stretch found and nothing to confirm it is `InPlaces`, and is NOT a finding: the one false alarm of the
// measurement campaign was exactly that — a 1024-sample block structure in two seconds of one programme, and
// 1024 is both AAC's hop and an ordinary buffer size.
//
// WHERE THE NUMBERS COME FROM. The rule, the stretch length (2 s), their number (8, at tenths of the programme)
// and the phase tolerance are the ones a measurement campaign fixed: 148 lossless programmes x 14 encoder settings
// x 2 conditions, then files from the wild. They are parameters here and their defaults are those numbers; a
// caller who changes one has left what was measured.
//
// THE RATES. 44.1 and 48 kHz programmes are examined. Any other rate is accepted and reported `NotExamined`
// (UnsupportedRate) — a measurement that cannot be made is a report, not a refused prepare().
//
// THE CHANNELS. The first two planes are read; a mono programme is read as two identical channels. The four
// signals the evidence is kept for are left, right, mid and side.
//
// STREAMING (law 11). prepare() takes the programme's TOTAL length, because where the stretches are is a function
// of it; it is the only allocation and what it asks for is published first (storageFor). process() takes the
// programme in any split and keeps a copy of the stretches only — the result does not depend on the split. The
// width is exact. A call past the prepared length is refused whole. finish() / finishStep() analyse what was
// seen: a programme that ended early is examined in the stretches that arrived complete.
//
// TIME. Nearly all of it is in finishStep(): one call does a bounded piece — a few dozen offsets of one scan, or
// one resampling, or one reading — and returns false while work remains. Nothing is allocated, locked or thrown
// after prepare(). A programme whose grid the rule confirms early costs two stretches; one without a grid costs
// all eight times seven, which is the price of saying so.
//
// DETERMINISM. Same input, same parameters, same binary and FFT backend: the same report whatever the split.
// Across platforms or backends the transforms' single-precision rounding may differ in the last bits, and so may
// a score in its sixth digit; nothing here is promised bit-identical across them.
//==============================================================================
namespace felitronics::codecgrid
{

enum class Family  : std::uint8_t { None, Mp3, Aac, Celt };
enum class Verdict : std::uint8_t { NotExamined, NoGrid, SeveralGrids, InPlaces, Confirmed };
enum class Ground  : std::uint8_t { None, TwoWindows, PhaseAgreement, Both };
enum class Reason  : std::uint8_t { None, NotFinished, UnsupportedRate, TooShort, Silent };

// How far the analysis goes once the rule has confirmed a grid.
enum class Depth : std::uint8_t
{
    Verdict,        // stop: the verdict is known
    Follow,         // read the remaining stretches with the confirmed family only — where in the programme it is
    Exhaustive      // every hypothesis in every stretch, always
};

constexpr Family familyOf (Transform t) noexcept
{
    return t == Transform::Mp3 ? Family::Mp3 : t == Transform::Celt ? Family::Celt : Family::Aac;
}

struct CodecGridParams
{
    GridRule rule;
    Depth depth = Depth::Follow;
    double windowSeconds = 2.0;         // the length of a stretch; 0.25 .. 8. The rule's numbers were measured at 2 s
    int phaseTolerance = 1;             // samples: two best offsets this close are one phase
    int phaseAgreeWindows = 4;          // this many stretches on one phase confirm a grid on their own
    double silentWindowDb = -70.0;      // a stretch at or under this level is not examined
};

struct Hypothesis
{
    Transform transform = Transform::Mp3;
    int codecRate = 0;
};

// One hypothesis in one stretch.
struct HypothesisReading
{
    bool scanned = false;
    bool found = false;                 // the rule finds a grid
    bool several = false;               // broad and deep but not unique
    int gridPhase = 0;                  // (offset + the stretch's start at the codec's rate) mod hop
    GridReading reading;
};

inline constexpr int kHypotheses = 7;
inline constexpr int kWindows = 8;

// One stretch of the programme.
struct WindowReading
{
    std::uint64_t startFrame = 0;       // of the 2 s that are read, in programme frames
    bool complete = false;              // every sample of it arrived
    bool silent = false;                // at or under silentWindowDb: not examined
    double rmsDb = 0.0;
    int best = -1;                      // the found hypothesis with the largest score, or -1
    HypothesisReading hypotheses[kHypotheses];
};

struct CodecGridResult
{
    Verdict verdict = Verdict::NotExamined;
    Reason reason = Reason::NotFinished;
    Ground ground = Ground::None;

    Family family = Family::None;       // of the grid reported (Confirmed, InPlaces, SeveralGrids)
    Transform transform = Transform::Mp3;   // which transform read it best — for AAC, which WINDOW won, not which encoder
    int codecRate = 0;                  // the rate the grid is at
    int hop = 0;                        // its frame hop in samples at that rate
    int gridPhase = 0;                  // 0 .. hop - 1

    int windows = 0;                    // stretches placed
    int windowsExamined = 0;            // complete and not silent
    int windowsFound = 0;               // in which the rule found the reported family
    int windowsAgreeing = 0;            // on the reported phase, of the reported hypothesis
    // The curve and the zero share that are published are those of ONE reading — the found one with the largest
    // score, or with nothing found the largest score — and it names itself: with two grids in a programme it
    // need not be the reported one.
    int bestWindow = -1;                // its stretch
    int bestHypothesis = -1;            // its hypothesis: CodecGridDetector::hypothesis (bestHypothesis)
    int bestOffset = 0;                 // its offset, 0 .. hop - 1
    double bestScore = 0.0;

    float zeroByBand[kZeroBands] {};    // the zero share at the found offset of bestWindow, 32 bands to the top of the transform
    float zeroShare = 0.0f;             // its mean over bands 2 .. 15 — coarser coding, more zeros
    std::uint64_t nonFiniteSamples = 0; // NaN or infinite samples met in the stretches kept (read as zero); one that lies in
                                        // two stretches' margins is met twice
};

template <ComplexFftBackend Fft = MixedRadixFft>
class CodecGridDetector
{
public:
    static constexpr double kMinWindowSeconds = 0.25, kMaxWindowSeconds = 8.0;
    static constexpr double kMarginSeconds = 0.25;       // kept each side of a stretch for the rate conversion
    static constexpr int kMaxZeroRows = 192;             // frames of the zero map kept: 2 s holds at most 167 granules
    static constexpr int kOffsetsPerStep = 32;           // one finishStep() of a scan
    static constexpr double kMaxSampleRate = 768000.0;
    static constexpr std::uint64_t kMaxFrames = std::uint64_t (1) << 40;

    //==============================================================================
    // WHAT prepare() ASKS THE HEAP FOR (law 11d), from the function prepare() sizes and validates itself with.
    struct Storage
    {
        bool ok = false;                    // false on exactly the arguments prepare() refuses
        bool examined = false;              // false when the programme cannot be examined (rate, length): nothing is allocated
        int windows = 0;
        std::size_t windowFloats = 0;       // per stretch: 2 channels x (2 s + two margins)
        std::size_t resampleFloats = 0;     // 2 channels x 2 s at the other rate
        std::size_t curveFloats = 0;        // one curve being filled, one kept
        std::size_t readFloats = 0, readDoubles = 0;
        std::size_t zeroMapBytes = 0;
        std::size_t resamplerFloats = 0;    // the taps of the one conversion this rate needs
        std::size_t resamplerDesignDoubles = 0;   // ...and the scratch they are designed in, freed before prepare() returns
        std::uint64_t scanBytes = 0;
        std::uint64_t bytes() const noexcept
        {
            return sizeof (float) * ((std::uint64_t) windows * windowFloats + resampleFloats + curveFloats + readFloats + resamplerFloats)
                 + sizeof (double) * ((std::uint64_t) readDoubles + resamplerDesignDoubles) + zeroMapBytes + scanBytes;
        }
    };

    // The geometry of a rate: the stretch, its margin, the grid its start is rounded to.
    struct Geometry
    {
        bool ok = false;
        int rate = 0, other = 0;
        int length = 0, margin = 0, grid = 0;       // at the programme's rate
        int otherLength = 0, otherMargin = 0;       // at the other rate
    };

    // The stretch is a whole number of grid steps (147 samples at 44.1 kHz are 160 at 48 kHz), so that it is a
    // whole number of samples at both rates.
    static constexpr Geometry geometryFor (double sampleRate, double windowSeconds) noexcept
    {
        Geometry g;
        int otherGrid = 0;
        if (sampleRate == 44100.0)      { g.rate = 44100; g.other = 48000; g.grid = 147; otherGrid = 160; }
        else if (sampleRate == 48000.0) { g.rate = 48000; g.other = 44100; g.grid = 160; otherGrid = 147; }
        else return g;
        if (! (windowSeconds >= kMinWindowSeconds && windowSeconds <= kMaxWindowSeconds)) return g;     // NaN fails
        const int steps = (int) (windowSeconds * (double) g.rate / (double) g.grid);
        g.length = steps * g.grid;      g.margin = g.rate / 4;
        g.otherLength = steps * otherGrid;  g.otherMargin = g.other / 4;
        g.ok = true;
        return g;
    }

    static constexpr bool rateExamined (double sampleRate) noexcept { return sampleRate == 44100.0 || sampleRate == 48000.0; }

    // Where stretch i (0 .. kWindows - 1) starts in a programme of `totalFrames`: at (i + 1) tenths of it,
    // rounded down to the grid, and kept a margin away from both ends.
    static constexpr std::uint64_t windowStart (const Geometry& g, std::uint64_t totalFrames, int i) noexcept
    {
        const std::uint64_t grid = (std::uint64_t) g.grid, need = (std::uint64_t) g.length + (std::uint64_t) g.margin;
        std::uint64_t start = (totalFrames / 10u * (std::uint64_t) (i + 1) + totalFrames % 10u * (std::uint64_t) (i + 1) / 10u) / grid * grid;
        const std::uint64_t last = (totalFrames - need) / grid * grid;
        if (start > last) start = last;
        const std::uint64_t first = ((std::uint64_t) g.margin + grid - 1u) / grid * grid;
        if (start < first) start = first;
        return start;
    }

    static Storage storageFor (double sampleRate, int channels, std::uint64_t totalFrames, const CodecGridParams& p) noexcept
    {
        Storage s;
        if (! (sampleRate >= core::kMinSampleRate && sampleRate <= kMaxSampleRate)) return s;      // NaN fails
        if (channels < 1 || channels > core::kMaxChannels) return s;
        if (totalFrames > kMaxFrames) return s;
        if (! (p.windowSeconds >= kMinWindowSeconds && p.windowSeconds <= kMaxWindowSeconds)) return s;
        if (p.phaseTolerance < 0 || p.phaseAgreeWindows < 2) return s;
        s.ok = true;
        const Geometry g = geometryFor (sampleRate, p.windowSeconds);
        if (! g.ok) return s;
        if (totalFrames < (std::uint64_t) g.length + 2u * (std::uint64_t) g.margin + (std::uint64_t) g.grid) return s;
        s.examined = true;
        s.windows = placedWindows (g, totalFrames, nullptr);
        s.windowFloats = 2u * ((std::size_t) g.length + 2u * (std::size_t) g.margin);
        s.resampleFloats = 2u * (std::size_t) g.otherLength;
        s.curveFloats = 2u * (std::size_t) GridScan<Fft>::kMaxLines * kCells;
        s.readFloats = ReadScratch::floatsFor (GridScan<Fft>::kMaxLines);
        s.readDoubles = ReadScratch::doublesFor (GridScan<Fft>::kMaxLines);
        s.zeroMapBytes = (std::size_t) kMaxZeroRows * kZeroBands;
        const auto ratio = BackResampler::ratioFor (g.rate, g.other);
        s.resamplerFloats = BackResampler::tapsFor (ratio.up, ratio.down);
        s.resamplerDesignDoubles = s.resamplerFloats;
        s.scanBytes = GridScan<Fft>::bytesFor (std::max (g.length, g.otherLength));
        return s;
    }

    void setParams (const CodecGridParams& p) noexcept { params_ = p; }           // takes effect at the next prepare()
    const CodecGridParams& params() const noexcept { return installed_; }

    [[nodiscard]] bool prepare (double sampleRate, int channels, std::uint64_t totalFrames)
    {
        // law 11b: disarm, validate, write
        prepared_ = false;
        finished_ = false;
        stage_ = Stage::Ready;
        channels_ = 0;
        total_ = 0;
        seen_ = 0;
        windows_ = 0;
        geometry_ = Geometry {};
        installed_ = CodecGridParams {};
        clearReport();
        const Storage st = storageFor (sampleRate, channels, totalFrames, params_);
        if (! st.ok) return false;

        installed_ = params_;
        channels_ = channels;
        total_ = totalFrames;
        examined_ = st.examined;
        notExamined_ = st.examined ? Reason::None : rateExamined (sampleRate) ? Reason::TooShort : Reason::UnsupportedRate;
        if (st.examined)
        {
            geometry_ = geometryFor (sampleRate, installed_.windowSeconds);
            windows_ = placedWindows (geometry_, totalFrames, starts_);
            for (int w = 0; w < kWindows; ++w) pcm_[w].assign (w < windows_ ? st.windowFloats : 0u, 0.0f);
            resampled_.assign (st.resampleFloats, 0.0f);
            curve_.assign (st.curveFloats / 2u, 0.0f);
            bestCurve_.assign (st.curveFloats / 2u, 0.0f);
            readFloats_.assign (st.readFloats, 0.0f);
            readDoubles_.assign (st.readDoubles, 0.0);
            zeroMap_.assign (st.zeroMapBytes, std::uint8_t (0));
            const auto ratio = BackResampler::ratioFor (geometry_.rate, geometry_.other);
            if (! resampler_.prepare (ratio.up, ratio.down)) return false;
            if (! scan_.prepare (std::max (geometry_.length, geometry_.otherLength))) return false;
            // the seven hypotheses, in the order they are tried
            const Transform order[3] { Transform::Mp3, Transform::AacSine, Transform::AacKbd };
            int h = 0;
            for (Transform t : order) hypotheses_[h++] = { t, geometry_.rate };
            for (Transform t : order) hypotheses_[h++] = { t, geometry_.other };
            hypotheses_[h] = { Transform::Celt, 48000 };
        }
        prepared_ = true;
        reset();
        return true;
    }

    // Back to the head of the same programme: nothing seen, nothing reported.
    void reset() noexcept
    {
        finished_ = false;
        stage_ = Stage::Ready;
        seen_ = 0;
        nonFinite_ = 0;
        for (int w = 0; w < kWindows; ++w) std::fill (pcm_[w].begin(), pcm_[w].end(), 0.0f);
        clearReport();
    }

    //==============================================================================
    // Law 11 order: malformed -> unprepared or finishing -> the width (EXACT) -> n == 0 -> a plane missing ->
    // past the prepared length -> run. A refused call consumes nothing.
    [[nodiscard]] bool process (const float* const* in, int numChannels, int n) noexcept
    {
        if (numChannels < 0 || n < 0) return false;
        if (numChannels > 0 && in == nullptr) return false;
        if (! prepared_ || stage_ != Stage::Ready) return false;
        if (numChannels != channels_) return false;
        if (n == 0) return true;
        for (int c = 0; c < channels_; ++c) if (in[c] == nullptr) return false;
        if ((std::uint64_t) n > total_ - seen_) return false;

        const std::uint64_t a = seen_, b = seen_ + (std::uint64_t) n;
        const float* left = in[0];
        const float* right = channels_ > 1 ? in[1] : in[0];
        for (int w = 0; w < windows_; ++w)
        {
            const std::uint64_t lo = starts_[w] - (std::uint64_t) geometry_.margin;
            const std::uint64_t hi = starts_[w] + (std::uint64_t) geometry_.length + (std::uint64_t) geometry_.margin;
            const std::uint64_t from = std::max (a, lo), to = std::min (b, hi);
            if (from >= to) continue;
            const std::size_t extended = (std::size_t) (hi - lo);
            float* dl = pcm_[w].data() + (std::size_t) (from - lo);
            float* dr = dl + extended;
            const float* sl = left + (std::size_t) (from - a);
            const float* sr = right + (std::size_t) (from - a);
            const std::size_t count = (std::size_t) (to - from);
            for (std::size_t i = 0; i < count; ++i)
            {
                const float l = sl[i], r = sr[i];
                const bool lf = std::isfinite (l), rf = std::isfinite (r);
                dl[i] = lf ? l : 0.0f;
                dr[i] = rf ? r : 0.0f;
                nonFinite_ += (lf ? 0u : 1u) + (rf ? 0u : 1u);
            }
        }
        seen_ = b;
        return true;
    }

    // End of the programme: drive the same bounded stages a session uses. Idempotent.
    [[nodiscard]] bool finish() noexcept
    {
        if (! prepared_) return false;
        while (! finishStep()) {}
        return true;
    }

    // One bounded unit of the analysis. False while work remains, and when the detector is not prepared.
    [[nodiscard]] bool finishStep() noexcept
    {
        if (! prepared_) return false;
        if (finished_) return true;
        switch (stage_)
        {
            case Stage::Ready:      beginAnalysis(); break;
            case Stage::NextScan:   nextScan(); break;
            case Stage::Resample:   resampleWindow(); break;
            case Stage::Scan:       if (scan_.step (kOffsetsPerStep)) stage_ = Stage::Read; break;
            case Stage::Read:       readScan(); break;
            case Stage::Finalize:   finalize(); break;
            case Stage::Done:       break;
        }
        return finished_;
    }

    bool finished() const noexcept { return finished_; }

    // 0 .. 1: the share of the planned scans that are done. The plan shrinks when the rule confirms a grid.
    double progress() const noexcept
    {
        if (finished_) return 1.0;
        if (! prepared_ || stage_ == Stage::Ready) return 0.0;
        const int planned = std::max (1, plannedScans_);
        const double within = scan_.offsets() > 0 && stage_ == Stage::Scan ? (double) scan_.offsetsDone() / (double) scan_.offsets() : 0.0;
        return std::min (1.0, ((double) doneScans_ + within) / (double) planned);
    }

    //==============================================================================
    // The report. Before the analysis has finished it says NotExamined / NotFinished and nothing else.
    const CodecGridResult& result() const noexcept { return result_; }
    int windows() const noexcept { return windows_; }
    const WindowReading& window (int i) const noexcept { return readings_[i >= 0 && i < kWindows ? i : 0]; }
    static constexpr int hypothesisCount() noexcept { return kHypotheses; }
    Hypothesis hypothesis (int h) const noexcept { return h >= 0 && h < kHypotheses ? hypotheses_[h] : Hypothesis {}; }

    // The curve of result().bestHypothesis in result().bestWindow: curveOffsets() x kCells floats,
    // [offset][signal][group], dB. Empty (0 offsets) when nothing was scanned.
    const float* curve() const noexcept { return bestCurve_.data(); }
    int curveOffsets() const noexcept { return curveOffsets_; }
    // The zero map of the same stretch at its found offset: zeroMapRows() frames x 32 bands, the share of zeros 0..255.
    const std::uint8_t* zeroMap() const noexcept { return zeroMap_.data(); }
    int zeroMapRows() const noexcept { return zeroMapRows_; }

private:
    enum class Stage : std::uint8_t { Ready, NextScan, Resample, Scan, Read, Finalize, Done };

    // The stretches of a programme that do not overlap one another, in order. Returns how many; fills `starts`.
    static constexpr int placedWindows (const Geometry& g, std::uint64_t totalFrames, std::uint64_t* starts) noexcept
    {
        int count = 0;
        std::uint64_t previous = 0;
        for (int i = 0; i < kWindows; ++i)
        {
            const std::uint64_t start = windowStart (g, totalFrames, i);
            if (count > 0 && start < previous + (std::uint64_t) g.length) continue;      // would share samples with the one before
            if (starts != nullptr) starts[count] = start;
            previous = start;
            ++count;
        }
        return count;
    }

    void clearReport() noexcept
    {
        result_ = CodecGridResult {};
        for (auto& r : readings_) r = WindowReading {};
        curveOffsets_ = 0;
        zeroMapRows_ = 0;
        bestHeld_ = false;
        bestHeldFound_ = false;
        bestHeldScore_ = 0.0;
        doneScans_ = 0;
        plannedScans_ = 0;
        confirmed_ = false;
    }

    int hopFor (int h) const noexcept { return hopOf (hypotheses_[h].transform); }

    // The start of stretch w in samples at a hypothesis's rate: whole, because the start is on the grid.
    std::uint64_t startAtRate (int w, int codecRate) const noexcept
    {
        return codecRate == geometry_.rate ? starts_[w]
                                           : starts_[w] / (std::uint64_t) geometry_.grid * (std::uint64_t) (geometry_.grid == 147 ? 160 : 147);
    }

    void beginAnalysis() noexcept
    {
        clearReport();
        result_.nonFiniteSamples = nonFinite_;
        result_.windows = windows_;
        if (! examined_)
        {
            finishWith (Verdict::NotExamined, notExamined_);
            return;
        }
        int usable = 0;
        for (int w = 0; w < windows_; ++w)
        {
            WindowReading& r = readings_[w];
            r.startFrame = starts_[w];
            const std::uint64_t end = starts_[w] + (std::uint64_t) geometry_.length + (std::uint64_t) geometry_.margin;
            r.complete = seen_ >= end;
            if (! r.complete) continue;
            // the level of the 2 s that are read, both channels
            const std::size_t extended = (std::size_t) geometry_.length + 2u * (std::size_t) geometry_.margin;
            const float* l = pcm_[w].data() + (std::size_t) geometry_.margin;
            const float* rr = l + extended;
            double power = 0.0;
            for (int i = 0; i < geometry_.length; ++i) power += (double) l[i] * (double) l[i] + (double) rr[i] * (double) rr[i];
            r.rmsDb = 10.0 * core::det::log10 (power / (2.0 * (double) geometry_.length) + 1.0e-30);
            r.silent = ! (r.rmsDb > installed_.silentWindowDb);
            if (! r.silent) ++usable;
        }
        result_.windowsExamined = usable;
        if (usable == 0)
        {
            bool any = false;
            for (int w = 0; w < windows_; ++w) any = any || readings_[w].complete;
            finishWith (Verdict::NotExamined, any ? Reason::Silent : Reason::TooShort);
            return;
        }
        plannedScans_ = usable * kHypotheses;
        window_ = -1;
        advanceWindow();
    }

    // Move to the next usable stretch and order its hypotheses: those the rule has already found somewhere first,
    // the best score among them first of all — so that the stretch which confirms a grid confirms it with the
    // hypothesis that reads it best — then the rest in their own order.
    void advanceWindow() noexcept
    {
        do { ++window_; } while (window_ < windows_ && (! readings_[window_].complete || readings_[window_].silent));
        if (window_ >= windows_) { stage_ = Stage::Finalize; return; }
        double top[kHypotheses] {};
        bool seen[kHypotheses] {};
        for (int h = 0; h < kHypotheses; ++h)
            for (int w = 0; w < window_; ++w)
                if (readings_[w].hypotheses[h].found)
                {
                    seen[h] = true;
                    top[h] = std::max (top[h], readings_[w].hypotheses[h].reading.score);
                }
        int n = 0;
        for (int h = 0; h < kHypotheses; ++h) if (seen[h]) order_[n++] = h;
        std::stable_sort (order_, order_ + n, [&top] (int a, int b) noexcept { return top[a] > top[b]; });
        for (int h = 0; h < kHypotheses; ++h) if (! seen[h]) order_[n++] = h;
        position_ = 0;
        resampledWindow_ = -1;
        stage_ = Stage::NextScan;
    }

    bool wanted (int h) const noexcept
    {
        if (! confirmed_ || installed_.depth == Depth::Exhaustive) return true;
        if (installed_.depth == Depth::Verdict) return false;
        return familyOf (hypotheses_[h].transform) == result_.family && hypotheses_[h].codecRate == result_.codecRate;
    }

    void nextScan() noexcept
    {
        while (position_ < kHypotheses && ! wanted (order_[position_])) ++position_;
        if (position_ >= kHypotheses) { advanceWindow(); return; }
        current_ = order_[position_];
        if (hypotheses_[current_].codecRate != geometry_.rate && resampledWindow_ != window_) { stage_ = Stage::Resample; return; }
        beginScan();
    }

    void resampleWindow() noexcept
    {
        const std::size_t extended = (std::size_t) geometry_.length + 2u * (std::size_t) geometry_.margin;
        const float* l = pcm_[window_].data();
        resampler_.resample (l, (int) extended, geometry_.otherMargin, geometry_.otherLength, resampled_.data());
        resampler_.resample (l + extended, (int) extended, geometry_.otherMargin, geometry_.otherLength, resampled_.data() + (std::size_t) geometry_.otherLength);
        resampledWindow_ = window_;
        beginScan();
    }

    void beginScan() noexcept
    {
        const Hypothesis hy = hypotheses_[current_];
        const std::size_t extended = (std::size_t) geometry_.length + 2u * (std::size_t) geometry_.margin;
        const float* l;
        const float* r;
        int n;
        if (hy.codecRate == geometry_.rate)
        {
            l = pcm_[window_].data() + (std::size_t) geometry_.margin;
            r = l + extended;
            n = geometry_.length;
        }
        else
        {
            l = resampled_.data();
            r = l + (std::size_t) geometry_.otherLength;
            n = geometry_.otherLength;
        }
        if (! scan_.begin (hy.transform, l, r, n, curve_.data()))
        {
            ++doneScans_;                    // cannot happen for the geometry prepare() accepted; counted, never looped on
            ++position_;
            stage_ = Stage::NextScan;
            return;
        }
        stage_ = Stage::Scan;
    }

    void readScan() noexcept
    {
        const int hop = hopFor (current_);
        ReadScratch scratch;
        scratch.local = readDoubles_.data();
        scratch.sort = readDoubles_.data() + hop;
        scratch.dip = readFloats_.data();
        HypothesisReading& hr = readings_[window_].hypotheses[current_];
        hr.reading = readCurve (curve_.data(), hop, scratch);
        hr.scanned = true;
        hr.found = installed_.rule.found (hr.reading);
        hr.several = installed_.rule.several (hr.reading);
        hr.gridPhase = (int) (((std::uint64_t) hr.reading.offset + startAtRate (window_, hypotheses_[current_].codecRate)) % (std::uint64_t) hop);

        WindowReading& wr = readings_[window_];
        if (hr.found && (wr.best < 0 || hr.reading.score > wr.hypotheses[wr.best].reading.score)) wr.best = current_;

        // the curve and the zero share that are published: a found reading beats any that is not, then the score
        const bool better = (hr.found && ! bestHeldFound_) || (hr.found == bestHeldFound_ && (! bestHeld_ || hr.reading.score > bestHeldScore_));
        if (better)
        {
            std::copy (curve_.begin(), curve_.begin() + (std::ptrdiff_t) ((std::size_t) hop * kCells), bestCurve_.begin());
            curveOffsets_ = hop;
            zeroMapRows_ = scan_.zeroProfile (hr.reading.offset, result_.zeroByBand, zeroMap_.data(), kMaxZeroRows);
            bestHeld_ = true;
            bestHeldFound_ = hr.found;
            bestHeldScore_ = hr.reading.score;
            heldWindow_ = window_;
            heldHypothesis_ = current_;
            heldOffset_ = hr.reading.offset;
        }

        ++doneScans_;
        ++position_;
        stage_ = Stage::NextScan;

        if (! confirmed_ && hr.found && twoWindows (current_, window_))
        {
            confirmed_ = true;
            result_.family = familyOf (hypotheses_[current_].transform);
            result_.codecRate = hypotheses_[current_].codecRate;
            if (installed_.depth == Depth::Verdict) { stage_ = Stage::Finalize; return; }
            if (installed_.depth == Depth::Follow)
            {
                // what is left: this stretch's and the later stretches' hypotheses of the confirmed family
                int left = 0;
                for (int w = window_; w < windows_; ++w)
                {
                    if (! readings_[w].complete || readings_[w].silent) continue;
                    for (int h = 0; h < kHypotheses; ++h)
                        if (wanted (h) && ! readings_[w].hypotheses[h].scanned) ++left;
                }
                plannedScans_ = doneScans_ + left;
            }
        }
    }

    // Is hypothesis h, found in stretch w, found by the same family at the same rate on exactly the same phase in
    // an earlier stretch?
    bool twoWindows (int h, int w) const noexcept
    {
        const Family f = familyOf (hypotheses_[h].transform);
        const int rate = hypotheses_[h].codecRate, phase = readings_[w].hypotheses[h].gridPhase;
        for (int v = 0; v < w; ++v)
            for (int k = 0; k < kHypotheses; ++k)
            {
                const HypothesisReading& other = readings_[v].hypotheses[k];
                if (other.found && familyOf (hypotheses_[k].transform) == f && hypotheses_[k].codecRate == rate && other.gridPhase == phase) return true;
            }
        return false;
    }

    // The largest number of stretches whose best offset of hypothesis h lands on one phase, and that phase.
    int agreement (int h, int& phase) const noexcept
    {
        const int hop = hopFor (h);
        int most = 0;
        phase = 0;
        for (int w = 0; w < windows_; ++w)
        {
            if (! readings_[w].hypotheses[h].scanned) continue;
            const int p = readings_[w].hypotheses[h].gridPhase;
            int count = 0;
            for (int v = 0; v < windows_; ++v)
                if (readings_[v].hypotheses[h].scanned && phasesAgree (readings_[v].hypotheses[h].gridPhase, p, hop, installed_.phaseTolerance)) ++count;
            if (count > most) { most = count; phase = p; }
        }
        return most;
    }

    bool sameGrid (int h, Family family, int rate) const noexcept
    {
        return familyOf (hypotheses_[h].transform) == family && hypotheses_[h].codecRate == rate;
    }

    void finalize() noexcept
    {
        // THE PHASE GROUND: per hypothesis, over the stretches it was read in. With a family the rule already
        // confirmed only its own hypotheses may speak; otherwise the most stretches wins, then the larger score.
        int agreeHyp = -1, agreeCount = 0, agreePhase = 0;
        double agreeScore = 0.0;
        for (int h = 0; h < kHypotheses; ++h)
        {
            if (confirmed_ && ! sameGrid (h, result_.family, result_.codecRate)) continue;
            int phase = 0;
            const int count = agreement (h, phase);
            if (count < installed_.phaseAgreeWindows) continue;
            double score = 0.0;
            for (int w = 0; w < windows_; ++w)
                if (readings_[w].hypotheses[h].scanned && phasesAgree (readings_[w].hypotheses[h].gridPhase, phase, hopFor (h), installed_.phaseTolerance))
                    score += readings_[w].hypotheses[h].reading.score;
            if (count > agreeCount || (count == agreeCount && score > agreeScore)) { agreeHyp = h; agreeCount = count; agreePhase = phase; agreeScore = score; }
        }

        Verdict verdict = Verdict::NoGrid;
        int reported = -1, phase = 0;
        if (confirmed_)
        {
            // THE RULE'S GROUND. The phase is the one most stretches of the family were found on; the hypothesis
            // reported is the one found most often there (for AAC: which window read it best), then the best score.
            verdict = Verdict::Confirmed;
            result_.ground = agreeHyp >= 0 ? Ground::Both : Ground::TwoWindows;
            int most = 0;
            for (int w = 0; w < windows_; ++w)
                for (int h = 0; h < kHypotheses; ++h)
                {
                    if (! sameGrid (h, result_.family, result_.codecRate) || ! readings_[w].hypotheses[h].found) continue;
                    const int p = readings_[w].hypotheses[h].gridPhase;
                    int count = 0;
                    for (int v = 0; v < windows_; ++v)
                    {
                        bool on = false;
                        for (int k = 0; k < kHypotheses; ++k)
                            on = on || (sameGrid (k, result_.family, result_.codecRate) && readings_[v].hypotheses[k].found && readings_[v].hypotheses[k].gridPhase == p);
                        if (on) ++count;
                    }
                    if (count > most) { most = count; phase = p; }
                }
            int best = -1;
            double top = 0.0;
            for (int h = 0; h < kHypotheses; ++h)
            {
                if (! sameGrid (h, result_.family, result_.codecRate)) continue;
                int found = 0;
                double score = 0.0;
                for (int w = 0; w < windows_; ++w)
                    if (readings_[w].hypotheses[h].found && readings_[w].hypotheses[h].gridPhase == phase)
                    {
                        ++found;
                        score = std::max (score, readings_[w].hypotheses[h].reading.score);
                    }
                if (found > best || (found == best && score > top)) { best = found; top = score; reported = h; }
            }
        }
        else if (agreeHyp >= 0)
        {
            verdict = Verdict::Confirmed;
            result_.ground = Ground::PhaseAgreement;
            reported = agreeHyp;
            phase = agreePhase;
        }
        else
        {
            // nothing confirmed: the best reading the rule found is `InPlaces`; failing that the best one that is
            // broad and deep but not unique is `SeveralGrids`
            double top = 0.0;
            for (int pass = 0; pass < 2 && reported < 0; ++pass)
                for (int w = 0; w < windows_; ++w)
                    for (int h = 0; h < kHypotheses; ++h)
                    {
                        const HypothesisReading& hr = readings_[w].hypotheses[h];
                        if (! (pass == 0 ? hr.found : hr.several)) continue;
                        if (reported < 0 || hr.reading.score > top)
                        {
                            verdict = pass == 0 ? Verdict::InPlaces : Verdict::SeveralGrids;
                            top = hr.reading.score;
                            reported = h;
                            phase = hr.gridPhase;
                        }
                    }
        }

        if (reported >= 0)
        {
            result_.family = familyOf (hypotheses_[reported].transform);
            result_.transform = hypotheses_[reported].transform;
            result_.codecRate = hypotheses_[reported].codecRate;
            result_.hop = hopFor (reported);
            result_.gridPhase = phase;
            for (int w = 0; w < windows_; ++w)
            {
                bool found = false, on = false;
                for (int h = 0; h < kHypotheses; ++h)
                {
                    if (! sameGrid (h, result_.family, result_.codecRate)) continue;
                    const HypothesisReading& hr = readings_[w].hypotheses[h];
                    found = found || hr.found;
                    on = on || (hr.scanned && phasesAgree (hr.gridPhase, phase, result_.hop, installed_.phaseTolerance));
                }
                if (found) ++result_.windowsFound;
                if (on) ++result_.windowsAgreeing;
            }
        }

        if (bestHeld_)
        {
            result_.bestWindow = heldWindow_;
            result_.bestHypothesis = heldHypothesis_;
            result_.bestOffset = heldOffset_;
            result_.bestScore = bestHeldScore_;
            double sum = 0.0;
            for (int b = 2; b < 16; ++b) sum += (double) result_.zeroByBand[b];
            result_.zeroShare = (float) (sum / 14.0);
        }
        finishWith (verdict, Reason::None);
    }

    void finishWith (Verdict v, Reason r) noexcept
    {
        result_.verdict = v;
        result_.reason = r;
        finished_ = true;
        stage_ = Stage::Done;
    }

    CodecGridParams params_ {}, installed_ {};
    bool prepared_ = false, finished_ = false, examined_ = false, confirmed_ = false;
    Reason notExamined_ = Reason::None;
    bool bestHeld_ = false, bestHeldFound_ = false;
    Stage stage_ = Stage::Ready;
    int channels_ = 0, windows_ = 0, window_ = 0, position_ = 0, current_ = 0, resampledWindow_ = -1;
    int order_[kHypotheses] {}, doneScans_ = 0, plannedScans_ = 0, curveOffsets_ = 0, zeroMapRows_ = 0;
    int heldWindow_ = -1, heldHypothesis_ = 0, heldOffset_ = 0;
    double bestHeldScore_ = 0.0;
    std::uint64_t total_ = 0, seen_ = 0, nonFinite_ = 0;
    std::uint64_t starts_[kWindows] {};
    Geometry geometry_ {};
    Hypothesis hypotheses_[kHypotheses] {};
    WindowReading readings_[kWindows] {};
    CodecGridResult result_ {};
    GridScan<Fft> scan_;
    BackResampler resampler_;
    std::vector<float> pcm_[kWindows], resampled_, curve_, bestCurve_, readFloats_;
    std::vector<double> readDoubles_;
    std::vector<std::uint8_t> zeroMap_;
};

} // namespace felitronics::codecgrid
