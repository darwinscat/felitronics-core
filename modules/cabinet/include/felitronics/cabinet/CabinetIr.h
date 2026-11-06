// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/core/Math.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/measurement/ReferenceUnity.h>

#include <felitronics/convolution/CabConvolver.h>
#include <felitronics/convolution/IrResampler.h>
#include <felitronics/convolution/MatrixConvolverNupc.h>

#include <algorithm>
#include <atomic>
#include <mutex>
#include <numeric>
#include <vector>

namespace felitronics::cabinet
{

/** The cabinet as a convolution: one IR at the end of the chain. Moved here from OrbitAmp (src/core/CabinetIr.h) —
    the DSP half, its arithmetic unchanged: the shot comes in as taps (load), and decoding a file is the caller's.

    The family's own convolver does the heavy lifting: felitronics-core's CabConvolver, the one
    OrbitCab plays — true sample-zero latency, a result that does not depend on the host's block
    size, and a new IR built on the message thread into the idle slot and crossfaded in over 50 ms
    on the audio thread, without a click. It replaced `juce::dsp::Convolution`, which loaded on a
    thread of its own at a pace nobody could see and did not agree across platforms about what a
    bad sample does to it (see `process`).

    WHAT A PLAYER DOES TO THE IR is baked into the IR, not run on the audio: a second-order
    high-pass and low-pass over the impulse itself, a trim of its tail, a flip of its sign — the
    raw shot is kept here and rebuilt whenever one of them moves, then handed to the convolver as
    a new IR. The convolution costs exactly what it cost before, whatever is switched on.

    Level: the family reference-unity contract (felitronics::measurement::referenceUnityGain),
    the same one OrbitCab ships — the IR passes the guitar-band reference at unity RMS, so a cab
    contributes TONE, not gain, and swapping cabinets compares voicings at matched loudness. Taken
    after the post, so a cut IR is as loud as the whole one was. A bass amp's cabinet passes the
    bass's reference instead (Reference::bass, bassReferenceUnityGain). JUCE's Normalise::yes is NOT
    usable here: it scales the IR to 0.125 total energy (−18 dB full-band), which lands a
    band-limited guitar signal ~12 dB low and drifts with the session rate. Never load a cab IR
    raw either — 61k taps sum to enormous gain; the reference gain is clamped ±30 dB. The bass's
    reference is felitronics::measurement::bassReferenceUnityGain.

    THREADS. `prepare`, `load`, `setPost`, `setReference` and `flushPending` are the producer's side and take one
    lock: the convolver allows one producer beside the audio thread, and a host may prepare from a
    thread of its own while the message thread's pump is loading. `process` and `idle` are the
    audio thread's and take none. */
class CabinetIr
{
public:
    struct Post
    {
        bool  trimOn       = false;
        float trimFraction = 1.0f;     // of the IR's length, kept
        bool  hpfOn        = false;
        float hpfHz        = 80.0f;
        int   hpfSlope     = 12;       // dB/oct: 6, 12, 18, 24 or 48 — the consoles' ladder
        bool  lpfOn        = false;
        float lpfHz        = 7000.0f;
        int   lpfSlope     = 12;
        bool  phase        = false;    // flipped

        bool operator== (const Post& o) const noexcept
        {
            return trimOn == o.trimOn && core::approximatelyEqual (trimFraction, o.trimFraction)
                && hpfOn == o.hpfOn && core::approximatelyEqual (hpfHz, o.hpfHz)
                && hpfSlope == o.hpfSlope
                && lpfOn == o.lpfOn && core::approximatelyEqual (lpfHz, o.lpfHz)
                && lpfSlope == o.lpfSlope
                && phase == o.phase;
        }
    };

    void prepare (double sampleRate, int blockSize, int numChannels)
    {
        const std::scoped_lock lock (producer);

        channels = std::clamp (numChannels, 1, 2);
        // Its schedule is fixed at prepare for the longest IR it will ever be handed — `load` refuses
        // longer — and the loudness is ours, set in `rebuild`, so the convolver applies none.
        prepared = conv.prepare (sampleRate, blockSize, channels, maxSeconds, false);
        spare.assign ((size_t) std::max (1, blockSize), 0.0f);
        silence.assign ((size_t) std::max (1, blockSize), 0.0f);
        hostRate    = sampleRate > 0.0 ? sampleRate : 48000.0;   // the convolver's own reading of a rate of nothing
        drainLeft   = 0;
        silenceLeft = 0;
        cleared     = true;
        handedOver.store (false, std::memory_order_relaxed);

        // A prepared convolver holds no cabinet, and the one kept here is forgotten with it: the chain
        // loads the chosen IR again straight after, and a cut that moves in between must not build a
        // cabinet of its own first — the convolver takes one handover at a time, and on a render with
        // no message loop the second one would wait for a retry that never comes.
        raw = {};
        atHost.clear();
    }

    /** The longest IR `load` takes, and the rates it takes one at — see `takes`. */
    static constexpr double maxSeconds = 5.0, minRate = 8000.0, maxRate = 768000.0;

    /** Whether `load` takes a shot this long at this rate. A header can say anything: a rate no converter ever ran at
        would let a tiny file claim billions of samples and still pass the duration test; and a cabinet with its room is
        a second or two — anything longer than `maxSeconds` is not a cabinet (millions of taps the convolution would try
        to run every block, and a rebuild over all of them on every move of a cut). */
    static bool takes (long long numSamples, double sampleRate) noexcept
    {
        return numSamples > 0 && sampleRate >= minRate && sampleRate <= maxRate
            && (double) numSamples <= maxSeconds * sampleRate;
    }

    /** The producer's side. ONE channel's taps — the FIRST of a file's (the convolution plays channel 0 of whatever it
        is handed, so a stereo IR kept whole would be measured across both and heard through one) — at the IR's own
        rate. Copied here: the caller's buffer may die. False, and nothing changed, for a shot `takes` refuses. */
    bool load (const float* taps, int numSamples, double sampleRate)
    {
        if (taps == nullptr || ! takes (numSamples, sampleRate))
            return false;

        std::vector<float> shot (taps, taps + numSamples);

        const std::scoped_lock lock (producer);
        raw     = std::move (shot);
        rawRate = sampleRate;
        atHost.clear();
        rebuild();
        return true;
    }

    /** Message thread. Rebuilds only when something actually moved. */
    void setPost (const Post& p)
    {
        const std::scoped_lock lock (producer);

        if (p == post)
            return;

        post = p;
        if (! raw.empty())
            rebuild();
    }

    const Post& currentPost() const noexcept { return post; }

    /** What the IR's level is taken against: the guitar's band, or the bass's (bassReferenceUnityGain). The producer's side. */
    enum class Reference { guitar, bass };

    void setReference (Reference r)
    {
        const std::scoped_lock lock (producer);

        if (r == reference)
            return;

        reference = r;
        if (! raw.empty())
            rebuild();
    }

    void process (float* const* io, int numChannels, int numSamples, bool on)
    {
        if (! on)
            return;

        const int nch = std::min (channels, numChannels);
        drainLeft = 0;   // it is sounding again: whatever was still being drained is its history now

        // Still flushing a poisoned history (see below): silence goes in, silence comes out.
        if (silenceLeft > 0)
        {
            flushWithSilence (numSamples);
            silenceLeft -= numSamples;
            for (int ch = 0; ch < nch; ++ch)
                std::fill_n (io[ch], numSamples, 0.0f);
            return;
        }

        // POISON. A sample that is not a number does not enter the convolution: a convolver's history
        // is a recursion of blocks, and the one this used to be kept a NaN for good on Windows — the
        // plugin went silent from that sample on. This link replaces the signal it is handed, so a
        // bad sample becomes a silent one here exactly as it would at the output door.
        for (int ch = 0; ch < nch; ++ch)
            for (int i = 0; i < numSamples; ++i)
                if (! std::isfinite (io[ch][i]))
                    io[ch][i] = 0.0f;

        cleared = false;   // there is history in the tail again

        // THE CONVOLVER'S WIDTH IS EXACT: prepared for the bus's two channels, it refuses a call with
        // one — a 2×2 operator needs both planes. The chain runs mono until its stereo seam, so a
        // one-channel call goes in with a silent second plane beside it: that channel's history
        // simply stays quiet, and nothing is torn when the chain's width changes.
        float* planes[2] { io[0], nch > 1 ? io[1] : spare.data() };

        if (nch < channels)
        {
            if (numSamples > (int) spare.size())
                planes[0] = nullptr;   // longer than prepared: refused below
            else
                std::fill_n (spare.data(), numSamples, 0.0f);
        }

        // A convolver that could not be prepared, or a block past what it was prepared for, plays
        // nothing rather than something wrong: the host promised this block size, and the chain
        // cuts longer blocks to it before they reach here.
        if (! prepared || planes[0] == nullptr || ! conv.process (planes, channels, numSamples))
        {
            for (int ch = 0; ch < nch; ++ch)
                std::fill_n (io[ch], numSamples, 0.0f);
            return;
        }

        // ...and one the convolver makes of its own is not kept either: the block goes silent and the
        // convolution starts again. `x - x` is 0 for every finite x and NaN for anything else.
        //
        // Silent for the length of the impulse playing, not of the whole schedule: that is what this
        // cabinet reads back, so after it the output is exact again. A longer IR picked before the
        // rest has been played out may read a poisoned sample from further back — and is caught here
        // the same way.
        float poison = 0.0f;
        for (int ch = 0; ch < nch; ++ch)
            for (int i = 0; i < numSamples; ++i)
                poison += io[ch][i] - io[ch][i];

        if (! std::isfinite (poison))
        {
            silenceLeft = impulseSamples();
            for (int ch = 0; ch < nch; ++ch)
                std::fill_n (io[ch], numSamples, 0.0f);
        }
    }

    /** The chain calls this every block the cabinet is out of the path: what it remembers is DRAINED
        — silence run through it for as long as the convolver can remember anything — and then it rests.

        Not the convolver's reset. That runs on this thread against a load the message thread may be
        handing over, which the convolver forbids, and it cancels a handover it has not picked up yet:
        the new cabinet would stay in the idle slot for good and the one playing would be silence.

        THE WHOLE SCHEDULE, NOT THE IMPULSE. The convolver keeps every stage of its history running
        whatever the IR playing reads of it, so that a longer IR finds a real past — which means a
        drain as long as the short IR that was playing leaves the rest of the phrase in there, and a
        longer cabinet picked while it stands down plays it back when it returns. Five seconds of a
        cabinet's work after each stand-down, then nothing.

        A cabinet handed over while it rests is run in on silence too, for as long as its crossfade
        takes: otherwise it would still be fading in from the old one when the chain brings it back. */
    void idle (int numSamples)
    {
        if (numSamples <= 0)
            return;

        const bool landing = handedOver.exchange (false, std::memory_order_relaxed);

        if (cleared)
        {
            if (! landing)
                return;

            cleared   = false;
            drainLeft = landingSamples();
        }
        else if (drainLeft <= 0)
            drainLeft = scheduleSamples();
        else if (landing)
            drainLeft = std::max (drainLeft, landingSamples());

        flushWithSilence (numSamples);
        drainLeft -= numSamples;

        if (drainLeft <= 0)
            cleared = true;
    }

    /** Message thread, every pump tick. An IR handed over while the previous one is still fading in
        waits in the convolver; this gives it the next chance. Dragging a cut rebuilds the IR on every
        tick, so without it the last position of a drag could be the one that never arrives. */
    void flushPending()
    {
        const std::scoped_lock lock (producer);

        if (conv.hasPending() && conv.flushPending())
            handedOver.store (true, std::memory_order_relaxed);
    }

    /** How long the loaded impulse is, in seconds — this block's whole tail. */
    float tailSeconds() const noexcept { return tailSec.load (std::memory_order_relaxed); }

private:
    using Nupc = felitronics::convolution::MatrixConvolverNupc<felitronics::convolution::CabConvFft>;

    bool cleared = true;   // nothing left in the convolver's history — see idle()
    int  drainLeft = 0, silenceLeft = 0;
    double hostRate = 48000.0;
    std::vector<float> silence;   // the second silent plane of a flush
    std::atomic<bool> handedOver { false };   // a cabinet was published since the audio thread last looked

    /** How much silence clears the whole of the convolver's memory: the schedule it was prepared for,
        which ends past that by up to one of its longest partitions, plus the input frame of one. */
    int scheduleSamples() const noexcept
    {
        return (int) std::ceil (maxSeconds * hostRate) + 2 * Nupc::kDefaultMaxBlock + 1;
    }

    /** How much of it the impulse playing reads back: the impulse, plus a block for the framing. */
    int impulseSamples() const noexcept
    {
        return core::roundToInt ((double) tailSec.load (std::memory_order_relaxed) * hostRate)
             + (int) spare.size() + 1;
    }

    /** How long a handover takes to land: the convolver's 50 ms crossfade, from the block that picks
        it up. */
    int landingSamples() const noexcept
    {
        return core::roundToInt (0.05 * hostRate) + (int) spare.size() + 1;
    }

    /** One block of silence through the convolver, and its answer thrown away. */
    void flushWithSilence (int numSamples)
    {
        if (! prepared || numSamples > (int) spare.size())
            return;

        std::fill_n (spare.data(), numSamples, 0.0f);
        std::fill_n (silence.data(), numSamples, 0.0f);
        float* planes[2] { spare.data(), silence.data() };
        (void) conv.process (planes, channels, numSamples);
    }

    /** The same rate, to within what a host's clock arithmetic leaves behind: 48000.0000001 is 48 kHz,
        and must not send a cabinet through a resampler for nothing. */
    static bool sameRate (double a, double b) noexcept
    {
        return std::abs (a - b) <= 1.0e-6 * std::max (a, b);
    }

    /** The shot at the session's rate, made once per IR and per rate — a cut dragged at 96 kHz used to
        pay for a resample on every move of the handle, a tenth of a second each.

        SILENCE AROUND IT, AND THE PART THAT LANDS ON IT THROWN AWAY. The resampler fills in what lies
        before the first sample and after the last from the samples it has, not with zeros — and a
        cabinet's first sample is its loudest edge, which that turns into a broadband floor 12 to 26 dB
        over the cabinet in the top octave. Handed real zeros on both sides, it has nothing to invent.
        The zeros in front come off again as a whole number of output samples, so the grid is the one
        the bare shot would have had: `lead` is a multiple of the rate's step whenever the two rates
        have one short enough, and within half a sample of it when they do not.

        Its level is the shot's: a resampled impulse carries the same response on more (or fewer)
        taps, and its gain grows with the ratio unless the ratio is taken out — +6 dB for a 48 kHz
        cabinet in a 96 kHz session. */
    std::vector<float> shotAtHost() const
    {
        const int len = (int) raw.size();
        const double ratio = hostRate / rawRate;
        const int reach = 2 * felitronics::convolution::IrResampleConfig{}.halfTaps;   // past the kernel's radius

        int lead = reach;
        if (core::exactlyEqual (std::floor (rawRate), rawRate) && core::exactlyEqual (std::floor (hostRate), hostRate)
            && hostRate <= 16.0e6)
        {
            const auto from = (long long) rawRate, to = (long long) hostRate;
            const auto step = from / std::gcd (from, to);   // input samples per whole output sample
            if (step <= 1 << 16)
                lead = (int) (step * ((reach + step - 1) / step));
        }

        std::vector<float> padded ((size_t) (lead + len + reach), 0.0f);
        std::copy_n (raw.data(), len, padded.begin() + lead);

        const auto out = felitronics::convolution::resampleIr (padded.data(), (int) padded.size(), rawRate, hostRate);
        const int drop = (int) std::llround ((double) lead * ratio);
        const int want = std::max (1, (int) std::llround ((double) len * ratio));
        const float rateGain = (float) (rawRate / hostRate);

        std::vector<float> taps ((size_t) want, 0.0f);
        for (int i = 0; i < want && drop + i < (int) out.size(); ++i)
            taps[(size_t) i] = out[(size_t) (drop + i)] * rateGain;
        return taps;
    }

    /** What the player does to the IR, over taps at `rate`: the trim, the cuts, the sign. */
    void bake (std::vector<float>& ir, double rate) const
    {
        const int full = (int) ir.size();
        const int keep = post.trimOn ? std::clamp (core::roundToInt ((float) full * post.trimFraction),
                                                   std::min (full, 64), full)
                                     : full;
        ir.resize ((size_t) keep);
        const float* const end = ir.data() + keep;
        float* const d = ir.data();

        // The trimmed edge leaves on a half-cosine, a few milliseconds wide: an abrupt cut rings
        // as spectral ripple, and the trim is meant to remove room, not to add zipper.
        if (keep < full)
        {
            const int fade = std::min (keep / 4, core::roundToInt (0.004 * rate));
            for (int i = 0; i < fade; ++i)
            {
                const float t = (float) (i + 1) / (float) fade;
                d[keep - fade + i] *= 0.5f + 0.5f * std::cos (t * (float) core::kPi);
            }
        }

        // The cuts, run over the impulse itself: matched-biquad Butterworth cascades — the same
        // ladder EqLink runs live (odd slopes lead with a first-order section), so the picture,
        // the consoles and the baked IR all mean the same steepness.
        const auto run = [&] (bool isHighpass, float hz, int slopeDb)
        {
            namespace m = felitronics::eq::matched;
            const double f = std::clamp ((double) hz, 10.0, 0.49 * rate);

            felitronics::eq::BiquadCoeffs sections[4];
            int count = 0;

            const auto fo = [&] { return isHighpass ? m::highpass1 (f, rate) : m::lowpass1 (f, rate); };
            const auto bq = [&] (double q) { return isHighpass ? m::highpass (f, rate, q) : m::lowpass (f, rate, q); };

            switch (slopeDb)
            {
                case 6:  sections[count++] = fo(); break;
                default:
                case 12: sections[count++] = bq (0.70710678); break;
                case 18: sections[count++] = fo();
                         sections[count++] = bq (1.0); break;
                case 24: sections[count++] = bq (0.54119610);
                         sections[count++] = bq (1.30656296); break;
                case 48: sections[count++] = bq (0.50979558);
                         sections[count++] = bq (0.60134489);
                         sections[count++] = bq (0.89997622);
                         sections[count++] = bq (2.56291545); break;
            }

            for (int s = 0; s < count; ++s)
            {
                felitronics::eq::Biquad biquad;
                biquad.setCoeffs (sections[s]);
                for (float* x = d; x != end; ++x)
                    *x = biquad.processSample (*x);
            }
        };

        if (post.hpfOn) run (true,  post.hpfHz, post.hpfSlope);
        if (post.lpfOn) run (false, post.lpfHz, post.lpfSlope);

        if (post.phase)
            for (float* x = d; x != end; ++x)
                *x = -*x;
    }

    void rebuild()
    {
        // THE LEVEL IS TAKEN AT THE IR'S OWN RATE, the shot with the player's hand on it: measured at
        // the session's, the reference's band would reach to that session's Nyquist and the same
        // cabinet would come out a tenth of a dB apart at 48 and 96 kHz.
        std::vector<float> own (raw);
        bake (own, rawRate);
        const float gain = reference == Reference::bass ? felitronics::measurement::bassReferenceUnityGain (own.data(), (int) own.size(), rawRate)
                                                        : felitronics::measurement::referenceUnityGain (own.data(), (int) own.size(), rawRate);

        // What this speaker goes on saying after the guitar stops — the plugin has to be able to
        // tell a host, and an impulse is exactly as long as it is. The silence cut below may shorten
        // it further at the tail; taking the length before that errs long, which is the side to be
        // wrong on. Stored where the load happens, read from wherever the host asks.
        tailSec.store ((float) ((double) own.size() / std::max (1.0, rawRate)), std::memory_order_relaxed);

        // ...and the taps at the session's rate, which is the only rate the convolver is ever handed:
        // it never resamples anything itself.
        std::vector<float> taps;
        if (sameRate (rawRate, hostRate))
            taps = std::move (own);
        else
        {
            if (atHost.empty())
                atHost = shotAtHost();
            taps = atHost;
            bake (taps, hostRate);
        }

        for (auto& t : taps)
            t *= gain;

        // SILENCE AT EITHER END IS NOT PART OF THE CABINET. Anything under -80 dBFS before the first
        // sound or after the last is cut before the convolver sees it — the rule the JUCE convolver
        // applied under `Trim::yes`, kept so the cabinet's timing is what it was: a vendor IR with a
        // few milliseconds of pre-roll would otherwise start that late.
        constexpr float threshold = 1.0e-4f;   // -80 dBFS
        const int count = (int) taps.size();
        int first = 0, last = count;

        while (first < count && std::abs (taps[(size_t) first]) < threshold) ++first;
        while (last > first && std::abs (taps[(size_t) (last - 1)]) < threshold) --last;

        // A shot with nothing in it — a silent file, or a trim that kept only its pre-roll — is a
        // cabinet that plays silence, and is loaded as one: skipping it would leave the last cabinet
        // playing under a pick that says otherwise.
        const float zero = 0.0f;
        const float* planes[1] { last > first ? taps.data() + first : &zero };
        conv.loadIR (planes, 1, last > first ? last - first : 1, hostRate);

        if (! conv.hasPending())
            handedOver.store (true, std::memory_order_relaxed);
    }

    std::atomic<float> tailSec { 0.0f };

    std::mutex producer;   // prepare, load, setPost, setReference, flushPending — see the class note
    felitronics::convolution::CabConvolver conv;
    bool prepared = false;
    int channels = 2;
    std::vector<float> spare;   // the silent second plane of a one-channel call — see process

    std::vector<float> raw;          // the shot, as it came — the post is baked into a copy
    double rawRate = 48000.0;
    std::vector<float> atHost;       // the shot at the session's rate — see shotAtHost
    Post post;
    Reference reference = Reference::guitar;
};

} // namespace felitronics::cabinet
