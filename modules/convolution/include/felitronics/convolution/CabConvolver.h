// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <felitronics/convolution/MatrixConvolverNupc.h>  // JUCE-free cab convolution backend (non-uniform/Gardner)
#include <felitronics/convolution/IrResampler.h>          // resample-on-load (message thread, one-shot)
#include <felitronics/core/Fft.h>                          // DefaultRealFft — analysis FFT + the scalar fallback backend

#if defined(FELITRONICS_WITH_PFFFT)
 #include <felitronics/fftpffft/PffftRealFft.h>            // SIMD (z-order) audio FFT — the desktop-speed backend
#endif

#include <algorithm>
#include <cmath>
#include <vector>

//==============================================================================
// felitronics::convolution::CabConvolver — the cab IR convolution. The backend is
// `MatrixConvolverNupc` (non-uniform / Gardner partition schedule, block-INDEPENDENT, TRUE
// sample-zero-latency), the JUCE-free cab DSP path. The public methods speak raw floats + sizes,
// so product adapters remain framework-free.
//
// PERF/PROVENANCE: this replaces the interim `juce::dsp::Convolution` backend (felitronics-core#21).
// NUPC keeps the head + doubling stages warm every sample so a live cabinet swap is click-free, yet
// its spectral MACs scale with the ACTUAL IR length (activeParts), so a short cab IR (<1 s of taps —
// the norm) stays cheap (<1% RT). The audio FFT is pffft (SIMD, BSD) when FELITRONICS_WITH_PFFFT is
// set, else the scalar reference. CabConvolver now LIVES in core; OrbitCab consumes this shared type.
//
// IR-LENGTH CAP: unlike juce (which sized its partitions per-IR), NUPC builds a FIXED schedule from
// maxIrSamples at prepare() time. We cap at maxIrSeconds of taps (kNupcHeadPartition head + one uniform
// tail stage): covers every real cabinet IR with headroom while keeping RAM sane. An IR longer than the cap
// convolves only its first maxIrSeconds — a deliberate product cap (cab IRs are <1 s; longer is off-label).
//
// LOUDNESS: every IR is normalized AT LOAD to reference-unity — unity RMS gain for a guitar-band
// reference signal (white noise through a one-pole low-pass at kIrRefShapeHz), computed in the
// frequency domain from the FINAL IR (post-trim, post-resample) and applied as ONE common gain
// across channels (stereo imaging untouched). Vendor cab IRs are peak-normalized bandpasses, so
// their passband sits +13…+18 dB above unity (measured across the factory set) — convolving them
// RAW made the wet path that much hotter than the dry, pushed the output past the ±24 dB trim range
// with auto-level off, and got worse the quieter the input. After normalization a cab contributes
// TONE, not gain: the auto-leveler shrinks to a small signal-dependent corrector and the plugin is
// usable with it off. The normalized taps are handed to the convolver verbatim (NUPC applies no gain
// of its own), so our reference-unity gain is the single, exact loudness contract.
//
// THE OTHER PATH HAS A LOUDNESS CONTRACT TOO, and it is not "no gain at all". normalize=false keeps the
// IR's authored level — and an IR's level is a CONVOLUTION GAIN, which is a sum over taps and therefore
// proportional to how many of them fit into a second. Resampling to the host rate changes that count, so
// keeping the authored level through a resample takes exactly one factor, irSr/hostSr (P68; the
// arithmetic is `convolutionRateGain` in IrResampler.h). It is 1 whenever the rates match. So: one of
// the two gains is applied at load, never both, and `irNormalizationGain()` reports whichever it was.
//
// THREADING: loadIR() runs on the message thread (resample + FFT-domain gain + the convolver's partition
// build all allocate). MatrixConvolverNupc::setIr()/setOperator() build synchronously into the INACTIVE
// slot then publish; process() picks up the new operator and crossfades it in over ~50 ms (click-free), so
// a live swap never touches the audio thread's active operator. While a crossfade is in flight the convolver
// REJECTS a new load (returns false) — we retain the staged taps and retry on the reload poll (flushPending),
// latest-wins (coalescing). Zero latency — latencySamples() == 0.
//==============================================================================
namespace felitronics::convolution
{

// Audio-path real-FFT backend for the convolver: pffft (SIMD, z-order) on desktop when opted in, else the
// header-only scalar reference. The choice MUST be uniform across every TU that sees this header (the two
// instantiations are distinct types) — it is a single target-wide compile definition.
#if defined(FELITRONICS_WITH_PFFFT)
using CabConvFft = felitronics::fftpffft::PffftRealFft;
#else
using CabConvFft = felitronics::core::fft::DefaultRealFft;
#endif

class CabConvolver
{
public:
    // The reference the IRs are normalized against: RMS gain for white noise shaped by a one-pole
    // low-pass at this corner — a deterministic stand-in for guitar-band program material. The
    // same 2 kHz/-18 dBFS reference anchors the NAM-stage level probes, so "unity" means the same
    // thing across the whole product family.
    static constexpr double kIrRefShapeHz = 2000.0;
    // Clamp: a pathological IR can't blast or vanish. ±30 dB clears every REAL library measured
    // (186 commercial IRs: reference gains +12…+24.1 dB — 96 kHz packs run the hottest) with
    // headroom, so the guard only ever bites genuine garbage.
    static constexpr float  kIrNormMinDb  = -30.0f;
    static constexpr float  kIrNormMaxDb  = 30.0f;
    static constexpr double kIrRefFloorDb = -60.0;    // near-silent IR ⇒ keep g = 1 (don't amplify garbage)
    // Two rates this close are ONE rate, and the IR is taken verbatim. A host that reports 48000.0000001 is
    // reporting 48 kHz; resampling for it band-limits a host-rate IR at 0.95 of Nyquist and moves every tap
    // for nothing. Relative, so it means the same thing at every rate. One part per million is far above
    // anything a rate's own arithmetic drifts by and far below any two rates that really differ, and it is
    // what orbit-amp's own `sameRate` already uses — the family agrees on what "the same rate" means.
    static constexpr double kRateMatchTolerance = 1.0e-6;

    // 🔴 THE CEILING ON A HOST RATE. Its job is not taste — 3 MHz is sixteen times the highest rate any
    // DAW offers — but to keep every `(int) f(sampleRate)` below INSIDE int and every buffer it sizes
    // inside memory. Two conversions in prepare() were undefined without it (the IR-sample count and the
    // crossfade length, both from a rate this class accepted unchecked), and `normalizationGain` sized an
    // analysis window from `lround(hostSr_)`.
    //
    // The VALUE is the house convention, not a new number: dynamics::Compressor::kMaxSampleRate and
    // limiter::TruePeakLimiter::kMaxSampleRate are both 3.0e6 with the same stated reason, eq::EqBand
    // refuses past the same figure as a literal, and rigplayer::RigPlayer spells it a fourth time while
    // naming the duplication as its own disease. This is the fifth spelling; consolidating them is still
    // its own task, and convolution links none of those modules, so it cannot simply ask.
    static constexpr double kMaxSampleRate = 3.0e6;

    // The IR budget prepare() will ask the backend for, as a pure function of the two arguments that
    // decide it — PUBLIC because it is the only way to pin the saturation without building the thing it
    // sizes: a budget of kMaxIrSamples is a 1.34 GB partition schedule (measured), which no test tier
    // should allocate to prove an arithmetic clamp. CLAMPED IN DOUBLE, BEFORE THE CAST — the ceiling used
    // to be applied to the result of `(long long) std::ceil(...)`, one step too late, since the conversion
    // is undefined for an argument outside long long's range: on arm64 a `maxIrSeconds` of +inf produced
    // LLONG_MAX and a one-sample crossfade (a hard switch, not a fade), on x86-64 it saturates the other
    // way. Saturating in double first leaves every representable request exactly where it was.
    // TOTAL, not merely documented: this is public, so it answers for every double a caller can type. A
    // comment is not a precondition — `std::min (NaN, ceiling)` returns the NaN (the comparison is false,
    // so it keeps the first argument) and converting that to int is undefined, which is the very defect
    // this function exists to have fixed one line further up. `! (want > 0.0)` is false for a NaN and for
    // everything at or below zero, so both leave through the same door.
    [[nodiscard]] static int maxIrSamplesFor (double sampleRate, double maxIrSeconds) noexcept
    {
        const double want = std::ceil (maxIrSeconds * sampleRate);
        if (! (want > 0.0)) return 0;
        return (int) std::min (want, (double) Conv::kMaxIrSamples);
    }

    // maxIrSeconds sizes the NUPC partition schedule (the fixed IR-length cap — see IR-LENGTH CAP above).
    // The default matches the historical value; cab IRs are far shorter, so the cap is generous headroom.
    // normalize=false skips the reference-unity RMS normalization (LOUDNESS above) — a REVERB IR is a
    // decay tail, not a tone-shaping cab bandpass, so RMS-normalizing its mostly-decayed length would
    // blow up the wet gain. The spring IRs are peak-normalized at bundle time; the Reverb Mix sets level.
    // It does NOT skip the rate factor: a resampled IR is scaled by irSr/hostSr so that "its authored
    // level" means the same loudness at 44.1, 48, 96 and 192 kHz (P68).
    // THAT FACTOR IS NOT CLAMPED, and the guard above it is — deliberately, and the asymmetry is the
    // point. The normalization gain is MEASURED from the IR's own content, so a pathological IR can make
    // it anything and +-30 dB stops that. The rate factor is arithmetic on two numbers the CALLER gave
    // us, and a clamped one would quietly deliver a different filter than the caller asked for. It spans
    // 4.7e-10 to 4.3e9 (-186 to +192.6 dB) — what `resampleIr`'s own length and position gates leave
    // reachable. A WAV's rate is a uint32 and nothing in this family validates it, so a garbage-but-
    // FINITE header is the one broken metadata a real file can carry, and it now plays LOUD where it
    // used to play quiet: a file claiming 352800 Hz on a 48 kHz host is +17.3 dB, 5e6 Hz is +40.4 dB.
    // Both are correct by this contract — those taps really would be that loud at the rate claimed — and
    // both are garbage. A consumer that loads UNTRUSTED files should bound the rate before it gets here,
    // the way orbit-amp's loader already refuses anything outside 8 kHz...768 kHz.
    // Law 11(b): `numChannels` is BINDING, and this convolver's ceiling is 2, not core::kMaxChannels.
    // It used to CLAMP — prepare(..., 4) succeeded silently as a stereo convolver, after which
    // process(io, 4, n) was a perfectly well-formed call that left planes 2 and 3 DRY. An observable
    // refusal in process() is worth nothing if prepare() already agreed to something it cannot do.
    [[nodiscard]] bool prepare (double sampleRate, int maxBlock, int numChannels, double maxIrSeconds = 4.0,
                                bool normalize = true)
    {
        prepared_ = false;                       // law 11: a REFUSED prepare leaves the object unusable,
        // ...AND UNUSABLE MEANS THE WHOLE OBJECT, NOT JUST THE FLAG. The pending-retry geometry used to be
        // cleared only on the way OUT of a successful prepare, so a refusal in between left `pendingRetry_`
        // true — after which `isBusy()` answered true for the life of the object and `flushPending()` could
        // never publish it, because it returns on `! prepared_`. Reachable in four calls: prepare, start a
        // fade, load again so the convolver rejects it, then re-prepare with anything it refuses. The width
        // refusal below could already do this; the rate and duration refusals would have widened it.
        pendingRetry_ = false;
        pendingLen_   = 0;
        pendingNch_   = 0;
        if (numChannels < 1 || numChannels > 2)  // the way Compressor and TruePeakLimiter already do —
            return false;                        // otherwise a rejected re-prepare silently keeps the old one
        // 🔴 THE RATE IS BINDING TOO, and it used to be GUESSED: `sampleRate > 0.0 ? sampleRate : 48000.0`
        // answered a rate it had not been given with the factory one, and then sized a convolver from the
        // answer — a host at 44.1 kHz that asked with a zero got a 48 kHz crossfade and a 48 kHz IR budget
        // and was told the object was ready. Law 11(b) says every argument prepare() takes is binding and a
        // value it cannot honour is refused HERE, which is what Compressor and TruePeakLimiter do.
        // Spelled as a RANGE and positively: `> 0.0` is false for a NaN and `<= kMaxSampleRate` is false for
        // an infinity, so between them they exclude everything std::isfinite would have — and an infinity is
        // what made `(long long) std::ceil (maxIrSeconds * hostSr_)` and `(int) std::lround (0.05 * hostSr_)`
        // below undefined. THE UNKNOWN-RATE RULE BELOW IS NOT THIS RULE: an IR file's broken metadata is
        // data, and loads as is; a host that cannot name its own clock is a broken call.
        if (! (sampleRate > 0.0 && sampleRate <= kMaxSampleRate)) return false;
        hostSr_    = sampleRate;
        channels_  = numChannels;
        maxBlock_  = std::max (1, maxBlock);                // retained for API parity — NUPC is block-independent
        normalize_ = normalize;

        // Fixed schedule: head kNupcHeadPartition + doubling + one uniform tail covering maxIrSamples.
        // CLAMPED IN DOUBLE, BEFORE THE CAST — the ceiling used to be applied to the result of
        // `(long long) std::ceil(...)`, which is one step too late: the conversion is undefined for an
        // argument outside long long's range, so `maxIrSeconds` of 1e300 (or, once the rate was unbounded,
        // an ordinary four seconds at an infinite rate) was undefined BEFORE anything got capped. Saturating
        // in double first leaves every representable request exactly where it was.
        // A DURATION HAS TO BE A NUMBER, AND NOT A NEGATIVE ONE. `std::max (0.0, NaN)` returns its first
        // operand, so a NaN duration used to prepare successfully with a zero IR budget — which is not
        // silence (the backend always allocates its 128-sample head, so the cab played its first 128 taps)
        // but is certainly not what the caller asked for either, and neither is the zero a negative request
        // was quietly turned into. Law 11(b) again: a value prepare() cannot honour is refused here. +inf
        // is honoured and means the backend's own ceiling, which is what every huge finite value already got.
        if (! (maxIrSeconds >= 0.0)) return false;
        const int maxIrSamples = maxIrSamplesFor (hostSr_, maxIrSeconds);
        // juce swapped IRs over a 50 ms crossfade — match it (the migration spec's verified parity value).
        const int crossfadeSamples = std::max (1, (int) std::lround (0.05 * hostSr_));

        prepared_     = convolution_.prepare (kNupcHeadPartition, maxIrSamples, crossfadeSamples, channels_);
        pendingRetry_ = false;
        pendingLen_   = 0;
        pendingNch_   = 0;
        return prepared_;
    }

    void reset() { convolution_.reset(); }

    // The history alone — the audio the caller fed — leaving a staged or fading IR swap exactly where
    // it is. reset() above ENDS a swap in flight by adopting the filter the caller last published (P88;
    // before that it dropped it, and a knob move was lost until the next knob move). Either verb keeps the
    // filter now, and neither races a loader in the convolver below (measured under ThreadSanitizer, not
    // yet the contract — law 11e). What separates them is law 11a:
    // reset() ends the fade, so a restart one block into it and a restart after it settled answer the next
    // programme identically; this verb leaves the fade running, so they do not — measured here, 2143 of
    // 5120 samples a channel differ (worst 4.186e-01) one block into a 50 ms fade. Reach for this one when the swap
    // must be left where it is.
    void clearAudioState() noexcept { convolution_.clearAudioState(); }

    // Load an IR (mono broadcasts to both channels — juce Stereo::yes parity) — normalized to reference-unity
    // (or, with normalize=false, scaled by the rate factor a resample costs — LOUDNESS above), resampled to
    // host rate unless within kRateMatchTolerance of it. Message thread: resample + gain + the convolver's
    // partition build all allocate.
    // AN UNUSABLE RATE IS AN UNKNOWN RATE, and an unknown rate loads the taps AS IS — NaN, zero, negative and
    // both infinities alike. The samples in the file are fine, only its metadata is broken: refusing would drop
    // the cabinet whole and play silence, where as-is at worst plays an impulse of the wrong length.
    // A load that stages nothing (no samples, a null plane, a known rate so far off that resampleIr cannot
    // address the result — its length or its last position past INT_MAX) is IGNORED whole: the playing IR, the
    // staged taps, their gain and any pending retry stay as they were.
    void loadIR (const float* const* samples, int numChannels, int numSamples, double irSampleRate)
    {
        buildAndStage (samples, numChannels, numSamples, irSampleRate);
    }

    // THE GAIN APPLIED to the last loaded IR — the exact linear factor (for references that must null
    // against the engine) and its dB reading (diagnostics). Message thread. It is the reference-unity
    // normalization on the normalized path and the rate factor irSr/hostSr on the verbatim one (P68);
    // the name is historical, the contract is "what was multiplied in", which is what a null needs.
    // Exactly 1.0f whenever nothing was applied, which on the VERBATIM path means the host's own rate
    // (within kRateMatchTolerance) or an unknown one. The normalized path always applies something: an
    // unknown rate there still loads the taps as is and still normalizes them.
    float irNormalizationGain()   const noexcept { return normGain_; }
    float irNormalizationGainDb() const noexcept { return normGainDb_; }

    // The staged taps of the last load that staged any (one that stages nothing leaves them) — resampled
    // to host rate + scaled by irNormalizationGain(), which is what the convolver plays once a pending
    // retry has published, up to the IR-LENGTH CAP: a load longer than maxIrSeconds stages whole and
    // convolves only its first maxIrSeconds, so past the cap these taps are staged, not audible.
    // (Retained anyway for the reject-retry coalescing.) Message thread;
    // input for offline blend analysis (auto-polarity / interference tint). NOTE: they persist after a
    // slot clear (the engine only gates the slot off) — callers gate on the slot's own loaded state, not
    // on non-emptiness here.
    const std::vector<std::vector<float>>& stagedTaps() const noexcept { return ir_; }

    // RT-safe in-place convolution of `numChannels` planar channels. NUPC processes the prepared channel
    // count and no-ops if handed fewer planes — callers pass the prepared width (wet buffer == bus width).
    [[nodiscard]] bool process (float* const* io, int numChannels, int numSamples)
    {
        if (! prepared_) return false;                                   // ...and process() has to READ that flag
        return convolution_.process (io, io, numChannels, numSamples);   // in-place (in == out), zero latency
    }

    // "A load is in flight": either the async crossfade hasn't finished, or a load is still queued because
    // the convolver rejected it mid-crossfade (retried by flushPending).
    bool isBusy() const noexcept { return pendingRetry_ || convolution_.isBusy(); }

    static constexpr int latencySamples() noexcept { return 0; }   // NUPC guarantees sample-zero latency

    // Message thread, driven by the reload poll. If a load was rejected while the convolver was
    // mid-crossfade, retry it now from the retained taps. Returns true when nothing is pending.
    bool flushPending() { return tryPublishPending(); }
    bool hasPending() const noexcept { return pendingRetry_; }

private:
    using Conv = felitronics::convolution::MatrixConvolverNupc<CabConvFft>;
    static constexpr int kNupcHeadPartition = 128;   // time-domain head P0 (pow2) — lineareq uses the same

    // STAGED IN LOCALS, COMMITTED WHOLE. The taps, their gain and the pending retry's geometry change together
    // or not at all. This used to overwrite `ir_` first and return on an empty result after — and with a
    // rejected load still pending, the retry then published the OLD length from an emptied or narrowed `ir_`:
    // a read past the end of a vector (a zero-length load; a mono load over a pending stereo one), or, with
    // a null data pointer, a retry refused forever — isBusy() stuck true, and neither the pending IR nor the
    // new one ever reached the convolver.
    void buildAndStage (const float* const* samples, int nch, int len, double irSr)
    {
        if (samples == nullptr || len <= 0) return;
        nch = std::clamp (nch, 1, 2);

        // Only a KNOWN rate — a positive finite number — that is off the host's by more than the tolerance.
        const bool resample = irSr > 0.0 && std::isfinite (irSr)
                           && std::fabs (irSr - hostSr_) > kRateMatchTolerance * std::max (irSr, hostSr_);
        std::vector<std::vector<float>> staged ((std::size_t) nch);
        for (int c = 0; c < nch; ++c)                                          // resample only off host rate
        {
            if (samples[c] == nullptr) return;
            auto& ch = staged[(std::size_t) c];
            if (resample)
                ch = felitronics::convolution::resampleIr (samples[c], len, irSr, hostSr_);
            else
                ch.assign (samples[c], samples[c] + len);
            if (ch.empty()) return;                                            // a result resampleIr cannot address
        }
        const int outLen = (int) staged[0].size();

        // Reference-unity normalization of the FINAL (resampled) IR — one common gain, all channels.
        // Skipped for a reverb IR (normalize_=false) — but a RESAMPLED reverb IR still gets the rate
        // factor, which is a different thing and not a normalization: the taps go in at their AUTHORED
        // level, and an IR's authored level is a convolution gain, not a tap amplitude. Without it the
        // spring was +6.02 dB on a 96 kHz host and -0.74 on a 44.1 one, for the same Mix knob and the
        // same file (P68; the arithmetic and the measurements live on convolutionRateGain). Exactly 1
        // when the rates match, so the verbatim path stays verbatim to the bit.
        const float g = normalize_ ? normalizationGain (staged, outLen)
                      : resample   ? (float) convolutionRateGain (irSr, hostSr_)
                                   : 1.0f;
        for (auto& ch : staged) for (float& v : ch) v *= g;

        ir_.resize ((std::size_t) nch);
        for (int c = 0; c < nch; ++c) ir_[(std::size_t) c].swap (staged[(std::size_t) c]);
        normGain_   = g;
        // THE FLOOR IS ONLY THERE SO A ZERO CANNOT READ AS -inf, and it has to sit below every gain the
        // loader can produce. It was 1.0e-6f, which was below the old minimum (the normalization clamps at
        // -30 dB, i.e. 0.0316) and ABOVE the new one: the rate factor's smallest value is about 4.7e-10 —
        // under that the output is longer than INT_MAX and resampleIr refuses the load — so an IR file
        // claiming 0.024 Hz on a 48 kHz host applies 5.0e-7 and this reported -120.0000 dB for a gain that
        // is -126.0206. The linear accessor was right throughout; only the diagnostic lied.
        // ...and a NaN gets its own answer, because `std::max` cannot give it one: max(a, b) is
        // (a < b) ? b : a, every comparison against a NaN is false, so max(floor, NaN) is the FLOOR —
        // a finite, plausible dB reading for a gain that is not a number. (A NaN tap in the IR makes
        // the normalization gain NaN; the taps are then NaN too, which is the real problem, but the
        // diagnostic must not be the thing that hides it.) Reporting the NaN is the honest answer and
        // matches the linear accessor, which has always returned it.
        normGainDb_ = std::isfinite (g) ? 20.0f * std::log10 (std::max (1.0e-20f, g)) : g;

        // Publish the normalized taps to the convolver. On rejection (mid-crossfade) the retry reads back
        // from ir_ — which always holds the LATEST staged IR — so there is no dangling snapshot: a newer
        // load simply overwrites ir_ and the pending state, and flushPending applies the latest (coalescing).
        pendingNch_   = nch;
        pendingLen_   = outLen;
        pendingRetry_ = true;
        tryPublishPending();
    }

    // Hand the currently-staged IR (in ir_) to the convolver. Cheap on rejection: the convolver's state
    // check happens BEFORE any partition build, so a mid-crossfade reject does no work. Message thread.
    bool tryPublishPending()
    {
        if (! pendingRetry_) return true;
        if (! prepared_)     return false;   // convolver failed to arm (impossible with our clamps) — keep pending

        bool ok;
        if (pendingNch_ <= 1)
            ok = convolution_.setIr (ir_[0].data(), pendingLen_);                 // mono → broadcast (Stereo::yes)
        else
        {
            const float* banks[2] { ir_[0].data(), ir_[1].data() };
            ok = convolution_.setOperator (Conv::Topology::LRDiag, banks, 2, pendingLen_);   // true stereo
        }
        if (ok) pendingRetry_ = false;
        return ! pendingRetry_;
    }

    // 1 / (reference RMS gain) of an IR about to be staged: G² = Σ w(f)·P(f) / Σ w(f) over the positive-
    // frequency bins (DC excluded), where P(f) is the channel-mean power response |H(f)|² and
    // w(f) = 1 / (1 + (f/kIrRefShapeHz)²) is the one-pole-shaped reference's power spectrum.
    // Frequency domain (one real FFT per channel on the message thread) — equals convolving the
    // reference noise through the IR and reading the RMS ratio, without needing a signal.
    float normalizationGain (const std::vector<std::vector<float>>& taps, int len) const
    {
        // Analysis window: the first second. An IR's tail past that carries so little band energy
        // that it moves the reference gain by < 0.1 dB (measured on the factory set down to 10%
        // trims), and the trim drag re-runs this per mouse move — the cap keeps the drag light.
        const int cap = std::min (len, (int) std::lround (hostSr_));
        int N = 256; while (N < cap * 2 && N < (1 << 21)) N <<= 1;   // dense DTFT sampling of the IR
        // AN ANALYSIS WINDOW CANNOT BE LONGER THAN THE TRANSFORM THAT CARRIES IT. N stops doubling at
        // 1<<21, so above a megasample of window the two part company and the copy below ran off the end of
        // `padded`: measured under ASan as a 12 MB heap-buffer-overflow WRITE for a 3 000 000-sample IR at a
        // 4 MHz host — reachable through the public loadIR, and still reachable inside the 3 MHz rate
        // ceiling. The window is the first second OR the transform, whichever is shorter; past that the
        // truncation is the cap doing its job, not a defect.
        const int used = std::min (cap, N);
        felitronics::core::fft::DefaultRealFft fft;
        if (! fft.prepare (N)) return 1.0f;

        std::vector<float> padded ((std::size_t) N, 0.0f);
        std::vector<float> spec ((std::size_t) felitronics::core::fft::DefaultRealFft::spectrumFloats (N));
        std::vector<double> power ((std::size_t) (N / 2 + 1), 0.0);
        for (const auto& ch : taps)
        {
            std::fill (padded.begin(), padded.end(), 0.0f);
            std::copy (ch.begin(), ch.begin() + std::min<std::ptrdiff_t> (used, (std::ptrdiff_t) ch.size()),
                       padded.begin());
            fft.forward (padded.data(), spec.data());
            power[0]                    += (double) spec[0] * spec[0];                 // DC (unused below)
            power[(std::size_t) (N / 2)] += (double) spec[1] * spec[1];               // Nyquist
            for (int k = 1; k < N / 2; ++k)
                power[(std::size_t) k] += (double) spec[(std::size_t) (2 * k)] * spec[(std::size_t) (2 * k)]
                                        + (double) spec[(std::size_t) (2 * k + 1)] * spec[(std::size_t) (2 * k + 1)];
        }

        double num = 0.0, den = 0.0;
        const double chInv = 1.0 / (double) std::max<std::size_t> (1, taps.size());
        for (int k = 1; k <= N / 2; ++k)                             // DC excluded: not audio
        {
            const double f = (double) k * hostSr_ / N;
            const double w = 1.0 / (1.0 + (f / kIrRefShapeHz) * (f / kIrRefShapeHz));
            num += w * power[(std::size_t) k] * chInv;
            den += w;
        }
        const double gSq = den > 0.0 ? num / den : 0.0;
        if (gSq < std::pow (10.0, kIrRefFloorDb / 10.0))             // near-silent IR: leave it alone
            return 1.0f;
        const float gDb = (float) (-10.0 * std::log10 (gSq));
        return std::pow (10.0f, std::clamp (gDb, kIrNormMinDb, kIrNormMaxDb) / 20.0f);
    }

    double hostSr_   = 48000.0;
    int    channels_ = 2;
    int    maxBlock_ = 512;
    float  normGain_ = 1.0f;
    float  normGainDb_ = 0.0f;
    bool   normalize_ = true;              // false = reverb IR (skip reference-unity RMS normalization)
    bool   prepared_ = false;
    Conv   convolution_;                   // non-uniform (Gardner), block-independent, zero-latency

    std::vector<std::vector<float>> ir_;   // staging + retained latest taps (message thread)
    int  pendingNch_   = 0;                // channels of the staged IR (1 = mono broadcast, 2 = true stereo)
    int  pendingLen_   = 0;                // length of the staged IR
    bool pendingRetry_ = false;            // a staged IR still needs publishing (rejected mid-crossfade)
};

} // namespace felitronics::convolution
