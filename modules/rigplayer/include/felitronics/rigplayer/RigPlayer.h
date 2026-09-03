// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
#pragma once
// felitronics::rigplayer — ONE DEVICE OF A PACK, PLAYING. Hand it the pack's structures
// (namz::rig — the manifest as the canonical reader returns it) and a way to get a file's bytes, tell it
// where the knobs stand, feed it audio. It does, in signal order:
//
//     in → dry copy → pre tone (bands, curve) → chain trim → ┬→ trim A → model A → align ─┐
//                                                            └→ trim B → model B → align ─┴→ mix
//        → post tone (bands, curve) → dry/wet blend → out
//
//   • WHICH files play is namz::rig's policy (a turned control is law) joined to ModelBlend.h's pair
//     along the gain dial — RigSelection.h, pure.
//   • WHEN a model may be replaced, and how loudly each is heard, is felitronics::nam::BlendLaw —
//     run here once per block on the audio thread, the only writer of the weight. The message thread
//     posts a request and answers load asks; two one-way mailboxes, no locks. A load itself — bytes
//     fetched, a network built and warmed, some twenty milliseconds for a WaveNet — is a JOB the
//     player hands out (takeLoadJob) and the host runs where it likes (run, any thread) before
//     bringing it back (deliver): no thread of its own, no policy, no hiccup on the drawing thread.
//   • A slot that has stood silent under an unchanged request for kColdAfterSeconds goes COLD and its
//     model is not run — a dial parked on a capture costs one network, not two. The law owns the flag
//     (BlendState::cold) and wakes the slot, warm-up first, the moment the request changes; the model
//     stays where it is throughout. slotCold() and coldBlocks() read it out.
//   • A tone knob is applied as the pack describes it — a curve becomes a minimum-phase FIR, bands
//     become biquads — on the side of the models its circuit puts it. ToneKnobs.h, pure.
//   • A blend knob mixes the DI the models were fed, through the dry path's own response. BlendKnob.h.
//   • The slot trims are the files' own `input_db` — a linked setting plays its neighbour softer — and
//     the chain trim is the extension past the ends of the captured range.
//   • The models' offsets are measured from the bytes (AlignmentTable.h) by the host, on any thread,
//     and handed in; each slot's delay travels with its model and lands when the slot is silent.
//
// NOTHING OF THE LIBRARY IS KNOWN HERE — no take, no database, no file system, no JUCE. The bytes come
// through ModelSource by `files[].id`, when the law asks and never before: a swept dial fetches each
// capture once. This is the unit that leaves for the plugin as one directory.
//
// Threads, exactly as NamStage: prepare() / load() / the knob setters / service() / takeLoadJob() /
// deliver() on the message thread; run() on any thread but the audio one; process() alone on the audio
// thread — it allocates nothing, locks nothing, touches no file. The read-outs are atomics.

#include <felitronics/rigplayer/AlignmentTable.h>
#include <felitronics/rigplayer/BlendKnob.h>
#include <felitronics/rigplayer/RigSelection.h>
#include <felitronics/rigplayer/ToneKnobs.h>

#include <felitronics/convolution/CabConvolver.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/lineareq/MagnitudeCurve.h>
#include <felitronics/nam/BlendLaw.h>
#include <felitronics/nam/NamStage.h>
#include <namz_rig.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace felitronics::rigplayer {

class RigPlayer {
public:
    static constexpr int kMaxChannels = 2;
    // Biquads per side. Every band of every `sections` knob on that side is one; a pack carrying more
    // than this on one side would be a first, and the rest are left out (bandsDropped() says so).
    static constexpr int kMaxBands = 24;
    // The curve form's FIR: 1024 taps designed through 8192 points, the numbers the bench auditions with.
    static constexpr int kFirTaps = 1024, kFirDesign = 8192;
    static constexpr int kMaxDelay = felitronics::nam::kBlendMaxDelay;
    // How long a slot may stand at exactly zero, its request unchanged, before its model stops being
    // run (BlendLaw.h, "and one economy"). Long against a hand — a pause this long IS the dial at rest
    // — and short against a session, so the saving is there whenever nobody is turning. A host may set
    // its own (setColdAfterSeconds); zero or less = never.
    static constexpr double kColdAfterSeconds = 2.0;

    RigPlayer() = default;
    RigPlayer(const RigPlayer&) = delete;
    RigPlayer& operator=(const RigPlayer&) = delete;

    // ---------------------------------------------------------------- message thread: setup ----

    // Sizes every buffer and (re)builds anything designed for a rate. Never while process() runs.
    void prepare(double sampleRate, int maxBlock, int numChannels) {
        fs_       = sampleRate > 0.0 ? sampleRate : 48000.0;
        maxBlock_ = std::max(1, maxBlock);
        channels_ = std::clamp(numChannels, 1, kMaxChannels);
        coldAfter_.store(coldSamples(coldSeconds_), std::memory_order_release);
        for (auto& n : nam_) n.prepare(fs_, maxBlock_);
        // Not normalised: a tone curve's broadband level is part of what the pack says, and a dry path's
        // level rides in its gain. Reference-unity RMS is for cabinets, which these are not.
        for (auto& f : fir_) f.prepare(fs_, maxBlock_, channels_, 0.1, false);
        dry_.prepare(fs_, maxBlock_, channels_, 0.1, false);
        for (int c = 0; c < kMaxChannels; ++c) {
            slotB_[c].assign((std::size_t) maxBlock_, 0.0f);
            dryBuf_[c].assign((std::size_t) maxBlock_, 0.0f);
            spare_[c].assign((std::size_t) maxBlock_, 0.0f);
            for (auto& t : lagTail_) t[c].fill(0.0f);
        }
        for (auto& side : bq_) for (auto& band : side) for (auto& b : band) b.reset();
        curChain_ = chainGain_.load(std::memory_order_relaxed);
        for (int i = 0; i < 2; ++i) curSlot_[i] = slotGain_[i].load(std::memory_order_relaxed);
        curDry_ = dryGain_.load(std::memory_order_relaxed);
        curWet_ = wetGain_.load(std::memory_order_relaxed);
        prepared_ = true;
        if (loaded_) {                                   // the FIRs and the bands were designed for the old rate
            for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
            rebuildDry();
        }
    }
    bool   prepared()   const { return prepared_; }
    double sampleRate() const { return fs_; }
    int    channels()   const { return channels_; }

    // The first stage of the chain this player can run. False = nothing to play (no NAM stage, or one
    // with no files) and the player is left empty.
    bool load(const namz::rig::Rig& rig, ModelSource source) {
        for (const auto& st : rig.chain)
            if (st.kind == namz::rig::StageKind::Nam) return load(st, std::move(source));
        unload();
        return false;
    }

    bool load(const namz::rig::Stage& stage, ModelSource source) {
        unload();
        if (stage.kind != namz::rig::StageKind::Nam || stage.device.files.empty()) return false;
        stage_  = stage;
        source_ = std::move(source);
        const auto& d = stage_.device;
        // One model per DISTINCT file: several settings pointing at one file — the pack's link — share
        // its bytes, its id and its delay; only their trims differ, and those belong to the slot.
        fileModel_.assign(d.files.size(), -1);
        for (std::size_t i = 0; i < d.files.size(); ++i) {
            const auto& id = d.files[i].id;
            int m = -1;
            for (std::size_t k = 0; k < models_.size(); ++k) if (models_[k].id == id) { m = (int) k; break; }
            if (m < 0) { models_.push_back({ id, {} }); m = (int) models_.size() - 1; }
            fileModel_[i] = m;
        }
        // Where every knob starts: the pack says, the player never invents (namz::rig::defaultSettings,
        // and `default`-else-`reference` for the linear knobs).
        axes_ = namz::rig::defaultSettings(d);
        if (const auto* dial = crossfadeDial(d)) {
            dial_ = dial->name;
            dialSweep_ = dial->sweep;
            const int deg = degreesOf(axes_[dial_]);
            dialDeg_ = deg >= 0 ? (double) deg : 0.0;
        }
        // …and the defaults are per CONTROL, not a promise that the combination was captured. If no
        // file stands there, pin the first control at its default and let resolve() find the closest
        // captured combination — the pack opens on something sounding, never on a silent panel.
        if (select(d, axes_, dial_, dialDeg_, shape_, topExtend_).fileA < 0)
            for (const auto& c : d.controls)
                if (namz::rig::resolve(d, axes_, c.name, axes_[c.name]) != nullptr) break;
        for (const auto& t : tones())      toneAt_[t.name]  = toneStart(t);
        for (const auto& b : stage_.blend) blendAt_[b.name] = blendStart(b);
        // The models' offsets, from the pack when it carries them — measured once, when it was built.
        // A pack without them leaves the table empty; the host may measure and setAlignment().
        align_ = AlignmentTable::fromDevice(d);
        loaded_ = true;
        if (prepared_) {
            for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
            rebuildDry();
        }
        apply();
        return true;
    }

    // Forget the device: both models cleared, the law told to start over, every filter a wire.
    void unload() {
        loaded_ = false;
        ++gen_;                                          // a job out for this pack comes back to be dropped
        setRequest(0, 0, 0.0f);
        // A landing (or a failure) published and not yet consumed dies with its pack — and it must die
        // BEFORE the forget is posted: consumed after it, it would land a stale model in a law that was
        // just wiped and a stage that was just emptied, and an empty stage carried is the raw DI.
        landFlag_.store(false, std::memory_order_release);
        landFail_.store(0, std::memory_order_release);
        forget_.store(true, std::memory_order_release);
        loadAsk_.store(0, std::memory_order_release);    // …and an ask for it is void with it
        for (auto& n : nam_) n.clearModel();
        models_.clear(); fileModel_.clear();
        axes_.clear(); dial_.clear(); dialSweep_ = 0; dialDeg_ = 0.0;
        toneAt_.clear(); blendAt_.clear(); toneOverride_.clear();
        sel_ = {}; plan_ = {};
        stage_ = {};
        source_ = nullptr;
        align_ = {};
        for (int s = 0; s < 2; ++s) {
            curveSum_[s].clear();
            firBypass_[s].store(true, std::memory_order_release);
            publishBands(s, BandSet {});
        }
        bandsDropped_ = 0;
        dryActive_.store(false, std::memory_order_release);
        dryFirBypass_.store(true, std::memory_order_release);
        dryGain_.store(0.0f, std::memory_order_release);
        wetGain_.store(1.0f, std::memory_order_release);
        chainGain_.store(1.0f, std::memory_order_release);
        for (auto& g : slotGain_) g.store(1.0f, std::memory_order_release);
        for (auto& p : pendDelay_) p.store(0, std::memory_order_release);
        for (auto& c : slotCold_) c.store(false, std::memory_order_release);
    }
    bool loaded() const { return loaded_; }
    const namz::rig::Stage& stage() const { return stage_; }

    // ---------------------------------------------------------------- message thread: the panel ----

    // A dial, in degrees of its rotation. The crossfade dial moves continuously; any other captured dial
    // selects its nearest captured position; a tone or blend dial moves its own filter. False = no such
    // dial on this device, or a value it cannot take.
    bool setDial(const std::string& control, double degrees) {
        if (! loaded_) return false;
        if (! dial_.empty() && control == dial_) {
            dialDeg_ = std::clamp(degrees, 0.0, (double) dialSweep_);
            apply();
            return true;
        }
        for (const auto& c : stage_.device.controls)
            if (c.name == control) {
                const auto v = nearestValue(c, degrees);
                return ! v.empty() && setSwitch(control, v);
            }
        for (const auto& t : tones())
            if (t.name == control && t.sweep > 0) {
                toneAt_[control] = degreeValue(degrees, t.sweep);
                if (prepared_) rebuildKnob(t);
                return true;
            }
        for (const auto& b : stage_.blend)
            if (b.name == control && b.sweep > 0) {
                blendAt_[control] = degreeValue(degrees, b.sweep);
                applyBlendGains();
                return true;
            }
        return false;
    }

    // A control set to one of its values. For a captured axis this is namz::rig::resolve: the turned
    // control is pinned, everything else stays where the hand left it, and the closest captured
    // combination is chosen — the crossfade dial keeps its angle and follows along the new knots. False
    // = no such control, or a value nothing was captured at (then nothing changes).
    bool setSwitch(const std::string& control, const std::string& value) {
        if (! loaded_) return false;
        for (const auto& c : stage_.device.controls)
            if (c.name == control) {
                if (! dial_.empty() && control == dial_) {
                    const int deg = degreesOf(value);
                    return deg >= 0 && setDial(control, (double) deg);
                }
                if (namz::rig::resolve(stage_.device, axes_, control, value) == nullptr) return false;
                apply();
                return true;
            }
        for (const auto& t : tones())
            if (t.name == control) {
                // A SWITCH ONLY HAS THE POSITIONS IT DECLARES. Accepting any word and storing it read
                // back as a setting the knob does not have, while the sound was the reference — the
                // caller was told `true` and heard nothing. A dial keeps taking any degree: every angle
                // of its travel is playable, swept or not.
                if (t.sweep <= 0 && ! t.positions.empty()) {
                    bool found = false;
                    for (const auto& p : t.positions) if (p.value == value) { found = true; break; }
                    if (! found) return false;
                }
                toneAt_[control] = value;
                if (prepared_) rebuildKnob(t);
                return true;
            }
        for (const auto& b : stage_.blend)
            if (b.name == control) {
                blendAt_[control] = value;
                applyBlendGains();
                return true;
            }
        return false;
    }

    // What the panel reads back. `settings()` is the captured combination — for the crossfade dial it
    // names the knot at or below the hand, the angle itself is `dialDegrees()`.
    const Settings&    settings()    const { return axes_; }
    const std::string& dialName()    const { return dial_; }
    int                dialSweep()   const { return dialSweep_; }
    double             dialDegrees() const { return dialDeg_; }
    // A tone or blend knob's position as the pack spells it ("150", or a switch's token). Empty = no
    // such knob.
    std::string knobValue(const std::string& control) const {
        if (const auto t = toneAt_.find(control); t != toneAt_.end()) return t->second;
        if (const auto b = blendAt_.find(control); b != blendAt_.end()) return b->second;
        if (const auto a = axes_.find(control); a != axes_.end()) return a->second;
        return {};
    }
    const Selection& selection() const { return sel_; }
    const SlotPlan&  plan()      const { return plan_; }

    // TONE HANDED IN BESIDE THE MANIFEST — the same structures the pack carries, from another source: a
    // bench fitting bands by ear hands them here and hears them at once, then packs and hears them from
    // the file. ONE path applies tone; only where the knob's description came from differs. A knob named
    // here replaces the pack's block of that name (its `positions` or `sections`, whichever it had); a
    // name the pack has no block for is a new knob, starting at its default. The pack's own blocks are
    // untouched and return with clearToneOverride().
    void setToneOverride(std::vector<namz::rig::Tone> tone) {
        toneOverride_ = std::move(tone);
        for (const auto& t : tones()) if (! toneAt_.count(t.name)) toneAt_[t.name] = toneStart(t);
        if (prepared_ && loaded_) for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
    }
    void clearToneOverride() {
        if (toneOverride_.empty()) return;
        toneOverride_.clear();
        if (prepared_ && loaded_) for (int s = 0; s < 2; ++s) { rebuildCurves(s); rebuildBands(s); }
    }
    bool toneOverridden() const { return ! toneOverride_.empty(); }
    // The tone knobs as they PLAY: the pack's, with any handed-in block of the same name in its place,
    // and the handed-in knobs the pack does not have after them.
    std::vector<namz::rig::Tone> tones() const {
        std::vector<namz::rig::Tone> out;
        for (const auto& t : stage_.tone) {
            const namz::rig::Tone* use = &t;
            for (const auto& o : toneOverride_) if (o.name == t.name) { use = &o; break; }
            out.push_back(*use);
        }
        for (const auto& o : toneOverride_)
            if (std::none_of(stage_.tone.begin(), stage_.tone.end(), [&o](const namz::rig::Tone& t) { return t.name == o.name; }))
                out.push_back(o);
        return out;
    }

    // THE SHAPE OF THE HANDOVER — where the 50/50 lands between two captures and how wide the fade is.
    // The default is the original law (midpoint, full span); the bench turns these, a plugin need not.
    void setBlendShape(BlendShape s) { shape_ = s; if (loaded_) apply(); }
    BlendShape blendShape() const { return shape_; }
    // How much harder the top capture is driven per degree past it (ModelBlend.h explains the guess).
    void setTopExtendDbPerDeg(double v) { topExtend_ = v; if (loaded_) apply(); }
    double topExtendDbPerDeg() const { return topExtend_; }

    // Apply each model's loudness tag. Off, the models play raw, as the hardware returned them.
    void setNormalize(bool on) { normalize_.store(on, std::memory_order_release); }
    bool normalize() const { return normalize_.load(std::memory_order_acquire); }

    // Apply the pack's per-file input trims (`input_db` — a linked setting plays its neighbour
    // softer). OFF feeds every capture at unity instead: for a library shot at one honest level
    // the stated attenuations are somebody else's story, and into a nonlinear model a few dB
    // less in is a lot less out. Read on the audio side, so the toggle lands on the next block.
    void setInputTrims(bool on) { inputTrims_.store(on, std::memory_order_release); }
    bool inputTrims() const { return inputTrims_.load(std::memory_order_acquire); }

    // The rest after which a silent slot's model is no longer run (kColdAfterSeconds). A slot that
    // sounds at all is never cold; between two captures both models run, on one they do not. Zero or
    // less = never. Takes effect from the next block; a slot already asleep stays asleep until woken.
    void   setColdAfterSeconds(double seconds) {
        coldSeconds_ = seconds;
        coldAfter_.store(coldSamples(seconds), std::memory_order_release);
    }
    double coldAfterSeconds() const { return coldSeconds_; }

    // The models' measured offsets (AlignmentTable.h), for a pack that does not carry its own. A delay
    // travels with a model and is applied at the one instant its slot is silent — weight exactly zero
    // — never as a splice on a live signal. So a table handed in BEFORE playing lands with the first
    // loads; one that arrives mid-mix waits until the dial visits a knot. A host that has the bytes in
    // hand should measure before it plays.
    void setAlignment(AlignmentTable table) {
        align_ = std::move(table);
        if (loaded_) stageDelays();
    }
    const AlignmentTable& alignment() const { return align_; }
    bool alignmentFromPack() const { return align_.fromPack && ! align_.empty(); }

    // Regular message-thread housekeeping: frees models the audio thread has stepped past, retries a
    // filter the convolver rejected mid-fade. Call from a timer, a few times a second at least — and
    // ask for the load job after it (takeLoadJob), because the law waits for that.
    void service() {
        for (auto& n : nam_) n.collectGarbage();
        for (auto& f : fir_) f.flushPending();
        dry_.flushPending();
    }

    // ---------------------------------------------------------------- message thread: loading ----

    // A MODEL THE LAW ASKED FOR, as work for the host to run where it likes — a worker thread, or right
    // here. The job carries everything the work needs (the source, the bytes if already fetched, the
    // numbers the stages were prepared with) and nothing of the player, so run() is a pure function of
    // it. At most ONE job is out at a time, and a job taken MUST come back through deliver() — with a
    // null model if it failed — or the law waits for it forever.
    struct LoadJob {
        int           slot = 0;
        std::uint64_t model = 0;                                   // BlendModelId: models_ index + 1
        std::string   fileId;                                      // the pack's `files[].id`
        std::shared_ptr<const std::vector<std::byte>> bytes;       // null = not fetched yet
        ModelSource   source;
        double        sampleRate = 48000.0;
        int           maxBlock = 512;
        std::uint64_t generation = 0;                              // the load() it belongs to
    };
    struct Loaded {
        LoadJob job;
        std::shared_ptr<const std::vector<std::byte>> bytes;       // the job's own, or what run() fetched
        bool    fetched = false;                                   // …and whether it had to
        felitronics::nam::NamStage::PreparedModel model;           // null = the bytes were not a model
    };

    std::optional<LoadJob> takeLoadJob() {
        if (jobOut_ || ! loaded_) return std::nullopt;
        const auto ask = loadAsk_.load(std::memory_order_acquire);
        if (ask == 0) return std::nullopt;
        const auto id = ask >> 1;
        if (id == 0 || id > models_.size()) {                      // an ask that outlived its pack
            loadAsk_.store(0, std::memory_order_release);
            landFail_.store((int) (ask & 1u) + 1, std::memory_order_release);
            return std::nullopt;
        }
        const auto& m = models_[(std::size_t) id - 1];
        LoadJob job;
        job.slot = (int) (ask & 1u); job.model = id; job.fileId = m.id; job.bytes = m.bytes;
        job.source = source_; job.sampleRate = fs_; job.maxBlock = maxBlock_; job.generation = gen_;
        jobOut_ = true;
        return job;
    }

    // Any thread but the audio one. The bytes through the source when the job has none, then the heavy
    // half of a load (NamStage::prepareModel). Touches nothing of any player.
    static Loaded run(LoadJob job) {
        Loaded out;
        out.bytes = job.bytes;
        if (out.bytes == nullptr && job.source) {
            out.bytes = std::make_shared<const std::vector<std::byte>>(job.source(job.fileId));
            out.fetched = true;
        }
        if (out.bytes != nullptr && ! out.bytes->empty())
            out.model = felitronics::nam::NamStage::prepareModel(out.bytes->data(), out.bytes->size(),
                                                                 job.sampleRate, job.maxBlock);
        out.job = std::move(job);
        return out;
    }

    // Message thread. The light half: the bytes kept (fetched once per file, however many times the
    // dial crosses it), the model installed in its slot, the law told it landed — or failed. A job from
    // a pack that is gone is dropped whole: the law was told to start over when the pack left.
    void deliver(Loaded loaded) {
        jobOut_ = false;
        auto& job = loaded.job;
        if (job.generation != gen_ || ! loaded_) return;
        if (loaded.fetched && job.model >= 1 && job.model <= models_.size()) {
            auto& m = models_[(std::size_t) job.model - 1];
            if (m.bytes == nullptr) { m.bytes = loaded.bytes; ++loads_; }
        }
        loadSlot(job.slot, std::move(loaded.model), job.model);
    }

    // service() with every job due run right here, on this thread — for a host without a worker, and
    // for a test. The same path on one thread, not a second one.
    void serviceHere() {
        service();
        while (auto job = takeLoadJob()) deliver(run(std::move(*job)));
    }

    // Host-rate latency: the models' own rate-matching (the filters are minimum-phase, the alignment
    // delays are relative and the reference model carries none).
    int latencySamples() const {
        return std::max(nam_[0].latencySamples(), nam_[1].latencySamples());
    }

    // ---- read-outs of the audio thread, for a strip or a dump ----
    float liveMix() const { return liveMix_.load(std::memory_order_acquire); }        // applied weight of slot 1
    // The file id held in a slot right now ("" = nothing yet).
    std::string heldFileId(int slot) const {
        const auto id = held_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed);
        return id == 0 || id > models_.size() ? std::string() : models_[(std::size_t) id - 1].id;
    }
    int appliedSlotDelay(int slot) const { return slotDelay_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    // The loudness tag of what slot 0 holds — during a crossfade each model carries its own, and slot
    // 0 is the one a panel names. Message thread.
    double modelLoudness()   const { return nam_[0].modelLoudness(); }
    bool   modelHasLoudness() const { return nam_[0].modelHasLoudness(); }
    // Instrumentation, for a bench's dump: blocks in which a slot was still warming, weight moves of
    // more than 2 % inside one block and the biggest of them — a whole swing inside one block is a step
    // in all but name. Counted since the last clearCounters().
    int   warmBlocks()  const { return warmBlocks_.load(std::memory_order_relaxed); }
    int   mixJumps()    const { return mixJumps_.load(std::memory_order_relaxed); }
    float biggestJump() const { return biggestJump_.load(std::memory_order_relaxed); }
    // Whether a slot is COLD right now — held, silent, its model not run (kColdAfterSeconds) — and the
    // blocks it has slept since the last clearCounters(), for the dump beside warmBlocks().
    bool  slotCold(int slot)   const { return slotCold_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    int   coldBlocks(int slot) const { return coldBlocks_[(std::size_t) (slot & 1)].load(std::memory_order_relaxed); }
    void  clearCounters() {
        warmBlocks_.store(0, std::memory_order_relaxed);
        mixJumps_.store(0, std::memory_order_relaxed);
        biggestJump_.store(0.0f, std::memory_order_relaxed);
        for (auto& c : coldBlocks_) c.store(0, std::memory_order_relaxed);
    }
    long long modelLoads() const { return loads_; }                                   // how many the source was asked for
    int bandsDropped() const { return bandsDropped_; }                                // sections past kMaxBands

    // The tone as it stands, for drawing: the summed curve of the curve-form knobs on a side (dB on
    // `commonGrid()`, empty = nothing on that side) and the bands of the section-form knobs.
    const std::vector<double>& commonGrid() const { return commonGrid_; }
    const std::vector<double>& curveDb(int side) const { return curveSum_[side & 1]; }
    bool curveActive(int side) const { return ! firBypass_[side & 1].load(std::memory_order_acquire); }   // …and whether it is a FIR or a wire
    const std::vector<felitronics::rigplayer::SectionBiquad>& bands(int side) const { return bandsShown_[side & 1]; }
    BlendGains blendGains() const { return { (double) dryGain_.load(std::memory_order_relaxed),
                                             (double) wetGain_.load(std::memory_order_relaxed) }; }

    // ---------------------------------------------------------------------- the audio thread ----

    // In place. `numChannels` up to what prepare() was given (fewer is fine: the missing planes are
    // silent). Any `numSamples`; longer than maxBlock is walked in pieces.
    //
    // THE MODELS RUN ON THE PLANES THEY ARE GIVEN, not on the width the player was prepared for.
    // A host prepared for a stereo bus that plays one plane — a mono chain — used to push the silent
    // second plane through both networks as well: four WaveNet passes a block for one channel of
    // sound, and a quarter of a core gone to nothing. The convolvers still take the prepared width
    // (the partitioned convolution no-ops when handed fewer planes), and the spare planes stay
    // zeroed for them; only the expensive part — the two models — shrinks to what is playing.
    void process(float* const* io, int numChannels, int numSamples) {
        if (! prepared_ || io == nullptr || numSamples <= 0) return;
        const int nch = std::clamp(numChannels, 0, channels_);
        if (nch == 0) return;
        const bool norm = normalize_.load(std::memory_order_acquire);
        float* a[kMaxChannels]; float* b[kMaxChannels]; float* d[kMaxChannels];
        for (int c = 0; c < kMaxChannels; ++c) { b[c] = slotB_[c].data(); d[c] = dryBuf_[c].data(); }

        int done = 0;
        while (done < numSamples) {
            const int count = std::min(maxBlock_, numSamples - done);
            for (int c = 0; c < channels_; ++c) {
                if (c < nch && io[c] != nullptr) a[c] = io[c] + done;
                else { a[c] = spare_[c].data(); std::fill(a[c], a[c] + count, 0.0f); }
            }
            for (int s = 0; s < 2; ++s) takeBands(s);

            // Keep the DRY block before anything touches it: the same DI the models are fed, which is
            // the whole economy of a blend. Taken ahead of the pre-model tone — the hardware's dry path
            // leaves at the input jack.
            const bool mixDry = dryActive_.load(std::memory_order_acquire);
            if (mixDry)
                for (int c = 0; c < channels_; ++c) std::copy(a[c], a[c] + count, d[c]);

            // The knobs that sit BEFORE the distortion in the hardware, shaping what gets distorted.
            runBands(0, a, count);
            if (! firBypass_[0].load(std::memory_order_acquire)) fir_[0].process(a, channels_, count);

            // The chain's trim: past the top capture it rises with the angle, and a fast hand moves
            // tens of degrees between two events — so it is smoothed like every other gain here.
            rampInto(a, nch, count, chainGain_.load(std::memory_order_acquire), curChain_);

            // THE LAW, once per block, and the only writer of the weight. Everything it needs arrives
            // through atomics; everything it asks for leaves the same way.
            if (forget_.exchange(false, std::memory_order_acquire)) { blend_ = {}; coldRt_[0] = coldRt_[1] = false; }
            if (const int f = landFail_.exchange(0, std::memory_order_acquire); f != 0)
                felitronics::nam::blendLoadFailed(blend_, f - 1);
            if (landFlag_.exchange(false, std::memory_order_acquire)) {
                const int ls = landSlot_.load(std::memory_order_relaxed);
                const int nd = landDelay_.load(std::memory_order_relaxed);
                slotDelay_[(std::size_t) (ls & 1)].store(nd, std::memory_order_relaxed);
                pendDelay_[(std::size_t) (ls & 1)].store(nd, std::memory_order_release);
                for (auto& t : lagTail_[(std::size_t) (ls & 1)]) t.fill(0.0f);
                felitronics::nam::blendLanded(blend_, ls, landModel_.load(std::memory_order_relaxed),
                                              landWarm_.load(std::memory_order_relaxed));
            }
            // ACQUIRE FIRST: the writer publishes the two ids and then releases on the target.
            felitronics::nam::BlendRequest req;
            req.targetB = (double) reqTarget_.load(std::memory_order_acquire);
            req.want[0] = reqWant_[0].load(std::memory_order_relaxed);
            req.want[1] = reqWant_[1].load(std::memory_order_relaxed);
            felitronics::nam::BlendPolicy policy;
            policy.coldAfterSamples = coldAfter_.load(std::memory_order_relaxed);
            const auto law = felitronics::nam::blendStep(blend_, req, count, policy);
            if (law.load.wanted && loadAsk_.load(std::memory_order_relaxed) == 0)
                loadAsk_.store((law.load.model << 1) | (std::uint64_t) (law.load.slot & 1), std::memory_order_release);
            for (int i = 0; i < 2; ++i) {
                held_[(std::size_t) i].store(blend_.held[i], std::memory_order_relaxed);
                slotCold_[(std::size_t) i].store(blend_.cold[i], std::memory_order_relaxed);
                // FALLING ASLEEP CLEARS THE DELAY LINE, exactly as a landing does (above): the tail would
                // otherwise hold the slot's last samples from before the rest, and a model that declares
                // no field (need == 0) is heard on the very block it wakes in — with those samples first.
                if (blend_.cold[i] && ! coldRt_[i])
                    for (auto& t : lagTail_[(std::size_t) i]) t.fill(0.0f);
                coldRt_[i] = blend_.cold[i];
                if (blend_.cold[i]) coldBlocks_[(std::size_t) i].fetch_add(1, std::memory_order_relaxed);
            }
            // …and a staged retime lands the same way a model does: only where the slot is silent.
            for (int i = 0; i < 2; ++i) {
                const double w = i == 0 ? 1.0 - law.endB : law.endB;
                const int want = pendDelay_[(std::size_t) i].load(std::memory_order_acquire);
                if (w <= 0.0 && want != slotDelay_[(std::size_t) i].load(std::memory_order_relaxed)) {   // exactly zero: the law clamps to [0, 1]
                    slotDelay_[(std::size_t) i].store(want, std::memory_order_relaxed);
                    for (auto& t : lagTail_[(std::size_t) i]) t.fill(0.0f);
                }
            }
            liveMix_.store((float) law.endB, std::memory_order_release);
            if (blend_.fed[0] < blend_.need[0] || blend_.fed[1] < blend_.need[1])
                warmBlocks_.fetch_add(1, std::memory_order_relaxed);
            if (const float jump = (float) std::abs(law.endB - law.beginB); jump > 0.02f) {
                mixJumps_.fetch_add(1, std::memory_order_relaxed);
                float prev = biggestJump_.load(std::memory_order_relaxed);
                while (jump > prev && ! biggestJump_.compare_exchange_weak(prev, jump, std::memory_order_relaxed)) {}
            }

            // THE SECOND CAPTURE'S COPY is taken here, after the pre-model tone and the chain trim, so
            // both models are fed the identical signal — and only now do the two part company, each
            // through its own trim: a linked setting plays its neighbour's model with less going in.
            // A COLD SLOT IS NOT RUN — not copied into, not trimmed, not modelled, not delayed, not even
            // mixed: the law holds its weight at exactly zero for as long as it sleeps, so the other
            // slot IS the sound, and its buffer — whatever it holds, however old — is left out rather
            // than multiplied by zero (a NaN times zero is a NaN, for ever). Its trim ramp resumes from
            // where it stopped and settles within its block, like any other gain change; its delay line
            // was cleared on the way to sleep. (The rail is checked as well as the flag: the flag says
            // what the law decided, the weights say what this block sounds like.)
            const bool run[2] { ! blend_.cold[0], ! blend_.cold[1] };
            const bool aAlone = ! run[1] && law.beginB <= 0.0 && law.endB <= 0.0;
            const bool bAlone = ! run[0] && law.beginB >= 1.0 && law.endB >= 1.0;
            const bool trims = inputTrims_.load(std::memory_order_acquire);
            if (run[1]) for (int c = 0; c < nch; ++c) std::copy(a[c], a[c] + count, b[c]);
            if (run[0]) rampInto(a, nch, count, trims ? slotGain_[0].load(std::memory_order_acquire) : 1.0f, curSlot_[0]);
            if (run[1]) rampInto(b, nch, count, trims ? slotGain_[1].load(std::memory_order_acquire) : 1.0f, curSlot_[1]);
            if (run[0]) nam_[0].process(a, nch, count, norm);
            if (run[1]) nam_[1].process(b, nch, count, norm);
            // Align BEFORE the weights: during a ramp the two gains must sum to one at the SAME instant.
            // Each slot carries its own history per channel, advanced every block including at zero.
            for (int c = 0; c < nch; ++c) {
                if (run[0]) felitronics::nam::blendDelay(a[c], lagTail_[0][(std::size_t) c].data(), kMaxDelay,
                                                         slotDelay_[0].load(std::memory_order_acquire), count);
                if (run[1]) felitronics::nam::blendDelay(b[c], lagTail_[1][(std::size_t) c].data(), kMaxDelay,
                                                         slotDelay_[1].load(std::memory_order_acquire), count);
                if (aAlone)      { /* slot 1 asleep at zero: `a` is the whole sound */ }
                else if (bAlone) std::copy(b[c], b[c] + count, a[c]);
                else             felitronics::nam::blendMix(a[c], b[c], count, law.beginB, law.endB);
                // …and the law's own gain, which is below one only while nothing has been fed yet.
                if (law.beginGain < 1.0 || law.endGain < 1.0) {
                    const float dg = (float) (law.endGain - law.beginGain) / (float) count;
                    float g = (float) law.beginGain;
                    for (int i = 0; i < count; ++i) { g += dg; a[c][i] *= g; }
                }
            }

            // The knobs AFTER the distortion, shaping what came out.
            runBands(1, a, count);
            if (! firBypass_[1].load(std::memory_order_acquire)) fir_[1].process(a, channels_, count);

            // The blend, as the hardware sums it: the dry path through its own response, then the two
            // gains the pack states for this position.
            if (mixDry) {
                if (! dryFirBypass_.load(std::memory_order_acquire)) dry_.process(d, channels_, count);
                const float wantDry = dryGain_.load(std::memory_order_acquire);
                const float wantWet = wetGain_.load(std::memory_order_acquire);
                const float decay   = std::exp(-(float) count / (float) (0.010 * fs_));
                const float endDry  = wantDry + (curDry_ - wantDry) * decay;
                const float endWet  = wantWet + (curWet_ - wantWet) * decay;
                const float stepDry = (endDry - curDry_) / (float) count;
                const float stepWet = (endWet - curWet_) / (float) count;
                for (int c = 0; c < channels_; ++c) {
                    float gd = curDry_, gw = curWet_;
                    for (int i = 0; i < count; ++i) {
                        gd += stepDry; gw += stepWet;
                        a[c][i] = a[c][i] * gw + d[c][i] * gd;
                    }
                }
                curDry_ = endDry; curWet_ = endWet;
            }
            done += count;
        }
    }

private:
    // ------------------------------------------------------------------ the decision, applied ----

    static std::string degreeValue(double degrees, int sweep) {
        return std::to_string((int) std::lround(std::clamp(degrees, 0.0, (double) sweep)));
    }

    static float linOf(double db) { return std::abs(db) < 0.0005 ? 1.0f : (float) std::pow(10.0, db / 20.0); }

    std::uint64_t modelIdOf(int file) const {
        return file < 0 || (std::size_t) file >= fileModel_.size() ? 0 : (std::uint64_t) fileModel_[(std::size_t) file] + 1;
    }

    const std::string& fileIdOfModel(std::uint64_t id) const {
        static const std::string none;
        return id == 0 || id > models_.size() ? none : models_[(std::size_t) id - 1].id;
    }

    // In HOST samples. A pack's numbers are in the model's own rate, which the stage knows once the
    // model is loaded; a measured table carries the rate it was measured at.
    int delayOfModel(std::uint64_t id, const felitronics::nam::NamStage& st) const {
        return align_.delayOf(fileIdOfModel(id), fs_, st.modelSampleRate());
    }

    // The whole decision for the panel as it stands, handed to the audio thread as a request. What the
    // law does with it — in what order, and whether a model may be replaced yet — is its business.
    void apply() {
        const auto& d = stage_.device;
        sel_  = select(d, axes_, dial_, dialDeg_, shape_, topExtend_);
        plan_ = slotPlan(sel_, d);
        if (plan_.file[0] < 0) { setRequest(0, 0, 0.0f); return; }
        setRequest(modelIdOf(plan_.file[0]), modelIdOf(plan_.file[1]), (float) plan_.targetB);
        for (int i = 0; i < 2; ++i) slotGain_[i].store(linOf(plan_.inputDb[i]), std::memory_order_release);
        chainGain_.store(linOf(sel_.extendDb), std::memory_order_release);
        stageDelays();
        // The dial's own entry names the knot at (or below) the hand, so a later resolve() on another
        // control scores "keep what the user had" against the capture actually sounding.
        if (! dial_.empty() && sel_.fileA >= 0)
            if (const auto v = d.files[(std::size_t) sel_.fileA].settings.find(dial_);
                v != d.files[(std::size_t) sel_.fileA].settings.end())
                axes_[dial_] = v->second;
    }

    // For the slots as loaded now; a model still to come brings its own delay when it lands.
    void stageDelays() {
        for (int i = 0; i < 2; ++i)
            pendDelay_[(std::size_t) i].store(delayOfModel(modelIdOf(plan_.file[i]), nam_[i]), std::memory_order_release);
    }

    void setRequest(std::uint64_t want0, std::uint64_t want1, float targetB) {
        reqWant_[0].store(want0, std::memory_order_relaxed);
        reqWant_[1].store(want1, std::memory_order_relaxed);
        reqTarget_.store(std::clamp(targetB, 0.0f, 1.0f), std::memory_order_release);
    }

    // Put this model in that slot and tell the law it landed, with the delay that travels with the
    // model — applied when it lands, the one instant the slot is guaranteed silent. Only the ask this
    // answers is cleared, and a failure carries its own slot: sharing either cell stranded the other
    // slot for the rest of a session.
    bool loadSlot(int slot, felitronics::nam::NamStage::PreparedModel prepared, std::uint64_t model) {
        auto& stage = nam_[(std::size_t) (slot & 1)];
        const bool ok = prepared != nullptr && stage.install(std::move(prepared));
        const int delay = ok ? delayOfModel(model, stage) : 0;    // after the load: the rate is the model's
        if (const auto ask = loadAsk_.load(std::memory_order_acquire); ask != 0 && (int) (ask & 1u) == slot)
            loadAsk_.store(0, std::memory_order_release);
        if (! ok) { landFail_.store(slot + 1, std::memory_order_release); return false; }
        landSlot_.store(slot, std::memory_order_relaxed);
        landModel_.store(model, std::memory_order_relaxed);
        landWarm_.store(warmFor(stage), std::memory_order_relaxed);
        landDelay_.store(std::clamp(delay, 0, kMaxDelay), std::memory_order_relaxed);
        landFlag_.store(true, std::memory_order_release);
        return true;
    }

    // How many samples this model owes before it may be heard, in THIS rate: its receptive field,
    // scaled — a 96 kHz host feeds twice as many to fill the same network — plus one block, so a slot
    // is never marked audible for a block it is still short in.
    long long warmFor(const felitronics::nam::NamStage& st) const {
        const int pre = st.prewarmSamples();
        if (pre <= 0) return 0;
        const double mr = st.modelSampleRate();
        const double scale = (mr > 0.0 && fs_ > 0.0) ? fs_ / mr : 1.0;
        return (long long) std::ceil((double) pre * scale) + maxBlock_;
    }

    // The rest before a slot goes cold, in this rate's samples, for the law. Zero = never.
    long long coldSamples(double seconds) const {
        return seconds > 0.0 ? (long long) std::llround(seconds * fs_) : 0;
    }

    // ------------------------------------------------------------------------------- the tone ----

    // Which side of the models a knob sits on: 0 before, 1 after, -1 for a placement this player does
    // not know — such a knob is skipped, as the format says.
    static int sideOf(const namz::rig::Tone& t) {
        return t.placement == "pre" ? 0 : t.placement == "post" ? 1 : -1;
    }

    void rebuildKnob(const namz::rig::Tone& t) {
        const int side = sideOf(t);
        if (side < 0) return;
        // Bands, whether they travel with a dial or stand still at a switch's position; else the curve.
        if (! t.sections.empty() || bandsPerPosition(t)) rebuildBands(side); else rebuildCurves(side);
    }

    // Every curve-form knob on one side, summed (cascaded linear filters multiply, which in dB is a
    // sum) on one common grid — each knob may ship its own — and designed as ONE minimum-phase FIR.
    // Flat, or nothing there, is a bypass rather than a convolution with an impulse.
    void rebuildCurves(int side) {
        std::vector<double> sum(commonGrid_.size(), 0.0);
        bool any = false;
        for (const auto& t : tones()) {
            // A knob whose positions carry BANDS ships no curve and no grid: it belongs to the other
            // path, and asking this one for its curve would ask an empty grid for a value.
            if (sideOf(t) != side || t.positions.empty() || bandsPerPosition(t)) continue;
            const auto grid  = gridOf(t);
            const auto curve = curveAt(t, grid, toneAt_[t.name]);
            if (curve.empty()) continue;
            any = true;
            for (std::size_t k = 0; k < sum.size(); ++k)
                sum[k] += felitronics::lineareq::curveDbAt(curve, grid, commonGrid_[k]);
        }
        curveSum_[side] = any ? sum : std::vector<double> {};
        auto taps = any ? felitronics::lineareq::magnitudeCurveToFir(sum, commonGrid_, fs_, kFirTaps, kFirDesign)
                        : std::vector<float> {};
        firTaps_[side] = std::move(taps);
        if (! firTaps_[side].empty() && prepared_) {
            const float* t[1] { firTaps_[side].data() };
            fir_[side].loadIR(t, 1, (int) firTaps_[side].size(), fs_);
            firBypass_[side].store(false, std::memory_order_release);
        } else
            firBypass_[side].store(true, std::memory_order_release);
    }

    // Every band of every section-form knob on one side, at its knob's position — in a fixed order,
    // so band k keeps its state across knob moves.
    void rebuildBands(int side) {
        BandSet set;
        bandsShown_[side].clear();
        int dropped = 0;
        for (const auto& t : tones()) {
            if (sideOf(t) != side) continue;
            // A DIAL's bands travel with its rotation; a SWITCH's stand whole at the position it is on.
            // One list either way, in a fixed order, so band k keeps its state across knob moves.
            const auto qs = ! t.sections.empty() ? sectionsAt(t, toneAt_[t.name], fs_)
                          : bandsPerPosition(t)  ? sectionsAtValue(t, toneAt_[t.name], fs_)
                                                 : std::vector<felitronics::rigplayer::SectionBiquad> {};
            for (const auto& q : qs) {
                if (set.count < kMaxBands) { set.c[set.count++] = q; bandsShown_[side].push_back(q); }
                else ++dropped;
            }
        }
        bandsDropped_ = dropped;
        publishBands(side, set);
    }

    // The first blend knob of the stage — a pedal has one. Its dry path's shape is designed once; the
    // gains move with the knob and never rebuild a filter.
    const namz::rig::Blend* blendKnob() const { return stage_.blend.empty() ? nullptr : &stage_.blend.front(); }

    void rebuildDry() {
        const auto* b = blendKnob();
        dryTaps_.clear();
        if (b != nullptr) {
            const auto grid = gridOf(*b);
            dryTaps_ = felitronics::lineareq::magnitudeCurveToFir(dryCurve(*b, grid), grid, fs_, kFirTaps, kFirDesign);
        }
        if (! dryTaps_.empty() && prepared_) {
            const float* t[1] { dryTaps_.data() };
            dry_.loadIR(t, 1, (int) dryTaps_.size(), fs_);
            dryFirBypass_.store(false, std::memory_order_release);
        } else
            dryFirBypass_.store(true, std::memory_order_release);
        applyBlendGains();
        dryActive_.store(b != nullptr, std::memory_order_release);
    }

    void applyBlendGains() {
        const auto* b = blendKnob();
        const BlendGains g = b != nullptr ? blendGainsAt(*b, blendAt_[b->name]) : BlendGains {};
        dryGain_.store((float) std::clamp(g.dry, -4.0, 4.0), std::memory_order_release);
        wetGain_.store((float) std::clamp(g.wet,  0.0, 4.0), std::memory_order_release);
    }

    // ------------------------------------------------------------ the bands, across the threads ----

    struct BandSet {
        int count = 0;
        felitronics::rigplayer::SectionBiquad c[kMaxBands];
    };
    // A small pool instead of a heap: the message thread fills a free set and publishes its index; the
    // audio thread takes the index, copies the set and frees it. At most one pending and one being
    // filled at a time, so four never run out — and nothing is allocated or freed on the audio thread.
    struct BandSlot { BandSet set; std::atomic<bool> busy { false }; };
    static constexpr int kBandPool = 4;

    void publishBands(int side, const BandSet& set) {
        auto& pool = bandPool_[side];
        for (int i = 0; i < kBandPool; ++i) {
            if (pool[i].busy.load(std::memory_order_acquire)) continue;
            pool[i].set = set;
            pool[i].busy.store(true, std::memory_order_release);
            const int old = bandNext_[side].exchange(i, std::memory_order_acq_rel);
            if (old >= 0) pool[old].busy.store(false, std::memory_order_release);   // never taken: superseded
            return;
        }
    }

    // Audio thread, at the top of every block.
    void takeBands(int side) {
        const int i = bandNext_[side].exchange(-1, std::memory_order_acq_rel);
        if (i < 0) return;
        auto& rt = bandRt_[side];
        const auto& in = bandPool_[side][i].set;
        // A band that appears (a load) starts AT its coefficients — there is nothing to ramp from; one
        // that is already running ramps to the new ones across the next block.
        for (int k = rt.count; k < in.count; ++k) {
            bandCur_[side][k] = in.c[k];
            for (auto& b : bq_[side][k]) b.reset();
        }
        rt = in;
        bandPool_[side][i].busy.store(false, std::memory_order_release);
    }

    // The same five numbers, bit for bit: a band that did not move ramps nothing.
    static bool same(const felitronics::rigplayer::SectionBiquad& x, const felitronics::rigplayer::SectionBiquad& y) {
        return std::memcmp(&x.b0, &y.b0, sizeof(double)) == 0 && std::memcmp(&x.b1, &y.b1, sizeof(double)) == 0
            && std::memcmp(&x.b2, &y.b2, sizeof(double)) == 0 && std::memcmp(&x.a1, &y.a1, sizeof(double)) == 0
            && std::memcmp(&x.a2, &y.a2, sizeof(double)) == 0;
    }

    // Every band of a side, on every channel. Coefficients that changed since the last block travel
    // linearly across this one — a knob dragged sixty times a second never steps a filter.
    void runBands(int side, float* const* planes, int count) {
        auto& rt = bandRt_[side];
        for (int k = 0; k < rt.count; ++k) {
            const auto& to = rt.c[k];
            auto& from = bandCur_[side][k];
            auto& bands = bq_[side][k];
            if (same(from, to)) {
                for (int c = 0; c < channels_; ++c) {
                    auto& bq = bands[(std::size_t) c];
                    bq.c = to;
                    float* x = planes[c];
                    for (int i = 0; i < count; ++i) x[i] = bq.processSample(x[i]);
                    bq.flushDenormals();
                }
            } else {
                const double inv = 1.0 / (double) count;
                for (int c = 0; c < channels_; ++c) {
                    auto& bq = bands[(std::size_t) c];
                    float* x = planes[c];
                    for (int i = 0; i < count; ++i) {
                        const double t = (double) (i + 1) * inv;
                        bq.c.b0 = from.b0 + (to.b0 - from.b0) * t;
                        bq.c.b1 = from.b1 + (to.b1 - from.b1) * t;
                        bq.c.b2 = from.b2 + (to.b2 - from.b2) * t;
                        bq.c.a1 = from.a1 + (to.a1 - from.a1) * t;
                        bq.c.a2 = from.a2 + (to.a2 - from.a2) * t;
                        x[i] = bq.processSample(x[i]);
                    }
                    bq.flushDenormals();
                }
                from = to;
            }
        }
    }

    // One linear ramp per block toward an exponential ~10 ms endpoint: smooth moves without an exp()
    // per sample. `current` is audio-thread-only.
    static bool isOne(float g) { const float one = 1.0f; return std::memcmp(&g, &one, sizeof(float)) == 0; }   // bit for bit, as same()
    void rampInto(float* const* planes, int nch, int count, float want, float& current) const {
        const float decay = std::exp(-(float) count / (float) (0.010 * fs_));
        const float end   = want + (current - want) * decay;
        const float step  = (end - current) / (float) count;
        if (! isOne(current) || ! isOne(want))
            for (int c = 0; c < nch; ++c) {
                float g = current;
                float* x = planes[c];
                for (int i = 0; i < count; ++i) { g += step; x[i] *= g; }
            }
        current = end;
    }

    // ------------------------------------------------------------------------------------ state ----

    // message thread
    struct Model { std::string id; std::shared_ptr<const std::vector<std::byte>> bytes; };   // null = not fetched yet
    std::vector<Model> models_;                        // BlendModelId = index + 1; 0 = nothing
    std::uint64_t      gen_ = 0;                       // bumped by every unload(): a job carries the one it was taken under
    bool               jobOut_ = false;                // a job is with the host; nothing else is handed out
    std::vector<int>   fileModel_;                     // Device::files index → models_ index
    namz::rig::Stage   stage_;
    ModelSource        source_;
    bool               loaded_ = false;
    Settings           axes_;
    std::string        dial_;
    int                dialSweep_ = 0;
    double             dialDeg_ = 0.0;
    std::map<std::string, std::string> toneAt_, blendAt_;
    std::vector<namz::rig::Tone> toneOverride_;        // tone handed in beside the manifest; empty = none
    Selection          sel_;
    SlotPlan           plan_;
    BlendShape         shape_ {};
    double             topExtend_ = kTopExtendDbPerDeg;
    AlignmentTable     align_;
    long long          loads_ = 0;
    int                bandsDropped_ = 0;
    double             coldSeconds_ = kColdAfterSeconds;
    // A 1/12-octave grid from 20 Hz to 20 kHz, on which every curve-form knob is summed whatever grid
    // it shipped on.
    std::vector<double> commonGrid_ = felitronics::lineareq::logFreqGrid(20.0, 20000.0, 121);
    std::vector<double> curveSum_[2];
    std::vector<float>  firTaps_[2], dryTaps_;
    std::vector<felitronics::rigplayer::SectionBiquad> bandsShown_[2];

    // both, by contract
    double fs_ = 48000.0;
    int    maxBlock_ = 0, channels_ = 1;
    bool   prepared_ = false;

    // the mailboxes (message → audio, audio → message)
    felitronics::nam::NamStage nam_[2];
    std::atomic<std::uint64_t> reqWant_[2] { 0, 0 };
    std::atomic<float>         reqTarget_ { 0.0f };
    std::atomic<std::uint64_t> loadAsk_ { 0 };          // (model << 1) | slot; 0 = nothing asked
    std::atomic<bool>          landFlag_ { false };
    std::atomic<int>           landSlot_ { 0 };
    std::atomic<std::uint64_t> landModel_ { 0 };
    std::atomic<long long>     landWarm_ { 0 };
    std::atomic<int>           landDelay_ { 0 };
    std::atomic<int>           landFail_ { 0 };          // slot + 1; 0 = none
    std::atomic<bool>          forget_ { false };
    std::atomic<int>           pendDelay_[2] { 0, 0 };
    std::atomic<int>           slotDelay_[2] { 0, 0 };
    std::atomic<std::uint64_t> held_[2] { 0, 0 };
    std::atomic<float>         liveMix_ { 0.0f };
    std::atomic<int>           warmBlocks_ { 0 }, mixJumps_ { 0 };
    std::atomic<float>         biggestJump_ { 0.0f };
    std::atomic<long long>     coldAfter_ { 0 };          // coldSeconds_ at fs_, for the law; 0 = never
    std::atomic<bool>          slotCold_[2] { false, false };
    std::atomic<int>           coldBlocks_[2] { 0, 0 };
    std::atomic<float>         chainGain_ { 1.0f };
    std::atomic<float>         slotGain_[2] { 1.0f, 1.0f };
    std::atomic<bool>          normalize_ { false };
    std::atomic<bool>          inputTrims_ { true };
    std::atomic<bool>          firBypass_[2] { true, true };
    std::atomic<bool>          dryActive_ { false }, dryFirBypass_ { true };
    std::atomic<float>         dryGain_ { 0.0f }, wetGain_ { 1.0f };
    BandSlot                   bandPool_[2][kBandPool];
    std::atomic<int>           bandNext_[2] { -1, -1 };

    // audio thread only
    felitronics::nam::BlendState blend_ {};
    bool coldRt_[2] {};                                // the law's cold flags as of the last block, to see a slot fall asleep
    std::array<float, (std::size_t) kMaxDelay> lagTail_[2][kMaxChannels] {};
    float curChain_ = 1.0f, curSlot_[2] { 1.0f, 1.0f }, curDry_ = 0.0f, curWet_ = 1.0f;
    BandSet                    bandRt_[2];
    felitronics::rigplayer::SectionBiquad   bandCur_[2][kMaxBands];
    felitronics::eq::Biquad    bq_[2][kMaxBands][kMaxChannels];
    felitronics::convolution::CabConvolver fir_[2], dry_;
    std::vector<float> slotB_[kMaxChannels], dryBuf_[kMaxChannels], spare_[kMaxChannels];
};

} // namespace felitronics::rigplayer