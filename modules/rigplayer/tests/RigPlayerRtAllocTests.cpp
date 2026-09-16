// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// 🔴 THE PLAYER'S RT GATE, AND UNTIL NOW THIS CLASS HAD NONE. `RigPlayer.h`'s own header says it in the
// first paragraph — "process() alone on the audio thread — it allocates nothing, locks nothing, touches
// no file" — and until P90 that sentence was checked by a PROBE somebody ran once (P85) rather than by
// anything that runs again. A claim with no gate is a claim that stops being true quietly, which is the
// whole lesson of this tree's law 9: the difference between a tier and a test suite is the gate.
//
// It is cheap now and was not before: `test_support/alloc_counter.h` (P52) replaces all eight `new` forms
// — plain, array, over-aligned, nothrow and their combinations — and proves itself before main(). This
// player's audio path runs THROUGH `core::AlignedVector`/`SeamAllocator<64>` (every `CabConvolver` here
// owns one), so the over-aligned form is not theory for this file: a two-form counter is exactly the
// blind instrument that header was written to retire.
//
// WHAT IS GATED, and it is deliberately wider than `process()`. The audio thread of a host running this
// player calls THREE things, and the class promises all three:
//   · `process()` — every block, on every shape the player has: bands, a curve FIR, a dry blend, two
//     models, a crossfade, a slot asleep, a slot warming.
//   · `reset()` — the stream restart P86 added, "callable from the audio thread — nothing here
//     allocates, locks, throws or touches a field the message thread owns".
//   · the read-outs a strip polls beside the callback.
// The message thread's verbs — `prepare()`, `load()`, `service()`, `deliver()` — allocate by design and
// are NOT gated; they are exercised only to reach the state the gated calls are measured in.
//
// ⚠️ THE CARVE-OUT IS NAM'S, WORD FOR WORD, and it is why the fixtures here are the shapes they are.
// `NamStage`'s own header records that some real captures allocate INSIDE `process()` — `wavenet_a2_max`
// makes 1024 allocations for one stereo 256-block — so a gate that drove an arbitrary capture would be
// measuring NAM's upstream, not this player. The fixtures are therefore the tree's own minimal Linear and
// WaveNet models, whose inference is allocation-free, so a non-zero delta names THIS file.
//
// ⚠️ AND THE CONTROL IS PART OF THE GATE. A no-allocation assertion is the kind that passes when the
// instrument is dead — the exact failure `alloc_counter.h`'s "no macro" note exists to prevent one form
// of. So the last group PLANTS an allocation on the same thread, through the same counter, and asserts
// the gate goes RED: it is built from the failure it is written for, rather than trusted because it is
// green. The planting goes through `::operator new` and not a `new T`, because a new-expression is
// elidable from -O1 and USING the result does not protect it — a working gate would then read zero and
// look blind.

#include <alloc_counter.h>
#include <felitronics_test.h>
#include <felitronics/rigplayer/RigPlayer.h>

#include <namz.h>

#include <cmath>
#include <functional>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

using felitronics::test::group;
using felitronics::test::ok;
using namespace felitronics::rigplayer;

namespace {

constexpr int kBlock = 256;

std::string gainModel(double w) {
    char buf[256];
    std::snprintf(buf, sizeof buf,
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[%.6f],"sample_rate":48000})", w);
    return buf;
}

// A capture with real memory, so the warm-up ledger and the per-slot alignment tail are live rather
// than skipped — the two places this player does per-block bookkeeping of its own.
std::string delayModel(int d) {
    std::string w = "[";
    for (int i = 0; i < d; ++i) w += "0.0,";
    w += "1.0]";
    return R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":)" + std::to_string(d + 1)
         + R"(,"bias":false,"implementation":"direct"},"weights":)" + w + R"(,"sample_rate":48000})";
}

// The minimal WaveNet the nam suite uses: one one-channel layer, allocation-free inference.
std::string waveNetModel(int dilation) {
    return std::string(R"({"version":"0.5.0","architecture":"WaveNet","config":{"layers":[{"input_size":1,)")
         + R"("condition_size":1,"head_size":1,"head_bias":false,"channels":1,"kernel_size":2,"dilations":[)"
         + std::to_string(dilation)
         + R"(],"activation":"Tanh","gated":false}],"head_scale":1.0},"weights":[1,1,0,0,1,0,0,1,1],)"
         + R"("sample_rate":48000})";
}

std::vector<std::byte> bytesOf(const std::string& s) {
    const auto* p = reinterpret_cast<const std::byte*>(s.data());
    return { p, p + s.size() };
}

// A device with a gain dial over three knots and a tone knob on each side — bands after the models,
// a curve before them — plus a blend knob, so every stage of the chain in the header's signal diagram
// is actually running while the counter is armed. A curve knob is what puts a `CabConvolver` (and with
// it the over-aligned allocator) on the audio path at all.
namz::rig::Rig rigOf() {
    namz::rig::Rig r;
    namz::rig::Stage st; st.kind = namz::rig::StageKind::Nam; st.rawKind = "nam"; st.slot = "pedal";
    namz::rig::Control g; g.name = "gain"; g.role = namz::rig::Role::Gain;
    g.values = { "60", "150", "240" }; g.sweep = 300;
    st.device.controls = { g };
    namz::rig::FileEntry fe; fe.id = "early"; fe.settings = { { "gain", "60" } };
    namz::rig::FileEntry fm; fm.id = "mid";   fm.settings = { { "gain", "150" } };
    namz::rig::FileEntry fl; fl.id = "late";  fl.settings = { { "gain", "240" } };
    st.device.files = { fe, fm, fl };

    namz::rig::Tone tone;
    tone.name = "tone"; tone.sweep = 300; tone.placement = "post"; tone.reference = "150"; tone.defaultValue = "150";
    namz::rig::Section hs; hs.kind = namz::rig::SectionKind::HighShelf; hs.hz = 3000.0; hs.q = 0.7;
    hs.dbAtMin = -6.0; hs.dbAtMax = 6.0;
    tone.sections = { hs };

    namz::rig::Tone bass;
    bass.name = "bass"; bass.sweep = 300; bass.placement = "pre"; bass.reference = "0"; bass.defaultValue = "0";
    bass.grid.fLo = 20.0; bass.grid.fHi = 20000.0; bass.grid.points = 25;
    bass.trusted.levels = 1;
    const auto grid = felitronics::lineareq::logFreqGrid(20.0, 20000.0, 25);
    namz::rig::TonePosition p0; p0.value = "0"; p0.norm = 0.0; p0.db.assign(25, 0.0);
    namz::rig::TonePosition p1; p1.value = "300"; p1.norm = 1.0;
    for (const double f : grid)
        p1.db.push_back(f <= 200.0 ? 6.0 : f >= 2000.0 ? 0.0 : 6.0 * (1.0 - std::log(f / 200.0) / std::log(10.0)));
    bass.positions = { p0, p1 };
    st.tone = { tone, bass };

    namz::rig::Blend mix;
    mix.name = "mix"; mix.sweep = 300; mix.reference = "300"; mix.dryEnd = "0"; mix.defaultValue = "300";
    mix.polarity = 1; mix.dryLevelDb = -6.0;
    mix.grid.fLo = 20.0; mix.grid.fHi = 20000.0; mix.grid.points = 25;
    mix.dryDb.assign(25, 0.0);
    namz::rig::BlendPosition dry; dry.value = "0";   dry.norm = 0.0; dry.dryDb = 0.0;    dry.wetDb = -120.0;
    namz::rig::BlendPosition wet; wet.value = "300"; wet.norm = 1.0; wet.dryDb = -120.0; wet.wetDb = 0.0;
    mix.positions = { dry, wet };
    st.blend = { mix };

    r.chain = { st };
    return r;
}

struct Bench {
    RigPlayer p;
    std::map<std::string, std::vector<std::byte>> files;
    double fs;
    int channels;
    double phase = 0.0;

    Bench(const std::string& nam, double sampleRate, int nch) : fs(sampleRate), channels(nch) {
        files = { { "early", bytesOf(nam) }, { "mid", bytesOf(nam) }, { "late", bytesOf(nam) } };
        felitronics::test::run (p.prepare(fs, kBlock, channels));
        p.load(rigOf(), [this](const std::string& id) {
            const auto it = files.find(id);
            return it == files.end() ? std::vector<std::byte> {} : it->second;
        });
    }
    // Message-thread work, run to completion OUTSIDE any armed window: this is where the player is
    // allowed to allocate, and where the models actually land.
    void settle(int blocks) {
        for (int b = 0; b < blocks; ++b) { one(0.3); p.serviceHere(); }
    }
    void one(double a) {
        std::vector<float> l((std::size_t) kBlock), r((std::size_t) kBlock);
        for (int i = 0; i < kBlock; ++i) {
            const float s = (float) (a * std::sin(phase));
            phase += 2.0 * 3.14159265358979323846 * 220.0 / fs;
            if (phase > 6.283185307179586) phase -= 6.283185307179586;
            l[(std::size_t) i] = s; r[(std::size_t) i] = -s;
        }
        float* io[2] { l.data(), r.data() };
        felitronics::test::run (p.process(io, channels, kBlock));
    }
};

// The delta the counter saw across `body`, with NOTHING else on this thread in between.
long long allocsAcross(const std::function<void()>& body) {
    const long long before = alloc::count.load(std::memory_order_relaxed);
    body();
    return alloc::count.load(std::memory_order_relaxed) - before;
}

// Every block of audio the gate drives, for one already-settled player. Kept apart from the counting so
// the buffers it needs are allocated BEFORE the window, not inside it.
struct Driver {
    std::vector<float> l, r;
    float* io[2];
    Driver() : l((std::size_t) kBlock), r((std::size_t) kBlock) { io[0] = l.data(); io[1] = r.data(); }
    void fill(double a, double& phase, double fs) {
        for (int i = 0; i < kBlock; ++i) {
            const float s = (float) (a * std::sin(phase));
            phase += 2.0 * 3.14159265358979323846 * 220.0 / fs;
            if (phase > 6.283185307179586) phase -= 6.283185307179586;
            l[(std::size_t) i] = s; r[(std::size_t) i] = -s;
        }
    }
};

} // namespace

int main() {
    // ------------------------------------------------------------------------------------------
    group("the counter itself is live in THIS binary before anything is claimed with it");
    {
        // The header proves its eight replacements before main(); what is proved here is that the
        // replacements are reachable from this translation unit's own numbers, so a zero below is a
        // measurement rather than a link-time accident.
        const long long seen = allocsAcross([] {
            auto* p = ::operator new(64);
            ::operator delete(p);
        });
        ok(seen == 1, "one ::operator new is counted as exactly one (" + std::to_string(seen) + ")");
        const long long none = allocsAcross([] {});
        ok(none == 0, "…and an empty window counts none (" + std::to_string(none) + ")");
    }

    // ------------------------------------------------------------------------------------------
    group("process() allocates nothing — every shape, mono and stereo");
    {
        struct Case { const char* name; std::string json; double fs; int nch; };
        const Case cases[] {
            { "a memoryless capture, mono, 48 kHz",           gainModel(0.5),    48000.0, 1 },
            { "a memoryless capture, stereo, 48 kHz",         gainModel(0.5),    48000.0, 2 },
            { "a 2001-tap capture, stereo, 48 kHz",           delayModel(2000),  48000.0, 2 },
            // 44.1 kHz puts the rate-matcher in: the model is tagged 48 kHz, so both resampler legs
            // run and the player's own dry aligner is holding a non-zero delay.
            { "a 2001-tap capture, stereo, 44.1 kHz (resampling)", delayModel(2000), 44100.0, 2 },
            { "a WaveNet, stereo, 96 kHz (resampling)",       waveNetModel(512), 96000.0, 2 },
        };
        for (const auto& c : cases) {
            Bench b(c.json, c.fs, c.nch);
            b.p.setDial("gain", 150.0);
            b.settle(40);                                   // land both models, off the clock
            Driver d;
            double ph = 0.0;
            d.fill(0.3, ph, c.fs);                          // …and fill once, also off the clock

            const long long steady = allocsAcross([&] {
                for (int k = 0; k < 24; ++k)
                    felitronics::test::run (b.p.process(d.io, c.nch, kBlock));
            });
            felitronics::test::okNoAlloc(steady == 0, std::string("24 steady blocks allocate nothing (") + std::to_string(steady)
                            + ") — " + c.name);

            // A TURN MID-BLOCK is the interesting one: the knobs are published by the message thread,
            // but the audio thread is what TAKES them — new bands, a new weight goal, a staged retime.
            b.p.setDial("gain", 240.0);
            const long long turning = allocsAcross([&] {
                for (int k = 0; k < 24; ++k)
                    felitronics::test::run (b.p.process(d.io, c.nch, kBlock));
            });
            felitronics::test::okNoAlloc(turning == 0, std::string("…and 24 blocks across a published turn allocate nothing (")
                             + std::to_string(turning) + ") — " + c.name);

            // DIGITAL SILENCE. What this does NOT do is put a slot to sleep — sleep is a property of a
            // WEIGHT at rest under an unchanged request, not of the input, and 24 blocks are a fraction of
            // the cold window anyway. An earlier version of this comment claimed it reached the cold path;
            // a drain-path allocation planted to check that claim passed the whole gate. The cold path and
            // the dropped-lane drain have their own group below.
            const long long hush = allocsAcross([&] {
                std::fill(d.l.begin(), d.l.end(), 0.0f);
                std::fill(d.r.begin(), d.r.end(), 0.0f);
                for (int k = 0; k < 24; ++k)
                    felitronics::test::run (b.p.process(d.io, c.nch, kBlock));
            });
            felitronics::test::okNoAlloc(hush == 0, std::string("…and 24 blocks of silence allocate nothing (") + std::to_string(hush)
                          + ") — " + c.name);
        }
    }

    // ------------------------------------------------------------------------------------------
    // 🔴 THE TWO PATHS THIS PLAYER RUNS ONLY WHEN SOMETHING STOPS, and they are the ones a gate built
    // from steady blocks never enters. (1) A slot the law puts to SLEEP is handed a width-zero call and
    // its stage feeds itself the silence it still owes (law 11c) — reachable only after a whole cold
    // window at rest under an unchanged request. (2) A host that drops from stereo to mono leaves lane 1
    // to drain on its falling edge, in every stage, the same way. Measured before this group existed: an
    // allocation planted inside NamStage's falling-edge drain passed all 29 checks of this file.
    group("the SLEEP transition and a DROPPED lane — the drain paths — allocate nothing");
    {
        struct Case { const char* name; std::string json; double fs; };
        const Case cases[] {
            { "a 2001-tap capture, 48 kHz",              delayModel(2000),  48000.0 },
            { "a 2001-tap capture, 44.1 kHz (resampling)", delayModel(2000), 44100.0 },
            { "a WaveNet, 96 kHz (resampling)",          waveNetModel(512), 96000.0 },
        };
        for (const auto& c : cases) {
            {   // (1) the sleep transition, with the drain that follows it inside the window
                Bench b(c.json, c.fs, 2);
                b.p.setBlendShape({ 0.5, 0.0 });           // STEP: the neighbour sits at exactly zero
                b.p.setDial("gain", 150.0);
                Driver d; double ph = 0.0; d.fill(0.3, ph, c.fs);
                const int rest = (int) std::ceil((RigPlayer::kColdAfterSeconds - 0.2) * c.fs / kBlock);
                for (int k = 0; k < rest; ++k) { felitronics::test::run (b.p.process(d.io, 2, kBlock)); b.p.serviceHere(); }
                const bool awakeBefore = ! b.p.slotCold(0) && ! b.p.slotCold(1);
                const int window = (int) std::ceil(0.6 * c.fs / kBlock);   // spans the 2 s line and the drain after it
                const long long seen = allocsAcross([&] {
                    for (int k = 0; k < window; ++k) felitronics::test::run (b.p.process(d.io, 2, kBlock));
                });
                const bool asleepAfter = b.p.slotCold(0) || b.p.slotCold(1);
                ok(awakeBefore && asleepAfter,
                   std::string("precondition: a slot fell asleep INSIDE the measured window — ") + c.name);
                felitronics::test::okNoAlloc(seen == 0,
                   std::string("a slot falling asleep, and draining at width zero, allocates nothing (")
                   + std::to_string(seen) + ") — " + c.name);
            }
            {   // (2) stereo -> mono -> nothing: lane 1 drains first, then lane 0
                Bench b(c.json, c.fs, 2);
                b.p.setDial("gain", 150.0);
                b.settle(40);
                Driver d; double ph = 0.0; d.fill(0.3, ph, c.fs);
                for (int k = 0; k < 8; ++k) felitronics::test::run (b.p.process(d.io, 2, kBlock));
                const long long mono = allocsAcross([&] {
                    for (int k = 0; k < 48; ++k) felitronics::test::run (b.p.process(d.io, 1, kBlock));
                });
                felitronics::test::okNoAlloc(mono == 0,
                   std::string("a host dropping to MONO — lane 1 draining — allocates nothing (")
                   + std::to_string(mono) + ") — " + c.name);
                const long long gap = allocsAcross([&] {
                    for (int k = 0; k < 48; ++k) felitronics::test::run (b.p.process(d.io, 0, kBlock));
                });
                felitronics::test::okNoAlloc(gap == 0,
                   std::string("…and a GAP — both lanes draining — allocates nothing (")
                   + std::to_string(gap) + ") — " + c.name);
            }
        }
    }

    // ------------------------------------------------------------------------------------------
    group("a call LONGER than the prepared block is chunked, not allocated for (law 11)");
    {
        Bench b(delayModel(2000), 48000.0, 2);
        b.p.setDial("gain", 150.0);
        b.settle(40);
        // Four times the prepared block: `process` walks it in `maxBlock` steps. Buffers made outside.
        std::vector<float> l((std::size_t) (4 * kBlock), 0.05f), r((std::size_t) (4 * kBlock), -0.05f);
        float* io[2] { l.data(), r.data() };
        const long long seen = allocsAcross([&] {
            for (int k = 0; k < 8; ++k)
                felitronics::test::run (b.p.process(io, 2, 4 * kBlock));
        });
        felitronics::test::okNoAlloc(seen == 0, "a 4x-maxBlock call allocates nothing (" + std::to_string(seen) + ")");
    }

    // ------------------------------------------------------------------------------------------
    group("reset() — the audio-thread stream restart — allocates nothing either");
    {
        struct Case { const char* name; std::string json; double fs; };
        const Case cases[] {
            { "a memoryless capture",             gainModel(0.5),    48000.0 },
            { "a 2001-tap capture",               delayModel(2000),  48000.0 },
            { "a 2001-tap capture, resampling",   delayModel(2000),  44100.0 },
            { "a WaveNet, resampling",            waveNetModel(512), 96000.0 },
        };
        for (const auto& c : cases) {
            Bench b(c.json, c.fs, 2);
            b.p.setDial("gain", 150.0);
            b.settle(40);
            Driver d; double ph = 0.0; d.fill(0.3, ph, c.fs);
            for (int k = 0; k < 4; ++k) felitronics::test::run (b.p.process(d.io, 2, kBlock));

            // The FIRST restart is the expensive one — both slots owe a whole drain — and it is the one
            // that must not allocate. The second is idempotent and costs nothing; both are measured.
            const long long first = allocsAcross([&] { b.p.reset(); });
            felitronics::test::okNoAlloc(first == 0, std::string("a restart with both lanes dirty allocates nothing (")
                           + std::to_string(first) + ") — " + c.name);
            const long long second = allocsAcross([&] { b.p.reset(); });
            felitronics::test::okNoAlloc(second == 0, std::string("…and a second, idempotent one too (") + std::to_string(second)
                            + ") — " + c.name);
        }
    }

    // ------------------------------------------------------------------------------------------
    group("the read-outs a strip polls beside the callback allocate nothing");
    {
        Bench b(delayModel(2000), 48000.0, 2);
        b.p.setDial("gain", 150.0);
        b.settle(40);
        // `heldFileId()` returns a std::string BY VALUE and is therefore expected to allocate for a name
        // longer than the small-string buffer; it is a message-thread read-out and is not claimed here.
        // What is claimed is the set a meter polls per block.
        const long long seen = allocsAcross([&] {
            volatile float sink = 0.0f;
            volatile int isink = 0;
            for (int k = 0; k < 64; ++k) {
                sink = sink + b.p.liveMix() + b.p.liveDry() + b.p.liveWet();
                isink = isink + b.p.warmBlocks() + b.p.mixJumps() + b.p.coldBlocks(0)
                      + b.p.appliedSlotDelay(0) + b.p.latencySamples() + (int) b.p.slotCold(1);
            }
            (void) sink; (void) isink;
        });
        felitronics::test::okNoAlloc(seen == 0, "64 rounds of the per-block read-outs allocate nothing (" + std::to_string(seen) + ")");
    }

    // ------------------------------------------------------------------------------------------
    // 🔴 THE CONTROL. Everything above is an assertion that a number is ZERO, and that is the shape of
    // assertion that keeps passing after the instrument dies. This group plants the failure the gate
    // exists to catch — an allocation on the audio thread, inside the measured window — and asserts the
    // gate SEES it. Without this group a linker that dropped the replacements, or a window opened around
    // the wrong call, would read green for ever.
    group("the gate is RED when an allocation is planted in the measured window");
    {
        Bench b(delayModel(2000), 48000.0, 2);
        b.p.setDial("gain", 150.0);
        b.settle(40);
        Driver d; double ph = 0.0; d.fill(0.3, ph, 48000.0);

        // Through `::operator new` and not `new T`: a new-expression is elidable from -O1 — and USING
        // the result does not protect it — so a `new float[…]` here can vanish and leave this control
        // silently vacuous, which is the same blindness in a new coat.
        //
        // AND IT IS A DIFFERENCE, NOT A COUNT, because this assertion is a hard `ok` on every stdlib
        // while the zeroes above are `okNoAlloc` — enforced only where the harness can separate our
        // allocations from the standard library's. Where a stdlib allocates inside `process()` itself,
        // "the planted window reads exactly one" would fail for a reason that has nothing to do with the
        // gate. The identical window without the plant is the baseline, so the control proves the
        // counter sees ONE MORE, whatever the library around it does.
        const auto window = [&](bool plant) {
            return allocsAcross([&] {
                for (int k = 0; k < 24; ++k) {
                    felitronics::test::run (b.p.process(d.io, 2, kBlock));
                    if (plant && k == 11) { void* p = ::operator new(sizeof(float) * (std::size_t) kBlock);
                                            ::operator delete(p); }
                }
            });
        };
        const long long base = window(false);
        const long long planted = window(true);
        ok(planted - base == 1, "one planted allocation inside 24 blocks reads as exactly one more ("
                                + std::to_string(planted) + " against " + std::to_string(base)
                                + ") — the zeroes above are measurements, not a dead gate");

        // …and the OVER-ALIGNED form separately, because that is the half fifty suites in this tree
        // could not see and the half this player's convolvers actually use.
        const long long alignedBefore = alloc::alignedCount.load(std::memory_order_relaxed);
        const long long over = allocsAcross([&] {
            void* p = ::operator new(256, std::align_val_t { 64 });
            ::operator delete(p, std::align_val_t { 64 });
        });
        const long long alignedDelta = alloc::alignedCount.load(std::memory_order_relaxed) - alignedBefore;
        ok(over == 1 && alignedDelta == 1,
           "…and an over-aligned allocation is counted, by both totals (" + std::to_string(over) + ", "
           + std::to_string(alignedDelta) + ")");
    }

    return felitronics::test::report();
}
