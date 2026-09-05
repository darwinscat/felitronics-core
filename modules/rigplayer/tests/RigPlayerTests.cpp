// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The pack player (src/rigplayer), driven the way a hand drives it and checked without a sound card.
// The pack is written by namz's writer and read back by its canonical reader, so what the player is
// handed is what a plugin will be handed. The models are tiny Linear NAMs — a gain, or a pure delay —
// whose output IDENTIFIES them, so "which capture is sounding, how loud, and how far apart" are numbers
// read off the audio, not off a variable.

#include <felitronics_test.h>
#include <felitronics/rigplayer/RigPlayer.h>

#include <namz.h>
#include <namz_rig_load.h>
#include <namz_rig_write.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

using felitronics::test::approx;
using felitronics::test::group;
using felitronics::test::ok;
using namespace felitronics::rigplayer;

namespace {

constexpr double kFs = 48000.0;
constexpr int    kBlock = 256;

// A NAM that is a gain: Linear, receptive field 1, one weight. What comes out says which one it is.
std::string gainModel(double w) {
    char buf[256];
    std::snprintf(buf, sizeof buf,
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[%.6f],"sample_rate":48000})", w);
    return buf;
}

// A NAM that is a gain AND AN OFFSET: y = w*x + b. The offset is what makes this fixture worth
// having — it is the cheapest thing that is not a pure scalar, and a pure scalar cannot tell the
// player's two levels apart. Feed the network half as much and the offset is untouched (w*x/2 + b);
// halve the STAGE's output and the offset halves with everything else ((w*x + b)/2). So the mean of
// the output says which of the two was applied, and the release's central claim — input is drive,
// output is volume — finally has a test that a swap would fail.
std::string biasModel(double w, double b) {
    char buf[320];
    std::snprintf(buf, sizeof buf,
        R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":true,"implementation":"direct"},"weights":[%.6f,%.6f],"sample_rate":48000})", w, b);
    return buf;
}

// A NAM that is a pure delay of `d` samples: the window's first tap is the oldest sample and its last
// the newest, so a weight on the first alone hands back what came in `d` samples ago.
std::string delayModel(int d) {
    std::string w = "[";
    for (int i = 0; i < d; ++i) w += "0.0,";
    w += "1.0]";
    return R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":)" + std::to_string(d + 1)
         + R"(,"bias":false,"implementation":"direct"},"weights":)" + w + R"(,"sample_rate":48000})";
}

std::vector<std::byte> bytesOf(const std::string& s) {
    const auto* p = reinterpret_cast<const std::byte*>(s.data());
    return { p, p + s.size() };
}

std::vector<std::byte> packed(const std::string& nam) {
    const auto z = namz::pack(nam.data(), nam.size());
    const auto* p = reinterpret_cast<const std::byte*>(z.data());
    return { p, p + z.size() };
}

// THE DEVICE UNDER TEST. Two axes: a channel switch and a gain dial captured at 60, 150 and 240 of a
// 300-degree sweep on the green channel, one capture on red. The bottom of the dial is a LINK — gain 0
// plays the 60 capture fed 6 dB softer — spelled the way the pack spells it: a second files[] entry
// pointing at the same file with an input_db. One tone knob as bands (a high shelf after the model),
// one as a curve (a bass lift before it), and a blend knob whose dry path is a wire 6 dB down.
namz::rig::Rig testRig(int polarity = 1) {
    namz::rig::Rig rig;
    rig.rigId = "test-rig"; rig.name = "Test Pedal"; rig.modeledBy = "the test";
    namz::rig::Stage st;
    st.kind = namz::rig::StageKind::Nam; st.rawKind = "nam"; st.slot = "pedal";
    st.make = "Darwin's Cat"; st.model = "Test Pedal"; st.gearType = "pedal";
    auto& d = st.device;
    d.family = "Test Pedal"; d.rigId = "test-rig"; d.slot = "pedal";
    namz::rig::Control ch; ch.name = "channel"; ch.role = namz::rig::Role::Channel; ch.values = { "green", "red" };
    namz::rig::Control g;  g.name = "gain";     g.role = namz::rig::Role::Gain;    g.values = { "60", "150", "240" }; g.sweep = 300;
    d.controls = { ch, g };
    const auto file = [](const char* id, const char* chan, const char* gain, double inDb = 0.0) {
        namz::rig::FileEntry f; f.id = id; f.settings = { { "channel", chan }, { "gain", gain } }; f.inputDb = inDb; return f;
    };
    d.files = { file("g60", "green", "60"), file("g150", "green", "150"), file("g240", "green", "240"),
                file("r150", "red", "150"), file("g60", "green", "0", -6.0) };

    // `tone`: one band after the model, +6 dB at the top of the travel, -6 at the bottom, zero at 150.
    namz::rig::Tone tone;
    tone.name = "tone"; tone.sweep = 300; tone.placement = "post"; tone.reference = "150"; tone.defaultValue = "150";
    namz::rig::Section hs; hs.kind = namz::rig::SectionKind::HighShelf; hs.hz = 3000.0; hs.q = 0.7;
    hs.dbAtMin = -6.0; hs.dbAtMax = 6.0;
    tone.sections = { hs };
    // `bass`: a curve before the model — flat at 0 (the reference), and at 300 a +6 dB lift that runs
    // out between 200 Hz and 2 kHz, on a 25-point grid of its own.
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

    // `mix`: wet at 300 (where the models were captured), dry at 0; the dry path is flat, 6 dB down.
    namz::rig::Blend mix;
    mix.name = "mix"; mix.sweep = 300; mix.reference = "300"; mix.dryEnd = "0"; mix.defaultValue = "300";
    mix.polarity = polarity; mix.dryLevelDb = -6.0;
    mix.grid.fLo = 20.0; mix.grid.fHi = 20000.0; mix.grid.points = 25;
    mix.dryDb.assign(25, 0.0);
    namz::rig::BlendPosition dry; dry.value = "0";   dry.norm = 0.0; dry.dryDb = 0.0;    dry.wetDb = -120.0;
    namz::rig::BlendPosition wet; wet.value = "300"; wet.norm = 1.0; wet.dryDb = -120.0; wet.wetDb = 0.0;
    mix.positions = { dry, wet };
    st.blend = { mix };

    rig.chain = { st };
    return rig;
}

// The same rig, as a plugin would meet it: written by the pack writer, read by the canonical reader.
namz::rig::Rig throughTheFormat(const namz::rig::Rig& rig, bool* okOut = nullptr) {
    return namz::rig::loadRigManifest(namz::rig::writeManifest(rig), okOut);
}

// The same device with the pack's own two levels written on the stage, through the format as well.
namz::rig::Rig levelled(double inDb, double outDb) {
    auto r = testRig();
    r.chain[0].inputDb = inDb; r.chain[0].outputDb = outDb;
    return throughTheFormat(r);
}

struct Bench {
    RigPlayer p;
    std::map<std::string, std::vector<std::byte>> files;
    int fetches = 0;
    double phase = 0.0;
    double fs = kFs;

    Bench(const namz::rig::Rig& rig, int channels = 1, double sampleRate = kFs) : fs(sampleRate) {
        files["g60"]  = packed(gainModel(0.25));           // packed, as the pack ships them…
        files["g150"] = packed(gainModel(0.5));
        files["g240"] = bytesOf(gainModel(1.0));           // …and raw, which the stage takes as well
        files["r150"] = bytesOf(gainModel(0.75));
        p.prepare(fs, kBlock, channels);
        load(rig);
    }
    void load(const namz::rig::Rig& rig) {
        p.load(rig, [this](const std::string& id) {
            ++fetches;
            const auto it = files.find(id);
            return it == files.end() ? std::vector<std::byte> {} : it->second;
        });
    }

    // Run the transport on a sine and read the RMS of channel `channel` over the last `measure` blocks
    // — RMS, because a sampled peak is off by up to a fraction of a sample. Channel 1, when there is
    // one, gets the same sine at `aR`. `serviceHere()` after every block — a host's timer, with the
    // load job run right here rather than on a worker.
    double rms(double hz, double a, int blocks = 48, int measure = 32, int channel = 0, double aR = 0.0) {
        const int nch = p.channels();
        std::vector<float> l((std::size_t) kBlock), r((std::size_t) kBlock);
        double sum = 0.0; long n = 0;
        for (int b = 0; b < blocks + measure; ++b) {
            for (int i = 0; i < kBlock; ++i) {
                const double s = std::sin(phase);
                phase += 2.0 * 3.14159265358979323846 * hz / fs;
                if (phase > 2.0 * 3.14159265358979323846) phase -= 2.0 * 3.14159265358979323846;
                l[(std::size_t) i] = (float) (a * s);
                r[(std::size_t) i] = (float) (aR * s);
            }
            float* io[2] { l.data(), r.data() };
            p.process(io, nch, kBlock);
            p.serviceHere();
            if (b >= blocks) {
                const auto& x = channel == 0 ? l : r;
                for (const float v : x) { sum += (double) v * v; ++n; }
            }
        }
        return std::sqrt(sum / (double) std::max(1L, n));
    }
    // A measurement AGAINST the input: the chain's gain at `hz`, linear.
    double gainAt(double hz, double a = 0.1, int blocks = 48, int measure = 32) {
        return rms(hz, a, blocks, measure) / (a / std::sqrt(2.0));
    }
};

double db(double lin) { return 20.0 * std::log10(std::max(1e-12, lin)); }

// What the test rig's high shelf (3 kHz, Q 0.7) reads between 10 kHz and 100 Hz at `gainDb`, by the
// format's own formula at this rate — the claim is "the player applies that formula", not "a shelf is
// exactly its nominal gain 1.7 octaves up", which it is not.
double shelfDb(double gainDb, double fs) {
    const auto q = felitronics::rigplayer::designSection(felitronics::rigplayer::SectionKind::HighShelf, 3000.0, gainDb, 0.7, fs);
    return felitronics::rigplayer::sectionMagnitudeDb(q, 10000.0, fs) - felitronics::rigplayer::sectionMagnitudeDb(q, 100.0, fs);
}

} // namespace

int main() {
    std::printf("felitronics::rigplayer::RigPlayer tests\n");

    group("the rig round-trips through the pack writer and the canonical reader");
    bool okManifest = false;
    const auto rig = throughTheFormat(testRig(), &okManifest);
    {
        ok(okManifest, "the writer's manifest is a manifest to the reader");
        ok(rig.chain.size() == 1 && rig.chain[0].kind == namz::rig::StageKind::Nam, "one NAM stage");
        const auto& st = rig.chain[0];
        ok(st.device.files.size() == 5, "five file entries, the link among them");
        ok(st.device.files.size() == 5 && st.device.files[4].id == "g60" && st.device.files[4].inputDb == -6.0,
           "the link is a second entry for the same file, with its input_db");
        ok(st.tone.size() == 2, "both tone knobs survive");
        ok(st.tone.size() == 2 && st.tone[0].sections.size() == 1 && st.tone[0].positions.empty(), "the band knob is bands");
        ok(st.tone.size() == 2 && st.tone[1].positions.size() == 2 && st.tone[1].sections.empty(), "the curve knob is a curve");
        ok(st.blend.size() == 1 && st.blend[0].positions.size() == 2, "the blend knob and its two ends");
        ok(crossfadeDial(st.device) != nullptr && crossfadeDial(st.device)->name == "gain", "the gain dial is the crossfade axis");
    }
    const auto& dev = rig.chain[0].device;

    group("selection: which files sound where (pure)");
    {
        Settings green { { "channel", "green" }, { "gain", "150" } };
        auto s = select(dev, green, "gain", 150.0);
        ok(s.knots.size() == 4, "four knots on green: the link at 0, then 60, 150, 240");
        ok(s.fileA == 1 && s.fileB == 2 && s.mixB == 0.0,
           "at 150 exactly the pair is 150 and 240 with nothing of 240: the neighbour is named, so it can be warm");
        s = select(dev, green, "gain", 200.0);
        ok(s.fileA == 1 && s.fileB == 2, "at 200 the pair is 150 and 240");
        approx(s.mixB, 50.0 / 90.0, 1e-9, "…mixed by angle: 50 of the 90 degrees between them");
        auto plan = slotPlan(s, dev);
        ok(plan.file[0] == 1 && plan.file[1] == 2, "150 is the third knot (even) so it keeps slot 0");
        approx(plan.targetB, 50.0 / 90.0, 1e-9, "…and slot 1's weight is the 240 share");
        s = select(dev, green, "gain", 30.0);
        ok(s.fileA == 4 && s.fileB == 0, "at 30 the pair is the link (a knot at 0) and the 60 capture");
        plan = slotPlan(s, dev);
        approx(plan.inputDb[0], -6.0, 1e-9, "…the link's slot is fed 6 dB softer");
        approx(plan.inputDb[1], 0.0, 1e-9, "…the capture's is not");
        s = select(dev, green, "gain", 300.0);
        ok(s.fileA == 2 && s.fileB == 2, "past the top capture it plays alone");
        approx(s.extendDb, 3.0, 1e-9, "…driven 0.05 dB per degree harder: 60 degrees, 3 dB");
        Settings red { { "channel", "red" }, { "gain", "150" } };
        s = select(dev, red, "gain", 240.0);
        ok(s.knots.size() == 1 && s.fileA == 3 && s.fileB == 3, "red has one knot; at 240 it plays alone");
        approx(s.extendDb, 4.5, 1e-9, "…90 degrees past it");
        s = select(dev, green, "", 150.0);
        ok(s.fileA == 1 && s.fileB == 1, "with no dial the panel's own file plays");
    }

    group("the law asks and service answers: silent until fed, then sounding");
    {
        Bench b(rig);
        ok(b.p.loaded(), "loaded");
        ok(b.p.settings().at("channel") == "green" && b.p.settings().at("gain") == "150",
           "the pack's defaults: first channel, the middle of the gain sweep");
        approx(b.p.dialDegrees(), 150.0, 1e-9, "…and the dial stands there");
        // No service: nothing can land, and the law keeps the unfed slots silent.
        std::vector<float> x((std::size_t) kBlock, 0.1f);
        float* io[1] { x.data() };
        for (int i = 0; i < 8; ++i) b.p.process(io, 1, kBlock);
        double e = 0.0; for (const float v : x) e += v * v;
        ok(e == 0.0, "before any model lands the output is silence, not the raw DI");
        ok(b.fetches == 0, "…and nothing was fetched: the law asks, the host answers");
        approx(b.gainAt(1000.0), 0.5, 0.01, "at 150 the 0.5 capture plays");
        ok(b.fetches == 2, "two fetches: the capture sounding, and its neighbour above, warm at zero weight");
        ok(b.p.heldFileId(0) == "g150" && b.p.heldFileId(1) == "g240", "slot 0 holds it, slot 1 the neighbour");
    }

    group("the load as a job: taken once, run anywhere, delivered back — and a stale one dropped");
    {
        Bench b(rig);
        std::vector<float> x((std::size_t) kBlock, 0.1f);
        float* io[1] { x.data() };
        b.p.process(io, 1, kBlock);                                  // the law asks
        b.p.service();                                               // …and service does not load
        ok(b.fetches == 0, "service() fetches nothing: the ask is work for the host");
        auto job = b.p.takeLoadJob();
        // Which slot first is the law's business (it frees the one at weight zero — slot 1, cold): the job
        // names one of the two wanted captures, in the slot the law freed for it.
        ok(job.has_value() && ((job->slot == 0 && job->fileId == "g150") || (job->slot == 1 && job->fileId == "g240")),
           "the job names a wanted capture, for the slot the law freed");
        ok(! b.p.takeLoadJob().has_value(), "one job out: no second one until it is back");
        auto loaded = RigPlayer::run(std::move(*job));                // any thread — here; a worker in the app
        ok(loaded.fetched && loaded.model != nullptr && b.fetches == 1, "run() fetched the bytes and built the model");
        b.p.deliver(std::move(loaded));
        ok(b.p.modelLoads() == 1, "delivered: the fetch is counted and the bytes kept");
        approx(b.gainAt(1000.0), 0.5, 0.01, "…and it sounds (the neighbour arrives the same way)");
        ok(b.fetches == 2, "…each file fetched once");

        // A job out while the pack changes comes back for a pack that is gone — and is dropped whole.
        b.p.setDial("gain", 0.0);                                    // wants the 60 capture
        std::optional<RigPlayer::LoadJob> late;
        for (int i = 0; i < 24 && ! late; ++i) { b.p.process(io, 1, kBlock); b.p.service(); late = b.p.takeLoadJob(); }
        ok(late.has_value() && late->fileId == "g60", "the job for the 60 capture is out");
        b.p.unload();
        b.load(rig);
        ok(! b.p.takeLoadJob().has_value(), "…and nothing new is handed out while it is out");
        b.p.deliver(RigPlayer::run(std::move(*late)));
        b.p.process(io, 1, kBlock);                                  // the new pack's law: nothing landed
        ok(b.p.heldFileId(0).empty() && b.p.heldFileId(1).empty(), "delivered to the new pack it is dropped: no slot holds it");
        approx(b.gainAt(1000.0), 0.5, 0.01, "…and the new pack loads its own and sounds");
    }

    group("a landing published before the pack leaves dies with it");
    {
        // deliver() and unload() on the message thread, with no audio block between them: the landing
        // sat published, the pack left, and the landing used to be consumed AFTER the forget — a stale
        // model in a wiped law, over an emptied stage.
        Bench b(rig);
        std::vector<float> x((std::size_t) kBlock, 0.1f);
        float* io[1] { x.data() };
        b.p.process(io, 1, kBlock);                                  // the law asks
        auto job = b.p.takeLoadJob();
        ok(job.has_value(), "a job is out");
        b.p.deliver(RigPlayer::run(std::move(*job)));                // …and lands, published for the audio thread
        b.p.unload();
        b.load(rig);
        b.p.process(io, 1, kBlock);
        ok(b.p.heldFileId(0).empty() && b.p.heldFileId(1).empty(),
           "the new pack's law holds nothing of the old landing");
        approx(b.gainAt(1000.0), 0.5, 0.01, "…and the new pack loads its own and sounds");
    }

    group("a file that cannot be a model is asked for once, not every block");
    {
        // THE STORM THIS PREVENTS: red's file is broken in this pack. The switch to red asks for it
        // once, fails once, and the captures the hand came from carry on; it used to be re-fetched and
        // re-parsed on every service tick, for ever, fixing nothing.
        Bench b(rig);
        b.files["r150"] = bytesOf("not a model at all");
        approx(b.gainAt(1000.0), 0.5, 0.01, "green 150 sounds");
        ok(b.p.setSwitch("channel", "red"), "the hand switches to red, whose file is broken");
        approx(b.gainAt(1000.0), 0.5, 0.01, "…and the captures it came from carry on");
        ok(b.fetches == 3, "the broken file was asked of the source once");
        ok(b.p.heldFileId(0) == "g150" && b.p.heldFileId(1) == "g240", "…and both slots keep their real models");
        std::vector<float> x((std::size_t) kBlock, 0.1f);
        float* io[1] { x.data() };
        int asks = 0;
        for (int k = 0; k < 100; ++k) {
            std::fill(x.begin(), x.end(), 0.1f);
            b.p.process(io, 1, kBlock);
            b.p.service();
            if (b.p.takeLoadJob().has_value()) ++asks;
        }
        ok(asks == 0 && b.fetches == 3, "a hundred blocks later it has not been asked for again");
        ok(b.p.setSwitch("channel", "green"), "the hand leaves for green…");
        b.p.process(io, 1, kBlock);                              // the law hears the wish change
        ok(b.p.setSwitch("channel", "red"), "…and asks for red again");
        std::optional<RigPlayer::LoadJob> again;
        for (int k = 0; k < 10 && ! again; ++k) { b.p.process(io, 1, kBlock); b.p.service(); again = b.p.takeLoadJob(); }
        ok(again.has_value() && again->fileId == "r150", "a new wish tries the file anew");
        b.p.deliver(RigPlayer::run(std::move(*again)));          // fails again; the slot is refused again
        approx(b.gainAt(1000.0), 0.5, 0.01, "…and the sound never blinked");
    }

    group("the dial: a pair mixed by angle, the extension past the top, each file fetched once");
    {
        Bench b(rig);
        approx(b.gainAt(1000.0), 0.5, 0.01, "150: the 0.5 capture");
        ok(b.p.setDial("gain", 200.0), "the dial turns to 200");
        approx(b.gainAt(1000.0), (40.0 / 90.0) * 0.5 + (50.0 / 90.0) * 1.0, 0.015,
               "200: the 150 and 240 captures, 4/9 and 5/9 of each");
        approx((double) b.p.liveMix(), 50.0 / 90.0, 0.01, "…and that is the weight the audio thread applied");
        b.p.setDial("gain", 300.0);
        approx(b.gainAt(1000.0), 1.0 * std::pow(10.0, 3.0 / 20.0), 0.03, "300: the top capture, 3 dB harder in");
        b.p.setDial("gain", 150.0);
        approx(b.gainAt(1000.0), 0.5, 0.01, "back at 150");
        b.p.setDial("gain", 240.0);
        approx(b.gainAt(1000.0), 1.0, 0.02, "240 alone");
        ok(b.fetches == 2, "two fetches for two files, however many times the dial crossed them");
        ok(b.p.modelLoads() == 2, "…which is what the player counts too");
        const auto& knots = b.p.selection().knots;
        ok(knots.size() == 4 && knots[0].deg == 0.0 && knots[3].deg == 240.0, "the ring's knots: 0 (the link) to 240");
    }

    group("the shape of the handover: where the 50/50 lands and how wide the fade is");
    {
        Bench b(rig);
        ok(b.p.blendShape().point == 0.5 && b.p.blendShape().width == 1.0, "the default is the original law: midpoint, full span");
        b.p.setDial("gain", 195.0);                                  // the middle of 150..240
        approx(b.p.selection().mixB, 0.5, 1e-9, "…so the 50/50 sits in the middle of the pair");
        b.p.setBlendShape({ 0.25, 1.0 });
        ok(b.p.blendShape().point == 0.25, "the point moves to a quarter of the span");
        b.p.setDial("gain", 172.5);                                  // 150 + 0.25 * 90
        approx(b.p.selection().mixB, 0.5, 1e-9, "…and the 50/50 sits there now");
        approx(b.gainAt(1000.0), 0.5 * 0.5 + 0.5 * 1.0, 0.015, "…which is what sounds: half of each capture");
        b.p.setDial("gain", 161.25);                                 // halfway up the near side
        approx(b.p.selection().mixB, 0.25, 1e-9, "the near side is a quarter of the way at its own half");
        b.p.setBlendShape({ 0.25, 0.0 });
        b.p.setDial("gain", 170.0);
        approx(b.p.selection().mixB, 0.0, 1e-9, "width zero: below the point the lower capture alone");
        b.p.setDial("gain", 175.0);
        approx(b.p.selection().mixB, 1.0, 1e-9, "…above it the upper alone — a step, where the bench asked for one");
        approx(b.gainAt(1000.0), 1.0, 0.02, "…and that is what sounds");
        b.p.setBlendShape({});
        b.p.setDial("gain", 195.0);
        approx(b.p.selection().mixB, 0.5, 1e-9, "back to the law");
    }

    group("every knob by name: a dial in degrees, a switch by value, its position read back");
    {
        Bench b(rig);
        ok(b.p.knobValue("tone") == "150" && b.p.knobValue("bass") == "0" && b.p.knobValue("mix") == "300",
           "the tone and blend knobs start where the pack says");
        ok(b.p.knobValue("gain") == "150" && b.p.knobValue("channel") == "green", "…and so do the captured axes");
        ok(b.p.setDial("tone", 300.0) && b.p.knobValue("tone") == "300", "a band knob turned by degrees reads back in degrees");
        ok(b.p.setDial("bass", 210.0) && b.p.knobValue("bass") == "210", "a curve knob the same");
        ok(b.p.setDial("mix", 0.0) && b.p.knobValue("mix") == "0", "the blend knob the same");
        ok(b.p.setSwitch("tone", "0") && b.p.knobValue("tone") == "0", "…or set to a value outright");
        ok(b.p.setDial("gain", 200.0) && b.p.knobValue("gain") == "150",
           "a captured dial reads back the knot at or below the hand; the angle is dialDegrees()");
        approx(b.p.dialDegrees(), 200.0, 1e-9, "…which is where the hand is");
        ok(! b.p.setDial("nope", 10.0) && b.p.knobValue("nope").empty(), "a knob the pack has not got: refused, and empty");
    }

    group("a linked setting plays its neighbour's weights, fed softer");
    {
        Bench b(rig);
        b.p.setDial("gain", 0.0);
        approx(b.gainAt(1000.0), 0.25 * std::pow(10.0, -6.0 / 20.0), 0.005,
               "at 0 the 60 capture sounds, 6 dB less going in");
        ok(b.p.heldFileId(0) == "g60", "…and it is that file in the slot");
        b.p.setDial("gain", 30.0);
        approx(b.gainAt(1000.0), 0.5 * 0.25 * std::pow(10.0, -6.0 / 20.0) + 0.5 * 0.25, 0.01,
               "at 30, halfway to 60, the same weights softer and louder are crossfaded");
    }

    group("a switch turn is namz::rig's resolve; the dial keeps its angle and follows the new knots");
    {
        Bench b(rig);
        b.p.setDial("gain", 240.0);
        ok(b.p.setSwitch("channel", "red"), "channel to red");
        ok(b.p.settings().at("channel") == "red", "…pinned");
        ok(b.p.settings().at("gain") == "150", "…and the combination is red's only capture");
        approx(b.p.dialDegrees(), 240.0, 1e-9, "…while the dial still stands at 240");
        approx(b.gainAt(1000.0), 0.75 * std::pow(10.0, 4.5 / 20.0), 0.03, "so the red capture plays, 4.5 dB harder in");
        ok(! b.p.setSwitch("channel", "blue"), "a value nothing was captured at is refused");
        ok(b.p.settings().at("channel") == "red", "…and changes nothing");
        ok(b.p.setSwitch("channel", "green"), "back to green");
        approx(b.gainAt(1000.0), 1.0, 0.02, "…and the dial, still at 240, finds its capture again");
        ok(b.p.setSwitch("gain", "60"), "the crossfade dial set by value");
        approx(b.p.dialDegrees(), 60.0, 1e-9, "…is the dial turned there");
    }

    group("tone as bands: the travel law at the reference and at the stops");
    {
        Bench b(rig);
        ok(b.p.bands(1).size() == 1 && b.p.bands(0).empty(), "one band after the model, none before");
        ok(b.p.knobValue("tone") == "150", "the knob starts at its default, the reference");
        const double flatLo = b.gainAt(100.0), flatHi = b.gainAt(10000.0);
        approx(db(flatHi / flatLo), 0.0, 0.1, "at the reference the band is flat");
        ok(b.p.setDial("tone", 300.0), "tone to the plus stop");
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)), shelfDb(6.0, kFs), 0.1, "+6 dB of high shelf at the plus stop, as the formula draws it");
        b.p.setDial("tone", 0.0);
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)), shelfDb(-6.0, kFs), 0.1, "-6 dB at the minus stop");
        b.p.setDial("tone", 225.0);
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)), shelfDb(3.0, kFs), 0.1, "halfway up from the reference: half the gain, in dB");
        ok(b.p.knobValue("tone") == "225", "…and the knob reads back in the pack's words");
    }

    group("tone as a curve: the pack's decibels at this frequency, as a FIR before the model");
    {
        Bench b(rig);
        ok(! b.p.curveActive(0), "at the reference the curve is flat and there is no FIR");
        const double refLo = b.gainAt(60.0), refHi = b.gainAt(10000.0);
        ok(b.p.setDial("bass", 300.0), "bass to the top");
        ok(b.p.curveActive(0) && ! b.p.curveDb(0).empty(), "…and now there is a FIR on the pre side");
        approx(db(b.gainAt(60.0) / refLo), 6.0, 0.6, "+6 dB at 60 Hz");
        approx(db(b.gainAt(10000.0) / refHi), 0.0, 0.3, "nothing at 10 kHz");
        b.p.setDial("bass", 150.0);
        approx(db(b.gainAt(60.0) / refLo), 3.0, 0.6, "halfway: the two curves interpolated, +3 dB");
    }

    group("tone handed in beside the manifest: the same structures, another source");
    {
        Bench b(rig);
        b.p.setDial("bass", 300.0);
        const double refLo = b.gainAt(60.0);                  // the pack's curve: +6 dB at 60 Hz
        // The bench rewrites `bass` as one band — a low shelf reaching +12 dB — and hears it at once.
        namz::rig::Tone bass;
        bass.name = "bass"; bass.sweep = 300; bass.placement = "pre"; bass.reference = "0"; bass.defaultValue = "0";
        namz::rig::Section ls; ls.kind = namz::rig::SectionKind::LowShelf; ls.hz = 200.0; ls.q = 0.7;
        ls.dbAtMin = 0.0; ls.dbAtMax = 12.0;
        bass.sections = { ls };
        b.p.setToneOverride({ bass });
        ok(b.p.toneOverridden(), "the override is in");
        ok(b.p.tones().size() == 2 && b.p.tones()[1].name == "bass" && b.p.tones()[1].positions.empty(),
           "…and the knob plays as the band, in the pack's place for it");
        ok(b.p.knobValue("bass") == "300", "the knob keeps its position across the swap");
        const auto q = felitronics::rigplayer::designSection(felitronics::rigplayer::SectionKind::LowShelf, 200.0, 12.0, 0.7, kFs);
        const double want = felitronics::rigplayer::sectionMagnitudeDb(q, 60.0, kFs);
        approx(db(b.gainAt(60.0) / refLo) + 6.0, want, 0.3, "at 60 Hz the band's own decibels, not the curve's");
        b.p.clearToneOverride();
        ok(! b.p.toneOverridden(), "cleared");
        approx(db(b.gainAt(60.0) / refLo), 0.0, 0.3, "…and the pack's curve is back");

        // A knob the pack has no block for is a new knob, at its default.
        namz::rig::Tone presence;
        presence.name = "presence"; presence.sweep = 300; presence.placement = "post"; presence.reference = "150";
        namz::rig::Section hs; hs.kind = namz::rig::SectionKind::HighShelf; hs.hz = 3000.0; hs.q = 0.7;
        hs.dbAtMin = -6.0; hs.dbAtMax = 6.0;
        presence.sections = { hs };
        b.p.setToneOverride({ presence });
        ok(b.p.tones().size() == 3 && b.p.knobValue("presence") == "150", "a new knob appears, at its reference");
        const double flat = db(b.gainAt(10000.0) / b.gainAt(100.0));
        ok(b.p.setDial("presence", 300.0), "…and turns");
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)) - flat, shelfDb(6.0, kFs), 0.1, "+6 dB of the new shelf");
    }

    // A KNOB THAT CLICKS. Its positions are words with an order and no angle, so there is no rotation
    // for a band's gain to travel on: each position states the filter it IS, and nothing is computed
    // between two of them. Before schema 4 such a knob could not say this at all — it was given an
    // evenly spaced `norm` nobody measured, and a bench that declared it as bands heard its measured
    // curve instead, silently.
    group("a knob that clicks: the filter is stated at the position, whole");
    {
        Bench b(rig);
        namz::rig::Tone edge;
        edge.name = "edge"; edge.placement = "post"; edge.reference = "sharp"; edge.defaultValue = "sharp";
        namz::rig::TonePosition sharp;  sharp.value  = "sharp";      // the anchor states nothing: flat by construction
        namz::rig::TonePosition smooth; smooth.value = "smooth";
        namz::rig::PositionSection hs;
        hs.kind = namz::rig::SectionKind::HighShelf; hs.hz = 3000.0; hs.q = 0.7; hs.gainDb = -6.0;
        smooth.sections = { hs };
        edge.positions = { sharp, smooth };
        b.p.setToneOverride({ edge });
        ok(b.p.knobValue("edge") == "sharp", "the switch opens at its anchor");
        const double flat = db(b.gainAt(10000.0) / b.gainAt(100.0));
        ok(b.p.setSwitch("edge", "smooth") && b.p.knobValue("edge") == "smooth", "…and clicks over");
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)) - flat, shelfDb(-6.0, kFs), 0.1,
               "what sounds is that position's own filter, at its own decibels");
        ok(b.p.setSwitch("edge", "sharp"), "back to the anchor");
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)) - flat, 0.0, 0.1,
               "…where the models play exactly as they were captured");
        ok(! b.p.setSwitch("edge", "smoothh"), "a position this knob does not declare is REFUSED");
        ok(b.p.knobValue("edge") == "sharp",
           "…and the knob stayed where it was: it used to accept the word, read it back, and play the anchor");
    }

    // LAW 8 (state, not timing). The blend ramps are `end = want + (current-want)*decay` — asymptotic, so
    // without a snap they never arrive. An EXACT-zero target is not hypothetical and needs no special pack:
    // BlendKnob::linOf returns 0.0 for any level at or below -120 dB, which is how this very rig spells
    // "off" (dry end: wetDb -120). Move the dial there after an audible position and the wet gain decays
    // to a subnormal fixed point (decay 0.587 at kBlock 256 => k <= 0.5/(1-decay) ~= 1.2, so 1 ulp) and
    // stays, while `mixDry` — gated on a dry IR being LOADED, not on the gains — keeps the per-sample mix
    // loop running over it for the life of the rig.
    //
    // Asserted twice, as in the Saturator's law-8 test: at 5 s the un-flushed value is 1 ulp (subnormal,
    // so a machine with hardware FTZ could read it as zero and hide a regression), while at 0.8 s it is
    // 1.8e-35 — a NORMAL float — and the flush has already fired (it does so at 0.69 s). The observable is
    // liveWet(), the applied gain the audio thread actually multiplied by.
    group("blend: law 8 — a gain aimed at exact zero ARRIVES (it used to park in the subnormals)");
    {
        Bench b(rig);
        b.rms(440.0, 0.2, 16, 0);                                   // default dial = wet end: curWet_ -> 1
        ok(b.p.liveWet() > 0.9f, "the wet gain is up before the move");
        ok(b.p.setDial("mix", 0.0), "mix to the dry end — the pack states the wet path at -120 dB, i.e. 0");
        b.rms(440.0, 0.2, 150, 0);                                  // 0.8 s: un-flushed would be 1.8e-35
        ok(b.p.liveWet() == 0.0f, "applied wet is EXACTLY 0 at 0.8 s (a normal 1.8e-35 without the snap)");
        b.rms(440.0, 0.2, 788, 0);                                  // out to 5 s
        ok(b.p.liveWet() == 0.0f, "still exactly 0 after 5 s");
    }

    group("blend: the dry path and the wet one as the pack states them");
    {
        Bench b(rig);
        approx(b.gainAt(1000.0), 0.5, 0.01, "wet end (the default): the model alone");
        ok(b.p.setDial("mix", 0.0), "mix to the dry end");
        approx(b.gainAt(1000.0), std::pow(10.0, -6.0 / 20.0), 0.01, "dry end: the DI through the dry path, 6 dB down");
        b.p.setDial("mix", 150.0);
        approx(b.gainAt(1000.0), 0.5 * std::pow(10.0, -6.0 / 20.0) + 0.5 * 0.5, 0.01,
               "halfway: half of each, in amplitude, summing in phase");
    }
    {
        Bench b(throughTheFormat(testRig(-1)));
        b.p.setSwitch("gain", "240");                          // a unity capture, so the two paths match
        b.p.setDial("mix", 150.0);
        const double half = b.gainAt(1000.0);
        ok(half < 0.5 * 0.501 + 0.5 - 0.4, "polarity -1: the dry path is subtracted, and the middle of the knob nearly nulls");
        approx(half, 0.5 - 0.5 * std::pow(10.0, -6.0 / 20.0), 0.01, "…to exactly the difference of the two");
    }

    group("the pack's own levels: the guitar in, the device out");
    {
        const double q6 = std::pow(10.0, -6.0 / 20.0);
        Bench b(levelled(-6.0, 0.0));
        ok(b.p.stageInputDb() == -6.0 && b.p.stageOutputDb() == 0.0, "the player reads both levels off the pack");
        approx(b.gainAt(1000.0), 0.5 * q6, 0.01, "at 150 the model is fed 6 dB softer");
        b.p.setDial("gain", 0.0);
        approx(b.gainAt(1000.0), 0.25 * q6 * q6, 0.005,
               "at the link the two ADD: the pack's 6 dB and the alias's own 6 dB");
    }
    {
        // The input level reaches the DRY side of a blend too. One guitar cannot arrive at the two ends
        // of a mix at two different levels, so both move together and the stated mix is untouched.
        const double q6 = std::pow(10.0, -6.0 / 20.0);
        Bench plain(rig), quiet(levelled(-6.0, 0.0));
        ok(plain.p.setDial("mix", 0.0), "the plain pack to the dry end");
        ok(quiet.p.setDial("mix", 0.0), "the levelled one too");
        approx(quiet.gainAt(1000.0), plain.gainAt(1000.0) * q6, 0.01, "the dry path is fed 6 dB softer as well");
        plain.p.setDial("mix", 150.0); quiet.p.setDial("mix", 150.0);
        approx(quiet.gainAt(1000.0), plain.gainAt(1000.0) * q6, 0.01,
               "…so halfway across the knob the mix is the same mix, 6 dB down");
    }
    {
        // The output level is a scalar on the whole stage, applied AFTER the mix — the same number at
        // every position of the blend knob, which is what "it cannot move the blend" means as a number.
        const double q6 = std::pow(10.0, -6.0 / 20.0);
        Bench plain(rig), quieter(levelled(0.0, -6.0));
        approx(quieter.gainAt(1000.0), plain.gainAt(1000.0) * q6, 0.005, "the wet end comes out 6 dB down");
        ok(plain.p.setDial("mix", 0.0), "the plain pack to the dry end");
        ok(quieter.p.setDial("mix", 0.0), "the levelled one too");
        approx(quieter.gainAt(1000.0), plain.gainAt(1000.0) * q6, 0.005, "so does the dry end — the same scalar");
        plain.p.setDial("mix", 150.0); quieter.p.setDial("mix", 150.0);
        approx(quieter.gainAt(1000.0), plain.gainAt(1000.0) * q6, 0.005,
               "and so does halfway across, which is the blend's ratio left exactly where the pack put it");
    }
    {
        // A pack level cannot put a step into a crossfade that had none: the ladder with one is the
        // ladder without one times the scalar, at every angle — the link at the bottom included.
        const double q3 = std::pow(10.0, -3.0 / 20.0);
        Bench plain(rig), lower(levelled(0.0, -3.0));
        for (const double deg : { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 }) {
            plain.p.setDial("gain", deg); lower.p.setDial("gain", deg);
            char what[96];
            std::snprintf(what, sizeof what, "the crossfade at %.0f deg is untouched but for the scalar", deg);
            approx(lower.gainAt(1000.0), plain.gainAt(1000.0) * q3, 0.01, what);
        }
    }
    {
        // The loudness tag is a contract, not a listener's option. A player that opens un-normalised
        // plays every capture at whatever level the hardware gave, and this default must not be
        // quietly turned back — hence a test on the default itself.
        RigPlayer fresh;
        ok(fresh.normalize(), "a player normalizes by default");
        Bench b(rig);
        ok(b.p.normalize(), "…and so does one with a pack in it");
        b.p.setNormalize(false);
        b.p.unload();
        b.load(rig);
        ok(! b.p.normalize(), "the switch is still the host's to throw, and survives a reload");
    }

    group("input is drive and output is volume: a model that is not a pure scalar tells them apart");
    {
        // The mean of the output IS the offset the network adds. Scale what goes IN and the offset
        // stands; scale what comes OUT and the offset scales too. Every other fixture in this file is
        // a pure gain, through which the two levels are indistinguishable — swap where the player
        // applies them and nothing else here would notice.
        const auto biased = [] (const std::string& nam) {
            return [nam] (const std::string&) { return bytesOf(nam); };
        };
        const auto meanOf = [] (Bench& b) {
            std::vector<float> l((std::size_t) kBlock);
            double sum = 0.0; long n = 0;
            for (int k = 0; k < 80; ++k) {
                std::fill(l.begin(), l.end(), 0.0f);          // silence in: what comes out is the offset
                float* io[2] { l.data(), nullptr };
                b.p.process(io, 1, kBlock);
                b.p.serviceHere();
                if (k >= 48) for (const float v : l) { sum += (double) v; ++n; }
            }
            return sum / (double) std::max(1L, n);
        };
        const auto model = biasModel(0.5, 0.25);

        Bench plain(rig);
        plain.p.load(throughTheFormat(testRig()).chain[0], biased(model));
        const double base = meanOf(plain);
        ok(std::abs(base - 0.25) < 0.02, "the fixture adds its offset: mean 0.25");

        auto inRig = testRig(); inRig.chain[0].inputDb = -6.0;
        Bench fedLess(rig);
        fedLess.p.load(throughTheFormat(inRig).chain[0], biased(model));
        ok(std::abs(meanOf(fedLess) - 0.25) < 0.02,
           "input_db does NOT touch the offset — it is drive, and drive is not volume");

        auto outRig = testRig(); outRig.chain[0].outputDb = -6.0;
        Bench quieter(rig);
        quieter.p.load(throughTheFormat(outRig).chain[0], biased(model));
        ok(std::abs(meanOf(quieter) - 0.25 * std::pow(10.0, -6.0 / 20.0)) < 0.02,
           "…while output_db takes the offset down with everything else: it is volume");
    }

    group("the host's own hand: the WHOLE number, and the player does the subtracting");
    {
        // A host with a fader states what it wants the device played at. It never sends a difference:
        // the difference needs to know which pack is loaded, the host reads that from its own document,
        // and the document and the pack disagree whenever an edit or a rebuild is in flight. Here there
        // is one number and one place that applies it.
        const double q6 = std::pow(10.0, -6.0 / 20.0);
        Bench b(levelled(-6.0, 0.0));
        ok(! b.p.hostInputDb() && ! b.p.hostOutputDb(), "a fresh player holds no hand: the pack's own stands");
        approx(b.gainAt(1000.0), 0.5 * q6, 0.01, "…and it is the pack's -6 dB that is heard");

        b.p.setHostInputDb(0.0);
        approx(b.gainAt(1000.0), 0.5, 0.01, "the hand at 0 REPLACES the pack's -6, it does not add to it");

        b.p.setHostInputDb(-12.0);
        approx(b.gainAt(1000.0), 0.5 * std::pow(10.0, -12.0 / 20.0), 0.01, "…and -12 is heard as -12, once");

        // THE WHOLE POINT: a new pack states something else, and the hand still wins without the host
        // saying a word. This is the case that was wrong when the host did the subtracting — the pack
        // changed under it and its arithmetic did not.
        b.load(levelled(-18.0, 0.0));
        approx(b.gainAt(1000.0), 0.5 * std::pow(10.0, -12.0 / 20.0), 0.01,
               "a pack swap does not move the hand, and does not double with it");

        b.p.setHostInputDb(std::nullopt);
        approx(b.gainAt(1000.0), 0.5 * std::pow(10.0, -18.0 / 20.0), 0.01,
               "letting go hands the level back to the pack");

        b.p.setHostOutputDb(-6.0);
        approx(b.gainAt(1000.0), 0.5 * std::pow(10.0, -18.0 / 20.0) * q6, 0.01, "the other end works the same way");
    }

    group("alignment: two captures that land apart are lined up before they mix");
    {
        // Two files: unity, and unity two samples late.
        namz::rig::Rig r;
        namz::rig::Stage st; st.kind = namz::rig::StageKind::Nam; st.rawKind = "nam";
        namz::rig::Control g; g.name = "gain"; g.role = namz::rig::Role::Gain; g.values = { "60", "240" }; g.sweep = 300;
        st.device.controls = { g };
        namz::rig::FileEntry a; a.id = "early"; a.settings = { { "gain", "60" } };
        namz::rig::FileEntry c; c.id = "late";  c.settings = { { "gain", "240" } };
        st.device.files = { a, c };
        r.chain = { st };
        std::map<std::string, std::vector<std::byte>> files {
            { "early", bytesOf(gainModel(1.0)) }, { "late", bytesOf(delayModel(2)) } };
        const ModelSource src = [&files](const std::string& id) { return files.at(id); };

        const auto table = measureAlignment(st.device, src, kFs);
        ok(table.lagByFile.size() == 2, "both files measured");
        ok(table.lagByFile.count("late") && table.lagByFile.at("late") == 2, "the late one reads two samples late");
        ok(table.delayOf("early") == 2 && table.delayOf("late") == 0, "so the early one is delayed by two, the late one not at all");

        // At 6 kHz two samples are a quarter turn: unaligned, a 50/50 sum of the two is 3 dB down.
        {
            Bench b(r);
            b.files = files;
            b.p.setDial("gain", 150.0);
            approx(b.gainAt(6000.0), std::sqrt(0.5), 0.02, "without the table the pair combs: 0.707 at 6 kHz");
        }
        {
            Bench b(r);
            b.files = files;
            b.p.setAlignment(table);                           // before anything lands: it travels with the loads
            b.p.setDial("gain", 150.0);
            approx(b.gainAt(6000.0), 1.0, 0.02, "with the table in hand before playing, the pair sums to one");
            ok(b.p.appliedSlotDelay(0) == 2 && b.p.appliedSlotDelay(1) == 0, "…the early slot delayed, the late one not");
        }
        {
            // A table that arrives MID-MIX lands at the one instant a slot is free — weight exactly zero —
            // never as a splice on a live signal. So the comb stands until the dial visits a knot, and is
            // gone once it has: the host that can measure before playing should.
            Bench b(r);
            b.files = files;
            b.p.setDial("gain", 150.0);
            approx(b.gainAt(6000.0), std::sqrt(0.5), 0.02, "mid-mix, before the table: the comb");
            b.p.setAlignment(table);
            approx(b.gainAt(6000.0), std::sqrt(0.5), 0.02, "…and still the comb: neither slot is silent, so nothing lands");
            b.p.setDial("gain", 240.0);
            b.gainAt(6000.0);                                  // the top: the early capture leaves its slot at zero
            b.p.setDial("gain", 150.0);
            approx(b.gainAt(6000.0), 1.0, 0.02, "back at the middle it lands again, delayed, and the pair sums to one");
        }
    }

    group("alignment from the pack: the lags written at pack time, no probe at load");
    {
        namz::rig::Rig r;
        namz::rig::Stage st; st.kind = namz::rig::StageKind::Nam; st.rawKind = "nam";
        namz::rig::Control g; g.name = "gain"; g.role = namz::rig::Role::Gain; g.values = { "60", "240" }; g.sweep = 300;
        st.device.controls = { g };
        namz::rig::FileEntry a; a.id = "early"; a.settings = { { "gain", "60" } };  a.lagSamples = 0;
        namz::rig::FileEntry c; c.id = "late";  c.settings = { { "gain", "240" } }; c.lagSamples = 2;
        st.device.files = { a, c };
        r.chain = { st };
        bool okM = false;
        const auto packed = throughTheFormat(r, &okM);
        ok(okM && packed.chain.size() == 1 && packed.chain[0].device.files.size() == 2
           && packed.chain[0].device.files[1].lagSamples == 2, "lag_samples round-trips through the writer and the reader");
        const auto table = AlignmentTable::fromDevice(packed.chain[0].device);
        ok(table.fromPack && ! table.empty() && table.delayOf("early") == 2 && table.delayOf("late") == 0,
           "the table comes from the pack: the early file delayed by two");
        approx(table.delayOf("early", 96000.0, 48000.0), 4.0, 0.0, "…and in host samples at 96 kHz, four");

        std::map<std::string, std::vector<std::byte>> files {
            { "early", bytesOf(gainModel(1.0)) }, { "late", bytesOf(delayModel(2)) } };
        Bench b(packed);
        b.files = files;
        ok(b.p.alignmentFromPack(), "the player took the pack's reading at load");
        b.p.setDial("gain", 150.0);
        approx(b.gainAt(6000.0), 1.0, 0.02, "no setAlignment, no probe: the pair sums to one from the pack's numbers");
        ok(b.p.appliedSlotDelay(0) == 2 && b.p.appliedSlotDelay(1) == 0, "…the early slot delayed by the pack's two");

        // Half a reading is no reading.
        auto half = packed;
        half.chain[0].device.files[1].lagSamples.reset();
        ok(AlignmentTable::fromDevice(half.chain[0].device).empty(), "a stage with one entry unmeasured is not measured");
        Bench h(half);
        h.files = files;
        ok(! h.p.alignmentFromPack(), "…and the player does not pretend it is");
    }

    group("stereo: each channel its own signal through the same decision");
    {
        Bench b(rig, 2);
        const double l = b.rms(1000.0, 0.1, 48, 32, 0, 0.2);
        const double rr = b.rms(1000.0, 0.1, 48, 32, 1, 0.2);
        approx(l / (0.1 / std::sqrt(2.0)), 0.5, 0.01, "left: the 0.5 capture on its own signal");
        approx(rr / (0.2 / std::sqrt(2.0)), 0.5, 0.01, "right: the same capture on the other");
    }

    group("another rate: the bands and the curves are designed for it");
    {
        Bench b(rig, 1, 96000.0);
        // A 48 kHz model at a 96 kHz host is rate-matched by the stage, and its resampler is not
        // perfectly flat to 10 kHz — so the shelf is read against what the chain does with it flat.
        const double base = db(b.gainAt(10000.0) / b.gainAt(100.0));
        b.p.setDial("tone", 300.0);
        approx(db(b.gainAt(10000.0) / b.gainAt(100.0)) - base, shelfDb(6.0, 96000.0), 0.1,
               "the shelf at 96 kHz, as the formula draws it there");
        b.p.setDial("tone", 150.0);
        const double refLo = b.gainAt(60.0, 0.1, 48, 64);
        b.p.setDial("bass", 300.0);
        approx(db(b.gainAt(60.0, 0.1, 48, 64) / refLo), 6.0, 0.6, "+6 dB of curve at 96 kHz");
    }

    group("a dial at rest: after kColdAfterSeconds the silent neighbour is not run — and stays loaded");
    {
        // THE ECONOMY, measured in the plugin: a dial parked on a capture kept its neighbour's network
        // running at weight zero, every block, for nothing — 4.5 % of a P-core, 14 % of an E-core.
        // With the handover a step (the plugin's STEP), the hand is always on a capture and the
        // neighbour is always at exactly zero: after the rest it sleeps, and one model plays.
        Bench b(rig);
        b.p.setBlendShape({ 0.5, 0.0 });
        ok(b.p.coldAfterSeconds() == RigPlayer::kColdAfterSeconds && RigPlayer::kColdAfterSeconds == 2.0,
           "the rest is the player's constant, two seconds");
        approx(b.gainAt(1000.0), 0.5, 0.01, "150: the 0.5 capture alone");
        ok(b.p.heldFileId(1) == "g240" && ! b.p.slotCold(1), "…its neighbour held at zero, awake: the hand only just arrived");
        ok(b.p.modelLoads() == 2, "two loads");
        b.p.clearCounters();
        const int rest = (int) std::ceil(2.0 * kFs / kBlock);
        approx(b.rms(1000.0, 0.1, rest, 32) / (0.1 / std::sqrt(2.0)), 0.5, 0.01, "two seconds later it sounds exactly the same");
        ok(b.p.slotCold(1), "…and the neighbour's slot is cold");
        ok(! b.p.slotCold(0), "the sounding slot is not");
        ok(b.p.heldFileId(1) == "g240", "the model is still in the cold slot");
        ok(b.p.modelLoads() == 2, "…nothing was loaded, nothing unloaded");
        // The rest began when the neighbour landed and settled — somewhere inside the first
        // measurement's 80 blocks — so of the 80 + 375 + 32 blocks since, it slept between 32 and 112.
        const int slept = b.p.coldBlocks(1);
        ok(slept >= 32 && slept <= 80 + 32, "it fell asleep at the two-second line, not before (" + std::to_string(slept)
                                            + " blocks slept of " + std::to_string(80 + rest + 32) + ")");
        ok(b.p.warmBlocks() == 0, "a sleeping slot is not a warming one");
        b.rms(1000.0, 0.1, 0, 20);
        ok(b.p.coldBlocks(1) == slept + 20, "twenty more blocks, twenty more slept: the model was not run once");
        ok(b.p.coldBlocks(0) == 0, "…and the sounding one never slept");

        // THE TURN AFTER THE REST. The hand moves to the next capture; in STEP that is all of the
        // neighbour at once. The cold slot wakes on the first block, is fed its field — a warm-up, not
        // a load — and only then does the weight travel, at the law's own pace, never in a step.
        b.p.clearCounters();
        ok(b.p.setDial("gain", 200.0), "the hand moves past the midpoint: the 240 capture, all of it");
        int firstMove = -1, arrived = -1;
        bool awakeAtOnce = false;
        std::vector<double> level;
        for (int k = 0; k < 40; ++k) {
            level.push_back(b.rms(1000.0, 0.1, 0, 1) / (0.1 / std::sqrt(2.0)));
            if (k == 0) awakeAtOnce = ! b.p.slotCold(1);
            const float m = b.p.liveMix();
            if (firstMove < 0 && m > 0.0f) firstMove = k;
            if (arrived < 0 && m >= 1.0f) arrived = k;
        }
        ok(awakeAtOnce, "the first block of the turn wakes the slot");
        // A Linear model DECLARES no field (NAM answers prewarm only for convnet and lstm; WaveNet's is
        // read from its config), so its need is zero and the law may move the weight in the block the
        // slot woke in. A model with a field waits it out first — the law's own test proves that with
        // a 6332-sample field; here the claim is only "no later than the warm-up plus one block".
        ok(firstMove >= 0 && firstMove <= 2, "the weight starts moving no later than the warm-up plus one block (block "
                                             + std::to_string(firstMove) + ")");
        ok(arrived == firstMove + 3, "…and arrives four blocks later, a quarter per block: the law's own slew, no step");
        ok(b.p.biggestJump() <= 0.25f + 1e-6f, "no block moved the weight more than the law allows");
        ok(b.p.modelLoads() == 2, "no load in the wake: the model never left");
        ok(b.p.heldFileId(1) == "g240" && ! b.p.slotCold(1), "…the same model, awake");
        // One block of a 1 kHz sine is 5⅓ cycles, so a block's RMS wanders by a percent or two with
        // its phase; a step would be the whole half of level in one block.
        bool monotone = true; double worstRise = 0.0;
        for (std::size_t k = 1; k < level.size(); ++k) {
            monotone = monotone && level[k] >= level[k - 1] - 0.03;
            worstRise = std::max(worstRise, level[k] - level[k - 1]);
        }
        approx(level.back(), 1.0, 0.02, "what sounds at the end is the 240 capture");
        ok(monotone && worstRise <= 0.5 * 0.25 * 1.3, "…reached by a rise, block on block, none of them a step (the worst "
                                                      + std::to_string(worstRise) + " of level in one block)");
        ok(b.p.warmBlocks() >= 0 && b.p.warmBlocks() <= 2, "the wake cost at most two blocks of warming");

        // Zero means never: a host that wants both networks running says so.
        b.p.setColdAfterSeconds(0.0);
        b.p.clearCounters();
        b.p.setDial("gain", 150.0);
        b.rms(1000.0, 0.1, rest + 40, 32);
        ok(! b.p.slotCold(0) && ! b.p.slotCold(1) && b.p.coldBlocks(0) == 0 && b.p.coldBlocks(1) == 0,
           "with the rest set to zero nothing sleeps, however long the hand rests");
        b.p.setColdAfterSeconds(0.5);
        b.rms(1000.0, 0.1, (int) std::ceil(0.5 * kFs / kBlock), 32);
        ok(b.p.slotCold(1), "…and half a second, once set, is a rest");
    }

    group("a player that sleeps sounds bit-identically to one that never does — through rests and turns");
    {
        // The same sine into two players, block for block: one sleeps after two seconds, one never.
        // A cold slot is left out of the mix, not multiplied by zero, and a wake re-lands a model that
        // has no memory of its own (a Linear gain, no alignment delay) — so not one sample may differ,
        // resting on the lower capture (slot 1 asleep) or on the upper (slot 0 — the live buffer —
        // asleep), nor across the turns between them.
        Bench sleeps(rig), never(rig);
        sleeps.p.setBlendShape({ 0.5, 0.0 }); never.p.setBlendShape({ 0.5, 0.0 });
        never.p.setColdAfterSeconds(0.0);
        const int rest = (int) std::ceil(2.0 * kFs / kBlock);
        std::vector<float> x((std::size_t) kBlock), y((std::size_t) kBlock);
        double phase = 0.0;
        const auto both = [&](int blocks) {
            float worst = 0.0f;
            for (int b = 0; b < blocks; ++b) {
                for (int i = 0; i < kBlock; ++i) {
                    x[(std::size_t) i] = y[(std::size_t) i] = (float) (0.1 * std::sin(phase));
                    phase += 2.0 * 3.14159265358979323846 * 1000.0 / kFs;
                }
                float* ix[1] { x.data() }; float* iy[1] { y.data() };
                sleeps.p.process(ix, 1, kBlock); sleeps.p.serviceHere();
                never.p.process(iy, 1, kBlock);  never.p.serviceHere();
                for (int i = 0; i < kBlock; ++i) worst = std::max(worst, std::abs(x[(std::size_t) i] - y[(std::size_t) i]));
            }
            return worst;
        };
        ok(both(rest + 60) == 0.0f, "two and a half seconds on 150: identical, and by then slot 1 sleeps in one player");
        ok(sleeps.p.slotCold(1) && ! never.p.slotCold(1), "…which it does");
        sleeps.p.setDial("gain", 200.0); never.p.setDial("gain", 200.0);
        ok(both(rest + 60) == 0.0f, "the turn to 240 and two and a half seconds there: identical, slot 0 now asleep");
        ok(sleeps.p.slotCold(0) && ! sleeps.p.slotCold(1), "…the live buffer's slot, and the other awake");
        sleeps.p.setDial("gain", 150.0); never.p.setDial("gain", 150.0);
        ok(both(60) == 0.0f, "…and the turn back: identical to the last sample");
        ok(sleeps.p.modelLoads() == 2 && never.p.modelLoads() == 2, "neither player loaded anything for a turn");
    }

    group("a slot's delay line is cleared on the way to sleep, as it is on a landing");
    {
        // THE BURST THIS CATCHES: with an alignment delay on the sleeping slot, its delay line held the
        // last hundred samples from before the rest; a model that declares no field is heard on the
        // block it wakes in, and those samples came out first — two seconds old, on a silent input.
        AlignmentTable table;
        table.lagByFile = { { "g60", 100 }, { "g150", 100 }, { "g240", 0 }, { "r150", 100 } };   // 240 is early: delayed by 100
        Bench b(rig);
        b.p.setBlendShape({ 0.5, 0.0 });
        b.p.setAlignment(table);
        const int rest = (int) std::ceil(2.0 * kFs / kBlock);
        b.rms(1000.0, 0.1, rest + 40, 8);
        ok(b.p.slotCold(1) && b.p.appliedSlotDelay(1) == 100, "the delayed neighbour is asleep");
        std::vector<float> z((std::size_t) kBlock);
        float* io[1] { z.data() };
        float peak = 0.0f;
        for (int k = 0; k < 8; ++k) { std::fill(z.begin(), z.end(), 0.0f); b.p.process(io, 1, kBlock); b.p.serviceHere(); }
        for (const float v : z) peak = std::max(peak, std::abs(v));
        ok(peak == 0.0f, "silence in, silence out, before the turn");
        b.p.setDial("gain", 200.0);
        for (int k = 0; k < 8; ++k) {
            std::fill(z.begin(), z.end(), 0.0f); b.p.process(io, 1, kBlock); b.p.serviceHere();
            for (const float v : z) peak = std::max(peak, std::abs(v));
        }
        ok(peak == 0.0f, "…and after it: nothing of the rest comes back out (peak " + std::to_string(peak) + ")");
        ok(! b.p.slotCold(1) && b.p.liveMix() >= 1.0f, "the slot woke and took the sound — of a silent input");
    }

    group("between two captures both models run, however long the hand rests");
    {
        // SMOOTH between two captures: both are heard, so neither can be cold — two passes, as it
        // must be. The saving is for the dial at rest on a capture, never for a mix.
        Bench b(rig);
        b.p.setDial("gain", 200.0);
        b.p.clearCounters();
        const int rest = (int) std::ceil(2.0 * kFs / kBlock);
        approx(b.rms(1000.0, 0.1, rest + 40, 32) / (0.1 / std::sqrt(2.0)), (40.0 / 90.0) * 0.5 + (50.0 / 90.0) * 1.0, 0.015,
               "two and a half seconds at 200: still 4/9 and 5/9 of each");
        ok(! b.p.slotCold(0) && ! b.p.slotCold(1), "neither slot is cold");
        ok(b.p.coldBlocks(0) == 0 && b.p.coldBlocks(1) == 0, "…and neither slept a single block");
    }

    group("unload: silence, and a second device loads clean");
    {
        Bench b(rig);
        approx(b.gainAt(1000.0), 0.5, 0.01, "sounding");
        b.p.unload();
        ok(! b.p.loaded() && b.p.settings().empty(), "nothing loaded");
        const double after = b.rms(1000.0, 0.1, 8, 8);
        ok(after < 1e-6, "…and nothing sounds");
        b.load(rig);
        approx(b.gainAt(1000.0), 0.5, 0.01, "loaded again, sounding again");
    }

    return felitronics::test::report();
}
