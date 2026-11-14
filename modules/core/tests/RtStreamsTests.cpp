// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// RtStreams under REAL threads: an audio thread opening Blocks against a controller publishing, switching,
// starting and harvesting. This is the suite ThreadSanitizer gates (the CI tsan job runs it with nothing
// suppressed): every sample, position and length behind the gate is PLAIN memory, so a missing
// happens-before edge in the gate is a data race TSan reports. The protocol's logic at every interleaving is
// the interleaving suite's job (RtStreamsInterleaveTests.cpp); here the checks are what a listener would hear.
//
// THE TEST'S OWN COORDINATION IS RELAXED. The hold flags, counters and stop flag below are relaxed atomics,
// which create no happens-before edge, so they cannot supply one the gate is missing. Every content check
// rests on the gate alone. The one exception is std::thread::join at the end of each section, after which
// only counters are read.
//
// Each section FORCES contention, so the refusal path runs on every row, not when the scheduler happens to
// cooperate: the audio thread is asked to hold a Block open, the controller tries while it holds (and must
// be refused), then lets it go and retries (and must land).

#include <felitronics/core/RtStreams.h>

// Samples here are MARKERS (clip ids, block numbers) copied, never computed: exact equality is the property.
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#include "felitronics_test.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

using felitronics::test::group;
using felitronics::test::ok;

namespace {

constexpr auto rlx = std::memory_order_relaxed;
using Clock = std::chrono::steady_clock;

// Spins until `done()` or the deadline; true when done.
template <class F> bool waitFor (F&& done, double seconds = 20.0) {
    const auto end = Clock::now() + std::chrono::duration_cast<Clock::duration> (std::chrono::duration<double> (seconds));
    while (! done()) {
        if (Clock::now() > end) return false;
        std::this_thread::yield();
    }
    return true;
}

// The hold handshake: the controller asks for round r; the audio thread, inside an admitted Block, reports
// it holds r and spins until released. Relaxed throughout (see the header).
struct Hold {
    std::atomic<long> request { 0 }, holding { 0 }, release { 0 };
    void insideBlock() {
        const long r = request.load (rlx);
        if (r == 0 || holding.load (rlx) == r) return;
        holding.store (r, rlx);
        (void) waitFor ([&] { return release.load (rlx) == r; });
    }
};

int lenFor (long id) { return (int) (id % 200) + 1; }

// ---- AuditionStream -----------------------------------------------------------------------------------
void audition()
{
    group ("AuditionStream: publishes over a playing clip (looping, then one-shot)");
    felitronics::core::AuditionStream a;
    ok (a.tryReserve (256), "reserve");
    a.loop.store (true, rlx);

    std::atomic<bool> go { true };
    std::atomic<long> blocks { 0 }, torn { 0 }, midway { 0 }, completed { 0 };
    Hold hold;
    std::thread reader ([&] {
        float last = -1.0f;
        while (go.load (rlx)) {
            auto b = a.beginBlock();
            if (! b) { std::this_thread::yield(); continue; }
            blocks.fetch_add (1, rlx);
            const float m = b.data()[0];
            if (b.len() != lenFor ((long) m)) torn.fetch_add (1, rlx);       // the clip's own length, never another's
            if (m != last && b.pos != 0) midway.fetch_add (1, rlx);          // a clip is first seen from its start
            last = m;
            for (int i = 0; i < 16; ++i) {
                if (b.pos >= b.len()) {
                    if (a.loop.load (rlx)) b.pos = 0;
                    else { b.finish(); completed.store ((long) m, rlx); break; }
                }
                if (b.data()[b.pos++] != m) torn.fetch_add (1, rlx);
            }
            hold.insideBlock();
        }
    });

    std::vector<float> clip (256);
    auto fill = [&] (long id) { for (int i = 0; i < lenFor (id); ++i) clip[(std::size_t) i] = (float) id; };
    long id = 0, refused = 0, forcedRounds = 0, forcedRefusals = 0, landedWhileHeld = 0, holdTimeouts = 0;

    // Phase 1, looping: a clip always plays, so every publish lands over a live reader.
    for (int k = 0; k < 4000; ++k) {
        fill (++id);
        if (k % 8 == 7) {                                                    // forced contention
            ++forcedRounds;
            hold.request.store (id, rlx);
            if (! waitFor ([&] { return hold.holding.load (rlx) == id; })) { ++holdTimeouts; hold.release.store (id, rlx); continue; }
            if (a.tryPublish (clip.data(), (std::size_t) lenFor (id))) ++landedWhileHeld; else ++forcedRefusals;
            hold.release.store (id, rlx);
        }
        if (k % 50 == 49) while (! a.tryReserve (256)) ++refused;            // resize under a live reader
        while (! a.tryPublish (clip.data(), (std::size_t) lenFor (id))) ++refused;
    }
    hold.request.store (0, rlx);

    // Phase 2, one-shot: every clip that lands must play to its END — a block of the previous clip that ends
    // it must not take the new one down behind its back (the reader's own flag-down vs the re-arm).
    a.loop.store (false, rlx);
    long lost = 0;
    for (int k = 0; k < 1500; ++k) {
        fill (++id);
        while (! a.tryPublish (clip.data(), (std::size_t) lenFor (id))) ++refused;
        if (k % 3 == 0) continue;                                            // publish over a clip still playing
        if (! waitFor ([&] { return completed.load (rlx) == id; }, 10.0)) ++lost;
    }
    go.store (false, rlx);
    reader.join();

    ok (blocks.load() > 0, "the audio thread was admitted to blocks");
    ok (torn.load() == 0, "no block saw a clip that was half one render and half another (" + std::to_string (torn.load()) + ")");
    ok (midway.load() == 0, "every clip was first heard from its start (" + std::to_string (midway.load()) + ")");
    ok (holdTimeouts == 0, "every forced hold engaged");
    ok (landedWhileHeld == 0, "no publish landed while a block was held open");
    ok (forcedRefusals == forcedRounds && forcedRounds == 500, "every forced round was refused, then landed on the retry ("
                                                                + std::to_string (forcedRefusals) + "/" + std::to_string (forcedRounds) + ")");
    ok (lost == 0, "every one-shot clip that landed played to its end (" + std::to_string (lost) + " lost)");
    std::printf ("    audition: %ld blocks, %ld unforced refusals\n", blocks.load(), refused);
}

// ---- ConvStream ---------------------------------------------------------------------------------------
void conv()
{
    group ("ConvStream: refill + clip, then live, racing the clip's own end");
    using C = felitronics::core::ConvStream;
    C c;
    ok (c.tryReserve (64), "reserve");

    std::atomic<bool> go { true };
    std::atomic<long> torn { 0 }, liveBlocks { 0 };
    Hold hold;
    std::thread reader ([&] {
        while (go.load (rlx)) {
            auto b = c.beginBlock();
            if (! b) { std::this_thread::yield(); continue; }
            if (b.mode() == C::kClip) {
                const float m = b.len() > 0 ? b.data()[0] : 0.0f;
                if (b.len() > 0 && b.len() != lenFor ((long) m) % 64 + 1) torn.fetch_add (1, rlx);
                for (int i = 0; i < 4 && b.pos < b.len(); ++i) if (b.data()[b.pos++] != m) torn.fetch_add (1, rlx);
                if (b.pos >= b.len()) b.finish();                            // the clip ends: one-shot
            } else {
                liveBlocks.fetch_add (1, rlx);
            }
            hold.insideBlock();
        }
    });

    std::vector<float> clip (64);
    long refused = 0, lostLive = 0, forcedRounds = 0, forcedRefusals = 0, landedWhileHeld = 0, notGuaranteed = 0, holdTimeouts = 0;
    for (long id = 1; id <= 3000; ++id) {
        const int n = lenFor (id) % 64 + 1;
        for (int i = 0; i < n; ++i) clip[(std::size_t) i] = (float) id;
        if (id % 50 == 0) while (! c.tryReserve (64)) ++refused;             // resize under a live reader
        while (! c.tryRefill (clip.data(), (std::size_t) n)) ++refused;
        if (! c.trySetMode (C::kClip)) ++notGuaranteed;                      // right after a true refill: never refused
        if (id % 10 == 0) {                                                  // forced contention on the switch
            ++forcedRounds;
            hold.request.store (id, rlx);
            if (! waitFor ([&] { return hold.holding.load (rlx) == id || c.mode() == 0; })) { ++holdTimeouts; }
            if (hold.holding.load (rlx) == id) { if (! c.trySetMode (C::kLive)) ++forcedRefusals; else ++landedWhileHeld; }
            else --forcedRounds;                                             // the clip ended before a block held it
            hold.release.store (id, rlx);
        }
        while (! c.trySetMode (C::kLive)) ++refused;
        const long seen = liveBlocks.load (rlx);
        (void) waitFor ([&] { return liveBlocks.load (rlx) >= seen + 2; }, 10.0);
        if (c.mode() != C::kLive) ++lostLive;                                // the old clip's end took live down
    }
    go.store (false, rlx);
    reader.join();

    ok (torn.load() == 0, "no clip block saw two clips (" + std::to_string (torn.load()) + ")");
    ok (notGuaranteed == 0, "trySetMode right after a true tryRefill was never refused");
    ok (holdTimeouts == 0, "every forced hold engaged or the clip had ended");
    ok (landedWhileHeld == 0, "no switch landed while a block was held open");
    ok (forcedRounds > 0 && forcedRefusals == forcedRounds, "every forced switch was refused while a block was held ("
                                                            + std::to_string (forcedRefusals) + "/" + std::to_string (forcedRounds) + ")");
    ok (lostLive == 0, "live monitoring was never taken down by the clip's end (" + std::to_string (lostLive) + ")");
    std::printf ("    conv: %ld live blocks, %ld unforced refusals\n", liveBlocks.load(), refused);
}

// ---- RecStream ----------------------------------------------------------------------------------------
// The audio thread numbers its admitted blocks and pushes the number 8 times: within a take the numbers
// rise in whole groups of 8, and a take begins with a block admitted AFTER its start, so its first number
// is above every number the previous take holds. A push of the old take landing after the rewind would
// put an old (smaller) number at the head.
void rec()
{
    group ("RecStream: takes, restarts at once, harvests, against a pusher");
    felitronics::core::RecStream r;
    ok (r.tryReserve (4096), "reserve");

    std::atomic<bool> go { true };
    Hold hold;
    std::thread pusher ([&] {
        long seq = 0;
        while (go.load (rlx)) {
            auto b = r.beginBlock();
            if (! b) { std::this_thread::yield(); continue; }
            const float v = (float) ++seq;
            for (int i = 0; i < 8; ++i) b.push (v);
            hold.insideBlock();
        }
    });

    long refusedStart = 0, refusedHarvest = 0, bad = 0, notRising = 0, shortTakes = 0, timeouts = 0;
    long forcedRounds = 0, forcedRefusals = 0;
    float lastOfPrevious = 0.0f;
    std::vector<float> take;
    for (long k = 1; k <= 600; ++k) {
        if (k % 50 == 0) while (! r.tryReserve (4096)) ++refusedStart;      // resize under a live pusher
        while (! r.tryStart()) ++refusedStart;
        if (! waitFor ([&] { return r.recorded() >= 64; })) ++timeouts;
        if (k % 2 == 0) {                                                    // restart at once, no harvest
            r.stop();
            while (! r.tryStart()) ++refusedStart;
            if (! waitFor ([&] { return r.recorded() >= 64; })) ++timeouts;
        }
        if (k % 5 == 0) {                                                    // forced: harvest/start while a block is held
            ++forcedRounds;
            hold.request.store (k, rlx);
            if (! waitFor ([&] { return hold.holding.load (rlx) == k; })) ++timeouts;
            std::vector<float> probe;
            const bool h = r.tryHarvest (probe), s = r.tryStart();
            if (! h && ! s) ++forcedRefusals;
            hold.release.store (k, rlx);
        }
        r.stop();
        while (! r.tryHarvest (take)) ++refusedHarvest;
        if (take.size() < 64) ++shortTakes;
        for (std::size_t i = 0; i < take.size(); ++i) {
            if (i % 8 != 0 && take[i] != take[i - 1]) ++bad;                // whole groups of 8
            if (i % 8 == 0 && i > 0 && ! (take[i] > take[i - 1])) ++notRising;
        }
        if (! take.empty()) {
            if (! (take.front() > lastOfPrevious)) ++notRising;              // the head belongs to THIS take
            lastOfPrevious = take.back();
        }
    }
    go.store (false, rlx);
    pusher.join();

    ok (timeouts == 0, "the pusher kept recording (no wait timed out)");
    ok (bad == 0, "every block's pushes are whole and in order (" + std::to_string (bad) + ")");
    ok (notRising == 0, "no take holds a sample of an earlier block than its own start (" + std::to_string (notRising) + ")");
    ok (shortTakes == 0, "every take kept what was recorded into it");
    ok (forcedRefusals == forcedRounds && forcedRounds == 120, "harvest and start were both refused while a block was held ("
                                                               + std::to_string (forcedRefusals) + "/" + std::to_string (forcedRounds) + ")");
    std::printf ("    rec: %ld start / %ld harvest unforced refusals\n", refusedStart, refusedHarvest);
}

// ---- forced ABA on the admitting CAS --------------------------------------------------------------------
// After the audio thread's load in enter() sees the mode up, the controller runs a WHOLE tryPublish (quiesce,
// copy, arm) before the audio thread's CAS: the CAS then succeeds on a word that holds the same mode as the
// load but belongs to a NEWER clip. That block must see the new clip whole, which only the CAS's acquire
// gives it. A scheduler produces this rarely, so it is forced through the Atom seam (the hook is relaxed:
// it supplies no edge of its own).
std::atomic<long> abaWant { 0 }, abaDone { 0 };
std::atomic<bool> abaGo { true };
thread_local bool abaAudioThread = false;
struct AbaU32 {
    std::atomic<std::uint32_t> a;
    AbaU32 (std::uint32_t v) noexcept : a (v) {}
    std::uint32_t load (std::memory_order o) const noexcept {
        const auto v = a.load (o);
        if (abaAudioThread && (v & 0x7fffffffu) != 0u && (v & 0x80000000u) == 0u) {   // about to CAS
            const long r = abaWant.load (rlx) + 1;
            abaWant.store (r, rlx);
            while (abaDone.load (rlx) < r && abaGo.load (rlx)) std::this_thread::yield();
        }
        return v;
    }
    bool compare_exchange_strong (std::uint32_t& e, std::uint32_t d, std::memory_order sc, std::memory_order f) noexcept {
        return a.compare_exchange_strong (e, d, sc, f);
    }
    std::uint32_t fetch_and (std::uint32_t m, std::memory_order o) noexcept { return a.fetch_and (m, o); }
};

void forcedAba()
{
    group ("AuditionStream: the admitting CAS succeeds on a re-armed word (forced ABA)");
    felitronics::core::BasicAuditionStream<AbaU32> s;
    std::vector<float> clip (64, 8.0f);
    ok (s.tryReserve (64) && s.tryPublish (clip.data(), 8), "reserve + first clip");
    std::atomic<long> admitted { 0 }, torn { 0 }, aba { 0 };
    std::thread audio ([&] {
        abaAudioThread = true;
        while (abaGo.load (rlx)) {
            auto b = s.beginBlock();
            if (! b) continue;
            admitted.fetch_add (1, rlx);
            const float m = b.data()[0];
            if (b.len() != (int) m || b.pos != 0) torn.fetch_add (1, rlx);
            for (int i = 0; i < b.len(); ++i) if (b.data()[i] != m) torn.fetch_add (1, rlx);
            if (m > 8.0f) aba.fetch_add (1, rlx);                            // admitted on a clip newer than its load
        }
    });
    long landed = 0, timeouts = 0;
    constexpr long kRounds = 2000;
    for (long r = 1; r <= kRounds; ++r) {
        if (! waitFor ([&] { return abaWant.load (rlx) >= r; }, 10.0)) { ++timeouts; break; }
        const int n = 9 + (int) (r % 50);
        for (int i = 0; i < n; ++i) clip[(std::size_t) i] = (float) n;
        if (s.tryPublish (clip.data(), (std::size_t) n)) ++landed;           // no block is open: the reader sits before its CAS
        abaDone.store (r, rlx);
    }
    abaGo.store (false, rlx);
    audio.join();
    ok (timeouts == 0 && landed == kRounds, "every publish between the audio thread's load and its CAS landed");
    ok (aba.load() >= kRounds - 1, "the blocks were admitted on the newer clip (" + std::to_string (aba.load()) + ")");
    ok (torn.load() == 0, "…and saw it whole, from 0 (" + std::to_string (torn.load()) + ")");
}

// ---- the old reader's pointer --------------------------------------------------------------------------
void storageStays()
{
    group ("storage never moves under publishes");
    felitronics::core::AuditionStream a;
    ok (a.tryReserve (16), "reserve");
    const std::vector<float> r1 (4, 1.0f), r2 (40, 2.0f);
    const float* base = nullptr;
    ok (a.tryPublish (r1.data(), r1.size()), "publish");
    { auto b = a.beginBlock(); base = b.data(); }
    bool all = true;
    for (int k = 0; k < 100; ++k) all = a.tryPublish ((k & 1) ? r1.data() : r2.data(), (k & 1) ? r1.size() : r2.size()) && all;
    ok (all, "100 republishes with no block open all land");
    { auto b = a.beginBlock(); ok (b.data() == base, "the audio thread's pointer is the same storage after 100 publishes"); }
}

} // namespace

int main()
{
    std::setvbuf (stdout, nullptr, _IONBF, 0);   // a protocol break can crash: keep every line printed before it
    std::printf ("felitronics::core rt-streams tests (threads)\n");
    storageStays();
    audition();
    conv();
    rec();
    forcedAba();
    return felitronics::test::report();
}
