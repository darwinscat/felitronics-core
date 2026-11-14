// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// RtStreams, deterministically: every interleaving of one controller call and one audio-thread block at
// the granularity that matters, the gate's atomic operations. No threads (so it runs in the single-threaded
// wasm tier too): the gate is instantiated on an atomic that calls a hook after each operation, and the hook
// runs the OTHER side's whole call there. A dry run counts the victim's operations, then every injection
// point 0 … K is tried on a fresh stream. What this proves is the protocol's LOGIC: no admitted block sees
// a half-written clip, no store a block makes at its close survives a re-arm, a refusal is followed by a
// retry that lands, and a rewind cannot be undone. The memory orders are the threaded suite's (TSan) job.

#include <felitronics/core/RtStreams.h>

// Samples here are MARKERS (1, 2, 3 …) copied, never computed: exact equality is the property under test.
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#include "felitronics_test.h"
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using felitronics::test::group;
using felitronics::test::ok;
namespace alloc = felitronics::test::alloc;

namespace {

// ---- the seam: an atomic that runs a hook after its `at`-th operation -------------------------------
struct Hook {
    int ops = 0, at = -1;
    bool fired = false, inside = false;
    void (*fn) (void*) = nullptr;
    void* ctx = nullptr;
};
Hook g;

void disarm() { g = Hook {}; }
template <class F> void armAt (int at, F* f) {
    g = Hook {};
    g.at = at;
    g.ctx = f;
    g.fn = [] (void* p) { (*(F*) p)(); };
}

struct HookedU32 {
    std::atomic<std::uint32_t> a;
    HookedU32 (std::uint32_t v) noexcept : a (v) {}
    void after() const noexcept {
        ++g.ops;
        if (g.fn != nullptr && ! g.inside && ! g.fired && g.ops == g.at) {
            g.fired = true; g.inside = true; g.fn (g.ctx); g.inside = false;
        }
    }
    std::uint32_t load (std::memory_order o) const noexcept { const auto v = a.load (o); after(); return v; }
    bool compare_exchange_strong (std::uint32_t& e, std::uint32_t d, std::memory_order s, std::memory_order f) noexcept {
        const bool r = a.compare_exchange_strong (e, d, s, f); after(); return r;
    }
    std::uint32_t fetch_and (std::uint32_t m, std::memory_order o) noexcept { const auto v = a.fetch_and (m, o); after(); return v; }
    std::uint32_t fetch_or  (std::uint32_t m, std::memory_order o) noexcept { const auto v = a.fetch_or (m, o); after(); return v; }
};

using Gate = felitronics::core::BasicRtGate<HookedU32>;
using Aud  = felitronics::core::BasicAuditionStream<HookedU32>;
using Conv = felitronics::core::BasicConvStream<HookedU32>;
using Rec  = felitronics::core::BasicRecStream<HookedU32>;

// Counts the gate operations `victim` performs: the injection points are 0 (before it) … K (after it).
template <class V> int countOps (V&& victim) { disarm(); victim(); const int k = g.ops; disarm(); return k; }

std::vector<float> clip (float marker, int n) { return std::vector<float> ((std::size_t) n, marker); }

// A clip a block saw: one marker throughout and the marker's own length — never two clips mixed.
bool whole (const float* d, int len, float marker, int wantLen) {
    if (len != wantLen) return false;
    for (int i = 0; i < len; ++i) if (d[i] != marker) return false;
    return true;
}

std::string at (const char* what, int k, int variant = -1) {
    std::string s = what;
    s += "  [inject at op " + std::to_string (k);
    if (variant >= 0) s += ", variant " + std::to_string (variant);
    return s + "]";
}

// ---- the gate on its own -----------------------------------------------------------------------------
void gateRules()
{
    group ("RtGate: the word's rules");
    Gate gate;
    ok (gate.enter() == 0u && ! gate.busy(), "a down gate admits nothing and stays idle");
    ok (gate.quiesce() && gate.arm (3u) && gate.mode() == 3u, "quiesce → arm raises the mode");
    ok (! gate.arm (1u) && gate.mode() == 3u, "arm on a raised word is refused (a second arm never ORs modes)");
    ok (gate.enter() == 3u && gate.busy(), "a block is admitted at the mode and holds BUSY");
    ok (gate.enter() == 0u, "a second block of the same gate is not admitted while one is open");
    ok (! gate.quiesce() && gate.mode() == 0u && gate.busy(), "quiesce with a block open: refused, mode down, BUSY kept");
    ok (! gate.arm (1u) && gate.mode() == 0u, "arm while a block is open is refused (the word is not exactly 0)");
    gate.leave (false);
    ok (! gate.busy() && gate.mode() == 0u, "leave(false) clears BUSY and keeps the (down) mode");
    ok (gate.quiesce(), "the retry after the block closed lands");
    ok (gate.arm (2u) && gate.enter() == 2u, "re-armed and admitted");
    gate.leave (true);
    ok (gate.mode() == 0u && ! gate.busy(), "leave(true) takes the mode down with BUSY");
    ok (! gate.arm (0u) && gate.mode() == 0u, "arm(0) raises nothing");
    ok (gate.quiesce() && gate.arm (0xffffffffu) && gate.mode() == Gate::kModeMask && ! gate.busy(),
        "a mode is masked to 31 bits: arming never fakes BUSY");
}

// ---- the streams' API, one thread ----------------------------------------------------------------------
void auditionApi()
{
    group ("AuditionStream: API");
    Aud a;
    ok (a.tryReserve (8) && a.capacity() == 8 && ! a.isPlaying(), "tryReserve sizes the fixed storage, stopped");
    const auto r1 = clip (1.0f, 3);
    ok (a.tryPublish (r1.data(), r1.size()) && a.isPlaying(), "a publish arms the stream");
    const float* base = nullptr;
    { auto b = a.beginBlock(); ok ((bool) b && b.pos == 0 && whole (b.data(), b.len(), 1.0f, 3), "a block sees the render from 0"); base = b.data(); }
    const auto big = clip (2.0f, 16);
    ok (a.tryPublish (big.data(), big.size()), "an oversized render is accepted");
    { auto b = a.beginBlock(); ok ((bool) b && whole (b.data(), b.len(), 2.0f, 8), "…and truncated to the storage");
      ok (b.data() == base, "a publish never reallocates"); }
    ok (a.tryPublish (big.data(), 0) && ! a.isPlaying(), "a zero-length publish lands and leaves the stream stopped");
    { auto b = a.beginBlock(); ok (! b, "…so no block is admitted (the loop-path OOB pin)"); }
    ok (a.tryPublish (nullptr, 5) && ! a.isPlaying(), "a null render is an empty one");

    ok (a.tryPublish (r1.data(), r1.size()), "republish");
    { auto b = a.beginBlock(); b.pos = -7; }
    { auto b = a.beginBlock(); ok (b.pos == 0, "a negative position is clamped to 0 on close"); b.pos = 99; }
    { auto b = a.beginBlock(); ok (b.pos == 3, "a position past the clip is clamped to its length on close"); b.finish(); }
    ok (! a.isPlaying(), "finish() takes the stream down on close");

    ok (a.tryPublish (r1.data(), r1.size()), "republish");
    {
        auto b = a.beginBlock();
        ok (! a.tryReserve (16) && a.capacity() == 8, "tryReserve is refused while a block is open, storage untouched");
        ok (! a.isPlaying(), "…and the refusal stopped the stream");
    }
    a.stop();
    ok (a.tryReserve (4) && a.capacity() == 4, "tryReserve lands once the block closed");

    const auto r2 = clip (3.0f, 4);
    ok (a.tryPublish (r2.data(), r2.size()), "publish");
    const long long c0 = alloc::count.load();
    for (int k = 0; k < 64; ++k) {
        { auto b = a.beginBlock(); if (b) b.pos += 1; }
        (void) a.tryPublish (r2.data(), r2.size());
        a.stop();
    }
    ok (alloc::count.load() == c0, "blocks, publishes and stops never allocate");
}

void convApi()
{
    group ("ConvStream: API");
    Conv c;
    ok (c.tryReserve (4) && c.capacity() == 4 && c.mode() == 0, "tryReserve sizes the fixed storage, off");
    const auto di = clip (7.0f, 5);
    ok (c.tryRefill (di.data(), di.size()) && c.mode() == 0, "tryRefill lands and leaves the route off");
    ok (c.trySetMode (Conv::kClip) && c.mode() == Conv::kClip, "trySetMode right after a true tryRefill lands");
    { auto b = c.beginBlock(); ok ((bool) b && b.mode() == Conv::kClip && whole (b.data(), b.len(), 7.0f, 4),
                                   "a clip block sees the truncated clip, length = the clip not the storage"); b.pos = 2; }
    ok (c.trySetMode (Conv::kLive) && c.mode() == Conv::kLive, "a mode switch lands with no block open");
    { auto b = c.beginBlock(); ok ((bool) b && b.mode() == Conv::kLive && b.pos == 2, "a mode switch keeps the clip position"); }
    {
        auto b = c.beginBlock();
        ok (! c.trySetMode (Conv::kClip), "a mode switch is refused while a block is open");
        ok (c.trySetMode (0) && c.mode() == 0, "mode 0 stops and is never refused, even with a block open");
        ok (c.trySetMode (-3), "a negative mode stops too");
        ok (! c.tryRefill (di.data(), di.size()), "tryRefill is refused while a block is open");
    }
    ok (c.tryRefill (di.data(), 2) && c.trySetMode (Conv::kClip), "refill + switch land after it closed");
    { auto b = c.beginBlock(); ok (b.pos == 0 && b.len() == 2, "a refill rewinds"); b.finish(); }
    ok (c.mode() == 0, "finish() ends the route on close");

    const long long c0 = alloc::count.load();
    for (int k = 0; k < 64; ++k) {
        (void) c.tryRefill (di.data(), di.size());
        (void) c.trySetMode (1 + (k & 1));
        { auto b = c.beginBlock(); if (b) b.pos += 1; }
        c.stop();
    }
    ok (alloc::count.load() == c0, "blocks, refills, mode switches and stops never allocate");
}

void recApi()
{
    group ("RecStream: API");
    Rec r;
    ok (r.tryReserve (6) && r.capacity() == 6 && ! r.isRecording() && r.recorded() == 0, "tryReserve sizes, stopped, empty");
    ok (r.tryStart() && r.isRecording(), "tryStart records");
    { auto b = r.beginBlock(); for (float v : { 1.f, 2.f, 3.f }) b.push (v); }
    ok (r.recorded() == 3, "recorded() reads the closed block's length");
    { auto b = r.beginBlock(); for (float v : { 4.f, 5.f, 6.f, 7.f, 8.f }) b.push (v); ok (b.full(), "the block sees the take full"); }
    ok (r.recorded() == 6 && r.full(), "capacity caps the take, no overrun");
    std::vector<float> take;
    ok (r.tryHarvest (take) && ! r.isRecording(), "tryHarvest lands and stops");
    ok (take == std::vector<float> ({ 1, 2, 3, 4, 5, 6 }), "the take is every pushed sample in order");
    ok (r.tryStart() && r.recorded() == 0, "a restart rewinds");
    {
        auto b = r.beginBlock();
        b.push (9.f);
        std::vector<float> t2 { 42.f };
        ok (! r.tryHarvest (t2) && t2.size() == 1 && t2[0] == 42.f, "tryHarvest is refused while a block is open, take untouched");
        ok (! r.tryStart(), "tryStart is refused while a block is open");
        ok (! r.tryReserve (12), "tryReserve is refused while a block is open");
        b.push (9.f);
    }
    ok (r.tryHarvest (take) && take == std::vector<float> ({ 9, 9 }), "the block's pushes after the stop are in the take");
    ok (r.tryReserve (3) && r.recorded() == 0 && r.tryHarvest (take) && take.empty(), "tryReserve discards the take");

    ok (r.tryStart(), "start");
    const long long c0 = alloc::count.load();
    for (int k = 0; k < 64; ++k) { auto b = r.beginBlock(); if (b) b.push ((float) k); }
    r.stop();
    ok (r.tryStart(), "restart");
    ok (alloc::count.load() == c0, "blocks, pushes, stop and start never allocate");
}

// ---- interleavings: AuditionStream -------------------------------------------------------------------
// Variant: 0 the block plays a little · 1 the block plays its clip to the end and finish()es · 2 the block
// opens and does nothing. The controller's call is a publish of a NEW clip (marker 2, length 3) over an OLD
// one (marker 1, length 4), or a stop.
void auditionInterleave()
{
    group ("AuditionStream: a publish at every point of a block");
    const auto oldc = clip (1.0f, 4), newc = clip (2.0f, 3);
    for (int variant = 0; variant < 3; ++variant) {
        struct Seen { bool admitted = false, sawOld = false, sawNew = false, intact = true; int startPos = -1; };
        auto block = [&] (Aud& a, Seen& s) {
            auto b = a.beginBlock();
            if (! b) return;
            s.admitted = true; s.startPos = b.pos;
            s.sawOld = b.len() > 0 && b.data()[0] == 1.0f; s.sawNew = b.len() > 0 && b.data()[0] == 2.0f;
            s.intact = s.sawOld ? whole (b.data(), b.len(), 1.0f, 4) : whole (b.data(), b.len(), 2.0f, 3);
            if (variant == 0) b.pos = 2;
            if (variant == 1) { b.pos = b.len(); b.finish(); }
        };
        int K = 0;
        { Aud a; (void) a.tryReserve (8); (void) a.tryPublish (oldc.data(), oldc.size()); Seen s; K = countOps ([&] { block (a, s); }); }
        for (int k = 0; k <= K; ++k) {
            Aud a; (void) a.tryReserve (8); (void) a.tryPublish (oldc.data(), oldc.size());
            Seen s; bool landed = false;
            auto inject = [&] { landed = a.tryPublish (newc.data(), newc.size()); };
            if (k == 0) inject(); else armAt (k, &inject);
            block (a, s);
            const bool fired = k == 0 || g.fired;
            disarm();
            ok (fired, at ("the injection point was reached", k, variant));
            ok (s.admitted, at ("the armed clip was admitted", k, variant));
            ok (s.intact, at ("an admitted block reads ONE whole clip", k, variant));
            ok (s.startPos == 0, at ("…from its start", k, variant));
            if (! landed) {
                ok (s.sawOld, at ("a refusal means the block held the OLD clip", k, variant));
                ok (a.tryPublish (newc.data(), newc.size()), at ("the refused publish lands on the retry", k, variant));
            }
            if (s.sawOld) {      // every store that block made at its close belongs to the OLD clip
                ok (a.isPlaying(), at ("the re-arm survives the old block's close", k, variant));
                auto b = a.beginBlock();
                ok ((bool) b && b.pos == 0 && whole (b.data(), b.len(), 2.0f, 3), at ("the new clip starts at 0, whole", k, variant));
            } else if (s.sawNew) { // admitted after the publish: its stores are the new clip's own
                if (variant == 1) ok (! a.isPlaying(), at ("the new clip's own end takes it down", k, variant));
                else {
                    auto b = a.beginBlock();
                    ok ((bool) b && b.pos == (variant == 0 ? 2 : 0) && b.len() == 3, at ("the new clip keeps its own position", k, variant));
                }
            }
        }
    }

    group ("AuditionStream: a stop at every point of a block");
    for (int variant = 0; variant < 3; ++variant) {
        auto block = [&] (Aud& a) {
            auto b = a.beginBlock();
            if (! b) return;
            if (variant == 0) b.pos = 2;
            if (variant == 1) { b.pos = b.len(); b.finish(); }
        };
        int K = 0;
        { Aud a; (void) a.tryReserve (8); (void) a.tryPublish (oldc.data(), oldc.size()); K = countOps ([&] { block (a); }); }
        for (int k = 0; k <= K; ++k) {
            Aud a; (void) a.tryReserve (8); (void) a.tryPublish (oldc.data(), oldc.size());
            auto inject = [&] { a.stop(); };
            if (k == 0) inject(); else armAt (k, &inject);
            block (a);
            disarm();
            ok (! a.isPlaying(), at ("a stop is never undone by a block's close", k, variant));
            auto b = a.beginBlock();
            ok (! b, at ("…and no block is admitted after it", k, variant));
        }
    }

    group ("AuditionStream: a block at every point of a publish");
    {
        Aud a0; (void) a0.tryReserve (8); (void) a0.tryPublish (oldc.data(), oldc.size());
        const int K = countOps ([&] { (void) a0.tryPublish (newc.data(), newc.size()); });
        for (int k = 0; k <= K; ++k) {
            Aud a; (void) a.tryReserve (8); (void) a.tryPublish (oldc.data(), oldc.size());
            bool admitted = false, good = true, sawNew = false;
            auto inject = [&] {
                auto b = a.beginBlock();
                if (! b) return;
                admitted = true;
                sawNew = b.len() > 0 && b.data()[0] == 2.0f;
                good = b.pos == 0 && (sawNew ? whole (b.data(), b.len(), 2.0f, 3) : whole (b.data(), b.len(), 1.0f, 4));
                b.pos = 1;
            };
            if (k == 0) inject(); else armAt (k, &inject);
            const bool landed = a.tryPublish (newc.data(), newc.size());
            disarm();
            ok (landed, at ("a publish with no block open lands", k));
            ok (good, at ("a block admitted around a publish sees one whole clip from 0", k));
            if (admitted && k > 0) ok (sawNew, at ("a block admitted INSIDE a publish sees the NEW clip", k));
            auto b = a.beginBlock();
            ok ((bool) b && b.len() == 3 && b.data()[0] == 2.0f, at ("the new clip plays after the publish", k));
            ok (b.pos == ((admitted && sawNew) ? 1 : 0), at ("…from 0, or from where its own block left it", k));
        }
    }
}

// ---- interleavings: ConvStream -----------------------------------------------------------------------
void convInterleave()
{
    group ("ConvStream: a mode switch at every point of a block that ends the clip");
    const auto oldc = clip (1.0f, 4), newc = clip (2.0f, 3);
    auto block = [] (Conv& c) {
        auto b = c.beginBlock();
        if (! b) return;
        if (b.mode() == Conv::kClip) { b.pos = b.len(); b.finish(); }       // the clip ended in this block
    };
    {
        int K = 0;
        { Conv c; (void) c.tryReserve (8); (void) c.tryRefill (oldc.data(), oldc.size()); (void) c.trySetMode (Conv::kClip);
          K = countOps ([&] { block (c); }); }
        for (int k = 0; k <= K; ++k) {
            Conv c; (void) c.tryReserve (8); (void) c.tryRefill (oldc.data(), oldc.size()); (void) c.trySetMode (Conv::kClip);
            bool landed = false;
            auto inject = [&] { landed = c.trySetMode (Conv::kLive); };
            if (k == 0) inject(); else armAt (k, &inject);
            block (c);
            disarm();
            if (! landed) ok (c.trySetMode (Conv::kLive), at ("a refused switch lands on the retry", k));
            ok (c.mode() == Conv::kLive, at ("live monitoring survives the clip's end (the old block cannot take it down)", k));
        }
    }

    group ("ConvStream: refill + switch at every point of a clip block");
    {
        auto play = [] (Conv& c, bool& sawOld, bool& good) {
            auto b = c.beginBlock();
            if (! b) return;
            sawOld = b.len() == 4;
            good = b.pos == 0 && (sawOld ? whole (b.data(), b.len(), 1.0f, 4) : whole (b.data(), b.len(), 2.0f, 3));
            b.pos = 1;
        };
        int K = 0;
        { Conv c; (void) c.tryReserve (8); (void) c.tryRefill (oldc.data(), oldc.size()); (void) c.trySetMode (Conv::kClip);
          bool so = false, gd = true; K = countOps ([&] { play (c, so, gd); }); }
        for (int k = 0; k <= K; ++k) {
            Conv c; (void) c.tryReserve (8); (void) c.tryRefill (oldc.data(), oldc.size()); (void) c.trySetMode (Conv::kClip);
            bool refilled = false, switched = false, sawOld = false, good = true;
            auto inject = [&] { refilled = c.tryRefill (newc.data(), newc.size()); if (refilled) switched = c.trySetMode (Conv::kClip); };
            if (k == 0) inject(); else armAt (k, &inject);
            play (c, sawOld, good);
            disarm();
            ok (good, at ("a clip block reads one whole clip from 0", k));
            if (refilled) ok (switched, at ("trySetMode right after a true tryRefill is never refused", k));
            else {
                ok (sawOld, at ("a refused refill means the block held the old clip", k));
                ok (c.tryRefill (newc.data(), newc.size()) && c.trySetMode (Conv::kClip), at ("the retry lands", k));
            }
            auto b = c.beginBlock();
            ok ((bool) b && b.mode() == Conv::kClip && b.len() == 3 && b.data()[0] == 2.0f, at ("the new clip plays", k));
            ok (b.pos == (sawOld ? 0 : 1), at ("…from 0, or from where its own block left it", k));
        }
    }

    group ("ConvStream: a block at every point of a refill");
    {
        Conv c0; (void) c0.tryReserve (8); (void) c0.tryRefill (oldc.data(), oldc.size()); (void) c0.trySetMode (Conv::kClip);
        const int K = countOps ([&] { (void) c0.tryRefill (newc.data(), newc.size()); });
        for (int k = 0; k <= K; ++k) {
            Conv c; (void) c.tryReserve (8); (void) c.tryRefill (oldc.data(), oldc.size()); (void) c.trySetMode (Conv::kClip);
            bool admitted = false, good = true;
            auto inject = [&] {
                auto b = c.beginBlock();
                if (! b) return;
                admitted = true;
                good = whole (b.data(), b.len(), 1.0f, 4);       // a refill never re-arms: only the OLD clip can be seen
            };
            if (k == 0) inject(); else armAt (k, &inject);
            const bool landed = c.tryRefill (newc.data(), newc.size());
            disarm();
            ok (landed, at ("a refill with no block open lands", k));
            ok (good, at ("a block around a refill sees only the old clip, whole", k));
            if (k > 0) ok (! admitted, at ("no block is admitted inside a refill", k));
            ok (c.mode() == 0, at ("the refill left the route off", k));
        }
    }
}

// ---- interleavings: RecStream ------------------------------------------------------------------------
void recInterleave()
{
    // The old take is marker 1; a push block pushes three samples of marker `m` (the take it was admitted
    // to). Injected: a harvest, a start (the new take, marker 2), a stop, a reserve.
    group ("RecStream: a controller call at every point of a push block");
    for (int what = 0; what < 4; ++what) {
        // `m` by reference: the audio thread pushes what it hears when it pushes, so a block admitted after an
        // injected start pushes the new take's marker.
        auto block = [] (Rec& r, const float& m) { auto b = r.beginBlock(); if (! b) return; b.push (m); b.push (m); b.push (m); };
        auto fresh = [] (Rec& r) { (void) r.tryReserve (16); (void) r.tryStart(); { auto b = r.beginBlock(); b.push (1.f); } };
        int K = 0;
        { Rec r; fresh (r); K = countOps ([&] { block (r, 1.f); }); }
        for (int k = 0; k <= K; ++k) {
            Rec r; fresh (r);
            bool landed = false; float marker = 1.f;
            std::vector<float> take;
            auto inject = [&] {
                switch (what) {
                    case 0: landed = r.tryHarvest (take); break;
                    case 1: landed = r.tryStart(); if (landed) marker = 2.f; break;
                    case 2: r.stop(); landed = true; break;
                    default: landed = r.tryReserve (16); break;
                }
            };
            if (k == 0) inject(); else armAt (k, &inject);
            block (r, marker);     // the marker is read when the block runs: a block admitted after a start pushes 2
            disarm();
            std::vector<float> fin;
            ok (r.tryHarvest (fin), at ("a harvest with no block open lands", k, what));
            for (float v : fin) if (v != fin.front()) { ok (false, at ("a take holds ONE take's samples", k, what)); break; }
            switch (what) {
                case 0:
                    if (landed) ok (take == std::vector<float> ({ 1.f }) || take == std::vector<float> ({ 1, 1, 1, 1 }),
                                    at ("a harvest that lands holds a whole take, before or after the block", k, what));
                    ok (fin.size() == take.size() || ! landed, at ("…and its length is final (no block appends after it)", k, what));
                    break;
                case 1:
                    if (landed) ok (fin.empty() || (fin.size() == 3 && fin[0] == 2.f),
                                    at ("a start that lands rewinds: the new take is empty or the block's own pushes", k, what));
                    else ok (fin.size() == 4 && fin[0] == 1.f, at ("a refused start leaves the old take whole", k, what));
                    break;
                case 2:
                    ok (! r.isRecording(), at ("a stop is never undone by a block's close", k, what));
                    ok (fin.size() == 1 || fin.size() == 4, at ("a block admitted before the stop finishes its pushes", k, what));
                    break;
                default:
                    if (landed) ok (fin.empty(), at ("a reserve that lands discards the take, no late push survives", k, what));
                    else ok (fin.size() == 4, at ("a refused reserve leaves the take", k, what));
                    break;
            }
        }
    }

    group ("RecStream: a push block at every point of a start");
    {
        Rec r0; (void) r0.tryReserve (16);
        const int K = countOps ([&] { (void) r0.tryStart(); });
        for (int k = 0; k <= K; ++k) {
            Rec r; (void) r.tryReserve (16); (void) r.tryStart(); { auto b = r.beginBlock(); b.push (1.f); } r.stop();
            bool admitted = false;
            auto inject = [&] { auto b = r.beginBlock(); if (! b) return; admitted = true; b.push (2.f); };
            if (k == 0) inject(); else armAt (k, &inject);
            const bool landed = r.tryStart();
            disarm();
            ok (landed, at ("a start with no block open lands", k));
            ok (k == 0 || admitted == (k >= K), at ("a block is admitted inside a start only after its arm", k));
            std::vector<float> fin;
            ok (r.tryHarvest (fin), at ("harvest", k));
            if (admitted && k > 0) ok (fin == std::vector<float> ({ 2.f }), at ("a block admitted after the arm pushes from 0", k));
            else ok (fin.empty(), at ("the rewind holds", k));
        }
    }
}

} // namespace

int main()
{
    std::printf ("felitronics::core rt-streams interleaving tests\n");
    gateRules();
    auditionApi();
    convApi();
    recApi();
    auditionInterleave();
    convInterleave();
    recInterleave();
    return felitronics::test::report();
}
