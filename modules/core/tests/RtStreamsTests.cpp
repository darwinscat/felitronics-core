// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The RT buffer-swap discipline types (promoted from OrbitCapture). The one property that matters
// and that no UI test can see: a publish/refill NEVER moves the buffer's memory once it's been
// reserved — the RT reader holds the pointer across swaps. Plus a discipline-respecting concurrent
// smoke (publisher vs reader / RT pusher vs message-thread stop) that gives ASan/UBSan something
// real to chew on.

#include <felitronics/core/RtStreams.h>

#include "felitronics_test.h"

#include <thread>
#include <vector>

using felitronics::test::group;
using felitronics::test::ok;

int main()
{
    std::printf ("felitronics::core rt-streams tests\n");

    group ("AuditionStream::publish");
    {
        felitronics::core::AuditionStream a;
        a.reserve (8);
        ok (a.buf.size() == 8, "reserve sizes the FIXED storage up front");
        const float* base = a.buf.data();
        const std::vector<float> r1 { 1, 2, 3 };
        a.publish (r1.data(), r1.size());
        ok (a.playing.load() && a.len.load() == 3 && a.pos.load() == 0, "publish arms the stream");
        ok (a.buf[0] == 1.0f && a.buf[2] == 3.0f, "buffer carries the render");

        const std::vector<float> big (16, 0.5f);                       // over the fixed size: truncate, never grow
        a.publish (big.data(), big.size());
        ok (a.len.load() == 8 && a.buf.size() == 8, "oversized render truncates to the fixed storage");
        ok (a.buf.data() == base, "publish never reallocates — the RT pointer stays valid");

        // Crew pin: a zero-length render must NOT arm — the reader's loop path resets pos to 0
        // and dereferences buf[0], which on an armed empty clip is an out-of-bounds read.
        a.publish (big.data(), 0);
        ok (! a.playing.load() && a.len.load() == 0 && a.pos.load() == 0,
            "zero-length publish rewinds but never arms (loop-path OOB pin)");

        a.stop();
        ok (! a.playing.load(), "stop takes the flag down");
    }

    group ("ConvStream::refill");
    {
        felitronics::core::ConvStream c;
        c.reserve (4);
        const float* base = c.buf.data();
        c.pos.store (99);
        const std::vector<float> di { 7, 8, 9, 10, 11 };
        c.refill (di.data(), di.size());                               // mode is down (0) — allowed
        ok (c.buf.size() == 4 && c.buf[0] == 7.0f, "refill truncates to capacity");
        ok (c.pos.load() == 0, "refill rewinds the read position");
        ok (c.buf.data() == base, "refill never reallocates");
        ok (c.mode.load() == 0, "refill itself never raises the mode — the caller arms it");

        c.setMode (1);
        ok (c.mode.load() == 1, "setMode arms the conv route");
        c.stop();
        ok (c.mode.load() == 0, "stop = mode 0");
    }

    group ("RecStream: the live-input recorder");
    {
        felitronics::core::RecStream r;
        r.reserve (6);
        const float* base = r.buf.data();
        r.start();
        ok (r.recording.load() && r.len.load() == 0, "start arms and rewinds");
        for (float v : { 1.f, 2.f, 3.f }) r.push (v);
        ok (r.len.load() == 3 && r.buf[0] == 1.f && r.buf[2] == 3.f, "pushes append in order");
        for (float v : { 4.f, 5.f, 6.f, 7.f, 8.f }) r.push (v);
        ok (r.len.load() == 6 && r.full(), "capacity caps the recording (no overrun)");
        ok (r.buf.data() == base, "recording never reallocates — the RT pointer stays valid");
        r.stop();
        ok (! r.recording.load() && r.len.load() == 6, "stop keeps the recorded length");
        r.start();
        ok (r.len.load() == 0, "restart rewinds");
    }

    group ("concurrent smoke (discipline respected; ASan/UBSan food)");
    {
        // Audition: a reader streaming buf[pos..len) whenever `playing` is up, while the writer
        // re-publishes renders. The reader mimics the RT callback: check the flag, then read.
        felitronics::core::AuditionStream a;
        a.reserve (1 << 12);
        std::atomic<bool> go { true };
        std::atomic<long> consumed { 0 };
        float sink = 0.0f;
        std::thread reader ([&] {
            while (go.load (std::memory_order_acquire)) {
                if (a.playing.load (std::memory_order_acquire)) {
                    const int n = a.len.load();
                    for (int i = 0; i < n; ++i) sink += a.buf[(std::size_t) i];
                    consumed.fetch_add (n);
                }
            }
        });
        std::vector<float> render (1 << 12, 0.25f);
        for (int k = 0; k < 2000; ++k)
            a.publish (render.data(), render.size());       // stop -> assign-in-capacity -> start
        go.store (false, std::memory_order_release);
        reader.join();
        ok (consumed.load() >= 0 && sink >= 0.0f, "publisher/reader ran to completion under the discipline");

        // Recorder: an "RT" pusher against message-thread stop/harvest cycles. The harvest reads
        // CONTENT, not just length — with a relaxed len publish, len=i+1 may become visible
        // before buf[i] and this check can catch a zero where a sample must be (crew finding;
        // push() now release-stores len).
        felitronics::core::RecStream r;
        r.reserve (1 << 12);
        std::atomic<bool> go2 { true };
        std::thread pusher ([&] {
            while (go2.load (std::memory_order_acquire))
                if (r.recording.load (std::memory_order_acquire)) r.push (1.0f);
        });
        long harvested = 0; bool contentOk = true;
        for (int k = 0; k < 500; ++k) {
            r.start();
            while (r.len.load() < 64 && ! r.full()) {}      // let the pusher work
            r.stop();
            const int n = r.len.load();                     // seq_cst load = acquire: pairs with push's release
            for (int i = 0; i < n; ++i) contentOk = contentOk && r.buf[(std::size_t) i] == 1.0f;
            harvested += n;                                 // message thread reads buf[0..len) after stop
        }
        go2.store (false, std::memory_order_release);
        pusher.join();
        ok (harvested >= 500L * 64L, "recorder harvested every cycle without overrun");
        ok (contentOk, "every harvested sample is the pushed value (len release-publishes the sample)");
    }

    return felitronics::test::report();
}
