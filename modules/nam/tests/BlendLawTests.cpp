// SPDX-License-Identifier: MIT
// The mixing law, driven the way a hand drives it — and checked without a sound card.
//
// This file exists because the same defect was chased seven times by ear in one day and every fix
// was judged by "does it still crackle". Each property below is one of those failures turned into a
// number. If any of them regresses, this goes red before anybody plays a note.

#include <felitronics/nam/BlendLaw.h>
#include <felitronics_test.h>

#include <cmath>
#include <vector>

using felitronics::test::ok;
using felitronics::test::group;
using namespace felitronics::nam;

namespace {

constexpr double kSr = 48000.0;
constexpr int    kBlock = 512;
constexpr long long kPrewarm = 6332;          // the receptive field of the captures here
constexpr double kMaxDelta = 0.25;

// Nine captures at thirty degrees apart, one model each, in the slots parity gives them.
constexpr double kFirst = 60.0, kStep = 30.0;
constexpr int    kKnots = 9;

int knotBelow(double deg) {
    const int i = (int) std::floor((deg - kFirst) / kStep);
    return std::clamp(i, 0, kKnots - 2);
}

// What the CALLER asks for at an angle: the pair either side, each in the slot its index's parity
// names, and the weight of slot 1.
BlendRequest requestAt(double deg) {
    const int lo = knotBelow(deg);
    const double t = std::clamp((deg - (kFirst + kStep * lo)) / kStep, 0.0, 1.0);
    BlendRequest r;
    const bool loEven = (lo % 2) == 0;
    r.want[loEven ? 0 : 1] = (BlendModelId) (lo + 1);          // ids are 1-based; 0 means "nothing"
    r.want[loEven ? 1 : 0] = (BlendModelId) (lo + 2);
    r.targetB = loEven ? t : 1.0 - t;
    return r;
}

// A run of the law with a scripted hand. `loadBlocks` is how many blocks a load takes to land.
struct Run {
    BlendState  s;
    int         loads = 0, gapBlocks = 0;
    double      worstStep = 0.0;
    double      worstSeam = 0.0;
    bool        unfedHeard = false, swapAboveZero = false, twoInFlight = false;
    double      prevEnd = 0.0;
    // A pending load, and how many blocks until it lands.
    bool pending = false; int pendSlot = 0; BlendModelId pendModel = 0; int pendLeft = 0;

    void block(const BlendRequest& r, int loadBlocks) {
        // THE MODEL ACTUALLY CHANGES HERE. That is the instant the weight must be zero — not the
        // instant the swap was asked for, which is a block or more earlier.
        if (pending && --pendLeft <= 0) {
            const double w = pendSlot == 0 ? 1.0 - s.x : s.x;
            if (w != 0.0) swapAboveZero = true;
            blendLanded(s, pendSlot, pendModel, kPrewarm);
            pending = false;
        }
        const bool wasInFlight = s.inFlight[0] || s.inFlight[1];
        const auto st = blendStep(s, r, kBlock);

        worstStep = std::max(worstStep, std::abs(st.endB - st.beginB));
        worstSeam = std::max(worstSeam, std::abs(st.beginB - prevEnd));
        prevEnd = st.endB;

        // NOTHING UNFED MAY BE HEARD. Checked on the weight the block ENDS with as well as begins,
        // because a slot may rise inside a block.
        for (int i = 0; i < kBlendSlots; ++i) {
            const double w = std::max(i == 0 ? 1.0 - st.beginB : st.beginB,
                                      i == 0 ? 1.0 - st.endB   : st.endB);
            if (s.held[i] != 0 && s.fed[i] < s.need[i] && w > 1.0e-9) unfedHeard = true;
        }
        // …and nothing may be REPLACED unless its weight is exactly zero.
        if (st.load.wanted) {
            // Asked for only once the slot has ARRIVED at zero — checked on the weight the block
            // ends with, since that is the weight the request was made against.
            const double w = st.load.slot == 0 ? 1.0 - st.endB : st.endB;
            if (w != 0.0) swapAboveZero = true;
            if (wasInFlight) twoInFlight = true;
            ++loads;
            pending = true; pendSlot = st.load.slot; pendModel = st.load.model; pendLeft = std::max(1, loadBlocks);
        }
        // A block where NEITHER slot is audible is a hole — the thing the law exists to avoid.
        if (! blendAudible(s, 0) && ! blendAudible(s, 1)) ++gapBlocks;
    }
};

// Turn from `from` to `to` over `seconds`, publishing the hand's position every block.
Run sweep(double from, double to, double seconds, int loadBlocks = 1, BlendState seed = {}) {
    Run run; run.s = seed; run.prevEnd = seed.x;
    const int blocks = std::max(1, (int) std::round(seconds * kSr / kBlock));
    for (int b = 0; b <= blocks; ++b)
        run.block(requestAt(from + (to - from) * (double) b / (double) blocks), loadBlocks);
    return run;
}

// Both slots already holding the pair at `deg`, fed and settled — a player that has been running.
BlendState settledAt(double deg) {
    BlendState s;
    const auto r = requestAt(deg);
    for (int i = 0; i < kBlendSlots; ++i) { s.held[i] = r.want[i]; s.need[i] = kPrewarm; s.fed[i] = kPrewarm; }
    s.x = r.targetB;
    return s;
}

} // namespace

int main() {
    std::printf("felitronics::nam blend-law tests\n");

    group("the weight never steps, at any speed");
    {
        // FAILURE 5 and 6 in one line: a warm gate that slammed the target to a rail produced 0.513
        // and then 1.000 swings inside a single block. Nothing here can, because one recurrence owns
        // the number and it is clamped.
        for (const double secs : { 8.0, 2.0, 0.5, 0.15 }) {
            const auto r = sweep(60.0, 300.0, secs, 1, settledAt(60.0));
            ok(r.worstStep <= kMaxDelta + 1e-9,
               "a " + std::to_string((int) (secs * 1000)) + " ms sweep moves at most one step per block");
            ok(r.worstSeam <= 1e-12, "…and consecutive blocks meet exactly, with no seam");
        }
    }

    group("a model is replaced only in a slot that is silent");
    {
        // FAILURE 4: loading whichever slot happened to be quiet, while the other still carried.
        // FAILURE 7 / the real one: a second writer swapping slot A at the midpoint, where its weight
        // is about a half.
        for (const double secs : { 8.0, 2.0, 0.5, 0.15 }) {
            const auto r = sweep(60.0, 300.0, secs, 3, settledAt(60.0));
            ok(! r.swapAboveZero, "no swap at " + std::to_string((int) (secs * 1000)) + " ms happens above zero");
            ok(! r.twoInFlight, "…and never two loads at once");
        }
    }

    group("a model nobody has fed is never heard");
    {
        // FAILURE 7: with the gate taken out entirely, a freshly loaded network became audible while
        // still describing the silence it was born into.
        for (const double secs : { 8.0, 0.5, 0.15 }) {
            const auto r = sweep(60.0, 300.0, secs, 3, settledAt(60.0));
            ok(! r.unfedHeard, "nothing unfed is audible on a " + std::to_string((int) (secs * 1000)) + " ms sweep");
        }
    }

    group("the dial is never silent");
    {
        for (const double secs : { 8.0, 0.5, 0.15 }) {
            const auto r = sweep(60.0, 300.0, secs, 3, settledAt(60.0));
            ok(r.gapBlocks == 0, "some real capture is audible in every block of a "
                                 + std::to_string((int) (secs * 1000)) + " ms sweep");
        }
    }

    group("the work is proportional to the captures crossed, not to the events");
    {
        // FAILURE: 62 loads and 401 parked retries on a nine-capture device, because a second path
        // reloaded a slot at every midpoint and the parking counted calls rather than waits.
        const auto slow = sweep(60.0, 300.0, 8.0, 1, settledAt(60.0));
        ok(slow.loads <= kKnots, "a full sweep loads at most nine models, not sixty-two");
        const auto fast = sweep(60.0, 300.0, 0.15, 1, settledAt(60.0));
        ok(fast.loads <= kKnots, "…and a fast sweep does not load MORE, because the last request wins");
    }

    group("a hand that does not cross a capture loads nothing");
    {
        // Dithering around the midpoint between two captures used to reload a slot on every wobble:
        // the NEAREST capture flips there, and something was following nearest.
        auto s = settledAt(105.0);
        Run run; run.s = s; run.prevEnd = s.x;
        for (int b = 0; b < 300; ++b)
            run.block(requestAt(105.0 + 2.0 * std::sin(0.1 * b)), 1);
        ok(run.loads == 0, "three seconds of wobbling across a midpoint loads nothing at all");
        ok(run.worstStep <= kMaxDelta + 1e-9, "…and the weight still never steps");
    }

    group("a jerk across the whole dial replaces both, one at a time");
    {
        // FAILURE 4 again, from the other side: the pair at the far end shares no model with the pair
        // being played, so BOTH slots must change — and they cannot both be silent at once.
        auto s = settledAt(75.0);
        Run run; run.s = s; run.prevEnd = s.x;
        const auto far = requestAt(285.0);
        for (int b = 0; b < 200; ++b) run.block(far, 3);
        ok(! run.swapAboveZero && ! run.unfedHeard, "both replacements happen silently and fed");
        ok(run.gapBlocks == 0, "…without a single silent block in between");
        ok(run.s.held[0] == far.want[0] && run.s.held[1] == far.want[1], "and it arrives where it was sent");
        ok(std::abs(run.s.x - far.targetB) < 1e-6, "…at the weight that was asked for");
        ok(run.loads == 2, "exactly two loads: no timeout, no retry storm");
    }

    group("an empty slot is never the one that sounds");
    {
        // A NamStage with no model is a PASSTHROUGH, not silence. Letting an empty slot carry would
        // put the raw DI on the output — wrong, and some ten decibels louder than any capture.
        BlendState s;                          // nothing held anywhere
        s.held[1] = 7; s.need[1] = kPrewarm; s.fed[1] = kPrewarm;   // …except a fed model in slot 1
        BlendRequest r; r.want[0] = 0; r.want[1] = 7; r.targetB = 0.5;
        for (int b = 0; b < 40; ++b) blendStep(s, r, kBlock);
        ok(std::abs(s.x - 1.0) < 1e-9, "the weight goes to the slot that actually holds a model");

        BlendState both;                       // and with neither holding anything, slot 0 by default
        BlendRequest none;
        for (int b = 0; b < 40; ++b) blendStep(both, none, kBlock);
        ok(both.x == 0.0, "…and with nothing anywhere it does not thrash");
    }

    group("a request repeated forever changes nothing");
    {
        auto s = settledAt(150.0);
        Run run; run.s = s; run.prevEnd = s.x;
        const auto same = requestAt(150.0);
        for (int b = 0; b < 500; ++b) run.block(same, 1);
        ok(run.loads == 0 && run.worstStep <= 1e-12, "no load, no movement - idle is idle");
    }

    group("the mix is the two slots, weighted — and nothing else");
    {
        // THE SLIP THIS CATCHES: one of the two model calls went missing in the host, so the mix was
        // the raw DI against a model rather than model against model. It was audible only as the
        // loudness rippling with the knob, loud at every even capture and quiet at every odd one,
        // because a DI is some ten decibels above a normalised capture. Here it is arithmetic.
        constexpr int n = 8;
        float a[n], b[n];
        for (int i = 0; i < n; ++i) { a[i] = 1.0f; b[i] = 3.0f; }
        blendMix(a, b, n, 0.0, 0.0);
        bool allA = true; for (const float v : a) allA = allA && std::abs(v - 1.0f) < 1e-6f;
        ok(allA, "at weight zero the first slot is the whole sound");

        for (int i = 0; i < n; ++i) { a[i] = 1.0f; b[i] = 3.0f; }
        blendMix(a, b, n, 1.0, 1.0);
        bool allB = true; for (const float v : a) allB = allB && std::abs(v - 3.0f) < 1e-6f;
        ok(allB, "…and at one, the second");

        for (int i = 0; i < n; ++i) { a[i] = 1.0f; b[i] = 3.0f; }
        blendMix(a, b, n, 0.5, 0.5);
        bool half = true; for (const float v : a) half = half && std::abs(v - 2.0f) < 1e-6f;
        ok(half, "halfway is the average - linear, so two coherent halves stay at unity");

        // Ramping across the block: the last sample must have arrived at the end weight.
        for (int i = 0; i < n; ++i) { a[i] = 0.0f; b[i] = 1.0f; }
        blendMix(a, b, n, 0.0, 1.0);
        ok(std::abs(a[n - 1] - 1.0f) < 1e-6f, "a ramp lands exactly on the weight the block ends with");
        bool rising = true;
        for (int i = 1; i < n; ++i) rising = rising && a[i] > a[i - 1];
        ok(rising, "…and gets there monotonically, without a step");
    }

    group("a slot's delay line carries between blocks");
    {
        // The correction for two captures sitting a few samples apart. It has to survive a block
        // boundary, and it has to keep advancing while the slot is silent — otherwise the first
        // audible samples read a line full of whatever was there before.
        constexpr std::size_t cap = 128;
        float tail[cap] {};
        float x[4] = { 1, 2, 3, 4 };
        blendDelay(x, tail, 2, 4);
        ok(x[0] == 0.0f && x[1] == 0.0f && x[2] == 1.0f && x[3] == 2.0f,
           "two samples of delay pushes the block back by two");
        float y[4] = { 5, 6, 7, 8 };
        blendDelay(y, tail, 2, 4);
        ok(y[0] == 3.0f && y[1] == 4.0f && y[2] == 5.0f && y[3] == 6.0f,
           "…and the next block continues from the tail, with nothing lost");

        float z[4] = { 1, 2, 3, 4 };
        float clean[cap] {};
        blendDelay(z, clean, 0, 4);
        ok(z[0] == 1.0f && z[3] == 4.0f, "no delay changes nothing at all");
    }

    return felitronics::test::report();
}
