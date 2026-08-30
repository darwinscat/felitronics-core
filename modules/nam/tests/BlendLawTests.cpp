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
    BlendPolicy policy;
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
        const auto st = blendStep(s, r, kBlock, policy);

        worstStep = std::max(worstStep, std::abs(st.endB - st.beginB));
        worstSeam = std::max(worstSeam, std::abs(st.beginB - prevEnd));
        prevEnd = st.endB;

        // NOTHING UNFED MAY BE HEARD — anywhere inside the block, and counted against the samples
        // that had been fed BEFORE it, since a counter crossing the line mid-block does not make the
        // first half of that block fed. The overall gain counts too: muted is not heard.
        for (int i = 0; i < kBlendSlots; ++i) {
            const double w = std::max(i == 0 ? 1.0 - st.beginB : st.beginB,
                                      i == 0 ? 1.0 - st.endB   : st.endB)
                           * std::max(st.beginGain, st.endGain);
            if (s.held[i] != 0 && (s.cold[i] || s.fed[i] - kBlock < s.need[i]) && w > 1.0e-9) unfedHeard = true;
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
        // A block where NEITHER slot is audible is a hole — tolerable only at a cold start.
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

        BlendState both;                       // and with neither holding anything at all
        BlendRequest none;
        BlendStep last;
        for (int b = 0; b < 40; ++b) last = blendStep(both, none, kBlock);
        ok(both.x == 0.0, "…and with nothing anywhere it does not thrash");
        ok(last.endGain == 0.0, "…and stays MUTED rather than letting a passthrough out");

        // A freshly landed model is silent for its whole receptive field, gain and all.
        BlendState cold; BlendRequest want; want.want[0] = 5; want.targetB = 0.0;
        blendLanded(cold, 0, 5, kPrewarm);
        bool heard = false;
        int blocks = 0;
        for (; blocks < 40; ++blocks) {
            const auto st = blendStep(cold, want, kBlock);
            const double w = std::max(st.beginGain, st.endGain) * std::max(1.0 - st.beginB, 1.0 - st.endB);
            if (cold.fed[0] - kBlock < kPrewarm && w > 1e-9) heard = true;
            if (st.endGain > 0.5) break;
        }
        ok(! heard, "a cold model is inaudible for its whole receptive field");
        ok(blocks >= kPrewarm / kBlock, "…which is at least " + std::to_string(kPrewarm / kBlock) + " blocks");
    }

    group("a request repeated forever changes nothing");
    {
        auto s = settledAt(150.0);
        Run run; run.s = s; run.prevEnd = s.x;
        const auto same = requestAt(150.0);
        for (int b = 0; b < 500; ++b) run.block(same, 1);
        ok(run.loads == 0 && run.worstStep <= 1e-12, "no load, no movement - idle is idle");
    }

    group("a slot at rest goes cold: at exactly zero, unchanged, and only once the rest is over");
    {
        // THE ECONOMY. A dial parked on a capture keeps the neighbour in the other slot at exactly
        // zero, and that neighbour's network used to run every block for nothing — a second WaveNet
        // pass for silence. The law may call the slot cold, and the host stop running it, only once
        // it has been at rest for as long as the host allows: a pause between two moves of the hand
        // is not a rest.
        BlendPolicy p; p.coldAfterSamples = 2 * 48000;          // two seconds, in the law's own unit
        const auto r = requestAt(150.0);                        // on a capture: by parity slot 1 carries
        const int rest = 0, carry = 1;                          // …and slot 0 waits at exactly zero
        auto s = settledAt(150.0);
        ok(s.x >= 1.0 && s.held[rest] == r.want[rest], "the waiting slot holds the wanted neighbour at zero");
        int blocks = 0;
        BlendStep last;
        while (! s.cold[rest] && blocks < 400) { last = blendStep(s, r, kBlock, p); ++blocks; }
        ok(s.cold[rest], "the waiting slot goes cold");
        ok((long long) (blocks - 1) * kBlock >= p.coldAfterSamples,
           "…only once the whole rest has been PLAYED (" + std::to_string(blocks) + " blocks)");
        ok((long long) blocks * kBlock < p.coldAfterSamples + 3 * kBlock, "…and not long after it");
        ok(! s.cold[carry], "the carrying slot does not: it has weight");
        ok(! blendAudible(s, rest) && blendAudible(s, carry), "cold reads as unfed — inaudible; the carrier still sounds");
        ok(s.x >= 1.0 && s.held[rest] == r.want[rest] && ! last.load.wanted, "nothing moved and nothing was replaced");
        for (int b = 0; b < 200; ++b) last = blendStep(s, r, kBlock, p);
        ok(s.cold[rest] && s.x >= 1.0 && ! last.load.wanted, "…and under the same request it sleeps on: no load, no move");
        ok(s.fed[rest] >= s.need[rest] && s.need[rest] == kPrewarm, "its counters are left as they were: the field is remembered for the wake");

        // WAKING. The hand moves to 160, wanting a third of the cold slot. It wakes on the first block
        // of the new request — re-landed, not reloaded — and is fed its whole field before a sample of
        // it is heard; then the weight travels as it always does. The first turn after a rest trails
        // the hand by one warm-up, the price already paid at every crossing.
        Run run; run.s = s; run.policy = p; run.prevEnd = s.x;
        const auto r2 = requestAt(160.0);
        ok(r2.want[rest] == r.want[rest] && r2.targetB < 1.0, "the new request wants the sleeping model, in its slot, with weight");
        run.block(r2, 1);
        ok(! run.s.cold[rest], "the first block of the new request wakes it");
        ok(run.s.fed[rest] == kBlock && run.s.need[rest] == kPrewarm, "…fed from zero, with the same field to serve");
        ok(run.s.held[rest] == r2.want[rest] && run.loads == 0, "…the same model, in place: no load was asked");
        int warm = 1;
        while (run.s.x >= 1.0 && warm < 100) { run.block(r2, 1); ++warm; }
        ok(run.s.x < 1.0, "the weight moves once the slot is fed");
        ok(warm > kPrewarm / kBlock, "…and not before: " + std::to_string(warm) + " blocks, the field is "
                                     + std::to_string(kPrewarm / kBlock));
        ok(warm <= kPrewarm / kBlock + 3, "…nor much after");
        for (int b = 0; b < 40; ++b) run.block(r2, 1);
        ok(run.loads == 0, "no load in the whole wake: the model never left");
        ok(! run.unfedHeard, "nothing unfed — and nothing cold — was heard on the way");
        ok(run.worstStep <= kMaxDelta + 1e-9, "…and the weight never stepped");
        ok(std::abs(run.s.x - r2.targetB) < 1e-6, "it arrives where it was sent");
        ok(! run.s.cold[0] && ! run.s.cold[1], "…and with both slots heard, neither can sleep");
    }

    group("the rest is counted as played: the block that completes it runs, the next sleeps");
    {
        // A rest of ONE block. The first block under the new request starts the count; the second is
        // the rest, and runs; the third finds the rest complete and is the first the host may skip.
        BlendPolicy p; p.coldAfterSamples = kBlock;
        auto s = settledAt(150.0);
        const auto r = requestAt(150.0);
        blendStep(s, r, kBlock, p);
        ok(! s.cold[0], "the first block of a request is never a rest");
        blendStep(s, r, kBlock, p);
        ok(! s.cold[0] && s.still[0] == kBlock, "the block that completes the rest is still played");
        blendStep(s, r, kBlock, p);
        ok(s.cold[0], "…and the one after it is the first asleep");
    }

    group("a slot asked for nothing is at rest too");
    {
        // The caller wants one capture and nothing beside it: the other slot keeps its old model (the
        // law never swaps for `want == 0`) at exactly zero — that is a rest, and it may sleep.
        BlendPolicy p; p.coldAfterSamples = 2 * 48000;
        auto s = settledAt(150.0);
        auto r = requestAt(150.0); r.want[0] = 0;
        BlendStep last;
        for (int b = 0; b < 400; ++b) last = blendStep(s, r, kBlock, p);
        ok(s.cold[0] && s.held[0] != 0 && ! last.load.wanted, "the slot wanted empty sleeps with its old model in place");
        ok(! s.cold[1] && s.x >= 1.0, "…and the one carrying does not");
    }

    group("a slot with any weight never goes cold; the law's own default never sleeps");
    {
        BlendPolicy p; p.coldAfterSamples = 2 * 48000;
        auto both = settledAt(165.0);                            // halfway between two captures
        const auto r = requestAt(165.0);
        for (int b = 0; b < 400; ++b) blendStep(both, r, kBlock, p);
        ok(! both.cold[0] && ! both.cold[1], "four seconds between two captures: both stay awake");

        auto s = settledAt(150.0);
        const auto on = requestAt(150.0);
        for (int b = 0; b < 400; ++b) blendStep(s, on, kBlock);     // the default policy
        ok(! s.cold[0] && ! s.cold[1], "four seconds on a capture under the default policy: nothing sleeps");
        BlendPolicy never; never.coldAfterSamples = 0;
        for (int b = 0; b < 400; ++b) blendStep(s, on, kBlock, never);
        ok(! s.cold[0] && ! s.cold[1], "…and zero means never, spelled out");

        // A hand that keeps moving, however little, keeps the request changing — and a changed request
        // is never a rest.
        auto wobble = settledAt(150.0); int coldBlocks = 0;
        for (int b = 0; b < 400; ++b) {
            blendStep(wobble, requestAt(150.0 + 0.01 * (double) (b % 2)), kBlock, p);
            if (wobble.cold[0] || wobble.cold[1]) ++coldBlocks;
        }
        ok(coldBlocks == 0, "a hand that wobbles a hundredth of a degree never lets a slot sleep");
    }

    group("a load that failed is refused: asked once, then not until the request changes");
    {
        // THE STORM THIS PREVENTS: a file that cannot be a model fails identically every time, and the
        // law used to ask for it again on every block — a fetch and a parse per timer tick, for ever.
        auto s = settledAt(150.0);
        auto r = requestAt(120.0);                            // wants model 3 in slot 0, all of it
        int asks = 0;
        BlendStep st = blendStep(s, r, kBlock);
        if (st.load.wanted) ++asks;
        for (int b = 0; b < 10 && asks == 0; ++b) { st = blendStep(s, r, kBlock); if (st.load.wanted) ++asks; }
        ok(asks == 1 && st.load.slot == 0 && st.load.model == 3, "the load is asked for");
        blendLoadFailed(s, 0);
        for (int b = 0; b < 400; ++b) { st = blendStep(s, r, kBlock); if (st.load.wanted) ++asks; }
        ok(asks == 1, "…and after the failure, never again under the same request");
        ok(s.held[0] == 5 && s.x >= 1.0, "the slot keeps its old capture, and the pair partner carries");
        BlendPolicy p; p.coldAfterSamples = 2 * 48000;
        for (int b = 0; b < 400; ++b) blendStep(s, r, kBlock, p);
        ok(s.cold[0], "a refused slot is at rest: it may sleep");

        auto back = requestAt(150.0);                          // slot 0's wish changes: 3 → 5, already held
        for (int b = 0; b < 10; ++b) { st = blendStep(s, back, kBlock, p); if (st.load.wanted) ++asks; }
        ok(asks == 1 && s.refused[0] == 0, "a new wish clears the refusal — here it is already held, so no load");
        for (int b = 0; b < 10 && asks == 1; ++b) { st = blendStep(s, r, kBlock); if (st.load.wanted) ++asks; }
        ok(asks == 2 && st.load.model == 3, "…and asking for the failed model AGAIN is a new question, tried anew");
    }

    group("a cold slot is replaced the way any silent slot is — and lands awake");
    {
        // The hand jumps to 120, wanting a different model in the sleeping slot. The change of request
        // wakes the slot, the law sees it wrong and silent and asks for the swap at once, and the
        // landing is a wake in its own right: fed from zero, then heard.
        BlendPolicy p; p.coldAfterSamples = 2 * 48000;
        auto s = settledAt(150.0);
        for (int b = 0; b < 400; ++b) blendStep(s, requestAt(150.0), kBlock, p);
        ok(s.cold[0], "asleep at 150");
        Run run; run.s = s; run.policy = p; run.prevEnd = s.x;
        const auto far = requestAt(120.0);
        ok(far.want[0] != s.held[0] && far.want[1] == s.held[1] && far.targetB <= 0.0,
           "120 wants another model in the sleeping slot, and all of it");
        for (int b = 0; b < 200; ++b) run.block(far, 3);
        ok(run.loads == 1, "exactly one load: the new model, into the slot that slept");
        ok(! run.swapAboveZero && ! run.unfedHeard, "…swapped silent, heard only fed");
        ok(run.gapBlocks == 0, "…and the old carrier sounded throughout");
        ok(run.s.held[0] == far.want[0] && run.s.x <= 0.0, "it arrives: the new model carries alone");
        for (int b = 0; b < 200; ++b) run.block(far, 3);
        ok(run.s.cold[1] && ! run.s.cold[0], "…and two seconds later the OTHER slot, now the silent one, sleeps");
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
        blendDelay(x, tail, (int) cap, 2, 4);
        ok(x[0] == 0.0f && x[1] == 0.0f && x[2] == 1.0f && x[3] == 2.0f,
           "two samples of delay pushes the block back by two");
        float y[4] = { 5, 6, 7, 8 };
        blendDelay(y, tail, (int) cap, 2, 4);
        ok(y[0] == 3.0f && y[1] == 4.0f && y[2] == 5.0f && y[3] == 6.0f,
           "…and the next block continues from the tail, with nothing lost");

        float z[4] = { 1, 2, 3, 4 };
        float clean[cap] {};
        blendDelay(z, clean, (int) cap, 0, 4);
        ok(z[0] == 1.0f && z[3] == 4.0f, "no delay changes nothing at all");
    }

    return felitronics::test::report();
}
