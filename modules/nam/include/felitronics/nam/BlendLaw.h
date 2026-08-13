// SPDX-License-Identifier: MIT
#pragma once
// felitronics::nam — WHO IS AUDIBLE, HOW LOUDLY, AND WHEN A MODEL MAY BE REPLACED.
//
// Two model slots play the same input and are mixed linearly; a knob between two captures asks for a
// fraction of each. That much is arithmetic. What is not arithmetic — and what cost a day of chasing
// crackle by ear — is WHO decides the number. This header exists because that answer must be "one
// thing, and it is testable without a sound card".
//
// THE FAILURE THIS PREVENTS. In the player this replaces, three separate places wrote the mix or
// swapped a model: the blend itself, the code that follows the nearest capture, and a warm-up gate
// bolted on later. None was wrong alone. Together they replaced models under a live gain and drove
// the weight across its whole range inside single 10 ms blocks. Measured on one sweep of a nine-
// capture device: 62 model loads, 401 parked retries, and a full 1.000 swing of the weight.
//
// SO: one writer, one recurrence, and every property below provable by construction rather than by
// listening.
//
//   1. The applied weight moves by at most `maxDeltaPerBlock` per block. There is no other writer,
//      so no path — a load, a device change, a give-up, a mode flip — can step it.
//   2. A model is replaced only in a slot whose applied weight is exactly zero.
//   3. A model that has not yet been fed its receptive field is never audible. Its weight is zero,
//      not small: a network that has not heard the last 132 ms does not produce a quiet version of
//      the right sound, it produces the wrong sound.
//   4. Progress is guaranteed without timeouts. Whenever a slot holds the wrong model, some slot's
//      goal is zero and the recurrence reaches it in a bounded number of blocks.
//
// WHAT GIVES, since something must: not (1) and not (3) — they constrain a value and a derivative
// and are compatible. What gives is INSTANTANEOUS ACCURACY. Crossing a capture, the incoming model
// arrives about 132 ms late and the outgoing capture carries alone until it does. That is a morph
// that trails the hand, which is what a person hears as a knob rather than as a fault.

#include <algorithm>
#include <cstdint>

namespace felitronics::nam {

using BlendModelId = std::uint64_t;      // 0 = this slot holds nothing
inline constexpr int kBlendSlots = 2;

struct BlendPolicy {
    // How far the weight may travel in one block. The full crossfade therefore takes at least
    // 1/maxDeltaPerBlock blocks — at 512 samples and 48 kHz, 0.25 gives 42.7 ms. This is the one
    // number a listener gets to argue with; everything else is derived from it.
    double maxDeltaPerBlock = 0.25;
};

struct BlendRequest {
    BlendModelId want[kBlendSlots] {};   // which capture belongs in each slot (assigned by the caller)
    double       targetB = 0.0;          // the weight slot 1 should end up with; slot 0 gets 1 - it
};

// What the law asks the message thread to do. At most one is outstanding at a time: a second swap
// while one is in flight would need a slot that is neither audible nor settled.
struct BlendLoad {
    bool         wanted = false;
    int          slot = 0;
    BlendModelId model = 0;
};

struct BlendState {
    BlendModelId held[kBlendSlots] {};
    long long    fed[kBlendSlots] {};    // samples of real signal since this model landed
    long long    need[kBlendSlots] {};   // …and how many it must have before it may be heard
    bool         inFlight[kBlendSlots] {};
    double       x = 0.0;                // THE applied weight of slot 1. The only state that sounds.
};

struct BlendStep {
    double    beginB = 0.0, endB = 0.0;  // the weight at the block's first and last sample
    BlendLoad load;
};

// A slot is SAFE TO HEAR when it holds a model that has been fed enough. Note this does not ask
// whether it holds the WANTED model: an old capture is still a real one, and letting it carry while
// its replacement warms is the whole reason a crossing sounds like a knob instead of a hole.
inline bool blendAudible(const BlendState& s, int i) {
    return s.held[i] != 0 && ! s.inFlight[i] && s.fed[i] >= s.need[i];
}

// One audio block. Advances the fed counters, decides the goal, moves the weight toward it by at
// most one step, and — only where a slot has arrived at exactly zero — asks for the swap.
inline BlendStep blendStep(BlendState& s, const BlendRequest& r, int blockSamples,
                           const BlendPolicy& p = {}) {
    BlendStep out;
    out.beginB = s.x;
    if (blockSamples <= 0) { out.endB = s.x; return out; }

    for (int i = 0; i < kBlendSlots; ++i)
        if (s.held[i] != 0 && ! s.inFlight[i]) s.fed[i] += blockSamples;

    // The goal, in priority order. Anything unsafe to hear outranks anything merely wrong, and being
    // wrong outranks the request — because a slot cannot take its new model until it is silent.
    const bool unsafe0 = ! blendAudible(s, 0), unsafe1 = ! blendAudible(s, 1);
    const bool wrong0 = s.held[0] != r.want[0], wrong1 = s.held[1] != r.want[1];
    double goal;
    // Nothing may be heard yet. One of them sounds anyway — 132 ms of a model converging on the
    // truth reads as a knob, while 132 ms of silence reads as a fault. But NEVER an empty slot: a
    // stage with no model passes its input straight through, and the raw DI is both wrong and about
    // ten decibels louder than any capture. An empty slot is silent or it is a disaster.
    if (unsafe0 && unsafe1) goal = (s.held[0] != 0 || s.held[1] == 0) ? 0.0 : 1.0;
    else if (unsafe1)       goal = 0.0;
    else if (unsafe0)       goal = 1.0;
    else if (wrong0 && wrong1) goal = s.x <= 0.5 ? 0.0 : 1.0;   // evacuate the LIGHTER one first and
                                                                // freeze the heavier: two swaps in a
                                                                // row, never a hole between them
    else if (wrong0)        goal = 1.0;
    else if (wrong1)        goal = 0.0;
    else                    goal = std::clamp(r.targetB, 0.0, 1.0);

    const double d = std::clamp(p.maxDeltaPerBlock, 1.0e-9, 1.0);
    s.x = std::clamp(s.x + std::clamp(goal - s.x, -d, d), 0.0, 1.0);
    out.endB = s.x;

    // A swap is asked for only at EXACTLY zero, and only one at a time. Everything above conspires to
    // make that reachable: a wrong slot's goal is a rail, and the rail is its own zero.
    const bool busy = s.inFlight[0] || s.inFlight[1];
    if (! busy)
        for (int i = 0; i < kBlendSlots; ++i) {
            const double w = i == 0 ? 1.0 - s.x : s.x;
            if (s.held[i] != r.want[i] && w == 0.0 && r.want[i] != 0) {
                out.load = { true, i, r.want[i] };
                s.inFlight[i] = true;
                break;
            }
        }
    return out;
}

// The message thread reporting back. `prewarm` is that model's own receptive field in HOST samples —
// see NamStage::prewarmSamples(), and convert it, because the model counts in its own rate.
inline void blendLanded(BlendState& s, int slot, BlendModelId model, long long prewarm) {
    if (slot < 0 || slot >= kBlendSlots) return;
    s.held[slot] = model;
    s.need[slot] = std::max(0LL, prewarm);
    s.fed[slot] = 0;
    s.inFlight[slot] = false;
}

// A load that could not be honoured (unreadable file, wrong rate). The slot keeps whatever it had,
// which is a real model, and the law will try again on the next block.
inline void blendLoadFailed(BlendState& s, int slot) {
    if (slot >= 0 && slot < kBlendSlots) s.inFlight[slot] = false;
}

} // namespace felitronics::nam
