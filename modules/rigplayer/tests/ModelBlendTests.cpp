// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Headless tests for rigplayer/ModelBlend.h — which two captures sound at a knob position. The rule is
// small on purpose; what is checked here is that it is small in the RIGHT places, because every one
// of these cases is a way for a dial to go dead or to jump.

#include <felitronics_test.h>
#include <felitronics/rigplayer/ModelBlend.h>

using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;
using namespace felitronics::rigplayer;

int main() {
    std::printf("felitronics::rigplayer::ModelBlend tests\n");

    // A nine-point dial as this library shoots them: 60 to 300 degrees, every 30.
    const std::vector<double> nine { 60, 90, 120, 150, 180, 210, 240, 270, 300 };

    group("on a capture, that capture is the whole sound");
    {
        // The PAIR may still name the neighbour above — deliberately. It is kept loaded at weight
        // zero so that leaving the knot does not have to start a model from an empty memory, which
        // is the one thing that makes a moving knob click. What matters here is only the weight.
        for (int i = 0; i < (int) nine.size(); ++i) {
            const auto p = pickBlend(nine, nine[(std::size_t) i]);
            ok(p.lo == i && p.mixHi == 0.0,
               "at " + std::to_string((int) nine[(std::size_t) i]) + " degrees it plays at full weight");
        }
    }

    group("between two captures the mix follows the ANGLE");
    {
        const auto mid = pickBlend(nine, 105.0);
        ok(mid.lo == 1 && mid.hi == 2, "the two neighbours are the ones either side");
        approx(mid.mixHi, 0.5, 1e-9, "…and halfway across is half of each");
        approx(pickBlend(nine, 100.0).mixHi, 1.0 / 3.0, 1e-9, "a third of the way is a third");
        approx(pickBlend(nine, 119.0).mixHi, 29.0 / 30.0, 1e-9, "…and just short of the next one is nearly all of it");
    }

    group("the mix never jumps as the knob crosses a capture");
    {
        // THE property the whole scheme rests on: approach a knot from below and the upper neighbour's
        // weight goes to 1; step past it and the pair re-forms with that same model at weight 0. If
        // these two did not meet, every capture would be a click.
        const auto below = pickBlend(nine, 149.999), above = pickBlend(nine, 150.001);
        ok(nine[(std::size_t) below.hi] == 150.0, "below, the upper neighbour is the capture at 150");
        ok(nine[(std::size_t) above.lo] == 150.0, "above, the LOWER neighbour is that same capture");
        approx(below.mixHi, 1.0, 1e-3, "…arriving at full weight");
        approx(above.mixHi, 0.0, 1e-3, "…and leaving from none");
    }

    group("below the lowest capture, the lowest one plays SOFTER");
    {
        // The dial turns its whole travel whatever was shot on it. Under the bottom capture there is
        // nothing to mix with, so that one plays alone with less going into it — the safe direction,
        // since a model fed less stays inside what it was trained on.
        const auto p = pickBlend(nine, 0.0);
        ok(p.lo == 0 && p.hi == 0, "no second model");
        approx(p.extendDb, -60.0 * kTopExtendDbPerDeg, 1e-9, "…and it is fed softer by the angle below it");
        ok(pickBlend(nine, 60.0).extendDb == 0.0, "exactly ON the bottom capture nothing is taken away");
        ok(pickBlend(nine, 30.0).extendDb > pickBlend(nine, 0.0).extendDb, "…and it deepens as the dial falls");
    }

    group("above the highest capture, that one is fed harder");
    {
        // The exception the bench asked for: a dial that turns to 17h when the last capture is at 15h
        // must not have a dead top third.
        const auto p = pickBlend(nine, 360.0);
        ok(p.lo == 8 && p.hi == 8, "the top capture, alone");
        approx(p.extendDb, 60.0 * kTopExtendDbPerDeg, 1e-9, "…driven by the angle past it");
        ok(pickBlend(nine, 300.0).extendDb == 0.0, "and exactly ON the top capture nothing is added");
    }

    group("a device with one capture still plays");
    {
        const std::vector<double> one { 150 };
        ok(pickBlend(one, 60.0).lo == 0 && pickBlend(one, 150.0).lo == 0, "everywhere below and at it");
        ok(pickBlend(one, 60.0).extendDb < 0.0, "…softer below");
        ok(pickBlend(one, 300.0).extendDb > 0.0, "…and above it the one model is driven");
    }

    group("no captures at all: nothing to play, and it says so");
    {
        ok(pickBlend({}, 150.0).lo < 0, "-1 rather than a slot that does not exist");
    }

    group("stepped() is the same choice with the crossfade taken away");
    {
        ok(stepped(pickBlend(nine, 100.0)).lo == 1, "a third of the way up still plays the lower one");
        ok(stepped(pickBlend(nine, 110.0)).lo == 2, "…and past halfway it plays the upper one");
        const auto s = stepped(pickBlend(nine, 110.0));
        ok(s.lo == s.hi && s.mixHi == 0.0, "one model, no mix - which is what it is for");
        approx(stepped(pickBlend(nine, 360.0)).extendDb, 60.0 * kTopExtendDbPerDeg, 1e-9,
               "…and the top extension survives stepping: it is not a crossfade");
    }

    group("an unevenly spaced dial is read by angle, not by index");
    {
        // The measured-control grids in this library are 7h/8h/9h/12h/15h/17h - deliberately uneven.
        // Halfway between two captures means halfway in ROTATION, which is what a hand feels.
        const std::vector<double> uneven { 0, 60, 120, 300 };
        approx(pickBlend(uneven, 210.0).mixHi, 0.5, 1e-9, "the middle of a wide gap is half of each");
        approx(pickBlend(uneven, 150.0).mixHi, 30.0 / 180.0, 1e-9,
               "…and a sixth into it is a sixth, not the half an index would have given");
    }

    return felitronics::test::report();
}
