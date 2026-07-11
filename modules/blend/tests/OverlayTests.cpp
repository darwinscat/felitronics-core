// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The mix-view overlay facade (Overlay.h, reuse audit N6). Two tiers of truth:
//   1. a NULL against the hand-composed pipeline (processedMic/blendIrs/logMagnitudeCurve/
//      interferenceDb) — the facade adds NO numerics of its own, machine-epsilon identical;
//   2. analytical pins of the baked-in display semantics — single audible mic → zero
//      interference; two in-phase copies → +3 dB constructive; a 180° flip → destructive;
//      the Master-|H| fade only ever SHRINKS the tint and never flips its sign.

#include <felitronics/blend/Overlay.h>

#include "felitronics_test.h"

#include <cmath>
#include <limits>

using felitronics::test::group;
using felitronics::test::ok;
namespace fb  = felitronics::blend;
namespace off = felitronics::analysis::offline;

// A deterministic "IR": a decaying multi-tone burst, distinct per seed.
static std::vector<float> makeIr (int n, double sr, int seed)
{
    std::vector<float> v ((std::size_t) n);
    for (int i = 0; i < n; ++i) {
        const double t = i / sr;
        v[(std::size_t) i] = (float) (std::exp (-t * 60.0)
            * (std::sin (2.0 * 3.14159265358979 * (200.0 + 130.0 * seed) * t)
               + 0.5 * std::sin (2.0 * 3.14159265358979 * (1700.0 + 90.0 * seed) * t)));
    }
    return v;
}

int main()
{
    std::printf ("felitronics::blend overlay tests\n");
    const double sr = 48000.0;
    const std::vector<std::vector<float>> irs { makeIr (4096, sr, 0), makeIr (4096, sr, 1) };

    group ("guards");
    {
        fb::MasterParams master;
        std::vector<fb::StripParams> one (1);
        ok (fb::makeOverlay ({}, {}, master, sr).blend.empty(), "no irs -> empty");
        ok (fb::makeOverlay (irs, one, master, sr).blend.empty(), "strip-count mismatch -> empty");
        ok (fb::makeOverlay (irs, std::vector<fb::StripParams> (2), master, 0.0).blend.empty(), "sr <= 0 -> empty");
        const std::vector<std::vector<float>> holey { {}, irs[1] };    // crew pin: no PARTIAL overlays
        const auto h = fb::makeOverlay (holey, std::vector<fb::StripParams> (2), master, sr);
        ok (h.blend.empty() && h.mic.empty() && h.interference.empty(), "any empty ir -> fully empty (never partial)");

        std::vector<fb::StripParams> sick (2);                         // crew pin: heal non-finite params
        sick[0].gainDb = std::nan ("");
        sick[1].phaseDeg = std::numeric_limits<double>::infinity();
        const auto healed = fb::makeOverlay (irs, sick, master, sr);
        const auto plain  = fb::makeOverlay (irs, std::vector<fb::StripParams> (2), master, sr);
        bool same = healed.blend.size() == plain.blend.size() && ! healed.blend.empty();
        for (std::size_t p = 0; same && p < healed.blend.size(); ++p)
            same = healed.blend[p] == plain.blend[p] && std::isfinite (healed.blend[p]);
        ok (same, "non-finite gain/phase heal to defaults (offline convention) — never a NaN picture");
    }

    group ("NULL vs the hand-composed pipeline");
    {
        std::vector<fb::StripParams> strips (2);
        strips[0].gainDb = -3.0;
        strips[1].phaseDeg = 45.0;
        strips[1].hpf = { true, 120.0, 24 };
        fb::MasterParams master;
        master.gainDb = -2.0;
        master.lpf = { true, 9000.0, 12 };

        const auto o = fb::makeOverlay (irs, strips, master, sr);

        off::LogCurveSpec cs; cs.normalize = false;
        bool micNull = true;
        for (std::size_t m = 0; m < irs.size(); ++m) {                 // per-mic: processedMic * gain -> curve
            auto gi = fb::processedMic (irs[m], strips[m], sr);
            const double g = std::pow (10.0, strips[m].gainDb / 20.0);
            for (auto& v : gi) v *= (float) g;                         // the app's exact expression (bit-parity)
            const auto c = off::logMagnitudeCurve (std::span<const float> (gi), sr, cs);
            micNull = micNull && c.size() == o.mic[m].size();
            for (std::size_t p = 0; micNull && p < c.size(); ++p) micNull = c[p] == o.mic[m][p];
        }
        ok (micNull, "per-mic curves NULL the primitives exactly");

        const auto post = fb::blendIrs (irs, strips, master, sr, true);
        const auto bc = off::logMagnitudeCurve (std::span<const float> (post), sr, cs);
        bool blendNull = bc.size() == o.blend.size();
        for (std::size_t p = 0; blendNull && p < bc.size(); ++p) blendNull = bc[p] == o.blend[p];
        ok (blendNull, "blend curve NULLs blendIrs(applyMaster=true) exactly");

        const auto pre = fb::blendIrs (irs, strips, master, sr, false);
        const auto rc = off::logMagnitudeCurve (std::span<const float> (pre), sr, cs);
        std::vector<off::MicCurveView> mv;
        for (std::size_t m = 0; m < o.mic.size(); ++m)
            mv.push_back ({ std::span<const double> (o.mic[m]), o.audible[m] });
        auto interf = off::interferenceDb (rc, mv);
        std::vector<float> imp (irs[0].size(), 0.0f); imp[0] = 1.0f;   // the |H| fade, recomputed by hand
        fb::applyFilters (imp, master.hpf, master.lpf, sr);
        const auto mc = off::logMagnitudeCurve (std::span<const float> (imp), sr, cs);
        for (std::size_t p = 0; p < interf.size() && p < mc.size(); ++p)
            interf[p] *= std::pow (10.0, std::min (0.0, mc[p]) / 20.0);
        bool intNull = interf.size() == o.interference.size();
        for (std::size_t p = 0; intNull && p < interf.size(); ++p) intNull = interf[p] == o.interference[p];
        ok (intNull, "interference NULLs the pre-Master basis + |H| fade exactly");
    }

    group ("analytical pins of the display semantics");
    {
        fb::MasterParams master;                                       // no master filters
        std::vector<fb::StripParams> strips (2);

        strips[1].mute = true;                                         // single audible mic -> no interference
        auto o = fb::makeOverlay (irs, strips, master, sr);
        ok (o.audible[0] && ! o.audible[1], "mute flags flow into audible[]");
        double worst = 0.0;
        for (double v : o.interference) worst = std::max (worst, std::abs (v));
        ok (worst < 1e-6, "one audible mic -> interference is identically ~0 dB");

        const std::vector<std::vector<float>> twin { irs[0], irs[0] }; // identical, in phase
        std::vector<fb::StripParams> both (2);
        o = fb::makeOverlay (twin, both, master, sr);
        double mid = 0.0; int cnt = 0;
        for (std::size_t p = o.interference.size() / 4; p < 3 * o.interference.size() / 4; ++p) { mid += o.interference[p]; ++cnt; }
        mid /= std::max (1, cnt);
        ok (std::abs (mid - 3.0103) < 0.05, "two identical in-phase mics -> +3.01 dB constructive everywhere");

        both[1].phaseDeg = 180.0;                                      // flip one -> destructive
        o = fb::makeOverlay (twin, both, master, sr);
        double deep = 0.0;
        for (double v : o.interference) deep = std::min (deep, v);
        ok (deep < -20.0, "a 180-degree flip of an identical copy digs a deep destructive notch");
    }

    group ("Master |H| fade: shrinks the tint, never flips or grows it");
    {
        std::vector<fb::StripParams> both (2);
        both[1].phaseDeg = 180.0;
        const std::vector<std::vector<float>> twin { irs[0], irs[0] };
        fb::MasterParams noF;                                          // reference: no master filter
        fb::MasterParams lpF; lpF.lpf = { true, 1000.0, 24 };
        const auto a = fb::makeOverlay (twin, both, noF, sr);
        const auto b = fb::makeOverlay (twin, both, lpF, sr);
        bool shrunk = a.interference.size() == b.interference.size();
        for (std::size_t p = 0; shrunk && p < a.interference.size(); ++p)
            shrunk = std::abs (b.interference[p]) <= std::abs (a.interference[p]) + 1e-9
                     && b.interference[p] * a.interference[p] >= -1e-12;   // same sign (or zero)
        ok (shrunk, "|H| in [0,1] only attenuates the interference column");

        std::vector<fb::StripParams> inPhase (2);                      // an AUDIBLE blend for the rolloff check
        const auto c = fb::makeOverlay (twin, inPhase, noF, sr);       // (the flipped pair above is near-silence)
        const auto d = fb::makeOverlay (twin, inPhase, lpF, sr);
        ok (d.blend.back() < c.blend.back() - 30.0, "the blend curve itself shows the Master rolloff");
    }

    return felitronics::test::report();
}
