// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// codecgrid::CodecGridDetector's VERDICTS, on programmes made to have exactly one property each: a grid at the
// programme's rate, a grid at the other rate, no grid, a grid too weak for the rule but on one phase everywhere,
// a grid in one stretch only, a grid in half the programme, material cut into pieces with a phase each. The
// programmes carry real transform grids (SyntheticCodec.h); no codec and no audio file is involved.

#include "DetectorHarness.h"
#include "felitronics_test.h"

#include <algorithm>
#include <cstdio>
#include <string>

using namespace harness;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
void print (const char* what, const Detector& d)
{
    const auto& r = d.result();
    std::printf ("    %-34s %-12s %-14s rate %5d phase %4d | found %d of %d, on the phase %d | best %.1f | %d scans\n", what, name (r.verdict), name (r.ground),
                 r.codecRate, r.gridPhase, r.windowsFound, r.windowsExamined, r.windowsAgreeing, r.bestScore, scans (d));
}
}

int main()
{
    std::printf ("felitronics::codecgrid detector verdict tests\n");
    constexpr int kSeconds = 3;
    const int n44 = 44100 * kSeconds, n48 = 48000 * kSeconds;
    const auto original44 = synthetic::programme (n44, 7u);
    const auto original48 = synthetic::programme (n48, 8u);

    group ("a grid at the programme's rate is confirmed, with its family, its rate and its phase");
    {
        struct Case { Transform t; int offset; double rate; Family family; int hop; };
        for (const Case c : { Case { Transform::AacSine, 480, 44100.0, Family::Aac, 1024 }, Case { Transform::Mp3, 133, 44100.0, Family::Mp3, 576 },
                              Case { Transform::Celt, 548, 48000.0, Family::Celt, 960 }, Case { Transform::AacKbd, 1023, 48000.0, Family::Aac, 1024 } })
        {
            const auto x = synthetic::coded (c.rate > 45000.0 ? original48 : original44, c.t, c.offset, 0.3);
            const auto d = analyse (x, c.rate, params());
            const std::string tag = std::string (c.t == Transform::Mp3 ? "MP3" : c.t == Transform::Celt ? "CELT" : c.t == Transform::AacSine ? "AAC sine" : "AAC KBD")
                                  + " at " + std::to_string (c.offset) + ": ";
            if (! felitronics::test::run (d != nullptr)) continue;
            print (tag.c_str(), *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::Confirmed && r.reason == Reason::None, tag + "confirmed");
            ok (r.ground == Ground::Both, tag + "by both grounds: the rule in two stretches, and the phase in every stretch that was followed");
            ok (r.family == c.family && r.codecRate == (int) c.rate && r.hop == c.hop, tag + "the family, the rate and the hop");
            ok (r.transform == c.t, tag + "the transform that read it best is the one it was coded with");
            ok (r.gridPhase == c.offset, tag + "the phase is the offset it was coded at, counted from sample 0 of the programme");
            ok (r.windows == 8 && r.windowsExamined == 8 && r.windowsFound == 8 && r.windowsAgreeing == 8, tag + "eight stretches placed, examined, found, and on the phase");
            ok (r.bestWindow >= 0 && r.bestHypothesis >= 0 && d->hypothesis (r.bestHypothesis).transform == c.t && r.bestOffset >= 0 && r.bestOffset < c.hop
                    && r.bestScore > 100.0, tag + "the published reading is one of this transform, with a score far over the rule's");
            ok (d->curveOffsets() == c.hop && d->zeroMapRows() > 5, tag + "its curve and its zero map are published");
            ok (r.zeroShare > 0.2f && r.zeroShare < 0.5f, tag + "the zero share is about the 30 % that were zeroed");
            // the curve really is that reading's: its deepest summed offset is bestOffset
            int deepest = 0;
            double least = 1.0e300;
            for (int o = 0; o < d->curveOffsets(); ++o)
            {
                double sum = 0.0;
                for (int k = 0; k < kCells; ++k) sum += (double) d->curve()[(std::size_t) o * kCells + (std::size_t) k];
                if (sum < least) { least = sum; deepest = o; }
            }
            ok (deepest == r.bestOffset, tag + "the published curve is deepest at the published offset");
            // Follow: the first stretch tried everything, the others only the confirmed family
            const int perWindow = c.family == Family::Aac ? 2 : 1;
            ok (scans (*d) == kHypotheses + 7 * perWindow, tag + "seven hypotheses in the first stretch, then only the family's in the other seven");
            bool phases = true;
            for (int w = 0; w < 8; ++w)
                for (int h = 0; h < kHypotheses; ++h)
                {
                    const auto& hr = d->window (w).hypotheses[h];
                    if (hr.found) phases = phases && hr.gridPhase == c.offset && familyOf (d->hypothesis (h).transform) == c.family;
                }
            ok (phases, tag + "every found reading is of the family, on the phase");
        }
    }

    group ("a grid at the OTHER rate is confirmed after the stretch is taken back to it");
    {
        struct Case { Transform t; int offset, codecRate, fileRate; Family family; const char* what; };
        for (const Case c : { Case { Transform::AacKbd, 300, 48000, 44100, Family::Aac, "AAC coded at 48 kHz, in a 44.1 kHz file" },
                              Case { Transform::Celt, 548, 48000, 44100, Family::Celt, "CELT (always 48 kHz), in a 44.1 kHz file" },
                              Case { Transform::Mp3, 133, 44100, 48000, Family::Mp3, "MP3 coded at 44.1 kHz, in a 48 kHz file" } })
        {
            const auto x = codedAtOtherRate (c.codecRate * kSeconds, 11u, c.t, c.offset, c.codecRate, c.fileRate);
            const auto d = analyse (x, (double) c.fileRate, params (Depth::Verdict));
            const std::string tag = std::string (c.what) + ": ";
            if (! felitronics::test::run (d != nullptr)) continue;
            print (c.what, *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::Confirmed && r.ground == Ground::TwoWindows, tag + "confirmed by the rule in two stretches");
            ok (r.family == c.family && r.codecRate == c.codecRate && r.transform == c.t, tag + "the family, the CODEC's rate and the transform");
            ok (r.gridPhase == c.offset, tag + "the phase is the offset at the codec's rate: the conversion back kept sample 0 on sample 0");
            ok (scans (*d) == kHypotheses + 1, tag + "depth Verdict stops at the second stretch's first hypothesis: eight scans");
            ok (r.windowsFound == 2 && r.windowsExamined == 8, tag + "found in the two stretches that were read, of the eight that could be");
        }
    }

    group ("no grid: every hypothesis in every stretch, and nothing found");
    {
        const auto d = analyse (original44, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("the programme uncoded", *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::NoGrid && r.reason == Reason::None && r.ground == Ground::None, "NoGrid");
            ok (r.family == Family::None && r.codecRate == 0 && r.hop == 0 && r.windowsFound == 0 && r.windowsAgreeing == 0, "no family, no rate, nothing found");
            ok (scans (*d) == 8 * kHypotheses, "all 56 scans were made: that is what saying `no grid` costs");
            ok (r.bestWindow >= 0 && r.bestScore > 0.0 && r.bestScore < 9.2 && d->curveOffsets() > 0, "the best reading there was is still published, under the rule's score");
            int most = 0;
            for (int h = 0; h < kHypotheses; ++h)
                for (int w = 0; w < 8; ++w)
                {
                    int count = 0;
                    for (int v = 0; v < 8; ++v)
                        if (phasesAgree (d->window (v).hypotheses[h].gridPhase, d->window (w).hypotheses[h].gridPhase, hopOf (d->hypothesis (h).transform), 1)) ++count;
                    most = std::max (most, count);
                }
            std::printf ("    the most stretches any hypothesis has on one phase: %d\n", most);
            ok (most < 4, "no hypothesis has four stretches on one phase");
        }
    }

    group ("depth Exhaustive reads everything even after the rule has confirmed");
    {
        const auto d = analyse (synthetic::coded (original44, Transform::AacSine, 480, 0.3), 44100.0, params (Depth::Exhaustive));
        if (felitronics::test::run (d != nullptr))
        {
            print ("AAC, exhaustive", *d);
            ok (d->result().verdict == Verdict::Confirmed && d->result().ground == Ground::Both && d->result().gridPhase == 480, "the same verdict");
            ok (scans (*d) == 8 * kHypotheses, "and all 56 scans");
            bool others = true;
            for (int w = 0; w < 8; ++w)
                for (int h = 0; h < kHypotheses; ++h)
                    if (! (familyOf (d->hypothesis (h).transform) == Family::Aac && d->hypothesis (h).codecRate == 44100)) others = others && ! d->window (w).hypotheses[h].found;
            ok (others, "no other family and no other rate finds anything in an AAC programme");
        }
    }

    group ("a grid too weak for the rule, on one phase in every stretch: confirmed by the phase alone");
    {
        // the coded programme under another, uncoded one at half its level
        auto x = synthetic::coded (original44, Transform::AacSine, 480, 0.3);
        const auto other = synthetic::programme (n44, 99u);
        for (int i = 0; i < n44; ++i)
        {
            x.left[(std::size_t) i] += 0.5f * other.left[(std::size_t) i];
            x.right[(std::size_t) i] += 0.5f * other.right[(std::size_t) i];
        }
        const auto d = analyse (x, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("AAC under an uncoded programme", *d);
            const auto& r = d->result();
            ok (r.windowsFound == 0, "the rule finds it in no stretch");
            ok (r.verdict == Verdict::Confirmed && r.ground == Ground::PhaseAgreement, "and it is confirmed, by the phase ground");
            ok (r.family == Family::Aac && r.codecRate == 44100 && r.gridPhase == 480 && r.windowsAgreeing >= 4, "the family, the rate and the phase are the grid's");
            ok (scans (*d) == 8 * kHypotheses, "it cost every scan: the phase ground has no early exit");
        }
    }

    group ("a grid in ONE stretch is `InPlaces`, not a finding");
    {
        const auto g = Detector::geometryFor (44100.0, kStretchSeconds);
        const std::uint64_t start = Detector::windowStart (g, (std::uint64_t) n44, 3);
        const auto coded = synthetic::coded (original44, Transform::AacSine, 480, 0.3);
        const auto x = spliced (original44, coded, (std::size_t) start, (std::size_t) start + (std::size_t) g.length);
        const auto d = analyse (x, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("AAC in the fourth stretch only", *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::InPlaces && r.ground == Ground::None, "InPlaces, with no ground");
            ok (r.family == Family::Aac && r.codecRate == 44100 && r.gridPhase == 480 && r.windowsFound == 1, "what was found, and that it was found once");
            ok (d->window (3).best >= 0 && d->window (2).best < 0 && d->window (4).best < 0, "it is the fourth stretch, and not its neighbours");
            ok (r.bestWindow == 3, "the published reading is that stretch's");
        }
    }

    group ("a grid in half the programme: confirmed, and the stretches say where");
    {
        const auto coded = synthetic::coded (original44, Transform::Mp3, 133, 0.3);
        const auto x = spliced (original44, coded, 0, (std::size_t) n44 / 2);
        const auto d = analyse (x, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("MP3 in the first half", *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::Confirmed && r.family == Family::Mp3 && r.gridPhase == 133, "confirmed, the family and the phase");
            bool where = true;
            for (int w = 0; w < 8; ++w)
            {
                const bool inside = d->window (w).startFrame + (std::uint64_t) Detector::geometryFor (44100.0, kStretchSeconds).length <= (std::uint64_t) n44 / 2u;
                const bool outside = d->window (w).startFrame >= (std::uint64_t) n44 / 2u;
                if (inside) where = where && d->window (w).best >= 0;
                if (outside) where = where && d->window (w).best < 0;
            }
            ok (where, "every stretch inside the coded half is found, none outside it");
            ok (r.windowsFound >= 3 && r.windowsFound <= 5 && r.windowsFound < r.windowsExamined, "found in about half of the stretches, not in all");
        }
    }

    group ("material cut into pieces, each with its own phase: several grids, not one");
    {
        synthetic::Stereo x = original44;
        const int piece = 3307;                                 // 75 ms
        const int offsets[5] { 0, 277, 554, 831, 84 };
        synthetic::Stereo codedAt[5];
        for (int i = 0; i < 5; ++i) codedAt[i] = synthetic::coded (original44, Transform::AacSine, offsets[i], 0.3);
        for (int k = 0; k * piece < n44; ++k) x = spliced (x, codedAt[k % 5], (std::size_t) (k * piece), (std::size_t) ((k + 1) * piece));
        const auto d = analyse (x, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("AAC in 75 ms pieces", *d);
            const auto& r = d->result();
            ok (r.verdict == Verdict::SeveralGrids && r.ground == Ground::None, "SeveralGrids: broad and deep in every stretch, unique in none");
            ok (r.family == Family::Aac && r.windowsFound == 0, "the family is named; the rule found nothing");
            bool several = true;
            for (int w = 0; w < 8; ++w)
            {
                bool any = false;
                for (int h = 0; h < kHypotheses; ++h) any = any || d->window (w).hypotheses[h].several;
                several = several && any;
            }
            ok (several, "every stretch has a reading that is `several`");
        }
    }

    group ("a mono programme is read as two identical channels");
    {
        synthetic::Stereo mono = synthetic::coded (original44, Transform::Mp3, 200, 0.3);
        mono.right = mono.left;
        const auto d = analyse (mono, 44100.0, params (Depth::Verdict), 1);
        const auto s = analyse (mono, 44100.0, params (Depth::Verdict), 2);
        if (felitronics::test::run (d != nullptr && s != nullptr))
        {
            print ("MP3, one channel", *d);
            ok (d->result().verdict == Verdict::Confirmed && d->result().family == Family::Mp3 && d->result().gridPhase == 200, "confirmed from one channel");
            ok (sameReport (*d, *s), "and the report is the one two identical channels give");
        }
    }

    group ("what cannot be examined is reported, not refused");
    {
        Detector d;
        d.setParams (params());
        const float* in[2] { original44.left.data(), original44.right.data() };
        // an unsupported rate
        ok (d.prepare (96000.0, 2, (std::uint64_t) n44), "96 kHz is accepted by prepare()");
        ok (d.windows() == 0 && d.process (in, 2, n44), "...no stretch is placed, and the programme is taken");
        ok (d.finish() && d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::UnsupportedRate, "NotExamined, UnsupportedRate");
        ok (d.progress() >= 1.0 && d.curveOffsets() == 0 && d.result().bestWindow < 0, "finished, with no curve");
        // too short: one sample under what one stretch and its margins need
        const auto g = Detector::geometryFor (44100.0, kStretchSeconds);
        const std::uint64_t least = (std::uint64_t) g.length + 2u * (std::uint64_t) g.margin + (std::uint64_t) g.grid;
        ok (d.prepare (44100.0, 2, least - 1) && d.windows() == 0, "a programme one sample too short is accepted and has no stretch");
        ok (d.process (in, 2, (int) least - 1) && d.finish() && d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::TooShort, "NotExamined, TooShort");
        ok (d.prepare (44100.0, 2, least) && d.windows() == 1, "the shortest programme that can be examined has one stretch");
        // an empty programme
        ok (d.prepare (44100.0, 2, 0) && d.finish() && d.result().verdict == Verdict::NotExamined && d.result().reason == Reason::TooShort, "an empty programme: TooShort");
        // silence
        synthetic::Stereo z { std::vector<float> ((std::size_t) n44, 0.0f), std::vector<float> ((std::size_t) n44, 0.0f) };
        const auto silent = analyse (z, 44100.0, params());
        if (felitronics::test::run (silent != nullptr))
        {
            ok (silent->result().verdict == Verdict::NotExamined && silent->result().reason == Reason::Silent, "digital silence: NotExamined, Silent");
            ok (silent->result().windows == 8 && silent->result().windowsExamined == 0 && scans (*silent) == 0, "eight stretches placed, none examined, no scan made");
        }
        // a programme whose only loud stretch is coded: one stretch examined
        const auto coded = synthetic::coded (original44, Transform::AacSine, 480, 0.3);
        const std::uint64_t start = Detector::windowStart (g, (std::uint64_t) n44, 5);
        const auto x = spliced (z, coded, (std::size_t) start - 2000u, (std::size_t) start + (std::size_t) g.length + 2000u);
        const auto one = analyse (x, 44100.0, params());
        if (felitronics::test::run (one != nullptr))
        {
            print ("one loud stretch in silence", *one);
            ok (one->result().windowsExamined == 1 && one->result().verdict == Verdict::InPlaces && one->window (5).best >= 0, "one stretch examined, found, and no more than InPlaces");
            ok (one->window (0).silent && ! one->window (5).silent, "the silent stretches are marked silent");
        }
    }

    group ("a short programme has fewer stretches, and they never share samples");
    {
        for (double seconds : { 0.8, 1.1, 1.6, 2.6, 3.0 })
        {
            const std::uint64_t total = (std::uint64_t) (seconds * 44100.0);
            const auto g = Detector::geometryFor (44100.0, kStretchSeconds);
            const auto st = Detector::storageFor (44100.0, 2, total, params());
            Detector d;
            d.setParams (params());
            if (! felitronics::test::run (d.prepare (44100.0, 2, total))) continue;
            // the stretches are those of the storage, in order, on the grid, inside the programme, a stretch apart.
            // One finish step is enough to read where they are: it places the report's stretches and scans nothing.
            bool sane = d.windows() == st.windows && d.windows() >= 1 && d.windows() <= 8;
            synthetic::Stereo x = synthetic::programme ((int) total, 5u);
            sane = sane && feed (d, x, 2) && ! d.finishStep();
            std::uint64_t previous = 0;
            for (int w = 0; w < d.windows(); ++w)
            {
                const std::uint64_t s = d.window (w).startFrame;
                sane = sane && s % (std::uint64_t) g.grid == 0 && s >= (std::uint64_t) g.margin && s + (std::uint64_t) g.length + (std::uint64_t) g.margin <= total;
                if (w > 0) sane = sane && s >= previous + (std::uint64_t) g.length;
                previous = s;
            }
            std::printf ("    %.1f s: %d stretches\n", seconds, d.windows());
            ok (sane, std::to_string (seconds) + " s: the stretches are on the grid, inside the programme, and do not overlap");
        }
        // 1.6 s coded: fewer than eight stretches, still confirmed by the rule
        const int n = (int) (1.6 * 44100.0);
        const auto x = synthetic::coded (synthetic::programme (n, 5u), Transform::AacSine, 480, 0.3);
        const auto d = analyse (x, 44100.0, params());
        if (felitronics::test::run (d != nullptr))
        {
            print ("AAC, 1.6 s", *d);
            ok (d->windows() < 8 && d->windows() >= 2 && d->result().verdict == Verdict::Confirmed && d->result().gridPhase == 480, "a short coded programme is confirmed from the stretches it has");
        }
    }

    return felitronics::test::report();
}
