// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// K3c — band bursts on Mid and Side, and the cross-reading that says which one carried the burst.
//
// THE ORACLE IS THE SCALED-COPY THEOREM, not a second detector. If R = k·L then M = (1+k)/2·L and
// S = (1−k)/2·L: BOTH axes are scaled copies of the same signal, so a detector that judges a hop against
// its OWN baseline must fire identically on both, at the same hops, whatever k is — while the cross
// reading between them must be exactly ((1+k)/(1−k))² in power. One fixture therefore pins the event
// logic (counts and positions must MATCH) and the cross arithmetic (a number that must DIFFER by a
// closed form) at the same time, and neither number is one this file measured and wrote down.

#include <felitronics_test.h>

#include <felitronics/analysis/StereoBandBursts.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
using felitronics::analysis::BandBursts;
using felitronics::analysis::BandBurstsHopTrace;
using felitronics::analysis::StereoBandBursts;
namespace test = felitronics::test;
using test::ok;
using test::approx;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 48000.0;

// A quiet continuous tone in the band so the baseline is eligible and non-zero, with ten raised-cosine
// bursts 20 dB over it, all after the 2 s baseline has filled.
std::vector<float> bandBursts (double seconds = 6.0, double hz = 7000.0)
{
    const std::size_t n = (std::size_t) (kFs * seconds);
    std::vector<float> x (n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / kFs;
        double a = 0.02;
        for (int b = 0; b < 10; ++b)
        {
            const double t0 = 3.0 + 0.25 * b, d = 0.04;
            if (t >= t0 && t < t0 + d)
            {
                const double u = (t - t0) / d;                       // a smooth envelope, so the burst's
                a += 0.2 * 0.5 * (1.0 - std::cos (2.0 * kPi * u));   // edges do not spray broadband energy
            }
        }
        x[i] = (float) (a * std::sin (2.0 * kPi * hz * t));
    }
    return x;
}

bool feed2 (StereoBandBursts& sb, const std::vector<float>& l, const std::vector<float>& r)
{
    const float* ch[2] { l.data(), r.data() };
    if (! sb.process (ch, 2, (int) l.size())) return false;
    return sb.finish();
}
} // namespace

int main()
{
    std::printf ("felitronics::analysis StereoBandBursts — K3c\n");

    test::group ("centred content: Side is exactly zero, and the instrument says so without a threshold");
    {
        const auto x = bandBursts();
        StereoBandBursts sb;
        if (test::run (sb.prepare (kFs, 2)) && test::run (feed2 (sb, x, x)))
        {
            ok (sb.mid().eventCount() == 10,
                "Mid carries all ten bursts (got " + std::to_string (sb.mid().eventCount()) + ")");
            ok (sb.side().eventCount() == 0, "Side carries none");
            // THE ABSENCE IS REPORTED AS A STRUCTURE, NOT A SMALL NUMBER. An exactly zero Side has no
            // eligible hop at all, which is the engine's own existing "there was nothing to judge".
            ok (! sb.side().eventsValid(), "…and Side's report refuses to call itself valid");
            const std::int64_t hops = sb.exactZeroHops (StereoBandBursts::kSide);
            ok (hops > 0 && hops == sb.side().hopCount(),
                "every Side hop was exactly zero (" + std::to_string (hops) + " of "
                    + std::to_string (sb.side().hopCount()) + ")");
            ok (sb.exactZeroHops (StereoBandBursts::kMid) == 0, "and no Mid hop was");
            // Two bit-identical channels are NOT the absent case: there are two of them.
            ok (! sb.sideAbsent(), "two identical channels are still two channels, not a missing Side");
        }
    }

    test::group ("a scaled copy: both axes must fire identically, and differ by a closed form");
    {
        // 0.18 dB of imbalance — small enough that a reader looking at levels would call it centred.
        const double k = std::pow (10.0, -0.18 / 20.0);
        const auto l = bandBursts();
        std::vector<float> r (l.size());
        for (std::size_t i = 0; i < l.size(); ++i) r[i] = (float) (k * (double) l[i]);

        StereoBandBursts sb;
        if (test::run (sb.prepare (kFs, 2)) && test::run (feed2 (sb, l, r)))
        {
            ok (sb.mid().eventCount() == 10 && sb.side().eventCount() == 10,
                "both axes fire on all ten (Mid " + std::to_string (sb.mid().eventCount()) + ", Side "
                    + std::to_string (sb.side().eventCount()) + ")");

            bool sameHops = sb.mid().storedEventCount() == sb.side().storedEventCount();
            for (std::int64_t e = 0; sameHops && e < sb.mid().storedEventCount(); ++e)
                sameHops = sb.mid().event (e).start == sb.side().event (e).start;
            ok (sameHops, "at the same hops — a relative detector cannot tell a scaled copy from its original");

            // THE CROSS ARITHMETIC. At Mid's peak hop, Side's band mean square is ((1-k)/2)^2 of the same
            // L band power that Mid holds ((1+k)/2)^2 of, so the ratio is fixed by k alone.
            const double want = 20.0 * std::log10 ((1.0 + k) / (1.0 - k));
            const auto cross = sb.crossAt (StereoBandBursts::kMid, 0);
            ok (cross.eligible && cross.hop >= 0, "the cross reading was captured and was eligible");
            if (cross.eligible && cross.power > 0.0)
            {
                const double got = 10.0 * std::log10 (sb.mid().event (0).peakPower / cross.power);
                approx (got, want, 1e-3,
                        "Mid stands " + std::to_string (want) + " dB over Side at its own peak");
            }
            // AND THE OTHER DIRECTION, which is where this class first shipped broken. Mid is traced
            // before Side within one input sample, so a join that only fills from a LATER trace answers
            // Mid and leaves every Side event reading power 0, eligible false, hop -1 — indistinguishable
            // from "Mid was silent there", and the exact opposite of the truth. Asserting only the Mid
            // direction is what let it pass.
            const auto back = sb.crossAt (StereoBandBursts::kSide, 0);
            ok (back.eligible && back.hop >= 0, "Side's own cross reading was captured too");
            if (back.eligible && sb.side().event (0).peakPower > 0.0)
            {
                const double got = 10.0 * std::log10 (back.power / sb.side().event (0).peakPower);
                approx (got, want, 1e-3,
                        "…and reads Mid standing the same " + std::to_string (want) + " dB over it");
            }
            // Structural, over EVERY stored event on both axes: a cross reading that was never filled is
            // the failure above, and one row of it is enough to make a page draw the wrong conclusion.
            std::int64_t unfilled = 0;
            for (int ax = 0; ax < StereoBandBursts::kAxes; ++ax)
                for (std::int64_t e = 0; e < sb.axis (ax).storedEventCount(); ++e)
                    if (sb.crossAt (ax, e).hop < 0) ++unfilled;
            ok (unfilled == 0, "no stored event on either axis is missing its cross hop (got "
                                   + std::to_string (unfilled) + " unfilled)");

            // The same fixture at 20 dB of imbalance must move that number and nothing else.
            ok (sb.exactZeroHops (StereoBandBursts::kSide) == 0,
                "and no Side hop is exactly zero any more — 0.18 dB is a real Side, not an absent one");
        }
    }

    test::group ("absence is one channel, and only one channel");
    {
        const auto x = bandBursts();
        StereoBandBursts sb;
        const float* ch[1] { x.data() };
        if (test::run (sb.prepare (kFs, 2)) && test::run (sb.process (ch, 1, (int) x.size()))
            && test::run (sb.finish()))
        {
            ok (sb.sideAbsent(), "a one-channel programme has no Side to report");
            ok (sb.mid().eventCount() == 10, "while Mid is the programme itself");
        }
    }

    // THE DOME, AGAINST THE FILTER AND NOT AGAINST ITS OWN FORMULA. The closed form b^8/(a^2+b^2)^4 is
    // derived from the LR4 power response; checking it by scanning that same response would only certify
    // the maximisation. So the tone is actually RUN through the detector and the band/wideband ratio is
    // read off the hop trace — an end-to-end measurement of the thing the number claims.
    test::group ("the dome: what the band can report at all, measured through the band");
    {
        StereoBandBursts probe;
        if (test::run (probe.prepare (kFs, 2)))
        {
            const double f = probe.domeHz(), want = 10.0 * std::log10 (probe.domeShare());
            ok (f > 0.0 && f > 5000.0 && f < 9000.0,
                "the dome sits inside the band at " + std::to_string (f) + " Hz");

            auto gainAt = [] (double hz)
            {
                struct Acc { double band = 0.0, wide = 0.0; } acc;
                BandBursts bb;
                if (! bb.prepare (kFs, 0, 1)) return -1000.0;
                bb.setObserver ([] (void* u, const BandBurstsHopTrace& t) noexcept
                                {
                                    auto* a = static_cast<Acc*> (u);
                                    if (t.full && t.startSample >= (std::int64_t) kFs)  // past the settling
                                    { a->band += t.energy; a->wide += t.wideEnergy; }
                                }, &acc);
                std::vector<float> s ((std::size_t) (kFs * 3.0));
                for (std::size_t i = 0; i < s.size(); ++i)
                    s[i] = (float) (0.5 * std::sin (2.0 * kPi * hz * (double) i / kFs));
                const float* c[1] { s.data() };
                if (! bb.process (c, 1, (int) s.size())) return -1000.0;
                bb.finish();
                return acc.wide > 0.0 ? 10.0 * std::log10 (acc.band / acc.wide) : -1000.0;
            };

            const double measured = gainAt (f);
            approx (measured, want, 1e-3,   // the end-to-end residual is under 1e-4 dB; a geometric-vs-arithmetic mean slips 0.072
                    "the band passes " + std::to_string (want) + " dB there, as the closed form says");
            // …and it is a MAXIMUM, which is the half a single-point check cannot see.
            ok (gainAt (f * 0.6) < measured - 1.0 && gainAt (f * 1.6) < measured - 1.0,
                "and less on either side of it");
        }
    }

    return test::report();
}
