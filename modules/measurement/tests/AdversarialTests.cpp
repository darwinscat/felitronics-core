// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// ADVERSARIAL self-tests for felitronics::measurement — one per bug the crew ("break it" consilium:
// deepseek + antigravity + Fable, sanitizer-confirmed) found at EXTREME/edge inputs. Each case fed a
// crash / UB / non-finite / silent-garbage input before the fix and now must run cleanly + return a
// sane, finite result. Run this suite under ASan/UBSan (FELITRONICS_ENABLE_SANITIZERS=ON) to also gate
// the heap-overflow (B2) and float→size_t (B6) fixes.

#include <felitronics_test.h>
#include <felitronics/measurement/Convolve.h>
#include <felitronics/measurement/Sweep.h>
#include <felitronics/measurement/CaptureGate.h>
#include <felitronics/measurement/IrPost.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

namespace
{
bool allFinite (const std::vector<double>& x)
{
    for (double v : x) if (! std::isfinite (v)) return false;
    return true;
}
const double kInf = std::numeric_limits<double>::infinity();
const double kNaN = std::numeric_limits<double>::quiet_NaN();
}

int main()
{
    std::printf ("felitronics::measurement adversarial tests\n");

    // B8 — nextPow2 must never hang; returns 0 past the largest representable power of two.
    group ("B8 nextPow2 overflow guard");
    {
        ok (measurement::nextPow2 (0) == 1, "nextPow2(0)=1");
        ok (measurement::nextPow2 (5) == 8, "nextPow2(5)=8");
        ok (measurement::nextPow2 (std::size_t (1) << 63) == (std::size_t (1) << 63), "nextPow2(2^63)=2^63");
        ok (measurement::nextPow2 ((std::size_t (1) << 63) + 1) == 0, "nextPow2(>2^63)=0 (no hang)");
    }

    // B2 — magSpectrum on a non-power-of-two nfft must NOT heap-overflow; it rounds nfft up.
    group ("B2 magSpectrum non-pow2 heals");
    {
        const std::vector<double> x { 1.0, 2.0, 3.0 };
        const auto m = measurement::magSpectrum (x, 3);       // was: heap-buffer-overflow in radix-2
        ok (m.size() == 2, "nfft 3 -> rounded to 4 -> 2 bins");
        ok (allFinite (m), "spectrum finite");
    }

    // B3 — makeSweep must sanitize EVERY parameter (no overflow / NaN / inf leaking through).
    group ("B3 makeSweep parameter sanitization");
    {
        measurement::SweepSpec s; s.durationSeconds = 1e18;   // was: llround overflow -> vector length_error
        const auto sw = measurement::makeSweep (s);
        ok (! sw.signal.empty() && sw.signal.size() < (std::size_t) (130.0 * sw.spec.sampleRate), "huge duration clamped");
        ok (allFinite (sw.signal) && allFinite (sw.inverse), "signal + inverse finite");

        measurement::SweepSpec t; t.tailSeconds = kNaN; t.fadeSeconds = kNaN;   // was: NaN through max/clamp
        const auto sw2 = measurement::makeSweep (t);
        ok (std::isfinite (sw2.spec.tailSeconds) && std::isfinite (sw2.spec.fadeSeconds), "NaN tail/fade healed");
        ok (allFinite (sw2.signal), "signal finite with NaN fade/tail");

        measurement::SweepSpec a; a.amplitude = kInf;         // was: inf>=0 passes -> non-finite signal
        ok (allFinite (measurement::makeSweep (a).signal), "amplitude +inf healed");

        measurement::SweepSpec z; z.amplitude = 0.0;          // was: 0 -> silent zero sweep
        const auto sz = measurement::makeSweep (z);
        double pk = 0.0; for (double v : sz.signal) pk = std::max (pk, std::fabs (v));
        ok (pk > 0.0, "amplitude 0 healed to a real sweep");
    }

    // B4 — f1 near Nyquist must NOT invert the f2 clamp (std::clamp lo>hi is UB) → ascending sweep, harmonicL>0.
    group ("B4 f1 near Nyquist: no clamp inversion");
    {
        measurement::SweepSpec s; s.sampleRate = 48000.0; s.f1 = 23990.0; s.f2 = 23995.0;
        const auto sw = measurement::makeSweep (s);           // was: f2<f1, R<0, harmonicL=-8565
        ok (sw.spec.f2 > sw.spec.f1, "still an ASCENDING sweep");
        ok (sw.harmonicL > 0.0 && std::isfinite (sw.harmonicL), "harmonicL positive + finite");
        ok (sw.spec.f2 <= 0.5 * sw.spec.sampleRate, "f2 <= Nyquist");
    }

    // B5 — a finite-but-huge amplitude must NOT overflow the Farina inverse to inf.
    group ("B5 huge finite amplitude");
    {
        measurement::SweepSpec s; s.amplitude = 1e308;        // was: inverse peak inf -> poisons everything
        ok (allFinite (measurement::makeSweep (s).inverse), "inverse stays finite");
    }

    // B6 — bandEnergy with negative / NaN / >Nyquist band edges must NOT do a UB float→size_t cast.
    group ("B6 bandEnergy degenerate band edges");
    {
        std::vector<double> x (8192, 0.0);
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::sin (0.05 * (double) i);
        ok (std::isfinite (measurement::bandEnergy (x, 48000.0, -100.0, 500.0)), "negative loHz -> finite");
        ok (std::isfinite (measurement::bandEnergy (x, 48000.0, 100.0, kNaN)),    "NaN hiHz -> finite");
        ok (std::isfinite (measurement::bandEnergy (x, 48000.0, 100.0, 1e9)),     "hiHz > Nyquist -> finite");
        ok (measurement::bandEnergy (x, 48000.0, 2000.0, 100.0) == 0.0,           "lo>hi -> 0 (no wrap)");
    }

    // B1 — gateRecording on a malformed sweep (empty inverse) must REJECT, not crash on d[0].
    group ("B1 gate on empty-inverse sweep rejects (no OOB)");
    {
        const std::vector<double> rec (1000, 0.1);
        const auto g = measurement::gateRecording (rec, measurement::Sweep {});   // was: OOB read d[0]
        ok (! g.ok && g.reason == measurement::GateReject::SweepNotDetected, "rejected, not crashed");
    }

    // B7 — snrValid distinguishes "too short to measure" from a real 0 dB.
    group ("B7 snrValid flag");
    {
        measurement::SweepSpec s; s.durationSeconds = 0.3; s.f1 = 200.0; s.f2 = 4000.0; s.tailSeconds = 0.05;
        const auto sw = measurement::makeSweep (s);           // signal ~16800 samples (>= 4096)
        std::vector<double> shortRec (3000, 0.0);             // < 4096 → SNR path skipped
        for (std::size_t i = 0; i < shortRec.size() && i < sw.signal.size(); ++i) shortRec[i] = sw.signal[i];
        ok (! measurement::gateRecording (shortRec, sw).snrValid, "short recording -> snrValid=false");
        ok (measurement::gateRecording (std::vector<double> (sw.signal.begin(), sw.signal.end()), sw).snrValid,
            "full recording -> snrValid=true");
    }

    // B9 — the double FFT (exact twiddles) is transparent: sweep⊛inverse floor FAR from the delta is tiny.
    group ("B9 double FFT numerical floor");
    {
        measurement::SweepSpec s; s.f1 = 100.0; s.f2 = 8000.0; s.durationSeconds = 0.3; s.tailSeconds = 0.05;
        const auto sw = measurement::makeSweep (s);
        std::vector<double> proper (sw.signal.begin(), sw.signal.begin() + (std::ptrdiff_t) sw.sweepLen);
        const auto d = measurement::convolve (proper, sw.inverse);
        std::size_t pk = 0; double peak = 0.0;
        for (std::size_t i = 0; i < d.size(); ++i) { const double a = std::fabs (d[i]); if (a > peak) { peak = a; pk = i; } }
        double floorFar = 0.0;                                // > 4000 samples from the delta = pure numerical floor
        for (std::size_t i = 0; i < d.size(); ++i)
            if ((i + 4000 < pk) || (i > pk + 4000)) floorFar = std::max (floorFar, std::fabs (d[i]));
        ok (floorFar < 1e-6 * peak, "numerical floor < -120 dB far from the delta");
    }

    return felitronics::test::report();
}
