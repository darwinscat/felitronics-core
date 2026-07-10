// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Headless property tests for felitronics::blend — ported from OrbitCapture's app-local
// BlendEngineTests.cpp + BlendKernelsTests.cpp (de-monolith steps 2-3 → promoted step 4). INDEPENDENT
// property checks (not a re-run of the engine): default params must equal a plain sum that bypasses the
// kernels; polarity must cancel; solo/mute/gain/master semantics; "filter parked at its extreme ≡ off";
// the kernels' own DC/−3dB/polarity/integer-shift behaviour; plus a new check for the odd-order
// applyBlendSlope guard added when promoting into core (no valid caller ever passed an odd order, so
// this is a rejects-to-no-op, not a behaviour change for any real caller).
//
// NOTE: BlendKernelsTests.cpp's "convolveDI: δ is identity" and "resampleIR DC preservation" groups are
// NOT ported — convolveDI/resampleIR are app-local audition/resample helpers, not among the 5 kernels
// promoted here (biquadInplace, onePoleInplace, applyBlendSlope, rotatePhase, shiftFrac).
#include <felitronics_test.h>
#include <felitronics/blend/Blend.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <vector>

using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace
{
std::vector<float> ramp (std::size_t n, float scale = 1.0f)
{
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = scale * std::sin (0.03f * (float) i + 0.5f);
    return x;
}
double maxAbsDiff (const std::vector<float>& a, const std::vector<float>& b)
{
    double d = 0.0; for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) d = std::max (d, (double) std::fabs (a[i] - b[i]));
    return d;
}
double rms (const std::vector<float>& x, std::size_t from, std::size_t to)
{
    double s = 0.0; std::size_t c = 0;
    for (std::size_t i = from; i < to && i < x.size(); ++i) { s += (double) x[i] * x[i]; ++c; }
    return c ? std::sqrt (s / (double) c) : 0.0;
}
std::vector<float> sine (double f, double sr, std::size_t n)
{
    std::vector<float> x (n);
    for (std::size_t i = 0; i < n; ++i) x[i] = (float) std::sin (2.0 * felitronics::core::kPi * f * (double) i / sr);
    return x;
}
}

int main()
{
    std::printf ("felitronics::blend tests\n");
    const double sr = 48000.0;
    using namespace felitronics::blend;

    const auto ir = ramp (2048);
    const MasterParams flat;   // gain 0 dB, filters off → a no-op master

    // Default strips + flat master == a plain Σ that never touches the kernels (independent).
    group ("default params → plain sum");
    {
        std::vector<std::vector<float>> irs { ir, ir, ir };
        std::vector<StripParams> strips (3);                 // all defaults: 0 dB, 0°, 0 ms, filters off
        const auto b = blendIrs (irs, strips, flat, sr);
        std::vector<float> want (ir.size());
        for (std::size_t i = 0; i < ir.size(); ++i) want[i] = 3.0f * ir[i];
        ok (maxAbsDiff (b, want) < 1e-5, "3 identical mics @ default → 3×ir (no kernel ops)");
    }

    // Polarity: two identical mics, one rotated 180° → cancellation.
    group ("phase 180° cancels");
    {
        std::vector<std::vector<float>> irs { ir, ir };
        std::vector<StripParams> strips (2);
        strips[1].phaseDeg = 180.0;
        const auto b = blendIrs (irs, strips, flat, sr);
        double pk = 0.0; for (float v : b) pk = std::max (pk, (double) std::fabs (v));
        ok (pk < 1e-3, "ir + (−ir) ≈ 0");
    }

    // Mute drops a strip; solo overrides mute.
    group ("solo / mute");
    {
        std::vector<std::vector<float>> irs { ir, ir };
        std::vector<StripParams> strips (2);
        strips[1].mute = true;
        ok (maxAbsDiff (blendIrs (irs, strips, flat, sr), ir) < 1e-5, "muted strip drops out → other mic alone");

        std::vector<std::vector<float>> irs3 { ir, ir, ir };
        std::vector<StripParams> s3 (3);
        s3[0].solo = true; s3[1].mute = true;                 // solo on 0; mute on 1 ignored (solo active)
        const auto b = blendIrs (irs3, s3, flat, sr);
        ok (maxAbsDiff (b, ir) < 1e-5, "any solo → only soloed audible (mute irrelevant)");
        ok (channelAudible (s3, 0) && ! channelAudible (s3, 1) && ! channelAudible (s3, 2), "audible mask: solo overrides");
    }

    // Gain: −20 dB → ×0.1 per sample (single mic, flat master).
    group ("gain scaling");
    {
        std::vector<std::vector<float>> irs { ir };
        std::vector<StripParams> strips (1);
        strips[0].gainDb = -20.0;
        const auto b = blendIrs (irs, strips, flat, sr);
        std::vector<float> want (ir.size()); for (std::size_t i = 0; i < ir.size(); ++i) want[i] = 0.1f * ir[i];
        ok (maxAbsDiff (b, want) < 1e-5, "gainDb=−20 → 0.1×ir");
    }

    // Master: −6 dB, filters off → whole blend ×10^(−6/20).
    group ("master gain");
    {
        std::vector<std::vector<float>> irs { ir };
        std::vector<StripParams> strips (1);
        MasterParams m; m.gainDb = -6.0;
        const auto post = blendIrs (irs, strips, m, sr, /*applyMaster*/ true);
        const auto pre  = blendIrs (irs, strips, m, sr, /*applyMaster*/ false);
        const double gm = std::pow (10.0, -6.0 / 20.0);
        std::vector<float> want (ir.size()); for (std::size_t i = 0; i < ir.size(); ++i) want[i] = (float) (gm * ir[i]);
        ok (maxAbsDiff (post, want) < 1e-5, "master −6 dB scales the whole blend");
        ok (maxAbsDiff (pre, ir) < 1e-5, "applyMaster=false → pre-master (unscaled)");
    }

    // "Filter parked at its extreme ≡ off" — hpf.on but hz at kHpfLo → no filtering.
    group ("parked filter ≡ off");
    {
        std::vector<std::vector<float>> irs { ir };
        std::vector<StripParams> parked (1), off (1);
        parked[0].hpf = { true, kHpfLo, 24 };                 // enabled but parked at the low extreme
        const auto bParked = blendIrs (irs, parked, flat, sr);
        const auto bOff    = blendIrs (irs, off,    flat, sr);
        ok (maxAbsDiff (bParked, bOff) < 1e-9, "hpf parked at kHpfLo bypasses (== filter off)");
        ok (hpActiveSlope (parked[0].hpf) == 0, "hpActiveSlope=0 when parked");
    }

    // Degenerate sets return empty.
    group ("degenerate → empty");
    {
        std::vector<std::vector<float>> none;
        std::vector<StripParams> s2 (2);
        ok (blendIrs (none, {}, flat, sr).empty(), "empty irs → {}");
        std::vector<std::vector<float>> one { ir };
        ok (blendIrs (one, s2, flat, sr).empty(), "size mismatch (1 ir, 2 strips) → {}");
    }

    // biquad HPF blocks DC; LPF passes it.
    group ("biquad DC response");
    {
        std::vector<float> hp (400, 1.0f); biquadInplace (hp, sr, 1000.0, true);
        ok (std::fabs (hp[399]) < 1e-3f, "HPF settles DC → 0");
        std::vector<float> lp (400, 1.0f); biquadInplace (lp, sr, 1000.0, false);
        ok (std::fabs (lp[399] - 1.0f) < 1e-2f, "LPF passes DC → ~1");
    }

    // applyBlendSlope is TRUE Butterworth: −3 dB (×0.708) at fc for EVERY (even-order) slope.
    group ("applyBlendSlope −3 dB at fc (all slopes)");
    {
        const double fc = 1000.0;
        for (int slope : { 12, 24, 48, 96 })
        {
            auto v = sine (fc, sr, 8192);
            const double in = rms (v, 4096, 8192);            // steady-state input RMS
            applyBlendSlope (v, sr, fc, /*hp*/ true, slope);
            const double out = rms (v, 4096, 8192);
            approx (out / in, 0.7079, 0.03, "HPF@fc ratio ≈ −3 dB, slope " + std::to_string (slope));
        }
    }

    // applyBlendSlope: an odd order (18 dB/oct → order 3) can't be factored into biquads — no-op guard.
    group ("applyBlendSlope odd-order guard");
    {
        const double fc = 1000.0;
        auto v = sine (fc, sr, 2048);
        const auto orig = v;
        applyBlendSlope (v, sr, fc, /*hp*/ true, 18);
        ok (maxAbsDiff (v, orig) < 1e-12, "slopeDbOct=18 (odd order 3) → no-op");
    }

    // rotatePhase: θ=180° is an exact polarity flip (sin180=0 kills the Hilbert term); θ=0 is a no-op.
    group ("rotatePhase 180°/0°");
    {
        auto x = sine (500.0, sr, 512);
        const auto orig = x;
        rotatePhase (x, 180.0);
        bool flip = true; for (std::size_t i = 0; i < x.size(); ++i) if (std::fabs (x[i] + orig[i]) > 1e-4f) { flip = false; break; }
        ok (flip, "θ=180° == polarity flip (y == −x)");
        auto y = orig; rotatePhase (y, 0.0);
        bool noop = true; for (std::size_t i = 0; i < y.size(); ++i) if (std::fabs (y[i] - orig[i]) > 1e-6f) { noop = false; break; }
        ok (noop, "θ=0 == no-op");
    }

    // shiftFrac by an integer moves a delta by exactly that many samples.
    group ("shiftFrac integer delay");
    {
        std::vector<float> d (64, 0.0f); d[10] = 1.0f;
        const auto y = shiftFrac (d, 5.0);
        ok (std::fabs (y[15] - 1.0f) < 1e-4f && std::fabs (y[10]) < 1e-4f, "δ@10 shifted +5 → δ@15");
    }

    return felitronics::test::report();
}
