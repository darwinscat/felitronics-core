// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::measurement — the IR-measurement pipeline (sweep + Farina
// deconvolution + IR post + capture gate + multi-mic align). Oracle/property tests (no golden files):
// the reference is recomputed independently (sweep⊛inverse≈δ; a known IR round-trips through deconv).
// Adversarial extremes + a scipy/numpy cross-NULL live in the wider ctest campaign.

#include <felitronics_test.h>
#include <felitronics/measurement/Sweep.h>
#include <felitronics/measurement/Deconvolve.h>
#include <felitronics/measurement/IrPost.h>
#include <felitronics/measurement/CaptureGate.h>
#include <felitronics/measurement/MicSetAlign.h>

#include <cmath>
#include <cstddef>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::approx;
using felitronics::test::group;

namespace
{
double peakAbs (const std::vector<double>& x, std::size_t& idxOut)
{
    double pk = 0.0; idxOut = 0;
    for (std::size_t i = 0; i < x.size(); ++i) { const double a = std::fabs (x[i]); if (a > pk) { pk = a; idxOut = i; } }
    return pk;
}
}

int main()
{
    std::printf ("felitronics::measurement tests\n");

    measurement::SweepSpec spec;
    spec.f1 = 100.0; spec.f2 = 8000.0; spec.durationSeconds = 0.3; spec.sampleRate = 48000.0;
    spec.tailSeconds = 0.1; spec.fadeSeconds = 0.01;
    const auto sw = measurement::makeSweep (spec);

    // --- Oracle 1: sweep_proper ⊛ inverse ≈ δ (unit peak at sweepLen-1, low sidelobes) ---
    group ("sweep ⊛ inverse ≈ δ");
    {
        std::vector<double> proper (sw.signal.begin(), sw.signal.begin() + (std::ptrdiff_t) sw.sweepLen);
        const auto d = measurement::convolve (proper, sw.inverse);
        std::size_t pk = 0; const double peak = peakAbs (d, pk);
        approx (peak, 1.0, 0.02, "delta unit peak");
        ok (pk == sw.sweepLen - 1, "delta lands at sweepLen-1");
        // sidelobe floor outside a small guard window
        double side = 0.0;
        for (std::size_t i = 0; i < d.size(); ++i)
            if (i + 64 < pk || i > pk + 64) side = std::max (side, std::fabs (d[i]));
        ok (side < 0.05 * peak, "sidelobes < -26 dB outside guard");
    }

    // --- Oracle 2: known-answer deconvolution + injected round-trip latency (deconv MATH; the linear
    // IR delta must land at sweepLen-1+tau and reproduce the known IR). Onset WALK-BACK is tested
    // separately on a realistic attack (Oracle 3) — a hard synthetic delta has non-physical pre-ring. ---
    group ("known-answer deconvolution (latency absorbed)");
    {
        const std::vector<double> knownIr { 1.0, 0.0, 0.5, 0.0, -0.3, 0.1 };   // peak tap at offset 0
        auto rec = measurement::convolve (
            std::vector<double> (sw.signal.begin(), sw.signal.end()), knownIr);
        const std::size_t tau = 137;                                   // fake round-trip latency
        rec.insert (rec.begin(), tau, 0.0);
        const auto dr = measurement::deconvolve (rec, sw, -40.0, /*preRoll*/ 0);
        const std::size_t lag = sw.sweepLen - 1 + tau;                 // where the linear IR must land
        std::size_t argmax = sw.sweepLen - 1; double pk = 0.0;
        for (std::size_t i = sw.sweepLen - 1; i < dr.full.size(); ++i)
        { const double a = std::fabs (dr.full[i]); if (a > pk) { pk = a; argmax = i; } }
        ok (argmax == lag, "linear-IR delta lands at sweepLen-1+tau (latency absorbed)");
        // The recovered transfer function matches the known system's IN-BAND (the sweep band-limits the
        // delta, so knownIr⊛delta smears in time — compare magnitude over [2*f1, f2/2], normalized).
        const std::size_t W = 1024, nfft = 2048;
        std::vector<double> ext (nfft, 0.0);
        for (long k = -(long) W; k < (long) W; ++k)
        { const long idx = (long) lag + k; if (idx >= 0 && idx < (long) dr.full.size()) ext[(std::size_t) (k + (long) W)] = dr.full[(std::size_t) idx]; }
        const auto extMag = measurement::magSpectrum (ext, nfft);
        const auto knMag  = measurement::magSpectrum (knownIr, nfft);
        const double binHz = spec.sampleRate / (double) nfft;
        const std::size_t b0 = (std::size_t) (2.0 * spec.f1 / binHz), b1 = (std::size_t) (0.5 * spec.f2 / binHz);
        auto inbandMax = [&] (const std::vector<double>& m) { double x = 0.0; for (std::size_t b = b0; b <= b1; ++b) x = std::max (x, m[b]); return x > 0.0 ? x : 1.0; };
        const double en = inbandMax (extMag), kn = inbandMax (knMag);
        double worst = 0.0;
        for (std::size_t b = b0; b <= b1; ++b) worst = std::max (worst, std::fabs (extMag[b] / en - knMag[b] / kn));
        ok (worst < 0.08, "recovered transfer function matches known IR in-band");
    }

    // --- Oracle 3: onset detection robust to a pre-onset floor ---
    group ("onset detection");
    {
        std::vector<double> ir (2000, 0.0);
        for (std::size_t i = 0; i < 500; ++i) ir[i] = 0.001 * std::sin (0.3 * (double) i); // pre-floor
        ir[900] = 1.0; ir[901] = -0.7; ir[902] = 0.4;                  // the real attack
        const auto on = measurement::detectOnset (ir, -40.0, 8, 0);
        ok (on.peakIndex == 900, "peak found");
        ok (on.onset >= 892 && on.onset <= 900, "onset near attack, not the pre-floor");
    }

    // --- Oracle 4: capture gate (clean / NaN / noise / clip) ---
    group ("capture gate");
    {
        std::vector<double> clean (sw.signal.begin(), sw.signal.end());   // record = the sweep itself
        const auto gc = measurement::gateRecording (clean, sw);
        ok (gc.ok && gc.reason == measurement::GateReject::None, "clean sweep passes");
        ok (gc.sweepPresent, "sweep detected");

        auto nan = clean; nan[1234] = std::nan ("");
        const auto gn = measurement::gateRecording (nan, sw);
        ok (! gn.ok && gn.reason == measurement::GateReject::NonFinite, "NaN → NonFinite (honest reason)");

        std::vector<double> noise (clean.size(), 0.0);
        for (std::size_t i = 0; i < noise.size(); ++i) noise[i] = 0.01 * std::sin (12.9898 * (double) i); // deterministic 'noise'
        const auto gnz = measurement::gateRecording (noise, sw);
        ok (! gnz.ok && gnz.reason == measurement::GateReject::SweepNotDetected, "no sweep → SweepNotDetected");

        auto clip = clean; for (int i = 0; i < 8; ++i) clip[(std::size_t) (2000 + i)] = 1.0;
        const auto gcl = measurement::gateRecording (clip, sw);
        ok (gcl.clipped && gcl.reason == measurement::GateReject::Clipped, "flat-top run → Clipped");
    }

    // --- Oracle 5: multi-mic common-onset preserves the inter-mic delay ---
    group ("multi-mic alignment preserves inter-mic delay");
    {
        const std::vector<double> ir { 1.0, 0.0, 0.5, -0.2 };
        const std::size_t deltaMic = 17;
        auto mk = [&] (std::size_t extraLatency)
        {
            auto rec = measurement::convolve (std::vector<double> (sw.signal.begin(), sw.signal.end()), ir);
            rec.insert (rec.begin(), extraLatency, 0.0);
            return measurement::deconvolve (rec, sw, -40.0, 0);
        };
        std::vector<measurement::DeconvResult> mics { mk (200), mk (200 + deltaMic) };
        const auto set = measurement::alignToCommonOnset (mics, 64);
        ok (set.delaySamples.size() == 2, "two mics aligned");
        ok (set.delaySamples[0] == 0, "earliest mic at delay 0");
        ok (set.delaySamples[1] == deltaMic, "second mic keeps the 17-sample inter-mic delay");
    }

    return felitronics::test::report();
}
