// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P57 — core::DeliveryResampler, THE BAR: stopband >= 140 dB and passband ripple <= 0.001 dB to
// 20 kHz, measured THROUGH the converter on all 30 directed delivery pairs.
//
// 🔴 WHY A SPECTRUM, AND WHY THIS ONE. The round-trip null in the sibling suite bottoms out at the
// float32 quantisation of its own API boundaries and cannot see a -140 dB stopband at all. Only a
// spectrum separates the through component from the spurs — and the first version of this spectrum
// was blind in three measured ways, each of which is why the instrument now looks the way it does:
//
//  1. IMAGES DO NOT LAND ON BIN CENTRES. The carrier was placed on an exact output bin so a
//     rectangular window leaked nothing, and the file claimed every alias landed on a bin too. Aliases
//     at m*fsOut +/- f do; IMAGES at m*fsIn +/- f do not, whenever fsIn*N/fsOut is not an integer
//     (44.1 <-> 48: 30105.6 bins). A rectangular FFT under-reads such a line by up to 3.9 dB, and an
//     injected -138 dB image read -141.3 and PASSED. So spurs are now read from a Kaiser-windowed
//     (beta 30), 4x zero-padded spectrum — every line within a fraction of a dB wherever it falls — and
//     the through gain, which must be read to 1e-6 dB, still comes from the rectangular bin of the
//     bin-centred carrier, where it is exact.
//  2. FIVE STOPBAND PROBES MISSED THE WORST LOBE. On 88.2 -> 44.1 they sat near a null (-181 dB) and
//     beside a lobe (-153), while a 0.1 Hz scan of the coefficients then shipped found -145.7 dB at
//     22130.3 Hz. A Kaiser stopband's worst lobe sits within a few lobe-widths of the edge, so every
//     decimating stage edge is now swept densely, at a small fraction of its own transition width.
//  3. IT NEVER ASSERTED THE LEVEL. Ripple was peak-to-peak among the probes; a converter at half gain,
//     inverted, or silent on every two-stage route passed all 127 checks — and its spurs READ BETTER,
//     being measured re 1.0 instead of re the carrier. Absolute gain is now asserted on every probe, and
//     a non-finite sample anywhere in an analysed segment is a failure rather than `max(1e-300, NaN)`.
//
//  4. A REAL LINE AT DC OR AT NYQUIST IS NOT SPLIT BETWEEN TWO BINS, so the single-sided `2|X|/sum(w)`
//     that is exact everywhere else reads it 6.02 dB HIGH. The first run of the fixed instrument
//     failed six pairs at exactly 0 Hz and exactly fsOut/2 — every one of them this factor, measured
//     on synthetic lines (-140.00 reads -133.98). Those two lines are now read by exact projections
//     onto 1 and (-1)^n, the scan skips a guard around them, and a tone placed EXACTLY on the output
//     Nyquist is probed in sine and cosine phase, because a sampled line there keeps only its cosine
//     part: the delivered level is sqrt(p^2 + q^2) of the two readings, not either one.
//
// Negative controls pin the instrument itself before any converter is measured: a 0 dB carrier ALONE
// must read no spur above -190 dB (the window's own floor); that carrier PLUS a -138 dB line placed
// 0.4626 bins off-centre must read -138 to within 0.2 dB — the gate must be able to FAIL; and -140 dB
// lines at exactly DC and exactly Nyquist must read -140, not -134.
//
// WHAT IS ASSERTED, AND WHERE THE SPEC IS DELIBERATELY SILENT:
//   * f <= passband edge    -> through gain within 0.001 dB of unity, and every spur at or under the bar
//                              relative to the carrier (IMAGING);
//   * f >= min(fsIn,fsOut)/2 -> every line at or under the bar re the input (ALIASING), including a
//                              full-scale tone at EXACTLY the input Nyquist;
//   * in between             -> nothing: that is the transition band, and a fold from it lands back in
//                              the same band, never in the protected passband.

#include <felitronics/core/DeliveryResampler.h>
#include <felitronics/core/OfflineFft.h>
#include <felitronics_test.h>

#include <cmath>
#include <complex>
#include <vector>

using felitronics::core::DeliveryResampler;
using felitronics::test::ok;

namespace
{
    constexpr double kPi      = 3.14159265358979323846;
    constexpr int    kSeg     = 16384;     // analysed output samples
    constexpr int    kPad     = 4;         // zero-padding factor for the spur spectrum
    constexpr double kBeta    = 30.0;      // analysis window
    constexpr int    kGuard   = 12;        // main-lobe exclusion (beta 30 lobe: 9.6 bins), in segment bins
    constexpr double kBarDb    = -140.0;
    constexpr double kRippleDb = 0.001;
    constexpr double kRates[] = { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

    // 🔴 THE CLAIM IS GATED, NOT ONLY THE REQUIREMENT. See the header block of DeliveryResampler.h: a
    // gate at -140 alone stayed green when a mutation stand dropped the design target to 140.
    // Delivered -146.58 dB. The gate sits 0.38 dB above it and below every regression measured on the way
    // here: no edge margin -142.01, A = 143 -142.32, a 0.03 margin -141.44, per-phase normalisation
    // -143.74, A = 145 -145.14 — and a Kaiser beta one dB low, -145.62, which a 1.1 dB gate let through
    // on a mutation stand. The margin is not fragile, and that is measured rather than argued: this suite
    // printed the same worst line, ripple, gain and instrument floor on five rows — arm64 macOS, x86-64
    // macOS (Apple libm), x86-64 Linux (gcc 14.2, glibc), x86-64 Windows (MSVC 19.44, UCRT) and wasm32
    // (emsdk 6.0.9, node). The coefficients differ between them only in libm's last places.
    constexpr double kClaimDb       = -146.2;
    constexpr double kClaimRippleDb = 0.00001;

    double besselI0 (double x)
    {
        double sum = 1.0, term = 1.0;
        const double y = x * x * 0.25;
        for (int k = 1; k < 400; ++k) { term *= y / ((double) k * (double) k); sum += term; if (term < 1e-18 * sum) break; }
        return sum;
    }

    struct Spectrum { double throughDb; double spurDb; double spurHz; double dcAmp; double nyqAmp; bool finite; };

    double db (double amp) { return 20.0 * std::log10 (std::max (1e-300, amp)); }

    // `y` is kSeg samples at fsOut. `carrierBin` < 0 means there is no carrier to exclude. With
    // `exactEnds`, the DC and Nyquist lines are read by exact projection and the scan skips their guard.
    Spectrum analyse (const std::vector<double>& y, double fsOut, int carrierBin, bool exactEnds)
    {
        Spectrum sp { -999.0, -999.0, 0.0, 0.0, 0.0, true };
        for (double v : y) if (! std::isfinite (v)) { sp.finite = false; return sp; }

        if (carrierBin >= 0)
        {
            std::vector<std::complex<double>> rect ((std::size_t) kSeg);
            for (int i = 0; i < kSeg; ++i) rect[(std::size_t) i] = { y[(std::size_t) i], 0.0 };
            felitronics::core::offline::fftInplace (rect, +1);
            sp.throughDb = 20.0 * std::log10 (std::max (1e-300, 2.0 * std::abs (rect[(std::size_t) carrierBin]) / (double) kSeg));
        }

        const int nfft = kSeg * kPad;
        std::vector<std::complex<double>> win ((std::size_t) nfft, { 0.0, 0.0 });
        const double i0b = besselI0 (kBeta);
        double wsum = 0.0;
        for (int i = 0; i < kSeg; ++i)
        {
            const double u = 2.0 * (double) i / (double) (kSeg - 1) - 1.0;
            const double w = besselI0 (kBeta * std::sqrt (std::max (0.0, 1.0 - u * u))) / i0b;
            win[(std::size_t) i] = { w * y[(std::size_t) i], 0.0 };
            wsum += w;
        }
        double dc = 0.0, nyq = 0.0;
        for (int i = 0; i < kSeg; ++i)
        {
            const double wy = win[(std::size_t) i].real();
            dc += wy;
            nyq += (i & 1) ? -wy : wy;
        }
        sp.dcAmp  = std::fabs (dc) / wsum;
        sp.nyqAmp = std::fabs (nyq) / wsum;

        felitronics::core::offline::fftInplace (win, +1);
        for (int i = 0; i <= nfft / 2; ++i)
        {
            if (carrierBin >= 0 && std::abs (i - carrierBin * kPad) <= kGuard * kPad) continue;
            if (exactEnds && (i <= kGuard * kPad || i >= nfft / 2 - kGuard * kPad)) continue;
            const double d = db (2.0 * std::abs (win[(std::size_t) i]) / wsum);
            if (d > sp.spurDb) { sp.spurDb = d; sp.spurHz = (double) i * fsOut / (double) nfft; }
        }
        if (exactEnds)
        {
            if (db (sp.dcAmp) > sp.spurDb)  { sp.spurDb = db (sp.dcAmp);  sp.spurHz = 0.0; }
            if (db (sp.nyqAmp) > sp.spurDb) { sp.spurDb = db (sp.nyqAmp); sp.spurHz = 0.5 * fsOut; }
        }
        return sp;
    }

    // One tone (or the literal input Nyquist, `fHz < 0`) through a prepared converter; the steady-state
    // segment is analysed. Returns false if the converter refused a call.
    bool probe (DeliveryResampler& r, double fsIn, double fsOut, double fHz, int carrierBin, Spectrum& out,
                bool exactEnds = false, double phase = 0.0)
    {
        r.reset();
        const int skipOut = 1024 + 4 * (int) std::ceil (r.latencyOutputSamples());
        const int wantOut = skipOut + kSeg;
        const int nIn = (int) std::ceil ((double) wantOut * fsIn / fsOut) + 2048;
        const int block = 4096;
        std::vector<float> in ((std::size_t) block);
        std::vector<float> buf ((std::size_t) r.maxOutputFor (block));
        std::vector<double> y;
        y.reserve ((std::size_t) wantOut + 8192);
        const double fade = 512.0;
        for (int i = 0; i < nIn && (int) y.size() < wantOut; i += block)
        {
            const int n = std::min (block, nIn - i);
            for (int j = 0; j < n; ++j)
            {
                const double t = (double) (i + j);
                const double w = (t < fade) ? (0.5 - 0.5 * std::cos (kPi * t / fade)) : 1.0;
                const double s = (fHz < 0.0) ? (((i + j) & 1) ? -1.0 : 1.0)
                                             : std::sin (2.0 * kPi * fHz * t / fsIn + phase);
                in[(std::size_t) j] = (float) (w * s);
            }
            const float* ip[1] = { in.data() };
            float* op[1] = { buf.data() };
            int got = 0;
            if (! r.process (ip, 1, n, op, (int) buf.size(), got)) return false;
            for (int k = 0; k < got; ++k) y.push_back ((double) buf[(std::size_t) k]);
        }
        if ((int) y.size() < wantOut) return false;
        std::vector<double> seg (y.begin() + skipOut, y.begin() + skipOut + kSeg);
        out = analyse (seg, fsOut, carrierBin, exactEnds);
        return true;
    }

    int binOf (double fHz, double fsOut) { return (int) std::llround (fHz * (double) kSeg / fsOut); }
    double snap (double wantHz, double fsOut)
    {
        const double step = fsOut / (double) kSeg;
        return std::max (step, std::round (wantHz / step) * step);
    }
}

//==============================================================================
static void testInstrument()
{
    std::printf ("  instrument negative controls\n");
    const double fs = 48000.0;
    std::vector<double> y ((std::size_t) kSeg);

    // A 0 dB carrier alone: whatever the window leaks is the floor under every reading below.
    for (int i = 0; i < kSeg; ++i) y[(std::size_t) i] = std::sin (2.0 * kPi * 20.0 * i / (double) kSeg);
    const Spectrum a = analyse (y, fs, 20, false);
    ok (a.finite, "instrument: finite");
    ok (std::fabs (a.throughDb) < 1e-9, "instrument: a bin-centred 0 dB carrier reads 0 dB");
    ok (a.spurDb < -190.0, "instrument: the window's own floor is under -190 dB");

    // The same carrier plus a -138 dB line 0.4626 bins off-centre — the image placement that fooled the
    // rectangular instrument into reading -141.3 and passing. It must read -138 and FAIL the bar.
    const double off = 2224.0 + 0.462585034;
    const double lvl = std::pow (10.0, -138.0 / 20.0);
    for (int i = 0; i < kSeg; ++i)
        y[(std::size_t) i] = std::sin (2.0 * kPi * 20.0 * i / (double) kSeg)
                           + lvl * std::sin (2.0 * kPi * off * i / (double) kSeg);
    const Spectrum b = analyse (y, fs, 20, false);
    ok (std::fabs (b.spurDb - (-138.0)) < 0.2, "instrument: an off-bin -138 dB line reads -138 dB");
    ok (b.spurDb > kBarDb, "instrument: …and would fail the bar, as it must");
    std::printf ("    carrier alone: floor %.1f dB;  off-bin -138 line reads %.3f dB\n", a.spurDb, b.spurDb);

    // A real line at exactly DC, and at exactly Nyquist: read -140, not the -134 a single-sided scan gives.
    const double l140 = std::pow (10.0, -140.0 / 20.0);
    for (int i = 0; i < kSeg; ++i) y[(std::size_t) i] = l140;
    const Spectrum c = analyse (y, fs, -1, true);
    for (int i = 0; i < kSeg; ++i) y[(std::size_t) i] = (i & 1) ? -l140 : l140;
    const Spectrum d = analyse (y, fs, -1, true);
    ok (std::fabs (c.spurDb - (-140.0)) < 0.01 && c.spurHz < 1.0, "instrument: a -140 dB DC line reads -140 dB at 0 Hz");
    ok (std::fabs (d.spurDb - (-140.0)) < 0.01 && d.spurHz > 0.4 * fs, "instrument: a -140 dB Nyquist line reads -140 dB at fs/2");
    std::printf ("    exact ends: DC line reads %.3f dB, Nyquist line reads %.3f dB\n", c.spurDb, d.spurDb);

    // A NaN in the segment is a failure, not a silence.
    y[100] = std::nan ("");
    ok (! analyse (y, fs, 20, false).finite, "instrument: a NaN in the segment is reported, not hidden");
}

int main()
{
    std::printf ("felitronics::core::DeliveryResampler — the bar, measured on all 30 pairs\n");
    testInstrument();

    std::printf ("  %-14s %-8s %11s %10s %10s %9s %7s\n", "pair", "route", "ripple/dB", "gain/dB", "worst/dB", "at Hz", "probes");
    double worstRipple = 0.0, worstGain = 0.0, worstSpur = -999.0;
    char worstRipplePair[24] = "", worstSpurPair[24] = "";
    int totalProbes = 0;

    for (double a : kRates)
        for (double b : kRates)
        {
            if (a == b) continue;
            char nm[24];
            std::snprintf (nm, sizeof (nm), "%g->%g", a / 1000.0, b / 1000.0);

            DeliveryResampler r;
            DeliveryResampler::Params p;
            p.inRate = a; p.outRate = b;
            if (! r.prepare (p, 1, 4096)) { ok (false, std::string ("prepare: ") + nm); continue; }
            const auto& pl = r.currentPlan();
            const double band = pl.stage[0].passbandHz;
            const double stopEdge = 0.5 * std::min (a, b);
            const double df = stopEdge - band;

            double lo = 1e9, hi = -1e9, gainErr = 0.0, spur = -999.0, spurHz = 0.0;
            int probes = 0;
            bool ran = true, finite = true;

            auto passband = [&] (double f)
            {
                f = snap (f, b);
                if (! (f < band)) return;
                Spectrum s;
                if (! probe (r, a, b, f, binOf (f, b), s)) { ran = false; return; }
                ++probes;
                if (! s.finite) { finite = false; return; }
                lo = std::min (lo, s.throughDb); hi = std::max (hi, s.throughDb);
                gainErr = std::max (gainErr, std::fabs (s.throughDb));
                const double rel = s.spurDb - s.throughDb;             // re the carrier, not re 1.0
                if (rel > spur) { spur = rel; spurHz = s.spurHz; }
            };
            auto stopband = [&] (double f)
            {
                f = snap (f, b);
                if (! (f >= stopEdge && f < 0.5 * a)) return;
                Spectrum s;
                if (! probe (r, a, b, f, -1, s)) { ran = false; return; }
                ++probes;
                if (! s.finite) { finite = false; return; }
                if (s.spurDb > spur) { spur = s.spurDb; spurHz = s.spurHz; }
            };
            const double guardHz = (double) (kGuard + 4) * b / (double) kSeg;

            // Passband: flatness and imaging, plus a denser walk just under the band edge, whose images
            // land nearest the interpolation filter's own stopband edge.
            for (double frac : { 0.02, 0.05, 0.2, 0.45, 0.7, 0.88, 0.97 }) passband (frac * band);
            for (int i = 0; i < 6; ++i) passband (band - df * 0.04 * (double) (i + 1));

            // Stopband: only a decimating pair has input content above the output Nyquist. Every
            // decimating stage's own edge, referred to the input, is swept at a small fraction of that
            // stage's transition width — that is where a Kaiser stopband keeps its worst lobe.
            if (a > b)
            {
                for (int s = 0; s < pl.count; ++s)
                {
                    const auto& st = pl.stage[s];
                    if (st.outRate >= st.inRate) continue;
                    // From the stage's NYQUIST, where aliasing begins — not from its design stop edge,
                    // which sits inside it by kStopEdgeMargin and bounds nothing that can fold.
                    const double e = st.nyquistHz, w = st.nyquistHz - st.passbandHz;
                    for (int i = 0; i < 40; ++i) stopband (e + guardHz + w * 0.015 * (double) i);
                }
                for (double mul : { 1.25, 1.6, 1.98 }) stopband (mul * stopEdge);

                // EXACTLY on the output Nyquist — the stopband's first point. Sine and cosine phase, read
                // by exact projection: the delivered line is sqrt(p^2 + q^2) of the two.
                if (stopEdge < 0.5 * a)
                {
                    Spectrum s0, s1;
                    const bool ok0 = probe (r, a, b, stopEdge, -1, s0, true, 0.0);
                    const bool ok1 = probe (r, a, b, stopEdge, -1, s1, true, 0.5 * kPi);
                    if (! ok0 || ! ok1) ran = false;
                    else
                    {
                        probes += 2;
                        if (! s0.finite || ! s1.finite) finite = false;
                        const double edgeDb = db (std::sqrt (s0.nyqAmp * s0.nyqAmp + s1.nyqAmp * s1.nyqAmp));
                        const double other  = std::max (s0.spurDb, s1.spurDb);
                        const double worst  = std::max (edgeDb, other);
                        if (worst > spur) { spur = worst; spurHz = (edgeDb >= other) ? stopEdge : s0.spurHz; }
                    }
                }
            }
            // The literal input Nyquist — a full-scale (+1, -1) alternation, which a sine cannot express.
            // It sits at or above every stage's stopband edge on every pair, so it must be gone. Its alias
            // lands exactly on DC or on the output Nyquist whenever the ratio is dyadic, so it is read with
            // exact ends.
            {
                Spectrum s;
                if (! probe (r, a, b, -1.0, -1, s, true)) ran = false;
                else
                {
                    ++probes;
                    if (! s.finite) finite = false;
                    if (s.spurDb > spur) { spur = s.spurDb; spurHz = s.spurHz; }
                }
            }

            totalProbes += probes;
            ok (ran, std::string ("every probe ran: ") + nm);
            ok (finite, std::string ("every analysed sample is finite: ") + nm);
            const double ripple = hi - lo;
            ok (ripple <= kRippleDb, std::string ("passband ripple <= 0.001 dB: ") + nm);
            ok (gainErr <= kRippleDb, std::string ("absolute gain within 0.001 dB of unity: ") + nm);
            ok (spur <= kBarDb, std::string ("no alias or image above -140 dB: ") + nm);

            if (ripple > worstRipple) { worstRipple = ripple; std::snprintf (worstRipplePair, sizeof (worstRipplePair), "%s", nm); }
            if (spur > worstSpur)     { worstSpur = spur;     std::snprintf (worstSpurPair, sizeof (worstSpurPair), "%s", nm); }
            worstGain = std::max (worstGain, gainErr);

            char route[16];
            std::snprintf (route, sizeof (route), "%d stage%s", pl.count, pl.count == 1 ? "" : "s");
            std::printf ("  %-14s %-8s %11.7f %10.7f %10.2f %9.0f %7d\n", nm, route, ripple, gainErr, spur, spurHz, probes);
        }

    std::printf ("\n  worst passband ripple : %.7f dB  (%s)   bar 0.001\n", worstRipple, worstRipplePair);
    std::printf ("  worst absolute gain   : %.7f dB             bar 0.001\n", worstGain);
    std::printf ("  worst spur            : %.2f dB     (%s)   bar -140.00   claim %.2f\n", worstSpur, worstSpurPair, kClaimDb);
    std::printf ("  probes                : %d\n", totalProbes);

    ok (worstSpur <= kClaimDb, "the matrix delivers the CLAIMED figure, not merely the -140 dB bar");
    ok (worstRipple <= kClaimRippleDb, "the matrix delivers the claimed ~2e-6 dB flatness");
    ok (totalProbes > 1000 && worstSpur > -300.0 && worstRipple > 0.0, "the sweep measured something (not a vacuous pass)");
    return felitronics::test::report();
}
