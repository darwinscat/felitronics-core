// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// EBU Tech 3341-2023 §2.6 / Table 1, true-peak test signals 15-19 — the FIRST external true-peak criterion
// in this repository. The core holds two independent true-peak paths (felitronics::analysis::TruePeakMeter
// and, in tools/, fcore::Probe), and they must both answer to the same published envelope, so the signals
// live here rather than in either suite: two copies of a specification drift, and a drifted fixture is the
// exact failure this finding already suffered three times.
//
// WHAT THIS IS AND IS NOT. Tests 15-19 are fully specified by Table 1 — frequency, amplitude, phase, taper —
// so they can be synthesized exactly, and that is what happens here. Tests 20-23 CANNOT: their definition
// contains an unspecified "lowpass (anti-aliasing) filter" applied before the 4*fs -> fs decimation, so
// synthesizing them would bake OUR filter choice into the test SIGNAL. They are therefore verified out of
// tree, against the official EBU Loudness Test Set, and are not part of ctest. (Measured 2026-09-04: both
// implementations pass all of 15-23 on the official files.)
//
// WHY SYNTHESIS IS TRUSTWORTHY HERE — measured, not assumed. Run against the official WAVs the same two
// meters read within 0.0005 dB of what they read on these synthesized signals: 400x below the +0.2/-0.4 dB
// tolerance. The 87 MB official set is third-party material that cannot enter a public AGPL repository and
// is a golden-file corpus besides; this fixture reproduces it to four decimal places instead.
//
// THE AMPLITUDE AUTHORITY IS NEITHER THE PRINTED FFS NOR THE FILES — IT IS THE dB TARGET. Three candidate
// values exist and they differ, so the choice has to be named rather than implied:
//   * what Table 1 PRINTS: 0.50 and 1.41 FFS — rounded, and 0.0206 / 0.0259 dB below their own targets;
//   * what the official WAVs CARRY, measured from their integer sample peaks: 0.5011720 and 1.414119;
//   * what makes the targets EXACT: 10^(-6/20) = 0.5011872336 and sqrt(2) = 1.4142135624 — used here.
// The third is the authority because it is derivable from the same table row as the target, so stimulus and
// oracle cannot drift apart; the printed FFS values would insert an undocumented 0.02-0.03 dB guard band
// between them, and a guard band belongs in an assertion where it can be read. The official files sit
// 0.00026 / 0.00058 dB below these constants — unexplained, 700x below tolerance, recorded not modelled —
// which is why readings here reproduce the official ones to 0.0005 dB.
//
// Note that this makes case 19's ANALYTIC peak +3.0103 dBTP while its Table target is +3.0: sqrt(2) is not
// 10^(3/20). The two assertions are deliberately separate for exactly this reason.
//
// TWO TRAPS, BOTH ALREADY PAID FOR:
//   * The taper is not decoration. A tone that stops at full scale gives a real edge-recovery overshoot,
//     which an honest meter finds — and then what is being measured is the edge of the drain, not the
//     filter's pass-band. Without a fade this fixture read +0.92 dB at 1 kHz. Table 1 specifies 10 ms
//     fade-in and fade-out for exactly this reason. The shape is not specified; raised-cosine is used here,
//     and it cannot matter, because the peak lives in the flat region between the fades.
//   * The phase is part of the signal. Table 1's phases are chosen so the sample grid MISSES the peak; at
//     phase 0 (case 15) it hits, and every meter then reads the truth regardless of interpolation quality.
//     Case 15 is the control, not the test.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace felitronics::test::ebu3341
{

struct TruePeakCase
{
    int         number;        // Table 1 test number
    const char* name;
    double      divisor;       // frequency = fs / divisor
    double      amplitude;     // linear; the EXACT value that makes targetDb analytic (see the note)
    double      phaseDeg;
    double      targetDb;      // Table 1's expected max true-peak level, tolerance +0.2 / -0.4 dB
};

inline constexpr double kTolAboveDb = 0.2;    // EBU Tech 3341-2023 §2.6: the envelope is asymmetric
inline constexpr double kTolBelowDb = 0.4;

// The regression gate, and the reason it is DERIVED rather than picked. The dominant error in tests 16-19 is
// not the prototype's pass-band at all — it is grid geometry. Both prototypes are symmetric with an even tap
// count (48 and 128), so the polyphase reconstruction is evaluated at the odd eighths of the input period,
// and Table 1's phases (45, 60, 67.5 degrees) place the crest exactly half an oversampled step away: the
// worst case the geometry allows. Half a step of an L-times grid costs
//
//     20 * log10 (cos (pi * f / (L * fs)))
//
// which for L = 4 is -0.16852 / -0.07463 / -0.04193 dB at fs/4, fs/6, fs/8 — against measured readings of
// -0.16875 / -0.07397 / -0.04206. The residual is at most 0.0007 dB. In other words the whole of the visible
// error is geometry and the pass-band contributes nothing measurable at these three frequencies.
//
// So the assertion is: the reading must sit within kBudgetSlackDb of that bound on the low side and of the
// analytic truth on the high side. It is not a golden pin — it is a function of the oversampling factor,
// which both meters expose — and it therefore catches a 0.02 dB pass-band droop where a flat +-0.25 dB
// budget would have needed 0.25. A redesign is free to move: an odd-length prototype puts the grid ON the
// crest and reads ~0.000, comfortably inside; a drop to 2x reads -0.688 at fs/4 and is caught by the EBU
// envelope above, not here. What the two-sided form does forbid is adopting a prototype whose pass-band
// ripple over-reads by more than 0.02 dB — the literal BS.1770 Annex 2 coefficient table is such a
// prototype. That is a deliberate change in the meter's declared accuracy and must be re-decided here
// rather than pass silently.
//
// (ITU-R BS.1770 Annex 2 tabulates a more conservative form of the same under-read, assuming a full step
// rather than the half step that a nearest-grid-point search actually costs; its 4x figure at f_norm 0.5 is
// 0.688 dB where the geometry above gives 0.169. We assert what we can derive and have measured.)
inline constexpr double kBudgetSlackDb = 0.02;

inline constexpr TruePeakCase kTruePeakCases[] {
    { 15, "Fs/4, 0.50 FFS, phase 0",      4.0, 0.5011872336272722,  0.0, -6.0 },
    { 16, "Fs/4, 0.50 FFS, phase 45",     4.0, 0.5011872336272722, 45.0, -6.0 },
    { 17, "Fs/6, 0.50 FFS, phase 60",     6.0, 0.5011872336272722, 60.0, -6.0 },
    { 18, "Fs/8, 0.50 FFS, phase 67.5",   8.0, 0.5011872336272722, 67.5, -6.0 },
    { 19, "Fs/4, 1.41 FFS, phase 45",     4.0, 1.4142135623730951, 45.0, +3.0 },
};

inline constexpr double kPi = 3.14159265358979323846;

// The analytic true peak: these are continuous sines, so the maximum of the reconstruction is the amplitude
// itself. Independently confirmed by ideal band-limited interpolation (32x and 64x spectral zero-padding) of
// the official WAVs, which lands within 0.0003 dB of this for all five.
inline double oracleDb (const TruePeakCase& c) { return 20.0 * std::log10 (c.amplitude); }

// The exact maximum the SAMPLE GRID can reach, in closed form: the phase advances by 2*pi/divisor per sample,
// so the grid repeats every `divisor` samples and the maximum is a scan of one cycle. Rate-independent by
// construction, which is itself worth asserting — a meter whose sample peak moved with the rate would be
// reading something other than the samples it was handed.
inline double sampleGridPeak (const TruePeakCase& c)
{
    const int    d   = (int) std::lround (c.divisor);
    const double phi = c.phaseDeg * kPi / 180.0;
    double mx = 0.0;
    for (int n = 0; n < d; ++n)
        mx = std::max (mx, std::fabs (std::sin (2.0 * kPi * (double) n / (double) d + phi)));
    return c.amplitude * mx;
}

inline double sampleGridPeakDb (const TruePeakCase& c) { return 20.0 * std::log10 (sampleGridPeak (c)); }

// The lower edge of the accuracy budget: how far below the analytic truth an L-times grid may legitimately
// land on this case, before any filter quality is involved. Takes the factor from the meter under test, so
// the assertion follows a redesign instead of pinning one.
inline double gridBoundDb (const TruePeakCase& c, int oversampleFactor)
{
    return 20.0 * std::log10 (std::cos (kPi / (c.divisor * (double) oversampleFactor)));
}

// The 10 ms raised-cosine taper, as its own function so it can be tested as itself. Indexed from the END on
// the way out, so the two ramps are exact mirrors: `cos(pi * (n - i) / fade)` looks equivalent and is not —
// it leaves the last sample at 0.5 - 0.5*cos(pi/fade) ~ 1.1e-5 of full amplitude instead of zero, with the
// whole ramp shifted a sample. Both forms span the same `fade` intervals, so nothing is shortened; only one
// of them actually ends in silence. That distinction changes no maximum here — the peak lives in the flat
// region — which is precisely why a fixture bug of this shape survives every assertion about maxima, and why
// the symmetry of THIS function is asserted directly rather than inferred from the signal (the carrier has a
// different phase at mirrored positions, so comparing samples would measure the sine, not the window).
inline double fadeWindow (long long i, long long n, long long fade)
{
    if (fade <= 0 || n <= 0) return 1.0;
    if (i < fade)         return 0.5 - 0.5 * std::cos (kPi * (double) i / (double) fade);
    if (i >= n - fade)    return 0.5 - 0.5 * std::cos (kPi * (double) (n - 1 - i) / (double) fade);
    return 1.0;
}

// ONE channel, `seconds` long, with the specified 10 ms raised-cosine fades at both ends. Table 1's signals
// are stereo in phase, which the callers do by handing this same buffer as both channels.
//
// Frequency is written as fs/divisor the way Table 1 phrases it; at 48 and 44.1 kHz that yields the same
// normalized signal, so running both rates tests that prepare() handles them, not that the signal differs.
//
// The fade-out is indexed FROM THE END so that it mirrors the fade-in exactly. Writing it as
// `cos(pi * (n - i) / fade)` looks equivalent and is not: at i = n-1 that leaves the last sample at
// 0.5 - 0.5*cos(pi/fade) ~ 1.1e-5 of full amplitude instead of zero, and the whole fade-out shifted one
// sample. It changes no maximum here — the peak lives in the flat region — which is exactly why a fixture
// bug of this shape survives every assertion about maxima. `endsAreSilent()` below is the check that does
// catch it.
inline void synthesize (const TruePeakCase& c, double sampleRate, double seconds, std::vector<float>& out)
{
    const long long n    = (long long) std::llround (sampleRate * seconds);
    const long long fade = (long long) std::llround (sampleRate * 0.010);      // Table 1: 10 ms
    const double    hz   = sampleRate / c.divisor;
    const double    phi  = c.phaseDeg * kPi / 180.0;

    out.assign ((std::size_t) n, 0.0f);
    for (long long i = 0; i < n; ++i)
        out[(std::size_t) i] = (float) (fadeWindow (i, n, fade) * c.amplitude
                                        * std::sin (2.0 * kPi * hz * (double) i / sampleRate + phi));
}

// The fixture self-check no assertion about a maximum can make: both ends must actually reach silence, and
// the fade must be symmetric. Everything else in this suite looks only at maxima, so a malformed envelope
// that leaves the file ending at -100 dBFS instead of zero would pass unnoticed — and a file that does not
// end in silence is precisely the condition that made the meter measure its own edge transient and produced
// the +0.92 dB absurdity this finding started from.
inline bool endsAreSilent (const std::vector<float>& x)
{
    return ! x.empty() && x.front() == 0.0f && x.back() == 0.0f;
}

// Worst mismatch between the two ramps of the taper itself. Zero when they mirror exactly.
inline double fadeAsymmetry (long long n, long long fade)
{
    double worst = 0.0;
    for (long long k = 0; k < fade; ++k)
        worst = std::max (worst, std::fabs (fadeWindow (k, n, fade) - fadeWindow (n - 1 - k, n, fade)));
    return worst;
}

} // namespace felitronics::test::ebu3341
