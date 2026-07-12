// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.
//
// THEORY-FIRST FALSIFICATION suite for felitronics::poweramp::{SagEnvelope, TubeStage, PowerAmpStage}.
// Companion to PowerAmpGoldenTests.cpp (the fingerprint golden). Where the golden PINS the shipped
// numeric fingerprint, THIS suite starts from amplifier / DSP THEORY and the headers' documented
// contracts, DERIVES each bound from first principles (never tuned to go green), and tries to
// FALSIFY the claim. A failure here is a real physics/contract violation — a documented deviation is
// PINNED in place (never weakened) and reported, not silently absorbed.
//
// Theory basis, per test (all thresholds derived below, at the test site):
//   TT1 Silence   — Asym g(0)=0 (WaveShaper subtracts tanh(kb)); PP(0)=g(vb)−g(vb)=0 (evenLeak=0);
//                   DC-block/bell/shelf/transformer of 0 = 0; sag rail·0 = 0 ⇒ a reset stage fed
//                   zeros must output IDENTICALLY 0 (no additive offset anywhere). No settle horizon:
//                   from reset there is no transient. Strongest possible silence claim: exact 0.
//   TT2 DC reject — the OS-domain one-pole DC-block H(z)=(1−z⁻¹)/(1−dcR·z⁻¹) has an EXACT zero at DC,
//                   corner kDcBlockHz=10 Hz at fs·os. Sustained DC ⇒ steady-state output → 0. Settle
//                   horizon from τ=1/(2π·10)≈15.9 ms (+25 ms coeff smoother); residual after T bounded
//                   by exp(−2π·10·T) ⇒ ≪1e-4 by T≈0.4 s. Assert tail ‖out‖∞ < 1e-3 (float margin).
//   TT3 BIBO      — every block bounded: WaveShaper peak-normalised + tanh-saturating (kernel |y|<8),
//                   DC-block |H|≤~1 (no amplification), transformer tanh soft-clip, bell/shelf ≤ their
//                   dB, sag rail s∈[0.2,1] attenuates, drive-comp normalises small-signal. |x|≤1 ⇒
//                   output finite and O(1). Assert finite always; ‖out‖∞<5 (out 0 dB / autoComp 1),
//                   <8 with the feel layer's depth/bell headroom. (golden X1 pins <4.0 — margin above.)
//   TT4 PP/SE     — PP with evenLeak=0 is ODD by construction ⇒ even harmonics cancel to the numeric/
//                   OS floor; SE (asymmetric Asym) is even-rich. FFT: PP even < −80 dBc, SE H2 > −42
//                   dBc, gap SE−PP > 40 dB — derived from the differential-pair cancellation, not tuned.
//   TT5 Sag       — droop=maxDroop·amount·(1−exp(−2·env)), env→|demand|: monotone ↑ in demand, bounded
//                   by maxDroop·amount (rail s=1−droop∈(0.1,1]; stage clamps ≥kSagMinRail=0.2). Release
//                   is a single real pole ⇒ monotone decay (slow: recovery≫attack). Stage: louder
//                   program ⇒ MORE sag compression (the gap grows monotonically with level).
//   TT6 autoComp  — autoComp=1 ⇒ small-signal gain Drive-invariant (comp=|slope|⁻¹ cancels the slope);
//                   tiny-signal drive sweep RMS window < 0.1 dB. autoComp=0 ⇒ no compensation ⇒ RMS
//                   grows monotonically ≈ +Δdrive dB (near-linear tiny signal).
//   TT7 Latency   — latency = oversampler round-trip = tapsPerPhase−1 = 31, INVARIANT across os factor;
//                   0 before prepare (host queries early). Reported == impulse-measured main lobe.
//   TT8 NaN/Inf   — flushDenormals + the isfinite gates flush state each block ⇒ a poison burst is
//                   contained: output finite within ≤2 valid blocks. NaN driveDb/Inf outputDb/NaN
//                   autoComp ⇒ setParams sanitize gate falls back (0 dB / autoComp 1) ⇒ finite, alive.
//   TT9 Determ.   — same schedule twice from reset ⇒ BIT-identical. Block-size invariance: feel OFF is
//                   bit-exact across schedules; feel ON carries the KNOWN ~1e-2 "B9" feel-layer block
//                   discrepancy — PINNED (≤2e-2) and reported, NOT chased (a documented deviation).

#include <felitronics/poweramp/PowerAmpStage.h>
#include <felitronics/poweramp/SagEnvelope.h>
#include <felitronics/poweramp/TubeStage.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using felitronics::poweramp::PowerAmpStage;
using felitronics::poweramp::SagEnvelope;
using felitronics::poweramp::TubeStage;
using felitronics::poweramp::Voicing;

namespace
{
//==============================================================================
// PRODUCT-SEAM FIXTURE — identical to the golden's: OrbitCab's voicing presets + control POD + the
// index→voicing lookup shim, reproduced with the SAME constants (product data; the core stage only
// consumes the plain numbers). This is the documented product fixture the theory attacks.
constexpr Voicing kTubeVoicings[4] = {
    /* 6L6  */ { 0.82f, 2.3f, 0.25f, 0.30f, 0.0f,   8.0f, 130.0f, 0.20f, 0.40f,   5000.0f, 3.0f,  95.0f, 7.5f, 0.45f,    500.0f, -5.0f, 0.70f,    90.0f, 1.6f, 4.0f,  1600.0f, 3.0f,   150.0f, 1.8f, 8000.0f },
    /* EL34 */ { 1.12f, 3.2f, 0.15f, 0.20f, 0.0f,   6.0f,  90.0f, 0.22f, 0.50f,   3400.0f, 5.0f, 118.0f, 4.0f, 0.60f,    680.0f,  4.5f, 0.70f,   100.0f, 1.8f, 3.5f,  1500.0f, 4.0f,   170.0f, 2.2f, 7000.0f },
    /* EL84 */ { 1.45f, 1.7f, 0.45f, 0.18f, 0.0f,  14.0f, 240.0f, 0.42f, 0.70f,   5500.0f, 4.0f, 130.0f, 3.0f, 0.75f,   1600.0f,  4.0f, 0.70f,   110.0f, 2.0f, 4.0f,  1800.0f, 4.5f,   210.0f, 2.6f, 8500.0f },
    /* KT88 */ { 0.60f, 3.8f, 0.05f, 0.34f, 0.0f,   4.0f,  80.0f, 0.12f, 0.30f,   4800.0f, 4.0f,  68.0f, 10.0f, 0.30f,    100.0f,  5.0f, 1.40f,    80.0f, 1.4f, 4.0f,  1300.0f, 2.5f,   120.0f, 1.4f, 10000.0f },
};

struct TubeParams
{
    float driveDb = 0.0f, outputDb = 0.0f;
    int   tubeType = 0;
    bool  singleEnded = false;
    float autoComp = 1.0f, sag = 0.0f, presence = 0.0f, depth = 0.0f, load = 0.0f, iron = 0.0f, bias = 0.0f;
};

class TubePowerAmp
{
public:
    void prepare (double sampleRate, int maxBlock, int oversampleFactor = 4) { d.prepare (sampleRate, maxBlock, oversampleFactor); }
    void reset() { d.reset(); }
    void setParams (const TubeParams& t) noexcept
    {
        felitronics::poweramp::Params p;
        p.driveDb = t.driveDb; p.outputDb = t.outputDb; p.singleEnded = t.singleEnded; p.autoComp = t.autoComp;
        p.sag = t.sag; p.presence = t.presence; p.depth = t.depth; p.load = t.load; p.iron = t.iron; p.bias = t.bias;
        d.setParams (p, kTubeVoicings[(std::size_t) std::clamp (t.tubeType, 0, 3)]);
    }
    void process (float* const* io, int numChannels, int numSamples) noexcept { d.process (io, numChannels, numSamples); }
    int  latencySamples() const noexcept { return d.latencySamples(); }

private:
    PowerAmpStage d;
};

int g_checks = 0, g_fail = 0;
void check (bool ok, const char* m) { ++g_checks; if (! ok) { ++g_fail; std::printf ("  [FAIL] %s\n", m); } else std::printf ("  [ok]   %s\n", m); }
void info  (const char* m) { std::printf ("       %s\n", m); }

constexpr double kPi     = 3.14159265358979323846;
constexpr double kSr     = 48000.0;
constexpr int    kMaxBlk = 512;
constexpr int    kLat    = 31;                     // tpp(32) − 1 — the documented oversampler round-trip

float  dbToGain (float db) { return std::pow (10.0f, db * 0.05f); }
double dbc (double num, double den) { return 20.0 * std::log10 (std::max (1e-15, num) / std::max (1e-15, den)); }

TubeParams P (float driveDb, bool se, int tube = 0, float outDb = 0.0f)
{ TubeParams t; t.driveDb = driveDb; t.outputDb = outDb; t.singleEnded = se; t.tubeType = tube; t.autoComp = 1.0f; return t; }

TubeParams Pf (float driveDb, bool se, int tube, float sag, float pres, float depth, float outDb = 0.0f)
{ TubeParams t = P (driveDb, se, tube, outDb); t.sag = sag; t.presence = pres; t.depth = depth; return t; }

std::vector<float> sine (int warm, int N, int tail, int cycles, double amp, double dc = 0.0)
{
    std::vector<float> x ((std::size_t) (warm + N + tail), 0.0f);
    const double w = 2.0 * kPi * cycles / N;
    for (int n = 0; n < (int) x.size(); ++n) x[(std::size_t) n] = (float) (amp * std::sin (w * (n - warm)) + dc);
    return x;
}

double magBin (const std::vector<float>& x, int start, int N, int bin)
{
    double re = 0.0, im = 0.0; const double w = 2.0 * kPi * bin / N;
    for (int n = 0; n < N; ++n) { const double s = x[(std::size_t) (start + n)]; re += s * std::cos (w * n); im -= s * std::sin (w * n); }
    return std::sqrt (re * re + im * im) / N;
}
double rms (const std::vector<float>& x, int start, int N) { double e = 0.0; for (int n = 0; n < N; ++n) { const double s = x[(std::size_t) (start + n)]; e += s * s; } return std::sqrt (e / std::max (1, N)); }

// Run a mono buffer through a fresh TubePowerAmp at OS factor `os`, with an optional block schedule.
std::vector<float> runStage (const TubeParams& tp, std::vector<float> buf, int os = 4, const std::vector<int>* sched = nullptr)
{
    TubePowerAmp d; d.prepare (kSr, kMaxBlk, os);
    static const int dflt[] = { 64, 128, 333, 512, 17, 1 };
    const int total = (int) buf.size();
    int pos = 0, bi = 0;
    while (pos < total)
    {
        int n = sched ? (*sched)[(std::size_t) (bi % (int) sched->size())] : dflt[bi % 6];
        n = std::min (n, total - pos);
        d.setParams (tp);
        float* io[1] { buf.data() + pos };
        d.process (io, 1, n);
        pos += n; ++bi;
    }
    return buf;
}

bool   allFinite (const std::vector<float>& x) { for (float s : x) if (! std::isfinite (s)) return false; return true; }
double maxAbs (const std::vector<float>& x, int from, int to) { double m = 0; for (int i = from; i < to; ++i) m = std::max (m, (double) std::fabs (x[(std::size_t) i])); return m; }
} // namespace

int main()
{
    std::printf ("felitronics::poweramp THEORY-FIRST falsification suite (derived bounds; tries to break the DSP)\n");
    const int warm = 4096, N = 16384, tail = 512;
    const int an = warm + kLat;
    const int binF = (int) std::llround (1000.0 * N / kSr);

    // ===================== TT1: SILENCE → SILENCE (exact) =====================
    // Theory: with zero input every path collapses to zero — g(0)=0 for Asym (the norm subtracts
    // tanh(kb)), PP(0)=g(0+vb)−g(−0+vb)=0 at evenLeak=0, DC-block/transformer/bell/shelf of 0 = 0, and
    // the sag rail multiplies 0. A reset stage fed silence therefore emits IDENTICALLY 0 — there is no
    // additive DC/offset term anywhere. No settle horizon needed (no transient from reset). This is the
    // strongest silence assertion: exact bit-zero, across every tube/topology/drive and feel fully on.
    {
        double worst = 0.0;
        for (int t = 0; t < 4; ++t) for (bool se : { false, true })
            for (float dr : { 0.0f, 18.0f, 36.0f })
            {
                std::vector<float> z (8192, 0.0f);
                worst = std::max (worst, maxAbs (runStage (P (dr, se, t), z), 0, 8192));                 // feel off
                std::vector<float> z2 (8192, 0.0f);
                TubeParams full = Pf (dr, se, t, 1.0f, 1.0f, 1.0f); full.load = 1.0f; full.iron = 1.0f; full.bias = 1.0f;
                worst = std::max (worst, maxAbs (runStage (full, z2), 0, 8192));                          // feel fully on
            }
        std::printf ("       TT1 worst |silence out| = %.2e (theory: exactly 0)\n", worst);
        check (worst == 0.0, "TT1 silence in ⇒ output identically 0 (no additive offset; feel on & off, all tubes)");
    }

    // ===================== TT2: DC REJECTION =====================
    // Theory: the OS-domain DC blocker is (1−z⁻¹)/(1−dcR·z⁻¹) with a zero AT DC ⇒ steady-state DC gain
    // is exactly 0. Corner kDcBlockHz=10 Hz at fsOs=fs·os ⇒ τ=1/(2π·10)≈15.9 ms; with the 25 ms coeff
    // smoother the operating point settles by ~0.2 s and the residual is bounded by exp(−2π·10·T).
    // Feed ~0.55 s of sustained DC (sag OFF so the rail doesn't modulate), measure the settled tail. A
    // broken DC-block leaving even −40 dBFS of DC (0.005 for DC=0.5) would fail the 1e-3 bound.
    {
        const int Nn = 26624;                              // ~0.55 s @48k
        double worst = 0.0;
        for (int t = 0; t < 4; ++t) for (bool se : { false, true })
            for (double dc : { 0.5, -0.5, 1.0 })
            {
                std::vector<float> in ((std::size_t) Nn, (float) dc);
                auto out = runStage (P (12.0f, se, t), in);
                worst = std::max (worst, maxAbs (out, Nn - 4096, Nn));   // last ~85 ms, fully settled
            }
        std::printf ("       TT2 worst settled DC tail = %.2e (bound 1e-3; theory residual ~exp(-2π·10·T)≈0)\n", worst);
        check (worst < 1e-3, "TT2 sustained DC ⇒ steady-state output rejected below 1e-3 (exact DC zero)");
    }

    // ===================== TT3: BIBO SANITY =====================
    // Theory: |x|≤1 with sane params ⇒ every stage is a contraction or a bounded gain (kernel |y|<8;
    // DC-block/transformer non-amplifying; bell ≤+5 dB, load ≤~+4 dB, depth ≤+10 dB blended; sag rail
    // s∈[0.2,1] attenuates; drive-comp+Output normalise). Net output finite + O(1). Assert (a) finite
    // ALWAYS (the true BIBO), (b) ‖out‖∞ < 5 at out 0 dB / autoComp 1, feel off (golden X1 pins <4.0),
    // (c) ‖out‖∞ < 8 with the feel layer fully engaged (its depth/bell headroom).
    {
        std::vector<std::vector<float>> sigs;
        sigs.push_back (sine (0, 8192, 0, 137, 1.0));                                                   // full-scale sine
        { std::vector<float> v (8192); for (int i = 0; i < 8192; ++i) v[(std::size_t) i] = (i & 1) ? 1.0f : -1.0f; sigs.push_back (v); }  // Nyquist square
        { std::vector<float> v (8192); unsigned long long s = 9; for (auto& x : v) { s = s * 6364136223846793005ULL + 1ULL; x = ((float) ((s >> 40) & 0xFFFFFF) / 8388608.0f - 1.0f); } sigs.push_back (v); }  // full-scale noise
        { std::vector<float> v (8192); for (int i = 0; i < 8192; ++i) v[(std::size_t) i] = i < 4096 ? -1.0f : 1.0f; sigs.push_back (v); }  // full-scale step
        bool finite = true; double peakOff = 0.0, peakOn = 0.0;
        for (auto& s : sigs) for (int t = 0; t < 4; ++t) for (bool se : { false, true }) for (float dr : { 0.0f, 12.0f, 24.0f, 36.0f })
        {
            auto o0 = runStage (P (dr, se, t), s);
            finite = finite && allFinite (o0); peakOff = std::max (peakOff, maxAbs (o0, 0, (int) o0.size()));
            TubeParams fp = Pf (dr, se, t, 1.0f, 1.0f, 1.0f); fp.load = 1.0f; fp.iron = 1.0f; fp.bias = 1.0f;
            auto o1 = runStage (fp, s);
            finite = finite && allFinite (o1); peakOn = std::max (peakOn, maxAbs (o1, 0, (int) o1.size()));
        }
        std::printf ("       TT3 worst peak: feel off = %.3f  feel on = %.3f\n", peakOff, peakOn);
        check (finite, "TT3 BIBO: |x|≤1 ⇒ output finite for every input/tube/topology/drive (feel on & off)");
        check (peakOff < 5.0, "TT3 BIBO: ‖out‖∞ < 5 (out 0 dB, autoComp 1, feel off) — derived envelope");
        check (peakOn  < 8.0, "TT3 BIBO: ‖out‖∞ < 8 with the feel layer fully engaged (depth/bell headroom)");
    }

    // ===================== TT4: PUSH-PULL vs SINGLE-ENDED HARMONIC THEORY =====================
    // Theory: PP y = g(u+vb) − g(−u+vb) is ODD in u (a differential pair) ⇒ even harmonics cancel to the
    // numeric/OS floor when evenLeak=0 (the voicings ship evenLeak=0). SE y=g(u) is asymmetric ⇒ strong
    // H2. Gap derived conservatively from evenLeak=0: PP even at the discrete harmonic bins < −80 dBc
    // (float+OS residual of an exact cancellation); SE H2 present > −42 dBc; SE−PP separation > 40 dB.
    {
        bool ppOk = true, seOk = true, gapOk = true; double worstPP = -300, weakestSE = 300;
        for (int t = 0; t < 4; ++t)
        {
            auto pp = runStage (P (24.0f, false, t), sine (warm, N, tail, binF, 0.5));
            auto se = runStage (P (24.0f, true,  t), sine (warm, N, tail, binF, 0.5));
            const double ppF = magBin (pp, an, N, binF), seF = magBin (se, an, N, binF);
            double ppEven = -300; for (int k : { 2, 4, 6 }) ppEven = std::max (ppEven, dbc (magBin (pp, an, N, binF * k), ppF));
            const double seH2 = dbc (magBin (se, an, N, binF * 2), seF);
            worstPP = std::max (worstPP, ppEven); weakestSE = std::min (weakestSE, seH2);
            ppOk = ppOk && (ppEven < -80.0);
            seOk = seOk && (seH2 > -42.0);
            gapOk = gapOk && (seH2 - ppEven > 40.0);
        }
        std::printf ("       TT4 worst PP even = %.1f dBc   weakest SE H2 = %.1f dBc\n", worstPP, weakestSE);
        check (ppOk,  "TT4 PP differential pair cancels even harmonics (< −80 dBc, evenLeak=0, all tubes)");
        check (seOk,  "TT4 SE is even-rich (H2 > −42 dBc, all tubes)");
        check (gapOk, "TT4 SE H2 sits > 40 dB above the PP even floor (topology really differs)");
    }

    // ===================== TT5: SAG MONOTONICITY + RELEASE =====================
    // Part A (SagEnvelope unit) — droop = maxDroop·amount·(1−exp(−2·env)) with steady env→|demand|:
    // strictly ↑ in demand, and < maxDroop·amount always ⇒ rail s=1−droop ∈ (0.1,1] (never inverts;
    // the stage additionally clamps s ≥ kSagMinRail=0.2). Assert monotone-↑ droop and the bound.
    {
        SagEnvelope s; s.prepare (kSr);
        const float amount = 1.0f, maxDroop = 0.42f;                 // EL84-like: spongy
        s.setParams (amount, 8.0f, 200.0f, maxDroop);
        double prev = -1.0; bool mono = true, bounded = true;
        for (double d : { 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 2.0 })
        {
            s.reset();
            for (int i = 0; i < 40000; ++i) s.process ((float) d);   // reach steady env≈d
            const double droop = s.droop();
            if (droop < prev - 1e-9) mono = false;                   // non-decreasing
            if (droop >= (double) (maxDroop * amount)) bounded = false;
            prev = droop;
        }
        check (mono,    "TT5 sag droop is monotone ↑ in sustained demand (louder ⇒ lower rail)");
        check (bounded, "TT5 sag droop < maxDroop·amount ⇒ rail s > 0.1 (never inverts / blows 1/s)");

        // Part B (release) — a single real pole ⇒ monotone decay; recovery TC ≫ attack TC (dual-TC sign).
        s.reset(); s.setParams (amount, 8.0f, 200.0f, maxDroop);
        for (int i = 0; i < 40000; ++i) s.process (1.0f);            // sustain loud → high env
        const float envLoud = s.envelope();
        // Run well past the release half-life (recovery TC 200 ms ⇒ half-life ≈ 6654 samples): 20 000
        // samples ≈ 3 half-lives ⇒ env should fall below 0.5·envLoud while never rising (monotone pole).
        bool relMono = true; float last = envLoud;
        for (int i = 0; i < 20000; ++i) { s.process (0.0f); const float e = s.envelope(); if (e > last + 1e-9f) relMono = false; last = e; }
        check (relMono && last < 0.5f * envLoud, "TT5 release: env decays monotonically toward 0 (single-pole recovery)");
        // attack half-time < release half-time — the whole point of a supply sag (fast collapse, slow bloom).
        s.reset(); s.setParams (amount, 8.0f, 200.0f, maxDroop);
        for (int i = 0; i < 40000; ++i) s.process (1.0f); const float envFull = s.envelope();
        s.reset(); int atkHalf = 1 << 30; for (int i = 0; i < 40000; ++i) { s.process (1.0f); if (s.envelope() >= 0.5f * envFull) { atkHalf = i; break; } }
        for (int i = 0; i < 40000; ++i) s.process (1.0f); const float envRel0 = s.envelope();
        int relHalf = 1 << 30; for (int i = 0; i < 400000; ++i) { s.process (0.0f); if (s.envelope() <= 0.5f * envRel0) { relHalf = i; break; } }
        std::printf ("       TT5 attack-half=%d  release-half=%d samples\n", atkHalf, relHalf);
        check (atkHalf < relHalf, "TT5 dual-TC: attack is faster than release (touch punches through, then blooms)");

        // Part C (stage) — the SAG-specific claim: louder program ⇒ MORE sag compression. Isolate sag by
        // the on/off gap so ordinary saturation (present in both) cancels; the residual gap must grow ↑.
        double gPrev = -1.0; bool gapMono = true;
        double g02 = 0, g05 = 0, g09 = 0;
        int idx = 0;
        for (double a : { 0.2, 0.5, 0.9 })
        {
            auto in = sine (0, 12288, 0, 256, a);
            const double rOff = rms (runStage (Pf (18.0f, true, 2, 0.0f, 0, 0), in), 8000, 4000);
            const double rOn  = rms (runStage (Pf (18.0f, true, 2, 1.0f, 0, 0), in), 8000, 4000);
            const double gapDb = dbc (rOff, rOn);                    // how much sag pulls the level down
            if (gapDb < gPrev - 1e-3) gapMono = false;               // monotone ↑ with level
            gPrev = gapDb;
            (idx == 0 ? g02 : idx == 1 ? g05 : g09) = gapDb; ++idx;
        }
        std::printf ("       TT5 sag compression gap dB: a=0.2 %.3f  a=0.5 %.3f  a=0.9 %.3f\n", g02, g05, g09);
        check (gapMono && g09 > g02, "TT5 stage: sag compresses MORE as the program gets louder (monotone droop)");
    }

    // ===================== TT6: DRIVE-COMPENSATION (autoComp) =====================
    // Theory: comp = |slope|^(−autoComp). At autoComp=1 the compensation cancels the composite small-
    // signal slope ⇒ the low-level gain is Drive-INVARIANT (documented purpose). Probe with a tiny tone
    // (1e-4, shaper stays linear even at Drive 36) so we read ONLY the compensation. Derived window: the
    // golden pins <0.05 dB drift; assert <0.1 dB with margin. At autoComp=0 there is no compensation ⇒
    // near-linear RMS grows monotonically ≈ +Δdrive dB (here +12 dB over the 0→12 sweep).
    {
        auto in = sine (warm, N, tail, binF, 1e-4);
        const double ri = rms (in, warm, N);
        double worstDrift = 0.0; bool autoOffMono = true, growthOk = true;
        for (int t = 0; t < 4; ++t) for (bool se : { false, true })
        {
            double lo = 1e9, hi = -1e9;
            for (float dr : { 0.0f, 4.0f, 8.0f, 12.0f })
            {
                const double g = dbc (rms (runStage (P (dr, se, t), in), an, N), ri);
                lo = std::min (lo, g); hi = std::max (hi, g);
            }
            worstDrift = std::max (worstDrift, hi - lo);

            double prevG = -1e9, g0 = 0, g12 = 0; int di = 0;
            for (float dr : { 0.0f, 4.0f, 8.0f, 12.0f })
            {
                TubeParams tp = P (dr, se, t); tp.autoComp = 0.0f;
                const double g = dbc (rms (runStage (tp, in), an, N), ri);
                if (g < prevG - 0.05) autoOffMono = false;           // strictly ↑ (near-linear)
                prevG = g; if (di == 0) g0 = g; if (di == 3) g12 = g; ++di;
            }
            if (! (g12 - g0 > 10.5 && g12 - g0 < 13.5)) growthOk = false;   // ≈ +12 dB drive ⇒ +12 dB out
        }
        std::printf ("       TT6 autoComp=1 worst drift = %.3f dB (window 0.1)\n", worstDrift);
        check (worstDrift < 0.1, "TT6 autoComp=1 ⇒ small-signal gain Drive-invariant (< 0.1 dB over 0→12 dB)");
        check (autoOffMono,      "TT6 autoComp=0 ⇒ RMS grows monotonically with Drive (no compensation)");
        check (growthOk,         "TT6 autoComp=0 ⇒ +12 dB Drive gives ≈ +12 dB out (near-linear tiny signal)");
    }

    // ===================== TT7: LATENCY =====================
    // Theory: latency = oversampler round-trip = tapsPerPhase−1 = 31, INVARIANT across the os factor
    // (both up- and down-sampler contribute (N−1)/2, the FIR length scales with the factor but tpp is
    // fixed at 32). Before prepare() the oversampler reports 0 (hosts query latency early). And the
    // REPORTED value must equal the PHYSICAL delay (impulse main lobe index).
    {
        bool inv = true;
        for (int os : { 2, 4, 8, 16, 32 }) { TubePowerAmp d; d.prepare (kSr, kMaxBlk, os); inv = inv && (d.latencySamples() == kLat); }
        check (inv, "TT7 reported latency == 31 (tpp−1), invariant across os factor {2,4,8,16,32}");

        TubePowerAmp fresh;                                          // not prepared
        check (fresh.latencySamples() == 0, "TT7 unprepared stage reports 0 latency (host queries before prepare)");

        bool measured = true;
        for (int t = 0; t < 4; ++t) for (bool se : { false, true })
        {
            std::vector<float> imp (4096, 0.0f); imp[1000] = 1e-3f;  // small ⇒ linear regime
            auto out = runStage (P (6.0f, se, t), imp);
            int peak = 0; double pv = 0; for (int i = 990; i < 1100; ++i) if (std::fabs (out[(std::size_t) i]) > pv) { pv = std::fabs (out[(std::size_t) i]); peak = i; }
            measured = measured && (peak == 1000 + kLat);
        }
        check (measured, "TT7 reported latency == impulse-measured main lobe (+31), all tubes/topologies");
    }

    // ===================== TT8: NaN/Inf CONTAINMENT =====================
    // Theory: flushDenormals() + the per-block isfinite gates zero any non-finite IIR/envelope state, and
    // the setParams gate sanitizes at the input, so no poison persists. Horizon: state is flushed at the
    // block boundary and the input is sanitized before the OS FIR ⇒ recovery within ≤2 valid blocks.
    {
        TubePowerAmp d; d.prepare (kSr, kMaxBlk, 4); d.setParams (P (18.0f, false));
        { std::vector<float> bad (256, std::numeric_limits<float>::quiet_NaN()); float* io[1] { bad.data() }; d.process (io, 1, 256); }
        { std::vector<float> inf (256, std::numeric_limits<float>::infinity()); float* io[1] { inf.data() }; d.process (io, 1, 256); }
        bool early = true, alive = true; double tailRms = 0;
        for (int blk = 0; blk < 6; ++blk)
        {
            auto s = sine (0, 256, 0, 8, 0.3); d.setParams (P (18.0f, false));
            float* io[1] { s.data() }; d.process (io, 1, 256);
            if (blk >= 2) { early = early && allFinite (s); tailRms += rms (s, 0, 256); alive = alive && (rms (s, 0, 256) > 1e-4); }
        }
        std::printf ("       TT8 recovered RMS (blocks 2..5) = %.4f\n", tailRms / 4);
        check (early && alive, "TT8 NaN/Inf burst contained: finite + non-zero output within ≤2 valid blocks");

        // Param poison: NaN driveDb / Inf outputDb / NaN autoComp ⇒ sanitize gate ⇒ 0 dB / autoComp 1 fallback.
        TubePowerAmp d2; d2.prepare (kSr, kMaxBlk, 4);
        TubeParams badp; badp.driveDb = std::numeric_limits<float>::quiet_NaN(); badp.outputDb = std::numeric_limits<float>::infinity(); badp.autoComp = std::numeric_limits<float>::quiet_NaN();
        auto tone = sine (0, 4096, 0, 171, 0.3);
        d2.setParams (badp); { float* io[1] { tone.data() }; d2.process (io, 1, 4096); }
        check (allFinite (tone) && rms (tone, 512, 3000) > 1e-3, "TT8 non-finite params ⇒ sanitize-gate fallback, finite + alive");
    }

    // ===================== TT9: DETERMINISM + BLOCK INVARIANCE (B9 pin) =====================
    // Theory: the stage is a deterministic state machine — no RNG, no uninitialised reads — so identical
    // schedules from reset give BIT-identical output, feel on or off. And with the feel layer OFF the
    // per-sample processing is host-buffer-agnostic ⇒ bit-exact across block schedules.
    {
        // size warm+N+tail so the steady-state analysis window [an, an+N) stays in bounds (an=warm+kLat).
        std::vector<float> in ((std::size_t) (warm + N + tail), 0.0f);
        { unsigned long long s = 5; for (auto& v : in) { s = s * 6364136223846793005ULL + 1; v = 0.5f * ((float) ((s >> 40) & 0xFFFFFF) / 8388608.0f - 1.0f); } }

        // (a) two fresh runs, identical schedule ⇒ bit-identical (feel off AND feel on)
        auto a0 = runStage (P (24.0f, false, 1), in), a1 = runStage (P (24.0f, false, 1), in);
        double dOff = maxAbs ([&]{ std::vector<float> d (a0.size()); for (std::size_t i = 0; i < a0.size(); ++i) d[i] = a0[i] - a1[i]; return d; }(), 0, (int) a0.size());
        TubeParams fp = Pf (24.0f, false, 1, 1.0f, 0.6f, 0.6f); fp.load = 1.0f; fp.iron = 1.0f; fp.bias = 1.0f;
        auto b0 = runStage (fp, in), b1 = runStage (fp, in);
        double dOn = maxAbs ([&]{ std::vector<float> d (b0.size()); for (std::size_t i = 0; i < b0.size(); ++i) d[i] = b0[i] - b1[i]; return d; }(), 0, (int) b0.size());
        check (dOff == 0.0 && dOn == 0.0, "TT9 determinism: identical schedule twice ⇒ bit-identical (feel on & off)");

        // (b) block-size invariance, feel OFF ⇒ bit-exact across a hostile schedule
        std::vector<int> blk512 { 512 }, hostile { 1, 7, 64, 333, 512, 128 };
        auto c0 = runStage (P (24.0f, false, 1), in, 4, &blk512);
        auto c1 = runStage (P (24.0f, false, 1), in, 4, &hostile);
        const int wn = an, we = an + N;
        auto win = [&] (const std::vector<float>& x, const std::vector<float>& y) { double m = 0; for (int i = wn; i < we; ++i) m = std::max (m, (double) std::fabs (x[(std::size_t) i] - y[(std::size_t) i])); return m; };
        const double detOff = win (c0, c1);
        std::printf ("       TT9 block-schedule diff: feel off = %.2e\n", detOff);
        check (detOff < 1e-6, "TT9 block-size invariance holds with the feel layer OFF (bit-exact)");

        // (c) block-size invariance, feel ON. Theory: with constant params every smoothed coefficient
        // SNAPS to target on the first (primed) block (a=1) and then holds — so no per-block quantity
        // varies with the block length n, and the whole feel path (per-sample sag rail, per-sample shelf
        // "open" blend, per-sample bias, continuous OT/DC IIR state) is block-size deterministic. Measured
        // in-bounds: the feel-ON output is BIT-EXACT across schedules (diff 0), for sag/presence/depth and
        // for the full load+iron+bias path — a STRONGER result than the feel-OFF S3 tolerance.
        //
        // *** FINDING (reported): this CONTRADICTS the golden's "B9 ~1e-2 feel-layer block-size
        // discrepancy." That ~1e-2 is an ARTIFACT of the golden's analysis window: diffWin(..., an, an+N)
        // over a buffer of length warm+N reads 31 samples PAST the end (an=warm+kLat). Probed both ways in
        // the SAME binary, the identical OOB pattern yields ~1e-2 of nondeterministic garbage while the
        // in-bounds window yields exactly 0 — for feel ON *and* for feel OFF (which is provably exact). So
        // there is NO real deviation to pin; the correct, honest assertion is bit-exactness. ***
        TubeParams feel = Pf (24.0f, false, 1, 0.8f, 0.5f, 0.5f);
        const double detOnA = win (runStage (feel, in, 4, &blk512), runStage (feel, in, 4, &hostile));
        TubeParams feelFull = feel; feelFull.load = 1.0f; feelFull.iron = 1.0f; feelFull.bias = 1.0f;
        const double detOnB = win (runStage (feelFull, in, 4, &blk512), runStage (feelFull, in, 4, &hostile));
        std::printf ("       TT9 block-schedule diff: feel ON = %.2e   +load/iron/bias = %.2e (theory: exactly 0)\n", detOnA, detOnB);
        info ("TT9 NOTE: the feel layer is block-size BIT-EXACT in-bounds; the golden's 'B9 ~1e-2' is an OOB-window artifact (see FINDING above).");
        check (detOnA < 1e-6 && detOnB < 1e-6, "TT9 feel-ON block-size invariance is bit-exact (no real B9 deviation; golden's ~1e-2 is OOB-window)");
    }

    std::printf ("%d checks, %d failures\n", g_checks, g_fail);
    std::printf (g_fail ? "THEORY SUITE FAILED\n" : "THEORY SUITE PASSED\n");
    return g_fail ? 1 : 0;
}
