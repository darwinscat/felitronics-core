// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics_chorus_render — a bass DI through chorus::BassChorus at several slider positions, written as
// 24-bit stereo WAVs for listening (the dry too, so the two can be A/B'd at the same level).
//
//   felitronics_chorus_render <outDir> [in.wav ...]
//
// With input files, each is run as given (a mono file is put on both sides of the stereo bus). Without any, a
// bass line is synthesised (Karplus–Strong: fingered notes from low B 30.9 Hz up to G2, a slap/pop stretch and
// long held notes, 24 s) and written as `chorus-dry.wav`. Outputs: `chorus-amount-025.wav` … `-100.wav` at the
// default menu (0.7 Hz, 200 Hz crossover), `chorus-amount-075-xover-500.wav`, and README.txt.
#include <felitronics/chorus/BassChorus.h>
#include <felitronics/chorus/ChorusLaw.h>
#include <felitronics/io/Wav.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
constexpr double kPi = 3.14159265358979323846;
using Planes = std::vector<std::vector<double>>;

struct Rng
{
    std::uint32_t s;
    double next() { s = s * 1664525u + 1013904223u; return (double) (s >> 8) / 8388608.0 - 1.0; }
};

enum class Style { finger, thumb, pop };

struct Note { double at, hz, dur, vel; Style style; };

/** One plucked string, Karplus–Strong: an excitation shaped by the pluck, a two-point loop filter, a release. */
void pluck (std::vector<double>& out, double fs, const Note& n, Rng& rng)
{
    const int P = std::max (2, (int) std::lround (fs / n.hz - 0.5));   // the loop filter adds about half a sample
    std::vector<double> line ((std::size_t) P, 0.0);

    // The excitation is the string's displacement at the pluck: a triangle peaking at the pluck point, whose
    // harmonics fall at 12 dB/oct (n-th as 1/n^2) as a real plucked string's do — the fundamental leads, as on a
    // bass DI. A little low-passed noise on top for the finger's grain; more, and brighter, for a thumb or a pop,
    // which also strike nearer the bridge.
    const double beta  = n.style == Style::finger ? 0.22 : (n.style == Style::thumb ? 0.10 : 0.06);
    const double grain = n.style == Style::finger ? 0.04 : (n.style == Style::thumb ? 0.12 : 0.20);
    const double lpHz  = n.style == Style::finger ? 1500.0 : 6000.0;
    const double a = 1.0 - std::exp (-2.0 * kPi * lpHz / fs);
    const int    peakAt = std::max (1, (int) (beta * P));
    double lp = 0.0, mean = 0.0;
    for (int i = 0; i < P; ++i)
    {
        lp += a * (rng.next() - lp);
        const double tri = i <= peakAt ? (double) i / peakAt : (double) (P - i) / (double) (P - peakAt);
        line[(std::size_t) i] = tri + grain * lp;
        mean += line[(std::size_t) i];
    }
    mean /= P;
    double peak = 1e-9;
    for (auto& v : line)
    {
        v -= mean;   // no DC in the string
        peak = std::max (peak, std::fabs (v));
    }
    for (auto& v : line)
        v *= n.vel / peak;

    // Sustain: a T60 of ~5 s on the low notes (shorter up the neck); a muted release of 70 ms after `dur`.
    const double t60   = n.style == Style::finger ? 5.0 * std::sqrt (41.2 / n.hz) : 2.5;
    const double gHold = std::pow (10.0, -3.0 / (t60 * n.hz));
    const double gMute = std::pow (10.0, -3.0 / (0.07 * n.hz));
    const double bright = n.style == Style::finger ? 0.5 : 0.35;   // a smaller second tap keeps more highs
    const int start = (int) (n.at * fs), hold = (int) (n.dur * fs), total = hold + (int) (0.4 * fs);

    int w = 0;
    for (int i = 0; i < total && start + i < (int) out.size(); ++i)
    {
        const int    w1 = w + 1 == P ? 0 : w + 1;
        const double y  = line[(std::size_t) w];
        out[(std::size_t) (start + i)] += y;
        const double g = i < hold ? gHold : gMute;
        line[(std::size_t) w] = g * ((1.0 - bright) * y + bright * line[(std::size_t) w1]);
        w = w1;
    }

    // A thumb thumps the fret, a pop snaps it: a short bright click on top of the string.
    if (n.style != Style::finger)
    {
        const int len = (int) (0.004 * fs);
        double hp = 0.0, prev = 0.0;
        for (int i = 0; i < len && start + i < (int) out.size(); ++i)
        {
            const double x = rng.next();
            hp = 0.9 * (hp + x - prev);
            prev = x;
            out[(std::size_t) (start + i)] += (n.style == Style::pop ? 0.5 : 0.3) * n.vel * hp * (1.0 - (double) i / len);
        }
    }
}

std::vector<double> synthBass (double fs)
{
    const double B0 = 30.87, D1 = 36.71, E1 = 41.20, Fs1 = 46.25, G1 = 49.00, A1 = 55.00, B1 = 61.74,
                 Cs2 = 69.30, D2 = 73.42, E2 = 82.41, Fs2 = 92.50, G2 = 98.00;
    const double e = 0.3;   // an eighth at 100 BPM
    std::vector<Note> notes;
    auto bar = [&] (double t0, std::initializer_list<double> hz, Style st = Style::finger)
    {
        int k = 0;
        for (double f : hz)
        {
            notes.push_back ({ t0 + k * e, f, e * 0.9, (k % 4 == 0) ? 1.0 : 0.8, st });
            ++k;
        }
    };
    // Fingered groove down to low B and up to G2.
    bar (0.0 * e, { B0, B0, B0, D1, E1, E1, Fs1, A1 });
    bar (8.0 * e, { E1, E1, E1, G1, A1, B1, A1, G1 });
    bar (16.0 * e, { A1, A1, E1, A1, D2, D2, Cs2, A1 });
    bar (24.0 * e, { G1, G1, D2, G1, G2, Fs2, E2, D2 });
    // Slap: thumb on the low notes, pops on the octaves.
    const struct { double hz; Style st; } slap[] = {
        { E1, Style::thumb }, { E2, Style::pop }, { E1, Style::thumb }, { E1, Style::thumb },
        { D2, Style::pop },   { E1, Style::thumb }, { G1, Style::thumb }, { G2, Style::pop },
        { A1, Style::thumb }, { G2, Style::pop }, { A1, Style::thumb }, { E2, Style::pop },
        { D2, Style::pop },   { B1, Style::thumb }, { A1, Style::thumb }, { G1, Style::thumb } };
    for (int k = 0; k < 16; ++k)
        notes.push_back ({ (32.0 + k) * e, slap[k].hz, e * 0.6, slap[k].st == Style::pop ? 1.0 : 0.95, slap[k].st });
    // Long held notes, where a chorus is heard best.
    const double held[] = { B0, E1, G2, D2 };
    for (int k = 0; k < 4; ++k)
        notes.push_back ({ (48.0 + 8.0 * k) * e, held[k], 7.5 * e, 1.0, Style::finger });

    const double len = (48.0 + 32.0) * e + 0.5;   // 24.5 s
    std::vector<double> out ((std::size_t) (len * fs), 0.0);
    Rng rng { 0x0BA55u };
    for (const auto& n : notes)
        pluck (out, fs, n, rng);

    // The pickup and the DI: a gentle 2-pole low pass at 5 kHz, a DC blocker, -12 dBFS peak.
    const double a = 1.0 - std::exp (-2.0 * kPi * 5000.0 / fs);
    double s1 = 0.0, s2 = 0.0, xPrev = 0.0, dc = 0.0, peak = 1e-9;
    for (auto& v : out)
    {
        s1 += a * (v - s1);
        s2 += a * (s1 - s2);
        dc = 0.9995 * dc + s2 - xPrev;
        xPrev = s2;
        v = dc;
        peak = std::max (peak, std::fabs (v));
    }
    const double g = std::pow (10.0, -12.0 / 20.0) / peak;
    for (auto& v : out)
        v *= g;
    return out;
}

double peakDb (const Planes& p)
{
    double pk = 1e-12;
    for (const auto& c : p)
        for (double v : c)
            pk = std::max (pk, std::fabs (v));
    return 20.0 * std::log10 (pk);
}

/** The stereo bus through the stage at `amount`, with `menu` as the right-click menu stands. */
Planes runChorus (const Planes& in, double fs, float amount, felitronics::chorus::BassChorusParams menu)
{
    felitronics::chorus::BassChorus s;
    s.setParams (felitronics::chorus::paramsForAmount (amount, menu));
    if (! s.prepare (fs, 512, 2))
        return {};
    const std::size_t n = in[0].size() + (std::size_t) (0.1 * fs);
    std::vector<float> l (n, 0.0f), r (n, 0.0f);
    for (std::size_t i = 0; i < in[0].size(); ++i)
    {
        l[i] = (float) in[0][i];
        r[i] = (float) in[in.size() > 1 ? 1 : 0][i];
    }
    for (std::size_t pos = 0; pos < n; pos += 512)
    {
        float* io[2] = { l.data() + pos, r.data() + pos };
        if (! s.process (io, 2, (int) std::min<std::size_t> (512, n - pos)))
            return {};
    }
    Planes out (2, std::vector<double> (n));
    for (std::size_t i = 0; i < n; ++i)
    {
        out[0][i] = l[i];
        out[1][i] = r[i];
    }
    return out;
}

bool write (const std::string& path, const Planes& p, double fs, std::string& readme, const std::string& what)
{
    const bool ok = felitronics::io::writeWav (path, p, fs, 24, false);
    const auto slash = path.find_last_of ('/');
    const std::string name = slash == std::string::npos ? path : path.substr (slash + 1);
    char pk[32];
    std::snprintf (pk, sizeof pk, "%.1f", peakDb (p));
    std::printf ("%s  %s  (peak %s dBFS)\n", ok ? "wrote" : "FAILED", name.c_str(), pk);
    readme += name + "\n    " + what + " Peak " + pk + " dBFS.\n";
    return ok;
}
} // namespace

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf (stderr, "usage: %s <outDir> [in.wav ...]\n", argv[0]);
        return 2;
    }
    const std::string dir = argv[1];
    std::string readme =
        "Bass chorus (felitronics::chorus::BassChorus) — listening renders, 24-bit stereo WAV.\n"
        "One slider, AMOUNT 0..1: depth 0..3 ms of head swing and the voice's level 0..1 (1 = the dry's level)\n"
        "rise together. Menu at its defaults unless named: rate 0.7 Hz, crossover 200 Hz, centre delay 8 ms.\n"
        "Below the crossover the voice falls at 24 dB/oct (-24.6 dB an octave down), so the lows stay the dry's,\n"
        "centred and in phase: attenuated in the voice, not excluded. The two sides carry\n"
        "two voices whose LFOs swing in quadrature (a quarter period apart), so a mono sum still swings.\n\n";

    struct Input { std::string stem; Planes ch; double fs; };
    std::vector<Input> inputs;
    for (int k = 2; k < argc; ++k)
    {
        const auto w = felitronics::io::readWav (argv[k]);
        if (! w.ok)
        {
            std::fprintf (stderr, "cannot read %s: %s\n", argv[k], w.error.c_str());
            return 1;
        }
        std::string stem = argv[k];
        stem = stem.substr (stem.find_last_of ('/') == std::string::npos ? 0 : stem.find_last_of ('/') + 1);
        stem = stem.substr (0, stem.rfind ('.'));
        inputs.push_back ({ "-" + stem, w.ch, w.sr });
    }
    if (inputs.empty())
    {
        const double fs = 48000.0;
        const auto bass = synthBass (fs);
        inputs.push_back ({ "", Planes { bass, bass }, fs });
    }

    bool allOk = true;
    for (const auto& in : inputs)
    {
        const std::string src = in.stem.empty() ? "the synthesised bass line" : "the input" + in.stem;
        Planes dry = in.ch;
        if (dry.size() == 1)
            dry.push_back (dry[0]);
        allOk &= write (dir + "/chorus" + in.stem + "-dry.wav", dry, in.fs, readme,
                        "BEFORE: " + src + " (Karplus-Strong: fingered notes from low B 30.9 Hz up to G2, a slap/pop"
                        " stretch from 9.6 s, long held B0 / E1 / G2 / D2 from 14.4 s), on both sides.");
        felitronics::chorus::BassChorusParams menu;
        for (const float amount : { 0.25f, 0.5f, 0.75f, 1.0f })
        {
            char tag[16];
            std::snprintf (tag, sizeof tag, "%03d", (int) std::lround (amount * 100.0f));
            const auto p = felitronics::chorus::paramsForAmount (amount, menu);
            char what[160];
            std::snprintf (what, sizeof what, "AFTER: amount %.2f (depth %.2f ms, mix %.2f, peak deviation %.1f cents).",
                           (double) amount, (double) p.depthMs, (double) p.mix,
                           (double) felitronics::chorus::pitchDeviationCents (p));
            const Planes out = runChorus (dry, in.fs, amount, menu);
            allOk &= write (dir + "/chorus" + in.stem + "-amount-" + tag + ".wav", out, in.fs, readme, what);
            if (amount > 0.99f && out.size() == 2)
            {
                // What a mono PA or a summed DI hears: (L + R) / 2 on both sides.
                Planes mono = out;
                for (std::size_t i = 0; i < mono[0].size(); ++i)
                    mono[0][i] = mono[1][i] = 0.5 * (out[0][i] + out[1][i]);
                allOk &= write (dir + "/chorus" + in.stem + "-amount-100-mono.wav", mono, in.fs, readme,
                                "AFTER: amount 1.00 folded to mono, (L + R) / 2 on both sides — the quadrature voices still"
                                " swing in the sum.");
            }
        }
        menu.crossoverHz = 500.0f;
        allOk &= write (dir + "/chorus" + in.stem + "-amount-075-xover-500.wav", runChorus (dry, in.fs, 0.75f, menu), in.fs, readme,
                        "AFTER: amount 0.75 with a HIGH crossover, 500 Hz — the voice only on the upper harmonics and the clank.");
    }

    if (std::FILE* f = std::fopen ((dir + "/README.txt").c_str(), "wb"))
    {
        std::fwrite (readme.data(), 1, readme.size(), f);
        std::fclose (f);
    }
    return allOk ? 0 : 1;
}
