// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// fcore_octave_render — offline listening renders for octave::Octaver.
//
//   fcore_octave_render <out-dir> [in.wav ...]
//
// Each input (a bass DI; channels are processed as the stage gets them, 1 or 2) is rendered through the one-slider
// law at amount 0.25 / 0.5 / 0.75 / 1 (the default menu: −1 up to 1.25, −2 off), through the −2 octave alone under
// the dry, and through each octave alone with no dry — 24-bit WAVs named `octave-<input>-<setting>.wav`, with the input itself as `-before`. With no
// input it synthesises a bass line first (modal plucked strings, fingered and slapped, low B 30.9 Hz up to G2
// 98 Hz, a long decay, a gap over a noise floor) and renders that, named `octave-<setting>.wav` — and, for the SUB HP
// question, four open-B plucks at the slider's top under SUB HP 25 Hz and 10 Hz (`octave-openB-*.wav`).

#include <felitronics/io/Wav.h>
#include <felitronics/octave/OctaveLaw.h>
#include <felitronics/octave/Octaver.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
namespace oc = felitronics::octave;
using Planes = std::vector<std::vector<double>>;

constexpr double kPi = 3.14159265358979323846;

struct Note { double at, len, hz; bool slap; double vel; };

/** One plucked note into `out` at `at` seconds: modal, plucked at `beta`, a velocity pickup near the bridge,
    upper modes damped faster, muted with a 25 ms fade at the note's end. A slap is plucked near the end of the
    string, brighter, louder, with a 3 ms click of the string hitting the frets. */
void addNote (std::vector<double>& out, double fs, const Note& n, std::uint32_t& seed)
{
    const double beta = n.slap ? 0.07 : 0.24, gamma = 0.12, bright = n.slap ? 0.985 : 0.93;
    const double tail = 0.025;
    const auto start = (std::size_t) std::lround (n.at * fs);
    const auto len   = (std::size_t) std::lround ((n.len + tail) * fs);
    struct Mode { double a, w, d; };
    std::vector<Mode> modes;
    for (int k = 1; k <= 120; ++k)
    {
        const double f = n.hz * k * (1.0 + 1.0e-4 * k * k);
        if (f > 0.45 * fs || f > 9000.0) break;
        const double a = std::sin (k * kPi * beta) * std::sin (k * kPi * gamma) / k * std::pow (bright, k - 1);
        const double d = (n.slap ? 0.5 + 0.03 * k * k * n.hz / 40.0 : 0.7 + 0.08 * k * k * n.hz / 40.0);
        modes.push_back ({ a, 2.0 * kPi * f / fs, d / fs });
    }
    const double muteAt = n.len * fs, muteLen = tail * fs;
    for (std::size_t i = 0; i < len && start + i < out.size(); ++i)
    {
        double s = 0.0;
        for (const auto& m : modes) s += m.a * std::exp (-m.d * (double) i) * std::sin (m.w * (double) i);
        if (n.slap && i < (std::size_t) (0.003 * fs))
        {
            seed = seed * 1664525u + 1013904223u;
            s += 0.15 * ((double) (std::int32_t) seed / 2147483648.0) * (1.0 - (double) i / (0.003 * fs));
        }
        const double mute = (double) i < muteAt ? 1.0 : std::max (0.0, 1.0 - ((double) i - muteAt) / muteLen);
        out[start + i] += n.vel * mute * s;
    }
}

/** ~26 s at 48 kHz: fingered lines across the range, a slap passage, a low B left to ring for four seconds, a gap
    over a -84 dBFS floor, and a last low E. Peak -10 dBFS. */
std::vector<double> synthLine (double fs)
{
    const double e = 0.3;   // an eighth at 100 BPM
    const double B0 = 30.868, Cs1 = 34.648, D1 = 36.708, E1 = 41.203, Fs1 = 46.249, G1 = 48.999, A1 = 55.0,
                 B1 = 61.735, Cs2 = 69.296, D2 = 73.416, E2 = 82.407, Fs2 = 92.499, G2 = 97.999;
    std::vector<Note> notes;
    double t = 0.2;
    const auto play = [&] (double hz, double eighths, bool slap = false, double vel = 1.0, double gate = 0.92)
    {
        notes.push_back ({ t, eighths * e * gate, hz, slap, vel });
        t += eighths * e;
    };
    // Fingered, low: a root-fifth figure on E and on B.
    for (int bar = 0; bar < 2; ++bar)
    {
        play (E1, 2); play (E1, 1, false, 0.8); play (B1, 1, false, 0.9);
        play (E2, 1, false, 0.85); play (D2, 1, false, 0.8); play (B1, 2);
    }
    play (B0, 3); play (B0, 1, false, 0.8); play (Fs1, 2); play (A1, 2);
    // A walk up to G2 and back.
    for (double hz : { G1, A1, B1, Cs2, D2, E2, Fs2, G2, Fs2, E2, D2, B1, A1, G1, Fs1, E1 })
        play (hz, 1, false, 0.85, 0.97);
    // Slap: thumb on the roots, pops an octave up (up to G2).
    for (int bar = 0; bar < 2; ++bar)
    {
        play (E1, 1, true, 1.3); play (E2, 1, true, 1.1); play (E1, 1, true, 1.2, 0.5); play (G1, 1, true, 1.2);
        play (A1, 1, true, 1.3); play (G2, 1, true, 1.1); play (D1, 1, true, 1.2); play (Cs1, 1, true, 1.2);
    }
    // A low B left to ring, then a gap.
    play (B0, 14, false, 1.1, 0.98);
    t += 1.5;
    play (E1, 6, false, 1.0, 0.95);

    std::vector<double> out ((std::size_t) std::lround ((t + 0.5) * fs), 0.0);
    std::uint32_t seed = 0x0C7A5EEDu;
    for (const auto& n : notes) addNote (out, fs, n, seed);

    double peak = 0.0;
    for (double v : out) peak = std::max (peak, std::fabs (v));
    const double floorAmp = std::pow (10.0, -84.0 / 20.0);
    for (double& v : out)
    {
        seed = seed * 1664525u + 1013904223u;
        v = v / peak * std::pow (10.0, -10.0 / 20.0) + floorAmp * ((double) (std::int32_t) seed / 2147483648.0);
    }
    return out;
}

/** Four open-B plucks (30.9 Hz), fingered, a second and a half each: what SUB HP does to the B string's −1. */
std::vector<double> openBLine (double fs)
{
    std::vector<double> out ((std::size_t) std::lround (6.5 * fs), 0.0);
    std::uint32_t seed = 0x0B0B0B0Bu;
    for (int k = 0; k < 4; ++k)
        addNote (out, fs, { 0.2 + 1.5 * k, 1.35, 30.868, false, 1.0 }, seed);
    double peak = 0.0;
    for (double v : out) peak = std::max (peak, std::fabs (v));
    for (double& v : out) v = v / peak * std::pow (10.0, -10.0 / 20.0);
    return out;
}

Planes run (const Planes& in, double fs, const oc::OctaverParams& p)
{
    const int nch = (int) std::min<std::size_t> (in.size(), 2), block = 256;
    oc::Octaver o;
    if (! o.prepare (fs, block, nch)) return {};
    o.setParams (p);
    std::vector<std::vector<float>> buf ((std::size_t) nch);
    for (int c = 0; c < nch; ++c)
        for (double v : in[(std::size_t) c]) buf[(std::size_t) c].push_back ((float) v);
    const std::size_t n = buf[0].size();
    for (std::size_t at = 0; at < n; at += (std::size_t) block)
    {
        float* io[2] = { buf[0].data() + at, nch > 1 ? buf[1].data() + at : nullptr };
        if (! o.process (io, nch, (int) std::min<std::size_t> ((std::size_t) block, n - at))) return {};
    }
    Planes out ((std::size_t) nch);
    for (int c = 0; c < nch; ++c) out[(std::size_t) c].assign (buf[(std::size_t) c].begin(), buf[(std::size_t) c].end());
    return out;
}

/** The sub ALONE: the stage has no dry level, so it is read from a stereo call whose right channel is silent —
    {2·mean, 0} tracks the mean, and the right channel comes back as the sub. Mono out. */
Planes alone (const Planes& in, double fs, const oc::OctaverParams& p)
{
    Planes two (2, std::vector<double> (in[0].size(), 0.0));
    for (const auto& c : in)
        for (std::size_t i = 0; i < c.size(); ++i) two[0][i] += 2.0 * c[i] / (double) in.size();
    const Planes out = run (two, fs, p);
    return out.empty() ? out : Planes { out[1] };
}

bool write (const std::string& path, const Planes& x, double fs)
{
    double peak = 0.0;
    for (const auto& c : x) for (double v : c) peak = std::max (peak, std::fabs (v));
    const bool okw = ! x.empty() && felitronics::io::writeWav (path, x, fs, 24, false);
    std::printf ("%s  %s  peak %.1f dBFS\n", okw ? "wrote" : "FAILED", path.c_str(), 20.0 * std::log10 (peak + 1e-12));
    return okw;
}

std::string stem (const std::string& path)
{
    const auto slash = path.find_last_of ("/\\");
    std::string s = slash == std::string::npos ? path : path.substr (slash + 1);
    const auto dot = s.find_last_of ('.');
    return dot == std::string::npos ? s : s.substr (0, dot);
}
} // namespace

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf (stderr, "usage: fcore_octave_render <out-dir> [in.wav ...]\n");
        return 2;
    }
    const std::string dir = argv[1];

    struct Input { std::string prefix; Planes x; double fs; };
    std::vector<Input> inputs;
    for (int k = 2; k < argc; ++k)
    {
        auto w = felitronics::io::readWav (argv[k]);
        if (! w.ok) { std::fprintf (stderr, "%s: %s\n", argv[k], w.error.c_str()); return 1; }
        if (w.ch.size() > 2) w.ch.resize (2);
        inputs.push_back ({ dir + "/octave-" + stem (argv[k]) + "-", w.ch, w.sr });
    }
    const bool synthesised = inputs.empty();
    if (synthesised)
        inputs.push_back ({ dir + "/octave-", Planes { synthLine (48000.0) }, 48000.0 });

    bool all = true;
    for (const auto& in : inputs)
    {
        all = write (in.prefix + "before.wav", in.x, in.fs) && all;
        for (float a : { 0.25f, 0.5f, 0.75f, 1.0f })
        {
            char name[32];
            std::snprintf (name, sizeof name, "amount-%.2f.wav", (double) a);
            all = write (in.prefix + name, run (in.x, in.fs, oc::paramsForAmount (a)), in.fs) && all;
        }
        oc::OctaverParams p;
        p.sub1 = 0.0f; p.sub2 = 1.0f;
        all = write (in.prefix + "sub2-under-dry.wav", run (in.x, in.fs, p), in.fs) && all;
        p.sub1 = 1.0f; p.sub2 = 0.0f;
        all = write (in.prefix + "sub1-alone.wav", alone (in.x, in.fs, p), in.fs) && all;
        p.sub1 = 0.0f; p.sub2 = 1.0f;
        all = write (in.prefix + "sub2-alone.wav", alone (in.x, in.fs, p), in.fs) && all;
    }
    if (synthesised)
    {
        // The B string's −1 under SUB HP 25 Hz (the default) and 10 Hz (the menu's lowest), at the slider's top.
        const Planes b { openBLine (48000.0) };
        all = write (dir + "/octave-openB-before.wav", b, 48000.0) && all;
        for (float hp : { 25.0f, 10.0f })
        {
            oc::OctaverParams menu;
            menu.subHpHz = hp;
            char name[64];
            std::snprintf (name, sizeof name, "/octave-openB-amount-1-subhp-%.0f.wav", (double) hp);
            all = write (dir + name, run (b, 48000.0, oc::paramsForAmount (1.0f, menu)), 48000.0) && all;
        }
    }
    return all ? 0 : 1;
}
