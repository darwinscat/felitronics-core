// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics_envfilter_render — a listening tool, not a test. Runs a bass DI through envfilter::EnvelopeFilter at
// several positions of the one slider (and the band-pass mode and the down direction) and writes 24-bit WAVs.
//
//   felitronics_envfilter_render <outDir> [inDir]
//
// Every *.wav in inDir is rendered (mono or stereo, at its own rate). Without one, a bass line is synthesised:
// plucked strings (a modal model — each partial its own decay, the high ones faster, a little inharmonicity) from
// low B (30.87 Hz) up to G2 (98 Hz), fingered notes, ghost notes and a few hard slap attacks with their click,
// 25 s at 48 kHz — and it is written beside the renders as the "before".
#include <felitronics/envfilter/EnvelopeFilter.h>
#include <felitronics/envfilter/EnvelopeFilterLaw.h>
#include <felitronics/io/Wav.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iterator>
#include <string>
#include <vector>

namespace ef = felitronics::envfilter;
namespace fs = std::filesystem;

namespace
{
constexpr double kPi = 3.14159265358979323846;

struct Rng { std::uint32_t s = 0x1234567u; double next() { s = s * 1664525u + 1013904223u; return (double) s / 4294967296.0 * 2.0 - 1.0; } };

double midiHz (int m) { return 440.0 * std::pow (2.0, (m - 69) / 12.0); }

/** One plucked note, added into `out` from sample `at`: partials with their own decays, a 30 ms damp at note
    off, and the attack's own noise (a soft thump for a finger, a bright click for a slap). */
void pluck (std::vector<double>& out, double sr, int at, double hz, double lenS, double vel, bool slap, Rng& rng)
{
    const double beta = slap ? 0.12 : 0.22;            // where the string is struck — nearer the end is brighter
    const double tilt = slap ? 1.25 : 1.9;             // the partials' fall
    const double inharm = 6.0e-5;
    const int    total = (int) ((lenS + 0.06) * sr);
    const int    offAt = (int) (lenS * sr);
    for (int h = 1; h * hz < 7000.0; ++h)
    {
        const double fh  = h * hz * std::sqrt (1.0 + inharm * h * h);
        const double ah  = vel * std::fabs (std::sin (kPi * h * beta)) / std::pow ((double) h, tilt);
        const double tau = 2.2 / (1.0 + std::pow (fh / 350.0, 1.4));
        const double ph  = rng.next() * 0.3;
        for (int i = 0; i < total && at + i < (int) out.size(); ++i)
        {
            const double t = i / sr;
            double env = std::exp (-t / tau) * std::min (1.0, t / 0.0015);
            if (i > offAt) env *= std::exp (-(i - offAt) / (0.012 * sr));
            out[(std::size_t) (at + i)] += ah * env * std::sin (2.0 * kPi * fh * t + ph);
        }
    }
    // the attack's noise: one-pole coloured, a few ms
    const int nn = (int) ((slap ? 0.004 : 0.003) * sr);
    double lp = 0.0, prev = 0.0;
    const double a = slap ? 0.15 : 0.92;
    for (int i = 0; i < nn && at + i < (int) out.size(); ++i)
    {
        const double w = rng.next();
        lp = a * lp + (1.0 - a) * w;
        const double v = slap ? (lp - prev) : lp;      // a slap's click is the noise's bright edge
        prev = lp;
        out[(std::size_t) (at + i)] += (slap ? 0.5 : 0.25) * vel * v * std::exp (-(double) i / (0.25 * nn));
    }
}

std::vector<double> synthBass (double sr)
{
    // 96 BPM, ten bars of sixteenths. Each bar: {sixteenth, interval above the root, length in sixteenths, kind}
    // kind 0 finger, 1 ghost, 2 slap.
    struct Hit { int at, iv, len, kind; };
    const Hit bar[] = { { 0, 0, 3, 0 }, { 3, 0, 1, 1 }, { 4, 12, 2, 2 }, { 6, 0, 2, 0 }, { 8, 7, 2, 0 },
                        { 10, 10, 1, 0 }, { 11, 12, 1, 2 }, { 12, 0, 3, 0 }, { 15, 5, 1, 1 } };
    const int roots[] = { 28, 28, 33, 23, 28, 38, 31, 23, 28, 33 };   // E1 E1 A1 B0 E1 D2 G1 B0 E1 A1
    const double sixteenth = 60.0 / 96.0 / 4.0;
    std::vector<double> out ((std::size_t) (25.0 * sr), 0.0);
    Rng rng;
    for (int b = 0; b < (int) std::size (roots); ++b)
        for (const auto& h : bar)
        {
            int m = roots[b] + h.iv;
            while (m > 43) m -= 12;                       // nothing above G2
            const double vel = h.kind == 2 ? 1.0 : (h.kind == 1 ? 0.22 : 0.5 + 0.12 * rng.next());
            const int at = (int) (((b * 16 + h.at) * sixteenth + 0.25) * sr);
            pluck (out, sr, at, midiHz (m), h.len * sixteenth * 0.92, vel, h.kind == 2, rng);
        }
    // a DC-free DI at a -6 dBFS peak
    double hp = 0.0, prevIn = 0.0, peak = 0.0;
    const double c = std::exp (-2.0 * kPi * 20.0 / sr);
    for (auto& v : out) { const double y = c * (hp + v - prevIn); prevIn = v; hp = y; v = y; peak = std::max (peak, std::fabs (y)); }
    for (auto& v : out) v *= 0.5 / peak;
    return out;
}

struct Setting { const char* name; float amount; ef::Mode mode; ef::Direction dir; };

std::vector<std::vector<double>> render (const std::vector<std::vector<double>>& in, double sr, const Setting& st)
{
    ef::EnvelopeFilter f;
    const int nch = (int) in.size(), n = (int) in[0].size();
    if (! f.prepare (sr, 512, nch)) { std::fprintf (stderr, "prepare refused\n"); return {}; }
    ef::EnvelopeFilterParams menu;
    menu.mode = st.mode;
    menu.direction = st.dir;
    f.setParams (ef::paramsForAmount (st.amount, menu));
    std::vector<std::vector<float>> io ((std::size_t) nch, std::vector<float> ((std::size_t) n));
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < n; ++i) io[(std::size_t) c][(std::size_t) i] = (float) in[(std::size_t) c][(std::size_t) i];
    for (int at = 0; at < n; at += 512)
    {
        float* p[2] = { io[0].data() + at, nch > 1 ? io[1].data() + at : nullptr };
        if (! f.process (p, nch, std::min (512, n - at))) { std::fprintf (stderr, "process refused\n"); return {}; }
    }
    std::vector<std::vector<double>> out ((std::size_t) nch, std::vector<double> ((std::size_t) n));
    for (int c = 0; c < nch; ++c)
        for (int i = 0; i < n; ++i) out[(std::size_t) c][(std::size_t) i] = (double) io[(std::size_t) c][(std::size_t) i];
    return out;
}

double peakDb (const std::vector<std::vector<double>>& x)
{
    double p = 0.0;
    for (const auto& c : x) for (double v : c) p = std::max (p, std::fabs (v));
    return 20.0 * std::log10 (std::max (p, 1e-12));
}
/** THE DOOR DEMO: one held A1, clean, and the slider in band pass going 0 -> 0.75 at 1.0 s and back to 0 at 3.5 s,
    so the ladder's way in and out (its small pitch bend included) can be judged by ear on a sustained note. */
std::vector<double> renderDoors (const std::vector<double>& note, double sr)
{
    ef::EnvelopeFilter f;
    const int n = (int) note.size();
    if (! f.prepare (sr, 64, 1)) return {};
    ef::EnvelopeFilterParams menu;
    menu.mode = ef::Mode::bandPass;
    f.setParams (ef::paramsForAmount (0.0f, menu));
    std::vector<float> io ((std::size_t) n);
    for (int i = 0; i < n; ++i) io[(std::size_t) i] = (float) note[(std::size_t) i];
    const int on = (int) (1.0 * sr) / 64 * 64, off = (int) (3.5 * sr) / 64 * 64;
    for (int at = 0; at < n; at += 64)
    {
        if (at == on)  f.setParams (ef::paramsForAmount (0.75f, menu));
        if (at == off) f.setParams (ef::paramsForAmount (0.0f, menu));
        float* p[1] = { io.data() + at };
        if (! f.process (p, 1, std::min (64, n - at))) return {};
    }
    return std::vector<double> (io.begin(), io.end());
}

std::vector<double> heldNote (double sr)
{
    std::vector<double> out ((std::size_t) (5.0 * sr), 0.0);
    Rng rng;
    pluck (out, sr, (int) (0.1 * sr), midiHz (33), 4.8, 0.7, false, rng);   // A1, 55 Hz, held
    double peak = 0.0;
    for (double v : out) peak = std::max (peak, std::fabs (v));
    for (auto& v : out) v *= 0.4 / peak;
    return out;
}
} // namespace

int main (int argc, char** argv)
{
    if (argc < 2) { std::fprintf (stderr, "usage: %s <outDir> [inDir]\n", argv[0]); return 2; }
    const fs::path outDir = argv[1];

    struct Input { std::string stem; std::vector<std::vector<double>> ch; double sr; };
    std::vector<Input> inputs;
    if (argc > 2 && fs::is_directory (argv[2]))
        for (const auto& e : fs::directory_iterator (argv[2]))
            if (e.path().extension() == ".wav")
            {
                auto w = felitronics::io::readWav (e.path().string());
                if (! w.ok) { std::fprintf (stderr, "skip %s: %s\n", e.path().string().c_str(), w.error.c_str()); continue; }
                if (w.ch.size() > 2) w.ch.resize (2);
                inputs.push_back ({ e.path().stem().string(), w.ch, w.sr });
            }
    const bool synth = inputs.empty();
    if (synth)
    {
        inputs.push_back ({ "", { synthBass (48000.0) }, 48000.0 });
        const auto path = (outDir / "envfilter-dry-synth-bass.wav").string();
        if (! felitronics::io::writeWav (path, inputs[0].ch, 48000.0, 24, false)) { std::fprintf (stderr, "cannot write %s\n", path.c_str()); return 1; }
        std::printf ("%-40s peak %6.2f dBFS  (the input)\n", "envfilter-dry-synth-bass.wav", peakDb (inputs[0].ch));
    }

    const Setting settings[] = {
        { "amount-025",          0.25f, ef::Mode::lowPass,  ef::Direction::up },
        { "amount-050",          0.50f, ef::Mode::lowPass,  ef::Direction::up },
        { "amount-075",          0.75f, ef::Mode::lowPass,  ef::Direction::up },
        { "amount-100",          1.00f, ef::Mode::lowPass,  ef::Direction::up },
        { "bandpass-amount-075", 0.75f, ef::Mode::bandPass, ef::Direction::up },
        { "down-amount-075",     0.75f, ef::Mode::lowPass,  ef::Direction::down },
    };
    for (const auto& in : inputs)
        for (const auto& st : settings)
        {
            const auto out = render (in.ch, in.sr, st);
            if (out.empty()) return 1;
            const std::string name = "envfilter-" + (in.stem.empty() ? std::string() : in.stem + "-") + st.name + ".wav";
            if (! felitronics::io::writeWav ((outDir / name).string(), out, in.sr, 24, false)) { std::fprintf (stderr, "cannot write %s\n", name.c_str()); return 1; }
            std::printf ("%-40s peak %6.2f dBFS\n", name.c_str(), peakDb (out));
        }
    // The door demo, on its own held note.
    {
        const auto note = heldNote (48000.0);
        const auto out  = renderDoors (note, 48000.0);
        if (out.empty()) return 1;
        if (! felitronics::io::writeWav ((outDir / "envfilter-dry-held-note.wav").string(), { note }, 48000.0, 24, false)
            || ! felitronics::io::writeWav ((outDir / "envfilter-bandpass-doors.wav").string(), { out }, 48000.0, 24, false))
            return 1;
        std::printf ("%-40s peak %6.2f dBFS  (held A1; band pass, slider 0 -> 0.75 at 1.0 s -> 0 at 3.5 s)\n",
                     "envfilter-bandpass-doors.wav", peakDb ({ out }));
    }
    return 0;
}
