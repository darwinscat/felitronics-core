// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The P77 ABI surface under adversarial arguments — natively, so it runs under ASan/UBSan on every row.
//
// WHY THIS EXISTS. The parity harnesses exercise these getters only with capacities the module itself
// handed out, which is exactly the case that cannot fail. Everything a real caller can do wrong — a null
// pointer, a capacity of zero, a capacity one short of a row, a buffer that is not a multiple of the
// stride, a getter read before any run or after a refused one — was checked by hand once and would not
// have been checked again. A safety contract with no test is a comment.
//
// Every buffer here is surrounded by canaries and every call is followed by a check that they are intact:
// a getter that writes one element past its capacity is the defect this suite exists to catch, and it is
// invisible to a diff of two successful runs.

#include <felitronics_test.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

extern "C"
{
    int           fc_probe_report_run      (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_report_count_rows (void);
    std::uint32_t fc_probe_report_value_rows (void);
    std::uint32_t fc_probe_report_counts   (double*, std::uint32_t);
    std::uint32_t fc_probe_report_values   (double*, std::uint32_t);
    std::uint32_t fc_probe_report_names    (char*, std::uint32_t);

    int           fc_probe_bursts_run      (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_bursts_scalars  (double*, std::uint32_t);
    std::uint32_t fc_probe_bursts_chan     (double*, std::uint32_t);
    std::uint32_t fc_probe_bursts_events   (double*, std::uint32_t);
    std::uint32_t fc_probe_bursts_ioi      (double*, std::uint32_t);
    std::uint32_t fc_probe_bursts_lag      (double*, std::uint32_t);

    int           fc_probe_hum_run         (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_hum_scalars     (double*, std::uint32_t);
    std::uint32_t fc_probe_hum_chan        (double*, std::uint32_t);
    std::uint32_t fc_probe_hum_cand        (double*, std::uint32_t);
    std::uint32_t fc_probe_hum_harm        (double*, std::uint32_t);
    std::uint32_t fc_probe_hum_stretch     (double*, std::uint32_t);

    int           fc_probe_forensics_run   (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_forensics_scalars (double*, std::uint32_t);
    std::uint32_t fc_probe_forensics_wall  (double*, std::uint32_t);
    std::uint32_t fc_probe_forensics_grid  (double*, std::uint32_t);
    std::uint32_t fc_probe_forensics_khist (double*, std::uint32_t);

    int           fc_probe_lowend_run      (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_lowend_scalars  (double*, std::uint32_t);
    std::uint32_t fc_probe_lowend_hist     (double*, std::uint32_t);
    std::uint32_t fc_probe_lowend_series   (double*, std::uint32_t);
    std::uint32_t fc_probe_lowend_bands    (double*, std::uint32_t);
    std::uint32_t fc_probe_lowend_note_name (char*, std::uint32_t);
}

using felitronics::test::ok;

namespace
{
using Getter = std::uint32_t (*) (double*, std::uint32_t);

constexpr double kCanary = 1234567.875;   // exactly representable, so a survivor is a survivor

// Calls `g` with `cap`, into a buffer of exactly `cap` doubles fenced by canaries, and reports whether
// the fence held and the return is within the capacity it was given.
bool fenced (Getter g, std::uint32_t cap, const std::string& name, std::uint32_t stride)
{
    std::vector<double> buf ((std::size_t) cap + 4, kCanary);
    const std::uint32_t rows = g (buf.data() + 2, cap);
    bool fence = true;
    for (int i = 0; i < 2; ++i)
    {
        if (buf[(std::size_t) i] != kCanary) fence = false;
        if (buf[buf.size() - 1 - (std::size_t) i] != kCanary) fence = false;
    }
    ok (fence, name + ": the fence around a " + std::to_string (cap) + "-double buffer is intact");
    ok (stride == 0 || (std::uint64_t) rows * stride <= (std::uint64_t) cap,
        name + ": returned " + std::to_string (rows) + " rows for a capacity of " + std::to_string (cap));
    return fence;
}

void hammer (Getter g, const std::string& name, std::uint32_t stride)
{
    ok (g (nullptr, 0) == 0, name + ": a null pointer with zero capacity writes nothing");
    ok (g (nullptr, 1000) == 0, name + ": a null pointer with a generous capacity writes nothing");
    for (std::uint32_t cap : { 0u, 1u, 2u, 3u, 7u, 64u })
        fenced (g, cap, name, stride);
    // a capacity that is NOT a whole number of rows: the getter must floor it, never round up
    if (stride > 1) fenced (g, stride - 1, name + " [one short of a row]", stride);
    if (stride > 1) fenced (g, stride * 3 - 1, name + " [three rows less one]", stride);

    // AND THE CAPACITY MUST ACTUALLY BIND, or nothing above tested the truncation path at all. Asking for
    // one row and getting one row does NOT prove that: it is what a list of length one returns anyway.
    // Learn the full length first, then ask for one row less — only then is the capacity the constraint.
    // Measured: without this, a `room + 1` mutant in bursts_events left the suite green.
    {
        std::vector<double> big (4096, kCanary);
        const std::uint32_t full = g (big.data(), (std::uint32_t) (big.size() / stride) * stride);
        if (full >= 2)
        {
            const std::uint32_t want = full - 1;
            std::vector<double> buf ((std::size_t) want * stride + 4, kCanary);
            const std::uint32_t got = g (buf.data() + 2, want * stride);
            bool fence = true;
            for (int i = 0; i < 2; ++i)
                if (buf[(std::size_t) i] != kCanary || buf[buf.size() - 1 - (std::size_t) i] != kCanary) fence = false;
            ok (got == want, name + ": a capacity of " + std::to_string (want) + " rows against a list of "
                             + std::to_string (full) + " returns exactly " + std::to_string (got));
            ok (fence, name + ": and does not write past it");
        }
        else
        {
            ok (true, name + ": only " + std::to_string (full) + " row(s) on this fixture — THE TRUNCATION "
                      "PATH IS NOT EXERCISED HERE, said out loud rather than passed in silence");
        }
    }
}

// White noise with PERIODIC BURSTS in the band BandBursts watches. Without them the event list holds
// fewer than two rows, its truncation path is never exercised, and a getter mutated to write one row past
// its capacity passes the whole suite — measured, before this fixture existed.
std::vector<float> burstyFixture (int frames, int channels)
{
    std::vector<float> v ((std::size_t) frames * (std::size_t) channels, 0.0f);
    std::uint64_t st = 0x9E3779B97F4A7C15ull;
    const int period = frames / 6, width = frames / 30;   // ~167 ms apart, ~33 ms wide: several hops each
    for (int i = 0; i < frames; ++i)
    {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        const double n = (double) (st >> 11) / 9007199254740992.0 * 2.0 - 1.0;
        const bool inBurst = (i % (period > 0 ? period : 1)) < (width > 0 ? width : 1);
        const double a = inBurst ? 0.7 : 0.05;
        for (int c = 0; c < channels; ++c) v[(std::size_t) c * (std::size_t) frames + (std::size_t) i] = (float) (a * n);
    }
    return v;
}

std::vector<float> fixture (int frames, int channels)
{
    std::vector<float> v ((std::size_t) frames * (std::size_t) channels, 0.0f);
    std::uint64_t st = 0x2545F4914F6CDD1Dull;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        v[i] = (float) ((double) (st >> 11) / 9007199254740992.0 * 1.6 - 0.8);
    }
    return v;
}
} // namespace

int main()
{
    const int frames = 48000, ch = 2;
    std::vector<float> planar = fixture (frames, ch);

    // ---------- before any run, every getter answers zero and touches nothing ----------
    {
        const std::pair<Getter, const char*> before[] = {
            { fc_probe_report_counts, "report_counts" }, { fc_probe_bursts_scalars, "bursts_scalars" },
            { fc_probe_hum_chan, "hum_chan" }, { fc_probe_forensics_wall, "forensics_wall" },
            { fc_probe_lowend_series, "lowend_series" },
        };
        for (auto [g, n] : before)
        {
            std::vector<double> buf (64, kCanary);
            ok (g (buf.data(), 60) == 0, std::string (n) + ": answers zero before any run");
            bool clean = true;
            for (double d : buf) if (d != kCanary) clean = false;
            ok (clean, std::string (n) + ": and writes nothing before any run");
        }
    }

    // ---------- a REFUSED run leaves the getters silent ----------
    {
        ok (fc_probe_report_run (planar.data(), (std::uint32_t) frames, 99u, 48000.0) == 0,
            "report_run: refuses 99 channels");
        std::vector<double> buf (64, kCanary);
        ok (fc_probe_report_counts (buf.data(), 60) == 0, "report_counts: silent after a refused run");
        ok (fc_probe_report_run (planar.data(), (std::uint32_t) frames, (std::uint32_t) ch, 0.0) == 0,
            "report_run: refuses a zero sample rate");
        ok (fc_probe_report_run (nullptr, (std::uint32_t) frames, (std::uint32_t) ch, 48000.0) == 0,
            "report_run: refuses a null input with a non-zero frame count");
    }

    // ---------- the real thing, then every getter hammered ----------
    struct Mode { const char* name; int (*run) (const float*, std::uint32_t, std::uint32_t, double); };
    const Mode modes[] = {
        { "report",    fc_probe_report_run },    { "bursts", fc_probe_bursts_run },
        { "hum",       fc_probe_hum_run },       { "forensics", fc_probe_forensics_run },
        { "lowend",    fc_probe_lowend_run },
    };
    for (const Mode& m : modes)
        ok (m.run (planar.data(), (std::uint32_t) frames, (std::uint32_t) ch, 48000.0) == 1,
            std::string (m.name) + "_run: accepts a real programme");
    // bursts gets the bursty one, so its event list is long enough for a capacity to bind against.
    {
        // Five seconds, not one: the detector's baseline is a long moving window, and a one-second
        // programme never establishes one, so a shorter fixture yields no events at all.
        const int burstFrames = 5 * 48000;
        const std::vector<float> bursty = burstyFixture (burstFrames, ch);
        ok (fc_probe_bursts_run (bursty.data(), (std::uint32_t) burstFrames, (std::uint32_t) ch, 48000.0) == 1,
            "bursts_run: accepts the bursty programme");
        std::vector<double> probe (4096, 0.0);
        const std::uint32_t events = fc_probe_bursts_events (probe.data(), 4092);
        ok (events >= 2, "bursts: the bursty fixture yields " + std::to_string (events)
                         + " events — enough for a capacity to bind against");
    }

    hammer (fc_probe_report_counts, "report_counts", 2);
    hammer (fc_probe_report_values, "report_values", 4);
    hammer (fc_probe_bursts_scalars, "bursts_scalars", 0);
    hammer (fc_probe_bursts_chan,   "bursts_chan",   4);
    hammer (fc_probe_bursts_events, "bursts_events", 12);
    hammer (fc_probe_bursts_ioi,    "bursts_ioi",    2);
    hammer (fc_probe_bursts_lag,    "bursts_lag",    2);
    hammer (fc_probe_hum_scalars,   "hum_scalars",   0);
    hammer (fc_probe_hum_chan,      "hum_chan",      21);
    hammer (fc_probe_hum_cand,      "hum_cand",      19);
    hammer (fc_probe_hum_harm,      "hum_harm",      8);
    hammer (fc_probe_hum_stretch,   "hum_stretch",   5);
    hammer (fc_probe_forensics_scalars, "forensics_scalars", 0);
    hammer (fc_probe_forensics_wall,    "forensics_wall",    39);
    hammer (fc_probe_forensics_grid,    "forensics_grid",    22);
    hammer (fc_probe_forensics_khist,   "forensics_khist",   0);
    hammer (fc_probe_lowend_scalars, "lowend_scalars", 0);
    hammer (fc_probe_lowend_hist,    "lowend_hist",    0);
    hammer (fc_probe_lowend_series,  "lowend_series",  6);
    hammer (fc_probe_lowend_bands,   "lowend_bands",   11);

    // ---------- the two char buffers, whose capacity is in BYTES ----------
    {
        const std::uint32_t need = fc_probe_report_names (nullptr, 0);
        ok (need > 0, "report_names: a zero capacity answers the size required (" + std::to_string (need) + " bytes)");
        std::vector<char> small ((std::size_t) need, '\0');
        ok (fc_probe_report_names (small.data(), need - 1) == 0,
            "report_names: a capacity one byte short writes NOTHING rather than a truncated blob");
        std::vector<char> exact ((std::size_t) need + 4, '\x7f');
        ok (fc_probe_report_names (exact.data() + 2, need) == need, "report_names: an exact capacity is filled");
        ok (exact[0] == '\x7f' && exact[1] == '\x7f' && exact[exact.size() - 1] == '\x7f',
            "report_names: and the fence around it is intact");
        std::uint32_t seps = 0;
        for (std::uint32_t i = 0; i < need; ++i) if (exact[(std::size_t) i + 2] == '\0') ++seps;
        ok (seps == fc_probe_report_count_rows() + fc_probe_report_value_rows(),
            "report_names: one NUL-terminated name per row, " + std::to_string (seps) + " of them");
    }
    {
        const std::uint32_t need = fc_probe_lowend_note_name (nullptr, 0);
        if (need > 0)
        {
            std::vector<char> exact ((std::size_t) need + 4, '\x7f');
            ok (fc_probe_lowend_note_name (exact.data() + 2, need) == need, "lowend_note_name: an exact capacity is filled");
            ok (exact[0] == '\x7f' && exact[exact.size() - 1] == '\x7f', "lowend_note_name: the fence is intact");
            ok (fc_probe_lowend_note_name (exact.data() + 2, need - 1) == 0,
                "lowend_note_name: a short capacity writes nothing");
        }
        else ok (true, "lowend_note_name: this fixture has no valid note, so there is no name — said out loud");
    }

    return felitronics::test::report();
}
