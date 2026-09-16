// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// The streaming surface of the wasm ABI (P120) — the fc_stream_* handles over fcore::StreamProbe — compiled
// natively, so ASan and UBSan see it.
//
// WHAT IT PROVES. A programme fed through the handles in pieces of 1 sample, of 4096, and of random sizes (zero
// included) ends with the same integrated, momentary and short-term loudness, bit for bit, and the same runs,
// field for field, as the SAME instruments — analysis::DeterministicLoudnessMeter and analysis::ClipDetector —
// handed the whole buffer in one call and read here directly, not through the class under test. The runs of the
// streams are the ones polled BETWEEN the pieces through `from` and concatenated, so the list a page assembles
// while the audio is still arriving is what is compared.
//
// The three streams are open AT ONCE and their pieces interleaved, at each rate: a handle that shared state with
// another would move a reading of both away from the oracle.
//
// 88200 IS IN THE RATES ON PURPOSE: it is where the system-libm meter's block energies differ between native and
// wasm (32 of 37 blocks of a 4 s generated fixture), and the reason this surface measures with the deterministic one.

#include <felitronics_test.h>

#include <felitronics/analysis/ClipDetector.h>
#include <felitronics/analysis/LoudnessMeter.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace felitronics;
using Planes = std::vector<std::vector<float>>;

// The C ABI, written out rather than included: fc_probe.cpp has no header, and a signature that drifts fails to
// link, which is the diagnosis.
extern "C" {
    std::uint32_t fc_stream_create (double sampleRate, std::uint32_t channels);
    int           fc_stream_process (std::uint32_t h, const float* planar, std::uint32_t frames);
    std::uint32_t fc_stream_loudness_fields (void);
    int           fc_stream_loudness (std::uint32_t h, double* out);
    double        fc_stream_clips_count (std::uint32_t h);
    std::uint32_t fc_stream_clips_stride (void);
    std::uint32_t fc_stream_clips (std::uint32_t h, std::uint32_t from, double* out, std::uint32_t cap);
    int           fc_stream_finish (std::uint32_t h);
    int           fc_stream_destroy (std::uint32_t h);
}

static std::uint64_t b64 (double d) noexcept { return std::bit_cast<std::uint64_t> (d); }

// Clamped, on the 16-bit lattice, with a loud passage, a passage 40 dB down (under the relative gate), digital
// silence, non-finite holes and a clamped tail whose runs exist only after finish(). Passages are fractions of the
// length, so every rate gets all of them.
static Planes programme (long long n, int channels, unsigned seed)
{
    Planes p ((std::size_t) channels, std::vector<float> ((std::size_t) n, 0.0f));
    std::mt19937 rng (seed);
    std::uniform_int_distribution<int> N (-40, 40);
    const int ceilCode = 29196;
    for (int c = 0; c < channels; ++c)
    {
        const long long period = 96 + 16 * c;
        for (long long i = 0; i < n; ++i)
        {
            const long long k = (i + 7 * c) % (2 * period), t = k % period;
            const double sign = k < period ? 1.0 : -1.0;
            const double u = (double) (2 * t - period) / (double) period;
            double code = sign * 1.25 * (double) ceilCode * (1.0 - u * u) + (double) N (rng);
            const long long part = i * 10 / n;                            // tenths of the programme
            if (part == 4 || part == 5) code *= 0.01;                     // -40 dB: a block the relative gate drops
            if (part == 7) code = 0.0;                                    // digital silence
            const double clamped = std::max (-(double) ceilCode, std::min ((double) ceilCode, std::round (code)));
            float v = (float) (clamped / 32768.0);
            if (part == 2 && (i % 977) == 13 * (c + 1)) v = std::numeric_limits<float>::quiet_NaN();
            p[(std::size_t) c][(std::size_t) i] = v;
        }
    }
    return p;
}

struct Oracle
{
    double momentary = 0, shortTerm = 0, integrated = 0;
    int dropped = 0;
    std::int64_t samples = 0, runCount = 0;
    std::vector<analysis::ClipRun> runs;
};

static Oracle wholeBuffer (const Planes& p, double sr)
{
    Oracle o;
    const int nc = (int) p.size();
    const int n = (int) p[0].size();
    const float* view[core::kMaxChannels] {};
    for (int c = 0; c < nc; ++c) view[c] = p[(std::size_t) c].data();

    analysis::DeterministicLoudnessMeter m;
    test::ok (m.prepare (sr, nc, 3600.0) && m.process (view, nc, n), "the oracle meter takes the whole buffer");
    o.momentary = m.momentaryLufs();
    o.shortTerm = m.shortTermLufs();
    o.integrated = m.integratedLufs();
    o.dropped = m.droppedBlocks();

    analysis::ClipDetector d;
    d.setParams (analysis::ClipDetectorParams {});
    test::ok (d.prepare (sr, 0, nc) && d.process (view, nc, n), "the oracle detector takes the whole buffer");
    d.finish();
    o.samples = d.samplesProcessed();
    o.runCount = d.runCount();
    for (std::int64_t i = 0; i < d.storedRunCount(); ++i) o.runs.push_back (d.run (i));
    return o;
}

struct Stream
{
    std::string name;
    std::vector<long long> pieces;
    std::size_t next = 0;
    long long at = 0;
    std::uint32_t h = 0;
    std::vector<double> runs;          // stride doubles a run, in the order they were polled
    bool processOk = true, loudnessOk = true;
};

// Polls the runs decided so far through a buffer of three runs, so `from` advances in short steps.
static void drain (Stream& s)
{
    const std::uint32_t stride = fc_stream_clips_stride();
    double buf[3 * 6] {};
    for (;;)
    {
        const std::uint32_t got = fc_stream_clips (s.h, (std::uint32_t) (s.runs.size() / stride), buf, 3u * stride);
        if (got == 0) break;
        s.runs.insert (s.runs.end(), buf, buf + (std::size_t) got * stride);
    }
}

static void piecewiseEqualsWhole (double sr)
{
    const std::string at = " @ " + std::to_string ((int) sr);
    test::group ("pieces of 1, 4096 and random sizes, three handles at once, against one call" + at);
    const int nc = 2;
    const long long n = (long long) (4.0 * sr);
    const Planes p = programme (n, nc, 20260916u);
    const Oracle o = wholeBuffer (p, sr);
    // A comparison of nothing proves nothing: the programme must be loud, clipped, and hold runs finish() decides.
    test::ok (o.integrated > -70.0, "the programme has an integrated loudness to compare" + at);
    test::ok (o.runs.size() > 20 && o.runCount == (std::int64_t) o.runs.size(), "and stored runs to compare (" + std::to_string (o.runs.size()) + ")" + at);

    std::vector<Stream> streams (3);
    streams[0].name = "1";
    streams[0].pieces.assign ((std::size_t) n, 1);
    streams[1].name = "4096";
    for (long long a = 0; a < n; a += 4096) streams[1].pieces.push_back (std::min<long long> (4096, n - a));
    streams[2].name = "random";
    std::mt19937 rng (7u);
    std::uniform_int_distribution<long long> U (1, 9000);
    for (long long a = 0; a < n; )
    {
        const long long k = std::min (n - a, U (rng));
        streams[2].pieces.push_back (k);
        if (streams[2].pieces.size() % 5 == 0) streams[2].pieces.push_back (0);   // a zero-length piece every fifth call
        a += k;
    }

    for (auto& s : streams)
    {
        s.h = fc_stream_create (sr, (std::uint32_t) nc);
        test::ok (s.h != 0, "fc_stream_create gives stream '" + s.name + "' a handle" + at);
    }
    test::ok (streams[0].h != streams[1].h && streams[1].h != streams[2].h && streams[0].h != streams[2].h,
              "three live streams, three handles" + at);

    const std::uint32_t fields = fc_stream_loudness_fields();
    test::ok (fields == 5, "five loudness fields");
    std::vector<float> planar;
    for (;;)
    {
        // The stream furthest behind goes next, so the pieces of all three interleave on one sample clock.
        Stream* s = nullptr;
        for (auto& c : streams)
            if (c.next < c.pieces.size() && (s == nullptr || c.at < s->at)) s = &c;
        if (s == nullptr) break;
        const long long k = s->pieces[s->next++];
        planar.assign ((std::size_t) std::max<long long> (1, k * nc), 0.0f);
        for (int c = 0; c < nc; ++c)
            std::copy_n (p[(std::size_t) c].data() + s->at, (std::size_t) k, planar.data() + (std::size_t) (c * k));
        s->processOk = s->processOk && fc_stream_process (s->h, planar.data(), (std::uint32_t) k) == 1;
        double lo[5] {};
        s->loudnessOk = s->loudnessOk && fc_stream_loudness (s->h, lo) == 1 && b64 (lo[3]) == b64 ((double) (s->at + k));
        s->at += k;
        drain (*s);
    }

    const std::uint32_t stride = fc_stream_clips_stride();
    for (auto& s : streams)
    {
        const std::string w = " (" + s.name + at + ")";
        test::ok (s.processOk, "every piece is accepted" + w);
        test::ok (s.loudnessOk, "and the loudness is readable after every piece, at the right sample count" + w);
        const std::size_t beforeFinish = s.runs.size() / stride;
        test::ok (fc_stream_finish (s.h) == 1, "finish" + w);
        drain (s);
        test::ok (s.runs.size() / stride > beforeFinish, "and finish() decides runs no piece could (" + std::to_string (beforeFinish) + " before it)" + w);

        double lo[5] {};
        test::ok (fc_stream_loudness (s.h, lo) == 1, "the loudness reads after finish" + w);
        test::ok (b64 (lo[2]) == b64 (o.integrated), "integrated LUFS, bit for bit" + w);
        test::ok (b64 (lo[0]) == b64 (o.momentary) && b64 (lo[1]) == b64 (o.shortTerm),
                  "momentary and short-term, bit for bit" + w);
        test::ok (b64 (lo[3]) == b64 ((double) o.samples) && b64 (lo[4]) == b64 ((double) o.dropped), "samples and dropped blocks" + w);
        test::ok (b64 (fc_stream_clips_count (s.h)) == b64 ((double) o.runCount), "the run count" + w);

        bool same = s.runs.size() == o.runs.size() * stride;
        for (std::size_t i = 0; same && i < o.runs.size(); ++i)
        {
            const analysis::ClipRun& r = o.runs[i];
            const double* f = s.runs.data() + i * stride;
            same = b64 (f[0]) == b64 ((double) r.start) && b64 (f[1]) == b64 ((double) r.length) && b64 (f[2]) == b64 (r.level)
                && b64 (f[3]) == b64 ((double) r.channel) && b64 (f[4]) == b64 ((double) r.sign)
                && b64 (f[5]) == b64 ((double) (int) r.evidence);
        }
        test::ok (same, "the runs polled between the pieces are the runs of the whole buffer, field for field ("
                        + std::to_string (s.runs.size() / stride) + " of " + std::to_string (o.runs.size()) + ")" + w);
    }
    for (auto& s : streams) test::ok (fc_stream_destroy (s.h) == 1, "destroy '" + s.name + "'" + at);
}

int main()
{
    std::printf ("felitronics_stream_abi_tests\n");
    for (double sr : { 44100.0, 48000.0, 88200.0 }) piecewiseEqualsWhole (sr);
    return test::report();
}
