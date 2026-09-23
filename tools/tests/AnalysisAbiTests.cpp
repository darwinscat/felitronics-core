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
#include <alloc_counter.h>   // installs the allocation counter: EVERY form of `new`, over-aligned included

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <vector>

// The allocation oracle, the repo's usual global-operator-new idiom, counting BYTES as well as calls: the
// P81 queries below claim a number, and the only witness that is not the same arithmetic saying so twice is
// what prepare() actually asks the heap for.
//
// IT COUNTS WHAT THE CONTAINER ASKED FOR, WHICH IS THE QUANTITY A BUDGET STATES — and on MSVC's STL that is
// not what operator new is handed. On x86/x64 it asks for sizeof(void*) + 31 more on a block of 4096 bytes
// or upwards (its hand alignment, <xmemory>), so the correction below takes that back off. It is the same
// correction, for the same reason, as ReferenceTruePeakMeterTests.cpp and LoudnessConformanceTests.cpp:
// spelled here rather than shared because `operator new` has to be defined once per PROGRAM and each of
// these suites is its own. MEASURED, not assumed: this suite was green on macOS and on deb and failed
// eleven checks on `win` — the calibration below caught it first and by itself, reporting a 4096-double
// P52: the allocation counter is `test_support/alloc_counter.h` — ONE definition per executable, and it
// installs EVERY replaceable form, the over-aligned ones included. The private copy that used to live here
// replaced only the two default-aligned forms, so an over-aligned request (core::AlignedVector, and every
// analyzer scratch that rides one) was invisible to the very oracle this file exists to be. The MSVC
// container-padding correction it carried is in the shared header too, as `alloc::containerBytes`.

// The independent side of the P81 oracle: the core's own storageFor(), called here with the arguments the
// shim is supposed to pass. Not an independent implementation of the formula — a second copy of THAT is
// exactly what law 11d exists to prevent — but independent of the SHIM, which is where this task's mutants
// live.
#include <felitronics/core/Config.h>
#include <felitronics/analysis/BandBursts.h>
#include <felitronics/analysis/HumDetector.h>
#include <felitronics/analysis/LowEnd.h>
#include <felitronics/analysis/SourceForensics.h>
#include <felitronics/analysis/ProgrammeReport.h>
#include "fcore_probe.h"

extern "C"
{
    int           fc_probe_report_run      (const float*, std::uint32_t, std::uint32_t, double);
    std::uint32_t fc_probe_report_count_rows (void);
    std::uint32_t fc_probe_report_value_rows (void);
    std::uint32_t fc_probe_report_counts   (double*, std::uint32_t);
    std::uint32_t fc_probe_report_values   (double*, std::uint32_t);
    std::uint32_t fc_probe_report_names    (char*, std::uint32_t);

    int           fc_probe_crest_run       (std::int32_t, const float*, std::uint32_t, std::uint32_t, double);
    int           fc_probe_crest_run_with  (std::int32_t, const float*, std::uint32_t, std::uint32_t, double,
                                            double, double, double, double, std::int32_t, double, double);
    std::uint32_t fc_probe_crest_scalars_len  (void);
    std::uint32_t fc_probe_crest_block_stride (void);
    std::uint32_t fc_probe_crest_loss_len     (void);
    std::uint32_t fc_probe_crest_scalars   (std::int32_t, double*, std::uint32_t);
    std::uint32_t fc_probe_crest_blocks    (std::int32_t, double*, std::uint32_t);
    std::uint32_t fc_probe_crest_loss      (std::int32_t, double*, std::uint32_t);
    double        fc_probe_crest_storage_bytes (std::uint32_t, double, std::uint32_t);

    int           fc_probe_bursts_run      (const float*, std::uint32_t, std::uint32_t, double);
    int           fc_probe_bursts_run_with (const float*, std::uint32_t, std::uint32_t, double,
                                            double, double, double, double, double, double);
    std::uint32_t fc_probe_bursts_scalars_len (void);
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
    std::uint32_t fc_probe_lowend_band_stride (void);

    // P81 — the price of a measurement, asked before it is paid. Geometry only: no audio pointer, because
    // the whole use is to ask BEFORE the input buffer exists.
    double fc_probe_report_storage_bytes    (std::uint32_t, double);
    double fc_probe_bursts_storage_bytes    (std::uint32_t, double);
    double fc_probe_bursts_storage_bytes_with (std::uint32_t, double, double, double, double, double,
                                               double, double);
    double fc_probe_hum_storage_bytes       (std::uint32_t, double);
    double fc_probe_forensics_storage_bytes (std::uint32_t, double);
    double fc_probe_lowend_storage_bytes    (std::uint32_t, double);
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
    // The guard is not cosmetic: stride 0 marks a getter that writes a fixed-length block rather than
    // rows, and dividing by it is UB — clang folded it away and gcc trapped with SIGFPE, which is how
    // this was found, on the second row and not the first.
    if (stride > 1)
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

//==================================================================================================
// P81 — `fc_probe_<mode>_storage_bytes`: the demand of a preparation, published before it is paid.
//
// WHAT IS BEING GUARDED. The module is built -fno-exceptions, where a failed allocation inside prepare()
// does not return and does not throw: it aborts. These five queries let a caller ask the price first, and
// everything below exists so that the price is the RIGHT one and agrees with the run about what it will
// even accept. Two oracles are used and they are different in kind:
//
//   (1) the core's own storageFor(), called from this file with the arguments the shim is supposed to
//       pass. That is NOT an independent formula — one formula is the whole point of law 11d — but it is
//       independent of the SHIM, and the shim is where the mutants of this task live: a swapped argument,
//       the wrong mode's budget, a hard-coded width, bytes divided by 1024, a narrowing through float.
//   (2) what the FIRST run of each mode actually asks the heap for, counted through a global
//       operator new. That one is outside the object entirely, and it is what would catch storageFor()
//       itself understating the allocation — the failure that would make the whole publication a lie.
//==================================================================================================
namespace
{
using Query = double (*) (std::uint32_t, double);
using Run   = int (*) (const float*, std::uint32_t, std::uint32_t, double);

// The core's budget for a geometry, per mode, with this ABI's own arguments: default parameters, and
// Probe::kChunk as the one maxBlock any of the five is sized by.
double oracleReport (std::uint32_t ch, double sr)
{
    const auto st = felitronics::analysis::ProgrammeReport::storageFor (
        sr, (int) fcore::Probe::kChunk, (int) ch, felitronics::analysis::ProgrammeReportParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}
double oracleBursts (std::uint32_t ch, double sr)
{
    const auto st = felitronics::analysis::BandBursts::storageFor (sr, (int) ch, felitronics::analysis::BandBurstsParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}
double oracleHum (std::uint32_t ch, double sr)
{
    const auto st = felitronics::analysis::HumDetector::storageFor (sr, (int) ch, felitronics::analysis::HumDetectorParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}
double oracleForensics (std::uint32_t ch, double sr)
{
    const auto st = felitronics::analysis::SourceForensics::storageFor (sr, (int) ch, felitronics::analysis::SourceForensicsParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}
double oracleLowEnd (std::uint32_t ch, double sr)
{
    const auto st = felitronics::analysis::LowEnd::storageFor (sr, (int) ch, felitronics::analysis::LowEndParams {});
    return st.ok ? (double) st.bytes() : 0.0;
}

// `acceptsEmpty` is the ONE deliberate asymmetry between the two roads: four of the five report on an
// empty programme and lowend refuses one, because `fcore_measure lowend` refuses it too. `silent` is that
// mode's first list getter, read to check that a refused run left nothing behind.
struct Priced { const char* name; Query query; Run run; double (*oracle) (std::uint32_t, double);
                bool acceptsEmpty; Getter silent; };

const Priced priced[] = {
    { "report",    fc_probe_report_storage_bytes,    fc_probe_report_run,    oracleReport,
      true,  fc_probe_report_counts },
    { "bursts",    fc_probe_bursts_storage_bytes,    fc_probe_bursts_run,    oracleBursts,
      true,  fc_probe_bursts_scalars },
    { "hum",       fc_probe_hum_storage_bytes,       fc_probe_hum_run,       oracleHum,
      true,  fc_probe_hum_scalars },
    { "forensics", fc_probe_forensics_storage_bytes, fc_probe_forensics_run, oracleForensics,
      true,  fc_probe_forensics_scalars },
    { "lowend",    fc_probe_lowend_storage_bytes,    fc_probe_lowend_run,    oracleLowEnd,
      false, fc_probe_lowend_scalars },
};

// THE GEOMETRY THE ALLOCATION ORACLE IS MEASURED AT, and the widest the ABI has rather than the ordinary
// one. Each analyzer here is a function-local static, so each mode gets exactly ONE first run per process
// and exactly one chance at an exact equality; spending it on 48 kHz stereo makes the assertion blind to
// everything a run might do CONDITIONALLY — half a megabyte extra above two channels, or at any rate but
// 48 kHz, were two mutants the crew's testing round wrote, and both survived while this was measured at
// (2, 48000). At the widest geometry both conditions are true, and it is also the row the whole task is
// about: hum's 352 688 184 bytes.
constexpr std::uint32_t kFirstCh = 16u;
constexpr double        kFirstSr = 768000.0;

// What `f` asked the heap for, in bytes.
template <typename F>
unsigned long long asked (F&& f)
{
    // A DELTA, not a reset: the shared counter is always armed and other TUs in this binary may be
    // allocating around us, so what this measures is what `f` added — which is what the budget is about.
    const long long before = felitronics::test::alloc::bytes.load();
    f();
    return (unsigned long long) (felitronics::test::alloc::bytes.load() - before);
}

// %.17g and not %.0f. A rate of 1000.5 printed through %.0f reads "1000" — which is exactly how a
// mislabelled number got into this task's own notes before the crew round caught it. A witness that
// renames its own coordinates is not a witness.
std::string num (double v) { char b[64]; std::snprintf (b, sizeof b, "%.17g", v); return b; }
std::string bytesOf (double v) { char b[64]; std::snprintf (b, sizeof b, "%.0f", v); return b; }

// ---- the counter is calibrated before it is believed ----------------------------------------------
// An exact equality against an allocation count is only as good as the count. MSVC's STL hand-aligns any
// block of 4096 bytes or more — a threshold, not a property of the element type; the vector below holds
// ordinary `double` — and an allocator that asked for more than the container wanted would make every
// assertion below fail for a reason that has nothing to do with this ABI. So the
// counter is asked a question whose answer is known first — the same calibration
// ReferenceTruePeakMeterTests makes before its own law-11d equality.
void theByteCounterCountsWhatWasAsked()
{
    felitronics::test::group ("storage_bytes — the allocation oracle is calibrated before it is trusted");
    unsigned long long counted = 0;
    counted = asked ([&]
    {
        std::vector<double> v;
        v.assign (4096, 0.0);                       // 32 768 bytes: a padded block on MSVC's STL
        volatile double* sink = v.data();
        sink[0] = 1.0;
    });
    ok (counted == 4096ull * sizeof (double),
        "the byte counter reports a 4096-double vector as " + std::to_string (counted)
        + " bytes, which is what its container asked for"
        + (felitronics::test::alloc::kStlBigPad != 0
               ? "  [after taking back this STL's " + std::to_string (felitronics::test::alloc::kStlBigPad)
                 + "-byte big-block padding]"
               : ""));
}

// ---- before a single measurement has run ---------------------------------------------------------
void pricesOnAVirginModule()
{
    felitronics::test::group ("storage_bytes — the price is askable before anything has run, and asking costs nothing");
    for (const Priced& m : priced)
    {
        double b = 0.0;
        const unsigned long long cost = asked ([&] { b = m.query (2u, 48000.0); });
        ok (b > 0.0, std::string (m.name) + "_storage_bytes: answers " + bytesOf (b)
                     + " bytes at 48 kHz stereo before any run");
        ok (cost == 0, std::string (m.name) + "_storage_bytes: and asked the heap for " + std::to_string (cost)
                       + " bytes doing it — a query that prepared anything would defeat its own purpose");
    }
}

// ---- the shim against the core's own budget -------------------------------------------------------
void theShimQuotesTheCoreBudget()
{
    felitronics::test::group ("storage_bytes — the quoted price IS the core's published demand");
    const std::uint32_t widths[] = { 1u, 2u, 16u };
    const double rates[] = { 44100.0, 48000.0, 96000.0, 768000.0 };
    for (const Priced& m : priced)
        for (std::uint32_t ch : widths)
            for (double sr : rates)
            {
                const double got = m.query (ch, sr), want = m.oracle (ch, sr);
                // `want > 0` beside the equality, or the row is `0 == 0` the day the core starts refusing
                // a rate it used to accept — green, and comparing nothing. Every one of these sixty rows
                // is a geometry all five modes price today.
                ok (got == want && want > 0.0,
                    std::string (m.name) + "_storage_bytes(" + std::to_string (ch) + ", " + num (sr)
                    + ") = " + bytesOf (got) + ", the core's own budget for those arguments, and positive");
            }
    // TWO MUTANTS THIS CANNOT KILL, measured and named rather than left for the next reader to rediscover:
    //   * masking the low three bits off every price. Every demand in the accepted domain is a multiple
    //     of 8 — 238 992 of 238 992 rows over a random sweep of the domain — so no input distinguishes it.
    //     Masking off four bits IS caught (142 448 of those rows are not multiples of 16), which is how
    //     the boundary was located. Every "exactly" below is therefore exact to the byte and pinned
    //     modulo 8.
    //   * hard-coding the WIDTH in the lowend query. LowEnd's demand does not depend on it: its spectrum
    //     is sized for two fixed mid/side axes, so 7 101 976 bytes at 48 kHz is the answer at every width
    //     from 1 to 16. The same mutant in any of the other four dies on the first row below.
    //
    // The headline number as a LITERAL, so a mutant that carried the byte count through a float — where
    // 352688184 is not representable — fails here even if every relative comparison above still held.
    ok (fc_probe_hum_storage_bytes (16u, 768000.0) == 352688184.0,
        "hum_storage_bytes(16, 768000) is exactly 352688184 bytes");
    // And the modes must not be interchangeable, or a copy-paste between the five would be invisible.
    // ALL TEN PAIRS, not a hand-picked three: "the five are distinct" is a claim about every pair, and three
    // inequalities do not make it. A query wired to the wrong analyzer has to be visible whichever two got
    // swapped.
    {
        int pairs = 0, equal = 0;
        for (std::size_t a = 0; a < 5; ++a)
            for (std::size_t b = a + 1; b < 5; ++b, ++pairs)
                if (priced[a].query (16u, 768000.0) == priced[b].query (16u, 768000.0)) ++equal;
        ok (equal == 0,
            "all " + std::to_string (pairs) + " pairs of the five prices differ at 16 x 768000 — a query "
            "wired to the wrong analyzer is visible whichever two were swapped");
    }

    // THE ARGUMENT ORDER, and exactly the part of it the rows above cannot see. A shim that reads the two
    // the other way round THROUGHOUT — guard included — still refuses everything the positive rows ask
    // about, so those rows catch it; (48000, 2) is what catches the other half, a shim that takes the two
    // along in the wrong order past a guard that was left alone. Measured, both ways: the swap-everything
    // mutant dies here, the swap-only-the-call mutant dies on `answers 0 bytes at 48 kHz stereo`.
    // Neither row checks the declared SIGNATURE: the prototypes in this file are a second copy of the
    // five, and on arm64 a (uint32_t, double) declared as (double, uint32_t) leaves both values in the
    // same register banks and would pass anyway. That is pinned where the call is actually positional —
    // from JavaScript, in tools/wasm/storage-probe.mjs.
    for (const Priced& m : priced)
        ok (m.query (48000u, 2.0) == 0.0,
            std::string (m.name) + "_storage_bytes(48000, 2) is refused — the first argument is the WIDTH");

    // THE SHARED FLOOR (P51): every mode reads the core's 8000 Hz, so these rows sit on it from both sides for each
    // of the five — a shim wired to the wrong analyzer, or an analyzer whose alias drifted, disagrees here first.
    // Bursts is priced at neither: its own band floor lies above the shared one (next row).
    const double floorHz = felitronics::core::kMinSampleRate, under = std::nextafter (floorHz, 0.0);
    ok (floorHz == 8000.0, "the shared floor is 8000 Hz — a literal pin");
    for (const Priced& m : priced)
    {
        const bool band = std::string (m.name) == "bursts";
        ok (m.query (2u, under) == 0.0 && (band || m.query (2u, floorHz) > 0.0) && m.query (2u, 44.1) == 0.0,
            std::string (m.name) + ": refused one ulp under 8000 Hz and at 44.1"
            + (band ? " (and its own band floor is higher still)" : ", priced at 8000"));
    }
    // Each mode may ALSO have an admission floor of its own, above the shared one or below it. Only bursts' is above
    // it at the parameters this ABI measures with; forensics' 2 kHz plateau span used to be the one below, and is
    // now subsumed. A grid that only visits 44.1/48 kHz cannot see either.
    ok (fc_probe_bursts_storage_bytes (2u, 18367.0) == 0.0 && fc_probe_bursts_storage_bytes (2u, 18368.0) > 0.0,
        "bursts: refused at 18367 Hz and priced at 18368 — its own 9 kHz band decides that, not a shared bound");
    ok (fc_probe_forensics_storage_bytes (2u, 1999.999) == 0.0 && fc_probe_forensics_storage_bytes (2u, 2000.0) == 0.0,
        "forensics: its 2 kHz plateau span no longer decides anything — the shared floor refuses both sides of it");
    // A FRACTIONAL RATE, which a shim narrowing the rate to an integer would lose — and chosen so that it loses it
    // whichever way it narrows. With N = 2^15 (the AUTO order at these rates) and bin = fs / N:
    //   bandHi    = ceil ((60.5 * 8 + 3 * bin + 10) / bin)   1995 at 8130, 1994 at 8130.5 and 8131
    //   stretchHi = ceil ((120 + 13 + bin) / bin)             538 at 8130 and 8130.5, 537 at 8131
    // so each half-hertz step drops one double (8 B) from a different row. (The witness was 1000 / 1000.5 Hz
    // before P51, below the floor now.)
    ok (fc_probe_hum_storage_bytes (1u, 8130.0) == 1110784.0,
        "hum(1, 8130) is exactly 1110784 bytes — the rate is read as given");
    ok (fc_probe_hum_storage_bytes (1u, 8130.5) == 1110776.0 && fc_probe_hum_storage_bytes (1u, 8131.0) == 1110768.0,
        "hum(1, 8130.5) is exactly 1110776 — neither 8130's price nor 8131's (1110768): a FRACTIONAL rate survives");
}

// ---- refused geometries ---------------------------------------------------------------------------
void refusedGeometriesQuoteACanonicalZero()
{
    felitronics::test::group ("storage_bytes — a refused geometry quotes +0.0, and 0 cannot be read as free");
    // THE TRAP THIS GATE EXISTS FOR, asserted rather than described: a refused Storage is not an empty
    // Storage. ProgrammeReport's carries DeterministicLoudnessMeter::Storage, whose bytes() has a constant
    // ring term, so a default-constructed one is worth 2400 bytes. A query that forwarded st.bytes()
    // without looking at st.ok would quote that for a measurement that cannot happen — and only for
    // `report`, so probing the other four would show nothing.
    ok (felitronics::analysis::ProgrammeReport::Storage {}.bytes() > 0,
        "a refused ProgrammeReport::Storage is worth "
        + std::to_string ((unsigned long long) felitronics::analysis::ProgrammeReport::Storage {}.bytes())
        + " bytes, not zero — which is what makes the st.ok gate load-bearing");

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const double badRates[] = { 0.0, -0.0, -1.0, 44.1, 999.0, 7999.0, std::nextafter (felitronics::core::kMinSampleRate, 0.0),
                                768000.5, inf, -inf, nan, 1.0e300 };
    const std::uint32_t badWidths[] = { 0u, 17u, 99u, 0x7FFFFFFFu, 0x80000000u, 0xFFFFFFFFu };
    for (const Priced& m : priced)
    {
        for (double sr : badRates)
        {
            const double b = m.query (2u, sr);
            ok (b == 0.0 && ! std::signbit (b),
                std::string (m.name) + "_storage_bytes: a sample rate of " + num (sr) + " quotes +0.0");
        }
        for (std::uint32_t ch : badWidths)
        {
            const double b = m.query (ch, 48000.0);
            ok (b == 0.0 && ! std::signbit (b),
                std::string (m.name) + "_storage_bytes: a width of " + std::to_string (ch) + " quotes +0.0");
        }
    }
}

// ---- what the run actually asks for ----------------------------------------------------------------
// THE FIRST RUN OF EACH MODE IN THIS PROCESS, and it has to be: `std::vector::assign` never gives capacity
// back, so a second run at any geometry asks for less than the demand — or nothing — and the comparison
// would pass on a budget that was far too small.
void theBudgetBoundsWhatTheRunAsks()
{
    felitronics::test::group ("storage_bytes — the published demand bounds what a FIRST run asks the heap for");
    // Its own fixture, at the full width — see kFirstCh / kFirstSr above for why the widest geometry and
    // not the ordinary one.
    const std::uint32_t frames = 2048;
    const std::vector<float> planar = fixture ((int) frames, felitronics::core::kMaxChannels);
    for (const Priced& m : priced)
    {
        const std::uint32_t ch = kFirstCh;
        const double sr = kFirstSr;
        const double budget = m.query (ch, sr);
        int accepted = 0;
        const unsigned long long spent = asked ([&] { accepted = m.run (planar.data(), frames, ch, sr); });
        ok (accepted == 1, std::string (m.name) + "_run: accepted the fixture at " + std::to_string (ch)
                           + " x " + num (sr));
        ok ((double) spent <= budget,
            std::string (m.name) + ": the first run asked for " + std::to_string (spent)
            + " bytes against a published " + bytesOf (budget));
        ok ((double) spent == budget,
            std::string (m.name) + ": and asked for EXACTLY that — the publication is the allocation, not an "
            "upper bound someone guessed. prepare(), process() and finish() together, so a byte asked for "
            "outside the budget in ANY of the three fails here");
    }
    // A SECOND run at the same geometry, so the sentence "a demand, not a footprint" is a measurement in
    // this file rather than a claim in a comment. Four of the five ask for nothing at all the second time
    // (vector::assign keeps its capacity); SourceForensics re-allocates part of its budget on every
    // prepare(). The assertion is the law — never more than was published — and the numbers are in the
    // message, so a change in that asymmetry is READ rather than silently locked in.
    felitronics::test::group ("storage_bytes — a second run at the same geometry stays inside the same demand");
    for (const Priced& m : priced)
    {
        const double budget = m.query (kFirstCh, kFirstSr);
        int accepted = 0;
        const unsigned long long spent = asked ([&] { accepted = m.run (planar.data(), frames, kFirstCh, kFirstSr); });
        // the status is READ, not discarded: a refused second run spends nothing and would satisfy the
        // budget assertion while proving the opposite of what it is here for
        ok (accepted == 1 && (double) spent <= budget,
            std::string (m.name) + ": the second run was accepted and asked for " + std::to_string (spent)
            + " of the same " + bytesOf (budget) + " bytes");
    }
}

// ---- what the query deliberately does NOT mirror -----------------------------------------------------
// Written down as assertions so that a later reader does not "fix" the asymmetry. The query answers about
// a GEOMETRY; `_run` also judges a BUFFER, and the two refusals are different things. A price above zero
// therefore never promises that this particular call will be accepted.
void theQueryPricesGeometryAndNotTheBuffer (const std::vector<float>& planar, std::uint32_t frames)
{
    felitronics::test::group ("storage_bytes — a price is not a promise about the buffer");
    ok (fc_probe_lowend_storage_bytes (2u, 48000.0) > 0.0
        && fc_probe_lowend_run (planar.data(), 0u, 2u, 48000.0) == 0,
        "lowend: the geometry is priced, and an EMPTY programme is still refused — lowend refuses one where "
        "the other four report on it, and that is an input contract, not a geometry");
    ok (fc_probe_hum_storage_bytes (2u, 48000.0) > 0.0
        && fc_probe_hum_run (nullptr, frames, 2u, 48000.0) == 0,
        "hum: the geometry is priced, and a null buffer is still refused");
    ok (fc_probe_report_storage_bytes (2u, 48000.0) > 0.0
        && fc_probe_report_run (planar.data(), 0xFFFFFFFFu, 2u, 48000.0) == 0,
        "report: the geometry is priced, and a frame count whose span passes 4 GiB is still refused — the "
        "query has no way to know how long the programme is, which is the point of it not asking");
}

// ---- the two roads refuse the same geometries -------------------------------------------------------
// Runs LAST: its final rows are refusals, so it leaves every analyzer with no readable result at all.
//
// IT BUILDS ITS OWN FIXTURE, AT THE FULL WIDTH, and that is not tidiness. This ABI is planar: `_run` reads
// plane c at `planar + c * frames`, and the caller's promise is frames * channels floats. Handing a
// two-channel buffer to a sixteen-channel row reads 14 planes past the end of it — 2.7 MB past, with the
// fixture this file already had. macOS let that through and every check stayed green; deb SEGFAULTED on
// it, which is the whole reason a change to this ABI is run on both rows. Short frames, full width: a
// grid of 1170 runs does not need a long programme, it needs an in-bounds one.
void theRefusalSetsAgree()
{
    felitronics::test::group ("storage_bytes — priced above zero exactly where _run is accepted");
    const std::uint32_t frames = 1024;
    const std::vector<float> planar = fixture ((int) frames, felitronics::core::kMaxChannels);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    // The rows that matter are the ones where the five DISAGREE — bursts' own admission floor, above the shared one
    // since P51 (999, 1000 and 2000 now refuse in all five) —
    // and the ones adjacent to a bound, where a `<` written for a `<=` lives. std::nextafter walks to the
    // neighbouring representable double, which is a finer step than any decimal literal can spell.
    // The low bound is the core's (P51); 999, 1000 and 2000 stay as rows — the old floor and forensics' own, both
    // below it now — and 44.1 is the unit error the floor exists for.
    const double lo = felitronics::core::kMinSampleRate, hi = 768000.0;
    const double rates[] = { 0.0, -1.0, -48000.0, 44.1, 999.0, 1000.0, 2000.0,
                             std::nextafter (lo, 0.0), lo, std::nextafter (lo, hi), 8130.5,
                             11025.0, 18367.0, 18368.0, 44100.0, 48000.0,
                             96000.0, 192000.0, 705600.0,
                             std::nextafter (hi, 0.0), hi, std::nextafter (hi, inf), 768000.5,
                             inf, -inf, nan };
    const std::uint32_t widths[] = { 0u, 1u, 2u, 6u, 16u, 17u, 0x7FFFFFFFu, 0x80000000u, 0xFFFFFFFFu };
    // AND EVERY PROGRAMME LENGTH. An EMPTY programme is the road on which the two can actually disagree, and
    // the grid did not walk it: `_run` skips process() when frames == 0, so a run that ignored what
    // prepare() returned still reached finish(), still set its `have*` flag and still answered 1 — at a
    // geometry priced at zero, serving the PREVIOUS programme's numbers out of its getters. Measured:
    // with that one line mutated in forensics, this suite stayed green at 574 checks while
    // `forensics_run(buf, 0, 2, 999.0)` returned 1 and published 29 scalars reporting 48 kHz. Four of
    // the five accept an empty programme; lowend refuses one, and `acceptsEmpty` is how that deliberate
    // asymmetry is expressed here rather than excused.
    // THREE FRAME COUNTS, and the short one is not decoration either: a run that refused programmes below
    // a hundred frames survived a grid of {0, 1024} with nothing to say. Empty, shorter than one analysis
    // window, and longer than one.
    const std::uint32_t counts[] = { 0u, 64u, 1024u };
    int rows = 0, disagreed = 0, spoke = 0;
    std::vector<double> after (64, 0.0);
    for (const Priced& m : priced)
        for (double sr : rates)
            for (std::uint32_t ch : widths)
                for (std::uint32_t n : counts)
                {
                    const bool priceable = m.query (ch, sr) > 0.0;
                    const bool expected  = priceable && (n != 0u || m.acceptsEmpty);
                    const bool runnable  = m.run (planar.data(), n, ch, sr) == 1;
                    ++rows;
                    if (runnable != expected)
                    {
                        ++disagreed;
                        ok (false, std::string (m.name) + ": priced " + (priceable ? "above zero" : "at zero")
                                   + " but _run " + (runnable ? "accepted" : "refused")
                                   + " at " + std::to_string (ch) + " x " + num (sr)
                                   + " with " + std::to_string (n) + " frames");
                    }
                    // AND A REFUSED RUN LEFT NOTHING READABLE. Half of the defect above is not the
                    // verdict, it is what the getters answer afterwards.
                    if (! runnable && m.silent (after.data(), 60u) != 0) ++spoke;
                }
    ok (disagreed == 0, "the price and the run agree about every one of " + std::to_string (rows)
                        + " rows — five modes, twenty-six rates, nine widths and three programme lengths; "
                        "bursts' own admission floor and both neighbours of both shared bounds among them");
    ok (spoke == 0, "and every refused run of those left its getters silent");
}

// ---- a query disturbs nothing ------------------------------------------------------------------------
std::vector<double> snapshot()
{
    const Getter all[] = { fc_probe_report_counts, fc_probe_report_values, fc_probe_bursts_scalars,
                           fc_probe_bursts_chan, fc_probe_bursts_events, fc_probe_bursts_ioi,
                           fc_probe_bursts_lag, fc_probe_hum_scalars, fc_probe_hum_chan, fc_probe_hum_cand,
                           fc_probe_hum_harm, fc_probe_hum_stretch, fc_probe_forensics_scalars,
                           fc_probe_forensics_wall, fc_probe_forensics_grid, fc_probe_forensics_khist,
                           fc_probe_lowend_scalars, fc_probe_lowend_hist, fc_probe_lowend_series,
                           fc_probe_lowend_bands };
    std::vector<double> out;
    for (Getter g : all)
    {
        std::vector<double> buf (8192, 0.0);
        const std::uint32_t n = g (buf.data(), (std::uint32_t) buf.size());
        out.push_back ((double) n);
        out.insert (out.end(), buf.begin(), buf.end());
    }
    return out;
}

// ---- the table, so the two TIERS can be compared and not only two builds of one ----------------------
// `felitronics_analysis_abi_tests --storage-table` prints nothing but this and exits. The wasm twin,
// tools/wasm/storage-probe.mjs, prints the same rows from the built module, and CI compares them with
// `cmp`. That is a CROSS-TIER oracle and not a self-comparison: wasm32 is a 32-bit target with its own
// `std::size_t`, its own struct layouts and its own libm, and a release-vs-debug diff of one tier would
// happily compare a wrong number with itself. Only integer sample rates are listed, so "48000" is
// "48000" on both sides and the comparison is of numbers rather than of two printf dialects.
//==================================================================================================
// K3 — THE BAND IS AN ARGUMENT. `fc_probe_bursts_run` answers one question (sibilance, 5-9 kHz);
// `fc_probe_bursts_run_with` is the SAME detector pointed somewhere else, which is what makes "how dense
// are the transients" answerable without a second analyzer.
//
// THREE THINGS HAVE TO BE TRUE, and only the first is obvious:
//   * the parameterised road at the documented defaults IS the default road — bit for bit, or there are
//     quietly two detectors;
//   * the scalars report the parameters THIS measurement ran at. They used to read a constant, which was
//     true while one road existed and would have become a plausible lie the moment a second one did: the
//     thresholds would have described 6/3 dB while the events came from the caller's numbers;
//   * a refused parameterised run leaves NOTHING half-swapped. The established contract here is that a
//     refusal goes silent, and the state to rule out is new parameters sitting beside old events, which
//     is the one combination a reader cannot detect from the outside.
void k3TheBandIsAnArgument (int ch)
{
    felitronics::test::group ("K3 — bursts with the caller's band: the same detector, pointed elsewhere");
    const int frames = 5 * 48000;
    const std::vector<float> bursty = burstyFixture (frames, ch);
    const auto n = (std::uint32_t) frames, nc = (std::uint32_t) ch;
    const std::uint32_t len = fc_probe_bursts_scalars_len();

    auto scalars = [&] ()
    {
        std::vector<double> v (len + 8u, kCanary);
        const std::uint32_t got = fc_probe_bursts_scalars (v.data(), len);
        v.resize (got);
        return v;
    };
    // The scalar block's own order, as fc_probe.cpp writes it and bursts-format.mjs reads it.
    enum { kSr = 0, kCh = 1, kHop = 2, kBase = 3, kLo = 4, kHi = 5, kEnter = 6, kExit = 7,
           kSamples = 9, kHops = 10, kEligible = 11, kOnsets = 26, kOnsetsPerSec = 31 };

    ok (len == 32u, "the scalar block is 32 doubles — 31 of them since P81, plus K3's onsets per second");

    // ---- the parameterised road at the documented defaults IS the default road ----
    const bool plain = fc_probe_bursts_run (bursty.data(), n, nc, 48000.0) == 1;
    const std::vector<double> a = scalars();
    const bool withDefaults = fc_probe_bursts_run_with (bursty.data(), n, nc, 48000.0,
                                                        5000.0, 9000.0, 10.0, 2000.0, 6.0, 3.0) == 1;
    const std::vector<double> b = scalars();
    std::size_t differing = 0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) if (! (a[i] == b[i])) ++differing;
    ok (plain && withDefaults && a.size() == b.size() && differing == 0,
        "run_with at the documented defaults answers the default road's scalars, bit for bit ("
        + std::to_string (differing) + " differ)");

    // ---- the thresholds in the scalars are the ones THIS run installed, not a constant ----
    const bool moved = fc_probe_bursts_run_with (bursty.data(), n, nc, 48000.0,
                                                 80.0, 8000.0, 20.0, 500.0, 9.0, 4.5) == 1;
    const std::vector<double> c = scalars();
    ok (moved && c[kEnter] == 9.0 && c[kExit] == 4.5,
        "the scalars report the thresholds this measurement ran at (" + std::to_string (c[kEnter]) + " / "
        + std::to_string (c[kExit]) + " dB), not the documented 6 / 3");
    ok (moved && c[kLo] == 80.0 && c[kHi] == 8000.0,
        "... and the band it ran at, read back out of the detector");
    // The HOP is the detector's own answer in SAMPLES, so it is evidence that the millisecond argument
    // reached the core rather than being stored beside it: 20 ms at 48 kHz is 960 samples, not 480.
    ok (moved && c[kHop] == 960.0 && c[kBase] == 25.0,
        "... and 20 ms of hop is 960 samples with a 25-hop baseline — the arguments reached the detector ("
        + std::to_string (c[kHop]) + " / " + std::to_string (c[kBase]) + ")");

    // ---- ONSETS PER SECOND, against the other published scalars. The analyzer owns the denominator; this
    // recomputes it from numbers the same block publishes, so a changed denominator shows up here.
    {
        const double sec = c[kEligible] * c[kHop] / c[kSr];
        const double want = sec > 0.0 ? c[kOnsets] / sec : 0.0;
        ok (std::fabs (c[kOnsetsPerSec] - want) <= 1.0e-12 * std::max (1.0, std::fabs (want)),
            "onsets per second is the onset count over the JUDGED programme, not the file ("
            + std::to_string (c[kOnsetsPerSec]) + "/s over " + std::to_string (sec) + " s)");
        // And the denominator is not the file: the first `baselineHops` are not judged, so the two differ.
        const double fileSec = c[kSamples] / c[kSr];
        ok (fileSec > sec, "... and the judged stretch is shorter than the file (" + std::to_string (sec)
                           + " s of " + std::to_string (fileSec) + " s), so the two denominators differ");
    }

    // ---- A REFUSED PARAMETERISED RUN LEAVES THE PREVIOUS MEASUREMENT WHOLE, parameters included ----
    {
        // 30 kHz cannot be the top corner at 48 kHz: 0.49*fs is 23520.
        const bool refused = fc_probe_bursts_run_with (bursty.data(), n, nc, 48000.0,
                                                       20000.0, 30000.0, 10.0, 2000.0, 6.0, 3.0) == 0;
        const std::vector<double> after = scalars();
        ok (refused && after.empty(), "a band the core cannot build is refused, and the getters go silent");
        // AND THE NEXT SUCCESSFUL RUN IS NOT CONTAMINATED BY IT. Silence alone does not rule out the half
        // swap: parameters recorded by the refused call would sit in the module and be read beside the NEXT
        // measurement's events. So run the defaults again and check the thresholds came back to 6 / 3.
        const bool again = fc_probe_bursts_run (bursty.data(), n, nc, 48000.0) == 1;
        const std::vector<double> back = scalars();
        ok (again && back.size() == a.size() && back[kEnter] == 6.0 && back[kExit] == 3.0
            && back[kLo] == 5000.0 && back[kHi] == 9000.0,
            "... and the run after it reports ITS OWN parameters, not the refused call's");
    }

    // ---- THE PRICE MOVES WITH THE PARAMETERS, which is why a second one had to exist: the baseline ring
    // is round(baselineMs / hopMs) hops of hopMs each, so a page sizing itself by the default figure would
    // be short exactly where it asked for a longer memory.
    {
        const double dflt = fc_probe_bursts_storage_bytes (nc, 48000.0);
        const double same = fc_probe_bursts_storage_bytes_with (nc, 48000.0, 5000.0, 9000.0, 10.0, 2000.0, 6.0, 3.0);
        ok (dflt > 0.0 && same == dflt, "the parameterised price at the documented defaults is the default price");
        const double longer = fc_probe_bursts_storage_bytes_with (nc, 48000.0, 5000.0, 9000.0, 10.0, 8000.0, 6.0, 3.0);
        ok (longer > dflt, "a four-times-longer baseline costs more (" + std::to_string (longer)
                           + " over " + std::to_string (dflt) + " bytes)");
        ok (fc_probe_bursts_storage_bytes_with (nc, 48000.0, 20000.0, 30000.0, 10.0, 2000.0, 6.0, 3.0) == 0.0,
            "and a band the run refuses is priced at the canonical zero, as the default query does");
    }
}

//==================================================================================================
// K1 — TWO SLOTS AND A COMPARISON. Every other analyzer here holds one result, because every other question is
// about one programme; this one is about the DIFFERENCE between two, so both live in the module and the loss
// is computed here rather than reassembled from printed columns on the other side of the ABI.
void k1TheTwoSlotsAndTheLoss (int ch)
{
    felitronics::test::group ("K1 — crest: two slots, a comparison, and what each refusal does");
    const int frames = 3 * 48000;
    const std::vector<float> src = burstyFixture (frames, ch);
    std::vector<float> quiet = src;
    for (auto& v : quiet) v = std::clamp (v, -0.2f, 0.2f);     // a "master" with its peaks taken off
    const auto n = (std::uint32_t) frames, nc = (std::uint32_t) ch;
    const std::uint32_t slen = fc_probe_crest_scalars_len(), stride = fc_probe_crest_block_stride();
    const std::uint32_t llen = fc_probe_crest_loss_len();
    ok (slen == 22u && stride == 15u && llen == 17u,
        "the shapes are 22 scalars (16, one accepted-block count per band, and the programme level), 15 "
        "doubles a block row and 17 loss fields");

    // A getter is silent until its OWN slot has run — the two are independent, which is the whole point.
    std::vector<double> v (64, kCanary);
    ok (fc_probe_crest_scalars (0, v.data(), slen) == 0 && fc_probe_crest_scalars (1, v.data(), slen) == 0,
        "both slots are silent before either has run");
    ok (fc_probe_crest_loss (0, v.data(), llen) == 0, "and the loss refuses while a slot is missing");

    ok (fc_probe_crest_run (0, src.data(), n, nc, 48000.0) == 1, "slot 0 accepts a real programme");
    ok (fc_probe_crest_loss (0, v.data(), llen) == 0,
        "the loss STILL refuses with only one slot filled — a comparison needs two, and answering from one "
        "would be a number about nothing");
    ok (fc_probe_crest_run (1, quiet.data(), n, nc, 48000.0) == 1, "slot 1 accepts the master");
    ok (fc_probe_crest_loss (0, v.data(), llen) == llen, "and now the loss answers");

    // THE SLOTS DO NOT SHARE STATE. Running slot 1 at different parameters must not move slot 0's numbers.
    std::vector<double> a0 (slen, 0.0), a1 (slen, 0.0);
    ok (fc_probe_crest_scalars (0, a0.data(), slen) == slen, "slot 0 answers its scalars");
    ok (fc_probe_crest_run_with (1, quiet.data(), n, nc, 48000.0, 200.0, 1500.0, 7000.0, 50.0, 8, -65.0, -35.0) == 1,
        "slot 1 runs again at other parameters");
    ok (fc_probe_crest_scalars (0, a1.data(), slen) == slen, "slot 0 still answers");
    std::size_t moved = 0;
    for (std::uint32_t i = 0; i < slen; ++i) if (! (a0[i] == a1[i])) ++moved;
    ok (moved == 0, "and slot 0's scalars did not move (" + std::to_string (moved) + " differ)");
    // ... and slot 1's DID, including the parameters it echoes back.
    std::vector<double> b1 (slen, 0.0);
    ok (fc_probe_crest_scalars (1, b1.data(), slen) == slen
        && b1[4] == 200.0 && b1[5] == 1500.0 && b1[6] == 7000.0 && b1[3] == 8.0,
        "while slot 1 reports the parameters THIS run installed, not the defaults");

    // THE BLOCK TABLE IS ALL OF IT OR NONE OF IT — a short capacity writes nothing rather than a prefix that
    // looks like a measurement of a shorter programme.
    {
        const auto rows = (std::uint32_t) a0[12];                 // blockCount
        std::vector<double> big ((std::size_t) rows * stride + 8u, kCanary);
        ok (rows > 0 && fc_probe_crest_blocks (0, big.data(), rows * stride) == rows,
            "the block table fills an exact capacity (" + std::to_string (rows) + " rows)");
        std::vector<double> small ((std::size_t) rows * stride, kCanary);
        ok (fc_probe_crest_blocks (0, small.data(), rows * stride - 1u) == 0,
            "and one double short writes NOTHING");
        bool clean = true;
        for (double d : small) clean = clean && (d == kCanary);
        ok (clean, "... leaving the buffer as it found it");
    }

    // A SLOT THAT DOES NOT EXIST, and a band that does not.
    ok (fc_probe_crest_run (2, src.data(), n, nc, 48000.0) == 0
        && fc_probe_crest_run (-1, src.data(), n, nc, 48000.0) == 0, "a slot outside 0..1 is refused");
    ok (fc_probe_crest_scalars (2, v.data(), slen) == 0 && fc_probe_crest_loss (5, v.data(), llen) == 0,
        "and so are a slot and a band outside their ranges on the getters");

    // A REFUSED RUN GOES SILENT, and does not leave the previous measurement readable under new parameters.
    ok (fc_probe_crest_run_with (0, src.data(), n, nc, 48000.0, 3000.0, 2000.0, 6000.0, 100.0, 4, -70.0, -40.0) == 0,
        "edges that do not rise are refused");
    ok (fc_probe_crest_scalars (0, v.data(), slen) == 0, "and that slot is silent afterwards");

    // THE PRICE MOVES WITH THE PROGRAMME, because the cell store is per hop — which is why this query takes a
    // length where its neighbours take only a geometry.
    const double oneSec = fc_probe_crest_storage_bytes (nc, 48000.0, 48000u);
    ok (oneSec > 0.0 && fc_probe_crest_storage_bytes (nc, 48000.0, 480000u) > oneSec,
        "ten times the programme costs more than one second (" + std::to_string (oneSec) + " B)");
}

void printStorageTable()
{
    const std::uint32_t widths[] = { 1u, 2u, 6u, 16u, 17u };
    // 999 and 1000 straddle the probe's pre-P51 floor and are both refused now; 7999 / 8000 straddle the core's.
    // 12000 and 16000 are priced rates below 22050 — where the analyzers' own geometry changes most — and they
    // replace what 1000 and 2000 used to price. storage-probe.mjs's TABLE_RATES is this list.
    const double rates[] = { 999.0, 1000.0, 7999.0, 8000.0, 11025.0, 12000.0, 16000.0, 22050.0, 44100.0, 48000.0,
                             88200.0, 96000.0, 176400.0, 192000.0, 352800.0, 384000.0, 705600.0, 768000.0, 768001.0 };
    for (const Priced& m : priced)
        for (std::uint32_t ch : widths)
            for (double sr : rates)
                std::printf ("%s %u %.0f %.0f\n", m.name, ch, sr, m.query (ch, sr));
}

void queriesAreStateless()
{
    felitronics::test::group ("storage_bytes — the query is a function of its arguments and of nothing else");
    const std::vector<double> before = snapshot();
    // NOT VACUOUS, AND PER MODE. If a result had already been cleared, both snapshots would agree about
    // its zeros and a query that wrongly cleared it would pass. A count over the WHOLE snapshot is not
    // enough either: four analyzers can carry it past any threshold while the fifth is already silent, so
    // each of the five is asked for itself.
    {
        int live = 0;
        // THE BUFFER IS SIZED FROM THE WIDEST MODE, NOT FROM A NUMBER. It held 64 doubles and asked for
        // 60, which was comfortable until K9 took lowend's scalars to 67 — and then a getter that refuses
        // a short capacity, correctly, read here as "this mode has no result", which is a different
        // statement entirely. A capacity this check is not about must never be the thing it measures.
        const std::uint32_t widest = 4096u;
        std::vector<double> probe ((std::size_t) widest + 8, 0.0);
        for (const Priced& m : priced) if (m.silent (probe.data(), widest) > 0) ++live;
        ok (live == 5, "all five modes have a readable result before the queries — "
                       + std::to_string (live) + " of five — so there is something for a query to disturb");
    }
    double first[5] {};
    for (int i = 0; i < 5; ++i) first[i] = priced[(std::size_t) i].query (2u, 48000.0);
    // ask about other geometries in between, including refused ones
    for (const Priced& m : priced)
    {
        (void) m.query (16u, 768000.0); (void) m.query (0u, 48000.0); (void) m.query (2u, -1.0);
    }
    const std::vector<double> after = snapshot();
    ok (before == after, "twenty getters answer identically across a round of queries — asking the price "
                         "neither clears a result nor reconfigures an analyzer");
    bool same = true;
    for (int i = 0; i < 5; ++i) if (priced[(std::size_t) i].query (2u, 48000.0) != first[i]) same = false;
    ok (same, "and each query repeats its own answer afterwards");
}
} // namespace

int main (int argc, char** argv)
{
    // The cross-tier table, and nothing else, when asked for it: CI compares this against the same rows
    // read out of the wasm module. Before any measurement runs, because these queries are pure.
    if (argc > 1 && std::strcmp (argv[1], "--storage-table") == 0) { printStorageTable(); return 0; }

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

    // ---------- P81: the price, on a module where nothing has run yet ----------
    theByteCounterCountsWhatWasAsked();
    pricesOnAVirginModule();
    theShimQuotesTheCoreBudget();
    refusedGeometriesQuoteACanonicalZero();
    // ...and THEN the first runs of this process, which is the only moment the allocation oracle works:
    // retained vector capacity makes every later run ask for less than the demand.
    theBudgetBoundsWhatTheRunAsks();

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
    // THE STRIDE IS SPELLED, AND CHECKED AGAINST THE MODULE. Spelling it makes this an oracle — a test
    // that read the subject's own stride would move with a stride bug instead of catching it. Checking it
    // makes the spelling survive: it was 11 until K9 widened the band row to 13, and the only symptom was
    // this truncation check quietly measuring a 46-row buffer at the wrong width and expecting 38.
    ok (fc_probe_lowend_band_stride() == 13, "the band row is 13 doubles, and the module agrees ("
        + std::to_string (fc_probe_lowend_band_stride()) + ")");
    hammer (fc_probe_lowend_bands,   "lowend_bands",   13);

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

    k3TheBandIsAnArgument (ch);
    k1TheTwoSlotsAndTheLoss (ch);

    // ---------- P81, last: asking the price disturbs nothing, and the two roads refuse the same set ----------
    // queriesAreStateless() first, while the analyzers still hold the fixtures the block above measured;
    // theRefusalSetsAgree() runs every mode at twenty-six rates and ends on refusals, so every result is
    // cleared by the time it returns and nothing may follow it.
    queriesAreStateless();
    theQueryPricesGeometryAndNotTheBuffer (planar, (std::uint32_t) frames);
    theRefusalSetsAgree();

    return felitronics::test::report();
}
