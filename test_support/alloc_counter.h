// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

// THE ALLOCATION COUNTER — one object, so that "process() does not allocate" means the same sentence in
// every suite that says it.
//
// WHAT IT IS FOR. A suite proves an RT claim by replacing the global allocation functions and reading a
// DELTA across the call under test. That is only as good as the set of functions it replaced, and the set
// is not two. C++17 split `new` by ALIGNMENT: an object whose `alignof` exceeds
// __STDCPP_DEFAULT_NEW_ALIGNMENT__ — 16 on every desktop row, and 8 on wasm32, both measured — is
// allocated through
// `operator new(std::size_t, std::align_val_t)`, a DIFFERENT function. A counter that replaces only
// `operator new(std::size_t)` and `new[]` never sees it.
//
// That was not hypothetical here. On `a9816e2`, 61 translation units replaced the global `operator new` to
// count; 11 installed the over-aligned form and 50 did not, while both kinds said "no heap allocation in
// the audio path" in the same words. The carriers are real and load-bearing: `core::AlignedVector` /
// `SeamAllocator<64>` (core/Fft.h), and through it every convolver and everything downstream of one, plus
// `eq::EqEngine`, whose `alignof` is 64 and whose own header says a two-form counter would not see its
// 331 KiB.
//
// AND THE HOLE WAS LATENT, which is worth writing down rather than leaving as an implication. Fixing the
// instrument and re-measuring found nothing that had fallen through it: 1741 over-aligned call stacks
// across all 128 test binaries, every one of them inside `prepare()`, `setIr()` or `Bank::build()`, not one
// naming `process`, `analyse` or `applyGain`. `modules/eq/tests/` is the sharpest case — the suite the
// blindness would have hurt most makes ZERO over-aligned allocations of its own, because it puts the
// engine on the STACK. A blind instrument is a defect whether or not something walked past it; it is just
// not the same defect as one that did.
//
// SO THIS HEADER REPLACES EVERY REPLACEABLE FORM — all eight `new`s (plain · array · over-aligned ·
// nothrow · and the combinations) and all twelve matching `delete`s. No route a `new` can take is left out,
// so the sentence "the counter saw none" has one meaning. (What a `new` is NOT the only route to is listed
// at the foot of this comment, because a gate's boundary is part of what it proves.)
//
// HOW TO USE IT. Include this header in EXACTLY ONE translation unit of an executable. Including it IS
// installing the counter — there is no macro to remember, deliberately:
//
//   * A macro that must be defined can be forgotten, and a forgotten one is SILENT — the counter reads 0
//     for ever and every "no allocation" assertion built on it passes vacuously. That is the very defect
//     this header exists to end, re-introduced in a new shape.
//   * Including it twice in one executable is not silent. The standard makes a program with two
//     replacements ill-formed, and every linker in this tree's matrix says so out loud — measured, not
//     assumed: `duplicate symbol 'operator new(unsigned long)' … ld: 2 duplicate symbols`.
//
//   #include <alloc_counter.h>     // installs; read the numbers through `alloc::` below
//
// WHAT THE NUMBERS MEAN.
//   alloc::count      every call to any replaced form. Read a DELTA; the self-check below spends eight of
//                     these before main(), so an ABSOLUTE count means nothing and never did.
//   alloc::bytes      requested bytes, corrected for what MSVC's STL adds on top of a CONTAINER's own
//                     request (see kStlBigPad). This is the quantity a `*Bytes()` budget states.
//   alloc::rawBytes   requested bytes exactly as they reached the allocator, no correction.
//   alloc::alignedCount   the over-aligned calls alone — the subset a two-form counter cannot see.
//
// Both totals are kept because the tree's suites mean both: a law-11d budget check (`prepareBytes() ==
// what prepare() asked for`) is a CONTAINER question, while a cap expressed in raw allocator traffic is
// not. Reading the wrong one changes an answer only on MSVC x86/x64, which is exactly where it is
// hardest to notice from here.
//
// A COUNT IS OF REQUESTS, not of successes: an allocation that is then made to fail (see `failNext`) has
// already been counted, and a zero-byte request counts one allocation and zero bytes though a byte is
// handed out. Both were true of the sixty-one counters this replaces, and the budgets calibrated against
// them are the reason to keep it so rather than an oversight.
//
// WHAT THIS INSTRUMENT STILL CANNOT SEE:
//   * `malloc` — it replaces the C++ allocation functions, not the C ones. `modules/fftpffft/pffft/pffft.c`
//     really does call `malloc` directly, so a counter reading zero says nothing about that route. Measured
//     rather than assumed: with `--wrap=malloc,calloc,realloc` on the gcc row, 200 `process()` calls through
//     `MatrixConvolver<PffftRealFft>` make ZERO C allocations (4 in the whole run, all inside `prepare`).
//   * an allocator a dependency reaches by its own syscall or pool — Eigen, for one, which a suite
//     covers with Eigen's own malloc gate instead.
//   * a path no suite drives. This counts what is executed, not what exists.

#include <felitronics_test.h>

#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>
#if defined(_MSC_VER)
 #include <malloc.h>   // _aligned_malloc / _aligned_free
#endif

namespace felitronics::test::alloc
{
    //==============================================================================
    // MSVC's STL asks operator new for more than the container did. On x86/x64 it hand-aligns any block of
    // 4096 bytes or more to 32 and asks for `_Non_user_size` extra — sizeof(void*) + 31, one word more
    // under _DEBUG. A budget states what the CONTAINER wanted and leaves that to the caller, so `bytes`
    // takes it back off. NOT MODELLED: iterator debugging (_ITERATOR_DEBUG_LEVEL), whose container proxy is
    // a separate ALLOCATION rather than padding and cannot be subtracted off a block — under it a
    // calibrated suite fails by name, which is what a calibration is for.
   #if defined(_MSVC_STL_VERSION) && (defined(_M_IX86) || defined(_M_X64))
    #if defined(_DEBUG)
    inline constexpr std::size_t kStlBigPad = 2 * sizeof (void*) + 31;
    #else
    inline constexpr std::size_t kStlBigPad = sizeof (void*) + 31;
    #endif
   #else
    inline constexpr std::size_t kStlBigPad = 0;
   #endif
    inline constexpr std::size_t kStlBigBlock = 4096;

    // ONE EXEMPTION from that correction, and it needs stating because it cannot be inferred from the call.
    // A plain `new T` IS NOT A CONTAINER: it asks for exactly `sizeof(T)`, so subtracting the STL's
    // container padding off it takes away bytes nobody added — 39 of them on the `win` row, which is where
    // this was found. An over-aligned `new T` is already safe (its bytes are counted raw), but an ordinary
    // one over the big-block threshold is indistinguishable from a vector BY SIZE ALONE. So a suite that
    // knows the size of such an object tells the counter, and requests of exactly it are left alone. Zero —
    // the default — exempts nothing.
    inline std::atomic<std::size_t> plainObjectSize { 0 };

    inline long long containerBytes (std::size_t s) noexcept
    {
        // `if constexpr`, not `if`: a runtime branch on a constant is C4127 under MSVC /W4, and this header
        // is read by a downstream that may well compile at /W4 even though this tree is /W3.
        if constexpr (kStlBigPad == 0)
        {
            return (long long) s;
        }
        else
        {
            if (s == plainObjectSize.load (std::memory_order_relaxed)) return (long long) s;
            return (long long) (s >= kStlBigBlock + kStlBigPad ? s - kStlBigPad : s);
        }
    }

    //==============================================================================
    // The totals. `inline` variables, so a multi-TU executable reads the same numbers from any of its TUs
    // while only one of them installs.
    inline std::atomic<long long> count    { 0 };
    inline std::atomic<long long> bytes    { 0 };
    inline std::atomic<long long> rawBytes { 0 };

    // The OVER-ALIGNED allocations alone — the subset fifty suites could not see. Kept as a number of its
    // own so a claim about them can be made directly, and so the question "does this path touch the
    // over-aligned allocator at all?" has an answer that is measured rather than reasoned about.
    inline std::atomic<long long> alignedCount { 0 };

    // FAULT INJECTION, for the suites that ask what an entry point does when the heap says no. Set it and
    // the NEXT allocation fails: an exception from the throwing forms where exceptions exist, a null
    // pointer from the nothrow ones — each the way its own contract spells failure. The flag clears itself,
    // so it arms exactly one allocation. Inert unless set.
    inline std::atomic<bool> failNext { false };

    // RE-ENTRANCY, for the suites that call back INTO the subject from inside one of its allocations. The
    // hook runs once, before the memory is handed back, then disarms itself. Inert unless set.
    // `noexcept` IS PART OF THE TYPE, and it is here on purpose: four of the eight replacements below are
    // themselves `noexcept` (the nothrow forms), so a hook that threw would not report a failed allocation —
    // it would call std::terminate. Making the pointer type noexcept moves that from a runtime abort in a
    // rare fixture to a compile error at the assignment.
    inline std::atomic<void (*)() noexcept> onNext { nullptr };

    //==============================================================================
    // The alignment the probe below uses. It is the tree's own — `core::kSeamAlignment` and
    // `alignof(eq::EqEngine)` are both 64 — and it must EXCEED the default or there is no second function
    // to reach and nothing to probe. A platform where the default were 64 would not merely make this
    // assert; it would make every over-aligned claim in the tree mean something else, so failing to
    // compile is the right answer rather than a quietly weaker check.
    inline constexpr std::size_t kProbeAlign = 64;
    static_assert (kProbeAlign > __STDCPP_DEFAULT_NEW_ALIGNMENT__,
                   "the probe must be over-aligned on this platform, or it proves nothing");

    // EIGHT probes, one per form. The first draft proved the two throwing scalar forms and inferred the rest,
    // and a crew seat showed what that costs: replace the aligned ARRAY form's body with one that allocates
    // WITHOUT counting, and the two-probe gate stays green while every allocation through it escapes
    // accounting. That is the "one anchor, one blind spot" mistake in miniature.
    //
    // THE DISTINCTION THAT MATTERS, because it is not the obvious one and a reader will otherwise "simplify"
    // this back to two probes. DELETING a form is harmless: the standard's default implementations forward —
    // `new[]` to `new`, nothrow to throwing, array-aligned to scalar-aligned — so an unreplaced form lands in
    // a replaced one and the allocation is still counted. Measured on all three libraries (libc++, libstdc++,
    // UCRT): deleting `operator new[](size_t, align_val_t)` still moves `alignedCount`. What is NOT harmless
    // is a form that is present and stops counting, and that is what these eight catch — one mutation per
    // form, eight of eight red, each naming itself.
    //
    // WHAT EACH PROBE WATCHES, and why the two answers differ. The nothrow and array forms are checked
    // against `count`, since their fallback route is counted too. The OVER-ALIGNED forms are checked against
    // `alignedCount`, which only this file's aligned path increments — so a default-aligned allocation
    // happening to land inside the window cannot stand in for the form under test, which is precisely how an
    // over-aligned blind spot hid for as long as it did.
    inline constexpr int kForms = 8;
    inline const char* const kFormNames[kForms] = {
        "new(size_t)", "new[](size_t)",
        "new(size_t, align_val_t)", "new[](size_t, align_val_t)",
        "new(size_t, nothrow_t)", "new[](size_t, nothrow_t)",
        "new(size_t, align_val_t, nothrow_t)", "new[](size_t, align_val_t, nothrow_t)"
    };
    inline long long probeSeen[kForms] = { -1, -1, -1, -1, -1, -1, -1, -1 };

    // PROVE THE INSTRUMENT. An uninstalled or half-installed counter reads 0 for ever, and a suite built on
    // it goes green while measuring nothing — the failure this whole header is an answer to. So before
    // main() runs, ask the counter for one allocation through each form and record whether the number moved,
    // then refuse to let the run continue if any shape went unseen.
    //
    // Direct `::operator new` CALLS, not `new` expressions: the standard lets an implementation elide the
    // allocation of a new-expression whose result it can account for, and this tree has already been bitten
    // by an optimizer removing an allocation a counter was watching. A call to the function itself is not
    // elidable, and the `volatile` sink keeps the memory observably used. (Checked at -O3 and with LTO: no
    // elision.)
    // CODEQL, and where the silence actually lives. `cpp/new-free-mismatch` traces the REPLACED
    // `::operator new` down to the `std::malloc` inside it and then calls the matching `::operator delete`
    // a malloc/delete mismatch. It is not one: these eight lines are the only way to exercise a global
    // replacement at all, because after the replacement `::operator new` IS the counter, and the pairs are
    // exactly what [new.delete] requires — sized with sized, aligned with aligned, nothrow with its plain
    // partner. The rule accepts the two sized plain forms and rejects the six aligned and nothrow ones,
    // which is the shape of a rule limitation and not of a defect.
    // ⚠ THE SUPPRESSION IS NOT HERE, and an inline `// codeql[...]` marker was REMOVED from these lines
    // rather than left in place: this repository's code-scanning setup does not honour them, so the marker
    // read as a working guard while doing nothing — the class of claim this file exists to stop. The six
    // alerts are dismissed as false positives on GitHub, with this reason, where the dismissal is audited
    // and a NEW alert on this construct would still be raised.
    // ⚠ THE PROBE DEALLOCATES THROUGH THE UNSIZED FORMS ON PURPOSE, and it is a portability rule, not a
    // preference. `probe()` sits ABOVE the replacements at the bottom of this header, so the only
    // declarations in scope here are the ones `<new>` provides — and the SIZED deallocation functions are
    // declared there only when sized deallocation is enabled. Clang does not enable it by default on every
    // row: this built on Apple clang and failed on ubuntu clang with "no matching function for call to
    // 'operator delete'" on exactly the four sized calls. What the probe measures is the `new` side, so the
    // deallocation spelling is free to be the portable one — aligned allocations still go back through the
    // ALIGNED delete, which is what routes them to alignedFree rather than std::free.
    inline void probe() noexcept
    {
        constexpr std::size_t n = 128;
        const std::align_val_t a { kProbeAlign };
        auto touch = [] (void* p) noexcept { if (p != nullptr) static_cast<volatile char*> (p)[0] = 1; };

        { const long long c = count.load();        void* p = ::operator new      (n);                 touch (p); ::operator delete   (p);       probeSeen[0] = count.load()        - c; }
        { const long long c = count.load();        void* p = ::operator new[]    (n);                 touch (p); ::operator delete[] (p);       probeSeen[1] = count.load()        - c; }
        { const long long c = alignedCount.load(); void* p = ::operator new      (n, a);              touch (p); ::operator delete   (p, a);    probeSeen[2] = alignedCount.load() - c; }
        { const long long c = alignedCount.load(); void* p = ::operator new[]    (n, a);              touch (p); ::operator delete[] (p, a);    probeSeen[3] = alignedCount.load() - c; }
        { const long long c = count.load();        void* p = ::operator new      (n, std::nothrow);   touch (p); ::operator delete   (p);       probeSeen[4] = count.load()        - c; }
        { const long long c = count.load();        void* p = ::operator new[]    (n, std::nothrow);   touch (p); ::operator delete[] (p);       probeSeen[5] = count.load()        - c; }
        { const long long c = alignedCount.load(); void* p = ::operator new      (n, a, std::nothrow); touch (p); ::operator delete  (p, a);    probeSeen[6] = alignedCount.load() - c; }
        { const long long c = alignedCount.load(); void* p = ::operator new[]    (n, a, std::nothrow); touch (p); ::operator delete[](p, a);    probeSeen[7] = alignedCount.load() - c; }
    }

    // WHERE THE GATE GOES RED. Not in a check the suite could forget to reach, and not in `report()` —
    // which is not universal: the teq-derived suites total in their own harness and never call it. An
    // instrument that cannot be believed makes every conclusion drawn from it void, so the run stops here,
    // before main(), naming the forms that went uncounted. There is no non-fatal version of this failure:
    // the counter is either installed for a form or it is not.
    inline void verify()
    {
        int missing = 0;
        for (int i = 0; i < kForms; ++i) if (probeSeen[i] < 1) ++missing;
        if (missing == 0) return;

        std::fprintf (stderr, "FATAL: the allocation counter is not installed for every form of `new`.\n");
        for (int i = 0; i < kForms; ++i)
            if (probeSeen[i] < 1)
                std::fprintf (stderr, "       UNCOUNTED: operator %s\n", kFormNames[i]);
        std::fprintf (stderr, "       every \"no heap allocation\" claim in this suite would be measuring nothing — P52.\n");
        std::fflush (stderr);
        std::abort();
    }

    // The standing evidence, printed by every counting suite before its own banner: the instrument was
    // proven in THIS run, not in a comment. It is emitted here rather than from report() because report()
    // is not universal — two suites in this tree total in a harness of their own and never call it, and a
    // line that is missing in exactly the places hardest to notice is not evidence.
    //
    // ⚠ STDERR, AND THAT IS LOAD-BEARING. A test binary's stdout is a DATA channel for whoever runs it:
    // `felitronics_analysis_abi_tests --storage-table` writes 425 rows of published demand there, and CI
    // compares them BYTE FOR BYTE against the same table off the wasm module (P81's cross-tier gate).
    // Announcing on stdout put this line at byte 1 of that file and the gate failed on the first byte —
    // the instrument's own proof of health corrupting the answer it was measuring. Diagnostics go to
    // stderr; the run's result goes to stdout. The verify() failure above follows the same rule.
    inline void announce()
    {
        std::fprintf (stderr, "  - allocation counter: all %d forms of `new` counted, over-aligned included "
                              "(probed one allocation through each)\n", kForms);
    }

    //==============================================================================
    // The one place any of this is decided.
    inline void* release (std::size_t s, std::size_t align, bool nothrow) noexcept (false)
    {
        count.fetch_add (1, std::memory_order_relaxed);
        if (align != 0) alignedCount.fetch_add (1, std::memory_order_relaxed);
        // RAW for an over-aligned request, corrected otherwise. kStlBigPad undoes what MSVC's STL adds ON
        // TOP of a container's request; an over-aligned `new` in this tree is an OBJECT whose size is
        // exactly `sizeof`, and taking the correction off it subtracts bytes nobody ever added. The `win`
        // row found that, on the mastering chain: 339 072 read as 339 033, and 88 budget rows with it.
        bytes.fetch_add (align != 0 ? (long long) s : containerBytes (s), std::memory_order_relaxed);
        rawBytes.fetch_add ((long long) s, std::memory_order_relaxed);

        // Claim the pending failure BEFORE the re-entry hook runs. The other order loses it: a hook that
        // allocates would consume the switch, and the allocation that was promised the failure would quietly
        // succeed instead. No fixture arms both today; the order is what keeps that from being load-bearing.
        const bool fail = failNext.exchange (false, std::memory_order_relaxed);

        if (auto hook = onNext.exchange (nullptr, std::memory_order_relaxed)) hook();

        if (fail)
        {
            if (nothrow) return nullptr;
           #if defined(__cpp_exceptions) || defined(_CPPUNWIND)
            throw std::bad_alloc();
           #else
            // The `wasm-audio` tier builds without exceptions, where the analogue of a throwing bad_alloc is
            // the abort itself and no suite outlives it. Returning null from a THROWING form would be a
            // contract violation, so refuse to pretend.
            std::abort();
           #endif
        }

        void* p = nullptr;
        if (align == 0) p = std::malloc (s ? s : 1);
       #if defined(_MSC_VER)
        else p = _aligned_malloc (s ? s : 1, align);
       #else
        else
        {
            // posix_memalign, not std::aligned_alloc: the latter requires the size to be a multiple of the
            // alignment (a C11 rule Apple's libc enforces), which a `sizeof` never has to be. Its alignment
            // must be a multiple of sizeof(void*), hence the floor.
            const std::size_t al = align < sizeof (void*) ? sizeof (void*) : align;
            if (posix_memalign (&p, al, s ? s : 1) != 0) p = nullptr;
        }
       #endif
        // A THROWING allocation function may not return null — it reports failure by throwing, and a caller
        // is entitled to dereference what it got back. Returning null there turns an out-of-memory into a
        // null dereference somewhere else entirely; SourceForensicsTests.cpp's own counter says so and this
        // keeps it. Where there are no exceptions (the `wasm-audio` tier) the honest answer is the abort.
        if (p == nullptr && ! nothrow)
        {
           #if defined(__cpp_exceptions) || defined(_CPPUNWIND)
            throw std::bad_alloc();
           #else
            std::abort();
           #endif
        }
        return p;
    }

    inline void alignedFree (void* p) noexcept
    {
       #if defined(_MSC_VER)
        _aligned_free (p);   // memory from _aligned_malloc; std::free on it corrupts the heap
       #else
        std::free (p);
       #endif
    }

    // Arm the probe before main(). A namespace-scope initializer, so a suite gets it by including the header
    // and cannot forget it.
    //
    // NOTHING ODR-USES `installed`, AND THAT IS NOT A HOLE — the question came up in review and the answer is
    // in [basic.start.dynamic]. A non-block inline variable's dynamic initialization may indeed be deferred
    // past the first statement of main; when it is, it strongly happens before any non-initialization odr-use
    // of a NON-INLINE function or variable defined in the same translation unit. The twenty replacements
    // below are exactly that — non-inline functions defined in this TU — so the first allocation anywhere in
    // the program has already forced this initializer to run, and so has the first call into the suite's own
    // non-inline code. The probe's own allocations are initialization odr-uses, which the same rule excludes,
    // so there is no circularity. (Observed to hold on Apple clang, gcc 14, MSVC 2022 and em++ 6.0.9 — but
    // the guarantee is the standard's, not the four rows'.)
    inline const bool installed = []
    {
        probe();
        verify();
        announce();
        return true;
    }();
}

//==============================================================================
// THE REPLACEMENTS. Eight `new`s, twelve `delete`s — every form the language can route an allocation
// through. The deletes are PAIRED with the news by allocator: the over-aligned ones came from
// _aligned_malloc / posix_memalign and are released the matching way, because std::free is not valid for
// the first of those.
void* operator new      (std::size_t s)                                              { return felitronics::test::alloc::release (s, 0, false); }
void* operator new[]    (std::size_t s)                                              { return felitronics::test::alloc::release (s, 0, false); }
void* operator new      (std::size_t s, std::align_val_t a)                          { return felitronics::test::alloc::release (s, (std::size_t) a, false); }
void* operator new[]    (std::size_t s, std::align_val_t a)                          { return felitronics::test::alloc::release (s, (std::size_t) a, false); }
void* operator new      (std::size_t s, const std::nothrow_t&) noexcept              { return felitronics::test::alloc::release (s, 0, true); }
void* operator new[]    (std::size_t s, const std::nothrow_t&) noexcept              { return felitronics::test::alloc::release (s, 0, true); }
void* operator new      (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return felitronics::test::alloc::release (s, (std::size_t) a, true); }
void* operator new[]    (std::size_t s, std::align_val_t a, const std::nothrow_t&) noexcept { return felitronics::test::alloc::release (s, (std::size_t) a, true); }

void  operator delete   (void* p)                                        noexcept { std::free (p); }
void  operator delete[] (void* p)                                        noexcept { std::free (p); }
void  operator delete   (void* p, std::size_t)                           noexcept { std::free (p); }
void  operator delete[] (void* p, std::size_t)                           noexcept { std::free (p); }
void  operator delete   (void* p, const std::nothrow_t&)                 noexcept { std::free (p); }
void  operator delete[] (void* p, const std::nothrow_t&)                 noexcept { std::free (p); }
void  operator delete   (void* p, std::align_val_t)                      noexcept { felitronics::test::alloc::alignedFree (p); }
void  operator delete[] (void* p, std::align_val_t)                      noexcept { felitronics::test::alloc::alignedFree (p); }
void  operator delete   (void* p, std::size_t, std::align_val_t)         noexcept { felitronics::test::alloc::alignedFree (p); }
void  operator delete[] (void* p, std::size_t, std::align_val_t)         noexcept { felitronics::test::alloc::alignedFree (p); }
void  operator delete   (void* p, std::align_val_t, const std::nothrow_t&) noexcept { felitronics::test::alloc::alignedFree (p); }
void  operator delete[] (void* p, std::align_val_t, const std::nothrow_t&) noexcept { felitronics::test::alloc::alignedFree (p); }

// The suites spell the numbers `alloc::count` / `alloc::bytes`. Declared here rather than in sixty-one
// files, for the same reason the counter itself is. It does claim the global name `alloc`, so a test that
// wanted its own `namespace alloc` would have to pick another — no such name exists in this tree, and the
// trade is one identifier against sixty-one repeated lines.
namespace alloc = felitronics::test::alloc;
