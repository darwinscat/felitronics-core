// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// P88. RANDOM INTERLEAVINGS of publish / process (random width, random length) / reset() /
// clearAudioState() on the swap-safe convolvers, against the promise law 11e states: after reset(),
// the live operator is the one the caller last PUBLISHED, the history is empty and the convolver is Idle.
//
// The oracle is a SHADOW of the same class brought to that state the long way — prepared, handed the
// operator, run on silence until its fade is over, and reset at Idle, where nothing can be lost — and then
// fed the stream the convolver under test is fed, compared bit for bit until the next accepted publication
// (whose fade the shadow does not share). A shadow of the same class is only half an oracle, so for the
// engine a reference convolution from OUTSIDE the class is compared as well (the matrix siblings are nulled
// against one in their own files).
// Measured on the bodies before P88: red on the first seed of every class. The fixed groups in each class's
// own test file pin the named cells; this is what finds the sequence nobody wrote down.

#include <felitronics_test.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>
#include <felitronics/core/Math.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

using namespace felitronics;

namespace
{
struct Rng
{
    unsigned long long s;
    unsigned next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return (unsigned) (s >> 33); }
    int   in (int a, int b) { return a + (int) (next() % (unsigned) (b - a + 1)); }
    float f() { return (float) (next() & 0xffffu) / 32768.0f - 1.0f; }
};

constexpr int kP = 64, kIrMax = 1024, kIrLen = 700, kXf = 96, kNumOps = 4;

std::vector<std::vector<float>> makeOperators()
{
    Rng r { 99 };
    std::vector<std::vector<float>> ops ((std::size_t) kNumOps, std::vector<float> ((std::size_t) kIrLen));
    for (int k = 0; k < kNumOps; ++k)
    {
        for (auto& v : ops[(std::size_t) k]) v = 0.05f * r.f();
        ops[(std::size_t) k][0] = 0.2f * (float) (k + 1) * ((k & 1) != 0 ? -1.0f : 1.0f);   // distinct heads
    }
    return ops;
}

struct Engine
{
    convolution::ConvolutionEngine<> e;
    int nch = 2;
    bool prepare (int c) { nch = c; return e.prepare (kP, kIrMax, kXf, c); }
    bool publish (const std::vector<float>& ir) { return e.setIr (ir.data(), kIrLen); }
    bool process (const float* const* in, float* const* out, int nc, int n) { return e.process (in, out, nc, n); }
    void reset() { e.reset(); }
    void clear() { e.clearAudioState(); }
    bool busy() const { return e.isBusy(); }
    static constexpr bool kWidthFree = true;    // may be called narrower than prepared
};

template <class C> struct Matrix
{
    C e;
    int nch = 2;
    bool prepare (int c) { nch = c; return e.prepare (kP, kIrMax, kXf, c); }
    bool publish (const std::vector<float>& ir)
    {
        const float* b[2] { ir.data(), ir.data() };
        return e.setOperator (C::Topology::LRDiag, b, nch == 1 ? 1 : 2, kIrLen);
    }
    bool process (const float* const* in, float* const* out, int nc, int n) { return e.process (in, out, nc, n); }
    void reset() { e.reset(); }
    void clear() { e.clearAudioState(); }
    bool busy() const { return e.isBusy(); }
    static constexpr bool kWidthFree = false;   // the matrix width is exact (law 11b)
};

// The uniform reference answers the engine within 1e-6 (0 measured); the matrix siblings are nulled against it
// in their own files, and Nupc's non-uniform schedule only to a relative 2e-4 — the shadow carries those.
template <class A> bool exactRef() { return std::is_same_v<A, Engine>; }

template <class A>
void fuzz (const std::string& name, int seeds, int opsPerSeed)
{
    const auto ops = makeOperators();
    int failures = 0, windows = 0, busyAfterReset = 0, clearMovedBusy = 0;
    long long checked = 0;
    for (int seed = 1; seed <= seeds && failures == 0; ++seed)
    {
        Rng r { (unsigned long long) seed * 7919ULL };
        const int nch = r.in (1, 2);
        A dut;
        felitronics::test::ok (dut.prepare (nch), name + " prepare");
        int published = -1;                          // the operator the caller last had ACCEPTED
        std::unique_ptr<A> shadow;                   // armed by reset(), disarmed by the next accepted publish
        int shadowNc = -1;
        convolution::PartitionedConvolver<> ref;
        bool refOn = false;
        std::vector<float> xl (700), xr (700), yl (700), yr (700), sl (700), sr (700), rr (700);
        for (int op = 0; op < opsPerSeed && failures == 0; ++op)
        {
            const int what = r.in (0, 9);
            if (what <= 1)                                                   // publish
            {
                const int k = r.in (0, kNumOps - 1);
                if (dut.publish (ops[(std::size_t) k])) { published = k; shadow.reset(); refOn = false; }
            }
            else if (what == 2)                                              // restart
            {
                dut.reset();
                if (dut.busy()) ++busyAfterReset;
                shadow.reset(); refOn = false;
                if (published >= 0)
                {
                    shadow = std::make_unique<A>();
                    felitronics::test::run (shadow->prepare (nch));
                    felitronics::test::ok (shadow->publish (ops[(std::size_t) published]), name + " shadow publish");
                    std::vector<float> z (512, 0.0f), zo (512, 0.0f);
                    const float* zi[2] { z.data(), z.data() };
                    float*       zz[2] { zo.data(), zo.data() };
                    for (int i = 0; i < 4 && shadow->busy(); ++i) felitronics::test::run (shadow->process (zi, zz, nch, 512));
                    felitronics::test::ok (! shadow->busy(), name + " shadow settled");
                    shadow->reset();
                    shadowNc = -1;
                    if (exactRef<A>())
                    {
                        ref = convolution::PartitionedConvolver<>();
                        ref.prepare (kP, kIrMax);
                        ref.setIr (ops[(std::size_t) published].data(), kIrLen);
                        refOn = true;
                    }
                    ++windows;
                }
            }
            else if (what == 3)                                              // history only
            {
                const bool before = dut.busy();
                dut.clear();
                if (dut.busy() != before) ++clearMovedBusy;
                shadow.reset(); refOn = false;       // a fade may still be running: the shadow no longer applies
            }
            else                                                             // process
            {
                int nc = A::kWidthFree ? r.in (0, nch) : nch;
                if (shadow != nullptr && nc == 0) nc = nch;                  // a zero-width call is its own falling edge
                if (shadow != nullptr && shadowNc >= 0 && A::kWidthFree) nc = shadowNc;   // keep the window's width
                const int n = r.in (1, 700);
                for (int i = 0; i < n; ++i) { xl[(std::size_t) i] = 0.5f * r.f(); xr[(std::size_t) i] = 0.5f * r.f(); }
                const float* in[2]  { xl.data(), xr.data() };
                float*       out[2] { yl.data(), yr.data() };
                felitronics::test::run (dut.process (in, out, nc, n));
                if (shadow != nullptr && nc > 0)
                {
                    shadowNc = nc;
                    float* so[2] { sl.data(), sr.data() };
                    felitronics::test::run (shadow->process (in, so, nc, n));
                    for (int i = 0; i < n && failures == 0; ++i)
                        if (! core::sameBits (yl[(std::size_t) i], sl[(std::size_t) i]) || (nc == 2 && ! core::sameBits (yr[(std::size_t) i], sr[(std::size_t) i])))
                        {
                            felitronics::test::ok (false, name + ": seed " + std::to_string (seed) + " op " + std::to_string (op)
                                                          + " — after reset() the output left the settled shadow at +" + std::to_string (i));
                            ++failures;
                        }
                    if (refOn)
                    {
                        felitronics::test::run (ref.process (xl.data(), rr.data(), n));
                        for (int i = 0; i < n && failures == 0; ++i)
                            if (! (std::fabs (yl[(std::size_t) i] - rr[(std::size_t) i]) <= 1e-6f))
                            {
                                felitronics::test::ok (false, name + ": seed " + std::to_string (seed) + " op " + std::to_string (op)
                                                              + " — after reset() the output left the REFERENCE convolution at +" + std::to_string (i));
                                ++failures;
                            }
                    }
                    checked += n;
                }
            }
        }
    }
    std::printf ("      %s: %d post-reset windows, %lld samples compared\n", name.c_str(), windows, checked);
    felitronics::test::ok (failures == 0, name + ": every post-reset window matches the settled shadow bit for bit");
    felitronics::test::ok (busyAfterReset == 0, name + ": Idle after EVERY reset() (" + std::to_string (busyAfterReset) + " were not)");
    felitronics::test::ok (clearMovedBusy == 0, name + ": clearAudioState() never ends or starts a swap (" + std::to_string (clearMovedBusy) + " did)");
    // A fuzz that never reached its own oracle certifies nothing.
    felitronics::test::ok (windows > 10 * seeds / 2 && checked > 100000, name + ": the oracle was actually exercised");
}
} // namespace

int main()
{
    std::printf ("felitronics::convolution restart fuzz (P88)\n");
    felitronics::test::group ("ConvolutionEngine — random publish / process / reset / clearAudioState");
    fuzz<Engine> ("ConvolutionEngine", 120, 120);
    felitronics::test::group ("MatrixConvolver — random publish / process / reset / clearAudioState");
    fuzz<Matrix<convolution::MatrixConvolver<>>> ("MatrixConvolver", 120, 120);
    return felitronics::test::report();
}
