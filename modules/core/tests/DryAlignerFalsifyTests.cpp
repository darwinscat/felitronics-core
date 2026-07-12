// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Adversarial falsification pass for DryAligner: small hand-derived schedules from the external
// case list plus reset/re-prepare boundary attacks.

#include <felitronics_test.h>
#include <felitronics/core/DryAligner.h>

#include <algorithm>
#include <vector>

using felitronics::core::DryAligner;
using felitronics::test::group;
using felitronics::test::ok;

namespace
{
std::vector<float> runMono (DryAligner& a, const std::vector<float>& in, int delay, int block)
{
    std::vector<float> out (in.size(), 0.0f);
    std::vector<float> tmp ((std::size_t) block);
    for (std::size_t off = 0; off < in.size(); off += (std::size_t) block)
    {
        const int n = (int) std::min<std::size_t> ((std::size_t) block, in.size() - off);
        for (int i = 0; i < n; ++i) tmp[(std::size_t) i] = in[off + (std::size_t) i];
        const float* io[1] { tmp.data() };
        a.advance (io, 1, n, delay);
        for (int i = 0; i < n; ++i) out[off + (std::size_t) i] = a.delayed (0)[i];
    }
    return out;
}

bool equals (const std::vector<float>& a, const std::vector<float>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}
} // namespace

int main()
{
    std::printf ("felitronics::core DryAligner FALSIFY tests\n");

    group ("DeepSeek 3.1 materialized: D=0 stages the input exactly");
    {
        DryAligner a; a.prepare (1, 5, 8);
        const std::vector<float> in { 1.0f, -2.0f, 3.5f, 4.0f, -5.0f, 6.0f, 7.0f };
        ok (equals (runMono (a, in, 0, 3), in), "D=0 is bit-exact identity across ragged blocks");
    }

    group ("DeepSeek 3.2 materialized: D=capacity-1 gives D leading zeros, then input prefix");
    {
        constexpr int cap = 8, D = cap - 1;
        DryAligner a; a.prepare (1, 4, cap);
        std::vector<float> in (20), want (20, 0.0f);
        for (int i = 0; i < 20; ++i) in[(std::size_t) i] = (float) (i + 10);
        for (int i = D; i < 20; ++i) want[(std::size_t) i] = in[(std::size_t) (i - D)];
        ok (equals (runMono (a, in, D, 4), want), "capacity-1 tap has exact warm-up zeros and delayed prefix");
    }

    group ("DeepSeek 3.3 materialized: current block's variable D reads the shared warm cursor");
    {
        DryAligner a; a.prepare (1, 3, 8);
        std::vector<float> b1 { 1.0f, 2.0f, 3.0f };
        const float* io1[1] { b1.data() };
        a.advance (io1, 1, 3, 0);
        bool first = a.delayed (0)[0] == 1.0f && a.delayed (0)[1] == 2.0f && a.delayed (0)[2] == 3.0f;

        std::vector<float> b2 { 4.0f, 5.0f, 6.0f };
        const float* io2[1] { b2.data() };
        a.advance (io2, 1, 3, 2);
        bool second = a.delayed (0)[0] == 2.0f && a.delayed (0)[1] == 3.0f && a.delayed (0)[2] == 4.0f;
        ok (first && second, "[A,B,C] D=0 then [D,E,F] D=2 stages [B,C,D]");
    }

    group ("DeepSeek 3.4 materialized: long wraparound through cap=8, D=3");
    {
        DryAligner a; a.prepare (1, 5, 8);
        std::vector<float> in (20), want (20, 0.0f);
        for (int i = 0; i < 20; ++i) in[(std::size_t) i] = (float) i;
        for (int i = 3; i < 20; ++i) want[(std::size_t) i] = in[(std::size_t) (i - 3)];
        ok (equals (runMono (a, in, 3, 5), want), "cap=8 wrap emits out[n] = n<3 ? 0 : in[n-3]");
    }

    group ("own attack: re-prepare to a smaller ring clears old history and clamps to the new capacity");
    {
        DryAligner a; a.prepare (1, 8, 64);
        std::vector<float> warm (32, 99.0f);
        (void) runMono (a, warm, 12, 8);

        a.prepare (1, 4, 5);
        std::vector<float> in { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f };
        std::vector<float> want { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 2.0f }; // requested 99 clamps to D=4
        ok (equals (runMono (a, in, 99, 4), want), "re-prepare clears stale samples and applies the new cap-1 clamp");
    }

    return felitronics::test::report();
}
