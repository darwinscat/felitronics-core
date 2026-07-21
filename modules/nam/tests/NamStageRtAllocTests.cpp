// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Eigen's own runtime malloc gate complements the house operator-new counter: Eigen can call its
// aligned allocator directly, so only this guard sees every dynamic Eigen allocation in NAM process().

#include <felitronics_test.h>
#include <felitronics/nam/NamStage.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
std::string linearModel()
{
    return R"({"version":"0.5.0","architecture":"Linear","config":{"receptive_field":1,"bias":false,"implementation":"direct"},"weights":[0.75],"sample_rate":48000})";
}

std::string waveNetModel()
{
    // One one-channel, one-layer stack. Weight order is rechannel; dilated conv+bias; condition
    // mixin; residual 1x1+bias; head rechannel; head scale (nine scalars total at the pinned SHA).
    return R"({"version":"0.5.0","architecture":"WaveNet","config":{"layers":[{"input_size":1,"condition_size":1,"head_size":1,"head_bias":false,"channels":1,"kernel_size":2,"dilations":[1],"activation":"Tanh","gated":false}],"head_scale":1.0},"weights":[1,0,0,0,1,0,0,1,1],"sample_rate":48000})";
}

std::string lstmModel()
{
    // A one-layer, input=hidden=1 LSTM consumes 8 matrix, 4 bias, one initial hidden, one
    // initial cell, one head-weight, and one head-bias scalar: sixteen deterministic zeros.
    return R"({"version":"0.5.0","architecture":"LSTM","config":{"num_layers":1,"input_size":1,"hidden_size":1},"weights":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"sample_rate":48000})";
}

std::string convNetModel()
{
    // One 1x1-channel block: two convolution taps + bias, then one head weight + bias.
    return R"({"version":"0.5.0","architecture":"ConvNet","config":{"channels":1,"dilations":[1],"batchnorm":false,"activation":"Tanh"},"weights":[0,0,0,1,0],"sample_rate":48000})";
}

bool finite (const std::vector<float>& values)
{
    for (float value : values)
        if (! std::isfinite (value))
            return false;
    return true;
}

void exerciseNoEigenMalloc (const std::string& json, const char* architecture)
{
    felitronics::nam::NamStage stage;
    stage.prepare (48000.0, 64);
    const bool loaded = stage.loadModelFromMemory (json.data(), json.size());
    felitronics::test::ok (loaded, std::string (architecture) + " minimal fixture loads");
    if (! loaded)
        return;

    std::vector<float> left (64), right (64);
    for (int i = 0; i < 64; ++i)
    {
        left[(std::size_t) i] = 0.2f * std::sin (0.09f * (float) i);
        right[(std::size_t) i] = -0.15f * std::cos (0.07f * (float) i);
    }
    float* io[2] { left.data(), right.data() };
    stage.process (io, 2, 64, false);       // warm every process-reachable path before closing the gate

#if defined(EIGEN_RUNTIME_NO_MALLOC)
    std::printf ("    Eigen malloc gate: %s\n", architecture);
    std::fflush (stdout);
    Eigen::internal::set_is_malloc_allowed (false);
    stage.process (io, 2, 64, false);
    stage.process (io, 2, 64, true);
    Eigen::internal::set_is_malloc_allowed (true);
    std::printf ("    Eigen malloc gate passed: %s\n", architecture);
    std::fflush (stdout);
    felitronics::test::ok (finite (left) && finite (right),
                           std::string (architecture) + " stays finite with Eigen malloc forbidden");
#else
    // The target defines this macro, but keep the source honest if the test is reused standalone.
    felitronics::test::ok (false, std::string (architecture) + " Eigen runtime malloc guard is enabled");
#endif
}
} // namespace

int main()
{
    using namespace felitronics;
    std::printf ("felitronics::nam Eigen RT-allocation tests\n");

    test::group ("pinned architectures that preallocate all process work");
    exerciseNoEigenMalloc (linearModel(), "Linear");
    exerciseNoEigenMalloc (waveNetModel(), "WaveNet");

    test::group ("upstream-allocating architectures remain load-compatible");
    {
        nam::NamStage lstm;
        lstm.prepare (48000.0, 64);
        const auto json = lstmModel();
        test::ok (lstm.loadModelFromMemory (json.data(), json.size()), "minimal LSTM fixture loads");
        std::vector<float> block (64, 0.125f);
        float* io[1] { block.data() };
        lstm.process (io, 1, 64, false);
        test::ok (finite (block), "minimal LSTM fixture processes finitely with allocation allowed");
        // Deliberately excluded from the closed Eigen-malloc gate: the pinned LSTM returns an
        // owning dynamic Eigen::VectorXf from get_hidden_state() for every processed sample.
    }
    {
        nam::NamStage stage;
        stage.prepare (48000.0, 64);
        const auto json = convNetModel();
        test::ok (stage.loadModelFromMemory (json.data(), json.size()), "minimal ConvNet fixture loads");
        // Deliberately excluded from the closed Eigen-malloc gate: ConvNet::process constructs
        // dynamic Eigen matrices at the pinned upstream SHA. OrbitCab inherits the same behavior;
        // it is an upstream issue, not a port regression, and load remains compatible by design.
    }

    return test::report();
}
