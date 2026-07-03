// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// HEADER HYGIENE — one TU that #includes EVERY public header (felitronics + the teq compat shims) and
// is compiled under a strict, downstream-grade warning set (-Wconversion, -Wfloat-equal, … + -Werror;
// see the target in the top-level CMakeLists). Consumers (TabbyEQ links its test targets against
// juce::juce_recommended_warning_flags) compile these headers under exactly this class of flags — this
// target makes header warning-cleanliness a CORE gate instead of a downstream surprise.
//
// Class templates are explicitly instantiated below (with their shipped default arguments) so warnings
// inside template member bodies fire HERE, not only when a consumer instantiates them.

#include <felitronics/analysis/CorrelationMeter.h>
#include <felitronics/analysis/KWeightingFilter.h>
#include <felitronics/analysis/LoudnessMeter.h>
#include <felitronics/analysis/SpectrumTap.h>
#include <felitronics/analysis/TruePeakMeter.h>
#include <felitronics/convolution/ConvolutionEngine.h>
#include <felitronics/convolution/IrResampler.h>
#include <felitronics/convolution/MatrixConvolver.h>
#include <felitronics/convolution/PartitionedConvolver.h>
#include <felitronics/core/Config.h>
#include <felitronics/core/DelayLine.h>
#include <felitronics/core/Fft.h>
#include <felitronics/core/FlushToZero.h>
#include <felitronics/core/Math.h>
#include <felitronics/core/Smoother.h>
#include <felitronics/deesser/DeEsser.h>
#include <felitronics/dither/Dither.h>
#include <felitronics/dynamiceq/DynamicEqBand.h>
#include <felitronics/dynamics/ChannelLinker.h>
#include <felitronics/dynamics/Compressor.h>
#include <felitronics/dynamics/EnvelopeFollower.h>
#include <felitronics/dynamics/GainComputer.h>
#include <felitronics/dynamics/GainReductionFollower.h>
#include <felitronics/dynamics/TransientShaper.h>
#include <felitronics/eq/Crossover2.h>
#include <felitronics/eq/EqBand.h>
#include <felitronics/eq/EqEngine.h>
#include <felitronics/eq/EqTypes.h>
#include <felitronics/eq/MatchedBiquad.h>
#include <felitronics/eq/MultibandSplitter.h>
#include <felitronics/eq/Svf.h>
#include <felitronics/limiter/TruePeakLimiter.h>
#include <felitronics/lineareq/LinearPhaseEq.h>
#include <felitronics/lineareq/MixedPhaseFir.h>
#include <felitronics/lineareq/NaturalPhaseEq.h>
#include <felitronics/multiband/MultibandCompressor.h>
#include <felitronics/multiband/MultibandProcessor.h>
#include <felitronics/multiband/MultibandWidth.h>
#include <felitronics/neural/Inference.h>
#include <felitronics/neural/NeuralStage.h>
#include <felitronics/oversampling/PolyphaseOversampler.h>
#include <felitronics/saturation/Saturator.h>
#include <felitronics/saturation/WaveShaper.h>
#include <felitronics/stereo/MidSide.h>
#include <felitronics/stereo/MonoBass.h>
#include <felitronics/stereo/StereoWidth.h>

// The teq compat shims are public too — they are what TabbyEQ/OrbitCab actually include.
#include <teq/EqBand.h>
#include <teq/EqEngine.h>
#include <teq/EqTypes.h>
#include <teq/MatchedBiquad.h>
#include <teq/Math.h>
#include <teq/Smoother.h>
#include <teq/SpectrumTap.h>
#include <teq/Svf.h>

namespace
{
    // Minimal Inference backend so NeuralStage's members instantiate.
    struct NullInference
    {
        void prepare (double, int, int) noexcept {}
        void process (float* const*, int, int) noexcept {}
        void reset() noexcept {}
        int  latencySamples() const noexcept { return 0; }
    };
    static_assert (felitronics::neural::Inference<NullInference>);
}

// Explicit instantiations (shipped defaults) — force every template member body through the compiler.
template class  felitronics::convolution::PartitionedConvolver<felitronics::core::fft::DefaultRealFft>;
template class  felitronics::convolution::ConvolutionEngine<felitronics::core::fft::DefaultRealFft, 2>;
template class  felitronics::convolution::MatrixConvolver<felitronics::core::fft::DefaultRealFft>;
template class  felitronics::lineareq::MixedPhaseFir<felitronics::core::fft::DefaultRealFft>;
template class  felitronics::eq::MultibandSplitter<4>;
template class  felitronics::multiband::MultibandProcessor<felitronics::dynamics::Compressor, 4>;
template class  felitronics::multiband::MultibandCompressor<4>;
template class  felitronics::multiband::MultibandWidth<4>;
template class  felitronics::neural::NeuralStage<NullInference>;
template struct felitronics::analysis::SpectrumTapT<felitronics::analysis::kSpectrumFftOrder>;

int main() { return 0; }
