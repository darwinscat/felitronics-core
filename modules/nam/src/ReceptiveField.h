// SPDX-License-Identifier: MIT
#pragma once
// HOW MUCH SIGNAL A MODEL NEEDS BEFORE ITS OUTPUT MEANS ANYTHING. A network started with empty
// buffers describes the silence it was born into for as long as its memory reaches back, so anything
// that fades a freshly loaded model in has to wait this long first.
//
// NAM answers this itself for convnet and lstm, and a container forwards to its active submodel —
// but WaveNet never overrides GetPrewarmSamples(), so for every WaveNet capture (which is nearly all
// of them) the honest upstream answer is zero. Hence this: the receptive field read straight off the
// config, one dilated convolution at a time.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>

namespace felitronics::nam::detail
{
inline int receptiveFieldFromConfig (const nlohmann::json& model);

inline int receptiveFieldOfLayers (const nlohmann::json& cfg)
{
    if (! cfg.contains ("layers") || ! cfg["layers"].is_array()) return 0;
    long long total = 0;
    for (const auto& grp : cfg["layers"])
    {
        if (! grp.contains ("kernel_sizes") || ! grp.contains ("dilations")) continue;
        const auto& ks = grp["kernel_sizes"];
        const auto& ds = grp["dilations"];
        if (! ks.is_array() || ! ds.is_array()) continue;
        for (std::size_t i = 0; i < ks.size() && i < ds.size(); ++i)
            total += (long long) ds[i].get<long long>() * ((long long) ks[i].get<long long>() - 1);
    }
    return total > 0 ? (int) std::min (total + 1, (long long) 1 << 20) : 0;
}

inline int receptiveFieldFromConfig (const nlohmann::json& model)
{
    if (! model.contains ("config")) return 0;
    const auto& cfg = model["config"];
    // A container holds several models and switches between them by level; any of them can be the
    // one speaking, so the longest memory is the one that has to be waited out.
    if (cfg.contains ("submodels") && cfg["submodels"].is_array())
    {
        int worst = 0;
        for (const auto& sub : cfg["submodels"])
            if (sub.contains ("model"))
                worst = std::max (worst, receptiveFieldFromConfig (sub["model"]));
        return worst;
    }
    return receptiveFieldOfLayers (cfg);
}
} // namespace felitronics::nam::detail
