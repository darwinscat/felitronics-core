// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <vector>

namespace felitronics::core
{

//==============================================================================
// felitronics::core::DelayLine — a fixed-size single-channel integer delay line (RT-safe ring). The
// generic building block for LOOKAHEAD (delay the signal by the lookahead while a detector reads ahead),
// PDC, and short delays. prepare() allocates; process()/reset are alloc-free. One per channel.
class DelayLine
{
public:
    void prepare (int maxDelaySamples)
    {
        cap_ = maxDelaySamples < 0 ? 0 : maxDelaySamples;
        buf_.assign ((std::size_t) cap_ + 1, 0.0f);   // +1 so a delay of `cap_` is addressable
        reset();
    }

    void reset() noexcept { std::fill (buf_.begin(), buf_.end(), 0.0f); pos_ = 0; }

    void setDelay (int d) noexcept { delay_ = d < 0 ? 0 : (d > cap_ ? cap_ : d); }
    int  delay()    const noexcept { return delay_; }
    int  capacity() const noexcept { return cap_; }

    // Push one sample, return the sample `delay()` samples ago (delay 0 = passthrough). RT-safe.
    inline float process (float x) noexcept
    {
        const int n = (int) buf_.size();
        buf_[(std::size_t) pos_] = x;
        int rd = pos_ - delay_;
        if (rd < 0) rd += n;
        const float y = buf_[(std::size_t) rd];
        if (++pos_ >= n) pos_ = 0;
        return y;
    }

private:
    std::vector<float> buf_;
    int cap_ = 0, delay_ = 0, pos_ = 0;
};

} // namespace felitronics::core
