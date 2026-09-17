// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>

namespace felitronics::mastering
{

struct SolvePassRecord;

enum class ProgressStage
{
    Convert       = 0,
    LoudnessRange = 1,
    SearchPass    = 2,
    FinalRender   = 3,
    Render        = 4
};

struct ProgressEvent
{
    ProgressStage stage = ProgressStage::Convert;
    int    pass = 0;
    int    maxPasses = 0;
    double fraction = 0.0;
    const SolvePassRecord* record = nullptr;
};

struct ProgressCallback
{
    bool (*fn) (void* context, const ProgressEvent& event) = nullptr;
    void* context = nullptr;
};

class ProgressClock
{
public:
    explicit ProgressClock (const ProgressCallback& callback) noexcept : cb_ (callback) {}

    bool stopped() const noexcept { return stopped_; }

    [[nodiscard]] bool begin (ProgressStage stage, int pass, int maxPasses, long long units, long long frames) noexcept
    {
        ev_ = ProgressEvent {};
        ev_.stage = stage; ev_.pass = pass; ev_.maxPasses = maxPasses;
        units_ = std::max (0LL, units);
        done_ = 0; last_ = 0;
        step_ = std::max (1LL, frames / 200);
        return emit (0.0, nullptr);
    }

    int piece (long long want) const noexcept
    {
        return (int) (cb_.fn == nullptr ? want : std::min (want, step_));
    }

    [[nodiscard]] bool advance (long long n) noexcept
    {
        done_ += n;
        if (cb_.fn == nullptr || stopped_ || done_ - last_ < step_ || done_ >= units_) return ! stopped_;
        last_ = done_;
        return emit ((double) done_ / (double) units_, nullptr);
    }

    [[nodiscard]] bool finish (const SolvePassRecord* record = nullptr) noexcept { return emit (1.0, record); }

private:
    bool emit (double fraction, const SolvePassRecord* record) noexcept
    {
        if (stopped_) return false;
        if (cb_.fn == nullptr) return true;
        ev_.fraction = fraction;
        ev_.record = record;
        stopped_ = ! cb_.fn (cb_.context, ev_);
        return ! stopped_;
    }

    ProgressCallback cb_;
    ProgressEvent    ev_ {};
    long long units_ = 0, done_ = 0, last_ = 0, step_ = 1;
    bool stopped_ = false;
};

} // namespace felitronics::mastering
