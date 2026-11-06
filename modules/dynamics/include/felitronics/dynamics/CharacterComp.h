// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>

namespace felitronics::dynamics
{

/** COMP'S OWN CHARACTERS (the owner, 01.11: «OPTO, FET и SQUASH» — ours, not copies of anybody's): one mono compressor
    whose curve, detector, release and gain element a character sets (Voice). Moved here from OrbitAmp
    (src/core/CompCharacters.h), its arithmetic unchanged — the bass amp's OPTO, FET and SQUASH.

      - THE CURVE hears the input (feed-forward) or the squeezed signal (FEEDBACK). Fed back, holding a ratio R takes
        R - 1 dB of reduction for every dB the OUTPUT stands over the threshold, so the static curve keeps the ratio far
        above the threshold and its knee is R times wider in the input's terms: it eases in, and works slower near the
        threshold (staticDeltaDb draws it). The loop's steps are divided by how much it gains where it works (1 under the
        knee, R above it), so ATTACK and RELEASE are what they say, squeezing or let go.
      - THE DETECTOR: the level's RMS over 10 ms, or its peaks.
      - THE RELEASE in one stage, or two (a MEMORY): a fast stage back to what the memory holds, the memory sinking in the
        longer the squeeze is held (chargeMs) and letting go slowly (memoryMs) — a short hit lets go at once, a long
        squeeze slowly.
      - THE GAIN ELEMENT'S COLOUR, growing with the squeeze (full at fullColourDb): GRAIN — odd harmonics, worked at twice
        the rate (a 31-tap half-band each way) and added as a residual without a fundamental of its own, so the note's
        level and phase stay the squeeze's; DARK — the top let down toward a one-pole low-pass at darkHz.

    Real-time safe, nothing allocated anywhere. */
class CharacterComp
{
public:
    struct Voice
    {
        bool   rms       = false;
        bool   feedback  = false;
        double attackMs  = 10.0;
        double releaseMs = 100.0;   // the (first) stage
        double memoryMs  = 0.0;     // the memory's let-go; 0 = one stage
        double chargeMs  = 800.0;   // how long a squeeze takes to sink into the memory
        double grain     = 0.0;     // at the full squeeze the third harmonic is about grain / 4 of the signal
        double darkHz    = 0.0;     // the low-pass the top is let down toward at the full squeeze; 0 = none
    };

    static constexpr double rmsMs = 10.0, rangeDb = 60.0, fullColourDb = 12.0, grainEnvelopeMs = 50.0;

    /** THE STATIC CURVE, dB (0 or below) at a steady level: feed-forward, core's GainComputer's (the soft knee centred on
        the threshold); fed back, out = in - (R - 1) * kneeOver (out - T) solved for the output. */
    static double staticDeltaDb (bool feedback, double levelDb, double thresholdDb, double ratio, double kneeDb) noexcept
    {
        const double d = levelDb - thresholdDb, h = 0.5 * std::max (0.0, kneeDb), r = std::max (1.0, ratio);
        double delta = 0.0;
        if (! feedback)
            delta = -kneeOver (d, h) * (1.0 - 1.0 / r);
        else
        {
            double z = d;                          // out - T
            if (h <= 0.0)
                z = d > 0.0 ? d / r : d;
            else if (d >= r * h)
                z = d / r;
            else if (d > -h)
            {
                // In the knee: z + (R - 1)(z + h)^2 / (4h) = d — with w = z + h, a w^2 + w - (d + h) = 0.
                const double a = (r - 1.0) / (4.0 * h);
                const double w = a > 0.0 ? (-1.0 + std::sqrt (1.0 + 4.0 * a * (d + h))) / (2.0 * a) : d + h;
                z = w - h;
            }
            delta = z - d;
        }
        return std::clamp (std::isfinite (delta) ? delta : 0.0, -rangeDb, 0.0);
    }

    void prepare (double sampleRate) noexcept
    {
        fs      = sampleRate > 0.0 ? sampleRate : 48000.0;
        rmsK    = 1.0 - std::exp (-1.0 / (rmsMs * 0.001 * fs));
        envK    = 1.0 - std::exp (-1.0 / (grainEnvelopeMs * 0.001 * 2.0 * fs));
        designHalfBand();
        reset();
    }

    /** Silence: nothing held, nothing in the filters. */
    void reset() noexcept { seed (0.0); }

    /** A fresh start holding `reductionDb` (0 or more) — taken over from another detector without a step. */
    void seed (double reductionDb) noexcept
    {
        red  = std::max (0.0, reductionDb);
        mem  = 0.0;
        ms   = 0.0;
        env  = 0.0;
        lp   = 0.0;
        up.fill (0.0f);
        down.fill (0.0f);
    }

    void set (const Voice& v, double thresholdDb, double ratio, double kneeDb) noexcept
    {
        voice = v;
        thr   = thresholdDb;
        ratio_ = std::max (1.0, ratio);
        half  = 0.5 * std::max (0.0, kneeDb);
        atkK    = coeff (v.attackMs);
        relK    = coeff (v.releaseMs);
        memK    = v.memoryMs > 0.0 ? coeff (v.memoryMs) : 0.0;
        chargeK = coeff (v.chargeMs);
        darkK   = v.darkHz > 0.0 ? 1.0 - std::exp (-2.0 * 3.141592653589793 * v.darkHz / fs) : 0.0;
    }

    /** One stretch, in place: `key` is what the detector hears (the input, filtered as the product likes), `grDb` the
        squeeze per sample, dB (0 or below) — before any make-up, as core's compressor taps its own. */
    void process (float* io, const float* key, int n, float* grDb) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const double k = (double) key[i];
            double levelDb;
            if (voice.rms)
            {
                ms += rmsK * (k * k - ms);
                levelDb = 10.0 * std::log10 (std::max (ms, 1.0e-14));
            }
            else
                levelDb = 20.0 * std::log10 (std::max (std::abs (k), 1.0e-7));

            // The reduction the curve asks for, dB — fed back, of the key under the gain it got a sample ago; the loop's
            // gain there divides its steps (1 under the knee, R above it).
            double want, loop = 1.0;
            if (voice.feedback)
            {
                const double over = levelDb - red - thr;
                want = std::min (rangeDb, kneeOver (over, half) * (ratio_ - 1.0));
                loop = 1.0 + (ratio_ - 1.0) * kneeSlope (over, half);
            }
            else
                want = std::min (rangeDb, kneeOver (levelDb - thr, half) * (1.0 - 1.0 / ratio_));

            if (want > red)
                red += atkK / loop * (want - red);
            else
            {
                const double rest = memK > 0.0 ? std::max (want, mem) : want;   // the fast stage lets go down to the memory
                if (rest < red)
                    red += relK / loop * (rest - red);
            }
            // The memory: while the squeeze is asked for it sinks in toward what is held; asked less than it holds, it lets
            // go — slowly.
            if (memK > 0.0)
            {
                if (want >= mem)
                    mem += chargeK * std::max (0.0, red - mem);
                else
                    mem += memK / loop * (want - mem);
            }

            double s = (double) io[i] * std::exp (-red * 0.11512925464970228);   // ln 10 / 20
            const double c = std::clamp (red / fullColourDb, 0.0, 1.0);

            if (voice.grain > 0.0)
                s += grainResidual ((float) s, voice.grain * c);

            if (darkK > 0.0)
            {
                lp += darkK * (s - lp);
                s  += c * (lp - s);
            }

            io[i] = (float) s;
            if (grDb != nullptr)
                grDb[i] = (float) -red;
        }
    }

    /** The reduction held now, dB (0 or more). */
    double reductionDb() const noexcept { return red; }

private:
    static double kneeOver (double x, double h) noexcept
    {
        if (h <= 0.0) return x > 0.0 ? x : 0.0;
        if (x <= -h)  return 0.0;
        if (x >= h)   return x;
        return (x + h) * (x + h) / (4.0 * h);
    }

    static double kneeSlope (double x, double h) noexcept
    {
        if (h <= 0.0) return x > 0.0 ? 1.0 : 0.0;
        return std::clamp ((x + h) / (2.0 * h), 0.0, 1.0);
    }

    double coeff (double timeMs) const noexcept
    {
        return timeMs > 0.0 ? 1.0 - std::exp (-1.0 / (timeMs * 0.001 * fs)) : 1.0;
    }

    /** The odd harmonics the gain element adds, at twice the rate: the squeezed signal up, its residual
        -amount * e * (u^3 - 3u/4) over a peak envelope e (u = x / e, so a steady note's residual is its third harmonic
        alone, no fundamental), and the residual down again — late by the two filters' 15 samples, harmonics only. */
    float grainResidual (float x, double amount) noexcept
    {
        // Up: the zero-stuffed stream's two phases — the even one the half-band's even taps, the odd one its centre.
        std::copy_backward (up.begin(), up.end() - 1, up.end());
        up[0] = x;
        float even = 0.0f;
        for (int k = 0; k < (int) up.size(); ++k)
            even += hEven[(size_t) k] * up[(size_t) k];
        const float odd = 2.0f * h[(size_t) centre] * up[(size_t) centreBase];

        for (const float v : { 2.0f * even, odd })
        {
            const double a = std::abs ((double) v);
            env = a > env ? a : env + envK * (a - env);
            double r = 0.0;
            if (env > 1.0e-9 && amount > 0.0)
            {
                const double u = (double) v / env;
                r = -amount * env * (u * u * u - 0.75 * u);
            }
            std::copy_backward (down.begin(), down.end() - 1, down.end());
            down[0] = (float) r;
        }

        float y = 0.0f;
        for (int k = 0; k < taps; ++k)
            y += h[(size_t) k] * down[(size_t) k];
        return y;
    }

    void designHalfBand() noexcept
    {
        // A Blackman-windowed half-band, each phase summing to one half: DC passes both phases alike.
        double sumEven = 0.0;
        for (int n = 0; n < taps; ++n)
        {
            const double t = (double) (n - centre) * 0.5;
            const double sinc = n == centre ? 1.0 : std::sin (3.141592653589793 * t) / (3.141592653589793 * t);
            const double w = 0.42 - 0.5 * std::cos (2.0 * 3.141592653589793 * n / (taps - 1))
                           + 0.08 * std::cos (4.0 * 3.141592653589793 * n / (taps - 1));
            hd[(size_t) n] = (n % 2 == 0) ? 0.5 * sinc * w : 0.0;
            sumEven += hd[(size_t) n];
        }
        for (int n = 0; n < taps; ++n)
            h[(size_t) n] = (float) (n == centre ? 0.5 : hd[(size_t) n] * 0.5 / sumEven);
        for (int k = 0; k < (int) hEven.size(); ++k)
            hEven[(size_t) k] = h[(size_t) (2 * k)];
    }

    static constexpr int taps = 31, centre = 15, centreBase = 7;   // the centre tap is odd: the odd phase is a pure delay

    Voice  voice;
    double fs = 48000.0, rmsK = 0.04, envK = 0.0002;
    double thr = -24.0, ratio_ = 4.0, half = 4.0;
    double atkK = 0.01, relK = 0.001, memK = 0.0, chargeK = 0.0001, darkK = 0.0;
    double red = 0.0, mem = 0.0, ms = 0.0, env = 0.0, lp = 0.0;
    std::array<double, taps> hd {};
    std::array<float, taps> h {};
    std::array<float, 16> hEven {};
    std::array<float, 16> up {};
    std::array<float, taps> down {};
};

} // namespace felitronics::dynamics
