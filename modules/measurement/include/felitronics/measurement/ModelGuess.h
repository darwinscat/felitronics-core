// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::measurement — gear-model detection in a free-form file name. PROMOTED verbatim
// from OrbitCapture (model/MicSetModel.h, reuse audit N5) so every capture product resolves
// imported-file names against a catalog the same way instead of reinventing the heuristic.
//
// "YA MES 412 TRAD 906-1.wav" + a catalog holding "Sennheiser e906" → "Sennheiser e906".
//
// The rule (deliberately conservative — a wrong guess poisons metadata forever):
//   • A catalog model's FINGERPRINT is its last word ("e906", "SM57", "C414", "R-121").
//   • A file-name token matches it exactly (case/punctuation-blind), or — weaker — matches its
//     bare digits when they're 3+ long ("906" → e906; a lone "57" stays ambiguous).
//   • Exact beats digits; ANY ambiguity (two models claiming a token) yields NO guess —
//     never a silent wrong pick. Whole-token matches only ("4140" is not "414").
//   • NB an EXACT fingerprint hit has no length floor: a catalog entry "Shure 57" is matched by
//     the token "57" (it IS that model's whole fingerprint) — only the digits FALLBACK needs 3+.
//
// ASCII-only by design: tokens are [A-Za-z0-9] runs, everything else is a separator — enforced
// with locale-FREE classifiers (std::isalnum is locale-sensitive: under a single-byte locale a
// high-bit byte could join a token and break tokenizing; crew finding). Gear catalogs are ASCII
// in practice; non-ASCII bytes simply act as separators.
//==============================================================================

#include <string>
#include <vector>

namespace felitronics::measurement {

namespace detail {
    inline bool asciiAlnum(char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    }
    inline char asciiLower(char c) { return c >= 'A' && c <= 'Z' ? (char)(c - 'A' + 'a') : c; }
}

// Lower-cased ASCII-alphanumeric tokens of `s`; every other byte separates (locale-independent).
inline std::vector<std::string> alnumTokens(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (detail::asciiAlnum(c)) cur += detail::asciiLower(c);
        else if (!cur.empty()) { out.push_back(cur); cur.clear(); }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// The unique catalog model whose fingerprint appears in `fileName`, or "" (none / ambiguous).
inline std::string guessModel(const std::string& fileName, const std::vector<std::string>& catalog) {
    const auto toks = alnumTokens(fileName);
    std::string exact, digits;
    int exactHits = 0, digitHits = 0;
    for (const auto& model : catalog) {
        const auto mt = alnumTokens(model);
        if (mt.empty()) continue;
        const std::string& fp = mt.back();                            // the fingerprint token
        std::string num;
        for (char c : fp) if (c >= '0' && c <= '9') num += c;
        bool hitExact = false, hitNum = false;
        for (const auto& t : toks) {
            if (t == fp) hitExact = true;
            else if (num.size() >= 3 && t == num) hitNum = true;
        }
        if (hitExact && exact != model) { exact = model; ++exactHits; }
        else if (hitNum && digits != model) { digits = model; ++digitHits; }
    }
    if (exactHits == 1) return exact;
    if (exactHits == 0 && digitHits == 1) return digits;
    return {};                                                         // none or ambiguous
}

} // namespace felitronics::measurement
