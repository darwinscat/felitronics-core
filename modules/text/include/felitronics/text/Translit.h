// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

#pragma once

//==============================================================================
// felitronics::text — JUCE-free UTF-8 → ASCII romanization for building filesystem-safe, ASCII-only
// slugs out of user text that may contain German umlauts / ß or Russian Cyrillic. The RAW text should
// still be stored elsewhere (metadata) — this is for FILENAMES and tokens, where cross-platform sync
// (macOS NFD ↔ Windows/Linux NFC) and legacy ASCII parsers demand plain ASCII and a never-empty result.
//
//   decodeUtf8(s)     → the code points of s (malformed bytes skipped, never guessed)
//   romanize(cp)      → one code point's ASCII slug piece: ASCII alnum passes through; German umlauts,
//                       ß/ẞ/Œ, common Latin-1 accents, and Russian Cyrillic romanize; else "".
//   transliterate(s)  → romanize() folded over every code point of s (a slug of ASCII alnum letters).
//
// Romanization is intentionally lossy and opinionated (BGN-ish Russian, ae/oe/ue/ss German). Distinct
// inputs can fold to the same slug (е and э both → "e"); the caller is expected to detect collisions.
// std-only, header-only, no allocations beyond the returned strings. Message-thread (allocates).
//==============================================================================

#include <cctype>
#include <string>
#include <string_view>
#include <vector>

namespace felitronics::text
{

// Minimal UTF-8 → code-point decoder. A stray lead/continuation byte or a truncated tail is skipped
// (the decoder resynchronises), never turned into a guessed character.
inline std::vector<char32_t> decodeUtf8 (std::string_view s)
{
    std::vector<char32_t> cps;
    for (std::size_t i = 0, n = s.size(); i < n;)
    {
        const unsigned char c = static_cast<unsigned char> (s[i]);
        char32_t cp;
        int len;
        if      (c < 0x80)         { cp = c;                                len = 1; }
        else if ((c >> 5) == 0x06) { cp = static_cast<char32_t> (c & 0x1F); len = 2; }
        else if ((c >> 4) == 0x0E) { cp = static_cast<char32_t> (c & 0x0F); len = 3; }
        else if ((c >> 3) == 0x1E) { cp = static_cast<char32_t> (c & 0x07); len = 4; }
        else { ++i; continue; }                                     // stray lead / continuation byte
        if (i + static_cast<std::size_t> (len) > n) break;          // truncated tail
        bool okSeq = true;
        for (int k = 1; k < len; ++k)
        {
            const unsigned char cc = static_cast<unsigned char> (s[i + static_cast<std::size_t> (k)]);
            if ((cc >> 6) != 0x02) { okSeq = false; break; }
            cp = (cp << 6) | static_cast<char32_t> (cc & 0x3F);
        }
        if (! okSeq) { ++i; continue; }
        cps.push_back (cp);
        i += static_cast<std::size_t> (len);
    }
    return cps;
}

// One code point → its ASCII slug piece. ASCII alnum passes through; German digraphs, Latin-1 accents,
// and Russian Cyrillic romanize; anything else (other scripts, ASCII punctuation) yields "".
inline std::string romanize (char32_t cp)
{
    if (cp < 0x80)
        return std::isalnum (static_cast<unsigned char> (cp)) ? std::string (1, static_cast<char> (cp))
                                                              : std::string();

    switch (cp)                                      // German digraphs take priority over accent-fold
    {
        case 0x00C4: return "Ae"; case 0x00E4: return "ae";   // Ä ä
        case 0x00D6: return "Oe"; case 0x00F6: return "oe";   // Ö ö
        case 0x00DC: return "Ue"; case 0x00FC: return "ue";   // Ü ü
        case 0x00DF: return "ss"; case 0x1E9E: return "SS";   // ß ẞ
        case 0x0152: return "OE"; case 0x0153: return "oe";   // Œ œ
        default: break;
    }
    switch (cp)                                      // Latin-1 accented letters → base letter
    {
        case 0x00C0: case 0x00C1: case 0x00C2: case 0x00C3: case 0x00C5: return "A";
        case 0x00E0: case 0x00E1: case 0x00E2: case 0x00E3: case 0x00E5: return "a";
        case 0x00C7: return "C"; case 0x00E7: return "c";
        case 0x00C8: case 0x00C9: case 0x00CA: case 0x00CB: return "E";
        case 0x00E8: case 0x00E9: case 0x00EA: case 0x00EB: return "e";
        case 0x00CC: case 0x00CD: case 0x00CE: case 0x00CF: return "I";
        case 0x00EC: case 0x00ED: case 0x00EE: case 0x00EF: return "i";
        case 0x00D1: return "N"; case 0x00F1: return "n";
        case 0x00D2: case 0x00D3: case 0x00D4: case 0x00D5: case 0x00D8: return "O";
        case 0x00F2: case 0x00F3: case 0x00F4: case 0x00F5: case 0x00F8: return "o";
        case 0x00D9: case 0x00DA: case 0x00DB: return "U";
        case 0x00F9: case 0x00FA: case 0x00FB: return "u";
        case 0x00DD: return "Y"; case 0x00FD: case 0x00FF: return "y";
        default: break;
    }
    if (cp == U'Ё') return "Yo";                     // Ё — outside the contiguous Cyrillic block
    if (cp == U'ё') return "yo";                     // ё
    static const char* const kUpper[] = {            // А…Я (U+0410…U+042F)
        "A","B","V","G","D","E","Zh","Z","I","Y","K","L","M","N","O","P",
        "R","S","T","U","F","Kh","Ts","Ch","Sh","Shch","","Y","","E","Yu","Ya" };
    static const char* const kLower[] = {            // а…я (U+0430…U+044F)
        "a","b","v","g","d","e","zh","z","i","y","k","l","m","n","o","p",
        "r","s","t","u","f","kh","ts","ch","sh","shch","","y","","e","yu","ya" };
    if (cp >= U'А' && cp <= U'Я') return kUpper[cp - U'А'];
    if (cp >= U'а' && cp <= U'я') return kLower[cp - U'а'];
    return {};                                       // out of scope — dropped
}

// romanize() folded over the whole string: an ASCII-alnum slug of the input's letters and digits.
inline std::string transliterate (std::string_view s)
{
    std::string out;
    for (char32_t cp : decodeUtf8 (s)) out += romanize (cp);
    return out;
}

} // namespace felitronics::text
