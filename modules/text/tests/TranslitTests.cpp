// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// JUCE-free self-tests for felitronics::text — the UTF-8 → ASCII romanization used to slug filenames.
// FALSIFICATION intent: try to prove a German umlaut / ß or a Cyrillic letter gets corrupted, dropped,
// or mis-decoded. Every mapping is asserted; malformed UTF-8 must be skipped, not guessed; distinct
// letters that intentionally fold to one slug are pinned so a caller knows the collision is real.
// (Source is UTF-8; U'…' char32_t literals ARE the code points, "…" literals ARE the bytes.)

#include <felitronics_test.h>
#include <felitronics/text/Translit.h>

#include <string>
#include <vector>

using felitronics::test::ok;
using felitronics::test::group;
using namespace felitronics::text;

int main()
{
    std::printf ("felitronics::text tests\n");

    group ("decodeUtf8 — widths, resync, and malformed bytes");
    {
        ok (decodeUtf8 ("") .empty(), "empty string decodes to no code points");
        ok (decodeUtf8 ("abc") == std::vector<char32_t> { U'a', U'b', U'c' }, "ASCII is 1 byte each");
        ok (decodeUtf8 ("ä") == std::vector<char32_t> { U'ä' }, "2-byte sequence (U+00E4)");
        ok (decodeUtf8 ("к") == std::vector<char32_t> { U'к' }, "2-byte Cyrillic (U+043A)");
        ok (decodeUtf8 ("€") == std::vector<char32_t> { U'€' }, "3-byte sequence (U+20AC)");
        ok (decodeUtf8 ("😀") == std::vector<char32_t> { U'😀' }, "4-byte sequence (U+1F600)");
        ok (decodeUtf8 (std::string ("a\xff" "b")) == std::vector<char32_t> { U'a', U'b' },
            "a stray 0xFF lead byte is skipped, the decoder resyncs");
        ok (decodeUtf8 (std::string ("a\x80" "b")) == std::vector<char32_t> { U'a', U'b' },
            "a stray continuation byte is skipped");
        ok (decodeUtf8 (std::string ("a\xc3")) == std::vector<char32_t> { U'a' },
            "a truncated 2-byte tail is dropped, not half-decoded");
    }

    group ("romanize — German + Latin-1");
    {
        ok (romanize (U'ä') == "ae" && romanize (U'Ä') == "Ae", "ä/Ä → ae/Ae");
        ok (romanize (U'ö') == "oe" && romanize (U'Ö') == "Oe", "ö/Ö → oe/Oe");
        ok (romanize (U'ü') == "ue" && romanize (U'Ü') == "Ue", "ü/Ü → ue/Ue");
        ok (romanize (U'ß') == "ss" && romanize (U'ẞ') == "SS", "ß/ẞ → ss/SS");
        ok (romanize (U'œ') == "oe" && romanize (U'Œ') == "OE", "œ/Œ → oe/OE");
        ok (romanize (U'é') == "e" && romanize (U'ñ') == "n" && romanize (U'ç') == "c",
            "common Latin-1 accents fold to the base letter");
        ok (romanize (U'å') == "a" && romanize (U'ø') == "o", "å/ø fold to a/o");
    }

    group ("romanize — full Russian Cyrillic, incl. the multi-letter cases");
    {
        ok (romanize (U'ж') == "zh" && romanize (U'Ж') == "Zh", "ж/Ж → zh/Zh");
        ok (romanize (U'х') == "kh" && romanize (U'Х') == "Kh", "х/Х → kh/Kh");
        ok (romanize (U'ц') == "ts" && romanize (U'Ц') == "Ts", "ц/Ц → ts/Ts");
        ok (romanize (U'ч') == "ch" && romanize (U'Ч') == "Ch", "ч/Ч → ch/Ch");
        ok (romanize (U'ш') == "sh" && romanize (U'Ш') == "Sh", "ш/Ш → sh/Sh");
        ok (romanize (U'щ') == "shch" && romanize (U'Щ') == "Shch", "щ/Щ → shch/Shch");
        ok (romanize (U'ё') == "yo" && romanize (U'Ё') == "Yo", "ё/Ё → yo/Yo");
        ok (romanize (U'ю') == "yu" && romanize (U'Ю') == "Yu", "ю/Ю → yu/Yu");
        ok (romanize (U'я') == "ya" && romanize (U'Я') == "Ya", "я/Я → ya/Ya");
        ok (romanize (U'э') == "e" && romanize (U'ы') == "y", "э→e, ы→y");
        ok (romanize (U'ъ').empty() && romanize (U'ь').empty(), "hard/soft signs romanize to nothing");
        ok (romanize (U'а') == "a" && romanize (U'я') != "" && romanize (U'п') == "p",
            "the contiguous block boundaries (а…я) are covered, none fall through");
    }

    group ("romanize — ASCII passthrough and out-of-scope drop");
    {
        ok (romanize (U'A') == "A" && romanize (U'z') == "z" && romanize (U'7') == "7",
            "ASCII alphanumerics pass through unchanged");
        ok (romanize (U'-').empty() && romanize (U' ').empty() && romanize (U'/').empty(),
            "ASCII punctuation is not slug material → \"\"");
        ok (romanize (U'😀').empty() && romanize (U'中').empty(),
            "out-of-scope scripts romanize to nothing (caller detects the empty)");
    }

    group ("transliterate — whole words fold to ASCII slugs");
    {
        ok (transliterate ("зелёный") == "zelyonyy", "зелёный → zelyonyy (ё→yo)");
        ok (transliterate ("красный") == "krasnyy", "красный → krasnyy");
        ok (transliterate ("Усилитель") == "Usilitel", "Усилитель → Usilitel (soft sign vanishes)");
        ok (transliterate ("Straße") == "Strasse" && transliterate ("Grün") == "Gruen"
            && transliterate ("Über") == "Ueber", "German words romanize");
        ok (transliterate ("Röhrenverstärker") == "Roehrenverstaerker", "long compound with two umlauts");
        ok (transliterate ("ёж") == "yozh", "leading ё");
        ok (transliterate ("a😀b") == "ab" && transliterate ("hi-cut 2") == "hicut2",
            "out-of-scope + ASCII punctuation are dropped, alnum survives");
        ok (transliterate ("ReVolt12h") == "ReVolt12h", "an all-ASCII slug is idempotent");
    }

    group ("transliterate — intentional collisions are real (caller must check)");
    {
        ok (transliterate ("е") == transliterate ("э") && transliterate ("е") == "e",
            "Cyrillic е and э both fold to 'e' — a documented, detectable collision");
        ok (transliterate ("й") == transliterate ("ы") && transliterate ("й") == "y",
            "й and ы both fold to 'y'");
    }

    return felitronics::test::report();
}
