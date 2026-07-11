// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// Gear-model detection in free-form file names (ModelGuess.h, promoted from OrbitCapture).
// Locks the conservative contract: exact fingerprint beats bare digits, whole tokens only,
// any ambiguity yields NO guess — a wrong guess poisons imported-take metadata forever.

#include <felitronics/measurement/ModelGuess.h>

#include "felitronics_test.h"

using felitronics::measurement::alnumTokens;
using felitronics::measurement::guessModel;
using felitronics::test::group;
using felitronics::test::ok;

int main()
{
    std::printf ("felitronics::measurement model-guess tests\n");

    group ("alnumTokens: lower-cased alnum runs, everything else separates");
    {
        const auto t = alnumTokens ("YA MES-412_TRAD 906-1.wav");
        ok (t.size() == 7 && t[0] == "ya" && t[2] == "412" && t[4] == "906" && t[6] == "wav",
            "tokens split on space/-/_/. and lower-case");
        ok (alnumTokens ("").empty(), "empty in, empty out");
        ok (alnumTokens ("---").empty(), "separators only -> no tokens");
    }

    group ("guessModel: known mics detected in imported file names (ported pins)");
    {
        const std::vector<std::string> cat { "Shure SM57", "Sennheiser e906", "Sennheiser e609",
                                             "AKG C414", "Royer R-121", "Beyerdynamic M160" };
        ok (guessModel ("YA MES 412 TRAD 906-1", cat) == "Sennheiser e906", "bare digits pick e906 (906-1)");
        ok (guessModel ("cab_SM57_capedge", cat) == "Shure SM57", "exact token, punctuation-blind");
        ok (guessModel ("Mesa e609 close", cat) == "Sennheiser e609", "e609 vs e906 stay distinct");
        ok (guessModel ("R-121 ribbon 2", cat) == "Royer R-121", "hyphenated fingerprint matches");
        ok (guessModel ("take 57", cat).empty(), "2-digit bare number is too ambiguous");
        ok (guessModel ("cab 4140 bright", cat).empty(), "'4140' is not '414' (whole-token only)");
        ok (guessModel ("sm57 or m160", cat).empty(), "two exact hits -> no guess (never a wrong pick)");
        ok (guessModel ("room mic", cat).empty(), "nothing recognizable -> empty");
    }

    group ("guessModel: contract edges");
    {
        const std::vector<std::string> cat { "Shure SM57", "Sennheiser e906" };
        ok (guessModel ("sm57 906 blend", cat) == "Shure SM57", "exact beats a simultaneous digit hit");
        ok (guessModel ("anything", {}).empty(), "empty catalog -> empty");
        ok (guessModel ("", cat).empty(), "empty file name -> empty");

        const std::vector<std::string> dup { "Shure SM57", "Shure SM57" };
        ok (guessModel ("sm57 take", dup) == "Shure SM57", "duplicate catalog entry is one model, not ambiguity");

        const std::vector<std::string> clash { "Shure SM57", "Fake SM57" };
        ok (guessModel ("sm57 take", clash).empty(), "two models sharing a fingerprint -> ambiguous -> empty");

        const std::vector<std::string> digitClash { "AKG C414", "Fake R414" };
        ok (guessModel ("mic 414 far", digitClash).empty(), "two models sharing bare digits -> ambiguous");

        const std::vector<std::string> odd { "   ", "Sennheiser e906" };
        ok (guessModel ("906 dry", odd) == "Sennheiser e906", "catalog entry with no tokens is skipped");

        const std::vector<std::string> bare { "Shure 57" };
        ok (guessModel ("take 57 close", bare) == "Shure 57",
            "a 2-digit token IS an exact hit when it is the model's whole fingerprint (doc'd)");
        ok (guessModel ("t\xc3\xa9st 906 caf\xc3\xa9", odd) == "Sennheiser e906",
            "UTF-8 bytes separate tokens (locale-free tokenizer)");
    }

    return felitronics::test::report();
}
