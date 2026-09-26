// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (c) 2026 Darwin's Cat — Oleh Tsymaienko & Alisa Lafoks. Part of felitronics-core — see LICENSE.

// felitronics::io and RIFF's word alignment. Every chunk of a RIFF file is padded to an even length with one zero
// byte that its own size field does NOT count, and the RIFF size DOES. The writer never wrote that byte, so a
// 24-bit data chunk of an odd length — an odd channel count times an odd frame count, 3 bytes a frame in mono —
// ended the file on an odd byte. Measured before the fix: 24-bit mono, 1 frame → a 47-byte file, RIFF 39, data 3;
// 1001 frames → 3047, RIFF 3039, data 3003. Forgiving readers (ffprobe, afinfo, ours) read those, but a chunk
// anything appends after one starts on an odd offset, where a reader that skips the pad — ours included — lands
// one byte into its header and rejects the file as corrupt.
//
//   (1) THE LAYOUT, across 16/24-bit PCM and 32-bit float × 1/2/3 channels × 1/2/3/1000/1001 frames: the file is
//       even, RIFF size = file − 8, data size = frames·blockAlign with the pad NOT counted, the pad is there iff
//       the data is odd and it is zero, and a strict chunk walk covers the file exactly; write∘read is the
//       identity on the samples and read∘write on the bytes;
//   (2) A CHUNK APPENDED after the image, as a tagger adds a LIST, reads: the audio is unchanged and the walk
//       finds the chunk;
//   (3) THE READER steps over a pad wherever it sits — an odd LIST before fmt, an odd chunk after an odd data
//       chunk, several in a row — and still reads the images the writer made before (odd data at EOF, no pad);
//   (4) THE PLANTED FAILURES: the old writer's layout (pad dropped, RIFF = 36 + data) fails (1)'s walk and (2)'s
//       append on exactly the odd cases and passes every even one; a walk that does NOT skip the pad loses every
//       image of (3) with an odd chunk before another — so a green run is the instruments seeing, not images
//       that would have read either way.

#include <felitronics_test.h>
#include <felitronics/io/Wav.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace felitronics;
using felitronics::test::ok;
using felitronics::test::group;

namespace
{
using Bytes = std::vector<std::uint8_t>;

std::uint32_t u32at (const Bytes& b, std::size_t p)
{
    return (std::uint32_t) b[p] | ((std::uint32_t) b[p + 1] << 8) | ((std::uint32_t) b[p + 2] << 16) | ((std::uint32_t) b[p + 3] << 24);
}
void putU32 (Bytes& b, std::size_t p, std::uint32_t v) { for (int i = 0; i < 4; ++i) b[p + i] = (std::uint8_t) ((v >> (8 * i)) & 0xFF); }
void tag (Bytes& o, const char* t) { o.insert (o.end(), t, t + 4); }
void u32 (Bytes& o, std::uint32_t v) { for (int i = 0; i < 4; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); }
void u16 (Bytes& o, std::uint16_t v) { for (int i = 0; i < 2; ++i) o.push_back ((std::uint8_t) ((v >> (8 * i)) & 0xFF)); }

// A chunk: its id, its size field, and a body of that many bytes (value b), then the pad when the size is odd.
void chunk (Bytes& o, const char* id, const Bytes& body, bool pad = true)
{
    tag (o, id); u32 (o, (std::uint32_t) body.size());
    o.insert (o.end(), body.begin(), body.end());
    if (pad && (body.size() & 1u)) o.push_back (0);
}
void fixRiff (Bytes& o) { putU32 (o, 4, (std::uint32_t) (o.size() - 8)); }

// RIFF read strictly: the container's size is the file's minus 8; every chunk header starts at an even offset;
// an odd chunk is followed by a ZERO pad its size does not count; the last padded chunk ends exactly at the end
// of the file. Returns why not (empty = conforming) and the ids it walked.
std::string conforms (const Bytes& img, std::vector<std::string>* ids = nullptr)
{
    if (img.size() < 12 || std::memcmp (img.data(), "RIFF", 4) != 0 || std::memcmp (img.data() + 8, "WAVE", 4) != 0) return "not RIFF/WAVE";
    if (u32at (img, 4) != img.size() - 8) return "RIFF size is not file size - 8";
    std::size_t pos = 12;
    while (pos < img.size())
    {
        if (pos % 2 != 0) return "a chunk header at an odd offset";
        if (img.size() - pos < 8) return "a truncated chunk header";
        const std::uint32_t size = u32at (img, pos + 4);
        const std::size_t padded = (std::size_t) size + (size & 1u);
        if (padded > img.size() - pos - 8) return "a chunk (with its pad) overruns the file";
        if ((size & 1u) && img[pos + 8 + size] != 0) return "a pad byte that is not zero";
        if (ids) ids->emplace_back ((const char*) img.data() + pos, 4);
        pos += 8 + padded;
    }
    return {};
}

// The chunk ids a READER walks to: readWavMemory's own loop — overrun rejects, a final chunk whose pad is
// missing ends the walk — with the pad step switchable. skipPad = false is THE PLANTED READER, advancing by
// the size field alone. Empty = the walk rejected the file.
std::vector<std::string> readerIds (const Bytes& img, bool skipPad)
{
    std::vector<std::string> ids;
    std::size_t pos = 12;
    while (pos + 8 <= img.size())
    {
        const std::uint32_t size = u32at (img, pos + 4);
        const std::size_t body = pos + 8;
        if (size > img.size() - body) return {};
        ids.emplace_back ((const char*) img.data() + pos, 4);
        const std::size_t advance = (std::size_t) size + (skipPad ? (size & 1u) : 0u);
        if (advance > img.size() - body) break;
        pos = body + advance;
    }
    return ids;
}

// THE PLANTED WRITER: the layout writeWavMemory emitted before the fix — the same bytes without the pad, the
// RIFF size 36 + data. An even image is its own old layout.
Bytes oldLayout (Bytes img)
{
    if (u32at (img, 40) & 1u) img.pop_back();
    fixRiff (img);
    return img;
}

// A LIST chunk appended at the end of the image, the way a tagger adds one, RIFF size moved to match.
Bytes withList (Bytes img)
{
    chunk (img, "LIST", { 'I', 'N', 'F', 'O' });
    fixRiff (img);
    return img;
}

// Samples that survive the trip exactly: PCM codes k / 2^(bits-1) over the whole range, float values that are
// floats. A fixed LCG, so every run writes the same images.
std::vector<std::vector<double>> programme (int bits, bool isFloat, int nch, std::size_t frames)
{
    std::vector<std::vector<double>> x ((std::size_t) nch, std::vector<double> (frames));
    unsigned long long s = 0x9E3779B97F4A7C15ull ^ (unsigned long long) (bits * 131 + nch * 7) ^ frames;
    const double full = (double) (1LL << (bits - 1));
    for (auto& c : x)
        for (auto& v : c)
        {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            const double u = (double) (s >> 11) * 0x1.0p-53;                              // [0, 1)
            v = isFloat ? (double) (float) (u * 2.0 - 1.0)
                        : std::floor (u * 2.0 * full) / full - 1.0;                     // k / 2^(bits-1), k in [−full, full)
        }
    return x;
}

bool sameAudio (const io::WavData& w, const std::vector<std::vector<double>>& x)
{
    if (! w.ok || w.ch.size() != x.size()) return false;
    for (std::size_t c = 0; c < x.size(); ++c)
        if (w.ch[c] != x[c]) return false;
    return true;
}

void note (std::string& list, const std::string& name) { list += (list.empty() ? "" : ", ") + name; }
} // namespace

int main()
{
    std::printf ("felitronics::io RIFF alignment tests\n");

    struct Format { int bits; bool isFloat; const char* name; };
    const Format formats[] = { { 16, false, "16-bit" }, { 24, false, "24-bit" }, { 32, true, "f32" } };
    const int channels[] = { 1, 2, 3 };
    const std::size_t frameCounts[] = { 1, 2, 3, 1000, 1001 };

    group ("the layout: word-aligned, RIFF counts the pad, data does not — every depth × channels × frames");
    {
        std::string odd, notEven, riff, dataSize, pad, walk, writeRead, readWrite, append;
        std::string plantWalkMissed, plantWalkWrong, plantAppendMissed, plantAppendWrong;
        int cases = 0, oddCases = 0;
        for (const Format& f : formats)
            for (int nch : channels)
                for (std::size_t n : frameCounts)
                {
                    const std::string name = std::string (f.name) + "×" + std::to_string (nch) + "ch×" + std::to_string (n);
                    const auto x = programme (f.bits, f.isFloat, nch, n);
                    const Bytes img = io::writeWavMemory (x, 48000.0, f.bits, f.isFloat);
                    ++cases;
                    if (img.size() < 44) { note (walk, name + " (no image)"); continue; }

                    const std::uint32_t block = (std::uint32_t) nch * (std::uint32_t) (f.bits / 8);
                    const std::uint32_t data = u32at (img, 40);
                    const bool isOdd = (data & 1u) != 0;
                    oddCases += isOdd;
                    if (isOdd != (f.bits == 24 && nch % 2 == 1 && n % 2 == 1)) note (odd, name);

                    if (img.size() % 2 != 0) note (notEven, name);
                    if (u32at (img, 4) != img.size() - 8) note (riff, name);
                    if (data != (std::uint32_t) n * block || (img[32] | (img[33] << 8)) != (int) block) note (dataSize, name);
                    if (img.size() != 44 + (std::size_t) data + (data & 1u) || (isOdd && img.back() != 0)) note (pad, name);
                    std::vector<std::string> ids;
                    if (! conforms (img, &ids).empty() || ids != std::vector<std::string> { "fmt ", "data" }) note (walk, name);

                    const auto w = io::readWavMemory (img.data(), img.size());
                    if (! sameAudio (w, x) || w.frames() != n) note (writeRead, name);
                    if (! w.ok || io::writeWavMemory (w.ch, w.sr, f.bits, f.isFloat) != img) note (readWrite, name);

                    std::vector<std::string> tagged;
                    const Bytes app = withList (img);
                    if (! sameAudio (io::readWavMemory (app.data(), app.size()), x) || ! conforms (app, &tagged).empty()
                        || tagged != std::vector<std::string> { "fmt ", "data", "LIST" })
                        note (append, name);

                    // THE PLANTED WRITER, through the same instruments: it must fail on exactly the odd cases.
                    const Bytes old = oldLayout (img);
                    const bool walkFails = ! conforms (old).empty();
                    if (isOdd && ! walkFails) note (plantWalkMissed, name);
                    if (! isOdd && walkFails) note (plantWalkWrong, name);
                    const Bytes oldApp = withList (old);
                    const bool appendFails = ! sameAudio (io::readWavMemory (oldApp.data(), oldApp.size()), x);
                    if (isOdd && ! appendFails) note (plantAppendMissed, name);
                    if (! isOdd && appendFails) note (plantAppendWrong, name);
                }

        ok (cases == 45 && oddCases == 6 && odd.empty(),
            "the data chunk is odd in exactly the 6 cases 24-bit × odd channels × odd frames makes (" + odd + ")");
        ok (notEven.empty(), "every file is an even number of bytes (" + notEven + ")");
        ok (riff.empty(), "RIFF size == file size − 8, the pad counted (" + riff + ")");
        ok (dataSize.empty(), "data size == frames · blockAlign, the pad NOT counted (" + dataSize + ")");
        ok (pad.empty(), "one zero pad byte after odd data, none after even (" + pad + ")");
        ok (walk.empty(), "a strict RIFF walk covers each file exactly: fmt , data (" + walk + ")");
        ok (writeRead.empty(), "write∘read gives every sample back exactly (" + writeRead + ")");
        ok (readWrite.empty(), "read∘write gives every image back byte for byte (" + readWrite + ")");
        ok (append.empty(), "a LIST appended after the image: same audio, and the walk finds fmt , data, LIST (" + append + ")");
        ok (plantWalkMissed.empty() && plantWalkWrong.empty(),
            "planted: the old writer's layout fails the walk on exactly the odd cases (missed: " + plantWalkMissed
                + "; wrongly failed: " + plantWalkWrong + ")");
        ok (plantAppendMissed.empty() && plantAppendWrong.empty(),
            "planted: …and with a LIST appended our own reader refuses exactly the odd cases (missed: " + plantAppendMissed
                + "; wrongly failed: " + plantAppendWrong + ")");

        // The two images measured on v0.53.0, so the plant is the writer that was, not a guess at it.
        const Bytes one  = io::writeWavMemory ({ std::vector<double> (1, 0.25) }, 48000.0, 24, false);
        const Bytes many = io::writeWavMemory ({ std::vector<double> (1001, 0.25) }, 48000.0, 24, false);
        ok (one.size() == 48 && u32at (one, 4) == 40 && u32at (one, 40) == 3 && one.back() == 0,
            "24-bit mono, 1 frame: 48 bytes, RIFF 40, data 3, a zero pad");
        ok (many.size() == 3048 && u32at (many, 4) == 3040 && u32at (many, 40) == 3003,
            "24-bit mono, 1001 frames: 3048 bytes, RIFF 3040, data 3003");
        const Bytes oldOne = oldLayout (one), oldMany = oldLayout (many);
        ok (oldOne.size() == 47 && u32at (oldOne, 4) == 39 && oldMany.size() == 3047 && u32at (oldMany, 4) == 3039,
            "planted: the old layout is what v0.53.0 wrote — 47 bytes RIFF 39, 3047 bytes RIFF 3039");
    }

    group ("the reader steps over a pad wherever it sits");
    {
        // 24-bit mono, three codes: 0.25, −0.5, one LSB. 9 bytes of data — odd.
        const Bytes pcm24 = { 0x00, 0x00, 0x20,  0x00, 0x00, 0xC0,  0x01, 0x00, 0x00 };
        const std::vector<double> want24 = { 0.25, -0.5, 1.0 / 8388608.0 };
        Bytes fmt24;
        u16 (fmt24, 1); u16 (fmt24, 1); u32 (fmt24, 48000); u32 (fmt24, 48000 * 3); u16 (fmt24, 3); u16 (fmt24, 24);
        // 16-bit stereo, two frames — even data, so only the chunk BEFORE fmt is odd.
        const Bytes pcm16 = { 0x00, 0x40,  0x00, 0xE0,  0xFF, 0x7F,  0x00, 0x80 };
        Bytes fmt16;
        u16 (fmt16, 1); u16 (fmt16, 2); u32 (fmt16, 44100); u32 (fmt16, 44100 * 4); u16 (fmt16, 4); u16 (fmt16, 16);
        const Bytes list5 = { 'I', 'N', 'F', 'O', 'x' }, id3 = { 'I', 'D', '3', 4, 0, 0, 0 }, junk3 = { 1, 2, 3 };

        struct Img { const char* name; Bytes bytes; std::vector<std::string> ids; bool is24; bool oddBeforeAnother; bool conforming; };
        std::vector<Img> imgs;
        auto start = [] { Bytes o; tag (o, "RIFF"); u32 (o, 0); tag (o, "WAVE"); return o; };
        {
            Bytes o = start(); chunk (o, "LIST", list5); chunk (o, "fmt ", fmt24); chunk (o, "data", pcm24); fixRiff (o);
            imgs.push_back ({ "an odd LIST before fmt, odd data + pad at the end", o, { "LIST", "fmt ", "data" }, true, true, true });
        }
        {
            Bytes o = start(); chunk (o, "fmt ", fmt24); chunk (o, "data", pcm24); chunk (o, "junk", junk3); chunk (o, "LIST", { 'I', 'N', 'F', 'O' }); fixRiff (o);
            imgs.push_back ({ "odd data, then an odd chunk, then a LIST", o, { "fmt ", "data", "junk", "LIST" }, true, true, true });
        }
        {
            Bytes o = start(); chunk (o, "LIST", list5); chunk (o, "id3 ", id3); chunk (o, "fmt ", fmt24); chunk (o, "data", pcm24); chunk (o, "LIST", list5); fixRiff (o);
            imgs.push_back ({ "odd, odd, fmt, odd data, odd LIST — a pad before and after everything", o, { "LIST", "id3 ", "fmt ", "data", "LIST" }, true, true, true });
        }
        {
            Bytes o = start(); chunk (o, "LIST", list5); chunk (o, "fmt ", fmt16); chunk (o, "data", pcm16); fixRiff (o);
            imgs.push_back ({ "16-bit stereo behind an odd LIST", o, { "LIST", "fmt ", "data" }, false, true, true });
        }
        {
            Bytes o = start(); chunk (o, "fmt ", fmt24); chunk (o, "data", pcm24, false); fixRiff (o);
            imgs.push_back ({ "odd data at EOF with NO pad — what the writer made before", o, { "fmt ", "data" }, true, false, false });
        }

        int discriminating = 0;
        for (const Img& m : imgs)
        {
            const std::string at = std::string (" — ") + m.name;
            ok (conforms (m.bytes).empty() == m.conforming,
                std::string (m.conforming ? "fixture: a conforming RIFF image" : "fixture: NOT conforming, on purpose") + at);
            ok (readerIds (m.bytes, true) == m.ids, "fixture: the pad-skipping walk finds exactly the planted chunks" + at);
            const auto w = io::readWavMemory (m.bytes.data(), m.bytes.size());
            if (m.is24)
                ok (w.ok && w.bits == 24 && w.ch.size() == 1 && w.ch[0] == want24, "readWavMemory reads the three codes exactly" + at);
            else
                ok (w.ok && w.bits == 16 && w.ch.size() == 2 && w.frames() == 2
                        && w.ch[0][0] == 0.5 && w.ch[1][0] == -0.25 && w.ch[0][1] == 32767.0 / 32768.0 && w.ch[1][1] == -1.0,
                    "readWavMemory reads both frames exactly" + at);
            // THE PLANTED READER: without the pad step it must lose every image with an odd chunk before another;
            // where nothing follows the odd chunk the pad step is moot and both walks agree.
            const bool planted = readerIds (m.bytes, false) != m.ids;
            ok (planted == m.oddBeforeAnother,
                std::string (m.oddBeforeAnother ? "planted: a reader that does not skip the pad misparses it"
                                                : "planted: nothing follows the odd chunk, so both walks agree") + at);
            discriminating += planted;
        }
        ok (discriminating == 4, "planted: 4 of the 5 images tell a pad-skipping reader from one that is not");
    }

    return felitronics::test::report();
}
