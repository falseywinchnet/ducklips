#include "ducklips/body_config.hpp"
#include "ducklips/nudges.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <sstream>
#include <vector>

namespace {
std::atomic<std::size_t> allocations{0};
std::size_t failures = 0;
void check(bool okay, const char* label) {
    std::printf("%s %s\n", okay ? "PASS" : "FAIL", label);
    if (!okay) { ++failures; }
}
}
void* operator new(std::size_t count) {
    ++allocations;
    void* result = std::malloc(std::max<std::size_t>(count, 1));
    if (!result) { throw std::bad_alloc(); }
    return result;
}
void* operator new[](std::size_t count) { return ::operator new(count); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }

int main() {
    ducklips::VoiceProfile body = ducklips::baritone_voice();
    std::string error{};
    std::istringstream good("preset=alto\nname=example\ntract_length=.15\nfolds.lower_mass=.00007\nspeaking_rate=1.2\n");
    check(ducklips::read_body_config(good, body, error) && body.name == "example" && body.tract_length == .15
          && body.folds.lower_mass == .00007 && body.base_pitch == ducklips::alto_voice().base_pitch,
          "body preset and anatomy overrides are applied");
    const char* bad_cases[] = {"tract_length=.18\nunknown=1", "base_pitch=nan", "lung_pressure=inf",
        "tract_length=-1", "base_pitch=100garbage", "base_pitch=110\nbase_pitch=120", "name=test\npreset=alto",
        "audible_breath=yes", "adducted_area=.00008\nabducted_area=.00001"};
    for (const char* text : bad_cases) {
        std::istringstream bad(text);
        check(!ducklips::read_body_config(bad, body, error) && body.name == "example" && body.tract_length == .15,
              "invalid body input fails transactionally");
    }
    ducklips::Lexicon lexicon{};
    std::size_t installed = lexicon.load_cmu_file(std::string(DUCKLIPS_SOURCE_DIR) + "/samples/nato.dict");
    ducklips::WordPronunciation papa = lexicon.pronounce("papa");
    ducklips::WordPronunciation quebec = lexicon.pronounce("quebec");
    check(installed == 29 && papa.phones[1].phone == ducklips::Phone::AX && papa.phones[3].stress == 1
          && quebec.phones[1].phone == ducklips::Phone::EH, "explicit NATO dictionary overrides built-in pronunciations");
    body = ducklips::baritone_voice();
    ducklips::NudgeStream whole{};
    ducklips::NudgeStream chunks{};
    whole.prepare(body);
    chunks.prepare(body);
    const char* lines[] = {"0 lungs 600 40", "0 glottis .000003 170", "0 vowel AA 75",
        ".1 diphthong AY .4 75", ".2 pitch 125 90", ".65 constrict S 1 120", ".65 glottis .000025 170",
        ".8 release tip 150", ".8 lungs 0 40"};
    for (const char* line : lines) {
        check(whole.submit(line, error) && chunks.submit(line, error), "valid overlapping nudge accepted");
    }
    std::size_t pending = whole.pending();
    const char* invalid[] = {"-1 pitch 100 40", "0 vowel S 40", "0 diphthong AA .4 60", "0 pitch nan 50",
        "0 release forehead 40", "0 lungs 600 0", "0 lungs 600 40 garbage", "0 constrict S inf 120"};
    for (const char* line : invalid) {
        check(!whole.submit(line, error) && whole.pending() == pending, "invalid nudge does not change schedule");
    }
    std::vector<float> a(48000 * 2);
    std::vector<float> b(a.size());
    std::size_t before = allocations.load();
    whole.render(a);
    check(allocations.load() == before, "render allocates no heap storage including cue delivery");
    for (std::size_t frame = 0; frame < 48000;) {
        std::size_t run = std::min<std::size_t>(17, 48000 - frame);
        chunks.render(std::span<float>(b.data() + frame * 2, run * 2));
        frame += run;
    }
    check(a == b, "audio is bit-identical for full buffer and 17-frame host chunks");
    double energy = 0;
    bool finite = true;
    for (float value : a) { finite = finite && std::isfinite(value); energy += value * value; }
    check(finite && energy > .001 && !whole.engine().statistics().nonfinite_samples
          && !whole.engine().statistics().clipped_samples, "stream produces finite non-silent unclipped speech audio");
    check(!whole.submit(".9 pitch 100 50", error), "late nudge is rejected");
    check(whole.submit("1 pitch 100 50", error), "next uncommitted block accepts nudge");
    whole.render(std::span<float>(a.data(), 2));
    check(!whole.submit("1.00001 pitch 100 50", error), "partially rendered control block rejects retroactive nudge");
    ducklips::NudgeStream bounded{};
    bounded.prepare(body, 2.0, 4);
    check(bounded.submit("0 diphthong AY .4 75", error) && bounded.submit("0 diphthong AW .4 75", error),
          "bounded stream accepts exactly its cue capacity");
    check(!bounded.submit("1 pitch 110 90", error) && bounded.pending() == 4, "bounded stream rejects overflow");
    std::vector<float> half_second(24000 * 2);
    bounded.render(half_second);
    check(bounded.pending() == 0, "speed-scaled diphthong targets are delivered");
    std::printf("%zu failures\n", failures);
    return failures ? 1 : 0;
}
