#include "ducklips/planner.hpp"
#include "ducklips/synthesizer.hpp"
#include <cstdio>
#include <string>

namespace {
int failures{};
void check(bool okay, const char* label) {
    std::printf("%s %s\n", okay ? "PASS" : "FAIL", label);
    if (!okay) { ++failures; }
}
bool same_phones(const ducklips::UtterancePlan& a, const ducklips::UtterancePlan& b) {
    if (a.phones.size() != b.phones.size()) { return false; }
    for (std::size_t i = 0; i < a.phones.size(); ++i) {
        if (a.phones[i].phone != b.phones[i].phone || a.phones[i].stress != b.phones[i].stress
            || a.phones[i].start != b.phones[i].start || a.phones[i].end != b.phones[i].end) { return false; }
    }
    return true;
}
}

int main() {
    ducklips::Synthesizer engine{};
    const ducklips::VoiceProfile body = ducklips::baritone_voice();
    const ducklips::SynthesisOptions options{};
    engine.prepare(body, options);
    ducklips::UtterancePlanner planner{};
    planner.prepare(body, engine.vowels(), options.rate);
    check(planner.load_dictionary(std::string(DUCKLIPS_SOURCE_DIR) + "/samples/nato.dict") == 30,
          "NATO pronunciation overrides installed");
    const ducklips::UtterancePlan xray = planner.plan_text("X-ray.");
    check(xray.guessed_words.empty() && xray.transcription == "EH1 K S R EY0",
          "official X-ray spelling retains its initial vowel and dictionary stress");
    check(same_phones(xray, planner.plan_text("XRAY.")), "hyphenated alias preserves phone timing");
    check(same_phones(planner.plan_text("x-RaY."), xray), "compound lookup is case insensitive");
    check(same_phones(planner.plan_text("unknown-compound."), planner.plan_text("unknown compound.")),
          "unknown hyphenated compounds remain separate words");
    check(same_phones(planner.plan_text("X--ray."), planner.plan_text("X ray.")),
          "double hyphens remain word boundaries");
    check(same_phones(planner.plan_text("-X-ray-."), xray), "leading and trailing hyphens do not alter a known compound");
    const ducklips::UtterancePlan nato = planner.plan_text(
        "Alfa. Bravo. Charlie. Delta. Echo. Foxtrot. Golf. Hotel. India. Juliett. Kilo. Lima. Mike. "
        "November. Oscar. Papa. Quebec. Romeo. Sierra. Tango. Uniform. Victor. Whiskey. X-ray. Yankee. Zulu.");
    check(nato.guessed_words.empty(), "all 26 NATO spellings use explicit dictionary entries");
    return failures ? 1 : 0;
}
