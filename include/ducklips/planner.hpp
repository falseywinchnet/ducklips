#pragma once
#include "ducklips/cues.hpp"
#include "ducklips/phonetics.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/vowel_library.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace ducklips {

struct PlannedPhone final {
    Phone phone{Phone::AX};
    int stress{-1};
    double start{}; // s
    double end{};   // s
};

struct UtterancePlan final {
    std::vector<PlannedPhone> phones{};
    std::vector<Cue> cues{};            // sorted by time
    double duration{};                  // s, until the finish cue
    std::string transcription{};        // ARPAbet with stress and phrase bars
    std::vector<std::string> guessed_words{}; // pronounced by letter-to-sound rules
    std::vector<std::string> transitions{};   // connected-speech events applied
};

// Turns text into a behavioural cue stream for one voice.
//
// Timing follows simple segmental duration rules (stress, tenseness, voicing
// of the following consonant, phrase-final lengthening, cluster shortening)
// with a small deterministic variation. Gestures are scheduled with
// anticipatory overlap, so coarticulation comes from the articulator dynamics
// rather than from transition tables. Breath groups get an inhale, a pressure
// rise and a trailing pressure fall. Intonation is declination plus accents on
// stressed syllables, with a final fall or, for questions, a rise.
class UtterancePlanner final {
public:
    void prepare(const VoiceProfile& voice, const VowelLibrary& vowels, SynthesisRate rate);
    // Adds a CMU-format dictionary file; returns entries added.
    std::size_t load_dictionary(const std::string& path);
    [[nodiscard]] UtterancePlan plan_text(const std::string& text) const;
    // Space-separated ARPAbet; "|" separates phrases; "AA1:0.8" forces a duration.
    [[nodiscard]] UtterancePlan plan_phonemes(const std::string& phonemes) const;

private:
    struct Phrase final {
        std::vector<PhoneToken> tokens{};
        bool question{false};
        double pause_after{0.25}; // s
    };
    [[nodiscard]] UtterancePlan plan_phrases(const std::vector<Phrase>& phrases) const;
    [[nodiscard]] double base_duration(const std::vector<PhoneToken>& tokens, std::size_t index) const;

    VoiceProfile voice_{};
    VowelLibrary vowels_{};
    Lexicon lexicon_{};
    SynthesisRate rate_{};
};

}
