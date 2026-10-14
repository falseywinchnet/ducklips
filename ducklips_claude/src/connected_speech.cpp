#include "ducklips/connected_speech.hpp"
#include <array>
#include <sstream>

namespace ducklips {
namespace {

struct WeakForm final {
    const char* word{};
    const char* phones{};
};

// Unstressed forms used when the word is not the last of its phrase.
constexpr std::array<WeakForm, 17> weak_forms{{
    {"A", "AX0"}, {"AN", "AX0 N"}, {"AND", "AX0 N D"}, {"AS", "AX0 Z"}, {"AT", "AX0 T"},
    {"BUT", "B AX0 T"}, {"CAN", "K AX0 N"}, {"FOR", "F ER0"}, {"FROM", "F R AX0 M"},
    {"OF", "AX0 V"}, {"OR", "ER0"}, {"THAN", "DH AX0 N"}, {"THE", "DH AX0"},
    {"TO", "T AX0"}, {"WAS", "W AX0 Z"}, {"YOU", "Y UW0"}, {"YOUR", "Y ER0"},
}};

bool weak_form(const std::string& word, std::vector<PhoneToken>& phones) {
    for (std::size_t index = 0; index < weak_forms.size(); ++index) {
        if (word == weak_forms[index].word) {
            phones.clear();
            std::istringstream stream(weak_forms[index].phones);
            std::string symbol{};
            while (stream >> symbol) {
                PhoneToken token{};
                if (parse_phone(symbol, token.phone, token.stress)) {
                    phones.push_back(token);
                }
            }
            return true;
        }
    }
    return false;
}

bool is_consonant(Phone phone) noexcept {
    bool consonant = !is_vowel(phone);
    return consonant;
}

bool vowel_like_before_flap(Phone phone) noexcept {
    bool result = is_vowel(phone) || phone == Phone::R;
    return result;
}

void note(std::vector<std::string>& events, const char* kind, const std::string& word) {
    std::string event = std::string(kind) + ": " + word;
    events.push_back(event);
}

}

std::vector<PhoneToken> join_words(const std::vector<SpokenWord>& words, std::vector<std::string>& events) {
    std::vector<std::vector<PhoneToken>> forms{};
    for (std::size_t index = 0; index < words.size(); ++index) {
        std::vector<PhoneToken> phones = words[index].phones;
        bool has_following = index + 1 < words.size();
        if (has_following && weak_form(words[index].word, phones)) {
            note(events, "weak form", words[index].word);
        }
        forms.push_back(phones);
    }
    // Word-boundary events, left to right.
    for (std::size_t index = 0; index + 1 < forms.size(); ++index) {
        std::vector<PhoneToken>& left = forms[index];
        std::vector<PhoneToken>& right = forms[index + 1];
        if (left.empty() || right.empty()) {
            continue;
        }
        Phone next = right[0].phone;
        // Cluster loss: "and the", "last night".
        if (left.size() >= 2 && is_consonant(next)) {
            Phone final_phone = left[left.size() - 1].phone;
            Phone before = left[left.size() - 2].phone;
            bool droppable = (final_phone == Phone::D && before == Phone::N)
                || (final_phone == Phone::T && (before == Phone::S || before == Phone::F || before == Phone::K));
            if (droppable) {
                left.pop_back();
                note(events, "cluster loss", words[index].word);
            }
        }
        PhoneToken& last = left[left.size() - 1];
        if (last.phone == Phone::N && (next == Phone::P || next == Phone::B || next == Phone::M)) {
            last.phone = Phone::M;
            note(events, "assimilation", words[index].word);
        } else if (last.phone == Phone::N && (next == Phone::K || next == Phone::G)) {
            last.phone = Phone::NG;
            note(events, "assimilation", words[index].word);
        } else if ((last.phone == Phone::T || last.phone == Phone::D) && next == Phone::Y) {
            last.phone = last.phone == Phone::T ? Phone::CH : Phone::JH;
            right.erase(right.begin());
            note(events, "palatalization", words[index].word);
        } else if (is_consonant(last.phone) && last.phone == next) {
            last.duration_scale *= 1.6;
            right.erase(right.begin());
            note(events, "geminate", words[index].word);
        }
    }
    std::vector<PhoneToken> stream{};
    std::vector<std::size_t> word_of{};
    for (std::size_t index = 0; index < forms.size(); ++index) {
        for (std::size_t phone = 0; phone < forms[index].size(); ++phone) {
            stream.push_back(forms[index][phone]);
            word_of.push_back(index);
        }
    }
    // Flapping across the whole phrase, inside words and across boundaries.
    for (std::size_t index = 1; index + 1 < stream.size(); ++index) {
        PhoneToken& token = stream[index];
        bool alveolar_stop = token.phone == Phone::T || token.phone == Phone::D;
        if (!alveolar_stop) {
            continue;
        }
        const PhoneToken& previous = stream[index - 1];
        const PhoneToken& following = stream[index + 1];
        if (vowel_like_before_flap(previous.phone) && is_vowel(following.phone) && following.stress == 0) {
            token.flap = true;
            note(events, "flap", words[word_of[index]].word);
        }
    }
    return stream;
}

}
