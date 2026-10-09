#pragma once
#include "ducklips/phonetics.hpp"
#include <string>
#include <vector>

namespace ducklips {

struct SpokenWord final {
    std::string word{};             // upper case
    std::vector<PhoneToken> phones{};
};

// Joins the words of one phrase into a phone stream, applying the transition
// events of connected English speech:
//   weak forms     function words before more speech lose their full vowels
//   assimilation   /n/ takes the place of a following labial or velar stop
//   palatalization /t d/ + /j/ become /tS dZ/ ("did you")
//   cluster loss   final /nd st ft kt/ lose the stop before a consonant
//   geminates      a consonant ending one word and starting the next is held once
//   flapping       /t d/ between a vowel and an unstressed vowel become a tap
// Each applied event is described in events ("flap: WATER").
[[nodiscard]] std::vector<PhoneToken> join_words(const std::vector<SpokenWord>& words,
                                                 std::vector<std::string>& events);

}
