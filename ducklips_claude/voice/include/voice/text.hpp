#pragma once

// From English text to phrased, stressed phones.
//
//   text -> words and punctuation -> phrases (with breath decisions)
//        -> phones with lexical stress -> syllables
//
// The phone alphabet is ARPAbet. Vowels carry stress 0 (none), 1 (primary)
// or 2 (secondary).

#include <cstdint>
#include <string>
#include <vector>

namespace voice {

enum class Boundary : std::uint8_t {
  none,       // inside a phrase
  minor,      // short pause, pitch continues (comma, conjunction split)
  major,      // sentence end
  question    // sentence end with a rising question
};

struct Phone {
  std::string name;        // ARPAbet without stress digit
  std::int32_t stress = -1;  // vowels: 0, 1, 2; consonants: -1
};

struct Word {
  std::string text;
  std::vector<Phone> phones;
  bool function = false;   // grammatical word, normally unaccented
  bool accent = false;     // carries a pitch accent
  bool nuclear = false;    // last accent in its phrase
};

struct Phrase {
  std::vector<Word> words;
  Boundary end = Boundary::minor;
  bool breath_before = false;
  double breath_seconds = 0.0;
};

std::vector<Phrase> parse_text(const std::string& text);

// Pronunciation of one lower-case word: the lexicon first, then
// letter-to-sound rules with a stress heuristic.
std::vector<Phone> pronounce(const std::string& word);
bool is_function_word(const std::string& word);
// The reduced pronunciation of a function word, if it has one.
bool weak_form(const std::string& word, std::vector<Phone>& phones);
bool is_vowel_phone(const std::string& name);

std::string describe_phrases(const std::vector<Phrase>& phrases);

}  // namespace voice
