// Text to phrases: tokens, numbers, phrase boundaries, breaths and accents.

#include "voice/text.hpp"

#include <cctype>
#include <cstdio>

namespace voice {

namespace {

constexpr const char* kOnes[] = {"zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
                                 "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen",
                                 "seventeen", "eighteen", "nineteen"};
constexpr const char* kTens[] = {"", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"};

// Words that may begin a new phrase when the current one is long.
constexpr const char* kPhraseOpeners[] = {"and", "but", "or", "because", "which", "that", "when", "while",
                                          "although", "though", "if", "so", "where", "who", "until", "since"};

void append_number_words(long value, std::vector<std::string>& words) {
  if (value < 0) {
    words.push_back("minus");
    value = -value;
  }
  if (value >= 1000000) {
    append_number_words(value / 1000000, words);
    words.push_back("million");
    value %= 1000000;
    if (value == 0) {
      return;
    }
  }
  if (value >= 1000) {
    append_number_words(value / 1000, words);
    words.push_back("thousand");
    value %= 1000;
    if (value == 0) {
      return;
    }
  }
  if (value >= 100) {
    words.push_back(kOnes[value / 100]);
    words.push_back("hundred");
    value %= 100;
    if (value == 0) {
      return;
    }
  }
  if (value >= 20) {
    words.push_back(kTens[value / 10]);
    if (value % 10 != 0) {
      words.push_back(kOnes[value % 10]);
    }
    return;
  }
  words.push_back(kOnes[value]);
}

struct Token {
  std::string text;
  char punctuation = 0;  // 0 for a word
};

std::vector<Token> tokenize(const std::string& text) {
  std::vector<Token> tokens;
  std::size_t index = 0;
  while (index < text.size()) {
    unsigned char letter = static_cast<unsigned char>(text[index]);
    if (std::isalpha(letter) != 0 || letter == '\'') {
      std::string word;
      while (index < text.size()) {
        unsigned char next = static_cast<unsigned char>(text[index]);
        if (std::isalpha(next) == 0 && next != '\'') {
          break;
        }
        word.push_back(static_cast<char>(std::tolower(next)));
        ++index;
      }
      while (!word.empty() && word.front() == '\'') {
        word.erase(word.begin());
      }
      while (!word.empty() && word.back() == '\'') {
        word.pop_back();
      }
      if (!word.empty()) {
        Token token;
        token.text = word;
        tokens.push_back(token);
      }
    } else if (std::isdigit(letter) != 0) {
      long value = 0;
      while (index < text.size() && std::isdigit(static_cast<unsigned char>(text[index])) != 0) {
        value = value * 10 + (text[index] - '0');
        ++index;
      }
      std::vector<std::string> words;
      append_number_words(value, words);
      for (const std::string& word : words) {
        Token token;
        token.text = word;
        tokens.push_back(token);
      }
    } else if (letter == ',' || letter == ';' || letter == ':' || letter == '.' || letter == '!' || letter == '?' ||
               letter == '-') {
      if (letter == '-' && (index + 1 >= text.size() || text[index + 1] != '-')) {
        ++index;  // a hyphen inside a compound joins words
        continue;
      }
      Token token;
      token.punctuation = static_cast<char>(letter);
      tokens.push_back(token);
      ++index;
      while (index < text.size() && text[index] == '-') {
        ++index;
      }
    } else {
      ++index;
    }
  }
  return tokens;
}

bool is_phrase_opener(const std::string& word) {
  for (const char* opener : kPhraseOpeners) {
    if (word == opener) {
      return true;
    }
  }
  return false;
}

int count_syllables(const Word& word) {
  int count = 0;
  for (const Phone& phone : word.phones) {
    if (is_vowel_phone(phone.name)) {
      ++count;
    }
  }
  return count;
}

// Unstressed AH is a schwa: the bdl labels have ax for 49 of our 71.
void reduce_unstressed_ah(std::vector<Phone>& phones) {
  for (Phone& phone : phones) {
    if (phone.name == "AH" && phone.stress == 0) {
      phone.name = "AX";
    }
  }
}

Word make_word(const std::string& text) {
  Word word;
  word.text = text;
  word.phones = pronounce(text);
  word.function = is_function_word(text);
  if (word.function) {
    // Function words lose their stress in running speech, and most reduce
    // their vowel.
    std::vector<Phone> weak;
    if (weak_form(text, weak)) {
      word.phones = weak;
    }
    for (Phone& phone : word.phones) {
      if (phone.stress > 0) {
        phone.stress = 0;
      }
    }
  }
  reduce_unstressed_ah(word.phones);
  return word;
}

// Pitch accents: content words carry them, except that two accents are not
// placed on adjacent short content words (the first is dropped); the last
// accent of a phrase is nuclear.
void place_accents(Phrase& phrase) {
  bool any = false;
  for (std::size_t index = 0; index < phrase.words.size(); ++index) {
    Word& word = phrase.words[index];
    word.accent = !word.function;
    any = any || word.accent;
  }
  for (std::size_t index = 0; index + 1 < phrase.words.size(); ++index) {
    Word& word = phrase.words[index];
    const Word& next = phrase.words[index + 1];
    if (word.accent && next.accent && count_syllables(word) == 1 && index + 2 < phrase.words.size()) {
      word.accent = false;
    }
  }
  if (!any && !phrase.words.empty()) {
    phrase.words.back().accent = true;
  }
  for (std::size_t index = phrase.words.size(); index > 0; --index) {
    if (phrase.words[index - 1].accent) {
      phrase.words[index - 1].nuclear = true;
      break;
    }
  }
}

}  // namespace

namespace {

// Weak forms depend on their neighbours: "the" and "to" keep a tense
// vowel before a vowel ("the end", "to all": DH IY, T UW), and a function
// word ending its phrase keeps its full form ("what is it for").
void settle_weak_forms(Phrase& phrase) {
  for (std::size_t index = 0; index < phrase.words.size(); ++index) {
    Word& word = phrase.words[index];
    if (!word.function) {
      continue;
    }
    bool last = index + 1 == phrase.words.size();
    if (last) {
      word.phones = pronounce(word.text);
      for (Phone& phone : word.phones) {
        if (phone.stress > 0) {
          phone.stress = 0;
        }
      }
      reduce_unstressed_ah(word.phones);
      continue;
    }
    const Word& next = phrase.words[index + 1];
    bool next_vowel = !next.phones.empty() && is_vowel_phone(next.phones.front().name);
    if (next_vowel && !word.phones.empty() && word.phones.back().name == "AX") {
      if (word.text == "the") {
        word.phones.back().name = "IY";
      } else if (word.text == "to") {
        word.phones.back().name = "UW";
      }
    }
  }
}

}  // namespace

std::vector<Phrase> parse_text(const std::string& text) {
  std::vector<Token> tokens = tokenize(text);
  std::vector<Phrase> phrases;
  Phrase current;
  int syllables_since_breath = 0;
  bool first_phrase = true;

  for (std::size_t index = 0; index < tokens.size(); ++index) {
    const Token& token = tokens[index];
    if (token.punctuation == 0) {
      // A long phrase may split before a conjunction or relative word.
      int current_syllables = 0;
      for (const Word& word : current.words) {
        current_syllables += count_syllables(word);
      }
      if (is_phrase_opener(token.text) && current.words.size() >= 4 && current_syllables >= 7) {
        current.end = Boundary::minor;
        phrases.push_back(current);
        current = Phrase();
      }
      current.words.push_back(make_word(token.text));
      continue;
    }
    if (current.words.empty()) {
      continue;
    }
    if (token.punctuation == '?') {
      current.end = Boundary::question;
    } else if (token.punctuation == '.' || token.punctuation == '!') {
      current.end = Boundary::major;
    } else {
      current.end = Boundary::minor;
    }
    phrases.push_back(current);
    current = Phrase();
  }
  if (!current.words.empty()) {
    current.end = Boundary::major;
    phrases.push_back(current);
  }
  for (Phrase& phrase : phrases) {
    settle_weak_forms(phrase);
  }

  // Breaths: always before speaking; then where the air would run low
  // (about 12 syllables at conversational pressure), preferring sentence
  // ends and taking a quick catch-breath at a minor boundary when needed.
  for (std::size_t index = 0; index < phrases.size(); ++index) {
    Phrase& phrase = phrases[index];
    int syllables = 0;
    for (const Word& word : phrase.words) {
      syllables += count_syllables(word);
    }
    if (first_phrase) {
      phrase.breath_before = true;
      phrase.breath_seconds = 0.45;
      first_phrase = false;
      syllables_since_breath = 0;
    } else {
      Boundary previous = phrases[index - 1].end;
      bool sentence = previous == Boundary::major || previous == Boundary::question;
      if (sentence && syllables_since_breath + syllables > 8) {
        phrase.breath_before = true;
        phrase.breath_seconds = 0.40;
      } else if (syllables_since_breath + syllables > 16) {
        phrase.breath_before = true;
        phrase.breath_seconds = sentence ? 0.40 : 0.28;
      }
      if (phrase.breath_before) {
        syllables_since_breath = 0;
      }
    }
    syllables_since_breath += syllables;
    place_accents(phrase);
  }
  return phrases;
}

std::string describe_phrases(const std::vector<Phrase>& phrases) {
  std::string text;
  for (const Phrase& phrase : phrases) {
    if (phrase.breath_before) {
      char line[64];
      std::snprintf(line, sizeof(line), "[breath %.2f] ", phrase.breath_seconds);
      text += line;
    }
    text += "{";
    for (const Word& word : phrase.words) {
      text += " ";
      if (word.nuclear) {
        text += "**";
      } else if (word.accent) {
        text += "*";
      }
      text += word.text;
      text += "(";
      for (std::size_t index = 0; index < word.phones.size(); ++index) {
        if (index > 0) {
          text += " ";
        }
        text += word.phones[index].name;
        if (word.phones[index].stress >= 0) {
          text += static_cast<char>('0' + word.phones[index].stress);
        }
      }
      text += ")";
    }
    text += " }";
    if (phrase.end == Boundary::major) {
      text += " .\n";
    } else if (phrase.end == Boundary::question) {
      text += " ?\n";
    } else {
      text += " ,\n";
    }
  }
  return text;
}

}  // namespace voice
