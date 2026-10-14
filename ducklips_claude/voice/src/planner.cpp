// Phones to postures: durations, intonation, breath and the gestures of
// every part of the body.

#include "voice/planner.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>

namespace voice {

namespace {

enum class Kind : std::uint8_t { vowel, diphthong, stop, fricative, affricate, nasal, liquid, glide, aspirate };
enum class Organ : std::uint8_t { none, lips, tip, dorsum, glottis };

struct Posture {
  double place = 10.5;
  double area = 1.85;
  double lips = 1.9;
  double protrusion = 0.0;
};

struct PhoneSpec {
  const char* name;
  Kind kind;
  bool voiced;
  double inherent_ms;
  double minimum_ms;
  Posture vowel;      // vowels: target; diphthongs: start
  Posture glide_to;   // diphthongs: end
  Organ organ;        // consonants: the part that makes the constriction
  double constriction_place;
  double constriction_area;  // 0 closes
  double floor;              // sublingual space under a raised blade, cm^3
  double pocket = 0.0;       // laterals: opening of the pocket behind the contact, cm^2
  // Vowels: a second shaping point, the front of the tongue raised toward
  // the palate (0 area: none). The body alone tops out near F2 1700 Hz, so
  // /ae/ and /eh/ came out as one vowel (560/1700 vs 580/1760 Hz); bdl has
  // 581/1885 vs 594/1785 and the owner heard "Alfa" as "elfa".
  double dorsum_place = 0.0;
  double dorsum_area = 0.0;
};

// Vowel postures: tools/vowel_fit bounded, against Hillenbrand 1995 men (UW
// toward the fronted American /u/, AX a connected-speech schwa), measured
// with the glottis slightly apart as in voicing, lip openings held to
// anthropometric ranges. Wider openings radiate F3 and F4 away: at 9 cm^2
// their bandwidths reach 250 Hz, and the 2-4 kHz band of running speech fell
// 8 dB below the human's. Front F2 tops out near 1700 Hz for EH/AE in this
// geometry. (REJECTED: a refit after the epilarynx was lengthened to
// 3.5 cm reached front F2 targets in the impulse fit but scored worse in
// running speech: phone test 30.0 vs 32.0%.)
constexpr Posture kIY = {13.00, 0.250, 3.500, 0.0};
constexpr Posture kIH = {11.25, 0.697, 4.500, 0.0};
constexpr Posture kEY = {9.56, 0.260, 5.000, 0.0};
constexpr Posture kEH = {9.50, 1.620, 6.000, 0.0};
constexpr Posture kAE = {9.50, 1.394, 7.134, 0.0};
constexpr Posture kAA = {6.25, 1.014, 7.000, 0.0};
constexpr Posture kAO = {6.00, 0.516, 1.832, 0.3};
constexpr Posture kOW = {5.00, 0.526, 0.800, 0.6};
constexpr Posture kUH = {8.62, 0.479, 1.000, 0.4};
constexpr Posture kUW = {13.88, 2.500, 0.600, 0.6};
constexpr Posture kAH = {6.50, 0.516, 3.664, 0.0};
constexpr Posture kAX = {11.00, 2.500, 2.000, 0.0};
constexpr Posture kER = {3.50, 1.200, 2.00, 0.3};  // plus a bunched dorsum, see below
constexpr Posture kNone = {10.75, 2.500, 1.62, 0.0};
// Tongue body behind an alveolar closure: fronted, as for /I/ (F2 locus
// near 1800 Hz).
constexpr Posture kAlveolarLocus = {11.25, 0.900, 4.0, 0.0};
// Dark /l/ around a tip contact at 16.7 cm (tools/lateral_fit, jaw raised):
// 380 920 2730 Hz against bdl's 380 950 2850.
constexpr Posture kDarkL = {8.25, 0.250, 1.00, 0.0};

constexpr PhoneSpec kPhones[] = {
    {"IY", Kind::vowel, true, 155, 55, kIY, kNone, Organ::none, 0, 0, 0.5},
    {"IH", Kind::vowel, true, 135, 40, kIH, kNone, Organ::none, 0, 0, 0.5},
    {"EY", Kind::diphthong, true, 190, 100, kEY, kIH, Organ::none, 0, 0, 0.5},
    {"EH", Kind::vowel, true, 150, 70, kEH, kNone, Organ::none, 0, 0, 0.5},
    {"AE", Kind::vowel, true, 230, 80, kAE, kNone, Organ::none, 0, 0, 0.5},
    {"AA", Kind::vowel, true, 240, 100, kAA, kNone, Organ::none, 0, 0, 0.5},
    {"AO", Kind::vowel, true, 240, 100, kAO, kNone, Organ::none, 0, 0, 0.5},
    {"OW", Kind::diphthong, true, 220, 80, kOW, kUH, Organ::none, 0, 0, 0.5},
    {"UH", Kind::vowel, true, 160, 60, kUH, kNone, Organ::none, 0, 0, 0.5},
    {"UW", Kind::vowel, true, 210, 70, kUW, kNone, Organ::none, 0, 0, 0.5},
    {"AH", Kind::vowel, true, 140, 60, kAH, kNone, Organ::none, 0, 0, 0.5},
    {"AX", Kind::vowel, true, 110, 40, kAX, kNone, Organ::none, 0, 0, 0.5},
    {"ER", Kind::vowel, true, 180, 80, kER, kNone, Organ::dorsum, 13.0, 0.30, 2.5},
    {"AY", Kind::diphthong, true, 250, 150, kAA, kIH, Organ::none, 0, 0, 0.5},
    {"AW", Kind::diphthong, true, 260, 100, kAA, kUH, Organ::none, 0, 0, 0.5},
    {"OY", Kind::diphthong, true, 280, 150, kAO, kIH, Organ::none, 0, 0, 0.5},
    {"P", Kind::stop, false, 90, 50, kNone, kNone, Organ::lips, 17.0, 0.0, 0.5},
    {"B", Kind::stop, true, 85, 60, kNone, kNone, Organ::lips, 17.0, 0.0, 0.5},
    {"T", Kind::stop, false, 75, 50, kNone, kNone, Organ::tip, 16.0, 0.0, 0.8},
    {"D", Kind::stop, true, 75, 50, kNone, kNone, Organ::tip, 16.0, 0.0, 0.8},
    {"DX", Kind::stop, true, 28, 20, kNone, kNone, Organ::tip, 16.0, 0.0, 0.8},
    {"K", Kind::stop, false, 80, 60, kNone, kNone, Organ::dorsum, 12.05, 0.0, 0.5},
    {"G", Kind::stop, true, 80, 60, kNone, kNone, Organ::dorsum, 12.05, 0.0, 0.5},
    {"M", Kind::nasal, true, 70, 60, kNone, kNone, Organ::lips, 17.0, 0.0, 0.5},
    {"N", Kind::nasal, true, 60, 50, kNone, kNone, Organ::tip, 16.0, 0.0, 0.8},
    {"NG", Kind::nasal, true, 95, 60, kNone, kNone, Organ::dorsum, 12.05, 0.0, 0.5},
    {"F", Kind::fricative, false, 100, 80, kNone, kNone, Organ::lips, 17.0, 0.12, 0.5},
    {"V", Kind::fricative, true, 60, 40, kNone, kNone, Organ::lips, 17.0, 0.14, 0.5},
    {"TH", Kind::fricative, false, 90, 60, kNone, kNone, Organ::tip, 16.7, 0.16, 0.6},
    {"DH", Kind::fricative, true, 50, 30, kNone, kNone, Organ::tip, 16.7, 0.18, 0.6},
    {"S", Kind::fricative, false, 105, 60, kNone, kNone, Organ::tip, 15.6, 0.06, 0.8},
    {"Z", Kind::fricative, true, 75, 40, kNone, kNone, Organ::tip, 15.6, 0.07, 0.8},
    {"SH", Kind::fricative, false, 105, 80, kNone, kNone, Organ::tip, 14.8, 0.16, 2.0},
    {"ZH", Kind::fricative, true, 70, 40, kNone, kNone, Organ::tip, 14.8, 0.17, 2.0},
    {"CH", Kind::affricate, false, 70, 50, kNone, kNone, Organ::tip, 15.0, 0.0, 2.0},
    {"JH", Kind::affricate, true, 70, 50, kNone, kNone, Organ::tip, 15.0, 0.0, 2.0},
    {"L", Kind::liquid, true, 80, 40, kDarkL, kNone, Organ::tip, 16.7, 0.18, 0.8},
    {"R", Kind::liquid, true, 80, 30, kER, kNone, Organ::dorsum, 13.0, 0.22, 2.5},
    {"W", Kind::glide, true, 80, 60, kUW, kNone, Organ::lips, 17.0, 0.20, 0.5},
    {"Y", Kind::glide, true, 80, 40, kIY, kNone, Organ::none, 0, 0, 0.5},
    {"HH", Kind::aspirate, false, 80, 20, kNone, kNone, Organ::glottis, 0, 0, 0.5},
};

// The phone inventory in use: the built-in table, optionally overridden
// from a data file (so a tuner, or another language, can change it without
// recompiling).
std::vector<PhoneSpec>& phone_table() {
  static std::vector<PhoneSpec> table(std::begin(kPhones), std::end(kPhones));
  return table;
}

const PhoneSpec& spec_of(const std::string& name) {
  const std::vector<PhoneSpec>& table = phone_table();
  for (const PhoneSpec& spec : table) {
    if (name == spec.name) {
      return spec;
    }
  }
  return table[11];  // AX
}

bool is_vowel_kind(Kind kind) {
  return kind == Kind::vowel || kind == Kind::diphthong;
}

bool is_front_vowel(const std::string& name) {
  return name == "IY" || name == "IH" || name == "EY" || name == "EH" || name == "AE" || name == "Y";
}

bool is_obstruent(Kind kind) {
  return kind == Kind::stop || kind == Kind::fricative || kind == Kind::affricate;
}

// A stop released into a sonorant bursts inside its own interval: the
// segment is closure, burst and aspiration, and the next phone begins at
// voice onset. Median burst-to-end in the bdl labels (tools/burst_position):
// p 24, t 18, k 43 ms; b 10, d 8, g 28 ms (the velar lag is longest, as in
// every English VOT survey). Released at the segment's end, our /k/ and /t/
// slices at 75% were 10-20 dB short of the human aspiration and the
// four-view identification heard them as /p/.
double post_burst_interval(Organ organ, bool voiced) {
  switch (organ) {
    case Organ::lips:
      return voiced ? 0.010 : 0.024;
    case Organ::tip:
      return voiced ? 0.008 : 0.018;
    case Organ::dorsum:
      return voiced ? 0.028 : 0.043;
    case Organ::glottis:
    case Organ::none:
      break;
  }
  return 0.0;
}

// Modal voicing (`PlanOptions::voicing_gap`): the folds held slightly
// pressed (negative rest gap) so they meet mid-swing and snap shut. At 0
// they once glided into the midline and left 2-4 kHz 8 dB short; since the
// junction-loss and source fixes, -0.015 was too pressed instead: the low
// harmonics sat 3-4 dB under bdl's long-term spectrum, H1-H2 -4.6 dB
// (bdl -3.7), shimmer 8.4% (7.3). -0.008 matches all three without
// costing identification (-0.004 and 0 did: 36.2 and 34.6%).
constexpr double kConnectedVowelScale = 0.5;
constexpr double kSealLead = 0.4;

double lung_pressure_for(double pitch, const PlanOptions& options) {
  double pressure = options.speech_pressure + 5.0 * (pitch - 125.0);
  return std::clamp(pressure, 560.0, 1000.0);
}
constexpr double kVoicelessGap = 0.14;
constexpr double kWaitingGap = 0.04;  // folds apart, not vibrating, before a phrase-initial voiced closure
constexpr double kVoiceBarLead = 0.08;
constexpr double kExpansionRelax = 0.10;
constexpr double kClosureVoicingGap = 0.0;
constexpr double kUnreleasedOverlap = 0.02;
constexpr double kLabialRelease = 0.04;
constexpr double kReleaseClear = 0.05;
constexpr double kBreathGap = 0.45;

struct Emitter {
  Score* score = nullptr;

  void part(double time, Part part, double a, double b, double move) const {
    PartEvent event;
    event.time = std::max(time, 0.0);
    event.part = part;
    event.a = static_cast<float>(a);
    event.b = static_cast<float>(b);
    event.move = static_cast<float>(move);
    score->parts.push_back(event);
  }

  void mark(double time, Mark mark, double a, const std::string& text) const {
    MarkEvent event;
    event.time = std::max(time, 0.0);
    event.mark = mark;
    event.a = static_cast<float>(a);
    event.text = text;
    score->marks.push_back(event);
  }

  void vowel_shape(double time, const Posture& posture, double move) const {
    part(time, Part::tongue, posture.place, posture.area, move);
    part(time, Part::lips, posture.lips, posture.protrusion, move);
  }
};

// Durations in the manner of Klatt (1979): each phone has an inherent and a
// minimum duration; rules scale the compressible part.
// The planner's own random sequence (xorshift64*, Box-Muller): seeded per
// utterance, so a given seed always speaks the same way.
struct Variation {
  std::uint64_t state = 0x9E3779B97F4A7C15ULL;

  explicit Variation(std::uint32_t seed) : state(0x9E3779B97F4A7C15ULL ^ (static_cast<std::uint64_t>(seed) << 17U)) {}

  double uniform() {
    state ^= state >> 12U;
    state ^= state << 25U;
    state ^= state >> 27U;
    std::uint64_t mixed = state * 0x2545F4914F6CDD1DULL;
    return (static_cast<double>(mixed >> 11U) + 0.5) * (1.0 / 9007199254740992.0);
  }

  double gauss() {
    double radius = std::sqrt(-2.0 * std::log(uniform()));
    return radius * std::cos(2.0 * kPi * uniform());
  }
};

double model_factor(const FactorModel& model, const std::string& name) {
  std::map<std::string, double>::const_iterator found = model.factors.find(name);
  return found == model.factors.end() ? 0.0 : found->second;
}

// The learned duration model's class for a neighbouring sound.
const char* duration_class(const std::string& name) {
  static const char* const kVoicedObstruents[] = {"B", "D", "G", "V", "DH", "Z", "ZH", "JH"};
  static const char* const kVoicelessObstruents[] = {"P", "T", "K", "F", "TH", "S", "SH", "CH", "HH"};
  for (const char* voiced : kVoicedObstruents) {
    if (name == voiced) {
      return "voiced";
    }
  }
  for (const char* voiceless : kVoicelessObstruents) {
    if (name == voiceless) {
      return "voiceless";
    }
  }
  return is_vowel_phone(name) ? "vowel" : "sonorant";
}

std::string lowercase(const std::string& text) {
  std::string result = text;
  for (char& letter : result) {
    letter = static_cast<char>(std::tolower(static_cast<unsigned char>(letter)));
  }
  return result;
}

// Seconds for one phone from the learned model (tools/duration_fit.py on
// rms: held-out log error 0.412 against the hand rules' 0.458). The
// factors and their names mirror duration_fit.py's context().
double learned_seconds(const std::vector<TimedPhone>& phones, std::size_t first, std::size_t index,
                       const Phrase& phrase, const FactorModel& model) {
  const TimedPhone& phone = phones[index];
  bool vowel = is_vowel_phone(phone.name);
  std::string boundary = phrase.end == Boundary::minor ? "comma" : "end";
  double sum = model.intercept + model_factor(model, "phone:" + lowercase(phone.name));
  if (phone.phrase_final) {
    sum += model_factor(model, (vowel ? "final_vowel:" : "final_consonant:") + boundary);
    if (!vowel && index + 1 == phones.size()) {
      sum += model_factor(model, "last_consonant:" + boundary);
    }
  }
  std::size_t word_first = index;
  while (word_first > first && !phones[word_first].word_start) {
    --word_first;
  }
  std::size_t word_last = index;
  while (word_last + 1 < phones.size() && !phones[word_last].word_end) {
    ++word_last;
  }
  int syllables = 0;
  bool earlier_vowel = false;
  bool later_vowel = false;
  for (std::size_t other = word_first; other <= word_last; ++other) {
    if (is_vowel_phone(phones[other].name)) {
      ++syllables;
      earlier_vowel = earlier_vowel || other < index;
      later_vowel = later_vowel || other > index;
    }
  }
  sum += model_factor(model, "syllables:" + std::to_string(std::min(syllables, 3)));
  if (vowel) {
    sum += model_factor(model, "stress:" + std::to_string(phone.stress));
    if (phone.nuclear) {
      sum += model_factor(model, "nuclear");
    } else if (phone.accent) {
      sum += model_factor(model, "accent");
    }
    const char* position = !earlier_vowel && !later_vowel ? "only"
                           : !earlier_vowel               ? "first"
                           : !later_vowel                 ? "last"
                                                          : "middle";
    sum += model_factor(model, std::string("syllable_position:") + position);
    std::string following = index < word_last ? duration_class(phones[index + 1].name) : "word_end";
    sum += model_factor(model, "vowel_before:" + following);
  } else {
    const char* position = phone.word_start ? "initial" : phone.word_end ? "final" : "medial";
    sum += model_factor(model, std::string("consonant_position:") + position);
    bool before = index > word_first && !is_vowel_phone(phones[index - 1].name);
    bool after = index < word_last && !is_vowel_phone(phones[index + 1].name);
    if (before || after) {
      sum += model_factor(model, "cluster");
    }
    if (index + 1 < phones.size() && is_vowel_phone(phones[index + 1].name) && phones[index + 1].stress == 1) {
      sum += model_factor(model, "before_stress");
    }
  }
  return std::exp(sum);
}

void time_phrase(const Phrase& phrase, std::vector<TimedPhone>& phones, double start, const PlanOptions& options,
                 Variation& variation) {
  std::size_t first = phones.size();
  for (const Word& word : phrase.words) {
    for (std::size_t index = 0; index < word.phones.size(); ++index) {
      TimedPhone phone;
      phone.name = word.phones[index].name;
      phone.stress = word.phones[index].stress;
      phone.word_start = index == 0;
      phone.word_end = index + 1 == word.phones.size();
      phone.accent = word.accent && phone.stress == 1;
      phone.nuclear = word.nuclear && phone.stress == 1;
      phone.word = word.text;
      phone.vary_duration = variation.gauss();
      phone.vary_pitch = variation.gauss();
      phone.vary_place = variation.gauss();
      phone.vary_area = variation.gauss();
      phone.vary_lips = variation.gauss();
      phone.vary_effort = variation.gauss();
      phones.push_back(phone);
    }
  }
  // American English flapping: /t/ or /d/ between a vowel and an unstressed
  // vowel becomes a quick voiced tap ("water", "at a").
  for (std::size_t index = first + 1; index + 1 < phones.size(); ++index) {
    TimedPhone& phone = phones[index];
    bool coronal = phone.name == "T" || phone.name == "D";
    bool after_vowel = is_vowel_phone(phones[index - 1].name) || phones[index - 1].name == "R";
    bool before_weak = is_vowel_phone(phones[index + 1].name) && phones[index + 1].stress == 0;
    if (coronal && after_vowel && before_weak) {
      phone.phoneme = phone.name;
      phone.name = "DX";
    }
  }
  // The last vowel of the phrase and everything after it is phrase-final.
  for (std::size_t index = phones.size(); index > first; --index) {
    phones[index - 1].phrase_final = true;
    if (is_vowel_phone(phones[index - 1].name)) {
      break;
    }
  }
  double time = start;
  for (std::size_t index = first; index < phones.size(); ++index) {
    TimedPhone& phone = phones[index];
    const PhoneSpec& spec = spec_of(phone.name);
    double factor = 1.0;
    bool previous_consonant = index > first && !is_vowel_phone(phones[index - 1].name) && !phone.word_start;
    bool next_consonant = index + 1 < phones.size() && !is_vowel_phone(phones[index + 1].name) &&
                          !phones[index + 1].word_start;
    if (is_vowel_kind(spec.kind)) {
      if (phone.stress == 0) {
        factor *= 0.65;
      } else if (phone.stress == 2) {
        factor *= 0.8;
      }
      if (phone.nuclear) {
        factor *= 1.25;
      } else if (phone.accent) {
        factor *= 1.1;
      }
      if (index + 1 < phones.size() && !phones[index + 1].word_start) {
        const PhoneSpec& next = spec_of(phones[index + 1].name);
        if (is_obstruent(next.kind)) {
          factor *= next.voiced ? 1.15 : 0.8;
        }
      }
      if (!phone.phrase_final) {
        bool word_has_later_vowel = false;
        for (std::size_t later = index + 1; later < phones.size() && !phones[later].word_start; ++later) {
          if (is_vowel_phone(phones[later].name)) {
            word_has_later_vowel = true;
          }
        }
        if (word_has_later_vowel) {
          factor *= 0.85;
        }
      }
    } else {
      if (previous_consonant || next_consonant) {
        factor *= 0.75;
      }
      bool before_stress = index + 1 < phones.size() && phones[index + 1].stress == 1;
      factor *= before_stress ? 1.1 : 0.9;
    }
    if (options.duration_model.loaded) {
      double seconds = learned_seconds(phones, first, index, phrase, options.duration_model);
      if (phone.name == "DX") {
        // The duration model has no flap (bdl's labels write a flap as /t/
        // or /d/, so its factor never existed): flaps took a generic stop's
        // 80-95 ms ("body", "water") where a tap lasts 20-30 ms. Its own
        // inherent length instead.
        seconds = spec.inherent_ms * 1e-3;
      }
      seconds = std::max(seconds * std::exp(options.variation * phone.vary_duration), 0.02);
      phone.start = time;
      phone.duration = seconds / options.rate;
      time += phone.duration;
      continue;
    }
    double milliseconds = spec.minimum_ms + (spec.inherent_ms - spec.minimum_ms) * factor;
    if (is_vowel_kind(spec.kind)) {
      // Klatt's inherent vowel durations describe careful citation speech;
      // read connected speech (CMU ARCTIC bdl) runs vowels at about 0.45x.
      milliseconds *= kConnectedVowelScale;
    }
    if (phone.phrase_final) {
      // Lengthening on the whole duration: applied to the stretchable part
      // only (Klatt), the last schwa of "Papa," lasted 48 ms (rms 85 before
      // a pause) and was cut off like a /t/ ("papat").
      // A list item or clause before a comma is held longest (continuation).
      double vowel_scale = phrase.end == Boundary::minor ? options.comma_vowel_scale : options.final_vowel_scale;
      milliseconds *= is_vowel_kind(spec.kind) ? vowel_scale : options.final_consonant_scale;
      // Short vowels stretch further: the final schwa of "Delta." held 70
      // ms behind its /t/ and the word lost its second syllable.
      if (is_vowel_kind(spec.kind)) {
        milliseconds = std::max(milliseconds, options.final_vowel_floor_ms);
      }
    }
    // No two tokens are the same length (log-normal, sd = variation).
    milliseconds = std::max(milliseconds * std::exp(options.variation * phone.vary_duration), spec.minimum_ms * 0.8);
    phone.start = time;
    phone.duration = milliseconds / 1000.0 / options.rate;
    time += phone.duration;
  }
}

// Intonation: a declining baseline per phrase, accent peaks on stressed
// syllables of accented words, and a boundary tone at the end.
// The learned model's name for the sound before a vowel.
const char* preceding_class(const std::vector<TimedPhone>& phones, std::size_t index) {
  if (index == 0) {
    return "start";
  }
  const std::string& name = phones[index - 1].name;
  if (is_vowel_phone(name)) {
    return "vowel";
  }
  const PhoneSpec& spec = spec_of(name);
  if (is_obstruent(spec.kind) || spec.kind == Kind::aspirate) {
    return spec.voiced ? "voiced" : "voiceless";
  }
  return "sonorant";
}

// Vowel pitch from the learned model: each vowel's semitones re the median
// are the sum of its factors. Fitted on rms (1131 sentences), it predicts
// held-out vowel pitch with median error 1.15 semitones (correlation 0.62);
// the hand rules below scored 1.50 (0.40), worse than a flat line (1.38):
// their accents added 22-26 Hz (3-3.5 semitones) where the speaker's add
// 0.4, and their declination ran by time instead of through the phrase.
void intone_phrase_learned(std::vector<TimedPhone>& phones, std::size_t first, const Phrase& phrase,
                           int phrase_in_sentence, const PlanOptions& options) {
  const FactorModel& model = options.pitch_model;
  std::vector<std::size_t> vowels;
  for (std::size_t index = first; index < phones.size(); ++index) {
    if (is_vowel_phone(phones[index].name)) {
      vowels.push_back(index);
    }
  }
  bool statement_end = phrase.end != Boundary::minor;
  for (std::size_t order = 0; order < vowels.size(); ++order) {
    TimedPhone& phone = phones[vowels[order]];
    double semitones = 0.0;
    semitones += model_factor(model, "stress:" + std::to_string(std::max(phone.stress, 0)));
    if (phone.nuclear) {
      semitones += model_factor(model, "nuclear");
    } else if (phone.accent) {
      semitones += model_factor(model, "accent");
    }
    double progress = static_cast<double>(order) / static_cast<double>(std::max<std::size_t>(vowels.size() - 1, 1));
    int quarter = std::min(static_cast<int>(progress * 4.0), 3);
    semitones += model_factor(model, "phrase_progress:" + std::to_string(quarter));
    semitones += model_factor(model, "phrase_number:" + std::to_string(std::min(phrase_in_sentence, 2)));
    if (order == 0) {
      semitones += model_factor(model, "phrase_first_vowel");
    }
    if (phone.phrase_final) {
      std::string boundary = statement_end ? "final:end" : "final:comma";
      semitones += model_factor(model, phone.nuclear ? boundary + ":nuclear" : boundary);
    }
    semitones += model_factor(model, std::string("after:") + preceding_class(phones, vowels[order]));
    semitones = model.intercept + options.pitch_spread * semitones;
    // No two syllables are pitched alike (sd about 0.5 semitone at 10%).
    semitones += 5.0 * options.variation * phone.vary_pitch;
    phone.pitch = std::max(options.pitch_median * std::pow(2.0, semitones / 12.0), options.lowest_pitch / options.final_fall);
    if (phone.phrase_final && phrase.end == Boundary::question) {
      phone.pitch *= 1.05;
    }
  }
}

void intone_phrase(std::vector<TimedPhone>& phones, std::size_t first, const Phrase& phrase, int phrase_in_sentence,
                   const PlanOptions& options) {
  if (first >= phones.size()) {
    return;
  }
  if (options.pitch_model.loaded) {
    intone_phrase_learned(phones, first, phrase, phrase_in_sentence, options);
    return;
  }
  double start = phones[first].start;
  double end = phones.back().start + phones.back().duration;
  double length = std::max(end - start, 0.2);
  double top = options.base_pitch * (1.0 - 0.05 * std::min(phrase_in_sentence, 3));
  for (std::size_t index = first; index < phones.size(); ++index) {
    TimedPhone& phone = phones[index];
    if (!is_vowel_phone(phone.name)) {
      continue;
    }
    double middle = phone.start + 0.5 * phone.duration;
    double progress = (middle - start) / length;
    double baseline = top * (1.0 - 0.14 * progress);
    double pitch = baseline;
    if (phone.accent) {
      // Accent heights differ from one accent to the next.
      double excursion = phone.nuclear ? options.accent * 1.2 : options.accent;
      pitch += excursion * (1.0 + 2.0 * options.variation * phone.vary_pitch);
    } else {
      pitch *= 1.0 + 0.3 * options.variation * phone.vary_pitch;
    }
    if (phone.phrase_final) {
      if (phrase.end == Boundary::major) {
        pitch = phone.nuclear ? baseline + 0.6 * options.accent : baseline * 0.95;
      } else if (phrase.end == Boundary::question) {
        pitch = baseline * 1.05;
      }
    }
    phone.pitch = pitch;
  }
}

double next_vowel_index(const std::vector<TimedPhone>& phones, std::size_t index, std::size_t end) {
  for (std::size_t later = index; later < end; ++later) {
    if (is_vowel_phone(phones[later].name)) {
      return static_cast<double>(later);
    }
  }
  return -1.0;
}

Posture vowel_posture_near(const std::vector<TimedPhone>& phones, std::size_t index, std::size_t first,
                           std::size_t end) {
  double next = next_vowel_index(phones, index, end);
  if (next >= 0.0) {
    return spec_of(phones[static_cast<std::size_t>(next)].name).vowel;
  }
  for (std::size_t earlier = index; earlier > first; --earlier) {
    if (is_vowel_phone(phones[earlier - 1].name)) {
      const PhoneSpec& spec = spec_of(phones[earlier - 1].name);
      return spec.kind == Kind::diphthong ? spec.glide_to : spec.vowel;
    }
  }
  return kAX;
}

// A stop's closure lets go in stages. The tissue first peels apart slowly
// (PlanOptions::release_peel_area over release_peel_time), then the channel
// opens, then the articulator clears. Opened at full speed from contact, the
// /t/ burst rang the tract at 0.8-1.2 kHz (identical with every noise
// source off); a 60 ms release lowered that ring 5 dB but smeared the burst.
void emit_stop_release(const Emitter& emit, Part part, double place, double release_time, double release_move,
                       const PlanOptions& options) {
  double channel_time = release_time;
  if (options.release_peel_area > 0.0) {
    emit.part(release_time, part, place, options.release_peel_area, options.release_peel_time);
    channel_time = release_time + options.release_peel_time;
  }
  emit.part(channel_time, part, place, options.release_channel, release_move);
  emit.part(channel_time + release_move, part, place, kReleased, kReleaseClear);
}

double velar_place(const std::vector<TimedPhone>& phones, std::size_t index, std::size_t first, std::size_t end,
                   const PlanOptions& options) {
  double next = next_vowel_index(phones, index, end);
  std::string neighbour;
  // A velar that closes a syllable (right after its vowel, not right
  // before one) is shaped by that vowel: in "quick brown" the next vowel
  // is /aU/, and placed by it the tongue slid back from /I/ into a closure
  // 0.85 cm behind its rest; F2 fell into the closure as for a labial
  // (the vowel before a velar raised F2-F3 apart, +25 Hz, where bdl's pinch
  // together, -46 Hz) and the owner heard "quip".
  bool after_vowel = index > first && is_vowel_phone(phones[index - 1].name);
  bool before_vowel = index + 1 < end && is_vowel_phone(phones[index + 1].name);
  if (after_vowel && !before_vowel) {
    neighbour = phones[index - 1].name;
  } else if (next >= 0.0) {
    neighbour = phones[static_cast<std::size_t>(next)].name;
  } else if (index > first) {
    neighbour = phones[index - 1].name;
  }
  // The dorsum contacts further forward next to front vowels: beside /i/
  // it closes nearly where the vowel narrows (13.75 cm). Fronted only
  // 0.85 cm, the tongue slid back 1.7 cm into the closure, F2 fell as for a
  // labial, and the owner heard "speak" as "speap".
  double centre = spec_of(phones[index].name).constriction_place;
  if (!is_front_vowel(neighbour)) {
    return centre - 0.85;
  }
  // Close behind the front vowel's own narrowing, as far forward as the
  // vowel goes: a fixed 1.8 cm put the /k/ of "quick" (after /I/ at
  // 11.25 cm) ahead of its vowel, and the tongue slid forward into it.
  double vowel_place = spec_of(neighbour).vowel.place;
  return std::clamp(vowel_place - 0.3, centre + 0.85, centre + options.velar_front);
}

// Writes the gestures of one phrase.
void gesture_phrase(const std::vector<TimedPhone>& phones, std::size_t first, std::size_t end, const Phrase& phrase,
                    const PlanOptions& options, const Emitter& emit) {
  if (first >= end) {
    return;
  }
  double phrase_start = phones[first].start;
  double phrase_end = phones[end - 1].start + phones[end - 1].duration;

  // Breath support: pressure rises just before speech, carries stressed
  // syllables a little higher, and is released at the end.
  // A phrase that opens on a vowel has no consonant to cover the pressure's
  // rise: with 40 ms of lead, "Echo" and "and" after a comma began with the
  // folds barely moving and were unvoiced through their 80 ms vowels.
  double lead = is_vowel_phone(phones[first].name) ? std::max(options.breath_lead, options.vowel_breath_lead)
                                                     : options.breath_lead;
  emit.part(phrase_start - lead, Part::lungs, options.speech_pressure, 0.0, options.breath_rise);
  emit.mark(phrase_start, Mark::phrase, 0.0, "");
  for (std::size_t look = first; look < end; ++look) {
    if (is_vowel_phone(phones[look].name)) {
      emit.part(phrase_start - 0.2, Part::pitch, phones[look].pitch, 0.0, 0.1);
      break;
    }
  }

  // Voicing posture before the first sound.
  const PhoneSpec& first_spec = spec_of(phones[first].name);
  bool first_voiceless = !first_spec.voiced;
  // A phrase-initial /b d g dZ/ is not prevoiced in English: the folds
  // wait slightly apart and voicing starts at the release. Voiced through
  // the closure, "Juliett" opened with 90 ms of hum before the word (NATO's
  // starts with the hiss).
  bool first_unvoiced_closure = first_spec.voiced && (first_spec.kind == Kind::stop || first_spec.kind == Kind::affricate);
  // Soft onset: voiced phrases adduct before the pressure rises.
  double onset_gap = first_voiceless ? kVoicelessGap : first_unvoiced_closure ? kWaitingGap : options.voicing_gap;
  emit.part(phrase_start - (first_voiceless ? 0.06 : 0.11), Part::glottis, onset_gap, 0.6, 0.05);

  // The first vowel's shape is set up during the breath.
  Posture opening = vowel_posture_near(phones, first, first, end);
  emit.vowel_shape(phrase_start - 0.15, opening, 0.10);
  emit.part(phrase_start - 0.15, Part::velum, 0.0, 0.0, 0.08);

  // Gliding intonation: anchor pitches per vowel. Stressed vowels and the
  // phrase's first and last keep their planned pitch; the others lie on the
  // line (in semitones, by time) between their kept neighbours, a little
  // lower.
  std::vector<double> glide_pitch(phones.size(), 0.0);
  double glide_from = phrase_start - 0.1;
  if (options.pitch_glide) {
    std::vector<std::size_t> kept;
    std::size_t last_vowel = end;
    for (std::size_t look = first; look < end; ++look) {
      if (is_vowel_phone(phones[look].name)) {
        last_vowel = look;
      }
    }
    for (std::size_t look = first; look < end; ++look) {
      if (!is_vowel_phone(phones[look].name)) {
        continue;
      }
      glide_pitch[look] = phones[look].pitch;
      if (phones[look].stress == 1 || kept.empty() || look == last_vowel) {
        kept.push_back(look);
      }
    }
    for (std::size_t pair = 0; pair + 1 < kept.size(); ++pair) {
      std::size_t left = kept[pair];
      std::size_t right = kept[pair + 1];
      double left_time = phones[left].start + 0.5 * phones[left].duration;
      double right_time = phones[right].start + 0.5 * phones[right].duration;
      double left_tone = std::log2(phones[left].pitch);
      double right_tone = std::log2(phones[right].pitch);
      for (std::size_t look = left + 1; look < right; ++look) {
        if (!is_vowel_phone(phones[look].name)) {
          continue;
        }
        double time = phones[look].start + 0.5 * phones[look].duration;
        double share = std::clamp((time - left_time) / std::max(right_time - left_time, 1.0e-3), 0.0, 1.0);
        double tone = left_tone + share * (right_tone - left_tone) - options.unstressed_dip / 12.0;
        glide_pitch[look] = std::pow(2.0, tone);
      }
    }
  }

  double tongue_free_at = phrase_start - 0.15;
  bool velum_open = false;
  for (std::size_t index = first; index < end; ++index) {
    const TimedPhone& phone = phones[index];
    const PhoneSpec& spec = spec_of(phone.name);
    double start = phone.start;
    double finish = phone.start + phone.duration;
    if (phone.word_start) {
      emit.mark(start, Mark::word, 0.0, phone.word);
    }
    bool next_exists = index + 1 < end;
    const PhoneSpec& next_spec = next_exists ? spec_of(phones[index + 1].name) : spec_of("AX");

    if (is_vowel_kind(spec.kind)) {
      // The tongue body heads for this vowel from the start of the syllable
      // onset (or as soon as a consonant releases it).
      std::size_t onset = index;
      while (onset > first && !is_vowel_phone(phones[onset - 1].name)) {
        --onset;
      }
      double begin = std::max(phones[onset].start - 0.02, tongue_free_at);
      // A tip fricative in the onset keeps the body out of its way: heading
      // for /i/ during /sh/, the tongue made a second, palatal constriction
      // behind the hiss channel, the flow fell to about 200 cm^3/s, and the
      // /sh/ of "She ..." sat 23 dB under its vowel (rms: 8 dB).
      for (std::size_t look = onset; look < index; ++look) {
        const PhoneSpec& onset_spec = spec_of(phones[look].name);
        if (onset_spec.kind == Kind::fricative && onset_spec.organ == Organ::tip) {
          begin = std::max(begin, phones[look].start + options.fricative_hold * phones[look].duration);
        }
      }
      // The tongue is still travelling when the vowel ends: it aims to
      // arrive late in the vowel, so formants move through the syllable.
      double reach = std::clamp(start + options.vowel_arrival * phone.duration - begin, 0.05, options.vowel_reach);
      // Each token aims a little differently (undershoot, overshoot).
      double aimed_place = spec.vowel.place + 3.0 * options.variation * phone.vary_place;
      double aimed_area = spec.vowel.area * std::exp(options.variation * phone.vary_area);
      double aimed_lips = spec.vowel.lips * std::exp(0.7 * options.variation * phone.vary_lips);
      emit.part(begin, Part::tongue, aimed_place, aimed_area, reach);
      bool lips_busy = false;
      for (std::size_t look = onset; look < index; ++look) {
        Organ organ = spec_of(phones[look].name).organ;
        if (organ == Organ::lips) {
          lips_busy = true;
        }
      }
      if (!lips_busy) {
        emit.part(begin, Part::lips, aimed_lips, spec.vowel.protrusion, reach);
      }
      if (spec.dorsum_area > 0.0) {
        // After a velar the dorsum is busy until the vowel itself begins.
        bool dorsum_busy = false;
        for (std::size_t look = onset; look < index; ++look) {
          if (spec_of(phones[look].name).organ == Organ::dorsum) {
            dorsum_busy = true;
          }
        }
        double dorsum_begin = dorsum_busy ? start : begin;
        emit.part(dorsum_begin, Part::dorsum, spec.dorsum_place, spec.dorsum_area,
                  std::max(start + 0.5 * phone.duration - dorsum_begin, 0.05));
        emit.part(finish - 0.03, Part::dorsum, spec.dorsum_place, kReleased, 0.08);
      }
      if (spec.kind == Kind::diphthong) {
        double glide_begin = start + 0.35 * phone.duration;
        double glide_move = std::max(0.6 * phone.duration, 0.08);
        emit.part(glide_begin, Part::tongue, spec.glide_to.place + 3.0 * options.variation * phone.vary_place,
                  spec.glide_to.area * std::exp(options.variation * phone.vary_area), glide_move);
        emit.part(glide_begin, Part::lips, spec.glide_to.lips, spec.glide_to.protrusion, glide_move);
      }
      // (REJECTED: a "relax" re-target of the tongue at 65% of each
      // monophthong. With no drift left it only re-issued the vowel's own
      // posture, and after a vowel shorter than about 57 ms it landed after
      // the next vowel had begun and pulled the tongue back: the /a/ of
      // "Papa" after the 56 ms /I/ of "is" stayed at F1 431, F2 2104 Hz.)
      if (phone.name == "ER") {
        emit.part(begin, Part::dorsum, spec.constriction_place, spec.constriction_area, reach);
        emit.part(begin, Part::floor, spec.floor, 0.0, reach);
        emit.part(finish - 0.03, Part::dorsum, 13.0, kReleased, 0.08);
        emit.part(finish - 0.03, Part::floor, 0.5, 0.0, 0.1);
      }
      // Pitch: a sloped syllable target the voice chases (after Xu's target
      // approximation). Level targets per syllable, reached in 50-150 ms,
      // gave a staircase of flat harmonics in the spectrogram: early
      // synthetic speech. A stressed syllable rises from just before its
      // vowel to a peak 40% into it and falls toward the next vowel; an
      // unstressed one drifts down. Phrase-final vowels keep their
      // boundary tones below.
      if (options.pitch_glide) {
        // Head for this vowel's anchor from the last one's (see glide_pitch).
        double anchor = start + (phone.stress == 1 ? 0.4 : 0.5) * phone.duration;
        emit.part(glide_from, Part::pitch, glide_pitch[index], 0.0, std::max(anchor - glide_from, 0.05));
        glide_from = anchor;
      } else if (options.accent_peak > 0.0) {
        double next_vowel = finish + 0.08;
        for (std::size_t look = index + 1; look < end; ++look) {
          if (is_vowel_phone(phones[look].name)) {
            next_vowel = phones[look].start;
            break;
          }
        }
        if (next_vowel - finish > 0.4) {
          next_vowel = finish + 0.08;
        }
        if (phone.stress == 0) {
          emit.part(start, Part::pitch, phone.pitch, 0.0, 0.6 * phone.duration);
          if (!phone.phrase_final) {
            emit.part(start + 0.6 * phone.duration, Part::pitch, phone.pitch * (1.0 - options.unstressed_drift), 0.0,
                      std::max(next_vowel - start - 0.6 * phone.duration, 0.05));
          }
        } else {
          double scale = phone.stress == 1 ? 1.0 : 0.5;
          double top = start + 0.4 * phone.duration;
          emit.part(start - 0.10, Part::pitch, phone.pitch * (1.0 + scale * options.accent_peak), 0.0, top - start + 0.10);
          if (!phone.phrase_final) {
            emit.part(top, Part::pitch, phone.pitch * (1.0 - scale * options.accent_fall), 0.0,
                      std::max(next_vowel - top, 0.05));
          }
        }
      } else {
        emit.part(start, Part::pitch, phone.pitch, 0.0, std::clamp(0.6 * phone.duration, 0.05, 0.15));
      }
      if (phone.phrase_final && phrase.end != Boundary::minor) {
        // Boundary tone: a statement falls to the floor through its last
        // vowel; a question rises.
        // (With the learned model the vowel already sits low; it falls a
        // little further, not to a fixed floor that could lie above it.)
        double statement_tone = options.pitch_model.loaded
                                    ? std::max(phone.pitch * options.final_fall, options.lowest_pitch)
                                    : options.pitch_floor;
        double tone = phrase.end == Boundary::question ? phone.pitch * 1.45 : statement_tone;
        emit.part(start + 0.3 * phone.duration, Part::pitch, tone, 0.0, 0.8 * phone.duration);
      } else if (phone.phrase_final) {
        // Before a comma the voice falls a little, not up: over the last
        // vowel before a comma pause rms falls 1.3 semitones (median; 52%
        // fall, 13% rise). A 6% rise made every comma a question.
        emit.part(start + 0.4 * phone.duration, Part::pitch, phone.pitch * options.comma_pitch, 0.0, 0.6 * phone.duration);
      }
      // Voice quality: firmer folds on stressed syllables, a relaxed and
      // breathier release through the last vowel of a phrase.
      bool after_voiceless = index > first && !spec_of(phones[index - 1].name).voiced;
      if (!after_voiceless) {
        double quality = phone.stress == 1   ? options.voicing_gap - options.stress_press
                         : phone.stress == 0 ? options.voicing_gap + options.unstressed_gap
                                             : options.voicing_gap;
        double compression = phone.stress == 0 ? options.unstressed_compression : options.vowel_compression;
        emit.part(start - 0.03, Part::glottis, quality, compression, 0.06);
      }
      if (phone.phrase_final) {
        emit.part(start + 0.5 * phone.duration, Part::glottis, options.final_relax_gap, 0.4, 0.5 * phone.duration);
      }
      // Breath support follows pitch: higher, stressed syllables are driven
      // harder, low ones gently (subglottal pressure and F0 rise together).
      // Each syllable is a pulse of breath: pressure swells into the vowel
      // and decays through it, so loudness rises and falls instead of
      // holding a plateau.
      // Breath effort differs from syllable to syllable.
      double peak_pressure = lung_pressure_for(phone.pitch, options) * (1.0 + 0.6 * options.variation * phone.vary_effort);
      if (phone.stress == 0) {
        peak_pressure *= options.unstressed_pressure;
      }
      emit.part(start - 0.02, Part::lungs, peak_pressure, 0.0, 0.35 * phone.duration + 0.03);
      emit.part(start + 0.4 * phone.duration, Part::lungs, 0.88 * peak_pressure, 0.0, 0.6 * phone.duration);
      // Nasalization by anticipation of a following nasal.
      if (next_exists && next_spec.kind == Kind::nasal && !velum_open) {
        emit.part(start + 0.5 * phone.duration, Part::velum, 0.45, 0.0, 0.10);
        velum_open = true;
      }
      continue;
    }

    // Consonants.
    Posture context = vowel_posture_near(phones, index, first, end);
    double place = spec.constriction_place;
    if (spec.organ == Organ::dorsum && phone.name != "R") {
      place = velar_place(phones, index, first, end, options);
    }
    bool next_same_organ = next_exists && next_spec.organ == spec.organ && !is_vowel_kind(next_spec.kind);
    // Tip and lip gestures are fast (about 25-30 ms); the dorsum is heavier.
    // Each aims slightly past its target so the constriction is formed by
    // the phone's start, as real articulators overshoot onto contact.
    double close_move = spec.organ == Organ::dorsum ? 0.04 : 0.028;
    double release_move = spec.kind == Kind::stop || spec.kind == Kind::affricate ? options.release_move : 0.04;
    // Voiced stops let go more gently: at the voiceless speed their bursts
    // stood 4-6.5 dB above bdl's and rms's (/b/ -0.3 dB re the vowel; theirs
    // -6.5 and -7.1).
    if (spec.kind == Kind::stop && spec.voiced) {
      release_move = options.voiced_release_move;
    }
    if (phone.name == "L") {
      // A lateral lets go abruptly: the tip snaps from the ridge and F1
      // jumps within 10-20 ms (the cue that /l/ began a syllable).
      release_move = 0.018;
    }
    if (spec.organ == Organ::lips && (spec.kind == Kind::stop || spec.kind == Kind::nasal)) {
      // Lips part more slowly than the tip: F1 and F2 rise out of a labial
      // over 30-50 ms. At 22 ms F2 started 1430 Hz into /a/ ("bog"), barely
      // below the alveolar 1570, and the owner heard "dog" as "bag".
      release_move = std::max(kLabialRelease, release_move);
    }
    // A closure must be sealed by the phone's start: contact happens before
    // the gesture's end (it aims past contact), so stops and nasals start a
    // full movement time early; constrictions that stay open start later.
    bool sealing = spec.kind == Kind::stop || spec.kind == Kind::nasal || spec.kind == Kind::affricate;
    double close_time = start - (sealing ? kSealLead : 0.5) * close_move;
    double aimed_area = spec.kind == Kind::fricative ? 0.6 * spec.constriction_area : spec.constriction_area;
    if (index == first) {
      // The first consonant of a phrase is in place before the air flows.
      close_time = start - 0.12;
    }
    // The moment the closure opens: a stop into a vowel, liquid or glide
    // bursts before its segment ends; everything else opens at the end.
    double opening = finish;
    bool released_into_sonorant = spec.kind == Kind::stop && next_exists &&
                                  (is_vowel_kind(next_spec.kind) || next_spec.kind == Kind::liquid ||
                                   next_spec.kind == Kind::glide);
    if (released_into_sonorant) {
      opening = finish - std::min(options.burst_scale * post_burst_interval(spec.organ, spec.voiced), 0.4 * phone.duration);
    }
    if (phone.name == "DX") {
      // A flap is a ballistic tap: the tip flicks onto the ridge and off
      // again, with full contact for its short duration. With a stop's
      // 28 ms approach and 50 ms voiced release, launched 16 ms apart, the
      // tip never touched (0.21 cm^2 at its closest) and "body" was heard
      // as "belly" (the owner).
      close_move = 0.012;
      close_time = start - 0.5 * close_move;
      release_move = 0.015;
      opening = finish;
    }
    double release_time = opening - 0.25 * release_move;
    // A stop before another stop or affricate made elsewhere ("quick brown",
    // "that dog") is unreleased in English: it opens only once the next
    // closure has formed. Released first, it burst into silence and the
    // owner heard "quick brown" as "uplep prawn".
    if (spec.kind == Kind::stop && next_exists && (next_spec.kind == Kind::stop || next_spec.kind == Kind::affricate) &&
        next_spec.organ != spec.organ) {
      release_time = finish + kUnreleasedOverlap;
    }

    switch (spec.organ) {
      case Organ::lips:
        emit.part(close_time, Part::lips, aimed_area, phone.name == "W" ? 0.6 : 0.0, close_move);
        if (!next_same_organ) {
          bool into_sibilant = next_exists && (phones[index + 1].name == "S" || phones[index + 1].name == "Z");
          double open_lips = into_sibilant ? std::max(context.lips, options.sibilant_lips) : context.lips;
          double open_protrusion = into_sibilant ? std::min(context.protrusion, options.sibilant_protrusion)
                                                 : context.protrusion;
          emit.part(release_time, Part::lips, open_lips, open_protrusion, release_move);
        }
        break;
      case Organ::tip:
        emit.part(close_time, Part::tip, place, aimed_area, close_move);
        if (spec.kind == Kind::affricate) {
          // The release opens to the /sh/ channel (its fricative aim), not
          // wider: at 0.165 cm^2 the jet was too slow to hiss, and /dZ/
          // went straight from its burst into the vowel.
          // The hiss is the /sh/ fricative's own channel and place (its aim:
          // 0.6 x area), so the two stay one gesture when /sh/ is refitted.
          // Copied once as fixed numbers (14.8 cm, 0.072 cm^2), they fell
          // behind the /sh/ fit and /dZ/ hissed flat to 8 kHz like /f/ (the
          // owner heard "jumps" as "femps").
          const PhoneSpec& hush = spec_of("SH");
          double channel = 0.6 * hush.constriction_area;
          // Closure about 40% of the affricate, hiss the rest (the NATO
          // /dZ/: a short stop, then 50-60 ms of voiced /sh/-like noise).
          emit.part(start + options.affricate_hiss_at * phone.duration, Part::tip, hush.constriction_place, channel, 0.02);
          if (spec.voiced) {
            // The folds part so the hiss dominates, as in English /dZ/.
            emit.part(start + (options.affricate_hiss_at - 0.05) * phone.duration, Part::glottis, options.affricate_hiss_gap, 0.4, 0.03);
            // Re-close only for a voiced sound: issued 10 ms before the end,
            // it overrode a following voiceless sound's opening (30 ms early).
            if (next_exists && next_spec.voiced) {
              emit.part(finish - 0.01, Part::glottis, options.voicing_gap, 0.6, 0.03);
            }
          }
        }
        if (!next_same_organ) {
          if (spec.kind == Kind::stop) {
            // A stop opens a narrow channel first and clears afterwards:
            // aimed straight at "released", the tip crossed contact at full
            // speed and the opening jumped (sharper bursts than bdl's).
            emit_stop_release(emit, Part::tip, place, release_time, release_move, options);
          } else {
            emit.part(release_time, Part::tip, place, kReleased, release_move);
          }
        }
        break;
      case Organ::dorsum:
        emit.part(close_time, Part::dorsum, place, aimed_area, close_move);
        if (!next_same_organ) {
          if (spec.kind == Kind::stop) {
            emit_stop_release(emit, Part::dorsum, place, release_time, release_move, options);
          } else {
            emit.part(release_time, Part::dorsum, place, kReleased, release_move);
          }
        }
        break;
      case Organ::glottis:
      case Organ::none:
        break;
    }

    // Alveolar stops and nasals: the tongue body is drawn forward behind the
    // tip contact (the alveolar F2 locus near 1800 Hz) and returns to the
    // vowel after the release.
    bool alveolar = spec.organ == Organ::tip && spec.constriction_place >= 15.0 && spec.constriction_area <= 0.0 &&
                    (spec.kind == Kind::stop || spec.kind == Kind::nasal);
    Posture closure_body = alveolar ? kAlveolarLocus : context;
    if (alveolar) {
      emit.part(close_time, Part::tongue, closure_body.place, closure_body.area, close_move + 0.02);
      if (!(spec.kind == Kind::stop && spec.voiced)) {
        emit.part(opening + 0.005, Part::tongue, context.place, context.area, 0.06);
      }
    }
    // Voiced stops: the speaker enlarges the pharynx during the closure
    // (tongue root forward, larynx down) so air keeps flowing through the
    // folds; the next vowel's tongue target is restored at the release. The
    // enlargement runs through the last 80 ms before the release, while
    // the lungs push: done at the closure itself, a phrase-initial /d/ had
    // finished expanding before any air flowed, and its closure was silent
    // where the NATO speaker's carries a 20-30 ms voice bar.
    if (spec.kind == Kind::stop && spec.voiced && !next_exists) {
      // A voiced stop ending its phrase devoices at the release: voiced
      // through it, "dog." released into a puff of schwa.
      emit.part(opening - 0.01, Part::glottis, kVoicelessGap, 0.4, 0.03);
    }
    if ((spec.kind == Kind::stop || spec.kind == Kind::affricate) && spec.voiced) {
      double expand_begin = std::max(close_time, opening - kVoiceBarLead);
      emit.part(expand_begin, Part::tongue, closure_body.place, std::min(closure_body.area * 1.8, 4.0),
                std::max(opening - expand_begin, 0.03));
      // The folds slacken to touching so air keeps passing them: held
      // pressed, the enlargement sucked the closed cavity below atmospheric
      // and the release imploded (air rushed in through the lips).
      if (index == first) {
        emit.part(opening - 0.005, Part::glottis, options.voicing_gap, 0.6, options.adduction_move);
      } else if (spec.kind == Kind::stop) {
        emit.part(expand_begin, Part::glottis, kClosureVoicingGap, 0.6, 0.03);
      }
      // The enlargement relaxes after the release, not with it: collapsed in
      // 60 ms it pushed a gust out through the opening stop (the largest
      // spikes, at /g/ releases).
      emit.part(opening + 0.01, Part::tongue, context.place, context.area, kExpansionRelax);
    }
    // While the tip holds a fricative channel, the body cannot form a
    // second tight one behind it: set up early for the /i/ of "She", the
    // body's palatal channel (0.18 cm^2 over 3 cm) choked the flow to about
    // 200 cm^3/s and the /sh/ sat 23 dB under its vowel (rms: 8 dB). The
    // body heads for the vowel, no narrower than fricative_body_area, and
    // takes the vowel's own aperture as the fricative ends.
    bool tip_channel = spec.organ == Organ::tip && (spec.kind == Kind::fricative || spec.kind == Kind::affricate);
    if (tip_channel && options.fricative_body_area > 0.0 && context.area < options.fricative_body_area) {
      emit.part(close_time, Part::tongue, context.place, options.fricative_body_area, close_move);
      if (next_exists && is_vowel_kind(next_spec.kind)) {
        emit.part(finish - 0.02, Part::tongue, context.place, context.area, 0.06);
      }
    }
    // Sibilants: the jaw rises so the jet strikes the incisors.
    bool alveolar_sibilant = phone.name == "S" || phone.name == "Z";
    bool postalveolar = phone.name == "SH" || phone.name == "ZH" || phone.name == "CH" || phone.name == "JH";
    if (alveolar_sibilant || postalveolar) {
      emit.part(close_time, Part::jaw, alveolar_sibilant ? 0.16 : 0.35, 1.0, 0.06);
      bool next_sibilant = next_exists && (phones[index + 1].name == "S" || phones[index + 1].name == "Z" ||
                                           phones[index + 1].name == "SH" || phones[index + 1].name == "ZH" ||
                                           phones[index + 1].name == "CH" || phones[index + 1].name == "JH");
      if (!next_sibilant) {
        emit.part(release_time, Part::jaw, kReleased, next_exists && is_vowel_kind(next_spec.kind) ? 0.0 : 1.0, 0.07);
      }
    } else if (spec.organ != Organ::glottis && spec.kind != Kind::glide && spec.kind != Kind::liquid) {
      // (Not the liquids: the tongue makes /l r/ with the jaw left open.
      // Raised for the /l/ of "Alfa", it closed the /ae/ before it to F1
      // 420-600 Hz where its posture gives 750, and the owner heard "elfa".)
      // Every other oral consonant raises the jaw; the next vowel drops it.
      // The lips close on their own, so a labial lifts the jaw only a
      // little: fully raised 40 ms ahead, the /p/ of "Papa" closed the jaw
      // over the 86 ms /a/ (F1 431, F2 2104 Hz: an /I/; "body" reached
      // 737/1336), and the owner heard "bipple".
      double jaw_raise = spec.organ == Organ::lips ? options.labial_jaw : 1.0;
      emit.part(start - options.jaw_lead, Part::jaw, kReleased, jaw_raise, 0.06);
      if (!next_exists || is_vowel_kind(next_spec.kind)) {
        emit.part(opening - 0.01, Part::jaw, kReleased, 0.0, 0.08);
      }
    }
    // Rounding for a following /w r/ starts with this consonant's closure
    // ("quick", "tree", "sweet"): rounded only as the /k/ let go, the burst
    // and the aspiration came out through open lips, bright to 5 kHz where
    // rms's /kw/ has almost nothing above 1.5 kHz.
    bool rounded_next = next_exists && (phones[index + 1].name == "W" || phones[index + 1].name == "R");
    if (options.round_ahead && rounded_next && spec.organ != Organ::lips) {
      double lips = phones[index + 1].name == "W" ? next_spec.constriction_area : next_spec.vowel.lips;
      double protrusion = phones[index + 1].name == "W" ? 0.6 : next_spec.vowel.protrusion;
      emit.part(close_time, Part::lips, lips, protrusion, close_move + 0.02);
    }
    // Lip rounding and the space under the blade for postalveolars and /r/.
    if (phone.name == "SH" || phone.name == "ZH" || phone.name == "CH" || phone.name == "JH") {
      emit.part(close_time, Part::lips, options.postalveolar_lips, options.postalveolar_protrusion, 0.06);
      emit.part(release_time, Part::lips, context.lips, context.protrusion, 0.06);
    }
    bool lips_free = index == first || spec_of(phones[index - 1].name).organ != Organ::lips;
    if ((phone.name == "S" || phone.name == "Z") && lips_free &&
        (context.lips < options.sibilant_lips || context.protrusion > options.sibilant_protrusion)) {
      emit.part(close_time, Part::lips, std::max(context.lips, options.sibilant_lips),
                std::min(context.protrusion, options.sibilant_protrusion), 0.06);
      emit.part(release_time, Part::lips, context.lips, context.protrusion, 0.06);
    }
    if (spec.floor > 1.0) {
      emit.part(close_time, Part::floor, spec.floor, 0.0, 0.06);
      emit.part(release_time, Part::floor, 0.5, 0.0, 0.08);
    }
    // Consonants that also shape the tongue body and lips: the liquids
    // (American /l/ is dark in every position: the back of the tongue
    // rises, and its F3 stays high above a deep F2-F3 valley, the cue that
    // keeps it from /r/) and the glides. Their postures come from the
    // phone table (tools/lateral_fit, tools/vowel_fit rhotic).
    if (phone.name == "R" || phone.name == "L" || phone.name == "W" || phone.name == "Y") {
      Posture body = spec.vowel;
      emit.part(start - 0.04, Part::tongue, body.place, body.area, 0.07);
      if (phone.name == "R" || phone.name == "L") {
        emit.part(start - 0.04, Part::lips, body.lips, body.protrusion, 0.07);
      }
      if (phone.name == "R") {
        emit.part(start - 0.04, Part::floor, std::max(spec.floor, 0.5), 0.0, 0.07);
        emit.part(release_time, Part::floor, 0.5, 0.0, 0.1);
      }
      if (phone.name == "Y") {
        emit.part(start - 0.04, Part::lips, 3.0, 0.0, 0.07);
      }
      // A lateral leaves a pocket over the blade behind its contact, which
      // the air bypasses through the side channels: a side resonator.
      if (spec.pocket > 0.0) {
        emit.part(start - 0.03, Part::lateral, spec.pocket, 0.0, 0.04);
        emit.part(release_time, Part::lateral, 0.0, 0.0, release_move);
      }
      tongue_free_at = start + 0.5 * phone.duration;
    }

    // Velum for nasals.
    if (spec.kind == Kind::nasal) {
      if (!velum_open) {
        emit.part(index == first ? start - 0.14 : start - 0.05, Part::velum, 0.45, 0.0, 0.07);
        velum_open = true;
      }
      if (!next_exists || next_spec.kind != Kind::nasal) {
        emit.part(finish - 0.01, Part::velum, 0.0, 0.0, 0.07);
        velum_open = false;
      }
    }

    // Voiced fricatives: the folds part slightly so enough air reaches the
    // constriction to hiss while they keep vibrating.
    if (spec.voiced && spec.kind == Kind::fricative && phone.name != "DH") {
      emit.part(start - 0.02, Part::glottis, options.voiced_fricative_gap, 0.4, 0.03);
      // Before a pause or a consonant, English /z v/ devoice through their
      // second half (bdl's final /z/: hiss to the end, voicing fading).
      bool vowel_next = next_exists && (is_vowel_kind(next_spec.kind) || next_spec.kind == Kind::liquid ||
                                        next_spec.kind == Kind::glide || next_spec.kind == Kind::nasal);
      if (options.final_devoicing > 0.0 && !vowel_next) {
        emit.part(start + 0.4 * phone.duration, Part::glottis, options.final_devoicing, 0.4, 0.5 * phone.duration);
        if (next_exists && next_spec.voiced) {
          emit.part(finish - 0.01, Part::glottis, options.voicing_gap, 0.6, 0.03);
        }
      } else if (next_exists && next_spec.voiced) {
        // Re-closed before a voiceless sound, this landed after that sound
        // had opened the folds (30 ms early) and voiced it: "sells sea" was
        // heard as "sells zea".
        emit.part(finish - 0.01, Part::glottis, options.voicing_gap, 0.6, 0.03);
      }
    }
    // Glottis: voiceless sounds open the folds; aspirated stops keep them
    // open past the release (voice onset time).
    if (!spec.voiced) {
      double gap = spec.kind == Kind::aspirate ? options.aspirate_gap : kVoicelessGap;
      double compression = spec.kind == Kind::aspirate ? 0.2 : 0.4;
      // Fricatives and /h/ need the air before their constriction forms. A
      // stop's folds part as the mouth seals: parted 30 ms early, the open
      // tract carried /h/-like flow (650 cm^3/s) and the tongue tip closed
      // into that stream with a burst of turbulence.
      double abduct = spec.kind == Kind::stop ? start : start - options.fricative_abduct_lead;
      emit.part(abduct, Part::glottis, gap, compression, options.abduction_move);
      emit.part(start - 0.03, Part::lungs, options.voiceless_pressure * options.speech_pressure, 0.0, 0.06);
      bool next_voiceless = next_exists && !next_spec.voiced;
      // Voicing resumes only if something voiced follows in the phrase; after
      // a phrase-final stop the folds stay apart for the pause.
      if (next_exists && !next_voiceless) {
        double voice_onset = finish;
        bool voiced_obstruent_next = next_exists && is_obstruent(next_spec.kind);
        // After /s/ in the same word, /p t k/ are unaspirated ("speak",
        // "stay", "sky"): the folds close in as the stop releases. Aspirated
        // like a word-initial /p/, the /i/ of "speak" was left weak and
        // breathy and the owner heard "speet".
        bool after_s = index > first && phones[index - 1].name == "S" && !phone.word_start;
        if (spec.kind == Kind::stop && after_s) {
          voice_onset = opening;
        } else if (spec.kind == Kind::stop && !voiced_obstruent_next) {
          bool stressed_next = next_exists && phones[index + 1].stress >= 1;
          // Released inside the segment, the aspiration is already there:
          // voicing begins at the segment's end (a little later into stress).
          voice_onset = released_into_sonorant ? finish + (stressed_next ? options.voice_lag : 0.0)
                                               : finish + (stressed_next ? 0.045 : 0.015);
          // Into /w r l j/ the aspiration runs on through the sonorant
          // ("quick", "play", "tree"): rms voices the /w/ of "quick" about
          // 80 ms after the burst. Voiced from the stop's end, our /w/ came
          // out a bright, voiced glide and the owner heard "clip".
          bool sonorant_next = next_spec.kind == Kind::liquid || next_spec.kind == Kind::glide;
          if (sonorant_next) {
            voice_onset = finish + options.sonorant_devoicing * phones[index + 1].duration;
          }
        } else if (spec.kind == Kind::stop) {
          // A voiced stop or fricative follows ("quick brown"): no
          // aspiration lag, or its closure goes unvoiced and /b/ is heard /p/.
          voice_onset = finish;
        } else if (spec.kind == Kind::fricative) {
          voice_onset = finish - 0.015;
        }
        // The folds close in quickly (20 ms). Over 70 ms the vowel after
        // every voiceless consonant began with chaotic, rough pulses
        // (jitter 5.8%, shimmer 35%; bdl 2.1%, 7%).
        if (options.onset_gap > options.voicing_gap) {
          // Voicing starts with the folds still a little apart, so the
          // first cycles do not seal: driven straight to the pressed rest,
          // the first closure cut 300 cm^3/s of aspiration flow inside a
          // millisecond and a quarter of the vowels after /p t k/ began
          // with a click (first pulse > 1.5x the vowel's later peak: ours
          // 23%, bdl none).
          emit.part(voice_onset - options.adduction_lead, Part::glottis, options.onset_gap, 0.6, options.adduction_move);
          emit.part(voice_onset - options.adduction_lead + options.onset_hold, Part::glottis, options.voicing_gap, 0.6,
                    options.onset_settle);
        } else {
          emit.part(voice_onset - options.adduction_lead, Part::glottis, options.voicing_gap, 0.6, options.adduction_move);
        }
      }
    }
  }

  // Phrase end: release pressure; open the glottis for the pause.
  // The phrase ends by running down, not by stopping: breath pressure
  // decays through the last syllable while the folds keep vibrating, and
  // they part only once the sound has faded.
  double last_vowel_start = phrase_end - 0.1;
  for (std::size_t look = end; look > first; --look) {
    if (is_vowel_phone(phones[look - 1].name)) {
      last_vowel_start = phones[look - 1].start;
      break;
    }
  }
  // The fade must finish with the last sound: carried 120 ms past it at
  // nearly full voice and then cut, the tail of "Bravo," and "Echo," was
  // heard as another vowel ("bravua", "echoa").
  // (From 45% into the last vowel, voicing stopped when the pressure fell
  // under the folds' threshold, near the vowel's middle: our "hello," was
  // voiced for 170 of its 308 ms where the owner's carries to its end.)
  double fade_start = last_vowel_start + options.phrase_fade_from * (phrase_end - last_vowel_start);
  double fade_time = phrase_end - fade_start + options.phrase_fade_tail;
  // A phrase that ends in a voiced stop or affricate keeps its breath until
  // the release: faded through the vowel, "dog." released no burst at all,
  // and its voiced closure was heard as /l/ ("doll"). A final voiceless
  // stop does not: held to its release, the /t/ of "Foxtrot." burst at
  // full pressure, 7x the level around it ("mechanical sound at end").
  const TimedPhone& final_phone = phones[end - 1];
  const PhoneSpec& final_spec = spec_of(final_phone.name);
  // (For the segment-module body, whose release drains the pressure through
  // the opening constriction in milliseconds, a final voiceless stop also
  // keeps its breath to the release: faded through the vowel, "Quebec."
  // released at 1400 dyn/cm^2 and no burst was heard, where the owner's
  // final /k/ bursts plainly; `final_stop_breath` keeps it.)
  bool keeps_breath = final_spec.voiced || options.final_stop_breath;
  if ((final_spec.kind == Kind::stop || final_spec.kind == Kind::affricate) && keeps_breath) {
    fade_start = phrase_end;
    fade_time = options.phrase_fade_tail + 0.04;
  }
  emit.part(fade_start, Part::lungs, phrase.end == Boundary::minor ? options.pause_pressure : 0.0, 0.0, fade_time);
  emit.part(phrase_end + options.phrase_open_delay, Part::glottis, 0.12, 0.3, 0.10);
  if (velum_open) {
    emit.part(phrase_end, Part::velum, 0.0, 0.0, 0.08);
  }
}

}  // namespace

bool load_phone_overrides(const std::string& path) {
  std::ifstream stream(path);
  if (!stream) {
    return false;
  }
  std::string line;
  while (std::getline(stream, line)) {
    std::istringstream fields(line);
    std::string name;
    if (!(fields >> name) || name[0] == '#') {
      continue;
    }
    PhoneSpec* target = nullptr;
    for (PhoneSpec& spec : phone_table()) {
      if (name == spec.name) {
        target = &spec;
      }
    }
    if (target == nullptr) {
      return false;
    }
    std::string key;
    double value = 0.0;
    while (fields >> key >> value) {
      if (key == "place") {
        target->constriction_place = value;
      } else if (key == "area") {
        target->constriction_area = value;
      } else if (key == "floor") {
        target->floor = value;
      } else if (key == "inherent") {
        target->inherent_ms = value;
      } else if (key == "minimum") {
        target->minimum_ms = value;
      } else if (key == "vplace") {
        target->vowel.place = value;
      } else if (key == "varea") {
        target->vowel.area = value;
      } else if (key == "vlips") {
        target->vowel.lips = value;
      } else if (key == "vprotrusion") {
        target->vowel.protrusion = value;
      } else if (key == "dplace") {
        target->dorsum_place = value;
      } else if (key == "darea") {
        target->dorsum_area = value;
      } else if (key == "pocket") {
        target->pocket = value;
      } else if (key == "gplace") {
        target->glide_to.place = value;
      } else if (key == "garea") {
        target->glide_to.area = value;
      } else if (key == "glips") {
        target->glide_to.lips = value;
      } else if (key == "gprotrusion") {
        target->glide_to.protrusion = value;
      } else {
        return false;
      }
    }
  }
  return true;
}

bool load_factor_model(const std::string& path, FactorModel& model) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream fields(line);
    std::string kind;
    fields >> kind;
    if (kind == "intercept" || kind == "tempo") {
      fields >> model.intercept;
    } else if (kind == "factor") {
      std::string name;
      double value = 0.0;
      fields >> name >> value;
      model.factors[name] = value;
    }
  }
  model.loaded = true;
  return true;
}

// Outside timing replaces the planned durations and pitches of one phrase's
// phones, which are then laid end to end again from the phrase's start.
void apply_timing(std::vector<TimedPhone>& phones, std::size_t first, std::size_t& cursor,
                  const PlanOptions& options) {
  if (first >= phones.size()) {
    return;
  }
  double time = phones[first].start;
  for (std::size_t index = first; index < phones.size(); ++index) {
    TimedPhone& phone = phones[index];
    if (cursor < options.timing.size()) {
      const PhoneTiming& given = options.timing[cursor];
      if (given.duration > 0.0) {
        phone.duration = given.duration;
      }
      if (given.pitch > 0.0 && is_vowel_phone(phone.name)) {
        phone.pitch = given.pitch;
      }
    }
    ++cursor;
    phone.start = time;
    time += phone.duration;
  }
}

Plan plan_speech(const std::vector<Phrase>& phrases, const PlanOptions& options) {
  Plan plan;
  Emitter emit;
  emit.score = &plan.score;

  // Rest posture before anything happens.
  emit.part(0.0, Part::lungs, 0.0, 0.0, 0.05);
  emit.part(0.0, Part::glottis, 0.15, 0.3, 0.05);
  emit.part(0.0, Part::pitch, options.base_pitch, 0.0, 0.05);
  emit.vowel_shape(0.0, kAX, 0.05);
  emit.part(0.0, Part::velum, 0.2, 0.0, 0.05);
  emit.part(0.0, Part::tip, 16.0, kReleased, 0.05);
  emit.part(0.0, Part::dorsum, 12.0, kReleased, 0.05);
  emit.part(0.0, Part::floor, 0.5, 0.0, 0.05);
  emit.part(0.0, Part::jaw, kReleased, 0.0, 0.05);

  double time = options.lead_in;
  int phrase_in_sentence = 0;
  std::size_t timing_cursor = 0;
  Variation variation(options.seed);
  for (std::size_t phrase_index = 0; phrase_index < phrases.size(); ++phrase_index) {
    const Phrase& phrase = phrases[phrase_index];
    if (phrase.breath_before) {
      // Inhalation: the folds part, the mouth opens a little and the lungs
      // draw air in; it is heard as a soft breath.
      emit.mark(time, Mark::breath, phrase.breath_seconds, "");
      emit.part(time, Part::glottis, kBreathGap, 0.0, 0.08);
      emit.part(time, Part::lungs, -110.0, 0.0, 0.12);
      emit.part(time, Part::lips, 2.6, 0.0, 0.10);
      emit.part(time + phrase.breath_seconds - 0.12, Part::lungs, 0.0, 0.0, 0.08);
      time += phrase.breath_seconds;
    }
    std::size_t first = plan.phones.size();
    time_phrase(phrase, plan.phones, time + 0.12, options, variation);
    intone_phrase(plan.phones, first, phrase, phrase_in_sentence, options);
    if (!options.timing.empty()) {
      apply_timing(plan.phones, first, timing_cursor, options);
    }
    gesture_phrase(plan.phones, first, plan.phones.size(), phrase, options, emit);
    const TimedPhone& last = plan.phones.back();
    time = last.start + last.duration;
    // The breath is let go as the phrase ends, not after its pause: held
    // through the pause, an open glottis hissed for the whole of it (the
    // final /k/ of "Mike" was followed by 220 ms of 300 cm^3/s breath
    // noise; the owner's final release lasts about 50 ms). The next phrase
    // sets its own breath at its onset.
    emit.part(time + 0.05, Part::lungs, 0.0, 0.0, 0.12);
    double pause = 0.0;
    if (phrase.end == Boundary::minor) {
      pause = 0.10;
      ++phrase_in_sentence;
    } else {
      pause = 0.22;
      phrase_in_sentence = 0;
    }
    if (phrase_index < options.pauses.size()) {
      pause = options.pauses[phrase_index];
    }
    time += pause;
  }
  emit.part(time, Part::lungs, 0.0, 0.0, 0.1);
  emit.part(time, Part::glottis, 0.15, 0.3, 0.1);
  emit.vowel_shape(time, kAX, 0.2);
  emit.mark(time + 0.35, Mark::stop, 0.0, "");
  plan.score.duration = time + 0.35;
  sort_score(plan.score);
  return plan;
}

}  // namespace voice
