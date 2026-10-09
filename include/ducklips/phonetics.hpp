#pragma once
#include "ducklips/articulation.hpp"
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace ducklips {

// ARPAbet phones (CMU dictionary inventory) plus schwa.
enum class Phone : std::size_t {
    AA, AE, AH, AO, AW, AY, EH, ER, EY, IH, IY, OW, OY, UH, UW, AX,
    B, CH, D, DH, F, G, HH, JH, K, L, M, N, NG, P, R, S, SH, T, TH, V, W, Y, Z, ZH,
    count
};
inline constexpr std::size_t phone_count = static_cast<std::size_t>(Phone::count);

enum class PhoneClass { vowel, stop, fricative, affricate, nasal, liquid, glide, aspirate };

// Where a consonant needs the tongue body. A tip or dorsum constriction is
// part of one tongue: to touch the alveolar ridge the whole body moves
// forward, and after release it travels back to the vowel at its own pace.
enum class TonguePlace : std::size_t { none, dental, alveolar, postalveolar, palatal, velar, count };
inline constexpr std::size_t tongue_place_count = static_cast<std::size_t>(TonguePlace::count);

// A constriction in reference units: centimetres and square centimetres.
struct ConstrictionSpec final {
    bool present{false};
    ConstrictionSlot slot{ConstrictionSlot::tongue_tip};
    double place_cm{};
    double width_cm{0.5};
    double area_cm2{};
    double obstacle_cm{-1.0};
    // Tongue-tip gestures: 0 raises the tip itself to the place (apical);
    // 1 leaves the tip down behind the lower incisors while the blade makes
    // the constriction (laminal), which fills the space under the tip.
    double tip_down{0.0};
};

struct PhoneSpec final {
    const char* name{};
    PhoneClass kind{PhoneClass::vowel};
    bool voiced{true};
    bool tense{false};
    // Vowels: formant targets (Hz) for the reference adult male speaker and
    // the lip posture in reference units. Values recalled from Hillenbrand et
    // al. (1995) male averages; not checked against the source in this study.
    double f1{};
    double f2{};
    double f3{};
    double lip_area_cm2{3.0};
    double lip_protrusion_cm{0.0};
    Phone glide_to{Phone::count}; // diphthong offglide, count when none
    TonguePlace tongue{TonguePlace::none}; // tongue-body locus the consonant pulls toward
    ConstrictionSpec primary{};
    ConstrictionSpec secondary{};
};

[[nodiscard]] const PhoneSpec& phone_spec(Phone phone) noexcept;
[[nodiscard]] bool is_vowel(Phone phone) noexcept;
// Parses an ARPAbet symbol with optional stress digit ("AH0", "T").
[[nodiscard]] bool parse_phone(const std::string& symbol, Phone& phone, int& stress);
[[nodiscard]] ConstrictionGesture gesture_from_spec(const ConstrictionSpec& spec) noexcept;

struct PhoneToken final {
    Phone phone{Phone::AX};
    int stress{-1};             // -1 consonant, 0 unstressed, 1 primary, 2 secondary
    double duration_override{}; // seconds; 0 uses the timing rules
    double duration_scale{1.0}; // connected-speech lengthening (geminates)
    bool flap{false};           // /t/ or /d/ realized as a brief voiced tap
};

struct WordPronunciation final {
    std::string word{};
    std::vector<PhoneToken> phones{};
    bool from_dictionary{true};
};

// Pronunciation lookup: a built-in core dictionary, an optional CMU-format
// dictionary file, then letter-to-sound rules as a last resort.
class Lexicon final {
public:
    Lexicon();
    // Reads "WORD  PH0 PH1 ..." lines (CMU dictionary format). Returns the
    // number of entries added; 0 when the file cannot be read.
    std::size_t load_cmu_file(const std::string& path);
    [[nodiscard]] WordPronunciation pronounce(const std::string& word) const;
    [[nodiscard]] std::size_t entry_count() const noexcept;

private:
    void add_line(const std::string& line);
    std::map<std::string, std::vector<PhoneToken>> entries_{};
};

[[nodiscard]] std::vector<PhoneToken> letter_to_sound(const std::string& word);

}
