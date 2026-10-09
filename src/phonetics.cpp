#include "ducklips/phonetics.hpp"
#include <array>
#include <utility>
#include <cctype>
#include <fstream>
#include <sstream>

namespace ducklips {
namespace {

ConstrictionSpec lips_closure() {
    return ConstrictionSpec{.present = true, .slot = ConstrictionSlot::lips, .place_cm = 0.0, .width_cm = 0.45,
                            .area_cm2 = 0.0, .obstacle_cm = -1.0};
}
ConstrictionSpec tip_closure(double place) {
    return ConstrictionSpec{.present = true, .slot = ConstrictionSlot::tongue_tip, .place_cm = place, .width_cm = 0.45,
                            .area_cm2 = 0.0, .obstacle_cm = -1.0};
}
// An alveolar stop: its release jet strikes the upper incisors like a
// sibilant's, which makes the sharp high-frequency burst of /t/ and /d/.
ConstrictionSpec alveolar_stop() {
    ConstrictionSpec spec = tip_closure(2.2);
    spec.obstacle_cm = 0.9;
    return spec;
}
ConstrictionSpec body_closure() {
    return ConstrictionSpec{.present = true, .slot = ConstrictionSlot::tongue_body, .place_cm = 6.0, .width_cm = 0.9,
                            .area_cm2 = 0.0, .obstacle_cm = -1.0};
}
ConstrictionSpec narrowing(ConstrictionSlot slot, double place, double width, double area, double obstacle) {
    return ConstrictionSpec{.present = true, .slot = slot, .place_cm = place, .width_cm = width,
                            .area_cm2 = area, .obstacle_cm = obstacle};
}

PhoneSpec vowel(const char* name, bool tense, double f1, double f2, double f3, double lips, double protrusion,
                Phone glide) {
    PhoneSpec spec{};
    spec.name = name;
    spec.kind = PhoneClass::vowel;
    spec.voiced = true;
    spec.tense = tense;
    spec.f1 = f1;
    spec.f2 = f2;
    spec.f3 = f3;
    spec.lip_area_cm2 = lips;
    spec.lip_protrusion_cm = protrusion;
    spec.glide_to = glide;
    return spec;
}

PhoneSpec consonant(const char* name, PhoneClass kind, bool voiced, ConstrictionSpec primary,
                    ConstrictionSpec secondary) {
    PhoneSpec spec{};
    spec.name = name;
    spec.kind = kind;
    spec.voiced = voiced;
    spec.primary = primary;
    spec.secondary = secondary;
    return spec;
}

std::array<PhoneSpec, phone_count> build_table() {
    const ConstrictionSpec none{};
    const ConstrictionSpec rounding = narrowing(ConstrictionSlot::lips, 0.0, 0.45, 1.4, -1.0);
    std::array<PhoneSpec, phone_count> table{};
    table[static_cast<std::size_t>(Phone::AA)] = vowel("AA", true, 768, 1333, 2522, 6.0, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::AE)] = vowel("AE", false, 588, 1952, 2601, 5.0, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::AH)] = vowel("AH", false, 623, 1200, 2550, 3.5, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::AO)] = vowel("AO", true, 652, 997, 2538, 1.6, 0.5, Phone::count);
    table[static_cast<std::size_t>(Phone::AW)] = vowel("AW", true, 768, 1333, 2522, 6.0, 0.0, Phone::UH);
    table[static_cast<std::size_t>(Phone::AY)] = vowel("AY", true, 768, 1333, 2522, 6.0, 0.0, Phone::IH);
    table[static_cast<std::size_t>(Phone::EH)] = vowel("EH", false, 580, 1799, 2605, 3.5, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::ER)] = vowel("ER", true, 474, 1379, 1710, 1.6, 0.4, Phone::count);
    table[static_cast<std::size_t>(Phone::EY)] = vowel("EY", true, 476, 2089, 2691, 3.2, 0.0, Phone::IH);
    table[static_cast<std::size_t>(Phone::IH)] = vowel("IH", false, 427, 2034, 2684, 3.0, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::IY)] = vowel("IY", true, 342, 2322, 3000, 2.6, 0.0, Phone::count);
    table[static_cast<std::size_t>(Phone::OW)] = vowel("OW", true, 497, 910, 2459, 1.3, 0.6, Phone::UH);
    table[static_cast<std::size_t>(Phone::OY)] = vowel("OY", true, 652, 997, 2538, 1.6, 0.5, Phone::IH);
    table[static_cast<std::size_t>(Phone::UH)] = vowel("UH", false, 469, 1122, 2434, 1.3, 0.5, Phone::count);
    table[static_cast<std::size_t>(Phone::UW)] = vowel("UW", true, 378, 997, 2343, 0.7, 0.9, Phone::count);
    table[static_cast<std::size_t>(Phone::AX)] = vowel("AX", false, 500, 1500, 2500, 3.0, 0.0, Phone::count);

    table[static_cast<std::size_t>(Phone::P)] = consonant("P", PhoneClass::stop, false, lips_closure(), none);
    table[static_cast<std::size_t>(Phone::B)] = consonant("B", PhoneClass::stop, true, lips_closure(), none);
    table[static_cast<std::size_t>(Phone::M)] = consonant("M", PhoneClass::nasal, true, lips_closure(), none);
    table[static_cast<std::size_t>(Phone::T)] = consonant("T", PhoneClass::stop, false, alveolar_stop(), none);
    table[static_cast<std::size_t>(Phone::D)] = consonant("D", PhoneClass::stop, true, alveolar_stop(), none);
    table[static_cast<std::size_t>(Phone::N)] = consonant("N", PhoneClass::nasal, true, tip_closure(2.2), none);
    table[static_cast<std::size_t>(Phone::K)] = consonant("K", PhoneClass::stop, false, body_closure(), none);
    table[static_cast<std::size_t>(Phone::G)] = consonant("G", PhoneClass::stop, true, body_closure(), none);
    table[static_cast<std::size_t>(Phone::NG)] = consonant("NG", PhoneClass::nasal, true, body_closure(), none);
    table[static_cast<std::size_t>(Phone::F)] = consonant("F", PhoneClass::fricative, false,
        narrowing(ConstrictionSlot::lips, 0.0, 0.35, 0.08, -1.0), none);
    table[static_cast<std::size_t>(Phone::V)] = consonant("V", PhoneClass::fricative, true,
        narrowing(ConstrictionSlot::lips, 0.0, 0.35, 0.08, -1.0), none);
    table[static_cast<std::size_t>(Phone::TH)] = consonant("TH", PhoneClass::fricative, false,
        narrowing(ConstrictionSlot::tongue_tip, 1.2, 0.35, 0.09, -1.0), none);
    table[static_cast<std::size_t>(Phone::DH)] = consonant("DH", PhoneClass::fricative, true,
        narrowing(ConstrictionSlot::tongue_tip, 1.2, 0.35, 0.09, -1.0), none);
    // English /s z/ are mostly laminal: the blade grooves against the ridge
    // while the tip rests behind the lower incisors, so no space opens under
    // the tongue. Raised apically (as here), the tip leaves a 1.25 cm^3
    // pocket under it with the floor, which pulls the /s/ noise from 4.3 kHz
    // down to 1.8-3.3 kHz (MEASURED). tip_down = 1 was REJECTED for now: the
    // blade is a straight segment, so with the tip down it cannot bulge up
    // behind the teeth and leaves a large front cavity instead (/s/ at 2.5
    // kHz with or without the floor, MEASURED). It needs a flexible blade.
    ConstrictionSpec alveolar_groove = narrowing(ConstrictionSlot::tongue_tip, 1.8, 0.35, 0.06, 0.9);
    alveolar_groove.tip_down = 0.0;
    table[static_cast<std::size_t>(Phone::S)] = consonant("S", PhoneClass::fricative, false, alveolar_groove, none);
    table[static_cast<std::size_t>(Phone::Z)] = consonant("Z", PhoneClass::fricative, true, alveolar_groove, none);
    // Postalveolars at 2.4 cm from the lips: the /sh/ noise peaks near
    // 2.9 kHz as in the human reading (3.1 kHz). At 3.1 cm it peaked at 1.9 kHz,
    // a dark "sh" heard as Irish.
    table[static_cast<std::size_t>(Phone::SH)] = consonant("SH", PhoneClass::fricative, false,
        narrowing(ConstrictionSlot::tongue_tip, 2.4, 0.55, 0.12, 0.9), rounding);
    table[static_cast<std::size_t>(Phone::ZH)] = consonant("ZH", PhoneClass::fricative, true,
        narrowing(ConstrictionSlot::tongue_tip, 2.4, 0.55, 0.12, 0.9), rounding);
    table[static_cast<std::size_t>(Phone::CH)] = consonant("CH", PhoneClass::affricate, false,
        narrowing(ConstrictionSlot::tongue_tip, 2.4, 0.55, 0.12, 0.9), rounding);
    table[static_cast<std::size_t>(Phone::JH)] = consonant("JH", PhoneClass::affricate, true,
        narrowing(ConstrictionSlot::tongue_tip, 2.4, 0.55, 0.12, 0.9), rounding);
    table[static_cast<std::size_t>(Phone::HH)] = consonant("HH", PhoneClass::aspirate, false, none, none);
    table[static_cast<std::size_t>(Phone::L)] = consonant("L", PhoneClass::liquid, true,
        narrowing(ConstrictionSlot::tongue_tip, 2.2, 0.45, 0.35, -1.0), none);
    // American /r/: a narrow constriction a little further back with firm lip
    // rounding brings F3 down to about 1.6 kHz (MEASURED in "shore" and
    // "aa-r-aa"; the earlier 3.6 cm, 0.30 cm^2, lips 1.1 cm^2 stopped near
    // 2.1 kHz after a back vowel and was heard as an Irish /r/).
    table[static_cast<std::size_t>(Phone::R)] = consonant("R", PhoneClass::liquid, true,
        narrowing(ConstrictionSlot::tongue_tip, 4.5, 0.7, 0.20, -1.0),
        narrowing(ConstrictionSlot::lips, 0.0, 0.45, 0.6, -1.0));
    table[static_cast<std::size_t>(Phone::W)] = consonant("W", PhoneClass::glide, true,
        narrowing(ConstrictionSlot::lips, 0.0, 0.45, 0.25, -1.0),
        narrowing(ConstrictionSlot::tongue_body, 6.0, 1.0, 0.7, -1.0));
    table[static_cast<std::size_t>(Phone::Y)] = consonant("Y", PhoneClass::glide, true,
        narrowing(ConstrictionSlot::tongue_body, 4.6, 1.2, 0.30, -1.0), none);
    const std::array<std::pair<Phone, TonguePlace>, 22> places{{
        {Phone::TH, TonguePlace::dental}, {Phone::DH, TonguePlace::dental},
        {Phone::T, TonguePlace::alveolar}, {Phone::D, TonguePlace::alveolar}, {Phone::N, TonguePlace::alveolar},
        {Phone::S, TonguePlace::alveolar}, {Phone::Z, TonguePlace::alveolar}, {Phone::L, TonguePlace::alveolar},
        {Phone::SH, TonguePlace::postalveolar}, {Phone::ZH, TonguePlace::postalveolar},
        {Phone::CH, TonguePlace::postalveolar}, {Phone::JH, TonguePlace::postalveolar},
        {Phone::R, TonguePlace::postalveolar}, {Phone::Y, TonguePlace::palatal},
        {Phone::K, TonguePlace::velar}, {Phone::G, TonguePlace::velar}, {Phone::NG, TonguePlace::velar},
        {Phone::W, TonguePlace::velar}, {Phone::P, TonguePlace::none}, {Phone::B, TonguePlace::none},
        {Phone::M, TonguePlace::none}, {Phone::F, TonguePlace::none},
    }};
    for (std::size_t index = 0; index < places.size(); ++index) {
        table[static_cast<std::size_t>(places[index].first)].tongue = places[index].second;
    }
    return table;
}

const std::array<PhoneSpec, phone_count>& table() {
    static const std::array<PhoneSpec, phone_count> specs = build_table();
    return specs;
}

// Core pronunciations in CMU dictionary notation. Hand-entered for this
// study; a full CMU dictionary file can be added with Lexicon::load_cmu_file.
constexpr const char* core_dictionary = R"(A AH0
ABOUT AH0 B AW1 T
AFTER AE1 F T ER0
AGAIN AH0 G EH1 N
AIR EH1 R
ALL AO1 L
ALPHA AE1 L F AH0
ALSO AO1 L S OW0
ALWAYS AO1 L W EY2 Z
AM AE1 M
AN AE1 N
AND AH0 N D
ANY EH1 N IY0
ARE AA1 R
AS AE1 Z
AT AE1 T
AWAY AH0 W EY1
BACK B AE1 K
BE B IY1
BEAUTIFUL B Y UW1 T AH0 F AH0 L
BECAUSE B IH0 K AO1 Z
BEEN B IH1 N
BEFORE B IH0 F AO1 R
BEING B IY1 IH0 NG
BETTER B EH1 T ER0
BIG B IH1 G
BLUE B L UW1
BODY B AA1 D IY0
BOOK B UH1 K
BOTH B OW1 TH
BOY B OY1
BRAVO B R AA1 V OW2
BREATH B R EH1 TH
BREATHE B R IY1 DH
BRING B R IH1 NG
BROWN B R AW1 N
BUT B AH1 T
BY B AY1
CALL K AO1 L
CAME K EY1 M
CAN K AE1 N
CHARLIE CH AA1 R L IY0
CHILD CH AY1 L D
CITY S IH1 T IY0
COME K AH1 M
COMPUTER K AH0 M P Y UW1 T ER0
COULD K UH1 D
DAY D EY1
DELTA D EH1 L T AH0
DID D IH1 D
DO D UW1
DOES D AH1 Z
DOG D AO1 G
DON'T D OW1 N T
DOWN D AW1 N
DUCK D AH1 K
DUCKLIPS D AH1 K L IH2 P S
EACH IY1 CH
EARTH ER1 TH
ECHO EH1 K OW0
EIGHT EY1 T
ENGINE EH1 N JH AH0 N
EVEN IY1 V IH0 N
EVERY EH1 V ER0 IY0
EYES AY1 Z
FACE F EY1 S
FAR F AA1 R
FAST F AE1 S T
FATHER F AA1 DH ER0
FEW F Y UW1
FIND F AY1 N D
FIRST F ER1 S T
FISH F IH1 SH
FIVE F AY1 V
FOR F AO1 R
FOUR F AO1 R
FOX F AA1 K S
FOXTROT F AA1 K S T R AA2 T
FROM F R AH1 M
GET G EH1 T
GIVE G IH1 V
GO G OW1
GOLF G AA1 L F
GOOD G UH1 D
GREAT G R EY1 T
GREEN G R IY1 N
HAD HH AE1 D
HAND HH AE1 N D
HAS HH AE1 Z
HAVE HH AE1 V
HE HH IY1
HEAR HH IY1 R
HELLO HH AH0 L OW1
HELP HH EH1 L P
HER HH ER1
HERE HH IY1 R
HIM HH IH1 M
HIS HH IH1 Z
HOME HH OW1 M
HOTEL HH OW0 T EH1 L
HOUSE HH AW1 S
HOW HH AW1
HUMAN HH Y UW1 M AH0 N
I AY1
I'M AY1 M
IF IH1 F
IN IH0 N
INDIA IH1 N D IY0 AH0
INTO IH0 N T UW1
IS IH1 Z
IT IH1 T
IT'S IH1 T S
ITS IH1 T S
JULIET JH UW1 L IY0 EH2 T
JUMPED JH AH1 M P T
JUST JH AH1 S T
KILO K IY1 L OW2
KNOW N OW1
LANGUAGE L AE1 NG G W AH0 JH
LARGE L AA1 R JH
LAST L AE1 S T
LAZY L EY1 Z IY0
LIGHT L AY1 T
LIKE L AY1 K
LIMA L IY1 M AH0
LIPS L IH1 P S
LISTEN L IH1 S AH0 N
LITTLE L IH1 T AH0 L
LIVE L IH1 V
LONG L AO1 NG
LOOK L UH1 K
MACHINE M AH0 SH IY1 N
MADE M EY1 D
MAKE M EY1 K
MAN M AE1 N
MANY M EH1 N IY0
MAY M EY1
ME M IY1
MIGHT M AY1 T
MIKE M AY1 K
MODEL M AA1 D AH0 L
MORE M AO1 R
MORNING M AO1 R N IH0 NG
MOST M OW1 S T
MOTHER M AH1 DH ER0
MOUTH M AW1 TH
MOVE M UW1 V
MOVES M UW1 V Z
MUCH M AH1 CH
MUSIC M Y UW1 Z IH0 K
MUST M AH1 S T
MY M AY1
NAME N EY1 M
NEED N IY1 D
NEVER N EH1 V ER0
NEW N UW1
NIGHT N AY1 T
NINE N AY1 N
NO N OW1
NOSE N OW1 Z
NOT N AA1 T
NOVEMBER N OW0 V EH1 M B ER0
NOW N AW1
NUMBER N AH1 M B ER0
OF AH1 V
OFF AO1 F
OLD OW1 L D
ON AA1 N
ONE W AH1 N
ONLY OW1 N L IY0
OPEN OW1 P AH0 N
OR AO1 R
OSCAR AA1 S K ER0
OTHER AH1 DH ER0
OUR AW1 ER0
OUT AW1 T
OVER OW1 V ER0
PAPA P AA1 P AH0
PEOPLE P IY1 P AH0 L
PERSON P ER1 S AH0 N
PHYSICAL F IH1 Z IH0 K AH0 L
PHYSICALLY F IH1 Z IH0 K L IY0
PHYSICS F IH1 Z IH0 K S
PLACE P L EY1 S
PLEASE P L IY1 Z
QUEBEC K W AH0 B EH1 K
QUICK K W IH1 K
QUIETLY K W AY1 AH0 T L IY0
RAIN R EY1 N
READ R IY1 D
REAL R IY1 L
RED R EH1 D
RIGHT R AY1 T
ROBOT R OW1 B AA2 T
ROMEO R OW1 M IY0 OW0
SAID S EH1 D
SAME S EY1 M
SAW S AO1
SAY S EY1
SEA S IY1
SEASHELL S IY1 SH EH2 L
SEASHELLS S IY1 SH EH2 L Z
SEASHORE S IY1 SH AO2 R
SEE S IY1
SELLS S EH1 L Z
SEVEN S EH1 V AH0 N
SHE SH IY1
SHELLS SH EH1 L Z
SHORE SH AO1 R
SHOULD SH UH1 D
SHUT SH AH1 T
SIERRA S IY0 EH1 R AH0
SIMULATED S IH1 M Y AH0 L EY2 T IH0 D
SING S IH1 NG
SIX S IH1 K S
SLEEPING S L IY1 P IH0 NG
SMALL S M AO1 L
SO S OW1
SOME S AH1 M
SOMETHING S AH1 M TH IH0 NG
SOUND S AW1 N D
SPEAK S P IY1 K
SPEAKER S P IY1 K ER0
SPEAKING S P IY1 K IH0 NG
SPEECH S P IY1 CH
STILL S T IH1 L
SUCH S AH1 CH
SUN S AH1 N
TAKE T EY1 K
TALK T AO1 K
TALKING T AO1 K IH0 NG
TANGO T AE1 NG G OW0
TEETH T IY1 TH
TELL T EH1 L
TEN T EH1 N
TEXT T EH1 K S T
THAN DH AE1 N
THANK TH AE1 NG K
THAT DH AE1 T
THE DH AH0
THEIR DH EH1 R
THEM DH EH1 M
THEN DH EH1 N
THERE DH EH1 R
THESE DH IY1 Z
THEY DH EY1
THING TH IH1 NG
THINK TH IH1 NG K
THIS DH IH1 S
THOSE DH OW1 Z
THREE TH R IY1
THROAT TH R OW1 T
THROUGH TH R UW1
TIME T AY1 M
TO T UW1
TODAY T AH0 D EY1
TONGUE T AH1 NG
TOO T UW1
TWO T UW1
UNDER AH1 N D ER0
UNIFORM Y UW1 N AH0 F AO2 R M
UP AH1 P
US AH1 S
USE Y UW1 Z
VERY V EH1 R IY0
VICTOR V IH1 K T ER0
VIRTUAL V ER1 CH UW0 AH0 L
VOICE V OY1 S
VOICES V OY1 S IH0 Z
WANT W AA1 N T
WAS W AA1 Z
WATER W AO1 T ER0
WAVE W EY1 V
WAY W EY1
WE W IY1
WELL W EH1 L
WERE W ER1
WHAT W AH1 T
WHEN W EH1 N
WHERE W EH1 R
WHICH W IH1 CH
WHILE W AY1 L
WHISKEY W IH1 S K IY0
WHO HH UW1
WHY W AY1
WILL W IH1 L
WITH W IH1 DH
WORD W ER1 D
WORDS W ER1 D Z
WORK W ER1 K
WORLD W ER1 L D
WOULD W UH1 D
XRAY EH1 K S R EY2
YANKEE Y AE1 NG K IY0
YEAR Y IH1 R
YES Y EH1 S
YOU Y UW1
YOUR Y AO1 R
ZERO Z IH1 R OW0
ZULU Z UW1 L UW0
)";

std::string uppercase(const std::string& text) {
    std::string result = text;
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = static_cast<char>(std::toupper(static_cast<unsigned char>(result[index])));
    }
    return result;
}

bool is_vowel_letter(char letter) noexcept {
    bool vowel_letter = letter == 'a' || letter == 'e' || letter == 'i' || letter == 'o' || letter == 'u';
    return vowel_letter;
}

void append(std::vector<PhoneToken>& output, Phone phone, int stress) {
    output.push_back(PhoneToken{.phone = phone, .stress = stress, .duration_override = 0.0});
}

struct Grapheme final {
    const char* letters{};
    std::array<Phone, 3> phones{};
    std::size_t phone_total{};
};

// Multi-letter spellings, longest first. A deliberately small rule set; its
// errors are expected and are why dictionary coverage matters.
constexpr std::array<Grapheme, 33> graphemes{{
    {"tion", {Phone::SH, Phone::AX, Phone::N}, 3},
    {"sion", {Phone::ZH, Phone::AX, Phone::N}, 3},
    {"eigh", {Phone::EY, Phone::count, Phone::count}, 1},
    {"ough", {Phone::AO, Phone::count, Phone::count}, 1},
    {"igh", {Phone::AY, Phone::count, Phone::count}, 1},
    {"tch", {Phone::CH, Phone::count, Phone::count}, 1},
    {"sch", {Phone::S, Phone::K, Phone::count}, 2},
    {"ch", {Phone::CH, Phone::count, Phone::count}, 1},
    {"sh", {Phone::SH, Phone::count, Phone::count}, 1},
    {"th", {Phone::TH, Phone::count, Phone::count}, 1},
    {"ph", {Phone::F, Phone::count, Phone::count}, 1},
    {"wh", {Phone::W, Phone::count, Phone::count}, 1},
    {"ck", {Phone::K, Phone::count, Phone::count}, 1},
    {"ng", {Phone::NG, Phone::count, Phone::count}, 1},
    {"qu", {Phone::K, Phone::W, Phone::count}, 2},
    {"ee", {Phone::IY, Phone::count, Phone::count}, 1},
    {"ea", {Phone::IY, Phone::count, Phone::count}, 1},
    {"oo", {Phone::UW, Phone::count, Phone::count}, 1},
    {"ou", {Phone::AW, Phone::count, Phone::count}, 1},
    {"ow", {Phone::OW, Phone::count, Phone::count}, 1},
    {"oa", {Phone::OW, Phone::count, Phone::count}, 1},
    {"ai", {Phone::EY, Phone::count, Phone::count}, 1},
    {"ay", {Phone::EY, Phone::count, Phone::count}, 1},
    {"oi", {Phone::OY, Phone::count, Phone::count}, 1},
    {"oy", {Phone::OY, Phone::count, Phone::count}, 1},
    {"au", {Phone::AO, Phone::count, Phone::count}, 1},
    {"aw", {Phone::AO, Phone::count, Phone::count}, 1},
    {"ew", {Phone::UW, Phone::count, Phone::count}, 1},
    {"ar", {Phone::AA, Phone::R, Phone::count}, 2},
    {"er", {Phone::ER, Phone::count, Phone::count}, 1},
    {"ir", {Phone::ER, Phone::count, Phone::count}, 1},
    {"ur", {Phone::ER, Phone::count, Phone::count}, 1},
    {"or", {Phone::AO, Phone::R, Phone::count}, 2},
}};

bool matches_at(const std::string& word, std::size_t position, const char* letters) noexcept {
    std::size_t index = 0;
    while (letters[index] != '\0') {
        if (position + index >= word.size() || word[position + index] != letters[index]) {
            return false;
        }
        index += 1;
    }
    return true;
}

std::size_t letter_count(const char* letters) noexcept {
    std::size_t count = 0;
    while (letters[count] != '\0') {
        count += 1;
    }
    return count;
}

Phone single_consonant(char letter, char next) noexcept {
    switch (letter) {
    case 'b': return Phone::B;
    case 'c':
        if (next == 'e' || next == 'i' || next == 'y') {
            return Phone::S;
        }
        return Phone::K;
    case 'd': return Phone::D;
    case 'f': return Phone::F;
    case 'g': return Phone::G;
    case 'h': return Phone::HH;
    case 'j': return Phone::JH;
    case 'k': return Phone::K;
    case 'l': return Phone::L;
    case 'm': return Phone::M;
    case 'n': return Phone::N;
    case 'p': return Phone::P;
    case 'q': return Phone::K;
    case 'r': return Phone::R;
    case 's': return Phone::S;
    case 't': return Phone::T;
    case 'v': return Phone::V;
    case 'w': return Phone::W;
    case 'z': return Phone::Z;
    default: return Phone::count;
    }
}

Phone vowel_letter_phone(char letter, bool long_vowel) noexcept {
    switch (letter) {
    case 'a': return long_vowel ? Phone::EY : Phone::AE;
    case 'e': return long_vowel ? Phone::IY : Phone::EH;
    case 'i': return long_vowel ? Phone::AY : Phone::IH;
    case 'o': return long_vowel ? Phone::OW : Phone::AA;
    case 'u': return long_vowel ? Phone::UW : Phone::AH;
    default: return Phone::AX;
    }
}

}

const PhoneSpec& phone_spec(Phone phone) noexcept {
    return table()[static_cast<std::size_t>(phone)];
}

bool is_vowel(Phone phone) noexcept {
    bool vowel_phone = phone_spec(phone).kind == PhoneClass::vowel;
    return vowel_phone;
}

bool parse_phone(const std::string& symbol, Phone& phone, int& stress) {
    std::string base = uppercase(symbol);
    int digit = -1;
    if (!base.empty() && std::isdigit(static_cast<unsigned char>(base.back()))) {
        digit = base.back() - '0';
        base.pop_back();
    }
    for (std::size_t index = 0; index < phone_count; ++index) {
        if (base == table()[index].name) {
            phone = static_cast<Phone>(index);
            stress = digit;
            if (is_vowel(phone) && stress < 0) {
                stress = 1;
            }
            if (!is_vowel(phone)) {
                stress = -1;
            }
            return true;
        }
    }
    return false;
}

ConstrictionGesture gesture_from_spec(const ConstrictionSpec& spec) noexcept {
    ConstrictionGesture gesture{};
    gesture.place_from_lips = spec.place_cm * centimetre;
    gesture.width = spec.width_cm * centimetre;
    gesture.target_area = spec.area_cm2 * square_centimetre;
    gesture.obstacle_from_lips = -1.0;
    if (spec.obstacle_cm >= 0.0) {
        gesture.obstacle_from_lips = spec.obstacle_cm * centimetre;
    }
    gesture.tip_down = spec.tip_down;
    return gesture;
}

Lexicon::Lexicon() {
    std::istringstream lines(core_dictionary);
    std::string line{};
    while (std::getline(lines, line)) {
        add_line(line);
    }
}

bool Lexicon::add_line(const std::string& line, bool replace) {
    if (line.empty() || line[0] == ';') {
        return false;
    }
    std::istringstream fields(line);
    std::string word{};
    fields >> word;
    // CMU alternates are written WORD(2); keep the first pronunciation only.
    if (word.empty() || word.find('(') != std::string::npos) {
        return false;
    }
    std::vector<PhoneToken> phones{};
    std::string symbol{};
    while (fields >> symbol) {
        Phone phone = Phone::AX;
        int stress = -1;
        if (!parse_phone(symbol, phone, stress)) {
            return false;
        }
        if (phone == Phone::AH && stress == 0) {
            phone = Phone::AX;
        }
        append(phones, phone, stress);
    }
    std::string key = uppercase(word);
    if (!phones.empty() && (replace || entries_.find(key) == entries_.end())) {
        entries_.insert_or_assign(key, phones);
        return true;
    }
    return false;
}

std::size_t Lexicon::load_cmu_file(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        return 0;
    }
    std::size_t installed = 0;
    std::string line{};
    while (std::getline(file, line)) {
        if (add_line(line, true)) { ++installed; }
    }
    return installed;
}

WordPronunciation Lexicon::pronounce(const std::string& word) const {
    WordPronunciation result{};
    result.word = uppercase(word);
    std::map<std::string, std::vector<PhoneToken>>::const_iterator found = entries_.find(result.word);
    if (found != entries_.end()) {
        result.phones = (*found).second;
        result.from_dictionary = true;
        return result;
    }
    std::string lower = word;
    for (std::size_t index = 0; index < lower.size(); ++index) {
        lower[index] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[index])));
    }
    result.phones = letter_to_sound(lower);
    result.from_dictionary = false;
    return result;
}

std::size_t Lexicon::entry_count() const noexcept {
    return entries_.size();
}

std::vector<PhoneToken> letter_to_sound(const std::string& raw) {
    std::string word{};
    for (std::size_t index = 0; index < raw.size(); ++index) {
        char letter = raw[index];
        if (letter >= 'a' && letter <= 'z') {
            word.push_back(letter);
        }
    }
    std::vector<PhoneToken> output{};
    if (word.empty()) {
        return output;
    }
    std::size_t position = 0;
    if (matches_at(word, 0, "kn") || matches_at(word, 0, "wr")) {
        position = 1;
    }
    // A final silent e lengthens the vowel before a single consonant.
    bool final_e = word.size() >= 3 && word.back() == 'e' && !is_vowel_letter(word[word.size() - 2])
        && is_vowel_letter(word[word.size() - 3]);
    std::size_t end = word.size();
    if (final_e) {
        end -= 1;
    }
    while (position < end) {
        bool matched = false;
        for (std::size_t rule = 0; rule < graphemes.size(); ++rule) {
            const Grapheme& grapheme = graphemes[rule];
            if (matches_at(word, position, grapheme.letters)) {
                for (std::size_t phone = 0; phone < grapheme.phone_total; ++phone) {
                    append(output, grapheme.phones[phone], is_vowel(grapheme.phones[phone]) ? 0 : -1);
                }
                position += letter_count(grapheme.letters);
                matched = true;
                break;
            }
        }
        if (matched) {
            continue;
        }
        char letter = word[position];
        char next = '\0';
        if (position + 1 < word.size()) {
            next = word[position + 1];
        }
        if (is_vowel_letter(letter)) {
            bool long_vowel = final_e && position == end - 2;
            append(output, vowel_letter_phone(letter, long_vowel), 0);
        } else if (letter == 'y') {
            if (position == 0) {
                append(output, Phone::Y, -1);
            } else if (word.size() <= 3) {
                append(output, Phone::AY, 0);
            } else {
                append(output, Phone::IY, 0);
            }
        } else if (letter == 'x') {
            append(output, Phone::K, -1);
            append(output, Phone::S, -1);
        } else if (letter == next) {
            // Doubled consonant letters spell one consonant.
        } else {
            Phone phone = single_consonant(letter, next);
            if (phone != Phone::count) {
                append(output, phone, -1);
            }
        }
        position += 1;
    }
    for (std::size_t index = 0; index < output.size(); ++index) {
        if (output[index].stress >= 0) {
            output[index].stress = 1;
            break;
        }
    }
    return output;
}

}
