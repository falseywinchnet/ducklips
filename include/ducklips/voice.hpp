#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/articulation_speeds.hpp"
#include "ducklips/sources.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace ducklips {

// A voice is a body, not a sound: the same behavioural cues drive every voice
// and the differences come from geometry, tissue and habit.
struct VoiceProfile final {
    std::string name{};
    // Geometry. The reference tract (17.5 cm, adult male) defines the neutral
    // area function and consonant places; other lengths scale places by
    // length ratio and areas by its square.
    double tract_length{0.175};       // m, glottis to lips at rest
    double nasal_length{0.115};       // m, velar port to nostrils
    double port_position{0.50};       // fraction of tract length from the glottis
    double nostril_area{1.4e-4};      // m^2, both nostrils
    bool side_cavities{true};         // piriform fossae and maxillary sinuses
    // Fraction of the airway-wall volume velocity that reaches the air through
    // the soft tissue of the neck and cheeks. Calibrated so a voiced closure
    // (voice bar) sits about 25-30 dB under the vowel; not measured.
    double tissue_transmission{0.04};
    // Vocal folds and breathing.
    FoldParameters folds{};
    bool phonation{true};             // false holds the folds still (whisper)
    // Pa, subglottal pressure in conversation (500-700 Pa measured in the
    // literature; 800 sounded pressed and blew the folds open on accented
    // phrase-final vowels).
    double lung_pressure{600.0};
    double adducted_area{3.0e-6};     // m^2 rest glottal area while voicing
    double abducted_area{2.5e-5};     // m^2 for voiceless consonants
    double static_wall_compliance{1.0e-8}; // m^3/Pa, cavity walls during closures
    // Pitch habit.
    double base_pitch{110.0};         // Hz
    // Fractional F0 rise on a primary accent (MEASURED: 0.18 gave half the
    // pitch motion of a human recording).
    double accent_rise{0.27};
    double declination{0.20};         // fractional F0 fall across a phrase
    // Turbulence source gains (dimensionless fraction of jet dynamic pressure).
    // Calibrated against vowel level (see FINDINGS.md), not measured.
    double aspiration_gain{0.0025};
    double frication_gain{0.0045};
    double obstacle_gain{1.5};        // multiplier when the jet strikes the teeth
    // Physiological variation.
    double jitter{0.006};             // rms fractional tension wander
    double shimmer{0.03};             // rms fractional lung pressure wander
    double tremor_rate{5.2};          // Hz
    double tremor_depth{0.004};       // fractional tension
    bool audible_breath{true};        // inhale before each phrase
    double breath_loudness{0.5};      // output gain while inhaling
    double speaking_rate{1.0};        // >1 is faster
    ArticulatorSpeeds speeds{};       // how fast each part of this body can move
    NetworkLosses losses{};
    std::uint64_t seed{1};
};

[[nodiscard]] VoiceProfile baritone_voice();
[[nodiscard]] VoiceProfile alto_voice();
[[nodiscard]] VoiceProfile child_voice();
[[nodiscard]] VoiceProfile gravel_voice();
[[nodiscard]] VoiceProfile whisper_voice();
[[nodiscard]] std::vector<VoiceProfile> all_voices();
// Returns false and leaves output unchanged when the name is unknown.
[[nodiscard]] bool voice_by_name(const std::string& name, VoiceProfile& output);

}
