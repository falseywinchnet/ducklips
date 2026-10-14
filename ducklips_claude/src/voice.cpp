#include "ducklips/voice.hpp"

namespace ducklips {
namespace {

// Length of the vibrating folds relative to the reference adult male's.
// Stiffness and coupling scale with length (force per unit length is kept),
// so the Bernoulli force per unit stiffness, and with it the phonation
// threshold pressure, stays the same.
void scale_fold_length(FoldParameters& folds, double ratio) noexcept {
    folds.length *= ratio;
    folds.lower_stiffness *= ratio;
    folds.upper_stiffness *= ratio;
    folds.coupling_stiffness *= ratio;
}

}

VoiceProfile baritone_voice() {
    VoiceProfile voice{};
    voice.name = "baritone";
    voice.seed = 11;
    return voice;
}

VoiceProfile alto_voice() {
    VoiceProfile voice{};
    voice.name = "alto";
    voice.tract_length = 0.148;
    voice.nasal_length = 0.100;
    voice.nostril_area = 1.1e-4;
    // Lighter vibrating tissue raises the natural pitch. Length and stiffness
    // are kept so the Bernoulli force per unit stiffness, and with it the
    // phonation threshold pressure, stays that of the reference folds.
    voice.folds.lower_mass *= 0.45;
    voice.folds.upper_mass *= 0.45;
    // Shorter folds (about 11 mm); stiffness scales with length, so the
    // phonation threshold pressure is kept.
    scale_fold_length(voice.folds, 0.79);
    voice.lung_pressure = 560.0;
    voice.adducted_area = 3.5e-6;
    voice.abducted_area = 2.0e-5;
    voice.base_pitch = 205.0;
    voice.accent_rise = 0.30;
    voice.aspiration_gain *= 1.6;
    voice.seed = 23;
    return voice;
}

VoiceProfile child_voice() {
    VoiceProfile voice{};
    voice.name = "child";
    // About 27% shorter than the alto's: at 12.2 cm the child was heard as a
    // woman.
    voice.tract_length = 0.108;
    voice.nasal_length = 0.072;
    voice.nostril_area = 0.8e-4;
    voice.folds.lower_mass *= 0.25;
    voice.folds.upper_mass *= 0.25;
    // A child's folds are about 8 mm long: with adult length the flow peaked
    // near 260 cm^3/s, twice a child's (MEASURED in the surface engine).
    scale_fold_length(voice.folds, 0.57);
    voice.lung_pressure = 520.0;
    voice.adducted_area = 2.0e-6;
    voice.abducted_area = 1.5e-5;
    voice.base_pitch = 265.0;
    voice.accent_rise = 0.25;
    voice.aspiration_gain *= 1.4;
    voice.speaking_rate = 0.92;
    // A child's articulators are lighter but less practised.
    voice.speeds.tongue_body = 65.0;
    voice.seed = 37;
    return voice;
}

VoiceProfile gravel_voice() {
    VoiceProfile voice{};
    voice.name = "gravel";
    voice.tract_length = 0.182;
    voice.nasal_length = 0.120;
    voice.lung_pressure = 490.0;
    voice.adducted_area = 0.5e-6;
    voice.base_pitch = 88.0;
    voice.accent_rise = 0.18;
    // Gravel is roughness at the scale of single cycles, from folds of
    // unequal tension, not a slow wobble: the slow wander and tremor at 2.5%
    // and 1% were heard as "queasy".
    voice.folds.asymmetry = 0.7;
    voice.folds.roughness = 0.03;
    voice.jitter = 0.008;
    voice.shimmer = 0.04;
    voice.tremor_depth = 0.003;
    voice.aspiration_gain *= 1.4;
    voice.speaking_rate = 0.88;
    voice.speeds.jaw = 35.0;
    voice.speeds.tongue_body = 60.0;
    voice.seed = 41;
    return voice;
}

VoiceProfile whisper_voice() {
    VoiceProfile voice{};
    voice.name = "whisper";
    voice.tract_length = 0.160;
    voice.nasal_length = 0.105;
    // The folds stay apart and still: there is no voicing, only glottal turbulence.
    voice.phonation = false;
    voice.adducted_area = 1.2e-5;
    voice.abducted_area = 3.0e-5;
    voice.lung_pressure = 400.0;
    voice.base_pitch = 150.0;
    voice.aspiration_gain *= 6.0;
    voice.audible_breath = false;
    voice.seed = 53;
    return voice;
}

std::vector<VoiceProfile> all_voices() {
    std::vector<VoiceProfile> voices{};
    voices.push_back(baritone_voice());
    voices.push_back(alto_voice());
    voices.push_back(child_voice());
    voices.push_back(gravel_voice());
    voices.push_back(whisper_voice());
    return voices;
}

bool voice_by_name(const std::string& name, VoiceProfile& output) {
    std::vector<VoiceProfile> voices = all_voices();
    for (std::size_t index = 0; index < voices.size(); ++index) {
        if (voices[index].name == name) {
            output = voices[index];
            return true;
        }
    }
    return false;
}

}
