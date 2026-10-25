#include "voice/plan_parameters.hpp"

namespace voice {

namespace {

// Where a declared path lives: a field of PlanOptions or of its consonant
// shapes. (Durations, "duration.<PHONE>", are handled by name.)
struct PlanSlot {
  const char* path = "";
  double PlanOptions::*option = nullptr;
  double ConsonantShapes::*shape = nullptr;
};

PlanSlot option_slot(const char* path, double PlanOptions::*field) {
  PlanSlot slot;
  slot.path = path;
  slot.option = field;
  return slot;
}

PlanSlot shape_slot(const char* path, double ConsonantShapes::*field) {
  PlanSlot slot;
  slot.path = path;
  slot.shape = field;
  return slot;
}

// The declared tree, block by block: timing, the glottis, breath, stop
// releases, intonation, the area planner's mouth, the articulatory body's
// consonant shapes.
const std::vector<PlanSlot>& declared_slots() {
  static const std::vector<PlanSlot> slots = {
      option_slot("timing.rate", &PlanOptions::rate),
      option_slot("timing.diphthong_glide_max", &PlanOptions::diphthong_glide_max),
      option_slot("timing.variation", &PlanOptions::variation),
      option_slot("timing.lead_in", &PlanOptions::lead_in),
      option_slot("timing.seal_lead", &PlanOptions::seal_lead),
      option_slot("timing.vowel_hold", &PlanOptions::vowel_hold),
      option_slot("timing.seal_settle", &PlanOptions::seal_settle),
      option_slot("timing.affricate_hiss_max", &PlanOptions::affricate_hiss_max),
      option_slot("timing.affricate_hiss_at", &PlanOptions::affricate_hiss_at),
      option_slot("timing.release_move", &PlanOptions::release_move),
      option_slot("timing.voiced_release_move", &PlanOptions::voiced_release_move),
      option_slot("timing.burst_scale", &PlanOptions::burst_scale),
      option_slot("timing.release_peel_time", &PlanOptions::release_peel_time),
      option_slot("timing.voice_lag", &PlanOptions::voice_lag),
      option_slot("timing.vowel_reach", &PlanOptions::vowel_reach),
      option_slot("timing.vowel_arrival", &PlanOptions::vowel_arrival),
      option_slot("timing.jaw_lead", &PlanOptions::jaw_lead),
      option_slot("timing.fricative_hold", &PlanOptions::fricative_hold),
      option_slot("timing.final_vowel_scale", &PlanOptions::final_vowel_scale),
      option_slot("timing.final_consonant_scale", &PlanOptions::final_consonant_scale),
      option_slot("timing.final_vowel_floor_ms", &PlanOptions::final_vowel_floor_ms),
      option_slot("timing.comma_vowel_scale", &PlanOptions::comma_vowel_scale),
      option_slot("timing.sonorant_devoicing", &PlanOptions::sonorant_devoicing),
      option_slot("glottis.voicing_gap", &PlanOptions::voicing_gap),
      option_slot("glottis.voiced_fricative_gap", &PlanOptions::voiced_fricative_gap),
      option_slot("glottis.aspirate_gap", &PlanOptions::aspirate_gap),
      option_slot("glottis.voiceless_gap", &PlanOptions::voiceless_gap),
      option_slot("glottis.unstressed_gap", &PlanOptions::unstressed_gap),
      option_slot("glottis.stress_press", &PlanOptions::stress_press),
      option_slot("glottis.final_relax_gap", &PlanOptions::final_relax_gap),
      option_slot("glottis.vowel_compression", &PlanOptions::vowel_compression),
      option_slot("glottis.unstressed_compression", &PlanOptions::unstressed_compression),
      option_slot("glottis.onset_gap", &PlanOptions::onset_gap),
      option_slot("glottis.onset_hold", &PlanOptions::onset_hold),
      option_slot("glottis.onset_settle", &PlanOptions::onset_settle),
      option_slot("glottis.fricative_abduct_lead", &PlanOptions::fricative_abduct_lead),
      option_slot("glottis.abduction_move", &PlanOptions::abduction_move),
      option_slot("glottis.adduction_move", &PlanOptions::adduction_move),
      option_slot("glottis.adduction_lead", &PlanOptions::adduction_lead),
      option_slot("glottis.affricate_voicing_move", &PlanOptions::affricate_voicing_move),
      option_slot("glottis.affricate_hiss_gap", &PlanOptions::affricate_hiss_gap),
      option_slot("glottis.affricate_adduction_move", &PlanOptions::affricate_adduction_move),
      option_slot("glottis.final_devoicing", &PlanOptions::final_devoicing),
      option_slot("breath.speech_pressure", &PlanOptions::speech_pressure),
      option_slot("breath.voiceless_pressure", &PlanOptions::voiceless_pressure),
      option_slot("breath.unstressed_pressure", &PlanOptions::unstressed_pressure),
      option_slot("breath.inhale_pressure", &PlanOptions::inhale_pressure),
      option_slot("breath.syllable_rise", &PlanOptions::syllable_rise),
      option_slot("breath.syllable_decay", &PlanOptions::syllable_decay),
      option_slot("timing.function_vowel_scale", &PlanOptions::function_vowel_scale),
      option_slot("nasal.velum", &PlanOptions::nasal_velum),
      option_slot("nasal.open_lead", &PlanOptions::velum_open_lead),
      option_slot("nasal.open_move", &PlanOptions::velum_open_move),
      option_slot("nasal.close_lead", &PlanOptions::velum_close_lead),
      option_slot("nasal.close_move", &PlanOptions::velum_close_move),
      option_slot("nasal.anticipation", &PlanOptions::nasal_anticipation),
      option_slot("nasal.anticipation_move", &PlanOptions::anticipation_move),
      option_slot("breath.pause_pressure", &PlanOptions::pause_pressure),
      option_slot("breath.breath_lead", &PlanOptions::breath_lead),
      option_slot("breath.breath_rise", &PlanOptions::breath_rise),
      option_slot("breath.vowel_breath_lead", &PlanOptions::vowel_breath_lead),
      option_slot("breath.phrase_fade_tail", &PlanOptions::phrase_fade_tail),
      option_slot("breath.phrase_fade_from", &PlanOptions::phrase_fade_from),
      option_slot("breath.phrase_open_delay", &PlanOptions::phrase_open_delay),
      option_slot("release.release_channel", &PlanOptions::release_channel),
      option_slot("release.release_peel_area", &PlanOptions::release_peel_area),
      option_slot("pitch.base_pitch", &PlanOptions::base_pitch),
      option_slot("pitch.pitch_floor", &PlanOptions::pitch_floor),
      option_slot("pitch.accent", &PlanOptions::accent),
      option_slot("pitch.accent_peak", &PlanOptions::accent_peak),
      option_slot("pitch.accent_fall", &PlanOptions::accent_fall),
      option_slot("pitch.unstressed_drift", &PlanOptions::unstressed_drift),
      option_slot("pitch.unstressed_dip", &PlanOptions::unstressed_dip),
      option_slot("pitch.pitch_spread", &PlanOptions::pitch_spread),
      option_slot("pitch.pitch_median", &PlanOptions::pitch_median),
      option_slot("pitch.lowest_pitch", &PlanOptions::lowest_pitch),
      option_slot("pitch.final_fall", &PlanOptions::final_fall),
      option_slot("pitch.comma_pitch", &PlanOptions::comma_pitch),
      option_slot("mouth.postalveolar_lips", &PlanOptions::postalveolar_lips),
      option_slot("mouth.postalveolar_protrusion", &PlanOptions::postalveolar_protrusion),
      option_slot("mouth.sibilant_lips", &PlanOptions::sibilant_lips),
      option_slot("mouth.sibilant_protrusion", &PlanOptions::sibilant_protrusion),
      option_slot("mouth.lip_press", &PlanOptions::lip_press),
      option_slot("mouth.labial_jaw", &PlanOptions::labial_jaw),
      option_slot("mouth.velar_front", &PlanOptions::velar_front),
      option_slot("mouth.fricative_body_area", &PlanOptions::fricative_body_area),
      shape_slot("articulation.labial_jaw", &ConsonantShapes::labial_jaw),
      shape_slot("articulation.labiodental_aperture", &ConsonantShapes::labiodental_aperture),
      shape_slot("articulation.labiodental_jaw", &ConsonantShapes::labiodental_jaw),
      shape_slot("articulation.dental_x", &ConsonantShapes::dental_x),
      shape_slot("articulation.dental_y", &ConsonantShapes::dental_y),
      shape_slot("articulation.dental_press", &ConsonantShapes::dental_press),
      shape_slot("articulation.dh_press", &ConsonantShapes::dh_press),
      shape_slot("articulation.dental_jaw", &ConsonantShapes::dental_jaw),
      shape_slot("articulation.alveolar_x", &ConsonantShapes::alveolar_x),
      shape_slot("articulation.alveolar_y", &ConsonantShapes::alveolar_y),
      shape_slot("articulation.alveolar_press", &ConsonantShapes::alveolar_press),
      shape_slot("articulation.alveolar_jaw", &ConsonantShapes::alveolar_jaw),
      shape_slot("articulation.tap_press", &ConsonantShapes::tap_press),
      shape_slot("articulation.sibilant_x", &ConsonantShapes::sibilant_x),
      shape_slot("articulation.sibilant_y", &ConsonantShapes::sibilant_y),
      shape_slot("articulation.sibilant_press", &ConsonantShapes::sibilant_press),
      shape_slot("articulation.sibilant_jaw", &ConsonantShapes::sibilant_jaw),
      shape_slot("articulation.sibilant_aperture", &ConsonantShapes::sibilant_aperture),
      shape_slot("articulation.sibilant_protrusion", &ConsonantShapes::sibilant_protrusion),
      shape_slot("articulation.sibilant_groove", &ConsonantShapes::sibilant_groove),
      shape_slot("articulation.postalveolar_x", &ConsonantShapes::postalveolar_x),
      shape_slot("articulation.postalveolar_y", &ConsonantShapes::postalveolar_y),
      shape_slot("articulation.postalveolar_press", &ConsonantShapes::postalveolar_press),
      shape_slot("articulation.postalveolar_hiss_press", &ConsonantShapes::postalveolar_hiss_press),
      shape_slot("articulation.postalveolar_jaw", &ConsonantShapes::postalveolar_jaw),
      shape_slot("articulation.velar_jaw", &ConsonantShapes::velar_jaw),
      shape_slot("articulation.lateral_jaw", &ConsonantShapes::lateral_jaw),
      shape_slot("articulation.rhotic_jaw", &ConsonantShapes::rhotic_jaw),
      shape_slot("articulation.liquid_yields_tip", &ConsonantShapes::liquid_yields_tip),
      shape_slot("articulation.lateral_hold", &ConsonantShapes::lateral_hold),
      shape_slot("articulation.postalveolar_aperture", &ConsonantShapes::postalveolar_aperture),
      shape_slot("articulation.postalveolar_protrusion", &ConsonantShapes::postalveolar_protrusion),
      shape_slot("articulation.postalveolar_groove", &ConsonantShapes::postalveolar_groove),
      shape_slot("articulation.affricate_closure", &ConsonantShapes::affricate_closure),
      shape_slot("articulation.lateral_x", &ConsonantShapes::lateral_x),
      shape_slot("articulation.lateral_y", &ConsonantShapes::lateral_y),
      shape_slot("articulation.lateral_press", &ConsonantShapes::lateral_press),
      shape_slot("articulation.lateral_retract", &ConsonantShapes::lateral_retract),
      shape_slot("articulation.dark_lateral_x", &ConsonantShapes::dark_lateral_x),
      shape_slot("articulation.dark_lateral_y", &ConsonantShapes::dark_lateral_y),
      shape_slot("articulation.rhotic_x", &ConsonantShapes::rhotic_x),
      shape_slot("articulation.rhotic_y", &ConsonantShapes::rhotic_y),
      shape_slot("articulation.rhotic_press", &ConsonantShapes::rhotic_press),
      shape_slot("articulation.rhotic_back", &ConsonantShapes::rhotic_back),
      shape_slot("articulation.rhotic_up", &ConsonantShapes::rhotic_up),
      shape_slot("articulation.rhotic_aperture", &ConsonantShapes::rhotic_aperture),
      shape_slot("articulation.rhotic_protrusion", &ConsonantShapes::rhotic_protrusion),
      shape_slot("articulation.velar_front", &ConsonantShapes::velar_front),
      shape_slot("articulation.velar_back", &ConsonantShapes::velar_back),
      shape_slot("articulation.velar_press", &ConsonantShapes::velar_press),
      shape_slot("articulation.w_aperture", &ConsonantShapes::w_aperture),
      shape_slot("articulation.w_protrusion", &ConsonantShapes::w_protrusion),
      shape_slot("articulation.w_spread", &ConsonantShapes::w_spread),
      shape_slot("articulation.w_round_share", &ConsonantShapes::w_round_share),
      shape_slot("articulation.y_raise", &ConsonantShapes::y_raise),
      shape_slot("articulation.y_round_share", &ConsonantShapes::y_round_share),
      shape_slot("articulation.lead", &ConsonantShapes::lead),
      shape_slot("articulation.tip_close_move", &ConsonantShapes::tip_close_move),
      shape_slot("articulation.dorsum_close_move", &ConsonantShapes::dorsum_close_move),
      shape_slot("articulation.tip_release_move", &ConsonantShapes::tip_release_move),
      shape_slot("articulation.dorsum_release_move", &ConsonantShapes::dorsum_release_move),
      shape_slot("articulation.liquid_move", &ConsonantShapes::liquid_move),
      shape_slot("articulation.affricate_release_move", &ConsonantShapes::affricate_release_move),
  };
  return slots;
}

// The phones a duration factor may name (the planner's phone table).
const std::vector<std::string>& duration_phones() {
  static const std::vector<std::string> phones = {"IY", "IH", "EY", "EH", "AE", "AA", "AO", "OW", "UH", "UW", "AH", "AX", "ER", "AY", "AW", "OY", "P", "B", "T", "D", "DX", "K", "G", "M", "N", "NG", "F", "V", "TH", "DH", "S", "Z", "SH", "ZH", "CH", "JH", "L", "R", "W", "Y", "HH"};
  return phones;
}

constexpr const char* kDurationBlock = "duration.";

bool is_duration_phone(const std::string& name) {
  for (const std::string& phone : duration_phones()) {
    if (phone == name) {
      return true;
    }
  }
  return false;
}

const PlanSlot* find_slot(const std::string& path) {
  for (const PlanSlot& slot : declared_slots()) {
    if (path == slot.path) {
      return &slot;
    }
  }
  return nullptr;
}

}  // namespace

bool apply_plan_parameters(const std::vector<ParameterEntry>& entries, PlanOptions& options, std::string& error) {
  const std::string duration_block = kDurationBlock;
  for (const ParameterEntry& entry : entries) {
    if (entry.path.compare(0, duration_block.size(), duration_block) == 0) {
      std::string phone = entry.path.substr(duration_block.size());
      if (!is_duration_phone(phone)) {
        error = "no such phone: " + entry.path;
        return false;
      }
      options.duration_scale[phone] = entry.value;
      continue;
    }
    const PlanSlot* slot = find_slot(entry.path);
    if (slot == nullptr) {
      error = "no such plan parameter: " + entry.path;
      return false;
    }
    if (slot->option != nullptr) {
      options.*(slot->option) = entry.value;
    } else {
      options.shapes.*(slot->shape) = entry.value;
    }
  }
  return true;
}

void collect_plan_parameters(const PlanOptions& options, std::vector<ParameterEntry>& entries) {
  entries.clear();
  for (const PlanSlot& slot : declared_slots()) {
    ParameterEntry entry;
    entry.path = slot.path;
    entry.value = slot.option != nullptr ? options.*(slot.option) : options.shapes.*(slot.shape);
    entries.push_back(entry);
  }
  for (const std::string& phone : duration_phones()) {
    ParameterEntry entry;
    entry.path = std::string(kDurationBlock) + phone;
    std::map<std::string, double>::const_iterator found = options.duration_scale.find(phone);
    entry.value = found == options.duration_scale.end() ? 1.0 : found->second;
    entries.push_back(entry);
  }
}

}  // namespace voice
