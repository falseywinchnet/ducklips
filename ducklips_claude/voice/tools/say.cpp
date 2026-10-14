// English text in, speech out.
//
//   say "text" OUT.wav [--score FILE] [--rate R] [--pitch HZ] [--seed N]
//   say --from-score FILE OUT.wav
//
// The planner writes the body language (printable with --score); the body
// receives each event when its time arrives, as a stream.

#include "voice/body.hpp"
#include "voice/planner.hpp"
#include "voice/room.hpp"
#include "voice/text.hpp"
#include "voice/wave_file.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

struct Arguments {
  std::string text;
  std::string output;
  std::string score_path;
  std::string from_score;
  voice::PlanOptions plan;
  voice::BodyOptions body;
  voice::RoomOptions room;
  std::string features_path;  // per-phone planning features, for prosody fitting
  std::string audit_path;     // the body's realized state at 4 kHz, for tools/audit.py
  bool plan_only = false;      // stop after planning (no synthesis)
  bool rule_pitch = false;     // ignore data/pitch_en.txt and use the hand rules
  bool rule_durations = false; // ignore data/durations_en.txt and use the hand rules
  bool dry = false;
  bool phones_loaded = false;
};

// One line per planned phone: everything prosody depends on.
void write_features(const std::string& path, const std::vector<voice::TimedPhone>& phones) {
  std::ofstream stream(path);
  stream << "# name phoneme stress start duration pitch word_start word_end accent nuclear phrase_final word\n";
  for (const voice::TimedPhone& phone : phones) {
    char line[256];
    std::snprintf(line, sizeof(line), "%s %s %d %.4f %.4f %.1f %d %d %d %d %d %s\n", phone.name.c_str(),
                  phone.phoneme.empty() ? phone.name.c_str() : phone.phoneme.c_str(), phone.stress, phone.start,
                  phone.duration, phone.pitch, phone.word_start ? 1 : 0, phone.word_end ? 1 : 0, phone.accent ? 1 : 0,
                  phone.nuclear ? 1 : 0, phone.phrase_final ? 1 : 0, phone.word.empty() ? "-" : phone.word.c_str());
    stream << line;
  }
}

// Outside prosody: one line per phone, "duration_seconds pitch_hz" (0 keeps
// the plan's), and "pause seconds" lines for the gaps after each phrase.
bool load_timing(const char* path, voice::PlanOptions& plan) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string first;
  while (input >> first) {
    if (first == "pause") {
      double seconds = 0.0;
      input >> seconds;
      plan.pauses.push_back(seconds);
      continue;
    }
    voice::PhoneTiming timing;
    timing.duration = std::strtod(first.c_str(), nullptr);
    input >> timing.pitch;
    plan.timing.push_back(timing);
  }
  return true;
}

bool parse_arguments(int count, char** values, Arguments& arguments) {
  std::vector<std::string> positional;
  for (int index = 1; index < count; ++index) {
    std::string value = values[index];
    bool has_next = index + 1 < count;
    if (value == "--score" && has_next) {
      arguments.score_path = values[++index];
    } else if (value == "--from-score" && has_next) {
      arguments.from_score = values[++index];
    } else if (value == "--rate" && has_next) {
      arguments.plan.rate = std::strtod(values[++index], nullptr);
    } else if (value == "--pitch" && has_next) {
      arguments.plan.base_pitch = std::strtod(values[++index], nullptr);
    } else if (value == "--seed" && has_next) {
      arguments.body.seed = static_cast<std::uint32_t>(std::strtoul(values[++index], nullptr, 10));
      arguments.plan.seed = arguments.body.seed;
    } else if (value == "--hold" && has_next) {
      arguments.body.hold_noise = std::strtod(values[++index], nullptr);
    } else if (value == "--aspiration" && has_next) {
      arguments.body.aspiration = std::strtod(values[++index], nullptr);
    } else if (value == "--frication" && has_next) {
      arguments.body.frication = std::strtod(values[++index], nullptr);
    } else if (value == "--reverb-db" && has_next) {
      arguments.room.reverb_level_db = std::strtod(values[++index], nullptr);
    } else if (value == "--noise-db" && has_next) {
      arguments.room.noise_level_db = std::strtod(values[++index], nullptr);
    } else if (value == "--skin" && has_next) {
      arguments.body.skin = std::strtod(values[++index], nullptr);
    } else if (value == "--release-channel" && has_next) {
      arguments.plan.release_channel = std::strtod(values[++index], nullptr);
    } else if (value == "--release-move" && has_next) {
      arguments.plan.release_move = std::strtod(values[++index], nullptr);
    } else if (value == "--burst-scale" && has_next) {
      arguments.plan.burst_scale = std::strtod(values[++index], nullptr);
    } else if (value == "--voice-lag" && has_next) {
      arguments.plan.voice_lag = std::strtod(values[++index], nullptr);
    } else if (value == "--unstressed-pressure" && has_next) {
      arguments.plan.unstressed_pressure = std::strtod(values[++index], nullptr);
    } else if (value == "--unstressed-gap" && has_next) {
      arguments.plan.unstressed_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--aspirate-gap" && has_next) {
      arguments.plan.aspirate_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--voiced-fricative-gap" && has_next) {
      arguments.plan.voiced_fricative_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--final-devoicing" && has_next) {
      arguments.plan.final_devoicing = std::strtod(values[++index], nullptr);
    } else if (value == "--lateral-volume" && has_next) {
      arguments.body.lateral_volume = std::strtod(values[++index], nullptr);
    } else if (value == "--lateral-neck" && has_next) {
      arguments.body.lateral_neck = std::strtod(values[++index], nullptr);
    } else if (value == "--velar-front" && has_next) {
      arguments.plan.velar_front = std::strtod(values[++index], nullptr);
    } else if (value == "--labial-jaw" && has_next) {
      arguments.plan.labial_jaw = std::strtod(values[++index], nullptr);
    } else if (value == "--jaw-lead" && has_next) {
      arguments.plan.jaw_lead = std::strtod(values[++index], nullptr);
    } else if (value == "--vowel-reach" && has_next) {
      arguments.plan.vowel_reach = std::strtod(values[++index], nullptr);
    } else if (value == "--vowel-arrival" && has_next) {
      arguments.plan.vowel_arrival = std::strtod(values[++index], nullptr);
    } else if (value == "--oral-loss" && has_next) {
      arguments.body.oral_loss = std::strtod(values[++index], nullptr);
    } else if (value == "--unstressed-compression" && has_next) {
      arguments.plan.unstressed_compression = std::strtod(values[++index], nullptr);
    } else if (value == "--proximity-db" && has_next) {
      arguments.room.proximity_db = std::strtod(values[++index], nullptr);
    } else if (value == "--presence-db" && has_next) {
      arguments.room.presence_db = std::strtod(values[++index], nullptr);
    } else if (value == "--voicing-gap" && has_next) {
      arguments.plan.voicing_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--adduction-move" && has_next) {
      arguments.plan.adduction_move = std::strtod(values[++index], nullptr);
    } else if (value == "--adduction-lead" && has_next) {
      arguments.plan.adduction_lead = std::strtod(values[++index], nullptr);
    } else if (value == "--release-peel-area" && has_next) {
      arguments.plan.release_peel_area = std::strtod(values[++index], nullptr);
    } else if (value == "--release-peel-time" && has_next) {
      arguments.plan.release_peel_time = std::strtod(values[++index], nullptr);
    } else if (value == "--wave-diffusion" && has_next) {
      arguments.body.wave_diffusion = std::strtod(values[++index], nullptr);
    } else if (value == "--phrase-fade-tail" && has_next) {
      arguments.plan.phrase_fade_tail = std::strtod(values[++index], nullptr);
    } else if (value == "--phrase-open-delay" && has_next) {
      arguments.plan.phrase_open_delay = std::strtod(values[++index], nullptr);
    } else if (value == "--pause-pressure" && has_next) {
      arguments.plan.pause_pressure = std::strtod(values[++index], nullptr);
    } else if (value == "--abduction-move" && has_next) {
      arguments.plan.abduction_move = std::strtod(values[++index], nullptr);
    } else if (value == "--obstacle-reach" && has_next) {
      arguments.body.obstacle_reach = std::strtod(values[++index], nullptr);
    } else if (value == "--affricate-hiss-at" && has_next) {
      arguments.plan.affricate_hiss_at = std::strtod(values[++index], nullptr);
    } else if (value == "--affricate-hiss-gap" && has_next) {
      arguments.plan.affricate_hiss_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--strouhal" && has_next) {
      arguments.body.strouhal = std::strtod(values[++index], nullptr);
    } else if (value == "--corner-ceiling" && has_next) {
      arguments.body.corner_ceiling = std::strtod(values[++index], nullptr);
    } else if (value == "--variation" && has_next) {
      arguments.plan.variation = std::strtod(values[++index], nullptr);
    } else if (value == "--frication-stages" && has_next) {
      arguments.body.frication_stages = static_cast<int>(std::strtol(values[++index], nullptr, 10));
    } else if (value == "--postalveolar-lips" && has_next) {
      arguments.plan.postalveolar_lips = std::strtod(values[++index], nullptr);
    } else if (value == "--postalveolar-protrusion" && has_next) {
      arguments.plan.postalveolar_protrusion = std::strtod(values[++index], nullptr);
    } else if (value == "--microphone-hz" && has_next) {
      arguments.body.microphone_hz = std::strtod(values[++index], nullptr);
    } else if (value == "--vowel-compression" && has_next) {
      arguments.plan.vowel_compression = std::strtod(values[++index], nullptr);
    } else if (value == "--voiceless-pressure" && has_next) {
      arguments.plan.voiceless_pressure = std::strtod(values[++index], nullptr);
    } else if (value == "--onset-gap" && has_next) {
      arguments.plan.onset_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--onset-hold" && has_next) {
      arguments.plan.onset_hold = std::strtod(values[++index], nullptr);
    } else if (value == "--onset-settle" && has_next) {
      arguments.plan.onset_settle = std::strtod(values[++index], nullptr);
    } else if (value == "--viscous-depth" && has_next) {
      arguments.body.viscous_depth = std::strtod(values[++index], nullptr);
    } else if (value == "--voiced-breath" && has_next) {
      arguments.body.voiced_breath = std::strtod(values[++index], nullptr);
    } else if (value == "--frication-band" && has_next) {
      arguments.body.frication_band = std::strtod(values[++index], nullptr);
    } else if (value == "--glide-ms" && has_next) {
      arguments.body.glide_ms = std::strtod(values[++index], nullptr);
    } else if (value == "--fricative-hold" && has_next) {
      arguments.plan.fricative_hold = std::strtod(values[++index], nullptr);
    } else if (value == "--fricative-body-area" && has_next) {
      arguments.plan.fricative_body_area = std::strtod(values[++index], nullptr);
    } else if (value == "--voiced-release-move" && has_next) {
      arguments.plan.voiced_release_move = std::strtod(values[++index], nullptr);
    } else if (value == "--fricative-abduct-lead" && has_next) {
      arguments.plan.fricative_abduct_lead = std::strtod(values[++index], nullptr);
    } else if (value == "--breath-lead" && has_next) {
      arguments.plan.breath_lead = std::strtod(values[++index], nullptr);
    } else if (value == "--breath-rise" && has_next) {
      arguments.plan.breath_rise = std::strtod(values[++index], nullptr);
    } else if (value == "--final-vowel-scale" && has_next) {
      arguments.plan.final_vowel_scale = std::strtod(values[++index], nullptr);
    } else if (value == "--final-consonant-scale" && has_next) {
      arguments.plan.final_consonant_scale = std::strtod(values[++index], nullptr);
    } else if (value == "--phrase-fade-from" && has_next) {
      arguments.plan.phrase_fade_from = std::strtod(values[++index], nullptr);
    } else if (value == "--vowel-breath-lead" && has_next) {
      arguments.plan.vowel_breath_lead = std::strtod(values[++index], nullptr);
    } else if (value == "--accent-peak" && has_next) {
      arguments.plan.accent_peak = std::strtod(values[++index], nullptr);
    } else if (value == "--accent-fall" && has_next) {
      arguments.plan.accent_fall = std::strtod(values[++index], nullptr);
    } else if (value == "--unstressed-drift" && has_next) {
      arguments.plan.unstressed_drift = std::strtod(values[++index], nullptr);
    } else if (value == "--pitch-glide" && has_next) {
      arguments.plan.pitch_glide = std::strtod(values[++index], nullptr) != 0.0;
    } else if (value == "--unstressed-dip" && has_next) {
      arguments.plan.unstressed_dip = std::strtod(values[++index], nullptr);
    } else if (value == "--comma-pitch" && has_next) {
      arguments.plan.comma_pitch = std::strtod(values[++index], nullptr);
    } else if (value == "--comma-vowel-scale" && has_next) {
      arguments.plan.comma_vowel_scale = std::strtod(values[++index], nullptr);
    } else if (value == "--final-vowel-floor" && has_next) {
      arguments.plan.final_vowel_floor_ms = std::strtod(values[++index], nullptr);
    } else if (value == "--sonorant-devoicing" && has_next) {
      arguments.plan.sonorant_devoicing = std::strtod(values[++index], nullptr);
    } else if (value == "--no-round-ahead") {
      arguments.plan.round_ahead = false;
    } else if (value == "--final-relax-gap" && has_next) {
      arguments.plan.final_relax_gap = std::strtod(values[++index], nullptr);
    } else if (value == "--pitch-floor" && has_next) {
      arguments.plan.pitch_floor = std::strtod(values[++index], nullptr);
    } else if (value == "--wall-coupling" && has_next) {
      arguments.body.wall_coupling = std::strtod(values[++index], nullptr);
    } else if (value == "--no-flow-carry") {
      arguments.body.flow_carry = false;
    } else if (value == "--sink-last" && has_next) {
      arguments.body.sink_last = static_cast<voice::Index>(std::strtol(values[++index], nullptr, 10));
    } else if (value == "--sink-always") {
      arguments.body.sink_closed_only = false;
    } else if (value == "--jet-growth" && has_next) {
      arguments.body.jet_growth = std::strtod(values[++index], nullptr);
    } else if (value == "--jet-decay" && has_next) {
      arguments.body.jet_decay = std::strtod(values[++index], nullptr);
    } else if (value == "--labiodental" && has_next) {
      arguments.body.labiodental = std::strtod(values[++index], nullptr);
    } else if (value == "--abduction-stiffening" && has_next) {
      arguments.body.abduction_stiffening = std::strtod(values[++index], nullptr);
    } else if (value == "--pressed-flattening" && has_next) {
      arguments.body.pressed_flattening = std::strtod(values[++index], nullptr);
    } else if (value == "--fold-damping" && has_next) {
      double damping = std::strtod(values[++index], nullptr);
      arguments.body.fold.lower_damping = damping;
      arguments.body.fold.upper_damping = damping;
      arguments.body.fold.body_damping = damping;
    } else if (value == "--contact-scale" && has_next) {
      arguments.body.fold.contact_scale = std::strtod(values[++index], nullptr);
    } else if (value == "--contact-damping" && has_next) {
      arguments.body.fold.contact_damping = std::strtod(values[++index], nullptr);
    } else if (value == "--abduction-damping" && has_next) {
      arguments.body.abduction_damping = std::strtod(values[++index], nullptr);
    } else if (value == "--abduction-onset" && has_next) {
      arguments.body.abduction_onset = std::strtod(values[++index], nullptr);
    } else if (value == "--chink" && has_next) {
      arguments.body.chink = std::strtod(values[++index], nullptr);
    } else if (value == "--free-front" && has_next) {
      arguments.body.free_front = std::strtod(values[++index], nullptr);
    } else if (value == "--jet-incoherent") {
      arguments.body.jet_incoherent = true;
    } else if (value == "--dry") {
      arguments.dry = true;
    } else if (value == "--phones" && has_next) {
      arguments.phones_loaded = true;
      if (!voice::load_phone_overrides(values[++index])) {
        std::fprintf(stderr, "cannot apply phone overrides\n");
        return false;
      }
    } else if (value == "--timing" && has_next) {
      if (!load_timing(values[++index], arguments.plan)) {
        std::fprintf(stderr, "cannot read timing file\n");
        return false;
      }
    } else if (value == "--audit" && has_next) {
      arguments.audit_path = values[++index];
    } else if (value == "--features" && has_next) {
      arguments.features_path = values[++index];
    } else if (value == "--stress-press" && has_next) {
      arguments.plan.stress_press = std::strtod(values[++index], nullptr);
    } else if (value == "--rule-durations") {
      arguments.rule_durations = true;
    } else if (value == "--rule-pitch") {
      arguments.rule_pitch = true;
    } else if (value == "--final-stop-breath") {
      arguments.plan.final_stop_breath = true;
    } else if (value == "--plan-only") {
      arguments.plan_only = true;
    } else if (value == "--pressure" && has_next) {
      arguments.plan.speech_pressure = std::strtod(values[++index], nullptr);
    } else {
      positional.push_back(value);
    }
  }
  if (!arguments.from_score.empty()) {
    if (positional.size() != 1) {
      return false;
    }
    arguments.output = positional[0];
    return true;
  }
  if (positional.size() != 2) {
    return false;
  }
  arguments.text = positional[0];
  arguments.output = positional[1];
  return true;
}

bool read_score(const std::string& path, voice::Score& score) {
  std::ifstream stream(path);
  if (!stream) {
    return false;
  }
  std::string line;
  while (std::getline(stream, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!voice::parse_score_line(line, score)) {
      std::fprintf(stderr, "bad score line: %s\n", line.c_str());
    }
  }
  voice::sort_score(score);
  return true;
}

}  // namespace

int main(int count, char** values) {
  Arguments arguments;
  if (!parse_arguments(count, values, arguments)) {
    std::fprintf(stderr, "usage: say \"text\" OUT.wav [--score FILE] [--rate R] [--pitch HZ] [--seed N]\n"
                         "       say --from-score FILE OUT.wav\n");
    return 2;
  }
  std::ifstream learned_pitch("data/pitch_en.txt");
  if (learned_pitch && !arguments.rule_pitch && !voice::load_factor_model("data/pitch_en.txt", arguments.plan.pitch_model)) {
    std::fprintf(stderr, "cannot read data/pitch_en.txt\n");
    return 1;
  }
  std::ifstream learned_durations("data/durations_en.txt");
  if (learned_durations && !arguments.rule_durations &&
      !voice::load_factor_model("data/durations_en.txt", arguments.plan.duration_model)) {
    std::fprintf(stderr, "cannot read data/durations_en.txt\n");
    return 1;
  }
  if (!arguments.phones_loaded) {
    // The tuned phone inventory, when present next to the working directory.
    std::ifstream tuned("data/phones_en.txt");
    if (tuned && !voice::load_phone_overrides("data/phones_en.txt")) {
      std::fprintf(stderr, "cannot apply data/phones_en.txt\n");
      return 1;
    }
  }
  voice::Score score;
  if (!arguments.from_score.empty()) {
    if (!read_score(arguments.from_score, score)) {
      std::fprintf(stderr, "cannot read %s\n", arguments.from_score.c_str());
      return 1;
    }
  } else {
    std::vector<voice::Phrase> phrases = voice::parse_text(arguments.text);
    std::fputs(voice::describe_phrases(phrases).c_str(), stdout);
    voice::Plan plan = voice::plan_speech(phrases, arguments.plan);
    score = plan.score;
    if (!arguments.score_path.empty()) {
      std::ofstream stream(arguments.score_path);
      stream << "# phones: name stress start duration pitch\n";
      for (const voice::TimedPhone& phone : plan.phones) {
        char line[128];
        std::snprintf(line, sizeof(line), "# %s%d %.3f %.3f %.0f\n",
                      phone.phoneme.empty() ? phone.name.c_str() : phone.phoneme.c_str(), phone.stress, phone.start,
                      phone.duration, phone.pitch);
        stream << line;
      }
      stream << voice::score_text(score);
    }
    if (!arguments.features_path.empty()) {
      write_features(arguments.features_path, plan.phones);
    }
    if (arguments.plan_only) {
      return 0;
    }
  }

  voice::Body body(arguments.body);
  std::size_t next_event = 0;
  while (next_event < score.parts.size() && score.parts[next_event].time <= 0.0) {
    body.apply(score.parts[next_event]);
    ++next_event;
  }
  body.settle();

  voice::Index total = static_cast<voice::Index>(std::ceil(score.duration * voice::kAudioRate));
  std::vector<float> samples(static_cast<std::size_t>(total), 0.0F);
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  std::ofstream audit;
  if (!arguments.audit_path.empty()) {
    audit.open(arguments.audit_path);
    audit << "# time glottal_area glottal_flow sub_pressure narrowest_area narrowest_place lip_area velum_area "
             "pitch_heard out_rms\n";
  }
  for (voice::Index start = 0; start < total; start += voice::kSamplesPerControl) {
    double now = start / voice::kAudioRate;
    while (next_event < score.parts.size() && score.parts[next_event].time <= now) {
      body.apply(score.parts[next_event]);
      ++next_event;
    }
    voice::Index block = std::min(voice::kSamplesPerControl, total - start);
    // The audit renders the control block in four parts and records what
    // the body actually did after each (the body counts its own control
    // ticks, so the sound is the same).
    if (audit.is_open()) {
      constexpr voice::Index kAuditParts = 4;
      voice::Index part_size = (block + kAuditParts - 1) / kAuditParts;
      for (voice::Index offset = 0; offset < block; offset += part_size) {
        voice::Index count = std::min(part_size, block - offset);
        body.render(samples.data() + start + offset, count);
        double energy = 0.0;
        for (voice::Index index = 0; index < count; ++index) {
          double sample = static_cast<double>(samples[static_cast<std::size_t>(start + offset + index)]);
          energy += sample * sample;
        }
        voice::BodyTrace trace = body.trace();
        char line[256];
        std::snprintf(line, sizeof(line), "%.5f %.5f %.2f %.1f %.5f %.3f %.4f %.4f %.2f %.6f\n",
                      now + (offset + count) / voice::kAudioRate, trace.glottal_area, trace.glottal_flow,
                      trace.sub_pressure, trace.narrowest_area, trace.narrowest_place, trace.lip_area,
                      trace.velum_area, trace.pitch_heard, std::sqrt(energy / static_cast<double>(count)));
        audit << line;
      }
      continue;
    }
    if (std::getenv("TRACEAT") != nullptr && std::abs(now - std::strtod(std::getenv("TRACEAT"), nullptr)) < 0.006) {
      for (voice::Index index = 0; index < block; ++index) {
        body.render(samples.data() + start + index, 1);
        voice::BodyTrace trace = body.trace();
        std::fprintf(stderr, "%.5f glottal %.4f flow %.0f ps %.0f lipflow %.1f out %.4f\n", now + index / voice::kAudioRate, trace.glottal_area, trace.glottal_flow, trace.sub_pressure / 10.0,
                     trace.lip_flow,
                     static_cast<double>(samples[static_cast<std::size_t>(start + index)]));
      }
      continue;
    }
    // FINE="t0:t1" renders that window one sample at a time and prints the
    // narrowest oral section's area, the flow through it, the lip flow and
    // the output (diagnostic for impulses).
    bool fine = false;
    if (std::getenv("FINE") != nullptr) {
      double from = std::strtod(std::getenv("FINE"), nullptr);
      const char* colon = std::strchr(std::getenv("FINE"), ':');
      double to = colon != nullptr ? std::strtod(colon + 1, nullptr) : from + 0.01;
      fine = now + block / voice::kAudioRate > from && now < to;
      if (fine) {
        for (voice::Index offset = 0; offset < block; ++offset) {
          body.render(samples.data() + start + offset, 1);
          const voice::Tube& oral = body.oral();
          voice::Index narrow = 20;
          for (voice::Index index = 20; index < oral.count - 1; ++index) {
            if (oral.area[static_cast<std::size_t>(index)] < oral.area[static_cast<std::size_t>(narrow)]) {
              narrow = index;
            }
          }
          double at = now + offset / voice::kAudioRate;
          if (at >= from && at < to) {
            voice::BodyTrace trace = body.trace();
            std::fprintf(stderr,
                         "%.5f section %d area %.6f flow %.2f pressure %.0f lips %.2f glottis %.4f gflow %.1f "
                         "sub %.0f out %.4f\n",
                         at, narrow, oral.area[static_cast<std::size_t>(narrow)],
                         oral.flow[static_cast<std::size_t>(narrow)],
                         oral.pressure[static_cast<std::size_t>(narrow) - 1], trace.lip_flow, trace.glottal_area,
                         trace.glottal_flow, trace.sub_pressure, static_cast<double>(samples[static_cast<std::size_t>(start + offset)]));
          }
        }
      }
    }
    if (!fine) {
      body.render(samples.data() + start, block);
    }
    if (std::getenv("WATCH") != nullptr && !std::isfinite(samples[static_cast<std::size_t>(start)])) {
      voice::BodyTrace trace = body.trace();
      std::fprintf(stderr, "non-finite at %.3f s: area %g flow %g ps %g tension %g\n", now, trace.glottal_area,
                   trace.glottal_flow, trace.sub_pressure, trace.tension);
      break;
    }
    if (std::getenv("AREAS") != nullptr) {
      double when = std::strtod(std::getenv("AREAS"), nullptr);
      if (now <= when && when < now + block / voice::kAudioRate) {
        const voice::Tube& oral = body.oral();
        std::fprintf(stderr, "areas");
        for (voice::Index index = 0; index < oral.count; ++index) {
          std::fprintf(stderr, " %.4f", oral.area[static_cast<std::size_t>(index)]);
        }
        std::fprintf(stderr, "\n");
      }
    }
    if (std::getenv("CONSTRICTION") != nullptr &&
        (start / voice::kSamplesPerControl) % std::max(std::atoi(std::getenv("CONSTRICTION")), 1) == 0) {
      voice::BodyTrace trace = body.trace();
      std::fprintf(stderr, "%.3f narrowest %.3f cm2 at %.1f cm  glottal %.3f\n", now, trace.narrowest_area,
                   trace.narrowest_place, trace.glottal_area);
    }
    if (std::getenv("WATCH") != nullptr && (start / voice::kSamplesPerControl) % 20 == 0) {
      voice::BodyTrace trace = body.trace();
      std::fprintf(stderr, "%.3f s: area %.4f flow %.0f ps %.0f tension %.3f heard %.0f out %.3g\n", now,
                   trace.glottal_area, trace.glottal_flow, trace.sub_pressure / 10.0, trace.tension,
                   trace.pitch_heard, static_cast<double>(samples[static_cast<std::size_t>(start)]));
    }
  }
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();

  double peak = 0.0;
  bool finite = true;
  for (float sample : samples) {
    if (!std::isfinite(sample)) {
      finite = false;
      continue;
    }
    peak = std::max(peak, static_cast<double>(std::abs(sample)));
  }
  if (!finite) {
    std::fprintf(stderr, "non-finite samples in render\n");
    return 1;
  }
  if (!arguments.dry) {
    voice::apply_room(samples, voice::kAudioRate, arguments.room);
    peak = 0.0;
    for (float sample : samples) {
      peak = std::max(peak, static_cast<double>(std::abs(sample)));
    }
  }
  double gain = peak > 0.0 ? 0.7 / peak : 1.0;
  for (float& sample : samples) {
    sample = static_cast<float>(sample * gain);
  }
  if (!voice::write_wave(arguments.output, samples, static_cast<int>(voice::kAudioRate))) {
    std::fprintf(stderr, "cannot write %s\n", arguments.output.c_str());
    return 1;
  }
  std::fprintf(stderr, "%.2f s of speech in %.2f s (%.1fx realtime), gain %.3g\n", score.duration, seconds,
               score.duration / seconds, gain);
  return 0;
}
