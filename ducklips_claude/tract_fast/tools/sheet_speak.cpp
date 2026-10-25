// Speak a score through the sheet tube.
// usage: sheet_speak SCORE OUT.wav [options]
//   SCORE is the body language written by the front end:
//     say "text" x.wav --score SCORE --plan-only      (in ../voice)
//   --seed N, --frication G, --aspiration G, --hold-noise N, --hold-common C,
//   --base-pitch HZ (the cord's pitch at its base stretch), --pressed M, --lead S
// Environment: GAIN=G scales by G instead of normalising; TRACE=FILE writes the speaker's state every 10 ms (TRACE_EVERY
// samples).
// Writes 16-bit at 48 kHz (the physics rate decimated), normalised to a 0.7 peak.

#include "tract/sheet_speaker.hpp"
#include "voice/room.hpp"
#include "voice/score.hpp"
#include "wave_out.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

bool read_score(const char* path, voice::Score& score) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!voice::parse_score_line(line, score)) {
      std::fprintf(stderr, "bad score line: %s\n", line.c_str());
      return false;
    }
  }
  voice::sort_score(score);
  return true;
}

bool read_options(int count, char** values, std::uint32_t& seed, tract::SheetVoice& voice, bool& room) {
  for (int i = 3; i < count; ++i) {
    const char* name = values[i];
    if (i + 1 >= count) {
      std::fprintf(stderr, "option %s needs a value\n", name);
      return false;
    }
    double value = std::strtod(values[i + 1], nullptr);
    if (std::strcmp(name, "--room") == 0) {
      room = value != 0.0;
    } else if (std::strcmp(name, "--seed") == 0) {
      seed = static_cast<std::uint32_t>(value);
    } else if (std::strcmp(name, "--obstacle") == 0) {
      voice.obstacle_gain = value;
    } else if (std::strcmp(name, "--monopole") == 0) {
      voice.jet_monopole = value != 0.0;
    } else if (std::strcmp(name, "--corner-ceiling") == 0) {
      voice.jet_corner_ceiling = value;
    } else if (std::strcmp(name, "--frication") == 0) {
      voice.frication = value;
    } else if (std::strcmp(name, "--aspiration") == 0) {
      voice.aspiration = value;
    } else if (std::strcmp(name, "--hold-noise") == 0) {
      voice.hold_noise = value;
    } else if (std::strcmp(name, "--hold-common") == 0) {
      voice.hold_common = value;
    } else if (std::strcmp(name, "--base-pitch") == 0) {
      voice.base_pitch = value;
    } else if (std::strcmp(name, "--cord-thickness") == 0) {
      voice.cord_thickness = value;
    } else if (std::strcmp(name, "--cord-mass") == 0) {
      voice.cord_mass = value;
    } else if (std::strcmp(name, "--cord-elasticity") == 0) {
      voice.cord_elasticity = value;
    } else if (std::strcmp(name, "--cord-support") == 0) {
      voice.cord_support = value;
    } else if (std::strcmp(name, "--lung-resistance") == 0) {
      voice.lung_resistance = value;
    } else if (std::strcmp(name, "--subglottal-narrowing") == 0) {
      voice.subglottal_narrowing = value;
    } else if (std::strcmp(name, "--lung-corner") == 0) {
      voice.lung_corner = value;
    } else if (std::strcmp(name, "--split-nose") == 0) {
      voice.split_nose = value != 0.0;
    } else if (std::strcmp(name, "--skin-corner") == 0) {
      voice.skin_corner = value;
    } else if (std::strcmp(name, "--skin") == 0) {
      voice.skin_radiation = value;
    } else if (std::strcmp(name, "--closure-volume") == 0) {
      voice.closure_volume = value;
    } else if (std::strcmp(name, "--expansion") == 0) {
      voice.closure_expansion = value;
    } else if (std::strcmp(name, "--tract-mass") == 0) {
      voice.tract_mass = value;
    } else if (std::strcmp(name, "--wall-loss") == 0) {
      voice.wall_loss = value;
    } else if (std::strcmp(name, "--subglottal-taps") == 0) {
      if (i + 3 >= count) {
        return false;
      }
      voice.subglottal_taps[0] = value;
      voice.subglottal_taps[1] = std::strtod(values[i + 2], nullptr);
      voice.subglottal_taps[2] = std::strtod(values[i + 3], nullptr);
      i += 2;
    } else if (std::strcmp(name, "--oral-taps") == 0) {
      if (i + 3 >= count) {
        return false;
      }
      voice.oral_taps[0] = value;
      voice.oral_taps[1] = std::strtod(values[i + 2], nullptr);
      voice.oral_taps[2] = std::strtod(values[i + 3], nullptr);
      i += 2;
    } else if (std::strcmp(name, "--nasal-taps") == 0) {
      if (i + 3 >= count) {
        return false;
      }
      voice.nasal_taps[0] = value;
      voice.nasal_taps[1] = std::strtod(values[i + 2], nullptr);
      voice.nasal_taps[2] = std::strtod(values[i + 3], nullptr);
      i += 2;
    } else if (std::strcmp(name, "--pitch-feedback") == 0) {
      voice.pitch_feedback = value;
    } else if (std::strcmp(name, "--pitch-limit") == 0) {
      voice.pitch_correction_limit = value;
    } else if (std::strcmp(name, "--nasal-length") == 0) {
      voice.nasal_length = value;
    } else if (std::strcmp(name, "--nasal-wall-loss") == 0) {
      voice.nasal_wall_loss = value;
    } else if (std::strcmp(name, "--nasal-mass") == 0) {
      voice.nasal_mass = value;
    } else if (std::strcmp(name, "--nasal-lowpass") == 0) {
      voice.nasal_lowpass = value;
    } else if (std::strcmp(name, "--release-speed") == 0) {
      voice.release_speed = value;
    } else if (std::strcmp(name, "--closure-release") == 0) {
      voice.closure_release = value;
    } else if (std::strcmp(name, "--closure-pull") == 0) {
      voice.closure_pull = value;
    } else if (std::strcmp(name, "--press") == 0) {
      voice.closure_press = value;
    } else if (std::strcmp(name, "--velum-gain") == 0) {
      voice.velum_gain = value;
    } else if (std::strcmp(name, "--under-at") == 0) {
      voice.under_tongue_at = value;
    } else if (std::strcmp(name, "--sinus-radius") == 0) {
      voice.sinus_radius = value;
    } else if (std::strcmp(name, "--release") == 0) {
      voice.breath_release_speed = value;
    } else if (std::strcmp(name, "--fastest") == 0) {
      voice.fastest = value;
    } else if (std::strcmp(name, "--shape-rate") == 0) {
      voice.shape_rate = value;
    } else if (std::strcmp(name, "--control-rate") == 0) {
      voice.control_rate = value;
    } else if (std::strcmp(name, "--listen") == 0) {
      voice.listen = static_cast<int>(value);
    } else if (std::strcmp(name, "--rate") == 0) {
      voice.rate = value;
    } else if (std::strcmp(name, "--modules") == 0) {
      voice.modules = static_cast<int>(value);
    } else if (std::strcmp(name, "--triangular") == 0) {
      voice.triangular_glottis = value != 0.0;
    } else if (std::strcmp(name, "--elongation") == 0) {
      voice.fold_elongation = value;
    } else if (std::strcmp(name, "--elongation-reference") == 0) {
      voice.fold_elongation_reference = value;
    } else if (std::strcmp(name, "--elongation-release") == 0) {
      voice.fold_elongation_release = value;
    } else if (std::strcmp(name, "--devoicing-raise") == 0) {
      voice.devoicing_raise = value;
    } else if (std::strcmp(name, "--devoicing-relax") == 0) {
      voice.devoicing_relax = value;
    } else if (std::strcmp(name, "--elongation-floor") == 0) {
      voice.fold_elongation_floor = value;
    } else if (std::strcmp(name, "--lateral") == 0) {
      voice.lateral_pocket = value != 0.0;
    } else if (std::strcmp(name, "--recesses") == 0) {
      voice.recesses = value != 0.0;
    } else if (std::strcmp(name, "--side-tubes") == 0) {
      voice.side_tubes = value != 0.0;
    } else if (std::strcmp(name, "--tissue") == 0) {
      voice.tissue_relaxation = value;
    } else if (std::strcmp(name, "--jet-loss") == 0) {
      voice.jet_loss = value != 0.0;
    } else if (std::strcmp(name, "--taper") == 0) {
      voice.cord_taper = value;
    } else if (std::strcmp(name, "--contact") == 0) {
      voice.cord_contact = value;
    } else if (std::strcmp(name, "--cord-loss") == 0) {
      voice.cord_loss = value;
    } else if (std::strcmp(name, "--lung-scale") == 0) {
      voice.lung_scale = value;
    } else if (std::strcmp(name, "--velum-speed") == 0) {
      voice.velum_speed = value;
    } else if (std::strcmp(name, "--glottis-speed") == 0) {
      voice.glottis_speed = value;
    } else if (std::strcmp(name, "--speed") == 0) {
      voice.articulator_speed = value;
    } else if (std::strcmp(name, "--voicing-lead") == 0) {
      voice.voicing_lead = value;
    } else if (std::strcmp(name, "--voicing-speed") == 0) {
      voice.voicing_speed = value;
    } else if (std::strcmp(name, "--laryngeal-lead") == 0) {
      voice.laryngeal_lead = value;
    } else if (std::strcmp(name, "--lead") == 0) {
      voice.lead = value;
    } else if (std::strcmp(name, "--voiced-open") == 0) {
      voice.voiced_open_gain = value;
    } else if (std::strcmp(name, "--voicing-lower") == 0) {
      voice.voicing_lower = value;
    } else if (std::strcmp(name, "--voicing-upper") == 0) {
      voice.voicing_upper = value;
    } else if (std::strcmp(name, "--pressed") == 0) {
      voice.pressed = value;
    } else {
      std::fprintf(stderr, "unknown option %s\n", name);
      return false;
    }
    ++i;
  }
  return true;
}

}  // namespace

void flush_subnormals() {
#if defined(__aarch64__)
  std::uint64_t control = 0;
  __asm__ volatile("mrs %0, fpcr" : "=r"(control));
  control |= (1ULL << 24);
  __asm__ volatile("msr fpcr, %0" : : "r"(control));
#endif
}

int main(int count, char** values) {
  flush_subnormals();
  if (count < 3) {
    std::fprintf(stderr, "usage: sheet_speak SCORE OUT.wav [options]\n");
    return 2;
  }
  voice::Score score;
  if (!read_score(values[1], score)) {
    std::fprintf(stderr, "cannot read score %s\n", values[1]);
    return 1;
  }
  std::uint32_t seed = 1U;
  tract::SheetVoice voice;
  const char* gain_text = std::getenv("GAIN");  // a fixed gain instead of normalising (comparisons)
  bool room = true;  // the front end's faint room (a recording is never dry); --room 0 leaves it dry
  const char* trace_path = std::getenv("TRACE");
  std::FILE* trace = trace_path != nullptr ? std::fopen(trace_path, "w") : nullptr;
  const char* every_text = std::getenv("TRACE_EVERY");
  std::size_t every = every_text != nullptr ? static_cast<std::size_t>(std::strtol(every_text, nullptr, 10)) : 480;
  if (!read_options(count, values, seed, voice, room)) {
    return 2;
  }
  tract::SheetSpeaker speaker(score, seed, voice);
  std::vector<float> samples;
  samples.reserve(static_cast<std::size_t>(score.duration * voice.rate) + 16);
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  while (!speaker.finished()) {
    samples.push_back(static_cast<float>(speaker.step()));
    if (trace != nullptr && samples.size() % every == 0) {
      speaker.trace(trace);
    }
    if (speaker.fault()) {
      std::fprintf(stderr, "fault at %.4f s\n", speaker.time());
      break;
    }
  }
  if (trace != nullptr) {
    std::fclose(trace);
  }
  double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  // Down to 48 kHz: two Butterworth biquads at 18 kHz, then every other sample.
  if (voice.rate > 60000.0) {
    Biquad first = lowpass(18000.0, voice.rate, 0.5412);
    Biquad second = lowpass(18000.0, voice.rate, 1.3066);
    std::vector<float> decimated;
    decimated.reserve(samples.size() / 2 + 1);
    for (std::size_t i = 0; i < samples.size(); ++i) {
      double filtered = run(second, run(first, static_cast<double>(samples[i])));
      if (i % 2 == 0) {
        decimated.push_back(static_cast<float>(filtered));
      }
    }
    samples.swap(decimated);
  }
  // A 40 Hz high-pass (two poles), as a recording chain: the walls' slow
  // holding wander pumps air below 50 Hz, physical but nearly inaudible,
  // and it set the normalisation's peak.
  {
    double corner = std::tan(3.14159265358979323846 * 40.0 / 48000.0);
    double norm = 1.0 / (1.0 + std::sqrt(2.0) * corner + corner * corner);
    double b0 = norm;
    double a1 = 2.0 * (corner * corner - 1.0) * norm;
    double a2 = (1.0 - std::sqrt(2.0) * corner + corner * corner) * norm;
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    for (float& sample : samples) {
      double x0 = static_cast<double>(sample);
      double y0 = b0 * (x0 - 2.0 * x1 + x2) - a1 * y1 - a2 * y2;
      x2 = x1;
      x1 = x0;
      y2 = y1;
      y1 = y0;
      sample = static_cast<float>(y0);
    }
  }
  double peak = 1e-30;
  for (float sample : samples) {
    peak = std::fmax(peak, std::fabs(static_cast<double>(sample)));
  }
  double gain = gain_text != nullptr ? std::strtod(gain_text, nullptr) : 0.7 / peak;
  for (float& sample : samples) {
    sample = static_cast<float>(gain * static_cast<double>(sample));
  }
  if (room) {
    // The space between the mouth and the microphone, as the front end's
    // recordings are made: a faint reverberant tail and background.
    voice::RoomOptions options;
    options.seed = seed;
    voice::apply_room(samples, 48000.0, options);
  }
  write_wave(values[2], samples, 48000);
  std::fprintf(stderr, "%s: %.2f s in %.2f s (%.1fx real time)\n", values[2], score.duration, elapsed,
               score.duration / elapsed);
  return 0;
}
