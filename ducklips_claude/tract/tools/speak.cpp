// Speak a score through the segment-module body.
// usage: speak SCORE OUT.wav [seed]
//   SCORE is the body language written by the front end:
//     say "text" x.wav --score SCORE --plan-only      (in ../voice)
// Environment: TRACE=ms prints the cord, breath and narrowest place;
// BRANCHES=pus keeps only the named side branches.

#include "tract/speaker.hpp"
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

}  // namespace

// Flush subnormal numbers to zero: decaying states (silent tubes, a released
// breath) otherwise pass through subnormals, which are slow on every CPU.
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
    std::fprintf(stderr, "usage: speak SCORE OUT.wav [seed]\n");
    return 2;
  }
  voice::Score score;
  if (!read_score(values[1], score)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  std::uint32_t seed = count > 3 ? static_cast<std::uint32_t>(std::strtoul(values[3], nullptr, 10)) : 73021U;
  double trace_ms = std::getenv("TRACE") != nullptr ? std::strtod(std::getenv("TRACE"), nullptr) : 0.0;
  tract::BodyLayout layout = tract::speech_layout();
  layout.seed = seed;
  // Diagnostic switches: FRICATION, ASPIRATION, WALLVAR, TREMOR gains.
  if (std::getenv("FRICATION") != nullptr) {
    layout.frication = std::strtod(std::getenv("FRICATION"), nullptr);
  }
  if (std::getenv("JETATTACK") != nullptr) {
    layout.jet_attack = std::strtod(std::getenv("JETATTACK"), nullptr);
  }
  if (std::getenv("ASPIRATION") != nullptr) {
    layout.aspiration = std::strtod(std::getenv("ASPIRATION"), nullptr);
  }
  if (std::getenv("FLUX") != nullptr) {
    layout.jet_flux = std::strtod(std::getenv("FLUX"), nullptr);
  }
  if (std::getenv("CLOSING") != nullptr) {
    layout.closing_turbulence = std::strtod(std::getenv("CLOSING"), nullptr);
  }
  if (std::getenv("WALLVAR") != nullptr) {
    layout.wall_variation = std::strtod(std::getenv("WALLVAR"), nullptr);
  }
  if (std::getenv("TREMOR") != nullptr) {
    layout.tremor = std::strtod(std::getenv("TREMOR"), nullptr);
    layout.adductor_tremor = 0.0;
    layout.breath_tremor = 0.0;
  }
  if (std::getenv("SINUSNECK") != nullptr) {
    layout.sinus_neck_scale = std::strtod(std::getenv("SINUSNECK"), nullptr);
  }
  if (std::getenv("SKIN") != nullptr) {
    layout.skin_gain = std::strtod(std::getenv("SKIN"), nullptr);
  }
  if (std::getenv("SKINCORNER") != nullptr) {
    layout.skin_corner = std::strtod(std::getenv("SKINCORNER"), nullptr);
  }
  if (std::getenv("WALLK") != nullptr || std::getenv("WALLLOSS") != nullptr) {
    // WALLK=factor scales every vocal wall's support stiffness; WALLLOSS=x
    // sets its support damping (fraction of critical).
    for (tract::SegmentMaterial& material : layout.vocal_material_list) {
      if (std::getenv("WALLK") != nullptr) {
        material.support_stiffness *= std::strtod(std::getenv("WALLK"), nullptr);
      }
      if (std::getenv("WALLLOSS") != nullptr) {
        material.support_loss = std::strtod(std::getenv("WALLLOSS"), nullptr);
      }
    }
  }
  if (std::getenv("BREATHVAR") != nullptr) {
    layout.breath_variation = std::strtod(std::getenv("BREATHVAR"), nullptr);
  }
  if (std::getenv("TENSIONVAR") != nullptr) {
    layout.tension_variation = std::strtod(std::getenv("TENSIONVAR"), nullptr);
  }
  if (std::getenv("POSTUREVAR") != nullptr) {
    layout.posture_variation = std::strtod(std::getenv("POSTUREVAR"), nullptr);
  }
  if (std::getenv("PITCHFB") != nullptr) {
    layout.pitch_feedback = std::strtod(std::getenv("PITCHFB"), nullptr);
  }
  if (std::getenv("BREATHTREMOR") != nullptr) {
    layout.breath_tremor = std::strtod(std::getenv("BREATHTREMOR"), nullptr);
  }
  if (std::getenv("CEILING") != nullptr) {
    layout.jet_corner_ceiling = std::strtod(std::getenv("CEILING"), nullptr);
  }
  if (std::getenv("LABIO") != nullptr) {
    layout.labiodental = std::strtod(std::getenv("LABIO"), nullptr);
  }
  if (std::getenv("LIPWIDTH") != nullptr) {
    layout.lip_slit_width = std::strtod(std::getenv("LIPWIDTH"), nullptr);
  }
  if (std::getenv("SPREAD") != nullptr) {
    layout.jet_spread = std::strtod(std::getenv("SPREAD"), nullptr);
  }
  if (std::getenv("ROUNDJET") != nullptr) {
    layout.glottal_round_jet = std::string(std::getenv("ROUNDJET")) == "1";
  }
  if (std::getenv("PIRILOSS") != nullptr) {
    layout.piriform_loss = std::strtod(std::getenv("PIRILOSS"), nullptr);
  }
  if (std::getenv("BRANCHES") != nullptr) {
    // BRANCHES=pus: keep piriform (p), under-tongue (u), sinuses (s).
    std::string keep = std::getenv("BRANCHES");
    layout.piriform = keep.find('p') != std::string::npos;
    layout.under_tongue = keep.find('u') != std::string::npos;
    layout.sinuses = keep.find('s') != std::string::npos;
  }
  tract::Speaker speaker(layout, score, seed);
  if (std::getenv("PULL") != nullptr) {
    // PULL=cm: outward pull on the pharyngeal walls during voiced closures.
    speaker.set_closure_pull(std::strtod(std::getenv("PULL"), nullptr));
  }
  if (std::getenv("CONSERVE") != nullptr && std::string(std::getenv("CONSERVE")) == "0") {
    // CONSERVE=0: the frame's air displacement is not volume-conserving.
    speaker.body().vocal().set_conserved_span(0, -1);
  }
  if (std::getenv("CONSERVE") != nullptr && std::string(std::getenv("CONSERVE")) == "all") {
    // CONSERVE=all: no net frame displacement anywhere (diagnostic).
    speaker.body().vocal().set_conserved_span(0, speaker.body().vocal().count() - 1);
  }
  if (std::getenv("VOICEPOST") != nullptr) {
    // VOICEPOST=lower,upper (cm): the voicing posture of the folds.
    speaker.set_voicing_posture(std::strtod(std::getenv("VOICEPOST"), nullptr),
                                std::strtod(std::strchr(std::getenv("VOICEPOST"), ',') + 1, nullptr));
  }
  if (std::getenv("NOSEHF") != nullptr) {
    speaker.body().nose().set_wave_smoothing(std::strtod(std::getenv("NOSEHF"), nullptr));
  }
  if (std::getenv("HFLOSS") != nullptr) {
    speaker.body().vocal().set_wave_smoothing(std::strtod(std::getenv("HFLOSS"), nullptr));
  }
  if (std::getenv("CORDDAMP") != nullptr) {
    // CORDDAMP=lower,upper: the fold edges' damping ratios (VocalTractLab's 0.1, 0.6).
    speaker.body().cord().set_damping(std::strtod(std::getenv("CORDDAMP"), nullptr),
                                      std::strtod(std::strchr(std::getenv("CORDDAMP"), ',') + 1, nullptr));
  }
  if (std::getenv("LUNGFOLLOW") != nullptr) {
    speaker.set_lung_follow(std::strtod(std::getenv("LUNGFOLLOW"), nullptr));
  }
  if (std::getenv("RELEASE") != nullptr) {
    speaker.set_release_speed(std::strtod(std::getenv("RELEASE"), nullptr));
  }
  if (std::getenv("PRESS") != nullptr) {
    // PRESS=cm: extra adduction during voiced closures.
    speaker.set_closure_press(std::strtod(std::getenv("PRESS"), nullptr));
  }
  if (std::getenv("WANDER") != nullptr) {
    // WANDER=gain: scales every motor's wander (0 holds every part steady).
    speaker.set_wander_gain(std::strtod(std::getenv("WANDER"), nullptr));
  }
  if (std::getenv("RELAX") != nullptr) {
    // RELAX=share: pharyngeal wall release during voiced closures.
    speaker.set_closure_relaxation(std::strtod(std::getenv("RELAX"), nullptr));
  }

  Biquad first = lowpass(18000.0, tract::kStepRate, 0.5412);
  Biquad second = lowpass(18000.0, tract::kStepRate, 1.3066);
  long trace_every = trace_ms > 0.0 ? static_cast<long>(trace_ms * 1e-3 * tract::kStepRate) : 0;
  long total = static_cast<long>(std::ceil(score.duration * tract::kStepRate));
  std::vector<double> raw;
  raw.reserve(static_cast<std::size_t>(total / 2 + 1));
  double peak = 0.0;
  // Diagnostics read once, not per step (getenv was 2% of the running time).
  double steps_from = 0.0;
  double steps_to = 0.0;
  if (std::getenv("STEPS") != nullptr) {
    steps_from = std::strtod(std::getenv("STEPS"), nullptr);
    steps_to = std::strtod(std::strchr(std::getenv("STEPS"), ',') + 1, nullptr);
  }
  long profile_step = std::getenv("PROFILE") != nullptr
                          ? static_cast<long>(std::strtod(std::getenv("PROFILE"), nullptr) * tract::kStepRate)
                          : -1;
  // CYCLES=path: one line per glottal cycle, timed from each opening (the
  // glottal area rising through 0.02 cm^2 after falling below 0.01): time
  // s, period ms, peak glottal flow cm^3/s, peak-to-peak of the sound,
  // intended pitch Hz, the ear's log correction, lung pressure Pa, and the
  // cycle's smallest glottal area cm^2 (the chink alone when the folds shut).
  std::FILE* cycles = std::getenv("CYCLES") != nullptr ? std::fopen(std::getenv("CYCLES"), "w") : nullptr;
  bool was_shut = false;
  long last_opening = -1;
  double cycle_flow = 0.0;
  double cycle_high = 0.0;
  double cycle_low = 0.0;
  double cycle_area = 1.0e9;
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  for (long step = 0; step < total; ++step) {
    double sound = speaker.step();
    double filtered = run(second, run(first, sound));
    if (!std::isfinite(filtered)) {
      std::fprintf(stderr, "non-finite output at %.4f s\n", step / tract::kStepRate);
      return 1;
    }
    if (cycles != nullptr) {
      tract::Cord& cord = speaker.body().cord();
      double area = cord.area();
      cycle_flow = std::fmax(cycle_flow, cord.flow());
      cycle_high = std::fmax(cycle_high, filtered);
      cycle_low = std::fmin(cycle_low, filtered);
      cycle_area = std::fmin(cycle_area, area);
      if (area < 0.01) {
        was_shut = true;
      }
      if (was_shut && area > 0.02) {
        if (last_opening >= 0) {
          std::fprintf(cycles, "%.5f %.4f %.2f %.6g %.2f %.4f %.1f %.5f\n", step / tract::kStepRate,
                       1000.0 * static_cast<double>(step - last_opening) / tract::kStepRate, cycle_flow,
                       cycle_high - cycle_low, speaker.body().intended_pitch(), speaker.body().pitch_correction(),
                       speaker.body().lung_pressure() / tract::kPascal, cycle_area);
        }
        last_opening = step;
        was_shut = false;
        cycle_flow = 0.0;
        cycle_high = 0.0;
        cycle_low = 0.0;
        cycle_area = 1.0e9;
      }
    }
    if (step % 2 == 0) {
      raw.push_back(filtered);
      peak = std::fmax(peak, std::fabs(filtered));
    }
    if (steps_from < steps_to) {
      // STEPS=from,to (seconds): every step's flows, for hunting artefacts.
      double from = steps_from;
      double to = steps_to;
      double now = step / tract::kStepRate;
      if (now >= from && now < to) {
        tract::VoiceBody& body = speaker.body();
        tract::Tube& vocal = body.vocal();
        double volume = 0.0;
        double frame_flow = 0.0;
        for (tract::Index node = 0; node < vocal.count(); ++node) {
          volume += vocal.area(node) * tract::kNodeLength;
          double radius = std::sqrt(vocal.rest_area(node) / tract::kPi);
          frame_flow += tract::kTwoPi * radius * tract::kNodeLength * vocal.frame_rate(node);
        }
        std::printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.3f %.3f\n", now, vocal.radiated_flow(),
                    body.nose().radiated_flow(), body.cord().flow(), body.cord().area(),
                    body.trachea().pressure(0), vocal.area(vocal.count() - 1), vocal.pressure(20),
                    body.nose().pressure(0), volume, frame_flow);
      }
    }
    if (step == profile_step) {
      // PROFILE=t: area, flow, pressure and frame speed along the vocal tube.
      tract::Tube& vocal = speaker.body().vocal();
      for (tract::Index node = 0; node < vocal.count(); ++node) {
        std::printf("x %5.2f area %6.3f rest %6.3f p %8.1f flow %8.1f frame %8.3f\n", (node + 0.5) * tract::kNodeLength,
                    vocal.area(node), vocal.rest_area(node), vocal.pressure(node),
                    node + 1 < vocal.count() ? vocal.link_flow(node) : vocal.radiated_flow(), vocal.frame_rate(node));
      }
    }
    if (trace_every > 0 && step % trace_every == 0) {
      tract::VoiceBody& body = speaker.body();
      tract::Tube& vocal = body.vocal();
      tract::Index narrowest = 6;
      for (tract::Index node = 6; node < vocal.count(); ++node) {
        if (vocal.area(node) < vocal.area(narrowest)) {
          narrowest = node;
        }
      }
      std::printf("%.4f glottis %.4f flow %7.1f sub %7.0f supra %7.0f mouth %7.0f narrowest %5.3f at %4.1f lips %6.3f\n",
                  step / tract::kStepRate, body.cord().area(), body.cord().flow(),
                  body.trachea().pressure(body.trachea().count() - 1), vocal.pressure(2), vocal.pressure(vocal.count() - 4),
                  vocal.area(narrowest),
                  (narrowest + 0.5) * tract::kNodeLength, vocal.area(vocal.count() - 1));
    }
  }
  if (cycles != nullptr) {
    std::fclose(cycles);
  }
  double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  std::vector<float> samples;
  samples.reserve(raw.size());
  for (double value : raw) {
    samples.push_back(static_cast<float>(peak > 0.0 ? 0.7 * value / peak : 0.0));
  }
  // The space between the mouth and a microphone: a faint reverberant tail
  // and background, as every recording has (the owner's reading of these
  // sentences sits on a soft floor; ours fell to digital silence between
  // words). ROOM=0 leaves it dry.
  if (std::getenv("ROOM") == nullptr || std::strtod(std::getenv("ROOM"), nullptr) != 0.0) {
    voice::RoomOptions room;
    room.seed = seed;
    voice::apply_room(samples, tract::kAudioRate, room);
  }
  write_wave(values[2], samples, static_cast<int>(tract::kAudioRate));
  std::fprintf(stderr, "%s: %.2f s of speech in %.2f s (%.1fx real time)\n", values[2], score.duration, elapsed,
               score.duration / elapsed);
  return 0;
}
