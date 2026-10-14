// Sustain one posture and report what the body does with it.
// usage: vowel_probe OUT.wav place area lip_area [pitch_hz] [lungs_pa] [gap_cm] [seconds]

#include "voice/body.hpp"
#include "voice/wave_file.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argument_count, char** arguments) {
  if (argument_count < 5) {
    std::fprintf(stderr, "usage: vowel_probe OUT.wav place area lip_area [pitch] [lungs] [gap] [seconds]\n");
    return 2;
  }
  float place = std::strtof(arguments[2], nullptr);
  float area = std::strtof(arguments[3], nullptr);
  float lip_area = std::strtof(arguments[4], nullptr);
  float pitch = argument_count > 5 ? std::strtof(arguments[5], nullptr) : 110.0F;
  float lungs = argument_count > 6 ? std::strtof(arguments[6], nullptr) : 800.0F;
  float gap = argument_count > 7 ? std::strtof(arguments[7], nullptr) : 0.0F;
  double seconds = argument_count > 8 ? std::strtod(arguments[8], nullptr) : 1.5;

  voice::BodyOptions options;
  if (std::getenv("ASP") != nullptr) {
    options.aspiration = std::strtod(std::getenv("ASP"), nullptr);
  }
  if (std::getenv("FRIC") != nullptr) {
    options.frication = std::strtod(std::getenv("FRIC"), nullptr);
  }
  if (std::getenv("SKIN") != nullptr) {
    options.skin = std::strtod(std::getenv("SKIN"), nullptr);
  }
  if (std::getenv("VDEPTH") != nullptr) {
    options.viscous_depth = std::strtod(std::getenv("VDEPTH"), nullptr);
  }
  if (std::getenv("ABD_K") != nullptr) {
    options.abduction_stiffening = std::strtod(std::getenv("ABD_K"), nullptr);
  }
  if (std::getenv("ABD_ON") != nullptr) {
    options.abduction_onset = std::strtod(std::getenv("ABD_ON"), nullptr);
  }
  if (std::getenv("FLATTEN") != nullptr) {
    options.pressed_flattening = std::strtod(std::getenv("FLATTEN"), nullptr);
  }
  if (std::getenv("F_CDAMP") != nullptr) {
    options.fold.contact_damping = std::strtod(std::getenv("F_CDAMP"), nullptr);
  }
  if (std::getenv("MASS_EXP") != nullptr) {
    options.fold_mass_exponent = std::strtod(std::getenv("MASS_EXP"), nullptr);
  }
  if (std::getenv("HOLD") != nullptr) {
    options.hold_noise = std::strtod(std::getenv("HOLD"), nullptr);
  }
  const char* fold_names[] = {"F_LM", "F_UM", "F_BM", "F_LK", "F_UK", "F_SK", "F_BK", "F_LD", "F_UD", "F_BD", "F_CONV",
                              "F_CONTACT", "CHINK", "THICK1", "THICK2"};
  double* fold_fields[] = {&options.fold.lower_mass, &options.fold.upper_mass, &options.fold.body_mass,
                           &options.fold.lower_stiffness, &options.fold.upper_stiffness, &options.fold.shear_stiffness,
                           &options.fold.body_stiffness, &options.fold.lower_damping, &options.fold.upper_damping,
                           &options.fold.body_damping, &options.fold.convergence, &options.fold.contact_scale,
                           &options.chink, &options.lower_thickness, &options.upper_thickness};
  for (int index = 0; index < 15; ++index) {
    if (std::getenv(fold_names[index]) != nullptr) {
      *fold_fields[index] = std::strtod(std::getenv(fold_names[index]), nullptr);
    }
  }
  std::FILE* trace_file = std::getenv("TRACE") != nullptr ? std::fopen(std::getenv("TRACE"), "w") : nullptr;
  voice::Body body(options);
  voice::PartEvent event;
  event.part = voice::Part::tongue;
  event.a = place;
  event.b = area;
  event.move = 0.05F;
  body.apply(event);
  event.part = voice::Part::lips;
  event.a = lip_area;
  event.b = 0.0F;
  body.apply(event);
  event.part = voice::Part::pitch;
  event.a = pitch;
  body.apply(event);
  event.part = voice::Part::glottis;
  event.a = gap;
  event.b = 0.6F;
  body.apply(event);
  body.settle();
  event.part = voice::Part::lungs;
  event.a = lungs;
  event.move = 0.08F;
  body.apply(event);

  voice::Index total = static_cast<voice::Index>(seconds * voice::kAudioRate);
  std::vector<float> samples(static_cast<std::size_t>(total), 0.0F);
  voice::Index block = std::getenv("BLOCK") ? std::atoi(std::getenv("BLOCK")) : 4800;
  double peak = 0.0;
  for (voice::Index start = 0; start < total; start += block) {
    voice::Index count = std::min(block, total - start);
    if (trace_file != nullptr) {
      for (voice::Index index = 0; index < count; ++index) {
        body.render(samples.data() + start + index, 1);
        voice::BodyTrace step_trace = body.trace();
        std::fprintf(trace_file, "%.6g %.6g %.6g %.6g\n", step_trace.glottal_area, step_trace.glottal_flow,
                     step_trace.sub_pressure, static_cast<double>(samples[static_cast<std::size_t>(start + index)]));
      }
    } else {
      body.render(samples.data() + start, count);
    }
    voice::BodyTrace trace = body.trace();
    double block_peak = 0.0;
    for (voice::Index index = start; index < start + count; ++index) {
      block_peak = std::max(block_peak, static_cast<double>(std::abs(samples[static_cast<std::size_t>(index)])));
    }
    peak = std::max(peak, block_peak);
    std::printf("t %.2f heard %.1f Hz tension %.3f ps %.0f Pa flow %.0f area %.4f peak %.4g\n",
                (start + count) / voice::kAudioRate, trace.pitch_heard, trace.tension,
                trace.sub_pressure / voice::kPascal, trace.glottal_flow, trace.glottal_area, block_peak);
  }
  if (std::getenv("AREAS") != nullptr) {
    const voice::Tube& oral = body.oral();
    std::fprintf(stderr, "areas");
    for (voice::Index index = 0; index < oral.count; ++index) {
      std::fprintf(stderr, " %.4f", oral.area[static_cast<std::size_t>(index)]);
    }
    std::fprintf(stderr, "\n");
  }
  if (peak > 0.0) {
    for (float& sample : samples) {
      sample = static_cast<float>(sample / peak * 0.7);
    }
  }
  if (!voice::write_wave(arguments[1], samples, static_cast<int>(voice::kAudioRate))) {
    std::fprintf(stderr, "cannot write %s\n", arguments[1]);
    return 1;
  }
  return 0;
}
