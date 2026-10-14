// Closed-glottis impulse response of one posture: the tract's own
// resonances and their damping, with walls, nose, floor and sinuses.
// usage: tract_probe OUT.wav place area lip_area [velum_area] [tip_place tip_area]

#include "voice/body.hpp"
#include "voice/wave_file.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argument_count, char** arguments) {
  if (argument_count < 5) {
    std::fprintf(stderr, "usage: tract_probe OUT.wav place area lip_area [velum] [tip_place tip_area]\n");
    return 2;
  }
  voice::BodyOptions options;
  options.hold_noise = 0.0;
  if (std::getenv("WALLS") != nullptr) {
    options.wall_coupling = std::strtod(std::getenv("WALLS"), nullptr);
  }
  voice::Body body(options);
  voice::PartEvent event;
  event.part = voice::Part::tongue;
  event.a = std::strtof(arguments[2], nullptr);
  event.b = std::strtof(arguments[3], nullptr);
  body.apply(event);
  event.part = voice::Part::lips;
  event.a = std::strtof(arguments[4], nullptr);
  event.b = 0.0F;
  body.apply(event);
  if (argument_count > 5) {
    event.part = voice::Part::velum;
    event.a = std::strtof(arguments[5], nullptr);
    body.apply(event);
  }
  if (argument_count > 7) {
    event.part = voice::Part::tip;
    event.a = std::strtof(arguments[6], nullptr);
    event.b = std::strtof(arguments[7], nullptr);
    body.apply(event);
  }
  event.part = voice::Part::glottis;
  event.a = -0.05F;
  event.b = 1.0F;
  body.apply(event);
  event.part = voice::Part::lungs;
  event.a = 0.0F;
  body.apply(event);
  body.settle();
  std::vector<float> samples(24000, 0.0F);
  body.render(samples.data(), 2400);
  body.kick(1000.0);
  body.render(samples.data() + 2400, 24000 - 2400);
  std::vector<float> response(samples.begin() + 2400, samples.end());
  for (std::size_t index = 0; index < response.size(); index += 1200) {
    std::fprintf(stderr, "%zu:%g ", index, static_cast<double>(response[index]));
  }
  std::fprintf(stderr, "\n");
  float peak = 0.0F;
  for (float value : response) {
    peak = std::max(peak, std::abs(value));
  }
  for (float& value : response) {
    value = value / peak * 0.9F;
  }
  if (!voice::write_wave(arguments[1], response, static_cast<int>(voice::kAudioRate))) {
    return 1;
  }
  return 0;
}
