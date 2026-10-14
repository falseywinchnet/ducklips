// Codebook from the engine itself: for each mouth shape, tap the closed-glottis
// tract with a pressure impulse and record the radiated response on the
// 1/24-octave grid (100 Hz - 7.6 kHz), in dB. Walls, piriform recesses,
// trachea and floor space all take part, exactly as in speech.
//
// usage: engine_codebook OUT_PREFIX SHARD SHARDS
// Writes OUT_PREFIX_SHARD_params.npy (n x 5: tongue place, tongue area, lip
// area, protrusion, jaw closure) and OUT_PREFIX_SHARD_response.npy.

#include "time_frequency.hpp"
#include "voice/body.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kResponseSamples = 4096;

struct Setting {
  double place = 0.0;
  double area = 0.0;
  double lips = 0.0;
  double protrusion = 0.0;
  double jaw = 0.0;
};

std::vector<Setting> all_settings() {
  std::vector<Setting> settings;
  for (double place = 3.0; place <= 15.01; place += 0.5) {
    for (int a = 0; a < 12; ++a) {
      for (int l = 0; l < 12; ++l) {
        for (double protrusion : {0.0, 0.6}) {
          for (double jaw : {0.0, 1.0}) {
            Setting setting;
            setting.place = place;
            setting.area = 0.15 * std::pow(4.5 / 0.15, a / 11.0);
            setting.lips = 0.1 * std::pow(12.0 / 0.1, l / 11.0);
            setting.protrusion = protrusion;
            setting.jaw = jaw;
            settings.push_back(setting);
          }
        }
      }
    }
  }
  return settings;
}

void apply_part(voice::Body& body, voice::Part part, double a, double b) {
  voice::PartEvent event;
  event.part = part;
  event.a = static_cast<float>(a);
  event.b = static_cast<float>(b);
  body.apply(event);
}

}  // namespace

int main(int count, char** values) {
  if (count < 4) {
    std::fprintf(stderr, "usage: engine_codebook OUT_PREFIX SHARD SHARDS\n");
    return 2;
  }
  int shard = std::atoi(values[2]);
  int shards = std::atoi(values[3]);
  std::vector<double> hz;
  for (int index = 0;; ++index) {
    double value = 100.0 * std::pow(2.0, index / 24.0);
    if (value > 7600.0) {
      break;
    }
    hz.push_back(value);
  }
  std::vector<Setting> settings = all_settings();
  std::vector<float> parameters;
  std::vector<float> responses;
  std::vector<float> samples(static_cast<std::size_t>(kResponseSamples), 0.0F);
  std::vector<float> settle(480, 0.0F);
  int written = 0;
  voice::BodyOptions options;
  options.hold_noise = 0.0;
  options.aspiration = 0.0;
  options.frication = 0.0;
  for (std::size_t index = static_cast<std::size_t>(shard); index < settings.size();
       index += static_cast<std::size_t>(shards)) {
    const Setting& setting = settings[index];
    voice::Body body(options);
    apply_part(body, voice::Part::tongue, setting.place, setting.area);
    apply_part(body, voice::Part::lips, setting.lips, setting.protrusion);
    apply_part(body, voice::Part::jaw, voice::kReleased, setting.jaw);
    apply_part(body, voice::Part::glottis, -0.05, 1.0);
    body.settle();
    body.render(settle.data(), 480);
    body.kick(1000.0);
    body.render(samples.data(), kResponseSamples);
    parameters.push_back(static_cast<float>(setting.place));
    parameters.push_back(static_cast<float>(setting.area));
    parameters.push_back(static_cast<float>(setting.lips));
    parameters.push_back(static_cast<float>(setting.protrusion));
    parameters.push_back(static_cast<float>(setting.jaw));
    // Skip the direct pulse (1.5 ms); DFT at each grid frequency.
    for (double frequency : hz) {
      double omega = 2.0 * kPi * frequency / voice::kAudioRate;
      double real = 0.0;
      double imaginary = 0.0;
      for (int sample = 72; sample < kResponseSamples; ++sample) {
        double value = samples[static_cast<std::size_t>(sample)];
        real += value * std::cos(omega * sample);
        imaginary -= value * std::sin(omega * sample);
      }
      responses.push_back(static_cast<float>(10.0 * std::log10(real * real + imaginary * imaginary + 1e-30)));
    }
    ++written;
  }
  std::string prefix = std::string(values[1]) + "_" + values[2];
  bool good = tf::write_npy(prefix + "_params.npy", parameters, written, 5) &&
              tf::write_npy(prefix + "_response.npy", responses, written, static_cast<std::int32_t>(hz.size()));
  return good ? 0 : 1;
}
