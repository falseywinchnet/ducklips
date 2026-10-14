// Closed-glottis transfer of one posture through the whole body: a pressure
// pulse enters above the glottis and the radiated response is read on a
// 10 Hz grid (150 Hz - 6 kHz). Prints each resonance peak with its level
// and -3 dB bandwidth, so tract losses can be checked formant by formant.
//
// usage: transfer_probe place area lip_area [protrusion] [gap]
// Environment: WALLS=0 makes every wall rigid; SPECTRUM=1 prints the grid;
// DIFFUSION sets the boundary-layer wave smoothing.

#include "voice/body.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

constexpr int kResponseSamples = 9600;
constexpr int kBins = 586;
constexpr double kFirstHz = 150.0;
constexpr double kStepHz = 10.0;

void apply_part(voice::Body& body, voice::Part part, double a, double b) {
  voice::PartEvent event;
  event.part = part;
  event.a = static_cast<float>(a);
  event.b = static_cast<float>(b);
  body.apply(event);
}

// Half-power width around a peak bin, from linear interpolation of the
// power grid in dB on each side.
double bandwidth_at(const std::vector<double>& level_db, int peak) {
  double half = level_db[static_cast<std::size_t>(peak)] - 3.0;
  double low = -1.0;
  for (int bin = peak; bin > 0; --bin) {
    double here = level_db[static_cast<std::size_t>(bin)];
    double below = level_db[static_cast<std::size_t>(bin - 1)];
    if (below < half) {
      low = bin - (here - half) / (here - below);
      break;
    }
  }
  double high = -1.0;
  for (int bin = peak; bin < kBins - 1; ++bin) {
    double here = level_db[static_cast<std::size_t>(bin)];
    double above = level_db[static_cast<std::size_t>(bin + 1)];
    if (above < half) {
      high = bin + (here - half) / (here - above);
      break;
    }
  }
  if (low < 0.0 || high < 0.0) {
    return -1.0;
  }
  return (high - low) * kStepHz;
}

}  // namespace

int main(int count, char** values) {
  if (count < 4) {
    std::fprintf(stderr, "usage: transfer_probe place area lip_area [protrusion] [gap]\n");
    return 2;
  }
  double place = std::strtod(values[1], nullptr);
  double area = std::strtod(values[2], nullptr);
  double lips = std::strtod(values[3], nullptr);
  double protrusion = count > 4 ? std::strtod(values[4], nullptr) : 0.0;
  double gap = count > 5 ? std::strtod(values[5], nullptr) : -0.05;
  voice::BodyOptions options;
  options.hold_noise = 0.0;
  options.aspiration = 0.0;
  options.frication = 0.0;
  if (std::getenv("WALLS") != nullptr) {
    options.wall_coupling = std::strtod(std::getenv("WALLS"), nullptr);
  }
  if (std::getenv("DIFFUSION") != nullptr) {
    options.wave_diffusion = std::strtod(std::getenv("DIFFUSION"), nullptr);
  }
  if (std::getenv("PIRIFORM") != nullptr) {
    options.piriform_scale = std::strtod(std::getenv("PIRIFORM"), nullptr);
  }
  if (std::getenv("SKIN") != nullptr) {
    options.skin = std::strtod(std::getenv("SKIN"), nullptr);
  }
  voice::Body body(options);
  apply_part(body, voice::Part::tongue, place, area);
  apply_part(body, voice::Part::lips, lips, protrusion);
  apply_part(body, voice::Part::jaw, voice::kReleased, 0.0);
  apply_part(body, voice::Part::glottis, gap, 1.0);
  body.settle();
  std::vector<float> samples(static_cast<std::size_t>(kResponseSamples), 0.0F);
  body.render(samples.data(), 480);
  body.kick(1000.0);
  body.render(samples.data(), kResponseSamples);
  // Divide out the radiation slope (+6 dB/octave) so peaks read as volume
  // velocity resonances.
  std::vector<double> level_db(static_cast<std::size_t>(kBins), 0.0);
  for (int bin = 0; bin < kBins; ++bin) {
    double frequency = kFirstHz + kStepHz * bin;
    double omega = 2.0 * voice::kPi * frequency / voice::kAudioRate;
    double real = 0.0;
    double imaginary = 0.0;
    for (int sample = 72; sample < kResponseSamples; ++sample) {
      double value = samples[static_cast<std::size_t>(sample)];
      real += value * std::cos(omega * sample);
      imaginary -= value * std::sin(omega * sample);
    }
    double power = (real * real + imaginary * imaginary) / (frequency * frequency);
    level_db[static_cast<std::size_t>(bin)] = 10.0 * std::log10(power + 1.0e-30);
  }
  double top = level_db[0];
  for (double value : level_db) {
    top = std::max(top, value);
  }
  bool print_grid = std::getenv("SPECTRUM") != nullptr;
  for (int bin = 2; bin < kBins - 2; ++bin) {
    double here = level_db[static_cast<std::size_t>(bin)];
    if (print_grid && bin % 5 == 0) {
      std::printf("  %5.0f Hz %6.1f dB\n", kFirstHz + kStepHz * bin, here - top);
    }
    if (here > level_db[static_cast<std::size_t>(bin - 1)] && here >= level_db[static_cast<std::size_t>(bin + 1)] &&
        here > level_db[static_cast<std::size_t>(bin - 2)] && here >= level_db[static_cast<std::size_t>(bin + 2)]) {
      std::printf("peak %5.0f Hz  level %6.1f dB  bandwidth %5.0f Hz\n", kFirstHz + kStepHz * bin, here - top,
                  bandwidth_at(level_db, bin));
    }
  }
  return 0;
}
